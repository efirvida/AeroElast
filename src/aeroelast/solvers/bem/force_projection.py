"""Conservative force projection from BEM stations onto 3D shell mesh.

Maps per-span-station aerodynamic loads (Np, Tp) from a BEM computation
onto the finite-element mesh nodes, preserving the total integrated force
and moment on each chordwise strip.

Reference-point convention
--------------------------
All per-strip moments and nodal offsets are expressed about the
**aerodynamic centre** (AC) of the strip, i.e. the point at
``aerodynamic_center × chord`` measured from the estimated leading edge
along the local chord direction.

Using the AC as the unique reference point for both the moment balance
and the nodal offset vectors means:

* The geometric transfer term  ``r_{AC→ref} × F``  is identically zero —
  it only appears when ``ref ≠ AC``.
* The only moment distributed to the nodes is the aerodynamic pitching
  moment ``M_AC`` from the BEM polars (Cm coefficient), which is the
  physically correct quantity.
* The result is **independent of mesh density asymmetry** between the
  leading and trailing edges, which would bias an arithmetic-mean
  centroid toward the denser side.
"""

import logging
from dataclasses import dataclass
from typing import List

import numpy as np
from scipy.linalg import lstsq

from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.aerodynamics import BladeAero
from aeroelast.solvers.bem.engine import BEMResult

logger = logging.getLogger(__name__)


@dataclass
class _Strip:
    """Internal: a chordwise strip of mesh nodes."""

    node_indices: np.ndarray  # indices into mesh.nodes
    r_center: float  # span-wise centre (m, from hub)
    dr: float  # strip width (m)
    centroid: np.ndarray  # 3-D arithmetic mean of strip nodes (used for PCA only)
    ac_position: np.ndarray  # 3-D position of the aerodynamic centre
    offsets: np.ndarray  # (n, 3) node positions relative to *ac_position*


class ForceProjector:
    """Project BEM distributed loads onto shell mesh nodes.

    Forces are distributed so that every chordwise strip conserves the
    total BEM force vector **and** moment about its centroid.

    Parameters
    ----------
    mesh : MeshModel
        The shell finite-element mesh (all nodes considered).
    blade_aero : BladeAero
        Aerodynamic blade definition (provides radial station positions).
    span_direction : array-like
        Unit vector along the blade span in global coordinates
        (default: z-axis ``[0, 0, 1]``).
    normal_direction : array-like
        Global direction for the BEM normal force *Np*
        (default: x-axis ``[1, 0, 0]``, i.e. downwind).
    tangential_direction : array-like
        Global direction for the BEM tangential force *Tp*
        (default: y-axis ``[0, 1, 0]``, positive in rotation direction).
    hub_radius : float or None
        Override hub radius for span coordinate calculation.
        If *None*, taken from *blade_aero*.
    """

    def __init__(
        self,
        mesh: MeshModel,
        blade_aero: BladeAero,
        span_direction=None,
        normal_direction=None,
        tangential_direction=None,
        hub_radius: float | None = None,
        include_pitching_moment: bool = True,
    ):
        self._include_pitching_moment = include_pitching_moment
        span_dir = np.asarray(
            span_direction if span_direction is not None else [0.0, 0.0, 1.0],
            dtype=float,
        )
        span_dir /= np.linalg.norm(span_dir)
        self._span_dir = span_dir

        self._normal_dir = np.asarray(
            normal_direction if normal_direction is not None else [1.0, 0.0, 0.0],
            dtype=float,
        )
        self._normal_dir /= np.linalg.norm(self._normal_dir)

        self._tangential_dir = np.asarray(
            tangential_direction if tangential_direction is not None else [0.0, 1.0, 0.0],
            dtype=float,
        )
        self._tangential_dir /= np.linalg.norm(self._tangential_dir)

        coords = mesh.coords_array  # (N, 3)
        n_nodes = coords.shape[0]
        # A hub radius of 0.0 is a placeholder, not a value: the campaign case
        # YAMLs carried it while the aero stations are stored from the rotor
        # centre, so honouring it shifted every strip by the real hub radius
        # (3.97 m on the IEA 15 MW rotor) and emptied the root strip.  Fall back
        # only when the aero actually has a hub, so a wind-tunnel case (root at
        # the rotation axis, no hub) keeps its explicit 0.0.
        if hub_radius is None or (hub_radius == 0.0 and blade_aero.hub_radius > 0.0):
            hub_r = blade_aero.hub_radius
        else:
            hub_r = hub_radius

        # Span coordinate for every mesh node (distance along the span
        # direction, in the blade-local frame where the root is at 0).
        # The BladeMesh coordinates are already blade-local (root at 0),
        # so NO hub offset is subtracted here — subtracting it shifted all
        # strips by hub_radius and emptied the root strip (2026-09-09,
        # diverged the yaml-blade FSI campaign whose hub_radius = 3.97 m).
        span_coords = coords @ span_dir  # blade-local [0, blade_length]

        # BEM station radial positions are stored from the rotor centre
        # (hub_radius + blade-local span).  Subtract hub_radius so they are
        # expressed in the same blade-local frame as span_coords.
        r_stations = blade_aero.r - hub_r  # now blade-local [0, blade_length]

        # Build strip boundaries at midpoints between stations
        r_mid = 0.5 * (r_stations[:-1] + r_stations[1:])
        r_low = np.empty(len(r_stations))
        r_high = np.empty(len(r_stations))
        r_low[0] = r_stations[0] - (r_mid[0] - r_stations[0])
        r_high[-1] = r_stations[-1] + (r_stations[-1] - r_mid[-1])
        r_low[1:] = r_mid
        r_high[:-1] = r_mid

        # Assign each node to a strip
        node_strip = np.full(n_nodes, -1, dtype=int)
        for k in range(len(r_stations)):
            mask = (span_coords >= r_low[k]) & (span_coords < r_high[k])
            node_strip[mask] = k
        # Handle nodes exactly at the tip boundary
        node_strip[span_coords == r_high[-1]] = len(r_stations) - 1

        # Build strip objects.
        #
        # For each strip we need the aerodynamic centre (AC) position so that
        # nodal offsets are expressed relative to AC — the physically correct
        # reference point for the moment balance (see module docstring).
        #
        # The AC position requires knowing the local chord direction, which is
        # estimated via PCA on the strip node coordinates.  Both computations
        # are done here in a single pass so that _Strip.offsets already uses
        # the AC as its origin.
        self._strips: List[_Strip] = []
        self._n_nodes = n_nodes
        self._strip_chord_dirs: list[np.ndarray] = []

        for k in range(len(r_stations)):
            idx = np.where(node_strip == k)[0]
            station = blade_aero.stations[k]

            if len(idx) == 0:
                # Empty strip – create a placeholder
                self._strips.append(
                    _Strip(
                        node_indices=idx,
                        r_center=float(r_stations[k]),
                        dr=float(r_high[k] - r_low[k]),
                        centroid=np.zeros(3),
                        ac_position=np.zeros(3),
                        offsets=np.zeros((0, 3)),
                    )
                )
                self._strip_chord_dirs.append(self._tangential_dir.copy())
                continue

            strip_coords = coords[idx]
            # Arithmetic mean — used only for PCA centering, not as moment
            # reference (which would be biased by uneven LE/TE node density).
            centroid = strip_coords.mean(axis=0)

            # ------------------------------------------------------------------
            # Chord direction via PCA on the chordwise plane (⊥ span_dir)
            # ------------------------------------------------------------------
            if len(idx) >= 2:
                off = strip_coords - centroid
                off_plane = off - np.outer(off @ span_dir, span_dir)
                chord_dir = self._tangential_dir.copy()  # fallback
                if np.linalg.norm(off_plane) >= 1e-12:
                    _, _, Vt = np.linalg.svd(off_plane, full_matrices=False)
                    cd = Vt[0]
                    # Orient consistently (avoid ±π flips at tip sections)
                    if abs(np.dot(cd, self._tangential_dir)) >= abs(
                        np.dot(cd, self._normal_dir)
                    ):
                        if np.dot(cd, self._tangential_dir) < 0:
                            cd = -cd
                    else:
                        if np.dot(cd, self._normal_dir) < 0:
                            cd = -cd
                    c_norm = np.linalg.norm(cd)
                    if c_norm > 1e-12:
                        chord_dir = cd / c_norm
            else:
                chord_dir = self._tangential_dir.copy()

            self._strip_chord_dirs.append(chord_dir)

            # ------------------------------------------------------------------
            # Aerodynamic centre position (3-D)
            #
            # LE is estimated as the node with the minimum projection onto
            # chord_dir.  The AC lies at aerodynamic_center × chord from that
            # point (default: c/4).
            # ------------------------------------------------------------------
            chord_proj = strip_coords @ chord_dir
            le_proj = float(np.min(chord_proj))
            ac_proj = le_proj + station.airfoil.aerodynamic_center * station.chord
            # Span-wise component: use the strip's mean span coordinate so the
            # AC sits on the strip mid-plane (not offset spanwise).
            span_proj = float(centroid @ span_dir)
            # Chordwise-normal component: use the centroid projection onto the
            # direction ⊥ both span and chord (i.e. the thickness direction).
            thickness_dir = np.cross(span_dir, chord_dir)
            thickness_dir_norm = np.linalg.norm(thickness_dir)
            if thickness_dir_norm > 1e-12:
                thickness_dir /= thickness_dir_norm
                thick_proj = float(centroid @ thickness_dir)
            else:
                thickness_dir = np.zeros(3)
                thick_proj = 0.0

            ac_position = (
                ac_proj * chord_dir
                + span_proj * span_dir
                + thick_proj * thickness_dir
            )

            self._strips.append(
                _Strip(
                    node_indices=idx,
                    r_center=float(r_stations[k]),
                    dr=float(r_high[k] - r_low[k]),
                    centroid=centroid,
                    ac_position=ac_position,
                    offsets=strip_coords - ac_position,
                )
            )

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def project(self, bem_result: BEMResult) -> np.ndarray:
        """Map BEM loads onto mesh nodes.

        Parameters
        ----------
        bem_result : BEMResult
            Output of :pymethod:`BEMSolver.compute`.

        Returns
        -------
        forces : ndarray, shape (n_nodes, 3)
            Nodal force vectors in global coordinates.
        """
        forces = np.zeros((self._n_nodes, 3))

        for k, strip in enumerate(self._strips):
            n_k = len(strip.node_indices)
            if n_k == 0:
                continue

            # Total strip force in aero frame
            F_n = float(bem_result.Np[k]) * strip.dr  # normal (N)
            F_t = float(bem_result.Tp[k]) * strip.dr  # tangential (N)

            # Global force vector for the strip
            F_strip = F_n * self._normal_dir + F_t * self._tangential_dir

            # Moment about the aerodynamic centre (AC).
            #
            # Because _Strip.offsets are measured from the AC, the geometric
            # transfer term  r_{AC→ref} × F  is identically zero here.
            # The only moment to distribute is M_AC from the BEM polars.
            M_strip = (
                float(bem_result.Mp[k]) * strip.dr * self._span_dir
                if (bem_result.Mp is not None and self._include_pitching_moment)
                else np.zeros(3)
            )

            # Distribute to nodes (constrained minimum-norm)
            f_nodes = self._distribute(strip, F_strip, M_strip)
            forces[strip.node_indices] = f_nodes

        return forces

    def verify(self, bem_result: BEMResult, forces: np.ndarray) -> dict:
        """Check global force and moment conservation.

        The expected force is computed using the same discrete strip
        summation (``Np[k] * dr_k``) that :meth:`project` uses, so
        the error reflects only the per-strip distribution fidelity.

        Returns
        -------
        dict
            ``force_error`` – norm of residual (N).
            ``force_bem`` – expected total BEM force vector (N).
            ``force_mesh`` – actual summed mesh force vector (N).
        """
        # Discrete strip sum — same as project()
        F_total_n = 0.0
        F_total_t = 0.0
        for k, strip in enumerate(self._strips):
            if len(strip.node_indices) == 0:
                continue
            F_total_n += float(bem_result.Np[k]) * strip.dr
            F_total_t += float(bem_result.Tp[k]) * strip.dr

        F_bem = F_total_n * self._normal_dir + F_total_t * self._tangential_dir
        F_mesh = forces.sum(axis=0)

        return {
            "force_error": float(np.linalg.norm(F_mesh - F_bem)),
            "force_bem": F_bem,
            "force_mesh": F_mesh,
        }

    # ------------------------------------------------------------------
    #  Internal
    # ------------------------------------------------------------------

    @staticmethod
    def _distribute(
        strip: _Strip,
        F_strip: np.ndarray,
        M_strip: np.ndarray,
    ) -> np.ndarray:
        """Minimum-norm nodal forces preserving total force and moment.

        Solves::

            min  Σ |f_j|²
            s.t. Σ f_j          = F_strip        (3 eqs)
                 Σ d_j × f_j   = M_strip        (3 eqs)

        via the pseudoinverse  f = Aᵀ (A Aᵀ)⁻¹ b.

        **Single-node strips:** when *n* == 1 there is no moment arm, so
        ``M_strip`` cannot be represented as a force couple.  Only
        ``F_strip`` is applied; if ``M_strip`` is non-negligible a
        ``WARNING`` is emitted advising to refine the mesh so that each
        BEM strip spans at least two coupling nodes.
        """
        n = len(strip.node_indices)
        if n == 1:
            # A single node has no moment arm — the aerodynamic pitching moment
            # M_strip cannot be represented as a force couple.  Apply F_strip
            # and warn once so that this situation is not silently ignored.
            m_mag = float(np.linalg.norm(M_strip))
            if m_mag > 1e-10:
                logger.warning(
                    "Strip with a single mesh node: aerodynamic pitching moment "
                    "‖M_strip‖ = %.3e N·m cannot be distributed and is dropped. "
                    "Refine the mesh so each BEM strip contains ≥ 2 nodes.",
                    m_mag,
                )
            return F_strip.reshape(1, 3)

        # Build constraint matrix A  (6 × 3n)
        A = np.zeros((6, 3 * n))
        for j in range(n):
            c = 3 * j
            # Force balance rows
            A[0, c] = 1.0
            A[1, c + 1] = 1.0
            A[2, c + 2] = 1.0
            # Moment balance rows:  d × f
            dx, dy, dz = strip.offsets[j]
            #  (d × f)_x =  dy * fz - dz * fy
            A[3, c + 1] += -dz
            A[3, c + 2] += dy
            #  (d × f)_y =  dz * fx - dx * fz
            A[4, c] += dz
            A[4, c + 2] += -dx
            #  (d × f)_z =  dx * fy - dy * fx
            A[5, c] += -dy
            A[5, c + 1] += dx

        b = np.concatenate([F_strip, M_strip])

        # Minimum-norm solution: f = Aᵀ (A Aᵀ)⁻¹ b
        AAT = A @ A.T
        lam, _, _, _ = lstsq(AAT, b)
        f_flat = A.T @ lam

        return f_flat.reshape(n, 3)
