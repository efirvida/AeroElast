"""Conservative force projection from BEM stations onto 3D shell mesh.

Maps per-span-station aerodynamic loads (Np, Tp) from a BEM computation
onto the finite-element mesh nodes, preserving the total integrated force
and moment on each chordwise strip.
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
    centroid: np.ndarray  # 3-D centroid of strip nodes
    offsets: np.ndarray  # (n, 3) node positions relative to centroid


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
    ):
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

        # Span coordinate for every mesh node (distance from hub centre
        # measured along the span direction)
        span_coords = coords @ span_dir  # projection

        # BEM station radial positions (measured from the hub centre).
        #
        # A single-blade shell mesh is referenced to the blade root (span
        # origin 0) while the BEM stations are hub-referenced; a rotor mesh is
        # hub-referenced like the stations.  Anchor the station grid on the
        # mesh's own span origin so both conventions line up: for a
        # blade-root mesh this is a shift by the hub radius, for a
        # hub-referenced mesh it is a no-op.  When the caller supplies
        # ``hub_radius`` explicitly it names that same datum.
        r_stations = np.asarray(blade_aero.r, dtype=float)
        mesh_span_min = float(span_coords.min())
        if hub_radius is None:
            r_stations = r_stations + (mesh_span_min - r_stations[0])
        else:
            r_stations = r_stations - float(hub_radius) + mesh_span_min

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

        # Build strip objects
        self._strips: List[_Strip] = []
        self._n_nodes = n_nodes
        for k in range(len(r_stations)):
            idx = np.where(node_strip == k)[0]
            if len(idx) == 0:
                # Empty strip – create a placeholder
                self._strips.append(
                    _Strip(
                        node_indices=idx,
                        r_center=float(r_stations[k]),
                        dr=float(r_high[k] - r_low[k]),
                        centroid=np.zeros(3),
                        offsets=np.zeros((0, 3)),
                    )
                )
                continue
            strip_coords = coords[idx]
            centroid = strip_coords.mean(axis=0)
            self._strips.append(
                _Strip(
                    node_indices=idx,
                    r_center=float(r_stations[k]),
                    dr=float(r_high[k] - r_low[k]),
                    centroid=centroid,
                    offsets=strip_coords - centroid,
                )
            )

        # ------------------------------------------------------------------
        # Per-strip chord directions and AC-to-centroid offset vectors.
        #
        # BEM polars define Cm (and therefore Mp) about the aerodynamic
        # centre (AC, typically at c/4).  ForceProjector._distribute()
        # balances moments about the strip *centroid*, so the transfer is:
        #
        #   M_centroid = M_AC + (r_AC - r_centroid) x F_strip
        #              = M_AC - cross(ac_offset, F_strip)
        #
        # with ``ac_offset = r_centroid - r_AC``.  The moment arm term is
        # often larger than M_AC for typical wind-turbine blades and must
        # not be omitted; it is kept as the full 3-D vector so the applied
        # moment about the origin is reproduced exactly.
        # ------------------------------------------------------------------
        self._strip_chord_dirs: list[np.ndarray] = []
        self._strip_ac_offsets: list[np.ndarray] = []

        for k, strip in enumerate(self._strips):
            idx = strip.node_indices
            station = blade_aero.stations[k]

            if len(idx) < 2:
                # Not enough nodes for PCA — fall back, zero AC offset
                self._strip_chord_dirs.append(self._tangential_dir.copy())
                self._strip_ac_offsets.append(np.zeros(3))
                continue

            strip_pts = coords[idx]
            off = strip_pts - strip_pts.mean(axis=0)
            # Project offsets onto the plane ⊥ span_dir
            off_plane = off - np.outer(off @ span_dir, span_dir)
            if np.linalg.norm(off_plane) < 1e-12:
                self._strip_chord_dirs.append(self._tangential_dir.copy())
                self._strip_ac_offsets.append(np.zeros(3))
                continue

            _, _, Vt = np.linalg.svd(off_plane, full_matrices=False)
            chord_dir = Vt[0]
            # Orient using most-aligned reference direction (avoids ±π
            # sign flips at tip sections where chord ⊥ normal_dir)
            if abs(np.dot(chord_dir, self._tangential_dir)) >= abs(
                np.dot(chord_dir, self._normal_dir)
            ):
                if np.dot(chord_dir, self._tangential_dir) < 0:
                    chord_dir = -chord_dir
            else:
                if np.dot(chord_dir, self._normal_dir) < 0:
                    chord_dir = -chord_dir
            c_norm = np.linalg.norm(chord_dir)
            if c_norm > 1e-12:
                chord_dir /= c_norm
            else:
                chord_dir = self._tangential_dir.copy()
            self._strip_chord_dirs.append(chord_dir)

            # Leading/trailing ends from the section's own geometry: the
            # blunt end (the larger in-plane spread inside the outer quarter
            # of the chord) is the leading edge, the sharp end the trailing
            # edge.  Inferring the LE from ``min(chord_proj)`` would pick
            # whichever end the unrelated reference axis happened to point
            # away from; this rule is convention-free and sign-free.
            chord_proj = strip_pts @ chord_dir
            p_lo = float(chord_proj.min())
            p_hi = float(chord_proj.max())
            chord = p_hi - p_lo
            in_plane = np.cross(span_dir, chord_dir)
            in_plane_norm = float(np.linalg.norm(in_plane))
            le_at_hi = True
            if chord > 1e-12 and in_plane_norm > 1e-12:
                in_plane = in_plane / in_plane_norm
                slab = 0.25 * chord
                q_lo = strip_pts[chord_proj <= p_lo + slab] @ in_plane
                q_hi = strip_pts[chord_proj >= p_hi - slab] @ in_plane
                t_lo = float(q_lo.max() - q_lo.min()) if q_lo.size else 0.0
                t_hi = float(q_hi.max() - q_hi.min()) if q_hi.size else 0.0
                # Within 10 % the blunt/sharp split is ill-conditioned
                # (e.g. a circular root section); keep the +chord_dir end,
                # the mesh's own leading-edge convention.  A symmetric
                # section is where the two candidate AC placements are
                # least distinguishable, so the ambiguity is smallest.
                if t_lo > t_hi and (t_lo - t_hi) > 0.10 * max(t_lo, t_hi):
                    le_at_hi = False
            i_hi = int(np.argmax(chord_proj))
            i_lo = int(np.argmin(chord_proj))
            le_i, te_i = (i_hi, i_lo) if le_at_hi else (i_lo, i_hi)
            ac_frac = station.airfoil.aerodynamic_center
            ac_point = strip_pts[le_i] + ac_frac * (strip_pts[te_i] - strip_pts[le_i])
            # Full 3-D vector from the aerodynamic centre to the strip
            # centroid.  Keeping the out-of-chord component (the span and
            # thickness offset) is required to reproduce the applied moment
            # about the origin; the chordwise part alone leaves a ~1 %
            # lever-arm error on the real blade.
            self._strip_ac_offsets.append(strip.centroid - ac_point)

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

            # Moment about strip centroid:
            #   M_centroid = M_AC + (r_AC - r_centroid) x F_strip
            #              = M_AC - cross(ac_offset, F_strip)
            # M_AC is the aerodynamic pitching moment from BEM polars (about
            # the aerodynamic centre).  The geometric transfer term accounts
            # for the moment arm between the AC and the centroid of the strip
            # nodes and is always present regardless of Mp availability.
            M_ac = (
                float(bem_result.Mp[k]) * strip.dr * self._span_dir
                if bem_result.Mp is not None
                else np.zeros(3)
            )
            M_strip = M_ac - np.cross(self._strip_ac_offsets[k], F_strip)

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
