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

    Conventions (declared by the model owner, pinned by the downwind
    thrust-direction and positive-power assertions in
    ``tests/test_force_projection_load_frame.py``):

    * span runs along ``+Z``;
    * the section chord runs along ``+X`` on this repo's IEA-15MW mesh,
      with the leading edge at positive ``x``;
    * the out-of-plane (flapwise/section-normal) direction is ``+Y``, which
      is the fluid / downwind direction;
    * the rotor turns clockwise viewed from behind, i.e.
      ``Omega = +omega * y``.

    The mesh axes were **measured** on this tree's IEA-15MW blade mesh: a
    ring's ``x`` extent is the section chord and its ``y`` extent the
    airfoil thickness, and the tip ring's mean ``y`` is the documented
    prebend ``BlCrvAC``.  The *sense* (which of the two opposite
directions along each axis the loads push) is the model owner's convention
    rather than something derivable from the geometry alone.

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
        Sense reference for the BEM normal force *Np* (default: y-axis
        ``[0, 1, 0]``, the fluid/downwind direction).  The *axis* *Np*
        rides is taken per strip from the section geometry (normal to the
        local chord); this vector only selects the sense, i.e. which of the
        two opposite in-plane directions is ``+normal_hat``.  It must be
        (nearly) parallel to the section normal - an orthogonal reference
        makes the sign a round-off decision.
    tangential_direction : array-like
        Sense reference for the BEM tangential force *Tp* (default: x-axis
        ``[1, 0, 0]``, the chordwise direction).  As with
        *normal_direction*, it selects the sense of the per-strip chord
        axis, not the axis itself, and must be (nearly) parallel to that
        axis.
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
            normal_direction if normal_direction is not None else [0.0, 1.0, 0.0],
            dtype=float,
        )
        self._normal_dir /= np.linalg.norm(self._normal_dir)

        self._tangential_dir = np.asarray(
            tangential_direction if tangential_direction is not None else [1.0, 0.0, 0.0],
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
        # Per-strip load frames and AC-to-centroid offset vectors.
        #
        # The load *axis* is the section's own chord (the principal in-plane
        # direction of the strip outline).  The absolute *sense* has no
        # source in this tree, so it is taken from the caller's configured
        # directions.  That sense is resolved **once for the whole blade**:
        # flipping each strip independently into a configured half-plane is
        # discontinuous wherever the twisting section rotates through the
        # reference direction, and the resulting 180 deg jump destroys
        # global force conservation - a uniform ``Np`` would no longer
        # integrate to ``sum_k Np_k dr_k`` (the guard
        # ``tests/test_force_projection_load_frame.py`` pins that).  A single
        # global flip merely mirrors the whole load, which changes no
        # magnitude, and leaves every station's normal parallel to its
        # neighbours.
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
        self._strip_normal_dirs: list[np.ndarray] = []
        self._strip_ac_offsets: list[np.ndarray] = []

        # Pass 1: the per-strip chord axes, made continuous along the span
        # (the SVD leaves a +/- 180 deg ambiguity per strip).
        raw_chords: list[np.ndarray | None] = []
        prev_chord: np.ndarray | None = None
        for strip in self._strips:
            chord_axis = self._strip_chord_axis(coords, strip, span_dir)
            if chord_axis is not None and prev_chord is not None:
                if float(chord_axis @ prev_chord) < 0.0:
                    chord_axis = -chord_axis
            if chord_axis is not None:
                prev_chord = chord_axis
            raw_chords.append(chord_axis)

        # Resolve the load senses once from the configured directions.
        weight = np.array([strip.dr for strip in self._strips])
        chord_sum = np.zeros(3)
        normal_sum = np.zeros(3)
        for k, axis in enumerate(raw_chords):
            if axis is None:
                continue
            chord_sum += weight[k] * axis
            normal_sum += weight[k] * np.cross(axis, span_dir)
        chord_sign = 1.0 if float(chord_sum @ self._tangential_dir) >= 0.0 else -1.0
        normal_sign = 1.0 if float(normal_sum @ self._normal_dir) >= 0.0 else -1.0

        # Pass 2: identify the ends on the resolved axis (so the AC datum
        # shares the frame's orientation), build the frame and keep the arm.
        for k, strip in enumerate(self._strips):
            station = blade_aero.stations[k]
            axis = raw_chords[k]

            if axis is None:
                # Not enough nodes for PCA — fall back to the configured
                # sense, zero AC offset.
                c_hat, n_hat = self._load_frame(self._tangential_dir)
                self._strip_chord_dirs.append(c_hat)
                self._strip_normal_dirs.append(n_hat)
                self._strip_ac_offsets.append(np.zeros(3))
                continue

            strip_pts = coords[strip.node_indices]
            chord_hat = chord_sign * axis
            le_i, te_i = self._section_ends(strip_pts, chord_hat, span_dir)
            normal_hat = normal_sign * np.cross(chord_hat, span_dir)
            n_norm = float(np.linalg.norm(normal_hat))
            if n_norm > 1e-12:
                normal_hat = normal_hat / n_norm
            else:
                normal_hat = self._normal_dir.copy()
            self._check_frame(chord_hat, normal_hat)
            self._strip_chord_dirs.append(chord_hat)
            self._strip_normal_dirs.append(normal_hat)

            ac_frac = station.airfoil.aerodynamic_center
            ac_point = strip_pts[le_i] + ac_frac * (strip_pts[te_i] - strip_pts[le_i])
            # Full 3-D vector from the aerodynamic centre to the strip
            # centroid.  Keeping the out-of-chord component (the span and
            # thickness offset) is required to reproduce the applied moment
            # about the origin; the chordwise part alone leaves a ~1 %
            # lever-arm error on the real blade.
            self._strip_ac_offsets.append(strip.centroid - ac_point)

    @staticmethod
    def _strip_chord_axis(
        coords: np.ndarray,
        strip: _Strip,
        span_dir: np.ndarray,
    ) -> np.ndarray | None:
        """Unit principal in-plane axis of a strip, or *None* if degenerate."""
        idx = strip.node_indices
        if len(idx) < 2:
            return None
        strip_pts = coords[idx]
        off = strip_pts - strip_pts.mean(axis=0)
        off_plane = off - np.outer(off @ span_dir, span_dir)
        if np.linalg.norm(off_plane) < 1e-12:
            return None
        _, _, Vt = np.linalg.svd(off_plane, full_matrices=False)
        axis = Vt[0]
        norm = float(np.linalg.norm(axis))
        return axis / norm if norm > 1e-12 else None

    @staticmethod
    def _section_ends(
        strip_pts: np.ndarray,
        chord_dir: np.ndarray,
        span_dir: np.ndarray,
    ) -> tuple[int, int]:
        """Leading/trailing node indices from the section's own geometry.

        The blunt end (the larger in-plane spread inside the outer quarter
        of the chord) is the leading edge, the sharp end the trailing edge.
        Inferring the LE from ``min(chord_proj)`` would pick whichever end
        the unrelated reference axis happened to point away from; this rule
        is convention-free and sign-free.  Within 10 % the blunt/sharp split
        is ill-conditioned (e.g. a circular root section) and the
        ``+chord_dir`` end is kept - a symmetric section is where the two
        candidate AC placements are least distinguishable.
        """
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
            if t_lo > t_hi and (t_lo - t_hi) > 0.10 * max(t_lo, t_hi):
                le_at_hi = False
        i_hi = int(np.argmax(chord_proj))
        i_lo = int(np.argmin(chord_proj))
        return (i_hi, i_lo) if le_at_hi else (i_lo, i_hi)

    def _check_frame(self, chord_hat: np.ndarray, normal_hat: np.ndarray) -> None:
        """Validate a per-strip frame is orthonormal and in the section plane."""
        tol = 1e-9
        residuals = (
            ("load frame is not orthogonal", abs(float(normal_hat @ chord_hat))),
            ("chord_hat is not unit", abs(float(np.linalg.norm(chord_hat)) - 1.0)),
            ("normal_hat is not unit", abs(float(np.linalg.norm(normal_hat)) - 1.0)),
            ("chord_hat leaves the section plane", abs(float(chord_hat @ self._span_dir))),
            ("normal_hat leaves the section plane", abs(float(normal_hat @ self._span_dir))),
        )
        for message, residual in residuals:
            if residual >= tol:
                raise ValueError(f"{message}: residual {residual:.3e} >= {tol}")

    def _load_frame(self, chord_axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Fallback frame for a strip with no usable outline geometry.

        There is no section to follow, so the configured axis itself is
        used (projected into the plane normal to the span) and the
        configured normal direction sets its sense.
        """
        c = np.asarray(chord_axis, dtype=float)
        c = c - float(c @ self._span_dir) * self._span_dir
        if np.linalg.norm(c) < 1e-12:
            # Configured tangential axis parallel to the span: any in-plane
            # axis perpendicular to the configured normal will do.
            c = np.cross(self._span_dir, self._normal_dir)
        if np.linalg.norm(c) < 1e-12:
            helper = np.array([1.0, 0.0, 0.0])
            if abs(float(helper @ self._span_dir)) > 0.9:
                helper = np.array([0.0, 1.0, 0.0])
            c = np.cross(self._span_dir, helper)
        c = c / np.linalg.norm(c)

        n = np.cross(c, self._span_dir)
        n = n / np.linalg.norm(n)
        if float(n @ self._normal_dir) < 0.0:
            n = -n
        self._check_frame(c, n)
        return c, n

    # ------------------------------------------------------------------
    #  Public API
    # ------------------------------------------------------------------

    def project(self, bem_result: BEMResult) -> np.ndarray:
        """Map BEM loads onto mesh nodes.

        Each strip's load frame is derived from that strip's own section
        geometry (:meth:`_load_frame`); the configured *normal_direction*
        and *tangential_direction* only select the sense of the two axes,
        they are no longer a fixed global direction.

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

            # Strip force on its own section axes (axis from geometry, sense
            # from the configured normal/tangential directions).
            F_strip = F_n * self._strip_normal_dirs[k] + F_t * self._strip_chord_dirs[k]

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
        """Check global force conservation across the discretised strips.

        The expected force is recomputed from the BEM loads and the *same*
        per-strip section frames :meth:`project` uses,
        ``sum_k (Np_k * normal_hat_k + Tp_k * chord_hat_k) * dr_k``.  The
        vector sum is therefore no longer parallel to any single configured
        axis: the direction is per strip.  Because the frame is shared with
        :meth:`project`, ``force_error`` measures only the fidelity of the
        per-strip nodal distribution (for example a dropped or clamped
        strip); it does **not** independently validate the section frame or
        the load sense - that is what
        ``tests/test_force_projection_load_frame.py`` does.

        Returns
        -------
        dict
            ``force_error`` – norm of residual (N).
            ``force_bem`` – expected total BEM force vector (N).
            ``force_mesh`` – actual summed mesh force vector (N).
        """
        # Discrete strip sum with the same per-strip frames as project().
        F_bem = np.zeros(3)
        for k, strip in enumerate(self._strips):
            if len(strip.node_indices) == 0:
                continue
            F_bem += (
                float(bem_result.Np[k]) * strip.dr * self._strip_normal_dirs[k]
                + float(bem_result.Tp[k]) * strip.dr * self._strip_chord_dirs[k]
            )

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
