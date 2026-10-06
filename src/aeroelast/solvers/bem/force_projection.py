"""Conservative force projection from BEM stations onto 3D shell mesh.

Maps per-span-station aerodynamic loads (Np, Tp) from a BEM computation
onto the finite-element mesh nodes, preserving the total integrated force
and moment on each chordwise strip.

The section frame and the aerodynamic-centre datum are taken from each strip's
**station ring** - the physical-ring group nearest ``strip.r_center`` - not from
the whole BEM band.  A band is ~2.4 m wide on the IEA-15MW blade and holds three
physical rings of a tapered, twisted, prebent blade, so a band-wide chord is not
any section's chord: measured on this tree's own blade, 12 of the 50 strips chose
their LE and TE on **different** physical rings (24%), and against the independent
AeroDyn chord the band datum is off by a median 1.517% / p90 5.538% (max 174.6%
at the tip) while the station-ring datum is off by a median 0.648% / p90 1.736%.

Pitching-moment sense
---------------------
The BEM polars define ``Cm`` (and therefore the per-unit-span ``Mp``) about the
aerofoil chord **from the leading edge to the trailing edge**, positive nose-up:
the aerodynamic moment vector rides ``chord_le_to_te x n_hat``.  The projector
resolves ``_strip_chord_dirs`` from the tangential-force sense, and that axis is
not guaranteed to run leading-to-trailing - measured against the deck's own
WindIO aerofoils (``tools/diagnose_sign_chain.py``) it runs *trailing-to-leading*
on **every real station** of the IEA-15MW blade (``LE . c_hat = -1.0``).  The
per-strip handedness is read back from :meth:`_section_ends` and stored as
``_strip_moment_axis_sign``: with the measured handedness ``c_hat x n_hat =
+span_hat``, a leading edge at the HIGH end of the chord projection means the
aerodynamic moment axis is ``-span_hat``.  Applying ``Mp`` on the raw
``+span_hat`` (the pre-fix behaviour) inverted the term: the requested section
torque came out NOSE-UP at ``+4.839888e+05 N.m``, the lever-arm transfer
``+7.607e+05`` dominating and opposing the polars' ``-2.766e+05`` (measured by
``tools/diagnose_applied_torsion_sign.py``).  The force terms and the lever-arm
transfer are frame-consistent and are left untouched; only the aerodynamic moment
term is applied on the axis that :meth:`_section_ends` selects.
"""

import logging
from dataclasses import dataclass
from typing import List

import numpy as np
from scipy.linalg import lstsq

from aeroelast.core.mesh.entities import ElementType
from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.aerodynamics import BladeAero
from aeroelast.solvers.bem.engine import BEMResult
from aeroelast.solvers.bem.section_contour import (
    SectionCell,
    cell_adjacency,
    contour_report,
    rings_in_band,
    section_cells,
)

logger = logging.getLogger(__name__)

#: Fraction of the **mesh's** span extent used as the ring-merging gap tolerance.
#: A strip is one BEM band (a few metres) and holds several physical mesh rings, each
#: of which is prebent: its nodes spread over ~1e-3 m of span instead of lying in one
#: plane.  Judging the ring gap against the strip's own extent puts the tolerance below
#: that spread and shears each physical ring into partial arcs - and an arc is still a
#: properly ordered, ``contour_report.usable`` polygon, so the gate below cannot catch
#: it and Bredt's flow is applied to a bogus, small-area "ring".  The tolerance must
#: therefore come from the mesh, whose extent it shares with the physical stations: on
#: the 117 m blade this is 0.0117 m, above the ~1e-3 m prebend spread and well below the
#: ~0.44 m station spacing.
RING_GAP_FRACTION = 1e-4


@dataclass
class _Strip:
    """Internal: a chordwise strip of mesh nodes."""

    node_indices: np.ndarray  # indices into mesh.nodes
    r_center: float  # span-wise centre (m, from hub)
    dr: float  # strip width (m)
    centroid: np.ndarray  # 3-D centroid of strip nodes
    offsets: np.ndarray  # (n, 3) node positions relative to centroid


@dataclass
class _RingSection:
    """One physical ring's wall graph, bound to the owning strip's local node space.

    ``cells`` / ``adjacency`` / ``edge_length`` / ``edge_shear_stiffness`` are the
    :func:`ring_section` output, whose node indices are **local to the ring**
    (``0 .. len(local_nodes) - 1``).  ``local_nodes`` maps each of those ring-local
    indices to the owning strip's own ``offsets`` index (a position in
    ``strip.node_indices``), so the multi-cell realisation can add a wall-flow half to
    ``f_shear[local_nodes[a]]``.  The two index spaces differ - a ring group is a
    permutation/subset of the strip's nodes - so the mapping is carried explicitly
    rather than assumed to be the identity.

    ``from_element_properties`` records whether the ring's wall stiffnesses came from
    the projector's ``element_properties`` map.  Without it :func:`ring_section`
    returns the uniform, geometric-only ``S = 1.0``; the multi-cell gate must refuse
    that case explicitly rather than infer it from a uniform stiffness.
    """

    cells: List[SectionCell]
    adjacency: dict
    edge_length: dict
    edge_shear_stiffness: dict
    local_nodes: np.ndarray
    from_element_properties: bool


def _cells_have_usable_area(
    cells: List[SectionCell], offsets: np.ndarray, span_dir: np.ndarray
) -> bool:
    """Are all of a ring's extracted cells large enough to solve on?

    :func:`section_cells` keeps every bounded face with a strictly positive shoelace
    area, but a projected drawing of a **deformed** wall graph can trace a
    self-overlapping loop whose area is a round-off residual (~1e-16 m^2).  Such a
    face carries a near-singular compliance and the multi-cell solve returns
    meaningless flows, so the ring is unusable.  The floor is a fraction of the
    ring's own in-plane size (``1e-12 * radius^2``, *radius* the largest in-plane
    distance from the ring's own mean), which is scale-free between meshes and far
    below any real wall cell.  A ring rejected here is recorded as a failed ring and
    its strip falls back to the minimum-norm realisation.
    """
    if not cells:
        return False
    pts = np.asarray(offsets, dtype=float)
    if pts.size == 0:
        return False
    rel = pts - pts.mean(axis=0)
    in_plane = rel - np.outer(rel @ span_dir, span_dir)
    radius = float(np.sqrt(np.max(np.sum(in_plane**2, axis=1))))
    floor = 1e-12 * radius**2
    return all(float(cell.area) > floor for cell in cells)


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
        element_properties: dict | None = None,
    ):
        span_dir = np.asarray(
            span_direction if span_direction is not None else [0.0, 0.0, 1.0],
            dtype=float,
        )
        span_dir /= np.linalg.norm(span_dir)
        self._span_dir = span_dir

        self._mesh = mesh
        #: ``{element_set_name: property}`` (the assembler's own map).  Fed to
        #: :func:`ring_section` so each wall's multi-cell flow sees its laminate's
        #: membrane shear stiffness ``S = G*t``.  When it is ``None``, or a ring's
        #: element set has no entry, every wall falls back to the uniform ``S = 1.0``
        #: and the cell-to-cell split is geometric-only; such a ring is marked
        #: non-physical and ``project()`` routes it to :meth:`_distribute` instead of
        #: claiming a physical split.
        self._element_properties = element_properties

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

        # Ring-merging tolerance from the mesh's global span extent, never a
        # strip's (see RING_GAP_FRACTION).  The same value separates every strip's
        # physical rings.
        self._ring_gap_tolerance = RING_GAP_FRACTION * float(span_coords.max() - span_coords.min())

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

        # Group every strip's nodes into physical rings **once**, with the
        # mesh-scaled tolerance.  ``project()`` hands these groups to
        # :func:`realise_section_load` so the strip's own (band-wide) extent can
        # never set a tolerance below its rings' prebend spread.
        self._strip_ring_groups: list[list[np.ndarray]] = []
        for strip in self._strips:
            if len(strip.node_indices) == 0:
                self._strip_ring_groups.append([])
                continue
            strip_span = (strip.centroid + strip.offsets) @ span_dir
            self._strip_ring_groups.append(rings_in_band(strip_span, self._ring_gap_tolerance))

        # Precompute each physical ring's multi-cell wall graph once, from the
        # mesh, so ``project()`` realises a strip's span moment through
        # :func:`multi_cell_shear_flow` instead of the minimum-norm field.  The
        # ring group holds **strip-local** offsets indices; :func:`ring_section`
        # wants **global** ``mesh.nodes`` indices, and its own cells index the ring
        # locally.  ``_RingSection.local_nodes`` keeps that local-vs-global mapping
        # explicit.  A ring that has fewer than three nodes, that bounds no closed
        # cell on this mesh (a degenerate/collinear outline, or a coupling mesh
        # stripped of its elements), or whose cells carry only a round-off-sized
        # area (see :func:`_cells_have_usable_area`) is recorded as ``None`` rather
        # than raising: :func:`realise_section_load` then falls back to
        # :meth:`_distribute` for the whole strip, exactly as before.
        self._strip_ring_sections: list[list[_RingSection | None]] = []
        for k, strip in enumerate(self._strips):
            if len(strip.node_indices) == 0:
                self._strip_ring_sections.append([])
                continue
            sections: list[_RingSection | None] = []
            for group in self._strip_ring_groups[k]:
                if len(group) < 3:
                    sections.append(None)
                    continue
                local = np.asarray(group, dtype=np.intp)
                try:
                    cells, adjacency, edge_length, edge_shear_stiffness = ring_section(
                        mesh, strip.node_indices[local], span_dir, element_properties
                    )
                except ValueError:
                    sections.append(None)
                    continue
                if not _cells_have_usable_area(cells, strip.offsets[local], span_dir):
                    # A face with a round-off-sized area (a self-overlapping loop
                    # the projected, deformed wall graph can trace) makes the
                    # multi-cell solve ill-conditioned and its flows meaningless.  Such
                    # a ring is a failure to extract a usable section, not a section.
                    sections.append(None)
                    continue
                sections.append(
                    _RingSection(
                        cells,
                        adjacency,
                        edge_length,
                        edge_shear_stiffness,
                        local,
                        from_element_properties=element_properties is not None,
                    )
                )
            self._strip_ring_sections.append(sections)

        # ------------------------------------------------------------------
        # Per-strip load frames and AC-to-centroid offset vectors.
        #
        # The load *axis* is the section's own chord: the principal in-plane
        # direction of the strip's **station ring** (the physical-ring group
        # nearest ``strip.r_center``), not of the whole band.  A band holds
        # several physical rings of a tapered, twisted, prebent blade - on
        # this tree's IEA-15MW blade 12 of the 50 strips pick their LE and TE
        # on different rings (24%), so a band-wide chord is not any section's
        # chord.  Against the independent AeroDyn chord the band datum is off
        # by a median 1.517% / p90 5.538% (max 174.6% at the tip); the
        # station-ring datum by a median 0.648% / p90 1.736%.  The absolute
        # *sense* has no source in this tree, so it is taken from the caller's
        # configured directions.  That sense is resolved **once for the whole
        # blade**:
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
        #: Per strip, the factor that carries the polars' nose-up ``Mp`` onto this
        #: frame's span axis: ``+1.0`` when leading-to-trailing runs along
        #: ``_strip_chord_dirs[k]``, ``-1.0`` when it runs against it.  See the module
        #: docstring's "Pitching-moment sense" and ``project()``.
        self._strip_moment_axis_sign: list[float] = []

        # Pass 1: the per-strip chord axes, made continuous along the span
        # (the SVD leaves a +/- 180 deg ambiguity per strip).  The axis of
        # strip k is measured on its station-ring points, the same points the
        # AC datum below uses, so the frame and the datum share one section.
        raw_chords: list[np.ndarray | None] = []
        prev_chord: np.ndarray | None = None
        for k, strip in enumerate(self._strips):
            ring_pts = self._station_ring_points(k, strip, span_dir)
            chord_axis = self._points_chord_axis(ring_pts, span_dir)
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
                # No section outline, hence no measured handedness to convert: keep
                # the pre-existing raw span axis rather than invent a flip.
                self._strip_moment_axis_sign.append(1.0)
                continue

            # The ends and the AC live on the strip's station ring, so the
            # chord the fraction is measured along is one physical section's
            # chord rather than a line between two different rings.  The arm
            # below still runs to the strip centroid, which is the point
            # ``_distribute`` balances moments about.
            ring_pts = self._station_ring_points(k, strip, span_dir)
            chord_hat = chord_sign * axis
            le_i, te_i = self._section_ends(ring_pts, chord_hat, span_dir)
            # The polars give Cm about the aerofoil chord ``chord_le_to_te`` (positive
            # nose-up).  ``_section_ends`` is sign-free and locates the leading edge;
            # comparing it with the HIGH end of the ring's projection on ``chord_hat``
            # reads the handedness of ``_strip_chord_dirs``.  ``le_at_hi`` means
            # leading-to-trailing runs AGAINST ``chord_hat``, so the aerodynamic moment
            # axis - normal to the chord, i.e. ``chord_le_to_te x n_hat`` - is
            # ``-(chord_hat x n_hat) = -span_hat`` (the measured handedness on every
            # real station of this blade; tools/diagnose_sign_chain.py).
            le_at_hi = le_i == int(np.argmax(ring_pts @ chord_hat))
            self._strip_moment_axis_sign.append(-1.0 if le_at_hi else 1.0)
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
            ac_point = ring_pts[le_i] + ac_frac * (ring_pts[te_i] - ring_pts[le_i])
            # Full 3-D vector from the aerodynamic centre to the strip
            # centroid.  Keeping the out-of-chord component (the span and
            # thickness offset) is required to reproduce the applied moment
            # about the origin; the chordwise part alone leaves a ~1 %
            # lever-arm error on the real blade.
            self._strip_ac_offsets.append(strip.centroid - ac_point)

    @staticmethod
    def _points_chord_axis(
        points: np.ndarray,
        span_dir: np.ndarray,
    ) -> np.ndarray | None:
        """Unit principal in-plane axis of a point set, or *None* if degenerate."""
        pts = np.asarray(points, dtype=float)
        if len(pts) < 2:
            return None
        off = pts - pts.mean(axis=0)
        off_plane = off - np.outer(off @ span_dir, span_dir)
        if np.linalg.norm(off_plane) < 1e-12:
            return None
        _, _, Vt = np.linalg.svd(off_plane, full_matrices=False)
        axis = Vt[0]
        norm = float(np.linalg.norm(axis))
        return axis / norm if norm > 1e-12 else None

    @classmethod
    def _strip_chord_axis(
        cls,
        coords: np.ndarray,
        strip: _Strip,
        span_dir: np.ndarray,
    ) -> np.ndarray | None:
        """Unit principal in-plane axis of a strip's nodes, or *None* if degenerate.

        Kept for existing callers; :meth:`_points_chord_axis` holds the body.
        """
        return cls._points_chord_axis(coords[strip.node_indices], span_dir)

    def _station_ring_points(
        self,
        k: int,
        strip: _Strip,
        span_dir: np.ndarray,
    ) -> np.ndarray:
        """The strip's station ring: the physical-ring group nearest ``r_center``.

        The aerodynamic centre is defined on a *station's section*, so it has to
        be read from one physical ring, not from the whole BEM band.  Of the
        groups in ``self._strip_ring_groups[k]`` pick the one whose mean span
        coordinate is nearest ``strip.r_center``, skipping any group with fewer
        than 3 nodes (a ring needs a plane outline for the chord axis and the
        blunt-end rule).  A strip with no such ring - an empty or a one-node
        strip - falls back to its whole node set, the pre-existing behaviour, so
        the caller keeps a defined datum.
        """
        ring_pts = strip.centroid + strip.offsets
        groups = self._strip_ring_groups[k]
        best: np.ndarray | None = None
        best_distance = np.inf
        for group in groups:
            if len(group) < 3:
                continue
            radius = float(np.mean(ring_pts[group] @ span_dir))
            distance = abs(radius - strip.r_center)
            if distance < best_distance:
                best_distance = distance
                best = group
        if best is None:
            return ring_pts
        return ring_pts[np.asarray(best, dtype=np.intp)]

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

        Each strip's load frame is derived from its own station ring
        (:meth:`_station_ring_points` and the section axes built on it); the
        configured *normal_direction* and *tangential_direction* only select
        the sense of the two axes, they are no longer a fixed global
        direction.

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
            # the aerodynamic centre).  The polars' Cm is positive nose-up about
            # ``chord_le_to_te x n_hat``; the raw span axis is only that chord axis
            # when leading-to-trailing runs along ``_strip_chord_dirs[k]``, so the
            # per-strip sign below carries the polars' sense onto this frame (see the
            # module docstring and ``_strip_moment_axis_sign``).  The geometric
            # transfer term accounts for the moment arm between the AC and the centroid
            # of the strip nodes and is always present regardless of Mp availability.
            M_ac = (
                float(bem_result.Mp[k])
                * strip.dr
                * (self._strip_moment_axis_sign[k] * self._span_dir)
                if bem_result.Mp is not None
                else np.zeros(3)
            )
            M_strip = M_ac - np.cross(self._strip_ac_offsets[k], F_strip)

            # Distribute to nodes.  The precomputed per-ring wall graphs carry the
            # span moment as Bredt-Batho multi-cell wall flows, one tributary share of
            # the strip's torsion per physical ring; whatever force/transverse moment
            # is left goes through the unchanged minimum-norm solve.  ``ring_groups``
            # and ``ring_sections`` are the strip's own precomputed physical rings, so
            # the multi-cell branch never regroups with a strip-scaled tolerance.  A
            # strip whose rings could not be extracted - or that was built without
            # ``element_properties``, so its walls are only the uniform ``S = 1.0`` -
            # falls back to ``_distribute`` inside ``realise_section_load``.
            f_nodes = realise_section_load(
                strip,
                F_strip,
                M_strip,
                self._span_dir,
                ring_groups=self._strip_ring_groups[k],
                ring_sections=self._strip_ring_sections[k],
            )
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
        # scipy's stubs type the return as optional although it never is for valid input;
        # narrow it explicitly instead of unpacking a possibly-None tuple.
        solution = lstsq(AAT, b)
        if solution is None:  # pragma: no cover - the stub's optional, not a real branch
            raise RuntimeError("scipy.linalg.lstsq returned no solution")
        lam = solution[0]
        f_flat = A.T @ lam

        return f_flat.reshape(n, 3)


def realise_section_load(
    strip: _Strip,
    F_strip: np.ndarray,
    M_strip: np.ndarray,
    span_dir,
    ring_groups: list[np.ndarray] | None = None,
    ring_sections: list[_RingSection | None] | None = None,
) -> np.ndarray:
    """Realise a strip load with its span moment carried by a wall shear flow.

    A closed thin-walled section under a torque ``M`` about its span axis carries the
    Bredt wall shear flow ``q = M / (2A)`` along its wall mid-line.  The constrained
    minimum-norm solve :meth:`ForceProjector._distribute` instead spreads a pure span moment
    as the circle-tangential field ``f_j = omega x d_j``, which is not a Saint-Venant wall
    traction.  This function realises the span component of ``M_strip`` as wall flows and
    hands whatever is left to :meth:`ForceProjector._distribute`.

    There are two wall-flow modes, selected by ``ring_sections``:

    * **Multi-cell (production):** when ``ring_sections`` is given, every physical ring of
      the strip is realised with the Bredt-Batho **multi-cell** flow of its own wall graph
      (:func:`multi_cell_shear_flow`), so an internal web shared by two cells carries
      ``q_i - q_j``.  This is what :meth:`ForceProjector.project` passes.
    * **Single-cell (direct callers):** when ``ring_sections`` is ``None``, the pre-existing
      single-ring :func:`contour_report` / Bredt ``q = torsion / (2 * area)`` route below is
      used unchanged - byte-for-byte.  This is what the thin-walled-tube test
      (``tests/validation/parity/test_thin_walled_tube_moment_realization.py``) calls
      directly, and it is all a caller without a mesh can get.

    Steps, in order:

    1. **Split the moment.**  ``span_hat = span_dir / ||span_dir||``;
       ``torsion = M_strip . span_hat`` is the part the wall flow can carry, and
       ``transverse = M_strip - torsion * span_hat`` is left to the fallback.
    2. **Fail fast on an unextractable or non-physical ring.**  In the multi-cell mode,
       ``ring_sections`` holds one entry per non-empty physical-ring group of the strip, in
       the strip's own ring order.  A ring is only realisable when it carries **real** wall
       stiffnesses (``_RingSection.from_element_properties``), bounds **at least one**
       closed cell, and every wall's ``edge_shear_stiffness`` is strictly positive.  If
       **any** entry fails any of those - an entry is ``None`` because the ring has fewer
       than three nodes, because :func:`ring_section` found no closed cell for it (a
       collinear/degenerate outline, or a coupling mesh stripped of its elements), or
       because a cell's shoelace area is a round-off residual rather than a wall cell
       (:func:`_cells_have_usable_area`); or the whole projector was built without
       ``element_properties``, so every wall fell back to the uniform, geometric-only
       ``S = 1.0`` - the strip falls back to :meth:`ForceProjector._distribute`; a
       half-realised or geometric-only strip is never returned.
    3. **Tributary per-ring flow.**  Otherwise the strip's torsion is split over its rings
       by their **tributary span weight**: sort the rings by their own mean span position
       ``t_0 < ... < t_{n-1}`` inside the strip and give ring ``i`` the trapezoidal
       length ``w_0 = (t_1 - t_0)/2``, ``w_i = (t_{i+1} - t_{i-1})/2``,
       ``w_{n-1} = (t_{n-1} - t_{n-2})/2``, normalised by ``sum w`` (a single ring takes
       all of ``torsion``; see :func:`_ring_tributary_weights`).  The trapezoid is the
       right quadrature here because the BEM ``Mp`` is a **per-unit-span** moment: a strip
       is one BEM band wider than one physical ring, so the band's moment per metre has to
       be integrated over each ring's own share of the span - an equal share would put the
       same moment on a short root ring and a long tip ring.  Each ring is then solved with
       :func:`multi_cell_shear_flow` on its own cells, adjacency, edge lengths and per-wall
       shear stiffness ``S = G*t``.
    4. **Per-cell walk.**  For every cell and every wall ``(a, b)`` of its boundary, in the
       cell's own (already ``+span``-CCW-normalised) traversal order, the load
       ``0.5 * q_cell * ell * tangent_hat`` is added to **both** endpoints.  A wall shared by
       cells ``i`` and ``j`` is traversed in opposite directions by the two cells, so it
       automatically receives ``q_i - q_j`` - this is exactly the wall flow the multi-cell
       system imposes, and it is why the per-cell walk is the right realisation and no
       shared-wall special case is needed.  Summing the cell moments, the realised span
       moment is ``sum_i 2 A_i q_i = torsion`` by the solver's own torque equilibrium.
    5. **Residual.**  The realised nodal moment ``M_shear = sum_j cross(offset_j, f_j)`` is
       read back about the strip centroid and the leftover force/moment
       ``F_strip - sum_j f_j`` and ``M_strip - M_shear`` are routed through
       :meth:`ForceProjector._distribute`.  The shear flow already delivers the span moment,
       so this keeps force and moment conservation exact (to the pseudoinverse) regardless
       of round-off in the shear-flow algebra.  The residual step and the signed span-moment
       guard below apply to the multi-cell mode as well: the realised span moment must still
       match the requested ``torsion`` or the call raises.

    In the single-cell mode the nodes are grouped as before: ``ring_groups`` (int index
    arrays into the strip's own ``offsets``) when given, otherwise :func:`rings_in_band`
    with ``1e-4 * (strip span extent)`` - which **assumes the strip is a single ring**, as it
    is for a tube cross-section.  A strip can hold several physical rings because the mesh's
    BEM band is wider than one ring; on a prebent blade one ring spreads over ~1e-3 m of span
    while stations are ~0.44 m apart, and a strip's own extent sets a tolerance below that
    spread, shearing the ring into arcs.  The single-cell mode is gated to **exactly one
    usable** group (the group is ordered and validated with :func:`contour_report`); anything
    else falls back to :meth:`ForceProjector._distribute`.

    Before adding the residual, the delivered span moment is checked against the requested one
    (when ``torsion != 0``).  The check is the **signed** difference
    ``abs(span_hat @ M_shear - torsion)``: a difference of magnitudes would accept an
    orientation bug that flips the delivered moment, and on the tube the signed difference is
    pure round-off because the flow's moment ``2 q A`` and the area ``A`` come from the same
    ordered polygon.  ``1e-9`` relative is therefore generous for the tube and hard for a sign
    or orientation error, which fails loudly instead of being absorbed by the residual.

    Pure over ``strip``: the ``_Strip`` and its arrays are read-only here.

    Raises
    ------
    ValueError
        If ``span_dir`` is zero-length, or the wall flow's realised span moment disagrees with
        the requested one by more than 1e-9 relative.
    """
    span_hat = np.asarray(span_dir, dtype=float).reshape(3)
    span_norm = float(np.linalg.norm(span_hat))
    if span_norm == 0.0:
        raise ValueError("span_dir must be non-zero")
    span_hat = span_hat / span_norm

    M_strip = np.asarray(M_strip, dtype=float)
    F_strip = np.asarray(F_strip, dtype=float)
    torsion = float(M_strip @ span_hat)

    if ring_sections is not None:
        return _realise_multi_cell_section_load(
            strip, F_strip, M_strip, span_hat, torsion, list(ring_sections)
        )

    if ring_groups is None:
        # Fallback for direct callers: the strip's own extent.  This assumes the strip is
        # one physical ring (true for a tube cross-section).  ``project()`` never takes this
        # path - it always passes the mesh-scaled groups precomputed in __init__.
        points = strip.centroid + strip.offsets
        span_coords = points @ span_hat
        if span_coords.size:
            span_extent = float(span_coords.max() - span_coords.min())
        else:
            span_extent = 0.0
        ring_groups = rings_in_band(span_coords, RING_GAP_FRACTION * span_extent)

    groups = [np.asarray(g, dtype=np.intp).ravel() for g in ring_groups if len(g) > 0]
    if len(groups) != 1:
        # Several physical rings in one BEM band.  The multi-cell realisation needs each
        # ring's own wall graph, and only ``ring_sections`` carries it; without it there is
        # no per-ring system to solve, so fall back rather than ship an unvalidated
        # single-cell flow on one arbitrary ring.
        return ForceProjector._distribute(strip, F_strip, M_strip)

    group = groups[0]
    report = contour_report(strip.centroid + strip.offsets[group], span_hat)
    if not report.usable:
        return ForceProjector._distribute(strip, F_strip, M_strip)

    n = len(strip.node_indices)
    f_shear = np.zeros((n, 3))
    ordered = strip.centroid + strip.offsets[group][report.order]
    local = group[report.order]
    q = torsion / (2.0 * report.area)
    m = len(local)
    for k in range(m):
        i = int(local[k])
        j = int(local[(k + 1) % m])
        edge = ordered[(k + 1) % m] - ordered[k]
        ell = float(np.linalg.norm(edge))
        if ell == 0.0:  # unusable rings are already excluded; keep the guard local
            continue
        half = 0.5 * (q * ell) * (edge / ell)
        f_shear[i] = f_shear[i] + half
        f_shear[j] = f_shear[j] + half

    M_shear = np.cross(strip.offsets, f_shear).sum(axis=0)
    if torsion != 0.0:
        delivered = float(M_shear @ span_hat)
        if abs(delivered - torsion) > 1e-9 * abs(torsion):
            raise ValueError(
                "wall shear flow realised a span moment of "
                f"{delivered:.6e} N.m but {torsion:.6e} N.m was requested: the ring order or "
                "the span axis is wrong (this is not round-off)"
            )

    residual_force = F_strip - f_shear.sum(axis=0)
    residual_moment = M_strip - M_shear
    return f_shear + ForceProjector._distribute(strip, residual_force, residual_moment)


def _ring_section_is_realisable(section: _RingSection) -> bool:
    """Can this ring's multi-cell system be solved on **real** wall stiffnesses?

    Three conditions, all required: the ring was built from ``element_properties``
    (:func:`ring_section` returns a uniform, geometric-only ``S = 1.0`` without them, so
    a no-properties projector must not claim a physical split), it bounds at least one
    closed cell, and every wall's ``S = G*t`` is strictly positive (a zero-stiffness wall
    makes the compliance singular).  Any failure routes the whole strip to
    :meth:`ForceProjector._distribute`; a half-realised or geometric-only strip is never
    returned.
    """
    if not section.from_element_properties:
        return False
    if not section.cells:
        return False
    return all(float(value) > 0.0 for value in section.edge_shear_stiffness.values())


def _ring_tributary_weights(positions: np.ndarray) -> np.ndarray:
    """Normalised trapezoidal tributary weights of rings at span ``positions``.

    A strip is one BEM band wider than one physical ring and the BEM ``Mp`` is a
    **per-unit-span** moment, so each ring carries the moment integrated over its own
    share of the span.  The shares are the composite trapezoid lengths of the rings' own
    span positions - ``w_0 = (t_1 - t_0)/2``, ``w_i = (t_{i+1} - t_{i-1})/2``,
    ``w_{n-1} = (t_{n-1} - t_{n-2})/2`` - normalised so they sum to one.  A single ring
    takes all of the torsion.  If the positions are coincident (a degenerate grouping)
    the weights fall back to an equal share so the row stays solvable rather than
    dividing by zero.
    """
    values = np.asarray(positions, dtype=float).ravel()
    count = values.size
    if count == 0:
        return np.zeros(0)
    if count == 1:
        return np.ones(1)
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    sorted_weights = np.empty(count, dtype=float)
    sorted_weights[0] = 0.5 * (sorted_values[1] - sorted_values[0])
    sorted_weights[-1] = 0.5 * (sorted_values[-1] - sorted_values[-2])
    sorted_weights[1:-1] = 0.5 * (sorted_values[2:] - sorted_values[:-2])
    total = float(sorted_weights.sum())
    if total > 0.0:
        sorted_weights = sorted_weights / total
    else:
        sorted_weights = np.full(count, 1.0 / count)
    weights = np.empty(count, dtype=float)
    weights[order] = sorted_weights
    return weights


def _realise_multi_cell_section_load(
    strip: _Strip,
    F_strip: np.ndarray,
    M_strip: np.ndarray,
    span_hat: np.ndarray,
    torsion: float,
    ring_sections: list[_RingSection | None],
) -> np.ndarray:
    """Multi-cell branch of :func:`realise_section_load`; see its docstring for the design.

    ``ring_sections`` is one entry per non-empty physical-ring group of the strip.  A
    ring is only realisable when it was built from ``element_properties`` (real wall
    stiffnesses), bounds at least one cell, and every wall stiffness is positive; any
    ``None`` entry, empty list, no-cell ring or non-positive wall fails the gate
    (:func:`_ring_section_is_realisable`) and the whole strip falls back to the unchanged
    :meth:`ForceProjector._distribute`, so a partially-realised or geometric-only strip
    is never returned.  Otherwise the torsion is split over the rings by their tributary
    span weight (:func:`_ring_tributary_weights`), ``multi_cell_shear_flow`` solves each
    ring's cells, and the per-cell wall walk below (no shared-wall special case: the two
    cells traverse a shared wall in opposite directions, so it collects ``q_i - q_j`` by
    construction) realises the field.
    """
    if not ring_sections:
        return ForceProjector._distribute(strip, F_strip, M_strip)
    sections: list[_RingSection] = []
    for section in ring_sections:
        if section is None or not _ring_section_is_realisable(section):
            return ForceProjector._distribute(strip, F_strip, M_strip)
        sections.append(section)

    n = len(strip.node_indices)
    f_shear = np.zeros((n, 3))
    # Tributary weight from each ring's own mean span position inside the strip.
    positions = np.array(
        [
            float((strip.centroid + strip.offsets[section.local_nodes]).mean(axis=0) @ span_hat)
            for section in sections
        ]
    )
    weights = _ring_tributary_weights(positions)
    for section, weight in zip(sections, weights, strict=True):
        ring_torsion = torsion * float(weight)
        flows, _theta_rate = multi_cell_shear_flow(
            section.cells,
            section.adjacency,
            section.edge_length,
            section.edge_shear_stiffness,
            ring_torsion,
        )
        local_nodes = section.local_nodes
        for cell, flow in zip(section.cells, flows, strict=True):
            for a, b in cell.edges:
                ia = int(local_nodes[int(a)])
                ib = int(local_nodes[int(b)])
                edge = strip.offsets[ib] - strip.offsets[ia]
                ell = float(np.linalg.norm(edge))
                if ell == 0.0:
                    continue
                half = 0.5 * (flow * ell) * (edge / ell)
                f_shear[ia] = f_shear[ia] + half
                f_shear[ib] = f_shear[ib] + half

    M_shear = np.cross(strip.offsets, f_shear).sum(axis=0)
    if torsion != 0.0:
        delivered = float(M_shear @ span_hat)
        if abs(delivered - torsion) > 1e-9 * abs(torsion):
            raise ValueError(
                "multi-cell wall shear flows realised a span moment of "
                f"{delivered:.6e} N.m but {torsion:.6e} N.m was requested: the cell "
                "order or the span axis is wrong (this is not round-off)"
            )

    residual_force = F_strip - f_shear.sum(axis=0)
    residual_moment = M_strip - M_shear
    return f_shear + ForceProjector._distribute(strip, residual_force, residual_moment)


def _membrane_shear_stiffness(property_object) -> float | None:
    """Integrated membrane shear stiffness ``A66`` [N/m] of a shell property.

    The shell Voigt layout is ``[sigma_xx, sigma_yy, tau_xy]``, so index ``(2, 2)``
    of the merged membrane stiffness ``Cm = C * h`` is the wall's in-plane shear
    stiffness ``A66`` in N/m.  The property objects this tree produces are the
    Rust ``_aeroelast.Laminate`` (which exposes the classical-lamination A matrix
    as ``a_matrix()``) and the isotropic / composite dicts the assembler accepts
    (``{"type": "isotropic", ...}`` or a flat ``"cm"`` A matrix); a Python
    property exposing ``Cm()`` -- the historical element API the stress-recovery
    docstring names -- is accepted as well.  Returns ``None`` when the property is
    absent or unrecognised, so the caller can fall back to a uniform stiffness.
    """
    if property_object is None:
        return None
    cm = None
    cm_method = getattr(property_object, "Cm", None)
    if callable(cm_method):
        cm = np.asarray(cm_method(), dtype=float)
    elif hasattr(property_object, "a_matrix"):
        cm = np.asarray(property_object.a_matrix(), dtype=float)
    elif isinstance(property_object, dict):
        if "cm" in property_object:
            cm = np.asarray(property_object["cm"], dtype=float).reshape(3, 3)
        elif property_object.get("type") == "isotropic":
            # Plane-stress shear modulus G = E / (2 (1 + nu)); A66 = G * t.  The
            # shear correction factor is for transverse shear, not the membrane A66.
            e = float(property_object["e"])
            nu = float(property_object["nu"])
            thickness = float(property_object["thickness"])
            return e / (2.0 * (1.0 + nu)) * thickness
    if cm is None or cm.shape != (3, 3):
        return None
    return float(cm[2, 2])


def _element_property_by_id(mesh: MeshModel, element_properties: dict | None) -> dict[int, object]:
    """Map element id to its property, mirroring ``_to_rust_mesh``'s set convention.

    ``element_properties`` is the ``{element_set_name: property}`` mapping returned
    by ``Blade.get_element_properties()`` (a Rust ``Laminate`` or an isotropic /
    composite dict).  An element's property is reached through the very membership
    the assembler uses -- ``mesh.element_sets[set_name].elements`` -- rather than a
    second, element-id-keyed convention.  When two sets claim the same element the
    first set in ``element_properties`` order wins, deterministically.
    """
    mapping: dict[int, object] = {}
    if not element_properties:
        return mapping
    for set_name, prop in element_properties.items():
        element_set = mesh.element_sets.get(set_name)
        if element_set is None:
            continue
        for element in element_set.elements:
            mapping.setdefault(int(element.id), prop)
    return mapping


def ring_section(
    mesh: MeshModel,
    ring_nodes,
    span_dir,
    element_properties: dict | None = None,
) -> tuple[list[SectionCell], dict, dict, dict]:
    """Cells, adjacency, edge lengths and edge shear stiffness of one physical ring.

    A closed thin-walled section's torsional flow depends on the wall graph of the
    section -- its outer skin plus any internal webs -- and on each wall's membrane
    shear stiffness.  This function recovers both from the finite-element mesh for
    **one physical ring** (a strip's ring group, or any constant-span node set).

    Parameters
    ----------
    mesh : MeshModel
        The shell mesh.
    ring_nodes : array-like of int
        Global node indices of one physical ring: positional indices into
        ``mesh.nodes`` / ``mesh.coords_array`` (the same indices
        ``MeshModel.node_id_to_index`` yields), not node ids.
    span_dir : array-like
        Unit span axis; the section plane is normal to it.
    element_properties : dict or None
        ``{element_set_name: property}`` as returned by
        ``Blade.get_element_properties()``.  The property is read through the
        ``mesh.element_sets[set_name].elements`` membership the assembler itself
        uses (see :func:`_element_property_by_id`).  When it is ``None``, or an
        element's set has no entry, every wall falls back to a **uniform
        ``S = 1.0`` N/m** and the returned split is **geometric only**; the caller
        can therefore tell the two cases apart by inspecting the stiffness values.

    The wall graph
    --------------
    For every mesh element the in-plane edges are the pairs of its nodes that
    **both** belong to ``ring_nodes``.  A shell element spans two rings, so exactly
    its chordwise edges lie in one ring (a quad ``[A_a, A_b, B_b, B_a]`` has the
    chordwise edges ``(A_a, A_b)`` and ``(B_b, B_a)``); its spanwise edges and
    diagonals do not.  Edges are enumerated along the element's node cycle, never
    across it, so a quad's ``(n0, n2)`` / ``(n1, n3)`` diagonals cannot arise, and
    only the corner cycle of the triangle / quad families is used so a quadratic
    element's mid-side nodes cannot manufacture a spurious pair.

    Per-edge shear stiffness
    ------------------------
    ``S = G * t`` (N/m) is read from the owning element's property via index
    ``(2, 2)`` of the merged membrane stiffness ``Cm = C * h`` (``a_matrix()[2, 2]``
    for the Rust laminate, the ``"cm"`` entry or ``G * t`` for the assembler dicts).
    A chordwise edge is shared by the element above and the element below the ring;
    adjacent elements across a ring normally share the same laminate, and the first
    element encountered in ``mesh.elements`` order is taken as the owner.

    Returns
    -------
    (cells, adjacency, edge_length, edge_shear_stiffness)
        ``cells`` and ``adjacency`` come from
        :func:`section_contour.section_cells` / :func:`section_contour.cell_adjacency`.
        ``edge_length`` (m) and ``edge_shear_stiffness`` (N/m) are keyed by the
        **undirected local ring edge** ``(min(a, b), max(a, b))`` with ``a, b`` local
        indices into ``ring_nodes`` -- the same indexing the cells use.

    Raises
    ------
    ValueError
        When the ring has fewer than three nodes, when no element joins two ring
        nodes, or when the wall graph bounds no cell.

    Pure: the mesh and ``element_properties`` are read, never written.
    """
    ring = np.asarray(ring_nodes, dtype=np.intp).ravel()
    if ring.size < 3:
        raise ValueError(f"ring_section needs at least 3 ring nodes, got {ring.size}")

    coords = mesh.coords_array
    ring_points = coords[ring]
    local_of_global = {int(global_index): local for local, global_index in enumerate(ring)}
    element_property = _element_property_by_id(mesh, element_properties)

    edge_owner: dict[tuple[int, int], int] = {}
    for element in mesh.elements:
        if element.element_type in (ElementType.triangle, ElementType.triangle6):
            corner_count = 3
        elif element.element_type in (
            ElementType.quad,
            ElementType.quad8,
            ElementType.quad9,
        ):
            corner_count = 4
        else:
            continue
        element_indices = [mesh.node_id_to_index[node_id] for node_id in element.node_ids]
        corners = element_indices[:corner_count]
        for a, b in zip(corners, corners[1:] + corners[:1], strict=True):
            if a not in local_of_global or b not in local_of_global or a == b:
                continue
            la, lb = local_of_global[a], local_of_global[b]
            key = (la, lb) if la < lb else (lb, la)
            edge_owner.setdefault(key, int(element.id))

    if not edge_owner:
        raise ValueError(
            "ring_section found no in-plane edge: no element chordwise edge joins two ring nodes"
        )

    cells = section_cells(ring_points, list(edge_owner), span_dir)
    if not cells:
        raise ValueError(
            "ring_section found no closed cell: the ring's in-plane edges bound no face"
        )

    adjacency = {key: tuple(owners) for key, owners in cell_adjacency(cells).items()}
    edge_length: dict[tuple[int, int], float] = {}
    edge_shear_stiffness: dict[tuple[int, int], float] = {}
    for key in edge_owner:
        a, b = key
        edge_length[key] = float(np.linalg.norm(ring_points[a] - ring_points[b]))
        stiffness = _membrane_shear_stiffness(element_property.get(edge_owner[key]))
        edge_shear_stiffness[key] = 1.0 if stiffness is None else float(stiffness)

    return cells, adjacency, edge_length, edge_shear_stiffness


def multi_cell_shear_flow(
    cells: list[SectionCell],
    adjacency: dict[tuple[int, int], tuple[int, ...]],
    edge_length: dict[tuple[int, int], float],
    edge_shear_stiffness: dict[tuple[int, int], float],
    torsion: float,
) -> tuple[list[float], float]:
    """Bredt-Batho shear flows of a multi-cell thin-walled section under torque ``T``.

    A singly closed section carries the one wall flow ``q = T / (2 A)``.  A
    section split by internal webs into several cells carries one flow per cell,
    coupled through the shared webs.  For every cell the solver enforces

        sum over the walls of cell i of  q_wall * ell_wall / S_wall = 2 * A_i * theta'

    where the wall flow is ``q_i`` on a wall that bounds only cell i and
    ``q_i - q_j`` on a wall shared by cells i and j, and ``S = G * t`` is the
    wall's **membrane shear stiffness** in N/m.  Together with the torque
    equilibrium ``sum_i 2 * A_i * q_i = T`` this is the ``(n + 1)`` linear system

        [ D       -2 a ] [ q   ]   [ 0 ]
        [ 2 a^T    0   ] [ theta'] = [ T ]

    with ``D[i, i] = sum ell / S`` over cell i's walls, ``D[i, j] = -ell_ij / S_ij``
    on a wall shared by cells i and j (subtracted once per owner), ``a_i = A_i``
    and ``theta'`` the last unknown.  It is solved in one call to
    :func:`numpy.linalg.solve` -- no iteration, no tolerance, no regularisation.

    Units.  ``q * ell / S`` is dimensionless (it is ``gamma`` integrated along the
    wall), so the compatibility equation puts ``q`` and ``2 A theta'`` on the same
    footing: ``q`` comes out in N/m and ``theta'`` in rad/m.  Using one ``S`` per
    wall -- rather than one ``G`` and one ``t`` -- is what lets the spar caps, the
    skin and the webs carry their own laminates and split the flow accordingly; the
    torsion constant ``J`` is then the caller's business, recovered as
    ``J = T / theta'`` for the section's effective ``G``.

    Parameters
    ----------
    cells : list[SectionCell]
        Bounded faces of the section's planar wall graph, as returned by
        :func:`section_contour.section_cells`.  ``cell.edges`` are the walls and
        ``cell.area`` the enclosed area in m^2.
    adjacency : dict[tuple[int, int], tuple[int, ...]]
        Wall-to-owner map from :func:`section_contour.cell_adjacency`: each
        undirected wall key maps to the one cell it bounds or to the two cells
        that share it.
    edge_length, edge_shear_stiffness : dict[tuple[int, int], float]
        Per-wall length (m) and membrane shear stiffness ``S = G * t`` (N/m),
        keyed by the undirected edge ``(min(a, b), max(a, b))``.
    torsion : float
        Torque about the section's span axis in N.m.

    Returns
    -------
    (flows, theta_rate)
        ``flows[i]`` is the shear flow of cell ``i`` in N/m, in ``cells`` order;
        ``theta_rate`` is the twist rate in rad/m (``= T / (G_eff * J)``).

    Raises
    ------
    ValueError
        For fewer than one cell, a non-positive cell area, a non-positive wall
        stiffness, a wall missing from ``edge_length``/``edge_shear_stiffness``, a
        wall owned by more than two cells, or a singular system.

    Pure: ``cells``, ``adjacency`` and both dicts are read, never written.
    """
    n_cells = len(cells)
    if n_cells < 1:
        raise ValueError(f"multi_cell_shear_flow needs at least one cell, got {n_cells}")

    matrix = np.zeros((n_cells + 1, n_cells + 1), dtype=float)
    rhs = np.zeros(n_cells + 1, dtype=float)

    for i, cell in enumerate(cells):
        area = float(cell.area)
        if not area > 0.0:
            raise ValueError(f"cell {i} has non-positive area {area!r}")
        for a, b in cell.edges:
            key = (a, b) if a < b else (b, a)
            if key not in edge_length:
                raise ValueError(f"cell {i} wall {key} is missing from edge_length")
            if key not in edge_shear_stiffness:
                raise ValueError(f"cell {i} wall {key} is missing from edge_shear_stiffness")
            stiffness = float(edge_shear_stiffness[key])
            if not stiffness > 0.0:
                raise ValueError(
                    f"cell {i} wall {key} has non-positive shear stiffness {stiffness!r}"
                )
            ratio = float(edge_length[key]) / stiffness
            matrix[i, i] += ratio
            owners = adjacency.get(key, (i,))
            if len(owners) > 2:
                raise ValueError(
                    f"cell {i} wall {key} is owned by {len(owners)} cells; a thin-walled "
                    "section wall is shared by at most two"
                )
            others = [int(owner) for owner in owners if int(owner) != i]
            if len(others) == 1:
                matrix[i, others[0]] -= ratio
        # Compatibility row i: sum_j D[i, j] q_j - 2 A_i theta' = 0, so the
        # `-2 A_i` coefficient lives in theta''s (last) column.
        matrix[i, n_cells] = -2.0 * area
        matrix[n_cells, i] = 2.0 * area

    # Equilibrium row: sum_i 2 A_i q_i = T.
    rhs[n_cells] = float(torsion)

    try:
        solution = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "multi-cell shear-flow system is singular: the cell/wall graph is degenerate"
        ) from exc

    flows = [float(solution[i]) for i in range(n_cells)]
    return flows, float(solution[n_cells])
