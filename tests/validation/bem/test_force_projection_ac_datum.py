"""Production-datum guard for the BEM -> shell force projection (issue #9).

The aerodynamic centre (AC) of every BEM strip must sit at the station's
``airfoil.aerodynamic_center`` fraction of the chord, measured from the
**true** leading edge.  ``ForceProjector`` used to infer the leading edge as
the node with the minimum chordwise projection, an orientation set by an
unrelated reference axis (``normal_direction``).  On the repo's own
IEA-15MW blade mesh the ``+x`` direction points toward the **leading** edge,
so ``min`` selected the trailing edge and the AC was applied at
``1 - ac_frac`` of the chord instead of ``ac_frac``.

The three assertions below use the real blade mesh and the real AeroDyn
``BladeAero`` (not a synthetic stub) and are deliberately independent of the
projector's own AC arithmetic:

1. every mesh node belongs to a strip (the strip grid must not be offset
   from the mesh datum by the hub radius);
2. the AC lands at the station's chord fraction from the geometrically
   identified leading edge (the **blunt** end of the section);
3. the projected nodal forces reproduce the analytic moment about the
   origin, ``sum_k [ r_ac_k x F_k + M_AC_k ]`` with a **3-D** geometric AC -
   the check ``ForceProjector.verify`` cannot make, because it balances
   force only and never inspects ``M_strip``.

Assertion 2 identifies the leading edge convention-free: within the outer
quarter of the chord, the blunt end has the larger in-plane thickness (the
spread of the nodes perpendicular to ``chord_dir``, inside the plane normal
to the span).  When the two thicknesses are within 10 % the blunt/sharp split
is ill-conditioned (a circular root section) and the ``+chord_dir`` end is
taken, which is the mesh's own leading-edge convention; a symmetric section is
where the two candidate AC placements are least distinguishable.

The 2 %-of-chord datum bound and the 1 % relative moment bound are the
task's bounds; they are not fitted to the measurement.
"""

from __future__ import annotations


import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

from tests.support.paths import DATA_DIR  # noqa: E402
YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AD_PRIMARY = (
    DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
)

#: Production defaults of the standalone / FSI projectors.
SPAN_DIR = np.array([0.0, 0.0, 1.0])
NORMAL_DIR = np.array([1.0, 0.0, 0.0])
TANGENTIAL_DIR = np.array([0.0, 1.0, 0.0])

DATUM_TOL = 0.02  # AC placement, as a fraction of the local chord
MOMENT_TOL = 0.01  # relative error on |M_expected|
END_SLAB = 0.25  # outer quarter of the chord, per end
THICKNESS_TIE_REL = 0.10  # <=10 % thickness difference is a tie


def _blunt_end(pts: np.ndarray, chord_dir: np.ndarray, span_dir: np.ndarray):
    """Identify the leading (blunt) and trailing (sharp) chordwise ends.

    Returns ``(le_index, te_index, thickness_lo, thickness_hi)`` where the
    thicknesses are the in-plane node spreads in the outer quarter of the
    chord nearest the ``min``/``max`` projection.  The larger spread is the
    blunt leading edge; a tie falls back to the ``+chord_dir`` end.
    """
    p = pts @ chord_dir
    lo, hi = float(p.min()), float(p.max())
    c = hi - lo
    w = np.cross(span_dir, chord_dir)
    w_norm = float(np.linalg.norm(w))
    if c <= 1e-12 or w_norm <= 1e-12:
        return int(np.argmax(p)), int(np.argmin(p)), 0.0, 0.0
    w = w / w_norm
    slab = END_SLAB * c
    q_lo = pts[p <= lo + slab] @ w
    q_hi = pts[p >= hi - slab] @ w
    t_lo = float(q_lo.max() - q_lo.min()) if q_lo.size else 0.0
    t_hi = float(q_hi.max() - q_hi.min()) if q_hi.size else 0.0
    if t_lo > t_hi and (t_lo - t_hi) > THICKNESS_TIE_REL * max(t_lo, t_hi):
        return int(np.argmin(p)), int(np.argmax(p)), t_lo, t_hi
    # Blunt at +chord_dir, or an (effectively) symmetric section: the
    # +chord_dir end is the mesh's leading-edge convention.
    return int(np.argmax(p)), int(np.argmin(p)), t_lo, t_hi


def _section_chord_axis(pts: np.ndarray, span_dir: np.ndarray) -> np.ndarray:
    """Principal in-plane axis of a strip outline, from the ring geometry alone.

    The section chord *axis* comes from an SVD of the strip's in-plane node
    offsets - never from ``ForceProjector`` - so a wrong production frame
    cannot make the expectation below agree with it by construction.  The
    sense is resolved once for the whole blade in
    :func:`_geometry_load_frames`.
    """
    off = pts - pts.mean(axis=0)
    off_plane = off - np.outer(off @ span_dir, span_dir)
    _, _, vt = np.linalg.svd(off_plane, full_matrices=False)
    axis = vt[0]
    return axis / np.linalg.norm(axis)


def _geometry_load_frames(strips, coords, span_dir, normal_dir, tangential_dir):
    """Per-strip ``(chord_hat, normal_hat)`` from the ring outlines alone.

    Each axis is the principal in-plane axis of that strip's own nodes
    (SVD), made continuous along the span, and the two senses are resolved
    **once for the whole blade** from the configured reference directions
    with the same rule the production code documents::

        chord_sign  = sign( sum_k dr_k * raw_chord_k  . tangential_dir )
        normal_sign = sign( sum_k dr_k * (raw_chord_k x span_dir) . normal_dir )

    ``projector._strip_chord_dirs`` / ``_strip_normal_dirs`` are deliberately
    **not** read: an expectation built from the implementation's own frame
    validates algebra, not physics (suite audit section 22.1).
    """
    raw = []
    prev = None
    for strip in strips:
        axis = _section_chord_axis(coords[strip.node_indices], span_dir)
        if prev is not None and float(axis @ prev) < 0.0:
            axis = -axis
        prev = axis
        raw.append(axis)

    weight = np.array([float(strip.dr) for strip in strips])
    chord_sum = sum(weight[k] * raw[k] for k in range(len(raw)))
    normal_sum = sum(
        weight[k] * np.cross(raw[k], span_dir) for k in range(len(raw))
    )
    chord_sign = 1.0 if float(chord_sum @ tangential_dir) >= 0.0 else -1.0
    normal_sign = 1.0 if float(normal_sum @ normal_dir) >= 0.0 else -1.0

    frames = []
    for axis in raw:
        chord_hat = chord_sign * axis
        chord_hat = chord_hat / np.linalg.norm(chord_hat)
        normal_hat = normal_sign * np.cross(chord_hat, span_dir)
        normal_hat = normal_hat / np.linalg.norm(normal_hat)
        frames.append((chord_hat, normal_hat))
    return frames


@pytest.fixture(scope="module")
def blade_case():
    """The repo's real IEA-15MW blade mesh, AeroDyn aero and projector."""
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None, "Blade.generate_mesh() produced no mesh"
    coords = mesh.coords_array
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    projector = ForceProjector(
        mesh,
        blade_aero,
        span_direction=SPAN_DIR,
        normal_direction=NORMAL_DIR,
        tangential_direction=TANGENTIAL_DIR,
    )
    return mesh, coords, blade_aero, projector


@pytest.fixture(scope="module")
def rated_bem(blade_case):
    """Bottom-up rated BEM result: Np, Tp and a physical (nose-down) Mp."""
    _, _, blade_aero, _ = blade_case
    solver = BEMSolver(
        blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0
    )
    return solver.compute(10.59, 7.56, 0.0)


def test_every_mesh_node_is_assigned(blade_case):
    """The strip grid must be aligned with the mesh span datum.

    The BEM stations are hub-referenced while the single-blade mesh is
    blade-root-referenced (span 0..117 m).  If the strip boundaries are not
    shifted onto the mesh datum, the root nodes below the first boundary are
    silently dropped (and the tip stations get no nodes at all).
    """
    _, coords, _, projector = blade_case
    assigned = np.zeros(coords.shape[0], dtype=bool)
    empty_strips = []
    for k, strip in enumerate(projector._strips):
        if len(strip.node_indices) == 0:
            empty_strips.append(k)
        assigned[strip.node_indices] = True
    unassigned = np.where(~assigned)[0]
    if len(unassigned):
        z = coords[unassigned, 2]
        print(
            f"\nunassigned nodes: {len(unassigned)}/{coords.shape[0]}, "
            f"z in [{z.min():.3f}, {z.max():.3f}] m"
        )
    print(f"empty strips: {empty_strips}")
    assert len(unassigned) == 0, (
        f"{len(unassigned)} mesh nodes were not assigned to any strip "
        f"(z in [{coords[unassigned, 2].min():.3f}, "
        f"{coords[unassigned, 2].max():.3f}] m)"
    )
    assert not empty_strips, f"strips without nodes: {empty_strips}"


def test_aerodynamic_centre_datum(blade_case):
    """The applied AC must sit at ``ac_frac`` of the chord from the true LE."""
    _, coords, blade_aero, projector = blade_case
    worst = 0.0
    print(f"\n{'k':>3} {'chord':>7} {'thk_lo':>8} {'thk_hi':>8} "
          f"{'frac_from_LE':>12} {'expected':>8}")
    for k, strip in enumerate(projector._strips):
        if len(strip.node_indices) < 2:
            continue
        idx = strip.node_indices
        chord_dir = projector._strip_chord_dirs[k]
        pts = coords[idx]
        le_i, te_i, t_lo, t_hi = _blunt_end(pts, chord_dir, SPAN_DIR)
        p = pts @ chord_dir
        le_proj = float(p[le_i])
        te_proj = float(p[te_i])
        chord = float(p.max() - p.min())
        ac_frac = blade_aero.stations[k].airfoil.aerodynamic_center
        ac_proj_expected = le_proj + ac_frac * (te_proj - le_proj)

        # Read the datum the projection actually uses from its own arm.
        centroid_proj = float(strip.centroid @ chord_dir)
        ac_proj_actual = centroid_proj - float(
            projector._strip_ac_offsets[k] @ chord_dir
        )
        # Fraction of the chord from the true leading edge toward the true
        # trailing edge: 0.25 is the physical value, ~0.75 the defect.
        frac_from_le = (le_proj - ac_proj_actual) / (le_proj - te_proj)
        error = abs(ac_proj_actual - ac_proj_expected) / chord
        worst = max(worst, error)
        if k % 5 == 0 or error > DATUM_TOL:
            print(
                f"{k:>3} {chord:>7.3f} {t_lo:>8.4f} {t_hi:>8.4f} "
                f"{frac_from_le:>12.3f} {ac_frac:>8.2f}"
            )
        assert error <= DATUM_TOL, (
            f"strip {k}: AC at {frac_from_le:.3f} of the chord from the true "
            f"leading edge, expected {ac_frac:.3f} (error {error:.3f} c)"
        )
    print(f"worst datum error: {worst:.4f} c")


def test_moment_conservation(blade_case, rated_bem):
    """The projected forces must reproduce the geometry-derived applied moment.

    ``ForceProjector.verify`` balances total force only, so a wrong section
    frame, a wrong AC or a wrong AC->centroid arm is invisible to it.  Here
    the applied moment ``sum_j r_j x f_j`` of the returned nodal forces is
    compared with

        sum_k [ r_ac_k x F_k + Mp_k dr_k span_dir ],
        F_k = Np_k dr_k n_hat_k + Tp_k dr_k c_hat_k,

    where **every frame** ``(c_hat_k, n_hat_k)``, the AC point ``r_ac_k`` and
    the strip ends are derived here from the ring outlines - the principal
    in-plane axis of each strip's nodes (:func:`_geometry_load_frames`) and
    the blunt-end LE rule (:func:`_blunt_end`) - with only the blade-wide
    *sign* convention taken from the configured directions.  Nothing is read
    from ``projector._strip_chord_dirs`` / ``_strip_normal_dirs``: an expected
    value built from the implementation's own frame validates algebra, not
    physics (suite audit section 22.1).  This geometry-derived expectation is
    what makes the 1 % bound meaningful - it can only hold if production
    follows each section, not a fixed global axis.
    """
    _, coords, blade_aero, projector = blade_case
    forces = projector.project(rated_bem)
    moment_applied = np.cross(coords, forces).sum(axis=0)

    frames = _geometry_load_frames(
        projector._strips, coords, SPAN_DIR, NORMAL_DIR, TANGENTIAL_DIR
    )
    moment_expected = np.zeros(3)
    mp_total = 0.0
    for k, strip in enumerate(projector._strips):
        idx = strip.node_indices
        assert len(idx) > 0, f"strip {k} has no nodes; cannot be loaded"
        pts = coords[idx]
        chord_hat, normal_hat = frames[k]
        le_i, te_i, _, _ = _blunt_end(pts, chord_hat, SPAN_DIR)
        ac_frac = blade_aero.stations[k].airfoil.aerodynamic_center
        r_ac = pts[le_i] + ac_frac * (pts[te_i] - pts[le_i])
        F = (
            float(rated_bem.Np[k]) * strip.dr * normal_hat
            + float(rated_bem.Tp[k]) * strip.dr * chord_hat
        )
        M_ac = float(rated_bem.Mp[k]) * strip.dr * SPAN_DIR
        moment_expected += np.cross(r_ac, F) + M_ac
        mp_total += abs(float(rated_bem.Mp[k]) * strip.dr)

    err = float(np.linalg.norm(moment_applied - moment_expected))
    rel = err / float(np.linalg.norm(moment_expected))
    print(
        f"\n|M_applied| = {np.linalg.norm(moment_applied):.6e} N*m\n"
        f"|M_expected| = {np.linalg.norm(moment_expected):.6e} N*m\n"
        f"|M_applied - M_expected| = {err:.6e} N*m ({rel:.4e} relative)\n"
        f"sum_k |Mp*dr| = {mp_total:.6e} N*m; error = {err / mp_total:.2f} x Mp"
    )
    assert rel <= MOMENT_TOL, (
        f"moment error {err:.6e} N*m is {rel:.4f} of |M_expected| "
        f"(bound {MOMENT_TOL})"
    )
