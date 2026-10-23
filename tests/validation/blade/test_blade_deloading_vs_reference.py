"""One de-loading table, one sign convention, six BEM feedbacks (issue #9, WU-D).

The production one-way path feeds the BEM two entangled feedbacks - the deformed strip radii
and the elastic section twist.  This module drives ``BEMFSIParticipant`` on the real IEA-15MW
blade mesh, the real AeroDyn ``BladeAero`` and the rated point (V = 10.59 m/s, 7.56 rpm, pitch
0), solves the shell **once** under the participant's own loads, and re-runs only the (cheap)
BEM on the feedback combinations - one table in one sign convention:

    thrust and power as (flexible - rigid) / rigid, negative = the rotor unloads; twist in
    degrees about +span under the declared convention fluid +Y / rotor clockwise from behind.

Rows: twist only (pure bend-twist unloading), the radii family as four rows - the radius
definition bias alone (no deformation), the deformation increment with the projection metric,
the same increment with a path metric of the deformed centroid line, and the two combined as
production feeds them - plus twist + radii (production).  The path metric is measured beside
the feedback and never fed to the BEM.

The finding: the participant's radius is the span projection of the deformed strip line - the
strip's own spanwise displacement plus the second-order shortening the linear solve cannot
produce (``radii_linear_kinematics_axial_term``) - while a path measure of the same
configuration lengthens with the transverse deflection instead.  The
table is **reported**; the only asserts are structural (rigid path invariance, twist sign,
``sum F`` conservation).  No magnitude is asserted against Zhou et al. 2025 Table 6
(-13.04% thrust / -8.38% power, flexible vs rigid).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_validation  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

from tests.support.paths import DATA_DIR  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"

V_RATED, RPM_RATED, PITCH_RATED = 10.59, 7.56, 0.0
BEM_CONFIG = {
    "wind_speed": V_RATED,
    "omega": RPM_RATED * 2.0 * np.pi / 60.0,
    "pitch": PITCH_RATED,
    "azimuth": 0.0,
    "air_density": 1.225,
    "dynamic_viscosity": 1.81206e-5,
    "hub_height": 150.0,
    "shear_exp": 0.0,
}
#: Declared sign convention, stated once and used by the header and the sign assertion.
#: ``tools/diagnose_sign_chain.py`` settles the sense: the deck's leading edge sits at +x,
#: the load-frame downwind (thrust) direction is +y, and a rigid +z rotation of the tip
#: ring moves the leading edge downwind, so a NOSE-DOWN section rotation is ``omega > 0``.
CONVENTION = (
    "fluid +Y / rotor clockwise viewed from behind; twist is the section rotation "
    "in degrees about +span, positive = nose-down toward feather (reduces the AoA)"
)
SENSE_THRUST_TOL = 0.02  # the applied-load invariant of verdict section 22.4 / the P5 guard
ZHOU_FLEX_THRUST_DELTA, ZHOU_FLEX_POWER_DELTA = -0.1304, -0.0838  # Zhou 2025 Table 6
#: Float64 accumulation bound for the span-projection identity (issue #13).  The test recomputes
#: the participant's increment from the same displacements with independent numpy code, so the
#: two must agree to rounding only; this is not a physical tolerance.
SPAN_PROJECTION_TOL = 1e-9
#: BeamDyn (OpenFAST v5.0.0, geometrically exact, static, dead loads) under this fixture's own
#: binned aero load, official IEA-15 BeamDyn deck: axial displacement at the tip-strip centroid
#: span (z0 = 116.403 m).  ``tools/diagnose_radii_beamdyn.py`` regenerates it.  The shell flaps
#: 9.15 % less than BeamDyn here (17.007 vs 18.720 m), so the second-order term, which scales with
#: the flap squared, is not expected to match it closely; only its sign is asserted.
BEAMDYN_TIP_AXIAL_M = -1.3146


@pytest.fixture(scope="module")
def deloading():
    """One shell solve under the participant's loads, then the baseline + three BEM feedbacks.

    Built as ``tests/test_bem_fsi_deformed_geometry.py`` builds it.  Every row shares one
    reference baseline (``_rebuild_bem_solver(_ref_r, _ref_twist)``), so the table's
    denominators are one construction.
    """
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=0.25)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None, "the blade mesh failed to generate"
    props = model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props),
        props,
        list(blade_validation.SPAN_DIRECTION),
        None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    free = np.array(
        [i for i in range(n) if i not in {6 * r + d for r in root for d in range(6)}],
        dtype=np.int64,
    )
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes], dtype=float)

    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    participant = BEMFSIParticipant(mesh, blade_aero, BEM_CONFIG)
    forces_rigid, bem_rigid = participant._compute_forces(np.zeros((len(coords), 3)))
    load = np.zeros(n)
    for dof in range(3):
        load[dof::6] = forces_rigid[:, dof]
    u = np.zeros(n)
    u[free] = np.asarray(spsolve(K[np.ix_(free, free)], load[free]))
    displacements = np.column_stack([u[0::6], u[1::6], u[2::6]])
    r_def, twist_def = participant._compute_deformed_geometry(displacements)
    # The same call at zero displacement: the participant's own radius definition on the
    # reference configuration.  Its offset from ``_ref_r`` was the mesh-vs-deck definition
    # bias (issue #13; removed for #33 L5, now pinned to zero by
    # ``test_zero_displacement_returns_the_deck_geometry``); ``dr_proj`` is the deformation
    # increment production feeds.
    r_def_rigid, _ = participant._compute_deformed_geometry(np.zeros_like(displacements))
    dr_proj = r_def - r_def_rigid
    # Strip centroids from the participant's own strip node sets, for the path metric.
    c_ref = np.stack(
        [participant._ref_coords[idx].mean(axis=0) for idx in participant._strip_node_indices]
    )
    c_def = np.stack(
        [
            (participant._ref_coords[idx] + displacements[idx]).mean(axis=0)
            for idx in participant._strip_node_indices
        ]
    )
    r_path_ref = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(c_ref, axis=0), axis=1))])
    r_path_def = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(c_def, axis=0), axis=1))])
    dr_path = r_path_def - r_path_ref
    L_ref, L_def = float(r_path_ref[-1]), float(r_path_def[-1])

    def bem(r: np.ndarray, twist: np.ndarray):
        solver, _ = participant._rebuild_bem_solver(r, twist)
        return solver.compute(V_RATED, RPM_RATED, PITCH_RATED)

    ref_r, ref_twist = participant._ref_r, participant._ref_twist
    # Back into the mesh's rotation sense about +span (fsi_participant._twist_mesh_to_bem).
    omega = participant._twist_mesh_to_bem * (twist_def - ref_twist)
    span = np.asarray(blade_validation.SPAN_DIRECTION, dtype=float)
    return {
        "participant": participant,
        "blade_aero": blade_aero,
        "coords": coords,
        "mesh": mesh,
        "props": props,
        "free": free,
        "Kff": K[np.ix_(free, free)],
        "n": n,
        "displacements": displacements,
        "r_def": r_def,
        "twist_def": twist_def,
        "ref_r": ref_r,
        "ref_twist": ref_twist,
        "omega": omega,
        "tip_k": int(np.argmax(ref_r)),
        "span": span,
        "tip_disp": displacements[int(np.argmax(coords @ span))],
        "forces_rigid": forces_rigid,
        "bem_rigid": bem_rigid,
        "base": bem(ref_r, ref_twist),
        "twist": bem(ref_r, twist_def),
        "radii": bem(r_def, ref_twist),
        "both": bem(r_def, twist_def),
        "r_def_rigid": r_def_rigid,
        "dr_proj": dr_proj,
        "dr_path": dr_path,
        "c_ref": c_ref,
        "c_def": c_def,
        "r_path_ref": r_path_ref,
        "r_path_def": r_path_def,
        "L_ref": L_ref,
        "L_def": L_def,
        "radii_dr_proj": bem(ref_r + dr_proj, ref_twist),
        "radii_dr_path": bem(ref_r + dr_path, ref_twist),
        "radii_bias": bem(r_def_rigid, ref_twist),
    }


def test_rigid_baseline_reproduces_the_standalone_bem(deloading):
    """Zero displacement must reproduce the standalone BEM exactly (path invariance).

    The 1e-9 relative bound is float64 accumulation across two CCBlade evaluations of the same
    rotor, not a physical tolerance.
    """
    _, bem = deloading["participant"]._compute_forces(np.zeros_like(deloading["displacements"]))
    reference = BEMSolver(
        deloading["blade_aero"], rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0
    ).compute(V_RATED, RPM_RATED, PITCH_RATED)
    thrust_rel = abs(bem.thrust - reference.thrust) / abs(reference.thrust)
    power_rel = abs(bem.power - reference.power) / abs(reference.power)
    print(
        f"\nrigid participant vs standalone BEM: thrust {bem.thrust / 1e6:.6f} MN "
        f"(rel {thrust_rel:.3e}), power {bem.power / 1e6:.6f} MW (rel {power_rel:.3e})"
    )
    assert thrust_rel < 1e-9
    assert power_rel < 1e-9


def test_twist_sign_matches_the_declared_convention(deloading):
    """The tip section rotation is nose-down about +span; the BEM twist increment must rise.

    A sign is an invariance, so nothing is bound.  ``tools/diagnose_sign_chain.py``
    settles the sense: the deck's leading edge sits at +x, the load-frame downwind
    (thrust) direction is +y, and a rigid +z rotation of the tip ring moves the leading
    edge downwind, so a NOSE-DOWN section rotation is ``omega > 0``.  With the pitching
    moment applied on the axis ``_section_ends`` selects
    (``ForceProjector._strip_moment_axis_sign``), the properties-less production path
    measures ``omega = +8.1048 deg``.  A nose-down twist reduces the angle of attack, and
    CCBlade's ``alpha = phi - theta`` needs ``theta`` to RISE for that, so the BEM twist
    increment must be positive - the SAME sense as the corrected section rotation.  It was
    pinned ``xfail(strict=True)`` while the applied pitching moment rode the raw span
    axis; that inversion is fixed and the marker is gone.
    """
    tip_k = deloading["tip_k"]
    omega_tip = float(deloading["omega"][tip_k])
    theta_bem_tip = float(deloading["twist_def"][tip_k] - deloading["ref_twist"][tip_k])
    print(f"\ndeclared convention: {CONVENTION}")
    print(
        f"tip strip r_ref {float(deloading['ref_r'][tip_k]):.4f} m: mesh section rotation "
        f"{np.rad2deg(omega_tip):+.4f} deg (nose-down), BEM twist increment "
        f"{np.rad2deg(theta_bem_tip):+.4f} deg (must be positive to reduce the AoA)"
    )
    assert omega_tip > 0.0, (
        f"the tip section rotation is {np.rad2deg(omega_tip):+.4f} deg, not nose-down "
        f"under the declared convention; on this frame nose-down is omega > 0 (deck "
        f"leading edge at +x, load-frame downwind at +y, a +z rotation moves the leading "
        f"edge downwind - tools/diagnose_sign_chain.py).  The corrected pitching-moment "
        f"axis gives +8.1048 deg"
    )
    # The physical claim: the BEM twist rises (``alpha = phi - theta``) so a nose-down
    # section rotation unloads the rotor.  The participant shares that sense
    # (``_twist_mesh_to_bem = +1``, fsi_participant.py); when it was -1 it compensated the
    # inverted pitching-moment axis fixed in 1146265 and turned the de-loading table's
    # twist-only row into a +21.88% RE-LOADING, which is what exposed the chain.
    assert theta_bem_tip > 0.0, (
        f"the BEM twist increment is {np.rad2deg(theta_bem_tip):+.4f} deg; a nose-down "
        f"section rotation (+{np.rad2deg(omega_tip):.4f} deg) must RISE the BEM twist "
        f"(alpha = phi - theta) so it unloads, so the increment must be positive"
    )


def test_multi_cell_twist_sign_matches_the_declared_convention(deloading):
    """The with-properties multi-cell path's measured sense: nose-down, like the fallback.

    The realisation itself is verified in ``tests/test_multicell_shear_flow.py``; this
    rebuilds the strip loads with the laminate map the fixture derives from the ``Blade``
    model (``get_element_properties``) and asserts the declared sense (nose-down is
    ``omega > 0``; tools/diagnose_sign_chain.py).  It was pinned ``xfail(strict=True)``
    while the applied pitching moment rode the raw span axis; once
    ``_strip_moment_axis_sign`` carries the polars' sense onto this frame the multi-cell
    path measures ``omega = +9.6669 deg`` - nose-down - so the marker is removed and the
    assertion kept.
    """
    participant = deloading["participant"]
    projector = ForceProjector(
        deloading["mesh"],
        deloading["blade_aero"],
        span_direction=participant._span_dir,
        element_properties=deloading["props"],
    )
    forces = projector.project(deloading["bem_rigid"])
    load = np.zeros(deloading["n"])
    for dof in range(3):
        load[dof::6] = forces[:, dof]
    u = np.zeros(deloading["n"])
    u[deloading["free"]] = np.asarray(spsolve(deloading["Kff"], load[deloading["free"]]))
    displacements = np.column_stack([u[0::6], u[1::6], u[2::6]])
    _, twist_def = participant._compute_deformed_geometry(displacements)
    omega = participant._twist_mesh_to_bem * (twist_def - deloading["ref_twist"])
    tip_k = deloading["tip_k"]
    omega_tip = float(omega[tip_k])
    theta_bem_tip = float(twist_def[tip_k] - deloading["ref_twist"][tip_k])
    print(
        f"\nmulti-cell tip strip r_ref {float(deloading['ref_r'][tip_k]):.4f} m: mesh "
        f"section rotation {np.rad2deg(omega_tip):+.4f} deg (nose-down, measured +9.6669 "
        f"deg), BEM twist increment {np.rad2deg(theta_bem_tip):+.4f} deg"
    )
    assert omega_tip > 0.0, (
        f"the multi-cell path gives a non-positive tip section rotation "
        f"{np.rad2deg(omega_tip):+.4f} deg; on this frame nose-down is omega > 0 "
        f"(tools/diagnose_sign_chain.py) and the corrected pitching-moment axis measures "
        f"+9.6669 deg, the declared sense"
    )


def test_de_loading_table_and_load_conservation(deloading):
    """Print the one de-loading table and assert only the applied-load ``sum F`` invariant.

    **Reported:** the three rows' thrust/power deltas (and Zhou's Table 6), the per-row tip
    axial-radius and radial-twist feedback, the tip nodal displacement, and the exact
    production-path cross-check.  **Asserted:** ``|sum(F)|`` against ``bem.thrust / n_blades``
    (the P5 guard's 2% bound).
    """
    base, tip_k = deloading["base"], deloading["tip_k"]
    dr_tip = float(deloading["r_def"][tip_k] - deloading["ref_r"][tip_k])
    dtheta_tip = float(np.rad2deg(deloading["omega"][tip_k]))
    tip_disp = deloading["tip_disp"]
    axial = float(tip_disp @ deloading["span"])
    radial = float(np.linalg.norm(tip_disp - axial * deloading["span"]))

    applied = deloading["forces_rigid"].sum(axis=0)
    n_blades = int(deloading["blade_aero"].n_blades)
    thrust_per_blade = float(deloading["bem_rigid"].thrust) / n_blades
    load_ratio = float(np.linalg.norm(applied)) / thrust_per_blade

    cases = (
        ("twist", "twist only, reference radii", "the pure bend-twist unloading", 0.0, dtheta_tip),
        (
            "radii",
            "deformed radii only, reference twist",
            "the radii feedback as production feeds it",
            dr_tip,
            0.0,
        ),
        (
            "both",
            "twist + radii (production path)",
            "what production does today",
            dr_tip,
            dtheta_tip,
        ),
        (
            "radii_dr_proj",
            "radii, deformation increment only (projection metric)",
            "the deformation increment production feeds, bias removed",
            float(deloading["dr_proj"][tip_k]),
            0.0,
        ),
        (
            "radii_dr_path",
            "radii, deformation increment only (path metric)",
            "what a path measure would add instead",
            float(deloading["dr_path"][tip_k]),
            0.0,
        ),
        (
            "radii_bias",
            "radii, definition bias only (zero displacement)",
            "the radius definition bias, no deformation (0 by construction)",
            float(deloading["r_def_rigid"][tip_k] - deloading["ref_r"][tip_k]),
            0.0,
        ),
    )
    print(f"\nDe-loading, one table, one sign convention: {CONVENTION}")
    print(
        "thrust/power = (flexible - rigid) / rigid, negative = the rotor unloads; twist [deg] "
        f"about +span; Zhou et al. 2025 Table 6: {ZHOU_FLEX_THRUST_DELTA:+.2%} thrust / "
        f"{ZHOU_FLEX_POWER_DELTA:+.2%} power"
    )
    print(
        f"rigid baseline (reference radii + reference twist): {base.thrust / 1e6:.6f} MN / "
        f"{base.power / 1e6:.6f} MW"
    )
    print(
        f"  {'feedback fed to the BEM':62} {'isolates':58} {'d thrust':>9} {'d power':>9} "
        f"{'tip dr [m]':>10} {'tip dtwist [deg]':>16}"
    )
    deltas = {}
    for key, feedback, isolates, dr, dtheta in cases:
        bem = deloading[key]
        dt = float(bem.thrust) / float(base.thrust) - 1.0
        dp = float(bem.power) / float(base.power) - 1.0
        deltas[key] = (dt, dp)
        print(
            f"  {feedback:62} {isolates:58} {dt:>+9.2%} {dp:>+9.2%} {dr:>+10.4f} {dtheta:>+16.4f}"
        )

    r_def_rigid = deloading["r_def_rigid"]
    print(
        "\nmetric block, the deformed-radius definition and the two increment metrics:\n"
        f"  definition bias r_def_rigid - _ref_r: strip 0 "
        f"{float(r_def_rigid[0] - deloading['ref_r'][0]):+.6f} m, tip "
        f"{float(r_def_rigid[tip_k] - deloading['ref_r'][tip_k]):+.6f} m\n"
        f"  tip increment: dr_proj {float(deloading['dr_proj'][tip_k]):+.6f} m, dr_path "
        f"{float(deloading['dr_path'][tip_k]):+.6f} m\n"
        f"  deformed centroid path length: L_ref {deloading['L_ref']:.6f} m, L_def "
        f"{deloading['L_def']:.6f} m, L_def/L_ref - 1 "
        f"{deloading['L_def'] / deloading['L_ref'] - 1.0:+.6f} "
        f"({deloading['L_def'] / deloading['L_ref'] - 1.0:+.4%})"
    )

    radii_power, both_power = deltas["radii"][1], deltas["both"][1]
    print(
        f"\npower sign, deformed radii only: {radii_power:+.2%} (positive = re-loads); "
        f"production path: {both_power:+.2%} (positive = re-loads); twist alone: "
        f"{deltas['twist'][1]:+.2%} -> "
        f"{'same sign' if np.sign(radii_power) == np.sign(both_power) else 'opposite signs'} "
        f"(radii only vs production)"
    )
    print(
        f"tip nodal displacement: axial (span) {axial:+.4f} m, radial (in-plane) {radial:+.4f} m "
        f"(u_tip = {np.array2string(tip_disp, precision=4)}); tip strip radius increment fed: "
        f"dr = {dr_tip:+.4f} m of {float(deloading['ref_r'][tip_k]):.4f} m"
    )
    _, bem_prod = deloading["participant"]._compute_forces(deloading["displacements"])
    print(
        f"cross-check, exact production path (deformed projector, own Rtip rule) vs the rigid "
        f"participant BEM: thrust {float(bem_prod.thrust) / float(deloading['bem_rigid'].thrust) - 1.0:+.2%}, "
        f"power {float(bem_prod.power) / float(deloading['bem_rigid'].power) - 1.0:+.2%}"
    )
    print(
        f"one-way rows production can feed (reported, not asserted): twist only "
        f"{deltas['twist'][0]:+.2%}, radii only {deltas['radii'][0]:+.2%}, twist + radii "
        f"{deltas['both'][0]:+.2%} thrust; Zhou {ZHOU_FLEX_THRUST_DELTA:+.2%}.  The path-metric "
        f"row printed above is a counterfactual that is never fed to the BEM; a converged fixed "
        f"point would deepen the twist."
    )
    print(
        f"\napplied-load invariant: |sum(F)| = {np.linalg.norm(applied):.4f} N vs bem.thrust/"
        f"{n_blades} = {thrust_per_blade:.4f} N -> ratio {load_ratio:.6f} "
        f"({load_ratio - 1.0:+.3%}, bound {SENSE_THRUST_TOL:.0%})"
    )
    assert abs(load_ratio - 1.0) < SENSE_THRUST_TOL, (
        f"the projected load |sum(F)| = {np.linalg.norm(applied):.4f} N is {load_ratio - 1.0:+.3%} "
        f"off bem.thrust/{n_blades} = {thrust_per_blade:.4f} N (the P5 guard's bound)"
    )


def test_zero_displacement_returns_the_deck_geometry(deloading):
    """At ``u = 0`` the deformed path is the rigid path: same annulus, same twist (#33 L5).

    The rigid short-circuit (``disp_max < 1e-12``) feeds the BEM the deck's ``_ref_r`` and
    ``_ref_twist``; the deformed path must reduce to exactly that, or its first iteration
    moves every station radius by a definition bias before any deformation (it was
    ``+0.397959 m`` at strip 0 and ``-0.596870 m`` at the tip strip while the radius was the
    node-cloud centroid; ``radii_datum_definition_bias``).  **Asserted:** exact equality, no
    tolerance: the radius is ``_ref_r[k] + mean(0)`` and the twist ``_ref_twist[k] + 0``.
    """
    participant = deloading["participant"]
    r0, twist0 = participant._compute_deformed_geometry(np.zeros_like(deloading["displacements"]))
    np.testing.assert_array_equal(r0, participant._ref_r)
    np.testing.assert_array_equal(twist0, participant._ref_twist)


def test_deformed_radii_are_the_second_order_span_projection(deloading):
    """The feedback's measuring stick is a second-order span projection, not a path length.

    The deformed radius (issue #13, #33 L5 ``radii_linear_kinematics_axial_term``) is

        r_def[k] = _ref_r[k] + mean(u[strip k] . s) - 1/2 int_root^{s_k} |u_perp'|^2 ds

    on the strip-centroid line, from the clamped root.  The first term is the linear solve's
    spanwise displacement (at the rated point almost entirely the prebent axis straightening);
    the second is the shortening of a near-inextensible axis under the same bending, which the
    linear solve cannot produce and which is about twice as large with the opposite sign.  The
    formula was checked against BeamDyn's own geometrically exact field under this load
    (``docs/validation_closures.md``): it reproduces BeamDyn's tip axial displacement to 0.3 %
    of the second-order term.  A path measure of the deformed centroid line lengthens instead;
    it is printed beside it and **never fed to the BEM**.

    **Reported:** the tip values of ``dr_proj`` (production's increment), ``dr_path``,
    ``L_def/L_ref - 1`` and BeamDyn's tip axial displacement.  **Asserted:** the identity above
    to float precision, recomputed independently, and the tip moving inward as BeamDyn's does.
    """
    participant = deloading["participant"]
    displacements = deloading["displacements"]
    span_dir = participant._span_dir
    idx_sets = participant._strip_node_indices
    u_mean = np.stack([displacements[idx].mean(axis=0) for idx in idx_sets])
    span = np.array([float(np.mean(participant._ref_coords[idx] @ span_dir)) for idx in idx_sets])
    order = np.argsort(span)
    u_perp = u_mean[order] - np.outer(u_mean[order] @ span_dir, span_dir)
    x = np.concatenate([[float(np.min(participant._ref_coords @ span_dir))], span[order]])
    slope2 = np.sum(np.gradient(np.vstack([np.zeros(3), u_perp]), x, axis=0) ** 2, axis=1)
    shortening = np.empty(len(idx_sets))
    shortening[order] = 0.5 * np.cumsum(0.5 * (slope2[1:] + slope2[:-1]) * np.diff(x))
    expected = u_mean @ span_dir - shortening

    dr_proj = deloading["dr_proj"]
    tip_k = deloading["tip_k"]
    ratio_minus_one = deloading["L_def"] / deloading["L_ref"] - 1.0
    print(
        f"\ntip increment: dr_proj {float(dr_proj[tip_k]):+.6f} m (fed to the BEM: first order "
        f"{float(u_mean[tip_k] @ span_dir):+.6f}, shortening {-float(shortening[tip_k]):+.6f}), "
        f"BeamDyn {BEAMDYN_TIP_AXIAL_M:+.4f} m; dr_path "
        f"{float(deloading['dr_path'][tip_k]):+.6f} m (path measure, not fed), "
        f"L_def/L_ref - 1 {ratio_minus_one:+.6f} ({ratio_minus_one:+.4%})"
    )
    deviation = float(np.max(np.abs(dr_proj - expected)))
    assert deviation < SPAN_PROJECTION_TOL, (
        f"the deformed-radius increment is {deviation:.3e} off mean(u_strip . s) minus the "
        f"second-order shortening (bound {SPAN_PROJECTION_TOL:.0e})"
    )
    assert np.sign(dr_proj[tip_k]) == np.sign(BEAMDYN_TIP_AXIAL_M), (
        f"tip radius increment {float(dr_proj[tip_k]):+.4f} m has the opposite sign of "
        f"BeamDyn's geometrically exact {BEAMDYN_TIP_AXIAL_M:+.4f} m"
    )
