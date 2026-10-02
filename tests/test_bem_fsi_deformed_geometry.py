"""Production guards for the BEM-FSI participant's deformed-geometry feedback.

Two confirmed defects in ``src/aeroelast/solvers/bem/fsi_participant.py``
(``BEMFSIParticipant``) live on the one-way de-loading path, which nothing in the suite
exercised through the production class:

**P6 - mixed datums.** ``_ref_r`` is hub-referenced (3.97 .. 120.97 m) while the deformed
radii were averaged from the mesh's span coordinate, blade-root-referenced (0 .. 117 m)
for a single-blade mesh; the comparison was minus the hub radius and the radii handed to
CCBlade started at ~0, so the deformed BEM ran on a corrupted annulus.

**P7 - the twist came from a node cloud that is not one section.** A strip spans
``dr ~ 2.39 m`` against ~0.63 m mesh stations, so a strip contains ~4 rings; on a blade
bent 16 m the PCA chord direction of the whole cloud follows the deformed arc.

The tests drive the production class on the real IEA-15MW blade mesh, the real AeroDyn
``BladeAero`` and the rated point (V = 10.59 m/s, 7.56 rpm, pitch 0), with the structural
response solved by the repo's own shell under the participant's own projected loads.
Reported, never asserted: the magnitude against Zhou et al. 2025 Table 6, and the shell's
nodal drilling rotation (measured, it is not the section rotation on this shell).
"""

from __future__ import annotations

from pathlib import Path

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
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402

import test_blade_iea15mw_validation as blade_validation  # noqa: E402
from _openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from test_blade_rated_twist import STATION_GAP_TOLERANCE, _physical_stations  # noqa: E402

YAML = Path(__file__).resolve().parent / "IEA-15-240-RWT.yaml"
AD_PRIMARY = (
    Path(__file__).resolve().parent
    / "reference"
    / "iea15mw_openfast"
    / "case"
    / "IEA-15-240-RWT_AeroDyn15.dat"
)

#: Rated point.  ``bem_config["omega"]`` is rad/s and ``_compute_forces`` converts back to
#: RPM, so the exact conversion is used: the rounded 0.79168 rad/s is 7.5599871 rpm
#: (-1.7e-6 relative), which is not invariance-scale.
V_RATED = 10.59
RPM_RATED = 7.56
OMEGA_RATED_RADS = RPM_RATED * 2.0 * np.pi / 60.0
PITCH_RATED = 0.0
BEM_CONFIG = {
    "wind_speed": V_RATED,
    "omega": OMEGA_RATED_RADS,
    "pitch": PITCH_RATED,
    "azimuth": 0.0,
    "air_density": 1.225,
    "dynamic_viscosity": 1.81206e-5,
    "hub_height": 150.0,
    "shear_exp": 0.0,
}

#: Zhou et al. 2025, *Energy* 336:138488, Table 6, flexible vs rigid.  Reported only.
ZHOU_FLEX_THRUST_DELTA = -0.1304
ZHOU_FLEX_POWER_DELTA = -0.0838

#: Physical order-of-magnitude bounds, justified in the test docstrings.
DATUM_SHIFT_MAX_M = 2.0
HUB_PROTECTION_FRACTION = 0.5
TWIST_AGREEMENT_DEG = 1e-6
OUTER_SPAN_FRACTION = 0.75


@pytest.fixture(scope="module")
def fsi():
    """Mesh + participant + one shell solve under the participant's projected loads.

    Built as the blade tests build it: ``Blade(..., element_size=1.0)``,
    ``PyMeshAssembler.from_model`` over ``_to_rust_mesh``, the root ring clamped in all
    6 DOF, and one ``spsolve``.
    """
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")

    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
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
    K = coo_matrix(
        (np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)
    ).tocsr()
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

    phys_stations = _physical_stations(coords)
    return {
        "participant": participant,
        "blade_aero": blade_aero,
        "coords": coords,
        "u": u,
        "displacements": np.column_stack([u[0::6], u[1::6], u[2::6]]),
        "bem_rigid": bem_rigid,
        "phys_stations": phys_stations,
        "phys_rings": {
            zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0]
            for zz in phys_stations
        },
    }


def _nearest_station(stations: list[float], target: float) -> float:
    """The merged physical station closest to ``target`` (both on the mesh datum)."""
    arr = np.asarray(stations, dtype=float)
    return float(arr[int(np.argmin(np.abs(arr - target)))])


def _ring_section_rotation(coords: np.ndarray, disp: np.ndarray, ring: np.ndarray) -> float:
    """Independent affine fit of one ring's in-plane field; returns its rotation [rad]."""
    pts = coords[ring]
    design = np.column_stack([pts[:, 0], pts[:, 1], np.ones(len(ring))])
    a12 = np.linalg.lstsq(design, disp[ring, 0], rcond=None)[0][1]
    a21 = np.linalg.lstsq(design, disp[ring, 1], rcond=None)[0][0]
    return 0.5 * (a21 - a12)


def test_rigid_path_matches_the_standalone_bem(fsi):
    """Zero displacement must reproduce the standalone BEM exactly.

    The short-circuit is an invariance of the code path, not a model, so the 1e-9
    relative bound is float64 accumulation across two CCBlade evaluations of the same
    rotor - not a fitted tolerance.
    """
    participant = fsi["participant"]
    _, bem = participant._compute_forces(np.zeros_like(fsi["displacements"]))
    reference = BEMSolver(
        fsi["blade_aero"], rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0
    ).compute(V_RATED, RPM_RATED, PITCH_RATED)

    thrust_rel = abs(bem.thrust - reference.thrust) / abs(reference.thrust)
    power_rel = abs(bem.power - reference.power) / abs(reference.power)
    print(f"\nrigid participant vs standalone BEM: thrust {bem.thrust / 1e6:.6f} MN "
          f"(rel {thrust_rel:.3e}), power {bem.power / 1e6:.6f} MW (rel {power_rel:.3e})")
    assert thrust_rel < 1e-9
    assert power_rel < 1e-9


def test_deformed_radii_stay_on_the_reference_datum(fsi):
    """The deformed radii must be an *elastic* change of the reference radii (P6).

    Both bounds are order-of-magnitude statements about a 117 m blade whose measured
    axial (span projection) change is +1.06 m at the tip and whose geometric shortening
    from a 16.4 m flapwise deflection is ``delta^2 / (2 r) ~ 1.1 m``; 2 m of elastic span
    change is already generous.  No hub-anchored station may end up inside
    ``0.5 * Rhub``.  A datum mismatch of minus the hub radius (-3.0 .. -4.2 m measured
    before the fix) is a rigid offset, not an elastic change, and fails both bounds.
    """
    participant = fsi["participant"]
    blade_aero = fsi["blade_aero"]
    r_def, _ = participant._compute_deformed_geometry(fsi["displacements"])
    shift = r_def - participant._ref_r

    print(f"\nP6 datum: r_ref {participant._ref_r.min():.4f} .. {participant._ref_r.max():.4f} m"
          f" | r_def {r_def.min():.4f} .. {r_def.max():.4f} m"
          f" | r_def - r_ref {shift.min():+.4f} .. {shift.max():+.4f} m"
          f" | max|shift| {np.max(np.abs(shift)):.4f} m (bound {DATUM_SHIFT_MAX_M})")
    assert float(np.max(np.abs(shift))) < DATUM_SHIFT_MAX_M, (
        f"the deformed radii are a rigid {shift.mean():+.4f} m off the reference datum "
        f"(max |r_def - r_ref| = {np.max(np.abs(shift)):.4f} m, bound {DATUM_SHIFT_MAX_M} m)"
    )
    assert float(r_def.min()) > HUB_PROTECTION_FRACTION * blade_aero.hub_radius, (
        f"the innermost deformed station is at {r_def.min():.4f} m, inside "
        f"{HUB_PROTECTION_FRACTION} * Rhub = "
        f"{HUB_PROTECTION_FRACTION * blade_aero.hub_radius:.4f} m"
    )


def test_elastic_twist_matches_the_nodal_section_rotation(fsi):
    """The elastic twist must be the section's rotation, not the strip cloud's (P7).

    The production twist is compared against an independent computation of the same
    physical quantity: the least-squares in-plane section rotation of the strip's own
    centre ring.  The 1e-6 deg bound is an exact invariance of the estimator, not a
    physical tolerance, so it fails the moment the twist is taken from another node cloud
    (measured before the fix: +2.3764 deg at the tip against this reference's -1.5112).

    Also asserted from physics: the outer-span section rotation is nose-down (negative
    about +span), the response to the nose-down pitching moment of every loaded station at
    rated (section 19).

    **Reported, not asserted:** the shell's mean nodal drilling rotation ``theta_z`` -
    measured, it is not the section rotation here (the tip ring reaches -3.8558 deg while
    its in-plane section rotation is -1.5112 deg, the same 2.34 deg gap sections 20 and
    22.4 record).  No translation-only estimator reaches -3.86 deg (rigid-rotation fit
    +0.325, chord-line rotation +0.381, Kabsch fit +0.451), so a bound between the two
    quantities would be fitted to that gap.
    """
    participant = fsi["participant"]
    blade_aero = fsi["blade_aero"]
    u = fsi["u"]
    _, twist_def = participant._compute_deformed_geometry(fsi["displacements"])
    # Back into the mesh's rotation sense about +span: the BEM twist angle and the mesh
    # rotation are opposite-sensed on this frame (participant._twist_mesh_to_bem).
    section = participant._twist_mesh_to_bem * (twist_def - participant._ref_twist)

    # The BEM stations are hub-referenced; the mesh's blade-root datum starts at 0.  This
    # only localises the ring - it is not the P6 guard.
    mesh_station = participant._ref_r - blade_aero.hub_radius
    outer = np.where(mesh_station >= OUTER_SPAN_FRACTION * blade_aero.blade_length)[0]

    rows = []
    for k in outer:
        ring = fsi["phys_rings"][_nearest_station(fsi["phys_stations"], mesh_station[k])]
        rows.append((
            k,
            participant._ref_r[k],
            len(ring),
            float(np.rad2deg(section[k])),
            float(np.rad2deg(
                _ring_section_rotation(fsi["coords"], fsi["displacements"], ring))),
            float(np.rad2deg(np.mean(u[6 * ring + 5]))),
        ))

    print("\nP7 twist, outer quarter of the span (degrees about +span):")
    print(f"  {'strip':>5} {'r_ref[m]':>9} {'nodes':>6} {'production':>11} "
          f"{'ring section':>13} {'nodal theta_z':>14}")
    for k, r_ref, count, prod, reference, nodal in rows:
        print(f"  {k:>5} {r_ref:>9.2f} {count:>6} {prod:>11.4f} {reference:>13.4f} "
              f"{nodal:>14.4f}")
    gap = np.array([abs(row[3] - row[4]) for row in rows])
    nodal_gap = np.array([abs(row[3] - row[5]) for row in rows])
    print(f"  production vs the ring's own section rotation: max {gap.max():.3e} deg; "
          f"vs the nodal drilling rotation: max {nodal_gap.max():.4f} deg (reported)")

    assert float(gap.max()) < TWIST_AGREEMENT_DEG, (
        f"the production elastic twist misses the centre ring's own section rotation by "
        f"{gap.max():.3e} deg: it is not the section's rotation"
    )
    assert all(row[3] < 0.0 for row in rows), (
        "an outer-span section rotates nose-up under the rated nose-down pitching moment"
    )


def test_one_way_de_loading_reduces_thrust_and_power(fsi):
    """The elastic twist must de-load the rotor (sign invariance).

    Asserted: with the elastic twist as the *only* feedback (reference station radii, same
    solver construction) both thrust and power fall - before the fix the inverted twist
    sense turned that into a re-loading.  On the full production path, thrust falls.

    **Reported, never asserted:** the magnitudes against Zhou et al. 2025 Table 6
    (-13.04% / -8.38%), and the production net **power**, which rises: the deformed radii
    grow by up to +0.93 m of span projection, adding +0.85% of power against the twist's
    -0.81%.  A one-way application of rigid-blade loads against a coupled solution has no
    defensible bound.
    """
    participant = fsi["participant"]
    _, bem_rigid = participant._compute_forces(np.zeros_like(fsi["displacements"]))
    _, bem_def = participant._compute_forces(fsi["displacements"])
    _, twist_def = participant._compute_deformed_geometry(fsi["displacements"])

    base_solver, _ = participant._rebuild_bem_solver(
        participant._ref_r, participant._ref_twist
    )
    twist_solver, _ = participant._rebuild_bem_solver(participant._ref_r, twist_def)
    bem_base = base_solver.compute(V_RATED, RPM_RATED, PITCH_RATED)
    bem_twist = twist_solver.compute(V_RATED, RPM_RATED, PITCH_RATED)

    net_thrust = bem_def.thrust / bem_rigid.thrust - 1.0
    net_power = bem_def.power / bem_rigid.power - 1.0
    twist_thrust = bem_twist.thrust / bem_base.thrust - 1.0
    twist_power = bem_twist.power / bem_base.power - 1.0
    print(f"\none-way de-loading (Zhou et al. Table 6 flexible vs rigid: "
          f"{ZHOU_FLEX_THRUST_DELTA:+.2%} thrust, {ZHOU_FLEX_POWER_DELTA:+.2%} power)")
    print(f"  elastic twist alone (reference radii): thrust {twist_thrust:+.2%}, "
          f"power {twist_power:+.2%}")
    print(f"  full production path (deformed radii): thrust {net_thrust:+.2%}, "
          f"power {net_power:+.2%} (reported: the radius growth re-loads)")

    assert twist_thrust < 0.0, "the elastic twist must reduce thrust"
    assert twist_power < 0.0, "the elastic twist must reduce power"
    assert net_thrust < 0.0, (
        f"the deformed blade does not de-load in thrust on the production path "
        f"({net_thrust:+.2%})"
    )
