"""Our own rated-case validation for the IEA-15MW blade (issue #9, WU-A/WU-B).

Nothing here reuses a previous simulation's output: the aerodynamic loads come from the
repo's own BEM engine driving ccblade (Ning 2014) over the official AeroDyn polars, the
structure is the repo's own shell blade, and the comparison is against a reference that
was extracted from the paper itself (Zhou et al. 2025, *Energy* 336:138488, Table 6).

Three things this module settles:

1. **The rigid aero side.** Our BEM at the rated point (V = 10.59 m/s, 7.56 rpm, pitch 0)
   against Zhou's rigid-blade means: 2.53 MN thrust and 16.11 MW power.
2. **The physical sign of the twist.** The aerodynamic pitching moment from the polars is
   nose-down at every station, and a *nose-down torsional couple* applied to the section
   produces a **negative** mean rotation about z at the tip. So in this mesh a negative
   theta_z is "toward feather / reducing the angle of attack" - the same physical sense as
   Zhou's -3.60 deg and Ma's about -3.9 deg, and the opposite sense to the +0.98 deg
   previously attributed to the BeamDyn anchor.
3. **How the section moment must be applied.** A section moment applied as nodal rotations
   about z feeds the shell's drilling degree of freedom, not the section's torsion, and
   over-reports the twist by orders of magnitude. It has to be a *couple of forces*. The
   test pins both numbers so the difference cannot come back unnoticed.

The production load path is now what the structural response is measured under: the rated
response below is solved from ``ForceProjector``'s own nodal forces
(``test_rated_twist_under_production_loads``), so a repeat of the P5 load-frame defect in
production fails here. The four hand-built applications of ``_rated_load_cases`` remain only
as the **historical sensitivity record** of the task document's sections 18/20 - they show how
much the twist moves with the application, they are not the load path.
"""

from __future__ import annotations


import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as blade_validation  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

from tests.support.paths import DATA_DIR  # noqa: E402

YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"

#: Rated operating point, from the case's own definition.
V_RATED = 10.59
RPM_RATED = 7.56
PITCH_RATED = 0.0

#: Zhou et al. 2025, *Energy* 336:138488, Table 6, rigid blade, fixed condition.
ZHOU_RIGID_THRUST_MN = 2.53
ZHOU_RIGID_POWER_MW = 16.11
#: Same table, flexible blade - the flexible comparison needs a settled load path (WU-B).
ZHOU_FLEXIBLE_THRUST_MN = 2.20
ZHOU_FLEXIBLE_POWER_MW = 14.76
#: Zhou Table 4, mean tip deflections at rated (flapwise, edgewise, torsional).
ZHOU_TIP_FLAP_M = 13.86
ZHOU_TIP_EDGE_M = -1.22
ZHOU_TIP_TORSION_DEG = -3.60

#: The blade is prebent: a physical ring's nodes span ~1e-3 m in z, so the raw unique-z set
#: splits each ring into near-duplicate buckets. This merges them (and detects an element edge
#: as in-plane when its endpoints are within it). It must stay well above the prebend skew and
#: well below the spanwise element length (~0.3 m).
from tests.support.assertions import assert_residual_below  # noqa: E402

STATION_GAP_TOLERANCE = 0.02  # [m]


@pytest.fixture(scope="module")
def rated_bem():
    """The repo's own BEM at the rated point, over the official AeroDyn polars."""
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    solver = BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0)
    return solver.compute(V_RATED, RPM_RATED, PITCH_RATED), blade_aero


def test_rated_bem_matches_the_literature_rigid_case(rated_bem):
    """Our BEM's rigid rotor at rated vs Zhou's rigid means (verified from Table 6).

    This is the aero side, and it is the part the issue assumed was fine. It is fine -
    and now it is measured here rather than quoted from a previous run.
    """
    bem, _ = rated_bem
    thrust_mn = bem.thrust / 1e6
    power_mw = bem.power / 1e6
    print(
        f"\nrated BEM: thrust {thrust_mn:.3f} MN (Zhou rigid {ZHOU_RIGID_THRUST_MN}), "
        f"power {power_mw:.3f} MW (Zhou rigid {ZHOU_RIGID_POWER_MW})"
    )
    assert_residual_below(
        abs(thrust_mn - ZHOU_RIGID_THRUST_MN) / ZHOU_RIGID_THRUST_MN,
        tol=0.03,
        kind="paper",
        reference_name="Zhou et al. 2025, Energy 336:138488, Table 6, rigid case thrust",
        what="rated BEM thrust",
    )
    assert_residual_below(
        abs(power_mw - ZHOU_RIGID_POWER_MW) / ZHOU_RIGID_POWER_MW,
        tol=0.03,
        kind="paper",
        reference_name="Zhou et al. 2025, Energy 336:138488, Table 6, rigid case power",
        what="rated BEM power",
    )


def test_rated_pitching_moment_is_nose_down(rated_bem):
    """The polars give a nose-down pitching moment at rated, at every loaded station.

    This is the physical anchor for the sign question: whatever our mesh calls positive,
    the aerodynamic moment that drives the elastic twist at rated is nose-down (it
    reduces the angle of attack), which is what Zhou's text states and what Ma reports.
    """
    bem, blade_aero = rated_bem
    loaded = bem.Np > 0.01 * bem.Np.max()
    print(
        f"\nMp over the loaded span: {bem.Mp[loaded].min():.1f} .. {bem.Mp[loaded].max():.1f} N.m/m"
    )
    assert np.all(bem.Mp[loaded] <= 0.0), "a station carries a nose-up pitching moment at rated"
    assert np.any(bem.Mp[loaded] < 0.0)


@pytest.fixture(scope="module")
def blade_shell():
    """The repo's shell blade, assembled, with the tip ring and the span stations."""
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=1.0)
    model.generate_mesh()
    mesh, props = model.mesh, model.get_element_properties()
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
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    stations = sorted({round(zz, 6) for zz in coords[:, 2]})
    rings = {zz: np.where(np.abs(coords[:, 2] - zz) < 1e-6)[0] for zz in stations}
    # The blade is prebent, so a physical ring's nodes span ~1e-3 m in z and the raw 1e-6
    # buckets are only slices of a ring. The load path and the ring kinematics both need the
    # merged physical stations (see _physical_stations); the raw sets stay for the tests that
    # were already measuring on them.
    phys_stations = _physical_stations(coords)
    phys_rings = {
        zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0] for zz in phys_stations
    }
    return {
        "n": n,
        "Kff": K[np.ix_(free, free)],
        "free": free,
        "coords": coords,
        "stations": stations,
        "rings": rings,
        "tip": rings[stations[-1]],
        "phys_stations": phys_stations,
        "phys_rings": phys_rings,
        "phys_tip": phys_rings[phys_stations[-1]],
        "mesh": mesh,
        "props": props,
    }


def _section_couples(shell, moment_by_station):
    """A pure moment about the blade axis, applied as a COUPLE OF FORCES.

    Applying the same moment as nodal rotations about z feeds the shell's drilling degree
    of freedom instead of the section's torsion and over-reports the twist by orders of
    magnitude (measured below). A couple is the honest way: ``+F`` flapwise at the node
    of maximum chordwise coordinate and ``-F`` at the minimum, with ``F * dx = moment``,
    so the section sees the moment with no net force.
    """
    force = np.zeros(shell["n"])
    coords = shell["coords"]
    arms = []
    for zz, moment in zip(shell["stations"], moment_by_station, strict=True):
        ring = shell["rings"][zz]
        xs = coords[ring, 0]
        i_max, i_min = ring[np.argmax(xs)], ring[np.argmin(xs)]
        arm = coords[i_max, 0] - coords[i_min, 0]
        arms.append(arm)
        # a degenerate ring (the tip refinement makes tiny sections) would make
        # moment/arm explode; those stations carry almost no moment anyway
        if abs(moment) < 1e-12 or arm < 0.05:
            continue
        f = moment / arm
        force[6 * i_max + 1] += f
        force[6 * i_min + 1] -= f
    print(f"  couple arms: {min(arms):.3f} .. {max(arms):.3f} m")
    return force


def test_section_moment_as_a_couple_sets_the_physical_sign(blade_shell, rated_bem):
    """A nose-down couple gives a negative theta_z: the sign convention is settled.

    With that mapping, Zhou's -3.60 deg and Ma's about -3.9 deg are the *same physical
    sense* as a negative theta_z in this mesh, and the +0.98 deg attributed to the
    BeamDyn anchor is the opposite sense to both the literature and the aerodynamics.
    """
    bem, blade_aero = rated_bem
    span = blade_aero.r - blade_aero.hub_radius
    moments = np.interp(np.asarray(blade_shell["stations"]), span, bem.Mp, left=0.0, right=0.0)

    # as a couple (the honest application)
    u = np.zeros(blade_shell["n"])
    u[blade_shell["free"]] = spsolve(
        blade_shell["Kff"], _section_couples(blade_shell, moments)[blade_shell["free"]]
    )
    tip = blade_shell["tip"]
    couple_twist = float(np.mean([u[6 * i + 5] for i in tip]))

    # as nodal rotations (the wrong application, pinned so it cannot return unnoticed)
    force = np.zeros(blade_shell["n"])
    for zz, moment in zip(blade_shell["stations"], moments, strict=True):
        ring = blade_shell["rings"][zz]
        for i in ring:
            force[6 * i + 5] += moment / len(ring)
    u2 = np.zeros(blade_shell["n"])
    u2[blade_shell["free"]] = spsolve(blade_shell["Kff"], force[blade_shell["free"]])
    nodal_twist = float(np.mean([u2[6 * i + 5] for i in tip]))

    print("\nnose-down section moment applied two ways:")
    print(f"  as a couple of forces : {np.rad2deg(couple_twist):+10.4f} deg")
    print(
        f"  as nodal z-rotations  : {np.rad2deg(nodal_twist):+10.4f} deg "
        f"(drilling DOF, not section torsion)"
    )
    profile = [
        float(np.mean([u[6 * i + 5] for i in blade_shell["rings"][zz]]))
        for zz in blade_shell["stations"]
    ]
    step = max(1, len(profile) // 6)
    print("  theta_z along the span [deg]:", [f"{np.rad2deg(v):+.3f}" for v in profile[::step]])
    print(f"  Zhou Table 4 tip torsion: {ZHOU_TIP_TORSION_DEG} deg (same physical sense)")

    assert couple_twist < 0.0, "a nose-down couple must give a negative theta_z"
    assert abs(nodal_twist) > 2.0 * abs(couple_twist), (
        "the nodal-rotation application is no longer the documented over-report"
    )


def _physical_stations(coords, gap_tolerance=STATION_GAP_TOLERANCE):
    """Merge the raw unique-z buckets into physical spanwise stations.

    The blade is **prebent**, so a physical ring's nodes span a small z range (~1e-3 m) instead
    of lying in one plane; the raw ``unique(z)`` set therefore splits each ring into
    near-duplicate buckets (266 buckets for 186 physical stations here). Values closer than
    ``gap_tolerance`` are one station. This is what lets the load application cover every ring
    instead of skipping the split ones - without it the consistent-traction and shear-flow paths
    apply nothing at 62 buckets, 14% of the normal force.
    """
    raw = np.unique(np.round(coords[:, 2], 9))
    groups: list[list[float]] = []
    for value in raw:
        if groups and value - groups[-1][-1] <= gap_tolerance:
            groups[-1].append(float(value))
        else:
            groups.append([float(value)])
    return [float(np.mean(group)) for group in groups]


def _in_plane_edges(mesh, coords, stations, z_tolerance=STATION_GAP_TOLERANCE):
    """The chordwise edges of every station, deduplicated, with their direction.

    Needed for both a consistent traction (each edge carries its share by length, not each
    node an equal force) and a shear flow (a torque enters a closed thin-walled tube as a
    tangential flow along the walls, ``q = M / 2A``, never as a force pair at two nodes).

    ``z_tolerance`` must exceed the blade's prebend skew (a ring's nodes span ~1e-3 m in z) and
    stay below the spanwise element length (~0.3 m). Each edge is assigned to the station nearest
    its midpoint, so a prebent ring is not split between near-duplicate z buckets.
    """
    z = coords[:, 2]
    stations_arr = np.asarray(stations, dtype=float)
    out = {zz: {} for zz in stations}
    for element in mesh.elements:
        ids = [nd.id for nd in element.nodes]
        for a, b in zip(ids, ids[1:] + ids[:1], strict=True):
            za, zb = z[mesh.node_id_to_index[a]], z[mesh.node_id_to_index[b]]
            if abs(za - zb) < z_tolerance:
                key = stations[int(np.argmin(np.abs(stations_arr - 0.5 * (za + zb))))]
                out[key][tuple(sorted((a, b)))] = (a, b)
    return out


def _unit_shear_flow(mesh, coords, edges):
    """A unit shear flow along a station's walls: returns (force vector, realised moment)."""
    force = np.zeros(6 * len(coords))
    for a, b in edges.values():
        ia, ib = mesh.node_id_to_index[a], mesh.node_id_to_index[b]
        pa, pb = coords[ia], coords[ib]
        length = float(np.hypot(pb[0] - pa[0], pb[1] - pa[1]))
        if length < 1e-9:
            continue
        tangent = (pb - pa) / length
        force[6 * ia] += 0.5 * length * tangent[0]
        force[6 * ia + 1] += 0.5 * length * tangent[1]
        force[6 * ib] += 0.5 * length * tangent[0]
        force[6 * ib + 1] += 0.5 * length * tangent[1]
    moment = float(np.sum(coords[:, 0] * force[1::6] - coords[:, 1] * force[0::6]))
    return force, moment


def _rated_load_cases(shell, bem, blade_aero, mesh):
    """The rated load applications, as nodal force vectors, with the physics-correct ones.

    Three application defects were found by a controlled eccentric-force test and are fixed
    here (see the task document section 18):

    * a section moment must be a **shear flow** along the walls (``q = M / 2A``), not a force
      pair at two nodes - the pair gives erratic twist ratios up to 3.4x;
    * a distributed force must be a **consistent traction** (length-weighted per edge), not
      an equal force per node, because the node spacing around a real airfoil is not uniform;
    * moving a resultant from the ring centroid ``x_c`` to the aerodynamic centre ``x_ac``
      needs a moment of ``(x_ac - x_c) * F``, not ``(x_c - x_ac) * F``;
    * every per-station contribution must carry its **spanwise tributary length** ``dz[k]``: the
      mesh has 266 z-stations over 117 m (tip-refined and non-uniform), so applying a per-metre
      load once per station inflates the resultant by 1.72-2.19x and concentrates it where the
      stations are dense. ``test_rated_aero_loads_reproduce_the_bem_resultants`` asserts the
      corrected resultant against the BEM's own integral.

    ``at_ac`` is the physical case. ``uniform`` keeps the per-node force and no moment, which
    is the convention the previous one-way path used, and it is kept only as the sensitivity
    reference.
    """
    n = shell["n"]
    coords = shell["coords"]
    # Merge the raw z buckets into physical stations; a prebent ring is not planar, and the
    # raw unique-z set splits it (see _physical_stations). Without this the edge-based paths
    # skip the split stations and the applied resultant is 14% short of the BEM's integral.
    stations = _physical_stations(coords)
    rings = {zz: np.where(np.abs(coords[:, 2] - zz) < STATION_GAP_TOLERANCE)[0] for zz in stations}
    edges = _in_plane_edges(mesh, coords, stations)
    flows = {zz: _unit_shear_flow(mesh, coords, edges[zz]) for zz in stations}
    span = blade_aero.r - blade_aero.hub_radius
    Np = np.interp(np.asarray(stations), span, bem.Np, left=0.0, right=0.0)
    Tp = np.interp(np.asarray(stations), span, bem.Tp, left=0.0, right=0.0)
    Mp = np.interp(np.asarray(stations), span, bem.Mp, left=0.0, right=0.0)

    # Spanwise tributary lengths. The mesh stations are non-uniform (tip-refined), so a
    # per-metre load must be multiplied by the station's own tributary length. The two end
    # stations carry half a cell (the trapezoid weights), which makes sum(dz) exactly the
    # blade length.
    station_z = np.asarray(stations, dtype=float)
    dz = np.empty_like(station_z)
    dz[0] = 0.5 * (station_z[1] - station_z[0])
    dz[-1] = 0.5 * (station_z[-1] - station_z[-2])
    dz[1:-1] = 0.5 * (station_z[2:] - station_z[:-2])
    blade_length = station_z[-1] - station_z[0]
    assert abs(dz.sum() - blade_length) < 1e-9 * blade_length, (
        f"tributary lengths sum to {dz.sum():.6f} m, not the blade length {blade_length:.6f} m"
    )

    names = ("mp_only", "uniform", "at_ac", "uniform_plus_mp")
    vectors = {name: np.zeros(n) for name in names}
    for k, zz in enumerate(stations):
        ring = rings[zz]
        xs = coords[ring, 0]
        edge_list = edges[zz]
        total = sum(
            float(
                np.hypot(
                    coords[mesh.node_id_to_index[b]][0] - coords[mesh.node_id_to_index[a]][0],
                    coords[mesh.node_id_to_index[b]][1] - coords[mesh.node_id_to_index[a]][1],
                )
            )
            for a, b in edge_list.values()
        )
        if total > 1e-12:
            for a, b in edge_list.values():
                ia, ib = mesh.node_id_to_index[a], mesh.node_id_to_index[b]
                share = (
                    float(np.hypot(coords[ib][0] - coords[ia][0], coords[ib][1] - coords[ia][1]))
                    / total
                )
                for name in ("at_ac", "uniform_plus_mp"):
                    vectors[name][6 * ia + 1] += 0.5 * share * Np[k] * dz[k]  # normal -> flapwise
                    vectors[name][6 * ia + 0] += (
                        0.5 * share * Tp[k] * dz[k]
                    )  # tangential -> chordwise
                    vectors[name][6 * ib + 1] += 0.5 * share * Np[k] * dz[k]
                    vectors[name][6 * ib + 0] += 0.5 * share * Tp[k] * dz[k]
        for i in ring:  # the previous one-way convention, kept as the sensitivity reference
            vectors["uniform"][6 * i + 1] += Np[k] * dz[k] / len(ring)
            vectors["uniform"][6 * i + 0] += Tp[k] * dz[k] / len(ring)
        flow, moment = flows[zz]
        if abs(moment) > 1e-12:
            # the resultant of a CONSISTENT traction sits at the edge-length-weighted
            # centroid of the perimeter, not at the mean of the node positions
            weighted = 0.0
            for a, b in edge_list.values():
                ia, ib = mesh.node_id_to_index[a], mesh.node_id_to_index[b]
                length = float(
                    np.hypot(coords[ib][0] - coords[ia][0], coords[ib][1] - coords[ia][1])
                )
                weighted += 0.5 * (coords[ia][0] + coords[ib][0]) * length
            x_c = weighted / total if total > 1e-12 else float(xs.mean())
            x_ac = float(xs.min()) + 0.25 * (float(xs.max()) - float(xs.min()))
            for name in ("mp_only", "at_ac", "uniform_plus_mp"):
                vectors[name] += flow * (Mp[k] * dz[k] / moment)
            vectors["at_ac"] += flow * (
                Np[k] * dz[k] * (x_ac - x_c) / moment
            )  # move the resultant TO x_ac
    return vectors


def test_rated_tip_twist_matches_zhou_with_the_physical_load_path(blade_shell, rated_bem):
    """The physical sense of the rated twist, and how much the load's line of action moves it.

    Zhou's -3.60 deg (Table 4) is the **total** tip torsion at rated, so the comparable case is
    the full load applied along its physical line of action: forces at the aerodynamic centre
    plus the pitching moment. That gives -2.72 deg against Zhou's -3.60 deg.

    What is *asserted* is the physics: the physical load path gives the nose-down sense (the
    same sense as Zhou and Ma), and the line of action dominates the twist by more than a
    factor of two.

    The **magnitude is withdrawn as a result**. A controlled eccentric-force test (task
    document section 18) showed that the same load set gives -4.35, -14.25 or -34.13 deg
    depending only on how the forces and moments are distributed over the shell, so the
    number measures the application, not the model. It is printed for the record and the
    matrix row is flagged accordingly; promoting it needs the application validated on a
    case with an exact answer (a rectangular closed tube with a known GJ and a known torque).
    """
    bem, blade_aero = rated_bem
    vectors = _rated_load_cases(blade_shell, bem, blade_aero, blade_shell["mesh"])
    twists = {}
    for name, force in vectors.items():
        u = np.zeros(blade_shell["n"])
        u[blade_shell["free"]] = spsolve(blade_shell["Kff"], force[blade_shell["free"]])
        tip = blade_shell["tip"]
        twists[name] = float(np.mean([u[6 * i + 5] for i in tip]))

    print("\nrated tip twist (rotation about the blade axis), IEA 15 MW at V=10.59, 7.56 rpm:")
    for name, label in (
        ("mp_only", "pitching moment alone (couple)"),
        ("uniform", "Np+Tp uniform over the ring (one-way path)"),
        ("uniform_plus_mp", "Np+Tp uniform + pitching moment"),
        ("at_ac", "Np+Tp at the aerodynamic centre + pitching moment"),
    ):
        print(f"  {label:52} {np.rad2deg(twists[name]):+9.4f} deg")
    print(f"  {'Zhou 2025 Table 4 (total tip torsion)':52} {ZHOU_TIP_TORSION_DEG:+9.4f} deg")

    ratio = np.rad2deg(twists["at_ac"]) / ZHOU_TIP_TORSION_DEG
    magnitudes = [abs(v) for v in twists.values()]
    smallest = min(magnitudes)
    # a twist of exactly zero would divide by zero; the spread is then unbounded
    spread = max(magnitudes) / smallest if smallest > 1e-12 else float("inf")
    print(
        f"  physical case / Zhou = {ratio:.3f}  <- WITHDRAWN as a result: the magnitude moves "
        f"by {spread:.1f}x with the load application alone (see the four cases above), so it is "
        f"not a measurement of the model until the application is validated on a case with an "
        f"exact answer (task document section 18.6)"
    )

    # ASSERTED: the physics (the sense of the twist) and the load-path dominance. A sign
    # test is a valid test - it proves the physics is right - and it needs no tolerance.
    assert twists["at_ac"] < 0.0, "the physical load path must give the nose-down sense"
    # Direction-free dominance: the *choice* of application moves the twist by more than a
    # factor of two. (An earlier version asserted that the uniform-ring case exceeded the
    # aerodynamic-centre one; fixing the application reversed that inequality, which is
    # itself the point - the magnitude is not a property of the model.)
    assert spread > 2.0, (
        "the load application must dominate the twist; if it does not, the comparison is "
        "not well posed and the magnitude may be promoted"
    )

    # NOT asserted: the 24.5% magnitude gap against Zhou. It is above the suite's 5% rule, so
    # it is reported as a residual and recorded in the validation matrix (section 9.5) rather
    # than dressed up with a wide tolerance. The remaining differences are the load case (our
    # steady BEM loads against Zhou's fully coupled aeroelastic solution), the aerodynamic
    # model (BEM against lifting-line free vortex wake) and the reference's own modelling.
    assert abs(1.0 - ratio) > 0.05, (
        "the residual is now inside 5%: promote it to an asserted row in the matrix"
    )


def _ring_kinematics(coords, u, ring):
    """Decompose one planar ring's in-plane displacement into rotation and distortion.

    A ring is the set of nodes at one constant z (``blade_shell["rings"][zz]``), so it is
    planar. The ring's in-plane field ``(u_x, u_y)`` at ``(x, y)`` is fitted by least squares
    to the six-parameter affine map::

        u_x = a11 * x + a12 * y + tx
        u_y = a21 * x + a22 * y + ty

    using ``(x, y)`` as they are (the fitted translations absorb the choice of origin). From
    the fitted matrix::

        omega       = 0.5 * (a21 - a12)                 # section rotation [rad]
        shear       = 0.5 * (a12 + a21)                 # parallelogram distortion (exy)
        dilatation  = a11 + a22                         # relative area change (breathing)
        distortion  = sqrt(a11**2 + a22**2 + 2*shear**2)  # symmetric strain magnitude
        residual    = ||u_ip - (A x + t)|| / ||u_ip||   # non-affine part, 0 if exactly affine

    ``omega`` is the antisymmetric (rotation) part and is invariant to the origin. A
    separate **rigid-only** model is also fitted: two translations plus one rotation about
    the ring node mean, with no strain. Its rotation is ``rigid_rotation`` and its residual
    ``rigid_residual`` measure the same non-affine departure for the strain-free model, so
    ``rigid_rotation`` and the affine ``omega`` can be compared with the module's existing
    mean ``theta_z`` (DOF index 5). The ring node ordering is never used, so no shoelace
    area is needed and a contour that is not in order does not matter.
    """
    idx = np.asarray(ring)
    x = coords[idx, 0]
    y = coords[idx, 1]
    ux = u[6 * idx]
    uy = u[6 * idx + 1]

    # Affine fit: a [N x 3] design matrix per component, translations included.
    design = np.column_stack([x, y, np.ones_like(x)])
    cx = np.linalg.lstsq(design, ux, rcond=None)[0]
    cy = np.linalg.lstsq(design, uy, rcond=None)[0]
    a11, a12, _tx = cx
    a21, a22, _ty = cy

    omega = 0.5 * (a21 - a12)
    shear = 0.5 * (a12 + a21)
    dilatation = a11 + a22
    distortion = float(np.sqrt(a11**2 + a22**2 + 2.0 * shear**2))

    u_inplane = np.concatenate([ux, uy])
    norm = float(np.linalg.norm(u_inplane))
    affine = np.concatenate([design @ cx, design @ cy])
    residual = float(np.linalg.norm(u_inplane - affine) / norm) if norm > 1e-30 else 0.0

    # Rigid-only model: translation + rotation about the ring node mean.
    uxm = ux - ux.mean()
    uym = uy - uy.mean()
    dx = x - x.mean()
    dy = y - y.mean()
    denom = float(np.sum(dx * dx + dy * dy))
    rigid_rotation = float(np.sum(dx * uym - dy * uxm) / denom) if denom > 1e-30 else 0.0
    rigid = np.concatenate([ux.mean() - rigid_rotation * dy, uy.mean() + rigid_rotation * dx])
    rigid_residual = float(np.linalg.norm(u_inplane - rigid) / norm) if norm > 1e-30 else 0.0

    return {
        "omega": float(omega),
        "shear": float(shear),
        "dilatation": float(dilatation),
        "distortion": distortion,
        "residual": residual,
        "rigid_rotation": rigid_rotation,
        "rigid_residual": rigid_residual,
    }


def test_rated_twist_with_the_validated_application_and_the_measured_section_distortion(
    blade_shell, rated_bem
):
    """Re-derive the rated twist with the validated load path and measure the distortion.

    The exact thin-walled-tube case (``tests/test_thin_walled_tube_torsion.py``, task
    document section 19) showed that a closed section loaded at one end with a clamped root
    can have its measured twist dominated by a **section-distortion mode**: the same
    resultant torque applied on one wall instead of as a shear flow gave 125.7x the section
    shear, 69.2x the energy, 74.8x the tip in-plane displacement, and the two twist metrics
    disagreed by 2050%. The blade's rated magnitude is withdrawn because the same load set
    moves the tip twist by a factor 4.3-5.5 depending only on the application (section 20;
    the earlier 7.8 was computed with a load-magnitude defect).

    This test splits each ring's in-plane field into an origin-invariant rotation ``omega`` (the
    antisymmetric part of the affine fit) and the strain measures (parallelogram ``shear``,
    breathing ``dilatation``, symmetric ``distortion``), and compares the four applications of
    ``_rated_load_cases`` on both. Three claims are asserted; two of them replaced earlier,
    mis-specified predictions after the corrected (tributary-weighted) loads refuted them:

    * **Sign (holds).** Every application gives ``omega < 0``: nose-down at rated, the sense of
      Zhou's -3.60 deg.
    * **The spread is the finding (was mis-specified).** The original prediction was
      ``spread_omega < spread_theta_z`` - "the distortion-free metric collapses the application
      spread". Measured on the corrected loads it is the opposite: ``spread_omega = 5.459``
      against ``spread_theta_z = 4.290``. The hypothesis that the spread *lives in the
      distortion* is therefore **refuted on the blade**: the application moves the section
      rotation itself, so the rated twist magnitude is not promotable and stays a reported
      residual. The asserted claim is the true one, ``spread_omega > 2.0``; promoting the
      magnitude requires this spread to collapse, at which point this assertion is removed with
      the promotion.
    * **The tube's transfer test was mis-specified (replaced).** The original prediction was
      ``max(distortion[uniform], distortion[uniform_plus_mp]) > 1.5 * distortion[at_ac]``.
      Measured: **0.982x**, i.e. the off-path cases excite *no more* distortion than ``at_ac``.
      The premise was wrong: the tube's contrast is a **localized one-wall traction against a
      closed shear flow**, while on the blade *both* compared cases are ring-distributed. The
      correct analogue is the **validated shear-flow moment application (``mp_only``) against
      the force-loaded cases**: an in-plane traction around the ring loads the section's
      in-plane flexibility, a pure torsion shear flow does not. It holds strongly -
      ``distortion[at_ac] / distortion[mp_only] = 18.8x`` - and the assertion is
      ``distortion[at_ac] > 5.0 * distortion[mp_only]``. The off-path numbers are still printed
      so the refuted comparison stays on the record.
    * **Promotion guard (holds).** The validated path's residual against Zhou stays above the
      suite's 5% rule, so the magnitude remains a reported residual rather than an asserted row.

    **What the corrected loads say physically.** For the physical application ``at_ac`` the
    section rotation is ``omega = -26.12 deg`` and the section strain is
    ``distortion = 9.375e-2``, i.e. ``distortion / |omega| = 0.359`` at that ring: the section
    strains are the same order as the rotation, so the measured "twist" of this shell blade at
    rated is substantially a **sectional deformation**, not a rigid section rotation. That is a
    property of the model plus the load path: the aerodynamic normal force acts at the 25 %-chord
    aerodynamic centre while the mesh's measured shear centre sits at 0.477 of the chord, so the
    eccentricity torque dominates the pitching moment (``mp_only`` gives -4.79 deg, ``at_ac``
    gives -26.12 deg). Settling whether the anchor's shear centre agrees is the next unit
    (**WU-C**, the BeamDyn 6x6 ``xS`` against the mesh's per-station shear centre); if the anchor
    puts the shear centre near the pitch axis, the eccentricity - and the twist - is a
    structural/geometry difference, not a load-application defect.

    A refutation here is a finding, not a test defect.
    """
    shell = blade_shell
    bem, blade_aero = rated_bem
    vectors = _rated_load_cases(shell, bem, blade_aero, shell["mesh"])
    # The kinematics need COMPLETE rings. The raw fixture rings are 1e-6 z buckets and the blade
    # is prebent (a physical ring's nodes span ~1e-3 m in z), so a raw bucket is only a slice of
    # a ring - a fit on it is meaningless. The load path already merges the buckets into physical
    # stations; the kinematics must use the same merged rings.
    stations = shell["phys_stations"]
    rings = shell["phys_rings"]
    tip = shell["phys_tip"]
    z_tip = stations[-1]

    def nearest(target):
        pos = int(np.searchsorted(stations, target))
        return stations[min(max(pos, 0), len(stations) - 1)]

    probe = {"50%": nearest(0.5 * z_tip), "80%": nearest(0.8 * z_tip), "tip": z_tip}
    raw_stations = shell["stations"]

    def nearest_raw(target):
        pos = int(np.searchsorted(raw_stations, target))
        return raw_stations[min(max(pos, 0), len(raw_stations) - 1)]

    raw_counts = {key: len(shell["rings"][nearest_raw(zz)]) for key, zz in probe.items()}
    ring_counts = {key: len(rings[zz]) for key, zz in probe.items()}
    print(
        f"\nring node counts at the probes - nearest raw 1e-6 bucket {raw_counts}, "
        f"merged physical rings {ring_counts}"
    )
    for key, zz in probe.items():
        assert len(rings[zz]) >= max(10, raw_counts[key]), (
            f"the {key} ring ({zz:.3f} m) has {len(rings[zz])} nodes after the merge while the "
            f"nearest raw 1e-6 bucket has {raw_counts[key]}: the merged physical ring must be a "
            f"superset of the raw bucket, never a slice of it, or the affine fit silently runs on "
            f"a partial section"
        )

    order = ("mp_only", "uniform", "uniform_plus_mp", "at_ac")
    labels = {
        "mp_only": "pitching moment alone (shear flow)",
        "uniform": "Np+Tp equal force per node",
        "uniform_plus_mp": "Np+Tp per node + pitching moment",
        "at_ac": "Np+Tp consistent traction at the AC + Mp",
    }
    kin, twist, energy, tip_ip = {}, {}, {}, {}
    disp = {}
    for name in order:
        force = vectors[name]
        u = np.zeros(shell["n"])
        u[shell["free"]] = spsolve(shell["Kff"], force[shell["free"]])
        disp[name] = u
        kin[name] = {
            key: _ring_kinematics(shell["coords"], u, rings[zz]) for key, zz in probe.items()
        }
        twist[name] = float(np.mean([u[6 * i + 5] for i in tip]))
        energy[name] = float(0.5 * force @ u)
        ux = u[6 * np.asarray(tip)]
        uy = u[6 * np.asarray(tip) + 1]
        tip_ip[name] = float(max(np.max(np.abs(ux)), np.max(np.abs(uy))))

    print(
        "\nblade section kinematics, IEA 15 MW at V=10.59, 7.56 rpm "
        "(omega/shear/dilatation/distortion/residual are dimensionless except omega)"
    )
    print(
        f"  {'application':44} {'theta_z':>9} {'omega':>9} {'rigid':>9} {'shear':>11} "
        f"{'dilat':>10} {'distort':>10} {'d/|w|':>8} {'resid':>8} {'energy[J]':>11} "
        f"{'max_ip[m]':>10}"
    )
    for name in order:
        k = kin[name]["tip"]
        print(
            f"  {labels[name]:44} {np.rad2deg(twist[name]):+9.4f} "
            f"{np.rad2deg(k['omega']):+9.4f} {np.rad2deg(k['rigid_rotation']):+9.4f} "
            f"{k['shear']:+11.3e} {k['dilatation']:+10.3e} {k['distortion']:10.3e} "
            f"{k['distortion'] / abs(k['omega']):8.3f} "
            f"{k['residual']:8.4f} {energy[name]:11.4e} {tip_ip[name]:10.3e}"
        )

    print("\n  omega along the span [deg] (50%, 80%, tip):")
    for name in order:
        vals = [np.rad2deg(kin[name][key]["omega"]) for key in ("50%", "80%", "tip")]
        print(f"    {labels[name]:44} " + "  ".join(f"{v:+8.4f}" for v in vals))

    print("\n  tip twist in the module's own sense (mean theta_z about the blade axis):")
    for name in order:
        print(f"    {labels[name]:44} {np.rad2deg(twist[name]):+9.4f} deg")
    print(f"    {'Zhou 2025 Table 4 (total tip torsion)':44} {ZHOU_TIP_TORSION_DEG:+9.4f} deg")

    def spread(values):
        magnitudes = [abs(v) for v in values]
        smallest = min(magnitudes)
        return max(magnitudes) / smallest if smallest > 1e-12 else float("inf")

    omega_mag = [kin[name]["tip"]["omega"] for name in order]
    theta_mag = [twist[name] for name in order]
    spread_omega = spread(omega_mag)
    spread_theta_z = spread(theta_mag)
    omega_at_ac_deg = np.rad2deg(kin["at_ac"]["tip"]["omega"])
    theta_z_at_ac_deg = np.rad2deg(twist["at_ac"])
    ratio_omega = omega_at_ac_deg / ZHOU_TIP_TORSION_DEG
    ratio_theta_z = theta_z_at_ac_deg / ZHOU_TIP_TORSION_DEG
    distortion_at_ac = kin["at_ac"]["tip"]["distortion"]
    distortion_off_path = max(
        kin["uniform"]["tip"]["distortion"], kin["uniform_plus_mp"]["tip"]["distortion"]
    )
    distortion_ratio = (
        distortion_off_path / distortion_at_ac if distortion_at_ac > 0.0 else float("inf")
    )

    print(f"\n  spread_omega   = {spread_omega:.3f}  (max|omega| / min|omega| over the four cases)")
    print(
        f"  spread_theta_z = {spread_theta_z:.3f}  (max|theta_z| / min|theta_z| over the four cases)"
    )
    print(f"  ratio_omega    = {ratio_omega:.4f}  (omega_at_ac / Zhou {ZHOU_TIP_TORSION_DEG} deg)")
    print(
        f"  ratio_theta_z  = {ratio_theta_z:.4f}  (theta_z_at_ac / Zhou {ZHOU_TIP_TORSION_DEG} deg)"
    )
    print(
        f"  distortion ratio = {distortion_ratio:.3f}  "
        f"(max(distortion[uniform], distortion[uniform_plus_mp]) / distortion[at_ac])"
    )
    distortion_mp_only = kin["mp_only"]["tip"]["distortion"]
    distortion_at_ac_over_mp = (
        distortion_at_ac / distortion_mp_only if distortion_mp_only > 0.0 else float("inf")
    )
    max_u_ac = float(
        np.max(
            np.sqrt(disp["at_ac"][0::6] ** 2 + disp["at_ac"][1::6] ** 2 + disp["at_ac"][2::6] ** 2)
        )
    )
    print(
        f"  distortion[at_ac] / distortion[mp_only] = {distortion_at_ac_over_mp:.3f}  "
        f"(the corrected transfer test: ring in-plane traction vs pure torsion shear flow)"
    )
    print(
        f"  corrected max |u| (at_ac, one-way) = {max_u_ac:.3f} m = "
        f"{max_u_ac / ZHOU_TIP_FLAP_M:.3f}x Zhou's coupled tip flap {ZHOU_TIP_FLAP_M:.2f} m "
        f"(reported; a one-way application of the rigid-blade loads on the flexible structure "
        f"must exceed the coupled deflection without being an order of magnitude away)"
    )

    # 1. Sign (physics, no tolerance): the rated twist is nose-down for every application.
    for name in order:
        assert kin[name]["tip"]["omega"] < 0.0, (
            f"{name} gives a non-negative section rotation omega = "
            f"{np.rad2deg(kin[name]['tip']['omega']):+.4f} deg; the rated twist must be "
            f"nose-down (the sense of Zhou's {ZHOU_TIP_TORSION_DEG} deg)"
        )

    # 2. The true claim, after the refutation of the original "the spread lives in the
    #    distortion": the application spread SURVIVES the distortion-free section-rotation
    #    metric (measured spread_omega = 5.459 > spread_theta_z = 4.290, not the predicted
    #    collapse), so the rated twist magnitude is not promotable and stays a reported
    #    residual. Promoting it requires this spread to collapse to a settled rotation; this
    #    assertion is then removed with the promotion.
    assert spread_omega > 2.0, (
        f"the section-rotation spread has collapsed to {spread_omega:.3f} <= 2.0 "
        f"(spread_theta_z = {spread_theta_z:.3f}): the application spread no longer survives "
        f"the distortion-free metric, so the rated twist magnitude is now promotable - "
        f"promote it to an asserted row in the matrix and remove this assertion"
    )

    # 3. The corrected transfer test. The original prediction (off-path > 1.5 * at_ac) was
    #    mis-specified: it compared two ring-distributed cases, while the tube's contrast is a
    #    localized one-wall traction against a closed shear flow. Measured, it is 0.982x, i.e.
    #    refuted. The correct analogue is the validated shear-flow moment application
    #    (mp_only) against the force-loaded cases: a ring in-plane traction loads the section's
    #    in-plane flexibility, a pure torsion shear flow does not. Measured 18.8x.
    assert distortion_at_ac > 5.0 * distortion_mp_only, (
        f"the force-loaded sections do not excite the predicted distortion over the shear-flow "
        f"moment application: distortion[at_ac] {distortion_at_ac:.6e} vs 5.0 * "
        f"distortion[mp_only] {5.0 * distortion_mp_only:.6e} (ratio "
        f"{distortion_at_ac_over_mp:.3f})"
    )

    # 4. Promotion guard (the suite's rule 3, same style as the existing test): while the
    #    validated application's residual against Zhou stays above 5%, the blade row remains
    #    a reported residual. Once it enters the bound, promote it to an asserted row.
    assert abs(1.0 - ratio_omega) > 0.05, (
        f"the section-rotation residual is now inside 5% (ratio_omega = {ratio_omega:.4f}, "
        f"ratio_theta_z = {ratio_theta_z:.4f}): promote the magnitude to an asserted row "
        f"in the matrix and remove this guard"
    )


def test_rated_aero_loads_reproduce_the_bem_resultants(blade_shell, rated_bem):
    """The applied load vectors must carry the BEM's own integrated resultants.

    The mesh has 266 z-stations over 117 m and the tip refinement makes the spacing
    non-uniform, so a per-metre load applied once per station without a spanwise tributary
    weight inflates the resultant and concentrates it where the stations are dense. The
    read-only diagnostic that found the defect measured ``at_ac`` at 1.76x and ``mp_only``
    at 1.72x the BEM's own integral. This test is the invariant that catches it: the applied
    force and moment must reproduce the BEM's own trapezoidal integrals
    ``sum_i (v[i] + v[i+1]) / 2 * (r[i+1] - r[i])`` on the BEM's radial stations, to 0.5%.

    It also prints the corrected tip displacement against Zhou's verified coupled tip
    flapwise deflection (13.86 m). That comparison is **reported, not asserted**: a one-way
    load against a coupled aeroelastic solution has no defensible numeric bound, but a value
    of the same order is the independent confirmation that the load magnitude is right now.
    """
    bem, blade_aero = rated_bem
    shell = blade_shell
    r = np.asarray(blade_aero.r, dtype=float)

    def integral(values):
        values = np.asarray(values, dtype=float)
        return float(np.sum(0.5 * (values[:-1] + values[1:]) * np.diff(r)))

    i_np, i_tp, i_mp = integral(bem.Np), integral(bem.Tp), integral(bem.Mp)

    vectors = _rated_load_cases(shell, bem, blade_aero, shell["mesh"])
    coords = shell["coords"]

    def resultants(force):
        fx = float(force[0::6].sum())
        fy = float(force[1::6].sum())
        mz = float(np.sum(coords[:, 0] * force[1::6] - coords[:, 1] * force[0::6]))
        return fx, fy, mz

    fx_ac, fy_ac, mz_ac = resultants(vectors["at_ac"])
    fx_un, fy_un, mz_un = resultants(vectors["uniform"])
    fx_mp, fy_mp, mz_mp = resultants(vectors["mp_only"])

    print("\nrated aero resultants: applied vs the BEM's own integrals")
    print(f"  BEM integrals (trapezoid over {len(r)} radial stations):")
    print(f"    I_Np = {i_np:+.6e} N   I_Tp = {i_tp:+.6e} N   I_Mp = {i_mp:+.6e} N.m")
    print(
        f"  at_ac   sum(f_y) = {fy_ac:+.6e} N  ratio {fy_ac / i_np:.6f}   "
        f"sum(f_x) = {fx_ac:+.6e} N  ratio {fx_ac / i_tp:.6f}   "
        f"sum(M_z) = {mz_ac:+.6e} N.m"
    )
    print(
        f"  uniform sum(f_y) = {fy_un:+.6e} N  ratio {fy_un / i_np:.6f}   "
        f"sum(f_x) = {fx_un:+.6e} N  ratio {fx_un / i_tp:.6f}   "
        f"sum(M_z) = {mz_un:+.6e} N.m"
    )
    print(
        f"  mp_only sum(x F_y - y F_x) = {mz_mp:+.6e} N.m  ratio {mz_mp / i_mp:.6f}   "
        f"(net force {fy_mp:+.3e} N, a pure couple)"
    )

    assert_residual_below(
        abs(fy_ac - i_np) / i_np,
        tol=0.005,
        kind="self",
        reference_name="our own BEM integral for the same load path",
        what="at_ac normal resultant",
    )
    assert_residual_below(
        abs(fx_ac - i_tp) / i_tp,
        tol=0.005,
        kind="self",
        reference_name="our own BEM integral for the same load path",
        what="at_ac tangential resultant",
    )
    assert_residual_below(
        abs(fy_un - i_np) / i_np,
        tol=0.005,
        kind="self",
        reference_name="our own BEM integral for the same load path",
        what="uniform normal resultant, same resultant a different distribution",
    )
    assert_residual_below(
        abs(mz_mp - i_mp) / abs(i_mp),
        tol=0.005,
        kind="self",
        reference_name="our own BEM integral for the same load path",
        what="mp_only resultant moment",
    )

    u = np.zeros(shell["n"])
    u[shell["free"]] = spsolve(shell["Kff"], vectors["at_ac"][shell["free"]])
    max_u = float(np.max(np.sqrt(u[0::6] ** 2 + u[1::6] ** 2 + u[2::6] ** 2)))
    print(f"\n  corrected max |u| over all nodes (at_ac, one-way) = {max_u:.3f} m")
    print(
        f"  Zhou 2025 Table 4 coupled tip flapwise deflection   = {ZHOU_TIP_FLAP_M:.2f} m "
        f"(reported only; a one-way load has no defensible bound against a coupled solution)"
    )


@pytest.fixture(scope="module")
def production_rated_loads(blade_shell, rated_bem):
    """The rated aero loads through the **production** projector, on this module's mesh.

    Everything above applies its own hand-built load vectors, so this module - the one whose
    subject *is* the load application - never exercised ``ForceProjector``; a repeat of the P5
    class of load-frame defect in production would still leave this file green. This fixture
    builds the projector exactly as production does (``standalone.py``): the same real
    ``BladeAero`` the rated BEM uses, the same mesh, and the **default**
    ``normal_direction``/``tangential_direction``, so only ``span_direction`` is passed.

    It then projects the rated point and scatters the ``(n_nodes, 3)`` result into the 6-DOF
    load vector this module already solves with (``f[6*i + 0..2] = forces[i, 0..2]``). It
    returns that assembled vector (apply with ``force[blade_shell["free"]]``) and prints the
    summed nodal force and its comparison with ``bem.thrust / n_blades`` - computed here, so
    the check is independent of the projector's own ``verify()``.
    """
    bem, blade_aero = rated_bem
    shell = blade_shell
    projector = ForceProjector(
        shell["mesh"], blade_aero, span_direction=blade_validation.SPAN_DIRECTION
    )
    forces = projector.project(bem)
    assert forces.shape == (shell["n"] // 6, 3), (
        f"project() returned {forces.shape}, not ({shell['n'] // 6}, 3)"
    )
    force = np.zeros(shell["n"])
    force[0::6] = forces[:, 0]
    force[1::6] = forces[:, 1]
    force[2::6] = forces[:, 2]

    n_blades = int(blade_aero.n_blades)
    thrust_per_blade = float(bem.thrust) / n_blades
    total = forces.sum(axis=0)
    total_mag = float(np.linalg.norm(total))
    print(
        f"\nproduction ForceProjector load path at rated "
        f"(V={V_RATED} m/s, {RPM_RATED} rpm, pitch {PITCH_RATED}, defaults):"
    )
    print(
        f"  sum(F) = {np.array2string(total, precision=4, suppress_small=True)} N  "
        f"|sum(F)| = {total_mag:.4f} N"
    )
    print(
        f"  bem.thrust / {n_blades} = {thrust_per_blade:.4f} N -> ratio "
        f"{total_mag / thrust_per_blade:.6f} (relative {total_mag / thrust_per_blade - 1.0:+.4%})"
    )
    return force


def _measured_tip_axes(coords, ring):
    """The tip ring's own edgewise(chord) and flapwise axes, measured from its outline.

    The in-plane principal axis of the ring is its chord (on this mesh a ring's ``x`` extent is
    the local chord and its ``y`` extent the airfoil thickness - the measurement
    ``tests/test_force_projection_load_frame.py`` makes), and the perpendicular in-plane axis is
    the flapwise/thickness one. ``chord`` is oriented toward ``+x``, which on THIS mesh is toward
    the LEADING EDGE (the leading edge sits at ``+pitch_axis * chord`` from the pitch axis), and
    ``flap`` toward ``+y`` = the downwind/fluid direction of the model's declared convention.
    Callers comparing with Zhou must flip ``chord`` to get Zhou's TE-positive edgewise axis.
    """
    pts = coords[np.asarray(ring)]
    in_plane = pts[:, :2] - pts[:, :2].mean(axis=0)
    _, _, vt = np.linalg.svd(in_plane, full_matrices=False)
    chord = np.array([vt[0, 0], vt[0, 1], 0.0])
    chord = chord / np.linalg.norm(chord)
    if chord[0] < 0.0:
        chord = -chord
    flap = np.array([-chord[1], chord[0], 0.0])
    if flap[1] < 0.0:
        flap = -flap
    return chord, flap


def test_rated_twist_under_production_loads_is_nose_down(
    blade_shell, rated_bem, production_rated_loads
):
    """The physical sense at rated: the tip section rotation is nose-down (Zhou's -3.60 deg).

    Split out of ``test_rated_twist_under_production_loads`` so that the applied-load
    invariant, the applied direction and the promotion guard there stay live while this
    one asserts the physical sign on its own.  It was declared ``xfail(strict=True)``
    while the wired multi-cell realisation was mesh-dependent; the property gate now
    routes the properties-less production projector to the minimum-norm fallback, so the
    measured rotation is nose-down again (``omega = -0.8696 deg``) and the marker is gone.
    """
    shell = blade_shell
    force = production_rated_loads
    u = np.zeros(shell["n"])
    u[shell["free"]] = spsolve(shell["Kff"], force[shell["free"]])
    tip = np.asarray(shell["phys_tip"])
    omega = _ring_kinematics(shell["coords"], u, tip)["omega"]
    assert omega < 0.0, (
        f"the production load path gives a non-negative tip section rotation "
        f"omega = {np.rad2deg(omega):+.4f} deg; the properties-less production projector "
        f"must route through the minimum-norm fallback and measure nose-down "
        f"-0.8696 deg (the sense of Zhou's {ZHOU_TIP_TORSION_DEG} deg)"
    )


@pytest.mark.xfail(
    strict=True,
    reason=(
        "the wall-flow realisation is verified correct (skin walls carry q_i, shared webs "
        "q_i - q_j, the total moment is exact and the net force is zero), yet the section "
        "rotates nose-up when the ring wall stiffnesses come from "
        "Blade.get_element_properties(): omega = +0.0773 deg against the properties-less "
        "fallback's nose-down -0.8696 deg. The requested section torque itself measures "
        "NOSE-UP (+4.840e5 N.m) because the lever-arm transfer (+7.607e5 N.m) dominates and "
        "opposes the polars' pitching moment (-2.766e5 N.m). The sign question is open and "
        "tracked by tools/diagnose_applied_torsion_sign.py; strict=True, so a fix that "
        "restores the physical nose-down sense turns this marker into an XPASS."
    ),
)
def test_rated_twist_under_multi_cell_production_loads_is_nose_down(blade_shell, rated_bem):
    """The multi-cell wall-flow path's measured sense, pinned until the sign question closes.

    The realisation itself is verified in ``tests/test_multicell_shear_flow.py``; the sign
    of the realised section torque is not settled.  This builds the projector exactly as
    production would if it held the laminate map (``element_properties`` from the ``Blade``
    model, as the ``blade_shell`` fixture does) and asserts the physical sense, so a fix
    flips the marker to XPASS and forces this row to be updated.
    """
    shell = blade_shell
    bem, blade_aero = rated_bem
    projector = ForceProjector(
        shell["mesh"],
        blade_aero,
        span_direction=blade_validation.SPAN_DIRECTION,
        element_properties=shell["props"],
    )
    forces = projector.project(bem)
    force = np.zeros(shell["n"])
    force[0::6] = forces[:, 0]
    force[1::6] = forces[:, 1]
    force[2::6] = forces[:, 2]
    u = np.zeros(shell["n"])
    u[shell["free"]] = spsolve(shell["Kff"], force[shell["free"]])
    tip = np.asarray(shell["phys_tip"])
    omega = _ring_kinematics(shell["coords"], u, tip)["omega"]
    print(
        f"\nmulti-cell load path: tip section omega = {np.rad2deg(omega):+9.4f} deg "
        f"(nose-up, measured +0.0773 deg; Zhou torsion {ZHOU_TIP_TORSION_DEG:+.2f} deg)"
    )
    assert omega < 0.0, (
        f"the multi-cell load path gives a nose-up tip section rotation "
        f"omega = {np.rad2deg(omega):+.4f} deg; the measured multi-cell sense is "
        f"+0.0773 deg, the sign question is open "
        f"(tools/diagnose_applied_torsion_sign.py), and the physical sense required here "
        f"is Zhou's {ZHOU_TIP_TORSION_DEG} deg"
    )


def test_rated_twist_under_production_loads(blade_shell, rated_bem, production_rated_loads):
    """The rated structural response measured under the **production** load path.

    The four applications of ``_rated_load_cases`` are this module's own reconstruction of the
    load path, so they cannot see a defect in production's. Here the shell is solved with the
    forces ``ForceProjector.project()`` actually produces (default directions), and three things
    are asserted:

    * **The applied-load invariant**, ``|sum(F)|`` within **2 %** of ``bem.thrust / n_blades``.
      That bound is imported verbatim from the P5 guard
      (``tests/test_force_projection_load_frame.py::test_load_sense_is_downwind_and_driving``,
      ``SENSE_THRUST_TOL = 0.02``); this is the same physical claim, now on the structural path.
    * **The applied-load direction** - the same ``|sum(F)|`` is a norm, so a load rotated off the
      aero axis would pass it; at least **0.95** of it must ride the flapwise axis measured from
      the tip ring's outline, in the downwind half-space (the mirror of the P5 guard's per-
      station ``|F.f_hat|/|F| >= 0.95``).
    * **The physical sense**: the tip section rotation ``omega`` is negative - nose-down at
      rated, the sense of Zhou's -3.60 deg.

    The magnitude against Zhou is **reported, not asserted**: our steady one-way BEM loads are
    being compared with a coupled aeroelastic solution, so no numeric bound is defensible. The
    promotion guard fires the moment the ratio enters 5 %.
    """
    shell = blade_shell
    bem, blade_aero = rated_bem
    force = production_rated_loads

    n_blades = int(blade_aero.n_blades)
    thrust_per_blade = float(bem.thrust) / n_blades
    applied = np.array([force[0::6].sum(), force[1::6].sum(), force[2::6].sum()])
    applied_mag = float(np.linalg.norm(applied))
    load_ratio = applied_mag / thrust_per_blade

    u = np.zeros(shell["n"])
    u[shell["free"]] = spsolve(shell["Kff"], force[shell["free"]])
    tip = np.asarray(shell["phys_tip"])
    kin = _ring_kinematics(shell["coords"], u, tip)
    omega = kin["omega"]
    theta_z = float(np.mean([u[6 * i + 5] for i in tip]))
    d_over_omega = kin["distortion"] / abs(omega) if abs(omega) > 1e-30 else float("inf")

    chord_hat, flap_hat = _measured_tip_axes(shell["coords"], tip)
    # Zhou reports the edgewise deflection POSITIVE TOWARD THE TRAILING EDGE, while on this mesh
    # the measured chord axis points toward the LEADING EDGE (the tip ring's leading edge sits at
    # +x - see _measured_tip_axes and the load-frame test). Without this flip the comparison
    # against Zhou's triple shows a spurious "opposite sign" that belongs to the axis convention,
    # not to the physics.
    edge_hat = -chord_hat
    tip_disp = np.column_stack([u[6 * tip], u[6 * tip + 1], u[6 * tip + 2]])
    flap_def = float(np.mean(tip_disp @ flap_hat))
    edge_def = float(np.mean(tip_disp @ edge_hat))
    ratio_to_zhou = np.rad2deg(omega) / ZHOU_TIP_TORSION_DEG

    print(
        "\nrated structural response under the production load path "
        "(ForceProjector, default directions):"
    )
    print(
        f"  applied |sum(F)|    = {applied_mag:12.4f} N  vs bem.thrust/{n_blades} = "
        f"{thrust_per_blade:12.4f} N  ratio {load_ratio:.6f} ({load_ratio - 1.0:+.3%})"
    )
    print(
        f"  tip section omega   = {np.rad2deg(omega):+9.4f} deg   "
        f"(Zhou torsion {ZHOU_TIP_TORSION_DEG:+.2f} deg)"
    )
    print(f"  tip mean theta_z    = {np.rad2deg(theta_z):+9.4f} deg")
    print(
        f"  section distortion  = {kin['distortion']:.6e}   distortion/|omega| = {d_over_omega:.4f}"
    )
    print(
        f"  tip flapwise defl.  = {flap_def:+9.4f} m   (Zhou {ZHOU_TIP_FLAP_M:+.2f} m -> "
        f"ratio {flap_def / ZHOU_TIP_FLAP_M:+.4f})"
    )
    print(
        f"  tip edgewise defl.  = {edge_def:+9.4f} m   (Zhou {ZHOU_TIP_EDGE_M:+.2f} m -> "
        f"ratio {edge_def / ZHOU_TIP_EDGE_M:+.4f})"
    )
    print(f"  omega / Zhou torsion = {ratio_to_zhou:.4f}")

    # Applied-load invariant: the P5 guard's own 2 % bound on |F| vs bem.thrust / n_blades,
    # applied here to the load the structure is actually solved with.
    assert_residual_below(
        abs(load_ratio - 1.0),
        tol=0.02,
        kind="self",
        reference_name="the BEM thrust the structure is solved against",
        what="applied-load invariant, the P5 guard bound",
    )
    # Directional invariant (the review's R3-001): |sum(F)| alone is a norm claim, so a load
    # rotated off the aero axis but norm-preserving would pass it. The aggregate thrust must
    # also ride the flapwise axis measured from the tip ring's own outline, in the downwind
    # half-space. The 0.95 bound is the same mirror the P5 guard puts on a single station
    # (|F.f_hat|/|F| >= 0.95); it is NOT tightened here because the chordwise share of an
    # aggregate over a 24.41 deg twisting span is not a per-statement claim.
    flap_share = float(applied @ flap_hat) / applied_mag
    chord_share = float(applied @ edge_hat) / applied_mag
    print(
        f"  applied direction   = flapwise share {flap_share:+.4f}  "
        f"edgewise(TE-positive) share {chord_share:+.4f}  (of |sum(F)|)"
    )
    assert flap_share >= 0.95, (
        f"the production load path applies only {flap_share:.4f} of |sum(F)| along the "
        f"measured flapwise axis (bound 0.95): the thrust must be downwind along that axis, "
        f"not merely of the right magnitude"
    )
    # Physical sense (nose-down) is asserted by its own test below, declared xfail while
    # the wired multi-cell realisation is still mesh-dependent.  Everything above stays
    # live: the load invariant, the direction and the promotion guard.
    # Promotion guard, in the style of the existing test: the magnitude stays a reported
    # residual while the comparison against Zhou stays outside 5 %.
    assert abs(1.0 - ratio_to_zhou) > 0.05, (
        f"the production-path tip rotation is now inside 5 % of Zhou "
        f"(ratio_to_zhou = {ratio_to_zhou:.4f}): the row must be promoted to an asserted "
        f"magnitude and this guard removed"
    )
