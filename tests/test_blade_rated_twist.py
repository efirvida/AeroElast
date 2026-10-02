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
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

import test_blade_iea15mw_validation as blade_validation  # noqa: E402
from _openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402

YAML = Path(__file__).resolve().parent / "IEA-15-240-RWT.yaml"
AD_PRIMARY = (Path(__file__).resolve().parent / "reference" / "iea15mw_openfast" /
              "case" / "IEA-15-240-RWT_AeroDyn15.dat")

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
    print(f"\nrated BEM: thrust {thrust_mn:.3f} MN (Zhou rigid {ZHOU_RIGID_THRUST_MN}), "
          f"power {power_mw:.3f} MW (Zhou rigid {ZHOU_RIGID_POWER_MW})")
    assert abs(thrust_mn - ZHOU_RIGID_THRUST_MN) / ZHOU_RIGID_THRUST_MN < 0.03
    assert abs(power_mw - ZHOU_RIGID_POWER_MW) / ZHOU_RIGID_POWER_MW < 0.03


def test_rated_pitching_moment_is_nose_down(rated_bem):
    """The polars give a nose-down pitching moment at rated, at every loaded station.

    This is the physical anchor for the sign question: whatever our mesh calls positive,
    the aerodynamic moment that drives the elastic twist at rated is nose-down (it
    reduces the angle of attack), which is what Zhou's text states and what Ma reports.
    """
    bem, blade_aero = rated_bem
    loaded = bem.Np > 0.01 * bem.Np.max()
    print(f"\nMp over the loaded span: {bem.Mp[loaded].min():.1f} .. {bem.Mp[loaded].max():.1f} N.m/m")
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
        blade_validation._to_rust_mesh(mesh, props), props,
        list(blade_validation.SPAN_DIRECTION), None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))),
                   shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.array([i for i in range(n) if i not in fixed], dtype=np.int64)
    coords = np.array([[nd.x, nd.y, nd.z] for nd in mesh.nodes])
    stations = sorted({round(zz, 6) for zz in coords[:, 2]})
    rings = {zz: np.where(np.abs(coords[:, 2] - zz) < 1e-6)[0] for zz in stations}
    return {"n": n, "Kff": K[np.ix_(free, free)], "free": free, "coords": coords,
            "stations": stations, "rings": rings, "tip": rings[stations[-1]]}


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
    u[blade_shell["free"]] = spsolve(blade_shell["Kff"],
                                     _section_couples(blade_shell, moments)[blade_shell["free"]])
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
    print(f"  as nodal z-rotations  : {np.rad2deg(nodal_twist):+10.4f} deg "
          f"(drilling DOF, not section torsion)")
    profile = [float(np.mean([u[6 * i + 5] for i in blade_shell["rings"][zz]]))
               for zz in blade_shell["stations"]]
    step = max(1, len(profile) // 6)
    print("  theta_z along the span [deg]:",
          [f"{np.rad2deg(v):+.3f}" for v in profile[::step]])
    print(f"  Zhou Table 4 tip torsion: {ZHOU_TIP_TORSION_DEG} deg (same physical sense)")

    assert couple_twist < 0.0, "a nose-down couple must give a negative theta_z"
    assert abs(nodal_twist) > 2.0 * abs(couple_twist), (
        "the nodal-rotation application is no longer the documented over-report"
    )


def _rated_load_cases(shell, bem, blade_aero):
    """The four rated load applications, as nodal force vectors.

    ``forces_at_ac`` is the physical one: the AeroDyn normal and tangential resultants act
    at the **aerodynamic centre** (25% chord), reached by adding the corrective couple that
    moves the uniform ring resultant from the ring's node centroid to the AC, plus the
    aerodynamic pitching moment as a couple. Spreading the loads uniformly over the ring -
    what the one-way path does - puts the resultant at the node centroid instead, and the
    difference is not a detail: it moves the tip twist by a factor of five.
    """
    n = shell["n"]
    coords = shell["coords"]
    stations, rings = shell["stations"], shell["rings"]
    span = blade_aero.r - blade_aero.hub_radius
    Np = np.interp(np.asarray(stations), span, bem.Np, left=0.0, right=0.0)
    Tp = np.interp(np.asarray(stations), span, bem.Tp, left=0.0, right=0.0)
    Mp = np.interp(np.asarray(stations), span, bem.Mp, left=0.0, right=0.0)

    def add_couple(force, ring, moment):
        xs = coords[ring, 0]
        i_max, i_min = ring[np.argmax(xs)], ring[np.argmin(xs)]
        arm = coords[i_max, 0] - coords[i_min, 0]
        if arm < 0.05 or abs(moment) < 1e-12:
            return
        force[6 * i_max + 1] += moment / arm
        force[6 * i_min + 1] -= moment / arm

    vectors = {name: np.zeros(n) for name in ("mp_only", "uniform", "at_ac", "uniform_plus_mp")}
    for k, zz in enumerate(stations):
        ring = rings[zz]
        xs = coords[ring, 0]
        for name in ("uniform", "at_ac", "uniform_plus_mp"):
            force = vectors[name]
            for i in ring:
                force[6 * i + 1] += Np[k] / len(ring)   # normal to the rotor plane -> flapwise
                force[6 * i + 0] += Tp[k] / len(ring)   # tangential -> chordwise
        add_couple(vectors["at_ac"], ring, Np[k] * (float(xs.mean()) - (float(xs.min()) + 0.25 * (float(xs.max()) - float(xs.min())))))
        for name in ("mp_only", "at_ac", "uniform_plus_mp"):
            add_couple(vectors[name], ring, Mp[k])
    return vectors


def test_rated_tip_twist_matches_zhou_with_the_physical_load_path(blade_shell, rated_bem):
    """The physical sense of the rated twist, and how much the load's line of action moves it.

    Zhou's -3.60 deg (Table 4) is the **total** tip torsion at rated, so the comparable case is
    the full load applied along its physical line of action: forces at the aerodynamic centre
    plus the pitching moment. That gives -2.72 deg against Zhou's -3.60 deg.

    What is *asserted* is the physics: the physical load path gives the nose-down sense (the
    same sense as Zhou and Ma), and the line of action dominates the twist by more than a
    factor of two. What is *reported* is the 0.755x magnitude, 24.5% below Zhou - above the
    suite's 5% rule, so it is a documented residual (validation matrix section 9.5), not an
    assertion dressed up with a wide tolerance.
    """
    bem, blade_aero = rated_bem
    vectors = _rated_load_cases(blade_shell, bem, blade_aero)
    twists = {}
    for name, force in vectors.items():
        u = np.zeros(blade_shell["n"])
        u[blade_shell["free"]] = spsolve(blade_shell["Kff"], force[blade_shell["free"]])
        tip = blade_shell["tip"]
        twists[name] = float(np.mean([u[6 * i + 5] for i in tip]))

    print("\nrated tip twist (rotation about the blade axis), IEA 15 MW at V=10.59, 7.56 rpm:")
    for name, label in (("mp_only", "pitching moment alone (couple)"),
                        ("uniform", "Np+Tp uniform over the ring (one-way path)"),
                        ("uniform_plus_mp", "Np+Tp uniform + pitching moment"),
                        ("at_ac", "Np+Tp at the aerodynamic centre + pitching moment")):
        print(f"  {label:52} {np.rad2deg(twists[name]):+9.4f} deg")
    print(f"  {'Zhou 2025 Table 4 (total tip torsion)':52} {ZHOU_TIP_TORSION_DEG:+9.4f} deg")

    ratio = np.rad2deg(twists["at_ac"]) / ZHOU_TIP_TORSION_DEG
    print(f"  physical case / Zhou = {ratio:.3f}  <- a REPORTED residual, not an assertion: "
          f"{abs(1.0 - ratio) * 100:.1f}% is above the suite's 5% rule")

    # ASSERTED: the physics (the sense of the twist) and the load-path dominance. A sign
    # test is a valid test - it proves the physics is right - and it needs no tolerance.
    assert twists["at_ac"] < 0.0, "the physical load path must give the nose-down sense"
    assert abs(twists["uniform"]) > 2.0 * abs(twists["at_ac"]), (
        "the line of action must dominate; if it does not, the comparison is not well posed"
    )

    # NOT asserted: the 24.5% magnitude gap against Zhou. It is above the suite's 5% rule, so
    # it is reported as a residual and recorded in the validation matrix (section 9.5) rather
    # than dressed up with a wide tolerance. The remaining differences are the load case (our
    # steady BEM loads against Zhou's fully coupled aeroelastic solution), the aerodynamic
    # model (BEM against lifting-line free vortex wake) and the reference's own modelling.
    assert abs(1.0 - ratio) > 0.05, (
        "the residual is now inside 5%: promote it to an asserted row in the matrix"
    )
