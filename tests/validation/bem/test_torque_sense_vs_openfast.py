"""Torsional load-chain sense of the BEM -> shell path, against OpenFAST v5 (L3 T5).

Mesh frame of the official blade (0.25 m): leading edge at ``+x``, suction side and
thrust at ``+y`` (downwind), span ``+z``.  A rotation about ``+z`` carries the leading
edge downwind, i.e. toward stall: **+z is nose-UP**.  AeroDyn's ``Cm`` is positive
nose-up about ``x_a x y_a`` (suction normal x chord), which is ``+z`` here, so the
cambered section's polar moment is negative about ``+z`` and the AC lever of the
downwind normal force (AC ahead of the reference axis, at ``+x``) is positive.  The two
oppose and the physical rated section torque is nose-DOWN.

Commit 1146265 inverted three links (``_strip_moment_axis_sign``, the guards' nose-down
convention, ``_twist_mesh_to_bem``): the production torque came out ``+4.57e5 N.m`` and
a rotation toward stall RAISED theta.  An OpenFAST v5 AeroDyn+BeamDyn rated run
(``logs/check/l3-bd-openfast/REPORT.md``) twists the tip 3.40 deg toward feather and
cuts thrust 13.2 %; the inverted chain gave +5.28 deg toward stall.

Two assertions through the public ``BEMFSIParticipant`` / ``ForceProjector`` path:

(a) **Applied torque.** Metric: ``Mz = sum_i [x_i f_y,i - (y_i - y_ac,k) f_x,i]`` of the
    participant's projected nodal forces, per strip ``k`` about the ``z`` line through
    the station-ring AC (so the rotor-torque lever of ``Tp`` is not counted).
    Reference: ``B - A = sum_k (Mp_k + x_ac,k Np_k) dr_k`` from the rigid BEM result
    and the AC read off the mesh station ring (``-9.52e4 N.m`` here; the
    interpolated-ring script ``logs/check/l3-bd-openfast/torque_sense.py`` gives
    ``-9.41e4``, 1.1 % apart).  Tolerance 2 % of ``|B - A|`` plus the sign: the
    projection realises the requested section moment exactly (measured applied
    ``-9.5173e4`` against ``-9.5173e4``), so the 2 % covers only the two ways of reading
    the AC off the mesh ring (1.1 % measured); the sign error fixed here (``+4.58e5``)
    is 580 % of ``|B - A|``, so 2 % is a statement about the arithmetic and not the sign.

(b) **Feedback sense.** A rigid ``+z`` (nose-up) rotation must LOWER the section pitch
    angle ``theta`` (CCBlade ``alpha = phi - theta``), so alpha rises and thrust
    rises.  Metric: ``theta_def - theta_ref`` on every strip.  Tolerance: the rigid
    field is an exact section rotation, so ``-angle`` to 1e-3 relative (estimator
    invariance, not a physical bound).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.mesh.entities import MeshElement, Node  # noqa: E402
from aeroelast.models.blade.model import Blade  # noqa: E402
from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.validation.bem.test_force_projection_ac_datum import (  # noqa: E402
    AD_PRIMARY,
    ELEMENT_SIZE,
    MESH_RING_GAP_ELEMENT_FRACTION,
    NORMAL_DIR,
    SPAN_DIR,
    TANGENTIAL_DIR,
    YAML,
    _blunt_end,
    _geometry_load_frames,
    _station_ring_points,
)
from tests.validation.rotor.test_bem_fsi_deformed_geometry import (  # noqa: E402
    BEM_CONFIG,
    PRODUCTION_FRAME,
    PITCH_RATED,
    RPM_RATED,
    V_RATED,
    _rigid_section_rotation,
)

#: Independent rated torque about +z [N.m]: Sum Cm term (-2.7676e5) + AC lever
#: (+1.8263e5), ``logs/check/l3-bd-openfast/torque_sense.log``.  Stated for the
#: sign/size context; the asserted reference is recomputed below from the mesh.
INDEPENDENT_B_MINUS_A = -9.4133e4
TORQUE_REL_TOL = 0.02
RIGID_NOSE_UP_DEG = 2.0
THETA_REL_TOL = 1e-3


@pytest.fixture(scope="module")
def case():
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(YAML), element_size=ELEMENT_SIZE)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None
    coords = mesh.coords_array
    blade_aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    participant = BEMFSIParticipant(mesh, blade_aero, BEM_CONFIG)
    participant_prod = BEMFSIParticipant(mesh, blade_aero, {**BEM_CONFIG, **PRODUCTION_FRAME})
    return mesh, coords, blade_aero, participant, participant_prod


def test_rated_section_torque_is_nose_down(case):
    _, coords, blade_aero, participant, _ = case
    projector = participant._projector
    forces, bem = participant._compute_forces(np.zeros((len(coords), 3)))

    tol = MESH_RING_GAP_ELEMENT_FRACTION * ELEMENT_SIZE
    rings = [
        _station_ring_points(coords, strip, blade_aero, k, SPAN_DIR, tol)
        for k, strip in enumerate(projector._strips)
    ]
    frames = _geometry_load_frames(
        rings, [float(s.dr) for s in projector._strips], SPAN_DIR, NORMAL_DIR, TANGENTIAL_DIR
    )

    applied = 0.0
    cm_term = 0.0
    lever = 0.0
    for k, strip in enumerate(projector._strips):
        pts = rings[k]
        le_i, te_i, _, _ = _blunt_end(pts, frames[k][0], SPAN_DIR)
        ac = pts[le_i] + blade_aero.stations[k].airfoil.aerodynamic_center * (pts[te_i] - pts[le_i])
        idx = strip.node_indices
        f = forces[idx]
        x = coords[idx]
        applied += float(np.sum(x[:, 0] * f[:, 1] - (x[:, 1] - ac[1]) * f[:, 0]))
        dr = float(strip.dr)
        cm_term += float(bem.Mp[k]) * dr
        lever += float(ac[0]) * float(bem.Np[k]) * dr
    expected = cm_term + lever
    print(
        f"\nCm term {cm_term:+.4e}  AC lever {lever:+.4e}  B-A {expected:+.4e} "
        f"(logged {INDEPENDENT_B_MINUS_A:+.4e})  applied {applied:+.4e} N.m about +z"
    )
    assert expected < 0.0, "independent torque must be nose-down (negative about +z)"
    assert applied < 0.0, (
        f"applied rated torque {applied:+.4e} N.m about +z is nose-UP; the physical "
        f"torque is nose-down ({expected:+.4e})"
    )
    assert applied == pytest.approx(expected, rel=TORQUE_REL_TOL)


@pytest.mark.parametrize("which", ["default", "production"])
def test_nose_up_rotation_lowers_theta(case, which):
    _, coords, _, participant, participant_prod = case
    p = participant if which == "default" else participant_prod
    angle = np.deg2rad(RIGID_NOSE_UP_DEG)
    rigid = _rigid_section_rotation(coords, p._span_dir, angle)  # +z rotation: nose-up
    _, twist_def = p._compute_deformed_geometry(rigid)
    dtheta = twist_def - p._ref_twist
    print(
        f"\n{which}: factor {p._twist_mesh_to_bem:+.0f}, dtheta [deg] "
        f"{np.rad2deg(dtheta.min()):+.4f} .. {np.rad2deg(dtheta.max()):+.4f}"
    )
    assert np.all(dtheta < 0.0), "a nose-up (+z) rotation must lower theta (alpha rises)"
    np.testing.assert_allclose(dtheta, -angle, rtol=THETA_REL_TOL)
    base, _ = p._rebuild_bem_solver(p._ref_r, p._ref_twist)
    deformed, _ = p._rebuild_bem_solver(p._ref_r, twist_def)
    assert (
        deformed.compute(V_RATED, RPM_RATED, PITCH_RATED).thrust
        > base.compute(V_RATED, RPM_RATED, PITCH_RATED).thrust
    )
