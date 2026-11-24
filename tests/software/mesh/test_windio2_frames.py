"""Independent section-frame checks through the public YAML reader (no FE solve)."""

from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.blade import Blade

ROOT = Path(__file__).resolve().parents[3]
OFFICIAL = ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
# windIO (suction, trailing-edge, span) -> NuMAD (-trailing-edge, suction, span).
TO_NUMAD = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])


def expected_frame(tangent, twist):
    """Construct basis columns geometrically, independently of production rotations.

    My @ Mx sends ez to t. Its first column lies in the x-z plane;
    its second column is t cross first. This fixes the forward-axis branch.
    """
    t = tangent / np.linalg.norm(tangent)
    first = np.array([t[2], 0., -t[0]])
    first /= np.linalg.norm(first)
    second = np.cross(t, first)
    correction = np.arctan2(second[0], first[0])
    # windIO twist is positive to feather (LE upwind): a NEGATIVE rotation about the section normal in the
    # windIO frame (x suction/downwind, y TE).  The correction defines the zero-twist reference and keeps its sign.
    angle = correction - np.deg2rad(twist)
    chord_normal = np.cos(angle) * first + np.sin(angle) * second
    chord = -np.sin(angle) * first + np.cos(angle) * second
    return np.column_stack((chord_normal, chord, t)), correction


def load_copy(tmp_path, source):
    path = tmp_path / "blade.yaml"
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    blade = Blade()
    blade.read_yaml(str(path))
    assert path.read_bytes() == before
    return blade


@pytest.fixture(scope="module", params=[0., 23.])
def oblique(request, tmp_path_factory):
    source = yaml.safe_load(OFFICIAL.read_text())
    model = source["components"]["blade"]
    # Constant analytic slope with both transverse components: correction != 0.
    origin = np.array([0.3, -0.4, 1.2])
    slope = np.array([2.1, -1.7, 117.])
    for i, name in enumerate("xyz"):
        grid = np.asarray(model["reference_axis"][name]["grid"])
        model["reference_axis"][name]["values"] = (origin[i] + slope[i] * grid).tolist()
    model["outer_shape"]["twist"]["values"] = [request.param] * len(
        model["outer_shape"]["twist"]["grid"]
    )
    blade = load_copy(tmp_path_factory.mktemp("oblique"), source)
    return blade, origin, slope, request.param


@pytest.fixture(scope="module", params=[10., -6.])
def oblique_free_straight(request, tmp_path_factory):
    """Straight axis along +z through the origin, constant twist: correction == 0."""
    source = yaml.safe_load(OFFICIAL.read_text())
    model = source["components"]["blade"]
    for name in "xyz":
        grid = np.asarray(model["reference_axis"][name]["grid"])
        model["reference_axis"][name]["values"] = (117. * grid if name == "z" else 0. * grid).tolist()
    model["outer_shape"]["twist"]["values"] = [request.param] * len(model["outer_shape"]["twist"]["grid"])
    return load_copy(tmp_path_factory.mktemp("straight"), source), request.param


def local_points(geometry, k):
    # Use the published generated section, not an independently re-lofted airfoil:
    # these tests isolate frames, not the still-pending surface interpolation.
    profile = geometry.profiles[:, :, k]
    return np.column_stack((
        profile[:, 1] * geometry.ichord[k],
        (profile[:, 0] - geometry.ichordoffset[k]) * geometry.ichord[k],
        np.zeros(len(profile)),
    ))


def fitted_frame(geometry, k):
    local = local_points(geometry, k)[:, :2]
    coordinates = geometry.coordinates[:, :, k]
    # Center before fitting: a ~118 m translation must not dominate the small
    # tip section in the least-squares solve for its dimensionless frame.
    local_mean = local.mean(axis=0)
    coordinate_mean = coordinates.mean(axis=0)
    fit = np.linalg.lstsq(local - local_mean, coordinates - coordinate_mean, rcond=None)[0]
    return fit.T, coordinate_mean - local_mean @ fit


def test_oblique_section_plane_follows_analytic_tangent(oblique):
    blade, _, slope, _ = oblique
    expected = TO_NUMAD @ (slope / np.linalg.norm(slope))
    for k in (0, 20, -1):
        basis, _ = fitted_frame(blade.geometry, k)
        normal = np.cross(basis[:, 0], basis[:, 1])
        np.testing.assert_allclose(normal, expected, atol=2e-13, rtol=0)


def test_oblique_twist_correction_and_rotation_order(oblique):
    blade, _, slope, twist = oblique
    frame, correction = expected_frame(slope, twist)
    assert abs(correction) > 1e-5  # distinguishes uncorrected My @ Mx @ Mz
    for k in (0, 20, -1):
        basis, _ = fitted_frame(blade.geometry, k)
        np.testing.assert_allclose(basis, (TO_NUMAD @ frame)[:, :2], atol=2e-13, rtol=0)
        # Correction cancels the first-row component of the rotated y column at zero twist.
        untilted = TO_NUMAD.T @ basis
        measured = np.arctan2(
            np.dot(untilted[:, 0], expected_frame(slope, 0)[0][:, 1]),
            np.dot(untilted[:, 0], expected_frame(slope, 0)[0][:, 0]),
        )
        # was +deg2rad(twist): the mirrored (wrong) pretwist sense, see test_chord_angle_is_minus_twist_*
        np.testing.assert_allclose(measured, -np.deg2rad(twist), atol=2e-13, rtol=0)


def test_oblique_coordinates_and_root_translation(oblique):
    blade, origin, slope, twist = oblique
    geometry = blade.geometry
    frame, _ = expected_frame(slope, twist)
    eta = geometry.ispan / blade.definition.span[-1]
    for k in (0, 20, -1):
        reference = TO_NUMAD @ (origin + eta[k] * slope)
        _, translation = fitted_frame(geometry, k)
        np.testing.assert_allclose(translation, reference, atol=2e-13, rtol=0)
        expected = local_points(geometry, k) @ frame.T @ TO_NUMAD.T + reference
        np.testing.assert_allclose(geometry.coordinates[:, :, k], expected, atol=2e-13, rtol=0)


def test_official_curved_axis_uses_native_pchip_derivative():
    source = yaml.safe_load(OFFICIAL.read_text())["components"]["blade"]
    blade = Blade()
    blade.read_yaml(str(OFFICIAL))
    geometry = blade.geometry
    eta = geometry.ispan / blade.definition.span[-1]
    curves = [PchipInterpolator(source["reference_axis"][c]["grid"],
                                source["reference_axis"][c]["values"]) for c in "xyz"]
    # Include the inserted chord knot .995 as well as native knots and the root.
    for k in (0, 30, np.flatnonzero(np.isclose(eta, .995))[0], len(eta) - 1):
        tangent = np.array([c.derivative()(eta[k]) for c in curves])
        basis, translation = fitted_frame(geometry, k)
        np.testing.assert_allclose(np.cross(basis[:, 0], basis[:, 1]),
                                   TO_NUMAD @ (tangent / np.linalg.norm(tangent)),
                                   atol=2e-13, rtol=0)
        np.testing.assert_allclose(translation, TO_NUMAD @ np.array([c(eta[k]) for c in curves]),
                                   atol=2e-13, rtol=0)


def chord_angle(geometry, k, tangent):
    """Chord angle psi (TE -> LE, right-handed about the axis) in the plane normal to the axis, in the mesh frame.

    Reference direction, built from the physical definition and not from the production rotations: the zero-twist
    LE direction is the unit vector normal to the axis with no windIO-x (suction) component, pointing to -y (LE);
    for a straight axis (t = +z) that is +x in the mesh frame, so psi = atan2(LE_y, LE_x).
    """
    t = tangent / np.linalg.norm(tangent)
    te_ref = np.array([0., t[2], -t[1]])
    le_ref = TO_NUMAD @ (-te_ref / np.linalg.norm(te_ref))
    t_mesh = TO_NUMAD @ t
    basis, _ = fitted_frame(geometry, k)
    le = -basis[:, 1]  # local y is the TE direction (profile x, TE at 1)
    return np.arctan2(le @ np.cross(t_mesh, le_ref), le @ le_ref)


def test_chord_angle_is_minus_twist_straight_axis(oblique_free_straight):
    """windIO/AeroDyn twist is positive to feather (LE toward upwind = -y in the mesh): psi = -twist, 1e-9."""
    blade, twist = oblique_free_straight
    for k in (0, 20, -1):
        psi = chord_angle(blade.geometry, k, np.array([0., 0., 1.]))
        np.testing.assert_allclose(psi, -np.deg2rad(twist), atol=1e-9, rtol=0)


def test_chord_angle_is_minus_twist_oblique_axis(oblique):
    """Same on a constant oblique axis (correction != 0): psi measured in the plane normal to the axis, 1e-9."""
    blade, _, slope, twist = oblique
    for k in (0, 20, -1):
        psi = chord_angle(blade.geometry, k, slope)
        np.testing.assert_allclose(psi, -np.deg2rad(twist), atol=1e-9, rtol=0)


def test_chord_angle_is_minus_twist_official_prebent_swept_axis():
    """Official blade (prebend + sweep, varying twist and tangent): psi = -twist(eta) at stations along the span.

    The identity is exact (rotation about the axis relative to the zero-twist reference); the tolerance 1e-9 rad
    only absorbs the least-squares frame fit."""
    source = yaml.safe_load(OFFICIAL.read_text())["components"]["blade"]
    blade = Blade()
    blade.read_yaml(str(OFFICIAL))
    g = blade.geometry
    eta = g.ispan / blade.definition.span[-1]
    curves = [PchipInterpolator(source["reference_axis"][c]["grid"],
                                source["reference_axis"][c]["values"]) for c in "xyz"]
    assert np.ptp(g.idegreestwist) > 10.  # a real twist range is exercised
    for k in (3, 15, 30, len(eta) // 2, len(eta) - 2):
        tangent = np.array([c.derivative()(eta[k]) for c in curves])
        assert abs(tangent[0]) + abs(tangent[1]) > 1e-4  # curved/prebent: the axis is not +z
        psi = chord_angle(g, k, tangent)
        np.testing.assert_allclose(psi, -np.deg2rad(g.idegreestwist[k]), atol=1e-9, rtol=0)


def test_positive_twist_moves_leading_edge_upwind(oblique_free_straight):
    """Mesh frame: +y is suction/downwind, LE at +x.  Twist t > 0 moves the LE toward -y (upwind) and the TE to +y."""
    blade, twist = oblique_free_straight
    g = blade.geometry
    k = 20
    local = local_points(g, k)
    # untwisted straight-axis section in the mesh frame: (x, y) = (-local_y, local_x), i.e. TO_NUMAD on (x_w, y_w)
    untwisted = local @ TO_NUMAD.T
    le, te = g.LEindex, 0
    d_le_y = g.coordinates[le, 1, k] - (untwisted[le, 1] + g.coordinates[:, 1, k].mean() - untwisted[:, 1].mean())
    d_te_y = g.coordinates[te, 1, k] - (untwisted[te, 1] + g.coordinates[:, 1, k].mean() - untwisted[:, 1].mean())
    assert untwisted[le, 0] > 0 > untwisted[te, 0]  # LE at +x, TE at -x
    assert (twist > 0) == (d_le_y < 0) and (twist > 0) == (d_te_y > 0)


def test_v2_blade_cache_key_follows_the_section_frame_source(tmp_path, monkeypatch):
    """The on-disk v2 mesh cache must not survive a change of geometry.py (which builds the section frames and twist):
    before the pretwist fix its digest did not cover that file and stale mirrored meshes would be reused silently."""
    import sys

    sys.path.insert(0, str(ROOT / "tests" / "support"))
    import v2_blade

    assert v2_blade.GEOMETRY.name == "geometry.py" and v2_blade.GEOMETRY.exists()
    for mesher in ("gmsh",):  # the windio2 mesher is retired (gmsh-official-mesher T12)
        before = v2_blade.cache_key(0.25, mesher=mesher)
        edited = tmp_path / f"geometry_{mesher}.py"
        edited.write_bytes(v2_blade.GEOMETRY.read_bytes() + b"\n# edit\n")
        monkeypatch.setattr(v2_blade, "GEOMETRY", edited)
        assert v2_blade.cache_key(0.25, mesher=mesher) != before
        monkeypatch.undo()
        assert v2_blade.cache_key(0.25, mesher=mesher) == before
    assert v2_blade.CACHE_VERSION >= 4  # the pretwist-fix bump: pre-fix cache files (v <= 3) are unreachable


def test_apply_ad_loads_builds_the_physical_airfoil_frame():
    """S-5 load mapping: Fn/Ft act in the frame whose chord sits at psi = -BlTwist (positive twist = LE upwind), the same
    frame as both meshes (numad v1 and windIO-v2/gmsh).  A pure Fn at twist th points along (+sin th, cos th) in the mesh
    x-y plane (not (-sin th, cos th), the frame of the former +twist mesh), a pure Ft along (cos th, -sin th)."""
    import sys
    from types import SimpleNamespace

    sys.path.insert(0, str(ROOT / "tools"))
    import run_s5_oneway as s5

    z = np.linspace(0.0, 100.0, 11)
    mesh = SimpleNamespace(node_count=len(z), nodes=[SimpleNamespace(coords=np.array([0., 0., v])) for v in z])
    asm = SimpleNamespace(dofs_count=6 * len(z), dofs_per_node=6)
    th = np.radians(12.0)
    r = np.array([0.0, 100.0])
    twist = np.full(2, 12.0)
    for fn, ft, expected in ((1.0, 0.0, [np.sin(th), np.cos(th)]), (0.0, 1.0, [np.cos(th), -np.sin(th)])):
        f = s5.apply_ad_loads(asm, mesh, r, np.full(2, fn), np.full(2, ft), twist).reshape(-1, 6)[:, :3]
        total = f.sum(axis=0)
        assert total[2] == pytest.approx(0.0, abs=1e-12)
        np.testing.assert_allclose(total[:2] / np.linalg.norm(total[:2]), expected, atol=1e-12, rtol=0)
        assert np.linalg.norm(total[:2]) == pytest.approx(100.0, rel=1e-9)  # fn*L conserved
