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
    angle = np.deg2rad(twist) + correction
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
        np.testing.assert_allclose(measured, np.deg2rad(twist), atol=2e-13, rtol=0)


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


def test_v1_public_reader_retains_parallel_frame_and_original_twist():
    blade = Blade()
    blade.read_yaml(str(ROOT / "tests/IEA-15-240-RWT.yaml"))
    g = blade.geometry
    for k in range(len(g.ispan)):
        p = g.profiles[:, :, k]
        x = -(p[:, 0] - g.ichordoffset[k]) * g.ichord[k]
        y = p[:, 1] * g.ichord[k]
        angle = np.deg2rad(-g.idegreestwist[k])
        expected = np.column_stack((
            np.cos(angle) * x - np.sin(angle) * y - g.isweep[k],
            np.sin(angle) * x + np.cos(angle) * y + g.iprebend[k],
            np.full(len(x), g.ispan[k]),
        ))
        np.testing.assert_array_equal(g.coordinates[:, :, k], expected)
