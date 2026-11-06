"""Independent source-coordinate loft and curvilinear station checks; no FE solve."""

from pathlib import Path

import numpy as np
import yaml
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.blade import Blade

ROOT = Path(__file__).resolve().parents[3]
OFFICIAL = ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"


def read_copy(tmp_path, source):
    path = tmp_path / "surface.yaml"
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    blade = Blade()
    blade.read_yaml(str(path))
    assert path.read_bytes() == before
    return blade


def source_curve(airfoil, fractions):
    """Independent whole-curve source oracle, with no sidewise normalization."""
    xy = np.array([airfoil["coordinates"][c] for c in "xy"]).T
    le = np.argmin(xy[:, 0])
    if xy[:le, 1].mean() > xy[le:, 1].mean():
        xy = xy[::-1]  # pressure first for the published source convention
    distance = np.r_[0., np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    return PchipInterpolator(distance / distance[-1], xy)(fractions)


def test_official_surface_uses_rthick_and_open_curve_pchip():
    before = OFFICIAL.read_bytes()
    source = yaml.safe_load(before)
    outer = source["components"]["blade"]["outer_shape"]
    blade = Blade()
    blade.read_yaml(str(OFFICIAL))
    g = blade.geometry
    masters = sorted(source["airfoils"], key=lambda a: a["rthick"])
    eta = g.ispan / g.ispan[-1]
    thick = PchipInterpolator(outer["rthick"]["grid"], outer["rthick"]["values"])(eta)
    # The representation publishes sample parameters, not an independently
    # trusted shape. Evaluate EVERY master at each SAME whole-curve fraction.
    expected, by_position = [], []
    for k, target in enumerate(thick):
        fractions = g.surface_fraction[:, k]
        assert fractions[0] == 0 and fractions[-1] == 1
        assert np.all(np.diff(fractions) > 0)
        curves = [source_curve(a, fractions) for a in masters]
        expected.append(PchipInterpolator([a["rthick"] for a in masters], curves, axis=0)(target))
        by_position.append(PchipInterpolator(
            [a["spanwise_position"] for a in outer["airfoils"]],
            [source_curve(next(m for m in masters if m["name"] == a["name"]), fractions)
             for a in outer["airfoils"]], axis=0,
        )(eta[k]))
    expected = np.array(expected)
    np.testing.assert_allclose(g.profiles[1:-1].transpose(2, 0, 1), expected,
                               atol=2e-13, rtol=0)
    midpoint = (expected[:, 0] + expected[:, -1]) / 2
    np.testing.assert_allclose(g.profiles[0].T, midpoint, atol=2e-13, rtol=0)
    np.testing.assert_array_equal(g.profiles[0], g.profiles[-1])
    np.testing.assert_array_equal(g.profiles[g.LEindex, 0], g.profiles[1:-1, 0].min(axis=0))
    # LE is a landmark on the completed loft, NOT an imposed common fraction .5.
    assert np.max(np.abs(g.surface_fraction[g.LEindex - 1] - .5)) > .001
    assert np.max(np.abs(g.profiles[1, 1] - g.profiles[-2, 1])) > .01
    # Thick-section x positions must follow arc fractions, not a common chord grid.
    assert np.max(np.abs(expected[4, :g.LEindex, 0] - expected[20, :g.LEindex, 0])) > .01
    # rthick dominates the simultaneously supplied, non-equivalent position grid.
    assert np.max(np.abs(np.array(by_position) - expected)) > .001
    assert OFFICIAL.read_bytes() == before


def test_oblique_stations_are_curvilinear_and_keep_native_field_grids(tmp_path):
    source = yaml.safe_load(OFFICIAL.read_text())
    model = source["components"]["blade"]
    origin, slope = np.array([.3, -.4, 1.2]), np.array([12., -5., 117.])
    for i, c in enumerate("xyz"):
        grid = np.asarray(model["reference_axis"][c]["grid"])
        model["reference_axis"][c]["values"] = (origin[i] + slope[i] * grid).tolist()
    outer = model["outer_shape"]
    outer["twist"] = {"grid": [0., .3, .7, 1.], "values": [20., 11., 3., -2.]}
    outer["section_offset_y"] = {"grid": [0., .4, 1.], "values": [2., 1., .2]}
    blade = read_copy(tmp_path, source)
    g, d = blade.geometry, blade.definition
    length = np.linalg.norm(slope)
    np.testing.assert_allclose(d.span[-1], length, atol=2e-13, rtol=0)
    eta = g.ispan / length
    np.testing.assert_allclose(d.chord_curve[0], np.array(outer["chord"]["grid"]) * length,
                               atol=2e-13, rtol=0)
    for field, actual in (("twist", g.idegreestwist), ("chord", g.ichord),
                          ("section_offset_y", g.ichordoffset * g.ichord)):
        curve = outer[field]
        np.testing.assert_allclose(actual, PchipInterpolator(curve["grid"], curve["values"])(eta),
                                   atol=2e-13, rtol=0)
    # Recover translations independently from the public profile/coordinates.
    tangent = slope / length
    first = np.array([tangent[2], 0., -tangent[0]])
    first /= np.linalg.norm(first)
    second = np.cross(tangent, first)
    correction = np.arctan2(second[0], first[0])
    for k in (0, 20, -1):
        angle = np.deg2rad(g.idegreestwist[k]) + correction
        basis = np.column_stack((np.cos(angle)*first + np.sin(angle)*second,
                                 -np.sin(angle)*first + np.cos(angle)*second))
        local = np.column_stack((g.profiles[:, 1, k] * g.ichord[k],
                                 (g.profiles[:, 0, k] - g.ichordoffset[k]) * g.ichord[k]))
        windio = local @ basis.T + origin + eta[k] * slope
        np.testing.assert_allclose(g.coordinates[:, :, k], windio[:, [1, 0, 2]] * [-1, 1, 1],
                                   atol=2e-13, rtol=0)


def test_v1_public_surface_preserves_chordwise_camber_thickness_recipe():
    blade = Blade()
    blade.read_yaml(str(ROOT / "tests/IEA-15-240-RWT.yaml"))
    g = blade.geometry
    for k in range(len(g.ispan)):
        thickness = g.ithickness[:, k]
        half = .5 * (thickness * (g.ipercentthick[k] / (thickness.max() * 100)))
        expected = np.column_stack((
            np.r_[g.ic[-1, k], g.ic[::-1, k], g.ic[1:, k], g.ic[-1, k]],
            np.r_[0., (g.icamber[:, k] - half)[::-1], (g.icamber[:, k] + half)[1:], 0.],
        ))
        np.testing.assert_array_equal(g.profiles[:, :, k], expected)
