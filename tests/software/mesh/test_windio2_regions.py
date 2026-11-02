"""Public placement contracts from the pinned windIO midpoint-origin arcs."""

import numpy as np
import pytest
import yaml
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.blade import Blade
from tests.software.mesh.test_windio2_official import OFFICIAL, ROOT, source_reference_length


@pytest.fixture(scope="module")
def official():
    source = yaml.safe_load(OFFICIAL.read_text())
    blade = Blade().read_yaml(str(OFFICIAL))
    return blade, source["components"]["blade"]


def source_arc(blade, source, anchor, handle):
    anchors = {a["name"]: a for a in source["structure"]["anchors"]}
    curve = anchors[anchor][handle]
    while "anchor" in curve:
        ref = curve["anchor"]
        curve = anchors[ref["name"]][ref["handle"]]
    eta = blade.geometry.ispan / source_reference_length(source["reference_axis"])
    g = blade.geometry
    return g.arclength[0] + (1 - np.interp(eta, curve["grid"], curve["values"])) * (
        g.arclength[-1] - g.arclength[0]
    )


def point_at(blade, arcs):
    g = blade.geometry
    result = []
    for k, a in enumerate(arcs):
        # Closed source TEs repeat the midpoint exactly. Removing repeated
        # interpolation abscissae does not move the source contour.
        arc, indices = np.unique(g.arclength[:, k], return_index=True)
        result.append(PchipInterpolator(arc, g.coordinates[indices, :, k])(a))
    return np.array(result).T


def test_web_midpoint_origin_and_physical_endpoints(official):
    blade, source = official
    # Source identities replace the obsolete fixed-12-row representation. The
    # independent oracle still uses source curves and the completed public loft.
    for name in ("web1", "web0"):
        for side, handle in enumerate(("start_nd_arc", "end_nd_arc")):
            expected = source_arc(blade, source, name, handle)
            xyz = point_at(blade, expected)
            for k in (4, 20, 40):
                eta = blade.geometry.ispan[k] / source_reference_length(source['reference_axis'])
                point = dict(blade.windio2_topology.section(eta).webs)[name][side]
                assert point.arc == pytest.approx(expected[k], abs=1e-12)
                np.testing.assert_allclose(point.xyz, xyz[:, k], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("anchor,handle", [
    ("LE_reinforcement", "end_nd_arc"),
    ("LE_reinforcement", "start_nd_arc"),
    ("TE_reinforcement_PS", "start_nd_arc"),
    ("TE_reinforcement_SS", "end_nd_arc"),
])
def test_reinforcement_boundary_arcs_and_coordinates(official, anchor, handle):
    blade, source = official
    expected = source_arc(blade, source, anchor, handle)
    xyz = point_at(blade, expected)
    for k in (4, 20, 40):
        eta = blade.geometry.ispan[k] / source_reference_length(source['reference_axis'])
        point, = [p for p in blade.windio2_topology.section(eta).points
                  if ('anchor', anchor, handle) in p.identities]
        assert point.arc == pytest.approx(expected[k], abs=1e-12)
        np.testing.assert_allclose(point.xyz, xyz[:, k], rtol=1e-12, atol=1e-12)


def test_v1_public_web_origin_remains_corner_based():
    blade = Blade().read_yaml(str(ROOT / "tests/IEA-15-240-RWT.yaml"))
    g = blade.geometry
    eta = (g.ispan - g.ispan[0]) / (g.ispan[-1] - g.ispan[0])
    for actual, (grid, start, end) in zip(blade.keypoints.web_attach_arcs,
                                         blade.definition.web_nd_arc, strict=True):
        expected = np.array([g.arclength[1] + (1 - np.interp(eta, grid, nd)) *
                             (g.arclength[-2] - g.arclength[1]) for nd in (end, start)])
        np.testing.assert_array_equal(actual, expected)
