"""Source-driven topology through the public reader, without legacy region rows."""
import copy
import hashlib

import numpy as np
import pytest
import yaml
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.blade import Blade
from tests.software.mesh.test_windio2_official import OFFICIAL, ROOT


@pytest.fixture(scope="module")
def official():
    before = OFFICIAL.read_bytes()
    blade = Blade().read_yaml(str(OFFICIAL))
    assert OFFICIAL.read_bytes() == before
    return blade, yaml.safe_load(before), hashlib.sha256(before).hexdigest()


def test_coincident_le_has_identities_but_no_finite_cell(official):
    blade, source, digest = official
    topology = blade.windio2_topology
    eta = 0.08163265306122448
    section = topology.section(eta)
    ids = {("anchor", "LE_reinforcement", h) for h in ("start_nd_arc", "end_nd_arc")}
    points = [p for p in section.points if ids.intersection(p.identities)]
    assert len(points) == 1
    assert ids.issubset(points[0].identities)
    anchor = next(a for a in source["components"]["blade"]["structure"]["anchors"]
                  if a["name"] == "LE_reinforcement")
    curve = anchor["start_nd_arc"]
    assert points[0].nd_arc == np.interp(eta, curve["grid"], curve["values"])
    assert points[0].arc > 0  # source coincidence is NOT the geometric LE
    assert len(section.intervals) == len(section.points) - 1
    assert all(a.nd_arc < b.nd_arc for a, b in section.intervals)
    assert topology.source_sha256 == digest
    assert tuple(a["name"] for a in topology.source["components"]["blade"]["structure"]["anchors"]) == tuple(a["name"] for a in source["components"]["blade"]["structure"]["anchors"])
    with pytest.raises(TypeError):
        topology.source["components"]["blade"]["structure"]["anchors"][0]["name"] = "changed"


def test_web_positions_use_completed_full_midpoint_curve(official):
    blade, source, _ = official
    topology, g = blade.windio2_topology, blade.geometry
    structure = source["components"]["blade"]["structure"]
    anchors = {a["name"]: a for a in structure["anchors"]}
    k = 20
    eta = g.ispan[k] / blade.definition.span[-1]
    before = g.coordinates.copy()
    section = topology.section(eta)
    for web in structure["webs"]:
        actual = dict(section.webs)[web["name"]]
        for handle, point in zip(("start_nd_arc", "end_nd_arc"), actual, strict=True):
            curve = web[handle]
            while "anchor" in curve:
                ref = curve["anchor"]
                curve = anchors[ref["name"]][ref["handle"]]
            nd = np.interp(eta, curve["grid"], curve["values"])
            arc = g.arclength[0, k] + (1 - nd) * np.ptp(g.arclength[:, k])
            xyz = PchipInterpolator(g.arclength[:, k], g.coordinates[:, :, k])(arc)
            assert point.nd_arc == nd
            np.testing.assert_allclose(point.xyz, xyz, atol=2e-12, rtol=0)
            assert abs(arc - (g.arclength[1, k] + (1-nd) *
                              (g.arclength[-2, k]-g.arclength[1, k]))) > 1e-6
    np.testing.assert_allclose(section.points[0].xyz, g.coordinates[-1, :, k], atol=2e-12)
    np.testing.assert_allclose(section.points[-1].xyz, g.coordinates[0, :, k], atol=2e-12)
    np.testing.assert_array_equal(blade.geometry.coordinates, before)


def test_crossing_events_native_grids_and_partial_domains(tmp_path, official):
    _, original, _ = official
    source = copy.deepcopy(original)
    anchors = source["components"]["blade"]["structure"]["anchors"]
    anchors.extend([
        {"name": "cross_a", "start_nd_arc": {"grid": [.2, .8], "values": [.2, .8]}},
        {"name": "cross_b", "start_nd_arc": {"grid": [.2, .4, .8], "values": [.5, .5, .5]}},
        {"name": "alias", "start_nd_arc": {"anchor": {"name": "cross_a", "handle": "start_nd_arc"}}},
    ])
    path = tmp_path / "cross.yaml"
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    blade = Blade().read_yaml(str(path))
    topology = blade.windio2_topology
    assert path.read_bytes() == before
    native = next(c for c in topology.curves if c.identity == ("anchor", "cross_b", "start_nd_arc"))
    assert native.grid == (.2, .4, .8) and native.values == (.5, .5, .5)
    assert np.min(np.abs(np.array(topology.event_stations) - .5)) < 1e-14
    for eta in (.2, .4, .8):
        assert eta in topology.event_stations
    a, b, alias = [("anchor", name, "start_nd_arc") for name in ("cross_a", "cross_b", "alias")]
    middle = topology.section(.5)
    point = next(p for p in middle.points if a in p.identities)
    assert {a, b, alias}.issubset(point.identities)
    for eta in (.3, .7):  # arbitrary eta, not a resampled geometry station
        section = topology.section(eta)
        pa = next(p for p in section.points if a in p.identities)
        pb = next(p for p in section.points if b in p.identities)
        assert (pa.nd_arc < pb.nd_arc) == (eta < .5)
        assert pa.nd_arc == pytest.approx(eta, abs=1e-15)
    assert not any(a in p.identities for p in topology.section(.1).points)
    with pytest.raises(ValueError, match="eta must be finite and in \\[0, 1\\]"):
        topology.section(float("nan"))


def event_topology(official, levels):
    """Native eta line and constant curves; crossing roots are the levels."""
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
    blade, original, _ = official
    source = copy.deepcopy(original)
    anchors = [{'name': 'diagonal', 'start_nd_arc':
                {'grid': [0., 1.], 'values': [0., 1.]}}]
    for name, grid, level in levels:
        anchors.append({'name': name, 'start_nd_arc':
                        {'grid': grid, 'values': [level]*len(grid)}})
    source['components']['blade']['structure'].update(anchors=anchors, webs=[], layers=[])
    before = yaml.safe_dump(source)
    topology = WindIO2Topology(source, hashlib.sha256(before.encode()).hexdigest(), blade.definition)
    return topology, source, before


def test_same_crossing_from_coincident_segments_has_one_event(official):
    topology, source, before = event_topology(official, [
        ('full', [0., 1.], .3), ('partial', [.1, .8], .3)])
    # The independent oracle is y=eta intersecting the unchanged constant .3;
    # neither the constant's domain nor a different pair may move its root.
    assert [e for e in topology.event_stations if .25 < e < .35] == [.3]
    identities = {('anchor', name, 'start_nd_arc') for name in ('diagonal', 'full', 'partial')}
    point = next(p for p in topology.section(.3).points if ('anchor', 'diagonal', 'start_nd_arc') in p.identities)
    assert identities.issubset(point.identities)
    assert point.nd_arc == .3
    assert yaml.safe_dump(source) == before
    assert topology.source_sha256 == hashlib.sha256(before.encode()).hexdigest()
    assert [a['name'] for a in topology.source['components']['blade']['structure']['anchors']] == ['diagonal', 'full', 'partial']


def test_truly_distinct_adjacent_crossings_are_not_canonicalized(official):
    upper = float(np.nextafter(.3, 1.))
    topology, source, before = event_topology(official, [
        ('lower', [0., 1.], .3), ('upper', [0., 1.], upper)])
    assert [e for e in topology.event_stations if .25 < e < .35] == [.3, upper]
    for eta, name, other in ((.3, 'lower', 'upper'), (upper, 'upper', 'lower')):
        point = next(p for p in topology.section(eta).points if ('anchor', 'diagonal', 'start_nd_arc') in p.identities)
        assert ('anchor', name, 'start_nd_arc') in point.identities
        assert ('anchor', other, 'start_nd_arc') not in point.identities
        assert point.nd_arc == eta
    assert yaml.safe_dump(source) == before


def test_distinct_exact_crossings_with_one_float_image_reject_without_mutation(tmp_path, official):
    """Post-implementation characterization of the public precision guard."""
    from fractions import Fraction
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology

    # Independent analytic roots: eta=.3 and .7*eta=.21. No production
    # intersection/evaluation helper is used by this oracle. Ratios describe
    # parsed binary64 inputs, not exact decimal or physical quantities.
    diagonal_root = Fraction(*(.3).as_integer_ratio())
    scaled_root = Fraction(*(.21).as_integer_ratio()) / Fraction(*(.7).as_integer_ratio())
    assert diagonal_root != scaled_root
    image = float(diagonal_root)
    assert float(scaled_root) == image

    _, source, _ = event_topology(official, [
        ('level', [0., 1.], .3), ('scaled_level', [0., 1.], .21)])
    source['components']['blade']['structure']['anchors'].append(
        {'name': 'scaled', 'start_nd_arc': {'grid': [0., 1.], 'values': [0., .7]}})
    path = tmp_path/'collision.yaml'
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    digest = hashlib.sha256(before).hexdigest()
    topology = WindIO2Topology(source, digest, official[0].definition)
    curves, events = topology.curves, topology.event_stations
    provenance = [(c.identity, c.grid, c.values) for c in curves]
    ids = {name: ('anchor', name, 'start_nd_arc')
           for name in ('diagonal', 'level', 'scaled', 'scaled_level')}
    expected = {'diagonal': (0., 1.), 'level': (.3, .3),
                'scaled': (0., .7), 'scaled_level': (.21, .21)}
    for name, values in expected.items():
        curve = next(c for c in curves if c.identity == ids[name])
        assert curve.grid == (0., 1.) and curve.values == values
    assert image in events
    neighbors = tuple(topology.section(eta) for eta in (.25, .35))
    for section, before_crossing in zip(neighbors, (True, False), strict=True):
        points = {name: next(p for p in section.points if identity in p.identities)
                  for name, identity in ids.items()}
        assert len({p.nd_arc for p in points.values()}) == 4
        assert (points['diagonal'].nd_arc < points['level'].nd_arc) == before_crossing
        assert (points['scaled'].nd_arc < points['scaled_level'].nd_arc) == before_crossing

    with pytest.raises(ValueError) as caught:
        topology.section(image)
    assert str(caught.value) == (
        'Unsupported float64 topology event precision: '
        'distinct exact events have the same float64 image; source events were not merged')
    # There are no public allocation/section counters in this API. Check all
    # relevant public state and successful neighboring queries instead.
    assert topology.curves == curves and topology.event_stations == events
    assert [(c.identity, c.grid, c.values) for c in topology.curves] == provenance
    assert tuple(topology.section(eta) for eta in (.25, .35)) == neighbors
    assert topology.source_sha256 == digest
    assert path.read_bytes() == before and yaml.safe_dump(source).encode() == before
    assert [a['name'] for a in topology.source['components']['blade']['structure']['anchors']] == [
        'diagonal', 'level', 'scaled_level', 'scaled']


def test_v1_has_no_topology():
    blade = Blade().read_yaml(str(ROOT / "tests/IEA-15-240-RWT.yaml"))
    assert blade.windio2_topology is None


def test_ulp_level_source_events_are_one_event_in_extended_view(official):
    """TE_reinforcement_SS end and Shell_skin_inner start are the same physical
    event computed by two double-precision paths (3e-18 apart at this eta).
    The extended view identifies them; legitimately separate events stay apart.
    """
    blade, _, _ = official
    topology = blade.windio2_topology.with_extended_precision()
    ids = {("anchor", "Shell_skin_inner", "start_nd_arc"),
           ("anchor", "TE_reinforcement_SS", "end_nd_arc")}
    eta = 0.2383414361472924
    points = [p for p in topology.section(eta).points if ids.intersection(p.identities)]
    assert len(points) == 1 and ids <= set(points[0].identities)
    # 1e-4 apart in nd at the neighbouring event station: still two events.
    near = [p for p in topology.section(0.23636263040030311).points
            if ids.intersection(p.identities)]
    assert len(near) == 2


def test_event_identification_is_anchored_and_snaps_to_the_unit_interval_ends():
    from fractions import Fraction as F
    from aeroelast.models.blade.numad.objects.windio2_topology import (
        EVENT_IDENTIFICATION_TOL as TOL, _identify_ulp_events)
    step = TOL*F(3, 4)
    # A chain of values each closer than TOL to the previous one must not drift
    # past TOL from its first member: the third value starts a new event.
    got = _identify_ulp_events({'a': F(1, 3), 'b': F(1, 3)+step, 'c': F(1, 3)+2*step})
    assert got['a'] == got['b'] == F(1, 3) and got['c'] == F(1, 3)+2*step
    # Just beyond the tolerance stays distinct (exact, no rounding).
    far = _identify_ulp_events({'a': F(1, 3), 'b': F(1, 3)+TOL+F(1, 10**30)})
    assert far['a'] != far['b']
    # Values within TOL of the interval ends join the end, so they cannot form
    # a micro-band with the knots 0 and 1 that the section adds separately.
    ends = _identify_ulp_events({'lo': TOL/2, 'hi': 1-TOL/2, 'mid': F(1, 2)})
    assert ends == {'lo': 0, 'hi': 1, 'mid': F(1, 2)}
