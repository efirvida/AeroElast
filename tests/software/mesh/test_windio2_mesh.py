"""Isolated skin mesh: source boundaries, not the legacy twelve-region path."""
import copy
from collections import Counter
import hashlib

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.objects.blade import Blade
from tests.software.mesh.test_windio2_topology import official  # noqa: F401
from tests.software.mesh.test_windio2_laminate import expected


def mesh_function():
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_skin_mesh
    return generate_skin_mesh


def check_mesh(mesh, source):
    assert mesh.elements.shape[1] == 4
    assert np.isfinite(mesh.nodes).all()
    edges = Counter()
    for element, cell, uv in zip(mesh.elements, mesh.cells, mesh.element_parameters, strict=True):
        assert len(set(element)) == 4
        xyz = mesh.nodes[element]
        normal = np.cross(xyz[1]-xyz[0], xyz[3]-xyz[0])
        assert np.linalg.norm(normal) > 0
        assert np.dot(normal, np.cross(xyz[2]-xyz[1], xyz[0]-xyz[1])) > 0
        eta, q = np.mean(uv, axis=0)
        want = expected(source, eta, 'skin', q, q)
        assert [(l.name, l.material_name) for l in cell.layers] == [(r[0], r[1]) for r in want]
        np.testing.assert_allclose([l.thickness_m for l in cell.layers], [r[2] for r in want], atol=1e-16, rtol=0)
        for a, b in zip(element, np.roll(element, -1), strict=True):
            edges[tuple(sorted((a, b)))] += 1
    assert set(edges.values()) <= {1, 2}
    # All singly occupied edges must be on the requested root/tip cut, never
    # interior cracks at crossing, coincidence, dropoff or the TE seam.
    for (a, b), count in edges.items():
        if count == 1:
            assert mesh.parameters[a, 0] == mesh.parameters[b, 0]
            assert mesh.parameters[a, 0] in mesh.eta_bounds


def test_short_source_crossing_coincidence_dropoff(tmp_path, official):
    generate = mesh_function()
    _, original, _ = official
    source = copy.deepcopy(original)
    blade = source['components']['blade']
    for axis in ('x', 'y', 'z'):
        field = blade['reference_axis'][axis]
        field['values'] = [0.5*q if axis == 'z' else 0. for q in field['grid']]
    shape = blade['outer_shape']
    for key, value in [('chord', .2), ('section_offset_y', 0.), ('twist', 0.)]:
        shape[key]['values'] = [value]*len(shape[key]['grid'])
    path = tmp_path/'geometry.yaml'
    path.write_text(yaml.safe_dump(source))
    definition = Blade().read_yaml(str(path)).definition
    structure = blade['structure']
    base = copy.deepcopy(structure['layers'][0])
    base.update(start_nd_arc={'grid': [0., 1.], 'values': [0., 0.]},
                end_nd_arc={'grid': [0., 1.], 'values': [1., 1.]})
    structure.update(webs=[], layers=[base])
    structure['anchors'] = [
        {'name': 'moving', 'start_nd_arc': {'grid': [0., 1.], 'values': [.25, .75]}},
        {'name': 'fixed', 'start_nd_arc': {'grid': [0., 1.], 'values': [.5, .5]}},
    ]
    layer = copy.deepcopy(structure['layers'][0])
    layer.update(name='drop', start_nd_grid=.25, end_nd_grid=.75,
                 start_nd_arc={'anchor': {'name': 'moving', 'handle': 'start_nd_arc'}},
                 end_nd_arc={'anchor': {'name': 'fixed', 'handle': 'start_nd_arc'}})
    structure['layers'].append(layer)
    path = tmp_path/'short.yaml'
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
    topology = WindIO2Topology(source, hashlib.sha256(before).hexdigest(), definition)
    mesh = generate(topology, element_size=.25)
    check_mesh(mesh, source)
    assert {.25, .5, .75}.issubset(set(mesh.parameters[:, 0]))
    assert any(l.name == 'drop' for c in mesh.cells for l in c.layers)
    # Independent physical oracle: native source endpoints, zero twist/offset,
    # straight 0.5 m axis. No production geometry/position helper is used.
    from scipy.interpolate import PchipInterpolator
    masters = sorted(source['airfoils'], key=lambda m: m['rthick'])
    endpoints = [np.array([m['coordinates'][c] for c in 'xy'])[:, [0, -1]]
                 for m in masters]
    seam = mesh.parameters[:, 1] == 0
    etas = mesh.parameters[seam, 0]
    thick = PchipInterpolator(shape['rthick']['grid'], shape['rthick']['values'])(etas)
    xy = PchipInterpolator([m['rthick'] for m in masters], endpoints, axis=0)(thick).mean(axis=2)
    np.testing.assert_allclose(mesh.nodes[seam], np.column_stack((-.2*xy[:, 0], .2*xy[:, 1], .5*etas)),
                               atol=2e-13, rtol=0)
    assert path.read_bytes() == before


def test_official_patch_midpoint_geometry_and_ownership(official):
    generate = mesh_function()
    blade, source, digest = official
    g = blade.geometry
    eta = .08163265306122448
    mesh = generate(blade.windio2_topology, element_size=.25, eta_bounds=(eta, eta+.0001))
    check_mesh(mesh, source)
    # Independent existing geometry TE midpoint, not the production mesh helper.
    k = np.argmin(abs(g.ispan/blade.definition.span[-1]-eta))
    assert abs(g.ispan[k]/blade.definition.span[-1]-eta) < 1e-14
    seam = np.where((mesh.parameters[:, 0] == eta) & (mesh.parameters[:, 1] == 0))[0]
    assert len(seam) == 1
    np.testing.assert_allclose(mesh.nodes[seam[0]], (g.coordinates[1, :, k]+g.coordinates[-2, :, k])/2, atol=2e-12, rtol=0)
    from tests.software.mesh.test_windio2_official import OFFICIAL
    assert hashlib.sha256(OFFICIAL.read_bytes()).hexdigest() == digest


def synthetic_chart(tmp_path, official, case):
    """Real source topology on a short loft; no test-local quadrangulator."""
    _, original, _ = official
    source = copy.deepcopy(original)
    blade = source['components']['blade']
    for axis in ('x', 'y', 'z'):
        field = blade['reference_axis'][axis]
        field['values'] = [.5*q if axis == 'z' else 0. for q in field['grid']]
    for key, value in [('chord', .2), ('section_offset_y', 0.), ('twist', 0.)]:
        field = blade['outer_shape'][key]
        field['values'] = [value]*len(field['grid'])
    path = tmp_path/'synthetic.yaml'
    path.write_text(yaml.safe_dump(source))
    definition = Blade().read_yaml(str(path)).definition
    structure = blade['structure']
    base = copy.deepcopy(structure['layers'][0])
    base.update(start_nd_arc={'grid': [0., 1.], 'values': [0., 0.]},
                end_nd_arc={'grid': [0., 1.], 'values': [1., 1.]})
    layer = copy.deepcopy(base)
    layer['name'] = 'probe'
    left, right = [.3, .3], [.7, .7]
    if case == 'crossing':
        left, right = [.25, .75], [.5, .5]
    elif case == 'finite':
        right = [.3+2**-12]*2
    elif case == 'zero':
        right = left
    elif case == 'arc_precision':
        # Narrower than float64 resolves for chart sizes, but wider than the
        # ulp-level event identification (windio2_topology.EVENT_IDENTIFICATION_TOL):
        # a 1-ulp band is by design one event with its neighbour.
        right = [.3+2**-48]*2
    layer.update(start_nd_arc={'grid': [0., 1.], 'values': left},
                 end_nd_arc={'grid': [0., 1.], 'values': right})
    structure.update(webs=[], anchors=[], layers=[base, layer])
    if case == 'span_precision':
        layer['thickness'] = {'grid': [0., .5, float(np.nextafter(.5, 1.)), 1.],
                              'values': [.001]*4}
    path.write_text(yaml.safe_dump(source))
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
    return WindIO2Topology(source, hashlib.sha256(path.read_bytes()).hexdigest(), definition), source, path


@pytest.mark.parametrize('case', ['rectangle', 'crossing', 'finite', 'zero'])
def test_general_chart_coverage_and_source_interfaces(tmp_path, official, case):
    topology, source, path = synthetic_chart(tmp_path, official, case)
    before = path.read_bytes()
    mesh = mesh_function()(topology)
    check_mesh(mesh, source)  # includes full source ownership/thickness and edge incidence
    uv = mesh.element_parameters
    # Translation-centred shoelace avoids cancellation. Analytic chart area is
    # one for the complete skin, regardless of cuts, crossing or zero-width band.
    local = uv-uv[:, :1]
    twice_area = np.sum(local[:, :, 0]*np.roll(local[:, :, 1], -1, axis=1)
                        - local[:, :, 1]*np.roll(local[:, :, 0], -1, axis=1), axis=1)
    assert np.all(twice_area > 0)
    np.testing.assert_allclose(twice_area.sum()/2, 1., atol=32*np.finfo(float).eps, rtol=0)
    for r in (-1., 0., 1.):
        for s in (-1., 0., 1.):
            dr = np.array([-(1-s), 1-s, 1+s, -(1+s)])/4
            ds = np.array([-(1-r), -(1+r), 1+r, 1-r])/4
            du, dv = dr@local, ds@local
            assert np.all(du[:, 0]*dv[:, 1]-du[:, 1]*dv[:, 0] > 0)
    probe = np.array([any(l.name == 'probe' for l in c.layers) for c in mesh.cells])
    if case in ('finite', 'zero', 'rectangle'):
        width = {'finite': 2**-12, 'zero': 0., 'rectangle': .4}[case]
        np.testing.assert_allclose(twice_area[probe].sum()/2, width,
                                   atol=32*np.finfo(float).eps, rtol=0)
        for boundary in (.3, .3+width):
            # Constants are source floats, so bitwise equality is required,
            # not an epsilon that could erase the finite band.
            assert np.any(mesh.parameters[:, 1] == boundary)
    assert path.read_bytes() == before


@pytest.mark.parametrize('case', ['arc_precision', 'span_precision'])
def test_general_chart_preserves_previously_unrepresentable_subdivision(tmp_path, official, case):
    topology, source, path = synthetic_chart(tmp_path, official, case)
    before = path.read_bytes()
    mesh = mesh_function()(topology)
    check_mesh(mesh, source)
    assert mesh.nodes.dtype == np.dtype(np.longdouble)
    assert path.read_bytes() == before


def test_hanging_cuts_keep_narrow_span_quads_aligned(tmp_path, official):  # noqa: F811
    topology, source, path = synthetic_chart(tmp_path, official, 'span_precision')
    before = path.read_bytes()
    mesh = mesh_function()(topology)
    check_mesh(mesh, source)
    lo, hi = np.longdouble(.5), np.longdouble(np.nextafter(.5, 1.))
    uv = mesh.element_parameters
    narrow = (uv[:, :, 0].min(axis=1) >= lo) & (uv[:, :, 0].max(axis=1) <= hi)
    assert narrow.any()
    for chart, indices in zip(uv[narrow], mesh.elements[narrow], strict=True):
        # No diagonal fan across different contour locations: both span levels
        # carry the same cuts, including the original hanging source cut.
        assert len(np.unique(chart[:, 0])) == 2
        assert len(np.unique(chart[:, 1])) == 2
        xyz = mesh.nodes[indices]-mesh.nodes[indices[0]]
        normal = np.cross(xyz[1], xyz[3])
        for r in (-1, 0, 1):
            for s in (-1, 0, 1):
                dr = np.array([-(1-s), 1-s, 1+s, -(1+s)], dtype=np.longdouble)/4
                ds = np.array([-(1-r), -(1+r), 1+r, 1-r], dtype=np.longdouble)/4
                assert np.dot(np.cross(dr@xyz, ds@xyz), normal) > 0
    assert path.read_bytes() == before


def test_duplicate_crossing_pairs_mesh_without_erasing_source_cuts(tmp_path, official):
    from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
    seed, source, path = synthetic_chart(tmp_path, official, 'rectangle')
    source['components']['blade']['structure']['anchors'] = [
        {'name': 'diagonal', 'start_nd_arc': {'grid': [0., 1.], 'values': [0., 1.]}},
        {'name': 'full', 'start_nd_arc': {'grid': [0., 1.], 'values': [.3, .3]}},
        {'name': 'partial', 'start_nd_arc': {'grid': [.1, .8], 'values': [.3, .3]}},
    ]
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    topology = WindIO2Topology(source, hashlib.sha256(before).hexdigest(), seed._definition)
    assert [e for e in topology.event_stations if .25 < e < .35] == [.3]
    mesh = mesh_function()(topology)
    check_mesh(mesh, source)
    uv = mesh.element_parameters-mesh.element_parameters[:, :1]
    twice_area = np.sum(uv[:, :, 0]*np.roll(uv[:, :, 1], -1, axis=1)
                        - uv[:, :, 1]*np.roll(uv[:, :, 0], -1, axis=1), axis=1)
    assert np.all(twice_area > 0)
    np.testing.assert_allclose(twice_area.sum()/2, 1., atol=32*np.finfo(float).eps, rtol=0)
    assert path.read_bytes() == before
    assert topology.source_sha256 == hashlib.sha256(before).hexdigest()


@pytest.mark.parametrize('size', [.1, .3, float('nan')])
def test_reject_nonvalidation_mesh_size(size, official):
    generate = mesh_function()
    with pytest.raises(ValueError, match='element_size must be in'):
        generate(official[0].windio2_topology, element_size=size)
