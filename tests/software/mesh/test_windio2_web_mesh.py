"""Public raw skin/web connectivity and explicit shell-coordinate contracts."""
from collections import Counter
import copy
import hashlib

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.objects.blade import Blade
from aeroelast.models.blade.numad.objects.windio2_topology import WindIO2Topology
from tests.software.mesh.test_windio2_laminate import expected
from tests.software.mesh.test_windio2_topology import official  # noqa: F401


def generate(*args, **kwargs):
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh
    return generate_mesh(*args, **kwargs)


@pytest.fixture(scope='module')
def short(tmp_path_factory, official):
    _, original, _ = official
    source = copy.deepcopy(original)
    blade = source['components']['blade']
    for axis in 'xyz':
        f = blade['reference_axis'][axis]
        f['values'] = [.5*q if axis == 'z' else 0. for q in f['grid']]
    for key, value in [('chord', 1.), ('section_offset_y', 0.), ('twist', 0.), ('rthick', 1.)]:
        f = blade['outer_shape'][key]
        f['values'] = [value]*len(f['grid'])
    path = tmp_path_factory.mktemp('web')/'geometry.yaml'
    path.write_text(yaml.safe_dump(source))
    definition = Blade().read_yaml(str(path)).definition
    def field(q):
        return {'grid': [0., 1.], 'values': [q, q]}
    base = copy.deepcopy(blade['structure']['layers'][0])
    base.update(start_nd_arc=field(0.), end_nd_arc=field(1.))
    # The independent source oracle reads each web's local anchor list.
    web = {'name': 'web', 'start_nd_arc': field(.25), 'end_nd_arc': field(.75),
           'anchors': []}
    # Two distinct webs share their start constraint exactly, not a spatial weld.
    second = dict(web, name='second', end_nd_arc=field(.625))
    layers = [base]
    for name in ('web', 'second'):
        for i, (lo, hi, angle) in enumerate([(0., 1., 0.), (.25, .75, 31.)]):
            layer = dict(base, name=f'{name}_{i}', web=name,
                         start_nd_arc=field(lo), end_nd_arc=field(hi),
                         fiber_orientation=field(angle), thickness=field(.001*(i+1)))
            if i:
                layer.update(start_nd_grid=.25, end_nd_grid=.75)
            layers.append(layer)
    blade['structure'].update(anchors=[], webs=[web, second], layers=layers)
    path = path.with_name('source.yaml')
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    return WindIO2Topology(source, hashlib.sha256(before).hexdigest(), definition), source, path, before


def edges(elements):
    return Counter(tuple(sorted((a, b))) for el in elements
                   for a, b in zip(el, np.roll(el, -1), strict=True))


def test_web_shared_edges_source_ownership_and_normals(short):
    topology, source, path, before = short
    mesh = generate(topology, element_size=.25)
    skin = edges(mesh.elements[np.array(mesh.surfaces) == 'skin'])
    for name in ('web', 'second'):
        selected = np.array(mesh.surfaces) == name
        web = edges(mesh.elements[selected])
        assert selected.any()
        for handle in ('start_nd_arc', 'end_nd_arc'):
            chain = mesh.attachments[(name, handle)]
            # span 0.5 m at es 0.25: one attachment node per row, edge ~ es
            assert len(chain) >= 3
            for a, b in zip(chain[:-1], chain[1:], strict=True):
                edge = tuple(sorted((a, b)))
                assert skin[edge] == 2 and web[edge] == 1
        for cell, uv in zip(np.array(mesh.cells, dtype=object)[selected], mesh.element_parameters[selected], strict=True):
            eta, q = uv.mean(axis=0)
            want = expected(source, eta, name, q, q)
            assert [(l.name, l.material_name, l.fiber_angle_deg) for l in cell.layers] == [(r[0], r[1], r[3]) for r in want]
            np.testing.assert_allclose([l.thickness_m for l in cell.layers], [r[2] for r in want], atol=1e-16, rtol=0)
    assert mesh.attachments[('web', 'start_nd_arc')] == mesh.attachments[('second', 'start_nd_arc')]
    for xyz, basis in zip(mesh.nodes[mesh.elements], mesh.bases, strict=True):
        np.testing.assert_allclose(basis@basis.T, np.eye(3), atol=1e-14)
        assert np.dot(np.cross(basis[0], basis[1]), basis[2]) > 0
        # Bilinear shell Jacobians, not just one triangle's area.
        for r in (-1., 0., 1.):
            for s in (-1., 0., 1.):
                dr = np.array([-(1-s), 1-s, 1+s, -(1+s)])/4
                ds = np.array([-(1-r), -(1+r), 1+r, 1-r])/4
                assert np.dot(np.cross(dr@xyz, ds@xyz), basis[2]) > 0
    assert mesh.material_angles is None  # geometry-only default; material mapping is opt-in
    assert path.read_bytes() == before


def test_explicit_physical_fiber_vectors_to_adapter_angles():
    from aeroelast.models.blade.numad.mesh_gen.windio2_basis import physical_fiber_angle
    # Rows are physical element x, y, normal. Independent rotated shell plane.
    basis = np.array([[0., 0., 1.], [0., 1., 0.], [-1., 0., 0.]])
    for angle in (0., 31., -47., 90.):
        a = np.deg2rad(angle)
        vector = np.array([0., np.sin(a), np.cos(a)])
        assert physical_fiber_angle(vector, basis) == pytest.approx(angle, abs=2e-14)


def test_official_patch_full_mesh(official):
    blade, source, digest = official
    eta = .08163265306122448
    mesh = generate(blade.windio2_topology, eta_bounds=(eta, eta+.0001))
    names = {w['name'] for w in source['components']['blade']['structure']['webs']
             if w.get('start_nd_grid', 0.) <= eta <= w.get('end_nd_grid', 1.)}
    assert {'skin', *names} == set(mesh.surfaces)
    from tests.software.mesh.test_windio2_official import OFFICIAL
    assert hashlib.sha256(OFFICIAL.read_bytes()).hexdigest() == digest


def test_web_is_meshed_only_where_its_laminate_has_thickness(short):
    """A web whose layers are zero thick on a span range does not exist there (no phantom web)."""
    topology, source, _, _ = short
    source = copy.deepcopy(source)
    for layer in source['components']['blade']['structure']['layers']:
        if layer.get('web') == 'second':
            layer['thickness'] = {'grid': [0., .4, 1.], 'values': [0., 0., layer['thickness']['values'][0]]}
    topology = WindIO2Topology(source, '0' * 64, topology._definition)
    mesh = generate(topology, element_size=.25)
    eta = mesh.element_parameters[:, :, 0].mean(axis=1)
    surfaces = np.array(mesh.surfaces)
    assert (eta[surfaces == 'second'] > .4 - 1e-12).all() and (surfaces == 'second').any()
    assert eta[surfaces == 'web'].min() < .25  # the untouched web still spans the root
