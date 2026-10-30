"""S13 public geometry-only boundary; no structural or physics acceptance."""
import hashlib

import numpy as np
import pytest

from aeroelast.core.mesh.entities import MeshElement, Node
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.model import build_rust_properties
from tests.software.mesh.test_windio2_official import OFFICIAL
from tests.software.mesh.test_windio2_topology import official  # noqa: F401
from tests.software.mesh.test_windio2_laminate import expected


PATCH = (.08163265306122448, .08173265306122448)


def test_default_v2_rejected_before_mesh_mutation(capsys, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('default v2 must reject before legacy station planning')
    monkeypatch.setattr(BladeMesh, '_insert_row_stations', forbidden)
    gen = BladeMesh(str(OFFICIAL), element_size=.25)
    counters = Node._id_counter, MeshElement._id_counter
    before = OFFICIAL.read_bytes()
    with pytest.raises(ValueError, match='windIO 2.0 requires geometry_only=True; physical FE offsets are not supported'):
        gen.generate(verbose=False)
    assert (Node._id_counter, MeshElement._id_counter) == counters
    assert gen._numad_mesh is None
    assert OFFICIAL.read_bytes() == before
    captured = capsys.readouterr()
    assert captured.out == '' and captured.err == ''


@pytest.mark.parametrize('renumber', [None, 'rcm'])
def test_public_geometry_retains_source_offsets_and_connectivity(official, monkeypatch, renumber):  # noqa: F811 -- imported pytest fixture
    from aeroelast.models.blade.numad.objects.blade import Blade
    import aeroelast.models.blade.numad.mesh_gen as legacy

    def forbidden(*args, **kwargs):
        pytest.fail('geometry-only v2 must not enter the legacy mesher or TE expansion')

    monkeypatch.setattr(Blade, 'expand_blade_geometry_te', forbidden)
    monkeypatch.setattr(legacy, 'get_shell_mesh', forbidden)
    _, source, digest = official
    gen = BladeMesh(str(OFFICIAL), element_size=.25)
    mesh = gen.generate(verbose=False, renumber=renumber, geometry_only=True, eta_bounds=PATCH)
    metadata = mesh.windio2_metadata
    assert metadata.structural_ready is False
    assert metadata.source_sha256 == digest == hashlib.sha256(OFFICIAL.read_bytes()).hexdigest()
    assert len(metadata.elements) == mesh.elements_count > 0
    assert gen.numad_mesh_data['structural_ready'] is False
    assert 'sections' not in gen.numad_mesh_data  # no centered-property surrogate
    # Object-keyed provenance must survive renumbering, without stale integer IDs.
    for element in mesh.elements:
        record = metadata.elements[element]
        uv, basis, cell, stack = record.parameters, record.basis, record.cell, record.stack
        eta, q = uv.mean(axis=0)
        want = expected(source, eta, cell.surface, q, q)
        assert [(l.name, l.material_name, l.fiber_angle_deg) for l in cell.layers] == [(r[0], r[1], r[3]) for r in want]
        np.testing.assert_allclose([l.thickness_m for l in cell.layers], [r[2] for r in want], rtol=0, atol=1e-16)
        assert {l.source.name for l in stack.layers} == {r[0] for r in want}
        xyz = np.array([node.coords for node in element.nodes])
        assert len(set(element.nodes)) == 4
        np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=1e-13)
        for r in (-1., 0., 1.):
            for s in (-1., 0., 1.):
                dr = np.array([-(1-s), 1-s, 1+s, -(1+s)])/4
                ds = np.array([-(1-r), -(1+r), 1+r, 1-r])/4
                assert np.cross(dr @ xyz, ds @ xyz) @ basis[2] > 0
        for layer in stack.layers:
            lo, hi = layer.z_bounds_m
            assert hi-lo == pytest.approx(layer.source.thickness_m, abs=1e-16)
            f = np.asarray(layer.fiber_vector)
            assert np.linalg.norm(f) == pytest.approx(1.)
            assert f @ basis[2] == pytest.approx(0., abs=1e-13)
            a = np.deg2rad(layer.angle_deg)
            np.testing.assert_allclose(f, np.cos(a)*basis[0]+np.sin(a)*basis[1], atol=1e-13)
        if cell.surface == 'skin':
            assert stack.reference == 'outer_mold_line'
            assert stack.midplane_offset_m == pytest.approx(-sum(r[2] for r in want)/2)
            assert max((l.z_bounds_m[1] for l in stack.layers), default=0.) == pytest.approx(0.)
    # Independent coordinates at a native source station on the first ring.
    # Check actual public nodes, not copied raw geometry arrays.
    from scipy.interpolate import PchipInterpolator
    blade, _, _ = official
    g = blade.geometry
    k = 4
    assert g.ispan[k]/blade.definition.span[-1] == pytest.approx(PATCH[0])
    arc_samples, indices = np.unique(g.arclength[:, k], return_index=True)
    interpolator = PchipInterpolator(arc_samples, g.coordinates[indices, :, k])
    checked = 0
    for element, record in metadata.elements.items():
        if record.cell.surface != 'skin':
            continue
        for node, (eta, q) in zip(element.nodes, record.parameters, strict=True):
            if eta == PATCH[0]:
                arc = g.arclength[0, k]+(1-q)*np.ptp(g.arclength[:, k])
                np.testing.assert_allclose(node.coords, interpolator(arc), atol=2e-12, rtol=0)
                checked += 1
    assert checked > 0
    assert np.isfinite(mesh.coords_array).all()
    assert set(metadata.elements) == set(mesh.elements)
    assert all(node in mesh.nodes for chain in metadata.attachments.values() for node in chain)
    before = mesh.coords_array.copy()
    counters = Node._id_counter, MeshElement._id_counter
    with pytest.raises(ValueError, match='geometry-only.*physical FE offsets'):
        build_rust_properties(gen.numad_mesh_data)
    np.testing.assert_array_equal(mesh.coords_array, before)
    assert counters == (Node._id_counter, MeshElement._id_counter)


def test_property_guard_before_material_lookup():
    with pytest.raises(ValueError, match='geometry-only.*physical FE offsets'):
        build_rust_properties({'structural_ready': False})


@pytest.mark.parametrize('options', [{'skin_reference': 'mid'}, {'include_webs': False}])
def test_unsupported_geometry_options_fail_before_entities(options):
    gen = BladeMesh(str(OFFICIAL), element_size=.25, **options)
    counters = Node._id_counter, MeshElement._id_counter
    with pytest.raises(ValueError, match='geometry-only'):
        gen.generate(verbose=False, geometry_only=True, eta_bounds=PATCH)
    assert counters == (Node._id_counter, MeshElement._id_counter)
    assert gen._numad_mesh is None


@pytest.mark.parametrize('options', [{}, {
    'n_samples': 150, 'refine_tip': False, 'span_grading': 'uniform',
    'airfoil_spacing': 'constant', 'max_bay_length': .2, 'laminate_per_row': False,
}, {'airfoil_spacing': 'linear'}])
def test_v1_public_legacy_deduplication_and_properties_preserved(monkeypatch, options):
    import aeroelast.models.blade.numad.mesh_gen as legacy
    from tests.software.mesh.test_windio2_official import ROOT

    # Small deterministic legacy output: exercise the unchanged public adapter,
    # including duplicate-node welding, sets and the centered property consumer.
    raw = {
        'nodes': np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.],
                           [0., 1., 0.], [0., 0., 0.]]),
        'elements': [[4, 1, 2, 3]],
        'sets': {'element': [{'name': 'skin', 'labels': [0]}],
                 'node': [{'name': 'root', 'labels': [0, 4]}]},
        'sections': [{'elementSet': 'skin', 'layup': [['iso', .002, 0.]]}],
        'materials': [{'name': 'iso', 'density': 1000.,
                       'elastic': {'E': 1e9, 'nu': .3}}],
    }
    calls = []
    def shell(*args, **kwargs):
        assert kwargs['spanGrading'] == options.get('span_grading', 'chord')
        return raw
    monkeypatch.setattr(legacy, 'get_shell_mesh', shell)
    monkeypatch.setattr(BladeMesh, '_insert_row_stations', lambda *a, **k: calls.append('row'))
    monkeypatch.setattr(BladeMesh, '_refine_high_gradient_sections', lambda *a, **k: calls.append('tip'))
    monkeypatch.setattr(BladeMesh, '_limit_bay_length', lambda self, length, **k: calls.append(length))
    gen = BladeMesh(str(ROOT/'tests/IEA-15-240-RWT.yaml'), element_size=.25, **options)
    effective = {'n_samples': 300, 'refine_tip': True, 'span_grading': 'chord',
                 'airfoil_spacing': 'cosine', 'max_bay_length': None, 'laminate_per_row': True}
    effective.update(options)
    assert {name: getattr(gen, name) for name in effective} == effective
    if options.get('airfoil_spacing') == 'linear':
        # V1 already rejects this spelling in its existing resampler; preserve
        # that error rather than inventing support while fixing the v2 boundary.
        counters = Node._id_counter, MeshElement._id_counter
        with pytest.raises(Exception) as error:
            gen.generate(verbose=False, renumber='simple')
        assert type(error.value) is Exception
        assert str(error.value) == 'Resampling algorithm specified is not an available option'
        assert counters == (Node._id_counter, MeshElement._id_counter)
        assert gen._numad_mesh is None and not calls
        return
    mesh = gen.generate(verbose=False, renumber='simple')
    assert calls == ([.2] if options else ['tip', 'row'])
    assert mesh.windio2_metadata is None
    assert mesh.node_count == 4 and mesh.elements_count == 1
    assert mesh.elements[0].node_ids == (0, 1, 2, 3)
    assert len(mesh.node_sets['root'].nodes) == 1
    assert gen.numad_mesh_data is raw
    assert build_rust_properties(raw)['skin'] == {
        'type': 'isotropic', 'e': 1e9, 'nu': .3, 'rho': 1000.,
        'thickness': .002, 'shear_correction': 5./6.,
    }


@pytest.mark.parametrize('option,value', [
    ('n_samples', 150), ('refine_tip', False), ('span_grading', 'uniform'),
    ('airfoil_spacing', 'linear'), ('max_bay_length', .2),
    ('laminate_per_row', False),
    # Explicit legacy defaults are still requests for unsupported controls.
    ('n_samples', 300), ('refine_tip', True), ('span_grading', 'chord'),
    ('airfoil_spacing', 'cosine'), ('max_bay_length', None),
    ('laminate_per_row', True),
])
def test_explicit_legacy_controls_rejected_before_raw(option, value, monkeypatch, capsys):
    import aeroelast.core.mesh.windio2 as adapter

    calls = []
    def forbidden(*args, **kwargs):
        calls.append(True)
        pytest.fail('unsupported control reached the raw mesher')
    monkeypatch.setattr(adapter, 'generate_mesh', forbidden)
    gen = BladeMesh(str(OFFICIAL), element_size=.25, **{option: value})
    counters = Node._id_counter, MeshElement._id_counter
    before = OFFICIAL.read_bytes()
    with pytest.raises(ValueError) as error:
        gen.generate(geometry_only=True, eta_bounds=PATCH, verbose=False)
    assert type(error.value) is ValueError
    assert str(error.value) == f"windIO 2.0 geometry-only does not support explicit option '{option}'"
    assert not calls
    assert gen._numad_mesh is None
    assert counters == (Node._id_counter, MeshElement._id_counter)
    assert OFFICIAL.read_bytes() == before
    output = capsys.readouterr()
    assert output.out == output.err == ''


@pytest.mark.parametrize('options,want', [({}, 300), ({'n_samples': 150}, 150)])
def test_v1_convenience_resolves_samples_without_sentinel(monkeypatch, options, want):
    from tests.software.mesh.test_windio2_official import ROOT

    def generate(self, **kwargs):
        return self.n_samples, self.refine_tip, self.span_grading, self.airfoil_spacing, self.max_bay_length, self.laminate_per_row, kwargs
    monkeypatch.setattr(BladeMesh, 'generate', generate)
    result = BladeMesh.create_from_yaml(str(ROOT/'tests/IEA-15-240-RWT.yaml'),
                                       element_size=.25, **options)
    assert result == (want, True, 'chord', 'cosine', None, True,
                      {'renumber': None, 'geometry_only': False, 'eta_bounds': (0., 1.)})


@pytest.mark.parametrize('n_samples', [150, 300])
def test_convenience_rejects_explicit_samples(n_samples, monkeypatch):
    import aeroelast.core.mesh.windio2 as adapter

    def forbidden(*args, **kwargs):
        pytest.fail('convenience entrypoint reached raw mesher with unsupported control')
    monkeypatch.setattr(adapter, 'generate_mesh', forbidden)
    counters = Node._id_counter, MeshElement._id_counter
    before = OFFICIAL.read_bytes()
    with pytest.raises(ValueError) as error:
        BladeMesh.create_from_yaml(str(OFFICIAL), element_size=.25,
                                   geometry_only=True, eta_bounds=PATCH, n_samples=n_samples)
    assert type(error.value) is ValueError
    assert str(error.value) == "windIO 2.0 geometry-only does not support explicit option 'n_samples'"
    assert counters == (Node._id_counter, MeshElement._id_counter)
    assert OFFICIAL.read_bytes() == before


def test_geometry_convenience_entrypoint():
    mesh = BladeMesh.create_from_yaml(str(OFFICIAL), element_size=.25,
                                     geometry_only=True, eta_bounds=PATCH)
    assert mesh.windio2_metadata.structural_ready is False
    assert mesh.elements_count == len(mesh.windio2_metadata.elements)
