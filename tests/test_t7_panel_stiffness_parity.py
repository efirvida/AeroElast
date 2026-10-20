"""Public mesh/deck invariants for the bounded T7(iv) diagnostic."""
import contextlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.model import build_rust_properties
from tools import t7_panel_stiffness_parity as t7


@pytest.fixture
def tmp_path():
    """Do not retain CCX deck scratch after these tests finish."""
    with TemporaryDirectory(prefix='t7-panel-test-') as directory:
        yield Path(directory)


def patch(nx=2, ny=2, blade_plane=False):
    mesh = MeshModel()
    grid = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            node = Node([i / nx, 0., j / ny] if blade_plane else [i / nx, j / ny, 0.])
            mesh.add_node(node)
            grid[i, j] = node
    for j in range(ny):
        for i in range(nx):
            mesh.add_element(MeshElement([grid[i, j], grid[i+1, j],
                                         grid[i+1, j+1], grid[i, j+1]], ElementType.quad))
    mesh.add_element_set(ElementSet('panel', mesh.elements))
    return mesh


def properties():
    return build_rust_properties({'materials': [{'name': 'steel', 'density': 7800.,
        'elastic': {'E': 210e9, 'nu': 0.}}],
        'sections': [{'elementSet': 'panel', 'layup': [['steel', .01, 0.]]}]})


def test_extraction_topology_fixation_and_force(tmp_path):
    source = patch()
    original = source.coords_array.copy()
    mesh, fixed, force = t7.prepare_patch(source, 'panel', span=(0., 1., 0.))
    assert len(mesh.elements) == 4
    np.testing.assert_array_equal(mesh.coords_array, original)
    assert all(a is not b for a, b in zip(source.elements, mesh.elements, strict=True))
    assert all(a is not b for a, b in zip(source.nodes, mesh.nodes, strict=True))
    assert len(t7.boundary_edges(mesh)) == 8
    assert len(fixed) == 3
    np.testing.assert_allclose(force.reshape(-1, 6)[:, :3].sum(axis=0), [0, 0, 1])
    np.testing.assert_array_equal(force.reshape(-1, 6)[fixed], 0.)
    deck = tmp_path / 'panel.inp'
    canonical = t7.write_deck(mesh, properties(), force, deck, span=(0., 1., 0.))
    # Canonical public deck loads are used unchanged in both static solvers.
    repeated = t7.write_deck(mesh, properties(), canonical, deck, span=(0., 1., 0.))
    np.testing.assert_array_equal(canonical, repeated)
    text = deck.read_text()
    assert 'NT7FIXED, 1, 6, 0.0' in text
    assert 'OFFSET' not in text.upper()
    assert '*DLOAD' not in text and 'CENTRIF' not in text
    # Read public deck artifacts, not a private quadratic-mesh storage oracle.
    labels = []
    active = False
    for line in deck.with_suffix('.nam').read_text().splitlines():
        if line.startswith('*'):
            active = line.upper().replace(' ', '').startswith('*NSET,NSET=NT7FIXED')
        elif active:
            labels.extend(int(x) for x in line.split(',') if x.strip())
    assert len(labels) == 5  # three corners and two inboard midsides
    loads = np.zeros((len(mesh.nodes), 3))
    active = False
    for line in text.splitlines():
        if line.startswith('*'):
            active = line.upper() == '*CLOAD'
        elif active and line.strip():
            node, dof, value = line.split(',')
            loads[int(node)-1, int(dof)-1] += float(value)
    np.testing.assert_array_equal(loads, canonical.reshape(-1, 6)[:, :3])
    np.testing.assert_allclose(loads, force.reshape(-1, 6)[:, :3], rtol=1e-6)
    exported = {}
    active = False
    for line in deck.with_suffix('.msh').read_text().splitlines():
        if line.startswith('*'):
            active = line.upper().startswith('*NODE')
        elif active:
            values = line.split(',')
            exported[int(values[0])] = np.array(values[1:], dtype=float)
    for i, node in enumerate(mesh.nodes, 1):
        np.testing.assert_array_equal(exported[i], node.coords)
    np.testing.assert_array_equal(source.coords_array, original)


def test_rejects_missing_or_disconnected_patch():
    with pytest.raises(ValueError, match='Element sets not found'):
        t7.prepare_patch(patch(), 'missing', span=(0, 1, 0))
    mesh = patch(1, 1)
    other = patch(1, 1)
    for n in other.nodes:
        mesh.add_node(n)
    for e in other.elements:
        mesh.add_element(e)
    mesh.add_element_set(ElementSet('both', mesh.elements))
    with pytest.raises(ValueError, match='connected'):
        t7.prepare_patch(mesh, 'both', span=(0, 1, 0))


def test_work_and_public_stiffness_energy():
    mesh, fixed, force = t7.prepare_patch(patch(), 'panel', span=(0, 1, 0))
    result, displacement = t7.solve_mitc4(mesh, properties(), fixed, force, (0, 1, 0))
    assert result['work_J'] > 0
    assert result['free_residual'] < 1e-9
    assert result['energy_relative_error'] < 1e-9
    assert result['force_balance_relative_error'] < 1e-9
    assert result['work_J'] == pytest.approx(.5 * force @ displacement)


def test_cli_rejected_output_is_clean(tmp_path):
    out = io.StringIO()
    err = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        with pytest.raises(SystemExit) as exit_info:
            t7.main(['--out', str(tmp_path / 'forbidden')])
    assert exit_info.value.code == 2
    assert out.getvalue() == ''
    assert err.getvalue().endswith('error: --out must be under logs/check/tipmodes/t7-iv\n')
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('unsupported', [False, True])
def test_panel_report_numpy_section_metadata(tmp_path, monkeypatch, capsys, unsupported):
    """Exercise panel's real report path; replace only expensive external inputs."""
    from aeroelast.core.mesh import generators

    source = patch(blade_plane=True)
    source.add_element_set(ElementSet(t7.PANEL, source.elements))
    section = {'elementSet': t7.PANEL, 'layup': [['steel', np.float64(.01), 0.]],
               'xDir': np.array([0., 0., 1.]),
               'metadata': {'station': np.int64(46), 'active': np.bool_(True),
                            'weight': np.float32(.25), 'grid': np.array([[1, 2], [3, 4]])}}
    if unsupported:
        section['unrelated'] = object()
    data = {'materials': [{'name': 'steel', 'density': 7800.,
                          'elastic': {'E': 210e9, 'nu': 0.}}], 'sections': [section]}

    def synthetic_blade(**kwargs):
        assert kwargs['element_size'] == .25
        assert kwargs['skin_reference'] == 'oml'
        return SimpleNamespace(numad_mesh_data=data, generate=lambda **kw: source)

    # Substitute the mesher factory, not mesh/assembler/writer production methods.
    monkeypatch.setattr(generators, 'BladeMesh', synthetic_blade)
    monkeypatch.setattr(t7, 'solve_ccx', lambda mesh, props, force, deck, ccx:
                        (1., np.zeros((len(mesh.nodes), 3))))
    if unsupported:
        with pytest.raises(TypeError, match='Object of type object is not JSON serializable'):
            t7.panel(tmp_path, 'unused-ccx', {'valid': True})
        assert not (tmp_path / 'result.json').exists()
        assert 'PANEL ' not in capsys.readouterr().out
    else:
        assert t7.panel(tmp_path, 'unused-ccx', {'valid': True}) is None
        report = json.loads((tmp_path / 'result.json').read_text())
        assert report['sections'] == [{'elementSet': t7.PANEL,
            'layup': [['steel', .01, 0.]], 'xDir': [0., 0., 1.],
            'metadata': {'station': 46, 'active': True, 'weight': .25, 'grid': [[1, 2], [3, 4]]}}]
        stdout = capsys.readouterr().out
        printed = stdout.split('PANEL ', 1)[1]
        assert json.loads(printed) == report
        assert (tmp_path / 'panel.npz').exists()
    # Reporting must not normalize or mutate the original section metadata.
    np.testing.assert_array_equal(section['xDir'], [0., 0., 1.])
    assert isinstance(section['metadata']['station'], np.int64)
    assert isinstance(section['metadata']['active'], np.bool_)


def test_cli_help():
    out = io.StringIO()
    err = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        with pytest.raises(SystemExit) as exit_info:
            t7.main(['--help'])
    assert exit_info.value.code == 0
    assert '--reference-only' in out.getvalue()
    assert 'No K_G' in out.getvalue()
    assert err.getvalue() == ''
