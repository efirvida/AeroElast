"""Public CLI usage: python tools/bd_k16_candidate.py --output ${SCRATCH}/tmp/new-k16.

All runtime products live in TemporaryDirectory, never in reference directories.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / 'tools/bd_k16_candidate.py'
GEOMETRY = ROOT / ('tests/IEA15MW/validation papers/github-IEA-15-240-RWT/'
                   'OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat')
PROPERTIES = ROOT / 'tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat'
# Exercise the public CLI while forbidding process execution and heavy imports.
GUARD = '''import sys, runpy
def audit(event, args):
    if event in ('subprocess.Popen', 'os.system', 'os.posix_spawn'):
        raise RuntimeError('binary execution forbidden')
    if event == 'import' and args[0].split('.')[0] in ('aeroelast', 'numpy', 'scipy'):
        raise RuntimeError('heavy import forbidden')
sys.addaudithook(audit)
sys.argv = sys.argv[1:]
sys.path.insert(0, str(__import__('pathlib').Path(sys.argv[0]).parent))
runpy.run_path(sys.argv[0], run_name='__main__')
'''


def invoke(*args):
    return subprocess.run([sys.executable, '-B', '-c', GUARD, str(CLI), *map(str, args)],
                          capture_output=True, text=True, timeout=10,
                          env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})


def parameter(text, key):
    return re.search(rf'^\s*(\S+)\s+{re.escape(key)}\s', text, re.M)[1]


def snapshot(directory):
    return {str(p.relative_to(directory)): p.read_bytes()
            for p in directory.rglob('*') if p.is_file()}


@pytest.mark.parametrize('case,gravity,tip,distributed', [
    ('LC1', (-.000981, 0, 0), (0, 0, 0), (0, 0, 0)),
    ('LC2', (0, 0, 0), (100, 0, 0), (0, 0, 0)),
    ('LC3', (0, 0, 0), (0, -100, 0), (0, 0, 0)),
    ('LC4', (0, 0, 0), (0, 0, 0), (.5, 0, 0)),
])
def test_candidate(case, gravity, tip, distributed):
    originals = (GEOMETRY.read_bytes(), PROPERTIES.read_bytes())
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / 'candidate'
        result = invoke('--output', target, '--case', case)
        assert (result.returncode, result.stdout, result.stderr) == (0, f'{target}\n', '')
        assert set(snapshot(target)) == {'beam.dat', 'blade.dat', 'driver.dvr', 'candidate.json'}
        beam = (target / 'beam.dat').read_text()
        assert parameter(beam, 'quadrature') == '1'
        assert parameter(beam, 'order_elem') == '60'
        assert parameter(beam, 'BldNd_BlOutNd') == '"All"'
        nodal = beam.split('BldNd_BlOutNd', 1)[1].split('OutList', 1)[1]
        assert re.findall(r'^"(\w+)"$', nodal, re.M) == [
            'TDxr', 'TDyr', 'TDzr', 'RDxr', 'RDyr', 'RDzr']
        assert nodal.strip().startswith('-')  # original schema comments retained
        assert 'PITCH ACTUATOR PARAMETERS' not in beam
        def geometry(text):
            return text.split('GEOMETRY PARAMETER', 1)[1].split('MESH PARAMETER')[0]

        assert geometry(beam) == geometry(GEOMETRY.read_text())
        driver = (target / 'driver.dvr').read_text()
        assert parameter(driver, 'DynamicSolve') == 'False'
        assert parameter(driver, 'GlbRotBladeT0') == 'True'
        assert '! global to root\n1 0 0\n0 1 0\n0 0 1\n' in driver
        for axis in 'xyz':
            assert parameter(driver, f'omega_{axis}') == '0'
        for names, values in [('gravity_', gravity), ('TipF', tip), ('DistrF', distributed)]:
            assert tuple(float(parameter(driver, names + a)) for a in 'xyz') == values
        for prefix in ('TipM', 'DistrM'):
            assert all(parameter(driver, prefix + a) == '0' for a in 'xyz')
        assert parameter(driver, 'NumPointLoads') == '0'
        assert f'"{target / "beam.dat"}" InputFile' in driver
        assert f'"{target / "blade.dat"}"' in beam
        blade = (target / 'blade.dat').read_text()
        assert parameter(blade, 'damp_type') == '0'
        assert '1 n_modes - unused\n0 zeta - unused\n' in blade
        assert blade.split('DISTRIBUTED PROPERTIES', 1)[1] == (
            PROPERTIES.read_text().split('DISTRIBUTED PROPERTIES', 1)[1])
        record = json.loads((target / 'candidate.json').read_text())
        assert record['g1_closed'] is False
        assert record['executed'] is False
        assert record['case'] == case
        assert record['missing_evidence'] == [
            'undeformed motion-node positions and directors',
            'verified source-content equivalence to installed binary',
            'verified node ordering, root/tip coverage, units and rotation conversion',
            'G2 geometric/load equivalence', 'G3 independent controls and uncertainty bounds',
            'G4 shared instrument']
    assert originals == (GEOMETRY.read_bytes(), PROPERTIES.read_bytes())


def test_reference_geometry_request_is_source_supported_not_acquired():
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / 'candidate'
        result = invoke('--output', target)
        assert (result.returncode, result.stdout, result.stderr) == (0, f'{target}\n', '')
        assert parameter((target / 'beam.dat').read_text(), 'SumPrint') == 'True'
        record = json.loads((target / 'candidate.json').read_text())
        assert record['reference_geometry_contract'] == {
            'request': {'SumPrint': True},
            'status': 'source-supported; not acquired or validated',
            'source_evidence': 'parent source review: BeamDyn_IO.f90:2049-2060',
            'summary_suffix': '.sum.yaml',
            'array_schema': 'Init_Nodes_E<n>',
            'columns_in_order': ['initial X', 'initial Y', 'initial Z',
                                 'initial WM parameter 1', 'initial WM parameter 2',
                                 'initial WM parameter 3'],
            'coordinate_system': 'IEC',
            'indexing': 'per element <n>, local FE/GLL nodes; no global indexing',
            'serialization': 'ES18.5',
            'caveats': [
                'ordering, root/tip coverage and units unverified',
                'WM-to-director conversion unverified',
                'finite serialization precision; independent propagated bound missing'],
        }
        assert record['g1_closed'] is False and record['executed'] is False
        assert 'undeformed motion-node positions and directors' in record['missing_evidence']
        assert set(snapshot(target)) == {'beam.dat', 'blade.dat', 'driver.dvr', 'candidate.json'}
    result = invoke('--help')
    assert result.returncode == 0 and result.stderr == ''
    for phrase in ('SumPrint=True', 'Init_Nodes_E<n>', 'ES18.5', '${SCRATCH}/tmp'):
        assert phrase in result.stdout


def test_overwrite_rejected_without_any_changes():
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / 'candidate'
        assert invoke('--output', target).returncode == 0
        before = snapshot(Path(scratch))
        result = invoke('--output', target)
        assert (result.returncode, result.stdout, result.stderr) == (
            2, '', f'error: output already exists: {target}\n')
        assert snapshot(Path(scratch)) == before


@pytest.mark.parametrize('args,message', [
    (['--order', '9'], 'order must be >=10'),
    (['--scale', '0'], 'scale must be positive and finite'),
    (['--scale', 'nan'], 'scale must be positive and finite'),
])
def test_invalid_inputs_leave_destination_absent(args, message):
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / 'candidate'
        result = invoke('--output', target, *args)
        assert (result.returncode, result.stdout, result.stderr) == (2, '', f'error: {message}\n')
        assert snapshot(Path(scratch)) == {}
        assert not target.exists()


def test_reference_tree_is_protected_even_through_symlink():
    with tempfile.TemporaryDirectory() as scratch:
        alias = Path(scratch) / 'alias'
        alias.symlink_to(ROOT / 'tests/IEA15MW', target_is_directory=True)
        result = invoke('--output', alias / 'bd_static' / 'k16-forbidden')
        assert (result.returncode, result.stdout, result.stderr) == (
            2, '', 'error: output must be outside tests/IEA15MW reference tree\n')
        assert not (alias / 'bd_static' / 'k16-forbidden').exists()


def test_help_and_required_destination():
    result = invoke('--help')
    assert result.returncode == 0 and result.stderr == ''
    for phrase in ('--output', 'generate-only', 'G1 is not closed', 'no overwrite'):
        assert phrase in result.stdout
    result = invoke()
    assert result.returncode == 2 and result.stdout == ''
    assert result.stderr.endswith('error: the following arguments are required: --output\n')
