"""Isolated OpenFAST v5 blade modal campaign (no aeroelast imports).

Usage (from repository root):
  python tools/of_bd_blade_modal.py build
  python tools/of_bd_blade_modal.py submit --case ed_0
  python tools/of_bd_blade_modal.py post --case ed_0

OpenFAST runs are serial SLURM jobs. Inputs live in tests/IEA15MW/of_bd_modal;
linearizations and unavoidable header-only OpenFAST output live in logs/check/of_bd.
No original deck is modified. Existing linearizations are never overwritten by
submit; rerun completed cases in a fresh checkout. Use --case to select a case,
and --attempt N to retry an initialization failure while preserving its logs.
post uses openfast_toolbox's MBC3, then extracts
its collective invariant block; it records off-block coupling as a diagnostic.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/scratch/leahk/eduardo.donestevez')
SOURCE = BASE / 'ofruns/OpenFAST'
BINARY = BASE / 'tools/openfast-env/bin/openfast'
PYTHON = BASE / 'venv/bin/python'
INPUTS = ROOT / 'tests/IEA15MW/of_bd_modal'
LOGS = ROOT / 'logs/check/of_bd'
CASES = [(f'{kind}_{rpm:g}', kind, rpm, n) for kind, n in
         [('ed', 0), ('bd15', 15), ('bd30', 30)] for rpm in [0.0, 7.56]]
# Separate numerical-stability attempts preserve the original failed decks.
CASES += [('bd15_7.56_dt001', 'bd15', 7.56, 15)]
CASES += [('bd15_7.56_loose', 'bd15', 7.56, 15)]


def replace(text, key, value):
    """Replace exactly one OpenFAST labelled parameter, preserving comments."""
    value_pattern = r'[0-9., \t]+' if key == 'LinTimes' else r'[^\s]+'
    pattern = rf'^[ \t]*{value_pattern}([ \t]+{re.escape(key)}[ \t]+[^\n]*)$'
    text, count = re.subn(pattern, lambda m: str(value) + m[1], text, flags=re.M)
    if count != 1:
        raise ValueError(f'{key}: expected exactly one parameter, got {count}')
    return text


def empty_outputs(text):
    return re.sub(r'(\s+OutList\s+[^\n]*\n).*?(^END[^\n]*)',
                  r'\1\2', text, flags=re.M | re.S)


def build(selected):
    import numpy as np
    from scipy.interpolate import CubicSpline

    mono = SOURCE / 'IEA-15-240-RWT-Monopile'
    blade = SOURCE / 'IEA-15-240-RWT'
    fst0 = (mono / 'IEA-15-240-RWT-Monopile.fst').read_text()
    ed0 = (mono / 'IEA-15-240-RWT-Monopile_ElastoDyn.dat').read_text()
    bd0 = (blade / 'IEA-15-240-RWT_BeamDyn.dat').read_text()
    geometry = bd0.split('(deg)\n', 1)[1].split('---------------------- MESH')[0]
    points = np.array([[float(v) for v in line.split()] for line in geometry.splitlines()
                       if line.strip()])
    LOGS.mkdir(parents=True, exist_ok=True)
    INPUTS.mkdir(parents=True, exist_ok=True)
    shared = {
        'ed_blade.dat': (blade / 'IEA-15-240-RWT_ElastoDyn_blade.dat').read_text(),
        'tower.dat': (mono / 'IEA-15-240-RWT-Monopile_ElastoDyn_tower.dat').read_text(),
        'bd_blade.dat': replace((ROOT / 'tests/IEA15MW/reference/'
                                'IEA-15-240-RWT_BeamDyn_blade.dat').read_text(),
                               'damp_type', 0),
    }
    # v5 reads modal-damping fields even when physical damping is disabled.
    shared['bd_blade.dat'] = shared['bd_blade.dat'].replace(
        ' ---------------------- DISTRIBUTED PROPERTIES',
        ' ---------------------- MODAL DAMPING --------------------------------------\n'
        '1 n_modes - Number of modal damping coefficients (unused)\n'
        '0 zeta - Modal damping coefficients (unused)\n'
        ' ---------------------- DISTRIBUTED PROPERTIES')
    for filename, text in shared.items():
        (INPUTS / filename).write_text(text)
    for name, kind, rpm, n in CASES:
        if selected and name != selected:
            continue
        case = INPUTS / name
        case.mkdir(exist_ok=True)
        out = LOGS / name
        out.mkdir(exist_ok=True)
        stable_attempt = name.endswith(('_dt001', '_loose'))
        t0 = 0.1 if stable_attempt else 30.0
        times = [1.0] if rpm == 0 else [t0 + k * 60 / rpm / 12 for k in range(12)]
        tmax = round(times[-1] + 0.01, 3)
        fst = fst0
        params = {'TMax': tmax, 'DT': 0.001 if stable_attempt else 0.005, 'CompElast': 1 if kind == 'ed' else 2,
                      'CompInflow': 0, 'CompAero': 0, 'CompServo': 0, 'CompSeaSt': 0,
                      'CompHydro': 0, 'CompSub': 0, 'CompMooring': 0, 'CompIce': 0, 'CompSoil': 0,
                      'Gravity': 0, 'Linearize': 'True', 'CalcSteady': 'False',
                      'NLinTimes': len(times), 'LinTimes': ', '.join(f'{t:.6f}' for t in times),
                      'LinInputs': 0, 'LinOutputs': 0, 'LinOutJac': 'False', 'LinOutMod': 'False',
                      'OutFileFmt': 1, 'TStart': tmax, 'DT_Out': 100, 'SumPrint': 'False',
                      'OutFmt': '"ES24.16E3"', 'EDFile': f'"{case / "ed.dat"}"'}
        if name.endswith('_loose'):
            params['ModCoupling'] = 1
        for key, value in params.items():
            fst = replace(fst, key, value)
        for i in range(1, 4):
            fst = replace(fst, f'BDBldFile({i})', f'"{case / "bd.dat"}"')
        # OpenFAST writes beside the primary .fst: keep generated files in logs.
        (case / 'case.fst').write_text(fst)
        (out / 'case.fst').write_text(fst)
        ed = empty_outputs(ed0)
        for key in re.findall(r'^\s*(?:True|False)\s+(\w*DOF\w*)\s', ed, re.M):
            ed = replace(ed, key, 'True' if key in ('FlapDOF1', 'FlapDOF2', 'EdgeDOF') else 'False')
        for key, value in {'RotSpeed': rpm, 'SumPrint': 'False',
                               'TwrFile': f'"{INPUTS / "tower.dat"}"'}.items():
            ed = replace(ed, key, value)
        for i in range(1, 4):
            ed = replace(ed, f'BlPitch({i})', 0)
            ed = replace(ed, f'BldFile{i}', f'"{INPUTS / "ed_blade.dat"}"')
        (case / 'ed.dat').write_text(ed)
        if n:
            bd = empty_outputs(bd0)
            # v5 removed the old BeamDyn pitch-actuator input section.
            bd = re.sub(r'---------------------- PITCH ACTUATOR PARAMETERS.*?'
                        r'(?=---------------------- OUTPUTS)', '', bd, flags=re.S)
            for key, value in {'rhoinf': 0.8 if stable_attempt else 0,
                               'order_elem': 4, 'quadrature': 1, 'member_total': n,
                                   'kp_total': 2*n+1, 'RotStates': 'True',
                                   'BldFile': f'"{INPUTS / "bd_blade.dat"}"'}.items():
                bd = replace(bd, key, value)
            newpoints = CubicSpline(points[:, 2], points, axis=0)(
                np.linspace(points[0, 2], points[-1, 2], 2*n+1))
            start = bd.index('     1     50')
            stop = bd.index('---------------------- MESH')
            geom = ''.join(f'{i+1} 3\n' for i in range(n))
            geom += 'kp_xr kp_yr kp_zr initial_twist\n(m) (m) (m) (deg)\n'
            geom += ''.join(' '.join(f'{v:.12e}' for v in row) + '\n' for row in newpoints)
            bd = bd[:start] + geom + bd[stop:]
            (case / 'bd.dat').write_text(bd)
        print(name, case / 'case.fst', 'output', out)


def submit(selected, attempt):
    for name, _, _, _ in CASES:
        if selected and name != selected:
            continue
        out = LOGS / name
        if any(out.glob('case.*.lin')):
            raise FileExistsError(f'{name}: linearizations already exist; rerun in a fresh checkout')
        jobfile = out / ('job.json' if attempt == 1 else f'job-{attempt}.json')
        logfile = out / ('openfast.log' if attempt == 1 else f'openfast-{attempt}.log')
        if jobfile.exists():
            raise FileExistsError(f'{name}: job already recorded; use a new explicit campaign for reruns')
        command = f'cd {shlex.quote(str(ROOT))} && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '
        command += f'{shlex.quote(str(BINARY))} {shlex.quote(str(out / "case.fst"))}'
        command += f' > {shlex.quote(str(logfile))} 2>&1'
        args = ['sbatch', '-p', 'sequana_cpu', '--ntasks=1', '--time=120',
                '-o', str(LOGS / f'{name}-%j.out'), '--wrap', command]
        result = subprocess.run(args, check=True, text=True, capture_output=True)
        record = {'case': name, 'job': result.stdout.strip(), 'command': args}
        jobfile.write_text(json.dumps(record, indent=2) + '\n')
        print(record)


def post(name, attempt):
    import numpy as np
    from openfast_toolbox.linearization.mbc import fx_mbc3

    out = LOGS / name
    files = sorted(out.glob('case.*.lin'), key=lambda p: int(p.stem.split('.')[-1]))
    expected = 1 if name.endswith('_0') else 12
    if len(files) != expected:
        raise RuntimeError(f'{name}: expected {expected} lin files, found {len(files)}')
    logfile = out / ('openfast.log' if attempt == 1 else f'openfast-{attempt}.log')
    if 'OpenFAST terminated normally.' not in logfile.read_text():
        raise RuntimeError(f'{name}: OpenFAST did not terminate normally')
    mbc, data = fx_mbc3([str(f) for f in files], verbose=False)
    if not mbc['performedTransformation']:
        raise RuntimeError('MBC3 transformation not performed')
    triplets = np.asarray(data['RotTripletIndicesStates2'], dtype=int)
    # fx_mbc3 stores transformed [collective, cosine, sine] in each original triplet.
    q = triplets[:, 0]
    if data['ndof1'] or 3*len(q) != data['ndof2']:
        raise RuntimeError('Expected blade-only second-order rotating states')
    idx = np.concatenate([q, q + data['ndof2']])
    other = np.setdiff1d(np.arange(mbc['AvgA'].shape[0]), idx)
    a = mbc['AvgA'][np.ix_(idx, idx)]
    eigenvalues, vectors = np.linalg.eig(a)
    modes = []
    for j in np.argsort(np.abs(eigenvalues)):
        ev = eigenvalues[j]
        if ev.imag <= 0:
            continue
        amplitudes = np.abs(vectors[:len(q), j])
        strongest = np.argsort(amplitudes)[-6:][::-1]
        modes.append({'hz': float(abs(ev)/(2*np.pi)), 'damped_hz': float(ev.imag/(2*np.pi)),
                          'damping': float(-ev.real/abs(ev)),
                          'dominant': [{'state': data['DescStates'][q[k]],
                                         'amplitude': float(amplitudes[k]/amplitudes.max())}
                                    for k in strongest]})
    record = {'case': name, 'binary': str(BINARY), 'deck': str(INPUTS / name / 'case.fst'),
                  'lin_files': [str(p) for p in files], 'rpm': float(mbc['RotSpeed_rpm']),
                  'azimuth_deg': np.asarray(data['Azimuth']).tolist(),
                  'collective_offblock_relative': float(np.linalg.norm(mbc['AvgA'][np.ix_(other, idx)]) /
                                                     np.linalg.norm(a)), 'modes': modes,
                  'all_mbc_hz': sorted(np.asarray(mbc['eigSol']['NaturalFreqs_Hz']).ravel().tolist())}
    (out / 'modes.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', choices=['build', 'submit', 'post'])
    parser.add_argument('--case', choices=[c[0] for c in CASES])
    parser.add_argument('--attempt', type=int, default=1,
                        help='Attempt number; preserves earlier job records and logs (default: 1)')
    args = parser.parse_args()
    if args.action == 'build':
        build(args.case)
    elif args.action == 'submit':
        submit(args.case, args.attempt)
    elif not args.case:
        parser.error('post requires --case')
    else:
        post(args.case, args.attempt)


if __name__ == '__main__':
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
    sys.exit(main())
