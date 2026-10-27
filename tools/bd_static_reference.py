"""BeamDyn v5 static S-2 reference, without importing aeroelast.

python tools/bd_static_reference.py build --case LC2 --order 60 --scale .0001
python tools/bd_static_reference.py run --case LC2 --order 60 --scale .0001
python tools/bd_static_reference.py post --case LC2 --order 60 --scale .0001

Run `run` on a compute node (sbatch -p sequana_cpu --ntasks=1). Build writes
inputs only under tests/IEA15MW/bd_static; run/post write logs/check/bd_static.
--blade iea34 uses the IEA-3.4 deck (tests/IEA34/bd_static, logs/check/bd_static/iea34).
Use distinct order/scale combinations to preserve prior evidence. No overwrite.

Frame: BD (x,y,z) -> shell (-y,x,z). The published prebend in BD x maps to
shell y; positive BD twist rotates chord y toward -x, matching NuMAD's
negative twist about shell z. Root DCM is identity, root motion fixed, omega=0.
Thus BD TipTDxr is shell tip uy (m), NOT TipTDyr. Loads are dead/global loads.
LC1 g=(-9.81,0,0), LC2 tip F=(1e6,0,0), LC3 tip F=(0,-1e6,0), LC4 q=(5e3,0,0),
scaled uniformly. LC3 is the shell edgewise tip load (+1 MN, 0, 0): its shell
tip ux is -TipTDyr. The summary file (SumPrint) is written for the deck mass.
Small-load displacement / scale estimates the linear response; compare at least
two scales before using it. Full LC2 is outside the small-deflection regime
(~82 m linear reference). p40..60 changed normalized tip uy by at most 0.13%
in the initial campaign; do not infer convergence from p15..20 alone.

--nodes 53,54,...,61 (at most 9) also writes the translations of those FE nodes (N#TDxr/yr/zr; with Gauss
quadrature BeamDyn's output nodes are the GLL FE nodes, node order_elem+1 is the tip) and the reference mesh
at init (VTK), whose positions give each node's initial z: the interior-station comparison of l3-closure S12.
Its directories carry a `-nodes` suffix, so tip-only evidence is never overwritten.

Limitations: BD tip is its reference axis, shell tip is the first maximum-z
node, not the same material point. BD q is per reference arc length whereas
shell LC4 is per projected z, lumped into 60 bins. EB endpoints are bin centres
rather than the true root/tip. Keep these differences visible in any verdict.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

sys.dont_write_bytecode = True
from of_bd_blade_modal import replace  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/scratch/leahk/eduardo.donestevez')
GEOMETRY = BASE / 'ofruns/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat'
PROPERTIES = ROOT / 'tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat'
BINARY = BASE / 'tools/openfast-env/bin/beamdyn_driver'
INPUTS = ROOT / 'tests/IEA15MW/bd_static'
LOGS = ROOT / 'logs/check/bd_static'
# Second official blade (l3-closure T10b): same OpenFAST v5 deck chain, its own inputs and logs.
BLADES = {
    'iea15': (GEOMETRY, PROPERTIES, INPUTS, LOGS),
    'iea34': (ROOT / 'tests/IEA34/reference/IEA-3.4-130-RWT_BeamDyn.dat',
              ROOT / 'tests/IEA34/reference/IEA-3.4-130-RWT_BeamDyn_blade.dat',
              ROOT / 'tests/IEA34/bd_static', ROOT / 'logs/check/bd_static/iea34'),
}
CHANNELS = ('TipTDxr', 'TipTDyr', 'TipTDzr')


def case_loads(case, scale):
    """Return gravity, tip force, and distributed force in BD GLOBAL axes."""
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError('scale must be positive and finite')
    index, vector = {'LC1': (0, (-9.81, 0, 0)), 'LC2': (1, (1e6, 0, 0)),
                     'LC3': (1, (0, -1e6, 0)), 'LC4': (2, (5e3, 0, 0))}[case]
    loads = [(0, 0, 0)] * 3
    loads[index] = tuple(v * scale for v in vector)
    return loads


def driver_text(case, scale, beam, vtk=0):
    gravity, tip, distributed = case_loads(case, scale)
    lines = ['BeamDyn driver v5', 'S-2 small-load static reference',
             '--- simulation', 'False DynamicSolve', '0 t_initial', '10 t_final',
             '1 dt', '--- gravity']
    lines += [f'{v:.12g} gravity_{axis}' for axis, v in zip('xyz', gravity, strict=True)]
    lines += ['--- frame', '0 x', '0 y', '0 z', '! DCM rows', '! global to root',
              '1 0 0', '0 1 0', '0 0 1', 'True GlbRotBladeT0', '--- velocity',
              '0 omega_x', '0 omega_y', '0 omega_z', '--- loads']
    for name, vector in [('Distr', distributed), ('Tip', tip)]:
        lines += [f'{v:.12g} {name}{component}' for component, v in
                  zip(('Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'), (*vector, 0, 0, 0), strict=True)]
    lines += ['0 NumPointLoads', '! point load table', '! units', '--- beam',
              f'"{beam}" InputFile', '--- outputs', f'{vtk} WrVTK', '1 VTK_fps']
    return '\n'.join(lines) + '\n'


def node_channels(nodes):
    """Translation channels of the requested output nodes, in OutNd order (N1 is the first listed node)."""
    return tuple(f'N{k}TD{a}r' for k in range(1, len(nodes) + 1) for a in 'xyz')


def parse_output(text, channels=CHANNELS):
    """Parse last finite driver row, requiring all the given translation channels (default: the tip)."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        headers = line.split()
        if headers and headers[0] == 'Time' and all(c in headers for c in channels):
            if i + 1 >= len(lines):
                break
            units = lines[i+1].split()
            if len(units) != len(headers) or any(units[headers.index(c)] != '(m)' for c in channels):
                raise ValueError('expected tip translations in metres')
            rows = [row.split() for row in lines[i+2:] if row.strip()]
            if not rows or len(rows[-1]) != len(headers):
                break
            try:
                values = [float(v.replace('D', 'E')) for v in rows[-1]]
            except ValueError as exc:
                raise ValueError('invalid BeamDyn output row') from exc
            if not all(math.isfinite(v) for v in values):
                raise ValueError('nonfinite BeamDyn output')
            if values[0] <= 0:
                raise ValueError('no solved static output after initialization')
            return {c: values[headers.index(c)] for c in channels}
    raise ValueError('missing BeamDyn translation channels or data')


def run_name(case, order, scale, nodes=()):
    return f'{case}-p{order}-s{scale:g}' + ('-nodes' if nodes else '')


def build(case, order, scale, nodes=()):
    name = run_name(case, order, scale, nodes)
    target = INPUTS / name
    target.mkdir(parents=True, exist_ok=False)
    prop = replace(PROPERTIES.read_text(), 'damp_type', 0)
    prop = prop.replace(' ---------------------- DISTRIBUTED PROPERTIES',
                        ' ---------------------- MODAL DAMPING ---------------------\n'
                        '1 n_modes - unused\n0 zeta - unused\n'
                        ' ---------------------- DISTRIBUTED PROPERTIES')
    (target / 'blade.dat').write_text(prop)
    beam = GEOMETRY.read_text()
    beam = re.sub(r'---------------------- PITCH ACTUATOR PARAMETERS.*?'
                  r'(?=---------------------- OUTPUTS)', '', beam, flags=re.S)
    for key, value in {'order_elem': order, 'quadrature': 1, 'NRMax': 50, 'SumPrint': 'True',
                       'BldFile': f'"{target / "blade.dat"}"',
                       'OutFmt': '"ES20.10E3"'}.items():
        beam = replace(beam, key, value)
    if nodes:
        beam = replace(beam, 'NNodeOuts', len(nodes))
        beam, count = re.subn(r'^[^\n]*?([ \t]+OutNd[ \t]+[^\n]*)$',
                              lambda m: ', '.join(map(str, nodes)) + m[1], beam, flags=re.M)
        if count != 1:
            raise ValueError(f'OutNd: expected exactly one parameter, got {count}')
        # first OutList only (the BldNd section after it has its own)
        beam = re.sub(r'(^[ \t]*OutList[^\n]*\n)', lambda m: m[1] + '"' + ', '.join(node_channels(nodes)) + '"\n',
                      beam, count=1, flags=re.M)
    (target / 'beam.dat').write_text(beam)
    (target / 'driver.dvr').write_text(driver_text(case, scale, target / 'beam.dat', vtk=1 if nodes else 0))
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', choices=('build', 'run', 'post'))
    parser.add_argument('--case', choices=('LC1', 'LC2', 'LC3', 'LC4'), required=True)
    parser.add_argument('--order', type=int, default=60, help='single member polynomial order, >=10 (default: 60)')
    parser.add_argument('--scale', type=float, default=.0001, help='positive load scale (default: .0001); check linearity across scales')
    parser.add_argument('--blade', choices=tuple(BLADES), default='iea15', help='official deck (default: iea15)')
    parser.add_argument('--nodes', default='', help='comma list of FE output nodes, at most 9 (interior station, S12)')
    args = parser.parse_args()
    try:
        nodes = tuple(int(n) for n in args.nodes.split(',') if n.strip())
    except ValueError:
        parser.error('--nodes must be a comma list of integers')
    if len(nodes) > 9 or any(not 1 <= n <= args.order + 1 for n in nodes):
        parser.error(f'--nodes: at most 9 nodes in 1..{args.order + 1}')
    global GEOMETRY, PROPERTIES, INPUTS, LOGS
    GEOMETRY, PROPERTIES, INPUTS, LOGS = BLADES[args.blade]
    case_loads(args.case, args.scale)
    if args.order < 10:
        parser.error('order must be >=10')
    name = run_name(args.case, args.order, args.scale, nodes)
    if args.action == 'build':
        print(build(args.case, args.order, args.scale, nodes))
        return
    output = LOGS / name
    if args.action == 'run':
        output.mkdir(parents=True, exist_ok=False)
        # Driver writes all products beside its .dvr file, not its BeamDyn input.
        driver = output / 'driver.dvr'
        driver.write_text((INPUTS / name / 'driver.dvr').read_text())
        with (output / 'driver.log').open('w') as log:
            subprocess.run([str(BINARY), str(driver)], cwd=output, stdout=log,
                           stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    else:
        log = (output / 'driver.log').read_text()
        if 'terminated normally' not in log:
            raise ValueError('driver did not terminate normally')
        tip = parse_output((output / 'driver.out').read_text())
        record = {'case': args.case, 'order': args.order, 'scale': args.scale,
                  'tip_bd_m': tip, 'linear_shell_uy_m': tip['TipTDxr'] / args.scale,
                  'linear_shell_ux_m': -tip['TipTDyr'] / args.scale,
                  'binary': str(BINARY), 'geometry': str(GEOMETRY), 'properties': str(PROPERTIES)}
        (output / 'result.json').write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps(record))


if __name__ == '__main__':
    main()
