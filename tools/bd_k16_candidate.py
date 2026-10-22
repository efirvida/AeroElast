"""Prepare a generate-only K16 G1 input candidate; never run BeamDyn or postprocess.

Usage: python tools/bd_k16_candidate.py --output ${SCRATCH}/tmp/new-k16 --case LC3
The destination must be new (no overwrite), outside tests/IEA15MW, including
symlink aliases. Inputs are official in-repository text, not an aeroelast mesh.
Gaussian quadrature and all-node TD/RD are requested; no fields are measured.
The existing static driver retains fixed root, identity root DCM, zero omega,
and scaled dead/global loads. LC1 includes gravity; LC2/LC3/LC4 do not.
SumPrint=True requests reference geometry: source-supported .sum.yaml arrays
Init_Nodes_E<n> contain initial XYZ and Wiener-Milenkovic (WM) parameters in
IEC coordinates, per element/local FE/GLL node, without global indexing.
No instrumentation is needed merely to request these arrays; none are acquired
or validated here. ES18.5 serialization has finite precision; its independently
propagated bound is missing. G1 is not closed: undeformed motion-node positions/
directors remain unacquired/unvalidated, source-content equivalence is unverified,
and G2-G4 evidence is missing. RD/WM conversion, node ordering, units and root/tip
coverage still require verification.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
from bd_static_reference import case_loads, driver_text  # noqa: E402
from of_bd_blade_modal import replace  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / 'tests/IEA15MW'
GEOMETRY = REFERENCES / ('validation papers/github-IEA-15-240-RWT/OpenFAST/'
                         'IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat')
PROPERTIES = REFERENCES / 'reference/IEA-15-240-RWT_BeamDyn_blade.dat'
MISSING = [
    'undeformed motion-node positions and directors',
    'verified source-content equivalence to installed binary',
    'verified node ordering, root/tip coverage, units and rotation conversion',
    'G2 geometric/load equivalence', 'G3 independent controls and uncertainty bounds',
    'G4 shared instrument',
]


def generate(target, case, order, scale):
    """Render first, then exclusively create the isolated candidate directory."""
    case_loads(case, scale)
    if order < 10:
        raise ValueError('order must be >=10')
    if target.is_relative_to(REFERENCES.resolve()):
        raise ValueError('output must be outside tests/IEA15MW reference tree')
    blade = replace(PROPERTIES.read_text(), 'damp_type', 0)
    blade = blade.replace(' ---------------------- DISTRIBUTED PROPERTIES',
                          ' ---------------------- MODAL DAMPING ---------------------\n'
                          '1 n_modes - unused\n0 zeta - unused\n'
                          ' ---------------------- DISTRIBUTED PROPERTIES')
    beam = re.sub(r'---------------------- PITCH ACTUATOR PARAMETERS.*?'
                  r'(?=---------------------- OUTPUTS)', '', GEOMETRY.read_text(), flags=re.S)
    for key, value in {'order_elem': order, 'quadrature': 1, 'NRMax': 50,
                       'SumPrint': 'True', 'BldFile': f'"{target / "blade.dat"}"',
                       'OutFmt': '"ES20.10E3"', 'BldNd_BlOutNd': '"All"'}.items():
        beam = replace(beam, key, value)
    # Retain the official all-node schema; add only the missing declared component.
    if beam.count('"TDyr"\n') != 1:
        raise ValueError('expected one all-node TDyr declaration')
    beam = beam.replace('"TDyr"\n', '"TDyr"\n"TDzr"\n')
    metadata = {'case': case, 'order': order, 'scale': scale,
                'geometry': str(GEOMETRY.relative_to(ROOT)),
                'properties': str(PROPERTIES.relative_to(ROOT)),
                'g1_closed': False, 'executed': False, 'missing_evidence': MISSING,
                'reference_geometry_contract': {
                    'request': {'SumPrint': True},
                    'status': 'source-supported; not acquired or validated',
                    'source_evidence': 'parent source review: BeamDyn_IO.f90:2049-2060',
                    'summary_suffix': '.sum.yaml',
                    'array_schema': 'Init_Nodes_E<n>',
                    # Semantic descriptions, not invented literal YAML column labels.
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
                }}
    files = {'blade.dat': blade, 'beam.dat': beam,
             'driver.dvr': driver_text(case, scale, target / 'beam.dat'),
             'candidate.json': json.dumps(metadata, indent=2) + '\n'}
    try:
        target.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise ValueError(f'output already exists: {target}') from exc
    for name, text in files.items():
        with (target / name).open('x') as stream:
            stream.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', type=Path, required=True,
                        help='new isolated directory outside tests/IEA15MW (no overwrite)')
    parser.add_argument('--case', choices=('LC1', 'LC2', 'LC3', 'LC4'), default='LC3')
    parser.add_argument('--order', type=int, default=60, help='polynomial order >=10 (default: 60)')
    parser.add_argument('--scale', type=float, default=.0001,
                        help='positive finite load scale (default: .0001)')
    args = parser.parse_args()
    target = args.output.resolve()
    try:
        generate(target, args.case, args.order, args.scale)
    except (ValueError, OSError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
    print(target)
    return 0


if __name__ == '__main__':
    sys.exit(main())
