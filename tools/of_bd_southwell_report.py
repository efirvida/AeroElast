"""Record per-code vs analytic errors for completed Southwell campaigns.

Usage: OPENBLAS_NUM_THREADS=1 python tools/of_bd_southwell_report.py
Reads only modes.json, generated parameters, and the published analytic box;
writes logs/check/of_bd_sw/report.json and report.md. No eigensolves or MBC3.
Errors are signed 100*(computed/reference-1); no cross-code judging.
Frequency tolerance: 1% parked, 5% rotating (S8 ceiling); K tolerance: 5%.
Box parked reference is the published Timoshenko value (EB also recorded);
box rotating reference is the published rotating EB value, since the file has
no rotating Timoshenko reference. This model difference remains explicit.
"""
from pathlib import Path
import json
import math

from of_bd_southwell import reference

ROOT = Path(__file__).resolve().parents[1]
LOGS = ROOT / 'logs/check/of_bd_sw'


def first_modes(data):
    """Identify directions through the public MBC3 dominant-state descriptions."""
    result = {}
    for mode in data['modes']:
        state = mode['dominant'][0]['state']
        if 'displacement in X,' in state or '1st flapwise' in state:
            result.setdefault('flap', mode)
        if 'displacement in Y,' in state or '1st edgewise' in state:
            result.setdefault('edge', mode)
    if set(result) != {'flap', 'edge'}:
        raise ValueError(f'Cannot classify first flap/edge: {data["case"]}')
    return result


def main():
    ref = reference()
    records = []
    for parked in ['sw_bd10_0', 'sw_bd20_0', 'sw_ed_0', 'sw_boxbd10d_0',
                   'sw_boxbd15d_0', 'sw_boxbd20d_0', 'sw_boxed_0']:
        box = parked.startswith('sw_box')
        rotating = parked[:-1] + ('rot' if box else '10')
        if not box and 'bd' in rotating:
            rotating = rotating.replace('_10', 'd_10')
        data0 = json.loads((LOGS / parked / 'modes.json').read_text())
        data1 = json.loads((LOGS / rotating / 'modes.json').read_text())
        modes0, modes1 = first_modes(data0), first_modes(data1)
        source = None
        if box:
            # Frozen parameters include the exact published analytic references.
            source = json.loads((ROOT / 'tests/IEA15MW/of_bd_modal' / parked /
                                 'parameters.json').read_text())['published_box']
        omega = source['Omega_rad_s'] if box else 10*math.pi/30
        for i, direction in enumerate(['flap', 'edge']):
            f0, f1 = modes0[direction]['hz'], modes1[direction]['hz']
            k = ((2*math.pi*f1)**2-(2*math.pi*f0)**2)/omega**2
            if box:
                analytic = source['analytic'][direction][0]
                f0ref = analytic['omega_0_Timoshenko_rad_s']/(2*math.pi)
                f0eb = analytic['omega_0_EB_rad_s']/(2*math.pi)
                f1ref = analytic['omega_rotating_rad_s']/(2*math.pi)
                kref = analytic['K1'] - (direction == 'edge')
            else:
                f0ref = f0eb = ref['parked_hz'][i]
                kref = ref[f'K_{direction}']
                f1ref = math.sqrt(f0ref**2+kref*(omega/(2*math.pi))**2)
            for data, frequency, frequency_ref, rpm in [
                    (data0, f0, f0ref, 0), (data1, f1, f1ref, omega*30/math.pi)]:
                records.append({'case': data['case'], 'direction': direction,
                                'rpm': rpm, 'hz': frequency, 'analytic_hz': frequency_ref,
                                'frequency_error_percent': 100*(frequency/frequency_ref-1),
                                'analytic_parked_EB_hz': f0eb,
                                'frequency_tolerance_percent': 5 if rpm else 1,
                                'K': k if rpm else None, 'analytic_K': kref if rpm else None,
                                'K_error_percent': 100*(k/kref-1) if rpm else None,
                                'K_tolerance_percent': 5,
                                'collective_offblock_relative': data['collective_offblock_relative']})
    for parked in ['sw_boxbd10d_0', 'sw_boxbd15d_0', 'sw_boxbd20d_0']:
        data = json.loads((LOGS / parked / 'modes.json').read_text())
        published = json.loads((ROOT / 'tests/IEA15MW/of_bd_modal' / parked /
                                'parameters.json').read_text())['published_box']
        for direction, marker, index in [('flap2', 'translational displacement in X,', 1),
                                          ('edge2', 'translational displacement in Y,', 1),
                                          ('torsion1', 'rotational displacement in Z,', 0)]:
            modes = [m for m in data['modes'] if marker in m['dominant'][0]['state']]
            mode = modes[index]
            if direction == 'torsion1':
                frequency_ref = published['analytic']['torsion']['omega_0_rad_s']/(2*math.pi)
                eb = frequency_ref
            else:
                analytic = published['analytic'][direction[:-1]][1]
                frequency_ref = analytic['omega_0_Timoshenko_rad_s']/(2*math.pi)
                eb = analytic['omega_0_EB_rad_s']/(2*math.pi)
            records.append({'case': parked, 'direction': direction, 'rpm': 0,
                            'hz': mode['hz'], 'analytic_hz': frequency_ref,
                            'frequency_error_percent': 100*(mode['hz']/frequency_ref-1),
                            'analytic_parked_EB_hz': eb, 'frequency_tolerance_percent': 1,
                            'K': None, 'analytic_K': None, 'K_error_percent': None,
                            'K_tolerance_percent': 5,
                            'collective_offblock_relative': data['collective_offblock_relative']})
    iea = {}
    for name in ['bp20_0', 'bp20_7.56']:
        data = json.loads((ROOT / 'logs/check/of_bd' / name / 'modes.json').read_text())
        iea[name] = first_modes(data)['edge']
    rotating_bd = [r for r in records if r['rpm'] and 'bd' in r['case']]
    artifact = all(abs(r['K_error_percent']) <= 5 for r in rotating_bd if r['direction'] == 'flap')
    artifact &= all(abs(r['K_error_percent']) > 5 for r in rotating_bd if r['direction'] == 'edge')
    report = {'references': ref, 'published_box_reference': source,
              'IEA_1E_dominant_states': iea, 'results': records,
              'verdict': ('BD rotating 1E is an artifact; only parked 1E usable' if artifact
                          else 'Canonical criterion not established; inspect raw errors'),
              'qualification': 'Canonical evidence identifies missing in-plane softening in the '
              'extracted collective linearized dynamics, not its precise OpenFAST source-code cause. '
              'Parked usable means unaffected by this artifact, not blanket validation of IEA physics.'}
    (LOGS / 'report.json').write_text(json.dumps(report, indent=2)+'\n')
    lines = ['# Southwell diagnostic: each code vs independent analytic reference', '',
             '| Case | mode | rpm | f Hz | analytic Hz | error % | K | analytic K | K error % |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in records:
        kcols = ('— | — | —' if r['K'] is None else
                 f'{r["K"]:.6f} | {r["analytic_K"]:.6f} | {r["K_error_percent"]:+.3f}')
        lines.append(f'| {r["case"]} | {r["direction"]} | {r["rpm"]:.6f} | '
                     f'{r["hz"]:.9f} | {r["analytic_hz"]:.9f} | '
                     f'{r["frequency_error_percent"]:+.4f} | {kcols} |')
    lines += ['', report['verdict'], '', report['qualification'], '', __doc__]
    (LOGS / 'report.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
