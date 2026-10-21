"""BeamDyn v5 rotating free vibration (K10, method (a)), without importing aeroelast.

python tools/bd_rotating_freevib.py build --beam uniform --rpm 10
python tools/bd_rotating_freevib.py run   --beam uniform --rpm 10   # on a compute node
python tools/bd_rotating_freevib.py post  --beam uniform --rpm 10
python tools/bd_rotating_freevib.py gate                             # K10.2 table

Why: the MBC3-extracted BeamDyn linearization has the in-plane centrifugal softening
missing (rotating edge Southwell K 1.19 vs 0.193, logs/check/of_bd_sw/report.md), so
it has no admissible rated 1E. This takes the frequency from the nonlinear time
integration instead: driver `DynamicSolve=True`, root turning at constant omega about
BeamDyn x (flap; shell +y), gravity 0, no structural damping, generalized-alpha
rhoinf 0.8. QuasiStaticInit is off: the beam starts undeformed with the rigid-body
velocity of the root, and the step tip load (global x flap, global y edge) plus the
centrifugal step excite the modes. The v5 driver applies loads unrotated in the
global frame, so the edge load turns at exactly 1P in the rotating frame; the fit
carries that line as a known harmonic.

Frequencies: tip TDxr/TDyr (root frame) after `--skip` s; joint nonlinear least
squares of a constant, the 1P harmonic and damped sinusoids shared by both
channels, modes added greedily from Hann-window spectrum peaks of the residual.
Beams: `uniform` (tools/of_bd_southwell.py constants), `box` (published uniform
box, logs/check/analytic_beam/cases.json), `iea` (official v1.1 deck, single
member, the deck's 50 keypoints, as the S-4 BeamDyn p20 arbiter). `--hub` puts the
root at z=hub on the rotation axis line (the shell's rated modal has hub 0).
Inputs under tests/IEA15MW/bd_freevib/<name>, outputs under logs/check/bd_freevib.
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

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from of_bd_blade_modal import replace  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE = Path('/scratch/leahk/eduardo.donestevez')
GEOMETRY = BASE / 'ofruns/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat'
PROPERTIES = ROOT / 'tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat'
BOX = ROOT / 'logs/check/analytic_beam/cases.json'
BINARY = BASE / 'tools/openfast-env/bin/beamdyn_driver'
INPUTS = ROOT / 'tests/IEA15MW/bd_freevib'
LOGS = ROOT / 'logs/check/bd_freevib'
CHANNELS = ('TipTDxr', 'TipTDyr')
# beam: default order, tip load (N), t_final (s), dt (s), fit band (Hz)
DEFAULTS = {'uniform': (10, 1.0, 300.0, 0.005, (0.1, 3.0)),
            'box': (10, 1e3, 150.0, 0.002, (0.2, 6.0)),
            'iea': (20, 1e3, 300.0, 0.005, (0.1, 3.0))}


# ----------------------------------------------------------------------------- fit
def _basis(t, one_p_hz, freqs, sigmas):
    w = 2 * np.pi * one_p_hz
    cols = [np.ones_like(t), np.cos(w * t), np.sin(w * t)]
    for f, s in zip(freqs, sigmas, strict=True):
        decay = np.exp(-s * t)
        cols += [decay * np.cos(2 * np.pi * f * t), decay * np.sin(2 * np.pi * f * t)]
    return np.column_stack(cols)


def _solve(t, y, one_p_hz, freqs, sigmas):
    b = _basis(t, one_p_hz, freqs, sigmas)
    coef = np.linalg.lstsq(b, y.T, rcond=None)[0]
    return coef, y.T - b @ coef


def _peak(t, residual, one_p_hz, band, taken):
    """Strongest combined-spectrum peak in band, away from 1P and from fitted modes."""
    n = len(t)
    window = np.hanning(n)
    nfft = 8 * (1 << (n - 1).bit_length())
    power = sum(np.abs(np.fft.rfft(residual[:, c] * window, nfft)) ** 2 for c in range(residual.shape[1]))
    f = np.fft.rfftfreq(nfft, t[1] - t[0])
    guard = 3.0 / (t[-1] - t[0])
    mask = (f > band[0]) & (f < band[1]) & (np.abs(f - one_p_hz) > guard)
    for g in taken:
        mask &= np.abs(f - g) > guard
    if not mask.any():
        return None
    i = np.flatnonzero(mask)[np.argmax(power[mask])]
    if 0 < i < len(f) - 1:  # parabolic refinement on the log spectrum
        a, b, c = np.log(power[i - 1:i + 2] + 1e-300)
        return f[i] + 0.5 * (a - c) / (a - 2 * b + c) * (f[1] - f[0])
    return f[i]


def band_limit(t, y, band):
    """Zero-phase low-pass (Butterworth 8, cutoff 1.5 band[1]) and decimation to >= 6 band[1].

    A damped sinusoid is an eigenfunction of the filter in both passes, so its frequency and
    decay survive (only its amplitude scales); the filter end transients are dropped.
    """
    from scipy.signal import butter, sosfiltfilt

    dt = t[1] - t[0]
    q = int(1 / (dt * 6 * band[1]))
    if q <= 1:
        return t, y
    cutoff = 1.5 * band[1]
    y = sosfiltfilt(butter(8, cutoff, fs=1 / dt, output='sos'), y, axis=1)
    margin = int(20 / (cutoff * dt))
    keep = slice(margin, len(t) - margin, q)
    return t[keep], y[:, keep]


def fit_frequencies(t, y, one_p_hz, band, max_modes=6, rel_amplitude=1e-3):
    """Modes shared by the channels of y (n_channels x n_samples); see the module docstring."""
    from scipy.optimize import least_squares

    t, y = band_limit(np.asarray(t, float), np.asarray(y, float), band)
    t = t - t[0]
    freqs, sigmas = [], []
    largest = 0.0
    scale = np.linalg.norm(y - y.mean(axis=1, keepdims=True))
    guard = 3.0 / t[-1]
    for _ in range(max_modes):
        _, residual = _solve(t, y, one_p_hz, freqs, sigmas)
        guess = _peak(t, residual, one_p_hz, band, freqs)
        if guess is None:
            break
        trial_f, trial_s = freqs + [guess], sigmas + [0.0]
        k = len(trial_f)
        # Bounded to the band: unbounded, the box parked fit wandered to -20.4 Hz (11621046).
        lower = [band[0]] * k + [-1.0] * k
        upper = [band[1]] * k + [1.0] * k
        fit = least_squares(
            lambda p: _solve(t, y, one_p_hz, p[:k], p[k:])[1].ravel() / scale,
            np.clip(np.array(trial_f + trial_s), lower, upper), bounds=(lower, upper),
            xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
        f = np.sort(fit.x[:k])
        if f[0] - band[0] < guard or band[1] - f[-1] < guard or (k > 1 and np.diff(f).min() < guard):
            break  # a mode pinned to the band edge or split in two: the trial mode is not real
        coef, _ = _solve(t, y, one_p_hz, fit.x[:k], fit.x[k:])
        amplitude = np.hypot(coef[3 + 2 * (k - 1)], coef[4 + 2 * (k - 1)]).max()
        largest = max(largest, amplitude)
        if amplitude < rel_amplitude * largest:
            break
        freqs, sigmas = list(fit.x[:k]), list(fit.x[k:])
    coef, residual = _solve(t, y, one_p_hz, freqs, sigmas)
    modes = [{'f_hz': float(f), 'sigma_1_s': float(s),
              'amplitudes': np.hypot(coef[3 + 2 * i], coef[4 + 2 * i]).tolist()}
             for i, (f, s) in enumerate(zip(freqs, sigmas, strict=True))]
    modes.sort(key=lambda m: m['f_hz'])
    return {'modes': modes, 'relative_residual': float(np.linalg.norm(residual) / scale)}


# --------------------------------------------------------------------------- inputs
def case_name(beam, rpm, order, hub, dt, load, profile=None):
    tag = f'-eie{Path(profile).stem}' if profile else ''
    return f'{beam}-p{order}-r{rpm:g}-h{hub:g}-dt{dt:g}-F{load:g}{tag}'


def _blade_properties(stations=None):
    """Official property file, undamped; `stations` replaces its distributed table."""
    prop = replace(PROPERTIES.read_text(), 'damp_type', 0)
    prop = prop.replace(' ---------------------- DISTRIBUTED PROPERTIES',
                        ' ---------------------- MODAL DAMPING ---------------------\n'
                        '1 n_modes - unused\n0 zeta - unused\n'
                        ' ---------------------- DISTRIBUTED PROPERTIES')
    if stations is None:
        return prop
    prop = replace(prop, 'station_total', len(stations))
    prop = prop[:prop.index(' ---------------------- DISTRIBUTED PROPERTIES')]
    prop += ' ---------------------- DISTRIBUTED PROPERTIES ----------------\n'
    for eta, stiffness, mass in stations:
        prop += f'{eta}\n'
        for matrix in (stiffness, mass):
            prop += '\n'.join(' '.join(f'{v:.16e}' for v in row) for row in matrix) + '\n\n'
    return prop


def _uniform_section(beam):
    """(length, stations) of the T7 uniform and box beams (BD X flap, Y edge)."""
    if beam == 'uniform':
        from of_bd_southwell import EDGE_EI, FLAP_EI, LENGTH, MASS
        stiffness = np.diag([1e10, 1e10, 1e10, EDGE_EI, FLAP_EI, 1e6])
        mass = np.diag([MASS] * 3 + [0.001, 0.001, 0.002])
        return LENGTH, [(0, stiffness, mass), (1, stiffness, mass)]
    box = json.loads(BOX.read_text())['uniform_box']
    b = box['beam']
    stiffness = np.diag([b['kGA_flap'], b['kGA_edge'], b['EA'], b['EI_edge'], b['EI_flap'], b['GJ']])
    mass = np.diag([b['m']] * 3 + [b['rotary_edge'], b['rotary_flap'], b['rotary_polar']])
    return box['geometry']['L'], [(0, stiffness, mass), (1, stiffness, mass)]


def rpm_of(beam, rpm):
    """The box rotates at its published Omega; any nonzero --rpm selects it."""
    if beam == 'box' and rpm:
        return json.loads(BOX.read_text())['uniform_box']['Omega_rad_s'] * 30 / math.pi
    return rpm


def _scaled_deck(profile):
    """Deck stations (eta, K, M) with K[3, 3] (edge EI) times the ratio profile `eta ratio` (B3-3 K13.C)."""
    from v05_deck_consistency import load_beamdyn
    table = np.loadtxt(profile)
    stations = []
    for eta, stiffness, mass in load_beamdyn(PROPERTIES):
        stiffness = stiffness.copy()
        stiffness[3, 3] *= np.interp(eta, table[:, 0], table[:, 1])
        stations.append((eta, stiffness, mass))
    return stations


def build(beam, rpm, order, hub, dt, load, t_final, profile=None):
    name = case_name(beam, rpm, order, hub, dt, load, profile)
    target = INPUTS / name
    target.mkdir(parents=True, exist_ok=False)
    primary = GEOMETRY.read_text()
    primary = re.sub(r'---------------------- PITCH ACTUATOR PARAMETERS.*?'
                     r'(?=---------------------- OUTPUTS)', '', primary, flags=re.S)
    if beam == 'iea':
        (target / 'blade.dat').write_text(_blade_properties(_scaled_deck(profile) if profile else None))
    else:
        length, stations = _uniform_section(beam)
        (target / 'blade.dat').write_text(_blade_properties(stations))
        primary = replace(primary, 'kp_total', 3)
        start = primary.index('     1     50')
        end = primary.index('---------------------- MESH')
        primary = primary[:start] + ('1 3\nkp_xr kp_yr kp_zr initial_twist\n(m) (m) (m) (deg)\n'
                                     f'0 0 0 0\n0 0 {length / 2} 0\n0 0 {length} 0\n') + primary[end:]
    for key, value in {'QuasiStaticInit': 'False', 'rhoinf': 0.8, 'order_elem': order,
                       'quadrature': 1, 'BldFile': f'"{target / "blade.dat"}"',
                       'OutFmt': '"ES20.10E3"'}.items():
        primary = replace(primary, key, value)
    primary = re.sub(r'(\s+OutList\s+[^\n]*\n).*', r'\1"TipTDxr, TipTDyr, TipTDzr"\nEND of input file\n',
                     primary, count=1, flags=re.S)
    (target / 'beam.dat').write_text(primary)
    omega = rpm_of(beam, rpm) * math.pi / 30
    lines = ['BeamDyn driver v5', 'K10 rotating free vibration', '--- simulation', 'True DynamicSolve',
             '0 t_initial', f'{t_final:g} t_final', f'{dt:g} dt', '--- gravity',
             '0 gravity_x', '0 gravity_y', '0 gravity_z', '--- frame', '0 x', '0 y', f'{hub:.12g} z',
             '! DCM rows', '! global to root', '1 0 0', '0 1 0', '0 0 1', 'True GlbRotBladeT0',
             '--- velocity', f'{omega:.15g} omega_x', '0 omega_y', '0 omega_z', '--- loads']
    for prefix, vector in (('Distr', (0, 0, 0)), ('Tip', (load, load, 0))):
        lines += [f'{v:.12g} {prefix}{c}' for c, v in
                  zip(('Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'), (*vector, 0, 0, 0), strict=True)]
    lines += ['0 NumPointLoads', '! point load table', '! units', '--- beam',
              f'"{target / "beam.dat"}" InputFile', '--- outputs', '0 WrVTK', '1 VTK_fps']
    (target / 'driver.dvr').write_text('\n'.join(lines) + '\n')
    return name


def read_output(path):
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        headers = line.split()
        if headers and headers[0] == 'Time' and all(c in headers for c in CHANNELS):
            rows = np.array([[float(v.replace('D', 'E')) for v in row.split()]
                             for row in lines[i + 2:] if row.strip()])
            if not np.isfinite(rows).all():
                raise ValueError('nonfinite BeamDyn output')
            return rows[:, 0], np.vstack([rows[:, headers.index(c)] for c in CHANNELS])
    raise ValueError('missing BeamDyn tip channels')


def post(name, beam, rpm, skip, band):
    output = LOGS / name
    if 'terminated normally' not in (output / 'driver.log').read_text():
        raise ValueError('driver did not terminate normally')
    t, y = read_output(output / 'driver.out')
    keep = t >= skip
    one_p = rpm_of(beam, rpm) / 60
    fitted = fit_frequencies(t[keep], y[:, keep], one_p_hz=one_p, band=band)
    for mode in fitted['modes']:
        mode['dominant'] = 'flap' if mode['amplitudes'][0] >= mode['amplitudes'][1] else 'edge'
    record = {'name': name, 'beam': beam, 'rpm': rpm_of(beam, rpm), 'one_p_hz': one_p,
              'skip_s': skip, 'band_hz': list(band), 'samples': int(keep.sum()),
              'channels': list(CHANNELS), **fitted, 'binary': str(BINARY)}
    (output / 'result.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


# ----------------------------------------------------------------------------- gate
def _first(record, family):
    return min((m for m in record['modes'] if m['dominant'] == family), key=lambda m: m['f_hz'])['f_hz']


def gate(names):
    """K10.2: parked within 1 % of the analytic, flap/edge K within 5 % (T7 references)."""
    reference = json.loads((ROOT / 'logs/check/of_bd_sw/reference.json').read_text())
    analytic = json.loads(BOX.read_text())['uniform_box']['analytic']
    refs = {'uniform': {'parked': dict(zip(('flap', 'edge'), reference['parked_hz'], strict=True)),
                        'K': {'flap': reference['K_flap'], 'edge': reference['K_edge']}},
            'box': {'parked': {f: analytic[f][0]['omega_0_Timoshenko_rad_s'] / (2 * math.pi) for f in ('flap', 'edge')},
                    'K': {'flap': analytic['flap'][0]['K1'], 'edge': analytic['edge'][0]['K1'] - 1}}}
    rows, ok = [], True
    for beam, (parked_name, rotating_name) in names.items():
        parked = json.loads((LOGS / parked_name / 'result.json').read_text())
        rotating = json.loads((LOGS / rotating_name / 'result.json').read_text())
        omega_hz = rotating['rpm'] / 60
        for family in ('flap', 'edge'):
            f0, fw = _first(parked, family), _first(rotating, family)
            k = (fw ** 2 - f0 ** 2) / omega_hz ** 2
            e0 = 100 * (f0 / refs[beam]['parked'][family] - 1)
            ek = 100 * (k / refs[beam]['K'][family] - 1)
            passed = abs(e0) <= 1 and abs(ek) <= 5
            ok &= passed
            rows.append(f'| {beam} | {family} | {f0:.9f} | {e0:+.4f} | {fw:.9f} | {k:.6f} | '
                        f'{refs[beam]["K"][family]:.6f} | {ek:+.3f} | {"PASS" if passed else "FAIL"} |')
    head = ('| beam | mode | f0 Hz | f0 err % | f_rot Hz | K | analytic K | K err % | K10.2 |\n'
            '|---|---|---:|---:|---:|---:|---:|---:|---|')
    return head + '\n' + '\n'.join(rows) + f'\n\nK10.2 {"PASS" if ok else "FAIL"}', ok


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', choices=('build', 'run', 'post', 'gate'))
    parser.add_argument('--beam', choices=tuple(DEFAULTS), default='uniform')
    parser.add_argument('--rpm', type=float, default=0.0, help='rotor speed (box: any nonzero = published Omega)')
    parser.add_argument('--order', type=int)
    parser.add_argument('--hub', type=float, default=0.0, help='root distance from the rotation axis (m)')
    parser.add_argument('--dt', type=float)
    parser.add_argument('--load', type=float, help='tip load per direction (N); halve it for the linearity check')
    parser.add_argument('--t-final', type=float)
    parser.add_argument('--skip', type=float, default=2.0, help='seconds discarded before the fit')
    parser.add_argument('--ei-edge-profile', help='iea only: `eta ratio` file scaling the deck K[3,3] (K13.C)')
    parser.add_argument('--gate', nargs=4, metavar=('UNI0', 'UNIW', 'BOX0', 'BOXW'), help='case names for gate')
    args = parser.parse_args()
    if args.action == 'gate':
        if not args.gate:
            parser.error('gate needs --gate UNI0 UNIW BOX0 BOXW')
        table, ok = gate({'uniform': args.gate[:2], 'box': args.gate[2:]})
        print(table)
        sys.exit(0 if ok else 1)
    order, load, t_final, dt, band = DEFAULTS[args.beam]
    order, load = args.order or order, args.load or load
    dt, t_final = args.dt or dt, args.t_final or t_final
    rpm = rpm_of(args.beam, args.rpm)
    name = case_name(args.beam, rpm, order, args.hub, dt, load, args.ei_edge_profile)
    if args.action == 'build':
        print(build(args.beam, rpm, order, args.hub, dt, load, t_final, args.ei_edge_profile))
    elif args.action == 'run':
        output = LOGS / name
        output.mkdir(parents=True, exist_ok=False)
        driver = output / 'driver.dvr'  # the driver writes its products beside the .dvr
        driver.write_text((INPUTS / name / 'driver.dvr').read_text())
        with (output / 'driver.log').open('w') as log:
            subprocess.run([str(BINARY), str(driver)], cwd=output, stdout=log, stderr=subprocess.STDOUT,
                           check=True, env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    else:
        record = post(name, args.beam, rpm, args.skip, band)
        print(json.dumps({k: record[k] for k in ('name', 'rpm', 'relative_residual', 'modes')}, indent=1))


if __name__ == '__main__':
    main()
