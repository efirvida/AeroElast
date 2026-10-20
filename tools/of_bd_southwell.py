"""Uniform rotating cantilever diagnostic, built by of_bd_blade_modal.py.

Usage: python tools/of_bd_blade_modal.py build --case sw_bd20_0 (also _10,
 sw_bd10_0/_10 and sw_ed_0/_10); submit/post use the same interface.
Published box cases: sw_boxbd10d_0/_rot, sw_boxbd15d_0/_rot,
sw_boxbd20d_0/_rot and sw_boxed_0/_rot read the authoritative SI parameters
from logs/check/analytic_beam/cases.json. BeamDyn local X is flap, local Y edge;
the independent source's global axes differ, but EI/m/shear/rotary values map
by the named bending direction. The ED fit and its residual are saved per case.
If undamped integration fails, sw_bd10d_10 and sw_bd20d_10 preserve those
attempts and use rhoinf=0.8, DT=0.001 (numerical damping only).
Outputs are isolated under logs/check/of_bd_sw. Postprocess on a compute node
with OPENBLAS_NUM_THREADS=1 /scratch/leahk/eduardo.donestevez/venv/bin/python.
Run this file directly to record the independent analytic reference.

Reference: Euler-Bernoulli rotating-cantilever energy equation, with centrifugal
 tension T(x)=m Omega^2 (L^2-x^2)/2. Its first-order eigenvalue perturbation is
 K_flap=int_0^1 (1-s^2)/2 phi'(s)^2 ds / int_0^1 phi(s)^2 ds;
 K_edge=K_flap-1 from the in-plane centrifugal potential -m Omega^2 u^2/2.
phi is the exact nonrotating clamped-free mode (cosh beta cos beta = -1).
This explicitly derived independent reference (not a fit to OpenFAST) gives
K_flap ~1.193 and K_edge ~0.193. A 5% relative tolerance on K is fixed before
execution, allowing finite-speed shape changes and numerical discretization.
Parked frequencies must also agree with beta^2 sqrt(EI/m/L^4)/(2pi) within 1%.
The rectangular/isotropic-principal-axis section is diagonal, with unequal
bending rigidities to separate flap/edge; no offsets, twist, damping or gravity.
"""
from pathlib import Path
import json
import math

LENGTH = 10.0
MASS = 10.0
FLAP_EI = 1e5
EDGE_EI = 2.25e5


def reference():
    import numpy as np
    from scipy.integrate import quad

    beta = 1.875104068711961
    sigma = (np.cosh(beta) + np.cos(beta)) / (np.sinh(beta) + np.sin(beta))

    def phi(s):
        z = beta*s
        return np.cosh(z)-np.cos(z)-sigma*(np.sinh(z)-np.sin(z))

    def derivative(s):
        z = beta*s
        return beta*(np.sinh(z)+np.sin(z)-sigma*(np.cosh(z)-np.cos(z)))

    k = quad(lambda s: (1-s*s)/2*derivative(s)**2, 0, 1, epsabs=1e-12)[0]
    k /= quad(lambda s: phi(s)**2, 0, 1, epsabs=1e-12)[0]
    return {'K_flap': k, 'K_edge': k-1, 'K_relative_tolerance': 0.05,
            'parked_relative_tolerance': 0.01,
            'parked_hz': [beta**2*math.sqrt(ei/MASS/LENGTH**4)/(2*math.pi)
                          for ei in [FLAP_EI, EDGE_EI]],
            'source': __doc__}


def build(api, selected):
    """Reuse the campaign parameter editor; never write its shared blade files."""
    import numpy as np

    template = api.INPUTS / 'bp20_0'
    for name, kind, rpm, order in api.CASES:
        if not name.startswith('sw_') or (selected and name != selected):
            continue
        box = None
        if kind.startswith('sw_box'):
            box = json.loads((api.ROOT / 'logs/check/analytic_beam/cases.json').read_text())['uniform_box']
        length = box['geometry']['L'] if box else LENGTH
        mass_per_length = box['beam']['m'] if box else MASS
        flap_ei = box['beam']['EI_flap'] if box else FLAP_EI
        edge_ei = box['beam']['EI_edge'] if box else EDGE_EI
        if box and rpm:
            rpm = box['Omega_rad_s']*30/math.pi
        case = api.INPUTS / name
        out = api.LOGS / name
        case.mkdir(parents=True, exist_ok=True)
        out.mkdir(parents=True, exist_ok=True)
        ed = (template / 'ed.dat').read_text()
        params = {'RotSpeed': rpm, 'TipRad': length, 'HubRad': 0,
                  'ShftTilt': 0, 'OverHang': 0, 'BldNodes': 100,
                  'FlapDOF2': 'False'}
        for i in range(1, 4):
            params[f'PreCone({i})'] = 0
            params[f'BldFile{i}'] = f'"{case / "ed_blade.dat"}"'
        for key, value in params.items():
            ed = api.replace(ed, key, value)
        (case / 'ed.dat').write_text(ed)
        (case / 'parameters.json').write_text(json.dumps(
            {'L': length, 'm': mass_per_length, 'EI_flap': flap_ei,
             'EI_edge': edge_ei, 'rpm': rpm, 'published_box': box}, indent=2)+'\n')
        # Fit the exact first cantilever shape to ElastoDyn's x^2..x^6 basis.
        s = np.linspace(0, 1, 1001)
        b = 1.875104068711961
        c = (np.cosh(b)+np.cos(b))/(np.sinh(b)+np.sin(b))
        phi = np.cosh(b*s)-np.cos(b*s)-c*(np.sinh(b*s)-np.sin(b*s))
        coefficients = np.linalg.lstsq(np.array([s**i for i in range(2, 7)]).T,
                                       phi/phi[-1], rcond=None)[0]
        coefficients /= coefficients.sum()
        residual = np.array([s**i for i in range(2, 7)]).T @ coefficients - phi/phi[-1]
        (case / 'mode_fit.json').write_text(json.dumps(
            {'basis': 's^2..s^6; least squares, 1001 equally spaced points; tip normalized',
             'coefficients': coefficients.tolist(),
             'max_absolute_residual': float(np.max(np.abs(residual))),
             'rms_residual': float(np.sqrt(np.mean(residual**2)))}, indent=2)+'\n')
        blade = (api.INPUTS / 'ed_blade.dat').read_text()
        for key, value in {'NBlInpSt': 2, 'BldFlDmp1': 0, 'BldFlDmp2': 0,
                           'BldEdDmp1': 0}.items():
            blade = api.replace(blade, key, value)
        start = blade.index(' 0.000000000000000e+00')
        end = blade.index('---------------------- BLADE MODE SHAPES')
        blade = blade[:start] + ''.join(
            f'{s} 0.5 0 {mass_per_length} {flap_ei} {edge_ei}\n' for s in [0, 1]) + blade[end:]
        for mode in ['BldFl1Sh', 'BldEdgSh']:
            for i, coefficient in enumerate(coefficients, 2):
                blade = api.replace(blade, f'{mode}({i})', f'{coefficient:.16e}')
        (case / 'ed_blade.dat').write_text(blade)
        bd = (template / 'bd.dat').read_text()
        bd = api.replace(bd, 'order_elem', order or 10)
        if kind.endswith('d'):
            bd = api.replace(bd, 'rhoinf', 0.8)
        bd = api.replace(bd, 'kp_total', 3)
        bd = api.replace(bd, 'BldFile', f'"{case / "bd_blade.dat"}"')
        start = bd.index('     1     50')
        end = bd.index('---------------------- MESH')
        bd = bd[:start] + ('1 3\nkp_xr kp_yr kp_zr initial_twist\n(m) (m) (m) (deg)\n'
                           f'0 0 0 0\n0 0 {length/2} 0\n0 0 {length} 0\n') + bd[end:]
        (case / 'bd.dat').write_text(bd)
        blade = (api.INPUTS / 'bd_blade.dat').read_text()
        blade = api.replace(blade, 'station_total', 2)
        blade = blade[:blade.index('DISTRIBUTED PROPERTIES')].rsplit('\n', 1)[0]
        blade += '\n ---------------------- DISTRIBUTED PROPERTIES ----------------\n'
        if box:
            beam = box['beam']
            stiffness = np.diag([beam['kGA_flap'], beam['kGA_edge'], beam['EA'],
                                 edge_ei, flap_ei, beam['GJ']])
            mass = np.diag([mass_per_length]*3 + [beam['rotary_edge'],
                           beam['rotary_flap'], beam['rotary_polar']])
        else:
            stiffness = np.diag([1e10, 1e10, 1e10, edge_ei, flap_ei, 1e6])
            mass = np.diag([mass_per_length]*3 + [0.001, 0.001, 0.002])
        for station in [0, 1]:
            blade += f'{station}\n'
            for matrix in [stiffness, mass]:
                blade += '\n'.join(' '.join(f'{v:.16e}' for v in row) for row in matrix)+'\n\n'
        (case / 'bd_blade.dat').write_text(blade)
        fst = (template / 'case.fst').read_text()
        times = [1.0] if rpm == 0 else [30+k*60/rpm/12 for k in range(12)]
        params = {'EDFile': f'"{case / "ed.dat"}"', 'CompElast': 1 if kind.endswith('ed') else 2,
                  'DT': 0.001 if kind.endswith('d') else 0.005,
                  'TMax': times[-1]+0.01, 'TStart': times[-1]+0.01,
                  'NLinTimes': len(times), 'LinTimes': ', '.join(f'{t:.6f}' for t in times)}
        for i in range(1, 4):
            params[f'BDBldFile({i})'] = f'"{case / "bd.dat"}"'
        for key, value in params.items():
            fst = api.replace(fst, key, value)
        (case / 'case.fst').write_text(fst)
        (out / 'case.fst').write_text(fst)
        print(name, case, out)


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1] / 'logs/check/of_bd_sw'
    root.mkdir(parents=True, exist_ok=True)
    result = reference()
    (root / 'reference.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
