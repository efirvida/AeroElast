"""L4 spanwise rows: is the AeroDyn 'dt dependence' an azimuth-sampling artefact?

Prediction (written before the run): per-node alpha/Cn taken from the LAST ROW (current parser) vary with dt by several pp of
max|alpha| because the last row sits at a dt-dependent azimuth and the 6 deg shaft tilt makes alpha azimuth-dependent; the
revolution MEAN of the same channels does not move with dt (<0.1 deg) and its distance to the BEM is the physical residual.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
import tests.support.openfast_bem as ob  # noqa: E402
import tests.validation.bem.test_bem_openfast_parity as T  # noqa: E402
from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402

driver = T._openfast_driver_or_skip()
blade = ob.parse_aerodyn_blade(T.BLADE)
r_ref = ob.IEA15MW_HUB_RAD + blade.bl_spn
interior = (r_ref > T.INTERIOR[0] * r_ref[-1]) & (r_ref < T.INTERIOR[1] * r_ref[-1])


def read_all(path):
    lines = Path(path).read_text(errors="replace").splitlines()
    h = next(i for i, ln in enumerate(lines) if ln.startswith(" Time") or "\tTime" in ln)
    header = lines[h].split()
    rows = []
    for ln in lines[h + 2:]:
        f = ln.split()
        if len(f) == len(header):
            try:
                rows.append([float(x) for x in f])
            except ValueError:
                pass
    return header, np.array(rows)


with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    primary = ob.write_bem_primary(T.PRIMARY, td / "AD_clean.dat")
    blade_aero = ob.build_blade_aero_from_aerodyn(primary)
    solver = BEMSolver(blade_aero, rho=1.225, mu=1.464e-5, precone=ob.IEA15MW_PRECONE_DEG,
                       tilt=ob.IEA15MW_SHAFT_TILT_DEG, yaw=0.0, hub_height=ob.IEA15MW_HUB_HEIGHT, shear_exp=0.0)
    comp = [solver.compute(*c) for c in T.CASES]
    print("V   dt     last-row max|da| deg  mean|dCn| %   |  rev-mean max|da| deg  mean|dCn| %   | az(last) deg  alpha spread over rev (interior, max) deg")
    for dt in (1.0, 0.5, 0.25):
        sub = td / f"dt{dt}"
        sub.mkdir()
        dvr = ob.write_aerodyn_dvr(sub / "a1.dvr", primary, T.CASES, dt=dt, tmax=60.0)
        outs = ob.run_aerodyn_driver(driver, dvr, sub)
        for k, case in enumerate(T.CASES):
            header, rows = read_all(outs[k])
            col = {n: j for j, n in enumerate(header)}
            nn = sum(1 for n in header if n.startswith("AB1N") and n.endswith("Alpha"))
            az = rows[:, col["Azimuth"]]
            wraps = np.where(np.diff(az) < 0)[0]
            seg = rows[wraps[-2] + 1:] if len(wraps) >= 2 else rows
            al = lambda r: np.array([r[col[f"AB1N{i:03d}Alpha"]] for i in range(1, nn + 1)])  # noqa: E731
            cn = lambda r: np.array([r[col[f"AB1N{i:03d}Cn"]] for i in range(1, nn + 1)])  # noqa: E731
            res = {}
            for tag, a, c in (("last", al(rows[-1]), cn(rows[-1])),
                              ("mean", np.mean([al(r) for r in seg], axis=0), np.mean([cn(r) for r in seg], axis=0))):
                da = np.abs(comp[k].alpha - a)[interior].max()
                dc = np.nanmean(np.abs((comp[k].cn - c) / np.where(c == 0, np.nan, c))[interior])
                res[tag] = (da, 100 * dc)
            spread = (np.array([al(r) for r in seg]).max(axis=0) - np.array([al(r) for r in seg]).min(axis=0))[interior].max()
            print(f"{case[0]:4.1f} {dt:5.3f}  {res['last'][0]:8.3f} {res['last'][1]:8.2f}      |  {res['mean'][0]:8.3f} {res['mean'][1]:8.2f}      | {az[-1]:7.1f}  {spread:6.2f}")
