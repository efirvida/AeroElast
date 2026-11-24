"""L4 T3 diagnostic: polar-handling variants (spline s=0.1 | exact spline | linear) at the pitched V-01 points,
plus Cl(spline)-Cl(table) at the AD alpha per station, plus tip-node contribution."""

import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t5_options as t5  # noqa: E402
import l4_t3_span as sp  # noqa: E402  (runs nothing at import: guarded below)

ob, p = t5.ob, t5.p


def exact_spline(solver, blade):
    by_id = {a.name: a for a in blade.airfoils}
    for af, st in zip(solver.rotor.af, blade.stations, strict=False):
        pol = by_id[st.airfoil.name].polars[0]
        a = np.asarray(pol.alpha)
        for nm, y in (("cl_spline", pol.cl), ("cd_spline", pol.cd)):
            Y = np.c_[y, y]
            setattr(af, nm, RectBivariateSpline(a, [1e1, 1e15], Y, kx=3, ky=1, s=0))


def run(case, variant, blade):
    s = p.BEMSolver(
        blade,
        rho=t5.RHO,
        mu=t5.MU_AD,
        precone=sp.CONE,
        tilt=sp.TILT,
        yaw=0.0,
        hub_height=ob.IEA15MW_HUB_HEIGHT,
        shear_exp=0.0,
    )
    if variant == "linear":
        t5._patch_linear(s, blade)
    elif variant == "exact":
        exact_spline(s, blade)
    r = s.compute(*case)
    return r.thrust, r.torque, s


for case in sp.CASES:
    blade, ref = t5.ad_run(case, tilt=sp.TILT, precone=sp.CONE)
    line = f"{case}: AD T={ref.thrust:.0f} Q={ref.torque:.0f} |"
    for v in ("spline", "exact", "linear"):
        T, Q, s = run(case, v, blade)
        line += f" {v}: T{T / ref.thrust - 1:+.2%} Q{Q / ref.torque - 1:+.2%}"
    print(line)
    # Cl error of the smoothed spline vs the table at the operating alpha, outboard stations
    s = run(case, "spline", blade)[2]
    by_id = {a.name: a for a in blade.airfoils}
    errs = []
    for i, st in enumerate(blade.stations):
        pol = by_id[st.airfoil.name].polars[0]
        al = np.deg2rad(ref.alpha_deg[i])
        cl_s = s.rotor.af[i].evaluate(al, 1e6)[0]
        cl_t = np.interp(al, pol.alpha, pol.cl)
        errs.append(cl_s - cl_t)
    errs = np.array(errs)
    print(
        "   Cl(spline)-Cl(table) at AD alpha, i=10..48 step 6:",
        np.round(errs[10:49:6], 4),
        "| mean rel over r>40:",
        f"{np.mean(errs[15:48]):+.4f}",
    )
