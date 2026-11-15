"""L4 T3 diagnostic: spanwise ours vs AeroDyn at the V-01 points (azimuth-averaged both sides)."""

import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t5_options as t5  # noqa: E402

ob, p = t5.ob, t5.p
TILT, CONE = ob.IEA15MW_SHAFT_TILT_DEG, ob.IEA15MW_PRECONE_DEG
CASES = [(10.659, 7.518, 0.0), (5.006, 5.0, 2.893), (12.259, 7.518, 6.76), (15.471, 7.518, 12.19)]


def ours_span(blade, case, linear=False, pitch_override=None, naz=36):
    s = p.BEMSolver(
        blade,
        rho=t5.RHO,
        mu=t5.MU_AD,
        precone=CONE,
        tilt=TILT,
        yaw=0.0,
        hub_height=ob.IEA15MW_HUB_HEIGHT,
        shear_exp=0.0,
    )
    if linear:
        t5._patch_linear(s, blade)
    pitch = case[2] if pitch_override is None else pitch_override
    acc = {}
    for az in np.arange(naz) * 360.0 / naz:
        l, _ = s.rotor.distributedAeroLoads(case[0], case[1], pitch, az)
        for k in ("alpha", "Cl", "Cd", "a", "ap", "Np", "Tp"):
            acc.setdefault(k, []).append(np.asarray(l[k]))
    m = {k: np.mean(v, axis=0) for k, v in acc.items()}
    r = s.compute(case[0], case[1], pitch)
    m["T"], m["Q"] = r.thrust, r.torque
    return m, blade.r


def main():
    out = {}
    blade = None
    for case in CASES:
        blade, ref = t5.ad_run(case, tilt=TILT, precone=CONE)
        m, r = ours_span(blade, case)
        ml, _ = ours_span(blade, case, linear=True)
        n = len(ref.alpha_deg)
        key = f"{case[0]}_{case[1]}_{case[2]}"
        print(
            f"\n=== {key}: nAD={n} nours={len(r)}  T ours/AD {m['T'] / ref.thrust - 1:+.4%} Q {m['Q'] / ref.torque - 1:+.4%}"
            f"  | linear-polar T {ml['T'] / ref.thrust - 1:+.4%} Q {ml['Q'] / ref.torque - 1:+.4%}"
        )
        rec = {
            "T": m["T"] / ref.thrust - 1,
            "Q": m["Q"] / ref.torque - 1,
            "Tlin": ml["T"] / ref.thrust - 1,
            "Qlin": ml["Q"] / ref.torque - 1,
        }
        if n == len(r):
            print(
                " i    r    a_AD  a_our  dα   Cl_AD Cl_our  Cd_AD  Cd_our  axI_AD axI_our tnI_AD tnI_our  dFn%   dFt%"
            )
            rows = []
            for i in range(n):
                fn_o, ft_o = m["Np"][i], m["Tp"][i]
                aa = np.deg2rad(ref.alpha_deg[i])
                cl_a = ref.cn[i] * np.cos(aa) + ref.ct_section[i] * np.sin(aa)
                cd_a = ref.cn[i] * np.sin(aa) - ref.ct_section[i] * np.cos(aa)
                _cn_o, _ct_o = cl_a, cd_a  # AD-derived Cl, Cd (columns reused)
                rows.append(
                    [
                        r[i],
                        ref.alpha_deg[i],
                        m["alpha"][i],
                        cl_a,
                        m["Cl"][i],
                        cd_a,
                        m["Cd"][i],
                        ref.ax_ind[i],
                        m["a"][i],
                        ref.tn_ind[i],
                        m["ap"][i],
                        ref.fn[i],
                        fn_o,
                        ref.ft[i],
                        ft_o,
                    ]
                )
                if i % 3 == 0 or i == n - 1:
                    print(
                        f"{i:2d} {r[i]:6.1f} {ref.alpha_deg[i]:6.2f} {m['alpha'][i]:6.2f} {m['alpha'][i] - ref.alpha_deg[i]:+5.2f} "
                        f"{cl_a:6.3f} {m['Cl'][i]:6.3f} {cd_a:.4f} {m['Cd'][i]:.4f} "
                        f"{ref.ax_ind[i]:.4f} {m['a'][i]:.4f} {ref.tn_ind[i]:.4f} {m['ap'][i]:.4f} "
                        ""
                    )
            rec["rows"] = rows
        out[key] = rec
    Path("logs/l4_t3diag/span.json").write_text(json.dumps(out, default=float))


if __name__ == "__main__":
    main()
