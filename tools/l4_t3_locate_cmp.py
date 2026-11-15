"""L4 T3 locate (step 1+2): per-station comparison ours (CCBlade, exact linear polars) vs AeroDyn nodal outputs, r > 60 m.
Reads logs/l4_t3locate/ad_<key>.json from l4_t3_locate_ad.py. Diagnostic only."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t3_span as sp  # noqa: E402
import l4_t5_options as t5  # noqa: E402

ob, p = t5.ob, t5.p
OUT = Path("logs/l4_t3locate")


def ours(blade, case, naz=36):
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
    t5._patch_linear(s, blade)
    rot = s.rotor
    case[1] * 2 * np.pi / 60
    acc = {}
    for az in np.arange(naz) * 360.0 / naz:
        l, _ = rot.distributedAeroLoads(case[0], case[1], case[2], az)
        vx, vy, *_ = rot._CCBlade__windComponents(case[0], case[1], np.deg2rad(az))
        l = dict(l)
        l["Vx"], l["Vy"] = vx, vy
        l["phi"] = np.rad2deg(np.deg2rad(l["alpha"]) + rot.theta + np.deg2rad(case[2]))
        for k in ("alpha", "a", "ap", "Vx", "Vy", "phi", "Cl", "Cd", "Np", "Tp"):
            acc.setdefault(k, []).append(np.asarray(l[k]))
    m = {k: np.mean(v, axis=0) for k, v in acc.items()}
    m["r_ccb"] = rot.r.copy()
    m["theta_deg"] = np.rad2deg(rot.theta)
    m["chord"] = rot.chord.copy()
    return m, s


def main():
    for case in sp.CASES[2:]:
        pass
    for case in [(15.471, 7.518, 12.19), (10.659, 7.518, 0.0)]:
        key = f"{case[0]}_{case[1]}_{case[2]}"
        ad = json.loads((OUT / f"ad_{key}.json").read_text())
        blade, _ = t5.ad_run(case, tilt=sp.TILT, precone=sp.CONE) if False else (None, None)
        driver_blade = ob.build_blade_aero_from_aerodyn(p.PRIMARY)
        m, s = ours(driver_blade, case)
        r_bl = driver_blade.r
        print(
            f"\n=== {key}   nodes: AD {len(ad['Alpha'])} ours {len(r_bl)}; r_blade[-3:]={r_bl[-3:]}  CCBlade r[-3:]={m['r_ccb'][-3:]} rotorR={s.rotor.rotorR:.4f}"
        )
        print(
            " i  r_blade  | alpha AD  our   d  | phi AD    our     d   | theta(AD) our_th+pitch | Vx AD   our   | Vy AD    our   | a AD   our   | ap AD  our | Curve  Chi"
        )
        for i in range(len(r_bl)):
            if r_bl[i] < 60:
                continue
            thAD = ad["Theta"][i]
            thO = m["theta_deg"][i] + case[2]
            print(
                f"{i:2d} {r_bl[i]:7.2f} | {ad['Alpha'][i]:6.3f} {m['alpha'][i]:6.3f} {m['alpha'][i] - ad['Alpha'][i]:+.3f} | "
                f"{ad['Phi'][i]:6.3f} {m['phi'][i]:6.3f} {m['phi'][i] - ad['Phi'][i]:+.3f} | {thAD:7.3f} {thO:7.3f} {thO - thAD:+.3f} | "
                f"{ad['Vx'][i]:6.2f} {m['Vx'][i]:6.2f} | {ad['Vy'][i]:6.2f} {m['Vy'][i]:6.2f} | "
                f"{ad['AxInd'][i]:.4f} {m['a'][i]:.4f} | {ad['TnInd'][i]:.4f} {m['ap'][i]:.4f} | {ad['Curve'][i]:+.3f} {ad['Chi'][i]:+.2f}"
            )


if __name__ == "__main__":
    main()
