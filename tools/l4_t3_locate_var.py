"""L4 T3 locate (steps 2+3): does adding the deck's prebend/sweep (BlCrvAC, BlSwpAC) to our CCBlade, and/or removing the
tip node on both sides, remove the outboard alpha offset and the thrust/torque error? Diagnostic only; solver untouched."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t3_locate_cmp as cmp  # noqa: E402
import l4_t3_span as sp  # noqa: E402
import l4_t5_options as t5  # noqa: E402

ob, p = t5.ob, t5.p
OUT = Path("logs/l4_t3locate")
CASES = [(15.471, 7.518, 12.19), (10.659, 7.518, 0.0)]


def deck_curve():
    lines = p.BLADE.read_text(errors="replace").splitlines()
    h = next(i for i, l in enumerate(lines) if "BlSpn" in l and "BlAFID" in l)
    rows = []
    for l in lines[h + 1 :]:
        f = l.split()
        if len(f) >= 7 and not f[0].startswith("("):
            try:
                rows.append([float(x) for x in f[:7]])
            except ValueError:
                continue
        if len(rows) == 50:
            break
    d = np.array(rows)
    return d[:, 1], d[:, 2], d[:, 3]  # BlCrvAC, BlSwpAC, BlCrvAng


def make(blade, curve_sign=0.0, sweep_sign=0.0, notip=False):
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
    crv, swp, _ = deck_curve()
    r_nd = (blade.r - blade.r[0]) / (blade.r[-1] - blade.r[0])
    nd_tip = max(0.99, r_nd[-2:].mean())
    if curve_sign:
        rot.precurve = curve_sign * crv.copy()
        rot.precurveTip = curve_sign * crv[-1]
        rot.precurve[-1] = np.interp(nd_tip, r_nd, rot.precurve)
        rot.rotorR = rot.Rtip * np.cos(rot.precone) + rot.precurveTip * np.sin(rot.precone)
    if sweep_sign:
        rot.presweep = sweep_sign * swp.copy()
        rot.presweepTip = sweep_sign * swp[-1]
        rot.presweep[-1] = np.interp(nd_tip, r_nd, rot.presweep)
    if notip:
        rot.chord = rot.chord.copy()
        rot.chord[-1] = 1e-9
    return s


def tq(s, case, d=0.0):
    r = s.compute(case[0], case[1], case[2] + d)
    return r.thrust, r.torque


def main():
    res = {}
    for case in CASES:
        key = f"{case[0]}_{case[1]}_{case[2]}"
        ad = json.loads((OUT / f"ad_{key}.json").read_text())
        blade = ob.build_blade_aero_from_aerodyn(p.PRIMARY)
        r = blade.r
        # AD tip-node share of thrust / torque: B * F * w_tip, w_tip = half spacing; rotor-frame Fx (axial), Fy (tangential)
        0.5 * (r[-1] - r[-2])
        B = 3
        print(f"\n=== {key}: AD RtFldFxh {ad['_thrust']:.1f} Mxh {ad['_torque']:.1f}")
        # reconstruct AD total from nodal Fx, Fy to validate the share estimate
        wi = np.gradient(r)  # not exact; use trapezoid-like widths
        wi = np.r_[0.5 * (r[1] - r[0]), 0.5 * (r[2:] - r[:-2]), 0.5 * (r[-1] - r[-2])]
        Fx, Fy, Fn, Ft = (np.array(ad[k]) for k in ("Fx", "Fy", "Fn", "Ft"))
        print(
            f"   AD nodal tip Fx,Fy,Fn,Ft = {Fx[-1]:.1f} {Fy[-1]:.1f} {Fn[-1]:.1f} {Ft[-1]:.1f} N/m; Fx[-2]={Fx[-2]:.1f}"
        )
        for nm, arr in (("Fx", Fx), ("Fn", Fn)):
            print(
                f"   B*sum({nm} w)= {B * np.sum(arr * wi):.1f}  vs RtFldFxh {ad['_thrust']:.1f} ({B * np.sum(arr * wi) / ad['_thrust'] - 1:+.3%})"
            )
        torq_arm = r * np.cos(np.deg2rad(sp.CONE))
        for nm, arr in (("Fy", Fy), ("Ft", Ft)):
            for sgn in (1, -1):
                v = sgn * B * np.sum(arr * wi * torq_arm)
                print(
                    f"   {sgn:+d}*B*sum({nm} w r cos(cone))={v:.1f} vs Mxh {ad['_torque']:.1f} ({v / ad['_torque'] - 1:+.3%})"
                )
        # tip-node share (use the Fx/Fy that best reconstruct, printed above)
        sh_T = B * Fx[-1] * wi[-1]
        sh_Q = B * Fy[-1] * wi[-1] * torq_arm[-1]
        T_ad, Q_ad = ad["_thrust"], ad["_torque"]
        T_ad_nt, Q_ad_nt = T_ad - sh_T, Q_ad - abs(sh_Q) * np.sign(Q_ad)
        print(
            f"   AD tip share  T {sh_T:.1f} ({sh_T / T_ad:+.3%})  Q {sh_Q:.1f} ({sh_Q / Q_ad:+.3%})"
        )
        for name, kw in [
            ("base(lin)", {}),
            ("crv-", {"curve_sign": -1}),
            ("crv+", {"curve_sign": 1}),
            ("crv-,swp-", {"curve_sign": -1, "sweep_sign": -1}),
            ("crv-,swp+", {"curve_sign": -1, "sweep_sign": 1}),
        ]:
            for notip in (False, True):
                s = make(blade, notip=notip, **kw)
                m, _ = cmp.ours(blade, case) if False else (None, None)
                T, Q = tq(s, case)
                tref, qref = (T_ad_nt, Q_ad_nt) if notip else (T_ad, Q_ad)
                # sensitivity -> pitch shift zeroing thrust
                ds = np.array([-0.3, 0.0, 0.3])
                eT = np.array([tq(s, case, d)[0] / tref - 1 for d in ds])
                zT = -eT[1] / np.polyfit(ds, eT, 1)[0]
                # outboard alpha offset r in 80..115 (azimuth mean)
                rot = s.rotor
                acc = []
                for az in np.arange(36) * 10.0:
                    l, _ = rot.distributedAeroLoads(case[0], case[1], case[2], az)
                    acc.append(l["alpha"])
                al = np.mean(acc, axis=0)
                msk = (r > 80) & (r < 116)
                da = (al - np.array(ad["Alpha"]))[msk]
                print(
                    f"   {name:10s} notip={notip!s:5}: T {T / tref - 1:+.3%} Q {Q / qref - 1:+.3%}  pitch-to-zero-T {zT:+.3f} deg  "
                    f"dalpha(80<r<116) mean {da.mean():+.3f} [{da.min():+.3f},{da.max():+.3f}]  dalpha(r=61..80) {np.mean((al - np.array(ad['Alpha']))[(r > 60) & (r <= 80)]):+.3f}"
                )
                res[f"{key}|{name}|{notip}"] = {
                    "T": T / tref - 1,
                    "Q": Q / qref - 1,
                    "zT": zT,
                    "da": float(da.mean()),
                }
    (OUT / "var.json").write_text(json.dumps(res))


if __name__ == "__main__":
    main()
