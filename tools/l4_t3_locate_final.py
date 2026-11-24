"""L4 T3 locate (final): our CCBlade with (a) deck prebend as precurve (local cant) and (b) AeroDyn's pitch-through-cant
effective twist  theta_eff = twist + pitch/cos(c)  [c = AD 'Curve' output], with/without the tip node, vs AeroDyn.
Diagnostic only; the solver is not changed. Output: logs/l4_t3locate/final.log."""

import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t3_locate_var as v  # noqa: E402

ob, p, sp, OUT = v.ob, v.p, v.sp, v.OUT


def main():
    for case in v.CASES:
        key = f"{case[0]}_{case[1]}_{case[2]}"
        ad = json.loads((OUT / f"ad_{key}.json").read_text())
        blade = ob.build_blade_aero_from_aerodyn(p.PRIMARY)
        r = blade.r
        cant = np.deg2rad(np.abs(np.array(ad["Curve"])))
        wi = np.r_[0.5 * (r[1] - r[0]), 0.5 * (r[2:] - r[:-2]), 0.5 * (r[-1] - r[-2])]
        Fx, Fy = np.array(ad["Fx"]), np.array(ad["Fy"])
        arm = r * np.cos(np.deg2rad(sp.CONE))
        shT = Fx[-1] * wi[-1] / np.sum(Fx * wi)
        shQ = Fy[-1] * wi[-1] * arm[-1] / np.sum(Fy * wi * arm)
        T_ad, Q_ad = ad["_thrust"], ad["_torque"]
        print(f"\n=== {key}  AD tip-node share (nodal recon): T {shT:+.3%} Q {shQ:+.3%}")
        for name, kw, cantfix in [
            ("base(lin)", {}, False),
            ("crv-", {"curve_sign": -1}, False),
            ("crv-+pitch/cos(cant)", {"curve_sign": -1}, True),
        ]:
            for notip in (False, True):
                s = v.make(blade, notip=notip, **kw)
                rot = s.rotor
                th0 = rot.theta.copy()

                def run(d):
                    pit = case[2] + d
                    rot.theta = th0 + (np.deg2rad(pit) * (1 / np.cos(cant) - 1) if cantfix else 0.0)
                    res = s.compute(case[0], case[1], pit)
                    return res.thrust, res.torque

                tref = T_ad * (1 - shT) if notip else T_ad
                qref = Q_ad * (1 - shQ) if notip else Q_ad
                T, Q = run(0.0)
                ds = np.array([-0.3, 0, 0.3])
                eT = np.array([run(d)[0] / tref - 1 for d in ds])
                eQ = np.array([run(d)[1] / qref - 1 for d in ds])
                zT = -eT[1] / np.polyfit(ds, eT, 1)[0]
                zQ = -eQ[1] / np.polyfit(ds, eQ, 1)[0]
                run(0.0)
                acc = []
                for az in np.arange(36) * 10.0:
                    l, _ = rot.distributedAeroLoads(case[0], case[1], case[2], az)
                    acc.append(l["alpha"])
                da = np.mean(acc, axis=0) - np.array(ad["Alpha"])
                mk = (r > 80) & (r < 116)
                print(
                    f"  {name:22s} notip={notip!s:5}: T {T / tref - 1:+.3%} Q {Q / qref - 1:+.3%} | pitch-to-zero T {zT:+.3f} Q {zQ:+.3f} deg | "
                    f"dalpha(80-116) mean {da[mk].mean():+.3f} max|.| {np.abs(da[mk]).max():.3f}; (61-80) {da[(r > 60) & (r <= 80)].mean():+.3f}"
                )
    print(
        "\nnote: base(lin) notip row differs slightly from l4_t3_locate_var.py because the AD tip share is here scaled by the nodal-sum reconstruction."
    )


if __name__ == "__main__":
    main()
