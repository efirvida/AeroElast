"""L4 T3 locate: per-station detail ours(crv-,[swp]) vs AD at one case; shows what the prebend model leaves. Diagnostic only."""

import json
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t3_locate_var as v  # noqa: E402

ob, p, sp = v.ob, v.p, v.sp
OUT = v.OUT
for case in v.CASES:
    key = f"{case[0]}_{case[1]}_{case[2]}"
    ad = json.loads((OUT / f"ad_{key}.json").read_text())
    blade = ob.build_blade_aero_from_aerodyn(p.PRIMARY)
    r = blade.r
    for name, kw in (
        ("crv-", {"curve_sign": -1}),
        ("crv-,swp+", {"curve_sign": -1, "sweep_sign": 1}),
    ):
        s = v.make(blade, **kw)
        rot = s.rotor
        om = case[1]
        acc = {}
        for az in np.arange(36) * 10.0:
            l, _ = rot.distributedAeroLoads(case[0], case[1], case[2], az)
            vx, vy, *_ = rot._CCBlade__windComponents(case[0], om, np.deg2rad(az))
            for k, x in (
                ("al", l["alpha"]),
                ("a", l["a"]),
                ("ap", l["ap"]),
                ("Vx", vx),
                ("Vy", vy),
            ):
                acc.setdefault(k, []).append(np.asarray(x))
        m = {k: np.mean(x, axis=0) for k, x in acc.items()}
        phi = np.rad2deg(np.deg2rad(m["al"]) + rot.theta + np.deg2rad(case[2]))
        print(
            f"\n=== {key} [{name}]   i  r  | dalpha  dphi | Vx AD/our  | Vy AD/our | a AD/our | ap AD/our"
        )
        for i in range(24, 49, 3):
            print(
                f"{i:2d} {r[i]:6.1f} | {m['al'][i] - ad['Alpha'][i]:+.3f} {phi[i] - ad['Phi'][i]:+.3f} | {ad['Vx'][i]:6.3f} {m['Vx'][i]:6.3f} | "
                f"{ad['Vy'][i]:7.3f} {m['Vy'][i]:7.3f} | {ad['AxInd'][i]:.4f} {m['a'][i]:.4f} | {ad['TnInd'][i]:.5f} {m['ap'][i]:.5f}"
            )
