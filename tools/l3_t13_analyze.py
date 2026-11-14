"""T13: analyse logs/chord_scaled/t13/band_*.json: class shares of the residual and band compliance FE/Bredt (rigid fit, line fit +-1.5 m).
python tools/l3_t13_analyze.py DIR
"""

import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
TQ = 1000.0


def load(f):
    r = json.load(open(d / f))
    return r


def theta_at(r, z, W=1.5):
    zs, th = np.array(r["slab_z"]), np.array(r["slab_th"])
    m = abs(zs - z) <= W
    return float(np.polyval(np.polyfit(zs[m], th[m], 1), z))


def band(r, z0, z1):
    zb, gb = np.array(r["bredt_z"]), np.array(r["bredt_off0137"])
    zz = np.arange(z0, z1 + 1e-9, 0.5)
    comp = float(np.trapezoid(1.0 / np.interp(zz, zb, gb), zz))
    return comp / ((theta_at(r, z1) - theta_at(r, z0)) / TQ)


names = sorted(p.name for p in d.glob("band_*.json"))
for f in names:
    r = load(f)
    print(f, "variant", r["variant"], "nodes", r["n_nodes"])
    print(
        "  band FE/Bredt compliance:",
        {
            f"{a}-{b}": round(band(r, a, b), 3)
            for a, b in ((80, 88), (88, 92), (90, 96), (92, 96), (96, 102), (100, 106))
        },
    )
    print(
        "  twist(z)/T [deg/kNm] at 88/92/96/100:",
        [round(float(np.degrees(theta_at(r, z))), 5) for z in (88, 92, 96, 100)],
    )
    if f.endswith("base.json"):
        print("  residual share by node class (res2_share / node_share), z: skin | panel | cap")
        for z, c in r["class_shares"].items():
            print(
                f"   {z}: "
                + " | ".join(
                    f"{k} {v['res2_share']:.2f}/{v['node_share']:.2f}" for k, v in c.items()
                )
            )
