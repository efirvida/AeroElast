"""L3 T10 (W7): band FE/Bredt compliance ratio vs N and extraction choices, Richardson-type extrapolation, from logs/chord_scaled/w7/tip_n<N>.json.

    python tools/l3_w7_analyze.py logs/chord_scaled/w7 15 20 30 40
Ratio = int(1/GJ_Bredt dz over the band) / ((theta(z1)-theta(z0))/T); <1 = FE more compliant than its Bredt.
Extraction options: theta(z) from a linear / quadratic fit of the slab rotations within +-W (W = 1.0, 1.5, 3.0 m) or from the nearest slab (no fit);
Bredt cut at z+0.0137 (as L24) or at z; rigid-fit rotation (S-7 estimator) or affine antisymmetric rotation.
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

BANDS = [(100, 106), (102, 108), (98, 104), (96, 102), (104, 110), (90, 96)]
T = 1000.0


def theta_at(d, z, w, order, key):
    zs, th = np.array(d["slab_z"]), np.array(d[key])
    if w == 0:
        return float(th[np.argmin(abs(zs - z))])
    m = (abs(zs - z) <= w) & (zs <= d["zmax"] - 0.25)
    return float(np.polyval(np.polyfit(zs[m], th[m], order), z))


def ratio(d, z0, z1, w=1.5, order=1, key="slab_th", off="bredt_off0137"):
    zb = np.array(d["bredt_z"])
    zz = np.arange(z0, z1 + 1e-9, 0.5)
    comp_b = np.trapezoid(1.0 / np.interp(zz, zb, d[off]), zz)
    return comp_b / ((theta_at(d, z1, w, order, key) - theta_at(d, z0, w, order, key)) / T)


def richardson(ns, rs):
    """r(N) = r_inf + a N^-p through the three finest levels; returns (r_inf, p) or None if non-monotone / not bracketed."""
    n, r = np.array(ns[-3:], float), np.array(rs[-3:])
    d1, d2 = r[1] - r[0], r[2] - r[1]
    if d1 * d2 <= 0:
        return None

    def f(p):
        return (n[1] ** -p - n[0] ** -p) * d2 - (n[2] ** -p - n[1] ** -p) * d1

    try:
        p = brentq(f, 0.01, 12.0)
    except ValueError:
        return None
    a = d2 / (n[2] ** -p - n[1] ** -p)
    return float(r[2] - a * n[2] ** -p), float(p)


def main():
    root, ns = Path(sys.argv[1]), [int(x) for x in sys.argv[2:]]
    D = {n: json.load(open(root / f"tip_n{n}.json")) for n in ns}
    print("nodes:", {n: D[n]["n_nodes"] for n in ns}, "ring:", {n: D[n]["ring_nodes"] for n in ns})
    opts = {
        "lin+-1.5 (L24)": {"w": 1.5},
        "lin+-3.0": {"w": 3.0},
        "lin+-1.0": {"w": 1.0},
        "quad+-3.0": {"w": 3.0, "order": 2},
        "nearest slab": {"w": 0},
        "lin+-1.5 affine": {"w": 1.5, "key": "slab_th_aff"},
        "lin+-1.5 Bredt@z": {"w": 1.5, "off": "bredt_off0"},
    }
    for z0, z1 in BANDS:
        print(f"\nband {z0}-{z1}:  " + "  ".join(f"N{n:>2d}" for n in ns) + "   Richardson(finest 3): r_inf, p")
        for name, kw in opts.items():
            rs = [ratio(D[n], z0, z1, **kw) for n in ns]
            rc = richardson(ns, rs) if len(ns) >= 3 else None
            tail = "non-monotone/none" if rc is None else f"{rc[0]:.3f}, p={rc[1]:.2f}"
            print(f"  {name:18s}" + "  ".join(f"{r:.3f}" for r in rs) + f"   {tail}")


if __name__ == "__main__":
    main()
