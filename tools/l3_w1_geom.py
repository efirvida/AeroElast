#!/usr/bin/env python
"""W1: dimensional comparison of the solid's 2D section (build_2d) with a SONATA mesh export (tools/l3_w1_sonata_export.py).
Registers the SONATA section on the solid's by rigid fit (rotation grid + Nelder-Mead on the IoU of the unions), then reports
per-material area, IoU, TE-zone content, web chart lines and writes a figure.
  python tools/l3_w1_geom.py ETA SONATA.json OUTPREFIX     (aeroenv, PYTHONPATH=$PWD)"""
import collections
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from shapely import affinity
from shapely.geometry import Polygon
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402


def solid_polys(eta, h=0.1, tmax=0.03, web_mode="bonded"):
    sec = T.blade_section(eta)
    m2 = T.build_2d(sec, h, tmax, web_mode)
    xy = m2["xy"]
    out = []
    for q in m2["quads"]:
        pg = Polygon(xy[q[0]]).buffer(0)
        if pg.area > 0:
            out.append((q[1], pg))
    return sec, m2, out


def sonata_polys(path):
    d = json.load(open(path))
    out = []
    for c in d["cells"]:
        pg = Polygon(c["xy"]).buffer(0)
        if pg.area > 0:
            out.append((d["names"][str(c["mat"])], pg))
    return d, out


def register(sol, son):
    U1 = unary_union([p for _, p in sol]).buffer(1e-4).buffer(-1e-4)
    U2 = unary_union([p for _, p in son]).buffer(1e-4).buffer(-1e-4)

    def iou(p, A=U2):
        a = affinity.translate(affinity.rotate(A, p[0], origin=(0, 0)), p[1], p[2])
        return a.intersection(U1).area / a.union(U1).area

    best = None
    c1 = U1.centroid
    for ang in np.arange(-180, 180, 2.0):
        a = affinity.rotate(U2, ang, origin=(0, 0))
        t = (c1.x - a.centroid.x, c1.y - a.centroid.y)
        v = iou([ang, *t])
        if best is None or v > best[0]:
            best = (v, [ang, *t])
    r = minimize(lambda p: -iou(p), best[1], method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 600})
    return r.x, -r.fun, U1, U2


def main():
    eta, jf, pre = float(sys.argv[1]), sys.argv[2], sys.argv[3]
    sec, m2, sol = solid_polys(eta)
    d, son = sonata_polys(jf)
    p, iou, U1, U2 = register(sol, son)
    tr = lambda g: affinity.translate(affinity.rotate(g, p[0], origin=(0, 0)), p[1], p[2])  # noqa: E731
    son = [(m, tr(g)) for m, g in son]
    U2 = tr(U2)
    print(f"eta {eta} {Path(jf).name}: rot {p[0]:.3f} deg  shift ({p[1]:.4f},{p[2]:.4f}) m  IoU {iou:.4f}")
    As, Ao = collections.defaultdict(float), collections.defaultdict(float)
    for m, g in sol:
        As[m] += g.area
    for m, g in son:
        Ao[m] += g.area
    print("material area [m2]   solid   SONATA   ratio S/solid")
    for m in sorted(set(As) | set(Ao)):
        print(f"  {m:22s} {As[m]:.6f} {Ao[m]:.6f} {Ao[m]/As[m] if As[m] else float('nan'):.3f}")
    print(f"  total {sum(As.values()):.6f} {sum(Ao.values()):.6f}")
    only_sol, only_son = U1.difference(U2), U2.difference(U1)
    print(f"area only in solid {only_sol.area:.3e}  only in SONATA {only_son.area:.3e}")
    # TE zone: x along chord -> solid frame: TE = sec.te_point; distance from the TE point
    te = np.array(sec.te_point)
    for r in (0.03, 0.06, 0.12):
        from shapely.geometry import Point

        zone = Point(*te).buffer(r)
        a1 = sum(g.intersection(zone).area for _, g in sol)
        a2 = collections.defaultdict(float)
        for m, g in son:
            a2[m] += g.intersection(zone).area
        print(f"  TE zone r={r}: solid {a1:.3e} | SONATA {sum(a2.values()):.3e} " + " ".join(f"{k}={v:.2e}" for k, v in a2.items()))
    # webs: foam cores (medium_density_foam) cells whose centroid lies between the caps; report solid chart lines
    for w in sec.webs:
        print(f"  solid {w['name']}: A {np.round(w['A'],4)} B {np.round(w['B'],4)} layers {[(m, round(t,5)) for m, t in w['layers']]}")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(16, 7))
    for a, (zx, zy, zr) in zip(ax, ((te[0], te[1], 0.12), (0, float(np.mean(m2["xy"][:, 1])), 1.6)), strict=False):
        for _, g in sol:
            a.plot(*g.exterior.xy, "b-", lw=0.4)
        for m, g in son:
            a.plot(*g.exterior.xy, "r-" if m != "adhesive" else "g-", lw=0.4)
        a.set_aspect("equal")
        if a is ax[0]:
            a.set_xlim(zx - zr, zx + zr)
            a.set_ylim(zy - zr, zy + zr)
        a.set_title(f"eta {eta}: blue solid, red SONATA ({Path(jf).name}), green adhesive")
    fig.savefig(pre + ".png", dpi=110)
    json.dump({"rot": p[0], "shift": [p[1], p[2]], "iou": iou}, open(pre + ".json", "w"))


if __name__ == "__main__":
    main()
