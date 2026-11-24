#!/usr/bin/env python
"""W8 / L3 T11: compare the warping of the 3D torque solution (middle plane) with the 2D St-Venant warping function on the SAME 2D mesh,
to locate where the 3D and 2D torsion stiffness part. Rebuilds the 2D mesh with the build arguments of the 3D case (read from CASE_DIR/meta.json
is not enough: pass ETA ARC H TMAX). Run: python tools/l3_w8_warp_compare.py CASE_DIR ETA ARC H TMAX"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402

d = Path(sys.argv[1])
eta, arc, h, tmax = (float(x) for x in sys.argv[2:6])
meta = json.load(open(d / "meta.json"))
sec = T.blade_section(eta)
m2 = T.build_2d(sec, h, tmax, "bonded", {"mat": "adhesive", "arc": arc} if arc > 0 else None)
ph = []
gj2 = T.torsion_2d(sec, m2, ph)
phi = ph[0]
blk = T.parse_dat(d / "case_T.dat", 1)[0]
xy, w, N2 = np.array(meta["corner_xy"]), np.array(meta["wnode"]), meta["N2"]
assert N2 == len(m2["xy"])
levels, zs = meta["levels"], meta["zs_levels"]
k_mid = levels[int(np.argmin([abs(z - 0.5 * meta["L"]) for z in zs]))]
k1, k2 = levels[2], levels[4]  # 0.4 L, 0.6 L
U = {k: np.array([blk[1 + k * N2 + i] for i in range(N2)]) for k in (k1, k_mid, k2)}
rate = np.array([T.fit_plane(xy, w, U[k][:, 0], U[k][:, 1], U[k][:, 2], meta["cx"], meta["cy"])[1] for k in (k1, k2)])
th = (rate[1] - rate[0]) / (zs[levels.index(k2)] - zs[levels.index(k1)])
X = np.column_stack([np.ones(N2), xy[:, 0], xy[:, 1]])
Wm = np.sqrt(w)[:, None]


def detrend(v):
    return v - X @ np.linalg.lstsq(X * Wm, v * Wm[:, 0], rcond=None)[0]


uz3 = detrend(U[k_mid][:, 2]) / th
uz2 = detrend(phi)
diff = uz3 - uz2
print(f"theta' 3D (0.4L-0.6L) {th:.4e} -> GJ_3D(T=1e6) {1e6 / th:.4e}; GJ_2D {gj2:.4e}; rms warping 2D {np.sqrt(np.sum(w * uz2**2) / w.sum()):.3e}, "
      f"rms(3D-2D) {np.sqrt(np.sum(w * diff**2) / w.sum()):.3e}")
sk = m2["node_sk"]
print("worst nodes (weight-scaled warping difference): idx  x  y  s  depth  uz3  uz2")
order = np.argsort(-np.abs(diff) * np.sqrt(w))[:15]
for i in order:
    s_, dpt = sk.get(i, (float("nan"), float("nan")))
    print(f"{i:5d} {xy[i,0]:8.4f} {xy[i,1]:8.4f}  s {s_:8.5f} depth {dpt:.4f}  {uz3[i]:+.4e} {uz2[i]:+.4e}")
