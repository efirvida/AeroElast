#!/usr/bin/env python
"""W8 / L3 T11: is the 3D torque case in uniform St-Venant torsion in the middle? Per plane: affine-fit twist, in-plane distortion rms
(displacement minus best affine fit, relative to the rotation displacement rms), and twist rate between consecutive planes.
Run: python tools/l3_w8_torsion_field.py CASE_DIR (needs result of `l3_solid_section_ccx.py run`)."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402

d = Path(sys.argv[1])
meta = json.load(open(d / "meta.json"))
blk = T.parse_dat(d / "case_T.dat", 1)[0]
xy, w, N2 = np.array(meta["corner_xy"]), np.array(meta["wnode"]), meta["N2"]
prev = None
for k, z in zip(meta["levels"], meta["zs_levels"], strict=False):
    U = np.array([blk[1 + k * N2 + i] for i in range(N2)])
    cw, th = T.fit_plane(xy, w, U[:, 0], U[:, 1], U[:, 2], meta["cx"], meta["cy"])
    X = np.column_stack([np.ones(N2), xy[:, 0], xy[:, 1]])
    res = []
    for c in (0, 1):
        co = np.linalg.lstsq(X * np.sqrt(w)[:, None], U[:, c] * np.sqrt(w), rcond=None)[0]
        res.append(U[:, c] - X @ co)
    dist = np.sqrt(np.sum(w * (res[0] ** 2 + res[1] ** 2)) / w.sum())
    rot = abs(th) * np.sqrt(np.sum(w * ((xy[:, 0] - meta["cx"]) ** 2 + (xy[:, 1] - meta["cy"]) ** 2)) / w.sum())
    rate = "" if prev is None else f"  twist rate vs previous plane {(th - prev[1]) / (z - prev[0]):.4e}"
    print(f"z {z:7.3f} twist {th:.4e} in-plane distortion/rotation rms {dist / rot:.4f}{rate}")
    prev = (z, th)

# where is the in-plane distortion at the plane closest to mid-length (worst nodes by weight-scaled magnitude)
k = meta["levels"][int(np.argmin([abs(z - 0.5 * meta["L"]) for z in meta["zs_levels"]]))]
U = np.array([blk[1 + k * N2 + i] for i in range(N2)])
X = np.column_stack([np.ones(N2), xy[:, 0], xy[:, 1]])
mag = np.zeros(N2)
for c in (0, 1):
    co = np.linalg.lstsq(X * np.sqrt(w)[:, None], U[:, c] * np.sqrt(w), rcond=None)[0]
    mag += (U[:, c] - X @ co) ** 2
mag = np.sqrt(mag)
print("worst in-plane-distortion nodes at mid plane: idx x y |u-affine|")
for i in np.argsort(-mag * np.sqrt(w))[:8]:
    print(f"{i:5d} {xy[i, 0]:8.4f} {xy[i, 1]:8.4f} {mag[i]:.3e}")
