"""T13 (L3 D3): in-plane section distortion ~0.12 at z 90-92 (L24), from stored W7/W4 data only (no solve).

dist_rel(z) = res_rigid / (rrms * |theta|) per 0.5 m slab (half width 0.25 m), tip-torque S-7 case, gmsh chord-scaled N15/20/30/40,
plus the W7 end-stiffened variant (G12 x1.28 for z >= 98) as a load/section probe.
    python tools/l3_t13_offline.py logs/chord_scaled/w7 [OUT.json]
"""

import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
out = {}
for tag, f in (
    ("N15", "tip_n15"),
    ("N20", "tip_n20"),
    ("N30", "tip_n30"),
    ("N40", "tip_n40"),
    ("N20_tipstiff1.28", "tip_n20_var"),
):
    r = json.load(open(d / f"{f}.json"))
    z, th, res, rr = (np.array(r[k]) for k in ("slab_z", "slab_th", "slab_res_rigid", "slab_rrms"))
    out[tag] = {
        "z": z.tolist(),
        "dist_rel": (res / (rr * abs(th))).tolist(),
        "res_abs": res.tolist(),
        "rrms": rr.tolist(),
        "th": th.tolist(),
    }
zs = np.arange(86.0, 96.01, 0.5)
print("dist_rel = rms residual of rigid fit / (rms r * theta), 0.5 m slabs")
print("z     " + "  ".join(f"{k:>9s}" for k in out))
for zq in zs:
    row = []
    for v in out.values():
        i = int(np.argmin(abs(np.array(v["z"]) - zq)))
        row.append(v["dist_rel"][i])
    print(f"{zq:5.1f} " + "  ".join(f"{x:9.3f}" for x in row))
print("peak (z, value) per run in 84-98:")
for k, v in out.items():
    z = np.array(v["z"])
    x = np.array(v["dist_rel"])
    m = (z >= 84) & (z <= 98)
    i = int(np.argmax(np.where(m, x, -1)))
    print(
        f"  {k}: z={z[i]:.2f} dist_rel={x[i]:.3f}  (80 m: {np.interp(80, z, x):.3f}, 96 m: {np.interp(96, z, x):.3f}, 100 m: {np.interp(100, z, x):.3f})"
    )
if len(sys.argv) > 2:
    json.dump(out, open(sys.argv[2], "w"))
