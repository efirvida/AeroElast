"""L3 T10 (W7): S-7 case (1000 N m wall-flow torque on TipNodes, root clamped) at gmsh chord-scaled N, AR 1; stores the raw slab twist
(rigid fit, 0.5 m slabs of half width 0.25 m, as tools/l3_w4_tip_compliance.py) and the mesh-cut Bredt GJ on a 0.5 m grid at two cut
offsets, so the band compliance ratio can be re-extracted offline (tools/l3_w7_analyze.py) with other fit windows / cut z.

    python tools/l3_w7_tip.py N OUT.json [--variant tipstiff:F:Z0:Z1]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l3_w4_tip_compliance as w4  # noqa: E402
from l3_w4_bredt import make  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade, load_raw  # noqa: E402


def main():
    a = sys.argv[1:]
    n, out = int(a[0]), a[1]
    variant = a[a.index("--variant") + 1] if "--variant" in a else None
    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", chord_divisions=n)
    print(f"mesh N={n} nodes={mesh.node_count} t={time.time() - t0:.0f}s", flush=True)
    if variant:
        import l3_w7_tipstiff as ts

        mesh, props = ts.from_spec(variant)(mesh, props)
    u, nring = w4.solve(mesh, props, None)
    t_solve = time.time() - t0
    coords = mesh.coords_array
    zmax = float(coords[:, 2].max())
    sl = [
        s
        for s in (w4.slab_stats(coords, u, 6, zc, 0.25) for zc in np.arange(58.0, zmax, 0.5))
        if s is not None
    ]
    raw = load_raw(0.25, mesher="gmsh", generate=False, chord_divisions=n)
    deck = str(REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    gj_at, _, _ = make(raw, deck)
    zb = np.arange(58.0, 116.01, 0.5)
    res = {
        "N": n,
        "variant": variant,
        "n_nodes": int(len(coords)),
        "zmax": zmax,
        "ring_nodes": int(nring),
        "torque": w4.TORQUE,
        "solve_s": round(t_solve, 1),
        "slab_z": [s["z"] for s in sl],
        "slab_th": [s["th"] for s in sl],
        "slab_th_aff": [s["th_aff"] for s in sl],
        "slab_res_rigid": [s["res_rigid"] for s in sl],
        "slab_rrms": [s["rrms"] for s in sl],
        "bredt_z": zb.tolist(),
        "bredt_off0137": [gj_at(z + 0.0137)[0] for z in zb],
        "bredt_off0": [gj_at(z + 0.0)[0] for z in zb],
        "total_s": round(time.time() - t0, 1),
    }
    json.dump(res, open(out, "w"))
    print("done", res["total_s"], flush=True)


if __name__ == "__main__":
    main()
