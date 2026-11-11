"""S-7 shell/beam twist ratio of the same IEA-15 gmsh mesh with and without the web0/web1 node weld (tools/weld_coincident_nodes.py).
usage: python tools/weld_effect_s7.py ORIGINAL.npz WELDED.npz   (both are raw v2 arrays; the mesh model is built with tests.support.v2_blade.mesh_from_raw)
Same load and estimator as tests/validation/blade/v2/test_v2_iea15mw_s7_torsion.py (1000 N.m by section_torque_wall_forces on the TipNodes ring,
section_twists_deg, ratio over 0.3-0.9 of the span against the BeamDyn K[5,5] profile)."""
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
os.environ["V2_BLADE_MESHER"] = "gmsh"
from run_s7_torsion import beam_twist_profile_deg, section_twists_deg  # noqa: E402
from tests.support.v2_blade import mesh_from_raw  # noqa: E402
from tests.validation.blade.v2.test_v2_iea15mw_s7_torsion import run_torsion_case_v2  # noqa: E402

res = {}
for label, path in (("original", sys.argv[1]), ("welded", sys.argv[2])):
    raw = dict(np.load(path, allow_pickle=False))
    raw["key"] = f"weld-effect-{label}"
    mesh, props = mesh_from_raw(raw, "mid")
    u = np.asarray(run_torsion_case_v2(mesh, props, 1000.0), float)
    zs, th = section_twists_deg(mesh, u)
    ref = beam_twist_profile_deg(1000.0, zs)
    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    res[label] = (zs, th)
    print(f"{label}: nodes {len(mesh.nodes)}  mean shell/beam over 0.3-0.9 span {np.nanmean(th[mid] / ref[mid]):.5f}", flush=True)
zs, a = res["original"]
_, b = res["welded"]
ok = np.isfinite(a) & np.isfinite(b) & (a > 0)
print(f"max |welded/original - 1| over the slices: {np.max(np.abs(b[ok] / a[ok] - 1)) * 100:.4f} %  (z of the max {zs[ok][np.argmax(np.abs(b[ok] / a[ok] - 1))]:.1f} m)")
print("tip slice:", f"original {a[ok][-1]:.6f} welded {b[ok][-1]:.6f} deg")
