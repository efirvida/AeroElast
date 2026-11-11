"""S-7 torsion on the gmsh IEA-15 mesh: the SAME exported mesh and the SAME tip-ring nodal force vector in the shell solver and in CalculiX S8R.

usage (env V2_BLADE_MESHER=gmsh is forced here; V2_BLADE_CACHE_DIR reuses a cached mesh):
    python tools/ccx_s7_torsion_gmsh.py [element_size=0.25] [out.json]

Three twist profiles over the span, each the best-fit rigid rotation of 40 z slices (tools/run_s7_torsion.section_twists_deg):
the shell solver (MITC4), CalculiX S8R (write_ccx_mesh, quadratic=True) and the BeamDyn deck K[5,5] integrated as a torsion beam.
Reported: mean ratio over 0.3-0.9 of the span shell/beam, ccx/beam and shell/ccx. Load: 1000 N.m by section_torque_wall_forces on
the TipNodes ring, computed once and given to both solvers, so the two see the same force vector.
"""
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))
os.environ["V2_BLADE_MESHER"] = "gmsh"

from aeroelast.core.mesh.io.writers import write_ccx_mesh  # noqa: E402
from aeroelast.solvers.bem.force_projection import section_torque_wall_forces  # noqa: E402
from run_s7_torsion import beam_twist_profile_deg, section_twists_deg  # noqa: E402
from tests.support.ccx_io import fail_ccx, parse_frd_disp, run_ccx  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402
from tests.validation.blade.v2.test_v2_iea15mw_s7_torsion import run_torsion_case_v2  # noqa: E402

MOMENT = 1000.0
L = 117.0
es = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "logs/check/gmsh_gap/ccx_s7_torsion_gmsh.json"
out.parent.mkdir(parents=True, exist_ok=True)

mesh, props = build_v2_blade(es, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.5))
coords = mesh.coords_array
n = len(coords)
idx_of = mesh.node_id_to_index
tip = np.asarray(sorted(idx_of[nid] for nid in mesh.get_node_set("TipNodes").nodes.keys()), dtype=np.intp)
tip_nodes, f_tip = section_torque_wall_forces(mesh, tip, (0.0, 0.0, 1.0), props, MOMENT, coords=coords)
print(f"mesh {n} nodes, {mesh.elements_count} elements, tip ring {len(tip_nodes)} nodes, |F| sum {np.abs(f_tip).sum():.3f} N", flush=True)

# shell solver (its own routine builds the same force vector with the same function)
u_sh = np.asarray(run_torsion_case_v2(mesh, props, MOMENT), dtype=float)
zs, th_sh = section_twists_deg(mesh, u_sh)
print("shell solved", flush=True)

if os.environ.get("SKIP_CCX") == "1":   # element-size sequences of the shell only (the CalculiX agreement was measured at 0.25)
    json.dump({"element_size": es, "nodes": n, "z": zs.tolist(), "shell_deg": th_sh.tolist(), "ccx_deg": None,
               "beam_deg": np.asarray(beam_twist_profile_deg(MOMENT, zs)).tolist()}, open(out, "w"), indent=1)
    print("SKIP_CCX: shell twist profile written", out)
    raise SystemExit(0)
# CalculiX: the force vector above, one *CLOAD per tip node
node_ids = [nd.id for nd in mesh.nodes]
field = {node_ids[int(t)]: [float(f[0]), float(f[1]), float(f[2])] for t, f in zip(tip_nodes, f_tip, strict=True)}
ccx = shutil.which("ccx") or os.environ.get("CCX_BIN")
if not ccx:
    raise SystemExit("ccx not found")
with tempfile.TemporaryDirectory(prefix="ccx_s7_gmsh_", dir=os.environ.get("TMPDIR")) as tmp:
    deck = Path(tmp) / "blade.inp"
    write_ccx_mesh(mesh, str(deck), properties=props, boundary_nodeset="RootNodes", solver_type="LinearStatic",
                   quadratic=True, span_direction=(0.0, 0.0, 1.0), load_field=field)
    print("ccx deck written; running", flush=True)
    proc = run_ccx(deck, ccx)
    if proc.returncode != 0:
        fail_ccx(proc, deck)
    disp = parse_frd_disp(deck.with_suffix(".frd"), list(range(1, n + 1)))
u_cx = np.zeros(6 * n)
miss = 0
for i in range(n):
    d = disp.get(i + 1)
    if d is None:
        miss += 1
        continue
    u_cx[6 * i:6 * i + 3] = np.asarray(d, float)[:3]
zs2, th_cx = section_twists_deg(mesh, u_cx)
th_bm = beam_twist_profile_deg(MOMENT, zs)
mid = (zs >= 0.3 * L) & (zs <= 0.9 * L)
r = {"shell/beam": float(np.nanmean(th_sh[mid] / th_bm[mid])), "ccx/beam": float(np.nanmean(th_cx[mid] / th_bm[mid])),
     "shell/ccx": float(np.nanmean(th_sh[mid] / th_cx[mid]))}
print("nodes without a CCX displacement:", miss)
print("mean ratio over 0.3-0.9 span:", {k: round(v, 4) for k, v in r.items()})
for z, a, b, c in list(zip(zs, th_sh, th_cx, th_bm, strict=True))[::4]:
    print(f"z {z:7.2f}  shell {a:9.5f}  ccx {b:9.5f}  beam {c:9.5f} deg")
json.dump({"element_size": es, "nodes": n, "ratios": r, "z": zs.tolist(), "shell_deg": th_sh.tolist(), "ccx_deg": th_cx.tolist(),
           "beam_deg": np.asarray(th_bm).tolist(), "ccx_missing_nodes": miss}, open(out, "w"), indent=1)
