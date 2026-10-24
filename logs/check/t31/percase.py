"""usage: percase.py <tree> <case> <element_size>   (cwd=<tree>, PYTHONPATH=<tree>.pp)
Prints one PERCASE line: tip-mean displacement components (ux edge, uy flap, uz axial, twist deg) and
W=1/2 sum f.u for AE and CCX on the common tip nodes, plus relative gaps.
tip_torsion = force couple of 5e6 N.m about z on the tip nodes (f_i = T (z x r_i)/sum|r|^2), same forces both sides."""
import sys, os, tempfile, numpy as np
from pathlib import Path
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
tree, case, es = sys.argv[1], sys.argv[2], float(sys.argv[3])
sys.path.insert(0, tree)
import tests.test_blade_ccx_parity as T
from aeroelast.core.mesh.io.writers import write_ccx_mesh
from tests.support.ccx_io import run_ccx
T.ELEMENT_SIZE = es
mesh, props, coords = T._blade()
tip = sorted(mesh.node_id_to_index[n.id] for n in mesh.get_node_set("tip").nodes.values())
if case == "tip_torsion":
    r = coords[tip][:, :2] - coords[tip][:, :2].mean(0); den = (r ** 2).sum()
    F = {i: np.array([-5e6 * r[k, 1] / den, 5e6 * r[k, 0] / den, 0.0]) for k, i in enumerate(tip)}
else:
    load = np.asarray(T.CASES[case][0], float); F = {i: load / len(tip) for i in tip}
asm = T.PyMeshAssembler.from_model(T._to_rust_mesh(mesh, props), props, list(T.SPAN_DIRECTION), None)
rows, cols, vals = asm.assemble_k(); n = asm.dofs_count
K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
f = np.zeros(n)
for i in tip: f[i * 6:i * 6 + 3] = F[i]
cl = {mesh.node_id_to_index[x.id] * 6 + d for x in mesh.get_node_set("RootNodes").nodes.values() for d in range(6)}
mask = np.ones(n, bool); mask[list(cl)] = False; free = np.where(mask)[0]
u = np.zeros(n); u[free] = spsolve(K[np.ix_(free, free)], f[free])
ae = {i: u[i * 6:i * 6 + 3] for i in tip}
d = Path(tempfile.mkdtemp(dir=os.environ["TMPDIR"])); inp = d / "b.inp"
write_ccx_mesh(mesh, str(inp), properties=props, boundary_nodeset="RootNodes", solver_type="LinearStatic",
    load_field={mesh.nodes[i].id: [float(c) for c in F[i]] for i in tip},
    span_direction=(0.0, 0.0, 1.0), quadratic=True)
p = run_ccx(inp, "ccx"); print("ccx rc", p.returncode, flush=True)
cd = T._ccx_tip_displacements(d / "b.frd", mesh, coords)
com = [i for i in tip if i in cd]
def comp(disp):
    D = np.array([disp[i] for i in com]); m = D.mean(0)
    r = coords[com][:, :2] - coords[com][:, :2].mean(0)
    tw = np.degrees((r[:, 0] * D[:, 1] - r[:, 1] * D[:, 0]).sum() / (r ** 2).sum())
    return m[0], m[1], m[2], tw
W = lambda disp: 0.5 * sum(F[i] @ np.asarray(disp[i]) for i in com)
a, c, wa, wc = comp(ae), comp(cd), W(ae), W(cd)
g = lambda x, y: (x / y - 1) * 100 if y else float("nan")
print(f"PERCASE tree={os.path.basename(tree)} case={case} es={es} matched={len(com)}/{len(tip)}")
print("PERCASE_AE  edge_ux=%.6e flap_uy=%.6e axial_uz=%.6e twist_deg=%.6e W=%.6e" % (*a, wa))
print("PERCASE_CCX edge_ux=%.6e flap_uy=%.6e axial_uz=%.6e twist_deg=%.6e W=%.6e" % (*c, wc))
print("PERCASE_GAP(AE/CCX-1,%%) edge=%.2f flap=%.2f axial=%.2f twist=%.2f W=%.2f" % (*[g(x, y) for x, y in zip(a, c)], g(wa, wc)))
