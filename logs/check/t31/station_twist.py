"""usage: station_twist.py <tree> <frd> <element_size> [case=tip_torsion]   (cwd=<tree>, PYTHONPATH=<tree>.pp)
T31 follow-up: is the tip_torsion W gap a GJ difference or loaded-end distortion?
Re-solves our side of tip_torsion (same forces as percase.py), reads the preserved CCX FRD of the same
case, and fits the rigid section rotation theta(z) = sum(r x u)_z / sum r^2 on slabs of nodes at interior
stations. Interior twist rate ratio dtheta/dz (AE/CCX) = GJ_CCX/GJ_AE away from the loaded tip section."""
import sys, numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
tree, frd, es = sys.argv[1], sys.argv[2], float(sys.argv[3])
case = sys.argv[4] if len(sys.argv) > 4 else "tip_torsion"
sys.path.insert(0, tree)
import tests.test_blade_ccx_parity as T
from tests.validation.parity.test_orthotropic_shell_parity import _parse_ccx_frd, _parse_ccx_frd_coords
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
UA = u.reshape(-1, 6)[:, :3]
cc = _parse_ccx_frd_coords(frd); ids = list(cc)
cd = _parse_ccx_frd(frd, ids)
dist, idx = cKDTree(np.array([cc[i] for i in ids])).query(coords)
print(f"corner match: max dist {dist.max():.2e} m, n {len(coords)}")
UC = np.array([cd[ids[j]] for j in idx])
z = coords[:, 2]; z0, z1 = z.min(), z.max()
def theta(U, sel):
    p = coords[sel][:, :2]; rr = p - p.mean(0)
    return np.degrees((rr[:, 0] * U[sel, 1] - rr[:, 1] * U[sel, 0]).sum() / (rr ** 2).sum())
print(f"span z {z0:.2f}..{z1:.2f} m, tip nodes {len(tip)}")
print(f"case {case}")
print("| z (m) | nodes | theta AE (deg) | theta CCX (deg) | AE/CCX-1 % | mean u AE/CCX-1 % (x, y, z) |")
print("|---|---|---|---|---|---|")
st = []
for frac in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99):
    zc = z0 + frac * (z1 - z0); sel = np.where(np.abs(z - zc) < 0.6 * es)[0]
    ta, tc = theta(UA, sel), theta(UC, sel); st.append((zc, ta, tc))
    ma, mc = UA[sel].mean(0), UC[sel].mean(0)
    print(f"| {zc:.1f} | {len(sel)} | {ta:.6e} | {tc:.6e} | {(ta / tc - 1) * 100:.2f} | " + ", ".join(f"{(a / c - 1) * 100:.2f}" for a, c in zip(ma, mc)) + " |")
print("| band (m) | dtheta/dz AE/CCX-1 % (= GJ_CCX/GJ_AE - 1) |")
print("|---|---|")
for (za, aa, ca), (zb, ab, cb) in zip(st, st[1:]):
    print(f"| {za:.1f}-{zb:.1f} | {((ab - aa) / (cb - ca) - 1) * 100:.2f} |")
