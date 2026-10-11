"""Counterfactual: the skin reference surface moved from the OML to the laminate mid-surface (E8).

Every outer-skin node is shifted in its section plane, inward along the averaged normal of its skin
elements, by the mean half-thickness of those elements. Before and after, on the S-7 production
mesh, it prints the S-7 band twist ratio against an integral-to-z beam reference (linear in 1/GJ
and linear in GJ) and the flap/edge compliance at the tip and at 0.9 L under a tip load, beside
the beam compliance at 0.9 L from the BeamDyn 6x6 rotated by the ElastoDyn structural twist.
Nothing in the solver changes: it measures what the reference-surface placement alone moves.

    scripts/aeroenv.sh python tools/diagnose_oml_reference_surface.py <element_size>

Recorded in odd/tasks/s7-torsion-ratio.md section E8.
"""

import sys
import numpy as np
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))
import run_s7_torsion as s7  # noqa: E402
from aeroelast.core.bc import DirichletCondition, NodalLoad  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver  # noqa: E402

es = float(sys.argv[1])
T = 1e6
L = 117.0

# beam reference, integral to z, linear in 1/GJ (S-7's own trapezoid) and linear in GJ
lines = [l.strip() for l in open(s7.BEAMDYN_FILE)]
i = next(k for k, l in enumerate(lines) if "DISTRIBUTED PROPERTIES" in l) + 1
zb = []
gj = []
while i < len(lines):
    p = lines[i].split()
    if len(p) == 1:
        try:
            e = float(p[0])
        except ValueError:
            i += 1
            continue
        j = i + 1
        while not lines[j]:
            j += 1
        zb.append(e * L)
        gj.append(float(lines[j + 5].split()[5]))
        i = j + 6
    else:
        i += 1
zb = np.array(zb)
gj = np.array(gj)
v = zb <= 111.2
zb, gj = zb[v], gj[v]
zz = np.linspace(0, 111.2, 200001)


def cum(f):
    return np.concatenate([[0], np.cumsum(np.diff(zz) * 0.5 * (f[1:] + f[:-1]))])


ref_inv = cum(np.interp(zz, zb, 1 / gj))
ref_lin = cum(1 / np.interp(zz, zb, gj))


def beam(zq, ref):
    return np.rad2deg(T * np.interp(np.minimum(zq, 111.2), zz, ref))


def solve(mesh, props, loads):
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0, 0, 1),
            "properties": props,
        },
    }
    s = StaticLinearSolver(mesh, cfg)
    dpn = 6
    root = sorted(mesh.get_node_set("RootNodes").nodes)
    idx = mesh.node_id_to_index
    s.add_dirichlet_conditions(
        [DirichletCondition(sorted(idx[n] * dpn + d for n in root for d in range(dpn)), 0.0)]
    )
    dofs, vals = [], []
    for n, f in loads:
        for d in range(3):
            dofs.append(int(n) * dpn + d)
            vals.append(float(f[d]))
    s.add_nodal_loads([NodalLoad(dofs, vals)])
    return s.solve()


def measure(mesh, props, label):
    X = mesh.coords_array
    tip, fy = s7.tip_couple(mesh, T)
    F = np.zeros((len(tip), 3))
    F[:, 1] = fy
    u = solve(mesh, props, zip(tip, F, strict=True))
    z7, th7 = s7.section_twists_deg(mesh, u)
    m = (z7 >= 0.3 * L) & (z7 <= 0.9 * L)
    r_inv = np.mean(th7[m] / beam(z7[m], ref_inv))
    r_lin = np.mean(th7[m] / beam(z7[m], ref_lin))
    out = [
        f"{label}: twist ratio 0.3-0.9 (1/GJ-linear ref) {r_inv:.4f}  (GJ-linear ref) {r_lin:.4f}"
    ]
    for name, d in (("flap(y)", 1), ("edge(x)", 0)):
        P = 1e5
        f = np.zeros(3)
        f[d] = P / len(tip)
        uu = solve(mesh, props, [(n, f) for n in tip]).reshape(-1, 6)
        band = np.abs(X[:, 2] - 0.9 * L) < 0.6
        out.append(
            f"  tip {name}: compliance at tip {uu[tip, d].mean() / P:.6e}  at 0.9L {uu[band, d].mean() / P:.6e} m/N"
        )
    print("\n".join(out))
    return r_inv


ed = [l.split() for l in open(s7.REF_DIR / "IEA-15-240-RWT_ElastoDyn_blade.dat")]
h = next(k for k, l in enumerate(ed) if l and l[0] == "BlFract")
rr = []
for l in ed[h + 2 :]:
    try:
        rr.append([float(x) for x in l[:6]])
    except ValueError:
        break
rr = np.array(rr)
Kb = []
i = next(k for k, l in enumerate(lines) if "DISTRIBUTED PROPERTIES" in l) + 1
while i < len(lines):
    p = lines[i].split()
    if len(p) == 1:
        try:
            e = float(p[0])
        except ValueError:
            i += 1
            continue
        j = i + 1
        while not lines[j]:
            j += 1
        K = [[float(x) for x in lines[j + k].split()] for k in range(6)]
        Kb.append((e, K[3][3], K[4][4]))
        i = j + 6
    else:
        i += 1
Kb = np.array(Kb)
z0 = 0.9 * L
zq = np.linspace(0, z0, 100001)
eq = zq / L
th = np.deg2rad(np.interp(eq, rr[:, 0], rr[:, 2]))
for lbl, f in (
    ("1/EI-linear", lambda a: np.interp(eq, Kb[:, 0], 1 / a)),
    ("EI-linear", lambda a: 1 / np.interp(eq, Kb[:, 0], a)),
):
    cf, ce = f(Kb[:, 2]), f(Kb[:, 1])
    w = (L - zq) * (z0 - zq)
    dy = np.trapezoid(w * (np.cos(th) ** 2 * cf + np.sin(th) ** 2 * ce), zq)
    dx = np.trapezoid(w * (np.sin(th) ** 2 * cf + np.cos(th) ** 2 * ce), zq)
    print(f"beam ({lbl}) compliance at 0.9L under tip load: flap(y) {dy:.6e}  edge(x) {dx:.6e} m/N")
mesh, props = s7.build_mesh_and_properties(es)
X0 = mesh.coords_array.copy()
idx = mesh.node_id_to_index
r0 = measure(mesh, props, "OML (as built)")

outer = {idx[n] for n in mesh.get_node_set("allOuterShellNods").nodes}
th_of = {}
for name in mesh.element_sets_names:
    lam = props.get(name)
    if lam is None or not hasattr(lam, "total_thickness"):
        continue
    for e in mesh.get_element_set(name).elements:
        th_of[getattr(e, "id", e)] = lam.total_thickness
nsum = defaultdict(lambda: np.zeros(3))
tsum = defaultdict(float)
cnt = defaultdict(int)
for el in mesh.elements:
    ids = (
        [idx[n] for n in el.node_ids] if hasattr(el, "node_ids") else [idx[n.id] for n in el.nodes]
    )
    if not all(a in outer for a in ids) or el.id not in th_of:
        continue
    P = X0[ids]
    nrm = (
        np.cross(P[2] - P[0], P[3 if len(ids) == 4 else 1] - P[1])
        if len(ids) == 4
        else np.cross(P[1] - P[0], P[2] - P[0])
    )
    nrm /= np.linalg.norm(nrm)
    for a in ids:
        nsum[a] += nrm
        tsum[a] += th_of[el.id]
        cnt[a] += 1
# orient: inward = towards the section's centroid at the node's z, by majority over the whole skin
zr = np.round(X0[:, 2], 6)
cz = {}
for a in outer:
    cz.setdefault(zr[a], []).append(a)
cent = {z: X0[ns, :2].mean(0) for z, ns in cz.items()}
votes = sum(np.sign(np.dot(nsum[a][:2], cent[zr[a]] - X0[a, :2])) for a in nsum)
sign = 1.0 if votes > 0 else -1.0
print(
    f"skin nodes shifted {len(nsum)}; inward-vote {votes:+.0f} of {len(nsum)} -> sign {sign:+.0f}"
)
keep_te_flat = "--keep-te-flat" in sys.argv
te_flat_nodes = set()
if keep_te_flat:
    te_ids = set()
    for name in mesh.element_sets_names:
        if name.endswith("_TE_FLAT"):
            te_ids.update(getattr(e, "id", e) for e in mesh.get_element_set(name).elements)
    for el in mesh.elements:
        if el.id in te_ids:
            te_flat_nodes.update(idx[n] for n in el.node_ids)
    print(f"--keep-te-flat: {len(te_flat_nodes)} TE_FLAT nodes stay on the OML")
X1 = X0.copy()
shift = []
for a in nsum:
    if a in te_flat_nodes:
        continue
    n = sign * nsum[a]
    n[2] = 0.0
    n /= np.linalg.norm(n)
    h = 0.5 * tsum[a] / cnt[a]
    X1[a] += h * n
    shift.append(h)
print(f"shift half-thickness: mean {np.mean(shift) * 1e3:.1f} mm  max {np.max(shift) * 1e3:.1f} mm")
# Geometry audit of the shift: who dissented from the inward vote, did any skin element invert,
# and how the shifted ring area compares with the first-order A_OML - 1/2 * contour of t ds.
set_of = {}
for name in mesh.element_sets_names:
    if not (len(name) > 6 and name[:2].isdigit() and name[3:5].isdigit()):
        continue  # only the region sets, e.g. 03_05_HP_SPAR
    for e in mesh.get_element_set(name).elements:
        set_of[getattr(e, "id", e)] = name
node_sets = defaultdict(set)
flips = 0
n_skin = 0
flip_where = defaultdict(int)
flip_eta = []
for el in mesh.elements:
    ids = (
        [idx[n] for n in el.node_ids] if hasattr(el, "node_ids") else [idx[n.id] for n in el.nodes]
    )
    if not all(a in outer for a in ids) or el.id not in th_of:
        continue
    n_skin += 1
    for a in ids:
        node_sets[a].add(set_of.get(el.id, "?").split("_", 2)[-1])

    def _n(P):
        return (
            np.cross(P[2] - P[0], P[3] - P[1])
            if len(P) == 4
            else np.cross(P[1] - P[0], P[2] - P[0])
        )

    if np.dot(_n(X0[ids]), _n(X1[ids])) <= 0.0:
        flips += 1
        flip_where[set_of.get(el.id, "?").split("_", 2)[-1]] += 1
        flip_eta.append(X0[ids, 2].mean() / L)
dissent = [a for a in nsum if sign * np.dot(nsum[a][:2], cent[zr[a]] - X0[a, :2]) < 0]
tally = defaultdict(int)
for a in dissent:
    for nm in node_sets[a]:
        tally[nm] += 1
print(
    f"dissenting nodes {len(dissent)}; by skin region: {dict(sorted(tally.items(), key=lambda kv: -kv[1]))}"
)
from matplotlib.path import Path as _MplPath  # noqa: E402

ring_edges = defaultdict(lambda: defaultdict(list))
for el in mesh.elements:
    ids = [idx[n] for n in el.node_ids]
    for a, b in zip(ids, ids[1:] + ids[:1], strict=True):
        if a in outer and b in outer and zr[a] == zr[b]:
            ring_edges[zr[a]][a].append(b)
            ring_edges[zr[a]][b].append(a)
outside = 0
checked = 0
for adj in ring_edges.values():
    adj = {k: sorted(set(v)) for k, v in adj.items()}
    if any(len(v) != 2 for v in adj.values()):
        continue
    start = next(iter(adj))
    order, prev, cur = [start], None, start
    while True:
        nxt = [n for n in adj[cur] if n != prev][0]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
    if len(order) != len(adj):
        continue
    moved = [a for a in order if a in nsum and a not in te_flat_nodes]
    inside = _MplPath(X0[order, :2]).contains_points(X1[moved, :2])
    outside += int((~inside).sum())
    checked += len(moved)
print(f"shifted ring nodes checked {checked}; landed OUTSIDE their original OML ring: {outside}")
print(
    f"skin elements whose normal flips after the shift: {flips} of {n_skin}; by region {dict(flip_where)}"
)
if flip_eta:
    print(
        f"  flipped elements span eta {min(flip_eta):.3f}..{max(flip_eta):.3f}; inside 0.3-0.9: {sum(0.3 <= e <= 0.9 for e in flip_eta)}"
    )
d_eta = [X0[a, 2] / L for a in dissent]
if d_eta:
    print(
        f"  dissenting nodes span eta {min(d_eta):.3f}..{max(d_eta):.3f}; inside 0.3-0.9: {sum(0.3 <= e <= 0.9 for e in d_eta)}"
    )
mesh.coords_array = X1
assert np.allclose(mesh.coords_array, X1)
r1 = measure(mesh, props, "mid-surface (inset t/2)")
print(
    f"GJ change implied by the inset: {r0 / r1:.4f} (shell GJ_after/GJ_before = ratio_before/ratio_after)"
)
