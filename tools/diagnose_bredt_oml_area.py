"""First-order Bredt estimate of the OML reference-surface effect on GJ (E8, 2026-10-11).

The blade mesh puts its skin nodes on the outer mould line and the laminate is centred on them,
so half of each skin laminate lies outside the aerofoil and the torsion cell encloses the OML area.
If the laminate grew inward from the OML, the shear-flow path would be the mid-surface, and to
first order A_mid = A_OML - 1/2 * contour integral of t ds. This prints (A_mid / A_OML)^2, the
single-cell Bredt factor on GJ, ring by ring over 0.3-0.9 R.

    scripts/aeroenv.sh python tools/diagnose_bredt_oml_area.py <element_size>

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

es = float(sys.argv[1])
L = 117.0
mesh, props = s7.build_mesh_and_properties(es)
X = mesh.coords_array
idx = mesh.node_id_to_index
outer = {idx[n] for n in mesh.get_node_set("allOuterShellNods").nodes}
zr = np.round(X[:, 2], 6)
# element thickness through its set
th_of = {}
for name in mesh.element_sets_names:
    lam = props.get(name)
    if lam is None:
        continue
    t = lam.total_thickness if hasattr(lam, "total_thickness") else float(lam.get("thickness", 0))
    for e in [getattr(x, "id", x) for x in mesh.get_element_set(name).elements]:
        th_of[e] = t
edges = defaultdict(list)  # ring edge -> thicknesses
for el in mesh.elements:
    ids = (
        [idx[n] for n in el.node_ids] if hasattr(el, "node_ids") else [idx[n.id] for n in el.nodes]
    )
    t = th_of.get(el.id)
    if t is None:
        continue
    for a, b in zip(ids, ids[1:] + ids[:1], strict=True):
        if a in outer and b in outer and zr[a] == zr[b]:
            edges[(min(a, b), max(a, b))].append(t)
rings = defaultdict(list)
for (a, b), ts in edges.items():
    rings[zr[a]].append((a, b, np.mean(ts)))
print(f"es={es}: {len(rings)} rings, {len(th_of)} elements with thickness")
rows = []
for z, E in sorted(rings.items()):
    if not (0.3 * L <= z <= 0.9 * L):
        continue
    adj = defaultdict(list)
    for a, b, _t in E:
        adj[a].append(b)
        adj[b].append(a)
    if any(len(v) != 2 for v in adj.values()):
        continue
    start = E[0][0]
    order = [start]
    prev = None
    cur = start
    while True:
        nxt = [n for n in adj[cur] if n != prev][0]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
    if len(order) != len(adj):
        continue
    P = X[order, :2]
    A = 0.5 * abs(np.dot(P[:, 0], np.roll(P[:, 1], -1)) - np.dot(P[:, 1], np.roll(P[:, 0], -1)))
    tmap = {(min(a, b), max(a, b)): t for a, b, t in E}
    half_t_ds = sum(
        0.5 * tmap[(min(a, b), max(a, b))] * np.linalg.norm(X[a, :2] - X[b, :2])
        for a, b in zip(order, order[1:] + order[:1], strict=True)
    )
    rows.append((z / L, A, half_t_ds / A, ((A - half_t_ds) / A) ** 2))
r = np.array(rows)
for q in r[:: max(1, len(r) // 12)]:
    print(f"  eta {q[0]:.2f}  A_OML {q[1]:.3f} m2  dA/A {q[2]:.3f}  GJ factor {q[3]:.3f}")
print(f"  band mean GJ factor (A_mid/A_OML)^2 = {r[:, 3].mean():.3f}  over {len(r)} rings")
