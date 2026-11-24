"""T13: layup events of the gmsh chord-scaled mesh at z 88-98 (no solve).

Per 0.5 m z bin of the element centroid: number of elements, shell thickness (laminate total_thickness) min/median/max, n_plies
min/max, number of distinct layup sets, and the thickness split by chord position of the element centroid (chord axis from the PCA
of the slab nodes, s = (x_c - x_min)/(x_max - x_min), 0 = one end of the chord): aft 0-0.15, mid 0.15-0.85, fore 0.85-1 (end labels
are PCA-arbitrary; both ends are reported by min/max over the two end regions).
    python tools/l3_t13_mesh_events.py N OUT.json
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.support.v2_blade import build_v2_blade  # noqa: E402

n, out = int(sys.argv[1]), sys.argv[2]
mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", chord_divisions=n)
pos = {
    nd.id: np.asarray(nd.coordinates if hasattr(nd, "coordinates") else nd.coords, float)
    for nd in mesh.nodes
}
rows = {}
for name, es in mesh.element_sets.items():
    if name not in props:
        continue
    lam = props[name]
    for e in es.elements:
        c = np.mean([pos[i] for i in e.node_ids], axis=0)
        rows.setdefault(int(np.floor(c[2] / 0.5)), []).append(
            (name, float(lam.total_thickness), int(lam.n_plies), c)
        )
res = []
print(
    " z     n_el  t_min   t_med   t_max  plies  nsets | end A t_med/n   end B t_med/n | mid t_med"
)
for b in sorted(rows):
    z = (b + 0.5) * 0.5
    if not 86 <= z <= 98:
        continue
    r = rows[b]
    t = np.array([x[1] for x in r])
    npl = np.array([x[2] for x in r])
    C = np.array([x[3][:2] for x in r])
    cen = C.mean(axis=0)
    w, v = np.linalg.eigh(np.cov((C - cen).T))
    s = (C - cen) @ v[:, 1]
    s = (s - s.min()) / (s.max() - s.min())
    ea, eb, md = s < 0.15, s > 0.85, (s >= 0.15) & (s <= 0.85)
    f = lambda m: (float(np.median(t[m])) if m.any() else float("nan"), int(m.sum()))
    d = {
        "z": z,
        "n": len(r),
        "t_min": float(t.min()),
        "t_med": float(np.median(t)),
        "t_max": float(t.max()),
        "p_min": int(npl.min()),
        "p_max": int(npl.max()),
        "nsets": len({x[0] for x in r}),
        "endA": f(ea),
        "endB": f(eb),
        "mid": f(md),
    }
    res.append(d)
    print(
        f"{z:5.2f} {d['n']:5d} {d['t_min'] * 1e3:6.2f} {d['t_med'] * 1e3:6.2f} {d['t_max'] * 1e3:7.2f} {d['p_min']:3d}-{d['p_max']:<3d} {d['nsets']:5d} | "
        f"{d['endA'][0] * 1e3:6.2f}/{d['endA'][1]:<4d}   {d['endB'][0] * 1e3:6.2f}/{d['endB'][1]:<4d} | {d['mid'][0] * 1e3:6.2f}"
    )
json.dump(res, open(out, "w"), indent=1)
