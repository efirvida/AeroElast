"""Weld nodes that sit at the same position (default tolerance 1e-9 m) in a cached v2 blade mesh (.npz of tests/support/v2_blade or
logs/check/gmsh_fem_webs/mesh): duplicates are replaced by the first node of their group, the node array is compacted and every
element index remapped. Element arrays (surfaces, params, layup_id, normals, offsets, areas) are untouched.
usage: python tools/weld_coincident_nodes.py IN.npz OUT.npz [tol]    prints the groups welded and the surfaces involved.
This is a measurement tool for the web0/web1 overlap of the IEA-15 gmsh mesh (odd/tasks/gmsh-gap-reduction.md); it does not change the mesher.
"""
import sys

import numpy as np

src, dst = sys.argv[1], sys.argv[2]
tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-9
z = dict(np.load(src, allow_pickle=False))
x = np.asarray(z["nodes"], float)
el = np.asarray(z["elements"])
key = np.round(x / tol).astype(np.int64)
_, first, inv, cnt = np.unique(key, axis=0, return_index=True, return_inverse=True, return_counts=True)
inv = np.ravel(inv)
rep = first[inv]                      # representative (lowest-index first occurrence) of every node
keep = np.unique(rep)
new_id = np.full(len(x), -1, np.int64)
new_id[keep] = np.arange(len(keep))
mapped = np.where(el >= 0, new_id[rep[np.where(el >= 0, el, 0)]], -1)
surf = np.asarray(z["surfaces"])
groups = np.nonzero(cnt > 1)[0]
pairs = set()
for g in groups:
    members = np.nonzero(inv == g)[0]
    own = set()
    for k, e in enumerate(el):
        if set(e[e >= 0]) & set(members.tolist()):
            own.add(str(surf[k]))
    pairs.add(tuple(sorted(own)))
print(f"welded {int((cnt - 1).sum())} nodes in {len(groups)} groups; surfaces {sorted(pairs)}; nodes {len(x)} -> {len(keep)}")
bad = [k for k, e in enumerate(mapped) if len(set(e[e >= 0])) != len(e[e >= 0])]
print("elements with a repeated node after welding:", len(bad))
z["nodes"] = x[keep]
z["elements"] = mapped
np.savez(dst, **z)
