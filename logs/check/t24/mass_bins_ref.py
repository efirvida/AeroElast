"""Assembled blade mass per 1 m z-bin, skin / webs, at es 0.25 (TestBladeMassIntegration path)."""
import json, sys
import numpy as np
from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.elements import ElementFamily
from aeroelast.models.blade.model import build_rust_properties

yaml_file, out = sys.argv[1], sys.argv[2]
g = BladeMesh(yaml_file=yaml_file, element_size=0.25, skin_reference=sys.argv[3])
mesh = g.generate(renumber="rcm")
props = build_rust_properties(g.numad_mesh_data)
asm = MeshAssembler(mesh, {"elements": {"element_family": ElementFamily.SHELL,
                                        "span_direction": (0.0, 0.0, 1.0), "properties": props}})
M = asm.assemble_mass_matrix(); M.destroy()
r, c, v = asm._rust.assemble_m()
total = float(np.asarray(v)[np.asarray(r) % 6 == 0].sum())

# per-set areal density rho*h from the NuMAD layup (same data the Rust laminates are built from)
mats = {m["name"]: float(m["density"]) for m in g.numad_mesh_data["materials"]}
arho = {s["elementSet"]: sum(mats[mn] * t for mn, t, a in s["layup"]) for s in g.numad_mesh_data["sections"]}
layup_of = {s["elementSet"]: s["layup"] for s in g.numad_mesh_data["sections"]}
nid = mesh.node_id_to_index
X = mesh.coords_array
bins = np.zeros(125); cls = {"skin": np.zeros(125), "web": np.zeros(125)}
setmass = {}
matbins = {}
for name, es in mesh.element_sets.items():
    if name not in arho: continue
    kind = "web" if ("_SW" in name or "web" in name.lower()) else "skin"
    for e in es.elements:
        P = X[[nid[n.id] for n in e.nodes]]
        if len(P) == 4:
            A = 0.5*np.linalg.norm(np.cross(P[2]-P[0], P[3]-P[1]))
        else:
            A = 0.5*np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
        m = A*arho[name]; b = min(int(P[:, 2].mean()), 124)
        for (mn, t, a) in layup_of[name]:
            if mn == 'glass_uni': matbins.setdefault('guni@' + name.split('_', 2)[2], np.zeros(125))[b] += A*mats[mn]*t
            matbins.setdefault(kind + ":" + mn, np.zeros(125))[b] += A*mats[mn]*t
        cls[kind][b] += m; setmass[name] = setmass.get(name, 0.0) + m
res = {"M_total": total, "skin_bins": cls["skin"].tolist(), "web_bins": cls["web"].tolist(),
       "skin": cls["skin"].sum(), "web": cls["web"].sum(), "sets": setmass,
       "matbins": {k: v.tolist() for k, v in matbins.items()}, "zmax": float(X[:, 2].max())}
json.dump(res, open(out, "w"))
print(out, "M=%.1f skin=%.1f web=%.1f zmax=%.3f" % (total, res["skin"], res["web"], res["zmax"]))
print(sorted(setmass)[:8], len(setmass))
