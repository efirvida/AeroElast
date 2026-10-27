"""T24 (c)/(b): per set, sum(rho t) from the NuMAD section layup (what mass_bins.py used) vs from the
laminate handed to the solver (Blade.get_element_properties().plies_data(), what the audit reads).
usage: aeroenv.sh python layup_vs_laminate.py <out.json>"""
import sys, json, numpy as np
from aeroelast.models.blade.model import Blade, build_rust_properties
b = Blade("tests/IEA-15-240-RWT.yaml", element_size=0.25); b.generate_mesh()
mats = {m["name"]: float(m["density"]) for m in b._numad_mesh["materials"]}
layups = {s["elementSet"]: s["layup"] for s in b._numad_mesh["sections"]}
props = b.get_element_properties()
rprops = build_rust_properties(b._numad_mesh)
print("same object as build_rust_properties?", type(props).__name__, type(rprops).__name__, len(props), len(rprops))
out = {}
for n, lam in props.items():
    lay = layups.get(n)
    if lay is None or isinstance(lam, dict): continue
    pd = lam.plies_data()
    a = sum(mats[m] * t for m, t, _ in lay); c = sum(p["material"]["rho"] * p["thickness"] for p in pd)
    out[n] = (a, c, [(m, t) for m, t, _ in lay], [(p["material"].get("name"), p["thickness"], p["material"]["rho"]) for p in pd])
json.dump(out, open(sys.argv[1], "w"))
sw = {n: v for n, v in out.items() if "_SW" in n}
print("web sets", len(sw), "sum arho_layup %.2f sum arho_laminate %.2f" % (sum(v[0] for v in sw.values()), sum(v[1] for v in sw.values())))
k = next(iter(sw)); print(k, sw[k][2], sw[k][3])
sk = {n: v for n, v in out.items() if "_SW" not in n}
print("skin sets", len(sk), "ratio laminate/layup: min %.4f max %.4f" % (min(v[1]/v[0] for v in sk.values() if v[0] > 0), max(v[1]/v[0] for v in sk.values() if v[0] > 0)))
