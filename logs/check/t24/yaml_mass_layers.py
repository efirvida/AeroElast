"""Yaml mass per 1 m z-bin from tools/audit_blade_input.py's YamlBlade (outer-surface arcs, no mesh)."""
import sys, json, numpy as np
sys.path.insert(0, "/scratch/leahk/eduardo.donestevez/fem-shell/tools")
from audit_blade_input import YamlBlade
yb = YamlBlade.from_file("/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA-15-240-RWT.yaml")
zg = np.array(yb.outer["reference_axis"]["z"]["grid"]); zv = np.array(yb.outer["reference_axis"]["z"]["values"])
L = yb.length; dz = 0.25
skin = {}; web = {}; webs_tot = np.zeros(125); skin_tot = np.zeros(125)
for z in np.arange(dz/2, L, dz):
    f = float(np.interp(z, zv, zg)); b = min(int(z), 124)
    per = yb.perimeter(f); o = yb.outer_shape(f)
    for name in yb._skin_layers():
        st = yb.layer_state(name, f)
        if st["thickness"] <= 1e-12: continue
        s, e = yb.layer_bounds(name, f, per)
        if s is None or e is None: continue
        m = (e - s) * per * st["thickness"] * yb.material(st["material"])["rho"] * dz
        skin.setdefault(name if st['material']=='glass_uni' else st["material"], np.zeros(125))[b] += m; skin_tot[b] += m
    poly, _, _ = yb.airfoil_at(f)
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]; seg /= seg[-1]
    for k in range(len(yb.struct["webs"])):
        w = yb.web_state(k, f)
        p = [np.array([np.interp(a, seg, poly[:, i]) for i in (0, 1)]) for a in (w["start"], w["end"])]
        h = np.linalg.norm(p[0] - p[1]) * o["chord"]
        for ly in yb.web_layers(k, f):
            m = h * ly["thickness"] * yb.material(ly["material"])["rho"] * dz
            web.setdefault(ly["material"], np.zeros(125))[b] += m; webs_tot[b] += m
json.dump({"skin_bins": skin_tot.tolist(), "web_bins": webs_tot.tolist(), "skin": skin_tot.sum(), "web": webs_tot.sum(),
           "skin_mat": {k: v.tolist() for k, v in skin.items()}, "web_mat": {k: v.tolist() for k, v in web.items()}}, open("bins-yaml-layers.json", "w"))
print("yaml skin %.1f web %.1f" % (skin_tot.sum(), webs_tot.sum()))
for k, v in skin.items(): print("skin", k, round(v.sum(), 1))
for k, v in web.items(): print("web", k, round(v.sum(), 1))
