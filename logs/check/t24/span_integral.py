"""T24 (b)(c): line mass rho*t*ds [kg/m] of triax, carbon, web (biax+foam) at z = 0.5,1.5,... mesh (planar cut)
vs yaml (outer arcs), then integrated over span. usage: aeroenv.sh python span_integral.py <mesh.pkl> <out.csv>"""
import sys, pickle, numpy as np
sys.path.insert(0, "tools")
import audit_blade_input as ab
mesh = pickle.loads(open(sys.argv[1], "rb").read()); ab.add_set_zrange(mesh)
yb = ab.YamlBlade.from_file("tests/IEA-15-240-RWT.yaml")
zg = np.array(yb.outer["reference_axis"]["z"]["grid"]); zv = np.array(yb.outer["reference_axis"]["z"]["values"])
skin = [n for n in mesh["sets"] if "_SW" not in n and n in mesh["stacks"]]
webs = [n for n in mesh["sets"] if "_SW" in n and n in mesh["stacks"]]
zr = mesh["set_zrange"]
rows = []
for z in np.arange(0.5, 117.0, 1.0):
    f = float(np.interp(z, zv, zg)); per = yb.perimeter(f); o = yb.outer_shape(f)
    mm = {"glass_triax": 0.0, "CarbonUD": 0.0, "web": 0.0}
    for n in skin:
        if not (zr[n][0] <= z <= zr[n][1]): continue
        for s in ab._slice(mesh, z, [n]):
            L = np.linalg.norm(s["pts"][0] - s["pts"][1])
            for p in mesh["stacks"][n]:
                if p["mat"] in mm: mm[p["mat"]] += p["rho"] * p["t"] * L
    hmesh = 0.0
    for n in webs:
        if not (zr[n][0] <= z <= zr[n][1]): continue
        for s in ab._slice(mesh, z, [n]):
            L = np.linalg.norm(s["pts"][0] - s["pts"][1]); hmesh += L
            mm["web"] += sum(p["rho"] * p["t"] for p in mesh["stacks"][n]) * L
    ym = {"glass_triax": 0.0, "CarbonUD": 0.0, "web": 0.0}
    for name in yb._skin_layers():
        st = yb.layer_state(name, f)
        if st["thickness"] <= 1e-12 or st["material"] not in ym: continue
        s0, e0 = yb.layer_bounds(name, f, per)
        if s0 is None or e0 is None: continue
        ym[st["material"]] += (e0 - s0) * per * st["thickness"] * yb.material(st["material"])["rho"]
    poly, _, _ = yb.airfoil_at(f)
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]; seg /= seg[-1]
    hy = 0.0
    for k in range(len(yb.struct["webs"])):
        w = yb.web_state(k, f)
        p = [np.array([np.interp(a, seg, poly[:, i]) for i in (0, 1)]) for a in (w["start"], w["end"])]
        h = np.linalg.norm(p[0] - p[1]) * o["chord"]; hy += h
        for ly in yb.web_layers(k, f): ym["web"] += h * ly["thickness"] * yb.material(ly["material"])["rho"]
    rows.append((z, mm["glass_triax"], ym["glass_triax"], mm["CarbonUD"], ym["CarbonUD"], mm["web"], ym["web"], hmesh, hy))
a = np.array(rows); np.savetxt(sys.argv[2], a, delimiter=",", header="z,triax_mesh,triax_yaml,carbon_mesh,carbon_yaml,web_mesh,web_yaml,hweb_mesh,hweb_yaml")
print("span integral (1 m midpoint rule, kg): mesh / yaml / diff / %")
for i, nm in enumerate(("triax", "carbon", "web")):
    m, y = a[:, 1+2*i].sum(), a[:, 2+2*i].sum(); print(f"{nm:7s} {m:9.1f} {y:9.1f} {m-y:8.1f} {100*(m-y)/y:6.2f}")
print("web height sum mesh/yaml:", a[:, 7].sum(), a[:, 8].sum())
