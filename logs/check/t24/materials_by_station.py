"""T24 (b): per-station material line-mass [kg/m] (rho t ds) mesh vs yaml, and TE_FLAT stiffness share.
usage: aeroenv.sh python materials_by_station.py <mesh.pkl>"""
import sys, pickle, numpy as np
sys.path.insert(0, "tools")
import audit_blade_input as ab
mesh = pickle.loads(open(sys.argv[1], "rb").read()); ab.add_set_zrange(mesh)
yb = ab.YamlBlade.from_file("tests/IEA-15-240-RWT.yaml")
zg = np.array(yb.outer["reference_axis"]["z"]["grid"]); zv = np.array(yb.outer["reference_axis"]["z"]["values"])
skin = [n for n in mesh["sets"] if "_SW" not in n and n in mesh["stacks"]]
print("z   perim_mid perim_yaml(open) | triax mesh yaml d% | carbon mesh yaml d% | flat_glass kg/m | EA_flat/EA_ring EIedge_flat/EIedge_ring (ring only, E1*t*ds about ring centroid)")
for z in (2, 5, 10, 15, 20, 30, 40, 60, 80, 100):
    z = float(z); f = float(np.interp(z, zv, zg)); per = yb.perimeter(f)
    at = [n for n in skin if mesh["set_zrange"][n][0] - 1e-9 <= z <= mesh["set_zrange"][n][1] + 1e-9]
    mm = {}; Lmid = 0; ea = 0; sx = 0; sxx = 0; eaf = 0; sxxf = 0; sxf = 0; gf = 0
    segs_all = []
    for n in at:
        for s in ab._slice(mesh, z, [n]):
            L = np.linalg.norm(s["pts"][0] - s["pts"][1]); xm = 0.5 * (s["pts"][0][0] + s["pts"][1][0])
            flat = n.endswith("TE_FLAT"); Lmid += 0 if flat else L
            for p in mesh["stacks"][n]:
                mm[p["mat"]] = mm.get(p["mat"], 0) + p["rho"] * p["t"] * L
                a = p["e1"] * p["t"] * L; ea += a; sx += a * xm
                segs_all.append((a, xm, flat))
                if flat and p["mat"] == "glass_uni": gf += p["rho"] * p["t"] * L
    xc = sx / ea
    EI = sum(a * (x - xc) ** 2 for a, x, fl in segs_all); EAr = ea
    # without the glass_uni of TE_FLAT: recompute excluding all TE_FLAT glass
    ea2 = 0; sx2 = 0; l2 = []
    for n in at:
        for s in ab._slice(mesh, z, [n]):
            L = np.linalg.norm(s["pts"][0] - s["pts"][1]); xm = 0.5 * (s["pts"][0][0] + s["pts"][1][0])
            for p in mesh["stacks"][n]:
                if n.endswith("TE_FLAT") and p["mat"] == "glass_uni": continue
                a = p["e1"] * p["t"] * L; ea2 += a; sx2 += a * xm; l2.append((a, xm))
    xc2 = sx2 / ea2; EI2 = sum(a * (x - xc2) ** 2 for a, x in l2)
    ym = {}
    for name in yb._skin_layers():
        st = yb.layer_state(name, f)
        if st["thickness"] <= 1e-12: continue
        s0, e0 = yb.layer_bounds(name, f, per)
        if s0 is None or e0 is None: continue
        ym[st["material"]] = ym.get(st["material"], 0) + (e0 - s0) * per * st["thickness"] * yb.material(st["material"])["rho"]
    d = lambda m: 100 * (mm.get(m, 0) - ym.get(m, 0)) / ym.get(m, 1)
    print(f"{z:5.0f} {Lmid:8.3f} {per:8.3f} | {mm.get('glass_triax',0):7.2f} {ym.get('glass_triax',0):7.2f} {d('glass_triax'):6.2f} | {mm.get('CarbonUD',0):7.2f} {ym.get('CarbonUD',0):7.2f} {d('CarbonUD'):6.2f} | {gf:6.2f} | dEA(glass flat)={100*(EAr-ea2)/EAr:5.2f}%  dEIedge={100*(EI-EI2)/EI:5.2f}%")
