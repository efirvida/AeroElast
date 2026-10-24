"""T24 (a)(c): per-station TE_FLAT / TE_REINF region geometry and web height, mesh vs yaml.
usage: aeroenv.sh python stations.py <mesh.pkl> """
import sys, pickle, numpy as np
sys.path.insert(0, "tools")
import audit_blade_input as ab
mesh = pickle.loads(open(sys.argv[1], "rb").read()); ab.add_set_zrange(mesh)
yb = ab.YamlBlade.from_file("tests/IEA-15-240-RWT.yaml")
zg = np.array(yb.outer["reference_axis"]["z"]["grid"]); zv = np.array(yb.outer["reference_axis"]["z"]["values"])
def sets_at(z, suffix):
    return [n for n, (a, b) in mesh["set_zrange"].items() if n.endswith(suffix) and a - 1e-9 <= z <= b + 1e-9]
def seglen(z, names):
    segs = ab._slice(mesh, z, names); return sum(np.linalg.norm(s["pts"][0]-s["pts"][1]) for s in segs), segs
def gl(name):  # glass_uni total thickness in the set
    return sum(p["t"] for p in mesh["stacks"][name] if p["mat"] == "glass_uni")
for z in (20.0, 40.0, 60.0):
    f = float(np.interp(z, zv, zg)); per = yb.perimeter(f); o = yb.outer_shape(f)
    print(f"\n=== z={z} f={f:.3f} chord={o['chord']:.3f} perimeter(open polyline)={per:.3f} m")
    for lay in ("TE_reinforcement_SS", "TE_reinforcement_PS"):
        st = yb.layer_state(lay, f); s, e = yb.layer_bounds(lay, f, per)
        print(f"yaml {lay}: t={st['thickness']*1e3:.2f} mm width={st.get('width')} nd_arc=[{s:.4f},{e:.4f}] -> {(e-s)*per:.3f} m of arc")
    for side in ("HP", "LP"):
        for reg in ("TE_FLAT", "TE_REINF", "TE_PANEL"):
            for n in sets_at(z, f"{side}_{reg}"):
                L, segs = seglen(z, [n]);
                if not segs: continue
                pts = np.array([p for s in segs for p in s["pts"]])
                print(f"mesh {n:22s} cut-length={L:.3f} m  glass_uni t={gl(n)*1e3:.2f} mm  x[{pts[:,0].min():.3f},{pts[:,0].max():.3f}] y[{pts[:,1].min():.3f},{pts[:,1].max():.3f}]")
    # webs
    names = [n for n in mesh["sets"] if "_SW" in n and mesh["set_zrange"].get(n, (1, 0))[0] - 1e-9 <= z <= mesh["set_zrange"][n][1] + 1e-9]
    poly, _, _ = yb.airfoil_at(f)
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]; seg /= seg[-1]
    groups = {}
    for n in names: groups.setdefault(n.split("_")[0], []).append(n)
    for k in range(len(yb.struct["webs"])):
        w = yb.web_state(k, f)
        p = [np.array([np.interp(a, seg, poly[:, i]) for i in (0, 1)]) for a in (w["start"], w["end"])]
        h = np.linalg.norm(p[0]-p[1]) * o["chord"]
        print(f"yaml web{k}: outer-point height {h:.4f} m")
    for g, ns in sorted(groups.items()):
        L, segs = seglen(z, ns); print(f"mesh web group {g}: height {L:.4f} m")
    for side in ("HP", "LP"):
        for n in sets_at(z, f"{side}_SPAR"):
            if ab._slice(mesh, z, [n]): print(f"  skin total t at spar {n}: {sum(p['t'] for p in mesh['stacks'][n])*1e3:.1f} mm")
