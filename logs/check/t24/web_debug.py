import sys, pickle, numpy as np
sys.path.insert(0, "tools")
import audit_blade_input as ab
mesh = pickle.loads(open(sys.argv[1], "rb").read()); ab.add_set_zrange(mesh)
X = mesh["nodes"]; els = mesh["elements"]
webs = [n for n in mesh["sets"] if "_SW" in n and n in mesh["stacks"]]
tot_area = 0; tot_m = 0
for n in webs:
    ar = 0
    for e in mesh["sets"][n]:
        P = X[els[e]]
        ar += 0.5*np.linalg.norm(np.cross(P[2]-P[0], P[3]-P[1])) if len(P)==4 else 0.5*np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
    rt = sum(p["rho"]*p["t"] for p in mesh["stacks"][n]); tot_m += ar*rt
print("n web sets", len(webs), "mass %.1f"%tot_m)
# one bay near z=50
for n in webs:
    zr = mesh["set_zrange"][n]
    if zr[0] <= 50 <= zr[1] and n.startswith("00"):
        P = X[np.unique(np.concatenate([els[e] for e in mesh["sets"][n]]))]
        print(n, "zrange", zr, "n_el", len(mesh["sets"][n]), "set mass/m(rho t)", sum(p["rho"]*p["t"] for p in mesh["stacks"][n]))
        # area vs h*dz
        ar = sum(0.5*np.linalg.norm(np.cross(X[els[e]][2]-X[els[e]][0], X[els[e]][3]-X[els[e]][1])) for e in mesh["sets"][n])
        L,_ = (lambda segs: (sum(np.linalg.norm(s["pts"][0]-s["pts"][1]) for s in segs), segs))(ab._slice(mesh, 0.5*(zr[0]+zr[1]), [n]))
        print("  area %.4f  cutlen_mid*dz %.4f  dz %.4f cutlen %.4f"%(ar, L*(zr[1]-zr[0]), zr[1]-zr[0], L))
