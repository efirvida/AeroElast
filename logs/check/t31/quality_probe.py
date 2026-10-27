import sys, re, numpy as np
tree, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, tree)
import tests.test_blade_ccx_parity as T
mesh, props, coords = T._blade()
idx = mesh.node_id_to_index
# isotropic dict sets: AE vs deck
txt = open(out + "/d.inp").read()
iso = {n: p for n, p in props.items() if isinstance(p, dict)}
n0 = next(iter(iso)); print("iso example", n0, {k: v for k, v in iso[n0].items()})
m = re.search(r"\*SHELL SECTION, ELSET=E%s,[^\n]*\n([^\n]*)" % n0.upper(), txt); print("deck section:", m.group(0) if m else None)
mm = re.search(r"\*MATERIAL, NAME=MAT_%s\n\*ELASTIC\n([^\n]*)" % n0, txt); print("deck mat:", mm.group(1) if mm else None)
# geometry quality
L = []; asp = []; area = []; zc = []
for e in mesh.elements:
    p = coords[[idx[n] for n in e.node_ids[:4]]]
    if len(p) < 4: continue
    ed = [np.linalg.norm(p[(i+1)%4]-p[i]) for i in range(4)]
    asp.append(max(ed)/min(ed)); area.append(0.5*np.linalg.norm(np.cross(p[2]-p[0], p[3]-p[1]))); zc.append(p[:,2].mean())
    # warp: distance of p3 from plane p0p1p2
    n = np.cross(p[1]-p[0], p[2]-p[0]); n /= max(np.linalg.norm(n),1e-30); L.append(abs((p[3]-p[0])@n)/max(ed))
asp=np.array(asp); area=np.array(area); zc=np.array(zc); W=np.array(L)
print("quads %d  aspect median %.2f p99 %.1f max %.1f  n(asp>10) %d n(asp>50) %d" % (len(asp), np.median(asp), np.percentile(asp,99), asp.max(), (asp>10).sum(), (asp>50).sum()))
print("area min %.3e n<1e-4 %d  warp max %.3f n(warp>0.05) %d" % (area.min(), (area<1e-4).sum(), W.max(), (W>0.05).sum()))
for lo,hi in [(0,20),(20,60),(60,100),(100,117.01)]:
    k=(zc>=lo)&(zc<hi); print("  z %3d-%3d n=%5d aspect p99 %.1f  n(asp>10)=%d" % (lo,hi,k.sum(),np.percentile(asp[k],99),(asp[k]>10).sum()))
# spanwise rows
zs=np.unique(np.round(coords[:,2],4)); dz=np.diff(zs); print("distinct z stations %d  dz min %.4f median %.4f max %.4f ; n(dz<0.01)=%d" % (len(zs), dz.min(), np.median(dz), dz.max(), (dz<0.01).sum()))
