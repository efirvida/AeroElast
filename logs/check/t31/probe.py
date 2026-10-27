import sys, time, numpy as np
sys.path.insert(0, "."); 
from tests.test_blade_ccx_parity import _blade, SPAN_DIRECTION, BLADE_YAML
from aeroelast.core.mesh.io.writers import _build_angle_bucket_sets
t=time.time()
mesh, props, coords = _blade()
print("mesh", len(mesh.nodes), len(mesh.elements), "SPAN", SPAN_DIRECTION, round(time.time()-t,1),"s")
sd=np.array(SPAN_DIRECTION,float); sd/=np.linalg.norm(sd)
idx=mesh.node_id_to_index
def frame(p):
    v1=p[1]-p[0]; v2=p[2]-p[0]
    n1=np.cross(v1,v2); n1/=np.linalg.norm(n1)
    if len(p)>=4:
        n2=np.cross(p[2]-p[0],p[3]-p[0]); n2/=np.linalg.norm(n2)
        e3=n1+n2; e3/=np.linalg.norm(e3)
    else: e3=n1
    e1=v1-v1@e3*e3; e1/=np.linalg.norm(e1); e2=np.cross(e3,e1); return e1,e2,e3
rows=[]
for el in mesh.elements:
    p=coords[[idx[n] for n in el.node_ids[:4]]]
    if len(p)<3: continue
    e1,e2,e3=frame(p)
    sp=sd-sd@e3*e3; nrm=np.linalg.norm(sp)
    if nrm<1e-10: rows.append((p.mean(0)[2],np.nan,np.nan,0)); continue
    off=np.degrees(np.arctan2(sp@e2/nrm,sp@e1/nrm))
    # AE: bucketed (0.1 deg) offset, fibre = cos e1 + sin e2
    offb=round(off*10)/10
    a=np.radians(offb); fib_ae=np.cos(a)*e1+np.sin(a)*e2
    fib_ccx=sp/nrm   # CCX: 1-axis=span projected on shell surface
    ang=np.degrees(np.arccos(np.clip(fib_ae@fib_ccx,-1,1)))
    # tilt of the fibre out of global span axis
    tilt=np.degrees(np.arccos(np.clip(abs(fib_ccx@sd),-1,1)))
    rows.append((p.mean(0)[2],ang,tilt,abs(e3@sd)))
R=np.array(rows); ok=~np.isnan(R[:,1])
print("elements",len(R),"degenerate(span||normal)",(~ok).sum())
a=R[ok,1]; print("AE-vs-CCX fibre angle deg: median %.4f p95 %.4f max %.4f"%(np.median(a),np.percentile(a,95),a.max()))
t=R[ok,2]; print("fibre tilt from global z (deg): median %.2f p95 %.2f max %.2f"%(np.median(t),np.percentile(t,95),t.max()))
z=R[ok,0]
for lo,hi in [(0,5),(5,20),(20,60),(60,100),(100,118),(118,121)]:
    m=(z>=lo)&(z<hi)
    if m.sum(): print(f"z {lo:>3}-{hi:<3} n={m.sum():5d} AE-CCX max {a[m].max():.4f}  tilt p95 {np.percentile(t[m],95):.2f}")
# ply angle sets + sign sensitivity
angs={}
for k,pr in props.items():
    try: angs[k]=sorted({round(p['angle'],1) for p in pr.plies_data()})
    except Exception as e: pass
allang=sorted({a for v in angs.values() for a in v}); print("ply angles present:",allang)
