import yaml,subprocess,numpy as np
P="/scratch/leahk/eduardo.donestevez/wt/input-fix-DECKY"
new=yaml.safe_load(open(P+"/tests/IEA-15-240-RWT.yaml"))
old=yaml.safe_load(subprocess.check_output(["git","-C",P,"show","ad5f3ac:tests/IEA-15-240-RWT.yaml"],text=True))
def layers(y):
    def walk(o):
        if isinstance(o,dict):
            if "name" in o and "thickness" in o and "material" in o: yield o
            for v in o.values(): yield from walk(v)
        elif isinstance(o,list):
            for v in o: yield from walk(v)
    return {l["name"]:l for l in walk(y)}
Lo,Ln=layers(old),layers(new)
g=np.linspace(0,1,50)
def f(l,key):
    d=l[key]; return np.interp(g,d["grid"],d["values"])
print("only old",sorted(set(Lo)-set(Ln)),"only new",sorted(set(Ln)-set(Lo)))
print("layers differing in thickness (max |dt| mm, at span nd, old-> new):")
for n in sorted(set(Lo)&set(Ln)):
    a,b=f(Lo[n],"thickness"),f(Ln[n],"thickness"); d=np.abs(a-b).max()
    if d>1e-6:
        i=np.abs(a-b).argmax(); print(f"  {n}: {d*1e3:.2f} mm at {g[i]:.2f} ({a[i]*1e3:.2f}->{b[i]*1e3:.2f}); mat {Lo[n]['material']}->{Ln[n]['material']}")
L=sorted(c["name"] for c in [])
for n in ("TE_reinforcement_SS","TE_reinforcement_PS"):
    a,b=f(Lo[n],"thickness"),f(Ln[n],"thickness")
    print(n,"material",Lo[n]["material"],Ln[n]["material"])
    for i in range(0,12): print(f"  r/R={g[i]:.3f} off={a[i]*1e3:.2f} mm deck={b[i]*1e3:.2f} mm")
    wa,wb=f(Lo[n],"width"),f(Ln[n],"width"); print("  width diff max",np.abs(wa-wb).max(), "grids old/new", len(Lo[n]['thickness']['grid']),len(Ln[n]['thickness']['grid']))
print("---- other differences")
def flat(o,p=""):
    if isinstance(o,dict):
        for k,v in o.items(): yield from flat(v,p+"/"+str(k))
    elif isinstance(o,list) and o and isinstance(o[0],(dict,list)):
        for i,v in enumerate(o): yield from flat(v,p+f"[{i}]")
    else: yield p,o
import math
A=dict(flat(old["components"]["blade"])); B=dict(flat(new["components"]["blade"]))
for k in sorted(set(A)|set(B)):
    a,b=A.get(k),B.get(k)
    try: same = a==b or np.allclose(np.interp(g,np.linspace(0,1,len(a)),a),np.interp(g,np.linspace(0,1,len(b)),b)) if isinstance(a,list) and isinstance(b,list) and a and b and not isinstance(a[0],str) else a==b
    except Exception: same=False
    if not same and "TE_reinforcement" not in k: print(" ",k[:110], str(a)[:50],"|",str(b)[:50])
print("materials equal:",old.get("materials")==new.get("materials"))
