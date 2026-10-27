import pickle, numpy as np
d=pickle.load(open("logs/check/gj/mesh_head.pkl","rb"))
zs=d["zs"]; print(zs.max(), len(zs))
i=np.argmin(abs(zs/zs.max()-0.70)); print(i, zs[i]/zs.max())
secs=d["numad_sections"]; print(len(secs)); 
for s in secs[:3]: print({k:(v if k!="layup" else v[:3]) for k,v in s.items()})
print([m["name"] for m in d["materials"]]); print(d["materials"][0])
print(list(d["esets"])[:60])
