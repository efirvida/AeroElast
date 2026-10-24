import pickle, numpy as np
d=pickle.load(open("logs/check/input-fix/capture-T15m.pkl","rb"))["fixtures"]
k=[k for k in d if k.startswith("s1_results")][0]; s=d[k]
z=np.asarray(s["z"]); gj=np.asarray(s["gj"]); bd=np.asarray(s["gj_bd"])
for zi,g,b in zip(z,gj,bd): print(f"r/R={zi/z.max():.3f} FE={g:.4e} BD={b:.4e} rel={g/b-1:+.1%}")
