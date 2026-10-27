import pickle
for n in ("DECKY","T14m"):
    d=pickle.load(open(f"logs/check/input-fix/capture-{n}.pkl","rb"))
    k=[k for k in d["fixtures"] if k.startswith("s2_results")][0]
    for c,v in d["fixtures"][k].items(): print(n,c,{a:(b if not hasattr(b,'__len__') or isinstance(b,str) else '<arr>') for a,b in v.items()} if isinstance(v,dict) else v)
