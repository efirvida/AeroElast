import yaml,subprocess,numpy as np
P="/scratch/leahk/eduardo.donestevez/wt/input-fix-DECKY"
new=yaml.safe_load(open(P+"/tests/IEA-15-240-RWT.yaml"))
old=yaml.safe_load(subprocess.check_output(["git","-C",P,"show","ad5f3ac:tests/IEA-15-240-RWT.yaml"],text=True))
ob,nb=old["components"]["blade"]["outer_shape_bem"],new["components"]["blade"]["outer_shape_bem"]
for k in ("chord","twist","reference_axis","airfoil_position"):
    if k not in ob: continue
    if k=="reference_axis":
        for ax in "xyz":
            a,b=ob[k][ax],nb[k][ax]; g=np.array(b["grid"]); d=np.interp(g,a["grid"],a["values"])-np.array(b["values"]); print(k,ax,"max|d|",np.abs(d).max())
    elif "values" in ob[k]:
        a,b=ob[k],nb[k]; g=np.array(b["grid"]); d=np.interp(g,a["grid"],a["values"])-np.array(b["values"]); i=np.abs(d).argmax()
        print(k,"grid pts",len(a["grid"]),len(nb[k]["grid"]),"max|d|",np.abs(d).max(),"at",g[i])
    else: print(k,ob[k]==nb[k])
print("chord tail old",list(zip(np.round(ob["chord"]["grid"][-4:],3),np.round(ob["chord"]["values"][-4:],3))))
print("chord tail new",list(zip(np.round(nb["chord"]["grid"][-4:],3),np.round(nb["chord"]["values"][-4:],3))))
mo,mn={m["name"]:m for m in old["materials"]},{m["name"]:m for m in new["materials"]}
print("mat names diff",set(mo)^set(mn))
for n in mo:
    if n in mn and mo[n]!=mn[n]:
        print("mat",n,{k:(mo[n].get(k),mn[n].get(k)) for k in set(mo[n])|set(mn[n]) if mo[n].get(k)!=mn[n].get(k)})
