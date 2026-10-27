import pickle
BD={"LC1_gravity":(1,-2.231483),"LC2_tip_flap":(1,82.273278),"LC3_tip_edge":(0,34.469441),"LC4_distributed":(1,9.062631)}
R={}
for n in ("DECKY","T14m"):
    d=pickle.load(open(f"logs/check/input-fix/capture-{n}.pkl","rb"))
    k=[k for k in d["fixtures"] if k.startswith("s2_results")][0]
    R[n]=d["fixtures"][k]
print("case comp shell_DECKY shell_T14m BD err%_DECKY err%_T14m d_abs d_pp ring_DECKY ring_T14m")
for c,(i,bd) in BD.items():
    a=R["DECKY"][c]["shell"][i]; b=R["T14m"][c]["shell"][i]
    ea=(a/bd-1)*100; eb=(b/bd-1)*100
    print(f"{c} {'xyz'[i]} {a:.6f} {b:.6f} {bd:.6f} {ea:+.3f} {eb:+.3f} {a-b:+.6f} {ea-eb:+.3f} {R['DECKY'][c]['ring'][i]:.6f} {R['T14m'][c]['ring'][i]:.6f}")
