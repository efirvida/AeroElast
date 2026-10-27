"""Independent reading of the windIO yaml: layers covering a nd_arc point at span eta -> ply stack."""
import numpy as np, yaml
Y = "/scratch/leahk/eduardo.donestevez/fem-shell/tests/IEA-15-240-RWT.yaml"
_d = yaml.safe_load(open(Y))
B = _d["components"]["blade"]["internal_structure_2d_fem"]
MATS = {m["name"]: m for m in _d["materials"]}
def at(spec, eta):
    return float(np.interp(eta, spec["grid"], spec["values"]))
def layers(eta):
    out = []
    for L in B["layers"]:
        if "thickness" not in L: continue
        t = at(L["thickness"], eta)
        fo = at(L["fiber_orientation"], eta) if "fiber_orientation" in L else 0.0
        ex = {k: (at(L[k], eta) if isinstance(L.get(k), dict) and "values" in L[k] else L.get(k))
              for k in ("start_nd_arc","end_nd_arc","midpoint_nd_arc","width","offset_y_pa")}
        fx = {k: L[k].get("fixed") for k in ("start_nd_arc","end_nd_arc","midpoint_nd_arc") if isinstance(L.get(k), dict) and "fixed" in L[k]}
        s0, s1 = ex["start_nd_arc"], ex["end_nd_arc"]
        out.append(dict(ex=ex, fx=fx, name=L["name"], mat=L["material"], t=t, ang=fo, s0=s0, s1=s1, web=L.get("web")))
    return out
if __name__ == "__main__":
    eta = 0.7002
    for L in layers(eta):
        print(f"{L['name']:22s} {L['mat']:18s} t={L['t']*1e3:8.3f} mm ang={L['ang']:+.3f} ex={ {k:(round(v,4) if isinstance(v,float) else v) for k,v in L['ex'].items() if v is not None} } fx={L['fx']} web={L['web']}")
