import pickle
from _aeroelast import Laminate
print([a for a in dir(Laminate) if not a.startswith("_")])
from aeroelast.models.blade.model import build_rust_properties
d=pickle.load(open("logs/check/gj/mesh_head.pkl","rb"))
secs=[s for s in d["numad_sections"] if s["elementSet"].startswith("03_412") or s["elementSet"].startswith("03_41")][:2]
p=build_rust_properties({"materials":d["materials"],"sections":secs})
for k,v in p.items():
    print(k, secs[0]["layup"][:3]); print({a:(getattr(v,a)() if callable(getattr(v,a)) else getattr(v,a)) for a in dir(v) if a in ("abd","a_matrix","Cm","thickness","total_thickness")})
