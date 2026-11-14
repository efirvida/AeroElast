"""T3: SONATA/ANBA4 on OUR v2 yaml (fem-shell tests/IEA15MW/reference/upstream, sha af51499e) at deck stations."""
import os, sys, json, hashlib
import numpy as np
from SONATA.classBlade import Blade
yaml_path = sys.argv[1]; out_dir = sys.argv[2]; os.makedirs(out_dir, exist_ok=True)
stations = [float(s) for s in sys.argv[3].split(",")]
sha = hashlib.sha256(open(yaml_path, "rb").read()).hexdigest()
flags = {"flag_wt_ontology": True, "flag_ref_axes_wt": True, "attribute_str": "MatID",
         "flag_plotDisplacement": False, "flag_plotTheta11": False, "flag_wf": False, "flag_lft": False,
         "flag_topo": True, "mesh_resolution": 400, "flag_recovery": False, "c2_axis": False}
job = Blade(name="IEA15_t3", filename=yaml_path, flags=flags, stations=stations)
job.blade_gen_section(topo_flag=True, mesh_flag=True)
job.blade_run_anbax()
res = {"yaml": yaml_path, "sha256": sha, "stations": stations, "TS": [], "MM": []}
for i, x in enumerate(stations):
    bp = job.beam_properties[i, 1]
    res["TS"].append(np.asarray(bp.TS).tolist()); res["MM"].append(np.asarray(bp.MM).tolist())
    print(f"eta {x:.3f}  TS diag (SONATA frame) =", " ".join(f"{v:.4e}" for v in np.diag(bp.TS)))
json.dump(res, open(os.path.join(out_dir, "sonata_t3.json"), "w"), indent=1)
try:
    job.blade_plot_sections(attribute="MatID", plotTheta11=False, plotDisplacement=False, savepath=out_dir + os.sep)
except Exception as e:
    print("plot failed:", e)
