"""W1: SONATA 2026 on ONE section with one factor changed after meshing (the cells are edited, then ANBA4 runs on the edited mesh).
  python tools/l3_w1_sonata_variants.py YAML ETA OUT.json [mesh_resolution]
variants: base | noadh (adhesive cells removed) | adh_x10 / adh_x0.1 (adhesive E and G scaled via a material copy is NOT done; see note)
Only structural edits of the mesh are applied (cell deletion / MatID remap), never a change of ply data."""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
from SONATA.classBlade import Blade  # noqa: E402

yaml, eta, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
res = int(sys.argv[4]) if len(sys.argv) > 4 else 400
flags = {"flag_wt_ontology": True, "flag_ref_axes_wt": True, "attribute_str": "MatID", "flag_plotDisplacement": False,
         "flag_plotTheta11": False, "flag_wf": False, "flag_lft": False, "flag_topo": True, "mesh_resolution": res,
         "flag_recovery": False}
res_out = {"yaml": yaml, "eta": eta, "res": res, "variants": {}}
names = None
for var in ("base", "noadh", "adh_as_foam"):
    job = Blade(name="w1", filename=yaml, flags=flags, stations=[eta])
    job.blade_gen_section(topo_flag=True, mesh_flag=True, split_quads=True)
    cs = job.sections[0][1]
    ids = {getattr(m, "name", k): k for k, m in job.materials.items()}
    adh, foam = ids["adhesive"], ids["medium_density_foam"]
    n0 = len(cs.mesh)
    if var == "noadh":
        cs.mesh = [c for c in cs.mesh if c.MatID != adh]
    elif var == "adh_as_foam":
        for c in cs.mesh:
            if c.MatID == adh:
                c.MatID = foam
    job.blade_run_anbax()
    bp = job.beam_properties[0, 1]
    res_out["variants"][var] = {"TS": np.asarray(bp.TS).tolist(), "cells": n0, "cells_after": len(cs.mesh)}
    print("VARIANT", var, "cells", n0, "->", len(cs.mesh), "TS diag", " ".join(f"{v:.4e}" for v in np.diag(bp.TS)), flush=True)
json.dump(res_out, open(out, "w"), indent=1)
