"""W1: export SONATA's section mesh (cells: MatID, area, theta_1, theta_3, node xy) at one eta. Runs in the sonata or sonata2021 env.
  python tools/l3_w1_sonata_export.py YAML ETA OUT.json [mesh_resolution]"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
from SONATA.classBlade import Blade  # noqa: E402

yaml, eta, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
res = int(sys.argv[4]) if len(sys.argv) > 4 else 400
flags = {"flag_wt_ontology": True, "flag_ref_axes_wt": True, "attribute_str": "MatID", "flag_plotDisplacement": False,
         "flag_plotTheta11": False, "flag_wf": False, "flag_lft": False, "flag_topo": True, "mesh_resolution": res,
         "flag_recovery": False}
job = Blade(name="w1", filename=yaml, flags=flags, stations=[eta])
job.blade_gen_section(topo_flag=True, mesh_flag=True, split_quads=True)
sec = job.sections[0][1]


def xy(n):
    for a in ("coordinates", "xy"):
        if hasattr(n, a):
            return [float(v) for v in getattr(n, a)[:2]]
    return [float(n.x), float(n.y)]


cells = [{"mat": int(c.MatID), "area": float(c.area), "th1": float(c.theta_1[0]), "th3": None if c.theta_3 is None else float(c.theta_3),
          "xy": [xy(n) for n in c.nodes]} for c in sec.mesh]
names = {int(k): getattr(m, "name", str(k)) for k, m in job.materials.items()}
json.dump({"yaml": yaml, "eta": eta, "res": res, "names": names, "cells": cells}, open(out, "w"))
print("cells", len(cells), "total area", sum(c["area"] for c in cells))
