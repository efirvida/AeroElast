#!/usr/bin/env python
"""L3 / T9 (W2): dump SONATA's section mesh (cells: MatID, area, centre) and its layer table, for cap width / thickness analysis.
Works with SONATA 2026 (env conda-envs/sonata) and 2021 (env $SCRATCH/tmp/t3b/bin/env.sh). args: yaml out.npz eta [mesh_resolution]
Layer table = segment.LayerLst (ID, S1, S2, thickness, MatID, name) as SONATA built it."""

import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
from SONATA.classBlade import Blade

yml, out, eta = sys.argv[1], sys.argv[2], float(sys.argv[3])
res = int(sys.argv[4]) if len(sys.argv) > 4 else 400
flags = {
    "flag_wt_ontology": True,
    "flag_ref_axes_wt": True,
    "attribute_str": "MatID",
    "flag_plotDisplacement": False,
    "flag_plotTheta11": False,
    "flag_wf": False,
    "flag_lft": False,
    "flag_topo": True,
    "mesh_resolution": res,
    "flag_recovery": False,
}
job = Blade(name="w2", filename=yml, flags=flags, stations=[eta])
job.blade_gen_section(topo_flag=True, mesh_flag=True, split_quads=True)
sec = job.sections[0][1]
names = {k: getattr(m, "name", str(k)) for k, m in job.materials.items()}
mat = np.array([c.MatID for c in sec.mesh])
area = np.array([c.area for c in sec.mesh])
ctr = np.array([np.mean([[n.Pnt2d.X(), n.Pnt2d.Y()] for n in c.nodes], axis=0) for c in sec.mesh])
lay = []
for si, seg in enumerate(sec.SegmentLst):
    for L in seg.LayerLst:
        lay.append(
            (si, int(L.ID), float(L.S1), float(L.S2), float(L.thickness), int(L.MatID), str(L.name))
        )
np.savez(
    out,
    mat=mat,
    area=area,
    ctr=ctr,
    names=np.array([f"{k}:{v}" for k, v in names.items()]),
    lay=np.array(lay, dtype=object),
    allow_pickle=True,
)
for k in sorted(set(mat)):
    print("MatID", k, names.get(k), "area %.5f" % area[mat == k].sum())
print("total %.5f" % area.sum())
for l in lay:
    print("LAYER seg%d id%d S=[%.4f,%.4f] t=%.5f mat%d %s" % l)
