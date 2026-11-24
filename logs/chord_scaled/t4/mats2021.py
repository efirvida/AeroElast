import sys, collections, matplotlib; matplotlib.use("Agg")
from SONATA.classBlade import Blade
Y = sys.argv[1]
flags = {"flag_wt_ontology": True, "flag_ref_axes_wt": True, "attribute_str": "MatID", "flag_plotDisplacement": False,
         "flag_plotTheta11": False, "flag_wf": False, "flag_lft": False, "flag_topo": True, "mesh_resolution": 400, "flag_recovery": False}
for x in map(float, sys.argv[2].split(",")):
    job = Blade(name="m", filename=Y, flags=flags, stations=[x])
    job.blade_gen_section(topo_flag=True, mesh_flag=True, split_quads=True)
    sec = job.sections[0][1]; A = collections.defaultdict(float)
    for c in sec.mesh: A[c.MatID] += c.area
    names = {k: getattr(m, "name", k) for k, m in job.materials.items()}
    print("eta", x)
    for k in sorted(A): print("  MatID", k, names.get(k), "area %.5e" % A[k])
    print("  total area %.5e" % sum(A.values()), flush=True)
