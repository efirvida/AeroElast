"""Re-mesh like mk_mesh.py and dump numad key arcs (region extents) + span stations. argv: out.pkl"""
import sys, pickle, numpy as np
from aeroelast.core.mesh.generators import BladeMesh
W = "/scratch/leahk/eduardo.donestevez/fem-shell/"
gen = BladeMesh(yaml_file=W + "tests/IEA-15-240-RWT.yaml", element_size=0.25, skin_reference="mid")
gen.generate(renumber="rcm")
b = gen._numad_blade
kp = b.keypoints
print([a for a in dir(kp) if not a.startswith("_")])
print([a for a in dir(b) if not a.startswith("_")])
print(np.shape(kp.key_arcs), np.shape(b.ispan))
pickle.dump(dict(key_arcs=np.asarray(kp.key_arcs), key_labels=getattr(kp,"key_labels",None), ispan=np.asarray(b.ispan)
                 ), open(sys.argv[1], "wb"))
