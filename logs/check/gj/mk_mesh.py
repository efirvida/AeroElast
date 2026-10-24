"""Mesh at ad5f3ac (overlay on PYTHONPATH), es 0.25 mid, dump the Bredt inputs (no solve). argv: out.pkl"""
import sys, pickle, numpy as np
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.model import build_rust_properties
from aeroelast.postprocess.sectional import SectionalExtractor
import aeroelast
print(aeroelast.__file__)
W = "/scratch/leahk/eduardo.donestevez/fem-shell/"
gen = BladeMesh(yaml_file=W + "tests/IEA-15-240-RWT.yaml", element_size=0.25, skin_reference="mid")
mesh = gen.generate(renumber="rcm")
props = build_rust_properties(gen.numad_mesh_data)
ext = SectionalExtractor(mesh, props)
mem = ext.station_membership()
sts = sorted(mem)
zs = np.array([mem[s]["z"] for s in sts])
coords = np.asarray([n.coords for n in mesh.nodes])
cents = np.array([coords[mem[s]["nodes"]].mean(axis=0) for s in sts])
tangs = np.gradient(cents, zs, axis=0)
tangs /= np.linalg.norm(tangs, axis=1)[:, None]
idx = mesh.node_id_to_index
elems = [(e.id, [idx[n] for n in e.node_ids]) for e in mesh.elements]
esets = {n: [e.id for e in es.elements] for n, es in mesh.element_sets.items()}
pickle.dump(dict(coords=coords, mem=dict(mem), sts=sts, zs=zs, tangs=tangs, elems=elems, esets=esets,
                 numad_sections=gen.numad_mesh_data["sections"], materials=gen.numad_mesh_data["materials"]),
            open(sys.argv[1], "wb"))
print("done", len(coords), len(elems), len(sts))
