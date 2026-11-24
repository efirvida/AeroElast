"""L3 T10 (W7): end-region shell GJ stiffening variant (D3 sensitivity) as a ``mesh_transform``.

G12 of every orthotropic ply (isotropic plies untouched, as W6 V1) times ``factor`` for elements whose centroid z >= Z1, linearly
ramped from 1 at Z0 to ``factor`` at Z1 (stepwise, 0.5 m z bins of the element centroid). Same mesh/property machinery as
``tools/l3_w6_scale.make_transform`` (one element set per layup set and factor bin).

    variant spec string: ``tipstiff:<factor>:<z0>:<z1>``   e.g. tipstiff:1.2821:96:98
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import l3_w6_scale as sc  # noqa: E402


def make_transform(factor, z0, z1):
    def fz(z):
        return 1.0 + (factor - 1.0) * float(np.clip((z - z0) / (z1 - z0), 0.0, 1.0))

    def transform(mesh, props):
        from aeroelast.core.mesh.entities import ElementSet
        from aeroelast.models.blade.model import build_rust_properties

        numad = mesh.numad_mesh_data
        mats = {m["name"]: m for m in numad["materials"]}
        layups = {s["elementSet"]: s["layup"] for s in numad["sections"]}
        pos = {n.id: n.coordinates if hasattr(n, "coordinates") else n.coords for n in mesh.nodes}
        new_props, new_mats, secs, nmoved = dict(props), {}, [], 0
        for name in [n for n in mesh.element_sets if n in layups]:
            groups: dict[float, list] = {}
            for e in mesh.element_sets[name].elements:
                zc = float(np.mean([pos[i][2] for i in e.node_ids]))
                zb = (np.floor(zc / 0.5) + 0.5) * 0.5
                groups.setdefault(round(fz(zb), 4), []).append(e)
            for f, g in groups.items():
                if f == 1.0:
                    continue
                sname = f"{name}__T{f:.4f}"
                mesh.element_sets[name].elements.difference_update(g)
                mesh.add_element_set(ElementSet(name=sname, elements=set(g)))
                lay = []
                for n, th, a in layups[name]:
                    nn = sc.key(n, f, 1.0)
                    new_mats.setdefault(
                        nn,
                        {
                            "name": nn,
                            "density": mats[n]["density"],
                            "elastic": sc.scale_material(mats[n]["elastic"], f, 1.0),
                        },
                    )
                    lay.append([nn, th, a])
                secs.append({"elementSet": sname, "layup": lay})
                nmoved += len(g)
        new_props.update(
            build_rust_properties({"materials": list(new_mats.values()), "sections": secs})
        )
        for name in [n for n in layups if n in mesh.element_sets and not mesh.element_sets[n].elements]:
            new_props.pop(name, None)
        print(f"w7 transform: {len(secs)} property sets for {nmoved} elements", flush=True)
        return mesh, new_props

    return transform


def from_spec(spec):
    _, f, z0, z1 = spec.split(":")
    return make_transform(float(f), float(z0), float(z1))
