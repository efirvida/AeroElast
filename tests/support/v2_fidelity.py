"""Material fidelity of a cached v2 skin mesh against the windIO YAML laminate (odd/tasks/gmsh-store-campaign.md, R2).

The assigned areal mass per material (the layup of each element) is compared with the YAML laminate sampled on a 3x3 grid
inside every skin element at the element's own chart points, so the reference does not depend on how the mesher snapped
or merged the layer boundaries.  Same measurement as ``logs/check/gmsh_merge/matbal.py`` (the numbers of
gmsh-validation-15-22 L23), packaged for the tests.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_GRID = [(a, b) for a in (1 / 6, 0.5, 5 / 6) for b in (1 / 6, 0.5, 5 / 6)]


def skin_material_balance(raw: dict, yaml_file) -> dict:
    """``{"mass_src_kg", "mass_asg_kg", "per_material": {name: (src_kg, asg_kg)}, "frac_gt5": fraction of skin elements whose
    assigned areal mass differs from the sampled source by more than 5 %}``."""
    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import lookup_cells
    from aeroelast.models.blade.numad.objects.blade import Blade

    blade = Blade()
    blade.read_yaml(str(Path(yaml_file)))
    topology = blade.windio2_topology
    tables = json.loads(bytes(raw["tables_json"]).decode())
    rho = {n: float(np.mean(np.asarray(m["rho"], float))) for n, m in tables["materials"].items()}
    mats = sorted(rho)
    skin = np.nonzero(np.asarray(raw["surfaces"]) == "skin")[0]
    pts, own = [], []
    for k in skin:
        e, p = raw["elements"][k], raw["params"][k]
        for a, b in _GRID:
            if e[3] >= 0:
                w = np.array([(1 - a) * (1 - b), a * (1 - b), a * b, (1 - a) * b])
                pts.append(w @ p[:4])
            else:
                w = np.array([1 - a * 0.5 - b * 0.5, a * 0.5, b * 0.5])
                pts.append(w @ p[:3])
            own.append(k)
    pts = np.array(pts)
    cells = lookup_cells(topology, "skin", pts[:, 0], np.clip(pts[:, 1], 0, 1))

    def areal(layers):
        v = dict.fromkeys(mats, 0.0)
        for name, thickness in layers:
            v[name] += rho[name] * thickness
        return v

    src = [areal([(s.material_name, s.thickness_m) for s in c.layers]) for c in cells]
    ne = len(skin)
    S = np.array([[s[m] for m in mats] for s in src]).reshape(ne, len(_GRID), len(mats)).mean(1)
    layups = tables["layups"]
    A = np.array([[areal([(n, th) for n, th, _ in layups[i]])[m] for m in mats] for i in raw["layup_id"][skin]])
    area = np.asarray(raw["areas"])[skin]
    tot_s, tot_a = (S * area[:, None]).sum(0), (A * area[:, None]).sum(0)
    dev = np.abs(A.sum(1) - S.sum(1)) / np.maximum(S.sum(1), 1e-9)
    return {"mass_src_kg": float(tot_s.sum()), "mass_asg_kg": float(tot_a.sum()),
            "per_material": {m: (float(a), float(b)) for m, a, b in zip(mats, tot_s, tot_a, strict=True)},
            "frac_gt5": float((dev > 0.05).mean())}
