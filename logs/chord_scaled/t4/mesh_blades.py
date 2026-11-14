"""Mesh one blade/case with the gmsh mesher and write a VTU (aspect, min_angle, jacobian, is_web) plus a counts JSON.

Usage: mesh_blades.py <blade> <case>   blade in iea34, iea15, iea22; case in uniform, n10, n20.
The uniform size is each blade's converged size (IEA-3.4 0.1, IEA-15 0.25, IEA-22 0.177 m).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "tests/support"))
BLADES = {
    "iea34": (REPO / "tests/IEA34/reference/upstream/IEA-3.4-130-RWT.yaml", 0.1),
    "iea15": (REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml", 0.25),
    "iea22": (REPO / "tests/IEA22MW/reference/upstream/IEA-22-280-RWT.yaml", 0.177),
}
CASES = {"uniform": {}, "n10": dict(chord_divisions=10, span_aspect=1.0, min_band_elements=2),
         "n20": dict(chord_divisions=20, span_aspect=1.0, min_band_elements=2)}


def main(blade, case):
    import meshio
    import v2_blade
    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import element_metrics, quality_gate

    yaml_file, h = BLADES[blade]
    kw = dict(element_size=h, **CASES[case])
    if blade == "iea34":
        kw["es_range"] = (0.05, 0.5)
    t0 = time.time()
    try:
        raw = v2_blade.load_raw(mesher="gmsh", yaml_file=yaml_file, **kw)
    except TypeError:
        kw.pop("es_range", None)
        raw = v2_blade.load_raw(mesher="gmsh", yaml_file=yaml_file, **kw)
    x = np.asarray(raw["nodes"], float)
    el = np.asarray(raw["elements"])
    skin = np.asarray(raw["surfaces"]) == "skin"
    m = element_metrics(x, el)
    rep = quality_gate(m, raise_on_fail=False)
    cells, cd = [], {}
    for t, sel in (("quad", el[:, 3] >= 0), ("triangle", el[:, 3] < 0)):
        if sel.any():
            cells.append((t, el[sel][:, :4] if t == "quad" else el[sel][:, :3]))
            for name, arr in (("aspect", m["aspect"]), ("min_angle", m["min_angle"]), ("jacobian", m["jacobian"]),
                              ("is_web", (~skin).astype(float))):
                cd.setdefault(name, []).append(arr[sel])
    out = HERE / f"{blade}_{case}.vtu"
    meshio.Mesh(x, cells, cell_data=cd).write(out)
    res = {"blade": blade, "case": case, "args": {k: v for k, v in kw.items()}, "nodes": len(x), "elements": len(el),
           "triangles": int((el[:, 3] < 0).sum()), "gate_passed": rep["passed"], "gate_violations": rep["violations"],
           "min_angle_p5": float(np.percentile(m["min_angle"], 5)), "n_min_angle_lt30": int((m["min_angle"] < 30).sum()),
           "n_min_angle_lt45": int((m["min_angle"] < 45).sum()), "aspect_p95": float(np.percentile(m["aspect"], 95)),
           "wall_s": time.time() - t0, "vtu": str(out)}
    (HERE / f"{blade}_{case}.json").write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({k: res[k] for k in ("blade", "case", "nodes", "elements", "gate_passed", "n_min_angle_lt45", "wall_s")}))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
