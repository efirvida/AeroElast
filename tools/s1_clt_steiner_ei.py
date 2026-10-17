"""Principal EI of the MESH laminate by a CLT Steiner integral (no solve), vs the yaml integral.

Partitions the S-1 EI gap: (a) yaml integral (tools/yaml_ei_independent.py, contour + yaml
angles) vs this number = mapping/geometry/angle of the mesh; (b) this number vs the static-solve
EI of tools/s1_reference_sweep.py = shell response. Weight per element Gauss point:
A11 * area / dz / 4 (A11 from the mesh section layup, same as the S-1 EA(CLT) path).

    scripts/aeroenv.sh python -B tools/s1_clt_steiner_ei.py <oml|mid> [element_size]
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from run_s1_sectional import laminate_a11  # noqa: E402
from yaml_ei_independent import principal_ei  # noqa: E402
from yaml_layup_independent import load  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402



def clt_ei_table(gen, mesh, ext):
    """Return rows (z, EI_flap_CLT, EI_edge_CLT) per station (min/max eigenvalue of the A11 block)."""
    membership = ext.station_membership()
    sts = sorted(membership)
    coords = np.asarray([n.coords for n in mesh.nodes])
    idx = mesh.node_id_to_index
    mats = {m["name"]: m for m in gen.numad_mesh_data["materials"]}
    sections = {s["elementSet"]: s for s in gen.numad_mesh_data["sections"]}
    elem_set = {}
    for name, es in mesh.element_sets.items():
        for e in es.elements:
            elem_set.setdefault(e.id, name)
    a11 = {eid: laminate_a11(sections[n]["layup"], mats) for eid, n in elem_set.items() if n in sections}
    gp = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3.0)
    out = []
    for st in sts:
        if st - 1 not in membership:
            continue
        z = membership[st]["z"]
        dz = max(z - membership[st - 1]["z"], 1e-9)
        ns = set(membership[st]["nodes"].tolist())
        pts, w = [], []
        for e in mesh.elements:
            a = a11.get(e.id)
            ids = [idx[n] for n in e.node_ids]
            if a is None or any(i not in ns for i in ids) or len(ids) != 4:
                continue
            c = coords[ids]
            area = 0.5 * np.linalg.norm(np.cross(c[2] - c[0], c[3] - c[1]))
            for xi, eta in gp:
                N = .25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
                pts.append(N @ c[:, :2]); w.append(a * area / 4 / dz)
        pts, w = np.array(pts), np.array(w)
        xc = (w[:, None] * pts).sum(0) / w.sum()
        d = pts - xc
        blk = np.array([[(w * d[:, 1] ** 2).sum(), -(w * d[:, 0] * d[:, 1]).sum()],
                        [-(w * d[:, 0] * d[:, 1]).sum(), (w * d[:, 0] ** 2).sum()]])
        ev = np.linalg.eigvalsh(blk)
        out.append((z, ev[0], ev[1]))
    return np.array(out)


def main():
    reference = sys.argv[1]
    size = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size,
                    skin_reference=reference)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ext = SectionalExtractor(mesh, props)
    table = clt_ei_table(gen, mesh, ext)
    data = load()
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    zmax = table[:, 0].max()
    rows = []
    for z, f, e in table:
        eta_s = z / zmax
        if not .2 < eta_s < .85:
            continue
        yf, ye, _ = principal_ei(data, eta_s)
        bf, be = np.interp(z, bd.z, bd.ei_flap_ec), np.interp(z, bd.z, bd.ei_edge_ec)
        rows.append((f / bf - 1, e / be - 1, f / yf - 1, e / ye - 1))
    r = 100 * np.mean(rows, axis=0)
    print(f"STEINER {reference} es={size} n={len(rows)}  mesh-CLT vs 6x6: flap {r[0]:+.2f} % edge {r[1]:+.2f} %"
          f" | mesh-CLT vs yaml-integral: flap {r[2]:+.2f} % edge {r[3]:+.2f} %")


if __name__ == "__main__":
    main()
