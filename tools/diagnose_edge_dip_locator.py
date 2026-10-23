"""B3-3 K13.A1/A2: where is the 63-100 m edge dip, and do MITC3 triangles carry it? (no solve)

Per station: hoop-free (1/a11) mesh-laminate edge/flap EI against the 6x6 deck, the FE 2x2
eigenvalue edge/flap EI read from a tools/diagnose_bending_compliance_2x2.py log, their gap, the
MITC3 triangle count and the share of the laminate edge integral they carry. Then the composite
triangle census by component and 10 m bin (odd/tasks/validation-l3-b3-audit.md, K13).

    scripts/aeroenv.sh python -u -B tools/diagnose_edge_dip_locator.py <es> logs/check/b3-3/m1-<es>-<job>.log
"""

import collections
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from diagnose_clt_compliance_steiner import section_integrals  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402

DIP_FROM = 62.0  # m, K13 window W = [62 m, 0.85 R]; inboard window I = [0.2 R, 62 m)


def component(name: str) -> str:
    tail = "_".join(name.split("_")[2:])
    for key, label in (("SW", "web"), ("TE_FLAT", "TE_flat"), ("TE_PANEL", "TE_panel"), ("LE_PANEL", "LE_panel"),
                       ("SPAR", "spar_cap"), ("CAP", "spar_cap")):
        if key in tail:
            return label
    return tail or name


def read_fe(path):
    rows = [tuple(map(float, l.split())) for l in Path(path).read_text().splitlines()
            if re.match(r"^\s+\d+\.\d+\s+[+-]", l)]
    return np.array(rows)  # z, single f, single e, 2x2 f, 2x2 e, angle (percent)


def main():
    size = float(sys.argv[1])
    fe = read_fe(sys.argv[2])
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    ext = SectionalExtractor(mesh, build_rust_properties(gen.numad_mesh_data))
    quad = section_integrals(gen, mesh, ext)[1]
    full = section_integrals(gen, mesh, ext, include_triangles=True)[1]
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    zmax = quad[:, 0].max()
    mem = ext.station_membership()
    coords = np.asarray([n.coords for n in mesh.nodes])
    idx = mesh.node_id_to_index
    es = {}
    for name, s in mesh.element_sets.items():
        for e in s.elements:
            es.setdefault(e.id, name)
    tris = [(coords[[idx[n] for n in e.node_ids]].mean(axis=0)[2], component(es.get(e.id, "?")))
            for e in mesh.elements if len(e.node_ids) == 3]
    stz = np.array([mem[st]["z"] for st in sorted(mem)])
    print(f"EDGE-DIP-LOCATOR mid es={size}  stations {len(quad)}  MITC3 triangles {len(tris)}")
    print("  z[m]   lam-edge  lam-flap  FE-edge  FE-flap  gap-edge[pp]  tri/station  tri-share-edge[%]")
    out = []
    for q, f in zip(quad, full):
        z = q[0]
        if not (0.2 < z / zmax < 0.85):
            continue
        j = int(np.argmin(np.abs(fe[:, 0] - z)))
        if abs(fe[j, 0] - z) > 0.02:
            continue
        le = q[3] / np.interp(z, bd.z, bd.ei_edge_ec) - 1
        lf = q[2] / np.interp(z, bd.z, bd.ei_flap_ec) - 1
        fee, fef = fe[j, 4] / 100, fe[j, 3] / 100
        dz = np.diff(stz)[np.argmin(np.abs(stz[1:] - z))] if len(stz) > 1 else 1.0
        nt = sum(1 for tz, _ in tris if abs(tz - z) <= dz / 2)
        share = 100 * (1 - q[3] / f[3])
        out.append((z, le, lf, fee, fef, 100 * (fee - le), nt, share))
        print(f"  {z:6.2f} {100 * le:+8.2f} {100 * lf:+9.2f} {100 * fee:+8.2f} {100 * fef:+8.2f} {100 * (fee - le):+11.2f}"
              f" {nt:9d} {share:12.2f}")
    r = np.array(out)
    for label, m in (("I inboard", r[:, 0] < DIP_FROM), ("W dip", r[:, 0] >= DIP_FROM), ("all", r[:, 0] > 0)):
        x = r[m]
        print(f"  {label:10s} n={m.sum():2d}  lam-edge {100 * x[:, 1].mean():+6.2f}  FE-edge {100 * x[:, 3].mean():+6.2f}"
              f"  gap {x[:, 5].mean():+6.2f} pp  lam-flap {100 * x[:, 2].mean():+6.2f}  FE-flap {100 * x[:, 4].mean():+6.2f}"
              f"  tri-share mean {x[:, 7].mean():.2f} max {x[:, 7].max():.2f} %")
    worst = r[r[:, 3] <= -0.10]
    print(f"  stations with FE edge <= -10 %: {len(worst)}; tri-share mean {worst[:, 7].mean() if len(worst) else float('nan'):.2f} %,"
          f" gap mean {worst[:, 5].mean() if len(worst) else float('nan'):+.2f} pp")
    print("  triangle census (10 m bins x component):")
    cnt = collections.Counter((int(tz // 10) * 10, c) for tz, c in tris)
    for (b, c), n in sorted(cnt.items()):
        print(f"    {b:3d}-{b + 10:3d} m  {c:10s} {n}")
    print("  triangle census by component:", dict(collections.Counter(c for _, c in tris)))


if __name__ == "__main__":
    main()
