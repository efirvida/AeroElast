"""B3-3 M2: mesh-laminate EA / EI integrals with A11 vs the free-hoop compliance modulus (no solve).

The S-1 EA row and the recorded Steiner EI partition (tools/s1_clt_steiner_ei.py) weight each
element by A11 = sum Qbar11 h, i.e. the span stiffness with the hoop and in-plane shear strains
held at zero. A beam section (BeamDyn's SONATA/ANBA4 6x6 included) lets the wall contract freely,
which gives 1/a11 with a = A^-1 (full 3x3 in-plane A). This script prints both integrals per
station against the 6x6 deck over 0.2-0.85 R (odd/tasks/validation-l3-b3-audit.md, "B3-3
rejected-rows analysis", K12.M2). The A11 columns must reproduce the recorded numbers (S-1 EA
+7.57 %, Steiner mid flap +7.74 % / edge -0.40 %) or the comparison is void.

    scripts/aeroenv.sh python -B tools/diagnose_clt_compliance_steiner.py <oml|mid> [element_size]
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def laminate_a(layup, mats) -> np.ndarray:
    """In-plane A [N/m] (3x3, order 1 = span reference direction, 2 = hoop, 6 = shear).

    Same ply data and reciprocal-Poisson Q as run_s1_sectional.laminate_a11, so A[0, 0] equals it.
    """
    a = np.zeros((3, 3))
    for mat_name, h, ang in layup:
        el = mats[mat_name]["elastic"]
        if isinstance(el["E"], (list, tuple)):
            E1, E2, nu12 = float(el["E"][0]), float(el["E"][1]), float(el["nu"][0])
            G12 = float(el["G"][0]) if el.get("G") else E1 / (2 * (1 + nu12))
        else:
            E1 = E2 = float(el["E"])
            nu12 = float(el["nu"])
            G12 = float(el["G"]) if el.get("G") else E1 / (2 * (1 + nu12))
        den = 1.0 - nu12 * nu12 * E2 / E1
        q11, q22, q12, q66 = E1 / den, E2 / den, nu12 * E2 / den, G12
        t = np.deg2rad(ang)
        c, s = np.cos(t), np.sin(t)
        b11 = q11 * c**4 + 2 * (q12 + 2 * q66) * s**2 * c**2 + q22 * s**4
        b22 = q11 * s**4 + 2 * (q12 + 2 * q66) * s**2 * c**2 + q22 * c**4
        b12 = (q11 + q22 - 4 * q66) * s**2 * c**2 + q12 * (s**4 + c**4)
        b66 = (q11 + q22 - 2 * q12 - 2 * q66) * s**2 * c**2 + q66 * (s**4 + c**4)
        b16 = (q11 - q12 - 2 * q66) * c**3 * s - (q22 - q12 - 2 * q66) * c * s**3
        b26 = (q11 - q12 - 2 * q66) * c * s**3 - (q22 - q12 - 2 * q66) * c**3 * s
        a += float(h) * np.array([[b11, b12, b16], [b12, b22, b26], [b16, b26, b66]])
    return a


def section_integrals(gen, mesh, ext):
    """Rows (z, EA, EI_lo, EI_hi) per station for weights A11 and 1/a11 (Steiner, Gauss points)."""
    mem = ext.station_membership()
    sts = sorted(mem)
    coords = np.asarray([n.coords for n in mesh.nodes])
    idx = mesh.node_id_to_index
    mats = {m["name"]: m for m in gen.numad_mesh_data["materials"]}
    sections = {s["elementSet"]: s for s in gen.numad_mesh_data["sections"]}
    elem_set = {}
    for name, es in mesh.element_sets.items():
        for e in es.elements:
            elem_set.setdefault(e.id, name)
    mods = {}
    for eid, n in elem_set.items():
        if n in sections:
            a = laminate_a(sections[n]["layup"], mats)
            mods[eid] = (a[0, 0], 1.0 / np.linalg.inv(a)[0, 0])
    gp = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3.0)
    out = {0: [], 1: []}
    for st in sts:
        if st - 1 not in mem:
            continue
        z = mem[st]["z"]
        dz = max(z - mem[st - 1]["z"], 1e-9)
        ns = set(mem[st]["nodes"].tolist())
        pts, w = [], []
        for e in mesh.elements:
            m = mods.get(e.id)
            ids = [idx[n] for n in e.node_ids]
            if m is None or any(i not in ns for i in ids) or len(ids) != 4:
                continue
            c = coords[ids]
            area = 0.5 * np.linalg.norm(np.cross(c[2] - c[0], c[3] - c[1]))
            for xi, eta in gp:
                N = .25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta),
                                    (1 - xi) * (1 + eta)])
                pts.append(N @ c[:, :2])
                w.append(np.array(m) * area / 4 / dz)
        pts, w = np.array(pts), np.array(w)
        for k in (0, 1):
            wk = w[:, k]
            xc = (wk[:, None] * pts).sum(0) / wk.sum()
            d = pts - xc
            blk = np.array([[(wk * d[:, 1] ** 2).sum(), -(wk * d[:, 0] * d[:, 1]).sum()],
                            [-(wk * d[:, 0] * d[:, 1]).sum(), (wk * d[:, 0] ** 2).sum()]])
            lo, hi = np.linalg.eigvalsh(blk)
            out[k].append((z, wk.sum(), lo, hi))
    return {k: np.array(v) for k, v in out.items()}


def main():
    reference = sys.argv[1]
    size = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=size,
                    skin_reference=reference)
    mesh = gen.generate(renumber="rcm")
    ext = SectionalExtractor(mesh, build_rust_properties(gen.numad_mesh_data))
    tables = section_integrals(gen, mesh, ext)
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    print(f"CLT-STEINER {reference} es={size} (mean signed / mean |.| vs 6x6 over 0.2-0.85 R)")
    for k, label in ((0, "A11 (hoop clamped)"), (1, "1/a11 (hoop free)")):
        t = tables[k]
        zmax = t[:, 0].max()
        m = (t[:, 0] / zmax > .2) & (t[:, 0] / zmax < .85)
        z = t[m, 0]
        e = np.column_stack([t[m, 1] / np.interp(z, bd.z, bd.ea) - 1,
                             t[m, 2] / np.interp(z, bd.z, bd.ei_flap_ec) - 1,
                             t[m, 3] / np.interp(z, bd.z, bd.ei_edge_ec) - 1])
        s, a = 100 * e.mean(axis=0), 100 * np.abs(e).mean(axis=0)
        print(f"  {label:20s} n={m.sum()}  EA {s[0]:+6.2f} / {a[0]:5.2f} %   EI flap {s[1]:+6.2f} / {a[1]:5.2f} %"
              f"   EI edge {s[2]:+6.2f} / {a[2]:5.2f} %")
    r = tables[1][:, 1:] / tables[0][:, 1:]
    print(f"  ratio (1/a11)/A11 over all stations: EA min {r[:, 0].min():.4f} median {np.median(r[:, 0]):.4f}"
          f"  EI flap median {np.median(r[:, 1]):.4f}  EI edge median {np.median(r[:, 2]):.4f}")


if __name__ == "__main__":
    main()
