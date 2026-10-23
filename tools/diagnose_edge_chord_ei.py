"""B3-3 K14: spline-free sectional EI (plane-section slope chords) vs the M1 spline 2x2 and the laminate.

M1 (diagnose_bending_compliance_2x2.py) reads curvature as the second derivative of a smoothing spline
through the section-mean displacement. Here each station's plane-section slope u_i' = -d(u_z)/d(x_i) is
fitted by least squares on the station nodes (u_z = a + b x + c y), and the compliance of a chord
[z_a, z_b] is S_ij = (u_i'^j(z_b) - u_i'^j(z_a)) / (F int (L - z) dz): no spline, no second derivative, exact moment.
Eigenvalue EIs of the symmetrised S are compared with the deck (harmonic mean over the chord).

    scripts/aeroenv.sh python -u -B tools/diagnose_edge_chord_ei.py <es> <chord stations: 1|2|3>
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from diagnose_bending_compliance_2x2 import F, solve  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def main():
    es = float(sys.argv[1])
    chords = [int(c) for c in sys.argv[2:]] or [1]
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=es, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    ext = SectionalExtractor(mesh, props)
    mem = ext.station_membership()
    sts = sorted(mem)
    zs = np.array([mem[st]["z"] for st in sts])
    L = zs.max()
    u = {j: solve(mesh, props, mem[sts[-1]]["nodes"], j) for j in (0, 1)}
    coords = mesh.coords_array
    slope = np.zeros((len(sts), 2, 2))  # [station, i (x,y), load j]
    resid = np.zeros(len(sts))
    for k, st in enumerate(sts):
        n = mem[st]["nodes"]
        xy = coords[n, :2] - coords[n, :2].mean(axis=0)
        A = np.column_stack([np.ones(len(n)), xy])
        for j in (0, 1):
            coef, *_ = np.linalg.lstsq(A, u[j][n, 2], rcond=None)
            slope[k, :, j] = -coef[1:3]
            resid[k] = max(resid[k], np.std(u[j][n, 2] - A @ coef) / (np.abs(u[j][n, 2]).max() + 1e-30))
    bd = load_beamdyn_blade(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    print(f"EDGE-CHORD-EI mid es={es} stations {len(sts)} L={L:.2f}; plane-fit residual median {np.median(resid):.3f} max {resid.max():.3f}")
    for c in chords:
        print(f"-- chord {c} station(s)   z_mid  edge-EI/deck-1 [%]  flap-EI/deck-1 [%]  coupling asym [%]")
        for a in range(0, len(sts) - c):
            b = a + c
            if not (0.2 * L < 0.5 * (zs[a] + zs[b]) < 0.85 * L):
                continue
            zz = np.linspace(zs[a], zs[b], 401)
            m_int = F * np.trapezoid(L - zz, zz)
            S = (slope[b] - slope[a]) / m_int
            ev = np.linalg.eigvalsh(0.5 * (S + S.T))
            ref = lambda arr: (zs[b] - zs[a]) / np.trapezoid(1 / np.interp(zz, bd.z, arr), zz)
            asym = abs(S[0, 1] - S[1, 0]) / np.abs(S).max()
            print(f"  {0.5 * (zs[a] + zs[b]):7.2f} {100 * (1 / ev[0] / ref(bd.ei_edge_ec) - 1):+10.2f} "
                  f"{100 * (1 / ev[1] / ref(bd.ei_flap_ec) - 1):+17.2f} {100 * asym:12.1f}")


if __name__ == "__main__":
    main()
