"""B3-3 K15: split the shell's tip-load deflection into plane-section bending and shear/distortion (no spline).

For each station the mean node displacement ubar_i(z) (i = loaded direction) and the plane-section slope
phi_i(z) = -d(u_z)/d(x_i) (least-squares plane on the station nodes, as K14) are taken from the same solve.
Per chord: total = ubar(z_b) - ubar(z_a); bending = int phi dz (trapezoid on the two stations); shear/distortion = total - bending.
BeamDyn reference from the deck 6x6 at the same stations: bending V... (L-z)^2/EI_eig, shear int dz/GA with GA = 1/C[k,k]
(C = K^-1, BD index 1 = edge = shell +x, index 0 = flap = shell +y), as tip compliance in m/MN.

    scripts/aeroenv.sh python -u -B tools/diagnose_edge_shear_split.py <es>
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT))

from beam_reference import load_beamdyn_blade  # noqa: E402
from diagnose_bending_compliance_2x2 import F, solve  # noqa: E402
from v05_deck_consistency import load_beamdyn  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from aeroelast.postprocess.sectional import SectionalExtractor  # noqa: E402


def main():
    es = float(sys.argv[1])
    gen = BladeMesh(yaml_file=str(ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=es, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    mem = SectionalExtractor(mesh, props).station_membership()
    sts = sorted(mem)
    zs = np.array([mem[st]["z"] for st in sts])
    L = zs.max()
    coords = mesh.coords_array
    deck = load_beamdyn(ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat")
    eta = np.array([d[0] for d in deck]) * 117.0
    ga = {0: np.array([1 / np.linalg.inv(K)[0, 0] for _, K, _ in deck]), 1: np.array([1 / np.linalg.inv(K)[1, 1] for _, K, _ in deck])}
    ei = {1: np.array([np.sort(np.linalg.eigvalsh(K[3:5, 3:5] - np.outer(K[2, 3:5], K[2, 3:5]) / K[2, 2]))[1] for _, K, _ in deck]),
          0: np.array([np.sort(np.linalg.eigvalsh(K[3:5, 3:5] - np.outer(K[2, 3:5], K[2, 3:5]) / K[2, 2]))[0] for _, K, _ in deck])}
    zz = np.linspace(0.0, L, 4001)
    for name, comp, bdk in (("EDGE (+x)", 0, 1), ("FLAP (+y)", 1, 0)):
        u = solve(mesh, props, mem[sts[-1]]["nodes"], comp)
        ubar, phi = np.zeros(len(sts)), np.zeros(len(sts))
        for k, st in enumerate(sts):
            n = mem[st]["nodes"]
            xy = coords[n, :2] - coords[n, :2].mean(axis=0)
            coef, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(n)), xy]), u[n, 2], rcond=None)
            ubar[k], phi[k] = u[n, comp].mean(), -coef[1 + comp]
        # root contribution: stations start above z = 0
        bend = np.concatenate([[0.0], np.cumsum(0.5 * (phi[1:] + phi[:-1]) * np.diff(zs))]) + 0.5 * phi[0] * zs[0]
        tot = ubar
        shear = tot - bend
        cb = np.trapezoid((L - zz) ** 2 / np.interp(zz, eta, ei[bdk]), zz) / 1e0
        cs = np.trapezoid(1.0 / np.interp(zz, eta, ga[bdk]), zz)
        print(f"== {name}: tip station mean ubar {1e6 * ubar[-1] / F:.3f} m/MN  plane-section bending {1e6 * bend[-1] / F:.3f}  shear+distortion {1e6 * shear[-1] / F:.3f}"
              f"   (shear share {100 * shear[-1] / ubar[-1]:.1f} %)")
        print(f"   BeamDyn deck (EB bending {1e6 * cb:.3f} + shear {1e6 * cs:.3f} = {1e6 * (cb + cs):.3f} m/MN; shear share {100 * cs / (cb + cs):.1f} %)")
        print("   z_mid   dBend/dTot-chord  shell GA_eff/deck_GA-1 [%]   shell EI_plane/deck-1 [%] (3-station chord)")
        for a in range(0, len(sts) - 3):
            b = a + 3
            zm = 0.5 * (zs[a] + zs[b])
            if not (0.2 * L < zm < 0.85 * L):
                continue
            ch = np.linspace(zs[a], zs[b], 201)
            dsh = (ubar[b] - ubar[a]) - (bend[b] - bend[a])
            gaeff = F * (zs[b] - zs[a]) / dsh if abs(dsh) > 1e-30 else np.nan
            m_int = F * np.trapezoid(L - ch, ch)
            eip = m_int / (phi[b] - phi[a])
            ref_ga = (zs[b] - zs[a]) / np.trapezoid(1 / np.interp(ch, eta, ga[bdk]), ch)
            ref_ei = (zs[b] - zs[a]) / np.trapezoid(1 / np.interp(ch, eta, ei[bdk]), ch)
            print(f"  {zm:7.2f} {(bend[b] - bend[a]) / (ubar[b] - ubar[a]):12.3f} {100 * (gaeff / ref_ga - 1):+22.1f} {100 * (eip / ref_ei - 1):+24.1f}")


if __name__ == "__main__":
    main()
