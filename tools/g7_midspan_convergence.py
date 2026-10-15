"""G7 paired mesh convergence, second metric: outer-fibre stress at mid-span, away from the clamp.

Why a second probe: tools/g7_stress_convergence.py showed the first metric (max von Mises over the
whole plate) does NOT converge in either code -- both rise past M c / I at the root as the mesh
refines (128x32: aero 61.8, ccx 65.9 MPa against 60.0) -- the clamped edge (six DOF held)
is a stress singularity that grows with refinement.

Predeclared BEFORE running this probe (2026-10-15):
  metric   outer-fibre sigma_xx at the plate centre line y = B/2, at the first station at or past
           x = L/2, normalised by the closed form M(x) c / I at that same x
           (AeroElast element value sits at its centroid x = (NX/2 - 1/2) h, CalculiX nodal value at x = L/2).
  meshes   NX x NY = 8x2, 16x4, 32x8, 64x16, 128x32 (same as the first probe; same load/material/BC).
  converged  a code is converged when the ratio changes < 1 % between the last two meshes.
  admissible pair  both converged AND the two ratios agree within 5 %; the analytical comparison
           then uses the converged ratio against 1.0 with the 5 % bound.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import tests.validation.parity.test_shell_stress_ccx_parity as m  # noqa: E402
from aeroelast.postprocess.stress_recovery import StressLocation, StressRecovery, StressType  # noqa: E402
from tests.support.ccx_io import parse_frd_stress  # noqa: E402


def frd_coords(path: Path) -> dict[int, np.ndarray]:
    out, on = {}, False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("    2C"):
            on = True
            continue
        if on and line.startswith(" -3"):
            break
        if on and line.startswith(" -1"):
            out[int(line[3:13])] = np.array([float(line[13 + 12 * j:25 + 12 * j]) for j in range(3)])
    return out


def m_over_i(x):
    return m.FORCE * (m.L - x) * (m.H / 2) / (m.B * m.H**3 / 12)


res = []
for nx, ny in [(8, 2), (16, 4), (32, 8), (64, 16), (128, 32)]:
    m.NX, m.NY = nx, ny
    mesh = m._build_plate()
    domain, u = m._solve_aeroelast(mesh)
    r = StressRecovery(domain, u).compute_element_stresses(
        location=StressLocation.TOP, stress_type=StressType.TOTAL)
    vm = np.asarray(r.von_mises)
    h = m.L / nx
    idx = (ny // 2) * nx + (nx // 2 - 1)
    xc = (nx // 2 - 0.5) * h
    aero = vm[idx] / m_over_i(xc)
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        m._run_ccx_outer_fibre(d, mesh)
        coords = frd_coords(d / "plate.frd")
        st = parse_frd_stress(d / "plate.frd")
    sel = [n for n, c in coords.items()
           if abs(c[0] - m.L / 2) < 1e-9 and abs(c[1] - m.B / 2) < 1e-9 and c[2] > 0 and n in st]
    if not sel:
        print(nx, ny, "no ccx node found at mid-span; coords sample", list(coords.items())[:3])
        continue
    sxx = abs(float(np.mean([st[n][0] for n in sel])))
    ccx = sxx / m_over_i(m.L / 2)
    res.append((aero, ccx))
    print(f"{nx}x{ny}: aero/MoverI {aero:.4f}  ccx/MoverI {ccx:.4f}  aero/ccx-1 {aero/ccx-1:+.2%}"
          f"  (ccx nodes {len(sel)})", flush=True)
for name, i in (("aero", 0), ("ccx", 1)):
    a, b = res[-2][i], res[-1][i]
    print(f"{name}: last-two-mesh change {abs(b-a)/b:.3%}")
