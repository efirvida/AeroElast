"""G7 paired mesh convergence: outer-fibre stress of the clamped plate, AeroElast vs CalculiX.

Predeclared BEFORE any run (2026-10-15):
  metric   max element von Mises at the outer fibre (the metric of
           tests/validation/parity/test_shell_stress_ccx_parity.py), same load, material, BC.
  meshes   NX x NY = 8x2, 16x4, 32x8, 64x16 (aspect ratio held constant); 128x32 added after the first read showed CalculiX still moving 2.56 %.
  converged  a code is converged when the change between the last two meshes is < 1 %.
  verdict  - both converge AND agree with each other within 5 %: the pair is admissible;
             the analytical M c / I is then compared at the converged value.
           - the analytical M c / I at x = 0 is not the quantity either code reports (an element
             value sits at x = h_x/2, where M c / I is (1 - h_x/(2L)) of the root value); the
             probe also prints that centroid-corrected reference, as a diagnosis, not a new bound.
No tolerance is changed by this script; it writes nothing to the store.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import tests.validation.parity.test_shell_stress_ccx_parity as m  # noqa: E402
from aeroelast.postprocess.stress_recovery import StressLocation, StressRecovery, StressType  # noqa: E402
import numpy as np  # noqa: E402

rows = []
for nx, ny in [(8, 2), (16, 4), (32, 8), (64, 16), (128, 32)]:
    m.NX, m.NY = nx, ny
    mesh = m._build_plate()
    domain, u = m._solve_aeroelast(mesh)
    rec = StressRecovery(domain, u)
    r = rec.compute_element_stresses(location=StressLocation.TOP, stress_type=StressType.TOTAL)
    aero = float(np.asarray(r.von_mises).max())
    with tempfile.TemporaryDirectory() as d:
        ccx = m._run_ccx_outer_fibre(Path(d), mesh)
    hx = m.L / nx
    ref_c = m.ANALYTICAL_STRESS * (1 - hx / (2 * m.L))
    rows.append((nx, ny, aero / 1e6, ccx / 1e6, ref_c / 1e6))
    print(f"{nx}x{ny}: aero {aero/1e6:.3f}  ccx {ccx/1e6:.3f}  aero/ccx-1 {aero/ccx-1:+.2%}  "
          f"analytic root {m.ANALYTICAL_STRESS/1e6:.2f}  centroid-corrected {ref_c/1e6:.3f}  "
          f"aero/analytic-root-1 {aero/m.ANALYTICAL_STRESS-1:+.2%}", flush=True)
for name, i in (("aero", 2), ("ccx", 3)):
    a, b = rows[-2][i], rows[-1][i]
    print(f"{name}: last-two-mesh change {abs(b-a)/b:.3%}")
