"""G36: paired convergence of the outer-fibre ply stress, AeroElast production path vs CalculiX S8R.

Predeclared criterion (written before the first run): meshes 8x2 (test), 16x4 (the row),
32x8, 64x16, both layups, free-field centre TOP sigma_xx. A code is converged when its value
changes < 1 % from 32x8 to 64x16. The CLT closed form (the test's own, 5 % bound) is the arbiter:
the CCX comparison is reported against it and a CalculiX value that converges away from the
closed form is a finding about the CalculiX deck, not a failure of the closed-form row.
Run through scripts/aeroenv.sh with PYTHONPATH=$SCRATCH/validenv.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_stress_ccx_parity as m  # noqa: E402
from aeroelast.core.assembler import MeshAssembler  # noqa: E402

MESHES = ((8, 2), (16, 4), (32, 8), (64, 16))
root = Path("logs/check/l2-audit")
StressLocation, StressType = m.StressLocation, m.StressType
for name, angles in m.LAYUPS.items():
    lam = m._laminate(angles)
    sigma_ref, _ = m._hand_clt_ply_stresses(angles, m.N_X)
    ref = float(sigma_ref[-1])
    res = {}
    for nx, ny in MESHES:
        mesh = m._build_mesh(nx, ny)
        cfg = m._model_cfg(lam)
        u = m._solve_production(mesh, cfg)
        dom = MeshAssembler(mesh=mesh, model=cfg)
        r = m.StressRecovery(dom, u).compute_element_stresses(
            location=StressLocation.TOP, stress_type=StressType.TOTAL
        )
        ours = float(r.sigma_xx[m._centre_element_index(dom)])
        with tempfile.TemporaryDirectory(dir=root) as td:
            ccx, _ = m._run_ccx(Path(td), mesh, lam)
        res[(nx, ny)] = (ours, ccx)
        print(f"{name} {nx:>2}x{ny:<2} ours {ours/1e6:.5f} ccx {ccx/1e6:.5f} CLT {ref/1e6:.5f} MPa | "
              f"ours-vs-CLT {100*(ours/ref-1):+.3f}% ccx-vs-CLT {100*(ccx/ref-1):+.3f}% gap {100*abs(ours-ccx)/ccx:.3f}%",
              flush=True)
    (o1, c1), (o2, c2) = res[MESHES[-2]], res[MESHES[-1]]
    do, dc = 100 * (o2 / o1 - 1), 100 * (c2 / c1 - 1)
    print(f"  => {name}: 32x8->64x16 ours {do:+.3f}% ccx {dc:+.3f}%  "
          f"{'both converged' if abs(do) < 1 and abs(dc) < 1 else 'NOT converged'}", flush=True)
