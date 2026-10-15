"""G26 CalculiX deck row: paired mesh convergence of element twist and CCX S8R twist.

Predeclared criterion (written before the first run). For each layup of
test_independent_ccx_deck_confirms_the_coupon ([45,0,0,45]s, [60,0,0,60]s, [45,-45,-45,45]s) the
free-edge least-squares twist is computed by the element (_run_coupon) and by the hand-authored
CCX S8R deck (ccx_tip_twist) on the same ladder 2x12, 4x24, 8x48, 16x96 (nx x ny). The test pairs
the element at 16x96 with CCX at 8x48 and calls that "matched refinement (converged)"; this tool
checks it. A layup is CLOSED when (a) both codes change < 1 % from 8x48 to 16x96 and (b) the gap
element/CCX at 16x96 vs 16x96 is <= 3 % (the test bound). If (a) fails the comparison is
discretization-limited and the row stays code parity. Run through scripts/aeroenv.sh with
PYTHONPATH=$SCRATCH/validenv.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, ".")
import tests.validation.blade.test_laminate_bend_twist as t  # noqa: E402

root = Path("logs/check/l2-audit")
ladder = ((2, 12), (4, 24), (8, 48), (16, 96))
for name in ("[45,0,0,45]s", "[60,0,0,60]s", "[45,-45,-45,45]s"):
    lam, half = t._laminate(t.LAYUPS[name]), t.LAYUPS[name]
    el, cc = {}, {}
    for nx, ny in ladder:
        coords, u, tips = t._run_coupon(lam, nx=nx, ny=ny)
        el[(nx, ny)] = t.twist_lsq(coords, u, tips)
        with tempfile.TemporaryDirectory(dir=root) as td:
            cc[(nx, ny)] = t.ccx_tip_twist(half, nx=nx, ny=ny, workdir=td)
        print(f"{name:>18} {nx:>2}x{ny:<3} element {np.rad2deg(el[(nx, ny)]):+.6f} deg  "
              f"ccx {np.rad2deg(cc[(nx, ny)]):+.6f} deg  el/ccx-1 {100*(el[(nx, ny)]/cc[(nx, ny)]-1):+.3f} %",
              flush=True)
    k = ladder
    de, dc = 100 * (el[k[3]] / el[k[2]] - 1), 100 * (cc[k[3]] / cc[k[2]] - 1)
    gap = 100 * abs(el[k[3]] / cc[k[3]] - 1)
    ok = abs(de) < 1 and abs(dc) < 1 and gap <= 3
    print(f"  => {name}: 8x48->16x96 element {de:+.3f} %, ccx {dc:+.3f} %, gap {gap:.3f} %: "
          f"{'CLOSED' if ok else 'OPEN'} under the predeclared criterion", flush=True)
