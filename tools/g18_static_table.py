"""Re-tabulate a g18_static_convergence log: error of each code against the hand CLT, per mesh."""
import re
import sys

import numpy as np

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_layup_parity as m  # noqa: E402
from test_composite_beam_parity import _hand_clt_abd  # noqa: E402


def analytic(angles):
    A, Bm, D = _hand_clt_abd(m.E1, m.E2, m.G12, m.NU12, list(angles), m.THICKNESS / len(angles))
    inv = np.linalg.inv(np.block([[A, Bm], [Bm, D]]))
    n_y = 1000.0 / m.B
    return {
        "axial": inv[1, 1] * n_y * m.L,
        "b_coupling": -0.5 * inv[4, 1] * n_y * m.L**2,  # CLT kappa = -w''
        "bending": inv[4, 4] * 100.0 * m.L**3 / (3.0 * m.B),
    }


mesh = None
for line in open(sys.argv[1]):
    if line.startswith("== mesh"):
        mesh = line.split()[2]
        print(f"== {mesh}   (error vs hand CLT, %: aero / ccx)")
        continue
    t = re.match(r"(\S+)\s+(axial|bending|b_coupling)\s+aero=\s*(\S+) ccx=\s*(\S+)", line)
    if not t:
        continue
    name, case, a, c = t.group(1), t.group(2), float(t.group(3)), float(t.group(4))
    ref = analytic(m.LAYUPS[name])[case]
    if abs(ref) < 1e-9:
        print(f"{name:10s} {case:10s} ref~0  aero={a: .2e} ccx={c: .2e}")
    else:
        print(f"{name:10s} {case:10s} {100 * (a / ref - 1):+7.3f} / {100 * (c / ref - 1):+7.3f}   ref={ref: .5e}")
