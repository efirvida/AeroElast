"""G9 mesh study: CalculiX S4/S8/S8R tip deflection vs Euler-Bernoulli and Timoshenko closed forms.

Re-measures the table in the docstring of test_ccx_shell_element_types_parity.py and adds
the Timoshenko shear-corrected closed form (kappa = 5/6), the analytical reason S8/S8R
settle above the Euler-Bernoulli value. Run through scripts/aeroenv.sh on a compute node.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")  # repo root first: tools/tests would otherwise shadow `tests`
import test_ccx_shell_element_types_parity as m  # noqa: E402

KAPPA = 5.0 / 6.0
G = m.E / (2.0 * (1.0 + m.NU))
shear = m.FORCE * m.L / (KAPPA * G * m.B * m.H)
timo = m.ANALYTICAL_UY + shear
print(f"EB {m.ANALYTICAL_UY:.7e}  Timoshenko {timo:.7e}  shear share {100 * shear / m.ANALYTICAL_UY:.3f} %")

ccx = m.ccx_bin_or_skip()
print("mesh   " + "  ".join(f"{t:>24}" for t in m.ELEMENT_TYPES))
for nx, ny in ((4, 2), (8, 4), (16, 8), (32, 16), (64, 32)):
    m.NX, m.NY = nx, ny
    cells = []
    for t in m.ELEMENT_TYPES:
        with tempfile.TemporaryDirectory(dir="logs/check/l2-audit") as d:
            uy = m._run_ccx_case(ccx, Path(d), t)
        cells.append(f"{uy:.6e} EB{100 * (uy / m.ANALYTICAL_UY - 1):+.2f}% T{100 * (uy / timo - 1):+.2f}%")
    print(f"{nx}x{ny:<3} " + "  ".join(cells), flush=True)
