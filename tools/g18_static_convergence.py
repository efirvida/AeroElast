"""G18 statics: AeroElast and CalculiX S8R vs an independent hand-CLT strip solution, by mesh.

For each layup/case it prints the free-edge mean response of both codes at 4x10, 8x20, 16x40 and
32x80, next to the analytical strip result from ``_hand_clt_abd`` (first-principles CLT of
``test_composite_beam_parity.py``, independent of the code's laminate class):

* strip free in width, uniform edge resultant N_y = F/b and zero moments: the 6x6 ABD inverse
  gives eps_y (axial) and kappa_y (b_coupling, w = kappa_y L^2 / 2);
* tip load F in z: m_y = F (L - y) / b, kappa_y = (ABD^-1)[kappa_y, M_y] m_y, w = d F L^3 / (3 b).

Fibre 0 deg lies along the width X (read off the first run, which printed both conventions). Run on a compute node through scripts/aeroenv.sh.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_layup_parity as m  # noqa: E402
from test_composite_beam_parity import _hand_clt_abd  # noqa: E402


def analytic(angles):
    """Hand-CLT strip result. Fibre 0 deg lies along the strip width X (read off the runs: the
    uni_0 axial response along Y is the matrix-dominated one), so the length axis is CLT y."""
    A, Bm, D = _hand_clt_abd(m.E1, m.E2, m.G12, m.NU12, list(angles), m.THICKNESS / len(angles))
    inv = np.linalg.inv(np.block([[A, Bm], [Bm, D]]))
    iy, ik = 1, 4
    n_y = 1000.0 / m.B
    return {
        "axial": inv[iy, iy] * n_y * m.L,
        "b_coupling": -0.5 * inv[ik, iy] * n_y * m.L**2,  # CLT kappa = -w''
        "bending": inv[ik, ik] * 100.0 * m.L**3 / (3.0 * m.B),
    }


ccx = m.ccx_bin_or_skip()
for nx, ny in ((4, 10), (8, 20), (16, 40), (32, 80)):
    mesh = m._build_mesh(nx, ny)
    face = [mesh.node_id_to_index[n.id] for n in mesh.get_node_set("free_face").nodes.values()]
    print(f"== mesh {nx}x{ny}", flush=True)
    for name, angles in m.LAYUPS.items():
        lam = m._laminate(angles)
        an = analytic(angles)
        for case, (dof, load, mdof) in m.CASES.items():
            u = m._aero_solution(mesh, lam, dof, load)
            aero = float(np.mean([u[6 * k + mdof] for k in face]))
            with tempfile.TemporaryDirectory(dir="logs/check/l2-audit") as td:
                disp, ids = m._ccx_solution(mesh, lam, dof, load, ccx, Path(td))
            cc = float(np.mean([disp[i][mdof] for i in ids]))
            print(
                f"{name:10s} {case:10s} aero={aero: .5e} ccx={cc: .5e} "
                f"analytic={an[case]: .5e}",
                flush=True,
            )
