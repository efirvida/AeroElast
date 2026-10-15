"""G18 modal: paired convergence of AeroElast and CalculiX S8R, plus beam closed-form anchors.

Predeclared criterion (odd/tasks/validation-l3-b3-audit.md, "Predeclared criterion - G18 modal
paired convergence"): the five modes the test matches are converged in a code when their
frequency changes < 1 % from 16x40 to 32x80; the comparison closes only if all five are converged
in both codes and the matched gap at 32x80 is <= 3 %. Mode identity is tracked by sorted index
within each code. Run on a compute node through scripts/aeroenv.sh.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "tests/validation/parity")
sys.path.insert(0, ".")
import test_composite_layup_parity as m  # noqa: E402
from test_composite_beam_parity import _hand_clt_abd  # noqa: E402

MESHES = ((8, 20), (16, 40), (32, 80))
BETA_L = (1.8751040687, 4.6940911330)


def beam_anchors(angles):
    A, Bm, D = _hand_clt_abd(m.E1, m.E2, m.G12, m.NU12, list(angles), m.THICKNESS / len(angles))
    inv = np.linalg.inv(np.block([[A, Bm], [Bm, D]]))
    mass_per_len = m.RHO * m.THICKNESS * m.B
    ei_out = m.B / inv[4, 4]
    ei_in = m.B**3 / (12.0 * inv[1, 1])
    f = lambda ei, bl: bl**2 / (2 * np.pi * m.L**2) * np.sqrt(ei / mass_per_len)
    return {
        "out1": f(ei_out, BETA_L[0]),
        "out2": f(ei_out, BETA_L[1]),
        "in1": f(ei_in, BETA_L[0]),
    }


def ccx_all(mesh, lam, ccx, workdir, name):
    inp = workdir / f"modal_{name}.inp"
    m.write_ccx_mesh(
        mesh,
        str(inp),
        properties={"plate": m.CompositeShellProperty(laminate=lam)},
        boundary_nodeset="clamped",
        solver_type="Modal",
        num_modes=m.N_SEARCH,
        quadratic=True,
    )
    res = m.run_ccx(inp, ccx)
    if res.returncode != 0:
        m.fail_ccx(res, inp)
    return np.sort(m.parse_ccx_frequencies(inp, n_modes=m.N_SEARCH))


ccx = m.ccx_bin_or_skip()
for name, angles in m.LAYUPS.items():
    lam = m._laminate(angles)
    runs = {}
    for nx, ny in MESHES:
        mesh = m._build_mesh(nx, ny)
        ae = np.sort(m._aero_frequencies(mesh, lam, m.N_SEARCH))
        with tempfile.TemporaryDirectory(dir="logs/check/l2-audit") as td:
            cc = ccx_all(mesh, lam, ccx, Path(td), name)
        runs[(nx, ny)] = (ae, cc)
    an = beam_anchors(angles)
    print(f"== {name}  beam anchors (Hz): out1 {an['out1']:.3f} out2 {an['out2']:.3f} in1 {an['in1']:.3f}")
    for k, (ae, cc) in runs.items():
        print(f"  {k[0]}x{k[1]:<3} aero " + " ".join(f"{v:8.3f}" for v in ae[:8]))
        print(f"  {'':6} ccx  " + " ".join(f"{v:8.3f}" for v in cc[:8]))
    ae_f, cc_f = runs[MESHES[-1]]
    ae_p, cc_p = runs[MESHES[-2]]
    cost = np.abs(ae_f[:, None] - cc_f[None, :]) / np.maximum(np.abs(cc_f[None, :]), 1e-14)
    ri, ci = linear_sum_assignment(cost)
    pairs = sorted(zip(ri, ci), key=lambda p: cost[p[0], p[1]])[: m.N_COMPARE]
    ok = True
    for i, j in pairs:
        d_ae = 100 * (ae_f[i] / ae_p[i] - 1)
        d_cc = 100 * (cc_f[j] / cc_p[j] - 1)
        gap = 100 * cost[i, j]
        conv = abs(d_ae) < 1 and abs(d_cc) < 1
        ok &= conv and gap <= 3
        print(
            f"  matched aero#{i} {ae_f[i]:8.3f} ccx#{j} {cc_f[j]:8.3f}  gap {gap:5.2f}%  "
            f"16->32 change aero {d_ae:+.3f}% ccx {d_cc:+.3f}%  {'converged' if conv else 'NOT converged'}"
        )
    print(f"  => {name}: comparison {'CLOSED' if ok else 'OPEN'} under the predeclared criterion", flush=True)
