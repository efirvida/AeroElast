"""S-3 — Newmark time-integration validation against analytic solutions.

Runs the same trapezoidal Newmark formulation used by the S-6 tool and the
FSI rotor loop on the thin-walled closed box (whose EI/GJ are already
validated analytically in test_box_torsion_bending_benchmark.py) and checks:

  1. Undamped free-vibration frequency vs the Euler-Bernoulli first mode.
  2. Energy conservation of the trapezoidal scheme (undamped, no loads).
  3. Rayleigh-damped decay ratio vs the prescribed damping coefficients.

Analytic reference (Euler-Bernoulli cantilever, mass per unit length m):
    f1 = 1.875^2 / (2 pi L^2) * sqrt(E I / m)
The box: I_x exact = w t h^2/2 + t h^3/6 (validated), m = 4 w t rho.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve, splu  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from test_box_torsion_bending_benchmark import (  # noqa: E402
    DOF,
    E,
    H,
    L,
    NU,
    W,
    _assemble,
    _box_mesh,
    _find_node,
    I_X_EXACT,
)

RHO = 7800.0            # steel-like density so the frequencies are comfortable
M_PER_LEN = 4.0 * W * 0.01 * RHO   # 4 walls * w * t * rho
F1_EXACT = 1.875**2 / (2 * np.pi * L**2) * np.sqrt(E * I_X_EXACT / M_PER_LEN)
ETA_K = 1.0e-4          # stiffness-proportional Rayleigh coefficient (s)
ETA_M = 0.0             # mass-proportional (1/s)
ZETA_1_TARGET = ETA_K * 2 * np.pi * F1_EXACT / 2.0   # zeta = eta_k * omega_1 / 2

DT = 1.0e-3
N_STEPS = 1500          # 1.5 s -> ~10 cycles of the first mode


def _fixed_dofs(mesh, m):
    fixed = []
    for n in mesh.nodes:
        if n.coords[2] < 1e-12:
            fixed.extend(m[n.id] * DOF + d for d in range(DOF))
    return np.array(fixed)


def _newmark_transient(K, M, C, u0, v0, dt, n_steps, fixed, hist_dof):
    """Trapezoidal Newmark (beta=0.25, gamma=0.5) — same coefficients as the
    S-6 tool / FSI loop.  Returns (hist, E_hist)."""
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    a0 = 4.0 / dt**2
    a1 = 2.0 / dt
    a2 = 4.0 / dt
    a3 = 1.0
    a4 = 1.0
    Keff = a0 * M + a1 * C + K
    Kf = splu(Keff[free][:, free].tocsc())  # factorize once

    u = u0.copy()
    v = v0.copy()
    Mff = M[free][:, free].tocsc()
    Cff = C[free][:, free].tocsc()
    a = np.zeros_like(u)
    a[free] = spsolve(Mff, -(K[free][:, free] @ u[free]) - Cff @ v[free])

    hist = np.zeros(n_steps + 1)
    hist[0] = u[hist_dof]
    E_hist = np.zeros(n_steps + 1)
    E_hist[0] = 0.5 * (u @ (K @ u)) + 0.5 * (v @ (M @ v))
    for step in range(n_steps):
        eff = a0 * (Mff @ u[free]) + a2 * (Mff @ v[free]) + a3 * (Mff @ a[free]) \
            + a1 * (Cff @ u[free]) + a4 * (Cff @ v[free])
        u_new = u.copy()
        u_new[free] = Kf.solve(eff)
        a_new = a.copy()
        a_new[free] = a0 * (u_new[free] - u[free]) - a2 * v[free] - a3 * a[free]
        v_new = v.copy()
        v_new[free] = v[free] + 0.5 * dt * (a[free] + a_new[free])
        u, v, a = u_new, v_new, a_new
        hist[step + 1] = u[hist_dof]
        E_hist[step + 1] = 0.5 * (u @ (K @ u)) + 0.5 * (v @ (M @ v))
    return hist, E_hist


pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def transient_case():
    mesh = _box_mesh(nx=8, ny=40)
    K, m = _assemble(mesh, rho=RHO)

    from _aeroelast import PyMeshAssembler  # type: ignore[import-not-found]
    nodes_sorted = sorted(mesh.nodes, key=lambda n: n.id)
    conn = [[m[n.id] for n in e.nodes] for e in mesh.elements]
    mat = {"type": "isotropic", "e": E, "nu": NU, "rho": RHO, "thickness": 0.01,
           "shear_correction": 5.0 / 6.0}
    asm = PyMeshAssembler(
        node_coords=np.asarray([n.coords[:3] for n in nodes_sorted], dtype=float),
        connectivity=conn,
        elem_types=[4] * len(conn),
        materials=[dict(mat) for _ in conn])
    mr, mc, mv = asm.assemble_m()
    M = coo_matrix((mv, (mr, mc)), shape=(K.shape[0], K.shape[0])).tocsr()

    fixed = _fixed_dofs(mesh, m)
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    tip_dof = m[_find_node(mesh, 0.0, H / 2, L).id] * DOF + 1

    # 1st bending mode (numeric reference — the pure-mode IC for the transient)
    from scipy.sparse.linalg import eigsh
    lam, vecs = eigsh(K[free][:, free], k=6, M=M[free][:, free],
                      sigma=1.0, which="LM")
    order = np.argsort(lam)
    # pick the lowest mode whose tip y-motion dominates (bending about x)
    for k in order:
        v = vecs[:, k]
        tip_pos = np.where(free == tip_dof)[0][0]
        if abs(v[tip_pos]) > 0.1 * np.abs(v).max():
            f1_num = np.sqrt(lam[k]) / (2 * np.pi)
            u_modal = np.zeros(K.shape[0])
            u_modal[free] = v
            break
    else:
        raise RuntimeError("no bending mode found")

    return {"mesh": mesh, "K": K, "M": M, "m": m, "fixed": fixed,
            "tip_dof": tip_dof, "f1_num": float(f1_num),
            "u_modal": u_modal}


def test_undamped_free_vibration_frequency(transient_case):
    """Released first-mode shape must oscillate at the analytic f1 (EB cantilever)."""
    K, M, fixed, tip_dof = (transient_case[k] for k in
                            ("K", "M", "fixed", "tip_dof"))
    u0 = transient_case["u_modal"]
    v0 = np.zeros_like(u0)

    hist, _ = _newmark_transient(K, M, 0.0 * M, u0, v0, DT, N_STEPS, fixed, tip_dof)

    s = hist - hist.mean()
    crossings = np.where(np.diff(np.sign(s)))[0]
    f_num = (len(crossings) / 2.0) / (N_STEPS * DT)

    rel = abs(f_num - F1_EXACT) / F1_EXACT
    assert rel < 0.05, (
        f"free-vibration frequency {f_num:.3f} Hz vs analytic {F1_EXACT:.3f} Hz "
        f"(rel {rel:.2%} > 5%; the shell's 1st mode is at "
        f"{transient_case['f1_num']:.3f} Hz)"
    )


def test_energy_conservation_undamped(transient_case):
    """The trapezoidal scheme must conserve energy in the undamped case."""
    K, M, fixed, tip_dof = (transient_case[k] for k in
                            ("K", "M", "fixed", "tip_dof"))
    u0 = transient_case["u_modal"]
    v0 = np.zeros_like(u0)

    _, E_hist = _newmark_transient(K, M, 0.0 * M, u0, v0, DT, N_STEPS, fixed, tip_dof)

    drift = np.abs(E_hist[1:] / E_hist[0] - 1.0).max()
    assert drift < 1e-6, (
        f"total-energy relative drift {drift:.2e} > 1e-6 (trapezoidal must conserve)"
    )


def test_rayleigh_damping_decay_ratio(transient_case):
    """The decay envelope must match the Rayleigh damping at the first mode."""
    K, M, fixed, tip_dof = (transient_case[k] for k in
                            ("K", "M", "fixed", "tip_dof"))
    C = ETA_K * K
    u0 = transient_case["u_modal"]
    v0 = np.zeros_like(u0)

    hist, _ = _newmark_transient(K, M, C, u0, v0, DT, N_STEPS, fixed, tip_dof)

    f1 = transient_case["f1_num"]
    zeta_target = ETA_K * 2 * np.pi * f1 / 2.0

    env = np.abs(hist)
    t = np.arange(len(hist)) * DT
    pk_idx = np.where((np.diff(np.sign(np.diff(env))) < 0))[0] + 1
    pk_idx = pk_idx[pk_idx < len(hist)][:12]
    peaks, pvals = t[pk_idx], env[pk_idx]
    if len(peaks) < 4:
        pytest.skip("not enough peaks for the decay fit")
    zeta_num = -np.polyfit(peaks, np.log(pvals), 1)[0] / (2 * np.pi * f1)

    rel = abs(zeta_num - zeta_target) / zeta_target
    assert rel < 0.10, (
        f"decay zeta {zeta_num:.5f} vs target {zeta_target:.5f} (rel {rel:.2%} > 10%)"
    )
