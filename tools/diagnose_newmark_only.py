"""Minimal Newmark stability test: yaml vs xlsx blade, no FSI, no rotation.

Assembles K (Rust COO), lumped M, Rayleigh C from the SLEPc modes (the same
calibration the rotor solver uses), clamps the root, and runs the plain
Newmark-beta (beta=0.25, gamma=0.5, dt=0.01) under the gravity body force.

If the yaml blade diverges here, the instability is purely structural
(K/M/C assembly), with no rotor machinery involved.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import spsolve

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

G = np.array([0.0, 0.0, -9.81])
DT = 0.01
BETA, GAMMA = 0.25, 0.5
ZETA = 0.03
N_STEPS = 150


def run(name, kw):
    import _aeroelast
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=0.5, **kw)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL,
                          "span_direction": (0.0, 0.0, 1.0),
                          "properties": props}}
    asm = MeshAssembler(mesh, model)
    n_dof = asm.dofs_count
    dpn = asm.dofs_per_node

    kr, kc, kv = asm._rust.assemble_k()
    K = coo_matrix((np.asarray(kv), (np.asarray(kr), np.asarray(kc))),
                   shape=(n_dof, n_dof)).tocsr()
    mr, mc, mv = asm._rust.assemble_m()
    M = coo_matrix((np.asarray(mv), (np.asarray(mr), np.asarray(mc))),
                   shape=(n_dof, n_dof)).tocsr()
    m_diag = np.asarray(M.sum(axis=1)).ravel()

    # Rayleigh calibration: the rotor solver uses the Rust SLEPc modes 1..2
    # (root clamped, like V-02)
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    root_dofs = {nid2idx[n] * dpn + d for n in root for d in range(dpn)}
    free = np.array(sorted(set(range(n_dof)) - root_dofs), dtype=np.int64)
    freqs, _ = _aeroelast.modal_solve_coo(
        np.asarray(kr, dtype=np.int64), np.asarray(kc, dtype=np.int64),
        np.asarray(kv, dtype=np.float64),
        np.asarray(mr, dtype=np.int64), np.asarray(mc, dtype=np.int64),
        np.asarray(mv, dtype=np.float64), n_dof, free, 3)
    f = np.sort(np.asarray(freqs, dtype=np.float64))
    w1, w2 = 2 * np.pi * f[0], 2 * np.pi * f[1]
    A = np.array([[1 / (2 * w1), w1 / 2], [1 / (2 * w2), w2 / 2]])
    eta_m, eta_k = np.linalg.solve(A, [ZETA, ZETA])
    C = eta_m * diags(m_diag) + eta_k * K

    # gravity body force: consistent nodal loads from the mass matrix
    F = np.zeros(n_dof)
    for i in range(mesh.node_count):
        for d in range(3):
            F[i * dpn + d] = m_diag[i * dpn + d] * G[d] if m_diag[i * dpn + d] > 0 else 0.0

    Kf = K[free][:, free].tocsr()
    Mf = M[free][:, free].tocsr()
    Cf = C[free][:, free].tocsr()

    a0 = 1.0 / (BETA * DT**2)
    a1 = GAMMA / (BETA * DT)
    a2 = 1.0 / (BETA * DT)
    a3 = 1.0 / (2 * BETA) - 1.0
    a4 = GAMMA / BETA - 1.0
    a5 = DT / 2 * (GAMMA / BETA - 2.0)

    Keff = Kf + a0 * Mf + a1 * Cf
    Keff = Keff.tocsc()
    u = np.zeros(len(free))
    v = np.zeros(len(free))
    acc = np.zeros(len(free))
    Ff = F[free]

    tip_nodes = np.nonzero(mesh.coords_array[:, 2] >= mesh.coords_array[:, 2].max() - 1e-3)[0]
    tip_free = [nid2idx[n] * dpn + 1 for n in tip_nodes if nid2idx[n] * dpn + 1 in free]

    history = []
    for step in range(N_STEPS):
        Feff = Ff + Mf @ (a0 * u + a2 * v + a3 * acc) + Cf @ (a1 * u + a4 * v + a5 * acc)
        u_new = spsolve(Keff, Feff)
        acc_new = a0 * (u_new - u) - a2 * v - a3 * acc
        v_new = v + DT * ((1 - GAMMA) * acc + GAMMA * acc_new)
        u, v, acc = u_new, v_new, acc_new
        if step % 10 == 0:
            history.append(float(np.mean(u[tip_free]) if tip_free else np.nan))
    print(f"{name}: tip uy(t) = {np.round(history, 3)}")
    return np.asarray(history)


if __name__ == "__main__":
    h_y = run("yaml", dict(yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml")))
    h_x = run("xlsx", dict(excel_file=str(REPO / "tests/NuMAD_utd_iea15mw.xlsx"),
                           airfoil_dir=str(REPO / "tests/airfoils")))
    print("max |uy|: yaml=%.3g  xlsx=%.3g" % (np.max(np.abs(h_y)), np.max(np.abs(h_x))))
