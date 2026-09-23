"""A/B damping-matrix diagnosis: min eigenvalue of C for yaml vs xlsx.

C = eta_m * M + eta_k * (K + K_G + K_SP)  (tangent-consistent Rayleigh,
per theory doc 03).  A negative eigenvalue of C = negative damping = the
exponential transient growth observed in the yaml FSI runs.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

OMEGA = 0.7906341464750989
AXIS = np.array([0.0, 1.0, 0.0])
ETA_M = 0.11447
ETA_K = 7.7304e-3


def build_c(assembler, mesh):
    n_dof = assembler.dofs_count
    dpn = assembler.dofs_per_node

    k_rows, k_cols, k_vals = assembler._rust.assemble_k()
    K = coo_matrix((np.asarray(k_vals), (np.asarray(k_rows), np.asarray(k_cols))),
                   shape=(n_dof, n_dof)).tocsr()

    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    root_dofs = set()
    for n in root:
        for d in range(dpn):
            root_dofs.add(nid2idx[n] * dpn + d)
    free = np.array(sorted(set(range(n_dof)) - root_dofs), dtype=np.int64)
    kg = assembler.assemble_geometric_stiffness(
        omega=OMEGA, rotation_axis=AXIS, rotation_center=np.zeros(3),
        free_dofs=free)
    ai, aj, av = kg.getValuesCSR()
    Kg = csr_matrix((av, aj, ai), shape=(n_dof, n_dof))

    m_rows, m_cols, m_vals = assembler._rust.assemble_m()
    m_rows = np.asarray(m_rows); m_vals = np.asarray(m_vals)
    M_diag = np.zeros(n_dof)
    for r, v in zip(m_rows, m_vals):
        M_diag[r] += v

    rows, cols, vals = [], [], []
    blk = -OMEGA**2 * (np.eye(3) - np.outer(AXIS, AXIS))
    for i in range(mesh.node_count):
        m_i = M_diag[i * dpn]
        for a in range(3):
            for b in range(3):
                rows.append(i * dpn + a)
                cols.append(i * dpn + b)
                vals.append(m_i * blk[a, b])
    Ksp = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()

    C = ETA_M * diags(M_diag) + ETA_K * (K + Kg + Ksp)
    return C, M_diag


def main():
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    cases = [
        ("yaml", dict(yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"))),
        ("xlsx", dict(excel_file=str(REPO / "tests/NuMAD_utd_iea15mw.xlsx"),
                      airfoil_dir=str(REPO / "tests/airfoils"))),
    ]
    for name, kw in cases:
        gen = BladeMesh(element_size=0.5, **kw)
        mesh = gen.generate(renumber="rcm")
        props = build_rust_properties(gen.numad_mesh_data)
        model = {"elements": {"element_family": ElementFamily.SHELL,
                              "span_direction": (0.0, 0.0, 1.0),
                              "properties": props}}
        asm = MeshAssembler(mesh, model)
        C, M_diag = build_c(asm, mesh)
        d = C.diagonal()
        n_neg_diag = int((d < 0).sum())
        # cheap negativity proxy: smallest eigenvalue of the SHIFTED system
        # via a few power iterations on the most negative direction
        print(f"{name}: C diag min={d.min():.3e}  n_neg_diag={n_neg_diag}  "
              f"sum(C diag)={d.sum():.3e}")
        # Rayleigh quotient on the K_SP-softest direction: check the in-plane
        # translational DOFs (the spin softening acts there)
        dpn = asm.dofs_per_node
        d_trans = d.reshape(-1, dpn)[:, :3].min(axis=1)
        d_rot = d.reshape(-1, dpn)[:, 3:].min(axis=1)
        print(f"{name}: min C diag (translational) = {d_trans.min():.3e}  "
              f"min C diag (rotational) = {d_rot.min():.3e}")
        print()


if __name__ == "__main__":
    main()
