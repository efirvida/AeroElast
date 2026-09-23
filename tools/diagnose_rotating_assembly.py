"""A/B rotating-assembly diagnosis: yaml vs xlsx blade (solid side).

Compares the rotating-frame effective stiffness (K + K_G + K_SP) of both
blade models via the Rust modal solver (root clamped, like V-02):

  - K_G : centrifugal geometric stiffness (Rust assembler, Python call)
  - K_SP: spin softening -w^2 m_i (I - n n^T) on translational DOFs

A collapsed/near-zero rotating modal of the yaml blade would explain the
FSI divergence (10x transient at t=0.06 s + Newmark collapse).

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/diagnose_rotating_assembly.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

OMEGA = 0.7906341464750989
AXIS = np.array([0.0, 1.0, 0.0])


def rotating_modes(assembler, mesh):
    import _aeroelast

    n_dof = assembler.dofs_count
    dpn = assembler.dofs_per_node

    k_rows, k_cols, k_vals = assembler._rust.assemble_k()
    K = coo_matrix((np.asarray(k_vals), (np.asarray(k_rows), np.asarray(k_cols))),
                   shape=(n_dof, n_dof)).tocsr()

    # free DOFs: everything except the clamped root nodes
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

    # K_SP on translational DOFs: -w^2 m_i (I - n n^T)
    m_rows, m_cols, m_vals = assembler._rust.assemble_m()
    m_rows = np.asarray(m_rows)
    nodal_m = np.zeros(mesh.node_count)
    for i in range(mesh.node_count):
        nodal_m[i] = np.asarray(m_vals)[m_rows == i * dpn].sum()
    rows, cols, vals = [], [], []
    blk = -OMEGA**2 * (np.eye(3) - np.outer(AXIS, AXIS))
    for i in range(mesh.node_count):
        for a in range(3):
            for b in range(3):
                rows.append(i * dpn + a)
                cols.append(i * dpn + b)
                vals.append(nodal_m[i] * blk[a, b])
    Ksp = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()

    Ksum = (K + Kg + Ksp).tocoo()

    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    root_dofs = set()
    for n in root:
        for d in range(dpn):
            root_dofs.add(nid2idx[n] * dpn + d)
    free = np.array(sorted(set(range(n_dof)) - root_dofs), dtype=np.int64)

    freqs, _ = _aeroelast.modal_solve_coo(
        Ksum.row.astype(np.int64), Ksum.col.astype(np.int64),
        Ksum.data.astype(np.float64),
        np.asarray(m_rows, dtype=np.int64), np.asarray(m_cols, dtype=np.int64),
        np.asarray(m_vals, dtype=np.float64),
        n_dof, free, 8,
    )
    return np.sort(np.asarray(freqs, dtype=np.float64))


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
    results = {}
    for name, kw in cases:
        gen = BladeMesh(element_size=0.5, **kw)
        mesh = gen.generate(renumber="rcm")
        props = build_rust_properties(gen.numad_mesh_data)
        model = {"elements": {"element_family": ElementFamily.SHELL,
                              "span_direction": (0.0, 0.0, 1.0),
                              "properties": props}}
        asm = MeshAssembler(mesh, model)
        freqs = rotating_modes(asm, mesh)
        results[name] = freqs
        print(f"{name}: rotating freqs [Hz] = {np.round(freqs, 3)}")

    print("\ncomparison:")
    for i in range(8):
        print(f"  mode {i}: yaml={results['yaml'][i]:.3f}  xlsx={results['xlsx'][i]:.3f}")


if __name__ == "__main__":
    main()
