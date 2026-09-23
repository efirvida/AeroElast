"""Debug: centrifugal K_G magnitudes (load, static solve, membrane forces, K_σ)."""
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

OMEGA = 0.7916808310986626
AXIS = np.array([0.0, 1.0, 0.0])


def main():
    import _aeroelast
    from petsc4py import PETSc

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=0.5, yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"))
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL,
                          "span_direction": (0.0, 0.0, 1.0),
                          "properties": props}}
    asm = MeshAssembler(mesh, model)
    n_dof = asm.dofs_count
    dpn = asm.dofs_per_node

    f_full = np.asarray(asm._rust.centrifugal_load(
        OMEGA, list(map(float, AXIS)), [0.0, 0.0, 0.0], asm._rho_per_elem))
    print(f"|f_cf| total = {np.abs(f_full).sum():.4e} N (F=0 expected)")
    f_nodes = f_full.reshape(-1, dpn)
    print(f"|f_cf| por nodo max = {np.abs(f_nodes).max():.4e} N")

    # analytic reference: w^2 * int_0^R r mu(r) dr
    xyz = np.asarray([n.coords for n in mesh.nodes])
    total_mass = 0.0
    # per-node mass from M diagonal
    m_rows, m_cols, m_vals = asm._rust.assemble_m()
    m_rows = np.asarray(m_rows)
    m_vals = np.asarray(m_vals)
    nodal_m = np.zeros(mesh.node_count)
    for r, v in zip(m_rows, m_vals):
        nodal_m[r // dpn] += v
    first_moment = (nodal_m * xyz[:, 2]).sum()
    print(f"blade mass = {nodal_m.sum():.4e} kg, first moment = {first_moment:.4e} kg·m")
    print(f"analytic root reaction = {OMEGA**2 * first_moment:.4e} N")

    # static solve
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    root_dofs = set()
    for n in root:
        for d in range(dpn):
            root_dofs.add(nid2idx[n] * dpn + d)
    free = np.array(sorted(set(range(n_dof)) - root_dofs), dtype=np.int64)

    k_rows, k_cols, k_vals = asm._rust.assemble_k()
    K_full = asm._coo_to_petsc(k_rows, k_cols, k_vals)
    is_free = PETSc.IS().createGeneral(free.astype(PETSc.IntType), comm=asm.comm)
    K_red = K_full.createSubMatrix(is_free, is_free)
    u_red = PETSc.Vec().createSeq(len(free), comm=asm.comm)
    f_red = PETSc.Vec().createSeq(len(free), comm=asm.comm)
    f_red.setArray(f_full[free])
    ksp = PETSc.KSP().create(asm.comm)
    ksp.setOperators(K_red)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setFromOptions()
    ksp.solve(f_red, u_red)
    u_full = np.zeros(n_dof)
    u_full[free] = u_red.getArray()
    u_tip = u_full.reshape(-1, dpn)[int(np.argmax(xyz[:, 2]))]
    print(f"u_tip (x,y,z) = {u_tip[:3]}  |u|={np.linalg.norm(u_tip[:3]):.4e} m")

    # membrane resultants
    (rows, cols, vals) = asm._rust.assemble_geometric_k_from_disp(u_full)
    print(f"K_G nnz = {len(vals)}, |max N| en K_G input…")

    # direct stress check via compute_stress_field
    sigma6, _ = asm._rust.compute_stress_field(u_full, 0.0, 0)
    sigma6 = np.asarray(sigma6)
    print(f"membrane σ (Pa): max={np.abs(sigma6[:, :3]).max():.3e}  mean|x|={np.abs(sigma6[:, 0]).mean():.3e}")

    # modal shift with uniform stress field (control): sigma=[1e6,0,0] N/m
    n_elems = len(asm.elements)
    s_arr = np.zeros((n_elems, 3))
    s_arr[:, 0] = 1.0e6
    r2, c2, v2 = asm._rust.assemble_geometric_k(s_arr)
    Kg_ctrl = coo_matrix((np.asarray(v2), (np.asarray(r2), np.asarray(c2))), shape=(n_dof, n_dof)).tocsr()

    K = coo_matrix((np.asarray(k_vals), (np.asarray(k_rows), np.asarray(k_cols))),
                   shape=(n_dof, n_dof)).tocsr()
    for name, Kg in [("K", coo_matrix((n_dof, n_dof))), ("K+Kg_ctrl(1e6)", Kg_ctrl)]:
        Ksum = (K + Kg).tocoo()
        freqs, _ = _aeroelast.modal_solve_coo(
            Ksum.row.astype(np.int64), Ksum.col.astype(np.int64), Ksum.data.astype(np.float64),
            np.asarray(m_rows, dtype=np.int64), np.asarray(m_cols, dtype=np.int64),
            np.asarray(m_vals, dtype=np.float64), n_dof, free, 6)
        print(f"{name}: f0..f5 = {np.round(np.sort(np.asarray(freqs)), 4)}")


if __name__ == "__main__":
    main()
