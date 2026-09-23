"""S-4: rotating blade modal (K + K_G + K_SP) vs OpenFAST MBC3 targets.

Computes the rotating-frame modal of the official WindIO blade (root clamped,
consistent mass) at the rated rotor speed and classifies each mode by its tip
motion: flap (u_x) vs edge (u_y) with span along +z, rotation about +y.

Comparison targets (OpenFAST v5.0.0 ElastoDyn, blade DOFs only, no aero,
RotSpeed=7.56 rpm, MBC3 rotating frame) from the S-4 linearization run:
  1st flap : 0.4405 (BW) / 0.5666 (collective) / 0.6926 (FW) Hz
  1st edge : 0.6200 (BW) / 0.7461 (collective) / 0.8721 (FW) Hz
  2nd flap : 1.5111 (BW) / 1.6372 (collective) / 1.7632 (FW) Hz

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s4_rotating_modal.py [--rpm 7.56] [--n-modes 20]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

RPM_DEFAULT = 7.56
AXIS = np.array([0.0, 1.0, 0.0])


def omega_from_rpm(rpm):
    return rpm * 2.0 * np.pi / 60.0


def assemble_rotating(assembler, mesh, omega, n_modes):
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

    if omega > 0.0:
        kg = assembler.assemble_geometric_stiffness(
            omega=omega, rotation_axis=AXIS, rotation_center=np.zeros(3),
            free_dofs=free)
        ai, aj, av = kg.getValuesCSR()
        Kg = csr_matrix((av, aj, ai), shape=(n_dof, n_dof))
    else:
        Kg = csr_matrix((n_dof, n_dof))

    m_rows, m_cols, m_vals = assembler._rust.assemble_m()
    m_rows = np.asarray(m_rows)

    if omega > 0.0:
        nodal_m = np.zeros(mesh.node_count)
        for i in range(mesh.node_count):
            nodal_m[i] = np.asarray(m_vals)[m_rows == i * dpn].sum()
        rows, cols, vals = [], [], []
        blk = -omega**2 * (np.eye(3) - np.outer(AXIS, AXIS))
        for i in range(mesh.node_count):
            for a in range(3):
                for b in range(3):
                    rows.append(i * dpn + a)
                    cols.append(i * dpn + b)
                    vals.append(nodal_m[i] * blk[a, b])
        Ksp = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()
    else:
        Ksp = csr_matrix((n_dof, n_dof))

    Ksum = (K + Kg + Ksp).tocoo()

    freqs, modes_flat = _aeroelast.modal_solve_coo(
        Ksum.row.astype(np.int64), Ksum.col.astype(np.int64),
        Ksum.data.astype(np.float64),
        m_rows.astype(np.int64), np.asarray(m_cols, dtype=np.int64),
        np.asarray(m_vals, dtype=np.float64),
        n_dof, free, n_modes,
    )
    freqs = np.asarray(freqs, dtype=np.float64)
    modes = np.asarray(modes_flat, dtype=np.float64).reshape(n_modes, len(free))
    return freqs, modes, free


def classify_modes(mesh, modes, free, dpn=6):
    """Classify each mode by tip-node translation direction (flap/edge/axial).

    Blade mesh convention: span along +z, chord along x (rotor plane), so
    out-of-plane flap motion = u_y and in-plane edge motion = u_x, with the
    rotation axis along y.
    """
    nid2idx = mesh.node_id_to_index
    xyz = np.asarray([n.coords for n in mesh.nodes])
    tip_idx = int(np.argmax(xyz[:, 2]))
    tip_dof0 = nid2idx[tip_idx] * dpn
    labels = []
    for k in range(modes.shape[0]):
        u = np.zeros(dpn)
        for d in range(3):
            dof = tip_dof0 + d
            pos = np.searchsorted(free, dof)
            if pos < len(free) and free[pos] == dof:
                u[d] = modes[k, pos]
        ux, uy, uz = abs(u[0]), abs(u[1]), abs(u[2])
        if uy >= 2.0 * ux and uy >= 2.0 * uz:
            labels.append("flap")
        elif ux >= 2.0 * uy and ux >= 2.0 * uz:
            labels.append("edge")
        elif uz >= 2.0 * ux and uz >= 2.0 * uy:
            labels.append("axial")
        else:
            labels.append("mixed")
    return labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpm", type=float, default=RPM_DEFAULT)
    ap.add_argument("--n-modes", type=int, default=20)
    ap.add_argument("--element-size", type=float, default=0.5)
    args = ap.parse_args()

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(element_size=args.element_size,
                    yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"))
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {"elements": {"element_family": ElementFamily.SHELL,
                          "span_direction": (0.0, 0.0, 1.0),
                          "properties": props}}
    asm = MeshAssembler(mesh, model)

    for rpm in (0.0, args.rpm):
        omega = omega_from_rpm(rpm)
        freqs, modes, free = assemble_rotating(asm, mesh, omega, args.n_modes)
        labels = classify_modes(mesh, modes, free)
        order = np.argsort(freqs)
        print(f"\nrpm={rpm} (omega={omega:.4f} rad/s), {args.n_modes} modos:")
        for i in order:
            print(f"  {i:2d}: f={freqs[i]:8.4f} Hz  {labels[i]}")


if __name__ == "__main__":
    main()
