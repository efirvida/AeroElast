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


def assemble_rotating(assembler, mesh, omega, n_modes, include_kg=True, include_ksp=True):
    """Rotating-frame modal of K + K_G + K_SP; ``include_kg``/``include_ksp`` drop a term for diagnosis."""
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

    if omega > 0.0 and include_kg:
        kg = assembler.assemble_geometric_stiffness(
            omega=omega, rotation_axis=AXIS, rotation_center=np.zeros(3),
            free_dofs=free)
        ai, aj, av = kg.getValuesCSR()
        Kg = csr_matrix((av, aj, ai), shape=(n_dof, n_dof))
    else:
        Kg = csr_matrix((n_dof, n_dof))

    m_rows, m_cols, m_vals = assembler._rust.assemble_m()
    m_rows = np.asarray(m_rows)

    if omega > 0.0 and include_ksp:
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


# Below this fraction of the mode explained by rigid-section (beam) motion the mode is "local".
# Measured on the IEA-15 blade, es 0.25, oml and mid, 0 and 7.56 rpm (logs/check/modecls/,
# odd/tasks/s4-mode-classifier.md T1): below 3.5 Hz the beam modes 1F-3F/1E-2E sit at >= 0.83
# and the local modes at <= 0.44; 0.6 splits that gap.
BEAM_FRACTION_MIN = 0.6
# A family wins when its beam amplitude is at least this factor above every other family.
FAMILY_DOMINANCE = 2.0
SPAN_BINS = 200


def nodal_translations(mode, free, n_nodes, dpn=6):
    """Scatter a free-DOF mode vector into an (n_nodes, 3) translation array."""
    u = np.zeros((n_nodes, 3))
    node = free // dpn
    comp = free % dpn
    m = comp < 3
    u[node[m], comp[m]] = mode[m]
    return u


def beam_projection(xyz, u, n_bins=SPAN_BINS):
    """Rigid-section fit u_i = t + theta x (p_i - c) of nodal translations, binned along +z.

    Returns the beam fraction sum|u_fit|^2 / sum|u|^2 (unweighted by mass), the energy per
    beam family (edge t_x, flap t_y, axial t_z, torsion theta_z x r) and the residual.
    """
    z = xyz[:, 2]
    edges = np.linspace(z.min(), z.max() + 1e-9, n_bins + 1)
    which = np.clip(np.searchsorted(edges, z, side="right") - 1, 0, n_bins - 1)
    fit = np.zeros_like(u)
    energy = {"flap": 0.0, "edge": 0.0, "axial": 0.0, "torsion": 0.0}
    for b in range(n_bins):
        idx = np.flatnonzero(which == b)
        if idx.size < 3:
            continue
        r = xyz[idx] - xyz[idx].mean(axis=0)
        # u_i = t + theta x r_i  ->  rows [I | -skew(r_i)] [t; theta]
        a = np.zeros((3 * idx.size, 6))
        a[0::3, 0] = a[1::3, 1] = a[2::3, 2] = 1.0
        a[0::3, 4], a[0::3, 5] = r[:, 2], -r[:, 1]
        a[1::3, 3], a[1::3, 5] = -r[:, 2], r[:, 0]
        a[2::3, 3], a[2::3, 4] = r[:, 1], -r[:, 0]
        sol, *_ = np.linalg.lstsq(a, u[idx].reshape(-1), rcond=None)
        fit[idx] = (a @ sol).reshape(-1, 3)
        t, th = sol[:3], sol[3:]
        energy["edge"] += idx.size * t[0] ** 2
        energy["flap"] += idx.size * t[1] ** 2
        energy["axial"] += idx.size * t[2] ** 2
        energy["torsion"] += th[2] ** 2 * float(np.sum(r[:, 0] ** 2 + r[:, 1] ** 2))
    total = float(np.sum(u**2))
    frac = float(np.sum(fit**2)) / total if total > 0 else 0.0
    return frac, energy, u - fit


def classify_modes(mesh, modes, free, dpn=6):
    """Label each mode flap/edge/axial/torsion/mixed/local from the whole blade.

    Blade mesh convention: span along +z, chord along x (rotor plane), so
    out-of-plane flap motion = u_y and in-plane edge motion = u_x, with the
    rotation axis along y.

    Every mode is projected on beam kinematics (rigid cross-sections per span bin,
    ``beam_projection``). A mode the projection explains below ``BEAM_FRACTION_MIN`` is
    "local" (wall/panel motion, not a beam mode). Otherwise the family whose beam amplitude
    (square root of its span-summed energy) is ``FAMILY_DOMINANCE`` times every other wins,
    else "mixed". This replaced a label read from the single highest-z node, which let a
    low-amplitude local mode into the flap family (odd/tasks/s4-mode-classifier.md).
    """
    xyz = np.asarray([n.coords for n in mesh.nodes], dtype=np.float64)
    labels = []
    for k in range(modes.shape[0]):
        u = nodal_translations(modes[k], free, len(xyz), dpn)
        frac, energy, _ = beam_projection(xyz, u)
        if frac < BEAM_FRACTION_MIN:
            labels.append("local")
            continue
        amp = {fam: np.sqrt(e) for fam, e in energy.items()}
        best = max(amp, key=amp.get)
        others = [a for fam, a in amp.items() if fam != best]
        labels.append(best if amp[best] >= FAMILY_DOMINANCE * max(others) else "mixed")
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
