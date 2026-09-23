"""S-6 — one-way DYNAMIC validation: AeroDyn load time series -> shell Newmark.

Extends S-5 to the time domain: the blade-1 AeroDyn distributed loads
(Fn/Ft per unit length, local chord frame) are sampled over a short window
of the OpenFAST blade-only rated run and applied to the official WindIO
shell blade with a Newmark-beta (trapezoidal) integration of the rotating
(K + K_G + K_SP) system, starting from the static equilibrium under the
window-mean loads.  The tip out-of-plane deflection time series is compared
against ElastoDyn's OoPDefl1 over the same window (mean + fluctuation).

Design notes:
  * The shell uses the twist-based section frames (WindIO sign, rotorspin
    fix) and the
    spanwise tributary weighting — identical to S-5's static mapping.
  * Damping: Rayleigh with the auto-computed campaign coefficients
    (eta_k = 7.73e-3 s on the stiffness, eta_m = 1.14e-1 1/s on the mass).
  * dt = 0.01 s (downsampled from the OpenFAST DT=0.005) — 185 samples per
    1st-flap cycle, well inside the Newmark accuracy band.
  * The K_G over-response documented in S-5 (~13% modal vs ~8.5% ElastoDyn)
    also affects the dynamic stiffness; the mean tolerance absorbs it.

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s6_oneway_dynamic.py --out <...>.out --ad-blade <...>.dat \
      --export tests/IEA15MW/reference/s6_ad_blade1_timeseries.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

RPM_RATED = 7.56
AXIS = np.array([0.0, 1.0, 0.0])
ETA_K = 7.73e-3   # s    (Rayleigh stiffness coefficient, campaign auto)
ETA_M = 1.14e-1   # 1/s  (Rayleigh mass coefficient, campaign auto)
DT = 0.01
WINDOW = (96.0, 100.0)


def extract_timeseries(out_path: Path, ad_blade_path: Path):
    """Sample the AD loads + the ElastoDyn tip channels over the window.

    Returns (t, fn, ft, twist, r_m, oop_ed):
      t      : (n_steps,) sampled at DT inside [t_lo, t_hi]
      fn/ft  : (n_steps, n_stations) N/m per unit length
      twist  : (n_stations,) geometric twist (deg) from the AD blade file
      r_m    : (n_stations,) spanwise stations (m)
      oop_ed : (n_steps,) ElastoDyn OoPDefl1 (m)
    """
    with open(out_path) as f:
        lines = f.readlines()
    names = lines[6].split()
    data = np.loadtxt(out_path, skiprows=8)
    t_all = data[:, 0]

    with open(ad_blade_path) as f:
        ad_lines = f.readlines()
    n_nodes = int(ad_lines[3].split()[0])
    blspn, twist = [], []
    for ln in ad_lines[6:6 + n_nodes]:
        parts = ln.split()
        if len(parts) >= 5:
            blspn.append(float(parts[0]))
            twist.append(float(parts[4]))
    blspn = np.asarray(blspn)
    twist = np.asarray(twist)

    t_lo, t_hi = WINDOW
    sel = (t_all >= t_lo) & (t_all <= t_hi)
    t = t_all[sel]
    keep = np.isclose(np.mod(t, DT), 0.0, atol=1e-9) | np.isclose(
        np.mod(t, DT), DT, atol=1e-9)
    t = t[keep]

    fn = np.zeros((len(t), n_nodes))
    ft = np.zeros((len(t), n_nodes))
    for n in range(1, n_nodes + 1):
        fn[:, n - 1] = data[sel, names.index(f"AB1N{n:03d}Fn")][keep]
        ft[:, n - 1] = data[sel, names.index(f"AB1N{n:03d}Ft")][keep]
    oop = data[sel, names.index("OoPDefl1")][keep]
    return t, fn, ft, twist, blspn, oop


def build_assembly(element_size=0.5):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(
        element_size=element_size,
        yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"),
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    model = {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": props,
        }
    }
    return MeshAssembler(mesh, model), mesh


def assemble_keff(assembler, mesh, omega):
    """K_eff = K + K_G + K_SP, M, and the clamped-root free DOFs."""
    from scipy.sparse import coo_matrix, csr_matrix

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
        omega=omega, rotation_axis=AXIS, rotation_center=np.zeros(3),
        free_dofs=free)
    ai, aj, av = kg.getValuesCSR()
    Kg = csr_matrix((av, aj, ai), shape=(n_dof, n_dof))

    m_rows, m_cols, m_vals = assembler._rust.assemble_m()
    nodal_m = np.zeros(mesh.node_count)
    for r, v in zip(np.asarray(m_rows), np.asarray(m_vals)):
        if r % dpn < 3:
            nodal_m[r // dpn] += v / 3.0
    rows, cols, vals = [], [], []
    blk = -omega**2 * (np.eye(3) - np.outer(AXIS, AXIS))
    for i in range(mesh.node_count):
        for a in range(3):
            for b in range(3):
                rows.append(i * dpn + a)
                cols.append(i * dpn + b)
                vals.append(nodal_m[i] * blk[a, b])
    Ksp = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()

    M = coo_matrix((np.asarray(m_vals), (np.asarray(m_rows), np.asarray(m_cols))),
                   shape=(n_dof, n_dof)).tocsr()

    return (K + Kg + Ksp).tocsr(), M, free


def nodal_load_series(assembler, mesh, r_m, fn_ts, ft_ts, twist_deg):
    """(n_steps, n_dof) nodal force matrix; frames/tributaries precomputed."""
    n_dof = assembler.dofs_count
    dpn = assembler.dofs_per_node
    xyz = np.asarray([n.coords for n in mesh.nodes])
    z = xyz[:, 2]
    n_nodes = mesh.node_count
    n_steps = fn_ts.shape[0]

    # per-station spanwise tributary split over the station's nodes
    # (fixes the 1-2 nodes per section artefact; see load_mapping.py)
    from load_mapping import spanwise_tributary

    trib = spanwise_tributary(z)

    theta = np.deg2rad(np.interp(z, r_m, twist_deg))
    n_hat = np.column_stack([-np.sin(theta), np.cos(theta), np.zeros_like(theta)])
    t_hat = np.column_stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)])

    fn_i = np.array([np.interp(z, r_m, row) for row in fn_ts])   # (n_steps, n_nodes)
    ft_i = np.array([np.interp(z, r_m, row) for row in ft_ts])

    F = np.zeros((n_steps, n_dof))
    for a in range(3):
        F[:, a::dpn] = (fn_i * n_hat[None, :, a]
                        + ft_i * t_hat[None, :, a]) * trib[None, :]
    return F


def newmark_tip_series(assembler, mesh, K_eff, M, free, F_ts, dt):
    """Trapezoidal Newmark from the static IC; returns (t_tip_xyz, steps)."""
    from petsc4py import PETSc

    n_dof = assembler.dofs_count
    dpn = assembler.dofs_per_node
    n_steps = F_ts.shape[0]

    C = ETA_K * K_eff + ETA_M * M
    # trapezoidal rule (beta=0.25, gamma=0.5) integration constants
    a0 = 4.0 / dt**2          # 1/(beta dt^2)
    a1 = 2.0 / dt             # gamma/(beta dt)
    a2 = 4.0 / dt             # 1/(beta dt)
    a3 = 1.0                  # 1/(2 beta) - 1
    a4 = 1.0                  # gamma/beta - 1
    Keff = a0 * M + a1 * C + K_eff

    Kp = assembler._coo_to_petsc(
        np.asarray(Keff.tocoo().row, dtype=np.int64),
        np.asarray(Keff.tocoo().col, dtype=np.int64),
        np.asarray(Keff.tocoo().data, dtype=np.float64))
    is_free = PETSc.IS().createGeneral(free.astype(PETSc.IntType), comm=assembler.comm)
    K_red = Kp.createSubMatrix(is_free, is_free)
    ksp = PETSc.KSP().create(assembler.comm)
    ksp.setOperators(K_red)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setFromOptions()

    # static equilibrium under the mean load -> IC (STATIC K_eff solve,
    # NOT the dynamic effective matrix)
    f_mean = F_ts.mean(axis=0)
    rhs = PETSc.Vec().createSeq(len(free), comm=assembler.comm)
    sol = PETSc.Vec().createSeq(len(free), comm=assembler.comm)
    Kp_static = assembler._coo_to_petsc(
        np.asarray(K_eff.tocoo().row, dtype=np.int64),
        np.asarray(K_eff.tocoo().col, dtype=np.int64),
        np.asarray(K_eff.tocoo().data, dtype=np.float64))
    K_red_static = Kp_static.createSubMatrix(is_free, is_free)
    ksp_static = PETSc.KSP().create(assembler.comm)
    ksp_static.setOperators(K_red_static)
    ksp_static.setType("preonly")
    pc_s = ksp_static.getPC()
    pc_s.setType("lu")
    pc_s.setFactorSolverType("mumps")
    ksp_static.setFromOptions()
    rhs.setArray(f_mean[free])
    ksp_static.solve(rhs, sol)
    u = np.zeros(n_dof)
    u[free] = sol.getArray()
    v = np.zeros(n_dof)
    a = np.zeros(n_dof)

    coords = np.asarray([n.coords for n in mesh.nodes])
    tip_idx = int(np.argmax(coords[:, 2]))
    tip_dof = tip_idx * dpn

    # M/C reduced row access via scipy slicing on the free partition
    M_ff = M[free][:, free].tocsr()
    C_ff = C[free][:, free].tocsr()

    hist = np.zeros((n_steps, 3))
    u_new = np.zeros(n_dof)
    v_new = np.zeros(n_dof)
    a_new = np.zeros(n_dof)
    rhs = PETSc.Vec().createSeq(len(free), comm=assembler.comm)
    sol = PETSc.Vec().createSeq(len(free), comm=assembler.comm)

    for step in range(n_steps):
        f = F_ts[step][free]
        eff = f \
            + (a0 * (M_ff @ u[free]) + a2 * (M_ff @ v[free]) + a3 * (M_ff @ a[free])) \
            + (a1 * (C_ff @ u[free]) + a4 * (C_ff @ v[free]))
        rhs.setArray(eff)
        ksp.solve(rhs, sol)
        u_new[free] = sol.getArray()
        a_new[free] = a0 * (u_new[free] - u[free]) - a2 * v[free] - a3 * a[free]
        v_new[free] = v[free] + 0.5 * dt * (a[free] + a_new[free])

        u, v, a = u_new.copy(), v_new.copy(), a_new.copy()
        hist[step] = u[tip_dof:tip_dof + 3]

    return hist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/"
                                 "IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.out"))
    ap.add_argument("--ad-blade", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/"
                                 "IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat"))
    ap.add_argument("--export", type=Path,
                    default=Path("tests/IEA15MW/reference/s6_ad_blade1_timeseries.csv"))
    ap.add_argument("--element-size", type=float, default=0.5)
    args = ap.parse_args()

    t, fn_ts, ft_ts, twist, r_m, oop_ed = extract_timeseries(args.out, args.ad_blade)
    print(f"window: {t[0]:.2f}..{t[-1]:.2f} s, {len(t)} steps, "
          f"mean Fn={fn_ts.mean():.4e} N/m")

    if args.export:
        header = "t,OoPDefl1," + ",".join(
            [f"AB1N{n:03d}Fn" for n in range(1, 51)]
            + [f"AB1N{n:03d}Ft" for n in range(1, 51)])
        np.savetxt(args.export, np.column_stack(
            [t, oop_ed, fn_ts, ft_ts]), header=header, delimiter=",", comments="")
        print(f"time series exported to {args.export}")

    asm, mesh = build_assembly(args.element_size)
    K_eff, M, free = assemble_keff(asm, mesh, RPM_RATED * 2 * np.pi / 60.0)
    F_ts = nodal_load_series(asm, mesh, r_m, fn_ts, ft_ts, twist)
    hist = newmark_tip_series(asm, mesh, K_eff, M, free, F_ts, DT)

    oop_shell = hist[:, 1]
    print(f"\nS-6 shell series (first 8): {oop_shell[:8].round(3)} m")
    print(f"S-6 shell series (last 4):  {oop_shell[-4:].round(3)} m")
    print(f"S-6 shell min={oop_shell.min():.3f} max={oop_shell.max():.3f} m")
    print(f"\nS-6 comparison over the window ({len(t)} steps):")
    print(f"  OoP mean: shell {oop_shell.mean():+.4f} m vs ElastoDyn "
          f"{oop_ed.mean():+.4f} m "
          f"-> {(abs(oop_shell.mean()) - abs(oop_ed.mean())) / abs(oop_ed.mean()):+.1%}")
    print(f"  OoP std : shell {oop_shell.std():.4f} m vs ElastoDyn {oop_ed.std():.4f} m "
          f"-> {(oop_shell.std() - oop_ed.std()) / oop_ed.std():+.1%}")


if __name__ == "__main__":
    main()
