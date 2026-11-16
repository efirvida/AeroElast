"""T6 verdict: settled mean of the rated one-way run vs the like-for-like static solve (odd/tasks/l5-one-way.md T6).

The T6 run (``cases/l5_oneway_rated``) carries the aero load and the centrifugal load on a constant
K_eff = K + K_G + K_SP (no gravity, no structural precone, K_G never updated). Its time mean over whole
revolutions must equal the static solution K_eff u = <f_aero> + f_cf; the T3 reference (t11/bem.json) is
aero only, so the centrifugal response is solved separately and reported as its own column.

Inputs: the ``results/<t>/{fields.vtu,state.npz}`` dumps of the run. Forces in the VTU are inertial-frame
(rotor.py ``to_inertial``) and are rotated back with the dump's ``theta``; the frame of ``U`` is detected
from its 1P content. Metrics reuse ``l5_oneway_bd_parity.shell_sections`` (ring affine fit at the
reference axis, equal-weight and skin-weight twist).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import l5_oneway_bd_parity as par  # noqa: E402
import twist_p1 as tp  # noqa: E402


def load_dumps(results: Path, n_rev: int, period: float):
    import meshio

    dirs = sorted((p for p in results.iterdir() if p.is_dir() and (p / "fields.vtu").exists()),
                  key=lambda p: float(p.name))
    t_end = float(dirs[-1].name)
    t0 = t_end - n_rev * period
    sel = [p for p in dirs if float(p.name) >= t0 - 1e-9]
    T, TH, U, FA, FI, X = [], [], [], [], [], None
    for p in sel:
        m = meshio.read(p / "fields.vtu")
        s = np.load(p / "state.npz")
        X = m.points
        T.append(float(s["t"])), TH.append(float(s["theta"]))
        U.append(np.asarray(m.point_data["U"])[:, :3])
        FA.append(np.asarray(m.point_data["F_AERO"])[:, :3])
        FI.append(np.asarray(m.point_data["F_INERTIAL"])[:, :3])
    return np.array(T), np.array(TH), np.array(U), np.array(FA), np.array(FI), X


def tmean(t, a):
    """Trapezoid time mean over the sampled window."""
    w = np.zeros_like(t)
    dt = np.diff(t)
    w[:-1] += dt / 2
    w[1:] += dt / 2
    return np.tensordot(w / w.sum(), a, axes=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, required=True, help="T6 results/ directory")
    ap.add_argument("--reference", type=Path, default=REPO / "logs/check/l5_oneway/t11/bem.json")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-rev", type=int, default=8, help="whole revolutions averaged, ending at the last dump")
    ap.add_argument("--hub-radius", type=float, default=3.97)
    args = ap.parse_args()
    t_start = time.time()
    from aeroelast.solvers.fsi.corotational import CoordinateTransforms, InertialForcesCalculator
    from run_s5_oneway import solve_static

    omega = par.OMEGA
    period = 2 * np.pi / omega
    t, th, U, FA, FI, Xv = load_dumps(args.results, args.n_rev, period)
    center = np.array([0.0, 0.0, -args.hub_radius])
    ct = CoordinateTransforms(rotation_axis=par.AXIS, rotation_center=center)
    FA_r = np.array([ct.to_rotating(f, a) for f, a in zip(FA, th, strict=True)])
    FI_r = np.array([ct.to_rotating(f, a) for f, a in zip(FI, th, strict=True)])

    mesh, props, asm = par.build_shell(20.0)
    dpn = asm.dofs_per_node
    f_bem, part = par.shell_loads("bem", asm, mesh, props, -6.0, -4.0)
    X = np.asarray(part._ref_coords, float)
    # the VTU writes the deformed configuration: points = X0 + U of that dump (same node order)
    if Xv.shape != X.shape or np.abs(Xv - U[-1] - X).max() > 1e-6:
        raise SystemExit(f"T6 mesh != rebuilt mesh: max |Xv - U - X| {np.abs(Xv - U[-1] - X).max():.3e}")

    # frame of U: a rotating-frame field has no 1P swing in the in-plane components at the tip
    tip = int(np.argmax(X[:, 2]))
    U_rot = U if U[:, tip, 0].std() < 0.1 else np.array([ct.to_rotating(u, a) for u, a in zip(U, th, strict=True)])
    u_frame = "rotating" if U_rot is U else "inertial->rotated"

    K, free = par.rotating_keff(asm, mesh, omega, args.hub_radius)
    mr, _, mv = asm._rust.assemble_m()
    nodal_m = np.zeros(mesh.node_count)
    for row, val in zip(np.asarray(mr), np.asarray(mv), strict=True):
        if row % dpn < 3:
            nodal_m[row // dpn] += val / 3.0
    calc = InertialForcesCalculator(rotation_axis=par.AXIS, rotation_center=center)
    f_cf = calc.compute_centrifugal_force(X, nodal_m, omega)  # the rotor solver's fallback path, X0

    # K_G the production Rust loop rebuilds every step (rotor_fsi.rs -> CentrifugalPrestress, T17/L33): the same
    # code is exposed by the PyO3 binding below, so static_t6_kgrust is the acceptance row (== static_t6, 1e-6 rel)
    from scipy.sparse import coo_matrix

    n_dof = asm.dofs_count
    r, c, v = asm._rust.assemble_k()
    K0 = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n_dof, n_dof)).tocsr()
    blk = -omega**2 * (np.eye(3) - np.outer(par.AXIS, par.AXIS))
    idx = np.arange(mesh.node_count) * dpn
    Ksp = coo_matrix((np.concatenate([nodal_m * blk[a, b] for a in range(3) for b in range(3)]),
                      (np.concatenate([idx + a for a in range(3) for b in range(3)]),
                       np.concatenate([idx + b for a in range(3) for b in range(3)]))), shape=(n_dof, n_dof)).tocsr()
    Kg_py = (K - K0 - Ksp).tocsr()
    r, c, v = asm._rust.assemble_centrifugal_geometric_k(
        float(omega), [float(x) for x in par.AXIS], [float(x) for x in center],
        np.asarray(free, dtype=np.int64))
    Kg_rs = coo_matrix((np.asarray(v), (np.asarray(r), np.asarray(c))), shape=(n_dof, n_dof)).tocsr()
    K_rust = (K0 + Ksp + Kg_rs).tocsr()

    def frob(A):
        return float(np.sqrt((A[free][:, free].power(2)).sum()))

    kg_norms = {"Kg_python_F": frob(Kg_py), "Kg_rust_F": frob(Kg_rs), "K_F": frob(K0)}
    print("K_G norms (free DOFs):", kg_norms, flush=True)

    def full(f3):
        f = np.zeros(asm.dofs_count)
        for d in range(3):
            f[d::dpn] = f3[:, d]
        return f

    def solve(f3, Kmat=None):
        u = solve_static(asm, K if Kmat is None else Kmat, free, full(f3))
        return np.column_stack([u[d::dpn] for d in range(3)])

    fa_mean, fi_mean, u_mean = tmean(t, FA_r), tmean(t, FI_r), tmean(t, U_rot)
    f_bem3 = np.column_stack([f_bem[d::dpn] for d in range(3)])
    sols = {"aero_t6": solve(fa_mean), "aero_t3": solve(f_bem3), "cf_calc": solve(f_cf), "inertial_t6": solve(fi_mean)}
    sols["static_t6"] = sols["aero_t6"] + sols["inertial_t6"]
    sols["static_t6_kgrust"] = solve(fa_mean + fi_mean, K_rust)
    sols["static_t6_nokg"] = solve(fa_mean + fi_mean, (K0 + Ksp).tocsr())
    sols["dyn_mean"] = u_mean

    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    ref = json.loads(args.reference.read_text())
    zq = np.array(ref["table"]["z"])
    table = {"z": zq.tolist()}
    for k, U3 in sols.items():
        sh = par.shell_sections(mesh, part, U3, kp_bd)
        table[k] = {"flap": np.interp(zq, sh["z"], sh["u"][:, 1]).tolist(),
                    "edge": np.interp(zq, sh["z"], sh["u"][:, 0]).tolist(),
                    "twist_eq_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_eq"])).tolist(),
                    "twist_skin_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_skin"])).tolist()}
    loads = {**kg_norms,"F_aero_t6_mean": fa_mean.sum(0).tolist(), "F_aero_t3": f_bem3.sum(0).tolist(),
             "F_inertial_t6_mean": fi_mean.sum(0).tolist(), "F_cf_calc": f_cf.sum(0).tolist(),
             "F_aero_rot_sum_std": FA_r.sum(1).std(0).tolist(),
             "tip_U_std_rot": U_rot[:, tip].std(0).tolist()}
    doc = {"results": str(args.results), "n_dumps": len(t), "window": [t[0], t[-1]], "n_rev": args.n_rev,
           "u_frame": u_frame, "tip_node": tip, "loads": loads, "table": table, "revision": tp.revision(REPO),
           "elapsed_s": time.time() - t_start}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(doc, indent=1))
    print(json.dumps({k: doc[k] for k in ("n_dumps", "window", "u_frame", "loads")}, indent=1))
    for k in sols:
        r = table[k]
        print(f"{k:12s} tip flap {r['flap'][-1]:8.4f} edge {r['edge'][-1]:7.4f} "
              f"tw_eq {r['twist_eq_deg'][-1]:7.3f} tw_skin {r['twist_skin_deg'][-1]:7.3f}")


if __name__ == "__main__":
    main()
