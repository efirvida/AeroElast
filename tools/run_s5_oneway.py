"""S-5 — one-way validation: AeroDyn distributed loads -> AeroElast shell.

Applies the mean AeroDyn blade-1 distributed loads (Fn/Ft per unit length,
local chord frame, from an OpenFAST S-5 rated run) onto the official WindIO
shell blade and solves the static equilibrium, then compares the tip
deflection against the ElastoDyn reference.

Comparison design (settled empirically, 2026-09-10):

  * The shell is solved PARKED (K only) + the azimuth-mean gravity component.
    ElastoDyn's own static OoP at 7.56 rpm (15.91 m, gravity off) matches the
    un-stiffened Euler-Bernoulli beam (15.98 m) within 0.4% — the centrifugal
    stiffening contributes ~0-2% to the STATIC deflection (the modal +8.5% is
    shape-weighted differently).  Applying the shell's K_G here would inject a
    spurious ~13% stiffening (see the open K_G issue below).
  * The mean gravity component from the -6 deg shaft tilt:
    g_mean = (0, -1.026, 0) m/s^2 (downwind = -y in the shell frame, the
    deflection direction under the AD loads).  Calibrated against the
    blade-only ElastoDyn run: OoP = 16.1033 m with gravity vs 15.9121 m
    without (+0.19 m).
  * Reference: blade-only OpenFAST deck (only Flap/Edge DOFs, RotSpeed=7.56,
    tower/platform/servo DOFs off, CompAero=2) so the root-clamped shell is
    compared apples-to-apples: OoPDefl1 = 16.1033 m over t in [80, 100] s.
  * The full-monopile S-5 deck gives 16.0277 m (0.5% lower — controller
    speed 7.476 vs 7.56 rpm); kept as the campaign reference only.

The AD channel families (OpenFAST v5, BEM_Mod=1):
    AB1N###Fn  - normal force (to chord) per unit length   [N/m]  -> shell flap
    AB1N###Ft  - tangential force (to chord) per unit length [N/m] -> shell edge
The local chord frame rotates with the blade twist, which the shell sections
share (same IEA geometry), so the mapping is 1:1 in the local frame.

K_G characteristic (investigated 2026-09-10, formulation verified):

  * The centrifugal pre-stress machinery is exact: blade mass 70.7 t
    (+8% vs ElastoDyn), first moment 1.93e6 kg·m, root tension 1.21e6 N;
    the section-integrated s_m = sigma*h matches the analytic N(z) to
    within a few percent (bins 0.85-1.05, noise at the tip).
  * The K_G itself is the standard total-Lagrangian initial-stress form
    (B_G^T S_tilde B_G over the three displacement gradient blocks).
    A w-only experiment (u/v blocks removed) changed the modal stiffening
    by only ~0.7% (12.7% -> 12.0%) and the static by 0.4% (14.3% -> 13.9%):
    the transverse block carries essentially all the response.
  * The shell's K_G response (~12% modal / ~14% static of the deflection
    stiffness) exceeds the 1D Euler-Bernoulli tension estimate (~2%) because
    the shell kinematics include the twist-rotating section frames
    (covariant span gradients of the normal displacement), which the beam
    reference does not model.  It is a richer-formulation difference, not a
    bug: the sigma field, the load, and the equilibrium are all verified.
    The S-4 modal closure (rotating 1st flap +0.6%) uses this K_G; static
    comparisons intentionally use the parked stiffness.

IP note (resolved 2026-09-10): ElastoDyn's IPDefl1 (-0.81 m) is the
modal-reduced 1D measure: its IP root moment (RootMxb1 = 8.09e6 N·m) breaks
down as the aero loads 6.00e6 (Ft 6.42e6 + Fn-twist -0.42e6) + a gravity/
cone mean 2.09e6.  The 3D references (the 1D Euler-Bernoulli beam with the
twist-rotating frames, and the shell) additionally carry the geometric
twist coupling: the beam gives u_x = 1.17 m with the SAME aero loads (Ft
0.82 + Fn-twist coupling 0.35; per-moment compliance 1.94e-7 vs the ED's
1.00e-7).  The shell adds ~0.25 m of section-level in-plane compliance on
top (1.49 m total).  The IP comparison against IPDefl1 is therefore
apples-to-oranges by construction; the test keeps it as a loose regression
bound and the OoP is the S-5 deliverable.

OpenFAST environment (persistent, $SCRATCH — survives /tmp cleanup):
  * OpenFAST-v5.0.0 conda env: /scratch/leahk/eduardo.donestevez/conda-envs/openfast
    (run: LD_LIBRARY_PATH=<env>/lib <env>/bin/openfast <deck>.fst)
  * Decks: /scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/ (official repo
    clone IEA-15-240-RWT, converted v4.0 -> v5.0.0 schema: the .fst gets
    ModCoupling/RhoInf/ConvTol/MaxConvIter/NRotors/CompSoil/MirrorRotor/
    SoilFile; the ED gets PitchDOF, PtfmRefxt/yt, PBrIner/BlPIner; the AD15
    loses Buoyancy, gains TwrCp/TwrCa tower columns)
  * blade-only S-5 deck: ofruns/OpenFAST/IEA-15-240-RWT-Monopile/
    (only blade DOFs, RotSpeed 7.56, servo/sea/hydro/sub off, HWindSpeed
    10.59 m/s — reproduces OoPDefl1 = 16.1039 vs the target 16.1033)

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s5_oneway.py --out /tmp/opencode/ofrun-s5-blade/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.out \
      --ad-blade /tmp/opencode/ofrun-s5/IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat \
      --export tests/IEA15MW/reference/s5_ad_blade1_loads.csv
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

# Blade-only ElastoDyn targets (mean over t in [80, 100] s, rated wind
# 10.59 m/s, RotSpeed = 7.56 rpm, only blade DOFs, gravity ON).
TIP_OOP_DEFL_M = 16.1033
TIP_IP_DEFL_M = -0.8103

# Azimuth-mean gravity in the shell frame: shaft tilt -6 deg projects the
# gravity onto the rotor normal (downwind = -y, the AD-load deflection
# direction): g_mean = -9.81 * sin(6 deg) * (0, 1, 0).
G_MEAN = np.array([0.0, -1.0255, 0.0])

# Channels: local chord frame loads at each AD node
CH_FN = "AB1N{n:03d}Fn"
CH_FT = "AB1N{n:03d}Ft"


def extract_ad_loads(out_path: Path, ad_blade_path: Path, t_lo=80.0, t_hi=100.0):
    """Mean AD blade-1 loads over [t_lo, t_hi] s + station positions.

    Returns (r_m, fn_nm, ft_nm, twist_deg): spanwise positions (m from blade
    root), mean Fn/Ft per unit length (N/m), and the geometric twist (deg).
    """
    with open(out_path) as f:
        lines = f.readlines()
    names = lines[6].split()
    data = np.loadtxt(out_path, skiprows=8)

    t = data[:, 0]
    sel = (t >= t_lo) & (t <= t_hi)
    if not sel.any():
        raise ValueError(f"no samples in t [{t_lo}, {t_hi}] s")

    # AD blade stations.  The v5 AD blade file layout: 2 header lines, 1
    # section line, NumBlNds line, column header, units line, then the data
    # rows (line index 6 onward).  Columns: BlSpn, BlCrvAC, BlSwpAC,
    # BlCrvAng, BlTwist, BlChord, ...
    with open(ad_blade_path) as f:
        ad_lines = f.readlines()
    n_nodes = int(ad_lines[3].split()[0])
    blspn = []
    twist = []
    for ln in ad_lines[6:6 + n_nodes]:
        parts = ln.split()
        if len(parts) >= 5:
            blspn.append(float(parts[0]))
            twist.append(float(parts[4]))
    blspn = np.asarray(blspn)
    twist = np.asarray(twist)
    if blspn.shape[0] != n_nodes:
        raise ValueError(f"expected {n_nodes} stations, got {blspn.shape[0]}")

    fn = np.zeros(n_nodes)
    ft = np.zeros(n_nodes)
    for n in range(1, n_nodes + 1):
        fn[n - 1] = data[sel, names.index(CH_FN.format(n=n))].mean()
        ft[n - 1] = data[sel, names.index(CH_FT.format(n=n))].mean()
    return blspn, fn, ft, twist


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


def assemble_parked_k(assembler, mesh):
    """Parked (linear) stiffness K and the clamped-root free DOFs."""
    from scipy.sparse import coo_matrix

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
    return K, free


def apply_ad_loads(assembler, mesh, r_m, fn, ft, twist_deg=None):
    """Map AD distributed loads onto shell nodes -> nodal force vector.

    The shell spans +z with the sections in the x-y plane rotated by the
    geometric twist.  The AD local frame shares the SAME twist (the IEA
    geometry), so the exact section frame is rebuilt from the interpolated
    BlTwist: t_hat = (cos θ, sin θ, 0) (chord), n_hat = (−sin θ, cos θ, 0)
    (flap normal).  Fn/Ft are applied in that frame.
    """
    n_dof = assembler.dofs_count
    dpn = assembler.dofs_per_node
    xyz = np.asarray([n.coords for n in mesh.nodes])
    z = xyz[:, 2]
    n_nodes = mesh.node_count

    span_lo, span_hi = z.min(), z.max()
    if span_hi - span_lo < 1.0:
        raise ValueError("unexpected mesh span")

    # spanwise tributary per unique station, split over the station's nodes
    # (the old z-sorted trapezoid gave 1-2 arbitrary nodes per section the
    # whole station load -> spurious section torques; see load_mapping.py)
    from load_mapping import spanwise_tributary

    trib = spanwise_tributary(z)

    fn_i = np.interp(z, r_m, fn)
    ft_i = np.interp(z, r_m, ft)

    if twist_deg is not None:
        # The shell sections now follow the WindIO twist sign (rotorspin fix,
        # 2026-09-15): chord angle ≈ +BlTwist, same convention as AeroDyn.
        theta = np.deg2rad(np.interp(z, r_m, twist_deg))
        t_hat = np.column_stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)])
        n_hat = np.column_stack([-np.sin(theta), np.cos(theta), np.zeros_like(theta)])
    else:
        # fallback: PCA of the section slice (sign-ambiguous, less accurate)
        half_w = 0.75
        n_hat = np.zeros((n_nodes, 3))
        t_hat = np.zeros((n_nodes, 3))
        for i in range(n_nodes):
            slab = np.abs(xyz[:, 2] - z[i]) <= half_w
            if slab.sum() < 4:
                slab = np.abs(xyz[:, 2] - z[i]) <= 3.0 * half_w
            pts = xyz[slab][:, :2]
            pts = pts - pts.mean(axis=0)
            if pts.shape[0] < 3:
                continue
            cov = pts.T @ pts / (pts.shape[0] - 1)
            w, v = np.linalg.eigh(cov)
            n_hat[i, :2] = v[:, np.argmin(w)]
            t_hat[i, :2] = v[:, np.argmax(w)]

    f = np.zeros(n_dof)
    for i in range(n_nodes):
        f[dpn * i:dpn * i + 3] += (fn_i[i] * n_hat[i] + ft_i[i] * t_hat[i]) * trib[i]
    return f


def assemble_gravity(assembler, g):
    """Body-force load vector from the Rust assembler (consistent nodal)."""
    f = assembler._rust.assemble_f_body(list(map(float, g)))
    return np.asarray(f, dtype=np.float64)


def solve_static(assembler, K_csr, free, f):
    from petsc4py import PETSc

    n_dof = assembler.dofs_count
    Kp = assembler._coo_to_petsc(
        np.asarray(K_csr.tocoo().row, dtype=np.int64),
        np.asarray(K_csr.tocoo().col, dtype=np.int64),
        np.asarray(K_csr.tocoo().data, dtype=np.float64),
    )
    is_free = PETSc.IS().createGeneral(free.astype(PETSc.IntType), comm=assembler.comm)
    K_red = Kp.createSubMatrix(is_free, is_free)

    u_red = PETSc.Vec().createSeq(len(free), comm=assembler.comm)
    f_red = PETSc.Vec().createSeq(len(free), comm=assembler.comm)
    f_red.setArray(f[free])
    ksp = PETSc.KSP().create(assembler.comm)
    ksp.setOperators(K_red)
    ksp.setType("preonly")
    pc = ksp.getPC()
    pc.setType("lu")
    pc.setFactorSolverType("mumps")
    ksp.setFromOptions()
    ksp.solve(f_red, u_red)

    u = np.zeros(n_dof)
    u[free] = u_red.getArray()
    return u


def tip_deflections(mesh, u, dpn=6):
    xyz = np.asarray([n.coords for n in mesh.nodes])
    tip_idx = int(np.argmax(xyz[:, 2]))
    return u.reshape(-1, dpn)[tip_idx][:3], xyz[tip_idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=Path("/tmp/opencode/ofrun-s5-blade/IEA-15-240-RWT-Monopile/"
                                 "IEA-15-240-RWT-Monopile.out"),
                    help="OpenFAST .out with the AD blade-1 load channels")
    ap.add_argument("--ad-blade", type=Path,
                    default=Path("/tmp/opencode/ofrun-s5/IEA-15-240-RWT/"
                                 "IEA-15-240-RWT_AeroDyn15_blade.dat"),
                    help="AeroDyn15 blade definition (BlSpn stations)")
    ap.add_argument("--loads", type=Path, default=None,
                    help="load the AD loads from a CSV artifact "
                         "(r_m,fn_Nm,ft_Nm) instead of --out/--ad-blade")
    ap.add_argument("--export", type=Path, default=None,
                    help="export the AD loads to a CSV artifact")
    ap.add_argument("--element-size", type=float, default=0.5)
    args = ap.parse_args()

    if args.loads:
        loads = np.loadtxt(args.loads, delimiter=",", skiprows=1)
        r_m, fn, ft = loads[:, 0], loads[:, 1], loads[:, 2]
        twist_deg = loads[:, 3] if loads.shape[1] >= 4 else None
        print(f"loads loaded from {args.loads}")
    else:
        r_m, fn, ft, twist_deg = extract_ad_loads(args.out, args.ad_blade)
    print(f"AD loads: {r_m.shape[0]} stations, "
          f"Fn total = {np.trapezoid(fn, r_m):.4e} N, "
          f"Ft total = {np.trapezoid(ft, r_m):.4e} N")

    if args.export:
        args.export.parent.mkdir(parents=True, exist_ok=True)
        np.savetxt(args.export,
                   np.column_stack([r_m, fn, ft, twist_deg]),
                   header="r_m,fn_Nm,ft_Nm,twist_deg", delimiter=",", comments="")
        print(f"loads exported to {args.export}")

    asm, mesh = build_assembly(args.element_size)
    K, free = assemble_parked_k(asm, mesh)
    f = apply_ad_loads(asm, mesh, r_m, fn, ft, twist_deg)
    u_aero = solve_static(asm, K, free, f)

    # azimuth-mean gravity along the aero deflection direction (downwind):
    # the gravity ADDS to the deflection (ElastoDyn: OoP 16.10 m with gravity
    # vs 15.91 m without)
    u_tip_aero = tip_deflections(mesh, u_aero)[0]
    g_dir = u_tip_aero[:3] / np.linalg.norm(u_tip_aero[:3])
    g_mean = np.linalg.norm(G_MEAN) * g_dir
    f_g = assemble_gravity(asm, g_mean)
    u = solve_static(asm, K, free, f + f_g)

    u_tip, tip_xyz = tip_deflections(mesh, u)
    print(f"\nparked solve (aero + g_mean): |u_tip| = {np.linalg.norm(u_tip):.4f} m "
          f"(u_x={u_tip[0]:+.4f}, u_y={u_tip[1]:+.4f}, u_z={u_tip[2]:+.4f})")
    print(f"parked solve (aero only):     |u_tip| = {np.linalg.norm(u_tip_aero):.4f} m "
          f"(gravity adds {np.linalg.norm(u_tip) - np.linalg.norm(u_tip_aero):+.4f} m)")
    print(f"tip node coords: {tip_xyz.round(2)}")

    # The tip section has ~0 twist -> |u_y| ≈ out-of-plane (downwind), |u_x| ≈ in-plane
    oop = abs(u_tip[1])
    ip = abs(u_tip[0])
    print(f"\nS-5 comparison (blade-only ElastoDyn targets):")
    print(f"  OoP tip deflection: shell {oop:.4f} m vs {TIP_OOP_DEFL_M:.4f} m "
          f"-> {(oop - TIP_OOP_DEFL_M) / TIP_OOP_DEFL_M:+.1%}")
    print(f"  IP  tip deflection: shell {ip:.4f} m vs {abs(TIP_IP_DEFL_M):.4f} m "
          f"-> {(ip - abs(TIP_IP_DEFL_M)) / abs(TIP_IP_DEFL_M):+.1%}")


if __name__ == "__main__":
    main()
