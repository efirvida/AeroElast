"""S-7 — torsional response validation: shell vs BeamDyn GJ tables.

Applies a torsional moment M_z at the blade tip via transverse force
couples on the tip-section nodes (NOT nodal moments on the θz drilling
DOF — for MITC shells the physical torsional response lives in the
section kinematics, and the drilling DOF is a numerical penalty of
0.15·E·t²). The tip section rotation is recovered from the in-plane
displacement field (best-fit rigid rotation, same construction as S-1)
and compared against the 1D beam integral

    θ(L) = M · ∫₀^L dz / GJ(z)

using the GJ(z) distribution from the official BeamDyn blade file
(``IEA-15-240-RWT_BeamDyn_blade.dat``, K[5,5] per station).

This is the global torsional response that feeds the BEM elastic-twist
feedback: S-1 showed the sectional shell GJ is ~17% above the BeamDyn
tables; S-7 checks whether the GLOBAL tip twist is consistent with that.

Output:
  docs/validation_data/generated/s7_torsion.csv

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/run_s7_torsion.py [--element-size 0.5] [--moment 1.0e6]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from tools.beam_reference import load_beamdyn_blade  # noqa: E402

REF_DIR = REPO_ROOT / "tests" / "IEA15MW" / "reference"
BEAMDYN_FILE = REF_DIR / "IEA-15-240-RWT_BeamDyn_blade.dat"
BLADE_YAML = REPO_ROOT / "tests" / "IEA-15-240-RWT.yaml"
OUT_DIR = REPO_ROOT / "docs" / "validation_data" / "generated"

OUT_DIR.mkdir(parents=True, exist_ok=True)


def build_mesh_and_properties(element_size: float):
    """Build the official WindIO blade mesh (same generator as S-5/S-2)."""
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    generator = BladeMesh(
        yaml_file=str(BLADE_YAML),
        airfoil_dir=str(REPO_ROOT / "tests" / "airfoils"),
        element_size=element_size,
        n_samples=300,
        airfoil_spacing="constant",
        span_grading="chord",
    )
    mesh = generator.generate(renumber="rcm")
    properties = build_rust_properties(generator.numad_mesh_data)
    return mesh, properties


def run_torsion_case(mesh, properties, torque: float) -> np.ndarray:
    """Clamp root, apply tip torque via force couples, return solution."""
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": properties,
        },
    }
    solver = StaticLinearSolver(mesh, cfg)

    dpn = solver.domain.dofs_per_node
    root_node_ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    node_id_to_idx = mesh.node_id_to_index
    root_dofs = sorted(
        node_id_to_idx[nid] * dpn + d for nid in root_node_ids for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    # Tip torque about the span via transverse force couples on the tip
    # section nodes (S-1 construction): F_y = torque · x_off / Σ(2·x_off²)
    coords = mesh.coords_array
    z = coords[:, 2]
    tip_nodes = np.nonzero(z >= z.max() - 1e-6)[0]
    tip_xyz = coords[tip_nodes]
    x_off = tip_xyz[:, 0] - tip_xyz[:, 0].mean()
    denom = max(2.0 * float(np.sum(x_off[x_off > 0] ** 2)), 1e-6)
    f_y = {int(t): [0.0, torque * xo / denom, 0.0] for t, xo in zip(tip_nodes, x_off)}

    dofs, vals = [], []
    for t, fvec in f_y.items():
        for d in range(3):
            dofs.append(int(t) * dpn + d)
            vals.append(float(fvec[d]))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])

    return solver.solve()


def section_twists_deg(mesh, u: np.ndarray, n_slices: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Section rotation about the span from the displacement field [deg].

    Same construction as S-1: best-fit rigid rotation
    θ = Σ(r_x·u_y − r_y·u_x) / Σ(r_x² + r_y²) over each section's nodes.
    """
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6
    coords = mesh.coords_array
    z = coords[:, 2]
    edges = np.linspace(z.min(), z.max(), n_slices + 1)
    zs, thetas = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        mask = (z >= a) & (z < b) if b < z.max() else (z >= a) & (z <= b)
        idx = np.nonzero(mask)[0]
        if len(idx) < 3:
            continue
        xyz = coords[idx]
        r = xyz - xyz.mean(axis=0)
        ux = np.asarray([u[int(i) * dpn + 0] for i in idx])
        uy = np.asarray([u[int(i) * dpn + 1] for i in idx])
        num = float((r[:, 0] * uy - r[:, 1] * ux).sum())
        den = float((r[:, 0] ** 2 + r[:, 1] ** 2).sum())
        if den < 1e-18:
            continue
        zs.append(float(xyz[:, 2].mean()))
        thetas.append(np.rad2deg(num / den))
    return np.asarray(zs), np.asarray(thetas)


def beam_twist_profile_deg(moment: float, z_query: np.ndarray, blade_length: float = 117.0) -> np.ndarray:
    """θ(z) = M·∫₀^z dz'/GJ(z') with the BeamDyn GJ distribution [deg].

    The BeamDyn table's last station (r/R=1.0, GJ ≈ 7e4 N·m²) is a
    degenerate razor-tip section — its contribution is excluded by never
    integrating past the last valid station (z=111.1).
    """
    bd = load_beamdyn_blade(BEAMDYN_FILE, blade_length)
    z = np.asarray(bd.z, dtype=np.float64)
    gj = np.asarray(bd.gj, dtype=np.float64)
    z_valid = z[z <= 111.2]
    gj_valid = gj[z <= 111.2]
    out = np.empty_like(z_query, dtype=np.float64)
    for i, zq in enumerate(z_query):
        m = z_valid <= min(zq, 111.2)
        theta_rad = moment * np.trapezoid(1.0 / gj_valid[m], z_valid[m])
        out[i] = np.rad2deg(theta_rad)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-size", type=float, default=0.5)
    ap.add_argument("--moment", type=float, default=1000.0, help="Tip torsional moment [N·m]")
    args = ap.parse_args()

    print(f"[S-7] Building blade mesh (element_size={args.element_size})...")
    mesh, properties = build_mesh_and_properties(args.element_size)
    print(f"[S-7] Mesh: {mesh.coords_array.shape[0]} nodes")

    print(f"[S-7] Solving torsion case (M_z = {args.moment:.1e} N·m via tip couples)...")
    u = run_torsion_case(mesh, properties, args.moment)
    u_arr = np.asarray(u, dtype=np.float64)

    zs, theta_s = section_twists_deg(mesh, u_arr)
    theta_beam = beam_twist_profile_deg(args.moment, zs)

    # Compare station by station; the degenerate BeamDyn tip station is
    # already excluded inside beam_twist_profile_deg.
    ratios = []
    print()
    print(f'{"z [m]":>8} {"shell [deg]":>12} {"beam [deg]":>12} {"ratio":>8}')
    for z, th_s, th_b in zip(zs, theta_s, theta_beam):
        if th_b < 1e-9:
            continue
        ratios.append(th_s / th_b)
        print(f'{z:8.1f} {th_s:12.4e} {th_b:12.4e} {th_s/th_b:8.3f}')
    ratios = np.asarray(ratios)

    # Mid-span stations (0.3-0.9 of span) — the reliable zone.
    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    ratio_mid = (theta_s[mid] / theta_beam[mid]).mean()
    print()
    print(f"  Mean shell/beam twist ratio (0.3-0.9 span): {ratio_mid:.3f}")
    print(f"  Implied shell GJ/beam GJ: {1.0/ratio_mid:.3f}")

    csv_path = OUT_DIR / "s7_torsion.csv"
    with csv_path.open("w") as f:
        f.write("z_m,theta_shell_deg,theta_beam_deg,ratio\n")
        for z, th_s, th_b in zip(zs, theta_s, theta_beam):
            f.write(f"{z:.6f},{th_s:.6e},{th_b:.6e},"
                    f"{th_s/th_b if th_b > 1e-9 else ''}\n")
    print(f"[S-7] Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
