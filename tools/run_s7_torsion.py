"""S-7 — torsional response validation: shell vs BeamDyn GJ tables.

Applies a torsional moment M_z at the blade tip as the Saint-Venant wall
shear flow of the tip ring (``section_torque_wall_forces``: multi-cell Bredt
on the walls' own G t; a ~x_off force couple excites section distortion that
the rigid fit reads as twist, blade-input-audit T39/T40) — NOT nodal
moments on the θz drilling DOF — for MITC shells the physical torsional response lives in the
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


def build_mesh_and_properties(element_size: float, skin_reference: str = "oml"):
    """Build the official WindIO blade mesh (same generator as S-5/S-2).

    ``skin_reference`` is passed to ``BladeMesh``: ``"oml"`` (the generator default, the
    CFD/FSI surface) or ``"mid"`` (nodes on the laminate mid-surface, as the S-1/S-2/V-02
    anchors build it).
    """
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    generator = BladeMesh(
        yaml_file=str(BLADE_YAML),
        airfoil_dir=str(REPO_ROOT / "tests" / "airfoils"),
        element_size=element_size,
        n_samples=300,
        airfoil_spacing="constant",
        span_grading="chord",
        skin_reference=skin_reference,
    )
    mesh = generator.generate(renumber="rcm")
    properties = build_rust_properties(generator.numad_mesh_data)
    return mesh, properties


def tip_couple(mesh, torque: float) -> tuple[np.ndarray, np.ndarray]:
    """Legacy ``F_y ~ x_off`` tip couple (no longer the S-7 load; kept for ``diagnose_oml_reference_surface``).

    Pure tip couple about the span axis, as transverse forces on the tip-ring nodes.

    Returns ``(tip_node_indices, F_y)`` with ``F_y = torque * x_off / denom``, where ``x_off`` is the
    chordwise node offset from the ring centroid.  By construction ``sum(x_off * F_y) = torque`` and
    ``sum(F_y) = 0``: a couple of the requested magnitude with no net force.

    The normaliser is the whole construction: it must be the ring's chordwise second moment,
    ``sum(x_off**2)``, because only then does ``sum(x_off * F_y) = torque`` hold for *any* node
    distribution.  ``2 * sum(x_off**2 for x_off > 0)`` equals it only for an x-symmetric ring, and the
    IEA 15 MW tip ring is not symmetric: that form under-delivered the applied moment by a
    mesh-dependent factor (0.8978 at ``element_size`` 0.25, 0.9154 on the coarser three) and carried the
    S-7 global twist ratio across 1.0, inverting the direction the S-7 module reports.  See
    ``odd/tasks/s7-torsion-ratio.md``.
    """
    coords = mesh.coords_array
    z = coords[:, 2]
    tip = np.nonzero(z >= z.max() - 1e-6)[0]
    x_off = coords[tip, 0] - coords[tip, 0].mean()
    denom = float(np.sum(x_off**2))
    if denom <= 0.0:
        raise ValueError("degenerate tip ring: zero chordwise second moment")
    return tip, torque * x_off / denom


def run_torsion_case(mesh, properties, torque: float) -> np.ndarray:
    """Clamp root, apply the tip torque as the end ring's wall shear flow, return solution."""
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
    root_dofs = sorted(node_id_to_idx[nid] * dpn + d for nid in root_node_ids for d in range(dpn))
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    # Tip torque about the span as the Saint-Venant wall shear flow of the end ring (multi-cell
    # Bredt on the walls' own G t, the BEM projector's machinery).  A force couple ~ x_off is not
    # that field: it excites section distortion the rigid-fit twist reads as rotation (T39).
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces

    coords = mesh.coords_array
    tip = np.nonzero(coords[:, 2] >= coords[:, 2].max() - 1e-6)[0]
    tip_nodes, f_tip = section_torque_wall_forces(
        mesh, tip, (0.0, 0.0, 1.0), properties, torque, coords=coords
    )
    dofs, vals = [], []
    for t, f in zip(tip_nodes, f_tip, strict=True):
        for d in range(3):
            dofs.append(int(t) * dpn + d)
            vals.append(float(f[d]))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])

    return solver.solve()


def section_twists_deg(
    mesh, u: np.ndarray, n_slices: int = 40, fit_nodes: str = "all"
) -> tuple[np.ndarray, np.ndarray]:
    """Section rotation about the span from the displacement field [deg].

    Same construction as S-1: best-fit rigid rotation
    θ = Σ(r_x·u_y − r_y·u_x) / Σ(r_x² + r_y²) over each section's nodes.

    ``fit_nodes="all"`` (default, unchanged) uses every node of the slice; ``"box"`` only the slice's box nodes
    (caps SPAR* and webs *_SW, ``sectional.box_node_mask``; l3-closure S15), so the in-plane distortion of the
    skin and TE panels is not read as twist. A slice with fewer than 6 box nodes gives NaN (kept, never dropped
    and never filled with all nodes): count them with ``np.isnan``.
    """
    from aeroelast.postprocess.sectional import BOX_MIN_NODES, _check_fit_nodes, box_node_mask

    _check_fit_nodes(fit_nodes)
    box = box_node_mask(mesh) if fit_nodes == "box" else None
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6
    coords = mesh.coords_array
    z = coords[:, 2]
    edges = np.linspace(z.min(), z.max(), n_slices + 1)
    zs, thetas = [], []
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        mask = (z >= a) & (z < b) if b < z.max() else (z >= a) & (z <= b)
        idx = np.nonzero(mask)[0]
        if len(idx) < 3:
            continue
        z_slice = float(coords[idx, 2].mean())  # the slice's own z, so the all and box fits share one grid
        if box is not None:
            idx = idx[box[idx]]
            if len(idx) < BOX_MIN_NODES:
                zs.append(z_slice)
                thetas.append(np.nan)
                continue
        xyz = coords[idx]
        r = xyz - xyz.mean(axis=0)
        ux = np.asarray([u[int(i) * dpn + 0] for i in idx])
        uy = np.asarray([u[int(i) * dpn + 1] for i in idx])
        num = float((r[:, 0] * uy - r[:, 1] * ux).sum())
        den = float((r[:, 0] ** 2 + r[:, 1] ** 2).sum())
        if den < 1e-18:
            continue
        zs.append(z_slice)
        thetas.append(np.rad2deg(num / den))
    return np.asarray(zs), np.asarray(thetas)


def beam_twist_profile_deg(
    moment: float, z_query: np.ndarray, blade_length: float = 117.0, beamdyn_file=None, z_cut: float = 111.2
) -> np.ndarray:
    """θ(z) = M·∫₀^z dz'/GJ(z') with the BeamDyn GJ distribution [deg].

    The integral runs to ``z`` itself: ``1/GJ`` is linear between stations (the trapezoid's own
    assumption), so the partial segment from the last station below ``z`` up to ``z`` is included.
    Until 2026-10-11 it was dropped, and that staircase under-stated the beam twist by 16 % on
    average over the 0.3-0.9 band (``odd/tasks/s7-torsion-ratio.md`` section E8).

    The BeamDyn table's last station (r/R=1.0, GJ ≈ 7e4 N·m²) is a
    degenerate razor-tip section — its contribution is excluded by never
    integrating past the last valid station (z=111.1).

    ``beamdyn_file`` (default: the IEA 15 MW deck) and ``z_cut`` (default 111.2 m = 0.9504 of the IEA 15 MW
    blade) let another blade use its own deck and razor-tip cut; the defaults are the historical reference.
    """
    bd = load_beamdyn_blade(BEAMDYN_FILE if beamdyn_file is None else beamdyn_file, blade_length)
    z = np.asarray(bd.z, dtype=np.float64)
    gj = np.asarray(bd.gj, dtype=np.float64)
    z_valid = z[z <= z_cut]
    inv_valid = 1.0 / gj[z <= z_cut]
    out = np.empty_like(z_query, dtype=np.float64)
    for i, zq in enumerate(z_query):
        z_end = min(float(zq), float(z_valid[-1]))
        m = z_valid <= z_end
        zz = z_valid[m]
        ii = inv_valid[m]
        if z_end > zz[-1]:
            zz = np.append(zz, z_end)
            ii = np.append(ii, np.interp(z_end, z_valid, inv_valid))
        theta_rad = moment * np.trapezoid(ii, zz)
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
    print(f"{'z [m]':>8} {'shell [deg]':>12} {'beam [deg]':>12} {'ratio':>8}")
    for z, th_s, th_b in zip(zs, theta_s, theta_beam, strict=True):
        if th_b < 1e-9:
            continue
        ratios.append(th_s / th_b)
        print(f"{z:8.1f} {th_s:12.4e} {th_b:12.4e} {th_s / th_b:8.3f}")
    ratios = np.asarray(ratios)

    # Mid-span stations (0.3-0.9 of span) — the reliable zone.
    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    ratio_mid = (theta_s[mid] / theta_beam[mid]).mean()
    print()
    print(f"  Mean shell/beam twist ratio (0.3-0.9 span): {ratio_mid:.3f}")
    print(f"  Implied shell GJ/beam GJ: {1.0 / ratio_mid:.3f}")

    csv_path = OUT_DIR / "s7_torsion.csv"
    with csv_path.open("w") as f:
        f.write("z_m,theta_shell_deg,theta_beam_deg,ratio\n")
        for z, th_s, th_b in zip(zs, theta_s, theta_beam, strict=True):
            f.write(f"{z:.6f},{th_s:.6e},{th_b:.6e},{th_s / th_b if th_b > 1e-9 else ''}\n")
    print(f"[S-7] Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
