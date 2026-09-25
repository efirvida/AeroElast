"""Blade structural matrix: the IEA 15 MW blade under every load case that
matters, AeroElast shell vs the beam reference, on the converged mesh.

Why a separate artifact instead of a plain test
-----------------------------------------------
The blade-level cases the S-2 already had (tip flap/edge, gravity, spanwise
distributed) are compared against a beam reference only.  This extends the
matrix with the cases the wind + rotation actually impose and that were
missing, and writes a CSV so the heavy side (full-blade CalculiX, ~250k DOF
with quadratic S8R on the 0.25 m mesh) can run separately on HPC and the test
stays fast -- the same artifact pattern the S-x suite uses.

Cases
-----
  B1 tip flap            (0, +F, 0)   at the tip
  B2 tip edge            (+F, 0, 0)   at the tip
  B3 surface distributed (0, +q, 0)   per unit span, over the whole blade
  B4 traction axial      (0, 0, +F)   at the tip
  B5 compression axial   (0, 0, -F)   at the tip
  B7 gravity             (0, -g, 0)   self weight

B6 (tip torsion) needs a couple, not a force; it stays in the S-7 torsion
test until this matrix grows a torque loader.
B8 (centrifugal) has no CalculiX counterpart -- the writer emits no
centrifugal *DLOAD -- so it is verified on the AeroElast side only (modal vs
OpenFAST, S-4).

Mesh: element_size defaults to 0.25 m, the converged size from the S-x study
(0.25 -> 0.125 moves f1f by 0.1% and the OoP deflection by 0.5%), which is
also the finest one full-blade CalculiX can realistically factor.

Usage:
  python tools/run_blade_structural_matrix.py --element-size 0.25 \
      --csv docs/validation_data/generated/blade_structural_matrix.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))
# The CalculiX helpers (_run_ccx, _frd_disp_at_point) live in the test tree.
sys.path.insert(0, str(REPO / "tests"))

FORCE_N = 1.0e6          # 1 MN tip load (matches the S-2 LC2/LC3 magnitude)
DISTRIB_Q = 5.0e3        # 5 kN/m spanwise distributed (matches the S-2 LC4)
GRAVITY = 9.81

CASES = (
    ("B1_tip_flap", {"type": "tip", "value": (0.0, FORCE_N, 0.0)}, "uy"),
    ("B2_tip_edge", {"type": "tip", "value": (FORCE_N, 0.0, 0.0)}, "ux"),
    ("B3_distributed", {"type": "distributed", "value": (0.0, DISTRIB_Q, 0.0),
                        "n_slices": 60}, "uy"),
    ("B4_traction_axial", {"type": "tip", "value": (0.0, 0.0, FORCE_N)}, "uz"),
    ("B5_compression_axial", {"type": "tip", "value": (0.0, 0.0, -FORCE_N)}, "uz"),
    ("B7_gravity", {"type": "gravity", "value": (0.0, -GRAVITY, 0.0)}, "uy"),
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-size", type=float, default=0.25)
    ap.add_argument("--blade", type=Path, default=REPO / "tests" / "IEA-15-240-RWT.yaml")
    ap.add_argument("--csv", type=Path,
                    default=REPO / "docs" / "validation_data" / "generated"
                    / "blade_structural_matrix.csv")
    ap.add_argument("--ccx", action="store_true",
                    help="also run full-blade CalculiX (heavy; intended for HPC)")
    ap.add_argument("--case-index", type=int, default=None,
                    help="run ONLY this case (0-based) and write --csv for it alone; "
                         "for a SLURM job array over the independent cases")
    ap.add_argument("--work-dir", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/tmp/opencode/blade_ccx"))
    args = ap.parse_args()

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from beam_reference import BeamReference
    from run_s2_static_cases import (
        beam_tip_displacement,
        run_shell_case,
        shell_tip_displacement,
    )

    print(f"building blade mesh (element_size={args.element_size} m) ...")
    gen = BladeMesh(element_size=args.element_size, yaml_file=str(args.blade))
    mesh = gen.generate(renumber="rcm", verbose=False)
    props = build_rust_properties(gen.numad_mesh_data)
    print(f"  mesh: {mesh.node_count} nodes, {mesh.elements_count} elements, "
          f"{len(props)} element sets")
    beam = BeamReference.from_mesh(mesh, *_load_refs())

    selected = list(enumerate(CASES))
    if args.case_index is not None:
        selected = [selected[args.case_index]]
        print(f"single-case mode: index {args.case_index} -> {selected[0][1][0]}")

    rows = []
    print(f"\n{'case':>22} {'shell [m]':>14} {'beam [m]':>14} {'rel err':>9}")
    for case_i, (name, load, comp) in selected:
        u = run_shell_case(mesh, props, load)
        d_shell = shell_tip_displacement(mesh, u)
        d_beam = beam_tip_displacement(beam, load)
        idx = {"ux": 0, "uy": 1, "uz": 2}[comp]
        s, b = float(d_shell[idx]), float(d_beam[idx])
        rel = abs(s - b) / max(abs(b), 1e-30)
        print(f"{name:>22} {s:14.6f} {b:14.6f} {rel*100:8.2f}%")
        rows.append({
            "case": name, "component": comp,
            "shell_m": s, "beam_m": b, "rel_err": rel,
            "shell_ux": d_shell[0], "shell_uy": d_shell[1], "shell_uz": d_shell[2],
            "beam_ux": d_beam[0], "beam_uy": d_beam[1], "beam_uz": d_beam[2],
            "element_size_m": args.element_size,
            "n_nodes": mesh.node_count,
        })

    if args.ccx:
        # Tip-load cases only; in single-case mode just this one.
        tip_cases = ("B1_tip_flap", "B2_tip_edge", "B4_traction_axial", "B5_compression_axial")
        if args.case_index is not None and rows[0]["case"] not in tip_cases:
            print(f"  {rows[0]['case']}: no CCX path (*DLOAD/couple not wired)")
        else:
            _run_ccx_column(mesh, props, rows, args)

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {args.csv}")


def _run_ccx_column(mesh, props, rows, args) -> None:
    """Full-blade CalculiX column for the tip-load cases.

    Only the tip-load cases are exported: the writer loads a named node set
    with a single vector, so the spanwise distributed (B3) and gravity (B7)
    cases have no *DLOAD path yet, and B6 needs a couple rather than a force.

    CalculiX needs quadratic S8R for *SHELL SECTION, COMPOSITE -- with linear
    S4 the writer falls back to a single-layer section and the laminate is
    silently lost.  Quadratic also means the FRD numbering includes mid-side
    nodes, so displacements are matched by COORDINATE, not by node id
    (_frd_disp_at_point, as test_orthotropic_shell_parity.py does).
    """
    import shutil

    from aeroelast.core.mesh.entities import NodeSet
    from aeroelast.core.mesh.io.writers import write_ccx_mesh
    from test_beam_4cases_parity import _run_ccx
    from test_orthotropic_shell_parity import _frd_disp_at_point

    ccx_bin = shutil.which("ccx") or shutil.which("CalculiX")
    if ccx_bin is None:
        raise SystemExit("CalculiX (ccx) not found in PATH")

    # Tip node set for the load (the blade mesh ships RootNodes only).
    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    z = coords[:, 2]
    tip_nodes = {n for n in mesh.nodes if n.z >= z.max() - 1e-6}
    if "TipNodes" not in mesh.node_sets:
        mesh.add_node_set(NodeSet("TipNodes", tip_nodes))
    tip_xyz = {mesh.node_id_to_index[n.id]: np.array([n.x, n.y, n.z]) for n in tip_nodes}

    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    print(f"\nCCX column (quadratic S8R, {len(tip_nodes)} tip nodes) -> {work}")

    for row in rows:
        name = row["case"]
        if name not in ("B1_tip_flap", "B2_tip_edge", "B4_traction_axial",
                        "B5_compression_axial"):
            print(f"  {name:>22}: skipped (needs *DLOAD / couple; not wired)")
            row["ccx_m"] = ""
            continue
        vec = {"B1_tip_flap": (0.0, FORCE_N, 0.0),
               "B2_tip_edge": (FORCE_N, 0.0, 0.0),
               "B4_traction_axial": (0.0, 0.0, FORCE_N),
               "B5_compression_axial": (0.0, 0.0, -FORCE_N)}[name]
        case_dir = work / name
        case_dir.mkdir(exist_ok=True)
        stem = f"blade_{name}"
        inp = case_dir / f"{stem}.inp"
        write_ccx_mesh(
            mesh, str(inp), properties=props,
            boundary_nodeset="RootNodes", solver_type="LinearStatic",
            load_nodeset="TipNodes", load_vector=list(vec),
            span_direction=(0.0, 0.0, 1.0), quadratic=True,
        )
        _run_ccx(ccx_bin, case_dir, stem)
        frd = case_dir / f"{stem}.frd"
        if not frd.exists():
            print(f"  {name:>22}: CCX produced no FRD")
            row["ccx_m"] = ""
            continue
        # average tip displacement over the tip section, same component as shell
        idx = {"ux": 0, "uy": 1, "uz": 2}[row["component"]]
        # The FRD holds the DEFORMED coordinates, so the tolerance must cover
        # the tip displacement (a 1 MN axial load moves it ~0.4 mm; the helper
        # defaults to 1e-9 and finds nothing).
        vals = [
            float(np.asarray(_frd_disp_at_point(frd, xyz, tol=5e-3))[idx])
            for xyz in tip_xyz.values()
        ]
        row["ccx_m"] = float(np.mean(vals))
        row["ccx_shell_rel_err"] = abs(row["shell_m"] - row["ccx_m"]) / max(abs(row["ccx_m"]), 1e-30)
        print(f"  {name:>22}: CCX={row['ccx_m']:14.6f} m  "
              f"(shell {row['shell_m']:.6f}, {row['ccx_shell_rel_err']*100:.2f}% apart)")


def _load_refs():
    """ElastoDyn + BeamDyn blade files used by the beam reference."""
    ref = REPO / "tests" / "IEA15MW" / "reference"
    from beam_reference import load_beamdyn_blade, load_elastodyn_blade

    ed = load_elastodyn_blade(ref / "IEA-15-240-RWT_ElastoDyn_blade.dat")
    bd = load_beamdyn_blade(ref / "IEA-15-240-RWT_BeamDyn_blade.dat")
    return ed, bd


if __name__ == "__main__":
    main()
