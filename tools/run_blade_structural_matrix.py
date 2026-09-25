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

    rows = []
    print(f"\n{'case':>22} {'shell [m]':>14} {'beam [m]':>14} {'rel err':>9}")
    for name, load, comp in CASES:
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
        print("\n--ccx requested: full-blade CalculiX is not wired here yet; "
              "see tests/test_composite_ccx_parity.py for the export path "
              "(properties=build_rust_properties(...), quadratic=True).")

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {args.csv}")


def _load_refs():
    """ElastoDyn + BeamDyn blade files used by the beam reference."""
    ref = REPO / "tests" / "IEA15MW" / "reference"
    from beam_reference import load_beamdyn_blade, load_elastodyn_blade

    ed = load_elastodyn_blade(ref / "IEA-15-240-RWT_ElastoDyn_blade.dat")
    bd = load_beamdyn_blade(ref / "IEA-15-240-RWT_BeamDyn_blade.dat")
    return ed, bd


if __name__ == "__main__":
    main()
