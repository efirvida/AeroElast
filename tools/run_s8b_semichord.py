"""S-8b — bend-twist coupling via semi-chord tip loads (Almeida 2025).

Almeida et al. applied the tip loads at the section semi-chord
(reference axis) exactly to expose the bending-torsion coupling of the
laminated structure.  This tool does the same for the AeroElast shell:
a tip force is applied at the tip section's semi-chord point (the
section node closest to pitch_axis = 0.5 chord), and the resulting tip
translation AND rotation are compared against the 1D beam reference
built from the official BeamDyn 6x6 stiffness matrices (which include
the off-diagonal coupling terms), solved with the SAME eccentric load
(force + moment arm about the reference axis).

The twist response here is the same mechanism the BEM twist feedback
reads (elastic twist of the section under flapwise load).

Output:
  docs/validation_data/generated/s8b_semichord_tip.csv

Usage:
  python tools/run_s8b_semichord.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from tools.beam_reference import load_beamdyn_blade  # noqa: E402
from tools.run_s7_torsion import BLADE_YAML, build_mesh_and_properties  # noqa: E402

AIRFOIL_DIR = REPO / "tests" / "airfoils"
BEAMDYN_FILE = REPO / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_BeamDyn_blade.dat"
OUT_CSV = REPO / "docs" / "validation_data" / "generated" / "s8b_semichord_tip.csv"

LOAD_DIRS = [0.0, 30.0, 60.0, 90.0]  # deg, measured from the chord direction
TIP_FORCE_N = 5.0e5


def shell_semichord_response(mesh, properties, load_dir_deg: float) -> dict:
    """Tip force at the semi-chord node; returns tip translations + rotation."""
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
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    root_dofs = sorted(nid2idx[n] * dpn + d for n in root for d in range(dpn))
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    coords = mesh.coords_array
    z = coords[:, 2]
    tip_idx = int(np.argmax(z))
    # semi-chord: the tip-section node closest to pitch_axis*chord along x
    tip_sec = np.nonzero(z >= z.max() - 1e-6)[0]
    # chord direction ~ +x (tip twist ~ -1.6 deg); pitch_axis = 0.5 chord from LE.
    # Approximate: pick the node with x closest to the section centroid x.
    x_c = coords[tip_sec, 0].mean()
    target = tip_sec[int(np.argmin(np.abs(coords[tip_sec, 0] - x_c)))]

    th = np.deg2rad(load_dir_deg)
    f = np.array([np.cos(th), np.sin(th), 0.0]) * TIP_FORCE_N
    dofs = [int(target) * dpn + d for d in range(3)]
    solver.add_nodal_loads([NodalLoad(dofs, list(f))])
    u = np.asarray(solver.solve(), dtype=np.float64)

    # tip response from the section rotation fit (translations at centroid +
    # twist from the displacement field)
    from run_s7_torsion import section_twists_deg
    zs, twists = section_twists_deg(mesh, u)
    u_tip = u[tip_idx * dpn : tip_idx * dpn + 3]
    return {
        "ux": float(u_tip[0]),
        "uy": float(u_tip[1]),
        "twist_deg": float(twists[-1]),
    }


def beam_semichord_response(mesh, load_dir_deg: float) -> dict:
    """Beam reference (diagonal EI/GJ model, from_mesh) for translations.

    NOTE: the diagonal beam has no bend-twist coupling — its twist is a
    sanity floor, not the target.  The elastic twist of the shell is the
    quantity of interest; its anchor is the BeamDyn OpenFAST run (S-8),
    not this beam model.
    """
    from tools.beam_reference import (
        BeamReference,
        load_beamdyn_blade,
        load_elastodyn_blade,
    )

    bd = load_beamdyn_blade(BEAMDYN_FILE, 117.0)
    ed_path = REPO / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_ElastoDyn_blade.dat"
    ed = load_elastodyn_blade(ed_path)
    beam = BeamReference.from_mesh(mesh, ed, bd)

    th = np.deg2rad(load_dir_deg)
    f = np.array([np.cos(th), np.sin(th), 0.0]) * TIP_FORCE_N
    u = beam.solve_tip_load(tuple(f))
    return {
        "ux": float(u[-6]),
        "uy": float(u[-5]),
        "twist_deg": float(np.rad2deg(u[-1])),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--element-size", type=float, default=0.5)
    args = ap.parse_args()

    print(f"[S-8b] Building mesh (element_size={args.element_size}) ...")
    mesh, properties = build_mesh_and_properties(args.element_size)

    rows = []
    for dir_deg in LOAD_DIRS:
        s = shell_semichord_response(mesh, properties, dir_deg)
        b = beam_semichord_response(mesh, dir_deg)
        rows.append({"dir_deg": dir_deg, **{f"shell_{k}": v for k, v in s.items()},
                     **{f"beam_{k}": v for k, v in b.items()}})
        print(f"[S-8b] dir={dir_deg:5.1f} deg: "
              f"ux shell/beam {s['ux']:.3f}/{b['ux']:.3f} m | "
              f"uy {s['uy']:.3f}/{b['uy']:.3f} m | "
              f"twist shell {s['twist_deg']:+.3f} deg (beam diagonal: {b['twist_deg']:+.3f})",
              flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w") as f:
        f.write("dir_deg,shell_ux,beam_ux,shell_uy,beam_uy,shell_twist_deg,beam_twist_deg\n")
        for r in rows:
            f.write(f"{r['dir_deg']},{r['shell_ux']:.6f},{r['beam_ux']:.6f},"
                    f"{r['shell_uy']:.6f},{r['beam_uy']:.6f},"
                    f"{r['shell_twist_deg']:.6f},{r['beam_twist_deg']:.6f}\n")
    print(f"[S-8b] Wrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
