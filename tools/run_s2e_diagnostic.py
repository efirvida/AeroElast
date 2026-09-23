"""S-2E — diagnostic of the 2nd edgewise mode vs mesh refinement.

The 2nd edgewise frequency (shell 1.9796 Hz vs IEA-15-240-RWT 2.167 Hz,
-8.65%) is the only V-02 mode outside the published beam-shell band
(Almeida 2025: -0.8%..+3.7%; Faccio 2019: up to -5.7% on high modes).
This tool tests whether the deficit is a discretization artifact: it
extracts the first 10 modes at element sizes 1.0 / 0.5 / 0.25 m on the
official WindIO blade, identifies the edgewise modes by their
directional effective-mass participation (Tx dominant) and reports the
frequency trend.  If f2e converges from below toward ~2.1 Hz with
refinement, the -8.65% is partly mesh-induced (the S-shape 2E mode is
more curvature-sensitive than 1E).

Output:
  docs/validation_data/generated/s2e_mesh_convergence.csv

Usage (login node or queue, GCC runtime required):
  python tools/run_s2e_diagnostic.py [--sizes 1.0 0.5 0.25]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

BLADE_YAML = REPO / "tests" / "IEA-15-240-RWT.yaml"
AIRFOIL_DIR = REPO / "tests" / "airfoils"
OUT_CSV = REPO / "docs" / "validation_data" / "generated" / "s2e_mesh_convergence.csv"

REF_2E = 2.1670  # IEA-15-240-RWT package


def run_size(element_size: float) -> dict:
    from aeroelast.core.bc import DirichletCondition
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.elements import ElementFamily
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.modal import ModalSolver

    gen = BladeMesh(
        yaml_file=str(BLADE_YAML),
        airfoil_dir=str(AIRFOIL_DIR),
        element_size=element_size,
        n_samples=300,
        airfoil_spacing="constant",
        span_grading="chord",
    )
    mesh = gen.generate(renumber="rcm")
    properties = build_rust_properties(gen.numad_mesh_data)

    cfg = {
        "solver": {"num_modes": 10},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": properties,
        },
    }
    solver = ModalSolver(mesh, cfg)
    root = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    nid2idx = mesh.node_id_to_index
    dpn = mesh.dofs_per_node if hasattr(mesh, "dofs_per_node") else 6
    root_dofs = sorted(nid2idx[n] * dpn + d for n in root for d in range(dpn))
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    freqs, _ = solver.solve()
    freqs = np.asarray(freqs, dtype=float)
    factors = np.asarray(solver.participation_factors)  # (n_modes, 6) in %
    tx = factors[:, 0]

    order = np.argsort(freqs)
    edge = [int(k) for k in order if tx[k] > 0.3 * float(tx.max()) and tx[k] > 5.0]
    f1e = freqs[edge[0]] if edge else np.nan
    f2e = freqs[edge[1]] if len(edge) > 1 else np.nan

    print(f"  nodes={mesh.coords_array.shape[0]}  n_modes={len(freqs)}  "
          f"1E={f1e:.4f} Hz  2E={f2e:.4f} Hz  "
          f"(delta vs ref {100*(f2e-REF_2E)/REF_2E:+.1f}%)", flush=True)
    return {
        "element_size_m": element_size,
        "n_nodes": int(mesh.coords_array.shape[0]),
        "f1e_Hz": f1e,
        "f2e_Hz": f2e,
        "delta_2e_pct": 100 * (f2e - REF_2E) / REF_2E if np.isfinite(f2e) else np.nan,
        "tx_2e_pct": float(tx[edge[1]]) if len(edge) > 1 else np.nan,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", nargs="+", type=float, default=[1.0, 0.5, 0.25])
    args = ap.parse_args()

    print(f"[S-2E] Official blade, sizes {args.sizes}, ref 2E = {REF_2E} Hz")
    results = []
    for es in args.sizes:
        print(f"[S-2E] element_size={es} m ...", flush=True)
        results.append(run_size(es))

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w") as f:
        f.write("element_size_m,n_nodes,f1e_Hz,f2e_Hz,delta_2e_pct,tx_2e_pct\n")
        for r in results:
            f.write(f"{r['element_size_m']},{r['n_nodes']},{r['f1e_Hz']:.6f},"
                    f"{r['f2e_Hz']:.6f},{r['delta_2e_pct']:.3f},{r['tx_2e_pct']:.2f}\n")
    print(f"[S-2E] Wrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
