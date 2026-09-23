"""S-x mesh convergence: S-4 rotating modal + S-5 parked static vs element size.

Runs the blade validations at multiple element sizes and reports the
convergence trend of:
  * S-4: the 1st flap / 1st edge natural frequencies (parked + rotating)
  * S-5: the parked static OoP tip deflection under the AD loads

Usage:
  python tools/run_sx_convergence.py --sizes 1.0 0.5 0.25 \
      --out <openfast .out> --ad-blade <AeroDyn15_blade.dat> \
      --csv docs/validation_data/generated/sx_mesh_convergence.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

RPM = 7.56


def _case_for(element_size, out_path, ad_blade_path):
    """One mesh generation shared by the S-4 and the S-5 measurements."""
    from run_s4_rotating_modal import assemble_rotating, classify_modes
    from run_s5_oneway import (
        apply_ad_loads,
        assemble_parked_k,
        build_assembly,
        extract_ad_loads,
        solve_static,
        tip_deflections,
    )

    asm, mesh = build_assembly(element_size)
    r_m, fn, ft, twist = extract_ad_loads(out_path, ad_blade_path)

    # S-4 modal
    f_park, modes_p, free_p = assemble_rotating(asm, mesh, 0.0, n_modes=6)
    f_rot, modes_r, free_r = assemble_rotating(asm, mesh, RPM * 2 * np.pi / 60.0, n_modes=6)
    lab_p = classify_modes(mesh, modes_p, free_p)
    lab_r = classify_modes(mesh, modes_r, free_r)
    f1f_park = min(f for f, l in zip(f_park, lab_p) if l == "flap")
    f1f_rot = min(f for f, l in zip(f_rot, lab_r) if l == "flap")
    f1e_park = min(f for f, l in zip(f_park, lab_p) if l == "edge")
    f1e_rot = min(f for f, l in zip(f_rot, lab_r) if l == "edge")

    # S-5 parked static
    K, free = assemble_parked_k(asm, mesh)
    f = apply_ad_loads(asm, mesh, r_m, fn, ft, twist)
    u = solve_static(asm, K, free, f)
    oop = abs(tip_deflections(mesh, u)[0][1])

    n_nodes = len(mesh.nodes)
    n_elems = len(asm.elements)
    return {
        "size": element_size,
        "n_nodes": n_nodes,
        "n_elems": n_elems,
        "f1f_park": f1f_park,
        "f1f_rot": f1f_rot,
        "f1e_park": f1e_park,
        "f1e_rot": f1e_rot,
        "oop_static": oop,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=float, nargs="+", default=[1.0, 0.5, 0.25])
    ap.add_argument("--out", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/"
                                 "IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.out"))
    ap.add_argument("--ad-blade", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/"
                                 "IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat"))
    ap.add_argument("--csv", type=Path,
                    default=Path("docs/validation_data/generated/sx_mesh_convergence.csv"))
    args = ap.parse_args()

    rows = []
    for size in args.sizes:
        print(f"\n=== element_size = {size} m ===", flush=True)
        r = _case_for(size, args.out, args.ad_blade)
        rows.append(r)
        print(f"  nodes={r['n_nodes']} elems={r['n_elems']} | "
              f"1F parked {r['f1f_park']:.4f} Hz, rotating {r['f1f_rot']:.4f} | "
              f"1E parked {r['f1e_park']:.4f} | OoP static {r['oop_static']:.4f} m",
              flush=True)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with open(args.csv, "w") as fh:
            fh.write("element_size_m,n_nodes,n_elems,f1f_park_Hz,f1f_rot_Hz,"
                     "f1e_park_Hz,f1e_rot_Hz,oop_static_m\n")
            for r in rows:
                fh.write(f"{r['size']},{r['n_nodes']},{r['n_elems']},"
                         f"{r['f1f_park']:.8f},{r['f1f_rot']:.8f},"
                         f"{r['f1e_park']:.8f},{r['f1e_rot']:.8f},"
                         f"{r['oop_static']:.8f}\n")
        print(f"\nwrote {args.csv}")

    if len(rows) >= 2:
        print("\nconvergence (relative change between consecutive sizes):")
        for key, label in [("f1f_park", "1F parked"), ("f1f_rot", "1F rotating"),
                           ("oop_static", "OoP static")]:
            deltas = [abs(rows[i][key] - rows[i + 1][key]) / rows[i + 1][key]
                      for i in range(len(rows) - 1)]
            print(f"  {label:12s}: " + " | ".join(f"{d:.2%}" for d in deltas))


if __name__ == "__main__":
    main()
