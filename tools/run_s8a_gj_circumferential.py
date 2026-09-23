"""S-8a — circumferential refinement of the torsional stiffness (Camarena 2025).

Camarena & Anderson (Sandia) traced the shell's -24% twist deficit (vs
solid) to the circumferential mesh resolution: fixing the number of
elements around the section contour per span station improved it.  This
tool reproduces that diagnostic on the AeroElast shell: the GJ is
measured at fixed spanwise element_size (0.5 m) while the airfoil
contour discretization (n_samples) is refined 120 -> 240 -> 480.  If
GJ_shell drops toward the BeamDyn 6x6 tables with refinement, the
documented S-1 +17% band is the same published shell artifact.

Method: S-7 construction — tip torque via transverse force couples,
section twist from the displacement field, station-by-station
comparison against theta(z) = M*int(dz/GJ_BeamDyn).

Output:
  docs/validation_data/generated/s8a_gj_circumferential.csv

Usage:
  python tools/run_s8a_gj_circumferential.py [--samples 120 240 480]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from run_s7_torsion import (  # noqa: E402
    BLADE_YAML,
    beam_twist_profile_deg,
    run_torsion_case,
    section_twists_deg,
)

AIRFOIL_DIR = REPO / "tests" / "airfoils"
OUT_CSV = REPO / "docs" / "validation_data" / "generated" / "s8a_gj_circumferential.csv"

MOMENT_NM = 1000.0
ELEMENT_SIZE = 0.5


def build_mesh_props(n_samples: int):
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(
        yaml_file=str(BLADE_YAML),
        airfoil_dir=str(AIRFOIL_DIR),
        element_size=ELEMENT_SIZE,
        n_samples=n_samples,
        airfoil_spacing="constant",
        span_grading="chord",
    )
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    return mesh, props


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", type=int, default=[120, 240, 480])
    args = ap.parse_args()

    print(f"[S-8a] GJ vs airfoil contour samples {args.samples} "
          f"(element_size={ELEMENT_SIZE}, M={MOMENT_NM:.0f} Nm)")
    results = []
    for ns in args.samples:
        print(f"[S-8a] n_samples={ns} ...", flush=True)
        mesh, props = build_mesh_props(ns)
        u = run_torsion_case(mesh, props, MOMENT_NM)
        u_arr = np.asarray(u, dtype=np.float64)
        zs, theta_s = section_twists_deg(mesh, u_arr)
        theta_b = beam_twist_profile_deg(MOMENT_NM, zs)
        mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
        ratio = float(np.mean(theta_s[mid] / theta_b[mid]))
        gj_ratio = 1.0 / ratio
        print(f"  nodes={mesh.coords_array.shape[0]:6d}  twist ratio shell/beam = "
              f"{ratio:.3f}  ->  GJ shell/beam = {gj_ratio:.3f}", flush=True)
        results.append({
            "n_samples": ns,
            "n_nodes": int(mesh.coords_array.shape[0]),
            "twist_ratio": ratio,
            "gj_ratio": gj_ratio,
        })

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w") as f:
        f.write("n_samples,n_nodes,twist_ratio_shell_beam,gj_ratio_shell_beam\n")
        for r in results:
            f.write(f"{r['n_samples']},{r['n_nodes']},{r['twist_ratio']:.5f},"
                    f"{r['gj_ratio']:.5f}\n")
    print(f"[S-8a] Wrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
