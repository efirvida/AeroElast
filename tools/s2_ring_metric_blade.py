#!/usr/bin/env python
"""S-2 tip metrics on an official blade: node, tip-ring mean, Saint-Venant tip (l3-closure T10b).

Runs the four S-2 cases (``tools/run_s2_static_cases.py`` loads, clamped root, mid-surface skin) on a windIO yaml and
writes, per case, the three tip estimators of the shell translation:
  node   ``argmax z``, the first tip-ring node in mesh order (the metric the S-2 test used)
  ring   mean over the tip-ring nodes (L12(4))
  sv     ``saint_venant_tip``: ring means 1-3 tip chords inboard, cubic in z, at the tip (T10a: within 0.43 %
         of the closed form on prismatic tubes, where the tip-ring mean missed by 1.24 %)
plus the ring node counts used, so the extrapolation's stations can be audited. The BeamDyn reference is not read
here; compare with ``tools/bd_static_reference.py post`` (``--blade iea34`` for the 3.4).

Usage:  scripts/aeroenv.sh python tools/s2_ring_metric_blade.py --yaml tests/IEA-15-240-RWT.yaml \\
            --element-size 0.25 --out logs/check/l3/t10b/iea15-es025.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from run_s2_static_cases import (  # noqa: E402
    DISTRIB_Q_NM,
    GRAVITY,
    TIP_FORCE_N,
    ring_mean_at,
    run_shell_case,
    saint_venant_tip,
    shell_tip_displacement,
    tip_ring_scale,
)

CASES = [
    ("LC1_gravity", {"type": "gravity", "value": GRAVITY}),
    ("LC2_tip_flap", {"type": "tip", "value": (0.0, TIP_FORCE_N, 0.0)}),
    ("LC3_tip_edge", {"type": "tip", "value": (TIP_FORCE_N, 0.0, 0.0)}),
    ("LC4_distributed", {"type": "distributed", "value": (0.0, DISTRIB_Q_NM, 0.0), "n_slices": 60}),
]
DEPTHS = (1.0, 1.5, 2.0, 2.5, 3.0)
DEPTHS_FAR = (3.0, 4.0, 5.0, 6.0)  # stability check: a longer distortional decay would move sv


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--element-size", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.element_size >= 0.3:
        raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties

    gen = BladeMesh(yaml_file=args.yaml, element_size=args.element_size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    z = mesh.coords_array[:, 2]
    scale = tip_ring_scale(mesh)
    z_tip = float(z.max())
    stations = [ring_mean_at(mesh, np.zeros(6 * len(z)), z_tip - d * scale) for d in DEPTHS]
    out = {
        "yaml": args.yaml,
        "element_size": args.element_size,
        "n_nodes": int(len(z)),
        "tip_z": z_tip,
        "tip_ring_nodes": int(np.sum(z >= z_tip - 1e-6)),
        "tip_scale_m": scale,
        "sv_stations": [{"z": s[0], "nodes": s[2]} for s in stations],
        "cases": {},
    }
    print(
        f"S-2 metrics {Path(args.yaml).name} es {args.element_size}: {out['n_nodes']} nodes, tip z {z_tip:.3f}, "
        f"tip ring {out['tip_ring_nodes']} nodes, scale {scale:.3f} m, stations "
        + ", ".join(f"{s[0]:.2f} ({s[2]})" for s in stations),
        flush=True,
    )

    for name, load in CASES:
        u = np.asarray(run_shell_case(mesh, props, load))
        node = shell_tip_displacement(mesh, u)
        _z, ring, _n = ring_mean_at(mesh, u, z_tip)
        sv = saint_venant_tip(mesh, u, scale, DEPTHS)
        sv_far = saint_venant_tip(mesh, u, scale, DEPTHS_FAR)
        out["cases"][name] = {
            "node": node.tolist(),
            "ring": ring.tolist(),
            "sv": sv.tolist(),
            "sv_far": sv_far.tolist(),
        }
        k = 0 if name == "LC3_tip_edge" else 1
        print(
            f"  {name:16s} {'xy'[k]}: node {node[k]:.6f}  ring {ring[k]:.6f}  sv {sv[k]:.6f}  "
            f"sv_far {sv_far[k]:.6f}  (ring/sv {ring[k] / sv[k] - 1:+.3%}, node/sv {node[k] / sv[k] - 1:+.3%}, "
            f"sv_far/sv {sv_far[k] / sv[k] - 1:+.3%})",
            flush=True,
        )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
