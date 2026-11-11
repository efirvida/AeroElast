"""T5b: per-station G-C elastic twist on the chord-scaled gmsh mesh (odd/tasks/gmsh-chord-scaled-mesh.md).

For one case (N<k> or es<h>) solves the rated one-way field (production wall flow) and writes, per BEM strip: radius,
the production ring twist (``BEMFSIParticipant._section_rotation``), the same fit with a per-ring bending-slope column
(``l3_torque_sensitivity._ring_rotation_with_slope``), the ring node count, its span spread (non-planarity) and the
strip-mean flap/edge displacement. Diagnostic only.

    scripts/aeroenv.sh python tools/diagnose_gc_twist_chord.py N20 --out logs/chord_scaled/t5b
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main() -> int:
    import l3_torque_sensitivity as l3

    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--out", default="logs/chord_scaled/t5b")
    a = ap.parse_args()
    if a.case.startswith("N"):
        es, nd = 0.25, float(a.case[1:])
    else:
        es, nd = float(a.case[2:]), None
    p, disp = l3.solve_rated(es, "gmsh", wall_flow=True, chord_divisions=nd)
    span = p._span_dir
    sense = p._twist_mesh_to_bem
    rows = []
    for k, ring in enumerate(p._strip_ring_indices):
        idx = p._strip_node_indices[k]
        if len(idx) == 0:
            continue
        w0 = p._section_rotation(disp, ring)
        w1 = l3._ring_rotation_with_slope(p, disp, ring)
        zs = p._ref_coords[ring] @ span
        um = disp[idx].mean(axis=0)
        rows.append({
            "k": k, "r": float(p._ref_r[k]), "n_ring": int(len(ring)), "ring_z_ptp": float(np.ptp(zs)),
            "twist_prod_deg": None if w0 is None else float(np.rad2deg(sense * w0)),
            "twist_slope_deg": None if w1 is None else float(np.rad2deg(sense * w1)),
            "u_mean": [float(v) for v in um],
        })
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{a.case}.json").write_text(json.dumps({"case": a.case, "n_nodes": int(disp.shape[0]), "strips": rows}, indent=1))
    for r in rows[-8:]:
        print(r["k"], round(r["r"], 2), r["n_ring"], round(r["ring_z_ptp"], 4), r["twist_prod_deg"], r["twist_slope_deg"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
