"""L3 T7 (W6): rated one-way G-C solve on gmsh chord-scaled N20 AR 1 with an outboard stiffness variant, then the BEM at the
deformed geometry (the two production feedback channels: strip radius, elastic twist). Same construction as
``tools/l3_torque_sensitivity.py`` (``solve_rated`` with the production wall flow, as ``diagnose_gc_twist_chord.py``).

    scripts/aeroenv.sh env PYTHONPATH=$PWD python tools/l3_w6_run.py {base,V1,V2,V3} logs/chord_scaled/w6/<variant>.json
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import l3_torque_sensitivity as l3  # noqa: E402
import l3_w6_scale as sc  # noqa: E402

FACTORS = Path(__file__).resolve().parents[1] / "logs/chord_scaled/w6/factors.json"


def main():
    variant, out = sys.argv[1], Path(sys.argv[2])
    do_g, do_e = variant in ("V1", "V3"), variant in ("V2", "V3")
    transform = sc.make_transform(FACTORS, do_g, do_e) if variant != "base" else None
    p, disp = l3.solve_rated(
        0.25, "gmsh", wall_flow=True, chord_divisions=20, mesh_transform=transform
    )
    r_ref, tw_ref = p._ref_r, p._ref_twist
    r_def, tw_def = p._compute_deformed_geometry(disp)
    mesh = p._mesh
    tip = np.array([mesh.node_id_to_index[i] for i in mesh.get_node_set("TipNodes").node_ids])
    cases = {"rigid": l3.evaluate(p, r_ref, tw_ref), "deformed": l3.evaluate(p, r_def, tw_def)}
    res = {
        "variant": variant,
        "cases": cases,
        "tip_x_edge_m": float(disp[tip, 0].mean()),
        "tip_y_flap_m": float(disp[tip, 1].mean()),
        "tip_z_m": float(disp[tip, 2].mean()),
        "tip_twist_deg": float(np.rad2deg((tw_def - tw_ref)[-1])),
        "min_dr_m": float((r_def - r_ref).min()),
        "max_dr_m": float((r_def - r_ref).max()),
        "twist_deg_strips": [float(v) for v in np.rad2deg(tw_def - tw_ref)],
        "r": [float(v) for v in r_ref],
        "n_tip_nodes": int(len(tip)),
    }
    out.write_text(json.dumps(res, indent=1))
    print(
        json.dumps({k: v for k, v in res.items() if k not in ("twist_deg_strips", "r")}, indent=1)
    )


if __name__ == "__main__":
    main()
