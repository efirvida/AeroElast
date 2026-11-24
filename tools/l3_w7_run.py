"""L3 T10 (W7): rated one-way G-C solve on gmsh chord-scaled N20 AR 1 with an optional end-region GJ variant (same as tools/l3_w6_run.py).

    scripts/aeroenv.sh env PYTHONPATH=$PWD python tools/l3_w7_run.py {base|tipstiff:F:Z0:Z1} OUT.json
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import l3_torque_sensitivity as l3  # noqa: E402
import l3_w7_tipstiff as ts  # noqa: E402


def main():
    variant, out = sys.argv[1], Path(sys.argv[2])
    transform = None if variant == "base" else ts.from_spec(variant)
    p, disp = l3.solve_rated(
        0.25, "gmsh", wall_flow=True, chord_divisions=20, mesh_transform=transform
    )
    r_ref, tw_ref = p._ref_r, p._ref_twist
    r_def, tw_def = p._compute_deformed_geometry(disp)
    mesh = p._mesh
    tip = np.array([mesh.node_id_to_index[i] for i in mesh.get_node_set("TipNodes").node_ids])
    res = {
        "variant": variant,
        "cases": {"rigid": l3.evaluate(p, r_ref, tw_ref), "deformed": l3.evaluate(p, r_def, tw_def)},
        "tip_x_edge_m": float(disp[tip, 0].mean()),
        "tip_y_flap_m": float(disp[tip, 1].mean()),
        "tip_z_m": float(disp[tip, 2].mean()),
        "tip_twist_deg": float(np.rad2deg((tw_def - tw_ref)[-1])),
        "min_dr_m": float((r_def - r_ref).min()),
        "twist_deg_strips": [float(v) for v in np.rad2deg(tw_def - tw_ref)],
        "r": [float(v) for v in r_ref],
    }
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k not in ("twist_deg_strips", "r")}, indent=1))


if __name__ == "__main__":
    main()
