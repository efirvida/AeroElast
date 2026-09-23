"""S-8c — elastic twist sign verification vs the BeamDyn anchor.

Applies the S-5 AeroDyn distributed loads (Fn/Ft per unit length, mean
over the rated window) to the shell and measures the elastic twist
along the span.  The BeamDyn anchor (official IEA-15-240-RWT deck,
CompElast=2, rated run) gives NOSE-UP twist: +0.55 deg mid-span,
+0.98 deg tip.  Before the rotorspin fix the shell gave nose-down —
the inverted bend-twist coupling that loaded the FSI by +6.6%.

Output:
  docs/validation_data/generated/s8c_twist_sign.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from run_s7_torsion import build_mesh_and_properties, section_twists_deg  # noqa: E402

AD_LOADS = REPO / "tests" / "IEA15MW" / "reference" / "s5_ad_blade1_loads.csv"
OUT_CSV = REPO / "docs" / "validation_data" / "generated" / "s8c_twist_sign.csv"

# BeamDyn anchor (steady-state rated run, blade 1)
ANCHOR_MID = +0.545
ANCHOR_TIP = +0.981


def apply_ad_loads(mesh, properties, r_m, fn, ft, twist_deg):
    """Apply the S-5 AD distributed loads (local chord frame, WindIO sign)
    and solve static."""
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
    n = len(coords)

    fn_i = np.interp(z, r_m, fn)
    ft_i = np.interp(z, r_m, ft)
    tw_i = np.deg2rad(np.interp(z, r_m, twist_deg))

    # local chord frame: Fn along local y (flap), Ft along local x (edge)
    # mesh twist now matches the yaml sign; section frames rotate with the
    # geometric twist: [Fx, Fy] = R(tw) @ [Ft, Fn]
    c = np.cos(tw_i)
    s = np.sin(tw_i)
    fx = c * ft_i - s * fn_i
    fy = s * ft_i + c * fn_i

    # per-station spanwise tributary split over the station's nodes
    # (fixes the 1-2 nodes per section artefact; see load_mapping.py)
    from load_mapping import spanwise_tributary

    dz = spanwise_tributary(z)

    dofs = []
    vals = []
    for i in range(n):
        w = dz[i]
        dofs += [int(i) * dpn + 0, int(i) * dpn + 1]
        vals += [fx[i] * w, fy[i] * w]
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return solver.solve()


def main() -> int:
    print("[S-8c] Building mesh ...")
    mesh, properties = build_mesh_and_properties(0.5)

    ad = np.loadtxt(AD_LOADS, delimiter=",", skiprows=1)
    r_m, fn, ft, twist_deg = ad[:, 0], ad[:, 1], ad[:, 2], ad[:, 3]

    print("[S-8c] Solving with AD rated loads ...")
    u = np.asarray(apply_ad_loads(mesh, properties, r_m, fn, ft, twist_deg),
                   dtype=np.float64)

    zs, twists = section_twists_deg(mesh, u)
    mid = (zs >= 0.2 * 117.0) & (zs <= 0.8 * 117.0)
    mid_twist = float(np.mean(twists[mid]))
    tip_twist = float(twists[-1])

    print(f"[S-8c] elastic twist: mid-span {mid_twist:+.3f} deg, tip {tip_twist:+.3f} deg")
    print(f"[S-8c] BeamDyn anchor: mid-span {ANCHOR_MID:+.3f} deg, tip {ANCHOR_TIP:+.3f} deg")
    sign_ok = (np.sign(mid_twist) == np.sign(ANCHOR_MID)) and (
        np.sign(tip_twist) == np.sign(ANCHOR_TIP)
    )
    print(f"[S-8c] sign matches anchor: {'YES' if sign_ok else 'NO'}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w") as f:
        f.write("z_m,twist_deg\n")
        for zz, tw in zip(zs, twists):
            f.write(f"{zz:.4f},{tw:.6f}\n")
    print(f"[S-8c] Wrote {OUT_CSV}")
    return 0 if sign_ok else 1


if __name__ == "__main__":
    sys.exit(main())
