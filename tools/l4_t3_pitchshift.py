"""L4 T3 diagnostic: equivalent pitch offset (linear polars) that zeroes ours-AD thrust / torque, per V-01 point;
plus tip-node load share on both sides."""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import l4_t5_options as t5  # noqa: E402
import l4_t3_span as sp  # noqa: E402

ob, p = t5.ob, t5.p
for case in sp.CASES:
    blade, ref = t5.ad_run(case, tilt=sp.TILT, precone=sp.CONE)
    s = p.BEMSolver(
        blade,
        rho=t5.RHO,
        mu=t5.MU_AD,
        precone=sp.CONE,
        tilt=sp.TILT,
        yaw=0.0,
        hub_height=ob.IEA15MW_HUB_HEIGHT,
        shear_exp=0.0,
    )
    t5._patch_linear(s, blade)
    ds = np.array([-0.3, -0.15, 0.0, 0.15, 0.3])
    TQ = np.array([[(r := s.compute(case[0], case[1], case[2] + d)).thrust, r.torque] for d in ds])
    eT, eQ = TQ[:, 0] / ref.thrust - 1, TQ[:, 1] / ref.torque - 1
    zT = -eT[2] / np.polyfit(ds, eT, 1)[0]
    zQ = -eQ[2] / np.polyfit(ds, eQ, 1)[0]
    print(
        f"{case}: linear-polar err T{eT[2]:+.2%} Q{eQ[2]:+.2%}; dT/dpitch {np.polyfit(ds, eT, 1)[0]:+.2%}/deg, dQ/dpitch {np.polyfit(ds, eQ, 1)[0]:+.2%}/deg; "
        f"pitch offset to zero T: {zT:+.3f} deg, Q: {zQ:+.3f} deg"
    )
    # tip node: share of ours thrust in last node (Np*local dr ~ last two nodes)
    l, _ = s.rotor.distributedAeroLoads(case[0], case[1], case[2], 0.0)
    r = blade.r
    print(
        f"   ours Np last 3 nodes {np.round(l['Np'][-3:], 1)} N/m at r={np.round(r[-3:], 2)}, Rtip={blade.rotor_radius}; AD Fn last 3 {np.round(ref.fn[-3:], 1)}"
    )
