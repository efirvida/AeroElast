"""T18 stage 2: static aero-only response (K_eff = K + K_G + K_SP at 7.56 rpm) of every candidate load of the T18 tools.

Same K_eff, same metric as ``tools/l5_t6_settled.py`` (ring-axis flap/edge, equal-weight and skin-weight twist), so the
``t6`` and ``t3_zero_az0`` rows reproduce the L31 ``aero_t6`` (-3.204 deg) and ``aero_t3`` (-2.253 deg) rows.
Reads ``loads.npz`` (l5_t18_loads.py) and ``rebuild.npz`` (l5_t18_rebuild.py); writes ``<dir>/twist.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

import l5_oneway_bd_parity as par  # noqa: E402
import twist_p1 as tp  # noqa: E402

# (label, source npz, key)
CASES = [
    ("t6_production_dump_mean", "loads", "t6"),
    ("t3_aerodyn_wall_zero_az0", "loads", "aerodyn_wall_az0"),
    ("aerodyn_wall_zero_azmean", "loads", "aerodyn_wall_azmean"),
    ("aerodyn_wall_rebuilt_azmean", "rebuild", "t3_rebuilt_azmean"),
    ("aerodyn_minnorm_zero_az0", "loads", "aerodyn_minnorm_az0"),
    ("windio_minnorm_zero_az0", "loads", "windio_minnorm_az0"),
    ("windio_minnorm_zero_azmean", "loads", "windio_minnorm_azmean"),
    ("windio_minnorm_rebuilt_az0", "rebuild", "prod_rebuilt_az0"),
    ("windio_minnorm_rebuilt_azmean_FULLMESH", "rebuild", "prod_rebuilt_azmean"),
    ("windio_minnorm_rebuilt_dumptheta_FULLMESH", "rebuild", "prod_rebuilt_dumptheta_mean"),
    ("windio_wall_zero_az0", "loads", "windio_wall_az0"),
    ("aerodyn_wall_zero_az0_noCm", "loads", "aerodyn_wall_az0_nocm"),
    ("windio_minnorm_zero_az0_noCm", "loads", "windio_minnorm_az0_nocm"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--reference", type=Path, default=REPO / "logs/check/l5_oneway/t11/bem.json")
    args = ap.parse_args()
    t0 = time.time()
    from run_s5_oneway import solve_static

    mesh, props, asm = par.build_shell(20.0)
    dpn = asm.dofs_per_node
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    part = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(par.AD_PRIMARY), par.bem_config(-6.0, -4.0),
                             element_properties=props)  # rings for the section estimators only
    K, free = par.rotating_keff(asm, mesh, par.OMEGA, 3.97)
    srcs = {"loads": np.load(args.dir / "loads.npz"), "rebuild": np.load(args.dir / "rebuild.npz")}
    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    zq = np.array(json.loads(args.reference.read_text())["table"]["z"])
    out = {"z": zq.tolist(), "cases": {}}
    for label, src, key in CASES:
        f3 = srcs[src][key]
        f = np.zeros(asm.dofs_count)
        for d in range(3):
            f[d::dpn] = f3[:, d]
        u = solve_static(asm, K, free, f)
        U = np.column_stack([u[d::dpn] for d in range(3)])
        sh = par.shell_sections(mesh, part, U, kp_bd)
        out["cases"][label] = {"flap": np.interp(zq, sh["z"], sh["u"][:, 1]).tolist(),
                               "edge": np.interp(zq, sh["z"], sh["u"][:, 0]).tolist(),
                               "twist_eq_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_eq"])).tolist(),
                               "twist_skin_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_skin"])).tolist(),
                               "W": 0.5 * float(f @ u), "F": f3.sum(0).tolist()}
        r = out["cases"][label]
        print(f"{label:44s} flap {r['flap'][-1]:8.4f} edge {r['edge'][-1]:7.4f} tw_eq {r['twist_eq_deg'][-1]:7.3f} "
              f"tw_skin {r['twist_skin_deg'][-1]:7.3f}  ({time.time() - t0:.0f} s)", flush=True)
    out["revision"] = tp.revision(REPO)
    (args.dir / "twist.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
