"""T18 follow-up: the non-zero-displacement path of ``BEMFSIParticipant._compute_forces``.

``_compute_forces`` evaluates ``self._bem_solver`` (full blade_aero incl. prebend) only when max|u| < 1e-12;
otherwise it calls ``_rebuild_bem_solver``, whose ``AeroStation(...)`` copies drop ``prebend``.  The production
one-way run always takes the second branch (the structure returns u != 0), the T3 loads the first.  This tool
evaluates both branches on the same participants (the displacement only triggers the branch: deformed_radius /
deformed_twist are false in the production config, so r_def = ref, twist_def = ref) and compares with the
production dumps.

Outputs ``<out>/rebuild.json`` and ``<out>/rebuild.npz``.  ``loads.npz`` of ``l5_t18_loads.py`` supplies the T6 mean.
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
import l5_t18_loads as lo  # noqa: E402
import l5_t6_settled as t6s  # noqa: E402
import twist_p1 as tp  # noqa: E402

TINY = 1e-6  # m; only selects the branch (max|u| >= 1e-12); geometry feedback is off in the production config


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--loads", type=Path, required=True, help="loads.npz of l5_t18_loads.py")
    args = ap.parse_args()
    t0 = time.time()
    omega = par.OMEGA
    mesh, props, asm = par.build_shell(20.0)
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant, build_from_config
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    prod = build_from_config(mesh, lo.production_cfg(), element_properties=None)
    t3 = BEMFSIParticipant(mesh, build_blade_aero_from_aerodyn(par.AD_PRIMARY), par.bem_config(-6.0, -4.0),
                           element_properties=props)
    X = np.asarray(prod._ref_coords, float)
    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    u_small = np.zeros((mesh.node_count, 3))
    u_small[:, 0] = TINY

    def f_at(part, theta, rebuilt):
        part._partner_theta = float(theta)
        return part._compute_forces(u_small if rebuilt else np.zeros_like(u_small))

    # dump reproduction (8 dumps over the last 8 revolutions)
    t, th, U, FA, FI, Xv = t6s.load_dumps(args.results, 8, 2 * np.pi / omega)
    from aeroelast.solvers.fsi.corotational import CoordinateTransforms

    ct = CoordinateTransforms(rotation_axis=par.AXIS, rotation_center=np.array([0.0, 0.0, -3.97]))
    FA_r = np.array([ct.to_rotating(f, a) for f, a in zip(FA, th, strict=True)])
    repro = []
    for i in np.linspace(0, len(t) - 1, 8).astype(int):
        row = {"t": float(t[i]), "theta_deg": float(np.degrees(th[i]) % 360), "sum_dump": FA_r[i].sum(0).tolist()}
        for label, rb in (("zero_disp_branch", False), ("rebuilt_branch", True)):
            f, _ = f_at(prod, th[i], rb)
            row[label] = {"sum": f.sum(0).tolist(), "rel_max_node": float(np.abs(f - FA_r[i]).max() / np.abs(FA_r[i]).max()),
                          "rel_norm": float(np.linalg.norm(f - FA_r[i]) / np.linalg.norm(FA_r[i]))}
        repro.append(row)
    print(f"repro done ({time.time() - t0:.0f} s)", flush=True)

    az = np.radians(np.arange(lo.N_AZ) * 360.0 / lo.N_AZ)
    fields = {"t6": np.load(args.loads)["t6"]}
    for name, part in (("prod", prod), ("t3", t3)):
        for lab, rb in (("zero", False), ("rebuilt", True)):
            fields[f"{name}_{lab}_az0"] = f_at(part, 0.0, rb)[0]
            fields[f"{name}_{lab}_azmean"] = np.mean([f_at(part, a, rb)[0] for a in az], axis=0)
        print(f"{name} done ({time.time() - t0:.0f} s)", flush=True)
    fields["prod_rebuilt_dumptheta_mean"] = t6s.tmean(t, np.array([f_at(prod, a, True)[0] for a in th]))
    # station-level effect of the branch on the BEM result (production polars, azimuth 0)
    bem = {}
    for lab, rb in (("zero", False), ("rebuilt", True)):
        res = f_at(prod, 0.0, rb)[1]
        bem[lab] = {k: np.asarray(getattr(res, k), float).tolist() for k in ("r", "Np", "Tp", "Mp", "alpha", "a")}
        bem[lab].update({"thrust": float(res.thrust), "torque": float(res.torque), "power": float(res.power)})
    prebend = np.asarray(prod._blade_aero.prebend, float)
    tables = {}
    for name, f3 in fields.items():
        rows, F, M, z = par.bin_rows(X, f3, kp_bd)
        tables[name] = {"z": z.tolist(), "Fx_edge": F[:, 0].tolist(), "Fy_flap": F[:, 1].tolist(), "Mz_torsion": M[:, 2].tolist(),
                        "total": f3.sum(0).tolist(), "Mz_total": float(M[:, 2].sum())}
    doc = {"tiny_u": TINY, "dump_reproduction": repro, "bem_level_prod_az0": bem,
           "prod_prebend_m": {"min": float(prebend.min()), "max": float(prebend.max()), "tip": float(prebend[-1])},
           "t3_prebend_m": {"min": float(np.min(t3._blade_aero.prebend)), "max": float(np.max(t3._blade_aero.prebend))},
           "tables": tables, "revision": tp.revision(REPO), "elapsed_s": time.time() - t0}
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out / "rebuild.npz", **fields)
    (args.out / "rebuild.json").write_text(json.dumps(doc, indent=1))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
