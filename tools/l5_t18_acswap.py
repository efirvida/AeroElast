"""T18: isolate the aerodynamic-centre (AC) fraction from the polars in the T3 -> production twist difference.

The T3 tool builds its blade from AeroDyn files with ``aerodynamic_center = 0.25`` everywhere; the production blade
(windIO yaml) carries per-airfoil values (0.5 / 0.316 inboard, 0.25 outboard).  The projector's moment arm
``centroid - AC`` reads it.  Two swaps, azimuth 0, zero-displacement branch, wall-flow off (production projector):

* ``windio_ac025``  windIO polars/chords/twist, AC forced to 0.25 everywhere;
* ``aerodyn_acwindio`` AeroDyn polars/chords/twist, AC per station = windIO's interpolated on r.

Compared with ``windio_minnorm_zero_az0`` / ``aerodyn_minnorm_zero_az0`` of ``l5_t18_twist.py`` (same K_eff, same metric).
Output ``<dir>/acswap.json``.
"""

from __future__ import annotations

import argparse
import dataclasses
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
import twist_p1 as tp  # noqa: E402


def with_ac(ba, ac):
    st = [dataclasses.replace(s, airfoil=dataclasses.replace(s.airfoil, aerodynamic_center=float(a))) for s, a in zip(ba.stations, ac, strict=True)]
    return dataclasses.replace(ba, stations=st)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dir", type=Path, required=True)
    args = ap.parse_args()
    t0 = time.time()
    from aeroelast.models.blade.aerodynamics import load_blade_aero
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant
    from run_s5_oneway import solve_static
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    mesh, props, asm = par.build_shell(20.0)
    dpn = asm.dofs_per_node
    cfg = lo.production_cfg()
    b = cfg["bem"]
    wio = load_blade_aero(cfg["blade_file"], default_re=b["default_re"], neuralfoil_model=b["neuralfoil_model"], hub_radius=b["hub_radius"],
                          n_blades=b["n_blades"], viterna_ar=b["viterna_ar"], viterna_confidence_threshold=b["viterna_confidence_threshold"])
    ad = build_blade_aero_from_aerodyn(par.AD_PRIMARY)
    r_w = np.array([s.r for s in wio.stations])
    ac_w = np.array([s.airfoil.aerodynamic_center for s in wio.stations])
    variants = {"windio_ac025": (with_ac(wio, np.full(len(wio.stations), 0.25)), cfg["bem"]),
                "aerodyn_acwindio": (with_ac(ad, np.interp([s.r for s in ad.stations], r_w, ac_w)), None)}
    K, free = par.rotating_keff(asm, mesh, par.OMEGA, 3.97)
    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    zq = np.array(json.loads((par.REPO / "logs/check/l5_oneway/t11/bem.json").read_text())["table"]["z"])
    out = {"z": zq.tolist(), "cases": {}}
    ref_part = None
    for label, (ba, _) in variants.items():
        bem_cfg = {**cfg["bem"], "omega": cfg["bem"]["omega"]["initial"]} if label.startswith("windio") else par.bem_config(-6.0, -4.0)
        part = BEMFSIParticipant(mesh, ba, bem_cfg, element_properties=None)
        ref_part = ref_part or part
        part._partner_theta = 0.0
        fa, _ = part._compute_forces(np.zeros((mesh.node_count, 3)))
        f = np.zeros(asm.dofs_count)
        for d in range(3):
            f[d::dpn] = fa[:, d]
        u = solve_static(asm, K, free, f)
        U = np.column_stack([u[d::dpn] for d in range(3)])
        sh = par.shell_sections(mesh, ref_part, U, kp_bd)
        out["cases"][label] = {"flap": np.interp(zq, sh["z"], sh["u"][:, 1]).tolist(), "edge": np.interp(zq, sh["z"], sh["u"][:, 0]).tolist(),
                               "twist_eq_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_eq"])).tolist(),
                               "twist_skin_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_skin"])).tolist(), "F": fa.sum(0).tolist()}
        r = out["cases"][label]
        print(f"{label:20s} F {np.round(fa.sum(0)/1e3,2)} flap {r['flap'][-1]:.4f} edge {r['edge'][-1]:.4f} tw_eq {r['twist_eq_deg'][-1]:.3f} ({time.time()-t0:.0f} s)", flush=True)
    out["revision"] = tp.revision(par.REPO)
    (args.dir / "acswap.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
