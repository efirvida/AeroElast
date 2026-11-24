"""T18 (odd/tasks/l5-one-way.md): why the T6 mean aero load and the T3 load give different twist on the same K_eff.

Stage 1 (this tool, no structural solve): nodal aero load vectors of

* ``t6``       time mean over whole revolutions of the production run's F_AERO (rotated to the rotating frame);
* ``t3``       T3 = ``BEMFSIParticipant._compute_forces(zeros)`` of ``l5_oneway_bd_parity.bem_config`` (AeroDyn
               polars, ``element_properties`` given, azimuth 0);
* ``prod``     the participant ``build_from_config`` makes from ``cases/l5_oneway_rated/fluid/fluid.yaml`` through the
               runner's own ``_effective_element_properties`` (windIO polars, wall-flow realisation OFF);
* the 2x2 grid {AeroDyn, windIO} x {wall-flow on, off}, each at azimuth 0, azimuth-mean (36 azimuths) and the
  production dumps' own thetas, and each with the pitching moment Cm removed (``Mp = None``).

Everything is binned with ``l5_oneway_bd_parity.bin_rows`` (0.5 m strips): resultant force (edge = mesh x,
flap = mesh y), strip moment about the reference axis (M_z = torsion), and the cumulative section torsion from the tip.
Output: ``<out>/loads.npz`` (nodal vectors) and ``<out>/loads.json`` (strip tables and totals).
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
import l5_t6_settled as t6s  # noqa: E402
import twist_p1 as tp  # noqa: E402

FLUID_YAML = par.FLUID_YAML
N_AZ = par.N_AZ
production_cfg = par.production_cfg


def participants(mesh, props):
    """{(polars, wall): participant}; ('windio', False) is the production one."""
    from aeroelast.cli.run_bem_fsi import _effective_element_properties
    from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant, build_from_config
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    cfg = production_cfg()
    out = {}
    for wall in (True, False):
        ep = props if wall else None
        out[("aerodyn", wall)] = BEMFSIParticipant(
            mesh, build_blade_aero_from_aerodyn(par.AD_PRIMARY), par.bem_config(-6.0, -4.0), element_properties=ep)
        out[("windio", wall)] = build_from_config(mesh, cfg, element_properties=ep)
    # the production runner's decision, asserted rather than assumed
    assert _effective_element_properties(cfg, props) is None, "production wall-flow realisation is ON"
    return out


def _arm_diff(pa, pb):
    """Max |ac_offset_a - ac_offset_b| [m] after interpolating b onto a's strip centres (stations differ: 53 vs 50)."""
    ca, cb = np.array([s.r_center for s in pa._projector._strips]), np.array([s.r_center for s in pb._projector._strips])
    oa, ob = np.array(pa._projector._strip_ac_offsets), np.array(pb._projector._strip_ac_offsets)
    ib = np.column_stack([np.interp(ca, cb, ob[:, k]) for k in range(3)])
    return float(np.abs(oa - ib).max())


def forces_at(part, theta, cm=True):
    part._partner_theta = float(theta)
    fa, res = part._compute_forces(np.zeros((part._n_nodes, 3)))
    if not cm:
        fa = part._projector.project(dataclasses.replace(res, Mp=None))
    return fa, res


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", type=Path, required=True, help="T6 results/ directory")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-rev", type=int, default=8)
    args = ap.parse_args()
    t0 = time.time()
    omega = par.OMEGA
    mesh, props, asm = par.build_shell(20.0)
    parts = participants(mesh, props)
    prod = parts[("windio", False)]
    X = np.asarray(prod._ref_coords, float)
    kp_bd = tp.read_keypoints((par.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]

    # configuration diff, production (parsed) vs the T3 tool
    diff = {"production_bem_solver_kwargs": prod._bem_solver_kwargs,
            "t3_bem_solver_kwargs": parts[("aerodyn", True)]._bem_solver_kwargs,
            "production_omega_rad_s": prod._current_omega, "t3_omega_rad_s": parts[("aerodyn", True)]._current_omega,
            "production_wind": prod._bem_cfg.get("wind_speed"), "t3_wind": parts[("aerodyn", True)]._bem_cfg.get("wind_speed"),
            "production_azimuth_offset": prod._azimuth_offset, "t3_azimuth_offset": parts[("aerodyn", True)]._azimuth_offset,
            "production_deformed_twist": prod._bem_cfg.get("deformed_twist"), "production_deformed_radius": prod._bem_cfg.get("deformed_radius"),
            "production_projector_element_properties_is_None": prod._projector._element_properties is None,
            "t3_projector_element_properties_is_None": parts[("aerodyn", True)]._projector._element_properties is None,
            "production_n_stations": int(len(prod._ref_r)), "t3_n_stations": int(len(parts[("aerodyn", True)]._ref_r)),
            "production_r_range": [float(prod._ref_r[0]), float(prod._ref_r[-1])],
            "t3_r_range": [float(parts[("aerodyn", True)]._ref_r[0]), float(parts[("aerodyn", True)]._ref_r[-1])],
            "projector_arm_prod_vs_t3_max_abs_m": _arm_diff(prod, parts[("aerodyn", True)]),
            "include_pitching_moment_read_by_participant": False}

    # T6 mean (rotated with the dumps' theta) and per-dump reproduction by the production participant
    t, th, U, FA, FI, Xv = t6s.load_dumps(args.results, args.n_rev, 2 * np.pi / omega)
    from aeroelast.solvers.fsi.corotational import CoordinateTransforms

    ct = CoordinateTransforms(rotation_axis=par.AXIS, rotation_center=np.array([0.0, 0.0, -3.97]))
    FA_r = np.array([ct.to_rotating(f, a) for f, a in zip(FA, th, strict=True)])
    fields = {"t6": t6s.tmean(t, FA_r)}
    print(f"dumps {len(t)} loaded ({time.time() - t0:.0f} s)", flush=True)

    repro = []
    for i in np.linspace(0, len(t) - 1, 8).astype(int):
        f, _ = forces_at(prod, th[i])
        repro.append({"t": float(t[i]), "theta_deg": float(np.degrees(th[i]) % 360),
                      "rel_diff_max_node": float(np.abs(f - FA_r[i]).max() / np.abs(FA_r[i]).max()),
                      "sum_prod": f.sum(0).tolist(), "sum_dump": FA_r[i].sum(0).tolist()})
    # azimuth offset scan on one dump: does the dump match theta + delta?
    i = len(t) // 2
    scan = []
    for dd in (-30, -15, -7.5, -3.75, 0.0, 3.75, 7.5, 15, 30):
        f, _ = forces_at(prod, th[i] + np.radians(dd))
        scan.append([dd, float(np.abs(f - FA_r[i]).max() / np.abs(FA_r[i]).max())])

    az = np.radians(np.arange(N_AZ) * 360.0 / N_AZ)
    az_curve = {}
    for key, part in parts.items():
        name = f"{key[0]}_{'wall' if key[1] else 'minnorm'}"
        fz, _ = forces_at(part, 0.0)
        fields[f"{name}_az0"] = fz
        fields[f"{name}_azmean"] = np.mean([forces_at(part, a)[0] for a in az], axis=0)
        fields[f"{name}_az0_nocm"] = forces_at(part, 0.0, cm=False)[0]
        fields[f"{name}_azmean_nocm"] = np.mean([forces_at(part, a, cm=False)[0] for a in az], axis=0)
        print(f"{name} done ({time.time() - t0:.0f} s)", flush=True)
    fields["prod_dumptheta_mean"] = t6s.tmean(t, np.array([forces_at(prod, a)[0] for a in th]))
    # production participant azimuth curve of the totals (rotating frame)
    az_curve = [{"azimuth_deg": float(np.degrees(a)), "F": forces_at(prod, a)[0].sum(0).tolist()} for a in az]
    # BEM-level Np/Tp/Mp, windIO vs AeroDyn at azimuth 0 (per station)
    bem_lvl = {}
    for key in (("aerodyn", True), ("windio", False)):
        _, res = forces_at(parts[key], 0.0)
        bem_lvl[key[0]] = {k: np.asarray(getattr(res, k), float).tolist() for k in ("r", "Np", "Tp", "Mp", "alpha", "cl", "cd", "a", "ap")
                           if getattr(res, k) is not None}
        bem_lvl[key[0]].update({"thrust": float(res.thrust), "torque": float(res.torque), "power": float(res.power)})

    # strip tables
    tables = {}
    root_ref = None
    for name, f3 in fields.items():
        rows, F, M, z = par.bin_rows(X, f3, kp_bd)
        r0 = par.ref_axis_mesh(kp_bd, z)
        # section torsion about the local reference axis, from the tip: M_z strip + (r0_k - r0(z)) x F_k, z component
        T = np.zeros(len(z))
        for j in range(len(z)):
            sel = np.arange(j, len(z))
            d = r0[sel] - r0[j]
            T[j] = float((M[sel, 2] + d[:, 0] * F[sel, 1] - d[:, 1] * F[sel, 0]).sum())
        tables[name] = {"z": z.tolist(), "Fx_edge": F[:, 0].tolist(), "Fy_flap": F[:, 1].tolist(), "Mx": M[:, 0].tolist(),
                        "My": M[:, 1].tolist(), "Mz_torsion": M[:, 2].tolist(), "T_cum_from_tip": T.tolist(),
                        "total": f3.sum(0).tolist(), "Mz_total": float(M[:, 2].sum()), "T_root": float(T[0]),
                        "Mroot_flap_about_root": float(np.cross(X - par.ref_axis_mesh(kp_bd, np.array([0.0]))[0], f3).sum(0)[0])}
    doc = {"results": str(args.results), "n_dumps": len(t), "window": [float(t[0]), float(t[-1])], "n_az": N_AZ,
           "config_diff": diff, "dump_reproduction": repro, "theta_offset_scan_max_rel": scan,
           "prod_az_curve": az_curve, "bem_level": bem_lvl, "tables": tables, "revision": tp.revision(REPO),
           "elapsed_s": time.time() - t0}
    args.out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.out / "loads.npz", **fields)
    (args.out / "loads.json").write_text(json.dumps(doc, indent=1))
    print("done", time.time() - t0)


if __name__ == "__main__":
    main()
