"""T5 (odd/tasks/gmsh-chord-scaled-mesh.md, S4): one metric group of the chord-scaled gmsh mesh convergence on the IEA-15.

Mesh: tests.support.v2_blade.build_v2_blade, gmsh mesher, windIO v2 YAML, mid-surface. ``--case N<k>`` = chord_divisions k,
span_aspect 1.0 (``N<k>A<ar>``: span_aspect ar), min_band_elements 2, min_band_fraction 1.0; ``--case es0.25`` = the uniform element size 0.25 reference.
Groups (each writes ``<out-dir>/<group>.json``; nothing is asserted, the numbers are the deliverable):
  mesh    build the mesh into the disk cache, report its size
  modal   S-4 rotating modal at 7.56 rpm: 1st/2nd flap, 1st edge, 1st torsion (reused tools/run_s4_rotating_modal.py)
  static  S-2 flap/edge tip cases LC1-LC4 (tip node, tip ring mean), W = 1/2 u.K.u = 1/2 f.u per case, S-7 torsion ratio
          (reused tests/validation/blade/v2 helpers, same code as the tests)
  gc      G-C one-way rated, tip elastic twist and rated thrust/torque/power (tools/l3_torque_sensitivity.solve_rated),
          with the deck property map (production, wall flow) and without it (minimum-norm)
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tools"))

RPM_RATED = 7.56
BD_ROTATING = {"1st flap": 0.532003414, "1st edge": 0.699413448, "2nd flap": 1.508057474}
BD_STATIC_UY = {"LC1_gravity": -2.231483, "LC2_tip_flap": 82.273278, "LC4_distributed": 9.062631}
BD_LC3_UX = 31.903817


def parse_case(case: str) -> dict:
    if case.startswith("es"):
        return {"element_size": float(case[2:])}
    if case.startswith("N"):  # N<k> or N<k>A<ar> (T6: span_aspect ar, default 1.0)
        n, _, ar = case[1:].partition("A")
        return {"element_size": 0.25, "chord_divisions": float(n), "span_aspect": float(ar) if ar else 1.0,
                "min_band_elements": 2, "min_band_fraction": 1.0}
    raise SystemExit(f"--case must be N<k>, N<k>A<ar> or es<h>, got {case!r}")


def build(case: str):
    from tests.support.v2_blade import build_v2_blade

    t0 = time.time()
    mesh, props = build_v2_blade(skin_reference="mid", mesher="gmsh", **parse_case(case))
    return mesh, props, time.time() - t0


def mesh_info(mesh, props, t_mesh):
    return {"nodes": int(mesh.node_count), "elements": int(mesh.elements_count), "dofs": int(6 * mesh.node_count),
            "mesh_s": round(t_mesh, 1)}


def group_mesh(case, args):
    mesh, props, t = build(case)
    return mesh_info(mesh, props, t)


def _assembler(mesh, props):
    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily

    return MeshAssembler(mesh, {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0),
                                             "properties": props}})


def group_modal(case, args):
    from run_s4_rotating_modal import assemble_rotating, classify_modes, omega_from_rpm

    mesh, props, t_mesh = build(case)
    asm = _assembler(mesh, props)
    t0 = time.time()
    freqs, modes, free = assemble_rotating(asm, mesh, omega_from_rpm(args.rpm), n_modes=args.n_modes)
    t_eig = time.time() - t0
    labels = classify_modes(mesh, modes, free)
    out = {"mesh": mesh_info(mesh, props, t_mesh), "rpm": args.rpm, "n_modes": args.n_modes, "eigen_s": round(t_eig, 1),
           "all_modes": [{"f_hz": float(f), "label": lab} for f, lab in zip(freqs, labels, strict=True)], "families": {}}
    for name, direction, order in [("1st flap", "flap", 1), ("1st edge", "edge", 1), ("2nd flap", "flap", 2),
                                   ("1st torsion", "torsion", 1)]:
        cand = sorted(float(f) for f, lab in zip(freqs, labels, strict=True) if lab == direction)
        f = cand[order - 1] if len(cand) >= order else None
        ref = BD_ROTATING.get(name)
        out["families"][name] = {"f_hz": f, "beamdyn_hz": ref, "rel_to_beamdyn": (f / ref - 1.0) if (f and ref) else None}
    return out


def group_static(case, args):
    from scipy.sparse import coo_matrix

    import tests.validation.blade.v2.test_v2_iea15mw_s2_static_prescribed as s2
    import tests.validation.blade.v2.test_v2_iea15mw_s7_torsion as s7
    from run_s7_torsion import beam_twist_profile_deg, section_twists_deg

    mesh, props, t_mesh = build(case)
    asm = _assembler(mesh, props)
    n = asm.dofs_count
    rows, cols, vals = asm._rust.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
    tip = s2._tip_idxs(mesh)
    out = {"mesh": mesh_info(mesh, props, t_mesh), "s2": {}, "s7": {}}
    for name, load in s2._TIP_CASES.items():
        t0 = time.time()
        u = np.asarray(s2._run(mesh, props, load), dtype=np.float64)
        w = 0.5 * float(u @ (K @ u))
        node = s2.shell_tip_displacement(mesh, u)
        ring = u.reshape(-1, 6)[tip, :3].mean(axis=0)
        rec = {"tip_node_u_m": [float(v) for v in node], "tip_ring_mean_u_m": [float(v) for v in ring],
               "W_joule": w, "solve_s": round(time.time() - t0, 1)}
        if load["type"] == "tip":  # W = 1/2 F . mean(u over the loaded ring): independent check of u.K.u
            rec["W_tip_force_check"] = 0.5 * float(np.dot(load["value"], ring))
        ref = BD_STATIC_UY.get(name)
        if ref is not None:
            rec["uy_beamdyn_m"], rec["uy_rel_to_beamdyn"] = ref, float(node[1] / ref - 1.0)
        if name == "LC3_tip_edge":
            rec["ux_beamdyn_m"], rec["ux_rel_to_beamdyn"] = BD_LC3_UX, float(node[0] / BD_LC3_UX - 1.0)
        out["s2"][name] = rec
    t0 = time.time()
    u = np.asarray(s7.run_torsion_case_v2(mesh, props, s7.MOMENT_NM), dtype=np.float64)
    zs, theta = section_twists_deg(mesh, u)
    beam = beam_twist_profile_deg(s7.MOMENT_NM, zs)
    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    out["s7"] = {"ratio_mean_0.3_0.9": float(np.nanmean(theta[mid] / beam[mid])), "moment_nm": s7.MOMENT_NM,
                 "W_joule": 0.5 * float(u @ (K @ u)), "solve_s": round(time.time() - t0, 1)}
    return out


def group_gc(case, args):
    import l3_torque_sensitivity as l3

    spec = parse_case(case)
    out = {}
    for label, wall_flow in (("wall_flow_production", True), ("min_norm", False)):
        t0 = time.time()
        participant, disp = l3.solve_rated(spec["element_size"], "gmsh", wall_flow=wall_flow,
                                           chord_divisions=spec.get("chord_divisions"))
        r_ref, tw_ref = participant._ref_r, participant._ref_twist
        r_def, tw_def = participant._compute_deformed_geometry(disp)
        d_twist = tw_def - tw_ref
        rigid, deformed = l3.evaluate(participant, r_ref, tw_ref), l3.evaluate(participant, r_def, tw_def)
        out[label] = {"tip_twist_deg": float(np.rad2deg(d_twist[-1])),
                      "max_abs_twist_deg": float(np.rad2deg(np.max(np.abs(d_twist)))),
                      "rigid": {k: float(v) for k, v in rigid.items()}, "deformed": {k: float(v) for k, v in deformed.items()},
                      "deformed_vs_rigid": {q: float(deformed[q] / rigid[q] - 1.0) for q in ("thrust", "torque", "power")},
                      "n_nodes": int(disp.shape[0]), "wall_s": round(time.time() - t0, 1)}
    return out


GROUPS = {"mesh": group_mesh, "modal": group_modal, "static": group_static, "gc": group_gc}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("group", choices=sorted(GROUPS))
    ap.add_argument("--case", required=True, help="N10, N15, ... (chord_divisions) or es0.25 (uniform)")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--rpm", type=float, default=RPM_RATED)
    ap.add_argument("--n-modes", type=int, default=20)
    args = ap.parse_args()
    t0 = time.time()
    result = GROUPS[args.group](args.case, args)
    result = {"case": args.case, "group": args.group, "spec": parse_case(args.case), **result,
              "wall_s": round(time.time() - t0, 1),
              "max_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6, 2)}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    path = args.out_dir / f"{args.group}.json"
    path.write_text(json.dumps(result, indent=2))
    print(f"wrote {path} wall {result['wall_s']} s maxrss {result['max_rss_gb']} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
