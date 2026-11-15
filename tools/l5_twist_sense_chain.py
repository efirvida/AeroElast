"""L5 T10: ONE sign chain for the pretwist sense of the shell blade (odd/tasks/l5-one-way.md, L10).

Question: BeamDyn places a section at -kp_twist about +z (BeamDyn_IO.f90:898); the shell mesh places the chord at +BlTwist.
(a) mirrored pretwist defect, (b) consistent mirror world, or (c) a mis-mapped BeamDyn frame?

The chain is built from physical definitions and measured in the generated gmsh N20 mesh (no production code is changed):

  1. physical blade: rotor axis, upwind, spin sense, LE-first motion, twist positive to feather (AeroDyn / BeamDyn docs)
  2. mesh frame: LE side, suction (cambered) side, span, prebend side, measured per ring
  3. handedness of (LE, suction, span) and of the spin/loads: proper or mirrored
  4. angle of attack the MESH chord sees in the BEM inflow vs the BEM's own alpha (phi - theta): the decisive number
  5. pitch-to-feather direction in the mesh, thrust and torque signs of the production BEM nodal loads
  6. AeroDyn Fn/Ft frame vs tools/run_s5_oneway.apply_ad_loads
  7. v1 (NuMAD) vs v2 (windIO frames) chord angle, geometry only

Pure helpers (unit-tested by tests/software/test_l5_twist_sense_chain.py) sit above ``# --- views``.

  scripts/aeroenv.sh python tools/l5_twist_sense_chain.py --out logs/check/l5_oneway/t10 [--work $SCRATCH/tmp/l5-t10]

Mesh frame (measured): +x leading edge, +y suction = downwind, +z span.  Angles are right-handed about +z; ``psi`` is the
direction of the chord measured from the trailing edge to the leading edge, from +x.
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

WIND = 10.59  # m/s, rated
RPM = 7.56
OMEGA = RPM * 2.0 * np.pi / 60.0
HUB_RADIUS = 3.97


# --- pure helpers (unit-tested) ----------------------------------------------------------------------------------


def handedness(a, b, c):
    """+1 when (a, b, c) is a right-handed triad, -1 when left-handed (sign of the triple product)."""
    return int(np.sign(np.linalg.det(np.array([a, b, c], float))))


def rotate_z(v, delta_deg):
    """Right-handed rotation of a 2-D (or 3-D) vector about +z by ``delta_deg``."""
    v = np.asarray(v, float)
    a = np.radians(delta_deg)
    c, s = np.cos(a), np.sin(a)
    out = v.copy()
    out[..., 0] = c * v[..., 0] - s * v[..., 1]
    out[..., 1] = s * v[..., 0] + c * v[..., 1]
    return out


def _farthest_pair(xy):
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2)
    i, j = np.unravel_index(np.argmax(d), d.shape)
    return xy[i], xy[j]


def section_from_ring(xy):
    """Chord geometry of a section ring (n, 2): LE, TE, psi (TE->LE angle from +x, deg), suction normal, camber/c.

    The chord is the farthest pair; the LE is the end whose outline is thicker 3-10 % of the chord from the end (airfoil
    nose against thin trailing edge); the suction side is the side the mean outline offset lies on (camber).
    """
    xy = np.asarray(xy, float)
    a, b = _farthest_pair(xy)
    c = float(np.linalg.norm(b - a))
    u = (b - a) / c
    n = np.array([-u[1], u[0]])
    s = (xy - a) @ u / c
    o = (xy - a) @ n / c

    def width(lo, hi):
        m = (s > lo) & (s < hi)
        return float(o[m].max() - o[m].min()) if m.sum() > 1 else 0.0

    wa, wb = width(0.03, 0.10), width(0.90, 0.97)
    le, te = (a, b) if wa > wb else (b, a)
    d = le - te
    psi = float(np.degrees(np.arctan2(d[1], d[0])))
    camber = float(np.mean(o))
    n_s = n if camber >= 0.0 else -n
    return {"le": le, "te": te, "chord": c, "psi_deg": psi, "suction_normal": n_s, "camber_over_c": abs(camber),
            "end_width_ratio": float(max(wa, wb) / max(min(wa, wb), 1e-12))}


def blade_velocity(omega_vec, r_rel):
    """Velocity of the blade point at ``r_rel`` from the hub centre: Omega x r."""
    return np.cross(np.asarray(omega_vec, float), np.asarray(r_rel, float))


def inflow(v_blade, v_wind, a=0.0, ap=0.0):
    """Relative velocity seen by the section: wind (1 - a) minus blade velocity (1 + ap)."""
    return (1.0 - a) * np.asarray(v_wind, float) - (1.0 + ap) * np.asarray(v_blade, float)


def angle_of_attack_deg(le_to_te, suction_normal, v_rel_xy):
    """Signed angle of attack: atan2(v . n_suction, v . c_hat); positive when the flow meets the pressure side."""
    return float(np.degrees(np.arctan2(np.dot(v_rel_xy, suction_normal), np.dot(v_rel_xy, le_to_te))))


def ad_to_rotor_plane(fn, ft, theta):
    """AeroDyn Fn/Ft (airfoil frame) -> (normal to the rotor plane, tangential in it); AeroDyn_IO.f90:486-487 inverted."""
    ct, st = np.cos(theta), np.sin(theta)
    return np.array([fn * ct - ft * st, fn * st + ft * ct])


def ad_force_in_mesh(fn, ft, theta, sense):
    """(f_x, f_y) in the mesh frame of an AeroDyn Fn/Ft pair for a section whose chord sits at ``sense * theta``.

    ``sense = -1`` is the physical blade (BlTwist positive to feather: TE->LE direction at -theta about +z);
    ``sense = +1`` is what ``tools/run_s5_oneway.apply_ad_loads`` builds (n_hat = (-sin, cos)).
    """
    a = sense * theta
    return fn * np.array([-np.sin(a), np.cos(a)]) + ft * np.array([np.cos(a), np.sin(a)])


def torque_about_y(X, f, hub_radius):
    """Torque about the rotor axis (+y through (0, 0, -hub)) of nodal forces f at X: sum (z + hub) f_x - x f_z."""
    X, f = np.asarray(X, float), np.asarray(f, float)
    return float(np.sum((X[:, 2] + hub_radius) * f[:, 0] - X[:, 0] * f[:, 2]))


# --- views (need the mesh and the production BEM; not unit-tested) -----------------------------------------------


def _work_mesh(work):
    """Mesh, properties, assembler of the reference gmsh N20 AR 1 mid-surface blade (same as l5_oneway_bd_parity)."""
    import l5_oneway_bd_parity as P

    return P.build_shell(20)


def _participant_and_loads(mesh, props, asm, tilt, precone, work):
    import l5_oneway_bd_parity as P

    tag = f"bem_t{tilt:g}_p{precone:g}"
    f_path = Path(work) / f"{tag}_f.npy"
    f, part = P.shell_loads("bem", asm, mesh, props, tilt, precone)
    np.save(f_path, f)
    return part, f


def view_geometry(mesh, part, out):
    """Links 2-3: LE/TE/psi/suction side/prebend per ring, handedness, vs BlTwist."""
    X = np.array([n.coords for n in mesh.nodes], float)
    import l5_oneway_bd_parity as P
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    aero = build_blade_aero_from_aerodyn(P.AD_PRIMARY)
    rows = []
    for k, ring in enumerate(part._strip_ring_indices):
        ring = np.asarray(ring)
        if len(ring) < 8:
            continue
        pts = X[ring]
        z = float(pts[:, 2].mean())
        sec = section_from_ring(pts[:, :2])
        tw = float(np.degrees(np.interp(z + aero.hub_radius, aero.r, aero.twist)))
        rows.append({"k": k, "z": z, "chord": sec["chord"], "le_x": float(sec["le"][0]), "le_y": float(sec["le"][1]),
                     "te_x": float(sec["te"][0]), "te_y": float(sec["te"][1]), "psi_deg": sec["psi_deg"],
                     "bl_twist_deg": tw, "camber_over_c": sec["camber_over_c"], "end_width_ratio": sec["end_width_ratio"],
                     "suction_nx": float(sec["suction_normal"][0]), "suction_ny": float(sec["suction_normal"][1]),
                     "centre_y": float(pts[:, 1].mean())})
    good = [r for r in rows if r["camber_over_c"] > 0.004 and r["end_width_ratio"] > 2.0 and 10.0 < r["z"] < 110.0]
    summary = {
        "n_rings": len(rows), "n_good": len(good),
        "le_on_plus_x": all(r["le_x"] > r["te_x"] for r in good),
        "suction_on_plus_y": all(r["suction_ny"] > 0.9 for r in good),
        "psi_minus_bltwist_deg_max_abs": max(abs(r["psi_deg"] - r["bl_twist_deg"]) for r in good),
        "psi_plus_bltwist_deg_max_abs": max(abs(r["psi_deg"] + r["bl_twist_deg"]) for r in good),
        "tip_ring_centre_y_m": rows[-1]["centre_y"],
        "handedness_le_suction_span": handedness([1, 0, 0], [0, 1, 0], [0, 0, 1]),
    }
    out["geometry"] = {"rows": rows, "summary": summary}
    return rows


def view_aoa(part, rows, out):
    """Link 4: the angle of attack the MESH chord sees in the BEM inflow vs the BEM's own alpha."""
    res = part._bem_solver.compute(WIND, RPM, 0.0, azimuth=0.0)
    spin = np.array([0.0, OMEGA, 0.0])
    wind = np.array([0.0, WIND, 0.0])
    out_rows = []
    for r in rows:
        k = r["k"]
        a_k, ap_k = float(res.a[k]), float(res.ap[k])
        r_rel = np.array([0.0, 0.0, r["z"] + HUB_RADIUS])  # section on the span axis, hub centre at z = -hub
        v = inflow(blade_velocity(spin, r_rel), wind, a_k, ap_k)
        sec_psi = r["psi_deg"]
        c_hat = -np.array([np.cos(np.radians(sec_psi)), np.sin(np.radians(sec_psi))])  # LE -> TE
        n_s = np.array([r["suction_nx"], r["suction_ny"]])
        alpha_mesh = angle_of_attack_deg(c_hat, n_s, v[:2])
        theta = float(res.twist_deg[k]) if getattr(res, "twist_deg", None) is not None else r["bl_twist_deg"]
        out_rows.append({"k": k, "z": r["z"], "psi_deg": sec_psi, "bl_twist_deg": theta,
                         "alpha_bem_deg": float(res.alpha[k]), "phi_bem_deg": float(res.alpha[k]) + theta,
                         "phi_geom_deg": float(np.degrees(np.arctan2(v[1], -v[0]))),
                         "alpha_mesh_chord_deg": alpha_mesh,
                         "alpha_if_chord_at_minus_twist_deg": angle_of_attack_deg(
                             -rotate_z(np.array([1.0, 0.0]), -theta), rotate_z(np.array([0.0, 1.0]), -theta), v[:2])})
        # pitch test: rotate the mesh section by +1 deg about +z
        c1 = rotate_z(c_hat, 1.0)
        n1 = rotate_z(n_s, 1.0)
        out_rows[-1]["d_alpha_per_plus_1deg_about_z"] = angle_of_attack_deg(c1, n1, v[:2]) - alpha_mesh
    out["aoa"] = {"rows": out_rows,
                  "max_abs_alpha_mesh_minus_bem_deg": max(abs(q["alpha_mesh_chord_deg"] - q["alpha_bem_deg"]) for q in out_rows
                                                          if 10.0 < q["z"] < 110.0),
                  "max_abs_alpha_minus_twist_minus_bem_deg": max(
                      abs(q["alpha_if_chord_at_minus_twist_deg"] - q["alpha_bem_deg"]) for q in out_rows
                      if 10.0 < q["z"] < 110.0),
                  "twist_mesh_to_bem": float(part._twist_mesh_to_bem)}


def view_loads(mesh, f, tag, out):
    """Link 5: thrust and torque signs of a production BEM nodal force vector."""
    X = np.array([n.coords for n in mesh.nodes], float)
    f3 = np.column_stack([f[d::6] for d in range(3)])
    F = f3.sum(0)
    torque = torque_about_y(X, f3, HUB_RADIUS)
    out.setdefault("loads", {})[tag] = {"F_total_N": F.tolist(), "torque_about_plus_y_Nm": torque,
                                        "power_per_blade_W": torque * OMEGA, "spin_vector_rad_s": [0.0, OMEGA, 0.0]}


def view_ad_frame(out):
    """Link 6: AeroDyn Fn/Ft -> mesh totals for the physical (-theta) and the apply_ad_loads (+theta) frames."""
    loads = np.loadtxt(REPO / "tests/IEA15MW/reference/s5_ad_blade1_loads.csv", delimiter=",", skiprows=1)
    r, fn, ft, tw = loads.T
    th = np.radians(tw)
    res = {}
    for name, sense in (("physical_minus_theta", -1), ("apply_ad_loads_plus_theta", 1)):
        fx, fy = np.array([ad_force_in_mesh(a, b, t, sense) for a, b, t in zip(fn, ft, th, strict=True)]).T
        res[name] = {"Fx_N": float(np.trapezoid(fx, r)), "Fy_N": float(np.trapezoid(fy, r)),
                     "torque_Nm": float(np.trapezoid(fx * (r + HUB_RADIUS), r))}
        res[name]["power_3_blades_MW"] = 3.0 * OMEGA * res[name]["torque_Nm"] / 1e6
    res["sum_fn_N"] = float(np.trapezoid(fn, r))
    res["sum_ft_N"] = float(np.trapezoid(ft, r))
    out["ad_frame"] = res


def view_v1_v2(rows, out):
    """Link 7: the legacy NuMAD (v1) and the windIO-frames (v2) chord angle, geometry only."""
    from aeroelast.models.blade.numad.io.yaml_to_blade import yaml_to_blade
    from aeroelast.models.blade.numad.objects.blade import Blade

    b = Blade()
    yaml_to_blade(b, str(REPO / "tests/IEA-15-240-RWT.yaml"))
    g = b.geometry
    v1 = []
    for k in np.linspace(5, len(g.ispan) - 5, 8).astype(int):
        xy = g.coordinates[:, :2, k]
        le = xy[g.LEindex]
        te = 0.5 * (xy[0] + xy[-1])
        d = le - te
        v1.append({"z": float(g.ispan[k]), "twist_deg": float(g.idegreestwist[k]),
                   "psi_deg": float(np.degrees(np.arctan2(d[1], d[0])))})
    out["v1_v2"] = {"v1_numad_rotorspin": int(b.definition.rotorspin), "v1": v1,
                    "v1_psi_equals_minus_twist": all(abs(q["psi_deg"] + q["twist_deg"]) < 0.05 for q in v1),
                    "v2_gmsh_psi_equals_plus_twist": all(abs(r["psi_deg"] - r["bl_twist_deg"]) < 0.4
                                                         for r in rows if 10.0 < r["z"] < 110.0)}


def flipped_yaml(src, dst):
    """Scratch copy of the windIO v2 deck with the twist values negated (the physical chord sense in the v2 frame code)."""
    import yaml

    deck = yaml.safe_load(Path(src).read_text())
    tw = deck["components"]["blade"]["outer_shape"]["twist"]
    tw["values"] = [-float(v) for v in tw["values"]]
    Path(dst).write_text(yaml.safe_dump(deck, sort_keys=False))
    return dst


def view_flipped(out_dir, work):
    """Link 8 (falsifier): the SAME parked AeroDyn loads (keep x / keep y) on the stock mesh and on a mesh whose chord sits at
    -BlTwist (scratch deck with the twist negated; the mesher is not touched), against BeamDyn lin recorded by L8.
    """
    import l5_oneway_bd_parity as P
    import twist_p1 as tp
    from run_s5_oneway import solve_static
    from tests.support.v2_blade import DEFAULT_YAML, build_v2_blade

    from aeroelast.core.assembler import MeshAssembler
    from aeroelast.elements import ElementFamily

    work = Path(work)
    alt = flipped_yaml(DEFAULT_YAML, work / "IEA-15-240-RWT_twist_negated.yaml")
    t0 = time.time()
    mesh, props = build_v2_blade(0.25, skin_reference="mid", mesher="gmsh", es_range=(0.05, 0.25), chord_divisions=20,
                                 yaml_file=Path(alt))
    model = {"elements": {"element_family": ElementFamily.SHELL, "span_direction": (0.0, 0.0, 1.0), "properties": props}}
    asm = MeshAssembler(mesh, model)
    print(f"flipped mesh {mesh.node_count} nodes ({time.time() - t0:.0f} s)", flush=True)
    dpn = asm.dofs_per_node
    f_all, part = P.shell_loads("ad", asm, mesh, props, -6.0, -4.0)
    K, free = P.rotating_keff(asm, mesh, 0.0, HUB_RADIUS)
    kp_bd = tp.read_keypoints((P.DECK_DIR / "IEA-15-240-RWT_BeamDyn.dat").read_text())[:, :3]
    res = {}
    for keep in ("x", "y"):
        f = f_all.copy()
        for d, a in enumerate("xyz"):
            if a not in keep:
                f[d::dpn] = 0.0
        u = solve_static(asm, K, free, f)
        U = np.column_stack([u[d::dpn] for d in range(3)])
        sh = P.shell_sections(mesh, part, U, kp_bd)
        old = json.loads((Path(out_dir).parent / f"ad-parked-{keep}.json").read_text())["table"]
        zq = np.array(old["z"])
        res[keep] = {"z": zq.tolist(),
                     "flipped_flap": np.interp(zq, sh["z"], sh["u"][:, 1]).tolist(),
                     "flipped_edge": np.interp(zq, sh["z"], sh["u"][:, 0]).tolist(),
                     "flipped_twist_skin_deg": np.rad2deg(np.interp(zq, sh["z"], sh["twist_skin"])).tolist(),
                     "stock_flap": old["shell"]["flap"], "stock_edge": old["shell"]["edge"],
                     "stock_twist_skin_deg": old["shell"]["twist_skin_deg"],
                     "bd_flap": old["beamdyn_x0.01"]["flap"], "bd_edge": old["beamdyn_x0.01"]["edge"],
                     "bd_twist_deg": old["beamdyn_x0.01"]["twist_deg"]}
        print(keep, "done", time.time() - t0, flush=True)
    (Path(out_dir) / "flipped.json").write_text(json.dumps(res, indent=1))
    lines = ["# T10 falsifier: parked AD loads, stock mesh vs chord-at--BlTwist mesh vs BeamDyn lin (L8 numbers)", ""]
    for keep, r in res.items():
        lines.append(f"## load kept: {keep} (mesh axis)\n\n| z | quantity | stock shell | twist-negated shell | BeamDyn lin |\n|---:|---|---:|---:|---:|")
        for i, z in enumerate(r["z"]):
            if z < 60:
                continue
            for q in ("flap", "edge"):
                lines.append(f"| {z:.1f} | {q} [m] | {r['stock_' + q][i]:+.3f} | {r['flipped_' + q][i]:+.3f} | {r['bd_' + q][i]:+.3f} |")
            lines.append(f"| {z:.1f} | twist [deg] | {r['stock_twist_skin_deg'][i]:+.2f} | {r['flipped_twist_skin_deg'][i]:+.2f} | "
                         f"{r['bd_twist_deg'][i]:+.2f} |")
        lines.append("")
    (Path(out_dir) / "flipped.md").write_text("\n".join(lines))
    print("\n".join(lines))


def markdown(out):
    g, a = out["geometry"]["summary"], out["aoa"]
    L = ["# L5 T10 sign chain: measured links", ""]
    L.append(f"- mesh rings (good) {g['n_good']}/{g['n_rings']}: LE on +x {g['le_on_plus_x']}, suction on +y {g['suction_on_plus_y']}, "
             f"handedness(LE, suction, span) {g['handedness_le_suction_span']:+d}, tip ring centre y {g['tip_ring_centre_y_m']:+.2f} m "
             "(prebend upwind = -y)")
    L.append(f"- chord psi vs BlTwist: max |psi - BlTwist| {g['psi_minus_bltwist_deg_max_abs']:.2f} deg, "
             f"max |psi + BlTwist| {g['psi_plus_bltwist_deg_max_abs']:.2f} deg")
    L.append(f"- alpha(mesh chord, BEM inflow) - alpha_BEM: max {a['max_abs_alpha_mesh_minus_bem_deg']:.2f} deg; "
             f"with the chord at -BlTwist: {a['max_abs_alpha_minus_twist_minus_bem_deg']:.3f} deg; "
             f"production _twist_mesh_to_bem {a['twist_mesh_to_bem']:+.0f}")
    L.append("")
    L.append("| z [m] | psi mesh [deg] | BlTwist [deg] | alpha BEM | alpha mesh chord | alpha chord at -tw | d alpha / +1 deg about +z |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|")
    for q in a["rows"][::4]:
        L.append(f"| {q['z']:.1f} | {q['psi_deg']:+.2f} | {q['bl_twist_deg']:+.2f} | {q['alpha_bem_deg']:+.2f} | "
                 f"{q['alpha_mesh_chord_deg']:+.2f} | {q['alpha_if_chord_at_minus_twist_deg']:+.2f} | "
                 f"{q['d_alpha_per_plus_1deg_about_z']:+.3f} |")
    L.append("")
    for tag, d in out["loads"].items():
        L.append(f"- loads {tag}: F = ({d['F_total_N'][0]:+.0f}, {d['F_total_N'][1]:+.0f}, {d['F_total_N'][2]:+.0f}) N, "
                 f"torque about +y {d['torque_about_plus_y_Nm']:+.4e} N m, spin +y -> power/blade {d['power_per_blade_W']:+.3e} W")
    ad = out["ad_frame"]
    for k in ("physical_minus_theta", "apply_ad_loads_plus_theta"):
        L.append(f"- AD Fn/Ft {k}: Fx {ad[k]['Fx_N']:+.0f} N, Fy {ad[k]['Fy_N']:+.0f} N, torque {ad[k]['torque_Nm']:.4e} N m, "
                 f"3-blade power {ad[k]['power_3_blades_MW']:.2f} MW")
    v = out["v1_v2"]
    L.append(f"- v1 NuMAD (rotorspin {v['v1_numad_rotorspin']}): psi = -twist {v['v1_psi_equals_minus_twist']}; "
             f"v2 gmsh/windIO frames: psi = +twist {v['v2_gmsh_psi_equals_plus_twist']}")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--work", type=Path, default=Path("/tmp") / "l5-t10-work")
    ap.add_argument("--view", choices=("chain", "flipped"), default="chain")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.work.mkdir(parents=True, exist_ok=True)
    if args.view == "flipped":
        view_flipped(args.out, args.work)
        return
    t0 = time.time()
    mesh, props, asm = _work_mesh(args.work)
    out = {"mesh_nodes": mesh.node_count}
    part, f_deck = _participant_and_loads(mesh, props, asm, -6.0, -4.0, args.work)
    rows = view_geometry(mesh, part, out)
    view_loads(mesh, f_deck, "bem_deck_tilt-6_precone-4", out)
    part0, f0 = _participant_and_loads(mesh, props, asm, 0.0, 0.0, args.work)
    view_loads(mesh, f0, "bem_tilt0_precone0", out)
    view_aoa(part0, rows, out)
    view_ad_frame(out)
    view_v1_v2(rows, out)
    out["elapsed_s"] = time.time() - t0
    (args.out / "chain.json").write_text(json.dumps(out, indent=1))
    (args.out / "chain.md").write_text(markdown(out))
    print(markdown(out))


if __name__ == "__main__":
    main()
