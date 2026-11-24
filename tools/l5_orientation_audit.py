"""L5 T12: measured orientation audit of the shell blade, the windIO yaml and the BeamDyn deck (odd/tasks/l5-one-way.md L16).

Every orientation is MEASURED from coordinates, never taken from a comment or a convention statement:

  A. physical frame (no convention): rotor axis and spin from the case solid.yaml, blade velocity Omega x r, wind direction
     from the sign of the total aerodynamic thrust of the production BEM loads on the mesh
  B. mesh, per ring: chord length, LE side, suction (cambered) side, chord angle psi (TE->LE from +x), reference-axis point
     (LE + section_offset_y along the chord) -> prebend (mesh y) and sweep (mesh x)
  C. yaml (windIO v2): chord, twist, section_offset_y, reference axis x/y/z at the same span
  D. BeamDyn deck (IEA-15-240-RWT_BeamDyn.dat): kp_xr/kp_yr/kp_zr, initial_twist, mapped into the mesh frame with the
     parity tool's MESH_TO_BD (tools/twist_p1.py:37); BeamDyn places a section at -initial_twist about +z
  E. physics-only angle of attack (no BEM): relative wind (-Omega r (1 + a'), V (1 - a)) against the measured chord, for
     the measured mesh chord and for the mirrored one (psi -> -psi), a in {0, 1/3}
  F. response direction of the T1/T3 parity runs: angle of the (edge, flap) deflection vector, shell vs BeamDyn, i.e. the
     rotation that would explain the edge difference as flap leakage

  scripts/aeroenv.sh python tools/l5_orientation_audit.py --out logs/check/l5_oneway/t12

Mesh frame: z span, rotor axis +y (case solid.yaml rotation_axis), Omega > 0.
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

YAML = REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
SOLID_YAML = REPO / "cases/l5_oneway_rated/solid/solid.yaml"
BD_PRIMARY = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat")
PARITY_DIR = REPO / "logs/check/l5_oneway/t11"
WIND = 10.59
RPM = 7.56
HUB_RADIUS = 3.97


# --- pure helpers (tests/software/test_l5_orientation_audit.py) --------------------------------------------------


def parse_bd_keypoints(text: str) -> np.ndarray:
    """(n, 4) kp_xr, kp_yr, kp_zr, initial_twist [deg] from a BeamDyn primary file (rows after the 'kp_xr' header)."""
    lines = text.splitlines()
    i = next(k for k, ln in enumerate(lines) if "kp_xr" in ln and "kp_yr" in ln)
    n = int(next(ln.split()[0] for ln in lines if "kp_total" in ln))
    rows = [list(map(float, lines[i + 2 + k].split()[:4])) for k in range(n)]
    return np.array(rows)


def bd_to_mesh(v_bd, mesh_to_bd):
    """BeamDyn blade-frame vector(s) -> mesh frame (v_mesh = v_bd @ MESH_TO_BD, as the parity tool)."""
    return np.asarray(v_bd, float) @ np.asarray(mesh_to_bd, float)


def bd_chord_angle_deg(initial_twist_deg, mesh_to_bd):
    """Mesh angle of the TE->LE direction of a BeamDyn section: the zero-twist TE->LE direction is -y_BD, rotated by
    -initial_twist about +z (BeamDyn_IO.f90:898), expressed in the mesh frame."""
    t = -np.deg2rad(np.asarray(initial_twist_deg, float))
    d_bd = np.stack([np.sin(t), -np.cos(t), np.zeros_like(t)], axis=-1)  # Rz(t) @ (0, -1, 0)
    d = bd_to_mesh(d_bd, mesh_to_bd)
    return np.degrees(np.arctan2(d[..., 1], d[..., 0]))


def physical_aoa_deg(psi_deg, suction_normal, r, wind_dir_y, omega, wind=WIND, a=0.0, ap=0.0):
    """Angle of attack from kinematics only: air relative to the section = wind (1 - a) along +/-y minus blade velocity
    (+x for Omega > 0 about +y at +z).  Positive when the flow comes from the pressure side (v . n_suction > 0), the
    convention of tools/l5_twist_sense_chain.angle_of_attack_deg."""
    v = np.array([-omega * r * (1.0 + ap), wind_dir_y * wind * (1.0 - a)])
    p = np.deg2rad(psi_deg)
    le_to_te = -np.array([np.cos(p), np.sin(p)])
    return float(np.degrees(np.arctan2(np.dot(v, suction_normal), np.dot(v, le_to_te))))


def rotate2(v, deg):
    """Right-handed rotation of a 2-D vector about +z (a mirrored pretwist is the section rotated by -2 psi, not reflected)."""
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])


def response_angle_deg(edge, flap):
    """Direction of the in-plane deflection vector (edge = mesh x, flap = mesh y) from +y, degrees."""
    return np.degrees(np.arctan2(np.asarray(edge, float), np.asarray(flap, float)))


# --- views --------------------------------------------------------------------------------------------------------


def yaml_blade():
    import yaml

    b = yaml.safe_load(YAML.read_text())["components"]["blade"]
    ra, os_ = b["reference_axis"], b["outer_shape"]
    return {"eta_ax": np.array(ra["z"]["grid"]), "x": np.array(ra["x"]["values"]), "y": np.array(ra["y"]["values"]),
            "z": np.array(ra["z"]["values"]), "x_grid": np.array(ra["x"]["grid"]), "y_grid": np.array(ra["y"]["grid"]),
            "chord": (np.array(os_["chord"]["grid"]), np.array(os_["chord"]["values"])),
            "twist": (np.array(os_["twist"]["grid"]), np.array(os_["twist"]["values"])),
            "offset": (np.array(os_["section_offset_y"]["grid"]), np.array(os_["section_offset_y"]["values"]))}


def yaml_at(yb, z):
    eta = np.interp(z, yb["z"], yb["eta_ax"])
    f = {k: float(np.interp(eta, *yb[k])) for k in ("chord", "twist", "offset")}
    f["ref_x"] = float(np.interp(eta, yb["x_grid"], yb["x"]))
    f["ref_y"] = float(np.interp(eta, yb["y_grid"], yb["y"]))
    f["eta"] = float(eta)
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    import yaml
    import l5_oneway_bd_parity as P
    import twist_p1 as tp
    from l5_twist_sense_chain import section_from_ring

    out = {}
    # A. physical frame
    rotor = yaml.safe_load(SOLID_YAML.read_text())["solver"]["rotor"]
    axis = np.array(rotor["rotation_axis"], float)
    omega = float(rotor["omega"])
    v_tip = np.cross(omega * axis, [0.0, 0.0, 1.0])
    mesh, props, asm = P.build_shell(20)
    f, part = P.shell_loads("bem", asm, mesh, props, -6.0, -4.0)
    f3 = np.asarray(f).reshape(-1, asm.dofs_per_node)[:, :3]
    thrust = f3.sum(0)
    wind_dir_y = float(np.sign(thrust[1]))
    out["frame"] = {"rotation_axis": axis.tolist(), "omega": omega, "blade_velocity_dir_at_+z": v_tip.tolist(),
                    "bem_total_force_N": thrust.tolist(), "wind_dir_y_from_thrust": wind_dir_y,
                    "upwind_is": "-y" if wind_dir_y > 0 else "+y"}

    # B-C. mesh rings vs yaml
    X = np.array([n.coords for n in mesh.nodes], float)
    yb = yaml_blade()
    rows = []
    for ring in part._strip_ring_indices:
        ring = np.asarray(ring)
        if len(ring) < 8:
            continue
        pts = X[ring]
        z = float(pts[:, 2].mean())
        sec = section_from_ring(pts[:, :2])
        y = yaml_at(yb, z)
        u_le_te = (sec["te"] - sec["le"]) / sec["chord"]
        ref = sec["le"] + y["offset"] * u_le_te
        aoa = {}
        for tag, psi, n_s in (("measured", sec["psi_deg"], sec["suction_normal"]),
                              ("mirrored", -sec["psi_deg"], rotate2(sec["suction_normal"], -2.0 * sec["psi_deg"]))):
            for a in (0.0, 1.0 / 3.0):
                aoa[f"{tag}_a{a:.2f}"] = physical_aoa_deg(psi, n_s, z + HUB_RADIUS, wind_dir_y, omega, a=a)
        rows.append({"z": z, "eta": y["eta"], "chord_mesh": sec["chord"], "chord_yaml": y["chord"],
                     "psi_mesh_deg": sec["psi_deg"], "twist_yaml_deg": y["twist"],
                     "le_dot_v": float(np.dot(sec["le"] - sec["te"], v_tip[:2])),
                     "suction_ny": float(sec["suction_normal"][1]), "camber_over_c": sec["camber_over_c"],
                     "end_width_ratio": sec["end_width_ratio"],
                     "ref_mesh_x": float(ref[0]), "ref_mesh_y": float(ref[1]),
                     "yaml_ref_x_prebend": y["ref_x"], "yaml_ref_y_sweep": y["ref_y"], "aoa_deg": aoa})
    good = [r for r in rows if r["camber_over_c"] > 0.004 and r["end_width_ratio"] > 2.0 and 10.0 < r["z"] < 115.0]

    # D. BeamDyn deck in the mesh frame
    kp = parse_bd_keypoints(BD_PRIMARY.read_text())
    kp_mesh = bd_to_mesh(kp[:, :3], tp.MESH_TO_BD)
    psi_bd = bd_chord_angle_deg(kp[:, 3], tp.MESH_TO_BD)
    for r in rows:
        r["bd_ref_mesh_x"] = float(np.interp(r["z"], kp_mesh[:, 2], kp_mesh[:, 0]))
        r["bd_ref_mesh_y"] = float(np.interp(r["z"], kp_mesh[:, 2], kp_mesh[:, 1]))
        r["bd_twist_deg"] = float(np.interp(r["z"], kp[:, 2], kp[:, 3]))
        r["psi_bd_deg"] = float(np.interp(r["z"], kp_mesh[:, 2], psi_bd))

    def mx(key):
        return float(max(abs(r[key]) for r in good))

    out["summary"] = {
        "n_rings": len(rows), "n_good": len(good),
        "LE_moves_with_blade_velocity_all": all(r["le_dot_v"] > 0 for r in good),
        "suction_faces_downwind_all": all(r["suction_ny"] * wind_dir_y > 0.9 for r in good),
        "max_abs_psi_plus_twist_deg": float(max(abs(r["psi_mesh_deg"] + r["twist_yaml_deg"]) for r in good)),
        "max_abs_psi_minus_twist_deg": float(max(abs(r["psi_mesh_deg"] - r["twist_yaml_deg"]) for r in good)),
        "max_abs_psi_mesh_minus_psi_bd_deg": float(max(abs(r["psi_mesh_deg"] - r["psi_bd_deg"]) for r in good)),
        "max_abs_bd_twist_minus_yaml_deg": float(max(abs(r["bd_twist_deg"] - r["twist_yaml_deg"]) for r in good)),
        "max_rel_chord_err": float(max(abs(r["chord_mesh"] / r["chord_yaml"] - 1) for r in good)),
        "tip_prebend_mesh_y_m": rows[-1]["ref_mesh_y"], "tip_prebend_yaml_x_m": rows[-1]["yaml_ref_x_prebend"],
        "tip_prebend_bd_mesh_y_m": rows[-1]["bd_ref_mesh_y"],
        "max_abs_prebend_mesh_minus_yaml_m": mx_diff(good, "ref_mesh_y", "yaml_ref_x_prebend"),
        "max_abs_prebend_mesh_minus_bd_m": mx_diff(good, "ref_mesh_y", "bd_ref_mesh_y"),
        "max_abs_sweep_mesh_x_m": mx("ref_mesh_x"),
        "aoa_measured_a0.33_range_deg": aoa_range(good, "measured_a0.33"),
        "aoa_mirrored_a0.33_range_deg": aoa_range(good, "mirrored_a0.33"),
    }
    out["rings"] = rows

    # F. response direction of the parity runs
    resp = {}
    for name in ("ad", "bem", "ad-parked-y"):
        p = PARITY_DIR / f"{name}.json"
        if not p.exists():
            continue
        t = json.loads(p.read_text())["table"]
        sh, bd = t["shell"], t["beamdyn_x0.01"]
        a_sh, a_bd = response_angle_deg(sh["edge"], sh["flap"]), response_angle_deg(bd["edge"], bd["flap"])
        resp[name] = {"z": t["z"], "angle_shell_deg": a_sh.tolist(), "angle_bd_deg": a_bd.tolist(),
                      "delta_deg": (a_sh - a_bd).tolist(),
                      "flap_ratio": (np.array(sh["flap"]) / np.array(bd["flap"])).tolist()}
    out["response"] = resp
    out["elapsed_s"] = time.time() - t0
    (args.out / "audit.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({"frame": out["frame"], "summary": out["summary"]}, indent=1, default=float))
    for name, r in resp.items():
        print(name, "delta response angle deg:", np.round(r["delta_deg"], 2).tolist())


def mx_diff(rows, a, b):
    return float(max(abs(r[a] - r[b]) for r in rows))


def aoa_range(rows, key):
    v = [r["aoa_deg"][key] for r in rows]
    return [float(min(v)), float(max(v))]


if __name__ == "__main__":
    main()
