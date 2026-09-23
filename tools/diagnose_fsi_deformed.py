"""Diagnose the BEM-FSI deformed-geometry pipeline: yaml blade vs xlsx blade.

Replicates the FSI participant's per-window reconstruction (strip node
assignment, chord-direction PCA, deformed r/twist, BEM rebuild) for both
blade inputs and applies a realistic first-window flapwise displacement
field, comparing the two paths station by station.

Usage (SDumont):
  module load glu gcc/14.2.0_sequana
  python tools/diagnose_fsi_deformed.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))


def chord_dirs(strip_indices, coords, span_dir, normal_dir, tangential_dir):
    out = []
    for idx in strip_indices:
        if len(idx) < 2:
            out.append(normal_dir.copy())
            continue
        pts = coords[idx]
        offsets = pts - pts.mean(axis=0)
        offsets_plane = offsets - np.outer(offsets @ span_dir, span_dir)
        _, _, Vt = np.linalg.svd(offsets_plane, full_matrices=False)
        c = Vt[0]
        if abs(np.dot(c, tangential_dir)) >= abs(np.dot(c, normal_dir)):
            if np.dot(c, tangential_dir) < 0:
                c = -c
        else:
            if np.dot(c, normal_dir) < 0:
                c = -c
        n = np.linalg.norm(c)
        out.append(c / n if n > 1e-12 else normal_dir.copy())
    return out


def run_pipeline(name, aero, mesh, span_dir, normal_dir, tangential_dir):
    from aeroelast.models.blade.aerodynamics import AeroStation, BladeAero
    from aeroelast.solvers.bem.engine import BEMSolver
    from aeroelast.solvers.bem.force_projection import ForceProjector

    proj = ForceProjector(mesh, aero, span_direction=span_dir,
                          normal_direction=normal_dir,
                          tangential_direction=tangential_dir)
    strips = proj._strips
    idx = [s.node_indices for s in strips]
    coords = mesh.coords_array

    ref_chord = chord_dirs(idx, coords, span_dir, normal_dir, tangential_dir)

    # realistic first-window deformation: rigid rotor rotation about the
    # Y axis (rotation_axis in the solid yaml) + flapwise bending
    z = coords[:, 2]
    L = z.max()
    theta = 0.055  # rad — blade angle at t≈0.07 s
    c, s_ = np.cos(theta), np.sin(theta)
    R = np.array([[c, 0.0, s_], [0.0, 1.0, 0.0], [-s_, 0.0, c]])
    disp = (R @ coords.T).T - coords
    disp[:, 1] += 0.5 * (z / L) ** 2
    def_coords = coords + disp

    def_chord = chord_dirs(idx, def_coords, span_dir, normal_dir, tangential_dir)
    for k in range(len(idx)):
        if np.dot(def_chord[k], ref_chord[k]) < 0:
            def_chord[k] = -def_chord[k]

    s = span_dir
    r_def = np.empty(len(idx))
    twist_def = np.empty(len(idx))
    for k in range(len(idx)):
        if len(idx[k]) == 0:
            r_def[k] = aero.r[k]
            twist_def[k] = aero.twist[k]
            continue
        r_def[k] = float(np.mean(def_coords[idx[k]] @ s)) + aero.hub_radius
        c_ref, c_def = ref_chord[k], def_chord[k]
        cos_a = float(np.clip(np.dot(c_ref, c_def), -1.0, 1.0))
        sin_a = float(np.dot(np.cross(c_ref, c_def), s))
        twist_def[k] = aero.twist[k] + np.arctan2(sin_a, cos_a)

    stations = [AeroStation(span_fraction=st.span_fraction, r=float(r_def[k]),
                            chord=st.chord, twist=float(twist_def[k]),
                            pitch_axis=st.pitch_axis, airfoil=st.airfoil)
                for k, st in enumerate(aero.stations)]
    ba = BladeAero(stations=stations, airfoils=aero.airfoils,
                   hub_radius=aero.hub_radius, rotor_radius=aero.rotor_radius,
                   n_blades=aero.n_blades, blade_length=aero.blade_length)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        bem = BEMSolver(ba, rho=1.225, mu=1.81e-05, precone=0.0, tilt=0.0,
                        yaw=0.0, hub_height=150.0, shear_exp=0.0)
        res = bem.compute(v_inf=10.59, omega=0.7906341464750989, pitch=0.0)
    return {
        "r_def": r_def, "twist_def": twist_def, "r_ref": aero.r,
        "twist_ref": aero.twist, "Np": res.Np, "alpha": res.alpha,
        "warnings": len(w),
    }


def main():
    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.aerodynamics import load_blade_aero

    s = np.array([0.0, 0.0, 1.0])
    n = np.array([0.0, 1.0, 0.0])
    t = np.array([-1.0, 0.0, 0.0])

    aero_y = load_blade_aero(str(REPO / "tests/IEA-15-240-RWT.yaml"),
                             default_re=1e7, neuralfoil_model="large",
                             viterna_ar=17.0, viterna_confidence_threshold=0.5)
    mesh_y = BladeMesh(yaml_file=str(REPO / "tests/IEA-15-240-RWT.yaml"),
                       element_size=0.25).generate(renumber="rcm")
    r_y = run_pipeline("yaml", aero_y, mesh_y, s, n, t)

    aero_x = load_blade_aero(str(REPO / "tests/NuMAD_utd_iea15mw.xlsx"),
                             default_re=1e7, neuralfoil_model="large",
                             airfoil_dir=str(REPO / "tests/airfoils"),
                             viterna_ar=17.0, viterna_confidence_threshold=0.5)
    mesh_x = BladeMesh(excel_file=str(REPO / "tests/NuMAD_utd_iea15mw.xlsx"),
                       airfoil_dir=str(REPO / "tests/airfoils"),
                       element_size=0.25).generate(renumber="rcm")
    r_x = run_pipeline("xlsx", aero_x, mesh_x, s, n, t)

    print(f"\n{'st':>3s} | {'r_def y/x':>16s} | {'dtwist y/x [deg]':>18s} | "
          f"{'Np y/x [kN/m]':>16s} | {'alpha y/x [deg]':>18s}")
    print("-" * 82)
    for k in range(len(r_y["r_ref"])):
        dty = np.rad2deg(r_y["twist_def"][k] - r_y["twist_ref"][k])
        dtx = np.rad2deg(r_x["twist_def"][k] - r_x["twist_ref"][k])
        flag = " <== " if abs(dty) > 30 else ""
        print(f"{k:3d} | {r_y['r_def'][k]-r_y['r_ref'][k]:8.3f} "
              f"{r_x['r_def'][k]-r_x['r_ref'][k]:7.3f} | {dty:9.2f} {dtx:8.2f} | "
              f"{r_y['Np'][k]/1e3:8.2f} {r_x['Np'][k]/1e3:7.2f} | "
              f"{r_y['alpha'][k]:9.2f} {r_x['alpha'][k]:8.2f}{flag}")
    print(f"\nyaml warnings: {r_y['warnings']}  | xlsx warnings: {r_x['warnings']}")


if __name__ == "__main__":
    main()
