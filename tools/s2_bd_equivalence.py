"""S-2 load-wrench equivalence: shell load sets vs the BeamDyn-driver load sets.

The BeamDyn static reference (tools/bd_static_reference.py) and the shell S-2 cases
(tools/run_s2_static_cases.py, used unchanged here) do not apply the same loads by
construction: BD q is per reference-axis arc length, shell LC4 is per projected z lumped
into 60 bins; BD gravity acts on the 6x6 mass, shell gravity on the shell mass; BD tip load
and output sit on the reference axis, shell tip load is spread over the max-z nodes and the
S-2 output is the first max-z node. This tool measures each difference instead of assuming it:

  * resultant force and moment about the root origin of the load the shell solver actually
    assembles (captured from StaticLinearSolver._build_f_ext) and of the BD load set
    integrated on the published reference axis (shell frame: (-yBD, xBD, zBD));
  * shell tip uy at the S-2 output node (first max-z node) and as the mean over the tip ring,
    next to the BD p60 linear tip uy (logs/check/bd_static/report.md).

Usage (skin reference forced through tools/with_skin_reference.py, one SLURM task):
    scripts/aeroenv.sh python tools/with_skin_reference.py oml tools/s2_bd_equivalence.py 0.25 out.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))

BD_GEOMETRY = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat")
BD_PROPERTIES = REPO_ROOT / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
# BeamDyn driver p60, scale 1e-4, linear tip uy in shell axes (m); p40-p60 within 0.13 %.
BD_TIP_UY = {"LC1_gravity": -2.2315, "LC2_tip_flap": 82.273, "LC4_distributed": 9.0626}
GRAVITY = np.array([0.0, -9.81, 0.0])
TIP_FORCE = np.array([0.0, 1.0e6, 0.0])
Q = np.array([0.0, 5.0e3, 0.0])


def bd_axis(n: int = 20001):
    """Reference axis in shell axes, sampled uniformly in arc length: (s/S, points, S)."""
    rows, lines = [], BD_GEOMETRY.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "kp_xr" in line) + 2
    for line in lines[start:start + 50]:
        rows.append([float(v) for v in line.split()[:3]])
    kp = np.asarray(rows)
    kp = np.column_stack([-kp[:, 1], kp[:, 0], kp[:, 2]])
    seg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(kp, axis=0), axis=1))]
    s = np.linspace(0.0, seg[-1], n)
    pts = np.column_stack([np.interp(s, seg, kp[:, k]) for k in range(3)])
    return s / seg[-1], pts, seg[-1]


def wrench(points: np.ndarray, forces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return forces.sum(axis=0), np.cross(points, forces).sum(axis=0)


def bd_wrenches() -> dict:
    from beam_reference import load_beamdyn_blade

    eta, pts, arc = bd_axis()
    w = np.full(eta.size, arc / (eta.size - 1))
    w[[0, -1]] *= 0.5  # trapezoid weights on the uniform arc-length grid
    bd = load_beamdyn_blade(BD_PROPERTIES, blade_length=1.0)
    m = np.interp(eta, bd.z, bd.mass_dens)
    out = {"arc_length_m": arc, "tip_point_m": pts[-1].tolist(), "mass_kg": float((m * w).sum())}
    for name, f in (("LC1_gravity", np.outer(m * w, GRAVITY)),
                    ("LC4_distributed", np.outer(w, Q))):
        F, M = wrench(pts, f)
        out[name] = {"F": F.tolist(), "M_root": M.tolist()}
    F, M = wrench(pts[-1:], TIP_FORCE[None, :])
    out["LC2_tip_flap"] = {"F": F.tolist(), "M_root": M.tolist()}
    return out


def shell_cases(es: float) -> dict:
    from run_s2_static_cases import run_shell_case

    from aeroelast.core.mesh.generators import BladeMesh
    from aeroelast.models.blade.model import build_rust_properties
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    captured = {}
    original = StaticLinearSolver._build_f_ext

    def capture(self):
        f = original(self)
        captured["f"] = f
        return f

    StaticLinearSolver._build_f_ext = capture
    gen = BladeMesh(yaml_file=str(REPO_ROOT / "tests/IEA-15-240-RWT.yaml"), element_size=es)
    mesh = gen.generate(renumber="rcm")
    props = build_rust_properties(gen.numad_mesh_data)
    xyz = mesh.coords_array
    z = xyz[:, 2]
    ring = np.nonzero(z >= z.max() - 1e-6)[0]
    out = {"skin_reference": gen.skin_reference if hasattr(gen, "skin_reference") else None,
           "n_nodes": int(xyz.shape[0]), "z_range_m": [float(z.min()), float(z.max())],
           "tip_ring_nodes": int(ring.size), "tip_ring_centroid_m": xyz[ring].mean(axis=0).tolist()}
    loads = {"LC1_gravity": {"type": "gravity", "value": tuple(GRAVITY)},
             "LC2_tip_flap": {"type": "tip", "value": tuple(TIP_FORCE)},
             "LC4_distributed": {"type": "distributed", "value": tuple(Q), "n_slices": 60}}
    for name, load in loads.items():
        u = run_shell_case(mesh, props, load).reshape(-1, 6)
        f = captured["f"].reshape(-1, 6)
        F, M = wrench(xyz, f[:, :3])
        M = M + f[:, 3:].sum(axis=0)
        out[name] = {"F": F.tolist(), "M_root": M.tolist(),
                     "tip_uy_s2_node": float(u[int(np.argmax(z)), 1]),
                     "tip_uy_ring_mean": float(u[ring, 1].mean())}
    return out


def main() -> int:
    es, path = float(sys.argv[1]), Path(sys.argv[2])
    shell, bd = shell_cases(es), bd_wrenches()
    path.write_text(json.dumps({"element_size": es, "shell": shell, "bd": bd}, indent=2) + "\n")
    print(f"shell z {shell['z_range_m']}, tip ring {shell['tip_ring_nodes']} nodes at "
          f"{np.round(shell['tip_ring_centroid_m'], 3)}; BD arc {bd['arc_length_m']:.4f} m, "
          f"tip {np.round(bd['tip_point_m'], 3)}, 6x6 mass {bd['mass_kg']:.0f} kg")
    print(f"{'case':16s} {'Fy shell':>12s} {'Fy BD':>12s} {'dF %':>7s} {'Mx shell':>13s} "
          f"{'Mx BD':>13s} {'dM %':>7s} {'uy node':>9s} {'uy ring':>9s} {'uy BD':>9s} "
          f"{'node %':>7s} {'ring %':>7s}")
    for name, uy_bd in BD_TIP_UY.items():
        s, b = shell[name], bd[name]
        dF = 100 * (s["F"][1] / b["F"][1] - 1)
        dM = 100 * (s["M_root"][0] / b["M_root"][0] - 1)
        print(f"{name:16s} {s['F'][1]:12.5g} {b['F'][1]:12.5g} {dF:7.2f} {s['M_root'][0]:13.6g} "
              f"{b['M_root'][0]:13.6g} {dM:7.2f} {s['tip_uy_s2_node']:9.4f} "
              f"{s['tip_uy_ring_mean']:9.4f} {uy_bd:9.4f} "
              f"{100 * (s['tip_uy_s2_node'] / uy_bd - 1):7.2f} "
              f"{100 * (s['tip_uy_ring_mean'] / uy_bd - 1):7.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
