"""#30 T3b part 2: do the two applied fields twist the *structure* differently, statically?

Part 1 (`tools/diagnose_wall_flow_pattern_ab.py`) measured, from the A/B runs' own saved fields,
that the wall-flow realisation delivers the *same-sign, +9 %* twist drive while the coupled tip
rotation moves by a factor ~2.6.  So the coupled movement is a **response** effect, not a load
effect - but two very different causes are still open:

* the structure answers the wall-flow **pattern** with more twist (a static, structural effect:
  section distortion / drilling), or
* the amplification lives in the **coupled loop** (aero-elastic torsional feedback, or the
  transient), and the static response is comparable for both fields.

This script settles it the cheapest way: apply each run's own frozen ``F_AERO`` field to the
production mesh in a **linear static** solve - no dynamics, no aero feedback, no corotational
update - and compare the section twist.  The mesh is regenerated with the case's own generator
parameters and is index-wise identical to the runs' ``solid_mesh.vtu``, so the saved field maps
node by node with no interpolation.

Reading the result:
* ``omega``        - the affine-fit section rotation, ``t._ring_kinematics`` (the repo's estimator);
* ``rigid``        - the rigid-only fit rotation, same helper, no strain allowed;
* ``distortion``   - the affine symmetric-strain magnitude, same helper.
The two rotation estimators are printed side by side because part 1 showed the tip reading is
estimator-sensitive.

Diagnostic only: no assertions, not collected by the suite, writes nothing.  Run with the
repository bootstrap::

    scripts/aeroenv.sh python tools/diagnose_wall_flow_static_ab.py
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import meshio  # noqa: E402
from _aeroelast import PyMeshAssembler  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

import tests.validation.blade.test_blade_iea15mw_validation as bv  # noqa: E402
import tests.validation.blade.test_blade_rated_twist as t  # noqa: E402
from aeroelast.core.mesh.generators import BladeMesh  # noqa: E402
from aeroelast.models.blade.model import build_rust_properties  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

RUNS = {"A (min-norm)": "A", "B (wall flow)": "B"}
YAML = DATA_DIR / "IEA-15-240-RWT.yaml"
AIRFOILS = DATA_DIR / "airfoils"
SPAN = np.array([0.0, 0.0, 1.0])
# The case's own generator parameters (case5s/solid_corotational_yaw_0.yaml).
GENERATOR = {
    "element_size": 0.25,
    "n_samples": 300,
    "airfoil_spacing": "constant",
    "span_grading": "chord",
}
N_STATIONS = 14


def load_field(root: Path, step: str) -> np.ndarray:
    """The run's applied nodal force field, shape ``(n_nodes, 3)``."""
    path = root / "corotational" / step / "fields.vtu"
    if not path.is_file():
        raise SystemExit(f"missing {path}")
    return np.asarray(meshio.read(path).point_data["F_AERO"], dtype=float)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="#30 T3b part 2: a linear static solve under each run's frozen field"
    )
    parser.add_argument(
        "--base",
        type=Path,
        default=Path("/scratch/leahk/eduardo.donestevez/bfs16"),
        help="directory holding the A/ and B/ run trees",
    )
    parser.add_argument("--step", default="5", help="time-step directory under corotational/")
    parser.add_argument("--stations", type=int, default=N_STATIONS)
    args = parser.parse_args()
    np.set_printoptions(precision=6, suppress=True, linewidth=160)

    # ── mesh: the case's own generator, checked against the runs' own solid_mesh.vtu ─────────
    t0 = time.time()
    generator = BladeMesh(yaml_file=str(YAML), airfoil_dir=str(AIRFOILS), **GENERATOR)
    mesh = generator.generate(renumber="rcm", verbose=False)
    props = build_rust_properties(generator.numad_mesh_data)
    coords = mesh.coords_array
    for key in RUNS:
        solid = args.base / RUNS[key] / "solid_mesh.vtu"
        if solid.is_file():
            drift = float(np.abs(np.asarray(meshio.read(solid).points) - coords).max())
            if drift > 1e-9:
                raise SystemExit(f"{solid}: mesh drift {drift:.3e} m - do not map node by node")
    print(
        f"mesh: {coords.shape[0]} nodes, {len(mesh.elements)} elements, "
        f"element_size={GENERATOR['element_size']} ({time.time() - t0:.1f} s); "
        f"index-wise identical to both runs' solid_mesh.vtu"
    )

    # ── assembly + clamped root (the case's own dirichlet: RootNodes, all DOFs) ─────────────
    t0 = time.time()
    k_ff, n_dof, free = _assemble(mesh, props)
    print(f"stiffness: {n_dof} dof, {free.size} free, assembly {time.time() - t0:.1f} s")

    fields = {key: load_field(args.base / RUNS[key], args.step) for key in RUNS}
    if any(f.shape[0] != coords.shape[0] for f in fields.values()):
        raise SystemExit("the saved field and the mesh disagree on the node count")

    stations = t._physical_stations(coords)
    tip = np.where(np.abs(coords[:, 2] - stations[-1]) < t.STATION_GAP_TOLERANCE)[0]
    # cumulative applied twist moment about the reference span axis, per station
    moment = {
        key: np.asarray(
            [
                _moment(coords[r] - coords[r].mean(axis=0), fields[key][r]) @ SPAN
                for r in [
                    np.where(np.abs(coords[:, 2] - z) < t.STATION_GAP_TOLERANCE)[0]
                    for z in stations
                ]
            ]
        )
        for key in RUNS
    }

    print(f"\nstations: {len(stations)}, tip ring {tip.size} nodes at z={stations[-1]:.3f} m")
    for key in RUNS:
        f = fields[key]
        print(
            f"  {RUNS[key]:16s} L2={np.linalg.norm(f):.6e} N  max|f|={np.abs(f).max():.6e} N  "
            f"sum F={np.array2string(f.sum(axis=0), precision=3)}"
        )

    solutions = {}
    reports = {}
    for key in RUNS:
        t0 = time.time()
        f = np.zeros(n_dof)
        f[0::6] = fields[key][:, 0]
        f[1::6] = fields[key][:, 1]
        f[2::6] = fields[key][:, 2]
        u = np.zeros(n_dof)
        u[free] = np.asarray(spsolve(k_ff, f[free]), dtype=float).ravel()
        solutions[key] = u
        reports[key] = _section_report(coords, u, stations)
        print(
            f"  {RUNS[key]:16s} static solve {time.time() - t0:.1f} s   max|u|={np.abs(u.reshape(-1, 6)[:, :3]).max():.4f} m"
        )

    print(f"\n=== section twist profile [deg] over {N_STATIONS} stations ===")
    print(
        f"{'z [m]':>7} | {'omega A':>9} {'omega B':>9} {'ratio':>7} | "
        f"{'rigid A':>9} {'rigid B':>9} | {'dist A':>9} {'dist B':>9} | {'Tcum A':>10} {'Tcum B':>10}"
    )
    pick = np.linspace(0, len(stations) - 1, args.stations).astype(int)
    for i in pick:
        zz = stations[i]
        oa, ob = (
            np.rad2deg(reports[list(RUNS)[0]]["omega"][i]),
            np.rad2deg(reports[list(RUNS)[1]]["omega"][i]),
        )
        ratio = ob / oa if abs(oa) > 1e-12 else float("nan")
        print(
            f"{zz:>7.2f} | {oa:>9.4f} {ob:>9.4f} {ratio:>7.3f} | "
            f"{np.rad2deg(reports[list(RUNS)[0]]['rigid'][i]):>9.4f} "
            f"{np.rad2deg(reports[list(RUNS)[1]]['rigid'][i]):>9.4f} | "
            f"{reports[list(RUNS)[0]]['distortion'][i]:>9.2e} "
            f"{reports[list(RUNS)[1]]['distortion'][i]:>9.2e} | "
            f"{moment[list(RUNS)[0]][i]:>10.3e} {moment[list(RUNS)[1]][i]:>10.3e}"
        )

    first, second = list(RUNS)
    tip_omega = {
        key: float(np.rad2deg(t._ring_kinematics(coords, solutions[key], tip)["omega"]))
        for key in RUNS
    }
    tip_rigid = {
        key: float(np.rad2deg(t._ring_kinematics(coords, solutions[key], tip)["rigid_rotation"]))
        for key in RUNS
    }
    tip_dist = {
        key: float(t._ring_kinematics(coords, solutions[key], tip)["distortion"]) for key in RUNS
    }
    print("\n=== tip section (the coupled runs moved this by a factor 2.6, same sign) ===")
    for key in RUNS:
        print(
            f"  {RUNS[key]:16s} omega {tip_omega[key]:+9.4f} deg   rigid {tip_rigid[key]:+9.4f} deg   "
            f"distortion {tip_dist[key]:.4e}   applied T(at tip ring) "
            f"{moment[key][-1]:.4e} N.m"
        )
    print(
        f"  static ratio B/A:  omega {tip_omega[second] / tip_omega[first]:+.4f}   "
        f"rigid {tip_rigid[second] / tip_rigid[first]:+.4f}   "
        f"distortion {tip_dist[second] / tip_dist[first]:+.4f}   "
        f"applied T {moment[second][-1] / moment[first][-1]:+.4f}"
    )

    # linear superposition: the response to the pattern difference alone
    u_diff = solutions[second] - solutions[first]
    diff_tip = t._ring_kinematics(coords, u_diff, tip)
    print(
        f"\n  response to the difference field alone (linear superposition): omega "
        f"{np.rad2deg(diff_tip['omega']):+.4f} deg, rigid "
        f"{np.rad2deg(diff_tip['rigid_rotation']):+.4f} deg, distortion {diff_tip['distortion']:.4e}"
    )

    # The very last ring is 16 nodes with distortion ~1.0 (fully distorted), so it cannot decide
    # the ratio.  The outer-span window does, and it is where the coupled movement lives.
    z_lo = float(coords[:, 2].min())
    z_hi = float(coords[:, 2].max())
    window = (z_lo + 0.80 * (z_hi - z_lo), z_lo + 0.96 * (z_hi - z_lo))
    pick_w = np.asarray([i for i, zz in enumerate(stations) if window[0] <= zz <= window[1]])
    omega_ratio = reports[second]["omega"][pick_w] / reports[first]["omega"][pick_w]
    rigid_ratio = reports[second]["rigid"][pick_w] / reports[first]["rigid"][pick_w]
    dist_ratio = reports[second]["distortion"][pick_w] / reports[first]["distortion"][pick_w]
    print(
        f"\n=== outer-span window z in [{window[0]:.1f}, {window[1]:.1f}] m ({pick_w.size} stations; "
        f"the coupled runs moved this by ~2.6x with the SAME sign) ==="
    )
    print(
        f"  static B/A: omega  median {np.median(omega_ratio):+.3f} "
        f"range [{omega_ratio.min():+.3f}, {omega_ratio.max():+.3f}]  |  "
        f"rigid median {np.median(rigid_ratio):+.3f} range [{rigid_ratio.min():+.3f}, "
        f"{rigid_ratio.max():+.3f}]  |  distortion median {np.median(dist_ratio):.3f}"
    )
    print(
        f"  static tip deflection: max|u| A {np.abs(solutions[first].reshape(-1, 6)[:, :3]).max():.4f} m "
        f"vs B {np.abs(solutions[second].reshape(-1, 6)[:, :3]).max():.4f} m "
        f"(the coupled runs read 24.10 / 23.63 m)"
    )
    return 0


def _assemble(mesh, props) -> tuple:
    """Assemble the material stiffness and the clamped-root free-DOF set (the case's own BC)."""
    assembler = PyMeshAssembler.from_model(
        bv._to_rust_mesh(mesh, props), props, list(bv.SPAN_DIRECTION), None
    )
    n_dof = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix(
        (np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n_dof, n_dof)
    ).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    fixed = {6 * i + d for i in root for d in range(6)}
    free = np.asarray([i for i in range(n_dof) if i not in fixed], dtype=np.int64)
    return K[np.ix_(free, free)].tocsc(), n_dof, free


def _moment(r: np.ndarray, f: np.ndarray) -> np.ndarray:
    return np.cross(r, f).sum(axis=0)


def _section_report(coords: np.ndarray, u: np.ndarray, stations: list[float]) -> dict:
    """Per-station affine rotation, rigid rotation and distortion from the displacement ``u``."""
    omega, rigid, distortion = [], [], []
    for zz in stations:
        ring = np.where(np.abs(coords[:, 2] - zz) < t.STATION_GAP_TOLERANCE)[0]
        kin = t._ring_kinematics(coords, u, ring)
        omega.append(kin["omega"])
        rigid.append(kin["rigid_rotation"])
        distortion.append(kin["distortion"])
    return {
        "omega": np.asarray(omega),
        "rigid": np.asarray(rigid),
        "distortion": np.asarray(distortion),
    }


if __name__ == "__main__":
    sys.exit(main())
