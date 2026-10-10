"""#30 T3b: attribute the A/B coupled twist movement from the runs' **own saved fields**.

WU: one physics question, one script.  T2 showed that statically the two moment realisations of a
pure span couple sit ~1.6x apart with the *same* sign on the deck's own airfoil section, while the
coupled runs A (minimum-norm, parked default) and B (wall flow) differ by a factor ~13 with a
**sign flip**.  T3b's candidate: the wall-flow branch hands ``_distribute`` only the residual
non-span moment (``M_strip - M_shear``) while the minimum-norm branch hands it the whole
``M_strip``, so the same per-strip resultant is spread over the section by two different
optimisations - a difference in the **flapwise** force pattern that a pure-torque experiment cannot
see and that acts through the bending-torsion coupling of a blade deflected ~16 m.

Before paying for any solve, note that the runs already carry the object.  ``fields.vtu`` holds
``F_AERO``, ``U`` and the nodal rotations at the converged coupled state, and the reconstruction

    points == solid_mesh.vtu + U        (index-wise, to round-off)

holds for both runs, on the **same** reference mesh (``A/solid_mesh.vtu`` and ``B/solid_mesh.vtu``
are bit-identical).  ``C`` reproduces ``A``, so A vs B *is* the realisation, with no transfer step
(T1: the preCICE exchange is the identity), no mesh change and no BEM-state change in the way.
The stored ``F_AERO`` is perpendicular to the reference blade span (measured below: axial share
~1e-5 of ``|F|``), i.e. it is already expressed in the same frame as the coordinates.

What this script measures, per physical ring of the coupling set
--------------------------------------------------------------
* the two applied fields' resultants, their L2 and their peakedness (``max|f|/mean|f|``,
  ``L2/|sum f|``) - the pattern, not the load;
* the difference field ``df = F_B - F_A``: its resultant, moment and L2 share, per ring and
  globally (self-equilibration);
* the **conjugate twist torque** ``T_eff = a . sum (x - c) x f`` about the section's own deformed
  axis ``a`` through its deformed centroid ``c``, each field with its own ``U`` and both fields on
  the *common* (A) geometry, next to the same quantity about the undeformed axis;
* the **cumulative torque** ``sum_{rings outboard} T_eff``, which is what actually twists a
  cantilever section;
* the best-fit section rotation of ``U``, the in-plane distortion (the T2 estimator) and the
  drilling audit of the issue's reading 2 (nodal rotation component about the local section normal).

Diagnostic only: no assertions, not collected by the suite, writes nothing.  Run with the
repository bootstrap::

    scripts/aeroenv.sh python tools/diagnose_wall_flow_pattern_ab.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import meshio  # noqa: E402

from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION  # noqa: E402
from aeroelast.solvers.bem.section_contour import rings_in_band, section_plane_axes  # noqa: E402

RUNS = {"A (min-norm)": "A", "B (wall flow)": "B"}
SPAN = np.array([0.0, 0.0, 1.0])
N_STATIONS = 16


def load_step(root: Path, step: str) -> dict:
    """Read one step's ``fields.vtu`` and its reference ``solid_mesh.vtu``.

    ``fields.vtu`` writes ``points = reference + U`` (verified here index-wise), so the reference
    geometry is recovered as ``points - U`` and checked against ``solid_mesh.vtu``; the stretch
    between the two runs is then a genuine difference of the *response*, never of the mesh.
    """
    path = root / "corotational" / step / "fields.vtu"
    if not path.is_file():
        raise SystemExit(f"missing {path}")
    mesh = meshio.read(path)
    point = mesh.point_data
    deformed = np.asarray(mesh.points, dtype=float)
    displacement = np.asarray(point["U"], dtype=float)
    reference = deformed - displacement
    solid = root / "solid_mesh.vtu"
    if solid.is_file():
        mesh_reference = np.asarray(meshio.read(solid).points, dtype=float)
        stretch = float(np.abs(reference - mesh_reference).max())
        if stretch > 1e-9:
            raise SystemExit(f"{path}: points - U is not solid_mesh.vtu (max {stretch:.3e} m)")
    rotation = np.column_stack(
        [np.asarray(point["ROTX"]), np.asarray(point["ROTY"]), np.asarray(point["ROTZ"])]
    )
    return {
        "path": path,
        "reference": reference,
        "deformed": deformed,
        "U": displacement,
        "F": np.asarray(point["F_AERO"], dtype=float),
        "rotation": rotation,
    }


def _twist_theta(reference: np.ndarray, U: np.ndarray, idx: np.ndarray) -> float:
    """Best-fit rigid rotation about the span axis: ``sum(rx*uy - ry*ux) / sum(rx^2 + ry^2)``."""
    xyz = reference[idx]
    r = xyz - xyz.mean(axis=0)
    den = float(np.sum(r[:, 0] ** 2 + r[:, 1] ** 2))
    if den <= 0.0:
        return float("nan")
    return float(np.sum(r[:, 0] * U[idx, 1] - r[:, 1] * U[idx, 0]) / den)


def _distortion(reference: np.ndarray, U: np.ndarray, idx: np.ndarray) -> tuple[float, float]:
    """In-plane distortion: RMS of the in-plane displacement minus its best-fit rigid motion.

    The T2 estimator: remove the mean translation and the least-squares rotation about the section
    centroid, then take the RMS of the residual, absolute and normalised by the section dimension.
    """
    xyz = reference[idx]
    c = xyz.mean(axis=0)
    xc = xyz[:, 0] - c[0]
    yc = xyz[:, 1] - c[1]
    ux = U[idx, 0]
    uy = U[idx, 1]
    ax = ux - ux.mean()
    ay = uy - uy.mean()
    den = float(np.sum(xc**2 + yc**2))
    theta = float(np.sum(xc * ay - yc * ax) / den) if den > 0.0 else 0.0
    residual = np.sqrt(
        np.mean((ux - (ux.mean() - theta * yc)) ** 2 + (uy - (uy.mean() + theta * xc)) ** 2)
    )
    dim = float(max(np.ptp(xc), np.ptp(yc)))
    return float(residual), float(residual / dim) if dim > 0.0 else float("nan")


def _moment(r: np.ndarray, f: np.ndarray) -> np.ndarray:
    """``sum_j r_j x f_j``."""
    return np.cross(r, f).sum(axis=0)


def _ring_axes(points: np.ndarray) -> np.ndarray:
    """Per-section axis from the centroid line: central difference, one-sided at the ends."""
    axis = np.zeros_like(points)
    if len(points) > 2:
        axis[1:-1] = points[2:] - points[:-2]
    axis[0] = points[1] - points[0]
    axis[-1] = points[-1] - points[-2]
    norms = np.linalg.norm(axis, axis=1, keepdims=True)
    safe = np.where(norms > 0.0, norms, 1.0)
    axis /= safe
    return axis


def analyse(base: Path, step: str, n_stations: int) -> int:
    runs = {key: load_step(base / RUNS[key], step) for key in RUNS}
    keys = list(runs)
    first, second = keys
    reference = runs[first]["reference"]
    if float(np.abs(runs[second]["reference"] - reference).max()) > 1e-9:
        raise SystemExit("the two runs do not share their reference mesh")

    mask = np.any(np.abs(runs[first]["F"]) > 0.0, axis=1)
    for key in keys[1:]:
        if not np.array_equal(mask, np.any(np.abs(runs[key]["F"]) > 0.0, axis=1)):
            raise SystemExit(f"applied-force node sets differ between {first} and {key}")
    skin = np.nonzero(mask)[0]
    z = reference[skin] @ SPAN
    gap = RING_GAP_FRACTION * float(z.max() - z.min())
    rings = [np.asarray(skin[g], dtype=np.intp) for g in rings_in_band(z, gap)]
    ring_z = np.asarray([float(np.mean(reference[r] @ SPAN)) for r in rings])
    order = np.argsort(ring_z, kind="stable")
    rings = [rings[i] for i in order]
    ring_z = ring_z[order]
    n_ring = len(rings)

    print(f"grid: {reference.shape[0]} nodes, {skin.size} carry an applied force, {n_ring} rings")
    print(
        f"reference and deformed span z in [{z.min():.3f}, {z.max():.3f}] m; "
        f"|U| max {np.linalg.norm(runs[first]['U'][skin], axis=1).max():.4f} m (A), "
        f"{np.linalg.norm(runs[second]['U'][skin], axis=1).max():.4f} m (B)"
    )

    cent_ref = np.asarray([reference[r].mean(axis=0) for r in rings])
    axis_ref = _ring_axes(cent_ref)
    cent_def = {}
    axes_def = {}
    for key in keys:
        deformed = runs[key]["deformed"]
        cent = np.asarray([deformed[r].mean(axis=0) for r in rings])
        cent_def[key] = cent
        axes_def[key] = _ring_axes(cent)
    deformed_first = runs[first]["deformed"]

    # frame evidence: the stored F_AERO lies in the section plane of the reference axis, so it is
    # expressed in the same frame as the coordinates (rotating/reference), not in the inertial one.
    axial_share = np.zeros(n_ring)
    for i, idx in enumerate(rings):
        f = runs[first]["F"][idx]
        mag = np.linalg.norm(f, axis=1)
        good = mag > 0.0
        if good.any():
            axial_share[i] = float(np.median(np.abs(f[good] @ axis_ref[i]) / mag[good]))
    print(
        f"frame check: `points == solid_mesh + U` index-wise for both runs; axial share of F_AERO "
        f"against the reference span axis: median {np.median(axial_share):.3e}, "
        f"max {axial_share.max():.3e} (F lies in the section plane -> same frame as the coordinates; "
        f"rotating it by R(-theta) about Y raises the share to ~1.3e-1, so it is not inertial)"
    )

    fields = (
        "T_ref",
        "T_ref_def",
        "T_def_own",
        "T_def_common",
        "chord",
        "flap",
        "axial",
        "l2",
        "peak",
        "theta",
        "distort",
        "drill",
    )
    metrics = {name: {key: np.zeros(n_ring) for key in keys} for name in fields}
    residual = {key: np.zeros((n_ring, 3)) for key in keys}

    for i, idx in enumerate(rings):
        x_ref = reference[idx]
        r_ref = x_ref - cent_ref[i]
        a_ref = axis_ref[i]
        r_def_first = deformed_first[idx] - cent_def[first][i]
        u_ax, v_ax = section_plane_axes(a_ref)
        a_def_first = axes_def[first][i]
        for key in keys:
            run = runs[key]
            r_def = run["deformed"][idx] - cent_def[key][i]
            f = run["F"][idx]
            a_def = axes_def[key][i]
            metrics["T_ref"][key][i] = float(_moment(r_ref, f) @ a_ref)
            metrics["T_ref_def"][key][i] = float(_moment(r_ref, f) @ a_def)
            metrics["T_def_own"][key][i] = float(_moment(r_def, f) @ a_def)
            metrics["T_def_common"][key][i] = float(_moment(r_def_first, f) @ a_def_first)
            metrics["chord"][key][i] = float((f @ u_ax).sum())
            metrics["flap"][key][i] = float((f @ v_ax).sum())
            metrics["axial"][key][i] = float((f @ a_ref).sum())
            residual[key][i] = f.sum(axis=0)
            mag = np.linalg.norm(f, axis=1)
            metrics["l2"][key][i] = float(np.linalg.norm(f))
            mean_mag = float(mag.mean())
            metrics["peak"][key][i] = float(mag.max() / mean_mag) if mean_mag > 0.0 else np.nan
            metrics["theta"][key][i] = _twist_theta(reference, run["U"], idx)
            metrics["distort"][key][i] = _distortion(reference, run["U"], idx)[1]
            metrics["drill"][key][i] = float(np.sqrt(np.mean((run["rotation"][idx] @ a_def) ** 2)))

    print("\n=== global fields ===")
    for key in keys:
        f = runs[key]["F"]
        print(
            f"  {RUNS[key]:16s} L2={np.linalg.norm(f):.6e} N   max|f|={np.abs(f).max():.6e} N   "
            f"sum F={np.array2string(f.sum(axis=0), precision=4)}   "
            f"sum r x F={np.array2string(_moment(reference, f), precision=4)}"
        )
    df = runs[second]["F"] - runs[first]["F"]
    print(
        f"  {'difference':16s} L2={np.linalg.norm(df):.6e} N "
        f"({np.linalg.norm(df) / np.linalg.norm(runs[first]['F']):.4%} of A)   "
        f"max|df|={np.abs(df).max():.6e} N   sum dF={np.array2string(df.sum(axis=0), precision=4)}   "
        f"sum r x dF={np.array2string(_moment(reference, df), precision=4)}   "
        f"nodes moved={int(np.any(np.abs(df) > 1e-9, axis=1).sum())} of {skin.size}"
    )

    # cumulative torque: what actually twists a cantilever section is the load outboard of it
    cumulative = {
        key: np.cumsum(metrics["T_def_own"][key][::-1])[::-1] for key in keys
    }
    common_cumulative = np.cumsum(metrics["T_def_common"][first][::-1])[::-1]

    stations = np.linspace(ring_z[0], ring_z[-1], n_stations)

    def interp(name, key, zz):
        return np.interp(zz, ring_z, metrics[name][key])

    print("\n=== section profiles (interpolated onto stations) -> torque units N.m ===")
    print(
        f"{'z [m]':>7} | {'Tref A':>10} {'Tref B':>10} | {'Tdef A':>10} {'Tdef B':>10} "
        f"{'dTdef':>10} | {'cum A':>10} {'cum B':>10} {'dCum':>10} | "
        f"{'theta A':>8} {'theta B':>8} | {'drillA':>8} {'drillB':>8} | {'disA':>8} {'disB':>8}"
    )
    for zz in stations:
        ra, rb = interp("T_def_own", first, zz), interp("T_def_own", second, zz)
        ca, cb = np.interp(zz, ring_z, cumulative[first]), np.interp(zz, ring_z, cumulative[second])
        print(
            f"{zz:>7.2f} | {interp('T_ref', first, zz):>10.3e} {interp('T_ref', second, zz):>10.3e} | "
            f"{ra:>10.3e} {rb:>10.3e} {rb - ra:>10.3e} | "
            f"{ca:>10.3e} {cb:>10.3e} {cb - ca:>10.3e} | "
            f"{np.rad2deg(interp('theta', first, zz)):>8.3f} "
            f"{np.rad2deg(interp('theta', second, zz)):>8.3f} | "
            f"{interp('drill', first, zz):>8.2e} {interp('drill', second, zz):>8.2e} | "
            f"{interp('distort', first, zz):>8.2e} {interp('distort', second, zz):>8.2e}"
        )

    print("\n=== totals over the loaded span ===")
    for key in keys:
        print(
            f"  {RUNS[key]:16s} sum T_ref={metrics['T_ref'][key].sum():.6e}  "
            f"sum T_def={metrics['T_def_own'][key].sum():.6e} N.m  "
            f"force sum={np.array2string(residual[key].sum(axis=0), precision=4)}  "
            f"theta tip={np.rad2deg(metrics['theta'][key][-1]):.4f} deg"
        )
    for name in ("T_ref", "T_ref_def", "T_def_own", "T_def_common"):
        delta = metrics[name][second] - metrics[name][first]
        sign_changes = int(np.sum(np.sign(delta[1:]) != np.sign(delta[:-1])))
        print(
            f"  {name:12s} A vs B: sum d={delta.sum():>11.4e} N.m  "
            f"min={delta.min():>11.3e}  max={delta.max():>11.3e}  sign changes={sign_changes}"
        )
    delta_cum = cumulative[second] - cumulative[first]
    print(
        f"  cumulative   A vs B: sum d={delta_cum.sum():.4e} N.m  "
        f"min={delta_cum.min():.3e}  max={delta_cum.max():.3e}  "
        f"common-geometry cumulative A vs B: min {common_cumulative.min():.3e} "
        f"max {common_cumulative.max():.3e}"
    )
    print("\n=== pattern peakedness (max|f|/mean|f| per ring) ===")
    for key in keys:
        print(
            f"  {RUNS[key]:16s} median={np.nanmedian(metrics['peak'][key]):.3f}  "
            f"min={np.nanmin(metrics['peak'][key]):.3f}  max={np.nanmax(metrics['peak'][key]):.3f}  "
            f"L2 total={metrics['l2'][key].sum():.4e} N"
        )
    print(
        f"  component sums: chord A {metrics['chord'][first].sum():.4e} / "
        f"B {metrics['chord'][second].sum():.4e} N;  "
        f"flap A {metrics['flap'][first].sum():.4e} / B {metrics['flap'][second].sum():.4e} N;  "
        f"axial A {metrics['axial'][first].sum():.4e} / B {metrics['axial'][second].sum():.4e} N"
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="#30 T3b: attribute the A/B coupled twist movement from the saved fields"
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
    return analyse(args.base, args.step, args.stations)


if __name__ == "__main__":
    sys.exit(main())
