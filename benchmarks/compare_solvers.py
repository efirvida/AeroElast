#!/usr/bin/env python3
"""
benchmarks/compare_solvers.py
==============================
A/B benchmark: Inertial solver vs. corotational-style full-rebuild per window.

Measures the performance impact of the Sprint 4 optimisation
(`update_node_coordinates` fast path) compared to the legacy approach
(full assembler rebuild at every rotation step).

Requires:
  - petsc4py
  - _aeroelast  (Rust backend)
  - aeroelast package on PYTHONPATH  (or installed in venv)

Usage::

    python benchmarks/compare_solvers.py
    python benchmarks/compare_solvers.py --nx 8 --ny 8 --n-windows 50
    python benchmarks/compare_solvers.py --help

Output:
  - Markdown summary table to stdout
  - PNG plots saved to benchmarks/results/ (if matplotlib available)
  - JSON results saved to benchmarks/results/compare_solvers_<timestamp>.json

Decision criterion (from IMPLEMENTATION_PLAN_INERTIAL.md Sprint 5):
  IF  speedup >= 1.5x  AND  K_eigenvalue_error < 1e-6
  THEN  update_node_coordinates fast path is production-ready.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import NamedTuple

# ---------------------------------------------------------------------------
# Path setup — allow running from repo root without install
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))

# ---------------------------------------------------------------------------
# Imports (with informative errors)
# ---------------------------------------------------------------------------
try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise SystemExit("numpy is required: pip install numpy") from exc

try:
    from petsc4py import PETSc
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "petsc4py is required — load the appropriate HPC module first."
    ) from exc

try:
    import _aeroelast  # noqa: F401
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        f"Rust backend (_aeroelast) not available: {exc}\n"
        "Build it first: cd crates/aeroelast-py && maturin develop --release"
    ) from exc

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.material import IsotropicMaterial, OrthotropicMaterial
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.elements import ElementFamily

# ---------------------------------------------------------------------------
# Constants / materials
# ---------------------------------------------------------------------------

_ISO_MAT = IsotropicMaterial(name="Al", E=70e9, nu=0.33, rho=2700.0)
_ORTHO_MAT = OrthotropicMaterial(
    name="CFRP",
    E=(150e9, 10e9, 10e9),
    G=(5e9, 5e9, 5e9),
    nu=(0.3, 0.3, 0.3),
    rho=1600.0,
)
_THICKNESS = 0.004  # 4 mm


# ---------------------------------------------------------------------------
# Mesh helpers
# ---------------------------------------------------------------------------

def _build_quad_plate(nx: int, ny: int, L: float = 1.0) -> MeshModel:
    """Flat MITC4 plate mesh (nx × ny quads)."""
    Node._id_counter = 0
    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, L, ny + 1)
    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            n = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(n)
            grid[(i, j)] = n
    for j in range(ny):
        for i in range(nx):
            mesh.add_element(
                MeshElement(
                    nodes=[
                        grid[(i, j)],
                        grid[(i + 1, j)],
                        grid[(i + 1, j + 1)],
                        grid[(i, j + 1)],
                    ],
                    element_type=ElementType.quad,
                )
            )
    return mesh


def _model_cfg(mat):
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "material": mat,
            "thickness": _THICKNESS,
        }
    }


def _rotate_coords(coords: np.ndarray, theta: float) -> np.ndarray:
    """Rotate (n,3) coords by theta around Z-axis."""
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return (R @ coords.T).T


def _petsc_to_dense(mat: PETSc.Mat) -> np.ndarray:
    n = mat.getSize()[0]
    dense = np.zeros((n, n))
    for i in range(n):
        cols, vals = mat.getRow(i)
        dense[i, cols] = vals
    return dense


# ---------------------------------------------------------------------------
# Benchmark result container
# ---------------------------------------------------------------------------

class BenchmarkResult(NamedTuple):
    label: str
    material: str
    n_nodes: int
    n_dofs: int
    n_windows: int
    time_init_s: float           # one-time assembler init
    time_per_window_s: float     # median window time (geometry update + assembly)
    time_total_s: float          # total for all windows
    eig_error_max: float         # max |λ(rotated) - λ(orig)| / λ_max(orig)
    speedup: float               # ratio vs reference (set later)


# ---------------------------------------------------------------------------
# Approach A: "corotational" — rebuild MeshAssembler from scratch each window
# ---------------------------------------------------------------------------

def _bench_rebuild(
    nx: int, ny: int, mat, n_windows: int, thetas: np.ndarray
) -> BenchmarkResult:
    """
    Corotational / legacy approach: re-initialise MeshAssembler at each window
    with a freshly built MeshModel containing rotated node coordinates.
    This is the baseline — it does maximum work every step.
    """
    label = "rebuild"
    mesh0 = _build_quad_plate(nx, ny)
    coords0 = np.stack([n.coords[:3] for n in mesh0.nodes], axis=0)

    # One-time cost: first assembly
    t0 = time.perf_counter()
    asm0 = MeshAssembler(mesh=_build_quad_plate(nx, ny), model=_model_cfg(mat))
    K0 = asm0.assemble_stiffness_matrix()
    t_init = time.perf_counter() - t0
    eigs0 = np.sort(np.linalg.eigvalsh(_petsc_to_dense(K0)))

    n_nodes = asm0.n_nodes
    n_dofs = n_nodes * 6

    # Per-window: rebuild mesh + assembler + assemble K
    window_times: list[float] = []
    eig_errors: list[float] = []

    for theta in thetas:
        coords_rot = _rotate_coords(coords0, theta)

        t0 = time.perf_counter()
        # Rebuild entire MeshModel + MeshAssembler
        Node._id_counter = 0
        mesh_rot = MeshModel()
        grid: dict[int, Node] = {}
        for idx, c in enumerate(coords_rot):
            n = Node(c.tolist(), geometric_node=False)
            mesh_rot.add_node(n)
            grid[idx] = n
        # Reconstruct element connectivity from original mesh
        orig_elements = list(asm0.elements)
        for elem in orig_elements:
            if hasattr(elem, "nodes"):
                # Get node indices and rebuild
                pass
        # Simpler: just build the same topology with rotated coords
        Node._id_counter = 0
        mesh_new = _build_quad_plate(nx, ny)
        for node, c_rot in zip(mesh_new.nodes, coords_rot):
            node.coords[:3] = c_rot
        asm_rot = MeshAssembler(mesh=mesh_new, model=_model_cfg(mat))
        K_rot = asm_rot.assemble_stiffness_matrix()
        dt = time.perf_counter() - t0

        window_times.append(dt)

        eigs_rot = np.sort(np.linalg.eigvalsh(_petsc_to_dense(K_rot)))
        eig_err = np.max(np.abs(eigs_rot - eigs0)) / (np.max(np.abs(eigs0)) + 1e-30)
        eig_errors.append(eig_err)

    t_total = sum(window_times)
    t_median = float(np.median(window_times))
    eig_error_max = float(np.max(eig_errors))

    return BenchmarkResult(
        label=label,
        material=mat.name,
        n_nodes=n_nodes,
        n_dofs=n_dofs,
        n_windows=n_windows,
        time_init_s=t_init,
        time_per_window_s=t_median,
        time_total_s=t_total,
        eig_error_max=eig_error_max,
        speedup=1.0,  # baseline
    )


# ---------------------------------------------------------------------------
# Approach B: "inertial" — update_node_coordinates fast path
# ---------------------------------------------------------------------------

def _bench_update(
    nx: int, ny: int, mat, n_windows: int, thetas: np.ndarray
) -> BenchmarkResult:
    """
    Inertial / Sprint 4 approach: build assembler once, then call
    update_node_coordinates() at each window — no Python-side mesh rebuild.
    """
    label = "update_coords"

    mesh0 = _build_quad_plate(nx, ny)
    coords0 = np.stack([n.coords[:3] for n in mesh0.nodes], axis=0)

    # One-time cost: first assembly
    t0 = time.perf_counter()
    asm = MeshAssembler(mesh=mesh0, model=_model_cfg(mat))
    K0 = asm.assemble_stiffness_matrix()
    t_init = time.perf_counter() - t0

    if asm._rust is None:
        raise RuntimeError("Rust assembler not available — cannot use update_node_coordinates")

    eigs0 = np.sort(np.linalg.eigvalsh(_petsc_to_dense(K0)))
    n_nodes = asm.n_nodes
    n_dofs = n_nodes * 6

    # Per-window: update coords in Rust + reassemble
    window_times: list[float] = []
    eig_errors: list[float] = []

    for theta in thetas:
        coords_rot = _rotate_coords(coords0, theta)

        t0 = time.perf_counter()
        asm._rust.update_node_coordinates(coords_rot)
        K_rot = asm.assemble_stiffness_matrix()
        dt = time.perf_counter() - t0

        window_times.append(dt)

        eigs_rot = np.sort(np.linalg.eigvalsh(_petsc_to_dense(K_rot)))
        eig_err = np.max(np.abs(eigs_rot - eigs0)) / (np.max(np.abs(eigs0)) + 1e-30)
        eig_errors.append(eig_err)

    # Restore original geometry
    asm._rust.update_node_coordinates(coords0)

    t_total = sum(window_times)
    t_median = float(np.median(window_times))
    eig_error_max = float(np.max(eig_errors))

    return BenchmarkResult(
        label=label,
        material=mat.name,
        n_nodes=n_nodes,
        n_dofs=n_dofs,
        n_windows=n_windows,
        time_init_s=t_init,
        time_per_window_s=t_median,
        time_total_s=t_total,
        eig_error_max=eig_error_max,
        speedup=1.0,  # filled in by caller
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def _print_markdown_table(results: list[tuple[BenchmarkResult, BenchmarkResult]]) -> None:
    print("\n## A/B Benchmark: inertial solver vs. full-rebuild")
    print()
    header = (
        "| Material | Nodes | DOFs | Approach | Init (s) | "
        "Median window (s) | Total (s) | K eig error | Speedup |"
    )
    sep = (
        "|----------|-------|------|----------|----------|"
        "-------------------|-----------|-------------|---------|"
    )
    print(header)
    print(sep)
    for rebuild, update in results:
        for r in (rebuild, update):
            print(
                f"| {r.material} | {r.n_nodes} | {r.n_dofs} "
                f"| `{r.label}` "
                f"| {r.time_init_s:.4f} "
                f"| {r.time_per_window_s:.4f} "
                f"| {r.time_total_s:.4f} "
                f"| {r.eig_error_max:.2e} "
                f"| {r.speedup:.2f}× |"
            )
    print()


def _decision(results: list[tuple[BenchmarkResult, BenchmarkResult]]) -> None:
    print("## Decision criterion")
    print()
    print("  speedup >= 1.5×  AND  K eigenvalue error < 1e-6")
    print()
    all_pass = True
    for rebuild, update in results:
        speedup_ok = update.speedup >= 1.5
        eig_ok = update.eig_error_max < 1e-6
        status = "✓ PROMOTE" if (speedup_ok and eig_ok) else "✗ KEEP EXPERIMENTAL"
        print(
            f"  {update.material}: speedup={update.speedup:.2f}×  "
            f"eig_err={update.eig_error_max:.2e}  → {status}"
        )
        if not (speedup_ok and eig_ok):
            all_pass = False
    print()
    if all_pass:
        print("**VERDICT: update_node_coordinates fast path is PRODUCTION-READY.**")
    else:
        print("**VERDICT: Fast path does not meet all criteria — keep as EXPERIMENTAL.**")
    print()


def _save_json(
    results: list[tuple[BenchmarkResult, BenchmarkResult]],
    out_dir: Path,
) -> None:
    import datetime
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = out_dir / f"compare_solvers_{ts}.json"
    data = []
    for rebuild, update in results:
        data.append({"rebuild": rebuild._asdict(), "update_coords": update._asdict()})
    out_path.write_text(json.dumps(data, indent=2))
    print(f"Results saved to: {out_path}")


def _save_plots(
    results: list[tuple[BenchmarkResult, BenchmarkResult]],
    out_dir: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available — skipping plots")
        return

    materials = [r.material for r, _ in results]
    t_rebuild = [r.time_per_window_s for r, _ in results]
    t_update  = [u.time_per_window_s for _, u in results]
    speedups  = [u.speedup for _, u in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    x = range(len(materials))
    ax1.bar([i - 0.2 for i in x], t_rebuild, width=0.4, label="rebuild (corotational)")
    ax1.bar([i + 0.2 for i in x], t_update, width=0.4, label="update_coords (inertial)")
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(materials)
    ax1.set_ylabel("Median time per window (s)")
    ax1.set_title("Assembly time per window")
    ax1.legend()
    ax1.axhline(0, color="k", linewidth=0.5)

    ax2.bar(list(x), speedups, color=["green" if s >= 1.5 else "orange" for s in speedups])
    ax2.axhline(1.5, color="red", linestyle="--", label="1.5× threshold")
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(materials)
    ax2.set_ylabel("Speedup (rebuild / update_coords)")
    ax2.set_title("Speedup: inertial vs. corotational rebuild")
    ax2.legend()

    fig.tight_layout()
    out_path = out_dir / "compare_solvers.png"
    fig.savefig(out_path, dpi=150)
    print(f"Plot saved to: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=6, help="Mesh divisions X (default 6)")
    parser.add_argument("--ny", type=int, default=6, help="Mesh divisions Y (default 6)")
    parser.add_argument(
        "--n-windows", type=int, default=20,
        help="Number of FSI windows (rotation steps) to time (default 20)"
    )
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip matplotlib plots"
    )
    args = parser.parse_args()

    nx, ny = args.nx, args.ny
    n_windows = args.n_windows

    # Rotation angles: evenly spaced over one full revolution
    thetas = np.linspace(0.0, 2 * np.pi, n_windows, endpoint=False)

    print(f"Benchmark: nx={nx}, ny={ny}, n_windows={n_windows}")
    print(f"  Mesh: {(nx+1)*(ny+1)} nodes, {(nx+1)*(ny+1)*6} DOFs")
    print()

    out_dir = _REPO_ROOT / "benchmarks" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_results: list[tuple[BenchmarkResult, BenchmarkResult]] = []

    for mat in (_ISO_MAT, _ORTHO_MAT):
        print(f"  Running {mat.name} …", end=" ", flush=True)

        r_rebuild = _bench_rebuild(nx, ny, mat, n_windows, thetas)
        r_update  = _bench_update(nx, ny, mat, n_windows, thetas)

        # Compute speedup
        speedup = (
            r_rebuild.time_per_window_s / r_update.time_per_window_s
            if r_update.time_per_window_s > 0
            else float("inf")
        )
        r_update = r_update._replace(speedup=speedup)

        all_results.append((r_rebuild, r_update))
        print(f"done  ({speedup:.2f}× speedup)")

    _print_markdown_table(all_results)
    _decision(all_results)
    _save_json(all_results, out_dir)

    if not args.no_plots:
        _save_plots(all_results, out_dir)


if __name__ == "__main__":
    main()
