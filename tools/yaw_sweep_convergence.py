#!/usr/bin/env python3
"""
yaw_sweep_convergence.py — preCICE coupling convergence analysis
================================================================

Reads preCICE log files (``precice-Solid-iterations.log`` and optionally
``precice-Solid-convergence.log``) for each yaw case and reports:

  * Iterations per time-window (mean, max, distribution)
  * Convergence rate evolution over simulation time
  * Residual decay statistics (Displacement and Force)
  * Cases where max-iterations was hit (convergence failures)

The time-window index is converted to simulation time using
``t = window * dt``, where ``dt`` is the preCICE time-window size.

Usage
-----
  python tools/yaw_sweep_convergence.py \\
        --results-dir /scratch/.../frontiersin_results \\
        [--dt 0.01] [--t-start 20] [--t-end 70] \\
        [--output convergence_summary.csv] \\
        [--verbose]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def parse_iterations_log(path: Path) -> list[dict]:
    """Parse precice-*-iterations.log into list of dicts."""
    rows = []
    with open(path) as f:
        lines = f.readlines()
    # Find header line (starts with whitespace + "TimeWindow")
    data_start = 0
    for i, line in enumerate(lines):
        if "TimeWindow" in line:
            data_start = i + 1
            break
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) >= 3:
            try:
                rows.append({
                    "time_window": int(parts[0]),
                    "total_iterations": int(parts[1]),
                    "iterations": int(parts[2]),
                    "convergence": int(parts[3]) if len(parts) > 3 else None,
                    "qn_columns": int(parts[4]) if len(parts) > 4 else None,
                })
            except (ValueError, IndexError):
                continue
    return rows


def parse_convergence_log(path: Path) -> list[dict]:
    """Parse precice-*-convergence.log into list of dicts."""
    rows = []
    with open(path) as f:
        lines = f.readlines()
    # Find header
    data_start = 0
    header_cols: list[str] = []
    for i, line in enumerate(lines):
        if "TimeWindow" in line and "Iteration" in line:
            # Extract field names from header
            header_cols = line.split()
            data_start = i + 1
            break
    if not header_cols:
        return rows
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) != len(header_cols):
            continue
        try:
            row = {"time_window": int(parts[0]), "iteration": int(parts[1])}
            for col, val in zip(header_cols[2:], parts[2:]):
                row[col] = float(val)
            rows.append(row)
        except (ValueError, IndexError):
            continue
    return rows


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------


def window_to_time(window: int, dt: float) -> float:
    return window * dt


def filter_windows(rows: list[dict], t_start: float, t_end: float, dt: float) -> list[dict]:
    return [r for r in rows if t_start < window_to_time(r["time_window"], dt) <= t_end]


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _std(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def analyze_iterations(rows: list[dict]) -> dict:
    iters = [r["iterations"] for r in rows]
    conv = [r["convergence"] for r in rows if r.get("convergence") is not None]
    return {
        "n_windows": len(rows),
        "iter_mean": _mean(iters),
        "iter_std": _std(iters),
        "iter_max": max(iters) if iters else 0,
        "iter_min": min(iters) if iters else 0,
        "iter_total": sum(iters),
        "no_conv_windows": sum(1 for c in conv if c == 0),
        "iter_hist": {i: iters.count(i) for i in sorted(set(iters))},
    }


def analyze_residuals(conv_rows: list[dict]) -> dict:
    """Compute final residual per window and decay stats."""
    # Group by time_window, take last iteration's residuals
    by_window: dict[int, list[dict]] = {}
    for r in conv_rows:
        by_window.setdefault(r["time_window"], []).append(r)

    # residual field names (all except time_window, iteration)
    if not conv_rows:
        return {}
    res_fields = [k for k in conv_rows[0] if k not in ("time_window", "iteration")]

    final_residuals: dict[str, list[float]] = {f: [] for f in res_fields}
    for tw_rows in by_window.values():
        last = max(tw_rows, key=lambda r: r["iteration"])
        for f in res_fields:
            if f in last:
                final_residuals[f].append(last[f])

    stats = {}
    for f, vals in final_residuals.items():
        stats[f] = {
            "mean_final": _mean(vals),
            "max_final": max(vals) if vals else float("nan"),
            "above_1e4": sum(1 for v in vals if v > 1e-4),
        }
    return stats


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def print_convergence_table(
    cases: list[tuple[int, dict, dict]],
    t_start: float,
    t_end: float,
) -> None:
    sep = "─" * 90

    print()
    print("  PRECICE COUPLING CONVERGENCE ANALYSIS")
    print(f"  Window: {t_start} < t <= {t_end} s")
    print(sep)

    print(
        f"  {'Yaw':>5}  {'Windows':>9}  {'Iter/W':>10}  {'Iter max':>9}  "
        f"{'Iter total':>11}  {'No-conv':>8}"
    )
    print(f"  {'[°]':>5}  {'in window':>9}  {'mean±std':>10}  {'':>9}  {'':>11}  {'windows':>8}")
    print(sep)

    for yaw_deg, iter_stats, res_stats in sorted(cases, key=lambda x: x[0]):
        i_mean = iter_stats["iter_mean"]
        i_std = iter_stats["iter_std"]
        i_max = iter_stats["iter_max"]
        i_tot = iter_stats["iter_total"]
        no_cv = iter_stats["no_conv_windows"]
        nc_str = f"{'⚠' if no_cv > 0 else ''}{no_cv}"
        print(
            f"  {yaw_deg:>5}  {iter_stats['n_windows']:>9}  "
            f"{i_mean:>6.2f}±{i_std:<4.2f}  {i_max:>9}  {i_tot:>11}  {nc_str:>8}"
        )

    print(sep)

    # Residual summary
    any_res = any(res for _, _, res in cases)
    if any_res:
        print()
        print("  FINAL RESIDUALS (last iteration per time-window, mean in ss window)")
        print(sep)
        # collect all field names
        all_fields: set[str] = set()
        for _, _, res in cases:
            all_fields.update(res.keys())
        for field in sorted(all_fields):
            print(f"  Field: {field}")
            print(f"  {'Yaw':>5}  {'Mean final':>14}  {'Max final':>14}  {'Above 1e-4':>12}")
            for yaw_deg, _, res in sorted(cases, key=lambda x: x[0]):
                if field not in res:
                    continue
                s = res[field]
                flag = "⚠" if s["above_1e4"] > 0 else " "
                print(
                    f"  {yaw_deg:>5}  {s['mean_final']:>14.3e}  "
                    f"{s['max_final']:>14.3e}  "
                    f"{flag}{s['above_1e4']:>10} windows"
                )
            print()

    # Iteration distribution
    print()
    print("  ITERATION DISTRIBUTION (fraction of time-windows)")
    print(sep)
    # Collect all seen iteration counts
    all_iter_counts: set[int] = set()
    for _, s, _ in cases:
        all_iter_counts.update(s["iter_hist"].keys())

    header = f"  {'Yaw':>5}"
    for k in sorted(all_iter_counts):
        header += f"  {'iter=' + str(k):>8}"
    print(header)

    for yaw_deg, s, _ in sorted(cases, key=lambda x: x[0]):
        n = s["n_windows"]
        row = f"  {yaw_deg:>5}"
        for k in sorted(all_iter_counts):
            cnt = s["iter_hist"].get(k, 0)
            pct = cnt / n * 100 if n > 0 else 0
            row += f"  {pct:>7.1f}%"
        print(row)

    print(sep)


def save_convergence_csv(cases: list[tuple[int, dict, dict]], output: Path) -> None:
    rows = []
    for yaw_deg, s, res in sorted(cases, key=lambda x: x[0]):
        row: dict = {
            "yaw_deg": yaw_deg,
            "n_windows": s["n_windows"],
            "iter_mean": round(s["iter_mean"], 3),
            "iter_std": round(s["iter_std"], 3),
            "iter_max": s["iter_max"],
            "iter_min": s["iter_min"],
            "iter_total": s["iter_total"],
            "no_conv_windows": s["no_conv_windows"],
        }
        for field, rs in res.items():
            safe = field.replace(":", "_").replace("(", "").replace(")", "").replace(" ", "_")
            row[f"res_{safe}_mean_final"] = f"{rs['mean_final']:.3e}"
            row[f"res_{safe}_max_final"] = f"{rs['max_final']:.3e}"
            row[f"res_{safe}_above_1e4"] = rs["above_1e4"]
        rows.append(row)

    if not rows:
        return

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Saved: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze preCICE coupling convergence for a yaw sweep.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--results-dir",
        "-d",
        default="/scratch/leahk/eduardo.donestevez/frontiersin_results",
        help="Root directory (default: %(default)s)",
    )
    parser.add_argument(
        "--dt", type=float, default=0.01, help="preCICE time-window size [s] (default: %(default)s)"
    )
    parser.add_argument(
        "--t-start",
        type=float,
        default=20.0,
        help="Analysis window start [s] (default: %(default)s)",
    )
    parser.add_argument(
        "--t-end", type=float, default=70.0, help="Analysis window end [s] (default: %(default)s)"
    )
    parser.add_argument("--output", "-o", default=None, help="Save summary as CSV")
    parser.add_argument("--verbose", action="store_true", help="Show per-window iteration counts")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        print(f"ERROR: {results_dir} not found", file=sys.stderr)
        sys.exit(1)

    print(f"\n  Scanning: {results_dir}")
    print(f"  dt={args.dt}s  |  window: {args.t_start} < t <= {args.t_end} s\n")

    cases: list[tuple[int, dict, dict]] = []

    for yaw_dir in sorted(results_dir.iterdir()):
        if not yaw_dir.is_dir() or not yaw_dir.name.startswith("yaw_"):
            continue
        try:
            yaw_deg = int(yaw_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue

        iter_log = yaw_dir / "precice-Solid-iterations.log"
        conv_log = yaw_dir / "precice-Solid-convergence.log"

        if not iter_log.exists():
            print(f"  [SKIP] yaw={yaw_deg}°: no iterations log found", file=sys.stderr)
            continue

        all_iter_rows = parse_iterations_log(iter_log)
        iter_rows = filter_windows(all_iter_rows, args.t_start, args.t_end, args.dt)

        if not iter_rows:
            print(f"  [SKIP] yaw={yaw_deg}°: no windows in analysis band", file=sys.stderr)
            continue

        t_span = (
            max(r["time_window"] for r in iter_rows) - min(r["time_window"] for r in iter_rows)
        ) * args.dt

        if t_span < 20.0:
            print(f"  [SKIP] yaw={yaw_deg}°: only {t_span:.1f}s of data in window", file=sys.stderr)
            continue

        i_stats = analyze_iterations(iter_rows)

        res_stats: dict = {}
        if conv_log.exists():
            all_conv_rows = parse_convergence_log(conv_log)
            conv_rows = filter_windows(all_conv_rows, args.t_start, args.t_end, args.dt)
            res_stats = analyze_residuals(conv_rows)

        cases.append((yaw_deg, i_stats, res_stats))
        print(
            f"  [OK]   yaw={yaw_deg:>3}°: {i_stats['n_windows']} windows, "
            f"iter mean={i_stats['iter_mean']:.2f}, max={i_stats['iter_max']}"
        )

        if args.verbose:
            print(f"         distribution: {i_stats['iter_hist']}")

    if not cases:
        print("\nNo cases processed.", file=sys.stderr)
        sys.exit(1)

    print_convergence_table(cases, args.t_start, args.t_end)

    if args.output:
        save_convergence_csv(cases, Path(args.output))


if __name__ == "__main__":
    main()
