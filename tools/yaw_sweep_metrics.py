#!/usr/bin/env python3
"""
yaw_sweep_metrics.py — Yaw-sweep FSI post-processing tool
==========================================================

Reads ``bem_report.csv`` files produced by the BEM-FSI fluid participant and
computes steady-state statistics for each yaw angle case.

Time window:
  Only rows with  ``t_start < t <= t_end``  are used.  Any case that has
  fewer than ``min_window`` seconds in that band is flagged and skipped so
  incomplete runs don't pollute the summary.

Usage
-----
  python tools/yaw_sweep_metrics.py \\
        --results-dir /scratch/.../frontiersin_results \\
        [--t-start 20] [--t-end 70] [--min-window 20] \\
        [--solver corotational] \\
        [--output metrics_summary.csv] \\
        [--no-table]

Output
------
  * Console table (unless --no-table)
  * CSV file with one row per case (optional, --output)

Published comparison targets (IEA 15 MW, rated, yaw sweep)
-----------------------------------------------------------
  Source: Ma et al. (2025), doi:10.3389/fenrg.2025.1571567
  Figures 11, 15a, 16a — corotational solver, yaw = 0–40°

  Yaw=0 flapwise reference:
    Zhou [R11] "Present":  13.86 m
    ALM-GEBT:              14.10 m
    Acceptance band:        9–15 m

  Power law approximation (cos^n):
    P(γ) ≈ P(0)·cos^n(γ)  with  n ≈ 3 for the IEA 15 MW at rated
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional

import numpy as np

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class CaseMetrics:
    yaw_deg: int
    solver: str
    t_start_actual: float
    t_end_actual: float
    n_samples: int
    # --- tip displacements ---
    flap_mean: float  # Tip Disp Y [m]  flapwise
    flap_std: float
    flap_max: float
    flap_min: float
    edge_mean: float  # Tip Disp X [m]  edgewise
    edge_std: float
    edge_max: float
    edge_min: float
    span_mean: float  # Tip Disp Z [m]  spanwise / chordwise
    span_std: float
    tip_mag_mean: float  # Tip Disp Mag [m]
    tip_mag_std: float
    tip_mag_max: float
    max_disp_mean: float  # Max structural displacement
    # --- aerodynamics ---
    thrust_mean: float  # Thrust [N]
    thrust_std: float
    torque_mean: float  # Torque [N.m]
    torque_std: float
    power_mean_mw: float  # Power [MW]
    power_std_mw: float
    power_max_mw: float
    ct_mean: float  # CT
    ct_std: float
    cq_mean: float  # CQ
    cp_mean: float  # CP
    cp_std: float
    # --- structural ---
    mb_mean: float  # Mb [N.m]  root bending moment
    mb_std: float
    max_force_mean: float  # Max Nodal Force [N]
    # --- oscillation amplitude (3P proxy) ---
    flap_p2p: float  # max - min of flapwise in ss window
    edge_p2p: float
    power_p2p: float
    # --- divergence flag ---
    diverged: bool = False
    # --- cos^n fit residual vs yaw=0 ---
    cp_cos3_predicted: Optional[float] = None  # filled in post-hoc
    cp_cos3_err_pct: Optional[float] = None


# ---------------------------------------------------------------------------
# Core helpers
# ---------------------------------------------------------------------------


_DIVERGE_THRUST_N = 1e9  # 1 GN — physically impossible for any wind turbine


def _mean(xs: list[float]) -> float:
    a = np.asarray(xs, dtype=float)
    return float(np.nanmean(a))


def _std(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    a = np.asarray(xs, dtype=float)
    if not np.all(np.isfinite(a)):
        return float("nan")
    return float(np.std(a, ddof=1))


def load_bem_csv(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def filter_window(rows: list[dict], t_start: float, t_end: float) -> list[dict]:
    return [r for r in rows if t_start < r["Time [s]"] <= t_end]


def compute_metrics(
    yaw_deg: int,
    solver: str,
    rows: list[dict],
) -> CaseMetrics:
    def col(name: str) -> list[float]:
        return [r[name] for r in rows]

    flap = col("Tip Disp Y [m]")
    edge = col("Tip Disp X [m]")
    span = col("Tip Disp Z [m]")
    mag = col("Tip Disp Mag [m]")
    mdisp = col("Max Disp [m]")
    thr = col("Thrust [N]")
    torq = col("Torque [N.m]")
    pwr = [v / 1e6 for v in col("Power [W]")]
    ct = col("CT")
    cq = col("CQ")
    cp = col("CP")
    mb = col("Mb [N.m]")
    mf = col("Max Nodal Force [N]")
    t = col("Time [s]")

    # Detect diverged cases by extreme thrust values (physically impossible)
    thr_arr = np.asarray(thr, dtype=float)
    diverged = bool(np.any(np.abs(thr_arr) > _DIVERGE_THRUST_N))
    if diverged:
        nan = float("nan")
        return CaseMetrics(
            yaw_deg=yaw_deg,
            solver=solver,
            t_start_actual=min(t),
            t_end_actual=max(t),
            n_samples=len(rows),
            flap_mean=nan, flap_std=nan, flap_max=nan, flap_min=nan,
            edge_mean=nan, edge_std=nan, edge_max=nan, edge_min=nan,
            span_mean=nan, span_std=nan,
            tip_mag_mean=nan, tip_mag_std=nan, tip_mag_max=nan,
            max_disp_mean=nan,
            thrust_mean=nan, thrust_std=nan,
            torque_mean=nan, torque_std=nan,
            power_mean_mw=nan, power_std_mw=nan, power_max_mw=nan,
            ct_mean=nan, ct_std=nan, cq_mean=nan, cp_mean=nan, cp_std=nan,
            mb_mean=nan, mb_std=nan, max_force_mean=nan,
            flap_p2p=nan, edge_p2p=nan, power_p2p=nan,
            diverged=True,
        )

    return CaseMetrics(
        yaw_deg=yaw_deg,
        solver=solver,
        t_start_actual=min(t),
        t_end_actual=max(t),
        n_samples=len(rows),
        flap_mean=_mean(flap),
        flap_std=_std(flap),
        flap_max=max(flap),
        flap_min=min(flap),
        edge_mean=_mean(edge),
        edge_std=_std(edge),
        edge_max=max(edge),
        edge_min=min(edge),
        span_mean=_mean(span),
        span_std=_std(span),
        tip_mag_mean=_mean(mag),
        tip_mag_std=_std(mag),
        tip_mag_max=max(mag),
        max_disp_mean=_mean(mdisp),
        thrust_mean=_mean(thr),
        thrust_std=_std(thr),
        torque_mean=_mean(torq),
        torque_std=_std(torq),
        power_mean_mw=_mean(pwr),
        power_std_mw=_std(pwr),
        power_max_mw=max(pwr),
        ct_mean=_mean(ct),
        ct_std=_std(ct),
        cq_mean=_mean(cq),
        cp_mean=_mean(cp),
        cp_std=_std(cp),
        mb_mean=_mean(mb),
        mb_std=_std(mb),
        max_force_mean=_mean(mf),
        flap_p2p=max(flap) - min(flap),
        edge_p2p=max(edge) - min(edge),
        power_p2p=max(pwr) - min(pwr),
        diverged=False,
    )


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def discover_cases(results_dir: Path, solver: str) -> list[tuple[int, str, Path]]:
    """Return list of (yaw_deg, solver_label, csv_path) sorted by yaw."""
    cases = []
    for yaw_dir in sorted(results_dir.iterdir()):
        if not yaw_dir.is_dir() or not yaw_dir.name.startswith("yaw_"):
            continue
        try:
            yaw_deg = int(yaw_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue

        # Prefer solver-specific subdir, fall back to bare fluid/
        solver_csv = yaw_dir / solver / "fluid" / "bem_report.csv"
        bare_csv = yaw_dir / "fluid" / "bem_report.csv"

        if solver_csv.exists():
            cases.append((yaw_deg, solver, solver_csv))
        elif bare_csv.exists():
            cases.append((yaw_deg, solver, bare_csv))
        else:
            print(f"  [SKIP] yaw={yaw_deg}°: no bem_report.csv found", file=sys.stderr)

    return cases


# ---------------------------------------------------------------------------
# cos^n fit
# ---------------------------------------------------------------------------


def annotate_cosn(metrics: list[CaseMetrics], n: float = 3.0) -> None:
    """Fill cp_cos3_predicted and cp_cos3_err_pct relative to yaw=0 case."""
    ref = next((m for m in metrics if m.yaw_deg == 0), None)
    if ref is None:
        return
    cp0 = ref.cp_mean
    for m in metrics:
        gamma = math.radians(m.yaw_deg)
        predicted = cp0 * math.cos(gamma) ** n
        m.cp_cos3_predicted = round(predicted, 4)
        if predicted > 0:
            m.cp_cos3_err_pct = round((m.cp_mean - predicted) / predicted * 100, 2)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

_YAW_REFS = {
    # yaw_deg: (flap_ref_m, ref_label)
    0: (13.86, "Zhou [R11]"),
}


def print_table(metrics: list[CaseMetrics]) -> None:
    sep = "─" * 100
    print()
    print("  YAW SWEEP — STEADY-STATE METRICS")
    print(f"  solver: {metrics[0].solver if metrics else '—'}")
    print(sep)

    # Header
    print(
        f"  {'Yaw':>5}  {'Window [s]':>12}  {'n':>6}  "
        f"{'Flap [m]':>14}  {'Edge [m]':>14}  {'Tip|m| [m]':>12}  "
        f"{'CP':>8}  {'CT':>8}  {'Power [MW]':>12}"
    )
    print(
        f"  {'[°]':>5}  {'t_start–t_end':>12}  {'':>6}  "
        f"{'mean ± std':>14}  {'mean ± std':>14}  {'mean':>12}  "
        f"{'mean':>8}  {'mean':>8}  {'mean ± std':>12}"
    )
    print(sep)

    for m in sorted(metrics, key=lambda x: x.yaw_deg):
        window = f"{m.t_start_actual:.1f}–{m.t_end_actual:.1f}"
        flap = f"{m.flap_mean:+.3f}±{m.flap_std:.3f}"
        edge = f"{m.edge_mean:+.3f}±{m.edge_std:.3f}"
        pwr = f"{m.power_mean_mw:.3f}±{m.power_std_mw:.3f}"
        print(
            f"  {m.yaw_deg:>5}  {window:>12}  {m.n_samples:>6}  "
            f"  {flap:>14}  {edge:>14}  {m.tip_mag_mean:>12.3f}  "
            f"{m.ct_mean:>8.4f}  {m.cp_mean:>8.4f}  {pwr:>12}"
        )

    print(sep)

    print()
    print("  OSCILLATION AMPLITUDES (peak-to-peak in steady-state window)")
    print(sep)
    print(f"  {'Yaw':>5}  {'Flap p2p [m]':>14}  {'Edge p2p [m]':>14}  {'Power p2p [MW]':>16}")
    for m in sorted(metrics, key=lambda x: x.yaw_deg):
        print(
            f"  {m.yaw_deg:>5}  {m.flap_p2p:>14.3f}  "
            f"{m.edge_p2p:>14.3f}  {m.power_p2p:>16.3f}"
        )
    print(sep)

    print()
    print("  AERODYNAMICS — THRUST & TORQUE")
    print(sep)
    print(
        f"  {'Yaw':>5}  {'Thrust [kN]':>14}  {'Torque [MN·m]':>16}  "
        f"{'CQ':>8}  {'Root Mb [MN·m]':>16}"
    )
    for m in sorted(metrics, key=lambda x: x.yaw_deg):
        print(
            f"  {m.yaw_deg:>5}  {m.thrust_mean / 1e3:>14.2f}  "
            f"{m.torque_mean / 1e6:>16.3f}  "
            f"{m.cq_mean:>8.5f}  {m.mb_mean / 1e6:>16.3f}"
        )
    print(sep)

    # cos^3 law check
    ref = next((m for m in metrics if m.yaw_deg == 0), None)
    if ref and any(m.cp_cos3_predicted is not None for m in metrics):
        print()
        print("  cos³(γ) LAW CHECK  — CP(γ) vs CP(0)·cos³(γ)")
        print(sep)
        print(f"  {'Yaw':>5}  {'CP_sim':>8}  {'CP_cos3':>8}  {'Err %':>8}")
        for m in sorted(metrics, key=lambda x: x.yaw_deg):
            pred = f"{m.cp_cos3_predicted:.4f}" if m.cp_cos3_predicted is not None else "—"
            err = f"{m.cp_cos3_err_pct:+.1f}%" if m.cp_cos3_err_pct is not None else "—"
            print(f"  {m.yaw_deg:>5}  {m.cp_mean:>8.4f}  {pred:>8}  {err:>8}")
        print(sep)

    # Validation references
    has_refs = any(m.yaw_deg in _YAW_REFS for m in metrics)
    if has_refs:
        print()
        print("  V-04 VALIDATION — FLAPWISE TIP DISPLACEMENT")
        print(sep)
        for m in sorted(metrics, key=lambda x: x.yaw_deg):
            if m.yaw_deg not in _YAW_REFS:
                continue
            ref_val, ref_label = _YAW_REFS[m.yaw_deg]
            err = (m.flap_mean - ref_val) / ref_val * 100
            status = "✓ PASS" if 9 <= m.flap_mean <= 15 else "✗ FAIL"
            print(
                f"  yaw={m.yaw_deg}°: sim={m.flap_mean:+.3f} m  "
                f"ref={ref_val} m ({ref_label})  err={err:+.1f}%  [{status}]"
            )
            print(f"           band: 9–15 m  |  max in ss: {m.flap_max:.3f} m")
        print(sep)


def save_csv(metrics: list[CaseMetrics], output: Path) -> None:
    field_names = [f.name for f in fields(CaseMetrics)]
    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()
        for m in sorted(metrics, key=lambda x: x.yaw_deg):
            writer.writerow({f.name: getattr(m, f.name) for f in fields(CaseMetrics)})
    print(f"\n  Saved: {output}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute steady-state metrics from a yaw-sweep FSI run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--results-dir",
        "-d",
        default="/scratch/leahk/eduardo.donestevez/frontiersin_results",
        help="Root directory containing yaw_N/ subdirectories (default: %(default)s)",
    )
    parser.add_argument(
        "--t-start",
        type=float,
        default=20.0,
        help="Lower bound of analysis window [s] (default: %(default)s)",
    )
    parser.add_argument(
        "--t-end",
        type=float,
        default=70.0,
        help="Upper bound of analysis window [s] (default: %(default)s)",
    )
    parser.add_argument(
        "--min-window",
        type=float,
        default=20.0,
        help="Minimum required time in window to include a case [s] (default: %(default)s)",
    )
    parser.add_argument(
        "--solver",
        default="corotational",
        help="Solver label used as subdirectory name (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Save summary as CSV to this path",
    )
    parser.add_argument(
        "--no-table",
        action="store_true",
        help="Suppress console table output",
    )
    parser.add_argument(
        "--cosn",
        type=float,
        default=3.0,
        help="Exponent n for cos^n power-law check (default: %(default)s)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        print(f"ERROR: results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"\n  Scanning: {results_dir}")
    print(f"  Window:   {args.t_start} < t <= {args.t_end} s  (min {args.min_window} s required)\n")

    cases = discover_cases(results_dir, args.solver)
    if not cases:
        print("No cases found.", file=sys.stderr)
        sys.exit(1)

    metrics: list[CaseMetrics] = []
    for yaw_deg, solver_label, csv_path in cases:
        all_rows = load_bem_csv(csv_path)
        rows = filter_window(all_rows, args.t_start, args.t_end)

        t_span = (rows[-1]["Time [s]"] - rows[0]["Time [s]"]) if rows else 0.0

        if not rows or t_span < args.min_window:
            total_t = all_rows[-1]["Time [s]"] if all_rows else 0.0
            print(
                f"  [SKIP] yaw={yaw_deg:>3}°: only {t_span:.1f}s in window "
                f"({args.t_start}–{min(args.t_end, total_t):.1f}s of "
                f"{total_t:.1f}s total). Need {args.min_window}s minimum.",
                file=sys.stderr,
            )
            continue

        m = compute_metrics(yaw_deg, solver_label, rows)
        metrics.append(m)
        print(
            f"  [OK]   yaw={yaw_deg:>3}°: {len(rows)} samples, "
            f"t={m.t_start_actual:.1f}–{m.t_end_actual:.1f}s"
        )

    if not metrics:
        print("\nNo cases passed the minimum-window filter.", file=sys.stderr)
        sys.exit(1)

    annotate_cosn(metrics, n=args.cosn)

    if not args.no_table:
        print_table(metrics)

    if args.output:
        save_csv(metrics, Path(args.output))


if __name__ == "__main__":
    main()
