#!/usr/bin/env python3
"""Campaign metrics: collect them, and compare a re-run against a recorded snapshot.

The FSI campaigns are compared against a "before" column that lived only in the
campaign directories on $SCRATCH.  `docs/validation_data/campaign_metrics_before_2026-09-30.csv`
is that column, taken before those directories were removed to free space, so the
re-run can be diffed against it.

Collect mode mirrors the recorded CSV's conventions: the settled window is the last
60% of each run, flap is `Tip Disp Y`, edge is `Tip Disp X`.

    python tools/campaign_metrics.py collect  <results-dir> [<results-dir> ...]
    python tools/campaign_metrics.py compare  <new.csv> [--baseline PATH]

A directory is searched for `<case>/fluid/bem_report.csv` and
`<case>/<sub>/bem_report.csv`, so it works for both the flat layout (yaw, parked,
ch6) and the nested one (the V-06 matrix).
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

WINDOW_FRACTION = 0.6  # the recorded runs used t in [40, 100] on a 100 s run
DEFAULT_BASELINE = (
    Path(__file__).resolve().parents[1]
    / "docs"
    / "validation_data"
    / "campaign_metrics_before_2026-09-30.csv"
)
METRICS = ("flap_mean", "flap_std", "edge_mean", "thrust_mean", "power_mean_mw", "ct_mean", "cp_mean")


def summarise(path: Path) -> dict | None:
    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        return None
    t = np.array([float(r["Time [s]"]) for r in rows])
    t0, t1 = float(t.min()), float(t.max())
    t_start = t0 + (t1 - t0) * (1.0 - WINDOW_FRACTION)
    mask = t >= t_start

    def col(name):
        return np.array([float(r[name]) for r in rows])[mask]

    flap, edge = col("Tip Disp Y [m]"), col("Tip Disp X [m]")
    thrust, torque, power = col("Thrust [N]"), col("Torque [N.m]"), col("Power [W]")
    return {
        "n_samples": int(mask.sum()),
        "t_start_actual": round(t_start, 2),
        "t_end_actual": round(t1, 2),
        "flap_mean": float(flap.mean()),
        "flap_std": float(flap.std()),
        "flap_max": float(flap.max()),
        "flap_min": float(flap.min()),
        "flap_p2p": float(flap.max() - flap.min()),
        "edge_mean": float(edge.mean()),
        "edge_std": float(edge.std()),
        "thrust_mean": float(thrust.mean()),
        "thrust_std": float(thrust.std()),
        "torque_mean": float(torque.mean()),
        "torque_std": float(torque.std()),
        "power_mean_mw": float(power.mean()) / 1e6,
        "power_std_mw": float(power.std()) / 1e6,
        "ct_mean": float(col("CT").mean()),
        "cq_mean": float(col("CQ").mean()),
        "cp_mean": float(col("CP").mean()),
        "diverged": bool(not np.isfinite(flap).all() or np.abs(flap).max() > 100.0),
    }


def collect(dirs: list[Path]) -> list[dict]:
    out, seen = [], set()
    for base in dirs:
        if not base.exists():
            print(f"(no such directory: {base})", file=sys.stderr)
            continue
        reports = sorted(base.glob("*/fluid/bem_report.csv")) + sorted(base.glob("*/*/bem_report.csv"))
        for report in reports:
            if report in seen:
                continue
            seen.add(report)
            summary = summarise(report)
            if summary is None:
                continue
            case = report.parent.parent.name
            out.append({"campaign": base.name, "case": case, **summary})
    return out


def compare(new_csv: Path, baseline: Path) -> int:
    new = {f"{r['campaign']}/{r['case']}": r for r in csv.DictReader(new_csv.open(newline=""))}
    old = {f"{r['campaign']}/{r['case']}": r for r in csv.DictReader(baseline.open(newline=""))}
    if not new:
        print(f"no rows in {new_csv}", file=sys.stderr)
        return 1
    print(f"{'case':>26} {'metric':>14} {'before':>12} {'after':>12} {'delta':>9}")
    matched = 0
    for key in sorted(new):
        if key not in old:
            continue
        matched += 1
        for metric in METRICS:
            try:
                a, b = float(old[key][metric]), float(new[key][metric])
            except (KeyError, TypeError, ValueError):
                continue
            if b == 0.0:
                continue
            print(f"{key:>26} {metric:>14} {a:12.4f} {b:12.4f} {100 * (b - a) / b:+8.2f}%")
    print(f"\nmatched {matched} of {len(new)} cases against {baseline.name}")
    if matched < len(new):
        only = sorted(set(new) - set(old))
        print(f"only in the re-run: {only}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("collect", help="summarise every campaign case under these dirs")
    c.add_argument("dirs", nargs="+", type=Path)
    c.add_argument("--csv", type=Path, default=Path("campaign_metrics.csv"))
    cmp_ = sub.add_parser("compare", help="diff a re-run against the recorded snapshot")
    cmp_.add_argument("csv", type=Path)
    cmp_.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    args = ap.parse_args()

    if args.mode == "collect":
        rows = collect(args.dirs)
        if not rows:
            print("nothing collected", file=sys.stderr)
            return 1
        with args.csv.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows)} cases -> {args.csv}")
        for row in rows:
            print(f"  {row['campaign']}/{row['case']}: flap {row['flap_mean']:.3f} m "
                  f"thrust {row['thrust_mean'] / 1e6:.3f} MN power {row['power_mean_mw']:.3f} MW")
        return 0
    return compare(args.csv, args.baseline)


if __name__ == "__main__":
    sys.exit(main())
