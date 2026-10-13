"""Contraction and whole-period statistics of a coupled BEM-FSI smoke run.

Usage: python tools/coupled_run_stats.py <run_dir> [<run_dir> ...] [--period 7.947] [--periods 2]

For each ``<run_dir>`` (a ``runs/smoke_fix/<RUN_NAME>`` folder) it reads

* ``precice-Solid-iterations.log``: windows, how many hit ``--max-iterations``
  (default 30, the smoke config's limit), mean and maximum sub-iterations;
* ``fluid/bem_report.csv``: means over the last ``--periods`` whole rotor periods
  ``(t_end - n T, t_end]``, never a fixed window that cuts a period.  ``--t-end``
  truncates every run at the same time, so a shorter run and a longer one are compared
  over the same span (the coupled solution is causal: truncating a run at ``t`` gives
  the run that stopped at ``t``).

The coupled solution is deterministic, so a difference between two runs is a model
difference, not noise; the tolerance that decides what a difference means belongs to
the caller.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

COLUMNS = {
    "thrust_MN": ("Thrust [N]", 1e-6),
    "power_MW": ("Power [W]", 1e-6),
    "Mb_MNm": ("Mb [N.m]", 1e-6),
    "tip_flap_Y_m": ("Tip Disp Y [m]", 1.0),
    "tip_edge_X_m": ("Tip Disp X [m]", 1.0),
}


def contraction(run: Path, max_iterations: int, t_end: float | None, dt: float) -> dict:
    rows = np.loadtxt(run / "precice-Solid-iterations.log", skiprows=1, ndmin=2)
    if t_end is not None:
        rows = rows[rows[:, 0] <= round(t_end / dt)]
    iterations = rows[:, 2]
    return {
        "windows": int(len(rows)),
        "at_max_iterations": int(np.sum(iterations >= max_iterations)),
        "mean_iterations": float(iterations.mean()),
        "max_iterations": int(iterations.max()),
    }


def whole_period_means(run: Path, period: float, n_periods: int, t_end: float | None) -> dict:
    with open(run / "fluid" / "bem_report.csv", newline="") as fh:
        reader = csv.DictReader(fh)
        data = {k: [] for k in ["Time [s]", *[c for c, _ in COLUMNS.values()]]}
        for row in reader:
            for k in data:
                data[k].append(float(row[k]))
    t = np.asarray(data["Time [s]"])
    if t_end is None:
        t_end = float(t[-1])
    elif t_end > t[-1] + 1e-9:
        raise ValueError(f"{run}: data end at {t[-1]}, before --t-end {t_end}")
    t0 = t_end - n_periods * period
    if t0 < t[0]:
        raise ValueError(f"{run}: {n_periods} periods need t >= {t0:.3f}, data start at {t[0]}")
    sel = (t > t0) & (t <= t_end + 1e-9)  # (t0, t_end]: n whole periods
    out = {"window_s": f"({t0:.3f}, {t_end:.3f}]", "samples": int(sel.sum())}
    for name, (col, scale) in COLUMNS.items():
        out[name] = float(np.mean(np.asarray(data[col])[sel])) * scale
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n")[0])
    ap.add_argument("runs", nargs="+", type=Path)
    ap.add_argument("--period", type=float, default=7.947)
    ap.add_argument("--periods", type=int, default=2)
    ap.add_argument("--max-iterations", type=int, default=30)
    ap.add_argument("--t-end", type=float, default=None)
    ap.add_argument("--dt", type=float, default=1e-2, help="preCICE time-window size")
    args = ap.parse_args()
    for run in args.runs:
        print(f"== {run}")
        for k, v in contraction(run, args.max_iterations, args.t_end, args.dt).items():
            print(f"  {k:20s} {v}")
        for k, v in whole_period_means(run, args.period, args.periods, args.t_end).items():
            print(f"  {k:20s} {v if isinstance(v, (str, int)) else f'{v:.4f}'}")


if __name__ == "__main__":
    main()
