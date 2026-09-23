"""A3 — Non-overlapping window bootstrap of FSI rotor signals.

Quantifies statistical stability over the definitive steady-state window used
by the corotational report. The default ``40 <= t <= 100`` s window gives six
non-overlapping 10 s sub-windows for the 100 s corotational campaign, removing
the residual transient contamination present in the previous 20-70 s cut.

Procedure (per solver, per yaw, per column):
  1. Read bem_report.csv, filter t ∈ [t_start, t_end].
  2. Split into K non-overlapping sub-windows of length T_w seconds.
  3. Within each sub-window compute mean, std, peak-to-peak.
  4. Across the K sub-windows compute mean of means and std of means,
     report IC95 = ±1.96 · std / sqrt(K).
  5. Optionally include FFT amplitude at characteristic frequencies if a
     constant-rate VTU/tip column is provided (out of scope here — keeps
     this script as a pure bem_report.csv reader for portability).

Output:
  docs/validation_data/generated/statistical_stability_bem_report.csv

Run from anywhere; uses SCRATCH env var to locate the campaign directories.
"""

from __future__ import annotations

import csv
import os
import pathlib
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SCRATCH = os.environ.get("SCRATCH") or str(pathlib.Path.home() / "scratch")


def _detect_corotational_root() -> pathlib.Path:
    override = os.environ.get("COROT_ROOT")
    if override:
        return pathlib.Path(override)
    candidates = [
        pathlib.Path(SCRATCH) / "frontiersin_results_corotational_100s",
        pathlib.Path(SCRATCH) / "frontiersin_results_corotational",
    ]
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


CAMPAIGNS = {
    "corotational": _detect_corotational_root(),
    "inertial": pathlib.Path(SCRATCH) / "frontiersin_results_inertial",
}

YAW_LIST = [0, 10, 20, 30, 40]
# Default window: [40, 100] s for the 100 s corotational campaign.
T_START = float(os.environ.get("T_START_S", "40.0"))
T_END = float(os.environ.get("T_END_S", "100.0"))
T_WINDOW = float(os.environ.get("T_WINDOW_S", "10.0"))

# Columns to evaluate. (csv_name, short_label, unit, statistic).
TARGETS = [
    ("Tip Disp Y [m]", "flap", "m", "mean"),
    ("Tip Disp Y [m]", "flap", "m", "std"),
    ("Tip Disp Y [m]", "flap", "m", "pp"),
    ("Tip Disp X [m]", "edge", "m", "mean"),
    ("Tip Disp X [m]", "edge", "m", "std"),
    ("Power [W]", "power", "MW", "mean"),
    ("Power [W]", "power", "MW", "pp"),
    ("Torque [N.m]", "torque", "MN·m", "mean"),
    ("Thrust [N]", "thrust", "MN", "mean"),
]


def load_signal(path: pathlib.Path) -> pd.DataFrame | None:
    if not path.exists():
        print(f"[A3] WARN: missing {path}", file=sys.stderr)
        return None
    df = pd.read_csv(path)
    return df[(df["Time [s]"] >= T_START) & (df["Time [s]"] <= T_END)].reset_index(drop=True)


def scale_for_label(unit: str) -> float:
    if unit == "MW":
        return 1e-6
    if unit == "MN":
        return 1e-6
    if unit == "MN·m":
        return 1e-6
    return 1.0


def stat_over_window(values: np.ndarray, kind: str) -> float:
    if values.size == 0:
        return float("nan")
    if kind == "mean":
        return float(np.mean(values))
    if kind == "std":
        return float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    if kind == "pp":
        return float(np.max(values) - np.min(values))
    raise ValueError(f"unknown kind {kind}")


def bootstrap_subwindows(df: pd.DataFrame, column: str, kind: str, scale: float) -> dict:
    """Compute per-subwindow statistic, return mean ± IC95 across sub-windows."""
    t = df["Time [s]"].to_numpy()
    y = df[column].to_numpy() * scale
    edges = np.arange(T_START, T_END + 1e-9, T_WINDOW)
    per_window = []
    for k in range(len(edges) - 1):
        mask = (t >= edges[k]) & (t < edges[k + 1])
        per_window.append(stat_over_window(y[mask], kind))
    per_window = np.asarray(per_window, dtype=float)
    k_eff = int(np.sum(np.isfinite(per_window)))
    if k_eff == 0:
        return dict(mean=float("nan"), std=float("nan"), ic95=float("nan"), k=0)
    mean_val = float(np.nanmean(per_window))
    std_val = float(np.nanstd(per_window, ddof=1)) if k_eff > 1 else 0.0
    ic95 = 1.96 * std_val / np.sqrt(k_eff) if k_eff > 1 else 0.0
    return dict(mean=mean_val, std=std_val, ic95=ic95, k=k_eff, per_window=per_window)


def main() -> int:
    out_rows = []
    for solver_label, root in CAMPAIGNS.items():
        if not root.exists():
            print(f"[A3] skip {solver_label}: directory not found ({root})")
            continue
        for yaw in YAW_LIST:
            csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
            df = load_signal(csv_path)
            if df is None or df.empty:
                continue
            for column, short, unit, kind in TARGETS:
                if column not in df.columns:
                    continue
                scale = scale_for_label(unit)
                res = bootstrap_subwindows(df, column, kind, scale)
                if not np.isfinite(res["mean"]):
                    continue
                out_rows.append({
                    "solver": solver_label,
                    "yaw_deg": yaw,
                    "signal": short,
                    "stat": kind,
                    "unit": unit,
                    "mean": round(res["mean"], 6),
                    "std": round(res["std"], 6),
                    "ic95": round(res["ic95"], 6),
                    "k_subwindows": res["k"],
                    "t_window_s": T_WINDOW,
                })

    if not out_rows:
        print("[A3] No data — check SCRATCH path and campaign directories.")
        return 1

    csv_out = OUT_DIR / "statistical_stability_bem_report.csv"
    with csv_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"[A3] Wrote {csv_out} ({len(out_rows)} rows)")

    # Side-by-side significance table for yaw=0° (main steady-window evidence).
    by_key = defaultdict(dict)
    for r in out_rows:
        if r["yaw_deg"] != 0:
            continue
        key = (r["signal"], r["stat"], r["unit"])
        by_key[key][r["solver"]] = r

    corot_k = next(
        (
            r["k_subwindows"]
            for r in out_rows
            if r["solver"] == "corotational" and r["yaw_deg"] == 0
        ),
        "n/a",
    )
    md_lines = [
        f"# A3 — Statistical stability (yaw=0°, T_w={T_WINDOW:g}s, corotational K={corot_k})",
        "",
        f"Window: [{T_START}, {T_END}] s, sub-window length = {T_WINDOW} s.",
        "Reports mean ± IC95 across non-overlapping sub-windows. The corotational",
        "rows are the evidence used by the standalone corotational report; any",
        "additional solver rows are kept only for cross-solver diagnostics.",
        "",
        "**Caveat para `std` y `pp`**: estas columnas reportan IC95 de la dispersión",
        "*dentro de una ventana de 10 s*. Eso es lo correcto si interesa saber si la",
        "amplitud típica observada en un fragmento de 10 s es robusta. La `std` global",
        "sobre la ventana completa es un objeto distinto: no se mide bien por",
        "promedio de stds locales, porque el promedio de stds locales subestima la std",
        "global cuando la señal contiene escalas de período mayor a 10 s. La columna",
        "`mean` no sufre este sesgo.",
        "",
        "| Signal | Stat | Unit | Corot mean ± IC95 | Inert mean ± IC95 | Δ (inert-corot) | Overlap? |",
        "|---|---|---|---|---|---:|:---:|",
    ]
    for (signal, kind, unit), per_solver in sorted(by_key.items()):
        c = per_solver.get("corotational")
        i = per_solver.get("inertial")
        if c is None or i is None:
            continue
        cell_c = f"{c['mean']:.4f} ± {c['ic95']:.4f}"
        cell_i = f"{i['mean']:.4f} ± {i['ic95']:.4f}"
        delta = i["mean"] - c["mean"]
        # Overlap test: bands overlap iff |Δ| <= IC95_c + IC95_i.
        overlap = abs(delta) <= (c["ic95"] + i["ic95"])
        flag = "sí" if overlap else "**no**"
        md_lines.append(
            f"| {signal} | {kind} | {unit} | {cell_c} | {cell_i} | {delta:+.4f} | {flag} |"
        )
    md_out = OUT_DIR / "statistical_stability_summary.md"
    md_out.write_text("\n".join(md_lines) + "\n")
    print(f"[A3] Wrote {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
