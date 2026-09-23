"""B1+B2 post — h/dt convergence plot for V-04 corotacional FSI.

Reads bem_report.csv from each convergence case + the baseline (frontiersin
yaw_0 corotacional) and produces:

  - Convergence plot: flap_mean, flap_std, flap_pp, power_mean vs h (B1) and vs dt (B2)
  - Richardson extrapolation estimate of asymptotic value and effective order
  - CSV summary with all metrics per case

Run AFTER `runConvergence.srm` completes (or as cases finish — partial plots
are produced if some cases are missing).

Outputs:
  docs/validation_data/generated/convergence_b1_b2_summary.csv
  docs/validation_plots/figures/fig_v04_convergence_b1_b2.png

Env vars:
  RESULTS_BASE   convergence results base (default: $SCRATCH/convergence_b1_b2_results_b1fix)
  BASELINE_PATH  baseline results path (default: $SCRATCH/frontiersin_results_corotational_100s/yaw_0)
"""

from __future__ import annotations

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT_DATA = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_FIG = REPO_ROOT / "docs" / "validation_plots" / "figures"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SCRATCH = os.environ.get("SCRATCH") or str(pathlib.Path.home() / "scratch")
RESULTS_BASE = pathlib.Path(os.environ.get("RESULTS_BASE")
                            or f"{SCRATCH}/convergence_b1_b2_results_b1fix")
BASELINE_PATH = pathlib.Path(os.environ.get("BASELINE_PATH")
                             or f"{SCRATCH}/frontiersin_results_corotational_100s/yaw_0")

T_STEADY_START = 20.0   # baseline ran to 70s; convergence runs to 30s — use [20,30] common window
T_STEADY_END = 30.0

# Case definitions: (tag, family, h, dt, csv_path)
CASES = [
    ("h_coarse",    "B1", 0.500, 0.01,  RESULTS_BASE / "h_coarse"  / "fluid" / "bem_report.csv"),
    ("h_baseline",  "B1", 0.250, 0.01,  BASELINE_PATH               / "fluid" / "bem_report.csv"),
    ("h_fine",      "B1", 0.125, 0.01,  RESULTS_BASE / "h_fine"    / "fluid" / "bem_report.csv"),
    ("dt_coarse",   "B2", 0.250, 0.02,  RESULTS_BASE / "dt_coarse" / "fluid" / "bem_report.csv"),
    ("dt_baseline", "B2", 0.250, 0.01,  BASELINE_PATH               / "fluid" / "bem_report.csv"),
    ("dt_fine",     "B2", 0.250, 0.005, RESULTS_BASE / "dt_fine"   / "fluid" / "bem_report.csv"),
]


def load_steady_metrics(csv_path: pathlib.Path) -> dict | None:
    if not csv_path.exists():
        return None
    df = pd.read_csv(csv_path)
    df = df[(df["Time [s]"] >= T_STEADY_START) & (df["Time [s]"] <= T_STEADY_END)].reset_index(drop=True)
    if df.empty:
        return None
    flap = df["Tip Disp Y [m]"].to_numpy()
    edge = df["Tip Disp X [m]"].to_numpy()
    power_mw = df["Power [W]"].to_numpy() / 1e6
    torque_mnm = df["Torque [N.m]"].to_numpy() / 1e6
    thrust_mn = df["Thrust [N]"].to_numpy() / 1e6
    return {
        "n_samples": len(df),
        "t_span_s": round(float(df["Time [s]"].iloc[-1] - df["Time [s]"].iloc[0]), 3),
        "flap_mean":  round(float(np.mean(flap)),  4),
        "flap_std":   round(float(np.std(flap, ddof=1)),  4),
        "flap_pp":    round(float(np.max(flap) - np.min(flap)), 4),
        "edge_mean":  round(float(np.mean(edge)),  4),
        "power_mean": round(float(np.mean(power_mw)), 4),
        "torque_mean": round(float(np.mean(torque_mnm)), 4),
        "thrust_mean": round(float(np.mean(thrust_mn)), 4),
    }


def richardson_extrapolate(coarse: float, baseline: float, fine: float,
                           r: float = 2.0) -> tuple[float, float] | None:
    """Estimate asymptotic value f_exact and effective order p.

    Uses three values at refinement levels h_c = r·h_b, h_b, h_f = h_b / r.
    f_h = f_exact + C h^p
    p = log[(f_coarse - f_baseline) / (f_baseline - f_fine)] / log(r)
    f_exact ≈ f_fine + (f_fine - f_baseline) / (r^p - 1)
    Returns None if the differences have the wrong sign for Richardson to apply.
    """
    num = coarse - baseline
    den = baseline - fine
    if abs(den) < 1e-12 or (num / den) <= 0:
        return None
    p = np.log(num / den) / np.log(r)
    denom = r ** p - 1.0
    if abs(p) < 1e-9 or abs(denom) < 1e-12:
        return None
    f_exact = fine + (fine - baseline) / denom
    return float(f_exact), float(p)


def main() -> int:
    rows = []
    metrics_by_tag = {}
    for tag, family, h, dt, csv_path in CASES:
        m = load_steady_metrics(csv_path)
        if m is None:
            print(f"[B1B2] WARN: missing or empty {csv_path}", file=sys.stderr)
            continue
        m["tag"] = tag
        m["family"] = family
        m["element_size_m"] = h
        m["time_window_size_s"] = dt
        rows.append(m)
        metrics_by_tag[tag] = m

    if not rows:
        print("[B1B2] no data — none of the cases completed yet.")
        return 1

    df = pd.DataFrame(rows)
    csv_out = OUT_DATA / "convergence_b1_b2_summary.csv"
    df.to_csv(csv_out, index=False)
    print(f"[B1B2] Wrote {csv_out}")

    # ── Richardson extrapolation (only if all 3 levels present) ───────────────
    summary_lines = ["# B1+B2 — Convergence summary (auto-generated)\n"]
    for family, triplet, x_label in [
        ("B1", ["h_coarse", "h_baseline", "h_fine"], "element_size"),
        ("B2", ["dt_coarse", "dt_baseline", "dt_fine"], "time-window-size"),
    ]:
        summary_lines.append(f"## {family} — {x_label}\n")
        if not all(t in metrics_by_tag for t in triplet):
            missing = [t for t in triplet if t not in metrics_by_tag]
            summary_lines.append(f"Missing: {missing}. Skip Richardson.\n")
            continue
        summary_lines.append("| Metric | Coarse | Baseline | Fine | Richardson f∞ | Order p |")
        summary_lines.append("|---|---:|---:|---:|---:|---:|")
        for metric in ["flap_mean", "flap_std", "flap_pp", "power_mean", "torque_mean"]:
            c = metrics_by_tag[triplet[0]][metric]
            b = metrics_by_tag[triplet[1]][metric]
            f = metrics_by_tag[triplet[2]][metric]
            re = richardson_extrapolate(c, b, f, r=2.0)
            if re is None:
                summary_lines.append(f"| {metric} | {c} | {b} | {f} | n/a | n/a |")
            else:
                f_exact, p = re
                summary_lines.append(f"| {metric} | {c} | {b} | {f} | {f_exact:.4f} | {p:.2f} |")
        summary_lines.append("")

    md_out = OUT_DATA / "convergence_b1_b2_summary.md"
    md_out.write_text("\n".join(summary_lines) + "\n")
    print(f"[B1B2] Wrote {md_out}")

    # ── Plot 2×2: flap_mean and power_mean vs h (left) and vs dt (right) ──────
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=False)
    metric_pairs = [("flap_mean", "Flap mean [m]"), ("power_mean", "Power mean [MW]")]

    for row, (metric, ylabel) in enumerate(metric_pairs):
        # B1: vs element_size (log x)
        ax = axes[row, 0]
        b1 = df[df["family"] == "B1"].sort_values("element_size_m")
        if not b1.empty:
            ax.plot(b1["element_size_m"], b1[metric], "o-", color="#c0392b", mfc="none", mew=1.5)
            for _, r in b1.iterrows():
                ax.annotate(r["tag"], xy=(r["element_size_m"], r[metric]),
                            xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("element_size [m]")
        ax.set_ylabel(ylabel)
        ax.set_title(f"B1 — h-convergence ({metric})")
        ax.grid(True, alpha=0.3, which="both")

        # B2: vs time-window-size (log x)
        ax = axes[row, 1]
        b2 = df[df["family"] == "B2"].sort_values("time_window_size_s")
        if not b2.empty:
            ax.plot(b2["time_window_size_s"], b2[metric], "s-", color="#2980b9", mfc="none", mew=1.5)
            for _, r in b2.iterrows():
                ax.annotate(r["tag"], xy=(r["time_window_size_s"], r[metric]),
                            xytext=(5, 5), textcoords="offset points", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("time-window-size [s]")
        ax.set_ylabel(ylabel)
        ax.set_title(f"B2 — dt-convergence ({metric})")
        ax.grid(True, alpha=0.3, which="both")

    fig.suptitle("B1+B2 — Convergence study, V-04 corotacional, yaw=0°, V=10.59 m/s, Ω=7.55 RPM",
                 fontsize=11)
    fig.tight_layout()
    png = OUT_FIG / "fig_v04_convergence_b1_b2.png"
    pdf = OUT_FIG / "fig_v04_convergence_b1_b2.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(f"[B1B2] Saved {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
