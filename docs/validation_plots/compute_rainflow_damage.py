"""A4 — Rainflow damage analysis on flapwise root bending moment Mb(t).

Computes a pseudo-DEL (Damage Equivalent Load) characterization across the yaw
sweep using the existing ``bem_report.csv`` time series. The standalone
corotational report uses the corotational rows; optional inertial rows are kept
only for cross-solver diagnostics when those campaign files are available.

DEL definition (Palmgren-Miner with Wöhler slope m):

    DEL_m = ( Σ_i n_i · ΔS_i^m / N_eq ) ^ (1/m)

where ΔS_i are the rainflow-counted stress (here moment) ranges with counts
n_i, m is the inverse Wöhler slope, and N_eq is the equivalent number of
cycles (typically 1 Hz × T_sim for time-normalised DEL). We use m=10 (typical
for glass-fibre composite blades per IEC 61400-1 / DNV-GL guidelines for
fatigue analysis of wind turbine rotor blades) and N_eq = T_sim (so DEL has
units of MN·m at 1 Hz equivalent reference cycle).

Window: t ∈ [40, 100] s for the 100 s corotational campaign. This is short for
absolute fatigue claims, but adequate for a comparative DEL trend across yaw.

Outputs:
  docs/validation_data/generated/rainflow_damage_summary.csv
  docs/validation_data/generated/rainflow_damage_summary.md
  docs/validation_plots/figures/fig_v04_rainflow_damage.png
"""

from __future__ import annotations

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rainflow

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT_DATA = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_FIG = REPO_ROOT / "docs" / "validation_plots" / "figures"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SCRATCH = os.environ.get("SCRATCH") or str(pathlib.Path.home() / "scratch")


def _detect_corotational_root() -> pathlib.Path:
    """Pick the most recent corotacional results directory."""
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

# Default window updated to [40, 100] s for the 100 s corotacional campaign
# (avoids initial transient and gives ~7.5 revolutions of steady state).
# Override via env vars T_START_S / T_END_S if needed.
T_START = float(os.environ.get("T_START_S", "40.0"))
T_END = float(os.environ.get("T_END_S", "100.0"))
WOHLER_M = 10.0  # composite (glass-epoxy) Wöhler slope (IEC/DNV).
SIGNAL_COL = "Mb [N.m]"  # root flapwise bending moment in bem_report.csv
OMEGA_RPM = 7.518


def compute_del(values: np.ndarray, m: float, n_eq: float) -> tuple[float, list]:
    """Return DEL at 1 Hz reference rate (units: input units of values).

    Uses ASTM E1049 rainflow counting via the ``rainflow`` library, which
    yields (range, mean, count, i_start, i_end) tuples per cycle. We discard
    mean and the index info because the m=10 Wöhler model for composite is
    typically mean-insensitive after Goodman correction (omitted here for
    a comparative claim).
    """
    cycles = list(rainflow.extract_cycles(values.tolist()))
    if not cycles:
        return float("nan"), []
    ranges_counts = [(c[0], c[2]) for c in cycles]  # (range, count)
    sum_ndm = sum(count * (rng**m) for rng, count in ranges_counts)
    del_val = (sum_ndm / n_eq) ** (1.0 / m)
    return float(del_val), ranges_counts


def main() -> int:
    rows = []
    histograms = {}
    for solver_label, root in CAMPAIGNS.items():
        if not root.exists():
            print(f"[A4] skip {solver_label}: directory not found ({root})")
            continue
        for yaw in YAW_LIST:
            csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
            if not csv_path.exists():
                print(f"[A4] skip yaw={yaw} solver={solver_label}: {csv_path} missing")
                continue
            df = pd.read_csv(csv_path)
            df = df[(df["Time [s]"] >= T_START) & (df["Time [s]"] <= T_END)].reset_index(drop=True)
            if SIGNAL_COL not in df.columns:
                print(f"[A4] WARN: column {SIGNAL_COL!r} not in {csv_path}", file=sys.stderr)
                continue
            t = df["Time [s]"].to_numpy()
            mb_nm = df[SIGNAL_COL].to_numpy()
            mb_mnm = mb_nm / 1e6
            t_sim = float(t[-1] - t[0])  # equivalent number of cycles at 1 Hz reference
            mean_mb = float(np.mean(mb_mnm))
            std_mb = float(np.std(mb_mnm, ddof=1))
            range_mb = float(np.max(mb_mnm) - np.min(mb_mnm))
            del_value, ranges_counts = compute_del(mb_mnm, WOHLER_M, n_eq=t_sim)
            rows.append({
                "solver": solver_label,
                "yaw_deg": yaw,
                "n_samples": len(mb_mnm),
                "t_span_s": round(t_sim, 3),
                "Mb_mean_MNm": round(mean_mb, 4),
                "Mb_std_MNm": round(std_mb, 4),
                "Mb_pp_MNm": round(range_mb, 4),
                "n_cycles": len(ranges_counts),
                f"DEL_m{int(WOHLER_M)}_MNm": round(del_value, 4),
            })
            histograms[(solver_label, yaw)] = ranges_counts

    if not rows:
        print("[A4] no data — check SCRATCH path.")
        return 1

    df_out = pd.DataFrame(rows)
    csv_path = OUT_DATA / "rainflow_damage_summary.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"[A4] Wrote {csv_path}")

    # ── Plot: DEL vs yaw per solver + cycle range histograms at yaw=0 ───────
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    for solver in df_out["solver"].unique():
        sub = df_out[df_out["solver"] == solver].sort_values("yaw_deg")
        marker = "o" if solver == "corotational" else "s"
        color = "#c0392b" if solver == "corotational" else "#2980b9"
        ax.plot(
            sub["yaw_deg"],
            sub[f"DEL_m{int(WOHLER_M)}_MNm"],
            marker=marker,
            color=color,
            lw=1.5,
            label=solver,
            mfc="none",
            mew=1.5,
        )
    ax.set_xlabel("Yaw [°]")
    ax.set_ylabel(f"DEL_m={int(WOHLER_M)} of $M_b$ [MN·m]")
    ax.set_title(f"Damage Equivalent Load — $M_b$ root, Wöhler slope m={int(WOHLER_M)}")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)

    ax = axes[1]
    bins = np.linspace(0, 8, 21)  # MN·m range bins
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    for solver in df_out["solver"].unique():
        if (solver, 0) not in histograms:
            continue
        rc = histograms[(solver, 0)]
        if not rc:
            continue
        ranges_arr = np.array([r for r, _ in rc])
        counts_arr = np.array([c for _, c in rc])
        hist, _ = np.histogram(ranges_arr, bins=bins, weights=counts_arr)
        color = "#c0392b" if solver == "corotational" else "#2980b9"
        ls = "-" if solver == "corotational" else "--"
        ax.plot(bin_centers, hist, color=color, ls=ls, lw=1.5, label=solver, marker=".")
    ax.set_xlabel("Cycle range $\\Delta M_b$ [MN·m]")
    ax.set_ylabel("Cycle count")
    ax.set_yscale("log")
    ax.set_title("Rainflow cycle histogram — yaw=0°")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(frameon=False)
    fig.tight_layout()

    png = OUT_FIG / "fig_v04_rainflow_damage.png"
    pdf = OUT_FIG / "fig_v04_rainflow_damage.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(f"[A4] Saved {png}")

    # ── Markdown summary table ──────────────────────────────────────────────
    md = [
        "# A4 — Rainflow / DEL analysis on $M_b$ root (flapwise root bending moment)",
        "",
        f"Source: `bem_report.csv` column `{SIGNAL_COL}` (units: N·m, converted to MN·m here).",
        f"Window: t ∈ [{T_START}, {T_END}] s "
        f"(target span ≈ {(T_END - T_START) * OMEGA_RPM / 60.0:.1f} revolutions at Ω = {OMEGA_RPM:.3f} RPM).",
        f"Wöhler slope m = {int(WOHLER_M)} (typical glass-fibre composite, IEC 61400-1 / DNV-GL).",
        "DEL normalised by N_eq = T_sim, so units are MN·m at 1 Hz equivalent cycle rate.",
        "",
        "## Comparison table",
        "",
        "| Solver | Yaw [°] | n samples | span [s] | Mb mean [MN·m] | Mb std [MN·m] | Mb p-p [MN·m] | n cycles | DEL_m10 [MN·m] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in df_out.iterrows():
        md.append(
            f"| {r['solver']} | {r['yaw_deg']} | {r['n_samples']} | {r['t_span_s']:.2f} | "
            f"{r['Mb_mean_MNm']:.4f} | {r['Mb_std_MNm']:.4f} | {r['Mb_pp_MNm']:.4f} | {r['n_cycles']} | "
            f"{r[f'DEL_m{int(WOHLER_M)}_MNm']:.4f} |"
        )

    md += [
        "",
        "## Lectura",
        "",
        "- **Tendencia con yaw**: en la rama corrotacional, el DEL del Mb raíz decrece",
        "  con yaw junto con la carga media y la amplitud cíclica absoluta. La métrica",
        "  captura una tendencia comparativa regular, sin picos aislados que dominen el daño.",
        "- **Alcance de las filas adicionales**: si existen filas de otro solver, deben leerse",
        "  solo como diagnóstico auxiliar y con atención a `span [s]`; el reporte corrotacional",
        "  usa las filas `corotational` de la campaña de 100 s.",
        "- **Limitación de tiempo simulado**: 60 s ≈ 7.5 revoluciones es corto para una",
        "  reclamación de fatiga absoluta. La DEL aquí es comparativa, no absoluta. Para",
        "  certificación IEC 61400-1 se requiere DLC 1.1 (NTM, 600 s × 6 semillas) que",
        "  excede el alcance de este dossier.",
        "- **Publicabilidad**: la combinación Rainflow + DEL sobre el barrido yaw de la",
        "  IEA 15 MW con BEM+shell es una caracterización dinámica original, pero todavía",
        "  no reemplaza una campaña DLC turbulenta.",
    ]
    md_path = OUT_DATA / "rainflow_damage_summary.md"
    md_path.write_text("\n".join(md) + "\n")
    print(f"[A4] Wrote {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
