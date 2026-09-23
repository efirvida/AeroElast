"""Sub-band spectral variance of rotor torque.

Quantifies where rotor-torque variance lives in frequency for the yaw sweep.
The standalone corotational report uses the corotational rows to show that the
torque signal remains concentrated in physically interpretable low-frequency
bands; optional inertial rows are retained only for cross-solver diagnostics.

Método: Parseval por bandas. Para cada solver/yaw se lee la señal de torque
Q(t) = `Torque [N.m]` en bem_report.csv sobre la ventana configurada, se centra,
se aplica ventana Hann con normalización Parseval-coherente, y se computa la PSD
unilateral. La varianza total se reparte en tres bandas físicamente significativas:

  B_low  : 0.05 ≤ f ≤ 1.0 Hz   — contiene 1P (0.125), 3P (0.376), f1_flap (0.554), f1_edge (0.629)
  B_mid  : 1.0  <  f ≤ 5.0 Hz  — contiene 2do flap (1.69), 2do edge (1.98), 1er torsional (4.54)
  B_high : 5.0  <  f ≤ 25 Hz   — contenido sin contraparte modal estructural

El reparto se reporta como σ_Q por banda (= sqrt(var_banda)) en kN·m y como
fracción de la varianza total por banda.

Sanity check Parseval: var_total ≈ Σ var_bins (verificado al final).

Outputs:
  docs/validation_data/generated/torque_subband_variance.csv
  docs/validation_data/generated/torque_subband_variance_summary.md
  docs/validation_plots/figures/fig_v04_torque_subband_variance.png
"""

from __future__ import annotations

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import windows

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT_DATA = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_FIG = REPO_ROOT / "docs" / "validation_plots" / "figures"
OUT_DATA.mkdir(parents=True, exist_ok=True)
OUT_FIG.mkdir(parents=True, exist_ok=True)

SCRATCH = os.environ.get("SCRATCH") or str(pathlib.Path.home() / "scratch")


def _detect_corotational_root() -> pathlib.Path:
    """Pick the most recent corotacional results directory."""
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

# Default window: [40, 100] s for the 100 s corotacional campaign.
T_START = float(os.environ.get("T_START_S", "40.0"))
T_END = float(os.environ.get("T_END_S", "100.0"))

# Frequency bands (Hz). Lower bound 0.05 Hz excludes DC/trend; upper 25 Hz is well
# above any modal frequency of interest. Edges chosen by structural meaning:
#   1.0 Hz separates fundamental flap/edge (~0.55-0.63) from higher modes.
#   5.0 Hz sits just above the 1st torsional (4.54 Hz).
BANDS = {
    "B_low": (0.05, 1.0),
    "B_mid": (1.0, 5.0),
    "B_high": (5.0, 25.0),
}

SIGNAL_COL = "Torque [N.m]"


def parseval_psd_bands(signal: np.ndarray, dt: float, bands: dict) -> tuple[dict, float, float]:
    """Compute one-sided PSD with Hann window (Parseval-coherent) and return
    variance per band (in (input units)^2) plus the total variance from the
    time series for a sanity check.

    Returns: (var_by_band, var_time_total, var_freq_total)
    """
    sig = np.asarray(signal, dtype=float)
    sig = sig - float(np.mean(sig))
    n = len(sig)
    win = windows.hann(n, sym=False)
    cw = float(np.sum(win**2) / n)  # power-correction factor for Hann
    x = sig * win
    spec = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=dt)
    # Variance per bin via Parseval (one-sided, window-corrected):
    #   Σ var_per_bin == var(x_centered) when properly normalised.
    var_per_bin = (np.abs(spec) ** 2) / (n**2 * cw)
    if n % 2 == 0:
        var_per_bin[1:-1] *= 2.0  # interior bins doubled for one-sided
    else:
        var_per_bin[1:] *= 2.0
    out = {}
    for label, (fmin, fmax) in bands.items():
        mask = (freqs >= fmin) & (freqs < fmax)
        out[label] = float(np.sum(var_per_bin[mask]))
    var_time_total = float(np.var(sig, ddof=1))
    var_freq_total = float(np.sum(var_per_bin))
    return out, var_time_total, var_freq_total


def main() -> int:
    rows = []
    for solver_label, root in CAMPAIGNS.items():
        if not root.exists():
            print(f"[SB] skip {solver_label}: missing {root}")
            continue
        for yaw in YAW_LIST:
            csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
            if not csv_path.exists():
                print(f"[SB] skip yaw={yaw} solver={solver_label}: {csv_path} missing")
                continue
            df = pd.read_csv(csv_path)
            df = df[(df["Time [s]"] >= T_START) & (df["Time [s]"] <= T_END)].reset_index(drop=True)
            if df.empty or SIGNAL_COL not in df.columns:
                continue
            t = df["Time [s]"].to_numpy()
            q_nm = df[SIGNAL_COL].to_numpy()
            q_knm = q_nm / 1e3
            dt = float(np.median(np.diff(t)))
            var_bands, var_time, var_freq = parseval_psd_bands(q_knm, dt, BANDS)
            sigma_bands = {k: float(np.sqrt(v)) for k, v in var_bands.items()}
            # Fraction of total variance per band
            sum_var = sum(var_bands.values())
            frac = {k: (v / sum_var if sum_var > 0 else float("nan")) for k, v in var_bands.items()}
            rows.append({
                "solver": solver_label,
                "yaw_deg": yaw,
                "dt_s": round(dt, 6),
                "n_samples": len(q_knm),
                "sigma_Q_total_kNm": round(float(np.sqrt(var_freq)), 4),
                "var_time_check": round(var_time, 4),
                "var_freq_check": round(var_freq, 4),
                "parseval_ratio": round(var_freq / var_time, 4) if var_time > 0 else float("nan"),
                **{f"sigma_{k}_kNm": round(sigma_bands[k], 4) for k in BANDS},
                **{f"frac_{k}": round(frac[k], 4) for k in BANDS},
            })

    if not rows:
        print("[SB] no data — check SCRATCH path.")
        return 1

    df_out = pd.DataFrame(rows)
    csv_path = OUT_DATA / "torque_subband_variance.csv"
    df_out.to_csv(csv_path, index=False)
    print(f"[SB] Wrote {csv_path}")

    # Sanity-check Parseval (var_freq / var_time should be ~1).
    par_min = df_out["parseval_ratio"].min()
    par_max = df_out["parseval_ratio"].max()
    print(
        f"[SB] Parseval check: ratio ∈ [{par_min:.3f}, {par_max:.3f}] (expected ~1; "
        f"slight deviation OK due to Hann + finite Δf)"
    )

    # ── Plot: two panels ─────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))

    # Panel A: σ per band, all yaws, both solvers (grouped bars)
    ax = axes[0]
    band_labels = list(BANDS.keys())
    n_bands = len(band_labels)
    n_yaws = len(YAW_LIST)
    bar_w = 0.35
    x_idx = np.arange(n_bands)
    # one set of grouped bars per yaw is too cluttered; show only yaw=0 here
    yaw_focus = 0
    corot_row = df_out[(df_out["solver"] == "corotational") & (df_out["yaw_deg"] == yaw_focus)]
    inert_row = df_out[(df_out["solver"] == "inertial") & (df_out["yaw_deg"] == yaw_focus)]
    if not corot_row.empty and not inert_row.empty:
        corot_vals = [corot_row[f"sigma_{b}_kNm"].iloc[0] for b in band_labels]
        inert_vals = [inert_row[f"sigma_{b}_kNm"].iloc[0] for b in band_labels]
        ax.bar(x_idx - bar_w / 2, corot_vals, bar_w, label="corotational", color="#c0392b")
        ax.bar(x_idx + bar_w / 2, inert_vals, bar_w, label="inertial", color="#2980b9")
        for i, (c, j) in enumerate(zip(corot_vals, inert_vals)):
            ratio = j / c if c > 0 else float("nan")
            ax.text(
                x_idx[i],
                max(c, j) * 1.05,
                f"×{ratio:.2f}",
                ha="center",
                fontsize=9,
                color="black",
                fontweight="bold",
            )
    ax.set_xticks(x_idx)
    ax.set_xticklabels([f"{b}\n{BANDS[b][0]:g}–{BANDS[b][1]:g} Hz" for b in band_labels])
    ax.set_ylabel("σ_Q per band [kN·m]")
    ax.set_title(f"Torque sub-band variance — yaw={yaw_focus}°")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, which="both", axis="y")
    ax.legend(frameon=False, loc="upper right")

    # Panel B: ratio inert/corot per band, across yaws
    ax = axes[1]
    for b, color in zip(band_labels, ["#27ae60", "#e67e22", "#8e44ad"]):
        ratios = []
        yaws_present = []
        for yaw in YAW_LIST:
            c = df_out[(df_out["solver"] == "corotational") & (df_out["yaw_deg"] == yaw)]
            i = df_out[(df_out["solver"] == "inertial") & (df_out["yaw_deg"] == yaw)]
            if c.empty or i.empty:
                continue
            cv = float(c[f"sigma_{b}_kNm"].iloc[0])
            iv = float(i[f"sigma_{b}_kNm"].iloc[0])
            if cv <= 0:
                continue
            ratios.append(iv / cv)
            yaws_present.append(yaw)
        if ratios:
            ax.plot(
                yaws_present,
                ratios,
                "o-",
                color=color,
                lw=1.5,
                mfc="none",
                mew=1.5,
                label=f"{b} ({BANDS[b][0]:g}–{BANDS[b][1]:g} Hz)",
            )
    ax.axhline(1.0, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("Yaw [°]")
    ax.set_ylabel("σ ratio: inertial / corotational")
    ax.set_title("Ratio inercial/corot por banda y por yaw")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()

    png = OUT_FIG / "fig_v04_torque_subband_variance.png"
    pdf = OUT_FIG / "fig_v04_torque_subband_variance.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(f"[SB] Saved {png}")

    # ── Markdown summary ────────────────────────────────────────────────────
    md = [
        "# Sub-band spectral variance of rotor torque",
        "",
        f"Source: `bem_report.csv` column `{SIGNAL_COL}` (units: N·m, converted to kN·m).",
        f"Window: t ∈ [{T_START}, {T_END}] s. Hann window, Parseval-coherent PSD.",
        "The standalone corotational report uses the corotational rows; optional inertial",
        "rows are retained only for cross-solver diagnostics when available.",
        "",
        "## Bands",
        "",
        "| Label | Range [Hz] | Physical content |",
        "|---|---|---|",
        "| `B_low`  | 0.05–1.0 | 1P (0.125), 3P (0.376), f1_flap (0.554), f1_edge (0.629) |",
        "| `B_mid`  | 1.0–5.0  | 2nd flap (1.69), 2nd edge (1.98), 1st torsional (4.54) |",
        "| `B_high` | 5.0–25.0 | sin contraparte modal estructural conocida |",
        "",
        "## σ_Q por banda (corotacional y diagnostico opcional)",
        "",
        "| Yaw [°] | Banda | σ corot [kN·m] | σ inert [kN·m] | ratio i/c |",
        "|---:|---|---:|---:|---:|",
    ]
    for yaw in YAW_LIST:
        for b in band_labels:
            c = df_out[(df_out["solver"] == "corotational") & (df_out["yaw_deg"] == yaw)]
            i = df_out[(df_out["solver"] == "inertial") & (df_out["yaw_deg"] == yaw)]
            if c.empty or i.empty:
                continue
            cv = float(c[f"sigma_{b}_kNm"].iloc[0])
            iv = float(i[f"sigma_{b}_kNm"].iloc[0])
            ratio = iv / cv if cv > 0 else float("nan")
            md.append(f"| {yaw} | {b} | {cv:.4f} | {iv:.4f} | {ratio:.2f} |")

    md += [
        "",
        "## Fracción de varianza total por banda (corotacional)",
        "",
        "| Yaw [°] | B_low | B_mid | B_high |",
        "|---:|---:|---:|---:|",
    ]
    for yaw in YAW_LIST:
        c = df_out[(df_out["solver"] == "corotational") & (df_out["yaw_deg"] == yaw)]
        if c.empty:
            continue
        md.append(
            f"| {yaw} | {float(c['frac_B_low'].iloc[0]):.4f} | "
            f"{float(c['frac_B_mid'].iloc[0]):.4f} | {float(c['frac_B_high'].iloc[0]):.4f} |"
        )

    md += [
        "",
        "## Fracción de varianza total por banda (inercial)",
        "",
        "| Yaw [°] | B_low | B_mid | B_high |",
        "|---:|---:|---:|---:|",
    ]
    for yaw in YAW_LIST:
        i = df_out[(df_out["solver"] == "inertial") & (df_out["yaw_deg"] == yaw)]
        if i.empty:
            continue
        md.append(
            f"| {yaw} | {float(i['frac_B_low'].iloc[0]):.4f} | "
            f"{float(i['frac_B_mid'].iloc[0]):.4f} | {float(i['frac_B_high'].iloc[0]):.4f} |"
        )

    md += [
        "",
        "## Parseval sanity check",
        "",
        f"Ratio var_freq / var_time ∈ [{par_min:.3f}, {par_max:.3f}] — expected ~1.0.",
        "Pequeñas desviaciones (<10 %) son aceptables y provienen de la ventana Hann + Δf finito.",
    ]
    md_path = OUT_DATA / "torque_subband_variance_summary.md"
    md_path.write_text("\n".join(md) + "\n")
    print(f"[SB] Wrote {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
