"""Paired corotational vs inertial campaign comparison.

This script assumes both campaigns are complete and comparable in the common
steady window used by the validation dossier: 40 <= t <= 100 s. It produces the
solver-to-solver evidence layer that complements the existing yaw sweep,
component influence, and FFT analyses.

Outputs
-------
docs/figures/fig_5_12_1_solver_pairwise_matrix.png
docs/figures/fig_5_12_1_solver_pairwise_matrix.pdf
docs/validation_data/generated/solver_pairwise_bem_timeseries.csv
docs/validation_data/generated/solver_pairwise_yaw_response.csv
docs/validation_data/generated/solver_pairwise_rotor_fft.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

BASE = Path("/scratch/leahk/eduardo.donestevez")
ROOTS = {
    "corotational": BASE / "frontiersin_results_corotational_100s",
    "inertial": BASE / "frontiersin_results_inertial",
}
REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / "docs" / "figures"
DATA_DIR = REPO_ROOT / "docs" / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
DT = 0.01
OMEGA_RAD_S = 7.55 * 2.0 * np.pi / 60.0
F_1P = OMEGA_RAD_S / (2.0 * np.pi)
F_3P = 3.0 * F_1P
F_FLAP = 0.5537
F_EDGE = 0.6290

BEM_SIGNALS = {
    "thrust": {"column": "Thrust [N]", "scale": 1.0e6, "unit": "MN", "label": "Thrust"},
    "torque": {"column": "Torque [N.m]", "scale": 1.0e6, "unit": "MNm", "label": "Torque"},
    "power": {"column": "Power [W]", "scale": 1.0e6, "unit": "MW", "label": "Power"},
    "edgewise": {"column": "Tip Disp X [m]", "scale": 1.0, "unit": "m", "label": "Edgewise"},
    "flapwise": {"column": "Tip Disp Y [m]", "scale": 1.0, "unit": "m", "label": "Flapwise"},
    "tip_mag": {"column": "Tip Disp Mag [m]", "scale": 1.0, "unit": "m", "label": "Tip magnitude"},
}

SUMMARY_METRICS = {
    "flap_mean": {"unit": "m", "label": "Flap mean"},
    "edge_mean": {"unit": "m", "label": "Edge mean"},
    "thrust_mean": {"unit": "N", "label": "Thrust mean"},
    "torque_mean": {"unit": "N.m", "label": "Torque mean"},
    "power_mean_mw": {"unit": "MW", "label": "Power mean"},
    "flap_std": {"unit": "m", "label": "Flap std"},
    "edge_std": {"unit": "m", "label": "Edge std"},
    "power_std_mw": {"unit": "MW", "label": "Power std"},
}

ROTOR_SIGNALS = {
    "power_rotor_eq": {
        "column": "Aero Power Rotor Equivalent [W]",
        "scale": 1.0e6,
        "unit": "MW",
        "label": "Power rotor-eq.",
    },
    "tilt_moment": {
        "column": "Aero Torque X [Nm]",
        "scale": 1.0e6,
        "unit": "MNm",
        "label": "Tilt moment",
    },
    "yaw_moment": {
        "column": "Aero Torque Z [Nm]",
        "scale": 1.0e6,
        "unit": "MNm",
        "label": "Yaw moment",
    },
    "max_displacement": {
        "column": "Max Displacement [m]",
        "scale": 1.0,
        "unit": "m",
        "label": "Max displacement",
    },
}

CHAR_FREQS = {
    "1P": F_1P,
    "3P": F_3P,
    "f1_flap": F_FLAP,
    "f1_edge": F_EDGE,
}


def steady_window(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df["Time [s]"] >= T_START) & (df["Time [s]"] <= T_END)].reset_index(drop=True)


def load_bem(solver: str, yaw: int) -> pd.DataFrame:
    path = ROOTS[solver] / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
    return steady_window(pd.read_csv(path))


def load_rotor_performance(solver: str, yaw: int) -> pd.DataFrame:
    path = ROOTS[solver] / f"yaw_{yaw}" / solver / "rotor_performance.csv"
    df = steady_window(pd.read_csv(path))
    if "Aero Power Rotor Equivalent [W]" not in df and "Aero Power [W]" in df:
        df["Aero Power Rotor Equivalent [W]"] = 3.0 * df["Aero Power [W]"]
    return df


def relative_pct(delta: float, reference: float) -> float:
    if not np.isfinite(reference) or abs(reference) < 1.0e-12:
        return float("nan")
    return 100.0 * delta / abs(reference)


def lagged_correlation(
    corotational: np.ndarray,
    inertial: np.ndarray,
    max_lag_s: float = 2.0,
) -> tuple[float, float]:
    x = np.asarray(corotational, dtype=float) - float(np.mean(corotational))
    y = np.asarray(inertial, dtype=float) - float(np.mean(inertial))
    max_lag = int(round(max_lag_s / DT))
    best_corr = -np.inf
    best_lag = 0
    for lag in range(-max_lag, max_lag + 1):
        if lag > 0:
            xs = x[:-lag]
            ys = y[lag:]
        elif lag < 0:
            xs = x[-lag:]
            ys = y[:lag]
        else:
            xs = x
            ys = y
        denom = float(np.linalg.norm(xs) * np.linalg.norm(ys))
        corr = float(np.dot(xs, ys) / denom) if denom > 0.0 else float("nan")
        if np.isfinite(corr) and corr > best_corr:
            best_corr = corr
            best_lag = lag
    return best_lag * DT, best_corr


def onesided_fft(signal: np.ndarray, dt: float = DT) -> tuple[np.ndarray, np.ndarray]:
    values = np.asarray(signal, dtype=float)
    n = len(values)
    window = np.hanning(n)
    normalized = window / float(window.mean())
    centered = (values - float(values.mean())) * normalized
    spectrum = np.fft.rfft(centered)
    freq = np.fft.rfftfreq(n, d=dt)
    amp = np.abs(spectrum) * 2.0 / n
    amp[0] = 0.0
    return freq, amp


def peak_near(freq: np.ndarray, amp: np.ndarray, target: float, bandwidth: float = 0.04) -> float:
    mask = np.abs(freq - target) <= bandwidth
    if not np.any(mask):
        return float("nan")
    return float(np.max(amp[mask]))


def compute_bem_pairwise_metrics() -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for yaw in YAWS:
        corot = load_bem("corotational", yaw)
        inert = load_bem("inertial", yaw)
        paired = pd.merge(corot, inert, on="Time [s]", suffixes=("_corot", "_inert"))
        for signal, spec in BEM_SIGNALS.items():
            col_corot = f"{spec['column']}_corot"
            col_inert = f"{spec['column']}_inert"
            scale = float(spec["scale"])
            c = paired[col_corot].to_numpy(dtype=float) / scale
            i = paired[col_inert].to_numpy(dtype=float) / scale
            c_center = c - float(np.mean(c))
            i_center = i - float(np.mean(i))
            bias = float(np.mean(i) - np.mean(c))
            rmse = float(np.sqrt(np.mean((i - c) ** 2)))
            centered_rmse = float(np.sqrt(np.mean((i_center - c_center) ** 2)))
            corr = float(np.corrcoef(c_center, i_center)[0, 1])
            lag_s, lag_corr = lagged_correlation(c, i)
            c_std = float(np.std(c, ddof=1))
            i_std = float(np.std(i, ddof=1))
            c_p2p = float(np.ptp(c))
            i_p2p = float(np.ptp(i))
            rows.append({
                "yaw_deg": yaw,
                "signal": signal,
                "label": spec["label"],
                "unit": spec["unit"],
                "n_samples": len(paired),
                "corot_mean": float(np.mean(c)),
                "inert_mean": float(np.mean(i)),
                "mean_bias_inert_minus_corot": bias,
                "mean_bias_pct": relative_pct(bias, float(np.mean(c))),
                "corot_std": c_std,
                "inert_std": i_std,
                "std_ratio_inert_over_corot": i_std / c_std if c_std else float("nan"),
                "corot_p2p": c_p2p,
                "inert_p2p": i_p2p,
                "p2p_ratio_inert_over_corot": i_p2p / c_p2p if c_p2p else float("nan"),
                "rmse": rmse,
                "centered_rmse": centered_rmse,
                "centered_rmse_over_corot_std": centered_rmse / c_std if c_std else float("nan"),
                "pearson_same_time": corr,
                "best_lag_s_positive_inert_lags": lag_s,
                "best_lag_corr": lag_corr,
            })
    return pd.DataFrame(rows)


def compute_yaw_response_metrics() -> pd.DataFrame:
    summary_path = DATA_DIR / "campaign_metrics_summary.csv"
    summary = pd.read_csv(summary_path)
    rows: list[dict[str, float | str]] = []
    for metric, spec in SUMMARY_METRICS.items():
        pivot = summary.pivot(index="yaw_deg", columns="solver", values=metric).loc[YAWS]
        diff = pivot["inertial"] - pivot["corotational"]
        diff_pct = 100.0 * diff / pivot["corotational"].abs()
        corot_drop_pct = relative_pct(
            float(pivot.loc[40, "corotational"] - pivot.loc[0, "corotational"]),
            float(pivot.loc[0, "corotational"]),
        )
        inert_drop_pct = relative_pct(
            float(pivot.loc[40, "inertial"] - pivot.loc[0, "inertial"]),
            float(pivot.loc[0, "inertial"]),
        )
        corot_slope = float(np.polyfit(YAWS, pivot["corotational"].to_numpy(dtype=float), 1)[0])
        inert_slope = float(np.polyfit(YAWS, pivot["inertial"].to_numpy(dtype=float), 1)[0])
        rows.append({
            "metric": metric,
            "label": spec["label"],
            "unit": spec["unit"],
            "corot_yaw0": float(pivot.loc[0, "corotational"]),
            "inert_yaw0": float(pivot.loc[0, "inertial"]),
            "yaw0_bias_pct": float(diff_pct.loc[0]),
            "corot_yaw40": float(pivot.loc[40, "corotational"]),
            "inert_yaw40": float(pivot.loc[40, "inertial"]),
            "yaw40_bias_pct": float(diff_pct.loc[40]),
            "mean_abs_bias_pct": float(diff_pct.abs().mean()),
            "max_abs_bias_pct": float(diff_pct.abs().max()),
            "corot_drop_0_to_40_pct": corot_drop_pct,
            "inert_drop_0_to_40_pct": inert_drop_pct,
            "drop_delta_pct_points": inert_drop_pct - corot_drop_pct,
            "corot_slope_per_yaw_deg": corot_slope,
            "inert_slope_per_yaw_deg": inert_slope,
        })

    power = summary.pivot(index="yaw_deg", columns="solver", values="power_mean_mw").loc[YAWS]
    x = np.log(np.cos(np.deg2rad(np.array(YAWS[1:], dtype=float))))
    for solver in ["corotational", "inertial"]:
        y = np.log(power.loc[YAWS[1:], solver].to_numpy(dtype=float) / float(power.loc[0, solver]))
        rows.append({
            "metric": "power_cos_exponent",
            "label": f"Power cos exponent ({solver})",
            "unit": "-",
            "corot_yaw0": float("nan"),
            "inert_yaw0": float("nan"),
            "yaw0_bias_pct": float("nan"),
            "corot_yaw40": float("nan"),
            "inert_yaw40": float("nan"),
            "yaw40_bias_pct": float("nan"),
            "mean_abs_bias_pct": float("nan"),
            "max_abs_bias_pct": float("nan"),
            "corot_drop_0_to_40_pct": float("nan"),
            "inert_drop_0_to_40_pct": float("nan"),
            "drop_delta_pct_points": float("nan"),
            "corot_slope_per_yaw_deg": float(np.polyfit(x, y, 1)[0])
            if solver == "corotational"
            else float("nan"),
            "inert_slope_per_yaw_deg": float(np.polyfit(x, y, 1)[0])
            if solver == "inertial"
            else float("nan"),
        })
    return pd.DataFrame(rows)


def compute_rotor_fft_metrics() -> pd.DataFrame:
    rows: list[dict[str, float | int | str]] = []
    for yaw in YAWS:
        amps_by_solver: dict[tuple[str, str, str], float] = {}
        for solver in ["corotational", "inertial"]:
            df = load_rotor_performance(solver, yaw)
            for signal, spec in ROTOR_SIGNALS.items():
                values = df[str(spec["column"])].to_numpy(dtype=float) / float(spec["scale"])
                freq, amp = onesided_fft(values)
                for freq_label, freq_target in CHAR_FREQS.items():
                    amps_by_solver[(solver, signal, freq_label)] = peak_near(freq, amp, freq_target)
        for signal, spec in ROTOR_SIGNALS.items():
            for freq_label, freq_target in CHAR_FREQS.items():
                corot_amp = amps_by_solver[("corotational", signal, freq_label)]
                inert_amp = amps_by_solver[("inertial", signal, freq_label)]
                rows.append({
                    "yaw_deg": yaw,
                    "signal": signal,
                    "label": spec["label"],
                    "unit": spec["unit"],
                    "frequency_label": freq_label,
                    "frequency_hz": freq_target,
                    "corot_amplitude": corot_amp,
                    "inert_amplitude": inert_amp,
                    "amplitude_delta_inert_minus_corot": inert_amp - corot_amp,
                    "amplitude_ratio_inert_over_corot": inert_amp / corot_amp
                    if corot_amp
                    else float("nan"),
                })
    return pd.DataFrame(rows)


def plot_pairwise_matrix(
    bem_metrics: pd.DataFrame,
    yaw_response: pd.DataFrame,
    rotor_fft: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), constrained_layout=True)
    fig.suptitle("Matriz pareada inercial vs corrotacional — ventana 40-100 s", fontsize=13)

    ax = axes[0, 0]
    mean_signals = ["flapwise", "edgewise", "thrust", "torque", "power"]
    for signal in mean_signals:
        subset = bem_metrics[bem_metrics["signal"] == signal]
        ax.plot(
            subset["yaw_deg"],
            subset["mean_bias_pct"],
            marker="o",
            label=str(subset["label"].iloc[0]),
        )
    ax.axhline(0.0, color="0.25", lw=0.8)
    ax.set_title("Sesgo medio: inercial - corrotacional")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("Diferencia media [%]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False, ncol=2)

    ax = axes[0, 1]
    dynamic_signals = ["flapwise", "edgewise", "power"]
    for signal in dynamic_signals:
        subset = bem_metrics[bem_metrics["signal"] == signal]
        ax.plot(
            subset["yaw_deg"],
            subset["std_ratio_inert_over_corot"],
            marker="s",
            label=str(subset["label"].iloc[0]),
        )
    ax.axhline(1.0, color="0.25", lw=0.8)
    ax.set_title("Relación dinámica en bem_report.csv")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("σ inercial / σ corrotacional")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 0]
    for signal in ["flapwise", "torque", "power"]:
        subset = bem_metrics[bem_metrics["signal"] == signal]
        ax.plot(
            subset["yaw_deg"],
            subset["pearson_same_time"],
            marker="^",
            label=str(subset["label"].iloc[0]),
        )
    ax.set_title("Correlación temporal pareada")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("r(corot, inert)")
    ax.set_ylim(-1.05, 1.05)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1, 1]
    subset = rotor_fft[rotor_fft["frequency_label"] == "1P"]
    for signal in ["power_rotor_eq", "tilt_moment", "yaw_moment", "max_displacement"]:
        signal_subset = subset[subset["signal"] == signal]
        ax.plot(
            signal_subset["yaw_deg"],
            signal_subset["amplitude_ratio_inert_over_corot"],
            marker="D",
            label=str(signal_subset["label"].iloc[0]),
        )
    ax.axhline(1.0, color="0.25", lw=0.8)
    ax.set_yscale("log")
    ax.set_title("Anomalía 1P en rotor_performance.csv")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("Amplitud 1P inercial / corrotacional")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8, frameon=False)

    png = OUT_DIR / "fig_5_12_1_solver_pairwise_matrix.png"
    pdf = OUT_DIR / "fig_5_12_1_solver_pairwise_matrix.pdf"
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")


def main() -> None:
    bem_metrics = compute_bem_pairwise_metrics()
    yaw_response = compute_yaw_response_metrics()
    rotor_fft = compute_rotor_fft_metrics()

    bem_path = DATA_DIR / "solver_pairwise_bem_timeseries.csv"
    yaw_path = DATA_DIR / "solver_pairwise_yaw_response.csv"
    fft_path = DATA_DIR / "solver_pairwise_rotor_fft.csv"
    bem_metrics.to_csv(bem_path, index=False)
    yaw_response.to_csv(yaw_path, index=False)
    rotor_fft.to_csv(fft_path, index=False)
    print(f"Saved: {bem_path}")
    print(f"Saved: {yaw_path}")
    print(f"Saved: {fft_path}")

    plot_pairwise_matrix(bem_metrics, yaw_response, rotor_fft)

    key_metrics = yaw_response[
        yaw_response["metric"].isin(["flap_mean", "power_mean_mw", "thrust_mean", "torque_mean"])
    ]
    print("\nMean-pairing summary")
    print(
        key_metrics[
            ["metric", "mean_abs_bias_pct", "max_abs_bias_pct", "drop_delta_pct_points"]
        ].to_string(index=False)
    )

    fft_1p_yaw0 = rotor_fft[(rotor_fft["yaw_deg"] == 0) & (rotor_fft["frequency_label"] == "1P")]
    print("\nRotor-performance 1P ratios at yaw=0")
    print(
        fft_1p_yaw0[
            [
                "signal",
                "corot_amplitude",
                "inert_amplitude",
                "amplitude_ratio_inert_over_corot",
                "unit",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
