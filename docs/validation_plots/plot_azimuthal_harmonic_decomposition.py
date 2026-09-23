"""Azimuthal and harmonic decomposition per revolution — IEA 15 MW, rated.

Reconstructs rotor azimuth from the steady-state time series under constant
rotational speed, folds the signals by revolution, and quantifies the 1P/2P/3P
harmonic budget across the yaw sweep for thrust, torque, power, and tip deformation.

Outputs
-------
docs/figures/fig_5_4_7f_azimuthal_folded_signals.png
docs/figures/fig_5_4_7f_azimuthal_folded_signals.pdf
docs/figures/fig_5_4_7g_harmonic_budget.png
docs/figures/fig_5_4_7g_harmonic_budget.pdf
docs/validation_data/generated/azimuthal_folded_signals.csv
docs/validation_data/generated/azimuthal_harmonics.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/scratch/leahk/eduardo.donestevez")
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

ROOTS = {
    "corot": BASE / "frontiersin_results_corotational_100s",
    "inert": BASE / "frontiersin_results_inertial",
}
YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
MIN_STEADY_WINDOW = 20.0
OMEGA_RPM = 7.55
OMEGA_RAD_S = OMEGA_RPM * 2.0 * np.pi / 60.0
AZIMUTH_BINS_DEG = np.linspace(0.0, 360.0, 73)
HARMONIC_ORDERS = (1, 2, 3)

SIGNALS = {
    "thrust": {"column": "Thrust [N]", "scale": 1e6, "unit": "MN", "label": "Thrust"},
    "torque": {"column": "Torque [N.m]", "scale": 1e6, "unit": "MN·m", "label": "Torque"},
    "power": {"column": "Power [W]", "scale": 1e6, "unit": "MW", "label": "Potencia"},
    "edgewise": {"column": "Tip Disp X [m]", "scale": 1.0, "unit": "m", "label": "Edgewise"},
    "flapwise": {"column": "Tip Disp Y [m]", "scale": 1.0, "unit": "m", "label": "Flapwise"},
    "resultante": {"column": "Tip Disp Mag [m]", "scale": 1.0, "unit": "m", "label": "Resultante"},
}

COLORS = {
    "thrust": "#2980b9",
    "torque": "#c0392b",
    "power": "#f39c12",
    "edgewise": "#8e44ad",
    "flapwise": "#16a085",
    "resultante": "#2c3e50",
}

HARMONIC_COLORS = {
    "1P": "#2980b9",
    "2P": "#16a085",
    "3P": "#f39c12",
    "Residual": "#bdc3c7",
}

SOLVER_TITLES = {
    "corot": "Corrotacional",
    "inert": "Inercial",
}


def load_case(csv_path: Path) -> tuple[pd.DataFrame | None, dict[str, float | str]]:
    cols = [cfg["column"] for cfg in SIGNALS.values()]
    df = pd.read_csv(csv_path)
    finite = np.asarray(np.isfinite(df[cols]).all(axis=1), dtype=bool)
    physical = np.ones(len(df), dtype=bool)
    physical &= df["Tip Disp Mag [m]"].abs().to_numpy() < 50.0
    physical &= df["Thrust [N]"].abs().to_numpy() < 1.0e8
    physical &= df["Torque [N.m]"].abs().to_numpy() < 1.0e8
    physical &= df["Power [W]"].abs().to_numpy() < 1.0e8
    valid = finite & physical

    first_bad = int(np.argmax(~valid)) if np.any(~valid) else len(df)
    truncated = df.iloc[:first_bad].copy() if first_bad < len(df) else df.copy()
    t_end_valid = float(truncated["Time [s]"].iloc[-1]) if not truncated.empty else float("nan")
    steady = truncated[
        (truncated["Time [s]"] >= T_START) & (truncated["Time [s]"] <= T_END)
    ].reset_index(drop=True)
    steady_window = min(t_end_valid, T_END) - T_START if not np.isnan(t_end_valid) else float("nan")

    status = "ok"
    if truncated.empty:
        status = "empty"
    elif steady.empty or steady_window < MIN_STEADY_WINDOW:
        status = "insufficient_steady_window"

    meta = {
        "rows_total": int(len(df)),
        "rows_valid": int(len(truncated)),
        "t_end_valid": t_end_valid,
        "steady_window": steady_window,
        "status": status,
    }
    if status != "ok":
        return None, meta
    return steady, meta


def reconstruct_azimuth_deg(time_s: np.ndarray) -> np.ndarray:
    return np.rad2deg(OMEGA_RAD_S * (time_s - float(time_s[0]))) % 360.0


def fit_single_harmonic(
    signal: np.ndarray, azimuth_deg: np.ndarray, order: int
) -> dict[str, float]:
    psi = np.deg2rad(azimuth_deg)
    signal = np.asarray(signal, dtype=float)
    centered = signal - float(np.mean(signal))
    design = np.column_stack([
        np.cos(order * psi),
        np.sin(order * psi),
    ])
    coeffs, *_ = np.linalg.lstsq(design, centered, rcond=None)
    a_cos, b_sin = coeffs
    component = a_cos * np.cos(order * psi) + b_sin * np.sin(order * psi)
    amplitude = float(np.hypot(a_cos, b_sin))
    phase_deg = float(np.rad2deg(np.arctan2(b_sin, a_cos)))
    peak_azimuth_deg = float((phase_deg / order) % 360.0)
    variance = float(np.var(centered, ddof=0))
    variance_fraction = float(np.var(component, ddof=0) / variance) if variance > 0.0 else 0.0
    return {
        "order": order,
        "amplitude": amplitude,
        "phase_deg": phase_deg,
        "peak_azimuth_deg": peak_azimuth_deg,
        "variance_fraction": variance_fraction,
    }


def harmonic_budget(
    signal: np.ndarray, azimuth_deg: np.ndarray, orders: tuple[int, ...]
) -> tuple[list[dict[str, float]], float]:
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    total_variance = float(np.var(centered, ddof=0))
    rows = []
    used = 0.0
    for order in orders:
        fit = fit_single_harmonic(signal, azimuth_deg, order)
        rows.append(fit)
        used += fit["variance_fraction"]
    residual_fraction = max(0.0, 1.0 - used) if total_variance > 0.0 else 0.0
    return rows, residual_fraction


def fold_signal(
    signal: np.ndarray, azimuth_deg: np.ndarray, azimuth_offset_deg: float
) -> pd.DataFrame:
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    std = float(np.std(centered, ddof=1))
    scaled = centered / std if std > 0.0 else np.zeros_like(centered)
    azimuth_aligned = (azimuth_deg - azimuth_offset_deg) % 360.0
    bin_ids = np.digitize(azimuth_aligned, AZIMUTH_BINS_DEG, right=False) - 1
    bin_ids = np.clip(bin_ids, 0, len(AZIMUTH_BINS_DEG) - 2)
    centers = 0.5 * (AZIMUTH_BINS_DEG[:-1] + AZIMUTH_BINS_DEG[1:])
    means = np.full_like(centers, np.nan, dtype=float)
    counts = np.zeros_like(centers, dtype=int)
    for idx in range(len(centers)):
        mask = bin_ids == idx
        if np.any(mask):
            means[idx] = float(np.mean(scaled[mask]))
            counts[idx] = int(np.sum(mask))
    return pd.DataFrame({
        "azimuth_deg": centers,
        "mean_norm": means,
        "count": counts,
    })


harmonic_rows: list[dict[str, float | int | str]] = []
fold_rows: list[dict[str, float | int | str]] = []
case_summaries: dict[tuple[str, int], dict[str, float | str]] = {}

for solver_name, root in ROOTS.items():
    for yaw in YAWS:
        csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
        if not csv_path.exists():
            continue
        steady_df, meta = load_case(csv_path)
        if steady_df is None:
            case_summaries[(solver_name, yaw)] = meta
            continue

        azimuth_deg = reconstruct_azimuth_deg(steady_df["Time [s]"].to_numpy(dtype=float))
        torque_fit_1p = fit_single_harmonic(
            steady_df[SIGNALS["torque"]["column"]].to_numpy(dtype=float), azimuth_deg, 1
        )
        torque_peak_deg = torque_fit_1p["peak_azimuth_deg"]

        case_summaries[(solver_name, yaw)] = {
            **meta,
            "torque_1p_peak_deg": torque_peak_deg,
        }

        for signal_name, signal_cfg in SIGNALS.items():
            values = steady_df[signal_cfg["column"]].to_numpy(dtype=float)
            fit_rows, residual_fraction = harmonic_budget(values, azimuth_deg, HARMONIC_ORDERS)
            signal_mean = float(np.mean(values)) / signal_cfg["scale"]
            signal_std = float(np.std(values, ddof=1)) / signal_cfg["scale"]
            for fit in fit_rows:
                harmonic_rows.append({
                    "solver": solver_name,
                    "yaw_deg": yaw,
                    "signal": signal_name,
                    "status": "ok",
                    "signal_mean": signal_mean,
                    "signal_std": signal_std,
                    "torque_1p_peak_deg": torque_peak_deg,
                    "harmonic": f"{fit['order']}P",
                    "order": fit["order"],
                    "amplitude": fit["amplitude"] / signal_cfg["scale"],
                    "phase_deg": fit["phase_deg"],
                    "peak_azimuth_deg": fit["peak_azimuth_deg"],
                    "variance_fraction": fit["variance_fraction"],
                    "unit": signal_cfg["unit"],
                })
            harmonic_rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                "signal": signal_name,
                "status": "ok",
                "signal_mean": signal_mean,
                "signal_std": signal_std,
                "torque_1p_peak_deg": torque_peak_deg,
                "harmonic": "Residual",
                "order": -1,
                "amplitude": float("nan"),
                "phase_deg": float("nan"),
                "peak_azimuth_deg": float("nan"),
                "variance_fraction": residual_fraction,
                "unit": signal_cfg["unit"],
            })

            folded = fold_signal(values, azimuth_deg, torque_peak_deg)
            for record in folded.to_dict(orient="records"):
                fold_rows.append({
                    "solver": solver_name,
                    "yaw_deg": yaw,
                    "signal": signal_name,
                    "torque_1p_peak_deg": float(torque_peak_deg),
                    "azimuth_deg": float(record["azimuth_deg"]),
                    "mean_norm": (
                        float(record["mean_norm"])
                        if np.isfinite(float(record["mean_norm"]))
                        else float("nan")
                    ),
                    "count": int(record["count"]),
                })

harmonics_df = pd.DataFrame(harmonic_rows)
folded_df = pd.DataFrame(fold_rows)

harmonics_out = DATA_DIR / "azimuthal_harmonics.csv"
folded_out = DATA_DIR / "azimuthal_folded_signals.csv"
harmonics_df.to_csv(harmonics_out, index=False)
folded_df.to_csv(folded_out, index=False)

print("\nAZIMUTHAL / HARMONIC ANALYSIS")
for solver_name in ("corot", "inert"):
    print(f"\nsolver={solver_name}")
    for yaw in YAWS:
        summary = case_summaries.get((solver_name, yaw), {})
        if summary.get("status") != "ok":
            print(f"  yaw={yaw:>2}: skipped ({summary.get('status', 'missing')})")
            continue
        torque_case = harmonics_df[
            (harmonics_df["solver"] == solver_name)
            & (harmonics_df["yaw_deg"] == yaw)
            & (harmonics_df["signal"] == "torque")
            & (harmonics_df["harmonic"].isin(["1P", "2P", "3P"]))
        ].copy()
        dominant = torque_case.iloc[torque_case["variance_fraction"].argmax()]
        print(
            f"  yaw={yaw:>2}: torque dominant={dominant['harmonic']} "
            f"(frac={dominant['variance_fraction']:.3f}, phase={dominant['phase_deg']:+.1f} deg)"
        )

fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.2), sharex=True, sharey=True)
fig.suptitle(
    "Senales plegadas por azimut alrededor del pico 1P del torque\n"
    "Campana rated, yaw=0°; azimut cero definido por el maximo 1P de Q",
    fontsize=11,
    y=1.02,
)

fold_yaw0 = folded_df[folded_df["yaw_deg"] == 0].copy()
panel_map = {
    (0, 0): ("corot", ["thrust", "torque", "power"], "(a) Corrotacional: cargas aerodinamicas"),
    (0, 1): (
        "corot",
        ["edgewise", "flapwise", "resultante"],
        "(b) Corrotacional: deformacion de punta",
    ),
    (1, 0): ("inert", ["thrust", "torque", "power"], "(c) Inercial: cargas aerodinamicas"),
    (1, 1): ("inert", ["edgewise", "flapwise", "resultante"], "(d) Inercial: deformacion de punta"),
}

for (i, j), (solver_name, signal_names, title) in panel_map.items():
    ax = axes[i, j]
    for signal_name in signal_names:
        subset = fold_yaw0[
            (fold_yaw0["solver"] == solver_name) & (fold_yaw0["signal"] == signal_name)
        ]
        ax.plot(
            subset["azimuth_deg"],
            subset["mean_norm"],
            color=COLORS[signal_name],
            lw=1.8,
            label=SIGNALS[signal_name]["label"],
        )
    ax.axvline(0.0, color="#666", lw=0.9, ls=":")
    ax.axhline(0.0, color="#666", lw=0.9, ls=":")
    ax.set_title(title)
    ax.set_ylabel("fluctuacion normalizada")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

for ax in axes[1, :]:
    ax.set_xlabel("Azimut alineado [deg]")

fig.tight_layout()
fold_png = OUT_DIR / "fig_5_4_7f_azimuthal_folded_signals.png"
fold_pdf = OUT_DIR / "fig_5_4_7f_azimuthal_folded_signals.pdf"
fig.savefig(fold_png, dpi=300, bbox_inches="tight")
fig.savefig(fold_pdf, bbox_inches="tight")

budget_signals = ["thrust", "torque", "edgewise", "flapwise"]
budget_labels = ["1P", "2P", "3P", "Residual"]

fig2, axes2 = plt.subplots(2, 4, figsize=(19.0, 8.0), sharex=True, sharey=True)
fig2.suptitle(
    "Presupuesto armonico por yaw para las fluctuaciones centradas\n"
    "Fracciones de varianza explicadas por 1P, 2P y 3P",
    fontsize=11,
    y=1.02,
)

for row_idx, solver_name in enumerate(("corot", "inert")):
    for col_idx, signal_name in enumerate(budget_signals):
        ax = axes2[row_idx, col_idx]
        subset = harmonics_df[
            (harmonics_df["solver"] == solver_name)
            & (harmonics_df["signal"] == signal_name)
            & (harmonics_df["yaw_deg"].isin(YAWS))
        ].copy()
        pivot = subset.pivot(
            index="yaw_deg", columns="harmonic", values="variance_fraction"
        ).reindex(YAWS)
        bottom = np.zeros(len(YAWS), dtype=float)
        for label in budget_labels:
            vals = (
                pivot[label].to_numpy(dtype=float)
                if label in pivot.columns
                else np.zeros(len(YAWS), dtype=float)
            )
            ax.bar(
                YAWS,
                vals,
                bottom=bottom,
                width=7.0,
                color=HARMONIC_COLORS[label],
                label=label,
            )
            bottom += vals
        for yaw, total in zip(YAWS, bottom):
            if np.isnan(total) or total <= 0.0:
                ax.text(yaw, 0.52, "N/A", ha="center", va="center", fontsize=8, color="#555")
        ax.set_title(
            f"({chr(97 + row_idx * len(budget_signals) + col_idx)}) {SOLVER_TITLES[solver_name]}: {SIGNALS[signal_name]['label']}"
        )
        ax.set_ylim(0.0, 1.02)
        ax.grid(axis="y", alpha=0.3)
        if col_idx == 0:
            ax.set_ylabel("fraccion de varianza")
        if row_idx == 1:
            ax.set_xlabel("Yaw [deg]")

handles, labels = axes2[0, 0].get_legend_handles_labels()
unique = dict(zip(labels, handles))
axes2[0, 3].legend(unique.values(), unique.keys(), frameon=False, fontsize=8)

fig2.tight_layout()
budget_png = OUT_DIR / "fig_5_4_7g_harmonic_budget.png"
budget_pdf = OUT_DIR / "fig_5_4_7g_harmonic_budget.pdf"
fig2.savefig(budget_png, dpi=300, bbox_inches="tight")
fig2.savefig(budget_pdf, bbox_inches="tight")

print(f"\nSaved: {fold_png}")
print(f"Saved: {fold_pdf}")
print(f"Saved: {budget_png}")
print(f"Saved: {budget_pdf}")
print(f"Saved: {folded_out}")
print(f"Saved: {harmonics_out}")
