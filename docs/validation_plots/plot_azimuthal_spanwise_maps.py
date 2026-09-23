r"""Azimuth-span modulation of sectional aerodynamic quantities under yaw.

Builds folded azimuth-span maps from the per-time-step ``bem_sectional.csv``
files of the corotational rated campaign. The absolute blade azimuth is not
exported by the current BEM reports, so the postprocess reconstructs azimuth
from time assuming constant rotor speed and aligns ``psi = 0`` with the 1P peak
of the blade-integrated tangential load ``\int T_p(r) r dr``.

Outputs
-------
docs/figures/fig_5_4_7j_azimuthal_spanwise_maps.png
docs/figures/fig_5_4_7j_azimuthal_spanwise_maps.pdf
docs/figures/fig_5_4_7k_spanwise_cyclic_intensity.png
docs/figures/fig_5_4_7k_spanwise_cyclic_intensity.pdf
docs/validation_data/generated/azimuthal_spanwise_folded.csv
docs/validation_data/generated/azimuthal_spanwise_metrics.csv
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, TypedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class SignalResult(TypedDict):
    centers: np.ndarray
    folded: np.ndarray
    matrix: np.ndarray


class SectionalWindow(TypedDict):
    time_s: np.ndarray
    torque_proxy: np.ndarray
    r_m: np.ndarray
    r_R: np.ndarray
    signals: dict[str, np.ndarray]


class CaseResult(TypedDict):
    r_m: np.ndarray
    r_R: np.ndarray
    align_deg: float
    signal_results: dict[str, SignalResult]


BASE = Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s"
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

YAWS = [0, 10, 20, 30, 40]
HEATMAP_YAWS = [0, 20, 40]
HEATMAP_SIGNALS = ["Np", "Tp", "alpha"]
T_START = 40.0
T_END = 100.0
DT_SAMPLE = 0.05
ROOT_CUTOFF_R_R = 0.16
OMEGA_RPM = 7.55
OMEGA_RAD_S = OMEGA_RPM * 2.0 * np.pi / 60.0
AZIMUTH_BINS_DEG = np.linspace(0.0, 360.0, 73)
EPS = 1.0e-12

SIGNALS = {
    "Np": {
        "column": "Np[N/m]",
        "label": r"$N_p'$",
        "title": "Fuerza normal",
        "unit": "kN/m",
        "scale": 1.0e3,
    },
    "Tp": {
        "column": "Tp[N/m]",
        "label": r"$T_p'$",
        "title": "Fuerza tangencial",
        "unit": "kN/m",
        "scale": 1.0e3,
    },
    "alpha": {
        "column": "alpha[deg]",
        "label": r"$\alpha'$",
        "title": "AoA",
        "unit": "deg",
        "scale": 1.0,
    },
    "a": {
        "column": "a",
        "label": r"$a$",
        "title": "Induccion axial",
        "unit": "-",
        "scale": 1.0,
    },
    "ap": {
        "column": "ap",
        "label": r"$a'$",
        "title": "Induccion tangencial",
        "unit": "-",
        "scale": 1.0,
    },
}


def list_time_dirs(root: Path) -> Iterable[Path]:
    dirs = [path for path in root.iterdir() if path.is_dir()]
    return sorted(dirs, key=lambda path: float(path.name))


def reconstruct_azimuth_deg(time_s: np.ndarray) -> np.ndarray:
    return np.rad2deg(OMEGA_RAD_S * (time_s - float(time_s[0]))) % 360.0


def fit_single_harmonic(
    signal: np.ndarray, azimuth_deg: np.ndarray, order: int = 1
) -> dict[str, float]:
    psi = np.deg2rad(np.asarray(azimuth_deg, dtype=float))
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    design = np.column_stack([np.cos(order * psi), np.sin(order * psi)])
    coeffs, *_ = np.linalg.lstsq(design, centered, rcond=None)
    a_cos, b_sin = coeffs
    component = a_cos * np.cos(order * psi) + b_sin * np.sin(order * psi)
    variance = float(np.var(centered, ddof=0))
    return {
        "amplitude": float(np.hypot(a_cos, b_sin)),
        "phase_deg": float(np.rad2deg(np.arctan2(b_sin, a_cos))),
        "peak_azimuth_deg": float((np.rad2deg(np.arctan2(b_sin, a_cos)) / order) % 360.0),
        "variance_fraction": float(np.var(component, ddof=0) / variance) if variance > 0.0 else 0.0,
    }


def fold_matrix(
    matrix: np.ndarray, azimuth_deg: np.ndarray, azimuth_offset_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    azimuth_aligned = (np.asarray(azimuth_deg, dtype=float) - azimuth_offset_deg) % 360.0
    bin_ids = np.digitize(azimuth_aligned, AZIMUTH_BINS_DEG, right=False) - 1
    bin_ids = np.clip(bin_ids, 0, len(AZIMUTH_BINS_DEG) - 2)
    centers = 0.5 * (AZIMUTH_BINS_DEG[:-1] + AZIMUTH_BINS_DEG[1:])
    folded = np.full((len(centers), matrix.shape[1]), np.nan, dtype=float)
    for idx in range(len(centers)):
        mask = bin_ids == idx
        if np.any(mask):
            folded[idx, :] = np.mean(matrix[mask, :], axis=0)
    return centers, folded


def cell_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mids = 0.5 * (values[:-1] + values[1:])
    start = values[0] - 0.5 * (values[1] - values[0])
    end = values[-1] + 0.5 * (values[-1] - values[-2])
    return np.concatenate([[start], mids, [end]])


def load_sectional_window(fluid_root: Path) -> SectionalWindow:
    times: list[float] = []
    torque_proxy: list[float] = []
    signal_samples: dict[str, list[np.ndarray]] = {name: [] for name in SIGNALS}
    r_ref: np.ndarray | None = None
    r_r_ref: np.ndarray | None = None

    last_t = -np.inf
    for time_dir in list_time_dirs(fluid_root):
        try:
            t = float(time_dir.name)
        except ValueError:
            continue
        if t < T_START or t > T_END:
            continue
        if t - last_t < DT_SAMPLE:
            continue

        csv_path = time_dir / "bem_sectional.csv"
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        r = df["r[m]"].to_numpy(dtype=float)
        if r_ref is None:
            r_over_r = r / float(np.max(r))
            mask = r_over_r >= ROOT_CUTOFF_R_R
            r_ref = r[mask]
            r_r_ref = r_over_r[mask]
        else:
            r_over_r = r / float(np.max(r))
            mask = r_over_r >= ROOT_CUTOFF_R_R

        for signal_name, cfg in SIGNALS.items():
            signal_samples[signal_name].append(df.loc[mask, cfg["column"]].to_numpy(dtype=float))

        torque_proxy.append(float(np.trapezoid(df["Tp[N/m]"].to_numpy(dtype=float) * r, r)))
        times.append(t)
        last_t = t

    if not times or r_ref is None or r_r_ref is None:
        raise RuntimeError(f"No bem_sectional.csv samples found in {fluid_root}")

    return {
        "time_s": np.asarray(times, dtype=float),
        "torque_proxy": np.asarray(torque_proxy, dtype=float),
        "r_m": r_ref,
        "r_R": r_r_ref,
        "signals": {name: np.vstack(rows) for name, rows in signal_samples.items()},
    }


def summarize_signal(
    yaw: int,
    signal_name: str,
    r_m: np.ndarray,
    r_r: np.ndarray,
    azimuth_deg: np.ndarray,
    aligned_azimuth_deg: float,
    matrix: np.ndarray,
    folded: np.ndarray,
) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    rows_metrics: list[dict[str, float | int | str]] = []
    rows_folded: list[dict[str, float | int | str]] = []
    mean_profile = np.mean(matrix, axis=0)
    std_profile = np.std(matrix, axis=0, ddof=1)
    centered_folded = folded - mean_profile[None, :]
    peak_to_peak = np.nanmax(centered_folded, axis=0) - np.nanmin(centered_folded, axis=0)

    for span_idx, (radius_m, radius_r, mean_val, std_val, ptp_val) in enumerate(
        zip(r_m, r_r, mean_profile, std_profile, peak_to_peak, strict=True)
    ):
        fit = fit_single_harmonic(matrix[:, span_idx], azimuth_deg, order=1)
        peak_aligned = (fit["peak_azimuth_deg"] - aligned_azimuth_deg) % 360.0
        rows_metrics.append({
            "yaw_deg": yaw,
            "signal": signal_name,
            "r_m": float(radius_m),
            "r_R": float(radius_r),
            "mean_value": float(mean_val),
            "std_value": float(std_val),
            "cyclic_intensity_pct": 100.0 * float(std_val) / max(abs(float(mean_val)), EPS),
            "peak_to_peak_value": float(ptp_val),
            "peak_to_peak_pct": 100.0 * float(ptp_val) / max(abs(float(mean_val)), EPS),
            "one_p_amplitude": fit["amplitude"],
            "one_p_amplitude_pct": 100.0 * fit["amplitude"] / max(abs(float(mean_val)), EPS),
            "one_p_variance_fraction": fit["variance_fraction"],
            "one_p_peak_azimuth_aligned_deg": peak_aligned,
        })

    centers = 0.5 * (AZIMUTH_BINS_DEG[:-1] + AZIMUTH_BINS_DEG[1:])
    for azimuth_bin, azimuth_center in enumerate(centers):
        for span_idx, (radius_m, radius_r) in enumerate(zip(r_m, r_r, strict=True)):
            folded_value = float(folded[azimuth_bin, span_idx])
            rows_folded.append({
                "yaw_deg": yaw,
                "signal": signal_name,
                "azimuth_deg": float(azimuth_center),
                "r_m": float(radius_m),
                "r_R": float(radius_r),
                "folded_value": folded_value,
                "folded_centered_value": folded_value - float(mean_profile[span_idx]),
            })

    return rows_metrics, rows_folded


results: dict[int, CaseResult] = {}
metrics_rows: list[dict[str, float | int | str]] = []
folded_rows: list[dict[str, float | int | str]] = []

for yaw in YAWS:
    fluid_root = CAMPAIGN_ROOT / f"yaw_{yaw}" / "fluid"
    if not fluid_root.exists():
        continue

    case = load_sectional_window(fluid_root)
    azimuth_deg = reconstruct_azimuth_deg(case["time_s"])
    torque_fit = fit_single_harmonic(case["torque_proxy"], azimuth_deg, order=1)
    align_deg = torque_fit["peak_azimuth_deg"]

    signal_results: dict[str, SignalResult] = {}
    for signal_name in SIGNALS:
        matrix = case["signals"][signal_name]
        centers, folded = fold_matrix(matrix, azimuth_deg, align_deg)
        signal_results[signal_name] = {"centers": centers, "folded": folded, "matrix": matrix}

        signal_metrics, signal_folded = summarize_signal(
            yaw=yaw,
            signal_name=signal_name,
            r_m=case["r_m"],
            r_r=case["r_R"],
            azimuth_deg=azimuth_deg,
            aligned_azimuth_deg=float(align_deg),
            matrix=matrix,
            folded=folded,
        )
        metrics_rows.extend(signal_metrics)
        folded_rows.extend(signal_folded)

    results[yaw] = {
        "r_m": case["r_m"],
        "r_R": case["r_R"],
        "align_deg": align_deg,
        "signal_results": signal_results,
    }

metrics_df = pd.DataFrame(metrics_rows)
folded_df = pd.DataFrame(folded_rows)
metrics_out = DATA_DIR / "azimuthal_spanwise_metrics.csv"
folded_out = DATA_DIR / "azimuthal_spanwise_folded.csv"
metrics_df.to_csv(metrics_out, index=False)
folded_df.to_csv(folded_out, index=False)

print("\nAZIMUTH-SPAN ANALYSIS")
for signal_name in SIGNALS:
    print(f"\nSignal: {signal_name}")
    subset = metrics_df[metrics_df["signal"] == signal_name].copy()
    for yaw in YAWS:
        yaw_subset = subset[subset["yaw_deg"] == yaw]
        if yaw_subset.empty:
            print(f"  yaw={yaw:>2}: missing")
            continue
        peak_row = yaw_subset.iloc[int(yaw_subset["cyclic_intensity_pct"].argmax())]
        print(
            f"  yaw={yaw:>2}: max cyclic intensity = {peak_row['cyclic_intensity_pct']:.1f}% "
            f"at r/R={peak_row['r_R']:.3f}; 1P frac mean = {yaw_subset['one_p_variance_fraction'].mean():.3f}"
        )

x_edges = AZIMUTH_BINS_DEG
y_edges = cell_edges(next(iter(results.values()))["r_R"])

fig, axes = plt.subplots(
    len(HEATMAP_SIGNALS), len(HEATMAP_YAWS), figsize=(13.8, 9.0), sharex=True, sharey=True
)
fig.suptitle(
    "Una pala desenrollada en fase de giro y posicion spanwise — AeroElast corrotacional\n"
    "Rojo: por encima de la media local de esa estacion; azul: por debajo. psi=0 es una referencia interna, no azimut geometrico absoluto",
    fontsize=11,
    y=1.02,
)

for row_idx, signal_name in enumerate(HEATMAP_SIGNALS):
    cfg = SIGNALS[signal_name]
    stacked = []
    for yaw in HEATMAP_YAWS:
        folded = results[yaw]["signal_results"][signal_name]["folded"]
        mean_profile = np.mean(results[yaw]["signal_results"][signal_name]["matrix"], axis=0)
        stacked.append((folded - mean_profile[None, :]) / cfg["scale"])
    limit = float(np.nanpercentile(np.abs(np.concatenate(stacked, axis=0)), 98.0))
    limit = max(limit, 1.0e-6)

    for col_idx, yaw in enumerate(HEATMAP_YAWS):
        ax = axes[row_idx, col_idx]
        folded = results[yaw]["signal_results"][signal_name]["folded"]
        mean_profile = np.mean(results[yaw]["signal_results"][signal_name]["matrix"], axis=0)
        centered = (folded - mean_profile[None, :]) / cfg["scale"]
        mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            centered.T,
            shading="auto",
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
        )
        if row_idx == 0:
            ax.set_title(f"yaw = {yaw}°")
        if col_idx == 0:
            ax.set_ylabel(rf"Posicion spanwise $r/R$ [-]\\n{cfg['title']}")
        ax.grid(alpha=0.15)
        if row_idx == len(HEATMAP_SIGNALS) - 1:
            ax.set_xlabel("Fase dentro de una vuelta [deg]")
        fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.02, label=cfg["unit"])

fig.tight_layout()
heatmap_png = OUT_DIR / "fig_5_4_7j_azimuthal_spanwise_maps.png"
heatmap_pdf = OUT_DIR / "fig_5_4_7j_azimuthal_spanwise_maps.pdf"
fig.savefig(heatmap_png, dpi=300, bbox_inches="tight")
fig.savefig(heatmap_pdf, bbox_inches="tight")

fig2, axes2 = plt.subplots(2, 3, figsize=(14.0, 8.0), sharex=True)
fig2.suptitle(
    "Resumen del mapa anterior: cuanto oscila cada estacion radial durante una vuelta\n"
    r"Metrica: $100\,\sigma_{\psi}(x)/|\mu_{\psi}(x)|$ sobre la ventana estacionaria",
    fontsize=11,
    y=1.02,
)

cmap = plt.get_cmap("viridis")
colors = {yaw: cmap(idx / max(1, len(YAWS) - 1)) for idx, yaw in enumerate(YAWS)}
signal_order = ["Np", "Tp", "alpha", "a", "ap"]

for panel_idx, signal_name in enumerate(signal_order):
    ax = axes2.flat[panel_idx]
    subset = metrics_df[metrics_df["signal"] == signal_name].copy()
    for yaw in YAWS:
        yaw_subset = subset[subset["yaw_deg"] == yaw].sort_values("r_R")
        if yaw_subset.empty:
            continue
        ax.plot(
            yaw_subset["r_R"],
            yaw_subset["cyclic_intensity_pct"],
            color=colors[yaw],
            lw=1.8,
            label=f"yaw = {yaw}°",
        )
    ax.set_title(SIGNALS[signal_name]["title"])
    ax.set_ylabel("intensidad ciclica [%]")
    ax.grid(alpha=0.3)

axes2.flat[-1].axis("off")
for ax in axes2[1, :2]:
    ax.set_xlabel(r"Posicion spanwise $r/R$ [-]")
axes2[0, 0].legend(frameon=False, fontsize=8)

fig2.tight_layout()
severity_png = OUT_DIR / "fig_5_4_7k_spanwise_cyclic_intensity.png"
severity_pdf = OUT_DIR / "fig_5_4_7k_spanwise_cyclic_intensity.pdf"
fig2.savefig(severity_png, dpi=300, bbox_inches="tight")
fig2.savefig(severity_pdf, bbox_inches="tight")

print(f"\nSaved: {heatmap_png}")
print(f"Saved: {heatmap_pdf}")
print(f"Saved: {severity_png}")
print(f"Saved: {severity_pdf}")
print(f"Saved: {folded_out}")
print(f"Saved: {metrics_out}")
