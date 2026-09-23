"""Aerodynamic-to-structural power exchange for the corotational yaw sweep.

This script uses checkpointed aerodynamic nodal forces and structural nodal
velocities to estimate the instantaneous power exchanged with the structural
motion channel,

    P_struct(t) = sum_i F_aero,i^global . v_i^global

where the stored VTU velocity field is rotated back from the corotating frame to
the inertial frame using the checkpoint angle ``theta``.

This is a direct aero-structural power channel, not an exact elastic-energy
budget. Splitting that exchange into elastic storage, kinetic storage, damping,
and returned shaft work would additionally require reconstructing the structural
energy functional offline.

Outputs
-------
docs/figures/fig_5_4_7l_structural_power_exchange.png
docs/figures/fig_5_4_7l_structural_power_exchange.pdf
docs/figures/fig_5_4_7m_structural_energy_channel.png
docs/figures/fig_5_4_7m_structural_energy_channel.pdf
docs/validation_data/generated/structural_power_exchange_timeseries.csv
docs/validation_data/generated/structural_power_exchange_summary.csv

By default, existing CSV outputs are reused to regenerate the figures quickly.
Set ``AEROELAST_REBUILD_STRUCTURAL_POWER=1`` to rebuild the CSV files from the
VTU checkpoints.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np
import pandas as pd

BASE = Path("/scratch/leahk/eduardo.donestevez")
ROOT = BASE / "frontiersin_results_corotational_100s"
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
CACHE_TIMESERIES = DATA_DIR / "structural_power_exchange_timeseries.csv"
CACHE_SUMMARY = DATA_DIR / "structural_power_exchange_summary.csv"
REBUILD_FROM_CHECKPOINTS = os.environ.get("AEROELAST_REBUILD_STRUCTURAL_POWER", "0") == "1"

YAWS = [0, 10, 20, 30, 40]
TIME_SERIES_YAWS = [0, 20, 40]
OMEGA_RPM = 7.55
ROTOR_PERIOD_S = 60.0 / OMEGA_RPM
LAST_REVOLUTIONS = 3

# The frontiersin corotational campaign uses a horizontal-axis turbine aligned
# with global Y. ``theta`` is stored per checkpoint, but the axis itself is not.
ROTATION_AXIS = np.array([0.0, 1.0, 0.0], dtype=float)
ROTATION_AXIS /= np.linalg.norm(ROTATION_AXIS)

COLORS = {
    0: "#1f3c88",
    10: "#3a73b4",
    20: "#2e8b57",
    30: "#d98e04",
    40: "#b33a3a",
}


def rotation_matrix(theta_rad: float) -> np.ndarray:
    axis = ROTATION_AXIS
    k = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ],
        dtype=float,
    )
    k2 = k @ k
    c = float(np.cos(theta_rad))
    s = float(np.sin(theta_rad))
    return np.eye(3) + s * k + (1.0 - c) * k2


def to_inertial(vec_local: np.ndarray, theta_rad: float) -> np.ndarray:
    return np.asarray(vec_local, dtype=float) @ rotation_matrix(theta_rad).T


def list_time_dirs(case_dir: Path) -> list[tuple[float, Path]]:
    time_dirs: list[tuple[float, Path]] = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            time_value = float(child.name)
        except ValueError:
            continue
        if not (child / "state.npz").exists() or not (child / "fields.vtu").exists():
            continue
        time_dirs.append((time_value, child))
    time_dirs.sort(key=lambda item: item[0])
    return time_dirs


def load_case_timeseries(case_dir: Path) -> pd.DataFrame:
    perf = pd.read_csv(case_dir / "rotor_performance.csv")
    perf_time = perf["Time [s]"].to_numpy(dtype=float)
    perf_power = perf["Aero Power [W]"].to_numpy(dtype=float)

    time_dirs = list_time_dirs(case_dir)
    if not time_dirs:
        return pd.DataFrame()

    end_time = time_dirs[-1][0]
    start_time = max(0.0, end_time - LAST_REVOLUTIONS * ROTOR_PERIOD_S)

    rows: list[dict[str, float]] = []
    for time_value, folder in time_dirs:
        if time_value < start_time:
            continue
        state = np.load(folder / "state.npz")
        theta_rad = float(state["theta"])
        mesh = meshio.read(folder / "fields.vtu")

        if "F_AERO" not in mesh.point_data or "VEL" not in mesh.point_data:
            continue

        force_global = np.asarray(mesh.point_data["F_AERO"], dtype=float)
        velocity_local = np.asarray(mesh.point_data["VEL"], dtype=float)
        velocity_global = to_inertial(velocity_local, theta_rad)
        structural_power_w = float(np.sum(force_global * velocity_global))
        aero_power_w = float(np.interp(time_value, perf_time, perf_power))

        rows.append({
            "time_s": time_value,
            "theta_rad": theta_rad,
            "aero_power_W": aero_power_w,
            "structural_power_W": structural_power_w,
        })

    return pd.DataFrame.from_records(rows)


def last_revolution_windows(end_time: float) -> list[tuple[int, float, float]]:
    windows: list[tuple[int, float, float]] = []
    for rev_index in range(LAST_REVOLUTIONS, 0, -1):
        start = end_time - rev_index * ROTOR_PERIOD_S
        end = end_time - (rev_index - 1) * ROTOR_PERIOD_S
        windows.append((LAST_REVOLUTIONS - rev_index + 1, start, end))
    return windows


def integrate_window(
    df: pd.DataFrame, yaw_deg: int, window_id: int, start: float, end: float
) -> dict[str, float | int]:
    window = df[(df["time_s"] >= start) & (df["time_s"] <= end)].copy()
    if len(window) < 2:
        return {
            "yaw_deg": yaw_deg,
            "window_id": window_id,
            "t_start_s": start,
            "t_end_s": end,
            "samples": int(len(window)),
            "status": "insufficient_samples",
        }

    time_s = window["time_s"].to_numpy(dtype=float)
    aero_power_w = window["aero_power_W"].to_numpy(dtype=float)
    structural_power_w = window["structural_power_W"].to_numpy(dtype=float)

    aero_energy_j = float(np.trapezoid(aero_power_w, time_s))
    into_structure_j = float(np.trapezoid(np.clip(structural_power_w, 0.0, None), time_s))
    returned_j = float(np.trapezoid(np.clip(-structural_power_w, 0.0, None), time_s))
    net_j = float(np.trapezoid(structural_power_w, time_s))
    exchange_j = 0.5 * (into_structure_j + returned_j)

    return {
        "yaw_deg": yaw_deg,
        "window_id": window_id,
        "t_start_s": float(time_s[0]),
        "t_end_s": float(time_s[-1]),
        "samples": int(len(window)),
        "status": "ok",
        "aero_energy_MJ": aero_energy_j / 1.0e6,
        "into_structure_kJ": into_structure_j / 1.0e3,
        "returned_to_aero_kJ": returned_j / 1.0e3,
        "net_structural_kJ": net_j / 1.0e3,
        "exchange_kJ": exchange_j / 1.0e3,
        "into_structure_pct_of_aero": 100.0 * into_structure_j / aero_energy_j,
        "returned_to_aero_pct_of_aero": 100.0 * returned_j / aero_energy_j,
        "net_structural_pct_of_aero": 100.0 * net_j / aero_energy_j,
        "exchange_pct_of_aero": 100.0 * exchange_j / aero_energy_j,
        "return_ratio_pct": 100.0 * returned_j / into_structure_j
        if into_structure_j > 0.0
        else float("nan"),
        "structural_power_rms_kW": float(np.sqrt(np.mean(structural_power_w**2)) / 1.0e3),
        "structural_power_mean_kW": float(np.mean(structural_power_w) / 1.0e3),
        "structural_power_min_kW": float(np.min(structural_power_w) / 1.0e3),
        "structural_power_max_kW": float(np.max(structural_power_w) / 1.0e3),
    }


def summarize_windows(window_df: pd.DataFrame) -> pd.DataFrame:
    ok = window_df[window_df["status"] == "ok"].copy()
    summary_rows: list[dict[str, float | int]] = []
    metric_names = [
        "aero_energy_MJ",
        "into_structure_kJ",
        "returned_to_aero_kJ",
        "net_structural_kJ",
        "exchange_kJ",
        "into_structure_pct_of_aero",
        "returned_to_aero_pct_of_aero",
        "net_structural_pct_of_aero",
        "exchange_pct_of_aero",
        "return_ratio_pct",
        "structural_power_rms_kW",
        "structural_power_mean_kW",
        "structural_power_min_kW",
        "structural_power_max_kW",
    ]
    for yaw_deg, group in ok.groupby("yaw_deg"):
        row: dict[str, float | int] = {
            "yaw_deg": int(yaw_deg),
            "window_count": int(len(group)),
        }
        for metric in metric_names:
            row[f"{metric}_mean"] = float(group[metric].mean())
            row[f"{metric}_std"] = float(group[metric].std(ddof=1)) if len(group) > 1 else 0.0
        summary_rows.append(row)
    return pd.DataFrame.from_records(summary_rows).sort_values("yaw_deg").reset_index(drop=True)


def plot_time_series(timeseries_map: dict[int, pd.DataFrame]) -> None:
    available_yaws = [yaw_deg for yaw_deg in TIME_SERIES_YAWS if yaw_deg in timeseries_map]
    if not available_yaws:
        return
    fig, axes = plt.subplots(len(available_yaws), 1, figsize=(10.5, 8.8), sharex=False)
    axes = np.atleast_1d(axes)
    for ax, yaw_deg in zip(axes, available_yaws):
        df = timeseries_map[yaw_deg].copy()
        end_time = float(df["time_s"].max())
        start_time = end_time - ROTOR_PERIOD_S
        final_rev = df[df["time_s"] >= start_time].copy()
        time_local = final_rev["time_s"].to_numpy(dtype=float) - start_time
        time_s = final_rev["time_s"].to_numpy(dtype=float)
        aero_power_w = final_rev["aero_power_W"].to_numpy(dtype=float)
        structural_power_w = final_rev["structural_power_W"].to_numpy(dtype=float)
        aero_power_mean_w = float(np.mean(aero_power_w))
        structural_power_pct = 100.0 * structural_power_w / aero_power_mean_w
        aero_energy_j = float(np.trapezoid(aero_power_w, time_s))
        into_pct = (
            100.0
            * float(np.trapezoid(np.clip(structural_power_w, 0.0, None), time_s))
            / aero_energy_j
        )
        returned_pct = (
            100.0
            * float(np.trapezoid(np.clip(-structural_power_w, 0.0, None), time_s))
            / aero_energy_j
        )
        net_pct = into_pct - returned_pct

        ax.fill_between(
            time_local,
            0.0,
            structural_power_pct,
            where=structural_power_pct >= 0.0,
            color="#287271",
            alpha=0.32,
            interpolate=True,
            label="Area +: hacia la estructura",
        )
        ax.fill_between(
            time_local,
            0.0,
            structural_power_pct,
            where=structural_power_pct < 0.0,
            color="#b7410e",
            alpha=0.32,
            interpolate=True,
            label="Area -: devuelta al canal aero",
        )
        ax.plot(
            time_local,
            structural_power_pct,
            color=COLORS[yaw_deg],
            lw=2.0,
            label="$P_{str}$ normalizada",
        )
        ax.axhline(0.0, color="#777777", lw=0.9, ls="--")
        limit = max(5.0, 1.18 * float(np.max(np.abs(structural_power_pct))))
        ax.set_ylim(-limit, limit)
        ax.set_ylabel("$P_{str}$ / $\\bar{P}_{aero}$ [%]")
        ax.set_title(f"Yaw = {yaw_deg} deg")
        ax.text(
            0.015,
            0.06,
            f"E+={into_pct:.2f}%  E-={returned_pct:.2f}%  net={net_pct:.2f}%",
            transform=ax.transAxes,
            fontsize=9.5,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#cccccc",
                "alpha": 0.9,
            },
        )
        ax.grid(alpha=0.25)
        if yaw_deg == available_yaws[0]:
            ax.legend(loc="upper right", frameon=False)
    axes[-1].set_xlabel("Tiempo dentro de la ultima vuelta [s]")
    fig.suptitle(
        "Intercambio aero-estructural instantaneo normalizado\n"
        "Areas positivas y negativas dentro de la ultima vuelta disponible"
    )
    fig.tight_layout(h_pad=2.6)
    fig.savefig(OUT_DIR / "fig_5_4_7l_structural_power_exchange.png", dpi=220)
    fig.savefig(OUT_DIR / "fig_5_4_7l_structural_power_exchange.pdf")
    plt.close(fig)


def plot_summary(summary_df: pd.DataFrame) -> None:
    yaw = summary_df["yaw_deg"].to_numpy(dtype=float)
    into_pct = summary_df["into_structure_pct_of_aero_mean"].to_numpy(dtype=float)
    into_std = summary_df["into_structure_pct_of_aero_std"].to_numpy(dtype=float)
    returned_pct = summary_df["returned_to_aero_pct_of_aero_mean"].to_numpy(dtype=float)
    returned_std = summary_df["returned_to_aero_pct_of_aero_std"].to_numpy(dtype=float)
    net_pct = summary_df["net_structural_pct_of_aero_mean"].to_numpy(dtype=float)
    net_std = summary_df["net_structural_pct_of_aero_std"].to_numpy(dtype=float)
    released_stored_pct = returned_pct - into_pct

    fig, axes = plt.subplots(2, 1, figsize=(9.8, 8.2), sharex=True)

    width = 3.0
    axes[0].bar(
        yaw - width / 2.0, into_pct, width=width, color="#287271", label="$E_+$ hacia la estructura"
    )
    axes[0].bar(
        yaw + width / 2.0,
        -returned_pct,
        width=width,
        color="#b7410e",
        label="-$E_-$ devuelta al canal aero",
    )
    axes[0].errorbar(
        yaw - width / 2.0, into_pct, yerr=into_std, fmt="none", ecolor="#124c4b", capsize=4
    )
    axes[0].errorbar(
        yaw + width / 2.0,
        -returned_pct,
        yerr=returned_std,
        fmt="none",
        ecolor="#7c2807",
        capsize=4,
    )
    axes[0].errorbar(
        yaw,
        net_pct,
        yerr=net_std,
        color="#222222",
        marker="D",
        lw=1.8,
        capsize=4,
        label="$E_{net}$",
    )
    axes[0].axhline(0.0, color="#777777", lw=0.9, ls="--")
    axes[0].set_ylabel("Energia / $E_{aero}$ [%]")
    axes[0].set_title("Balance firmado del canal estructural en las ultimas 3 vueltas")
    axes[0].grid(alpha=0.25, axis="y")
    axes[0].legend(frameon=False, ncol=3, loc="upper center")

    lower_band = released_stored_pct - net_std
    upper_band = released_stored_pct + net_std
    axes[1].fill_between(
        yaw,
        lower_band,
        upper_band,
        color="#bbbbbb",
        alpha=0.35,
        label="Variabilidad entre vueltas (±1σ)",
        zorder=1,
    )
    axes[1].plot(
        yaw,
        released_stored_pct,
        color="#222222",
        marker="o",
        lw=2.0,
        label="Media últimas 3 vueltas",
        zorder=3,
    )
    for yaw_value, release_value in zip(yaw, released_stored_pct):
        axes[1].text(
            yaw_value,
            release_value + 0.025,
            f"{release_value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            zorder=4,
        )
    axes[1].axhline(0.0, color="#777777", lw=0.9, ls="--")
    axes[1].set_ylim(float(np.min(lower_band)) - 0.04, float(np.max(upper_band)) + 0.08)
    axes[1].set_xlabel("Angulo de yaw [deg]")
    axes[1].set_ylabel("$(E_- - E_+) / E_{aero}$ [%]")
    axes[1].set_title(
        "Energia almacenada liberada dentro de la ventana: escala real sobre $E_{aero}$", pad=14
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_5_4_7m_structural_energy_channel.png", dpi=220)
    fig.savefig(OUT_DIR / "fig_5_4_7m_structural_energy_channel.pdf")
    plt.close(fig)


def build_from_checkpoints() -> tuple[pd.DataFrame, pd.DataFrame, dict[int, pd.DataFrame]]:
    timeseries_rows: list[pd.DataFrame] = []
    window_rows: list[dict[str, float | int]] = []
    timeseries_map: dict[int, pd.DataFrame] = {}

    for yaw_deg in YAWS:
        case_dir = ROOT / f"yaw_{yaw_deg}" / "corotational"
        if not case_dir.exists():
            continue
        case_df = load_case_timeseries(case_dir)
        if case_df.empty:
            continue
        case_df.insert(0, "yaw_deg", yaw_deg)
        timeseries_rows.append(case_df)
        timeseries_map[yaw_deg] = case_df

        end_time = float(case_df["time_s"].max())
        for window_id, start_time, end_time_window in last_revolution_windows(end_time):
            window_rows.append(
                integrate_window(case_df, yaw_deg, window_id, start_time, end_time_window)
            )

    timeseries_df = pd.concat(timeseries_rows, ignore_index=True)
    window_df = pd.DataFrame.from_records(window_rows)
    summary_df = summarize_windows(window_df)

    timeseries_df.to_csv(CACHE_TIMESERIES, index=False)
    summary_df.to_csv(CACHE_SUMMARY, index=False)
    return timeseries_df, summary_df, timeseries_map


def load_cached_results() -> tuple[pd.DataFrame, pd.DataFrame, dict[int, pd.DataFrame]]:
    timeseries_df = pd.read_csv(CACHE_TIMESERIES)
    summary_df = pd.read_csv(CACHE_SUMMARY)
    timeseries_map = {
        yaw_deg: timeseries_df[timeseries_df["yaw_deg"] == yaw_deg].copy()
        for yaw_deg in TIME_SERIES_YAWS
        if not timeseries_df[timeseries_df["yaw_deg"] == yaw_deg].empty
    }
    return timeseries_df, summary_df, timeseries_map


def main() -> None:
    use_cache = (
        CACHE_TIMESERIES.exists() and CACHE_SUMMARY.exists() and not REBUILD_FROM_CHECKPOINTS
    )
    if use_cache:
        _, summary_df, timeseries_map = load_cached_results()
        source = "cached CSV"
    else:
        _, summary_df, timeseries_map = build_from_checkpoints()
        source = "VTU checkpoints"

    plot_time_series(timeseries_map)
    plot_summary(summary_df)

    print("AERO-STRUCTURAL POWER EXCHANGE")
    print(f"source: {source}")
    print()
    for _, row in summary_df.iterrows():
        yaw_deg = int(row["yaw_deg"])
        print(
            f"yaw={yaw_deg:>2}: "
            f"into={row['into_structure_pct_of_aero_mean']:.3f}% +/- {row['into_structure_pct_of_aero_std']:.3f}% of aero energy, "
            f"returned={row['returned_to_aero_pct_of_aero_mean']:.3f}% +/- {row['returned_to_aero_pct_of_aero_std']:.3f}%, "
            f"net={row['net_structural_pct_of_aero_mean']:.3f}% +/- {row['net_structural_pct_of_aero_std']:.3f}%"
        )

    print()
    print(f"Saved: {OUT_DIR / 'fig_5_4_7l_structural_power_exchange.png'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7l_structural_power_exchange.pdf'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7m_structural_energy_channel.png'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7m_structural_energy_channel.pdf'}")
    print(f"Saved: {CACHE_TIMESERIES}")
    print(f"Saved: {CACHE_SUMMARY}")


if __name__ == "__main__":
    main()
