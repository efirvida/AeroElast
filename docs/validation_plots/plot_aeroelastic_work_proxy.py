"""Aeroelastic work and damping proxies per revolution — IEA 15 MW, rated.

Quantifies signed hysteresis loops between aerodynamic load fluctuations and the
flapwise tip motion, plus cycle-averaged quadrature proxies against flapwise
velocity. The goal is not to claim exact modal work, but to compare whether the
two structural formulations inject or extract flapwise-coupled energy-like
content with the exported rotor-equivalent signals.

Outputs
-------
docs/figures/fig_5_4_7h_hysteresis_work_loops.png
docs/figures/fig_5_4_7h_hysteresis_work_loops.pdf
docs/figures/fig_5_4_7i_cycle_work_proxy.png
docs/figures/fig_5_4_7i_cycle_work_proxy.pdf
docs/validation_data/generated/aeroelastic_work_proxy.csv
docs/validation_data/generated/aeroelastic_work_loops.csv
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
AZIMUTH_BINS_DEG = np.linspace(0.0, 360.0, 97)

COLORS = {
    "corot": "#c0392b",
    "inert": "#2980b9",
}

SIGNAL_COLUMNS = {
    "time": "Time [s]",
    "thrust": "Thrust [N]",
    "torque": "Torque [N.m]",
    "power": "Power [W]",
    "flapwise": "Tip Disp Y [m]",
    "resultante": "Tip Disp Mag [m]",
}


def load_case(csv_path: Path) -> tuple[pd.DataFrame | None, dict[str, float | str]]:
    cols = [
        SIGNAL_COLUMNS["time"],
        SIGNAL_COLUMNS["thrust"],
        SIGNAL_COLUMNS["torque"],
        SIGNAL_COLUMNS["power"],
        SIGNAL_COLUMNS["flapwise"],
        SIGNAL_COLUMNS["resultante"],
    ]
    df = pd.read_csv(csv_path)
    finite = np.asarray(np.isfinite(df[cols]).all(axis=1), dtype=bool)
    physical = np.ones(len(df), dtype=bool)
    physical &= df[SIGNAL_COLUMNS["resultante"]].abs().to_numpy() < 50.0
    physical &= df[SIGNAL_COLUMNS["thrust"]].abs().to_numpy() < 1.0e8
    physical &= df[SIGNAL_COLUMNS["torque"]].abs().to_numpy() < 1.0e8
    physical &= df[SIGNAL_COLUMNS["power"]].abs().to_numpy() < 1.0e8
    valid = finite & physical

    first_bad = int(np.argmax(~valid)) if np.any(~valid) else len(df)
    truncated = df.iloc[:first_bad].copy() if first_bad < len(df) else df.copy()
    t_end_valid = (
        float(truncated[SIGNAL_COLUMNS["time"]].iloc[-1]) if not truncated.empty else float("nan")
    )
    steady = truncated[
        (truncated[SIGNAL_COLUMNS["time"]] >= T_START)
        & (truncated[SIGNAL_COLUMNS["time"]] <= T_END)
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


def reconstruct_theta_rad(time_s: np.ndarray) -> np.ndarray:
    return OMEGA_RAD_S * (time_s - float(time_s[0]))


def fit_single_harmonic(
    signal: np.ndarray, azimuth_deg: np.ndarray, order: int
) -> dict[str, float]:
    psi = np.deg2rad(azimuth_deg)
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    design = np.column_stack([
        np.cos(order * psi),
        np.sin(order * psi),
    ])
    coeffs, *_ = np.linalg.lstsq(design, centered, rcond=None)
    a_cos, b_sin = coeffs
    amplitude = float(np.hypot(a_cos, b_sin))
    phase_deg = float(np.rad2deg(np.arctan2(b_sin, a_cos)))
    peak_azimuth_deg = float((phase_deg / order) % 360.0)
    return {
        "order": order,
        "amplitude": amplitude,
        "phase_deg": phase_deg,
        "peak_azimuth_deg": peak_azimuth_deg,
    }


def full_cycle_masks(theta_rad: np.ndarray) -> list[np.ndarray]:
    cycle_ids = np.floor((theta_rad - float(theta_rad[0])) / (2.0 * np.pi)).astype(int)
    masks: list[np.ndarray] = []
    for cycle_id in range(int(cycle_ids.min()), int(cycle_ids.max())):
        mask = cycle_ids == cycle_id
        if int(np.sum(mask)) < 8:
            continue
        theta_span = float(theta_rad[mask][-1] - theta_rad[mask][0])
        if theta_span < 1.8 * np.pi:
            continue
        masks.append(mask)
    return masks


def mean_cycle_loop(
    azimuth_deg: np.ndarray,
    torque_prime: np.ndarray,
    thrust_prime: np.ndarray,
    flapwise_prime: np.ndarray,
    flapwise_velocity: np.ndarray,
    azimuth_offset_deg: float,
) -> pd.DataFrame:
    azimuth_aligned = (azimuth_deg - azimuth_offset_deg) % 360.0
    bin_ids = np.digitize(azimuth_aligned, AZIMUTH_BINS_DEG, right=False) - 1
    bin_ids = np.clip(bin_ids, 0, len(AZIMUTH_BINS_DEG) - 2)
    centers = 0.5 * (AZIMUTH_BINS_DEG[:-1] + AZIMUTH_BINS_DEG[1:])

    records: list[dict[str, float]] = []
    for idx, center in enumerate(centers):
        mask = bin_ids == idx
        if not np.any(mask):
            records.append({
                "azimuth_deg": float(center),
                "torque_prime_MNm": float("nan"),
                "thrust_prime_MN": float("nan"),
                "flapwise_prime_m": float("nan"),
                "flapwise_velocity_m_s": float("nan"),
            })
            continue
        records.append({
            "azimuth_deg": float(center),
            "torque_prime_MNm": float(np.mean(torque_prime[mask]) / 1.0e6),
            "thrust_prime_MN": float(np.mean(thrust_prime[mask]) / 1.0e6),
            "flapwise_prime_m": float(np.mean(flapwise_prime[mask])),
            "flapwise_velocity_m_s": float(np.mean(flapwise_velocity[mask])),
        })
    return pd.DataFrame.from_records(records)


def summarize_case(df: pd.DataFrame) -> tuple[dict[str, float], pd.DataFrame]:
    t = df[SIGNAL_COLUMNS["time"]].to_numpy(dtype=float)
    azimuth_deg = reconstruct_azimuth_deg(t)
    theta_rad = reconstruct_theta_rad(t)
    thrust = df[SIGNAL_COLUMNS["thrust"]].to_numpy(dtype=float)
    torque = df[SIGNAL_COLUMNS["torque"]].to_numpy(dtype=float)
    flapwise = df[SIGNAL_COLUMNS["flapwise"]].to_numpy(dtype=float)

    thrust_prime = thrust - float(np.mean(thrust))
    torque_prime = torque - float(np.mean(torque))
    flapwise_prime = flapwise - float(np.mean(flapwise))
    flapwise_velocity = np.gradient(flapwise_prime, t)

    torque_1p = fit_single_harmonic(torque, azimuth_deg, 1)
    azimuth_offset_deg = torque_1p["peak_azimuth_deg"]
    cycle_masks = full_cycle_masks(theta_rad)

    qdy_proxy: list[float] = []
    qv_proxy: list[float] = []
    tdy_work: list[float] = []
    tv_power: list[float] = []
    rotor_work: list[float] = []

    for mask in cycle_masks:
        tm = t[mask]
        theta_m = theta_rad[mask]
        yy = flapwise_prime[mask]
        vv = flapwise_velocity[mask]
        qq = torque_prime[mask]
        tt = thrust_prime[mask]
        qq_abs = torque[mask]
        duration = float(tm[-1] - tm[0])
        qdy_proxy.append(float(np.trapezoid(qq, yy)))
        tdy_work.append(float(np.trapezoid(tt, yy)))
        qv_proxy.append(float(np.trapezoid(qq * vv, tm) / duration))
        tv_power.append(float(np.trapezoid(tt * vv, tm) / duration))
        rotor_work.append(float(np.trapezoid(qq_abs, theta_m)))

    summary = {
        "cycle_count": float(len(cycle_masks)),
        "flapwise_mean_m": float(np.mean(flapwise)),
        "flapwise_std_m": float(np.std(flapwise, ddof=1)),
        "torque_mean_MNm": float(np.mean(torque) / 1.0e6),
        "thrust_mean_MN": float(np.mean(thrust) / 1.0e6),
        "torque_1p_peak_deg": float(azimuth_offset_deg),
        "qdy_proxy_mean_MJm": float(np.mean(qdy_proxy) / 1.0e6),
        "qdy_proxy_std_MJm": float(np.std(qdy_proxy, ddof=1) / 1.0e6)
        if len(qdy_proxy) > 1
        else 0.0,
        "qv_proxy_mean_MWm": float(np.mean(qv_proxy) / 1.0e6),
        "qv_proxy_std_MWm": float(np.std(qv_proxy, ddof=1) / 1.0e6) if len(qv_proxy) > 1 else 0.0,
        "tdy_work_mean_MJ": float(np.mean(tdy_work) / 1.0e6),
        "tdy_work_std_MJ": float(np.std(tdy_work, ddof=1) / 1.0e6) if len(tdy_work) > 1 else 0.0,
        "tv_power_mean_MW": float(np.mean(tv_power) / 1.0e6),
        "tv_power_std_MW": float(np.std(tv_power, ddof=1) / 1.0e6) if len(tv_power) > 1 else 0.0,
        "rotor_work_mean_MJ": float(np.mean(rotor_work) / 1.0e6),
        "rotor_work_std_MJ": float(np.std(rotor_work, ddof=1) / 1.0e6)
        if len(rotor_work) > 1
        else 0.0,
    }

    loops = mean_cycle_loop(
        azimuth_deg,
        torque_prime,
        thrust_prime,
        flapwise_prime,
        flapwise_velocity,
        azimuth_offset_deg,
    )
    return summary, loops


summary_rows: list[dict[str, float | int | str]] = []
loop_rows: list[dict[str, float | int | str]] = []

for solver_name, root in ROOTS.items():
    for yaw in YAWS:
        csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
        if not csv_path.exists():
            summary_rows.append({"solver": solver_name, "yaw_deg": yaw, "status": "missing_csv"})
            continue
        steady_df, meta = load_case(csv_path)
        if steady_df is None:
            summary_rows.append({"solver": solver_name, "yaw_deg": yaw, **meta})
            continue

        summary, loops = summarize_case(steady_df)
        summary_rows.append({"solver": solver_name, "yaw_deg": yaw, "status": "ok", **summary})
        for record in loops.to_dict(orient="records"):
            loop_rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                **record,
            })

summary_df = pd.DataFrame(summary_rows)
loops_df = pd.DataFrame(loop_rows)

summary_out = DATA_DIR / "aeroelastic_work_proxy.csv"
loops_out = DATA_DIR / "aeroelastic_work_loops.csv"
summary_df.to_csv(summary_out, index=False)
loops_df.to_csv(loops_out, index=False)

print("\nAEROELASTIC WORK / DAMPING PROXY")
for solver_name in ("corot", "inert"):
    print(f"\nsolver={solver_name}")
    for yaw in YAWS:
        row = summary_df[(summary_df["solver"] == solver_name) & (summary_df["yaw_deg"] == yaw)]
        if row.empty:
            continue
        status = row.iloc[0].get("status", "missing")
        if status != "ok":
            print(f"  yaw={yaw:>2}: skipped ({status})")
            continue
        r = row.iloc[0]
        print(
            f"  yaw={yaw:>2}: QdY={r['qdy_proxy_mean_MJm']:+.3f} MJ·m, "
            f"TdY={r['tdy_work_mean_MJ']:+.3f} MJ, "
            f"<Qv>={r['qv_proxy_mean_MWm']:+.3f} MW·m, "
            f"<Tv>={r['tv_power_mean_MW']:+.3f} MW, "
            f"Wrev={r['rotor_work_mean_MJ']:.3f} MJ/rev"
        )

fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.0), sharex=False, sharey=False)
fig.suptitle(
    "Lazos de histéresis medios por vuelta en yaw=0°\nAzimut alineado al máximo 1P del torque",
    fontsize=11,
    y=1.02,
)

panel_specs = {
    (0, 0): (
        "corot",
        "torque_prime_MNm",
        "flapwise_prime_m",
        "(a) Corrotacional: Q' vs Y'",
        "Torque fluctuation [MN·m]",
    ),
    (0, 1): (
        "inert",
        "torque_prime_MNm",
        "flapwise_prime_m",
        "(b) Inercial: Q' vs Y'",
        "Torque fluctuation [MN·m]",
    ),
    (1, 0): (
        "corot",
        "thrust_prime_MN",
        "flapwise_prime_m",
        "(c) Corrotacional: T' vs Y'",
        "Thrust fluctuation [MN]",
    ),
    (1, 1): (
        "inert",
        "thrust_prime_MN",
        "flapwise_prime_m",
        "(d) Inercial: T' vs Y'",
        "Thrust fluctuation [MN]",
    ),
}

for (i, j), (solver_name, y_col, x_col, title, ylabel) in panel_specs.items():
    ax = axes[i, j]
    subset = loops_df[(loops_df["solver"] == solver_name) & (loops_df["yaw_deg"] == 0)].copy()
    ax.plot(subset[x_col], subset[y_col], color=COLORS[solver_name], lw=2.0)
    ax.scatter(
        subset[x_col].iloc[0], subset[y_col].iloc[0], color=COLORS[solver_name], s=28, zorder=3
    )
    arrow_idx = min(12, len(subset) - 2)
    ax.annotate(
        "",
        xy=(subset[x_col].iloc[arrow_idx + 1], subset[y_col].iloc[arrow_idx + 1]),
        xytext=(subset[x_col].iloc[arrow_idx], subset[y_col].iloc[arrow_idx]),
        arrowprops=dict(arrowstyle="->", color=COLORS[solver_name], lw=1.5),
    )
    metric_row = summary_df[
        (summary_df["solver"] == solver_name) & (summary_df["yaw_deg"] == 0)
    ].iloc[0]
    metric_text = (
        f"∮Q' dY = {metric_row['qdy_proxy_mean_MJm']:+.3f} MJ·m"
        if "torque" in y_col
        else f"∮T' dY = {metric_row['tdy_work_mean_MJ']:+.3f} MJ"
    )
    ax.text(0.03, 0.95, metric_text, transform=ax.transAxes, va="top", fontsize=8)
    ax.axhline(0.0, color="#777", lw=0.9, ls=":")
    ax.axvline(0.0, color="#777", lw=0.9, ls=":")
    ax.set_title(title)
    ax.set_xlabel("Flapwise fluctuation Y' [m]")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)

fig.tight_layout()
loops_png = OUT_DIR / "fig_5_4_7h_hysteresis_work_loops.png"
loops_pdf = OUT_DIR / "fig_5_4_7h_hysteresis_work_loops.pdf"
fig.savefig(loops_png, dpi=300, bbox_inches="tight")
fig.savefig(loops_pdf, bbox_inches="tight")

fig2, axes2 = plt.subplots(2, 2, figsize=(12.8, 8.0), sharex=True)
fig2.suptitle(
    "Proxies energéticos por vuelta para el acoplamiento flapwise\n"
    "Lazo firmado, cuadratura velocidad-carga y control de trabajo del rotor",
    fontsize=11,
    y=1.02,
)

ok_df = summary_df[summary_df["status"] == "ok"].copy()

summary_panels = [
    ((0, 0), "qdy_proxy_mean_MJm", "(a) Proxy de lazo ∮Q' dY", "MJ·m"),
    ((0, 1), "tdy_work_mean_MJ", "(b) Trabajo axial proxy ∮T' dY", "MJ"),
    ((1, 0), "qv_proxy_mean_MWm", "(c) Proxy cuadratura <Q'·Ydot>", "MW·m"),
    ((1, 1), "tv_power_mean_MW", "(d) Potencia axial proxy <T'·Ydot>", "MW"),
]

for (i, j), column, title, ylabel in summary_panels:
    ax = axes2[i, j]
    for solver_name in ("corot", "inert"):
        subset = ok_df[ok_df["solver"] == solver_name]
        ax.plot(
            subset["yaw_deg"],
            subset[column],
            marker="o" if solver_name == "corot" else "s",
            color=COLORS[solver_name],
            lw=1.8,
            label="Corrotacional" if solver_name == "corot" else "Inercial",
        )
    ax.axhline(0.0, color="#777", lw=0.9, ls=":")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    if i == 1:
        ax.set_xlabel("Yaw [deg]")

handles, labels = axes2[0, 0].get_legend_handles_labels()
unique = dict(zip(labels, handles))
axes2[0, 1].legend(unique.values(), unique.keys(), frameon=False, fontsize=8)

fig2.tight_layout()
summary_png = OUT_DIR / "fig_5_4_7i_cycle_work_proxy.png"
summary_pdf = OUT_DIR / "fig_5_4_7i_cycle_work_proxy.pdf"
fig2.savefig(summary_png, dpi=300, bbox_inches="tight")
fig2.savefig(summary_pdf, bbox_inches="tight")

print(f"\nSaved: {loops_png}")
print(f"Saved: {loops_pdf}")
print(f"Saved: {summary_png}")
print(f"Saved: {summary_pdf}")
print(f"Saved: {summary_out}")
print(f"Saved: {loops_out}")
