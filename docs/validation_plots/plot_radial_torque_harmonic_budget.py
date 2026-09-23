"""Radial harmonic budget of torque from sectional BEM loads.

This postprocess locates where the single-blade sectional torque modulation is
generated along the blade span. It reconstructs a rotor-equivalent sectional
torque from the corotational rated yaw campaign using

    Q_band(t) = 3 * integral_band T_p(r, t) * r dr

and decomposes each radial band into 0P, 1P, 2P, and 3P contributions. The
factor 3 converts the exported single-blade sectional BEM loads into the rotor
equivalent convention used by ``bem_report.csv`` for the mean torque. Under yaw,
the global rotor harmonics additionally include multi-blade phase cancellation,
so the sectional harmonic budget is interpreted as local blade activity and is
shown together with a global BEM closure diagnostic.

Outputs
-------
docs/figures/fig_5_4_7p_radial_torque_harmonic_budget.png
docs/figures/fig_5_4_7p_radial_torque_harmonic_budget.pdf
docs/validation_data/generated/radial_torque_harmonic_budget.csv
docs/validation_data/generated/radial_torque_global_closure.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s"
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
DT_SAMPLE = 0.05
ROOT_CUTOFF_R_R = 0.16
OMEGA_RPM = 7.55
OMEGA_RAD_S = OMEGA_RPM * 2.0 * np.pi / 60.0
HARMONIC_ORDERS = [1, 2, 3]

BANDS = [
    {"name": "inboard", "label": "Inboard root-0.40R", "min": 0.0, "max": 0.40},
    {"name": "midspan", "label": "Midspan 0.40-0.70R", "min": 0.40, "max": 0.70},
    {"name": "outboard", "label": "Outboard 0.70-1.00R", "min": 0.70, "max": 1.01},
]

COLORS = {
    "inboard": "#3a73b4",
    "midspan": "#2e8b57",
    "outboard": "#d98e04",
    "1P": "#1f3c88",
    "2P": "#b7410e",
    "3P": "#6f4aa8",
}


def list_time_dirs(root: Path) -> list[Path]:
    dirs = [path for path in root.iterdir() if path.is_dir()]
    numeric_dirs: list[tuple[float, Path]] = []
    for path in dirs:
        try:
            numeric_dirs.append((float(path.name), path))
        except ValueError:
            continue
    numeric_dirs.sort(key=lambda item: item[0])
    return [path for _, path in numeric_dirs]


def reconstruct_azimuth_deg(time_s: np.ndarray) -> np.ndarray:
    return np.rad2deg(OMEGA_RAD_S * (time_s - float(time_s[0]))) % 360.0


def harmonic_coefficients(
    signal: np.ndarray, azimuth_deg: np.ndarray, order: int
) -> dict[str, float]:
    psi = np.deg2rad(np.asarray(azimuth_deg, dtype=float))
    centered = np.asarray(signal, dtype=float) - float(np.mean(signal))
    design = np.column_stack([np.cos(order * psi), np.sin(order * psi)])
    coeffs, *_ = np.linalg.lstsq(design, centered, rcond=None)
    a_cos, b_sin = coeffs
    component = a_cos * np.cos(order * psi) + b_sin * np.sin(order * psi)
    variance = float(np.var(centered, ddof=0))
    phase_deg = float(np.rad2deg(np.arctan2(b_sin, a_cos)))
    return {
        "vector_a": float(a_cos),
        "vector_b": float(b_sin),
        "a_cos": float(a_cos),
        "b_sin": float(b_sin),
        "amplitude": float(np.hypot(a_cos, b_sin)),
        "phase_deg": phase_deg,
        "peak_azimuth_deg": float((phase_deg / order) % 360.0),
        "variance_fraction": float(np.var(component, ddof=0) / variance) if variance > 0.0 else 0.0,
    }


def integrate_band_trapezoid(
    r_m: np.ndarray,
    r_r: np.ndarray,
    integrand: np.ndarray,
    band_min: float,
    band_max: float,
) -> float:
    rotor_radius_m = float(np.max(r_m))
    start_m = max(float(np.min(r_m)), band_min * rotor_radius_m)
    end_m = min(float(np.max(r_m)), band_max * rotor_radius_m)
    if end_m <= start_m:
        return 0.0

    mask = (r_r > band_min) & (r_r < band_max)
    band_r = np.concatenate([[start_m], r_m[mask], [end_m]])
    band_y = np.concatenate([
        [float(np.interp(start_m, r_m, integrand))],
        integrand[mask],
        [float(np.interp(end_m, r_m, integrand))],
    ])
    order = np.argsort(band_r)
    return float(np.trapezoid(band_y[order], band_r[order]))


def load_case_band_torque(fluid_root: Path) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    last_t = -np.inf

    for time_dir in list_time_dirs(fluid_root):
        try:
            time_s = float(time_dir.name)
        except ValueError:
            continue
        if time_s < T_START or time_s > T_END:
            continue
        if time_s - last_t < DT_SAMPLE:
            continue

        csv_path = time_dir / "bem_sectional.csv"
        if not csv_path.exists():
            continue

        df = pd.read_csv(csv_path)
        r_m = df["r[m]"].to_numpy(dtype=float)
        r_r = r_m / float(np.max(r_m))
        tp_n_m = df["Tp[N/m]"].to_numpy(dtype=float)
        integrand = tp_n_m * r_m

        row: dict[str, float] = {"time_s": time_s}
        total = 0.0
        for band in BANDS:
            value = float(
                3.0
                * integrate_band_trapezoid(
                    r_m,
                    r_r,
                    integrand,
                    float(band["min"]),
                    float(band["max"]),
                )
                / 1.0e6
            )
            row[str(band["name"])] = value
            total += value
        row["total"] = total
        rows.append(row)
        last_t = time_s

    if not rows:
        raise RuntimeError(f"No bem_sectional.csv samples found in {fluid_root}")

    return pd.DataFrame.from_records(rows)


def load_global_torque_at_times(fluid_root: Path, time_s: np.ndarray) -> np.ndarray:
    report = pd.read_csv(fluid_root / "bem_report.csv")
    return np.interp(
        time_s,
        report["Time [s]"].to_numpy(dtype=float),
        report["Torque [N.m]"].to_numpy(dtype=float) / 1.0e6,
    )


def analyze_case(yaw: int, band_df: pd.DataFrame) -> list[dict[str, float | int | str]]:
    azimuth_deg = reconstruct_azimuth_deg(band_df["time_s"].to_numpy(dtype=float))
    rows: list[dict[str, float | int | str]] = []

    total_mean = float(band_df["total"].mean())
    for band in BANDS:
        band_name = str(band["name"])
        band_signal = band_df[band_name].to_numpy(dtype=float)
        rows.append({
            "yaw_deg": yaw,
            "band": band_name,
            "band_label": str(band["label"]),
            "quantity": "mean",
            "order": 0,
            "mean_torque_MNm": float(np.mean(band_signal)),
            "mean_share_pct": 100.0 * float(np.mean(band_signal)) / total_mean,
            "harmonic_amplitude_MNm": np.nan,
            "harmonic_activity_share_pct": np.nan,
            "coherent_harmonic_contribution_MNm": np.nan,
            "coherent_harmonic_share_pct": np.nan,
            "phase_deg": np.nan,
            "variance_fraction": np.nan,
            "total_harmonic_amplitude_MNm": np.nan,
        })

    for order in HARMONIC_ORDERS:
        total_fit = harmonic_coefficients(
            band_df["total"].to_numpy(dtype=float), azimuth_deg, order=order
        )
        band_fits = {
            str(band["name"]): harmonic_coefficients(
                band_df[str(band["name"])].to_numpy(dtype=float), azimuth_deg, order=order
            )
            for band in BANDS
        }
        amplitude_sum = sum(fit["amplitude"] for fit in band_fits.values())
        total_vector = np.array([total_fit["vector_a"], total_fit["vector_b"]], dtype=float)
        total_amplitude = float(total_fit["amplitude"])
        total_direction = total_vector / total_amplitude if total_amplitude > 0.0 else total_vector
        for band in BANDS:
            band_name = str(band["name"])
            fit = band_fits[band_name]
            band_vector = np.array([fit["vector_a"], fit["vector_b"]], dtype=float)
            coherent_projection = (
                float(np.dot(band_vector, total_direction)) if total_amplitude > 0.0 else np.nan
            )
            rows.append({
                "yaw_deg": yaw,
                "band": band_name,
                "band_label": str(band["label"]),
                "quantity": f"{order}P",
                "order": order,
                "mean_torque_MNm": np.nan,
                "mean_share_pct": np.nan,
                "harmonic_amplitude_MNm": fit["amplitude"],
                "harmonic_activity_share_pct": 100.0 * fit["amplitude"] / amplitude_sum
                if amplitude_sum > 0.0
                else np.nan,
                "coherent_harmonic_contribution_MNm": coherent_projection,
                "coherent_harmonic_share_pct": 100.0 * coherent_projection / total_amplitude
                if total_amplitude > 0.0
                else np.nan,
                "phase_deg": fit["phase_deg"],
                "variance_fraction": fit["variance_fraction"],
                "total_harmonic_amplitude_MNm": total_fit["amplitude"],
            })

    return rows


def analyze_global_closure(
    yaw: int, band_df: pd.DataFrame, global_torque_mnm: np.ndarray
) -> list[dict[str, float | int | str]]:
    time_s = band_df["time_s"].to_numpy(dtype=float)
    azimuth_deg = reconstruct_azimuth_deg(time_s)
    reconstructed = band_df["total"].to_numpy(dtype=float)

    rows: list[dict[str, float | int | str]] = [
        {
            "yaw_deg": yaw,
            "quantity": "mean",
            "order": 0,
            "sectional_equivalent_MNm": float(np.mean(reconstructed)),
            "global_bem_MNm": float(np.mean(global_torque_mnm)),
            "ratio_sectional_to_global": float(np.mean(reconstructed) / np.mean(global_torque_mnm)),
            "correlation_sectional_global": float(
                np.corrcoef(reconstructed, global_torque_mnm)[0, 1]
            ),
        }
    ]

    for order in HARMONIC_ORDERS:
        sectional_fit = harmonic_coefficients(reconstructed, azimuth_deg, order=order)
        global_fit = harmonic_coefficients(global_torque_mnm, azimuth_deg, order=order)
        rows.append({
            "yaw_deg": yaw,
            "quantity": f"{order}P",
            "order": order,
            "sectional_equivalent_MNm": sectional_fit["amplitude"],
            "global_bem_MNm": global_fit["amplitude"],
            "ratio_sectional_to_global": sectional_fit["amplitude"] / global_fit["amplitude"]
            if global_fit["amplitude"] > 0.0
            else np.nan,
            "correlation_sectional_global": np.nan,
        })

    return rows


def plot_budget(summary_df: pd.DataFrame, closure_df: pd.DataFrame) -> None:
    yaw_values = np.asarray(YAWS, dtype=float)
    figure, axes = plt.subplots(2, 2, figsize=(12.4, 8.4), sharex=False)
    figure.suptitle(
        "Presupuesto radial del torque armónico: campaña corrotacional",
        fontsize=12,
        y=1.01,
    )

    ax_harm = axes[0, 0]
    for order in [1, 2]:
        subset = summary_df[
            (summary_df["quantity"] == f"{order}P") & (summary_df["band"] == "inboard")
        ]
        closure_subset = closure_df[closure_df["quantity"] == f"{order}P"].sort_values("yaw_deg")
        ax_harm.plot(
            subset["yaw_deg"],
            subset["total_harmonic_amplitude_MNm"],
            marker="o",
            lw=2.0,
            color=COLORS[f"{order}P"],
            label=f"{order}P seccional equiv.",
        )
        ax_harm.plot(
            closure_subset["yaw_deg"],
            closure_subset["global_bem_MNm"],
            marker="s",
            lw=1.7,
            ls="--",
            color=COLORS[f"{order}P"],
            alpha=0.85,
            label=f"{order}P rotor global",
        )
    ax_harm.set_title("Amplitud armónica: pala seccional vs rotor global")
    ax_harm.set_xlabel("Yaw [deg]")
    ax_harm.set_ylabel("Amplitud [MNm]")
    ax_harm.grid(alpha=0.25)
    ax_harm.legend(frameon=False)

    panels = [
        (axes[0, 1], "mean", "Torque medio por banda", "Participación media [%]"),
        (axes[1, 0], "1P", "Origen radial de la actividad 1P", "Participación de amplitud [%]"),
        (axes[1, 1], "2P", "Origen radial de la actividad 2P", "Participación de amplitud [%]"),
    ]

    for ax, quantity, title, ylabel in panels:
        bottom = np.zeros(len(yaw_values), dtype=float)
        for band in BANDS:
            band_name = str(band["name"])
            subset = summary_df[
                (summary_df["quantity"] == quantity) & (summary_df["band"] == band_name)
            ].sort_values("yaw_deg")
            values = (
                subset["mean_share_pct"].to_numpy(dtype=float)
                if quantity == "mean"
                else subset["harmonic_activity_share_pct"].to_numpy(dtype=float)
            )
            ax.bar(
                yaw_values,
                values,
                bottom=bottom,
                width=7.2,
                color=COLORS[band_name],
                label=str(band["label"]),
            )
            bottom += values
        ax.set_title(title)
        ax.set_xlabel("Yaw [deg]")
        ax.set_ylabel(ylabel)
        ax.set_xticks(YAWS)
        ax.set_ylim(0.0, 105.0)
        ax.grid(axis="y", alpha=0.25)

    axes[0, 1].legend(frameon=False, fontsize=8.0, loc="upper left")
    figure.tight_layout()
    figure.savefig(
        OUT_DIR / "fig_5_4_7p_radial_torque_harmonic_budget.png",
        dpi=280,
        bbox_inches="tight",
    )
    figure.savefig(
        OUT_DIR / "fig_5_4_7p_radial_torque_harmonic_budget.pdf",
        bbox_inches="tight",
    )
    plt.close(figure)


def main() -> None:
    rows: list[dict[str, float | int | str]] = []
    closure_rows: list[dict[str, float | int | str]] = []
    for yaw in YAWS:
        fluid_root = CAMPAIGN_ROOT / f"yaw_{yaw}" / "fluid"
        band_df = load_case_band_torque(fluid_root)
        global_torque_mnm = load_global_torque_at_times(
            fluid_root, band_df["time_s"].to_numpy(dtype=float)
        )
        rows.extend(analyze_case(yaw, band_df))
        closure_rows.extend(analyze_global_closure(yaw, band_df, global_torque_mnm))

    summary_df = pd.DataFrame.from_records(rows)
    closure_df = pd.DataFrame.from_records(closure_rows)
    out_csv = DATA_DIR / "radial_torque_harmonic_budget.csv"
    closure_csv = DATA_DIR / "radial_torque_global_closure.csv"
    summary_df.to_csv(out_csv, index=False)
    closure_df.to_csv(closure_csv, index=False)
    plot_budget(summary_df, closure_df)

    print("\nRADIAL TORQUE HARMONIC BUDGET")
    for yaw in YAWS:
        yaw_df = summary_df[summary_df["yaw_deg"] == yaw]
        mean_row = yaw_df[yaw_df["quantity"] == "mean"].sort_values("mean_share_pct")
        one_p = yaw_df[yaw_df["quantity"] == "1P"].sort_values(
            "harmonic_activity_share_pct", ascending=False
        )
        two_p = yaw_df[yaw_df["quantity"] == "2P"].sort_values(
            "harmonic_activity_share_pct", ascending=False
        )
        amp_1p = float(one_p["total_harmonic_amplitude_MNm"].iloc[0])
        amp_2p = float(two_p["total_harmonic_amplitude_MNm"].iloc[0])
        closure_yaw = closure_df[closure_df["yaw_deg"] == yaw]
        mean_ratio = float(
            closure_yaw[closure_yaw["quantity"] == "mean"]["ratio_sectional_to_global"].iloc[0]
        )
        one_p_ratio = float(
            closure_yaw[closure_yaw["quantity"] == "1P"]["ratio_sectional_to_global"].iloc[0]
        )
        print(
            f"yaw={yaw:>2}: mean dominant={mean_row.iloc[-1]['band']} "
            f"({mean_row.iloc[-1]['mean_share_pct']:.1f}%), "
            f"1P amp={amp_1p:.3f} MNm dominated by {one_p.iloc[0]['band']} "
            f"({one_p.iloc[0]['harmonic_activity_share_pct']:.1f}%), "
            f"2P amp={amp_2p:.3f} MNm dominated by {two_p.iloc[0]['band']} "
            f"({two_p.iloc[0]['harmonic_activity_share_pct']:.1f}%), "
            f"mean closure={mean_ratio:.3f}, 1P sectional/global={one_p_ratio:.2f}"
        )

    print(f"\nSaved: {OUT_DIR / 'fig_5_4_7p_radial_torque_harmonic_budget.png'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7p_radial_torque_harmonic_budget.pdf'}")
    print(f"Saved: {out_csv}")
    print(f"Saved: {closure_csv}")


if __name__ == "__main__":
    main()
