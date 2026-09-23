"""Parked/flag V50 response and frequency characterization.

This postprocess documents the ``bem_90_50_S`` simulation used in the
validation report section 5.8. It produces one deformation-vs-time figure with
the NREL DLC 6.x parked deflection scale, plus one frequency-domain figure for
the exported tip displacement channels.

Outputs
-------
docs/figures/fig_5_8_1_parked_v50_deformation.png
docs/figures/fig_5_8_1_parked_v50_deformation.pdf
docs/figures/fig_5_8_2_parked_v50_frequency.png
docs/figures/fig_5_8_2_parked_v50_frequency.pdf
docs/validation_data/generated/parked_flag_v50_response_summary.csv
docs/validation_data/generated/parked_flag_v50_frequency_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, windows

BASE = Path("/scratch/leahk/eduardo.donestevez")
CASE_DIR = BASE / "simulations" / "bem_90_50_S"
SOLID_CSV = CASE_DIR / "solid" / "results" / "structural_report.csv"
FLUID_CSV = CASE_DIR / "fluid" / "bem_fsi_results" / "bem_report.csv"

OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

DLC6_PARKED_DEFLECTION_M = 8.0
FORCE_RAMP_TIME_S = 2.0
STEADY_START_S = 35.0
FFT_START_S = FORCE_RAMP_TIME_S
MAX_FREQUENCY_HZ = 2.5

# V-02 natural frequencies for the unloaded stationary blade. These are shown
# as reference markers only; the FFT peaks below are observed response peaks.
F1_FLAP_HZ = 0.5537166320
F1_EDGE_HZ = 0.6289939601

COLORS = {
    "max_disp": "#1f3c88",
    "tip_mag": "#2e8b57",
    "dlc": "#8a2d2d",
    "steady": "#a6a6a6",
    "edgewise": "#d98e04",
    "flapwise": "#287271",
    "reference": "#444444",
}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    solid = pd.read_csv(SOLID_CSV)
    fluid = pd.read_csv(FLUID_CSV)
    return solid, fluid


def summarize_response(solid: pd.DataFrame, fluid: pd.DataFrame) -> pd.DataFrame:
    steady_solid = solid[solid["Time [s]"] >= STEADY_START_S]
    steady_fluid = fluid[fluid["Time [s]"] >= STEADY_START_S]
    peak_index = solid["Max Disp [m]"].idxmax()
    rows = [
        {
            "case": "bem_90_50_S",
            "dlc6_parked_reference_m": DLC6_PARKED_DEFLECTION_M,
            "peak_max_disp_m": float(solid.loc[peak_index, "Max Disp [m]"]),
            "peak_time_s": float(solid.loc[peak_index, "Time [s]"]),
            "steady_start_s": STEADY_START_S,
            "steady_samples": int(len(steady_solid)),
            "steady_max_disp_mean_m": float(steady_solid["Max Disp [m]"].mean()),
            "steady_max_disp_std_m": float(steady_solid["Max Disp [m]"].std(ddof=1)),
            "steady_max_disp_min_m": float(steady_solid["Max Disp [m]"].min()),
            "steady_max_disp_max_m": float(steady_solid["Max Disp [m]"].max()),
            "steady_tip_mag_mean_m": float(steady_fluid["Tip Disp Mag [m]"].mean()),
            "steady_tip_mag_std_m": float(steady_fluid["Tip Disp Mag [m]"].std(ddof=1)),
        }
    ]
    return pd.DataFrame.from_records(rows)


def one_sided_fft(
    fluid: pd.DataFrame, column: str, start_time_s: float = FFT_START_S
) -> tuple[np.ndarray, np.ndarray, float, float]:
    signal = fluid[fluid["Time [s]"] >= start_time_s][["Time [s]", column]].dropna().copy()
    times = signal["Time [s]"].to_numpy(dtype=float)
    values = signal[column].to_numpy(dtype=float)
    time_step_s = float(np.median(np.diff(times)))
    centered = values - values.mean()
    window = windows.hann(len(centered), sym=False)
    windowed = centered * (window / window.mean())

    frequencies_hz = np.fft.rfftfreq(len(windowed), d=time_step_s)
    amplitudes_m = np.abs(np.fft.rfft(windowed)) * 2.0 / len(windowed)
    amplitudes_m[0] = 0.0
    resolution_hz = 1.0 / (len(windowed) * time_step_s)
    rms_m = float(np.sqrt(np.mean(centered**2)))
    return frequencies_hz, amplitudes_m, resolution_hz, rms_m


def dominant_frequency(frequencies_hz: np.ndarray, amplitudes_m: np.ndarray) -> tuple[float, float]:
    useful = (frequencies_hz >= 0.05) & (frequencies_hz <= MAX_FREQUENCY_HZ)
    useful_frequencies = frequencies_hz[useful]
    useful_amplitudes = amplitudes_m[useful]
    if len(useful_amplitudes) == 0:
        return float("nan"), float("nan")

    prominence = max(float(useful_amplitudes.max()) * 0.02, 1.0e-12)
    bin_spacing = float(np.median(np.diff(useful_frequencies)))
    min_distance_bins = max(1, int(round(0.08 / bin_spacing)))
    peak_indices, _ = find_peaks(
        useful_amplitudes,
        prominence=prominence,
        distance=min_distance_bins,
    )
    if len(peak_indices) == 0:
        peak_index = int(np.argmax(useful_amplitudes))
    else:
        peak_index = int(peak_indices[np.argmax(useful_amplitudes[peak_indices])])

    return float(useful_frequencies[peak_index]), float(useful_amplitudes[peak_index])


def plot_deformation(solid: pd.DataFrame, fluid: pd.DataFrame) -> None:
    summary = summarize_response(solid, fluid).iloc[0]
    steady_solid = solid[solid["Time [s]"] >= STEADY_START_S]
    steady_mean_m = float(summary["steady_max_disp_mean_m"])
    steady_std_m = float(summary["steady_max_disp_std_m"])
    peak_time_s = float(summary["peak_time_s"])
    peak_disp_m = float(summary["peak_max_disp_m"])

    figure, (full_axis, steady_axis) = plt.subplots(
        2,
        1,
        figsize=(10.4, 7.4),
        sharex=False,
        gridspec_kw={"height_ratios": [2.1, 1.0]},
    )

    full_axis.plot(
        solid["Time [s]"],
        solid["Max Disp [m]"],
        color=COLORS["max_disp"],
        lw=2.0,
        label="AeroElast max. estructural",
    )
    full_axis.plot(
        fluid["Time [s]"],
        fluid["Tip Disp Mag [m]"],
        color=COLORS["tip_mag"],
        lw=1.4,
        ls="--",
        alpha=0.85,
        label="BEM punta, magnitud",
    )
    full_axis.axhline(
        DLC6_PARKED_DEFLECTION_M,
        color=COLORS["dlc"],
        lw=1.5,
        ls=":",
        label="NREL DLC 6.x parqueado (~8 m)",
    )
    full_axis.axvline(
        FORCE_RAMP_TIME_S,
        color=COLORS["reference"],
        lw=1.0,
        ls="--",
        alpha=0.75,
        label="fin rampa de fuerza",
    )
    full_axis.axvline(
        STEADY_START_S,
        color=COLORS["steady"],
        lw=1.0,
        ls="-.",
        alpha=0.85,
        label="inicio ventana final",
    )
    full_axis.scatter(
        [peak_time_s],
        [peak_disp_m],
        s=48,
        color=COLORS["dlc"],
        zorder=5,
    )
    full_axis.annotate(
        f"pico {peak_disp_m:.2f} m\nt = {peak_time_s:.2f} s",
        xy=(peak_time_s, peak_disp_m),
        xytext=(4.8, peak_disp_m - 0.25),
        arrowprops={"arrowstyle": "->", "color": COLORS["dlc"], "lw": 1.0},
        fontsize=9,
    )
    full_axis.set_ylabel("Deflexion [m]")
    full_axis.set_title("Caso parqueado V50: deformacion temporal y escala DLC 6.x")
    full_axis.grid(alpha=0.25)
    full_axis.legend(frameon=False, ncol=2, fontsize=8.5)

    steady_axis.plot(
        steady_solid["Time [s]"],
        steady_solid["Max Disp [m]"],
        color=COLORS["max_disp"],
        lw=1.6,
    )
    steady_axis.axhline(steady_mean_m, color=COLORS["reference"], lw=1.2, label="media final")
    steady_axis.fill_between(
        steady_solid["Time [s]"],
        steady_mean_m - steady_std_m,
        steady_mean_m + steady_std_m,
        color=COLORS["steady"],
        alpha=0.25,
        label="+/- 1 sigma",
    )
    steady_axis.set_xlabel("Tiempo [s]")
    steady_axis.set_ylabel("Deflexion [m]")
    steady_axis.set_title(f"Ventana final: {steady_mean_m:.3f} +/- {steady_std_m:.3f} m")
    steady_axis.grid(alpha=0.25)
    steady_axis.legend(frameon=False, fontsize=8.5)

    figure.tight_layout()
    figure.savefig(OUT_DIR / "fig_5_8_1_parked_v50_deformation.png", dpi=260, bbox_inches="tight")
    figure.savefig(OUT_DIR / "fig_5_8_1_parked_v50_deformation.pdf", bbox_inches="tight")
    plt.close(figure)


def plot_frequency(fluid: pd.DataFrame) -> pd.DataFrame:
    channels = [
        {
            "name": "edgewise_proxy_x",
            "label": "Canal X (edgewise proxy)",
            "column": "Tip Disp X [m]",
            "color": COLORS["edgewise"],
            "reference_frequency_hz": F1_EDGE_HZ,
        },
        {
            "name": "flapwise_proxy_y",
            "label": "Canal Y (flapwise proxy)",
            "column": "Tip Disp Y [m]",
            "color": COLORS["flapwise"],
            "reference_frequency_hz": F1_FLAP_HZ,
        },
    ]

    analysis_window = fluid[fluid["Time [s]"] >= FFT_START_S].copy()
    local_time_s = analysis_window["Time [s]"] - FFT_START_S

    figure, (time_axis, spectrum_axis) = plt.subplots(1, 2, figsize=(12.0, 4.9))
    summary_rows: list[dict[str, float | str]] = []
    resolution_hz = float("nan")

    for channel in channels:
        column = str(channel["column"])
        centered_m = analysis_window[column] - analysis_window[column].mean()
        time_axis.plot(
            local_time_s,
            centered_m,
            lw=1.2,
            color=str(channel["color"]),
            label=str(channel["label"]),
        )

        frequencies_hz, amplitudes_m, resolution_hz, rms_m = one_sided_fft(fluid, column)
        useful = frequencies_hz <= MAX_FREQUENCY_HZ
        spectrum_axis.plot(
            frequencies_hz[useful],
            amplitudes_m[useful],
            lw=1.8,
            color=str(channel["color"]),
            label=str(channel["label"]),
        )
        peak_frequency_hz, peak_amplitude_m = dominant_frequency(frequencies_hz, amplitudes_m)
        reference_frequency_hz = float(channel["reference_frequency_hz"])

        spectrum_axis.scatter(
            [peak_frequency_hz],
            [peak_amplitude_m],
            s=35,
            color=str(channel["color"]),
            zorder=5,
        )
        y_offset = -0.035 if peak_amplitude_m > 0.2 else 0.01
        spectrum_axis.annotate(
            f"{peak_frequency_hz:.3f} Hz",
            xy=(peak_frequency_hz, peak_amplitude_m),
            xytext=(peak_frequency_hz + 0.08, peak_amplitude_m + y_offset),
            fontsize=8,
            color=str(channel["color"]),
        )
        summary_rows.append({
            "channel": str(channel["name"]),
            "column": column,
            "fft_start_s": FFT_START_S,
            "frequency_resolution_hz": resolution_hz,
            "dominant_frequency_hz": peak_frequency_hz,
            "dominant_amplitude_m": peak_amplitude_m,
            "rms_after_ramp_m": rms_m,
            "v02_unloaded_reference_hz": reference_frequency_hz,
            "shift_vs_v02_pct": 100.0 * (peak_frequency_hz / reference_frequency_hz - 1.0),
        })

    time_axis.set_title("Oscilacion residual despues de la rampa")
    time_axis.set_xlabel("Tiempo desde t = 2 s [s]")
    time_axis.set_ylabel("Desplazamiento centrado [m]")
    time_axis.grid(alpha=0.25)
    time_axis.legend(frameon=False, fontsize=8.5)

    spectrum_axis.axvline(
        F1_FLAP_HZ,
        color=COLORS["reference"],
        lw=1.2,
        ls=":",
        alpha=0.8,
        label=f"V-02 f1 flap {F1_FLAP_HZ:.3f} Hz",
    )
    spectrum_axis.axvline(
        F1_EDGE_HZ,
        color=COLORS["reference"],
        lw=1.2,
        ls="--",
        alpha=0.8,
        label=f"V-02 f1 edge {F1_EDGE_HZ:.3f} Hz",
    )
    spectrum_axis.set_title(f"FFT one-sided, Hann, Delta f = {resolution_hz:.4f} Hz")
    spectrum_axis.set_xlabel("Frecuencia [Hz]")
    spectrum_axis.set_ylabel("Amplitud [m]")
    spectrum_axis.set_xlim(0.0, MAX_FREQUENCY_HZ)
    spectrum_axis.grid(alpha=0.25)
    spectrum_axis.legend(frameon=False, fontsize=8.0)

    figure.suptitle("Caso parqueado V50: caracterizacion flapwise/edgewise exportada")
    figure.tight_layout()
    figure.savefig(OUT_DIR / "fig_5_8_2_parked_v50_frequency.png", dpi=260, bbox_inches="tight")
    figure.savefig(OUT_DIR / "fig_5_8_2_parked_v50_frequency.pdf", bbox_inches="tight")
    plt.close(figure)

    return pd.DataFrame.from_records(summary_rows)


def main() -> None:
    solid, fluid = load_data()
    response_summary = summarize_response(solid, fluid)
    frequency_summary = plot_frequency(fluid)
    plot_deformation(solid, fluid)

    response_summary.to_csv(DATA_DIR / "parked_flag_v50_response_summary.csv", index=False)
    frequency_summary.to_csv(DATA_DIR / "parked_flag_v50_frequency_summary.csv", index=False)

    print("Saved figures:")
    print(OUT_DIR / "fig_5_8_1_parked_v50_deformation.png")
    print(OUT_DIR / "fig_5_8_2_parked_v50_frequency.png")
    print("\nResponse summary:")
    print(response_summary.to_string(index=False))
    print("\nFrequency summary:")
    print(frequency_summary.to_string(index=False))


if __name__ == "__main__":
    main()
