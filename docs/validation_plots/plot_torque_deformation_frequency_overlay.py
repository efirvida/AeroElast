"""Torque/deformation frequency overlay for the corotational validation report.

The goal is to compare the spectral peaks of rotor-equivalent torque against
tip deformation in the definitive steady window, and to avoid over-reading
every peak as a blade-count harmonic. The script exports:

* docs/figures/fig_5_4_7q_torque_deformation_fft_overlay.{png,pdf}
* docs/figures/fig_5_4_7r_torque_deformation_yaw_marker_map.{png,pdf}
* docs/validation_plots/figures/fig_v04_torque_deformation_fft_overlay.{png,pdf}
* docs/validation_plots/figures/fig_v04_torque_deformation_yaw_marker_map.{png,pdf}
* docs/validation_data/generated/torque_deformation_frequency_peaks.csv
* docs/validation_data/generated/torque_deformation_frequency_markers.csv
* docs/validation_data/generated/torque_deformation_frequency_summary.md
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.signal import find_peaks, windows

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS_DIR = REPO_ROOT / "docs"
OUT_FIG = DOCS_DIR / "figures"
OUT_PLOTS = DOCS_DIR / "validation_plots" / "figures"
OUT_DATA = DOCS_DIR / "validation_data" / "generated"
OUT_FIG.mkdir(parents=True, exist_ok=True)
OUT_PLOTS.mkdir(parents=True, exist_ok=True)
OUT_DATA.mkdir(parents=True, exist_ok=True)

COROT_ROOT = BASE / "frontiersin_results_corotational_100s"
YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
OMEGA_RPM = 7.518
F_1P = OMEGA_RPM / 60.0
N_BLADES = 3
MODE_MATCH_TOL_HZ = 0.025

MODAL_FREQUENCIES = {
    "1F": 0.5537,
    "1E": 0.6290,
    "2F": 1.6946,
    "2E": 1.9796,
    "1T": 4.5347,
}

SIGNALS = {
    "torque": {
        "column": "Torque [N.m]",
        "label": "Torque Q",
        "unit": "MN m",
        "scale": 1.0e6,
        "color": "#c0392b",
    },
    "flapwise": {
        "column": "Tip Disp Y [m]",
        "label": "Tip flapwise Y",
        "unit": "m",
        "scale": 1.0,
        "color": "#16a085",
    },
    "edgewise": {
        "column": "Tip Disp X [m]",
        "label": "Tip edgewise X",
        "unit": "m",
        "scale": 1.0,
        "color": "#8e44ad",
    },
}


@dataclass(frozen=True)
class Spectrum:
    freq: np.ndarray
    amp: np.ndarray
    amp_norm: np.ndarray


def load_window(yaw_deg: int) -> pd.DataFrame:
    csv_path = COROT_ROOT / f"yaw_{yaw_deg}" / "fluid" / "bem_report.csv"
    df = pd.read_csv(csv_path)
    window = df[(df["Time [s]"] >= T_START) & (df["Time [s]"] <= T_END)].copy()
    return window.reset_index(drop=True)


def one_sided_fft(values: np.ndarray, dt: float, nfft_factor: int = 8) -> Spectrum:
    signal = np.asarray(values, dtype=float)
    signal = signal[np.isfinite(signal)]
    n = len(signal)
    window = windows.hann(n)
    centered = signal - float(np.mean(signal))
    tapered = centered * (window / window.mean())
    nfft = int(2 ** np.ceil(np.log2(max(n * nfft_factor, n))))
    spectrum = np.fft.rfft(tapered, n=nfft)
    freq = np.fft.rfftfreq(nfft, d=dt)
    amp = np.abs(spectrum) * 2.0 / n
    amp[0] = 0.0
    max_amp = float(np.max(amp)) if len(amp) else 0.0
    amp_norm = amp / max_amp if max_amp > 0.0 else np.zeros_like(amp)
    return Spectrum(freq=freq, amp=amp, amp_norm=amp_norm)


def nearest_modal(freq: float) -> tuple[str, float, float]:
    label, modal_freq = min(MODAL_FREQUENCIES.items(), key=lambda item: abs(freq - item[1]))
    delta = freq - modal_freq
    return label, modal_freq, delta


def modal_match_label(freq: float) -> tuple[str, float, float, str]:
    label, modal_freq, delta = nearest_modal(freq)
    if abs(delta) <= MODE_MATCH_TOL_HZ:
        return label, modal_freq, delta, f"{label} ({modal_freq:.4f} Hz)"
    return "none", float("nan"), delta, "-"


def classify_blade_filter(order: int) -> str:
    if order == 0:
        return "DC/mean"
    if order % N_BLADES == 0:
        return "survives ideal 3-blade sum"
    return "single-blade/local; cancels in ideal 3-blade sum"


def extract_peaks(
    signal_key: str,
    cfg: dict[str, str | float],
    spectrum: Spectrum,
    yaw_deg: int,
    fmax: float = 2.2,
    max_peaks: int = 6,
) -> list[dict[str, float | int | str]]:
    mask = (spectrum.freq >= 0.04) & (spectrum.freq <= fmax)
    freq = spectrum.freq[mask]
    amp = spectrum.amp[mask]
    amp_norm = spectrum.amp_norm[mask]
    if len(freq) == 0:
        return []
    peak_idx, props = find_peaks(amp_norm, prominence=0.025, distance=8)
    if len(peak_idx) == 0:
        return []

    order_by_height = peak_idx[np.argsort(amp_norm[peak_idx])[::-1]][:max_peaks]
    rows = []
    for idx in order_by_height:
        peak_freq = float(freq[idx])
        p_order = peak_freq / F_1P
        nearest_order = int(round(p_order))
        nearest_p_freq = nearest_order * F_1P
        modal_label, modal_freq, modal_delta, modal_reading = modal_match_label(peak_freq)
        rows.append({
            "yaw_deg": yaw_deg,
            "signal": signal_key,
            "label": str(cfg["label"]),
            "frequency_hz": peak_freq,
            "amplitude": float(amp[idx]),
            "amplitude_unit": str(cfg["unit"]),
            "relative_amplitude": float(amp_norm[idx]),
            "order_p": p_order,
            "nearest_p": nearest_order,
            "nearest_p_hz": nearest_p_freq,
            "delta_nearest_p_hz": peak_freq - nearest_p_freq,
            "nearest_mode": modal_label,
            "nearest_mode_hz": modal_freq,
            "delta_nearest_mode_hz": modal_delta,
            "modal_reading": modal_reading,
            "three_blade_filter": classify_blade_filter(nearest_order),
        })
    return rows


def marker_rows(spectra: dict[str, Spectrum], yaw_deg: int) -> pd.DataFrame:
    markers: list[tuple[str, float, str]] = []
    for order in range(1, 10):
        markers.append((f"{order}P", order * F_1P, classify_blade_filter(order)))
    for label, freq in MODAL_FREQUENCIES.items():
        markers.append((label, freq, "structural modal marker"))

    rows = []
    for signal_key, spectrum in spectra.items():
        cfg = SIGNALS[signal_key]
        for marker, marker_freq, reading in markers:
            idx = int(np.argmin(np.abs(spectrum.freq - marker_freq)))
            rows.append({
                "yaw_deg": yaw_deg,
                "signal": signal_key,
                "label": str(cfg["label"]),
                "marker": marker,
                "marker_frequency_hz": marker_freq,
                "nearest_fft_frequency_hz": float(spectrum.freq[idx]),
                "amplitude": float(spectrum.amp[idx]),
                "amplitude_unit": str(cfg["unit"]),
                "relative_amplitude": float(spectrum.amp_norm[idx]),
                "reading": reading,
            })
    return pd.DataFrame(rows)


def write_outputs(df: pd.DataFrame, spectra: dict[str, Spectrum], dt: float) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.2), sharex=True)
    fig.suptitle(
        "Torque-deformation spectral overlay - corotational yaw=0",
        fontsize=11,
    )

    ax = axes[0]
    for key, cfg in SIGNALS.items():
        spectrum = spectra[key]
        mask = spectrum.freq <= 2.2
        ax.plot(
            spectrum.freq[mask],
            spectrum.amp_norm[mask],
            color=str(cfg["color"]),
            lw=1.2,
            label=str(cfg["label"]),
        )
    for order in range(1, 10):
        freq = order * F_1P
        if freq > 2.2:
            break
        color = "#2c3e50" if order % N_BLADES == 0 else "#95a5a6"
        ax.axvline(
            freq, color=color, lw=0.85, ls=":" if order % N_BLADES == 0 else "--", alpha=0.75
        )
        ax.text(
            freq,
            1.03,
            f"{order}P",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=6.5,
            color=color,
        )
    for label, freq in MODAL_FREQUENCIES.items():
        if freq <= 2.2:
            ax.axvline(freq, color="#f39c12", lw=1.0, ls="-.", alpha=0.85)
            ax.text(
                freq + 0.008,
                0.72,
                label,
                ha="left",
                va="center",
                rotation=90,
                fontsize=6.5,
                color="#b9770e",
            )
    ax.set_ylabel("Normalized FFT amplitude")
    ax.set_ylim(0.0, 1.12)
    ax.set_title("(a) Same-frequency overlay, each signal normalized by its own peak")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, ncol=3, loc="upper right")

    ax = axes[1]
    top = df.sort_values(["signal", "relative_amplitude"], ascending=[True, False])
    for key, cfg in SIGNALS.items():
        rows = top[top["signal"] == key].head(4)
        ax.scatter(
            rows["order_p"],
            rows["relative_amplitude"],
            s=48,
            color=str(cfg["color"]),
            label=str(cfg["label"]),
            zorder=3,
        )
        for _, row in rows.iterrows():
            ax.text(
                float(row["order_p"]) + 0.04,
                float(row["relative_amplitude"]) + 0.015,
                f"{row['frequency_hz']:.3f} Hz\n{int(row['nearest_p'])}P"
                if row["nearest_mode"] == "none"
                else f"{row['frequency_hz']:.3f} Hz\n{row['nearest_mode']}",
                fontsize=6.3,
                color=str(cfg["color"]),
            )
    for order in range(1, 10):
        color = "#2c3e50" if order % N_BLADES == 0 else "#bdc3c7"
        ax.axvline(order, color=color, lw=0.8, ls=":" if order % N_BLADES == 0 else "--", alpha=0.6)
    ax.set_xlabel("Frequency order f / 1P")
    ax.set_ylabel("Relative peak amplitude")
    ax.set_title("(b) Peak locations as P-order; multiples of 3 survive ideal 3-blade summation")
    ax.set_xlim(0.0, 10.0)
    ax.set_ylim(0.0, 1.12)
    ax.grid(alpha=0.25, lw=0.5)

    handles = [
        Line2D([0], [0], color="#2c3e50", ls=":", lw=1.2, label="3P, 6P, 9P: ideal 3-blade global"),
        Line2D(
            [0],
            [0],
            color="#95a5a6",
            ls="--",
            lw=1.2,
            label="Other mP: local/single-blade unless symmetry breaks",
        ),
        Line2D([0], [0], color="#f39c12", ls="-.", lw=1.2, label="Structural modal frequency"),
    ]
    axes[0].legend(
        handles=axes[0].get_legend_handles_labels()[0] + handles,
        labels=axes[0].get_legend_handles_labels()[1] + [handle.get_label() for handle in handles],
        frameon=False,
        fontsize=7,
        ncol=2,
        loc="upper right",
    )

    fig.text(
        0.01,
        0.01,
        f"Window: {T_START:.0f}-{T_END:.0f} s, dt={dt:.3f} s, display uses zero-padding; physical resolution ~{1 / (T_END - T_START):.4f} Hz.",
        fontsize=7,
        color="#555",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))

    for out_dir, basename in [
        (OUT_PLOTS, "fig_v04_torque_deformation_fft_overlay"),
        (OUT_FIG, "fig_5_4_7q_torque_deformation_fft_overlay"),
    ]:
        png = out_dir / f"{basename}.png"
        pdf = out_dir / f"{basename}.pdf"
        fig.savefig(png, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        print(f"Wrote {png}")
        print(f"Wrote {pdf}")
    plt.close(fig)


def write_yaw_marker_map(markers: pd.DataFrame) -> None:
    marker_order = ["1P", "2P", "3P", "4P", "1F", "1E"]
    fig = plt.figure(figsize=(13.2, 4.6))
    grid = fig.add_gridspec(
        1,
        4,
        width_ratios=[1.0, 1.0, 1.0, 0.045],
        left=0.06,
        right=0.94,
        top=0.80,
        bottom=0.24,
        wspace=0.18,
    )
    axes = [fig.add_subplot(grid[0, 0])]
    axes.append(fig.add_subplot(grid[0, 1], sharey=axes[0]))
    axes.append(fig.add_subplot(grid[0, 2], sharey=axes[0]))
    cax = fig.add_subplot(grid[0, 3])
    fig.suptitle(
        "Yaw sweep of torque-deformation spectral markers - corotational solver",
        fontsize=11,
        y=0.94,
    )
    for ax, signal_key in zip(axes, ["torque", "flapwise", "edgewise"]):
        subset = markers[(markers["signal"] == signal_key) & (markers["marker"].isin(marker_order))]
        pivot = (
            subset
            .pivot_table(
                index="yaw_deg",
                columns="marker",
                values="relative_amplitude",
                aggfunc="first",
            )
            .reindex(index=YAWS, columns=marker_order)
            .fillna(0.0)
        )
        im = ax.imshow(
            pivot.to_numpy(dtype=float), vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto"
        )
        ax.set_title(str(SIGNALS[signal_key]["label"]))
        ax.set_xticks(np.arange(len(marker_order)), labels=marker_order, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(YAWS)), labels=[str(yaw) for yaw in YAWS])
        ax.tick_params(axis="y", labelleft=signal_key == "torque")
        ax.set_xlabel("Marker")
        ax.grid(False)
        for row_idx, yaw in enumerate(YAWS):
            for col_idx, marker in enumerate(marker_order):
                value = float(pivot.loc[yaw, marker])
                text_color = "white" if value > 0.55 else "black"
                ax.text(
                    col_idx,
                    row_idx,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=7,
                    color=text_color,
                )
    axes[0].set_ylabel("Yaw [deg]")
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Relative FFT amplitude")
    fig.text(
        0.06,
        0.08,
        "Values are normalized by the maximum FFT amplitude of each signal at each yaw; compare marker hierarchy, not absolute amplitudes.",
        fontsize=7,
        color="#555",
    )
    for out_dir, basename in [
        (OUT_PLOTS, "fig_v04_torque_deformation_yaw_marker_map"),
        (OUT_FIG, "fig_5_4_7r_torque_deformation_yaw_marker_map"),
    ]:
        png = out_dir / f"{basename}.png"
        pdf = out_dir / f"{basename}.pdf"
        fig.savefig(png, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        print(f"Wrote {png}")
        print(f"Wrote {pdf}")
    plt.close(fig)


def write_summary(peaks: pd.DataFrame, markers: pd.DataFrame, dt: float) -> None:
    csv_path = OUT_DATA / "torque_deformation_frequency_peaks.csv"
    peaks.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path}")

    markers_path = OUT_DATA / "torque_deformation_frequency_markers.csv"
    markers.to_csv(markers_path, index=False)
    print(f"Wrote {markers_path}")

    yaw0_peaks = peaks[peaks["yaw_deg"] == 0]
    yaw0_markers = markers[markers["yaw_deg"] == 0]
    dominant = (
        yaw0_peaks
        .sort_values(["signal", "relative_amplitude"], ascending=[True, False])
        .groupby("signal")
        .head(3)
    )
    lines = [
        "# Torque-deformation frequency overlay",
        "",
        f"Source: `frontiersin_results_corotational_100s/yaw_*/fluid/bem_report.csv`, window `{T_START:.0f}-{T_END:.0f} s`.",
        f"Rotation: `1P = {F_1P:.5f} Hz`; for a three-blade ideal rotor, `3P`, `6P`, `9P`, ... survive exact phase summation, while non-multiples of 3 are local/single-blade harmonics unless symmetry is broken.",
        f"Frequency-bin resolution is approximately `{1 / (T_END - T_START):.5f} Hz`; zero-padding is used only to make peak locations easier to read in the figure.",
        "",
        "## Dominant detected peaks at yaw=0",
        "",
        "| Signal | f [Hz] | f/1P | nearest P | modal match | rel. amp. | 3-blade reading |",
        "|---|---:|---:|---:|---|---:|---|",
    ]
    for _, row in dominant.iterrows():
        lines.append(
            f"| {row['label']} | {row['frequency_hz']:.4f} | {row['order_p']:.2f} | "
            f"{int(row['nearest_p'])}P | {row['modal_reading']} | "
            f"{row['relative_amplitude']:.3f} | {row['three_blade_filter']} |"
        )
    lines += [
        "",
        "## Marker amplitudes",
        "",
        "| Marker | Torque rel. amp. | Flapwise rel. amp. | Edgewise rel. amp. | Reading |",
        "|---|---:|---:|---:|---|",
    ]
    marker_table = yaw0_markers[yaw0_markers["marker"].isin(["1P", "2P", "3P", "4P", "1F", "1E"])]
    for marker in ["1P", "2P", "3P", "4P", "1F", "1E"]:
        subset = marker_table[marker_table["marker"] == marker]
        values = {row["signal"]: float(row["relative_amplitude"]) for _, row in subset.iterrows()}
        reading = str(subset["reading"].iloc[0]) if not subset.empty else ""
        lines.append(
            f"| {marker} | {values.get('torque', float('nan')):.3f} | "
            f"{values.get('flapwise', float('nan')):.3f} | "
            f"{values.get('edgewise', float('nan')):.3f} | {reading} |"
        )
    lines += [
        "",
        "## Yaw-sweep marker hierarchy",
        "",
        "| Yaw [deg] | Torque dominant marker | Torque 1P/2P/3P | Flapwise dominant marker | Flapwise 1P/3P/4P | Edgewise dominant marker | Edgewise 1P/2P |",
        "|---:|---|---:|---|---:|---|---:|",
    ]
    marker_subset = markers[markers["marker"].isin(["1P", "2P", "3P", "4P", "1F", "1E"])]
    for yaw in YAWS:
        yaw_rows = marker_subset[marker_subset["yaw_deg"] == yaw]

        def value(signal: str, marker: str) -> float:
            row = yaw_rows[(yaw_rows["signal"] == signal) & (yaw_rows["marker"] == marker)]
            if row.empty:
                return float("nan")
            return float(row["relative_amplitude"].iloc[0])

        def dominant_marker(signal: str) -> str:
            signal_rows = yaw_rows[yaw_rows["signal"] == signal]
            if signal_rows.empty:
                return "-"
            row = signal_rows.sort_values("relative_amplitude", ascending=False).iloc[0]
            return f"{row['marker']} ({row['relative_amplitude']:.2f})"

        lines.append(
            f"| {yaw} | {dominant_marker('torque')} | "
            f"{value('torque', '1P'):.2f}/{value('torque', '2P'):.2f}/{value('torque', '3P'):.2f} | "
            f"{dominant_marker('flapwise')} | "
            f"{value('flapwise', '1P'):.2f}/{value('flapwise', '3P'):.2f}/{value('flapwise', '4P'):.2f} | "
            f"{dominant_marker('edgewise')} | "
            f"{value('edgewise', '1P'):.2f}/{value('edgewise', '2P'):.2f} |"
        )
    lines += [
        "",
        "## Reading",
        "",
        "- At yaw=0, torque and edgewise deformation share the strongest low-order rotational content near 1P, with visible 2P and weaker 3P content.",
        "- Across yaw, the behaviour is not identical: torque remains dominated by low-order rotational content, but flapwise dominance moves between 1P and 4P, with significant 3P content at low yaw.",
        "- The yaw dependence is physical information, not a failure of the method: yaw reorganizes inflow azimuthally, changes phase, and changes how local harmonics survive in the rotor-equivalent signal.",
        "- The strongest yaw=0 flapwise peak falls near 4P, not on the static 1F marker; it should be read as rotational-order content in the modal neighbourhood, not as a closed modal-identification claim.",
        "- Peaks near integer `mP` but not at multiples of 3 should not be interpreted as true rotor-global harmonics without a phase-summed three-blade structural response.",
        "- The clearly modal-matched peaks near 1E are weak in all three signals; they are much smaller than the low-order rotational content.",
    ]
    md_path = OUT_DATA / "torque_deformation_frequency_summary.md"
    md_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {md_path}")


def main() -> None:
    peak_rows: list[dict[str, float | int | str]] = []
    marker_frames: list[pd.DataFrame] = []
    yaw0_spectra: dict[str, Spectrum] | None = None
    yaw0_dt = float("nan")

    for yaw in YAWS:
        df = load_window(yaw)
        time = df["Time [s]"].to_numpy(dtype=float)
        dt = float(np.median(np.diff(time)))
        spectra: dict[str, Spectrum] = {}
        for key, cfg in SIGNALS.items():
            values = df[str(cfg["column"])].to_numpy(dtype=float) / float(cfg["scale"])
            spectrum = one_sided_fft(values, dt)
            spectra[key] = spectrum
            peak_rows.extend(extract_peaks(key, cfg, spectrum, yaw))
        marker_frames.append(marker_rows(spectra, yaw))
        if yaw == 0:
            yaw0_spectra = spectra
            yaw0_dt = dt

    peaks = pd.DataFrame(peak_rows)
    markers = pd.concat(marker_frames, ignore_index=True)
    if yaw0_spectra is None:
        raise RuntimeError("yaw=0 data not found; cannot create overlay figure")
    write_outputs(peaks[peaks["yaw_deg"] == 0], yaw0_spectra, yaw0_dt)
    write_yaw_marker_map(markers)
    write_summary(peaks, markers, yaw0_dt)


if __name__ == "__main__":
    main()
