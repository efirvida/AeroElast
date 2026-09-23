#!/usr/bin/env python3
"""
yaw_sweep_plots.py — Yaw-sweep FSI results vs. published reference data
========================================================================

Generates publication-quality comparison plots between the aeroelast FSI
corotational solver results and the reference data from:

  Ma et al. (2025) — "Aeroelastic Analysis of Offshore Wind Turbine Blades
  under Yaw Misalignment", Frontiers in Energy Research
  doi: 10.3389/fenrg.2025.1571567
  Note: Ma et al. uses LL-FVW (Lifting-Line Free Vortex Wake) + GEBT.
  Their "reference" comparison data comes from Zhou et al. [R11] (ALM-GEBT).

  Zhou et al. (2025) [R11] — "Unsteady aeroelastic performance of the 15 MW
  floating offshore wind turbine under surge condition", *Energy* **336**, 138488
  doi: 10.1016/j.energy.2025.138488
  Note: U=10.59 m/s, Ω=7.55 rpm, rated, fixed-platform condition (no surge).
  Table 4: flapwise tip = 13.86 m, edgewise = -1.22 m, torsion = -3.60°.
  Table 6: mean power = 14.76 MW, mean thrust = 2.20 MN.
  Fig. 21 (FFT, fixed): dominant frequency = 1P = 0.126 Hz (gravity excitation).

  Gaertner et al. (2020) [R1] — "Definition of the IEA 15-Megawatt Offshore
  Reference Wind Turbine", NREL/TP-5000-75698
  doi: 10.2172/1603478
  Note: worst-case out-of-plane tip deflection = 22.8 m (DLC analysis,
  used here as design-envelope reference line on flapwise plots)

Plots generated:
  1.  time_series_flap_yaw_N.png   — flapwise time series per yaw angle
  2.  time_series_edge_yaw_N.png   — edgewise time series per yaw angle
  3.  time_series_power_yaw_N.png  — power time series per yaw angle
  4.  yaw_sweep_cp.png             — CP vs yaw: sim + cos³ + article reference
  5.  yaw_sweep_ct.png             — CT vs yaw: sim + article reference
  6.  yaw_sweep_flap.png           — flapwise tip vs yaw: sim + article refs
  7.  yaw_sweep_power.png          — power vs yaw: sim + article reference
  8.  yaw_sweep_oscillation.png    — peak-to-peak oscillation amplitudes
  9.  tip_fft_yaw_N.png           — FFT of flapwise tip displacement per yaw angle

  10. convergence_iterations.png   — preCICE iterations per case
  python tools/yaw_sweep_plots.py \\
        --results-dir /scratch/.../frontiersin_results \\
        [--output-dir ./plots] \\
        [--t-start 20] [--t-end 70] \\
        [--format png]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------------------------------------------------------------------------
# Published reference data
# (All values from the papers cited in docs/iea15mw_validation_reference.md)
# ---------------------------------------------------------------------------

# --- Ma et al. 2025 [R14] / Zhou et al. 2025 [R11] — yaw sweep ---
# NOTE: Ma et al. (2025) [R14] (doi:10.3389/fenrg.2025.1571567) reports power,
# thrust, tilt and yaw moments vs yaw angle (Fig. 15), and spanwise deflection
# *distributions* (Fig. 16) — NOT numeric tip deflection values vs yaw angle.
# Zhou et al. [R11] covers floating platform motion, not yaw sweep.
# No numeric CP, CT, or flapwise tip values were extractable from those PDFs.
# => All previously listed MA2025_* arrays have been removed to maintain
#    scientific integrity. Qualitative trends (decrease with yaw) are cited
#    in text; no numeric reference lines are plotted from these papers.
MA2025_YAW_DEG = [0, 10, 20, 30, 40]  # yaw angles studied in [R14]

# --- Gaertner et al. 2020 [R1] — DLC design envelope ---
# Worst-case out-of-plane (flapwise) tip deflection from IEC 61400-3 DLC analysis
# Used as design-limit reference line; blade-tower clearance is checked against this
DLC_MAX_FLAP = 22.8  # m  (Gaertner et al. 2020 [R1], NREL/TP-5000-75698)

# --- Bernardi et al. 2025 [R17] — LES + ALM + modal CSD, IEA 15 MW ---
# "Large Eddy Simulation of the IEA 15-MW Wind Turbine Using a Two-Way Coupled
# Fluid-Structure Interaction Model", IEA WIND TASK 47 TURBINIA, preprint wes-2025-120
# Config: U_inf=10 m/s (slightly sub-rated), TSR=9, yaw=0°, tower+nacelle included
# Flapwise tip deflection: ~16 m (time-averaged, TN config)
# Note: LES/CFD-CSD consistently predicts higher deflection than BEM-based methods
BERNARDI2025_YAW_DEG = 0  # only yaw=0° in this study
BERNARDI2025_FLAP = 16.0  # m, LES+ALM+CSD, U_inf=10 m/s

# --- V-04 acceptance band (docs/iea15mw_validation_reference.md §5) ---
# Upper bound: Bernardi [R17] LES+CSD = 16 m (BEM-FSI should be below LES)
# Lower bound: ~9 m (physically unreasonable to be below this at rated, yaw=0°)
# Note: no numeric flapwise reference vs yaw is available from BEM/LL-FVW papers;
#       bounds are set from the verified LES reference and physical reasoning.
V04_FLAP_UPPER = BERNARDI2025_FLAP  # BEM-FSI result must be < LES-CSD upper bound
V04_FLAP_BAND = (9.0, BERNARDI2025_FLAP)  # (lower_physical, upper_lescsd)

# --- Zhou et al. 2025 [R11] — LL-FVW+GEBT, U=10.59 m/s, rated, fixed platform ---
# "Unsteady aeroelastic performance of the 15 MW FOWT under surge condition"
# doi: 10.1016/j.energy.2025.138488
# Table 4 (mean values, flexible blade, fixed condition):
ZHOU2025_FLAP_FIXED = 13.86  # m   flapwise tip deflection, LL-FVW+GEBT
ZHOU2025_FLAP_ALM = 14.10  # m   ALM-GEBT reference used in [R11]
ZHOU2025_EDGE_FIXED = -1.22  # m   edgewise tip deflection (trailing-edge side)
ZHOU2025_TORS_FIXED = -3.60  # deg torsional tip rotation
# Table 6 (mean values, flexible blade, fixed condition):
ZHOU2025_POWER_MW = 14.76  # MW  mean power
ZHOU2025_THRUST_MN = 2.20  # MN  mean rotor thrust
# Fig. 21 (FFT of flapwise tip deflection, fixed condition):
ZHOU2025_FFT_1P_HZ = 0.126  # Hz  dominant frequency = 1P (gravity excitation)

# --- Rated conditions ---
P_RATED_MW = 15.0  # MW
CP_RATED = 0.4618  # [R3] aerodynamic

# ---------------------------------------------------------------------------
# CSV loaders (same logic as yaw_sweep_metrics.py — no shared import needed)
# ---------------------------------------------------------------------------


def load_bem_csv(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({k: float(v) for k, v in row.items()})
    return rows


def filter_window(rows: list[dict], t_start: float, t_end: float) -> list[dict]:
    return [r for r in rows if t_start < r["Time [s]"] <= t_end]


def _mean(xs: list[float]) -> float:
    return sum(xs) / len(xs) if xs else float("nan")


def _std(xs: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    m = _mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def discover_csvs(results_dir: Path) -> list[tuple[int, Path]]:
    cases = []
    for yaw_dir in sorted(results_dir.iterdir()):
        if not yaw_dir.is_dir() or not yaw_dir.name.startswith("yaw_"):
            continue
        try:
            yaw_deg = int(yaw_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue
        csv_path = yaw_dir / "fluid" / "bem_report.csv"
        if csv_path.exists():
            cases.append((yaw_deg, csv_path))
    return sorted(cases)


def parse_iterations_log(path: Path) -> list[tuple[float, int]]:
    """Return list of (sim_time_s, n_iterations), dt=0.01."""
    rows = []
    with open(path) as f:
        lines = f.readlines()
    data_start = 0
    for i, line in enumerate(lines):
        if "TimeWindow" in line:
            data_start = i + 1
            break
    dt = 0.01
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) >= 3:
            try:
                tw = int(parts[0])
                itr = int(parts[2])
                rows.append((tw * dt, itr))
            except (ValueError, IndexError):
                continue
    return rows


# ---------------------------------------------------------------------------
# Plot style helpers
# ---------------------------------------------------------------------------

PALETTE = {
    0: "#1f77b4",  # blue
    10: "#ff7f0e",  # orange
    20: "#2ca02c",  # green
    30: "#d62728",  # red
    40: "#9467bd",  # purple
}

REF_COLOR_PRESENT = "#e377c2"  # pink — Ma2025 "Present"
REF_COLOR_ALM = "#8c564b"  # brown — Ma2025 ALM-GEBT

ARTICLE_MARKER = "s"  # square for article data
SIM_MARKER = "o"  # circle for simulation


def setup_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "lines.linewidth": 1.2,
        "axes.grid": True,
        "grid.alpha": 0.35,
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    })


def add_ss_band(ax, t_start: float, t_end: float):
    """Shade the steady-state analysis window."""
    ax.axvspan(t_start, t_end, color="gold", alpha=0.08, label=f"SS window ({t_start}–{t_end}s)")


def add_mean_line(ax, mean_val: float, color: str, label: str):
    ax.axhline(mean_val, color=color, linestyle="--", linewidth=1.0, alpha=0.7, label=label)


# ---------------------------------------------------------------------------
# 1–3. Time series plots per yaw angle
# ---------------------------------------------------------------------------


def plot_time_series(
    cases: list[tuple[int, list[dict]]],
    col: str,
    ylabel: str,
    title: str,
    output_dir: Path,
    fmt: str,
    t_start: float,
    t_end: float,
    scale: float = 1.0,
    ref_yaw0: float | None = None,
    ref_label: str = "",
    dlc_max: float | None = None,
    dlc_label: str = "",
):
    """One time series per yaw angle in a grid of subplots."""
    n = len(cases)
    ncols = min(n, 2)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(7 * ncols, 3.5 * nrows), sharex=False, sharey=False
    )
    axes = np.array(axes).flatten()

    for i, (yaw_deg, rows) in enumerate(cases):
        ax = axes[i]
        t = [r["Time [s]"] for r in rows]
        val = [r[col] * scale for r in rows]
        color = PALETTE.get(yaw_deg, "#333333")

        ax.plot(t, val, color=color, linewidth=0.7, alpha=0.85, label=f"Sim yaw={yaw_deg}°")
        add_ss_band(ax, t_start, min(t_end, max(t)))

        ss_rows = filter_window(rows, t_start, t_end)
        if ss_rows:
            m = _mean([r[col] * scale for r in ss_rows])
            add_mean_line(ax, m, color, f"Mean ss = {m:.2f}")

        if yaw_deg == 0 and ref_yaw0 is not None:
            ax.axhline(
                ref_yaw0,
                color=REF_COLOR_PRESENT,
                linestyle=":",
                linewidth=1.4,
                label=f"{ref_label} = {ref_yaw0:.2f}",
            )
        if yaw_deg == 0 and dlc_max is not None:
            ax.axhline(
                dlc_max,
                color="firebrick",
                linestyle="--",
                linewidth=1.4,
                label=f"{dlc_label} = {dlc_max} m",
            )

        ax.set_title(f"yaw = {yaw_deg}°")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=8, loc="lower right")

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(title, fontsize=13, y=1.01)
    plt.tight_layout()
    safe_title = col.replace(" ", "_").replace("[", "").replace("]", "").replace("/", "_")
    path = output_dir / f"time_series_{safe_title}.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 4. CP vs yaw
# ---------------------------------------------------------------------------


def plot_cp_vs_yaw(
    sim_yaw: list[int],
    sim_cp: list[float],
    sim_cp_std: list[float],
    output_dir: Path,
    fmt: str,
):
    fig, ax = plt.subplots(figsize=(7, 5))

    # cos^3 curve
    gamma_fine = np.linspace(0, 45, 200)
    cos3 = sim_cp[0] * np.cos(np.radians(gamma_fine)) ** 3
    ax.plot(gamma_fine, cos3, "--", color="gray", linewidth=1.2, label="cos³(γ) fit from yaw=0°")

    # Rated CP line
    ax.axhline(
        CP_RATED,
        color="lightblue",
        linestyle=":",
        linewidth=1.0,
        label=f"CP rated [R3] = {CP_RATED:.4f}",
    )

    # No numeric CP reference values are available from [R14]/[R11] PDFs;
    # cos³ theoretical curve and rated CP line serve as the only references.

    # Simulation
    ax.errorbar(
        sim_yaw,
        sim_cp,
        yerr=sim_cp_std,
        fmt=SIM_MARKER,
        color="#1f77b4",
        markersize=8,
        linewidth=1.5,
        capsize=4,
        label="aeroelast (this work)",
    )

    ax.set_xlabel("Yaw angle γ [°]")
    ax.set_ylabel("CP  [–]")
    ax.set_title("Power coefficient vs. yaw angle\nIEA 15 MW RWT — rated conditions")
    ax.set_xlim(-2, 47)
    ax.legend()
    plt.tight_layout()
    path = output_dir / f"yaw_sweep_cp.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 5. CT vs yaw
# ---------------------------------------------------------------------------


def plot_ct_vs_yaw(
    sim_yaw: list[int],
    sim_ct: list[float],
    sim_ct_std: list[float],
    output_dir: Path,
    fmt: str,
):
    fig, ax = plt.subplots(figsize=(7, 5))

    # cos²(γ) theoretical curve — standard yaw correction for CT
    # (see e.g. Madsen et al. 2020; physically well-established)
    gamma_fine = np.linspace(0, 45, 200)
    cos2 = sim_ct[0] * np.cos(np.radians(gamma_fine)) ** 2
    ax.plot(gamma_fine, cos2, "--", color="gray", linewidth=1.2, label="cos²(γ) fit from yaw=0°")

    # No numeric CT reference values available from [R14]/[R11] PDFs.

    ax.errorbar(
        sim_yaw,
        sim_ct,
        yerr=sim_ct_std,
        fmt=SIM_MARKER,
        color="#1f77b4",
        markersize=8,
        linewidth=1.5,
        capsize=4,
        label="aeroelast (this work)",
    )

    ax.set_xlabel("Yaw angle γ [°]")
    ax.set_ylabel("CT  [–]")
    ax.set_title("Thrust coefficient vs. yaw angle\nIEA 15 MW RWT — rated conditions")
    ax.set_xlim(-2, 47)
    ax.legend()
    plt.tight_layout()
    path = output_dir / f"yaw_sweep_ct.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 6. Flapwise tip vs yaw
# ---------------------------------------------------------------------------


def plot_flap_vs_yaw(
    sim_yaw: list[int],
    sim_flap: list[float],
    sim_flap_std: list[float],
    sim_flap_max: list[float],
    output_dir: Path,
    fmt: str,
):
    fig, ax = plt.subplots(figsize=(7, 5))

    # V-04 acceptance band
    ax.axhspan(*V04_FLAP_BAND, color="lightgreen", alpha=0.15, label="V-04 band (9–15 m)")

    # DLC design-envelope limit from the IEA 15 MW definition report [R1]
    ax.axhline(
        DLC_MAX_FLAP,
        color="firebrick",
        linestyle="--",
        linewidth=1.4,
        label=f"DLC worst-case [R1] = {DLC_MAX_FLAP} m",
    )

    # No numeric flapwise tip vs yaw data available from [R14]/[R11] PDFs
    # (Fig. 16 in [R14] shows spanwise distributions, not tip values vs yaw).
    # Only verified external reference: Bernardi [R17] LES+CSD at yaw=0°.

    # Bernardi 2025 [R17]: single point at yaw=0°, LES+CSD (high-fidelity)
    # U_inf=10 m/s (slightly sub-rated); star marker to distinguish from BEM/LL-FVW references
    ax.plot(
        BERNARDI2025_YAW_DEG,
        BERNARDI2025_FLAP,
        "*",
        color="#17becf",
        markersize=12,
        zorder=5,
        label="Bernardi2025 LES+CSD [R17]\n(U=10 m/s, yaw=0°)",
    )

    # Zhou 2025 [R11]: LL-FVW+GEBT at yaw=0°, rated, fixed platform
    # Table 4: flapwise = 13.86 m (Present LL-FVW), 14.10 m (ALM-GEBT reference)
    # Note: LL-FVW is higher-fidelity than BEM; our BEM result should be below LL-FVW
    ax.errorbar(
        [0],
        [ZHOU2025_FLAP_FIXED],
        yerr=[
            [ZHOU2025_FLAP_FIXED - ZHOU2025_FLAP_FIXED],
            [ZHOU2025_FLAP_ALM - ZHOU2025_FLAP_FIXED],
        ],
        fmt="D",
        color="#8c564b",
        markersize=8,
        capsize=5,
        linewidth=1.2,
        zorder=4,
        label=f"Zhou2025 LL-FVW [R11]\n(Table 4: {ZHOU2025_FLAP_FIXED} m ± ALM {ZHOU2025_FLAP_ALM} m)",
    )

    # Simulation mean ± std
    ax.errorbar(
        sim_yaw,
        sim_flap,
        yerr=sim_flap_std,
        fmt=SIM_MARKER,
        color="#1f77b4",
        markersize=8,
        linewidth=1.5,
        capsize=4,
        label="aeroelast mean ± std",
    )

    # Simulation max (envelope)
    ax.plot(
        sim_yaw,
        sim_flap_max,
        "^",
        color="#1f77b4",
        markersize=6,
        linestyle=":",
        linewidth=0.8,
        alpha=0.6,
        label="aeroelast max",
    )

    ax.set_xlabel("Yaw angle γ [°]")
    ax.set_ylabel("Tip flapwise displacement [m]")
    ax.set_title("Flapwise tip deflection vs. yaw angle\nIEA 15 MW RWT — rated conditions")
    ax.set_xlim(-2, 47)
    ax.set_ylim(8, 25)
    ax.legend()
    plt.tight_layout()
    path = output_dir / f"yaw_sweep_flap.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 7. Power vs yaw
# ---------------------------------------------------------------------------


def plot_power_vs_yaw(
    sim_yaw: list[int],
    sim_power: list[float],
    sim_power_std: list[float],
    output_dir: Path,
    fmt: str,
):
    fig, ax = plt.subplots(figsize=(7, 5))

    # cos^3 curve
    gamma_fine = np.linspace(0, 45, 200)
    cos3 = sim_power[0] * np.cos(np.radians(gamma_fine)) ** 3
    ax.plot(gamma_fine, cos3, "--", color="gray", linewidth=1.2, label="cos³(γ) fit from yaw=0°")

    # Rated power line
    ax.axhline(
        P_RATED_MW,
        color="lightblue",
        linestyle=":",
        linewidth=1.0,
        label=f"P rated = {P_RATED_MW} MW",
    )

    ax.errorbar(
        sim_yaw,
        sim_power,
        yerr=sim_power_std,
        fmt=SIM_MARKER,
        color="#1f77b4",
        markersize=8,
        linewidth=1.5,
        capsize=4,
        label="aeroelast (this work)",
    )

    ax.set_xlabel("Yaw angle γ [°]")
    ax.set_ylabel("Power [MW]")
    ax.set_title("Power vs. yaw angle\nIEA 15 MW RWT — rated conditions")
    ax.set_xlim(-2, 47)
    ax.legend()
    plt.tight_layout()
    path = output_dir / f"yaw_sweep_power.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 8. Oscillation amplitudes
# ---------------------------------------------------------------------------


def plot_oscillation(
    sim_yaw: list[int],
    flap_p2p: list[float],
    edge_p2p: list[float],
    power_p2p: list[float],
    output_dir: Path,
    fmt: str,
):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    data = [
        (axes[0], flap_p2p, "Flapwise p2p [m]", "#1f77b4"),
        (axes[1], edge_p2p, "Edgewise p2p [m]", "#ff7f0e"),
        (axes[2], power_p2p, "Power p2p [MW]", "#2ca02c"),
    ]
    for ax, vals, ylabel, color in data:
        ax.bar(
            [str(y) + "°" for y in sim_yaw],
            vals,
            color=color,
            alpha=0.8,
            edgecolor="k",
            linewidth=0.5,
        )
        ax.set_xlabel("Yaw angle")
        ax.set_ylabel(ylabel)

    axes[0].set_title("Flapwise tip oscillation amplitude\n(peak-to-peak in SS window)")
    axes[1].set_title("Edgewise tip oscillation amplitude\n(peak-to-peak in SS window)")
    axes[2].set_title("Power oscillation amplitude\n(peak-to-peak in SS window)")

    fig.suptitle("3P oscillation envelopes — IEA 15 MW RWT rated", fontsize=12)
    plt.tight_layout()
    path = output_dir / f"yaw_sweep_oscillation.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 9. FFT of flapwise tip displacement — V-07 operational frequency analysis
# ---------------------------------------------------------------------------

# Nominal rotation rate and harmonic markers
_OMEGA_RPM = 7.55  # rpm
_1P_HZ = _OMEGA_RPM / 60  # 0.12583 Hz

# Theoretical nP frequencies
_NP_FREQS = {
    "1P": 1 * _1P_HZ,
    "3P": 3 * _1P_HZ,
    "6P": 6 * _1P_HZ,
    "9P": 9 * _1P_HZ,
    "12P": 12 * _1P_HZ,
}

# First flapwise parked natural frequency (V-02 / [R11] Table 3)
_F1_FLAP_HZ = 0.555  # Hz  [R1]/[R16] reference (parked)


def plot_tip_fft(
    cases: list[tuple[int, list[dict]]],
    output_dir: Path,
    fmt: str,
    t_start: float,
    t_end: float,
):
    """
    FFT of flapwise (Tip Disp Y) and edgewise (Tip Disp X) tip displacement
    in the SS window.  One row per yaw case; left column = flapwise,
    right column = edgewise.

    Reference (V-07): Zhou et al. [R11] Fig. 21 (fixed condition):
      - Flapwise: 1P = 0.126 Hz dominant (aerodynamic excitation).
      - Edgewise: 1P = 0.126 Hz dominant (gravity excitation).
    Also: Bernardi [R17] Fig. 13 (3P, 6P visible in LES+CSD).

    Windowing: Hann window to reduce spectral leakage.
    """
    n = len(cases)
    if n == 0:
        return

    _COLS = [
        ("Tip Disp Y [m]", "Flapwise"),
        ("Tip Disp X [m]", "Edgewise"),
    ]
    _NP_COLORS = {
        "1P": "tab:red",
        "3P": "tab:orange",
        "6P": "tab:green",
        "9P": "tab:purple",
        "12P": "tab:brown",
    }

    fig, axes = plt.subplots(n, 2, figsize=(14, 4 * n), sharex=True, sharey=False)
    # Always 2-D array
    if n == 1:
        axes = axes[np.newaxis, :]

    def _plot_one(ax, signal_raw, freqs, color, yaw_deg, label):
        signal = signal_raw - signal_raw.mean()
        window = np.hanning(len(signal))
        fft_vals = np.fft.rfft(signal * window)
        var = signal.var() if signal.var() > 0 else 1.0
        psd = (np.abs(fft_vals) ** 2) / (var * len(signal) / 2)

        ax.semilogy(freqs, psd, color=color, linewidth=0.8, alpha=0.85)

        for nlabel, fq in _NP_FREQS.items():
            ax.axvline(
                fq,
                color=_NP_COLORS.get(nlabel, "gray"),
                linestyle="--",
                linewidth=0.9,
                alpha=0.75,
                label=f"{nlabel}={fq:.3f} Hz",
            )

        # First flapwise eigenfrequency (parked reference)
        ax.axvline(
            _F1_FLAP_HZ,
            color="black",
            linestyle=":",
            linewidth=1.0,
            alpha=0.6,
            label=f"f₁_flap={_F1_FLAP_HZ:.3f} Hz [R1]",
        )

        # Zhou [R11] 1P reference — solid thin red
        ax.axvline(
            ZHOU2025_FFT_1P_HZ,
            color="tab:red",
            linestyle="-",
            linewidth=1.4,
            alpha=0.35,
        )

        df = freqs[1] - freqs[0] if len(freqs) > 1 else 0
        ax.set_title(
            f"{label} — yaw={yaw_deg}°  (Δf={df:.4f} Hz, N={len(signal)})",
            fontsize=9,
        )
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Norm. PSD [–]")
        ax.set_xlim(0, min(2.0, freqs[-1]))
        ax.legend(fontsize=7, loc="upper right", ncol=2)

    for i, (yaw_deg, rows) in enumerate(cases):
        ss = [r for r in rows if t_start < r["Time [s]"] <= t_end]
        if len(ss) < 64:
            axes[i, 0].set_visible(False)
            axes[i, 1].set_visible(False)
            continue

        dt = ss[1]["Time [s]"] - ss[0]["Time [s]"]
        n_sig = len(ss)
        freqs = np.fft.rfftfreq(n_sig, d=dt)
        color = PALETTE.get(yaw_deg, "#333333")

        for col_idx, (col_key, col_label) in enumerate(_COLS):
            signal_raw = np.array([r[col_key] for r in ss])
            _plot_one(axes[i, col_idx], signal_raw, freqs, color, yaw_deg, col_label)

    fig.suptitle(
        "FFT — Flapwise & Edgewise tip displacement (V-07)\n"
        "IEA 15 MW RWT — rated, SS window\n"
        f"Ref: Zhou [R11] 1P={ZHOU2025_FFT_1P_HZ:.3f} Hz dominant (flap & edge); "
        "Bernardi [R17] Fig. 13",
        fontsize=11,
    )
    plt.tight_layout()
    path = output_dir / f"tip_fft.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 10. BEM spanwise profiles
# ---------------------------------------------------------------------------

# IEA 15 MW blade arc length (hub_radius=0 in BEM config, r goes 0→117 m).
# Do NOT use 120 m (rotor radius incl. hub) — BEM r values top out at ~117 m.
_BLADE_R = 117.0  # m  blade tip radius from YAML reference axis

# --- Rated operating conditions (Ma et al. 2025 [R14] / our YAML) ---
_U_WIND = 10.59  # m/s
_OMEGA = 0.7906  # rad/s  (7.55 rpm)
_PITCH = 0.0  # deg  (from fluid_yaw_0.yaml)

# IEA 15 MW YAML file (tests directory, checked in)
_IEA_YAML = Path(__file__).parent.parent / "tests" / "IEA-15-240-RWT.yaml"


def _load_iea_twist() -> tuple[np.ndarray, np.ndarray] | None:
    """
    Load the design twist from tests/IEA-15-240-RWT.yaml.

    Returns (r_norm, twist_deg) where r_norm is the normalized blade span
    [0, 1] and twist_deg is the design twist in degrees.  Returns None when
    the YAML or the PyYAML package is unavailable.
    """
    try:
        import yaml  # type: ignore
    except ImportError:
        return None
    if not _IEA_YAML.exists():
        return None
    with open(_IEA_YAML) as fh:
        d = yaml.safe_load(fh)
    blade = d["components"]["blade"]["outer_shape_bem"]
    grid = np.asarray(blade["twist"]["grid"], dtype=float)  # 0→1, 50 pts
    twist_rad = np.asarray(blade["twist"]["values"], dtype=float)
    return grid, np.degrees(twist_rad)


def _load_bem_spanwise_mean(
    fluid_dir: Path, t_start: float, t_end: float, stride: int = 10
) -> dict | None:
    """
    Read all fluid/<t>/bem_sectional.csv files in the SS window, subsample by
    *stride*, and return time-averaged spanwise profiles as a dict.

    Keys: r, r_R, alpha, cl, cd, Np, Tp, a, ap  (numpy arrays, shape=(n_stations,))
    """
    # Collect timestep directories in the SS window
    ts_dirs: list[tuple[float, Path]] = []
    for d in fluid_dir.iterdir():
        if not d.is_dir():
            continue
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t_start < t <= t_end:
            ts_dirs.append((t, d))

    ts_dirs.sort(key=lambda x: x[0])
    ts_dirs = ts_dirs[::stride]
    if not ts_dirs:
        return None

    # Column names from first file
    first_path = ts_dirs[0][1] / "bem_sectional.csv"
    if not first_path.exists():
        return None
    with open(first_path) as fh:
        header = fh.readline().strip().split(",")
    col_idx = {h: i for i, h in enumerate(header)}

    # Read all files into a list of 2-D arrays (n_stations × n_cols)
    arrays: list[np.ndarray] = []
    for _, d in ts_dirs:
        p = d / "bem_sectional.csv"
        if p.exists():
            arrays.append(np.loadtxt(p, delimiter=",", skiprows=1))

    if not arrays:
        return None

    stacked = np.stack(arrays, axis=0)  # (n_t, n_stations, n_cols)
    mean = stacked.mean(axis=0)  # (n_stations, n_cols)

    def _col(name: str) -> np.ndarray:
        return mean[:, col_idx[name]]

    r = _col("r[m]")
    return {
        "r": r,
        "r_R": r / _BLADE_R,
        "alpha": _col("alpha[deg]"),
        "cl": _col("cl"),
        "cd": _col("cd"),
        "Np": _col("Np[N/m]"),
        "Tp": _col("Tp[N/m]"),
        "a": _col("a"),
        "ap": _col("ap"),
    }


def plot_bem_spanwise(
    results_dir: Path,
    output_dir: Path,
    fmt: str,
    t_start: float,
    t_end: float,
    stride: int = 10,
):
    """
    Time-averaged BEM spanwise profiles for each yaw case (V-08/V-09).

    Three subplots, one line per yaw case:
      - α(r/R)   — angle of attack [deg]
      - Np(r/R)  — normal force per unit span [kN/m]
      - Tp(r/R)  — tangential force per unit span [kN/m]

    Reference lines added to α subplot:
      - Undisturbed inflow AoA (no induction): arctan(U/(Ω·r)) − pitch − twist_IEA
        Shows the geometric upper bound; BEM values are lower due to axial induction.

    Reference annotation on Np subplot:
      - Zhou et al. [R11] Table 6: thrust 2.20 MN / 3 blades = 733 kN per blade
        compared against ∫Np·dr from our yaw=0° simulation.

    Physical annotation on α subplot:
      - Shaded band 7°–10° for r/R ∈ [0.3, 0.9]: design operating range [R1]
      - Grey fill for r/R < 0.2: inner root (3-D flow dominates BEM assumptions)

    References:
      - [R1]  Gaertner et al. 2020, NREL/TP-5000-75698
      - [R11] Zhou et al. 2025, Energy 336:138488, Fig. 10–11
      - [R14] Ma et al. 2025, Frontiers in Energy Research
    Subsampling by *stride* (default 10) reduces I/O without affecting the mean.
    """
    fig, (ax_alpha, ax_Np, ax_Tp) = plt.subplots(1, 3, figsize=(15, 5))

    # ------------------------------------------------------------------
    # Compute undisturbed inflow AoA reference on a fine r/R grid.
    # alpha_ref(r) = arctan(U_wind / (Ω·r)) − pitch − twist_IEA(r/R_blade)
    # Only valid for r > 0; we skip the innermost hub region (r < 2 m).
    # ------------------------------------------------------------------
    iea_twist = _load_iea_twist()  # (r_norm_0to1, twist_deg) or None
    ref_rR = np.linspace(0.05, 1.0, 200)
    ref_r = ref_rR * _BLADE_R
    ref_phi = np.degrees(np.arctan2(_U_WIND, _OMEGA * ref_r))  # inflow angle
    if iea_twist is not None:
        r_norm_iea, twist_iea_deg = iea_twist
        twist_at_ref = np.interp(ref_rR, r_norm_iea, twist_iea_deg)
        ref_alpha_geo = ref_phi - _PITCH - twist_at_ref
    else:
        ref_alpha_geo = None

    any_data = False
    yaw0_prof: dict | None = None  # keep yaw=0° profile for thrust annotation

    for yaw_dir in sorted(results_dir.iterdir()):
        if not yaw_dir.is_dir() or not yaw_dir.name.startswith("yaw_"):
            continue
        try:
            yaw_deg = int(yaw_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue

        fluid_dir = yaw_dir / "fluid"
        if not fluid_dir.exists():
            continue

        print(f"    Loading BEM sectional: yaw={yaw_deg}° ...", end=" ", flush=True)
        prof = _load_bem_spanwise_mean(fluid_dir, t_start, t_end, stride)
        if prof is None:
            print("no SS data — skipped", file=sys.stderr)
            continue
        print(f"{len(prof['r'])} stations, stride={stride}")
        any_data = True
        if yaw_deg == 0:
            yaw0_prof = prof

        color = PALETTE.get(yaw_deg, "#333333")
        kw = dict(color=color, linewidth=1.4)
        label = f"yaw={yaw_deg}°"
        r_R = prof["r_R"]

        ax_alpha.plot(r_R, prof["alpha"], label=label, **kw)
        ax_Np.plot(r_R, prof["Np"] / 1e3, label=label, **kw)  # N/m → kN/m
        ax_Tp.plot(r_R, prof["Tp"] / 1e3, label=label, **kw)

    if not any_data:
        plt.close()
        return

    # --- α subplot: reference lines and annotations ---
    if ref_alpha_geo is not None:
        ax_alpha.plot(
            ref_rR,
            ref_alpha_geo,
            color="black",
            linestyle="--",
            linewidth=1.2,
            label="Undisturbed inflow\n(no induction) [R14]",
            zorder=3,
        )

    # Design operating band 7–10° for r/R ∈ [0.3, 0.9] (IEA 15 MW, rated [R1])
    ax_alpha.axhspan(
        7.0,
        10.0,
        alpha=0.08,
        color="green",
        label="Design AoA range 7–10° [R1]",
        zorder=1,
    )
    # Root region: BEM assumptions break down due to 3-D flow (r/R < 0.2)
    ax_alpha.axvspan(
        0.0,
        0.2,
        alpha=0.08,
        color="grey",
        label="Root 3-D region",
        zorder=1,
    )
    # Tip region: tip-loss corrections active (r/R > 0.9)
    ax_alpha.axvspan(
        0.9,
        1.0,
        alpha=0.08,
        color="grey",
        zorder=1,
    )

    # --- Np subplot: thrust integral annotation vs [R11] ---
    if yaw0_prof is not None:
        r_arr = yaw0_prof["r"]
        Np_arr = yaw0_prof["Np"]
        T_per_blade_kN = float(np.trapezoid(Np_arr, r_arr)) / 1e3  # N/m·m → kN
        T_rotor_MN = T_per_blade_kN * 3 / 1e3  # 3 blades → MN
        ref_T_kN = ZHOU2025_THRUST_MN * 1e3 / 3.0  # 2.20 MN / 3 blades
        ax_Np.text(
            0.97,
            0.97,
            f"∫Np dr (yaw=0°):\n"
            f"  Sim:   {T_per_blade_kN:.0f} kN/blade ({T_rotor_MN:.2f} MN)\n"
            f"  [R11]: {ref_T_kN:.0f} kN/blade ({ZHOU2025_THRUST_MN:.2f} MN)",
            transform=ax_Np.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7, ec="grey"),
        )

    for ax, ylabel, title in (
        (ax_alpha, "α [deg]", "V-08: Angle of attack"),
        (ax_Np, "Np [kN/m]", "V-09: Normal force per unit span"),
        (ax_Tp, "Tp [kN/m]", "V-09: Tangential force per unit span"),
    ):
        ax.set_xlabel("r/R [–]")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, 1)
        ax.legend(fontsize=7.5)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Time-averaged BEM spanwise profiles — IEA 15 MW RWT, rated\n"
        "pitch = 0°, U = 10.59 m/s, Ω = 7.55 rpm  "
        "|  SS window averaged  |  Ref: [R11] Fig. 10–11, [R14]",
        fontsize=10,
    )
    plt.tight_layout()
    path = output_dir / f"bem_spanwise.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# 11. preCICE convergence
# ---------------------------------------------------------------------------


def plot_convergence(
    results_dir: Path,
    output_dir: Path,
    fmt: str,
    t_start: float,
    t_end: float,
):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    axes = axes.flatten()
    i = 0

    for yaw_dir in sorted(results_dir.iterdir()):
        if not yaw_dir.is_dir() or not yaw_dir.name.startswith("yaw_"):
            continue
        try:
            yaw_deg = int(yaw_dir.name.split("_")[1])
        except (IndexError, ValueError):
            continue

        log = yaw_dir / "precice-Solid-iterations.log"
        if not log.exists() or i >= len(axes):
            continue

        data = parse_iterations_log(log)
        data = [(t, n) for t, n in data if t_start < t <= t_end]
        if not data:
            continue

        times, iters = zip(*data)
        color = PALETTE.get(yaw_deg, "#333333")
        ax = axes[i]

        # rolling average over 50 windows
        iters_arr = np.array(iters, dtype=float)
        window_size = 50
        kernel = np.ones(window_size) / window_size
        iters_smooth = np.convolve(iters_arr, kernel, mode="same")

        ax.plot(times, iters, color=color, alpha=0.25, linewidth=0.5)
        ax.plot(
            times,
            iters_smooth,
            color=color,
            linewidth=1.5,
            label=f"Rolling avg ({window_size} windows)",
        )
        ax.axhline(
            np.mean(iters_arr),
            color="red",
            linestyle="--",
            linewidth=1.0,
            alpha=0.7,
            label=f"Mean = {np.mean(iters_arr):.2f}",
        )
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Coupling iterations")
        ax.set_title(f"yaw = {yaw_deg}°")
        ax.legend(fontsize=8)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        i += 1

    for j in range(i, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("preCICE coupling iterations per time-window\n(SS window only)", fontsize=12)
    plt.tight_layout()
    path = output_dir / f"convergence_iterations.{fmt}"
    plt.savefig(path)
    plt.close()
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Generate yaw-sweep comparison plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--results-dir",
        "-d",
        default="/scratch/leahk/eduardo.donestevez/frontiersin_results",
        help="Root results directory (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        default="/scratch/leahk/eduardo.donestevez/frontiersin_results/plots",
        help="Output directory for plots (default: %(default)s)",
    )
    parser.add_argument(
        "--t-start", type=float, default=20.0, help="SS window start [s] (default: %(default)s)"
    )
    parser.add_argument(
        "--t-end", type=float, default=70.0, help="SS window end [s] (default: %(default)s)"
    )
    parser.add_argument(
        "--min-window",
        type=float,
        default=20.0,
        help="Min seconds in SS window to include case (default: %(default)s)",
    )
    parser.add_argument(
        "--format",
        "-f",
        default="png",
        choices=["png", "pdf", "svg"],
        help="Output image format (default: %(default)s)",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_style()

    print(f"\n  Scanning:   {results_dir}")
    print(f"  Output dir: {output_dir}")
    print(f"  SS window:  {args.t_start} < t <= {args.t_end} s\n")

    cases_raw = discover_csvs(results_dir)
    if not cases_raw:
        print("No bem_report.csv files found.", file=sys.stderr)
        sys.exit(1)

    # Load full time series
    all_cases: list[tuple[int, list[dict]]] = []
    for yaw_deg, csv_path in cases_raw:
        rows = load_bem_csv(csv_path)
        all_cases.append((yaw_deg, rows))
        print(
            f"  [OK] yaw={yaw_deg:>3}°: {len(rows)} rows, "
            f"t={rows[0]['Time [s]']:.1f}–{rows[-1]['Time [s]']:.1f}s"
        )

    # Compute SS metrics
    sim_yaw, sim_flap, sim_flap_std, sim_flap_max = [], [], [], []
    sim_edge, sim_edge_std = [], []
    sim_cp, sim_cp_std, sim_ct, sim_ct_std = [], [], [], []
    sim_power, sim_power_std = [], []
    sim_flap_p2p, sim_edge_p2p, sim_power_p2p = [], [], []

    ss_cases: list[tuple[int, list[dict]]] = []  # only cases with enough SS data

    for yaw_deg, rows in all_cases:
        ss = filter_window(rows, args.t_start, args.t_end)
        if not ss:
            print(f"  [SKIP] yaw={yaw_deg}°: no data in SS window", file=sys.stderr)
            continue
        t_span = ss[-1]["Time [s]"] - ss[0]["Time [s]"]
        if t_span < args.min_window:
            print(f"  [SKIP] yaw={yaw_deg}°: only {t_span:.1f}s in SS window", file=sys.stderr)
            continue

        flap = [r["Tip Disp Y [m]"] for r in ss]
        edge = [r["Tip Disp X [m]"] for r in ss]
        pwr = [r["Power [W]"] / 1e6 for r in ss]
        cp = [r["CP"] for r in ss]
        ct = [r["CT"] for r in ss]

        sim_yaw.append(yaw_deg)
        sim_flap.append(_mean(flap))
        sim_flap_std.append(_std(flap))
        sim_flap_max.append(max(flap))
        sim_edge.append(_mean(edge))
        sim_edge_std.append(_std(edge))
        sim_cp.append(_mean(cp))
        sim_cp_std.append(_std(cp))
        sim_ct.append(_mean(ct))
        sim_ct_std.append(_std(ct))
        sim_power.append(_mean(pwr))
        sim_power_std.append(_std(pwr))
        sim_flap_p2p.append(max(flap) - min(flap))
        sim_edge_p2p.append(max(edge) - min(edge))
        sim_power_p2p.append(max(pwr) - min(pwr))
        ss_cases.append((yaw_deg, rows))

    print(f"\n  Cases in SS metrics: {sim_yaw}\n")

    # --- Time series plots ---
    print("  Generating time series plots ...")
    plot_time_series(
        all_cases,
        "Tip Disp Y [m]",
        "Flapwise tip [m]",
        "Flapwise tip displacement — IEA 15 MW RWT rated, corotational",
        output_dir,
        args.format,
        args.t_start,
        args.t_end,
        ref_yaw0=BERNARDI2025_FLAP,
        ref_label="Bernardi2025 LES+CSD [R17] (U=10 m/s)",
        dlc_max=DLC_MAX_FLAP,
        dlc_label="DLC worst-case [R1]",
    )
    plot_time_series(
        all_cases,
        "Tip Disp X [m]",
        "Edgewise tip [m]",
        "Edgewise tip displacement — IEA 15 MW RWT rated, corotational",
        output_dir,
        args.format,
        args.t_start,
        args.t_end,
    )
    plot_time_series(
        all_cases,
        "Power [W]",
        "Power [MW]",
        "Power — IEA 15 MW RWT rated, corotational",
        output_dir,
        args.format,
        args.t_start,
        args.t_end,
        scale=1e-6,
        ref_yaw0=P_RATED_MW,
        ref_label="P_rated",
    )

    # --- Sweep plots ---
    print("  Generating yaw-sweep summary plots ...")
    if sim_cp:
        plot_cp_vs_yaw(sim_yaw, sim_cp, sim_cp_std, output_dir, args.format)
    if sim_ct:
        plot_ct_vs_yaw(sim_yaw, sim_ct, sim_ct_std, output_dir, args.format)
    if sim_flap:
        plot_flap_vs_yaw(sim_yaw, sim_flap, sim_flap_std, sim_flap_max, output_dir, args.format)
    if sim_power:
        plot_power_vs_yaw(sim_yaw, sim_power, sim_power_std, output_dir, args.format)
    if sim_flap_p2p:
        plot_oscillation(
            sim_yaw, sim_flap_p2p, sim_edge_p2p, sim_power_p2p, output_dir, args.format
        )

    # --- V-07: FFT of flapwise & edgewise tip displacement ---
    print("  Generating FFT plot (V-07) ...")
    plot_tip_fft(ss_cases, output_dir, args.format, args.t_start, args.t_end)

    # --- V-08/V-09: BEM spanwise profiles ---
    print("  Generating BEM spanwise profiles (V-08/V-09) ...")
    plot_bem_spanwise(results_dir, output_dir, args.format, args.t_start, args.t_end)

    # --- Convergence ---
    print("  Generating convergence plot ...")
    plot_convergence(results_dir, output_dir, args.format, args.t_start, args.t_end)

    print(f"\n  All plots saved to: {output_dir}")


if __name__ == "__main__":
    main()
