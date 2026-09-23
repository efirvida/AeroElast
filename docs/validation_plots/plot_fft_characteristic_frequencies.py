"""
FFT analysis of characteristic frequencies — IEA 15 MW, yaw=0°, rated wind.

Reproduces the style of Zhou et al. 2025 Energy §4.2 (Fig. 17-18):
  - Time series and one-sided FFT spectra for key loads
  - Vertical markers at 1P, 3P, f₁_flap, f₁_edge
  - Side-by-side corrotacional vs inercial comparison

Outputs
-------
docs/validation_plots/figures/fig_fft_characteristic_freq.png  (300 dpi)
docs/validation_plots/figures/fig_fft_characteristic_freq.pdf
stdout: amplitude table at each characteristic frequency
"""
import pathlib
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.signal import windows

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
COROT_CSV = BASE / "frontiersin_results_corotational_100s/yaw_0/corotational/rotor_performance.csv"
INERT_CSV = BASE / "frontiersin_results_inertial/yaw_0/inertial/rotor_performance.csv"
OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical parameters — IEA 15 MW rated condition
# ---------------------------------------------------------------------------
OMEGA_RAD = 0.7872          # rad/s  (7.518 RPM)
F_1P   = OMEGA_RAD / (2 * np.pi)   # = 0.12527 Hz
F_3P   = 3 * F_1P                  # = 0.37581 Hz
F_FLAP = 0.5537                     # Hz  (1st flapwise — V-02 AeroElast shell)
F_EDGE = 0.6290                     # Hz  (1st edgewise — V-02 AeroElast shell)
T_START = 40.0                      # s  — skip initial transient (definitive window)
T_END = 100.0                       # s  — end of 100 s campaign

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_trim(csv_path, t_start=T_START):
    df = pd.read_csv(csv_path)
    if "Aero Power Rotor Equivalent [W]" not in df and "Aero Power [W]" in df:
        df["Aero Power Rotor Equivalent [W]"] = 3.0 * df["Aero Power [W]"]
    return df[(df["Time [s]"] >= t_start) & (df["Time [s]"] <= T_END)].reset_index(drop=True)


def onesided_fft(signal, dt):
    """Return (freqs [Hz], amplitude [same units as signal]) using Hann window."""
    n = len(signal)
    win = windows.hann(n)
    # Normalize window so amplitude is preserved
    win_norm = win / win.mean()
    x = (signal.values - signal.mean()) * win_norm
    sp = np.fft.rfft(x)
    freq = np.fft.rfftfreq(n, d=dt)
    amp = np.abs(sp) * 2.0 / n   # two-sided → one-sided; /n for amplitude
    amp[0] = 0.0                  # DC component forced to zero (we subtracted mean)
    return freq, amp


def peak_near(freq, amp, f_target, bw=0.04):
    """Return amplitude at the FFT bin closest to f_target (within bandwidth bw)."""
    mask = np.abs(freq - f_target) <= bw
    if not mask.any():
        return 0.0
    return amp[mask].max()


def fmt_amp(val, scale):
    """Format amplitude with given scale and unit."""
    return f"{val / scale:.3f}"


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
corot = load_trim(COROT_CSV)
inert = load_trim(INERT_CSV)

dt = 0.01  # s

print(f"\nCorrotacional: {len(corot)} samples, "
      f"t=[{corot['Time [s]'].min():.1f}, {corot['Time [s]'].max():.1f}] s, "
      f"Δf={1/(len(corot)*dt):.4f} Hz")
print(f"Inercial:      {len(inert)} samples, "
      f"t=[{inert['Time [s]'].min():.1f}, {inert['Time [s]'].max():.1f}] s, "
      f"Δf={1/(len(inert)*dt):.4f} Hz\n")

# ---------------------------------------------------------------------------
# Signals to analyse
# (label, column, scale, unit)
# ---------------------------------------------------------------------------
SIGNALS = [
    ("Aero Thrust",       "Aero Thrust [N]",      1e3,  "kN"),
    ("Aero Power rotor-eq.", "Aero Power Rotor Equivalent [W]", 1e6, "MW"),
    ("Tilt Moment (T_x)", "Aero Torque X [Nm]",    1e6,  "MNm"),
    ("Yaw Moment (T_z)",  "Aero Torque Z [Nm]",    1e6,  "MNm"),
    ("Tip Displacement",  "Max Displacement [m]",  1.0,  "m"),
]

# Char. frequencies (name, Hz, marker colour)
CHAR_FREQS = [
    ("1P",        F_1P,   "#e63946"),
    ("3P",        F_3P,   "#457b9d"),
    ("f_1,flap",  F_FLAP, "#2a9d8f"),
    ("f_1,edge",  F_EDGE, "#e9c46a"),
]

# ---------------------------------------------------------------------------
# Figure — 5 rows × 2 cols (time-domain | FFT)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(len(SIGNALS), 2, figsize=(14, 3.2 * len(SIGNALS)))
fig.suptitle(
    "Characteristic-frequency FFT analysis — IEA 15 MW, yaw=0°, rated condition\n"
    "Corotational (—) vs Inertial (- -)",
    fontsize=12, y=1.01
)

print("=" * 80)
print(f"{'Signal':<22} | {'Freq':>8} | {'Amp Corrot':>12} | {'Amp Inert':>12} | {'Unit':>4}")
print("-" * 80)

for row, (label, col, scale, unit) in enumerate(SIGNALS):
    ax_t = axes[row, 0]  # time domain
    ax_f = axes[row, 1]  # frequency domain

    # Time series
    tc, yc = corot["Time [s]"], corot[col] / scale
    ti, yi = inert["Time [s]"], inert[col] / scale
    ax_t.plot(tc, yc, lw=0.9, color="#e63946", label="Corotational")
    ax_t.plot(ti, yi, lw=0.9, color="#457b9d", ls="--", alpha=0.8, label="Inertial")
    ax_t.set_ylabel(f"{label}\n[{unit}]", fontsize=8)
    ax_t.set_xlim(T_START, max(tc.max(), ti.max()))
    ax_t.grid(alpha=0.3)
    if row == 0:
        ax_t.legend(fontsize=7, frameon=False)
        ax_t.set_title("Time series", fontsize=8)
    if row == len(SIGNALS) - 1:
        ax_t.set_xlabel("Time [s]", fontsize=8)

    # FFT
    fc, ac = onesided_fft(corot[col] / scale, dt)
    fi, ai = onesided_fft(inert[col] / scale, dt)

    ax_f.semilogy(fc, np.clip(ac, 1e-6, None), lw=0.9, color="#e63946")
    ax_f.semilogy(fi, np.clip(ai, 1e-6, None), lw=0.9, color="#457b9d", ls="--", alpha=0.8)
    ax_f.set_xlim(0, 1.5)
    ax_f.set_ylabel(f"Amplitude [{unit}]", fontsize=8)
    ax_f.grid(alpha=0.3, which="both")

    # Vertical lines for char. frequencies
    for fname, fval, fcolor in CHAR_FREQS:
        ax_f.axvline(fval, color=fcolor, lw=1.2, ls=":", alpha=0.9)
        ax_f.text(fval, ax_f.get_ylim()[1] if ax_f.get_ylim()[1] > 0 else 1.0,
                  fname, color=fcolor, fontsize=6, ha="center", va="bottom",
                  rotation=90)

    if row == len(SIGNALS) - 1:
        ax_f.set_xlabel("Frequency [Hz]", fontsize=8)
    if row == 0:
        ax_f.set_title("FFT spectrum (Hann window, log-y)", fontsize=8)

    # Print amplitude table
    for fname, fval, _ in CHAR_FREQS:
        amp_c = peak_near(fc, ac, fval)
        amp_i = peak_near(fi, ai, fval)
        print(f"{label:<22} | {fname:>3} {fval:.4f} Hz | {amp_c:>12.4f} | {amp_i:>12.4f} | {unit}")

print("=" * 80)

ax_t0 = axes[0, 0]
ax_f0 = axes[0, 1]
# Add characteristic-frequency legend to the first FFT subplot
handles = [
    Line2D([0], [0], color=c, ls=":", lw=1.5, label=f"{n} ({v:.4f} Hz)")
    for n, v, c in CHAR_FREQS
]
handles += [
    Line2D([0], [0], color="#e63946", lw=1.2, label="Corotational"),
    Line2D([0], [0], color="#457b9d", lw=1.2, ls="--", label="Inertial"),
]
ax_f0.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper right")

fig.tight_layout()
png = OUT_DIR / "fig_fft_characteristic_freq.png"
pdf = OUT_DIR / "fig_fft_characteristic_freq.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"\nSaved: {png}")
print(f"Saved: {pdf}")

# Also save a copy under the name referenced from the validation document
report_png = pathlib.Path(__file__).resolve().parent.parent / "figures" / "fig_5_5_1_fft_cargas.png"
report_pdf = pathlib.Path(__file__).resolve().parent.parent / "figures" / "fig_5_5_1_fft_cargas.pdf"
report_png.parent.mkdir(parents=True, exist_ok=True)
shutil.copy(png, report_png)
shutil.copy(pdf, report_pdf)
print(f"Saved (mirror): {report_png}")
print(f"Saved (mirror): {report_pdf}")
