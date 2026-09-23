"""
Spanwise distribution of angle of attack and induction factors — IEA 15 MW.

Reproduces the comparison style of Zhou et al. 2025 Figure 10
("Comparison of spanwise distributions of the angle of attack along the
IEA-15 MW blade") and the standard BEM induction factor plots.

Two panels per quantity (AoA, axial induction a, tangential induction a'):
  - All yaw angles in the campaign (0°, 10°, 20°, 30°, 40°)
  - Time-averaged in the steady-state window

Outputs
-------
docs/figures/fig_5_4_8_aoa_induction.png  (300 dpi)
docs/figures/fig_5_4_8_aoa_induction.pdf
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths — update for the new campaign
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s"
YAW_ANGLES = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
DT_PLOT = 0.5
ROOT_CUTOFF_R_R = 0.15
OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "validation_data"
ZHOU_FIG10 = DATA_DIR / "zhou_2025_fig10_aoa.csv"
MA_FIG17 = DATA_DIR / "ma_2025_fig17_aoa_flap_velocity.csv"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def time_averaged_sectional(
    fluid_root: pathlib.Path,
    t_start: float = T_START,
    t_end: float = T_END,
    dt_plot: float = DT_PLOT,
) -> pd.DataFrame:
    """Average bem_sectional.csv across the steady-state window."""
    dirs = [p for p in fluid_root.iterdir() if p.is_dir()]
    dirs.sort(key=lambda p: float(p.name))
    last_t = -np.inf
    accum = None
    n = 0
    for d in dirs:
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t < t_start or t > t_end:
            continue
        if t - last_t < dt_plot:
            continue
        csv = d / "bem_sectional.csv"
        if not csv.exists():
            continue
        df = pd.read_csv(csv)
        if accum is None:
            accum = df.copy()
        else:
            accum.iloc[:, 1:] = accum.iloc[:, 1:] + df.iloc[:, 1:].values
        last_t = t
        n += 1
    if accum is None or n == 0:
        return pd.DataFrame()
    accum.iloc[:, 1:] = accum.iloc[:, 1:] / n
    accum.attrs["n_samples"] = n
    return accum


# ---------------------------------------------------------------------------
# Read all yaw angles
# ---------------------------------------------------------------------------
data = {}
for yaw in YAW_ANGLES:
    fluid_dir = CAMPAIGN_ROOT / f"yaw_{yaw}" / "fluid"
    if not fluid_dir.exists():
        print(f"Skipping yaw={yaw} — no fluid dir")
        continue
    df = time_averaged_sectional(fluid_dir)
    if df.empty:
        print(f"Skipping yaw={yaw} — no data in window")
        continue
    data[yaw] = df
    print(f"yaw={yaw}: {df.attrs['n_samples']} samples averaged")

zhou_aoa = pd.read_csv(ZHOU_FIG10) if ZHOU_FIG10.exists() else pd.DataFrame()
ma_aoa = pd.read_csv(MA_FIG17) if MA_FIG17.exists() else pd.DataFrame()

# ---------------------------------------------------------------------------
# Plot 1×3: AoA(r), a(r), a'(r)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
fig.suptitle(
    "Spanwise distribution of AoA and induction factors — IEA 15 MW, rated\n"
    f"AeroElast BEM, time-averaged on t = [{T_START}, {T_END}] s, "
    f"yaw sweep 0°–40°, r/R ≥ {ROOT_CUTOFF_R_R:.2f}",
    fontsize=11,
    y=1.02,
)

cmap = plt.get_cmap("viridis")
colors = {yaw: cmap(i / max(1, len(YAW_ANGLES) - 1)) for i, yaw in enumerate(YAW_ANGLES)}

# (a) Angle of attack
ax = axes[0]
for yaw, df in data.items():
    r_R = df["r[m]"] / df["r[m]"].max()
    mask = r_R >= ROOT_CUTOFF_R_R
    ax.plot(
        r_R[mask],
        df.loc[mask, "alpha[deg]"],
        "o-",
        color=colors[yaw],
        lw=1.4,
        ms=3,
        label=f"yaw = {yaw}°",
    )
if not zhou_aoa.empty:
    ax.plot(
        zhou_aoa["r_R"],
        zhou_aoa["angle_of_attack_deg"],
        "s--",
        color="#1d3557",
        lw=1.3,
        ms=3,
        label="Zhou 2025 Fig. 10",
    )
if not ma_aoa.empty:
    for yaw in [10, 20, 30, 40]:
        ax.plot(
            ma_aoa["r_R"],
            ma_aoa[f"aoa_{yaw}deg_deg"],
            "--",
            color=colors[yaw],
            lw=1.0,
            alpha=0.65,
            label=f"Ma 2025 Fig. 17 yaw={yaw}°",
        )
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$\alpha$ [deg]")
ax.set_title("(a) Angle of attack")
ax.set_xlim(ROOT_CUTOFF_R_R, 1.0)
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

# (b) Axial induction
ax = axes[1]
for yaw, df in data.items():
    r_R = df["r[m]"] / df["r[m]"].max()
    mask = r_R >= ROOT_CUTOFF_R_R
    ax.plot(
        r_R[mask], df.loc[mask, "a"], "o-", color=colors[yaw], lw=1.4, ms=3, label=f"yaw = {yaw}°"
    )
ax.axhline(1 / 3, color="black", ls="--", lw=0.8, alpha=0.5)
ax.text(0.50, 1 / 3 + 0.01, "Betz limit a=1/3", fontsize=7, color="black", alpha=0.7)
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$a$ [-]")
ax.set_title("(b) Axial induction factor")
ax.set_xlim(ROOT_CUTOFF_R_R, 1.0)
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

# (c) Tangential induction
ax = axes[2]
for yaw, df in data.items():
    r_R = df["r[m]"] / df["r[m]"].max()
    mask = r_R >= ROOT_CUTOFF_R_R
    ax.plot(
        r_R[mask], df.loc[mask, "ap"], "o-", color=colors[yaw], lw=1.4, ms=3, label=f"yaw = {yaw}°"
    )
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$a'$ [-]")
ax.set_title("(c) Tangential induction factor")
ax.set_xlim(ROOT_CUTOFF_R_R, 1.0)
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_8_aoa_induction.png"
pdf = OUT_DIR / "fig_5_4_8_aoa_induction.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
