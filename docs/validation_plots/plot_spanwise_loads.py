"""
Spanwise distribution of BEM aerodynamic loads — IEA 15 MW, rated.

Reproduces the comparison style of Zhou et al. 2025 Figure 11
("Comparison of spanwise distributions of the loads along the IEA-15 MW blade
under rated condition") using the per-time-step ``bem_sectional.csv`` files
produced by the AeroElast BEM fluid participant.

Outputs
-------
docs/figures/fig_5_4_6_spanwise_loads.png  (300 dpi)
docs/figures/fig_5_4_6_spanwise_loads.pdf

The figure shows the time-averaged normal force Np(r) and tangential force Tp(r)
distributions along the blade span, computed from the rated-condition campaign
(yaw=0°, t in [40, 100] s window). When a new campaign is available, only the BEM_DIR
path below needs to be updated.
"""

import pathlib
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths — update these when the rerun completes
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s" / "yaw_0" / "fluid"
T_START = 40.0  # s — definitive steady-state window starts here
T_END = 100.0  # s — end of 100 s campaign
DT_PLOT = 0.5  # s — sample every 0.5 s for time-mean (smooth enough)

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "validation_data"
ZHOU_FIG11 = DATA_DIR / "zhou_2025_fig11_loads.csv"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _list_time_dirs(root: pathlib.Path) -> Iterable[pathlib.Path]:
    """Return time-directory paths sorted by their numeric timestamp name."""
    dirs = [p for p in root.iterdir() if p.is_dir()]
    return sorted(dirs, key=lambda p: float(p.name))


def time_averaged_sectional(
    root: pathlib.Path, t_start: float, t_end: float, dt_plot: float
) -> pd.DataFrame:
    """Average bem_sectional.csv across the steady-state window."""
    dirs = list(_list_time_dirs(root))
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
            r_ref = df["r[m]"].values
        else:
            accum.iloc[:, 1:] = accum.iloc[:, 1:] + df.iloc[:, 1:].values
        last_t = t
        n += 1
    if accum is None or n == 0:
        raise RuntimeError(
            f"No bem_sectional.csv samples found in {root} between t={t_start} and t={t_end}"
        )
    accum.iloc[:, 1:] = accum.iloc[:, 1:] / n
    accum.attrs["n_samples"] = n
    return accum


# ---------------------------------------------------------------------------
# Compute and plot
# ---------------------------------------------------------------------------
mean = time_averaged_sectional(CAMPAIGN_ROOT, T_START, T_END, DT_PLOT)
n_samples = mean.attrs["n_samples"]

r = mean["r[m]"].values
r_R = r / r.max()
Np = mean["Np[N/m]"].values
Tp = mean["Tp[N/m]"].values
Mp = mean["Mp[N.m/m]"].values
zhou_loads = pd.read_csv(ZHOU_FIG11) if ZHOU_FIG11.exists() else pd.DataFrame()

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
fig.suptitle(
    "Spanwise aerodynamic load distributions — IEA 15 MW, rated, yaw=0°\n"
    f"AeroElast BEM, time-averaged on t = [{T_START}, {T_END}] s "
    f"({n_samples} samples)",
    fontsize=11,
    y=1.02,
)

# (a) Normal force Np
ax = axes[0]
ax.plot(r_R, Np / 1e3, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM")
if not zhou_loads.empty:
    ax.plot(
        zhou_loads["r_R"],
        zhou_loads["normal_force_kN_m"],
        "s--",
        color="#1d3557",
        lw=1.3,
        ms=3,
        label="Zhou 2025 Fig. 11",
    )
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$N_p$ [kN/m]")
ax.set_title("(a) Normal force (out-of-plane)")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

# (b) Tangential force Tp
ax = axes[1]
ax.plot(r_R, Tp / 1e3, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM")
if not zhou_loads.empty:
    ax.plot(
        zhou_loads["r_R"],
        zhou_loads["tangential_force_kN_m"],
        "s--",
        color="#1d3557",
        lw=1.3,
        ms=3,
        label="Zhou 2025 Fig. 11",
    )
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$T_p$ [kN/m]")
ax.set_title("(b) Tangential force (in-plane)")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

# (c) Pitching moment Mp
ax = axes[2]
ax.plot(r_R, Mp / 1e3, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM")
ax.set_xlabel(r"Spanwise position $r/R$ [-]")
ax.set_ylabel(r"$M_p$ [kN·m/m]")
ax.set_title("(c) Sectional pitching moment")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_6_spanwise_loads.png"
pdf = OUT_DIR / "fig_5_4_6_spanwise_loads.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
print(f"Time-averaged over {n_samples} samples between t = {T_START} and {T_END} s")
print("Spanwise integrals (single blade):")
print(f"  ∫Np dr = {np.trapezoid(Np, r) / 1e3:.1f} kN")
print(f"  ∫Tp dr = {np.trapezoid(Tp, r) / 1e3:.1f} kN")
print(f"  ∫(Tp · r) dr = {np.trapezoid(Tp * r, r) / 1e6:.2f} MN·m  (root bending edgewise)")
print(f"  ∫(Np · r) dr = {np.trapezoid(Np * r, r) / 1e6:.2f} MN·m  (root bending flapwise)")
