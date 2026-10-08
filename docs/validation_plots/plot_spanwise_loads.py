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

**Frame, declared before the curves are compared (issue #15).** The two codes do not
publish their loads in the same frame, and this figure used to overlay them as if they
did:

* `BEMSolver`/`bem_sectional.csv` emit the **rotor-plane** pair - ccblade computes
  `cn = cl cos(phi) + cd sin(phi)` with `phi` the inflow angle from the rotor plane, so
  `atan2(Tp, Np) + atan2(Cd, Cl) - twist == alpha` (measured to 1.4e-14 deg).
* Zhou et al. 2025's Fig. 11 is a **section (chord) frame** pair: their published
  `(Np, Tp)` read with `alpha = atan2(Tp, Np) + atan2(Cd, Cl)` reproduces their own
  Fig. 10 angle of attack to 0.37 deg over the blade core, where the rotor-plane reading
  is off by 6.83 deg. `tools/diagnose_zhou_tp_frame.py` measures both.

The rotor-plane pair is rotated into the section frame (`theta` is the local twist in
radians, taken from the campaign's own `twist[deg]` column) before the Zhou overlay is
drawn, so the two tangential curves are the same quantity:

    Np_sec = Np_rot cos(theta) + Tp_rot sin(theta)
    Tp_sec = Tp_rot cos(theta) - Np_rot sin(theta)

The panel (c) moment is unaffected: `Mp` is about the span axis in both frames. The
remaining difference between the curves is the inflow-angle bias the report already
carries in section 5.11, magnified on `Tp` because it is the small difference of two
large terms. Two further caveats travel with Fig. 11 and are stated in section 5.9: it
is the **flexible** case (its digitized thrust is 2.20 MN, the paper's flexible value,
not the rigid 2.53 MN), and its trapezoid carries only 11.68 MW of the paper's own
14.76 MW, so no integral claim is taken through the digitization.
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

r = np.asarray(mean["r[m]"].values, dtype=float)
r_R = r / r.max()
Np = np.asarray(mean["Np[N/m]"].values, dtype=float)
Tp = np.asarray(mean["Tp[N/m]"].values, dtype=float)
Mp = np.asarray(mean["Mp[N.m/m]"].values, dtype=float)

# The campaign's own twist is the chord angle from the rotor plane, so it is the angle
# that takes the rotor-plane pair (what ccblade emits) into the section frame that
# Zhou's Fig. 11 is published in. See the module docstring for the measurement.
theta = np.deg2rad(np.asarray(mean["twist[deg]"].values, dtype=float))
Np_section = Np * np.cos(theta) + Tp * np.sin(theta)
Tp_section = Tp * np.cos(theta) - Np * np.sin(theta)

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

# (b) Tangential force Tp, in the section frame Zhou publishes (see the docstring)
ax = axes[1]
ax.plot(r_R, Tp / 1e3, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM (rotor plane)")
ax.plot(
    r_R,
    Tp_section / 1e3,
    "^-",
    color="#e76f51",
    lw=1.2,
    ms=3,
    label="AeroElast BEM (section frame, comparable)",
)
if not zhou_loads.empty:
    ax.plot(
        zhou_loads["r_R"],
        zhou_loads["tangential_force_kN_m"],
        "s--",
        color="#1d3557",
        lw=1.3,
        ms=3,
        label="Zhou 2025 Fig. 11 (section frame)",
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
print("  frame note: rotor plane as emitted; the section-frame conversion is the")
print("  comparable one against Zhou's Fig. 11 (see the module docstring, issue #15).")
print(f"  ∫Np dr = {np.trapezoid(Np, r) / 1e3:.1f} kN")
print(f"  ∫Tp dr = {np.trapezoid(Tp, r) / 1e3:.1f} kN  (rotor plane)")
print(f"  ∫Tp dr = {np.trapezoid(Tp_section, r) / 1e3:.1f} kN  (section frame, comparable)")
print(
    f"  ∫(Tp · r) dr = {np.trapezoid(Tp * r, r) / 1e6:.2f} MN·m"
    "  (edgewise root moment, rotor plane)"
)
print(
    f"  ∫(Tp · r) dr = {np.trapezoid(Tp_section * r, r) / 1e6:.2f} MN·m"
    "  (edgewise root moment, section frame, comparable)"
)
print(f"  ∫(Np · r) dr = {np.trapezoid(Np * r, r) / 1e6:.2f} MN·m  (root bending flapwise)")
