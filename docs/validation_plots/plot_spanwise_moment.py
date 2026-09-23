"""
Spanwise distribution of blade root bending moment — IEA 15 MW, rated.

Computes the cumulative spanwise bending moment from the BEM sectional loads:

    M_flap(r)   = ∫_r^{R_tip} (r' - r) · Np(r') dr'      [N·m]
    M_edge(r)   = ∫_r^{R_tip} (r' - r) · Tp(r') dr'      [N·m]

These distributions are nominally comparable with the BeamDyn outputs of the
IEA-15-240-RWT package (M(r) along the blade span) for shell-vs-beam validation
of the structural load path.

Outputs
-------
docs/figures/fig_5_4_7_spanwise_moments.png  (300 dpi)
docs/figures/fig_5_4_7_spanwise_moments.pdf
"""
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Reuse the time-averaging helper from plot_spanwise_loads
import importlib.util
_loads_spec = importlib.util.spec_from_file_location(
    "loads_helper",
    pathlib.Path(__file__).resolve().parent / "plot_spanwise_loads.py",
)
# Avoid executing the plot module — just borrow the helper if needed.

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s" / "yaw_0" / "fluid"
T_START = 40.0
T_END = 100.0
DT_PLOT = 0.5
OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def time_averaged_sectional(root, t_start, t_end, dt_plot):
    dirs = [p for p in root.iterdir() if p.is_dir()]
    dirs.sort(key=lambda p: float(p.name))
    last_t = -np.inf
    accum = None
    n = 0
    for d in dirs:
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t < t_start or t > t_end or t - last_t < dt_plot:
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
    accum.iloc[:, 1:] = accum.iloc[:, 1:] / n
    accum.attrs["n_samples"] = n
    return accum


# ---------------------------------------------------------------------------
# Compute moment distributions
# ---------------------------------------------------------------------------
mean = time_averaged_sectional(CAMPAIGN_ROOT, T_START, T_END, DT_PLOT)
r = mean["r[m]"].values
Np = mean["Np[N/m]"].values
Tp = mean["Tp[N/m]"].values

# Cumulative bending: M(r) = ∫_r^{R} (r' - r) · F(r') dr'
# Compute by integrating from the tip back to root.
R_tip = r[-1]


def cumulative_bending(r, F):
    """Return M(r_i) = ∫_{r_i}^{R} (r' - r_i) · F(r') dr', vectorised."""
    M = np.zeros_like(r)
    for i in range(len(r)):
        M[i] = np.trapezoid((r[i:] - r[i]) * F[i:], r[i:])
    return M


M_flap = cumulative_bending(r, Np)
M_edge = cumulative_bending(r, Tp)

# Shear (cumulative integral from tip)
V_flap = np.array([np.trapezoid(Np[i:], r[i:]) for i in range(len(r))])
V_edge = np.array([np.trapezoid(Tp[i:], r[i:]) for i in range(len(r))])

print(f"Root flapwise bending moment M_flap(r=0) = {M_flap[0] / 1e6:.2f} MN·m")
print(f"Root edgewise bending moment M_edge(r=0) = {M_edge[0] / 1e6:.2f} MN·m")
print(f"Root flapwise shear V_flap(r=0)          = {V_flap[0] / 1e3:.1f} kN")
print(f"Root edgewise shear V_edge(r=0)          = {V_edge[0] / 1e3:.1f} kN")

# ---------------------------------------------------------------------------
# Plot 1×2: M_flap(r), M_edge(r)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
fig.suptitle(
    "Spanwise distribution of bending moment — IEA 15 MW, rated, yaw=0°\n"
    f"AeroElast BEM (time-avg t = [{T_START}, {T_END}] s) "
    f"vs. BeamDyn package output [R1] (placeholder)",
    fontsize=11, y=1.02,
)

ax = axes[0]
ax.plot(r, M_flap / 1e6, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM")
# TODO: when BeamDyn M(r) distribution is digitized from the IEA-15-240-RWT
# package outputs, plot it here as a dashed reference line.
# ax.plot(r_bd, M_flap_bd / 1e6, "s--", color="#457b9d", label="BeamDyn [R1]")
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel(r"$M_{\rm flap}$ [MN·m]")
ax.set_title("(a) Flapwise (out-of-plane) bending moment")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

ax = axes[1]
ax.plot(r, M_edge / 1e6, "o-", color="#e63946", lw=1.4, ms=4, label="AeroElast BEM")
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel(r"$M_{\rm edge}$ [MN·m]")
ax.set_title("(b) Edgewise (in-plane) bending moment")
ax.grid(alpha=0.3)
ax.legend(frameon=False)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_7_spanwise_moments.png"
pdf = OUT_DIR / "fig_5_4_7_spanwise_moments.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
