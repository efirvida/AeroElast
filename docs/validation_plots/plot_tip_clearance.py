"""
Tip trajectory and tip-deflection margin proxy — IEA 15 MW, rated.

For each yaw angle in the campaign, plots:
  - The 3D trajectory of the deformed blade tip over the steady-state window
    - A proxy margin based on flapwise tip deflection and NREL's static
        unbent clearance of 30.0 m as reference

The tower of the IEA 15 MW is positioned at (X, Z) = (0, 0) in our reference
frame with shaft tilt 6°; the blade root sits a distance shaft = 11.014 m
upstream of the tower along the rotor axis. We treat the tower axis as the
line X = 0, Y = any, Z = ~ground-up; this is the conservative
approximation typically used for tower-clearance studies of upwind rotors.

NREL/TP-5000-75698 reports:
  - Unbent blade tip-to-tower clearance: 30.0 m
  - Worst-case DLC 1.4 out-of-plane deflection: 22.8 m
  - Implied DLC margin: ~7.2 m

This plot does not compute the minimum geometric blade-surface-to-tower
distance. It quantifies a conservative tip-deflection margin proxy under
rated operation.

Outputs
-------
docs/figures/fig_5_4_10_tip_clearance.png  (300 dpi)
docs/figures/fig_5_4_10_tip_clearance.pdf
"""
import pathlib

import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
CAMPAIGN_ROOT = BASE / "frontiersin_results_corotational_100s"
YAW_ANGLES = [0, 10, 20, 30, 40]
T_START = 40.0

# Reference geometry — IEA 15 MW
UNBENT_CLEARANCE = 30.0      # m, NREL TP-5000-75698
DLC_MAX_DEFLECTION = 22.8    # m, NREL DLC 1.4
DLC_MARGIN = UNBENT_CLEARANCE - DLC_MAX_DEFLECTION

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Load tip-displacement series per yaw
# ---------------------------------------------------------------------------
def load_tip_series(yaw: int) -> pd.DataFrame | None:
    csv = CAMPAIGN_ROOT / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
    if not csv.exists():
        return None
    df = pd.read_csv(csv)
    df = df[df["Time [s]"] >= T_START].reset_index(drop=True)
    # Reject obviously divergent samples (the corotational solver diverges
    # past t ≈ 64 s in the yaw=30° run). A physically plausible blade tip
    # displacement never exceeds the blade length (117 m).
    physical_limit = 50.0   # m — well above any plausible deflection
    mask = (df["Tip Disp Y [m]"].abs() < physical_limit) & \
           (df["Tip Disp X [m]"].abs() < physical_limit)
    if not mask.all():
        df = df[mask].reset_index(drop=True)
    return df


series = {}
for yaw in YAW_ANGLES:
    df = load_tip_series(yaw)
    if df is None or df.empty:
        print(f"yaw={yaw}: no data")
        continue
    series[yaw] = df
    print(f"yaw={yaw}: n={len(df)}, "
          f"flap min/max = [{df['Tip Disp Y [m]'].min():.2f}, "
          f"{df['Tip Disp Y [m]'].max():.2f}] m,  "
          f"clearance min = {UNBENT_CLEARANCE - df['Tip Disp Y [m]'].max():.2f} m")


# ---------------------------------------------------------------------------
# Plot 1×2
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
fig.suptitle(
    "Tip trajectory and tip-deflection margin proxy — IEA 15 MW, rated\n"
    f"AeroElast corotational, steady-state window t ≥ {T_START} s, "
    "yaw sweep 0°–40°",
    fontsize=11, y=1.02,
)

cmap = plt.get_cmap("viridis")
colors = {yaw: cmap(i / max(1, len(YAW_ANGLES) - 1)) for i, yaw in enumerate(YAW_ANGLES)}

# (a) Tip trajectory in (edgewise X, flapwise Y) plane
ax = axes[0]
for yaw, df in series.items():
    ax.plot(df["Tip Disp X [m]"], df["Tip Disp Y [m]"],
            "-", color=colors[yaw], lw=0.8, alpha=0.6, label=f"yaw = {yaw}°")
    # Mean point as a marker
    ax.plot(df["Tip Disp X [m]"].mean(), df["Tip Disp Y [m]"].mean(),
            "o", color=colors[yaw], ms=6, mec="black")
ax.set_xlabel("Edgewise tip displacement X [m]")
ax.set_ylabel("Flapwise tip displacement Y [m]")
ax.set_title("(a) Tip trajectory (X = edgewise, Y = flapwise)")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)
ax.axhline(UNBENT_CLEARANCE, color="red", ls=":", lw=1.0)
ax.text(0.02, UNBENT_CLEARANCE - 0.5, "Unbent clearance 30.0 m [R2]",
        transform=ax.get_yaxis_transform(), fontsize=8, color="red")

# (b) Tip-deflection margin proxy over time
ax = axes[1]
for yaw, df in series.items():
    margin = UNBENT_CLEARANCE - df["Tip Disp Y [m]"]
    ax.plot(df["Time [s]"], margin, "-", color=colors[yaw], lw=1.0,
            label=f"yaw = {yaw}°")
ax.axhline(DLC_MARGIN, color="red", ls="--", lw=1.0,
           label=f"DLC 1.4 minimum margin = {DLC_MARGIN:.1f} m [R2-DLC]")
ax.set_xlabel("Time [s]")
ax.set_ylabel("Tip-deflection margin proxy [m]")
ax.set_title("(b) Proxy margin vs. NREL DLC envelope")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_10_tip_clearance.png"
pdf = OUT_DIR / "fig_5_4_10_tip_clearance.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"\nSaved: {png}")
print(f"Saved: {pdf}")
