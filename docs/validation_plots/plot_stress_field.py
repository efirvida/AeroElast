"""
Layer strain-field visualization on the blade shell — IEA 15 MW, rated.

Reads a single steady-state ``fields.vtu`` snapshot from the AeroElast solver
and plots the recovered layer strains on the top and bottom laminate surfaces.
This avoids comparing orthotropic laminate stress scalars against references
that do not share the same constitutive/material-axis convention.

Outputs
-------
docs/figures/fig_5_4_9_layer_strain_field.png  (300 dpi)
docs/figures/fig_5_4_9_layer_strain_field.pdf
stdout: peak recovered strain values per laminate side and per spanwise band

The figure shows the blade unrolled in a (spanwise z, chordwise x') projection
for the spar cap region, color-coded by principal strain, alongside the
spanwise envelope of recovered layer strains.
"""
import pathlib

import matplotlib.pyplot as plt
import meshio
import numpy as np

# ---------------------------------------------------------------------------
# Paths and time selection
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
# Pick a steady-state snapshot. The exact value of T_SNAPSHOT must match an
# existing fields.vtu time-directory name in the run.
CASE_ROOT = BASE / "frontiersin_results_corotational_100s" / "yaw_0" / "corotational"
T_SNAPSHOT = 70.05     # s — steady-state snapshot inside [40, 100] window (VTU export every 0.05 s)
N_BANDS = 30            # spanwise bands for envelope plot

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Load mesh + fields
# ---------------------------------------------------------------------------
vtu_path = CASE_ROOT / f"{T_SNAPSHOT}" / "fields.vtu"
print(f"Reading {vtu_path}")
m = meshio.read(vtu_path)
print(f"  points={m.points.shape}, fields={len(m.point_data)}")

pts = m.points        # shape (N, 3); IEA 15 MW blade: Z is spanwise (0..117 m)
disp = m.point_data["U"]
top_e1 = m.point_data["TOP_epsilon_1"] * 1e6
mid_e1 = m.point_data["MID_epsilon_1"] * 1e6
bot_e1 = m.point_data["BOT_epsilon_1"] * 1e6
top_gmax = m.point_data["TOP_gamma_max"] * 1e6
mid_gmax = m.point_data["MID_gamma_max"] * 1e6
bot_gmax = m.point_data["BOT_gamma_max"] * 1e6

# Spanwise axis: the Z axis of the (undeformed reference) mesh points.
# The mesh stored in fields.vtu is the reference geometry — displacements are
# in the U field and not added to the coordinates.
r_span = pts[:, 2]
chord_x = pts[:, 0]
r_min, r_max = r_span.min(), r_span.max()
print(f"  spanwise range: [{r_min:.2f}, {r_max:.2f}] m")

# ---------------------------------------------------------------------------
# Spanwise envelope of peak recovered strain (top/mid/bottom)
# ---------------------------------------------------------------------------
r_edges = np.linspace(r_min, r_max, N_BANDS + 1)
r_mid = 0.5 * (r_edges[:-1] + r_edges[1:])
layer_e1 = {"TOP": top_e1, "MID": mid_e1, "BOT": bot_e1}
layer_gmax = {"TOP": top_gmax, "MID": mid_gmax, "BOT": bot_gmax}
e1_envelope = {name: np.zeros(N_BANDS) for name in layer_e1}
gmax_envelope = {name: np.zeros(N_BANDS) for name in layer_gmax}

for i in range(N_BANDS):
    mask = (r_span >= r_edges[i]) & (r_span < r_edges[i + 1])
    if mask.sum() == 0:
        for name in layer_e1:
            e1_envelope[name][i] = np.nan
            gmax_envelope[name][i] = np.nan
        continue
    for name, values in layer_e1.items():
        e1_envelope[name][i] = np.nanmax(np.abs(values[mask]))
    for name, values in layer_gmax.items():
        gmax_envelope[name][i] = np.nanmax(np.abs(values[mask]))

print()
print("Peak recovered layer strain:")
print(f"  TOP |epsilon_1|: {np.max(np.abs(top_e1)):.1f} microstrain")
print(f"  MID |epsilon_1|: {np.max(np.abs(mid_e1)):.1f} microstrain")
print(f"  BOT |epsilon_1|: {np.max(np.abs(bot_e1)):.1f} microstrain")
print()
for name in ("TOP", "MID", "BOT"):
    peak_i = np.nanargmax(e1_envelope[name])
    print(
        f"Spanwise band with maximum {name} |epsilon_1|: "
        f"r ≈ {r_mid[peak_i]:.1f} m ({e1_envelope[name][peak_i]:.1f} microstrain)"
    )

# ---------------------------------------------------------------------------
# Plot: 2 rows × 2 cols
#   (a) TOP scatter on (r, x),  (b) BOT scatter on (r, x)
#   (c) |epsilon_1| envelope, (d) |gamma_max| envelope
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(
    f"Recovered laminate layer strains — IEA 15 MW, yaw=0°, rated, t = {T_SNAPSHOT} s\n"
    "(AeroElast shell-3D, MITC4+ + composite layup)",
    fontsize=11, y=1.01,
)

# (a) TOP scatter: spanwise vs. chordwise projection, colour = principal strain
ax = axes[0, 0]
top_lim = np.percentile(np.abs(top_e1), 99)
sc = ax.scatter(r_span, chord_x, c=top_e1, s=1.0, cmap="coolwarm",
                vmin=-top_lim, vmax=top_lim)
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel(r"Chordwise (deformed) $x$ [m]")
ax.set_title("(a) TOP layer recovered ε₁")
plt.colorbar(sc, ax=ax, label="ε₁ [με]")
ax.grid(alpha=0.3)

# (b) BOT scatter
ax = axes[0, 1]
bot_lim = np.percentile(np.abs(bot_e1), 99)
sc = ax.scatter(r_span, chord_x, c=bot_e1, s=1.0, cmap="coolwarm",
                vmin=-bot_lim, vmax=bot_lim)
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel(r"Chordwise (deformed) $x$ [m]")
ax.set_title("(b) BOT layer recovered ε₁")
plt.colorbar(sc, ax=ax, label="ε₁ [με]")
ax.grid(alpha=0.3)

# (c) Spanwise envelope epsilon_1
ax = axes[1, 0]
for name, color in (("TOP", "#e63946"), ("MID", "#1d3557"), ("BOT", "#457b9d")):
    ax.plot(r_mid, e1_envelope[name], "-", color=color, lw=1.4, label=f"{name} peak |ε₁|")
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel("|ε₁| [με]")
ax.set_title(f"(c) Spanwise peak principal strain ({N_BANDS} bands)")
ax.legend(frameon=False)
ax.grid(alpha=0.3)

# (d) Spanwise envelope gamma_max
ax = axes[1, 1]
for name, color in (("TOP", "#e63946"), ("MID", "#1d3557"), ("BOT", "#457b9d")):
    ax.plot(r_mid, gmax_envelope[name], "-", color=color, lw=1.4, label=f"{name} peak |γmax|")
ax.set_xlabel("Spanwise position r [m]")
ax.set_ylabel("|γmax| [με]")
ax.set_title(f"(d) Spanwise peak shear strain ({N_BANDS} bands)")
ax.legend(frameon=False)
ax.grid(alpha=0.3)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_9_layer_strain_field.png"
pdf = OUT_DIR / "fig_5_4_9_layer_strain_field.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"\nSaved: {png}")
print(f"Saved: {pdf}")
