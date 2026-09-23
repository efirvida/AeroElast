"""
3D rotor load map — IEA 15 MW at rated, all three blades.

Reads a single ``fields.vtu`` snapshot (one blade) and synthesises the
three-bladed rotor by rigidly rotating the mesh ±120°. The aerodynamic
force vector F_AERO is colour-coded by magnitude and projected onto the
deformed configuration.

This view shows immediately:
  - How the load magnitude distributes along the span
  - The location of the high-load region (outboard, where Np peaks)
  - For yaw cases, the azimuthal asymmetry across the three blades

Outputs
-------
docs/figures/fig_5_4_12_rotor_load_map.png  (300 dpi)
docs/figures/fig_5_4_12_rotor_load_map.pdf
"""
import pathlib

import matplotlib.pyplot as plt
import meshio
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
CASE_ROOT = BASE / "frontiersin_results_corotational_100s" / "yaw_0" / "corotational"
T_SNAPSHOT = 70.05
DEFORM_SCALE = 1.0       # 1.0 = true deformation; increase to exaggerate

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Load one blade and rotate to form a three-bladed rotor
# ---------------------------------------------------------------------------
vtu = CASE_ROOT / f"{T_SNAPSHOT}" / "fields.vtu"
print(f"Reading {vtu}")
m = meshio.read(vtu)
pts = m.points       # (N, 3) — Z is spanwise, X chordwise, Y vertical (rotor plane)
U = m.point_data["U"]
F = m.point_data["F_AERO"]
F_mag = np.linalg.norm(F, axis=1)
print(f"  Force magnitude range: [{F_mag.min():.2f}, {F_mag.max():.2f}] N")
print(f"  Spanwise extent: [{pts[:,2].min():.1f}, {pts[:,2].max():.1f}] m")

# Deformed coordinates (one blade, in the rotating frame).
xyz_def = pts + DEFORM_SCALE * U

# Build the three-bladed rotor by rotating the deformed shape ±120° about Z.
# In the local frame Z is the blade span and the rotor axis is along X
# (perpendicular to blade). To put the blade radially in the global rotor
# plane we map (local x → tangential, local y → upstream, local z → radial).
# A standard rotor view places X as upstream and the (Y, Z) plane as the
# rotor disk.  For visualisation we use the simple convention:
#   blade frame:  Xb = chord, Yb = flap, Zb = span
#   rotor frame:  X = upstream (out-of-disk), Yr,Zr = disk
# i.e. rotate the blade so its span aligns with the disk radial direction.
def to_rotor_frame(pts_blade: np.ndarray, azimuth_deg: float) -> np.ndarray:
    """Rotate blade-frame points to rotor-frame at the given azimuth."""
    # First reorient: blade-frame (Xb, Yb, Zb) → (Yr, Zr, X) so that span Zb
    # becomes radial Zr.
    R0 = np.array([
        [0.0, 1.0, 0.0],   # X (upstream) ← Yb (flap)
        [1.0, 0.0, 0.0],   # Yr (disk-1)   ← Xb (chord)
        [0.0, 0.0, 1.0],   # Zr (radial)   ← Zb (span)
    ])
    out = pts_blade @ R0.T
    # Then spin around X (rotor axis) by azimuth
    c, s = np.cos(np.deg2rad(azimuth_deg)), np.sin(np.deg2rad(azimuth_deg))
    Raz = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return out @ Raz.T


def rotate_force(F_blade: np.ndarray, azimuth_deg: float) -> np.ndarray:
    R0 = np.array([
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    out = F_blade @ R0.T
    c, s = np.cos(np.deg2rad(azimuth_deg)), np.sin(np.deg2rad(azimuth_deg))
    Raz = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    return out @ Raz.T


# ---------------------------------------------------------------------------
# Plot the three blades
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(111, projection="3d")
fig.suptitle(
    "Three-bladed rotor with aerodynamic load magnitude\n"
    f"IEA 15 MW, yaw=0°, rated, t = {T_SNAPSHOT} s — AeroElast + BEM",
    fontsize=11, y=0.97,
)

vmax = np.percentile(F_mag, 99)
azimuths = [0.0, 120.0, 240.0]
cmap = plt.get_cmap("viridis")

for az in azimuths:
    xyz_r = to_rotor_frame(xyz_def, az)
    sc = ax.scatter(xyz_r[:, 0], xyz_r[:, 1], xyz_r[:, 2],
                    c=F_mag, s=2.0, cmap=cmap, vmin=0, vmax=vmax,
                    alpha=0.8)

cb = fig.colorbar(sc, ax=ax, label=r"$|F_{\rm aero}|$ [N/node]",
                  shrink=0.5, pad=0.1)
# Add a faint tower line at X = -shaft_distance and any Y
ax.plot([-12, -12], [0, 0], [-130, 5], color="grey", lw=1.5, alpha=0.5,
        label="Tower axis (illustrative)")
# Rotor axis
ax.plot([-12, 5], [0, 0], [0, 0], color="black", ls="--", lw=1.0, alpha=0.5,
        label="Rotor axis")

ax.set_xlabel("X (upstream) [m]")
ax.set_ylabel("Y (in-plane) [m]")
ax.set_zlabel("Z (vertical) [m]")
ax.view_init(elev=20, azim=-65)
ax.legend(loc="upper left", fontsize=8, frameon=False)

# Match physical scale
all_xyz = np.vstack([to_rotor_frame(xyz_def, az) for az in azimuths])
half = max(np.abs(all_xyz).max(), 130)
ax.set_xlim(-half * 0.3, half * 0.3)
ax.set_ylim(-half, half)
ax.set_zlim(-half, half)
ax.set_box_aspect((0.6, 2.0, 2.0))

png = OUT_DIR / "fig_5_4_12_rotor_load_map.png"
pdf = OUT_DIR / "fig_5_4_12_rotor_load_map.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
