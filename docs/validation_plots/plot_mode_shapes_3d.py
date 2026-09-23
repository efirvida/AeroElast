"""
3D visualisation of the first four blade mode shapes — IEA 15 MW.

This script is a TEMPLATE. To produce the actual figure, first run the modal
solver with VTU export enabled so the eigenvectors are written as point data:

    pytest tests/test_iea15mw_v02_natural_frequencies.py
    # then add to the test or solver call:
    #     modal_solver.export_modes_vtu("modal_modes.vtu", n_modes=10)

When the VTU file is present, this script reads it and produces a 4-panel
figure showing the displacement field of modes 1-4. Mode shapes are colour-
coded by displacement magnitude on the deformed configuration.

Inputs
------
MODES_VTU  pathlib.Path to the modal VTU file (point data must include
           UMAG_MODE_1, UMAG_MODE_2, ... and U_MODE_1, U_MODE_2, ...)

Outputs
-------
docs/figures/fig_5_3_3_mode_shapes_3d.png  (300 dpi)
docs/figures/fig_5_3_3_mode_shapes_3d.pdf
"""
import pathlib

import matplotlib.pyplot as plt
import meshio
import numpy as np

# ---------------------------------------------------------------------------
# Paths — TODO: update once modal VTU export is in place
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
MODES_VTU = BASE / "fem-shell" / "tests" / "NuMAD_utd" / "modal_modes.vtu"
# Expected naming convention in point_data: U_MODE_i (Nx3), UMAG_MODE_i (N,)

N_MODES_TO_PLOT = 4
EXAGGERATION = 30.0   # display scaling for mode shapes (mass-normalised)
MODE_NAMES = ["1st flapwise (0.554 Hz)",
              "1st edgewise (0.629 Hz)",
              "2nd flapwise (1.695 Hz)",
              "2nd edgewise / 3rd flap (mixed, 1.98 Hz)"]

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Guard: file not yet exported — emit a placeholder and exit cleanly
# ---------------------------------------------------------------------------
if not MODES_VTU.exists():
    print(f"[INFO] Modal VTU not found at {MODES_VTU}")
    print("       To regenerate this figure, run the modal solver with VTU export "
          "enabled and re-run this script. See module docstring.")
    # Write a minimal placeholder figure so the doc has *something* until the
    # real eigenvectors are exported.
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5,
            "Mode-shape figure pending —\n"
            "modal solver VTU export not yet generated.\n"
            "See docs/validation_plots/plot_mode_shapes_3d.py",
            ha="center", va="center", fontsize=11, color="#444")
    ax.set_axis_off()
    fig.savefig(OUT_DIR / "fig_5_3_3_mode_shapes_3d.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "fig_5_3_3_mode_shapes_3d.pdf", bbox_inches="tight")
    print(f"Placeholder saved at {OUT_DIR}")
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# Real plot — runs when the modal VTU is available
# ---------------------------------------------------------------------------
m = meshio.read(MODES_VTU)
pts = m.points

fig = plt.figure(figsize=(13, 10))
fig.suptitle("First four blade mode shapes — IEA 15 MW, AeroElast shell-3D",
             fontsize=12, y=0.96)

for i in range(N_MODES_TO_PLOT):
    ax = fig.add_subplot(2, 2, i + 1, projection="3d")
    U_key = f"U_MODE_{i+1}"
    UMAG_key = f"UMAG_MODE_{i+1}"
    if U_key not in m.point_data:
        ax.text(0.5, 0.5, 0.5, f"{U_key} missing", transform=ax.transAxes)
        continue
    U = m.point_data[U_key]
    UMAG = m.point_data[UMAG_key]
    xyz = pts + EXAGGERATION * U / np.abs(U).max()
    sc = ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=UMAG, s=1.0,
                    cmap="viridis")
    ax.set_title(MODE_NAMES[i], fontsize=10)
    ax.view_init(elev=15, azim=-60)
    ax.set_box_aspect((1, 1, 4))
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(sc, ax=ax, shrink=0.5, label="|U_mode|")

fig.tight_layout()
fig.savefig(OUT_DIR / "fig_5_3_3_mode_shapes_3d.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "fig_5_3_3_mode_shapes_3d.pdf", bbox_inches="tight")
print(f"Saved: {OUT_DIR / 'fig_5_3_3_mode_shapes_3d.png'}")
