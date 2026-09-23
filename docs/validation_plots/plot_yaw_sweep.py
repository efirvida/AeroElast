import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CSV = pathlib.Path(__file__).resolve().parents[1] / "validation_data" / "generated" / "campaign_metrics_summary.csv"
DATA_DIR = pathlib.Path(__file__).resolve().parents[1] / "validation_data"
MA_FIG16 = DATA_DIR / "ma_2025_fig16_spanwise_deflections.csv"
OUT_DIR = pathlib.Path(__file__).resolve().parents[1] / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(CSV)
df = df[df["solver"].eq("corotational")].sort_values("yaw_deg")
yaw = df["yaw_deg"].to_numpy()
cp = df["cp_mean"].to_numpy()
flap = df["flap_mean"].to_numpy()
cp0 = cp[0]
cp_cos3 = cp0 * np.cos(np.deg2rad(yaw)) ** 3
ma_fig16 = pd.read_csv(MA_FIG16)
ma_yaw = np.array([10, 20, 30, 40], dtype=float)
ma_tip = ma_fig16.loc[np.isclose(ma_fig16["r_R"], 1.0)].iloc[0]
ma_flap = np.array([ma_tip[f"flap_{int(angle)}deg_m"] for angle in ma_yaw])

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(yaw, cp, "o-", label="AeroElast corrotacional")
axes[0].plot(yaw, cp_cos3, "--", label=r"$C_{P,0}\cos^3(\gamma)$")
axes[0].set_xlabel("Yaw [deg]")
axes[0].set_ylabel("Cp [-]")
axes[0].grid(alpha=0.3)
axes[0].legend(frameon=False)

axes[1].plot(yaw, flap, "o-", label="AeroElast corrotacional")
axes[1].plot(
    ma_yaw,
    ma_flap,
    "D-.",
    color="#16a085",
    label="Ma 2025 Fig. 16 tip r/R=1",
)
axes[1].annotate(
    "Ma Fig. 16: medias de punta r/R=1\n(no incluye yaw=0)",
    xy=(20, ma_flap[1]),
    xytext=(18, ma_flap[1] + 0.35),
    fontsize=8,
    color="#16a085",
    arrowprops={"arrowstyle": "->", "color": "#16a085", "lw": 0.8},
)
axes[1].set_xlabel("Yaw [deg]")
axes[1].set_ylabel("Flapwise tip medio [m]")
axes[1].grid(alpha=0.3)
axes[1].legend(frameon=False)

fig.suptitle("Barrido en yaw: desempeño aerodinámico y respuesta estructural")
fig.tight_layout()

png = OUT_DIR / "fig_4_8_4_9_yaw_sweep.png"
pdf = OUT_DIR / "fig_4_8_4_9_yaw_sweep.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
