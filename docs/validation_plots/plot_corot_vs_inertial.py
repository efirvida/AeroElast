import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

METRICS = pathlib.Path(__file__).resolve().parents[1] / "validation_data" / "generated" / "campaign_metrics_summary.csv"
OUT_DIR = pathlib.Path(__file__).resolve().parents[1] / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(METRICS)
c = df[(df["yaw_deg"] == 0) & (df["solver"] == "corotational")].iloc[0]
i = df[(df["yaw_deg"] == 0) & (df["solver"] == "inertial")].iloc[0]

labels = ["Flap mean [m]", "Cp mean [-]", "Power mean [MW]"]
corot_vals = np.array([c["flap_mean"], c["cp_mean"], c["power_mean_mw"]])
inert_vals = np.array([i["flap_mean"], i["cp_mean"], i["power_mean_mw"]])

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(x - width / 2, corot_vals, width, label="Corrotacional")
ax.bar(x + width / 2, inert_vals, width, label="Inercial")
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_title("Comparación corrotacional vs inercial (yaw=0°)")
ax.grid(axis="y", alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()

png = OUT_DIR / "fig_4_10_corot_vs_inertial.png"
pdf = OUT_DIR / "fig_4_10_corot_vs_inertial.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
