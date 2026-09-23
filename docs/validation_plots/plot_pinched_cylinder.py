import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N = np.array([4, 8, 16, 32, 64])

regular = {
    "MITC3+": np.array([0.4074, 0.7681, 0.9308, 0.9847, 1.0020]),
    "S4": np.array([0.3882, 0.7544, 0.9328, 0.9902, 1.0080]),
    "MITC4": np.array([0.3788, 0.7469, 0.9286, 0.9871, 1.0050]),
    "MITC4+": np.array([0.3904, 0.7548, 0.9313, 0.9878, 1.0050]),
}

distorted = {
    "MITC3+": np.array([0.4088, 0.7237, 0.8986, 0.9734, 0.9992]),
    "S4": np.array([0.4010, 0.7647, 0.9375, 0.9924, 1.0090]),
    "MITC4": np.array([0.3751, 0.7498, 0.9308, 0.9882, 1.0050]),
    "MITC4+": np.array([0.3793, 0.7535, 0.9321, 0.9886, 1.0050]),
}

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for label, values in regular.items():
    axes[0].plot(N, values, marker="o", label=label)
for label, values in distorted.items():
    axes[1].plot(N, values, marker="s", label=label)

for ax, title in zip(axes, ["Malla regular", "Malla distorsionada"]):
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("N")
    ax.set_title(title)
    ax.grid(alpha=0.3)

axes[0].set_ylabel("w / w_ref")
axes[0].legend(frameon=False, fontsize=8)
axes[1].legend(frameon=False, fontsize=8)
fig.suptitle("Cilindro pinchado (Ko et al. 2017)")
fig.tight_layout()

png = OUT_DIR / "fig_4_2_pinched_cylinder.png"
pdf = OUT_DIR / "fig_4_2_pinched_cylinder.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
