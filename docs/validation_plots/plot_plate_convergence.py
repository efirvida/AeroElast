import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N = np.array([4, 8, 16, 32, 64])

regular = {
    "t/L=1/100": np.array([0.9321, 0.9817, 0.9947, 0.9980, 0.9988]),
    "t/L=1/1000": np.array([0.9316, 0.9813, 0.9943, 0.9976, 0.9984]),
    "t/L=1/10000": np.array([0.9230, 0.9806, 0.9942, 0.9976, 0.9984]),
}

distorted = {
    "t/L=1/100": np.array([1.0100, 1.0040, 1.0000, 0.9995, 0.9992]),
    "t/L=1/1000": np.array([1.0090, 1.0030, 1.0000, 0.9991, 0.9988]),
    "t/L=1/10000": np.array([1.0040, 1.0030, 1.0000, 0.9990, 0.9988]),
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
axes[0].legend(frameon=False)
axes[1].legend(frameon=False)
fig.suptitle("Placa cuadrada empotrada (MITC4+, Ko et al. 2017)")
fig.tight_layout()

png = OUT_DIR / "fig_4_1_plate_convergence.png"
pdf = OUT_DIR / "fig_4_1_plate_convergence.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
