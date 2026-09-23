import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N = np.array([4, 8, 16, 32, 64])

scordelis = {
    "MITC3+": np.array([0.7982, 0.9223, 0.9757, 0.9909, 0.9948]),
    "MITC4+": np.array([1.0400, 0.9973, 0.9942, 0.9948, 0.9957]),
}
hypar = {
    "t/L=1/1000": np.array([0.5685, 0.8359, 0.9483, 0.9816, 0.9936]),
    "t/L=1/10000": np.array([0.2176, 0.6945, 0.9607, 0.9948, 0.9990]),
}

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for label, values in scordelis.items():
    axes[0].plot(N, values, marker="o", label=label)
for label, values in hypar.items():
    axes[1].plot(N, values, marker="s", label=label)

axes[0].set_title("Scordelis-Lo (malla distorsionada)")
axes[1].set_title("Paraboloide hiperbólico (distorsionada)")
for ax in axes:
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_xscale("log", base=2)
    ax.set_xlabel("N")
    ax.grid(alpha=0.3)

axes[0].set_ylabel("w / w_ref")
axes[0].legend(frameon=False)
axes[1].legend(frameon=False)
fig.suptitle("Convergencia en problemas de cascarón complejos")
fig.tight_layout()

png = OUT_DIR / "fig_4_3_scordelis_hypar.png"
pdf = OUT_DIR / "fig_4_3_scordelis_hypar.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
