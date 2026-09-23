import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

L = 10.0
cases = [
    (r"$\lambda=\pi/2$", -3.6338, 6.3662),
    (r"$\lambda=\pi$", -10.0, 6.3662),
    (r"$\lambda=3\pi/2$", -12.1220, 2.1221),
    (r"$\lambda=2\pi$", -10.0, 0.0),
]

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot([0.0], [0.0], "ko", label="Raiz")

for label, u_tip, w_tip in cases:
    t = np.linspace(0.0, 1.0, 200)
    x = L * t + u_tip * t**2
    z = w_tip * t**2
    ax.plot(x, z, linewidth=2, label=label)
    ax.plot([x[-1]], [z[-1]], "x", color="black")

ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("u [m]")
ax.set_ylabel("w [m]")
ax.set_title("Cantilever bajo momento terminal (Simo-Vu Quoc)")
ax.grid(alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()

png = OUT_DIR / "fig_4_4_large_rotation.png"
pdf = OUT_DIR / "fig_4_4_large_rotation.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
