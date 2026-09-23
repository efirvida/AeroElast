import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

wind = np.array([5.006, 7.159, 9.027, 10.659, 12.259, 15.471])
cp_ref = np.array([0.4164, 0.4616, 0.4616, 0.4618, 0.3036, 0.1511])
ct_ref = np.array([0.7842, 0.7783, 0.7783, 0.7718, 0.3875, 0.1796])

# AeroElast values extracted with module load gcc + pytest/standalone run.
cp_model = np.array([
    0.4234428151,
    0.4736153282,
    0.4736192030,
    0.4738699002,
    0.3078547571,
    0.1537507601,
])
ct_model = np.array([
    0.7947163376,
    0.7932954185,
    0.7932174409,
    0.7865740151,
    0.3918674008,
    0.1808992925,
])

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(wind, cp_ref, "o-", label="Referencia NREL/IEA")
axes[0].plot(wind, cp_model, "s--", label="AeroElast")
axes[0].set_xlabel("Velocidad de viento [m/s]")
axes[0].set_ylabel("Cp [-]")
axes[0].grid(alpha=0.3)
axes[0].legend(frameon=False, fontsize=8)

axes[1].plot(wind, ct_ref, "o-", label="Referencia NREL/IEA")
axes[1].plot(wind, ct_model, "s--", label="AeroElast")
axes[1].set_xlabel("Velocidad de viento [m/s]")
axes[1].set_ylabel("Ct [-]")
axes[1].grid(alpha=0.3)
axes[1].legend(frameon=False, fontsize=8)

fig.suptitle("Curvas BEM de referencia para IEA 15 MW")
fig.tight_layout()

png = OUT_DIR / "fig_4_5_4_6_bem_curves.png"
pdf = OUT_DIR / "fig_4_5_4_6_bem_curves.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
