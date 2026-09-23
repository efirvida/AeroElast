import pathlib

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

modes = ["1F", "1E", "2F", "2E"]
f_ref = np.array([0.5585, 0.6406, 1.6590, 2.1670])

# AeroElast values extracted with module load gcc + modal solve using tests/NuMAD_utd_iea15mw.xlsx.
f_model = np.array([0.5537166320, 0.6289939601, 1.6946080421, 1.9795958844])

x = np.arange(len(modes))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(x - width / 2, f_ref, width, label="Referencia NREL")
ax.bar(x + width / 2, f_model, width, label="AeroElast")
ax.set_xticks(x)
ax.set_xticklabels(modes)
ax.set_ylabel("Frecuencia [Hz]")
ax.set_title("Frecuencias naturales de pala estacionaria")
ax.grid(axis="y", alpha=0.3)
ax.legend(frameon=False)
fig.tight_layout()

png = OUT_DIR / "fig_4_6_natural_frequencies.png"
pdf = OUT_DIR / "fig_4_6_natural_frequencies.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
