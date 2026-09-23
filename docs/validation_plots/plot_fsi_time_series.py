import pathlib

import matplotlib.pyplot as plt
import pandas as pd

BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
COROT_CSV = BASE / "frontiersin_results_corotational_100s/yaw_0/corotational/rotor_performance.csv"
INERTIAL_CSV = BASE / "frontiersin_results_inertial/yaw_0/inertial/rotor_performance.csv"
OUT_DIR = pathlib.Path(__file__).resolve().parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
T_END = 100.0

corot = pd.read_csv(COROT_CSV)
inert = pd.read_csv(INERTIAL_CSV)
corot = corot[corot["Time [s]"] <= T_END].reset_index(drop=True)
inert = inert[inert["Time [s]"] <= T_END].reset_index(drop=True)


def rotor_equivalent_power(df: pd.DataFrame) -> pd.Series:
    """Return full-rotor-equivalent aerodynamic power.

    New rotor logs expose this convention explicitly. Older campaign logs store
    the single-blade value in ``Aero Power [W]``; use 3x to compare with rotor
    literature values.
    """
    if "Aero Power Rotor Equivalent [W]" in df:
        return df["Aero Power Rotor Equivalent [W]"]
    return 3.0 * df["Aero Power [W]"]


t1 = corot["Time [s]"]
flap1 = corot["Max Displacement [m]"]
t2 = inert["Time [s]"]
flap2 = inert["Max Displacement [m]"]

fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
axes[0].plot(t1, flap1, label="Corrotacional")
axes[0].plot(t2, flap2, label="Inercial", alpha=0.8)
axes[0].axhline(13.86, color="black", linestyle="--", linewidth=1.0, label="Zhou 2025 (flap media)")
axes[0].set_ylabel("Desplazamiento máximo [m]")
axes[0].set_title("Serie temporal FSI, yaw=0°")
axes[0].grid(alpha=0.3)
axes[0].legend(frameon=False)

axes[1].plot(t1, rotor_equivalent_power(corot) / 1e6, label="Corrotacional")
axes[1].plot(t2, rotor_equivalent_power(inert) / 1e6, label="Inercial", alpha=0.8)
axes[1].set_xlabel("Tiempo [s]")
axes[1].set_ylabel("Potencia aerodinámica rotor-eq. [MW]")
axes[1].grid(alpha=0.3)
axes[1].legend(frameon=False)

fig.tight_layout()

png = OUT_DIR / "fig_4_7_fsi_time_series.png"
pdf = OUT_DIR / "fig_4_7_fsi_time_series.pdf"
fig.savefig(png, dpi=300)
fig.savefig(pdf)
print(f"Saved: {png}")
print(f"Saved: {pdf}")
