"""Torque and power coupling analysis — IEA 15 MW, rated, yaw=0°.

Compares the aeroelastic rotor torque signal against a rigid BEM-only
reference at the same operating point. The objective is to preserve the
temporal characterization of torque while also quantifying, component by
component, how edgewise, flapwise, and resultant tip deformation relate to
both torque and power.

Outputs
-------
docs/figures/fig_5_4_7_torque_signal_analysis.png  (300 dpi)
docs/figures/fig_5_4_7_torque_signal_analysis.pdf
docs/figures/fig_5_4_7b_component_influence.pdf

stdout
------
Three summary tables:
    - rigid BEM reference and mean aeroelastic torque deficit
    - centered torque-flapwise coupling metrics in steady state
    - component-wise coupling metrics for torque and power
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import windows

BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
COROT_CSV = BASE / "frontiersin_results_corotational_100s/yaw_0/fluid/bem_report.csv"
INERT_CSV = BASE / "frontiersin_results_inertial/yaw_0/fluid/bem_report.csv"
OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

T_START = 40.0
T_END = 100.0
RATED_WIND_SPEED = 10.59
RATED_OMEGA_RPM = 7.55
RATED_PITCH_DEG = 0.0
RATED_YAW_DEG = 0.0
RATED_HUB_HEIGHT = 150.0
RATED_SHEAR_EXP = 0.2
RATED_RHO = 1.225
RATED_MU = 1.81206e-5
OMEGA_RAD_S = RATED_OMEGA_RPM * 2.0 * np.pi / 60.0

# Cached fallback computed from the same BEM-only operating point used here.
RIGID_BEM_TORQUE_FALLBACK = 20.081384915993045e6
RIGID_BEM_POWER_FALLBACK = RIGID_BEM_TORQUE_FALLBACK * OMEGA_RAD_S

TARGETS = {
    "torque": {"column": "Torque [N.m]", "scale": 1e6, "unit": "MN·m", "label": "Torque"},
    "power": {"column": "Power [W]", "scale": 1e6, "unit": "MW", "label": "Potencia"},
}

COMPONENTS = {
    "edgewise": {"column": "Tip Disp X [m]", "short": "X", "label": "Edgewise"},
    "flapwise": {"column": "Tip Disp Y [m]", "short": "Y", "label": "Flapwise"},
    "resultante": {"column": "Tip Disp Mag [m]", "short": "|U|", "label": "Resultante"},
}

COLORS = {
    "corot": "#c0392b",
    "inert": "#2980b9",
    "rigid": "#2c3e50",
    "tip": "#16a085",
    "edgewise": "#8e44ad",
    "flapwise": "#16a085",
    "resultante": "#f39c12",
}


def load_signal(path: pathlib.Path, t_start=None) -> pd.DataFrame:
    df = pd.read_csv(path)
    if t_start is not None:
        df = df[(df["Time [s]"] >= t_start) & (df["Time [s]"] <= T_END)].reset_index(drop=True)
    return df[
        [
            "Time [s]",
            "Torque [N.m]",
            "Power [W]",
            "Tip Disp X [m]",
            "Tip Disp Y [m]",
            "Tip Disp Mag [m]",
        ]
    ].copy()


def compute_rigid_bem_reference() -> tuple[dict[str, float], str]:
    """Return rigid BEM torque and power for the same nominal operating point.

    Falls back to a cached value when the local node cannot import the hybrid
    aeroelast runtime, which can happen on SDumont without the GCC14 runtime
    libraries preloaded.
    """

    try:
        from aeroelast.models.blade.aerodynamics import load_blade_aero
        from aeroelast.solvers.bem.engine import BEMSolver

        blade_yaml = REPO_ROOT / "tests/IEA-15-240-RWT.yaml"
        if not blade_yaml.exists():
            raise FileNotFoundError(blade_yaml)

        blade = load_blade_aero(str(blade_yaml))
        bem = BEMSolver(
            blade,
            rho=RATED_RHO,
            mu=RATED_MU,
            yaw=RATED_YAW_DEG,
            hub_height=RATED_HUB_HEIGHT,
            shear_exp=RATED_SHEAR_EXP,
        )
        result = bem.compute(
            v_inf=RATED_WIND_SPEED,
            omega=RATED_OMEGA_RPM,
            pitch=RATED_PITCH_DEG,
        )
        return {"torque": float(result.torque), "power": float(result.power)}, "computed"
    except Exception as exc:
        print(
            "Rigid BEM reference: using cached value because direct BEM execution "
            f"is unavailable on this node ({exc.__class__.__name__}: {exc})."
        )
        return {"torque": RIGID_BEM_TORQUE_FALLBACK, "power": RIGID_BEM_POWER_FALLBACK}, "cached"


def zscore(values: np.ndarray) -> np.ndarray:
    std = float(np.std(values, ddof=1))
    if std == 0.0:
        return np.zeros_like(values)
    return (values - float(np.mean(values))) / std


def onesided_fft(signal: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    signal = np.asarray(signal, dtype=float)
    n = len(signal)
    win = windows.hann(n)
    x = (signal - float(np.mean(signal))) * (win / win.mean())
    sp = np.fft.rfft(x)
    freq = np.fft.rfftfreq(n, d=dt)
    amp = np.abs(sp) * 2.0 / n
    amp[0] = 0.0
    return freq, amp


def dominant_frequency(signal: np.ndarray, dt: float, fmin=0.05, fmax=1.0) -> float:
    freq, amp = onesided_fft(signal, dt)
    mask = (freq >= fmin) & (freq <= fmax)
    if not np.any(mask):
        return float("nan")
    idx = int(np.argmax(amp[mask]))
    return float(freq[mask][idx])


def steady_stats(df: pd.DataFrame, q_rigid: float) -> dict[str, float]:
    q = df["Torque [N.m]"].to_numpy(dtype=float)
    y = df["Tip Disp Y [m]"].to_numpy(dtype=float)
    q_center = q - float(np.mean(q))
    y_center = y - float(np.mean(y))
    delta_q = q - q_rigid
    slope, intercept = np.polyfit(y_center, q_center, 1)
    return {
        "q_mean": float(np.mean(q)),
        "q_std": float(np.std(q, ddof=1)),
        "tip_mean": float(np.mean(y)),
        "tip_std": float(np.std(y, ddof=1)),
        "delta_q_mean": float(np.mean(delta_q)),
        "delta_q_pct": float(100.0 * np.mean(delta_q) / q_rigid),
        "delta_q_rms": float(np.sqrt(np.mean(delta_q**2))),
        "pearson": float(np.corrcoef(y_center, q_center)[0, 1]),
        "slope": float(slope),
        "intercept": float(intercept),
    }


def coupling_metrics(
    df: pd.DataFrame,
    target_column: str,
    displacement_column: str,
    target_scale: float,
    dt: float,
) -> dict[str, float]:
    target = df[target_column].to_numpy(dtype=float)
    disp = df[displacement_column].to_numpy(dtype=float)
    target_center = target - float(np.mean(target))
    disp_center = disp - float(np.mean(disp))
    slope = float(np.polyfit(disp_center, target_center, 1)[0])
    return {
        "pearson": float(np.corrcoef(disp_center, target_center)[0, 1]),
        "slope": slope,
        "slope_scaled": slope / target_scale,
        "target_dom_freq": dominant_frequency(target_center, dt),
        "disp_dom_freq": dominant_frequency(disp_center, dt),
    }


def strongest_component(metrics: dict[str, dict[str, float]], key: str) -> str:
    return max(metrics.items(), key=lambda item: abs(item[1][key]))[0]


def fit_line(x: np.ndarray, slope: float, intercept: float) -> tuple[np.ndarray, np.ndarray]:
    x_fit = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    y_fit = slope * x_fit + intercept
    return x_fit, y_fit


full = {"corot": load_signal(COROT_CSV), "inert": load_signal(INERT_CSV)}
steady = {
    name: load_signal(path, T_START)
    for name, path in {"corot": COROT_CSV, "inert": INERT_CSV}.items()
}

rigid_reference, rigid_source = compute_rigid_bem_reference()
q_rigid = rigid_reference["torque"]
stats = {name: steady_stats(df, q_rigid) for name, df in steady.items()}
dt = float(np.median(np.diff(steady["corot"]["Time [s]"].to_numpy(dtype=float))))
component_metrics: dict[str, dict[str, dict[str, dict[str, float]]]] = {}
for solver_name, df in steady.items():
    component_metrics[solver_name] = {}
    for target_name, target_cfg in TARGETS.items():
        component_metrics[solver_name][target_name] = {}
        for component_name, component_cfg in COMPONENTS.items():
            component_metrics[solver_name][target_name][component_name] = coupling_metrics(
                df,
                target_cfg["column"],
                component_cfg["column"],
                target_cfg["scale"],
                dt,
            )

print("\nRIGID BEM-ONLY REFERENCE (same Frontiers operating point)")
print(
    f"Torque = {rigid_reference['torque'] / 1e6:.3f} MN·m | "
    f"Power = {rigid_reference['power'] / 1e6:.3f} MW   source = {rigid_source}"
)

print("\nAEROELASTIC TORQUE VS RIGID BEM (t in [40, 100] s)")
print(
    f"{'Solver':<8} {'Q_mean':>9} {'Tip_mean':>10} {'ΔQ_mean':>10} {'ΔQ/Qrig':>10} {'Q_std':>9} {'Tip_std':>9}"
)
for name in ("corot", "inert"):
    item = stats[name]
    print(
        f"{name:<8} {item['q_mean'] / 1e6:>9.3f} {item['tip_mean']:>10.3f} "
        f"{item['delta_q_mean'] / 1e6:>10.3f} {item['delta_q_pct']:>9.2f}% "
        f"{item['q_std'] / 1e6:>9.3f} {item['tip_std']:>9.3f}"
    )

print("\nCENTERED TORQUE–FLAPWISE COUPLING (t in [40, 100] s)")
print(f"{'Solver':<8} {'Pearson r':>10} {'dQ/dy [MNm/m]':>16}")
for name in ("corot", "inert"):
    item = stats[name]
    print(f"{name:<8} {item['pearson']:>10.3f} {item['slope'] / 1e6:>16.3f}")

print("\nCOMPONENT-WISE COUPLING METRICS (t in [40, 100] s)")
print(
    f"{'Solver':<8} {'Target':<8} {'Component':<11} {'r':>8} {'d(target)/dU':>14} {'f_dom,target':>13} {'f_dom,U':>10}"
)
for solver_name in ("corot", "inert"):
    for target_name in ("torque", "power"):
        for component_name in ("edgewise", "flapwise", "resultante"):
            item = component_metrics[solver_name][target_name][component_name]
            unit = TARGETS[target_name]["unit"]
            print(
                f"{solver_name:<8} {target_name:<8} {component_name:<11} {item['pearson']:>8.3f} "
                f"{item['slope_scaled']:>14.3f} {item['target_dom_freq']:>13.3f} {item['disp_dom_freq']:>10.3f}"
            )

print("\nSTRONGEST COMPONENT BY METRIC")
for solver_name in ("corot", "inert"):
    for target_name in ("torque", "power"):
        metrics = component_metrics[solver_name][target_name]
        strongest_r = strongest_component(metrics, "pearson")
        strongest_slope = strongest_component(metrics, "slope_scaled")
        print(
            f"{solver_name} {target_name}: strongest |r| = {strongest_r}, strongest |slope| = {strongest_slope}"
        )


fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.0))
fig.suptitle(
    "Acoplamiento torque-deformación contra una referencia BEM rígida — IEA 15 MW, nominal, yaw=0°\n"
    "Respuesta aeroelástica corrotacional vs inercial en la ventana estacionaria t ∈ [40, 100] s",
    fontsize=11,
    y=1.02,
)

# (a) Full signal vs rigid BEM reference
ax = axes[0, 0]
for name, df in full.items():
    ax.plot(df["Time [s]"], df["Torque [N.m]"] / 1e6, lw=1.0, color=COLORS[name], label=name)
ax.axhline(q_rigid / 1e6, color=COLORS["rigid"], lw=1.2, ls=":", label="BEM rígido")
ax.axvline(T_START, color="#666", lw=1.0, ls="--")
ax.text(T_START + 0.6, ax.get_ylim()[1] - 0.12, "ventana estacionaria", fontsize=7, color="#555")
ax.set_title("(a) Torque aeroelástico vs línea rígida BEM")
ax.set_xlabel("Tiempo [s]")
ax.set_ylabel("Torque [MN·m]")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

# (b) Corotational coupling in z-scores
ax = axes[0, 1]
df = steady["corot"]
ax.plot(
    df["Time [s]"],
    zscore(df["Torque [N.m]"].to_numpy(float)),
    lw=1.0,
    color=COLORS["corot"],
    label="Torque",
)
ax.plot(
    df["Time [s]"],
    zscore(df["Tip Disp Y [m]"].to_numpy(float)),
    lw=1.0,
    color=COLORS["tip"],
    label="Punta flapwise Y",
)
ax.axhline(0.0, color="#666", lw=0.8, ls=":")
ax.set_title("(b) Corrotacional: torque y deformación normalizados")
ax.set_xlabel("Tiempo [s]")
ax.set_ylabel("z-score")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

# (c) Inertial coupling in z-scores
ax = axes[1, 0]
df = steady["inert"]
ax.plot(
    df["Time [s]"],
    zscore(df["Torque [N.m]"].to_numpy(float)),
    lw=1.0,
    color=COLORS["inert"],
    label="Torque",
)
ax.plot(
    df["Time [s]"],
    zscore(df["Tip Disp Y [m]"].to_numpy(float)),
    lw=1.0,
    color=COLORS["tip"],
    label="Punta flapwise Y",
)
ax.axhline(0.0, color="#666", lw=0.8, ls=":")
ax.set_title("(c) Inercial: torque y deformación normalizados")
ax.set_xlabel("Tiempo [s]")
ax.set_ylabel("z-score")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)

# (d) Centered torque-deflection coupling
ax = axes[1, 1]
summary_lines = []
for name in ("corot", "inert"):
    df = steady[name]
    y_center = df["Tip Disp Y [m]"].to_numpy(float) - stats[name]["tip_mean"]
    q_center = (df["Torque [N.m]"] - stats[name]["q_mean"]).to_numpy(float)
    ax.scatter(y_center, q_center / 1e6, s=8, alpha=0.10, color=COLORS[name])
    x_fit, y_fit = fit_line(y_center, stats[name]["slope"], stats[name]["intercept"])
    ax.plot(x_fit, y_fit / 1e6, lw=1.6, color=COLORS[name], label=name)
    summary_lines.append(
        f"{name}: ΔQ̄/Qrig = {stats[name]['delta_q_pct']:.2f}% | "
        f"r = {stats[name]['pearson']:.2f} | dQ/dy = {stats[name]['slope'] / 1e6:.3f}"
    )
ax.axhline(0.0, color="#666", lw=0.8, ls=":")
ax.axvline(0.0, color="#666", lw=0.8, ls=":")
ax.set_title("(d) Fluctuación de torque vs fluctuación flapwise")
ax.set_xlabel(r"$Y_{tip} - \bar{Y}_{tip}$ [m]")
ax.set_ylabel(r"$Q - \bar{Q}$ [MN·m]")
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=8)
ax.text(
    0.02,
    0.98,
    "\n".join(summary_lines),
    transform=ax.transAxes,
    va="top",
    ha="left",
    fontsize=8,
    bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#bbbbbb", "alpha": 0.92},
)

fig.tight_layout()
png = OUT_DIR / "fig_5_4_7_torque_signal_analysis.png"
pdf = OUT_DIR / "fig_5_4_7_torque_signal_analysis.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"\nSaved: {png}")
print(f"Saved: {pdf}")


fig2, axes2 = plt.subplots(2, 2, figsize=(13.2, 8.0))
fig2.suptitle(
    "Influencia separada de edgewise, flapwise y resultante sobre torque y potencia\n"
    "Comparación corrotacional vs inercial en la ventana estacionaria t ∈ [40, 100] s",
    fontsize=11,
    y=1.02,
)

component_names = list(COMPONENTS.keys())
component_labels = [COMPONENTS[name]["label"] for name in component_names]
x = np.arange(len(component_names))
width = 0.36


def add_grouped_bars(ax, target_name: str, metric_key: str, ylabel: str, title: str) -> None:
    corot_values = [
        component_metrics["corot"][target_name][name][metric_key] for name in component_names
    ]
    inert_values = [
        component_metrics["inert"][target_name][name][metric_key] for name in component_names
    ]
    ax.bar(x - width / 2, corot_values, width=width, color=COLORS["corot"], label="corot")
    ax.bar(x + width / 2, inert_values, width=width, color=COLORS["inert"], label="inert")
    ax.set_xticks(x, component_labels)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3, axis="y")


add_grouped_bars(
    axes2[0, 0],
    target_name="torque",
    metric_key="pearson",
    ylabel="r",
    title="(a) Correlación con torque",
)
add_grouped_bars(
    axes2[0, 1],
    target_name="torque",
    metric_key="slope_scaled",
    ylabel="dQ/dU [MN·m/m]",
    title="(b) Sensibilidad de torque",
)
add_grouped_bars(
    axes2[1, 0],
    target_name="power",
    metric_key="pearson",
    ylabel="r",
    title="(c) Correlación con potencia",
)
add_grouped_bars(
    axes2[1, 1],
    target_name="power",
    metric_key="slope_scaled",
    ylabel="dP/dU [MW/m]",
    title="(d) Sensibilidad de potencia",
)

for ax in axes2.flat:
    ax.legend(frameon=False, fontsize=8)

fig2.tight_layout()
png2 = OUT_DIR / "fig_5_4_7b_component_influence.png"
pdf2 = OUT_DIR / "fig_5_4_7b_component_influence.pdf"
fig2.savefig(png2, dpi=300, bbox_inches="tight")
fig2.savefig(pdf2, bbox_inches="tight")
print(f"Saved: {png2}")
print(f"Saved: {pdf2}")
