"""Yaw-sweep component influence analysis — IEA 15 MW, rated.

Quantifies how edgewise, flapwise, and resultant tip deformation influence
torque and power across the yaw sweep for the corotational and inertial
solvers. The script uses the same steady-state contract as the validation
report, truncates cases at the first non-physical sample, and exports a CSV so
the analysis can be refreshed when a new inertial campaign becomes available.

Outputs
-------
docs/figures/fig_5_4_7c_yaw_component_consistency.png
docs/figures/fig_5_4_7c_yaw_component_consistency.pdf
docs/figures/fig_5_4_7d_yaw_modal_closure.png
docs/figures/fig_5_4_7d_yaw_modal_closure.pdf
docs/figures/fig_5_4_7e_power_chain_closure.png
docs/figures/fig_5_4_7e_power_chain_closure.pdf
docs/validation_data/generated/yaw_component_influence.csv
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import coherence, csd, windows

BASE = Path("/scratch/leahk/eduardo.donestevez")
REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

ROOTS = {
    "corot": BASE / "frontiersin_results_corotational_100s",
    "inert": BASE / "frontiersin_results_inertial",
}
YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
MIN_STEADY_WINDOW = 20.0

RATED_WIND_SPEED = 10.59
RATED_OMEGA_RPM = 7.55
RATED_PITCH_DEG = 0.0
RATED_HUB_HEIGHT = 150.0
RATED_SHEAR_EXP = 0.2
RATED_RHO = 1.225
RATED_MU = 1.81206e-5
OMEGA_RAD_S = RATED_OMEGA_RPM * 2.0 * np.pi / 60.0

COMPONENTS = {
    "edgewise": {"column": "Tip Disp X [m]", "label": "Edgewise"},
    "flapwise": {"column": "Tip Disp Y [m]", "label": "Flapwise"},
    "resultante": {"column": "Tip Disp Mag [m]", "label": "Resultante"},
}

TARGETS = {
    "torque": {"column": "Torque [N.m]", "scale": 1e6, "unit": "MN·m"},
    "power": {"column": "Power [W]", "scale": 1e6, "unit": "MW"},
}

COLORS = {
    "edgewise": "#8e44ad",
    "flapwise": "#16a085",
    "resultante": "#f39c12",
}
MARKERS = {"corot": "o", "inert": "s"}


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


def spectral_coupling_metrics(
    target: np.ndarray,
    disp: np.ndarray,
    dt: float,
    fmin: float = 0.05,
    fmax: float = 1.0,
) -> dict[str, float]:
    target = np.asarray(target, dtype=float)
    disp = np.asarray(disp, dtype=float)
    n = len(target)
    if n < 32:
        return {
            "coupling_freq": float("nan"),
            "coherence_peak": float("nan"),
            "phase_deg": float("nan"),
            "cross_power": float("nan"),
        }

    fs = 1.0 / dt
    nperseg = min(4096, n)
    noverlap = nperseg // 2
    disp_center = disp - float(np.mean(disp))
    target_center = target - float(np.mean(target))
    freq, coh = coherence(
        disp_center,
        target_center,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
    )
    _, cross = csd(
        disp_center,
        target_center,
        fs=fs,
        window="hann",
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="constant",
    )
    mask = (freq >= fmin) & (freq <= fmax)
    if not np.any(mask):
        return {
            "coupling_freq": float("nan"),
            "coherence_peak": float("nan"),
            "phase_deg": float("nan"),
            "cross_power": float("nan"),
        }

    cross_abs = np.abs(cross[mask])
    idx = int(np.argmax(cross_abs))
    freq_sel = float(freq[mask][idx])
    coherence_sel = float(coh[mask][idx])
    phase_deg = float(np.rad2deg(np.angle(cross[mask][idx])))
    return {
        "coupling_freq": freq_sel,
        "coherence_peak": coherence_sel,
        "phase_deg": phase_deg,
        "cross_power": float(cross_abs[idx]),
    }


def compute_rigid_bem_reference(yaw_deg: float) -> tuple[dict[str, float], str]:
    try:
        from aeroelast.models.blade.aerodynamics import load_blade_aero
        from aeroelast.solvers.bem.engine import BEMSolver

        blade_yaml = REPO_ROOT / "tests/IEA-15-240-RWT.yaml"
        blade = load_blade_aero(str(blade_yaml))
        bem = BEMSolver(
            blade,
            rho=RATED_RHO,
            mu=RATED_MU,
            yaw=float(yaw_deg),
            hub_height=RATED_HUB_HEIGHT,
            shear_exp=RATED_SHEAR_EXP,
        )
        result = bem.compute(v_inf=RATED_WIND_SPEED, omega=RATED_OMEGA_RPM, pitch=RATED_PITCH_DEG)
        return {"torque": float(result.torque), "power": float(result.power)}, "computed"
    except Exception as exc:
        print(
            f"yaw={yaw_deg}: using cached rigid reference because direct BEM execution is unavailable "
            f"({exc.__class__.__name__}: {exc})"
        )
        torque = float("nan")
        power = float("nan")
        return {"torque": torque, "power": power}, "unavailable"


def load_case(csv_path: Path) -> tuple[pd.DataFrame | None, dict[str, float | str]]:
    cols = [
        "Torque [N.m]",
        "Power [W]",
        "Tip Disp X [m]",
        "Tip Disp Y [m]",
        "Tip Disp Mag [m]",
    ]
    df = pd.read_csv(csv_path)
    finite = np.isfinite(df[cols]).all(axis=1).to_numpy()
    physical = np.ones(len(df), dtype=bool)
    physical &= df["Tip Disp Mag [m]"].abs().to_numpy() < 50.0
    physical &= df["Torque [N.m]"].abs().to_numpy() < 1.0e8
    physical &= df["Power [W]"].abs().to_numpy() < 1.0e8
    valid = finite & physical

    first_bad = int(np.argmax(~valid)) if np.any(~valid) else len(df)
    truncated = df.iloc[:first_bad].copy() if first_bad < len(df) else df.copy()
    t_end_valid = float(truncated["Time [s]"].iloc[-1]) if not truncated.empty else float("nan")

    steady = truncated[
        (truncated["Time [s]"] >= T_START) & (truncated["Time [s]"] <= T_END)
    ].reset_index(drop=True)
    steady_window = min(t_end_valid, T_END) - T_START if not np.isnan(t_end_valid) else float("nan")
    status = "ok"
    if truncated.empty:
        status = "empty"
    elif steady.empty or steady_window < MIN_STEADY_WINDOW:
        status = "insufficient_steady_window"

    meta = {
        "rows_total": int(len(df)),
        "rows_valid": int(len(truncated)),
        "t_end_valid": t_end_valid,
        "steady_window": steady_window,
        "status": status,
    }
    if status != "ok":
        return None, meta
    return steady, meta


def coupling_metrics(
    df: pd.DataFrame, target_column: str, displacement_column: str, target_scale: float, dt: float
) -> dict[str, float]:
    target = df[target_column].to_numpy(dtype=float)
    disp = df[displacement_column].to_numpy(dtype=float)
    target_center = target - float(np.mean(target))
    disp_center = disp - float(np.mean(disp))
    slope = float(np.polyfit(disp_center, target_center, 1)[0])
    metrics = {
        "pearson": float(np.corrcoef(disp_center, target_center)[0, 1]),
        "slope": slope,
        "slope_scaled": slope / target_scale,
        "target_mean": float(np.mean(target)) / target_scale,
        "target_std": float(np.std(target, ddof=1)) / target_scale,
        "disp_mean": float(np.mean(disp)),
        "disp_std": float(np.std(disp, ddof=1)),
        "target_dom_freq": dominant_frequency(target_center, dt),
        "disp_dom_freq": dominant_frequency(disp_center, dt),
    }
    metrics.update(spectral_coupling_metrics(target, disp, dt))
    return metrics


rows: list[dict[str, float | int | str]] = []
for solver_name, root in ROOTS.items():
    for yaw in YAWS:
        csv_path = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
        if not csv_path.exists():
            rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                "target": "torque",
                "component": "edgewise",
                "status": "missing_csv",
            })
            continue

        steady_df, meta = load_case(csv_path)
        rigid_ref, rigid_source = compute_rigid_bem_reference(yaw)
        if steady_df is None:
            rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                "target": "torque",
                "component": "edgewise",
                "status": meta["status"],
                "rows_total": meta["rows_total"],
                "rows_valid": meta["rows_valid"],
                "t_end_valid": meta["t_end_valid"],
            })
            continue

        dt = float(np.median(np.diff(steady_df["Time [s]"].to_numpy(dtype=float))))
        for target_name, target_cfg in TARGETS.items():
            target_mean_scaled = (
                float(np.mean(steady_df[target_cfg["column"]].to_numpy(dtype=float)))
                / target_cfg["scale"]
            )
            rigid_scaled = (
                rigid_ref[target_name] / target_cfg["scale"]
                if np.isfinite(rigid_ref[target_name])
                else float("nan")
            )
            delta_mean_scaled = (
                target_mean_scaled - rigid_scaled if np.isfinite(rigid_scaled) else float("nan")
            )
            delta_mean_pct = (
                100.0 * delta_mean_scaled / rigid_scaled
                if np.isfinite(rigid_scaled) and rigid_scaled != 0.0
                else float("nan")
            )

            for component_name, component_cfg in COMPONENTS.items():
                metrics = coupling_metrics(
                    steady_df,
                    target_cfg["column"],
                    component_cfg["column"],
                    target_cfg["scale"],
                    dt,
                )
                rows.append({
                    "solver": solver_name,
                    "yaw_deg": yaw,
                    "target": target_name,
                    "component": component_name,
                    "status": meta["status"],
                    "rows_total": meta["rows_total"],
                    "rows_valid": meta["rows_valid"],
                    "steady_window": meta["steady_window"],
                    "t_end_valid": meta["t_end_valid"],
                    "rigid_source": rigid_source,
                    "rigid_reference": rigid_scaled,
                    "target_mean": target_mean_scaled,
                    "delta_mean": delta_mean_scaled,
                    "delta_mean_pct": delta_mean_pct,
                    **metrics,
                })

results = pd.DataFrame(rows)
csv_out = DATA_DIR / "yaw_component_influence.csv"
results.to_csv(csv_out, index=False)

valid = results[results["status"] == "ok"].copy()


def wrapped_phase_diff_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return (a - b + 180.0) % 360.0 - 180.0


print("\nYAW-SWEEP COMPONENT INFLUENCE")
for solver_name in ("corot", "inert"):
    print(f"\nsolver={solver_name}")
    solver_df = valid[(valid["solver"] == solver_name) & (valid["target"] == "torque")]
    for yaw in sorted(solver_df["yaw_deg"].unique()):
        case = solver_df[solver_df["yaw_deg"] == yaw]
        best_r = case.iloc[case["pearson"].abs().argmax()]
        best_s = case.iloc[case["slope_scaled"].abs().argmax()]
        print(
            f"  yaw={yaw:>2}: best |r|={best_r['component']} ({best_r['pearson']:+.3f}) | "
            f"best |dQ/dU|={best_s['component']} ({best_s['slope_scaled']:+.3f} MNm/m)"
        )

print("\nYAW-SWEEP MODAL COUPLING (TORQUE vs DEFORMATION)")
for solver_name in ("corot", "inert"):
    print(f"\nsolver={solver_name}")
    solver_df = valid[(valid["solver"] == solver_name) & (valid["target"] == "torque")]
    for yaw in sorted(solver_df["yaw_deg"].unique()):
        case = solver_df[solver_df["yaw_deg"] == yaw]
        best_c = case.iloc[case["cross_power"].argmax()]
        print(
            f"  yaw={yaw:>2}: strongest copower={best_c['component']} (C={best_c['coherence_peak']:.3f}, "
            f"f={best_c['coupling_freq']:.3f} Hz, phase={best_c['phase_deg']:+.1f} deg)"
        )

skipped = results[results["status"] != "ok"]
if not skipped.empty:
    print("\nSKIPPED CASES")
    print(
        skipped[["solver", "yaw_deg", "status", "t_end_valid"]]
        .drop_duplicates()
        .to_string(index=False)
    )

torque_df = valid[valid["target"] == "torque"].copy()
power_df = valid[valid["target"] == "power"].copy()
closure = torque_df.merge(
    power_df,
    on=["solver", "yaw_deg", "component"],
    suffixes=("_torque", "_power"),
)
closure["slope_ratio"] = closure["slope_scaled_power"] / closure["slope_scaled_torque"]
closure["pearson_diff"] = closure["pearson_power"] - closure["pearson_torque"]
closure["freq_diff"] = closure["coupling_freq_power"] - closure["coupling_freq_torque"]
closure["phase_diff"] = wrapped_phase_diff_deg(
    closure["phase_deg_power"].to_numpy(dtype=float),
    closure["phase_deg_torque"].to_numpy(dtype=float),
)

print("\nPOWER-CHAIN CLOSURE (POWER vs TORQUE METRICS)")
print(
    f"dP/dU over dQ/dU: mean={closure['slope_ratio'].mean():.6f} rad/s | "
    f"std={closure['slope_ratio'].std(ddof=1):.3e} | expected Omega={OMEGA_RAD_S:.6f} rad/s"
)
print(
    f"max |Δpearson|={closure['pearson_diff'].abs().max():.3e} | "
    f"max |Δf_c|={closure['freq_diff'].abs().max():.3e} Hz | "
    f"max |Δphase|={np.abs(closure['phase_diff']).max():.3e} deg"
)

fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.4), sharex=True)
fig.suptitle(
    "Consistencia del acoplamiento por componente a lo largo del barrido yaw\n"
    "Torque: correlación y sensibilidad para edgewise, flapwise y resultante",
    fontsize=11,
    y=1.02,
)

panel_map = {
    (0, 0): ("corot", "pearson", "(a) Corrotacional: correlación con torque", "r"),
    (0, 1): ("inert", "pearson", "(b) Inercial: correlación con torque", "r"),
    (1, 0): (
        "corot",
        "slope_scaled",
        "(c) Corrotacional: sensibilidad de torque",
        "dQ/dU [MN·m/m]",
    ),
    (1, 1): ("inert", "slope_scaled", "(d) Inercial: sensibilidad de torque", "dQ/dU [MN·m/m]"),
}

for (i, j), (solver_name, metric_key, title, ylabel) in panel_map.items():
    ax = axes[i, j]
    case_df = valid[(valid["solver"] == solver_name) & (valid["target"] == "torque")]
    for component_name, component_cfg in COMPONENTS.items():
        comp_df = case_df[case_df["component"] == component_name].sort_values("yaw_deg")
        ax.plot(
            comp_df["yaw_deg"],
            comp_df[metric_key],
            marker="o",
            lw=1.5,
            color=COLORS[component_name],
            label=component_cfg["label"],
        )
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

for ax in axes[1, :]:
    ax.set_xlabel("Yaw [°]")

fig.tight_layout()
png_out = OUT_DIR / "fig_5_4_7c_yaw_component_consistency.png"
pdf_out = OUT_DIR / "fig_5_4_7c_yaw_component_consistency.pdf"
fig.savefig(png_out, dpi=300, bbox_inches="tight")
fig.savefig(pdf_out, bbox_inches="tight")

fig2, axes2 = plt.subplots(2, 2, figsize=(13.5, 8.4), sharex=True)
fig2.suptitle(
    "Cierre modal del acoplamiento torque-deformacion en el barrido yaw\n"
    "Frecuencia y fase en la frecuencia de maximo copoder dentro de 0.05-1.0 Hz",
    fontsize=11,
    y=1.02,
)

phase_panel_map = {
    (0, 0): (
        "corot",
        "coupling_freq",
        "(a) Corrotacional: frecuencia de maximo copoder",
        "f_c [Hz]",
    ),
    (0, 1): ("inert", "coupling_freq", "(b) Inercial: frecuencia de maximo copoder", "f_c [Hz]"),
    (1, 0): (
        "corot",
        "phase_deg",
        "(c) Corrotacional: fase torque respecto a deformacion",
        "Fase [deg]",
    ),
    (1, 1): (
        "inert",
        "phase_deg",
        "(d) Inercial: fase torque respecto a deformacion",
        "Fase [deg]",
    ),
}

for (i, j), (solver_name, metric_key, title, ylabel) in phase_panel_map.items():
    ax = axes2[i, j]
    case_df = valid[(valid["solver"] == solver_name) & (valid["target"] == "torque")]
    for component_name, component_cfg in COMPONENTS.items():
        comp_df = case_df[case_df["component"] == component_name].sort_values("yaw_deg")
        ax.plot(
            comp_df["yaw_deg"],
            comp_df[metric_key],
            marker="o",
            lw=1.5,
            color=COLORS[component_name],
            label=component_cfg["label"],
        )
    if metric_key == "phase_deg":
        ax.axhline(0.0, color="#666", lw=0.9, ls=":")
    if metric_key == "coupling_freq":
        ax.axhline(0.122, color="#666", lw=0.9, ls=":")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, fontsize=8)

for ax in axes2[1, :]:
    ax.set_xlabel("Yaw [°]")

fig2.tight_layout()
png_out_phase = OUT_DIR / "fig_5_4_7d_yaw_modal_closure.png"
pdf_out_phase = OUT_DIR / "fig_5_4_7d_yaw_modal_closure.pdf"
fig2.savefig(png_out_phase, dpi=300, bbox_inches="tight")
fig2.savefig(pdf_out_phase, bbox_inches="tight")

print(f"\nSaved: {png_out}")
print(f"Saved: {pdf_out}")
print(f"Saved: {png_out_phase}")
print(f"Saved: {pdf_out_phase}")

fig3, axes3 = plt.subplots(2, 2, figsize=(11.8, 9.0))
fig3.suptitle(
    "Cierre de la cadena deformacion -> torque -> potencia\n"
    "Verificacion numerica de que la potencia hereda la misma estructura de acoplamiento cuando Omega es constante",
    fontsize=11,
    y=1.02,
)

scatter_specs = {
    (0, 0): (
        "slope_scaled_torque",
        "slope_scaled_power",
        "(a) Sensibilidad: dP/dU vs dQ/dU",
        "dQ/dU [MN·m/m]",
        "dP/dU [MW/m]",
        "omega",
    ),
    (0, 1): (
        "pearson_torque",
        "pearson_power",
        "(b) Correlacion: potencia vs torque",
        "r con torque",
        "r con potencia",
        "identity",
    ),
    (1, 0): (
        "coupling_freq_torque",
        "coupling_freq_power",
        "(c) Frecuencia de acoplamiento: potencia vs torque",
        "f_c con torque [Hz]",
        "f_c con potencia [Hz]",
        "identity",
    ),
    (1, 1): (
        "phase_deg_torque",
        "phase_deg_power",
        "(d) Fase: potencia vs torque",
        "Fase con torque [deg]",
        "Fase con potencia [deg]",
        "identity",
    ),
}

for (i, j), (xkey, ykey, title, xlabel, ylabel, ref_kind) in scatter_specs.items():
    ax = axes3[i, j]
    for solver_name in ("corot", "inert"):
        for component_name, component_cfg in COMPONENTS.items():
            subset = closure[
                (closure["solver"] == solver_name) & (closure["component"] == component_name)
            ].sort_values("yaw_deg")
            ax.scatter(
                subset[xkey],
                subset[ykey],
                s=42,
                marker=MARKERS[solver_name],
                color=COLORS[component_name],
                alpha=0.9,
                label=f"{solver_name}-{component_cfg['label']}",
            )

    xvals = closure[xkey].to_numpy(dtype=float)
    yvals = closure[ykey].to_numpy(dtype=float)
    lo = float(np.nanmin(np.concatenate([xvals, yvals])))
    hi = float(np.nanmax(np.concatenate([xvals, yvals])))
    if ref_kind == "identity":
        ref_x = np.linspace(lo, hi, 200)
        ref_y = ref_x
    else:
        ref_x = np.linspace(float(np.nanmin(xvals)), float(np.nanmax(xvals)), 200)
        ref_y = OMEGA_RAD_S * ref_x
    ax.plot(ref_x, ref_y, color="#2c3e50", lw=1.1, ls="--")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)

handles, labels = axes3[0, 0].get_legend_handles_labels()
unique = dict(zip(labels, handles))
axes3[0, 0].legend(unique.values(), unique.keys(), frameon=False, fontsize=7, ncol=2)

fig3.tight_layout()
png_out_power = OUT_DIR / "fig_5_4_7e_power_chain_closure.png"
pdf_out_power = OUT_DIR / "fig_5_4_7e_power_chain_closure.pdf"
fig3.savefig(png_out_power, dpi=300, bbox_inches="tight")
fig3.savefig(pdf_out_power, bbox_inches="tight")

print(f"Saved: {png_out_power}")
print(f"Saved: {pdf_out_power}")
print(f"Saved: {csv_out}")
