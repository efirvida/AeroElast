"""Generate corotational-only figures for the standalone rotor validation report.

This script avoids touching the mixed comparison figures used by the master
validation dossier. It reuses the current corotational campaign plus the
generated CSVs already exported by the mechanism-analysis pipeline.
"""

from __future__ import annotations

import base64
import mmap
import pathlib
import re
import zlib

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
DOCS_DIR = pathlib.Path(__file__).resolve().parents[1]
OUT_DIR = DOCS_DIR / "figures"
DATA_DIR = DOCS_DIR / "validation_data"
GENERATED_DIR = DATA_DIR / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COROT_ROOT = BASE / "frontiersin_results_corotational_100s"
MA_FIG15 = DATA_DIR / "ma_2025_fig15_yaw_metrics.csv"
MA_FIG16 = DATA_DIR / "ma_2025_fig16_spanwise_deflections.csv"

T_START = 40.0
OMEGA_RPM = 7.55
OMEGA_RAD_S = OMEGA_RPM * 2.0 * np.pi / 60.0
RIGID_BEM_TORQUE_MNM = 20.081384915993045
RIGID_BEM_POWER_MW = RIGID_BEM_TORQUE_MNM * OMEGA_RAD_S

F_1P = 7.518 / 60.0
F_3P = 3.0 * F_1P
F_FLAP = 0.5537
F_EDGE = 0.6290

PALETTE = {
    "corot": "#c0392b",
    "ma": "#16a085",
    "zhou": "#8e44ad",
    "ref": "#2c3e50",
    "neutral": "#7f8c8d",
    "edgewise": "#8e44ad",
    "flapwise": "#16a085",
    "resultante": "#f39c12",
    "1P": "#e74c3c",
    "3P": "#3498db",
    "flap": "#27ae60",
    "edge": "#f39c12",
    "residual": "#d0d3d4",
}

COMPONENT_ORDER = ["edgewise", "flapwise", "resultante"]
COMPONENT_LABELS = {
    "edgewise": "Edgewise",
    "flapwise": "Flapwise",
    "resultante": "Resultante",
}
SIGNAL_LABELS = {
    "thrust": "Thrust",
    "torque": "Torque",
    "power": "Potencia",
    "edgewise": "Edgewise",
    "flapwise": "Flapwise",
    "resultante": "Resultante",
}
YAWS = [0, 10, 20, 30, 40]


def save_figure(fig: plt.Figure, basename: str) -> None:
    png = OUT_DIR / f"{basename}.png"
    pdf = OUT_DIR / f"{basename}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> {png.name}")


def load_trim(csv_path: pathlib.Path, t_start: float = T_START) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    return df[df["Time [s]"] >= t_start].reset_index(drop=True)


def load_full(csv_path: pathlib.Path) -> pd.DataFrame:
    return pd.read_csv(csv_path)


def onesided_fft(values: pd.Series | np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    arr = arr - float(np.nanmean(arr))
    win = np.hanning(len(arr))
    arr = arr * (win / win.mean())
    spectrum = np.fft.rfft(arr)
    freq = np.fft.rfftfreq(len(arr), d=dt)
    amp = np.abs(spectrum) * 2.0 / len(arr)
    amp[0] = 0.0
    return freq, amp


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    if isinstance(values, np.ndarray):
        arr = values.astype(float, copy=False)
    else:
        arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    std = float(np.std(arr, ddof=1))
    if std == 0.0:
        return np.zeros_like(arr)
    return (arr - float(np.mean(arr))) / std


def save_df_plot_data(df: pd.DataFrame, path: pathlib.Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    raise FileNotFoundError(path)


def summarize_corotational_campaign() -> pd.DataFrame:
    rows = []
    for yaw in YAWS:
        csv_path = COROT_ROOT / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
        df = load_trim(csv_path)
        rows.append({
            "yaw_deg": yaw,
            "n_samples": int(len(df)),
            "t_end_s": float(df["Time [s]"].max()),
            "cp": float(df["CP"].mean()),
            "ct": float(df["CT"].mean()),
            "power_mw": float(df["Power [W]"].mean() / 1.0e6),
            "power_std_mw": float(df["Power [W]"].std() / 1.0e6),
            "thrust_mn": float(df["Thrust [N]"].mean() / 1.0e6),
            "torque_mnm": float(df["Torque [N.m]"].mean() / 1.0e6),
            "flap_m": float(df["Tip Disp Y [m]"].mean()),
            "flap_std_m": float(df["Tip Disp Y [m]"].std()),
            "flap_p2p_m": float(df["Tip Disp Y [m]"].max() - df["Tip Disp Y [m]"].min()),
            "edge_m": float(df["Tip Disp X [m]"].mean()),
            "edge_std_m": float(df["Tip Disp X [m]"].std()),
        })
    return pd.DataFrame(rows)


def component_color(component: str) -> str:
    return PALETTE[component]


def component_label(component: str) -> str:
    return COMPONENT_LABELS[component]


def extract_ma_tip_series(column_prefix: str) -> tuple[np.ndarray, np.ndarray]:
    ma = pd.read_csv(MA_FIG16)
    tip = ma.loc[np.isclose(ma["r_R"], 1.0)].iloc[0]
    yaws = np.array([10, 20, 30, 40], dtype=float)
    values = np.array([tip[f"{column_prefix}_{int(yaw)}deg_m"] for yaw in yaws], dtype=float)
    return yaws, values


def extract_ma_yaw_metrics() -> pd.DataFrame:
    return pd.read_csv(MA_FIG15)


def _time_value(path: pathlib.Path) -> float:
    return float(path.parent.name)


def _field_vtu_paths(case_dir: pathlib.Path) -> list[pathlib.Path]:
    paths = [path for path in case_dir.glob("*/fields.vtu") if path.parent.name != "fluid"]
    return sorted(paths, key=_time_value)


def _decode_vtu_payload(encoded: bytes, dtype: np.dtype) -> np.ndarray:
    payload = base64.b64decode(encoded)
    if len(payload) >= 12:
        n_blocks, _, _ = np.frombuffer(payload, dtype="<u4", count=3)
        header_size = 4 * (3 + int(n_blocks))
        if 0 < n_blocks < 10000 and header_size <= len(payload):
            compressed_sizes = np.frombuffer(
                payload, dtype="<u4", count=int(n_blocks), offset=12
            ).astype(int)
            if header_size + int(compressed_sizes.sum()) == len(payload):
                chunks = []
                offset = header_size
                for size in compressed_sizes:
                    chunks.append(zlib.decompress(payload[offset : offset + int(size)]))
                    offset += int(size)
                return np.frombuffer(b"".join(chunks), dtype=dtype)

    n_bytes = int(np.frombuffer(payload, dtype="<u4", count=1)[0])
    if 4 + n_bytes <= len(payload):
        return np.frombuffer(payload[4 : 4 + n_bytes], dtype=dtype)
    return np.frombuffer(payload, dtype=dtype)


def _vtu_data_array(path: pathlib.Path, name: str) -> np.ndarray:
    name_pattern = f'Name="{name}"'.encode()
    with path.open("rb") as handle:
        with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            name_start = mapped.find(name_pattern)
            if name_start == -1:
                raise KeyError(f"VTU {path} does not contain DataArray {name!r}")
            tag_start = mapped.rfind(b"<DataArray", 0, name_start)
            tag_end = mapped.find(b">", name_start)
            data_end = mapped.find(b"</DataArray>", tag_end)
            if tag_start == -1 or tag_end == -1 or data_end == -1:
                raise ValueError(f"Malformed DataArray {name!r} in {path}")
            tag = mapped[tag_start : tag_end + 1].decode("ascii", errors="ignore")
            encoded = mapped[tag_end + 1 : data_end].strip()

    vtk_type = re.search(r'type="([^"]+)"', tag)
    if vtk_type is None or vtk_type.group(1) != "Float64":
        raise TypeError(f"Unsupported VTU array type for {name!r}")
    n_components_match = re.search(r'NumberOfComponents="(\d+)"', tag)
    n_components = int(n_components_match.group(1)) if n_components_match else 1
    values = _decode_vtu_payload(encoded, np.dtype("<f8"))
    if n_components > 1:
        values = values.reshape((-1, n_components))
    return values.astype(float, copy=False)


def _select_tip_node_from_vtu(path: pathlib.Path) -> tuple[int, float, int]:
    points = _vtu_data_array(path, "Points")
    displacement = _vtu_data_array(path, "U")
    points0 = points - displacement
    z0 = points0[:, 2]
    z_tip = float(z0.max())
    candidates = np.where(np.isclose(z0, z_tip, rtol=0.0, atol=1e-8))[0]
    if len(candidates) == 1:
        return int(candidates[0]), z_tip, 1
    centroid_xy = points0[candidates, :2].mean(axis=0)
    distances = np.sum((points0[candidates, :2] - centroid_xy) ** 2, axis=1)
    return int(candidates[int(np.argmin(distances))]), z_tip, int(len(candidates))


def load_tip_vtu_series(case_dir: pathlib.Path, t_start: float = T_START) -> pd.DataFrame:
    paths = [path for path in _field_vtu_paths(case_dir) if _time_value(path) >= t_start]
    if not paths:
        raise FileNotFoundError(f"No fields.vtu files found under {case_dir}")
    tip_node, z_tip, n_candidates = _select_tip_node_from_vtu(paths[0])
    rows = []
    for path in paths:
        displacement = _vtu_data_array(path, "U")[tip_node]
        rows.append({
            "Time [s]": _time_value(path),
            "Tip Disp X [m]": displacement[0],
            "Tip Disp Y [m]": displacement[1],
            "Tip Disp Z [m]": displacement[2],
            "tip_node": tip_node,
            "tip_z0": z_tip,
            "tip_z0_candidates": n_candidates,
        })
    return pd.DataFrame(rows).sort_values("Time [s]").reset_index(drop=True)


def _representative_station_nodes(path: pathlib.Path, r_targets: np.ndarray) -> np.ndarray:
    points = _vtu_data_array(path, "Points")
    displacement = _vtu_data_array(path, "U")
    points0 = points - displacement
    z0 = points0[:, 2]
    z0 = z0 - float(z0.min())
    span = float(z0.max())
    r = z0 / span
    station_nodes = []
    for target in np.asarray(r_targets, dtype=float):
        distances = np.abs(r - target)
        tolerance = float(distances.min()) + 1e-12
        candidates = np.where(distances <= tolerance)[0]
        centroid_xy = points0[candidates, :2].mean(axis=0)
        distances_xy = np.sum((points0[candidates, :2] - centroid_xy) ** 2, axis=1)
        station_nodes.append(int(candidates[int(np.argmin(distances_xy))]))
    return np.asarray(station_nodes, dtype=int)


def _spanwise_time_upper_bound(case_dir: pathlib.Path, t_start: float = T_START) -> float:
    report_path = case_dir / "structural_report.csv"
    if not report_path.exists():
        return float("inf")
    report = pd.read_csv(report_path)
    if "Time [s]" not in report.columns or "Max Disp Y [m]" not in report.columns:
        return float("inf")
    report = report[report["Time [s]"] >= t_start].reset_index(drop=True)
    if report.empty:
        return float("inf")
    flap = pd.to_numeric(report["Max Disp Y [m]"], errors="coerce")
    valid = np.isfinite(flap) & (np.abs(flap) <= 20.0)
    if valid.all():
        return float("inf")
    invalid_idx = np.flatnonzero((~valid).to_numpy())
    if len(invalid_idx) == 0 or invalid_idx[0] == 0:
        return float("inf")
    return float(report.loc[int(invalid_idx[0] - 1), "Time [s]"])


def load_spanwise_vtu_profile(
    case_dir: pathlib.Path, r_targets: np.ndarray, t_start: float = T_START, stride: int = 50
) -> dict[str, np.ndarray | float | int]:
    t_end = _spanwise_time_upper_bound(case_dir, t_start)
    paths = [
        path
        for path in _field_vtu_paths(case_dir)
        if _time_value(path) >= t_start and _time_value(path) <= t_end + 1e-12
    ]
    if not paths:
        raise FileNotFoundError(f"No stationary fields.vtu files found under {case_dir}")
    sampled_paths = paths[::stride] if stride > 1 else list(paths)
    if sampled_paths[-1] != paths[-1]:
        sampled_paths.append(paths[-1])
    station_nodes = _representative_station_nodes(paths[0], r_targets)
    flap_sum = np.zeros(len(r_targets), dtype=float)
    edge_sum = np.zeros(len(r_targets), dtype=float)
    for path in sampled_paths:
        displacement = _vtu_data_array(path, "U")
        flap_sum += displacement[station_nodes, 1]
        edge_sum += displacement[station_nodes, 0]
    flap_mean = flap_sum / len(sampled_paths)
    edge_mean = edge_sum / len(sampled_paths)
    bem_path = case_dir.parent / "fluid" / "bem_report.csv"
    if bem_path.exists():
        bem_trim = load_trim(bem_path)
        tip_edge_bem = float(bem_trim["Tip Disp X [m]"].mean())
        tip_edge_vtu = float(edge_mean[-1])
        bem_sign = float(np.sign(tip_edge_bem))
        vtu_sign = float(np.sign(tip_edge_vtu))
        if bem_sign != 0.0 and vtu_sign != 0.0:
            edge_mean = edge_mean * (bem_sign * vtu_sign)
    return {
        "r_R": np.asarray(r_targets, dtype=float),
        "flap": flap_mean,
        "edge": edge_mean,
        "n_frames": int(len(sampled_paths)),
    }


def plot_timeseries_deflection(summary: pd.DataFrame) -> None:
    csv_path = COROT_ROOT / "yaw_0" / "fluid" / "bem_report.csv"
    full = load_full(csv_path)
    steady = full[full["Time [s]"] >= T_START].reset_index(drop=True)
    zhou_ref = {"Tip Disp Y [m]": 13.86, "Tip Disp X [m]": -1.22}
    labels = {
        "Tip Disp Y [m]": "(a) Flapwise tip deflection",
        "Tip Disp X [m]": "(b) Edgewise tip deflection",
    }
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.4), sharex=True)
    for ax, column in zip(axes, ["Tip Disp Y [m]", "Tip Disp X [m]"]):
        ax.axvspan(T_START, float(full["Time [s]"].max()), color="#000000", alpha=0.04)
        ax.plot(
            full["Time [s]"], full[column], color=PALETTE["corot"], lw=1.1, label="Corrotacional"
        )
        ax.axhline(
            float(steady[column].mean()),
            color=PALETTE["corot"],
            lw=1.0,
            ls="--",
            label="Media t>=40 s",
        )
        ax.axhline(zhou_ref[column], color=PALETTE["zhou"], lw=1.0, ls=":", label="Zhou 2025")
        ax.axvline(T_START, color=PALETTE["ref"], lw=0.9, ls="-.")
        ax.set_ylabel("Displacement [m]")
        ax.set_title(labels[column])
        ax.grid(alpha=0.3, lw=0.5)
    handles = [
        Line2D([0], [0], color=PALETTE["corot"], lw=1.1, label="Corrotacional"),
        Line2D([0], [0], color=PALETTE["corot"], lw=1.0, ls="--", label="Media t>=40 s"),
        Line2D([0], [0], color=PALETTE["zhou"], lw=1.0, ls=":", label="Zhou 2025"),
        Line2D([0], [0], color=PALETTE["ref"], lw=0.9, ls="-.", label="Inicio ventana estadistica"),
    ]
    axes[0].legend(handles=handles, frameon=False, ncol=4, loc="upper center")
    axes[-1].set_xlabel("Time [s]")
    fig.suptitle("Tip-deflection history - corotational campaign, yaw=0", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    save_figure(fig, "fig_5_1_timeseries_deflection_corotational")


def plot_yaw_sweep(summary: pd.DataFrame) -> None:
    ma_metrics = extract_ma_yaw_metrics()
    _, ma_tip_flap = extract_ma_tip_series("flap")
    yaw = summary["yaw_deg"].to_numpy(dtype=float)
    yaw_fine = np.linspace(0.0, 40.0, 200)
    cos3 = np.cos(np.radians(yaw_fine)) ** 3
    cp0 = float(summary.loc[summary["yaw_deg"] == 0, "cp"].iloc[0])
    p0 = float(summary.loc[summary["yaw_deg"] == 0, "power_mw"].iloc[0])
    ma_p0 = float(ma_metrics["mean_power_MW"].iloc[0] / (np.cos(np.radians(10.0)) ** 1.5))
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))

    ax = axes[0]
    ax.plot(
        yaw_fine, cp0 * cos3, color=PALETTE["neutral"], lw=1.4, ls="--", label="CP0 cos^3(gamma)"
    )
    ax.plot(yaw, summary["cp"], color=PALETTE["corot"], lw=1.2, marker="o", label="Corrotacional")
    ax.set_title("(a) CP vs yaw")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("CP [-]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(yaw_fine, p0 * cos3, color=PALETTE["neutral"], lw=1.3, ls="--", label="P0 cos^3(gamma)")
    ax.plot(
        yaw_fine,
        ma_p0 * np.cos(np.radians(yaw_fine)) ** 1.5,
        color=PALETTE["ma"],
        lw=1.1,
        ls="-.",
        label="Ma fit cos^1.5",
    )
    ax.plot(
        yaw, summary["power_mw"], color=PALETTE["corot"], lw=1.2, marker="o", label="Corrotacional"
    )
    ax.scatter(
        ma_metrics["yaw_angle_deg"],
        ma_metrics["mean_power_MW"],
        color=PALETTE["ma"],
        s=48,
        marker="D",
        zorder=4,
        label="Ma 2025",
    )
    ax.set_title("(b) Mean power")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("Power [MW]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    ax = axes[2]
    ax.errorbar(
        yaw,
        summary["flap_m"],
        yerr=summary["flap_p2p_m"] / 2.0,
        color=PALETTE["corot"],
        lw=1.1,
        marker="o",
        capsize=3,
        label="Corrotacional",
    )
    ax.scatter(
        np.array([10, 20, 30, 40], dtype=float),
        ma_tip_flap,
        color=PALETTE["ma"],
        marker="D",
        s=48,
        label="Ma 2025",
    )
    ax.set_title("(c) Tip flapwise")
    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("Tip flapwise [m]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    for ax in axes:
        ax.set_xlim(-1.5, 41.5)
    fig.suptitle("Corotational yaw sweep with literature references", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_2_yaw_sweep_corotational")


def plot_literature_comparison(summary: pd.DataFrame) -> None:
    ma_metrics = extract_ma_yaw_metrics()
    ma_yaw = ma_metrics["yaw_angle_deg"].to_numpy(dtype=float)
    _, ma_tip_flap = extract_ma_tip_series("flap")
    _, ma_tip_edge = extract_ma_tip_series("edge")
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))

    series = [
        ("flap_m", ma_tip_flap, 13.86, "Tip flapwise [m]", "(a) Flapwise tip"),
        ("edge_m", ma_tip_edge, -1.22, "Tip edgewise [m]", "(b) Edgewise tip"),
        (
            "thrust_mn",
            ma_metrics["mean_thrust_1e6N"].to_numpy(dtype=float),
            2.20,
            "Thrust [MN]",
            "(c) Mean thrust",
        ),
    ]
    for ax, (column, ma_values, zhou_y0, ylabel, title) in zip(axes, series):
        ax.plot(
            summary["yaw_deg"],
            summary[column],
            color=PALETTE["corot"],
            lw=1.2,
            marker="o",
            label="Corrotacional",
        )
        ax.scatter(ma_yaw, ma_values, color=PALETTE["ma"], marker="D", s=48, label="Ma 2025")
        ax.scatter([0.0], [zhou_y0], color=PALETTE["zhou"], marker="s", s=54, label="Zhou 2025")
        ax.set_xlim(-1.5, 41.5)
        ax.set_xlabel("Yaw [deg]")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.3, lw=0.5)
    axes[0].legend(frameon=False, ncol=3, loc="upper center")
    fig.suptitle("Corotational mean values vs literature", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_2b_literature_comparison_corotational")


def plot_ma_spanwise_comparison() -> None:
    ma = pd.read_csv(MA_FIG16)
    r_targets = ma["r_R"].to_numpy(dtype=float)
    fig, axes = plt.subplots(2, 4, figsize=(13.4, 6.0), sharex=True)
    for column_index, yaw in enumerate([10, 20, 30, 40]):
        profile = load_spanwise_vtu_profile(COROT_ROOT / f"yaw_{yaw}" / "corotational", r_targets)
        ax_flap = axes[0, column_index]
        ax_edge = axes[1, column_index]

        ax_flap.plot(
            r_targets, ma[f"flap_{yaw}deg_m"], color=PALETTE["ma"], lw=1.1, label="Ma 2025"
        )
        ax_flap.plot(
            profile["r_R"], profile["flap"], color=PALETTE["corot"], lw=1.2, label="Corrotacional"
        )
        ax_flap.set_title(f"Yaw={yaw} deg")
        ax_flap.grid(alpha=0.3, lw=0.5)

        ax_edge.plot(
            r_targets, ma[f"edge_{yaw}deg_m"], color=PALETTE["ma"], lw=1.1, label="Ma 2025"
        )
        ax_edge.plot(
            profile["r_R"], profile["edge"], color=PALETTE["corot"], lw=1.2, label="Corrotacional"
        )
        ax_edge.grid(alpha=0.3, lw=0.5)

    axes[0, 0].set_ylabel("Flapwise [m]")
    axes[1, 0].set_ylabel("Edgewise [m]")
    for ax in axes[1, :]:
        ax.set_xlabel("r/R [-]")
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Spanwise deflection profiles - corotational vs Ma 2025", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_2c_ma_fig16_digitization_corotational")


def plot_torque_signal_analysis() -> None:
    csv_path = COROT_ROOT / "yaw_0" / "fluid" / "bem_report.csv"
    full = load_full(csv_path)
    steady = full[full["Time [s]"] >= T_START].reset_index(drop=True)
    q_center = steady["Torque [N.m]"] - float(steady["Torque [N.m]"].mean())
    p_center = steady["Power [W]"] - float(steady["Power [W]"].mean())
    y_center = steady["Tip Disp Y [m]"] - float(steady["Tip Disp Y [m]"].mean())
    slope, intercept = np.polyfit(
        y_center.to_numpy(dtype=float), q_center.to_numpy(dtype=float) / 1.0e6, 1
    )
    period_s = 60.0 / OMEGA_RPM
    recent = steady[
        steady["Time [s]"] >= float(steady["Time [s]"].max() - 3.0 * period_s)
    ].reset_index(drop=True)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.0))
    ax = axes[0, 0]
    ax.plot(
        full["Time [s]"],
        full["Torque [N.m]"] / 1.0e6,
        color=PALETTE["corot"],
        lw=1.1,
        label="Corrotacional",
    )
    ax.axhline(RIGID_BEM_TORQUE_MNM, color=PALETTE["ref"], lw=1.0, ls=":", label="Rotor rigido BEM")
    ax.axhline(
        float(steady["Torque [N.m]"].mean()) / 1.0e6,
        color=PALETTE["corot"],
        lw=1.0,
        ls="--",
        label="Media t>=40 s",
    )
    ax.axvline(T_START, color=PALETTE["neutral"], lw=0.9, ls="-.")
    ax.set_title("(a) Torque history and rigid reference")
    ax.set_ylabel("Torque [MNm]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(
        recent["Time [s]"],
        zscore(recent["Torque [N.m]"].to_numpy(dtype=float)),
        color=PALETTE["corot"],
        lw=1.1,
        label="Torque",
    )
    ax.plot(
        recent["Time [s]"],
        zscore(recent["Tip Disp Y [m]"].to_numpy(dtype=float)),
        color=PALETTE["flapwise"],
        lw=1.1,
        label="Flapwise",
    )
    ax.set_title("(b) Steady-state normalized torque and flapwise")
    ax.set_ylabel("z-score [-]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(
        recent["Time [s]"],
        zscore(recent["Power [W]"].to_numpy(dtype=float)),
        color=PALETTE["resultante"],
        lw=1.1,
        label="Power",
    )
    ax.plot(
        recent["Time [s]"],
        zscore(recent["Tip Disp Y [m]"].to_numpy(dtype=float)),
        color=PALETTE["flapwise"],
        lw=1.1,
        label="Flapwise",
    )
    ax.set_title("(c) Steady-state normalized power and flapwise")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("z-score [-]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    ax = axes[1, 1]
    x = y_center.to_numpy(dtype=float)
    y = q_center.to_numpy(dtype=float) / 1.0e6
    ax.scatter(x, y, s=10, alpha=0.20, color=PALETTE["corot"], edgecolors="none")
    x_fit = np.linspace(float(np.min(x)), float(np.max(x)), 200)
    ax.plot(
        x_fit,
        slope * x_fit + intercept,
        color=PALETTE["ref"],
        lw=1.2,
        label=f"dQ/dY={slope:.3f} MNm/m",
    )
    ax.set_title("(d) Centered torque vs centered flapwise")
    ax.set_xlabel("Y - Y_mean [m]")
    ax.set_ylabel("Q - Q_mean [MNm]")
    ax.grid(alpha=0.3, lw=0.5)
    ax.legend(frameon=False)

    deficit_pct = (
        100.0
        * ((float(steady["Torque [N.m]"].mean()) / 1.0e6) - RIGID_BEM_TORQUE_MNM)
        / RIGID_BEM_TORQUE_MNM
    )
    fig.suptitle(
        f"Torque-signal analysis for the corotational rotor (mean deficit {deficit_pct:.2f}%)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7_torque_signal_analysis_corotational")


def plot_component_influence() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "yaw_component_influence.csv")
    df = df[(df["solver"] == "corot") & (df["yaw_deg"] == 0) & (df["status"] == "ok")].copy()
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 6.4))
    configs = [
        ("torque", "pearson", "r_Q [-]", "(a) Torque correlation"),
        ("torque", "slope_scaled", "dQ/dU [MNm/m]", "(b) Torque sensitivity"),
        ("power", "pearson", "r_P [-]", "(c) Power correlation"),
        ("power", "slope_scaled", "dP/dU [MW/m]", "(d) Power sensitivity"),
    ]
    x = np.arange(len(COMPONENT_ORDER))
    for ax, (target, metric, ylabel, title) in zip(axes.ravel(), configs):
        sub = df[df["target"] == target].set_index("component").loc[COMPONENT_ORDER].reset_index()
        values = sub[metric].to_numpy(dtype=float)
        bars = ax.bar(
            x,
            values,
            color=[component_color(component) for component in COMPONENT_ORDER],
            width=0.65,
        )
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                value,
                f"{value:.3f}",
                ha="center",
                va="bottom",
                fontsize=7,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([component_label(component) for component in COMPONENT_ORDER])
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.3, axis="y", lw=0.5)
    fig.suptitle("Component-wise influence at yaw=0 - corotational solver", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7b_component_influence_corotational")


def plot_yaw_component_consistency() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "yaw_component_influence.csv")
    df = df[(df["solver"] == "corot") & (df["target"] == "torque") & (df["status"] == "ok")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), sharex=True)
    for component in COMPONENT_ORDER:
        sub = df[df["component"] == component].sort_values("yaw_deg")
        axes[0].plot(
            sub["yaw_deg"],
            sub["pearson"],
            color=component_color(component),
            marker="o",
            lw=1.2,
            label=component_label(component),
        )
        axes[1].plot(
            sub["yaw_deg"],
            sub["slope_scaled"],
            color=component_color(component),
            marker="o",
            lw=1.2,
            label=component_label(component),
        )
    axes[0].set_title("(a) Correlation with torque")
    axes[0].set_ylabel("Pearson r [-]")
    axes[1].set_title("(b) Sensitivity with torque")
    axes[1].set_ylabel("dQ/dU [MNm/m]")
    for ax in axes:
        ax.set_xlabel("Yaw [deg]")
        ax.set_xlim(-1.5, 41.5)
        ax.grid(alpha=0.3, lw=0.5)
    axes[0].legend(frameon=False)
    fig.suptitle("Yaw dependence of the torque-deformation coupling", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7c_yaw_component_consistency_corotational")


def plot_yaw_modal_closure() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "yaw_component_influence.csv")
    df = df[(df["solver"] == "corot") & (df["target"] == "torque") & (df["status"] == "ok")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2), sharex=True)
    for component in COMPONENT_ORDER:
        sub = df[df["component"] == component].sort_values("yaw_deg")
        axes[0].plot(
            sub["yaw_deg"],
            sub["coupling_freq"],
            color=component_color(component),
            marker="o",
            lw=1.2,
            label=component_label(component),
        )
        axes[1].plot(
            sub["yaw_deg"],
            sub["phase_deg"],
            color=component_color(component),
            marker="o",
            lw=1.2,
            label=component_label(component),
        )
    axes[0].axhline(F_1P, color=PALETTE["1P"], lw=1.0, ls=":", label="1P")
    axes[0].axhline(F_FLAP, color=PALETTE["flap"], lw=1.0, ls="--", label="f1 flap")
    axes[0].set_title("(a) Dominant coupling frequency")
    axes[0].set_ylabel("f_c [Hz]")
    axes[1].set_title("(b) Torque phase relative to deformation")
    axes[1].set_ylabel("Phase [deg]")
    for ax in axes:
        ax.set_xlabel("Yaw [deg]")
        ax.set_xlim(-1.5, 41.5)
        ax.grid(alpha=0.3, lw=0.5)
    axes[0].legend(frameon=False)
    fig.suptitle("Modal closure of the torque-deformation coupling", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7d_yaw_modal_closure_corotational")


def plot_power_chain_closure() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "yaw_component_influence.csv")
    df = df[(df["solver"] == "corot") & (df["status"] == "ok")].copy()
    torque = df[df["target"] == "torque"].copy()
    power = df[df["target"] == "power"].copy()
    merged = torque.merge(
        power, on=["solver", "yaw_deg", "component", "status"], suffixes=("_torque", "_power")
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.0))
    panels = [
        (
            "slope_scaled_torque",
            "slope_scaled_power",
            OMEGA_RAD_S,
            "dQ/dU [MNm/m]",
            "dP/dU [MW/m]",
            "(a) Sensitivity closure",
        ),
        ("pearson_torque", "pearson_power", 1.0, "r_Q [-]", "r_P [-]", "(b) Correlation closure"),
        (
            "coupling_freq_torque",
            "coupling_freq_power",
            1.0,
            "f_Q [Hz]",
            "f_P [Hz]",
            "(c) Frequency closure",
        ),
        (
            "phase_deg_torque",
            "phase_deg_power",
            1.0,
            "phi_Q [deg]",
            "phi_P [deg]",
            "(d) Phase closure",
        ),
    ]
    for ax, (xcol, ycol, factor, xlabel, ylabel, title) in zip(axes.ravel(), panels):
        for component in COMPONENT_ORDER:
            sub = merged[merged["component"] == component]
            x = sub[xcol].to_numpy(dtype=float)
            y = sub[ycol].to_numpy(dtype=float)
            if xcol == "slope_scaled_torque":
                x = x * factor
            ax.scatter(
                x, y, color=component_color(component), s=42, label=component_label(component)
            )
            for x_value, y_value, yaw in zip(x, y, sub["yaw_deg"].to_numpy(dtype=int)):
                ax.text(x_value, y_value, str(yaw), fontsize=6, ha="left", va="bottom")
        values = np.concatenate([
            ax.collections[0].get_offsets()[:, 0],
            ax.collections[0].get_offsets()[:, 1],
        ])
        finite = values[np.isfinite(values)]
        if len(finite) > 0:
            low = float(np.min(finite))
            high = float(np.max(finite))
            line = np.linspace(low, high, 100)
            ax.plot(line, line, color=PALETTE["neutral"], lw=1.0, ls="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.3, lw=0.5)
    axes[0, 0].legend(frameon=False)
    fig.suptitle("Numerical closure of deformation -> torque -> power", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7e_power_chain_closure_corotational")


def plot_azimuthal_folded_signals() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "azimuthal_folded_signals.csv")
    df = df[(df["solver"] == "corot") & (df["yaw_deg"] == 0)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharex=True)
    grouped = [
        (axes[0], ["thrust", "torque", "power"], "(a) Rotor-level loads"),
        (axes[1], ["edgewise", "flapwise", "resultante"], "(b) Tip-deformation channels"),
    ]
    for ax, signals, title in grouped:
        for signal in signals:
            sub = df[df["signal"] == signal].sort_values("azimuth_deg")
            ax.plot(
                sub["azimuth_deg"],
                sub["mean_norm"],
                color=PALETTE.get(signal, PALETTE["ref"]),
                lw=1.2,
                label=SIGNAL_LABELS[signal],
            )
        ax.set_title(title)
        ax.set_xlabel("Azimuth aligned with torque 1P peak [deg]")
        ax.set_ylabel("Normalized centered signal [-]")
        ax.set_xlim(0, 360)
        ax.grid(alpha=0.3, lw=0.5)
        ax.legend(frameon=False)
    fig.suptitle("Azimuthal folding around the torque 1P peak - corotational solver", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7f_azimuthal_folded_signals_corotational")


def plot_harmonic_budget() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "azimuthal_harmonics.csv")
    df = df[(df["solver"] == "corot") & (df["status"] == "ok")].copy()
    signals = ["thrust", "torque", "edgewise", "flapwise"]
    harmonics = ["1P", "2P", "3P", "Residual"]
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 6.8), sharex=True, sharey=True)
    for ax, signal in zip(axes.ravel(), signals):
        pivot = (
            df[df["signal"] == signal]
            .pivot_table(
                index="yaw_deg", columns="harmonic", values="variance_fraction", aggfunc="first"
            )
            .reindex(index=YAWS)
            .fillna(0.0)
        )
        bottom = np.zeros(len(pivot), dtype=float)
        for harmonic in harmonics:
            values = pivot[harmonic].to_numpy(dtype=float)
            ax.bar(
                pivot.index.to_numpy(dtype=float),
                values,
                bottom=bottom,
                width=7.0,
                color=PALETTE.get(
                    {"1P": "1P", "2P": "flapwise", "3P": "resultante", "Residual": "residual"}[
                        harmonic
                    ],
                    PALETTE["neutral"],
                ),
                label=harmonic,
            )
            bottom += values
        ax.set_title(SIGNAL_LABELS[signal])
        ax.set_ylim(0.0, 1.0)
        ax.grid(alpha=0.25, axis="y", lw=0.5)
    for ax in axes[1, :]:
        ax.set_xlabel("Yaw [deg]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Variance fraction [-]")
    handles = [
        Line2D([0], [0], color=PALETTE["1P"], lw=6, label="1P"),
        Line2D([0], [0], color=PALETTE["flapwise"], lw=6, label="2P"),
        Line2D([0], [0], color=PALETTE["resultante"], lw=6, label="3P"),
        Line2D([0], [0], color=PALETTE["residual"], lw=6, label="Residual"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=4)
    fig.suptitle("Harmonic budget across the corotational yaw sweep", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save_figure(fig, "fig_5_4_7g_harmonic_budget_corotational")


def plot_hysteresis_work_loops() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "aeroelastic_work_loops.csv")
    df = df[(df["solver"] == "corot") & (df["yaw_deg"] == 0)].sort_values("azimuth_deg")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    axes[0].plot(df["flapwise_prime_m"], df["torque_prime_MNm"], color=PALETTE["corot"], lw=1.2)
    axes[0].scatter(
        df["flapwise_prime_m"].iloc[0],
        df["torque_prime_MNm"].iloc[0],
        color=PALETTE["ref"],
        s=28,
        zorder=4,
    )
    axes[0].set_title("(a) Torque prime vs flapwise prime")
    axes[0].set_xlabel("Flapwise prime [m]")
    axes[0].set_ylabel("Torque prime [MNm]")
    axes[0].grid(alpha=0.3, lw=0.5)

    axes[1].plot(df["flapwise_prime_m"], df["thrust_prime_MN"], color=PALETTE["ma"], lw=1.2)
    axes[1].scatter(
        df["flapwise_prime_m"].iloc[0],
        df["thrust_prime_MN"].iloc[0],
        color=PALETTE["ref"],
        s=28,
        zorder=4,
    )
    axes[1].set_title("(b) Thrust prime vs flapwise prime")
    axes[1].set_xlabel("Flapwise prime [m]")
    axes[1].set_ylabel("Thrust prime [MN]")
    axes[1].grid(alpha=0.3, lw=0.5)

    fig.suptitle("Cycle-averaged hysteresis proxies at yaw=0 - corotational solver", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7h_hysteresis_work_loops_corotational")


def plot_cycle_work_proxy() -> None:
    df = save_df_plot_data(pd.DataFrame(), GENERATED_DIR / "aeroelastic_work_proxy.csv")
    df = df[(df["solver"] == "corot") & (df["status"] == "ok")].sort_values("yaw_deg")
    metrics = [
        ("qdy_proxy_mean_MJm", "Integral Q' dY [MJ m]", "(a) Tangential loop proxy"),
        ("tdy_work_mean_MJ", "Integral T' dY [MJ]", "(b) Axial loop proxy"),
        ("qv_proxy_mean_MWm", "<Q' Ydot> [MW m]", "(c) Tangential quadrature proxy"),
        ("tv_power_mean_MW", "<T' Ydot> [MW]", "(d) Axial quadrature proxy"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.8), sharex=True)
    for ax, (column, ylabel, title) in zip(axes.ravel(), metrics):
        ax.axhline(0.0, color=PALETTE["neutral"], lw=0.9, ls="--")
        ax.plot(df["yaw_deg"], df[column], color=PALETTE["corot"], marker="o", lw=1.2)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3, lw=0.5)
    for ax in axes[1, :]:
        ax.set_xlabel("Yaw [deg]")
        ax.set_xlim(-1.5, 41.5)
    fig.suptitle("Cycle-to-cycle aeroelastic work proxies - corotational solver", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_4_7i_cycle_work_proxy_corotational")


def plot_fft_tip_deflection() -> None:
    tip = load_tip_vtu_series(COROT_ROOT / "yaw_0" / "corotational")
    dt = float(tip["Time [s]"].diff().dropna().median())
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0))
    for ax, column, title in zip(
        axes,
        ["Tip Disp Y [m]", "Tip Disp X [m]"],
        ["(a) Flapwise tip FFT", "(b) Edgewise tip FFT"],
    ):
        freq, amp = onesided_fft(tip[column], dt)
        mask = freq <= 0.75
        ax.plot(freq[mask], amp[mask], color=PALETTE["corot"], lw=1.2)
        for marker, label, color in [
            (F_1P, "1P", PALETTE["1P"]),
            (F_3P, "3P", PALETTE["3P"]),
            (F_FLAP, "f1 flap", PALETTE["flap"]),
            (F_EDGE, "f1 edge", PALETTE["edge"]),
        ]:
            ax.axvline(marker, color=color, lw=1.0, ls=":")
            ax.text(
                marker + 0.004,
                0.92 * ax.get_ylim()[1],
                label,
                color=color,
                fontsize=6.5,
                rotation=90,
                va="top",
            )
        ax.set_title(title)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("Amplitude [m]")
        ax.set_xlim(0.0, 0.75)
        ax.grid(alpha=0.3, lw=0.5)
    fig.suptitle(
        f"Tip-node FFT spectra for the corotational solver (node={int(tip['tip_node'].iloc[0])})",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_figure(fig, "fig_5_5_1_fft_cargas_corotational")


def main() -> None:
    from plot_torque_deformation_frequency_overlay import main as plot_frequency_overlay

    summary = summarize_corotational_campaign()
    print("Generating corotational-only standalone figures ...")
    plot_timeseries_deflection(summary)
    plot_yaw_sweep(summary)
    plot_literature_comparison(summary)
    plot_ma_spanwise_comparison()
    plot_torque_signal_analysis()
    plot_component_influence()
    plot_yaw_component_consistency()
    plot_yaw_modal_closure()
    plot_power_chain_closure()
    plot_azimuthal_folded_signals()
    plot_harmonic_budget()
    plot_hysteresis_work_loops()
    plot_cycle_work_proxy()
    plot_fft_tip_deflection()
    plot_frequency_overlay()
    print("Done.")


if __name__ == "__main__":
    main()
