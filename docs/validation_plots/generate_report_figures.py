"""
Generate all report figures for validation_results_v01_v02_v04.md

Outputs (docs/figures/):
  fig_5_5_1_fft_cargas.png            — §5.5.1  FFT deflexión flapwise/edgewise (estilo Zhou)
  fig_v02_modal_comparison.png        — V-02     Comparación modal multi-método (Zhou Table 3 + AeroElast)
  fig_5_2_yaw_sweep.png               — §5.2     Barrido de yaw corrotacional vs inercial (Cp, P, flap)
  fig_5_1_timeseries_deflection.png   — §5.1     Series temporales tip deflexión flapwise/edgewise
    fig_5_2b_literature_comparison.png  — §5.4     Comparación literatura: flap/edge tip + thrust vs yaw
    fig_5_2c_ma_fig16_digitization.png  — §5.4     Comparación spanwise AeroElast vs Ma et al. 2025 Figure 16
"""

import base64
import mmap
import pathlib
import re
import zlib

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from numpy import trapezoid as _trapz
from scipy.signal import windows

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
OUT_DIR = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA_DIR = pathlib.Path(__file__).resolve().parents[1] / "validation_data"
GENERATED_DIR = DATA_DIR / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
GENERATED_DIR.mkdir(parents=True, exist_ok=True)

COROT_ROOT = BASE / "frontiersin_results_corotational_100s"
INERT_ROOT = BASE / "frontiersin_results_inertial"
METRICS_ALL = GENERATED_DIR / "campaign_metrics_summary.csv"
METRICS_INERT = GENERATED_DIR / "inertial_metrics_summary.csv"
COROT_TIP = COROT_ROOT / "yaw_0/fluid/bem_report.csv"
INERT_TIP = INERT_ROOT / "yaw_0/fluid/bem_report.csv"
MA_FIG15 = DATA_DIR / "ma_2025_fig15_yaw_metrics.csv"
MA_FIG16 = DATA_DIR / "ma_2025_fig16_spanwise_deflections.csv"

T_START = 40.0
T_END = 100.0
DT = 0.01
SPANWISE_PROFILE_STRIDE = 50

F_1P = 7.518 / 60  # 0.12530 Hz
F_3P = 3 * F_1P  # 0.37590 Hz
F_FLAP = 0.5537  # Hz  V-02
F_EDGE = 0.6290  # Hz  V-02

PALETTE = {
    "corrot": "#c0392b",  # dark red
    "inert": "#2980b9",  # steel blue
    "1P": "#e74c3c",
    "3P": "#3498db",
    "flap": "#27ae60",
    "edge": "#f39c12",
    "zhou": "#8e44ad",
    "ref": "#2c3e50",
    "ma2025": "#16a085",  # teal — Ma et al. 2025 LL-FVW+GEBT
}


# ============================================================
# Helper: FFT with Hann window, one-sided, amplitude-preserving
# ============================================================


def onesided_fft(series: pd.Series, dt: float):
    n = len(series)
    win = windows.hann(n)
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    x = (values - np.nanmean(values)) * (win / win.mean())
    sp = np.fft.rfft(x)
    freq = np.fft.rfftfreq(n, d=dt)
    amp = np.abs(sp) * 2.0 / n
    amp[0] = 0.0
    return freq, amp


def load_trim(path, t0=T_START, t1=T_END):
    df = pd.read_csv(path)
    return df[(df["Time [s]"] >= t0) & (df["Time [s]"] <= t1)].reset_index(drop=True)


def steady_window(df: pd.DataFrame, t0=T_START, t1=T_END) -> pd.DataFrame:
    return df[(df["Time [s]"] >= t0) & (df["Time [s]"] <= t1)].reset_index(drop=True)


def _time_value(path: pathlib.Path) -> float:
    return float(path.parent.name)


def _field_vtu_paths(case_dir: pathlib.Path) -> list[pathlib.Path]:
    paths = [path for path in case_dir.glob("*/fields.vtu") if path.parent.name != "fluid"]
    return sorted(paths, key=_time_value)


def _decode_vtu_payload(encoded: bytes, dtype: np.dtype) -> np.ndarray:
    payload = base64.b64decode(encoded)
    if len(payload) >= 12:
        n_blocks, block_size, _ = np.frombuffer(payload, dtype="<u4", count=3)
        header_size = 4 * (3 + int(n_blocks))
        if 0 < n_blocks < 10000 and block_size > 0 and header_size <= len(payload):
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
        raise TypeError(
            f"Unsupported VTU array type for {name!r}: {vtk_type.group(1) if vtk_type else 'missing'}"
        )
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


def load_tip_vtu_series(case_dir: pathlib.Path, t0=T_START, t1=T_END) -> pd.DataFrame:
    paths = _field_vtu_paths(case_dir)
    if not paths:
        raise FileNotFoundError(f"No fields.vtu files found under {case_dir}")
    paths = [path for path in paths if t0 <= _time_value(path) <= t1]
    tip_node, z_tip, n_candidates = _select_tip_node_from_vtu(paths[0])
    rows = []
    for path in paths:
        time = _time_value(path)
        displacement = _vtu_data_array(path, "U")[tip_node]
        rows.append({
            "Time [s]": time,
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
    if span <= 0.0:
        raise ValueError(f"Invalid span reconstructed from {path}")
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


def _spanwise_time_upper_bound(
    case_dir: pathlib.Path, t0=T_START, t1=T_END, flap_limit=20.0
) -> float:
    report_path = case_dir / "structural_report.csv"
    if not report_path.exists():
        return t1

    report = pd.read_csv(report_path)
    if "Time [s]" not in report.columns or "Max Disp Y [m]" not in report.columns:
        return t1

    report = report[(report["Time [s]"] >= t0) & (report["Time [s]"] <= t1)].reset_index(drop=True)
    if report.empty:
        return t1

    flap = pd.to_numeric(report["Max Disp Y [m]"], errors="coerce")
    valid = np.asarray(np.isfinite(flap) & (np.abs(flap) <= flap_limit), dtype=bool)
    if valid.all():
        return t1

    invalid_idx = np.flatnonzero(~valid)
    if len(invalid_idx) == 0 or invalid_idx[0] == 0:
        return t1
    return float(report["Time [s]"].iloc[int(invalid_idx[0] - 1)])


def load_spanwise_vtu_profile(
    case_dir: pathlib.Path,
    r_targets: np.ndarray,
    t0=T_START,
    t1=T_END,
    stride=SPANWISE_PROFILE_STRIDE,
) -> dict[str, np.ndarray | float | int]:
    t_end = _spanwise_time_upper_bound(case_dir, t0, t1)
    paths = [
        path
        for path in _field_vtu_paths(case_dir)
        if _time_value(path) >= t0 and _time_value(path) <= t_end + 1e-12
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

    n_frames = len(sampled_paths)
    flap_mean = flap_sum / n_frames
    edge_mean = edge_sum / n_frames

    edge_sign = 1.0
    bem_path = case_dir.parent / "fluid/bem_report.csv"
    if bem_path.exists():
        bem_trim = load_trim(bem_path, t0)
        if not bem_trim.empty:
            tip_edge_bem = float(bem_trim["Tip Disp X [m]"].mean())
            tip_edge_vtu = float(edge_mean[-1])
            if np.isfinite(tip_edge_bem) and np.isfinite(tip_edge_vtu):
                bem_sign = float(np.sign(tip_edge_bem))
                vtu_sign = float(np.sign(tip_edge_vtu))
                if bem_sign != 0.0 and vtu_sign != 0.0:
                    edge_sign = bem_sign * vtu_sign
                    edge_mean = edge_mean * edge_sign

    return {
        "r_R": np.asarray(r_targets, dtype=float),
        "flap": flap_mean,
        "edge": edge_mean,
        "edge_sign": float(edge_sign),
        "n_frames": int(n_frames),
        "station_nodes": station_nodes,
        "t_end": float(t_end),
    }


def vline_linear(ax, x, label, color, yrel=0.88):
    ax.axvline(x, color=color, lw=1.0, ls=":", alpha=0.9)
    ylims = ax.get_ylim()
    ypos = ylims[0] + yrel * (ylims[1] - ylims[0])
    ax.text(x + 0.004, ypos, label, color=color, fontsize=6.5, va="center", rotation=90)


# ============================================================
# Figure 1 — FFT deflexión flapwise/edgewise  §5.5.1
#   Estilo Zhou (2024): escala lineal, 2 paneles, 0–0.5 Hz
# ============================================================

print("Generating fig_5_5_1_fft_cargas.png …")

corot_s = load_tip_vtu_series(COROT_ROOT / "yaw_0/corotational")
inert_s = load_tip_vtu_series(INERT_ROOT / "yaw_0/inertial")
print(
    "    VTU tip nodes: "
    f"corot={int(corot_s['tip_node'].iloc[0])} "
    f"(Z0={corot_s['tip_z0'].iloc[0]:.3f}, candidates={int(corot_s['tip_z0_candidates'].iloc[0])}), "
    f"inert={int(inert_s['tip_node'].iloc[0])} "
    f"(Z0={inert_s['tip_z0'].iloc[0]:.3f}, candidates={int(inert_s['tip_z0_candidates'].iloc[0])})"
)

# Compute dt from time column (robust to varying step)
dt_c = float(corot_s["Time [s]"].diff().dropna().median())
dt_i = float(inert_s["Time [s]"].diff().dropna().median())

PANELS_FFT = [
    # (panel_title, col, ylabel)
    ("(a) Flapwise tip deflection", "Tip Disp Y [m]", "Amplitude [m]"),
    ("(b) Edgewise tip deflection", "Tip Disp X [m]", "Amplitude [m]"),
]

fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=False)

for ax, (title, col, ylabel) in zip(axes, PANELS_FFT):
    fc, ac = onesided_fft(corot_s[col], dt_c)
    fi, ai = onesided_fft(inert_s[col], dt_i)

    # Restrict to 0–0.75 Hz (so f₁_flap and f₁_edge are visible)
    mask_c = fc <= 0.75
    mask_i = fi <= 0.75

    ax.plot(fc[mask_c], ac[mask_c], lw=1.3, color=PALETTE["corrot"], label="Corotational")
    ax.plot(
        fi[mask_i],
        ai[mask_i],
        lw=1.3,
        color=PALETTE["inert"],
        label="Inertial",
        ls="--",
        alpha=0.85,
    )

    ax.set_xlim(0, 0.75)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.grid(alpha=0.3, lw=0.5)

    # Vertical frequency markers — placed AFTER ylim is auto-set
    ax.autoscale(enable=True, axis="y", tight=False)
    for fx, lbl, col_ in [
        (F_1P, "1P", PALETTE["1P"]),
        (F_3P, "3P", PALETTE["3P"]),
        (F_FLAP, "$f_1^{fl}$", PALETTE["flap"]),
        (F_EDGE, "$f_1^{ed}$", PALETTE["edge"]),
    ]:
        vline_linear(ax, fx, lbl, col_)

    ax.legend(frameon=False, loc="upper right", fontsize=7.5)

# Common frequency legend at bottom
freq_handles = [
    Line2D([0], [0], color=PALETTE["1P"], ls=":", lw=1.2, label=f"1P = {F_1P:.4f} Hz"),
    Line2D([0], [0], color=PALETTE["3P"], ls=":", lw=1.2, label=f"3P = {F_3P:.4f} Hz"),
    Line2D(
        [0], [0], color=PALETTE["flap"], ls=":", lw=1.2, label=f"$f_1^{{fl}}$ = {F_FLAP:.4f} Hz"
    ),
    Line2D(
        [0], [0], color=PALETTE["edge"], ls=":", lw=1.2, label=f"$f_1^{{ed}}$ = {F_EDGE:.4f} Hz"
    ),
]
fig.legend(
    handles=freq_handles,
    loc="lower center",
    ncol=4,
    frameon=False,
    bbox_to_anchor=(0.5, -0.04),
    fontsize=7.5,
)

fig.suptitle(
    "Frequency spectra (Hann window) — IEA 15 MW, yaw=0°, rated condition",
    fontsize=9.5,
)
fig.tight_layout(rect=(0, 0.06, 1, 1))
out1 = OUT_DIR / "fig_5_5_1_fft_cargas.png"
fig.savefig(out1, bbox_inches="tight")
fig.savefig(out1.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"  → {out1}")


# ============================================================
# Figure 2 — Comparación modal multi-método  V-02
# ============================================================

print("Generating fig_v02_modal_comparison.png …")

# Data from Zhou Energy 2025 Table 3 + AeroElast V-02 (Hz)
# NaN = not reported in that source
MODAL_DATA = {
    "HWAC2\n[74]": [0.512, 0.692, 1.509, 2.120, 4.314],
    "HWAC2\n[26]": [0.521, 0.619, 1.559, 1.933, 4.475],
    "3D-FEM\n[55]": [0.465, 0.547, 1.452, 1.671, 4.295],
    "3D-FEM\n[75]": [0.555, 0.642, 1.598, 1.925, 3.911],
    "Tech.Rep.\n[56]": [0.555, 0.642, np.nan, np.nan, np.nan],
    "Zhou\nGEBT": [0.523, 0.665, 1.475, 2.124, 4.072],
    "AeroElast\nShell": [0.5537, 0.6290, 1.6946, np.nan, np.nan],
}

MODES = ["1er Flapwise", "1er Edgewise", "2do Flapwise", "2do Edgewise", "1er Torsional"]
methods = list(MODAL_DATA.keys())
n_methods = len(methods)
n_modes = len(MODES)

x = np.arange(n_modes)
w = 0.11
gap = 0.01

# Colours: grey for literature, highlight for Zhou and AeroElast
col_list = [
    "#95a5a6",
    "#7f8c8d",
    "#bdc3c7",
    "#aab7b8",
    "#808b96",
    PALETTE["zhou"],
    PALETTE["corrot"],
]

fig, ax = plt.subplots(figsize=(11, 4.5))

for j, (method, vals) in enumerate(MODAL_DATA.items()):
    offsets = x + j * (w + gap) - (n_methods - 1) * (w + gap) / 2
    color = col_list[j]
    lw = 1.8 if method in ("Zhou\nGEBT", "AeroElast\nShell") else 0.6
    ec = color if method in ("Zhou\nGEBT", "AeroElast\nShell") else "#555"
    for xi, v in zip(offsets, vals):
        if not np.isnan(v):
            ax.bar(xi, v, width=w, color=color, edgecolor=ec, linewidth=lw, zorder=3)
        else:
            ax.bar(xi, 0, width=w, color="none", zorder=3)

ax.set_xticks(x)
ax.set_xticklabels(MODES, fontsize=8.5)
ax.set_ylabel("Frecuencia natural [Hz]")
ax.set_title(
    "Frecuencias naturales de la pala IEA 15 MW — comparación multi-método\n"
    "(Fuente: Zhou et al. Energy 2025, Tabla 3)",
    fontsize=9.5,
)
ax.grid(axis="y", alpha=0.3, zorder=0)
ax.set_ylim(0, 5.2)

# Legend
handles = [
    mpatches.Patch(color=col_list[j], label=m.replace("\n", " ")) for j, m in enumerate(methods)
]
ax.legend(
    handles=handles, ncol=7, loc="upper left", frameon=False, fontsize=7.5, bbox_to_anchor=(0, 1.0)
)

fig.tight_layout()
out2 = OUT_DIR / "fig_v02_modal_comparison.png"
fig.savefig(out2, bbox_inches="tight")
fig.savefig(out2.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"  → {out2}")


# ============================================================
# Figure 3 — Barrido de yaw  §5.2
# ============================================================

print("Generating fig_5_2_yaw_sweep.png …")

N_BLADES = 3
RHO = 1.225
_R = 120.0
_A = np.pi * _R**2
_V = 10.659


def _bem_rotor_ts(fluid_base: pathlib.Path, omega_series: pd.Series) -> pd.DataFrame:
    """Integrate per-timestep BEM sectional CSV → rotor Cp/P/Ct time series."""
    times, Cp_v, P_v, Ct_v = [], [], [], []
    dirs = sorted([float(d.name) for d in fluid_base.iterdir() if d.is_dir()])
    for t in dirs:
        if t < T_START or t > T_END:
            continue
        sec_path = fluid_base / f"{t:.2g}" / "bem_sectional.csv"
        if not sec_path.exists():
            sec_path = fluid_base / str(t) / "bem_sectional.csv"
        if not sec_path.exists():
            continue
        sec = pd.read_csv(sec_path)
        thrust = N_BLADES * _trapz(sec["Np[N/m]"], sec["r[m]"])
        torque = N_BLADES * _trapz(sec["Tp[N/m]"] * sec["r[m]"], sec["r[m]"])
        idx = omega_series.index.get_indexer(pd.Index([t]), method="nearest")[0]
        omega = omega_series.iloc[idx]
        power = torque * omega
        times.append(t)
        Cp_v.append(power / (0.5 * RHO * _A * _V**3))
        P_v.append(power / 1e6)
        Ct_v.append(thrust / (0.5 * RHO * _A * _V**2))
    return pd.DataFrame({"t": times, "Cp": Cp_v, "P_MW": P_v, "Ct": Ct_v})


def summarize_corotational_campaign(
    root: pathlib.Path, yaw_angles=(0, 10, 20, 30, 40)
) -> pd.DataFrame:
    """Summarize the updated corotational campaign from signed BEM report files."""
    return summarize_fluid_campaign(root, "corotational", yaw_angles)


def summarize_inertial_campaign(root: pathlib.Path, yaw_angles=(0, 10, 20, 30, 40)) -> pd.DataFrame:
    """Summarize the completed inertial campaign from signed BEM report files."""
    return summarize_fluid_campaign(root, "inertial", yaw_angles)


def summarize_fluid_campaign(
    root: pathlib.Path, solver_name: str, yaw_angles=(0, 10, 20, 30, 40)
) -> pd.DataFrame:
    """Summarize one campaign from raw fluid/bem_report.csv files."""
    rows = []
    for yaw in yaw_angles:
        report_csv = root / f"yaw_{yaw}" / "fluid" / "bem_report.csv"
        if not report_csv.exists():
            rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                "status": "missing_csv",
                "source": str(report_csv),
                "diverged": True,
            })
            continue
        df = pd.read_csv(report_csv)
        df = steady_window(df)
        if df.empty:
            rows.append({
                "solver": solver_name,
                "yaw_deg": yaw,
                "status": "empty_window",
                "source": str(report_csv),
                "diverged": True,
            })
            continue
        thrust = pd.to_numeric(df["Thrust [N]"], errors="coerce")
        torque = pd.to_numeric(df["Torque [N.m]"], errors="coerce")
        power = pd.to_numeric(df["Power [W]"], errors="coerce")
        finite = bool(
            np.isfinite(df.select_dtypes(include=[np.number]).to_numpy(dtype=float)).all()
        )
        diverged = bool(
            (not finite)
            or np.any(np.abs(thrust.to_numpy(dtype=float)) > 1.0e8)
            or np.any(np.abs(torque.to_numpy(dtype=float)) > 1.0e8)
            or np.any(np.abs(power.to_numpy(dtype=float)) > 1.0e8)
        )
        rows.append({
            "solver": solver_name,
            "yaw_deg": yaw,
            "status": "ok" if not diverged else "invalid_physical_bounds",
            "source": str(report_csv),
            "n_samples": len(df),
            "t_start_actual": float(df["Time [s]"].min()),
            "t_end_actual": float(df["Time [s]"].max()),
            "flap_mean": float(df["Tip Disp Y [m]"].mean()),
            "flap_std": float(df["Tip Disp Y [m]"].std()),
            "flap_max": float(df["Tip Disp Y [m]"].max()),
            "flap_min": float(df["Tip Disp Y [m]"].min()),
            "flap_p2p": float(df["Tip Disp Y [m]"].max() - df["Tip Disp Y [m]"].min()),
            "edge_mean": float(df["Tip Disp X [m]"].mean()),
            "edge_std": float(df["Tip Disp X [m]"].std()),
            "edge_max": float(df["Tip Disp X [m]"].max()),
            "edge_min": float(df["Tip Disp X [m]"].min()),
            "edge_p2p": float(df["Tip Disp X [m]"].max() - df["Tip Disp X [m]"].min()),
            "thrust_mean": float(df["Thrust [N]"].mean()),
            "thrust_std": float(df["Thrust [N]"].std()),
            "torque_mean": float(df["Torque [N.m]"].mean()),
            "torque_std": float(df["Torque [N.m]"].std()),
            "power_mean_mw": float(df["Power [W]"].mean() / 1e6),
            "power_std_mw": float(df["Power [W]"].std() / 1e6),
            "power_max_mw": float(df["Power [W]"].max() / 1e6),
            "power_p2p": float((df["Power [W]"].max() - df["Power [W]"].min()) / 1e6),
            "ct_mean": float(df["CT"].mean()),
            "ct_std": float(df["CT"].std()),
            "cq_mean": float(df["CQ"].mean()),
            "cp_mean": float(df["CP"].mean()),
            "cp_std": float(df["CP"].std()),
            "diverged": diverged,
        })
    return pd.DataFrame(rows).sort_values("yaw_deg").reset_index(drop=True)


def load_metrics_corot(root: pathlib.Path) -> dict:
    """Corotational: reads aerodynamic and signed tip metrics from the new campaign."""
    df = summarize_corotational_campaign(root)
    return {
        "yaw": df["yaw_deg"].values.astype(float),
        "Cp": df["cp_mean"].values,
        "P_MW": df["power_mean_mw"].values,
        "flap_m": df["flap_mean"].values,
        "flap_pp": df["flap_p2p"].values,
        "valid": ~df["diverged"].to_numpy(dtype=bool),
    }


def load_metrics_inert(base: pathlib.Path, yaw_angles=(0, 10, 20, 30, 40)) -> dict:
    """Inertial: metrics from the completed signed BEM report campaign."""
    df = summarize_inertial_campaign(base, yaw_angles)
    ok = df["status"].eq("ok")
    return {
        "yaw": df["yaw_deg"].values.astype(float),
        "Cp": df["cp_mean"].values,
        "P_MW": df["power_mean_mw"].values,
        "flap_m": df["flap_mean"].values,
        "flap_pp": df["flap_p2p"].values,
        "valid": ok.to_numpy(dtype=bool),
    }


mc = load_metrics_corot(COROT_ROOT)
print("  Loading inertial metrics from completed BEM reports …")
mi = load_metrics_inert(INERT_ROOT)

_mc_df_export = summarize_corotational_campaign(COROT_ROOT)
_mi_df_export = summarize_inertial_campaign(INERT_ROOT)
_metrics_export = pd.concat([_mc_df_export, _mi_df_export], ignore_index=True)
_metrics_export.to_csv(METRICS_ALL, index=False)
_mi_df_export.to_csv(METRICS_INERT, index=False)
print(f"  → {METRICS_ALL}")

yaw_arr_c = mc["yaw"]
yaw_arr_i = mi["yaw"]
yaw_fine = np.linspace(0, 40, 200)
cos3_fine = np.cos(np.radians(yaw_fine)) ** 3

fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))

# Panel 1: Cp + cos³ fit (corotational only — inertial Cp is anomalous)
ax = axes[0]
ax.plot(yaw_fine, mc["Cp"][0] * cos3_fine, color="#95a5a6", lw=1.5, ls="--", label="cos³(γ) fit")
ax.scatter(
    yaw_arr_c[mc["valid"]],
    mc["Cp"][mc["valid"]],
    color=PALETTE["corrot"],
    s=60,
    zorder=5,
    label="Corrotacional",
)
ax.scatter(
    yaw_arr_c[~mc["valid"]],
    mc["Cp"][~mc["valid"]],
    color=PALETTE["corrot"],
    s=60,
    zorder=5,
    marker="x",
    linewidths=2,
)
# Inertial: correct values from BEM fluid integration
ax.scatter(
    yaw_arr_i[mi["valid"]],
    mi["Cp"][mi["valid"]],
    color=PALETTE["inert"],
    s=60,
    zorder=5,
    marker="^",
    label="Inercial",
    alpha=0.85,
)
ax.scatter(
    yaw_arr_i[~mi["valid"]],
    mi["Cp"][~mi["valid"]],
    color=PALETTE["inert"],
    s=60,
    zorder=5,
    marker="x",
    linewidths=2,
    alpha=0.85,
)
ax.set_xlabel("Ángulo de yaw [°]")
ax.set_ylabel("C$_P$")
ax.set_title("Coeficiente de potencia")
ax.legend(frameon=False, fontsize=7)
ax.grid(alpha=0.3)
ax.set_xlim(-2, 43)

# Panel 2: Potencia [MW]
ax = axes[1]
ax.plot(yaw_fine, mc["P_MW"][0] * cos3_fine, color="#95a5a6", lw=1.5, ls="--", label="P₀ · cos³(γ)")
ax.scatter(
    yaw_arr_c[mc["valid"]],
    mc["P_MW"][mc["valid"]],
    color=PALETTE["corrot"],
    s=60,
    zorder=5,
    label="Corrotacional",
)
ax.scatter(
    yaw_arr_c[~mc["valid"]],
    mc["P_MW"][~mc["valid"]],
    color=PALETTE["corrot"],
    s=60,
    zorder=5,
    marker="x",
    linewidths=2,
)
ax.scatter(
    yaw_arr_i[mi["valid"]],
    mi["P_MW"][mi["valid"]],
    color=PALETTE["inert"],
    s=60,
    zorder=5,
    marker="^",
    label="Inercial",
    alpha=0.85,
)
ax.scatter(
    yaw_arr_i[~mi["valid"]],
    mi["P_MW"][~mi["valid"]],
    color=PALETTE["inert"],
    s=60,
    zorder=5,
    marker="x",
    linewidths=2,
    alpha=0.85,
)
ma_fig15 = pd.read_csv(MA_FIG15)
ma_fig16 = pd.read_csv(MA_FIG16)
ma_fig16_tip = ma_fig16.loc[np.isclose(ma_fig16["r_R"], 1.0)].iloc[0]

# Ma et al. 2025 reference points (LL-FVW+GEBT, digitized from Figure 15)
_ma_yaw = ma_fig15["yaw_angle_deg"].to_numpy(dtype=float)
_ma_P = ma_fig15["mean_power_MW"].to_numpy(dtype=float)
_ma_P0 = _ma_P[0] / np.cos(np.radians(_ma_yaw[0])) ** 1.5  # ~14.52 MW fitted anchor
ax.plot(
    yaw_fine,
    _ma_P0 * np.cos(np.radians(yaw_fine)) ** 1.5,
    color=PALETTE["ma2025"],
    lw=1.3,
    ls="-.",
    label="P₀·cos^1.5 (Ma fit)",
)
ax.scatter(
    _ma_yaw,
    _ma_P,
    color=PALETTE["ma2025"],
    s=55,
    zorder=6,
    marker="D",
    label="Ma et al. 2025 [LL-FVW]",
)
ax.set_xlabel("Ángulo de yaw [°]")
ax.set_ylabel("Potencia [MW]")
ax.set_title("Potencia aerodinámica")
ax.legend(frameon=False, fontsize=7)
ax.grid(alpha=0.3)
ax.set_xlim(-2, 43)

# Panel 3: Deflexión flapwise (media ± p-p/2) — both solvers are reliable here
ax = axes[2]
ax.errorbar(
    yaw_arr_c[mc["valid"]],
    mc["flap_m"][mc["valid"]],
    yerr=mc["flap_pp"][mc["valid"]] / 2,
    fmt="o",
    color=PALETTE["corrot"],
    capsize=4,
    lw=1.2,
    ms=6,
    label="Corrotacional media ± p-p/2",
)
ax.errorbar(
    yaw_arr_c[~mc["valid"]],
    mc["flap_m"][~mc["valid"]],
    yerr=mc["flap_pp"][~mc["valid"]] / 2,
    fmt="x",
    color=PALETTE["corrot"],
    capsize=4,
    lw=1.5,
    ms=8,
    markeredgewidth=2,
)
ax.errorbar(
    yaw_arr_i[mi["valid"]],
    mi["flap_m"][mi["valid"]],
    yerr=mi["flap_pp"][mi["valid"]] / 2,
    fmt="^",
    color=PALETTE["inert"],
    capsize=4,
    lw=1.2,
    ms=6,
    label="Inercial media ± p-p/2",
    alpha=0.85,
)
ax.errorbar(
    yaw_arr_i[~mi["valid"]],
    mi["flap_m"][~mi["valid"]],
    yerr=mi["flap_pp"][~mi["valid"]] / 2,
    fmt="x",
    color=PALETTE["inert"],
    capsize=4,
    lw=1.5,
    ms=8,
    markeredgewidth=2,
    alpha=0.85,
)
# Ma et al. 2025 reference series — Figure 16 mean blade-tip values at r/R=1
_ma_flap = np.array([ma_fig16_tip[f"flap_{int(yaw)}deg_m"] for yaw in _ma_yaw])
ax.plot(
    _ma_yaw,
    _ma_flap,
    color=PALETTE["ma2025"],
    lw=1.2,
    ls="-.",
    alpha=0.85,
    zorder=5,
)
ax.scatter(
    _ma_yaw,
    _ma_flap,
    color=PALETTE["ma2025"],
    s=55,
    zorder=6,
    marker="D",
    label="Ma Fig. 16 tip r/R=1",
)
ax.annotate(
    "Ma Fig. 16: valores medios de punta\nr/R=1, sin yaw=0",
    xy=(_ma_yaw[1], _ma_flap[1]),
    xytext=(16, 14.35),
    textcoords="data",
    fontsize=6.5,
    color=PALETTE["ma2025"],
    arrowprops={"arrowstyle": "->", "color": PALETTE["ma2025"], "lw": 0.8},
)
ax.axhline(13.86, color=PALETTE["zhou"], lw=1.3, ls="--", label="Zhou yaw=0° [13.86 m]")
ax.axhline(
    mc["flap_m"][0],
    color=PALETTE["ref"],
    lw=1.0,
    ls=":",
    label=f"Corrot. yaw=0° ({mc['flap_m'][0]:.2f} m)",
)
ax.set_xlabel("Ángulo de yaw [°]")
ax.set_ylabel("Deflexión flapwise tip media [m]")
ax.set_title("Deflexión flapwise tip — medias")
ax.legend(frameon=False, fontsize=7)
ax.grid(alpha=0.3)
ax.set_xlim(-2, 43)

fig.suptitle(
    "Barrido de yaw — Corrotacional vs Inercial — IEA 15 MW, condición rated\n"
    "△ = inercial  ○ = corrotacional  ◇ = Ma et al. 2025 [LL-FVW+GEBT]",
    fontsize=9,
    y=1.02,
)
fig.tight_layout()
out3 = OUT_DIR / "fig_5_2_yaw_sweep.png"
fig.savefig(out3, bbox_inches="tight")
fig.savefig(out3.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"  → {out3}")

# ============================================================
# Figure 4 — Series temporales de deflexión tip  §5.1
#   Eje X: tiempo [s] desde t=0 (incluye transitorio inicial)
#   Referencia literatura: Zhou LL-FVW+GEBT y ALM-GEBT (Tabla 4)
#   Referencia diseño: NREL DLC 1.4 pico máximo 22.8 m [R2-DLC]
# ============================================================

print("Generating fig_5_1_timeseries_deflection.png …")

# Cargar datos completos desde t=0 (sin trim)
cs_full = pd.read_csv(COROT_TIP)
is_full = pd.read_csv(INERT_TIP)
cs_full = cs_full[cs_full["Time [s]"] <= T_END].reset_index(drop=True)
is_full = is_full[is_full["Time [s]"] <= T_END].reset_index(drop=True)

# Ventana de estadísticas (régimen permanente)
cs_trim = steady_window(cs_full)
is_trim = steady_window(is_full)

# Medias de régimen permanente
mean_flap_c = float(cs_trim["Tip Disp Y [m]"].mean())
mean_flap_i = float(is_trim["Tip Disp Y [m]"].mean())
mean_edge_c = float(cs_trim["Tip Disp X [m]"].mean())
mean_edge_i = float(is_trim["Tip Disp X [m]"].mean())

T_MAX = max(cs_full["Time [s]"].max(), is_full["Time [s]"].max())

# ── Valores de referencia (Zhou Energy 2025, Tabla 4) ──────────────────────────
REF_FLAP = [
    (13.86, "Zhou LL-FVW+GEBT [13.86 m]", PALETTE["zhou"]),
    (14.10, "ALM-GEBT [14.10 m]", "#7f8c8d"),
]
REF_EDGE = [
    (-1.22, "Zhou LL-FVW+GEBT [−1.22 m]", PALETTE["zhou"]),
    (-1.27, "ALM-GEBT [−1.27 m]", "#7f8c8d"),
]
# Referencia de diseño (NREL DLC 1.4 — pico extremo, no condición comparable)
DLC_FLAP_MAX = 22.8  # m — NREL/TP-5000-75698 §6

fig, (ax_fl, ax_ed) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

# ─── Panel flapwise ─────────────────────────────────────────────────────────────
ax_fl.plot(
    cs_full["Time [s]"],
    cs_full["Tip Disp Y [m]"],
    lw=0.8,
    color=PALETTE["corrot"],
    label="AeroElast — Corotacional",
    zorder=3,
)
ax_fl.plot(
    is_full["Time [s]"],
    is_full["Tip Disp Y [m]"],
    lw=0.8,
    color=PALETTE["inert"],
    ls="--",
    label="AeroElast — Inercial",
    alpha=0.85,
    zorder=3,
)

# Línea de inicio de ventana estadística
ax_fl.axvline(T_START, color="#888", lw=1.0, ls=":", zorder=2)
ax_fl.text(
    T_START + 0.5,
    ax_fl.get_ylim()[0] if ax_fl.get_ylim()[0] > 0 else 1.0,
    f"t={T_START:.0f} s\nventana\nestadística",
    fontsize=6.5,
    color="#666",
    va="bottom",
)

# Referencia de diseño DLC 1.4
ax_fl.axhline(
    DLC_FLAP_MAX,
    color="#e67e22",
    lw=1.3,
    ls=(0, (4, 3)),
    zorder=2,
    label=f"NREL DLC 1.4 pico máx. [{DLC_FLAP_MAX} m] — no comparable",
)

# Medias de régimen permanente
ax_fl.axhline(
    mean_flap_c,
    color=PALETTE["corrot"],
    lw=1.1,
    ls=":",
    label=f"Media corrot. ({mean_flap_c:.2f} m)",
)
ax_fl.axhline(
    mean_flap_i,
    color=PALETTE["inert"],
    lw=1.1,
    ls=":",
    label=f"Media inercial ({mean_flap_i:.2f} m)",
)

# Referencias de literatura (régimen permanente)
for val, lbl, col in REF_FLAP:
    ax_fl.axhline(val, color=col, lw=1.4, ls=(0, (5, 3)), zorder=2)
    ax_fl.annotate(
        lbl,
        xy=(T_MAX * 0.99, val),
        xytext=(-4, 4),
        textcoords="offset points",
        fontsize=7,
        color=col,
        ha="right",
        va="bottom",
    )

ax_fl.set_ylabel("Flapwise tip deflection [m]")
ax_fl.set_title("(a) Deflexión flapwise — tip de pala, yaw=0°, condición rated", fontsize=9)
ax_fl.legend(frameon=False, loc="upper right", fontsize=7.0, ncol=2)
ax_fl.grid(alpha=0.3, lw=0.5)

# ─── Panel edgewise ─────────────────────────────────────────────────────────────
ax_ed.plot(
    cs_full["Time [s]"],
    cs_full["Tip Disp X [m]"],
    lw=0.8,
    color=PALETTE["corrot"],
    label="AeroElast — Corotacional",
    zorder=3,
)
ax_ed.plot(
    is_full["Time [s]"],
    is_full["Tip Disp X [m]"],
    lw=0.8,
    color=PALETTE["inert"],
    ls="--",
    label="AeroElast — Inercial",
    alpha=0.85,
    zorder=3,
)

# Línea de inicio de ventana estadística
ax_ed.axvline(T_START, color="#888", lw=1.0, ls=":", zorder=2)

# Medias de régimen permanente
ax_ed.axhline(
    mean_edge_c,
    color=PALETTE["corrot"],
    lw=1.1,
    ls=":",
    label=f"Media corrot. ({mean_edge_c:.2f} m)",
)
ax_ed.axhline(
    mean_edge_i,
    color=PALETTE["inert"],
    lw=1.1,
    ls=":",
    label=f"Media inercial ({mean_edge_i:.2f} m)",
)

# Referencias de literatura
for val, lbl, col in REF_EDGE:
    ax_ed.axhline(val, color=col, lw=1.4, ls=(0, (5, 3)), zorder=2)
    ax_ed.annotate(
        lbl,
        xy=(T_MAX * 0.99, val),
        xytext=(-4, 4),
        textcoords="offset points",
        fontsize=7,
        color=col,
        ha="right",
        va="bottom",
    )

ax_ed.set_xlabel("Tiempo [s]")
ax_ed.set_ylabel("Edgewise tip deflection [m]")
ax_ed.set_title("(b) Deflexión edgewise — tip de pala, yaw=0°, condición rated", fontsize=9)
ax_ed.legend(frameon=False, loc="lower right", fontsize=7.0, ncol=2)
ax_ed.grid(alpha=0.3, lw=0.5)

ax_ed.set_xlim(0, T_MAX)

# Anotación: diferencia edgewise con literatura
ax_ed.annotate(
    "Componente tangencial\naero. subestimada por BEM",
    xy=(T_MAX * 0.55, -1.22),
    xytext=(T_MAX * 0.48, -0.15),
    fontsize=7,
    color=PALETTE["zhou"],
    arrowprops={"arrowstyle": "->", "color": PALETTE["zhou"], "lw": 0.8},
    ha="center",
)

fig.suptitle(
    "Series temporales de deflexión de punta — IEA 15 MW, yaw=0°, condición rated\n"
    "Línea punteada vertical: inicio ventana estadística (t=40 s) · "
    "Líneas de puntos horizontales: medias y referencias Zhou 2025",
    fontsize=9.5,
)
fig.tight_layout()
out4 = OUT_DIR / "fig_5_1_timeseries_deflection.png"
fig.savefig(out4, bbox_inches="tight")
fig.savefig(out4.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig)
print(f"  → {out4}")


# ============================================================
# Figure 5 — Comparación con literatura
#   Flapwise tip | Edgewise tip | Mean thrust
# ============================================================

print("Generating fig_5_2b_literature_comparison.png …")

# Leer métricas de campaña (todos los campos estáticos y dinámicos)
_mc_df = summarize_corotational_campaign(COROT_ROOT)
_mi_df = pd.read_csv(METRICS_INERT).sort_values("yaw_deg").reset_index(drop=True)

# ---- Datos externos ----
# Ma et al. 2025 — digitized Figure 16, tip row r/R=1
_ma_yaw = np.array([10.0, 20.0, 30.0, 40.0])
_ma_flap_m = np.array([ma_fig16_tip[f"flap_{int(yaw)}deg_m"] for yaw in _ma_yaw])
_ma_edge_m = np.array([ma_fig16_tip[f"edge_{int(yaw)}deg_m"] for yaw in _ma_yaw])
_ma_thrust_mn = ma_fig15["mean_thrust_1e6N"].to_numpy(dtype=float)

# Zhou 2025 — solo yaw=0°, Tabla 4 (pala flexible, sin surge)
_zhou_flap_m0 = 13.86  # m
_zhou_edge_m0 = -1.22  # m
_zhou_thrust_mn0 = 2.20  # MN, Table 6 flexible blade


def _robust_mean(values) -> float:
    series = pd.to_numeric(pd.Series(values), errors="coerce").replace([np.inf, -np.inf], np.nan)
    series = series.dropna()
    if series.empty:
        return float("nan")
    median = series.median()
    mad = (series - median).abs().median()
    if mad > 0:
        clipped = series[(series - median).abs() <= 8.0 * 1.4826 * mad]
        if len(clipped) >= max(10, int(0.2 * len(series))):
            series = clipped
    return float(series.mean())


def _load_rotor_thrust_mn(root: pathlib.Path, solver: str, yaws: np.ndarray) -> np.ndarray:
    values = []
    for yaw in yaws:
        perf_path = root / f"yaw_{int(yaw)}" / solver / "rotor_performance.csv"
        perf = pd.read_csv(perf_path)
        perf = steady_window(perf)
        values.append(N_BLADES * _robust_mean(perf["Aero Thrust [N]"]) / 1e6)
    return np.array(values, dtype=float)


_mc_thrust_mn = _mc_df["thrust_mean"].to_numpy(dtype=float) / 1e6
_mi_thrust_mn = _load_rotor_thrust_mn(INERT_ROOT, "inertial", _mi_df["yaw_deg"].to_numpy())

fig5, axes5 = plt.subplots(1, 3, figsize=(13, 4.5))

# [0] Flapwise tip mean [m]
ax = axes5[0]
ax.plot(
    _mc_df["yaw_deg"],
    _mc_df["flap_mean"],
    "o-",
    color=PALETTE["corrot"],
    lw=1.5,
    ms=5,
    label="AeroElast corrot.",
)
ax.plot(
    _mi_df["yaw_deg"],
    _mi_df["flap_mean"],
    "^--",
    color=PALETTE["inert"],
    lw=1.5,
    ms=5,
    alpha=0.85,
    label="AeroElast inercial",
)
ax.scatter(
    _ma_yaw,
    _ma_flap_m,
    color=PALETTE["ma2025"],
    s=60,
    zorder=6,
    marker="D",
    label="Ma Fig. 16 tip r/R=1",
)
ax.plot(_ma_yaw, _ma_flap_m, color=PALETTE["ma2025"], lw=1.0, ls="-.", alpha=0.8)
ax.scatter(
    [0],
    [_zhou_flap_m0],
    color=PALETTE["zhou"],
    s=70,
    zorder=7,
    marker="s",
    label="Zhou 2025 (yaw=0°)",
)
for yaw, flap in zip(_ma_yaw, _ma_flap_m):
    ax.annotate(
        f"{flap:.1f}",
        xy=(yaw, flap),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=6.5,
        color=PALETTE["ma2025"],
    )
ax.set_title("(a) Deflexión flapwise tip r/R=1 [m]", fontsize=9)
ax.set_ylabel("Flapwise tip [m]")
ax.set_xlabel("Ángulo de yaw [°]")
ax.legend(frameon=False, fontsize=6.5)
ax.grid(alpha=0.3)
ax.set_xlim(-3, 44)

# [1] Signed edgewise tip mean [m]
ax = axes5[1]
ax.plot(
    _mc_df["yaw_deg"],
    _mc_df["edge_mean"],
    "o-",
    color=PALETTE["corrot"],
    lw=1.5,
    ms=5,
    label="AeroElast corrot.",
)
ax.plot(
    _mi_df["yaw_deg"],
    _mi_df["edge_mean"],
    "^--",
    color=PALETTE["inert"],
    lw=1.5,
    ms=5,
    alpha=0.85,
    label="AeroElast inercial",
)
ax.scatter(
    _ma_yaw,
    _ma_edge_m,
    color=PALETTE["ma2025"],
    s=60,
    zorder=6,
    marker="D",
    label="Ma Fig. 16 tip r/R=1",
)
ax.plot(_ma_yaw, _ma_edge_m, color=PALETTE["ma2025"], lw=1.0, ls="-.", alpha=0.8)
ax.scatter(
    [0],
    [_zhou_edge_m0],
    color=PALETTE["zhou"],
    s=70,
    zorder=7,
    marker="s",
    label="Zhou 2025 (yaw=0°)",
)
for yaw, edge in zip(_ma_yaw, _ma_edge_m):
    ax.annotate(
        f"{edge:.2f}",
        xy=(yaw, edge),
        xytext=(0, -12),
        textcoords="offset points",
        ha="center",
        fontsize=6.5,
        color=PALETTE["ma2025"],
    )
ax.set_title("(b) Deflexión edgewise tip r/R=1 [m]", fontsize=9)
ax.set_ylabel("Edgewise tip [m]")
ax.set_xlabel("Ángulo de yaw [°]")
ax.legend(frameon=False, fontsize=6.5)
ax.grid(alpha=0.3)
ax.set_xlim(-3, 44)
edge_min = min(
    float(_ma_edge_m.min()),
    float(_mc_df["edge_mean"].min()),
    float(_mi_df["edge_mean"].min()),
    _zhou_edge_m0,
)
edge_max = max(
    float(_ma_edge_m.max()),
    float(_mc_df["edge_mean"].max()),
    float(_mi_df["edge_mean"].max()),
    _zhou_edge_m0,
)
ax.set_ylim(edge_min - 0.2, edge_max + 0.2)
ax.annotate(
    "Valores Ma conservan el signo\nde Fig. 16",
    xy=(0.02, 0.04),
    xycoords="axes fraction",
    fontsize=6.5,
    color="#777",
)

# [2] Mean aerodynamic thrust [MN]
ax = axes5[2]
ax.plot(
    _mc_df["yaw_deg"],
    _mc_thrust_mn,
    "o-",
    color=PALETTE["corrot"],
    lw=1.5,
    ms=5,
    label="AeroElast corrot.",
)
ax.plot(
    _mi_df["yaw_deg"],
    _mi_thrust_mn,
    "^--",
    color=PALETTE["inert"],
    lw=1.5,
    ms=5,
    alpha=0.85,
    label="AeroElast inercial",
)
ax.scatter(
    _ma_yaw,
    _ma_thrust_mn,
    color=PALETTE["ma2025"],
    s=60,
    zorder=6,
    marker="D",
    label="Ma Fig. 15 thrust",
)
ax.plot(_ma_yaw, _ma_thrust_mn, color=PALETTE["ma2025"], lw=1.0, ls="-.", alpha=0.8)
ax.scatter(
    [0],
    [_zhou_thrust_mn0],
    color=PALETTE["zhou"],
    s=70,
    zorder=7,
    marker="s",
    label="Zhou 2025 (yaw=0°)",
)
for yaw, thrust in zip(_ma_yaw, _ma_thrust_mn):
    ax.annotate(
        f"{thrust:.2f}",
        xy=(yaw, thrust),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center",
        fontsize=6.5,
        color=PALETTE["ma2025"],
    )
ax.set_title("(c) Thrust medio del rotor [MN]", fontsize=9)
ax.set_ylabel("Thrust [MN]")
ax.set_xlabel("Ángulo de yaw [°]")
ax.legend(frameon=False, fontsize=6.5)
ax.grid(alpha=0.3)
ax.set_xlim(-3, 44)
thrust_min = min(float(_mc_thrust_mn.min()), float(_mi_thrust_mn.min()), float(_ma_thrust_mn.min()))
thrust_max = max(
    float(_mc_thrust_mn.max()),
    float(_mi_thrust_mn.max()),
    float(_ma_thrust_mn.max()),
    _zhou_thrust_mn0,
)
ax.set_ylim(thrust_min - 0.1, thrust_max + 0.15)
ax.annotate(
    "AeroElast: rotor_performance\nMa: Fig. 15",
    xy=(0.02, 0.04),
    xycoords="axes fraction",
    fontsize=6.5,
    color="#777",
)

fig5.suptitle(
    "IEA 15 MW — Respuesta de punta r/R=1 y thrust medio contra yaw\n"
    "◇ = Ma Fig. 16 en (a,b), Ma Fig. 15 en (c) · □ = Zhou 2025 (solo yaw=0°)",
    fontsize=9.5,
    y=1.03,
)
fig5.tight_layout()
out5 = OUT_DIR / "fig_5_2b_literature_comparison.png"
fig5.savefig(out5, bbox_inches="tight")
fig5.savefig(out5.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig5)
print(f"  → {out5}")


# ============================================================
# Figure 5c — Ma et al. 2025 Figure 16 digitization (spanwise)
#   Overlays AeroElast spanwise mean profiles only for the structural
#   quantities that are directly comparable from the VTU fields.
# ============================================================

print("Generating fig_5_2c_ma_fig16_digitization.png …")

ma_yaws = [10, 20, 30, 40]
ma_colors = {10: "#e74c3c", 20: "#2980b9", 30: "#f1c40f", 40: "#2ecc71"}
ma_markers = {10: "s", 20: "o", 30: "^", 40: "v"}

fig6, axes6 = plt.subplots(1, 2, figsize=(11.6, 4.6), sharex=True)
ma_r = ma_fig16["r_R"]
corot_profiles = {}
inert_profiles = {}

for yaw in ma_yaws:
    corot_profiles[yaw] = load_spanwise_vtu_profile(
        COROT_ROOT / f"yaw_{yaw}/corotational",
        ma_r.to_numpy(dtype=float),
    )
    inert_profiles[yaw] = load_spanwise_vtu_profile(
        INERT_ROOT / f"yaw_{yaw}/inertial",
        ma_r.to_numpy(dtype=float),
    )
    print(
        f"    yaw={yaw:02d}° spanwise profiles: "
        f"corot frames={corot_profiles[yaw]['n_frames']}, edge_sign={corot_profiles[yaw]['edge_sign']:+.0f}; "
        f"inert frames={inert_profiles[yaw]['n_frames']}, edge_sign={inert_profiles[yaw]['edge_sign']:+.0f}, "
        f"t_end={inert_profiles[yaw]['t_end']:.2f}s"
    )

for yaw in ma_yaws:
    ref_style = {
        "color": ma_colors[yaw],
        "lw": 1.3,
        "ls": "--",
        "marker": ma_markers[yaw],
        "markevery": 2,
        "ms": 4,
    }
    corot_style = {"color": ma_colors[yaw], "lw": 1.6, "ls": "-"}
    inert_style = {"color": ma_colors[yaw], "lw": 1.3, "ls": "-."}

    axes6[0].plot(ma_r, ma_fig16[f"flap_{yaw}deg_m"], **ref_style)
    axes6[0].plot(ma_r, corot_profiles[yaw]["flap"], **corot_style)
    axes6[0].plot(ma_r, inert_profiles[yaw]["flap"], **inert_style)

    axes6[1].plot(ma_r, ma_fig16[f"edge_{yaw}deg_m"], **ref_style)
    axes6[1].plot(ma_r, corot_profiles[yaw]["edge"], **corot_style)
    axes6[1].plot(ma_r, inert_profiles[yaw]["edge"], **inert_style)

axes6[0].set_title("(a) Deflexión flapwise")
axes6[0].set_ylabel("Deflexión flapwise [m]")
axes6[0].set_ylim(0, 15)

axes6[1].set_title("(b) Deflexión edgewise")
axes6[1].set_ylabel("Deflexión edgewise [m]")
axes6[1].set_ylim(-2.5, 0.5)

for ax in axes6:
    ax.set_xlabel("Posición spanwise r/R")
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)

source_handles = [
    Line2D(
        [0],
        [0],
        color="#555",
        lw=1.3,
        ls="--",
        marker="o",
        ms=4,
        label="Ma et al. 2025 (digitalizado)",
    ),
    Line2D([0], [0], color="#555", lw=1.6, ls="-", label="AeroElast corrotacional (VTU)"),
    Line2D([0], [0], color="#555", lw=1.3, ls="-.", label="AeroElast inercial (VTU)"),
]
yaw_handles = [
    Line2D([0], [0], color=ma_colors[yaw], lw=1.5, label=f"yaw={yaw}°") for yaw in ma_yaws
]
source_legend = fig6.legend(
    handles=source_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.95),
    ncol=3,
    frameon=False,
)
fig6.add_artist(source_legend)
fig6.legend(
    handles=yaw_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.88),
    ncol=4,
    frameon=False,
)

fig6.suptitle(
    "Ma et al. 2025 Figure 16 vs AeroElast — respuesta spanwise media\n"
    "Solo se muestran flapwise y edgewise, que sí tienen una contraparte VTU comparable; la curva yaw=0 no existe en la figura original",
    fontsize=9.5,
    y=1.02,
)
fig6.tight_layout(rect=(0, 0, 1, 0.82))
out6 = OUT_DIR / "fig_5_2c_ma_fig16_digitization.png"
fig6.savefig(out6, bbox_inches="tight")
fig6.savefig(out6.with_suffix(".pdf"), bbox_inches="tight")
plt.close(fig6)
print(f"  → {out6}")


print("\nDone. All figures written to:", OUT_DIR)
