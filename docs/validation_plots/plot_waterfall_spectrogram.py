"""
Short-Time Fourier Transform (waterfall) spectrogram — IEA 15 MW, rated, yaw=0°.

Visualises how the dominant frequencies of the flapwise tip displacement
evolve from the initial transient (t < 20 s) to the steady-state regime
(t ∈ [40, 100] s). The plot complements the static FFT of §5.5: the FFT collapses
the time axis, while this waterfall shows the transient settling and the
persistence of the 1P, 3P and structural-mode peaks.

Outputs
-------
docs/figures/fig_5_4_11_waterfall_spectrogram.png  (300 dpi)
docs/figures/fig_5_4_11_waterfall_spectrogram.pdf
"""
import base64
import mmap
import pathlib
import re
import zlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import spectrogram, windows

# ---------------------------------------------------------------------------
# Paths and parameters
# ---------------------------------------------------------------------------
BASE = pathlib.Path("/scratch/leahk/eduardo.donestevez")
COROT_CASE = BASE / "frontiersin_results_corotational_100s/yaw_0/corotational"
INERT_CASE = BASE / "frontiersin_results_inertial/yaw_0/inertial"

WINDOW_SEC = 8.0         # s, STFT window length (~1 P period × 1 s safety)
OVERLAP_FRAC = 0.875     # 87.5 % overlap → high time resolution
FMAX = 1.5               # Hz, frequency axis upper limit
T_END = 100.0            # s, end of 100 s corotacional campaign

# Characteristic frequencies (Hz) for vertical reference
F_1P = 0.12527
F_3P = 0.37581
F_FLAP = 0.5537
F_EDGE = 0.6290

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _time_value(path: pathlib.Path) -> float:
    return float(path.parent.name)


def _field_vtu_paths(case_dir: pathlib.Path) -> list[pathlib.Path]:
    return sorted(case_dir.glob("*/fields.vtu"), key=_time_value)


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
        raise TypeError(f"Unsupported VTU array type for {name!r}")
    n_components_match = re.search(r'NumberOfComponents="(\d+)"', tag)
    n_components = int(n_components_match.group(1)) if n_components_match else 1
    values = _decode_vtu_payload(encoded, np.dtype("<f8"))
    if n_components > 1:
        values = values.reshape((-1, n_components))
    return values.astype(float, copy=False)


def _select_tip_node(path: pathlib.Path) -> tuple[int, float, int]:
    points = _vtu_data_array(path, "Points")
    displacement = _vtu_data_array(path, "U")
    points0 = points - displacement
    z0 = points0[:, 2]
    z_tip = float(z0.max())
    candidates = np.where(np.isclose(z0, z_tip, rtol=0.0, atol=1e-8))[0]
    centroid_xy = points0[candidates, :2].mean(axis=0)
    distances = np.sum((points0[candidates, :2] - centroid_xy) ** 2, axis=1)
    return int(candidates[int(np.argmin(distances))]), z_tip, int(len(candidates))


def load_tip_vtu(case_dir: pathlib.Path) -> pd.DataFrame:
    paths = _field_vtu_paths(case_dir)
    if not paths:
        raise FileNotFoundError(f"No fields.vtu files found under {case_dir}")
    paths = [path for path in paths if _time_value(path) <= T_END]
    tip_node, z_tip, n_candidates = _select_tip_node(paths[0])
    rows = []
    for path in paths:
        displacement = _vtu_data_array(path, "U")[tip_node]
        rows.append(
            {
                "Time [s]": _time_value(path),
                "Tip Disp Y [m]": displacement[1],
                "tip_node": tip_node,
                "tip_z0": z_tip,
                "tip_z0_candidates": n_candidates,
            }
        )
    return pd.DataFrame(rows).sort_values("Time [s]").reset_index(drop=True)


def _spectrogram(t, y, win_sec=WINDOW_SEC, overlap_frac=OVERLAP_FRAC):
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt
    nperseg = int(win_sec * fs)
    noverlap = int(nperseg * overlap_frac)
    f, tspec, Sxx = spectrogram(
        y - y.mean(),
        fs=fs,
        window=windows.hann(nperseg),
        nperseg=nperseg,
        noverlap=noverlap,
        scaling="spectrum",
        mode="magnitude",
    )
    return f, tspec + float(t[0]), Sxx, dt


# ---------------------------------------------------------------------------
# Compute spectrograms
# ---------------------------------------------------------------------------
corot = load_tip_vtu(COROT_CASE)
inert = load_tip_vtu(INERT_CASE)
print(
    "Fixed VTU tip node: "
    f"corot={int(corot['tip_node'].iloc[0])} "
    f"(Z0={corot['tip_z0'].iloc[0]:.3f}, candidates={int(corot['tip_z0_candidates'].iloc[0])}), "
    f"inert={int(inert['tip_node'].iloc[0])} "
    f"(Z0={inert['tip_z0'].iloc[0]:.3f}, candidates={int(inert['tip_z0_candidates'].iloc[0])})"
)

f_c, t_c, S_c, dt_c = _spectrogram(corot["Time [s]"].values,
                                   corot["Tip Disp Y [m]"].values)
f_i, t_i, S_i, dt_i = _spectrogram(inert["Time [s]"].values,
                                   inert["Tip Disp Y [m]"].values)

# ---------------------------------------------------------------------------
# Plot 2×1 — corotational top, inertial bottom
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=False)
fig.suptitle(
    "Fixed-tip flapwise displacement spectrogram — IEA 15 MW, yaw=0°, rated\n"
    f"STFT (Hann, {WINDOW_SEC:.0f} s window, {int(OVERLAP_FRAC*100)}% overlap)",
    fontsize=11, y=1.00,
)

for ax, label, f, t, S, cmap in [
    (axes[0], "Corotational", f_c, t_c, S_c, "magma"),
    (axes[1], "Inertial",     f_i, t_i, S_i, "viridis"),
]:
    f_mask = f <= FMAX
    f_lim = f[f_mask]
    S_lim = S[f_mask, :]
    pcm = ax.imshow(
        100 * S_lim,
        origin="lower",
        aspect="auto",
        extent=(float(t[0]), float(t[-1]), float(f_lim[0]), FMAX),
        interpolation="bilinear",
        cmap=cmap,
    )
    cb = plt.colorbar(pcm, ax=ax, label="Magnitude [cm]")
    ax.set_ylabel(f"Frequency [Hz]\n({label})")
    ax.set_ylim(0, FMAX)
    ax.set_xlim(float(t[0]), float(t[-1]))
    for fname, fval, color in [("1P", F_1P, "white"),
                                ("3P", F_3P, "white"),
                                ("f_1,flap", F_FLAP, "cyan"),
                                ("f_1,edge", F_EDGE, "yellow")]:
        ax.axhline(fval, color=color, ls=":", lw=0.8, alpha=0.8)
        ax.text(t[-1] * 0.99, fval + 0.02, fname, color=color, fontsize=7,
                ha="right", va="bottom")

axes[0].set_xlabel("Time [s]")
axes[1].set_xlabel("Time [s]")
fig.tight_layout()
png = OUT_DIR / "fig_5_4_11_waterfall_spectrogram.png"
pdf = OUT_DIR / "fig_5_4_11_waterfall_spectrogram.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
print(f"Saved: {png}")
print(f"Saved: {pdf}")
print()
print(f"Corotational: time range [{t_c.min():.1f}, {t_c.max():.1f}] s, "
      f"freq resolution = {f_c[1]-f_c[0]:.3f} Hz, time resolution = {t_c[1]-t_c[0]:.2f} s")
print(f"Inertial:     time range [{t_i.min():.1f}, {t_i.max():.1f}] s")
print(f"VTU sampling: corot dt={dt_c:.3f} s, inert dt={dt_i:.3f} s")
