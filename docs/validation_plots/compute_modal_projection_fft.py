"""Project corotational checkpoint response onto a reduced modal basis and FFT it.

This script is the *viable first-pass* version of the modal-response analysis
for the Frontiers corotational campaign:

- build one reduced modal basis from the exact structural model used by the
  campaign (same YAML / mesh / BCs);
- project checkpointed reduced displacements ``u_red(t)`` onto that basis;
- compute one-sided FFTs of the modal coordinates ``q_i(t)`` for each yaw case.

It does **not** solve a new eigenproblem at every checkpoint. That stronger
approach would require mode tracking (MAC / sign / ordering continuity) across
time, which is a separate analysis phase. The present script is intended for
modal-family attribution: *which spatial mode family carries the peaks already
seen in the blade-tip FFT?*

Basis options
-------------
``unloaded``
    Solve modes from the reduced elastic stiffness and lumped mass matrices.
``kg_static``
    Add a single geometric-stiffness matrix ``K_G(omega_ref)`` built from the
    campaign rotor speed before solving the basis once. This is closer to the
    rotating operating point, but still not a full loaded tangent basis because
    it does not include per-checkpoint mode updates or explicit spin softening
    as a standalone exported matrix.

Outputs
-------
docs/validation_data/generated/modal_projection_fft_<tag>_timeseries.csv
docs/validation_data/generated/modal_projection_fft_<tag>_peaks.csv
docs/validation_data/generated/modal_projection_fft_<tag>_markers.csv
docs/validation_data/generated/modal_projection_fft_<tag>_summary.md
docs/figures/fig_modal_projection_fft_<tag>_yaw0.png
docs/figures/fig_modal_projection_fft_<tag>_yaw0.pdf
docs/figures/fig_modal_projection_fft_<tag>_peak_tracking.png
docs/figures/fig_modal_projection_fft_<tag>_peak_tracking.pdf
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_ROOT = Path("/scratch/leahk/eduardo.donestevez/frontiersin_results_corotational_100s")
DEFAULT_CONFIG_PATH = (
    REPO_ROOT / "tests" / "IEA15MW" / "frontiersin_2025_yaw" / "solid_corotational_yaw_0.yaml"
)
DEFAULT_OUT_DATA = REPO_ROOT / "docs" / "validation_data" / "generated"
DEFAULT_OUT_FIG = REPO_ROOT / "docs" / "figures"

DEFAULT_YAWS = (0, 10, 20, 30, 40)
DEFAULT_NUM_MODES = 6
DEFAULT_START_TIME = 40.0
DEFAULT_END_TIME = 100.0
DEFAULT_MAX_FREQ = 2.5
MIN_PEAK_FREQ = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=DEFAULT_RESULTS_ROOT,
        help="Campaign root containing yaw_*/corotational/<time>/state.npz.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Structural YAML used to assemble the campaign model.",
    )
    parser.add_argument(
        "--yaws",
        type=int,
        nargs="+",
        default=list(DEFAULT_YAWS),
        help="Yaw angles to analyse.",
    )
    parser.add_argument(
        "--num-modes",
        type=int,
        default=DEFAULT_NUM_MODES,
        help="Number of modes to retain in the reduced basis.",
    )
    parser.add_argument(
        "--basis-mode",
        choices=("unloaded", "kg_static"),
        default="unloaded",
        help="Basis used for the modal projection.",
    )
    parser.add_argument("--start-time", type=float, default=DEFAULT_START_TIME)
    parser.add_argument("--end-time", type=float, default=DEFAULT_END_TIME)
    parser.add_argument(
        "--checkpoint-stride",
        type=int,
        default=1,
        help="Keep every N-th checkpoint after time filtering.",
    )
    parser.add_argument(
        "--max-checkpoints",
        type=int,
        default=None,
        help="Optional cap used mainly for smoke tests; samples the filtered window evenly.",
    )
    parser.add_argument(
        "--plot-modes",
        type=int,
        default=4,
        help="Number of leading modes to show in the figures.",
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=DEFAULT_MAX_FREQ,
        help="Maximum frequency shown in the FFT figures [Hz].",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip figure generation and write only tabular outputs.",
    )
    return parser.parse_args()


def import_runtime() -> tuple[type, object]:
    sys.path.insert(0, str((REPO_ROOT / "src").resolve()))
    try:
        from aeroelast.solvers.fsi.runner import FSIRunner
        from _aeroelast import modal_solve_coo
    except ImportError as exc:  # pragma: no cover - cluster/runtime specific
        raise RuntimeError(
            "Could not import the AeroElast runtime. Load the SDumont GCC/GLU modules "
            "and ensure the compiled _aeroelast extension is on LD_LIBRARY_PATH."
        ) from exc
    return FSIRunner, modal_solve_coo


def petsc_to_coo(mat) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ai, aj, av = mat.getValuesCSR()
    counts = np.diff(ai.astype(np.int64))
    rows = np.repeat(np.arange(len(counts), dtype=np.int32), counts)
    cols = aj.astype(np.int32)
    vals = av.astype(np.float64)
    return rows, cols, vals


def csr_from_petsc(mat) -> sp.csr_matrix:
    indptr, indices, data = mat.getValuesCSR()
    return sp.csr_matrix((data, indices, indptr), shape=mat.getSize())


def compute_participation_factors(
    modes_red: np.ndarray,
    mass_red: sp.csr_matrix,
    free_dofs_global: np.ndarray,
    coords: np.ndarray,
    dofs_per_node: int,
    n_dof_total: int,
) -> tuple[np.ndarray, list[str]]:
    direction_labels = ["Tx", "Ty", "Tz"][: min(dofs_per_node, 3)]
    rigid_vectors: list[np.ndarray] = []

    for direction in range(min(dofs_per_node, 3)):
        r_full = np.zeros(n_dof_total, dtype=np.float64)
        r_full[direction::dofs_per_node] = 1.0
        rigid_vectors.append(r_full[free_dofs_global])

    if dofs_per_node >= 6:
        direction_labels += ["Rx", "Ry", "Rz"]
        for direction in range(3):
            r_full = np.zeros(n_dof_total, dtype=np.float64)
            if direction == 0:
                r_full[1::dofs_per_node] = coords[:, 2]
                r_full[2::dofs_per_node] = -coords[:, 1]
                r_full[3::dofs_per_node] = 1.0
            elif direction == 1:
                r_full[0::dofs_per_node] = -coords[:, 2]
                r_full[2::dofs_per_node] = coords[:, 0]
                r_full[4::dofs_per_node] = 1.0
            else:
                r_full[0::dofs_per_node] = -coords[:, 1]
                r_full[1::dofs_per_node] = coords[:, 0]
                r_full[5::dofs_per_node] = 1.0
            rigid_vectors.append(r_full[free_dofs_global])

    total_directional_mass = np.array([float(r @ (mass_red @ r)) for r in rigid_vectors])
    factors = np.zeros((modes_red.shape[1], len(rigid_vectors)), dtype=np.float64)

    for mode_idx in range(modes_red.shape[1]):
        phi = modes_red[:, mode_idx]
        m_phi = mass_red @ phi
        generalized_mass = float(phi @ m_phi)
        for dir_idx, rigid_vec in enumerate(rigid_vectors):
            directional_mass = total_directional_mass[dir_idx]
            if generalized_mass <= 0.0 or directional_mass <= 0.0:
                continue
            gamma = float(phi @ (mass_red @ rigid_vec))
            factors[mode_idx, dir_idx] = 100.0 * (gamma * gamma) / (generalized_mass * directional_mass)

    return factors, direction_labels


def classify_modes(factors: np.ndarray, direction_labels: Sequence[str]) -> tuple[list[str], list[str]]:
    counters = {"flapwise": 0, "edgewise": 0, "torsion": 0, "axial": 0}
    labels: list[str] = []
    characters: list[str] = []

    for row in factors:
        values = {label: row[idx] for idx, label in enumerate(direction_labels)}
        dominant = max(values, key=values.get)
        if dominant in ("Rx", "Ty"):
            character = "flapwise"
            prefix = "F"
        elif dominant in ("Ry", "Tx"):
            character = "edgewise"
            prefix = "E"
        elif dominant == "Rz":
            character = "torsion"
            prefix = "T"
        else:
            character = "axial"
            prefix = "A"

        counters[character] += 1
        labels.append(f"{counters[character]}{prefix}")
        characters.append(character)

    return labels, characters


def canonicalise_modes(modes_red: np.ndarray, mass_red: sp.csr_matrix) -> np.ndarray:
    modes = np.array(modes_red, dtype=np.float64, copy=True)
    for idx in range(modes.shape[1]):
        vec = modes[:, idx]
        anchor = int(np.argmax(np.abs(vec)))
        if vec[anchor] < 0.0:
            vec *= -1.0
        norm = math.sqrt(float(vec @ (mass_red @ vec)))
        if norm <= 0.0:
            raise RuntimeError(f"Mode {idx + 1} has non-positive generalized mass.")
        modes[:, idx] = vec / norm
    return modes


def build_modal_basis(
    config_path: Path,
    num_modes: int,
    basis_mode: str,
) -> dict[str, object]:
    FSIRunner, modal_solve_coo = import_runtime()

    runner = FSIRunner(config_path, working_dir=config_path.parent)
    runner.mesh = runner._setup_mesh()
    runner._material = runner._create_material()
    model_config = runner._build_model_config()
    runner.solver = runner._create_solver(model_config)
    runner._apply_boundary_conditions()

    solver = runner.solver
    omega_ref = float(solver.solver_params["rotor"]["omega"])
    omega_rpm = omega_ref * 60.0 / (2.0 * np.pi)

    damping_enabled = getattr(solver, "_damping_enabled", None)
    damping_auto = getattr(solver, "_damping_auto", None)
    if damping_enabled is not None:
        solver._damping_enabled = False
    if damping_auto is not None:
        solver._damping_auto = False

    try:
        (stiffness, mass), bc_manager, k_geom = solver._assemble_system_matrices(omega_ref)
    finally:
        if damping_enabled is not None:
            solver._damping_enabled = damping_enabled
        if damping_auto is not None:
            solver._damping_auto = damping_auto

    k_red = bc_manager.reduce_matrix(stiffness)
    m_red = bc_manager.reduce_matrix(mass)

    if basis_mode == "kg_static":
        if k_geom is None:
            raise RuntimeError("basis_mode=kg_static requested, but K_G was not assembled.")
        k_geom_red = bc_manager.reduce_matrix(k_geom)
        k_red.axpy(1.0, k_geom_red)

    mass_red_sp = csr_from_petsc(m_red)
    k_rows, k_cols, k_vals = petsc_to_coo(k_red)
    m_rows, m_cols, m_vals = petsc_to_coo(m_red)

    n_free = k_red.getSize()[0]
    reduced_dofs = np.arange(n_free, dtype=np.int64)
    frequencies_hz, modes_flat = modal_solve_coo(
        k_rows.astype(np.int64),
        k_cols.astype(np.int64),
        k_vals.astype(np.float64),
        m_rows.astype(np.int64),
        m_cols.astype(np.int64),
        m_vals.astype(np.float64),
        n_free,
        reduced_dofs,
        num_modes,
    )

    frequencies_hz = np.asarray(frequencies_hz, dtype=np.float64)
    if frequencies_hz.size < num_modes:
        raise RuntimeError(
            f"Requested {num_modes} modes, but only {frequencies_hz.size} converged."
        )

    modes_red = np.asarray(modes_flat, dtype=np.float64).reshape(frequencies_hz.size, n_free).T
    modes_red = canonicalise_modes(modes_red, mass_red_sp)

    free_dofs_global = np.asarray(sorted(bc_manager.free_dofs), dtype=np.int64)
    factors, direction_labels = compute_participation_factors(
        modes_red,
        mass_red_sp,
        free_dofs_global,
        runner.mesh.coords_array,
        solver.domain.dofs_per_node,
        solver.domain.dofs_count,
    )
    mode_labels, mode_characters = classify_modes(factors, direction_labels)

    modal_projection_operator = mass_red_sp @ modes_red

    return {
        "frequencies_hz": frequencies_hz,
        "modes_red": modes_red,
        "mass_red": mass_red_sp,
        "projection_operator": modal_projection_operator,
        "participation_factors": factors,
        "direction_labels": direction_labels,
        "mode_labels": mode_labels,
        "mode_characters": mode_characters,
        "omega_ref_rad_s": omega_ref,
        "omega_ref_rpm": omega_rpm,
        "config_path": config_path,
        "basis_mode": basis_mode,
    }


def list_time_dirs(case_dir: Path) -> list[tuple[float, Path]]:
    time_dirs: list[tuple[float, Path]] = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            time_value = float(child.name)
        except ValueError:
            continue
        if (child / "state.npz").exists():
            time_dirs.append((time_value, child))
    time_dirs.sort(key=lambda item: item[0])
    return time_dirs


def evenly_sample(items: Sequence[tuple[float, Path]], max_items: int | None) -> list[tuple[float, Path]]:
    if max_items is None or len(items) <= max_items:
        return list(items)
    if max_items <= 1:
        return [items[-1]]
    raw_idx = np.linspace(0, len(items) - 1, max_items)
    indices = np.unique(np.round(raw_idx).astype(int))
    return [items[idx] for idx in indices]


def select_time_dirs(
    case_dir: Path,
    start_time: float,
    end_time: float,
    stride: int,
    max_checkpoints: int | None,
) -> list[tuple[float, Path]]:
    filtered = [
        item
        for item in list_time_dirs(case_dir)
        if start_time <= item[0] <= end_time
    ]
    if stride > 1:
        filtered = filtered[::stride]
    filtered = evenly_sample(filtered, max_checkpoints)
    if not filtered:
        raise RuntimeError(
            f"No checkpoint states found in {case_dir} within [{start_time}, {end_time}] s."
        )
    return filtered


def load_u_red(state_path: Path) -> np.ndarray:
    with np.load(state_path) as state:
        return np.asarray(state["u_red"], dtype=np.float64)


def compute_mean_displacement(selected_dirs: Sequence[tuple[float, Path]]) -> np.ndarray:
    accum: np.ndarray | None = None
    count = 0
    for _, folder in selected_dirs:
        disp = load_u_red(folder / "state.npz")
        if accum is None:
            accum = np.zeros_like(disp, dtype=np.float64)
        accum += disp
        count += 1
    assert accum is not None
    return accum / float(count)


def project_case_timeseries(
    selected_dirs: Sequence[tuple[float, Path]],
    mean_disp: np.ndarray,
    projection_operator: np.ndarray | sp.spmatrix,
) -> tuple[np.ndarray, np.ndarray]:
    times = np.zeros(len(selected_dirs), dtype=np.float64)
    q_series = np.zeros((len(selected_dirs), projection_operator.shape[1]), dtype=np.float64)
    if sp.issparse(projection_operator):
        psi = projection_operator.toarray()
    else:
        psi = np.asarray(projection_operator, dtype=np.float64)

    for idx, (time_value, folder) in enumerate(selected_dirs):
        disp = load_u_red(folder / "state.npz") - mean_disp
        times[idx] = time_value
        q_series[idx, :] = disp @ psi

    return times, q_series


def onesided_fft(signal: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(signal, dtype=np.float64)
    if x.ndim != 1:
        raise ValueError("onesided_fft expects a 1D signal")
    x = x - x.mean()
    n = x.size
    if n < 2:
        raise ValueError("Need at least two samples to compute an FFT")

    window = np.hanning(n)
    spectrum = np.fft.rfft(x * window)
    freq = np.fft.rfftfreq(n, d=dt)

    amplitude = 2.0 * np.abs(spectrum) / window.sum()
    amplitude[0] *= 0.5
    if n % 2 == 0 and amplitude.size > 1:
        amplitude[-1] *= 0.5
    return freq, amplitude


def nearest_bin_amplitude(freq: np.ndarray, amp: np.ndarray, target_hz: float) -> float:
    idx = int(np.argmin(np.abs(freq - target_hz)))
    return float(amp[idx])


def dominant_peak(freq: np.ndarray, amp: np.ndarray) -> tuple[float, float]:
    mask = freq >= MIN_PEAK_FREQ
    if not np.any(mask):
        idx = int(np.argmax(amp))
    else:
        idx_local = int(np.argmax(amp[mask]))
        idx = int(np.flatnonzero(mask)[idx_local])
    return float(freq[idx]), float(amp[idx])


def nearest_p_label(freq_hz: float, one_p_hz: float) -> str:
    if one_p_hz <= 0.0:
        return "-"
    harmonic = int(round(freq_hz / one_p_hz))
    if harmonic < 1:
        return "-"
    return f"{harmonic}P"


def make_output_paths(tag: str) -> dict[str, Path]:
    DEFAULT_OUT_DATA.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT_FIG.mkdir(parents=True, exist_ok=True)
    return {
        "timeseries_csv": DEFAULT_OUT_DATA / f"modal_projection_fft_{tag}_timeseries.csv",
        "peaks_csv": DEFAULT_OUT_DATA / f"modal_projection_fft_{tag}_peaks.csv",
        "markers_csv": DEFAULT_OUT_DATA / f"modal_projection_fft_{tag}_markers.csv",
        "summary_md": DEFAULT_OUT_DATA / f"modal_projection_fft_{tag}_summary.md",
        "yaw0_png": DEFAULT_OUT_FIG / f"fig_modal_projection_fft_{tag}_yaw0.png",
        "yaw0_pdf": DEFAULT_OUT_FIG / f"fig_modal_projection_fft_{tag}_yaw0.pdf",
        "tracking_png": DEFAULT_OUT_FIG / f"fig_modal_projection_fft_{tag}_peak_tracking.png",
        "tracking_pdf": DEFAULT_OUT_FIG / f"fig_modal_projection_fft_{tag}_peak_tracking.pdf",
    }


def plot_yaw0_fft(
    freq_map: dict[int, np.ndarray],
    amp_map: dict[int, np.ndarray],
    basis_info: dict[str, object],
    max_freq: float,
    plot_modes: int,
    output_png: Path,
    output_pdf: Path,
) -> None:
    yaw0 = 0
    if yaw0 not in freq_map:
        return

    frequencies_hz = np.asarray(basis_info["frequencies_hz"], dtype=np.float64)
    mode_labels = list(basis_info["mode_labels"])
    mode_characters = list(basis_info["mode_characters"])
    one_p_hz = float(basis_info["omega_ref_rpm"]) / 60.0
    markers = {
        "1P": one_p_hz,
        "2P": 2.0 * one_p_hz,
        "3P": 3.0 * one_p_hz,
        "4P": 4.0 * one_p_hz,
    }

    n_plot = min(plot_modes, len(mode_labels))
    fig, axes = plt.subplots(n_plot, 1, figsize=(10.4, 2.4 * n_plot), sharex=True)
    if n_plot == 1:
        axes = [axes]

    freq = freq_map[yaw0]
    amp = amp_map[yaw0]

    for mode_idx, ax in enumerate(axes):
        ax.plot(freq, amp[:, mode_idx], color="#1f3c88", lw=1.8)
        for marker_label, marker_freq in markers.items():
            ax.axvline(marker_freq, color="#7f8c8d", lw=0.9, ls="--", alpha=0.7)
            if marker_freq <= max_freq:
                ax.text(marker_freq, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0 else 0.0, marker_label,
                        color="#7f8c8d", fontsize=7, ha="center", va="top")
        ax.axvline(frequencies_hz[mode_idx], color="#d35400", lw=1.0, ls="-.")
        ax.set_ylabel("|Q(f)|")
        ax.set_xlim(0.0, max_freq)
        ax.grid(alpha=0.25)
        ax.set_title(
            f"yaw=0° — {mode_labels[mode_idx]} ({mode_characters[mode_idx]})  "
            f"basis={frequencies_hz[mode_idx]:.3f} Hz"
        )

    axes[-1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_peak_tracking(
    peaks_df: pd.DataFrame,
    basis_info: dict[str, object],
    plot_modes: int,
    output_png: Path,
    output_pdf: Path,
) -> None:
    frequencies_hz = np.asarray(basis_info["frequencies_hz"], dtype=np.float64)
    mode_labels = list(basis_info["mode_labels"])
    mode_characters = list(basis_info["mode_characters"])

    n_plot = min(plot_modes, len(mode_labels))
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    palette = plt.cm.tab10(np.linspace(0, 1, n_plot))

    for mode_idx in range(n_plot):
        subset = peaks_df[peaks_df["mode_index"] == (mode_idx + 1)].sort_values("yaw_deg")
        color = palette[mode_idx]
        label = f"{mode_labels[mode_idx]} ({mode_characters[mode_idx]})"
        ax.plot(
            subset["yaw_deg"],
            subset["dominant_peak_hz"],
            marker="o",
            color=color,
            lw=1.8,
            label=label,
        )
        ax.axhline(frequencies_hz[mode_idx], color=color, lw=1.0, ls="--", alpha=0.55)

    ax.set_xlabel("Yaw [deg]")
    ax.set_ylabel("Dominant modal-coordinate peak [Hz]")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    ax.set_title("Dominant FFT peak of projected modal coordinates across yaw")
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    summary_path: Path,
    basis_info: dict[str, object],
    args: argparse.Namespace,
    sample_counts: dict[int, int],
    peaks_df: pd.DataFrame,
    markers_df: pd.DataFrame,
) -> None:
    frequencies_hz = np.asarray(basis_info["frequencies_hz"], dtype=np.float64)
    direction_labels = list(basis_info["direction_labels"])
    factors = np.asarray(basis_info["participation_factors"], dtype=np.float64)
    mode_labels = list(basis_info["mode_labels"])
    mode_characters = list(basis_info["mode_characters"])
    one_p_hz = float(basis_info["omega_ref_rpm"]) / 60.0

    lines = [
        f"# Modal projection FFT — {basis_info['basis_mode']}",
        "",
        "This analysis uses a **fixed reduced modal basis** built once from the exact campaign model,",
        "then projects checkpointed `u_red(t)` onto that basis before computing FFTs.",
        "It is a modal-family attribution tool, not a per-checkpoint loaded modal identification.",
        "",
        f"Results root: `{args.results_root}`",
        f"Config: `{Path(basis_info['config_path'])}`",
        f"Window: `{args.start_time:.3f} ≤ t ≤ {args.end_time:.3f}` s",
        f"Basis mode: `{basis_info['basis_mode']}`",
        f"Reference rotor speed: `{basis_info['omega_ref_rpm']:.5f} RPM` → `1P = {one_p_hz:.5f} Hz`",
        "",
        "## Basis frequencies and dominant directions",
        "",
        "| Mode | Character | f [Hz] | " + " | ".join(direction_labels) + " |",
        "|---|---|---:|" + "---:|" * len(direction_labels),
    ]

    for idx, label in enumerate(mode_labels):
        row = [label, mode_characters[idx], f"{frequencies_hz[idx]:.4f}"]
        row.extend(f"{factors[idx, j]:.2f}" for j in range(len(direction_labels)))
        lines.append("| " + " | ".join(row) + " |")

    lines += ["", "## Checkpoint counts", "", "| Yaw [deg] | Samples |", "|---:|---:|"]
    for yaw_deg in sorted(sample_counts):
        lines.append(f"| {yaw_deg} | {sample_counts[yaw_deg]} |")

    lines += [
        "",
        "## Dominant FFT peak per modal coordinate",
        "",
        "| Yaw [deg] | Mode | Character | Basis f [Hz] | Dominant peak [Hz] | f/1P | Nearest P | Peak rel. amp. |",
        "|---:|---|---|---:|---:|---:|---|---:|",
    ]
    for _, row in peaks_df.sort_values(["yaw_deg", "mode_index"]).iterrows():
        lines.append(
            f"| {int(row['yaw_deg'])} | {row['mode_label']} | {row['mode_character']} | "
            f"{row['basis_frequency_hz']:.4f} | {row['dominant_peak_hz']:.4f} | "
            f"{row['dominant_peak_hz'] / one_p_hz:.2f} | {row['nearest_p']} | {row['dominant_peak_rel_amp']:.3f} |"
        )

    marker_names = ["1P", "2P", "3P", "4P", "basis"]
    lines += [
        "",
        "## Marker amplitudes (relative to the dominant FFT amplitude of each modal coordinate)",
        "",
        "| Yaw [deg] | Mode | " + " | ".join(marker_names) + " |",
        "|---:|---|" + "---:|" * len(marker_names),
    ]
    for _, row in markers_df.sort_values(["yaw_deg", "mode_index"]).iterrows():
        values = [f"{row[f'marker_{name.lower()}_rel_amp']:.3f}" for name in marker_names]
        lines.append(f"| {int(row['yaw_deg'])} | {row['mode_label']} | " + " | ".join(values) + " |")

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.start_time >= args.end_time:
        raise RuntimeError("start-time must be smaller than end-time")
    if args.checkpoint_stride < 1:
        raise RuntimeError("checkpoint-stride must be >= 1")

    tag = args.basis_mode
    out = make_output_paths(tag)
    basis_info = build_modal_basis(args.config.resolve(), args.num_modes, args.basis_mode)

    one_p_hz = float(basis_info["omega_ref_rpm"]) / 60.0
    marker_targets = {
        "1p": one_p_hz,
        "2p": 2.0 * one_p_hz,
        "3p": 3.0 * one_p_hz,
        "4p": 4.0 * one_p_hz,
    }

    projection_operator = basis_info["projection_operator"]
    basis_frequencies = np.asarray(basis_info["frequencies_hz"], dtype=np.float64)
    mode_labels = list(basis_info["mode_labels"])
    mode_characters = list(basis_info["mode_characters"])

    peaks_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []
    timeseries_frames: list[pd.DataFrame] = []
    sample_counts: dict[int, int] = {}
    freq_map: dict[int, np.ndarray] = {}
    amp_map: dict[int, np.ndarray] = {}

    for yaw_deg in args.yaws:
        case_dir = args.results_root / f"yaw_{yaw_deg}" / "corotational"
        selected_dirs = select_time_dirs(
            case_dir,
            args.start_time,
            args.end_time,
            args.checkpoint_stride,
            args.max_checkpoints,
        )
        sample_counts[yaw_deg] = len(selected_dirs)

        mean_disp = compute_mean_displacement(selected_dirs)
        times, q_series = project_case_timeseries(selected_dirs, mean_disp, projection_operator)

        if times.size < 2:
            raise RuntimeError(f"Need at least two checkpoints for yaw={yaw_deg}")
        dt = float(np.median(np.diff(times)))

        timeseries_df = pd.DataFrame({"time_s": times, "yaw_deg": yaw_deg})
        for mode_idx, label in enumerate(mode_labels):
            timeseries_df[label] = q_series[:, mode_idx]
        timeseries_frames.append(timeseries_df)

        mode_freq = None
        mode_amp = None
        for mode_idx, label in enumerate(mode_labels):
            freq, amp = onesided_fft(q_series[:, mode_idx], dt)
            if mode_freq is None:
                mode_freq = freq
                mode_amp = np.zeros((freq.size, len(mode_labels)), dtype=np.float64)
            mode_amp[:, mode_idx] = amp

            peak_freq_hz, peak_amp = dominant_peak(freq, amp)
            peak_rel_amp = peak_amp / float(np.max(amp)) if np.max(amp) > 0.0 else 0.0
            peaks_rows.append(
                {
                    "yaw_deg": yaw_deg,
                    "mode_index": mode_idx + 1,
                    "mode_label": label,
                    "mode_character": mode_characters[mode_idx],
                    "basis_frequency_hz": basis_frequencies[mode_idx],
                    "dominant_peak_hz": peak_freq_hz,
                    "dominant_peak_amp": peak_amp,
                    "dominant_peak_rel_amp": peak_rel_amp,
                    "nearest_p": nearest_p_label(peak_freq_hz, one_p_hz),
                }
            )

            amp_max = float(np.max(amp)) if np.max(amp) > 0.0 else 1.0
            marker_row = {
                "yaw_deg": yaw_deg,
                "mode_index": mode_idx + 1,
                "mode_label": label,
                "mode_character": mode_characters[mode_idx],
                "basis_frequency_hz": basis_frequencies[mode_idx],
            }
            for marker_name, marker_freq in marker_targets.items():
                marker_row[f"marker_{marker_name}_freq_hz"] = marker_freq
                marker_row[f"marker_{marker_name}_amp"] = nearest_bin_amplitude(freq, amp, marker_freq)
                marker_row[f"marker_{marker_name}_rel_amp"] = marker_row[f"marker_{marker_name}_amp"] / amp_max
            marker_row["marker_basis_freq_hz"] = basis_frequencies[mode_idx]
            marker_row["marker_basis_amp"] = nearest_bin_amplitude(freq, amp, basis_frequencies[mode_idx])
            marker_row["marker_basis_rel_amp"] = marker_row["marker_basis_amp"] / amp_max
            marker_rows.append(marker_row)

        assert mode_freq is not None and mode_amp is not None
        freq_map[yaw_deg] = mode_freq
        amp_map[yaw_deg] = mode_amp

    peaks_df = pd.DataFrame.from_records(peaks_rows)
    markers_df = pd.DataFrame.from_records(marker_rows)
    timeseries_df = pd.concat(timeseries_frames, ignore_index=True)

    timeseries_df.to_csv(out["timeseries_csv"], index=False)
    peaks_df.to_csv(out["peaks_csv"], index=False)
    markers_df.to_csv(out["markers_csv"], index=False)
    write_summary(out["summary_md"], basis_info, args, sample_counts, peaks_df, markers_df)

    if not args.skip_plots:
        plot_yaw0_fft(
            freq_map,
            amp_map,
            basis_info,
            args.max_freq,
            args.plot_modes,
            out["yaw0_png"],
            out["yaw0_pdf"],
        )
        plot_peak_tracking(
            peaks_df,
            basis_info,
            args.plot_modes,
            out["tracking_png"],
            out["tracking_pdf"],
        )

    print(f"Wrote {out['timeseries_csv']}")
    print(f"Wrote {out['peaks_csv']}")
    print(f"Wrote {out['markers_csv']}")
    print(f"Wrote {out['summary_md']}")
    if not args.skip_plots:
        print(f"Wrote {out['yaw0_png']}")
        print(f"Wrote {out['tracking_png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())