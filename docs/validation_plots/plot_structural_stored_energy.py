"""Dynamic elastic and kinetic energy stored in the blade motion.

This postprocess reconstructs the reduced stiffness and lumped mass matrices
offline from the Frontiers solid YAML, then evaluates dynamic stored energy from
checkpointed reduced states over the final three revolutions of the
corotational yaw sweep.

Definitions
-----------
Let ``u'(t) = u(t) - <u>`` and ``v'(t) = v(t) - <v>`` over the analysis window.
The dynamic stored energies are evaluated as

    U_el(t) = 0.5 * u'(t)^T K_red u'(t)
    T_kin(t) = 0.5 * v'(t)^T M_red v'(t)

where ``K_red`` is the reduced elastic stiffness matrix after Dirichlet BCs and
``M_red`` is the reduced lumped mass matrix. This is a dynamic-storage metric,
not a full rotating-frame energy balance: it excludes spin-softening and uses
the linear elastic stiffness as the material strain-energy measure.

Runtime requirements
--------------------
On SDumont this script needs the same runtime libraries as the compiled
``_aeroelast`` module. Use the cluster pattern already adopted by the repo:

    module purge
    module load glu gcc/14.2.0_sequana
    export LD_LIBRARY_PATH=${GCC14_LIB}:$VENV_DIR/lib:$VENV_DIR/lib64:/usr/lib64/psm2-compat:/scratch/app/openmpi/4.1.4_gnu/lib

Outputs
-------
docs/figures/fig_5_4_7n_dynamic_structural_energy.png
docs/figures/fig_5_4_7n_dynamic_structural_energy.pdf
docs/figures/fig_5_4_7o_dynamic_energy_summary.png
docs/figures/fig_5_4_7o_dynamic_energy_summary.pdf
docs/validation_data/generated/structural_stored_energy_timeseries.csv
docs/validation_data/generated/structural_stored_energy_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

BASE = Path("/scratch/leahk/eduardo.donestevez")
ROOT = BASE / "frontiersin_results_corotational_100s"
OUT_DIR = Path(__file__).resolve().parent.parent / "figures"
DATA_DIR = Path(__file__).resolve().parent.parent / "validation_data" / "generated"
OUT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

YAWS = [0, 10, 20, 30, 40]
TIME_SERIES_YAWS = [0, 20, 40]
OMEGA_RPM = 7.55
ROTOR_PERIOD_S = 60.0 / OMEGA_RPM
LAST_REVOLUTIONS = 3
CONFIG_PATH = Path("tests/IEA15MW/frontiersin_2025_yaw/solid_corotational_yaw_0.yaml").resolve()

COLORS = {
    "elastic": "#1f3c88",
    "kinetic": "#b7410e",
    "total": "#287271",
}


def import_runner() -> type:
    sys.path.insert(0, str(Path("src").resolve()))
    try:
        from aeroelast.solvers.fsi.runner import FSIRunner
    except ImportError as exc:  # pragma: no cover - runtime environment specific
        raise RuntimeError(
            "Could not import aeroelast runtime. On SDumont, load GCC 14 and add "
            "OpenMPI + venv libs to LD_LIBRARY_PATH before running this script."
        ) from exc
    return FSIRunner


def assemble_reduced_matrices() -> tuple[sp.csr_matrix, np.ndarray, float, float]:
    FSIRunner = import_runner()
    runner = FSIRunner(CONFIG_PATH, working_dir=CONFIG_PATH.parent)
    runner.mesh = runner._setup_mesh()
    runner._material = runner._create_material()
    model_config = runner._build_model_config()
    runner.solver = runner._create_solver(model_config)
    runner._apply_boundary_conditions()

    solver = runner.solver
    omega0 = float(solver.solver_params["rotor"]["omega"])
    (stiffness, mass), bc_manager, _ = solver._assemble_system_matrices(omega0)
    k_red = bc_manager.reduce_matrix(stiffness)
    m_red = bc_manager.reduce_matrix(mass)

    indptr, indices, data = k_red.getValuesCSR()
    k_csr = sp.csr_matrix((data, indices, indptr), shape=k_red.getSize())
    m_diag = m_red.getDiagonal().getArray(readonly=True).copy()

    return k_csr, m_diag, float(solver._eta_k), float(solver._eta_m)


def list_time_dirs(case_dir: Path) -> list[tuple[float, Path]]:
    time_dirs: list[tuple[float, Path]] = []
    for child in case_dir.iterdir():
        if not child.is_dir():
            continue
        try:
            time_value = float(child.name)
        except ValueError:
            continue
        state_file = child / "state.npz"
        if state_file.exists():
            time_dirs.append((time_value, child))
    time_dirs.sort(key=lambda item: item[0])
    return time_dirs


def load_case_states(case_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    time_dirs = list_time_dirs(case_dir)
    if not time_dirs:
        raise RuntimeError(f"No checkpoint states found in {case_dir}")

    end_time = time_dirs[-1][0]
    start_time = max(0.0, end_time - LAST_REVOLUTIONS * ROTOR_PERIOD_S)

    times: list[float] = []
    disp: list[np.ndarray] = []
    vel: list[np.ndarray] = []
    for time_value, folder in time_dirs:
        if time_value < start_time:
            continue
        state = np.load(folder / "state.npz")
        times.append(time_value)
        disp.append(np.asarray(state["u_red"], dtype=float))
        vel.append(np.asarray(state["v_red"], dtype=float))

    return np.asarray(times, dtype=float), np.vstack(disp), np.vstack(vel)


def last_revolution_windows(end_time: float) -> list[tuple[int, float, float]]:
    windows: list[tuple[int, float, float]] = []
    for rev_index in range(LAST_REVOLUTIONS, 0, -1):
        start = end_time - rev_index * ROTOR_PERIOD_S
        end = end_time - (rev_index - 1) * ROTOR_PERIOD_S
        windows.append((LAST_REVOLUTIONS - rev_index + 1, start, end))
    return windows


def compute_case_energy(
    k_csr: sp.csr_matrix,
    m_diag: np.ndarray,
    times: np.ndarray,
    u_red: np.ndarray,
    v_red: np.ndarray,
) -> pd.DataFrame:
    u_dyn = u_red - u_red.mean(axis=0)
    v_dyn = v_red - v_red.mean(axis=0)
    elastic = 0.5 * np.sum(u_dyn.T * (k_csr @ u_dyn.T), axis=0)
    kinetic = 0.5 * np.sum((v_dyn**2) * m_diag[None, :], axis=1)
    total = elastic + kinetic
    return pd.DataFrame({
        "time_s": times,
        "elastic_energy_kJ": elastic / 1.0e3,
        "kinetic_energy_kJ": kinetic / 1.0e3,
        "total_dynamic_energy_kJ": total / 1.0e3,
    })


def summarize_case(
    yaw_deg: int, energy_df: pd.DataFrame, aero_energy_mj: float, into_structure_kj: float
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    end_time = float(energy_df["time_s"].max())
    for window_id, start, end in last_revolution_windows(end_time):
        window = energy_df[(energy_df["time_s"] >= start) & (energy_df["time_s"] <= end)].copy()
        if len(window) < 2:
            continue
        row = {
            "yaw_deg": yaw_deg,
            "window_id": window_id,
            "elastic_energy_mean_kJ": float(window["elastic_energy_kJ"].mean()),
            "elastic_energy_max_kJ": float(window["elastic_energy_kJ"].max()),
            "kinetic_energy_mean_kJ": float(window["kinetic_energy_kJ"].mean()),
            "kinetic_energy_max_kJ": float(window["kinetic_energy_kJ"].max()),
            "total_energy_mean_kJ": float(window["total_dynamic_energy_kJ"].mean()),
            "total_energy_max_kJ": float(window["total_dynamic_energy_kJ"].max()),
            "elastic_max_pct_of_aero": 100.0
            * float(window["elastic_energy_kJ"].max())
            / (aero_energy_mj * 1.0e3),
            "total_max_pct_of_aero": 100.0
            * float(window["total_dynamic_energy_kJ"].max())
            / (aero_energy_mj * 1.0e3),
            "elastic_max_pct_of_into_structure": 100.0
            * float(window["elastic_energy_kJ"].max())
            / into_structure_kj,
        }
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def plot_time_series(timeseries_map: dict[int, pd.DataFrame]) -> None:
    fig, axes = plt.subplots(len(TIME_SERIES_YAWS), 1, figsize=(10.2, 8.8), sharex=False)
    for ax, yaw_deg in zip(axes, TIME_SERIES_YAWS):
        df = timeseries_map[yaw_deg].copy()
        end_time = float(df["time_s"].max())
        start_time = end_time - ROTOR_PERIOD_S
        final_rev = df[df["time_s"] >= start_time].copy()
        time_local = final_rev["time_s"].to_numpy(dtype=float) - start_time

        ax.plot(
            time_local,
            final_rev["elastic_energy_kJ"],
            color=COLORS["elastic"],
            lw=2.0,
            label="Energia elastica dinamica",
        )
        ax.plot(
            time_local,
            final_rev["kinetic_energy_kJ"],
            color=COLORS["kinetic"],
            lw=2.0,
            label="Energia cinetica dinamica",
        )
        ax.plot(
            time_local,
            final_rev["total_dynamic_energy_kJ"],
            color=COLORS["total"],
            lw=2.0,
            label="Energia almacenada total",
        )
        ax.set_ylabel("Energia [kJ]")
        ax.set_title(f"Yaw = {yaw_deg} deg")
        ax.grid(alpha=0.25)
        ax.legend(loc="upper right", frameon=False)

    axes[-1].set_xlabel("Tiempo dentro de la ultima vuelta [s]")
    fig.suptitle("Energia dinamica almacenada en la pala durante la ultima vuelta")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_5_4_7n_dynamic_structural_energy.png", dpi=220)
    fig.savefig(OUT_DIR / "fig_5_4_7n_dynamic_structural_energy.pdf")
    plt.close(fig)


def plot_summary(summary_df: pd.DataFrame) -> None:
    grouped = summary_df.groupby("yaw_deg", as_index=False).mean(numeric_only=True)
    yaw = grouped["yaw_deg"].to_numpy(dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 8.2), sharex=True)

    width = 3.0
    axes[0].bar(
        yaw - width / 2.0,
        grouped["elastic_max_pct_of_aero"],
        width=width,
        color=COLORS["elastic"],
        label="Pico elastico / energia aero",
    )
    axes[0].bar(
        yaw + width / 2.0,
        grouped["total_max_pct_of_aero"],
        width=width,
        color=COLORS["total"],
        label="Pico total / energia aero",
    )
    axes[0].set_ylabel("Pico de energia almacenada [% de energia aero por vuelta]")
    axes[0].set_title(
        "Escala de energia almacenada dinamica frente a la energia aerodinamica extraida"
    )
    axes[0].grid(alpha=0.25, axis="y")
    axes[0].legend(frameon=False)

    axes[1].plot(
        yaw, grouped["elastic_max_pct_of_into_structure"], color="#222222", marker="o", lw=2.0
    )
    axes[1].set_xlabel("Angulo de yaw [deg]")
    axes[1].set_ylabel("Pico elastico / energia que entra a la estructura [%]")
    axes[1].grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_5_4_7o_dynamic_energy_summary.png", dpi=220)
    fig.savefig(OUT_DIR / "fig_5_4_7o_dynamic_energy_summary.pdf")
    plt.close(fig)


def main() -> None:
    k_csr, m_diag, eta_k, eta_m = assemble_reduced_matrices()
    channel_summary = pd.read_csv(DATA_DIR / "structural_power_exchange_summary.csv").set_index(
        "yaw_deg"
    )

    timeseries_map: dict[int, pd.DataFrame] = {}
    timeseries_rows: list[pd.DataFrame] = []
    summary_rows: list[pd.DataFrame] = []

    for yaw_deg in YAWS:
        case_dir = ROOT / f"yaw_{yaw_deg}" / "corotational"
        if not case_dir.exists():
            continue
        times, u_red, v_red = load_case_states(case_dir)
        energy_df = compute_case_energy(k_csr, m_diag, times, u_red, v_red)
        energy_df.insert(0, "yaw_deg", yaw_deg)
        timeseries_rows.append(energy_df)
        timeseries_map[yaw_deg] = energy_df

        summary_rows.append(
            summarize_case(
                yaw_deg=yaw_deg,
                energy_df=energy_df,
                aero_energy_mj=float(channel_summary.loc[yaw_deg, "aero_energy_MJ_mean"]),
                into_structure_kj=float(channel_summary.loc[yaw_deg, "into_structure_kJ_mean"]),
            )
        )

    timeseries_df = pd.concat(timeseries_rows, ignore_index=True)
    summary_df = pd.concat(summary_rows, ignore_index=True)
    timeseries_df.to_csv(DATA_DIR / "structural_stored_energy_timeseries.csv", index=False)
    summary_df.to_csv(DATA_DIR / "structural_stored_energy_summary.csv", index=False)

    plot_time_series(timeseries_map)
    plot_summary(summary_df)

    grouped = summary_df.groupby("yaw_deg", as_index=False).mean(numeric_only=True)
    print("DYNAMIC STORED STRUCTURAL ENERGY")
    print(f"Rayleigh coefficients used offline: eta_k={eta_k:.6e} s, eta_m={eta_m:.6e} 1/s")
    print()
    for _, row in grouped.iterrows():
        print(
            f"yaw={int(row['yaw_deg']):>2}: "
            f"Uel_max={row['elastic_energy_max_kJ']:.3f} kJ, "
            f"Tkin_max={row['kinetic_energy_max_kJ']:.3f} kJ, "
            f"Etot_max={row['total_energy_max_kJ']:.3f} kJ, "
            f"Uel_max/Eaero={row['elastic_max_pct_of_aero']:.3f}%, "
            f"Uel_max/Ein={row['elastic_max_pct_of_into_structure']:.2f}%"
        )

    print()
    print(f"Saved: {OUT_DIR / 'fig_5_4_7n_dynamic_structural_energy.png'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7n_dynamic_structural_energy.pdf'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7o_dynamic_energy_summary.png'}")
    print(f"Saved: {OUT_DIR / 'fig_5_4_7o_dynamic_energy_summary.pdf'}")
    print(f"Saved: {DATA_DIR / 'structural_stored_energy_timeseries.csv'}")
    print(f"Saved: {DATA_DIR / 'structural_stored_energy_summary.csv'}")


if __name__ == "__main__":
    main()
