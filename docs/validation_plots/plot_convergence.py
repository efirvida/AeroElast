"""Plot preCICE convergence diagnostics for the corotational yaw campaign."""

from __future__ import annotations

import pathlib

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
DOCS_DIR = pathlib.Path(__file__).resolve().parents[1]
OUT_DIR = DOCS_DIR / "figures"
SUMMARY_CSV = DOCS_DIR / "validation_data" / "generated" / "corotational_convergence_summary.csv"
COROT_ROOT = BASE / "frontiersin_results_corotational_100s"

YAWS = [0, 10, 20, 30, 40]
T_START = 40.0
T_END = 100.0
DT = 0.01
ROLLING_WINDOWS = 100
RESIDUAL_THRESHOLD = 1e-4

COLORS = {
    0: "#c0392b",
    10: "#d35400",
    20: "#16a085",
    30: "#2980b9",
    40: "#8e44ad",
    "threshold": "#2c3e50",
    "disp": "#8e44ad",
    "force": "#16a085",
}

OUT_DIR.mkdir(parents=True, exist_ok=True)


def _window_time(time_window: int) -> float:
    return time_window * DT


def _in_window(time_window: int) -> bool:
    time_s = _window_time(time_window)
    return T_START < time_s <= T_END


def parse_iterations(path: pathlib.Path, yaw: int) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            if "TimeWindow" in line:
                break
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                time_window = int(parts[0])
                iterations = int(parts[2])
                convergence = int(parts[3])
            except ValueError:
                continue
            if _in_window(time_window):
                rows.append({
                    "yaw_deg": yaw,
                    "time_window": time_window,
                    "time_s": _window_time(time_window),
                    "iterations": iterations,
                    "convergence": convergence,
                })
    return pd.DataFrame(rows)


def parse_final_residuals(path: pathlib.Path, yaw: int) -> pd.DataFrame:
    rows = []
    with path.open() as handle:
        for line in handle:
            if "TimeWindow" in line and "Iteration" in line:
                break
        for line in handle:
            parts = line.split()
            if len(parts) != 4:
                continue
            try:
                time_window = int(parts[0])
                iteration = int(parts[1])
                displacement = float(parts[2])
                force = float(parts[3])
            except ValueError:
                continue
            if _in_window(time_window):
                rows.append({
                    "yaw_deg": yaw,
                    "time_window": time_window,
                    "time_s": _window_time(time_window),
                    "iteration": iteration,
                    "displacement": displacement,
                    "force": force,
                })
    if not rows:
        return pd.DataFrame(columns=["yaw_deg", "time_window", "time_s", "displacement", "force"])
    df = pd.DataFrame(rows)
    idx = df.groupby(["yaw_deg", "time_window"])["iteration"].idxmax()
    return df.loc[idx, ["yaw_deg", "time_window", "time_s", "displacement", "force"]].reset_index(
        drop=True
    )


def load_iteration_data() -> pd.DataFrame:
    frames = []
    for yaw in YAWS:
        frames.append(
            parse_iterations(COROT_ROOT / f"yaw_{yaw}" / "precice-Solid-iterations.log", yaw)
        )
    return pd.concat(frames, ignore_index=True)


def load_residual_data() -> pd.DataFrame:
    frames = []
    for yaw in YAWS:
        frames.append(
            parse_final_residuals(COROT_ROOT / f"yaw_{yaw}" / "precice-Solid-convergence.log", yaw)
        )
    return pd.concat(frames, ignore_index=True)


def build_summary(iterations: pd.DataFrame, residuals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for yaw in YAWS:
        iter_case = iterations[iterations["yaw_deg"] == yaw]
        res_case = residuals[residuals["yaw_deg"] == yaw]
        if iter_case.empty or res_case.empty:
            continue
        rows.append({
            "yaw_deg": yaw,
            "n_windows": int(len(iter_case)),
            "t_start_s": T_START,
            "t_end_s": T_END,
            "iter_mean": round(float(iter_case["iterations"].mean()), 3),
            "iter_std": round(float(iter_case["iterations"].std(ddof=1)), 3),
            "iter_max": int(iter_case["iterations"].max()),
            "iter_min": int(iter_case["iterations"].min()),
            "iter_total": int(iter_case["iterations"].sum()),
            "no_conv_windows": int((iter_case["convergence"] == 0).sum()),
            "res_ResRelSolid-Mesh_Displacement_mean_final": float(
                f"{res_case['displacement'].mean():.3e}"
            ),
            "res_ResRelSolid-Mesh_Displacement_max_final": float(
                f"{res_case['displacement'].max():.3e}"
            ),
            "res_ResRelSolid-Mesh_Displacement_above_1e-4": int(
                (res_case["displacement"] > RESIDUAL_THRESHOLD).sum()
            ),
            "res_ResRelSolid-Mesh_Force_mean_final": float(f"{res_case['force'].mean():.3e}"),
            "res_ResRelSolid-Mesh_Force_max_final": float(f"{res_case['force'].max():.3e}"),
            "res_ResRelSolid-Mesh_Force_above_1e-4": int(
                (res_case["force"] > RESIDUAL_THRESHOLD).sum()
            ),
        })
    return pd.DataFrame(rows)


def save_figure(fig: plt.Figure, basename: str) -> None:
    for extension in ("png", "pdf"):
        output = OUT_DIR / f"{basename}.{extension}"
        fig.savefig(output, bbox_inches="tight")
        print(f"Saved: {output}")
    plt.close(fig)


def plot_iterations(iterations: pd.DataFrame, summary: pd.DataFrame) -> None:
    fig, (ax_summary, ax_histogram) = plt.subplots(1, 2, figsize=(10.4, 4.0))

    ax_summary.errorbar(
        summary["yaw_deg"],
        summary["iter_mean"],
        yerr=summary["iter_std"],
        marker="o",
        color="#c0392b",
        linewidth=1.5,
        capsize=3,
        label="media +/- std",
    )
    ax_summary.plot(
        summary["yaw_deg"],
        summary["iter_max"],
        "s--",
        color="#2c3e50",
        linewidth=1.2,
        label="maximo observado",
    )
    ax_summary.set_xlabel("Yaw [deg]")
    ax_summary.set_ylabel("Iteraciones por ventana")
    ax_summary.set_title("Nivel y extremos de convergencia")
    ax_summary.set_xticks(YAWS)
    ax_summary.set_ylim(0, float(summary["iter_max"].max()) + 2)
    ax_summary.grid(alpha=0.3)
    ax_summary.legend(frameon=False, loc="upper right")

    categories = ["2", "3", "4", ">=5"]
    category_colors = ["#d5dbdb", "#7fb3d5", "#2874a6", "#922b21"]
    bottoms = np.zeros(len(YAWS))
    for category, color in zip(categories, category_colors):
        values = []
        for yaw in YAWS:
            case = iterations[iterations["yaw_deg"] == yaw]
            if category == ">=5":
                count = int((case["iterations"] >= 5).sum())
            else:
                count = int((case["iterations"] == int(category)).sum())
            values.append(100.0 * count / len(case))
        ax_histogram.bar(YAWS, values, bottom=bottoms, width=7.0, color=color, label=category)
        bottoms += np.array(values)

    for _, row in summary.iterrows():
        ax_histogram.text(
            row["yaw_deg"],
            102.0,
            f"max {int(row['iter_max'])}",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    ax_histogram.set_xlabel("Yaw [deg]")
    ax_histogram.set_ylabel("Ventanas [%]")
    ax_histogram.set_title("Distribucion de subiteraciones")
    ax_histogram.set_xticks(YAWS)
    ax_histogram.set_ylim(0, 110)
    ax_histogram.grid(axis="y", alpha=0.3)
    ax_histogram.legend(title="Iter.", frameon=False, loc="center right")

    fig.suptitle("Convergencia iterativa preCICE - campana corrotacional", y=1.02, fontsize=10)
    fig.tight_layout()
    save_figure(fig, "fig_5_4_5a_precice_iterations_corotational")


def plot_residuals(residuals: pd.DataFrame) -> None:
    fig, (ax_disp, ax_force) = plt.subplots(1, 2, figsize=(10.4, 4.0), sharey=True)

    for ax, column, color, title in [
        (ax_disp, "displacement", COLORS["disp"], "Residual de desplazamiento"),
        (ax_force, "force", COLORS["force"], "Residual de fuerza"),
    ]:
        values = [residuals.loc[residuals["yaw_deg"] == yaw, column].to_numpy() for yaw in YAWS]
        box = ax.boxplot(
            values,
            positions=YAWS,
            widths=6.0,
            patch_artist=True,
            showfliers=False,
            whis=(5, 95),
        )
        for patch in box["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.32)
            patch.set_edgecolor(color)
        for key in ("whiskers", "caps", "medians"):
            for artist in box[key]:
                artist.set_color(color)
                artist.set_linewidth(1.2)
        max_values = [float(np.max(v)) for v in values]
        ax.plot(YAWS, max_values, "o--", color="#2c3e50", linewidth=1.0, label="maximo")
        ax.axhline(
            1e-4, color=COLORS["threshold"], linewidth=1.0, linestyle=":", label="umbral 1e-4"
        )
        ax.set_xlabel("Yaw [deg]")
        ax.set_title(title)
        ax.set_xticks(YAWS)
        ax.set_yscale("log")
        ax.grid(alpha=0.3, which="both")
        ax.legend(frameon=False, loc="lower left")

    ax_disp.set_ylabel("Residual relativo final por ventana")
    fig.suptitle("Residuales finales preCICE - campana corrotacional", y=1.02, fontsize=10)
    fig.tight_layout()
    save_figure(fig, "fig_5_4_5b_precice_residuals_corotational")


def main() -> None:
    iterations = load_iteration_data()
    residuals = load_residual_data()
    summary = build_summary(iterations, residuals).sort_values("yaw_deg")
    summary.to_csv(SUMMARY_CSV, index=False)
    print(f"Saved: {SUMMARY_CSV}")
    plot_iterations(iterations, summary)
    plot_residuals(residuals)


if __name__ == "__main__":
    main()
