#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from yaw_sweep_convergence import (
    analyze_iterations,
    analyze_residuals,
    filter_windows,
    parse_convergence_log,
    parse_iterations_log,
)

T_START = 20.0
T_END = 70.0


def _pick_residual_field(stats: dict, token: str) -> dict:
    for key, value in stats.items():
        if token in key:
            return value
    return {}


def read_bem_report(path: Path, t_start: float = T_START, t_end: float = T_END) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[(df["Time [s]"] >= t_start) & (df["Time [s]"] <= t_end)].reset_index(drop=True)


def read_structural_report(
    path: Path, t_start: float = T_START, t_end: float = T_END
) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df[(df["Time [s]"] >= t_start) & (df["Time [s]"] <= t_end)].reset_index(drop=True)


def summarize_bem(df: pd.DataFrame) -> dict:
    return {
        "n_samples": int(len(df)),
        "t_end": float(df["Time [s]"].iloc[-1]),
        "flap_mean": float(df["Tip Disp Y [m]"].mean()),
        "flap_std": float(df["Tip Disp Y [m]"].std()),
        "flap_max": float(df["Tip Disp Y [m]"].max()),
        "flap_p2p": float(df["Tip Disp Y [m]"].max() - df["Tip Disp Y [m]"].min()),
        "edge_mean": float(df["Tip Disp X [m]"].mean()),
        "power_mean_mw": float(df["Power [W]"].mean() / 1e6),
        "power_std_mw": float(df["Power [W]"].std() / 1e6),
        "power_p2p_mw": float((df["Power [W]"].max() - df["Power [W]"].min()) / 1e6),
        "power_max_mw": float(df["Power [W]"].max() / 1e6),
        "cp_mean": float(df["CP"].mean()),
        "ct_mean": float(df["CT"].mean()),
        "thrust_mean_mn": float(df["Thrust [N]"].mean() / 1e6),
        "torque_mean_mnm": float(df["Torque [N.m]"].mean() / 1e6),
    }


def summarize_structural_flap(df: pd.DataFrame) -> dict:
    flap = df["Max Disp Y [m]"]
    return {
        "n_samples": int(len(df)),
        "t_end": float(df["Time [s]"].iloc[-1]),
        "flap_mean": float(flap.mean()),
        "flap_std": float(flap.std()),
        "flap_p2p": float(flap.max() - flap.min()),
    }


def summarize_corot_campaign(root: Path) -> list[dict]:
    rows = []
    for yaw in (0, 10, 20, 30, 40):
        df = read_bem_report(root / f"yaw_{yaw}" / "fluid" / "bem_report.csv")
        stats = summarize_bem(df)
        rows.append({"yaw_deg": yaw, **stats})
    return rows


def summarize_inertial_yaw(root: Path) -> list[dict]:
    rows = []
    for yaw in (0, 10, 20, 30, 40):
        bem = read_bem_report(root / f"yaw_{yaw}" / "fluid" / "bem_report.csv")
        bem_stats = summarize_bem(bem)
        rows.append({
            "yaw_deg": yaw,
            "source": "bem_report.csv",
            "note": "completed campaign, common 20-70 s window",
            **{
                k: bem_stats[k] for k in ("n_samples", "t_end", "flap_mean", "flap_std", "flap_p2p")
            },
        })
    return rows


def summarize_yaw0_pair(corot_root: Path, inert_root: Path) -> dict:
    corot = summarize_bem(read_bem_report(corot_root / "yaw_0" / "fluid" / "bem_report.csv"))
    inert = summarize_bem(read_bem_report(inert_root / "yaw_0" / "fluid" / "bem_report.csv"))
    return {"corot": corot, "inert": inert}


def summarize_convergence(root: Path, dt: float = 0.01) -> list[dict]:
    rows = []
    for yaw in (0, 10, 20, 30, 40):
        yaw_dir = root / f"yaw_{yaw}"
        iterations = parse_iterations_log(yaw_dir / "precice-Solid-iterations.log")
        convergence = parse_convergence_log(yaw_dir / "precice-Solid-convergence.log")
        iterations = filter_windows(iterations, T_START, T_END, dt)
        convergence = filter_windows(convergence, T_START, T_END, dt)
        iter_stats = analyze_iterations(iterations)
        res_stats = analyze_residuals(convergence)
        disp_stats = _pick_residual_field(res_stats, "Displacement")
        force_stats = _pick_residual_field(res_stats, "Force")
        rows.append({
            "yaw_deg": yaw,
            "n_windows": int(iter_stats["n_windows"]),
            "iter_mean": float(iter_stats["iter_mean"]),
            "iter_std": float(iter_stats["iter_std"]),
            "iter_min": int(iter_stats["iter_min"]),
            "iter_max": int(iter_stats["iter_max"]),
            "iter_total": int(iter_stats["iter_total"]),
            "no_conv_windows": int(iter_stats["no_conv_windows"]),
            "disp_res_mean": float(disp_stats.get("mean_final", float("nan"))),
            "disp_res_max": float(disp_stats.get("max_final", float("nan"))),
            "force_res_mean": float(force_stats.get("mean_final", float("nan"))),
            "force_res_max": float(force_stats.get("max_final", float("nan"))),
        })
    return rows


def build_payload(corot_root: Path, inert_root: Path) -> dict:
    return {
        "corot_yaw_sweep": summarize_corot_campaign(corot_root),
        "inertial_yaw_sweep": summarize_inertial_yaw(inert_root),
        "yaw0_pair": summarize_yaw0_pair(corot_root, inert_root),
        "convergence": summarize_convergence(corot_root),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rebuild V-04 validation data from current results."
    )
    parser.add_argument(
        "--corot-root",
        type=Path,
        default=Path("/scratch/leahk/eduardo.donestevez/frontiersin_results_corotational"),
    )
    parser.add_argument(
        "--inert-root",
        type=Path,
        default=Path("/scratch/leahk/eduardo.donestevez/frontiersin_results_inertial"),
    )
    parser.add_argument("--indent", type=int, default=2)
    args = parser.parse_args()

    payload = build_payload(args.corot_root, args.inert_root)
    print(json.dumps(payload, indent=args.indent, ensure_ascii=False))


if __name__ == "__main__":
    main()
