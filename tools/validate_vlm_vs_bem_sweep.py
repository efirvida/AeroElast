"""Run repeatable IEA-15 MW reduced-order VLM versus BEM validation sweeps."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from aeroelast.solvers.aero.validation import (
    IEA15MW_REFERENCE_OPERATING_POINTS,
    IEA15MW_ZERO_PITCH_TSR_SWEEP,
    RotorModelComparison,
    RotorValidationObjectiveMetrics,
    compare_vlm_experimental_gated_lift_to_bem_sweep,
    compare_vlm_experimental_gated_tsr_lift_to_bem_sweep,
    compare_vlm_experimental_lift_to_bem_sweep,
    compare_vlm_polar_corrected_to_bem_sweep,
    compare_vlm_to_bem_sweep,
    comparison_as_dict,
    compute_validation_objective_metrics,
)


def _resolve_default_blade_file(repo_root: Path) -> Path:
    candidates = [
        repo_root / "tests" / "IEA-15-240-RWT.yaml",
        repo_root / "examples" / "reference_turbines" / "yamls" / "IEA-15-240-RWT.yaml",
        repo_root.parent / "simulations" / "blade" / "solid" / "IEA-15-240-RWT.yaml",
        repo_root.parent
        / "simulations"
        / "IEA-15-240-RWT"
        / "fsi"
        / "NoTower"
        / "solid"
        / "IEA-15-240-RWT.yaml",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "IEA-15-240-RWT.yaml not found. Use --blade-file to point at the WindIO blade definition."
    )


def _format_ratio(value: float) -> str:
    return "nan" if value != value else f"{value:6.3f}"


def _print_table(title: str, comparisons: list[RotorModelComparison]) -> None:
    print(title)
    print("case          V [m/s]  rpm   pitch    T_ratio  Q_ratio  CP_vlm  CP_bem  CT_vlm  CT_bem")
    for row in comparisons:
        point = row.operating_point
        print(
            f"{point.name:<12}"
            f"{point.wind_speed:7.3f}"
            f" {point.omega_rpm:6.3f}"
            f"{point.pitch_deg:8.3f}"
            f"{_format_ratio(row.thrust_ratio)}"
            f"{_format_ratio(row.torque_ratio)}"
            f"{row.vlm.cp:8.3f}"
            f"{row.bem.cp:8.3f}"
            f"{row.vlm.ct:8.3f}"
            f"{row.bem.ct:8.3f}"
        )
    print()


def _print_summary(
    model_label: str,
    reference_rows: list[RotorModelComparison],
    tsr_rows: list[RotorModelComparison],
) -> RotorValidationObjectiveMetrics:
    metrics = compute_validation_objective_metrics(reference_rows, tsr_rows)
    thrust_ratios = [row.thrust_ratio for row in reference_rows]
    torque_ratios = [row.torque_ratio for row in reference_rows]
    zero_pitch_rows = [
        row for row in reference_rows if abs(row.operating_point.pitch_deg) < 1.0e-12
    ]
    zero_pitch_torque_ratios = [row.torque_ratio for row in zero_pitch_rows]

    bem_cp_peak = max(tsr_rows, key=lambda row: row.bem.cp)
    vlm_cp_peak = max(tsr_rows, key=lambda row: row.vlm.cp)

    print(f"Summary ({model_label})")
    print(
        "- Reference-curve thrust ratios stay between "
        f"{min(thrust_ratios):.3f} and {max(thrust_ratios):.3f}."
    )
    print(
        "- Reference-curve torque ratios stay between "
        f"{min(torque_ratios):.3f} and {max(torque_ratios):.3f}."
    )
    print(
        "- Zero-pitch torque ratios collapse to "
        f"{min(zero_pitch_torque_ratios):.3f}–{max(zero_pitch_torque_ratios):.3f},"
        " which remains the main power-prediction gap of this reduced-order model.",
    )
    print(
        "- In the zero-pitch TSR sweep, BEM peaks at "
        f"{bem_cp_peak.operating_point.omega_rpm:.3f} rpm with CP={bem_cp_peak.bem.cp:.3f},"
        " while the reduced VLM keeps increasing up to "
        f"{vlm_cp_peak.operating_point.omega_rpm:.3f} rpm with CP={vlm_cp_peak.vlm.cp:.3f}."
    )
    print("- Objective error metrics:")
    print(
        "  MAE(|1-T_ratio|)="
        f"{metrics.mean_abs_thrust_ratio_error:.3f}, "
        "MAE(|1-Q_ratio|)="
        f"{metrics.mean_abs_torque_ratio_error:.3f}, "
        "MAE(|1-P_ratio|)="
        f"{metrics.mean_abs_power_ratio_error:.3f}."
    )
    print(
        "  Zero-pitch MAE(|1-Q_ratio|)="
        f"{metrics.zero_pitch_mean_abs_torque_ratio_error:.3f}, "
        "CP_RMSE="
        f"{metrics.cp_curve_rmse:.3f}."
    )
    print(
        "  CP peak mismatch: Δrpm="
        f"{metrics.cp_peak_rpm_error:.3f}, "
        "ΔCP="
        f"{metrics.cp_peak_value_error:.3f}."
    )
    return metrics


def _print_objective_scoreboard(
    objective_rows: list[tuple[str, RotorValidationObjectiveMetrics]],
) -> None:
    print("\nObjective Metrics Comparison")
    print("model                         MAE|1-Q|  zeroQ_MAE  CP_RMSE  dRPM_peak  dCP_peak")
    for model_label, metrics in objective_rows:
        print(
            f"{model_label:<28}"
            f"{metrics.mean_abs_torque_ratio_error:9.3f}"
            f"{metrics.zero_pitch_mean_abs_torque_ratio_error:11.3f}"
            f"{metrics.cp_curve_rmse:9.3f}"
            f"{metrics.cp_peak_rpm_error:11.3f}"
            f"{metrics.cp_peak_value_error:10.3f}"
        )


def _write_csv(path: Path, rows: list[RotorModelComparison]) -> None:
    flat_rows = [comparison_as_dict(row) for row in rows]
    if not flat_rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat_rows[0].keys()))
        writer.writeheader()
        writer.writerows(flat_rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blade-file", type=Path, help="Path to the IEA-15 MW WindIO YAML.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Optional directory where the reference-curve and TSR-sweep CSV files are written.",
    )
    parser.add_argument("--n-spanwise-panels", type=int, default=12)
    parser.add_argument("--n-chordwise-panels", type=int, default=2)
    parser.add_argument(
        "--model",
        choices=(
            "vlm",
            "hybrid",
            "experimental",
            "experimental-gated",
            "experimental-gated-tsr",
            "both",
            "all",
        ),
        default="both",
        help="Which reduced-order model to compare against BEM.",
    )
    parser.add_argument(
        "--experimental-lift-gain",
        type=float,
        default=0.15,
        help="Gain for experimental delta-Cl trailing lift closure (model=experimental/all).",
    )
    parser.add_argument(
        "--experimental-alpha-gate-deg",
        type=float,
        default=6.0,
        help="Incident-angle gate for experimental-gated model (degrees).",
    )
    parser.add_argument(
        "--experimental-tsr-cutoff",
        type=float,
        default=9.5,
        help="TSR cutoff above which TSR-weighted experimental mode damps lift correction.",
    )
    parser.add_argument(
        "--experimental-tsr-falloff",
        type=float,
        default=3.0,
        help="TSR falloff width for TSR-weighted experimental mode.",
    )
    parser.add_argument(
        "--experimental-radial-tip-scale",
        type=float,
        default=1.0,
        help="Tip-side scaling factor for radial weighting in TSR-weighted mode (1.0 disables radial weighting).",
    )
    parser.add_argument(
        "--experimental-radial-power",
        type=float,
        default=1.0,
        help="Exponent for radial weighting in TSR-weighted mode.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    blade_file = (
        args.blade_file.resolve() if args.blade_file else _resolve_default_blade_file(repo_root)
    )

    model_specs = {
        "vlm": ("Reduced VLM", compare_vlm_to_bem_sweep),
        "hybrid": (
            "VLM + Incident Polar Drag Correction",
            compare_vlm_polar_corrected_to_bem_sweep,
        ),
        "experimental": (
            f"VLM + Experimental Delta-Cl Trailing Lift (gain={args.experimental_lift_gain:.3f})",
            lambda blade, points, **kwargs: compare_vlm_experimental_lift_to_bem_sweep(
                blade,
                points,
                experimental_lift_gain=args.experimental_lift_gain,
                **kwargs,
            ),
        ),
        "experimental-gated": (
            "VLM + Experimental Gated Positive Delta-Cl Trailing Lift "
            f"(gain={args.experimental_lift_gain:.3f}, gate={args.experimental_alpha_gate_deg:.1f} deg)",
            lambda blade, points, **kwargs: compare_vlm_experimental_gated_lift_to_bem_sweep(
                blade,
                points,
                experimental_lift_gain=args.experimental_lift_gain,
                experimental_alpha_gate_deg=args.experimental_alpha_gate_deg,
                **kwargs,
            ),
        ),
        "experimental-gated-tsr": (
            "VLM + Experimental Gated Positive Delta-Cl Trailing Lift + TSR Damping "
            "("
            f"gain={args.experimental_lift_gain:.3f}, "
            f"gate={args.experimental_alpha_gate_deg:.1f} deg, "
            f"tsr_cutoff={args.experimental_tsr_cutoff:.2f}, "
            f"tsr_falloff={args.experimental_tsr_falloff:.2f}, "
            f"tip_scale={args.experimental_radial_tip_scale:.2f}, "
            f"radial_power={args.experimental_radial_power:.2f}"
            ")",
            lambda blade, points, **kwargs: compare_vlm_experimental_gated_tsr_lift_to_bem_sweep(
                blade,
                points,
                experimental_lift_gain=args.experimental_lift_gain,
                experimental_alpha_gate_deg=args.experimental_alpha_gate_deg,
                experimental_tsr_cutoff=args.experimental_tsr_cutoff,
                experimental_tsr_falloff=args.experimental_tsr_falloff,
                experimental_radial_tip_scale=args.experimental_radial_tip_scale,
                experimental_radial_power=args.experimental_radial_power,
                **kwargs,
            ),
        ),
    }
    if args.model == "both":
        requested_models = ("vlm", "hybrid")
    elif args.model == "all":
        requested_models = (
            "vlm",
            "hybrid",
            "experimental",
            "experimental-gated",
            "experimental-gated-tsr",
        )
    else:
        requested_models = (args.model,)

    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=True)

    objective_rows: list[tuple[str, RotorValidationObjectiveMetrics]] = []

    for model_key in requested_models:
        model_label, compare_fn = model_specs[model_key]
        reference_rows = compare_fn(
            blade_file,
            IEA15MW_REFERENCE_OPERATING_POINTS,
            n_spanwise_panels=args.n_spanwise_panels,
            n_chordwise_panels=args.n_chordwise_panels,
        )
        tsr_rows = compare_fn(
            blade_file,
            IEA15MW_ZERO_PITCH_TSR_SWEEP,
            n_spanwise_panels=args.n_spanwise_panels,
            n_chordwise_panels=args.n_chordwise_panels,
        )

        _print_table(f"IEA-15 MW Reference Operating Points ({model_label})", reference_rows)
        _print_table(f"IEA-15 MW Zero-Pitch TSR Sweep ({model_label})", tsr_rows)
        objective_rows.append((model_label, _print_summary(model_label, reference_rows, tsr_rows)))

        if args.output_dir:
            _write_csv(
                args.output_dir / f"iea15mw_reference_curve_comparison_{model_key}.csv",
                reference_rows,
            )
            _write_csv(
                args.output_dir / f"iea15mw_zero_pitch_tsr_sweep_{model_key}.csv",
                tsr_rows,
            )

        if len(requested_models) > 1 and model_key != requested_models[-1]:
            print()

    if len(objective_rows) > 1:
        _print_objective_scoreboard(objective_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
