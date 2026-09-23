from __future__ import annotations

import math

import numpy as np
import pytest

ccblade = pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")

from aeroelast.models.blade.aerodynamics import load_blade_aero
from aeroelast.solvers.aero.validation import (
    IEA15MW_REFERENCE_OPERATING_POINTS,
    IEA15MW_ZERO_PITCH_TSR_SWEEP,
    build_reduced_vlm_rotor_mesh,
    compare_vlm_experimental_gated_lift_to_bem_sweep,
    compare_vlm_experimental_gated_tsr_lift_to_bem_sweep,
    compare_vlm_experimental_lift_to_bem_sweep,
    compare_vlm_polar_corrected_to_bem_sweep,
    compare_vlm_to_bem_sweep,
    compute_validation_objective_metrics,
    compute_vlm_polar_corrected_reduced_rotor_metrics,
    compute_vlm_sectional_diagnostics,
)
from aeroelast.solvers.bem.engine import BEMSolver


def test_build_reduced_vlm_rotor_mesh_has_expected_topology(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    mesh = build_reduced_vlm_rotor_mesh(
        blade_aero,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert mesh.node_count == blade_aero.n_blades * (12 + 1) * (2 + 1)
    assert mesh.elements_count == blade_aero.n_blades * 12 * 2
    assert sorted(mesh.node_sets_names) == ["rotor_blade_1", "rotor_blade_2", "rotor_blade_3"]
    assert sorted(mesh.element_sets_names) == ["rotor_blade_1", "rotor_blade_2", "rotor_blade_3"]


def test_build_reduced_vlm_rotor_mesh_preserves_prebend_and_sweep_reference_line(
    iea_blade_yaml: str,
):
    blade_aero = load_blade_aero(iea_blade_yaml)
    blade_aero.n_blades = 1

    for station in blade_aero.stations:
        eta = float(station.span_fraction)
        station.prebend = 1.2 * eta
        station.sweep = -0.6 * eta
        station.pitch_axis = 0.5

    n_spanwise_panels = 8
    mesh = build_reduced_vlm_rotor_mesh(
        blade_aero,
        n_spanwise_panels=n_spanwise_panels,
        n_chordwise_panels=1,
    )

    coords = mesh.coords_array
    radii_measured = coords[:, 2]

    span_grid = np.linspace(
        blade_aero.stations[0].span_fraction,
        blade_aero.stations[-1].span_fraction,
        n_spanwise_panels + 1,
        dtype=float,
    )
    station_span = np.array(
        [station.span_fraction for station in blade_aero.stations],
        dtype=float,
    )
    expected_radii = np.interp(span_grid, station_span, blade_aero.r)
    expected_prebend = np.interp(span_grid, station_span, blade_aero.prebend)
    expected_sweep = np.interp(span_grid, station_span, blade_aero.sweep)
    expected_prebend = expected_prebend - expected_prebend[0]
    expected_sweep = expected_sweep - expected_sweep[0]

    prebend_measured = []
    sweep_measured = []
    for expected_radius in expected_radii:
        station_mask = np.isclose(radii_measured, expected_radius, atol=1.0e-9)
        station_coords = coords[station_mask]
        assert station_coords.shape[0] == 2
        prebend_measured.append(float(np.mean(station_coords[:, 1])))
        sweep_measured.append(float(np.mean(station_coords[:, 0])))

    np.testing.assert_allclose(np.asarray(prebend_measured), expected_prebend, atol=1.0e-8)
    np.testing.assert_allclose(np.asarray(sweep_measured), expected_sweep, atol=1.0e-8)
    assert np.max(np.abs(prebend_measured)) > 1.0e-3
    assert np.max(np.abs(sweep_measured)) > 1.0e-3


def test_iea15mw_reference_curve_validation_returns_finite_metrics(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert len(rows) == len(IEA15MW_REFERENCE_OPERATING_POINTS)
    assert {row.panel_count for row in rows} == {72}

    thrust_ratios = []
    for row in rows:
        thrust_ratios.append(row.thrust_ratio)
        for value in (
            row.vlm.thrust,
            row.vlm.torque,
            row.vlm.power,
            row.vlm.cp,
            row.vlm.ct,
            row.bem.thrust,
            row.bem.torque,
            row.bem.power,
            row.bem.cp,
            row.bem.ct,
        ):
            assert math.isfinite(value)

    assert min(thrust_ratios) > 0.5
    assert max(thrust_ratios) < 1.5
    assert rows[1].vlm.ct > rows[-1].vlm.ct
    assert rows[2].vlm.ct > rows[-1].vlm.ct


def test_vlm_polar_correction_improves_pitch_controlled_torque_cases(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    baseline_rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    corrected_rows = compare_vlm_polar_corrected_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    baseline_by_case = {row.operating_point.name: row for row in baseline_rows}
    pitch_controlled_rows = [row for row in corrected_rows if row.operating_point.pitch_deg > 1.0]
    zero_pitch_rows = [
        row for row in corrected_rows if abs(row.operating_point.pitch_deg) < 1.0e-12
    ]

    assert {row.panel_count for row in corrected_rows} == {72}
    assert pitch_controlled_rows
    assert zero_pitch_rows

    for row in corrected_rows:
        for value in (
            row.vlm.thrust,
            row.vlm.torque,
            row.vlm.power,
            row.vlm.cp,
            row.vlm.ct,
        ):
            assert math.isfinite(value)

    for row in pitch_controlled_rows:
        baseline_row = baseline_by_case[row.operating_point.name]
        assert row.torque_ratio > baseline_row.torque_ratio

    for row in zero_pitch_rows:
        baseline_row = baseline_by_case[row.operating_point.name]
        assert row.torque_ratio > baseline_row.torque_ratio


def test_vlm_polar_correction_breaks_monotonic_zero_pitch_cp_ramp(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    baseline_rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    corrected_rows = compare_vlm_polar_corrected_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    baseline_peak = max(baseline_rows, key=lambda row: row.vlm.cp)
    corrected_peak = max(corrected_rows, key=lambda row: row.vlm.cp)

    assert baseline_peak.operating_point.omega_rpm == max(
        row.operating_point.omega_rpm for row in baseline_rows
    )
    assert corrected_peak.operating_point.omega_rpm < max(
        row.operating_point.omega_rpm for row in corrected_rows
    )


def test_validation_objective_metrics_are_finite_and_consistent(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    reference_rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    tsr_rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    metrics = compute_validation_objective_metrics(reference_rows, tsr_rows)

    for value in (
        metrics.mean_abs_thrust_ratio_error,
        metrics.mean_abs_torque_ratio_error,
        metrics.mean_abs_power_ratio_error,
        metrics.mean_abs_cp_error,
        metrics.zero_pitch_mean_abs_torque_ratio_error,
        metrics.cp_curve_rmse,
        metrics.cp_peak_rpm_error,
        metrics.cp_peak_value_error,
        metrics.bem_cp_peak_rpm,
        metrics.model_cp_peak_rpm,
        metrics.bem_cp_peak,
        metrics.model_cp_peak,
    ):
        assert math.isfinite(value)

    assert metrics.reference_sample_count == len(IEA15MW_REFERENCE_OPERATING_POINTS)
    assert metrics.tsr_sample_count == len(IEA15MW_ZERO_PITCH_TSR_SWEEP)


def test_vlm_sectional_diagnostics_have_expected_strip_shapes(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)
    diagnostics = compute_vlm_sectional_diagnostics(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS[0],
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    strip_count = 12
    arrays = (
        diagnostics.radii,
        diagnostics.alpha_incident_deg,
        diagnostics.alpha_total_deg,
        diagnostics.cl_polar,
        diagnostics.cd_polar,
        diagnostics.cl_vlm_trailing,
        diagnostics.gamma_leading,
        diagnostics.gamma_trailing,
        diagnostics.strip_base_thrust,
        diagnostics.strip_base_torque,
    )
    for array in arrays:
        assert len(array) == strip_count
        assert all(math.isfinite(value) for value in array)

    assert all(
        diagnostics.radii[index + 1] > diagnostics.radii[index] for index in range(strip_count - 1)
    )
    assert (
        max(
            abs(inc - total)
            for inc, total in zip(diagnostics.alpha_incident_deg, diagnostics.alpha_total_deg)
        )
        > 1.0e-6
    )


def test_experimental_lift_mode_runs_and_returns_finite_metrics(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    baseline_hybrid_rows = compare_vlm_polar_corrected_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    experimental_rows = compare_vlm_experimental_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        experimental_lift_gain=0.15,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert len(experimental_rows) == len(IEA15MW_REFERENCE_OPERATING_POINTS)
    assert {row.panel_count for row in experimental_rows} == {72}

    for row in experimental_rows:
        for value in (
            row.vlm.thrust,
            row.vlm.torque,
            row.vlm.power,
            row.vlm.cp,
            row.vlm.ct,
            row.vlm.cq,
            row.thrust_ratio,
            row.torque_ratio,
            row.power_ratio,
        ):
            assert math.isfinite(value)

    assert any(
        abs(exp.vlm.torque - base.vlm.torque) > 1.0e-9
        for exp, base in zip(experimental_rows, baseline_hybrid_rows)
    )


def test_drag_only_mode_is_default_and_matches_explicit_selection(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)
    operating_point = IEA15MW_REFERENCE_OPERATING_POINTS[0]

    default_metrics, default_panels = compute_vlm_polar_corrected_reduced_rotor_metrics(
        blade_aero,
        operating_point,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    explicit_metrics, explicit_panels = compute_vlm_polar_corrected_reduced_rotor_metrics(
        blade_aero,
        operating_point,
        lift_mode="drag_only",
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert default_panels == explicit_panels
    assert default_metrics == explicit_metrics


def test_reduced_vlm_inboard_negative_torque_extent_matches_bem(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)
    operating_point = next(
        point for point in IEA15MW_ZERO_PITCH_TSR_SWEEP if abs(point.omega_rpm - 7.0) < 1.0e-12
    )

    diagnostics = compute_vlm_sectional_diagnostics(
        blade_aero,
        operating_point,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    bem = BEMSolver(blade_aero).compute(
        v_inf=operating_point.wind_speed,
        omega=operating_point.omega_rpm,
        pitch=operating_point.pitch_deg,
    )

    bem_negative_radii = [radius for radius, tp in zip(bem.r, bem.Tp) if tp < 0.0]
    vlm_negative_radii = [
        radius for radius, dq in zip(diagnostics.radii, diagnostics.strip_base_torque) if dq < 0.0
    ]

    assert bem_negative_radii
    assert vlm_negative_radii

    bem_last_negative_radius = max(bem_negative_radii)
    vlm_last_negative_radius = max(vlm_negative_radii)

    # Reduced VLM should not keep braking deep into mid-span once formulation is corrected.
    assert vlm_last_negative_radius <= bem_last_negative_radius + 10.0


def test_experimental_gated_mode_improves_zero_pitch_torque_error(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    hybrid_reference = compare_vlm_polar_corrected_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    hybrid_tsr = compare_vlm_polar_corrected_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    gated_reference = compare_vlm_experimental_gated_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    gated_tsr = compare_vlm_experimental_gated_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    hybrid_metrics = compute_validation_objective_metrics(hybrid_reference, hybrid_tsr)
    gated_metrics = compute_validation_objective_metrics(gated_reference, gated_tsr)

    assert gated_metrics.zero_pitch_mean_abs_torque_ratio_error < (
        hybrid_metrics.zero_pitch_mean_abs_torque_ratio_error - 0.05
    )


def test_experimental_gated_tsr_mode_improves_cp_peak_rpm_mismatch(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)

    gated_reference = compare_vlm_experimental_gated_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    gated_tsr = compare_vlm_experimental_gated_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    gated_tsr_weighted_reference = compare_vlm_experimental_gated_tsr_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        experimental_tsr_cutoff=9.5,
        experimental_tsr_falloff=3.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    gated_tsr_weighted_tsr = compare_vlm_experimental_gated_tsr_lift_to_bem_sweep(
        blade_aero,
        IEA15MW_ZERO_PITCH_TSR_SWEEP,
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        experimental_tsr_cutoff=9.5,
        experimental_tsr_falloff=3.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    gated_metrics = compute_validation_objective_metrics(gated_reference, gated_tsr)
    weighted_metrics = compute_validation_objective_metrics(
        gated_tsr_weighted_reference,
        gated_tsr_weighted_tsr,
    )

    assert weighted_metrics.cp_peak_rpm_error < gated_metrics.cp_peak_rpm_error


def test_experimental_gated_tsr_mode_radial_defaults_are_neutral(iea_blade_yaml: str):
    blade_aero = load_blade_aero(iea_blade_yaml)
    operating_point = IEA15MW_REFERENCE_OPERATING_POINTS[3]

    default_metrics, default_panels = compute_vlm_polar_corrected_reduced_rotor_metrics(
        blade_aero,
        operating_point,
        lift_mode="delta_cl_trailing_positive_gated_tsr_weighted",
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        experimental_tsr_cutoff=9.5,
        experimental_tsr_falloff=3.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )
    explicit_metrics, explicit_panels = compute_vlm_polar_corrected_reduced_rotor_metrics(
        blade_aero,
        operating_point,
        lift_mode="delta_cl_trailing_positive_gated_tsr_weighted",
        experimental_lift_gain=1.0,
        experimental_alpha_gate_deg=6.0,
        experimental_tsr_cutoff=9.5,
        experimental_tsr_falloff=3.0,
        experimental_radial_tip_scale=1.0,
        experimental_radial_power=1.0,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert default_panels == explicit_panels
    assert default_metrics == explicit_metrics
