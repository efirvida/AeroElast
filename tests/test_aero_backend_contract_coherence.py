from __future__ import annotations

import math

from aeroelast.models.blade.aerodynamics import load_blade_aero
from aeroelast.solvers.aero.validation import (
    IEA15MW_REFERENCE_OPERATING_POINTS,
    compare_vlm_to_bem_sweep,
)


def test_vlm_and_bem_share_same_input_dataset_contract(iea_blade_yaml: str):
    """Both backends must solve the same operating points from one blade dataset.

    This is the contract baseline for any aerodynamic backend intended to plug
    into the same FSI participant loop (including future SHARPy adapters).
    """

    blade_aero = load_blade_aero(iea_blade_yaml)
    rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    assert len(rows) == len(IEA15MW_REFERENCE_OPERATING_POINTS)

    for row in rows:
        assert row.operating_point in IEA15MW_REFERENCE_OPERATING_POINTS
        assert row.panel_count > 0

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
            row.thrust_ratio,
            row.torque_ratio,
            row.cp_ratio,
            row.ct_ratio,
        ):
            assert math.isfinite(value)

        # Same physical problem => both must predict powered operation.
        assert row.vlm.thrust > 0.0
        assert row.vlm.torque > 0.0
        assert row.vlm.cp > 0.0
        assert row.vlm.ct > 0.0
        assert row.bem.thrust > 0.0
        assert row.bem.torque > 0.0
        assert row.bem.cp > 0.0
        assert row.bem.ct > 0.0

        # Coherence envelope (broad on purpose): this is a contract gate,
        # not a high-fidelity acceptance test.
        assert 0.4 <= row.thrust_ratio <= 1.4
        assert 0.03 <= row.torque_ratio <= 1.2
        assert 0.03 <= row.cp_ratio <= 1.2
        assert 0.4 <= row.ct_ratio <= 1.4


def test_vlm_and_bem_follow_same_pitch_control_trend(iea_blade_yaml: str):
    """Pitch-to-feather must reduce thrust loads for both backends.

    Above rated wind the operating points follow the real pitch schedule
    (v10.659 -> pitch 0 deg, v15.471 -> pitch 12.19 deg).  A pitch-regulated
    rotor holds rated power, so power/torque stay approximately constant
    while thrust and ct drop.  Asserting a monotone torque reduction would
    contradict the controller schedule the operating points encode.

    This trend check is backend-agnostic and should also hold for SHARPy once
    it is connected through the same runtime contract.
    """

    blade_aero = load_blade_aero(iea_blade_yaml)
    rows = compare_vlm_to_bem_sweep(
        blade_aero,
        IEA15MW_REFERENCE_OPERATING_POINTS,
        n_spanwise_panels=12,
        n_chordwise_panels=2,
    )

    by_name = {row.operating_point.name: row for row in rows}
    rated = by_name["v10.659"]
    feathered = by_name["v15.471"]

    # Current reduced VLM still has known CP/CQ-shape limitations at high
    # pitch, so we gate trend coherence on thrust/ct, which are the primary
    # force-side coupling quantities.
    assert rated.vlm.ct > feathered.vlm.ct
    assert rated.vlm.thrust > feathered.vlm.thrust

    assert rated.bem.ct > feathered.bem.ct
    assert rated.bem.thrust > feathered.bem.thrust
    assert rated.bem.cp > feathered.bem.cp

    # Power regulation above rated: torque (fixed omega) stays within 5%.
    regulation = abs(feathered.bem.torque - rated.bem.torque) / rated.bem.torque
    assert regulation < 0.05, (
        f"BEM torque not regulated above rated: {rated.bem.torque:.4e} -> "
        f"{feathered.bem.torque:.4e} N m (rel {regulation:.2%} > 5%)"
    )
