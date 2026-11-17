"""tools/l5_t9_verdict.py: model-vs-solid comparison and the written psi agreement tolerance (L5 T9c)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import l5_t9_verdict as V  # noqa: E402

REF = {"psi_deg": 2.5, "ratio_fe_ff": 0.02, "GJ": 2.0e7, "EI_flap": 4.0e8, "EI_edge": 8.0e8}


def test_compare_gives_differences_for_angles_and_ratios_for_stiffness():
    m = {"psi_deg": 2.9, "ratio_fe_ff": 0.026, "GJ": 3.0e7, "EI_flap": 4.0e8, "EI_edge": 7.2e8}
    c = V.compare(m, REF)
    assert c["dpsi"] == pytest.approx(0.4)
    assert c["dfe"] == pytest.approx(0.006)
    assert (c["GJ"], c["EI_flap"], c["EI_edge"]) == pytest.approx((1.5, 1.0, 0.9))
    assert c["agrees_psi"] is True


def test_psi_agreement_band_is_one_degree_and_is_not_widened_by_the_sign():
    assert V.TOL_PSI_AGREE == 1.0
    assert V.compare(REF | {"psi_deg": 3.6}, REF)["agrees_psi"] is False
    assert V.compare(REF | {"psi_deg": 1.4}, REF)["agrees_psi"] is False
    assert V.compare(REF | {"psi_deg": 1.6}, REF)["agrees_psi"] is True
