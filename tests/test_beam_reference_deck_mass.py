"""Deck-integrated blade mass (`tools/beam_reference.py`, K11.2): closed forms, no aeroelast."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from beam_reference import beamdyn_deck_mass, reference_axis_length  # noqa: E402


def test_straight_axis_length_is_exact():
    kp = np.column_stack([np.zeros(5), np.zeros(5), np.linspace(0.0, 117.0, 5)])
    assert reference_axis_length(kp) == pytest.approx(117.0, rel=1e-12)


def test_quarter_circle_axis_length_converges_to_pi_r_over_2():
    theta = np.linspace(0.0, np.pi / 2, 9)
    kp = np.column_stack([10.0 * (1 - np.cos(theta)), np.zeros(9), 10.0 * np.sin(theta)])
    assert reference_axis_length(kp, refine=64) == pytest.approx(5.0 * np.pi, rel=1e-4)


def test_linear_mass_density_integrates_exactly():
    # BeamDyn interpolates sectional properties linearly in eta: the trapezoid is exact for that.
    eta = np.array([0.0, 0.3, 1.0])
    m = np.array([1000.0, 400.0, 100.0])
    expected = 50.0 * (0.3 * (1000 + 400) / 2 + 0.7 * (400 + 100) / 2)
    assert beamdyn_deck_mass(eta, m, 50.0) == pytest.approx(expected, rel=1e-12)


def test_rejects_eta_outside_unit_interval_or_unsorted():
    with pytest.raises(ValueError):
        beamdyn_deck_mass(np.array([0.0, 0.6, 0.5]), np.ones(3), 1.0)
    with pytest.raises(ValueError):
        beamdyn_deck_mass(np.array([0.1, 0.5, 1.0]), np.ones(3), 1.0)
