"""Mesh chord follows the yaml chord on the yaml's own grid.

Candidate C / audit finding G2 of the blade input audit (``odd/tasks/blade-input-audit.md``, T8):
the yaml chord has 53 points (extra tip knots 0.985, 0.99, 0.995) but was resampled linearly onto
the 50-point reference_axis grid and then pchip-interpolated, so the tip knots were lost (+0.77 %
at r/R 0.97, -7..-35 % at the tip knots). Reference: linear interpolation on the yaml grid.
"""

from pathlib import Path

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.objects.blade import Blade

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
TOL = 1e-3


@pytest.fixture(scope="module")
def blade_and_chord():
    osb = yaml.safe_load(YAML.read_text())["components"]["blade"]["outer_shape_bem"]
    grid = np.array(osb["chord"]["grid"], dtype=float)
    values = np.array(osb["chord"]["values"], dtype=float)
    blade = Blade(str(YAML))
    length = float(blade.definition.span[-1])
    return blade, grid, values, length


def _check(blade, grid, values, length, fractions):
    span = np.asarray(blade.geometry.ispan, dtype=float)
    got = np.interp(fractions * length, span, np.asarray(blade.geometry.ichord, dtype=float))
    want = np.interp(fractions, grid, values)
    np.testing.assert_allclose(got, want, rtol=TOL)


def test_chord_at_yaml_tip_knots(blade_and_chord):
    blade, grid, values, length = blade_and_chord
    knots = np.array([0.985, 0.99, 0.995])
    assert np.all(np.isin(knots, grid))
    _check(blade, grid, values, length, knots)


def test_chord_at_0p97(blade_and_chord):
    blade, grid, values, length = blade_and_chord
    _check(blade, grid, values, length, np.array([0.97]))


def test_chord_follows_yaml_at_inserted_stations(blade_and_chord):
    _, grid, values, length = blade_and_chord
    blade = Blade(str(YAML))
    new = np.array([0.5, 0.9, 0.97, 0.9725, 0.9875, 0.9925]) * length
    blade.add_interpolated_stations(new)
    span = np.asarray(blade.geometry.ispan, dtype=float)
    idx = [int(np.argmin(abs(span - s))) for s in new]
    got = np.asarray(blade.geometry.ichord, dtype=float)[idx]
    want = np.interp(span[idx] / length, grid, values)
    np.testing.assert_allclose(got, want, rtol=TOL)
