"""A11 of the S-1 EA(CLT) diagnostic against closed forms (tools/run_s1_sectional.laminate_a11)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))

pytest.importorskip("scipy")
from run_s1_sectional import laminate_a11  # noqa: E402

_MATS = {
    "ud": {"elastic": {"E": [44.6e9, 14.3e9, 14.3e9], "nu": [0.3, 0.3, 0.3], "G": [3.0e9]}},
    "iso": {"elastic": {"E": 3.44e9, "nu": 0.3}},
}


def test_unidirectional_ply_uses_reciprocal_poisson_ratio():
    """A 0 deg UD ply: A11 = E1 h / (1 - nu12 nu21), nu21 = nu12 E2 / E1 (not 1 - nu12^2)."""
    h = 0.01
    nu12, e1, e2 = 0.3, 44.6e9, 14.3e9
    expected = e1 * h / (1.0 - nu12 * (nu12 * e2 / e1))
    assert laminate_a11([("ud", h, 0.0)], _MATS) == pytest.approx(expected, rel=1e-12)


def test_isotropic_ply_is_e_over_one_minus_nu_squared():
    h = 0.02
    assert laminate_a11([("iso", h, 37.0)], _MATS) == pytest.approx(3.44e9 * h / (1 - 0.3**2), rel=1e-12)


def test_ninety_degree_ply_carries_the_transverse_modulus():
    h = 0.01
    nu12, e1, e2 = 0.3, 44.6e9, 14.3e9
    expected = e2 * h / (1.0 - nu12 * (nu12 * e2 / e1))
    assert laminate_a11([("ud", h, 90.0)], _MATS) == pytest.approx(expected, rel=1e-12)
