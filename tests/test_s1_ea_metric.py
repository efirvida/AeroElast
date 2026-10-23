"""S-1 EA(CLT) metric: hoop-stress-free axial stiffness 1/(A^-1)[0,0] (tools/run_s1_sectional)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))
sys.path.insert(0, str(REPO_ROOT))

pytest.importorskip("scipy")
from run_s1_sectional import laminate_a11_eff  # noqa: E402

_MATS = {
    "ud": {"elastic": {"E": [44.6e9, 14.3e9, 14.3e9], "nu": [0.3, 0.3, 0.3], "G": [3.0e9]}},
    "iso": {"elastic": {"E": 3.44e9, "nu": 0.3}},
    "iso0": {"elastic": {"E": 3.44e9, "nu": 0.0}},
}


def test_single_orthotropic_ply_is_e1_h():
    """One 0 deg ply, hoop free: 1/(A^-1)11 = E1 h exactly, for any nu12 (A11 itself is E1 h / D)."""
    h = 0.01
    assert laminate_a11_eff([("ud", h, 0.0)], _MATS) == pytest.approx(44.6e9 * h, rel=1e-12)


def test_symmetric_two_material_stack_closed_form():
    """[iso / ud(0) / iso]: A is diagonal-coupled only through Q12; compare with the 2x2 Schur form."""
    h_i, h_u = 0.004, 0.012
    e1, e2, nu12 = 44.6e9, 14.3e9, 0.3
    ei, nui = 3.44e9, 0.3
    d_u = 1.0 - nu12 * nu12 * e2 / e1
    d_i = 1.0 - nui**2
    a11 = 2 * h_i * ei / d_i + h_u * e1 / d_u
    a12 = 2 * h_i * nui * ei / d_i + h_u * nu12 * e2 / d_u
    a22 = 2 * h_i * ei / d_i + h_u * e2 / d_u
    expected = a11 - a12**2 / a22  # 0 deg plies: no shear coupling, Schur complement of the hoop row
    got = laminate_a11_eff([("iso", h_i, 0.0), ("ud", h_u, 0.0), ("iso", h_i, 0.0)], _MATS)
    assert got == pytest.approx(expected, rel=1e-12)
    assert got < a11


def test_isotropic_nu_zero_is_unchanged():
    """nu = 0 decouples axial and hoop: the effective value equals the old A11 = E h."""
    h = 0.02
    assert laminate_a11_eff([("iso0", h, 37.0)], _MATS) == pytest.approx(3.44e9 * h, rel=1e-12)
    assert laminate_a11_eff([("iso0", h, 0.0)], _MATS) == pytest.approx(3.44e9 * h, rel=1e-12)


def test_matches_inverse_of_full_a_matrix():
    """Generic angled stack: equals 1/inv(A)[0,0] built independently from the ply Q."""
    layup = [("ud", 0.003, 45.0), ("ud", 0.005, -45.0), ("ud", 0.004, 0.0)]
    a = np.zeros((3, 3))
    e1, e2, nu12, g12 = 44.6e9, 14.3e9, 0.3, 3.0e9
    den = 1.0 - nu12 * nu12 * e2 / e1
    q = np.array([[e1 / den, nu12 * e2 / den, 0], [nu12 * e2 / den, e2 / den, 0], [0, 0, g12]])
    for _, h, ang in layup:
        t = np.deg2rad(ang)
        c, s = np.cos(t), np.sin(t)
        tm = np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])
        r = np.diag([1.0, 1.0, 2.0])
        a += h * np.linalg.inv(tm) @ q @ r @ tm @ np.linalg.inv(r)
    assert laminate_a11_eff(layup, _MATS) == pytest.approx(1.0 / np.linalg.inv(a)[0, 0], rel=1e-9)
