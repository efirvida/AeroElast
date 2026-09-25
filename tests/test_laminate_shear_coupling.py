"""Analytic verification of the laminate shear-coupling (the "16" terms).

The blade's bend-twist coupling — the mechanism that produces the large
aeroelastic twist measured in the FSI (`bem_diag.csv`: +7.4 deg at r/R=0.9)
— needs the extension-shear terms Qbar_16 / A16 to be right.  Those terms
had NO test at any level: the composite suite covers the B matrix
(extension-bending) and CalculiX parity covers axial->bending only.

This locks in the closed-form behaviour of `compute_Qbar` and the A matrix:

  * Qbar_16 / Qbar_26 against the textbook transformation,
  * a balanced angle-ply (equal +-45) gives A16 = A26 = 0,
  * an unbalanced angle-ply gives the analytic non-zero value, with the
    correct sign.

No external dependencies, runs in milliseconds.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.laminate import compute_Q, compute_Qbar, create_laminate_from_angles
from aeroelast.core.material import Material

E1, E2, G12, NU12 = 44.6e9, 17.0e9, 3.27e9, 0.262
PLY_T = 5.0e-3


@pytest.fixture(scope="module")
def mat() -> Material:
    return Material(
        name="glass_uni_test",
        E=(E1, E2, E2),
        G=(G12, G12, G12),
        nu=(NU12, NU12, 0.0),
        rho=1940.0,
    )


def _closed_form_16_26(mat: Material, theta_deg: float) -> tuple[float, float]:
    """Textbook Qbar_16 and Qbar_26 for a rotated lamina."""
    Q = compute_Q(mat)
    q11, q12, q22, q66 = Q[0, 0], Q[0, 1], Q[1, 1], Q[2, 2]
    c = np.cos(np.radians(theta_deg))
    s = np.sin(np.radians(theta_deg))
    q16 = (q11 - q12 - 2 * q66) * c**3 * s - (q22 - q12 - 2 * q66) * c * s**3
    q26 = (q11 - q12 - 2 * q66) * c * s**3 - (q22 - q12 - 2 * q66) * c**3 * s
    return float(q16), float(q26)


@pytest.mark.parametrize("theta", [30.0, 45.0, 60.0, -45.0])
def test_qbar_16_26_match_closed_form(mat: Material, theta: float) -> None:
    """The transformed shear-coupling terms must match the textbook form."""
    q16, q26 = _closed_form_16_26(mat, theta)
    Qb = compute_Qbar(mat, theta)
    assert Qb[0, 2] == pytest.approx(q16, rel=1e-10), f"Qbar_16 at {theta} deg"
    assert Qb[1, 2] == pytest.approx(q26, rel=1e-10), f"Qbar_26 at {theta} deg"


def test_qbar_16_flips_sign_with_plying_angle(mat: Material) -> None:
    """+45 and -45 must give opposite shear coupling."""
    q16_p = compute_Qbar(mat, 45.0)[0, 2]
    q16_m = compute_Qbar(mat, -45.0)[0, 2]
    assert q16_p > 0.0 > q16_m
    assert q16_p == pytest.approx(-q16_m, rel=1e-10)


def test_balanced_angle_ply_has_zero_shear_coupling(mat: Material) -> None:
    """Equal +-45 counts -> A16 = A26 = 0 (the code's own documented invariant)."""
    lam = create_laminate_from_angles(mat, PLY_T, [45.0, -45.0, -45.0, 45.0])
    A = lam.get_ABD_matrix()[:3, :3]
    # A is [A11 A12 A16; A12 A22 A26; A16 A26 A66]
    assert np.abs(A[0, 2]) < 1e-3 * np.abs(A[0, 0])
    assert np.abs(A[1, 2]) < 1e-3 * np.abs(A[1, 1])


def test_unbalanced_angle_ply_shear_coupling_matches_closed_form(mat: Material) -> None:
    """An unbalanced stack's A16 must equal sum(Qbar_16 * t) over the plies."""
    angles = [45.0, 0.0, 0.0, 45.0]
    lam = create_laminate_from_angles(mat, PLY_T, angles)
    A = lam.get_ABD_matrix()[:3, :3]
    expected = sum(compute_Qbar(mat, a)[0, 2] * PLY_T for a in angles)
    assert A[0, 2] == pytest.approx(expected, rel=1e-9)
    assert A[0, 2] > 0.0  # two +45 plies, no -45: net positive shear coupling


def test_pure_angle_ply_coupling_is_four_plies(mat: Material) -> None:
    """[45]s: A16 = 4 * Qbar_16(45) * t and it is 2x the two-ply unbalanced case."""
    lam4 = create_laminate_from_angles(mat, PLY_T, [45.0] * 4)
    lam2 = create_laminate_from_angles(mat, PLY_T, [45.0, 0.0, 0.0, 45.0])
    a4 = lam4.get_ABD_matrix()[0, 2]
    a2 = lam2.get_ABD_matrix()[0, 2]
    assert a4 == pytest.approx(4.0 * compute_Qbar(mat, 45.0)[0, 2] * PLY_T, rel=1e-9)
    assert a4 == pytest.approx(2.0 * a2, rel=1e-9)
