"""The [0/90/45/-45]s cantilever against its classical-lamination closed form.

The laminate comparison with CalculiX (``test_shell_convergence``) has a CalculiX sequence that does
not converge (self-convergence orders -0.07 and -0.09), so that code is not an independent
arbiter of the converged answer. This test supplies one that is: the in-plane bending of the strip
(load in the width direction, depth ``B``) as a Timoshenko beam whose axial and shear moduli come
from the laminate's membrane matrix ``A`` (built here with plain numpy, not with the element code):

    w_tip = P L^3 / (3 E_x t B^3 / 12) + P L / (kappa G_xy t B),   kappa = 5/6

with ``E_x = 1 / (a_11 t)`` and ``G_xy = 1 / (a_33 t)``, ``a = A^-1``.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.assertions import assert_relative_error  # noqa: E402

# Imported lazily by the parity module itself guards petsc4py; this file needs only its helpers.
parity = pytest.importorskip(
    "tests.validation.parity.test_composite_beam_parity",
    reason="composite parity helpers unavailable",
)
convergence = pytest.importorskip(
    "tests.validation.parity.test_shell_convergence",
    reason="shell convergence helpers unavailable",
)

PLIES = [0.0, 90.0, 45.0, -45.0, -45.0, 45.0, 90.0, 0.0]
KAPPA = 5.0 / 6.0
# Timoshenko beam theory of a strip with B / L = 0.1 is good to a percent or two (clamped-end
# warping and the anticlastic effect are neglected); the bound was fixed from that reasoning.
TOL = 0.02


def _q_bar(q: np.ndarray, theta_deg: float) -> np.ndarray:
    c, s = np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg))
    t = np.array(
        [[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]]
    )
    r = np.diag([1.0, 1.0, 2.0])
    return np.linalg.inv(t) @ q @ r @ t @ np.linalg.inv(r)


def laminate_tip_closed_form(load: float) -> float:
    e1, e2, g12, nu12 = parity.E1, parity.E2, parity.G12, parity.nu12
    nu21 = nu12 * e2 / e1
    den = 1.0 - nu12 * nu21
    q = np.array(
        [[e1 / den, nu12 * e2 / den, 0.0], [nu12 * e2 / den, e2 / den, 0.0], [0.0, 0.0, g12]]
    )
    t_total = parity.thickness
    a_mat = sum(_q_bar(q, th) * (t_total / len(PLIES)) for th in PLIES)
    a_inv = np.linalg.inv(a_mat)
    e_x = 1.0 / (a_inv[0, 0] * t_total)
    g_xy = 1.0 / (a_inv[2, 2] * t_total)
    length, depth = convergence.L, convergence.B
    bending = load * length**3 / (3.0 * e_x * t_total * depth**3 / 12.0)
    shear = load * length / (KAPPA * g_xy * t_total * depth)
    return bending + shear


def test_laminate_tip_matches_the_classical_lamination_closed_form():
    nx, ny = convergence.LAMINATE_MESHES[-1]
    tip = convergence._laminate_aero_tip_x(parity, nx, ny)
    expected = laminate_tip_closed_form(convergence._LAMINATE_LOAD_X)
    assert_relative_error(
        tip,
        expected,
        tol=TOL,
        reference_name="classical lamination theory, Timoshenko strip in in-plane bending",
        kind="analytical",
        what=f"[0/90/45/-45]s laminate tip displacement, {nx} x {ny} elements",
    )
