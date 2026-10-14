"""StaticNonlinearSolver ``scheme: newton`` (total-Lagrangian Newton) against the elastica.

The legacy updated-Lagrangian scheme (``scheme: ul``) integrates ``K_e(x_ref) du = dlam F`` with the linear
elastic stiffness on the current geometry: no equilibrium iteration and no initial-stress
stiffness.  Its limit as the load stepping refines is not the dead-load equilibrium.  Measured
2026-10-14 on the cantilever strip of ``test_static_nonlinear_ul_elastica.py`` (alpha = P L^2 / EI
= 2, mesh and step converged): w/L is +5.29 % above the exact elastica and the tip shortening
u/L is 13 % too large; the gap grows as the strip narrows (+4.7, +8.8, +10.7 % for b = 0.1, 0.05,
0.025).  A 2-D corotational beam with the same scheme reaches +11.3 % at alpha = 2, while the same
beam with Newton iteration on the consistent tangent reproduces the elastica to 4 digits.

``scheme: newton`` (now the default) iterates the residual ``f_int(u) - lam F`` on the Green-Lagrange internal
force and its tangent.  Its objectivity for rotations about any axis is checked in
``test_twisted_beam_nlgeom_calculix.py``.

Reference: the exact large-deflection solution of a cantilever under a dead tip load
(Bisshopp & Drucker 1945; tabulated via Bathe & Bolourchi 1979):

    alpha = P L^2 / EI     w/L (tip, load direction)    u/L (tip, axial)
        1                    0.30172                       -0.05643
        2                    0.49346                       -0.16064

The elastica is a beam.  A shell strip is a plate: a wide strip is stiffer (anticlastic and
Ashwell effects, about 5 % at b/L = 0.1 and h/L = 1e-3).  The comparison against the
elastica is therefore made on a narrow strip (b = 0.025), where the plate effect is below 0.5 %.
"""

from __future__ import annotations

import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.assertions import assert_residual_below
from tests.support.strip_cantilever import ELASTICA, ELASTICA_REFERENCE, solve_strip


@pytest.mark.parametrize("alpha", [1.0, 2.0])
def test_newton_scheme_narrow_strip_matches_the_elastica(alpha):
    """Narrow strip, mesh and step converged: tip deflection and shortening vs the elastica."""
    w_ref, u_ref = ELASTICA[alpha]
    w, u = solve_strip(alpha, b=0.025, nx=64, ny=2, steps=8, scheme="newton")
    assert_residual_below(
        abs(w - w_ref) / w_ref,
        tol=0.01,
        kind="paper",
        reference_name=ELASTICA_REFERENCE,
        what=f"alpha={alpha} tip w/L, narrow strip, newton scheme",
    )
    assert_residual_below(
        abs(u - u_ref) / abs(u_ref),
        tol=0.01,
        kind="paper",
        reference_name=ELASTICA_REFERENCE,
        what=f"alpha={alpha} tip u/L, narrow strip, newton scheme",
    )
