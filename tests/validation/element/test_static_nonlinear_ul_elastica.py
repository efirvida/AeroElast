"""StaticNonlinearSolver (default scheme) against the exact elastica.

The exact large-deflection solution of a cantilever under a dead tip load
(Bisshopp & Drucker 1945; tabulated via Bathe & Bolourchi 1979):

    alpha = P·L²/EI     w/L (tip, load direction)    u/L (tip, axial)
        1                 0.30172                       -0.05643
        2                 0.49346                       -0.16064

The elastica is a beam.  The strip is a plate, and a wide one (b = 0.1, h = 1e-3) is 3-5 % stiffer
than the elastica at these loads (anticlastic and Ashwell effects: CalculiX S8R NLGEOM reads it
-3.3 / -4.3 / -3.9 % at alpha = 4/7, 1, 2; see ``test_strip_elastica_calculix.py``), so the
comparison here is made on a narrow strip (b = 0.025), where the plate effect is below 0.5 %.

**History (2026-09-30 to 2026-10-14).**  This file used to run the strip b = 0.1 on an 8 x 4 mesh
through the updated-Lagrangian loop and read 3.63 % (alpha = 1) and 4.27 % (alpha = 2) from the
elastica, which the file explained as a stiff element converging "about 4 % below".  That was
wrong twice over.  On the 8 x 4 mesh the alpha = 1 gap is span discretisation (nx = 8 / 16 / 32 / 64:
-3.83 / -2.29 / -0.83 / -0.11 % at 64 steps), and its agreement at nx >= 64 was a compensation of
the plate effect (-4.3 %) by the scheme's bias (+4.4 %).  The scheme was the real defect: the
updated-Lagrangian loop integrates ``K_e(x_ref) du = dlam F`` with no equilibrium iteration and no
initial-stress stiffness, so refining the load steps does not converge to the dead-load
equilibrium (alpha = 2, mesh and step converged: w/L +5.29 %, u/L +13 %; narrow strip +10.7 %; a
2-D corotational beam with the same scheme +11.3 %, with Newton on the consistent tangent exact).
A pure end moment is exact to 0.1-0.6 % up to 115 deg, so the element kinematics were never at
fault.  ``StaticNonlinearSolver`` now defaults to a total-Lagrangian Newton (``scheme: newton``);
the old loop stays as ``scheme: ul`` and its bias is pinned in
``tests/software/test_static_nonlinear_scheme.py``.
"""

from __future__ import annotations

import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.assertions import assert_residual_below
from tests.support.strip_cantilever import ELASTICA, solve_strip

B_NARROW = 0.025


def _solve(alpha: float, *, steps: int):
    """Tip (w/L, u/L) of the narrow strip with the default scheme."""
    return solve_strip(alpha, b=B_NARROW, nx=64, ny=2, steps=steps, scheme=None)


@pytest.mark.parametrize("alpha,steps", [(1.0, 8), (2.0, 8)])
def test_default_scheme_tip_deflection_matches_the_elastica(alpha, steps):
    """Narrow strip, mesh and step converged: tip deflection and shortening vs the elastica."""
    w_ref, u_ref = ELASTICA[alpha]
    w, u = _solve(alpha, steps=steps)
    reference_name = (
        "the elastica reference table of Bisshopp & Drucker (1945), "
        "tabulated by Bathe & Bolourchi (1979)"
    )
    assert_residual_below(
        abs(w - w_ref) / w_ref,
        tol=0.01,
        kind="paper",
        reference_name=reference_name,
        what=f"alpha={alpha} tip w/L, default scheme",
    )
    assert_residual_below(
        abs(u - u_ref) / abs(u_ref),
        tol=0.01,
        kind="paper",
        reference_name=reference_name,
        what=f"alpha={alpha} tip u/L, default scheme",
    )
