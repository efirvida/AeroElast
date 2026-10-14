"""Software checks of the ``scheme`` option of StaticNonlinearSolver (no physical reference)."""

from __future__ import annotations

import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.validation.element.test_static_nonlinear_newton_elastica import solve_strip


def test_newton_scheme_does_not_depend_on_the_load_stepping():
    """A consistent Newton result is a property of the equilibrium, not of the path."""
    w4, u4 = solve_strip(2.0, b=0.025, nx=32, ny=2, steps=4, scheme="newton")
    w16, u16 = solve_strip(2.0, b=0.025, nx=32, ny=2, steps=16, scheme="newton")
    assert abs(w4 - w16) / w16 < 1e-3
    assert abs(u4 - u16) / abs(u16) < 1e-3


def test_default_scheme_is_the_updated_lagrangian_one_and_unchanged():
    """No ``scheme`` key keeps the recorded UL reading (w/L = 0.29075, 32 steps, 8x4 mesh)."""
    w, _ = solve_strip(1.0, b=0.1, nx=8, ny=4, steps=32, scheme=None)
    assert w == pytest.approx(0.29075, abs=2e-5)
    w_explicit, _ = solve_strip(1.0, b=0.1, nx=8, ny=4, steps=32, scheme="ul")
    assert w_explicit == pytest.approx(w, abs=1e-12)


def test_unknown_scheme_is_rejected():
    with pytest.raises(ValueError, match="scheme"):
        solve_strip(1.0, b=0.1, nx=8, ny=4, steps=4, scheme="bogus")
