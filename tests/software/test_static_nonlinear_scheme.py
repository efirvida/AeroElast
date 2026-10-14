"""Software checks of the ``scheme`` option of StaticNonlinearSolver (no physical reference)."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.conftest import ccx_bin_or_skip
from tests.validation.benchmarks.test_twisted_beam_nlgeom_calculix import _ccx_tip, _ours_tip
from tests.support.strip_cantilever import solve_strip


def test_newton_scheme_does_not_depend_on_the_load_stepping():
    """A consistent Newton result is a property of the equilibrium, not of the path."""
    w4, u4 = solve_strip(2.0, b=0.025, nx=32, ny=2, steps=4, scheme="newton")
    w16, u16 = solve_strip(2.0, b=0.025, nx=32, ny=2, steps=16, scheme="newton")
    assert abs(w4 - w16) / w16 < 1e-3
    assert abs(u4 - u16) / abs(u16) < 1e-3


def test_default_scheme_is_newton():
    w, u = solve_strip(1.0, b=0.1, nx=16, ny=2, steps=8, scheme=None)
    w_explicit, u_explicit = solve_strip(1.0, b=0.1, nx=16, ny=2, steps=8, scheme="newton")
    assert w == pytest.approx(w_explicit, abs=1e-12)
    assert u == pytest.approx(u_explicit, abs=1e-12)


def test_legacy_ul_scheme_keeps_its_recorded_reading():
    """``scheme: ul`` still reproduces the old number (w/L = 0.29075, 8x4 mesh, 32 steps)."""
    w, _ = solve_strip(1.0, b=0.1, nx=8, ny=4, steps=32, scheme="ul")
    assert w == pytest.approx(0.29075, abs=2e-5)


def test_legacy_ul_scheme_bias_is_pinned():
    """Known defect of ``scheme: ul``: 23 % from CalculiX on the pre-twisted beam, mesh independent.

    Pinned so that a change of the legacy loop shows up; if it is ever fixed, delete the scheme
    and this test instead of loosening the bound.
    """
    tmp = pytest.importorskip("tempfile").mkdtemp()
    from pathlib import Path

    ccx = _ccx_tip(Path(tmp), ccx_bin_or_skip(), "in-plane", 10.0, 24, 4)
    ours = _ours_tip("in-plane", 10.0, 24, 4, scheme="ul")
    gap = float(np.linalg.norm(ours - ccx)) / float(np.linalg.norm(ccx))
    assert 0.15 < gap < 0.30


def test_unknown_scheme_is_rejected():
    with pytest.raises(ValueError, match="scheme"):
        solve_strip(1.0, b=0.1, nx=8, ny=4, steps=4, scheme="bogus")
