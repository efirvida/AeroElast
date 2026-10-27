"""The deformed-geometry BEM must not hand CCBlade a station inside the hub.

``BEMFSIParticipant._compute_deformed_geometry`` returns ``r_def = _ref_r + u_span - shortening``.
On the rated IEA-15 one-way solve the root station comes back at ``Rhub - 8e-6 m`` (the
second-order span shortening exceeds the root's first-order stretch).  CCBlade's hub-loss
factor is ``exp(-B/2 (r - Rhub) / (Rhub sin phi))``; with ``r < Rhub`` the exponent flips
sign, the residual is NaN at the lower phi bracket, ``brentq`` raises, CCBlade prints
``NaNs at 0/50: 1e-06 0.0 1.5707963267948966`` and zeroes the root strip's load.  The
rigid-geometry root strip (``r = Rhub``) carries a finite load, so the deformed one must too.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.solvers.bem.fsi_participant import BEMFSIParticipant  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
V_RATED, RPM_RATED, PITCH_RATED = 10.59, 7.56, 0.0
BEM_KWARGS = {"hub_height": 150.0, "shear_exp": 0.0}
#: Measured on the rated one-way solve (job 11625411): root station inside the hub.
ROOT_BELOW_HUB_M = 8.0e-6


@pytest.fixture(scope="module")
def stub():
    """A participant carrying only the state ``_rebuild_bem_solver`` reads (no mesh)."""
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    ba = build_blade_aero_from_aerodyn(AD_PRIMARY)
    p = BEMFSIParticipant.__new__(BEMFSIParticipant)
    p._blade_aero = ba
    p._ref_r = np.array([s.r for s in ba.stations])
    p._bem_solver_kwargs = BEM_KWARGS
    return p


def _loads(p, r, capsys):
    solver, _ = p._rebuild_bem_solver(r, np.array([s.twist for s in p._blade_aero.stations]))
    capsys.readouterr()
    res = solver.compute(V_RATED, RPM_RATED, PITCH_RATED)
    return res, capsys.readouterr().out


def test_rigid_root_station_is_finite(stub, capsys):
    res, out = _loads(stub, stub._ref_r.copy(), capsys)
    assert "NaNs" not in out
    assert np.isfinite(res.Np).all() and res.Np[0] > 0.0


def test_root_station_slightly_inside_hub_keeps_its_load(stub, capsys):
    rigid, _ = _loads(stub, stub._ref_r.copy(), capsys)
    r = stub._ref_r.copy()
    r[0] -= ROOT_BELOW_HUB_M
    res, out = _loads(stub, r, capsys)
    assert "NaNs" not in out
    assert np.isfinite(res.Np).all() and np.isfinite(res.Tp).all()
    np.testing.assert_allclose(res.Np[0], rigid.Np[0], rtol=1e-3)
    np.testing.assert_allclose(res.Tp[0], rigid.Tp[0], rtol=1e-3)
