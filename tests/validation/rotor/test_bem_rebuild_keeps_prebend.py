"""Rebuilding the BEM solver at the reference geometry must be the identity.

``BEMFSIParticipant._compute_forces`` uses the solver built at construction only when the
displacement is exactly zero; any other step goes through ``_rebuild_bem_solver``.  That
rebuild dropped the stations' ``prebend`` and ``sweep`` (L5 T18, odd/tasks/l5-one-way.md L37), so
every coupled run evaluated a straight blade: thrust -1.1 %, torque -1.9 % at rated on the
IEA-15.  With ``deformed_twist``/``deformed_radius`` off the rebuild receives the reference
radii and twist, and must reproduce the construction-time solver exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.solvers.bem import fsi_participant as fp  # noqa: E402

from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
V_RATED, RPM_RATED, PITCH_RATED = 10.59, 7.56, 0.0
BEM_KWARGS = {"hub_height": 150.0, "shear_exp": 0.0, "precone": -4.0, "tilt": -6.0}


@pytest.fixture(scope="module")
def stub():
    """A participant carrying only the state ``_rebuild_bem_solver`` reads (no mesh)."""
    if not AD_PRIMARY.exists():
        pytest.skip(f"AeroDyn reference not present: {AD_PRIMARY}")
    ba = build_blade_aero_from_aerodyn(AD_PRIMARY)
    p = fp.BEMFSIParticipant.__new__(fp.BEMFSIParticipant)
    p._blade_aero = ba
    p._ref_r = np.array([s.r for s in ba.stations])
    p._ref_twist = np.array([s.twist for s in ba.stations])
    p._bem_solver_kwargs = BEM_KWARGS
    return p


def test_rebuild_keeps_station_prebend_and_sweep(stub):
    ba = stub._blade_aero
    assert np.any(np.abs(ba.prebend) > 0.0), "fixture blade must carry prebend for this test to mean anything"
    _, rebuilt = stub._rebuild_bem_solver(stub._ref_r, stub._ref_twist)
    np.testing.assert_array_equal([s.prebend for s in rebuilt.stations], [s.prebend for s in ba.stations])
    np.testing.assert_array_equal([s.sweep for s in rebuilt.stations], [s.sweep for s in ba.stations])


def test_rebuild_at_reference_reproduces_construction_loads(stub, capsys):
    base = fp.BEMSolver(stub._blade_aero, **BEM_KWARGS).compute(V_RATED, RPM_RATED, PITCH_RATED)
    solver, _ = stub._rebuild_bem_solver(stub._ref_r, stub._ref_twist)
    again = solver.compute(V_RATED, RPM_RATED, PITCH_RATED)
    capsys.readouterr()
    np.testing.assert_allclose(again.Np, base.Np, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(again.Tp, base.Tp, rtol=1e-12, atol=1e-9)
