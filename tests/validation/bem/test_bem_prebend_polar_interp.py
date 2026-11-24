"""BEMSolver model inputs added by L4 T7 (odd/tasks/l4-bem-rigid-closure.md, L8-L10).

Two rules, each tested through the public ``BEMSolver`` interface:

* **Exact polar tables.**  Cl/Cd/Cm evaluated at a table node equal the table
  value (CCBlade's own ``CCAirfoil`` smooths Cl with ``s = 0.1``, which moves
  Cl by 0.007-0.013 outboard and is not what AeroDyn does).
* **Blade prebend.**  A blade that carries out-of-plane prebend
  (``AeroStation.prebend``, negative = upwind, the windIO ``reference_axis.x``
  and AeroDyn ``BlCrvAC`` convention) is passed to CCBlade as ``precurve``
  by default, with ``use_prebend=False`` as the straight-blade switch.  CCBlade's
  ``precurve`` is positive upwind for negative precone (``x_az = -r sin(precone) +
  precurve cos(precone)``, ``ccblade/src/*.f90`` ``defineCurvature``), so
  ``precurve = -prebend``.  The independent confirmation against AeroDyn is
  ``test_bem_openfast_parity.py::test_prebend_outboard_alpha_matches_aerodyn``.
"""

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")

from aeroelast.models.blade.aerodynamics import (
    AeroStation,
    AirfoilAero,
    BladeAero,
    PolarData,
    load_blade_aero,
)
from aeroelast.solvers.bem.engine import BEMSolver

from tests.support.paths import V2_BLADE_YAML  # noqa: E402

IEA_YAML = V2_BLADE_YAML  # windIO 2.x official file (gmsh-official-mesher T14)


def _wiggly_airfoil() -> AirfoilAero:
    """A polar with deliberate node-to-node noise: a smoothing spline cannot reproduce it."""
    alpha = np.deg2rad(np.arange(-180.0, 181.0, 2.0))
    n = alpha.size
    sign = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    cl = 2.0 * np.pi * np.sin(alpha) * np.cos(alpha) + 0.03 * sign
    cd = 0.01 + 1.2 * np.sin(alpha) ** 2 + 0.004 * sign
    cm = -0.05 + 0.01 * sign
    return AirfoilAero(
        name="wiggly",
        coordinates=np.empty((0, 2)),
        relative_thickness=0.2,
        aerodynamic_center=0.25,
        polars=[PolarData(alpha=alpha, cl=cl, cd=cd, cm=cm, re=1.0e6)],
    )


def _blade(prebend_tip: float = 0.0, airfoil: AirfoilAero | None = None) -> BladeAero:
    af = airfoil or _wiggly_airfoil()
    r = np.linspace(4.0, 60.0, 12)
    eta = (r - r[0]) / (r[-1] - r[0])
    stations = [
        AeroStation(
            span_fraction=float(e),
            r=float(ri),
            chord=float(3.0 - 2.0 * e),
            twist=float(np.deg2rad(12.0 * (1.0 - e))),
            pitch_axis=0.25,
            airfoil=af,
            prebend=float(prebend_tip * e**2),
        )
        for e, ri in zip(eta, r, strict=True)
    ]
    return BladeAero(
        airfoils=[af],
        stations=stations,
        blade_length=float(r[-1] - r[0]),
        hub_radius=float(r[0]),
        rotor_radius=float(r[-1]),
        n_blades=3,
    )


def _solver(blade: BladeAero, **kw) -> BEMSolver:
    return BEMSolver(blade, precone=-4.0, tilt=-6.0, shear_exp=0.0, hub_height=100.0, **kw)


# ---------------------------------------------------------------- polars


def test_polar_evaluation_equals_table_at_every_node():
    blade = _blade()
    solver = _solver(blade)
    pol = blade.airfoils[0].polars[0]
    cl, cd, cm = solver.polar(0, pol.alpha, re=1.0e6)
    np.testing.assert_allclose(cl, pol.cl, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(cd, pol.cd, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(cm, pol.cm, rtol=0.0, atol=1e-12)


def test_smoothed_polars_remain_available_and_differ_from_the_table():
    blade = _blade()
    pol = blade.airfoils[0].polars[0]
    smoothed = _solver(blade, polar_interp="smoothed")
    cl, _, _ = smoothed.polar(0, pol.alpha, re=1.0e6)
    assert np.max(np.abs(cl - pol.cl)) > 1e-3  # the old behaviour, opt-in only


def test_solution_cl_is_the_table_value_at_the_solved_alpha():
    """End to end: the Cl the BEM reports is the tabulated Cl at its own alpha."""
    blade = _blade()
    pol = blade.airfoils[0].polars[0]
    res = _solver(blade, polar_interp="linear").compute(9.0, 8.0, 0.0)
    cl_table = np.interp(np.deg2rad(res.alpha), pol.alpha, pol.cl)
    np.testing.assert_allclose(res.cl, cl_table, rtol=0.0, atol=1e-6)


# ---------------------------------------------------------------- prebend


def test_prebend_is_passed_to_ccblade_as_precurve_with_flipped_sign():
    blade = _blade(prebend_tip=-4.0)  # upwind
    solver = _solver(blade)
    np.testing.assert_allclose(solver.precurve, -blade.prebend)
    assert solver.precurve[-1] == pytest.approx(4.0)
    assert solver.rotor.precurveTip == pytest.approx(4.0)


def test_prebend_changes_outboard_inflow_and_loads():
    blade = _blade(prebend_tip=-4.0)
    with_pb = _solver(blade).compute(9.0, 8.0, 0.0)
    straight = _solver(blade, use_prebend=False).compute(9.0, 8.0, 0.0)
    outboard = slice(-4, None)
    # Outboard relative inflow changes (local cone), root stays put.
    assert np.max(np.abs(with_pb.W[outboard] - straight.W[outboard])) > 0.02
    assert abs(with_pb.W[0] - straight.W[0]) < 0.1 * np.max(np.abs(with_pb.W - straight.W))
    assert with_pb.thrust != pytest.approx(straight.thrust, rel=1e-4)


def test_use_prebend_false_equals_a_blade_without_prebend():
    straight_blade = _blade(prebend_tip=0.0)
    pb_blade = _blade(prebend_tip=-4.0)
    a = _solver(pb_blade, use_prebend=False).compute(9.0, 8.0, 0.0)
    b = _solver(straight_blade).compute(9.0, 8.0, 0.0)
    assert a.thrust == b.thrust
    assert a.torque == b.torque
    np.testing.assert_array_equal(a.alpha, b.alpha)


def test_blade_without_prebend_is_unchanged_by_the_default():
    blade = _blade(prebend_tip=0.0)
    a = _solver(blade).compute(9.0, 8.0, 0.0)
    b = _solver(blade, use_prebend=False).compute(9.0, 8.0, 0.0)
    assert a.thrust == b.thrust and a.torque == b.torque


def test_official_yaml_blade_carries_prebend_into_the_solver():
    blade = load_blade_aero(str(IEA_YAML))
    assert blade.prebend[-1] == pytest.approx(-4.0)
    solver = BEMSolver(blade)
    assert solver.rotor.precurveTip == pytest.approx(4.0)
    assert BEMSolver(blade, use_prebend=False).rotor.precurveTip == 0.0
