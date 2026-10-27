"""C8 D2: a driving aerodynamic torque reports positive power, CP and CQ.

Since T7c the tangential load Tp points to the leading edge (+x) and drives the
rotor with ``tau_aero > 0`` about ``+axis`` (see
``tests/validation/bem/test_force_projection_tangential_sense.py``). Power
extracted from the wind is ``tau_aero * omega`` > 0 at omega > 0.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.support.assertions import assert_residual_below

try:
    from aeroelast.solvers.fsi.rotor import LinearDynamicFSIRotorCorotationalSolver

    _HAS_ROTOR = True
except (ImportError, OSError):
    _HAS_ROTOR = False

pytestmark = pytest.mark.skipif(not _HAS_ROTOR, reason="rotor solver import unavailable")

RHO, V, R = 1.225, 10.0, 60.0
TAU, OMEGA = 2.0e6, 1.0


def _solver():
    solver = object.__new__(LinearDynamicFSIRotorCorotationalSolver)
    solver._fluid_density = RHO
    solver._flow_velocity = V
    return solver


def test_cp_cq_positive_for_driving_torque():
    _, cp, cq, _ = _solver()._compute_performance_coefficients(
        thrust=1.0e5, torque_aero=TAU, omega=OMEGA, radius=R
    )
    area = np.pi * R * R
    q_dyn = 0.5 * RHO * V * V
    for name, val, ref in (
        ("CP", cp, TAU * OMEGA / (q_dyn * area * V)),
        ("CQ", cq, TAU / (q_dyn * area * R)),
    ):
        assert val > 0.0, f"{name} = {val}: a driving torque must report a positive coefficient"
        assert_residual_below(
            abs(val - ref) / ref,
            tol=1e-12,
            kind="analytical",
            reference_name="tau * omega / (q A V) and tau / (q A R)",
            what=name,
        )


def test_aerodynamic_power_positive_for_driving_torque():
    power_aero, power_total = _solver()._aerodynamic_powers(TAU, 0.75 * TAU, OMEGA)
    assert_residual_below(
        abs(power_aero - TAU * OMEGA) / (TAU * OMEGA),
        tol=1e-12,
        kind="analytical",
        reference_name="tau_aero * omega",
        what="aerodynamic power",
    )
    assert_residual_below(
        abs(power_total - 0.75 * TAU * OMEGA) / (TAU * OMEGA),
        tol=1e-12,
        kind="analytical",
        reference_name="tau_total * omega",
        what="total power",
    )
