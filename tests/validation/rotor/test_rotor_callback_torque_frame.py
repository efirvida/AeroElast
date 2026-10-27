"""C8 D3: the step callback's torques use the rotating-frame coordinates unrotated.

The mesh coordinates, nodal forces and displacements the corotational solver
holds are all rotating-frame quantities (Rust uses them unrotated,
``rotor_fsi.rs`` converged-window torque). The torque about the axis is
therefore ``sum (X + u - c) x F . axis`` at every azimuth theta; it must not be
modulated by theta.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.support.assertions import assert_residual_below

try:
    from aeroelast.solvers.fsi.corotational import CoordinateTransforms
    from aeroelast.solvers.fsi.rotor import LinearDynamicFSIRotorCorotationalSolver

    _HAS_ROTOR = True
except (ImportError, OSError):
    _HAS_ROTOR = False

pytestmark = pytest.mark.skipif(not _HAS_ROTOR, reason="rotor solver import unavailable")

AXIS = np.array([0.0, 1.0, 0.0])
X = np.array([[3.0, 0.0, 0.0], [5.0, 0.0, 1.0], [8.0, 0.5, -0.5]])
U = np.array([[0.0, 0.1, 0.0], [0.0, 0.2, 0.05], [0.1, 0.3, 0.0]])
F = np.array([[10.0, 40.0, 0.0], [-5.0, 60.0, 2.0], [7.0, 80.0, -3.0]])


def _reference_torque() -> float:
    return float(np.sum(np.cross(X + U, F), axis=0) @ AXIS)


def _solver():
    solver = object.__new__(LinearDynamicFSIRotorCorotationalSolver)
    solver._coord_transforms = CoordinateTransforms(AXIS, np.zeros(3))
    return solver


@pytest.mark.parametrize("theta_deg", [0.0, 90.0, 180.0])
def test_window_torques_not_modulated_by_theta(theta_deg):
    tau_ref = _reference_torque()
    zeros = np.zeros_like(F)
    out = _solver()._window_torques(X, U, F, X, U, zeros, zeros, F, np.deg2rad(theta_deg))
    tau_aero, tau_total = out[1], out[5]
    for name, tau in (("aero", tau_aero), ("total", tau_total)):
        assert_residual_below(
            abs(tau - tau_ref) / abs(tau_ref),
            tol=1e-12,
            kind="analytical",
            reference_name="sum (X + u) x F . axis on the rotating-frame data (as Rust)",
            what=f"{name} torque at theta = {theta_deg} deg",
        )


def test_window_torques_inertial_and_gravity_not_modulated_by_theta():
    tau_ref = _reference_torque()
    zeros = np.zeros_like(F)
    theta = np.pi / 2
    out = _solver()._window_torques(X, U, F, X, U, F, F, zeros, theta)
    for name, tau in (("inertial", out[2]), ("gravity", out[3])):
        assert_residual_below(
            abs(tau - tau_ref) / abs(tau_ref),
            tol=1e-12,
            kind="analytical",
            reference_name="sum (X + u) x F . axis on the rotating-frame data (as Rust)",
            what=f"{name} torque at theta = 90 deg",
        )
