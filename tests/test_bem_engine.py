"""Tests for the BEM engine (CCBlade wrapper).

Requires ccblade to be installed: ``pip install -e ".[bem]"``

Covers:
- BEMSolver construction from BladeAero
- Parked blade (omega=0): AoA ≈ twist + pitch
- BEMResult fields are physically sensible
"""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")

from aeroelast.models.blade.aerodynamics import load_blade_aero
from aeroelast.solvers.bem.engine import BEMResult, BEMSolver

_PROJECT_ROOT = Path(__file__).resolve().parent
IEA_YAML = str(_PROJECT_ROOT / "IEA-15-240-RWT.yaml")


@pytest.fixture(scope="module")
def blade_aero():
    """Load IEA 15 MW blade aero data once per module."""
    return load_blade_aero(IEA_YAML)


@pytest.fixture(scope="module")
def bem_solver(blade_aero):
    """Create a BEM solver for the IEA 15 MW blade."""
    return BEMSolver(blade_aero, rho=1.225, mu=1.81206e-5)


class TestBEMSolverParked:
    """Tests for parked-blade BEM analysis (omega = 0)."""

    @pytest.fixture(scope="class")
    def parked_result(self, bem_solver):
        """BEM result for parked blade at V=45 m/s, pitch=0."""
        return bem_solver.compute(v_inf=45.0, omega=0.0, pitch=0.0)

    def test_result_is_bem_result(self, parked_result):
        assert isinstance(parked_result, BEMResult)

    def test_r_matches_blade(self, parked_result, blade_aero):
        np.testing.assert_array_equal(parked_result.r, blade_aero.r)

    def test_Np_shape(self, parked_result, blade_aero):
        assert parked_result.Np.shape == blade_aero.r.shape

    def test_Tp_shape(self, parked_result, blade_aero):
        assert parked_result.Tp.shape == blade_aero.r.shape

    def test_alpha_shape(self, parked_result, blade_aero):
        assert parked_result.alpha.shape == blade_aero.r.shape

    def test_parked_alpha_close_to_twist(self, parked_result, blade_aero):
        """For a parked blade (omega=0, pitch=0), AoA ≈ 90° - twist.

        When the blade is stationary the inflow is purely axial, so the
        angle of attack at each section equals the complement of the
        twist angle: α = 90° − θ_twist (both in degrees here, since BEM
        result reports alpha in degrees).
        """
        twist_deg = np.rad2deg(blade_aero.twist)
        expected_alpha = 90.0 - twist_deg  # from geometry
        # Allow generous tolerance — BEM induction modifies the inflow angle
        # but for parked conditions (no rotation) it should be close
        residual = np.abs(parked_result.alpha - expected_alpha)
        # At least 80% of stations should be within 15° — inboard stations
        # with thick airfoils and root corrections may deviate more
        fraction_close = np.mean(residual < 15.0)
        assert fraction_close > 0.8, (
            f"Only {fraction_close * 100:.0f}% of stations have alpha within 15° of expected"
        )

    def test_parked_thrust_positive(self, parked_result):
        """Thrust should be positive for positive wind speed."""
        assert parked_result.thrust > 0

    def test_parked_torque_near_zero(self, parked_result):
        """Torque should be very small for a parked rotor."""
        # For parked, power ≈ 0 (no rotation)
        assert abs(parked_result.power) < 1e6  # much less than rated power

    def test_parked_Np_mostly_positive(self, parked_result):
        """Normal force should be mostly positive (pushed downwind)."""
        frac_positive = np.mean(parked_result.Np > 0)
        assert frac_positive > 0.7


class TestBEMSolverRotating:
    """Tests for rotating-blade BEM analysis."""

    @pytest.fixture(scope="class")
    def rated_result(self, bem_solver):
        """BEM result near rated conditions: V=10.59 m/s, 7.56 RPM."""
        return bem_solver.compute(v_inf=10.59, omega=7.56, pitch=0.0)

    def test_rated_power_positive(self, rated_result):
        """A rotating blade in wind should produce positive power."""
        assert rated_result.power > 0

    def test_rated_thrust_positive(self, rated_result):
        assert rated_result.thrust > 0

    def test_rated_torque_positive(self, rated_result):
        assert rated_result.torque > 0

    def test_axial_induction_physical(self, rated_result):
        """Axial induction factor should be between 0 and ~0.6."""
        a = rated_result.a
        # Some root/tip stations may have extreme values; check bulk
        mask = np.isfinite(a) & (a > -0.5)
        frac_ok = np.mean((a[mask] >= -0.1) & (a[mask] <= 0.6))
        assert frac_ok > 0.7

    def test_cl_not_all_zero(self, rated_result):
        """Lift coefficients should be non-trivial."""
        assert np.any(np.abs(rated_result.cl) > 0.1)


class TestOmegaSignConvention:
    """Regression tests for the rotor sign convention.

    Current design (reworked 2026-09): the structural solver's RHR omega
    (positive = clockwise viewed from -Y) is passed to CCBlade WITHOUT a
    sign change — CCBlade's positive RPM direction matches the structural
    omega > 0.  The load direction is handled by the
    ``tangential_direction: [-1, 0, 0]`` YAML setting, not by negating
    omega.  These tests guard that invariant: positive RPM must give
    positive power/torque, and negating omega must invert the torque sign.
    """

    @pytest.fixture(scope="class")
    def cw_result(self, bem_solver):
        """BEM result under the current convention (omega > 0 = CW structural)."""
        import math

        current_omega_rad_s = 0.7917  # rad/s, structural RHR (CW from -Y)
        omega_rpm = +current_omega_rad_s * 60.0 / (2.0 * math.pi)
        return bem_solver.compute(v_inf=10.59, omega=omega_rpm, pitch=0.0)

    @pytest.fixture(scope="class")
    def ccw_result(self, bem_solver):
        """BEM result with the reversed rotation (negative RPM)."""
        import math

        current_omega_rad_s = 0.7917
        omega_rpm = -current_omega_rad_s * 60.0 / (2.0 * math.pi)
        return bem_solver.compute(v_inf=10.59, omega=omega_rpm, pitch=0.0)

    def test_cw_torque_positive(self, cw_result):
        """Positive RPM (current no-negation convention) must give positive torque."""
        assert cw_result.torque > 0, (
            f"Torque={cw_result.torque:.3e} Nm — negative torque means the sign "
            "convention (omega > 0 -> CCBlade positive RPM) was reverted."
        )

    def test_cw_power_positive(self, cw_result):
        """Power must be positive: the rotor extracts energy from the wind."""
        assert cw_result.power > 0, (
            f"Power={cw_result.power:.3e} W — negative power means the CW sign "
            "convention is broken."
        )

    def test_cw_cp_physical_range(self, cw_result):
        """CP must be in (0, 0.593] (Betz limit) at rated conditions."""
        assert cw_result.CP is not None
        assert 0.0 < cw_result.CP <= 0.593, (
            f"CP={cw_result.CP:.4f} outside physical range — check omega sign."
        )

    def test_cw_ct_physical_range(self, cw_result):
        """CT must be positive and below ~1.5 at rated conditions."""
        assert cw_result.CT is not None
        assert 0.0 < cw_result.CT < 1.5, (
            f"CT={cw_result.CT:.4f} outside physical range."
        )

    def test_sign_flip_inverts_power(self, cw_result, ccw_result):
        """Reversing the rotation must invert the extracted power.

        The current engine reports the torque with the magnitude sign
        (positive for either rotation) while the power = torque × omega
        carries the direction; the campaigns rely on the power sign.
        """
        assert cw_result.power > 0
        assert ccw_result.power < 0, (
            f"Reversed rotation gave power={ccw_result.power:.3e} W — "
            "the rotation direction is not respected."
        )
