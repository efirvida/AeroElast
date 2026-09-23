"""V-01 — IEA 15 MW RWT: rotor aerodynamic performance at rated conditions.

Validates that the BEM solver reproduces the reference Cp, Ct, Thrust and
Torque from the IEA 15 MW definition at rated wind speed.

Reference: NREL/TP-5000-75698, Table 3-1; IEA-15-240-RWT_tabular.xlsx sheet
"Rotor Performance" (auto-generated from WISDEM quasi-static).

Tolerances:
    - Thrust (rotor total): ±5 % of reference
    - Torque (rotor total): ±5 % of reference
  - Cp_aero             : ±0.02 absolute
  - Ct                  : ±0.05 absolute
"""

from __future__ import annotations

import pytest

ccblade = pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")

from aeroelast.models.blade.aerodynamics import load_blade_aero
from aeroelast.solvers.bem.engine import BEMSolver

# ---------------------------------------------------------------------------
# IEA 15 MW rated-condition reference values
# Source: IEA-15-240-RWT_tabular.xlsx, sheet "Rotor Performance", row V=10.659 m/s
# [R3] IEA-15-240-RWT_tabular.xlsx
# ---------------------------------------------------------------------------
_RATED_WIND_SPEED = 10.659  # m/s   [R3]
_RATED_OMEGA_RPM = 7.518  # RPM   [R3]
_RATED_OMEGA = _RATED_OMEGA_RPM  # RPM (BEMSolver.compute expects RPM)

_REF_THRUST_ROTOR = 2.457e6  # N     [R3] — total rotor
_REF_TORQUE_ROTOR = 19.91e6  # N·m   [R3] — total rotor
_REF_CT = 0.7718  # [-]   [R3]
_REF_CP_AERO = 0.4618  # [-]   [R3]

_RHO = 1.225  # kg/m³ — air density at sea level
_ROTOR_RADIUS = 120.675  # m     [R3] — half of Ø 241.35 m

# Reference curve at selected wind speeds (V, Ct_ref) for curve-shape check
# Source: [R3] Rotor Performance sheet
_PERF_CURVE = [
    # (V m/s, Cp_aero,  Ct)
    (5.006, 0.4164, 0.7842),
    (7.159, 0.4616, 0.7783),
    (9.027, 0.4616, 0.7783),
    (10.659, 0.4618, 0.7718),
    (12.259, 0.3036, 0.3875),
    (15.471, 0.1511, 0.1796),
]


@pytest.fixture(scope="module")
def blade_aero(iea_blade_yaml):
    """IEA 15 MW blade aero model."""
    return load_blade_aero(iea_blade_yaml)


@pytest.fixture(scope="module")
def bem(blade_aero):
    return BEMSolver(blade_aero, rho=_RHO, mu=1.81206e-5)


class TestRatedConditions:
    """V-01a — rated-point validation."""

    @pytest.fixture(scope="class")
    def rated(self, bem):
        return bem.compute(v_inf=_RATED_WIND_SPEED, omega=_RATED_OMEGA, pitch=0.0)

    def test_thrust_within_5pct(self, rated):
        """Rotor thrust must be within ±5 % of IEA reference [R3]."""
        rel_err = abs(rated.thrust - _REF_THRUST_ROTOR) / _REF_THRUST_ROTOR
        assert rel_err < 0.05, (
            f"Thrust {rated.thrust / 1e3:.1f} kN vs ref {_REF_THRUST_ROTOR / 1e3:.1f} kN "
            f"(err={rel_err:.1%})"
        )

    def test_torque_within_5pct(self, rated):
        """Rotor torque must be within ±5 % of IEA reference [R3]."""
        rel_err = abs(rated.torque - _REF_TORQUE_ROTOR) / _REF_TORQUE_ROTOR
        assert rel_err < 0.05, (
            f"Torque {rated.torque / 1e6:.3f} MN·m vs ref {_REF_TORQUE_ROTOR / 1e6:.3f} MN·m "
            f"(err={rel_err:.1%})"
        )

    def test_cp_aero_within_tolerance(self, rated):
        """Cp_aero must be within ±0.02 of IEA reference [R3]."""
        assert rated.CP is not None, "BEMSolver did not return CP coefficient"
        assert abs(rated.CP - _REF_CP_AERO) < 0.02, (
            f"Cp_aero={rated.CP:.4f} vs ref={_REF_CP_AERO:.4f} "
            f"(diff={rated.CP - _REF_CP_AERO:+.4f})"
        )

    def test_ct_within_tolerance(self, rated):
        """Ct must be within ±0.05 of IEA reference [R3]."""
        assert rated.CT is not None, "BEMSolver did not return CT coefficient"
        assert abs(rated.CT - _REF_CT) < 0.05, (
            f"Ct={rated.CT:.4f} vs ref={_REF_CT:.4f} (diff={rated.CT - _REF_CT:+.4f})"
        )


class TestPowerCurveShape:
    """V-01b — verify Cp and Ct trend across wind speeds [R3]."""

    @pytest.mark.parametrize("v_inf, cp_ref, ct_ref", _PERF_CURVE)
    def test_cp_ct_at_wind_speed(self, bem, v_inf, cp_ref, ct_ref):
        """Cp and Ct must each be within ±0.05 of the WISDEM reference [R3]."""
        pitch_map = {
            5.006: 2.893,  # [R3] below rated: pitch follows minimum pitch schedule
            7.159: 0.0,
            9.027: 0.0,
            10.659: 0.0,
            12.259: 6.76,
            15.471: 12.19,
        }
        rpm_map = {
            5.006: 5.0,
            7.159: 5.10,
            9.027: 6.43,
            10.659: 7.518,
            12.259: 7.518,
            15.471: 7.518,
        }
        pitch = pitch_map.get(v_inf, 0.0)
        rpm = rpm_map.get(v_inf, 7.518)

        result = bem.compute(v_inf=v_inf, omega=rpm, pitch=pitch)

        assert result.CP is not None, "BEMSolver did not return CP"
        assert result.CT is not None, "BEMSolver did not return CT"

        assert abs(result.CP - cp_ref) < 0.05, (
            f"V={v_inf} m/s: Cp={result.CP:.4f} vs ref={cp_ref:.4f}"
        )
        assert abs(result.CT - ct_ref) < 0.05, (
            f"V={v_inf} m/s: Ct={result.CT:.4f} vs ref={ct_ref:.4f}"
        )
