"""Tests for the BEM aerodynamic data layer (aerodynamics.py).

Covers:
- PolarData interpolation and periodic wrapping
- AirfoilAero Reynolds-number selection
- BladeAero loading from WindIO YAML (IEA 15 MW)
- BladeAero property accessors (r, chord, twist)
"""


import numpy as np
import pytest

from aeroelast.models.blade.aerodynamics import (
    AirfoilAero,
    PolarData,
    load_blade_aero,
)

from tests.support.paths import DATA_DIR  # noqa: E402

IEA_YAML = str(DATA_DIR / "IEA-15-240-RWT.yaml")


# =====================================================================
# PolarData
# =====================================================================


class TestPolarData:
    """Unit tests for the PolarData dataclass."""

    @pytest.fixture
    def naca0012_polar(self):
        """Simple NACA 0012-like polar (symmetric)."""
        alpha = np.linspace(-np.pi, np.pi, 361)
        cl = 2 * np.pi * np.sin(alpha)  # thin-airfoil Cl(alpha) ≈ 2π sin(α)
        cd = 0.01 + 0.1 * np.sin(alpha) ** 2
        cm = np.zeros_like(alpha)
        return PolarData(alpha=alpha, cl=cl, cd=cd, cm=cm, re=1e6)

    def test_evaluate_at_zero(self, naca0012_polar):
        """Cl should be ~0 at alpha = 0 for a symmetric airfoil."""
        cl, cd, cm = naca0012_polar.evaluate(np.array([0.0]))
        assert abs(cl[0]) < 1e-3
        assert cd[0] > 0

    def test_evaluate_at_known_alpha(self, naca0012_polar):
        """The polar interpolates linearly between tabulated alpha nodes.

        This checks the interpolator itself (the midpoint of two adjacent nodes is
        their mean), not the formula that generated the table, so a formula error
        cannot pass and an interpolation error cannot hide.
        """
        alpha = naca0012_polar.alpha
        cl = naca0012_polar.cl
        i = int(np.argmin(np.abs(np.rad2deg(alpha) - 5.0)))
        a_mid = 0.5 * (alpha[i] + alpha[i + 1])
        cl_mid, _, _ = naca0012_polar.evaluate(np.array([a_mid]))
        np.testing.assert_allclose(cl_mid[0], 0.5 * (cl[i] + cl[i + 1]), atol=1e-12)

    def test_evaluate_periodic_wrapping(self, naca0012_polar):
        """Alpha outside [-π, π] should wrap correctly."""
        alpha_pos = np.deg2rad(10.0)
        alpha_wrapped = alpha_pos + 2 * np.pi  # 370° should == 10°
        cl_ref, _, _ = naca0012_polar.evaluate(np.array([alpha_pos]))
        cl_wrap, _, _ = naca0012_polar.evaluate(np.array([alpha_wrapped]))
        np.testing.assert_allclose(cl_wrap, cl_ref, atol=1e-6)

    def test_evaluate_vectorised(self, naca0012_polar):
        """evaluate() should handle arrays."""
        alphas = np.deg2rad(np.array([-10.0, 0.0, 5.0, 10.0]))
        cl, cd, cm = naca0012_polar.evaluate(alphas)
        assert cl.shape == (4,)
        assert cd.shape == (4,)
        # Cl should increase from -10° to +10°
        assert cl[0] < cl[1] < cl[3]


# =====================================================================
# AirfoilAero
# =====================================================================


class TestAirfoilAero:
    """Reynolds number selection tests."""

    @pytest.fixture
    def multi_re_airfoil(self):
        """Airfoil with polars at Re = 3e6, 6e6, 9e6."""
        alpha = np.linspace(-np.pi, np.pi, 37)
        polars = []
        for re in [3e6, 6e6, 9e6]:
            cl = 2 * np.pi * np.sin(alpha) * (1 + 0.1 * re / 9e6)
            polars.append(
                PolarData(
                    alpha=alpha,
                    cl=cl,
                    cd=np.full_like(alpha, 0.01),
                    cm=np.zeros_like(alpha),
                    re=re,
                )
            )
        return AirfoilAero(
            name="test_af",
            coordinates=np.zeros((10, 2)),
            relative_thickness=0.21,
            aerodynamic_center=0.25,
            polars=polars,
        )

    def test_get_polar_exact_match(self, multi_re_airfoil):
        """Exact Re match should return the right polar."""
        polar = multi_re_airfoil.get_polar(6e6)
        assert polar.re == 6e6

    def test_get_polar_closest(self, multi_re_airfoil):
        """Nearest-Re selection should work."""
        polar = multi_re_airfoil.get_polar(7e6)
        assert polar.re == 6e6  # 7e6 is closer to 6e6 than 9e6

    def test_get_polar_single(self):
        """Single polar should be returned regardless of Re."""
        alpha = np.linspace(-np.pi, np.pi, 37)
        polar = PolarData(alpha=alpha, cl=np.zeros(37), cd=np.zeros(37), cm=np.zeros(37), re=1e6)
        af = AirfoilAero(
            name="single",
            coordinates=np.zeros((5, 2)),
            relative_thickness=0.1,
            aerodynamic_center=0.25,
            polars=[polar],
        )
        result = af.get_polar(9e6)
        assert result.re == 1e6


# =====================================================================
# BladeAero — YAML loading (IEA 15 MW)
# =====================================================================


class TestBladeAeroYAML:
    """Integration tests that load polars from the IEA 15 MW WindIO YAML."""

    @pytest.fixture(scope="class")
    def blade_aero(self):
        """Load IEA 15 MW blade aero data once per class."""
        return load_blade_aero(IEA_YAML)

    def test_blade_length_positive(self, blade_aero):
        assert blade_aero.blade_length > 100.0  # IEA 15 MW ~ 117 m

    def test_hub_radius_positive(self, blade_aero):
        assert blade_aero.hub_radius > 0.0

    def test_rotor_radius(self, blade_aero):
        assert blade_aero.rotor_radius > blade_aero.hub_radius

    def test_n_blades(self, blade_aero):
        assert blade_aero.n_blades == 3

    def test_stations_ordered(self, blade_aero):
        """Stations should be root-to-tip ordered."""
        r = blade_aero.r
        assert len(r) > 5
        assert np.all(np.diff(r) >= 0)

    def test_chord_physical(self, blade_aero):
        """Chord should be positive and < 10 m for IEA 15 MW."""
        chord = blade_aero.chord
        assert np.all(chord > 0)
        assert np.all(chord < 10.0)

    def test_twist_range(self, blade_aero):
        """Twist (rad) should be within ±30° for a utility-scale blade."""
        twist = blade_aero.twist
        assert np.all(np.abs(twist) < np.deg2rad(30))

    def test_airfoils_have_polars(self, blade_aero):
        """Every airfoil should have at least one polar table."""
        for af in blade_aero.airfoils:
            assert len(af.polars) >= 1, f"Airfoil {af.name} has no polars"

    def test_polar_alpha_covers_range(self, blade_aero):
        """Polars should cover a wide alpha range (at least ±90°)."""
        af = blade_aero.airfoils[0]
        polar = af.polars[0]
        alpha_range = polar.alpha.max() - polar.alpha.min()
        assert alpha_range >= np.pi  # at least 180°

    def test_polar_cl_not_constant(self, blade_aero):
        """The mid-span airfoil has a physical lift curve, not a flat table."""
        # Skip the root cylinder (airfoils[0]) — pick a mid-span airfoil
        af = blade_aero.airfoils[len(blade_aero.airfoils) // 2]
        polar = af.polars[0]
        assert np.std(polar.cl) > 0.1

        # Lift-curve slope near zero incidence [per radian], and near-zero lift
        # at zero incidence.  A flat or sign-flipped table fails here.
        alpha_deg = np.rad2deg(polar.alpha)
        m = np.abs(alpha_deg) <= 5.0
        slope = float(np.polyfit(polar.alpha[m], polar.cl[m], 1)[0])
        cl0 = float(np.interp(0.0, polar.alpha, polar.cl))
        print(f"  lift-curve slope = {slope:.3f} /rad, Cl(0) = {cl0:.4f}")
        assert 3.0 < slope < 8.0, f"non-physical lift-curve slope {slope:.3f} /rad"
        # A cambered section lifts at zero incidence; a symmetric one is ~0.
        assert -0.5 < cl0 < 1.0, f"Cl at zero incidence is {cl0:.4f}"

    def test_r_property_shape(self, blade_aero):
        """r property should match number of stations."""
        assert blade_aero.r.shape == (len(blade_aero.stations),)

    def test_chord_property_shape(self, blade_aero):
        assert blade_aero.chord.shape == (len(blade_aero.stations),)

    def test_twist_property_shape(self, blade_aero):
        assert blade_aero.twist.shape == (len(blade_aero.stations),)


def test_viterna_matches_aerodyn_theory_manual():
    """_viterna_extrapolation implements the AeroDyn manual eqs [98]-[102].

    The AeroDyn Theory Manual (Moriarty & Hansen 2005, NREL/TP-500-36881,
    p. 22) prints the Viterna equations and attributes them to Viterna &
    Janetzke (1982, NASA TM-82944).  This hand-codes them independently, so a
    wrong ``A2``/``B2`` cannot hide behind the model difference of the
    NeuralFoil-vs-official parity test, and checks the manual's two named
    limits at 90 deg (CL = 0, CD = Cd_max).
    """
    from aeroelast.models.blade.aerodynamics import _viterna_extrapolation

    ar = 17.0
    cd_max = 1.11 + 0.018 * ar
    alpha_s = np.deg2rad(15.0)
    cl_s, cd_s = 1.70, 0.050

    # Monotone attached samples whose last point on each side is the matching
    # point the extrapolation must use.
    alpha_attach = np.deg2rad(np.linspace(-15.0, 15.0, 61))
    cl_attach = cl_s * np.sin(alpha_attach) / np.sin(alpha_s)
    cd_attach = np.full_like(alpha_attach, cd_s)
    cm_attach = np.zeros_like(alpha_attach)

    alpha_full, cl_full, cd_full, _ = _viterna_extrapolation(
        alpha_attach, cl_attach, cd_attach, cm_attach, ar=ar
    )

    sin_s, cos_s = np.sin(alpha_s), np.cos(alpha_s)
    a2 = (cl_s - cd_max * sin_s * cos_s) * sin_s / cos_s**2
    b2 = (cd_s - cd_max * sin_s**2) / cos_s
    for a_deg in (20.0, 30.0, 45.0, 60.0, 90.0):
        a = np.deg2rad(a_deg)
        cl_ref = cd_max / 2.0 * np.sin(2.0 * a) + a2 * np.cos(a) ** 2 / np.sin(a)
        cd_ref = cd_max * np.sin(a) ** 2 + b2 * np.cos(a)
        cl = float(np.interp(a, alpha_full, cl_full))
        cd = float(np.interp(a, alpha_full, cd_full))
        assert abs(cl - cl_ref) < 1e-10, f"alpha={a_deg}: CL {cl} != manual {cl_ref}"
        assert abs(cd - cd_ref) < 1e-10, f"alpha={a_deg}: CD {cd} != manual {cd_ref}"

    # The manual's two named limits at 90 deg: CL = 0 and CD = Cd_max.
    assert abs(float(np.interp(np.pi / 2, alpha_full, cl_full))) < 1e-10
    assert abs(float(np.interp(np.pi / 2, alpha_full, cd_full)) - cd_max) < 1e-10
