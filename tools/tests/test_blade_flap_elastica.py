import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from tools.blade_flap_elastica_reference import linear_deflection, solve_elastica

LENGTH, EI0 = 2.0, 5.0e3


def _uniform(n=60):
    s = np.linspace(0.0, LENGTH, n)
    return s, np.full(n, EI0)


def _shooting_tip(alpha_load):
    """Independent path: shoot on the root curvature of a uniform-EI cantilever (tip load)."""

    def tip_moment(m0):
        def f(s, u):
            return [u[1] / EI0, -np.cos(u[0]) * alpha_load, np.cos(u[0]), np.sin(u[0])]

        return solve_ivp(f, [0, LENGTH], [0, m0, 0, 0], rtol=1e-11, atol=1e-13).y[:, -1]

    m0 = brentq(lambda m: tip_moment(m)[1], 0.0, 4.0 * alpha_load * LENGTH)
    return tip_moment(m0)


def test_small_load_recovers_double_integrated_euler_bernoulli_with_variable_ei():
    s = np.linspace(0.0, LENGTH, 40)
    ei = EI0 * (1.0 - 0.6 * s / LENGTH)
    w = 3.0 * (1.0 - 0.3 * s / LENGTH)
    sol = solve_elastica(s, ei, tip_load=1e-3, weight_tab=w * 1e-3)
    lin = linear_deflection(s, ei, tip_load=1e-3, weight_tab=w * 1e-3)
    assert sol["y"][-1] == pytest.approx(lin[-1], rel=2e-3)


def test_uniform_ei_tip_load_matches_independent_shooting():
    load = 1.0 * EI0 / LENGTH**2  # PL^2/EI = 1
    s, ei = _uniform()
    sol = solve_elastica(s, ei, tip_load=load)
    tip = _shooting_tip(load)
    assert sol["x"][-1] == pytest.approx(tip[2], rel=1e-6)
    assert sol["y"][-1] == pytest.approx(tip[3], rel=1e-6)


def test_uniform_ei_matches_published_bisshopp_drucker_value():
    s, ei = _uniform()
    sol = solve_elastica(s, ei, tip_load=EI0 / LENGTH**2)
    assert sol["y"][-1] / LENGTH == pytest.approx(0.3017, abs=2e-4)
