"""Planar large-deflection reference (inextensible elastica) for a cantilever with variable ``EI(s)``.

Dead loads only: a tip force ``P`` and an optional distributed weight ``w(s)``, both acting along
the fixed transverse direction ``y`` (positive along the load). Equations, with ``Q(s)`` the load
resultant beyond ``s``::

    theta' = M / EI(s),   M' = -cos(theta) Q(s),   x' = cos(theta),   y' = sin(theta)
    theta(0) = 0, x(0) = y(0) = 0, M(L) = 0

``Q`` is known in closed form from the loads, so the BVP has no integral unknowns. Solved with
``scipy.integrate.solve_bvp``. Reference checks live in ``tools/tests/test_blade_flap_elastica.py``:
the uniform-``EI`` tip-load elastica against the Bisshopp-Drucker (1945) tabulated values and the
small-load limit against the double-integrated Euler-Bernoulli deflection.

Usage: python tools/blade_flap_elastica_reference.py --tip-load 1e5 [--gravity]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "tools"))


def _cumulative_beyond(s: np.ndarray, w: np.ndarray) -> np.ndarray:
    """``int_s^L w dt`` on the grid ``s`` (trapezoid)."""
    seg = 0.5 * (w[1:] + w[:-1]) * np.diff(s)
    return np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])


def solve_elastica(
    s_tab: np.ndarray,
    ei_tab: np.ndarray,
    tip_load: float = 0.0,
    weight_tab: np.ndarray | None = None,
    n: int = 2001,
) -> dict[str, np.ndarray]:
    """Return ``s, theta, x, y`` of the dead-loaded elastica; ``y`` is positive along the load."""
    length = float(s_tab[-1])
    s_grid = np.linspace(0.0, length, 4001)
    w = np.zeros_like(s_grid) if weight_tab is None else np.interp(s_grid, s_tab, weight_tab)
    q_tab = _cumulative_beyond(s_grid, w) + tip_load

    def ei(s):
        return np.interp(s, s_tab, ei_tab)

    def rhs(s, u):
        theta, moment = u[0], u[1]
        return np.vstack(
            [
                moment / ei(s),
                -np.cos(theta) * np.interp(s, s_grid, q_tab),
                np.cos(theta),
                np.sin(theta),
            ]
        )

    def bc(ua, ub):
        return np.array([ua[0], ua[2], ua[3], ub[1]])

    s0 = np.linspace(0.0, length, n)
    guess = np.zeros((4, n))
    guess[2] = s0
    sol = solve_bvp(rhs, bc, s0, guess, tol=1e-9, max_nodes=200000)
    if not sol.success:
        raise RuntimeError(f"elastica BVP did not converge: {sol.message}")
    s_out = np.linspace(0.0, length, n)
    theta, _, x, y = sol.sol(s_out)
    return {"s": s_out, "theta": theta, "x": x, "y": y}


def linear_deflection(
    s_tab: np.ndarray,
    ei_tab: np.ndarray,
    tip_load: float = 0.0,
    weight_tab: np.ndarray | None = None,
) -> np.ndarray:
    """Euler-Bernoulli deflection by double integration of ``M(s) = int_s^L (t-s) q dt`` (small load)."""
    s = np.linspace(0.0, float(s_tab[-1]), 4001)
    w = np.zeros_like(s) if weight_tab is None else np.interp(s, s_tab, weight_tab)
    q_beyond = _cumulative_beyond(s, w) + tip_load
    m = np.concatenate(
        [
            [0.0],
            np.cumsum(0.5 * (q_beyond[1:] + q_beyond[:-1]) * np.diff(s)),
        ]
    )
    m = m[-1] - m  # M(s) = int_s^L Q(t) dt
    kappa = m / np.interp(s, s_tab, ei_tab)
    theta = np.concatenate([[0.0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * np.diff(s))])
    return np.concatenate([[0.0], np.cumsum(0.5 * (theta[1:] + theta[:-1]) * np.diff(s))])


def main() -> None:
    from beam_reference import load_elastodyn_blade

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tip-load", type=float, default=0.0, help="N, dead, along +y")
    ap.add_argument("--gravity", action="store_true", help="add the blade self-weight (9.81 m/s^2)")
    args = ap.parse_args()
    ed = load_elastodyn_blade(
        REPO_ROOT / "tests" / "IEA15MW" / "reference" / "IEA-15-240-RWT_ElastoDyn_blade.dat"
    )
    weight = ed.mass_dens * 9.81 if args.gravity else None
    sol = solve_elastica(ed.z, ed.ei_flap, args.tip_load, weight)
    lin = linear_deflection(ed.z, ed.ei_flap, args.tip_load, weight)
    print(f"tip y  elastica {sol['y'][-1]:+.5f} m  linear {lin[-1]:+.5f} m")
    print(f"tip x  shortening {sol['x'][-1] - ed.z[-1]:+.5f} m")


if __name__ == "__main__":
    main()
