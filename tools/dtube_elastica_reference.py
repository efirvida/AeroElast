"""Independent large-deflection reference for the D-Tube cantilever under a vertical dead load.

Inextensible planar elastica of the cantilever (``dM/ds = -V cos(theta)``, ``theta' = M/EI(kappa)``)
solved with ``scipy.integrate.solve_bvp``. Two section laws:

* ``rigid``: constant ``EI = E I_c`` (no ovalization). Gives the beam answer, which stiffens slightly.
* ``brazier``: ``EI`` is reduced by the secant factor ``M_B(kappa) / (E I_B kappa)`` of the Brazier
  law of a CIRCULAR tube of the same radius and wall (energy minimum over a ``zeta cos 2theta``
  ovalization). The D-Tube web closes the arc against ovalization, so this factor is an UPPER bound
  of the softening the real section can show: it tests whether ovalization can explain a gap at all.

Usage: python tools/dtube_elastica_reference.py [--scale 0.5]
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.integrate import solve_bvp

from tests.validation.benchmarks import test_dtube_camarena_bending as dtube


def brazier_secant_factor(kappa: float, radius: float, thickness: float, nu: float) -> float:
    """``M(kappa) / (E I kappa)`` of a circular tube with a cos(2 theta) ovalization (Brazier)."""
    if kappa == 0.0:
        return 1.0
    # U = 1/2 E t k^2 [pi R^3 + pi R^2 z + pi/2 R z^2] + 1/2 E t^3/(12(1-nu^2)) 9 z^2/R^4 pi R
    bend = 9.0 * math.pi * thickness**3 / (12.0 * (1.0 - nu**2) * radius**3)
    zeta = (
        -0.5
        * math.pi
        * thickness
        * kappa**2
        * radius**2
        / (0.5 * math.pi * thickness * kappa**2 * radius + bend)
    )
    moment_over_e = (
        thickness
        * kappa
        * (math.pi * radius**3 + math.pi * radius**2 * zeta + 0.5 * math.pi * radius * zeta**2)
    )
    return moment_over_e / (math.pi * radius**3 * thickness * kappa)


def tip_deflection(scale: float, section: str = "rigid", n: int = 401) -> float:
    """Tip y of the dead-loaded elastica; ``scale`` is the fraction of the 63 kPa resultant."""
    w, length = scale * dtube.W_RES, dtube.L
    ei0 = dtube.E * dtube.I_C

    def stiffness(kappa: float) -> float:
        if section == "rigid":
            return ei0
        return ei0 * brazier_secant_factor(abs(kappa), dtube.R, dtube.T, dtube.NU)

    def rhs(s, u):
        theta, moment, _, _ = u
        # solve_bvp is vectorised: kappa solves kappa * EI(kappa) = M by fixed-point iteration
        kappa = moment / ei0
        for _ in range(6):
            kappa = np.array([m / stiffness(k) for m, k in zip(moment, kappa, strict=True)])
        return np.vstack([kappa, -w * (length - s) * np.cos(theta), np.cos(theta), np.sin(theta)])

    s = np.linspace(0.0, length, n)
    guess = np.zeros((4, n))
    guess[3] = w * length**4 / (8.0 * ei0) * (s / length) ** 2
    sol = solve_bvp(rhs, lambda a, b: [a[0], b[1], a[2], a[3]], s, guess, tol=1e-7, max_nodes=20000)
    if sol.status != 0:
        raise RuntimeError(f"solve_bvp did not converge: {sol.message}")
    return float(sol.y[3, -1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", type=float, default=0.5)
    args = parser.parse_args()
    linear = args.scale * dtube.DELTA_REF
    for section in ("rigid", "brazier"):
        tip = tip_deflection(args.scale, section)
        print(
            f"{section:8s} tip y = {tip:.4f} m  ({tip / linear - 1:+.2%} vs linear {linear:.4f} m)"
        )


if __name__ == "__main__":
    main()
