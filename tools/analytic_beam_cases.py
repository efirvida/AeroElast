"""Write the shared, versioned analytic beam deck once; refuse overwrites.

Usage: scripts/aeroenv.sh python tools/analytic_beam_cases.py
No shell solve is needed; all derived references are independent beam formulas.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.validation.benchmarks.analytic_beam_helpers import (  # noqa: E402
    E, HEIGHT, LENGTH, NU, RHO, THICKNESS, WIDTH,
    parked_reference, rotating_beam_ritz, section, shear_areas, timoshenko_reference,
)


def definition():
    area, ix, iy, j = section()
    g = E/(2*(1+NU))
    refs, timo = parked_reference(), timoshenko_reference()
    omega = .4 * 2*np.pi*refs["flap"][0]
    analytic = {}
    for family in ("flap", "edge"):
        modes = []
        for order, f in enumerate(refs[family]):
            record = {"order": order+1, "omega_0_EB_rad_s": 2*np.pi*f,
                      "omega_0_Timoshenko_rad_s": 2*np.pi*timo[family][order]}
            k1, k0 = rotating_beam_ritz(omega/(2*np.pi*f), order=order)
            record.update(K1=k1, K1_low_speed=k0,
                          omega_rotating_rad_s=float(np.sqrt((2*np.pi*f)**2 +
                              (k1-(family == "edge"))*omega**2)))
            modes.append(record)
        analytic[family] = modes
    analytic["torsion"] = {"order": 1, "omega_0_rad_s": 2*np.pi*refs["torsion"][0]}
    return {"version": 2, "units": "SI; flap=y, edge=x, span=z, rotation axis=y",
            "uniform_box": {
                "geometry": {"section": "four non-overlapping shell wall strips",
                             "w": WIDTH, "h": HEIGHT, "t": THICKNESS, "L": LENGTH,
                             "dimensions_reference": "wall midplanes", "r0": 0, "precone": 0},
                "material": {"E": E, "nu": NU, "rho": RHO},
                "Omega_rad_s": omega, "parked_Omega_rad_s": 0,
                "beam": {"EA": E*area, "EI_flap": E*ix, "EI_edge": E*iy, "GJ": g*j,
                         "kGA_flap": g*shear_areas()[0], "kGA_edge": g*shear_areas()[1],
                         "kappa_flap": shear_areas()[0]/area,
                         "kappa_edge": shear_areas()[1]/area,
                         "shear_method": "closed-cell thin-wall shear energy, integral(q/t ds)=0",
                         "m": RHO*area, "rotary_flap": RHO*ix, "rotary_edge": RHO*iy,
                         "rotary_polar": RHO*(ix+iy),
                         "elastic_centre_xyz": [0, 0, 0], "shear_centre_xyz": [0, 0, 0],
                         "mass_centre_xyz": [0, 0, 0]},
                "analytic": analytic,
                "formulas_sources": [
                    "EB omega_n=beta_n^2 sqrt(EI/(m L^4)); beta=1.8751040687,4.694091133; Blevins",
                    "torsion omega=pi/(2L) sqrt(GJ/(rho(Ix+Iy))); J=2w^2h^2t/(w+h), Bredt-Batho",
                    "I_flap=t(w h^2/2+h^3/6)+w t^3/6; edge swap w,h; m=rho 2(w+h)t",
                    "rotating EB y''''-(mu^2(1-x^2)y'/2)'=lambda*y, Legendre Ritz degree12",
                    "K1=(lambda(mu)-lambda(0))/mu^2, lead-lag K=K1-1; finite-speed, not constant",
                    "Wright, Smith, Thresher & Wang 1982, J Appl Mech 49:197-202; K1(0)~1.193",
                    "Timoshenko energy EI theta'^2+GAs(w'-theta)^2; mass m w^2+rho I theta^2"
                ]}}


if __name__ == "__main__":
    target = Path("logs/check/analytic_beam/cases.json")
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("x") as stream:
        json.dump(definition(), stream, indent=2)
        stream.write("\n")
    print(target)
