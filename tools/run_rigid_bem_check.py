"""Rigid (undeformed) BEM check — isolate the BEM model from the FSI coupling.

The FSI campaign mixes three things: the BEM model (CCBlade + NeuralFoil
polars), the yaw model, and the aeroelastic coupling (deformed radius +
delta twist from the chord direction).  This tool removes the coupling by
running the *undeformed* blade through the same ``BEMSolver`` configuration
the FSI participant builds, so the remaining discrepancy can be assigned:

  * rigid vs the IEA tabular ``Rotor Performance`` sheet  -> BEM/polar model
  * rigid vs the FSI campaign                             -> the coupling
  * rigid vs Ma 2025 / Zhou 2025 (beam FSI)               -> the total gap

Two sweeps:

``--yaw-sweep``  yaw 0..40 deg (step 5) at the S-6 operating point
                 (wind 10.59 m/s, pitch 0 deg, 7.55 rpm) — the FSI basis.
``--curve``      the V-06 wind list at yaw 0 with the IEA reference schedule
                 (tabulated pitch and rpm) — the power/thrust curve check.

No preCICE and no structural solve, so the whole thing runs in minutes.

Usage:
  python tools/run_rigid_bem_check.py --yaw-sweep --curve \
      --csv docs/validation_data/generated/rigid_bem_check.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

BLADE_FILE = REPO / "tests" / "IEA-15-240-RWT.yaml"
REF_CSV = REPO / "docs" / "validation_data" / "reference_iea15mw_rotor_performance.csv"

# S-6 / FSI yaw campaign operating point
WIND_YAW = 10.59
PITCH_YAW = 0.0
OMEGA_RAD_S = 0.7906341464750989          # 7.55 rpm, rated
R_TIP = 120.97                            # rotor radius [m]
RHO = 1.225

# BEM config mirrored from tests/IEA15MW/frontiersin_2025_yaw/fluid_yaw_0.yaml
BEM_CFG = dict(
    default_re=1.0e7,
    hub_radius=0.0,
    n_blades=3,
    viterna_ar=17.0,
    viterna_confidence_threshold=0.5,
)


def _interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            t = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


def _cp(power: float, wind: float) -> float:
    return power / (0.5 * RHO * math.pi * R_TIP**2 * wind**3)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--yaw-sweep", action="store_true")
    ap.add_argument("--curve", action="store_true")
    ap.add_argument("--neuralfoil-model", default="large")
    ap.add_argument("--csv", type=Path,
                    default=REPO / "docs" / "validation_data" / "generated" / "rigid_bem_check.csv")
    args = ap.parse_args()
    if not (args.yaw_sweep or args.curve):
        ap.error("select at least one of --yaw-sweep / --curve")

    from aeroelast.models.blade.aerodynamics import load_blade_aero
    from aeroelast.solvers.bem.engine import BEMSolver

    print(f"loading blade aerodynamics from {BLADE_FILE.name} "
          f"(neuralfoil={args.neuralfoil_model}) ...")
    blade = load_blade_aero(
        str(BLADE_FILE),
        neuralfoil_model=args.neuralfoil_model,
        **BEM_CFG,
    )
    omega_rpm = OMEGA_RAD_S * 60.0 / (2.0 * math.pi)
    print(f"rotor rpm = {omega_rpm:.4f}  (rated)")

    rows = []

    if args.yaw_sweep:
        print("\n=== rigid yaw sweep (wind 10.59, pitch 0, 7.55 rpm) ===")
        print(f"{'yaw':>4} {'P [MW]':>8} {'T [MN]':>8} {'Q [MNm]':>9} {'CP':>7} {'CT':>7}")
        fs = {}
        for yaw in range(0, 45, 5):
            solver = BEMSolver(blade, rho=RHO, mu=1.81e-5, precone=0.0, tilt=0.0,
                               yaw=float(yaw), hub_height=150.0, shear_exp=0.0)
            res = solver.compute(WIND_YAW, omega_rpm, PITCH_YAW, azimuth=0.0)
            P, T, Q = res.power / 1e6, res.thrust / 1e6, res.torque / 1e6
            fs[yaw] = (P, T)
            print(f"{yaw:>4} {P:8.3f} {T:8.3f} {Q:9.3f} {_cp(res.power, WIND_YAW):7.4f} "
                  f"{T * 1e6 / (0.5 * RHO * math.pi * R_TIP**2 * WIND_YAW**2):7.4f}")
            rows.append(["yaw", yaw, WIND_YAW, PITCH_YAW, omega_rpm, P, T,
                         _cp(res.power, WIND_YAW), ""])
        # effective cosine exponent of the power degradation
        p0 = fs[0][0]
        print("\ncos^n fit of P(gamma)/P(0):")
        for yaw in (10, 20, 30, 40):
            if yaw in fs:
                n = math.log(fs[yaw][0] / p0) / math.log(math.cos(math.radians(yaw)))
                print(f"  yaw {yaw:2d}: P/P0 = {fs[yaw][0]/p0:6.3f}  -> n = {n:5.2f}")

    if args.curve:
        ref = list(csv.DictReader(REF_CSV.open()))
        rw = [float(r["wind_mps"]) for r in ref]
        rp = [float(r["power_MW"]) for r in ref]
        rt = [float(r["thrust_MN"]) for r in ref]
        rr = [float(r["rpm"]) for r in ref]
        rpit = [float(r["pitch_deg"]) for r in ref]
        print("\n=== rigid power curve at yaw 0 (IEA reference schedule) ===")
        print(f"{'V':>5} {'pitch':>6} {'rpm':>6} | {'P_rig':>7} {'P_ref':>7} {'dP%':>7} | "
              f"{'T_rig':>6} {'T_ref':>6} {'dT%':>7}")
        solver0 = BEMSolver(blade, rho=RHO, mu=1.81e-5, precone=0.0, tilt=0.0,
                            yaw=0.0, hub_height=150.0, shear_exp=0.0)
        for V in (3, 5, 7, 8, 9, 10, 10.59, 11, 12, 13, 15, 17, 19, 21, 23, 25):
            rpm = _interp(V, rw, rr)
            pitch = _interp(V, rw, rpit)
            res = solver0.compute(V, rpm, pitch, azimuth=0.0)
            P, T = res.power / 1e6, res.thrust / 1e6
            Pr, Tr = _interp(V, rw, rp), _interp(V, rw, rt)
            print(f"{V:5.2f} {pitch:6.2f} {rpm:6.2f} | {P:7.3f} {Pr:7.3f} "
                  f"{100*(P-Pr)/Pr:+7.2f} | {T:6.3f} {Tr:6.3f} {100*(T-Tr)/Tr:+7.2f}")
            rows.append(["curve", 0.0, V, pitch, rpm, P, T, _cp(res.power, V), Pr])

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sweep", "yaw_deg", "wind_mps", "pitch_deg", "rpm",
                    "power_MW", "thrust_MN", "CP", "power_ref_MW"])
        w.writerows(rows)
    print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
