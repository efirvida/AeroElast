"""Which frame is Zhou et al. 2025's Fig. 11 in, and does the pair hold together?

Issue #15 (roadmap item 3) states that our tangential load integral runs `+31%` above
Zhou's while the normal force and the pitching moment match to 1%, and that the
difference "has not been attributed: neither to the polar used, nor to the tangential
induction, nor to the definition of the tangential direction in their frame". This tool
runs the third of those, on the production BEM path, from what they publish.

The two frames
--------------
A section's aerodynamic force can be split two ways, and they differ by the local twist
angle `theta` (chord angle measured from the rotor plane):

    rotor plane (what ccblade emits, and what `BEMSolver` returns):
        Np = qc (CL cos(phi) + CD sin(phi))
        Tp = qc (CL sin(phi) - CD cos(phi))          phi = alpha + theta
    section / chord frame (what AeroDyn publishes):
        Np = qc (CL cos(alpha) + CD sin(alpha))
        Tp = qc (CL sin(alpha) - CD cos(alpha))

Both give the same resultant magnitude, so the frames are told apart by the resultant's
*direction*, `psi := atan2(Tp, Np)`. With `gamma := atan2(CD, CL)`:

    section frame:    psi = alpha - gamma      =>  alpha_implied = psi + gamma
    rotor plane:      psi = phi   - gamma      =>  alpha_implied = psi + gamma - theta

Zhou publish the pair (`Fig. 11` `Np`, `Tp`) and the angle of attack (`Fig. 10`), so the
two hypotheses are separable **without any of their unpublished inputs** (pitch, section
stiffness, torsion definition - none of which enter: `theta` is the official deck's own
twist, and only the rotor-plane reading needs it).

The other two terms of the issue's list
---------------------------------------
* **polar `Cd`/`Cl` projection**: carried explicitly as `gamma`, evaluated from the
  official AeroDyn polar at their own published `alpha`. It enters both readings
  identically, so it cannot by itself explain a frame-shaped difference.
* **tangential induction `a'`**: not published, and not needed for the frame question -
  the reading that fits does so without it. It is what is left once the frame is fixed.

Quadrature, stated before it is used
------------------------------------
The digitized curves are read through a trapezoid on their own 18 points. That is not
the paper's quadrature, so each curve's implied rotor thrust and torque are checked
against the paper's own Table 6 numbers before any integral claim is made.

Run it with the repository root on the path (the tool's own directory is what a bare
`python tools/...` puts on `sys.path`)::

    PYTHONPATH=$PWD python tools/diagnose_zhou_tp_frame.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from aeroelast.solvers.bem.engine import BEMSolver  # noqa: E402
from tests.support.openfast_bem import build_blade_aero_from_aerodyn  # noqa: E402
from tests.support.paths import DATA_DIR  # noqa: E402

#: Zhou et al. 2025, *Energy* 336:138488, rated condition.
V_RATED = 10.59
RPM_ZHOU = 7.55
PITCH_RATED = 0.0

#: Table 6, fixed condition: (power [MW], thrust [MN]) for the rigid and flexible blade.
TABLE6_RIGID = (16.11, 2.53)
TABLE6_FLEXIBLE = (14.76, 2.20)

#: The reference's own digitization uncertainty, printed rather than assumed away.
FIG11_TP_BAND_KN_M = 0.05  # +/- 3-5 % of level, per the traceability card in the report
FIG10_AOA_BAND_DEG = 0.5

#: The blade core, where Fig. 10's angle of attack is well resolved and the polar is not
#: being extrapolated past its table. The digitized tip stations are excluded separately.
CORE_LO = 0.26
CORE_HI = 0.80

AD_PRIMARY = DATA_DIR / "reference" / "iea15mw_openfast" / "case" / "IEA-15-240-RWT_AeroDyn15.dat"
ZHOU_LOADS = _REPO_ROOT / "docs" / "validation_data" / "zhou_2025_fig11_loads.csv"
ZHOU_AOA = _REPO_ROOT / "docs" / "validation_data" / "zhou_2025_fig10_aoa.csv"


def _polar_at(aero, index: int, alpha_deg: float) -> tuple[float, float, float]:
    """The official polar at one station, interpolated at ``alpha_deg`` (as of #12)."""
    polar = aero.stations[index].airfoil.polars[0]
    alpha = np.degrees(polar.alpha)
    return (
        float(np.interp(alpha_deg, alpha, polar.cl)),
        float(np.interp(alpha_deg, alpha, polar.cd)),
        float(np.interp(alpha_deg, alpha, polar.cm)),
    )


def main() -> int:
    if not AD_PRIMARY.exists():
        print(f"AeroDyn deck not present: {AD_PRIMARY}")
        return 2

    aero = build_blade_aero_from_aerodyn(AD_PRIMARY)
    r_hub = np.asarray(aero.r, dtype=float)
    twist_deg = np.degrees([st.twist for st in aero.stations])
    rotor_radius = float(aero.rotor_radius)

    bem = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.2).compute(
        V_RATED, RPM_ZHOU, PITCH_RATED
    )

    # ---------------------------------------------------------------- our side
    # The production output is rotor-plane by construction. Pin it before trusting the
    # same identity on their data: if this ever stops holding, the wrapper changed frame.
    psi_ours = np.degrees(np.arctan2(bem.Tp, bem.Np))
    gamma_ours = np.degrees(np.arctan2(bem.cd, bem.cl))
    identity_rotor = psi_ours + gamma_ours - twist_deg - bem.alpha
    identity_section = psi_ours + gamma_ours - bem.alpha

    print("our production BEM at Zhou's rated point (V=10.59 m/s, 7.55 rpm, pitch 0)")
    print(f"  thrust {bem.thrust / 1e6:.3f} MN   power {bem.power / 1e6:.3f} MW"
          f"   torque {bem.torque / 1e6:.3f} MN.m")
    print(f"  int Np dr = {np.trapezoid(bem.Np, r_hub) / 1e3:.1f} kN"
          f"   int Tp dr = {np.trapezoid(bem.Tp, r_hub) / 1e3:.1f} kN")
    print("  frame identity, residual over the 50 stations [deg]:")
    print(f"    rotor-plane reading (psi + gamma - theta - alpha): max "
          f"{np.max(np.abs(identity_rotor)):.3e}")
    print(f"    section reading     (psi + gamma - alpha)        : max "
          f"{np.max(np.abs(identity_section)):.3e}")
    print("    => ccblade/cl.BEMSolver emit the ROTOR-PLANE pair (it is qc(CL cos(phi), ...))")

    # -------------------------------------------------------------- their side
    loads = pd.read_csv(ZHOU_LOADS)
    aoa = pd.read_csv(ZHOU_AOA)
    r_R = loads["r_R"].to_numpy(dtype=float)
    np_zhou = loads["normal_force_kN_m"].to_numpy(dtype=float) * 1.0e3
    tp_zhou = loads["tangential_force_kN_m"].to_numpy(dtype=float) * 1.0e3
    alpha_zhou = np.interp(
        r_R, aoa["r_R"].to_numpy(dtype=float), aoa["angle_of_attack_deg"].to_numpy(dtype=float)
    )

    print()
    print("their published pair (Fig. 11 Np,Tp) against their published alpha (Fig. 10),")
    print("official AeroDyn polar at that alpha, deck twist for the rotor-plane reading:")
    header = (
        f"{'r/R':>5} {'alpha_10':>9} {'CL':>6} {'CD':>7} {'gamma':>6} {'psi':>6} "
        f"{'a_section':>10} {'d_sec':>7} {'a_rotor':>8} {'d_rot':>7}"
    )
    print(header)
    d_section: list[float] = []
    d_rotor: list[float] = []
    r_R_used: list[float] = []
    for i, rr in enumerate(r_R):
        r_abs = rr * rotor_radius
        if rr < 0.15 or r_abs <= r_hub[0] or r_abs >= r_hub[-1]:
            continue
        j = int(np.argmin(np.abs(r_hub - r_abs)))
        cl, cd, _ = _polar_at(aero, j, alpha_zhou[i])
        gamma = np.degrees(np.arctan2(cd, cl))
        psi = np.degrees(np.arctan2(tp_zhou[i], np_zhou[i]))
        a_section = psi + gamma
        a_rotor = psi + gamma - twist_deg[j]
        d_section.append(a_section - alpha_zhou[i])
        d_rotor.append(a_rotor - alpha_zhou[i])
        r_R_used.append(rr)
        print(f"{rr:5.2f} {alpha_zhou[i]:9.2f} {cl:6.3f} {cd:7.4f} {gamma:6.2f} {psi:6.2f} "
              f"{a_section:10.2f} {a_section - alpha_zhou[i]:+7.2f} "
              f"{a_rotor:8.2f} {a_rotor - alpha_zhou[i]:+7.2f}")

    d_section_a = np.asarray(d_section)
    d_rotor_a = np.asarray(d_rotor)
    r_R_used_a = np.asarray(r_R_used)
    core = (r_R_used_a >= CORE_LO) & (r_R_used_a <= CORE_HI)
    print()
    print(f"  section reading residual: max |d| {np.max(np.abs(d_section_a)):.2f} deg, "
          f"rms {np.sqrt(np.mean(d_section_a**2)):.2f} deg")
    print(f"  rotor   reading residual: max |d| {np.max(np.abs(d_rotor_a)):.2f} deg, "
          f"rms {np.sqrt(np.mean(d_rotor_a**2)):.2f} deg")
    print(f"  blade core r/R [{CORE_LO}, {CORE_HI}] ({int(core.sum())} stations): "
          f"section max {np.max(np.abs(d_section_a[core])):.2f} deg, "
          f"rotor max {np.max(np.abs(d_rotor_a[core])):.2f} deg "
          f"(Fig. 10 band +/-{FIG10_AOA_BAND_DEG} deg)")
    print("  reading: over the blade core the section reading is within a few tenths of a")
    print("  degree of their own published alpha, and the rotor-plane reading is not; the")
    print("  section reading's worst points are the digitized tip, where Fig. 10 falls")
    print("  from 4.25 to 1.22 deg between two digitized stations.")

    # ---------------------------------------------------- quadrature of the curve
    r_abs_z = r_R * rotor_radius
    blade_thrust_mn = np.trapezoid(np_zhou, r_abs_z) * 3.0 / 1e6
    blade_torque_mnm = np.trapezoid(tp_zhou * r_abs_z, r_abs_z) * 3.0 / 1e6
    omega_rad = RPM_ZHOU * 2.0 * np.pi / 60.0
    power_mw = blade_torque_mnm * omega_rad
    print()
    print("the digitized curve against the paper's own Table 6 (fixed condition):")
    print(f"  their curve gives thrust {blade_thrust_mn:.2f} MN  (Table 6: rigid "
          f"{TABLE6_RIGID[1]:.2f}, flexible {TABLE6_FLEXIBLE[1]:.2f})")
    print(f"  their curve gives power  {power_mw:.2f} MW  (Table 6: rigid "
          f"{TABLE6_RIGID[0]:.2f}, flexible {TABLE6_FLEXIBLE[0]:.2f})")
    print("  => Fig. 11 is the flexible case, and its trapezoid does not carry the")
    print("     paper's own torque; integral comparisons through it carry that band.")
    print(f"  the pointwise Tp digitization band alone is +/-{FIG11_TP_BAND_KN_M:.2f} kN/m")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
