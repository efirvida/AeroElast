"""L4 T5 matched-options measurement: sensitivity of rated rigid thrust/torque/power to every
option that differs between our BEM (CCBlade wrapper) and the AeroDyn deck used as arbiter.

Measurement instrument only: no band or anchor is touched. Usage:
    scripts/aeroenv.sh env PYTHONPATH=$PWD OPENFAST_BIN=... python tools/l4_t5_options.py --output logs/l4_t4t5/t5_effects.json
Helpers (`ad_run`, `ours`) are reused by tools/l4_t4_zhou.py.
"""

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.validation.bem import test_bem_openfast_parity as p  # noqa: E402

ob = p.ob
RHO, MU_AD = 1.225, 1.464e-5
MU_DEFAULT = 1.81206e-5  # BEMSolver / production default


def ad_run(case, *, overrides=None, tilt, precone, dt=0.25, tmax=40.0):
    """One AeroDyn driver run (quasi-steady, tower off); returns (blade_aero, AeroDynResult)."""
    driver = p._openfast_driver_or_skip()
    with tempfile.TemporaryDirectory(prefix="l4t5-") as td:
        td = Path(td)
        primary = ob.write_bem_primary(p.PRIMARY, td / "AD_clean.dat", overrides)
        blade = ob.build_blade_aero_from_aerodyn(primary)
        dvr = ob.write_aerodyn_dvr(
            td / "t.dvr",
            primary,
            [case],
            dt=dt,
            tmax=tmax,
            precone_deg=precone,
            shaft_tilt_deg=tilt,
        )
        ref = ob.parse_aerodyn_outs(ob.run_aerodyn_driver(driver, dvr, td))[0]
    return blade, ref


def ours(blade, case, *, tilt, precone, mu=MU_AD, nsector=None, linear_polars=False):
    s = p.BEMSolver(
        blade,
        rho=RHO,
        mu=mu,
        precone=precone,
        tilt=tilt,
        yaw=0.0,
        hub_height=ob.IEA15MW_HUB_HEIGHT,
        shear_exp=0.0,
    )
    if nsector:
        s.rotor.nSector = nsector
    if linear_polars:
        _patch_linear(s, blade)
    r = s.compute(case[0], case[1], case[2])
    om = case[1] * 2 * np.pi / 60
    return dict(thrust=r.thrust, torque=r.torque, power=r.torque * om, ccblade_power=r.power)


def _patch_linear(solver, blade):
    """Replace CCAirfoil's smoothed cubic spline by AeroDyn-style linear interpolation of the table."""
    seen = {}
    for st in blade.stations:
        seen[st.airfoil.name] = st.airfoil.polars[0]
    by_id = {a.name: a for a in blade.airfoils}
    for af, st in zip(solver.rotor.af, blade.stations):
        pol = by_id[st.airfoil.name].polars[0]
        a, cl, cd = pol.alpha, pol.cl, pol.cd

        def ev(alpha, Re, a=a, cl=cl, cd=cd):
            return np.interp(alpha, a, cl), np.interp(alpha, a, cd)

        af.evaluate = ev


def rel(x, ref):
    return {k: (x[k] - ref[k]) / ref[k] for k in ("thrust", "torque", "power")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    tilt, cone = ob.IEA15MW_SHAFT_TILT_DEG, ob.IEA15MW_PRECONE_DEG
    case = (10.59, ob.IEA15MW_RATED_RPM, 0.0)
    om = case[1] * 2 * np.pi / 60
    out = {"case": case, "tilt": tilt, "precone": cone}

    def adpack(ref):
        return dict(thrust=ref.thrust, torque=ref.torque, power=ref.torque * om)

    blade, ref_deck = ad_run(
        case, tilt=tilt, precone=cone
    )  # deck as the tests write it (Skew_Mod=1)
    _, ref_noskew = ad_run(case, overrides={"Skew_Mod": "0"}, tilt=tilt, precone=cone)
    out["AD_deck_skew1"] = adpack(ref_deck)
    out["AD_skew0"] = adpack(ref_noskew)
    out["AD_skew_effect_on_AD"] = rel(adpack(ref_deck), adpack(ref_noskew))
    base = ours(blade, case, tilt=tilt, precone=cone)
    out["ours_base"] = base
    out["ours_vs_AD_skew1"] = rel(base, adpack(ref_deck))
    out["ours_vs_AD_skew0"] = rel(base, adpack(ref_noskew))
    var = {
        "mu_1.81206e-5_(production default)": dict(mu=MU_DEFAULT),
        "nSector_24": dict(nsector=24),
        "linear_polar_interp": dict(linear_polars=True),
    }
    for k, kw in var.items():
        r = ours(blade, case, tilt=tilt, precone=cone, **kw)
        out[f"ours[{k}]"] = r
        out[f"ours[{k}]_vs_ours_base"] = rel(r, base)
        out[f"ours[{k}]_vs_AD_skew1"] = rel(r, adpack(ref_deck))
    r0 = ours(blade, case, tilt=0.0, precone=0.0)
    out["ours[tilt0_precone0]"] = r0
    out["ours[tilt0_precone0]_vs_ours_base"] = rel(r0, base)
    # YAML (production) blade: Rtip inconsistency and polars
    yb = p.load_blade_aero(str(p.YAML_BLADE))
    out["yaml_blade"] = dict(
        Rhub=yb.hub_radius, Rtip=yb.rotor_radius, r_tip_station=float(yb.r[-1])
    )
    ry = ours(yb, case, tilt=tilt, precone=cone)
    out["ours[yaml_blade]"] = ry
    out["ours[yaml_blade]_vs_AD_skew1"] = rel(ry, adpack(ref_deck))
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, default=float) + "\n")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
