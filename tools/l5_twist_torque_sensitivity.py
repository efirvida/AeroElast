"""Rated torque/power/thrust sensitivity to a spanwise elastic-twist difference (L5 T8, S2 of odd/tasks/l5-one-way.md).

A T1/T3 twist red (shell vs BeamDyn elastic twist, ``tools/l5_oneway_bd_parity.py``) is declared with its measured effect on the
thesis quantity: the production BEM (``BEMSolver``, T7 defaults: prebend on, exact polar interpolation) is run on the AeroDyn-deck
blade at the rated point (10.59 m/s, 7.56 rpm, pitch 0, tilt -6, precone -4, all explicit flags) with three twist distributions

* rigid: the deck twist;
* BeamDyn: deck twist + BeamDyn elastic twist (``table["beamdyn_x0.01"]["twist_deg"]``, the like-for-like ``lin`` reference);
* shell: deck twist + shell elastic twist (``table["shell"]["twist_<estimator>_deg"]``, ``--estimator skin|eq``, default skin).

The elastic twist is linearly interpolated onto the BEM stations (span z = BlSpn = r - hub radius), exactly at the table nodes, zero
inboard of the first z and held at the last value outboard.  The output is each case plus delta(shell - BeamDyn) in absolute terms
and in % of the rigid value.  One-way, quasi-static: a sensitivity of the rotor integrals to twist, not a coupled answer
(same spirit and constants as ``tools/l3_torque_sensitivity.py``, whose rated operating point is reused).

Sign applied.  The table twists are the mesh/BeamDyn elastic twist: right-hand rotation about +z, nose-up positive (R1 of
``tools/twist_p1.py``; ``tools/l5_oneway_bd_parity.py``, docstring).  The production feedback converts it to the BEM twist with
``theta_BEM = theta_ref + _twist_mesh_to_bem * omega`` (``src/aeroelast/solvers/bem/fsi_participant.py:824``) and the factor is -1
for the deck frame (leading edge at +x, downwind +y; derivation ``fsi_participant.py:494-503`` and ``_resolve_twist_sense``,
``:863``): the mesh rotation and CCBlade's theta are opposite-sensed (``alpha = phi - theta``), so a positive (nose-up) elastic twist
LOWERS theta and raises the angle of attack toward stall, a negative one raises theta toward feather.  This tool applies that
production factor (``MESH_TO_BEM_SENSE = -1``) to both the shell and the BeamDyn twist, never flipped to make a delta look better.
Wording caveat: ``odd/tasks/rated-twist-sign-convention.md`` (E2, "nose-down is omega > 0") and
``tests/validation/blade/test_blade_deloading_vs_reference.py:222-228`` ("nose-down is omega < 0") word the geometric label
differently from the code comment; all of them agree on the operative rule used here, ``d theta_BEM = -omega_mesh``.

Run (BEM only, seconds):
    scripts/aeroenv.sh python tools/l5_twist_torque_sensitivity.py --json logs/check/l5/bd_parity_bem.json --out sens.json
    scripts/aeroenv.sh python tools/l5_twist_torque_sensitivity.py --z 60 90 117 --dtwist-deg 0 -0.5 -1.0 --out sens.json
``--z/--dtwist-deg`` is the shell twist; BeamDyn is then the rigid (zero) twist unless ``--bd-dtwist-deg`` is given.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tools"))

from l3_torque_sensitivity import PITCH_RATED, RPM_RATED, V_RATED  # noqa: E402  (rated operating point, reused)

#: d theta_BEM = MESH_TO_BEM_SENSE * omega_mesh (fsi_participant.py:824, factor from :863; -1 on the deck frame).
MESH_TO_BEM_SENSE = -1.0
QUANTITIES = ("thrust", "torque", "power")


def interpolate_dtwist(z_query, z_nodes, dtwist) -> np.ndarray:
    """Linear interpolation of ``dtwist`` given at ``z_nodes`` onto ``z_query``: exact at nodes, 0 inboard of the first node,
    last value held outboard."""
    z_nodes = np.asarray(z_nodes, dtype=float)
    dtwist = np.asarray(dtwist, dtype=float)
    if z_nodes.shape != dtwist.shape or z_nodes.ndim != 1 or z_nodes.size == 0:
        raise ValueError("z and dtwist must be 1-D lists of equal, non-zero length")
    if np.any(np.diff(z_nodes) <= 0.0):
        raise ValueError("z must be strictly increasing")
    return np.interp(np.asarray(z_query, dtype=float), z_nodes, dtwist, left=0.0)


def mesh_twist_to_bem(omega_mesh) -> np.ndarray:
    """Mesh/BeamDyn elastic twist (nose-up positive about +z) -> BEM theta increment, same unit (fsi_participant.py:824)."""
    return MESH_TO_BEM_SENSE * np.asarray(omega_mesh, dtype=float)


def bem_twist(ref_twist_rad, z_stations, z_nodes, dtwist_deg_mesh) -> np.ndarray:
    """BEM twist [rad] per station: deck twist + converted, interpolated elastic twist (mesh sense, deg)."""
    om = np.deg2rad(interpolate_dtwist(z_stations, z_nodes, dtwist_deg_mesh))
    return np.asarray(ref_twist_rad, dtype=float) + mesh_twist_to_bem(om)


def read_parity_json(path, estimator: str):
    """``(z, shell twist [deg], BeamDyn twist [deg])`` from a ``tools/l5_oneway_bd_parity.py`` JSON."""
    if estimator not in ("skin", "eq"):
        raise ValueError(f"estimator must be 'skin' or 'eq', got {estimator!r}")
    t = json.loads(Path(path).read_text())["table"]
    return (np.asarray(t["z"], float), np.asarray(t["shell"][f"twist_{estimator}_deg"], float),
            np.asarray(t["beamdyn_x0.01"]["twist_deg"], float))


def delta_shell_minus_beamdyn(cases: dict) -> dict:
    """delta(shell - BeamDyn) per quantity: absolute and in % of the rigid value."""
    return {q: {"abs": cases["shell"][q] - cases["beamdyn"][q],
                "pct_of_rigid": 100.0 * (cases["shell"][q] - cases["beamdyn"][q]) / cases["rigid"][q]} for q in QUANTITIES}


def build_blade_aero():
    from tests.support.openfast_bem import DEFAULT_DECK, build_blade_aero_from_aerodyn

    return build_blade_aero_from_aerodyn(DEFAULT_DECK / "case" / "IEA-15-240-RWT_AeroDyn15.dat")


def evaluate(blade_aero, theta_rad, tilt, precone, mu, wind=V_RATED, rpm=RPM_RATED, pitch=PITCH_RATED) -> dict[str, float]:
    """Production BEM (prebend on, exact polar interpolation = engine defaults) with the station twist replaced by ``theta_rad``."""
    import dataclasses

    from aeroelast.solvers.bem.engine import BEMSolver

    stations = [dataclasses.replace(s, twist=float(t)) for s, t in zip(blade_aero.stations, theta_rad, strict=True)]
    solver = BEMSolver(dataclasses.replace(blade_aero, stations=stations), rho=1.225, mu=mu, precone=precone, tilt=tilt, yaw=0.0,
                       hub_height=150.0, shear_exp=0.0)
    res = solver.compute(wind, rpm, pitch)
    return {"thrust": float(res.thrust), "torque": float(res.torque), "power": float(res.power)}


def run(z, shell_deg, bd_deg, tilt=-6.0, precone=-4.0, mu=1.81206e-5, blade_aero=None) -> dict:
    ba = blade_aero or build_blade_aero()
    z_st = ba.r - ba.hub_radius
    ref = ba.twist
    cases = {"rigid": evaluate(ba, ref, tilt, precone, mu),
             "beamdyn": evaluate(ba, bem_twist(ref, z_st, z, bd_deg), tilt, precone, mu),
             "shell": evaluate(ba, bem_twist(ref, z_st, z, shell_deg), tilt, precone, mu)}
    return {"cases": cases, "delta_shell_minus_beamdyn": delta_shell_minus_beamdyn(cases),
            "bem_twist_increment_deg": {"z_station": z_st.tolist(),
                                        "beamdyn": np.rad2deg(bem_twist(ref, z_st, z, bd_deg) - ref).tolist(),
                                        "shell": np.rad2deg(bem_twist(ref, z_st, z, shell_deg) - ref).tolist()}}


def print_table(rep: dict) -> None:
    c, d = rep["cases"], rep["delta_shell_minus_beamdyn"]
    print(f"{'case':10s} {'thrust [N]':>14s} {'torque [N m]':>14s} {'power [W]':>14s}")
    for k in ("rigid", "beamdyn", "shell"):
        print(f"{k:10s} " + " ".join(f"{c[k][q]:14.1f}" for q in QUANTITIES))
    print(f"{'d(sh-BD)':10s} " + " ".join(f"{d[q]['abs']:+14.1f}" for q in QUANTITIES))
    print(f"{'% rigid':10s} " + " ".join(f"{d[q]['pct_of_rigid']:+13.4f}%" for q in QUANTITIES))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--json", type=Path, help="output of tools/l5_oneway_bd_parity.py")
    ap.add_argument("--estimator", choices=("skin", "eq"), default="skin", help="shell twist estimator of the JSON table")
    ap.add_argument("--z", type=float, nargs="+", help="explicit span stations [m] (with --dtwist-deg)")
    ap.add_argument("--dtwist-deg", type=float, nargs="+", help="explicit shell elastic twist [deg], mesh sense (nose-up +)")
    ap.add_argument("--bd-dtwist-deg", type=float, nargs="+", help="explicit BeamDyn elastic twist [deg] at --z (default 0)")
    ap.add_argument("--tilt", type=float, default=-6.0, help="shaft tilt [deg] (deck)")
    ap.add_argument("--precone", type=float, default=-4.0, help="precone [deg] (deck)")
    ap.add_argument("--mu", type=float, default=1.81206e-5, help="dynamic viscosity, as bem_config of l5_oneway_bd_parity")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    if (a.json is None) == (a.z is None):
        ap.error("give exactly one of --json or --z/--dtwist-deg")
    if a.json is not None:
        z, shell, bd = read_parity_json(a.json, a.estimator)
    else:
        if a.dtwist_deg is None:
            ap.error("--z needs --dtwist-deg")
        z, shell = np.asarray(a.z, float), np.asarray(a.dtwist_deg, float)
        bd = np.asarray(a.bd_dtwist_deg, float) if a.bd_dtwist_deg else np.zeros_like(z)
    rep = run(z, shell, bd, a.tilt, a.precone, a.mu)
    rep["inputs"] = {"wind": V_RATED, "rpm": RPM_RATED, "pitch": PITCH_RATED, "tilt": a.tilt, "precone": a.precone, "mu": a.mu,
                     "estimator": a.estimator if a.json else None, "json": str(a.json) if a.json else None,
                     "mesh_to_bem_sense": MESH_TO_BEM_SENSE, "z": z.tolist(), "shell_dtwist_deg": shell.tolist(),
                     "beamdyn_dtwist_deg": bd.tolist()}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rep, indent=1) + "\n")
    print_table(rep)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
