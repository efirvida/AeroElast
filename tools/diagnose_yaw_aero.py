"""Isolate the yaw degradation one knob at a time (aero side).

The FSI campaign mixes four things: the aero model, the yaw/wake model, the
aeroelastic coupling and the structural solver.  This tool attacks the aero
half with a ladder in which each rung changes exactly ONE thing, all at the
same operating point on the same rotor:

  R1  CCBlade (ours, in-process)          no skew model exists at all
  R2  AeroDyn driver, Skew_Mod = 0        skew explicitly off
  R3  AeroDyn driver, Skew_Mod = 1        skew on (Glauert/Pitt/Peters)

Reading:

  R2 - R1   the code/config difference with NO skew on either side
  R3 - R2   the skew model alone (same code, same polars, same settings)
  R3 - R1   the total aero-model gap

The UA and DBEMT models are switched OFF in R2/R3 so the comparison is
steady BEM vs steady BEM (CCBlade is steady too); that keeps the skew as the
single deliberate difference.  Turn them back on later as a separate rung.

Everything here is configuration + orchestration: no solver code is touched.

Usage:
  python tools/diagnose_yaw_aero.py --work-dir $SCRATCH/tmp/opencode/yawdiag
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

IEA_OPENFAST = (
    REPO / "tests" / "IEA15MW" / "validation papers"
    / "github-IEA-15-240-RWT" / "OpenFAST"
)
# The vendored IEA deck is the v4.2 schema; the CONVERTED v5.0 deck (with
# NacelleDrag / TwrAero / AA_InputFile) survived the /tmp wipe in $SCRATCH.
OFRUNS_DECK = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST")
DECK_AD = OFRUNS_DECK / "IEA-15-240-RWT-Monopile" / "IEA-15-240-RWT-Monopile_AeroDyn15.dat"
BLADE_SRC = OFRUNS_DECK / "IEA-15-240-RWT" / "IEA-15-240-RWT_AeroDyn15_blade.dat"
AIRFOIL_DIR = OFRUNS_DECK / "IEA-15-240-RWT" / "Airfoils"

OPENFAST_BIN = Path("/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/aerodyn_driver")

WIND = 10.59
RPM = 7.55
PITCH = 0.0
RHO = 1.225
KINVISC = 1.81e-5
R_TIP = 120.97

# fields rewritten in the AeroDyn deck for the like-for-like rungs
AD_EDIT = {
    "TwrPotent": "0",     # no tower influence (CCBlade has none)
    "TwrShadow": "0",
    "TwrAero": "False",   # no tower drag
    "UA_Mod": "0",        # quasi-steady polars only (no unsteady aero)
    "DBEMT_Mod": "0",     # no dynamic inflow
}
OUT_LIST = [
    "RtAeroPwr", "RtAeroFxh", "RtAeroFyh", "RtAeroMxh", "RtAeroMyh", "RtAeroMzh",
    "RtAeroCp", "RtAeroCt", "RtSpeed", "RtTSR", "RtArea",
]


def _set_field(text: str, field: str, value: str) -> str:
    """Replace the leading value on the line whose comment names ``field``."""
    pat = re.compile(rf"^(\s*)(\S+)(\s+{re.escape(field)}\b)", re.MULTILINE)

    def repl(m):
        return f"{m.group(1)}{value}{m.group(3)}"

    new, n = pat.subn(repl, text)
    if n == 0:
        print(f"  [warn] AeroDyn field {field!r} not present in this deck — skipped")
        return text
    if n > 1:
        raise RuntimeError(f"AeroDyn field {field!r} matched {n} times (expected 1)")
    return new


def _set_outlist(text: str) -> str:
    """Replace the first OutList block with our channel list."""
    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines) if "OutList" in ln)
    end = next(i for i in range(start, len(lines)) if lines[i].strip().startswith("END"))
    block = [f'"{c}"' for c in OUT_LIST] + ["END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)"]
    return "\n".join(lines[:start + 1] + block + lines[end + 1:]) + "\n"


def build_case(work: Path, tag: str, overrides: dict, yaws, tmax: float, dt: float) -> Path:
    """Write AeroDyn_<tag>.dat + driver_<tag>.dvr in ``work``. Returns the .dvr."""
    work.mkdir(parents=True, exist_ok=True)
    # inputs next to the case (no spaces in paths)
    shutil.copy(BLADE_SRC, work / "blade.dat")
    if not (work / "Airfoils").exists():
        shutil.copytree(AIRFOIL_DIR, work / "Airfoils")

    text = DECK_AD.read_text()
    edits = dict(AD_EDIT)
    edits.update(overrides)
    for f, v in edits.items():
        text = _set_field(text, f, v)
    # paths -> local copies
    text = re.sub(r'"[^"]*_AeroDyn15_blade\.dat"', '"blade.dat"', text)
    text = re.sub(r'"[^"]*Airfoils/([^"/]+)"', r'"Airfoils/\1"', text)
    text = _set_outlist(text)
    ad_file = work / f"AeroDyn_{tag}.dat"
    ad_file.write_text(text)

    rows = "\n".join(
        f"{WIND:10.4f} {0.0:8.3f} {RPM:9.3f} {PITCH:8.3f} {y:8.3f} {dt:8.4f} {tmax:8.2f} "
        f"{0:4d} {0.0:10.4f} {0.0:10.4f}"
        for y in yaws
    )
    dvr = f"""----- AeroDyn Driver Input File ---------------------------------------------------------
Yaw-model diagnostic: same rotor and conditions as the CCBlade rung
----- Input Configuration ---------------------------------------------------------------
False           Echo                   - Echo input parameters to "<rootname>.ech"?
          0     MHK                    - MHK turbine type
          3     AnalysisType           - {{1: multiple turbines, one sim, 2: one turbine, 3: combined cases}}
         60     TMax                   - Total run time [used only when AnalysisType/=3] (s)
       0.01     DT                     - Simulation time step [used only when AnalysisType/=3] (s)
"{ad_file.name}" AeroFile - Name of the primary AeroDyn input file
----- Environmental Conditions ----------------------------------------------------------
    {RHO:8.4f}     FldDens                - Density of working fluid (kg/m^3)
 {KINVISC:9.2E}     KinVisc                - Kinematic viscosity of working fluid (m^2/s)
      340.0     SpdSound               - Speed of sound in working fluid (m/s)
     101325     Patm                   - Atmospheric pressure (Pa) [MHK only]
       2500     Pvap                   - Vapour pressure of working fluid (Pa) [MHK only]
         50     WtrDpth                - Water depth (m)
----- Inflow Data -----------------------------------------------------------------------
          0     CompInflow             - Compute inflow wind velocities (0=Steady Wind; 1=InflowWind)
"unused"        InflowFile             - InflowWind input file [used only when CompInflow=1]
   {WIND:8.4f}     HWindSpeed             - Horizontal wind speed [used only when CompInflow=0 and AnalysisType=1] (m/s)
        150     RefHt                  - Reference height for horizontal wind speed (m)
        0.0     PLExp                  - Power law exponent (uniform flow) (-)
----- SeaState Data [MHK only] ----------------------------------------------------------
          0     CompSeaSt              - Compute wave velocities (0=No Waves; 1=SeaState)
"unused"        SeaStFile              - SeaState input file [used only when CompSeaSt=1]
----- Turbine Data ----------------------------------------------------------------------
          1     NumTurbines            - Number of turbines
----- Turbine(1) Geometry ---------------------------------------------------------------
True            BasicHAWTFormat(1)     - Flag basic vs generic input format
      0,0,0     BaseOriginInit(1)      - Coordinate of tower base in base coordinates (m)
          3     NumBlades(1)           - Number of blades (-)
        1.5     HubRad(1)              - Hub radius (m)
      150.0     HubHt(1)               - Hub height (m)
      -11.5     Overhang(1)            - Overhang (m)
        0.0     ShftTilt(1)            - Shaft tilt (deg)
        0.0     Precone(1)             - Blade precone (deg)
        0.0     Twr2Shft(1)            - Vertical distance from tower-top to rotor shaft (m)
----- Turbine(1) Motion [used only when AnalysisType=1] ---------------------------------
          0     BaseMotionType(1)      - Base motion type
          1     DegreeOfFreedom(1)     - DOF for sinusoidal motion
          0     Amplitude(1)           - Amplitude
          0     Frequency(1)           - Frequency
"unused"        BaseMotionFileName(1)  - Arbitrary base motion file
          0     NacYaw(1)              - Yaw angle about z_t (deg)
     {RPM:7.3f}     RotSpeed(1)            - Rotational speed of rotor (rpm)
     {PITCH:7.3f}     BldPitch(1)            - Blade 1 pitch (deg)
----- Time-dependent Analysis [used only when AnalysisType=2] ---------------------------
"unused"        TimeAnalysisFileName   - Time series file
----- Combined-Case Analysis [used only when AnalysisType=3] ----------------------------
{len(list(yaws)):11d}     NumCases               - Number of cases to run
HWndSpeed     PLExp     RotSpd     Pitch     Yaw     dT     Tmax     DOF     Amplitude     Frequency
(m/s)         (-)       (rpm)      (deg)     (deg)   (s)    (s)      (-)     (-)           (Hz)
{rows}
----- Output Settings -------------------------------------------------------------------
"ES15.8E2"      OutFmt                 - Format for text output
          1     OutFileFmt             - Output file format (1=text)
          0     WrVTK                  - VTK output (0=none)
          1     WrVTK_Type             - VTK type
          0     VTKHubRad              - Hub radius for VTK
0,0,0,0,0,0     VTKNacDim              - Nacelle dimension for VTK
EOF
"""
    dvr_path = work / f"driver_{tag}.dvr"
    dvr_path.write_text(dvr)
    return dvr_path


def run_driver(work: Path, dvr: Path) -> Path:
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "4")
    proc = subprocess.run(
        [str(OPENFAST_BIN), dvr.name], cwd=work, env=env,
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"aerodyn_driver failed for {dvr.name} (rc={proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-3000:]}\n--- stderr ---\n{proc.stderr[-2000:]}"
        )
    cands = sorted(work.glob(dvr.stem + ".*.out"), key=lambda p: int(p.stem.split(".")[-1]))
    if not cands:
        raise RuntimeError(f"no per-case .out produced for {dvr.name}: {proc.stdout[-1500:]}")
    return cands


def read_out(path: Path) -> tuple[list[str], np.ndarray]:
    """Read an AeroDyn driver text output (banner + header + units + rows)."""
    lines = path.read_text().splitlines()
    hdr_i = next(i for i, ln in enumerate(lines) if ln.strip().startswith("Time"))
    names = lines[hdr_i].split()
    rows = []
    for ln in lines[hdr_i + 1:]:
        toks = ln.split()
        if len(toks) < 2:
            continue
        try:
            rows.append([float(t) for t in toks])
        except ValueError:
            continue  # units line / banner
    return names, np.array(rows)


def aero_stats(names, data, tail_frac: float = 0.5):
    """Tail-average the requested rotor channels of one case."""
    out = {}
    for key in ("RtAeroPwr", "RtAeroFxh", "RtAeroMxh", "RtAeroCp", "RtAeroCt"):
        if key not in names:
            out[key] = float("nan")
            continue
        block = data[:, names.index(key)]
        out[key] = float(block[int(len(block) * (1 - tail_frac)):].mean())
    return out


def ccblade_rung(yaws):
    """R1: our in-process CCBlade, identical settings to the FSI participant."""
    from aeroelast.models.blade.aerodynamics import load_blade_aero
    from aeroelast.solvers.bem.engine import BEMSolver

    blade = load_blade_aero(
        str(REPO / "tests" / "IEA-15-240-RWT.yaml"),
        default_re=1.0e7, neuralfoil_model="large", hub_radius=0.0, n_blades=3,
        viterna_ar=17.0, viterna_confidence_threshold=0.5,
    )
    res = {}
    for y in yaws:
        s = BEMSolver(blade, rho=RHO, mu=KINVISC, precone=0.0, tilt=0.0, yaw=float(y),
                      hub_height=150.0, shear_exp=0.0)
        r = s.compute(WIND, RPM, PITCH, azimuth=0.0)
        res[y] = dict(P=r.power / 1e6, T=r.thrust / 1e6, Q=r.torque / 1e6,
                      CP=r.power / (0.5 * RHO * math.pi * R_TIP**2 * WIND**3),
                      CT=r.thrust / (0.5 * RHO * math.pi * R_TIP**2 * WIND**2))
    return res


VARIANTS = (
    ("skew0", "R2", {"Skew_Mod": "0"}),
    ("skew1", "R3", {"Skew_Mod": "1"}),
    ("deck", "R4", {"Skew_Mod": "1", "DBEMT_Mod": "2", "UA_Mod": "3"}),
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--yaws", type=float, nargs="+", default=[0, 5, 10, 15, 20, 25, 30, 35, 40])
    ap.add_argument("--tmax", type=float, default=40.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--skip-ccblade", action="store_true")
    ap.add_argument("--variants", nargs="*", default=[v[0] for v in VARIANTS],
                    help="AeroDyn variants to run (default: all)")
    ap.add_argument("--csv", type=Path, default=None)
    args = ap.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    yaws = [float(y) for y in args.yaws]

    print("=== R1: CCBlade (ours, no skew model exists) ===")
    r1 = {} if args.skip_ccblade else ccblade_rung(yaws)
    if r1:
        for y in yaws:
            print(f"  yaw {y:4.0f}:  P={r1[y]['P']:7.3f} MW  T={r1[y]['T']:6.3f} MN  "
                  f"Q={r1[y]['Q']:7.3f} MNm  CP={r1[y]['CP']:.4f}")

    runners: dict[str, dict] = {}
    order: list[str] = []
    for tag, label, overrides in VARIANTS:
        if tag not in args.variants:
            continue
        print(f"\n=== {label}: AeroDyn driver [{tag}] settings={overrides} ===")
        dvr = build_case(args.work_dir, tag, overrides, yaws, args.tmax, args.dt)
        outs = run_driver(args.work_dir, dvr)
        if len(outs) != len(yaws):
            raise RuntimeError(f"{len(outs)} case files for {len(yaws)} yaws")
        runners[label] = {}
        for y, out in zip(yaws, outs):
            names, data = read_out(out)
            runners[label][y] = aero_stats(names, data)
        order.append(label)
        print(f"  outputs: {len(outs)} case files, {len(names)} channels, {len(data)} rows each")
        for y in yaws:
            s = runners[label][y]
            print(f"  yaw {y:4.0f}:  P={s['RtAeroPwr']/1e6:7.3f} MW  "
                  f"T={s['RtAeroFxh']/1e6:6.3f} MN  Q={s['RtAeroMxh']/1e6:7.3f} MNm  "
                  f"CP={s['RtAeroCp']:.4f}")

    print("\n" + "=" * 104)
    print("DIAGNOSTIC TABLE — P [MW]")
    print("=" * 104)
    head = f"{'yaw':>4} " + " ".join(f"{lab:>11}" for lab in ["R1"] + order)
    head += " " + " ".join(f"{o+'-prev%':>9}" for o in order)
    print(head)
    p0 = {lab: runners[lab][0.0]["RtAeroPwr"] / 1e6 for lab in order}
    r1_p0 = r1[0.0]["P"] if r1 else float("nan")
    rows = []
    for y in yaws:
        vals = [r1[y]["P"] if r1 else float("nan")]
        vals += [runners[lab][y]["RtAeroPwr"] / 1e6 for lab in order]
        line = f"{y:4.0f} " + " ".join(f"{v:11.3f}" for v in vals)
        line += " " + " ".join(f"{100*(vals[i+1]-vals[i])/vals[i]:9.2f}" for i in range(len(order)))
        print(line)
        rows.append([y] + vals)

    base_of = {"R1": (r1_p0, lambda y: r1[y]["P"] if r1 else float("nan"))}
    for lab in order:
        base_of[lab] = (p0[lab], lambda y, lab=lab: runners[lab][y]["RtAeroPwr"] / 1e6)
    for lab in ["R1"] + order:
        base, get = base_of[lab]
        if not (base == base) or base <= 0:
            continue
        ns = [math.log(get(y) / base) / math.log(math.cos(math.radians(y)))
              for y in (10, 20, 30, 40) if y in yaws]
        print(f"  cos^n [{lab}]: " + "  ".join(f"yaw {y:2.0f}->{n:5.2f}" for y, n in zip((10,20,30,40), ns)))

    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["yaw_deg"] + ["R1_ccblade_MW"] + [f"{lab}_MW" for lab in order])
            w.writerows(rows)
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
