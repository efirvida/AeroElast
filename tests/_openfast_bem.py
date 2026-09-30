"""OpenFAST/AeroDyn helpers for BEM code-to-code validation.

Builds an AeroElast :class:`~aeroelast.models.blade.aerodynamics.BladeAero` and a
comparable OpenFAST ``aerodyn_driver`` case from the **same** AeroDyn deck, so
that the BEM comparison uses identical geometry and identical airfoil polar
tables.  This is the A1 work unit of the wind-turbine validation campaign.

Why both sides are built from the AeroDyn deck
----------------------------------------------
The AeroElast WindIO loader computes the station radius as
``r = hub_radius + eta * blade_length`` while setting
``rotor_radius = rotor_diameter / 2``.  For the IEA 15 MW YAML those two
disagree: ``hub_radius + blade_length = 3.97 + 117.0 = 120.97 m`` but
``rotor_diameter / 2 = 121.1189 m``.  The 0.149 m gap shifts the tip and the
tip-loss region.  This module therefore takes the geometry from the OpenFAST
deck (``HubRad = 3.97``, ``TipRad = 120.97``) and makes ``r = HubRad + BlSpn``
consistent with ``Rtip`` on both codes.

Only the pieces needed for the BEM comparison are parsed: the AeroDyn primary
file (blade file and airfoil list), the blade file, and the AirfoilInfo polar
tables.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from aeroelast.models.blade.aerodynamics import (
    AeroStation,
    AirfoilAero,
    BladeAero,
    PolarData,
)

#: Vendored IEA 15 MW OpenFAST deck, shipped as a single archive so the
#: repository carries one file instead of ~100 (Apache-2.0; see the NOTICE).
VENDORED_DECK_ARCHIVE = Path(__file__).resolve().parent / "reference/iea15mw_openfast.tar.gz"
#: Locally fetched copy of the same deck (`.sources` is gitignored).
FETCHED_DECK = Path(".sources/openfast/iea15mw")
#: Cache root and the deck directory the archive unpacks to.
_CACHE_ROOT = Path(tempfile.gettempdir()) / "aeroelast-openfast-deck"
_CACHED_DECK = _CACHE_ROOT / "iea15mw_openfast"


def _ensure_vendored_deck() -> Optional[Path]:
    """Return the vendored deck directory, unpacking the archive once if needed.

    The archive is unpacked under the temporary directory, never inside the
    repository, so the working tree stays clean.  Returns ``None`` when the
    archive is absent so the caller can fall back to the fetched deck.
    """
    if _CACHED_DECK.joinpath("NOTICE").is_file():
        return _CACHED_DECK
    if not VENDORED_DECK_ARCHIVE.is_file():
        return None
    _CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    with tarfile.open(VENDORED_DECK_ARCHIVE) as archive:
        archive.extractall(_CACHE_ROOT, filter="data")
    return _CACHED_DECK


#: Deck used by the tests: the vendored one when present, else the fetched one.
DEFAULT_DECK = _ensure_vendored_deck() or FETCHED_DECK

#: Official IEA 15 MW rotor geometry, read from the OpenFAST ElastoDyn file
#: (`HubRad`, `TipRad`, `PreCone`, `ShftTilt`, `OverHang`, `Twr2Shft`,
#: `HubHt`).  Kept here so the comparison does not silently inherit the YAML
#: rotor-diameter inconsistency described in the module docstring.
IEA15MW_HUB_RAD = 3.97
IEA15MW_TIP_RAD = 120.97
IEA15MW_PRECONE_DEG = -4.0
IEA15MW_SHAFT_TILT_DEG = -6.0
IEA15MW_OVERHANG = -12.097571763912535
IEA15MW_TWR2SHFT = 4.349459414248071
IEA15MW_HUB_HEIGHT = 150.0
IEA15MW_RATED_RPM = 7.56


def _strip_comment(line: str) -> str:
    return line.split("!", 1)[0].strip()


def _value_lines(text: str) -> List[str]:
    """Return the non-empty, comment-stripped lines of an AeroDyn file."""
    return [v for v in (_strip_comment(ln) for ln in text.splitlines()) if v]


def parse_aerodyn_airfoil(path: str | Path) -> AirfoilAero:
    """Parse an AirfoilInfo (AeroDyn 15) file into an :class:`AirfoilAero`.

    Handles ``NumTabs`` tables, each with ``Re`` (millions), ``Ctrl`` and
    ``InclUAdata``; unsteady-aerodynamics tables (``InclUAdata = True``) are
    rejected because the BEM comparison runs the quasi-steady model.
    """
    lines = _value_lines(Path(path).read_text(errors="replace"))

    def label(line: str) -> str:
        toks = line.split()
        return toks[1] if len(toks) > 1 else ""

    num_tabs = next(int(float(ln.split()[0])) for ln in lines if label(ln) == "NumTabs")

    cursor = 0
    polars: List[PolarData] = []
    for _ in range(num_tabs):
        while cursor < len(lines) and label(lines[cursor]) != "Re":
            cursor += 1
        re_millions = float(lines[cursor].split()[0])
        cursor += 1

        # Skip any unsteady-aerodynamics coefficients (`InclUAdata=True`): the
        # BEM comparison runs the quasi-steady tables, so only the Alpha/Cl/Cd/Cm
        # table that follows `NumAlf` is used.
        while cursor < len(lines) and label(lines[cursor]) != "NumAlf":
            cursor += 1
        n_alf = int(float(lines[cursor].split()[0]))
        cursor += 1

        rows: List[List[float]] = []
        while len(rows) < n_alf:
            toks = lines[cursor].split()
            cursor += 1
            if len(toks) >= 4:
                try:
                    rows.append([float(t) for t in toks[:4]])
                except ValueError:
                    continue
        data = np.array(rows)
        polars.append(
            PolarData(
                alpha=np.deg2rad(data[:, 0]),
                cl=data[:, 1],
                cd=data[:, 2],
                cm=data[:, 3],
                re=re_millions * 1.0e6,
            )
        )

    return AirfoilAero(
        name=Path(path).stem,
        coordinates=np.empty((0, 2), dtype=float),
        relative_thickness=0.0,
        aerodynamic_center=0.25,
        polars=polars,
    )


@dataclass
class AeroDynBlade:
    """Distributed properties parsed from an AeroDyn 15 blade file."""

    bl_spn: np.ndarray  # span from the blade root [m]
    chord: np.ndarray  # [m]
    twist_deg: np.ndarray  # [deg]
    af_id: np.ndarray  # 1-based airfoil table id


def parse_aerodyn_blade(path: str | Path) -> AeroDynBlade:
    """Parse an AeroDyn 15 blade definition file."""
    lines = Path(path).read_text(errors="replace").splitlines()

    header = None
    for idx, ln in enumerate(lines):
        if "BlSpn" in ln and "BlAFID" in ln:
            header = idx
            break
    if header is None:  # pragma: no cover - malformed deck
        raise ValueError(f"{path}: no BlSpn/BlAFID header row found")

    # The blade-node count is on the first non-comment line above the header.
    n_nodes = None
    for ln in lines[:header]:
        v = _strip_comment(ln)
        if v and "NumBlNds" in ln:
            n_nodes = int(float(v.split()[0]))
            break

    rows = []
    for ln in lines[header + 1 :]:
        v = _strip_comment(ln)
        if not v:
            continue
        if v.startswith("("):  # units row
            continue
        fields = v.split()
        if len(fields) < 7:
            continue
        rows.append([float(x) for x in fields[:7]])
        if n_nodes is not None and len(rows) >= n_nodes:
            break

    data = np.array(rows)
    return AeroDynBlade(
        bl_spn=data[:, 0],
        chord=data[:, 5],
        twist_deg=data[:, 4],
        af_id=data[:, 6].astype(int),
    )


def parse_aerodyn_primary(path: str | Path) -> Tuple[Path, List[Path]]:
    """Return ``(blade_file, airfoil_files)`` from an AeroDyn primary file.

    Paths are resolved relative to the primary file, as AeroDyn does.
    """
    primary = Path(path)
    root = primary.parent
    blade: Optional[Path] = None
    airfoils: List[Path] = []

    in_af_block = False
    for ln in primary.read_text(errors="replace").splitlines():
        raw = _strip_comment(ln)
        if not raw:
            continue
        if "ADBlFile(1)" in raw and blade is None:
            blade = root / raw.split()[0].strip('"')
        if "AFNames" in raw:
            in_af_block = True
            # The first filename is on the same line as the `AFNames` label.
            if raw.lstrip().startswith('"'):
                airfoils.append(root / raw.split('"')[1])
            continue
        if in_af_block:
            if raw.startswith('"'):
                airfoils.append(root / raw.strip().strip('"'))
            else:
                in_af_block = False

    if blade is None:  # pragma: no cover - malformed deck
        raise ValueError(f"{path}: no ADBlFile(1) found")
    if not airfoils:  # pragma: no cover - malformed deck
        raise ValueError(f"{path}: no AFNames block found")
    return blade, airfoils


def write_bem_primary(
    src: str | Path,
    dst: str | Path,
    overrides: Optional[Dict[str, str]] = None,
) -> Path:
    """Write a BEM-clean copy of an AeroDyn primary file.

    The OpenFAST IEA 15 MW primary file ships with unsteady aerodynamics
    (``UA_Mod=3``, Beddoes-Leishman) and tower influence/shadow enabled.  The
    AeroElast BEM is quasi-steady and has no tower, so the comparison must turn
    those off on the OpenFAST side too.  This helper copies the primary file,
    rewrites its relative ``"../`` paths to absolute ones (so it can live in a
    temporary directory), and overrides the listed switches.

    The default overrides are the ones that make the OpenFAST BEM match a
    quasi-steady, tower-free BEM:

    ==============  =======  =============================================
    switch          value    why
    ==============  =======  =============================================
    ``UA_Mod``      ``0``    quasi-steady (AeroElast has no unsteady model)
    ``TwrPotent``   ``0``    no tower potential-flow correction
    ``TwrShadow``   ``0``    no tower shadow
    ``TwrAero``     ``False``  no tower aerodynamics
    ``Wake_Mod``    ``1``    BEMT (not OLAF/dynamic)
    ``BEM_Mod``     ``1``    legacy BEM formulation
    ==============  =======  =============================================
    """
    defaults = {
        "UA_Mod": "0",
        "TwrPotent": "0",
        "TwrShadow": "0",
        "TwrAero": "False",
        "Wake_Mod": "1",
        "BEM_Mod": "1",
    }
    if overrides:
        defaults.update(overrides)

    src = Path(src).resolve()
    text = src.read_text(errors="replace")
    base = src.parent.parent
    text = text.replace('"../', f'"{base}/')
    for label, value in defaults.items():
        text = re.sub(
            rf"^(\s*)\S+(\s+{label}\b)", rf"\g<1>{value}\g<2>", text, flags=re.M
        )
    dst = Path(dst)
    dst.write_text(text)
    return dst


def build_blade_aero_from_aerodyn(
    primary_file: str | Path,
    hub_rad: float = IEA15MW_HUB_RAD,
) -> BladeAero:
    """Build the AeroElast :class:`BladeAero` from an AeroDyn primary file.

    The station radius is ``r = hub_rad + BlSpn`` and ``rotor_radius`` is the
    tip radius implied by the same sum, so ``Rhub``/``Rtip`` and the station
    grid are internally consistent on both codes.
    """
    blade_file, airfoil_files = parse_aerodyn_primary(primary_file)
    blade = parse_aerodyn_blade(blade_file)
    airfoils = [parse_aerodyn_airfoil(p) for p in airfoil_files]

    stations: List[AeroStation] = []
    for i in range(len(blade.bl_spn)):
        af = airfoils[int(blade.af_id[i]) - 1]
        stations.append(
            AeroStation(
                span_fraction=float(blade.bl_spn[i] / blade.bl_spn[-1]),
                r=float(hub_rad + blade.bl_spn[i]),
                chord=float(blade.chord[i]),
                twist=float(np.deg2rad(blade.twist_deg[i])),
                pitch_axis=0.25,
                airfoil=af,
            )
        )

    return BladeAero(
        airfoils=airfoils,
        stations=stations,
        blade_length=float(blade.bl_spn[-1]),
        hub_radius=float(hub_rad),
        rotor_radius=float(hub_rad + blade.bl_spn[-1]),
        n_blades=3,
    )


def find_openfast_bin() -> Optional[Path]:
    """Resolve the ``openfast`` binary following ``OPENFAST_BIN`` → PATH → env."""
    env = os.environ.get("OPENFAST_BIN")
    if env and Path(env).exists():
        return Path(env)
    found = shutil.which("openfast")
    if found:
        return Path(found)
    candidate = Path.home() / "miniconda3/envs/openfast/bin/openfast"
    if candidate.exists():
        return candidate
    return None


def aerodyn_driver_of(openfast_bin: Path) -> Optional[Path]:
    """Return the ``aerodyn_driver`` next to ``openfast``, if present."""
    driver = openfast_bin.parent / "aerodyn_driver"
    return driver if driver.exists() else None


def write_aerodyn_dvr(
    path: str | Path,
    primary_file: str | Path,
    cases: Sequence[Tuple[float, ...]],
    *,
    hub_rad: float = IEA15MW_HUB_RAD,
    precone_deg: float = IEA15MW_PRECONE_DEG,
    shaft_tilt_deg: float = IEA15MW_SHAFT_TILT_DEG,
    overhang: float = IEA15MW_OVERHANG,
    twr2shft: float = IEA15MW_TWR2SHFT,
    hub_height: float = IEA15MW_HUB_HEIGHT,
    dt: float = 1.0,
    tmax: float = 60.0,
    shear_exp: float = 0.0,
    yaw_deg: float = 0.0,
) -> Path:
    """Write an AeroDyn driver input file (OpenFAST 4.x format).

    ``cases`` is a sequence of ``(wind_speed, rpm, pitch)`` or
    ``(wind_speed, rpm, pitch, shear_exp, yaw_deg)``; the 3-tuple form uses the
    ``shear_exp`` and ``yaw_deg`` keyword defaults.  The combined-case mode
    (``AnalysisType=3``) is used so one run covers the whole sweep.
    """

    def _case_row(case: Tuple[float, ...]) -> str:
        if len(case) == 3:
            v, rpm, pitch = case
            shear, yaw = shear_exp, yaw_deg
        else:
            v, rpm, pitch, shear, yaw = case
        return (
            f"{v:8.3f}   {shear:<6.3f} {rpm:8.3f} {pitch:8.3f} {yaw:6.2f} "
            f"{dt:6.2f} {tmax:6.1f} 0    0          0"
        )

    rows = "\n".join(_case_row(c) for c in cases)
    aero_file = Path(primary_file).resolve()
    text = f"""----- AeroDyn Driver Input File ---------------------------------------------------------
AeroElast vs AeroDyn BEM parity, cases generated programmatically
----- Input Configuration ---------------------------------------------------------------
False           Echo         - Echo input parameters to "<rootname>.ech"?
        0       MHK          - MHK turbine type (switch) {{0: not an MHK turbine, 1: fixed MHK turbine, 2: floating MHK turbine}}
        3       AnalysisType - {{1: multiple turbines, one simulation, 2: one turbine, one time-dependent simulation, 3: one turbine, combined cases}}
       10.0     TMax         - Total run time [used only when AnalysisType/=3] (s)
        0.5     DT           - Simulation time step [used only when AnalysisType/=3] (s)
"{aero_file}"      AeroFile - Name of the primary AeroDyn input file
----- Environmental Conditions ----------------------------------------------------------
       1.225    FldDens      - Density of working fluid (kg/m^3)
 1.464E-05      KinVisc      - Kinematic viscosity of working fluid (m^2/s)
      340.0     SpdSound     - Speed of sound in working fluid (m/s)
     101325     Patm         - Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
       2000     Pvap         - Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
        150     WtrDpth      - Water depth (m)
----- Inflow Data -----------------------------------------------------------------------
          0      CompInflow  - Compute inflow wind velocities (switch) {{0=Steady Wind; 1=InflowWind}}
"unused"         InflowFile  - Name of the InflowWind input file [used only when CompInflow=1]
        9.0      HWindSpeed  - Horizontal wind speed   [used only when CompInflow=0 and AnalysisType=1] (m/s)
      {hub_height:8.1f}   RefHt       - Reference height for horizontal wind speed [used only when CompInflow=0]  (m)
       {shear_exp:6.3f}   PLExp       - Power law exponent   [used only when CompInflow=0 and AnalysisType=1] (-)
----- Turbine Data ----------------------------------------------------------------------
1               NumTurbines  - Number of turbines
----- Turbine(1) Geometry ---------------------------------------------------------------
        True    BasicHAWTFormat(1) - Flag to switch between basic or generic input format {{True: next 7 lines are basic inputs, False: Base/Twr/Nac/Hub/Bld geometry and motion must follow}}
       0,0,0    BaseOriginInit(1) - Coordinate of tower base in base coordinates (m)
           3    NumBlades(1)    - Number of blades (-)
      {hub_rad:8.4f}    HubRad(1)       - Hub radius (m)
        {hub_height:6.1f}    HubHt(1)        - Hub height (m)
      {overhang:8.4f}    Overhang(1)     - Overhang (m)
        {shaft_tilt_deg:6.2f}    ShftTilt(1)     - Shaft tilt (deg)
        {precone_deg:6.2f}    Precone(1)      - Blade precone (deg)
      {twr2shft:8.4f}    Twr2Shft(1)     - Vertical distance from the tower-top to the rotor shaft (m)
----- Turbine(1) Motion [used only when AnalysisType=1] ---------------------------------
0               BaseMotionType(1)      - Type of motion prescribed for this base {{0: fixed, 1: Sinusoidal motion, 2: arbitrary motion}} (flag)
1               DegreeOfFreedom(1)     - {{1:xt, 2:yt, 3:zt, 4:theta_xt, 5:theta_yt, 6:theta_zt}} [used only when BaseMotionType=1] (flag)
5.0             Amplitude(1)           - Amplitude of sinusoidal motion   [used only when BaseMotionType=1] (m or rad)
0.1             Frequency(1)           - Frequency of sinusoidal motion   [used only when BaseMotionType=1] (Hz)
""              BaseMotionFileName(1)  - Filename containing arbitrary base motion (19 columns)  [used only when BaseMotionType=2]
0               NacYaw(1)              - Yaw angle (about z_t) of the nacelle (deg)
{IEA15MW_RATED_RPM:8.2f}            RotSpeed(1)            - Rotational speed of rotor in rotor coordinates (rpm)
0               BldPitch(1)            - Blade 1 pitch (deg)
----- Time-dependent Analysis [used only when AnalysisType=2, numTurbines=1] ------------
"unused"         TimeAnalysisFileName - Filename containing time series (6 column: Time, HWndSpeed, PLExp, RotSpd, Pitch, Yaw).
-----  Combined-Case Analysis [used only when AnalysisType=3, numTurbines=1 -------------
{len(cases):10d}  NumCases     - Number of cases to run
HWndSpeed  PLExp  RotSpd  Pitch   Yaw   dT    Tmax  DOF  Amplitude Frequency
(m/s)      (-)    (rpm)   (deg)  (deg)  (s)   (s)   (-)   (-)       (Hz)
{rows}
----- Output Settings -------------------------------------------------------------------
"ES15.8E2"       OutFmt      - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)
1                OutFileFmt  - Format for tabular (time-marching) output file (switch) {{1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}}
0                WrVTK       - VTK visualization data output: (switch) {{0=none; 1=init; 2=animation}}
1                WrVTK_Type  - VTK visualization data type: (switch) {{1=surfaces; 2=lines; 3=both}}
2                VTKHubRad   - HubRadius for VTK visualization (m)
-1,-1,-1,2,2,2   VTKNacDim   - Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)
"""
    path = Path(path).resolve()
    if not path.parent.is_dir() or path.name in ("", ".", ".."):
        raise ValueError(f"refusing to write a driver file to {path}")
    path.write_text(text)
    return path


@dataclass
class AeroDynResult:
    """Parsed ``aerodyn_driver`` output for one operating point."""

    time: float
    wind_speed: float
    rpm: float
    pitch_deg: float
    thrust: float  # RtFldFxh [N]
    torque: float  # RtFldMxh [N*m] (about the rotor axis)
    cp: float
    ct: float
    r: np.ndarray  # radial station [m]
    alpha_deg: np.ndarray
    cl: np.ndarray
    cd: np.ndarray
    cn: np.ndarray
    ct_section: np.ndarray
    fn: np.ndarray  # normal force per length [N/m]
    ft: np.ndarray  # tangential force per length [N/m]
    ax_ind: np.ndarray  # axial induction a
    tn_ind: np.ndarray  # tangential induction a'


def parse_aerodyn_out(path: str | Path) -> List[AeroDynResult]:
    """Parse a whitespace AeroDyn driver ``.out`` file into per-case results.

    The last row of each combined case is used (the driver writes a converged
    steady row at the end of every case; DBEMT transients are discarded).
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    header_idx = next(
        i for i, ln in enumerate(lines) if ln.startswith(" Time") or "\tTime" in ln
    )
    header = lines[header_idx].split()

    data_lines = []
    for ln in lines[header_idx + 2 :]:
        s = ln.strip()
        if not s or s.startswith("-"):
            continue
        fields = s.split()
        if len(fields) == len(header):
            try:
                data_lines.append([float(x) for x in fields])
            except ValueError:
                continue
    data = np.array(data_lines)

    col = {name: j for j, name in enumerate(header)}
    n_nodes = sum(1 for name in header if name.startswith("AB1N") and name.endswith("Alpha"))

    def _revolution_mean(rows: np.ndarray, names: Sequence[str]) -> np.ndarray:
        """Mean over the last full rotor revolution (azimuth wrap to wrap).

        CCBlade ``evaluate`` integrates across azimuth, so a yawed or sheared
        case must be compared against an azimuth-averaged AeroDyn value, not
        against one instantaneous row.
        """
        az = rows[:, col["Azimuth"]]
        wraps = np.where(np.diff(az) < 0)[0]
        segment = rows[wraps[-2] + 1 :] if len(wraps) >= 2 else rows
        return segment[:, [col[n] for n in names]].mean(axis=0)

    results: List[AeroDynResult] = []
    for case_id in np.unique(data[:, col["Case"]]):
        rows = data[data[:, col["Case"]] == case_id]
        last = rows[-1]
        thrust, torque, cp, ct = _revolution_mean(
            rows, ["RtFldFxh", "RtFldMxh", "RtFldCp", "RtFldCt"]
        )
        alpha = np.array([last[col[f"AB1N{i:03d}Alpha"]] for i in range(1, n_nodes + 1)])

        # Radial station: HubRad + BlSpn.  The driver does not print BlSpn, so
        # the caller supplies the station radii from the blade file; here we
        # return the raw per-node arrays and leave r to the caller.
        def arr(suffix: str) -> np.ndarray:
            return np.array([last[col[f"AB1N{i:03d}{suffix}"]] for i in range(1, n_nodes + 1)])

        results.append(
            AeroDynResult(
                time=float(last[col["Time"]]),
                wind_speed=float(last[col["HWindSpeedX"]]),
                rpm=float(last[col["RotSpeed"]]),
                pitch_deg=float(last[col["BldPitch1"]]),
                thrust=float(thrust),
                torque=float(torque),
                cp=float(cp),
                ct=float(ct),
                r=np.array([]),
                alpha_deg=alpha,
                cl=arr("Cl") if "AB1N001Cl" in col else np.array([]),
                cd=arr("Cd") if "AB1N001Cd" in col else np.array([]),
                cn=arr("Cn"),
                ct_section=arr("Ct"),
                fn=arr("Fn"),
                ft=arr("Ft"),
                ax_ind=arr("AxInd"),
                tn_ind=arr("TnInd"),
            )
        )
    return results


def run_aerodyn_driver(driver: Path, dvr: Path, workdir: Path) -> List[Path]:
    """Run ``aerodyn_driver`` in ``workdir`` and return the produced ``.out`` files.

    In combined-case mode (``AnalysisType=3``) the driver writes one file per
    case, named ``<root>.<case>.out``; every file is returned, sorted by case.
    """
    proc = subprocess.run(
        [str(driver), dvr.name],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=1800,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"aerodyn_driver failed (exit {proc.returncode})\n"
            f"--- stdout ---\n{proc.stdout[-4000:]}\n--- stderr ---\n{proc.stderr[-2000:]}"
        )
    outs = sorted(workdir.glob("*.out"))
    if not outs:
        raise RuntimeError(f"aerodyn_driver produced no .out in {workdir}")
    return outs


def parse_aerodyn_outs(paths: Sequence[str | Path]) -> List[AeroDynResult]:
    """Parse every ``.out`` file returned by :func:`run_aerodyn_driver`."""
    results: List[AeroDynResult] = []
    for p in paths:
        results.extend(parse_aerodyn_out(p))
    return results
