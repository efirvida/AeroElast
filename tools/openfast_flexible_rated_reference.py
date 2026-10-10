#!/usr/bin/env python3
"""Regenerate the durable OpenFAST flexible-torsion rated reference (P2b / T3).

Work item: ``odd/tasks/openfast-flexible-coupled-arbiter.md`` T3.

The ElastoDyn deck in ``ofruns/`` has no blade torsion degree of freedom
(``CompElast = 1``: two flapwise and one edgewise mode, no twist), so it cannot
arbitrate a twist in either direction.  Blade torsion needs BeamDyn
(``CompElast = 2``).  This script materialises a durable, self-contained BeamDyn
case under ``$SCRATCH/bfs16/openfast-flexible/``, forces the BeamDyn per-node
rotation outputs on in a *case-local* copy of ``IEA-15-240-RWT_BeamDyn.dat``,
runs OpenFAST, and writes the blade-1 spanwise torsion profile.

What is copied, never edited in place
-------------------------------------
``$SCRATCH/tmp/opencode/ofrun-s6-bd/`` (the 2026-09-15 Beamdyn run) and
``$SCRATCH/IEA-15-240-RWT/`` are read-only inputs.  The script copies out of
them; it never moves or rewrites them.

Model tree choice
-----------------
The BeamDyn case reaches its shared model tree as ``../IEA-15-240-RWT/``.  The
sibling next to the source case is used, because it carries ``HWindSpeed =
10.59 m/s`` -- the rated inflow.  The ``$SCRATCH/IEA-15-240-RWT/OpenFAST/``
clone carries ``HWindSpeed = 10.0``, which would move the operating point off
rated and away from ``bfs16/case5s`` (the A/B harness of issue #16, wind 10.59).

What is measured
----------------
Blade 1, BeamDyn nodal channel ``B1N###_RDxr`` -- the rotation about the local
blade (span) axis.  The BeamDyn registry documents this channel as
"rotational displacement in X, rad" (string found in
``conda-envs/openfast/lib/libbeamdynlib.a``); the OpenFAST ``.out`` unit row
prints ``-`` for it.  Values are reported raw: no sign, frame or unit
conversion is applied here.  Sign conventions belong to
``odd/tasks/rated-twist-sign-convention.md``.

Station definition
------------------
BeamDyn reports its output nodes as *quadrature points*
(``Output_nodes_location: "quadrature points"`` in
``<RootName>.BD.R1.B1.sum.yaml``); this mesh is 51 points for the
IEA-15-240-RWT BeamDyn blade.  Their undisplaced reference positions are read
from the ``vtk/<RootName>.BD_Blade_R1B1_Reference.vtp`` initialization mesh
(``WrVTK = 1``, ``VTK_type = 3``), and ``span_m`` is the cumulative straight-line
distance along that node polyline, origin at output node 1.  The total
(117.1487 m) reproduces the ``Length: 117.149 m`` the BeamDyn summary reports
for the same blade.

Usage
-----
    scripts/aeroenv.sh python tools/openfast_flexible_rated_reference.py both
    scripts/aeroenv.sh python tools/openfast_flexible_rated_reference.py short
    scripts/aeroenv.sh python tools/openfast_flexible_rated_reference.py full --skip-run
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

SCRATCH = Path("/scratch/leahk/eduardo.donestevez")
OPENFAST_BIN = SCRATCH / "conda-envs" / "openfast" / "bin" / "openfast"
SOURCE_CASE_DIR = SCRATCH / "tmp" / "opencode" / "ofrun-s6-bd" / "IEA-15-240-RWT-Monopile"
SOURCE_MODEL_DIR = SCRATCH / "tmp" / "opencode" / "ofrun-s6-bd" / "IEA-15-240-RWT"
DURABLE_ROOT = SCRATCH / "bfs16" / "openfast-flexible"

MODEL_DIRNAME = "IEA-15-240-RWT"
ROOTNAME = "IEA-15-240-RWT-Monopile"
# The case dir must keep the deck's own root name: ElastoDyn reaches its tower
# properties through the self-referential `../IEA-15-240-RWT-Monopile/...` path.
CASE_DIRNAME = ROOTNAME
FST_NAME = f"{ROOTNAME}.fst"
LOCAL_BEAMDYN_NAME = f"{ROOTNAME}_BeamDyn.dat"
BEAMDYN_BLADE_NAME = "IEA-15-240-RWT_BeamDyn_blade.dat"
# OpenFAST writes VTK initialization meshes under <case>/vtk/; this one carries
# the BeamDyn blade-1 output-node reference positions.
VTK_REFERENCE_RELPATH = f"vtk/{ROOTNAME}.BD_Blade_R1B1_Reference.vtp"
PROFILE_NAME = "blade1_torsion_profile.csv"
PROFILE_NAME_SHORT = "blade1_torsion_profile_short.csv"

# BeamDyn nodal output channels requested by this reference.  RDxr (rotation
# about the local blade axis) is the torsion channel; the rest are kept so the
# written file is self-describing and the probe can be re-read.
NODAL_OUT_CHANNELS = ("TDxr", "TDyr", "RDxr", "RDyr", "RDzr")
TORSION_CHANNEL = "RDxr"
# See configure_case(): rejected by this OpenFAST build and absent from the
# recorded output of the source run.
UNAVAILABLE_AERODYN_CHANNELS = ("B1Mp",)

WINDOW_START = 90.0
WINDOW_END = 100.0
SHORT_TMAX = 3.0
FULL_TMAX = 100.0

SKIP_SUFFIXES = {".out", ".outb", ".ech", ".sum", ".log", ".yaml", ".IN", ".pyc"}
SKIP_NAME_PREFIXES = ("bd_smoke", "s6bd")

_QUOTED = re.compile(r'^"([^"]*)"$')


def _log(message: str) -> None:
    print(message, flush=True)


# --------------------------------------------------------------------------- #
# materialisation
# --------------------------------------------------------------------------- #
def _copy_inputs(src: Path, dst: Path) -> None:
    """Copy every input file under ``src`` into ``dst`` (never deletes)."""
    for path in sorted(src.rglob("*")):
        rel = path.relative_to(src)
        if path.is_dir():
            (dst / rel).mkdir(parents=True, exist_ok=True)
            continue
        if path.suffix in SKIP_SUFFIXES:
            continue
        if path.name.startswith(SKIP_NAME_PREFIXES):
            continue
        (dst / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dst / rel)


def materialize(root: Path) -> tuple[Path, Path]:
    """Create the durable case dir and its ``../IEA-15-240-RWT/`` sibling."""
    case_dir = root / CASE_DIRNAME
    model_dir = root / MODEL_DIRNAME
    case_dir.mkdir(parents=True, exist_ok=True)
    model_dir.mkdir(parents=True, exist_ok=True)
    _copy_inputs(SOURCE_CASE_DIR, case_dir)
    _copy_inputs(SOURCE_MODEL_DIR, model_dir)
    (root / "runs").mkdir(parents=True, exist_ok=True)
    return case_dir, model_dir


# --------------------------------------------------------------------------- #
# deck patching (case-local copies only)
# --------------------------------------------------------------------------- #
def _set_scalar(text: str, key: str, value: str) -> str:
    """Replace the value token of a single ``value  Key  - comment`` line."""
    pattern = re.compile(
        rf"(?m)^(?P<indent>[ \t]*)(?P<val>\S+)(?P<tail>[ \t]+{re.escape(key)}(?![A-Za-z0-9_]).*)$"
    )

    def _sub(match: re.Match[str]) -> str:
        return f"{match.group('indent')}{value:<12}{match.group('tail')}"

    patched, count = pattern.subn(_sub, text)
    if count != 1:
        raise RuntimeError(f"expected exactly one {key!r} line, found {count}")
    return patched


def _ensure_nodal_outlist(text: str, channels: tuple[str, ...]) -> str:
    """Add any missing channel to the optional BeamDyn nodal-out section."""
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if "Outputs for all blade stations" in line)
    end = next(
        i for i in range(start + 1, len(lines)) if lines[i].strip().upper().startswith("END")
    )
    present = [
        match.group(1).strip()
        for line in lines[start:end]
        if (match := _QUOTED.match(line.strip()))
    ]
    missing = [channel for channel in channels if channel not in present]
    if not missing:
        return text
    patched = lines[:end] + [f'"{channel}"' for channel in missing] + lines[end:]
    return "\n".join(patched) + "\n"


def _drop_outlist_channels(text: str, channels: tuple[str, ...]) -> str:
    """Remove ``"channel"`` lines from an OpenFAST OutList section."""
    wanted = {f'"{channel}"' for channel in channels}
    kept = [line for line in text.splitlines() if line.strip() not in wanted]
    return "\n".join(kept) + "\n"


def configure_case(case_dir: Path, model_dir: Path, n_node_outs: int) -> Path:
    """Write the case-local module decks and repoint the ``.fst`` at them."""
    source_beamdyn = model_dir / "IEA-15-240-RWT_BeamDyn.dat"
    text = source_beamdyn.read_text()
    # OpenFAST needs NNodeOuts >= 1 before the optional nodal section is written;
    # with BldNd_BlOutNd = "All" the run emits every blade station (51 for the
    # IEA-15-240-RWT BeamDyn mesh), not just NNodeOuts of them.
    text = _set_scalar(text, "NNodeOuts", str(n_node_outs))
    text = _set_scalar(text, "OutNd", "1")
    text = _set_scalar(text, "BldNd_BlOutNd", '"All"')
    text = _set_scalar(text, "SumPrint", "True")
    text = _ensure_nodal_outlist(text, NODAL_OUT_CHANNELS)

    local_beamdyn = case_dir / LOCAL_BEAMDYN_NAME
    local_beamdyn.write_text(text)
    # The local deck names its BldFile with a bare relative path, so the blade
    # properties file has to sit beside it.
    shutil.copy2(model_dir / BEAMDYN_BLADE_NAME, case_dir / BEAMDYN_BLADE_NAME)

    fst_path = case_dir / FST_NAME
    fst = fst_path.read_text()
    for blade in (1, 2, 3):
        fst = _set_scalar(fst, f"BDBldFile({blade})", f'"{LOCAL_BEAMDYN_NAME}"')
    fst = _set_scalar(fst, "SumPrint", "True")
    # WrVTK = 1 writes the initialization meshes only; VTK_type = 3 (all meshes)
    # is what emits the BeamDyn blade mesh, which is how the 51 output-node
    # reference positions are obtained (the <RootName>.BD.R1.B1.sum.yaml position
    # tables are emitted empty by this build).  Output selection only -- no effect
    # on the solution.
    fst = _set_scalar(fst, "WrVTK", "1")
    fst = _set_scalar(fst, "VTK_type", "3")
    fst_path.write_text(fst)

    # The copied source deck carries an AeroDyn OutList entry this OpenFAST build
    # rejects ("B1Mp is not an available output channel").  It was added to the
    # source after the recorded run: the archived <RootName>.out channel list has
    # no B1Mp channel.  Dropping it in the case-local copy restores a runnable
    # deck without touching the source or the physics.
    aerodyn_path = case_dir / f"{ROOTNAME}_AeroDyn15.dat"
    aerodyn_path.write_text(
        _drop_outlist_channels(aerodyn_path.read_text(), UNAVAILABLE_AERODYN_CHANNELS)
    )
    return local_beamdyn


def set_tmax(case_dir: Path, tmax: float) -> None:
    fst_path = case_dir / FST_NAME
    fst_path.write_text(_set_scalar(fst_path.read_text(), "TMax", f"{tmax}"))


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
def run_openfast(case_dir: Path, log_path: Path) -> float:
    command = [str(OPENFAST_BIN), FST_NAME]
    _log(f"$ cd {case_dir} && {' '.join(command)}")
    start = time.monotonic()
    with log_path.open("w") as handle:
        proc = subprocess.run(
            command,
            cwd=case_dir,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    wall = time.monotonic() - start
    log = log_path.read_text(errors="replace")
    if "terminated normally" not in log:
        _log(f"--- tail of {log_path} ---")
        _log("\n".join(log.splitlines()[-15:]))
        raise RuntimeError(f"OpenFAST did not terminate normally (exit {proc.returncode})")
    return wall


# --------------------------------------------------------------------------- #
# extraction
# --------------------------------------------------------------------------- #
def _load_out(outb_path: Path) -> tuple[list[str], list[str], np.ndarray]:
    from openfast_toolbox.io import FASTOutputFile

    out = FASTOutputFile(str(outb_path))
    return list(out.channels), list(out.units), np.asarray(out.data)


def _index(names: list[str], channel: str) -> int:
    if channel in names:
        return names.index(channel)
    raise KeyError(f"channel {channel!r} is not in the written output file")


def _node_number(channel: str) -> int:
    match = re.search(r"B1N(\d+)", channel)
    if match is None:
        raise ValueError(f"cannot read a node index out of {channel!r}")
    return int(match.group(1))


def _torsion_channels(names: list[str]) -> list[str]:
    found = [
        name
        for name in names
        if re.fullmatch(rf"B1N(\d+)_?{TORSION_CHANNEL}", name)
        or re.fullmatch(rf"B1N(\d+){TORSION_CHANNEL}", name)
    ]
    if not found:
        raise RuntimeError(
            f"no blade-1 {TORSION_CHANNEL} nodal channel in the written file "
            f"(found {len(names)} channels) -- refusing to report a reference"
        )
    return sorted(found, key=_node_number)


def _invalid(units: list[str], names: list[str], channel: str) -> bool:
    if channel not in names:
        return True
    label = units[names.index(channel)].upper()
    return "NVALI" in label or "INVALID" in label


@dataclass
class Extraction:
    """Raw channel data read back out of a written OpenFAST output file."""

    names: list[str]
    units: list[str]
    data: np.ndarray
    mask: np.ndarray
    channels: list[str]
    samples: int
    window: tuple[float, float]
    stations: list[tuple[float, float, str, str]]
    station_label: str
    total_span: float


def _spans_from_vtk(vtk_path: Path, n_nodes: int) -> tuple[list[float], float, str]:
    """Read the BeamDyn blade-1 output-node reference positions from a ``.vtp``.

    OpenFAST writes the initialization meshes (``WrVTK = 1``); ``VTK_type = 3``
    is what emits the BeamDyn blade mesh.  Its point list is exactly the 51
    quadrature points that carry the ``B1N###_`` channels, in node order.

    Returns per-node span [m] (cumulative straight-line distance along the node
    reference polyline, origin at output node 1), the total arc length, and a
    label describing the source.
    """
    if not vtk_path.exists():
        raise RuntimeError(f"missing {vtk_path}; cannot define the blade stations")
    block = re.search(
        r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>",
        vtk_path.read_text(errors="replace"),
        re.DOTALL,
    )
    if block is None:
        raise RuntimeError(f"no <Points> array in {vtk_path}")
    values = [float(token) for token in block.group(1).split()]
    points = [values[i : i + 3] for i in range(0, len(values), 3)]
    if len(points) != n_nodes:
        raise RuntimeError(f"{vtk_path.name}: has {len(points)} points, expected {n_nodes}")
    spans = [0.0]
    for previous, current in zip(points[:-1], points[1:], strict=True):
        step = sum((b - a) ** 2 for a, b in zip(previous, current, strict=True)) ** 0.5
        spans.append(spans[-1] + step)
    label = (
        f"BeamDyn blade-1 reference mesh {vtk_path.name} ({len(points)} output nodes); "
        "span = cumulative arc length along the node reference polyline, origin at B1N001"
    )
    return spans, spans[-1], label


def _fallback_spans(n_nodes: int) -> tuple[list[float], float, str]:
    blade_span_m = 117.0
    spacing = blade_span_m / (n_nodes - 1)
    label = (
        f"uniform spacing over the {blade_span_m:.1f} m blade span "
        "(no BeamDyn node-position mesh available)"
    )
    return [index * spacing for index in range(n_nodes)], blade_span_m, label


def write_profile(
    outb_path: Path,
    vtk_path: Path,
    csv_path: Path,
    window: tuple[float, float],
    tmax: float,
) -> "Extraction":
    names, units, data = _load_out(outb_path)
    time_index = _index(names, "Time")
    times = data[:, time_index]
    mask = (times >= window[0]) & (times <= window[1])
    samples = int(mask.sum())
    if samples == 0:
        raise RuntimeError(f"window {window} contains no samples (t max {times.max()})")

    channels = _torsion_channels(names)
    try:
        spans, total_span, station_label = _spans_from_vtk(vtk_path, len(channels))
    except RuntimeError as error:
        _log(f"station mesh unavailable ({error}); falling back to uniform spacing")
        spans, total_span, station_label = _fallback_spans(len(channels))

    header = [
        "# OpenFAST flexible-torsion rated reference -- blade 1, BeamDyn nodal rotation",
        f"# source deck : {DURABLE_ROOT / CASE_DIRNAME / FST_NAME}",
        f"# openfast    : {OPENFAST_BIN} (v5.0.0, single precision)",
        f"# run length  : TMax = {tmax} s, DT_Out = 0.05 s",
        f"# window      : {window[0]:.3f} .. {window[1]:.3f} s ({samples} samples)",
        f"# channel     : B1N###_{TORSION_CHANNEL} (BeamDyn: rotational displacement in X, rad)",
        f"# stations    : {len(channels)} nodes, total arc {total_span:.4f} m; {station_label}",
        "# raw values  : no sign, frame or unit conversion applied",
        "# units       : 'rad' per the BeamDyn registry; the .out unit row prints '-'",
    ]
    stations = [
        (span, float(data[mask, _index(names, channel)].mean()), "rad", channel)
        for span, channel in zip(spans, channels, strict=True)
    ]

    with csv_path.open("w", newline="") as handle:
        handle.write("\n".join(header))
        handle.write("\n")
        handle.write("span_m,mean_value,unit,channel\n")
        for span, mean_value, unit, channel in stations:
            handle.write(f"{span:.6f},{mean_value:.10f},{unit},{channel}\n")

    return Extraction(
        names=names,
        units=units,
        data=data,
        mask=mask,
        channels=channels,
        samples=samples,
        window=window,
        stations=stations,
        station_label=station_label,
        total_span=total_span,
    )


def report_sanity(extracted: "Extraction") -> dict[str, float | None]:
    names = extracted.names
    units = extracted.units
    data = extracted.data
    mask = extracted.mask

    def mean_of(channel: str) -> float | None:
        if _invalid(units, names, channel):
            return None
        return float(data[mask, names.index(channel)].mean())

    rotor_speed = mean_of("RotSpeed")
    gen_speed = mean_of("GenSpeed")
    pitch = mean_of("BldPitch1")
    torque = mean_of("RotTorq")
    power_mw = None
    if torque is not None and rotor_speed is not None:
        power_mw = torque * rotor_speed * 2.0 * math.pi / 60.0 / 1000.0

    tip = mean_of("TipDxc1")
    tip_state = "n/a (INVALID unit with CompElast=2)" if tip is None else f"{tip:.3f} m"

    _log("")
    _log("--- sanity block -------------------------------------------------")
    _log(f"rotor speed RotSpeed        : {rotor_speed} rpm")
    _log(f"generator speed GenSpeed    : {gen_speed} rpm")
    _log(f"blade pitch BldPitch1       : {pitch} deg")
    _log(
        "rotor aero power RotTorq*RotSpeed : "
        f"{power_mw:.3f} MW  (CompServo=0 -> no generator model, no GenPwr channel)"
    )
    _log(f"TipDxc1 (S-5 record 16.03 m) : {tip_state}")
    _log("------------------------------------------------------------------")
    return {"rotor_speed": rotor_speed, "pitch": pitch, "power_mw": power_mw}


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def _archive(case_dir: Path, root: Path, mode: str) -> None:
    dest = root / "runs" / mode
    dest.mkdir(parents=True, exist_ok=True)
    for name in (f"{ROOTNAME}.outb", f"{ROOTNAME}.BD.R1.B1.sum.yaml", VTK_REFERENCE_RELPATH):
        src = case_dir / name
        if src.exists():
            (dest / src.name).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest / src.name)


def run_mode(root: Path, mode: str, window: tuple[float, float], skip_run: bool) -> None:
    case_dir = root / CASE_DIRNAME
    tmax = SHORT_TMAX if mode == "short" else FULL_TMAX
    if mode == "short":
        window = (0.0, tmax)
    set_tmax(case_dir, tmax)
    log_path = root / f"run_{mode}.log"
    wall = None
    if skip_run:
        _log(f"[{mode}] --skip-run: reusing the existing outputs in {case_dir}")
    else:
        wall = run_openfast(case_dir, log_path)
        _log(f"[{mode}] wall time: {wall:.1f} s ({wall / 60.0:.2f} min)")

    csv_path = root / (PROFILE_NAME_SHORT if mode == "short" else PROFILE_NAME)
    extracted = write_profile(
        case_dir / f"{ROOTNAME}.outb",
        case_dir / VTK_REFERENCE_RELPATH,
        csv_path,
        window,
        tmax,
    )
    _log(f"[{mode}] wrote {csv_path}")
    _log(f"[{mode}] torsion channels: {len(extracted.channels)} blade-1 nodes")
    _log(f"[{mode}] stations: total arc {extracted.total_span:.4f} m; {extracted.station_label}")
    _log(f"[{mode}] first station: {extracted.stations[0]}")
    _log(f"[{mode}] last  station: {extracted.stations[-1]}")
    report_sanity(extracted)
    if not skip_run:
        _archive(case_dir, root, mode)


def main(argv: list[str] | None = None) -> int:
    description = (__doc__ or "OpenFAST flexible-torsion reference").splitlines()[0]
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "mode",
        nargs="?",
        default="both",
        choices=("short", "full", "both"),
        help="short validation run, full 100 s run, or both (default)",
    )
    parser.add_argument(
        "--window",
        type=float,
        nargs=2,
        metavar=("T0", "T1"),
        default=(WINDOW_START, WINDOW_END),
    )
    parser.add_argument("--n-node-outs", type=int, default=1)
    parser.add_argument("--skip-run", action="store_true", help="reuse existing outputs")
    parser.add_argument("--root", type=Path, default=DURABLE_ROOT)
    args = parser.parse_args(argv)

    if not OPENFAST_BIN.exists():
        _log(f"BLOCKED: OpenFAST binary missing at {OPENFAST_BIN}")
        return 2
    if not SOURCE_CASE_DIR.is_dir():
        _log(f"BLOCKED: source BeamDyn case missing at {SOURCE_CASE_DIR}")
        return 2
    if not SOURCE_MODEL_DIR.is_dir():
        _log(f"BLOCKED: source model tree missing at {SOURCE_MODEL_DIR}")
        return 2

    case_dir, model_dir = materialize(args.root)
    local_beamdyn = configure_case(case_dir, model_dir, args.n_node_outs)
    _log(f"durable root : {args.root}")
    _log(f"case dir     : {case_dir}")
    _log(f"model tree   : {model_dir}")
    _log(f"case-local BeamDyn deck: {local_beamdyn}")

    window = (args.window[0], args.window[1])
    modes = {"short": ["short"], "full": ["full"], "both": ["short", "full"]}[args.mode]
    for mode in modes:
        run_mode(args.root, mode, window, args.skip_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
