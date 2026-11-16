#!/usr/bin/env python3
"""L5 two-way (thesis plan step 5) OpenFAST v5.0.0 reference decks: flexible blades and rigid twin.

Generates, under ``--out`` (default ``$SCRATCH/l5_twoway/openfast``)::

    common/         copy of the IEA-15-240-RWT deck files shared by both decks (airfoils, AD blade, ED blade,
                    ED tower, BD main file, BD blade file) -- inputs are copied, never edited in place
    flex/           ElastoDyn + BeamDyn + AeroDyn + InflowWind, flexible blades (BeamDyn 6x6 as in the deck)
    rigid/          identical deck, BeamDyn stiffness 6x6 multiplied by --rigid-scale (default 1e3): same blade
                    geometry (prebend, cant, twist), same mass, negligible deformation (the flexible-minus-rigid twin)
    manifest.json   sha256 of every source file and every generated file, option table, generator version

Source deck (override with --deck-root): ``ofruns/OpenFAST/{IEA-15-240-RWT,IEA-15-240-RWT-Monopile}``, the
v4.2 -> v5.0.0 conversion that S-4/S-5 and the L4 closure use (v5 AeroDyn layout with TwrCp/TwrCa already
applied; see Engram ``config/openfast-5-0-0-corriendo-iea-15-mw-deck-v4-2-v5-0-convertido``).

Options (all matched to the L4 closure, odd/tasks/l4-bem-rigid-closure.md L5 "T5 matched options", and to the
production two-way case cases/l5_twoway_rated): rated 10.59 m/s uniform (PLexp 0), 7.56 rpm fixed (GenDOF/DrTrDOF
off, ElastoDyn RotSpeed), pitch 0, tilt -6, precone -4 (deck ElastoDyn), tower influence/shadow/aero off, UA_Mod 0,
DBEMT_Mod 0, Wake_Mod 1, BEM_Mod 1, Skew_Mod 1 (deck), tip/hub loss on, gravity off, CompServo 0.

No numbers are produced here: the smoke (``--smoke``) only proves the decks run.  Usage::

    scripts/aeroenv.sh python tools/l5_twoway_openfast_deck.py                 # generate
    scripts/aeroenv.sh python tools/l5_twoway_openfast_deck.py --smoke 2.0     # generate + 2 s run of each deck
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRATCH = Path(os.environ.get("SCRATCH", "/scratch/leahk/eduardo.donestevez"))
DEFAULT_DECK_ROOT = Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST")
DEFAULT_OPENFAST = Path("/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast")
GENERATOR_VERSION = "1"
DECK_FILES = {"run.fst", "ElastoDyn.dat", "AeroDyn15.dat", "BeamDyn.dat", "IEA-15-240-RWT_BeamDyn_blade.dat"}

RPM_RATED = 7.56
WIND = 10.59

COMMON_SRC = "IEA-15-240-RWT"
MONO_SRC = "IEA-15-240-RWT-Monopile"
# files copied verbatim into common/ (rel. to deck root) -- hashes recorded as the deck inputs
COMMON_FILES = [
    f"{COMMON_SRC}/IEA-15-240-RWT_AeroDyn15_blade.dat",
    f"{COMMON_SRC}/IEA-15-240-RWT_ElastoDyn_blade.dat",
    f"{COMMON_SRC}/IEA-15-240-RWT_BeamDyn.dat",
    f"{COMMON_SRC}/IEA-15-240-RWT_BeamDyn_blade.dat",
    f"{COMMON_SRC}/IEA-15-240-RWT_InflowFile.dat",
    f"{MONO_SRC}/IEA-15-240-RWT-Monopile.fst",
    f"{MONO_SRC}/IEA-15-240-RWT-Monopile_ElastoDyn.dat",
    f"{MONO_SRC}/IEA-15-240-RWT-Monopile_ElastoDyn_tower.dat",
    f"{MONO_SRC}/IEA-15-240-RWT-Monopile_AeroDyn15.dat",
]

AD_OVERRIDES = {  # same switches as tests/support/openfast_bem.py write_bem_primary (L4 harness) + explicit pins
    "UA_Mod": "0",
    "TwrPotent": "0",
    "TwrShadow": "0",
    "TwrAero": "False",
    "NacelleDrag": "False",
    "Wake_Mod": "1",
    "BEM_Mod": "1",
    "DBEMT_Mod": "0",
}
AD_EXTRA_OUT = ["RtAeroFxh", "RtAeroFyh", "RtAeroFzh", "RtAeroMxh", "RtAeroMyh", "RtAeroMzh", "RtAeroPwr", "RtAeroCp", "RtAeroCt"]
ED_OUT = ["Azimuth", "RotSpeed", "RotThrust", "RotTorq", "RotPwr", "LSShftFxa", "LSShftMxa"]


def sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def set_key(text: str, key: str, value: str) -> str:
    """Replace the value token(s) on the line whose label (2nd field) is ``key``; keep the description."""
    pat = re.compile(r"^(\s*)(\"[^\"]*\"|\S+)(\s+)(" + re.escape(key) + r")((?=\s).*|)$", re.M)
    new, n = pat.subn(lambda m: f"{m.group(1)}{value}{m.group(3)}{m.group(4)}{m.group(5)}", text, count=1)
    if n != 1:
        raise KeyError(f"key {key!r} not found")
    return new


def replace_outlist(text: str, channels: list[str]) -> str:
    """Replace the first OutList block (up to its END line) by ``channels`` (one quoted name per line)."""
    lines = text.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if re.match(r"^\s*OutList\b", ln))
    i1 = next(i for i in range(i0 + 1, len(lines)) if lines[i].startswith("END"))
    lines[i0 + 1 : i1] = [f'"{c}"' for c in channels]
    return "\n".join(lines) + "\n"


def extend_outlist(text: str, channels: list[str]) -> str:
    lines = text.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if re.match(r"^\s*OutList\b", ln))
    i1 = next(i for i in range(i0 + 1, len(lines)) if lines[i].startswith("END"))
    present = {ln.strip().strip('"') for ln in lines[i0 + 1 : i1]}
    lines[i1:i1] = [f'"{c}"' for c in channels if c not in present]
    return "\n".join(lines) + "\n"


def bd_blade_v5(text: str) -> str:
    """v4.2 -> v5.0.0 BeamDyn blade-file header: add the modal-damping block (BeamDyn_IO.f90:1054-1080).

    damp_type stays 1 (stiffness-proportional, deck mu); the inserted n_modes/zeta are unused unless damp_type = 2.
    """
    lines = text.splitlines()
    i = next(k for k, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln)
    extra = ["------ Modal Damping [used only if damp_type=2] --------------------------------",
             "1                        n_modes    - Number of modal damping coefficients (-)", "0.0"]
    return "\n".join(lines[:i] + extra + lines[i:]) + "\n"


def scale_bd_stiffness(text: str, scale: float) -> str:
    """Multiply every 6x6 stiffness matrix of a BeamDyn blade file; mass matrices and damping mu unchanged."""
    lines = text.splitlines()
    out, i = [], 0
    while i < len(lines):
        out.append(lines[i])
        if re.fullmatch(r"\s*[-+0-9.eE]+\s*", lines[i]) and i + 6 < len(lines) and len(lines[i + 1].split()) == 6:
            # station fraction line followed by 6 stiffness rows
            for r in range(1, 7):
                vals = [float(v) * scale for v in lines[i + r].split()]
                out.append("\t " + " \t ".join(f"{v:.16e}" for v in vals))
            i += 7
            continue
        i += 1
    return "\n".join(out) + "\n"


def write_deck(out: Path, deck_root: Path, flavour: str, tmax: float, dt: float, rigid_scale: float, gravity: float):
    d = out / flavour
    d.mkdir(parents=True, exist_ok=True)
    com = "../common"

    def rd(rel):
        return (out / "common" / Path(rel).name).read_text()

    # ---- main fst
    fst = rd(f"{MONO_SRC}/IEA-15-240-RWT-Monopile.fst")
    for k, v in {
        "TMax": f"{tmax}",
        "DT": f"{dt}",
        "CompElast": "2",
        "CompInflow": "1",
        "CompAero": "2",
        "CompServo": "0",
        "Gravity": f"{gravity}",
        "EDFile": '"ElastoDyn.dat"',
        "BDBldFile(1)": '"BeamDyn.dat"',
        "BDBldFile(2)": '"BeamDyn.dat"',
        "BDBldFile(3)": '"BeamDyn.dat"',
        "InflowFile": f'"{com}/IEA-15-240-RWT_InflowFile.dat"',
        "AeroFile": '"AeroDyn15.dat"',
    }.items():
        fst = set_key(fst, k, v)
    (d / "run.fst").write_text(fst)

    # ---- ElastoDyn: rotor speed fixed (GenDOF/DrTrDOF false in the deck), blade DOFs irrelevant under BeamDyn
    ed = rd(f"{MONO_SRC}/IEA-15-240-RWT-Monopile_ElastoDyn.dat")
    ed = set_key(ed, "RotSpeed", f"{RPM_RATED}")
    for k in ("GenDOF", "DrTrDOF", "YawDOF", "TwFADOF1", "TwSSDOF1", "FlapDOF1", "FlapDOF2", "EdgeDOF"):
        ed = set_key(ed, k, "False")
    for k in ("BldFile1", "BldFile2", "BldFile3"):
        ed = set_key(ed, k, f'"{com}/IEA-15-240-RWT_ElastoDyn_blade.dat"')
    ed = set_key(ed, "TwrFile", f'"{com}/IEA-15-240-RWT-Monopile_ElastoDyn_tower.dat"')
    ed = set_key(ed, "SumPrint", "False")
    ed = set_key(ed, "NTwGages", "0")
    ed = set_key(ed, "NBlGages", "0")
    ed = replace_outlist(ed, ED_OUT)
    (d / "ElastoDyn.dat").write_text(ed)

    # ---- AeroDyn
    ad = rd(f"{MONO_SRC}/IEA-15-240-RWT-Monopile_AeroDyn15.dat")
    for k, v in AD_OVERRIDES.items():
        ad = set_key(ad, k, v)
    ad = ad.replace("../IEA-15-240-RWT/Airfoils/", f"{com}/Airfoils/")
    ad = ad.replace("../IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat", f"{com}/IEA-15-240-RWT_AeroDyn15_blade.dat")
    ad = extend_outlist(ad, AD_EXTRA_OUT)
    (d / "AeroDyn15.dat").write_text(ad)

    # ---- BeamDyn (blade file local so the rigid twin can differ)
    bd = rd(f"{COMMON_SRC}/IEA-15-240-RWT_BeamDyn.dat")
    # v5.0.0 BeamDyn dropped the PITCH ACTUATOR section (header + UsePitchAct, PitchJ, PitchK, PitchC) of the
    # v4.2 deck layout (build/openfast-v5.0.0/modules/beamdyn/src/BeamDyn_IO.f90:905-928); rest unchanged
    bl = bd.splitlines()
    i = next(k for k, ln in enumerate(bl) if "PITCH ACTUATOR" in ln)
    j = next(k for k, ln in enumerate(bl) if "OUTPUTS" in ln and ln.startswith("-"))
    bd = "\n".join(bl[:i] + bl[j:]) + "\n"
    (d / "BeamDyn.dat").write_text(bd)
    blade = bd_blade_v5(rd(f"{COMMON_SRC}/IEA-15-240-RWT_BeamDyn_blade.dat"))
    if flavour == "rigid":
        blade = scale_bd_stiffness(blade, rigid_scale)
    (d / "IEA-15-240-RWT_BeamDyn_blade.dat").write_text(blade)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--deck-root", type=Path, default=DEFAULT_DECK_ROOT)
    ap.add_argument("--out", type=Path, default=SCRATCH / "l5_twoway" / "openfast")
    ap.add_argument("--openfast", type=Path, default=Path(os.environ.get("OPENFAST_BIN", DEFAULT_OPENFAST)))
    ap.add_argument("--tmax", type=float, default=60.0, help="s; mean over whole revolutions (T_rev 7.937 s) of the tail")
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--rigid-scale", type=float, default=1.0e3)
    ap.add_argument("--gravity", type=float, default=0.0, help="m/s^2; 0 matches cases/l5_twoway_rated (gravity off)")
    ap.add_argument("--smoke", type=float, default=0.0, help="if > 0: run each deck for this many s in <deck>/smoke/")
    a = ap.parse_args()

    out = a.out
    (out / "common" / "Airfoils").mkdir(parents=True, exist_ok=True)
    manifest = {"generator": Path(__file__).name, "version": GENERATOR_VERSION, "deck_root": str(a.deck_root), "sources": {}, "generated": {}}
    for rel in COMMON_FILES:
        src = a.deck_root / rel
        shutil.copy2(src, out / "common" / src.name)
        manifest["sources"][rel] = sha256(src)
    for src in sorted((a.deck_root / COMMON_SRC / "Airfoils").glob("*")):
        shutil.copy2(src, out / "common" / "Airfoils" / src.name)
        manifest["sources"][f"{COMMON_SRC}/Airfoils/{src.name}"] = sha256(src)

    # steady uniform inflow: the deck InflowWind has a 0.12 power-law shear; L4 and our BEM case use none
    inflow = (out / "common" / "IEA-15-240-RWT_InflowFile.dat")
    t = set_key(inflow.read_text(), "PLexp", "0.0")
    t = set_key(t, "HWindSpeed", f"{WIND}")
    inflow.write_text(t)

    for flavour in ("flex", "rigid"):
        d = write_deck(out, a.deck_root, flavour, a.tmax, a.dt, a.rigid_scale, a.gravity)
        for f in sorted(d.glob("*")):
            if f.name in DECK_FILES:
                manifest["generated"][f"{flavour}/{f.name}"] = sha256(f)
    manifest["generated"]["common/IEA-15-240-RWT_InflowFile.dat"] = sha256(inflow)
    manifest["options"] = {"tmax": a.tmax, "dt": a.dt, "rigid_scale": a.rigid_scale, "gravity": a.gravity, "rpm": RPM_RATED, "wind": WIND,
                           "aerodyn_overrides": AD_OVERRIDES}
    manifest["generator_sha256"] = sha256(Path(__file__))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
    print(f"decks written under {out}; manifest {out / 'manifest.json'}")

    if a.smoke > 0:
        rc = 0
        for flavour in ("flex", "rigid"):
            sm = out / flavour
            fst = set_key((sm / "run.fst").read_text(), "TMax", f"{a.smoke}")
            (sm / "smoke.fst").write_text(fst)
            p = subprocess.run([str(a.openfast), "smoke.fst"], cwd=sm, capture_output=True, text=True)
            (sm / "smoke.run.log").write_text(p.stdout + p.stderr)
            tail = [ln for ln in (p.stdout + p.stderr).splitlines() if ln.strip()][-4:]
            print(f"[smoke {flavour}] rc={p.returncode}: " + " | ".join(tail))
            rc |= p.returncode
        sys.exit(rc)


if __name__ == "__main__":
    main()
