"""OpenFAST Campbell reference: MBC3 linearizations at the sweep speeds.

Runs the blade-only OpenFAST deck at each rotor speed with the linearization
enabled, post-processes the .lin files with MBC3 (openfast_toolbox) and
writes the Campbell reference CSV (rpm, f1f, f2f, f1e) to compare against
the shell's diagram (tools/run_sx_campbell.py).

Usage (needs the OpenFAST env + the blade-only deck):
  python tools/run_sx_of_campbell.py --rpms 3 5 9 11 \
      --deck-dir /scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT-Monopile \
      --csv docs/validation_data/generated/sx_campbell_openfast.csv
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

OPENFAST_BIN = "/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast"
OPENFAST_LIB = "/scratch/leahk/eduardo.donestevez/conda-envs/openfast/lib"


def _edit_deck(deck_dir: Path, rpm: float) -> None:
    """Set RotSpeed + the linearization block (TMax=35, Linearize at t=30).

    Also disables the aero (CompAero=0) — the structural-only linearization
    matching the S-4 procedure.  The MBC3 periodic transformation needs the
    azimuth steps: NLinTimes=36 linearly over one revolution around t=30 s.
    """
    import re

    fst = deck_dir / "IEA-15-240-RWT-Monopile.fst"
    s = fst.read_text()
    t0 = 30.0
    period = 60.0 / rpm
    n_lin = 12   # the .fst LinTimes line must stay under the parser's
                 # character limit (~26 values max); 12 azimuth steps suffice
    lin_times = ", ".join(f"{t0 + k * period / n_lin:.6f}" for k in range(n_lin))
    tmax = t0 + period + 2.0
    s = re.sub(r"[0-9.]+(\s+)TMax", f"{tmax:<20}TMax", s, count=1)
    s = re.sub(r"False\s+Linearize", "True                   Linearize", s, count=1)
    s = re.sub(r"30\.000000, 60\.000000\s+LinTimes", f"{lin_times}   LinTimes", s, count=1)
    s = re.sub(r"^\s*2\s+NLinTimes", f"{n_lin}                  NLinTimes", s,
               count=1, flags=re.M)
    s = re.sub(r"^\s*2\s+CompAero", "0                      CompAero", s,
               count=1, flags=re.M)
    fst.write_text(s)

    ed = deck_dir / "IEA-15-240-RWT-Monopile_ElastoDyn.dat"
    s = ed.read_text()
    s = re.sub(r"^\s*\d+(\.\d+)?\s+RotSpeed", f"{rpm:<22} RotSpeed", s,
               count=1, flags=re.M)
    ed.write_text(s)


def run_one(deck_dir: Path, workdir: Path, rpm: float) -> list[Path]:
    import os
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = f"{OPENFAST_LIB}:{env.get('LD_LIBRARY_PATH', '')}"
    fst = deck_dir / "IEA-15-240-RWT-Monopile.fst"
    log_path = workdir / "openfast.log"
    with open(log_path, "w") as log:
        res = subprocess.run(
            [OPENFAST_BIN, str(fst)], cwd=str(workdir), env=env,
            stdout=log, stderr=subprocess.STDOUT, timeout=1500)
    if res.returncode != 0:
        tail = log_path.read_text()[-1200:]
        raise RuntimeError(f"openfast failed for rpm={rpm}:\n{tail}")
    # the .lin files land NEXT TO THE .fst (not the CWD)
    lin_files = sorted(deck_dir.glob("*.lin"))
    return lin_files


def mbc3_frequencies(lin_files: list[Path]):
    """MBC3 collective natural frequencies (Hz) from the azimuth-averaged A."""
    from openfast_toolbox.linearization import mbc

    mbc_out, _ = mbc.fx_mbc3([str(f) for f in lin_files], verbose=False)
    eig = mbc_out["eigSol"]
    freq_rads = np.asarray(eig["NaturalFrequencies"], dtype=float)
    return freq_rads / (2 * np.pi)


def classify_modes(freqs):
    """Keep the positive real frequencies, sorted ascending."""
    f = np.asarray(freqs)
    keep = f > 0.05
    return np.sort(f[keep])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rpms", type=float, nargs="+", default=[3.0, 5.0, 9.0, 11.0])
    ap.add_argument("--deck-dir", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/"
                                 "IEA-15-240-RWT-Monopile"))
    ap.add_argument("--workdir", type=Path,
                    default=Path("/scratch/leahk/eduardo.donestevez/of_campbell"))
    ap.add_argument("--csv", type=Path,
                    default=Path("docs/validation_data/generated/sx_campbell_openfast.csv"))
    args = ap.parse_args()

    args.workdir.mkdir(parents=True, exist_ok=True)

    # the deck edits are destructive — back up and restore afterwards
    fst = args.deck_dir / "IEA-15-240-RWT-Monopile.fst"
    ed = args.deck_dir / "IEA-15-240-RWT-Monopile_ElastoDyn.dat"
    fst_backup = fst.read_text()
    ed_backup = ed.read_text()

    rows = []
    try:
        for rpm in args.rpms:
            print(f"\n=== rpm={rpm} ===", flush=True)
            _edit_deck(args.deck_dir, rpm)
            (args.workdir / f"rpm_{rpm:g}").mkdir(exist_ok=True)
            lin_files = run_one(args.deck_dir, args.workdir / f"rpm_{rpm:g}", rpm)
            freqs = classify_modes(mbc3_frequencies(lin_files))
            print(f"  freqs: {np.round(freqs[:6], 4)}", flush=True)
            rows.append((rpm, freqs))
            # move the .lin files into the per-speed workdir and clear the deck dir
            for f in lin_files:
                shutil.move(str(f), str(args.workdir / f"rpm_{rpm:g}" / f.name))
    finally:
        fst.write_text(fst_backup)
        ed.write_text(ed_backup)
        print("\ndeck restored", flush=True)

    # write the raw frequencies (the classification happens in the test)
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with open(args.csv, "w") as fh:
            fh.write("rpm,freq_hz\n")
            for rpm, freqs in rows:
                for f in freqs:
                    fh.write(f"{rpm},{f:.6f}\n")
        print(f"\nwrote {args.csv}")


if __name__ == "__main__":
    main()
