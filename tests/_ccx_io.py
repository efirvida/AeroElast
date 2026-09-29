"""Shared CalculiX (CCX) helpers for the parity tests.

These three functions used to be private copies inside
``test_composite_beam_parity.py``.  They live here so that new CCX-backed tests
(the composite layup parity and the blade mechanical validation) do not have to
duplicate the fixed-width FRD reader, which is easy to get wrong: CalculiX
writes nodal records as ``(1X,'-1',I10,3E12.5)``, so a negative component abuts
the previous one with no separator and a greedy regex silently drops every node
that has one.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

#: One row of the CalculiX "E I G E N V A L U E   O U T P U T" table:
#: mode number then four numbers (frequency, and three more / eigenvalue / etc.).
_FREQ_ROW = re.compile(
    r"^\s*(\d+)\s+"
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s+"
    r"([+\-]?\d+(?:\.\d+)?(?:[EeDd][+\-]?\d+)?)\s*$"
)


def run_ccx(inp_path: Path, ccx_bin: str) -> subprocess.CompletedProcess:
    """Run CalculiX on ``inp_path`` inside its own directory.

    ``cwd`` is passed to the subprocess instead of mutating the process-wide
    working directory with ``os.chdir`` (which would leak into later tests).
    """
    return subprocess.run(
        [ccx_bin, inp_path.stem],
        cwd=inp_path.parent,
        capture_output=True,
        text=True,
    )


def fail_ccx(result: subprocess.CompletedProcess, inp_path: Path) -> None:
    """Hard-fail on a non-zero CCX exit; ccx writes its diagnostics to stdout."""
    pytest.fail(
        f"CCX failed for {inp_path.name} (rc={result.returncode}):\n"
        f"STDOUT tail:\n{result.stdout[-2000:]}\n"
        f"STDERR tail:\n{result.stderr[-2000:]}"
    )


def parse_frd_disp(frd_file: Path, node_ids: list[int]) -> dict[int, np.ndarray]:
    """Parse the last FRD displacement block, returning only requested nodes.

    Node numbers are read from columns 4-13 and each component from the
    following 12-column fields, so a negative component does not break the
    parse.
    """
    wanted = set(node_ids)
    disps: dict[int, np.ndarray] = {}

    with open(frd_file, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    last_disp_start = -1
    for i, line in enumerate(lines):
        if "-4" in line and "DISP" in line.upper():
            last_disp_start = i

    if last_disp_start == -1:
        return disps

    for line in lines[last_disp_start + 1 :]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("-3") or stripped.startswith("*") or "STEP" in stripped.upper():
            break
        if not stripped.startswith("-1"):
            continue
        try:
            node_id = int(line[3:13])
            if node_id not in wanted:
                continue
            disps[node_id] = np.array(
                [float(line[13:25]), float(line[25:37]), float(line[37:49])],
                dtype=float,
            )
        except ValueError:
            continue

    return disps


def parse_ccx_frequencies(dat_path: Path, n_modes: int = 5) -> np.ndarray:
    """Parse eigenfrequencies (Hz) from a CalculiX ``.dat`` file."""
    dat_file = dat_path.with_suffix(".dat")
    if not dat_file.exists():
        raise RuntimeError(f"DAT file not found: {dat_file}")

    out: dict[int, float] = {}
    in_eig = False
    for line in dat_file.read_text(errors="replace").splitlines():
        up = line.upper()
        if "E I G E N V A L U E" in up and "O U T P U T" in up:
            in_eig = True
            continue
        if not in_eig:
            continue
        m = _FREQ_ROW.match(line)
        if m:
            out[int(m.group(1))] = float(m.group(4).replace("D", "E"))

    if not out:
        raise RuntimeError(f"No frequencies parsed from {dat_file}")
    return np.array([out[k] for k in sorted(out)][:n_modes], dtype=float)
