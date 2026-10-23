#!/usr/bin/env python3
"""Fail when a coarse blade mesh enters the repository (AGENTS.md rule 9).

Scans the tracked and untracked (non-ignored) tests, tools, src and scripts for a blade element
size >= THRESHOLD: an assignment to any name containing ``element_size`` (``GATE_ELEMENT_SIZE``),
a positional size passed to a mesh builder (``build_assembly(0.5)``) or a tuple of sizes
(``MESHES = (...)``) and compares it with ``scripts/mesh-size-debt.txt``.

    scripts/check_mesh_size.py          # exit 1 on a NEW coarse mesh
    scripts/check_mesh_size.py --list   # print the registry with the findings

A line of the registry is ``path  v1,v2  category  # note``. A file may carry the listed values
and no other >= THRESHOLD value. The registry is debt, not permission: new entries need the reason
in ``docs/blade_mesh_size_rule.md``. A listed file that no longer has a coarse value is reported
as cleared so the entry can be deleted.

On a node without ``git`` (SLURM compute nodes) the file list comes from the snapshot that
``scripts/check_sbatch.sh`` writes on the login node (``$CHECK_GIT_SNAPSHOT/tracked``); with
neither, the script exits 2 rather than give a verdict. ``CHECK_GIT`` names the git executable.
"""

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

THRESHOLD = 0.3
ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "scripts" / "mesh-size-debt.txt"
SUFFIXES = (".py", ".yaml", ".yml", ".srm", ".sh")
NUM = r"([0-9]*\.?[0-9]+)"
# A size assigned to any name containing element_size/elem_size (GATE_ELEMENT_SIZE, _ELEMENT_SIZE).
PATTERN = re.compile(r"\b\w*(?:element_size|elem_size)\b\s*[=:]\s*" + NUM, re.I)
# A size passed positionally to a blade mesh builder (build_assembly(0.5)).
POSITIONAL = re.compile(
    r"\b(?:_?build_assembly|build_mesh_and_properties|_modal_couple)\(\s*" + NUM)
# A flat tuple or list of sizes (MESHES = (2.0, 1.0, 0.5), ELEMENT_SIZES = [...]); a nested
# sequence of integer division counts (MESHES = [(40, 2), ...]) is not a size and is skipped.
SEQUENCE = re.compile(
    r"\b\w*(?:element_sizes|elem_sizes|meshes|mesh_sizes)\b\s*=\s*[\(\[]([^\(\)\[\]]*)[\)\]]",
    re.I)


def coarse_values(text: str) -> set[float]:
    """Every blade element size >= THRESHOLD written in ``text``, in any of the three spellings."""
    values = [float(m.group(1)) for m in PATTERN.finditer(text)]
    values += [float(m.group(1)) for m in POSITIONAL.finditer(text)]
    for m in SEQUENCE.finditer(text):
        values += [float(v) for v in re.findall(r"[0-9]*\.[0-9]+", m.group(1))]
    return {v for v in values if v >= THRESHOLD}


def tracked_files() -> list[str]:
    git = shutil.which(os.environ.get("CHECK_GIT", "git"))
    snapshot = os.environ.get("CHECK_GIT_SNAPSHOT")
    if git:
        out = subprocess.check_output(
            [git, "ls-files", "--cached", "--others", "--exclude-standard",
             "tests", "tools", "src", "scripts"], cwd=ROOT, text=True)
    elif snapshot and (Path(snapshot) / "tracked").is_file():
        out = (Path(snapshot) / "tracked").read_text()
    else:
        print("mesh-size: ERROR - git is not on PATH and CHECK_GIT_SNAPSHOT names no snapshot "
              "with a 'tracked' file; on a compute node submit with scripts/check_sbatch.sh",
              file=sys.stderr)
        sys.exit(2)
    own = {"scripts/check_mesh_size.py", "tests/software/test_check_mesh_size_patterns.py"}
    return [f for f in out.split() if f.endswith(SUFFIXES) and f not in own]


def scan() -> dict[str, set[float]]:
    found: dict[str, set[float]] = {}
    for rel in tracked_files():
        try:
            text = (ROOT / rel).read_text(errors="replace")
        except OSError:
            continue
        values = coarse_values(text)
        if values:
            found.setdefault(rel, set()).update(values)
    return found


def registry() -> dict[str, tuple[set[float], str]]:
    reg = {}
    for line in REGISTRY.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            path, values, category = line.split()[:3]
            reg[path] = ({float(v) for v in values.split(",")}, category)
    return reg


def main() -> int:
    found, reg = scan(), registry()
    new = {f: v - reg.get(f, (set(), ""))[0] for f, v in found.items() if v - reg.get(f, (set(), ""))[0]}
    cleared = [f for f in reg if f not in found]
    if "--list" in sys.argv:
        for f, (vals, cat) in sorted(reg.items()):
            print(f"{cat:<18} {f}  {sorted(vals)}")
    for f in cleared:
        print(f"mesh-size: cleared, remove from scripts/mesh-size-debt.txt: {f}")
    if new:
        for f, v in sorted(new.items()):
            print(f"mesh-size: NEW coarse mesh {sorted(v)} m in {f}")
        print(f"mesh-size: FAIL - the blade is meshed at 0.15-0.25 m (AGENTS.md rule 9, "
              f"docs/blade_mesh_size_rule.md); {len(new)} file(s) outside the debt registry")
        return 1
    print(f"mesh-size: OK ({len(found)} registered file(s) with element_size >= {THRESHOLD}, no new)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
