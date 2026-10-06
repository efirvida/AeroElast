"""Wiring guard: every rotor FSI case on disk must reach the PRODUCTION solver.

The production rotor solver is ``LinearDynamicFSIRotorCorotationalSolver``
(``solvers/fsi/rotor.py``). The legacy enum value ``LinearDynamicFSIRotor`` is an
alias for it (``rotor.py`` ends with ``LinearDynamicFSIRotorSolver =
LinearDynamicFSIRotorCorotationalSolver``), so both names dispatch to the same
class through ``runner.py::_create_solver``.

``LinearDynamicFSIRotorInertialSolver`` is in development and deliberately out of
scope: a case that silently switches to it would be validating a half-built
solver. Cases whose filename names the inertial variant may declare it; every
other case must declare the production solver.

This guard is fast: no mesh, no preCICE, no solver instantiation.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
TESTS_DIR = REPO_ROOT / "tests"

ALIAS = "LinearDynamicFSIRotor"
EXPLICIT = "LinearDynamicFSIRotorCorotational"
INERTIAL = "LinearDynamicFSIRotorInertial"

PRODUCTION_NAMES = {ALIAS, EXPLICIT}

_SOLVER_BLOCK = re.compile(r"^solver:[ \t]*$", re.MULTILINE)
_TYPE_LINE = re.compile(r"^[ \t]+type:[ \t]*(\S+)[ \t]*$", re.MULTILINE)


def _declared_solver_type(text: str) -> str | None:
    """Return the ``solver.type`` value of a YAML document, or None if absent."""
    block = _SOLVER_BLOCK.search(text)
    if block is None:
        return None
    match = _TYPE_LINE.search(text, block.end())
    if match is None:
        return None
    # The first indented `type:` after `solver:` must belong to that block: stop
    # if a new top-level key starts first.
    between = text[block.end() : match.start()]
    if any(line and not line[0].isspace() for line in between.splitlines()):
        return None
    return match.group(1)


def _case_declarations() -> list[tuple[Path, str]]:
    """(path, solver type) for every YAML under tests/ that names a rotor solver."""
    found: list[tuple[Path, str]] = []
    for path in sorted(TESTS_DIR.rglob("*.yaml")):
        try:
            text = path.read_text(errors="replace")
        except OSError:  # pragma: no cover - unreadable file is not this test's business
            continue
        if "LinearDynamicFSIRotor" not in text:
            continue
        declared = _declared_solver_type(text)
        if declared is not None:
            found.append((path, declared))
    return found


def test_alias_is_the_production_corotational_class():
    from aeroelast.solvers.fsi.rotor import (
        LinearDynamicFSIRotorCorotationalSolver,
        LinearDynamicFSIRotorSolver,
    )

    assert LinearDynamicFSIRotorSolver is LinearDynamicFSIRotorCorotationalSolver, (
        "the legacy alias no longer points at the corotational solver: every case "
        "declaring it would silently change physics"
    )


def test_enum_values_match_the_case_yaml_names():
    from aeroelast.core.config import SolverType

    assert SolverType.LINEAR_DYNAMIC_FSI_ROTOR.value == ALIAS
    assert SolverType.LINEAR_DYNAMIC_FSI_ROTOR_COROTATIONAL.value == EXPLICIT
    assert SolverType.LINEAR_DYNAMIC_FSI_ROTOR_INERTIAL.value == INERTIAL


def test_every_case_declares_the_production_solver():
    declarations = _case_declarations()
    assert declarations, "no rotor FSI case found under tests/ - the scan is broken"

    off_path = []
    inertial_cases = []
    for path, declared in declarations:
        rel = path.relative_to(REPO_ROOT).as_posix()
        if declared in PRODUCTION_NAMES:
            continue
        if declared == INERTIAL and "inertial" in path.name:
            inertial_cases.append(rel)
            continue
        off_path.append(f"{rel} -> {declared}")

    assert not off_path, (
        "these rotor cases declare a solver outside production "
        f"({sorted(PRODUCTION_NAMES)}) without naming it in the filename: " + "; ".join(off_path)
    )
    # The guard is only meaningful if it actually saw the corpus.
    assert len(declarations) >= 60, (
        f"expected the full case corpus, scanned only {len(declarations)} declarations"
    )
