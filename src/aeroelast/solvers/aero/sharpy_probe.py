"""SHARPy integration feasibility probe.

This module intentionally does not import SHARPy at module import time.  It is a
small, optional diagnostic layer used to decide whether SHARPy can be driven as
an in-process UVLM backend or only as an external benchmark.
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from dataclasses import asdict, dataclass
from typing import Sequence


@dataclass(frozen=True)
class SharpyModuleProbe:
    """Import result for one SHARPy module and expected public attributes."""

    label: str
    module: str
    importable: bool
    attributes_found: tuple[str, ...] = ()
    attributes_missing: tuple[str, ...] = ()
    error: str | None = None


@dataclass(frozen=True)
class SharpyFeasibilityReport:
    """Summary of the SHARPy import/API surface needed by a future backend."""

    python_executable: str
    python_version: str
    sharpy_available: bool
    sharpy_version: str | None
    in_process_uvlm_candidate: bool
    module_probes: tuple[SharpyModuleProbe, ...]
    blockers: tuple[str, ...]
    next_steps: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-serializable representation."""

        return asdict(self)


_MODULE_SPECS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("core package", "sharpy", ()),
    ("single-step UVLM", "sharpy.solvers.stepuvlm", ("StepUvlm",)),
    ("static UVLM", "sharpy.solvers.staticuvlm", ("StaticUvlm",)),
    ("dynamic UVLM", "sharpy.solvers.dynamicuvlm", ("DynamicUVLM",)),
    ("prescribed UVLM", "sharpy.solvers.prescribeduvlm", ("PrescribedUvlm",)),
    ("dynamic coupled loop", "sharpy.solvers.dynamiccoupled", ("DynamicCoupled",)),
    ("aero grid loader", "sharpy.solvers.aerogridloader", ("AerogridLoader",)),
    ("beam loader", "sharpy.solvers.beamloader", ("BeamLoader",)),
)


def _probe_module(label: str, module_name: str, expected_attrs: Sequence[str]) -> SharpyModuleProbe:
    try:
        module = importlib.import_module(module_name)
    except Exception as exc:  # pragma: no cover - exact import errors depend on host env
        return SharpyModuleProbe(
            label=label,
            module=module_name,
            importable=False,
            attributes_missing=tuple(expected_attrs),
            error=f"{type(exc).__name__}: {exc}",
        )

    found = tuple(attr for attr in expected_attrs if hasattr(module, attr))
    missing = tuple(attr for attr in expected_attrs if not hasattr(module, attr))
    return SharpyModuleProbe(
        label=label,
        module=module_name,
        importable=True,
        attributes_found=found,
        attributes_missing=missing,
    )


def probe_sharpy_installation() -> SharpyFeasibilityReport:
    """Inspect whether SHARPy exposes the minimum modules for a UVLM backend."""

    probes = tuple(_probe_module(label, module, attrs) for label, module, attrs in _MODULE_SPECS)
    probe_by_module = {probe.module: probe for probe in probes}
    core_probe = probe_by_module["sharpy"]
    step_probe = probe_by_module["sharpy.solvers.stepuvlm"]
    aero_loader_probe = probe_by_module["sharpy.solvers.aerogridloader"]

    sharpy_version = None
    if core_probe.importable:
        core_module = importlib.import_module("sharpy")
        sharpy_version = str(getattr(core_module, "__version__", "unknown"))

    blockers: list[str] = []
    if not core_probe.importable:
        blockers.append("SHARPy is not importable in the active Python environment.")
    if not step_probe.importable or step_probe.attributes_missing:
        blockers.append("StepUvlm is not available as an importable Python solver class.")
    if not aero_loader_probe.importable or aero_loader_probe.attributes_missing:
        blockers.append("AerogridLoader is not available for in-process aerodynamic grid setup.")

    in_process_candidate = not blockers
    next_steps = [
        "Run a SHARPy NREL-5MW example in an isolated environment before adding dependencies.",
        "Inspect StepUvlm initialisation/run signatures and required Data containers.",
        "Verify whether wake/circulation state can be captured and restored for preCICE rollback.",
    ]
    if in_process_candidate:
        next_steps.append("Build a one-step rigid rotor prototype outside the preCICE participant.")
    else:
        next_steps.append(
            "Keep SHARPy as an external benchmark until the missing API surface is resolved."
        )

    return SharpyFeasibilityReport(
        python_executable=sys.executable,
        python_version=sys.version.split()[0],
        sharpy_available=core_probe.importable,
        sharpy_version=sharpy_version,
        in_process_uvlm_candidate=in_process_candidate,
        module_probes=probes,
        blockers=tuple(blockers),
        next_steps=tuple(next_steps),
    )


def format_sharpy_feasibility_report(report: SharpyFeasibilityReport) -> str:
    """Format a human-readable SHARPy feasibility report."""

    status = "yes" if report.sharpy_available else "no"
    candidate = "yes" if report.in_process_uvlm_candidate else "no"
    lines = [
        "SHARPy Feasibility Probe",
        f"Python: {report.python_version} ({report.python_executable})",
        f"SHARPy importable: {status}",
        f"SHARPy version: {report.sharpy_version or 'n/a'}",
        f"In-process UVLM candidate: {candidate}",
        "",
        "Module surface:",
    ]

    for probe in report.module_probes:
        mark = "ok" if probe.importable and not probe.attributes_missing else "missing"
        detail = probe.module
        if probe.attributes_found:
            detail += f" attrs={','.join(probe.attributes_found)}"
        if probe.attributes_missing:
            detail += f" missing={','.join(probe.attributes_missing)}"
        if probe.error:
            detail += f" error={probe.error}"
        lines.append(f"- {mark}: {probe.label} -> {detail}")

    if report.blockers:
        lines.extend(["", "Blockers:"])
        lines.extend(f"- {blocker}" for blocker in report.blockers)

    lines.extend(["", "Next steps:"])
    lines.extend(f"- {step}" for step in report.next_steps)
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero unless SHARPy looks usable as an in-process UVLM backend.",
    )
    args = parser.parse_args(argv)

    report = probe_sharpy_installation()
    if args.json:
        print(json.dumps(report.to_dict(), indent=2, sort_keys=True))
    else:
        print(format_sharpy_feasibility_report(report))

    if args.strict and not report.in_process_uvlm_candidate:
        return 2
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
