#!/usr/bin/env python3
"""Run the SHARPy integration feasibility probe."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_probe_main():
    repo_root = Path(__file__).resolve().parents[1]
    module_path = repo_root / "src" / "aeroelast" / "solvers" / "aero" / "sharpy_probe.py"
    spec = importlib.util.spec_from_file_location("_aeroelast_sharpy_probe", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load SHARPy probe module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.main


if __name__ == "__main__":
    raise SystemExit(_load_probe_main()())
