from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path


def _load_sharpy_probe_module():
    repo_root = Path(__file__).resolve().parents[1]
    module_path = repo_root / "src" / "aeroelast" / "solvers" / "aero" / "sharpy_probe.py"
    spec = importlib.util.spec_from_file_location("_test_sharpy_probe", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load SHARPy probe module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


sharpy_probe = _load_sharpy_probe_module()


def test_probe_reports_missing_sharpy_without_raising(monkeypatch):
    real_import_module = sharpy_probe.importlib.import_module

    def fake_import_module(name: str):
        if name == "sharpy" or name.startswith("sharpy."):
            raise ModuleNotFoundError("No module named 'sharpy'")
        return real_import_module(name)

    monkeypatch.setattr(sharpy_probe.importlib, "import_module", fake_import_module)

    report = sharpy_probe.probe_sharpy_installation()

    assert not report.sharpy_available
    assert not report.in_process_uvlm_candidate
    assert any("not importable" in blocker for blocker in report.blockers)
    assert "SHARPy importable: no" in sharpy_probe.format_sharpy_feasibility_report(report)


def test_probe_detects_minimum_in_process_uvlm_surface(monkeypatch):
    module_attrs = {
        "sharpy": {"__version__": "test-sharpy"},
        "sharpy.solvers.stepuvlm": {"StepUvlm": object},
        "sharpy.solvers.staticuvlm": {"StaticUvlm": object},
        "sharpy.solvers.dynamicuvlm": {"DynamicUVLM": object},
        "sharpy.solvers.prescribeduvlm": {"PrescribedUvlm": object},
        "sharpy.solvers.dynamiccoupled": {"DynamicCoupled": object},
        "sharpy.solvers.aerogridloader": {"AerogridLoader": object},
        "sharpy.solvers.beamloader": {"BeamLoader": object},
    }

    def fake_import_module(name: str):
        attrs = module_attrs[name]
        return types.SimpleNamespace(**attrs)

    monkeypatch.setattr(sharpy_probe.importlib, "import_module", fake_import_module)

    report = sharpy_probe.probe_sharpy_installation()

    assert report.sharpy_available
    assert report.sharpy_version == "test-sharpy"
    assert report.in_process_uvlm_candidate
    assert not report.blockers
    assert report.to_dict()["sharpy_version"] == "test-sharpy"
