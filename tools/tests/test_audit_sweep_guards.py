"""Exercise the diagnostic scripts with deterministic, isolated solver doubles."""
from pathlib import Path
import runpy
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(autouse=True)
def isolate_script_import_paths(monkeypatch):
    """Scripts prepend import paths; do not leak those changes into other tests."""
    monkeypatch.setattr(sys, "path", list(sys.path))


def run_g4(monkeypatch, source, displacement=1.0):
    module = ModuleType("test_orthotropic_shell_parity")
    module._build_plate = lambda: None
    module._ORTHO_MAT = object()
    module.thickness = 0.01
    module.B = module.L = 1.0
    module.create_laminate_from_angles = lambda *args: SimpleNamespace(A=np.eye(3))
    def bending():
        module._build_plate()
        module.assert_relative_error(displacement, 1.0)

    module.test_orthotropic_bending = bending
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setattr("inspect.getsource", lambda obj: source)
    return runpy.run_path(str(ROOT / "tools/g4_bending_sweep.py"))


def test_g4_rejects_a_patch_that_no_longer_matches(monkeypatch):
    with pytest.raises(RuntimeError, match="exactly one"):
        run_g4(monkeypatch, "def _build_plate():\n    return None\n")


def test_g4_keeps_closed_form_mismatch_unattributed(monkeypatch, capsys):
    run_g4(monkeypatch, "def _build_plate():\n    nx, nz = 4, 10\n    return nx, nz\n")
    assert "UNATTRIBUTED" in capsys.readouterr().out


def run_g18(monkeypatch, control_gap):
    module = ModuleType("test_composite_layup_parity")
    module.E1 = module.E2 = module.G12 = module.NU12 = module.THICKNESS = module.L = 1.0
    module.B = 0.1
    module.LAYUPS = {"quasi_iso": [45, -45], "sym_0_90s": [0, 90]}
    module.CASES = {"bending": (0, 100.0, 0)}
    module._laminate = lambda angles: angles
    face = SimpleNamespace(nodes={0: SimpleNamespace(id=0)})
    mesh = SimpleNamespace(node_id_to_index={0: 0}, get_node_set=lambda name: face)
    module._build_mesh = lambda nx, ny: mesh

    def solution(mesh, angles, dof, load):
        gap = -1.1 * module.B / 0.1 if angles[0] == 45 else control_gap
        ref = 100.0 / (3 * module.B)
        return np.array([ref * (1 + gap / 100), 0, 0, 0, 0, 0])

    module._aero_solution = solution
    clt = ModuleType("test_composite_beam_parity")
    clt._hand_clt_abd = lambda *args: (np.eye(3), np.zeros((3, 3)), np.eye(3))
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setitem(sys.modules, clt.__name__, clt)
    runpy.run_path(str(ROOT / "tools/g18_quasi_iso_root_sweep.py"))


def test_g18_control_failure_prevents_attribution(monkeypatch, capsys):
    run_g18(monkeypatch, 0.4)
    assert "UNATTRIBUTED" in capsys.readouterr().out


def test_g18_passing_control_preserves_attribution(monkeypatch, capsys):
    run_g18(monkeypatch, 0.1)
    output = capsys.readouterr().out
    assert "=> ATTRIBUTED" in output
    assert "UNATTRIBUTED" not in output
