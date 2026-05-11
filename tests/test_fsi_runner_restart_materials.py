"""Regression tests for FSIRunner material handling on restart."""

from __future__ import annotations

import io
import sys
from enum import Enum
from types import SimpleNamespace

try:
    import _aeroelast  # noqa: F401
except (ImportError, OSError):

    class _DummyElementFamily(Enum):
        PLANE = "plane"
        SHELL = "shell"
        SOLID = "solid"

    sys.modules["_aeroelast"] = SimpleNamespace(ElementFamily=_DummyElementFamily)

from rich.console import Console

from aeroelast.core.config import MeshGeneratorType, MeshSource
from aeroelast.solvers.fsi.runner import FSIRunner


def test_create_material_uses_blade_properties_on_checkpoint_restart(monkeypatch):
    runner = object.__new__(FSIRunner)
    runner._console = Console(file=io.StringIO(), force_terminal=False)
    runner._blade_properties = None
    runner._mesh_generator = None
    runner.config = SimpleNamespace(
        material=None,
        mesh=SimpleNamespace(
            source=MeshSource.GENERATOR.value,
            generator=SimpleNamespace(type=MeshGeneratorType.BLADE.value, params={}),
        ),
    )

    expected_properties = {"section": object()}
    calls = []

    def _fake_extract(self):
        calls.append(True)
        return expected_properties

    monkeypatch.setattr(FSIRunner, "_extract_blade_properties", _fake_extract)

    material = FSIRunner._create_material(runner)

    assert material is None
    assert calls == [True]
    assert runner._blade_properties is expected_properties
