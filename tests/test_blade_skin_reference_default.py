"""The production blade path builds its shell on the mid-surface unless told otherwise.

Every structural anchor of L3 (S-1, S-2, S-4, S-7 and the gmsh campaign) is measured on ``skin_reference="mid"``, while ``BladeMesh``
defaults to ``"oml"``.  At rated load the same blade deflects 16.2 m (oml) against 18.7 m (mid) at the tip
(odd/tasks/gmsh-official-mesher.md L5), so the FSI must not silently run another configuration than the one validated.
``BladeMesh`` itself keeps its default; the production entry points (the ``Blade`` wrapper and the two ``BladeMesh`` sites of
``FSIRunner``) pass ``skin_reference`` explicitly: ``mid`` by default, a generator key ``skin_reference`` overrides it.
"""

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

import pytest
from rich.console import Console

from aeroelast.core.config import MeshGeneratorType
from aeroelast.solvers.fsi.runner import FSIRunner


class _SpyBladeMesh:
    """Stands in for BladeMesh: records the constructor kwargs, builds nothing."""

    calls: list = []

    def __init__(self, **kwargs):
        type(self).calls.append(kwargs)
        self.numad_mesh_data = {}
        self.numad_blade = None

    def generate(self, **_kwargs):
        return SimpleNamespace(tag="mesh")


@pytest.fixture
def spy(monkeypatch):
    _SpyBladeMesh.calls = []
    monkeypatch.setattr("aeroelast.solvers.fsi.runner.BladeMesh", _SpyBladeMesh)
    monkeypatch.setattr("aeroelast.models.blade.model.BladeMesh", _SpyBladeMesh)
    return _SpyBladeMesh


def _runner(params):
    runner = object.__new__(FSIRunner)
    runner._console = Console(file=io.StringIO(), force_terminal=False)
    runner._mesh_generator = None
    runner.config_path = None
    runner.config = SimpleNamespace(
        mesh=SimpleNamespace(generator=SimpleNamespace(type=MeshGeneratorType.BLADE.value, params=params))
    )
    return runner


def test_blade_wrapper_uses_the_mid_surface_by_default(spy):
    from aeroelast.models.blade.model import Blade

    Blade("blade.yaml", element_size=0.25).generate_mesh()
    assert spy.calls[-1]["skin_reference"] == "mid"


def test_blade_wrapper_accepts_an_explicit_oml(spy):
    from aeroelast.models.blade.model import Blade

    Blade("blade.yaml", element_size=0.25, skin_reference="oml").generate_mesh()
    assert spy.calls[-1]["skin_reference"] == "oml"


def test_runner_generate_mesh_uses_the_mid_surface_by_default(spy):
    runner = _runner({"yaml_file": "blade.yaml", "element_size": 0.25})
    FSIRunner._generate_mesh(runner)
    assert spy.calls[-1]["skin_reference"] == "mid"


def test_runner_generate_mesh_honours_the_generator_key(spy):
    runner = _runner({"yaml_file": "blade.yaml", "element_size": 0.25, "skin_reference": "oml"})
    FSIRunner._generate_mesh(runner)
    assert spy.calls[-1]["skin_reference"] == "oml"


def test_runner_rehydration_uses_the_same_surface(spy, monkeypatch):
    runner = _runner({"yaml_file": "blade.yaml", "element_size": 0.25})
    monkeypatch.setattr(FSIRunner, "_uses_generator_composite_materials", property(lambda self: True))
    FSIRunner._rehydrate_composite_generator(runner)
    assert spy.calls[-1]["skin_reference"] == "mid"


def test_an_invalid_surface_is_rejected_before_any_mesh_is_built(spy):
    runner = _runner({"yaml_file": "blade.yaml", "skin_reference": "inner"})
    with pytest.raises(ValueError, match="skin_reference"):
        FSIRunner._generate_mesh(runner)
    assert spy.calls == []
