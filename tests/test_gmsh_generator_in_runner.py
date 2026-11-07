"""The runner builds the blade with the gmsh generator when the config asks for ``GmshBladeMesh`` (odd/tasks/gmsh-official-mesher.md, R1)."""

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

from aeroelast.core.config import MeshGeneratorType, MeshGeneratorConfig
from aeroelast.solvers.fsi.runner import FSIRunner


class _SpyGmsh:
    calls: list = []

    def __init__(self, **kwargs):
        type(self).calls.append(kwargs)
        self.numad_mesh_data = {"materials": [], "sections": []}
        self.numad_blade = None

    def generate(self, **_kwargs):
        return SimpleNamespace(tag="gmsh-mesh")


@pytest.fixture
def spy(monkeypatch):
    _SpyGmsh.calls = []
    monkeypatch.setattr("aeroelast.solvers.fsi.runner.GmshBladeMesh", _SpyGmsh)
    return _SpyGmsh


def _runner(params):
    runner = object.__new__(FSIRunner)
    runner._console = Console(file=io.StringIO(), force_terminal=False)
    runner._mesh_generator = None
    runner.config_path = None
    runner.config = SimpleNamespace(
        material=None,
        mesh=SimpleNamespace(source="generator", generator=SimpleNamespace(type="GmshBladeMesh", params=params)),
    )
    return runner


def test_the_generator_type_exists_and_has_typed_params():
    assert MeshGeneratorType.GMSH_BLADE.value == "GmshBladeMesh"
    cfg = MeshGeneratorConfig(type="GmshBladeMesh", params={"yaml_file": "blade.yaml", "element_size": 0.25})
    params = cfg.get_typed_params()
    assert params.yaml_file == "blade.yaml" and params.element_size == 0.25 and params.skin_reference == "mid"


def test_runner_generates_with_the_gmsh_class_and_the_mid_surface_by_default(spy):
    mesh = FSIRunner._generate_mesh(_runner({"yaml_file": "blade.yaml", "element_size": 0.177}))
    assert mesh.tag == "gmsh-mesh"
    assert spy.calls[-1] == {"yaml_file": "blade.yaml", "element_size": 0.177, "skin_reference": "mid"}


def test_runner_honours_the_surface_key_for_the_gmsh_generator(spy):
    FSIRunner._generate_mesh(_runner({"yaml_file": "blade.yaml", "element_size": 0.25, "skin_reference": "oml"}))
    assert spy.calls[-1]["skin_reference"] == "oml"


def test_composite_materials_come_from_the_gmsh_generator(spy, monkeypatch):
    runner = _runner({"yaml_file": "blade.yaml", "element_size": 0.25})
    assert FSIRunner._uses_generator_composite_materials.fget(runner) is True
    FSIRunner._rehydrate_composite_generator(runner)
    assert spy.calls and isinstance(runner._mesh_generator, _SpyGmsh)
