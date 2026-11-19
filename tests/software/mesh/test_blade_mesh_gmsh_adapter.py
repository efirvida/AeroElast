"""BladeMesh / Blade build with GmshBladeMesh (odd/tasks/gmsh-official-mesher.md, T6a, S7 (a)-(e)).

One test per S7 rule plus the equivalence of ``BladeMesh`` with ``GmshBladeMesh`` on the official windIO v2 IEA-15 yaml and the
check that the runner's own GMSH_BLADE dispatch is untouched.  The NuMAD meshing modules stay on disk (T9 deletes them); BladeMesh
simply no longer calls them.
"""

from __future__ import annotations

import io
from pathlib import Path
from types import SimpleNamespace

import pytest

from aeroelast.core.config import MeshGeneratorType
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.core.mesh.gmsh_blade import GmshBladeMesh

ROOT = Path(__file__).resolve().parents[3]
V2_YAML = ROOT / "tests" / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"
V1_YAML = ROOT / "tests" / "IEA-15-240-RWT.yaml"
V2_NAME = "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"


def test_s7a_windio_v1_yaml_is_rejected_naming_input_and_v2_file():
    with pytest.raises(ValueError) as exc:
        BladeMesh(yaml_file=str(V1_YAML)).generate()
    msg = str(exc.value)
    assert str(V1_YAML) in msg and V2_NAME in msg and "windIO" in msg


def test_s7a_excel_file_is_rejected_naming_input_and_v2_file():
    with pytest.raises(ValueError) as exc:
        BladeMesh(excel_file="blade.xlsx").generate()
    msg = str(exc.value)
    assert "blade.xlsx" in msg and V2_NAME in msg and "excel_file" in msg


@pytest.mark.parametrize(
    "key, value",
    [("n_samples", 300), ("airfoil_spacing", "cosine"), ("span_grading", "chord"), ("airfoil_dir", "airfoils"),
     ("refine_tip", True), ("max_bay_length", 1.0), ("laminate_per_row", True)],
)
def test_s7b_each_numad_only_key_is_rejected_by_name_even_at_its_old_default(key, value):
    with pytest.raises(ValueError) as exc:
        BladeMesh(yaml_file=str(V2_YAML), **{key: value}).generate()
    assert key in str(exc.value) and "NuMAD" in str(exc.value)


def test_s7b_all_offending_keys_are_listed_together():
    with pytest.raises(ValueError) as exc:
        BladeMesh(yaml_file=str(V2_YAML), n_samples=300, span_grading="chord", refine_tip=False).generate()
    msg = str(exc.value)
    assert all(k in msg for k in ("n_samples", "span_grading", "refine_tip"))


def test_s7c_no_size_means_the_reference_mesh_and_never_a_coarse_default(monkeypatch):
    seen = {}

    class _Spy:
        def __init__(self, **kwargs):
            seen.update(kwargs)
            self.numad_mesh_data, self.numad_blade = {"k": 1}, object()

        def generate(self, renumber=None, verbose=False):
            return "mesh"

    monkeypatch.setattr("aeroelast.core.mesh.gmsh_blade.GmshBladeMesh", _Spy)
    mesh = BladeMesh(yaml_file=str(V2_YAML)).generate()
    assert mesh == "mesh"
    assert seen["chord_divisions"] == 20 and seen["span_aspect"] == 1
    # an explicit element size keeps the uniform path of GmshBladeMesh
    seen.clear()
    BladeMesh(yaml_file=str(V2_YAML), element_size=0.2).generate()
    assert seen["element_size"] == 0.2 and seen["chord_divisions"] is None
    import inspect

    default = inspect.signature(BladeMesh.__init__).parameters["element_size"].default
    assert default is None


def test_s7d_adapter_exposes_the_gmsh_contract(monkeypatch):
    blade, data = object(), {"materials": {}}

    class _Spy:
        def __init__(self, **kwargs):
            self.numad_mesh_data, self.numad_blade = data, blade

        def generate(self, renumber=None, verbose=False):
            return "mesh"

    monkeypatch.setattr("aeroelast.core.mesh.gmsh_blade.GmshBladeMesh", _Spy)
    gen = BladeMesh(yaml_file=str(V2_YAML), chord_divisions=20)
    with pytest.raises(RuntimeError):
        gen.numad_mesh_data  # noqa: B018
    gen.generate()
    assert gen.numad_mesh_data is data and gen.numad_blade is blade


def test_s7e_no_webs_is_t7():
    with pytest.raises(NotImplementedError, match="T7"):
        BladeMesh(yaml_file=str(V2_YAML), include_webs=False).generate()


def test_s7e_geometry_only_is_t7():
    with pytest.raises(NotImplementedError, match="T7"):
        BladeMesh(yaml_file=str(V2_YAML)).generate(geometry_only=True)


def test_blade_wrapper_forwards_no_numad_key_and_has_no_blocks_sections(monkeypatch):
    from aeroelast.models.blade.model import Blade

    calls = []

    class _Spy:
        def __init__(self, **kwargs):
            calls.append(kwargs)
            self.numad_mesh_data, self.numad_blade = {}, None

        def generate(self, **_kw):
            return "mesh"

    monkeypatch.setattr("aeroelast.models.blade.model.BladeMesh", _Spy)
    blade = Blade(str(V2_YAML))
    blade.generate_mesh()
    assert "n_samples" not in calls[-1] and calls[-1]["element_size"] is None
    assert not hasattr(blade, "generate_blocks_sections")  # NuMAD block meshing deleted (T11, S9)


def test_runner_gmsh_blade_dispatch_is_unchanged(monkeypatch):
    from rich.console import Console

    from aeroelast.solvers.fsi import runner as runner_mod

    seen = []

    class _SpyGmsh:
        def __init__(self, **kwargs):
            seen.append(kwargs)

        def generate(self, renumber=None, verbose=False):
            return SimpleNamespace(tag="mesh")

    def _no_numad(**_kw):
        raise AssertionError("the GMSH_BLADE branch must not build a NuMAD BladeMesh")

    monkeypatch.setattr(runner_mod, "GmshBladeMesh", _SpyGmsh)
    monkeypatch.setattr(runner_mod, "BladeMesh", _no_numad)
    params = {"yaml_file": "blade.yaml", "element_size": 0.2, "chord_divisions": 20, "span_aspect": 1.0}
    runner = object.__new__(runner_mod.FSIRunner)
    runner._console = Console(file=io.StringIO(), force_terminal=False)
    runner.config_path = None
    runner.config = SimpleNamespace(
        mesh=SimpleNamespace(generator=SimpleNamespace(type=MeshGeneratorType.GMSH_BLADE.value, params=params)))
    runner._generate_mesh()
    assert seen == [{"yaml_file": "blade.yaml", "element_size": 0.2, "skin_reference": "mid",
                     "chord_divisions": 20, "span_aspect": 1.0}]


@pytest.mark.slow
def test_adapter_matches_gmsh_blade_mesh_counts_and_set_names():
    """Reference mesh (chord_divisions 20, span_aspect 1), both built for real from the same v2 yaml (~100 s CPU each)."""
    direct = GmshBladeMesh(str(V2_YAML), skin_reference="oml", chord_divisions=20, span_aspect=1).generate()
    adapter = BladeMesh(yaml_file=str(V2_YAML), skin_reference="oml", chord_divisions=20, span_aspect=1)
    mesh = adapter.generate()
    assert len(mesh.nodes) == len(direct.nodes) and mesh.elements_count == direct.elements_count
    assert set(mesh.element_sets) == set(direct.element_sets) and set(mesh.node_sets) == set(direct.node_sets)
    assert {"RootNodes", "TipNodes", "allOuterShellNods"} <= set(mesh.node_sets)
    assert set(adapter.numad_mesh_data) == set(direct.numad_mesh_data)


# --- T6c: callers forward only what the config / caller gave (S7b, S7c, AGENTS rule 9) -------------------------------------


def _runner_with(gen_type, params):
    from rich.console import Console

    from aeroelast.solvers.fsi import runner as runner_mod

    runner = object.__new__(runner_mod.FSIRunner)
    runner._console = Console(file=io.StringIO(), force_terminal=False)
    runner.config_path = None
    runner.config = SimpleNamespace(mesh=SimpleNamespace(generator=SimpleNamespace(type=gen_type, params=params)))
    return runner_mod, runner


def _spy(seen):
    class _Spy:
        def __init__(self, **kwargs):
            seen.append(kwargs)

        def generate(self, renumber=None, verbose=False):
            return SimpleNamespace(tag="mesh")

    return _Spy


def test_t6c_runner_blade_branch_forwards_only_the_keys_the_config_gives(monkeypatch):
    seen = []
    runner_mod, runner = _runner_with(MeshGeneratorType.BLADE.value, {"yaml_file": "blade.yaml", "element_size": 0.2})
    monkeypatch.setattr(runner_mod, "BladeMesh", _spy(seen))
    runner._generate_mesh()
    assert seen == [{"yaml_file": "blade.yaml", "skin_reference": "mid", "element_size": 0.2}]


def test_t6c_runner_blade_branch_without_size_injects_no_size_and_no_numad_default(monkeypatch):
    seen = []
    runner_mod, runner = _runner_with(MeshGeneratorType.BLADE.value, {"yaml_file": "blade.yaml"})
    monkeypatch.setattr(runner_mod, "BladeMesh", _spy(seen))
    runner._generate_mesh()
    assert seen == [{"yaml_file": "blade.yaml", "skin_reference": "mid"}]


def test_t6c_runner_blade_branch_with_numad_keys_gets_the_s7b_error():
    runner_mod, runner = _runner_with(
        MeshGeneratorType.BLADE.value,
        {"yaml_file": str(V2_YAML), "n_samples": 300, "span_grading": "chord"},
    )
    with pytest.raises(ValueError) as exc:
        runner._generate_mesh()
    assert "n_samples" in str(exc.value) and "span_grading" in str(exc.value) and "NuMAD" in str(exc.value)


def test_t6c_runner_rotor_branch_has_no_coarse_default_and_forwards_only_given_keys(monkeypatch):
    seen = []
    runner_mod, runner = _runner_with(MeshGeneratorType.ROTOR.value, {"yaml_file": "blade.yaml", "n_blades": 3})
    monkeypatch.setattr(runner_mod, "RotorMesh", _spy(seen))
    runner._generate_mesh()
    assert len(seen) == 1
    kw = seen[0]
    assert "element_size" not in kw and "n_samples" not in kw and "airfoil_spacing" not in kw
    assert kw["yaml_file"] == "blade.yaml" and kw["n_blades"] == 3


def test_t6c_rotormesh_forwards_only_given_keys_to_blademesh_and_has_no_coarse_default(monkeypatch):
    from aeroelast.core.mesh import generators
    from aeroelast.core.mesh.generators import RotorMesh

    seen = []

    class _Spy:
        numad_blade = SimpleNamespace(definition=SimpleNamespace(resolve_hub_radius=lambda override=None: (0.0, "x")))

        def __init__(self, **kwargs):
            seen.append(kwargs)

        def generate(self, renumber=None, verbose=False):
            from aeroelast.core.mesh.model import MeshModel

            return MeshModel()

    monkeypatch.setattr(generators, "BladeMesh", _Spy)
    RotorMesh(yaml_file="blade.yaml", n_blades=1).generate(verbose=False)
    assert seen == [{"yaml_file": "blade.yaml"}]
    seen.clear()
    RotorMesh(yaml_file="blade.yaml", n_blades=1, element_size=0.2, chord_divisions=12).generate(verbose=False)
    assert seen == [{"yaml_file": "blade.yaml", "element_size": 0.2, "chord_divisions": 12}]
    seen.clear()
    RotorMesh(yaml_file="blade.yaml", n_blades=1, n_samples=60).generate(verbose=False)
    assert seen == [{"yaml_file": "blade.yaml", "n_samples": 60}]
    assert inspect_default(RotorMesh, "element_size") is None


def inspect_default(cls, name):
    import inspect

    return inspect.signature(cls.__init__).parameters[name].default


def test_t6c_turbinemesh_forwards_only_given_keys_to_blademesh(monkeypatch):
    from aeroelast.core.mesh import turbine

    seen = []
    monkeypatch.setattr(turbine, "BladeMesh", _spy(seen))
    turbine.TurbineMesh("t.yaml", include_webs=True)._build_blade(verbose=False)
    assert seen == [{"yaml_file": "t.yaml", "include_webs": True}]
    assert inspect_default(turbine.TurbineMesh, "element_size") is None


def test_t6c_cli_mesh_blade_forwards_only_given_options(monkeypatch):
    from aeroelast.cli import mesh as cli_mesh
    from aeroelast.core.mesh import generators

    seen = []
    monkeypatch.setattr(generators, "BladeMesh", _spy(seen))
    monkeypatch.setattr(cli_mesh, "_write_mesh", lambda *a, **k: None)
    args = cli_mesh._build_parser().parse_args(["blade", "t.yaml", "--out", "b.stl", "--quiet"])
    args.func(args)
    assert seen == [{"yaml_file": "t.yaml"}]


def test_rotor_and_blade_model_defaults_forward_no_numad_key(monkeypatch):
    """Rotor and BladeModel defaults must not inject NuMAD-only keys (S7b) or a coarse size (S7c)."""
    import aeroelast.models.blade.model as model

    captured = {}

    class _FakeRotorMesh:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def generate(self, **_):
            return None

    monkeypatch.setattr(model, "RotorMesh", _FakeRotorMesh)
    model.Rotor("blade.yaml", 3).generate_mesh()
    assert "n_samples" not in captured
    assert captured.get("element_size") is None
    import inspect
    params = inspect.signature(model.BladeModel.__init__).parameters
    assert "n_samples" not in params  # NuMAD sampling key deleted (T11, S9)
    assert all(p.default is None for p in params.values() if p.default is not inspect.Parameter.empty)


def test_mesh_param_dataclasses_carry_no_numad_default():
    from aeroelast.core.config import BladeMeshParams, RotorMeshParams

    for cls in (BladeMeshParams, RotorMeshParams):
        params = cls(yaml_file="blade.yaml")
        assert not hasattr(params, "n_samples") and not hasattr(params, "airfoil_spacing")
        assert params.element_size is None
