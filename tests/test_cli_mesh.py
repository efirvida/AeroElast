"""Tests for the ``aeroelast mesh`` CLI."""

import os

import pytest

from aeroelast.cli.aeroelast import main as cli_main
from aeroelast.cli.mesh import main as mesh_main

_HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(_HERE, "IEA-15-240-RWT.yaml")

_SURFACE = ["--n-circ", "8"]


def test_help_exits_zero():
    with pytest.raises(SystemExit) as exc:
        mesh_main(["--help"])
    assert exc.value.code == 0


def test_missing_subcommand_errors():
    with pytest.raises(SystemExit) as exc:
        mesh_main([])
    assert exc.value.code == 2


def test_tower_from_params_writes_obj(tmp_path):
    out = tmp_path / "tower.obj"
    rc = mesh_main(
        [
            "tower",
            "--height",
            "50",
            "--base-diameter",
            "6",
            "--top-diameter",
            "3",
            "--n-axial",
            "3",
            *_SURFACE,
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.stat().st_size > 0
    assert "v " in out.read_text()


def test_tower_from_windio(tmp_path):
    out = tmp_path / "tower.stl"
    rc = mesh_main(["tower", _YAML, "--n-axial", "4", *_SURFACE, "--quiet", "--out", str(out)])
    assert rc == 0
    assert out.stat().st_size > 0


def test_hub_from_windio(tmp_path):
    out = tmp_path / "hub.stl"
    rc = mesh_main(["hub", _YAML, *_SURFACE, "--n-merid", "6", "--quiet", "--out", str(out)])
    assert rc == 0
    assert out.stat().st_size > 0


def test_nacelle_from_windio(tmp_path):
    out = tmp_path / "nacelle.vtk"
    rc = mesh_main(
        [
            "nacelle",
            _YAML,
            *_SURFACE,
            "--n-axial",
            "4",
            "--n-nose",
            "3",
            "--quiet",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.stat().st_size > 0


def test_hub_without_input_or_diameter_fails(tmp_path):
    rc = mesh_main(["hub", "--out", str(tmp_path / "hub.stl")])
    assert rc == 1


def test_tower_height_without_base_diameter_fails(tmp_path):
    rc = mesh_main(["tower", "--height", "10", "--out", str(tmp_path / "tower.stl")])
    assert rc == 1


def test_output_extension_fallback(tmp_path):
    stem = tmp_path / "tower"
    rc = mesh_main(
        [
            "tower",
            "--height",
            "30",
            "--base-diameter",
            "4",
            "--n-axial",
            "2",
            *_SURFACE,
            "--out",
            str(stem),
            "--format",
            "vtk",
        ]
    )
    assert rc == 0
    assert (tmp_path / "tower.vtk").stat().st_size > 0


def test_blade_no_webs(tmp_path):
    out = tmp_path / "blade.stl"
    rc = mesh_main(
        [
            "blade",
            _YAML,
            "--element-size",
            "15",
            "--no-webs",
            "--quiet",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.stat().st_size > 0


def test_rotor(tmp_path):
    out = tmp_path / "rotor.vtk"
    rc = mesh_main(
        [
            "rotor",
            _YAML,
            "--n-blades",
            "2",
            "--element-size",
            "20",
            "--no-webs",
            "--quiet",
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.stat().st_size > 0


def test_turbine_writes_every_component(tmp_path):
    out_dir = tmp_path / "meshes"
    rc = mesh_main(
        [
            "turbine",
            _YAML,
            "--out-dir",
            str(out_dir),
            "--format",
            "stl",
            "--n-blades",
            "1",
            "--element-size",
            "20",
            "--tower-n-axial",
            "4",
            "--tower-n-circ",
            "8",
            "--hub-n-circ",
            "8",
            "--hub-n-merid",
            "6",
            "--nacelle-n-circ",
            "8",
            "--quiet",
        ]
    )
    assert rc == 0
    assert sorted(p.name for p in out_dir.iterdir()) == [
        "blade_1.stl",
        "hub.stl",
        "nacelle.stl",
        "tower.stl",
    ]


def test_dispatch_through_aeroelast_main(tmp_path):
    out = tmp_path / "tower.obj"
    rc = cli_main(
        [
            "mesh",
            "tower",
            "--height",
            "25",
            "--base-diameter",
            "3",
            "--n-axial",
            "2",
            *_SURFACE,
            "--out",
            str(out),
        ]
    )
    assert rc == 0
    assert out.stat().st_size > 0
