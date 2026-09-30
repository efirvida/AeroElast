from pathlib import Path

import pytest

from aeroelast.cli.aeroelast import main


def test_aeroelast_cli_exports_blade_mesh_without_solver_yaml(tmp_path: Path, iea_blade_yaml: str):
    output_path = tmp_path / "blade.vtk"

    exit_code = main([
        "--export-mesh",
        str(output_path),
        "--mesh-generator",
        "BladeMesh",
        "--blade-yaml",
        iea_blade_yaml,
        "--element-size",
        "12.0",
        "--n-samples",
        "60",
    ])

    assert exit_code == 0
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_aeroelast_cli_exports_rotor_mesh_without_solver_yaml(tmp_path: Path, iea_blade_yaml: str):
    output_path = tmp_path / "rotor.vtk"

    exit_code = main([
        "--export-mesh",
        str(output_path),
        "--mesh-generator",
        "RotorMesh",
        "--blade-yaml",
        iea_blade_yaml,
        "--n-blades",
        "2",
        "--element-size",
        "12.0",
        "--n-samples",
        "60",
    ])

    assert exit_code == 0
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_aeroelast_cli_exports_rotor_hub_mesh_without_solver_yaml(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    import aeroelast.core.mesh as mesh_module

    output_path = tmp_path / "hub.stl"
    calls: dict[str, object] = {}

    class FakeMesh:
        def write_mesh(self, filename, **kwargs):
            Path(filename).write_text("solid fake\n")

    class FakeRotorHubMesh:
        def __init__(self, **kwargs):
            calls["kwargs"] = kwargs

        def generate(self, renumber=None, verbose=True):
            calls["renumber"] = renumber
            calls["verbose"] = verbose
            return FakeMesh()

    monkeypatch.setattr(mesh_module, "RotorHubMesh", FakeRotorHubMesh)

    exit_code = main([
        "--export-mesh",
        str(output_path),
        "--mesh-generator",
        "RotorHubMesh",
        "--blade-yaml",
        "tests/IEA-15-240-RWT.yaml",
        "--n-blades",
        "3",
        "--hub-length",
        "8.0",
        "--hub-connector-radius",
        "1.4",
        "--hub-nose-radius",
        "4.2",
    ])

    assert exit_code == 0
    assert output_path.exists()
    assert calls["kwargs"] == {
        "yaml_file": str(Path.cwd() / "tests/IEA-15-240-RWT.yaml"),
        "excel_file": None,
        "airfoil_dir": None,
        "n_blades": 3,
        "hub_radius": None,
        "hub_diameter": None,
        "rotor_diameter": None,
        "element_size": 0.1,
        "n_samples": 300,
        "airfoil_spacing": "cosine",
        "hub_length": 8.0,
        "connector_radius": 1.4,
        "nose_radius": 4.2,
    }
    assert calls["renumber"] is None
    assert calls["verbose"] is False


def test_aeroelast_cli_export_mesh_requires_generator(tmp_path: Path, iea_blade_yaml: str):
    output_path = tmp_path / "mesh.vtk"

    exit_code = main([
        "--export-mesh",
        str(output_path),
        "--blade-yaml",
        iea_blade_yaml,
    ])

    assert exit_code == 1
    assert not output_path.exists()


def test_aeroelast_cli_exports_rotor_blades_separately(tmp_path: Path, iea_blade_yaml: str):
    parts_dir = tmp_path / "parts"

    exit_code = main([
        "--export-parts-dir",
        str(parts_dir),
        "--mesh-generator",
        "RotorMesh",
        "--blade-yaml",
        iea_blade_yaml,
        "--n-blades",
        "3",
        "--element-size",
        "12.0",
        "--n-samples",
        "60",
    ])

    assert exit_code == 0
    assert (parts_dir / "rotor_blade_1.stl").exists()
    assert (parts_dir / "rotor_blade_2.stl").exists()
    assert (parts_dir / "rotor_blade_3.stl").exists()


def test_aeroelast_cli_exports_rotor_blades_from_excel_separately(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    import aeroelast.core.mesh as mesh_module

    parts_dir = tmp_path / "parts"
    excel_file = tmp_path / "blade.xlsx"
    airfoil_dir = tmp_path / "airfoils"
    excel_file.write_text("placeholder")
    airfoil_dir.mkdir()

    calls: dict[str, object] = {}

    class FakeMesh:
        element_sets_names = ["rotor_blade_1", "rotor_blade_2", "rotor_blade_3"]

        def write_element_sets(self, output_dir, set_names, file_format):
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            calls["set_names"] = list(set_names)
            calls["file_format"] = file_format
            written = []
            for set_name in set_names:
                output_file = output_dir / f"{set_name}.{file_format}"
                output_file.write_text("solid fake\n")
                written.append(output_file)
            return written

    class FakeRotorMesh:
        BLADE_PART_SET_PREFIX = "rotor_blade_"

        def __init__(self, **kwargs):
            calls["kwargs"] = kwargs

        def generate(self, renumber=None, verbose=True):
            calls["renumber"] = renumber
            calls["verbose"] = verbose
            return FakeMesh()

    monkeypatch.setattr(mesh_module, "RotorMesh", FakeRotorMesh)

    exit_code = main([
        "--export-parts-dir",
        str(parts_dir),
        "--mesh-generator",
        "RotorMesh",
        "--excel-file",
        str(excel_file),
        "--airfoil-dir",
        str(airfoil_dir),
        "--rotor-diameter",
        "242.23775645",
        "--n-blades",
        "3",
    ])

    assert exit_code == 0
    assert calls["kwargs"] == {
        "yaml_file": None,
        "excel_file": str(excel_file),
        "airfoil_dir": str(airfoil_dir),
        "n_blades": 3,
        "hub_radius": None,
        "hub_diameter": None,
        "rotor_diameter": 242.23775645,
        "element_size": 0.1,
        "n_samples": 300,
        "airfoil_spacing": "cosine",
    }
    assert calls["renumber"] is None
    assert calls["verbose"] is False
    assert calls["set_names"] == ["rotor_blade_1", "rotor_blade_2", "rotor_blade_3"]
    assert calls["file_format"] == "stl"
    assert (parts_dir / "rotor_blade_1.stl").exists()
    assert (parts_dir / "rotor_blade_2.stl").exists()
    assert (parts_dir / "rotor_blade_3.stl").exists()
