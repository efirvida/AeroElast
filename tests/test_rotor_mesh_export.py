from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aeroelast.core.config import FSISimulationConfig
from aeroelast.core.mesh import MeshModel, RotorHubMesh, RotorMesh
from aeroelast.core.mesh.utils import detect_open_boundaries, get_open_boundary_loops
from aeroelast.models.blade.numad.objects.definition import Definition


def _make_definition(
    *,
    hub_diameter: float | None = None,
    rotor_diameter: float | None = None,
    blade_length: float | None = None,
) -> Definition:
    definition = Definition()
    definition.hub_diameter = hub_diameter
    definition.rotor_diameter = rotor_diameter
    if blade_length is not None:
        definition.span = np.array([0.0, blade_length], dtype=float)
    return definition


def test_definition_resolve_hub_radius_explicit_override_wins():
    definition = _make_definition(hub_diameter=8.0, rotor_diameter=28.0, blade_length=10.0)

    radius, source = definition.resolve_hub_radius(override=2.5)

    assert radius == pytest.approx(2.5)
    assert source == "explicit"


def test_definition_resolve_hub_radius_uses_hub_diameter_when_available():
    definition = _make_definition(hub_diameter=8.0, rotor_diameter=28.0, blade_length=10.0)

    radius, source = definition.resolve_hub_radius()

    assert radius == pytest.approx(4.0)
    assert source == "hub_diameter"


def test_definition_resolve_hub_radius_falls_back_to_rotor_geometry():
    definition = _make_definition(rotor_diameter=28.0, blade_length=10.0)

    radius, source = definition.resolve_hub_radius()

    assert radius == pytest.approx(4.0)
    assert source == "rotor_diameter_minus_blade_length"


def test_definition_resolve_hub_radius_rejects_negative_results():
    definition = _make_definition(rotor_diameter=12.0, blade_length=10.0)

    with pytest.raises(ValueError, match="negative"):
        definition.resolve_hub_radius()


def test_rotor_hub_mesh_estimates_first_profile_diameter_from_xy_footprint():
    profile = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )

    diameter = RotorHubMesh._estimate_first_profile_diameter(profile)

    assert diameter == pytest.approx(2.0)


def test_rotor_hub_mesh_default_connector_radius_matches_root_profile_scale():
    hub = RotorHubMesh(yaml_file="dummy.yaml")

    _, _, connector_radius, _ = hub._resolve_hub_geometry_parameters(
        hub_radius=4.0,
        root_diameter=5.5,
    )

    assert connector_radius == pytest.approx(2.75)


def test_rotor_hub_mesh_default_connector_length_extends_past_hub_radius():
    hub = RotorHubMesh(yaml_file="dummy.yaml")

    _, connector_length, _, _ = hub._resolve_hub_geometry_parameters(
        hub_radius=4.0,
        root_diameter=5.5,
    )

    assert connector_length == pytest.approx(5.375)


def test_rotor_hub_mesh_default_connector_radius_is_bounded_by_hub_radius():
    hub = RotorHubMesh(yaml_file="dummy.yaml")

    _, _, connector_radius, _ = hub._resolve_hub_geometry_parameters(
        hub_radius=2.0,
        root_diameter=12.0,
    )

    assert connector_radius == pytest.approx(1.9)


def test_fsisimulationconfig_parses_mesh_export_parts(tmp_path: Path):
    config = FSISimulationConfig.from_dict(
        {
            "mesh": {
                "source": "generator",
                "generator": {
                    "type": "RotorMesh",
                    "params": {
                        "yaml_file": "blade.yaml",
                        "n_blades": 3,
                    },
                },
                "export_parts": {
                    "enabled": True,
                    "format": ".STL",
                    "output_dir": "blade_parts",
                },
            },
            "solver": {"type": "LinearStatic"},
            "boundary_conditions": {},
        },
        base_path=tmp_path,
    )

    assert config.mesh.export_parts is not None
    assert config.mesh.export_parts.enabled is True
    assert config.mesh.export_parts.format == "stl"
    assert config.mesh.export_parts.output_dir == str(tmp_path / "blade_parts")


@pytest.fixture(scope="module")
def coarse_rotor_mesh(iea_blade_yaml: str):
    rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=2,
        hub_radius=None,
        element_size=12.0,
        n_samples=60,
    )
    return rotor.generate(renumber=None, verbose=False)


def test_rotor_mesh_auto_hub_radius_applies_expected_offset(iea_blade_yaml: str):
    auto_rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=1,
        hub_radius=None,
        element_size=12.0,
        n_samples=60,
    )
    zero_rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=1,
        hub_radius=0.0,
        element_size=12.0,
        n_samples=60,
    )

    auto_mesh = auto_rotor.generate(renumber=None, verbose=False)
    zero_mesh = zero_rotor.generate(renumber=None, verbose=False)

    resolved_hub_radius, _ = auto_rotor.numad_blade.definition.resolve_hub_radius()

    np.testing.assert_allclose(auto_mesh.coords_array[:, 0], zero_mesh.coords_array[:, 0])
    np.testing.assert_allclose(auto_mesh.coords_array[:, 1], zero_mesh.coords_array[:, 1])
    np.testing.assert_allclose(
        auto_mesh.coords_array[:, 2] - zero_mesh.coords_array[:, 2],
        resolved_hub_radius,
    )


def test_rotor_mesh_creates_aggregate_blade_sets(coarse_rotor_mesh):
    blade_set_names = [RotorMesh.blade_part_set_name(0), RotorMesh.blade_part_set_name(1)]

    for set_name in blade_set_names:
        assert set_name in coarse_rotor_mesh.element_sets_names
        assert set_name in coarse_rotor_mesh.node_sets_names

    blade_meshes = [coarse_rotor_mesh.extract_submesh(name) for name in blade_set_names]
    assert blade_meshes[0].elements_count == blade_meshes[1].elements_count
    assert blade_meshes[0].elements_count * 2 == coarse_rotor_mesh.elements_count


def test_rotor_mesh_write_element_sets_exports_one_file_per_blade(
    tmp_path: Path, coarse_rotor_mesh
):
    blade_set_names = [RotorMesh.blade_part_set_name(0), RotorMesh.blade_part_set_name(1)]

    written_files = coarse_rotor_mesh.write_element_sets(tmp_path, blade_set_names, "stl")

    assert [path.name for path in written_files] == [
        "rotor_blade_1.stl",
        "rotor_blade_2.stl",
    ]
    for path in written_files:
        assert path.exists()
        assert path.stat().st_size > 0


def test_rotor_mesh_write_element_sets_closes_exported_blade_end_loops(
    tmp_path: Path, iea_blade_yaml: str
):
    rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=3,
        hub_radius=None,
        element_size=12.0,
        n_samples=60,
    )
    mesh = rotor.generate(renumber=None, verbose=False)
    blade_set_names = [RotorMesh.blade_part_set_name(i) for i in range(3)]

    written_files = mesh.write_element_sets(tmp_path, blade_set_names, "stl")

    for path in written_files:
        loaded = MeshModel.load(str(path))
        assert get_open_boundary_loops(loaded) == [], path.name


def test_rotor_hub_mesh_generates_closed_surface(iea_blade_yaml: str):
    hub = RotorHubMesh(
        yaml_file=iea_blade_yaml,
        n_blades=3,
        element_size=12.0,
        n_samples=60,
    )

    mesh = hub.generate(renumber=None, verbose=False)

    assert mesh.node_count > 0
    assert mesh.elements_count > 0
    assert RotorHubMesh.SURFACE_SET_NAME in mesh.node_sets_names
    assert RotorHubMesh.SURFACE_SET_NAME in mesh.element_sets_names
    assert detect_open_boundaries(mesh) is False


def test_rotor_mesh_resolves_hub_radius_from_explicit_rotor_diameter_override():
    rotor = RotorMesh(
        excel_file="blade.xlsx",
        airfoil_dir="airfoils",
        n_blades=3,
        rotor_diameter=242.23775645,
    )
    definition = _make_definition(blade_length=117.0)
    rotor._blade_generator = SimpleNamespace(numad_blade=SimpleNamespace(definition=definition))

    radius, source = rotor._resolve_hub_radius()

    assert radius == pytest.approx(4.118878225)
    assert source == "rotor_diameter_minus_blade_length"
