"""Turbine mesh assembly: separate blade/hub/nacelle/tower meshes in one frame."""

import os

import numpy as np
import pytest

from aeroelast.core.mesh.turbine import TurbineMesh

_HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(_HERE, "IEA-15-240-RWT.yaml")

_HUB_CENTER = np.array([12.0313, 0.0, 135.0])
_HUB_RADIUS = 7.94 / 2.0
_ROTOR_RADIUS = 242.23775645 / 2.0
_UPTILT = 0.10471975511965977
_TOWER_TOP_Z = 129.386


@pytest.fixture(scope="module")
def turbine():
    return TurbineMesh(_YAML, element_size=3.0, n_blades=3).generate(verbose=False)


@pytest.fixture(scope="module")
def overridden():
    return TurbineMesh(
        _YAML,
        element_size=8.0,
        n_blades=1,
        hub_diameter=10.0,
        tower_height=80.0,
        tower_base_diameter=8.0,
        tower_top_diameter=4.0,
        nacelle_length=20.0,
        nacelle_body_diameter=4.0,
        nacelle_nose_diameter=2.0,
    ).generate(verbose=False)


def test_component_names(turbine):
    assert set(turbine.meshes) == {"blade_1", "blade_2", "blade_3", "hub", "nacelle", "tower"}
    for name, mesh in turbine.meshes.items():
        assert mesh.node_count > 0, name
        assert mesh.elements_count > 0, name


def test_shared_frame(turbine):
    tower = turbine.tower
    assert tower is not None
    coords = tower.coords_array
    assert coords[:, 2].min() == pytest.approx(0.0)
    assert turbine.tower_top[2] == pytest.approx(_TOWER_TOP_Z)

    np.testing.assert_allclose(turbine.hub_center, _HUB_CENTER, atol=1e-6)
    assert turbine.hub_radius == pytest.approx(_HUB_RADIUS)
    np.testing.assert_allclose(
        turbine.rotor_axis,
        [np.cos(_UPTILT), 0.0, np.sin(_UPTILT)],
        atol=1e-12,
    )


def test_blades_root_at_hub_surface_and_reach_the_rotor_radius(turbine):
    for name, blade in turbine.blade_meshes.items():
        assert f"RootNodes_{name}" in blade.node_sets, name
        root_nodes = blade.node_sets[f"RootNodes_{name}"].nodes.values()
        root_center = np.mean([node.coords for node in root_nodes], axis=0)
        assert np.linalg.norm(root_center - turbine.hub_center) == pytest.approx(
            turbine.hub_radius, rel=0.02
        ), name

        distances = np.linalg.norm(blade.coords_array - turbine.hub_center, axis=1)
        assert distances.max() == pytest.approx(_ROTOR_RADIUS, rel=0.02), name


def test_blades_are_azimuthally_distributed(turbine):
    axis = turbine.rotor_axis
    # Build an orthonormal basis of the rotor plane.
    e1 = np.cross(axis, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    azimuths = []
    for blade in turbine.blade_meshes.values():
        offsets = blade.coords_array - turbine.hub_center
        distances = np.linalg.norm(offsets, axis=1)
        tip = offsets[int(np.argmax(distances))]
        planar = tip - np.dot(tip, axis) * axis
        azimuths.append(np.arctan2(np.dot(planar, e2), np.dot(planar, e1)))
    azimuths.sort()

    gaps = np.diff(azimuths + [azimuths[0] + 2.0 * np.pi])
    np.testing.assert_allclose(gaps, 2.0 * np.pi / 3.0, atol=0.05)


def test_no_webs_by_default(turbine):
    for blade in turbine.blade_meshes.values():
        assert not any("web" in name.lower() for name in blade.element_sets)
        assert not any("web" in name.lower() for name in blade.node_sets)


def test_include_webs_is_forwarded():
    result = TurbineMesh(_YAML, element_size=8.0, n_blades=1, include_webs=True).generate(
        verbose=False
    )
    blade = result.blade_meshes["blade_1"]
    assert "allShearWebEls_blade_1" in blade.element_sets


def test_params_override_missing_components(overridden):
    assert overridden.hub_radius == pytest.approx(5.0)
    assert overridden.tower_top[2] == pytest.approx(80.0)
    tower = overridden.tower
    assert tower is not None
    assert tower.coords_array[:, 2].max() == pytest.approx(80.0)
    # With a zero base offset the hub keeps the imported hub height.
    assert overridden.hub_center[2] == pytest.approx(150.0)
    assert overridden.hub_center[0] == pytest.approx(12.0313)


@pytest.mark.parametrize("suffix", ["stl", "vtk", "obj"])
def test_write_directory(turbine, tmp_path, suffix):
    written = turbine.write(tmp_path, format=suffix)
    assert set(written) == {"blade_1", "blade_2", "blade_3", "hub", "nacelle", "tower"}
    for name, path in written.items():
        assert os.path.getsize(path) > 0, name
