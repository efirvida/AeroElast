"""Turbine assembly: rotor centred at the origin, tower displaced behind it."""

import os

import numpy as np
import pytest

from aeroelast.core.mesh.turbine import TurbineMesh

_HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(_HERE, "IEA-15-240-RWT.yaml")

_TOWER_TOP = np.array([0.0, 12.0313, 0.0])
_NACELLE_RADIUS = 7.94 / 2.0
_ROTOR_RADIUS = 242.23775645 / 2.0


@pytest.fixture(scope="module")
def turbine():
    return TurbineMesh(_YAML, element_size=5.0, n_blades=3).generate(verbose=False)


@pytest.fixture(scope="module")
def overridden():
    return TurbineMesh(
        _YAML,
        element_size=10.0,
        n_blades=1,
        rotor_axis=(1.0, 0.0, 0.0),
        tower_offset=(-20.0, 0.0, -30.0),
        tower_height=80.0,
        tower_base_diameter=8.0,
        tower_top_diameter=4.0,
        nacelle_diameter=10.0,
    ).generate(verbose=False)


def test_component_names(turbine):
    assert set(turbine.meshes) == {"blade_1", "blade_2", "blade_3", "nacelle", "tower"}
    for name, mesh in turbine.meshes.items():
        assert mesh.node_count > 0, name
        assert mesh.elements_count > 0, name


def test_rotor_is_centred_at_the_origin(turbine):
    np.testing.assert_allclose(turbine.rotor_axis, [0.0, 1.0, 0.0])
    assert turbine.nacelle_radius == pytest.approx(_NACELLE_RADIUS)

    # One body: hemispherical hub cap centred on the rotor origin, a cylinder
    # extending backwards along +Y, and a rounded tail.
    nacelle = turbine.nacelle
    assert nacelle is not None
    assert "surface" in nacelle.element_sets
    coords = nacelle.coords_array
    assert coords[:, 1].min() == pytest.approx(-_NACELLE_RADIUS)
    # length = 2 * overhang + radius, so the body is symmetric about the tower.
    assert coords[:, 1].max() == pytest.approx(
        2.0 * np.linalg.norm(turbine.tower_top) + _NACELLE_RADIUS
    )


def test_tower_is_behind_the_rotor_plane(turbine):
    np.testing.assert_allclose(turbine.tower_top, _TOWER_TOP, atol=1e-6)
    # Behind == positive along the rotor axis (+Y); it rises to the nacelle.
    assert turbine.tower_top[1] > 0.0
    assert turbine.tower_top[0] == pytest.approx(0.0)
    assert turbine.tower_top[2] == pytest.approx(0.0)

    tower = turbine.tower
    assert tower is not None
    coords = tower.coords_array
    # The tower is vertical: its top sits at the tower_top height, base below.
    assert coords[:, 2].max() == pytest.approx(_TOWER_TOP[2])
    assert coords[:, 2].min() < _TOWER_TOP[2]
    # And it stays entirely behind the rotor plane.
    assert coords[:, 1].min() > 0.0


def test_nacelle_is_symmetric_about_the_tower_and_clears_it(turbine):
    nacelle = turbine.nacelle
    assert nacelle is not None
    y = nacelle.coords_array[:, 1]
    # The tower axis is the nacelle midpoint.
    assert (y.min() + y.max()) / 2.0 == pytest.approx(turbine.tower_top[1], abs=1e-6)
    # The tail sticks out behind the tower by more than the tower base diameter.
    assert y.max() - turbine.tower_top[1] > 10.0  # IEA-15 tower base diameter


def test_blades_root_at_the_body_surface_and_reach_the_rotor_radius(turbine):
    for name, blade in turbine.blade_meshes.items():
        assert f"RootNodes_{name}" in blade.node_sets, name
        root_nodes = blade.node_sets[f"RootNodes_{name}"].nodes.values()
        root_center = np.mean([node.coords for node in root_nodes], axis=0)
        assert np.linalg.norm(root_center) == pytest.approx(turbine.nacelle_radius, rel=0.02), name
        assert np.linalg.norm(blade.coords_array, axis=1).max() == pytest.approx(
            _ROTOR_RADIUS, rel=0.02
        ), name


def test_rotation_about_y_preserves_the_y_profile(turbine):
    # A rotation about the rotor axis (+Y) must leave every blade's y
    # coordinates identical up to reordering.  The base blade's prebend/sweep
    # lives in x/z, so y is the invariant that proves the rotation axis.
    reference = np.sort(turbine.blade_meshes["blade_1"].coords_array[:, 1])
    for name, blade in turbine.blade_meshes.items():
        np.testing.assert_allclose(
            np.sort(blade.coords_array[:, 1]), reference, atol=1e-9, err_msg=name
        )


def test_coning_tilts_the_tips_away_from_the_tower(turbine):
    # The tower sits at positive y (behind the rotor plane); coning must move
    # the blade tips the other way.
    assert turbine.tower_top[1] > 0.0
    for name, blade in turbine.blade_meshes.items():
        distances = np.linalg.norm(blade.coords_array, axis=1)
        tip = blade.coords_array[int(np.argmax(distances))]
        assert tip[1] < turbine.tower_top[1], name


def test_blades_are_azimuthally_distributed_about_y(turbine):
    azimuths = []
    for blade in turbine.blade_meshes.values():
        distances = np.linalg.norm(blade.coords_array, axis=1)
        tip = blade.coords_array[int(np.argmax(distances))]
        # Rotor plane is X-Z when rotor_axis is +Y.
        azimuths.append(np.arctan2(tip[2], tip[0]))
    azimuths.sort()
    gaps = np.diff(azimuths + [azimuths[0] + 2.0 * np.pi])
    np.testing.assert_allclose(gaps, 2.0 * np.pi / 3.0, atol=0.05)


def test_no_webs_by_default(turbine):
    for blade in turbine.blade_meshes.values():
        assert not any("web" in name.lower() for name in blade.element_sets)


def test_include_webs_is_forwarded():
    result = TurbineMesh(_YAML, element_size=10.0, n_blades=1, include_webs=True).generate(
        verbose=False
    )
    assert "allShearWebEls_blade_1" in result.blade_meshes["blade_1"].element_sets


def test_rotor_axis_is_configurable(overridden):
    np.testing.assert_allclose(overridden.rotor_axis, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(overridden.tower_top, [-20.0, 0.0, -30.0], atol=1e-9)
    assert overridden.nacelle_radius == pytest.approx(5.0)
    # The tower follows the requested base height.
    tower = overridden.tower
    assert tower is not None
    assert tower.coords_array[:, 2].min() == pytest.approx(-30.0 - 80.0)


def test_tower_offset_is_used_verbatim():
    from aeroelast.core.mesh.components import read_windio_components

    mesh = TurbineMesh(_YAML, tower_offset=(1.0, 2.0, 3.0), element_size=10.0, n_blades=1)
    definition = read_windio_components(_YAML)
    np.testing.assert_allclose(mesh._resolve_tower_top(definition), [1.0, 2.0, 3.0])


@pytest.mark.parametrize("suffix", ["stl", "vtk", "obj"])
def test_write_directory(turbine, tmp_path, suffix):
    written = turbine.write(tmp_path, format=suffix)
    assert set(written) == {"blade_1", "blade_2", "blade_3", "nacelle", "tower"}
    for name, path in written.items():
        assert os.path.getsize(path) > 0, name
