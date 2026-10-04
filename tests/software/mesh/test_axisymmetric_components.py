"""Axisymmetric turbine component meshes: tower and the nacelle body."""

import os
from collections import Counter

import numpy as np
import pytest

from aeroelast.core.mesh import (
    NacelleMesh,
    TowerMesh,
    build_revolved_shell,
    read_windio_components,
)
from aeroelast.core.mesh.io.writers import write_mesh

_HERE = os.path.dirname(os.path.abspath(__file__))
from tests.support.paths import DATA_DIR  # noqa: E402
_YAML = str(DATA_DIR / "IEA-15-240-RWT.yaml")

# IEA-15-240-RWT reference values (tests/IEA-15-240-RWT.yaml).
_TOWER_Z_BASE = 15.0
_TOWER_Z_TOP = 144.386
_TOWER_DIAMETER_BASE = 10.0
_TOWER_DIAMETER_TOP = 6.5
_HUB_DIAMETER = 7.94
_BODY_RADIUS = _HUB_DIAMETER / 2.0
_BODY_LENGTH = 12.0313  # overhang (yaw axis -> hub)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _closed_volume(mesh) -> float:
    total = 0.0
    for element in mesh.elements:
        pts = [node.coords for node in element.nodes]
        if len(pts) == 3:
            triangles = [pts]
        elif len(pts) == 4:
            triangles = [[pts[0], pts[1], pts[2]], [pts[0], pts[2], pts[3]]]
        else:
            raise AssertionError(f"unexpected element with {len(pts)} nodes")
        for a, b, c in triangles:
            total += float(np.dot(a, np.cross(b, c)))
    return total / 6.0


def _assert_closed(mesh) -> None:
    edges: Counter = Counter()
    for element in mesh.elements:
        ids = list(element.node_ids)
        for i, node_id in enumerate(ids):
            other = ids[(i + 1) % len(ids)]
            edges[tuple(sorted((node_id, other)))] += 1
    bad = {edge: c for edge, c in edges.items() if c != 2}
    assert not bad, f"surface is not watertight: {len(bad)} bad edges"


def _assert_outward(mesh) -> None:
    volume = _closed_volume(mesh)
    assert volume > 0.0, f"orientation is inward (signed volume {volume})"


# ---------------------------------------------------------------------------
# WindIO reader
# ---------------------------------------------------------------------------


def test_read_windio_components_parses_all_blocks():
    definition = read_windio_components(_YAML)

    assert definition.number_of_blades == 3
    assert definition.hub_height == pytest.approx(150.0)
    assert definition.rotor_diameter == pytest.approx(242.23775645)

    tower = definition.tower
    assert tower is not None
    assert tower.z[0] == pytest.approx(_TOWER_Z_BASE)
    assert tower.z[-1] == pytest.approx(_TOWER_Z_TOP)
    assert tower.diameter[0] == pytest.approx(_TOWER_DIAMETER_BASE)
    assert tower.diameter[-1] == pytest.approx(_TOWER_DIAMETER_TOP)

    hub = definition.hub
    assert hub is not None
    assert hub.diameter == pytest.approx(_HUB_DIAMETER)
    assert hub.cone_angle == pytest.approx(0.06981317007977318)

    nacelle = definition.nacelle
    assert nacelle is not None
    assert nacelle.body_diameter == pytest.approx(3.0)
    assert nacelle.nose_diameter == pytest.approx(2.2)
    assert nacelle.length == pytest.approx(_BODY_LENGTH)
    assert nacelle.distance_tt_hub == pytest.approx(5.614)


def test_from_windio_raises_when_component_missing(tmp_path):
    trimmed = tmp_path / "hub_only.yaml"
    trimmed.write_text("components:\n  hub:\n    diameter: 4.0\nassembly:\n  number_of_blades: 2\n")

    with pytest.raises(ValueError, match="no tower definition"):
        TowerMesh.from_windio(trimmed)
    # The body needs both the hub diameter and a nacelle length.
    with pytest.raises(ValueError, match="no nacelle length"):
        NacelleMesh.from_windio(trimmed)

    definition = read_windio_components(trimmed)
    assert definition.tower is None
    assert definition.nacelle is None
    assert definition.hub is not None
    assert definition.hub.diameter == 4.0


# ---------------------------------------------------------------------------
# Generic builder
# ---------------------------------------------------------------------------


def test_revolved_cylinder_is_closed_and_outward():
    centers = [[0, 0, 0], [0, 0, 1]]
    mesh = build_revolved_shell(centers, [0.5, 0.5], n_circ=8, cap_start=True, cap_end=True)
    assert mesh.node_count == 2 * 8 + 2
    assert mesh.elements_count == 8 + 2 * 8
    _assert_closed(mesh)
    _assert_outward(mesh)

    refined = build_revolved_shell(centers, [0.5, 0.5], n_circ=64, cap_start=True, cap_end=True)
    assert _closed_volume(refined) == pytest.approx(np.pi * 0.25 * 1.0, rel=1e-2)


def test_revolved_shell_rejects_bad_input():
    with pytest.raises(ValueError, match="n_circ"):
        build_revolved_shell([[0, 0, 0], [0, 0, 1]], [1.0, 1.0], n_circ=2)
    with pytest.raises(ValueError, match="same length"):
        build_revolved_shell([[0, 0, 0], [0, 0, 1]], [1.0])
    with pytest.raises(ValueError, match="non-zero"):
        build_revolved_shell([[0, 0, 0], [0, 0, 1]], [1.0, 1.0], axis=(0, 0, 0))
    with pytest.raises(ValueError, match="two consecutive"):
        build_revolved_shell([[0, 0, 0], [0, 0, 1]], [0.0, 0.0])


# ---------------------------------------------------------------------------
# Tower
# ---------------------------------------------------------------------------


def test_tower_from_windio_geometry():
    tower = TowerMesh.from_windio(_YAML, n_circ=8, n_axial=20)
    mesh = tower.generate()

    assert mesh.node_count == 21 * 8 + 2
    assert mesh.elements_count == 20 * 8 + 2 * 8
    assert set(mesh.element_sets) == {"surface", "cap_start", "cap_end"}
    assert {"base", "top", "all"} <= set(mesh.node_sets)

    coords = mesh.coords_array
    assert coords[:, 2].min() == pytest.approx(_TOWER_Z_BASE)
    assert coords[:, 2].max() == pytest.approx(_TOWER_Z_TOP)

    radius = np.linalg.norm(coords[:, :2], axis=1)
    base_ring = coords[np.isclose(coords[:, 2], _TOWER_Z_BASE)]
    assert np.linalg.norm(base_ring[:, :2], axis=1).max() == pytest.approx(
        _TOWER_DIAMETER_BASE / 2.0
    )
    assert radius.max() == pytest.approx(_TOWER_DIAMETER_BASE / 2.0)

    _assert_closed(mesh)
    _assert_outward(mesh)


def test_tower_element_size_controls_axial_count():
    tower = TowerMesh.from_windio(_YAML, n_circ=8, element_size=2.0)
    # ceil((144.386 - 15) / 2) = 65 axial segments.
    assert tower.stations.shape[0] == 66


def test_tower_from_params_is_a_linear_frustum():
    tower = TowerMesh.from_params(
        height=100.0, base_diameter=10.0, top_diameter=5.0, n_axial=10, n_circ=12
    )
    mesh = tower.generate()

    assert mesh.node_count == 11 * 12 + 2
    assert mesh.elements_count == 10 * 12 + 2 * 12
    coords = mesh.coords_array
    assert coords[:, 2].min() == pytest.approx(0.0)
    assert coords[:, 2].max() == pytest.approx(100.0)
    top_ring = coords[np.isclose(coords[:, 2], 100.0)]
    assert np.linalg.norm(top_ring[:, :2], axis=1).max() == pytest.approx(2.5)
    _assert_closed(mesh)
    _assert_outward(mesh)


# ---------------------------------------------------------------------------
# Hub + nacelle body
# ---------------------------------------------------------------------------


def test_body_from_windio_is_two_hemispheres_and_a_cylinder():
    body = NacelleMesh.from_windio(_YAML, n_circ=8, n_axial=6, n_tip=4, n_tail=4)
    mesh = body.generate()

    # Default axis is +Y: the hub cap spans [-r, 0], the cylinder [0, L-r] and
    # the drag-reducing tail hemisphere [L-r, L].
    coords = mesh.coords_array
    assert coords[:, 1].min() == pytest.approx(-_BODY_RADIUS)
    assert coords[:, 1].max() == pytest.approx(_BODY_LENGTH)

    radial = np.hypot(coords[:, 0], coords[:, 2])
    assert radial.max() == pytest.approx(_BODY_RADIUS)

    # Constant radius only on the cylinder, between the two hemispheres.
    cylinder = coords[
        (coords[:, 1] >= -1e-9) & (coords[:, 1] <= _BODY_LENGTH - _BODY_RADIUS + 1e-9)
    ]
    assert cylinder.shape[0] > 0
    assert np.allclose(np.hypot(cylinder[:, 0], cylinder[:, 2]), _BODY_RADIUS)

    # Both ends close in a single pole on the axis.
    for y_pole in (-_BODY_RADIUS, _BODY_LENGTH):
        pole = coords[np.isclose(coords[:, 1], y_pole)]
        assert pole.shape[0] == 1
        assert np.hypot(pole[0, 0], pole[0, 2]) == pytest.approx(0.0, abs=1e-12)

    assert {"tail", "nose", "all"} <= set(mesh.node_sets)
    assert "cap_end" not in mesh.element_sets  # rounded tail, no flat cap
    _assert_closed(mesh)
    _assert_outward(mesh)


def test_body_flat_tail_is_supported():
    body = NacelleMesh.from_params(length=12.0, radius=2.0, n_circ=8, n_axial=4, rear_tip=False)
    mesh = body.generate()
    assert "cap_end" in mesh.element_sets
    _assert_closed(mesh)
    _assert_outward(mesh)


def test_body_volume_matches_the_analytic_cylinder_plus_two_hemispheres():
    length, radius = 20.0, 3.0
    body = NacelleMesh(length=length, radius=radius, n_circ=48, n_axial=8, n_tip=8, n_tail=8)
    mesh = body.generate()
    exact = np.pi * radius**2 * (length - radius) + 2.0 * (2.0 / 3.0) * np.pi * radius**3
    assert _closed_volume(mesh) == pytest.approx(exact, rel=1e-2)


def test_body_from_params_axis_and_center():
    body = NacelleMesh.from_params(
        length=10.0, diameter=2.0, n_circ=8, n_axial=4, n_tip=3, center=(1.0, 2.0, 3.0)
    )
    mesh = body.generate()
    coords = mesh.coords_array
    # Default axis is +Y, so the body runs along y: hemisphere at 2-1, cylinder
    # from 2 to 12.
    assert coords[:, 1].min() == pytest.approx(1.0)
    assert coords[:, 1].max() == pytest.approx(12.0)
    assert np.hypot(coords[:, 0] - 1.0, coords[:, 2] - 3.0).max() == pytest.approx(1.0)
    _assert_closed(mesh)


def test_body_rejects_bad_parameters():
    with pytest.raises(ValueError, match="either radius or diameter"):
        NacelleMesh(length=10.0)
    with pytest.raises(ValueError, match="not both"):
        NacelleMesh(length=10.0, radius=1.0, diameter=2.0)
    with pytest.raises(ValueError, match="must be positive"):
        NacelleMesh(length=-1.0, radius=1.0)


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("suffix", [".stl", ".obj", ".vtk"])
def test_components_export(tmp_path, suffix):
    meshes = {
        "tower": TowerMesh.from_params(height=50.0, base_diameter=6.0, n_axial=4).generate(),
        "nacelle": NacelleMesh.from_params(
            length=8.0, diameter=4.0, n_circ=8, n_axial=4
        ).generate(),
    }
    for name, mesh in meshes.items():
        target = tmp_path / f"{name}{suffix}"
        write_mesh(mesh, str(target))
        assert target.exists()
        assert target.stat().st_size > 0
