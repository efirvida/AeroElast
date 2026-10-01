"""Axisymmetric turbine component meshes: tower, hub and nacelle."""

import os
from collections import Counter

import numpy as np
import pytest

from aeroelast.core.mesh import (
    HubMesh,
    NacelleMesh,
    TowerMesh,
    build_revolved_shell,
    read_windio_components,
)
from aeroelast.core.mesh.io.writers import write_mesh

_HERE = os.path.dirname(os.path.abspath(__file__))
_YAML = os.path.join(_HERE, "IEA-15-240-RWT.yaml")

# IEA-15-240-RWT reference values (tests/IEA-15-240-RWT.yaml).
_TOWER_Z_BASE = 15.0
_TOWER_Z_TOP = 144.386
_TOWER_DIAMETER_BASE = 10.0
_TOWER_DIAMETER_TOP = 6.5
_HUB_DIAMETER = 7.94
_NACELLE_BODY_DIAMETER = 3.0
_NACELLE_NOSE_DIAMETER = 2.2
_NACELLE_LENGTH = 5.614
_UPTILT = 0.10471975511965977
_OVERHANG = 12.0313


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


def _edge_counts(mesh) -> Counter:
    edges: Counter = Counter()
    for element in mesh.elements:
        ids = list(element.node_ids)
        for i, node_id in enumerate(ids):
            other = ids[(i + 1) % len(ids)]
            edges[tuple(sorted((node_id, other)))] += 1
    return edges


def _assert_closed(mesh) -> None:
    counts = _edge_counts(mesh)
    non_manifold = {edge: c for edge, c in counts.items() if c != 2}
    assert not non_manifold, f"surface is not watertight: {len(non_manifold)} bad edges"


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
    assert nacelle.body_diameter == pytest.approx(_NACELLE_BODY_DIAMETER)
    assert nacelle.nose_diameter == pytest.approx(_NACELLE_NOSE_DIAMETER)
    assert nacelle.length == pytest.approx(_NACELLE_LENGTH)
    assert nacelle.uptilt == pytest.approx(_UPTILT)
    assert nacelle.overhang == pytest.approx(_OVERHANG)


def test_from_windio_raises_when_component_missing(tmp_path):
    trimmed = tmp_path / "hub_only.yaml"
    trimmed.write_text("components:\n  hub:\n    diameter: 4.0\nassembly:\n  number_of_blades: 2\n")

    with pytest.raises(ValueError, match="no tower definition"):
        TowerMesh.from_windio(trimmed)
    with pytest.raises(ValueError, match="no nacelle definition"):
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

    # A discretised cylinder inscribes the exact one; refine to compare volumes.
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
    base_radius = np.linalg.norm(base_ring[:, :2], axis=1)
    assert base_radius.max() == pytest.approx(_TOWER_DIAMETER_BASE / 2.0)
    assert radius.max() == pytest.approx(_TOWER_DIAMETER_BASE / 2.0)

    _assert_closed(mesh)
    _assert_outward(mesh)


def test_tower_element_size_controls_axial_count():
    tower = TowerMesh.from_windio(_YAML, n_circ=8, element_size=2.0)
    # ceil((144.386 - 15) / 2) = 65 axial segments.
    assert tower.stations.shape[0] == 66


def test_tower_from_params_is_a_linear_frustum():
    tower = TowerMesh.from_params(
        height=100.0,
        base_diameter=10.0,
        top_diameter=5.0,
        n_axial=10,
        n_circ=12,
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
# Hub
# ---------------------------------------------------------------------------


def test_hub_from_windio_is_the_expected_sphere():
    hub = HubMesh.from_windio(_YAML, n_circ=8, n_merid=16)
    mesh = hub.generate()

    # 17 stations, two of them poles.
    assert mesh.node_count == 15 * 8 + 2
    assert mesh.elements_count == 16 * 8
    assert set(mesh.element_sets) == {"surface"}
    assert {"tail", "nose", "all"} <= set(mesh.node_sets)

    coords = mesh.coords_array
    assert coords.max() == pytest.approx(_HUB_DIAMETER / 2.0)
    assert coords.min() == pytest.approx(-_HUB_DIAMETER / 2.0)

    _assert_closed(mesh)
    _assert_outward(mesh)


def test_hub_volume_matches_the_exact_sphere():
    hub = HubMesh(diameter=_HUB_DIAMETER, n_circ=64, n_merid=32)
    mesh = hub.generate()
    exact = 4.0 / 3.0 * np.pi * (_HUB_DIAMETER / 2.0) ** 3
    assert _closed_volume(mesh) == pytest.approx(exact, rel=1e-2)


def test_hub_from_params_axis_and_center():
    hub = HubMesh.from_params(diameter=2.0, n_circ=8, n_merid=8, center=(1.0, 2.0, 3.0))
    mesh = hub.generate()
    coords = mesh.coords_array
    assert coords[:, 0].min() == pytest.approx(0.0)
    assert coords[:, 0].max() == pytest.approx(2.0)
    assert coords[:, 1].min() == pytest.approx(1.0)
    assert coords[:, 2].min() == pytest.approx(2.0)
    assert coords[:, 2].max() == pytest.approx(4.0)
    _assert_closed(mesh)


# ---------------------------------------------------------------------------
# Nacelle
# ---------------------------------------------------------------------------


def test_nacelle_from_windio_envelope():
    nacelle = NacelleMesh.from_windio(_YAML, n_circ=8, n_axial=6, n_nose=4)
    mesh = nacelle.generate()

    coords = mesh.coords_array
    assert coords[:, 0].min() == pytest.approx(0.0)
    assert coords[:, 0].max() == pytest.approx(_NACELLE_LENGTH)
    assert coords[:, 1].min() == pytest.approx(-_NACELLE_BODY_DIAMETER / 2.0)
    assert coords[:, 1].max() == pytest.approx(_NACELLE_BODY_DIAMETER / 2.0)
    assert coords[:, 2].min() == pytest.approx(-_NACELLE_BODY_DIAMETER / 2.0)
    assert coords[:, 2].max() == pytest.approx(_NACELLE_BODY_DIAMETER / 2.0)

    assert {"tail", "nose", "all"} <= set(mesh.node_sets)
    assert "cap_start" in mesh.element_sets
    _assert_closed(mesh)
    _assert_outward(mesh)


def test_nacelle_from_params_axis():
    nacelle = NacelleMesh.from_params(
        length=10.0,
        body_diameter=2.0,
        nose_diameter=1.0,
        n_circ=8,
        n_axial=4,
        n_nose=3,
        axis=(0.0, 1.0, 0.0),
    )
    mesh = nacelle.generate()
    coords = mesh.coords_array
    assert coords[:, 1].min() == pytest.approx(0.0)
    assert coords[:, 1].max() == pytest.approx(10.0)
    _assert_closed(mesh)
    _assert_outward(mesh)


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("suffix", [".stl", ".obj", ".vtk"])
def test_components_export(tmp_path, suffix):
    meshes = {
        "tower": TowerMesh.from_params(height=50.0, base_diameter=6.0, n_axial=4).generate(),
        "hub": HubMesh(diameter=4.0, n_circ=8, n_merid=8).generate(),
        "nacelle": NacelleMesh.from_params(
            length=6.0, body_diameter=2.4, n_axial=4, n_nose=3
        ).generate(),
    }
    for name, mesh in meshes.items():
        target = tmp_path / f"{name}{suffix}"
        write_mesh(mesh, str(target))
        assert target.exists()
        assert target.stat().st_size > 0
