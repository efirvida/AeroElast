"""Tests for reduced aerodynamic surface extraction on a full rotor."""

from __future__ import annotations

import numpy as np

from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.core.mesh.generators import RotorMesh
from aeroelast.solvers.aero import build_full_rotor_aero_geometry, build_full_rotor_aero_surface


def _build_surface_toy_rotor_mesh() -> MeshModel:
    mesh = MeshModel()

    blade_1_nodes = [
        Node([0.0, 0.0, 1.0]),
        Node([1.0, 0.0, 1.0]),
        Node([1.0, 0.0, 2.0]),
        Node([0.0, 0.0, 2.0]),
    ]
    blade_2_nodes = [
        Node([0.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -2.0]),
        Node([0.0, 0.0, -2.0]),
    ]

    for node in blade_1_nodes + blade_2_nodes:
        mesh.add_node(node)

    blade_1_elem = MeshElement(blade_1_nodes, ElementType.quad)
    blade_2_elem = MeshElement(blade_2_nodes, ElementType.quad)
    mesh.add_element(blade_1_elem)
    mesh.add_element(blade_2_elem)

    mesh.add_node_set(NodeSet("rotor_blade_1", set(blade_1_nodes)))
    mesh.add_node_set(NodeSet("rotor_blade_2", set(blade_2_nodes)))
    mesh.add_element_set(ElementSet("rotor_blade_1", {blade_1_elem}))
    mesh.add_element_set(ElementSet("rotor_blade_2", {blade_2_elem}))

    return mesh


def test_build_full_rotor_aero_surface_extracts_panel_metrics():
    geometry = build_full_rotor_aero_geometry(_build_surface_toy_rotor_mesh(), hub_radius=1.0)

    surface = build_full_rotor_aero_surface(geometry)

    assert surface.n_blades == 2
    assert surface.n_panels == 2

    blade_1 = surface.get_blade(0)
    blade_2 = surface.get_blade(1)

    assert blade_1.n_panels == 1
    assert blade_2.n_panels == 1
    np.testing.assert_allclose(blade_1.panel_areas, [1.0])
    np.testing.assert_allclose(blade_2.panel_areas, [1.0])
    np.testing.assert_allclose(blade_1.panel_centroids[0], [0.5, 0.0, 1.5])
    np.testing.assert_allclose(blade_2.panel_centroids[0], [-0.5, 0.0, -1.5])
    np.testing.assert_allclose(blade_1.span_direction, [0.0, 0.0, 1.0])
    np.testing.assert_allclose(blade_1.chord_direction, [1.0, 0.0, 0.0])
    np.testing.assert_allclose(blade_1.thickness_direction, [0.0, 1.0, 0.0])
    np.testing.assert_allclose(blade_1.panel_normals[0], [0.0, 1.0, 0.0])
    np.testing.assert_allclose(blade_2.panel_normals[0], [0.0, 1.0, 0.0])
    np.testing.assert_allclose(blade_1.root_reference_radius, 1.0)
    np.testing.assert_allclose(blade_2.root_reference_radius, 1.0)
    np.testing.assert_allclose(blade_1.node_reference_radii, [1.0, 1.0, 2.0, 2.0])
    np.testing.assert_allclose(blade_2.node_reference_radii, [1.0, 1.0, 2.0, 2.0])
    np.testing.assert_allclose(blade_1.panel_reference_radii, [1.5])
    np.testing.assert_allclose(blade_2.panel_reference_radii, [1.5])
    np.testing.assert_allclose(blade_1.panel_spanwise_coordinates, [0.5])
    np.testing.assert_allclose(blade_2.panel_spanwise_coordinates, [0.5])
    np.testing.assert_allclose(blade_1.panel_span_lengths, [1.0])
    np.testing.assert_allclose(blade_2.panel_span_lengths, [1.0])
    np.testing.assert_allclose(blade_1.panel_chord_lengths, [1.0])
    np.testing.assert_allclose(blade_2.panel_chord_lengths, [1.0])
    np.testing.assert_allclose(blade_1.panel_collocation_points[0], [0.75, 0.0, 1.5])
    np.testing.assert_allclose(
        blade_1.panel_bound_vortex_points[0],
        [[0.25, 0.0, 1.0], [0.25, 0.0, 2.0]],
    )
    np.testing.assert_array_equal(blade_1.parent_element_indices, [0])
    np.testing.assert_array_equal(blade_2.parent_element_indices, [1])


def test_build_full_rotor_aero_surface_tracks_parent_panel_nodes():
    geometry = build_full_rotor_aero_geometry(_build_surface_toy_rotor_mesh())

    surface = build_full_rotor_aero_surface(geometry)

    blade_1 = surface.get_blade(0)
    blade_2 = surface.get_blade(1)
    np.testing.assert_array_equal(
        blade_1.panel_parent_node_indices[0],
        geometry.blades[0].parent_node_indices[blade_1.panel_node_indices[0]],
    )
    np.testing.assert_array_equal(
        blade_2.panel_parent_node_indices[0],
        geometry.blades[1].parent_node_indices[blade_2.panel_node_indices[0]],
    )


def test_full_rotor_aero_surface_assembles_global_panel_forces():
    geometry = build_full_rotor_aero_geometry(_build_surface_toy_rotor_mesh(), hub_radius=1.0)
    surface = build_full_rotor_aero_surface(geometry)

    global_forces = surface.assemble_global_nodal_forces([
        np.array([[0.0, 8.0, 0.0]]),
        np.array([[0.0, 4.0, 0.0]]),
    ])

    expected = np.zeros((geometry.rotor_mesh.node_count, 3), dtype=float)
    expected[geometry.blades[0].parent_node_indices, 1] = 2.0
    expected[geometry.blades[1].parent_node_indices, 1] = 1.0
    np.testing.assert_allclose(global_forces, expected)


def test_full_rotor_aero_surface_assembles_global_forces_from_pressure():
    geometry = build_full_rotor_aero_geometry(_build_surface_toy_rotor_mesh(), hub_radius=1.0)
    surface = build_full_rotor_aero_surface(geometry)

    global_forces = surface.assemble_global_nodal_forces_from_pressure([
        np.array([2.0]),
        np.array([1.0]),
    ])

    expected = np.zeros((geometry.rotor_mesh.node_count, 3), dtype=float)
    expected[geometry.blades[0].parent_node_indices, 1] = 0.5
    expected[geometry.blades[1].parent_node_indices, 1] = 0.25
    np.testing.assert_allclose(global_forces, expected)


def test_full_rotor_aero_surface_builds_deformed_state_from_global_displacements():
    geometry = build_full_rotor_aero_geometry(_build_surface_toy_rotor_mesh(), hub_radius=1.0)
    surface = build_full_rotor_aero_surface(geometry)

    displacements = np.zeros((geometry.rotor_mesh.node_count, 3), dtype=float)
    displacements[geometry.blades[0].parent_node_indices] = [0.1, 0.2, 0.3]

    deformed = surface.build_deformed_state(displacements)

    blade_1 = deformed.get_blade(0)
    blade_2 = deformed.get_blade(1)
    np.testing.assert_allclose(
        blade_1.node_coordinates,
        surface.blades[0].reference_node_coordinates + [0.1, 0.2, 0.3],
    )
    np.testing.assert_allclose(blade_1.panel_centroids[0], [0.6, 0.2, 1.8])
    np.testing.assert_allclose(blade_1.panel_collocation_points[0], [0.85, 0.2, 2.1])
    np.testing.assert_allclose(
        blade_1.panel_bound_vortex_points[0],
        [[0.35, 0.2, 1.6], [0.35, 0.2, 2.6]],
    )
    np.testing.assert_allclose(blade_1.panel_normals[0], [0.0, 1.0, 0.0])
    np.testing.assert_allclose(blade_2.panel_centroids[0], surface.blades[1].panel_centroids[0])


def test_build_full_rotor_aero_surface_from_rotormesh(iea_blade_yaml: str):
    rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=3,
        hub_radius=None,
        element_size=12.0,
        n_samples=60,
    )
    mesh = rotor.generate(renumber=None, verbose=False)

    geometry = build_full_rotor_aero_geometry(mesh)
    surface = build_full_rotor_aero_surface(geometry)

    assert surface.n_blades == 3
    assert surface.n_panels == mesh.elements_count
    for blade_surface, blade_geometry in zip(surface.blades, geometry.blades):
        assert blade_surface.n_panels == blade_geometry.element_count
        assert np.all(blade_surface.panel_areas > 0.0)
        assert blade_surface.root_reference_radius >= 0.0
        assert np.all(blade_surface.panel_reference_radii >= blade_surface.root_reference_radius)
