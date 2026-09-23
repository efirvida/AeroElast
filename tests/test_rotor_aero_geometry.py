"""Tests for full-rotor aerodynamic geometry extraction helpers."""

from __future__ import annotations

import numpy as np

from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.core.mesh.generators import RotorMesh
from aeroelast.solvers.aero import build_full_rotor_aero_geometry


def _build_toy_rotor_mesh() -> MeshModel:
    mesh = MeshModel()

    blade_1_nodes = [
        Node([0.0, 0.0, 1.0]),
        Node([1.0, 0.0, 1.0]),
        Node([1.0, 1.0, 1.0]),
        Node([0.0, 1.0, 1.0]),
    ]
    blade_2_nodes = [
        Node([0.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -1.0]),
        Node([-1.0, 1.0, -1.0]),
        Node([0.0, 1.0, -1.0]),
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


def test_build_full_rotor_aero_geometry_detects_blades_and_azimuths():
    mesh = _build_toy_rotor_mesh()

    geometry = build_full_rotor_aero_geometry(mesh)

    assert geometry.n_blades == 2
    np.testing.assert_allclose(geometry.rotation_axis, [0.0, 1.0, 0.0])
    np.testing.assert_allclose(geometry.rotation_center, [0.0, 0.0, 0.0])
    assert geometry.hub_radius is None
    assert [blade.element_set_name for blade in geometry.blades] == [
        "rotor_blade_1",
        "rotor_blade_2",
    ]
    assert [blade.azimuth_deg for blade in geometry.blades] == [0.0, 180.0]
    assert geometry.blades[0].mesh.elements_count == 1
    assert geometry.blades[1].mesh.elements_count == 1


def test_blade_deformation_uses_parent_node_order_and_full_rotor_field():
    mesh = _build_toy_rotor_mesh()
    geometry = build_full_rotor_aero_geometry(mesh)

    displacements = np.zeros((mesh.node_count, 6), dtype=float)
    blade = geometry.get_blade(1)
    displacements[blade.parent_node_indices, 0] = 0.25
    displacements[blade.parent_node_indices, 2] = -0.10

    deformed_mesh = blade.get_deformed_mesh(
        displacements,
        dofs_per_node=6,
        displacement_indices=(0, 1, 2),
    )

    expected_coords = blade.reference_coords.copy()
    expected_coords[:, 0] += 0.25
    expected_coords[:, 2] -= 0.10

    np.testing.assert_allclose(deformed_mesh.coords_array, expected_coords)


def test_full_rotor_aero_geometry_splits_blade_fields():
    mesh = _build_toy_rotor_mesh()
    geometry = build_full_rotor_aero_geometry(mesh)

    field = np.arange(mesh.node_count * 3, dtype=float).reshape(mesh.node_count, 3)
    blade_fields = geometry.extract_blade_fields(field)

    assert len(blade_fields) == 2
    np.testing.assert_array_equal(
        blade_fields[0],
        field[geometry.blades[0].parent_node_indices],
    )
    np.testing.assert_array_equal(
        blade_fields[1],
        field[geometry.blades[1].parent_node_indices],
    )


def test_full_rotor_aero_geometry_assembles_global_fields():
    mesh = _build_toy_rotor_mesh()
    geometry = build_full_rotor_aero_geometry(mesh)

    blade_fields = tuple(
        np.full((blade.node_count, 3), blade.blade_index + 1.0, dtype=float)
        for blade in geometry.blades
    )

    assembled = geometry.assemble_global_nodal_field(blade_fields)

    expected = np.zeros((mesh.node_count, 3), dtype=float)
    for blade, blade_field in zip(geometry.blades, blade_fields):
        expected[blade.parent_node_indices] = blade_field

    np.testing.assert_array_equal(assembled, expected)


def test_build_full_rotor_aero_geometry_stores_rotation_center_and_hub_radius():
    mesh = _build_toy_rotor_mesh()

    geometry = build_full_rotor_aero_geometry(
        mesh,
        rotation_center=(1.0, 2.0, 3.0),
        hub_radius=4.5,
    )

    np.testing.assert_allclose(geometry.rotation_center, [1.0, 2.0, 3.0])
    assert geometry.hub_radius == 4.5


def test_build_full_rotor_aero_geometry_from_rotormesh(iea_blade_yaml: str):
    rotor = RotorMesh(
        yaml_file=iea_blade_yaml,
        n_blades=3,
        hub_radius=None,
        element_size=12.0,
        n_samples=60,
    )
    mesh = rotor.generate(renumber=None, verbose=False)

    geometry = build_full_rotor_aero_geometry(mesh)

    assert geometry.n_blades == 3
    assert sum(blade.mesh.elements_count for blade in geometry.blades) == mesh.elements_count
    for blade_index, blade in enumerate(geometry.blades):
        assert blade.element_set_name == RotorMesh.blade_part_set_name(blade_index)
        assert blade.node_count > 0
