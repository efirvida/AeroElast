"""Strip-edge round-off must not split a closed physical shell ring."""

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.solvers.bem.force_projection import ForceProjector

from tests.validation.bem.test_force_projection import _make_simple_blade_aero


def _tube(edge_jitter):
    """Seven elliptical rings, including the interior edge z=1 and tip edge z=3."""
    nodes = []
    n_ring = 8
    for z in (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        for i in range(n_ring):
            angle = 2 * np.pi * i / n_ring
            jitter = edge_jitter * (-1 if i % 2 else 1) if z == 1.0 else 0.0
            nodes.append(Node([0.2 * np.cos(angle), 0.1 * np.sin(angle), z + jitter]))
    elements = []
    for ring in range(6):
        for i in range(n_ring):
            j = (i + 1) % n_ring
            elements.append(
                MeshElement(
                    [
                        nodes[ring * n_ring + i],
                        nodes[ring * n_ring + j],
                        nodes[(ring + 1) * n_ring + j],
                        nodes[(ring + 1) * n_ring + i],
                    ],
                    ElementType.quad,
                )
            )
    return MeshModel(nodes=nodes, elements=elements)


@pytest.mark.parametrize("edge_jitter", [1e-12, 0.0, 1e-4])
def test_edge_ring_is_whole_and_realisable(edge_jitter):
    mesh = _tube(edge_jitter)
    # Stations 0 and 2 give bands [-1, 1) and [1, 3], including the tip.
    aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=2.0)
    projector = ForceProjector(mesh, aero)
    edge_ring = np.arange(16, 24)
    overlap = [np.intersect1d(edge_ring, strip.node_indices).size for strip in projector._strips]
    assert overlap == [0, 8], f"strip-edge ring split across strips: {overlap}"
    for k, _strip in enumerate(projector._strips):
        for group, section in zip(
            projector._strip_ring_groups[k], projector._strip_ring_sections[k], strict=True
        ):
            assert section is not None, f"strip {k} has an unrealisable ring"
            assert len(group) == 8
    np.testing.assert_array_equal([strip.dr for strip in projector._strips], [2.0, 2.0])


def test_non_edge_nodes_and_exact_tip_keep_their_strip():
    """Away from the edge tolerance, preserve the original half-open assignment."""
    mesh = _tube(1e-12)
    aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=2.0)
    projector = ForceProjector(mesh, aero)
    span = mesh.coords_array[:, 2]
    away = np.abs(span - 1.0) > projector._ring_gap_tolerance
    expected = np.where(span < 1.0, 0, 1)
    actual = np.full(len(span), -1)
    for k, strip in enumerate(projector._strips):
        actual[strip.node_indices] = k
    np.testing.assert_array_equal(actual[away], expected[away])
    np.testing.assert_array_equal(np.bincount(actual[away]), [16, 32])
    np.testing.assert_array_equal(actual[-8:], np.ones(8, dtype=int))


def _tilted_ring_tube(tilt):
    """The tube with its z=1 ring tilted by ``tilt * cos(angle)``: one physical ring that straddles the strip edge.

    The unstructured gmsh skin mesh has such rings (odd/tasks/gmsh-official-mesher.md L-S7, L-S14): the spread is far
    beyond the edge-snapping tolerance, but the in-ring edges keep it one ring (|dz| / length < RING_EDGE_SLOPE).
    """
    mesh = _tube(0.0)
    for i, node in enumerate(mesh.nodes[16:24]):
        node.coords[2] = 1.0 + tilt * np.cos(2 * np.pi * i / 8)
    return mesh


def test_a_ring_straddling_a_strip_edge_goes_whole_to_one_strip():
    mesh = _tilted_ring_tube(0.02)
    aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=2.0)
    projector = ForceProjector(mesh, aero)
    edge_ring = np.arange(16, 24)
    overlap = [np.intersect1d(edge_ring, strip.node_indices).size for strip in projector._strips]
    # mean span 1.0 lies on the interior edge, which belongs to the higher half-open band
    assert overlap == [0, 8], f"straddling ring split across strips: {overlap}"
    for k in range(2):
        assert all(len(group) == 8 for group in projector._strip_ring_groups[k])
        assert all(section is not None for section in projector._strip_wall_flow_sections[k])


def test_a_stub_ring_does_not_take_the_wall_flow_from_its_strip():
    """A 1-2 node group (no section) is left out of the wall-flow rings instead of failing the whole strip."""
    mesh = _tube(0.0)
    mesh.nodes.append(Node([0.15, 0.0, 2.25]))  # isolated node at its own span: a one-node "ring"
    mesh = MeshModel(nodes=mesh.nodes, elements=mesh.elements)
    aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=2.0)
    projector = ForceProjector(mesh, aero)
    groups, sections = projector._strip_ring_groups[1], projector._strip_wall_flow_sections[1]
    assert sorted(len(g) for g in groups) == [1, 8, 8, 8, 8, 8]
    assert len(sections) == 5 and all(section is not None for section in sections)
