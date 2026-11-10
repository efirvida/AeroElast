"""Strip-edge round-off must not split a closed physical shell ring."""

import dataclasses

import numpy as np
import pytest

from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.solvers.bem.force_projection import ForceProjector, realise_section_load

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


def test_an_unrealisable_ring_is_skipped_and_the_others_keep_the_wall_flow():
    """A real (>= 3 node) ring with no section is skipped; its strip is not sent whole to the minimum-norm field."""
    mesh = _tube(0.0)
    aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=2.0)
    projector = ForceProjector(mesh, aero)
    strip = projector._strips[1]
    # a projector built without element_properties marks every section non-physical; flag them as physical here
    good = [dataclasses.replace(s, from_element_properties=True) for s in projector._strip_wall_flow_sections[1]]
    assert len(good) == 5
    with_bad = [good[0], good[1], None, good[3], good[4]]
    without_bad = [good[0], good[1], good[3], good[4]]
    span = np.array([0.0, 0.0, 1.0])
    F = np.array([3.0, -2.0, 5.0])
    M = np.array([1.5, -0.7, 4.0])  # span component 4 N.m carried by the wall flow, the rest transverse

    f = realise_section_load(strip, F, M, span, ring_sections=with_bad)
    mn = ForceProjector._distribute(strip, F, M)

    np.testing.assert_allclose(f.sum(axis=0), F, atol=1e-9)
    np.testing.assert_allclose(np.cross(strip.offsets, f).sum(axis=0), M, atol=1e-9)
    assert not np.allclose(f, mn, atol=1e-6), "the strip fell back to the minimum-norm field"
    # skipping the bad ring is the same as not listing it: weights renormalised over the realisable rings, the rest minimum-norm
    np.testing.assert_allclose(f, realise_section_load(strip, F, M, span, ring_sections=without_bad), atol=1e-12)

    # a strip with no realisable ring at all still falls back to the minimum-norm solve
    f_none = realise_section_load(strip, F, M, span, ring_sections=[None, None])
    np.testing.assert_allclose(f_none, mn, atol=1e-12)
