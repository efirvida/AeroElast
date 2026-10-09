"""Unit tests for ``MeshModel.subset_to_nodes`` (node-subset view keeping elements).

The mesh is built by hand so this module stays dependency-free (no ccblade, no
Rust extension, no file IO).
"""

import warnings

import pytest

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel

# Node layout: a 3x3 skin grid at z=0 plus one web-only node at index 9.
_SKIN_COORDS = [
    (0.0, 0.0, 0.0),  # 0
    (1.0, 0.0, 0.0),  # 1
    (2.0, 0.0, 0.0),  # 2
    (0.0, 1.0, 0.0),  # 3
    (1.0, 1.0, 0.0),  # 4
    (2.0, 1.0, 0.0),  # 5
    (0.0, 2.0, 0.0),  # 6
    (1.0, 2.0, 0.0),  # 7
    (2.0, 2.0, 0.0),  # 8
    (0.5, 0.5, 0.5),  # 9 (web-only, never in the requested set)
]

# A 2x2 quad skin (four quads, source order) then a web triangle using node 9.
_QUADS = [(0, 1, 4, 3), (1, 2, 5, 4), (3, 4, 7, 6), (4, 5, 8, 7)]
_WEB = (0, 1, 9)

# Request every skin node in an order that deliberately differs from source order.
_REQUESTED_INDEX_ORDER = [2, 0, 1, 3, 4, 5, 6, 7, 8]


def _build_mesh() -> MeshModel:
    """Build the synthetic skin + web mesh used by every test."""
    mesh = MeshModel()
    nodes = []
    for coords in _SKIN_COORDS:
        node = Node(list(coords))
        mesh.add_node(node)
        nodes.append(node)

    quads = []
    for conn in _QUADS:
        quad = MeshElement([nodes[i] for i in conn], ElementType.quad)
        mesh.add_element(quad)
        quads.append(quad)
    web = MeshElement([nodes[i] for i in _WEB], ElementType.triangle)
    mesh.add_element(web)

    mesh.add_node_set(NodeSet("corner_and_web", {nodes[0], nodes[4], nodes[9]}))
    mesh.add_node_set(NodeSet("web_nodes", {nodes[9]}))
    mesh.add_element_set(ElementSet("skin_and_web", {quads[0], quads[1], web}))
    mesh.add_element_set(ElementSet("web_only", {web}))
    return mesh


def _requested(mesh: MeshModel) -> list[int]:
    """The requested node IDs: all skin nodes, in ``_REQUESTED_INDEX_ORDER``."""
    return [mesh.nodes[i].id for i in _REQUESTED_INDEX_ORDER]


def _coords_index(mesh: MeshModel) -> dict:
    """Map node coordinates to their position (== source index) in the mesh."""
    return {tuple(node.coords.tolist()): idx for idx, node in enumerate(mesh.nodes)}


def _source_index_connectivity(mesh: MeshModel, result: MeshModel) -> list:
    """Result element connectivity expressed in source-mesh node positions."""
    lookup = _coords_index(mesh)
    return [
        tuple(lookup[tuple(node.coords.tolist())] for node in element.nodes)
        for element in result.elements
    ]


def _expected_connectivity(mesh: MeshModel, requested_indices: list) -> list:
    """Source elements fully contained in ``requested_indices``, in source order."""
    id_to_index = {node.id: idx for idx, node in enumerate(mesh.nodes)}
    requested_set = set(requested_indices)
    expected = []
    for element in mesh.elements:
        positions = tuple(id_to_index[nid] for nid in element.node_ids)
        if all(pos in requested_set for pos in positions):
            expected.append(positions)
    return expected


def test_requested_node_order_and_coords_are_preserved():
    mesh = _build_mesh()
    requested = _requested(mesh)
    result = mesh.subset_to_nodes(requested)

    assert len(result.nodes) == len(_REQUESTED_INDEX_ORDER)
    assert [tuple(node.coords.tolist()) for node in result.nodes] == [
        _SKIN_COORDS[i] for i in _REQUESTED_INDEX_ORDER
    ]
    # Ids are the selection key and stay those of the source nodes.
    assert [node.id for node in result.nodes] == requested
    # A repeated id is honoured once, in its first position.
    deduped = mesh.subset_to_nodes([requested[2], *requested])
    assert [node.id for node in deduped.nodes] == [requested[2], *requested[:2], *requested[3:]]
    assert len(deduped.nodes) == len(requested)


def test_only_fully_contained_elements_survive_in_source_order():
    mesh = _build_mesh()
    result = mesh.subset_to_nodes(_requested(mesh))

    assert len(result.elements) == 4
    assert _source_index_connectivity(mesh, result) == _expected_connectivity(
        mesh, _REQUESTED_INDEX_ORDER
    )
    assert _WEB not in _source_index_connectivity(mesh, result)


def test_element_sets_are_restricted_and_empty_ones_dropped():
    mesh = _build_mesh()
    result = mesh.subset_to_nodes(_requested(mesh))

    assert "skin_and_web" in result.element_sets
    lookup = _coords_index(mesh)
    surviving = {
        tuple(lookup[tuple(node.coords.tolist())] for node in element.nodes)
        for element in result.element_sets["skin_and_web"].elements
    }
    assert surviving == {_QUADS[0], _QUADS[1]}
    assert "web_only" not in result.element_sets


def test_node_sets_are_restricted_and_empty_ones_dropped():
    mesh = _build_mesh()
    result = mesh.subset_to_nodes(_requested(mesh))

    assert "corner_and_web" in result.node_sets
    lookup = _coords_index(mesh)
    surviving = {
        lookup[tuple(node.coords.tolist())]
        for node in result.node_sets["corner_and_web"].nodes.values()
    }
    assert surviving == {0, 4}
    assert "web_nodes" not in result.node_sets


def test_deprecated_thickness_is_not_propagated():
    """The deprecated per-element thickness must not flood warnings into the view.

    The property map is the supported channel for shell thickness, and the
    ``thickness`` setter emits one ``DeprecationWarning`` per call, so the view
    deliberately drops the field the way :meth:`MeshModel.extract_submesh` does.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mesh = _build_mesh()
        mesh.elements[0].thickness = 0.1
        mesh.elements[2].thickness = 0.2
        mesh.elements[4].thickness = 0.3  # web element, excluded by the subset
        with warnings.catch_warnings():
            warnings.simplefilter("error", DeprecationWarning)
            result = mesh.subset_to_nodes(_requested(mesh))

    # Surviving elements are the four quads, in source order, all without thickness.
    assert [element.thickness for element in result.elements] == [None] * 4


def test_unknown_node_id_raises_value_error():
    mesh = _build_mesh()
    missing = [nid for nid in range(0, 1_000_000) if nid not in mesh.node_map][:3]
    with pytest.raises(ValueError) as excinfo:
        mesh.subset_to_nodes([mesh.nodes[0].id, *missing])

    message = str(excinfo.value)
    for nid in missing:
        assert str(nid) in message, f"missing id {nid} not named in: {message}"


def test_request_without_a_full_element_returns_nodes_only_mesh():
    mesh = _build_mesh()
    result = mesh.subset_to_nodes([mesh.nodes[0].id, mesh.nodes[1].id])

    assert result.node_count == 2
    assert result.elements_count == 0


def test_source_mesh_is_not_mutated():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        mesh = _build_mesh()
        mesh.elements[0].thickness = 0.1

        def snapshot():
            return (
                [node.id for node in mesh.nodes],
                [tuple(node.coords.tolist()) for node in mesh.nodes],
                [element.id for element in mesh.elements],
                [element.thickness for element in mesh.elements],
                {name: len(eset.elements) for name, eset in mesh.element_sets.items()},
                {name: len(nset.nodes) for name, nset in mesh.node_sets.items()},
            )

        before = snapshot()
        mesh.subset_to_nodes(_requested(mesh))
        after = snapshot()

    assert after == before
