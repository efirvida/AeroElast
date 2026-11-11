"""Geometry-only public windIO 2 adapter; no centered shell properties.

Object-keyed records survive MeshModel renumbering. The immutable source and
raw records retain physical reference offsets, not an FE approximation.
"""
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel
from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh


@dataclass(frozen=True)
class ElementMetadata:
    source_index: int
    parameters: object
    basis: object
    cell: object
    stack: object


@dataclass(frozen=True)
class WindIO2Metadata:
    source_sha256: str
    source: object
    raw_mesh: object
    elements: object
    attachments: object
    structural_ready: bool = False


def geometry_mesh(topology, element_size, eta_bounds, renumber=None, min_strip_m=1e-6):
    """Adapt without welding, dropping cells, centering or material averaging."""
    if renumber not in (None, 'simple', 'rcm'):
        raise ValueError(f'Unknown renumbering algorithm: {renumber}')
    # Complete every fallible raw geometry/material check before entity creation.
    raw = generate_mesh(topology, element_size=element_size, eta_bounds=eta_bounds,
                        with_materials=True, min_strip_m=min_strip_m)
    mesh = MeshModel()
    for xyz in raw.nodes:
        mesh.add_node(Node(xyz, coordinate_dtype=np.longdouble))
    records, surfaces = {}, {}
    for i, connectivity in enumerate(raw.elements):
        live = [int(j) for j in connectivity if j >= 0]
        element = MeshElement([mesh.nodes[j] for j in live],
                              ElementType.quad if len(live) == 4 else ElementType.triangle)
        mesh.add_element(element)
        records[element] = ElementMetadata(i, raw.element_parameters[i][:len(live)], raw.bases[i],
                                          raw.cells[i], raw.material_stacks[i])
        surfaces.setdefault(raw.surfaces[i], set()).add(element)
    for name, elements in surfaces.items():
        mesh.add_element_set(ElementSet(name=name, elements=elements))
    attachments = {key: tuple(mesh.nodes[i] for i in chain)
                   for key, chain in raw.attachments.items()}
    for (name, handle), nodes in attachments.items():
        mesh.add_node_set(NodeSet(name=f'{name}:{handle}', nodes=set(nodes)))
    mesh.windio2_metadata = WindIO2Metadata(
        topology.source_sha256, topology.source, raw,
        MappingProxyType(records), MappingProxyType(attachments))
    if renumber is not None:
        mesh.renumber_mesh(algorithm=renumber)
    return mesh
