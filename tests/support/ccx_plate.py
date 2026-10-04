"""The 2x2 quad-plate fixture shared by the CCX id-scheme tests.

``build_mesh`` and ``write_deck`` define the model and the deck that both the
software contract tests (``tests/software/contracts/test_ccx_writer_ids.py``)
and the physical parity test (``tests/validation/parity/test_ccx_writer_ids.py``)
operate on. They live here so the two halves cannot drift onto different
meshes: the parity test only means something if CalculiX solves the same deck
the contract tests describe.

The entity id counters ``Node._id_counter`` and ``MeshElement._id_counter`` are
process-global and never reset by ``MeshModel`` itself, so ``build_mesh(offset)``
sets them explicitly to place the ids at a chosen origin.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.io.writers import write_ccx_mesh
from aeroelast.core.mesh.model import MeshModel

OFFSET = 1000  # arbitrary: ids must not be assumed 0-based

ISOTROPIC = {"type": "isotropic", "e": 2.1e11, "nu": 0.3, "rho": 7800.0, "thickness": 0.01}


def build_mesh(offset: int) -> MeshModel:
    """A 2x2 quad plate strip; ids start at ``offset`` when it is non-zero."""
    Node._id_counter = offset
    MeshElement._id_counter = offset

    mesh = MeshModel()
    grid: dict[tuple[int, int], Node] = {}
    for j in range(3):
        for i in range(3):
            node = Node([float(i), float(j), 0.0], geometric_node=False)
            mesh.add_node(node)
            grid[(i, j)] = node

    for j in range(2):
        for i in range(2):
            mesh.add_element(
                MeshElement(
                    nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                    element_type=ElementType.quad,
                )
            )

    clamped = {n for n in mesh.nodes if np.isclose(n.y, 0.0)}
    free_face = {n for n in mesh.nodes if np.isclose(n.y, 2.0)}
    mesh.add_node_set(NodeSet("clamped", clamped))
    mesh.add_node_set(NodeSet("free_face", free_face))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def write_deck(tmp_path: Path, mesh: MeshModel, stem: str) -> Path:
    """Write the CCX deck for ``mesh`` and return the ``.inp`` path."""
    inp = tmp_path / f"{stem}.inp"
    write_ccx_mesh(
        mesh,
        str(inp),
        properties={"plate": ISOTROPIC},
        boundary_nodeset="clamped",
        load_nodeset="free_face",
        load_vector=[0.0, 0.0, 100.0],
        solver_type="LinearStatic",
    )
    return inp
