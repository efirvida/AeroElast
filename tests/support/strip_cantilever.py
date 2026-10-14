"""Shared helpers for the dead-tip-loaded shell cantilever strip (elastica and CalculiX pairs)."""

from __future__ import annotations

import numpy as np

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

L, B, H = 1.0, 0.1, 0.001
E = 2.1e11
EI = E * B * H**3 / 12.0

# The same steel the suite's other cantilever plate uses
# (tests/validation/benchmarks/test_shell_static_and_modal.py:66) and the same one the
# deleted `test_shell_comprehensive` defined: E = 2.1e11, nu = 0.3, rho = 7800.
STEEL = IsotropicMaterial(name="Steel", E=E, nu=0.3, rho=7800.0)


def build_cantilever_mesh(L: float = 1.0, b: float = 0.1, nx: int = 8, ny: int = 4) -> MeshModel:
    """The rectangular MITC4 cantilever strip the deleted helper built.

    Reproduced with the old dimensions (``nx=8``, ``ny=4``, ``L=1.0``, ``b=0.1``) so the
    recorded elastica table stays comparable: a ``clamped`` node set at ``x=0`` and a
    ``free_edge`` node set at ``x=L``.
    """
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, b, ny + 1)

    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            node = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(node)
            grid[(i, j)] = node

    for j in range(ny):
        for i in range(nx):
            mesh.add_element(
                MeshElement(
                    nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                    element_type=ElementType.quad,
                )
            )

    clamped = {n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)}
    free = {n for n in mesh.nodes if np.isclose(n.x, L, atol=1e-12)}
    mesh.add_node_set(NodeSet("clamped", clamped))
    mesh.add_node_set(NodeSet("free_edge", free))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def get_center_node(mesh: MeshModel, L: float, b: float) -> Node:
    """The node at the centre of the free edge."""
    target_y = 0.5 * b
    candidates = [n for n in mesh.nodes if np.isclose(n.x, L, atol=1e-9)]
    return min(candidates, key=lambda n: abs(float(n.y) - target_y))


def apply_edge_load(
    mesh: MeshModel, solver, nodeset: str, direction: str, total_load: float
) -> None:
    """Distribute the total load ``P`` over the nodes of ``nodeset`` in ``direction``."""
    nodes = sorted(mesh.get_node_set(nodeset).nodes.values(), key=lambda n: n.id)
    load_per_node = total_load / len(nodes)
    comp = {"x": 0, "y": 1, "z": 2}[direction]

    for node in nodes:
        idx = mesh.node_id_to_index[node.id]
        dpn = solver.domain.dofs_per_node
        force = [0.0] * dpn
        force[comp] = load_per_node
        solver.add_nodal_loads([NodalLoad(list(range(idx * dpn, idx * dpn + dpn)), force)])


ELASTICA = {1.0: (0.30172, -0.05643), 2.0: (0.49346, -0.16064)}
ELASTICA_REFERENCE = (
    "the elastica reference table of Bisshopp & Drucker (1945), "
    "tabulated by Bathe & Bolourchi (1979)"
)


def solve_strip(alpha: float, *, b: float, nx: int, ny: int, steps: int, scheme: str | None):
    """Tip (w/L, u/L) of the dead-tip-loaded strip of width ``b``."""
    solver_cfg: dict = {"continuation_steps": steps, "continuation_max_steps": 4 * steps}
    if scheme is not None:
        solver_cfg["scheme"] = scheme
    mesh = build_cantilever_mesh(L=L, b=b, nx=nx, ny=ny)
    cfg = {
        "solver": solver_cfg,
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": ShellProperty(material=STEEL, thickness=H)},
            "span_direction": (1.0, 0.0, 0.0),
        },
    }
    solver = StaticNonlinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamped = []
    for node in mesh.get_node_set("clamped").nodes.values():
        i0 = mesh.node_id_to_index[node.id] * dpn
        clamped.extend(range(i0, i0 + dpn))
    solver.add_dirichlet_conditions([DirichletCondition(clamped, 0.0)])
    ei = E * b * H**3 / 12.0
    apply_edge_load(mesh, solver, "free_edge", "z", alpha * ei / L**2)
    u = np.asarray(solver.solve(), dtype=np.float64)
    idx = mesh.node_id_to_index[get_center_node(mesh, L, b).id]
    return abs(u[idx * dpn + 2]) / L, u[idx * dpn] / L
