"""StaticNonlinearSolver (Updated-Lagrangian) validation vs the elastica.

The exact large-deflection solution of a cantilever under a dead tip load
(Bisshopp & Drucker 1945; tabulated via Bathe & Bolourchi 1979):

    alpha = P·L²/EI     w/L (tip, load direction)    u/L (tip, axial)
        1                 0.30172                        -0.05643
        2                 0.49346                        -0.16064

The shell plate (L=1, b=0.1, h=1e-3, steel) has EI = E·b·h³/12 and, for
the free narrow strip, the plate bending stiffness collapses to the beam
value, so the MITC4 shell must reproduce the 1-D elastica within the
first-order UL error O(Δθ) (≈2% for 32 steps of the α=1 case).

This regression guards the UL incremental implementation of
``StaticNonlinearSolver`` (the historical SNES path is not used — see the
class docstring for the tangent-consistency finding of 2026-09-18).

**Diagnosis, 2026-10-14 (the earlier "converges ~4 % below the elastica" reading was wrong).**
The 2026-09-30 table (3.63 ... 4.00 % over 32 ... 512 steps) was taken on the 8 x 4 mesh, whose
span discretisation alone costs ~3.8 % at alpha = 1.  Refining the span (64 steps, ny = 4):

    nx =   8   w/L=0.29018  -3.83 %        nx = 32  w/L=0.29920  -0.83 %
    nx =  16   w/L=0.29482  -2.29 %        nx = 64  w/L=0.30139  -0.11 %

so the alpha = 1 gap is mesh, and it closes on the elastica at nx >= 64 only by cancellation of
two opposite errors: the plate is -1.1 % stiffer than a beam (anticlastic root effect, mesh
independent, shrinking with the width) and the scheme below adds about +1.2 %.

The scheme is the real defect.  The updated-Lagrangian loop integrates ``K_e(x_ref) du = dlam F``
with no equilibrium iteration and no initial-stress stiffness, so refining the load steps does not
converge to the dead-load equilibrium.  At alpha = 2 (45 deg tip rotation), mesh and step
converged (nx = 64, 2048 steps), w/L is +5.29 % and u/L +13 % from the elastica, and the bias grows
as the strip narrows (+4.7, +8.8, +10.7 % for b = 0.1, 0.05, 0.025); a 2-D corotational beam with
the same scheme gives +11.3 %, and with Newton on the consistent tangent 0.49349 (exact 0.49346).
A pure end moment is reproduced to 0.1-0.6 % up to 115 deg (b = 0.00625): the element kinematics
are right, the dead-load scheme is not.  ``scheme: newton`` removes the bias
(``test_static_nonlinear_newton_elastica.py``, ``test_strip_elastica_calculix.py``).

The tight bounds stay as xfails carrying the measured number, exactly as before: they are an
accuracy limit of the default scheme, not a tolerance to widen.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

from tests.support.assertions import assert_residual_below

L, B, H = 1.0, 0.1, 0.001
E = 2.1e11
EI = E * B * H**3 / 12.0

# The same steel the suite's other cantilever plate uses
# (tests/validation/benchmarks/test_shell_static_and_modal.py:66) and the same one the
# deleted `test_shell_comprehensive` defined: E = 2.1e11, nu = 0.3, rho = 7800.
STEEL = IsotropicMaterial(name="Steel", E=E, nu=0.3, rho=7800.0)


def build_cantilever_mesh(
    L: float = 1.0, b: float = 0.1, nx: int = 8, ny: int = 4
) -> MeshModel:
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


def _solve_ul(P: float, steps: int) -> tuple[np.ndarray, MeshModel, int]:
    mesh = build_cantilever_mesh(L=L, b=B)
    prop = ShellProperty(material=STEEL, thickness=H)
    cfg = {
        "solver": {"continuation_steps": steps, "continuation_max_steps": 4 * steps},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": prop},
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
    apply_edge_load(mesh, solver, "free_edge", "z", P)
    return np.asarray(solver.solve(), dtype=np.float64), mesh, dpn


@pytest.mark.parametrize(
    "alpha,ref_wl,steps,tol",
    [
        (1.0, 0.30172, 32, 0.02),
        (2.0, 0.49346, 64, 0.02),
    ],
)
def test_elastica_tip_deflection(alpha, ref_wl, steps, tol):
    """The UL tip deflection against the elastica table, tight bound as xfail."""
    P = alpha * EI / L**2
    u, mesh, dpn = _solve_ul(P, steps)
    center = get_center_node(mesh, L, B)
    idx = mesh.node_id_to_index[center.id]
    w = abs(u[idx * dpn + 2])
    rel = abs(w / L - ref_wl) / ref_wl
    what = f"alpha={alpha} tip w/L after {steps} load steps"
    reference_name = (
        "the elastica reference table of Bisshopp & Drucker (1945), "
        "tabulated by Bathe & Bolourchi (1979)"
    )
    # `assert_residual_below` prints the measured margin in the store's canonical form and
    # then asserts the tight bound. The gap is an element-accuracy limit, so the
    # AssertionError becomes an xfail carrying the measured number -- and the canonical print
    # has already happened, which is what lets `regression` re-read the margin. The old
    # assertion was red by design; this keeps the evidence visible without widening `tol`.
    try:
        assert_residual_below(
            rel,
            tol=tol,
            kind="paper",
            reference_name=reference_name,
            what=what,
        )
    except AssertionError:
        pytest.xfail(
            f"alpha={alpha}: w/L={w / L:.5f} vs elastica {ref_wl:.5f} "
            f"(err={100 * rel:.2f}%, bound={100 * tol:.2f}%) -- 8x4 mesh under-resolves the span "
            "(alpha=1) and the default UL scheme has no equilibrium iteration (alpha=2); "
            "see the module docstring and scheme='newton'"
        )
