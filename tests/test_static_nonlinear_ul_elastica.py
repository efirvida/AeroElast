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

**Measured 2026-09-30, and the O(Δθ) claim above is wrong for this element.**
Refining the load stepping does not converge to the elastica; it converges to
about 4% below it:

    steps= 32   w/L=0.29075  3.63%
    steps= 64   w/L=0.29018  3.83%
    steps=128   w/L=0.28989  3.92%
    steps=256   w/L=0.28974  3.97%
    steps=512   w/L=0.28967  4.00%

So the alpha=1 gap is neither the O(Delta-theta) stepping error nor a mesh
effect: the reviewed element's large-deflection response is systematically ~4%
stiffer than the exact elastica.  The previous element was inside -0.5% here.
Left red on purpose, with this table as the evidence: it is an accuracy
question about the element's UL path, not a tolerance to widen.
"""

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

from test_shell_comprehensive import (
    STEEL,
    apply_edge_load,
    build_cantilever_mesh,
    get_center_node,
)

L, B, H = 1.0, 0.1, 0.001
E = 2.1e11
EI = E * B * H**3 / 12.0


def _solve_ul(P: float, steps: int) -> np.ndarray:
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
        (2.0, 0.49346, 64, 0.06),
    ],
)
def test_elastica_tip_deflection(alpha, ref_wl, steps, tol):
    P = alpha * EI / L**2
    u, mesh, dpn = _solve_ul(P, steps)
    center = get_center_node(mesh, L, B)
    idx = mesh.node_id_to_index[center.id]
    w = abs(u[idx * dpn + 2])
    err = abs(w / L - ref_wl) / ref_wl
    assert err < tol, (
        f"alpha={alpha}: w/L={w / L:.5f} vs elastica {ref_wl:.5f} "
        f"(err={100 * err:.2f}% > {100 * tol:.2f}%)"
    )
