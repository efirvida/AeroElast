"""Multi-cell closed-section torsion vs Bredt-Batho (isotropic box).

Verification case for the MITC4+ torsional path on a *multi-cell* section:
a rectangular hollow box (1 x 1 x 10, t = 0.01, steel-like isotropic) with a
central web splitting it into two cells.

Analytic reference (multi-cell Bredt-Batho shear-flow system):
    cell i:  q_i*(p_i/t) + (q_i - q_j)*(h/t_m) = 2*A_i*G*theta'
    torque:  T = 2*sum(q_i*A_i)      ->  J = T/(G*theta')
For the symmetric two-cell box the shared web carries zero NET flow, so the
result equals the single-cell value J = 2*w^2*h^2*t/(w+h) = 1.0e-2 m^4.

LOADING CONVENTION (important): the tip torque must be applied as *corner
force couples*, not as tangential forces spread over every tip-section node.
A distributed perimeter couple is not a physical load for a thin-walled
section and excites the very soft in-plane shear (parallelogram) mode of the
box (wall bending ~ E t^3 vs membrane torsion ~ G t), which collapses the
apparent J by ~40 %.  The corner couple (the same convention as
``test_box_torsion_bending_benchmark.py``) reproduces Bredt to <0.1 %.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import ElementSet, MeshModel
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

W, H, L, T = 1.0, 1.0, 10.0, 0.01
E, NU = 70.0e9, 0.33
G = E / (2.0 * (1.0 + NU))
J_EXACT = 2.0 * W**2 * H**2 * T / (W + H)  # single cell == two-cell (web inert)
NX, NY = 8, 40
TORQUE = 1.0e6


def _box_mesh(with_web: bool) -> MeshModel:
    mesh = MeshModel()
    xs = np.linspace(-W / 2, W / 2, NX + 1)
    ys = np.linspace(-H / 2, H / 2, NX + 1)
    zs = np.linspace(0.0, L, NY + 1)
    grid: dict[tuple[float, float, float], Node] = {}

    def node(x, y, z):
        key = (round(x, 12), round(y, 12), round(z, 12))
        if key not in grid:
            grid[key] = Node([float(x), float(y), float(z)], geometric_node=False)
            mesh.add_node(grid[key])
        return grid[key]

    def wall(xs_w, ys_w):
        for j in range(NY):
            for i in range(NX):
                n00 = node(xs_w[i], ys_w[i], zs[j])
                n10 = node(xs_w[i + 1], ys_w[i + 1], zs[j])
                n11 = node(xs_w[i + 1], ys_w[i + 1], zs[j + 1])
                n01 = node(xs_w[i], ys_w[i], zs[j + 1])
                mesh.add_element(MeshElement(nodes=[n00, n10, n11, n01],
                                             element_type=ElementType.quad))

    for y in (-H / 2, H / 2):
        wall(xs, np.full(NX + 1, y))
    for x in (-W / 2, W / 2):
        wall(np.full(NX + 1, x), ys)
    if with_web:
        wall(np.zeros(NX + 1), ys)
    mesh.add_element_set(ElementSet("box", set(mesh.elements)))
    return mesh


def _gj(mesh: MeshModel) -> float:
    prop = ShellProperty(
        material=IsotropicMaterial(name="steel", E=E, nu=NU, rho=0.0),
        thickness=T,
    )
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": {"box": prop},
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    nid2i = mesh.node_id_to_index
    coords = mesh.coords_array
    root = [n for n in mesh.nodes if abs(n.coords[2]) < 1e-12]
    solver.add_dirichlet_conditions([DirichletCondition(
        sorted(nid2i[n.id] * dpn + d for n in root for d in range(dpn)), 0.0)])

    # corner force couples (physical convention for a thin-walled section)
    corners = [(-W / 2, -H / 2), (W / 2, -H / 2), (W / 2, H / 2), (-W / 2, H / 2)]
    dofs, vals = [], []
    for cx, cy in corners:
        n = next(nd for nd in mesh.nodes
                 if abs(nd.coords[0] - cx) < 1e-9 and abs(nd.coords[1] - cy) < 1e-9
                 and abs(nd.coords[2] - L) < 1e-9)
        r = np.hypot(cx, cy)
        f = TORQUE / (4.0 * r)
        dofs += [nid2i[n.id] * dpn, nid2i[n.id] * dpn + 1]
        vals += [f * (-cy / r), f * (cx / r)]
    solver.add_nodal_loads([NodalLoad(dofs, vals)])

    u = np.asarray(solver.solve(), dtype=np.float64)
    zs = np.linspace(0.2 * L, 0.8 * L, 13)
    theta = []
    for z in zs:
        idx = np.nonzero(np.abs(coords[:, 2] - z) < 1e-9)[0]
        p = coords[idx][:, :2] - coords[idx][:, :2].mean(axis=0)
        ux = u[idx * dpn + 0]
        uy = u[idx * dpn + 1]
        theta.append(float((p[:, 0] * uy - p[:, 1] * ux).sum() / (p**2).sum()))
    kappa = abs(np.polyfit(zs, theta, 1)[0])
    return TORQUE / kappa


@pytest.mark.parametrize("with_web", [False, True])
def test_closed_section_torsion_matches_bredt(with_web):
    """Single- and two-cell boxes both reproduce J = 1.0e-2 m^4 within 2 %."""
    gj = _gj(_box_mesh(with_web))
    rel = abs(gj - G * J_EXACT) / (G * J_EXACT)
    label = "two-cell (central web)" if with_web else "single-cell"
    assert rel < 0.02, (
        f"{label}: GJ={gj:.4e} vs Bredt {G * J_EXACT:.4e} N m^2 "
        f"(rel {rel:.2%} > 2%)"
    )
