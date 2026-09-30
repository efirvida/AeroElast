"""D-Tube benchmark vs beam theory (Camarena & Anderson 2025).

Reproduces the verification case of Camarena & Anderson, "A critical
verification of beam and shell models of wind turbine blades",
Composite Structures 360 (2025) 118999, Section 3.1.1:

    Prismatic D-section cantilever, L = 100 m.  Semicircular shell
    (outer radius R = 2.25 m = D/2) closed by a flat web on the
    diameter, wall thickness t = 80 mm, steel, uniform 63 kPa traction
    applied normally on the semicircular portion (non-follower).

The analytical anchor is the classical prismatic-cantilever solution
about the section centroid:

    A      = pi*R*t + 2*R*t
    y_c    = -(2R/pi) * (pi*R*t) / A          (centroid below the web)
    I_y0   = (pi/2)*R^3*t                     (about the web line)
    I_c    = I_y0 - A*y_c^2                   (parallel axis, SUBTRACTED)
    w      = p * 2R                           (resultant per unit span)
    delta  = w*L^4 / (8*E*I_c)

Transverse shear adds w*L^2/(2*G*A) ~ 0.02 m (~0.1%), inside the tolerance.
TIP WARNING: adding the parallel-axis term instead of subtracting it
(y_c is 0.875 m below the web) inflates I_c and produced a bogus 2.6x
"shell/beam" ratio when this benchmark was first transcribed.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh import ElementSet, MeshElement, MeshModel, Node, NodeSet
from aeroelast.core.mesh.entities import ElementType
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

R = 2.25
T = 0.08
L = 100.0
P_TRACTION = 63.0e3
E = 2.0e11
NU = 0.3
N_ARC = 32
N_WEB = 16
# Span resolution from the measured convergence of the tip deflection
# (beam anchor delta = 24.5254 m, the test bounds the relative gap at 1%):
#
#   n_arc/n_web/n_z = 32/16/100   23.3594 m   4.75%   <- the earlier value
#   n_arc/n_web/n_z = 64/32/100   23.3862 m   4.65%   (section alone: no help)
#   n_arc/n_web/n_z = 64/32/200   24.2549 m   1.10%   (span: converged)
#   n_arc/n_web/n_z = 96/48/200   24.2561 m   1.10%   (more section: no change)
#
# The 4.75% was span discretisation and is gone.  What is left is a converged
# 1.10% and it is left red on purpose: it is about ten times the transverse
# shear estimate in the module docstring (w*L^2/(2GA) ~ 0.1%), so the residual
# is an open question about the element's bending-shear response, not a mesh
# artefact and not something to absorb by widening the bound.
NZ = 200

A_SEC = np.pi * R * T + 2.0 * R * T
Y_C = -(2.0 * R / np.pi) * (np.pi * R * T) / A_SEC
I_Y0 = 0.5 * np.pi * R**3 * T
I_C = I_Y0 - A_SEC * Y_C**2
W_RES = P_TRACTION * 2.0 * R
DELTA_REF = W_RES * L**4 / (8.0 * E * I_C)

A_ENCLOSED = 0.5 * np.pi * R**2
PERIM_OVER_T = (np.pi * R + 2.0 * R) / T
GJ_EXACT = (E / (2.0 * (1.0 + NU))) * 4.0 * A_ENCLOSED**2 / PERIM_OVER_T
TORQUE = 1.0e8


def _build_dtube_mesh(n_arc: int = N_ARC, n_web: int = N_WEB,
                      n_z: int = NZ) -> tuple[MeshModel, list[list[Node | None]], np.ndarray]:
    """Semicircle (below y=0) + flat web, extruded over L.  Shared nodes."""
    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()

    arc_thetas = np.linspace(0.0, np.pi, n_arc + 1)
    web_x = np.linspace(-R, R, n_web + 1)
    z_layers = np.linspace(0.0, L, n_z + 1)

    arc_nodes: list[list[Node | None]] = [[None] * (n_arc + 1) for _ in range(n_z + 1)]
    web_nodes: list[list[Node]] = [[None] * (n_web + 1) for _ in range(n_z + 1)]  # type: ignore[list-item]
    for k, zz in enumerate(z_layers):
        for j, xx in enumerate(web_x):
            n = Node([float(xx), 0.0, float(zz)], geometric_node=False)
            mesh.add_node(n)
            web_nodes[k][j] = n
        for i, th in enumerate(arc_thetas):
            if i == 0:
                arc_nodes[k][i] = web_nodes[k][n_web]
                continue
            if i == n_arc:
                arc_nodes[k][i] = web_nodes[k][0]
                continue
            n = Node([float(R * np.cos(th)), float(-R * np.sin(th)), float(zz)],
                     geometric_node=False)
            mesh.add_node(n)
            arc_nodes[k][i] = n

    for k in range(n_z):
        for i in range(n_arc):
            mesh.add_element(MeshElement(
                nodes=[arc_nodes[k][i], arc_nodes[k][i + 1],
                       arc_nodes[k + 1][i + 1], arc_nodes[k + 1][i]],
                element_type=ElementType.quad))
        for j in range(n_web):
            mesh.add_element(MeshElement(
                nodes=[web_nodes[k][j], web_nodes[k][j + 1],
                       web_nodes[k + 1][j + 1], web_nodes[k + 1][j]],
                element_type=ElementType.quad))

    clamped = {n for n in mesh.nodes if np.isclose(n.coords[2], 0.0)}
    mesh.add_node_set(NodeSet("clamped", clamped))
    mesh.add_element_set(ElementSet("tube", set(mesh.elements)))
    return mesh, arc_nodes, arc_thetas


def _apply_traction(mesh: MeshModel, solver: StaticLinearSolver,
                    arc_nodes: list[list[Node | None]], arc_thetas: np.ndarray,
                    n_arc: int = N_ARC, n_z: int = NZ) -> float:
    """Inward normal traction on the semicircle; returns the y-resultant."""
    dpn = solver.domain.dofs_per_node
    dth = np.pi / n_arc
    dz = L / n_z
    dofs: list[int] = []
    vals: list[float] = []
    f_y_total = 0.0
    for k in range(n_z + 1):
        w_z = dz if 0 < k < n_z else 0.5 * dz
        for i in range(n_arc + 1):
            th = arc_thetas[i]
            nrm = np.array([-np.cos(th), np.sin(th), 0.0])
            ds = dth * R if 0 < i < n_arc else 0.5 * dth * R
            f = P_TRACTION * ds * w_z * nrm
            node = arc_nodes[k][i]
            assert node is not None
            base = mesh.node_id_to_index[node.id] * dpn
            dofs += [base, base + 1]
            vals += [float(f[0]), float(f[1])]
            f_y_total += float(f[1])
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return f_y_total


@pytest.fixture(scope="module")
def dtube_linear():
    mesh, arc_nodes, arc_thetas = _build_dtube_mesh()
    prop = ShellProperty(
        material=IsotropicMaterial(name="Steel", E=E, nu=NU, rho=7850.0),
        thickness=T,
    )
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"tube": prop},
            "span_direction": (0.0, 0.0, 1.0),
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamped_dofs = sorted(
        mesh.node_id_to_index[n.id] * dpn + d
        for n in mesh.get_node_set("clamped").nodes.values()
        for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(clamped_dofs, 0.0)])
    f_y_total = _apply_traction(mesh, solver, arc_nodes, arc_thetas)
    u = np.asarray(solver.solve(), dtype=np.float64)
    tip_arc = [n for n in arc_nodes[NZ] if n is not None]
    tip_uy = float(np.mean([u[mesh.node_id_to_index[n.id] * dpn + 1] for n in tip_arc]))
    return {"u": u, "tip_uy": tip_uy, "f_y_total": f_y_total, "dpn": dpn, "mesh": mesh}


def test_applied_traction_resultant(dtube_linear):
    """The discretized traction integrates to w*L = p*D*L (vertical)."""
    assert dtube_linear["f_y_total"] == pytest.approx(W_RES * L, rel=1e-3)


def test_dtube_tip_deflection_matches_beam_theory(dtube_linear):
    """Linear MITC4 shell vs the prismatic-cantilever analytical deflection.

    The linear shell must reproduce beam theory (bending about the
    centroidal axis) well inside 1% for a prismatic closed section.
    """
    tip_uy = dtube_linear["tip_uy"]
    rel = abs(tip_uy - DELTA_REF) / DELTA_REF
    assert rel < 0.01, (
        f"D-Tube tip {tip_uy:.3f} m vs beam theory {DELTA_REF:.3f} m "
        f"(I_c={I_C:.4f} m^4, rel {rel:.2%} > 1%)"
    )


def _find_node(mesh: MeshModel, x: float, y: float, z: float) -> Node:
    return next(
        n for n in mesh.nodes
        if abs(n.coords[0] - x) < 1e-9 and abs(n.coords[1] - y) < 1e-9
        and abs(n.coords[2] - z) < 1e-9
    )


@pytest.fixture(scope="module")
def dtube_torsion():
    """Tip torque as a statically equivalent Bredt shear-flow load.

    Counterclockwise tangential forces along the wall (weights = trapezoid
    segment lengths, direction via the local wall tangent), scaled so the
    net moment about z is exactly TORQUE.  No nodal moments on the
    drilling DOF are used.

    Circumferential resolution is doubled with respect to the bending mesh:
    the torsional error is controlled by the wall discretization
    (32 arc elements -> 2.10%, 64 -> 0.92%, span refinement has no effect).
    """
    mesh, _arc_nodes, _thetas = _build_dtube_mesh(n_arc=64, n_web=32)
    prop = ShellProperty(
        material=IsotropicMaterial(name="Steel", E=E, nu=NU, rho=7850.0),
        thickness=T,
    )
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"tube": prop},
            "span_direction": (0.0, 0.0, 1.0),
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamped_dofs = sorted(
        mesh.node_id_to_index[n.id] * dpn + d
        for n in mesh.get_node_set("clamped").nodes.values()
        for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(clamped_dofs, 0.0)])

    dth = np.pi / N_ARC
    dweb = 2.0 * R / N_WEB
    items: list[tuple[Node, float, float, float]] = []
    torque_unit = 0.0
    for node in (n for n in mesh.nodes if abs(n.coords[2] - L) < 1e-9):
        x, y = node.coords[0], node.coords[1]
        if y < -1e-9:
            w = dth * R
            tx, ty = -y / R, x / R
        else:
            w = dweb
            tx, ty = -1.0, 0.0
        if abs(y) < 1e-9 and abs(abs(x) - R) < 1e-9:
            w = 0.5 * dth * R + 0.5 * dweb
        items.append((node, w, tx, ty))
        torque_unit += w * (x * ty - y * tx)

    scale = TORQUE / torque_unit
    dofs: list[int] = []
    vals: list[float] = []
    for node, w, tx, ty in items:
        base = mesh.node_id_to_index[node.id] * dpn
        dofs += [base, base + 1]
        vals += [scale * w * tx, scale * w * ty]
    solver.add_nodal_loads([NodalLoad(dofs, vals)])

    u = np.asarray(solver.solve(), dtype=np.float64)
    zs = np.linspace(0.2 * L, 0.8 * L, 13)
    theta = []
    for z in zs:
        plus = _find_node(mesh, R, 0.0, z)
        minus = _find_node(mesh, -R, 0.0, z)
        theta.append(
            (u[mesh.node_id_to_index[plus.id] * dpn + 1]
             - u[mesh.node_id_to_index[minus.id] * dpn + 1]) / (2.0 * R)
        )
    kappa = float(np.polyfit(zs, theta, 1)[0])
    return {"gj_num": TORQUE / kappa, "torque_unit": torque_unit}


def test_dtube_torsion_matches_bredt(dtube_torsion):
    """Linear MITC4 shell vs Bredt-Batho closed-cell torsional stiffness.

    Single-cell D-section: GJ = G * 4A^2 / oint(ds/t) with A = pi*R^2/2.
    The twist rate is measured in the St. Venant window [0.2, 0.8] L from
    the relative transverse displacement of the web ends.  With the 64
    arc-element mesh the error is 0.92% (converged in the span).
    """
    gj_num = dtube_torsion["gj_num"]
    rel = abs(gj_num - GJ_EXACT) / GJ_EXACT
    assert rel < 0.02, (
        f"D-Tube GJ {gj_num:.4e} vs Bredt {GJ_EXACT:.4e} N m^2 "
        f"(rel {rel:.2%} > 2%)"
    )
