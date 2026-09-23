"""Torsion + bending analytic benchmarks: thin-walled closed box (Bredt-Batho).

A rectangular hollow-section cantilever (four thin walls, isotropic) with
EXACT closed-form solutions:

  Torsion (Bredt-Batho, constant wall thickness t):
      GJ_exact = 4 A^2 / oint(ds/t) = 2 w^2 h^2 t / (w + h)
      twist rate: kappa = T / GJ

  Bending about x (flap, if x = the chord direction):
      I_x = w t h^2 / 2 + t h^3 / 6
      tip deflection under a tip load P: delta = P L^3 / (3 E I_x)

These are the absolute judges for the shell's torsional (membrane-shear)
and bending response — independent of any published beam table.  The
torsion case directly exercises the MITC4 membrane-shear path (the SRI
scheme), which is the open suspect for the GJ band vs the yaml's 6x6.

Reference magnitudes (w = h = 1.0 m, t = 0.01 m, L = 10 m, E = 70 GPa,
nu = 0.33):  G = 26.3 GPa,  GJ = 263.2e6 N m^2,  I = 6.667e-3 m^4.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve  # noqa: E402

DOF = 6

# Box geometry + material (SI)
W, H, L, T = 1.0, 1.0, 10.0, 0.01
E, NU = 70.0e9, 0.33
G = E / (2.0 * (1.0 + NU))

GJ_EXACT = G * 2.0 * W**2 * H**2 * T / (W + H)
I_X_EXACT = W * T * H**2 / 2.0 + T * H**3 / 6.0


def _box_mesh(nx: int, ny: int) -> MeshModel:
    """Four walls of a rectangular hollow section, corner nodes shared."""
    mesh = MeshModel()
    xs = np.linspace(-W / 2, W / 2, nx + 1)
    ys = np.linspace(-H / 2, H / 2, nx + 1)
    zs = np.linspace(0.0, L, ny + 1)
    grid: dict[tuple[float, float, float], Node] = {}

    def node(x, y, z):
        key = (round(x, 12), round(y, 12), round(z, 12))
        if key not in grid:
            grid[key] = Node([float(x), float(y), float(z)], geometric_node=False)
            mesh.add_node(grid[key])
        return grid[key]

    walls = [
        # (u-range, v-range, position function)
        (xs, [-H / 2], lambda u, v: (u, v, None)),   # bottom (y=-H/2)
        (xs, [+H / 2], lambda u, v: (u, v, None)),   # top    (y=+H/2)
        (ys, [-W / 2], lambda u, v: (v, u, None)),   # left   (x=-W/2)
        (ys, [+W / 2], lambda u, v: (v, u, None)),   # right  (x=+W/2)
    ]
    for uu, vv, fn in walls:
        for j in range(ny):
            for i in range(nx):
                n00 = node(*fn(uu[i], vv[0])[0:2], zs[j])
                n10 = node(*fn(uu[i + 1], vv[0])[0:2], zs[j])
                n11 = node(*fn(uu[i + 1], vv[0])[0:2], zs[j + 1])
                n01 = node(*fn(uu[i], vv[0])[0:2], zs[j + 1])
                mesh.add_element(
                    MeshElement(nodes=[n00, n10, n11, n01],
                                element_type=ElementType.quad))
    return mesh


def _assemble(mesh: MeshModel, rho: float = 0.0):
    from _aeroelast import PyMeshAssembler  # type: ignore[import-not-found]

    nodes_sorted = sorted(mesh.nodes, key=lambda n: n.id)
    node_id_to_idx = {n.id: i for i, n in enumerate(nodes_sorted)}
    ndof = len(nodes_sorted) * DOF
    mat = {"type": "isotropic", "e": E, "nu": NU, "rho": rho, "thickness": T,
           "shear_correction": 5.0 / 6.0}
    connectivity = [[node_id_to_idx[n.id] for n in e.nodes] for e in mesh.elements]
    asm = PyMeshAssembler(
        node_coords=np.asarray([n.coords[:3] for n in nodes_sorted], dtype=float),
        connectivity=connectivity,
        elem_types=[4] * len(connectivity),
        materials=[dict(mat) for _ in connectivity],
    )
    rows, cols, vals = asm.assemble_k()
    return coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr(), node_id_to_idx


def _solve(K, F, fixed):
    free = np.setdiff1d(np.arange(K.shape[0]), fixed)
    u = np.zeros(K.shape[0])
    u[free] = spsolve(K[free][:, free].tocsc(), F[free])
    return u


def _section_twist(mesh, m, u, z):
    """Best-fit rotation about the span (z) of the section's nodes at level z."""
    pts = np.array([n.coords for n in mesh.nodes if abs(n.coords[2] - z) < 1e-9])
    if len(pts) < 3:
        return np.nan
    c = pts.mean(axis=0)
    r = pts - c
    ux = np.array([u[m[_nid(mesh, p)][0]] if False else _ux(mesh, m, u, p) for p in pts])
    uy = np.array([_uy(mesh, m, u, p) for p in pts])
    num = (r[:, 0] * uy - r[:, 1] * ux).sum()
    den = (r[:, 0] ** 2 + r[:, 1] ** 2).sum()
    return num / den if den > 1e-30 else np.nan


def _nid(mesh, p):
    for n in mesh.nodes:
        if (abs(n.coords[0] - p[0]) < 1e-9 and abs(n.coords[1] - p[1]) < 1e-9
                and abs(n.coords[2] - p[2]) < 1e-9):
            return n.id
    raise KeyError(p)


def _ux(mesh, m, u, p):
    return u[m[_nid(mesh, p)] * DOF + 0]


def _uy(mesh, m, u, p):
    return u[m[_nid(mesh, p)] * DOF + 1]


def _find_node(mesh, x, y, z):
    return next(n for n in mesh.nodes
                if abs(n.coords[0] - x) < 1e-9 and abs(n.coords[1] - y) < 1e-9
                and abs(n.coords[2] - z) < 1e-9)


@pytest.fixture(scope="module")
def box_case():
    mesh = _box_mesh(nx=8, ny=40)
    K, m = _assemble(mesh)
    fixed = []
    for n in mesh.nodes:
        if n.coords[2] < 1e-12:
            fixed.extend(m[n.id] * DOF + d for d in range(DOF))
    return {"mesh": mesh, "K": K, "m": m, "fixed": np.array(fixed)}


def test_box_gj_matches_bredt(box_case):
    """Torsional stiffness vs the exact Bredt-Batho (closed thin-walled box)."""
    mesh, K, m, fixed = box_case["mesh"], box_case["K"], box_case["m"], box_case["fixed"]

    T_torque = 1.0e6  # N m
    corners = [(-W / 2, -H / 2), (W / 2, -H / 2), (W / 2, H / 2), (-W / 2, H / 2)]
    F = np.zeros(K.shape[0])
    for cx, cy in corners:
        nid = _find_node(mesh, cx, cy, L).id
        r = np.hypot(cx, cy)
        fx, fy = -cy / r, cx / r          # tangential unit vector
        fmag = T_torque / (4.0 * r)       # 4 corners, total torque = T
        F[m[nid] * DOF + 0] += fmag * fx
        F[m[nid] * DOF + 1] += fmag * fy

    u = _solve(K, F, fixed)

    zs = np.linspace(0.2 * L, 0.8 * L, 13)
    theta = np.array([_section_twist(mesh, m, u, z) for z in zs])
    kappa = np.polyfit(zs, theta, 1)[0]
    gj_num = T_torque / kappa

    rel = abs(gj_num - GJ_EXACT) / GJ_EXACT
    assert rel < 0.02, (
        f"GJ numeric {gj_num:.4e} vs exact {GJ_EXACT:.4e} N m^2 (rel {rel:.2%} > 2%)"
    )


def test_box_ei_matches_analytic(box_case):
    """Bending stiffness vs the analytic thin-wall section moment of inertia.

    Uses the mid-span curvature method (EI = M/kappa) — the tip deflection
    is contaminated by the local distortion under the discrete tip loads
    (the corners read 0.733 m vs 0.714 m analytic, +2.7%; the mid-width
    nodes bulge under the load).  The curvature estimate is immune to that.
    """
    mesh, K, m, fixed = box_case["mesh"], box_case["K"], box_case["m"], box_case["fixed"]

    P = 1.0e6  # tip load in +y (flap about x)
    F = np.zeros(K.shape[0])
    tip_ys = [-H / 2, H / 2]
    tip_xs = np.linspace(-W / 2, W / 2, 5)
    for y in tip_ys:
        for x in tip_xs:
            nid = _find_node(mesh, x, y, L).id
            F[m[nid] * DOF + 1] += P / (len(tip_ys) * len(tip_xs))

    u = _solve(K, F, fixed)
    zs = np.linspace(0.2 * L, 0.8 * L, 13)
    uys = np.array([u[m[_find_node(mesh, 0.0, H / 2, z).id] * DOF + 1] for z in zs])
    p = np.polyfit(zs, uys, 3)
    kappa = abs(np.polyval(np.polyder(np.polyder(p)), zs.mean()))
    M = P * (L - zs.mean())
    ei_num = M / kappa

    rel = abs(ei_num - E * I_X_EXACT) / (E * I_X_EXACT)
    assert rel < 0.02, (
        f"EI numeric {ei_num:.4e} vs exact {E * I_X_EXACT:.4e} N m^2 (rel {rel:.2%} > 2%)"
    )


def test_box_ea_matches_analytic(box_case):
    """Axial stiffness vs the analytic E * (structural area) = E * 4 w t.

    The clean discriminator for the blade's EA extraction bias: without the
    pre-bend/twist, the static-response extraction must match exactly.
    """
    mesh, K, m, fixed = box_case["mesh"], box_case["K"], box_case["m"], box_case["fixed"]

    F_ax = 1.0e6
    F = np.zeros(K.shape[0])
    tip_nodes = [n for n in mesh.nodes if n.coords[2] >= L - 1e-9]
    for n in tip_nodes:
        F[m[n.id] * DOF + 2] += F_ax / len(tip_nodes)

    u = _solve(K, F, fixed)
    zs = np.linspace(0.2 * L, 0.8 * L, 13)
    uz = np.array([u[m[_find_node(mesh, 0.0, H / 2, z).id] * DOF + 2] for z in zs])
    eps = np.polyfit(zs, uz, 1)[0]
    ea_num = F_ax / eps

    ea_exact = E * (4.0 * W * T)
    rel = abs(ea_num - ea_exact) / ea_exact
    assert rel < 0.02, (
        f"EA numeric {ea_num:.4e} vs exact {ea_exact:.4e} N (rel {rel:.2%} > 2%)"
    )
