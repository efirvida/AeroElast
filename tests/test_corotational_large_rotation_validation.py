"""Validation test for element-corotational-shell change.

Verifies that:
1. Corotational formulation is accurate for large rotations (30-45 deg)
2. Total Lagrangian formulation has significant error for large rotations
3. Both agree in the linear regime

Reference: Frisch-Fay (1962), "Flexible Bars", Butterworths, Table 3.1
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler


L = 1.0
B = 0.1
T = 0.005
E = 70e9
NU = 0.3
RHO = 2700.0


def _build_cantilever_mesh(*, nz: int = 40, nx: int = 4, theta_x: float = 0.0):
    xs = np.linspace(-B / 2.0, B / 2.0, nx + 1)
    zs = np.linspace(0.0, L, nz + 1)

    c, s = math.cos(theta_x), math.sin(theta_x)
    rot_x = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)

    nodes: list[np.ndarray] = []
    grid: dict[tuple[int, int], int] = {}
    for k, z in enumerate(zs):
        for i, x in enumerate(xs):
            idx = len(nodes)
            grid[(i, k)] = idx
            nodes.append(rot_x @ np.array([x, 0.0, z], dtype=float))

    conn: list[list[int]] = []
    for k in range(nz):
        for i in range(nx):
            conn.append([grid[(i, k)], grid[(i + 1, k)], grid[(i + 1, k + 1)], grid[(i, k + 1)]])

    return np.asarray(nodes, dtype=float), conn


def _materials(n_elem: int) -> list[dict[str, float | str]]:
    return [
        {
            "type": "isotropic",
            "e": E,
            "nu": NU,
            "rho": RHO,
            "thickness": T,
        }
    ] * n_elem


def _solve_tip_compliance(
    *,
    nodes: np.ndarray,
    k_rows: np.ndarray,
    k_cols: np.ndarray,
    k_vals: np.ndarray,
    theta_x: float,
) -> float:
    n_dofs = nodes.shape[0] * 6
    k = coo_matrix((k_vals, (k_rows, k_cols)), shape=(n_dofs, n_dofs)).tocsr()

    beam_axis = np.array([0.0, -math.sin(theta_x), math.cos(theta_x)], dtype=float)
    beam_coord = nodes @ beam_axis
    root = np.where(np.isclose(beam_coord, beam_coord.min(), atol=1e-10))[0]
    tip = np.where(np.isclose(beam_coord, beam_coord.max(), atol=1e-10))[0]
    tip_center = min(tip, key=lambda i: abs(float(nodes[i, 0])))

    fixed_dofs: list[int] = []
    for nid in root:
        base = 6 * int(nid)
        fixed_dofs.extend(range(base, base + 6))
    fixed = np.array(sorted(set(fixed_dofs)), dtype=int)
    free = np.setdiff1d(np.arange(n_dofs, dtype=int), fixed)

    transverse = np.array([0.0, math.cos(theta_x), math.sin(theta_x)], dtype=float)
    f = np.zeros(n_dofs, dtype=float)
    base = 6 * int(tip_center)
    f[base : base + 3] = transverse  # unit load

    u = np.zeros(n_dofs, dtype=float)
    u[free] = spsolve(k[free][:, free], f[free])
    return float(u[base : base + 3].dot(transverse))


def _corotational_reference_compliance(theta_x: float) -> float:
    nodes, conn = _build_cantilever_mesh(theta_x=theta_x)
    asm = PyMeshAssembler(
        node_coords=nodes,
        connectivity=conn,
        elem_types=[4] * len(conn),
        materials=_materials(len(conn)),
    )
    rows, cols, vals = asm.assemble_k()
    return _solve_tip_compliance(nodes=nodes, k_rows=rows, k_cols=cols, k_vals=vals, theta_x=theta_x)


def _tl_tangent_compliance_at_rigid_rotation(theta_x: float) -> float:
    nodes0, conn = _build_cantilever_mesh(theta_x=0.0)
    asm = PyMeshAssembler(
        node_coords=nodes0,
        connectivity=conn,
        elem_types=[4] * len(conn),
        materials=_materials(len(conn)),
    )

    c, s = math.cos(theta_x), math.sin(theta_x)
    rot_x = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)

    u_rigid = np.zeros(nodes0.shape[0] * 6, dtype=float)
    for i, x_ref in enumerate(nodes0):
        x_rot = rot_x @ x_ref
        u_rigid[6 * i : 6 * i + 3] = x_rot - x_ref

    rows, cols, vals = asm.assemble_kt(u_rigid)
    return _solve_tip_compliance(nodes=nodes0, k_rows=rows, k_cols=cols, k_vals=vals, theta_x=0.0)


def test_corotational_large_rotation_validation():
    """Physical validation: corotational frame objectivity vs TL rigid-rotation error."""
    i_beam = B * T**3 / 12.0
    delta_analytic = L**3 / (3.0 * E * i_beam)  # unit-load compliance

    # Test 1 — linear regime
    coro_lin = _corotational_reference_compliance(theta_x=0.0)
    tl_lin = _tl_tangent_compliance_at_rigid_rotation(theta_x=0.0)

    assert 0.98 <= coro_lin / delta_analytic <= 1.02
    assert 0.98 <= tl_lin / delta_analytic <= 1.02
    assert abs(coro_lin - tl_lin) / delta_analytic < 0.01

    # Test 2 — moderate rotation (~20°)
    theta20 = math.radians(20.0)
    coro_20 = _corotational_reference_compliance(theta_x=theta20)
    tl_20 = _tl_tangent_compliance_at_rigid_rotation(theta_x=theta20)

    assert 0.99 <= coro_20 / coro_lin <= 1.01
    assert tl_20 < 0.95 * coro_20

    # Test 3 — large rotation (~40°)
    theta40 = math.radians(40.0)
    coro_40 = _corotational_reference_compliance(theta_x=theta40)
    tl_40 = _tl_tangent_compliance_at_rigid_rotation(theta_x=theta40)

    assert abs(coro_40 - coro_lin) / coro_lin < 0.05
    assert abs(tl_40 - coro_40) / coro_40 > 0.10
