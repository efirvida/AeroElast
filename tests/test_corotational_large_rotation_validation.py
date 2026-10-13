"""Validation tests for corotational tangent frame objectivity.

These tests compare total-Lagrangian tangent (assemble_kt) against
corotational tangent (assemble_kt_corotational) directly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.sparse import coo_matrix

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler


L = 1.0
B = 0.1
T = 0.005
E = 70e9
NU = 0.3
RHO = 2700.0


def _build_cantilever_mesh(*, nx: int = 4, ny: int = 1):
    """4 MITC4 shell elements (5x2 nodes) in XY plane."""
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(-B / 2.0, B / 2.0, ny + 1)

    nodes: list[np.ndarray] = []
    grid: dict[tuple[int, int], int] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            idx = len(nodes)
            grid[(i, j)] = idx
            nodes.append(np.array([x, y, 0.0], dtype=float))

    conn: list[list[int]] = []
    for j in range(ny):
        for i in range(nx):
            conn.append([grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]])

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


def _build_assembler() -> tuple[np.ndarray, PyMeshAssembler]:
    nodes, conn = _build_cantilever_mesh()
    asm = PyMeshAssembler(
        node_coords=nodes,
        connectivity=conn,
        elem_types=[4] * len(conn),
        materials=_materials(len(conn)),
    )
    return nodes, asm


def _coo_to_dense(rows: np.ndarray, cols: np.ndarray, vals: np.ndarray, n: int) -> np.ndarray:
    return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr().toarray()


def _translation_dof_indices(n_nodes: int) -> np.ndarray:
    idx: list[int] = []
    for i in range(n_nodes):
        base = 6 * i
        idx.extend([base, base + 1, base + 2])
    return np.asarray(idx, dtype=int)


def _rigid_rotation_displacement(nodes: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    c, s = math.cos(theta), math.sin(theta)
    r = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)

    u_rigid = np.zeros(nodes.shape[0] * 6, dtype=float)
    for i, x_ref in enumerate(nodes):
        x_rot = r @ x_ref
        u_rigid[6 * i : 6 * i + 3] = x_rot - x_ref
    return u_rigid, r


def _transform_translation_submatrix(k_tt: np.ndarray, r: np.ndarray, n_nodes: int) -> np.ndarray:
    t = np.kron(np.eye(n_nodes), r)
    return t @ k_tt @ t.T


def _rel_err(a: np.ndarray, b: np.ndarray, ref: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(ref))



def test_kt_coro_equals_kt_tl_at_zero_displacement():
    nodes, asm = _build_assembler()
    u0 = np.zeros(asm.dofs_count, dtype=float)
    rows_tl, cols_tl, vals_tl = asm.assemble_kt(u0)
    rows_c, cols_c, vals_c = asm.assemble_kt_corotational(u0)

    k_tl = _coo_to_dense(rows_tl, cols_tl, vals_tl, asm.dofs_count)
    k_c = _coo_to_dense(rows_c, cols_c, vals_c, asm.dofs_count)
    err = _rel_err(k_c, k_tl, k_tl)
    assert err < 1e-10


# Round-off bound: the measured error is 1.7e-12..2.6e-12 for 1..30 deg (SLURM 11614644,
# tools/diagnose_frame_objectivity.py), so 1e-9 leaves three orders of margin and still
# fails for any frame-dependent tangent (errors of 1e-3 or more).
TOL_OBJECTIVITY = 1e-9


@pytest.mark.parametrize("assemble", ["assemble_kt", "assemble_kt_corotational"])
@pytest.mark.parametrize("angle_deg", [1.0, 5.0, 30.0])
def test_tangent_is_frame_objective_under_rigid_rotation(assemble, angle_deg):
    """K_t(u_rigid) == T K0 T^T on the translational block, for TL and corotational paths.

    Replaces the former ``..._corotational_is_frame_objective_tl_is_not``, whose premise
    (the TL tangent is not objective) is false for MITC4: both tangents are objective to
    round-off.
    """
    nodes, asm = _build_assembler()
    n_nodes = nodes.shape[0]
    t_idx = _translation_dof_indices(n_nodes)

    rows0, cols0, vals0 = asm.assemble_kt(np.zeros(asm.dofs_count, dtype=float))
    k0_tt = _coo_to_dense(rows0, cols0, vals0, asm.dofs_count)[np.ix_(t_idx, t_idx)]

    u_rigid, r = _rigid_rotation_displacement(nodes, math.radians(angle_deg))
    rows, cols, vals = getattr(asm, assemble)(u_rigid)
    k_tt = _coo_to_dense(rows, cols, vals, asm.dofs_count)[np.ix_(t_idx, t_idx)]

    k0_rot = _transform_translation_submatrix(k0_tt, r, n_nodes)
    assert _rel_err(k_tt, k0_rot, k0_tt) < TOL_OBJECTIVITY


@pytest.mark.parametrize("angle_deg", [1.0, 5.0, 30.0])
def test_rigid_rotation_produces_no_internal_force(angle_deg):
    """A finite rigid rotation about z must not strain the shell: f_int(u_rigid) ~ 0."""
    nodes, asm = _build_assembler()
    u_rigid, _ = _rigid_rotation_displacement(nodes, math.radians(angle_deg))

    f = np.asarray(asm.assemble_fint(u_rigid, True))
    # Scale: the force of a 1e-3 uniform stretch of the undeformed structure.
    rows0, cols0, vals0 = asm.assemble_kt(np.zeros(asm.dofs_count, dtype=float))
    k0 = _coo_to_dense(rows0, cols0, vals0, asm.dofs_count)
    scale = np.linalg.norm(k0 @ np.full(asm.dofs_count, 1e-3))
    assert np.linalg.norm(f) / scale < TOL_OBJECTIVITY


def test_small_displacement_kt_coro_agrees_with_tl():
    nodes, asm = _build_assembler()

    u_small = np.zeros(asm.dofs_count, dtype=float)
    for i, x in enumerate(nodes):
        # Deflexión suave: ~0.05% de L en el extremo libre
        u_small[6 * i + 2] = 5e-4 * (x[0] / L) ** 2

    rows_tl, cols_tl, vals_tl = asm.assemble_kt(u_small)
    rows_c, cols_c, vals_c = asm.assemble_kt_corotational(u_small)

    k_tl = _coo_to_dense(rows_tl, cols_tl, vals_tl, asm.dofs_count)
    k_c = _coo_to_dense(rows_c, cols_c, vals_c, asm.dofs_count)
    err = _rel_err(k_c, k_tl, k_tl)
    assert err < 1e-2
