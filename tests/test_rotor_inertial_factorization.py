"""Regression tests for PETSc factorization reuse helpers in the inertial rotor solver."""

from __future__ import annotations

import sys
from enum import Enum
from types import ModuleType
from typing import Any, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("petsc4py", reason="PETSc not available")
from petsc4py import PETSc  # noqa: E402

try:
    import precice  # noqa: F401
except (ImportError, OSError):

    class _DummyParticipant:
        def __init__(self, *args, **kwargs):
            pass

    precice_stub = cast(Any, ModuleType("precice"))
    precice_stub.Participant = _DummyParticipant
    sys.modules["precice"] = precice_stub

try:
    import _aeroelast  # noqa: F401
except (ImportError, OSError):

    class _DummyElementFamily(Enum):
        PLANE = "plane"
        SHELL = "shell"
        SOLID = "solid"

    aeroelast_stub = cast(Any, ModuleType("_aeroelast"))
    aeroelast_stub.ElementFamily = _DummyElementFamily
    sys.modules["_aeroelast"] = aeroelast_stub

try:
    from aeroelast.solvers.fsi.rotor_inertial import LinearDynamicFSIRotorInertialSolver
except (ImportError, OSError):
    pytest.skip("Rotor inertial solver not importable in this environment", allow_module_level=True)


def _create_sparse_mat(values: np.ndarray) -> PETSc.Mat:
    """Create a small sequential AIJ matrix from a dense numpy array."""
    mat = PETSc.Mat().createAIJ(values.shape, comm=PETSc.COMM_SELF)
    mat.setUp()
    for row_idx in range(values.shape[0]):
        cols = np.flatnonzero(values[row_idx]).astype(PETSc.IntType)
        if cols.size > 0:
            mat.setValues([row_idx], cols.tolist(), values[row_idx, cols])
    mat.assemble()
    return mat


def _petsc_to_dense(mat: PETSc.Mat) -> np.ndarray:
    dense = np.zeros(mat.getSize(), dtype=np.float64)
    for row_idx in range(mat.getSize()[0]):
        cols, vals = mat.getRow(row_idx)
        dense[row_idx, cols] = vals
    return dense


def test_assemble_inertial_effective_system_reuses_result_matrix() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
    a0 = 1.25
    a1 = 0.75

    k_first = np.array(
        [
            [4.0, 1.0, 0.0],
            [1.0, 5.0, 2.0],
            [0.0, 2.0, 6.0],
        ],
        dtype=np.float64,
    )
    k_second = np.array(
        [
            [10.0, 3.0, 0.0],
            [3.0, 11.0, 4.0],
            [0.0, 4.0, 12.0],
        ],
        dtype=np.float64,
    )
    mass = np.diag([1.0, 2.0, 3.0]).astype(np.float64)
    c_first = np.array(
        [
            [0.5, 0.1, 0.0],
            [0.1, 0.6, 0.2],
            [0.0, 0.2, 0.7],
        ],
        dtype=np.float64,
    )
    c_second = np.array(
        [
            [0.7, 0.2, 0.0],
            [0.2, 0.9, 0.3],
            [0.0, 0.3, 1.1],
        ],
        dtype=np.float64,
    )

    mats: list[PETSc.Mat] = []
    try:
        K_first = _create_sparse_mat(k_first)
        K_second = _create_sparse_mat(k_second)
        M = _create_sparse_mat(mass)
        C_first = _create_sparse_mat(c_first)
        C_second = _create_sparse_mat(c_second)
        mats.extend([K_first, K_second, M, C_first, C_second])

        K_eff = solver._assemble_inertial_effective_system(K_first, M, C_first, a0, a1)
        mats.append(K_eff)
        expected_first = k_first + a0 * mass + a1 * c_first
        assert_allclose(_petsc_to_dense(K_eff), expected_first)

        K_eff_reused = solver._assemble_inertial_effective_system(
            K_second,
            M,
            C_second,
            a0,
            a1,
            result=K_eff,
        )

        assert K_eff_reused is K_eff
        expected_second = k_second + a0 * mass + a1 * c_second
        assert_allclose(_petsc_to_dense(K_eff_reused), expected_second)
    finally:
        for mat in reversed(mats):
            mat.destroy()


def test_assemble_inertial_effective_system_reuses_matrix_after_zero_rows() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
    a0 = 1.25
    a1 = 0.75

    k_first = np.array(
        [
            [4.0, 1.0, 0.0],
            [1.0, 5.0, 2.0],
            [0.0, 2.0, 6.0],
        ],
        dtype=np.float64,
    )
    k_second = np.array(
        [
            [10.0, 3.0, 0.0],
            [3.0, 11.0, 4.0],
            [0.0, 4.0, 12.0],
        ],
        dtype=np.float64,
    )
    mass = np.diag([1.0, 2.0, 3.0]).astype(np.float64)
    c_first = np.array(
        [
            [0.5, 0.1, 0.0],
            [0.1, 0.6, 0.2],
            [0.0, 0.2, 0.7],
        ],
        dtype=np.float64,
    )
    c_second = np.array(
        [
            [0.7, 0.2, 0.0],
            [0.2, 0.9, 0.3],
            [0.0, 0.3, 1.1],
        ],
        dtype=np.float64,
    )

    mats: list[PETSc.Mat] = []
    vecs: list[PETSc.Vec] = []
    try:
        K_first = _create_sparse_mat(k_first)
        K_second = _create_sparse_mat(k_second)
        M = _create_sparse_mat(mass)
        C_first = _create_sparse_mat(c_first)
        C_second = _create_sparse_mat(c_second)
        mats.extend([K_first, K_second, M, C_first, C_second])

        K_eff = solver._assemble_inertial_effective_system(K_first, M, C_first, a0, a1)
        mats.append(K_eff)

        fixed_idx = np.array([1], dtype=PETSc.IntType)
        fixed_values = np.array([0.0], dtype=np.float64)
        U_fixed = K_eff.createVecRight()
        F_eff = K_eff.createVecRight()
        vecs.extend([U_fixed, F_eff])
        U_fixed.setValues(fixed_idx, fixed_values)
        U_fixed.assemble()
        F_eff.assemble()
        K_eff.zeroRows(fixed_idx, 1.0, U_fixed, F_eff)
        K_eff.assemble()

        K_eff_reused = solver._assemble_inertial_effective_system(
            K_second,
            M,
            C_second,
            a0,
            a1,
            result=K_eff,
        )

        assert K_eff_reused is K_eff
        expected_second = k_second + a0 * mass + a1 * c_second
        assert_allclose(_petsc_to_dense(K_eff_reused), expected_second)
    finally:
        for vec in reversed(vecs):
            vec.destroy()
        for mat in reversed(mats):
            mat.destroy()


def test_assemble_inertial_effective_system_handles_pattern_change() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
    a0 = 0.5
    a1 = 0.25

    k_first = np.array(
        [
            [4.0, 1.0, 0.0],
            [1.0, 5.0, 2.0],
            [0.0, 2.0, 6.0],
        ],
        dtype=np.float64,
    )
    k_second = np.array(
        [
            [4.0, 1.0, 9.0],
            [1.0, 5.0, 2.0],
            [9.0, 2.0, 6.0],
        ],
        dtype=np.float64,
    )
    mass = np.diag([1.0, 1.0, 1.0]).astype(np.float64)
    damping = np.zeros((3, 3), dtype=np.float64)

    mats: list[PETSc.Mat] = []
    try:
        K_first = _create_sparse_mat(k_first)
        K_second = _create_sparse_mat(k_second)
        M = _create_sparse_mat(mass)
        C = _create_sparse_mat(damping)
        mats.extend([K_first, K_second, M, C])

        K_eff = solver._assemble_inertial_effective_system(K_first, M, C, a0, a1)
        mats.append(K_eff)

        recreated = solver._assemble_inertial_effective_system(
            K_second,
            M,
            C,
            a0,
            a1,
            result=K_eff,
        )

        if recreated is not K_eff:
            mats.append(recreated)
        expected = k_second + a0 * mass
        assert_allclose(_petsc_to_dense(recreated), expected)
    finally:
        for mat in reversed(mats):
            mat.destroy()


def test_configure_reusable_factorization_ksp_sets_petsc_options() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
    object.__setattr__(solver, "comm", PETSc.COMM_SELF)
    solver._petsc_factorization_type = "lu"
    solver._petsc_factor_options_prefix = "rotor_inertial_test_"
    solver._petsc_factor_reuse_ordering = True
    solver._petsc_factor_reuse_fill = True
    solver._petsc_factor_mat_ordering_type = "natural"

    ksp = PETSc.KSP().create(comm=PETSc.COMM_SELF)
    opts = PETSc.Options()
    try:
        actual = solver._configure_reusable_factorization_ksp(ksp)
        configured = opts.getAll()
        assert actual == "lu"
        assert ksp.getPC().getType() == "lu"
        assert configured.get("rotor_inertial_test_pc_factor_reuse_ordering") == "true"
        assert configured.get("rotor_inertial_test_pc_factor_reuse_fill") == "true"
        assert configured.get("rotor_inertial_test_pc_factor_mat_ordering_type") == "natural"
    finally:
        ksp.destroy()
        opts.delValue("rotor_inertial_test_pc_factor_reuse_ordering")
        opts.delValue("rotor_inertial_test_pc_factor_reuse_fill")
        opts.delValue("rotor_inertial_test_pc_factor_mat_ordering_type")


def test_configure_reusable_factorization_ksp_supports_cholesky() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)
    object.__setattr__(solver, "comm", PETSc.COMM_SELF)
    solver._petsc_factorization_type = "cholesky"
    solver._petsc_factor_options_prefix = ""
    solver._petsc_factor_reuse_ordering = True
    solver._petsc_factor_reuse_fill = True
    solver._petsc_factor_mat_ordering_type = None

    ksp = PETSc.KSP().create(comm=PETSc.COMM_SELF)
    try:
        actual = solver._configure_reusable_factorization_ksp(ksp)
        assert actual == "cholesky"
        assert ksp.getPC().getType() == "cholesky"
    finally:
        ksp.destroy()


def test_fallback_factorization_type_uses_lu_after_cholesky() -> None:
    solver = object.__new__(LinearDynamicFSIRotorInertialSolver)

    assert solver._fallback_factorization_type("cholesky") == "lu"
    assert solver._fallback_factorization_type("lu") is None
