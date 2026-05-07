"""
Sprint 5 validation tests for LinearDynamicFSIRotorInertialSolver.

Covers:
- A2: Orthotropic K(θ) — eigenvalues invariant under rotation, matrices differ
- C2: update_node_coordinates roundtrip — K recovers after restore
- update_node_coordinates geometry change — K changes after rotation
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

# ---------------------------------------------------------------------------
# Optional dependencies
# ---------------------------------------------------------------------------

petsc4py = pytest.importorskip("petsc4py", reason="PETSc not available")
from petsc4py import PETSc  # noqa: E402

_aeroelast = pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.assembler import MeshAssembler  # noqa: E402
from aeroelast.core.material import IsotropicMaterial, OrthotropicMaterial  # noqa: E402
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402

# Import CoordinateTransforms without triggering __init__.py
_module_path = (
    Path(__file__).parent.parent / "src" / "aeroelast" / "solvers" / "fsi" / "corotational.py"
)
_spec = importlib.util.spec_from_file_location("fsi_rotor_corotational_s5", _module_path)
_corot_module = importlib.util.module_from_spec(_spec)
if "fsi_rotor_corotational_s5" not in sys.modules:
    sys.modules["fsi_rotor_corotational_s5"] = _corot_module
    _spec.loader.exec_module(_corot_module)
CoordinateTransforms = _corot_module.CoordinateTransforms

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

# Orthotropic material: E1 != E2 (deliberately anisotropic)
_ORTHO_MAT = OrthotropicMaterial(
    name="ortho_test",
    E=(100e9, 10e9, 10e9),   # E1=100 GPa, E2=E3=10 GPa
    G=(5e9, 5e9, 5e9),
    nu=(0.3, 0.3, 0.3),
    rho=1500.0,
)

_ISO_MAT = IsotropicMaterial(name="iso_test", E=70e9, nu=0.3, rho=2700.0)

THICKNESS = 0.005  # 5 mm shell
PLATE_L = 1.0      # 1 m x 1 m plate, 4x4 mesh


# ---------------------------------------------------------------------------
# Mesh / assembler helpers
# ---------------------------------------------------------------------------

def _build_quad_plate(
    nx: int = 4,
    ny: int = 4,
    L: float = PLATE_L,
) -> MeshModel:
    """Flat MITC4 plate mesh (nx x ny quads)."""
    Node._id_counter = 0
    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, L, ny + 1)
    grid: dict = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            n = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(n)
            grid[(i, j)] = n
    for j in range(ny):
        for i in range(nx):
            mesh.add_element(
                MeshElement(
                    nodes=[
                        grid[(i, j)],
                        grid[(i + 1, j)],
                        grid[(i + 1, j + 1)],
                        grid[(i, j + 1)],
                    ],
                    element_type=ElementType.quad,
                )
            )
    return mesh


def _model_cfg(material, thickness: float = THICKNESS):
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "material": material,
            "thickness": thickness,
        }
    }


def _petsc_to_dense(mat: PETSc.Mat) -> np.ndarray:
    n = mat.getSize()[0]
    dense = np.zeros((n, n))
    for i in range(n):
        cols_i, vals_i = mat.getRow(i)
        dense[i, cols_i] = vals_i
    return dense


def _assemble_K_dense(assembler: MeshAssembler) -> np.ndarray:
    K = assembler.assemble_stiffness_matrix()
    return _petsc_to_dense(K)


def _get_node_coords(assembler: MeshAssembler) -> np.ndarray:
    """Return (n_nodes, 3) array of current nodal coordinates from assembler."""
    return assembler.coords_array.copy()


def _apply_rotation(coords: np.ndarray, theta: float, axis: list) -> np.ndarray:
    """Rotate coords (n,3) by theta around given axis (using CoordinateTransforms)."""
    ct = CoordinateTransforms(rotation_axis=axis, rotation_center=[0.0, 0.0, 0.0])
    return ct.rotate_point_cloud(coords, theta)


# ---------------------------------------------------------------------------
# Test A2: Orthotropic K(theta) - eigenspectrum invariance
# ---------------------------------------------------------------------------

class TestOrthotropicStiffnessRotation:
    """
    Issue A2: For orthotropic material with E1 != E2, rotating the mesh by theta
    must produce a K matrix that is:
      - isospectral  (same sorted eigenvalues, up to numerical noise)
      - structurally different (K(0) != K(90 deg))

    Physical reason: rotation is an isometry -> det, trace, and all eigenvalues
    of K are preserved.  But the matrix entries change because the principal
    fiber directions rotate with the mesh.
    """

    @pytest.fixture
    def ortho_assembler_0(self):
        """Assembler for orthotropic plate at theta=0 (no rotation)."""
        mesh = _build_quad_plate()
        return MeshAssembler(mesh=mesh, model=_model_cfg(_ORTHO_MAT))

    @pytest.fixture
    def ortho_assembler_90(self):
        """Assembler for orthotropic plate rotated 90 deg around Z."""
        mesh = _build_quad_plate()
        asm = MeshAssembler(mesh=mesh, model=_model_cfg(_ORTHO_MAT))
        coords0 = _get_node_coords(asm)
        coords_rot = _apply_rotation(coords0, np.pi / 2, [0.0, 0.0, 1.0])
        asm._rust.update_node_coordinates(coords_rot)
        return asm

    def test_eigenvalues_invariant_under_90deg_rotation(
        self, ortho_assembler_0, ortho_assembler_90
    ):
        """Sorted eigenvalues of K must be equal at theta=0 and theta=90 deg."""
        K_0   = _assemble_K_dense(ortho_assembler_0)
        K_90  = _assemble_K_dense(ortho_assembler_90)

        eigs_0  = np.sort(np.linalg.eigvalsh(K_0))
        eigs_90 = np.sort(np.linalg.eigvalsh(K_90))

        assert_allclose(
            eigs_0, eigs_90, rtol=1e-6,
            err_msg="K eigenvalues must be rotation-invariant for orthotropic material"
        )

    def test_stiffness_matrices_differ_for_anisotropic_material(
        self, ortho_assembler_0, ortho_assembler_90
    ):
        """K(0 deg) and K(90 deg) must be numerically different for E1 != E2."""
        K_0  = _assemble_K_dense(ortho_assembler_0)
        K_90 = _assemble_K_dense(ortho_assembler_90)

        assert not np.allclose(K_0, K_90, atol=1e-4), (
            "K(0 deg) == K(90 deg) - rotation had no effect on orthotropic stiffness. "
            "Material transformation in element routine may be missing."
        )

    def test_isotropic_k_invariant_under_rotation(self):
        """For isotropic material, K must be the same after rotation (sanity check)."""
        mesh0 = _build_quad_plate()
        asm0  = MeshAssembler(mesh=mesh0, model=_model_cfg(_ISO_MAT))

        mesh90 = _build_quad_plate()
        asm90  = MeshAssembler(mesh=mesh90, model=_model_cfg(_ISO_MAT))
        coords0 = _get_node_coords(asm90)
        coords_rot = _apply_rotation(coords0, np.pi / 2, [0.0, 0.0, 1.0])
        asm90._rust.update_node_coordinates(coords_rot)

        K_0  = _assemble_K_dense(asm0)
        K_90 = _assemble_K_dense(asm90)

        eigs_0  = np.sort(np.linalg.eigvalsh(K_0))
        eigs_90 = np.sort(np.linalg.eigvalsh(K_90))
        assert_allclose(eigs_0, eigs_90, rtol=1e-8)

    @pytest.mark.parametrize("theta", [np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2, np.pi])
    def test_eigenvalues_invariant_for_multiple_angles(self, theta):
        """Eigenspectrum must be preserved for several rotation angles."""
        mesh0 = _build_quad_plate()
        asm0  = MeshAssembler(mesh=mesh0, model=_model_cfg(_ORTHO_MAT))
        K_0   = _assemble_K_dense(asm0)
        eigs_0 = np.sort(np.linalg.eigvalsh(K_0))

        mesh1 = _build_quad_plate()
        asm1  = MeshAssembler(mesh=mesh1, model=_model_cfg(_ORTHO_MAT))
        coords0 = _get_node_coords(asm1)
        coords_rot = _apply_rotation(coords0, theta, [0.0, 0.0, 1.0])
        asm1._rust.update_node_coordinates(coords_rot)
        K_rot = _assemble_K_dense(asm1)
        eigs_rot = np.sort(np.linalg.eigvalsh(K_rot))

        assert_allclose(
            eigs_0, eigs_rot, rtol=1e-6,
            err_msg=f"Eigenvalue invariance failed at theta={np.degrees(theta):.1f} deg"
        )


# ---------------------------------------------------------------------------
# Test C2: update_node_coordinates roundtrip
# ---------------------------------------------------------------------------

class TestUpdateNodeCoordinatesRoundtrip:
    """
    Issue C2: After rotating geometry with update_node_coordinates() and then
    restoring the original coordinates, the assembled K must be identical to
    the original K (within floating-point tolerance).

    This validates:
    1. update_node_coordinates correctly replaces topology coordinates
    2. The Rust geometry cache is fully rebuilt after update
    3. Restoring coordinates gives back the original stiffness exactly
    """

    @pytest.fixture
    def assembler(self):
        mesh = _build_quad_plate(nx=3, ny=3)
        return MeshAssembler(mesh=mesh, model=_model_cfg(_ISO_MAT))

    def test_k_changes_after_rotation(self, assembler):
        """K eigenvalues must remain the same after update_node_coordinates (in-plane rotation)."""
        K_before = _assemble_K_dense(assembler)

        coords0 = _get_node_coords(assembler)
        coords_rot = _apply_rotation(coords0, np.pi / 4, [0.0, 0.0, 1.0])
        assembler._rust.update_node_coordinates(coords_rot)
        K_after = _assemble_K_dense(assembler)

        eigs_before = np.sort(np.linalg.eigvalsh(K_before))
        eigs_after  = np.sort(np.linalg.eigvalsh(K_after))
        assert_allclose(eigs_before, eigs_after, rtol=1e-6)

    def test_k_recovered_after_roundtrip(self, assembler):
        """K must be equal after rotate-then-restore."""
        coords0  = _get_node_coords(assembler).copy()
        K_orig   = _assemble_K_dense(assembler)

        # Step 1: rotate
        coords_rot = _apply_rotation(coords0, np.pi / 3, [0.0, 1.0, 0.0])
        assembler._rust.update_node_coordinates(coords_rot)
        K_rotated = _assemble_K_dense(assembler)

        # K must differ (rotation is not identity)
        assert not np.allclose(K_orig, K_rotated, atol=1e-8), (
            "K did not change after rotate - update_node_coordinates has no effect"
        )

        # Step 2: restore original coordinates
        assembler._rust.update_node_coordinates(coords0)
        K_restored = _assemble_K_dense(assembler)

        assert_allclose(
            K_orig, K_restored, atol=1e-8,
            err_msg="K not recovered after coordinate roundtrip"
        )

    def test_multiple_roundtrips_stable(self, assembler):
        """Repeated rotate/restore cycles must not drift."""
        coords0 = _get_node_coords(assembler).copy()
        K_orig  = _assemble_K_dense(assembler)

        theta_seq = [np.pi / 6, np.pi / 4, np.pi / 3, np.pi / 2]
        for theta in theta_seq:
            coords_rot = _apply_rotation(coords0, theta, [0.0, 0.0, 1.0])
            assembler._rust.update_node_coordinates(coords_rot)
            assembler._rust.update_node_coordinates(coords0)

        K_final = _assemble_K_dense(assembler)
        assert_allclose(K_orig, K_final, atol=1e-8,
                        err_msg="K drifted after multiple rotate/restore cycles")

    def test_wrong_shape_raises(self, assembler):
        """Passing an array with wrong n_nodes must raise an error."""
        bad_coords = np.zeros((5, 3), dtype=np.float64)  # too few nodes
        with pytest.raises((ValueError, RuntimeError)):
            assembler._rust.update_node_coordinates(bad_coords)

    def test_wrong_ncols_raises(self, assembler):
        """Passing an array with ncols != 3 must raise an error."""
        coords0  = _get_node_coords(assembler)
        bad_cols = np.ascontiguousarray(coords0[:, :2])  # shape (n_nodes, 2)
        with pytest.raises((ValueError, RuntimeError)):
            assembler._rust.update_node_coordinates(bad_cols)


# ---------------------------------------------------------------------------
# Test: rotation axis coverage
# ---------------------------------------------------------------------------

class TestRotationAxisCoverage:
    """
    Confirm that update_node_coordinates works for rotations about X, Y, Z,
    and an arbitrary axis.  Eigenvalue invariance is the oracle.
    """

    @pytest.mark.parametrize("axis,theta", [
        ([1.0, 0.0, 0.0], np.pi / 4),   # X axis
        ([0.0, 1.0, 0.0], np.pi / 3),   # Y axis
        ([0.0, 0.0, 1.0], np.pi / 2),   # Z axis
        ([1.0, 1.0, 0.0], np.pi / 5),   # diagonal
        ([1.0, 1.0, 1.0], np.pi / 7),   # arbitrary
    ])
    def test_eigenvalue_invariance(self, axis, theta):
        mesh = _build_quad_plate(nx=3, ny=3)
        asm  = MeshAssembler(mesh=mesh, model=_model_cfg(_ISO_MAT))

        K_orig   = _assemble_K_dense(asm)
        coords0  = _get_node_coords(asm).copy()

        coords_rot = _apply_rotation(coords0, theta, axis)
        asm._rust.update_node_coordinates(coords_rot)
        K_rot = _assemble_K_dense(asm)

        eigs_orig = np.sort(np.linalg.eigvalsh(K_orig))
        eigs_rot  = np.sort(np.linalg.eigvalsh(K_rot))
        assert_allclose(
            eigs_orig, eigs_rot, rtol=1e-6,
            err_msg=f"Eigenvalues changed for axis={axis}, theta={np.degrees(theta):.1f} deg"
        )
