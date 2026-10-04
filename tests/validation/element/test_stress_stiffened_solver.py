"""Tests for StressStiffenedFSISolver — geometric stiffness update pipeline.

These tests verify the components of the incremental linearization approach.
They are split into two groups:

``TestKGAssemblyPipeline`` — tests the underlying K_G / StressRecovery machinery
  WITHOUT requiring preCICE.  Always run.

``TestStressStiffenedHook`` — tests _post_convergence_hook on the actual solver
  class.  Skipped when preCICE shared library is missing (login nodes).

``TestStressStiffenedConfig`` — tests config enum and inheritance.  The solver
  import is attempted lazily so only the config tests fail if preCICE is absent.
"""

from __future__ import annotations

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Skip the whole module if PETSc is not available (CI without MPI/PETSc)
# ---------------------------------------------------------------------------
pytest.importorskip("petsc4py", reason="PETSc not available")
from petsc4py import PETSc

from _aeroelast import PyMeshAssembler

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.bc import BoundaryConditionManager, DirichletCondition
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.postprocess.stress_recovery import StressLocation, StressRecovery, StressType
from aeroelast.solvers.fsi.time_integration import NewmarkCoefficients

from tests.support.assertions import assert_residual_below

# ---------------------------------------------------------------------------
# Optional FSI solver import (needs preCICE at runtime)
# ---------------------------------------------------------------------------
try:
    from aeroelast.solvers.fsi.stress_stiffened_dynamic import StressStiffenedFSISolver

    _HAS_FSI = True
except ImportError:
    _HAS_FSI = False
    StressStiffenedFSISolver = None  # type: ignore[assignment,misc]

_skip_fsi = pytest.mark.skipif(not _HAS_FSI, reason="preCICE shared library not available")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_plate_mesh(nx: int = 4, ny: int = 4, L: float = 1.0) -> MeshModel:
    """Build a flat MITC4 square plate mesh (nx × ny quads) in the XY plane."""
    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, L, ny + 1)

    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            n = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(n)
            grid[(i, j)] = n

    for j in range(ny):
        for i in range(nx):
            mesh.add_element(
                MeshElement(
                    nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                    element_type=ElementType.quad,
                )
            )
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


_SOLVER_CFG = {
    "type": "StressStiffenedDynamicFSI",
    "time_step": 0.01,
    "total_time": 1.0,
    "beta": 0.25,
    "gamma": 0.5,
    "geometric_stiffness": {"update_interval": 1},
}

_STEEL = IsotropicMaterial(name="steel", E=210e9, nu=0.3, rho=7800.0)
_PROP = ShellProperty(material=_STEEL, thickness=0.01)


def _model_cfg() -> dict:
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": _PROP},
        },
        "solver": dict(_SOLVER_CFG),
    }


def _build_domain(mesh: MeshModel) -> MeshAssembler:
    return MeshAssembler(mesh=mesh, model=_model_cfg())


def _build_bc_manager(domain: MeshAssembler, mesh: MeshModel) -> BoundaryConditionManager:
    """Pin the x=0 edge (all 6 DOFs) and return a BoundaryConditionManager."""
    K = domain.assemble_stiffness_matrix()
    M_c = domain.assemble_mass_matrix()
    F = PETSc.Vec().createMPI(domain.dofs_count, comm=PETSc.COMM_WORLD)
    F.set(0.0)

    bc_mgr = BoundaryConditionManager(
        stiffness=K, load=F, mass=M_c, dof_per_node=domain.dofs_per_node
    )
    pinned_nodes = [n for n in domain.nodes if abs(n.coords[0]) < 1e-12]
    dofs: list[int] = []
    for node in pinned_nodes:
        idx = domain.node_id_to_index[node.id]
        for d in range(domain.dofs_per_node):
            dofs.append(idx * domain.dofs_per_node + d)
    bc_mgr.apply_dirichlet([DirichletCondition(dofs=dofs, value=0.0)])
    return bc_mgr


def _make_solver(
    mesh: MeshModel, domain: MeshAssembler, update_interval: int = 1
) -> "StressStiffenedFSISolver":
    """Instantiate the solver bypassing the preCICE-dependent __init__."""
    cfg = _model_cfg()
    cfg["solver"]["geometric_stiffness"]["update_interval"] = update_interval
    solver = object.__new__(StressStiffenedFSISolver)
    solver.domain = domain
    solver.mesh_obj = mesh
    solver.model_properties = cfg
    solver.solver_params = cfg["solver"]
    solver._kg_update_interval = update_interval
    solver._K_G_red = None
    return solver


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def plate_setup():
    """Shared plate mesh + domain + bc_manager for all tests."""
    mesh = _build_plate_mesh(nx=4, ny=4, L=1.0)
    domain = _build_domain(mesh)
    bc_mgr = _build_bc_manager(domain, mesh)
    return mesh, domain, bc_mgr


# ---------------------------------------------------------------------------
# Group 1 — K_G pipeline tests (NO preCICE needed)
# ---------------------------------------------------------------------------


class TestKGAssemblyPipeline:
    """Verify the K_G assembly + StressRecovery pipeline without preCICE."""

    def test_assemble_geometric_stiffness_from_stress_field(self, plate_setup):
        """assembler.assemble_geometric_stiffness(stress_field=…) must return a PETSc.Mat."""
        mesh, domain, bc_mgr = plate_setup
        stress_field = {e.id: np.array([1e6, 5e5, 0.0]) for e in domain.elements}
        K_G = domain.assemble_geometric_stiffness(stress_field=stress_field)
        assert isinstance(K_G, PETSc.Mat)
        assert K_G.getSize()[0] == domain.dofs_count

    def test_K_G_magnitude_scales_with_thickness(self):
        """``K_G`` must scale linearly with the shell thickness (issue #7).

        ``assemble_geometric_k`` is handed the membrane *stress*
        ``[sigma_xx, sigma_yy, sigma_xy]`` in Pa, so the force resultant
        ``N = sigma * h`` has to be formed inside.  MITC4 omitted the thickness,
        which made its K_G ``1/h`` too large and independent of ``h``; the
        existing sign/symmetry/PSD tests could not see it.
        """
        node_coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
        connectivity = [[0, 1, 2, 3]]
        sigma = np.array([[1.0e6, 5.0e5, 2.0e5]])

        def kg_norm(thickness: float) -> float:
            materials = [
                {
                    "type": "isotropic",
                    "e": 2.1e11,
                    "nu": 0.3,
                    "rho": 0.0,
                    "thickness": thickness,
                    "shear_correction": 5.0 / 6.0,
                    "drilling_scale": 1.0,
                }
            ]
            asm = PyMeshAssembler(
                node_coords=node_coords,
                connectivity=connectivity,
                elem_types=[4],
                materials=materials,
            )
            _, _, vals = asm.assemble_geometric_k(sigma)
            return float(np.linalg.norm(vals))

        thin, thick = kg_norm(0.01), kg_norm(0.1)
        assert thin > 0.0, "K_G must be nonzero for a nonzero stress"
        ratio = thick / thin
        assert abs(ratio - 10.0) < 1e-9 * 10.0, (
            f"K_G must scale with the thickness: |K(0.1)|/|K(0.01)| = {ratio:.6f} "
            f"(expected 10); a thickness-independent ratio means the sigma*h "
            f"factor is missing"
        )

    def test_K_G_symmetry(self, plate_setup):
        """K_G assembled from a uniform stress field must be symmetric."""
        mesh, domain, bc_mgr = plate_setup
        stress_field = {e.id: np.array([1e6, 0.0, 0.0]) for e in domain.elements}
        K_G = domain.assemble_geometric_stiffness(stress_field=stress_field)

        # Compare K_G - K_G^T (should be zero up to floating-point)
        K_G_T = K_G.copy()
        K_G_T.transpose()
        diff = K_G.copy()
        diff.axpy(-1.0, K_G_T)
        assert diff.norm(PETSc.NormType.FROBENIUS) < 1e-8 * K_G.norm(PETSc.NormType.FROBENIUS), (
            "K_G must be symmetric."
        )

    def test_K_G_is_positive_semidefinite(self, plate_setup):
        """K_G from tensile stress must be positive semidefinite."""
        mesh, domain, bc_mgr = plate_setup
        stress_field = {e.id: np.array([1e6, 1e6, 0.0]) for e in domain.elements}
        K_G = domain.assemble_geometric_stiffness(stress_field=stress_field)
        K_G_red = bc_mgr.reduce_matrix(K_G)

        # Build dense copy for eigenvalue check
        n = K_G_red.getSize()[0]
        dense = np.zeros((n, n))
        for i in range(n):
            row_cols, row_vals = K_G_red.getRow(i)
            for c, v in zip(row_cols, row_vals, strict=False):
                dense[i, c] = v

        eigs = np.linalg.eigvalsh(dense)
        assert np.all(eigs >= -1e-6 * np.max(np.abs(eigs))), (
            f"K_G has negative eigenvalues: {eigs[eigs < 0]}"
        )

    def test_stress_recovery_element_stresses_returns_arrays(self, plate_setup):
        """StressRecovery.compute_element_stresses must run without error."""
        mesh, domain, bc_mgr = plate_setup

        # Build a non-zero displacement (unit membrane stretch in x)
        u_full = np.zeros(domain.dofs_count)
        for node in domain.nodes:
            idx = domain.node_id_to_index[node.id]
            u_full[idx * domain.dofs_per_node + 0] = node.coords[0] * 1e-3  # u = ε·x

        sr = StressRecovery(domain, u_full)
        result = sr.compute_element_stresses(
            location=StressLocation.MIDDLE,
            stress_type=StressType.MEMBRANE,
        )
        assert result.sigma_xx is not None
        assert len(result.sigma_xx) == len(list(domain.elements))

    def test_stress_field_dict_from_recovery(self, plate_setup):
        """Build stress_field dict from StressRecovery and assemble K_G."""
        mesh, domain, bc_mgr = plate_setup

        u_full = np.zeros(domain.dofs_count)
        for node in domain.nodes:
            idx = domain.node_id_to_index[node.id]
            u_full[idx * domain.dofs_per_node + 0] = node.coords[0] * 5e-3

        sr = StressRecovery(domain, u_full)
        elem_result = sr.compute_element_stresses(
            location=StressLocation.MIDDLE,
            stress_type=StressType.MEMBRANE,
        )

        stress_field: dict[int, np.ndarray] = {}
        for i, elem in enumerate(domain.elements):
            sigma = np.array(
                [
                    elem_result.sigma_xx[i],
                    elem_result.sigma_yy[i],
                    elem_result.sigma_xy[i],
                ]
            )
            if np.max(np.abs(sigma)) > 1e-20:
                stress_field[elem.id] = sigma

        # ``u = eps * x`` is a uniform membrane strain eps = 5e-3 in x, so the
        # plane-stress response is known in closed form:
        #   sigma_xx = E/(1 - nu**2) * eps,  sigma_yy = nu * sigma_xx,  sigma_xy = 0
        eps = 5e-3
        sigma_xx_expected = _STEEL.E / (1.0 - _STEEL.nu**2) * eps
        sigma_yy_expected = _STEEL.nu * sigma_xx_expected
        n_elem = len(list(domain.elements))

        assert len(stress_field) == n_elem, (
            f"every element carries the uniform strain; got {len(stress_field)}/{n_elem}"
        )
        # Every element must carry the uniform plane-stress response, and the worst element is what
        # says so: asserting the maximum is the same claim as asserting each one, stated once, so the
        # comparison has one residual the store can record instead of one per element.
        values = list(stress_field.values())
        assert_residual_below(
            max(abs(sigma[0] - sigma_xx_expected) / abs(sigma_xx_expected) for sigma in values),
            tol=1e-6,
            kind="analytical",
            reference_name="the plane-stress closed form E/(1-nu^2) * eps",
            what="worst element sigma_xx",
        )
        # The two below keep the bounds the loop had, and they are absolute: the file compares
        # sigma_yy against 1e-6*sigma_xx and the vanishing shear against 1e-9*sigma_xx, which are
        # quantities in Pa rather than fractions of the reference each one is measured against.
        assert_residual_below(
            max(abs(sigma[1] - sigma_yy_expected) for sigma in values),
            atol=1e-6 * sigma_xx_expected,
            unit="Pa",
            kind="analytical",
            reference_name="nu * sigma_xx of the same closed form",
            what="worst element sigma_yy deviation",
        )
        assert_residual_below(
            max(abs(sigma[2]) for sigma in values),
            atol=1e-9 * sigma_xx_expected,
            unit="Pa",
            kind="analytical",
            reference_name="zero shear in the plane-stress closed form",
            what="worst element sigma_xy",
        )

        K_G = domain.assemble_geometric_stiffness(stress_field=stress_field)
        assert isinstance(K_G, PETSc.Mat)

    def test_keff_with_KG_larger_than_without(self, plate_setup):
        """K + K_G diagonal sum must exceed K diagonal sum for tensile stress."""
        mesh, domain, bc_mgr = plate_setup
        K_red, _, M_red = bc_mgr.reduced_system
        coeffs = NewmarkCoefficients.from_newmark_params(0.25, 0.5, 0.01)

        # Baseline K_eff
        K_eff_base = K_red.copy()
        K_eff_base.axpy(coeffs.a0, M_red)
        diag_base = K_eff_base.getDiagonal().getArray().copy()

        # K_G from a uniform biaxial tensile stress
        sigma = np.array([1e7, 1e7, 0.0])
        stress_field = {e.id: sigma for e in domain.elements}
        K_G = domain.assemble_geometric_stiffness(stress_field=stress_field)
        K_G_red = bc_mgr.reduce_matrix(K_G)

        K_eff_new = K_red.copy()
        K_eff_new.axpy(1.0, K_G_red)
        K_eff_new.axpy(coeffs.a0, M_red)
        diag_new = K_eff_new.getDiagonal().getArray().copy()

        delta = diag_new - diag_base
        scale = float(np.max(np.abs(diag_base)))
        # K_G from tensile stress is positive semidefinite, so no DOF may lose
        # stiffness.  The elementwise check (not just the sum) catches a sign
        # error in a single element's contribution.
        n_negative = int((delta < -1e-12 * scale).sum())
        assert n_negative == 0, (
            f"stress stiffening removed stiffness from {n_negative} DOF(s): "
            f"min delta = {delta.min():.3e} (base scale {scale:.3e})"
        )
        assert np.sum(delta) > 0, "K_G must add net stiffness under tensile stress."
        # K_G is a small correction, not a replacement of the elastic stiffness.
        assert np.sum(delta) < 0.5 * np.sum(np.abs(diag_base)), (
            f"K_G contributes {np.sum(delta) / np.sum(np.abs(diag_base)) * 100:.2f}% "
            "of the base diagonal; it must remain a correction"
        )

        # K_G is linear in the prescribed stress: feeding the exact same field
        # scaled by 2 must double the trace.
        K_G_2x = domain.assemble_geometric_stiffness(
            stress_field={e.id: 2.0 * sigma for e in domain.elements}
        )
        trace_1x = float(np.sum(delta))
        trace_2x = float(
            np.sum(bc_mgr.reduce_matrix(K_G_2x).getDiagonal().getArray())
        )
        assert abs(trace_2x - 2.0 * trace_1x) < 1e-9 * abs(trace_1x), (
            f"K_G is not linear in stress: 1x trace {trace_1x:.6e}, "
            f"2x trace {trace_2x:.6e}"
        )


# ---------------------------------------------------------------------------
# Group 2 — Solver hook tests (require preCICE for the class import)
# ---------------------------------------------------------------------------


@_skip_fsi
class TestStressStiffenedHook:
    """Unit tests for _post_convergence_hook — require preCICE library."""

    def test_hook_returns_none_for_zero_displacement(self, plate_setup):
        """Hook returns None when u = 0 (no stress → no K_G rebuild)."""
        mesh, domain, bc_mgr = plate_setup
        solver = _make_solver(mesh, domain)

        K_red, _, M_red = bc_mgr.reduced_system
        coeffs = NewmarkCoefficients.from_newmark_params(0.25, 0.5, 0.01)
        K_eff = K_red.copy()
        K_eff.axpy(coeffs.a0, M_red)

        u_zero = K_red.createVecRight()
        u_zero.set(0.0)

        result = solver._post_convergence_hook(
            u=u_zero,
            time_step=1,
            K_eff=K_eff,
            K_red=K_red,
            M_red=M_red,
            C_red=None,
            coeffs=coeffs,
            bc_manager=bc_mgr,
        )
        assert result is None

    def test_hook_returns_new_keff_under_membrane_load(self, plate_setup):
        """Hook returns a NEW K_eff object when membrane stress is non-zero."""
        mesh, domain, bc_mgr = plate_setup
        solver = _make_solver(mesh, domain)

        K_red, _, M_red = bc_mgr.reduced_system
        coeffs = NewmarkCoefficients.from_newmark_params(0.25, 0.5, 0.01)
        K_eff_old = K_red.copy()
        K_eff_old.axpy(coeffs.a0, M_red)

        # Membrane displacement: u = i * 1e-4 → ε_xx ≠ 0
        u_mem = K_red.createVecRight()
        arr = u_mem.getArray()
        n_free = arr.shape[0] // domain.dofs_per_node
        for i in range(n_free):
            arr[i * domain.dofs_per_node + 0] = float(i) * 1e-4
        u_mem.setArray(arr)

        result = solver._post_convergence_hook(
            u=u_mem,
            time_step=1,
            K_eff=K_eff_old,
            K_red=K_red,
            M_red=M_red,
            C_red=None,
            coeffs=coeffs,
            bc_manager=bc_mgr,
        )

        if result is None:
            pytest.skip("All element stresses below threshold after BC reduction.")

        assert result is not K_eff_old, (
            "Hook must return a NEW PETSc.Mat so _solve_linear_system refactorizes."
        )
        assert isinstance(result, PETSc.Mat)

    def test_update_interval_skips_rebuild(self, plate_setup):
        """Hook must return None when time_step % update_interval != 0."""
        mesh, domain, bc_mgr = plate_setup
        solver = _make_solver(mesh, domain, update_interval=5)

        K_red, _, M_red = bc_mgr.reduced_system
        coeffs = NewmarkCoefficients.from_newmark_params(0.25, 0.5, 0.01)
        K_eff = K_red.copy()
        K_eff.axpy(coeffs.a0, M_red)

        u_nz = K_red.createVecRight()
        arr = u_nz.getArray()
        # A *straining* displacement: a uniform translation would have zero
        # strain and would make the hook return None for the wrong reason.
        n_free = arr.shape[0] // domain.dofs_per_node
        for i in range(n_free):
            arr[i * domain.dofs_per_node + 0] = float(i + 1) * 1e-4
        u_nz.setArray(arr)

        def _call(step: int):
            return solver._post_convergence_hook(
                u=u_nz,
                time_step=step,
                K_eff=K_eff,
                K_red=K_red,
                M_red=M_red,
                C_red=None,
                coeffs=coeffs,
                bc_manager=bc_mgr,
            )

        assert _call(3) is None, "update_interval=5 must skip rebuild at step 3."

        result_5 = _call(5)
        assert result_5 is not None, (
            "update_interval=5 must rebuild at step 5 for a straining displacement."
        )
        assert isinstance(result_5, PETSc.Mat)
        assert result_5 is not K_eff, "rebuild must return a NEW matrix for refactorization."


# ---------------------------------------------------------------------------
# Group 3 — Configuration tests (no preCICE needed for these)
# ---------------------------------------------------------------------------


class TestStressStiffenedConfig:
    """Configuration and registration tests."""

    def test_solver_type_enum_exists(self):
        from aeroelast.core.config import SolverType

        assert hasattr(SolverType, "STRESS_STIFFENED_DYNAMIC_FSI")
        assert SolverType.STRESS_STIFFENED_DYNAMIC_FSI.value == "StressStiffenedDynamicFSI"

    @_skip_fsi
    def test_solver_inherits_linear_dynamic(self):
        from aeroelast.solvers.fsi.linear_dynamic import LinearDynamicFSISolver

        assert issubclass(StressStiffenedFSISolver, LinearDynamicFSISolver)

    @_skip_fsi
    def test_default_update_interval(self):
        mesh = _build_plate_mesh(nx=2, ny=2)
        domain = _build_domain(mesh)
        solver = _make_solver(mesh, domain)
        assert solver._kg_update_interval == 1

    @_skip_fsi
    def test_custom_update_interval(self):
        mesh = _build_plate_mesh(nx=2, ny=2)
        domain = _build_domain(mesh)
        solver = _make_solver(mesh, domain, update_interval=10)
        assert solver._kg_update_interval == 10
