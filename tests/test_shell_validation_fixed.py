"""Comprehensive Shell Element Validation Tests - Fixed Version.

Uses exact same approach as working CCX parity tests.
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
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver
from aeroelast.solvers.modal import ModalSolver


# =============================================================================
# REFERENCE DATA
# =============================================================================

L, b, h = 1.0, 0.1, 0.001
E, nu, rho = 2.1e11, 0.3, 7800.0
STEEL = IsotropicMaterial(name="Steel", E=E, nu=nu, rho=rho)

# Expected reference values from validation tests
EXPECTED = {
    "ux": 600.0 * L / (E * b * h),
    "uy": 600.0 * L**3 / (3.0 * E * (h * b**3 / 12.0)),
    "uz": 600.0 * L**3 / (3.0 * E * (b * h**3 / 12.0)),
    "modal_1": (1.875104**2)
    * np.sqrt(E * (b * h**3 / 12.0) / (rho * b * h * L**4))
    / (2.0 * np.pi),
}

EXPECTED_RATIO = EXPECTED["uy"] / EXPECTED["ux"]


# Measured relative errors of the 8x4 mesh against the analytical references:
# FX (axial) 1.15%, FY (in-plane bending) 1.19%, FZ (out-of-plane bending)
# 1.92%.  The window is 3%: it covers those discretisation errors while still
# failing loudly for a wrong formulation, which is off by tens of percent.  Mesh
# convergence is asserted separately by
# test_axial_load_converges_to_the_analytical_solution.  The previous window was
# a flat 5% on a percentage, applied to three different physical regimes.
TOL_STATIC = 3.0  # percent

# Mesh sequence for the axial convergence study.  (2,1) is excluded because it
# is not on the asymptotic branch (measured 1.35%, below the 2.70% of (4,2)).
CONVERGENCE_MESHES = [(4, 2), (8, 4), (16, 8), (32, 16)]


def _solve_static(nx: int, ny: int, load: tuple, dof: int) -> float:
    """Solve the cantilever and return |displacement| at the free-edge centre node."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = _build_cantilever_mesh(nx=nx, ny=ny)
    prop = ShellProperty(material=STEEL, thickness=h)
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": prop},
            "span_direction": (1.0, 0.0, 0.0),
        },
    }

    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    solver.add_dirichlet_conditions([DirichletCondition(_clamped_dofs(mesh, dpn), 0.0)])
    solver.add_nodal_loads(_load_as_nodal(mesh, dpn, load))

    u = solver.solve()
    center = _center_free_edge_node(mesh)
    idx = mesh.node_id_to_index[center.id]
    return abs(u[idx * dpn + dof])


# =============================================================================
# EXACT COPY OF WORKING HELPERS FROM CCX PARITY TEST
# =============================================================================


def _build_cantilever_mesh(nx=8, ny=4):
    """Build cantilever mesh - EXACT copy."""
    Node._id_counter = 0
    MeshElement._id_counter = 0

    mesh = MeshModel()
    xs = np.linspace(0.0, L, nx + 1)
    ys = np.linspace(0.0, b, ny + 1)

    grid = {}
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

    clamped_nodes = {n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)}
    free_nodes = {n for n in mesh.nodes if np.isclose(n.x, L, atol=1e-12)}

    mesh.add_node_set(NodeSet("clamped", clamped_nodes))
    mesh.add_node_set(NodeSet("free_edge", free_nodes))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _center_free_edge_node(mesh):
    """Get center node of free edge - EXACT copy."""
    target_y = 0.5 * b
    candidates = [n for n in mesh.nodes if np.isclose(n.x, L, atol=1e-12)]
    return min(candidates, key=lambda n: abs(float(n.y) - target_y))


def _clamped_dofs(mesh, dofs_per_node):
    """Get clamped DOFs - EXACT copy."""
    dofs = []
    m = mesh.node_id_to_index
    for node in mesh.get_node_set("clamped").nodes.values():
        i0 = m[node.id] * dofs_per_node
        dofs.extend(range(i0, i0 + dofs_per_node))
    return sorted(set(dofs))


def _load_as_nodal(mesh, dofs_per_node, load6):
    """Create nodal loads - EXACT copy."""
    nodes = sorted(mesh.get_node_set("free_edge").nodes.values(), key=lambda n: n.id)
    n = max(len(nodes), 1)
    per_node = np.asarray(load6, dtype=float) / float(n)
    idx = mesh.node_id_to_index
    out = []
    for node in nodes:
        i0 = idx[node.id] * dofs_per_node
        dofs = list(range(i0, i0 + dofs_per_node))
        out.append(NodalLoad(dofs, per_node.tolist()))
    return out


# =============================================================================
# TEST CASES
# =============================================================================


class TestLinearStatic:
    """Linear static analysis tests."""

    def test_fx(self):
        """FX in-plane loading."""
        ux = _solve_static(8, 4, (600.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0)

        error = abs(ux - EXPECTED["ux"]) / EXPECTED["ux"] * 100
        print(f"\nFX: {ux * 1000:.4f} mm (ref: {EXPECTED['ux'] * 1000:.4f} mm, err {error:.2f}%)")

        assert error < TOL_STATIC, f"FX: error {error:.2f}% > {TOL_STATIC}%"

    def test_fy(self):
        """FY in-plane loading."""
        uy = _solve_static(8, 4, (0.0, 600.0, 0.0, 0.0, 0.0, 0.0), 1)

        error = abs(uy - EXPECTED["uy"]) / EXPECTED["uy"] * 100
        print(f"\nFY: {uy * 1000:.4f} mm (ref: {EXPECTED['uy'] * 1000:.4f} mm, err {error:.2f}%)")

        assert error < TOL_STATIC, f"FY: error {error:.2f}% > {TOL_STATIC}%"

    def test_fz(self):
        """FZ out-of-plane loading."""
        uz = _solve_static(8, 4, (0.0, 0.0, 600.0, 0.0, 0.0, 0.0), 2)

        error = abs(uz - EXPECTED["uz"]) / EXPECTED["uz"] * 100
        print(f"\nFZ: {uz * 1000:.4f} mm (ref: {EXPECTED['uz'] * 1000:.4f} mm, err {error:.2f}%)")

        assert error < TOL_STATIC, f"FZ: error {error:.2f}% > {TOL_STATIC}%"

    def test_ratio_physical(self):
        """UY must dominate UX: the strip is far more flexible in bending.

        Beam theory puts the ratio at 400.0 (UY/Ux = L^2 * A / I).  The previous
        window was -50%/+20%, which no factor-of-two error could fail; the
        measured ratio is 399.84, i.e. 0.04% off.
        """
        ux = _solve_static(8, 4, (600.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0)
        uy = _solve_static(8, 4, (0.0, 600.0, 0.0, 0.0, 0.0, 0.0), 1)
        ratio = uy / ux

        print(f"\nPhysical ratio: uY/uX = {ratio:.2f} (beam theory {EXPECTED_RATIO:.2f})")

        assert 0.98 * EXPECTED_RATIO <= ratio <= 1.02 * EXPECTED_RATIO, (
            f"ratio {ratio:.2f} outside 2% of the beam-theory {EXPECTED_RATIO:.2f}"
        )

    def test_axial_load_converges_to_the_analytical_solution(self):
        """The axial case must converge to P*L/(E*A), not merely land near it.

        A single-mesh comparison can only bound the discretisation error; this
        one shows that error is discretisation and not a wrong reference.  The
        measured relative errors are 2.70% (4,2), 1.15% (8,4), 0.75% (16,8) and
        0.49% (32,16): monotonically decreasing from (4,2) on, with the coarsest
        mesh excluded because it is not on the asymptotic branch.
        """
        errors = []
        for nx, ny in CONVERGENCE_MESHES:
            ux = _solve_static(nx, ny, (600.0, 0.0, 0.0, 0.0, 0.0, 0.0), 0)
            errors.append(abs(ux - EXPECTED["ux"]) / EXPECTED["ux"])

        print("\nAxial convergence (rel err vs P*L/(E*A)):")
        for (nx, ny), err in zip(CONVERGENCE_MESHES, errors, strict=True):
            print(f"  ({nx:>2},{ny:>2}): {err * 100:.3f}%")

        assert all(errors[i + 1] < errors[i] for i in range(len(errors) - 1)), (
            f"axial error must decrease under refinement, got {[f'{e * 100:.3f}%' for e in errors]}"
        )
        assert errors[-1] < 0.01, f"finest axial error {errors[-1] * 100:.3f}% must be below 1%"


class TestNonlinearStatic:
    """Nonlinear static tests."""

    def test_geometric_nonlinearity(self):
        """Large-deflection UL solve converges to a physical deflection."""
        mesh = _build_cantilever_mesh()
        prop = ShellProperty(material=STEEL, thickness=h)

        # Linear solution
        cfg_lin = {
            "solver": {},
            "elements": {
                "element_family": ElementFamily.SHELL,
                "properties": {"plate": prop},
                "span_direction": (1.0, 0.0, 0.0),
            },
        }

        solver_lin = StaticLinearSolver(mesh, cfg_lin)
        dpn = solver_lin.domain.dofs_per_node
        solver_lin.add_dirichlet_conditions([DirichletCondition(_clamped_dofs(mesh, dpn), 0.0)])
        solver_lin.add_nodal_loads(_load_as_nodal(mesh, dpn, (0.0, 0.0, 17.5, 0.0, 0.0, 0.0)))

        u_lin = solver_lin.solve()
        center = _center_free_edge_node(mesh)
        idx = mesh.node_id_to_index[center.id]
        dz_lin = abs(u_lin[idx * dpn + 2])

        # Nonlinear (UL incremental)
        cfg_nl = {
            "solver": {"continuation_steps": 32, "continuation_max_steps": 128},
            "elements": {
                "element_family": ElementFamily.SHELL,
                "properties": {"plate": prop},
                "span_direction": (1.0, 0.0, 0.0),
            },
        }

        mesh = _build_cantilever_mesh()
        solver_nl = StaticNonlinearSolver(mesh, cfg_nl)
        solver_nl.add_dirichlet_conditions([DirichletCondition(_clamped_dofs(mesh, dpn), 0.0)])
        solver_nl.add_nodal_loads(_load_as_nodal(mesh, dpn, (0.0, 0.0, 17.5, 0.0, 0.0, 0.0)))

        print(f"\nLinear estimate before nonlinear solve: {dz_lin:.3e}")
        assert dz_lin > L, "Benchmark should be strongly nonlinear"

        u_nl = np.asarray(solver_nl.solve(), dtype=np.float64)
        dz_nl = abs(u_nl[idx * dpn + 2])
        print(f"Nonlinear tip deflection: {dz_nl:.3e}")
        # Exact elastica for alpha=PL²/EI=10: w/L = 0.8034.  The first-order
        # UL gives ~0.9·L; the linear estimate is 3.3·L.
        assert dz_nl < 0.5 * dz_lin, "Nonlinear deflection must be far below the linear one"
        assert 0.4 * L < dz_nl < 1.1 * L, (
            f"Tip deflection {dz_nl:.3f} outside the physical O(L) range"
        )


class TestModal:
    """Modal tests."""

    def test_first_mode(self):
        """First modal frequency."""
        mesh = _build_cantilever_mesh()
        prop = ShellProperty(material=STEEL, thickness=h)

        cfg = {
            "solver": {"num_modes": 3},
            "elements": {
                "element_family": ElementFamily.SHELL,
                "properties": {"plate": prop},
                "span_direction": (1.0, 0.0, 0.0),
            },
        }

        solver = ModalSolver(mesh, cfg)
        dpn = solver.domain.dofs_per_node
        solver.add_dirichlet_conditions([DirichletCondition(_clamped_dofs(mesh, dpn), 0.0)])

        freqs, _ = solver.solve()
        f1 = freqs[0]

        print(f"\nModal mode 1: {f1:.3f} Hz (ref: {EXPECTED['modal_1']:.3f} Hz)")

        error = abs(f1 - EXPECTED["modal_1"]) / EXPECTED["modal_1"] * 100
        assert error < 2.0


if __name__ == "__main__":
    print("=" * 60)
    print("SHELL VALIDATION TESTS")
    print("=" * 60)

    for cls in [TestLinearStatic, TestNonlinearStatic, TestModal]:
        print(f"\n{cls.__name__}")
        print("-" * 40)

        instance = cls()
        for m in dir(instance):
            if m.startswith("test_"):
                try:
                    getattr(instance, m)()
                except AssertionError as e:
                    print(f"  FAILED: {e}")
                except Exception as e:
                    print(f"  ERROR: {e}")
