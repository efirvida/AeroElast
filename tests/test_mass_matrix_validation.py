"""Mass matrix validation tests.

This module tests the mass matrix assembly against analytical solutions
and compares different mass computation methods to detect bugs.

Tests:
1. Element mass vs total mass - direct comparison
2. Lumped mass matrix trace validation
3. Consistent mass matrix validation
4. Mass matrix diagonal (lumped) vs full consistency
5. Modal frequency convergence with mesh refinement
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from aeroelast.core.bc import DirichletCondition
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.mesh.entities import ElementType, MeshElement, Node, NodeSet
from aeroelast.core.mesh.model import MeshModel
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
from aeroelast.solvers.modal import ModalSolver

logger = logging.getLogger(__name__)


# =============================================================================
# ANALYTICAL MASS FORMULAS
# =============================================================================


def quadrilateral_area(
    x1: float, y1: float, x2: float, y2: float, x3: float, y3: float, x4: float, y4: float
) -> float:
    """Area of quadrilateral via shoelace formula."""
    return 0.5 * abs(x1 * y2 + x2 * y3 + x3 * y4 + x4 * y1 - x2 * y1 - x3 * y2 - x4 * y3 - x1 * y4)


def triangle_area(x1: float, y1: float, x2: float, y2: float, x3: float, y3: float) -> float:
    """Area of triangle via determinant."""
    return 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))


# =============================================================================
# HELPERS: exact physical mass extraction from mass-matrix COO data
# =============================================================================
#
# A consistent mass matrix M_ij = rho * integral(N_i * N_j) satisfies, for every
# translational direction d,
#
#     sum_ij M[(i,d),(j,d)] = rho * integral((sum_i N_i) * (sum_j N_j))
#                           = rho * integral(1) = rho * h * A
#
# because the shape functions are a partition of unity.  That identity holds for
# every element type and every quadrature that integrates the shape functions
# exactly, so it is the correct way to read the physical mass out of M.
#
# A trace-based shortcut does NOT generalise: tr(M) = (4/3) * m holds only for a
# bilinear quad.  For a linear triangle tr(M) = (3/2) * m, which is exactly the
# 12.5% artefact this module used to assert against.

# Exact consistent-mass coefficients for a single element, per translational
# direction, as a fraction of ``rho * h * A``.  The key is the cyclic node-index
# distance ``(i - j) % n_nodes``.
#
#   tri3  (linear triangle):  M_ii = 1/6, M_ij = 1/12
#   quad4 (bilinear quad, rectangle or parallelogram): M_ii = 1/9,
#         adjacent = 1/18, opposite = 1/36
_EXACT_TRANSLATIONAL_COEFFICIENTS = {
    "tri3": {0: 2.0 / 12.0, 1: 1.0 / 12.0, 2: 1.0 / 12.0},
    "quad4": {0: 4.0 / 36.0, 1: 2.0 / 36.0, 2: 1.0 / 36.0, 3: 2.0 / 36.0},
}


def mass_matrix_dense(m_rows, m_cols, m_vals, n_dofs: int) -> np.ndarray:
    """Return the assembled mass matrix as a dense array."""
    from scipy.sparse import coo_matrix

    return coo_matrix((m_vals, (m_rows, m_cols)), shape=(n_dofs, n_dofs)).toarray()


def _direction_slice(d: int, n_dofs: int, dofs_per_node: int) -> np.ndarray:
    return np.arange(d, n_dofs, dofs_per_node)


def translational_mass_per_direction(M: np.ndarray, dofs_per_node: int) -> np.ndarray:
    """Total translational mass per direction, summed over the whole model.

    Each entry equals ``rho * h * A`` for a consistent mass matrix built from
    partition-of-unity shape functions.
    """
    n_dofs = M.shape[0]
    return np.array(
        [
            M[
                np.ix_(
                    _direction_slice(d, n_dofs, dofs_per_node),
                    _direction_slice(d, n_dofs, dofs_per_node),
                )
            ].sum()
            for d in range(3)
        ]
    )


def rotational_mass_per_direction(M: np.ndarray, dofs_per_node: int) -> np.ndarray:
    """Total rotary-inertia mass per rotation direction, summed over the model.

    Each entry equals ``rho * h**3 / 12 * A``: the shell mass matrices here
    integrate ``rho * h**3 / 12`` with the same shape-function products as the
    translational block, so the two blocks share one normalisation.
    """
    if dofs_per_node != 6:
        raise ValueError(f"rotational_mass_per_direction assumes 6 DOF/node, got {dofs_per_node}")
    n_dofs = M.shape[0]
    return np.array(
        [
            M[
                np.ix_(
                    _direction_slice(3 + d, n_dofs, dofs_per_node),
                    _direction_slice(3 + d, n_dofs, dofs_per_node),
                )
            ].sum()
            for d in range(3)
        ]
    )


def expected_translational_block(element_type: str, n_nodes: int, mass: float) -> np.ndarray:
    """Exact single-element translational mass block, per direction.

    ``mass`` is ``rho * h * A`` for that element.
    """
    coeff = _EXACT_TRANSLATIONAL_COEFFICIENTS[element_type]
    return np.array(
        [[coeff[(i - j) % n_nodes] * mass for j in range(n_nodes)] for i in range(n_nodes)]
    )


def single_element_solver(coords, element_type, material, thickness):
    """Build a one-element shell model and return ``(solver, mesh, element)``."""
    mesh = MeshModel()
    nodes = []
    for x, y in coords:
        node = Node([float(x), float(y), 0.0])
        mesh.add_node(node)
        nodes.append(node)

    element = MeshElement(nodes=nodes, element_type=element_type)
    mesh.add_element(element)
    mesh.add_node_set(NodeSet("fixed", {nodes[0]}))

    properties = {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "material": material,
            "thickness": thickness,
        },
        "solver": {"num_modes": 3},
    }
    return StaticLinearSolver(mesh, properties), mesh, element


def build_grid_mesh(element_type: str, nx: int, ny: int, length: float, width: float):
    """Build a structured shell grid, quad4 or split into tri3 pairs."""
    mesh = MeshModel()
    xs = np.linspace(0.0, length, nx + 1)
    ys = np.linspace(0.0, width, ny + 1)

    grid = {}
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            node = Node([float(x), float(y), 0.0])
            mesh.add_node(node)
            grid[(i, j)] = node

    for j in range(ny):
        for i in range(nx):
            if element_type == "tri3":
                mesh.add_element(
                    MeshElement(
                        nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)]],
                        element_type=ElementType.triangle,
                    )
                )
                mesh.add_element(
                    MeshElement(
                        nodes=[grid[(i, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                        element_type=ElementType.triangle,
                    )
                )
            else:
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


def mass_matrix_from_solver(solver, mesh) -> tuple[np.ndarray, int]:
    """Assemble the consistent mass matrix, or skip when Rust is unavailable."""
    domain = solver.domain
    if domain._rust is None:
        pytest.skip("Rust assembler not available")
    dofs_per_node = domain.dofs_per_node
    n_dofs = dofs_per_node * len(list(mesh.nodes))
    rows, cols, vals = domain._rust.assemble_m()
    return mass_matrix_dense(rows, cols, vals, n_dofs), dofs_per_node


def tributary_areas(mesh) -> dict[int, float]:
    """Integral of each node's shape function: its tributary area.

    For a partition-of-unity basis, ``integral(N_i)`` is the node's equal share
    of every element it belongs to.  The row sum of a consistent mass matrix is
    ``rho * h * integral(N_i)`` for the matching translational DOF.
    """
    areas = {node.id: 0.0 for node in mesh.nodes}
    for element in mesh.elements:
        corners = [(float(n.x), float(n.y)) for n in element.nodes]
        if len(corners) == 3:
            element_area = triangle_area(*corners[0], *corners[1], *corners[2])
        else:
            element_area = quadrilateral_area(*corners[0], *corners[1], *corners[2], *corners[3])
        share = element_area / len(element.nodes)
        for node in element.nodes:
            areas[node.id] += share
    return areas


# =============================================================================
# TEST FIXTURES
# =============================================================================


@pytest.fixture
def material_steel() -> IsotropicMaterial:
    return IsotropicMaterial(name="Steel", E=210e9, nu=0.3, rho=7850.0)


@pytest.fixture
def fem_properties() -> dict:
    return {
        "elements": {
            "element_family": ElementFamily.SHELL,
            "material": None,  # Will be set in tests
            "thickness": 0.01,
        },
        "solver": {"num_modes": 3},
    }


# =============================================================================
# TEST 1: ELEMENT MASS VS TOTAL MASS
# =============================================================================


class TestElementMassVsTotalMass:
    """Compare element-by-element mass sum vs assembled matrix trace.

    This is the primary test for mass matrix bugs.
    The total_mass should equal sum(element_masses).
    """

    @pytest.mark.parametrize("element_type", ["tri3", "quad4"])
    def test_mass_consistency(self, material_steel, fem_properties, element_type):
        """Verify assembled mass equals direct element mass sum."""
        # Parameters
        L, b, h = 1.0, 0.1, 0.01  # m, m, m
        rho = material_steel.rho

        # Build mesh based on element type
        mesh = MeshModel()

        if element_type == "tri3":
            nx, ny = 4, 2
            xs = np.linspace(0, L, nx + 1)
            ys = np.linspace(0, b, ny + 1)

            # Create triangular mesh
            nodes = {}
            for j, y in enumerate(ys):
                for i, x in enumerate(xs):
                    n = Node([float(x), float(y), 0.0])
                    mesh.add_node(n)
                    nodes[(i, j)] = n

            # Triangular elements (split quads)
            for j in range(ny):
                for i in range(nx):
                    # Lower triangle
                    mesh.add_element(
                        MeshElement(
                            nodes=[
                                nodes[(i, j)],
                                nodes[(i + 1, j)],
                                nodes[(i + 1, j + 1)],
                            ],
                            element_type=ElementType.triangle,
                        )
                    )
                    # Upper triangle
                    mesh.add_element(
                        MeshElement(
                            nodes=[
                                nodes[(i, j)],
                                nodes[(i + 1, j + 1)],
                                nodes[(i, j + 1)],
                            ],
                            element_type=ElementType.triangle,
                        )
                    )

        elif element_type == "quad4":
            nx, ny = 4, 2
            xs = np.linspace(0, L, nx + 1)
            ys = np.linspace(0, b, ny + 1)

            nodes = {}
            for j, y in enumerate(ys):
                for i, x in enumerate(xs):
                    n = Node([float(x), float(y), 0.0])
                    mesh.add_node(n)
                    nodes[(i, j)] = n

            for j in range(ny):
                for i in range(nx):
                    mesh.add_element(
                        MeshElement(
                            nodes=[
                                nodes[(i, j)],
                                nodes[(i + 1, j)],
                                nodes[(i + 1, j + 1)],
                                nodes[(i, j + 1)],
                            ],
                            element_type=ElementType.quad,
                        )
                    )

        # Fix properties for solver
        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = h

        # Solve
        solver = StaticLinearSolver(mesh, fem_properties)

        # Add boundary (fixed at one end)
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        try:
            # Get mass matrix from rust assembler
            domain = solver.domain
            if domain._rust is None:
                pytest.skip("Rust assembler not available")

            m_rows, m_cols, m_vals = domain._rust.assemble_m()
            dofs_per_node = domain.dofs_per_node
            n_dofs = dofs_per_node * len(list(mesh.nodes))
            M = mass_matrix_dense(m_rows, m_cols, m_vals, n_dofs)

            # Analytical mass per translational direction
            analytical_mass = L * b * rho * h
            per_direction = translational_mass_per_direction(M, dofs_per_node)

            logger.info(
                "Element %s: per-direction mass=%s, analytical=%.6f",
                element_type,
                per_direction,
                analytical_mass,
            )

            # Exact identity: the translational block of a consistent mass matrix
            # sums to rho*h*A in every direction, for every element type.  No
            # element-specific correction factor is involved.
            np.testing.assert_allclose(per_direction, analytical_mass, rtol=1e-12)

        except ImportError:
            pytest.skip("scipy not available")


# =============================================================================
# TEST 2: LUMPED MASS MATRIX VALIDATION
# =============================================================================


class TestLumpedMassMatrix:
    """Validate lumped mass matrix (diagonal) approach.

    Lumped mass should equal total analytical mass.
    """

    def test_rust_lumped_binding_matches_domain_matrix(self, material_steel, fem_properties):
        """High-level lumped mass assembly should use the Rust binding consistently."""
        mesh = MeshModel()
        nx, ny = 2, 1
        xs = np.linspace(0.0, 1.0, nx + 1)
        ys = np.linspace(0.0, 0.1, ny + 1)

        nodes_map = {}
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                n = Node([float(x), float(y), 0.0])
                mesh.add_node(n)
                nodes_map[(i, j)] = n

        for j in range(ny):
            for i in range(nx):
                mesh.add_element(
                    MeshElement(
                        nodes=[
                            nodes_map[(i, j)],
                            nodes_map[(i + 1, j)],
                            nodes_map[(i + 1, j + 1)],
                            nodes_map[(i, j + 1)],
                        ],
                        element_type=ElementType.quad,
                    )
                )

        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = 0.01

        solver = StaticLinearSolver(mesh, fem_properties)
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        domain = solver.domain
        if domain._rust is None:
            pytest.skip("Rust assembler not available")

        diag_rust = np.asarray(domain._rust.assemble_m_lumped(), dtype=float)

        m_lumped = domain.assemble_mass_matrix_lumped()
        diag_vec = m_lumped.getDiagonal()
        try:
            diag_domain = np.array(diag_vec.getArray(readonly=True), copy=True)
        finally:
            diag_vec.destroy()
            m_lumped.destroy()

        np.testing.assert_allclose(diag_domain, diag_rust)

    def test_lumped_vs_analytical(self, material_steel, fem_properties):
        """Lumped mass should match analytical."""
        L, b, h = 1.0, 0.1, 0.01
        rho = material_steel.rho

        # Analytical mass
        m_expected = L * b * rho * h

        # Build mesh
        mesh = MeshModel()
        nx, ny = 4, 2
        xs = np.linspace(0, L, nx + 1)
        ys = np.linspace(0, b, ny + 1)

        nodes_map = {}
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                n = Node([float(x), float(y), 0.0])
                mesh.add_node(n)
                nodes_map[(i, j)] = n

        for j in range(ny):
            for i in range(nx):
                mesh.add_element(
                    MeshElement(
                        nodes=[
                            nodes_map[(i, j)],
                            nodes_map[(i + 1, j)],
                            nodes_map[(i + 1, j + 1)],
                            nodes_map[(i, j + 1)],
                        ],
                        element_type=ElementType.quad,
                    )
                )

        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = h

        solver = StaticLinearSolver(mesh, fem_properties)

        # Boundary
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        domain = solver.domain

        if domain._rust is None:
            pytest.skip("Rust assembler not available")

        # Get the row-sum lumped diagonal from the Rust assembler
        dofs_per_node = domain.dofs_per_node
        n_dofs = dofs_per_node * len(list(mesh.nodes))
        lumped = np.asarray(domain._rust.assemble_m_lumped(), dtype=float)
        assert lumped.shape == (n_dofs,)

        # Row-sum lumping conserves the translational mass: the diagonal entries
        # of each translational direction sum to rho*h*A.  (The assembler also
        # applies a floor to massless rotational rows; that only affects the
        # rotational block, so the translational sum stays exact.)
        per_direction = np.array(
            [lumped[_direction_slice(d, n_dofs, dofs_per_node)].sum() for d in range(3)]
        )

        logger.info(
            "Lumped mass per direction: %s kg, analytical: %.6f kg",
            per_direction,
            m_expected,
        )

        np.testing.assert_allclose(per_direction, m_expected, rtol=1e-10)


# =============================================================================
# TEST 3: MODAL MASS CONVERGENCE
# =============================================================================


class TestModalMassConvergence:
    """Verify modal mass converges with mesh refinement.

    For a free-free cantilever, the first bending mode mass
    should approach the analytical Rayleigh mass.
    """

    @pytest.mark.parametrize("mesh_density", [(2, 1), (4, 2), (8, 4)])
    def test_first_mode_mass(self, material_steel, fem_properties, mesh_density):
        """First mode effective mass should converge."""
        nx, ny = mesh_density
        L, b, h = 1.0, 0.1, 0.01
        E = material_steel.E
        rho = material_steel.rho

        # Build mesh
        mesh = MeshModel()
        xs = np.linspace(0, L, nx + 1)
        ys = np.linspace(0, b, ny + 1)

        nodes_map: dict[tuple[int, int], Node] = {}
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                n = Node([float(x), float(y), 0.0])
                mesh.add_node(n)
                nodes_map[(i, j)] = n

        for j in range(ny):
            for i in range(nx):
                mesh.add_element(
                    MeshElement(
                        nodes=[
                            nodes_map[(i, j)],
                            nodes_map[(i + 1, j)],
                            nodes_map[(i + 1, j + 1)],
                            nodes_map[(i, j + 1)],
                        ],
                        element_type=ElementType.quad,
                    )
                )

        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = h

        # Modal solve
        solver = ModalSolver(mesh, fem_properties)

        # Fixed at root
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        try:
            freqs, modes = solver.solve()

            # First mode for cantilever beam:
            # ω₁ = 1.875² √(EI/ρAL⁴) = 3.516² √(EI/ρAL⁴)
            # f₁ = ω₁/(2π) = 3.516/(2π) √(EI/ρAL⁴)
            I = b * h**3 / 12
            A = b * h
            f_analytical = (3.516 / (2 * np.pi)) * np.sqrt(E * I / (rho * A * L**4))

            freq = freqs[0]
            error = abs(freq - f_analytical) / f_analytical

            logger.info(
                "Mesh %sx%s: f1=%.1f Hz, analytical=%.1f Hz, error=%.1f%%",
                nx,
                ny,
                freq,
                f_analytical,
                error * 100,
            )

            # Both branches of the original `0.05 if nx <= 2 else 0.05` were the
            # same value, so it never branched; this is the 5% it always applied.
            tol = 0.05
            assert error < tol, f"Frequency error: {error * 100:.1f}% (tol {tol * 100:.0f}%)"

        except AssertionError:
            raise  # Re-raise assert errors, don't skip
        except Exception as e:
            pytest.skip(f"Modal solve failed: {e}")


# =============================================================================
# TEST 4: CONSISTENT MASS TOTAL
# =============================================================================


class TestConsistentMassTotal:
    """Validate the assembled consistent mass against the total model mass.

    The translational block of a consistent mass matrix sums to ``rho*h*A`` in
    every direction.  The trace does NOT equal the mass: it equals ``4/3 * m``
    for a bilinear quad and ``3/2 * m`` for a linear triangle, so a single
    trace-based factor is wrong for at least one of them.
    """

    def test_consistent_mass_total_equals_analytical(self, material_steel, fem_properties):
        """The translational block must sum to the analytical mass."""
        L, b, h = 1.0, 0.1, 0.01
        rho = material_steel.rho

        # Analytical mass
        m_total = L * b * rho * h

        # Build mesh
        mesh = MeshModel()
        nx, ny = 4, 2
        xs = np.linspace(0, L, nx + 1)
        ys = np.linspace(0, b, ny + 1)

        nodes_map = {}
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                n = Node([float(x), float(y), 0.0])
                mesh.add_node(n)
                nodes_map[(i, j)] = n

        for j in range(ny):
            for i in range(nx):
                mesh.add_element(
                    MeshElement(
                        nodes=[
                            nodes_map[(i, j)],
                            nodes_map[(i + 1, j)],
                            nodes_map[(i + 1, j + 1)],
                            nodes_map[(i, j + 1)],
                        ],
                        element_type=ElementType.quad,
                    )
                )

        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = h

        solver = StaticLinearSolver(mesh, fem_properties)

        # Boundary
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        domain = solver.domain

        if domain._rust is None:
            pytest.skip("Rust assembler not available")

        m_rows, m_cols, m_vals = domain._rust.assemble_m()
        dofs_per_node = domain.dofs_per_node
        n_dofs = dofs_per_node * len(list(mesh.nodes))
        M = mass_matrix_dense(m_rows, m_cols, m_vals, n_dofs)

        # Exact total mass per direction, read from the translational block
        per_direction = translational_mass_per_direction(M, dofs_per_node)

        logger.info("Consistent mass per direction: %s, analytical: %.6f", per_direction, m_total)

        np.testing.assert_allclose(per_direction, m_total, rtol=1e-12)


# =============================================================================
# TEST 5: MASS MATRIX SYMMETRY
# =============================================================================


class TestMassMatrixSymmetry:
    """Verify mass matrix is symmetric (M = M^T)."""

    def test_symmetry(self, material_steel, fem_properties):
        """Mass matrix should be symmetric."""
        L, b, h = 1.0, 0.1, 0.01

        mesh = MeshModel()
        nx, ny = 4, 2
        xs = np.linspace(0, L, nx + 1)
        ys = np.linspace(0, b, ny + 1)

        nodes = {}
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                n = Node([float(x), float(y), 0.0])
                mesh.add_node(n)
                nodes[(i, j)] = n

        for j in range(ny):
            for i in range(nx):
                mesh.add_element(
                    MeshElement(
                        nodes=[
                            nodes[(i, j)],
                            nodes[(i + 1, j)],
                            nodes[(i + 1, j + 1)],
                            nodes[(i, j + 1)],
                        ],
                        element_type=ElementType.quad,
                    )
                )

        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = h

        solver = StaticLinearSolver(mesh, fem_properties)

        # Boundary
        fixed = [n for n in mesh.nodes if np.isclose(n.x, 0.0, atol=1e-12)]
        mesh.add_node_set(NodeSet("fixed", set(fixed)))
        fixed_dofs = solver.get_dofs_by_nodeset_name("fixed")
        solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

        domain = solver.domain

        try:
            from scipy.sparse import coo_matrix

            domain.assemble_mass_matrix()
            m_rows, m_cols, m_vals = domain._rust.assemble_m()

            dofs = domain.dofs_per_node * len(mesh.nodes)
            M_sparse = coo_matrix((m_vals, (m_rows, m_cols)), shape=(dofs, dofs))

            # Check M - M^T should be ~zero
            M_diff = M_sparse - M_sparse.T
            max_diff = np.abs(M_diff.data).max() if M_diff.nnz > 0 else 0.0

            logger.info("Symmetry check: max_diff=%.2e", max_diff)

            assert max_diff < 1e-10, f"Mass matrix not symmetric: {max_diff}"

        except ImportError:
            pytest.skip("scipy not available")


# =============================================================================
# TEST 6: EXACT CONSISTENT-MASS COEFFICIENTS (single element)
# =============================================================================

_REFERENCE_ELEMENTS = [
    ("tri3", ElementType.triangle, [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], 0.5),
    ("quad4", ElementType.quad, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], 1.0),
]


class TestExactConsistentMassCoefficients:
    """The single-element mass matrix must match the closed-form integral.

    For one element the translational block carries no discretisation error, so
    the tolerance is machine precision.  A global mass sum cannot catch a wrong
    shape-function product or a wrong quadrature weight: a wrong distribution
    still conserves the total.
    """

    @pytest.mark.parametrize(("name", "element_type", "coords", "area"), _REFERENCE_ELEMENTS)
    def test_translational_block(self, material_steel, name, element_type, coords, area):
        """M_ij must equal the exact integral of rho*h*N_i*N_j."""
        thickness = 0.01
        solver, mesh, element = single_element_solver(
            coords, element_type, material_steel, thickness
        )
        M, dofs_per_node = mass_matrix_from_solver(solver, mesh)

        expected = expected_translational_block(
            name, len(element.nodes), material_steel.rho * thickness * area
        )
        indices = [mesh.node_id_to_index[node.id] for node in element.nodes]

        for direction in range(3):
            dofs = [index * dofs_per_node + direction for index in indices]
            np.testing.assert_allclose(M[np.ix_(dofs, dofs)], expected, rtol=1e-12, atol=1e-18)

    @pytest.mark.parametrize(("name", "element_type", "coords", "area"), _REFERENCE_ELEMENTS)
    def test_rotary_inertia_block(self, material_steel, name, element_type, coords, area):
        """Each rotational direction must carry rho*h**3/12*A of rotary inertia."""
        thickness = 0.01
        solver, mesh, element = single_element_solver(
            coords, element_type, material_steel, thickness
        )
        M, dofs_per_node = mass_matrix_from_solver(solver, mesh)

        expected = material_steel.rho * thickness**3 / 12.0 * area
        np.testing.assert_allclose(
            rotational_mass_per_direction(M, dofs_per_node), expected, rtol=1e-12
        )


# =============================================================================
# TEST 7: CONSISTENT-MASS DISTRIBUTION
# =============================================================================


class TestConsistentMassDistribution:
    """The mass distribution, not only its total, must be physical.

    The row sum of a consistent mass matrix is the lumped mass distribution:
    for a translational DOF it equals ``rho * h * integral(N_i)``, the node's
    tributary area.  Putting all mass on one node would conserve the total and
    still pass a sum-only check, so this is the test that pins the distribution.
    """

    @pytest.mark.parametrize("element_type", ["tri3", "quad4"])
    def test_row_sums_match_tributary_areas(self, material_steel, fem_properties, element_type):
        length, width, thickness = 1.0, 0.1, 0.01
        mesh = build_grid_mesh(element_type, 4, 2, length, width)
        fem_properties["elements"]["material"] = material_steel
        fem_properties["elements"]["thickness"] = thickness
        solver = StaticLinearSolver(mesh, fem_properties)
        M, dofs_per_node = mass_matrix_from_solver(solver, mesh)

        areas = tributary_areas(mesh)
        for node in mesh.nodes:
            expected = material_steel.rho * thickness * areas[node.id]
            base = mesh.node_id_to_index[node.id] * dofs_per_node
            for direction in range(3):
                row_sum = M[base + direction, :].sum()
                assert row_sum == pytest.approx(expected, rel=1e-12), (
                    f"node {node.id} direction {direction}: row sum {row_sum} "
                    f"!= rho*h*tributary_area {expected}"
                )
