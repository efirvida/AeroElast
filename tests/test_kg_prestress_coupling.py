"""K_G uses the local CLT resultant A eps + B kappa, not output membrane stress."""

import numpy as np
import pytest

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.laminate import Laminate, Ply
from aeroelast.core.material import OrthotropicMaterial
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import CompositeShellProperty
from aeroelast.elements import ElementFamily


def _plate(symmetric):
    def material(name, e):
        g = e / 2.6
        return OrthotropicMaterial(name=name, E=(e, e, e), G=(g, g, g),
                                   nu=(0.3, 0.3, 0.3), rho=1000.0)

    stiff, soft = material("stiff", 100e9), material("soft", 10e9)
    plies = ([Ply(stiff, 0.0025, 0), Ply(soft, 0.005, 0), Ply(stiff, 0.0025, 0)]
             if symmetric else [Ply(soft, 0.005, 0), Ply(stiff, 0.005, 0)])
    lam = Laminate(plies)
    mesh = MeshModel()
    nodes = {}
    for j in range(2):
        for i in range(3):
            nodes[i, j] = Node([float(i), float(j), 0.0], geometric_node=False)
            mesh.add_node(nodes[i, j])
    for i in range(2):
        mesh.add_element(MeshElement(
            nodes=[nodes[i, 0], nodes[i + 1, 0], nodes[i + 1, 1], nodes[i, 1]],
            element_type=ElementType.quad))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    asm = MeshAssembler(mesh=mesh, model={"elements": {
        "span_direction": [1, 0, 0],
        "element_family": ElementFamily.SHELL,
        "properties": {"plate": CompositeShellProperty(laminate=lam)}}})
    return asm, lam


def _dense(mat):
    indptr, indices, data = mat.getValuesCSR()
    from scipy.sparse import csr_matrix
    return csr_matrix((data, indices, indptr), shape=mat.getSize()).toarray()


@pytest.mark.parametrize("symmetric", [False, True])
@pytest.mark.parametrize("state", ["bending", "mixed", "moment_solve"])
def test_prestress_clt_and_geometric_stiffness(symmetric, state):
    asm, lam = _plate(symmetric)
    u = np.zeros(asm.dofs_count)
    fixed = []
    f = np.zeros_like(u)
    for node in asm.nodes:
        i = asm.node_id_to_index[node.id]
        x, y, _ = node.coords
        if x == 0:
            fixed.extend(range(6 * i, 6 * i + 6))
        if x == 2:
            f[6 * i + 4] = 1.0  # pure end moment, no in-plane force
        u[6 * i + 2] = -0.001 * x**2
        u[6 * i + 4] = 0.002 * x
        if state == "mixed":
            u[6 * i] = 0.001 * x
            u[6 * i + 1] = 0.0003 * x * y
    if state == "moment_solve":
        free = np.setdiff1d(np.arange(len(u)), fixed)
        k = _dense(asm.assemble_stiffness_matrix())
        u[free] = np.linalg.solve(k[np.ix_(free, free)], f[free])
        u[fixed] = 0
        np.testing.assert_allclose((k @ u)[free], f[free], atol=1e-10)
    old, eps = asm._rust.compute_stress_field(u, 0.0, 0)
    _, top = asm._rust.compute_stress_field(u, 0.5, 1)
    _, bottom = asm._rust.compute_stress_field(u, -0.5, 1)
    slots = [0, 1, 3]
    kappa = (top[:, slots] - bottom[:, slots]) / lam.total_thickness
    expected = (eps[:, slots] @ lam.A.T + kappa @ lam.B.T) / lam.total_thickness
    actual = asm._rust.compute_prestress_field(u)
    assert actual.shape == (len(asm.elements), 3)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-8)
    # Postprocessing remains A eps/h: the new field is separate.
    np.testing.assert_allclose(old[:, slots], eps[:, slots] @ lam.A.T / lam.total_thickness,
                               rtol=1e-14, atol=1e-8)
    if state == "bending":
        np.testing.assert_allclose(old, 0, atol=1e-10)
        if not symmetric:
            assert np.linalg.norm(actual) > 1
    def geometric(stress):
        return _dense(asm.assemble_geometric_stiffness(
            stress_field={el.id: stress[i] for i, el in enumerate(asm.elements)}))
    new_kg, old_kg = geometric(actual), geometric(old[:, slots])
    rows, cols, vals = asm._rust.assemble_geometric_k_from_disp(u)
    from scipy.sparse import coo_matrix
    from_disp = coo_matrix((vals, (rows, cols)), shape=new_kg.shape).toarray()
    np.testing.assert_allclose(from_disp, new_kg, rtol=1e-13, atol=1e-10)
    if symmetric:
        np.testing.assert_allclose(actual, old[:, slots], rtol=1e-14, atol=1e-8)
        np.testing.assert_allclose(new_kg, old_kg, rtol=1e-14, atol=1e-10)
    else:
        assert np.linalg.norm(new_kg - old_kg) > 1e-3 * max(np.linalg.norm(old_kg), 1)


# ---------------------------------------------------------------------------
# T17 (L33): one implementation of the centrifugal static-solve prestress K_G.
# ``_rust.assemble_centrifugal_geometric_k`` must equal the former Python
# orchestration (load -> MUMPS solve on the free DOFs -> from_disp), which is
# written out below as the reference.
# ---------------------------------------------------------------------------

def _isotropic_hub_plate(nx=4, ny=2, hub=1.0, thickness=0.02):
    """Clamped (x == hub) isotropic plate in the x-y plane, rotating about z."""
    from aeroelast.core.material import IsotropicMaterial
    from aeroelast.core.properties import ShellProperty

    mesh = MeshModel()
    nodes = {}
    for j in range(ny + 1):
        for i in range(nx + 1):
            nodes[i, j] = Node([hub + float(i), float(j), 0.0], geometric_node=False)
            mesh.add_node(nodes[i, j])
    for j in range(ny):
        for i in range(nx):
            mesh.add_element(MeshElement(
                nodes=[nodes[i, j], nodes[i + 1, j], nodes[i + 1, j + 1], nodes[i, j + 1]],
                element_type=ElementType.quad))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    material = IsotropicMaterial(name="steel", E=2.0e11, nu=0.3, rho=7800.0)
    asm = MeshAssembler(mesh=mesh, model={"elements": {
        "span_direction": [1, 0, 0],
        "element_family": ElementFamily.SHELL,
        "properties": {"plate": ShellProperty(material, thickness)}}})
    return asm, hub


def _clamped_free_dofs(asm, x_clamp):
    fixed = set()
    for node in asm.nodes:
        if node.coords[0] == x_clamp:
            i = asm.node_id_to_index[node.id]
            fixed.update(range(6 * i, 6 * i + 6))
    return np.array(sorted(set(range(asm.dofs_count)) - fixed), dtype=np.int64)


def _reference_centrifugal_kg(asm, omega, axis, center, free_dofs, thickness_scale=1.0):
    """The pre-T17 Python orchestration, verbatim in substance.

    ``thickness_scale`` is the old isotropic rescale loop: the Rust ``rho_per_elem`` is
    volumetric for isotropic shells and is multiplied by the thickness; laminates pass 1.
    """
    from petsc4py import PETSc
    from scipy.sparse import coo_matrix

    rho = np.asarray(asm._rust.rho_per_elem(), dtype=float) * thickness_scale
    f_full = np.asarray(asm._rust.centrifugal_load(
        float(omega), list(map(float, axis)), list(map(float, center)), rho))
    rows, cols, vals = asm._rust.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    Kr = K[free_dofs][:, free_dofs].tocsr()
    A = PETSc.Mat().createAIJ(size=Kr.shape, csr=(Kr.indptr.astype(PETSc.IntType),
                                                  Kr.indices.astype(PETSc.IntType), Kr.data))
    A.assemble()
    u_red, f_red = A.createVecs()
    f_red.setArray(f_full[free_dofs])
    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    ksp.setType("preonly")
    ksp.getPC().setType("lu")
    ksp.getPC().setFactorSolverType("mumps")
    ksp.setFromOptions()
    ksp.solve(f_red, u_red)
    u_full = np.zeros(n)
    u_full[free_dofs] = u_red.getArray()
    r, c, v = asm._rust.assemble_geometric_k_from_disp(u_full)
    return coo_matrix((v, (r, c)), shape=(n, n)).toarray()


def _binding_kg(asm, omega, axis, center, free_dofs):
    from scipy.sparse import coo_matrix

    r, c, v = asm._rust.assemble_centrifugal_geometric_k(
        float(omega), list(map(float, axis)), list(map(float, center)), free_dofs)
    n = asm.dofs_count
    return coo_matrix((v, (r, c)), shape=(n, n)).toarray()


def _composite_case():
    asm, _ = _plate(symmetric=False)
    return asm, _clamped_free_dofs(asm, 0.0), 1.0


def _isotropic_case():
    asm, hub = _isotropic_hub_plate(thickness=0.02)
    return asm, _clamped_free_dofs(asm, hub), 0.02


@pytest.mark.parametrize("case", [_composite_case, _isotropic_case],
                         ids=["composite_plate", "isotropic_hub_offset"])
def test_centrifugal_geometric_k_binding_equals_python_orchestration(case):
    asm, free, t_scale = case()
    axis, center, omega = [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 12.0
    ref = _reference_centrifugal_kg(asm, omega, axis, center, free, t_scale)
    assert np.linalg.norm(ref) > 0.0  # a non-trivial membrane state
    got = _binding_kg(asm, omega, axis, center, free)
    np.testing.assert_allclose(got, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())
    # The public entry point routes through the same code (no second implementation).
    pub = _dense(asm.assemble_geometric_stiffness(
        omega=omega, rotation_axis=np.array(axis), rotation_center=np.array(center),
        free_dofs=free))
    np.testing.assert_allclose(pub, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())


def test_mass_per_area_binding_and_no_centrifugal_k_binding():
    asm, hub = _isotropic_hub_plate(thickness=0.02)
    np.testing.assert_allclose(np.asarray(asm._rust.mass_per_area()), 7800.0 * 0.02, rtol=1e-14)
    np.testing.assert_allclose(asm._rho_per_elem, 7800.0 * 0.02, rtol=1e-14)
    casm, lam = _plate(symmetric=False)
    np.testing.assert_allclose(np.asarray(casm._rust.mass_per_area()),
                               sum(p.material.rho * p.thickness for p in lam.plies), rtol=1e-12)
    # L33: the local-sigma kernel is deleted, no fallback.
    assert not hasattr(asm._rust, "assemble_centrifugal_k")


def test_centrifugal_geometric_k_requires_free_dofs():
    asm, _ = _isotropic_hub_plate()
    with pytest.raises(ValueError, match="free_dofs"):
        asm.assemble_geometric_stiffness(omega=3.0, rotation_axis=np.array([0., 0., 1.]),
                                         rotation_center=np.zeros(3))
