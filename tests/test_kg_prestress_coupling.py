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
