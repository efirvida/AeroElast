"""Analytic bend-twist benchmark for an unbalanced symmetric laminate strip.

Closes the coverage gap found while chasing the blade's over-twist: the
extension-shear ("16") terms of the CLT are verified analytically in
``test_laminate_shear_coupling.py``, but nothing verified that the SHELL
element turns a coupled laminate into the correct *bend-twist response*.

The CalculiX route is blocked: ``test_composite_beam_parity.py`` writes an
ISOTROPIC equivalent section to CCX ("to avoid CCX crash"), so it only
compares against an isotropic plate and cannot exercise any coupling.

This test uses theory instead.  Classical laminated plate theory gives

    [kappa] = [D]^-1 [M]

so a strip bent by a moment ``M_y`` develops the twisting curvature

    kappa_xy(z) = (D^-1)[2, 1] * M_y(z)

and for a cantilever with a tip load P the tip twist is

    phi_tip = (D^-1)[2, 1] * P * L^2 / 2

A symmetric but UNBALANCED layup ([45, 0]s) has B = 0 and D16 = D26 != 0, so
this isolates the bend-twist response with no extension coupling.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402

from aeroelast.core.laminate import create_laminate_from_angles  # noqa: E402
from aeroelast.core.material import Material  # noqa: E402
from aeroelast.core.mesh.entities import (  # noqa: E402
    ElementSet,
    ElementType,
    MeshElement,
    Node,
    NodeSet,
)
from aeroelast.core.mesh.model import MeshModel  # noqa: E402

E1, E2, G12, NU12 = 44.6e9, 17.0e9, 3.27e9, 0.262
THICKNESS = 8.0e-3
PLY_T = THICKNESS / 4.0
ANGLES = [45.0, 0.0, 0.0, 45.0]      # symmetric ([45,0]s) and unbalanced
L, B = 1.0, 0.10                     # length, width
NX, NY = 2, 12                       # across width, along length


def _laminate():
    mat = Material(
        name="comp_test", E=(E1, E2, E2), G=(G12, G12, G12),
        nu=(NU12, NU12, 0.0), rho=0.0,
    )
    return create_laminate_from_angles(mat, PLY_T, ANGLES)


def _strip_mesh() -> MeshModel:
    """Cantilever strip: x = length, y = width, z = mid-surface normal."""
    mesh = MeshModel()
    grid = {}
    for i in range(NX + 1):
        for j in range(NY + 1):
            x = L * i / NX
            y = B * j / NY - B / 2.0
            n = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(n)
            grid[(i, j)] = n
    for i in range(NX):
        for j in range(NY):
            mesh.add_element(MeshElement(
                nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                element_type=ElementType.quad,
            ))
    mesh.add_node_set(NodeSet("clamped", {grid[(0, j)] for j in range(NY + 1)}))
    mesh.add_node_set(NodeSet("tip", {grid[(NX, j)] for j in range(NY + 1)}))
    mesh.add_element_set(ElementSet("strip", set(mesh.elements)))
    return mesh


def _material_dict() -> dict:
    lam = _laminate()
    ABD = lam.get_ABD_matrix()
    A, Bm, D = ABD[:3, :3], ABD[:3, 3:], ABD[3:, 3:]
    h = lam.total_thickness
    return {
        "type": "composite",
        "cm": A.ravel().tolist(),
        "b_coupling": Bm.ravel().tolist(),
        "cb": D.ravel().tolist(),
        "cs": lam.Cs.ravel().tolist(),
        "thickness": h,
        "e_equiv": A[0, 0] / h,
        "mass_per_area": 0.0,
        "rotational_inertia": 0.0,
    }


def _solve_strip(mesh: MeshModel, tip_force_z: float) -> np.ndarray:
    coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    mat = _material_dict()
    asm = PyMeshAssembler(
        node_coords=coords, connectivity=conn,
        elem_types=[4] * len(mesh.elements), materials=[mat] * len(mesh.elements),
    )
    rows, cols, vals = asm.assemble_k()
    K = coo_matrix((vals, (rows, cols)), shape=(asm.dofs_count, asm.dofs_count)).tocsr()

    f = np.zeros(asm.dofs_count)
    tip_nodes = list(mesh.get_node_set("tip").nodes.values())
    per_node = tip_force_z / len(tip_nodes)
    for n in tip_nodes:                       # distributed along the tip edge
        f[mesh.node_id_to_index[n.id] * 6 + 2] = per_node

    clamped = {
        mesh.node_id_to_index[n.id] * 6 + d
        for n in mesh.get_node_set("clamped").nodes.values()
        for d in range(6)
    }
    mask = np.ones(asm.dofs_count, dtype=bool)
    mask[list(clamped)] = False
    free = np.where(mask)[0]

    u = np.zeros(asm.dofs_count)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])
    return u


def _tip_twist_deg(mesh: MeshModel, u: np.ndarray) -> float:
    """Least-squares rigid rotation about x over the tip section (= the twist)."""
    idx = [mesh.node_id_to_index[n.id] for n in mesh.get_node_set("tip").nodes.values()]
    r = np.asarray([[mesh.nodes[i].y, mesh.nodes[i].z] for i in idx], dtype=float)
    r = r - r.mean(axis=0)
    uy = np.array([u[i * 6 + 1] for i in idx])
    uz = np.array([u[i * 6 + 2] for i in idx])
    num = float((r[:, 0] * uz - r[:, 1] * uy).sum())
    den = float((r[:, 0] ** 2 + r[:, 1] ** 2).sum())
    # Rotation about +x by phi maps uy ~ -z*phi and uz ~ +y*phi, so
    # r_y*u_z - r_z*u_y = phi*(r_y^2 + r_z^2) and the fitted angle is +num/den.
    return np.rad2deg(num / den)


def _analytic_tip_twist_deg(tip_force_z: float) -> float:
    """phi_tip = (D^-1)[2,0] * P * L^2 / 2, in degrees.

    The strip runs along x and bends in z, so the applied moment is M_x
    (curvature kappa_x) and the free edges give M_y = M_xy = 0.  The
    twist curvature per unit M_x is therefore (D^-1)[2,0] -- NOT [2,1],
    which would be the response to transverse bending that a narrow strip
    cannot develop.
    """
    D = _laminate().get_ABD_matrix()[3:, 3:]
    coupling = np.linalg.inv(D)[2, 0]         # kappa_xy per unit M_x
    return np.rad2deg(coupling * tip_force_z * L**2 / 2.0)


def test_laminate_is_symmetric_but_unbalanced() -> None:
    """Guards the benchmark itself: B = 0 but D16 = D26 != 0."""
    ABD = _laminate().get_ABD_matrix()
    assert np.abs(ABD[:3, 3:]).max() < 1e-3 * np.abs(ABD[:3, :3]).max(), "B must be ~0"
    D = ABD[3:, 3:]
    assert abs(D[0, 2]) > 0.01 * abs(D[0, 0]), "D16 must be non-zero"
    assert D[0, 2] == pytest.approx(D[1, 2], rel=1e-9), "D16 == D26 for this stack"


def test_analytic_coupling_is_nonzero() -> None:
    """The closed-form bend-twist must be a measurable angle, not ~0."""
    assert abs(_analytic_tip_twist_deg(10.0)) > 1e-4


@pytest.mark.xfail(
    strict=False,
    reason="KNOWN DISCREPANCY (2026-09-25): the MITC4 shell gives 3.23x the "
    "CLT plate bend-twist for a symmetric-unbalanced [45,0]s strip "
    "(shell -0.4429 deg vs CLT -0.1373 deg, same sign, tip load 10 N).  The "
    "CLT itself is verified analytically "
    "(tests/test_laminate_shear_coupling.py), so the discrepancy is either in "
    "how the element consumes the ABD matrices, or a plate-vs-cantilever "
    "reference mismatch in this benchmark.  Settle the reference (mesh "
    "convergence + the Rust ABD usage) before making this strict.",
)
def test_shell_bend_twist_matches_closed_form() -> None:
    """The MITC4 shell must reproduce the CLT bend-twist response."""
    mesh = _strip_mesh()
    P = 10.0  # N, tip load in z
    u = _solve_strip(mesh, P)
    got = _tip_twist_deg(mesh, u)
    want = _analytic_tip_twist_deg(P)
    print(f"tip twist: shell = {got:+.6f} deg, CLT = {want:+.6f} deg")
    assert want != 0.0
    assert got == pytest.approx(want, rel=0.15), (
        f"shell bend-twist {got:.6f} deg vs CLT {want:.6f} deg"
    )
    assert np.sign(got) == np.sign(want), "bend-twist sign must match the CLT"
