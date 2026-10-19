"""Analytic plate-buckling anchor for the MITC4 geometric stiffness (T7 (i), rotating-local-tip-modes).

Why: the rotating IEA-15 blade buckles its tip TE panels at lam ~0.5 with MITC4, where CalculiX S8R
on the same mesh and prestress gives ~2.3 (``odd/tasks/rotating-local-tip-modes.md`` L13-L15). No test
validated out-of-plane plate buckling of MITC4 against a closed form; this one does, for the simplest
case: an isotropic, simply supported, thin rectangular plate under uniform uniaxial compression.

Reference (Timoshenko & Gere, *Theory of Elastic Stability*, 2nd ed., §9.2, Kirchhoff plate):

    N_cr = k pi^2 D / b^2,   k = min_m (m b/a + a/(m b))^2,   D = E h^3 / (12 (1 - nu^2))

with ``a`` the loaded length (along the compression) and ``b`` the loaded edge. Before the FE result is
read, the reference is checked twice (standing rule L20 of ``odd/tasks/static-dtube-before-blade.md``):

  A. against Timoshenko & Gere's tabulated k (square plate k = 4.00, a/b = 1.5 k = 4.34, ...);
  B. against an independent finite-difference solve of D lap^2 w + N w_xx = 0 with w = lap w = 0 on
     the edges, which never uses the closed form.

FE model: steel plate, b/h = 100 (the Mindlin shear correction is pi^2 D (1/a^2 + 1/b^2)/(kappa G h)
= 0.06 %, below every tolerance here). Hard simple support: w = 0 on every edge and the rotation about
the edge tangent fixed (theta_x on x = const, theta_y on y = const); u = 0 on x = 0, v = 0 on y = 0 and
theta_z everywhere remove the rigid and drilling modes without loading the bending problem. The
membrane stress is prescribed (not solved): sigma = -sigma0 (Pa) in every element through the public
``assemble_geometric_stiffness(stress_field=...)``, so K phi = lam (-K_G) phi gives N_cr = lam sigma0 h.

Accepted tolerances, fixed before the measurement:
  * the relative error decreases monotonically over 8 -> 16 -> 32 elements per side;
  * the observed order from 8/16/32 lies in [1.5, 2.5] (bilinear element, eigenvalue error O(h^2));
  * the Richardson extrapolation from 16/32 is within 0.3 % of N_cr; the 32 x 32 value within 1 %.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.linalg as sla
import scipy.sparse as sp

pytest.importorskip("petsc4py", reason="PETSc not available")
from scipy.sparse.linalg import LinearOperator, eigs, splu  # noqa: E402

from aeroelast.core.assembler import MeshAssembler  # noqa: E402
from aeroelast.core.material import IsotropicMaterial  # noqa: E402
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.core.properties import ShellProperty  # noqa: E402
from aeroelast.elements import ElementFamily  # noqa: E402

E, NU, B_EDGE, H = 210e9, 0.3, 1.0, 0.01
D_PLATE = E * H**3 / (12.0 * (1.0 - NU**2))
SIGMA0 = 1.0e6  # Pa, the prescribed compressive membrane stress


# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------


def ss_buckling_coefficient(aspect: float, m_max: int = 20) -> tuple[float, int]:
    """k and the number of half-waves m along the load for a/b = ``aspect``."""
    m = np.arange(1, m_max + 1)
    k = (m / aspect + aspect / m) ** 2
    i = int(np.argmin(k))
    return float(k[i]), int(m[i])


def ss_plate_ncr(a: float, b: float, d: float) -> float:
    """Kirchhoff SS plate, uniaxial compression along ``a``: N_cr (N/m)."""
    return ss_buckling_coefficient(a / b)[0] * np.pi**2 * d / b**2


def fd_plate_ncr(a: float, b: float, d: float, n: int) -> float:
    """Independent check B: FD solve of D lap^2 w = N (-w_xx), w = lap w = 0 on the edges.

    The squared Dirichlet 5-point Laplacian imposes exactly the SS conditions w = lap w = 0, so the
    smallest generalized eigenvalue converges to N_cr at O(dx^2) without using the closed form.
    """
    nx = max(int(round(n * a / b)), 2)
    dx, dy = a / (nx + 1), b / (n + 1)

    def second_difference(k, step):
        return sp.diags([np.ones(k - 1), -2.0 * np.ones(k), np.ones(k - 1)], [-1, 0, 1]) / step**2

    lxx = sp.kron(sp.eye(n), second_difference(nx, dx))
    lap = lxx + sp.kron(second_difference(n, dy), sp.eye(nx))
    stiff = (d * (lap @ lap)).toarray()
    geo = (-lxx).toarray()
    return float(sla.eigh(stiff, geo, eigvals_only=True, subset_by_index=[0, 0])[0])


# ---------------------------------------------------------------------------
# FE model
# ---------------------------------------------------------------------------


def _plate_mesh(nx: int, ny: int, a: float, b: float) -> MeshModel:
    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()
    grid: dict[tuple[int, int], Node] = {}
    for j, y in enumerate(np.linspace(0.0, b, ny + 1)):
        for i, x in enumerate(np.linspace(0.0, a, nx + 1)):
            grid[(i, j)] = Node([float(x), float(y), 0.0], geometric_node=False)
            mesh.add_node(grid[(i, j)])
    for j in range(ny):
        for i in range(nx):
            mesh.add_element(MeshElement(
                nodes=[grid[(i, j)], grid[(i + 1, j)], grid[(i + 1, j + 1)], grid[(i, j + 1)]],
                element_type=ElementType.quad,
            ))
    mesh.add_element_set(ElementSet("plate", set(mesh.elements)))
    return mesh


def _to_scipy(mat) -> sp.csr_matrix:
    indptr, indices, data = mat.getValuesCSR()
    return sp.csr_matrix((data, indices, indptr), shape=mat.getSize())


def mitc4_plate_buckling(nx: int, ny: int, a: float, b: float, load: str = "xx", n_modes: int = 4):
    """Lowest buckling factors and mode shapes (w at the nodes) of the SS plate, N_cr = lam sigma0 h."""
    mesh = _plate_mesh(nx, ny, a, b)
    prop = ShellProperty(material=IsotropicMaterial(name="steel", E=E, nu=NU, rho=7800.0), thickness=H)
    asm = MeshAssembler(mesh=mesh, model={"elements": {"element_family": ElementFamily.SHELL,
                                                       "properties": {"plate": prop}}})
    dpn, tol = asm.dofs_per_node, 1e-9 * max(a, b)
    if dpn != 6:
        raise AssertionError(f"SHELL family expected 6 dofs per node, got {dpn}")

    fixed = set()
    xyz = np.zeros((mesh.node_count, 3))
    for node in asm.nodes:
        i = asm.node_id_to_index[node.id]
        x, y, _ = xyz[i] = node.coords
        on_x = abs(x) < tol or abs(x - a) < tol
        on_y = abs(y) < tol or abs(y - b) < tol
        dofs = {5}  # theta_z: drilling
        if on_x or on_y:
            dofs.add(2)  # w
        if on_x:
            dofs.add(3)  # theta_x = rotation about the edge tangent y
        if on_y:
            dofs.add(4)  # theta_y = rotation about the edge tangent x
        if abs(x) < tol:
            dofs.add(0)
        if abs(y) < tol:
            dofs.add(1)
        fixed.update(i * dpn + d for d in dofs)
    free = np.array(sorted(set(range(asm.dofs_count)) - fixed), dtype=np.int64)

    voigt = {"xx": [-SIGMA0, 0.0, 0.0], "yy": [0.0, -SIGMA0, 0.0]}[load]
    stress = {el.id: np.array(voigt) for el in asm.elements}
    K = _to_scipy(asm.assemble_stiffness_matrix())[free][:, free].tocsc()
    Gf = -_to_scipy(asm.assemble_geometric_stiffness(stress_field=stress))[free][:, free].tocsr()

    # K phi = lam Gf phi; Gf is PSD under pure compression, so 1/lam are the largest eigenvalues of
    # K^-1 Gf (the Arnoldi route of tools/diagnose_rotating_prestress.py buckle, shift 0).
    lu = splu(K)
    nu, vec = eigs(LinearOperator(K.shape, matvec=lambda v: lu.solve(Gf @ v)), k=n_modes, which="LR",
                   tol=1e-12, v0=np.random.default_rng(0).normal(size=len(free)))
    order = np.argsort(-nu.real)
    lam, vec = 1.0 / nu.real[order], np.real(vec[:, order])
    residual = [np.linalg.norm(K @ v - lv * (Gf @ v)) / np.linalg.norm(K @ v) for lv, v in zip(lam, vec.T, strict=True)]
    w = np.zeros((n_modes, mesh.node_count))
    for k in range(n_modes):
        full = np.zeros(asm.dofs_count)
        full[free] = vec[:, k]
        w[k] = full[2::dpn]
    return lam, np.asarray(residual), w, xyz, (K, Gf)


def half_waves(w: np.ndarray, xyz: np.ndarray, axis: int, at: float) -> int:
    """Half-waves of w along ``axis`` on the grid line through the other coordinate ``at``."""
    other = 1 - axis
    line = np.abs(xyz[:, other] - at) < 1e-9
    prof = w[line][np.argsort(xyz[line, axis])][1:-1]
    prof = prof[np.abs(prof) > 1e-6 * np.abs(prof).max()]
    return int(np.count_nonzero(np.diff(np.sign(prof)))) + 1


# ---------------------------------------------------------------------------
# Reference checks (L20: twice, before any FE comparison)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("aspect", "k_table", "m_table"),
    [(0.5, 6.25, 1), (0.6, 5.14, 1), (0.8, 4.20, 1), (1.0, 4.00, 1), (1.5, 4.34, 2), (2.0, 4.00, 2)],
)
def test_reference_check_a_timoshenko_gere_table(aspect, k_table, m_table):
    k, m = ss_buckling_coefficient(aspect)
    assert m == m_table
    assert k == pytest.approx(k_table, abs=0.005)  # the table is printed to 3 significant figures


@pytest.mark.parametrize("aspect", [1.0, 1.5])
def test_reference_check_b_finite_difference_pde(aspect):
    exact = ss_plate_ncr(aspect * B_EDGE, B_EDGE, D_PLATE)
    errs = [fd_plate_ncr(aspect * B_EDGE, B_EDGE, D_PLATE, n) / exact - 1.0 for n in (15, 31)]
    assert abs(errs[1]) < 2e-3
    assert abs(errs[1]) < abs(errs[0]) / 3.0  # O(dx^2): halving dx cuts the error ~4x


# ---------------------------------------------------------------------------
# MITC4 vs the reference
# ---------------------------------------------------------------------------


def test_mitc4_square_plate_converges_to_4pi2d_over_b2():
    exact = ss_plate_ncr(B_EDGE, B_EDGE, D_PLATE)
    ncr, rows = {}, []
    for n in (4, 8, 16, 32):
        lam, res, w, xyz, _ = mitc4_plate_buckling(n, n, B_EDGE, B_EDGE)
        ncr[n] = lam[0] * SIGMA0 * H
        rows.append(f"n={n:2d} N_cr={ncr[n]:.6e} err={ncr[n] / exact - 1.0:+.4%} res={res[0]:.1e}")
        assert res[0] < 1e-8
        if n >= 8:
            assert half_waves(w[0], xyz, 0, B_EDGE / 2) == 1
            assert half_waves(w[0], xyz, 1, B_EDGE / 2) == 1
    print(f"\nexact N_cr = {exact:.6e} N/m\n" + "\n".join(rows))
    err = {n: abs(ncr[n] / exact - 1.0) for n in ncr}
    assert err[32] < err[16] < err[8]
    order = np.log2(abs(ncr[8] - ncr[16]) / abs(ncr[16] - ncr[32]))
    richardson = ncr[32] + (ncr[32] - ncr[16]) / (2.0**order - 1.0)
    print(f"observed order {order:.3f}, Richardson {richardson:.6e} ({richardson / exact - 1.0:+.4%})")
    assert 1.5 <= order <= 2.5
    assert abs(richardson / exact - 1.0) < 3e-3
    assert err[32] < 1e-2


def test_mitc4_square_plate_dense_eigensolve_agrees_with_arnoldi():
    lam, _, _, _, (K, Gf) = mitc4_plate_buckling(8, 8, B_EDGE, B_EDGE)
    dense = sla.eigh(Gf.toarray(), K.toarray(), eigvals_only=True)
    assert 1.0 / dense.max() == pytest.approx(lam[0], rel=1e-9)


def test_mitc4_square_plate_load_direction_is_frame_independent():
    lam_x = mitc4_plate_buckling(16, 16, B_EDGE, B_EDGE, "xx")[0][0]
    lam_y = mitc4_plate_buckling(16, 16, B_EDGE, B_EDGE, "yy")[0][0]
    assert lam_y == pytest.approx(lam_x, rel=1e-8)


def test_mitc4_rectangular_plate_picks_two_half_waves():
    a = 1.5 * B_EDGE
    exact = ss_plate_ncr(a, B_EDGE, D_PLATE)
    lam, res, w, xyz, _ = mitc4_plate_buckling(48, 32, a, B_EDGE)
    ncr = lam[0] * SIGMA0 * H
    print(f"\na/b=1.5 exact {exact:.6e}, MITC4 48x32 {ncr:.6e} ({ncr / exact - 1.0:+.4%})")
    assert res[0] < 1e-8
    assert half_waves(w[0], xyz, 0, B_EDGE / 2) == 2
    assert half_waves(w[0], xyz, 1, a / 4) == 1  # x = a/2 is a nodal line of the m = 2 mode
    assert abs(ncr / exact - 1.0) < 1e-2
