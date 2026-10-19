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
from aeroelast.core.laminate import Laminate, Ply  # noqa: E402
from aeroelast.core.material import IsotropicMaterial, OrthotropicMaterial  # noqa: E402
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node  # noqa: E402
from aeroelast.core.mesh.model import MeshModel  # noqa: E402
from aeroelast.core.properties import CompositeShellProperty, ShellProperty  # noqa: E402
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


def mitc4_plate_buckling(nx: int, ny: int, a: float, b: float, load: str = "xx", n_modes: int = 4, prop=None):
    """Lowest buckling factors and mode shapes (w at the nodes) of the SS plate, N_cr = lam sigma0 h.

    ``prop`` defaults to the isotropic steel plate of thickness ``H``; a ``CompositeShellProperty``
    gives the laminate case, with h its total thickness.
    """
    mesh = _plate_mesh(nx, ny, a, b)
    if prop is None:
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


# ---------------------------------------------------------------------------
# T7 (iii): soft-core sandwich, first-order shear deformation (FSDT) reference
# ---------------------------------------------------------------------------
#
# Symmetric sandwich, isotropic faces and core with one Poisson ratio, so D11 = D22 = D, D12 = nu D and
# the section is an isotropic FSDT plate with bending stiffness D (CLT) and transverse-shear stiffness S.
# SS plate, uniaxial N along x, w = sin(kx x) sin(ky y), kx = m pi/a, ky = pi/b, K^2 = kx^2 + ky^2:
#
#     N_cr = min_m  D K^4 / kx^2 / (1 + D K^2 / S)
#
# (Engesser form; for a wide plate it is Allen's strut 1/N = 1/N_E + 1/S). Checked twice: A, a numerical
# Navier solve of the 3x3 FSDT energy Hessian in (w, psi_x, psi_y); B, its limits (S -> inf gives
# Kirchhoff, b -> inf gives Allen) and S itself from two definitions, the energy-equivalent
# S = D^2 / int Phi^2/G dz (own quadrature, the definition MITC4 consumes from the laminate) and Allen's
# antiplane-core S = G_c d^2/t_c. Accepted FE tolerances (fixed before measuring): b = 1.5 m (m = 1,
# shear lowers N 25 %) monotone over 8/16/32, order in [1.5, 2.5], Richardson within 0.5 %, 32x32 within
# 1 %; b = 0.75 m (m = 2, shear lowers N 2.25x) within 1 % at 48x48 with m = 2 picked. The S spread
# (0.16 %) moves N by < 0.1 %.

E_FACE, T_FACE, E_CORE, T_CORE = 20e9, 0.002, 60e6, 0.030
SANDWICH = ((E_FACE, T_FACE), (E_CORE, T_CORE), (E_FACE, T_FACE))
H_SANDWICH = sum(t for _, t in SANDWICH)


def sandwich_section(n_quad: int = 200001) -> tuple[float, float, float]:
    """(D, S_energy, S_allen) of SANDWICH, own implementation (independent of aeroelast.core.laminate)."""
    q = lambda e: e / (1.0 - NU**2)  # noqa: E731
    zb = np.cumsum([0.0] + [t for _, t in SANDWICH]) - H_SANDWICH / 2.0
    d = sum(q(e) * (zb[i + 1] ** 3 - zb[i] ** 3) / 3.0 for i, (e, _) in enumerate(SANDWICH))
    z = np.linspace(-H_SANDWICH / 2.0, H_SANDWICH / 2.0, n_quad)
    e_z = np.where(np.abs(z) < T_CORE / 2.0, E_CORE, E_FACE)
    integrand = q(e_z) * z
    phi = -np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(z))])
    s_energy = d**2 / np.trapezoid(phi**2 / (e_z / (2.0 * (1.0 + NU))), z)
    g_core = E_CORE / (2.0 * (1.0 + NU))
    return float(d), float(s_energy), g_core * (T_CORE + T_FACE) ** 2 / T_CORE


def fsdt_ncr_mode(a: float, b: float, d: float, s: float, m: int) -> float:
    kx, ky = m * np.pi / a, np.pi / b
    k2 = kx**2 + ky**2
    return d * k2**2 / kx**2 / (1.0 + d * k2 / s)


def fsdt_ncr(a: float, b: float, d: float, s: float, m_max: int = 20) -> tuple[float, int]:
    values = [fsdt_ncr_mode(a, b, d, s, m) for m in range(1, m_max + 1)]
    i = int(np.argmin(values))
    return values[i], i + 1


def fsdt_navier_ncr(a: float, b: float, d: float, s: float, m: int, n: float = 1.0) -> float:
    """Check A: smallest N of the FSDT energy Hessian in (W, X, Y) for one (m, n) Navier term."""
    kx, ky = m * np.pi / a, n * np.pi / b
    c = (1.0 - NU) / 2.0
    k_bend = d * np.array([[0.0, 0.0, 0.0],
                           [0.0, kx**2 + c * ky**2, (NU + c) * kx * ky],
                           [0.0, (NU + c) * kx * ky, ky**2 + c * kx**2]])
    k_shear = s * np.array([[kx**2 + ky**2, kx, ky], [kx, 1.0, 0.0], [ky, 0.0, 1.0]])
    k_geo = np.diag([kx**2, 0.0, 0.0])
    return float(1.0 / np.linalg.eigvals(np.linalg.solve(k_bend + k_shear, k_geo)).real.max())


def sandwich_property() -> CompositeShellProperty:
    def iso(name, e):
        g = e / (2.0 * (1.0 + NU))
        return OrthotropicMaterial(name=name, E=(e, e, e), G=(g, g, g), nu=(NU, NU, NU), rho=1000.0)

    face, core = iso("face", E_FACE), iso("core", E_CORE)
    plies = [Ply(face, T_FACE, 0.0), Ply(core, T_CORE, 0.0), Ply(face, T_FACE, 0.0)]
    return CompositeShellProperty(laminate=Laminate(plies))


@pytest.mark.parametrize("b", [0.75, 1.5])
def test_sandwich_reference_check_a_navier_fsdt(b):
    d, s, _ = sandwich_section()
    for m in range(1, 5):
        assert fsdt_navier_ncr(b, b, d, s, m) == pytest.approx(fsdt_ncr_mode(b, b, d, s, m), rel=1e-12)


def test_sandwich_reference_check_b_limits_and_shear_stiffness():
    d, s_energy, s_allen = sandwich_section()
    assert s_energy / s_allen - 1.0 == pytest.approx(0.0, abs=5e-3)  # measured +0.16 %
    for b in (0.75, 1.5):  # S -> inf: Kirchhoff (the reference already checked in A/B above)
        k, _ = ss_buckling_coefficient(1.0)
        assert fsdt_ncr(b, b, d, 1e20)[0] == pytest.approx(k * np.pi**2 * d / b**2, rel=1e-10)
    kx = np.pi  # a = 1, b -> inf: Allen's sandwich strut
    allen = 1.0 / (1.0 / (d * kx**2) + 1.0 / s_allen)
    assert fsdt_navier_ncr(1.0, 1e6, d, s_allen, 1, n=1e-9) == pytest.approx(allen, rel=1e-9)
    assert fsdt_ncr(1.5, 1.5, d, s_energy)[1] == 1
    assert fsdt_ncr(0.75, 0.75, d, s_energy)[1] == 2


def test_sandwich_laminate_section_matches_reference():
    """The input MITC4 consumes: aeroelast's CLT D and energy-equivalent Cs vs the own implementation."""
    d, s_energy, _ = sandwich_section()
    lam = sandwich_property().laminate
    assert lam.total_thickness == pytest.approx(H_SANDWICH, rel=1e-12)
    assert lam.D[0, 0] == pytest.approx(d, rel=1e-6)
    assert lam.Cs[0, 0] == pytest.approx(s_energy, rel=1e-4)
    assert lam.Cs[1, 1] == pytest.approx(s_energy, rel=1e-4)


_SECTION_SHEAR_DEFECT = pytest.mark.xfail(
    strict=True,
    reason="T7(iii) L17: the MITC4 composite shear block consumes A55 = sum G_k t_k "
    "(crates/aeroelast-py/src/assembler.rs, ADR-1 amended, 8cbfc0b), 40x the sandwich S; "
    "MITC4 then equals FSDT with S = A55 (test_mitc4_sandwich_equals_fsdt_with_the_shear_it_consumes)",
)


@_SECTION_SHEAR_DEFECT
def test_mitc4_sandwich_plate_m1_converges_to_fsdt():
    b = 1.5
    d, s, _ = sandwich_section()
    exact, _ = fsdt_ncr(b, b, d, s)
    ncr, rows = {}, []
    for n in (8, 16, 32):
        lam, res, w, xyz, _ = mitc4_plate_buckling(n, n, b, b, prop=sandwich_property())
        ncr[n] = lam[0] * SIGMA0 * H_SANDWICH
        rows.append(f"n={n:2d} N_cr={ncr[n]:.6e} err={ncr[n] / exact - 1.0:+.4%} res={res[0]:.1e} "
                    f"m={half_waves(w[0], xyz, 0, b / 2)}")
        assert res[0] < 1e-8
    print(f"\nsandwich b={b}: FSDT N_cr = {exact:.6e} N/m (Kirchhoff {4 * np.pi**2 * d / b**2:.6e})\n"
          + "\n".join(rows))
    assert half_waves(w[0], xyz, 0, b / 2) == 1
    err = {n: abs(ncr[n] / exact - 1.0) for n in ncr}
    assert err[32] < err[16] < err[8]
    order = np.log2(abs(ncr[8] - ncr[16]) / abs(ncr[16] - ncr[32]))
    richardson = ncr[32] + (ncr[32] - ncr[16]) / (2.0**order - 1.0)
    print(f"observed order {order:.3f}, Richardson {richardson:.6e} ({richardson / exact - 1.0:+.4%})")
    assert 1.5 <= order <= 2.5
    assert abs(richardson / exact - 1.0) < 5e-3
    assert err[32] < 1e-2


@_SECTION_SHEAR_DEFECT
def test_mitc4_sandwich_plate_m2_shear_dominated():
    b = 0.75
    d, s, _ = sandwich_section()
    exact, m = fsdt_ncr(b, b, d, s)
    lam, res, w, xyz, _ = mitc4_plate_buckling(48, 48, b, b, prop=sandwich_property())
    ncr = lam[0] * SIGMA0 * H_SANDWICH
    print(f"\nsandwich b={b}: FSDT N_cr = {exact:.6e} (m={m}), MITC4 48x48 {ncr:.6e} "
          f"({ncr / exact - 1.0:+.4%}), m_FE={half_waves(w[0], xyz, 0, b / 4)}, lam2/lam1={lam[1] / lam[0]:.4f}")
    assert res[0] < 1e-8
    assert half_waves(w[0], xyz, 0, b / 4) == 2
    assert abs(ncr / exact - 1.0) < 1e-2


@pytest.mark.parametrize(("b", "n", "m_expected"), [(1.5, 32, 1), (0.75, 48, 1)])
def test_mitc4_sandwich_equals_fsdt_with_the_shear_it_consumes(b, n, m_expected):
    """Characterization of the L17 defect: the element is a correct FSDT plate for the S it is given.

    With S = A55 (``Laminate.shear_stiffness_uncorrected``) the same reference matches MITC4 (Richardson
    +0.002 % at b = 1.5, +0.06 % at 48x48 for b = 0.75), so the error is the section input, not the
    element. Delete this test when the composite shear block receives the energy-equivalent ``Cs``: the
    two xfail tests above then pass and this one fails.
    """
    d, _, _ = sandwich_section()
    a55 = float(sandwich_property().laminate.shear_stiffness_uncorrected()[0, 0])
    expected, m = fsdt_ncr(b, b, d, a55)
    lam, res, w, xyz, _ = mitc4_plate_buckling(n, n, b, b, prop=sandwich_property())
    assert m == m_expected
    assert res[0] < 1e-8
    assert half_waves(w[0], xyz, 0, b / 4) == m_expected
    assert lam[0] * SIGMA0 * H_SANDWICH / expected - 1.0 == pytest.approx(0.0, abs=5e-3)


# ---------------------------------------------------------------------------
# T7 (ii): specially orthotropic cross-ply plate, Kirchhoff (CLT) reference
# ---------------------------------------------------------------------------
#
# Symmetric cross-ply [0/90/90/0] (B = 0, D16 = D26 = 0), SS plate, uniaxial N along the a edge,
# w = sin(m pi x/a) sin(pi y/b) solves D11 w_xxxx + 2 (D12 + 2 D66) w_xxyy + D22 w_yyyy = -N w_xx exactly:
#
#     N_cr = min_m  pi^2/b^2 [D11 (m b/a)^2 + 2 (D12 + 2 D66) + D22 (a/(m b))^2]
#
# (Whitney, Structural Analysis of Laminated Anisotropic Plates, 1987; Jones, Mechanics of Composite
# Materials, 2nd ed., 1999, specially orthotropic plate buckling). The form is checked twice before the FE
# is read: A, it reduces to the isotropic T&G k already checked against the table, its m -> m+1 switch
# sits at a/b = (D11/D22)^(1/4) sqrt(m (m+1)) and b -> inf gives the wide column pi^2 D11/a^2; B, an
# independent FD solve of the orthotropic PDE (squared 1D Dirichlet second differences impose w = 0 and
# M_n = 0 on every edge), which never uses the closed form. The section D is checked against an own CLT.
#
# Carbon/epoxy plies, h = 4 mm, b/h = 250: D11/D22 = 4.7, so the square plate buckles with m = 1 when
# compressed along the outer fibres (x) and with m = 2 along y (switch at a/b = 0.96); the two loads
# therefore pin the ply-angle convention (0 deg = element e1 = global x here). Transverse shear lowers
# N_cr by < 0.1 % for any S the element may consume (bound asserted below with the smaller of the
# energy-equivalent and uncorrected section shear), so the Kirchhoff value is the target.
#
# Accepted FE tolerances, fixed before the measurement (same basis as the isotropic case above):
#   * x load (m = 1): relative error monotone over 8 -> 16 -> 32, observed order in [1.5, 2.5],
#     Richardson from 16/32 within 0.3 % of N_cr, 32 x 32 within 1 %, one half-wave each way;
#   * y load (m = 2): 48 x 48 within 1 % with two half-waves along y picked;
#   * the [90/0/0/90] plate under y load equals the [0/90/90/0] plate under x load to 1e-8.

E1_CF, E2_CF, G12_CF, G23_CF, NU12_CF, T_PLY = 140e9, 10e9, 5e9, 3.5e9, 0.3, 0.001
CROSS_PLY = (0.0, 90.0, 90.0, 0.0)
H_CROSS = T_PLY * len(CROSS_PLY)


def cross_ply_d(angles=CROSS_PLY) -> np.ndarray:
    """(D11, D12, D22, D66) of the cross-ply, own CLT (independent of aeroelast.core.laminate)."""
    den = 1.0 - NU12_CF**2 * E2_CF / E1_CF
    q11, q22, q12 = E1_CF / den, E2_CF / den, NU12_CF * E2_CF / den
    zb = np.arange(len(angles) + 1) * T_PLY - H_CROSS / 2.0
    d = np.zeros(4)
    for i, ang in enumerate(angles):
        a, b = (q11, q22) if ang == 0.0 else (q22, q11)
        d += np.array([a, q12, b, G12_CF]) * (zb[i + 1] ** 3 - zb[i] ** 3) / 3.0
    return d


def ortho_ncr_mode(a: float, b: float, d: np.ndarray, m: int) -> float:
    d11, d12, d22, d66 = d
    r = m * b / a
    return np.pi**2 / b**2 * (d11 * r**2 + 2.0 * (d12 + 2.0 * d66) + d22 / r**2)


def ortho_ncr(a: float, b: float, d: np.ndarray, m_max: int = 20) -> tuple[float, int]:
    values = [ortho_ncr_mode(a, b, d, m) for m in range(1, m_max + 1)]
    i = int(np.argmin(values))
    return values[i], i + 1


def fd_ortho_ncr(a: float, b: float, d: np.ndarray, n: int) -> float:
    """Check B: FD solve of the orthotropic plate PDE with SS edges (w = w_nn = 0)."""
    d11, d12, d22, d66 = d
    nx = max(int(round(n * a / b)), 2)
    dx, dy = a / (nx + 1), b / (n + 1)

    def second_difference(k, step):
        return sp.diags([np.ones(k - 1), -2.0 * np.ones(k), np.ones(k - 1)], [-1, 0, 1]) / step**2

    lxx = sp.kron(sp.eye(n), second_difference(nx, dx))
    lyy = sp.kron(second_difference(n, dy), sp.eye(nx))
    stiff = (d11 * (lxx @ lxx) + 2.0 * (d12 + 2.0 * d66) * (lxx @ lyy) + d22 * (lyy @ lyy)).toarray()
    return float(sla.eigh(stiff, (-lxx).toarray(), eigvals_only=True, subset_by_index=[0, 0])[0])


def fsdt_ortho_ncr_mode(a: float, b: float, d: np.ndarray, s: float, m: int) -> float:
    """Navier FSDT of the specially orthotropic plate, both transverse shear stiffnesses equal to s."""
    d11, d12, d22, d66 = d
    kx, ky = m * np.pi / a, np.pi / b
    k_bend = np.array([[0.0, 0.0, 0.0],
                       [0.0, d11 * kx**2 + d66 * ky**2, (d12 + d66) * kx * ky],
                       [0.0, (d12 + d66) * kx * ky, d22 * ky**2 + d66 * kx**2]])
    k_shear = s * np.array([[kx**2 + ky**2, kx, ky], [kx, 1.0, 0.0], [ky, 0.0, 1.0]])
    k_geo = np.diag([kx**2, 0.0, 0.0])
    return float(1.0 / np.linalg.eigvals(np.linalg.solve(k_bend + k_shear, k_geo)).real.max())


def cross_ply_property(angles=CROSS_PLY) -> CompositeShellProperty:
    cf = OrthotropicMaterial(name="cf", E=(E1_CF, E2_CF, E2_CF), G=(G12_CF, G12_CF, G23_CF),
                             nu=(NU12_CF, NU12_CF, 0.4), rho=1500.0)
    return CompositeShellProperty(laminate=Laminate([Ply(cf, T_PLY, ang) for ang in angles]))


def test_ortho_reference_check_a_limits_and_mode_switch():
    iso = D_PLATE * np.array([1.0, NU, 1.0, (1.0 - NU) / 2.0])
    for aspect in (0.5, 1.0, 1.5, 2.0):
        assert ortho_ncr(aspect, 1.0, iso)[0] == pytest.approx(ss_plate_ncr(aspect, 1.0, D_PLATE), rel=1e-12)
    d = cross_ply_d()
    for m in (1, 2, 3):  # N_m = N_m+1 at a/b = (D11/D22)^(1/4) sqrt(m (m+1))
        switch = (d[0] / d[2]) ** 0.25 * np.sqrt(m * (m + 1))
        assert ortho_ncr_mode(switch, 1.0, d, m) == pytest.approx(ortho_ncr_mode(switch, 1.0, d, m + 1), rel=1e-12)
        assert ortho_ncr(0.99 * switch, 1.0, d)[1] == m
        assert ortho_ncr(1.01 * switch, 1.0, d)[1] == m + 1
    assert ortho_ncr(1.0, 1e4, d)[0] == pytest.approx(np.pi**2 * d[0], rel=1e-6)  # wide column, a = 1
    # The Navier FSDT Hessian with S -> inf is the same closed form. S = 1e11: shear term ~1e-7 and
    # round-off ~1e-8; S = 1e20 makes the 3x3 solve singular in float64.
    for m in (1, 2, 3):
        assert fsdt_ortho_ncr_mode(1.0, 1.0, d, 1e11, m) == pytest.approx(ortho_ncr_mode(1.0, 1.0, d, m), rel=1e-6)
    assert ortho_ncr(1.0, 1.0, d)[1] == 1
    assert ortho_ncr(1.0, 1.0, d[[2, 1, 0, 3]])[1] == 2  # load along y: D11 and D22 swap roles


@pytest.mark.parametrize("load", ["xx", "yy"])
def test_ortho_reference_check_b_finite_difference_pde(load):
    d = cross_ply_d() if load == "xx" else cross_ply_d()[[2, 1, 0, 3]]
    exact = ortho_ncr(1.0, 1.0, d)[0]
    errs = [fd_ortho_ncr(1.0, 1.0, d, n) / exact - 1.0 for n in (15, 31)]
    assert abs(errs[1]) < 2e-3
    assert abs(errs[1]) < abs(errs[0]) / 3.0


def test_cross_ply_laminate_section_matches_reference_and_shear_is_negligible():
    lam = cross_ply_property().laminate
    d11, d12, d22, d66 = cross_ply_d()
    scale = d11
    assert lam.total_thickness == pytest.approx(H_CROSS, rel=1e-12)
    np.testing.assert_allclose(lam.D, [[d11, d12, 0.0], [d12, d22, 0.0], [0.0, 0.0, d66]], rtol=0, atol=1e-9 * scale)
    np.testing.assert_allclose(lam.B, 0.0, atol=1e-9 * float(np.abs(lam.A).max()) * H_CROSS)
    s_min = min(float(np.diag(lam.Cs).min()), float(np.diag(lam.shear_stiffness_uncorrected()).min()))
    d = cross_ply_d()
    for dd in (d, d[[2, 1, 0, 3]]):
        kirchhoff, m = ortho_ncr(1.0, 1.0, dd)
        assert 0.0 < 1.0 - fsdt_ortho_ncr_mode(1.0, 1.0, dd, s_min, m) / kirchhoff < 1e-3


def test_mitc4_cross_ply_square_plate_x_load_converges():
    exact = ortho_ncr(1.0, 1.0, cross_ply_d())[0]
    ncr, rows = {}, []
    for n in (8, 16, 32):
        lam, res, w, xyz, _ = mitc4_plate_buckling(n, n, 1.0, 1.0, "xx", prop=cross_ply_property())
        ncr[n] = lam[0] * SIGMA0 * H_CROSS
        rows.append(f"n={n:2d} N_cr={ncr[n]:.6e} err={ncr[n] / exact - 1.0:+.4%} res={res[0]:.1e}")
        assert res[0] < 1e-8
        assert half_waves(w[0], xyz, 0, 0.5) == 1
        assert half_waves(w[0], xyz, 1, 0.5) == 1
    print(f"\ncross-ply x load: CLT N_cr = {exact:.6e} N/m\n" + "\n".join(rows))
    err = {n: abs(ncr[n] / exact - 1.0) for n in ncr}
    assert err[32] < err[16] < err[8]
    order = np.log2(abs(ncr[8] - ncr[16]) / abs(ncr[16] - ncr[32]))
    richardson = ncr[32] + (ncr[32] - ncr[16]) / (2.0**order - 1.0)
    print(f"observed order {order:.3f}, Richardson {richardson:.6e} ({richardson / exact - 1.0:+.4%})")
    assert 1.5 <= order <= 2.5
    assert abs(richardson / exact - 1.0) < 3e-3
    assert err[32] < 1e-2


def test_mitc4_cross_ply_square_plate_y_load_picks_two_half_waves():
    exact, m = ortho_ncr(1.0, 1.0, cross_ply_d()[[2, 1, 0, 3]])
    lam, res, w, xyz, _ = mitc4_plate_buckling(48, 48, 1.0, 1.0, "yy", prop=cross_ply_property())
    ncr = lam[0] * SIGMA0 * H_CROSS
    print(f"\ncross-ply y load: CLT N_cr = {exact:.6e} (m={m}), MITC4 48x48 {ncr:.6e} ({ncr / exact - 1.0:+.4%}), "
          f"lam2/lam1={lam[1] / lam[0]:.4f}")
    assert m == 2
    assert res[0] < 1e-8
    assert half_waves(w[0], xyz, 1, 0.5) == 2
    assert half_waves(w[0], xyz, 0, 0.25) == 1
    assert abs(ncr / exact - 1.0) < 1e-2


def test_mitc4_cross_ply_ply_angle_follows_the_element_frame():
    lam_x = mitc4_plate_buckling(16, 16, 1.0, 1.0, "xx", prop=cross_ply_property())[0][0]
    lam_y = mitc4_plate_buckling(16, 16, 1.0, 1.0, "yy", prop=cross_ply_property((90.0, 0.0, 0.0, 90.0)))[0][0]
    assert lam_y == pytest.approx(lam_x, rel=1e-8)
