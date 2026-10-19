"""T8 anchor: centrifugal K_G prestress (``compute_prestress_field``) vs a closed form.

Problem. A flat strip in the xy plane, root clamped at x = r0 = 2 m, free tip at R = 10 m,
width b = 0.5 m, spins about the global y axis through the origin. The centrifugal body force
is m omega^2 x along +x (m = mass per area); it lies in the strip plane, so nothing loads the
strip transversely.

Closed form. Equilibrium alone fixes the membrane resultant (statically determinate):

    N_x(x) = m omega^2 (R^2 - x^2) / 2,    N_y = N_xy = 0,    M_x = 0 (no transverse load).

With nu = 0 the uniaxial state is compatible (eps_y = gamma_xy = 0), so it is the exact plane
elasticity solution, not a beam approximation; this is why every material here has nu = 0.
Because M_x = 0, an unsymmetric laminate (B != 0) bends, kappa_x = -(B11/D11) eps_x, yet
N_x is unchanged. ``compute_prestress_field`` returns (A eps + B kappa)/h = N/h, so it must
equal N_x/h for any layup, while the output membrane stress A eps/h (``compute_stress_field``)
exceeds it by exactly A11 D11 / (A11 D11 - B11^2). That ratio is the B kappa term fixed in
41e9a16; for the unsymmetric case below it is 2.0083, far above every tolerance, so the
pre-41e9a16 A eps/h recovery would read the tension 100 % high.

Discretisation, derived before measuring. ``centrifugal_load`` puts each element's load
m omega^2 x_c A_e at its centroid radius x_c, split equally over its nodes. Summing the loads
outboard of element i (midpoint rule, exact for the linear integrand) gives the element
resultant N_i = N_x(x_c) - m omega^2 L^2 / 8, L the element length. Predeclared tolerances:

1. discrete identity: |sigma_FE + m omega^2 L^2/(8h) - sigma(x_c)| <= 1e-6 sigma_root;
2. physical metric: max |sigma_FE - sigma(x_c)| / sigma_root <= L^2 / (4 (R^2 - r0^2))
   (1.04e-4 at 40 elements), and second-order convergence: error ratio 40 vs 20 elements
   within 4 +/- 1 %;
3. sigma_yy and sigma_xy: |.| <= 1e-6 sigma_root;
4. B kappa: sigma_out / sigma_prestress = A11 D11 / (A11 D11 - B11^2) within 1e-6 relative;
   it is 1 for the symmetric layups.

Coverage: MITC4 quads only. MITC3 triangles split loads unevenly across the width, so the 1D
closed form is not exact for them; that path is not anchored here.
"""

import numpy as np
import pytest

from aeroelast.core.assembler import MeshAssembler
from aeroelast.core.laminate import Laminate, Ply
from aeroelast.core.material import IsotropicMaterial, OrthotropicMaterial
from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.core.properties import CompositeShellProperty, ShellProperty
from aeroelast.elements import ElementFamily

R0, R_TIP, WIDTH, OMEGA = 2.0, 10.0, 0.5, 1.2
AXIS, CENTER = np.array([0.0, 1.0, 0.0]), np.zeros(3)
IDENTITY_TOL = 1e-6
TRANSVERSE_TOL = 1e-6
RATIO_TOL = 1e-6


def _orthotropic(name, e, rho):
    g = e / 2.0
    return OrthotropicMaterial(name=name, E=(e, e, e), G=(g, g, g), nu=(0.0, 0.0, 0.0), rho=rho)


def _section(layup):
    """Return (property, thickness, mass per area, A11, B11, D11 or None)."""
    if layup == "isotropic":
        h, rho = 0.01, 7850.0
        prop = ShellProperty(material=IsotropicMaterial(name="steel", E=210e9, nu=0.0, rho=rho),
                             thickness=h)
        return prop, h, rho * h, None
    stiff, soft = _orthotropic("stiff", 100e9, 2000.0), _orthotropic("soft", 10e9, 1200.0)
    plies = ([Ply(stiff, 0.0025, 0), Ply(soft, 0.005, 0), Ply(stiff, 0.0025, 0)]
             if layup == "symmetric" else [Ply(soft, 0.005, 0), Ply(stiff, 0.005, 0)])
    lam = Laminate(plies)
    mass = sum(p.material.rho * p.thickness for p in plies)
    return CompositeShellProperty(laminate=lam), lam.total_thickness, mass, lam


def _strip(layup, n_span, n_width=2):
    prop, h, mass, lam = _section(layup)
    mesh = MeshModel()
    nodes = {}
    for j in range(n_width + 1):
        for i in range(n_span + 1):
            x = R0 + (R_TIP - R0) * i / n_span
            y = WIDTH * (j / n_width - 0.5)
            nodes[i, j] = Node([x, y, 0.0], geometric_node=False)
            mesh.add_node(nodes[i, j])
    for j in range(n_width):
        for i in range(n_span):
            mesh.add_element(MeshElement(
                nodes=[nodes[i, j], nodes[i + 1, j], nodes[i + 1, j + 1], nodes[i, j + 1]],
                element_type=ElementType.quad))
    mesh.add_element_set(ElementSet("strip", set(mesh.elements)))
    asm = MeshAssembler(mesh=mesh, model={"span_direction": [1, 0, 0], "elements": {
        "element_family": ElementFamily.SHELL, "properties": {"strip": prop}}})
    x_c = {el.id: float(np.mean([n.coords[0] for n in el.nodes])) for el in mesh.elements}
    return asm, h, mass, lam, x_c


def _dense(mat):
    from scipy.sparse import csr_matrix

    indptr, indices, data = mat.getValuesCSR()
    return csr_matrix((data, indices, indptr), shape=mat.getSize()).toarray()


def _solve(layup, n_span):
    asm, h, mass, lam, x_c_by_id = _strip(layup, n_span)
    fixed = []
    for node in asm.nodes:
        if abs(node.coords[0] - R0) < 1e-12:
            i = asm.node_id_to_index[node.id]
            fixed.extend(range(6 * i, 6 * i + 6))
    free = np.setdiff1d(np.arange(asm.dofs_count), fixed)
    k = _dense(asm.assemble_stiffness_matrix())
    f = asm.assemble_centrifugal_load(OMEGA, AXIS, CENTER)
    u = np.zeros(asm.dofs_count)
    u[free] = np.linalg.solve(k[np.ix_(free, free)], f[free])
    x_c = np.array([x_c_by_id[el.id] for el in asm.elements])
    prestress = np.asarray(asm._rust.compute_prestress_field(u))
    output, _ = asm._rust.compute_stress_field(u, 0.0, 0)
    return {"x_c": x_c, "prestress": prestress, "output": np.asarray(output)[:, [0, 1, 3]],
            "h": h, "mass": mass, "lam": lam, "length": (R_TIP - R0) / n_span, "f": f}


def _sigma_exact(x, mass, h):
    return mass * OMEGA**2 * (R_TIP**2 - x**2) / (2.0 * h)


@pytest.mark.parametrize("layup", ["isotropic", "symmetric", "unsymmetric"])
def test_prestress_matches_rotating_strip_closed_form(layup):
    errors = {}
    for n_span in (20, 40):
        r = _solve(layup, n_span)
        sig_root = _sigma_exact(R0, r["mass"], r["h"])
        exact = _sigma_exact(r["x_c"], r["mass"], r["h"])
        sxx = r["prestress"][:, 0]
        # The load really is the closed-form body force: total = m w^2 b (R^2 - r0^2)/2.
        total = r["f"][0::6].sum()
        assert total == pytest.approx(r["mass"] * OMEGA**2 * WIDTH * (R_TIP**2 - R0**2) / 2,
                                      rel=1e-12)
        bias = r["mass"] * OMEGA**2 * r["length"] ** 2 / (8.0 * r["h"])
        identity = np.max(np.abs(sxx + bias - exact)) / sig_root
        assert identity <= IDENTITY_TOL, f"{layup} n={n_span}: discrete identity {identity:.3e}"
        errors[n_span] = np.max(np.abs(sxx - exact)) / sig_root
        bound = r["length"] ** 2 / (4.0 * (R_TIP**2 - R0**2))
        assert errors[n_span] <= bound * (1 + IDENTITY_TOL), (
            f"{layup} n={n_span}: error {errors[n_span]:.3e} > bound {bound:.3e}")
        transverse = np.max(np.abs(r["prestress"][:, 1:])) / sig_root
        assert transverse <= TRANSVERSE_TOL, f"{layup}: sigma_yy/xy {transverse:.3e}"
    assert errors[20] / errors[40] == pytest.approx(4.0, rel=0.01)


@pytest.mark.parametrize("layup", ["symmetric", "unsymmetric"])
def test_output_stress_differs_from_prestress_by_the_b_kappa_factor(layup):
    r = _solve(layup, 40)
    lam = r["lam"]
    a11, b11, d11 = lam.A[0, 0], lam.B[0, 0], lam.D[0, 0]
    factor = a11 * d11 / (a11 * d11 - b11**2)
    ratio = r["output"][:, 0] / r["prestress"][:, 0]
    np.testing.assert_allclose(ratio, factor, rtol=RATIO_TOL)
    if layup == "symmetric":
        assert factor == pytest.approx(1.0, abs=1e-12)
    else:
        assert factor - 1.0 > 1e3 * RATIO_TOL  # B kappa is not a rounding-level effect
