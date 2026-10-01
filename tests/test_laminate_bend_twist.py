"""Laminate bend-twist coupon — MITC4 composite shell vs a CLT plate reference.

Rebuilt 2026-10-01 from GitHub issue #9 (*"bend-twist: MITC4 matches CalculiX but
is 3.4-4.8x the CLT coupon — the coupon reference needs settling"*). The file the
issue cites did not exist in this tree nor anywhere in its git history, so this is
an authored reconstruction that pins down the definitions the original report left
implicit.

Three design rules, each one the failure mode of the original report:

1. **One laminate object.** The ABD the element receives and the ABD the reference
   uses come from the same ``Laminate`` instance, and the property dict handed to
   the assembler is *asserted* against ``lam.A/lam.B/lam.D``
   (``test_abd_is_the_single_source_of_truth``). The original sweep took the
   reference from ``t._laminate()`` and the displacement from
   ``t._solve_strip(...)``; any divergence between those two paths turns the
   reported ratio into a difference between two *laminates*.
2. **A mesh-independent resultant.** The tip load is distributed over every free
   edge node, so refining the mesh does not change the applied work.
3. **The metric is explicit.** Three independent twist metrics are computed from
   the same displacement field (``test_metric_sensitivity``), because a
   least-squares fit over the tip nodes of a 2-element cantilever has no reason to
   equal the twist angle of the strip.

Coupon geometry (the issue's): ``L = 1.0 m`` along X, ``b = 0.1 m`` along Y,
8 plies of 1 mm (``[45,0,0,45]s``) -> ``h = 8 mm``, tip load ``P = 10 N`` in +Z,
mesh ``NX = 2 x NY = 12``. The element is the composite MITC4 (code 44).

Scope of this file: WU1-WU3 of ``odd/tasks/composite-bend-twist-verdict.md``
(the measurement, its mesh convergence and its metric sensitivity). The reference
verdict (WU4) is a separate step: the CLT formulas below are reproduced *to
characterise the issue's reference*, and are deliberately not asserted as truth.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")
from _aeroelast import PyMeshAssembler

from aeroelast.core.laminate import Laminate, create_laminate_from_angles
from aeroelast.core.material import OrthotropicMaterial

# ─────────────────────────────────────────────────────────────────────────────
# The coupon, as the issue defines it
# ─────────────────────────────────────────────────────────────────────────────

L = 1.0  # strip length along X [m]
B = 0.1  # strip width along Y [m]
PLY_T = 1.0e-3  # ply thickness [m]  (8 plies -> h = 8 mm, the issue's thickness)
P_TIP = 10.0  # total tip load in +Z [N]

NX0, NY0 = 2, 12  # the issue's mesh

# CFRP, the repo's standard coupon material (tests/test_material_suite.py)
_CFRP = OrthotropicMaterial(
    "CFRP", E=(120e9, 10e9, 10e9), G=(5e9, 3e9, 5e9), nu=(0.3, 0.3, 0.3), rho=1500.0
)
# A second, more anisotropic reference material: a conclusion that only holds for
# one material is not a conclusion about the element.
_T300 = OrthotropicMaterial(
    "T300/5208", E=(181e9, 10.3e9, 10.3e9), G=(7.17e9, 3.0e9, 7.17e9),
    nu=(0.28, 0.28, 0.28), rho=1600.0,
)

# The issue's sweep, plus the D16 = 0 control. Each entry is the *half* stack; the
# laminate is its symmetric expansion (the issue's trailing "s").
LAYUPS: dict[str, list[float]] = {
    "[45,0,0,45]s": [45.0, 0.0, 0.0, 45.0],
    "[45,-45,-45,45]s": [45.0, -45.0, -45.0, 45.0],
    "[0,45,45,0]s": [0.0, 45.0, 45.0, 0.0],
    "[45,45,45,45]s": [45.0, 45.0, 45.0, 45.0],
    "[30,0,0,30]s": [30.0, 0.0, 0.0, 30.0],
    "[60,0,0,60]s": [60.0, 0.0, 0.0, 60.0],
    "[0,0,0,0]": [0.0, 0.0, 0.0, 0.0],
}


def _laminate(half: list[float], material: OrthotropicMaterial = _CFRP,
              ply_t: float = PLY_T) -> Laminate:
    """The single source of truth: a symmetric expansion of ``half``.

    ``half + reversed(half)`` is the explicit symmetric stack, so the ply count
    and the stacking order are visible at the call site instead of implied by a
    suffix.
    """
    angles = [float(a) for a in half] + [float(a) for a in reversed(half)]
    return create_laminate_from_angles(material, ply_t, angles)


def _lam_prop(lam: Laminate) -> dict:
    """Composite property dict for ``PyMeshAssembler`` (same keys as the suite).

    ``cs_uncorrected`` is what the element consumes (ADR-1): the section integral
    without the 5/6 of a homogeneous stack.
    """
    h = lam.total_thickness
    mpa = sum(p.material.rho * p.thickness for p in lam.plies)
    ri = sum(p.material.rho * (p.z_top**3 - p.z_bottom**3) / 3 for p in lam.plies)
    return {
        "type": "composite",
        "cm": lam.A.ravel().tolist(),
        "b_coupling": lam.B.ravel().tolist(),
        "cb": lam.D.ravel().tolist(),
        "cs": lam.Cs.ravel().tolist(),
        "cs_uncorrected": lam.shear_stiffness_uncorrected().ravel().tolist(),
        "thickness": h,
        "e_equiv": float(np.trace(lam.A) / (3.0 * h)),
        "mass_per_area": mpa,
        "rotational_inertia": ri,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Mesh, solve, metrics
# ─────────────────────────────────────────────────────────────────────────────


def _strip_mesh(nx: int = NX0, ny: int = NY0, length: float = L, width: float = B):
    """Flat MITC4 strip in the X-Y plane, clamped at ``x = 0``, free at ``x = L``.

    X is the length, Y the width (centred on 0 so the twist metrics are
    sign-symmetric), Z the out-of-plane direction the tip load acts along.
    """
    xs = np.linspace(0.0, length, nx + 1)
    ys = np.linspace(-0.5 * width, 0.5 * width, ny + 1)
    coords = np.array([[x, y, 0.0] for y in ys for x in xs], dtype=float)

    def nid(i: int, j: int) -> int:
        return j * (nx + 1) + i

    conn = [
        [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
        for j in range(ny)
        for i in range(nx)
    ]
    clamped = [nid(0, j) * 6 + d for j in range(ny + 1) for d in range(6)]
    tip_nodes = [nid(nx, j) for j in range(ny + 1)]
    return coords, conn, clamped, tip_nodes


def _solve(asm: PyMeshAssembler, f: np.ndarray, clamped: list[int]) -> np.ndarray:
    """Sparse direct solve with Dirichlet BC."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    free_mask = np.ones(n, dtype=bool)
    free_mask[clamped] = False
    free = np.where(free_mask)[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u_free = spsolve(K[np.ix_(free, free)], f[free])
    u = np.zeros(n)
    u[free] = u_free
    return u


def _run_coupon(lam, nx=NX0, ny=NY0, load=P_TIP, elem_type=44):
    """Solve the coupon. Returns ``(coords, u, tip_nodes)``."""
    prop = _lam_prop(lam)
    coords, conn, clamped, tip_nodes = _strip_mesh(nx, ny)
    asm = PyMeshAssembler(coords, conn, [elem_type] * len(conn), [prop] * len(conn))
    f = np.zeros(asm.dofs_count)
    for nd in tip_nodes:  # mesh-independent resultant
        f[nd * 6 + 2] = load / len(tip_nodes)
    return coords, _solve(asm, f, clamped), tip_nodes


# --- the three twist metrics ------------------------------------------------
# The strip's twist is the rotation of the end cross-section about the X axis, i.e.
# ``dw/dy`` at the free edge. The three metrics below read it from three different
# places of the same displacement field.


def twist_lsq(coords, u, tip_nodes) -> float:
    """Least squares slope of the tip edge in the (y, w) plane [rad]."""
    y = np.array([coords[nd][1] for nd in tip_nodes])
    w = np.array([u[nd * 6 + 2] for nd in tip_nodes])
    return float(np.polyfit(y, w, 1)[0])


def twist_edge(coords, u, tip_nodes) -> float:
    """Exact slope through the two corner nodes of the free edge [rad]."""
    y = np.array([coords[nd][1] for nd in tip_nodes])
    w = np.array([u[nd * 6 + 2] for nd in tip_nodes])
    order = np.argsort(y)
    return float((w[order[-1]] - w[order[0]]) / (y[order[-1]] - y[order[0]]))


def twist_rotation(u, tip_nodes) -> float:
    """Mean nodal rotation about the length axis (dof 3) over the free edge [rad]."""
    return float(np.mean([u[nd * 6 + 3] for nd in tip_nodes]))


def twist_all(coords, u, tip_nodes) -> dict[str, float]:
    return {
        "lsq": twist_lsq(coords, u, tip_nodes),
        "edge": twist_edge(coords, u, tip_nodes),
        "rot_x": twist_rotation(u, tip_nodes),
    }


def twist_profile(coords, u, nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Twist angle ``theta(x) = dw/dy`` at every span station [rad].

    The strip is bent along X and twisted about X, so ``dw/dy`` evaluated at each
    spanwise station ``x_i`` is the rotation of that section. ``theta'(x)`` is the
    interior twist rate: for a laminate in the free-edge plate state it equals
    ``kappa_xy / 2`` (the engineering curvature ``kappa_xy`` is twice the geometric
    twist rate).
    """
    xs, thetas = [], []
    for i in range(nx + 1):
        nds = [j * (nx + 1) + i for j in range(ny + 1)]
        y = np.array([coords[nd][1] for nd in nds])
        w = np.array([u[nd * 6 + 2] for nd in nds])
        xs.append(coords[nds[0]][0])
        thetas.append(float(np.polyfit(y, w, 1)[0]))
    return np.array(xs), np.array(thetas)


def _run_coupon_moment(lam, moment: float, nx: int = NX0, ny: int = NY0, elem_type=44):
    """Solve the coupon under a **constant** tip moment about Y (no tip force).

    A constant end moment is the one load case whose plate solution is exactly the
    free-edge state ``M = (Mx, 0, 0)`` with constant curvatures, so the interior
    curvatures can be compared with ``D^-1 (Mx, 0, 0)`` *directly*, without any
    integration of a varying moment and without a closed-form formula in between.
    """
    prop = _lam_prop(lam)
    coords, conn, clamped, tip_nodes = _strip_mesh(nx, ny)
    asm = PyMeshAssembler(coords, conn, [elem_type] * len(conn), [prop] * len(conn))
    f = np.zeros(asm.dofs_count)
    for nd in tip_nodes:  # total moment about Y, distributed over the free edge
        f[nd * 6 + 4] = moment / len(tip_nodes)
    return coords, _solve(asm, f, clamped), tip_nodes


# ─────────────────────────────────────────────────────────────────────────────
# The issue's CLT reference, reproduced exactly (characterisation, not truth)
# ─────────────────────────────────────────────────────────────────────────────


def clt_twist_issue(lam: Laminate, p: float = P_TIP, length: float = L, index=(2, 0)):
    """``rad2deg( (D^-1)[index] * P L^2 / 2 )`` — the issue's own formula [deg].

    ``index`` distinguishes the two columns the issue itself reports as ambiguous:
    ``(2, 0)`` is what ``_analytic_tip_twist_deg`` used, ``(2, 1)`` is what the
    module docstring stated. Both are reproduced so the ambiguity is measurable.
    """
    D_inv = np.linalg.inv(lam.D)
    return float(np.rad2deg(D_inv[index] * (p * length**2 / 2.0)))


def clt_twist_per_width(lam: Laminate, p: float = P_TIP, length: float = L,
                        width: float = B, index=(2, 0)):
    """The same formula with the moment per unit width, ``P L^2 / (2 b)`` [deg].

    Dimensional check: ``kappa = D^-1 M`` needs a moment *per unit width*, so the
    reference above is a factor ``1/b = 10`` off unless the load is already a force
    per unit width. Reproduced to show which convention the numbers came from.
    """
    D_inv = np.linalg.inv(lam.D)
    return float(np.rad2deg(D_inv[index] * (p * length**2 / (2.0 * width))))


# ─────────────────────────────────────────────────────────────────────────────
# WU1 — the coupon is deterministic, and its laminate is a single object
# ─────────────────────────────────────────────────────────────────────────────


def test_abd_is_the_single_source_of_truth():
    """The ABD the element receives is the ABD the reference uses.

    This is the guard against the original report's structural risk: a reference
    computed from one laminate and a displacement computed from another.
    """
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    prop = _lam_prop(lam)
    assert prop["cm"] == pytest.approx(lam.A.ravel().tolist())
    assert prop["cb"] == pytest.approx(lam.D.ravel().tolist())
    assert prop["b_coupling"] == pytest.approx(lam.B.ravel().tolist())
    assert prop["thickness"] == pytest.approx(lam.total_thickness)
    # [45,0,0,45]s is symmetric, so B must vanish exactly
    assert np.abs(lam.B).max() < 1e-6 * np.abs(lam.D).max()
    # and it is unbalanced in the bending sense: D16 != 0
    assert abs(lam.D[0, 2]) > 1.0


def test_d16_zero_control_gives_zero_twist():
    """``[0,0,0,0]`` (D16 = D26 = 0): the coupling path must produce no twist."""
    lam = _laminate(LAYUPS["[0,0,0,0]"])
    assert lam.D[0, 2] == pytest.approx(0.0, abs=1e-9)
    assert lam.D[1, 2] == pytest.approx(0.0, abs=1e-9)

    coords, u, tips = _run_coupon(lam)
    metrics = twist_all(coords, u, tips)
    for name, value in metrics.items():
        assert abs(value) < 1e-9, f"control layup twisted via metric {name}: {value}"


def test_unbalanced_layup_matches_the_free_edge_reference():
    """The tip-force coupon must match the free-edge reference, 1:1.

    The reference carries the two unit factors the issue's formula omits: the
    moment per unit width (``1/b``) and the geometric twist rate (``kappa_xy / 2``).
    With both of them the element and CLT agree, which is what turns the original
    3.4-4.8x into a units/convention error *of the reference* rather than an
    element defect.
    """
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    coords, u, tips = _run_coupon(lam, nx=8, ny=48)
    measured = twist_lsq(coords, u, tips)
    expected = tip_force_free_edge_reference(lam)
    ratio = measured / expected

    naive = np.deg2rad(clt_twist_issue(lam))
    print("\ntip-force coupon [45,0,0,45]s (mesh 8x48):")
    print(f"  measured twist            = {measured:+.6e} rad")
    print(f"  issue formula (no 1/b,2x) = {naive:+.6e} rad   ratio = {measured / naive:+.3f}")
    print(f"  free-edge reference       = {expected:+.6e} rad   ratio = {ratio:+.4f}")
    assert abs(ratio - 1.0) < 0.05, (
        f"measured {measured:.6e} rad vs free-edge reference {expected:.6e} rad "
        f"(ratio {ratio:.4f})"
    )


# ─────────────────────────────────────────────────────────────────────────────
# The constitutive coupling, measured inside a real strip (constant tip moment)
# ─────────────────────────────────────────────────────────────────────────────

#: Consistent tip moment [N·m] for the constant-moment probe. Chosen so the
#: resulting curvature is of the same order as the tip-force case.
M_TIP = 1.0


def moment_reference(lam: Laminate, moment: float = M_TIP, length: float = L,
                     width: float = B) -> tuple[float, float]:
    """The **free-edge** twist reference, with both unit factors explicit.

    For a laminate strip in the free-edge plate state (``M = (Mx, 0, 0)``,
    constant curvatures) the exact CLT result is ``kappa = D^-1 (Mx, 0, 0)``. Two
    unit factors separate that from the formula in issue #9:

    * the moment resultant is **per unit width**, ``Mx = m / b`` -- not the total
      moment ``m``;
    * the *geometric* twist rate of the section is ``theta' = kappa_xy / 2``, because
      ``kappa_xy`` is the **engineering** twist curvature (``kappa_xy = 2 d2w/dxdy``).

    Returns ``(theta_prime, theta_tip)`` for a uniform ``Mx`` over the length.
    """
    D_inv = np.linalg.inv(lam.D)
    kappa_xy = D_inv[2, 0] * (moment / width)
    theta_prime = kappa_xy / 2.0
    return float(theta_prime), float(theta_prime * length)


def tip_force_free_edge_reference(lam: Laminate, p: float = P_TIP, length: float = L,
                                  width: float = B) -> float:
    """Free-edge strip twist under a **tip force**, in radians.

    The moment resultant per unit width is ``Mx(x) = p (L - x) / b`` and the
    geometric twist rate is ``theta'(x) = (D^-1)[2,0] * Mx(x) / 2``, so

        theta(L) = (D^-1)[2,0] * p * L^2 / (4 b).

    This is the free-edge (cylindrical-bending) plate state integrated along the
    strip. It is exact for a *constant* moment; for a tip force it neglects the
    y-redistribution that a varying moment forces through compatibility, so it is a
    first-order reference whose residual is what WU4 (a converged Ritz solution of
    the same BVP) has to settle.
    """
    D_inv = np.linalg.inv(lam.D)
    return float(D_inv[2, 0] * p * length**2 / (4.0 * width))


def test_constant_moment_interior_state_is_the_free_edge_state():
    """The interior of a constant-moment strip must satisfy ``M = D kappa``.

    This isolates the element's coupling from the boundary-value problem. For a
    constant end moment the exact plate solution has constant curvatures equal to
    ``D^-1 (Mx, 0, 0)`` with ``Mx = m / b``, and it satisfies the free-edge
    conditions exactly. So the *measured* interior twist rate must be
    ``kappa_xy / 2 = (D^-1)[2,0] * m / (2 b)`` if the element implements the CLT
    coupling -- whatever the tip-force reference formula turns out to be.

    Sign convention (the repo's, tests/test_mitc3_benchmarks.py:253): a positive
    nodal moment about +y rotates the tip section from +x toward -z, so ``w_tip``
    is negative and the *effective* ``Mx`` is ``-m/b``.
    """
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    coords, u, tips = _run_coupon_moment(lam, M_TIP, nx=8, ny=24)
    xs, thetas = twist_profile(coords, u, 8, 24)

    # interior window, away from the clamped edge and from the loaded edge
    mid = (xs > 0.25 * L) & (xs < 0.75 * L)
    rate = float(np.polyfit(xs[mid], thetas[mid], 1)[0])  # theta'(x) [rad/m]

    # lock the repo's documented tip-moment sign convention on the same solution:
    # a positive moment about +y must push the tip down in -z.
    w_tip = float(np.mean([u[nd * 6 + 2] for nd in tips]))
    assert w_tip < 0.0, f"positive +y tip moment must give w_tip < 0, got {w_tip}"

    mxx = -M_TIP / B  # effective moment resultant per unit width [N]
    D_inv = np.linalg.inv(lam.D)
    kappa = D_inv @ np.array([mxx, 0.0, 0.0])
    predicted_rate = kappa[2] / 2.0  # geometric twist rate 1/m

    ratio = rate / predicted_rate
    print(f"\nconstant tip moment {M_TIP} N.m -> Mx = {mxx:.1f} N per width "
          f"(mesh 8x24, interior 25-75%)")
    print(f"  w_tip                = {w_tip:+.6e} m")
    print(f"  CLT D^-1 (Mx,0,0)    = {kappa}   (curvatures 1/m)")
    print(f"  measured theta'      = {rate:+.6e} rad/m")
    print(f"  CLT kappa_xy / 2     = {predicted_rate:+.6e} rad/m   ratio = {ratio:+.4f}")
    print(f"  CLT kappa_xy         = {kappa[2]:+.6e} rad/m   ratio = {rate / kappa[2]:+.4f}")

    # The element must reproduce the free-edge state: this is the guard that makes
    # the 3.4-4.8x claim falsifiable.
    assert abs(ratio - 1.0) < 0.05, (
        f"interior twist rate {rate:.4e} rad/m is {ratio:.3f}x the free-edge CLT "
        f"value {predicted_rate:.4e} rad/m"
    )


def test_constant_moment_tip_twist_matches_the_free_edge_state():
    """Tip twist under a constant moment equals the free-edge state, 1:1."""
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    coords, u, tips = _run_coupon_moment(lam, M_TIP, nx=8, ny=24)
    measured = twist_lsq(coords, u, tips)
    _, expected = moment_reference(lam, moment=-M_TIP)
    ratio = measured / expected
    print(f"\nconstant tip moment {M_TIP} N.m: measured tip twist = {measured:+.6e} rad, "
          f"free-edge CLT = {expected:+.6e} rad, ratio = {ratio:+.4f}")
    assert abs(ratio - 1.0) < 0.05, f"tip twist ratio {ratio:.3f} is not 1"


# ─────────────────────────────────────────────────────────────────────────────
# WU3 — metric sensitivity: the ratio must not be an artefact of the metric
# ─────────────────────────────────────────────────────────────────────────────


def test_metric_sensitivity():
    """The three twist metrics must agree on a refined mesh.

    On the issue's coarse mesh (2 elements along a 1 m cantilever) the free edge
    carries ``NY + 1 = 13`` nodes but only two element lengths of curvature, so the
    metrics are expected to disagree. On a refined mesh they must converge to each
    other; if they do not, any single-metric ratio is a metric artefact (H2).
    """
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    for nx, ny in ((NX0, NY0), (8, 48)):
        coords, u, tips = _run_coupon(lam, nx=nx, ny=ny)
        m = twist_all(coords, u, tips)
        spread = (max(m.values()) - min(m.values())) / max(abs(m["lsq"]), 1e-30)
        print(f"  mesh {nx}x{ny}: lsq={m['lsq']:+.8e} edge={m['edge']:+.8e} "
              f"rot_x={m['rot_x']:+.8e}  spread={spread:.1%}")
    # On the refined mesh the metrics must agree to a few percent.
    coords, u, tips = _run_coupon(lam, nx=8, ny=48)
    m = twist_all(coords, u, tips)
    assert abs(m["edge"] - m["lsq"]) < 0.05 * abs(m["lsq"]), m
    assert abs(m["rot_x"] - m["lsq"]) < 0.05 * abs(m["lsq"]), m


def test_ratio_deviation_shrinks_with_width_refinement():
    """The residual left of the free-edge reference must behave like a mesh effect.

    The free-edge reference is exact for a constant moment. Under a tip force the
    measured/expected ratio ranges over the sweep (roughly 0.87-1.32), and the
    question is whether that residual is discretisation or physics. Refining *across
    the width* is the discriminating knob: the y-structure of the strip is what a
    free-edge state suppresses, so if the deviation is a coarse-width artefact it
    must shrink with NY, and if it is real plate physics it must converge to a
    non-unity limit.
    """
    for name in ("[45,0,0,45]s", "[60,0,0,60]s"):
        lam = _laminate(LAYUPS[name])
        expected = tip_force_free_edge_reference(lam)
        print(f"\n{name}: measured / free-edge reference vs NY (NX = 8)")
        ratios = []
        for ny in (12, 24, 48, 96):
            coords, u, tips = _run_coupon(lam, nx=8, ny=ny)
            ratio = twist_lsq(coords, u, tips) / expected
            ratios.append(ratio)
            print(f"  NY={ny:3d}: ratio = {ratio:+.4f}")
        # the sequence must settle: the last two refinements closer than the first two
        assert abs(ratios[3] - ratios[2]) < abs(ratios[1] - ratios[0]), (name, ratios)


# ─────────────────────────────────────────────────────────────────────────────
# WU2 — mesh convergence: the ratio must not be a coarse-mesh artefact
# ─────────────────────────────────────────────────────────────────────────────


def test_mesh_convergence_of_the_twist_ratio():
    """Refine NX x NY and report shell / reference (H3).

    The issue reports one mesh (2 x 12). If the ratio moves with refinement, the
    reported 3.4-4.8x was a discretisation artefact of a two-element cantilever.
    """
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    ref20 = clt_twist_issue(lam, index=(2, 0))
    ref21 = clt_twist_issue(lam, index=(2, 1))
    ref_pw = clt_twist_per_width(lam)

    rows = []
    for nx, ny in ((2, 12), (4, 24), (8, 48), (16, 96)):
        coords, u, tips = _run_coupon(lam, nx=nx, ny=ny)
        twist = np.rad2deg(twist_lsq(coords, u, tips))
        rows.append((nx, ny, twist, twist / ref20, twist / ref21, twist / ref_pw))
        print(f"  {nx:2d}x{ny:3d}: shell={twist:+.6f} deg  "
              f"sh/issue[2,0]={twist / ref20:+.3f}  sh/issue[2,1]={twist / ref21:+.3f}  "
              f"sh/issue[2,0],1/b={twist / ref_pw:+.3f}  "
              f"sh/free-edge={twist / np.rad2deg(tip_force_free_edge_reference(lam)):+.4f}")

    # The measured twist must be mesh-convergent (in the Cauchy sense): the last two
    # refinements must agree far better than the first pair. Which side the limit
    # lands on is the WU4 verdict, not this test's business.
    d_coarse = abs(rows[1][2] - rows[0][2])
    d_fine = abs(rows[3][2] - rows[2][2])
    assert d_fine < 0.5 * d_coarse, (
        "twist is not converging with mesh refinement: "
        f"|d(4x24 - 2x12)|={d_coarse:.3e} vs |d(16x96 - 8x48)|={d_fine:.3e}"
    )


# ─────────────────────────────────────────────────────────────────────────────
# WU4 route (a) — a converged Rayleigh-Ritz solution of the same plate BVP
# ─────────────────────────────────────────────────────────────────────────────
# Independent of `mitc4.rs`, of the CCX writer and of any closed-form formula: the
# CLT bending energy of the same strip is minimised over a truncated global basis.

def kirchhoff_clt_energy_matrix(lam: Laminate, p_ord: int = 6, q_ord: int = 16,
                                length: float = L, width: float = B,
                                nq_x: int = 40, nq_y: int = 64):
    """Ritz stiffness and load for a clamped-free CLT strip.

    Basis: ``w(x,y) = sum_{p>=2, q} a_pq (x/L)^p P_q(2y/b)`` with ``P_q`` the
    Legendre polynomials. ``p >= 2`` enforces the clamped edge ``w = w,x = 0``
    exactly, and ``w,y`` is left free there (the classical plate clamped condition).
    The basis spans the free-edge state, so in the interior the Ritz solution must
    approach ``D^-1 (Mx, 0, 0)`` whenever the exact solution is constant-curvature.
    """
    half = 0.5 * width
    gx, wx = np.polynomial.legendre.leggauss(nq_x)
    gy, wy = np.polynomial.legendre.leggauss(nq_y)
    X = 0.5 * length * (gx + 1.0)
    WX = 0.5 * length * wx
    S = gy  # y / half
    WY = half * wy
    XX = X[:, None]
    SS = S[None, :]
    WW = np.outer(WX, WY)

    def legendre_basis(q: int, s: np.ndarray):
        coef = np.zeros(q + 1)
        coef[q] = 1.0
        d1 = np.polynomial.legendre.legder(coef, 1) if q >= 1 else np.array([0.0])
        d2 = np.polynomial.legendre.legder(coef, 2) if q >= 2 else np.array([0.0])
        return (
            np.polynomial.legendre.legval(s, coef),
            np.polynomial.legendre.legval(s, d1),
            np.polynomial.legendre.legval(s, d2),
        )

    psi = [legendre_basis(q, SS) for q in range(q_ord + 1)]

    modes = []  # (p, q, kappa array (3, ny, nx))
    for p in range(2, p_ord + 1):
        phi = (XX / length) ** p
        dphi = p * (XX / length) ** (p - 1) / length
        d2phi = p * (p - 1) * (XX / length) ** (p - 2) / length**2
        for q in range(q_ord + 1):
            ps, d1, d2 = psi[q]
            kx = d2phi * ps / 1.0
            ky = phi * (d2 / half**2)
            kxy = 2.0 * dphi * (d1 / half)
            modes.append((p, q, np.stack([kx, ky, kxy], axis=0)))

    D = lam.D
    n = len(modes)
    # kappa of every mode, flattened over the quadrature grid: (n, 3, N) -> (n, 3N)
    A = np.stack([m[2].reshape(3, -1) for m in modes])          # (n, 3, N)
    B = np.einsum("ab,ibp->iap", D, A)                          # (n, 3, N)
    Wf = WW.reshape(-1)                                         # (N,)
    K = (A * Wf).reshape(n, -1) @ B.reshape(n, -1).T
    return K, modes, (X, S, WX, WY)


def ritz_tip_twist(lam: Laminate, load: str = "force", magnitude: float = P_TIP,
                   p_ord: int = 6, q_ord: int = 16, n_fit: int = 201) -> float:
    """Tip twist of the Ritz solution, read with the *element's* metric (LSQ dw/dy).

    ``load='force'``  : uniform line load ``P/b`` over the free edge.
    ``load='moment'`` : uniform line moment ``m/b`` about Y at the free edge.
    """
    if load not in ("force", "moment"):
        raise ValueError(load)
    K, modes, (X, S, WX, WY) = kirchhoff_clt_energy_matrix(
        lam, p_ord=p_ord, q_ord=q_ord
    )
    half = 0.5 * B
    f = np.zeros(len(modes))
    for i, (p, q, _) in enumerate(modes):
        # integral over y of P_q(2y/b) dy, by the same quadrature as the energy
        coef_q = np.zeros(q + 1)
        coef_q[q] = 1.0
        iy = float(np.sum(WY * np.polynomial.legendre.legval(S, coef_q)))
        if load == "force":
            f[i] = (magnitude / B) * 1.0 * iy        # phi_p(L) = 1
        else:
            f[i] = (magnitude / B) * (p / L) * iy    # phi_p'(L) = p/L
    c = np.linalg.solve(K, f)

    # w(L, y) and the LSQ slope dw/dy over the free edge
    y = np.linspace(-half, half, n_fit)
    s = y / half
    w = np.zeros_like(y)
    for (_, q, _), ci in zip(modes, c, strict=True):
        coef_q = np.zeros(q + 1)
        coef_q[q] = 1.0
        w += ci * 1.0 * np.polynomial.legendre.legval(s, coef_q)
    return float(np.polyfit(y, w, 1)[0])


def test_ritz_reference_converges_and_matches_the_element():
    """The independent Ritz reference vs the element, per layup.

    This is WU4's verdict surface. The Ritz solution is an upper bound on the
    compliance of the same CLT BVP and shares no code with the element, so if the
    element converges to it, the reported 3.4-4.8x is fully accounted for by the
    reference's missing unit factors.
    """
    print("\nRitz convergence of the tip twist [deg], P order 6, Q order swept:")
    lam0 = _laminate(LAYUPS["[45,0,0,45]s"])
    for q_ord in (6, 8, 12, 16, 20, 24):
        print(f"  Q={q_ord:2d}: {np.rad2deg(ritz_tip_twist(lam0, q_ord=q_ord)):+.6f}")

    print("\nelement (8x48) vs Ritz (Q=20) vs free-edge formula, tip force:")
    print(f"  {'layup':>18} {'shell[deg]':>11} {'ritz[deg]':>11} {'freeedge':>10} "
          f"{'shell/ritz':>10} {'shell/free':>10}")
    for name, half_stack in LAYUPS.items():
        lam = _laminate(half_stack)
        coords, u, tips = _run_coupon(lam, nx=8, ny=48)
        shell = np.rad2deg(twist_lsq(coords, u, tips))
        ritz = np.rad2deg(ritz_tip_twist(lam, q_ord=20))
        free = np.rad2deg(tip_force_free_edge_reference(lam))
        r_shell_ritz = shell / ritz if abs(ritz) > 1e-12 else float("nan")
        r_shell_free = shell / free if abs(free) > 1e-12 else float("nan")
        print(f"  {name:>18} {shell:11.6f} {ritz:11.6f} {free:10.6f} "
              f"{r_shell_ritz:10.4f} {r_shell_free:10.4f}")

    # The machinery control: with a constant tip moment the exact interior state is
    # the free-edge state, so the *tip* twist must land within a few percent of it.
    # (The Ritz applies ``Mx = +m/b``; ``moment_reference`` follows the element's
    # documented sign convention, so both are taken with the same sign of ``m``.)
    _, expected_moment = moment_reference(lam0, moment=M_TIP)
    ritz_moment = ritz_tip_twist(lam0, load="moment", magnitude=M_TIP, q_ord=20)
    print(f"\nconstant moment control: ritz = {ritz_moment:+.6e} rad, "
          f"free-edge = {expected_moment:+.6e} rad, ratio = {ritz_moment / expected_moment:+.4f}")
    assert abs(ritz_moment / expected_moment - 1.0) < 0.05, ritz_moment

    # The element must land on the Ritz reference for the issue's coupon.
    lam = _laminate(LAYUPS["[45,0,0,45]s"])
    coords, u, tips = _run_coupon(lam, nx=8, ny=48)
    shell = twist_lsq(coords, u, tips)
    ritz = ritz_tip_twist(lam, q_ord=20)
    assert abs(shell / ritz - 1.0) < 0.05, (
        f"element {shell:.6e} rad vs Ritz {ritz:.6e} rad (ratio {shell / ritz:.4f})"
    )


# ─────────────────────────────────────────────────────────────────────────────
# WU4 route (b) — an independent judge: a hand-authored CCX S8R deck
# ─────────────────────────────────────────────────────────────────────────────
# `write_ccx_mesh` is deliberately **not** used. An earlier audit localised a
# defect in its ``*SHELL SECTION, COMPOSITE`` block (per-ply orientations realised
# so that ``[0/90]2s`` came out 34.5% off an independent CLT while a smeared
# ``MATERIAL=`` section matched to 0.6%), so a deck built by that writer cannot
# judge the element. The deck below is authored from the CalculiX manual's own
# syntax, and validated by a decoupled control and a sign-flip control.

_CCX_E1, _CCX_E2 = 120e9, 10e9
_CCX_NU12, _CCX_G12, _CCX_G23 = 0.3, 5e9, 3e9
#: Wall-clock bound for one CalculiX run [s]; a hung solver must fail fast (R4-ccx-timeout).
_CCX_TIMEOUT_S = 300


def write_ccx_coupon(path, angles, ply_t, nx, ny, load=P_TIP):
    """Write the S8R composite coupon deck; return the free-edge node ids."""
    n2x, n2y = 2 * nx, 2 * ny

    def nid(i, j):
        return j * (n2x + 1) + i + 1

    coords = [
        (L * i / n2x, -0.5 * B + B * j / n2y)
        for j in range(n2y + 1)
        for i in range(n2x + 1)
    ]
    tip = [nid(n2x, j) for j in range(n2y + 1)]
    clamp = [nid(0, j) for j in range(n2y + 1)]

    def wset(handle, name, ids):
        handle.write(f"*NSET, NSET={name}\n")
        for k in range(0, len(ids), 16):  # CalculiX: at most 16 entries per line
            handle.write(", ".join(str(c) for c in ids[k:k + 16]) + "\n")

    with open(path, "w", encoding="ascii") as fh:
        fh.write("*HEADING\nbend-twist coupon: hand-authored S8R composite deck\n")
        fh.write("*NODE, NSET=NALL\n")
        for k, (x, y) in enumerate(coords, start=1):
            fh.write(f"{k}, {x:.10e}, {y:.10e}, 0.0\n")
        fh.write("*ELEMENT, TYPE=S8R, ELSET=EALL\n")
        k = 0
        for j in range(ny):
            for i in range(nx):
                k += 1
                i0, j0 = 2 * i, 2 * j
                conn = [nid(i0, j0), nid(i0 + 2, j0), nid(i0 + 2, j0 + 2), nid(i0, j0 + 2),
                        nid(i0 + 1, j0), nid(i0 + 2, j0 + 1), nid(i0 + 1, j0 + 2),
                        nid(i0, j0 + 1)]
                fh.write(f"{k}, " + ", ".join(str(c) for c in conn) + "\n")
        wset(fh, "CLAMP", clamp)
        wset(fh, "TIP", tip)
        fh.write("*MATERIAL, NAME=PLY\n*ELASTIC, TYPE=ENGINEERING CONSTANTS\n")
        fh.write(f"{_CCX_E1}, {_CCX_E2}, {_CCX_E2}, {_CCX_NU12}, {_CCX_NU12}, {_CCX_NU12}, "
                 f"{_CCX_G12}, {_CCX_G12}\n{_CCX_G23}, 0.0\n")
        fh.write("*DENSITY\n1500.0\n")
        fh.write("*SHELL SECTION, ELSET=EALL, COMPOSITE\n")
        for ang in angles:
            tag = f"{abs(ang):g}_{'M' if ang < 0 else 'P'}"
            fh.write(f"{ply_t:.10e}, , PLY, ORI_{tag}\n")
        for Ang in sorted({abs(float(a)) for a in angles}):
            for sgn, tag in ((1.0, "P"), (-1.0, "M")):
                fh.write(f"*ORIENTATION, NAME=ORI_{Ang:g}_{tag}, SYSTEM=RECTANGULAR\n")
                fh.write("1., 0., 0., 0., 1., 0.\n")
                fh.write(f"3, {sgn * Ang:.4f}\n")
        fh.write("*STEP\n*STATIC\n")
        fh.write("*BOUNDARY\nCLAMP, 1, 6, 0.0\n")
        fh.write("*CLOAD\n")
        for nd in tip:
            fh.write(f"{nd}, 3, {load / len(tip):.10e}\n")
        fh.write("*NODE FILE, OUTPUT=2D\nU\n")
        fh.write("*EL FILE, OUTPUT=2D\nS\n")
        fh.write("*END STEP\n")
    return coords, tip


def ccx_tip_twist(half, ply_t=PLY_T, nx=4, ny=24, workdir=None):
    """Run the hand-authored deck and return the free-edge LSQ twist [rad]."""
    import subprocess
    from pathlib import Path

    from conftest import ccx_bin_or_skip
    from _ccx_io import parse_frd_disp

    ccx = ccx_bin_or_skip()
    work = Path(workdir) if workdir else Path("/tmp") / "ccx_bend_twist"
    work.mkdir(parents=True, exist_ok=True)
    angles = [float(a) for a in half] + [float(a) for a in reversed(half)]
    coords, tip = write_ccx_coupon(work / "c.inp", angles, ply_t, nx, ny)
    try:
        result = subprocess.run(
            [ccx, "c"], cwd=work, capture_output=True, text=True, timeout=_CCX_TIMEOUT_S
        )
    except subprocess.TimeoutExpired:
        pytest.fail(f"CalculiX did not finish within {_CCX_TIMEOUT_S}s on the coupon deck")
    if result.returncode != 0:
        pytest.fail(f"CCX failed (rc={result.returncode}):\n{result.stdout[-2000:]}")
    disp = parse_frd_disp(work / "c.frd", tip)
    if not disp:
        pytest.fail("no displacements parsed from the CCX FRD")
    u = np.zeros(6 * len(coords))
    for node, d in disp.items():
        u[(node - 1) * 6:(node - 1) * 6 + 3] = d
    return twist_lsq(np.array([[x, y, 0.0] for x, y in coords]), u, [c - 1 for c in tip])


def test_independent_ccx_deck_confirms_the_coupon(tmp_path):
    """The element against an independently authored CalculiX S8R model.

    Controls first: the decoupled layup must give zero twist in CCX too (this
    validates the deck, the section, the load path and the parser), and the
    ``+/-45`` pair must give opposite signs (this validates that the per-ply
    orientations are actually realised). Then the issue's coupon must agree with
    the element at matched refinement.
    """
    control = ccx_tip_twist([0.0, 0.0, 0.0, 0.0], nx=2, ny=12, workdir=tmp_path / "ctrl")
    print(f"\nCCX deck control [0,0,0,0]: twist = {np.rad2deg(control):+.8f} deg")
    assert abs(control) < 1e-9, "decoupled control twisted: the deck or section is wrong"

    plus = ccx_tip_twist([45.0, 0.0, 0.0, 45.0], nx=4, ny=24, workdir=tmp_path / "plus")
    minus = ccx_tip_twist([-45.0, 0.0, 0.0, -45.0], nx=4, ny=24, workdir=tmp_path / "minus")
    assert plus * minus < 0.0, "ply orientation sign is not realised by the deck"
    assert abs(plus + minus) < 0.02 * abs(plus), (plus, minus)

    print("element vs hand-authored CCX S8R, matched refinement (converged):")
    for name in ("[45,0,0,45]s", "[60,0,0,60]s", "[45,-45,-45,45]s"):
        lam = _laminate(LAYUPS[name])
        ccx = ccx_tip_twist(LAYUPS[name], nx=8, ny=48, workdir=tmp_path / name[1:-1])
        coords, u, tips = _run_coupon(lam, nx=16, ny=96)
        elem = twist_lsq(coords, u, tips)
        ratio = elem / ccx
        print(f"  {name:>18}: ccx={np.rad2deg(ccx):+.6f} deg  element={np.rad2deg(elem):+.6f} deg  "
              f"element/ccx = {ratio:.4f}")
        assert abs(ratio - 1.0) < 0.03, (name, ratio)


# ─────────────────────────────────────────────────────────────────────────────
# Report: the issue's sweep table, regenerated (run with -s)
# ─────────────────────────────────────────────────────────────────────────────


def test_report_bend_twist_table():
    """Print the issue's sweep with this tree's own numbers (``pytest -s``).

    This is a measurement, not an assertion: it exists so the tables in
    ``odd/tasks/composite-bend-twist-verdict.md`` can be regenerated.
    """
    coords0, _, _, _ = _strip_mesh(NX0, NY0)
    print(f"\ncoupon: L={L} b={B} h={8 * PLY_T * 1e3:.0f} mm  P={P_TIP} N  mesh {NX0}x{NY0}")
    header = (f"{'layup':>18} {'D16':>9} {'D26':>9} {'shell[deg]':>11} "
              f"{'issue20':>9} {'free-edge':>10} {'sh/issue':>9} {'sh/freeedge':>12}")
    print(header)
    for name, half in LAYUPS.items():
        lam = _laminate(half)
        _, u, tips = _run_coupon(lam)
        shell = np.rad2deg(twist_lsq(coords0, u, tips))
        c20 = clt_twist_issue(lam, index=(2, 0))
        c21 = clt_twist_issue(lam, index=(2, 1))
        cpw = clt_twist_per_width(lam)
        cfe = np.rad2deg(tip_force_free_edge_reference(lam))
        ratio = shell / c20 if abs(c20) > 1e-12 else float("nan")
        rfe = shell / cfe if abs(cfe) > 1e-12 else float("nan")
        print(f"{name:>18} {lam.D[0, 2]:9.2f} {lam.D[1, 2]:9.2f} {shell:11.6f} "
              f"{c20:9.5f} {cfe:10.5f} {ratio:9.3f} {rfe:12.4f}   (CLT21={c21:+.5f}, CLTpw={cpw:+.5f})")

    print("\nsame sweep, T300/5208 plies:")
    for name, half in LAYUPS.items():
        lam = _laminate(half, _T300)
        _, u, tips = _run_coupon(lam)
        shell = np.rad2deg(twist_lsq(coords0, u, tips))
        c20 = clt_twist_issue(lam, index=(2, 0))
        cfe = np.rad2deg(tip_force_free_edge_reference(lam))
        rfe = shell / cfe if abs(cfe) > 1e-12 else float("nan")
        print(f"{name:>18} {lam.D[0, 2]:9.2f} {shell:11.6f} {c20:9.5f} {cfe:10.5f} {rfe:12.4f}")
