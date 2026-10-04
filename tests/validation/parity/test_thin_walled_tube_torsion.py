"""Saint-Venant torsion of a closed rectangular thin-walled tube — the load application and the metric.

This module is **WU-B2 step 2** of ``odd/tasks/composite-bend-twist-verdict.md`` (§18.5/§18.6/§19).
The blade's rated twist magnitude was withdrawn because the same load set moved the twist by a factor
of 7.8 depending only on how the forces and moments were distributed over the shell. This file settles
two questions on a case whose answer is exact - the Bredt-Batho torsion of a closed thin-walled tube,

    theta' = T / GJ,        GJ = 4 A^2 A66 / perimeter,

with ``A66`` the wall's in-plane shear stiffness per unit width (``G t`` isotropic,
``t_total * Qbar66(45 deg)`` for a laminated wall), ``A`` the **mid-line** enclosed area and the
integral over the **mid-line** perimeter:

1. **The load application.** A tip shear flow ``q = T/(2A)`` must reproduce ``T/GJ``. It does
   (§ ``_self_equilibrated``). A statically equivalent torque applied as a traction on one wall does
   not, and the mechanism is a section-distortion mode (the last test).
2. **Which twist metric is the section rotation.** Two metrics are computed side by side in every
   case: ``theta_fit`` (in-plane rigid-ring fit) and ``theta_z`` (mean DOF-5 rotation about z). In
   the exact solution both must equal the section rotation; the self-equilibrated case shows they do,
   and the clamped-root convergence table shows ``theta_z`` is the one that does not depend on the
   end boundary layer.

The boundary-layer-free experiment is the primary validation. The original clamped-root case is kept
as the quantified engineering comparison, and its two metrics are resolved with an L x t convergence
table rather than an assertion. The one-wall traction is recast to assert what the physics pins
(the mirror, the torque ruler, the distortion amplitude) and to report the rest.

Geometry: ``b = 1.0 m`` (x), ``h = 0.6 m`` (y), axis along z, ``L`` from 6 to 24 m, wall
``t = 0.01 m`` (isotropic) or 0.5 mm (laminate), ``n_seg = 4`` segments per side with corner nodes
shared (16-node closed ring), 0.2 m elements. ``T = 1.0e4 N.m`` throughout. The reference is
hand-written and independent of everything the repository computes.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler  # noqa: E402
from aeroelast.core.laminate import create_laminate_from_angles  # noqa: E402
from aeroelast.core.material import OrthotropicMaterial  # noqa: E402
from scipy.sparse import block_diag, bmat, coo_matrix, csr_matrix  # noqa: E402
from scipy.sparse.linalg import splu, spsolve

from tests.support.assertions import assert_relative_error  # noqa: E402  # noqa: E402

# ─────────────────────────────────────────────────────────────────────────────
# The exact case (fixed spec)
# ─────────────────────────────────────────────────────────────────────────────

B = 1.0  # mid-line width along x [m]
H = 0.6  # mid-line height along y [m]
L = 6.0  # reference tube length along z [m]
TORQUE = 1.0e4  # applied torque [N.m]
N_SEG = 4  # segments per side (corners shared -> ring of 4*n_seg nodes)
N_Z = 30  # elements along z for L = 6
TOL = 0.05  # the suite's 5% rule
SELF_CONSISTENCY = 0.01  # two sub-windows must agree within 1%
METRIC_AGREEMENT = 0.02  # theta_fit and theta_z must agree within 2% in the exact solution

ENCLOSED_AREA = B * H  # mid-line enclosed area [m^2]
PERIMETER = 2.0 * (B + H)  # mid-line perimeter [m]
DZ = L / N_Z  # element length along z [m]

E_ISO = 2.1e11  # [Pa]
NU_ISO = 0.3
RHO_ISO = 7800.0  # [kg/m^3]
THICKNESS = 0.01  # isotropic wall thickness [m]
G_ISO = E_ISO / (2.0 * (1.0 + NU_ISO))

E1 = 120e9  # [Pa]
E2 = 10e9  # [Pa]
G12 = 5e9  # [Pa]  (G13 = G23 = G12)
NU12 = 0.3
RHO_PLY = 1500.0  # [kg/m^3]
PLY_T = 0.125e-3  # ply thickness [m]
LAM_PLIES = 4  # [45, -45]s -> 4 plies
LAM_TOTAL_THICKNESS = LAM_PLIES * PLY_T  # [m]

#: Penalty scale for the documented rigid-mode penalty. The assembled K only annihilates the
#: analytic rigid modes to round-off, so the penalty is not exactly solution-neutral: measured
#: against the exact constraint it drifts the twist rate by 2.0e-5 relative at the plan's
#: ``1e-6 * median(diag K)`` and by more for the laminate. ``1e-14`` is below the measured drift
#: floor of the isotropic tube while remaining solvable.
_PENALTY_MEDIAN = 1e-14
#: The plan's prescribed penalty scale, kept for the measured drift report.
_PENALTY_PRESCRIBED = 1e-6


# ─────────────────────────────────────────────────────────────────────────────
# Mesh
# ─────────────────────────────────────────────────────────────────────────────


def _midline_ring(b: float = B, h: float = H, n_seg: int = N_SEG) -> list[tuple[float, float]]:
    """Closed mid-line rectangle, corners shared: 4*n_seg nodes around the perimeter."""
    pts: list[tuple[float, float]] = []
    for j in range(n_seg):  # bottom, y = -h/2, x from -b/2 to +b/2
        pts.append((-b / 2.0 + j * b / n_seg, -h / 2.0))
    for j in range(n_seg):  # right wall, x = +b/2, y from -h/2 to +h/2
        pts.append((b / 2.0, -h / 2.0 + j * h / n_seg))
    for j in range(n_seg):  # top, y = +h/2, x from +b/2 down to -b/2
        pts.append((b / 2.0 - j * b / n_seg, h / 2.0))
    for j in range(n_seg):  # left wall, x = -b/2, y from +h/2 down to -h/2
        pts.append((-b / 2.0, h / 2.0 - j * h / n_seg))
    return pts


def _tube_mesh(length: float = L, n_z: int = N_Z, n_seg: int = N_SEG):
    """Hand-built quad grid: `(coords, conn, rings, n_ring)`."""
    dz = length / n_z
    ring = _midline_ring(n_seg=n_seg)
    n_ring = len(ring)
    coords: list[list[float]] = []
    rings: list[list[int]] = []
    for k in range(n_z + 1):
        ids = []
        for x, y in ring:
            ids.append(len(coords))
            coords.append([x, y, k * dz])
        rings.append(ids)
    coords_arr = np.asarray(coords, dtype=float)
    conn = []
    for k in range(n_z):
        for i in range(n_ring):
            j = (i + 1) % n_ring
            conn.append([rings[k][i], rings[k][j], rings[k + 1][j], rings[k + 1][i]])
    return coords_arr, conn, rings, n_ring


# ─────────────────────────────────────────────────────────────────────────────
# Properties (the mechanics of tests/test_laminate_bend_twist.py, copied not imported)
# ─────────────────────────────────────────────────────────────────────────────


def _iso_prop(thickness: float = THICKNESS) -> dict:
    return {
        "type": "isotropic",
        "e": E_ISO,
        "nu": NU_ISO,
        "rho": RHO_ISO,
        "thickness": thickness,
        "shear_correction": 5.0 / 6.0,
    }


def _laminate_and_prop():
    """The [45, -45]s CFRP laminate and the composite property dict the element consumes (code 44).

    Balanced symmetric is the right wall here: ``A16 = A26 = 0`` and ``B = 0``, so the wall has no
    extension-shear or bend-twist coupling and its only job under a shear flow is the membrane shear
    stiffness ``A66``, for which ``GJ = 4 A^2 A66 / perimeter`` is exact.
    """
    material = OrthotropicMaterial(
        "CFRP", E=(E1, E2, E2), G=(G12, G12, G12), nu=(NU12, NU12, NU12), rho=RHO_PLY
    )
    angles = [45.0, -45.0, -45.0, 45.0]  # half [45, -45] plus its reversed symmetric expansion
    lam = create_laminate_from_angles(material, PLY_T, angles)
    # A property of the fixture, checked where the fixture is built: four plies of PLY_T must give
    # LAM_TOTAL_THICKNESS. It used to sit in the test body, where the extractor counted it as a
    # comparison of the test against a reference, which is what the contract says it is not.
    assert abs(lam.total_thickness - LAM_TOTAL_THICKNESS) < 1e-15, (
        f"laminate thickness {lam.total_thickness} is not the {LAM_TOTAL_THICKNESS} of four plies"
    )
    h = lam.total_thickness
    mpa = sum(p.material.rho * p.thickness for p in lam.plies)
    ri = sum(p.material.rho * (p.z_top**3 - p.z_bottom**3) / 3.0 for p in lam.plies)
    prop = {
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
    return lam, prop


# ─────────────────────────────────────────────────────────────────────────────
# Assembly, rigid modes, solvers
# ─────────────────────────────────────────────────────────────────────────────


def _assemble(coords, conn, elem_type, prop):
    """`(assembler, K)` for a uniform element type and property."""
    asm = PyMeshAssembler(coords, conn, [elem_type] * len(conn), [prop] * len(conn))
    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    return asm, K


def _rigid_modes(coords: np.ndarray) -> np.ndarray:
    """Six orthonormal rigid-body modes: 3 unit translations, 3 rotations ``e_a x r``.

    Returns a ``(6, ndof)`` array; each mode also carries the corresponding nodal rotation.
    """
    n = 6 * len(coords)
    modes = []
    for a in range(3):
        v = np.zeros(n)
        v[a::6] = 1.0
        modes.append(v)
    rot_fields = [
        (0.0, -coords[:, 2], coords[:, 1]),
        (coords[:, 2], 0.0, -coords[:, 0]),
        (-coords[:, 1], coords[:, 0], 0.0),
    ]
    for a, field in enumerate(rot_fields):
        v = np.zeros(n)
        v[0::6] = field[0]
        v[1::6] = field[1]
        v[2::6] = field[2]
        v[3 + a :: 6] = 1.0
        modes.append(v)
    out = []
    for v in modes:  # Gram-Schmidt so the penalty is a clean projector
        for u in out:
            v = v - u * (u @ v)
        out.append(v / np.linalg.norm(v))
    return np.asarray(out)


def _low_rank_penalty(modes: np.ndarray):
    """The penalty matrix ``sum_a phi_a phi_a^T``, node-block sparse."""
    blocks = [
        csr_matrix(modes[:, 6 * i : 6 * i + 6].T @ modes[:, 6 * i : 6 * i + 6])
        for i in range(modes.shape[1] // 6)
    ]
    return block_diag(blocks, format="csr")


def _solve_rigid_removed(K, f: np.ndarray, modes: np.ndarray) -> np.ndarray:
    """Solve ``K u = f`` with the rigid-body modes removed by the **exact** constraint ``Phi^T u = 0``.

    The plan's penalty ``K + k Phi Phi^T`` is the ``k -> 0`` limit of this saddle-point constraint.
    It was measured here that the penalty is *not* exactly solution-neutral: the assembled K
    annihilates the analytic rigid modes only to round-off, so the twist rate drifts by 2.0e-5
    relative at ``k = 1e-6 * median(diag K)`` (isotropic) and by up to 3e-7 for the laminate, above
    the ``1e-9`` the plan asks for. The exact constraint is k-free and is used for the physics; the
    penalty is still run and reported by the tests so the drift is on the record.
    """
    n = K.shape[0]
    m = modes.shape[0]
    augmented = bmat(
        [[K, coo_matrix(modes.T)], [coo_matrix(modes), coo_matrix((m, m))]], format="csc"
    )
    rhs = np.concatenate([f, np.zeros(m)])
    sol = splu(augmented).solve(rhs)
    return np.asarray(sol[:n])


def _solve_penalty(K, f: np.ndarray, modes: np.ndarray, k: float) -> np.ndarray:
    """The plan's rigid-mode penalty: ``(K + k Phi Phi^T) u = f + k Phi (Phi^T f)``."""
    penalised = (K + k * _low_rank_penalty(modes)).tocsr()
    rhs = f + k * (modes.T @ (modes @ f))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return spsolve(penalised, rhs)


def _solve_clamped(K, f: np.ndarray, clamped: list[int]) -> np.ndarray:
    """Sparse direct solve with a Dirichlet BC (the clamped-root engineering case)."""
    n = K.shape[0]
    free_mask = np.ones(n, dtype=bool)
    free_mask[clamped] = False
    free = np.where(free_mask)[0]
    u = np.zeros(n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u[free] = spsolve(K[np.ix_(free, free)], f[free])
    return u


def _clamped_dofs(rings: list[list[int]]) -> list[int]:
    return [6 * nd + d for nd in rings[0] for d in range(6)]


# ─────────────────────────────────────────────────────────────────────────────
# Loads
# ─────────────────────────────────────────────────────────────────────────────


def _realised_torque(coords: np.ndarray, f: np.ndarray) -> float:
    """The torque ruler: ``sum(x_i F_yi - y_i F_xi)`` from the assembled load vector."""
    return float(np.sum(coords[:, 0] * f[1::6] - coords[:, 1] * f[0::6]))


def _ring_shear_flow(coords: np.ndarray, ring: list[int], torque: float) -> np.ndarray:
    """Uniform shear flow ``q = T / 2A`` along one ring's edges (the Saint-Venant torsion)."""
    f = np.zeros(6 * len(coords))
    q = torque / (2.0 * ENCLOSED_AREA)
    n = len(ring)
    for i in range(n):
        j = (i + 1) % n
        a, b = ring[i], ring[j]
        pa, pb = coords[a], coords[b]
        ell = float(np.hypot(pb[0] - pa[0], pb[1] - pa[1]))
        tangent = (pb[:2] - pa[:2]) / ell
        for nd in (a, b):
            f[6 * nd] += 0.5 * q * ell * tangent[0]
            f[6 * nd + 1] += 0.5 * q * ell * tangent[1]
    return f


def _self_equilibrated_load(coords: np.ndarray, rings: list[list[int]]) -> np.ndarray:
    """``+T`` at the tip ring and ``-T`` at the root ring: a fully self-equilibrated torsion."""
    return _ring_shear_flow(coords, rings[-1], TORQUE) + _ring_shear_flow(coords, rings[0], -TORQUE)


def _wall_traction(coords: np.ndarray, rings: list[list[int]], n_ring: int,
                   x_wall: float, traction: float) -> np.ndarray:
    """Uniform traction ``traction`` [N/m] in +y along the tip edges of the wall at ``x = x_wall``."""
    f = np.zeros(6 * len(coords))
    tip = rings[-1]
    for i in range(n_ring):
        j = (i + 1) % n_ring
        a, b = tip[i], tip[j]
        if abs(coords[a, 0] - x_wall) < 1e-12 and abs(coords[b, 0] - x_wall) < 1e-12:
            ell = float(np.hypot(coords[b, 0] - coords[a, 0], coords[b, 1] - coords[a, 1]))
            for nd in (a, b):
                f[6 * nd + 1] += 0.5 * traction * ell
    return f


def _generalised_resultants(modes: np.ndarray, f: np.ndarray) -> np.ndarray:
    """The six rigid-mode resultants of ``f``, normalised by ``|f|`` (a dimensionless number)."""
    return np.asarray([abs(modes[a] @ f) / np.linalg.norm(f) for a in range(modes.shape[0])])


# ─────────────────────────────────────────────────────────────────────────────
# Metrics, distortion, and the hand-written reference
# ─────────────────────────────────────────────────────────────────────────────


def _theta_fit(coords: np.ndarray, u: np.ndarray, ring: list[int]) -> float:
    """In-plane rigid rotation of the ring from its translations [rad]."""
    xs = coords[ring, 0]
    ys = coords[ring, 1]
    ux = u[6 * np.asarray(ring)]
    uy = u[6 * np.asarray(ring) + 1]
    ux = ux - ux.mean()
    uy = uy - uy.mean()
    return float(np.sum(xs * uy - ys * ux) / np.sum(xs * xs + ys * ys))


def _theta_z(u: np.ndarray, ring: list[int]) -> float:
    """Mean nodal rotation about the tube axis (DOF 5) over the ring [rad]."""
    return float(np.mean([u[6 * nd + 5] for nd in ring]))


def _parallelogram(coords: np.ndarray, u: np.ndarray, ring: list[int]) -> tuple[float, float, float]:
    """Section-distortion amplitude of one ring: ``(alpha, beta, deviation)``.

    ``alpha = du_x/dy`` and ``beta = du_y/dx`` are least-squares fits of the mean-removed ring field
    to ``u_x = alpha y`` and ``u_y = beta x`` (the section parallelogram); ``deviation`` is the
    relative norm of the field after subtracting the best-fit rigid rotation.
    """
    xs = coords[ring, 0]
    ys = coords[ring, 1]
    ux = u[6 * np.asarray(ring)]
    uy = u[6 * np.asarray(ring) + 1]
    ux = ux - ux.mean()
    uy = uy - uy.mean()
    alpha = float(np.linalg.lstsq(ys[:, None], ux, rcond=None)[0][0])
    beta = float(np.linalg.lstsq(xs[:, None], uy, rcond=None)[0][0])
    theta = float(np.sum(xs * uy - ys * ux) / np.sum(xs * xs + ys * ys))
    residual = np.sqrt(np.sum((ux + theta * ys) ** 2 + (uy - theta * xs) ** 2))
    norm = float(np.sqrt(np.sum(ux**2 + uy**2)))
    return alpha, beta, residual / norm if norm > 0.0 else 0.0


def _lsq_slope(z: np.ndarray, y: np.ndarray, z_lo: float, z_hi: float) -> float:
    m = (z >= z_lo - 1e-9) & (z <= z_hi + 1e-9)
    return float(np.polyfit(z[m], y[m], 1)[0])


def _bredt_isotropic(thickness: float = THICKNESS) -> tuple[float, float]:
    """Hand-written Bredt: ``GJ = 4 A^2 (G t) / perimeter``, ``theta' = T / GJ``."""
    a66 = G_ISO * thickness
    gj = 4.0 * ENCLOSED_AREA**2 * a66 / PERIMETER
    return gj, TORQUE / gj


def _bredt_laminate() -> tuple[float, float, float, float]:
    """Hand-written ``Qbar66`` and Bredt for the balanced symmetric laminate (independent CLT)."""
    nu21 = NU12 * E2 / E1
    den = 1.0 - NU12 * nu21
    q11 = E1 / den
    q22 = E2 / den
    q12 = NU12 * E2 / den
    q66 = G12
    th = np.deg2rad(45.0)
    s, c = float(np.sin(th)), float(np.cos(th))
    qbar66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * s**2 * c**2 + q66 * (c**4 + s**4)
    a66 = LAM_TOTAL_THICKNESS * qbar66
    gj = 4.0 * ENCLOSED_AREA**2 * a66 / PERIMETER
    return gj, a66, qbar66, TORQUE / gj


def _rate_report(label: str, coords: np.ndarray, u: np.ndarray, rings: list[list[int]],
                 ref_rate: float, window: tuple[float, float]) -> dict[str, float]:
    """Print both rate metrics over ``window`` and return them."""
    z = np.asarray([coords[r[0], 2] for r in rings])
    theta_fit = np.asarray([_theta_fit(coords, u, r) for r in rings])
    theta_z = np.asarray([_theta_z(u, r) for r in rings])
    m = {
        "slope_fit": _lsq_slope(z, theta_fit, *window),
        "slope_z": _lsq_slope(z, theta_z, *window),
    }
    m["ratio_fit"] = m["slope_fit"] / ref_rate
    m["ratio_z"] = m["slope_z"] / ref_rate
    m["metric_diff"] = abs(m["slope_fit"] - m["slope_z"]) / abs(m["slope_z"])
    print(f"\n[{label}] hand-written reference theta' = T/GJ = {ref_rate:.6e} rad/m")
    print(f"  theta_fit rate = {m['slope_fit']:.6e} rad/m  ratio = {m['ratio_fit']:.6f}")
    print(f"  theta_z   rate = {m['slope_z']:.6e} rad/m  ratio = {m['ratio_z']:.6f}")
    print(f"  metric difference |fit - z| / |z| = {m['metric_diff']:.4%}")
    return m


# ─────────────────────────────────────────────────────────────────────────────
# 1-2. Self-equilibrated (boundary-layer-free) torsion - the primary validation
# ─────────────────────────────────────────────────────────────────────────────


def _self_equilibrated_validation(
    label: str, elem_type: int, prop: dict, ref_rate: float
) -> dict[str, float]:
    """Solve and report the two self-equilibrated cases, and return their metrics.

    The machinery and the properties stay here and keep asserting exactly what they asserted before:
    the load rulers, the self-equilibrium of the load, the removal of the rigid modes, the agreement
    of the two metrics and the linearity of the profile. None of them is a comparison against an
    independent reference, so none of them is a store row, and the store reads none of them.

    The comparison against the Bredt rate is *not* here: it is written in the two tests, where the
    store can read it. See `docs/adding-validation-tests.md`.
    """
    coords, conn, rings, n_ring = _tube_mesh()
    asm, K = _assemble(coords, conn, elem_type, prop)

    f_tip = _ring_shear_flow(coords, rings[-1], TORQUE)
    f_root = _ring_shear_flow(coords, rings[0], -TORQUE)
    f = f_tip + f_root
    assert abs(_realised_torque(coords, f_tip) - TORQUE) / TORQUE < 1e-9, "tip shear-flow torque"
    assert abs(_realised_torque(coords, f_root) + TORQUE) / TORQUE < 1e-9, "root shear-flow torque"

    modes = _rigid_modes(coords)
    resultants = _generalised_resultants(modes, f)
    print(f"\n[{label}] self-equilibrated: generalised resultants (normalised) = "
          f"{' '.join(f'{r:.2e}' for r in resultants)}")
    assert resultants.max() < 1e-9, f"load is not self-equilibrated: {resultants}"

    u = _solve_rigid_removed(K, f, modes)
    assert np.linalg.norm(modes @ u) / np.linalg.norm(u) < 1e-10, "rigid modes not removed"

    z = np.asarray([coords[r[0], 2] for r in rings])
    theta_fit = np.asarray([_theta_fit(coords, u, r) for r in rings])
    theta_z = np.asarray([_theta_z(u, r) for r in rings])
    m = _rate_report(label, coords, u, rings, ref_rate, (0.4 * L, 0.9 * L))

    # Linearity: the exact solution is uniform torsion, so the two halves must have the same slope.
    half_a = _lsq_slope(z, theta_fit, 0.1 * L, 0.5 * L)
    half_b = _lsq_slope(z, theta_fit, 0.5 * L, 0.9 * L)
    linearity = abs(half_a - half_b) / abs(half_b)
    print(f"  linearity 0.1-0.5L {half_a:.6e} vs 0.5-0.9L {half_b:.6e} -> {linearity:.4%}")

    # The plan's prescribed penalty, run and reported (see _solve_rigid_removed).
    k_bad = _PENALTY_PRESCRIBED * np.median(K.diagonal())
    u_k = _solve_penalty(K, f, modes, k_bad)
    u_1000k = _solve_penalty(K, f, modes, 1000.0 * k_bad)
    rate_k = _lsq_slope(z, np.asarray([_theta_fit(coords, u_k, r) for r in rings]), 0.4 * L, 0.9 * L)
    rate_1000k = _lsq_slope(
        z, np.asarray([_theta_fit(coords, u_1000k, r) for r in rings]), 0.4 * L, 0.9 * L
    )
    print(f"  prescribed penalty k=1e-6*median(diagK)={k_bad:.3e}: rate drifts "
          f"{abs(rate_k - rate_1000k) / abs(rate_k):.2e} relative to 1000k "
          f"(the exact constraint is used instead)")

    print("  theta_fit(z) profile:")
    for k in range(0, len(z), max(1, len(z) // 8)):
        print(f"    z={z[k]:.2f}  theta_fit={theta_fit[k]: .6e}  theta_z={theta_z[k]: .6e}")

    assert m["metric_diff"] < METRIC_AGREEMENT, (
        f"{label} self-equilibrated: theta_fit and theta_z differ by {m['metric_diff']:.4%}; in the "
        f"exact Saint-Venant solution they must be the same section rotation"
    )
    assert linearity < SELF_CONSISTENCY, (
        f"{label} self-equilibrated: the 0.1-0.5L and 0.5-0.9L slopes differ by {linearity:.4%}"
    )
    return m


def test_self_equilibrated_isotropic_reproduces_bredt():
    """Boundary-layer-free torsion, isotropic MITC4: both metrics equal ``T/GJ``.

    ``+T`` at the tip ring and ``-T`` at the root ring (both as ``q = T/(2A)``) is fully
    self-equilibrated - all six rigid-mode resultants vanish - so the solution is uniform torsion
    with no clamped root and no boundary layer. In the exact solution the section does not distort
    in plane, so ``theta_fit`` and ``theta_z`` must agree and both must equal the section rotation.
    """
    _, ref = _bredt_isotropic()
    m = _self_equilibrated_validation("iso / self-equilibrated", 4, _iso_prop(), ref)
    assert_relative_error(
        m["ratio_fit"],
        1.0,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt T L / GJ, the closed form for a closed thin-walled tube",
        what="isotropic theta_fit rate on the Bredt rate",
    )
    assert_relative_error(
        m["ratio_z"],
        1.0,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt T L / GJ, the closed form for a closed thin-walled tube",
        what="isotropic theta_z rate on the Bredt rate",
    )


def test_self_equilibrated_laminate_reproduces_bredt():
    """Boundary-layer-free torsion, MITC4Composite ``[45,-45]s``: both metrics equal ``T/GJ``."""
    _, prop = _laminate_and_prop()
    _, _, _, ref = _bredt_laminate()
    m = _self_equilibrated_validation("laminate / self-equilibrated", 44, prop, ref)
    assert_relative_error(
        m["ratio_fit"],
        1.0,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt T L / GJ, the closed form for a closed thin-walled tube",
        what="laminate theta_fit rate on the Bredt rate",
    )
    assert_relative_error(
        m["ratio_z"],
        1.0,
        tol=TOL,
        kind="analytical",
        reference_name="Bredt T L / GJ, the closed form for a closed thin-walled tube",
        what="laminate theta_z rate on the Bredt rate",
    )


# ─────────────────────────────────────────────────────────────────────────────
# 3. Clamped-root case as the quantified boundary-layer comparison
# ─────────────────────────────────────────────────────────────────────────────


def test_clamped_root_metric_convergence_identifies_the_section_rotation():
    """Clamped root + tip shear flow: which metric is the section rotation, and how long is the layer.

    The two metrics are reported for ``L = 6, 12, 24 m`` and ``t = 10, 0.5 mm`` over the proportional
    window ``0.4L-0.9L`` and the fixed window ``2.4-5.4 m``, with the local secant slope along z.
    Only ``theta_z`` converges: it is the section rotation. ``theta_fit`` is polluted by the
    clamped-root boundary layer and by the thickness-dependent section in-plane shear, so it is
    reported, not asserted.
    """
    ref_10, ref_05 = _bredt_isotropic(0.01)[1], _bredt_isotropic(0.0005)[1]
    print("\n================ clamped-root convergence: theta_fit vs theta_z ================")
    for thickness, ref_rate in ((0.01, ref_10), (0.0005, ref_05)):
        prop = _iso_prop(thickness)
        print(f"\n--- t = {thickness * 1e3:.1f} mm   T/GJ = {ref_rate:.6e} rad/m")
        for length in (6.0, 12.0, 24.0):
            n_z = int(round(length / DZ))
            coords, conn, rings, _ = _tube_mesh(length=length, n_z=n_z)
            _, K = _assemble(coords, conn, 4, prop)
            f = _ring_shear_flow(coords, rings[-1], TORQUE)
            u = _solve_clamped(K, f, _clamped_dofs(rings))
            z = np.asarray([coords[r[0], 2] for r in rings])
            theta_fit = np.asarray([_theta_fit(coords, u, r) for r in rings])
            theta_z = np.asarray([_theta_z(u, r) for r in rings])
            fit_prop = _lsq_slope(z, theta_fit, 0.4 * length, 0.9 * length) / ref_rate
            z_prop = _lsq_slope(z, theta_z, 0.4 * length, 0.9 * length) / ref_rate
            fit_abs = _lsq_slope(z, theta_fit, 2.4, 5.4) / ref_rate
            z_abs = _lsq_slope(z, theta_z, 2.4, 5.4) / ref_rate
            print(f"  L={length:4.0f} m: 0.4L-0.9L fit={fit_prop:.4f} z={z_prop:.4f}   "
                  f"| 2.4-5.4 m fit={fit_abs:.4f} z={z_abs:.4f}")
            local = np.asarray([(theta_fit[k + 1] - theta_fit[k]) / (length / n_z) / ref_rate
                                for k in range(n_z)])
            stations = np.arange(0, n_z, max(1, n_z // 8))
            print("        local theta_fit secant slope/ref: "
                  + " ".join(f"{((k + 0.5) * length / n_z):.1f}:{(local[k]):.3f}" for k in stations))
            if length == 24.0:
                asymptotic_z = z_prop
        assert_relative_error(
            asymptotic_z,
            1.0,
            tol=TOL,
        kind="analytical",
            reference_name=(
                "the idealised limit: at L=24 m the section rotation rate is T/GJ, so the "
                "theta_z rate on it must reach 1"
            ),
            what=f"t={thickness * 1e3:.1f} mm asymptotic theta_z rate at L=24 m",
        )
    print("=================================================================================")


# ─────────────────────────────────────────────────────────────────────────────
# 4. One-wall traction - a statically equivalent load that is not equivalent
# ─────────────────────────────────────────────────────────────────────────────


def test_wall_traction_excites_a_section_distortion():
    """One wall's tip traction: same torque, but the response is a section-distortion mode.

    The plan's case 2. The resultant torque equals the shear-flow case's, but the load is applied
    on a single wall, so it is not the Saint-Venant traction: it excites a soft in-plane section
    parallelogram (``u_x = alpha y``, ``u_y = beta x``) that the shear flow does not. The test
    asserts what the physics pins - the mirror (equal magnitude, opposite sign), the torque ruler,
    and that the distortion amplitude is an order of magnitude larger than the shear flow's - and
    reports the twist-rate ratios by both metrics, the energy ratio and the tip-displacement ratio.
    """
    tau = 2.0 * TORQUE / (B * H)  # so that T_eff = tau * h * b / 2 = TORQUE
    intended = tau * H * B / 2.0
    gj, _ = _bredt_isotropic()

    coords, conn, rings, n_ring = _tube_mesh()
    _, K = _assemble(coords, conn, 4, _iso_prop())
    clamped = _clamped_dofs(rings)

    f_right = _wall_traction(coords, rings, n_ring, +B / 2.0, tau)
    f_left = _wall_traction(coords, rings, n_ring, -B / 2.0, tau)
    t_right = _realised_torque(coords, f_right)
    t_left = _realised_torque(coords, f_left)
    assert abs(t_right - intended) / abs(intended) < 1e-9, "right-wall torque ruler"
    assert abs(t_left + intended) / abs(intended) < 1e-9, "left-wall torque ruler (mirror)"

    # Clamped root, as in the plan's case 2: the engineering comparison, not the boundary-layer-free one.
    u_right = _solve_clamped(K, f_right, clamped)
    u_left = _solve_clamped(K, f_left, clamped)
    u_flow = _solve_clamped(K, _ring_shear_flow(coords, rings[-1], TORQUE), clamped)

    z = np.asarray([coords[r[0], 2] for r in rings])
    m_right = _rate_report("wall traction / right wall", coords, u_right, rings, t_right / gj,
                           (0.4 * L, 0.9 * L))
    m_left = _rate_report("wall traction / left wall (mirror)", coords, u_left, rings, t_left / gj,
                          (0.4 * L, 0.9 * L))
    m_flow = _rate_report("shear flow reference (tip only)", coords, u_flow, rings, TORQUE / gj,
                          (0.4 * L, 0.9 * L))

    # Mirror: a physics claim, so it is checked before any magnitude one.
    slope_fit_right = m_right["slope_fit"]
    slope_fit_left = m_left["slope_fit"]
    assert abs(abs(slope_fit_left) - abs(slope_fit_right)) / abs(slope_fit_right) < 0.01, (
        f"mirror magnitude: {abs(slope_fit_left):.6e} vs {abs(slope_fit_right):.6e}"
    )
    assert np.sign(slope_fit_left) != np.sign(slope_fit_right), "mirror sign"

    # Section-distortion amplitude in the interior, compared with the shear flow.
    station = len(rings) // 2
    alpha_w, beta_w, dev_w = _parallelogram(coords, u_right, rings[station])
    alpha_f, beta_f, dev_f = _parallelogram(coords, u_flow, rings[station])
    shear_w = 0.5 * (alpha_w + beta_w)
    shear_f = 0.5 * (alpha_f + beta_f)
    distortion_ratio = abs(shear_w) / abs(shear_f)
    print(f"\n[wall traction] at z={z[station]:.2f}: alpha={alpha_w:+.4e} beta={beta_w:+.4e} "
          f"shear={shear_w:+.4e} deviation-from-rigid={dev_w:.3f}")
    print(f"[shear flow]   at z={z[station]:.2f}: alpha={alpha_f:+.4e} beta={beta_f:+.4e} "
          f"shear={shear_f:+.4e} deviation-from-rigid={dev_f:.3f}")
    print(f"  distortion amplitude ratio (wall / shear-flow) = {distortion_ratio:.1f}x")

    energy_ratio = float(f_right @ u_right) / float(_ring_shear_flow(coords, rings[-1], TORQUE)
                                                     @ u_flow)
    tip_right = float(np.abs(u_right[1::6]).max())
    tip_flow = float(np.abs(u_flow[1::6]).max())
    print(f"  elastic-energy ratio = {energy_ratio:.2f}x   tip in-plane displacement ratio = "
          f"{tip_right / tip_flow:.2f}x")
    print(f"  reported twist-rate ratios: right wall fit={m_right['ratio_fit']:.4f} "
          f"z={m_right['ratio_z']:.4f}; shear flow fit={m_flow['ratio_fit']:.4f} "
          f"z={m_flow['ratio_z']:.4f}")

    assert distortion_ratio > 10.0, (
        f"the one-wall traction should excite a section distortion at least an order of magnitude "
        f"larger than the shear flow; measured {distortion_ratio:.1f}x"
    )
