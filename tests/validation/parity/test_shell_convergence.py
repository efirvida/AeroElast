"""MITC4 shell convergence-order study.

This module holds two studies:

1. ``test_in_plane_bending_convergence`` -- the observed convergence order of a
   single, clean, mesh-independent MITC4 problem (CalculiX-free).
2. ``test_composite_laminate_gap_mesh_study`` -- whether the composite laminate
   mismatch against CalculiX S8R is a mesh artifact or a formulation/ABD
   difference (needs CalculiX).

Study 1 answers a question the rest of the suite cannot: is a tip-displacement
mismatch caused by a wrong element formulation, or simply by a mesh that is too
coarse?  It measures the observed convergence order of a single, clean,
mesh-independent problem.

Geometry (plate in the XY plane, mid-surface at Z = 0):

- length      ``L = 1.0``   m along Y (cantilever axis)
- width       ``B = 0.1``   m along X
- thickness   ``t = 0.005`` m along Z
- clamped on the whole ``y = 0`` edge (all 6 DOFs)
- total load  ``P = 100``   N in ``+X``, split equally over every node of the
  free ``y = L`` edge, so the total force is mesh-independent

Analytical reference -- in-plane bending of the strip, load along X:

    ``delta = P * L**3 / (3 * E * I)``,  with ``I = t * B**3 / 12``

For the constants below ``delta ~= 1.2308e-3`` m (about 1230.8 um).

The load is deliberately spread over the free edge: a single-node point load
makes the local displacement non-convergent and would corrupt the rate.

The FE solution converges to the *Timoshenko* beam value, not the
Euler-Bernoulli one: Euler-Bernoulli omits transverse-shear deformation, so the
residual against it has a physical floor and stops decreasing once the
discretization error drops below that floor.  Measuring the order from those raw
residuals is therefore meaningless -- the pairwise order collapses to ~0 as both
errors sit on the floor.  The order is instead measured by *self-convergence*
(Richardson): for a sequence with ``h`` halving and ``u_h = u_exact + C h**p``
the successive differences satisfy
``(u_n - u_{n-1}) / (u_{n-1} - u_{n-2}) = 2**-p``, hence
``p = log2((u_{n-1} - u_{n-2}) / (u_n - u_{n-1}))``, which cancels the unknown
reference floor.  The Euler-Bernoulli value is then used only as an independent
cross-check on the Richardson-extrapolated limit.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from _aeroelast import PyMeshAssembler

from tests.conftest import ccx_bin_or_skip

from tests.support.assertions import assert_residual_below  # noqa: E402

# ---------------------------------------------------------------------------
# Geometry and material
# ---------------------------------------------------------------------------
L = 1.0  # length along Y (m)
B = 0.1  # width along X (m)
T = 0.005  # shell thickness along Z (m)

E = 65e9  # isotropic equivalent modulus (the same e_equiv the parity tests use)
NU = 0.3
P = 100.0  # total force in +X (N)

# Analytical cantilever reference: delta = P * L**3 / (3 * E * I)
I = T * B**3 / 12.0
DELTA = P * L**3 / (3.0 * E * I)  # ~= 1.2308e-3 m

# Element size h = L / ny halves on each refinement.
MESHES = [(2, 10), (4, 20), (8, 40), (16, 80)]

_MITC4 = 4  # Rust element type code for the 6-DOF-per-node quad shell
_DOFS_PER_NODE = 6

# MITC4 displacement should be ~O(h**2); 1.5 leaves headroom for the coarse
# meshes.
MIN_ORDER = 1.5

# The Richardson-extrapolated limit must be within this relative error of the
# Euler-Bernoulli reference.  The 2% is justified physics, not convenience:
# Euler-Bernoulli omits the transverse-shear deflection.  For this strip the
# Timoshenko shear contribution is delta_shear = P*L/(G*A_s) with
# G = E/(2*(1+nu)) and A_s = (5/6)*B*t:
#   G           = 65e9 / (2 * 1.3)                 = 2.5000e10 Pa
#   A_s         = (5/6) * 0.1 * 0.005              = 4.1667e-4 m^2
#   delta_shear = 100 * 1.0 / (2.5e10 * 4.1667e-4) = 9.6e-6 m ~= 9.6 um
# on ~1230 um of bending, i.e. ~0.8%.  A 2% window therefore contains the
# expected Timoshenko-vs-Euler-Bernoulli physics while still failing loudly if
# the FE solution converges to something else (e.g. a wrong formulation).
EXTRAPOLATED_TOL = 0.02


def _build_strip(nx: int, ny: int):
    """Build a structured MITC4 grid for the plate strip.

    Returns ``(node_coords, connectivity, n_nodes_x)``.  The node index is
    ``j * (nx + 1) + i``, with ``i`` along X and ``j`` along Y.
    """
    n_nodes_x = nx + 1
    xs = np.linspace(0.0, B, n_nodes_x)
    ys = np.linspace(0.0, L, ny + 1)
    node_coords = np.array([[x, y, 0.0] for y in ys for x in xs], dtype=float)
    connectivity = []
    for j in range(ny):
        for i in range(nx):
            n00 = j * n_nodes_x + i
            connectivity.append([n00, n00 + 1, n00 + 1 + n_nodes_x, n00 + n_nodes_x])
    return node_coords, connectivity, n_nodes_x


def _tip_displacement(nx: int, ny: int) -> float:
    """Solve the clamped strip and return the +X tip displacement (m)."""
    node_coords, connectivity, n_nodes_x = _build_strip(nx, ny)
    n_elems = len(connectivity)
    material = {"type": "isotropic", "e": E, "nu": NU, "rho": 0.0, "thickness": T}

    asm = PyMeshAssembler(
        node_coords=node_coords,
        connectivity=connectivity,
        elem_types=[_MITC4] * n_elems,
        materials=[material] * n_elems,
    )
    rows, cols, vals = asm.assemble_k()
    n_dof = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n_dof, n_dof)).tocsr()

    # Total force P in +X, shared equally by every node on the free edge.
    free_edge = [ny * n_nodes_x + i for i in range(n_nodes_x)]
    f = np.zeros(n_dof, dtype=float)
    for node in free_edge:
        f[node * _DOFS_PER_NODE] = P / len(free_edge)

    # Clamp all 6 DOFs on the y = 0 edge.
    fixed = np.array(
        [i * _DOFS_PER_NODE + d for i in range(n_nodes_x) for d in range(_DOFS_PER_NODE)]
    )
    free = np.setdiff1d(np.arange(n_dof), fixed)

    u = np.zeros(n_dof, dtype=float)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])

    # Tip displacement at the mid-width loaded node (beam-axis point).
    tip_node = ny * n_nodes_x + nx // 2
    return float(u[tip_node * _DOFS_PER_NODE])


def _fit_order(h: np.ndarray, error: np.ndarray) -> float:
    """Observed order = least-squares slope of log(error) vs log(h)."""
    return float(np.polyfit(np.log(h), np.log(error), 1)[0])


def test_in_plane_bending_convergence():
    """The in-plane tip displacement converges with observed order >= 1.5.

    The order is measured by self-convergence (Richardson) so the physical
    Timoshenko-vs-Euler-Bernoulli floor in the raw residuals does not corrupt
    it.  The Richardson-extrapolated limit is then cross-checked against the
    Euler-Bernoulli reference within the physical shear window.
    """
    h = np.array([L / ny for _, ny in MESHES], dtype=float)
    displacements = np.array([_tip_displacement(nx, ny) for nx, ny in MESHES])
    errors = np.abs(displacements - DELTA) / DELTA  # vs Euler-Bernoulli, info only

    # Self-convergence: u_h = u_exact + C h**p with h halving gives
    #   (u_n - u_{n-1}) / (u_{n-1} - u_{n-2}) = 2**-p.
    diffs = np.diff(displacements)
    orders = [float(np.log2(diffs[k - 1] / diffs[k])) for k in range(1, len(diffs))]

    # Richardson-extrapolated h -> 0 limit from the finest pair and the finest
    # order estimate.
    p_fine = orders[-1]
    u_extrap = float(
        displacements[-1] + (displacements[-1] - displacements[-2]) / (2.0**p_fine - 1.0)
    )
    extrap_error = abs(u_extrap - DELTA) / DELTA

    print()
    print("In-plane bending convergence (MITC4, Euler-Bernoulli reference for info)")
    print(f"  analytical delta = {DELTA * 1e6:.4f} um  (P*L^3/(3*E*I), I = t*B^3/12)")
    print(
        f"  {'mesh':>10} {'h [m]':>10} {'u_tip [um]':>12} "
        f"{'rel err [%]':>12} {'pairwise order':>15}"
    )
    for k, (nx, ny) in enumerate(MESHES):
        pairwise = ""
        if k > 0:
            pairwise = f"{np.log(errors[k] / errors[k - 1]) / np.log(h[k] / h[k - 1]):.3f}"
        print(
            f"  {f'({nx},{ny})':>10} {h[k]:>10.5f} {displacements[k] * 1e6:>12.4f} "
            f"{errors[k] * 100:>12.4f} {pairwise:>15}"
        )
    print(
        "  (the pairwise order above is floor-limited and shown for information "
        "only; it is NOT asserted)"
    )
    print(f"  raw-error log-log fit (informational): {_fit_order(h, errors):.3f}")
    print(f"  successive differences [um]: {[f'{d * 1e6:.4f}' for d in diffs]}")
    print(f"  self-convergence order, triple 1: {orders[0]:.4f}")
    print(f"  self-convergence order, triple 2: {orders[1]:.4f}")
    agree = abs(orders[0] - orders[1]) <= 0.1
    print(
        f"  the two independent order estimates {'AGREE' if agree else 'DISAGREE'} "
        f"(|delta| = {abs(orders[0] - orders[1]):.4f})"
    )
    print(f"  Richardson-extrapolated limit = {u_extrap * 1e6:.4f} um")
    print(
        f"  extrapolated limit vs Euler-Bernoulli: {extrap_error * 100:.4f}% "
        f"(physical shear window {EXTRAPOLATED_TOL * 100:.1f}%)"
    )

    min_order = min(orders)
    assert min_order >= MIN_ORDER, (
        f"self-convergence order {min_order:.3f} < {MIN_ORDER}: "
        f"estimates={[f'{p:.3f}' for p in orders]}"
    )
    assert_residual_below(
        extrap_error,
        tol=EXTRAPOLATED_TOL,
        kind="analytical",
        reference_name="the Euler-Bernoulli tip deflection of the same cantilever",
        what="Richardson-extrapolated tip deflection",
    )


# ---------------------------------------------------------------------------
# Composite laminate mesh-dependence study (needs CalculiX)
# ---------------------------------------------------------------------------
# h halves between consecutive meshes (ny: 5 -> 10 -> 20 -> 40).
LAMINATE_MESHES = [(2, 5), (4, 10), (8, 20), (16, 40)]
_LAMINATE_LOAD_X = 100.0  # N, +X at the free-edge centre node (same as parity test)

# Project honesty bound (docs/validation-policy.md): any AeroElast-vs-reference
# structural comparison must stay within 5%.  The coarsest laminate mesh is the
# binding case here at 4.18%.
LAMINATE_GAP_TOL = 0.05
# Tight diagnostic bound on the mesh-independent (finest) comparison.  A 2%
# bound turns the printed plateau into a real, falsifiable assertion: it
# localises the residual MITC4+CLT-vs-S8R difference to a bounded band and
# would fail if the gap regressed.  Measured: 1.729% (bound 2.0%).
LAMINATE_FINEST_GAP_TOL = 0.02
# AeroElast's own MITC4+CLT sequence must still converge; observed self-
# convergence orders are 1.573 and 1.555.
LAMINATE_MIN_ORDER = 1.5


def _laminate_aero_tip_x(parity, nx: int, ny: int) -> float:
    """AeroElast MITC4 + CLT ABD tip displacement (+X) for the laminate."""
    mesh = parity._build_composite_plate_mesh(nx=nx, ny=ny)
    node_coords = np.asarray([[n.x, n.y, n.z] for n in mesh.nodes], dtype=float)
    conn = [[mesh.node_id_to_index[nid] for nid in el.node_ids] for el in mesh.elements]
    elem_types = [4] * len(mesh.elements)
    mat_dict = parity._make_laminate_mat(
        parity.E1, parity.E2, parity.G12, parity.nu12, parity.thickness
    )
    mats = [mat_dict] * len(mesh.elements)

    asm = PyMeshAssembler(
        node_coords=node_coords, connectivity=conn, elem_types=elem_types, materials=mats
    )
    rows, cols, vals = asm.assemble_k()
    n = asm.dofs_count
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    f = np.zeros(n, dtype=float)
    center = next(iter(mesh.get_node_set("free_center").nodes.values()))
    i0 = mesh.node_id_to_index[center.id] * 6
    f[i0] = _LAMINATE_LOAD_X

    clamped = {
        mesh.node_id_to_index[n.id] * 6 + d
        for n in mesh.get_node_set("clamped").nodes.values()
        for d in range(6)
    }
    free_mask = np.ones(n, dtype=bool)
    for dof in clamped:
        free_mask[dof] = False
    free = np.where(free_mask)[0]

    u = np.zeros(n, dtype=float)
    u[free] = spsolve(K[np.ix_(free, free)], f[free])
    return float(u[i0])


def _laminate_ccx_tip_x(parity, ccx_bin: str, inp_path, nx: int, ny: int) -> float:
    """CalculiX S8R + *SHELL SECTION, COMPOSITE tip displacement (+X)."""
    mesh = parity._build_composite_plate_mesh(nx=nx, ny=ny)
    parity._write_composite_ccx_inp(inp_path, mesh, (_LAMINATE_LOAD_X, 0.0, 0.0))
    node_ids = [
        mesh.node_id_to_index[n.id] + 1 for n in mesh.get_node_set("free_center").nodes.values()
    ]
    result = parity._run_ccx(inp_path, ccx_bin)
    if result.returncode != 0:
        parity._fail_ccx(result, inp_path)
    frd_path = inp_path.with_suffix(".frd")
    if not frd_path.exists():
        pytest.fail(f"No FRD output from CCX for {inp_path.name}")
    ccx_disp = parity._parse_frd_disp(frd_path, node_ids)
    if not ccx_disp:
        pytest.fail(f"Could not parse CCX displacements from {frd_path.name}")
    return float(ccx_disp[node_ids[0]][0])


def _richardson_orders(values: np.ndarray) -> tuple[list[float], np.ndarray]:
    """Self-convergence orders from an h-halving sequence (Richardson)."""
    diffs = np.diff(values)
    orders = [float(np.log2(diffs[k - 1] / diffs[k])) for k in range(1, len(diffs))]
    return orders, diffs


def _richardson_limit(values: np.ndarray, order: float) -> float:
    """Richardson-extrapolated h -> 0 limit: u_n + (u_n - u_{n-1})/(2**p - 1)."""
    return float(values[-1] + (values[-1] - values[-2]) / (2.0**order - 1.0))


def _extrapolatable(diffs: np.ndarray, orders: list[float]) -> bool:
    """True when a monotone sequence supports a Richardson extrapolation."""
    return (
        bool(np.all(diffs > 0))
        and len(orders) >= 1
        and all(np.isfinite(o) and o > 0.5 for o in orders)
    )


def test_composite_laminate_gap_mesh_study(tmp_path):
    """Is the AeroElast-vs-CalculiX laminate gap a mesh artifact or a
    formulation/ABD difference?

    The gap is measured on the SAME 8-ply [0/90/45/-45]s laminate on both
    sides, across four meshes with h halving.  If the gap shrinks toward zero
    with refinement it is a discretization artifact; if it plateaus it is a
    formulation / ABD-handling difference.

    The measurement is now an assertion, not a print: AeroElast's own sequence
    must keep converging (order >= ``LAMINATE_MIN_ORDER``), every mesh gap must
    stay inside the 5% project bound, the mesh-independent (finest) gap must
    stay inside the tight 2% diagnostic bound, and refinement must not make the
    coarse mismatch worse.  The plateau/shrink verdict is printed as the
    diagnosis that justifies the bound; it is deliberately not asserted, so a
    future formulation fix that shrinks the gap cannot fail this test.
    """
    ccx_bin = ccx_bin_or_skip()
    # Imported lazily: the parity module owns the laminate construction, and
    # its petsc4py guard must not skip this module's CalculiX-free test.
    parity = pytest.importorskip(
        "tests.validation.parity.test_composite_beam_parity",
        reason="composite parity helpers unavailable",
    )

    h = np.array([L / ny for _, ny in LAMINATE_MESHES], dtype=float)
    aero = np.array([_laminate_aero_tip_x(parity, nx, ny) for nx, ny in LAMINATE_MESHES])
    ccx = np.array(
        [
            _laminate_ccx_tip_x(parity, ccx_bin, tmp_path / f"laminate_{nx}x{ny}.inp", nx, ny)
            for nx, ny in LAMINATE_MESHES
        ]
    )
    gaps = (aero - ccx) / ccx

    print()
    print("Composite laminate [0/90/45/-45]s bending: AeroElast vs CalculiX S8R")
    print(
        f"  {'mesh':>10} {'h [m]':>10} {'aero [um]':>12} {'ccx [um]':>12} "
        f"{'gap [%]':>10} {'gap ratio':>10}"
    )
    for k, (nx, ny) in enumerate(LAMINATE_MESHES):
        ratio = "" if k == 0 else f"{gaps[k] / gaps[k - 1]:.3f}"
        print(
            f"  {f'({nx},{ny})':>10} {h[k]:>10.5f} {aero[k] * 1e6:>12.4f} "
            f"{ccx[k] * 1e6:>12.4f} {gaps[k] * 100:>10.4f} {ratio:>10}"
        )

    aero_orders, aero_diffs = _richardson_orders(aero)
    ccx_orders, ccx_diffs = _richardson_orders(ccx)
    print(f"  AeroElast self-convergence orders: {[f'{p:.3f}' for p in aero_orders]}")
    print(f"  CalculiX  self-convergence orders: {[f'{p:.3f}' for p in ccx_orders]}")

    if _extrapolatable(aero_diffs, aero_orders) and _extrapolatable(ccx_diffs, ccx_orders):
        aero_lim = _richardson_limit(aero, aero_orders[-1])
        ccx_lim = _richardson_limit(ccx, ccx_orders[-1])
        gap_lim = (aero_lim - ccx_lim) / ccx_lim
        print(f"  Richardson limit  aero = {aero_lim * 1e6:.4f} um, ccx = {ccx_lim * 1e6:.4f} um")
        print(
            f"  Richardson-extrapolated gap = {gap_lim * 100:.4f}% "
            f"(finest-mesh gap = {gaps[-1] * 100:.4f}%)"
        )
        trend = (
            "SHRINKS toward zero -> mesh artifact"
            if abs(gap_lim) < 0.5 * abs(gaps[-1])
            else "PLATEAUS -> formulation / ABD"
        )
    else:
        print(
            "  Richardson extrapolation not supported by this data "
            "(non-monotone or noisy sequence); reporting raw trend only."
        )
        # The raw trend still answers the decisive question: a gap that keeps
        # shrinking toward zero is a mesh artifact; a gap whose successive
        # ratios approach 1 has plateaued, i.e. a formulation/ABD difference.
        ratio_fine = abs(gaps[-1] / gaps[-2]) if len(gaps) >= 2 else float("inf")
        trend = (
            "SHRINKS toward zero -> mesh artifact"
            if ratio_fine < 0.75
            else "PLATEAUS -> formulation / ABD"
        )

    print(f"  gap trend with refinement: {[f'{g * 100:.4f}%' for g in gaps]}")
    print(f"  DECISIVE ANSWER: {trend}")

    assert np.all(np.isfinite(aero)) and np.all(aero > 0), "non-physical AeroElast tip"
    assert np.all(np.isfinite(ccx)) and np.all(ccx > 0), "non-physical CalculiX tip"

    min_aero_order = min(aero_orders)
    assert min_aero_order >= LAMINATE_MIN_ORDER, (
        f"AeroElast laminate self-convergence order {min_aero_order:.3f} < "
        f"{LAMINATE_MIN_ORDER}: estimates={[f'{p:.3f}' for p in aero_orders]}"
    )

    assert np.all(np.abs(gaps) < LAMINATE_GAP_TOL), (
        f"AeroElast-vs-CalculiX laminate gap exceeds the "
        f"{LAMINATE_GAP_TOL * 100:.0f}% project bound: "
        f"{[f'{g * 100:.4f}%' for g in gaps]}"
    )

    finest_gap = abs(gaps[-1])
    assert_residual_below(
        finest_gap,
        tol=LAMINATE_FINEST_GAP_TOL,
        kind="code",
        reference_name="CalculiX 2.23, the same laminate on the same meshes",
        what="finest-mesh laminate gap",
    )

    # Refinement must not make the coarse-to-fine mismatch worse.  The 4.18%
    # coarse gap is largely a mesh artifact; the finest gap is the residual.
    # (Satisfied both when the residual plateaus and when a future fix shrinks
    # it toward zero.)
    assert abs(gaps[-1]) < abs(gaps[0]), (
        f"refinement increased the laminate gap: coarsest {gaps[0] * 100:.4f}% -> "
        f"finest {gaps[-1] * 100:.4f}%"
    )
