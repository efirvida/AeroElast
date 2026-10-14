"""CalculiX S8R as an independent code for the Ko et al. (2017) hook and hyperbolic paraboloid.

The references of both benchmarks are Ko's own MITC9 solutions (hook 4.82482 at N = 64; paraboloid
2.8780e-4 and 2.3856e-4 on 72 x 144 elements). CalculiX S8R is a different code and a different
shell formulation, so a refined S8R answer on the *same model the live tests build* says whether
those references are the converged answers of that model.
"""

from __future__ import annotations

import math

import pytest

from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_relative_error
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck

TOL = 0.05  # suite rule

# --- Raasch hook -----------------------------------------------------------------------------
R1, R2, ANGLE1, ANGLE2 = 14.0, 46.0, math.radians(60.0), math.radians(150.0)
HOOK_WIDTH = 20.0
HOOK_THICKNESS = 2.0
HOOK_E = 3.3e3
HOOK_WREF = 4.82482  # Ko et al. 2017, MITC9 N = 64
_L2, _L1 = R2 * ANGLE2, R1 * ANGLE1


def _hook_centre_line(s: float) -> tuple[float, float]:
    """(x, z) of the hook mid-line at arc length ``s``, main arc first (RaaschHookMesh)."""
    if s <= _L2:
        phi = s / R2
        return -R2 + R2 * math.cos(phi), R2 * math.sin(phi)
    x1, z1 = -R2 + R2 * math.cos(ANGLE2), R2 * math.sin(ANGLE2)
    cx, cz = x1 + (x1 + R2) / R2 * R1, z1 + z1 / R2 * R1
    psi = (s - _L2) / R1
    ang = ANGLE2 + math.pi - psi
    return cx + R1 * math.cos(ang), cz + R1 * math.sin(ang)


def _hook_point(u: float, v: float):
    x, z = _hook_centre_line(u * (_L1 + _L2))
    return (x, v * HOOK_WIDTH, z)


def _hook_displacement(tmp_path, ccx_bin: str, n_width: int, poisson: float) -> float:
    """Mean Y displacement of the loaded edge under a unit total Y shear.

    ``RaaschHookMesh`` names its node sets the other way round from the geometry it draws: its
    ``root`` set (clamped by the live test) is the end of the small R = 14 arc, at
    ``(-97.96, 0, 16)``, and its ``tip`` set (loaded) is the origin where the large R = 46 arc
    starts. That is the benchmark's hook (clamped at the small arc, as in the Abaqus Benchmarks
    Manual), so this deck clamps ``u = 1`` and loads ``u = 0``.
    """
    n_curve = 6 * n_width
    tip_nodes = [node_id(0, j, n_curve) for j in range(2 * n_width + 1)]
    # Simpson weights of a quadratic edge: 1/6, 4/6, 2/6, ... per unit element length.
    loads = {}
    for j, nid in enumerate(tip_nodes):
        w = 1.0 if j in (0, 2 * n_width) else (4.0 if j % 2 else 2.0)
        loads[nid] = (2, w / (6.0 * n_width))
    inp = tmp_path / f"hook_{n_width}_{int(poisson * 100)}.inp"
    write_s8r_surface_deck(
        inp,
        _hook_point,
        n_curve,
        n_width,
        young=HOOK_E,
        poisson=poisson,
        thickness=HOOK_THICKNESS,
        boundary=lambda x, y, z: (
            [1, 2, 3, 4, 5, 6] if abs(x + 97.96152423) < 1e-6 and abs(z - 16.0) < 1e-6 else []
        ),
        nodal_loads=loads,
        print_nodes=tip_nodes,
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    return sum(abs(float(disp[n][1])) for n in tip_nodes) / len(tip_nodes)


@pytest.mark.parametrize("n_width", [8, 16])
def test_hook_calculix_s8r_reaches_ko_reference(tmp_path, n_width):
    """Same model as the live hook test (E = 3300, nu = 0.3, t = 2, unit tip shear)."""
    wy = _hook_displacement(tmp_path, ccx_bin_or_skip(), n_width, poisson=0.3)
    assert_relative_error(
        wy,
        HOOK_WREF,
        tol=TOL,
        reference_name="Ko et al. 2017 hook reference 4.82482 (MITC9, N = 64)",
        kind="code",
        what=f"CalculiX S8R hook, {n_width} x {6 * n_width} elements",
    )


# --- Hyperbolic paraboloid -------------------------------------------------------------------
HP_E, HP_NU = 2.0e11, 0.3


def _hp_point(u: float, v: float):
    x, y = -0.5 + u, -0.5 + v
    return (x, y, y * y - x * x)


def _hp_free_edge_uz(tmp_path, ccx_bin: str, n: int, t_over_l: float, rho: float) -> float:
    centre = node_id(2 * n, n, n)  # x = +0.5, y = 0
    inp = tmp_path / f"hp_{n}_{int(1 / t_over_l)}.inp"
    write_s8r_surface_deck(
        inp,
        _hp_point,
        n,
        n,
        young=HP_E,
        poisson=HP_NU,
        thickness=t_over_l,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(x + 0.5) < 1e-9 else [],
        density=rho,
        gravity=(0.0, 0.0, -1.0),
        print_nodes=[centre],
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    return abs(float(read_dat_displacements(inp.with_suffix(".dat"))[centre][2]))


@pytest.mark.parametrize("n", [48, 64])
def test_hyperbolic_paraboloid_calculix_s8r_reaches_ko_reference(tmp_path, n):
    """t/L = 1/1000 only.

    Measured S8R sweep, w / w_ref over 12 / 24 / 48 / 64 elements per side: 0.814 / 0.953 / 0.989 /
    0.994 (t/L = 1/1000, monotone). At t/L = 1/10000 the same sweep is not monotone (0.28 / 0.77 /
    0.985 / 0.940, and 0.88 / 0.95 / 0.94 / 0.99 after rescaling E), so CalculiX S8R is not a
    converged arbiter there and that case has no second source.
    """
    uz = _hp_free_edge_uz(tmp_path, ccx_bin_or_skip(), n, 1 / 1000, 360.0)
    assert_relative_error(
        uz,
        2.8780e-4,
        tol=TOL,
        reference_name="Ko et al. 2017 hypar reference (MITC9, 72 x 144), t/L = 1/1000",
        kind="code",
        what=f"CalculiX S8R hypar t/L = 1/1000, {n} x {n} elements",
    )


# --- Thin twisted beam, Simo et al. 1989 -----------------------------------------------------
TB_LENGTH, TB_WIDTH, TB_TWIST = 12.0, 1.1, math.radians(90.0)
TB_E, TB_NU, TB_THICKNESS = 29.0e6, 0.22, 0.05
SIMO_THIN = {"In-plane": 1.390, "Out-of-plane": 0.3431}


def _twisted_point(u: float, v: float):
    """Same mapping as the live test: width along y at the root, rotated by 90 deg at the tip."""
    phi = TB_TWIST * u
    y = (v - 0.5) * TB_WIDTH
    return (u * TB_LENGTH, y * math.cos(phi), y * math.sin(phi))


@pytest.mark.parametrize("load_case", ["In-plane", "Out-of-plane"])
def test_thin_twisted_beam_calculix_s8r_matches_simo(tmp_path, load_case):
    n_width = 8
    n_length = 6 * n_width
    tip_centre = node_id(2 * n_length, n_width, n_length)
    load_dof = 3 if load_case == "In-plane" else 2  # 1-based: z for in-plane, y for out-of-plane
    inp = tmp_path / f"twisted_{load_case[0]}.inp"
    write_s8r_surface_deck(
        inp,
        _twisted_point,
        n_length,
        n_width,
        young=TB_E,
        poisson=TB_NU,
        thickness=TB_THICKNESS,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(x) < 1e-9 else [],
        nodal_loads={tip_centre: (load_dof, 1.0)},
        print_nodes=[tip_centre],
    )
    result = run_ccx(inp, ccx_bin_or_skip())
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))[tip_centre]
    assert_relative_error(
        abs(float(disp[load_dof - 1])),
        SIMO_THIN[load_case],
        tol=TOL,
        reference_name="Simo, Fox & Rifai 1989 thin twisted beam (b = 0.05 in)",
        kind="paper",
        what=f"CalculiX S8R thin twisted beam {load_case}",
    )
