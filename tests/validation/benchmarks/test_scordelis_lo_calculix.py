"""Scordelis-Lo roof against an independent code: CalculiX S8R, refined.

The reference ``w_ref = 3.0240e-1`` (MacNeal & Harder 1985) is itself numerical, and the published
cells the smoothed-MITC3 tests use (Lee & Lee 2019) are finite-element results. CalculiX S8R is a
third, independent code, so its refined value says whether ``w_ref`` is the converged answer of
this model before our element is judged against it.
"""

from __future__ import annotations

import math

import pytest

from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_relative_error
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck

R = 25.0
L_HALF = 25.0
ANGLE = math.radians(40.0)
THICKNESS = 0.25
WREF = 3.0240e-1
TOL = 0.05


def _point(u: float, v: float):
    th = u * ANGLE
    return (R * math.cos(th), R * math.sin(th), v * L_HALF)


def _boundary(x: float, y: float, z: float):
    """Quarter model: symmetry in length (z = 0), rigid diaphragm (z = L/2), symmetry in y = 0.

    Only the bending rotation of each symmetry plane is constrained. A symmetry plane also asks
    for the drilling rotation to vanish, but CalculiX S8R expands the shell and does not carry a
    drilling dof: constraining both rotations at z = 0 stiffens the roof by a factor 7.7
    (w/w_ref 0.13 instead of 1.00 at 8x8, measured), while either rotation alone gives
    0.9985 and 0.9985 -- the answer is insensitive to which one is kept.
    """
    tol = 1e-6
    dofs: list[int] = []
    if abs(z) < tol:  # w = 0, bending rotation about x
        dofs += [3, 4]
    if abs(z - L_HALF) < tol:  # diaphragm: u, v
        dofs += [1, 2]
    if abs(y) < tol and x > 0.0:  # v = 0, rotation about z
        dofs += [2, 6]
    return dofs


def _normalized_tip(tmp_path, ccx_bin: str, n: int) -> float:
    inp = tmp_path / f"scordelis_s8r_{n}.inp"
    write_s8r_surface_deck(
        inp,
        _point,
        n,
        n,
        young=4.32e8,
        poisson=0.0,
        thickness=THICKNESS,
        boundary=_boundary,
        density=360.0,
        gravity=(-1.0, 0.0, 0.0),
        print_nodes=[node_id(2 * n, 0, n)],
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    tip = node_id(2 * n, 0, n)  # theta = 40 deg, z = 0
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    return abs(float(disp[tip][0])) / WREF


@pytest.mark.parametrize("n", [8, 16, 32])
def test_scordelis_lo_calculix_s8r_reaches_macneal_harder(tmp_path, n):
    norm = _normalized_tip(tmp_path, ccx_bin_or_skip(), n)
    assert_relative_error(
        norm,
        1.0,
        tol=TOL,
        reference_name="MacNeal & Harder 1985 Scordelis-Lo w_ref = 0.3024",
        kind="code",
        what=f"CalculiX S8R Scordelis-Lo {n}x{n}",
    )
