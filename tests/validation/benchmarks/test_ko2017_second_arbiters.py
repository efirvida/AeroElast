"""Second, independent arbiters for Ko et al. (2017) cases whose reference is one paper's own cell.

Each test here compares our MITC4+ shell with a source other than the Ko table the live test in
``test_ko2017_performance.py`` uses:

* Pinched cylinder: the MacNeal-Harder closed-form (series) tip displacement ``w_ref = 1.8248e-5``,
  approached by mesh refinement. At N = 16 (the live test's mesh) the shell is below it, which is a
  discretisation effect that this table measures instead of assuming.
* Thin twisted beam: the reference of Simo, Fox & Rifai (1989) for the thin case ``b = 0.05 in``
  (in-plane 1.390 in, out-of-plane 0.3431 in), as tabulated in the Abaqus Benchmarks Manual.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("_aeroelast", reason="Rust backend not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.validation.benchmarks.test_ko2017_performance import (  # noqa: E402
    DOF,
    MAT_CYL,
    MAT_TB,
    _build_cylindrical_patch,
    _build_twisted_beam_mesh,
    _Case,
    _find_node_by_xyz,
    _pinched_cylinder_fixed,
    _pinched_cylinder_load,
    _point_load_at_coords,
    _run_case,
    _twisted_beam_fixed,
)

TOL = 0.05  # suite rule

# --- Pinched cylinder -------------------------------------------------------------------------
CYL_R = 300.0
CYL_L_HALF = 300.0
CYL_T = 3.0
CYL_WREF = 1.8248e-5  # MacNeal & Harder 1985


def _pinched_cylinder_norm(n: int) -> float:
    case = _Case(
        name=f"pinched_cylinder_reg_N{n}",
        element="MITC4",
        build_mesh=lambda: _build_cylindrical_patch(
            radius=CYL_R, length=CYL_L_HALF, angle_deg=90.0, nx=n, ny=n, triangular=False
        ),
        material=MAT_CYL,
        thickness=CYL_T,
        load_vector=lambda mesh, m: _pinched_cylinder_load(mesh, m, radius=CYL_R),
        fixed_dofs=lambda mesh, m: _pinched_cylinder_fixed(mesh, m, length=CYL_L_HALF),
        measure=lambda mesh, m, u: float(
            u[m[_find_node_by_xyz(mesh, CYL_R, 0.0, 0.0).id] * DOF + 0]
        ),
        expected_normalized=1.0,
        wref=CYL_WREF,
    )
    return _run_case(case)


@pytest.mark.parametrize("n", [32])
def test_pinched_cylinder_refined_mesh_reaches_macneal_harder(n):
    """N = 32 is the finest mesh the dense global assembly of this module solves in memory."""
    norm = _pinched_cylinder_norm(n)
    assert_relative_error(
        norm,
        1.0,
        tol=TOL,
        reference_name="MacNeal & Harder 1985 pinched cylinder w_ref = 1.8248e-5",
        kind="analytical",
        what=f"pinched cylinder regular mesh N={n}",
    )


def test_pinched_cylinder_converges_monotonically_towards_the_reference():
    """The N = 16 gap to the reference is discretisation: the series approaches 1 from below."""
    norms = [_pinched_cylinder_norm(n) for n in (8, 16, 32)]
    print(f"pinched cylinder w/w_ref at N=8/16/32: {norms}")
    assert all(a < b for a, b in zip(norms, norms[1:])), norms
    assert norms[-1] < 1.0, norms


# --- Thin twisted beam, Simo 1989 --------------------------------------------------------------
TB_LENGTH = 12.0
TB_WIDTH = 1.1
TB_THIN = 0.05  # in, Simo b = 0.05
SIMO_THIN = {"In-plane": 1.390, "Out-of-plane": 0.3431}


@pytest.mark.parametrize("load_case", ["In-plane", "Out-of-plane"])
def test_thin_twisted_beam_against_simo_1989(load_case):
    n_width = 8
    load_dof = 2 if load_case == "In-plane" else 1
    case = _Case(
        name=f"twisted_beam_thin_{load_case}_b0.05",
        element="MITC4",
        build_mesh=lambda: _build_twisted_beam_mesh(
            length=TB_LENGTH,
            width=TB_WIDTH,
            twist_deg=90.0,
            nx=6 * n_width,
            ny=n_width,
            triangular=False,
        ),
        material=MAT_TB,
        thickness=TB_THIN,
        load_vector=lambda mesh, m: _point_load_at_coords(
            mesh, m, x=TB_LENGTH, y=0.0, z=0.0, load=1.0, dof=load_dof
        ),
        fixed_dofs=_twisted_beam_fixed,
        measure=lambda mesh, m, u: float(
            u[m[_find_node_by_xyz(mesh, TB_LENGTH, 0.0, 0.0).id] * DOF + load_dof]
        ),
        expected_normalized=1.0,
        wref=SIMO_THIN[load_case],
    )
    norm = _run_case(case)
    assert_relative_error(
        norm,
        1.0,
        tol=TOL,
        reference_name="Simo, Fox & Rifai 1989 thin twisted beam (b = 0.05 in)",
        kind="paper",
        what=case.name,
    )
