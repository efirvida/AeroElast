"""D-Tube cantilever under a dead body load (self weight), AeroElast and CalculiX S8R.

The pressure tests (``test_dtube_camarena_bending.py``, ``test_dtube_calculix.py``) load the arc by
a surface traction. Here the same tube carries a vertical body acceleration ``a`` chosen so that the
resultant per unit span equals the pressure one, ``rho a A = p D``: a uniform dead load through the
section centroid, so the closed form is the same prismatic cantilever ``w L^4 / (8 E I_c)``. This
exercises the body-force and mass path of the shell solvers, which the blade gravity cases depend on.
A body load is not follower: it is the dead load a nonlinear step must also keep fixed.

CalculiX point loads (``*CLOAD``) spread over the shell nodes lose 25-44 % of the resultant on this
deck and are not used; ``*DLOAD ... GRAV`` is.

Measured 2026-10-14: AeroElast -16.8 / -4.70 / -1.02 / -0.05 % at 50 / 100 / 200 / 400 layers,
CalculiX S8R +0.11 % at both 50 and 100 layers.
"""

from __future__ import annotations

import pytest

from aeroelast.core.bc import BodyForce, DirichletCondition
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_relative_error
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck
from tests.validation.benchmarks.test_dtube_camarena_bending import (
    A_SEC,
    DELTA_REF,
    E,
    L,
    N_ARC,
    N_WEB,
    NU,
    R,
    T,
    W_RES,
    _build_dtube_mesh,
)

RHO = 7850.0
ACCEL = W_RES / (RHO * A_SEC)  # m/s^2, +y: rho * a * A = p * D


def _aeroelast_tip(n_z: int) -> float:
    mesh, arc_nodes, _ = _build_dtube_mesh(n_z=n_z)
    prop = ShellProperty(material=IsotropicMaterial(name="Steel", E=E, nu=NU, rho=RHO), thickness=T)
    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"tube": prop},
            "span_direction": (0.0, 0.0, 1.0),
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamped = sorted(
        mesh.node_id_to_index[n.id] * dpn + d
        for n in mesh.get_node_set("clamped").nodes.values()
        for d in range(dpn)
    )
    solver.add_dirichlet_conditions([DirichletCondition(clamped, 0.0)])
    solver.add_body_forces([BodyForce([0.0, ACCEL, 0.0])])
    u = solver.solve()
    tip = [n for n in arc_nodes[n_z] if n is not None]
    return float(sum(u[mesh.node_id_to_index[n.id] * dpn + 1] for n in tip) / len(tip))


def _calculix_tip(tmp_path, ccx_bin: str, n_z: int) -> float:
    n_arc, n_web = N_ARC, N_WEB
    n_u = n_arc + n_web

    def point(u: float, v: float):
        import math

        s = u * n_u
        if s <= n_arc:
            th = math.pi * s / n_arc
            return (R * math.cos(th), -R * math.sin(th), v * L)
        return (-R + 2.0 * R * (s - n_arc) / n_web, 0.0, v * L)

    tip_ids = [node_id(i, 2 * n_z, n_u, periodic=True) for i in range(0, 2 * n_arc + 1)]
    inp = tmp_path / f"dtube_gravity_{n_z}.inp"
    write_s8r_surface_deck(
        inp,
        point,
        n_u,
        n_z,
        young=E,
        poisson=NU,
        thickness=T,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(z) < 1e-9 else [],
        periodic_u=True,
        density=RHO,
        gravity=(0.0, ACCEL, 0.0),
        print_nodes=tip_ids,
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))
    return float(sum(disp[n][1] for n in tip_ids) / len(tip_ids))


@pytest.fixture(scope="module")
def aeroelast_tips():
    return {n_z: _aeroelast_tip(n_z) for n_z in (100, 200, 400)}


@pytest.fixture(scope="module")
def calculix_tips(tmp_path_factory):
    ccx = ccx_bin_or_skip()
    root = tmp_path_factory.mktemp("dtube_gravity")
    return {n_z: _calculix_tip(root, ccx, n_z) for n_z in (50, 100)}


@pytest.mark.slow
def test_dtube_gravity_tip_deflection_aeroelast_against_beam_theory(aeroelast_tips):
    assert_relative_error(
        aeroelast_tips[400],
        DELTA_REF,
        tol=0.01,
        kind="analytical",
        reference_name="the prismatic cantilever w L^4 / (8 E I_c) about the section centroid",
        what="D-Tube tip deflection under self weight, AeroElast MITC4, 400 layers",
    )


@pytest.mark.slow
def test_dtube_gravity_aeroelast_converges_with_span_refinement(aeroelast_tips):
    errors = {n_z: abs(tip - DELTA_REF) / DELTA_REF for n_z, tip in aeroelast_tips.items()}
    assert errors[100] > errors[200] > errors[400], errors
    assert_relative_error(
        aeroelast_tips[200],
        DELTA_REF,
        tol=0.05,
        kind="analytical",
        reference_name="the prismatic cantilever w L^4 / (8 E I_c) about the section centroid",
        what="D-Tube tip deflection under self weight, AeroElast MITC4, 200 layers",
    )


@pytest.mark.slow
def test_dtube_gravity_tip_deflection_calculix_s8r_against_beam_theory(calculix_tips):
    assert_relative_error(
        calculix_tips[100],
        DELTA_REF,
        tol=0.05,
        kind="analytical",
        reference_name="the prismatic cantilever w L^4 / (8 E I_c) about the section centroid",
        what="D-Tube tip deflection under self weight, CalculiX S8R, 100 layers",
    )


@pytest.mark.slow
def test_dtube_gravity_calculix_is_span_converged(calculix_tips):
    assert_relative_error(
        calculix_tips[50],
        calculix_tips[100],
        tol=0.01,
        kind="self",
        reference_name="CalculiX S8R at 100 layers",
        what="D-Tube tip deflection under self weight, CalculiX S8R, 50 layers",
    )


@pytest.mark.slow
def test_dtube_gravity_aeroelast_and_calculix_agree_when_converged(aeroelast_tips, calculix_tips):
    assert_relative_error(
        aeroelast_tips[400],
        calculix_tips[100],
        tol=0.01,
        kind="code",
        reference_name="CalculiX 2.20 S8R at 100 layers",
        what="D-Tube tip deflection under self weight, AeroElast 400 layers",
    )
