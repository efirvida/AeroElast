"""Pre-twisted beam under a large tip load: the production nonlinear static solver against CalculiX.

Simo, Fox & Rifai (1989) thin twisted beam: L = 12 in, b = 1.1 in, t = 0.05 in, 90 deg of twist
from root to tip (mapping of ``test_ko2017_calculix_arbiters.py``), E = 29e6, nu = 0.22.  The
published tip loads (1.0) are linear; loads of 5 and 10 move the tip by 15-65 % of the length, so
the section rotates about more than one axis while it bends and twists, which the planar elastica
cannot exercise.  The arbiter is CalculiX 2.20 S8R with ``NLGEOM``
(20 increments), same geometry, material, clamp and dead tip load.

Measured 2026-10-14 (tip displacement vector, relative norm of the difference to CalculiX):

    case                    mesh   Newton (default)   updated-Lagrangian (scheme ul)
    in-plane, P = 10        24x4        0.52 %             23.30 %
    in-plane, P = 10        48x8        0.19 %             23.28 %
    out-of-plane, P = 5     24x4        2.36 %              2.71 %
    out-of-plane, P = 5     48x8        0.81 %              1.98 %

Newton converges with the mesh; the legacy scheme's 23 % does not move with it, so it is a
property of the scheme and not of the discretisation (pinned in
``tests/software/test_static_nonlinear_scheme.py``).  Out-of-plane loads bend about the weak
axis of the twisted section and are the more mesh-sensitive case in both codes.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.bc import DirichletCondition, NodalLoad
from aeroelast.core.material import IsotropicMaterial
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_residual_below
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck
from tests.support.strip_cantilever import build_cantilever_mesh

LENGTH, WIDTH, TWIST = 12.0, 1.1, math.radians(90.0)
YOUNG, POISSON, THICKNESS = 29.0e6, 0.22, 0.05
CCX_NAME = "CalculiX 2.20 S8R, NLGEOM, dead CLOAD, 20 increments"
# 0-based global dof of the tip load: z is the in-plane case, y the out-of-plane one.
LOAD_DOF = {"in-plane": 2, "out-of-plane": 1}


def _point(u: float, v: float):
    phi = TWIST * u
    y = (v - 0.5) * WIDTH
    return (u * LENGTH, y * math.cos(phi), y * math.sin(phi))


def _ccx_tip(tmp_path, ccx_bin: str, case: str, load: float, nl: int, nw: int) -> np.ndarray:
    tip = node_id(2 * nl, nw, nl)
    inp = tmp_path / f"twisted_{case}_{nl}x{nw}.inp"
    write_s8r_surface_deck(
        inp,
        _point,
        nl,
        nw,
        young=YOUNG,
        poisson=POISSON,
        thickness=THICKNESS,
        boundary=lambda x, y, z: [1, 2, 3, 4, 5, 6] if abs(x) < 1e-9 else [],
        nodal_loads={tip: (LOAD_DOF[case] + 1, load)},
        print_nodes=[tip],
        nlgeom=20,
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    return read_dat_displacements(inp.with_suffix(".dat"))[tip]


def _ours_tip(case: str, load: float, nl: int, nw: int, scheme: str | None) -> np.ndarray:
    mesh = build_cantilever_mesh(L=LENGTH, b=WIDTH, nx=nl, ny=nw)
    centre = min(
        (n for n in mesh.nodes if abs(n.x - LENGTH) < 1e-9), key=lambda n: abs(n.y - WIDTH / 2)
    )
    for n in mesh.nodes:  # the flat strip becomes the twisted one; node sets keep their nodes
        phi = TWIST * n.x / LENGTH
        y = n.y - WIDTH / 2
        n.coords[0], n.coords[1], n.coords[2] = (
            float(n.x),
            float(y * math.cos(phi)),
            float(y * math.sin(phi)),
        )
    solver_cfg: dict = {"continuation_steps": 8, "continuation_max_steps": 64}
    if scheme is not None:
        solver_cfg["scheme"] = scheme
    if scheme == "ul":
        solver_cfg = {"scheme": "ul", "continuation_steps": 64, "continuation_max_steps": 512}
    cfg = {
        "solver": solver_cfg,
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {
                "plate": ShellProperty(
                    material=IsotropicMaterial(name="Simo", E=YOUNG, nu=POISSON, rho=1.0),
                    thickness=THICKNESS,
                )
            },
            "span_direction": (1.0, 0.0, 0.0),
        },
    }
    solver = StaticNonlinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    clamped = []
    for node in mesh.get_node_set("clamped").nodes.values():
        i0 = mesh.node_id_to_index[node.id] * dpn
        clamped.extend(range(i0, i0 + dpn))
    solver.add_dirichlet_conditions([DirichletCondition(clamped, 0.0)])
    i0 = mesh.node_id_to_index[centre.id] * dpn
    force = [0.0] * dpn
    force[LOAD_DOF[case]] = load
    solver.add_nodal_loads([NodalLoad(list(range(i0, i0 + dpn)), force)])
    return np.asarray(solver.solve(), dtype=np.float64)[i0 : i0 + 3]


@pytest.mark.parametrize(
    "case,load,nl,nw,tol",
    [
        ("in-plane", 10.0, 24, 4, 0.01),
        ("in-plane", 10.0, 48, 8, 0.01),
        ("out-of-plane", 5.0, 48, 8, 0.015),
    ],
)
def test_twisted_beam_default_scheme_matches_calculix_nlgeom(tmp_path, case, load, nl, nw, tol):
    ccx = _ccx_tip(tmp_path, ccx_bin_or_skip(), case, load, nl, nw)
    ours = _ours_tip(case, load, nl, nw, scheme=None)
    assert_residual_below(
        float(np.linalg.norm(ours - ccx)) / float(np.linalg.norm(ccx)),
        tol=tol,
        kind="code",
        reference_name=CCX_NAME,
        what=f"twisted beam {case} P={load:g} tip displacement vector, {nl}x{nw}",
    )
