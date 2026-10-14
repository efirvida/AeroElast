"""Large-deflection cantilever shell strip against CalculiX S8R (geometrically nonlinear).

The strip of ``test_static_nonlinear_ul_elastica.py`` (L = 1 m, b = 0.1 m, h = 1e-3, steel, dead
tip load alpha = P L^2 / EI of 1 and 2) is a plate, not the beam the elastica describes: at
b/L = 0.1 and h/L = 1e-3 it is about 4-5 % stiffer than the exact elastica (anticlastic and
Ashwell effects, which vanish as b -> 0: measured -5.0, -1.2, -0.36 % for b = 0.1, 0.05, 0.025).
The elastica is therefore the wrong arbiter for this strip, and CalculiX S8R with ``NLGEOM`` is
the right one: a different code and a different (quadratic) element, same geometry, material,
clamp and dead load.

Measured 2026-10-14 (CalculiX 2.20, 20 increments; 10 and 40 increments agree to 4 digits):

    alpha  mesh (u x v)   w/L CalculiX   w/L ours (newton, ny=4)   u/L CalculiX   u/L ours
    4/7    32 x 2           0.17731          0.17712                 -0.01916       -0.01909
    4/7    64 x 2           0.17778          0.17753                 -0.01925       -0.01918
    4/7    32 x 2           0.17731          0.17712                 -0.01916       -0.01909
    4/7    64 x 2           0.17778          0.17753                 -0.01925       -0.01918
      1    16 x 2           0.28664             -                    -0.05111          -
      1    32 x 2           0.28811          0.28745                 -0.05154       -0.05121
      1    64 x 2           0.28880          0.28806                 -0.05175       -0.05144
      2    32 x 2           0.47329          0.47252                 -0.14712       -0.14630
      2    64 x 2           0.47424          0.47300                 -0.14761       -0.14665

Both codes move monotonically with the span mesh and are still drifting 0.1-0.3 % per doubling, so
the paired bound is 1 % on the tip deflection and 2 % on the shortening, not a tighter number the
sequences cannot yet support.  The default updated-Lagrangian scheme is not compared here: it sits
away from CalculiX by +2.0, +4.4 and +9.6 % at alpha = 4/7, 1 and 2 (see the module docstring of
``test_static_nonlinear_newton_elastica.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.core.bc import DirichletCondition
from aeroelast.core.properties import ShellProperty
from aeroelast.elements import ElementFamily
from aeroelast.solvers.elasticity.static_nonlinear import StaticNonlinearSolver

from tests.conftest import ccx_bin_or_skip
from tests.support.assertions import assert_relative_error
from tests.support.ccx_io import fail_ccx, run_ccx
from tests.support.ccx_surface import node_id, read_dat_displacements, write_s8r_surface_deck
from tests.validation.element.test_static_nonlinear_ul_elastica import (
    E,
    H,
    L,
    STEEL,
    apply_edge_load,
    build_cantilever_mesh,
    get_center_node,
)

B, NU = 0.1, 0.3
EI = E * B * H**3 / 12.0
CCX_NAME = "CalculiX 2.20 S8R, NLGEOM, dead CLOAD, 20 increments"


def _ccx_tip(tmp_path, ccx_bin: str, alpha: float, n: int) -> tuple[float, float]:
    """Tip (w/L, u/L) of the strip in CalculiX: ``n x 2`` S8R elements."""
    nv = 2
    load = alpha * EI / L**2
    weights = np.array(
        [1.0 if j in (0, 2 * nv) else (4.0 if j % 2 else 2.0) for j in range(2 * nv + 1)]
    )
    weights /= weights.sum()  # consistent nodal loads of a quadratic edge
    tip = node_id(2 * n, nv, n)
    inp = tmp_path / f"strip_{alpha:g}_{n}.inp"
    write_s8r_surface_deck(
        inp,
        lambda u, v: (L * u, B * v, 0.0),
        n,
        nv,
        young=E,
        poisson=NU,
        thickness=H,
        boundary=lambda x, y, z: range(1, 7) if abs(x) < 1e-12 else (),
        nodal_loads={node_id(2 * n, j, n): (3, load * weights[j]) for j in range(2 * nv + 1)},
        print_nodes=[tip],
        nlgeom=20,
    )
    result = run_ccx(inp, ccx_bin)
    if result.returncode != 0:
        fail_ccx(result, inp)
    disp = read_dat_displacements(inp.with_suffix(".dat"))[tip]
    return abs(float(disp[2])) / L, float(disp[0]) / L


def _ours_tip(alpha: float, n: int) -> tuple[float, float]:
    """Tip (w/L, u/L) with ``scheme: newton`` on an ``n x 4`` MITC4 mesh."""
    mesh = build_cantilever_mesh(L=L, b=B, nx=n, ny=4)
    cfg = {
        "solver": {"scheme": "newton", "continuation_steps": 8, "continuation_max_steps": 64},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "properties": {"plate": ShellProperty(material=STEEL, thickness=H)},
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
    apply_edge_load(mesh, solver, "free_edge", "z", alpha * EI / L**2)
    u = np.asarray(solver.solve(), dtype=np.float64)
    idx = mesh.node_id_to_index[get_center_node(mesh, L, B).id]
    return abs(u[idx * dpn + 2]) / L, u[idx * dpn] / L


@pytest.mark.parametrize("n", [32, 64])
# 4/7 is the load of ``TestNonlinearStatic`` in ``test_shell_static_and_modal.py`` (P = 1 N, EI = 1.75).
@pytest.mark.parametrize("alpha", [4.0 / 7.0, 1.0, 2.0])
def test_newton_strip_matches_calculix_nlgeom(tmp_path, alpha, n):
    w_ccx, u_ccx = _ccx_tip(tmp_path, ccx_bin_or_skip(), alpha, n)
    w, u = _ours_tip(alpha, n)
    assert_relative_error(
        w,
        w_ccx,
        tol=0.01,
        kind="code",
        reference_name=CCX_NAME,
        what=f"alpha={alpha:.6g} tip w/L, {n} span elements",
    )
    assert_relative_error(
        u,
        u_ccx,
        tol=0.02,
        kind="code",
        reference_name=CCX_NAME,
        what=f"alpha={alpha:.6g} tip u/L, {n} span elements",
    )
