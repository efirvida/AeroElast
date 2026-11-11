"""S-7 on the windIO v2 blade: shell tip torque vs BeamDyn K[5,5] (v2 counterpart of group 34).

Same load (tip ring Saint-Venant wall shear flow, 1000 N.m), same section-rotation estimator
(``section_twists_deg``), same reference (BeamDyn deck integrated to each z, ``beam_twist_profile_deg``),
same metric (mean shell/beam twist ratio over 0.3-0.9 span), same tolerance (``_RATIO_TOL = 0.05``,
the L3 ceiling) as ``tests/validation/blade/test_iea15mw_s7_torsion.py``.  Only the mesh differs:
``build_v2_blade(0.25, skin_reference="mid")`` instead of ``BladeMesh(...).generate(renumber="rcm")``.

Not repeated here: the v1 file's two machinery tests (reference integral exactness, applied moment
matches request) depend on neither mesher.  See docs/validation/_pending/62-v2_blade_s7_torsion/NOTES.md.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

from tests.support.assertions import assert_relative_error  # noqa: E402
from tests.support.paths import REPO_ROOT  # noqa: E402
from tests.support.v2_blade import build_v2_blade  # noqa: E402

sys.path.insert(0, str(REPO_ROOT / "tools"))

from run_s7_torsion import (  # noqa: E402
    beam_twist_profile_deg,
    section_twists_deg,
)

MOMENT_NM = 1000.0
ELEMENT_SIZE = 0.25
#: L3 ceiling, identical to the v1 test: band [0.95, 1.05] around exact agreement.
_RATIO_TOL = 0.05
SKIN_REFERENCE = "mid"


def run_torsion_case_v2(mesh, properties, torque: float) -> np.ndarray:
    """``tools.run_s7_torsion.run_torsion_case`` with the tip ring taken from the ``TipNodes`` set.

    The only change: v1 picks the ring as ``z >= z.max() - 1e-6``.  The v2 tip section is not planar in z
    (30 ring nodes spread over 1.0 cm of z because of pre-bend/sweep and the mid-surface inset), so that
    mask selects a single node and ``ring_section`` raises.  Everything else is copied verbatim.
    """
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.elements import ElementFamily
    from aeroelast.solvers.bem.force_projection import section_torque_wall_forces
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    cfg = {
        "solver": {},
        "elements": {
            "element_family": ElementFamily.SHELL,
            "span_direction": (0.0, 0.0, 1.0),
            "properties": properties,
        },
    }
    solver = StaticLinearSolver(mesh, cfg)
    dpn = solver.domain.dofs_per_node
    idx_of = mesh.node_id_to_index
    root_ids = sorted(mesh.get_node_set("RootNodes").nodes.keys())
    root_dofs = sorted(idx_of[nid] * dpn + d for nid in root_ids for d in range(dpn))
    solver.add_dirichlet_conditions([DirichletCondition(root_dofs, 0.0)])

    tip = np.asarray(
        sorted(idx_of[nid] for nid in mesh.get_node_set("TipNodes").nodes.keys()), dtype=np.intp
    )
    tip_nodes, f_tip = section_torque_wall_forces(
        mesh, tip, (0.0, 0.0, 1.0), properties, torque, coords=mesh.coords_array
    )
    dofs, vals = [], []
    for t, f in zip(tip_nodes, f_tip, strict=True):
        for d in range(3):
            dofs.append(int(t) * dpn + d)
            vals.append(float(f[d]))
    solver.add_nodal_loads([NodalLoad(dofs, vals)])
    return solver.solve()


@pytest.mark.slow
def test_v2_s7_global_torsion_ratio_within_band():
    mesh, properties = build_v2_blade(ELEMENT_SIZE, skin_reference=SKIN_REFERENCE)
    u = run_torsion_case_v2(mesh, properties, MOMENT_NM)
    u_arr = np.asarray(u, dtype=np.float64)

    zs, theta_s = section_twists_deg(mesh, u_arr)
    theta_beam = beam_twist_profile_deg(MOMENT_NM, zs)

    mid = (zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)
    ratio_mean = float(np.nanmean(theta_s[mid] / theta_beam[mid]))

    assert_relative_error(
        ratio_mean,
        1.0,
        tol=_RATIO_TOL,
        kind="code",
        reference_name=(
            "BeamDyn blade deck of the official IEA 15 MW reference turbine "
            "(IEA-15-240-RWT_BeamDyn_blade.dat, K[5,5] per station), integrated "
            "as a Timoshenko torsion beam"
        ),
        what="shell/beam tip twist ratio over the 0.3-0.9 span band (windIO v2 mesh)",
    )
