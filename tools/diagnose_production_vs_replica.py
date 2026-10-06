"""WU-3a prototype: is the production static solver an arbiter of the scipy replica?

Every test that calls ``spsolve`` assembles ``K`` with the production Rust
assembler and then solves with scipy. None of them ever drives a production
solver, so the production solve path is currently unasserted.

This probe runs the production solver on the SAME problem as an existing
replica test -- same mesh, same ``model_config``, same BCs, same loads -- and
reports the gap. The fixture is loaded from the test file by path rather than
copied, so the two sides cannot drift.

The number it prints is the acceptance gap that WU-3b's shared helper would
assert on; the tolerance is a maintainer decision.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE = REPO_ROOT / "tests" / "validation" / "parity" / "test_shell_stress_ccx_parity.py"


def _load_fixture(path: Path):
    """Load the reference test module from its file path (no package, no path mutation)."""
    spec = importlib.util.spec_from_file_location("_arbiter_fixture", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load the reference fixture at {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    from aeroelast.core.bc import DirichletCondition, NodalLoad
    from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver

    ref = _load_fixture(FIXTURE)
    mesh = ref._build_plate()
    cfg = ref._model_cfg()
    n_dof = 6 * len(list(mesh.nodes))

    # ---- replica side: the test's own solve, production K + scipy ----------
    _domain, u_replica = ref._solve_aeroelast(mesh)

    # ---- production side: same problem, driven by the production solver ----
    solver = StaticLinearSolver(mesh, cfg)

    index = mesh.node_id_to_index
    root = {index[nid] for nid in mesh.get_node_set("clamped").node_ids}
    fixed_dofs = sorted({6 * i + d for i in root for d in range(6)})
    solver.add_dirichlet_conditions([DirichletCondition(fixed_dofs, 0.0)])

    free_nodes = sorted(
        (nd for nd in mesh.nodes if np.isclose(nd.x, ref.L, atol=1e-12)),
        key=lambda nd: nd.y,
    )
    per_node = ref.FORCE / len(free_nodes)
    load_dofs = [6 * index[nd.id] + 2 for nd in free_nodes]
    solver.add_nodal_loads([NodalLoad(load_dofs, [-per_node] * len(load_dofs))])

    u_prod = np.asarray(solver.solve(), dtype=np.float64)

    # ---- the gap ----------------------------------------------------------
    if u_prod.shape != (n_dof,) or u_replica.shape != (n_dof,):
        raise RuntimeError(
            f"shape mismatch: production {u_prod.shape}, replica {u_replica.shape}, expected ({n_dof},)"
        )
    du = u_prod - u_replica
    scale = float(np.max(np.abs(u_replica)))
    rel_inf = float(np.max(np.abs(du)) / scale) if scale else float("nan")

    tip = 6 * index[free_nodes[len(free_nodes) // 2].id] + 2
    rel_tip = abs(float(u_prod[tip]) - float(u_replica[tip])) / abs(float(u_replica[tip]))

    print("=" * 78)
    print("WU-3a  production static solver vs scipy replica  (same K, same f, same BCs)")
    print("=" * 78)
    print(f"  case                     : shell_stress_ccx_parity plate (8x2 quads, {n_dof} DOF)")
    print("  replica  source          : ref._solve_aeroelast -> spsolve on the Rust K")
    print("  production source        : StaticLinearSolver -> _aeroelast.linear_static_solve_coo")
    print(f"  max |u|                  : {scale:.6e}")
    print(f"  max |u_prod - u_replica| : {float(np.max(np.abs(du))):.6e}")
    print(f"  relative (inf, vs max)   : {rel_inf:.3e}")
    print(f"  tip w  replica           : {float(u_replica[tip]):.10e}")
    print(f"  tip w  production        : {float(u_prod[tip]):.10e}")
    print(f"  tip w  relative gap      : {rel_tip:.3e}")
    print(f"  the test's tolerance vs CCX : {ref.TOL_CCX}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
