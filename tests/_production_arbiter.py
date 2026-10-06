"""Shared arbiter: drive the PRODUCTION solver on a replica test's own problem.

Almost every validation test assembles ``K`` with the production Rust assembler
and then solves with scipy (``spsolve`` / ``splu``). The reference comparison
proves the *model*; nothing proves the production *solver*. This module closes
that: given the mesh, the ``model_config`` and the boundary conditions a replica
test already builds, it runs the production solver on the same problem and
returns the gap.

Contract, deliberately split so the validation store still sees the comparison:

* this module solves and MEASURES -- it returns numbers;
* the test body writes the ``assert_relative_error(...)`` call on the physical
  scalar, because the store's extractor only reads assertions written in the
  test body (see ``tests/support/assertions.py``).

The ``1e-6`` default bound comes from the first measurement (WU-3a, a 162 DOF
plate: 1.85e-10 relative between PETSc KSP CG+GAMQ and scipy), which leaves
~5000x of margin while still being five orders tighter than the 5% the same test
uses against CalculiX. It is a per-case bound: measure the gap on a new case
before adopting it, because a well-conditioned toy does not predict a blade mesh.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pytest

#: Default acceptance bound for the production-vs-replica gap. See the module docstring.
DEFAULT_RTOL = 1e-6


def fixed_dofs_from_node_set(mesh, node_set: str, dofs_per_node: int = 6) -> list[int]:
    """All DOFs of the nodes in ``node_set``, sorted (the production BC convention)."""
    index = mesh.node_id_to_index
    nodes = {index[nid] for nid in mesh.get_node_set(node_set).node_ids}
    return sorted({dofs_per_node * i + d for i in nodes for d in range(dofs_per_node)})


def nodes_at_coordinate(mesh, value: float, axis: int = 0, atol: float = 1e-12) -> list[int]:
    """Node indices (``node_id_to_index``, as the tests use) at ``coordinate[axis] == value``.

    Ordered by the two remaining coordinates so a load spread over the returned nodes is
    deterministic. Indices are the ones the DOF numbering uses, not row positions.
    """
    index = mesh.node_id_to_index
    picked = []
    for node in mesh.nodes:
        coords = (node.x, node.y, node.z)
        if abs(float(coords[axis]) - value) <= atol:
            picked.append((index[node.id], coords))
    others = [a for a in (0, 1, 2) if a != axis]
    picked.sort(key=lambda item: (item[1][others[0]], item[1][others[1]]))
    return [i for i, _ in picked]


def nodal_loads(node_indices: Sequence[int], value: float, dof: int, dofs_per_node: int = 6):
    """A single ``NodalLoad`` spreading ``value`` equally over ``dof`` of every node."""
    from aeroelast.core.bc import NodalLoad

    per_node = value / len(node_indices)
    dofs = [dofs_per_node * i + dof for i in node_indices]
    return NodalLoad(dofs, [per_node] * len(dofs))


def nodal_loads_from_vector(force: np.ndarray, dofs_per_node: int = 6):
    """A single ``NodalLoad`` carrying every non-zero entry of a full force vector.

    Lets an arbiter test reuse the force vector its replica already builds instead of rebuilding
    the same load through a different API, which would be a second chance to get it wrong.
    """
    from aeroelast.core.bc import NodalLoad

    dofs = [i for i in range(force.size) if force[i] != 0.0]
    return NodalLoad(dofs, [float(force[i]) for i in dofs])


def solve_static(
    mesh,
    model_cfg: dict,
    *,
    fixed_dofs: Iterable[int],
    nodal_loads_list: Sequence | None = None,
    body_forces: Sequence | None = None,
) -> np.ndarray:
    """Run the production static solve (PETSc KSP via ``_aeroelast``) and return ``u``.

    Skips -- rather than fails -- when the Rust extension is unavailable, mirroring
    ``tests/_ccx_io.ccx_bin_or_skip``: the suite must still run without PETSc.
    """
    try:
        from aeroelast.core.bc import DirichletCondition
        from aeroelast.solvers.elasticity.static_linear import StaticLinearSolver
    except ImportError as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"production static solver unavailable: {exc}")

    solver = StaticLinearSolver(mesh, model_cfg)
    solver.add_dirichlet_conditions([DirichletCondition(sorted(fixed_dofs), 0.0)])
    if nodal_loads_list:
        solver.add_nodal_loads(list(nodal_loads_list))
    if body_forces:
        solver.add_body_forces(list(body_forces))

    try:
        u = np.asarray(solver.solve(), dtype=np.float64)
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"production static solve did not run: {exc}")
    return u


def displacement_gap(u_prod: np.ndarray, u_replica: np.ndarray) -> dict:
    """Absolute and relative gap between two displacement fields, for the record."""
    if u_prod.shape != u_replica.shape:
        raise ValueError(f"shape mismatch: production {u_prod.shape} vs replica {u_replica.shape}")
    du = np.abs(u_prod - u_replica)
    scale = float(np.max(np.abs(u_replica)))
    return {
        "max_abs": float(np.max(du)),
        "max_ref": scale,
        "rel_inf": float(np.max(du) / scale) if scale else float("nan"),
        "dofs": int(u_prod.size),
    }


def require_gap(gap: dict, *, rtol: float = DEFAULT_RTOL, label: str = "") -> None:
    """Assert the field-level gap stays inside ``rtol``, reporting the measured numbers."""
    print(
        f"  production vs replica{(' [' + label + ']') if label else ''}: "
        f"max|du|={gap['max_abs']:.3e} rel={gap['rel_inf']:.3e} bound={rtol:.1e} "
        f"({gap['dofs']} DOF)"
    )
    assert gap["rel_inf"] < rtol, (
        f"the production solver and the scipy replica disagree by {gap['rel_inf']:.3e} "
        f"(bound {rtol:.1e}) on {label or 'this case'}: max|du|={gap['max_abs']:.6e} "
        f"against max|u|={gap['max_ref']:.6e}. That is a production-path finding, "
        f"not noise: measure the case before relaxing the bound."
    )
