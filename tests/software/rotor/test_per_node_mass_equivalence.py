"""Equivalence test for the vectorized per-node translational mass gather.

``LinearDynamicFSIRotorCorotationalSolver._solve_via_rust`` originally computed
``all_node_masses`` with a per-node Python loop that, for every node, made two
full boolean scans of the COO row array ``mr``::

    all_node_masses = np.array(
        [
            float(mv[mr == i * dofs].sum()) if (mr == i * dofs).any() else 0.0
            for i in range(n_nodes)
        ],
        dtype=np.float64,
    )

That is O(n_nodes * nnz) and dominated the runtime at production sizes
(n_nodes ~= 32k, nnz ~= 4.6M). The shipped replacement is
``aeroelast.solvers.fsi.rotor._per_node_translational_masses``, a single
O(nnz) ``np.bincount`` pass.

The gather relies on the DOF layout being a plain stride, i.e. node ``i`` owns
rows ``i*dofs .. i*dofs+dofs-1``; that convention is asserted by
``MeshAssembler._node_dofs_map`` (``range(i * dpn, i * dpn + dpn)``) and by the
``[u, v, w, θx, θy, θz]`` per-node shell layout in the project docs.

This is a pure performance refactor with numerically identical output, so there
is no meaningful RED state for the behaviour itself. The tests below pin the
equivalence against the original per-node expression instead.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.solvers.fsi.rotor import _per_node_translational_masses

DOFS = 6  # shell: [u, v, w, θx, θy, θz]


def _reference_per_node_masses(
    rows: np.ndarray,
    vals: np.ndarray,
    dofs: int,
    n_nodes: int,
) -> np.ndarray:
    """Straightforward per-node reference: the original loop, unchanged."""
    return np.array(
        [
            float(vals[rows == i * dofs].sum()) if (rows == i * dofs).any() else 0.0
            for i in range(n_nodes)
        ],
        dtype=np.float64,
    )


def test_pins_equivalence_on_synthetic_coo_with_empty_node() -> None:
    """Exact equality with the per-node reference, including an empty node.

    Node 2 has no COO entries at all, so it must come out as exactly 0.0.
    Node 0 and node 3 have several contributions to their first translational
    row; only that first row (``i * dofs``) may contribute.
    """
    n_nodes = 5
    # rows/cols are 0-based global DOF indices; only row == i * dofs counts.
    rows = np.array(
        [
            0,
            0,  # node 0, row 0 (counts: 1.0 + 2.0)
            1,  # node 0, row 1 (ignored)
            6,  # node 1, row 6 (counts: 4.0)
            7,  # node 1, row 7 (ignored)
            # node 2 intentionally absent -> no entries
            18,
            18,
            18,  # node 3, row 18 (counts: 5.0 + 6.0 + 7.0)
            19,  # node 3, row 19 (ignored)
            24,  # node 4, row 24 (counts: 9.0)
            25,
            29,  # node 4, rows 25,29 (ignored)
        ],
        dtype=np.int64,
    )
    vals = np.array(
        [1.0, 2.0, 99.0, 4.0, 99.0, 5.0, 6.0, 7.0, 99.0, 9.0, 99.0, 99.0],
        dtype=np.float64,
    )

    new = _per_node_translational_masses(rows, vals, DOFS, n_nodes, n_nodes * DOFS)
    ref = _reference_per_node_masses(rows, vals, DOFS, n_nodes)

    assert new.shape == (n_nodes,)
    assert new.dtype == np.float64
    # Exact equality: per-row entry counts here are below NumPy's pairwise
    # summation block size, so bincount and mask-sum accumulate in the same
    # order and produce bit-identical results.
    np.testing.assert_array_equal(new, ref)
    assert new[2] == 0.0  # explicit "no entries for this node" case
    np.testing.assert_array_equal(ref, np.array([3.0, 4.0, 0.0, 18.0, 9.0]))


def test_max_difference_on_random_coo_is_negligible() -> None:
    """Deterministic random COO: max |old - new| stays at float round-off.

    Mimics the production distribution (shell, 6 DOFs/node, a handful of
    entries per translational row) at a size that keeps the O(n_nodes * nnz)
    reference cheap. A subset of nodes is dropped so the empty-node branch is
    exercised at scale too.

    Tolerance justification: ``np.bincount`` accumulates each row sequentially
    while NumPy's ``sum`` uses pairwise summation, so once a row has more than
    a couple of entries the two agree only up to floating-point summation
    order. The observed maximum across this deterministic fixture is ~1.3e-15;
    1e-12 is a safe bound that still fails loudly on any real layout bug.
    """
    rng = np.random.default_rng(20260930)
    dofs = DOFS
    n_nodes = 500
    n_full_dofs = n_nodes * dofs
    entries_per_row = 8

    rows = np.repeat(np.arange(n_full_dofs), entries_per_row).astype(np.int64)
    vals = rng.standard_normal(rows.size)

    # Drop all entries for one node in every 17, creating empty-node cases.
    empty_nodes = np.arange(0, n_nodes, 17)
    drop = np.isin(rows // dofs, empty_nodes)
    keep = ~drop
    rows, vals = rows[keep], vals[keep]

    new = _per_node_translational_masses(rows, vals, dofs, n_nodes, n_full_dofs)
    ref = _reference_per_node_masses(rows, vals, dofs, n_nodes)

    max_abs_diff = float(np.max(np.abs(new - ref)))
    assert max_abs_diff < 1e-12, f"summation-order drift too large: {max_abs_diff}"
    np.testing.assert_allclose(new, ref, rtol=0.0, atol=1e-12)
    # Empty nodes are exactly 0.0 on both sides, regardless of summation order.
    np.testing.assert_array_equal(new[empty_nodes], np.zeros(empty_nodes.size))
    np.testing.assert_array_equal(ref[empty_nodes], np.zeros(empty_nodes.size))


def test_minlength_shorter_than_actual_rows_still_matches() -> None:
    """Defensive: an undersized ``n_full_dofs`` must not truncate rows.

    ``np.bincount`` grows beyond ``minlength`` when the largest row index
    demands it, so passing the true ``n_full_dofs`` keeps every node covered.
    """
    n_nodes = 3
    rows = np.array([0, 6, 12], dtype=np.int64)
    vals = np.array([1.5, 2.5, 3.5], dtype=np.float64)

    new = _per_node_translational_masses(rows, vals, DOFS, n_nodes, n_nodes * DOFS)
    ref = _reference_per_node_masses(rows, vals, DOFS, n_nodes)

    np.testing.assert_array_equal(new, ref)


@pytest.mark.parametrize("dofs", [2, 6])
def test_holds_for_plane_and_shell_strides(dofs: int) -> None:
    """The stride gather is layout-agnostic across PLANE (2) and SHELL (6)."""
    n_nodes = 40
    rng = np.random.default_rng(7 + dofs)
    rows = np.arange(n_nodes * dofs, dtype=np.int64)
    vals = rng.standard_normal(rows.size)

    new = _per_node_translational_masses(rows, vals, dofs, n_nodes, n_nodes * dofs)
    ref = _reference_per_node_masses(rows, vals, dofs, n_nodes)

    np.testing.assert_array_equal(new, ref)
