"""``_matched_pairs`` keeps each row on one physical mode (synthetic, no blade).

Group 23 asserts one row per reference mode.  The helper used to sort the Hungarian pairs
by gap and return the ``n`` smallest, so entry ``k`` was the k-th best-matching pair rather
than reference mode ``k + 1``: the rows were the best of the matched set, a selection
biased towards green (the defect T29 removed from group 5, ``odd/tasks/l3-closure.md`` T13).
"""

from __future__ import annotations

import numpy as np

from tests.validation.blade.test_blade_iea15mw_validation import _matched_pairs


def test_entry_j_belongs_to_reference_mode_j() -> None:
    """Entry ``j`` pairs reference frequency ``j`` whatever the size of its gap."""
    ae = np.array([1.0, 2.0, 3.0, 4.0])
    ref = np.array([1.0001, 2.2, 3.0002, 4.4])  # gaps: tiny, 10 %, tiny, 10 %
    pairs = _matched_pairs(ae, ref, 4)
    assert [p[2] for p in pairs] == ref.tolist()
    assert [p[1] for p in pairs] == ae.tolist()
    np.testing.assert_allclose([p[0] for p in pairs], np.abs(ae - ref) / ref)


def test_first_n_reference_modes_not_the_n_best() -> None:
    """With ``n`` below the number of references, the first ``n`` modes are returned."""
    ae = np.array([1.0, 2.0, 3.0, 4.0])
    ref = np.array([1.0, 2.2, 3.0, 4.0])  # mode 1 is the worst, modes 2, 3 are exact
    pairs = _matched_pairs(ae, ref, 2)
    assert [p[2] for p in pairs] == [1.0, 2.2]
    assert pairs[1][0] > 0.05
