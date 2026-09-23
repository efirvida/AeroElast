"""Shared nodal load mapping for the one-way AD-load tools (S-5/S-6/S-8c).

Fixes the "1-2 nodes per section" bug: the old code computed the spanwise
tributary on the raw z-sorted node list.  All nodes of a cross-section share
the same z, so in the sorted array only the block's first/last entries got a
non-zero weight and each station's whole distributed load landed on 1-2
arbitrary nodes (mesh-order dependent), creating spurious section torques.

``spanwise_tributary`` groups nodes by (rounded) z, assigns each UNIQUE
station its trapezoidal span slice, and splits that slice evenly over the
nodes of the station.  The weights sum to the span and each station's total
load equals fn/ft * slice (verified: 831.5 kN vs 832.9 kN integral for the
S-5 rated loads).
"""

from __future__ import annotations

import numpy as np


def spanwise_tributary(z: np.ndarray, decimals: int = 9) -> np.ndarray:
    """Per-node spanwise tributary weights summing to (z.max() - z.min())."""
    z = np.asarray(z, dtype=float)
    uniq, inv = np.unique(np.round(z, decimals), return_inverse=True)
    spread = np.zeros(len(uniq))
    for k in range(len(uniq)):
        lo = uniq[k - 1] if k > 0 else uniq[k]
        hi = uniq[k + 1] if k < len(uniq) - 1 else uniq[k]
        spread[k] = 0.5 * (hi - lo)
    cnt = np.bincount(inv, minlength=len(uniq))
    return spread[inv] / cnt[inv]
