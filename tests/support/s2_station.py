"""Interior-station helpers of the S-2 static comparison (l3-closure S12, T10c-station).

Pure numpy / json, no mesh and no solver, so they unit-test in milliseconds
(``tests/support/test_s2_station.py``).  The rule they encode is the one predeclared in
``odd/tasks/l3-closure.md`` ("T10c-station predeclaration"): the BeamDyn p60 output node used is the inner of the
two GLL nodes nearest the tip whose tip distance is at least ``MIN_CHORDS`` tip chords, its z is read from the
BeamDyn reference output (never the GLL estimate), and the shell value at that z is the arc-weighted ring mean
interpolated linearly in z between the two bracketing rings.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

MIN_CHORDS = 2.0
BD_SCALES = (1e-4, 2e-4)
_CHORD_EPS = 1e-9  # a node exactly MIN_CHORDS from the tip counts


def component_index(case: str) -> int:
    """Dominant shell translation component: ux (0) for the edgewise tip load LC3, uy (1) otherwise."""
    return 0 if case == "LC3" else 1


def bd_station_nodes(root: Path, case: str, order: int = 60, scales=BD_SCALES) -> tuple[list[dict], float]:
    """BeamDyn output nodes of ``<case>-p<order>-s<scale>-nodes/result.json`` and the scale-to-scale spread.

    The nodes (shell axes: ``linear_shell_uy_m`` = N#TDxr / scale, ``linear_shell_ux_m`` = -N#TDyr / scale) are
    those of the first scale; the spread is the largest relative difference of the dominant component between
    the scales, the linearity gate of the reference.
    """
    runs = [
        json.loads((Path(root) / f"{case}-p{order}-s{scale:g}-nodes" / "result.json").read_text())["nodes"]
        for scale in scales
    ]
    key = f"linear_shell_u{'xy'[component_index(case)]}_m"
    spread = max(
        abs(a[key] / b[key] - 1.0) for other in runs[1:] for a, b in zip(runs[0], other, strict=True)
    )
    return runs[0], spread


def pick_station_node(nodes: list[dict], z_tip: float, chord: float, min_chords: float = MIN_CHORDS) -> dict:
    """Inner of the two nodes nearest the tip that lie at least ``min_chords`` chords from it.

    ``z_tip`` is the shell tip and ``chord`` its tip chord (``tip_ring_scale``); each node dict carries the
    ``z`` of the BeamDyn reference output.
    """
    far = [n for n in nodes if z_tip - n["z"] >= min_chords * chord * (1.0 - _CHORD_EPS)]
    if len(far) < 2:
        raise ValueError(f"fewer than two BeamDyn nodes at least {min_chords} chords ({min_chords * chord:.3f} m) inboard")
    return sorted(far, key=lambda n: -n["z"])[1]


def ring_mean_at_z(profile, z: float) -> np.ndarray:
    """Arc ring-mean translation (ux, uy, uz) linear in z between the two bracketing rings.

    ``profile`` rows are ``[ring z, ux, uy, uz, node count]`` (``tools/s2_ring_metric_blade.profile_means``).
    """
    prof = np.array(sorted(profile, key=lambda r: r[0]), dtype=float)
    if not prof[0, 0] <= z <= prof[-1, 0]:
        raise ValueError(f"z {z:.3f} outside the ring profile [{prof[0, 0]:.3f}, {prof[-1, 0]:.3f}]")
    return np.array([np.interp(z, prof[:, 0], prof[:, 1 + k]) for k in range(3)])
