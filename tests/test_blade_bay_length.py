"""Span bays are short enough for the per-bay laminate to follow a tapering layer.

Candidate B of the blade input audit (``odd/tasks/blade-input-audit.md``, T7): each NuMAD (region,
bay) holds one bay-mean laminate and the bays were the yaml grid step (about 2.4 m), so the spar
caps at r/R 0.90 were 11 % thin (audit F5 mass -6.86 %). Default 0.25 m (user L5).
"""

from pathlib import Path

import numpy as np

from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.models.blade.numad.objects.blade import Blade

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"


def test_default_bay_length_is_a_quarter_metre():
    assert BladeMesh(str(YAML), element_size=0.25).max_bay_length == 0.25


def test_limit_bay_length_splits_bays_and_keeps_the_yaml_stations():
    mesh = BladeMesh(str(YAML), element_size=0.25)
    mesh._numad_blade = Blade(str(YAML))
    before = np.array(mesh._numad_blade.ispan, dtype=float)
    mesh._limit_bay_length(0.5, verbose=False)
    after = np.asarray(mesh._numad_blade.ispan, dtype=float)
    assert np.all(np.isin(before, after))
    assert np.all(np.diff(after) <= 0.5 + 1e-9)
    assert mesh._numad_blade.keypoints.key_arcs.shape[1] == after.size
