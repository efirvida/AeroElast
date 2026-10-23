"""Small NuMAD mesher defects of the blade input audit (``odd/tasks/blade-input-audit.md``, T15).

One test per item: G8 web control weights, G9 exact tip span, never-raised ``ValueError``s.
"""

import importlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.mesh_gen import mesh_gen

y2b = importlib.import_module("aeroelast.models.blade.numad.io.yaml_to_blade")
YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"


def test_web_edge_control_points_are_exact_thirds():
    """G8: 0.6666/0.3333 summed to 0.9999 and pulled every web control point 1e-4 |P| to the origin."""
    rng = np.random.default_rng(0)
    kp = np.zeros((16, 3))
    kp[[0, 1, 2, 3, 6, 7, 10, 11]] = 100.0 * rng.random((8, 3))
    mesh_gen._web_control_points(kp)
    third = lambda a, b: (2.0 * kp[a] + kp[b]) / 3.0  # noqa: E731
    for lo, hi, a, b in [(4, 5, 0, 1), (8, 9, 2, 3), (12, 13, 11, 6), (14, 15, 7, 10)]:
        np.testing.assert_allclose(kp[lo], third(a, b), rtol=1e-14, atol=0)
        np.testing.assert_allclose(kp[hi], third(b, a), rtol=1e-14, atol=0)


