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


def _scaled_yaml(tmp_path, factor):
    data = yaml.safe_load(YAML.read_text())
    ref = data["components"]["blade"]["outer_shape_bem"]["reference_axis"]
    for ax in ("x", "y", "z"):
        ref[ax]["values"] = [v * factor for v in ref[ax]["values"]]
    path = tmp_path / "scaled.yaml"
    path.write_text(yaml.safe_dump(data))
    return path, ref["z"]["values"][-1]


def test_span_ends_at_the_yaml_tip_when_the_tip_is_not_an_integer(tmp_path):
    """G9: ``L = ceil(z_tip)`` stretched the whole span by up to 1 m."""
    from aeroelast.models.blade.numad.objects.blade import Blade

    path, tip = _scaled_yaml(tmp_path, 0.995)
    assert tip != np.ceil(tip)
    blade = Blade(str(path))
    assert blade.definition.span[-1] == pytest.approx(tip, rel=0, abs=1e-9)


def _layer(name):
    return {
        "name": name,
        "material": "m",
        "fiber_orientation": {"values": [0.0]},
        "thickness": {"grid": [0.0, 1.0], "values": [0.001, 0.001]},
    }


def _components(extra):
    names = ["spar_cap_ps", "spar_cap_ss", "uv", "shell_a", "shell_b", *extra]
    structure = {"layers": [_layer(n) for n in names]}
    definition = SimpleNamespace(materials={"m": SimpleNamespace(layerthickness=1.0)})
    return definition, structure


@pytest.mark.parametrize(
    "extra, match",
    [
        (["te_reinf_a", "te_reinf_b"], "te reinf"),
        (["le_reinf_a", "le_reinf_b"], "le reinf"),
    ],
)
def test_two_reinforcements_without_ss_ps_pair_raise(extra, match):
    """Four ``ValueError(...)`` statements were built and dropped, leaving the band unassigned."""
    valid = ["te_reinf_ss", "te_reinf_ps"] if "le" in extra[0] else ["le_reinf"]
    definition, structure = _components(valid + extra)
    with pytest.raises(ValueError, match=match):
        y2b._add_components(definition, structure, {})
