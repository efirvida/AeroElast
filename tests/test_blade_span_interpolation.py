"""Layer thickness is interpolated along the span linearly, the same scheme for every layer.

Candidate G10 of the blade input audit (``odd/tasks/blade-input-audit.md``, T9, user L15). NuMAD
interpolated a component with pchip when its name contained "spar" and linearly otherwise, and
resampled every layer onto the outer-shape grid with pchip. WISDEM uses one scheme for all layers
(``np.interp`` up to v3.2, the era of the 2020 IEA 15 MW BeamDyn deck); windIO prescribes none.
The audit reference (spec S2) is linear on the yaml's own grids. A pchip option is planned, not
implemented.
"""

from pathlib import Path

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.objects.blade import Blade

YAML = Path(__file__).parent / "IEA-15-240-RWT.yaml"
ETA = np.linspace(0.0, 1.0, 2001)


def _coarse_cap_yaml(tmp_path):
    """The official yaml with the spar caps on every other grid point (a grid unlike the outer shape)."""
    data = yaml.safe_load(YAML.read_text())
    for layer in data["components"]["blade"]["internal_structure_2d_fem"]["layers"]:
        if layer["name"].startswith("Spar_Cap"):
            for key in ("thickness",):
                g, v = layer[key]["grid"], layer[key]["values"]
                keep = sorted(set(range(0, len(g), 2)) | {len(g) - 1})
                layer[key]["grid"] = [g[k] for k in keep]
                layer[key]["values"] = [v[k] for k in keep]
    path = tmp_path / "coarse_caps.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def _check_linear(path):
    layers = {
        layer["name"]: layer
        for layer in yaml.safe_load(Path(path).read_text())["components"]["blade"][
            "internal_structure_2d_fem"
        ]["layers"]
    }
    blade = Blade()
    blade.read_yaml(str(path))
    definition = blade.definition
    components = definition.components
    components = components if isinstance(components, list) else list(components.values())
    checked = []
    for comp in components:
        layer = layers.get(comp.name)
        if layer is None or len(layer["thickness"]["grid"]) <= 2:
            continue
        ply = definition.materials[comp.materialid].layerthickness * 1e-3  # m
        got = comp.get_num_layers(ETA) * ply
        want = np.interp(ETA, layer["thickness"]["grid"], layer["thickness"]["values"])
        mask = want > 1e-6
        rel = np.abs(got[mask] - want[mask]) / want[mask]
        assert rel.max() < 1e-6, f"{comp.name}: {comp.imethod}, max rel {rel.max():.3e}"
        checked.append(comp.name)
    return checked


def test_every_layer_follows_the_yaml_linearly():
    checked = _check_linear(YAML)
    assert {"Spar_Cap_SS", "Spar_Cap_PS"} <= set(checked)
    assert len(checked) >= 10


def test_a_layer_grid_unlike_the_outer_shape_grid_stays_linear(tmp_path):
    checked = _check_linear(_coarse_cap_yaml(tmp_path))
    assert {"Spar_Cap_SS", "Spar_Cap_PS"} <= set(checked)


@pytest.mark.parametrize("name", ["Spar_Cap_SS", "Spar_Cap_PS", "Shell_skin"])
def test_component_interpolation_method_is_linear(name):
    blade = Blade()
    blade.read_yaml(str(YAML))
    components = blade.definition.components
    components = components if isinstance(components, list) else list(components.values())
    assert next(c for c in components if c.name == name).imethod == "linear"
