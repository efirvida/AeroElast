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


def _skin_bounds_yaml(tmp_path):
    """The official yaml with Shell_skin_inner's chordwise bounds varying non-linearly on a coarse grid."""
    data = yaml.safe_load(YAML.read_text())
    for layer in data["components"]["blade"]["internal_structure_2d_fem"]["layers"]:
        if layer["name"] == "Shell_skin_inner":
            start = [0.0, 0.02, 0.2, 0.2]
            layer["start_nd_arc"] = {"grid": [0.0, 0.3, 0.6, 1.0], "values": start}
            layer["end_nd_arc"] = {"grid": [0.0, 0.3, 0.6, 1.0], "values": [1.0 - s for s in start]}
    path = tmp_path / "skin_bounds.yaml"
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def _check_nd_arc_linear(path, monkeypatch):
    """Every nd_arc bound the stack database evaluates along the span follows the yaml linearly (G10b)."""
    from aeroelast.models.blade.numad.objects import stackdb

    layers = {
        layer["name"]: layer
        for layer in yaml.safe_load(Path(path).read_text())["components"]["blade"][
            "internal_structure_2d_fem"
        ]["layers"]
    }
    real = stackdb.interpolator_wrap
    calls = []

    def spy(x, v, xq, method="linear", axis=0):
        out = real(x, v, xq, method, axis)
        calls.append(
            (np.asarray(x, float), np.asarray(v, float), np.asarray(xq, float), out, method)
        )
        return out

    monkeypatch.setattr(stackdb, "interpolator_wrap", spy)
    blade = Blade()
    blade.read_yaml(str(path))
    components = blade.definition.components.values()
    shells = [c for c in components if c.nd_arc_start is not None]
    assert {"Shell_skin", "Shell_skin_inner"} <= {c.name for c in shells}
    arcs = [c for c in calls if c[2].ndim == 1 and c[4] in ("pchip", "linear") and c[0].size > 2]
    assert len(arcs) >= 2 * len(shells)
    for comp in shells:
        for key in ("start_nd_arc", "end_nd_arc"):
            want_grid = layers[comp.name][key]
            hit = [
                c
                for c in calls
                if np.allclose(
                    c[0],
                    comp.nd_arc_start["grid"] if key == "start_nd_arc" else comp.nd_arc_end["grid"],
                )
                and c[4] is not None
                and c[2].size > 1
                and np.allclose(
                    c[1],
                    (comp.nd_arc_start if key == "start_nd_arc" else comp.nd_arc_end)["values"],
                )
            ]
            assert hit, f"{comp.name} {key}: never evaluated"
            _, _, xq, out, method = hit[0]
            ref = np.interp(xq, want_grid["grid"], want_grid["values"])
            assert method == "linear", f"{comp.name} {key}: {method}"
            assert np.abs(out - ref).max() < 1e-9, (
                f"{comp.name} {key}: {np.abs(out - ref).max():.3e}"
            )


def test_official_skin_bounds_follow_the_yaml_linearly(monkeypatch):
    _check_nd_arc_linear(YAML, monkeypatch)


def test_nonlinear_skin_bounds_on_a_coarse_grid_stay_linear(tmp_path, monkeypatch):
    _check_nd_arc_linear(_skin_bounds_yaml(tmp_path), monkeypatch)
