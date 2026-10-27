"""Composite blade cases without ``elements.span_direction`` must still orient fibres along the span.

On the production blade mesh the quad edge 0->1 is chordwise (about 90 deg from global z), so a
missing ``span_direction`` would put every ply fibre across the chord. The loader therefore defaults
it to the span axis whenever the properties come from a blade/rotor generator; an explicit value
(and the isotropic path) is left untouched.
"""

import pytest

from aeroelast.core.config import FSISimulationConfig


def _data(generator="BladeMesh", elements=None, material=None):
    data = {
        "mesh": {
            "source": "generator",
            "generator": {"type": generator, "params": {"yaml_file": "blade.yaml"}},
        },
        "solver": {"type": "LinearDynamicFSIRotor", "total_time": 1.0, "time_step": 0.01},
        "boundary_conditions": {},
        "elements": {"family": "SHELL", **(elements or {})},
    }
    if material is not None:
        data["material"] = material
    return data


def _load(data):
    return FSISimulationConfig.from_dict(data)


@pytest.mark.parametrize("generator", ["BladeMesh", "RotorMesh"])
def test_composite_without_span_direction_defaults_to_z(generator):
    cfg = _load(_data(generator))
    assert list(cfg.elements.span_direction) == [0.0, 0.0, 1.0]


def test_explicit_span_direction_is_unchanged():
    cfg = _load(_data(elements={"span_direction": [1, 0, 0]}))
    assert list(cfg.elements.span_direction) == [1, 0, 0]


def test_isotropic_material_keeps_span_direction_unset():
    mat = {"type": "isotropic", "E": 1e9, "nu": 0.3, "rho": 1000.0}
    cfg = _load(_data(material=mat))
    assert cfg.elements.span_direction is None
