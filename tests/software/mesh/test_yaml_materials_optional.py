"""windIO v1 material parsing: strengths are optional, an isotropic shear modulus may be derived.

The IEA-3.4-130-RWT yaml (824f8d6) lists ``resin`` with only E, nu, rho and ``steel`` without G or S;
neither is used by the blade layup, but the reader parses every material and raised ``KeyError: 'Xt'``
(odd/tasks/l3-closure.md, T9b). Strengths do not enter the stiffness or the mass, so a missing one is
NaN; an isotropic material without G gets G = E / (2 (1 + nu)).
"""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aeroelast.models.blade.numad.io.yaml_to_blade import _add_materials


def _parse(materials):
    definition = SimpleNamespace()
    _add_materials(definition, materials)
    return definition.materials


def test_isotropic_material_without_strengths_or_shear_modulus():
    mats = _parse([{"name": "resin", "orth": 0, "E": 3.86e9, "nu": 0.35, "rho": 1150.0}])
    resin = mats["resin"]
    assert math.isnan(resin.uts) and math.isnan(resin.ucs) and math.isnan(resin.uss)
    g = 3.86e9 / (2.0 * 1.35)
    assert resin.gxy == pytest.approx(g) and resin.gxz == pytest.approx(g) and resin.gyz == pytest.approx(g)
    assert resin.ex == resin.ey == resin.ez == 3.86e9
    assert resin.density == 1150.0


def test_complete_material_is_read_as_before():
    mats = _parse([{
        "name": "glass_uni", "orth": 1, "rho": 1940.0, "ply_t": 0.001,
        "E": [4.46e10, 1.70e10, 1.67e10], "G": [3.27e9, 3.48e9, 3.50e9], "nu": [0.262, 0.35, 0.264],
        "Xt": [6.09e8, 3.81e7, 3.81e7], "Xc": [4.49e8, 1.04e8, 1.04e8], "S": [3.0e7, 3.0e7, 3.0e7],
    }])
    m = mats["glass_uni"]
    assert m.type == "orthotropic" and m.layerthickness == 1.0
    assert (m.ex, m.ey, m.ez) == (4.46e10, 1.70e10, 1.67e10)
    assert (m.gxy, m.gxz, m.gyz) == (3.27e9, 3.48e9, 3.50e9)
    assert (m.prxy, m.prxz, m.pryz) == (0.262, 0.35, 0.264)
    assert list(m.uts) == [6.09e8, 3.81e7, 3.81e7]
    assert list(m.ucs) == [-4.49e8, -1.04e8, -1.04e8]


def test_orthotropic_material_without_shear_modulus_still_fails():
    with pytest.raises(KeyError):
        _parse([{"name": "bad", "orth": 1, "rho": 1.0, "E": [1.0, 1.0, 1.0], "nu": [0.3, 0.3, 0.3]}])


IEA34_YAML = Path(__file__).resolve().parents[2] / "IEA-3.4-130-RWT.yaml"


def test_iea34_blade_is_read_with_every_layer():
    """The second validation blade reads end to end (fiber_orientation is optional in windIO v1)."""
    import yaml

    from aeroelast.models.blade.numad.objects.blade import Blade

    blade = Blade()
    blade.read_yaml(str(IEA34_YAML))
    layers = yaml.safe_load(IEA34_YAML.read_text())["components"]["blade"]["internal_structure_2d_fem"]["layers"]
    names = {c.name for c in blade.definition.components} if isinstance(blade.definition.components, list) \
        else set(blade.definition.components)
    assert names == {layer["name"] for layer in layers}
    assert all(float(np.mean(c.fabricangle)) == 0.0 for c in (
        blade.definition.components if isinstance(blade.definition.components, list)
        else blade.definition.components.values()))
