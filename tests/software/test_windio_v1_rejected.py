"""windIO 1.x blades are no longer read (odd/tasks/gmsh-official-mesher.md, S8/S9, T16).

Both readers that consume a blade yaml (the structural ``yaml_to_blade`` through ``normalize_windio`` and the
aerodynamic ``load_blade_aero``) must refuse a 1.x file with an explicit message instead of silently reading it
through the removed v1 branches.  The official v2 file keeps loading.
"""

import pytest

from aeroelast.models.blade.windio import normalize_windio
from tests.support.paths import V2_BLADE_YAML

V1_MESSAGE = "windIO 1.x blade files are no longer read"


@pytest.mark.parametrize("version", [None, "1.0", "1.1"])
def test_normalize_windio_rejects_v1(version):
    data = {"components": {"blade": {}}}
    if version is not None:
        data["windIO_version"] = version
    with pytest.raises(ValueError, match=V1_MESSAGE):
        normalize_windio(data)


def test_load_blade_aero_rejects_v1(tmp_path):
    from aeroelast.models.blade.aerodynamics import load_blade_aero

    path = tmp_path / "v1.yaml"
    path.write_text("windIO_version: '1.0'\ncomponents: {blade: {}}\n")
    with pytest.raises(ValueError, match=V1_MESSAGE):
        load_blade_aero(str(path), hub_radius=0.0)


def test_official_v2_still_normalizes():
    import yaml

    data = yaml.safe_load(V2_BLADE_YAML.read_text())
    out = normalize_windio(data)
    assert str(out["windIO_version"]).startswith("2.")
    assert out is not data
