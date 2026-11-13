"""Laminate ownership from original YAML, via public reader section queries."""
import copy

import numpy as np
import pytest
import yaml

from aeroelast.models.blade.numad.objects.blade import Blade
from tests.software.mesh.test_windio2_official import ROOT


def expected(source, eta, surface, lo, hi):
    structure = source["components"]["blade"]["structure"]
    anchors = {a["name"]: a for a in structure["anchors"]}
    anchors.update({a["name"]: a for w in structure["webs"] for a in w["anchors"]})

    def arc(field):
        while "anchor" in field:
            ref = field["anchor"]
            field = anchors[ref["name"]][ref["handle"]]
        if not field["grid"][0] <= eta <= field["grid"][-1]:
            return None
        return np.interp(eta, field["grid"], field["values"])

    result = []
    for layer in structure["layers"]:
        if layer.get("web", "skin") != surface:
            continue
        a, b = arc(layer["start_nd_arc"]), arc(layer["end_nd_arc"])
        if a is None or b is None or a == b:
            continue
        # Interior test independent of implementation's interval containment.
        q = (lo + hi) / 2
        if not ((a < q < b) if a < b else (q > a or q < b)):
            continue
        if not layer.get("start_nd_grid", 0) <= eta <= layer.get("end_nd_grid", 1):
            continue
        fields = [layer["thickness"], layer["fiber_orientation"]]
        if any(not f["grid"][0] <= eta <= f["grid"][-1] for f in fields):
            continue
        t, angle = [np.interp(eta, f["grid"], f["values"]) for f in fields]
        if t > 0:
            result.append((layer["name"], layer["material"], t, angle))
    return result


@pytest.mark.parametrize("eta", [.08163265306122448, .37, .95])
def test_official_all_cells_exact_source_stacks(official, eta):
    blade, source, _ = official
    section = blade.windio2_topology.section(eta)
    materials = {m["name"]: m for m in source["materials"]}
    assert len(section.skin_cells) == len(section.intervals)
    assert {c.surface for c in section.web_cells} == {w[0] for w in section.webs}
    for cell in (*section.skin_cells, *section.web_cells):
        lo, hi = cell.nd_bounds
        assert 0 <= lo < hi <= 1
        want = expected(source, eta, cell.surface, lo, hi)
        assert [(l.name, l.material_name) for l in cell.layers] == [(r[0], r[1]) for r in want]
        np.testing.assert_allclose([l.thickness_m for l in cell.layers], [r[2] for r in want], rtol=0, atol=1e-16)
        assert [l.fiber_angle_deg for l in cell.layers] == [r[3] for r in want]
        assert cell.thickness_m == pytest.approx(sum(r[2] for r in want), abs=1e-16)
        for layer in cell.layers:
            # Full source material record, including anisotropic E/G/nu and ply_t.
            for key, value in materials[layer.material_name].items():
                if isinstance(value, list):
                    assert layer.material[key] == tuple(value)
                else:
                    assert layer.material[key] == value
            with pytest.raises(TypeError):
                layer.material["rho"] = 0
    if eta == .08163265306122448:
        assert all(l.name != "LE_reinforcement" for c in section.skin_cells for l in c.layers)


def test_native_dropoff_web_local_and_wrap(tmp_path, official):
    _, original, _ = official
    source = copy.deepcopy(original)
    structure = source["components"]["blade"]["structure"]
    for name, surface, bounds in [("local", "web0", (.2, .6)), ("wrap", "skin", (.9, .1)),
                                  ("zero", "skin", (.3, .3))]:
        layer = copy.deepcopy(structure["layers"][0])
        layer.update(name=name, thickness={"grid": [.21, .43, .73], "values": [.0037, .0041, 0.]},
                     fiber_orientation={"grid": [.21, .73], "values": [-23., 41.]},
                     start_nd_grid=.27, end_nd_grid=.67)
        if surface != "skin":
            layer["web"] = surface
            # Explicit counts are metadata, never a reason to quantize thickness.
            layer["n_plies"] = {"grid": [.21, .73], "values": [3., 7.]}
        for key, value in zip(("start_nd_arc", "end_nd_arc"), bounds, strict=True):
            layer[key] = {"grid": [0., 1.], "values": [value, value]}
        structure["layers"].append(layer)
    path = tmp_path / "layup.yaml"
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    blade = Blade().read_yaml(str(path))
    topology = blade.windio2_topology
    coords = blade.geometry.coordinates.copy()
    assert {.21, .27, .43, .67, .73}.issubset(topology.event_stations)
    for eta in (.2, .27, .4, .67, .7, .74):
        section = topology.section(eta)
        for cell in (*section.skin_cells, *section.web_cells):
            want = expected(source, eta, cell.surface, *cell.nd_bounds)
            assert [l.name for l in cell.layers] == [r[0] for r in want]
            np.testing.assert_allclose([l.thickness_m for l in cell.layers], [r[2] for r in want], atol=1e-16, rtol=0)
            assert [l.fiber_angle_deg for l in cell.layers] == [r[3] for r in want]
        if eta == .4:
            local = [c for c in section.web_cells if any(l.name == "local" for l in c.layers)]
            assert len(local) == 1 and local[0].nd_bounds == (.2, .6)
            sample = next(l for l in local[0].layers if l.name == "local")
            assert sample.n_plies == np.interp(eta, [.21, .73], [3., 7.])
            assert sample.thickness_m != round(sample.thickness_m / sample.material["ply_t"]) * sample.material["ply_t"]
            assert {c.nd_bounds for c in section.skin_cells if any(l.name == "wrap" for l in c.layers)}
            assert all(not any(l.name == "local" for l in c.layers) for c in section.skin_cells)
            ends = dict(section.webs)["web0"]
            np.testing.assert_allclose(local[0].xyz, [np.array(ends[0].xyz) * (1-q) + np.array(ends[1].xyz) * q for q in (.2, .6)])
    accepted = topology.section(.4)
    for eta in (-.1, 1.1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match=r"^eta must be finite and in \[0, 1\]$"):
            topology.section(eta)
    assert topology.section(.4) == accepted
    assert path.read_bytes() == before
    np.testing.assert_array_equal(blade.geometry.coordinates, coords)
    assert topology.source["components"]["blade"]["structure"]["layers"][-1]["name"] == "zero"


def test_v1_none_preserved():
    blade = Blade().read_yaml(str(ROOT / "tests/IEA-15-240-RWT.yaml"))
    assert blade.windio2_topology is None
