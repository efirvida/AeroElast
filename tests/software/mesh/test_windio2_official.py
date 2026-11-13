"""Load the pinned official windIO 2.0 turbine without rewriting its inputs."""

import hashlib
from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

from aeroelast.models.blade.numad.objects.blade import Blade

ROOT = Path(__file__).resolve().parents[3]
OFFICIAL = ROOT / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
SHA256 = "af51499e4d4b35293f9c6bca7701a9af6df0b3884c532891b1cc133efd420978"


def source_reference_length(axis):
    """Independent schema reference_axis integral using only source data."""
    curves = [PchipInterpolator(axis[c]["grid"], axis[c]["values"]) for c in "xyz"]
    knots = sorted({v for c in "xyz" for v in axis[c]["grid"]})
    def speed(eta):
        return np.sqrt(sum(curve.derivative()(eta)**2 for curve in curves))
    return sum(quad(speed, a, b)[0] for a, b in zip(knots[:-1], knots[1:], strict=True))


def test_official_geometry_materials_and_input_identity():
    before = OFFICIAL.read_bytes()
    source = yaml.safe_load(before)
    blade = Blade()
    blade.read_yaml(str(OFFICIAL))
    definition = blade.definition
    geometry = source["components"]["blade"]
    z = definition.span  # Native distributions; ispan additionally holds inserted chord knots.
    length = source_reference_length(geometry["reference_axis"])
    np.testing.assert_allclose(z, np.array(geometry["reference_axis"]["z"]["grid"]) * length,
                               rtol=1e-12, atol=1e-12)
    twist = geometry["outer_shape"]["twist"]
    np.testing.assert_allclose(
        definition.degreestwist,
        PchipInterpolator(twist["grid"], twist["values"])(z / length),
        rtol=1e-12, atol=1e-12,
    )
    chord = geometry["outer_shape"]["chord"]
    np.testing.assert_array_equal(
        definition.chord_curve[0], np.array(chord["grid"]) * definition.span[-1]
    )
    np.testing.assert_array_equal(definition.chord_curve[1], chord["values"])
    # The knot at eta=.995 must not disappear during format conversion.
    assert definition.chord_curve[1][-2] == 1.472482673397968
    assert definition.materials["adhesive"].ex == 4.56e9
    assert definition.materials["adhesive"].gxy == 1.52e9
    assert {c.name for c in definition.components.values()} == {
        layer["name"] for layer in geometry["structure"]["layers"]
    }
    assert OFFICIAL.read_bytes() == before
    assert hashlib.sha256(before).hexdigest() == SHA256


def test_official_anchor_arcs_widths_and_layer_thickness():
    source = yaml.safe_load(OFFICIAL.read_text())
    structure = source["components"]["blade"]["structure"]
    anchors = {a["name"]: a for a in structure["anchors"]}
    blade = Blade()
    blade.read_yaml(str(OFFICIAL))
    definition = blade.definition
    # NuMAD group 1 is the fore web1, group 2 the aft web0.
    for row, name in zip(definition.web_nd_arc, ["web1", "web0"], strict=True):
        for idx, key in [(1, "start_nd_arc"), (2, "end_nd_arc")]:
            expected = anchors[name][key]
            np.testing.assert_array_equal(row[0], expected["grid"])
            np.testing.assert_array_equal(row[idx], expected["values"])
    np.testing.assert_allclose(
        definition.sparcapwidth_lp,
        1000 * np.interp(
            definition.ispan / source_reference_length(source["components"]["blade"]["reference_axis"]),
            anchors["Spar_Cap_SS"]["width"]["grid"],
            anchors["Spar_Cap_SS"]["width"]["values"],
        ),
    )
    # windIO 63dcad6 schema nd_arc_position and guide: both endpoints are the
    # TE midpoint, including both legs of an open TE (not just its corners).
    geometry = blade.geometry
    eta = geometry.ispan / geometry.ispan[-1]
    hp_te_leg = geometry.arclength[1] - geometry.arclength[0]
    lp_te_leg = geometry.arclength[-1] - geometry.arclength[-2]
    # The open TE is closed to its midpoint when the profiles are read, so both
    # legs vanish; the arcs still start at that midpoint.
    assert np.all(hp_te_leg == 0) and np.all(lp_te_leg == 0), "TE must be closed at read"
    arc0 = geometry.arclength[0]
    perimeter = geometry.arclength[-1] - arc0
    for name, rows in [("Spar_Cap_PS", (4, 3)), ("Spar_Cap_SS", (9, 8))]:
        for key, row in zip(("start_nd_arc", "end_nd_arc"), rows, strict=True):
            source_arc = anchors[name][key]
            expected = arc0 + (1 - np.interp(eta, source_arc["grid"], source_arc["values"])) * perimeter
            np.testing.assert_allclose(blade.keypoints.key_arcs[row], expected)
    inner = definition.components["Shell_skin_inner"]
    np.testing.assert_array_equal(
        inner.nd_arc_start["values"], anchors["Shell_skin_inner"]["start_nd_arc"]["values"]
    )
    for layer in structure["layers"]:
        component = definition.components[layer["name"]]
        material = definition.materials[layer["material"]]
        thickness = np.interp(
            component.control_points[:, 0],
            layer["thickness"]["grid"], layer["thickness"]["values"],
            left=0, right=0,
        )
        np.testing.assert_allclose(
            component.control_points[:, 1] * material.layerthickness / 1000,
            thickness,
            # Conversion/interpolation roundoff only, not a physics pass band.
            atol=np.finfo(float).eps * max(layer["thickness"]["values"]),
        )


def test_invalid_anchor_does_not_replace_blade_definition_or_input(tmp_path):
    source = yaml.safe_load(OFFICIAL.read_text())
    source["components"]["blade"]["structure"]["layers"][0]["start_nd_arc"] = {
        "anchor": {"name": "missing", "handle": "start_nd_arc"}
    }
    path = tmp_path / "invalid.yaml"
    path.write_text(yaml.safe_dump(source))
    before = path.read_bytes()
    blade = Blade()
    definition = blade.definition
    with pytest.raises(ValueError, match="Unknown windIO anchor: missing"):
        blade.read_yaml(str(path))
    assert blade.definition is definition
    assert path.read_bytes() == before


def test_official_aerodynamic_geometry_and_polar_units():
    from aeroelast.models.blade.aerodynamics import load_blade_aero

    source = yaml.safe_load(OFFICIAL.read_text())
    before = OFFICIAL.read_bytes()
    aero = load_blade_aero(str(OFFICIAL))
    outer = source["components"]["blade"]["outer_shape"]
    twist = outer["twist"]
    np.testing.assert_allclose(aero.twist, np.deg2rad(np.interp(
        outer["chord"]["grid"], twist["grid"], twist["values"],
    )))
    np.testing.assert_array_equal(aero.chord, outer["chord"]["values"])
    assert len(aero.airfoils[0].polars) == 1
    np.testing.assert_allclose(aero.airfoils[0].polars[0].alpha[[0, -1]], [-np.pi, np.pi])
    assert OFFICIAL.read_bytes() == before


@pytest.mark.parametrize("name", ["IEA-15-240-RWT.yaml", "IEA-3.4-130-RWT.yaml"])
def test_legacy_official_inputs_still_load_without_modification(name):
    path = ROOT / "tests" / name
    before = path.read_bytes()
    blade = Blade()
    blade.read_yaml(str(path))
    assert blade.definition.ispan[-1] == (117 if "15-240" in name else 63)
    assert blade.definition.web_nd_arc is not None
    assert path.read_bytes() == before
