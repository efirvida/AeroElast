"""Hand-computed per-length mass of the YAML inner skin, before meshing.

Regression target for the confirmed full-perimeter mapping defect. No solver,
mesh refinement, reference deck, or fitted tolerance is involved.
"""

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
import yaml

from aeroelast.models.blade.numad import Blade


def test_inner_skin_per_length_mass_respects_yaml_arc_interval():
    blade = Blade()
    blade.read_yaml(str(Path(__file__).parent / "IEA-15-240-RWT.yaml"))
    # At station 25 (eta=25/49), the official layer is exactly 2 mm thick,
    # rho=1940 kg/m^3, and occupies ND [0.08, 0.92]: 84% of the perimeter.
    station = 25
    arcs = blade.keypoints.key_arcs[:, station]
    perimeter = arcs[-1]-arcs[0]
    expected_kg_m = 1940.0 * 0.002 * (0.92-0.08) * perimeter
    actual_kg_m = 0.0
    for segment, stack in enumerate(blade.stackdb.stacks[:, station]):
        width = arcs[segment+1]-arcs[segment]
        for group in stack.plygroups:
            if group.component == "Shell_skin_inner":
                actual_kg_m += 1940.0 * .001 * group.thickness * group.nPlies * width
    assert actual_kg_m == pytest.approx(expected_kg_m, rel=1e-12)


def test_legacy_integer_ply_mapping_is_preserved():
    blade = Blade()
    blade.read_yaml(str(Path(__file__).parent / "IEA-15-240-RWT.yaml"))
    # Exercise the legacy NuMAD/Excel component contract with the same geometry.
    # Its 2.5 stock plies still round to two; continuous WindIO is opt-in metadata.
    for component in blade.definition.components.values():
        component.continuous_layers = False
    skin = blade.definition.components["Shell_skin"]
    skin.control_points[:, 1] = 2.5
    blade.update_blade()
    station = 10
    arcs = blade.keypoints.key_arcs[:, station]
    actual_kg_m = sum(
        1940. * .001 * group.thickness * group.nPlies * (arcs[j+1]-arcs[j])
        for j, stack in enumerate(blade.stackdb.stacks[:, station])
        for group in stack.plygroups if group.component == "Shell_skin"
    )
    assert actual_kg_m == pytest.approx(1940. * .002 * (arcs[-1]-arcs[0]), rel=1e-12)


def test_yaml_layer_mass_integrates_thickness_over_the_span_bay():
    blade = _load_skin([.001, .002])
    station = 10
    arcs = blade.keypoints.key_arcs[:, station]
    # Integral of h(eta)=0.001+0.001*eta across [10/49, 11/49]
    # per unit span is its midpoint value, NOT its inboard value.
    mean_h = .001 + .001 * (10.5/49)
    expected_kg_m = 1940. * mean_h * (arcs[-1]-arcs[0])
    actual_kg_m = sum(
        1940. * .001 * group.thickness * group.nPlies * (arcs[j+1]-arcs[j])
        for j, stack in enumerate(blade.stackdb.stacks[:, station])
        for group in stack.plygroups if group.component == "Shell_skin"
    )
    assert actual_kg_m == pytest.approx(expected_kg_m, rel=1e-12)


def _load_skin(thickness):
    source = Path(__file__).parent / "IEA-15-240-RWT.yaml"
    data = yaml.safe_load(source.read_text())
    layers = data["components"]["blade"]["internal_structure_2d_fem"]["layers"]
    layer = next(l for l in layers if l["name"] == "Shell_skin")
    layer["thickness"] = {"grid": [0., 1.], "values": thickness}
    layer["start_nd_arc"] = {"grid": [0., 1.], "values": [0., 0.]}
    layer["end_nd_arc"] = {"grid": [0., 1.], "values": [1., 1.]}
    with TemporaryDirectory(prefix="blade-mass-") as scratch:
        path = Path(scratch) / "blade.yaml"
        path.write_text(yaml.safe_dump(data))
        blade = Blade()
        blade.read_yaml(str(path))
    return blade


def test_yaml_layer_mass_preserves_non_integer_ply_thickness():
    # A uniform 2.5 mm full-perimeter layer of 1 mm stock has exactly
    # rho * 0.0025 * perimeter kg/m, not 2 or 3 rounded stock plies.
    blade = _load_skin([.0025, .0025])
    station = 10
    arcs = blade.keypoints.key_arcs[:, station]
    expected_kg_m = 1940. * .0025 * (arcs[-1]-arcs[0])
    actual_kg_m = sum(
        1940. * .001 * group.thickness * group.nPlies * (arcs[j+1]-arcs[j])
        for j, stack in enumerate(blade.stackdb.stacks[:, station])
        for group in stack.plygroups if group.component == "Shell_skin"
    )
    assert actual_kg_m == pytest.approx(expected_kg_m, rel=1e-12)
