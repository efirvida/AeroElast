"""Pure functions of ``tools/audit_blade_input.py`` (odd/tasks/blade-input-audit.md, T2).

The audit compares the 0.25 m shell mesh with ``tests/IEA-15-240-RWT.yaml`` through a standalone
yaml reader. Every function tested here is checked against a case whose answer is known in closed
form, before the harness is trusted on the real blade. The mesh side (``build_mesh`` and the
station slicing) needs the 60 s production mesh and is exercised by the tool run, not here.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import audit_blade_input as ab  # noqa: E402


def mini_yaml():
    """Two-airfoil, four-layer, one-web blade with round numbers."""
    lin = lambda g, v: {"grid": list(g), "values": list(v)}  # noqa: E731
    return {
        "airfoils": [{"name": "A", "coordinates": {"x": [1, 0.5, 0, 0.5, 1], "y": [0.01, 0.1, 0, -0.1, -0.01]}}],
        "materials": [
            {"name": "iso", "orth": 0, "E": 200e9, "G": 76.9e9, "nu": 0.3, "rho": 1000.0},
            {"name": "uni", "orth": 1, "E": [40e9, 10e9, 10e9], "G": [4e9, 4e9, 3e9],
             "nu": [0.3, 0.3, 0.4], "rho": 2000.0},
        ],
        "components": {"blade": {
            "outer_shape_bem": {
                "chord": lin([0, 1], [4.0, 2.0]),
                "twist": lin([0, 0.5, 1], [0.2, 0.1, 0.0]),
                "pitch_axis": lin([0, 1], [0.5, 0.4]),
                "airfoil_position": {"grid": [0, 1], "labels": ["A", "A"]},
                "reference_axis": {"x": lin([0, 1], [0.0, -4.0]), "y": lin([0, 1], [0, 0]),
                                   "z": lin([0, 1], [0.0, 100.0])},
            },
            "internal_structure_2d_fem": {
                "layers": [
                    {"name": "skin", "material": "iso", "thickness": lin([0, 1], [0.004, 0.002]),
                     "start_nd_arc": lin([0, 1], [0.0, 0.0]), "end_nd_arc": lin([0, 1], [1.0, 1.0]),
                     "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                    {"name": "cap", "material": "uni", "thickness": lin([0, 1], [0.02, 0.0]),
                     "start_nd_arc": lin([0, 1], [0.2, 0.3]), "end_nd_arc": lin([0, 1], [0.4, 0.4]),
                     "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                    {"name": "gone", "material": "iso", "thickness": lin([0, 1], [0.0, 0.0]),
                     "start_nd_arc": lin([0, 1], [0.0, 0.0]), "end_nd_arc": lin([0, 1], [1.0, 1.0]),
                     "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                    {"name": "te_ss", "material": "iso", "thickness": lin([0, 1], [0.01, 0.01]),
                     "start_nd_arc": lin([0, 1], [0.0, 0.0]), "width": lin([0, 1], [0.5, 0.5]),
                     "side": "suction", "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                    {"name": "te_ps", "material": "iso", "thickness": lin([0, 1], [0.01, 0.01]),
                     "end_nd_arc": lin([0, 1], [1.0, 1.0]), "width": lin([0, 1], [0.5, 0.5]),
                     "side": "pressure", "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                    {"name": "web0_skin", "material": "uni", "web": "web0",
                     "thickness": lin([0, 1], [0.003, 0.001]),
                     "fiber_orientation": lin([0, 1], [0.0, 0.0])},
                ],
                "webs": [{"name": "web0", "start_nd_arc": lin([0, 1], [0.2, 0.3]),
                          "end_nd_arc": lin([0, 1], [0.7, 0.6])}],
            },
        }},
    }


# --- span interpolation on the yaml's own grids ----------------------------------------------
def test_interp_is_linear_on_the_given_grid_and_clamped():
    g, v = [0.0, 0.5, 1.0], [0.0, 2.0, 2.0]
    assert ab.interp(g, v, 0.25) == pytest.approx(1.0, abs=1e-15)
    assert ab.interp(g, v, 0.75) == pytest.approx(2.0, abs=1e-15)
    assert ab.interp(g, v, -1.0) == 0.0 and ab.interp(g, v, 2.0) == 2.0


def test_outer_shape_uses_each_quantity_own_grid_and_converts_twist_to_degrees():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    o = yb.outer_shape(0.25)
    assert o["chord"] == pytest.approx(3.5, abs=1e-12)  # 4 -> 2 linear
    assert o["twist_deg"] == pytest.approx(np.degrees(0.15), abs=1e-12)  # 0.2 -> 0.1 on [0, .5]
    assert o["pitch_axis"] == pytest.approx(0.475, abs=1e-12)
    assert o["prebend"] == pytest.approx(-1.0, abs=1e-12)  # reference_axis.x
    assert o["sweep"] == 0.0  # reference_axis.y
    assert o["z"] == pytest.approx(25.0, abs=1e-12)


def test_blade_length_is_the_reference_axis_z_end():
    assert ab.YamlBlade.from_dict(mini_yaml()).length == pytest.approx(100.0)


def test_layer_state_interpolates_thickness_and_arcs_linearly():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    cap = yb.layer_state("cap", 0.5)
    assert cap["thickness"] == pytest.approx(0.01, abs=1e-15)
    assert cap["start"] == pytest.approx(0.25, abs=1e-15)
    assert cap["end"] == pytest.approx(0.4, abs=1e-15)
    assert cap["material"] == "uni" and cap["angle_deg"] == 0.0


# --- nd_arc convention: 0 at the suction-side trailing edge, 1 at the pressure-side one -------
def test_ring_arcs_is_cumulative_length_over_total():
    pts = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 2.0]])  # lengths 1 and 2
    arcs, total = ab.ring_arcs(pts)
    assert total == pytest.approx(3.0)
    assert arcs == pytest.approx([0.0, 1.0 / 3.0, 1.0])


def test_orient_ring_starts_on_the_suction_side():
    t = np.linspace(0.0, 2 * np.pi, 81)
    # ellipse, TE (suction, y > 0 first) -> LE -> TE (pressure), as the yaml airfoils are ordered
    prof = np.c_[0.5 * (1 + np.cos(t)), 0.06 * np.sin(t)]
    ring = ab.profile_to_mesh_ring(prof, chord=2.0, pitch_axis=0.5)
    oriented = ab.orient_ring(ring[::-1])
    assert oriented == pytest.approx(ring)  # the reversed ring is flipped back


def test_active_layers_follow_yaml_order_and_drop_zero_thickness_and_out_of_range():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    perimeter = 10.0
    names = lambda r, a: [x["name"] for x in yb.active_layers(r, a, perimeter)]  # noqa: E731
    # r = 0: cap spans 0.2-0.4 with thickness 0.02; "gone" has zero thickness; te_ss width 0.5 m
    # of a 10 m perimeter -> arc 0..0.05; te_ps -> 0.95..1.
    assert names(0.0, 0.3) == ["skin", "cap"]
    assert names(0.0, 0.03) == ["skin", "te_ss"]
    assert names(0.0, 0.97) == ["skin", "te_ps"]
    assert names(0.0, 0.5) == ["skin"]
    # r = 1: the cap has thickness 0 -> inactive even inside its arc
    assert names(1.0, 0.35) == ["skin"]


def test_missing_bound_is_filled_from_width_over_perimeter():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    s, e = yb.layer_bounds("te_ss", 0.0, perimeter=10.0)
    assert (s, e) == (pytest.approx(0.0), pytest.approx(0.05))
    s, e = yb.layer_bounds("te_ps", 0.0, perimeter=10.0)
    assert (s, e) == (pytest.approx(0.95), pytest.approx(1.0))


def test_web_state_and_web_layers():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    w = yb.web_state(0, 0.5)
    assert w["start"] == pytest.approx(0.25) and w["end"] == pytest.approx(0.65)
    layers = yb.web_layers(0, 0.5)
    assert [x["name"] for x in layers] == ["web0_skin"]
    assert layers[0]["thickness"] == pytest.approx(0.002)


def test_material_constants_for_isotropic_and_orthotropic_entries():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    iso = yb.material("iso")
    assert (iso["e1"], iso["e2"], iso["g12"], iso["nu12"], iso["rho"]) == (200e9, 200e9, 76.9e9, 0.3, 1000.0)
    uni = yb.material("uni")
    assert (uni["e1"], uni["e2"], uni["g12"], uni["nu12"], uni["rho"]) == (40e9, 10e9, 4e9, 0.3, 2000.0)
    # the other windIO slots are kept under their windIO names (xy, xz, yz)
    assert (uni["g13"], uni["g23"], uni["nu13"], uni["nu23"]) == (4e9, 3e9, 0.3, 0.4)


# --- mass and axial integrals on a section with a closed form --------------------------------
def test_section_integrals_closed_form():
    # ply A: rho 1000, t 2 mm, E1 = 200 GPa, fibre along the span; ply B: rho 500, t 10 mm, E1 = 10 GPa
    plies = [{"rho": 1000.0, "t": 0.002, "e1": 200e9, "e2": 200e9, "g12": 76.9e9, "nu12": 0.3, "angle_deg": 0.0},
             {"rho": 500.0, "t": 0.010, "e1": 10e9, "e2": 10e9, "g12": 4e9, "nu12": 0.25, "angle_deg": 0.0}]
    mass, axial = ab.section_integrals([(plies, 2.0)])  # one element, ds = 2 m
    assert mass == pytest.approx(2.0 * (1000 * 0.002 + 500 * 0.010), rel=1e-14)  # 14 kg/m
    assert axial == pytest.approx(2.0 * (200e9 * 0.002 + 10e9 * 0.010), rel=1e-14)  # 1000 MN


def test_axial_modulus_of_a_ply_rotated_off_the_span():
    e1, e2, g12, nu12 = 40e9, 10e9, 4e9, 0.3
    assert ab.axial_modulus(e1, e2, g12, nu12, 0.0) == pytest.approx(e1, rel=1e-14)
    assert ab.axial_modulus(e1, e2, g12, nu12, 90.0) == pytest.approx(e2, rel=1e-12)
    # an isotropic ply has Ex = E at any angle
    e, nu = 200e9, 0.3
    assert ab.axial_modulus(e, e, e / (2 * (1 + nu)), nu, 37.0) == pytest.approx(e, rel=1e-12)
    # 45 degrees, closed form 1/(1/4 (1/E1 + 1/E2) + 1/4 (1/G12 - 2 nu12/E1))
    expect = 1.0 / (0.25 * (1 / e1 + 1 / e2) + 0.25 * (1 / g12 - 2 * nu12 / e1))
    assert ab.axial_modulus(e1, e2, g12, nu12, 45.0) == pytest.approx(expect, rel=1e-12)


# --- outer shape recovered from a synthetic ring ---------------------------------------------
def test_measure_outer_recovers_chord_twist_and_pitch_point():
    chord, pa, twist_deg, prebend, sweep = 3.0, 0.4, 7.5, -1.25, 0.0
    t = np.linspace(0.0, 2 * np.pi, 201)
    # ellipse airfoil in profile coordinates: x in [0, 1] LE -> TE, y thickness 12 %
    x = 0.5 * (1 + np.cos(t))  # starts at TE (1), via LE (0), back to TE
    y = 0.06 * np.sin(t)
    prof = np.c_[x, y]
    ring = ab.profile_to_mesh_ring(prof, chord=chord, pitch_axis=pa, twist_deg=twist_deg,
                                   prebend=prebend, sweep=sweep)
    got = ab.measure_outer(ab.orient_ring(ring), pitch_axis=pa)
    assert got["chord"] == pytest.approx(chord, rel=1e-3)
    assert got["twist_deg"] == pytest.approx(twist_deg, abs=0.05)
    assert got["px"] == pytest.approx(-sweep, abs=1e-2)
    assert got["py"] == pytest.approx(prebend, abs=1e-2)


# --- airfoil shape distance ------------------------------------------------------------------
def test_max_normal_distance_against_a_known_offset_polyline():
    poly = np.c_[np.linspace(0.0, 1.0, 50), np.zeros(50)]
    pts = np.c_[np.linspace(0.1, 0.9, 7), np.full(7, 0.003)]
    pts[3, 1] = -0.004  # the worst point, on the other side
    d_max, d_mean = ab.normal_distance(pts, poly)
    assert d_max == pytest.approx(0.004, abs=1e-12)
    assert d_mean == pytest.approx((6 * 0.003 + 0.004) / 7, abs=1e-12)  # unsigned mean


# --- classification and the boundary rule -----------------------------------------------------
@pytest.mark.parametrize(
    ("res", "tol", "boundary", "gated", "expected"),
    [(0.5, 1.0, False, True, "pass"), (1.0, 1.0, False, True, "pass"), (1.01, 1.0, False, True, "fail"),
     (9.0, 1.0, True, True, "boundary"), (9.0, 1.0, False, False, "not-gated"),
     (9.0, 1.0, True, False, "not-gated")],
)
def test_classify(res, tol, boundary, gated, expected):
    assert ab.classify(res, tol, boundary=boundary, gated=gated) == expected


def test_is_boundary_uses_half_an_element_in_metres():
    bounds = [0.2, 0.4]
    # perimeter 10 m, es 0.25 -> half element 0.125 m = 0.0125 in nd_arc
    assert ab.is_boundary(0.2 + 0.012, bounds, perimeter=10.0, element_size=0.25)
    assert not ab.is_boundary(0.2 + 0.013, bounds, perimeter=10.0, element_size=0.25)
    assert ab.is_boundary(0.4 - 0.0125, bounds, perimeter=10.0, element_size=0.25)


# --- independence and CLI contract -------------------------------------------------------------
def test_reader_imports_nothing_from_the_production_blade_path():
    tree = ast.parse((TOOLS / "audit_blade_input.py").read_text())
    top = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
    mods = [a.name for n in top if isinstance(n, ast.Import) for a in n.names]
    mods += [n.module for n in top if isinstance(n, ast.ImportFrom)]
    assert not any(m and m.startswith("aeroelast") for m in mods)
    # the reader class and its helpers never mention the production layer at all
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "YamlBlade":
            assert "aeroelast" not in ast.get_source_segment((TOOLS / "audit_blade_input.py").read_text(), node)


def test_cli_defaults_match_the_spec():
    args = ab.parse_args([])
    assert args.yaml == "tests/IEA-15-240-RWT.yaml"
    assert args.element_size == 0.25
    assert args.output_dir == "logs/check/input-audit"
    assert ab.STATIONS == (0.30, 0.60, 0.90, 0.97)


def test_bay_mean_thickness_and_deck_triangle_indices():
    yb = ab.YamlBlade.from_dict(mini_yaml())
    # cap thickness 0.02 -> 0.0 linearly over [0, 100] m: mean over [0, 50] is 0.015
    assert yb.bay_mean_thickness("cap", 0.0, 50.0) == pytest.approx(0.015, rel=1e-12)
    data = mini_yaml()
    vals = [[float(i) for i in range(21)], [2.0 * i for i in range(21)]]
    data["components"]["blade"]["elastic_properties_mb"] = {"six_x_six": {
        "stiff_matrix": {"grid": [0.0, 1.0], "values": vals},
        "inertia_matrix": {"grid": [0.0, 1.0], "values": vals}}}
    m11, k11, k33 = ab.YamlBlade.from_dict(data).deck_m11_k11_k33(0.5)
    assert (m11, k11, k33) == (pytest.approx(0.0 * 1.5), pytest.approx(0.0), pytest.approx(11 * 1.5))


# --- partial-span layers (IEA-3.4: caps and webs on 0.1-0.94, fillers on 0.1-0.8) ------------
def partial_span_yaml():
    """mini_yaml with the cap and the web skin defined on a span grid of 0.2-0.6 only."""
    data = mini_yaml()
    lin = lambda g, v: {"grid": list(g), "values": list(v)}  # noqa: E731
    layers = {x["name"]: x for x in data["components"]["blade"]["internal_structure_2d_fem"]["layers"]}
    layers["cap"]["thickness"] = lin([0.2, 0.6], [0.02, 0.01])
    layers["web0_skin"]["thickness"] = lin([0.2, 0.6], [0.003, 0.001])
    return data


def test_partial_span_layer_is_zero_outside_its_grid_and_linear_inside():
    yb = ab.YamlBlade.from_dict(partial_span_yaml())
    t = lambda r: yb.layer_state("cap", r)["thickness"]  # noqa: E731
    assert t(0.1) == 0.0 and t(0.7) == 0.0 and t(1.0) == 0.0
    assert t(0.2) == pytest.approx(0.02, abs=1e-15) and t(0.6) == pytest.approx(0.01, abs=1e-15)
    assert t(0.4) == pytest.approx(0.015, abs=1e-15)
    # the arcs of the layer keep their own (full-span) grid
    assert yb.layer_state("cap", 0.7)["start"] == pytest.approx(0.27, abs=1e-15)


def test_partial_span_layer_is_inactive_outside_its_grid():
    yb = ab.YamlBlade.from_dict(partial_span_yaml())
    names = lambda r, a: [x["name"] for x in yb.active_layers(r, a, 10.0)]  # noqa: E731
    assert names(0.1, 0.3) == ["skin"]
    assert names(0.4, 0.3) == ["skin", "cap"]
    assert names(0.7, 0.35) == ["skin"]
    assert 0.2 not in yb.all_bounds(0.1, 10.0)


def test_web_is_absent_outside_its_layer_grid():
    yb = ab.YamlBlade.from_dict(partial_span_yaml())
    assert yb.web_layers(0, 0.1) == [] and yb.web_layers(0, 0.7) == []
    assert yb.web_layers(0, 0.4)[0]["thickness"] == pytest.approx(0.002, abs=1e-15)


def test_bay_mean_thickness_counts_zero_outside_the_grid():
    yb = ab.YamlBlade.from_dict(partial_span_yaml())
    # bay 50-70 m: 0.0125 -> 0.01 linear over 50-60 m, zero over 60-70 m
    assert yb.bay_mean_thickness("cap", 50.0, 70.0) == pytest.approx(0.5 * 0.01125, rel=1e-2)
    assert yb.bay_mean_thickness("cap", 0.0, 15.0) == 0.0
