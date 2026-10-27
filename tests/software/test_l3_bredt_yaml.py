"""``tools/l3_bredt_yaml.py``: thin-wall multi-cell Bredt GJ of a windIO v1 yaml section (odd/tasks/l3-closure.md, T9).

The Bredt core is checked against closed forms (circular tube, two-cell rectangle), the 6x6 shear
condensation on a hand matrix, and the partial-span rule on a tiny yaml dict.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np
import pytest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import l3_bredt_yaml as lb  # noqa: E402
from audit_blade_input import YamlBlade  # noqa: E402
from test_audit_blade_input import mini_yaml  # noqa: E402


def test_single_cell_circular_tube_matches_2_pi_R3_Gt():
    R, gt, n = 0.7, 3.0e6, 400
    th = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    pts = R * np.c_[np.cos(th), np.sin(th)]
    out = lb.multicell_bredt(pts, np.full(n, gt), np.zeros(n, int), [])
    assert out["gj"] == pytest.approx(2 * np.pi * R**3 * gt, rel=5e-3)
    assert out["n_cells"] == 1


def test_two_cell_rectangle_matches_closed_form():
    a1, a2, h = 0.6, 0.9, 0.4
    L = a1 + a2
    pts = np.array([[L, 0], [a1, 0], [0, 0], [0, h], [a1, h], [L, h]], float)
    g = [2.0e6, 3.0e6, 5.0e6, 7.0e6, 1.1e6, 1.7e6]  # Gt of segments v_i -> v_i+1
    cells = [0, 1, 1, 1, 0, 0]
    gw = 4.0e6
    out = lb.multicell_bredt(pts, np.array(g), np.array(cells), [{"i_s": 1, "i_p": 4, "gt": gw}])
    # cell 0 = right (walls 0, 4, 5 + web), cell 1 = left (walls 1, 2, 3 + web)
    f11 = (a2 / g[0] + a2 / g[4] + h / g[5]) + h / gw
    f22 = (a1 / g[1] + h / g[2] + a1 / g[3]) + h / gw
    f12 = -h / gw
    b1, b2 = 2 * a2 * h, 2 * a1 * h
    det = f11 * f22 - f12**2
    q1 = (b1 * f22 - f12 * b2) / det
    q2 = (f11 * b2 - f12 * b1) / det
    assert out["gj"] == pytest.approx(b1 * q1 + b2 * q2, rel=1e-12)
    assert out["areas"] == pytest.approx([a2 * h, a1 * h], rel=1e-12)


def test_shear_condensation_on_a_hand_matrix():
    k = np.diag([4.0, 2.0, 9.0, 9.0, 9.0, 10.0])
    k[0, 5] = k[5, 0] = 2.0
    k[1, 5] = k[5, 1] = 1.0
    k[0, 1] = k[1, 0] = 0.5
    s = np.array([[4.0, 0.5], [0.5, 2.0]])
    c = np.array([2.0, 1.0])
    assert lb.condense_shear(k) == pytest.approx(10.0 - c @ np.linalg.solve(s, c), rel=1e-14)


def test_beamdyn_interpolation_is_linear_in_eta():
    etas = np.array([0.0, 1.0])
    ks = np.array([np.eye(6), 3 * np.eye(6)])
    assert lb.interp_deck(etas, ks, 0.25) == pytest.approx(1.5 * np.eye(6))


def test_web_absent_outside_its_layer_grid_gives_one_cell():
    d = copy.deepcopy(mini_yaml())
    for layer in d["components"]["blade"]["internal_structure_2d_fem"]["layers"]:
        if layer["name"] == "web0_skin":
            layer["thickness"] = {"grid": [0.0, 0.5], "values": [0.003, 0.001]}
    yb = YamlBlade.from_dict(d)
    inside = lb.yaml_section(yb, 0.25, n=200, offset=False)
    outside = lb.yaml_section(yb, 0.8, n=200, offset=False)
    assert inside["n_cells"] == 2 and outside["n_cells"] == 1
    assert outside["gj"] > 0.0
