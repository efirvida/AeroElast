"""Pure parts of ``tools/l5_twist_torque_sensitivity.py`` (L5 T8, S2): twist interpolation and the sign conversion.

No BEM evaluation here: the convention is asserted through the tool's constant and the way it is applied, not through a torque number.
"""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[2] / "tools" / "l5_twist_torque_sensitivity.py"
_spec = importlib.util.spec_from_file_location("l5_twist_torque_sensitivity", _PATH)
tool = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(tool)


def test_interpolation_is_exact_at_nodes():
    z = np.array([60.0, 90.0, 117.0])
    d = np.array([0.0, -0.5, -1.0])
    assert np.array_equal(tool.interpolate_dtwist(z, z, d), d)


def test_interpolation_is_zero_inboard_of_first_z_and_linear_between():
    z = np.array([60.0, 90.0, 117.0])
    d = np.array([0.3, -0.5, -1.0])
    out = tool.interpolate_dtwist(np.array([0.0, 30.0, 59.999, 75.0]), z, d)
    assert np.all(out[:3] == 0.0)
    assert out[3] == pytest.approx(0.5 * (0.3 - 0.5))


def test_interpolation_rejects_unsorted_or_mismatched_input():
    with pytest.raises(ValueError):
        tool.interpolate_dtwist(np.array([1.0]), np.array([2.0, 1.0]), np.array([0.0, 0.0]))
    with pytest.raises(ValueError):
        tool.interpolate_dtwist(np.array([1.0]), np.array([1.0, 2.0]), np.array([0.0]))


def test_sign_convention_nose_up_mesh_twist_lowers_the_bem_twist():
    # fsi_participant.py:494-503: a positive mesh rotation about +span (nose-up, R1 of twist_p1) is opposite-sensed to the CCBlade
    # twist theta (alpha = phi - theta), so it must DECREASE theta, i.e. raise the angle of attack toward stall.
    assert tool.MESH_TO_BEM_SENSE == -1.0
    nose_up = np.deg2rad([0.5, 1.0])
    dtheta = tool.mesh_twist_to_bem(nose_up)
    assert np.all(dtheta < 0.0)
    assert np.array_equal(tool.mesh_twist_to_bem(-nose_up), -dtheta)  # nose-down raises theta, toward feather
    assert np.allclose(dtheta, -nose_up)


def test_bem_twist_is_reference_plus_signed_elastic_twist():
    ref = np.deg2rad(np.array([5.0, 3.0, 1.0]))
    z = np.array([10.0, 20.0, 30.0])
    out = tool.bem_twist(ref, z, z[1:], np.array([1.0, 2.0]))  # deg, mesh sense, zero at z=10
    assert out[0] == ref[0]
    assert out[1:] == pytest.approx(ref[1:] - np.deg2rad([1.0, 2.0]))


def test_table_extraction_picks_the_estimator(tmp_path):
    doc = {"table": {"z": [1.0, 2.0], "shell": {"twist_skin_deg": [0.1, 0.2], "twist_eq_deg": [0.3, 0.4]},
                     "beamdyn_x0.01": {"twist_deg": [0.5, 0.6]}}}
    f = tmp_path / "t.json"
    f.write_text(json.dumps(doc))
    z, shell, bd = tool.read_parity_json(f, "skin")
    assert (z.tolist(), shell.tolist(), bd.tolist()) == ([1.0, 2.0], [0.1, 0.2], [0.5, 0.6])
    assert tool.read_parity_json(f, "eq")[1].tolist() == [0.3, 0.4]
    with pytest.raises(ValueError):
        tool.read_parity_json(f, "nope")


def test_relative_deltas_are_against_rigid():
    cases = {"rigid": {"thrust": 100.0, "torque": 10.0, "power": 1.0},
             "beamdyn": {"thrust": 101.0, "torque": 10.2, "power": 1.02},
             "shell": {"thrust": 103.0, "torque": 10.5, "power": 1.05}}
    d = tool.delta_shell_minus_beamdyn(cases)
    assert d["torque"]["abs"] == pytest.approx(0.3)
    assert d["torque"]["pct_of_rigid"] == pytest.approx(3.0)
    assert d["thrust"]["pct_of_rigid"] == pytest.approx(2.0)
