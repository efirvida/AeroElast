"""Pure helpers of tools/l5_orientation_audit.py (L5 T12): BeamDyn key points, frame map, physical AoA sign."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import l5_orientation_audit as A  # noqa: E402
import twist_p1 as tp  # noqa: E402

BD_TEXT = """\
          1   member_total    - Total number of members (-)
          3   kp_total        - Total number of key points (-) [must be at least 3]
     1     3                 - Member number; Number of key points in this member
\t\t kp_xr \t\t\t kp_yr \t\t\t kp_zr \t\t initial_twist
\t\t  (m)  \t\t\t  (m)  \t\t\t  (m)  \t\t   (deg)
\t 0.00000e+00 \t 0.00000e+00 \t 0.00000e+00 \t 1.55946e+01
\t 1.00000e-01 \t 0.00000e+00 \t 2.00000e+00 \t 1.00000e+01
\t-4.00000e+00 \t 0.00000e+00 \t 1.17000e+02 \t-1.24000e+00
-------------------------------------------------------------------------------
"""


def test_parse_bd_keypoints_reads_the_key_point_table():
    kp = A.parse_bd_keypoints(BD_TEXT)
    assert kp.shape == (3, 4)
    np.testing.assert_allclose(kp[:, 0], [0.0, 0.1, -4.0])
    np.testing.assert_allclose(kp[:, 3], [15.5946, 10.0, -1.24])


def test_bd_frame_maps_flap_x_to_mesh_y_and_keeps_span():
    np.testing.assert_allclose(A.bd_to_mesh([1.0, 0.0, 0.0], tp.MESH_TO_BD), [0.0, 1.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(A.bd_to_mesh([0.0, 1.0, 0.0], tp.MESH_TO_BD), [-1.0, 0.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(A.bd_to_mesh([0.0, 0.0, 1.0], tp.MESH_TO_BD), [0.0, 0.0, 1.0], atol=1e-15)


@pytest.mark.parametrize("twist", [-6.0, 0.0, 2.5, 15.6])
def test_bd_chord_angle_is_minus_initial_twist_in_the_mesh_frame(twist):
    assert A.bd_chord_angle_deg(twist, tp.MESH_TO_BD) == pytest.approx(-twist, abs=1e-12)


def test_physical_aoa_flat_section_wind_from_minus_y_is_positive_atan_v_over_omega_r():
    # chord along +x (LE +x), suction +y (downwind), wind +y, blade moving +x: flow from the pressure side
    aoa = A.physical_aoa_deg(0.0, np.array([0.0, 1.0]), r=100.0, wind_dir_y=1.0, omega=0.8, wind=10.0)
    assert aoa == pytest.approx(np.degrees(np.arctan2(10.0, 80.0)), abs=1e-12)


def test_physical_aoa_drops_by_the_nose_up_rotation_toward_the_wind():
    # LE rotated toward upwind (-y) by 5 deg (psi = -5, positive feather twist): the AoA falls by 5 deg
    n = A.rotate2(np.array([0.0, 1.0]), -5.0)
    base = A.physical_aoa_deg(0.0, np.array([0.0, 1.0]), 100.0, 1.0, 0.8, wind=10.0)
    assert A.physical_aoa_deg(-5.0, n, 100.0, 1.0, 0.8, wind=10.0) == pytest.approx(base - 5.0, abs=1e-12)


def test_response_angle_is_measured_from_flap_toward_edge():
    np.testing.assert_allclose(A.response_angle_deg([0.0, 1.0], [1.0, 1.0]), [0.0, 45.0])
