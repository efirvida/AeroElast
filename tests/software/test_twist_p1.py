"""Pure functions of ``tools/twist_p1.py`` (P1 re-pre-registration R1-R3, odd/tasks/twist-discrepancy-closure.md).

They are the arithmetic behind the controls: Wiener-Milenkovic (WM) component to angle (R2), the
swing-twist split, the mesh <-> BeamDyn frame map (R1), the 6x6 compliance integration of C1 and
the rotating root-axial-force reference of C2a.  Each is tested against a case whose answer is
known by construction, before the harness is used.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import twist_p1 as tp  # noqa: E402


def rot(axis, angle):
    """Rodrigues rotation matrix about ``axis`` (right-handed)."""
    a = np.asarray(axis, float)
    a = a / np.linalg.norm(a)
    k = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


# --- WM component -> angle (R2 primary) and WM -> matrix ---------------------------------------


@pytest.mark.parametrize("angle", [-0.5, -0.05, 0.0, 0.01, 0.3, 1.2])
def test_wm_z_component_to_angle(angle):
    c = 4.0 * np.tan(angle / 4.0)  # WM component of a pure z rotation
    assert tp.wm_z_angle(c) == pytest.approx(angle, rel=1e-12, abs=1e-14)


def test_wm_z_angle_is_not_the_small_angle_component():
    # c_z = 4 tan(phi/4) ~ phi only to first order: the primary estimator must use 4 atan(c/4).
    angle = 0.8
    c = 4.0 * np.tan(angle / 4.0)
    assert abs(c - angle) > 1e-2
    assert tp.wm_z_angle(c) == pytest.approx(angle, rel=1e-12)


def test_wm_to_matrix_matches_rodrigues():
    axis, angle = np.array([0.3, -0.5, 0.8]), 0.7
    c = 4.0 * np.tan(angle / 4.0) * axis / np.linalg.norm(axis)
    assert np.allclose(tp.wm_to_matrix(c), rot(axis, angle), atol=1e-13)
    assert np.allclose(tp.wm_to_matrix(np.zeros(3)), np.eye(3))


def test_wm_to_matrix_is_proper_rotation():
    R = tp.wm_to_matrix([0.2, -0.1, 0.35])
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-13)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-13)


# --- swing-twist split (R2 secondary) ----------------------------------------------------------


def test_swing_twist_pure_z_rotation():
    twist, swing = tp.swing_twist(rot([0, 0, 1], 0.4), [0, 0, 1])
    assert twist == pytest.approx(0.4, abs=1e-13)
    assert swing == pytest.approx(0.0, abs=1e-13)


def test_swing_twist_pure_swing_has_zero_twist():
    # the minimal rotation taking the tangent to a tilted one carries no twist by definition
    t0 = np.array([0.0, 0.0, 1.0])
    R = tp.minimal_rotation(t0, [0.2, -0.1, 1.0])
    twist, swing = tp.swing_twist(R, t0)
    assert twist == pytest.approx(0.0, abs=1e-13)
    assert swing == pytest.approx(np.arctan(np.hypot(0.2, 0.1)), rel=1e-12)


def test_swing_twist_recovers_composed_twist():
    t0 = np.array([0.0, 0.0, 1.0])
    S = tp.minimal_rotation(t0, [0.15, 0.05, 1.0])
    t_def = S @ t0
    for tw in (-0.3, 0.02, 0.5):
        R = rot(t_def, tw) @ S  # R = R_twist(deformed tangent) . R_swing
        twist, swing = tp.swing_twist(R, t0)
        assert twist == pytest.approx(tw, abs=1e-12)
        assert swing == pytest.approx(np.arctan(np.hypot(0.15, 0.05)), rel=1e-12)


def test_twist_and_wm_z_agree_to_second_order_in_the_flap_slope():
    # a flap slope tilts the tangent: the z-component estimator and the tangent twist differ at O(slope^2)
    t0 = np.array([0.0, 0.0, 1.0])
    S = rot([0, 1, 0], 0.05)
    R = rot(S @ t0, 0.02) @ S
    twist, _ = tp.swing_twist(R, t0)
    q = tp.rot_to_quat(R)
    phi_z = 2.0 * np.arctan2(q[3], q[0])  # the z-only estimator
    assert twist == pytest.approx(0.02, abs=1e-12)
    assert abs(phi_z - twist) < 2e-3


# --- frame map (R1) ----------------------------------------------------------------------------


def test_frame_map_is_a_proper_rotation():
    Q = tp.MESH_TO_BD
    assert np.linalg.det(Q) == pytest.approx(1.0)
    assert np.allclose(Q @ Q.T, np.eye(3))


def test_frame_map_axes():
    # x_BD = mesh +y, y_BD = -mesh x, z_BD = z: components of a mesh vector in the BD basis
    assert np.allclose(tp.MESH_TO_BD @ [0, 1, 0], [1, 0, 0])
    assert np.allclose(tp.MESH_TO_BD @ [1, 0, 0], [0, -1, 0])
    assert np.allclose(tp.MESH_TO_BD @ [0, 0, 1], [0, 0, 1])


@pytest.mark.parametrize("angle", [-0.3, 0.2])
def test_frame_map_preserves_the_z_rotation_sense(angle):
    R_mesh = rot([0, 0, 1], angle)
    R_bd = tp.MESH_TO_BD @ R_mesh @ tp.MESH_TO_BD.T
    assert np.allclose(R_bd, rot([0, 0, 1], angle), atol=1e-14)
    assert tp.wm_z_angle(4 * np.tan(angle / 4)) == pytest.approx(angle)


# --- 6x6 compliance integration (C1) -----------------------------------------------------------


def _prismatic(n, GJ):
    K = np.tile(np.diag([1e9, 1e9, 5e9, 1e8, 1e8, GJ]), (n, 1, 1))
    return np.linspace(0.0, 40.0, n), K, np.tile(np.eye(3), (n, 1, 1))


def test_compliance_integration_prismatic_isotropic_is_TL_over_GJ():
    GJ, T = 3.0e8, 1.0e4
    arc, K, R = _prismatic(21, GJ)
    twist = tp.integrate_twist(arc, K, R, [0, 0, T])
    assert twist[0] == 0.0
    assert twist[-1] == pytest.approx(T * 40.0 / GJ, rel=1e-12)
    assert twist[10] == pytest.approx(T * 20.0 / GJ, rel=1e-12)


def test_compliance_integration_sees_the_k16_coupling():
    GJ, T = 3.0e8, 1.0e4
    arc, K, R = _prismatic(21, GJ)
    K[:, 0, 5] = K[:, 5, 0] = 4.0e7  # K16 = K61, positive definite still
    twist = tp.integrate_twist(arc, K, R, [0, 0, T])
    naive = T * 40.0 / GJ
    exact = T * np.linalg.inv(K[0])[5, 5] * 40.0
    assert twist[-1] == pytest.approx(exact, rel=1e-12)
    assert abs(twist[-1] / naive - 1.0) > 1e-3


def test_compliance_integration_transports_the_couple_to_the_local_frame():
    # a section frame tilted by 90 deg about x: the global z couple is local +y (Rx(90) maps local y onto z), hence resisted by K55
    arc = np.linspace(0.0, 10.0, 5)
    K = np.tile(np.diag([1e9, 1e9, 1e9, 2e8, 4e8, 8e8]), (5, 1, 1))
    R = np.tile(rot([1, 0, 0], np.pi / 2), (5, 1, 1))
    local = tp.transport_moment(R, [0, 0, 1.0e4])
    assert np.allclose(local[0], [0.0, 1.0e4, 0.0], atol=1e-9)
    assert tp.integrate_twist(arc, K, R, [0, 0, 1.0e4])[-1] == pytest.approx(0.0, abs=1e-15)


def test_section_frames_untwisted_straight_blade_is_the_identity():
    kp = np.column_stack([np.zeros(5), np.zeros(5), np.linspace(0, 8, 5)])
    R = tp.section_frames(kp, np.zeros(5))
    assert np.allclose(R, np.eye(3))


# --- rotating root axial force (C2a) -----------------------------------------------------------


def test_rotating_root_force_uniform_straight_beam_closed_form():
    m, L, Rh, Om = 3127.4, 117.0, 3.97, 0.791681
    s = np.linspace(0.0, L, 401)
    F = tp.centrifugal_root_force(s, np.full_like(s, m), Rh + s, Om)
    assert F == pytest.approx(m * Om**2 * (L * Rh + L**2 / 2), rel=1e-12)


def test_rotor_axis_radius_is_the_distance_from_the_x_axis():
    kp = np.array([[0.1, 0.0, 0.0], [-0.5, 0.0, 10.0], [-4.0, 3.0, 20.0]])
    rho = tp.rotor_axis_radius(kp, 3.97)
    assert np.allclose(rho, [3.97, 13.97, np.hypot(3.0, 23.97)])


def test_relative_error_and_pass_flag():
    assert tp.within(1.005, 1.0, 0.01)
    assert not tp.within(1.02, 1.0, 0.01)
    assert not tp.within(-1.0, 1.0, 0.01)  # opposite sign never passes
