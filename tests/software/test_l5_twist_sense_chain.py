"""Pure helpers of ``tools/l5_twist_sense_chain.py`` (odd/tasks/l5-one-way.md T10): each against a case with a known answer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import l5_twist_sense_chain as ch  # noqa: E402


def airfoil_ring(chord, psi_deg, camber=0.04, n=80):
    """Cambered section in the mesh frame: LE at +x, suction (camber) side at +y, TE->LE direction at psi about +z."""
    s = np.linspace(0.0, 1.0, n)  # 0 = LE, 1 = TE
    thick = 0.12 * (0.2969 * np.sqrt(s) - 0.126 * s - 0.3516 * s**2 + 0.2843 * s**3 - 0.1036 * s**4)
    mid = camber * 4.0 * s * (1.0 - s)
    upper = np.column_stack([-s, mid + thick])  # suction side, x = -s (LE at x = 0, TE at x = -1)
    lower = np.column_stack([-s[::-1], (mid - thick)[::-1]])
    xy = np.vstack([upper, lower[1:]]) * chord
    a = np.radians(psi_deg)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return xy @ rot.T


@pytest.mark.parametrize("psi", [-14.0, 0.0, 3.0, 15.0])
def test_section_from_ring_recovers_le_psi_and_suction_side(psi):
    xy = airfoil_ring(4.0, psi) + np.array([0.3, -2.0])
    sec = ch.section_from_ring(xy)
    a = np.radians(psi)
    np.testing.assert_allclose(sec["psi_deg"], psi, atol=0.3)
    np.testing.assert_allclose(sec["chord"], 4.0, rtol=2e-2)
    # the LE is the +x end (rotated by psi), the suction normal is +y rotated by psi
    np.testing.assert_allclose(sec["suction_normal"], [-np.sin(a), np.cos(a)], atol=5e-2)
    assert (sec["le"] - sec["te"]) @ np.array([np.cos(a), np.sin(a)]) > 3.9


def test_handedness_of_the_mesh_triad_is_proper():
    le, suction, span = np.eye(3)  # +x LE, +y suction, +z span (the measured mesh triad)
    assert ch.handedness(le, suction, span) == 1
    assert ch.handedness(-le, suction, span) == -1  # one mirrored axis flips it


def test_blade_velocity_points_to_the_leading_edge_for_spin_about_plus_y():
    v = ch.blade_velocity(np.array([0.0, 0.8, 0.0]), np.array([0.0, 0.0, 50.0]))
    np.testing.assert_allclose(v, [40.0, 0.0, 0.0], atol=1e-12)


def test_angle_of_attack_is_phi_plus_psi_in_the_mesh_frame():
    # LE at +x, suction +y, blade moving toward +x at u, wind +y at w: alpha = phi + psi (psi: TE->LE angle about +z)
    u, w = 60.0, 10.0
    phi = np.degrees(np.arctan2(w, u))
    v_rel = ch.inflow(np.array([u, 0.0, 0.0]), np.array([0.0, w, 0.0]))
    for psi in (-12.0, 0.0, 12.0):
        a = np.radians(psi)
        c_hat = -np.array([np.cos(a), np.sin(a)])  # LE -> TE
        n_s = np.array([-np.sin(a), np.cos(a)])
        assert ch.angle_of_attack_deg(c_hat, n_s, v_rel[:2]) == pytest.approx(phi + psi, abs=1e-9)


def test_inflow_adds_induction_to_both_components():
    v = ch.inflow(np.array([60.0, 0.0, 0.0]), np.array([0.0, 10.0, 0.0]), a=0.2, ap=0.1)
    np.testing.assert_allclose(v, [-66.0, 8.0, 0.0], atol=1e-12)


def test_rotate_about_z_is_right_handed():
    np.testing.assert_allclose(ch.rotate_z(np.array([1.0, 0.0]), 90.0), [0.0, 1.0], atol=1e-15)


def test_ad_normal_and_tangential_follow_the_aerodyn_source_rotation():
    # AeroDyn_IO.f90:486-487: Fn = X ct - Y st, Ft = -X st - Y ct with BNFx = X, BNFy = -Y
    th = 0.2
    X, Yneg = 3.0, 1.0
    fn = X * np.cos(th) + Yneg * np.sin(th)
    ft = -X * np.sin(th) + Yneg * np.cos(th)
    np.testing.assert_allclose(ch.ad_to_rotor_plane(fn, ft, th), [X, Yneg], atol=1e-14)


def test_physical_ad_force_in_mesh_has_the_airfoil_normal_at_minus_theta():
    th = np.radians(10.0)
    n_only = ch.ad_force_in_mesh(1.0, 0.0, th, sense=-1)
    t_only = ch.ad_force_in_mesh(0.0, 1.0, th, sense=-1)
    np.testing.assert_allclose(n_only, [np.sin(th), np.cos(th)], atol=1e-15)  # +y rotated by -theta
    np.testing.assert_allclose(t_only, [np.cos(th), -np.sin(th)], atol=1e-15)  # +x rotated by -theta
    mirrored = ch.ad_force_in_mesh(1.0, 0.0, th, sense=+1)  # what tools/run_s5_oneway.apply_ad_loads builds
    np.testing.assert_allclose(mirrored, [-np.sin(th), np.cos(th)], atol=1e-15)


def test_ad_force_in_mesh_equals_rotor_plane_components():
    th = 0.3
    fx, fy = ch.ad_force_in_mesh(5.0, 2.0, th, sense=-1)
    x_n, x_t = ch.ad_to_rotor_plane(5.0, 2.0, th)
    assert fy == pytest.approx(x_n) and fx == pytest.approx(x_t)


def test_rotor_torque_about_plus_y_through_the_hub_centre():
    X = np.array([[0.0, 0.0, 10.0], [1.0, 0.0, 20.0]])
    f = np.array([[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    # arm z + hub = 10 + 4 = 14, force +x: torque about +y is +28
    assert ch.torque_about_y(X, f, hub_radius=4.0) == pytest.approx(28.0)
