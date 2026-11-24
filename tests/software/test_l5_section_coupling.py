"""Pure helpers of ``tools/l5_section_coupling.py`` (odd/tasks/l5-one-way.md T9): each against a case with a known answer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import l5_section_coupling as sc  # noqa: E402
import twist_p1 as tp  # noqa: E402


def rz(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def test_wm_to_rotvec_matches_rotation_matrix_and_small_angle():
    n = np.array([1.0, 2.0, -2.0]) / 3.0
    phi = 0.7
    c = 4.0 * np.tan(phi / 4.0) * n
    np.testing.assert_allclose(sc.wm_to_rotvec(c), phi * n, atol=1e-14)
    np.testing.assert_allclose(sc.wm_to_rotvec([[0.0, 0.0, 0.0], c])[1], phi * n, atol=1e-14)
    np.testing.assert_array_equal(sc.wm_to_rotvec([0.0, 0.0, 0.0]), np.zeros(3))
    np.testing.assert_allclose(sc.wm_to_rotvec([1e-6, 0.0, 0.0]), [1e-6, 0.0, 0.0], rtol=1e-9)


def test_frame_map_is_the_twist_p1_one_and_invertible():
    v = np.array([0.3, -1.2, 0.5])
    # x_BD = mesh +y, y_BD = -mesh x: a mesh +x vector is -y in BD
    np.testing.assert_allclose(sc.mesh_to_bd([1.0, 0.0, 0.0]), [0.0, -1.0, 0.0])
    np.testing.assert_allclose(sc.bd_to_mesh([0.0, 1.0, 0.0]), [-1.0, 0.0, 0.0])
    np.testing.assert_allclose(sc.bd_to_mesh(sc.mesh_to_bd(v)), v, atol=1e-15)
    np.testing.assert_allclose(sc.mesh_to_bd(v), tp.MESH_TO_BD @ v, atol=1e-15)


def test_tip_couple_forces_have_the_requested_moment_and_no_net_force():
    rng = np.random.default_rng(3)
    xy = rng.normal(size=(40, 2)) * [1.5, 0.2] + [0.3, -0.1]
    for mx, my in ((1.0e4, 0.0), (0.0, 2.5e3), (1.0e4, -7.0e3)):
        f = sc.tip_couple_forces(xy, mx, my)
        assert abs(f.sum()) < 1e-9 * max(abs(mx), abs(my))
        assert (xy[:, 1] * f).sum() == pytest.approx(mx, abs=1e-8)
        assert (-xy[:, 0] * f).sum() == pytest.approx(my, abs=1e-8)


def test_fit_rotation_curvature_recovers_a_synthetic_field():
    rng = np.random.default_rng(1)
    X = np.column_stack([rng.uniform(-1.5, 1.5, 600), rng.uniform(-0.3, 0.3, 600), rng.uniform(40.0, 50.0, 600)])
    z0, th, ka = 45.0, np.array([2e-3, -1e-3]), np.array([3e-4, 5e-5])
    dz = X[:, 2] - z0
    uz = 7e-3 + 2e-4 * dz + (th[0] + ka[0] * dz) * X[:, 1] - (th[1] + ka[1] * dz) * X[:, 0]
    fit = sc.fit_rotation_curvature(X, uz, z0, 5.0)
    np.testing.assert_allclose(fit["kappa"], ka, rtol=1e-10)
    # the fit centres x, y on the slab centroid: that shift moves only c and e, never the rotation or its rate
    np.testing.assert_allclose(fit["theta"], th, rtol=1e-9)
    assert fit["rms"] < 1e-12 and fit["n"] == 600


def test_fit_rotation_curvature_window_and_too_few_nodes():
    X = np.column_stack([np.linspace(-1, 1, 30), np.linspace(-0.2, 0.2, 30), np.linspace(0.0, 30.0, 30)])
    with pytest.raises(ValueError, match="nodes"):
        sc.fit_rotation_curvature(X, np.zeros(30), 15.0, 1.0)


def test_slope_fit_is_exact_for_a_linear_field_and_centred_for_a_quadratic():
    z = np.linspace(0.0, 20.0, 41)
    y = np.column_stack([2.0 + 3.0 * z, -1.0 * z])
    np.testing.assert_allclose(sc.slope_fit(z, y, 10.0, 3.0), [3.0, -1.0], atol=1e-12)
    q = 5.0 + 0.5 * (z - 10.0) ** 2  # symmetric window: slope at the centre is 0
    assert sc.slope_fit(z, q, 10.0, 3.0) == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError, match="samples"):
        sc.slope_fit(z, q, 10.0, 0.2)


def test_compliance_and_principal_angle_of_a_rotated_section():
    ee, ff = 1.0e9, 1.0e8  # edge stiff, flap soft
    for psi in (0.0, np.radians(-2.0), np.radians(20.0), np.radians(-60.0)):
        # soft (flap) axis at angle psi from +x: K = R diag(ff, ee) R^T, C = K^-1
        Rm = np.array([[np.cos(psi), -np.sin(psi)], [np.sin(psi), np.cos(psi)]])
        K = Rm @ np.diag([ff, ee]) @ Rm.T
        C = np.linalg.inv(K)
        assert sc.principal_angle_deg(C) == pytest.approx(np.degrees(psi), abs=1e-9)
        # responses to Mx, My reproduce C
        M = 1e4
        k_mx, k_my = C @ [M, 0.0], C @ [0.0, M]
        np.testing.assert_allclose(sc.compliance_from_responses(k_mx, k_my, M), C, rtol=1e-12)
    row = sc.table_row(np.linalg.inv(K))
    assert row["C_fe"] == row["C"][0][1] and row["ratio_fe_ff"] == pytest.approx(row["C_fe"] / row["C_ff"])
    assert sc.principal_angle_deg(np.eye(2) * 2.0 + [[0, 1e-3], [1e-3, 0]]) in (45.0, -45.0)


def test_deck_bending_compliance_frames_and_signs():
    # decoupled section K: flap-like K44 (about section x_BD) soft, K55 stiff; no twist, straight: mesh x ~ -y_BD
    K = np.diag([1e9, 1e9, 5e10, 1e8, 1e9, 4e9])
    C = sc.deck_bending_compliance(K, np.eye(3))
    # section x_BD = mesh y: curvature about x_BD (soft, 1/1e8) is about mesh y; mesh x = -y_BD (stiff)
    np.testing.assert_allclose(C, np.diag([1.0 / 1e9, 1.0 / 1e8]), rtol=1e-12)
    assert sc.principal_angle_deg(C) == pytest.approx(90.0)
    # a coupling term K45 and the force-bending condensation do enter through the full inverse
    K2 = K.copy()
    K2[2, 3] = K2[3, 2] = 1e9  # extension-bending coupling K34: pure couple condenses it
    S = np.linalg.inv(K2)[3:, 3:]
    np.testing.assert_allclose(sc.deck_bending_compliance(K2, np.eye(3))[1, 1], S[0, 0], rtol=1e-12)
    assert S[0, 0] == pytest.approx(1.0 / (1e8 - 1e18 / 5e10))  # S44 = 1/(K44 - K34^2/K33): softer than 1/K44
    # per-z conversion divides by t_z
    np.testing.assert_allclose(sc.deck_bending_compliance(K, np.eye(3), 0.5), 2.0 * C, rtol=1e-14)


def test_deck_frame_rotation_sign_matters_and_twist_sense_is_explicit():
    # a section rotated by +a about z_BD: its soft axis (x_BD of the section) points at +a from x_BD
    K = np.diag([1e9, 1e9, 5e10, 1e8, 1e9, 4e9])
    a = np.radians(10.0)
    Cp, Cm = (sc.deck_bending_compliance(K, rz(s * a)) for s in (+1.0, -1.0))
    assert np.sign(Cp[0, 1]) == -np.sign(Cm[0, 1]) != 0
    assert Cp[0, 0] == pytest.approx(Cm[0, 0]) and Cp[1, 1] == pytest.approx(Cm[1, 1])
    # soft axis x_BD = mesh y rotated by +a about z (BD and mesh share z and handedness): mesh angle 90 + a
    assert sc.principal_angle_deg(Cp) == pytest.approx(90.0 + 10.0 - 180.0, abs=1e-9)


def test_thin_wall_block_two_flanges_and_cross_term():
    # two equal flanges at (+-1, 0): Iyy = 2 w, Ixx = 0 (no wall own bending), no cross term
    s = np.array([[1.0, 1.0, 0.0], [1.0, -1.0, 0.0]])
    np.testing.assert_allclose(sc.thin_wall_block(s), [[0.0, 0.0], [0.0, 2.0]], atol=1e-14)
    # rotating the section by 30 deg about z gives a cross term -Ixy and keeps the invariants
    a = np.radians(30.0)
    A = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    B = sc.thin_wall_block(s, A)
    assert B[0, 1] == pytest.approx(-(2.0 * np.cos(a) * np.sin(a)))
    assert np.trace(B) == pytest.approx(2.0) and np.linalg.det(B) == pytest.approx(0.0, abs=1e-14)


def test_alignment_recovers_rotation_translation_and_reflection():
    th = np.linspace(0.0, 2.0 * np.pi, 1200, endpoint=False)
    contour = np.column_stack([1.3 * np.cos(th), 0.2 * np.sin(th) + 0.08 * np.cos(th) ** 2])  # thin cambered ellipse
    a = np.radians(-7.0)
    A = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    for refl in (False, True):
        T = A @ (np.diag([1.0, -1.0]) if refl else np.eye(2))
        ring = (contour @ T.T + [4.0, -2.0])[::7]
        best = sc.align_candidates(contour, ring)[0]
        np.testing.assert_allclose(best["A"], T, atol=2e-3)
        assert best["residual"] < 2e-3
    assert sc.chord_angle_deg(ring) == pytest.approx(-7.0, abs=0.3)


def test_verdict_applies_the_registered_tolerances_without_widening():
    ref = [{"z": 1.0, "C_ff": 1.0, "C_ee": 0.1, "psi_deg": 0.0}]
    ok = [{"z": 1.0, "C_ff": 1.04, "C_ee": 0.096, "psi_deg": 0.9}]
    bad_c = [{"z": 1.0, "C_ff": 1.06, "C_ee": 0.1, "psi_deg": 0.0}]
    bad_a = [{"z": 1.0, "C_ff": 1.0, "C_ee": 0.1, "psi_deg": -1.2}]
    assert sc.verdict(ok, ref)["all_ok"]
    assert not sc.verdict(bad_c, ref)["all_ok"]
    assert not sc.verdict(bad_a, ref)["all_ok"]


def test_tip_deflections_match_the_cantilever_closed_forms():
    L, EI, q = 100.0, 2.0e9, 50.0
    z = np.linspace(0.0, L, 4001)
    C = np.broadcast_to(np.diag([1.0 / EI, 1.0 / (3.0 * EI)]), (len(z), 2, 2))
    one, zero = np.full_like(z, q), np.zeros_like(z)
    ux, uy = sc.tip_deflections(z, C, zero, one)  # flap load f_y: u_y = q L^4 / (8 EI_flap), no edge response
    assert uy == pytest.approx(q * L**4 / (8 * EI), rel=1e-5) and ux == pytest.approx(0.0, abs=1e-12)
    ux, uy = sc.tip_deflections(z, C, one, zero)  # edge load f_x: u_x = q L^4 / (8 EI_edge), EI_edge = 3 EI
    assert ux == pytest.approx(q * L**4 / (8 * 3.0 * EI), rel=1e-5) and uy == pytest.approx(0.0, abs=1e-12)
    Cx = C.copy()
    Cx[:, 0, 1] = Cx[:, 1, 0] = 0.1 / EI  # a cross term couples them: kappa_y = C_ef Mx
    ux, _ = sc.tip_deflections(z, Cx, zero, one)
    assert ux == pytest.approx(-0.1 * q * L**4 / (8 * EI), rel=1e-5)  # Mx = -q (L - z)^2 / 2 -> negative edge deflection
