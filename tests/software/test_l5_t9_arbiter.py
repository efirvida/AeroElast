"""Pure helpers of tools/l5_t9_arbiter.py (L5 T9): section compliance, principal angle, frame rotation, kinematic fit."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import l5_t9_arbiter as A  # noqa: E402
import twist_p1 as tp  # noqa: E402


def test_compliance_from_responses_is_response_per_unit_load():
    C0 = np.array([[2.0, 0.3, 0.0], [0.3, 1.0, 0.1], [0.0, 0.1, 5.0]]) * 1e-9
    loads = np.array([1e6, 2e6, 4e6])
    resp = [C0[:, j] * loads[j] for j in range(3)]  # response vector under load j (linear)
    np.testing.assert_allclose(A.compliance_from_responses(resp, loads), C0, rtol=1e-12)


def test_rot_z3_is_right_handed_about_plus_z():
    np.testing.assert_allclose(A.rot_z3(90.0) @ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(A.rot_z3(90.0) @ [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], atol=1e-15)


def test_rotate_compliance_swaps_axes_at_90_deg_and_keeps_torsion():
    C = np.diag([2.0, 1.0, 5.0])
    out = A.rotate_compliance(C, A.rot_z3(90.0))
    np.testing.assert_allclose(out, np.diag([1.0, 2.0, 5.0]), atol=1e-15)


def test_principal_angle_follows_the_rotation_and_the_cross_term_sign():
    C = np.diag([2.0, 1.0, 5.0])
    assert A.principal_angle_deg(C) == pytest.approx(0.0, abs=1e-12)
    assert A.principal_angle_deg(A.rotate_compliance(C, A.rot_z3(30.0))) == pytest.approx(30.0, abs=1e-9)
    assert A.principal_angle_deg(A.rotate_compliance(C, A.rot_z3(100.0))) == pytest.approx(-80.0, abs=1e-9)
    assert A.principal_angle_deg(np.diag([1.0, 2.0, 5.0])) == pytest.approx(90.0)
    cross = np.array([[2.0, 0.1], [0.1, 1.0]])  # analytic: 0.5 atan2(2 c12, c11 - c22)
    assert A.principal_angle_deg(cross) == pytest.approx(0.5 * np.degrees(np.arctan2(0.2, 1.0)), abs=1e-9)
    assert A.principal_angle_deg(np.array([[2.0, -0.1], [-0.1, 1.0]])) < 0.0


def test_section_metrics_of_a_diagonal_and_of_a_rotated_section():
    EIf, EIe, GJ = 4.0e9, 9.0e9, 1.0e8
    C = np.diag([1 / EIf, 1 / EIe, 1 / GJ])
    m = A.section_metrics(C)
    assert m["EI_flap"] == pytest.approx(EIf) and m["EI_edge"] == pytest.approx(EIe) and m["GJ"] == pytest.approx(GJ)
    assert m["ratio_fe_ff"] == 0.0 and m["psi_deg"] == pytest.approx(0.0, abs=1e-12)
    Cr = A.rotate_compliance(C, A.rot_z3(10.0))
    mr = A.section_metrics(Cr)
    assert mr["psi_deg"] == pytest.approx(10.0, abs=1e-9)
    assert mr["EI_min"] == pytest.approx(EIf) and mr["EI_max"] == pytest.approx(EIe)  # principal values are frame invariant
    assert mr["GJ"] == pytest.approx(GJ)
    assert mr["ratio_fe_ff"] == pytest.approx(Cr[0, 1] / Cr[0, 0]) and mr["ratio_ef_ee"] == pytest.approx(Cr[1, 0] / Cr[1, 1])
    assert abs(mr["ratio_fe_ff"]) > 0.0


def test_deck_compliance_mesh_maps_axes_and_divides_by_tz():
    K = np.diag([1.0, 1.0, 1.0, 3.0e9, 5.0e9, 2.0e8])
    C = A.deck_compliance_mesh(K, np.eye(3), tz=1.0)
    # BeamDyn x_BD = mesh +y, y_BD = -mesh x: mesh flap (about x) is BeamDyn y
    np.testing.assert_allclose(np.diag(C), [1 / 5.0e9, 1 / 3.0e9, 1 / 2.0e8], rtol=1e-12)
    np.testing.assert_allclose(A.deck_compliance_mesh(K, np.eye(3), tz=0.5), 2.0 * C, rtol=1e-12)


def test_deck_compliance_mesh_cross_term_sign_and_frame_rotation():
    K = np.diag([1.0, 1.0, 1.0, 3.0e9, 5.0e9, 2.0e8])
    K[3, 4] = K[4, 3] = 1.0e9
    S = np.linalg.inv(K)[3:, 3:]
    C = A.deck_compliance_mesh(K, np.eye(3), 1.0)
    P = tp.MESH_TO_BD  # v_mesh = v_bd @ P
    # explicit vector algebra: curvature (mesh) under a mesh couple e_j
    for j in range(3):
        m_bd = P @ np.eye(3)[j]  # mesh -> BD components
        k_mesh = (S @ m_bd) @ P
        np.testing.assert_allclose(C[:, j], k_mesh, rtol=1e-12, atol=1e-30)
    # section rotated by Rz(-tw): the 6x6 is expressed in the section frame, R maps it into the blade frame
    tw = 7.0
    Cr = A.deck_compliance_mesh(K, A.rot_z3(-tw), 1.0)
    Cb = A.rot_z3(-tw) @ S @ A.rot_z3(-tw).T
    np.testing.assert_allclose(Cr, P.T @ Cb @ P, rtol=1e-12)
    # a z rotation in BD is the same rotation about mesh z: the principal angle moves by the rotation of the MESH axes
    Kd = np.diag([1.0, 1.0, 1.0, 3.0e9, 5.0e9, 2.0e8])
    a0 = A.principal_angle_deg(A.deck_compliance_mesh(Kd, np.eye(3), 1.0)[:2, :2])
    a1 = A.principal_angle_deg(A.deck_compliance_mesh(Kd, A.rot_z3(-tw), 1.0)[:2, :2])
    assert (a1 - a0 + 90.0) % 180.0 - 90.0 == pytest.approx(-tw, abs=1e-9)


def _synthetic_field(X, kappa, twist_rate, z0):
    """Displacements of a section moving rigidly: translation cubic in dz, rotations quadratic, known rates at dz = 0."""
    dz = X[:, 2] - z0
    xc, yc = X[:, 0] - X[:, 0].mean(), X[:, 1] - X[:, 1].mean()
    thx = 1e-3 + kappa[0] * dz + 0.2 * kappa[0] * dz**2
    thy = -2e-3 + kappa[1] * dz - 0.1 * kappa[1] * dz**2
    thz = 5e-4 + twist_rate * dz + 0.05 * twist_rate * dz**2
    ux = 1e-2 + 2e-3 * dz + 1e-4 * dz**2 + 1e-5 * dz**3 - thz * yc
    uy = -3e-2 + 1e-3 * dz - 2e-4 * dz**2 + thz * xc
    uz = 4e-3 + 1e-3 * dz + thx * yc - thy * xc
    return np.column_stack([ux, uy, uz])


@pytest.mark.parametrize("weighted", [False, True])
def test_fit_section_kinematics_recovers_curvatures_and_twist_rate(weighted):
    rng = np.random.default_rng(3)
    X = np.column_stack([rng.uniform(-1.5, 1.5, 400), rng.uniform(-0.3, 0.3, 400), rng.uniform(-3.0, 3.0, 400)]) + [0, 0, 80.0]
    kappa, tr = (1.3e-3, -4.0e-4), 2.2e-4
    U = _synthetic_field(X, kappa, tr, 80.0)
    w = rng.uniform(0.5, 2.0, len(X)) if weighted else None
    fit = A.fit_section_kinematics(X, U, 80.0, half=3.0, weights=w)
    np.testing.assert_allclose(fit["kappa"], [kappa[0], kappa[1], tr], rtol=1e-8)
    assert fit["rms"] < 1e-12


def test_fit_section_kinematics_rejects_a_thin_slab():
    X = np.column_stack([np.linspace(-1, 1, 10), np.zeros(10), np.zeros(10)])
    with pytest.raises(ValueError, match="slab"):
        A.fit_section_kinematics(X, np.zeros((10, 3)), 0.0, half=1.0)


def test_rigid_end_map_is_u_ref_plus_theta_cross_r():
    pts = np.array([[1.0, 0.0, 5.0], [0.0, 2.0, 5.0]])
    ref = np.array([0.0, 0.0, 5.0])
    T = A.rigid_end_map(pts, ref).toarray()  # (6 n, 6): slave dofs from [u_ref, theta_ref]
    assert T.shape == (12, 6)
    q = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03])
    slave = T @ q
    for i, p in enumerate(pts):
        u = q[:3] + np.cross(q[3:], p - ref)
        np.testing.assert_allclose(slave[6 * i:6 * i + 3], u, atol=1e-15)
        np.testing.assert_allclose(slave[6 * i + 3:6 * i + 6], q[3:], atol=1e-15)


def test_beam_displacements_through_the_fit_recover_a_prismatic_compliance_and_a_tapered_one_at_the_centre():
    rng = np.random.default_rng(5)
    X = np.column_stack([rng.uniform(-1.5, 1.5, 600), rng.uniform(-0.3, 0.3, 600), rng.uniform(60.0, 90.0, 600)])
    C0 = np.array([[6.5e-10, 7e-12, 3e-13], [7e-12, 3.4e-10, -7.4e-10], [3e-13, -7.4e-10, 1.3e-8]])
    ref = X[:, :2].mean(axis=0)
    loads = [A.M0 * np.eye(3)[j] for j in range(3)]
    # prismatic: exact
    Us = [A.beam_displacements(X, lambda z: C0, 60.0, 90.0, m, ref) for m in loads]
    fits = [A.fit_section_kinematics(X, U, 75.0, half=5.0) for U in Us]
    np.testing.assert_allclose(A.compliance_from_responses([f["kappa"] for f in fits], [A.M0] * 3), C0, rtol=1e-4, atol=1e-4 * 6.5e-10)  # beam integration is piecewise linear on 801 points
    # tapered (compliance exponential in z, 4.6 %/m): the fit at the centre is within 1 % of the local value
    taper = lambda z: C0 * np.exp(0.046 * (z - 75.0))  # noqa: E731
    Us = [A.beam_displacements(X, taper, 60.0, 90.0, m, ref) for m in loads]
    fits = [A.fit_section_kinematics(X, U, 75.0, half=3.0) for U in Us]
    Ct = A.compliance_from_responses([f["kappa"] for f in fits], [A.M0] * 3)
    np.testing.assert_allclose(Ct, C0, rtol=1e-2, atol=1e-2 * 6.5e-10)


def test_fit_section_kinematics_is_not_biased_by_in_plane_affine_strain():
    rng = np.random.default_rng(11)
    X = np.column_stack([rng.uniform(-1.5, 1.5, 500), rng.uniform(-0.3, 0.3, 500), rng.uniform(-3.0, 3.0, 500)])
    kappa, tr = (1.0e-3, 2.0e-4), 3.0e-4
    U = _synthetic_field(X, kappa, tr, 0.0)
    dz = X[:, 2]
    xc, yc = X[:, 0] - X[:, 0].mean(), X[:, 1] - X[:, 1].mean()
    U[:, 0] += (-0.3 * 1e-3 * dz) * xc + 4e-4 * dz * yc  # Poisson-type e_xx(dz) and a shear gamma(dz)/2
    U[:, 1] += (-0.3 * 1e-3 * dz) * yc + 4e-4 * dz * xc
    fit = A.fit_section_kinematics(X, U, 0.0, half=3.0)
    np.testing.assert_allclose(fit["kappa"], [kappa[0], kappa[1], tr], rtol=1e-8)


# --- L5 T9c: shell measurement convergence ---------------------------------------------------------------------------------


def _row(z, n, psi, gj, ef=1.0e9, ee=2.0e9, frac=0.25):
    return {"z": z, "n_chord": n, "frac": frac, "corrected": {"psi_deg": psi, "ratio_fe_ff": 0.0, "GJ": gj, "EI_flap": ef, "EI_edge": ee}}


def test_bias_corrected_is_identity_when_the_deck_reads_back_exactly_and_removes_the_fit_bias_otherwise():
    sh = {"psi_deg": 1.0, "ratio_fe_ff": 0.01, "GJ": 5.0e7, "EI_flap": 1.0e9, "EI_edge": 2.0e9}
    ex = {"psi_deg": 0.5, "ratio_fe_ff": 0.005, "GJ": 6.0e7, "EI_flap": 1.1e9, "EI_edge": 2.2e9}
    same = A.bias_corrected(sh, ex, ex)
    for k in sh:
        assert same[k] == pytest.approx(sh[k])
    fit = {"psi_deg": 0.8, "ratio_fe_ff": 0.007, "GJ": 5.4e7, "EI_flap": 1.0e9, "EI_edge": 2.0e9}  # the fit reads the deck 0.3 deg / 10 % off
    c = A.bias_corrected(sh, fit, ex)
    assert c["psi_deg"] == pytest.approx(0.7)
    assert c["ratio_fe_ff"] == pytest.approx(0.008)
    assert c["GJ"] == pytest.approx(5.0e7 * 6.0 / 5.4)
    assert c["EI_flap"] == pytest.approx(1.1e9)


def test_converged_setting_returns_the_smallest_length_whose_later_steps_all_stay_inside_the_tolerances():
    rows = [_row(80, 4, -0.7, 8.0e7), _row(80, 6, 1.0, 8.2e7), _row(80, 8, 1.84, 8.3e7), _row(80, 10, 2.29, 8.4e7), _row(80, 12, 2.45, 8.45e7)]
    got = A.converged_setting(rows, tol_psi=0.5, tol_stiff=0.02)
    assert got == {80.0: 8.0}  # 6 -> 8 moves 0.84 deg, 8 -> 10 0.45, 10 -> 12 0.16
    # a stiffness step above the tolerance after n=8 pushes the answer to n=10
    rows[3]["corrected"]["GJ"] = 9.0e7
    rows[4]["corrected"]["GJ"] = 9.05e7
    assert A.converged_setting(rows, tol_psi=0.5, tol_stiff=0.02) == {80.0: 10.0}


def test_converged_setting_ignores_other_windows_reports_none_when_the_last_step_fails_and_keeps_stations_apart():
    rows = [_row(80, 4, 0.0, 1.0e7), _row(80, 6, 1.0, 1.0e7), _row(95, 4, 0.0, 1.0e7), _row(95, 6, 0.1, 1.0e7), _row(95, 6, 9.0, 1.0e7, frac=0.4)]
    got = A.converged_setting(rows, tol_psi=0.5, tol_stiff=0.02, frac=0.25)
    assert got[95.0] == 4.0
    assert got[80.0] is None  # the only step (4 -> 6) is 1 deg
