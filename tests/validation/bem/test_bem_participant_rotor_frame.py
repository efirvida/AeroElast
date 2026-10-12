"""Frame contract at the BEM participant's preCICE boundary (T7a).

The Rust rotor solver (``crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs``)
writes its interface displacements as ``disps_to_inertial(u, theta_target)``,
i.e. ``R(theta) . u_rot`` with the right-handed Rodrigues matrix
``R = I + sin(theta) K + (1 - cos(theta)) K^2`` (``rotor_physics.rs``), and it
reads the fluid forces as global, applying ``forces_to_rotating`` = ``R^T``.
``theta_target`` is the END-of-window angle.

The BEM participant works in the rotating (blade) frame, so with an inertial
partner it must un-rotate what it reads and rotate back what it writes.  These
tests pin the rotation sense against Rust, the round trip in both directions,
and that the boundary helpers act only when ``partner_frame == "inertial"``.
No preCICE, CCBlade or Rust backend is needed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aeroelast.solvers.bem import fsi_participant as fp

AXIS_Y = np.array([0.0, 1.0, 0.0])
THETAS = [0.0, 0.3, math.pi / 2, 2.0, math.pi, 4.5, 2 * math.pi - 0.1]
# A direction off every coordinate axis, plus an oblique rotation axis.
VECS = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.3, -1.2, 2.5], [-4.0, 0.7, 0.1]])
OBLIQUE_AXIS = np.array([0.2, 0.9, -0.4])


def _rust_rotation(axis, theta):
    """Independent re-statement of ``RotorTransforms::rotation_matrix``."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    k = np.array([[0.0, -n[2], n[1]], [n[2], 0.0, -n[0]], [-n[1], n[0], 0.0]])
    return np.eye(3) + math.sin(theta) * k + (1.0 - math.cos(theta)) * (k @ k)


def test_rodrigues_sense_matches_rust_right_handed():
    r = fp._rodrigues(AXIS_Y, math.pi / 2)
    np.testing.assert_allclose(r @ [1.0, 0.0, 0.0], [0.0, 0.0, -1.0], atol=1e-14)
    np.testing.assert_allclose(r @ [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], atol=1e-14)


@pytest.mark.parametrize("axis", [AXIS_Y, OBLIQUE_AXIS, 3.0 * AXIS_Y])
@pytest.mark.parametrize("theta", THETAS)
def test_rodrigues_equals_rust_matrix(axis, theta):
    np.testing.assert_allclose(fp._rodrigues(axis, theta), _rust_rotation(axis, theta), atol=1e-14)


@pytest.mark.parametrize("axis", [AXIS_Y, OBLIQUE_AXIS])
@pytest.mark.parametrize("theta", THETAS)
def test_displacement_round_trip(axis, theta):
    """Rust writes R.u; un-rotating what the participant reads returns u."""
    written_by_rust = VECS @ _rust_rotation(axis, theta).T  # row-wise R.u
    np.testing.assert_allclose(fp._to_rotating(written_by_rust, axis, theta), VECS, atol=1e-12)


@pytest.mark.parametrize("axis", [AXIS_Y, OBLIQUE_AXIS])
@pytest.mark.parametrize("theta", THETAS)
def test_force_round_trip(axis, theta):
    """The participant writes R.F_rot; Rust's R^T returns F_rot."""
    written = fp._to_inertial(VECS, axis, theta)
    seen_by_rust = written @ _rust_rotation(axis, theta)  # row-wise R^T.F
    np.testing.assert_allclose(seen_by_rust, VECS, atol=1e-12)


def _bare_participant(partner_frame, omega=0.8, frame_theta=1.1, azimuth_deg=200.0):
    """A participant with only the frame state set (no mesh, preCICE or CCBlade)."""
    p = object.__new__(fp.BEMFSIParticipant)
    p._partner_frame = fp._parse_partner_frame({"partner_frame": partner_frame})
    p._rotation_axis = fp._parse_rotation_axis({})
    p._current_omega = omega
    p._frame_theta = frame_theta
    p._azimuth = azimuth_deg
    return p


def test_boundary_helpers_transform_only_for_inertial_partner():
    dt = 0.01
    omega, theta0 = 0.8, 1.1
    theta_end = theta0 + omega * dt  # Rust's theta_target

    inertial = _bare_participant("inertial", omega=omega, frame_theta=theta0)
    u_read = VECS @ _rust_rotation(AXIS_Y, theta_end).T
    np.testing.assert_allclose(inertial._displacements_from_partner(u_read, dt), VECS, atol=1e-12)
    f_written = inertial._forces_to_partner(VECS, dt)
    np.testing.assert_allclose(f_written @ _rust_rotation(AXIS_Y, theta_end), VECS, atol=1e-12)

    rotating = _bare_participant("rotating", omega=omega, frame_theta=theta0)
    np.testing.assert_array_equal(rotating._displacements_from_partner(VECS, dt), VECS)
    np.testing.assert_array_equal(rotating._forces_to_partner(VECS, dt), VECS)


def test_frame_angle_ignores_bem_azimuth_offset():
    """The frame angle follows Rust's theta (starts at 0), not the BEM azimuth."""
    p = _bare_participant("inertial", omega=0.5, frame_theta=0.0, azimuth_deg=120.0)
    assert p._frame_theta_end(0.02) == pytest.approx(0.01)


def test_partner_frame_default_and_validation():
    assert fp._parse_partner_frame({}) == "rotating"
    assert fp._parse_partner_frame({"partner_frame": "Inertial"}) == "inertial"
    with pytest.raises(ValueError, match="partner_frame"):
        fp._parse_partner_frame({"partner_frame": "global"})
    np.testing.assert_allclose(fp._parse_rotation_axis({}), AXIS_Y)
    np.testing.assert_allclose(fp._parse_rotation_axis({"rotation_axis": [0, 0, 2]}), [0, 0, 1])
    with pytest.raises(ValueError, match="rotation_axis"):
        fp._parse_rotation_axis({"rotation_axis": [0, 0, 0]})
