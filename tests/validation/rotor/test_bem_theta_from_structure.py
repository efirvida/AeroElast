"""The structure publishes the rotor angle theta; the BEM participant never integrates its own.

Decision D4 (``odd/tasks/core-torsion.md``, S8): the corotational Rust rotor solver writes
``theta_target`` (the end-of-window angle it rotates forces and displacements with) on the omega
mesh next to omega.  ``BEMFSIParticipant`` reads it and uses that one value (a) for the frame
rotation of the read displacements and the written forces and (b) as the CCBlade azimuth:

    azimuth_deg = (bem.azimuth + degrees(theta)) mod 360

Why that mapping.  CCBlade (``bem.f90 windComponents``) measures the azimuth as a right-handed
rotation about the hub +x (shaft) axis starting from the blade pointing up (+z hub):
``heightFromHub = y_az sin(az) + z_az cos(az)``, and the rotation wind is ``-Omega x r`` with
Omega along +x, so a positive rotor speed raises the azimuth.  The structure rotates by
``R(axis, theta)`` with ``theta`` increasing for positive omega.  Both are the same positive
rotation, so ``d azimuth = + d theta`` independently of the coordinate labels; ``bem.azimuth`` is
the azimuth at ``theta = 0``.  No preCICE, CCBlade or Rust backend is needed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aeroelast.solvers.bem import fsi_participant as fp

AXIS_Y = np.array([0.0, 1.0, 0.0])
F_ROT = np.array([[120.0, -35.0, 800.0], [3.0, 40.0, -60.0], [-7.0, 1.5, 2500.0]])
U_ROT = np.array([[0.02, -0.01, 0.0], [0.5, 0.1, -0.3], [1.2, 0.0, 0.04]])


def _rust_rotation(axis, theta):
    """Independent re-statement of ``RotorTransforms::rotation_matrix``."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    k = np.array([[0.0, -n[2], n[1]], [n[2], 0.0, -n[0]], [-n[1], n[0], 0.0]])
    return np.eye(3) + math.sin(theta) * k + (1.0 - math.cos(theta)) * (k @ k)


class _FakeAdapter:
    """Serves the omega mesh the way the Rust solver writes it: omega and theta scalars."""

    def __init__(self, omega, theta, mesh="GlobalFluidMesh"):
        self.mesh = mesh
        self.omega = omega
        self.theta = theta
        self.reads = []

    def read_data(self, mesh, data):
        self.reads.append((mesh, data))
        assert mesh == self.mesh
        return np.array([{"AngularVelocity": self.omega, "Theta": self.theta}[data]])


def _participant(azimuth_offset=0.0, omega=0.0, theta_data="Theta"):
    p = object.__new__(fp.BEMFSIParticipant)
    p._rotation_axis = fp._parse_rotation_axis({})
    p._partner_frame = "inertial"
    p._omega_mesh = "GlobalFluidMesh"
    p._omega_data = "AngularVelocity"
    p._theta_data = theta_data
    p._current_omega = omega
    p._azimuth_offset = azimuth_offset
    p._partner_theta = 0.0
    p._window_count = 0
    p._iteration_count = 0
    return p


@pytest.mark.parametrize("theta", [0.0, 0.3, math.pi / 2, 2.5, math.pi])
def test_frame_rotation_uses_exactly_the_published_theta(theta):
    p = _participant(omega=0.9)
    p._read_rotor_state(_FakeAdapter(omega=0.9, theta=theta))
    assert p._partner_theta == theta
    written = p._forces_to_partner(F_ROT)
    np.testing.assert_allclose(written, F_ROT @ _rust_rotation(AXIS_Y, theta).T, atol=1e-12)
    np.testing.assert_allclose(
        p._displacements_from_partner(U_ROT @ _rust_rotation(AXIS_Y, theta).T), U_ROT, atol=1e-12
    )


def test_restart_with_nonzero_initial_theta_is_followed_from_the_first_read():
    """A Rust restart starts at theta0 != 0; the first iteration already rotates with it."""
    theta0 = 4.2
    p = _participant(azimuth_offset=10.0, omega=0.5)
    p._read_rotor_state(_FakeAdapter(omega=0.5, theta=theta0))
    np.testing.assert_allclose(
        p._forces_to_partner(F_ROT), F_ROT @ _rust_rotation(AXIS_Y, theta0).T, atol=1e-12
    )
    assert p._azimuth_deg() == pytest.approx((10.0 + math.degrees(theta0)) % 360.0)


@pytest.mark.parametrize(
    "offset, theta, expected",
    [
        (0.0, 0.0, 0.0),
        (0.0, math.pi / 2, 90.0),
        (30.0, math.pi / 2, 120.0),
        (350.0, math.radians(20.0), 10.0),
        (0.0, 2.0 * math.pi + math.radians(15.0), 15.0),
        (0.0, -math.radians(30.0), 330.0),
    ],
)
def test_ccblade_azimuth_is_offset_plus_theta_in_degrees(offset, theta, expected):
    p = _participant(azimuth_offset=offset)
    p._read_rotor_state(_FakeAdapter(omega=0.4, theta=theta))
    assert p._azimuth_deg() == pytest.approx(expected)


def test_azimuth_increases_with_positive_rotation():
    """Positive omega: theta and the CCBlade azimuth both increase (same rotation sense)."""
    p = _participant()
    azimuths = []
    for theta in (0.0, 0.01, 0.02, 0.03):
        p._read_rotor_state(_FakeAdapter(omega=0.8, theta=theta))
        azimuths.append(p._azimuth_deg())
    assert np.all(np.diff(azimuths) > 0.0)


def test_varying_omega_does_not_move_the_angle():
    """The participant does not integrate omega: theta from the partner is the only source."""
    p = _participant()
    out = []
    for omega in (0.1, 5.0, 123.0, -2.0):
        p._read_rotor_state(_FakeAdapter(omega=omega, theta=1.234))
        assert p._current_omega == omega
        out.append((p._forces_to_partner(F_ROT), p._azimuth_deg()))
    for forces, azimuth in out[1:]:
        np.testing.assert_array_equal(forces, out[0][0])
        assert azimuth == out[0][1]


def test_theta_and_omega_are_read_from_the_omega_mesh():
    p = _participant(theta_data="Angle")
    adapter = _FakeAdapter(omega=0.7, theta=0.2)
    adapter.read_data = lambda mesh, data: np.array(  # noqa: E731
        [{"AngularVelocity": 0.7, "Angle": 0.2}[data]]
    )
    p._read_rotor_state(adapter)
    assert p._partner_theta == pytest.approx(0.2)


def test_missing_theta_data_fails_at_setup_with_a_clear_message():
    class _Precice:
        def get_data_dimensions(self, mesh, data):
            if data == "Theta":
                raise RuntimeError("unknown data")
            return 1

    class _Adapter:
        precice = _Precice()

    p = _participant()
    with pytest.raises(RuntimeError) as exc:
        p._require_partner_theta(_Adapter())
    msg = str(exc.value)
    assert "Theta" in msg and "GlobalFluidMesh" in msg and "precice" in msg.lower()


def test_theta_data_present_passes_setup_check():
    class _Precice:
        def get_data_dimensions(self, mesh, data):
            return 1

    class _Adapter:
        precice = _Precice()

    _participant()._require_partner_theta(_Adapter())


def test_nonzero_static_omega_without_omega_mesh_is_rejected():
    """Without the omega mesh there is no theta; a spinning rotor cannot be followed."""
    with pytest.raises(ValueError, match="theta"):
        fp._check_theta_source(omega_mesh=None, omega_initial=0.79)
    fp._check_theta_source(omega_mesh=None, omega_initial=0.0)  # parked: theta stays 0
    fp._check_theta_source(omega_mesh="GlobalFluidMesh", omega_initial=0.79)
