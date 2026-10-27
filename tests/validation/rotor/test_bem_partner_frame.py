"""The BEM participant always exchanges inertial-axes data with its structural partner.

The structural partner is always the corotational Rust rotor solver
(``crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs``): it writes Displacement as
``R(theta_target) u`` (elastic, inertial axes) and reads Force as inertial, applying
``R^T(theta_target)``.  So the participant must un-rotate what it reads and rotate what
it writes whether or not the YAML declares ``partner_frame``; the old ``rotating``
default (no transform) fed the solver forces in the wrong axes.  Decision recorded in
``odd/tasks/core-torsion.md`` (C8, D1).  No preCICE, CCBlade or Rust backend is needed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from aeroelast.solvers.bem import fsi_participant as fp

AXIS_Y = np.array([0.0, 1.0, 0.0])
F_ROT = np.array([[120.0, -35.0, 800.0], [3.0, 40.0, -60.0], [-7.0, 1.5, 2500.0]])
U_ROT = np.array([[0.02, -0.01, 0.0], [0.5, 0.1, -0.3], [1.2, 0.0, 0.04]])
DT = 0.01


def _rust_rotation(axis, theta):
    """Independent re-statement of ``RotorTransforms::rotation_matrix``."""
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    k = np.array([[0.0, -n[2], n[1]], [n[2], 0.0, -n[0]], [-n[1], n[0], 0.0]])
    return np.eye(3) + math.sin(theta) * k + (1.0 - math.cos(theta)) * (k @ k)


def _participant(bem_config: dict, theta_end: float):
    """Participant carrying only frame state, configured through the loader's parsers,
    with the end-of-window angle ``theta_end`` the structure publishes."""
    p = object.__new__(fp.BEMFSIParticipant)
    p._partner_frame = fp._parse_partner_frame(bem_config)
    p._rotation_axis = fp._parse_rotation_axis(bem_config)
    p._current_omega = 0.5
    p._partner_theta = theta_end  # published by the structure (D4)
    return p


@pytest.mark.parametrize("theta_deg", [90.0, 180.0])
@pytest.mark.parametrize("cfg", [{}, {"partner_frame": "inertial"}])
def test_forces_written_inertial_by_default(theta_deg, cfg):
    theta = math.radians(theta_deg)
    p = _participant(cfg, theta)
    written = p._forces_to_partner(F_ROT)
    np.testing.assert_allclose(written, F_ROT @ _rust_rotation(AXIS_Y, theta).T, atol=1e-9)
    assert not np.allclose(written, F_ROT)
    # Rust applies R^T(theta_target) on read: it must recover the rotating-frame force.
    np.testing.assert_allclose(written @ _rust_rotation(AXIS_Y, theta), F_ROT, atol=1e-9)


@pytest.mark.parametrize("theta_deg", [90.0, 180.0])
@pytest.mark.parametrize("cfg", [{}, {"partner_frame": "inertial"}])
def test_displacements_read_unrotated_by_default(theta_deg, cfg):
    theta = math.radians(theta_deg)
    p = _participant(cfg, theta)
    written_by_rust = U_ROT @ _rust_rotation(AXIS_Y, theta).T  # R u, row-wise
    np.testing.assert_allclose(p._displacements_from_partner(written_by_rust), U_ROT, atol=1e-9)


def test_rotating_partner_frame_is_rejected():
    with pytest.raises(ValueError, match="corotational") as exc:
        fp._parse_partner_frame({"partner_frame": "rotating"})
    assert "inertial" in str(exc.value)
    assert "partner_frame" in str(exc.value)


def test_inertial_partner_frame_still_accepted():
    assert fp._parse_partner_frame({"partner_frame": "inertial"}) == "inertial"
    assert fp._parse_partner_frame({"partner_frame": " Inertial "}) == "inertial"
    assert fp._parse_partner_frame({}) == "inertial"
