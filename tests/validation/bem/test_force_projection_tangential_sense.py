"""Guard: the tangential load sense follows the measured leading edge (T7c).

``Tp`` is the in-plane force that drives the rotor, so for ``Omega = +omega * y`` and
a blade along ``+z`` (``omega x r = +x``) it must point toward the measured leading
edge.  On this mesh the leading edge is at ``+x`` at every station.  A configured
``tangential_direction`` is therefore not a free choice: one that contradicts the
measured edge reverses the torque (static FSI: +5.55 % thrust with ``-x``, -21.6 %
with ``+x``), so the projector refuses it instead of trusting it.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
pytest.importorskip("_aeroelast", reason="Rust backend not available")

from aeroelast.solvers.bem.force_projection import ForceProjector  # noqa: E402

from tests.validation.bem.test_force_projection_load_frame import (  # noqa: E402,F401
    NORMAL_DIR,
    RATED_OMEGA_RPM,
    SPAN_DIR,
    _max_extent_direction,
    _orient_leading_to_trailing,
    _physical_stations,
    Z_MERGE_TOL,
    blade_case,
    rated_bem,
)

MID_SPAN_FRACTIONS = (0.25, 0.50, 0.75)


def _measured_leading_edge(coords):
    """Mean leading-edge direction from the ring outlines alone (mid-span rings).

    ``_orient_leading_to_trailing`` returns the axis pointing *toward* the blunt
    (leading) end, despite its name.
    """
    stations = _physical_stations(coords)
    z0, z1 = stations[0], stations[-1]
    total = np.zeros(3)
    for frac in MID_SPAN_FRACTIONS:
        z = min(stations, key=lambda s: abs(s - (z0 + frac * (z1 - z0))))
        pts = coords[np.abs(coords[:, 2] - z) < Z_MERGE_TOL]
        axis, _ = _max_extent_direction(pts)
        total += _orient_leading_to_trailing(pts, axis, SPAN_DIR)
    return total / np.linalg.norm(total)


def test_contradicting_tangential_direction_is_refused(blade_case):
    mesh, _, blade_aero, _, _ = blade_case
    with pytest.raises(ValueError, match="tangential_direction") as err:
        ForceProjector(
            mesh,
            blade_aero,
            span_direction=SPAN_DIR,
            normal_direction=NORMAL_DIR,
            tangential_direction=[-1.0, 0.0, 0.0],
        )
    assert "leading edge" in str(err.value)
    assert "drop the key" in str(err.value).lower()


def test_tangential_force_points_to_measured_leading_edge(blade_case, rated_bem):
    mesh, coords, blade_aero, _, _ = blade_case
    projector = ForceProjector(mesh, blade_aero, span_direction=SPAN_DIR)
    le = _measured_leading_edge(coords)
    assert float(le @ np.array([1.0, 0.0, 0.0])) > 0.9, "mesh convention changed: LE not at +x"
    assert float(projector._rotor_tangential_dir @ le) > 0.9


def test_aerodynamic_torque_about_rotation_axis_is_positive(blade_case, rated_bem):
    mesh, coords, blade_aero, _, _ = blade_case
    projector = ForceProjector(mesh, blade_aero, span_direction=SPAN_DIR)
    assert float(np.min(rated_bem.Tp)) > -1e-9 or float(np.sum(rated_bem.Tp)) > 0.0
    forces = projector.project(rated_bem)
    r = coords + blade_aero.hub_radius * SPAN_DIR
    torque_y = float(np.cross(r, forces)[:, 1].sum())
    assert RATED_OMEGA_RPM > 0.0
    assert torque_y > 0.0, f"the structure is braked: torque about +y = {torque_y:.4e} N m"
