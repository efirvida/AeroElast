"""Tests for rotor_performance.csv logging in the rotor FSI solver."""

from __future__ import annotations

import csv
import sys
from enum import Enum
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("petsc4py", reason="PETSc not available")

try:
    import precice  # noqa: F401
except (ImportError, OSError):

    class _DummyParticipant:
        def __init__(self, *args, **kwargs):
            pass

    sys.modules["precice"] = SimpleNamespace(Participant=_DummyParticipant)

try:
    import _aeroelast  # noqa: F401
except (ImportError, OSError):

    class _DummyElementFamily(Enum):
        PLANE = "plane"
        SHELL = "shell"
        SOLID = "solid"

    sys.modules["_aeroelast"] = SimpleNamespace(ElementFamily=_DummyElementFamily)

try:
    from aeroelast.solvers.fsi.rotor import LinearDynamicFSIRotorSolver

    _HAS_FSI = True
except (ImportError, OSError):
    _HAS_FSI = False
    LinearDynamicFSIRotorSolver = None  # type: ignore[assignment,misc]

_skip_fsi = pytest.mark.skipif(not _HAS_FSI, reason="preCICE shared library not available")


@_skip_fsi
def test_rotor_performance_report_creates_output_folder_and_csv(tmp_path):
    solver = object.__new__(LinearDynamicFSIRotorSolver)
    solver.solver_params = {"output_folder": str(tmp_path / "nested" / "results")}
    solver._n_blades = 3
    solver._is_primary_rank = lambda: True

    solver._log_rotor_performance(
        t=0.25,
        omega_rpm=12.5,
        omega_rad=1.308996939,
        alpha=0.0,
        angle_deg=15.0,
        thrust=100.0,
        torque_aero=20.0,
        torque_non_aero=-2.0,
        torque_inertial=1.0,
        torque_gravity=1.0,
        torque_total=18.0,
        power_aero=26.17993878,
        power_total=23.561944902,
        structural_efficiency=0.1,
        cp=0.45,
        cq=0.12,
        ct=0.8,
        tsr=6.5,
        torque_aero_global=np.array([1.0, 2.0, 3.0], dtype=np.float64),
        torque_total_global=np.array([4.0, 5.0, 6.0], dtype=np.float64),
        max_displacement=0.02,
        deformed_radius=7.5,
    )

    csv_path = tmp_path / "nested" / "results" / "rotor_performance.csv"
    assert csv_path.exists()

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader)

    assert reader.fieldnames is not None
    assert "Aero Thrust [N]" in reader.fieldnames
    assert "Aero Power Single Blade [W]" in reader.fieldnames
    assert "Aero Power Rotor Equivalent [W]" in reader.fieldnames
    assert "Max Displacement [m]" in reader.fieldnames
    assert float(row["Time [s]"]) == pytest.approx(0.25)
    assert float(row["Aero Torque [Nm]"]) == pytest.approx(20.0)
    assert float(row["Aero Power Single Blade [W]"]) == pytest.approx(26.17993878)
    assert float(row["Aero Power Rotor Equivalent [W]"]) == pytest.approx(78.53981634)
    assert float(row["Blade Count [-]"]) == pytest.approx(3.0)
    assert float(row["Total Torque Z [Nm]"]) == pytest.approx(6.0)
    assert float(row["Max Displacement [m]"]) == pytest.approx(0.02)
    assert float(row["Deformed Radius [m]"]) == pytest.approx(7.5)


@_skip_fsi
def test_rotor_power_phase_sum_matches_analytic_three_phase_total():
    solver = object.__new__(LinearDynamicFSIRotorSolver)
    solver._n_blades = 3
    solver._init_state_tracking()

    omega = 0.7906341464750989  # 7.55 RPM [rad/s]
    period = 2.0 * np.pi / omega
    dt = 0.01
    n_steps = 5000  # ~6 revolutions
    t = np.arange(n_steps) * dt

    # Synthetic single-blade power: mean + 1P + 2P azimuthal harmonics.
    mean_p = 5.3e6
    single = mean_p * (1.0 + 0.3 * np.sin(omega * t) + 0.1 * np.cos(2.0 * omega * t))

    rotor_aero = []
    rotor_total = []
    for i, ti in enumerate(t):
        pa = float(single[i])
        pt = float(single[i] * 0.9)  # structural losses keep ratio constant
        aero_rotor, total_rotor = solver._rotor_power_phase_sum(
            t=ti, power_aero=pa, power_total=pt, omega_rad=omega, n_blades=3
        )
        rotor_aero.append(aero_rotor)
        rotor_total.append(total_rotor)

    rotor_aero = np.array(rotor_aero)
    rotor_total = np.array(rotor_total)

    # Mean of the rotor signal must equal 3× the single-blade mean once the
    # history covers a full revolution (before that, the legacy 3× fallback
    # samples a partial window whose harmonic mean is biased).
    steady_mask = t >= 2.0 * period
    assert rotor_aero[steady_mask].mean() == pytest.approx(3.0 * mean_p, rel=1e-3)
    assert rotor_total[steady_mask].mean() == pytest.approx(3.0 * 0.9 * mean_p, rel=1e-3)

    # Instantaneous check in steady state: the phase sum reproduces the
    # analytic 3-blade total — the sum over the three equally-spaced phases
    # of the 1P+2P harmonics cancels exactly, so the rotor signal is flat
    # ≈ 3·mean_p.
    steady = rotor_aero[int(4.0 * period / dt):]
    assert steady.std() / steady.mean() < 1e-3

    # Legacy fallback: without history (fresh solver), the returned value is
    # 3× the single-blade sample.
    fresh = object.__new__(LinearDynamicFSIRotorSolver)
    fresh._n_blades = 3
    aero_r, total_r = fresh._rotor_power_phase_sum(
        t=0.25, power_aero=2.0, power_total=1.0, omega_rad=omega, n_blades=3
    )
    assert aero_r == pytest.approx(6.0)
    assert total_r == pytest.approx(3.0)


@_skip_fsi
def test_rotor_power_phase_sum_interpolates_between_logged_samples():
    solver = object.__new__(LinearDynamicFSIRotorSolver)
    solver._n_blades = 3
    solver._init_state_tracking()

    omega = 0.7906341464750989
    period = 2.0 * np.pi / omega
    dt = 0.02
    n_steps = int(np.ceil(4.0 * period / dt))
    t = np.arange(n_steps) * dt

    single = 5.0e6 * (1.0 + 0.3 * np.sin(omega * t))
    for i, ti in enumerate(t):
        solver._rotor_power_phase_sum(
            t=ti, power_aero=float(single[i]), power_total=float(single[i]), omega_rad=omega, n_blades=3
        )

    # Request a total at a time BETWEEN logged samples: the interpolation
    # path must not raise and must stay within the plausible envelope.
    t_query = t[int(2.5 * period / dt)] + 0.007
    aero_r, total_r = solver._rotor_power_phase_sum(
        t=t_query, power_aero=float(single[0]), power_total=float(single[0]), omega_rad=omega, n_blades=3
    )
    envelope = 3.0 * 5.0e6 * (1.0 + 0.4)
    assert 0.0 < aero_r < envelope
    assert 0.0 < total_r < envelope
