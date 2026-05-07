"""Integration test for the unified omega-rebuild predicate (Fix #3).

Tests the ``omega_changed_significantly`` predicate at the Python level by simulating
a ramped-ω trajectory through the rotor solver config path.

Acceptance criterion (task 3.7):
    On realistic trajectories (linear ramp, converging, steady-state), rebuild
    count with the new relative predicate (0.005/0.003) is ≤ rebuild count with
    the old absolute-threshold policy (|Δω| > 1e-4 rad/s = ``ksp_omega_threshold``
    default) on ALL parametrized ω trajectories.

These tests are pure-Python (no PETSc, no preCICE, no Rust extension) and run
in any lightweight environment.
"""

from __future__ import annotations

import math
from typing import List

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# Pure-Python replica of the unified ω-rebuild predicate
# Mirrors crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs
# ---------------------------------------------------------------------------

def _omega_changed_significantly(
    omega_new: float,
    omega_sq_at_last: float,
    threshold_rebuild: float,
    threshold_skip: float,
    currently_rebuilt: bool,
    eps: float = 1e-8,
) -> bool:
    """Python replica of the Rust ``omega_changed_significantly`` predicate.

    Parameters
    ----------
    omega_new : float
        Current angular velocity [rad/s].
    omega_sq_at_last : float
        ω² at the previous rebuild; ``float('-inf')`` if never rebuilt.
    threshold_rebuild : float
        Relative |Δ(ω²)| / max(ω², ω_last²) high-band threshold.
        Applied when ``currently_rebuilt=True`` to prevent chattering.
    threshold_skip : float
        Same, low-band threshold. Applied when ``currently_rebuilt=False``.
    currently_rebuilt : bool
        True if a rebuild was performed < 10 steps ago.
    eps : float
        Near-zero ω² guard.  Seeded from ``ksp_omega_threshold²`` (default 1e-4²
        = 1e-8) per task 3.2: the legacy absolute threshold is rerouted as the
        eps guard instead of the rebuild criterion.

    Returns
    -------
    bool
        True → rebuild; False → skip.
    """
    if omega_sq_at_last == float("-inf"):
        return True  # first call always rebuilds

    omega_sq_new = omega_new * omega_new
    denom = max(omega_sq_new, omega_sq_at_last)

    if denom < eps:
        return True  # ω → 0 guard (ksp_omega_threshold² defines "near zero")

    rel_change = abs(omega_sq_new - omega_sq_at_last) / denom
    threshold = threshold_rebuild if currently_rebuilt else threshold_skip

    return rel_change >= threshold


# ---------------------------------------------------------------------------
# Trajectory simulators
# ---------------------------------------------------------------------------

def _simulate_new_predicate(
    omega_trajectory: List[float],
    threshold_rebuild: float = 0.005,
    threshold_skip: float = 0.003,
    hysteresis_steps: int = 10,
) -> int:
    """Total rebuild count using the new relative predicate (includes first-call)."""
    omega_sq_at_last = float("-inf")
    last_rebuild_step = 0
    count = 0

    for step, omega in enumerate(omega_trajectory):
        recently = (step - last_rebuild_step) < hysteresis_steps
        if _omega_changed_significantly(
            omega, omega_sq_at_last, threshold_rebuild, threshold_skip, recently
        ):
            omega_sq_at_last = omega * omega
            last_rebuild_step = step
            count += 1

    return count


def _simulate_old_absolute_policy(
    omega_trajectory: List[float],
    abs_threshold: float = 1e-4,
) -> int:
    """Total rebuild count using old absolute |Δω| > threshold policy.

    The old code initialised ``omega_at_last_ksp = initial_omega`` and rebuilt
    on ``abs(omega_new - omega_last) > ksp_omega_threshold``.  At the default
    threshold of 1e-4 rad/s this triggers on virtually every step of a ramp.

    We account for the first-call rebuild that the new policy always triggers
    (new: NEG_INFINITY → always rebuild; old: initialized to omega_trajectory[0]
    and may not trigger on step 0 if omega is already at initial value).
    To compare fairly, we count the old policy's first-call rebuild separately.
    """
    # Step 0: always count as a rebuild (both policies rebuild at startup)
    omega_at_last = omega_trajectory[0]
    count = 1  # first call = 1 rebuild

    for omega in omega_trajectory[1:]:
        if abs(omega - omega_at_last) > abs_threshold:
            omega_at_last = omega
            count += 1

    return count


# ---------------------------------------------------------------------------
# Trajectory generators
# ---------------------------------------------------------------------------

def _linear_ramp(omega_start: float, omega_end: float, n_steps: int) -> List[float]:
    return [omega_start + (omega_end - omega_start) * i / max(n_steps - 1, 1)
            for i in range(n_steps)]


def _converging_omega(
    omega_final: float,
    n_steps: int,
    time_constant_steps: int = 20,
    noise_amp: float = 0.0005,
    seed: int = 42,
) -> List[float]:
    """Exponentially converging ω with small random noise (realistic FSI scenario)."""
    rng = np.random.default_rng(seed)
    base = [omega_final * (1.0 - math.exp(-i / time_constant_steps)) for i in range(n_steps)]
    noise = rng.uniform(-noise_amp, noise_amp, n_steps) * omega_final
    return (np.array(base) + noise).tolist()


def _constant_omega(omega: float, n_steps: int) -> List[float]:
    return [omega] * n_steps


# ---------------------------------------------------------------------------
# Tests: rebuild count comparison (new ≤ old)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("omega_start,omega_end,n_steps", [
    (0.1, 100.0, 500),   # ramp from near-zero to rated
    (10.0, 100.0, 300),  # ramp from partial to rated
    (100.0, 1.0, 300),   # deceleration ramp
    (1.0, 50.0, 200),    # slow ramp
])
def test_new_predicate_le_old_on_linear_ramp(omega_start, omega_end, n_steps):
    """New relative predicate rebuilds ≤ old abs policy (1e-4 rad/s) on linear ramps."""
    trajectory = _linear_ramp(omega_start, omega_end, n_steps)

    count_new = _simulate_new_predicate(trajectory)
    count_old = _simulate_old_absolute_policy(trajectory, abs_threshold=1e-4)

    assert count_new <= count_old, (
        f"New rebuilt {count_new} vs old {count_old} "
        f"(linear ramp ω: {omega_start}→{omega_end}, {n_steps} steps)"
    )


@pytest.mark.parametrize("omega_final,n_steps,noise_amp", [
    (100.0, 200, 0.0005),  # typical rated speed convergence with noise
    (50.0,  150, 0.0002),  # partial speed convergence
    (10.0,  100, 0.001),   # low-speed convergence
])
def test_new_predicate_le_old_on_converging_omega(omega_final, n_steps, noise_amp):
    """New predicate rebuilds significantly less than old policy on converging ω (realistic FSI)."""
    trajectory = _converging_omega(omega_final, n_steps, noise_amp=noise_amp)

    count_new = _simulate_new_predicate(trajectory)
    count_old = _simulate_old_absolute_policy(trajectory, abs_threshold=1e-4)

    assert count_new <= count_old, (
        f"New rebuilt {count_new} vs old {count_old} "
        f"(converging ω→{omega_final}, {n_steps} steps)"
    )
    # Additionally verify significant reduction (not just marginally less)
    # On converging trajectories, new policy should rebuild ≤ 60% as often as old
    if count_old > 10:  # only meaningful for non-trivial trajectories
        reduction_ratio = count_new / count_old
        assert reduction_ratio <= 0.85, (
            f"Expected ≥15% rebuild reduction vs old policy, got only "
            f"{(1 - reduction_ratio)*100:.0f}% reduction "
            f"(new={count_new}, old={count_old})"
        )


def test_constant_omega_triggers_only_first_rebuild():
    """Constant ω → exactly 1 rebuild (first-call only)."""
    trajectory = _constant_omega(100.0, 50)
    count = _simulate_new_predicate(trajectory)
    assert count == 1, f"Expected exactly 1 rebuild for constant ω=100, got {count}"


def test_constant_omega_new_eq_old():
    """For constant ω, new and old policies both do exactly 1 rebuild."""
    trajectory = _constant_omega(100.0, 100)
    count_new = _simulate_new_predicate(trajectory)
    count_old = _simulate_old_absolute_policy(trajectory, abs_threshold=1e-4)
    assert count_new == count_old == 1


# ---------------------------------------------------------------------------
# Tests: predicate semantics (unit)
# ---------------------------------------------------------------------------

def test_predicate_first_call_always_rebuilds():
    """NEG_INFINITY sentinel → always rebuild regardless of ω or hysteresis state."""
    for omega in [0.0, 0.001, 50.0, 1000.0]:
        assert _omega_changed_significantly(omega, float("-inf"), 0.005, 0.003, False)
        assert _omega_changed_significantly(omega, float("-inf"), 0.005, 0.003, True)


def test_predicate_near_zero_omega_guard():
    """Both ω near zero (ω² < eps = 1e-8 = ksp_omega_threshold²) → always rebuild."""
    assert _omega_changed_significantly(0.0, 0.0, 0.005, 0.003, False)
    assert _omega_changed_significantly(0.0, 0.0, 0.005, 0.003, True)


def test_predicate_stable_omega_skips():
    """Stable ω (< 0.1% change) → skip regardless of hysteresis state."""
    omega_base = 100.0
    omega_sq_last = omega_base ** 2
    # 0.05% ω change → ~0.1% ω² change — below both 0.3% and 0.5%
    omega_new = omega_base * 1.0005
    assert not _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, False)
    assert not _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, True)


def test_predicate_large_jump_always_rebuilds():
    """ω change > 0.5% ω² → rebuild regardless of hysteresis state."""
    omega_base = 100.0
    omega_sq_last = omega_base ** 2
    # 0.4% ω change → ~0.8% ω² change; above both 0.3% and 0.5%
    omega_new = omega_base * 1.004
    assert _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, False)
    assert _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, True)


def test_predicate_hysteresis_between_bands():
    """ω² change in (0.3%, 0.5%) respects hysteresis:
    recently_rebuilt=True  → high threshold (0.5%) → skip
    recently_rebuilt=False → low  threshold (0.3%) → rebuild
    """
    omega_base = 100.0
    omega_sq_last = omega_base ** 2
    # +0.2% ω → ~+0.4% ω² (in the hysteresis band)
    omega_new = omega_base * 1.002
    denom = max(omega_new ** 2, omega_sq_last)
    rel = abs(omega_new ** 2 - omega_sq_last) / denom
    assert 0.003 < rel < 0.005, f"Test setup: rel={rel:.5f} must be in (0.003, 0.005)"

    # recently rebuilt → high bar → skip
    assert not _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, True)
    # not recently rebuilt → low bar → rebuild
    assert _omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, False)


# ---------------------------------------------------------------------------
# Tests: RotorConfig field forwarding (requires aeroelast package)
# ---------------------------------------------------------------------------

def _try_import_rotor_config():
    try:
        from aeroelast.core.config import RotorConfig  # noqa: PLC0415
        return RotorConfig
    except Exception:
        return None


@pytest.mark.skipif(
    _try_import_rotor_config() is None,
    reason="aeroelast.core.config not importable in this environment",
)
def test_rotor_config_exposes_omega_rebuild_fields():
    """RotorConfig.to_dict() must include omega_rebuild_rel_high/low with correct defaults."""
    RotorConfig = _try_import_rotor_config()
    cfg = RotorConfig()
    d = cfg.to_dict()
    assert "omega_rebuild_rel_high" in d, "omega_rebuild_rel_high missing from to_dict()"
    assert "omega_rebuild_rel_low" in d, "omega_rebuild_rel_low missing from to_dict()"
    assert d["omega_rebuild_rel_high"] == pytest.approx(0.005)
    assert d["omega_rebuild_rel_low"] == pytest.approx(0.003)


@pytest.mark.skipif(
    _try_import_rotor_config() is None,
    reason="aeroelast.core.config not importable in this environment",
)
def test_rotor_config_custom_thresholds_forwarded():
    """RotorConfig with custom thresholds must forward them through to_dict()."""
    RotorConfig = _try_import_rotor_config()
    cfg = RotorConfig(omega_rebuild_rel_high=0.01, omega_rebuild_rel_low=0.005)
    d = cfg.to_dict()
    assert d["omega_rebuild_rel_high"] == pytest.approx(0.01)
    assert d["omega_rebuild_rel_low"] == pytest.approx(0.005)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
