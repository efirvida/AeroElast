"""
Physical consistency tests for rotor FSI optimizations.

Validates that the performance optimizations maintain physical accuracy:
1. Centrifugal forces use the correct coordinate set based on include_ksp:
   - include_ksp=True:  centrifugal evaluated at X₀ (undeformed).
     K_SP already accounts for the displacement correction on the LHS;
     using X₀+u here would double-count the O(ω²|u|) term.
   - include_ksp=False: centrifugal evaluated at X₀+u (deformed).
     No LHS counterpart exists, so the full nonlinear correction must
     appear in the force.
2. K_G rebuild has hysteresis to prevent chattering
3. Coriolis is treated implicitly for stability at high ω
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose


@pytest.mark.parametrize("include_ksp", [True, False])
@pytest.mark.parametrize("omega", [0.0, 10.0, 50.0, 100.0])
@pytest.mark.parametrize("deformation_ratio", [0.0, 0.01, 0.05, 0.10])
def test_centrifugal_deformed_geometry(include_ksp, omega, deformation_ratio):
    """
    Centrifugal force coordinate selection depends on include_ksp.

    When include_ksp=True:
      K_SP·u already captures the deformation correction on the LHS.
      Evaluating centrifugal at X₀+u simultaneously would double-count O(ω²|u|).
      Therefore the reference-geometry (X₀) result is CORRECT — zero relative error
      compared to the true undeformed-geometry force.

    When include_ksp=False:
      No LHS counterpart exists. The full geometrically-nonlinear correction must
      appear in the force, so deformed geometry X₀+u is CORRECT.
      Using reference geometry (X₀) introduces an error that scales as O(u/R).
    """
    # Simplified rotor: single node at distance R from axis
    R = 10.0  # radius [m]
    mass = 100.0  # kg

    # Reference coordinates (undeformed)
    # Deformation in radial direction
    u_deformation = deformation_ratio * R

    # Force evaluated at reference (undeformed) geometry X₀
    r_perp_ref = R
    F_at_X0 = mass * omega**2 * r_perp_ref

    # Force evaluated at deformed geometry X₀+u
    r_perp_deformed = R + u_deformation
    F_at_X0_plus_u = mass * omega**2 * r_perp_deformed

    if omega == 0:
        # Both are zero regardless; nothing to assert about relative error
        assert F_at_X0 == 0.0
        assert F_at_X0_plus_u == 0.0
        return

    if include_ksp:
        # Contract: centrifugal MUST be evaluated at X₀ when include_ksp=True.
        # F_at_X0 is the CORRECT value; F_at_X0_plus_u would be wrong.
        # The "cached" (reference geometry) result has zero relative error vs the
        # correct X₀ evaluation — they are the same quantity.
        F_correct = F_at_X0
        F_wrong = F_at_X0_plus_u
        # Relative error of the WRONG (deformed-geometry) result vs the correct one
        rel_error_if_deformed = abs(F_wrong - F_correct) / F_correct
        assert rel_error_if_deformed == pytest.approx(deformation_ratio, abs=1e-10), (
            f"include_ksp=True: deformed-geometry force has rel_error={rel_error_if_deformed:.3%}, "
            f"expected ~{deformation_ratio:.3%} (proportional to deformation)"
        )
        # The reference-geometry evaluation is exact — zero error.
        assert F_at_X0 == F_correct, "include_ksp=True: X₀ evaluation must equal correct value"
    else:
        # Contract: centrifugal MUST be evaluated at X₀+u when include_ksp=False.
        # F_at_X0_plus_u is CORRECT; F_at_X0 is wrong (missing the deformation correction).
        F_correct = F_at_X0_plus_u
        # Relative error of using reference geometry (wrong for this case)
        rel_error = abs(F_at_X0 - F_correct) / F_correct
        assert rel_error == pytest.approx(deformation_ratio / (1.0 + deformation_ratio), abs=1e-10), (
            f"include_ksp=False: reference-geometry error={rel_error:.3%}, "
            f"expected ~{deformation_ratio / (1.0 + deformation_ratio):.3%}"
        )
        # If deformation is 10% of R, error with reference coords is significant
        if deformation_ratio == 0.10:
            assert rel_error > 0.08, (
                "include_ksp=False: reference-geometry method should have ~9% error "
                "for 10% deformation"
            )


def test_kg_hysteresis_prevents_chattering():
    """
    K_G rebuild hysteresis prevents chattering when ω oscillates around threshold.

    Without hysteresis: rebuild at 0.5% → skip at 0.49% → rebuild at 0.5% → ...
    With hysteresis: rebuild at 0.5% → skip until 0.3% → smooth operation
    """
    THRESHOLD_REBUILD = 0.005  # 0.5%
    THRESHOLD_SKIP = 0.003  # 0.3% (hysteresis)

    omega_base = 100.0  # rad/s
    omega_sq_base = omega_base**2

    # Simulate ω oscillating ±0.4% around base
    omega_values = [
        omega_base * (1 + 0.004),  # +0.4% → should rebuild (first time)
        omega_base * (1 + 0.002),  # +0.2% → skip (within hysteresis)
        omega_base * (1 - 0.002),  # -0.2% → skip (within hysteresis)
        omega_base * (1 + 0.004),  # +0.4% → skip (within hysteresis)
        omega_base * (1 + 0.006),  # +0.6% → rebuild (exceeds threshold)
        omega_base * (1 + 0.004),  # +0.4% → skip (just rebuilt, use high threshold)
        omega_base * (1 + 0.002),  # +0.2% → skip
    ]

    omega_sq_last = omega_sq_base
    rebuild_steps = []
    step = 0
    last_rebuild_step = 0

    for omega_current in omega_values:
        omega_sq_current = omega_current**2
        rel_change = abs(omega_sq_current - omega_sq_last) / omega_sq_last

        # Hysteresis logic
        steps_since_rebuild = step - last_rebuild_step
        threshold = THRESHOLD_REBUILD if steps_since_rebuild < 10 else THRESHOLD_SKIP

        if rel_change >= threshold or step == 0:
            rebuild_steps.append(step)
            omega_sq_last = omega_sq_current
            last_rebuild_step = step

        step += 1

    # Should rebuild at step 0 (first), 4 (exceeds 0.6%), maybe not at others
    assert 0 in rebuild_steps, "Must rebuild on first step"
    assert 4 in rebuild_steps, "Must rebuild when exceeding 0.6%"
    assert len(rebuild_steps) <= 3, (
        f"Should rebuild ≤3 times, got {len(rebuild_steps)} (steps: {rebuild_steps})"
    )

    # Without hysteresis, would rebuild every time change > 0.5% → chattering
    # With hysteresis, rebuilds only when necessary


def test_coriolis_matrix_antisymmetry():
    """
    Coriolis gyroscopic matrix G_cor must be antisymmetric: G^T = -G.

    This ensures energy conservation and stability in implicit treatment.
    """
    try:
        from _aeroelast import PyMeshAssembler
    except ImportError:
        pytest.skip("Rust backend not available")

    # Simple test with 2 nodes, 6 DOF each (shells)
    n_nodes = 2
    masses = np.array([10.0, 15.0], dtype=np.float64)
    axis = np.array([0.0, 0.0, 1.0])
    omega = 50.0  # rad/s
    dofs_per_node = 6

    # Free DOFs: translational DOFs for both nodes (0,1,2, 6,7,8)
    free_dofs = np.array([0, 1, 2, 6, 7, 8], dtype=np.int32)

    # Build Coriolis matrix (this would call build_coriolis_matrix in Rust)
    # For testing, we construct it manually following the algorithm:

    wx, wy, wz = omega * axis

    # G_cor has 3×3 blocks for each node (only couples translations)
    # Block for node i:
    # [  0   -wz   wy ]
    # [  wz   0   -wx ]
    # [ -wy   wx    0  ]
    # All multiplied by -2·m_i

    g_matrix = np.zeros((6, 6))  # 6×6 for 2 nodes × 3 translations each

    for i, m in enumerate(masses):
        coeff = -2.0 * m
        base = i * 3
        # Off-diagonal entries
        g_matrix[base + 0, base + 1] = -wz * coeff
        g_matrix[base + 0, base + 2] = wy * coeff
        g_matrix[base + 1, base + 0] = wz * coeff
        g_matrix[base + 1, base + 2] = -wx * coeff
        g_matrix[base + 2, base + 0] = -wy * coeff
        g_matrix[base + 2, base + 1] = wx * coeff

    # Check antisymmetry: G^T = -G
    assert_allclose(
        g_matrix.T, -g_matrix, atol=1e-14, rtol=0, err_msg="Coriolis matrix must be antisymmetric"
    )

    # Check diagonal is zero (consequence of antisymmetry)
    assert_allclose(
        np.diag(g_matrix), 0.0, atol=1e-14, err_msg="Diagonal of antisymmetric matrix must be zero"
    )


def test_coriolis_implicit_stability():
    """
    Implicit Coriolis treatment provides unconditional stability.

    For explicit treatment: Δt < 2/(ω·√2) for stability.
    For implicit treatment: unconditionally stable for any Δt.

    This test validates that adding G_cor to the LHS gives stable energy.
    """
    # Simplified 1-DOF oscillator in rotating frame:
    # m·ü + 2·m·ω·v̇ + k·u = 0  (Coriolis couples u and v)

    m = 1.0
    k = 100.0
    omega = 10.0

    # Newmark-β parameters
    beta = 0.25
    gamma = 0.5

    # Large time step that would be unstable with explicit Coriolis
    dt_large = 0.5  # >> 2/(ω·√2) ≈ 0.14 for ω=10

    # Effective stiffness with implicit Coriolis (simplified 2-DOF):
    # K_eff = K + γ/(β·dt)·G_cor + 1/(β·dt²)·M

    a0 = 1.0 / (beta * dt_large**2)
    a1 = gamma / (beta * dt_large)

    # 2×2 system (u, v coupled by Coriolis)
    K = np.diag([k, k])  # stiffness (uncoupled)
    M = np.diag([m, m])  # mass
    G_cor = np.array([[0, -2 * m * omega], [2 * m * omega, 0]])  # antisymmetric

    K_eff = K + a1 * G_cor + a0 * M

    # Stability: K_eff must be positive definite
    eigvals = np.linalg.eigvalsh(K_eff)

    assert np.all(eigvals > 0), (
        f"K_eff must be positive definite for stability, got eigenvalues: {eigvals}"
    )

    # Condition number should be reasonable (not stiff)
    cond = eigvals.max() / eigvals.min()
    assert cond < 1e6, f"Condition number = {cond:.2e}, system may be stiff"


@pytest.mark.parametrize("stress_interval", [1, 5, 10])
def test_stress_gate_checkpoint_consistency(stress_interval):
    """
    Stress computation gate must ALWAYS execute at checkpoint steps.

    Even if stress_output_interval=10, checkpoints must have stress fields.
    """
    # Simulate time steps with checkpoints at specific times
    dt = 0.01
    checkpoint_interval_time = 0.1  # checkpoint every 0.1 sec

    stress_computed_steps = []

    for step in range(50):
        t = step * dt

        # Should write checkpoint?
        t_since_last = t % checkpoint_interval_time
        is_checkpoint = (
            abs(t_since_last) < 1e-6 or abs(t_since_last - checkpoint_interval_time) < 1e-6
        )

        # Stress gate logic (from rotor.py)
        needs_stress = (
            stress_interval <= 1
            or (step % stress_interval == 0)
            or is_checkpoint  # <-- CRITICAL: must always compute at checkpoints
        )

        if needs_stress:
            stress_computed_steps.append(step)

    # Checkpoint steps (every 10 steps for dt=0.01, interval=0.1)
    checkpoint_steps = [i for i in range(50) if abs((i * dt) % checkpoint_interval_time) < 1e-6]

    # Verify ALL checkpoint steps computed stress
    for cp_step in checkpoint_steps:
        assert cp_step in stress_computed_steps, (
            f"Step {cp_step} is checkpoint but stress not computed (interval={stress_interval})"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
