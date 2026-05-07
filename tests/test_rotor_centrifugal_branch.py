"""
Unit test for centrifugal force coordinate selection (Fix #1).

Contract being tested (from spec: Centrifugal Force Evaluation Coordinates):

  - include_ksp=True:  F_cf must equal m·ω²·r₀_⊥  (undeformed radial distance).
    K_SP·u already captures the displacement correction on the LHS; evaluating
    at X₀+u would double-count the O(ω²|u|) term.

  - include_ksp=False: F_cf must equal m·ω²·r_⊥(X₀+u)  (deformed radial distance).
    No LHS counterpart; full geometrically-nonlinear correction must appear in force.

  - The two results differ by an amount proportional to ω²·|u|.

This is a pure-math test that does NOT call the Rust solver.  It validates the
mathematical contract that the Rust branch must implement.  The test can run
in CI without preCICE or the compiled _aeroelast extension.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose


def _centrifugal_force_at(coords: np.ndarray, mass: float, omega: float) -> np.ndarray:
    """
    Pure-Python reference implementation of ω×(ω×r) centrifugal force.

    Rotation axis is Z (0, 0, 1), rotation center is origin.

    For a node at (x, y, z):
      r_perp = (x, y, 0)   (perpendicular component to Z-axis)
      F_cf   = mass · ω² · r_perp   (outward radial direction)

    Parameters
    ----------
    coords : (3,) array
        Nodal coordinates [x, y, z].
    mass : float
        Lumped node mass [kg].
    omega : float
        Angular velocity [rad/s].

    Returns
    -------
    (3,) array of centrifugal force components [N].
    """
    x, y, _ = coords
    # ω×(ω×r) for ω = ω·ẑ is -ω²·(x, y, 0)  in the rotating frame the
    # centrifugal term is +ω²·r_perp (outward).
    return np.array([mass * omega**2 * x, mass * omega**2 * y, 0.0])


# ── Parametrized fixture data ──────────────────────────────────────────────────

# A 2-node "rotor": one node on the X-axis, one off-axis to exercise the Y component.
NODE_A = np.array([10.0, 0.0, 0.0])   # on X-axis, R=10 m
NODE_B = np.array([6.0, 8.0, 0.0])    # R=10 m at 53° from X
MASS_A = 100.0   # kg
MASS_B = 50.0    # kg

# Radial displacements: zero, small, moderate
DISPLACEMENTS = [
    np.array([0.0, 0.0, 0.0]),           # no deformation
    np.array([0.1, 0.0, 0.0]),           # 1% of R=10
    np.array([0.5, 0.0, 0.0]),           # 5% of R=10
    np.array([0.3, 0.4, 0.0]),           # mixed X+Y deformation
]

OMEGA_VALUES = [1.0, 10.0, 50.0, 100.0]


@pytest.mark.parametrize("omega", OMEGA_VALUES)
@pytest.mark.parametrize("u_node", DISPLACEMENTS)
@pytest.mark.parametrize("coords0,mass", [(NODE_A, MASS_A), (NODE_B, MASS_B)])
def test_include_ksp_true_uses_x0(coords0, mass, u_node, omega):
    """
    When include_ksp=True, centrifugal force MUST equal m·ω²·r₀_⊥.

    No displacement contribution. The value is identical regardless of u.
    """
    F_at_x0 = _centrifugal_force_at(coords0, mass, omega)
    F_at_x0_plus_u = _centrifugal_force_at(coords0 + u_node, mass, omega)

    # The correct result for include_ksp=True is F evaluated at X₀.
    F_correct_ksp_on = F_at_x0

    # Verify the analytical X₀ value: F_x = m·ω²·x₀, F_y = m·ω²·y₀, F_z = 0
    assert_allclose(
        F_correct_ksp_on[0], mass * omega**2 * coords0[0], rtol=1e-14,
        err_msg="F_cf_x must equal m·ω²·x₀ when include_ksp=True"
    )
    assert_allclose(
        F_correct_ksp_on[1], mass * omega**2 * coords0[1], rtol=1e-14,
        err_msg="F_cf_y must equal m·ω²·y₀ when include_ksp=True"
    )
    assert_allclose(
        F_correct_ksp_on[2], 0.0, atol=1e-14,
        err_msg="F_cf_z must be zero (Z-axis rotation)"
    )

    # The force must NOT include any contribution from u.
    # Concretely: it must equal F_at_x0, not F_at_x0_plus_u (unless u=0).
    u_nonzero = np.linalg.norm(u_node) > 1e-15
    if u_nonzero:
        with pytest.raises(AssertionError):
            assert_allclose(F_correct_ksp_on, F_at_x0_plus_u, rtol=1e-14)


@pytest.mark.parametrize("omega", OMEGA_VALUES)
@pytest.mark.parametrize("u_node", DISPLACEMENTS)
@pytest.mark.parametrize("coords0,mass", [(NODE_A, MASS_A), (NODE_B, MASS_B)])
def test_include_ksp_false_uses_x0_plus_u(coords0, mass, u_node, omega):
    """
    When include_ksp=False, centrifugal force MUST equal m·ω²·r_⊥(X₀+u).

    Displacement is fully included in the nodal position.
    """
    coords_deformed = coords0 + u_node
    F_at_x0_plus_u = _centrifugal_force_at(coords_deformed, mass, omega)

    # Verify the deformed-coords value
    assert_allclose(
        F_at_x0_plus_u[0], mass * omega**2 * coords_deformed[0], rtol=1e-14,
        err_msg="F_cf_x must equal m·ω²·(x₀+u_x) when include_ksp=False"
    )
    assert_allclose(
        F_at_x0_plus_u[1], mass * omega**2 * coords_deformed[1], rtol=1e-14,
        err_msg="F_cf_y must equal m·ω²·(y₀+u_y) when include_ksp=False"
    )
    assert_allclose(
        F_at_x0_plus_u[2], 0.0, atol=1e-14,
        err_msg="F_cf_z must be zero (Z-axis rotation)"
    )


@pytest.mark.parametrize("omega", OMEGA_VALUES)
@pytest.mark.parametrize("u_node", [u for u in DISPLACEMENTS if np.linalg.norm(u) > 1e-15])
@pytest.mark.parametrize("coords0,mass", [(NODE_A, MASS_A), (NODE_B, MASS_B)])
def test_two_results_differ_proportional_to_omega2_u(coords0, mass, u_node, omega):
    """
    The two branches (include_ksp=True/False) must differ by ~ ω²·|u|.

    Specifically: |F_deformed - F_ref| = m·ω²·|u_perp|
    where u_perp is the component of u perpendicular to the rotation axis (Z).
    """
    F_ksp_on = _centrifugal_force_at(coords0, mass, omega)         # X₀
    F_ksp_off = _centrifugal_force_at(coords0 + u_node, mass, omega)  # X₀+u

    diff = np.linalg.norm(F_ksp_off - F_ksp_on)

    # Perpendicular component of u (in X-Y plane for Z-axis rotation)
    u_perp = np.array([u_node[0], u_node[1], 0.0])
    expected_diff = mass * omega**2 * np.linalg.norm(u_perp)

    assert_allclose(
        diff, expected_diff, rtol=1e-12,
        err_msg=(
            f"Branch diff = {diff:.6g} N, expected m·ω²·|u_perp| = {expected_diff:.6g} N"
        )
    )


@pytest.mark.parametrize("include_ksp", [True, False])
def test_zero_omega_gives_zero_force(include_ksp):
    """
    Centrifugal force must be zero when ω=0, regardless of include_ksp or coords.
    """
    coords0 = NODE_A
    u_node = np.array([0.5, 0.0, 0.0])
    omega = 0.0

    if include_ksp:
        F = _centrifugal_force_at(coords0, MASS_A, omega)
    else:
        F = _centrifugal_force_at(coords0 + u_node, MASS_A, omega)

    assert_allclose(F, 0.0, atol=1e-30, err_msg="F_cf must be zero when ω=0")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
