"""
Regression tests for PETSc factorization reuse helpers in the inertial rotor solver.

DEPRECATED: The Python methods tested here (_assemble_inertial_effective_system,
_configure_reusable_factorization_ksp, _fallback_factorization_type) were removed
in the Rust migration (inertial-rotor-rust-migration). All FSI logic now lives in
crates/aeroelast-solvers/src/petsc/fsi/rotor_inertial.rs.
This file is kept as a placeholder; delete it once the Rust integration tests
(Phase 6 in the SDD task list) provide equivalent coverage.
"""

import pytest

pytestmark = pytest.mark.skip(
    reason=(
        "Dead tests: _assemble_inertial_effective_system, "
        "_configure_reusable_factorization_ksp and _fallback_factorization_type "
        "were deleted in the Rust migration. "
        "Replace with Rust integration tests in Phase 6."
    )
)
