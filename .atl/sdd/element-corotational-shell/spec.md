# Spec: element-corotational-shell

## Requirements

### R-01 — SVD polar decomposition
Given a deformation gradient F, the system MUST compute polar decomposition F = R·S
via SVD with det(R) = +1 and error ||F − R·S||_F < 1e-12.

### R-02 — Corotational frame for MITC4
Given nodal displacements u (24 DOFs), MITC4 MUST compute the deformed mid-plane
frame {ê₁_def, ê₂_def, ê₃_def} and the 24×24 transformation matrix T.

### R-03 — Corotational frame for MITC3
Same as R-02 for MITC3 (18 DOFs).

### R-04 — Corotational tangent stiffness MITC4
MITC4 MUST expose `compute_kt_corotational(u_elem, coords_elem) → Mat24` implementing
K_T^coro = Tᵀ · K_L · T + K_σ (geometric stiffness in deformed frame).

### R-05 — Corotational tangent stiffness MITC3
Same as R-04 for MITC3 → Mat18.

### R-06 — Assembler integration
`MeshAssembler` MUST expose `assemble_kt_corotational(u: &[f64]) → CooMatrix`
using the same parallelism as `assemble_kt()`.

### R-07 — Config flags
`RotorConfig` MUST accept `use_corotational_kt: bool` (default false) and
`kt_update_freq: usize` (default 1). These MUST be passed through to `rotor_fsi.rs`
and `rotor_inertial.rs`.

### R-08 — Backward compatibility
With `use_corotational_kt: false` (or absent), behavior MUST be identical to the
current TL implementation.

### R-09 — Solver wiring
When `use_corotational_kt: true`, the FSI rotor solver MUST call
`assemble_kt_corotational()` every `kt_update_freq` windows instead of `assemble_kt()`.

---

## Scenarios

### SC-001 — SVD decomposes identity
- Given: F = I₃ₓ₃
- When: `polar_decomposition_svd(F)` is called
- Then: R = I, S = I, error < 1e-12

### SC-002 — SVD decomposes pure rotation 30°
- Given: F = Rz(30°) (no stretch)
- When: `polar_decomposition_svd(F)` is called
- Then: R = Rz(30°), S = I, error < 1e-12, det(R) = +1

### SC-003 — SVD decomposes pure rotation 90°
- Same as SC-002 for θ = 90°

### SC-004 — SVD decomposes pure rotation 180°
- Same as SC-002 for θ = 180°

### SC-005 — SVD decomposes combined rotation+stretch
- Given: F = Rz(45°) · diag(1.1, 0.9, 1.0)
- When: `polar_decomposition_svd(F)` is called
- Then: R is orthogonal, det(R) = +1, ||F − R·S||_F < 1e-12

### SC-006 — Old Denman-Beavers test removed
- Given: test `test_polar_decomposition` using tol=0.05
- When: replaced by SVD tests
- Then: no test references the old DB iteration

### SC-010 — T24 frame for undeformed element
- Given: MITC4 with u = 0 (all zeros)
- When: `build_t24_deformed(frame_def)` is called
- Then: T = I₂₄ₓ₂₄ (identity, frame equals reference)

### SC-011 — T24 frame for translated element (rigid body)
- Given: MITC4 with uniform translation u_t (no rotation)
- When: `build_t24_deformed(frame_def)` is called
- Then: T = I₂₄ₓ₂₄ (pure translation doesn't change frame)

### SC-012 — Invariance under rigid rotation (CRITICAL)
- Given: MITC4 element rotated rigidly by angle θ (any axis)
- When: `compute_kt_corotational(u_elem, ...)` is called
- Then: K_T_coro = K_T_ref (same eigenvalues as reference undeformed stiffness)

### SC-013 — Same as SC-012 for MITC3

### SC-020 — Patch test MITC4 corotational
- Given: 4-element patch under uniform membrane strain (plane stress)
- When: `assemble_kt_corotational(u)` is called
- Then: equilibrium residual < 1e-10

### SC-021 — Patch test MITC3 corotational
- Same as SC-020 for a MITC3 patch

### SC-030 — Small-angle agreement (θ < 5°)
- Given: MITC4 element with rotation θ = 5°
- When: K_T_coro and K_T_TL are computed
- Then: max element-wise relative diff < 0.1%

### SC-031 — Regression: TL unchanged
- Given: existing test suite runs with `use_corotational_kt=false` (default)
- When: full test suite executes
- Then: all previously passing tests pass, no change in output

### SC-032 — Large-angle benefit (θ = 20°)
- Given: MITC4 element with rotation θ = 20°
- When: K_T_coro and K_T_TL are compared
- Then: max element-wise relative diff > 5% (confirms corotational benefit)

### SC-040 — Config parsing: use_corotational_kt
- Given: YAML with `use_corotational_kt: true`
- When: `RotorConfig.from_yaml()` parses it
- Then: `config.use_corotational_kt == True`

### SC-041 — Config parsing: kt_update_freq
- Given: YAML with `kt_update_freq: 3`
- When: parsed
- Then: `config.kt_update_freq == 3`

### SC-042 — Config defaults
- Given: YAML without `use_corotational_kt` or `kt_update_freq`
- When: parsed
- Then: `use_corotational_kt == False`, `kt_update_freq == 1`
