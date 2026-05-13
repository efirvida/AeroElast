# Design: element-corotational-shell

## Architecture Overview

```
Python layer
  src/aeroelast/core/config.py
    └─ RotorConfig: +use_corotational_kt, +kt_update_freq

Rust layer
  crates/aeroelast-core/src/elements/
    ├─ mitc4.rs    ← SVD polar decomp, build_t24_deformed, compute_kt_corotational
    └─ mitc3.rs    ← build_t18_deformed, compute_kt_corotational

  crates/aeroelast-core/src/assembly/
    └─ assembler.rs ← assemble_kt_corotational(u: &[f64])

  crates/aeroelast-solvers/src/petsc/fsi/
    ├─ rotor_fsi.rs      ← use_corotational_kt, kt_update_freq wiring
    └─ rotor_inertial.rs ← same wiring
```

## Design Decisions

### D-01 — SVD replaces Denman-Beavers
Use `nalgebra::Matrix3::svd(true, true)` for polar decomposition.
R = U·Vᵀ; if det(R) < 0 flip sign of last column of U.
Rationale: exact in one call, no iteration, handles degenerate F.

### D-02 — Frame = deformed mid-plane
The corotational frame is extracted from deformed nodal coordinates (translations only).
No director interpolation. Consistent with Mindlin thin-shell assumption.
Algorithm: g_r = Σ dN/dr · x_def; g_s = Σ dN/ds · x_def; ê₃ = g_r × g_s / ||...||.

### D-03 — T is block-diagonal 3×3
T24 (MITC4) and T18 (MITC3) are assembled as block-diagonal with the same 3×3
rotation R_frame applied to all translation AND rotation triplets.
This is valid under the Mindlin assumption that rotations transform like vectors.

### D-04 — K_σ ≡ K_σ_TL(frame_def)
K_spin term is dropped (see proposal §K_σ^coro). K_σ is the geometric stiffness
evaluated at deformed coordinates, equal to the existing `compute_kg_local()` call
expressed in global frame via Tᵀ K_g T.

### D-05 — Assembler: no new data structures
`assemble_kt_corotational` takes `u: &[f64]` (full DOF vector), slices it per-element,
and calls `element.compute_kt_corotational(u_elem, coords_elem)` inside the existing
rayon parallel loop. Returns `CooMatrix` same as `assemble_kt`.

### D-06 — Rust ↔ Python flag passing
`RotorFsiConfig` (Rust struct in rotor_fsi.rs) gains two fields:
```rust
pub use_corotational_kt: bool,
pub kt_update_freq: usize,
```
Python `RotorConfig.to_dict()` serializes them; `run_rotor_fsi_solver()` PyO3 binding
reads them from the dict.

## Files Changed

| File | Change |
|------|--------|
| `crates/aeroelast-core/src/elements/mitc4.rs` | Add SVD polar_decomp, build_t24_deformed, compute_kt_corotational |
| `crates/aeroelast-core/src/elements/mitc3.rs` | Add build_t18_deformed, compute_kt_corotational |
| `crates/aeroelast-core/src/assembly/assembler.rs` | Add assemble_kt_corotational |
| `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` | Add use_corotational_kt, kt_update_freq; conditional assembly |
| `crates/aeroelast-solvers/src/petsc/fsi/rotor_inertial.rs` | Same as rotor_fsi.rs |
| `src/aeroelast/core/config.py` | RotorConfig: add two fields + to_dict() |
| `crates/aeroelast-core/tests/test_mitc4_corotational.rs` | New: SC-001..SC-013, SC-020, SC-030..SC-032 |
| `crates/aeroelast-core/tests/test_mitc3_corotational.rs` | New: SC-013, SC-021 |
| `tests/test_config_corotational.py` | New: SC-040..SC-042 |

## Key Invariants

1. `use_corotational_kt=false` (default): zero code-path change — same as today.
2. `compute_kt_corotational(u=0)` must return the same matrix as `compute_kt()`.
3. `build_t24_deformed` never panics — fallback to identity frame if cross product degenerates.
4. `assemble_kt_corotational` is thread-safe (rayon, no shared mut state beyond CooMatrix builder).

## Dependency Order

```
D-01 (SVD) → D-02 (frame) → D-03 (T matrix) → D-04 (K_sigma) → D-05 (assembler)
                                                                → D-06 (config flags) → assembler dispatch
```
