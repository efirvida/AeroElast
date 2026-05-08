# Apply Progress: rotor-corotational-consistency-fixes (Fix #3 batch)

## Status: done

## Tasks completed

- [x] 3.1 Rust unit tests for `omega_changed_significantly` predicate (6 tests)
- [x] 3.2 Extract predicate into `rotor_fsi.rs` with `eps` parameter seeded from `ksp_omega_threshold²`
- [x] 3.3 Replace K_G rebuild gate with unified predicate call
- [x] 3.4 Replace K_SP rebuild gate with unified predicate call
- [x] 3.5 Add `omega_rebuild_rel_high`/`omega_rebuild_rel_low` to `RotorConfig`
- [x] 3.6 Wire new fields through PyO3 binding (`lib.rs`) and Python call site (`rotor.py`)
- [x] 3.7 Add `tests/test_omega_rebuild_predicate.py` — 15 passed, 1 skipped (env)

## Files changed

- `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs`
  - `RotorFsiConfig`: added `omega_rebuild_rel_high`, `omega_rebuild_rel_low`; doc-updated `ksp_omega_threshold`
  - Added `omega_changed_significantly(omega_new, omega_sq_at_last, threshold_rebuild, threshold_skip, currently_rebuilt, eps)` free function
  - `omega_sq_at_last_ksp: f64` field (was `omega_at_last_ksp`), initialized to `NEG_INFINITY`
  - `last_ksp_rebuild_step: usize` field
  - K_G gate: replaced hardcoded relative check with `omega_changed_significantly` + `ksp_omega_threshold²` eps
  - K_SP gate: replaced absolute `|Δω| > threshold` with `omega_changed_significantly` + `ksp_omega_threshold²` eps
  - `apply_ksp()`: stores `omega_sq_at_last_ksp = omega * omega`
  - `dummy_config()`: added new fields
  - 6 new `#[cfg(test)]` predicate tests

- `crates/aeroelast-py/src/lib.rs`
  - `run_rotor_fsi_solver`: added `omega_rebuild_rel_high: f64`, `omega_rebuild_rel_low: f64` params

- `src/aeroelast/core/config.py`
  - `RotorConfig`: added `omega_rebuild_rel_high: float = 0.005`, `omega_rebuild_rel_low: float = 0.003`
  - `to_dict()`: forwarded both fields

- `src/aeroelast/solvers/fsi/rotor.py`
  - `_init_rotor_config()`: stored `_omega_rebuild_rel_high`, `_omega_rebuild_rel_low`
  - `_aeroelast.run_rotor_fsi_solver()` call: passes both new kwargs

- `tests/test_rotor_rust_parity.py`
  - `_call_binding()` and `test_mismatched_dofs_raises`: added two new positional args

- `tests/test_omega_rebuild_predicate.py` (new)
  - Pure-Python replica of predicate + 15 tests

## Key design decision resolved (task 3.2)

`ksp_omega_threshold` is rerouted as the `eps` guard in `omega_changed_significantly`.
Both K_G and K_SP call sites pass `eps = ksp_omega_threshold * ksp_omega_threshold`.
The Python replica mirrors this: `eps: float = 1e-8` (= (1e-4)²).

## Behavior change (documented risk)

K_SP first-call behavior changed: old code initialized `omega_at_last_ksp = initial_omega`
so constant-ω cases never applied K_SP. New code initializes to `NEG_INFINITY` so K_SP
is applied once at startup even for constant ω. This is a latent bug fix — K_SP was never
being applied when ω started at its steady-state value.

## Test verification

- Rust: `cargo check -p aeroelast-solvers` → Finished dev (9 warnings, no errors)
- Rust: `cargo check -p aeroelast-py` → Finished dev (8 warnings, no errors)
- Rust tests: `cargo test -- omega_pred` → 6 passed, 0 failed
- Python: `pytest tests/test_omega_rebuild_predicate.py` → 15 passed, 1 skipped
- Python parity: `pytest tests/test_rotor_rust_parity.py` → 15 passed, 8 pre-existing failures (TestMapOmegaProvider), 30 skipped

---

# Apply Progress: rotor-corotational-consistency-fixes (Fix #1 batch)

## Status: done

## Tasks completed

- [x] 1.1 Rewrote `tests/test_rotor_physical_consistency.py::test_centrifugal_deformed_geometry`:
  added `include_ksp` to `@pytest.mark.parametrize`; when `include_ksp=True`, the reference-geometry
  (X₀) result is CORRECT — asserted to have zero relative error vs the true X₀ force, and the
  deformed-geometry result is WRONG (error proportional to `deformation_ratio`); when
  `include_ksp=False`, the deformed-geometry result is CORRECT (current behaviour preserved).
  Module docstring updated to reflect the corrected contract.

- [x] 1.2 Implemented centrifugal coordinate branch in
  `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` (previously around line 696-722, now
  lines 695-741 after edit): the `if self.config.include_centrifugal` block now branches on
  `self.config.include_ksp`. When true, `compute_centrifugal_force` is called with
  `&self.all_node_coords` (X₀), skipping the per-step `deformed_coords` Vec build entirely.
  When false, the existing deformed-coords build and call is retained verbatim. The Euler block
  immediately following is untouched per design. Comment updated to explain the double-count risk.

- [x] 1.3 Created `tests/test_rotor_centrifugal_branch.py`: 2-node mass-rotor (nodes at R=10 m on
  X-axis and at 53° off-axis) with prescribed displacements. Four test functions validate:
  (a) include_ksp=True → F_cf at X₀ matches m·ω²·r₀ within 1e-14; displacement contribution absent;
  (b) include_ksp=False → F_cf at X₀+u matches m·ω²·r_⊥(X₀+u) within 1e-14;
  (c) branch diff equals m·ω²·|u_perp| within 1e-12 for all ω and u combinations;
  (d) zero ω gives zero force in both branches. 90 parametrized cases, all passed.
  Pure-math, no _aeroelast import required.

- [x] 1.4 Test runs completed:
  - `test_rotor_physical_consistency.py`: 36 passed, 1 pre-existing failure
    (`test_kg_hysteresis_prevents_chattering` fails on baseline too — confirmed by stash check),
    1 skipped. No new failures.
  - `test_rotor_rust_parity.py`: 8 pre-existing failures (TestMapOmegaProvider — same as baseline),
    29 skipped. No new failures.
  - Rust binary not rebuilt; Rust change verified by inspection only — the binary in venv
    reflects the pre-Fix#1 code. Numerical validation of the branch at runtime is deferred to
    the orchestrator's rebuild step.
  - No RMS displacement shift measurable in this environment (Rust not rebuilt). Spec tolerance
    of ≤1% RMS shift on include_ksp=True cases to be validated after orchestrator rebuilds.

- [x] 1.5 Read-only check of `docs/rotor_inertial_solver_design.md`: the inertial solver design
  explicitly states the first version will NOT include a K_G equivalent and does not implement
  centrifugal force directly (`prestress centrifugo` section, line 193). There is therefore NO
  analogous include_ksp centrifugal double-count bug in `LinearDynamicFSIRotorInertialSolver` —
  the K_SP mechanism is absent from that solver by design. No follow-up issue required.

## Files changed

- `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs`
  - Centrifugal force block: added `include_ksp` branch; X₀ path skips deformed-coords Vec build
  - Old comment "computed at DEFORMED geometry" replaced with accurate two-branch explanation

- `tests/test_rotor_physical_consistency.py`
  - `test_centrifugal_deformed_geometry`: added `include_ksp` param, flipped assertions per contract
  - Module docstring: updated item 1 to describe the correct two-branch contract

- `tests/test_rotor_centrifugal_branch.py` (new)
  - 4 test functions, 90 parametrized cases, pure-math, no Rust import needed

## Pre-existing failures confirmed (not introduced by Fix #1)

- `test_kg_hysteresis_prevents_chattering`: fails on baseline (stash-verified)
- `TestMapOmegaProvider::test_constant_omega` and 8 others: fail on baseline (stash-verified)
- All `test_rotor_rust_parity.py::TestRotorAutoInertia` tests: SKIP (no Rust binary in this env)

---

# Apply Progress: rotor-corotational-consistency-fixes (Fix #2 batch)

## Status: done

## Tasks completed

- [x] 2.1 Added `test_gcor_rhs_history_secondorder` Rust unit test inside `dynamic_newmark.rs`
  `#[cfg(test)]` block. NOTE: task specified `tests/test_newmark_coriolis_history.py` but was
  implemented as a Rust in-crate test per the invocation prompt directive (PyO3 extension not
  rebuildable in this environment; Python test would exercise the OLD binary). The Rust test
  uses the public API (`update_spin_softening_and_gyroscopic`, `step`) and measures convergence
  order via L2 error at three halved dt values. CONFIRMED FAIL before fix (slope = -1.153,
  err_coarse=7.318e-1 but err_medium=1.627e0 — error grows as dt shrinks, confirming
  unbounded drift from missing G_cor history); PASS after fix (slope_cm=1.9996, slope_mf=1.9999).

- [x] 2.2 Modified `refactorize()` in `dynamic_newmark.rs`:
  - Added `assemble_union_aij()` helper function (assembles non-symmetric AIJ matrix
    from two COO sets via `ADD_VALUES`; does NOT set MAT_SYMMETRIC).
  - Added `MATAIJ_STR` constant (duplicated from `assembler.rs` for local use).
  - Added `mat_c_rhs: Option<PetscMat>` field to `NewmarkStepper` struct.
  - In `refactorize()`: when `g_cor_vals` is non-empty, builds `mat_c_rhs = C ⊕ G_cor`
    by concatenating COO triplets and calling `assemble_union_aij`. When empty, sets
    `mat_c_rhs = None`. Initialized to `None` in `new()`.

- [x] 2.3 Modified `step()` in `dynamic_newmark.rs`:
  - The C-history MatMult now uses `self.mat_c_rhs.as_ref().unwrap_or(&self.mat_c)`.
  - When G_cor is absent: `mat_c_rhs` is `None`, falls back to `mat_c` — zero-cost
    for all non-rotor callers (LinearDynamicSolver, FSI base, etc.).
  - When G_cor is active: uses the pre-built `mat_c_rhs` — single MatMult, zero new
    allocations per step.

- [x] 2.4 Python test runs:
  - `test_rotor_physical_consistency.py`: same pre-existing 1 failure, no new failures.
  - `test_rotor_rust_parity.py`: same 9 pre-existing failures, no new failures.
  - NOTE: Python tests exercise OLD binary (_aeroelast not rebuildable). Numerical
    validation at runtime deferred to orchestrator's rebuild step. Convergence-order
    gate is the Rust unit test (2.1), which passed cleanly.

## Files changed

- `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs`
  - Added `MATAIJ_STR` constant
  - Added `assemble_union_aij()` free function (non-symmetric union of two COO sets)
  - `NewmarkStepper` struct: added `mat_c_rhs: Option<PetscMat>` field
  - `new()`: added `mat_c_rhs: None` initialization
  - `refactorize()`: builds `mat_c_rhs` when `g_cor_vals` non-empty, else `None`
  - `step()`: routes C-history MatMult through `mat_c_rhs` (falls back to `mat_c`)
  - `#[cfg(test)]`: added `make_2dof_stepper()`, `run_coriolis_and_get_error()`,
    `test_gcor_rhs_history_secondorder()` in existing tests module

## Test verification

- Rust: 62 library unit tests — 62 passed, 0 failed (sequential: `--test-threads=1`)
- Rust convergence gate: `test_gcor_rhs_history_secondorder` — PASS (slope_cm=1.9996, slope_mf=1.9999; measured empirically, not assumed)
- `fsi_mock_loop.rs` integration tests: 4 pre-existing failures (user's in-progress inertial
  solver work, NOT touched by Fix #2)
- Python rotor tests: 9 pre-existing failures, no new failures

## Analytical solution verification

System: M·ü + G_cor·u̇ = 0, M=I, G_cor=[[0,-1],[1,0]]
Reduces to: ṗ₁=p₂, ṗ₂=-p₁ with p(0)=[1,0]
Solution: p₁(t)=cos(t), p₂(t)=-sin(t)
Displacement: u₁(t)=sin(t), u₂(t)=cos(t)-1
Initial acceleration: a₀ = M⁻¹·(-G_cor·v₀) = [0,-1]

## Risk note

`mat_c_rhs` is built as a non-symmetric matrix. PETSc LU handles non-symmetric matrices
correctly. The symmetric MAT_SYMMETRIC hint is intentionally omitted from `assemble_union_aij`.
The field is `Option<PetscMat>`, so non-rotor callers (where `g_cor_vals` is always empty)
never allocate or use it — zero overhead on the existing code path.
