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
