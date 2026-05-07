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
