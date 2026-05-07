# Tasks: Rotor Corotational Solver Consistency Fixes

> **Implementation order**: Fix #3 → Fix #1 → Fix #2 → Fix #5 → Fix #4.
> Design (`design.md §Implementation Order`) puts Fix #3 first because it is a pure
> behaviour-equivalent refactor that establishes a stable baseline before any numerical
> change lands. The invocation prompt lists a different order; the design is the
> authoritative source.
>
> Each task maps to one focused commit and can be independently reverted with
> `git revert <sha>`. Tasks within the same fix group are sequential unless noted.

---

## Fix #3: Unified ω-rebuild policy (refactor, no numerical change)

Covers spec requirements: **Omega-Rebuild Policy for K_G and K_SP** and
**Configurable Omega-Rebuild Thresholds on RotorConfig**.

- [x] 3.1 Add in-crate Rust unit tests for `omega_changed_significantly` predicate
  (`crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs`, new `#[cfg(test)]` block).
  Cover: stable ω → skip; large jump > high threshold → rebuild; ω→0 → rebuild;
  first-call (`omega_sq_at_last` unset) → rebuild; edge case where both ω² < `eps=1e-12`
  → rebuild (guard against divide-by-zero). Tests MUST fail before the predicate exists;
  acceptance: all pass after 3.2.

- [x] 3.2 Extract `omega_changed_significantly(omega_new, omega_last, threshold_rebuild,
  threshold_skip, currently_rebuilt) -> bool` into `rotor_fsi.rs`.
  Semantics (from design §Decision): relative `|Δ(ω²)|/ω²` predicate with hysteresis;
  when `max(ω², ω_last²) < eps` fall back to `true` (rebuild). The existing absolute
  `ksp_omega_threshold` config field is KEPT but rerouted: it seeds the relative test
  only as the `eps` guard, not as the rebuild criterion. Acceptance: unit tests in 3.1
  all pass; `cargo test` green in `aeroelast-solvers`.

- [x] 3.3 Replace K_G rebuild gate (`rotor_fsi.rs:439-457`) with a call to
  `omega_changed_significantly` using the new predicate. No threshold value change —
  defaults remain `0.005`/`0.003`. Acceptance: behaviour-equivalent; existing
  `test_rotor_rust_parity.py` passes unchanged.

- [x] 3.4 Replace K_SP rebuild gate (`rotor_fsi.rs:808`) with a call to the same
  `omega_changed_significantly` predicate. Acceptance: same as 3.3; additionally verify
  with `test_rotor_physical_consistency.py::test_kg_hysteresis_prevents_chattering`
  (already exercises the hysteresis path).

- [x] 3.5 Add `omega_rebuild_rel_high: float` (default `0.005`) and
  `omega_rebuild_rel_low: float` (default `0.003`) to `RotorConfig` in
  `src/aeroelast/core/config.py`. Forward through `RotorConfig.to_dict()` and thread
  through the Rust `RotorFSIConfig` struct in
  `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs`. Acceptance: a YAML without
  these fields uses the defaults; a YAML that sets them overrides the predicate
  thresholds; `test_rotor_rust_parity.py` still passes.

- [x] 3.6 Wire new threshold fields through the PyO3 binding in
  `crates/aeroelast-py/src/lib.rs` so `run_rotor_fsi_solver` accepts and passes them.
  Also thread from the Python call site in `rotor.py` (~`:2150` neighbourhood).
  Acceptance: smoke test — instantiate `RotorConfig` with explicit thresholds and
  confirm the Rust solver receives them (assert via debug log or a new test parameter).

- [x] 3.7 Add `tests/test_omega_rebuild_predicate.py` (Python integration layer).
  Verify rebuild count on a ramped-ω synthetic case matches the expected reduction vs
  the old absolute-threshold policy. Acceptance: rebuild count with new predicate ≤
  rebuild count with old policy on all parametrized ω ramps.

---

## Fix #1: Centrifugal force evaluated at X₀ when K_SP active

Covers spec requirement: **Centrifugal Force Evaluation Coordinates**.

- [ ] 1.1 Rewrite `tests/test_rotor_physical_consistency.py::test_centrifugal_deformed_geometry`
  to reflect the new contract. **Read the test body first** — it is a pure-math test that
  does not call the Rust solver; it computes `F_exact` (deformed) and `F_cached` (reference)
  analytically and asserts the cached method has error ≈ `deformation_ratio`. The minimal
  correct rewrite (option c) is: keep the pure-math structure but flip the parametrized
  assertion so that when `include_ksp = true`, the reference-geometry value is CORRECT
  (zero error expected) and when `include_ksp = false`, the deformed-geometry value is
  correct (current behaviour). Do NOT merge with 1.3 unless you also want to add a real
  solver call here; 1.3 is the definitive Rust-level gate. This and task 1.2 are in the
  same commit.

- [ ] 1.2 (same commit as 1.1) Implement the centrifugal branch in `rotor_fsi.rs:622-648`:
  when `self.config.include_ksp`, pass `&self.all_node_coords` (`X₀`) to
  `compute_centrifugal_force`; otherwise pass the per-step deformed vector as today.
  The Euler block (`:651-677`) is unchanged — no LHS counterpart, deformed coords correct
  in both cases. Acceptance: rewritten test in 1.1 passes; `test_rotor_rust_parity.py`
  passes; `cargo test` green.

- [ ] 1.3 Add `tests/test_rotor_centrifugal_branch.py`: 2-node mass-rotor at fixed ω
  and prescribed u. For `include_ksp = true` verify `|F_cf|` matches `m·ω²·r_0,⊥`
  within 1e-12; for `include_ksp = false` verify it matches `m·ω²·r_⊥(X₀+u)` within
  1e-12; verify the two differ by an amount ∝ `ω²·|u|`. Acceptance: all assertions
  pass; this test is independent of 1.1/1.2 and can run in CI without preCICE.

- [ ] 1.4 Run `test_rotor_physical_consistency.py` in full; run `test_rotor_rust_parity.py`.
  Confirm no new failures beyond the intentional replacement of the old centrifugal
  assertion. Document any numerical shifts observed as the expected delta for the
  regression baseline (≤ 1% RMS displacement shift on `include_ksp=true` cases is
  the acceptance threshold from the spec). Acceptance: written diff is committed to the
  PR as a benchmark annotation.

- [ ] 1.5 Read-only cross-check: skim `docs/rotor_inertial_solver_design.md` to
  confirm whether `LinearDynamicFSIRotorInertialSolver` has an analogous centrifugal
  bug. Flag a follow-up issue if yes; do NOT change that solver in this task. No code
  changes. Acceptance: a single comment in the PR description recording the finding.

---

## Fix #2: Newmark RHS history term for G_cor

Covers spec requirement: **Newmark RHS History Term for Gyroscopic Matrix**.

- [ ] 2.1 Add `tests/test_newmark_coriolis_history.py`: synthetic 2-DOF rotating point
  mass, no elastic stiffness or structural damping; known analytical precession solution.
  Integrate at three successively halved time steps and verify:
  - With the current code: `O(Δt)` error (first-order convergence — the bug).
  - After the fix in 2.2: `O(Δt²)` convergence (trapezoidal).
  This test MUST fail on the convergence-order assertion before 2.2 lands; it is the
  regression gate. Acceptance: all assertions pass after 2.2.

- [ ] 2.2 Modify `dynamic_newmark.rs::refactorize`
  (`crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs:722-748`):
  when `g_cor_vals` is non-empty, build `mat_c_rhs` whose triplets are the union of
  `c_vals` and `g_cor_vals` entries (summed at shared `(row, col)` pairs). Store
  `mat_c_rhs` on the struct. Acceptance: `cargo test` green; `test_newmark_coriolis_history.py`
  shows `O(Δt²)` convergence.

- [ ] 2.3 Update `dynamic_newmark.rs::step()` (`:914-1018`): route the C-history
  MatMult (`C·(a1·u + a4·v + a5·a)`) through `mat_c_rhs` instead of `mat_c`. When
  `mat_c_rhs` is absent (G_cor not used), fall back to `mat_c`. Zero new allocations
  per step. Acceptance: `test_newmark_coriolis_history.py` full pass; convergence
  order is `O(Δt²)`; `test_rotor_rust_parity.py` passes within existing 1e-9 tolerance.

- [ ] 2.4 Run `test_rotor_physical_consistency.py` and `test_rotor_performance_report.py`.
  Record RMS displacement shift attributable to the G_cor RHS fix (separate from Fix #1
  shift already recorded in 1.4). Acceptance: combined shift from Fix #1 + Fix #2 ≤ 1%
  RMS on `include_ksp=true` cases.

---

## Fix #5: Per-window inertial force passthrough via PyO3

Covers spec requirements: **Per-Window Inertial Force Field in Rust Result Struct**.

- [ ] 5.1 Add `applied_inertial_forces: Vec<f64>` to `FsiResult` in
  `crates/aeroelast-solvers/src/petsc/fsi/linear_elastic.rs`. Length = `n_full_dofs`.
  Rust accumulates centrifugal + Coriolis + Euler into an `inertial_scratch: Vec<f64>`
  per sub-iteration (buffer reused, zero per-step alloc) and clones it into `FsiResult`
  on the per-window callback path. Acceptance: `cargo build` green; field accessible
  from Rust tests.

- [ ] 5.2 Expose `applied_inertial_forces` through the PyO3 return tuple in
  `crates/aeroelast-py/src/lib.rs`. The new field MUST be additive: existing callers
  that do not reference it continue to work. No existing field renamed or removed.
  Acceptance: `python -c "import _aeroelast"` succeeds; the tuple's new slot is a
  numpy array or `None` and is accessible by index.

- [ ] 5.3 Add a Python fallback in `rotor.py::_step_cb`: if the passthrough field is
  absent or `None` (old pickled results, old Rust binary), fall back to the legacy
  recomputation block. Acceptance: running `_step_cb` with a mock result that lacks the
  field does not raise; the fallback produces values identical to the legacy path.

- [ ] 5.4 Update `rotor.py::_step_cb` (`:1909-1927`) to consume `applied_inertial_forces`
  from the Rust result instead of recomputing. Remove the legacy recomputation block
  (leave the fallback guard from 5.3 in place). Acceptance: diagnostic CSV columns
  are numerically equivalent within 1e-12 on a recorded checkpoint run compared to
  the pre-change output.

---

## Fix #4: Docstring synchronization

Covers spec requirement: **Coriolis Treatment Docstring in rotor.py**.

- [ ] 4.1 Update `rotor.py:284-289` "Theoretical Limitations & Risks" section:
  replace "explicit lagged Coriolis force on RHS" with a description of the implicit
  `G_cor` LHS placement (`a1·G_cor` contribution to `K_eff`) and the (now confirmed
  and fixed) RHS history term `G_cor·(a1·u_n + a4·v_n + a5·a_n)`. Acceptance: manual
  read confirms no contradiction with the Rust implementation.

- [ ] 4.2 Update the LHS/RHS treatment table at `rotor.py:31-43`: move `F_cor` out of
  "explicit lagged" column and into "implicit via G_cor·v_{n+1}" column. Acceptance:
  the table is internally consistent and matches `dynamic_newmark.rs`'s treatment.

- [ ] 4.3 Read-only cross-check: skim `docs/teoria_formulacion_fsi_rotor.md` and
  `docs/validez_teorica_fsi_rotor_corotational.md` for any sections that still describe
  Coriolis as explicit-lagged. Flag a follow-up note in the PR if yes; do NOT edit the
  docs in this change (out of scope per proposal). Acceptance: a single comment
  recording the finding.

---

## Final Integration

- [ ] I.1 Run the full test suite (excluding known-stale modules per `CLAUDE.md`):
  ```
  python -m pytest tests/ -q --tb=short \
      --ignore=tests/test_blade_mesh.py \
      --ignore=tests/test_rotor_inertial.py
  ```
  Acceptance: zero new failures beyond the intentional replacement in task 1.1.
  Confirm `tests/test_ko2017_performance.py` has the same 8 pre-existing failures and
  no new ones.

- [ ] I.2 Run benchmark suite and capture wall-clock deltas vs pre-change baseline:
  ```
  pytest tests/ -m "benchmark" --tb=short
  ```
  Acceptance: wall-clock time does NOT increase by more than 5% relative to the
  pre-change baseline on the same hardware. Record the delta as a benchmark annotation
  in the PR. Open the HPC follow-up question from `design.md §Open Questions` (verify
  `mat_c_rhs` allocation cost on ~10⁶ DOF mesh) as a non-blocking note.

- [ ] I.3 Update `CLAUDE.md` if any architectural fact changed (e.g. `LinearDynamicFSIRotorCorotationalSolver`
  theoretical limitations, new `RotorConfig` fields, new test files to exclude).
  Acceptance: `CLAUDE.md` is consistent with the post-change codebase.

---

## Dependency Graph

```
3.1 → 3.2 → 3.3 → 3.4 → 3.5 → 3.6 → 3.7
                                        ↓
                                   1.1+1.2 → 1.3 → 1.4 → 1.5
                                                     ↓
                                                  2.1 → 2.2 → 2.3 → 2.4
                                                                     ↓
                                                                  5.1 → 5.2 → 5.3 → 5.4
                                                                                     ↓
                                                                                  4.1 → 4.2 → 4.3
                                                                                               ↓
                                                                                            I.1 → I.2 → I.3
```

Tasks 1.1 and 1.2 are one commit. Tasks 1.3 and 1.4 can run in parallel after 1.2.
All other chains, including 3.1 → 3.2 → 3.3 → 3.4, are strictly sequential.

---

## Acceptance Summary

| Spec Requirement | Covered By | Gate |
|-----------------|------------|------|
| Centrifugal at X₀ when K_SP on | 1.1+1.2, 1.3 | `test_rotor_centrifugal_branch.py` all pass |
| Centrifugal at X₀+u when K_SP off | 1.1+1.2, 1.3 | same |
| Regression: tip displacement ≤ 1% RMS shift | 1.4, 2.4 | benchmark annotation in PR |
| Newmark G_cor history term (O(Δt²)) | 2.1, 2.2, 2.3 | `test_newmark_coriolis_history.py` convergence-order assertion |
| ω-rebuild policy coherence (K_G and K_SP) | 3.1-3.4 | `test_omega_rebuild_predicate.py`, parity tests |
| Configurable ω-rebuild thresholds on RotorConfig | 3.5-3.6 | smoke test, parity tests |
| Coriolis docstring aligned with implementation | 4.1-4.2 | manual review |
| Per-window inertial force passthrough (additive) | 5.1-5.4 | bit-for-bit CSV equivalence within 1e-12 |
| Validated benchmarks must not regress | I.1 | full test suite, zero new failures |
| Wall-clock ≤ +5% vs baseline | I.2 | benchmark annotation |
