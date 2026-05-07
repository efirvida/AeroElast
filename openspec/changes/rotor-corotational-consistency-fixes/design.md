# Design: Rotor Corotational Solver Consistency Fixes

## Finding #2 verification: CONFIRMED

`dynamic_newmark.rs::refactorize` (`crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs:722-748`) adds `a1·G_cor` to `K_eff` on the LHS via `MatSetValues(... ADD_VALUES)`. The author's intent is explicit at `:725` — `let g_scale = a1; // γ/(β·dt) — same coefficient as C in K_eff`. However `step()` (`:914-1018`) builds RHS as `F_ext + M·(a0·u + a2·v + a3·a) + C·(a1·u + a4·v + a5·a)` only; `g_cor_*` fields are read solely in `refactorize()`. The implicit Newmark substitution `v_{n+1} = a1·(u_{n+1}-u_n) - a4·v_n - a5·a_n` makes `a1·G_cor` on the LHS inseparable from a `G_cor·(a1·u_n + a4·v_n + a5·a_n)` RHS history term. The current code is solving a non-physical equation that lags Coriolis by exactly the missing history contribution — not a stable antisymmetric damper, just a partial implicitization. Fix #2 is therefore in scope as written in the proposal.

## Technical Approach

Five surgical, independently revertable edits restore correctness without touching the KSP cache, scratch buffers, Rodrigues transforms, or the symmetric-LHS contract.

## Architecture Decisions

### Decision: Centrifugal coordinate selection lives at the call site
**Choice**: Branch on `self.config.include_ksp` in `rotor_fsi.rs::run` and pass either `&self.all_node_coords` (`X_0`) or the per-step deformed vector to the existing `compute_centrifugal_force` helper.
**Alternatives considered**: (a) make `compute_centrifugal_force` itself K_SP-aware; (b) precompute one branch into a struct field.
**Rationale**: The helper is a pure $\omega \times (\omega \times r)$ kernel; pushing policy into it leaks coupling. The call-site branch is two lines, allocates nothing in the `include_ksp=true` path (skips the per-step `deformed_coords` Vec build at `:626-640`), and is trivial to revert. *Caveat*: the Euler block at `:651-677` keeps its own deformed-coord allocation — this is not a bug (no LHS counterpart) and is out of scope.

### Decision: G_cor history term enters via merged `mat_c_plus_gcor` built in `refactorize`
**Choice**: When `g_cor_vals` is non-empty, allocate a `mat_c_rhs` PETSc matrix at refactorization time whose triplets are the union of `c_vals` and `g_cor_vals` at the matching `(row, col)` pairs. `step()` then performs exactly one MatMult on this combined operator with the existing `(a1·u + a4·v + a5·a)` work vector.
**Alternatives considered**: (a) add a second MatMult on `g_cor_vals` inside `step()`; (b) augment `mat_c` itself in place.
**Rationale**: G_cor's sparsity pattern is a superset of C's (translational 3-DOF cross-coupling per node vs Rayleigh damping which inherits K's pattern). Using `mat_c` directly would require expanding K's sparsity, which ripples through `assemble_keff`. A separate `mat_c_rhs` avoids touching the LHS pipeline and keeps `step()`'s hot path at one MatMult — zero new allocations per step. The helper `matvec_add` (`:375`) is already in the right shape; we extend it conceptually by swapping the matrix only.

### Decision: Single relative-$\omega^2$ rebuild predicate, two thresholds promoted to config
**Choice**: Extract `omega_changed_significantly(omega_new, omega_last, threshold_rebuild, threshold_skip, currently_rebuilt) -> bool` into `rotor_fsi.rs`. Use it for both K_G (`:439-457`) and K_SP (`:808`). Promote `THRESHOLD_REBUILD = 0.005` and `THRESHOLD_SKIP = 0.003` to two new `RotorFsiConfig` fields `omega_rebuild_rel_high` / `omega_rebuild_rel_low` with defaults matching the current K_G constants. The K_SP magic absolute `ksp_omega_threshold` is kept as a config field but rerouted: it now seeds the relative test only when both $\omega^2$ values are below `eps = 1e-12`.
**Alternatives considered**: (a) separate fields per matrix; (b) keep absolute Δω.
**Rationale**: K_SP $= -\omega^2 M(I-\hat n \hat n^T)$ and K_G's centrifugal prestress both scale as $\omega^2$, so spin-softening sensitivity does not differ from stress-stiffening sensitivity. Two thresholds (rebuild high / skip low) are kept to preserve the existing hysteresis. **First-call semantics**: when `omega_sq_at_last_*` is not yet set OR `max(ω², ω_last²) < eps`, the predicate returns `true` (rebuild). This matches the current K_G behaviour at startup and avoids a divide-by-zero in the relative form.

### Decision: Per-window inertial passthrough — single summed vector
**Choice**: Add `applied_inertial_forces: Vec<f64>` to `FsiResult` (length = full-DOF, matching `n_full_dofs`). Rust accumulates centrifugal + Coriolis + Euler into a single struct-resident `inertial_scratch: Vec<f64>` during the sub-iteration loop (zero per-step alloc — buffer reused) and clones it into the result on the per-window callback path. PyO3 binding gains an additional positional return slot.
**Alternatives considered**: (a) three separate vectors (centrifugal/Coriolis/Euler); (b) leave `_step_cb` recomputing.
**Rationale**: The Python `_step_cb` (`rotor.py:1909-1927`) only needs the sum to feed `_compute_axis_torque(... f_inertial_nodes_local ...)` and the diagnostic CSV columns — splitting offers no caller benefit and would require the spec's bit-equivalence guarantee on three separate vectors. The Python side defaults the new field to `None`; if absent, `_step_cb` falls back to the legacy recomputation path (additive, non-breaking). Diagnostic CSV column order is therefore preserved exactly.

### Decision: Docstring rewrite is purely textual
**Choice**: Rewrite bullet "1. Explicit Coriolis Force" in `rotor.py:284-289` to describe implicit $G_{cor}$ on the LHS plus the (then-fixed) RHS history term; update the LHS/RHS treatment table at `rotor.py:31-43` to move `{F_cor}` out of "explicit lagged" and into "implicit via $G_{cor}\cdot v_{n+1}$".

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` | Modify | Centrifugal call branches on `include_ksp`; extract `omega_changed_significantly`; replace K_G (`:439-457`) and K_SP (`:808`) gates with the unified predicate; accumulate into `inertial_scratch` and copy into `FsiResult` per converged window. |
| `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs` | Modify | In `refactorize()`, build a `mat_c_rhs` when `g_cor_vals` is non-empty and store it on the struct. In `step()`, route the C-history MatMult through `mat_c_rhs` (falls back to `mat_c` when G_cor is absent). |
| `crates/aeroelast-solvers/src/petsc/fsi/linear_elastic.rs` | Modify | Add `applied_inertial_forces: Vec<f64>` to `FsiResult`. |
| `crates/aeroelast-py/src/lib.rs` | Modify | Expose new field through the PyO3 return tuple; thread `omega_rebuild_rel_high/_low` through `run_rotor_fsi_solver`. |
| `src/aeroelast/core/config.py` | Modify | Add `omega_rebuild_rel_high` (default `0.005`), `omega_rebuild_rel_low` (default `0.003`) to `RotorConfig.to_dict()`. |
| `src/aeroelast/solvers/fsi/rotor.py` | Modify | Replace `_step_cb` recomputation block (`:1909-1927`) with a passthrough read; rewrite docstring sections at `:31-43` and `:284-289`; thread new config keys into the `run_rotor_fsi_solver` call (`:2150` neighbourhood). |
| `tests/test_rotor_centrifugal_branch.py` | New | 2-node mass-rotor at high $\omega$ asserts $|F_{cf}|$ at `X_0` when `include_ksp=true` and at `X_0+u` when `False`. |
| `tests/test_newmark_coriolis_history.py` | New | Synthetic Coriolis-only system, exact ODE benchmark, $O(\Delta t^2)$ convergence. |
| `tests/test_omega_rebuild_predicate.py` | New | Pure-Rust unit (in-crate `#[cfg(test)]`) covering hysteresis, ω=0, and first-call branches. |

## Data Flow

```
preCICE forces (global)
   │
   ▼
forces_to_rotating(θ)
   │
   ├──[include_ksp]──┐
   │                 ▼
   │         compute_centrifugal_force(X₀, ...)        ← Fix #1: NEW path
   │                 │
   │  [!include_ksp] │
   │                 ▼
   │         compute_centrifugal_force(X₀+u, ...)      ← legacy path
   │                 │
   ▼                 ▼
   F_red ←── scatter(centrifugal) ←── scatter(coriolis) ←── scatter(euler) ←── scatter(gravity)
                                          │
                                          ▼
                              inertial_scratch[i] += each contribution    ← Fix #5
                                          │
                                          ▼
   stepper.step(F_red, dt)
        │  in refactorize: K_eff += a1·G_cor      (existing)
        │                  mat_c_rhs := mat_c ⊕ G_cor                      ← Fix #2
        │  in step():     RHS += mat_c_rhs · (a1·u + a4·v + a5·a)          ← Fix #2
        ▼
   converged window:
        FsiResult.applied_inertial_forces ← inertial_scratch.clone()        ← Fix #5
        omega_changed_significantly(...) gates K_G and K_SP rebuilds        ← Fix #3
```

## Testing Strategy

| Layer | What to Test | Approach |
|-------|--------------|----------|
| Unit (Rust) | Centrifugal at $X_0$ vs $X_0+u$ | 2-node mass-rotor, fix $\omega$ and $u$, assert $\|F_{cf}\|$ matches the expected branch within 1e-12. |
| Unit (Rust) | $\omega$-rebuild predicate | Table-driven: stable $\omega$ → skip; jump $> $ high → rebuild; ω→0 → rebuild; first call → rebuild. |
| Unit (Rust) | Newmark Coriolis-only | 2-DOF rotating point mass, no elastic stiffness; compare against analytic precession; assert $O(\Delta t^2)$ convergence by halving $\Delta t$. |
| Integration | Rust↔Python parity | `tests/test_rotor_rust_parity.py` — must pass unchanged; the parity tolerance was already 1e-9, so the new G_cor RHS term cannot drift Python and Rust apart. |
| Integration | Physical consistency | `tests/test_rotor_physical_consistency.py` — assertions are about energy/torque balance; should tighten slightly because the spurious $\omega^2 u$ ghost force is gone. |
| Integration | Rotor performance | `tests/test_rotor_performance_report.py` — track $C_T, C_P$ at converged $\omega$; expected RMS shift $\le 1\%$ on `include_ksp=true` cases; document the new value. |
| Reference | Ko et al. (2017) | `tests/test_ko2017_performance.py` — 8 pre-existing failures per `CLAUDE.md` are NOT regressions; this change must not change their failure modes (no new failures, no spurious passes). |
| Excluded | Rotor inertial / blade-mesh | `tests/test_rotor_inertial*.py` and `tests/test_blade_mesh.py` are excluded per `CLAUDE.md`; do not add new dependencies on them. |

## Implementation Order

1. **Fix #3** (predicate + config plumbing) — refactor only, behaviour-equivalent under `0.005`/`0.003` defaults; lands first as a safe baseline.
2. **Fix #1** (centrifugal at $X_0$) — smallest surgical change, biggest correctness payoff. Validate against the new branch unit test, then re-run `test_rotor_physical_consistency.py`.
3. **Fix #2** (G_cor RHS history) — touches `dynamic_newmark.rs::refactorize` and `step()`. Validate against the new Coriolis-only convergence test.
4. **Fix #5** (inertial passthrough) — extends `FsiResult` and the PyO3 return tuple; no numerical change. Confirm bit-for-bit equivalence on a recorded checkpoint.
5. **Fix #4** (docstring) — last, after all numerical behaviour is locked in.

Each fix is a single commit; `git revert <sha>` independently rolls back any one of them without touching the others. Risk #1 in the proposal (validated cases shift) is bounded by the integration-test suite above; the acceptable threshold is **interface-displacement RMS shift $\le 1\%$ at the converged FSI state with `include_ksp=true`**, with the shift directionally consistent with removing the spurious extra $\omega^2 u$ contribution (i.e. larger steady tip displacement, lower steady $C_P$ at the same $\omega$).

## Open Questions

- [ ] Should the new `omega_rebuild_rel_high/_low` fields appear in the published rotor performance docs (`docs/mejoras_rendimiento_fsi_rotor_corotational.md`)? Out of this change's scope but worth a follow-up note.
- [ ] Verify with one HPC test that `mat_c_rhs` allocation in `refactorize` does not measurably regress the rebuild path on a large mesh (~10⁶ DOF). Refactorizations are rare (ω-threshold gated), so the worst case is negligible — but record the number on a real case.
