# Proposal: Rotor Corotational Solver Consistency Fixes

## Intent

`LinearDynamicFSIRotorCorotationalSolver` is the production rotor FSI solver in this codebase. A fresh code review surfaced five inconsistencies between the documented theoretical contract (Python solver header) and the Rust hot-loop implementation. The most severe is a **silent double-count of the spin-softening term**: centrifugal force is evaluated at deformed coordinates $X_0 + u$ while $K_{SP}$ is simultaneously assembled on the LHS, producing an $O(\omega^{2}\,|u|)$ error in the equilibrium residual that grows precisely in the high-RPM, large-displacement regime that motivates the rotor solver in the first place. The remaining findings are a likely Newmark RHS asymmetry, an inconsistent $\omega$-rebuild policy across $K_G$ and $K_{SP}$, a stale docstring, and an inertial-force recomputation in the per-window Python callback. This change addresses them as a single coherent pass to restore correctness and remove implementation/documentation drift.

## Scope

### In Scope
- **Centrifugal/$K_{SP}$ coherence (CRITICAL)**: Branch the centrifugal evaluation in `rotor_fsi.rs:622-648` on `config.include_ksp`. Use $X_0$ (undeformed) when `include_ksp = true`; only use $X_0 + u$ when $K_{SP}$ is disabled. The Euler block at `:651-677` is left as-is (no LHS counterpart exists, evaluation at deformed coords is correct per the design contract).
- **Newmark RHS consistency for $G_{cor}$ (HIGH, conditional on design-phase verification)**: `dynamic_newmark.rs::refactorize` adds $a_1\,G_{cor}$ to the LHS, but `step()` does not appear to mirror the corresponding $G_{cor}\cdot(a_1 u_n + a_4 v_n + a_5 a_n)$ history term on the RHS as it does for $C$. Design phase will line-verify and, if confirmed, add the missing term.
- **Unified $\omega$-rebuild policy (MEDIUM)**: Replace the absolute $\Delta\omega$ test for $K_{SP}$ rebuild (`rotor_fsi.rs:800-813`) with the same relative $\Delta(\omega^2)/\omega^2$ predicate already used for $K_G$ (`:425-498`). Promote the hysteresis thresholds (currently magic numbers `0.005`/`0.003`) to fields on `RotorFSIConfig`.
- **Docstring alignment (LOW)**: Update `rotor.py:281-313` "Theoretical Limitations & Risks" so the Coriolis treatment matches the Rust implementation (implicit $G_{cor}$ on LHS, no longer "explicit lagged").
- **Per-window inertial-force reuse (EFFICIENCY)**: Have the Rust per-window result struct return the last-applied inertial force vector; replace the recomputation in `rotor.py::_step_cb` (`:1909-1927`) with a direct read.
- **Regression coverage**: Add a targeted test that toggles `include_ksp` and asserts the centrifugal branch selects the correct coordinates; protect the existing validated rotor cases (blade-only static, dynamic FSI) against regression.

### Out of Scope
- Non-symmetric solver paths (full antisymmetric $[G]$ treatment beyond the implicit symmetric-on-LHS form already in place).
- Full geometric nonlinearity (corotational element kinematics, large-rotation update of $K$).
- Migration to consistent (non-lumped) mass.
- Replacing the constant/ramped $\omega$ providers with a higher-order ($\textit{e.g.}$ RK4) integrator.
- Any change to the KSP factorization cache, scratch buffers, or Rodrigues helpers — those are correct and load-bearing.

The "Theoretical Limitations & Risks" section in `rotor.py` defines the boundary between this change and a future nonlinear/non-symmetric formulation; we stay strictly inside it.

## Capabilities

### New Capabilities
- None.

### Modified Capabilities
- `rotor-corotational-solver`: hot-loop centrifugal evaluation, Newmark RHS history term for $G_{cor}$, unified $\omega$-rebuild policy, per-window inertial-force passthrough, and aligned docstring.

## Approach

1. **Centrifugal branch**: Pass `self.all_node_coords` ($X_0$) into `compute_centrifugal_force` when `self.config.include_ksp` is true; pass the existing deformed vector only when $K_{SP}$ is off. Single branch, no allocation in the $K_{SP}$-on path (avoids the per-step `deformed_coords` build).
2. **Newmark RHS**: After design-phase line verification of `step()`, if $G_{cor}$'s history term is indeed missing, mirror the $C$ contribution exactly: add $G_{cor}\cdot(a_1 u_n + a_4 v_n + a_5 a_n)$ using the same scratch buffers.
3. **$\omega$ predicate**: Extract a single `omega_changed_significantly(omega_prev, omega_new)` helper; use it for both $K_G$ and $K_{SP}$ rebuild gates. Thresholds become `RotorFSIConfig::omega_rebuild_rel_high` / `_rel_low`.
4. **Inertial force passthrough**: Extend the Rust per-window result struct with `last_inertial_force: Vec<f64>` (already computed each sub-iteration); expose via PyO3 binding; consume in `_step_cb`.
5. **Docstring**: Rewrite the Coriolis bullet to describe the implicit $G_{cor}$ LHS placement; cross-reference the (then-fixed) RHS history term.

Detailed numerical formulation, signatures, and call-site changes are deferred to the design phase.

## Affected Areas

| Area | Impact | Description |
|------|--------|-------------|
| `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` | Modified | Centrifugal branch on `include_ksp`; unified $\omega$-rebuild predicate; per-window result struct extension. |
| `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs` | Modified (conditional) | Add $G_{cor}$ history term to RHS in `step()` if confirmed missing. |
| `crates/aeroelast-py/` | Modified | Expose new field on the per-window result struct. |
| `src/aeroelast/solvers/fsi/rotor.py` | Modified | Update theoretical-limitations docstring; replace inertial-force recomputation in `_step_cb` with passthrough read. |
| `src/aeroelast/core/config.py` | Modified | Add `omega_rebuild_rel_high`/`_rel_low` to `RotorConfig` (with defaults matching current $K_G$ thresholds). |
| `tests/` | New | Regression test for the centrifugal $X_0$-vs-deformed branch under both `include_ksp` settings. |

## Risks

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Fix #1 changes computed forces enough to shift validated FSI results | Medium | Re-run validated cases (blade-only static, dynamic FSI benchmarks) with `include_ksp=true`; expected change is within numerical tolerance because the previous formulation was inconsistent, not "calibrated". |
| Finding #2 turns out to be a misread of `step()` | Medium | Design phase MUST line-verify `step()` end-to-end before any edit; if the term is already there under a different form, reduce this finding to a docstring/comment fix. |
| Unified $\omega$-rebuild thresholds trigger more refactorizations on some cases | Low | Defaults match current $K_G$ behaviour; users can override via `RotorConfig`; benchmark before merging. |
| Per-window inertial passthrough exposes ABI mismatch in PyO3 binding | Low | Keep field optional / additive; old callers ignore it. |

## Rollback Plan

Each fix is independent and gated behind a small surface (centrifugal branch, RHS term, predicate helper, struct field, docstring). Revert per-finding via `git revert` on the corresponding commit; production behaviour returns to the pre-change state. No schema or persisted artifact changes.

## Dependencies

- Existing fresh code review on this branch (findings 1–5 above). No new exploration required.
- `docs/teoria_formulacion_fsi_rotor.md` and `docs/validez_teorica_fsi_rotor_corotational.md` for the theoretical contract.

## Success Criteria

- [ ] `rotor_fsi.rs` centrifugal call selects $X_0$ when `include_ksp=true`, $X_0 + u$ otherwise; covered by a regression test.
- [ ] Newmark RHS in `step()` is line-verified; $G_{cor}$ history term present (added if missing, documented if already there).
- [ ] $K_G$ and $K_{SP}$ rebuild on the same relative $\omega^2$ predicate, with thresholds exposed on `RotorConfig`.
- [ ] `rotor.py` "Theoretical Limitations & Risks" reflects the implicit Coriolis treatment.
- [ ] `_step_cb` reads inertial force from the Rust per-window result struct; the recomputation block is removed.
- [ ] Validated rotor benchmarks (blade-only static, dynamic FSI) show no regression beyond expected numerical drift from fix #1.
