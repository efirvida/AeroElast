# Feature: Test-suite physical correctness

## Objective

Make the AeroElast test suite tell the truth about the physics. Two red tests are
test defects over correct code; several green tests cannot fail; and the SHELL
family (MITC3/MITC4, the elements under active validation) has none of the
classical physical invariants. After this feature the suite must fail when the
physics is wrong and pass only when a measured physical invariant holds.

## Baseline (measured, 2026-09-22)

| Check | Result |
| --- | --- |
| `pytest -m "not slow"` | 1 failed / 338 passed / 4 skipped |
| `cargo test -p aeroelast-core` | 72 passed / 1 failed |
| `cargo test -p aeroelast-mesh` | blocked: `hdf5-metno-sys` build script fails |
| `cargo test -p aeroelast-solvers` | blocked: `precice` build script fails |
| `ccx` on PATH | absent -> every CalculiX parity test skips |

Run commands:

```text
~/miniconda3/envs/aeroelast-dev/bin/python -m pytest -m "not slow" -q
cd crates && cargo test -p aeroelast-core -q
```

## Findings that drive the work

### A. The two red tests are test bugs (proved by measurement)

- **A1** `tests/test_mass_matrix_validation.py::test_mass_consistency[tri3]`.
  `compute_physical_mass_from_matrix()` returns `0.75 * trace(M)` and claims
  `tr(M) = 4/3 * m`. Measured on the test's own mesh: `tr(M)/m = 1.500012`
  (tri3) and `1.333344` (quad4); the per-direction translational block sum is
  `m` exactly (`1.000000`) for both. `0.75 * 1.5 = 1.125` is the observed
  12.5%. The constant is Q4-specific and was applied to a triangle.
- **A2** Rust `test_drill_membrane_operator_detects_drill_gradient`. It calls
  `b_md_mitc4_plus(pre, 0.0, 0.0)`. At `(xi, eta) = (0, 0)` every entry of
  `drill_midside_shape_derivatives` is identically zero by construction (each
  carries a factor of `xi` or `eta`), so the operator is the zero matrix and the
  assertion holds for no displacement. The companion zero-test passes for the
  same trivial reason. `b_md_mitc4_plus` is also not wired into any production
  path; only the standalone `b_drill` penalty (mitc4.rs:1126-1140) is.

### B. Validation that cannot fail

- `tests/test_beam_shell_4cases_parity.py::test_linear_static_with_analytical`
  computes `anal_disp`, uses it only in prints, asserts only AE-vs-CCX, and
  skips entirely when CalculiX is absent (this machine).
- `tests/test_shell_validation_fixed.py::test_ratio_physical`: window
  `-50% / +20%`; measured error `0.04%`.
- `tests/test_shell_validation_fixed.py` FX/FY/FZ: `error < 5.0` percent;
  measured `3.4% / 1.2% / 1.9%`.
- `tests/test_ko2017_performance.py::assert_relative_error` hard-fails above 5%
  regardless of its `tol` argument, so every per-test `tol` is decorative.
  Hook Table 14: measured `0.9927` vs reference `1.12` (11.4% off), asserted
  only as `assert norm > 0.1`.
- Rust `test_ke_local_eigenvalues_nonsymmetric` feeds `SymmetricEigen::new(ke)`,
  which reads one triangle of a nonsymmetric input, so its PSD assertion is
  invalid under its own premise.

### C. Missing invariants for the SHELL family

Plane quads have rigid-body-mode tests (`tests/test_quad_elements.py`).
MITC3/MITC4 have no rigid-body-mode test, no K symmetry/PSD test outside the
flat-plate case, no patch test, and no mass test. MITC3 has 7 unit tests total.

### D. Open physics questions

- FX (axial) case is 3.4% above the exact analytical `P*L/(E*A)`: discretization,
  or the clamped edge suppressing Poisson contraction (a real physical effect
  that makes the reference formula wrong)?
- Hook Table 14 is 11.4% off the published reference: mesh (`n_width = 8` vs the
  paper's `N = 16`) or formulation?

## Decisions taken (user, 2026-09-22)

1. Scope: **A + B + C + D** (full).
2. `b_md_mitc4_plus` and its two tests are **deleted** as dead code rather than
   repaired.

## Work units

Each unit ends with one Conventional Commit on `test/physical-correctness`.
Tests and docs travel with the behaviour they describe.

### T1 - Mass tests: honest extraction and exact coefficients (A1)

- Replace `compute_physical_mass_from_matrix` with an element-agnostic extraction
  of the physical mass: the sum of the translational block for each direction,
  which equals `rho*h*A` exactly for any consistent mass matrix.
- Keep the global conservation test for `tri3` and `quad4`; it must pass on
  physics, not on a tuned constant.
- Add exact per-element coefficient tests on a one-element mesh, per
  translational direction:
  - tri3: `M_ii = rho*h*A/6`, `M_ij = rho*h*A/12`;
  - quad4 (rectangle/parallelogram): `M_ii = rho*h*A/9`, adjacent `rho*h*A/18`,
    opposite `rho*h*A/36`.
- Add a row-sum check: the row sum of the consistent mass matrix must equal the
  lumped mass distribution (`rho*h*A/n_nodes` per translational DOF).
- Files: `tests/test_mass_matrix_validation.py`.
- Acceptance: the full mass file passes; the exact coefficients are asserted,
  not approximated.

### T2 - Delete the dead drill-membrane operator (A2)

- Delete `b_md_mitc4_plus`, `drill_midside_shape_derivatives`, and the two tests
  that only exercise them.
- Remove any now-unused helper they were the sole caller of (verify, do not
  assume).
- Files: `crates/aeroelast-core/src/elements/mitc4.rs`.
- Acceptance: `cargo test -p aeroelast-core` is green with 71 tests; `cargo
  build` shows no new dead-code warning for the removed names.

### T3 - Shell K invariants: symmetry, PSD, rigid body, patch tests (C1-C3)

- In `mitc3.rs` and `mitc4.rs` test modules, on `compute_ke_global`:
  - symmetry: `|K - K^T|` below a scaled tolerance;
  - positive semi-definiteness: no eigenvalue below a scaled tolerance;
  - rigid-body modes: for the three translations and the three rotations built
    from the element's own coordinates, `|K * u_rigid|` is below a scaled
    tolerance. Record which modes are penalized if the drill penalty does not
    vanish for rigid in-plane rotation;
  - membrane patch test: a linear in-plane displacement field must reproduce the
    imposed constant strain at every Gauss point, exactly (to round-off);
  - bending patch test: a constant-curvature field must reproduce the imposed
    curvature at every Gauss point, exactly (to round-off).
- Replace `test_ke_local_eigenvalues_nonsymmetric` with a test whose name matches
  what it proves.
- Files: `crates/aeroelast-core/src/elements/mitc3.rs`,
  `crates/aeroelast-core/src/elements/mitc4.rs`.
- Acceptance: `cargo test -p aeroelast-core` green; each new assertion is
  exercised by a mutation check (documented in the task file), not assumed.

### T4 - Shell mass invariants (C4)

- In `mitc3.rs` and `mitc4.rs` test modules, on `compute_me_global` for a single
  element: symmetry, PSD, total translational mass `= rho*h*A` per direction,
  the exact coefficients from T1, and the rotary-inertia block `= rho*h^3*A/12`
  scaled per direction.
- Files: same as T3.
- Acceptance: green Rust suite; the rotational block is asserted, not ignored.

### T5 - Restore the real analytical assertion in the shell parity suite (B1)

- Split `_run_static_case` so the AeroElast solve and the analytical comparison
  happen before, and independently of, CalculiX availability: assert AE vs
  analytical first, then run CCX only when it is present.
- Files: `tests/test_beam_shell_4cases_parity.py`.
- Acceptance: with no `ccx` installed the four static load cases assert against
  the analytical solution instead of skipping; with `ccx` present the CCX
  comparison still runs.

### T6 - Settle the axial 3.4% with a convergence study (D1)

- Measure the FX case error against `P*L/(E*A)` across a mesh sequence.
- If the error decreases with refinement, keep the analytical reference and
  assert the converged error with a justified bound.
- If it plateaus, the clamped-edge Poisson effect is physical: document it, use
  the physically correct reference (or a constrained-edge reference), and assert
  the plateau with a justified bound.
- Files: `tests/test_shell_validation_fixed.py` and, if the study needs its own
  home, `tests/test_shell_convergence.py`.
- Acceptance: a printed convergence table in the test output plus an assertion
  whose bound is derived from the measured physics and stated in a comment.

### T7 - Settle the Hook 11.4% with a convergence study (D2)

- Measure the Hook Table 14 normalized value for `n_width` in a refinement
  sequence up to the paper's `N = 16`.
- Assert the finest value within a bound justified by the measured convergence,
  replacing `assert norm > 0.1`.
- Files: `tests/test_ko2017_performance.py`.
- Acceptance: the assertion can fail; the bound is justified in a comment with
  the measured numbers.

### T8 - Remove the decorative tolerances (B2-B4)

- `assert_relative_error`: delete the hard-coded 5% `pytest.fail` so `tol` is
  authoritative; remove the duplicated `np.isclose(..., rtol=0.05)` call sites.
- `test_ratio_physical`: tighten to a bound the measured 0.04% justifies.
- FX/FY/FZ: replace the flat 5% percentage window with the bounds T6 justifies.
- Files: `tests/test_ko2017_performance.py`, `tests/test_shell_validation_fixed.py`.
- Acceptance: green suite; every remaining tolerance has a stated physical
  justification in a comment or docstring.

## Honest error accounting

| | Before | After (target) |
| --- | --- | --- |
| Python non-slow failures | 1 | **0** |
| Python non-slow skips caused by a missing binary | 4 | 0 for analytical-only paths |
| Rust `aeroelast-core` failures | 1 | **0** |
| Rust `aeroelast-core` tests | 73 | 71 + new invariants |
| Tests whose assertion cannot fail | 5 identified | **0** |

## Out of scope

- `mitc3.rs` / `mitc4.rs` production formulation changes (the nonlinear D1-D7 plan
  in `docs/improvement_plan_mitc4_vs_s4r.md` stays separate).
- Wiring `b_md_mitc4_plus` into K: deleted instead (user decision).
- The vendored NuMAD tree and its lint debt.
- Installing CalculiX; the suite must be honest without it.
