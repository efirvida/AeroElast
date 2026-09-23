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
## Progress log

| Task | State | Evidence |
| --- | --- | --- |
| T1 mass tests | **done** | commit `55d65e7`; 15 tests in the file pass (was 9 with 1 red); non-slow suite 345 passed / 4 skipped / 0 failed; liveness proved by two mutations |
| T2 delete dead drill operator | **done** | commit `4a46af2`; `cargo test -p aeroelast-core` 73/1-failed -> 71/0-failed; no reference to either symbol remains |
| T3a MITC4 invariants | **done** | commit `785a48d`; 71 -> 76 Rust tests, 0 failed; all five assertions proved live by mutation |
| T3b MITC3 invariants | **blocked on a production defect** | the MITC3 file carries the five mirrored tests; four pass, `test_ke_global_leaves_all_six_rigid_body_modes_free` fails by design. See "Finding 1". |
| T4-T8 | not started | designs and measurements are recorded below |

## Finding 1 (critical): MITC3 and MITC4 disagree on the rotational DOF sign

Discovered while adding T3b. MITC3's global stiffness does not leave the physical rigid-body mode free, and the two element families interpret the shared rotational DOFs with opposite signs.

### Measured evidence

1. **Rigid-body invariance.** For the physical field `u = Omega x (x - x_c)`, `theta = Omega`: MITC3 residual `|K u| / (|K| |u|) = 7.887e-5`, MITC4 `6.5e-18`. The small number is a normalisation artifact (see 3); the penalty itself is O(1).
2. **The free mode.** MITC3 annihilates `(u = Omega x (x - x_c), theta_x,theta_y = -Omega_x,-Omega_y, theta_z = +Omega_z)` to machine precision (1e-17). Negating all three components breaks the `theta_z` case (4.3e-4). So MITC3's three rotational DOFs are not a rotation vector: the out-of-plane components are negated relative to the in-plane one.
3. **The penalty is O(1) on the relevant block.** `|K| = 3.885e9` is membrane-dominated (translational sub-block 3.885e9, rotational 5.44e5, bending/shear 7.47e5). The physical rigid rotation gives `|K u| = 5.87e5`, i.e. `|K u| / |K_bs| = 0.786`.
4. **Mixed meshes are grossly wrong.** Cantilever, 4x20 cells, transverse tip load, analytical `w = 3.047619e-01` m:

   | mesh | w_tip [m] | error | MITC3 elements |
   | --- | --- | --- | --- |
   | all quads (MITC4) | 3.026836e-01 | 0.682% | 0/80 |
   | all triangles (MITC3) | 3.016701e-01 | 1.015% | 160/160 |
   | triangles on the last row only | 2.205342e-01 | **27.637%** | 8/84 |
   | triangles on rows k>=15 | 4.676760e-02 | **84.654%** | 40/100 |
   | triangles on alternating rows | 5.690370e-04 | **99.813%** | 80/120 |

   Each family alone is accurate; any mixture is not. `MeshTopology::new` stores per-element connectivity lengths and per-element types with `dofs_per_node = max`, so mixed connectivity is supported by construction — this is not an assembler artifact.
5. **Composite B-coupling sign.** `Laminate([0/90], h=4mm)`, cantilever 2x20, axial load: MITC4Comp `uy = -2.627972e-03`, MITC3Comp `uy = +2.627691e-03`. Same magnitude to 0.01%, opposite sign. `test_b_coupling_analytical_mitc4comp` validates MITC4Comp against the Reddy CLT formula, so MITC3 is the wrong one.

### Root cause

MITC3's `eval_covariant_shear_ext` builds `e_rt = dw/dr + V3·g_r` with the director `V3 = e3 - theta_y·e1 + theta_x·e2` (its own comment). The Reissner-Mindlin director is `V3 = e3 + Omega x e3 = e3 + theta_y·e1 - theta_x·e2` for `theta = Omega`. MITC4's `b_gamma_mitc4` G-matrix implements the physical form. MITC3's curvature operator, its drilling operator (`theta_z` is the physical in-plane rotation) and MITC4 all use the physical convention; MITC3's shear does not.

Commit `b136ce5` ("fix(mitc3): correct covariant shear sign convention and add update_reference") introduced exactly this inversion, justified by the inverted director. Its parent (`b136ce5^`) had the convention that matches MITC4. The same commit also added `update_reference` (updated-Lagrangian nonlinear support), which independently accounts for 7 of the 8 tests its message claims: `test_cantilever_large_rotation_*`, `test_equilibrium_path[*]` and `test_simo_vu_quoc_rollup_360` are all nonlinear.

### Why the suite was green

- `tests/test_material_suite.py::test_axial_produces_bending_mitc3comp` asserts only `abs(uy) > 1e-8` — sign-blind, so a wrong-signed coupling passes.
- The other MITC3 composite tests use `[0]` or `[0/90/90/0]`, whose `B` matrix is zero.
- No test builds a mesh that mixes MITC3 and MITC4.
- `tests/test_mitc3_benchmarks.py::test_linear_tip_moment_sign` asserts `w_tip > 0` for a positive tip moment on `theta_y` — a statement about the convention, not the physics, which is why the inversion could be called a fix.

### Production reachability

`src/aeroelast/core/mesh/generators.py:1247` emits `ElementType.triangle` (MITC3) for degenerate NuMAD blade elements while the rest of the blade is quads (MITC4). One degenerate element makes the whole blade a mixed mesh.

### Disposition (pending user decision)

- **A** keep the strict failing test as the finding, suite reports 80 passed / 1 failed;
- **B** `#[ignore]` the test with a reason and a bug reference so the suite stays green while the defect is documented;
- **C** adapt the test to MITC3's internal convention — rejected: it hides a cross-element mismatch.

## Finding 2: the Hook benchmark expectation (1.12) is not in the paper

`tests/test_ko2017_performance.py::test_3_6_hook_table_14_minimal_fix` builds the paper's geometry exactly (R1=14, R2=46, width=20, angles 60/150, t=2.0, E=3.3e3, nu=0.3, P=1.0, an N x 6N mesh, load as a uniformly distributed traction at the tip) and normalises by `wref = 4.82482`. The paper (Ko, Lee & Bathe, *Performance of the MITC3+ and MITC4+ shell elements in widely used benchmark problems*, section 3.6) states that `wref = 4.82482` is the reference displacement obtained with the MITC9 element at N=64, and its Table 14 gives for MITC4+:

| N | 2 | 4 | 8 | 16 | 32 |
| --- | --- | --- | --- | --- | --- |
| MITC4+ | 0.9531 | 0.9635 | 0.9782 | 0.9911 | 0.9973 |

The table converges to 1.0, and the test's `expected_norm = 1.12` appears nowhere in it. Measured on the refinement sequence (this repo): N=2 0.96034, N=4 0.98175, N=6 0.98813, N=8 0.99268, N=12 0.99844, N=16 1.00164 — the same convergence to 1.0. The test's only assertion is `assert norm > 0.1`.

Also: the module-level `REFERENCE_VALUES` table (which does contain the correct 4.82482) is never read, which is how each test came to hard-code its own expectation.

## Finding 3: seven sign-blind assertions

`assert abs(x) > tiny` can never detect a sign error — the class of defect in Finding 1. Seven sites:

- `tests/test_material_suite.py:580`, `:593`, `:640`
- `tests/test_composite_b_coupling.py:202`, `:365`
- `tests/test_composite_b_coupling.py:122` (a matrix-norm check)
- `tests/test_force_projection.py:288`

## Measurements already taken for the remaining tasks

- **T5** `tests/test_beam_shell_4cases_parity.py` (slow-marked). The test named `..._with_analytical` computes the analytical value, prints it, asserts only AeroElast-vs-CCX, and calls `_ccx_bin()` first, so without CalculiX it skips entirely. Measured AeroElast vs analytical, using the loaded tip node vs the mean over the free face:

  | case | loaded node error | free-face mean error |
  | --- | --- | --- |
  | tension_axial | 5.093% | **0.560%** |
  | compression_axial | 5.093% | **0.560%** |
  | bending_fx | 0.890% | 0.890% |
  | transverse_fy | 0.681% | 0.682% |

  The axial 5.09% is a point-load artifact, not a formulation error: distributing the load over the free face instead gives 0.079%. Measuring the cross-sectional mean keeps the model identical to the CCX comparison and recovers the analytical agreement.
- **T6** `tests/test_shell_validation_fixed.py` axial case, error vs `P*L/(E*A)` by mesh: (2,1) 1.35%, (4,2) 2.70%, (8,4) 1.15%, (16,8) 0.75%, (32,16) 0.49% — monotone from (4,2) on, so the reference is right and the error is discretisation. FY: 14.11, 4.67, 1.19, 0.012, 0.39%. FZ: 9.12, 3.58, 1.92, 1.32, 1.12%. At the test's (8,4) mesh: FX 1.15%, FY 1.19%, FZ 1.92%.
- **T8** the measured values above justify a 3% window on FX/FY/FZ (was a flat 5%) and ±2% on `test_ratio_physical` (measured 0.04% inside a -50%/+20% window). `assert_relative_error` hard-fails above 5% regardless of its `tol` argument, and every call site duplicates `np.isclose(..., rtol=0.05)`.

## New work unit

- **T9** replace the seven sign-blind assertions with signed comparisons against a reference (or an explicit, documented sign assertion).

## Progress log (final state of this session)

Branch `test/physical-correctness`, 12 commits, nothing pushed.

| Task | State | Commit / evidence |
| --- | --- | --- |
| T1 mass tests | **done** | `55d65e7` — exact extraction + consistent-mass coefficients; 15 tests pass (was 9 with 1 red); two mutations prove liveness |
| T2 delete dead drill operator | **done** | `4a46af2` — Rust 73/1F -> 71/0F |
| T3a MITC4 invariants | **done** | `785a48d` — 5 invariants, all mutation-checked |
| T3b MITC3 invariants | **done, with the rigid-body test red by design** | `21c5fbe` — the disposition chosen was option A: keep the strict failing test as the evidence rather than `#[ignore]` it or adapt it to the element's convention |
| T4 shell mass invariants | **done** | `ae34315` (MITC4), `21c5fbe` (MITC3) — symmetry/PSD, total mass `rho*h*A`, exact coefficients, rotary inertia; three mutations per element prove liveness |
| T5 analytical assertion in the parity suite | **done** | `e252665` — 4 analytical cases now run instead of skipping; free-face mean instead of the loaded node (5.09% -> 0.56% axial); window 2% |
| T6 axial convergence study | **done** | `e0bdbd3` — 2.702% (4,2), 1.154% (8,4), 0.748% (16,8), 0.487% (32,16), monotone, finest below 1% |
| T7 Hook reference | **done** | `b611ad7` — the paper's Table 14 N=8 value 0.9782 replaces 1.12; 3% window; liveness proved by restoring 1.12 (fails at 11.37%) |
| T8 decorative tolerances | **done** | `e0bdbd3` (FX/FY/FZ 5% -> 3%, ratio -50%/+20% -> +-2%) and `9bb85a8` (`assert_relative_error`'s hidden 5% ceiling removed, 7 duplicate `np.isclose` calls dropped) |
| T9 sign-blind assertions | in flight | 6 sites across 3 files |

### Verified suite state after T1-T8

| Check | Before | After |
| --- | --- | --- |
| `pytest -m "not slow"` | 1 failed / 338 passed / 4 skipped | **345 passed / 4 skipped / 0 failed** (before T3b's red test was committed) |
| `pytest tests/test_shell_validation_fixed.py` | 6 passed, one window that could not fail | 7 passed, four windows mutation-proved |
| `pytest tests/test_beam_shell_4cases_parity.py` | 4 skipped without CalculiX | 4 passed (analytical) + 5 skipped (CCX) |
| `cargo test -p aeroelast-core` | 73 tests / 1 failed | 89 tests / 88 passed, the single failure being the MITC3 rigid-body finding |

Note: `cargo test -p aeroelast-mesh` and `-p aeroelast-solvers` still cannot run on this machine (the `hdf5-metno-sys` and `precice` build scripts fail), so their tests are unmeasured.

### Residual work

- **T9**: the six sign-blind sites. Two of them (`test_composite_b_coupling.py`) already state the physical sign in their own docstrings while checking only the magnitude; one of those is literally named `test_b_coupling_sign` and says "we check magnitude only" — the exact excuse that let Finding 1 through. `test_material_suite.py::test_axial_produces_bending_mitc3comp` becomes the signed CLT comparison and is expected to fail, so it is marked `xfail(strict=True)` with the defect reference, keeping the suite's red count at one while still turning into a failure the moment the production sign is fixed.
- `test_material_suite.py::test_bending_produces_extension_mitc4comp` and `test_force_projection.py:288`: signed only if the sign is determinate; otherwise the magnitude check stays with an explicit note.
- **Finding 1's production fix** remains unauthorized. The decisive experiment (revert the shear sign in `eval_covariant_shear_ext`) measured:
  - Rust suite 81/81 (the rigid-body test passes);
  - mixed-mesh errors drop to the pure-case level (0.682% and 1.009%);
  - composite and coupling suites: 48 passed;
  - but five MITC3-only tests that apply a moment to the `theta_y` DOF fail with a purely negated tip deflection (`w_tip = -6.3740` against a reference `6.3662`): `test_linear_tip_moment_sign`, `test_cantilever_large_rotation_half_circle[10]` and three `test_equilibrium_path` cases.

  So the sign is load-bearing in two directions: reverting it fixes the cross-element convention and rigid-body invariance but breaks the tests that encode the flipped convention. `test_linear_tip_moment_sign` asserts `w_tip > 0` for a positive moment on the `theta_y` DOF; under the physical convention (where `theta_y` is the rotation vector component and `theta_y = -w_,x` in the thin limit) a positive moment about +y gives `w_tip < 0`, so that expectation is convention-dependent rather than physical. The remaining four are large-rotation paths whose reported sign follows the same convention. The fix therefore has two halves and needs an explicit decision.
