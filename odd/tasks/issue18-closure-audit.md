# Audit of the closed items of #18 (2026-10-10)

Purpose: decide, issue by issue, whether each closure recorded in #18 still holds on the current
tree, as a precondition for citing the solver in scientific work. This audits **closures**, not
the solver: "closure holds" means the claim, its guard and its band survive; it does not mean
the quantity is validated against experiment (`experimental_validation` stays an open gap).

## Method

Tree measured: `integrate/origin-main-2026-09-30` at `2467dbc` **plus the uncommitted working
tree** of the shared checkout (`src/aeroelast/cli/run_bem_fsi.py`,
`tests/validation/bem/test_wall_flow_activation.py`, `tools/validation_matrix.py`, ... belong to
another session). The Rust extension (`_aeroelast`, built 2026-10-04 13:24) is newer than the
last `crates/` commit (`26ffe6e`, 13:12), so no stale `.so`.

Five checks per issue:

1. **Claim** - the closure statement, with metric, reference and tolerance.
2. **Fix on branch** - every SHA cited by #18 is an ancestor of HEAD.
3. **Guard alive** - the guard test exists and passes on the tree above.
4. **Band untouched** - no tolerance widened in the guard after the closure date.
5. **Not contradicted since** - later closure sections, moved anchors, later findings.

Verdicts: **holds** / **holds, with recorded gap** / **reopen candidate**. Reopening is the
maintainer's decision; nothing was changed on GitHub by this audit.

## Global results

- Check 2: all 16 SHAs cited by #18 (`6765633 e3278ba 7a84da1 208f220 5260233 756a3b2 81a1e15
  790b80e a9769bb f3edb65..ae25e03`) are ancestors of HEAD.
- Check 4: `git log -p --since=<closedAt>` over every guard file found **no widened band**.
  Every post-closure tolerance edit is a refactor into `assert_rel_within(tol=...)` with the same
  value (`b4f88b9`, `b6187fc`, `0a01423`, `9d284b3`), or a tightening (`ed79a92`: the MITC3
  symmetry check lost its `rtol`).
- Store gate: `tools/validation_matrix.py check` = 211 rows / 285 comparisons / 0 errors /
  0 warnings; `docs/validation` untouched afterwards (#32 not triggered).
- Guard run: 24 guard files, **188 passed, 0 failed, 0 xfail**, 529 s
  (`logs/check/issue18_audit_guards.log`). A pass proves only that each guard's own band holds;
  it does not prove the band can detect the effect it names (S-7's ±30 % band could not).

## Per issue

| # | closed | claim (short) | guard | verdict |
| --- | --- | --- | --- | --- |
| 2 | 10-01 | smoothed MITC3+ carries the relative frame; Scordelis-Lo vs Lee & Lee T6 within 5 % (0.99/1.00 at N=8/16) | `benchmarks/test_mitc3_smoothed.py` | holds |
| 3 | 10-01 | composite stiffness MITC4+CLT vs CCX: plateau explained, ABD within 1.5 % (group 18) | `parity/test_composite_layup_parity.py`, `element/test_composite_b_coupling.py` | holds (stress half -> #27) |
| 4 | 09-30 | blade modal shell-vs-CCX gap is discretisation (falls with mesh) | `parity/test_blade_iea15mw_mesh_convergence.py` | holds, with recorded gap: #8 corrected the closing claim - modes 5-8 plateau at ~2 % (element order, not mesh); shell-vs-beam is bounded only by the beam references' own scatter (5.8 % / 11.8 %) +3 %, which is a code-vs-code statement |
| 5 | 09-30 | Viterna is implemented exactly (AeroDyn manual eqs 98-102); the ~40 % post-stall gap is NeuralFoil input | `bem/test_bem_polars.py` | holds, with recorded gap: the post-stall gap is attributed, not removed. Production with the official YAML parses the deck's polars (`aerodynamics.py:726`); NeuralFoil+Viterna is only a fallback / xlsx path. Any run on a NeuralFoil polar inherits the gap |
| 6 | 09-30 | nonlinear static converges on defaults (`max_it` 100); cantilever 0.172820 m | `benchmarks/test_large_rotation_benchmarks.py` | holds, with envelope: P >= ~2 is outside the validated envelope (solver docstring) |
| 7 | 10-01 | MITC4 `K_G` applies the thickness factor | `element/test_stress_stiffened_solver.py` | holds (S-4 rotating modes vs OpenFAST are the downstream anchor: +0.6 / -4.5 / -1.5 %) |
| 8 | 10-01 | per-mode convergence claim restated per group (LOW 1 %, HIGH 5 % element-order) | same as #4 | holds |
| 9 | 10-04 | MITC4 exonerated vs three references (1-3 %); coupon reference was wrong; D16 = D26 = 0 on the blade | `blade/test_laminate_bend_twist*.py`, `blade/test_blade_deloading_vs_reference.py` | holds; the physics gap moved to #11/#12 |
| 11 | 10-06 | projector realised a section moment as circle-tangential field | `parity/test_thin_walled_tube_projection.py` | holds, with recorded gap `moment_realization_over_delivers` (production still runs the minimum-norm field; #16 open) |
| 12 | 10-07 | twist / de-loading 2.4x / 3-4x short of reference: split into items | (split) | holds as a split, **not as a physics closure**: the coupled twist against literature was never closed - #14 declared it non-transferable. Its numbers are pre-#26/#19 and uncitable |
| 13 | 10-07 | the path-measure artefact contributes zero to the radii feedback's +1.24 % thrust | `blade/test_blade_deloading_vs_reference.py` (span-projection identity, 1e-9) | holds as a measurement-identity statement; the +1.24 % itself is a pre-#26/#19 coupled number and uncitable |
| 14 | 10-08 | Zhou 2025 twist not transferable; CCX S8R arbiter with the same nodal vector | `parity/test_thin_walled_tube_projection.py` + `tools/ccx_blade_twist_arbitration.py` | holds as structural parity (both codes get the same vector, so a load-frame error cancels); absolute twist magnitudes are pre-#26 (`7a84da1` landed 2 h after the closure) |
| 15 | 10-08 | the +31 % tangential is the reference's frame plus our inflow bias, not a polar defect | `bem/test_bem_load_frame.py` | holds: the claim is on BEM `Np/Tp` before projection, which #26 did not touch. Residual gap `zhou_spanwise_load_frame` |
| 19 | 10-09 (10-08 in #18) | projection keeps reference geometry; coupling contracts, 500/500 windows, mean 2.60 sub-iterations | `rotor/test_bem_fsi_deformed_geometry.py` | holds (convergence claim, not accuracy); residual gap `reference_projection_geometry_approximation` |
| 20 | 10-08 | S-7: shell is 6.4 % **stiffer** than BeamDyn; the "softer" table was the artefact | `blade/test_iea15mw_s7_torsion.py` | **reopen candidate.** `72f64de` (2026-10-10) found the tip couple applied `kappa*M`, `kappa = 0.8978`; corrected ratio converges to 1.039 => shell ~3.8 % **softer**. The closure's direction is inverted, and S-1 sectional (+23.9 % stiffer) vs S-7 global (softer) now disagree in sign with nothing attributing it. GitHub still shows #20 CLOSED, and `AGENTS.md` still carries the superseded S-7 anchor `0.9396 -> 0.9336` |
| 21 | 10-08 | `regression` honours `non_validation_tests` | tooling | holds (gate integrity) |
| 26 | 10-09 | BEM `Np/Tp` ride the rotor plane, not the section | `bem/test_force_projection_load_frame.py`, `bem/test_force_projection_ac_datum.py` | holds (twisted known-answer case) |
| 27 | 10-09 | composite ply stress vs CalculiX S8R COMPOSITE and a closed form (group 36) | `parity/test_composite_ply_stress_parity.py`, `parity/test_composite_stress_ccx_parity.py` | holds, with recorded gap `composite_stress_recovery` (two limits) |
| 28 | 10-09 | recovery returned the thickness mean, now resolves the ply | same as #27 | holds |
| 29 | 10-10 | `regression` runs multi-file groups; one comparison per printed residual | tooling | holds (`check` green above); its limits are #31/#32 |
| 30 | 10-10 | wall-flow coupled movement is structural response to the pattern | `parity/test_thin_walled_tube_moment_realization.py` | holds as attribution; the activation decision stays in #16 |

## Findings on #18's own bookkeeping

1. **#20 is an orphan in the sense of #18's rule**: the artifact (closure section, corrected
   test) says the closure was inverted, the issue is closed. #18 must list it under *State to
   settle* and the maintainer must reopen or open a successor (the S-1 vs S-7 sign question).
2. **#30 is missing from the closed index** of #18, and its closing paragraph sits under the
   #27 entry in *Entry and closure per item*.
3. **#19 / #26 closure dates**: GitHub says 2026-10-09, #18's index says 2026-10-08.
4. `AGENTS.md` S-7 row (`superseded 2026-10-07: 0.9396 -> 0.9336`) is stale after `72f64de`;
   it should point at 1.0466 (fixture) / 1.039 (finest). Not edited mid-session per the
   AGENTS.md rule.
5. `AGENTS.md` says `test_rotor_inertial.py` was removed, but
   `tests/validation/rotor/test_rotor_inertial.py` exists.

## What this means for scientific use

- Element-level and code-vs-code claims (#2, #3, #4/#8, #7, #9, #26, #27/#28) hold, with no
  band widened. They are verification, not validation.
- Coupled-path numbers older than `7a84da1` (2026-10-08 16:01) and `6765633` (2026-10-09) are
  uncitable; #12, #13 and #14's absolute magnitudes are in that set.
- **Torsional stiffness of the blade is not settled**: after `72f64de` the global shell is ~3.8 %
  softer than BeamDyn while the sectional S-1 says +23.9 % stiffer. Any twist result should be
  quoted with that open question until #20 (or its successor) closes.
- Nothing in the list is validated against experiment.
