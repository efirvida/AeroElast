# Feature: composite outer-fibre stress recovery is not delivered yet (#27, item 4 of #18)

Status: T1-T6 done 2026-10-09 (`95bfbab`, `208f220`, `06f8be0`, `bd4cc05`, `62d806a`,
`5260233`); T7 (the isotropic guard, already green) is the only item left, and the work is
published (push + #27 comment/close + #28 + #29)
Owner: this session (2026-10-09)
Related: issue #27 (roadmap item 4, `P1`, successor of #3), issue #18 (roadmap and
order of record), `docs/validation/gaps.yaml` id `composite_stress_recovery`,
`docs/validation/rows/18-composite_layup.yaml` and `groups.yaml` group 18 (the
stiffness half, closed), `odd/tasks/composite-bend-twist-verdict.md` (the #3 record).
New issue required by the roadmap rule: the smeared recovery found here is an unknown
found mid-item, so it gets its own issue and its own row in #18 instead of being
absorbed into #27.

> **Handoff 2026-10-09 - read this block first.**
>
> **Done and committed** on `integrate/origin-main-2026-09-30` (local, nothing pushed):
> `95bfbab` T1 RED (6 failed, 1 passed: the smeared 2.0000 MPa against the closed form
> 3.6959 / 0.3041, i.e. -45.9 % / +557.7 %); `208f220` T2 + T2b GREEN (14 passed;
> isotropic bit-identical, 0 of 108 captured arrays differ; the Rust `_aeroelast.Laminate`
> production form resolves the same stack as `CompositeShellProperty` with max |diff|
> 0.000e+00 Pa; a raw ABD dict keeps the thickness mean and warns once per element set);
> `06f8be0` this note; `bd4cc05` T3 coverage (20 passed in 2.95 s: pure bending
> `kappa = D^-1 M` within 0.0026 % / 0.72 % and station-independent to 9.9e-06, the +/-45
> discrimination lives in `sigma_xy` because a balanced laminate has `A16 = A26 = 0`, and
> the `MIDDLE` tie-break is pinned by ply index).
>
> **Done in the resumed session: T4, the CalculiX row test, `62d806a`** (the first writer
> attempt failed on 2026-10-09 and created no file; the relaunch with the same specification
> produced 578 lines, 7 passed in 2.9 s). `TOL_CCX = 0.05` as specified, two inline
> references (CalculiX `kind="code"`, closed form `kind="analytical"`), each code on its own
> 8x2/16x4 sequence, read at the free-field centre.
>
> **Its anchors moved and that is the finding.** At 16x4 ours is 3.6959 / 0.3041 MPa (the
> closed form the same to <= 0.0006 %) while CalculiX gives 3.6853 / 0.3148 MPa, i.e. 0.29 %
> and 3.41 % - not the 0.07 % / 1.6 % the E2 probe had quoted. The cause is a deck/route
> difference that could not be recovered (E2 correction below), and the `[90/0]s` gap does
> not close with refinement, so `TOL_CCX` is recorded as a bound and not as a convergence
> claim. Two further corrections came out of the verification (E4 attribution, expanded FRD).
>
> **Next step: T7 only (the isotropic guard, green and unchanged), then this note closes.**
>
> **Published 2026-10-09 with the user's go-ahead:** push `adce7ab..a198ce0` to
> `origin/integrate/origin-main-2026-09-30`; issue **#28** filed for the smeared-recovery defect and
> closed the same day with the fix pointer; issue **#29** filed for the two `validation_matrix` gate
> defects found while running T5; a closing comment on **#27** (comment 6088796371) plus close; the
> **#18** body updated (item 4 out of the open table, `Closed 2026-10-09` block, `#28` in the closed
> index, `#29` at P4) and the closure record written as the `#27` section of
> `docs/validation_closures.md`.
>
> **Environment:** every command needs
> `bash -lc 'module load glu gcc/14.2.0_sequana; export LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64:$LD_LIBRARY_PATH; ...'`
> or the import dies on `CXXABI_1.3.15`.
>
> **The working tree is not ours right now:** a staged `CLAUDE.md -> AGENTS.md` rename plus
> modified `docs/blade_input_divergence_utd_vs_official.md`,
> `tests/test_iea15mw_v05_structural_properties.py`, `odd/tasks/token-efficiency.md` and
> `scripts/---` appeared during this session and belong to another session. Never
> `git add -A` here.
>
> **Review workload:** `assess` flagged `reviewDue: slice_budget_reached` - `208f220` is 611
> insertions against a ~400-line soft budget. T5 is the place to split.

## Why this exists

Issue #27 was written on the assumption that the missing thing is the **reference**:
"`*SHELL SECTION, COMPOSITE` ignores `OUTPUT=3D`, so CalculiX cannot act as the
independent stress judge it is for isotropic shells". Measured on 2026-10-09, both
halves of that assumption are wrong, in opposite directions:

* CalculiX **can** judge the outer fibre (E2), it just reports one value per node on
  the midsurface, and `OUTPUT=3D` is inert for composites (E3);
* the recovery in the tree **does not produce a ply stress at all** (E1). It produces
  the homogenised, thickness-averaged equivalent stress of the section.

So the item cannot close on a row today: there is nothing on our side to compare
against the reference. The user's decision (2026-10-09) is to build the capability
test-first and then close #27 with the CalculiX + closed-form row, instead of closing
item 4 on a weaker claim.

## Evidence as of 2026-10-09

Reproduce with `module load glu gcc/14.2.0_sequana` plus
`LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64` and the probes in
`/tmp/ccx_composite_probe/` (`probe.py`, `probe2.py`, `probe3.py`); see the instrument
notes at the bottom before re-deriving any of this.

### E1 - the delivered recovery is smeared, not ply-level

Membrane tension (`L=1.0`, `B=0.1`, four plies of `0.00125 m`, `E1=120 GPa`,
`E2=10 GPa`, `G12=5 GPa`, `G23=3 GPa`, `nu12=0.3`, `F=1000 N`), centre element, free
field:

| quantity | `[0/90]s` | `[90/0]s` |
| --- | --- | --- |
| CLT closed form, outer ply | **3.696 MPa** | **0.304 MPa** |
| CLT closed form, second ply | 0.304 MPa | 3.696 MPa |
| CLT thickness mean `A/h . eps0` | 2.000 MPa | 2.000 MPa |
| `StressRecovery`, TOP / MIDDLE / BOTTOM | **2.0000 / 2.0000 / 2.0000** | **2.0000 / 2.0000 / 2.0000** |

Stable across `8x2` and `16x4`, and independent of the stack: no `z` dependence in
membrane, no dependence on the ply whose fibre is outermost. Error against the real
outer fibre: **-45.9 %** for `[0/90]s`, **+558 %** for `[90/0]s`. In bending the error is
smaller by accident: our 616.10 MPa against CalculiX 595.59 MPa (-3.33 %) on
`[0/90/90/0]`.

Cause, not a scale slip: `src/aeroelast/postprocess/stress_recovery.py:24-32` defines
`sigma_m = (C / h) . (B_m . u_e)` with `C = Cm() / h`, and for the composite element
`Cm()` returns the laminate's integrated `A`, so `C / h` is the **homogenised section
modulus**. The documented `h` trap is the mechanism of the smear.

### E2 - CalculiX's composite-shell FRD value is the outer-fibre ply stress

Same deck through `write_ccx_mesh(..., quadratic=True)` with
`*SHELL SECTION, COMPOSITE` (S8R), membrane tension, **free-field** node at the
mid-length:

| laminate | CCX 8x2 | CCX 16x4 | CLT outer ply |
| --- | --- | --- | --- |
| `[0/90]s` | 3.6984 MPa | 3.7005 MPa | 3.696 MPa (+0.07 %) |
| `[90/0]s` | 0.2991 MPa | 0.2994 MPa | 0.304 MPa (-1.6 %) |

The value follows the closure ply, so it is the outer fibre and not a thickness mean
(which would be 2.000 MPa for both stacks).

**E2 mechanism confirmed by the T4 row, exact splits not, and the reason is now known.**
The composite card **expands the section in the FRD**: the file carries a midsurface node
(exactly 0.00000 MPa) plus nodes at ply boundaries and ply midpoints - nine z-levels in
all, ±0.0025, ±0.001875, ±0.00125, ±0.000625 and 0 - with interface nodes duplicated per
ply. The T4 row reads the topmost surface node (`z = +H/2`), which is the outer fibre and
carries the closure ply's stress: 3.6820 / 3.6853 MPa for `[0/90]s` and 0.3150 / 0.3148 MPa
for `[90/0]s` at 8x2 / 16x4. No node choice on the T4 deck can produce the probe's
3.6984 / 0.2991 MPa (the centre field holds only 0.0000, 0.3150 and 3.6820 MPa), so the
probe ran a different deck; what differed (load weighting, mesh aspect, section offset) is
unrecoverable because `/tmp/ccx_composite_probe/` no longer exists. Recorded rather than
smoothed over: the table above is the probe's measurement, the row's own numbers are
0.2-1.8 percentage points away from it, and the row states both.

### E3 - `OUTPUT=3D` is inert for composite sections

Bending-dominated `[0/90/90/0]` cantilever, same deck, `OUTPUT=2D` against `OUTPUT=3D`:
max `|SZZ|` 595.588 MPa both, max von Mises 574.521 MPa both, **relative difference
0.000e+00**. This confirms and sharpens the note already committed in
`tests/validation/parity/test_shell_stress_ccx_parity.py` (899.67 MPa unchanged);
that note was measured on a membrane-dominated case, this one on pure bending.
Re-confirmed by the 2026-10-09 T4 verifier on the T4 deck: `OUTPUT=3D` reproduces the centre
(3.68528 MPa) and the maximum (6.5803 MPa) value-for-value on `[0/90]s`; only the 233
zero-stress midsurface reference nodes disappear from the file (2437 -> 2204 stress nodes).

### E4 - the loaded free edge is a load-introduction singularity, not the answer

On the membrane case the maximum interior `SXX` **grows with refinement** (4.4413 MPa at
8x2, 5.0905 MPa at 16x4 on the original probe; 5.8838 -> 6.5803 MPa for `[0/90]s` on the T4
deck, whose maximum sits on the node at `(1.0, 0.0, -H/2)`) while the free field stays at
3.700 MPa. **Correction to this section as first written:** the singularity is at the
**loaded free edge** (`x = L`), not at the clamped restraint - measured by the T4 verifier on
2026-10-09, the clamped-edge value is not singular (3.69 -> 3.71 MPa). The conclusion is
unchanged, and it is what the row asserts: read the centre of the patch, never the global
maximum, or the row validates a mesh-dependent effect.

## Contract note for #27

The issue's `closes-when` ("a converged-vs-converged comparison against CalculiX S8R
with `*SHELL SECTION, COMPOSITE` ... plus a row with its reference and its bound") stays
as written; only its premise changes: the CalculiX route exists and the missing half is
ours. The premise in the issue body and in the `composite_stress_recovery` gap reason
must be corrected in the same push as the closure.

## Tasks

* **T1 - RED: ply-level outer-fibre stress against an independent closed form.**
  New `tests/validation/parity/test_composite_ply_stress_parity.py`. The reference must
  be re-implemented in the test (the `_hand_clt_abd` pattern of
  `test_composite_layup_parity.py:375-438`: own `Q`, own `Qbar(theta)`, own `A`), not
  imported from `aeroelast.core.laminate`. Cases: `[0/90]s` and `[90/0]s` membrane,
  asserted against 3.696 MPa and 0.304 MPa. RED today because the smeared value misses
  the bound by 46 % / 558 %. Solve through the **production** path (the standing rule:
  measure through the production path, never through a test-only construction), i.e.
  the Rust assembly plus the production static solver as
  `test_shell_stress_ccx_parity.py::test_production_solver_matches_the_scipy_replica`
  does, feeding `StressRecovery`.
  **Done 2026-10-09: `95bfbab`.** RED observed, verbatim: `2e+06 vs 3.69592e+06 (error:
  45.8863%)` for `[0/90]s` and `2e+06 vs 304080 (error: 557.7223%)` for `[90/0]s`, both
  fibres, plus the shape guard (`TOP` equals the thickness mean to 0.0000 %). 6 failed,
  1 passed in 2.23 s. Case: 8x2 QUAD4, 27 nodes, 96 DOF, mean tip `u_x = 3.111571e-05 m`.
* **T2 - GREEN: ply-resolved location, not a parallel API.**
  **Design decision (2026-10-09): fix the existing API instead of adding a second one.**
  `compute_element_stresses(location=...)` and `compute_nodal_stresses(location=...)`
  promise the stress at a through-thickness location; for a composite the stress at
  `TOP` **is** the outer ply's stress in its own axes, so the smeared answer is a bug in
  that promise, not a missing feature. Resolve the ply that contains the requested `z`
  and evaluate `sigma = Qbar(theta_ply) . eps(z)` with `eps(z)` from the existing strain
  recovery (check the Voigt shear convention against the Rust `compute_stress_field`
  before trusting it). The ply stack comes from the element's
  `CompositeShellProperty.laminate`, because the Rust element only holds ABD. For an
  isotropic section the new path must reproduce the old numbers to machine precision:
  that is the no-regression guard (T7).
  **Done 2026-10-09: `208f220`** (with T2b, below). 14 passed; centre-element sigma_xx
  `[0/90]s` 3.6960 / 0.3040 / 3.6960 MPa and `[90/0]s` 0.3040 / 3.6960 / 0.3040 MPa,
  errors 0.0017 % / 0.0210 % against the closed form. Isotropic bit-identity: 0 of 108
  captured arrays differ for the 52.46 MPa plate. Parity suite 68 passed, 1 xfailed.
  Known limit, now loud instead of silent: a raw homogenised ABD dict
  (`{"type": "composite", ... cm/cb/cs}`) carries no ply stack, so those elements keep the
  thickness-mean value and emit one `logger.warning` per element set.
  Review workload: 611 insertions in one slice, over the ~400-line soft budget that
  `assess` flagged as `reviewDue: slice_budget_reached` on 2026-10-09. Split in T5.
* **T2b - the production section form (done with T2, `208f220`).** T2 alone only covered
  `CompositeShellProperty`, which production does not use: `runner._extract_blade_properties`
  (`src/aeroelast/solvers/fsi/runner.py:1072-1090`) hands the solver the map built by
  `build_rust_properties` (`src/aeroelast/models/blade/model.py:721`), whose values are
  Rust `_aeroelast.Laminate`. That object does carry the stack (`plies` getter,
  `crates/aeroelast-py/src/materials.rs:288`: dicts with `thickness`, `angle` in degrees,
  `z_bottom`, `z_top`, `material`). The three section forms now share one normalisation, the
  Rust form and `CompositeShellProperty` agree bitwise (max |diff| 0.000e+00 Pa over the 16
  elements), and the production-shaped case is built in the test through
  `build_rust_properties` itself. Trap found here: the ply span is
  `plies[-1].z_top - plies[0].z_bottom`, not `z_top` (plies are centred on `z = 0`); the
  wrong span halves `h` and sends `BOTTOM` to the wrong ply.
* **T3 - coverage beyond the discriminator.** Pure bending of a symmetric laminate
  (`sigma(z)` linear within each ply, TOP and BOTTOM of opposite sign) against
  `kappa = D^-1 M`; an angle-ply case (`[+-45]s`) where the flight to the ply axes in
  `Qbar(theta)` is what is under test; and the `MIDDLE` tie-break convention pinned
  explicitly (at a ply interface the half-open rule picks the first ply whose
  `z_bottom <= z`; either neighbour is defensible, so the rule must be asserted, not
  assumed). **Out of scope on purpose**: an additive all-plies API
  (`compute_ply_stresses`). Production reads TOP/MIDDLE/BOTTOM through
  `compute_nodal_stresses_all_layers_dict` (`src/aeroelast/solvers/fsi/rotor.py:2535`), so a
  ply-wise field is a new capability with no consumer yet; it gets its own issue if it is
  ever wanted, and it would double this item's review workload for no claim.
* **T3 done 2026-10-09: `bd4cc05`.** 20 passed in 2.95 s in
  `tests/validation/parity/test_composite_ply_stress_parity.py`. Two premises of this plan
  did not survive measurement and the test records why: the outer plies of a symmetric
  `[0/90]s` stack share the same angle, so `BOTTOM` is exactly `-TOP` at the outer fibres
  (the non-mirror pair is the internal `z = +/-h/4` pair, +0.1356 MPa at 0 deg against
  -0.0109 MPa at 90 deg), and under uniaxial load the `+/-45` plies share the same global
  `sigma_xx`, so the angle-ply discriminator is `sigma_xy` (+8.1619e+05 against -8.0127e+05
  against a thickness mean of 0), on a 16x4 mesh where the 8x2 clamped-edge boundary layer
  had faked a 2.4 % `sigma_xx` split that refinement removes.
* **T4 - CalculiX row. Done 2026-10-09: `62d806a`** (the writer attempt of the earlier
  session failed and left no file; the relaunch of the same specification produced
  `tests/validation/parity/test_composite_stress_ccx_parity.py`, 578 lines). Converged-vs-
  converged outer-fibre comparison against S8R with `*SHELL SECTION, COMPOSITE` on each
  code's own 8x2/16x4 sequence, read in the free field (E4), two inline references
  (CalculiX `kind="code"`, closed form `kind="analytical"`), `TOL_CCX = 0.05`. 7 passed in
  2.9 s; the isotropic sibling unchanged (3 passed, 1 xfailed); ruff clean. The row's notes
  carry the measured E2/E3/E4 characterisation, including the two corrections found by
  verification (loaded-edge singularity; an expanded FRD whose midsurface node carries
  zero). Review workload: 578 insertions, the same order as `208f220` (611), so T5 carries
  the store row and T6 the prose and issue payload instead of folding more into this file.
* **T5 - store. Done 2026-10-09.** New group `36` (`composite_ply_stress`, `groups.yaml`)
  covering both test files, chosen by the user over "one group per file" and over extending
  group 18; the store before this change did not cover either file (`status`: two validation
  files neither grouped nor declared). `rows/36-composite_ply_stress.yaml` holds 13 rows /
  18 comparisons, 12 measured; `non_validation_tests` declares the 8 own-field tests with a
  reason each; no group-level `reference_kind` (each comparison declares its own, policy rule
  7). `gaps.yaml` `composite_stress_recovery` **sharpened, not removed**: the false premise
  ("CalculiX cannot judge composites") is gone and only two genuinely uncitable residuals
  remain (a raw homogenised ABD dict keeps the thickness mean and only warns; the
  `[0/90/90/0]` bending measurement 616.10 vs CalculiX 595.59 MPa has no row).
  `residual-patterns.json` gained the group-36 pattern (without it `regression` cannot read
  the group's prints). Gates: `check` `211 row(s), 277 comparison(s), 0 error(s), 0 warning(s)`;
  `status` `neither grouped nor declared: 0`; both modules unchanged at 27 passed.
  Do not run the `coherence` write path: it rewrites the whole row file and wipes the
  hand-written prose and `measured` blocks (issue #24, still open).
  **Two tool limits found here, both real and both outside this unit (in the T6 payload):**
  (a) `regression` cannot run a multi-file group in one call - `capture_prints`
  (`tools/validation_matrix.py:3589-3592`) passes the space-joined `source_files` as one
  pytest argv element while `triage` (`:3436`) splats it, so `regression --group 36` fails
  with "no node output captured"; the per-file `--scope` form is the working equivalent and
  is what T5 used; (b) `regression --scope <file>` without `--group` silently falls back to
  the parser default `--group 3` (`:4076`) and prints misleading `source digest: moved`
  lines - always pass `--group` explicitly. A third limit, inside the extractor: a canonical
  assertion inside a loop body is deduped to one comparison per site while the run prints
  one residual per iteration, so 3 rows / 6 comparisons stay `not_measured` with a `notes`
  reason instead of a guessed baseline.
  **The `regression` exit code is not green for this group (3 `unmapped`) and is also not
  green for the committed group 7 (1 `unclaimed`, 1 `unmapped`)**: a non-zero `regression`
  exit is pre-existing store practice, not introduced here, and "regression passes" cannot
  be claimed for either group today.
* **T6 - closure and publication. Done 2026-10-09.** Issue **#28** for the smeared-recovery defect
  (E1 table, mechanism, impact/scope, reproduction), filed and closed the same day as fixed; issue
  **#29** for the two `validation_matrix` gate defects plus the loop-dedup limitation (separate
  issue: tooling, not physics, so a different root class); closing comment on **#27** with the CCX
  evidence, the satisfied `closes-when`, the moved anchor and the bound caveat, then close; **#18**
  body updated per its own maintenance rule (open-table row removed, closed-index rows for #27 and
  #28, `#29` at P4 with its Entry/closes-when, the order section and "State to settle" refreshed);
  closure record `docs/validation_closures.md` § `#27` (Spanish, like its neighbours); push
  `adce7ab..a198ce0`.
* **T7 - guards. Verified green during T4-T6 and not a remaining task:** the isotropic outer-fibre
  test (`test_shell_stress_ccx_parity.py`) stays green and unchanged in its numbers - measured
  `3 passed, 1 xfailed` on every run of this work unit.

## Instrument notes (they cost time, do not re-derive)

* `import aeroelast` fails on SDumont without
  `module load glu gcc/14.2.0_sequana` and its `libstdc++` on `LD_LIBRARY_PATH`; a
  verifier subagent burned 30 min on a hanging bash call for exactly this.
* The workable composite CalculiX deck recipe is
  `tests/validation/parity/test_orthotropic_shell_parity.py` (`_write_ccx_inp` +
  `write_ccx_mesh(..., quadratic=True)`); S8R is required, S4 does not accept
  `COMPOSITE`.
* `tests/support/ccx_io.py` gives `parse_frd_stress`, `von_mises_from_voigt`;
  `ccx_bin_or_skip` lives in `tests/conftest.py` (`CCX_BIN` env var, then `PATH`).
* The FRD coordinates of an **isotropic** shell deck are midsurface nodes. That is not true
  for a composite section: the card expands the section, so the file carries nine z-levels
  (a midsurface node at exactly zero stress plus ply boundaries and ply midpoints, interface
  nodes duplicated per ply) and the through-thickness location of the reported value *is*
  readable off the file. The first version of this note said otherwise; the T4 row's ability
  to name the outer-fibre node depends on the correction.
* The validation store's `regression` gate takes `--group` explicitly, always: without it the
  parser default is `--group 3` (`tools/validation_matrix.py:4076`) and the output looks like
  a source-digest move. And a group whose `source_files` has more than one entry cannot be
  run in one call at all (`capture_prints` joins the paths with a space, `:3589-3592`); run it
  once per file with `--scope`.
* Unrelated discrepancy noticed on 2026-10-09 while running the isotropic sibling: its
  docstring records "CalculiX OUTPUT=3D 57.41 MPa" while the local ccx 2.20 reports
  `51.64` (the live assertion is the 5 % comparison, so the test passes at 52.46 against
  51.64). Likely a CalculiX-version or route difference, not a regression from this
  change; the stale prose figure belongs with the #23/#25 class of store and prose debt.
