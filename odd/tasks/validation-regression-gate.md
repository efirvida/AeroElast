# Feature: the store's `regression` gate runs a multi-file group, and a comparison is one printed residual

Issue: **#29** — `[validation] regression cannot run a multi-file group, and --group defaults silently to group 3`
Roadmap: item **P4** of **#18**. Feature branch: `integrate/origin-main-2026-09-30` (the integration branch;
this repository does not work on `main`). Related open items: **#24** (a `coherence` run rewrites the whole
row file and wipes hand-written prose and recorded measurements), **#23**/**#25** (a stale declared citation
site and pre-fix prose — the three pre-existing red tests of `tools/tests`).

## Objective

Three defects in `tools/validation_matrix.py`, measured in #29, plus the prerequisite that turned out to sit
under them:

1. **`regression` cannot run a multi-file group.** `command_regression` joins `source_files` with spaces and
   `capture_prints` (`:3585`) hands that single string to pytest as **one** argv element, so pytest tries to
   collect a file literally named `a.py b.py`. `command_triage` already splats (`*scope.split()`, `:3436`); the
   fix is to mirror it. Only group **36** declares more than one `source_files` entry today
   (`test_composite_ply_stress_parity.py`, `test_composite_stress_ccx_parity.py`), so this is its blast radius
   for now and every future multi-file group's.
2. **`--group` defaults silently to 3.** `regression --scope <file>` runs with the parser default `--group 3`
   (`:4076`) and prints `source digest: moved` against an unrelated group. `extract` carries the same default
   (`:4063`).
3. **A canonical assertion inside a `for` body is deduped.** The extractor emits one comparison per distinct
   call site (`tolerance_sites` dedups by `(line, kind, value)`, `:2479-2482`), while the run prints one
   residual per iteration. `compare_row` (`:3750`) pairs the Nth print with the Nth asserted comparison and
   refuses a count mismatch with `unmapped`, because attaching a margin to the wrong comparison would
   manufacture a baseline. Measured in group 36: 4 printed residuals against 2, 1 and 3 asserted comparisons
   → 3 rows / 6 comparisons stuck at `not_measured`.
4. **Prerequisite found while proving the fix (scope decision A of this session):** `command_extract` (`:2868`)
   and `command_coherence` (`:3027`) read **only `source_files[0]`**. A multi-file group therefore cannot be
   re-derived, and `extract --write` **deletes** the rows of every other declared file — that is how group 36
   lost its two CalculiX-parity rows (`test_outer_fibre_stress_matches_ccx_and.sym_*`) and why
   `coherence --group 36` reports `the extraction itself did not finish cleanly` (its
   `non_validation_tests` include tests that live in the second file).

## Decisions taken (this session, with the maintainer)

- **Row model = one comparison per printed residual (option A).** A comparison stays scalar (`measured` is one
  residual), so `margin_pct`, slack, `list --unit comparison` and the evidence digest keep their meaning, and
  policy rule 2 — *the comparison is the unit: one result against one independent reference* — is honoured
  literally. The rejected alternative was a row-declared print→comparison mapping, which would make
  `measured` multi-valued and force a redefinition of slack and the digest.
- **`--group` is required** on `regression` and `extract`; no default is assumed. When `--scope` is given and
  no group is, the error names the groups that declare that scope, so the fix is one token away.
- **The multi-file derivation of `extract`/`coherence` is in scope** for this work unit: it is the same defect
  ("a group is not really a set of files") and it is a hard prerequisite for group 36's loop rows to reach the
  store without deleting the CCX rows.
- **The measurement-key half of #24 stays out of scope** and is *not* fixed here: `preserve_measurements`
  keys by `tolerance.source` = `"<file>:<line>"`, so a moved line loses its `measured` block. That is why
  `test_every_group_re_derives_to_the_rows_on_disk` will still be red after this work unit (groups 6, 20, 31,
  34, 35 are stale for exactly that reason). Recorded, not absorbed.

## Row model A — how it is implemented

A canonical call site may execute more than once inside one test node. Only a `for` loop over a **literal
sequence** (a tuple/list of literals, or a module constant holding one — the same `resolve_literal` the
tolerance and parametrisation readers already use) gives a statically knowable multiplicity. Every measured
case in group 36 is of that kind:

| row | source shape | sites × iterations | prints |
| --- | --- | --- | --- |
| `test_bending_top_and_bottom_match_clt_at_two_stations` | `for station in ("centre","three_quarter")` + 2 asserts | 2 × 2 | 4 |
| `test_bending_kappa_matches_D_inverse_M` | nested `for station in (...)`, `for comp in (0,1)` + 1 assert | 1 × 4 | 4 |
| `test_angle_ply_rotates_Qbar_per_ply` | 2 asserts + `for label, value, ref in ((..),(..))` + 1 assert | 3 + 1 × 2 | 4 |

Rules:

1. The extractor walks the test function's statement tree **in execution order**, binding loop targets to the
   literal elements, and emits one `ToleranceSite` per execution.
2. **Ordering is the contract.** `compare_row` pairs the Nth print with the Nth comparison, and prints come out
   in execution order. A per-site `for site in sites: for iteration in ...` expansion would NOT reproduce it:
   for two sites in one loop body the run prints `[site1/s0, site2/s0, site1/s1, site2/s1]`, not
   `[site1/s0, site1/s1, site2/s0, site2/s1]`. The walk must therefore re-enter the loop body once per
   iteration, so the emitted order *is* the runtime order.
3. A site whose multiplicity is **not** statically knowable (a loop over a computed sequence, a `while`, an
   `if` whose test is not a literal, or a site inside a declared helper) is emitted **once** and reported as
   `dynamic_multiplicity` in the `extract` payload. Its row still ends `unmapped` at `regression` time —
   refusing to guess is unchanged; what changes is that the reason is now named by the tool instead of by a
   hand-written `notes` line.
4. `tolerance.source` stays `"<file>:<line>"` — three readers parse it that way (`preserve_measurements`,
   the `triage` review sheet at `:3505`, the row-schema validator). The iteration is distinguished in
   `label` (`rtol at line 873 (station=centre)`), and `preserve_measurements` consumes several comparisons
   that share a source **in order** instead of collapsing them into a dict.
5. A single-execution site emits byte-identical rows to today: no `label` suffix, one comparison, same
   `source`. Every group except 36 must re-derive without a diff.

## Tasks

- [x] **T1 — `regression` runs the group's whole scope in one pytest call.** `capture_prints` splats
      the scope, mirroring `command_triage`. Commit `f3edb65`;
      `test_capture_prints_runs_every_scope_argument` stubs `subprocess.run` and pins the argv.
- [x] **T2 — `--group` is never silently assumed.** `require_group` replaced the `default="3"` in
      `regression` and `extract`; the error names the groups that declare the `--scope`. Commit
      `99bd123`. The one test that relied on the default now states the group.
- [x] **T3 — a group's scope is every declared source file.** `command_extract` derives each
      `source_files` entry and merges them into one report (`used_ids` keeps row ids unique; the stale
      classification runs over the union); `command_coherence` passes the whole set. Commit `f199860`.
      `coherence --group 36` is now `coherent`.
- [x] **T4 — row model A.** `execution_contexts` walks the body in execution order and emits one site
      per iteration of a `for` over a tuple/list literal (a literal states its length even when its
      elements are computed); `binding_label` names the iteration; `dynamic_multiplicity` reports a
      multiplicity the code does not state; `preserve_measurements` consumes duplicate sources in
      order. Commit `d49400c`. Verified store-wide on a copy: only groups 26 and 36 change, no other
      group's derivation moves a byte.
- [x] **T5 — the tool's own tests.** 11 new tests, all on fixture stores or the real store read-only;
      `-m "not slow"` gives 133 passed and exactly the three pre-existing #23/#25 reference failures.
- [x] **T6 — group 36 (and 26) re-derived and measured.** `extract --group 36 --write` → 13 rows / 24
      comparisons with all 12 prior measurements carried; `regression --group 36 --write` → 12 new
      baselines; a second run → `24 same`, exit 0. `check`: 211 rows / 285 comparisons / 0 errors / 0
      warnings. Commits `fc24808` (36) and `9cd0c02` (26, the same loop shape: 3 coupons, 0.90% /
      0.51% / 0.24%). Group 26's gate stays red on 6 pre-existing `unclaimed` tests.
- [x] **T7 — the docs that state the old model.** `odd/tasks/validation-matrix-store.md` §8,
      `docs/validation-policy.md` rule 2, and the `#29` section in `docs/validation_closures.md`.
- [x] **T8 — published (2026-10-10, authorized).** Push `a9769bb..ae25e03`. Comment on #29
      (`6098795177`) with the measured outcome, and #29 closed. #18's body updated: the #29 row left the
      open table, **#31** and **#32** entered P4 with their Entry/closes-when, the `Closed 2026-10-10 —
      #29` paragraph, a row in the closed index and a line in "State to settle". Two new issues filed
      instead of absorbed: **#31** (group 26's six tests that are neither rows nor declared) and **#32**
      (the tool suite writing the real store — `#24` is what the rewrite loses, `#32` is the suite
      performing it). #24 got the measured record of the six damaged row files (comment `6098797801`)
      and keeps its status.

**Note on the push.** It carried 23 commits: this work unit's seven plus sixteen that were already
committed locally on `integrate/origin-main-2026-09-30` (the #30 session's series) and had not been
pushed. No unrelated working-tree change was committed: the shared tree's staged `CLAUDE.md ->
AGENTS.md` rename and the other session's four modified files are still uncommitted, untracked or
staged exactly as they were.

## Verification (independent, `gentle-ai-verify`, 2026-10-10)

All eight checks held: `check` 211/285/0/0; `status` 0 ungrouped; `coherence` on a **copy** of the store
stale for exactly `6-`, `20-`, `31-`, `34-`, `35-` with no extraction failure; `tools/tests -m "not slow"`
= 3 failed / 133 passed / 1 deselected, the three being the #23/#25 citation-site failures;
`regression --group 36` = `24 same`, exit 0; `check.sh quick` exit 0 (only the documented-red
`test_corotational_is_frame_objective_tl_is_not`); `git status --porcelain -- docs/validation` empty after
every run.

## Traps

- **Never run `coherence` against the real store** while working: it writes. `tools/tests`'s
  `test_every_group_re_derives_to_the_rows_on_disk` calls `run(REAL_STORE, "coherence")` and *does* write to
  `docs/validation` — this is #24, not something to reproduce. Copy the store to a scratch directory first.
- **The commit is pathspec-scoped.** The tree is shared with another session
  (`src/aeroelast/cli/run_bem_fsi.py`, `tests/validation/bem/test_wall_flow_activation.py`,
  `docs/blade_input_divergence_utd_vs_official.md`, `scripts/`, `.pi/` are dirty and not ours). Never
  `git add -A`.
- **Bootstrap every command that imports `aeroelast`:** `scripts/aeroenv.sh <cmd>`. `tools/validation_matrix.py`
  itself needs only PyYAML, but `collect_nodes` shells out to pytest.
- **This repository's output filter (`rtk`) drops and repeats lines.** Verify counts with `grep -c`, `wc -l`
  or a `Counter`; never trust a long listing, and always pass the explicit path to `git status`/`git diff`.
- The working tree's row files were restored from HEAD at the start of this work unit (the 6 files that a
  previous `coherence` run had rewritten: `6-`, `20-`, `31-`, `34-`, `35-`, `36-`).

## Handoff

- Store: `docs/validation/{groups,rows,gaps,sources,residual-patterns}.yaml|json`; policy in
  `docs/validation-policy.md`; design in `odd/tasks/validation-matrix-store.md`.
- Tool: `tools/validation_matrix.py` — `capture_prints:3585`, `command_regression:3879`, `command_extract:2863`,
  `coherence:3009`, `preserve_measurements:2800`, `compare_row:3750`, `tolerance_sites:2422`,
  `assertion_calls:2329`, `resolve_literal:2222`.
- Tests: `tools/tests/test_validation_matrix.py` (`_load_tool_module`, `make_row`, `write_store`, `run`,
  `REAL_STORE`).
- The 6 unmeasurable comparisons live in `docs/validation/rows/36-composite_ply_stress.yaml`; their present
  `notes` lines quote the old `unmapped` refusal and go away when the group is re-derived.
