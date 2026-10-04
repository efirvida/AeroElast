# The validation store: workflow for an agent

Read this file to work on validation evidence in this repository. It is written to be executed,
not admired: every command below runs, and every number is the one the tool prints today.

**The one rule to internalise first.** `check` validates the store's **data** and never runs a
test. A green `check` means the store agrees with itself; it never means the physics passes. The
suite is the authority on the physics, and `triage` is the bridge between the two.

## 1. The model

- A **row** is one own result (`result`, optional) plus N **comparisons**. Never one cell with a
  slash in it.
- A **comparison** is a canonical call in a test:

  ```python
  assert_relative_error(value, reference, tol=TOL, kind="paper", reference_name="...", what="...")
  assert_residual_below(residual, tol=TOL, kind="self", reference_name="...", what="...")
  assert_residual_below(residual, atol=BOUND, unit="deg", kind="code", reference_name="...", what="...")
  ```

  `tol` is a fraction (a ratio tolerance must be `<= 1.0`); `atol` is a quantity in `unit`. Both
  helpers print in one shape: `{what} vs {reference_name}: {error} (bound: {bound})`.
- **Reference kinds** (`kind`): `paper`, `analytical`, `code`, `self`. `paper` requires a
  `citation` that resolves in `references.yaml`; `citation` is refused on every other kind. The
  kind is read from the comparison, never inferred from the file or its prose.
- **The site criterion**: only a canonical call is a site. A bare `assert err < 0.05` is **not**
  read: it names no reference, so the store cannot attribute it. Convert it or leave it alone --
  do not declare it.
- **Digest vs delta**: the digest identifies *who changed* (it hashes the printed `value`,
  `error`, `expected` strings); the stored margin measures *by how much*.

## 2. File map

| path | holds |
| --- | --- |
| `docs/validation/groups.yaml` | the registries: one group per validation file, plus `out_of_scope`, `non_validation_tests`, `non_reference_asserts`, `validation_helpers`, `reference_kind`, `citation` |
| `docs/validation/rows/<group>-<slug>.yaml` | one file per group, the rows |
| `docs/validation/references.yaml` | the bibliography, 66 entries; `doi_status` per entry |
| `docs/validation/gaps.yaml` | the 6 declared gaps: what the suite does **not** validate |
| `docs/validation/residual-patterns.json` | per-group regexes (`asserted`, `unasserted`) under `patterns` |
| `docs/validation/sources.json` | per-group source digests |
| `docs/validation/flags.yaml` | stored flags |
| `docs/validation/adjudications/3-ko2017.yaml` | the pilot's adjudication record |
| `docs/validation/policy` docs | `validation-policy.md` (rules), `adding-validation-tests.md` (the contract), `validation-diagnostics.md` (symptom index), `validation-environment.md` (pinned versions) |
| `tools/validation_matrix.py` | the tool. `tools/tests/test_validation_matrix.py` is its suite |

## 3. Commands

| command | answers | notes |
| --- | --- | --- |
| `check [--group N]` | is the store's data consistent? | the only command that needs nothing but PyYAML |
| `status [--json]` | what does the store cover, per group? | one row per group; see §5 for the columns |
| `gaps` | what the suite does not validate | prints the 6 declared gaps |
| `references (list\|gaps\|get\|check\|bibtex\|where-used\|render)` | the bibliography | `references check` = 66 entries, 0 errors, 13 warnings (DOI debt) |
| `find` / `list` / `get` / `headline` | query rows | `list --gt5`, `list --group N` |
| `extract --group N [--scope PATH] [--write]` | derive rows from the code | needs the built env; without `--write` it is a dry run |
| `regression --group N [--scope] [--write]` | did the printed numbers move? | needs `residual-patterns.json` for the group; `--write` records margins |
| `coherence [--group N]` | do the row files still match the code? | re-derives every group and diffs in memory |
| `triage --scope PATH [--first] [--json]` | which test failed, in which row, by how much | runs the scope |
| `kinds` | undeclared reference kinds | review sheet |
| `set` | edit one row in place | refuses an edit that introduces an error outside `rows/` |

Exit codes: `0` ok, `1` findings, `2` error.

## 4. Invariants `check` enforces

- `reference.kind` is mandatory and read from the comparison.
- `paper` needs a resolving `citation`; `citation` is refused on other kinds.
- `bibliographic_gaps` may only contain `authors, pages, title, venue, volume, year`. **A DOI is
  not declarable**: an unknown DOI is `doi_status: to_verify`.
- A ratio tolerance (`rtol`, `rel_err`) above `RATIO_TOLERANCE_MAX = 1.0` is refused.
- `expected` must not be a numeric string; `expected` is a number only when the reference side is
  one value.
- A row may only claim a test under `tests/validation/`.
- A declaration that matches nothing is a finding (the `extract` reports it, `check` cannot).
- Constants: `SUITE_TOLERANCE_RULE = 0.05`, `RATIO_TOLERANCE_KINDS = {rtol, rel_err}`,
  `RATIO_TOLERANCE_MAX = 1.0`. Declared pytest markers, for `-m`: `slow`, `benchmark`,
  `composite` (there is no `fast`).

## 5. Reading `status`

```text
group      src  rows   cmp measured near  gt5
3            1    31    37       37    4    0
10        hand     1     4        4    0    0
```

- `src` = how many source files the group declares; **`hand`** = none, its rows were driven by
  hand and nothing can re-derive them. Not a gap.
- `rows` / `cmp` = what the store holds for the group.
- `measured` = comparisons carrying a margin. Zero until `regression --write` runs.
- **`near` is about the margin**: the margin is within 1 point of the bound.
- **`gt5` is about the tolerance**: the declared bound is looser than `5 * SUITE_TOLERANCE_RULE`,
  i.e. above 5%. It does **not** mean "past the bound five times over". It fires on comparisons
  with no measurement at all: a **normalised residual** legitimately carries `tol=1.0` because the
  number is already a fraction of its own allowance, so a `gt5` there is the shape, not a loose
  window.

## 6. Gaps

`gaps.yaml` holds claims the store **cannot derive**: absence is not in the data, so a quantity
nobody validates is invisible unless declared. Two statuses:

- `not_validated` -- nothing in the suite validates this quantity.
- `open_defect` -- a defect is open, measured, with its evidence and its consequence.

They are **not generated automatically**. Today 5 of the 6 were migrated by hand from the old
`docs/validation-matrix.md` not-citable list, and 1 was found during the sweep. The gap registry
is therefore **pilot-era**: the sweep that gave 28 groups their rows mined rows and never asked
the complementary question, so most domains have no gap recorded. Signals worth mining by hand:

- an `xfail` marker (the suite admits a defect): 12 exist today, none in `gaps.yaml`;
- a test that prints a verdict it never asserts (group 3 prints `[x] Norm vs Paper 3D: ...`
  under a comment saying "informational only");
- a comparison the store reads whose bound is loose (`list --gt5`).

Use `gaps` to read them; edit `gaps.yaml` to add one; no command derives them.

## 7. Recipes

### A. Migrate a validation file into a group

```bash
# 1. read the file: what does it compare, against what, with which bounds
# 2. convert its comparisons to the canonical calls (recipe B)
# 3. add the group to groups.yaml: id, title, slug, source_files, provenance_note
# 4. add a residual pattern (recipe C)
# 5. extract --group N --scope <file> --write
python tools/validation_matrix.py check                      # must be 0 errors
python tools/validation_matrix.py extract --group N --scope <file>   # dejad: nothing undeclared
python -m pytest -o addopts= -q <file>                       # the file still passes
python tools/validation_matrix.py regression --group N --write        # record the margins
python tools/validation_matrix.py coherence --group N        # the rows match the code
```

Rules that decide whether an assertion becomes a row:

- A comparison the criterion reads, against an independent reference -> a row.
- A property of our own matrix (a symmetry, a numerical zero, a sign the theory predicts) ->
  **not** a row: write it as a bare `assert` so it stops looking like a comparison. If it is a
  call that carries a bound while both sides are ours (`assert_allclose`), the store still reads
  it and it needs `non_reference_asserts`.
- A whole test that produces no row (mapping, plumbing, binding) -> `non_validation_tests`, one
  entry per test name with a reason.
- A test whose comparison lives in a same-module helper -> `validation_helpers`.

### B. Convert a comparison

| before | after |
| --- | --- |
| `assert err < 0.05` | `assert_residual_below(err, tol=0.05, kind=..., reference_name=..., what=...)` |
| `assert_allclose(a, b, rtol=r)` | `assert_residual_below(max(abs(a-b)/abs(b)), tol=r, ...)` (exact when `b != 0`) |
| `assert_allclose(a, b, rtol=r, atol=atol)` | `resid = max(abs(a-b) / (atol + r*abs(b)))` with `tol=1.0`; below 1 means inside both bounds |
| `assert_allclose(a, b, atol=atol)` with zero entries | `max(abs(a-b)) / max(abs(b))`, `tol` from the same ratio |
| a bound that is an absolute quantity | `atol=BOUND, unit="deg"` (a `%` is not a ratio; it needs `atol`) |
| a class-level bound `self.TOL` | hoist it to a module constant: the extractor resolves module and class level, but not `self` |

The frequency case that needed `unit="%"`: a bound of 5 percent is `atol=5.0, unit="%"`.

### C. Add a residual pattern

`regression` reads the test's stdout through `patterns[<group>]`, keyed by group id:

```json
"7": {
  "asserted": "(?P<what>.+?) vs (?P<reference_name>.+?): (?P<error>\\S+)% \\(bound: (?P<bound>\\S+)%\\)",
  "unasserted": "..."
}
```

- **Measure the print before writing the regex.** `capture_prints(scope)` in the tool returns the
  real lines; a guessed pattern silently matches nothing and the comparison stays silent.
- The canonical shape above matches every comparison written through the suite's helpers, so one
  pattern usually serves a whole group.
- A print that matches the `unasserted` pattern becomes an `INFO ... (never asserted)` line, and
  its stored text is the printed line verbatim.

### D. Record margins

```bash
python tools/validation_matrix.py regression --group N --write
```

- Needs the built environment (it re-runs the scope) and a pattern for the group.
- Writes `measured` per comparison plus a source digest. A comparison whose print the pattern did
  not read is left unmeasured and reported; it never stores a null (`check` would reject it).
- `measured` therefore means *this machine, this date, this revision*.

### E. Verify

```bash
python tools/validation_matrix.py check                  # data consistency
python tools/validation_matrix.py coherence              # row files vs the code, byte for byte
python tools/validation_matrix.py extract --group N      # nothing undeclared, no stale declaration
python -m pytest -o addopts= -q -m "not slow"            # the suite (27 min)
python -m pytest -o addopts= -q -m "not slow" tools/tests/test_validation_matrix.py   # the tool's suite
```

- Compare **bytes**, not counts: an earlier count-based audit passed while row files cited lines
  that had moved. `coherence` is that audit as a command, and `extract --write` + `git diff` is
  the same thing by hand.
- A green suite and a green `check` can coexist with 39 errors, which happened: the tests had
  anchored data paths to each file's own depth and the `tests/` reorganisation moved them.
  **Running the suite is the only instrument that sees this class.**

### F. Declarations, by level

| level | mechanism | when |
| --- | --- | --- |
| file | `out_of_scope` in `groups.yaml` (files + a reason) | a validation-tree file that is not a validation |
| test | `non_validation_tests` (test names + a reason) | a test that produces no row, and *why* |
| comparison | `non_reference_asserts` (`<path>:<line>` + a reason) | a call-shaped comparison with no independent reference |
| helper | `validation_helpers` | the comparison lives in a same-module helper |
| absence | `gaps.yaml` | nothing validates this, or a defect is open |

Ask **which level** before writing: a file-level or test-level declaration is often the right
answer where a list of line exclusions would be, and it says more.

## 8. Traps (each one measured, not theorised)

1. **`check` never runs tests.** Green `check` + red pytest is not a contradiction.
2. **A false skip is worse than a failure.** A skip that says "dependency absent" while the
   dependency is present turns off evidence. A missing tool skips; a tool that is present and
   refuses its input must stay an error.
3. **Do not infer a tool's behaviour from one artefact.** `assert_allclose` *is* read (as a
   comparison with no reference); `extract` without `--group` covers one scope, not all groups;
   the `gt5` column is about the tolerance. Each of those was inferred wrongly first. Run it.
4. **The editors' summaries strip indentation.** Read with `repr`, `read`, or `git show` before
   reasoning about whitespace.
5. **A moved line leaves a row file stale.** Nothing else sees it: `check` cannot read the code.
   Run `coherence` after editing a test.
6. **Declarations and data belong to different commits.** `git add tools docs` sweeps row files
   into a code commit. Use exact paths and read `git show --stat`.
7. **The extractor resolves constants, not `self`.** A bound on a class attribute is not read.

## 9. Pending work, recorded

- **The Rust crates.** `crates/aeroelast-{core,mesh,py,solvers}` hold tests validated against
  papers (`elements/{mitc4,mitc3,quad,smoothing}.rs`, `materials/laminate.rs`, `py/src/{materials,
  elements}.rs`). The store covers `tests/` only, so they are outside it. Including them is a
  design decision: the extractor reads pytest nodes today.
- **Gaps per domain.** See §6: the registry is pilot-era and the sweep produced none.
- **Comparisons the patterns do not match.** `status` shows them per group: the `measured` column
  is below `cmp`. Group 5 prints a list of gaps, group 11 prints the OpenFAST lines, group 14
  prints `physical shear window` instead of a bound.
