# Adding a validation test

The contract for adding a physical test to the validation scheme. The *principles* are
`CONTRIBUTING.md` "Test rules" 1-7; this file is the procedure that satisfies them, and
the one place they are applied. Read `docs/reading-sources.md` before quoting a value
from a source: it is the input half of this contract.

Applies to every test whose result becomes a matrix row, not only to "validation" tests
by name.

## A row has two halves, and both must be real

| half | what it must be |
| --- | --- |
| the reference | independent of the implementation under test, read from the source with vision, with its page and locator |
| the measurement | a residual the test actually prints, from a real run, at a named revision |

A row whose reference is the formula re-implemented inside the test is not a row
(rule 5). A row whose margin nobody printed is `measured.status: not_printed`, not a
number.

## What belongs to this scheme

A row exists for a test that validates a **physical quantity against an independent
reference**: a displacement, a stress, a natural frequency, a mass, a buckling or
stability limit, a constitutive response, a load. If the test asserts that a mesh loads,
that the CLI writes a deck, that a schema rejects a bad input, or that a logger records a
field, it is a **software test**: it belongs to the suite, it must pass, and it does not
belong here. Declare its file out of scope in `docs/validation/groups.yaml` with a reason,
rather than giving it a row.

The one-sentence test: name the physical quantity and its independent reference. If you
cannot, it is a software test.

Why this boundary is strict: a row that cannot be wrong about physics dilutes the rows
that can. The matrix's own history is the argument — rows comparing a value against a
schema, or asserting nothing that could fail, had to be flagged one by one and finally
excluded, and the flag set that marks them (`unfailable`, `tautology`, `info_only`) exists
only because they were counted as evidence first.

A declared-software file needs no row, but the declaration is required: `check` must be
able to tell *out of scope* from *forgotten*, and a collected file that is in neither a
group nor the out-of-scope list is an error, not a gap.

## Steps

1. **Pick the reference and read it.** Read it from the rendered page with vision and
   record the source key, the printed page and the locator
   (`docs/reading-sources.md`). If the work is not in
   `docs/validation/references.yaml` yet, add it there first; the key is what a row
   cites.
2. **Choose the tolerance before running.** It comes from the reference's own scatter or
   from a stated theoretical error, never from what the result turns out to be (rule 2).
   Anything above 5% needs `tolerance.justified: true` and a written
   `tolerance.justification` (rule 1).
3. **Write the test as restrictively as the physics allows.** When a magnitude is not yet
   inside the bound, assert what *is* pinned — the sign, an invariant, a dominance
   relation, a self-consistency between two of our own models — and report the magnitude
   as a residual (rule 3). A sign or physical-sense assertion needs no tolerance and is
   the strongest claim available (rule 4).
4. **Print the residual.** The measured margin is transcribed from the test's own
   output. An assertion that cannot print one yields `not_printed`, and the row says so.
5. **Declare the group, or declare the file out of scope.** A new test file needs an entry
   in `docs/validation/groups.yaml` with its `id`, `slug`, `source_files` and the
   `citation` key its paper references use — or, if it is a software test, an out-of-scope
   entry with the reason. A collected file in neither list is reported as drift by
   `check`.
6. **Seed the rows from the code:**
   `python tools/validation_matrix.py extract --group <id> --write`. It claims every
   collected node exactly once and fills what the code states — the tolerances, their
   sources, and the references it can label. It writes nothing it cannot read, and
   reports such a node instead of inventing a value.
7. **Complete the comparisons by hand where the code cannot state them.** Each
   independent reference gets its own entry: `reference.label` carries the locator and
   the cell, `expected` the reference value, `tolerance.source` where the code sets the
   tolerance, `measured.text` the verbatim printed residual. One test measuring against
   CCX and against a closed form is two comparisons, not one cell with a slash in it.
8. **Record the measurement** with `measured.status: measured`, `margin_pct`, `run` (the
   revision) and `date`; `not_measured` when the run did not happen. A comparison the
   test prints and never asserts is `asserted: false`, and the row then carries the
   `info_only` flag so it is never quoted as a result.
9. **Prove the reconciliation:** `python tools/validation_matrix.py check --group <id>`
   must exit 0, and the row count must match
   `pytest -o addopts="" --collect-only -q <file>` (rule 6). A zero-node collection is a
   failure, not a clean run: it usually means the module skipped for a missing
   dependency.
10. **Do not hand-edit the generated views.** `docs/references.md` is generated from the
    store (`references render --check` in CI). The matrix view becomes generated at T7;
    until then its Markdown is still hand-maintained, and the store is the rector for
    section 3.

## What `check` will refuse

These are rejections, not warnings, and each one is a rule above made executable:

- a stored derived flag (`near`, `gt5`) — they are computed per comparison;
- a comparison with `asserted: false` without the row's `info_only` flag;
- `measured.status: measured` without a numeric `margin_pct`, or `margin_pct` present
  when the status is not `measured`;
- `reference.kind: paper` without a `citation` that resolves in the store;
- a missing required bibliographic field that is not declared in `bibliographic_gaps`,
  or a declared gap whose field is actually present;
- a tolerance above 5% with no written justification;
- a comparison whose `reference.kind` is `schema`: a row with no independent reference is a
  software test, and there is no such row here;
- a `tests[]` entry that is not an exact collected node id, or names a missing file.

## The Rust side

The crates hold roughly 201 `#[test]` functions that validate against published work and
that the matrix does not model at all. Until a Rust group exists, this contract applies
to them in principle and cannot be applied in practice: a Rust row needs a node-id scheme
that expresses a cargo test, and the reconciliation step above is pytest-specific.
Treat a Rust validation as documented-but-not-inventoried, and never as a row.

## Do not

- Re-implement the reference inside the test (rule 5).
- Pick a tolerance so that today's number passes (rule 2).
- Count an assertion that cannot fail, or a value printed and never asserted (rule 6).
- Quote a number you did not see on a rendered page (`docs/reading-sources.md`).
- Widen a window instead of reporting the finding. The deliverable is the finding.
