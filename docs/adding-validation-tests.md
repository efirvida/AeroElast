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
belong here. Put it under `tests/software/`, rather than giving it a row or a declaration.

The one-sentence test: name the physical quantity and its independent reference. If you
cannot, it is a software test.

**The directory is that answer, and it is enforced.** Write a validation test under
`tests/validation/<domain>/` and a software test under `tests/software/<domain>/`; shared
helpers go in `tests/support/`, which is a package so any test can import them by name. A row
may only claim a test under `tests/validation/`, so `check` fails on a row for a software
test, and the two cannot drift apart. There is no list of excluded files to maintain: the path
is the classification.

Why this boundary is strict: a row that cannot be wrong about physics dilutes the rows
that can. The matrix's own history is the argument — rows comparing a value against a
schema, or asserting nothing that could fail, had to be flagged one by one and finally
excluded, and the flag set that marks them (`unfailable`, `tautology`, `info_only`) exists
only because they were counted as evidence first.

A declared-software file needs no row, but the declaration is required: `check` must be
able to tell *out of scope* from *forgotten*, and a collected file that is in neither a
group nor the out-of-scope list is an error, not a gap.

## The shape the store can read

A validation test stays an ordinary pytest test: it runs, it asserts, it fails when the physics is
wrong. What is specified here is the *form* of one assertion, so that the store reads the comparison
instead of inferring it.

**A comparison goes through the suite's assertion helper, and the helper is told what the reference
is.**

```python
assert_relative_error(tip_disp, REF_TABLE_6, tol=TOL_REFERENCE, reference="Ko et al. 2017 Table 6")
```

The reference argument is not decoration: it is the only thing that says whether the two quantities
in the comparison are our result and an independent reference, or two of our own numbers. No
structural rule can tell those apart, which is why the store grew a list of per-test declarations
instead of reading the tests. With the reference named, it does not need any:

- **the extractor reads the helper calls.** A test that makes one is a comparison and becomes a
  comparison in a row. A test that makes none is not a row, and that is a fact about its shape, not
  something a maintainer has to declare.
- **every other assertion in a validation test is a physics property or machinery.** A symmetry, an
  invariant, a dominance relation, a torque ruler, a setup constant: all of them stay in the test
  and keep asserting exactly what they assert today. None of them is a comparison against an
  independent reference, so none of them is a row.
- **assertions in a local helper are invisible.** Put the comparison where the test is; a helper may
  build, solve and print, but the assertion that becomes a row is written in the test body.
- **one statement per reference.** Two references are two helper calls, never one cell with a slash
  in it.

The store's per-test declarations -- `non_reference_asserts`, `validation_helpers`,
`non_validation_tests` -- exist because the tests did not say which of these they were. They are
scaffolding for the files written before this contract, and they go away as those files conform.
Twenty `non_reference_asserts` entries went away on their own when the criterion stopped reading
bare relational asserts: `assert err < 0.05` names no reference, so it is not a site and needs no
declaration. The declaration is still the right one for a single case, a comparison written as a
call that carries a bound while both sides are our own -- `assert_allclose` of a quantity against a
literal -- which the store does read and has to be told not to attribute.

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
5. **Put the file in the right directory, and declare its group if it is validation.**
   `tests/validation/<domain>/` needs an entry in `docs/validation/groups.yaml` with its `id`,
   `slug`, `source_files` and the `citation` key its paper references use. A test under
   `tests/software/<domain>/` needs none of that: the directory already says it is not
   evidence.
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
10. **Do not hand-edit the generated views.** `docs/validation/references.yaml` is generated from the
    store (`references render --check` in CI). Both Markdown views were deleted at T7, so there is no
    longer a hand-maintained copy of anything: the store is the rector, and the only file that still
    renders is the bibliography.

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
