# Contributing to AeroElast

Contributions are allowed and welcome even though the project itself is not
open-source in the OSI sense.

## What is welcome

- bug reports
- documentation improvements
- tests and reproducible cases
- solver fixes
- performance improvements
- new examples and tooling

## Before contributing

Please keep in mind:

- the project license allows public use only for non-commercial research,
  academic, educational, or evaluation purposes
- submitting a contribution does not grant you commercial usage rights to the
  project as a whole
- you should only submit code, data, or documentation that you have the right
  to contribute

## Contribution terms

By submitting any contribution to this repository, you agree that:

- you are legally allowed to submit the contribution
- the contribution may be reviewed, modified, rejected, or merged by the
  maintainer
- if accepted, the contribution becomes part of the project and may be used,
  distributed, relicensed, or sublicensed by the maintainer under the project
  licensing model
- public versions of accepted contributions will normally be distributed under
  the repository license unless explicitly agreed otherwise in writing

## Test rules (inviolable)

These govern every test in this repository. They exist because a test that cannot fail is
worse than no test: it certifies a result nobody measured.

1. **No tolerance above 5%.** A test comparing against a reference asserts inside 5% unless
   the reference's own scatter justifies more, and then the reason is written next to the
   tolerance. The suite's tolerance audit lives in `docs/validation-matrix.md` section 13.1.
2. **A tolerance is never chosen to accommodate the measured result.** It comes from the
   reference's own scatter or from a stated theoretical error, and it is written *before* the
   run. If a comparison misses the bound, the deliverable is the **finding** - a flagged
   residual, a documented limitation, an issue - never a wider window. A window picked so
   that today's number passes is exactly the failure mode this rule exists to stop.
3. **A test must be as restrictive as the physics allows.** When a magnitude comparison is not
   yet inside the bound, assert what the physics *does* pin - the sign, an invariant, a
   dominance relation, a self-consistency between two of our own models - and report the
   magnitude as a residual. Carry a **guard** that fails once the residual enters the bound, so
   the claim gets promoted to an asserted row instead of being left loose:
   `assert abs(1.0 - ratio) > 0.05, "the residual is now inside 5%: promote it to an asserted row"`.
4. **A sign or physical-sense assertion needs no tolerance**, and it is the strongest
   statement available while a magnitude is unsettled: "the twist is nose-down at rated" is a
   physics claim, and it is either true or false.
5. **The reference must be independent of the implementation under test**: a different code, a
   published table extracted from the paper itself, or a closed form. A formula re-implemented
   inside the test is not a reference.
6. **Every physics test row appears in the validation store** (`docs/validation/`) with its
   reference, tolerance, measured margin and flags, and the row counts must reconcile with
   `pytest -o addopts="" --collect-only`. A row that cannot fail is flagged there, never
   counted as evidence. The store is the rector: `docs/validation-policy.md` says what may be
   cited and `docs/adding-validation-tests.md` says how a row is added.
7. **A quoted source value is read, not extracted.** Every number, equation, table cell
   or figure value used as a reference is read from the rendered page with vision and
   recorded with its source key, printed page and locator. `pdftotext` is never the
   device for a value you will quote, and a citation without provenance is not evidence.
   The protocol is `docs/reading-sources.md`.

Rules 1-6 are the contract a test must satisfy; `docs/adding-validation-tests.md` is the
procedure that satisfies it, and `docs/reading-sources.md` is how a source value is read.

## Practical guidance

- keep changes focused
- include tests when fixing behavior
- document non-obvious assumptions
- avoid adding new dependencies without justification

## Legal note

If you need commercial use rights, a separate written agreement is required.