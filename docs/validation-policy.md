# Validation policy

The rules that decide what this project may cite as validation evidence, and where every
number lives. The rows are in `docs/validation/`; this file is the policy that reads them.
Read `docs/adding-validation-tests.md` before adding a test and `docs/reading-sources.md`
before quoting a source.

This file replaces the prose that used to wrap the validation matrix. The matrix is gone
because it was a second copy of the numbers; a policy is not a copy of anything, so it can
stay.

## What counts as evidence

1. **A row in the store is evidence. Nothing else is.** A green suite is not evidence: a
   test that cannot fail certifies a result nobody measured, and `check` refuses the shapes
   that make that possible.
2. **The comparison is the unit.** One result against one independent reference, with its
   own tolerance and its own margin. A result measured against two references is two
   comparisons, never one cell with a slash in it.
3. **A number nobody printed is not a number.** `measured.status` says `not_printed` or
   `not_measured`; only `measured` carries a margin, and it carries the revision it was
   measured at.
4. **A comparison printed and never asserted is not a result.** Section 3 prints one
   against the 3D-exact reference values and asserts nothing; `validation_matrix regression`
   reports those as `informational`, with values that run to `622932%` and `6.2e8%`. It may
   be cited as a limitation, never as a result, and it is deliberately not stored: an
   unasserted comparison bounds nothing, so it has no tolerance to record and no slack to
   report.
5. **A gap is not citable.** Run `validation_matrix gaps` before citing anything that smells
   unverified. Every entry there is a claim the suite does not support, or a defect whose
   numbers must not be cited. The list is data and it is validated, because it is the part
   whose loss would otherwise be silent: the suite stays green and the caveat disappears.
6. **A reference must be independent.** A different code, a published cell, or a closed
   form. A formula re-implemented inside the test is not a reference
   (`CONTRIBUTING.md` rule 5).

## How tight a comparison has to be

- **Same method** — our code against another implementation of the same formulation, or
  against a closed form — must be tight. The suite rule is 5% (`CONTRIBUTING.md` rule 1),
  and a tolerance above it needs a written justification.
- **Different method or different author** — our shell against a beam model, our polars
  against another generator — is tight **on purpose**, so the difference gets measured,
  named and analysed instead of absorbed by a wide tolerance. Where the difference exceeds
  the bound it is declared as a validity bound, and the bound is quoted next to the number.
- **A green suite never unflags a row.** A row that cannot fail stays flagged, and a passed
  assertion does not promote it.

This is the whole of the validity envelope. What the suite may be trusted for is a
consequence of which comparisons exist, at which margins, and which gaps are declared — so
it is read from the store, never restated here where it would go stale.

## Where the numbers are

| what you want | the command |
| --- | --- |
| one row, verbatim | `validation_matrix get <id>` |
| the derived summary per row or comparison | `validation_matrix list --unit row\|comparison` |
| the comparisons closest to failing | `validation_matrix list --unit comparison --sort slack --near` |
| tolerances above the 5% rule | `validation_matrix list --gt5` |
| what is not validated | `validation_matrix gaps` |
| the headline per group | `validation_matrix headline` |
| did anything move | `validation_matrix regression` |
| the bibliography and its DOI audit | `validation_matrix references check` |

Do not hand-maintain a table of any of these. Two copies of the same number is how a
document manufactures a false defect, which is exactly what the matrix did: its first
version recorded an 8.5% deviation on a thin twisted beam that no longer existed.

## Before you quote a number

1. `validation_matrix get <id>` and read the tolerance next to the margin, never the margin
   alone.
2. `validation_matrix gaps` for the scope you are writing about.
3. If the number came from a paper, check its provenance in
   `docs/validation/references.yaml`; a citation without a page and a locator is not
   evidence (`docs/reading-sources.md`).
4. If a number is not in the store, it is not a result yet: add the test
   (`docs/adding-validation-tests.md`) or report the gap.
