# Diagnosing a validation result

A symptom-indexed catalogue of the mistakes this project has already paid for, each with
the check that catches it. Read it when a number looks wrong, when a finding looks too
good, or before believing a clean run.

Rules of use:

- **A finding you cannot substantiate is not a finding.** If the check below cannot run,
  the verdict is `unverifiable`, never "stale", "wrong" or "fabricated".
- **Prefer a runnable check to an argument.** Every entry names one where one exists.
- **Never fix a symptom by widening a bound.** The deliverable is the finding
  (`CONTRIBUTING.md` rule 2).

## Before believing a clean run: what each command answers

`check` validates the store's **data** -- schemas, references, citations, literals, and the links
between rows and groups. It never runs a test and never reads the code, so it can be green while the
suite is red, and both are telling the truth: they answer different questions. The suite is the
authority on the physics; `check` is the authority on the store.

| the question | the command | what a clean answer means |
| --- | --- | --- |
| Is the store's data consistent? | `validation_matrix.py check` | the rows, the bibliography and the registries agree with each other |
| Does the physics pass? | `pytest`, on `tests/` | the comparisons held on this machine, in this run |
| Which row does this failure belong to? | `validation_matrix.py triage --scope <file>` | nothing: it names the failing test, the row it feeds, and how far past the bound it went |
| Did a stored number move? | `validation_matrix.py regression` | the recorded margins still match a fresh run |
| Do the row files still describe the code? | `validation_matrix.py coherence` | every group re-derives to the bytes on disk |
| What does the store cover, and what not yet? | `validation_matrix.py status` | every validation file is grouped, declared out of scope, or was deleted as software |

**A skip is not a pass, and not every error is a physics failure.** The suite skips when a
dependency or a reference deck is absent, and the message names which. When a skip claims a file is
missing, confirm it is not looking in the wrong place: thirteen tests skipped with "AeroDyn
reference not present" while the deck sat in `tests/reference/iea15mw_openfast`, because their
guard resolved the path from the test file's own depth and the file had moved. For the same reason
an `ERROR` is usually setup -- a file that is not where it is looked for -- and not a comparison
that failed.

**`status` renders one row per group.** `src` is how many source files the group declares, and it
reads `hand` where it declares none, because those rows were driven by hand and nothing can
re-derive them. `rows` and `cmp` are what the store holds for that group; `measured`, `near` and
`gt5` are about recorded margins and stay at zero until `regression --write` records them, which is
a run of the suite and not of `status`.

## Reading the source

### The residual is large, or a cell moved from the value the prose records

**Likely cause:** the wrong cell to begin with — MITC4 vs MITC4+ vs another code in the
same table, the wrong mesh (`N=8` vs `N=16`), or a different normalization
(Kirchhoff-normalized vs 3D-exact).
**Check:** re-read the cell with vision naming table, column, row and N
(`docs/reading-sources.md`), then diff it against the test's `expected` literal:
`python tools/validation_matrix.py extract --group <id>` prints the code's literals.
**Paid for:** section 3's twisted-beam and hemisphere-cutout cells were sourced from the
wrong tables and had to be corrected at `2fd847d`.

### An equation in the source does not match the implementation

**Likely cause:** the equation was read by `pdftotext`, which drops and reorders the
mathematics of these scans.
**Check:** re-render the page and read it with vision; validate the transcription against
a limiting case or a tabulated quantity.
**Paid for:** `pdftotext` placed Eq. (16) of Ko, Bathe & Zhang 2025 on a page where it is
not (`docs/formulations/mitc4plus-2017-extract.md`).

### A page citation points at the wrong page

**Likely cause:** PDF page number vs printed page number, an off-by-one from a
`pdftotext` page-finding, or a site that moved when the file was edited.
**Check:** `git log -L <line>,<line>:<file>` for a stale site; re-render the printed page.

### A citation site is reported stale

**Likely cause, in order:** (1) the line moved and the prose kept the old number;
(2) the mention is inside the declared range but not on its first line; (3) the entry has
no `code_mentions`, so the scan had nothing to search for and the verdict is vacuous.
**Check:** `python tools/validation_matrix.py references check`. An entry that declares
sites with an empty `code_mentions` is reported `unverifiable`, not stale: declare the
literal strings the code uses (`Ko2017`, `"Ko, Lee & Bathe (2017)"`).
**Paid for:** an early version of this checker turned an empty needle set into "the claim
is false" for ten entries.

### A work looks uncited, and the citation is in a manifest or a document

**Likely cause:** the scan reads `.py`, `.rs` and `.toml`. A site in a `.md` file is not
declarable today, because the bibliography itself names every work and would match every
entry trivially.
**Check:** grep the file by hand and record it as a finding for the maintainer.

## Choosing the reference and the tolerance

### A tolerance is above 5%, or was chosen after the run

**Check:** `python tools/validation_matrix.py list --gt5` lists every comparison above the
suite rule; each needs `tolerance.justified: true` and a written justification.
**Paid for:** a window existed for an 8.5% deviation on the thin twisted beam. The element
now measures inside the bound and the window is gone.

### The comparison cannot fail

**Symptoms:** both sides of the assertion come from the same formula; the assertion holds
whatever the solver returns; the value is printed and never asserted.
**Check:** `python tools/validation_matrix.py list --flag tautology --flag unfailable`,
and `--unit comparison` for comparisons with `asserted: false`. Those rows may be cited as
limitations, never as results.
**Paid for:** 28 tests that could never fail were removed from the suite, and the flags in
`docs/validation/flags.yaml` exist because such rows were counted as evidence first.

### A reference is the same formula re-implemented in the test

**Check:** read the test's reference against `CONTRIBUTING.md` rule 5: it must be a
different code, a published cell, or a closed form.

## Measuring

### The margin came from memory or from an older revision

**Check:** `measured.run` names the revision and `measured.date` the day; re-run the exact
node and compare. A margin nobody printed is `measured.status: not_printed`, not a number.

### The measurement is being compared against the wrong half of itself

**Likely cause:** one side is a node, the other is a maximum over all nodes, or the two are
normalized differently.
**Check:** the two sides must come from the same point of the model.

## The store and the checker

### `collect-only` returns 0 nodes

**Cause:** a module-level `pytest.importorskip` skipped the module, usually a missing
Rust extension or an unactivated environment.
**Check:**
`source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev`, then
`python -m pytest -o addopts="" --collect-only -q <file> | grep -c '::'`.
A zero is a failure, never a clean run: it reads as "zero drift".
**Paid for:** section 3 collects 31 nodes inside the pinned environment and 0 outside it.

### A capture command reports nothing where output clearly exists

**Symptom:** `validation_matrix regression` says `no node output captured` while the tests
pass.
**Likely cause:** the parser was written for the assumed output shape instead of the emitted
one. pytest prints a verbose node id with **no newline**, so the test's first print lands on
that same line and the status word comes after on a line of its own: `nodeid PASSED` on one
line never occurs.
**Check:** run one node with `-s -v` and read the last lines verbatim *before* writing the
parser, then pin the shape in a unit test with a synthetic output string
(`test_parse_prints_reads_the_shape_pytest_actually_emits`).
**Paid for:** the first regression parser matched zero nodes out of 31.

### A collected node is not claimed by any row

**Cause, in order:** a new test file with no group and no out-of-scope declaration; a
renamed test; a row whose `tests[]` names a node that no longer exists.
**Check:** `python tools/validation_matrix.py check`, which exits 1.

### A row claims a node that no longer collects

**Check:** `check` reports `tests[]` entries that name a missing file. Node ids are opaque
strings: never retype one, copy it from `--collect-only`.

### A derived flag (`near`, `gt5`) is stored in a row

**Cause:** someone wrote the flag instead of the numbers it is computed from.
**Check:** `check` refuses it. Delete it; the renderer recomputes it from tolerance and
margin.

### A tolerance is `null`, or a row was silently skipped

**Cause:** the tolerance is not a literal, so the extractor could not read it — it comes
from a parametrize parameter, a loop variable, or a helper that builds the case list.
**Check:** `validation_matrix extract` reports the node as `UNRESOLVED` and writes no row.
Map the parameter by position against the node id (`param_value_by_name`), or add the
missing link.

### An audit script accuses a source of fabrication

**Check first:** the script's join, not the source. Re-verify with an unambiguous key — a
DOI, a node id — before reporting anything.
**Paid for:** a comparison script paired entries with the wrong prose bullet and reported
"24 DOIs the source does not print". All 30 were present; the script was wrong.

### A pytest node id is being parsed with a split or a regex

**Cause:** ids join parameters with `-`, an exponent carries its own minus (`1e-06`), a
value may too (`In-plane`), and a parameter *name* can contain a number
(`expected_mitc40`).
**Check:** `param_tokens(params, expected)` in `tools/validation_matrix.py`, where
`expected` is the number of parameters the function declares; the tests in
`tools/tests/test_validation_matrix.py` pin each of those cases.

### A field required by the store is absent from the source

**Cause:** the source genuinely does not state it — a software version instead of a year,
a report with no author, a work with no printed title.
**Check:** declare it in `bibliographic_gaps` and leave the value empty. **Never** invent a
year, an author, a title or a DOI to satisfy the validator: an absence is recorded, never
filled in.
**Paid for:** a migration blocked on exactly this, and the tempting fix was an "audit year"
for software that states none. The validator was over-constrained, not the source.

### A DOI claim describes a copy that is not held

**Check:** `doi_status: printed_on_pdf` requires `held: true`. When the DOI was read from
a publisher page or a registry with no copy held, the state is `verified_externally` with
`verification: external_record`.
**Paid for:** `hughes1977` was stored as `printed_on_pdf` with `held: false`, which claims
a DOI printed on a copy nobody holds.

## Environment and tooling

### A row behaves differently on another machine

**Check:** `docs/validation-environment.md` pins the versions the rows were measured
under. CalculiX, OpenFAST and the polars are version-sensitive; a missing tool makes a row
`skip`, and a skip is not a pass.

### A generated view and the store disagree

**Check:** the store is the rector, always. `references render --check` exits 1 when the
export is stale. Never edit a generated file.
