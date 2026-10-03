# Feature: a machine-queryable store behind `docs/validation-matrix.md`

## Objective

Replace `docs/validation-matrix.md` and `docs/references.md` as the *storage* of validation
evidence with structured, schema-validated stores, keeping both Markdown files as **generated**
artifacts so that every existing cross-reference (CONTRIBUTING rules 1 and 6, `docs/*.md`,
`openspec/specs/`, the `tests/test_shell_convergence.py` comment, the traceability gate in
`tests/test_mitc4plusd_traceability.py`) keeps resolving against the same paths.

Success is measured in three ways:

1. An LLM asked to update one row reads **one record (~150 tokens)** instead of scanning a
   196 KB / ~50k-token document, and never edits the same number in two places.
2. `pytest --collect-only` reconciliation (CONTRIBUTING rule 6) becomes a **command** with a
   non-zero exit code, not a hand-maintained inventory table.
3. §2.1, §13, §13.1, the index, the flag views and the reference bibliography are **derived**,
   so they cannot disagree with the per-row data.

## 1. Status and decisions taken

- **Phase: design only.** No data migration, no code, no edit to `docs/validation-matrix.md`
  or `docs/references.md` until the pilot in §12 is authorized.
- Maintainer decisions, this session:
  1. Store format: **YAML, one file per group.** Comments are needed to carry tolerance
     justifications and git diffs must stay readable.
  2. **References get their own store.** `references.yaml` is the rector for structured,
     auditable, queryable reference fields (including the DOI); `docs/references.md` becomes
     generated from it.
  3. **Build from the tests, not from the Markdown.** The store is populated by static
     extraction from `tests/` plus a live run, then **cross-validated against** the current
     Markdown. Every disagreement is adjudicated by re-running, and the **new schema wins**:
     it becomes the authority for that group, and the authority for everything when the full
     migration lands.
  4. **Pilot group: §3, `tests/test_ko2017_performance.py`** (31 tests, 23 Markdown rows,
     2132 lines, one file).
  5. **The references migration ships with the pilot** (WU-2), not as a separate wave.
  6. **The pilot is the minimal production test, not a spike.** Nothing in it is throwaway and
     nothing is marked temporary: the schemas, the CLI subcommands (including the extractor),
     `references.yaml`, `rows/3-ko2017.yaml` and the adjudication log are the artifacts WU-3,
     WU-4 and WU-5 will reuse unchanged. The pilot's job is to prove the pipeline end to end on
     one group and leave the remaining groups as mechanical repetition of it.
  7. **Eliminate the data-shaped Markdown.** The end state has no generated `.md` at all:
     `docs/validation-matrix.md` and `docs/references.md` are deleted, and the store plus its
     query verbs replace them. Only *theoretical* documents (`docs/formulations/*`) and
     *general* documentation (README, CONTRIBUTING, `docs/reading-sources.md`,
     `docs/adding-validation-tests.md`, `docs/validation-environment.md`) survive. The policy
     prose the matrix carries — what may be cited as evidence, and the validity envelope —
     becomes a hand-written `docs/validation-policy.md`; everything else it contains is either
     already data in the store or derived by the tool. `references render` survives as an
     export/debug command only.
     **Neither file is edited once deletion is decided.** A pointer added to a document that is
     about to be deleted is churn: it lands in the review candidate, costs a diff, and then
     disappears with the file. The matrix is *read* (for T5's diff) and deleted; consolidation
     edits target only documents that survive. In particular the reading rules that the matrix
     restated already live in `docs/reading-sources.md`, so removing the restatement loses
     nothing — including the record of the `pdftotext -layout` failure, which survives there.
  8. **The scheme covers physics only.** A row exists for a test that validates a physical
     quantity against an independent reference. General software tests — CLI, mesh I/O,
     logging, schema and contract tests — are declared *out of scope* with a reason, not given
     rows. `reference.kind: schema` therefore leaves the model, and the `un_inventoried` list
     stops being a debt bucket and becomes a classification.
  9. **The surviving documents are prescriptive, and they carry their errors.** The point of
     `docs/reading-sources.md`, `docs/adding-validation-tests.md` and
     `docs/validation-diagnostics.md` is to stop an agent hallucinating and to stop it paying
     for the same mistake twice. Each one states the forbidden move, the runnable check that
     catches it, and the symptom that reveals it — not a description of the current state,
     which belongs in the store where it cannot go stale. `docs/validation-policy.md`, written
     at T7, holds the citation policy and the validity envelope and points at the diagnostic
     index. The test for any sentence in these files: does it change what the reader *does*?
  10. **`tests/` is organized by role, then domain.** `tests/validation/<domain>/` holds a test
     that validates a physical quantity against an independent reference;
     `tests/software/<domain>/` holds everything else; `tests/support/` holds `conftest.py` and
     the shared helpers, so that imports across directories resolve (today several test modules
     import each other by bare name, which works only while every file sits in one directory).
     The location **is** the classification, which is why the hand-written `un_inventoried`
     bucket is deleted when the move lands rather than refined now: a classification written
     today would be replaced by structure tomorrow. Measured cost: 56 files moved and 142
     `tests/…` path references outside `tests/`. The move happens **after** T7 deletes the
     Markdown, because `diff-against-md` finds its table by the heading that names the source
     file and would break on a rename.

## 2. Why the current files are not maintainable

Measured on the working tree at the time of writing (`docs/validation-matrix.md`: 196 288 bytes,
1757 lines; `docs/references.md`: ~586 lines).

The same facts are stored four times in four different shapes:

| representation | content | relation to the others |
| --- | --- | --- |
| §3-§8, 15 tables, ~240 rows | `test \| what it validates \| reference \| tolerance \| measured margin \| notes` | row data |
| §13, 252 rows | `§ \| test \| tolerance \| measured margin \| slack to fail \| flags` | **duplicate** of tolerance + margin |
| §2.1, 108 rows | group headline + worst measured result | aggregate of the same numbers |
| §9, 6 tables, ~54 rows | flags by defect class, with evidence prose | metadata of the row data |

Concrete defects this duplication has already produced, all admitted inside the document:

- Stale cells. "A stale row here is how this document can manufacture a false defect: the
  first version recorded an 8.5% deviation on the thin twisted beam that no longer exists."
  §4.8 also carries "Per-file drift, corrected rather than left implicit" for eight files.
- A standing drift table of 7 files / 68 tests that are collected but have no row, reconciled
  by hand with two shell commands.
- Row identity is prose: `...[dist SS]`, `(2 cases)`,
  `test_3_1_square_plate_tables_2_to_5[reg clamped, t/L=1/100..1/10000] (3 cases)`. No key
  joins against `pytest --collect-only`, which is exactly what CONTRIBUTING rule 6 requires.
  For §3 the mismatch is measurable today: **31 collected nodes in 23 Markdown rows.**
- Nothing is checkable: no schema, no validator, no duplicate detection, no
  "flag without a live test" detection.
- Updating one measured margin means editing the §3-§8 row *and* the §13 row, and re-deriving
  `slack to fail` and the flags by hand.

`docs/references.md` has the same shape of problem on a smaller scale. Entries are unkeyed
prose bullets, and two of their fields are claims about the code that nothing verifies:

- `DOI: to verify` vs `DOI: 10.1016/j.compstruc.2016.11.004` — the audit status is a sentence
  fragment inside a bullet, so there is no way to list the DOI gaps or fail CI when a
  verifiable DOI is left unverified.
- `*Cited by the code: src/aeroelast/core/assembler.py:681, ...*` — a hand-written claim that a
  given file:line cites the work. Nothing checks that line 681 still cites it, and the
  convention distinguishes "Verified against the recovered PDF" from "Source: repository
  citation" purely in prose.
- The matrix's §14 map points at references.md **by section number** (`§5`), so inserting an
  entry renumbers the target and silently breaks the map.

The root cause is one sentence: **the Markdown is simultaneously the store and the view.**

## 3. Target architecture

```text
docs/validation/
  references.yaml             # canonical bibliography: structured + auditable
  groups.yaml                 # group registry: id, slug, citation, source files
  flags.yaml                  # declarative flag registry (the section 9 content)
  gaps.yaml                   # what the suite does NOT validate (not derivable)
  residual-patterns.json      # how each group's tests print their residual
  schemas/
    validation-row.schema.json
    reference.schema.json
  adjudications/
    3-ko2017.yaml             # the audit trail of a reconciled group
  rows/
    3-ko2017.yaml             # the pilot's rows: one per collected node
    10-out-of-band.yaml       # evidence with no live test
tools/validation_matrix.py
docs/validation-policy.md     # the surviving policy; the Markdown views are deleted
```

Rationale for one file per group: an LLM updating §4.8 opens `rows/4-ccx-parity.yaml` and still
sees 87 rows. If that is still too much, the CLI is the intended access path (§9) and the file
layout is only a fallback. A per-row-file layout is rejected: 240 tiny files make `git log` and
review worse, not better.

## 4. References store

`docs/validation/references.yaml` is the rector. One record per work. `docs/references.md` is
generated from it, grouped by `section`, with the verification preamble kept as prose.

```yaml
- key: ko2017_perf                     # stable, ASCII, used by matrix rows and by tests
  section: "1. Shell element formulations"
  kind: journal                        # journal | conference | report | manual | software | thesis
  authors: ["Ko, Y.", "Lee, Y.", "Lee, P.-S.", "Bathe, K.-J."]
  title: "Performance of the MITC3+ and MITC4+ shell elements in widely-used benchmark problems"
  venue: "Computers and Structures"
  volume: "193"
  pages: "187-206"
  year: 2017
  doi: "10.1016/j.compstruc.2017.08.003"
  doi_status: verified                 # verified | printed_on_pdf | to_verify | not_applicable
  held: true
  held_files:
    - ".sources/papers/1-s2.0-S0045794917309550-main.pdf"
    - ".sources/papers/Performance_of_the_MITC3+_and_MITC4+_shell_elements_in_widely_used_benchmark_problems.pdf"
  verification: verified_pdf           # verified_pdf | repository_citation | unverified
  verification_note: "authors, title, journal, volume, year; DOI printed as http://dx.doi.org/10.1016/j.compstruc.2017.08.003"
  cited_by_declared:                   # the AUDITED claim, carried over from the prose today
    - "tests/test_ko2017_performance.py:4-7"
    - "src/aeroelast/core/mesh/generators.py:65"
  code_mentions:                       # literal strings the code uses to cite this work
    - "Ko2017"
  notes: "the README entry omitted Bathe"
```

Field rules:

- `key` is the join key: a matrix row's `reference.citation` holds it, and the §14 map is
  generated from it. `check` fails if any row cites an unresolvable key.
- `code_mentions` lists the **literal strings the code uses** when citing the work. It is not
  decoration: the code cites author-year text (`Apply Ko2017 ratio-based mesh distortion` at
  `src/aeroelast/core/mesh/generators.py:65`) and never the store's key, so a key-only scan would
  report every declared site as a false mismatch. With `code_mentions` the scan is exact, and a
  declared site that contains none of them is a genuine finding: the code stopped citing the work
  at that line.
- `doi_status` replaces the prose marker `DOI: to verify`. `verified` and `printed_on_pdf`
  require a conforming DOI (`^10\.\d{4,9}/\S+$`) **and** `held: true` **and** a non-empty
  `verification_note`. `not_applicable` covers NAFEMS-style report numbers, proprietary
  manuals and journals that register no DOIs, and requires a `notes` explaining why.
- `cited_by_declared` keeps today's audited claim verbatim. It is deliberately *not* the
  same field as what the code actually cites: see `cited_by_found` below. Both exist so the
  difference is reportable.
- `cited_by_found` is **derived**, never stored: the CLI scans `src/` and `tests/` for the
  `key` (and, during migration, for author-year strings) and reports `file:line`.
- `held_files` lives under `.sources/`, which is gitignored. `check` must therefore warn, not
  fail, when a fresh clone has no `.sources/`: absence is expected, presence is a bonus.

Query surface:

```bash
validation_matrix references list [--doi-status to_verify] [--kind journal] [--held|--not-held]
validation_matrix references gaps                     # regenerates the "Gap list" section
validation_matrix references get ko2017_perf
validation_matrix references check                    # DOI format, verification consistency,
                                                      # declared-vs-found citations, orphans
validation_matrix references bibtex [--out refs.bib]  # for the paper
validation_matrix references where-used ko2017_perf   # matrix rows + code sites
```

What this buys over the current file: the DOI gap list becomes a command output instead of a
hand-edited section; a declared `Cited by the code: file:line` that no longer cites the work
becomes a `check` finding; and §14's reference map stops pointing at renumbered sections.

## 5. Row schema

One record per matrix row. A row is a coherent **group** of parametrisations (the current
"(3 cases)" granularity), which is also the granularity §13 uses. A row in turn carries **one
or more comparisons**: a single test routinely measures against several independent references
at once, each with its own tolerance and its own margin.

That is not hypothetical. Three rows of the current document already carry it, crammed into
slash-separated cells:

- §4.3 `test_composite_axial_tension` compares against **CCX 2.23, S8R** under
  `CCX_MEMBRANE_TOL = 0.015` *and* against an independent CLT bar under
  `CLT_ANALYTICAL_TOL = 0.02`, reporting two margins in one cell: `AE 41.67 um, CCX 41.94 um
  (0.65%); CLT bar 42.19 um (1.25%)`.
- §4.10 `test_outer_fibre_stress_matches_ccx_and_analytical` carries tolerance `15% / 20%` and
  margin `52.46 vs 51.64 MPa (1.58%); vs analytical 12.57%` — two tolerances, two margins, one
  row, one of them above the 5% suite rule and the other not.
- Every §3 node compares against the paper's N=16 cell *and* against the Kirchhoff closed form
  used for normalization, and §3 additionally reaches the `PAPER_REFS` 3D-exact values in a
  print that is never asserted.

A scalar `reference`/`tolerance`/`measured` triple cannot represent that without a
slash-concatenated string, which is unqueryable and uncheckable. Hence `comparisons`.

```yaml
- id: ko2017.twisted_beam.thin_inplane
  group: "3"                              # key into groups.yaml
  title: "test_3_5_twisted_beam_tables_12_to_13[0.0002667, In-plane]"
  tests:                                  # exact collected node ids, from collect-only
    - "tests/test_ko2017_performance.py::test_3_5_twisted_beam_tables_12_to_13[0.0002667-In-plane-1e-06-0.005256-0.001294-0.9978-0.01]"
  validates: "thin twisted beam, in-plane, N=16x96"
  result:                                 # our own measured value, when the test has one
    text: "AE 0.9972 normalized"
    raw: 0.9972
    unit: "w / w_Kirchhoff"
  flags: []                               # test-level judgements; near/gt5 are DERIVED
  notes: "hardcoded 0.92 replaced by the paper cell"
  history:
    - rev: e879eba
      note: "the first version recorded an 8.5% deviation that no longer exists"
      evidence: "pilot adjudication log 2026-10-02, entry ADJ-0007"
  comparisons:
    - label: "paper cell, N=16"           # which part of the test this comparison covers
      asserted: true
      reference:
        kind: paper                       # paper | code | analytical | self
        label: "Table 12, MITC4+ N=16: 0.9978"
        citation: ko2017_perf             # key into references.yaml when kind == paper
      tolerance:
        kind: rtol                        # rtol | atol | rel_err | sign | exact | subset
        value: 0.01
        source: "_TWISTED_BEAM_CASES field 8 (0.01 per case)"
        justified: true
        justification: null
      expected: "0.9978"
      measured:
        status: measured                  # measured | not_printed | not_measured | sign_only | n/a
        run: "34e2328"                    # revision the numbers were measured at
        date: 2026-10-02
        raw: 0.9972
        margin_pct: 0.06                  # the number §13 sorts and flags on
        text: "0.9972 (0.06%)"            # verbatim residual, for the paper table
    - label: "Kirchhoff closed form used for normalization"
      asserted: true
      reference: { kind: analytical, label: "w_ref = alpha p L^4 / D, alpha = 1.267e-3" }
      tolerance: { kind: rtol, value: 0.05, source: "assert_relative_error", justified: true }
      measured: { status: measured, raw: 1.0, margin_pct: 0.0, text: "1.0000 (0.00%)" }
    - label: "PAPER_REFS 3D-exact value, printed only"
      asserted: false                     # a printed comparison that no assert consumes
      reference: { kind: paper, label: "3D exact w_ref = 5.256e-3", citation: ko2017_perf }
      tolerance: { kind: rtol, value: 0.05, source: "never enforced" }
      measured: { status: n/a, margin_pct: null, text: "printed errors 523% to 6.2e8%" }
```

`asserted` and `measured.status` are orthogonal on purpose, and both are needed:

| case | `asserted` | `status` | example |
| --- | --- | --- | --- |
| asserted and printed | `true` | `measured` | the ordinary row |
| asserted, no residual printed | `true` | `not_printed` | the Markdown's "not printed" cells |
| printed, never asserted | `false` | `n/a` | §9.5; requires the `info_only` flag |

Field rules:

- `id` is stable and immutable once published. Scheme `<group-slug>.<file-slug>.<case-slug>`,
  ASCII, lowercase, dots as separators. Renaming requires a redirect entry so prose anchors and
  old links do not break.
- `tests` holds **exact collected node ids**, never prose. A row whose reference is out-of-band
  (no live test, e.g. §10.1) carries `tests: []` and `evidence: out_of_band`.
- `result` is **our** measured value, optional and only meaningful when the test produces one.
  It exists so that "what did AeroElast measure?" is a queryable field rather than something to
  parse out of a comparison's `text`. Each comparison then says what that same result was
  compared against, under which tolerance, with which margin. Where a test normalizes the same
  raw quantity against different bases (every §3 case), the shared `result` is the raw quantity
  and each comparison's `measured.text` carries its own normalized value.
- `comparisons` holds **at least one** entry. `label` names which part of the test the
  comparison covers ("modes 1-5", "the CCX leg", "the analytical leg", "N=16 cell"), because
  that is the only thing that tells two references of one test apart.
- `asserted: false` marks a comparison the test prints and never asserts. `check` requires the
  row to carry the `info_only` flag in that case, which is how §9.5 becomes derivable instead of
  hand-listed.
- `measured.status != measured` suppresses the numeric fields; `margin_pct: null` propagates to
  `slack: n/a`, exactly as today. An unasserted comparison uses `status: not_asserted`.
- `flags` stores only **human judgements** (`tautology`, `unfailable`, `overpromise`,
  `unjustified_tolerance`, `info_only`, `out_of_band`). `near` (slack < 1 percentage point) and
  `gt5` (tolerance > 5%) are **derived per comparison** and rejected if present in the file. A
  row is reported `near` or `gt5` when any of its comparisons is.
- Every flag id must exist in `flags.yaml`, which carries the evidence prose currently in §9.
  This turns §9 into a generated view of data instead of a hand-synced list.
- `expected` is per comparison: one gap can have a single-number reference side while another
  gap of the same test compares a whole mode table.

## 6. Group and flag registries

```yaml
# groups.yaml
- id: "3"
  title: '`tests/test_ko2017_performance.py` (Ko, Lee, Lee & Bathe 2017)'
  source_files: [tests/test_ko2017_performance.py]
  common_tolerance: { kind: rtol, value: 0.05, source: "assert_relative_error" }
  headline: "31 benchmark cases; worst 2.73% (Scordelis-Lo regular)"
  provenance_note: "measured at e879eba and reproduced digit for digit at 34e2328"
```

```yaml
# flags.yaml
- id: near
  derived: "slack_pp < 1.0"
  user_facing: true
- id: tautology
  section: "9.1"
  label: "Tautological references (the arithmetic under test re-implemented in the test)"
  legend: "May be cited as a limitation or a negative finding, never as a validation result."
- id: unfailable
  section: "9.2"
  # ...
```

## 7. Derived views (removed from hand maintenance)

| today | becomes |
| --- | --- |
| §13 consolidated matrix (252 rows) | generated from `rows/*.yaml`: one line **per comparison** (`tolerance`, `measured.text`, computed `slack_pp`, stored + derived flags), so the slash-concatenated cells disappear |
| §13.1 tolerances above 5% | generated: every row with `gt5`, requiring `justified: true` and a `justification` |
| §2.1 results at a glance | generated from `groups.yaml.headline` + computed worst margin per group |
| Index group table | generated from row counts and `len(tests)`, cross-checked against collect-only |
| §2 file-by-file inventory and the known-drift table | generated diff of `tests[]` against `pytest --collect-only -q` |
| §9 flag summary | generated from `flags.yaml` + rows carrying the flag |
| §14 test -> reference map | generated from `comparisons[].reference.kind == paper`, keyed by `citation`; one line per (row, citation) pair |
| `docs/references.md` entries, the gap list, "Cited by the code" | generated from `references.yaml`; `cited_by_found` computed from a code scan |

`slack_pp = tolerance_pct - margin_pct` when both are percentages; `near` iff `slack_pp < 1.0`.
For `kind: sign`, `exact`, `subset` and `status != measured`, `slack` is `n/a`, matching the
current table. Both flags are computed **per comparison**; a row is shown flagged when any of
its comparisons is, and §13 lists the comparison that triggered it.

## 8. CLI contract (`tools/validation_matrix.py`)

Read paths are the token-saving surface; write paths must be surgical and validated.

```bash
validation_matrix check  [--group ID] [--json]
validation_matrix find   SUBSTRING [--group ID] [--file F] [--limit N] [--json]
validation_matrix list   [--group ID] [--file F] [--id ID] [--flag FLAG] [--near] [--gt5] \
                         [--unit row|comparison] [--sort id|slack] [--fields A,B] [--json]
validation_matrix get    ID [--comparisons] [--json]
validation_matrix headline [--json]
validation_matrix set    ID PATH=VALUE... [--unset PATH]... [--add-flag F]... [--remove-flag F]...
validation_matrix regression [--group ID] [--scope S] [--write] [--json]

validation_matrix extract [--group ID] [--scope S] [--citation KEY] [--write] [--json]
validation_matrix references (list|gaps|get|check|bibtex|where-used|render)
validation_matrix gaps   [--json]
validation_matrix status [--json]
```

`regression` is the answer to "the values must serve as a regression, without re-reading the
document": it re-runs the group's scope, reads each node's printed residual with the pattern
the group declares in `groups.yaml`, pairs the Nth print with the Nth asserted comparison, and
diffs it against the stored `measured`. It prints only what is not `same` and exits 1 on drift,
so it is a gate. `--write` records the baseline (`status`, `raw`, `margin_pct`, `text`, `run`,
`date`), and is how a new group's rows gain their first measurements. A count mismatch between
prints and asserted comparisons is reported `unmapped`, never guessed: attaching a margin to the
wrong comparison would manufacture a baseline.

Read verbs, and what makes them worth preferring to the Markdown:

- `find` searches id, title, validates, group and node ids, and prints one line per row: id,
  section, title, and the test/comparison counts. It never dumps the document.
- `list --unit row` prints the derived summary: comparison count, **worst margin**, **tightest
  slack**, and the union of stored and derived flags. `--unit comparison` prints one line per
  comparison (label, tolerance, margin, slack, asserted, flags), which is what makes "my result
  against several references" queryable. `--near` (slack < 1 pp), `--gt5`, `--flag` and
  `--sort slack` narrow it; with `--unit comparison` they filter comparisons, not rows.
- `get ID` prints the stored record verbatim as YAML, so it can be edited and pasted back;
  `--json` is the machine form and `--comparisons` prints the derived table for that row.
- `headline` renders the derived §2.1 table: each group's headline, its row and comparison
  counts, the worst margin it measures, and how many comparisons sit above the 5% rule.

Write rules:

- `set` accepts dotted paths over an explicit schema (`SET_SCHEMA`), refuses unknown fields,
  refuses derived names (`near`, `gt5`, `slack_pp`), refuses `id`, and requires an index for a
  list (`comparisons[0].measured.margin_pct`, never `comparisons.measured...`).
- Values parse as bool, null, int, float, JSON (`[...]`/`{...}`) or a string; a value wrapped in
  matching quotes is taken literally.
- `set` is atomic and validated: it re-validates the edited row and **refuses to write** when the
  result would break an invariant, leaving the file untouched. It also refuses to write when the
  store already has errors, so one bad row cannot be extended.
- Writes are per-file (one YAML file), so concurrent group edits do not conflict. Note that YAML
  comments in a row file are not preserved across a `set`; keep comments in `groups.yaml`,
  `flags.yaml` and `prose/`.
- `render --check` asserts each committed Markdown file is byte-identical to a fresh render; CI
  runs it (T7).
- Output defaults to compact text (one line per unit); `--json` is for machine consumption.

## 9. `check` invariants (CONTRIBUTING rule 6, made executable)

**Matrix invariants.**

1. Schema-valid rows; unique `id`; no comparison with `reference.kind: schema`, because a row
   with no independent reference is a software test and belongs out of scope.
2. Every row has at least one comparison, and every comparison has a `label`, a `reference` and
   a `tolerance`. This is the invariant that keeps "my result against several references"
   expressible without a concatenated string.
3. Every collected node is either claimed by exactly one row, or belongs to a file declared
   `software_only: true` in `groups.yaml` with a reason. Nothing is silently unaccounted: the
   old `un_inventoried` debt bucket becomes an explicit classification, because the scheme
   covers physics tests only. **Until T8 lands the declaration is still hand-written**; after it,
   the path is the classification (`tests/software/**` never needs a row) and the bucket is
   deleted rather than maintained.
4. Every entry of `tests[]` exists in the collected node set.
5. **A collect-only run that returns 0 nodes is a failure, not a clean result.** This is not
   hypothetical: `python -m pytest -o addopts="" --collect-only -q tests/test_ko2017_performance.py`
   outside the conda environment prints `no tests collected` because the module-level
   `pytest.importorskip("_aeroelast")` skips the whole module. A silent 0 would read as
   "zero drift".
6. `measured.status == measured` implies `margin_pct` is a number; otherwise `margin_pct` is
   null. Evaluated per comparison.
7. Every `flags[]` id exists in `flags.yaml`; derived flags (`near`, `gt5`) are absent from the
   files and computed per comparison.
8. A comparison with `asserted: false` implies the row carries the `info_only` flag (§9.5), so
   the flag is derivable from the data instead of hand-listed.
9. Every comparison whose ratio tolerance exceeds 5% has `tolerance.justified: true` and a
   non-empty `justification`.
10. When a row carries `result`, its `raw` is numeric and its `unit` is a non-empty string.
11. Every prose anchor referenced from `groups.yaml` exists in the generated Markdown.

**Reference invariants.**

1. Every `comparisons[].reference.citation` used by any row resolves to a key in
   `references.yaml`.
2. `doi_status` is consistent with `held`, `doi` format and `verification_note`.
3. Every `cited_by_declared` entry resolves to an existing file and line, and
   `cited_by_found` covers it. A declared site that no longer cites the work is a finding.
4. No orphan keys: a reference neither cited by a comparison nor present in `cited_by_declared`
   must set `orphan_ok: true` with a reason (the bibliography legitimately holds sources the
   code does not cite yet).
5. `held_files` under `.sources/` are reported as `warn` when absent (gitignored), never as a
   failure.
6. `where-used` reports which (row, comparison label) pairs use a key, so a reference's blast
   radius is visible before it is edited.

## 10. Method: extract from the tests, cross-validate against the matrix

The migration does **not** trust the Markdown as input. It derives rows from the code and the
run, then uses the Markdown as an independent witness.

| stage | what it does | output |
| --- | --- | --- |
| **E1 Collect** | `pytest -o addopts="" --collect-only -q` in the pinned environment (`conda activate aeroelast-dev`, `_aeroelast` importable) | authoritative node set: `file::node[param...]` |
| **E2 Static extract (AST)** | per test function: decorators and parametrize tuples with their field names, `assert` expressions, `rtol`/`atol`/`tol` literals or resolved module constants, expected-value dicts/tuples, docstring and adjacent-comment citation text, `xfail` marks (including conditional and `strict=`), rule-3 promotion guards (`assert abs(1.0 - ratio) > 0.05, "promote it to an asserted row"`), and `print()` calls that no assertion consumes | candidate rows **without** `measured.*` |
| **E3 Param-id channel** | pytest already encodes some arguments in the node id. In §3 the id of the thin in-plane twisted-beam case is `...[0.0002667-In-plane-1e-06-0.005256-0.001294-0.9978-0.01]`, carrying the expected value *and* the tolerance; other ids hide them (`expected_table2_30`, `expected0`, `expected_mitc40`). Where present, the id is a free second witness for E2 | per-node agreement/divergence between id and AST |
| **E4 Normalize the matrix** | parse §3-§8 tables into the same key space by resolving each prose row title to node ids | md-projected rows |
| **E5 Diff** | field-by-field comparison: `match`, `only_in_code`, `only_in_md`, `value_conflict` | machine-readable adjudication report |
| **E6 Adjudicate** | for each conflict, **re-run the exact nodes** (`-s`, and `-rA`) and decide from the run, the test source and the published reference. Every verdict is recorded with its evidence | adjudication log entries |
| **E7 Freeze** | the new rows become the authority for that group; the Markdown section is regenerated and its diff reviewed | rendered section + `check` green |

**Authority rule, stated once.** The test code is the authority for *the tolerance as enforced*
and for *the reference as asserted*. The published paper/report is the authority for *the
reference as published*. The run is the authority for *the measured margin*. Where those
disagree, the disagreement is a **finding** to record, never a value to pick silently — which
is exactly the class of defect §3 currently hides (§12).

## 11. Adjudication record

The log is validated data, not a document: `check` reads
`docs/validation/adjudications/<group>.yaml` and fails on an entry that names a row that no
longer exists, states no verdict, resolves no conflict, or carries a malformed date. An audit
trail nobody reads is exactly what this migration exists to remove, so the file that records
the adjudications is held to the same standard as the store it describes.

A group with no conflict still gets a file: it records the reconciliation and the findings the
probe produced. Section 3's is `docs/validation/adjudications/3-ko2017.yaml`, and it carries the
slot identity (`23 Markdown rows covering 37 case slots` against `31 nodes and 37 comparisons`),
the two phrasings the document uses for a row's case count, and the `tol=` to `rtol` mapping.

Every `value_conflict`, `only_in_code` and `only_in_md` produces one entry, stored next to the
group it belongs to and referenced from the affected row's `history`:

```yaml
# docs/validation/adjudications/3-ko2017.yaml
- id: ADJ-0007
  group: "3"
  row: ko2017.twisted_beam.thin_inplane
  field: reference.label
  candidate_code: "3D exact w_ref = 5.2560e-3 from PAPER_REFS, normalized against Kirchhoff"
  candidate_md:   "Table 12, MITC4+ N=16: 0.9978"
  check: "assert np.isclose(norm, expected_normalized, rtol=0.05) at tests/test_ko2017_performance.py:1488"
  rerun: "pytest -o addopts=\"\" -q -s 'tests/test_ko2017_performance.py::test_3_5_twisted_beam_tables_12_to_13[0.0002667-In-plane-1e-06-0.005256-0.001294-0.9978-0.01]'"
  rerun_result: "Norm vs Kirchhoff: 0.9972 (expected: 0.9978, error: 0.06%)"
  verdict: "md_was_correct_but_incomplete"
  resolution: "reference.label kept; add reference.also_asserts to record the unasserted 3D print"
  decided_by: maintainer
  date: 2026-10-02
```

The log is the audit trail for the migration and stays in the repository: it is the evidence
that the new schema was not just a re-typing of the old prose.

## 12. Pilot: §3, `tests/test_ko2017_performance.py`

Chosen because it concentrates the hardest cases and the widest mismatch. Facts measured now:

- 2132 lines, **31 collected nodes**, **23 Markdown rows**, 9 test functions (`test_3_1` ..
  `test_3_9`), 13 `parametrize` decorators (several stacked, so one function yields 6 nodes).
- Tolerances are **three different mechanisms** in one file: `rtol=0.05` inside
  `assert np.isclose(...)` (`:821`, `:845`), `tol=0.05` via the local `assert_relative_error`
  helper (`:105`, called at `:940`, `:1129`, `:1287`, `:1920`, `:2023`, `:2130`), `tol=0.01`
  per case in the `_TWISTED_BEAM_CASES` loop (`:1363`, `:1393`) and a literal
  `rel_err < 0.03` for the Raasch hook (`:1698`).
- Conflicts already visible from static inspection, before any run — these are the pilot's
  real payload:
  - **Two reference semantics in one file.** `PAPER_REFS` (`:56`) holds **3D exact** values
    (`clamped_square_plate = -2.2137e-1`, `twisted_beam_thin_inplane = 5.256e-3`, ...), while
    the Markdown rows cite **normalized N=16 paper cells** (`0.9984 / 0.9980 / 0.9979`,
    `0.9978`). Both reach the assertion indirectly. The extractor must record both and the
    schema must say which one a row cites.
  - **The informational-only print.** `:671-676` prints `[x] Norm vs Paper 3D` by dividing a
    displacement by `PAPER_REFS[...] * pressure` and comparing against a *normalized* value,
    printing errors from 523% to 6.2e8% that no test fails on. §9.5 flags this. The extractor
    can emit it as a candidate `info_only` flag **with evidence**, not as a verdict.
  - **A resolved xfail window.** `_TWISTED_BEAM_CASES` carries a per-case `xfail` field and a
    comment explaining a window that "existed to accommodate an 8.5% deviation". The Markdown
    records this as RESOLVED. The extractor asserts `xfail is None` for all four cases.
  - **An unreachable code path.** `_assemble_global` selects triangles from the node count, but
    every live call site passes `use_triangular=False`; the Markdown notes it as out-of-band.
- The 23-vs-31 gap is the row-granularity rule to settle in the pilot: one row per
  function-and-case-class, with `tests[]` listing every node it claims.
- The multi-reference rule is not hypothetical in this group: every §3 node measures against the
  paper's N=16 cell **and** against the Kirchhoff closed form used for normalization, and three
  of the nine functions additionally reach the `PAPER_REFS` 3D-exact values in a print that is
  never asserted. So the pilot is where `comparisons` proves it can hold one own result against
  several references without a concatenated string.

Pilot acceptance criteria:

1. 31/31 nodes claimed by exactly one row; each row's `tests[]` resolves.
2. `only_in_code` = 0 and `only_in_md` = 0 after E5, or each is listed in the adjudication log.
3. Every `value_conflict` has an `ADJ-*` entry with a re-run command and its result.
4. The regenerated §3 renders byte-identical to the current section **except** for the
   adjudicated corrections, which are listed explicitly in the log.
5. `validation_matrix check --group 3` exits 0, and exits 1 when a synthetic node is removed.
6. `references.yaml` holds the §3 keys (`ko2017_perf`, `ko2017_new`, `ko2017_nonlinear`,
   `knight1997`, `dvorkin1984`, `lee2014`, ...), `references get ko2017_perf` shows
   `doi: 10.1016/j.compstruc.2017.08.003` with `doi_status: verified`, and
   `references check` reports the declared-vs-found citation mismatch for
   `src/aeroelast/core/mesh/generators.py:65`.

Pilot run environment (required; without it collect-only silently returns 0):

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev
cd ~/Desktop/dev/fem-shell
python -m pytest -o addopts="" --collect-only -q tests/test_ko2017_performance.py | grep -c '::'   # 31
```

## 13. Migration plan

### 13.1 Pilot task list (the minimal production test)

The pilot is the first production slice, so its deliverables are the real artifacts. Two
consequences:

- The extractor and the matrix normalizer are **CLI subcommands** of `tools/validation_matrix.py`
  (`extract`, `diff-against-md`), not a scratch script. WU-3..WU-5 run them on the remaining
  groups without rewriting them.
- The pilot writes `rows/3-ko2017.yaml` as the finished §3 store, and the regenerated §3 of
  `docs/validation-matrix.md` is the proof that rendering is lossless.

| # | task | deliverable | acceptance criteria |
| --- | --- | --- | --- |
| T1 | Schemas and registries (WU-0) | `docs/validation/schemas/validation-row.schema.json`, `schemas/reference.schema.json`, `docs/validation/groups.yaml` (seeded with §3), `docs/validation/flags.yaml` (the §9.1-§9.5 registry plus the derived `near`/`gt5`) | a 3-row fixture validates; `check` exits 1 on a row with a stored derived flag, an unknown flag id, or `measured.status == measured` without `margin_pct` |
| T2 | CLI read side (WU-1) | `validation_matrix.py` with `find`, `get`, `list`, `headline`, `--json`, schema loading | each query over the fixture returns the expected records; unknown fields in `set` are refused |
| T3a | Bibliography store and verbs | `docs/validation/references.yaml` seeded with the section 1 entries; `references list\|gaps\|get\|check\|bibtex\|where-used`; `code_mentions` in the schema; hand-rolled validation of `doi_status`, held state and citation sites | `references get ko2017_perf` shows `doi: 10.1016/j.compstruc.2017.08.003` with `doi_status: verified`; a declared site that no longer contains a `code_mentions` string is reported and a live one is not; a row citing an unresolvable key fails `check` |
| T3b | Full bibliography migration and the generator | every entry of the current `docs/references.md` in the store, with its `cited_by_declared` sites and their `code_mentions`; `references render`; `docs/references.md` regenerated | the generated `references.md` diff against the current file is reviewed line by line for content loss; the gate in `tests/test_mitc4plusd_traceability.py` still resolves |
| T4 | Code extractor (WU-P1) | `validation_matrix extract --scope tests/test_ko2017_performance.py`: collect-only node set + AST rows + param-id cross-check, no `measured.*` | 31 nodes, every one claimed, 0 unexplained extraction misses; byte-identical output across two runs |
| T5 | Matrix normalizer and diff (WU-P2) | `validation_matrix diff-against-md --group 3`: §3 tables projected into the same key space, compared field by field | every conflict classified as `match`, `only_in_code`, `only_in_md` or `value_conflict`; the 23-vs-31 gap explained per row |
| T6 | Adjudication and freeze (WU-P3) | re-run each conflicting node (`-s`, `-rA`), `docs/validation/adjudications/3-ko2017.yaml`, complete `rows/3-ko2017.yaml` with `measured.*` | the six acceptance criteria in §12; no conflict resolved without a recorded command and its output |
| T7 | Elimination, classification and CI | repoint every live reference at the store (`docs/formulations/{mitc4plus-2017-extract,shell-elements,materials,solvers}.md`, `docs/validation-environment.md`, `tests/test_shell_convergence.py`, `openspec/specs/mitc4plusd-element/spec.md`, `tests/test_mitc4plusd_traceability.py`), delete `docs/validation-matrix.md` and `docs/references.md`, rewrite CONTRIBUTING rule 6, write `docs/validation-policy.md` with the surviving prose, and move the physics-only classification into `groups.yaml` | no live document points at a deleted path; archived `openspec/changes/archive/**` left untouched as history; `check` and `regression` exit 0 in CI; a synthetic unclassified node exits 1 |

### 13.2 Progress

| task | status | commit | evidence |
| --- | --- | --- | --- |
| design | done | `7a011f3` | this document |
| T1 | done | `d947535`, `7b2bcee` | `python -m pytest tools/tests` -> 22 passed; `check` -> 0 errors; ruff clean. Reviewed as `review-4ecd3ccb79863462` (tier low, no lenses, `non_executable_only`) |
| T2 | done | - | `python -m pytest tools/tests` -> 43 passed (21 new); `ruff check` 0.16.0 clean; `check` -> 0 errors; `headline` renders the group 3 row |
| T3a | done | - | `python -m pytest tools/tests` -> 62 passed (19 new); `ruff check tools/` clean; matrix `check` 0 errors; `references check` -> 6 entr(ies), 0 errors, 1 warning |
| T3b | pending | - | - |
| T4 | pending | - | - |
| T5 | done | - | `diff-against-md --group 3` -> `23 Markdown row(s) covering 37 case slot(s) against 31 collected node(s); 0 function(s) need adjudication`, exit 0; 80 tool tests. **The command was retired** once both views were deleted (see T7): its only input is gone, so its tests, its helpers and its CLI entry went with it |
| T6 | done | - | `check` validates `docs/validation/adjudications/3-ko2017.yaml`; no conflict needed adjudication; 90 tool tests |
| T7 | pending | - | - |
| T8 | done | `93c7384` (T8a), below (T8b) | T8a: `tests/` is a package, helpers in `tests/support/`, 40 import sites rewritten, 527 nodes collected with zero import errors. T8b: 44 validation files in six domains and 9 software files in three, the drift bucket deleted, and the path made the classification |

| T7 detail | deliverable | acceptance criteria |
| --- | --- | --- |
| Eliminate the generated views | (1) migrate the not-citable list into `docs/validation/gaps.yaml` — **done at T7a**, because it is the one part of the document that is not derivable: a gap is a claim about absence, and a row that does not exist cannot be queried; (2) write `docs/validation-policy.md` with the citation policy, the validity envelope's criteria, and pointers to `gaps`, `check` and the diagnostic index; (3) repoint the live references (`docs/formulations/{mitc4plus-2017-extract,shell-elements,materials,solvers}.md`, `docs/validation-environment.md`, `tests/test_shell_convergence.py`, `openspec/specs/mitc4plusd-element/spec.md`) and rewrite CONTRIBUTING rule 6; (4) repoint `tests/test_mitc4plusd_traceability.py:36` at the store; (5) delete `docs/validation-matrix.md` and `docs/references.md`; (6) a Make target running `check` and `regression` | no live document points at a deleted path; archived `openspec/changes/archive/**` untouched; `check`, `regression` and `gaps` pass; `diff-against-md` reports a clear error once its input is gone instead of crashing — and it is a migration-time command, not a CI one, because CI cannot diff against a deleted file. **As built, T7 went further and removed it**: a command whose only input no longer exists has no reason to stay in the CLI, where its presence invites a rerun that can only fail. Its history is in git, and what it proved is recorded in `docs/validation/adjudications/3-ko2017.yaml` |

| T8 detail | deliverable | acceptance criteria |
| --- | --- | --- |
| Reorganize `tests/` by role then domain | `git mv` the 56 test files into `tests/validation/<domain>/` and `tests/software/<domain>/`; promote `conftest.py`, `_ccx_io.py` and `_openfast_bem.py` into an importable `tests/support/`; rewrite the store's node ids and `source_files` by path prefix so the 37 baselines survive without a re-measure; update the 142 `tests/…` references outside `tests/`; rewrite the authoring rules (`CONTRIBUTING.md` test rules and `docs/adding-validation-tests.md`) with where a new test goes and what its location commits it to; delete the `un_inventoried` bucket | `pytest --collect-only` reports the same node count as before the move; `check`, `regression` and `diff-against-md` pass; a cross-directory import still resolves; no hand-written classification remains in `groups.yaml` |

Branch: `feat/validation-matrix-store`.

### 13.3 Review candidates must fit the provider budget

The review provider refused the accumulated candidate at preflight with
`lens_context_budget_exceeded` (13 paths / 4436 lines, base `ac69c8d4`): "the candidate's complete
reviewer evidence exceeds the native context budget... immutable candidate evidence is never
truncated and retrying this exact candidate cannot succeed". `mutation_outcome: not_started`, so no
authority or lineage was created and nothing needs repairing. The maintainer chose to leave the
candidate unreviewed and continue, so no further transaction is started for it.

Consequences for this work:

- **A review candidate is one work unit, not the branch.** The pilot's own commits are already
  separable: `7a011f3` (design), `d947535` (T1), `100b5a8` (T2). Reviewing them as a chain needs a
  distinct base ref per candidate.
- **The provider's base ref may include work that is not ours.** `ac69c8d4` predates four
  unrelated maintainer commits, so its window swept in the blade tests. A base ref that isolates
  this branch is `d322f73` (the merge-base with `main`).
- **A docs-only candidate is comfortably admitted**: the first review closed immediately as tier
  `low`, no lenses, `non_executable_only`. A design document of ~650 lines is not the problem;
  the tool plus its tests plus the accumulated commits are.
- **The facade START contract**, learned by three attempts: `mode` is mandatory, a committed range
  needs `baseRef` paired with `committedOnly: true` in the same JSON string, and the workspace
  projection sees nothing once everything is committed (`empty_candidate_base_ref_required`).

| WU | deliverable | acceptance evidence |
| --- | --- | --- |
| WU-0 | `docs/validation/schemas/validation-row.schema.json`, `schemas/reference.schema.json`, `groups.yaml`, `flags.yaml`, empty `rows/` | `check` passes on an empty store |
| WU-1 | `tools/validation_matrix.py` read side (`find`, `get`, `list`, `headline`) over a 3-row fixture | fixture queries return the expected records |
| WU-2 | `references.yaml` seeded for the references the pilot group uses; `references.md` generator | `references check` findings match a hand audit of `docs/references.md` |
| WU-P1 | E1-E3 over §3: node set + AST rows + param-id cross-check, no `measured.*` | 31 nodes, 0 unexplained extraction misses |
| WU-P2 | E4-E5: normalize the §3 Markdown rows and diff | adjudication report with every conflict classified |
| WU-P3 | E6-E7: adjudicate (re-runs), freeze rows, regenerate §3 | the six pilot acceptance criteria above |
| WU-3 | §4 migration using the E1-E7 pipeline the pilot proved | same criteria for 87 rows |
| WU-4 | §5 + §6 migration (195 rows) | same |
| WU-5 | §7 + §8 migration (139 rows) + §10 out-of-band | same |
| WU-6 | prose moved to `prose/`, both Markdown files fully generated | `render --check` clean; diffs reviewed line by line for content loss |
| WU-7 | retire the known-drift table: populate `un_inventoried:` and close the 7-file gap | `check` green with every collected node claimed or explicitly excluded |
| WU-8 | CI/Make target `validation-matrix-check` | fails on synthetic drift, passes on the tree |

Each work unit is a separate, reviewable commit on a feature branch; the generated Markdown is
regenerated in the same commit as the data it renders.

## 14. Risks

| risk | mitigation |
| --- | --- |
| The extractor encodes a wrong guess about intent | E5 keeps the Markdown as an independent witness; E3 gives a second in-file witness; E6 requires a re-run per conflict |
| Adjudication turns into a judgement call that silently widens a tolerance | the authority rule (§10) plus the adjudication log (§11), which records the command and its output |
| ~240 rows migrated by hand drift from the prose | §3 is the pilot; the pipeline, not hand editing, carries the remaining groups |
| Row `id`s churn and break anchors | immutable scheme (§5), redirect map, `check` matrix invariant 11 |
| The store duplicates the generated Markdown instead of replacing it | `render --check` in CI (WU-8); both Markdown files declare themselves generated |
| Reviewers stop reading a 196 KB file and read only YAML | the generated Markdown remains the paper-facing artifact; the YAML is the authoring artifact |
| `margin_pct` still hand-typed and therefore stale | accepted in this phase; §15 addresses it |
| Tolerating absence of `.sources/` hides a missing PDF | `warn`, not `fail`, and `references gaps` lists every unverified DOI explicitly |
| A store adds a dependency to a suite that must run anywhere | PyYAML is already a runtime dependency (`pyproject.toml`); `jsonschema` is optional and `check` degrades to hand-rolled validation without it |

## 15. Phase 3 (out of scope here): measured margins from the run

The margins are currently transcribed from `pytest -o addopts="" -q -s` output by reading
free-form `print()` calls; the suite has no shared residual emitter, so this cannot be automated
today. The follow-on change is:

1. A shared helper (`tests/_matrix.py`) that emits one JSON object per asserted residual on a
   machine-readable channel, e.g. `PYTEST_MATRIX_RESIDUALS=<path>`.
2. `validation_matrix sync <residuals.jsonl>` that updates `measured.*` for every row whose
   `tests[]` intersection produced a residual, and reports rows with no residual.
3. `measured.run` set to the tree revision of the run, so a stale margin is detectable, not just
   suspicious.

This is a `tests/` change with its own review and its own work units; it is listed here only so
that the schema in §5 does not have to be redesigned later. The same channel would let
`references check` verify `cited_by_declared` against AST-found citations rather than grep.

## 16. Open decisions for the maintainer

Decision 6 is settled: the references migration ships with the pilot. Decisions 1-5 are settled
by default inside the pilot, and the pilot's output is what they are judged against:

1. **Row granularity** — **settled at T5, and it is not the question it looked like.** The store's
   row is the *collected test*, and the Markdown's row is a *comparison slot*: §3 has 23 Markdown
   rows covering **37 case slots**, and the store has 31 nodes with **37 comparisons**. The two
   counts differ because a test that asserts two benchmark tables is one node and two slots
   (`test_3_1` is 4 Markdown rows over 6 nodes; `test_3_2` is 2 over 6; `test_3_9` is 2 over 4;
   the other six functions are 1 node per slot). `diff-against-md` computes that identity and
   fails when it does not hold, so the reconciliation is checked, not asserted.
   Mapping Markdown row to node *per row* is deliberately not attempted: the distinguishing
   parameter is not always in the assertion (`test_3_1`'s four rows split on `distorted`, which
   appears in the parametrisation and not in the assert), so a per-row mapping would be a guess.
   Function-level comparison is what the evidence supports, and it is what the command reports.
2. **Generated output** — **settled: eliminated.** Both Markdown views are deleted once the
   store covers them, so the anchor-redirect question disappears with them. `references render`
   survives as an export/debug command writing to a throwaway path, never as a committed
   artifact, and CI runs `check` rather than `render --check`.
3. **Prose ownership** — default: `prose/*.md` is maintainer-edited directly and only tables are
   generated.
4. **`jsonschema` dependency** — default: use it if it is already importable in the pinned
   environment, otherwise hand-rolled validation in `check` with no new dependency.
5. **Generated Markdown under version control** — default: committed, because existing
   cross-references and the `tests/test_mitc4plusd_traceability.py` gate read the file, and
   `render --check` keeps it honest.
6. **Scope of the references migration** — **settled: in the pilot** (T3).

## 17. Non-goals

- Not changing any test, tolerance, reference or measured number. This change moves data and
  adjudicates conflicts; it does not judge the physics.
- Not rewriting archived `openspec/changes/archive/**`: those documents are the record of what
  was true when they were written, and repointing them at paths that did not exist yet would
  falsify history.
- Not resolving the §9 flags. Migrating a flag into `flags.yaml` records it; it does not clear it.
- ~~Not deleting either Markdown file~~ — **reversed at T7**. Both were deleted once the store
  could answer everything asked of them: most of their bulk was data, and data belongs in the
  store. What must be *read* fell from 193 KB to a 5.1 KB policy.
- ~~Not migrating the `tests/test_mitc4plusd_traceability.py` gate~~ — **done at T7**, and the
  gate moved to `tests/software/contracts/` at T8b because it validates a contract between
  documents, not a physical quantity. It reads `docs/validation/references.yaml` directly.
- Not filling the DOI gaps. `references.yaml` makes them listable; obtaining the DOIs is
  bibliography work.

## 18. T11: the exclusion rule (proposal, awaiting one decision)

### 18.1 What is measured

Run over the 42 validation files that still have no group, against a throwaway copy of the store
so the real one is untouched:

| | count |
| --- | --- |
| collected tests | 432 |
| rows the extractor derives | 229 |
| tests with no row | 203 |

Every one of the 203 is `no_comparison`: no comparison site in the test. None is `not_in_ast` (the
class-based node-id fix closed that) and none is `tolerance_unresolved` (an unreadable comparison
now costs its comparison, not its row).

They cluster into 54 blocks of file + class, and the reasons repeat across those blocks: matrix
structural properties (`TestStiffnessMatrix::test_symmetry`), plumbing (`TestCoordinateTransforms`,
`TestOmegaProviders`, `TestUseRustFlag`), configuration loading (`TestBladeAeroYAML`), sanity
(`TestCompositeSanity`). A small vocabulary covers all of them.

The split that matters:

| | count |
| --- | --- |
| compare against another code (Rust, CCX, OpenFAST) | 50 |
| the rest | 153 |

### 18.2 The rule is already in the policy

`docs/validation-policy.md` rule 6: a reference must be independent -- a different code, a
published cell, or a closed form; a formula re-implemented inside the test is not one. Read
backwards, a test that compares one of our models with itself (`K` against `K.T`) or with nothing
(YAML loading, an API shape, a flag's plumbing) has no independent reference and no physical
quantity, so it is not a row. That is not a defect: it is a test doing its job in a suite that is
not only a validation suite.

It also means the 50 are **not** excludable on the same grounds. A Rust-against-Python parity test
compares against a different code, which rule 6 admits, so by the rule it *should* be a row -- and
the extractor cannot read it, because the assertion is `allclose(a, b)` rather than a residual
against a tolerance.

### 18.3 The mechanism

The declaration belongs in `gaps.yaml`, which already holds claims about absence, and `check`
enforces closure:

1. For a group whose scope is a test file, every collected node in that file is either claimed by
   a row or covered by a `not_a_validation` entry naming it.
2. A node that is neither is an error. A new test with no comparison site therefore cannot join a
   grouped file silently.

That is the drift bucket inverted: the classification is declared, validated, and impossible to
leave implicit. An entry needs `tests` (node ids or a class prefix), `status: not_a_validation`,
and `reason` drawn from the shared vocabulary.

### 18.4 The decision

The 50 comparisons against another code need an answer that rule 6 does not give:

- **Hand-write them** as out-of-band rows, the way group 10 works. Honest, and 50 rows of one-off
  work.
- **Teach the extractor to read them**: an `allclose(a, b, ...)` against a second code is a
  comparison with a known shape, and a residual printed next to it is already handled. More work
  once, and it covers every future parity test.
- **Exclude them** as software. Contradicts rule 6, which is the maintainer's own text.


### 18.5 A parity file mixes kinds, so the kind is decided per comparison

Reading the comparison sites of one clean parity file -- `tests/validation/parity/test_orthotropic_shell_parity.py`,
three tests, one assertion each, nothing unclaimed -- gives three different references:

| line | the comparison | honest kind |
| --- | --- | --- |
| 342 | the AeroElast displacement against `ccx_disp` | `code` |
| 404 | the AeroElast displacement against the Euler-Bernoulli value the test prints | `analytical` |
| 485 | the composite stiffness matrix against the isotropic one | `self` |

The module titles itself "AeroElast vs CalculiX", and declaring that would have been false for two of
its three comparisons. `reference_kind` is a group default and is right only where every comparison
shares one kind -- group 3, where all of them are cells of the paper. For a mixed file the group
declares none, `check` fails on the null kinds, and each comparison is decided. The workflow is
extract, declare, check, commit, and never the intermediate state.

So T11 cannot be filled mechanically, and the decision is not 226 rows one at a time either. It is
one criterion applied per comparison: **the reference the comparison names decides the kind**, and
the label the extractor derived is already the evidence for it. `rel_error` against `ccx_uy` is a code
comparison; the same `rel_error` against a value the test prints as "Analytical" is a closed form;
`||K_comp - K_iso|| / ||K_iso||` compares our own two models and is `self`. Approving the criterion
is what unblocks the domain, because the criterion is what makes 47 declarations written once.

