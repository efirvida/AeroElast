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
  groups.yaml                 # group registry: id, §, title, source files, headline
  flags.yaml                  # declarative flag registry (the §9 content, as data)
  rows/
    3-ko2017.yaml             # ~31 rows
    4-ccx-parity.yaml         # ~87
    5-analytical.yaml         # ~40
    6-element-invariants.yaml # ~155
    7-rotor-fsi.yaml          # ~67
    8-bem-aero-mesh.yaml      # ~72
    10-out-of-band.yaml       # §10.1-§10.3, evidence with no live test
  prose/
    validation-matrix/
      guide.md  validity-envelope.md  01-provenance.md  02-suite-snapshot.md
      09-flag-summary.md  11-cross-references.md  12-reproducing.md
      13-1-tolerances-above-5.md  14-reference-map.md
    references/
      preamble.md             # the recovered-PDF story and the verification convention
schemas/
  validation-row.schema.json
  reference.schema.json
tools/validation_matrix.py
docs/validation-matrix.md     # generated, deterministic
docs/references.md            # generated, deterministic
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
  notes: "the README entry omitted Bathe"
```

Field rules:

- `key` is the join key: a matrix row's `reference.citation` holds it, and the §14 map is
  generated from it. `check` fails if any row cites an unresolvable key.
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
        kind: paper                       # paper | code | analytical | self | schema
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
validation_matrix find <substring> [--group G] [--file F]     # ids + one-line summary
validation_matrix get <id> [--json|--yaml|--md]               # one record
validation_matrix list [--group G] [--file F] [--flagged|--flag X] \
                       [--near] [--gt5] [--sort slack] [--fields a,b,c]
validation_matrix headline                                    # derived §2.1 table
validation_matrix check [--collect-only-out FILE]             # schema + reconciliation, exit 1 on drift
validation_matrix render [--out docs/validation-matrix.md] [--check]
validation_matrix references <list|gaps|get|check|bibtex|where-used>
```

```bash
# surgical writes: exactly one record, validated before and after
validation_matrix set <id> measured.margin_pct=0.06 measured.run=HEAD
validation_matrix set <id> tolerance.value=0.02 tolerance.justification="CCX scatter"
validation_matrix set <id> --add-flag unjustified_tolerance
validation_matrix add-row --test <node-id>                    # seed a row from the code, see §10
validation_matrix rm <id>
```

Write rules:

- `set` accepts dotted paths, refuses unknown fields, refuses derived fields (`flags.near`,
  `slack`), and re-validates the whole row against the schema before writing.
- Writes are per-file (one YAML file), so concurrent group edits do not conflict.
- `render --check` asserts each committed Markdown file is byte-identical to a fresh render;
  CI runs it.
- Output defaults to compact text (one line per row); `--json` is for machine consumption.

## 9. `check` invariants (CONTRIBUTING rule 6, made executable)

**Matrix invariants.**

1. Schema-valid rows; unique `id`.
2. Every row has at least one comparison, and every comparison has a `label`, a `reference` and
   a `tolerance`. This is the invariant that keeps "my result against several references"
   expressible without a concatenated string.
3. Every collected node is claimed by exactly one row, or its file is listed in an explicit
   `un_inventoried:` allowlist in `groups.yaml` with a reason. This retires the hand-written
   known-drift table.
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
| T1 | Schemas and registries (WU-0) | `schemas/validation-row.schema.json`, `schemas/reference.schema.json`, `docs/validation/groups.yaml` (seeded with §3), `docs/validation/flags.yaml` (the §9.1-§9.5 registry plus the derived `near`/`gt5`) | a 3-row fixture validates; `check` exits 1 on a row with a stored derived flag, an unknown flag id, or `measured.status == measured` without `margin_pct` |
| T2 | CLI read side (WU-1) | `validation_matrix.py` with `find`, `get`, `list`, `headline`, `--json`, schema loading | each query over the fixture returns the expected records; unknown fields in `set` are refused |
| T3 | References store and generator (WU-2) | `docs/validation/references.yaml` seeded with the keys §3 uses; `references list\|gaps\|get\|check\|bibtex\|where-used`; `docs/references.md` generated | `references get ko2017_perf` shows `doi: 10.1016/j.compstruc.2017.08.003`; `references check` reports the declared-vs-found mismatch at `src/aeroelast/core/mesh/generators.py:65`; the generated `references.md` diff against the current file is reviewed for content loss |
| T4 | Code extractor (WU-P1) | `validation_matrix extract --scope tests/test_ko2017_performance.py`: collect-only node set + AST rows + param-id cross-check, no `measured.*` | 31 nodes, every one claimed, 0 unexplained extraction misses; byte-identical output across two runs |
| T5 | Matrix normalizer and diff (WU-P2) | `validation_matrix diff-against-md --group 3`: §3 tables projected into the same key space, compared field by field | every conflict classified as `match`, `only_in_code`, `only_in_md` or `value_conflict`; the 23-vs-31 gap explained per row |
| T6 | Adjudication and freeze (WU-P3) | re-run each conflicting node (`-s`, `-rA`), `docs/validation/adjudications/3-ko2017.yaml`, complete `rows/3-ko2017.yaml` with `measured.*` | the six acceptance criteria in §12; no conflict resolved without a recorded command and its output |
| T7 | Render, check and CI wiring (WU-8 for the pilot) | `render --check` for the pilot group, `check --group 3`, a Make target | exits 0 on the tree; exits 1 on a synthetic drift (a removed node, a tampered generated file), and exits 1 on a 0-node collect-only |

| WU | deliverable | acceptance evidence |
| --- | --- | --- |
| WU-0 | `schemas/validation-row.schema.json`, `schemas/reference.schema.json`, `groups.yaml`, `flags.yaml`, empty `rows/` | `check` passes on an empty store |
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

1. **Row granularity** — default: one row per function-and-case-class, with `tests[]` listing
   every node the row claims and `comparisons[]` listing every independent reference it measures
   against. The Markdown's 23 rows for 31 §3 nodes implies the row rule; the slash-concatenated
   tolerance cells (§4.3, §4.10) imply the comparison rule. Revisit at T5 if the 23-vs-31
   mapping is not expressible as a stable rule.
2. **Generated output** — default: keep single `docs/validation-matrix.md` and
   `docs/references.md`, because that preserves every existing anchor without a redirect map.
   Splitting is a later change if the single-file render diff becomes unreviewable.
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
- Not resolving the §9 flags. Migrating a flag into `flags.yaml` records it; it does not clear it.
- Not deleting either Markdown file; both remain the paper-facing views.
- Not migrating the `tests/test_mitc4plusd_traceability.py` gate to `references.yaml`. It can
  resolve against the new keys later; changing that test is its own work unit, and until then
  `references.md` is still generated, so the gate keeps passing.
- Not filling the DOI gaps. `references.yaml` makes them listable; obtaining the DOIs is
  bibliography work.
