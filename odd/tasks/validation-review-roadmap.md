# [validation] general review roadmap: what to fix, review and build before the solver is citable

*Local mirror of issue #33, which replaced #18 on 2026-10-10. The issue is the order of record; keep this file in step with it.*

## What this issue is

The order of record for making AeroElast citable as an FSI framework for wind turbines with a
**shell structural solver** coupled through **preCICE** to any fluid solver. The production fluid
is meant to be OpenFOAM; the built-in CCBlade BEM is a cheap fluid for development and coupling
checks.

#18 tracked defects one by one. Its closures were real work, but "closed" ended up meaning "we
have an explanation", and an issue could be closed by one session and reopened by the next. The
clearest case is #20. Its closure rested on a ±30 % band, five times the effect it named, and
that band could not see that the direction was inverted (`72f64de`). This issue replaces
narrative closure with **measured verdicts that expire on their own**.

## The acceptance rule

- A **physics comparison** passes at **≤ 5 % relative error** against an independent reference
  (a paper, a closed form, or another code), at a converged mesh or with a convergence trend, on
  a global scalar where one exists.
- A wider bound is an **exception**. It is valid only if:
  - an independent source puts a number on the model difference;
  - that number was fixed **before** measuring;
  - the band is smaller than the effect it claims to detect.

  "The models differ" with no number attached is not an exception.
- **Arbiters.** Data to validate against is limited, so the arbiter depends on what exists. It
  can be:
  - experimental data;
  - a closed form or theory;
  - published numerical results;
  - numerical results we produce with another code, such as CalculiX or OpenFAST.

  Experimental data or a closed form settles a comparison alone. **Numerical references settle
  it only in pairs:** agreement with a single model is not validation, and at least two
  independent ones are needed. With no arbiter available, the quantity is **not validated**,
  and the record says so. Experimental validation is a kind of arbiter, available or not on any
  level. It is not a level of its own.
- **Identity** checks (round-off or consistency bounds of our own computation) and
  **regression pins** (no independent reference) are outside the rule. They count toward no
  level.
- The maintainer decides whether to close or reopen. Agents measure and propose.

## The instrument

The store is the cheap way to see where validation stands. Since `a20a78b` it carries a contract
layer that `extract` and `coherence` never write:

- `docs/validation/levels.yaml`
  - places every group on one level of the ladder;
  - sets the 5 % ceiling;
  - lists the evidence that lives outside the store.
- `docs/validation/contract.yaml`
  - holds one verdict per audited row: class, verdict, audited revision and the files it
    depends on;
  - a verdict turns **stale** automatically when anything it rests on changes after its audited
    revision: its test, its listed reference files, or the code and data its level declares in
    `depends_on` (element sources for L1, the solver for L2, the mesh generator and the decks for
    L3, ...), including the levels it requires.
- `python tools/validation_matrix.py contract [--level Ln] [--rows]` reports, per level, what is
  trusted, stale, undecided, unarbitrated, unaudited, absent or over the ceiling. A level is **closed** only
  when nothing in it is open **and every level it requires is closed**. The tool enforces this
  and prints `blocked by Lk`. A gap blocks its level when it is `not_validated` or an
  `open_defect`. A `bounded` gap is a stated limit and shows as a caveat.

Day 0 (`a20a78b`): every row is **unaudited**. 280 of the 285 comparisons carry a
`justified: true` that extraction filled in, so that flag is not a judgement.

| Level | Scope | Store rows | Absent | Over ceiling | Single or no arbiter | Blocking gaps |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| L1 | MITC shell elements vs published | 103 | 2 | 3 | 55 | 0 |
| L2 | Linear solver vs published / CalculiX | 61 | 0 | 1 | 27 | 1 |
| L3 | IEA 15 MW blade structure | 24 | 9 | 2 | 32 | 1 |
| L4 | BEM vs OpenFAST AeroDyn | 12 | 3 | 0 | 15 | 0 |
| L5 | Coupled IEA 15 MW FSI with the built-in BEM (through preCICE) | 11 | 8 | 0 | 9 | 3 |
| L6 | Coupling with OpenFOAM through preCICE (planned; OpenFOAM side in progress) | 0 | 1 | 0 | 1 | 0 |

The "single or no arbiter" column is a heuristic. It reads the store's reference kinds and
treats a `paper` as published numerical data until a verdict says the paper reports a
measurement.

**Each level conditions the next**, and the tool enforces it through `requires`. The default
is every lower level; L4 (pure BEM against AeroDyn) requires none. A finding is argued on the lowest level where it appears,
and it suspends every conclusion above it that uses the same element, solver or path.

## How an item is worked

- One item per session, started from a clean context.
- **Entry** says what to read first. **Closes when** is a measurable state of the contract, not
  an activity.
- A closure ships three things in one push:
  1. the verdict in `contract.yaml`;
  2. the closure section;
  3. the comment on the item.
- An unknown found mid-item becomes a new item. It is never absorbed into the current one.

## Workstreams, in order

### A. Make the instrument trustworthy (first; cheap; it blocks nothing physical)

| # | Item | Closes when |
| --- | --- | --- |
| A1 | #32: the tool suite writes the real store | the suite never writes `docs/validation`; a test proves it |
| A2 | #24: `coherence` wipes prose and `measured` blocks | a rewrite preserves both; *stale* means a citation mismatch |
| A3 | #23 and the three failing `references` tests (line drift) | `references check` and the non-slow tool suite are green |
| A4 | #31: group 26 has six undeclared tests | `regression --group 26` exits 0 |
| A5 | Rows with `rtol 1.0` (groups 24, 25), probably extraction misreads | each is corrected or classed as an identity |
| A6 | Import the absent anchors (S-0..S-6, V-01/03/05, the #19/#26 guards, the BEM frame and polars) | `contract` reports 0 absent rows that exist as tests. A path-only absent entry cannot take a verdict; only importing it clears it |
| A7 | `AGENTS.md`: stale S-7 anchor; `test_rotor_inertial.py` is said to be removed but exists | the text matches the tree |

### B. Review, level by level (bottom-up; produces verdicts, not fixes)

| # | Item | Closes when |
| --- | --- | --- |
| B1 | Audit L1 | every L1 row has a verdict at one revision, with its arbiters named; the documented reds (D-Tube, UL elastica, ko2017's 8 failures, frame objectivity) are each an exception with a source, or rejected; every `atol` comparison (23 in the store) is classed or given a relative equivalent, because the ceiling cannot judge an absolute bound |
| B2 | Audit L2 | as B1; each CalculiX comparison is shown converged **in both codes** |
| B3 | Audit L3 | as B1; S-2 (+8.1 %) and V-02 first edge (10 %) carry a quantified exception or are rejected |
| B4 | Audit L4 | as B1 |
| B5 | Audit L5 | as B1, once the D-items for L5 exist |
| B6 | Audit L6 | as B1, once the OpenFOAM coupling exists |

### C. Physics questions already open (each attacked only once its level's audit has reached it)

| Level | Item | Issue |
| --- | --- | --- |
| L3 | Blade torsional stiffness: S-1 sectional is +23.9 % stiffer, S-7 global is ~3.8 % softer (after `72f64de`); the sign disagreement is unattributed | #20 (reopened) |
| L3 | S-7's old table, the residual between 1.014 and 1.286 | #22 |
| L3 | The rated-twist docstring's shear centre (0.477c) and the pre-fix application figures | #25 |
| L5 | Wall-flow moment realisation: the activation decision | #16 |
| L5 | Per-process node ordering in the mesh generator | #17 |
| L5 | Parked V50 divergence (2168 m against ~8 m) | #10 |

### D. Missing evidence to build (after B1–B4)

| Level | Item |
| --- | --- |
| L5 | Conservation of force and work (`Σf`, `Σf·u`) across the preCICE mapping, independent of the fluid participant |
| L6 | A published FSI benchmark with OpenFOAM. Candidate: the preCICE perpendicular flap. Feasibility is unverified: there is no SOLID family, so the flap must run as a shell or PLANE strip |
| L5 | A coupled rated run against OpenFAST coupled (ElastoDyn/BeamDyn + AeroDyn) at the same operating point |
| L5 | A published IEA 15 MW coupled reference to replace Zhou 2025 (declared non-transferable in #14). It is the second, independent arbiter next to OpenFAST coupled |

## Decisions recorded (2026-10-10)

- **L5 and L6 were swapped.** The coupled IEA 15 MW rotor with the built-in BEM has references
  today (OpenFAST coupled, published FSI) and couples through preCICE, so it is L5. Coupling
  with OpenFOAM is separate work in progress, so it is L6.
- **Experimental validation is an arbiter, not a level.** It left L6. The arbiter rule above
  replaces it.
- **Only unbounded gaps block.** A gap that is `not_validated` or an `open_defect` keeps its
  level open. A `bounded` gap is shown as a caveat.

## Standing caveats (carried over from #18)

- Coupled numbers older than `7a84da1` (load frame) and `6765633` (projection geometry) cannot
  be cited.
- Where no experimental data exists, validation rests on closed forms or on at least two
  independent numerical arbiters. The store gap `experimental_validation` records that no
  experiment arbitrates the coupled rotor.
- A digitized figure is a reference with a quadrature: check its trapezoid against the paper's
  own integrals.

## What #18 hands over

| #18 item | Here |
| --- | --- |
| #16 | C (L5) |
| #10 | C (L5) |
| #17 | C (L5) |
| #24, #23, #31, #32 | A2, A3, A4, A1 |
| #25, #22 | C (L3) |
| #20 (closed in #18, reopened since) | C (L3) |
| Closed items #2–#30 | audited on 2026-10-10 (`odd/tasks/issue18-closure-audit.md`); they enter the contract as unaudited rows and get verdicts in B1–B6 |

## Closures

Each closure records what was measured and how to re-check it. A tooling item carries no
`contract.yaml` verdict: it touches no store row.

### A1 — #32: the tool suite writes the real store (proposed for closing 2026-10-10)

- **Cause.** Two tests ran `coherence` on `docs/validation`, not one:
  `test_every_group_re_derives_to_the_rows_on_disk` (slow) and `test_the_real_store_is_coherent`
  (group 4, *not* slow, so the `-m "not slow"` rule in #32 was not safe either).
- **Fix.** Both run on `real_store_copy` (a `copytree` into `tmp_path`). The autouse fixture
  `store_is_never_written` in `tools/tests/conftest.py` snapshots every file of the shipped store
  before each test. After the test it fails the test if any file was added, removed or rewritten,
  and puts the bytes back. Every current and future test in `tools/tests` is held to it.
- **Proof.** `tools/tests/test_store_guard.py` drives the guard against a scratch tree: identical
  bytes are no change; an added, a removed and a modified file are each named and restored. A
  green suite alone proves nothing here.
- **Measured.** At `f5405b8` plus this change, `scripts/aeroenv.sh python -m pytest -o addopts=""
  -q tools/tests` (slow included, 4 min 52 s): 156 passed, 4 failed. `git status --porcelain --
  docs/validation` was empty before and after.
- **Failures, none from A1:**
  - Three `references` tests (`cited_by_stale`, `check_reports_a_stale_declared_site`,
    `where_used_lists_rows_and_code`): line drift, item A3.
  - The coherence sweep, now red honestly. Five row files cite moved lines: `6-tube_torsion`,
    `20-blade_anchor_beam`, `31-tube_moment_realization`, `34-blade_s7_torsion` and
    `35-tube_projection`. Before this change the test rewrote them in place and dropped their
    `measured` blocks, which is the damage #32 describes. Refreshing them is the maintainer's
    explicit `coherence`. It should wait for A2 (#24), so that the rewrite keeps prose and
    measurements.

### A2 — #24: `coherence` wipes prose and `measured` blocks (proposed for closing 2026-10-10)

- **Cause, two parts.**
  - `coherence` re-ran `extract --write`, which rebuilds the whole file from the extractor's output
    and copies back only `measured`, `expected`, `flags` and `history`. Hand-written `validates`,
    `notes`, `tolerance.justification` and `reference.label` became placeholders. Stale meant
    "the bytes differ", so prose alone made a file stale.
  - `preserve_measurements` paired comparisons by `<file>:<line>`. An edit above the assertions
    shifts every line, so groups 6, 20 and 31 lost every `measured` block on a pure shift.
- **Stale is now a citation mismatch.** `refresh_citations` starts from the file on disk and
  overwrites only what the extractor owns: `id`, `group`, `title`, `tests`, and per comparison
  `label`, `asserted`, `tolerance.kind/value/source` and `reference.kind/citation`. Everything else
  is kept. The file is written only when an owned field changed. A comparison is carried only
  when `pair_comparisons` pairs it with a stored one.
- **Pairing contract (changed, stated here on purpose).** The old docstring refused to carry a
  measurement across a moved line. The new rule:
  - **No source moved in the row:** pair by `<file>:<line>`, in order among comparisons that share
    a line. A pair whose assertion changed is refused.
  - **A source moved:** pair by a line-free fingerprint (label without `at line N`, keeping the loop
    binding; tolerance kind and bound; `asserted`; reference kind). Pairing is in order, and only
    for a fingerprint that occurs as many times on disk as in the code. If the count changed,
    nothing is paired, because nothing says which identical comparison is the new one. So a
    margin is never slid onto a neighbour's reference.
  - A carried `measured.run` is left untouched, so it still names the revision that produced the
    number.
- **Not changed:** `extract --write` still rebuilds a row file from the code. It now pairs
  measurements the new way, but still replaces prose. That is its documented job ("a better
  description").
- **Proof.** Six tests in `tools/tests/test_validation_matrix.py`:
  - a moved line keeps its measurement;
  - an ambiguous move carries nothing;
  - `refresh_citations` moves citations, keeps prose and margins, and refuses a changed bound;
  - `coherence` on a store copy leaves hand-written prose byte-identical and reports *coherent*;
  - `coherence` refreshes a moved citation in place and keeps the prose and the margin;
  - the "no counterpart" message test now uses a changed bound, because a moved line alone is
    carried.
- **Measured.** At `d117877` plus this change:
  - Sweep on a copy, then on `docs/validation`: **4 of 35** stale (6, 20, 31, 34). The diff holds
    only `label`/`source` lines, plus one real structural change in group 6. There, line 614 now
    asserts inside a loop over thickness, so one stored comparison (never measured) became two
    new ones.
  - `status: measured` counts are unchanged: group 6 has 8, group 20 has 4, group 31 has 2, group 34
    has 1. S-7's hand-written prose (group 34) and group 6's hand-written `reference.label` survive.
  - **Group 35 is not stale.** A1 counted it only because its P1a prose differed from the
    extractor's placeholders. Its citations were right.
  - `check`: 211 rows, 286 comparisons (+1 from the group-6 split), 0 errors.
  - `pytest tools/tests -m slow`: the coherence sweep **passes** (red since A1).
  - `-m "not slow"`: 161 passed, 3 failed. All three are the A3 `references` tests.
  - `scripts/check.sh quick` is OK. `git status --porcelain -- docs/validation` shows only the
    four refreshed files.
- **Checked against the files' own history, not only HEAD.** The hashes cited in #24 do not touch
  these files (they are probably from before a rebase). So each file's prose and `measured` count
  were compared through its `git log`:
  - Group 34 matches `19cf8be`, the S-7 re-measure: same prose, 1 measured.
  - Group 20 matches `9d71c94`: same prose, 4 measured.
  - Group 6 matches `c0c5587`: 8 measured. The only row that differs is the group-6 split.
  - Group 35 is untouched since `6780835`.
- **Contract.** The only verdict on these rows is S-7's `undecided` (row 34). `contract` still
  reports it as undecided, not stale, and no level changed state. No verdict is needed: only
  citation fields changed. The store now holds 286 comparisons, not #33's 285.
- **Known limit.** `reference.label` is kept from disk, because a row does not record whether the
  code or a person wrote it. Editing a test's `reference_name=` is therefore invisible to
  `coherence`.
