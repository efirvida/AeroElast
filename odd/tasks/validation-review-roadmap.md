# [validation] general review roadmap: what to fix, review and build before the solver is citable

*Draft of the issue that replaces #18 (2026-10-10). This is a plan, not a set of results. Nothing
below is fixed by writing it down.*

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
  trusted, stale, undecided, unaudited, absent or over the ceiling. A level is **closed** only
  when nothing in it is open **and every level it requires is closed**. The tool enforces this
  and prints `blocked by Lk`.

Day 0 (`a20a78b`): every row is **unaudited**. 280 of the 285 comparisons carry a
`justified: true` that extraction filled in, so that flag is not a judgement.

| Level | Scope | Store rows | Absent | Over ceiling | Gaps |
| --- | --- | ---: | ---: | ---: | ---: |
| L1 | MITC shell elements vs published | 103 | 2 | 3 | 0 |
| L2 | Linear solver vs published / CalculiX | 61 | 0 | 1 | 1 |
| L3 | IEA 15 MW blade structure | 24 | 9 | 2 | 1 |
| L4 | BEM vs OpenFAST AeroDyn | 12 | 3 | 0 | 0 |
| L5 | Coupling between domains (preCICE) | 11 | 5 | 0 | 2 |
| L6 | Coupled IEA 15 MW FSI | 0 | 4 | 0 | 3 |

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
| B1 | Audit L1 | every L1 row has a verdict at one revision; the documented reds (D-Tube, UL elastica, ko2017's 8 failures, frame objectivity) are each an exception with a source, or rejected; every `atol` comparison (23 in the store) is classed or given a relative equivalent, because the ceiling cannot judge an absolute bound |
| B2 | Audit L2 | as B1; each CalculiX comparison is shown converged **in both codes** |
| B3 | Audit L3 | as B1; S-2 (+8.1 %) and V-02 first edge (10 %) carry a quantified exception or are rejected |
| B4 | Audit L4 | as B1 |
| B5 | Audit L5 | as B1, once C-items for L5 exist |
| B6 | Audit L6 | as B1, once D-items exist |

### C. Physics questions already open (each attacked only once its level's audit has reached it)

| Level | Item | Issue |
| --- | --- | --- |
| L3 | Blade torsional stiffness: S-1 sectional is +23.9 % stiffer, S-7 global is ~3.8 % softer (after `72f64de`); the sign disagreement is unattributed | #20 (reopened) |
| L3 | S-7's old table, the residual between 1.014 and 1.286 | #22 |
| L3 | The rated-twist docstring's shear centre (0.477c) and the pre-fix application figures | #25 |
| L5 | Wall-flow moment realisation: the activation decision | #16 |
| L5 | Per-process node ordering in the mesh generator | #17 |
| L6 | Parked V50 divergence (2168 m against ~8 m) | #10 |

### D. Missing evidence to build (after B1–B4)

| Level | Item |
| --- | --- |
| L5 | Conservation of force and work (`Σf`, `Σf·u`) across the preCICE mapping, independent of the fluid participant |
| L5 | A published FSI benchmark with OpenFOAM. Candidate: the preCICE perpendicular flap. Feasibility is unverified: there is no SOLID family, so the flap must run as a shell or PLANE strip |
| L6 | A coupled rated run against OpenFAST coupled (ElastoDyn/BeamDyn + AeroDyn) at the same operating point |
| L6 | A published IEA 15 MW coupled reference to replace Zhou 2025 (declared non-transferable in #14) |

## Open scope decision (maintainer)

A level closes only with no gap listed, whatever the gap's status. Two consequences need a
decision before B-audits close levels:

- `experimental_validation` sits on L6, so **L6 can never close**. The stated scope is
  published results, theory and other codes, not experiment.
- `bounded` gaps, such as L2's `composite_stress_recovery`, block the same way as
  `not_validated` ones.

Options: a gap blocks its level, or it is shown as a non-blocking caveat. This can be decided
per status or per gap.

## Standing caveats (carried over from #18)

- Coupled numbers older than `7a84da1` (load frame) and `6765633` (projection geometry) cannot
  be cited.
- Verification here is against models, codes and literature. `experimental_validation` is an
  open gap that this list does not reach.
- A digitized figure is a reference with a quadrature: check its trapezoid against the paper's
  own integrals.

## What #18 hands over

| #18 item | Here |
| --- | --- |
| #16 | C (L5) |
| #10 | C (L6) |
| #17 | C (L5) |
| #24, #23, #31, #32 | A2, A3, A4, A1 |
| #25, #22 | C (L3) |
| #20 (closed in #18, reopened since) | C (L3) |
| Closed items #2–#30 | audited on 2026-10-10 (`odd/tasks/issue18-closure-audit.md`); they enter the contract as unaudited rows and get verdicts in B1–B6 |
