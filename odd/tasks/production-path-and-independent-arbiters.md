# Feature: #14 re-planned — measure on the production path, and add independent arbiters

Status: planned 2026-10-08 (P1/P2 open; nothing measured yet on these lines)
Owner: next session
Related: issue #14 (roadmap item 2 of #18), `odd/tasks/rated-twist-sign-convention.md`,
`odd/tasks/shear-centre-arbitration.md`, `odd/tasks/application-exact-case-validation.md` (the units
that came before), `odd/tasks/rotor-fsi-production-wiring.md` (the same rule, applied to the rotor).

## Why this exists: the measurements so far are test-local

The issue's own routes are exhausted for what can be read off the paper, and every *unexplained*
term has been measured. But the magnitudes that carried those measurements were produced by
**constructions written inside the test file**, not by the production path:

| measured quantity | how it was produced | production? |
| --- | --- | --- |
| the 6.17x spread across "the four applications" | `_rated_load_cases`, four hand-built load vectors | no |
| `mp_only` `+4.8337`, `at_ac` `+15.4921` deg | the same, solved with `spsolve` on the fixture's K | no |
| the "nose-down couple" sign test | `_section_couples`, a `+F`/`-F` pair | no |
| the tube validation (`0.23%` / `0.11%`) | `_ring_traction` + `_distributed_pattern_load`, written in that test | no |

**And the production path has exactly one application rule.** So:

- the **quotable** rated twist is the one `ForceProjector` produces, and it is a *single* number, not
  a spread: `test_rated_twist_under_production_loads` already asserts its sense, and the production
  participant is the path the campaigns use;
- the 6.17x is a spread across **test-local alternatives**, i.e. a sensitivity study, and quoting it
  as "the application's spread" overstates what was measured;
- the tube case validates **the pattern this repository's test wrote**, not the projector that the
  blade path actually uses. The verifier said exactly this ("the test validates the pattern, not the
  production path") and the response was to narrow the prose, not to fix the measurement. That was
  the wrong trade.

The repository already has the rule, in `odd/tasks/force-projection-shear-flow.md`'s constraints:
"`src/aeroelast/solvers/bem/force_projection.py` is production; this issue authorises changing it,
with a test that exercises the production entry point (not a test-local copy of the rule)." And
`odd/tasks/rotor-fsi-production-wiring.md` exists because 23 tests built `K` with the production
assembler and none invoked a production solver. This unit applies the same discipline to this issue's
measurements.

## P1 — measure on the production path

- [ ] P1a **The tube case through `ForceProjector`.** The projector takes `(mesh, BladeAero, span
  direction, element properties)`, so a tube-like `BladeAero` (stations with chord and a polar per
  station) makes the rectangular tube projectable. Then the same closed-form comparison runs through
  **the production entry point**, and the `0.23%` / `0.11%` claim attaches to production rather than
  to a test-local pattern. Reuse the module's existing tube mesh, properties, Bredt reference and
  metrics; write only the adapter.
- [ ] P1b **The rated twist through the production path, as the headline number.** State the rated
  tip twist that `ForceProjector` produces and what it is against Zhou's `-3.60 deg`. That is the
  quantity the write-up would quote, and it is one number rather than four.
- [ ] P1c **Re-frame the sensitivity honestly.** Keep the four applications, but say what they are: a
  spread across test-local constructions, measured to show what an *inconsistent* distribution would
  cost. The exact-case validation is what bounds the production path's pattern, not the spread.
- [ ] P1d Re-check the sign work's claims under the same lens: the cross-path test did use the
  projector (as the arbiter, and the test-local side was the defective one), so those results stand —
  but say so explicitly rather than leaving it implied.

## P2 — independent arbiters against the paper

Two of the three route-(c) references the issue names are available in this repository and unused for
this question.

- [ ] P2a **CalculiX S8R shell of the same blade under the same loads.** The writer, the decks and the
  S8R parity machinery exist (`tests/validation/parity/`, `docs/validation/rows/17-beam_shell_4cases.yaml`,
  `tools/ccx_blade_twist_arbitration.py`). This is the issue's own third route: a reference whose
  model *is* published, on the same geometry and loads, so the rated twist gets a structural arbiter
  that does not depend on our shell implementation at all.
- [ ] P2b **OpenFAST's flexible rated response.** The deck already arbitrates the *stiffness* (S-7 /
  #20). Running the rated case flexibly gives the *coupled* response — the quantity the paper reports
  and the one our one-way application approximates — so the `1.61x` gets a second arbiter on the
  aerodynamic-plus-structural side. Note the standing constraint from S-7: no coupled number is
  citable while item 0 (#19) is broken; an OpenFAST comparison is a third-party code's own coupled
  solution, so it is not blocked by our coupling, but our side's numbers would still be.
- [ ] P2c **Say which arbiter decides what.** Three different references answer three different
  questions (stiffness, structural response, coupled response). The write-up needs the assignment
  stated, or a reader will compare incomparable numbers — which is exactly what the S-7 table (#20)
  and the `origin_main` row did.

## What changes in the record

- Anything that presented a test-local construction's number as "the application" gets re-stated with
  the path it came from. The sign defects do not change (production was the arbiter there).
- The tube validation's claim narrows to what P1a makes true, or stays explicitly pattern-scoped.
- The 6.17x becomes "spread across test-local applications", which is what it is.
- Comments already posted on #14 that overstated the scope get a correction; the repository's own
  practice in this session has been to correct such statements rather than leave them.

## Risks

- **A tube is not a blade.** The projector expects a `BladeAero` (chords, polars, a hub radius); a
  rectangular tube can be dressed as one, but the adapter must be honest about the artificial inputs
  it invents (a flat plate polar, a chord equal to the tube's width) and the test must say so, or the
  production-path claim inherits a fake input.
- **CalculiX is a big step.** The blade deck writes and runs today for the static parity cases; the
  rated twist with the torsional load path may need the load realisation carried across, which is a
  unit of its own. Establish the shell-only torsional response first if the full load path is heavy.
- **OpenFAST answers a different question than the paper does.** Their torsional quantity is
  definitionally different (issue #14's own finding), so an OpenFAST comparison arbitrates *our*
  response, not the paper's number. That is still the useful part: it separates "our model is wrong"
  from "the quantities are not the same".
