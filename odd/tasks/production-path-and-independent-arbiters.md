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
  a spread. **Corrected 2026-10-08: it is one number per realization, two in total.** The coupled
  campaigns run through `BEMFSIParticipant`, which passes no `element_properties`, so a strip takes
  the minimum-norm `_distribute` fallback and measures `+8.1048 deg`; `standalone.py` passes the
  config's `elements.properties`, takes the multi-cell wall flow and measures `+9.6669 deg`. Both are
  production; which one applies is a property of the caller, not of the projector. See the Work
  units section;
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

- [x] P1a **The tube case through `ForceProjector`.** The projector takes `(mesh, BladeAero, span
  direction, element properties)`, so a tube-like `BladeAero` (stations with chord and a polar per
  station) makes the rectangular tube projectable. Then the same closed-form comparison runs through
  **the production entry point**, and the `0.23%` / `0.11%` claim attaches to production rather than
  to a test-local pattern. Reuse the module's existing tube mesh, properties, Bredt reference and
  metrics; write only the adapter. **Probed 2026-10-08** (`$SCRATCH/p1a_diag/probe.py`, non-committed
  diagnostic); the adapter works and the measurement is in the Work units section below. **Done**:
  `tests/validation/parity/test_thin_walled_tube_projection.py`, store group 35.
- [x] P1b **The rated twist through the production path, as the headline number.** DONE: the number
  is stated in the Work units section below - `+8.1048 deg` for the coupled path the campaigns use
  (2.25x Zhou's `-3.60 deg` in magnitude, opposite sign by frame convention), `+9.6669 deg` for the
  standalone multi-cell caller. The comparison carries the assignment Zhou is a **beam** (GEBT), so
  a shell magnitude against their beam number is not transferable; the arbitrable counterpart is
  beam-vs-beam. This unit also corrected two record defects: the "one number" premise above and the
  `production_rated_loads` fixture citing `standalone.py` while configuring the participant.
- [x] P1c **Re-frame the sensitivity honestly.** DONE 2026-10-08: the four applications stay, now
  named as **test-local constructions written in `test_blade_rated_twist.py`**, i.e. a sensitivity
  study of what an *inconsistent* distribution costs. Every live place that read as a property of
  the load path was re-framed (the module docstring, both four-application test docstrings, the
  withdrawn-magnitude print label and the `spread_omega > 2.0` comment), each pointing at store
  group 35 as what bounds the production pattern; the historical note carries a dated re-frame and
  the closure log the same qualifier. **No assertion, tolerance or numeric value changed** -
  `ruff` clean, the two affected tests `2 passed`, `regression --group 27 --write` `7 same` with
  only the call-site lines and the source digest moved, `check` 198 rows / 259 comparisons / 0
  errors. The stale mesh shear-centre `0.477` paragraph in the same docstring was **left alone**:
  it is the separately tracked follow-up, not this unit.
- [x] P1d Re-check the sign work's claims under the same lens. **Done 2026-10-08**: the results
  stand - the cross-path test `test_the_two_moment_applications_agree_in_sign` did exercise the
  production projector, and the defective side was the test-local one - and both the test docstring
  and the sign section of `docs/validation_closures.md` now say so explicitly, naming the projector
  **the arbiter** and its chain as anchored to measured geometry
  (`tools/diagnose_leading_edge.py`: leading edge at `+x`, downwind at `+y`, `+z` moves the leading
  edge downwind) rather than to the cross-check it participates in. Prose only; `regression --group
  27 --write` re-measured the same seven margins.

## P2 — independent arbiters against the paper

Two of the three route-(c) references the issue names are available in this repository and unused for
this question.

- [x] P2a **CalculiX S8R shell of the same blade under the same loads.** **Already measured
  2026-10-06; recorded here 2026-10-08.** `tools/ccx_blade_twist_arbitration.py` feeds an
  independent FE code **the exact same nodal force vector** the production `ForceProjector`
  produces, clamped the same way, and converges each code on its own mesh sequence (same-mesh
  MITC4-vs-S8R is invalid): at 0.25 m, AeroElast `+8.0980` vs CalculiX `+8.4396` = **4.22%**; the
  23.77% / 27.91% at 1.00 / 0.50 m was discretisation. Evidence:
  `odd/tasks/force-projection-shear-flow.md:202-204`, `tools/ccx_blade_twist_arbitration.py`
  (header + `--element-size`). So the structural arbiter **exists**; what P2a still owes is its
  registration as *the* arbiter for the rated twist in the assignment below.
- [ ] P2b **OpenFAST's flexible rated response.** The deck already arbitrates the *stiffness* (S-7 /
  #20). Running the rated case flexibly gives the *coupled* response — the quantity the paper reports
  and the one our one-way application approximates — so the `1.61x` gets a second arbiter on the
  aerodynamic-plus-structural side. Note the standing constraint from S-7: no coupled number is
  citable while item 0 (#19) is broken; an OpenFAST comparison is a third-party code's own coupled
  solution, so it is not blocked by our coupling, but our side's numbers would still be.
- [x] P2c **Say which arbiter decides what.** DONE 2026-10-08: the assignment is stated in the
  section at the end of this document, every row carrying its live evidence and, as importantly,
  what the arbiter does **not** decide. It is stated here first because it gates every quotation;
  propagating it to the write-up (the article draft / the validation docs) is a separate,
  maintainer-owned edit.

## What changes in the record

- Anything that presented a test-local construction's number as "the application" gets re-stated with
  the path it came from. The sign defects do not change (production was the arbiter there).
- The tube validation's claim is now **production-scoped**: group 35 runs the same closed-form
  comparison through `ForceProjector.project()`, so the pattern claim no longer rests on a
  construction written in the test. The pattern test's own docstring stays accurate for itself (it is
  the test-local comparator) and was left untouched; the production companion is the new module and
  group 35.
- The 6.17x becomes "spread across test-local applications", which is what it is. **Done
  2026-10-08**: the live texts say it in
  `tests/validation/blade/test_blade_rated_twist.py` (module docstring, both four-application test
  docstrings, the withdrawn-magnitude print and the `spread_omega > 2.0` comment), in
  `docs/validation_closures.md` and in `odd/tasks/application-exact-case-validation.md`, each
  pointing at store group 35 as what bounds the production pattern.
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

## Work units — P1a (status 2026-10-08)

Measured on `$SCRATCH/p1a_diag/probe.py` before the first source write, so the bound is fixed and
the adapter is known to work. The probe builds a `MeshModel` from `_tube_mesh`, a tube `BladeAero`
(chord `B`, `aerodynamic_center=0.5`, flat fake polar - the artificial inputs are declared), a
`BEMResult` with `Np = FORCE_PER_LENGTH`, `Tp = 0`, `Mp = e * FORCE_PER_LENGTH`, and passes them
through `ForceProjector(mesh, blade_aero, hub_radius=0.0, element_properties={"shell":_iso_prop(t)})`,
then `.project()`. Clamped-free, `f = 1e5 N/m`, `e = 0.2 m`:

| wall | `theta_z` | vs `e f L^2/(2 GJ)` | vs the discrete Saint-Venant response | torque ruler |
| ---: | ---: | ---: | ---: | ---: |
| 10 mm | `-1.021242e-3` | `1.0311` (3.11% high) | **`0.997801`** (0.220% low) | `-1.24e5` |
| 0.5 mm | `-2.045078e-2` | `1.0324` (3.24% high) | **`0.999069`** (0.093% low) | `-1.24e5` |

Two findings, both measured:

1. **The production pattern is as accurate as the test-local one.** Against the analytic response to
   the projector's *own* discrete station load - couples `e f dr_k` at `z_k`, `theta(L) = (e f/GJ)
   sum_k dr_k (L - z_k)` - the production path lands at 0.220% / 0.093%, the same numbers the
   test-local pattern gave (0.23% / 0.11%). The application is validated through the production entry
   point, not through a test-local copy.
2. **The 3.1-3.2% continuum excess is the projector's band quadrature, not the pattern.**
   `ForceProjector.__init__` extends the outer strip bands half a spacing past the first and last
   station (`force_projection.py:259-263`), so `sum_k dr_k = L + DZ = 6.2` here, and the projected
   load's own torque ruler reads exactly `e f (L + DZ) = 1.24e5`, not `e f L = 1.2e5`. The remainder of
   the 3.11 vs 3.33 difference is the clamped-end boundary layer.

Control: with `Mp = 0` (force through the section centre) the production load twists the tube by
`1.97e-15` rad (10 mm) / `-8.34e-12` rad (0.5 mm) - round-off - and `verify().force_error = 2.3e-10 N`.

- [x] WU-P1a-1 **The adapter and the production test.** DONE:
  `tests/validation/parity/test_thin_walled_tube_projection.py` (new file, its own group so
  `extract --write` cannot regenerate group 6's rows, the same reason group 31 exists). Four
  `assert_relative_error` call sites, `TOL = 0.05` fixed before the run: signed `theta_z` against
  both the continuum closed form and the discrete Saint-Venant response, per wall. `1 passed`,
  `ruff` clean.
- [x] WU-P1a-2 **Store rows.** DONE: group 35 registered in `docs/validation/groups.yaml`, its
  residual pattern declared in `docs/validation/residual-patterns.json`, `extract --group 35
  --write` + `regression --group 35 --write` (4 measured), `check` 198 rows / 259 comparisons / 0
  errors, `status` 0 files neither grouped nor declared.
- [x] WU-P1a-3 **Re-scope the record.** DONE: this document's "What changes in the record" now
  states the tube claim is production-scoped; the new module's docstring carries the two findings and
  the artificial inputs. `test_thin_walled_tube_torsion.py` was deliberately **not** edited: its
  pattern test's docstring claims nothing false, and editing it would move group 6's call-site line
  numbers for no claim correction.
- [x] WU-P1a-4 **Verify and commit.** DONE: independent read-only verifier confirmed all nine
  claims (production path, bound, the four numbers, the quadrature attribution re-derived from the
  code, reference independence, declared artificial inputs, sign pinning, store integrity, scope);
  work-unit commit `6780835` on `integrate/origin-main-2026-09-30`. The only caveat it raised - the
  row file's stale `notes` boilerplate - was fixed before the commit.

Environment: `export LD_LIBRARY_PATH=/petrobr/app_sequana/gcc/14.2.0/lib64:$LD_LIBRARY_PATH`
(without it `_aeroelast` fails to import; `module load` is not required for the import, only the
libstdc++ with `CXXABI_1.3.15`).

## Work units — P1b (status 2026-10-08)

Measured read-only at HEAD `5e2832c` with the three production tests of
`tests/validation/blade/test_blade_rated_twist.py` (`3 passed in 33.75 s`, `-k production`). The
rated tip twist the **production** load path produces, and what it is against Zhou et al. 2025
Table 4 (`-3.60 deg`, a GEBT **beam** model):

| production call site | realization | `omega` (tip section) | vs Zhou | `theta_z` | distortion/`\|omega\|` |
| --- | --- | ---: | ---: | ---: | ---: |
| `BEMFSIParticipant` (the coupled campaigns) | minimum-norm `_distribute` fallback | **`+8.1048 deg`** | **2.2513x** | `+10.6994 deg` | 1.6728 |
| `standalone.py` with `elements.properties` | multi-cell wall flow | `+9.6669 deg` | 2.6853x | - | - |

Other measured quantities on the coupled path: applied `|sum(F)|` = `852562.70 N` against
`bem.thrust/3 = 841688.79 N` (`+1.292%`); flapwise tip deflection `+17.2256 m` (1.2428x Zhou's
`13.86 m`); edgewise `-1.9926 m` (1.6333x Zhou's `-1.22 m`). All three tests pass and the
promotion guard (`abs(1.0 - ratio_to_zhou) > 0.05`) is live and far from firing.

**Premise corrected.** This plan said the production twist is "a *single* number, not a spread". It
is **one number per realization, two in total**: a strip only takes the multi-cell wall flow when
the caller hands `ForceProjector` an `element_properties` map. `standalone.py:91-100` does;
`fsi_participant.py:357-371` and its `_rebuild_projector` at `:848-858` explicitly do not, because
the fluid side never builds the laminate map - so the coupled campaigns, the path this issue is
about, are the fallback and the headline is `+8.1048 deg`.

**The assignment, so the number is not read as a model error.** Zhou's quantity is a **beam's**
rotation about its reference axis; ours is a **shell's** tip section rotation, which carries the
section distortion the beam cannot represent (`distortion/|omega| = 1.6728`). A shell magnitude
against their `-3.60 deg` is therefore **not transferable**, and the live repo statement of that is
`docs/validation_closures.md:244` ("Twist del shell ~5x la viga: fisica del modelo ... Las
referencias beam-based (Zhou/Ma/BeamDyn) son cota inferior"). The arbitrable counterpart is
beam-vs-beam: under Zhou's own Fig. 11 loads the shell gives `+5.7822 deg` (1.61x) and **our beam
`-2.0790 deg` against their `-3.60 deg`**, a `0.58x` twist (`tools/diagnose_zhou_loads_reverse.py`);
the `1.61x` is that tool's shell figure and must never be quoted as the production rated number. The
"~2.03x torque deficit" that an earlier session note carries is **not in the tree** and is not cited
here.

**Record defect fixed.** The `production_rated_loads` fixture docstring cited `standalone.py` as the
production construction it mirrors, while configuring the projector the way the participant does
(default directions, no `element_properties`). Corrected to cite `fsi_participant.py`, the coupled
path, at `tests/validation/blade/test_blade_rated_twist.py:1048` - one line, no line-number shift, so
group 27's call-site rows did not move.

**Prose audit.** Every quoted production figure (`+8.1048`, `+9.6669`) in the test docstrings
(`:234`, `:1128`, `:1141`, `:1155`, `:1177`, `:1183`), `odd/tasks/rated-twist-sign-convention.md:41`
and `docs/validation_closures.md:27` matches the live run; nothing was stale. The `1.61x` figures are
correctly scoped to the Fig.-11 tool run, not the production rated number.

**Not a row.** The magnitude stays withdrawn from the store: it is compared with a beam number under
unknown stiffness source and pitch, so no defensible bound exists. Group 27 keeps the two Table-6
paper comparisons, the four self-resultant invariants and the applied-load invariant.

## P2c — which arbiter decides what (2026-10-08)

Settled **before** any magnitude is quoted, because three of the references answer three different
questions and reading one as another is exactly what the S-7 table (#20) and the
`origin_main_integration` row did. Each row carries its live evidence and, as importantly, what that
arbiter does **not** decide.

| question | arbiter | what it decides | live evidence | what it does **not** decide |
| --- | --- | --- | --- | --- |
| Is the load **application pattern** right? | Bredt closed form on the closed tube, driven through the production `ForceProjector.project()` | the pattern (transverse force + its transfer moment): `0.22%` / `0.09%` of the exact discrete response, `0.03` quadrature on the continuum comparand | store **group 35**, `tests/validation/parity/test_thin_walled_tube_projection.py`, commit `6780835` | the blade's geometry, any magnitude against a paper |
| Is our **sectional stiffness** right? | OpenFAST MBC3 / BeamDyn deck, rotating modal (S-7) | `GJ` and the sectional stiffness: shell `+6.4%` global, `+23.9%` sectional | group 34, issue #20, `docs/validation_closures.md` S-7 | the twist magnitude, the load application |
| Is our **structural response** right for a given load? | CalculiX S8R, same blade, **the same nodal load** | the shell's response to that load: MITC4 `+8.0980` vs S8R `+8.4396` = **4.22%** at the converged 0.25 m mesh (each code converged on its own sequence) | `tools/ccx_blade_twist_arbitration.py`, `odd/tasks/force-projection-shear-flow.md:202-204` | the load itself (it is fed our nodal vector) and the coupled response |
| Is our **coupled** response right? | OpenFAST flexible rated (**P2b, not yet run**) | the aeroelastic response - the quantity Zhou approximates | pending | the sectional stiffness (that is S-7's), and the paper's definition of "torsion" |
| Is our magnitude **comparable with the paper**? | Zhou et al. 2025 Table 4 | **nothing about our model**: it is a GEBT beam number and our shell carries the section distortion a beam cannot (`distortion/\|omega\| = 1.6728`) | `docs/validation_closures.md:244` | any shell magnitude; the arbitrable pair is **beam-vs-beam** (our `-2.0790` vs their `-3.60`) |

Rules that follow, and that the write-up must keep:

1. **Every quoted number names its arbiter.** A magnitude without one is a number without a question.
2. **The shell rated magnitude is cited against CalculiX** (same load) and **never against Zhou
   directly.** The production number is `omega = +8.1048 deg` on the coupled path the campaigns use
   (fallback realization; `+9.6669 deg` multi-cell), and the arbiter that makes it a *measurement* is
   the 4.22% agreement with S8R under the identical load - not the paper.
3. **Zhou is cited for the beam-vs-beam comparison only**: under their own Fig. 11 loads our beam gives
   `-2.0790 deg` against their `-3.60 deg`, a `0.58x` twist (`tools/diagnose_zhou_loads_reverse.py`).
   The `1.61x` is that tool's *shell* figure (`+5.7822 deg`) and is not the production rated number.
4. **A missing arbiter is stated as missing.** P2b has not run, so no coupled number of *ours* is
   citable; OpenFAST's own coupled solution is a third-party result and stands on its own, but it
   arbitrates our side only once the coupled gates (#19 item 0) clear.
