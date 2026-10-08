# Feature: S-7's torsion ratio says both directions at once (#20)

Status: delivered 2026-10-07 (T1-T6 done)
Owner: this session (2026-10-07)
Related: issue #20 (roadmap item 8 of #18), roadmap item 2 (#14) which this gates,
roadmap item 1 (#13, closed) whose twist-only row this magnitude sets, issue #19
(item 0, unaffected: nothing here needs a coupled run).
Deliberately independent of the coupled path: the whole feature is offline, 81 s of
PETSc per run, no preCICE.

## Why this exists

S-7 (`tests/validation/blade/test_iea15mw_s7_torsion.py`) compares the shell's global
torsional response against the BeamDyn `GJ(z)` tables of the official IEA-15-240-RWT
deck, integrated as a Timoshenko torsion beam. It carries two numbers that point in
opposite directions: the live assertion measures the shell **stiffer** than the deck,
and the convergence table committed in the same file says it is **softer**.

The roadmap placed this item after item 0 and **before items 2 and 3**, not by weight
(it is the cheapest item on the list) but by influence: the direction of the blade's
torsional stiffness is upstream of both torsion items. Item 2 (#14) asks whether the
`1.61x` twist against Zhou et al. 2025 is structural, and names "compare against a
reference whose model is published" as its third option; S-7 is that comparison and
bounds the structural side. A 6 % stiffness difference cannot produce a 61-125 % twist
difference; a 30-40 % torsional softness is a different conversation. Which bound
applies is decided here, and nothing torsion-shaped is citable from the tree until it
is.

## Evidence as of 2026-10-07

### E1 - the definitions are the same (refutes the issue's likely origin)

Issue #20 hypothesises that `tools/run_s7_torsion.py` measures a different band or
slice count than the test, because it prints its own `ratio_mid`. It does not:

| | test | tool |
| --- | --- | --- |
| band | `(zs >= 0.3 * 117.0) & (zs <= 0.9 * 117.0)` | identical |
| ratio | `theta_s[mid] / theta_beam[mid]` | identical |
| reducer | `float(np.nanmean(ratios))` | `(...).mean()` |
| sections | `section_twists_deg(mesh, u)` -> `n_slices=40` | same call |
| reference | `beam_twist_profile_deg(MOMENT_NM, zs)` | same call |

`nanmean` and `mean` agree unless a station's beam twist is exactly zero, and the
band's upper edge is `105.3 m`, inside the `z <= 111.2` validity cut of
`beam_twist_profile_deg`. The ratio expression `float(np.nanmean(ratios))` has not
changed since the test was created in `5234b48` (2026-09-23): `git log -p --follow`
shows exactly one definition, added once.

=> "A different band or slice count" is not the explanation.

### E2 - the committed table, and its shape

Added by `9d2bcee` (2026-09-30 18:17), verbatim:

```text
#   es=2.000  1460 nodes  ratio=1.386
#   es=1.000  3040 nodes  ratio=1.406   <- the earlier fixture
#   es=0.500  9271 nodes  ratio=1.303
#   es=0.250 32325 nodes  ratio=1.286   <- now
#   es=0.125 120352 nodes ratio=1.273
```

All five are `> 1` (the shell softer than the deck), the committed reading is
"the 1.406 was coarse-mesh discretisation, not the element's torsional response, and
the band itself does not move".

### E3 - the table does not reproduce at HEAD, under either winding convention

Measured 2026-10-07 on HEAD, same mesh generator and same `ELEMENT_SIZE` per row. The
control column is `build_mesh_and_properties` as HEAD produces it; the treatment column
applies the `canonicalize_windings(mesh_model, span_axis=2)` call that `22d3ccc`
(2026-10-04) removed from `BladeMesh.generate`:

| es | nodes | HEAD | HEAD + canonicalised | committed table |
| --- | ---: | ---: | ---: | ---: |
| 2.000 | 1460 | **1.0237** | 1.1530 | 1.386 |
| 1.000 | 3040 | **1.0219** | 1.1640 | 1.406 |
| 0.500 | 9271 | **0.9790** | 1.0763 | 1.303 |
| 0.250 | 32325 | **0.939603** | 1.014223 | 1.286 |

The `es=0.250` control reproduces the live test's `0.939603` exactly. Restoring the
removed canonicalisation flips **1732 of the 33462** elements and moves the ratio
`0.939603 -> 1.014223`, a factor **1.079** - the same *direction* as the table (softer)
but only about a quarter of the way to `1.286` (which would need 1.37).

=> The committed table is **unreproducible at HEAD at every mesh size, with and without
the code it was measured against**. It is a record from a tree state that no longer
exists, i.e. the 2026-09-30 -> 2026-10-04 window, which holds the `origin/main` merge of
2026-09-30, `22d3ccc` (the winding removal), and the nuMAD geometry changes
`680cf81` ("drop the rotorspin override that mirrored the blade in x") and `6b2cbd6`
("give the generated sweep/prebend the input span shape"). The winding removal is a
measured contributor and not the whole cause; the residual is not bisected.

### E4 - the direction is settled: the shell is stiffer

Four independent sources agree, and only the stale table disagrees:

| source | value | direction |
| --- | --- | --- |
| live S-7 at HEAD (`es=0.250`) | `0.939603` -> `1/ratio = 1.064` | **stiffer +6.4 %** |
| S-1 sectional, per `docs/validation_closures.md` G2 | shell `GJ` **+23.9 %** | stiffer |
| `docs/model_parity_audit.md:125,295` | ratio `0.84`, "+17-19 % stiffer" | stiffer |
| pre-merge closure line `docs/validation_closures.md:84` | `shell/viga GJ = 1.080` -> `1/1.080 = 0.926` | stiffer |
| the committed table + `docs/origin_main_integration_2026-09-30.md:156` | `1.273 .. 1.406` | softer |

The S-1 sectional `GJ` and the S-7 global twist measure the same physical quantity two
ways, and they agree on the sign. The pre-merge `1.080` (`GJ` ratio) sits 1.5 % from the
live `1/0.939603 = 1.064` - and the `origin_main_integration` row compares that `1.080`
`GJ` ratio against a `1.406` twist ratio, which are reciprocals.

=> The side to burn is the table's. The live `0.939603` stands, and it also shows that
the *sectional*-versus-*global* tension the closure table carries (S-1's `+23.9 %`
sectional against `+6.4 %` global here, `+7.1 %` at the converged mesh) is real and not a
direction question.

### E5 - a second, new finding: the convergence record is replaced by measured rows

The committed table claimed the `0.5 -> 0.25` step as `1.303 -> 1.286` (**1.3 %**) and
concluded "the band itself does not move". The HEAD record, `es = 0.125` measured to
settle it:

| es | nodes | ratio | step | cost |
| --- | ---: | ---: | ---: | --- |
| 2.000 | 1460 | 1.0237 | | |
| 1.000 | 3040 | 1.0219 | -0.2 % | |
| 0.500 | 9271 | 0.9790 | -4.2 % | |
| 0.250 | 32325 | 0.939603 | -4.0 % | 74 s for both tests |
| 0.125 | 120352 | **0.933641** | **-0.63 %** | 42 s mesh + 230 s solve, **12 GB peak RSS** |

The `step` column is each value against the coarser mesh above it; the last step's two
directions are `0.635 %` below and `0.6385 %` above, so `0.250` sits **0.64 %** above
`0.125`.

=> The fixture **is** converged, so `ELEMENT_SIZE = 0.25` stands and its justification is
now reproducible. The 4 % steps belong to the two coarsest meshes, not to the fixture.
What the old table got wrong is the **value and the direction** (1.286 against 0.934), not
the existence of a plateau. The `12 GB` peak RSS is recorded so the next person does not
run `0.125` casually.

### E6 - the independent verification, and the one claim it corrected

A read-only verifier re-checked every numeric claim against live output and against git at
`73a0736`. All of them hold except one: "within 0.6 %" is false unrounded - `0.939603`
against `0.933641` is `0.6385 %` above the converged value (`0.635 %` in the other
direction) - so the four sites that said `0.6 %` now say `0.64 %`. It independently
confirmed the store state (`check` 0 errors; the S-7 module moving from the only "neither
grouped nor declared" file to none; the reference error count unchanged at 1, and that
error pre-existing), the code facts (`22d3ccc` removes the canonicalisation call, the ratio
expression created once in `5234b48` and never changed, the deck md5), and that no site in
`docs/` or `tests/` still asserts the softer direction.

## Tasks

- [x] T1 Measure `es = 0.125` (120352 nodes) at HEAD and rewrite the convergence record
  with HEAD-measured rows, so `ELEMENT_SIZE = 0.25` is justified by a record that
  reproduces. **Done 2026-10-07**: `0.933641`, `0.64 %` above the `0.250` fixture, and the
  decision is to keep `0.25`; numbers and cost in E5.
- [x] T2 Burn the stale side in the test and land the WIP that was in the working tree.
  **Done 2026-10-07**, commit `4258ba4`: the convergence comment carries the HEAD rows and
  names the window the old one came from, the module docstring states the measured
  direction instead of a derived `0.84`, the comparison uses the suite's
  `assert_relative_error` with `kind="code"`, and `tests/run_step1_revalidate.srm` points
  at the module's current path. Bound untouched. `2 passed`; `ruff check` and
  `ruff format --check` clean.
- [x] T3 Correct the documents that read the stale side. **Done 2026-10-07**, commit
  `9a09772`: the anchor row of `docs/origin_main_integration_2026-09-30.md` and its
  next-steps list are marked superseded (the row had also compared a `GJ` ratio against a
  twist ratio, which are reciprocals); `docs/model_parity_audit.md` no longer attributes
  S-1's sectional `+17-19 %` to S-7 nor quotes the derived `0.84` as the measurement;
  `docs/validation_closures.md` carries the measured record, the converged fixture and the
  re-run trigger, and its G2 S-7 row moves to the re-measured `GJ = 1.071`. The sites that
  already named the stiffer side were re-checked and left as they were.
- [x] T4 Give S-7 its store home. **Done 2026-10-07**, commit `73a0736`: group 34 in
  `docs/validation/groups.yaml` with its `source_files` and provenance note, the rows file
  from `extract --group 34 --write`, the comparison measured against the BeamDyn deck with
  `rtol 0.30` and a written justification (the store demands one above 5 %, and this is the
  widest bound in it), the residual pattern for the group in
  `docs/validation/residual-patterns.json`, and `iea15mw_deck` in `references.yaml`.
  Correction to the issue's own wording: the applied-moment test has **no comparison at
  all** (the extractor reads it as `no_comparison`), so it is declared through the group's
  `non_validation_tests`, not through `non_reference_asserts` - the latter is for a bare
  assertion inside a row that otherwise has one. The citation is `iea15mw_deck` (v1.1,
  Apache-2.0, byte-identical to the deck the test reads) rather than an OpenFAST `r-test`
  key: the deck is vendored here, so its provenance is recorded from the vendored tree.
- [x] T5 Verify. **Done 2026-10-07**: the module green (`2 passed`), `check` 0 errors,
  `status` 0 files neither grouped nor declared, `references check` at its one pre-existing
  error, `ruff` clean, and the independent read-only verifier of E6.
- [x] T6 Close #20 with the measured evidence and update the roadmap comment on #18.
  **Done 2026-10-07** with the maintainer's authorization: #20 is closed as completed with
  the measurement (comment 6049335952) and the map in #18 carries the updated order
  (comment 6049338529). Three follow-ups were opened from it: #21 (the `regression`
  declaration gap), #22 (the unattributed residual) and #23 (the pre-existing
  `references check` error).

## Delivered, with the evidence

| commit | what |
| --- | --- |
| `4258ba4` | the test: the convergence record re-measured at HEAD, the direction stated once, `kind="code"`, the WIP landed |
| `216dbc2` | this feature note, opened with the measurement |
| `9a09772` | the documents that read the stale side, and the closure log's measured record with its re-run trigger |
| `73a0736` | the store home: group 34, the rows file, the residual pattern, the `iea15mw_deck` citation |
| (this note) | the verifier's corrections, the `0.6 %` -> `0.64 %` fix and the refreshed group digest |

## Measurements kept outside the repo

`$SCRATCH/s7_diag/` (diagnostics are not repo artifacts):

- `winding_old.py` - `src/aeroelast/core/mesh/winding.py` as of `22d3ccc^`.
- `diag_winding_ratio.py` - control against treatment at `es = 0.250`.
- `winding_ratio.log` - the raw output of the second, re-run on 2026-10-07 so the
  canonicalised column and the `1732` flip count have a retained source.
- `diag_convergence.py` - the four mesh sizes, both conventions.
- `convergence.log` - the raw output of the second.
- `diag_es0125.py` and `es0125.log` - the `es = 0.125` run, with its wall time and RSS.

## Follow-ups, out of scope here

- The residual between the canonicalised reconstruction (`1.014`) and the table
  (`1.286`) is not attributed. Bisecting the 2026-09-30 -> 2026-10-04 window
  (`680cf81`, `6b2cbd6`, the merge's element line) would close it, but it does not change
  what #20 decides and it is not needed to burn the table. **Now #22.**
- **`regression` does not honour `non_validation_tests`.** It reports a declared test as
  `unclaimed` and counts that in its failing set, while `extract` honours the same
  declaration. Pre-existing and store-wide: `regression --group 30` reports 29 the same
  way. The fix is the one `extract` already applies - route the leftover nodes through
  `classify_unclaimed` instead of comparing against `store.rows` directly. It costs a few
  lines but it is store-tooling work, not #20's, so it is not done here. **Now #21.**
- **`references check` carries one pre-existing error**, `ko2017_nonlinear` at
  `src/aeroelast/core/assembler.py:681` (the declared citation site no longer mentions the
  work). Unrelated to this feature; recorded so the next session does not read it as new.
  **Now #23.**
- The row files under `docs/validation/rows/` are written by the store with one long line
  per scalar, so the prose linter flags all 31 of them; the fields are machine-owned and a
  re-write re-collapses any attempt to reflow them.
- `docs/validation_closures.md:84` (`GJ = 1.080`) is now marked as the invalidated
  pre-merge reading and the G2 row carries the re-measured `1.071`; nothing further owed.
- The S-1-versus-S-7 tension (`+23.9 %` sectional against `+6.4 %` global) is already in
  the closure table and is not opened here.
