# Feature: S-7's torsion ratio says both directions at once (#20)

Status: in progress (T1-T5 open, measurement done 2026-10-07)
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
the *sectional*-versus-*global* tension the closure table already carried (`+23.9 %`
against `+8.0 %`) is real and not a direction question.

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
| 0.125 | 120352 | **0.933641** | **-0.6 %** | 42 s mesh + 230 s solve, **12 GB peak RSS** |

=> The fixture **is** converged: `0.250` is within **0.6 %** of `0.125`, so
`ELEMENT_SIZE = 0.25` stands and its justification is now reproducible. The 4 % steps
belong to the two coarsest meshes, not to the fixture. What the old table got wrong is
the **value and the direction** (1.286 against 0.934), not the existence of a plateau.
The `12 GB` peak RSS is recorded so the next person does not run `0.125` casually.

## Tasks

- [x] T1 Measure `es = 0.125` (120352 nodes) at HEAD and rewrite the convergence record
  with HEAD-measured rows, so `ELEMENT_SIZE = 0.25` is justified by a record that
  reproduces. **Done 2026-10-07**: `0.933641`, `0.6 %` from the `0.250` fixture, and the
  decision is to keep `0.25`; numbers and cost in E5.
- [x] T2 Burn the stale side in the test and land the WIP that was in the working tree.
  **Done 2026-10-07**, commit `4258ba4`: the convergence comment carries the HEAD rows and
  names the window the old one came from, the module docstring states the measured
  direction instead of a derived `0.84`, the comparison uses the suite's
  `assert_relative_error` with `kind="code"`, and `tests/run_step1_revalidate.srm` points
  at the module's current path. Bound untouched. `2 passed` in 74 s; `ruff check` and
  `ruff format --check` clean.
- [ ] T3 Correct the documents that read the stale side:
  `docs/origin_main_integration_2026-09-30.md:156` (stale value **and** a `GJ` ratio
  compared against a twist ratio) and any other site quoting the table. Sites that quote
  the *stiffer* side (`docs/model_parity_audit.md:125,295`,
  `docs/validation_closures.md:84`, `docs/shell_vs_beam_sectional_validation.md:167`) are
  re-checked, not rewritten.
- [ ] T4 Give S-7 its store home, which it does not have today (`status` reports the
  module as the only "neither grouped nor declared" file):
  a group entry in `docs/validation/groups.yaml` with `source_files`, the rows file
  generated by `extract --group <id> --write`, the `ratio` comparison measured against
  the BeamDyn deck with `rtol 0.30`, the applied-moment assertion declared through the
  group's `non_reference_asserts` (the mechanism exists: `groups.yaml` validates the key
  and group 31 declares a bare assertion the same way), and the BeamDyn citation resolved
  in `docs/validation/references.yaml` (`gaertner2020` is the definition report;
  the `K[5,5]` tables come from the OpenFAST `r-test` deck, so the reference may need its
  own key).
- [ ] T5 Verify: the S-7 module green at the chosen `ELEMENT_SIZE`;
  `tools/validation_matrix.py check` with 0 errors; `status` reporting 0 files "neither
  grouped nor declared"; `ruff check` and `ruff format --check` on every touched Python
  file; both YAML files parse. An independent read-only verifier re-checks every numeric
  claim in the store text and the rewritten comments against the live output.
- [ ] T6 Close #20 with the measured evidence and update the roadmap comment on #18
  (item 8 done, and item 2 unblocked with the bound it now has). Publishing and closing
  are the maintainer's decision.

## Measurements kept outside the repo

`$SCRATCH/s7_diag/` (diagnostics are not repo artifacts):

- `winding_old.py` - `src/aeroelast/core/mesh/winding.py` as of `22d3ccc^`.
- `diag_winding_ratio.py` - control against treatment at `es = 0.250`.
- `diag_convergence.py` - the four mesh sizes, both conventions.
- `convergence.log` - the raw output of the second.

## Follow-ups, out of scope here

- The residual between the canonicalised reconstruction (`1.014`) and the table
  (`1.286`) is not attributed. Bisecting the 2026-09-30 -> 2026-10-04 window
  (`680cf81`, `6b2cbd6`, the merge's element line) would close it, but it does not change
  what #20 decides and it is not needed to burn the table.
- `docs/validation_closures.md:84` (`GJ = 1.080`) sits inside the G1/G2 block the
  2026-09-30 merge invalidated; T2 decides whether it is re-anchored or left as a
  historical record.
- The S-1-versus-S-7 tension (`+23.9 %` sectional against `+8.0 %` global) is already in
  the closure table and is not opened here.
