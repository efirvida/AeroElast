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

> **Reversed by E7 (2026-10-10).** The direction below is an artifact of the tip couple's
> normalisation, not a measurement: with the load honoured the shell is `~3.8 %` softer globally. The
> stale table's side was right about the direction. Read E7 before using anything in this section.

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

> **Superseded by E7 (2026-10-10)** for the ratios in the table below, every one of which is a
> stiffness ratio times that mesh's `kappa`. The *plateau* result stands; the values do not.

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

### E7 - the direction was an artifact: the tip couple under-delivered its moment (2026-10-10)

> **Direction reversed again by E8 (2026-10-11).** The `kappa` correction below stands, but the
> beam reference it was compared against is a staircase that under-states the beam twist; against
> the integral to `z` the shell is stiffer, on S-1's side. Read E8 before using the direction here.

**This reverses E4.** `run_torsion_case` built its tip couple with
`denom = 2 * sum(x_off**2 for x_off > 0)`, which gives `sum(x_off * F_y) = torque` only for an
x-symmetric ring. The IEA 15 MW tip ring has 16 nodes and is not symmetric, so the applied moment was
`kappa * torque`, with `kappa` a function of the ring geometry alone:

| es | nodes | `kappa` | recorded ratio | ratio / `kappa` | re-measured with the fix |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2.000 | 1460 | 0.915365 | 1.0237 | 1.118351 | **1.118319** |
| 1.000 | 3040 | 0.915365 | 1.0219 | 1.116385 | **1.116392** |
| 0.500 | 9271 | 0.915365 | 0.9790 | 1.069518 | **1.069529** |
| 0.250 | 32325 | 0.897777 | 0.939603 | 1.046588 | **1.046588** |
| 0.125 | 120352 | 0.898586 | 0.933641 | **1.039011** | (derived, see below) |

`kappa` was measured two ways -- `sum(x_off * F_y)` over the tip ring, and the moment read by the
instrument at a station 35 m inboard -- and dividing by it is exact: the response is linear in the
applied moment. The analytic correction predicts the four re-measured rows to 5-6 digits, which is
what licenses the `0.125` entry without re-running its 230 s / 12 GB solve. (Cross-check that the
defect moved the *recorded* numbers and not the solver: T8 had already shown `StaticLinearSolver` and
an independent `_assemble` + `spsolve` give the same band mean to six decimals.)

Two consequences:

1. **The trend never crosses 1.0.** The recorded series crossed between `0.500` and `0.250`, and that
   crossing is the entire "the shell is stiffer" reading. Corrected, the series is monotone,
   `1.118 -> 1.039`, and the shell's global `GJ` is **~3.8 % BELOW** the beam's at convergence.
2. **E4's argument loses its strongest leg.** Three of E4's four "stiffer" sources --
   `1/0.939603`, the pre-merge `GJ = 1.080`, and `model_parity_audit`'s `0.84` -- all derive from
   this same ratio. Only S-1's *sectional* `+23.9 %` is independent, and it measures a different
   quantity (St-Venant section stiffness, not a global tip twist). E4 did not have four independent
   sources agreeing; it had one measurement and three restatements of it, against the table.

**What this does *not* settle.** Sectional and global now disagree in sign: S-1's sectional shell
`GJ` is `+17-24 %` above the deck, S-7's global twist says the shell is `~3.8 %` softer. Those are
different quantities, and a shell can legitimately be sectionally stiff and globally soft (shear lag,
warping restraint, the clamped root) -- but nothing here demonstrates that it is. That is the open
question, and it is exactly the question #20 believed it was closing.

**Burn as settled**: the direction claim in the fixture docstring, in E4 above, and in the S-7 rows.
**Keep**: `_RATIO_TOL = 0.30` -- this is a load-construction defect, not a physics bound, and the
bound was never the defect. **Note**: the older committed table (`1.386 .. 1.273`) had the *right
direction* and still does not reproduce at HEAD at any mesh size; E3's factor `1.079` from the
element-winding canonicalisation does not close it, and neither does `kappa` (which moves those
numbers the other way, to `1.51 .. 1.42`).

### E8 - the sign disagreement was a second comparison error; the stiffness excess tracks the OML reference surface (2026-10-11)

Measured at `a9d3454`, sdumont login node, `scripts/aeroenv.sh`, official blade
`tests/IEA-15-240-RWT.yaml` (S-1's fixture resolves to the same file). Nothing in the solver, the
S-7 module, row 34 or the S-1 docstring was changed. The three scripts below reproduce every number.

**1. S-7's beam reference is a staircase (comparison error).** `beam_twist_profile_deg` integrates
`1/GJ` only up to the last BeamDyn station at or below each query `z` and drops the partial segment.
BeamDyn stations are about 5.85 m apart in the band, and `GJ` falls fast. So the reference under-states
the beam twist: over 0.3-0.9 R it gives **0.84 of the integral to `z` on average, and 0.63 at worst**,
and the shell/beam ratio comes out inflated. `tools/diagnose_s7_beam_reference.py`, same mesh, same
load, same twist estimator:

| es | nodes | staircase (published) | integral to z, K66 | integral to z, 1/C66 |
| --- | ---: | ---: | ---: | ---: |
| 1.000 | 3040 | 1.1164 | 0.9227 | 0.9110 |
| 0.500 | 9271 | 1.0695 | 0.8874 | 0.8763 |
| 0.250 | 32325 | 1.0466 | **0.8670** | 0.8561 |

The staircase column reproduces E7's re-measured rows to four digits, so the reference is the only
thing that moved. Dividing E7's derived `0.125` row by the same correction gives about `0.861`. The
trend is monotone and still moving about 2 % per halving.

Corrected, the shell's global `GJ` is **about 12-15 % ABOVE the deck**, not 3.8 % below it, and on
the same side as S-1. The spread depends on how BeamDyn interpolates between key points: +15 %
against a reference linear in `1/GJ`, +12 % against one linear in `GJ` (es = 0.25: 0.8670 and
0.8920). The BeamDyn user guide is where that gets settled. This is the second time S-7's direction reverses, and the cause is independent of
E7's `kappa`.

Uncertainty carried by the reference itself:
- Interpolating `GJ` linearly instead of `1/GJ` changes the beam twist by `0.972`, and log-linear by
  `0.986`. Which convention BeamDyn uses between key points is **not verified**. Under linear-`GJ`
  the es = 0.5 ratio is `0.913`.
- `K66` against the effective `1/C66` (`C = K^-1`, which frees the shear-torsion coupling of the
  6x6): `1/C66` sits 0.5-5 % below `K66` over 0.3-0.95 R, and bending-torsion coupling is zero.

**2. Sectional and global agree once both are honest.** I re-ran S-1's own torsion path
(`ea_gj_from_static_response`: local-tangent couple, station estimator, spline `s = 1e-3`) at its own
mesh (es = 0.5, default generator):
- the signed mid-span GJ is **+10.8 %** above `K66`;
- S-1's metric, mean `|rel|`, gives **12.7 %**, not the documented **23.9 %**: an anchor that moved,
  reported here and not touched. This was confirmed with the production function itself:
  `ea_gj_from_static_response` with the test's mask gives 12.73 % (+10.76 % signed, 32 stations).
  `z_static` and `z_eg` are both built from `station_membership`, so they are the same grid;
- the cumulative twist of that same solve gives a band ratio of 0.897.

A factorial over mesh (S-1 generator against S-7 generator) x load (local-tangent couple against
global-z `F_y` couple) x estimator (stations against 40 z-slices), all against S-7's staircase
reference, spans only **1.0695-1.0931**. None of the three instrument differences carries the sign.
The reference did.

The couple applied by S-1 is exact: `M . t = T` and the net force is zero, by construction.

**3. Where the excess plausibly comes from: the skin's reference surface.** Two premises, each
checked:
- **The skin nodes lie on the OML.** At eta 0.3 / 0.5 / 0.7 the largest node-to-node extent of the
  skin ring equals the yaml chord to within `3e-4` (5.3739 against 5.3724 m, 4.1545 against 4.1548 m,
  3.2255 against 3.2255 m).
- **The laminate is centred on the nodes.** `Laminate::new` sets `z` from `-h/2` to `+h/2`. Half of
  every skin laminate therefore sits outside the aerofoil, and the torsion cell encloses the OML area.

The tools that produce beam 6x6 tables from a layup conventionally stack plies inward from the OML.
For **this** deck that is **not verified**: the tool that generated `elastic_properties_mb` and its
layer convention remain to be read from the IEA 15 MW report or the yaml's provenance.

Prediction and counterfactual:
- **Bredt prediction** (`tools/diagnose_bredt_oml_area.py`): `A_mid = A_OML - 1/2 * contour integral
  of t ds`, giving a single-cell factor `(A_mid/A_OML)^2` of **0.864** over 169 rings at es = 0.5
  (0.85-0.90 across the band) and 0.876 at es = 1.0.
- **Counterfactual** (`tools/diagnose_oml_reference_surface.py`): every skin node moves in its section
  plane, inward by the mean half-thickness of its skin elements (mean 21.8 mm, maximum 54.5 mm at
  es = 0.5). The solver, laminates and loads are unchanged.

| es | case | twist ratio, 1/GJ-linear ref | twist ratio, GJ-linear ref | flap at 0.9 L [m/N] | edge at 0.9 L [m/N] |
| --- | --- | ---: | ---: | ---: | ---: |
| beam | 6x6 rotated by ElastoDyn `StrcTwst` | 1 | 1 | 4.62e-5 to 4.74e-5 | 1.95e-5 to 1.99e-5 |
| 1.0 | OML (as built) | 0.9227 | 0.9491 | 4.198e-5 | 2.098e-5 |
| 1.0 | mid-surface | 1.0210 | 1.0500 | 4.666e-5 | 2.177e-5 |
| 0.5 | OML (as built) | 0.8874 | 0.9130 | 3.864e-5 | 2.138e-5 |
| 0.5 | mid-surface | 0.9944 | 1.0228 | 4.448e-5 | 2.211e-5 |
| 0.25 | OML (as built) | 0.8670 | 0.8920 | 3.788e-5 | 2.144e-5 |
| 0.25 | mid-surface | **0.9732** | 1.0010 | 4.415e-5 | 2.216e-5 |

Geometry audit of the shift at es = 0.5 (printed by the script):
- 175 nodes in `HP_TE_REINF`/`HP_TE_PANEL` disagree with the centroid vote. That is the concave
  aft pressure side, where the centroid test is unreliable. A point-in-polygon test against each
  original OML ring finds **4 of 7005** shifted ring nodes outside it.
- **90 `TE_FLAT` elements invert**, 42 of them in the band, because the HP and LP shifts cross in
  the thin trailing-edge band.
- The `--keep-te-flat` variant leaves those 1050 nodes on the OML. It gives **0 outside and 0
  inverted**, with GJ factor 0.8934 against 0.8924, identical flap, and edge 2.172e-5 against
  2.211e-5.

The table therefore stands, and the edge shift is the only entry the TE band moves (+1.6 % rather
than +3.4 %).

The inset changes the shell's global GJ by a factor of 0.891-0.904 at all three meshes, which
accounts for **the bulk** of the excess. The mid-surface ratio still falls about 2-3 % per halving
(1.021, 0.994, 0.973). At convergence it would sit a few percent below 1 against the
`1/GJ`-linear reference, and near 1 against the `GJ`-linear one. That residual, roughly 0-5 % too
stiff, is **unattributed**: its size is inside the stacked uncertainty of the reference (±3 %), the
first-order counterfactual, and the unconverged trend.

The same change moves the static flap compliance at 0.9 L from 17-19 % too stiff to 4-6 % too
stiff, so it closes flap as well. Edge goes from 7.5-10 % too soft to 11-13.5 % too soft, so it
opens edge. Both beam compliances are Euler-Bernoulli from the 6x6 with no shear deformation. The
tip compliance cannot serve as a reference: the degenerate `eta = 1.0` station changes it by a
factor of 2 between interpolations.

**What this does not settle.**
- The dynamic flap is already in tension at HEAD. V-02's 1st flap is -3.8 % and S-4's is +0.6 %,
  while the static flap above is about 18 % too stiff. Moving the reference surface would lower the
  flap frequency by about 6 %. Mass, the modal set-up, or the frequency references are therefore
  part of the question, not only stiffness.
- Edge moves the wrong way.
- The reference convention of the deck is unverified, so the OML account is still a hypothesis.

**What follows for the store.** The fix (a reference-surface offset in the laminate, or nodes moved
to the mid-surface) would live in `crates/aeroelast-core/src/materials/**` or the mesh generator.
Neither appears in L3's `depends_on`, which is exactly #33 A8. A8 has to land first, or the fix will
expire no verdict.

**A new item, not absorbed here (A8-adjacent).** `crates/aeroelast-py/src/assembler.rs` computes the
span-angle offset only for 4-node composites (`node_ids.len() == 4`). Composite MITC3 triangles get
offset 0, so their fibres follow their own edge 01. The blade has **43 of them at es = 0.5 and 88 at
es = 0.25** (about 0.3-0.4 % of the composite elements). Also unverified: that the offset's `e1`
(edge 01 projected on the averaged normal) is the frame in which the MITC4 element consumes its
constitutive matrix.

Proposed for the maintainer, nothing applied:
- (a) fix S-7's reference to integrate to `z`, re-measure row 34, and re-state #20's direction;
- (b) open the reference-surface question as a modelling issue on L3, with this counterfactual as its
  entry evidence;
- (c) record that S-1's 23.9 % no longer reproduces.

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
- [x] T7 Fix the couple normalisation and re-derive the direction. **Done 2026-10-10**, commit
  `72f64de`: `tip_couple()` extracted as the single construction, `denom = sum(x_off**2)`, the
  applied-moment test tightened from a `+-30 %` window to exactness (`1e-12` on the forces that reach
  the solver, plus zero net force) after observing RED at `102.22 N.m` of `1000` missing, and the
  module docstring, convergence record and assertion comment corrected to the measured direction;
  E7 records the reversal and the sign disagreement it leaves open. `2 passed`; `ruff check` clean on
  both changed files. **Two follow-ups belong to the maintainer**: store row 34 still carries the old
  direction in `measured`/`notes`/`justification`, and #20's closure rests on E4.

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
