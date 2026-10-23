# Blade input audit: mesh vs yaml, field by field

Does the 0.25 m shell mesh carry what `tests/IEA-15-240-RWT.yaml` says? This audits the yaml →
NuMAD → section → element translation, the one cause category that neither L1/L2 nor the CalculiX
parity can exclude (CCX shares the geometry writer and the CLT). Diagnostic only: no fix and no
store row in this feature. Context: `l3-physical-validation.md`, `twist-discrepancy-closure.md` (P1 J1).

## Specs

- S1 (user, verbatim): "en lo que esperamos vamos a hacer la auditoria campo por campo de la malla contra el yaml"
- S2 Independent reference. The yaml is read by a standalone reader (PyYAML + numpy) that imports
  nothing from `aeroelast.models.blade` or its NuMAD layer. Span interpolation is linear on the
  yaml's own grids. The `nd_arc` convention is taken from the WindIO ontology and cited, not
  assumed.
- S3 Stations: r/R = 0.30, 0.60, 0.90 (the P3a stations) and 0.97 (where K66 collapses, P1 J1).
  Mesh: official yaml, `element_size = 0.25`, the production `Blade` path.
- S4 Fields and pre-registered tolerances (gated unless marked):

  | id | field | mesh side | yaml side | tolerance |
  |---|---|---|---|---|
  | F1 | chord, twist, pitch-axis position, prebend, sweep | the ring of outer nodes at the station | `outer_shape_bem` | chord 0.5 %; twist 0.1°; positions 1 cm or 0.5 % c |
  | F2 | airfoil shape | ring polygon | the yaml airfoils interpolated as the yaml defines (`airfoil_position`) | max normal distance 0.5 % c. If the yaml leaves the blending undefined, the scheme is reported and the field is not gated |
  | F3 | laminate stack per ring element: material, ply thickness, fiber angle, order | the element's section properties | the layers active at (r, `nd_arc` of the element centroid) | total and per-ply thickness 1 %; angle 0.1°; material and order exact |
  | F4 | webs: attachment `nd_arc` on both sides, and web stack | web elements at the station | `webs`, `web*` layers | position within es/2 = 0.125 m; stack as F3 |
  | F5 | integrals per unit length: mass Σρt·ds, axial Σ(E₁-projected)t·ds | ring elements | the yaml stack on the yaml geometry | 1 %. Deck M11/K11 printed as context, not gated |
  | F6 | material constants passed to the element (E1, E2, G12, ν12, ρ) | element material | `materials` | 1e-9 relative |
  | F7 | reference surface: where the nodes sit (OML or mid), and the offset used | mesh nodes and element offset | yaml OML plus the stack thickness | reported, not gated |

  An element whose centroid lies within es/2 of a layer's `start_nd_arc`/`end_nd_arc` is a
  "boundary" element. It is reported apart and not gated, because the boundary discretisation is
  ±es/2 by construction.
- S5 Decision, fixed now. Any gated field outside its tolerance on a non-boundary element is an
  **input-translation defect candidate**. It is reported with its location, then estimated for its
  effect on S-1 (EA, EI, GJ at that station) in a later, separately pre-registered step. It is not
  fixed here. If all gated fields pass, the input translation is cleared at these four stations, and
  the S-1/S-2 gaps point to model class or reference quality (P3, C1).
- S6 Outputs: `tools/audit_blade_input.py` and its unit tests; `logs/check/input-audit/report.json`
  and `stations.csv`; the summary table committed here.

## Tasks

- T1 | S1-S4 | parent | done | pre-registration (this file), committed before any tool edit
- T2 | S2-S4, S6 | worker | done | `tools/audit_blade_input.py`, 24 unit tests (RED: collection error; GREEN 24/24), ruff clean
- T3 | S4-S6 | parent | done | login-node run (~25 s), table below; S5: defect candidates A-E, none fixed
- T4 | S5 | - | superseded by T5-T9 (user L2): correct first, measure the effect of each correction
- T5 | S5 | parent | done fb273f2 | D, out-of-plane constants: order fixed; measured null (see Corrections)
- T6 | S5 | parent + explore | done 0c110d3 | A: skin smearing justified (mass-conserving); TE band clamp fixed
- T7 | S5 | parent | f5f860f; base aa62f91 (T7a 49299ca, T7b c9fad3b, T7c aa62f91) | B: bays limited to 0.5 m; S-2 measured; S-1 blocked by the extractor (T7a)
- T7a | S5 | - | pending | `postprocess/sectional.py:325-333`: the tangent-frame local fit uses a window of 2x the median station spacing, which assumes uniform stations; with B's stations one station is alone (`polyfit` degree 0, IndexError). Fix: at least the 3 nearest stations (identical when the window already holds 3+). Then remeasure T7 (`submit_wait.sh T7`, worktree pinned to the fix)
- T7b | S5 | - | decided L5 | default of `max_bay_length`: global in `BladeMesh` (now) or production `Blade` only until the full suite has run -> user L5: global default 0.25 m
- T10 | S5 | worker | committed 8db6c19 (fix/input-T10), measuring | F1: S-1 EA(CLT) metric with hoop stress free, 1/(A^-1)11 instead of A11 (`tools/run_s1_sectional.py:55-86`); verify the elem_areas index (`:218-235`)
- T8 | S5 | worker | G2 committed c5ad881 (fix/input-T8), measuring; G3: option 3 committed c157c62 (fix/input-T8w), measuring | C+E as G2+G3: chord read on its native yaml grid (tip knots); webs placed by the yaml web definition, not the spar-cap edges
- T11 | S5 | worker | committed 7242e92 (fix/input-T11), measuring | G5: LE band clamp at 0.1 side arc (`keypoints.py:216-217`, `:264-265`), same treatment as A
- T12 | S5 | worker | committed ccc0987 (fix/input-T12), measuring | F2: ply stack order so the outermost yaml layer sits on the outer face (`build_rust_properties`, `laminate.rs:77-84`)
- T13 | S5 | explore | done: not established, no change | G1: establish the reference lofting of relative thickness first; change only if the reference is established, else record
- T14 | S5 | worker | committed (fix/input-T14), integrates after T9 | F4: `ElementConfig.span_direction` default None (`core/config.py:661`), fibre along edge 0->1 in FSI configs that omit it; the key is misplaced in `tests/test_kg_prestress_coupling.py:36` and `tests/validation/element/test_prestress_centrifugal_analytic.py:104,214`
- T15 | S5 | worker | committed fix/input-T15 (63760fb G8, b487763 G9, 7453c0c ValueErrors, 4ada981 d from corner); G6/G7 reported; integrates after T9 |
- T16 | S5 | - | deferred to the second mesh-correction iteration (user L8), after T9/T14/T15 | multi-blade `RotorMesh` fibre direction: blades i>0 are rotated by i*2pi/n about y (`generators.py:1569`) but `span_direction` is one global vector (`assembler.py:457`), so on the 3-blade rotor (e.g. `tests/IEA15MW/frontiersin_2025_yaw/solid_corotational_sharpy_yaw_0.yaml`, explicit [0,0,1]) blades 2 and 3 get the fibre ~60 deg off their span (estimate). Pre-existing; found while reviewing T14. Likely fix: per-element radial span direction about the rotor axis (assembler or Rust). G6-G9: planar z sections ignore the prebend tilt; flat-TE half-face cubic with repeated control points (fold at >= 4 elements); web control weights 0.6666/0.3333; `L = ceil(z_tip)`; also the open `d` measured from the flat-back face midpoint and the never-raised `ValueError`s in `yaml_to_blade.py`
- T17 | - | explore | inventory running; separate feature after the second iteration (user L9) | multi-turbine cross-check: inventory the IEAWindSystems RWT repos (windIO version, BeamDyn/6x6 and the tool that produced it, OpenFAST tests, docs), primary goal (user L10): close every possible meshing error, so run the input audit and the mesh checks on each ingestible blade first (clamps and guards the IEA 15 MW never binds); then pick a second blade for S-1, pre-register which hypotheses it discriminates (EA/EI_edge bias systematic or turbine-specific; deck missing TE band mass; G1 lofting; mesher clamps not exercised by IEA 15 MW)
- T18 | - | - | future, separate feature (user L11) | generic windIO importer for any version: detect the version, normalise to one canonical internal form (prefer the official windIO converters/schema if they cover v1->v2; to check), validate against the windIO schema, one reader downstream. Acceptance idea: the upstream IEA-15 windIO v2 file and our v1 file give the same mesh and the same audit. Today `yaml_to_blade.py` reads v1 only. Prerequisite for T17 on v2 turbines
Order fixed by user L5: T7a, T7b, T10, T8, T11, T12, T13, T9, then T14, T15 (user L6: not in the order does not mean not fixed).
Parallel execution (user L6): each candidate in its own worktree/branch off the T7a/T7b base, measured against that base; the parent integrates in order; T9 measures the combination.
- T8 | S5 | parent | pending | C and E, chord at 0.97 and the web attachment at 0.30: diagnose, fix if they are defects
- T9 | S5 | parent | pending | re-audit, every gated field passes or is justified; re-measure the L3 rows; restamp

Rules for T5-T9, fixed now:
- One candidate per commit, in the order D, A, B, C/E: cheapest and clearest first.
- Each fix is test-first: the audit's own gated field for that candidate is the RED test.
- After each fix, S-1 (EA, EI_edge, EI_flap, GJ, mass) and S-2 LC1-LC4 are re-measured at es 0.25 on
  SLURM. The effect of each candidate is then a measured difference, not an inference.
- No band and no anchor is touched. A row that moves is a finding to record.
- A candidate that turns out to be a deliberate NuMAD modelling choice is recorded with its
  justification instead of being "fixed".

## Result (2026-10-23, revision of the tool commit; `logs/check/input-audit/report.json`, `stations.csv`, 1568 checks)

Counts are pass/fail/boundary/not-gated.

| r/R | F1 shape | F2 airfoil | F3 stack | F4 webs | F5 integrals | F6 materials | F7 surface |
|---|---|---|---|---|---|---|---|
| 0.30 | 5/0/0/0 | 0/0/0/1 | 294/72/106/52 | 20/1/0/0 | 2/0/0/6 | 25/0/0/25 | 0/0/0/3 |
| 0.60 | 5/0/0/0 | 0/0/0/1 | 199/17/104/36 | 21/0/0/0 | 0/2/0/6 | 25/0/0/25 | 0/0/0/3 |
| 0.90 | 5/0/0/0 | 1/0/0/0 | 80/44/68/24 | 17/4/0/0 | 0/2/0/6 | 15/0/0/15 | 0/0/0/3 |
| 0.97 | 4/1/0/0 | 1/0/0/0 | 96/8/68/24 | 1/0/0/0 | 0/2/0/6 | 15/0/0/15 | 0/0/0/3 |

F5, mesh vs the yaml point value:

| r/R | mass | axial |
|---|---|---|
| 0.30 | -0.10 % | -0.48 % |
| 0.60 | -5.15 % | -2.28 % |
| 0.90 | -6.86 % | -10.2 % |
| 0.97 | +2.73 % | +2.03 % |

Fields that pass:
- Outer shape: twist within 0.035°, positions within tolerance.
- Airfoil: where gated, within 0.29 % c.
- In-plane materials: exact.
- Fibre angles: all 0 in the yaml and in the mesh.
- Material and ply order: they differ only on boundary elements.
- Nodes: they sit on the OML with the laminate centred on them (F7, reported only).

**S5 verdict: input-translation defect candidates exist.** None was fixed. They are listed with the
worker's numbers, and only D was spot-checked by the parent:
- **A. Partly covered layers are smeared over the whole NuMAD region** (`stackdb._map_continuous_layers`):
  - `Shell_skin` in `TE_REINF` is -15.5 % at 0.30 and -13.3 % at 0.60;
  - `Shell_skin_inner` in `TE_PANEL` is -9 % at 0.30, -23 % at 0.90 and -30 % at 0.97;
  - TE-reinforcement glass mass at 0.60 is 47.1 vs 61.0 kg/m (-23 %), which is the whole F5 mass miss there.

  This bears on EI_edge, which S-1 measures as -7.5 % signed (softer).
- **B. Thickness is held per NuMAD bay**, not evaluated at the element. At 0.90 the spar caps are 10.18 vs
  11.43 mm (-11.0 %). Bay-averaging the yaml brings F5 at 0.90 to -0.9 % mass and -0.6 % axial, so
  B accounts for the 0.90 miss. Mesh refinement cannot remove it: the bay comes from the NuMAD
  stations, not from `element_size`. Rule 9 convergence studies therefore do not see it.
- **C. Chord at 0.97 is +0.77 %** (1.8897 vs 1.8752 m). Pchip versus linear explains part of it, and
  0.5 % is unexplained.
- **D. Out-of-plane constants** (`build_rust_properties`, model.py:721):
  - the Rust order is (g12, g23, g13; nu12, nu23, nu31), and windIO (xy, xz, yz) is passed in;
  - glass_uni g13/g23 come out swapped (0.57 %);
  - nu23 is 0.35 vs 0.264, and 0.0 vs 0.17 for glass_triax.

  How the shell uses these constants (transverse shear) has not been checked.
- **E. The web1 suction-side attachment at 0.30 is 0.157 m off** (tolerance 0.125 m).

Literal-application notes, each recorded in the report meta:
- The ring is the planar cut at z.
- The TE flat-back elements are excluded from F3/F5 (0.001-1.6 % of the ring mass).
- r/R is the blade-length fraction.
- The `nd_arc` convention is inferred from the airfoil order and the `side` tags. The windIO schema
  text is not vendored, so it is not cited, which departs from S2.
- The position tolerance is the larger of the two values.

**Reading.** "Implementation is clean" no longer holds for the input chain. A and B reach
the stiffness of the TE and the spar caps, the members that set EI_edge and EA. Whether they explain
S-1/S-2 is a quantitative question for T4. The answer is not assumed here.

## Corrections (T5-T7, 2026-10-23)

Measurement: `logs/check/input-fix/measure.sbatch <TAG>` runs S-1, S-2 (official) and V-03 at es 0.25 from a
detached worktree pinned to the candidate commit (`$SCRATCH/wt/input-fix-<TAG>`, with a PYTHONPATH shim that
exposes only `aeroelast` so the site-packages `_aeroelast` stays in use). `compare.py` applies the test
metrics (S-1 mid-span mean |rel| and signed; S-2 tip vs the BeamDyn p60 constants of the S-2 test). The
baseline is the 9f72c6b capture `logs/check/restamp-9f72c6b/capture-s1s2v3.pkl`; `src/` and `crates/` did
not change between 9f72c6b and 96901a3. Each row is cumulative (T6 includes T5).

| metric | 9f72c6b | T5 fb273f2 | T6 0c110d3 | T7 f5f860f |
|---|---|---|---|---|
| S-1 EA mean / signed % | 7.594 / +7.569 | 7.594 / +7.569 | 8.030 / +8.029 | error (see below) |
| S-1 EI_edge mean / signed % | 8.821 / -7.506 | 8.821 / -7.506 | 7.663 / -3.766 | error |
| S-1 EI_flap mean / signed % | 5.697 / +3.786 | 5.697 / +3.786 | 5.727 / +3.766 | error |
| S-1 GJ mean / signed % | 5.051 / -2.759 | 5.051 / -2.759 | 6.454 / -1.366 | error |
| S-1 mass mean / signed % | 2.235 / +0.727 | 2.235 / +0.727 | 3.194 / +1.740 | error |
| S-2 LC1 / LC2 / LC3 / LC4 vs BD % | -7.44 / -3.47 / +15.57 / -6.27 | same | -6.44 / -3.49 / +12.00 / -6.28 | -6.54 / -3.26 / +12.16 / -6.30 |
| V-03 mass kg | 66682.3 | 66682.3 | 66903.4 | 66902.8 |

**D (T5, fb273f2).** Confirmed: windIO/NuMAD give `G = [G12, G13, G23]` and `nu = [nu12, nu13, nu23]`;
`OrthotropicMaterial` (Rust and Python) takes `(g12, g23, g13)` and `(nu12, nu23, nu31)`. Both
`build_rust_properties` and `material_factory` passed the windIO order through. Rust reads g13/g23 only in
the transverse shear (`orthotropic.rs:116-118`, c44/c55) and never reads nu23/nu31. RED: the new
`tests/test_blade_material_constant_order.py` (g13 got G23). After the fix the audit's not-gated extras match
the yaml (glass_uni g13 3.48, g23 3.50 GPa, nu23 0.264; glass_triax nu23 0.17). The gated summary is identical.
**Measured null.** S-1, S-2 and V-03 agree to 4 decimals. At full precision GJ moves 1.7e-5, EI 8e-7,
S-2 at most 4e-7, and EA and mass by exactly 0. This confirms the fix is active and that a thin shell does not
see a 0.57 % G13/G23 swap.

**A (T6, 0c110d3).** The explore subagent's diagnosis, checked by the parent with
`logs/check/input-fix/t6/kp_dump.py`, splits A in two:
- **Skin smearing (`Shell_skin` in TE_REINF, `Shell_skin_inner` in TE_PANEL): a deliberate NuMAD choice, kept.**
  - The arc-overlap fraction times the region arc equals the yaml overlap arc, so mass is conserved per
    component. It only moves a strip of at most 0.05-0.2 m inward. The `stackdb.py:147-157` docstring states
    this is deliberate.
  - The fraction normaliser excludes the flat-back face, as the audit does; the perimeters agree to 0.1 %.
  - The -9/-23/-30 % F3 rows at 0.30/0.90/0.97 are this redistribution, seen element by element.
- **TE reinforcement band: a defect, fixed.**
  - Keypoint `d` was capped at 85 % of each side's arc (`keypoints.py:225`, `:272`), so the band was at most
    15 % of the side. The cap exists to keep `d` aft of `c` (`c` is clamped at 80 %).
  - Measured, as the band from the TE reference point: 0.59 m at r/R 0.60 and 0.43 m at 0.80, against 0.75 m
    in the yaml. The clip runs from r/R 0.47 to 0.82. 0.59/0.75 × 0.977 = 0.77, the audit's 47.1/61.0 kg/m.
  - The TE panel (c to d) is at least 0.8 m wide there, so the cap protected nothing.
  - Fix: `d` yields only to stay 5 % of the side arc aft of `c`. When `c` sits at its 80 % clamp this equals
    the old cap.
  - RED: `tests/test_blade_te_band_width.py` (19 of 45 stations clipped, at most 0.315 m).
  - Audit at 0.60: F5 goes from 0/2 to 2/0 (mass and axial within 1 %); F3 failures go from 17 to 4, the
    worst now total thickness -1.7 %, which is B. 0.30, 0.90 and 0.97 are unchanged.
- **Left open (not fixed):** `d` is measured from the flat-back face midpoint (`arclength[0]`), so the band
  starts h/2 short of the corner inboard, 0.732 vs 0.75 m at 0.60. That part goes to the TE_FLAT stacks:
  the mass moves rather than disappears. The yaml does not say whether the flat face carries the
  reinforcement.

**B (T7, f5f860f).** The stations are the `outer_shape_bem` z grid (50 points, bays of about 2.4 m), and each
(region, bay) holds one bay-mean laminate. `element_size` cannot reach it.
- Fix:
  - `BladeMesh(max_bay_length=0.5)` inserts interpolated stations through the new
    `Blade.add_interpolated_stations` (one regeneration); `add_interpolated_station` now delegates to it.
  - The value is fixed in metres, so it is a laminate resolution independent of `element_size`.
  - RED: `tests/test_blade_bay_length.py` (the default was None).
- Span study with the audit (`logs/check/input-fix/t7/`):

  | bay limit | stations added | F5 mass 0.90 | F3 fails 0.90 | chord 0.97 | audit wall |
  |---|---|---|---|---|---|
  | native ~2.4 m | 0 | -6.86 % | 44 | +0.77 % | ~25 s |
  | 0.5 m | 196 | -0.84 % | 8 | +0.70 % | 40 s |
  | 0.25 m (= es) | 441 | -1.73 % | 20 | +0.70 % | 66 s |

- **Finding: not monotone.** At bay = es the 0.90 integrals get worse. A probable cause is that element rows
  and bays fall out of alignment. This was not diagnosed. 0.5 m is the only value the evidence supports.

**Prediction check (A and B add the material the mesh was missing, so EA, already +7.6 % stiff, should get
worse).**

T6 (job 11623114) confirms it for A: EA goes from +7.57 to +8.03 %. It also brings findings to record:
- **EI_edge improves sharply**, from -7.51 to -3.77 % signed and 8.82 to 7.66 % mean. S-2 LC3 (edge) follows,
  from +15.57 to +12.00 %. The clipped TE glass was half of the edgewise softness.
- **Mass and GJ get worse in the mean.**
  - Mass goes from +0.73 to +1.74 % signed, and from 2.23 to 3.19 % mean; V-03 rises by 221 kg.
  - GJ goes from 5.05 to 6.45 % mean, although the signed GJ improves from -2.76 to -1.37 %, so the error is
    now spread both ways along the span.
  - The mesh moved toward the yaml (audit F5 at 0.60) and away from the BeamDyn deck mass. That suggests the
    deck itself may not carry the full 0.75 m TE band. Not checked.
- **The flap rows do not move**: EI_flap, LC2 and LC4 change by at most 0.03 points. LC1 improves from -7.44 to
  -6.44 %, through the added gravity load.
- No band was touched. The S-1 EA, GJ and mass rows moving outward are a finding, not a regression to
  undo.

**T7 measurement (job 11623121).** S-2 LC1 / LC2 / LC3 / LC4 moved from T6 by -0.10 / +0.23 / +0.16 / -0.02
points, V-03 mass by -0.6 kg. B barely moves the global deflections, as expected for a bay-mean
(mass-conserving) laminate. All five S-1 tests errored at setup (T7a), so EA and EI after B are unmeasured.

**The 0.25 m bay result was a sampling artefact, not a mesh effect (geometry audit G4, inferred, not run).**
Element rows align with the bays by construction (`mesh_gen.py:300-330`). The audit's `_slice` moves the
cut to z + 1e-9 when a node lies on the plane (`tools/audit_blade_input.py:451`), so it samples the outboard
bay. At a 0.25 m limit, r/R 0.90 falls exactly on a station; at 0.5 m it falls at a bay mid-point. Check
before revisiting the 0.5 m choice: rerun the audit at 0.25 m with cuts at bay mid-points.

## NuMAD mesher theoretical audit (2026-10-23, read-only, two explore subagents, estimates only)

User L4: NuMAD targets ANSYS/Abaqus layered sections; the shell adaptation for aeroelast is the user's. Full
lists are in Engram `fem-shell/blade-input-audit/numad-theory-audit`. The new candidates, ranked by expected
effect, are not yet in the plan:

| id | finding | evidence | expected effect |
|---|---|---|---|
| F1 | S-1 EA(CLT) metric uses the hoop-constrained A11, not 1/(A^-1)11 | `tools/run_s1_sectional.py:55-86` | EA overstated by ~1-3 % (part of +8 %); metric only |
| G1 | relative thickness linear between airfoil positions; the reference lofting may be pchip | `yaml_to_blade.py` `_add_stations` | +/-4 % thickness at r/R 0.05-0.25, EI_flap +/-8 % locally |
| G2 | chord tip knots 0.985/0.99/0.995 dropped (resampled to the 50-pt grid) | `yaml_to_blade.py` `_add_stations`, `geometry.py:225` | explains C fully; tip chord -7..-35 % |
| G3 | webs always on the spar-cap edges, yaml web `offset_y_pa` ignored | `mesh_gen.py:461-463`, `:557-559` | explains ~0.10 m of E; GJ <~1-2 % |
| G5 | LE band clamp at 0.1 side arc clips the 0.8 m LE reinforcement, r/R ~0.59-0.80 | `keypoints.py:216-217`, `:264-265` | same class as A; EI_edge -0.3..-0.5 % |
| F2 | ply stack inside-out: ply 0 (UV/gelcoat) at -h/2 with outward normals, so it sits on the inner face (likely the ANSYS layered-section convention) | `stackdb.py`, `laminate.rs:77-84` | EI/GJ +0.1-0.4 %; ply stress on the wrong face |
| F4 | `ElementConfig.span_direction` defaults to None: FSI configs without it get the fibre along edge 0->1; two tests set the key at the top level, where it is ignored | `core/config.py:661`; `tests/test_kg_prestress_coupling.py:36`; `tests/validation/element/test_prestress_centrifugal_analytic.py:104,214` | none on S-1/S-2; production FSI risk |
| G6-G9, F3 | prebend tilt of planar sections, flat-TE cubic overshoot, web weights 0.9999, `ceil(z_tip)`, fibre projection under prebend | Engram | <0.3 % or latent |

Found correct: spar-cap placement (`swtwisted=0` matches the windIO `rotation: twist`), fibre frame (projected
global z, `assembler.rs:257-335`), Q/Qbar, thickness units, mass integration, web stacks, shear correction.

## Second batch (2026-10-23, user L5/L6: everything this session, in parallel)

**T7c, the root cause behind T7a (aa62f91).** `SectionalExtractor.station_membership` read the station as
`name[3:5]`. Past station 99 the NuMAD set name has three digits, so `08_147_LP_SPAR` counted as station 14.
S-1 merged stations 14 and 140-149: the T7 IndexError, the B0 `x must be increasing` (job 11623128), and EI
computed without error on mixed stations. The fix parses the whole field (`_station_of`). On the es 0.25
mesh this gives 490 stations, z increasing, dz 0.234-0.244 m. The T7a window widening stays (harmless,
robust) but was a symptom. Also: the first B0 (11623128) ran in the main repo, because its submission
raced the worktree creation; the sbatch now aborts when `cd $W` fails or `aeroelast` is not from the shim.

**Parallel candidates.** Each is one commit on its own branch off aa62f91 and is measured against B0 (aa62f91).
The parent integrates them in L5 order in T9.
- T10 F1 (8db6c19): `laminate_a11_eff = 1/(A^-1)[0,0]` feeds `ea_gj_from_static_response`. Analytic UD
  ply: 0.971 x the old value (= E1 h exactly). `elem_areas` id vs position was checked: ids are contiguous,
  so not a bug, only fragile.
- T8 G2 (c5ad881): `Definition.chord_curve` holds the yaml chord on its own grid. Its knots become
  stations, and inserted stations interpolate it linearly. `ipercentthick` is pchip of `percentthick`.
  Audit F1 at 0.97: 1 fail (+0.70 %) -> 5 pass. Only chord had the grid loss. Left open:
  `definition.chord` (used by `aerodynamics.py`) is still on the coarse grid, so the BEM chord differs from
  the structural chord (added to T15's follow-ups).
- T11 G5 (7242e92): keypoint `a` yields only to stay 5 % of the side arc forward of `b`. Before, the clip ran
  r/R 0.59-0.80 at 0.289-0.393 m per side against the yaml's 0.4 m (up to 28 % of the LE glass). Audit F3/F5
  change only at 0.60 (F5 axial -2.19e-3 -> -2.18e-3).
- T12 F2 (ccc0987): skin stacks are reversed (inner -> outer). Normals are 13506/13506 outward by
  point-in-polygon; the audit's 92-96 % was its centroid heuristic. Stress recovery resolves by z and
  follows. `tools/audit_blade_input.py` now reads stacks outer face first; its F3 counts equal the base
  (182/2/108 at 0.60). Ply-stress parity: 20 passed.
- T13 G1: **not established, no change.** WISDEM and windIO v2 (`outer_shape.rthick`) use pchip. SONATA,
  which per the IEA-15 v1.1.4 release notes produced the BeamDyn deck S-1 uses, blends linearly in span.
  Rule 2 forbids picking the scheme that closes a number. Linear minus pchip t/c: -6.5 % at s 0.06,
  +4.5 % at s 0.20, < 1 % past s 0.25. Revisit when the deck's SONATA input path is pinned.
- T8 G3 webs: **stopped, needs the user.** Every yaml web attachment lies inside a spar cap (for example
  web1 SS at 0.30: nd 0.2878, cap 0.2794-0.3558). The mesh places webs on the cap edges. Only web1 SS at
  0.30 fails (-0.157 m against 0.125); three more rows sit at 0.113-0.122 m. The 0.157 m splits into
  ~0.06 m (mesh cap edge `b` vs yaml cap `start_nd_arc`) and ~0.10 m (web inside the cap). Options:
  (1) document as an approximation; (2) derive the cap edges from the yaml cap nd_arc (~1 day, removes
  the 0.06 m, shifts every cap region); (3) split the cap regions at the yaml web arcs (~2-3 days,
  faithful, 8 extra regions per bay, new set naming).

**T15 (worker, four commits on fix/input-T15, integrated after T9):**
- **G8:** web control weights exact (`_web_control_points`).
- **G9:** the span ends at the yaml tip z. The official yaml is bit-identical, since its tip is exactly 117.0.
- **ValueErrors:** the four errors for reinforcements without an ss/ps pair now raise; the LE ones no longer say "te reinf".
- **d from the corner:** the TE band is measured from the flat-back corner, not the face midpoint, so flat-back bands grow by the
  half-face height, at most 0.06 m. `tests/test_blade_te_band_width.py` moved its reference point accordingly, with the
  tolerance unchanged.
- **Audit:** every changed row moves toward the yaml (for example, F3 at 0.60: -6.07e-3 -> -4.14e-3), and no row gets worse.
- **G7 reported only, does not bind:** the largest half-face is 0.062 m, so there is one element across it. None of the 1898
  TE_FLAT elements is degenerate or inverted.
- **G6 reported only:** prebend tilt is 1.6 deg at mid-span and 5.6 deg at the tip; cos phi is 0.48 % at most.
  - Options: (a) leave as is; (b) re-enable the NuMAD `prebend_rot`, which moves every node and the FSI surface;
    (c) a cos phi scale.
  - Candidate for the second iteration.
- **Untraced side warning:** `element_utils.py:78` (`rotVec = (theta/mag)*cp`, invalid divide).

**Measurements of the second batch** (S-1, S-2 official, V-03 at es 0.25, sequana_cpu; `logs/check/input-fix/capture-<TAG>.pkl`).
Signed %, each candidate against B0 = aa62f91 (T6 + bays 0.25 m + the S-1 extractor fixes):

| metric | T6 0c110d3 | B0 aa62f91 | T10 F1 | T8 G2 | T11 G5 | T12 F2 |
|---|---|---|---|---|---|---|
| S-1 EA | +8.03 | +8.10 | **+6.28** | +8.10 | +8.13 | +8.10 |
| S-1 EI_edge | -3.77 | -2.92 | -2.92 | -2.93 | -2.85 | -2.92 |
| S-1 EI_flap | +3.77 | +5.55 | +5.55 | +5.53 | +5.57 | **+5.12** |
| S-1 GJ | -1.37 | -2.16 | -2.16 | -2.23 | -2.12 | -2.01 |
| S-1 mass | +1.74 | +1.62 | +1.62 | +1.62 | +1.66 | +1.62 |
| S-2 LC1 / LC2 | -6.44 / -3.49 | -6.45 / -2.82 | same as B0 | -6.26 / -3.30 | -6.41 / -2.83 | -6.19 / **-2.18** |
| S-2 LC3 / LC4 | +12.00 / -6.28 | +11.75 / -6.25 | same as B0 | +11.21 / -6.17 | +11.68 / -6.26 | +11.68 / -5.89 |
| V-03 mass kg | 66903.4 | 66904.1 | 66904.1 | 66914.5 | 66910.2 | 66904.1 |

Readings:
- **F1:** EA drops 1.8 points from the metric alone; the S-2 rows match B0, which confirms it is metric-only.
- **G2:** the tip chord moves S-2 (LC2 -0.48, LC3 -0.53 points); S-1 is about unchanged.
- **G5:** small, as expected for a narrow band; it adds material, so EA moves outward by +0.03.
- **F2:** the predicted sign holds. EI_flap -0.44 points (softer), and S-2 flap moves toward BeamDyn (LC2 +0.64, LC4
  +0.36). The magnitude is a little above the 0.1-0.4 % estimate. EA and mass are unchanged, as they must be. GJ +0.15
  points is not interpreted.
- **T6 -> B0 is confounded, not attributable (finding):**
  - The bays change the spanwise element rows: rows are `ceil(bay/size)`, about 0.12 m at mid-span instead of about 0.18.
    S-2 LC2 moved +0.23 points from 2.4 to 0.5 m bays and +0.44 from 0.5 to 0.25 m.
  - The S-1 extractor depends on station density: the EI local-fit window is twice the median station spacing (4.8
    -> 0.48 m), and `UnivariateSpline(s=1e-1)` saw 490 instead of 50 points.
  - Candidate-vs-B0 deltas are valid, because they share one extractor. B0 against 9f72c6b is not, so T7d comes before T9.
- T7d | S5 | worker | pending | make S-1 extraction independent of station density (reproduce the 9f72c6b method on any
  station set), and measure the spanwise row size and element count at 9f72c6b-geometry vs B0.

**T7d (133edc5, fix/input-T7d).** S-1 now fits on the mesh stations nearest to the native bay centres (51,
built from the yaml with the generator's tip refinement). On the native mesh this is the identity; on the 0.25 m mesh
each selected station is within 0.126 m of its native centre.
- **Measured (same mesh as B0; S-2 and V-03 identical, as they must be).** Signed %, B0 -> T7d:
  - EA +8.10 -> +8.21
  - EI_edge -2.92 -> -3.42
  - EI_flap +5.55 -> +4.85
  - GJ -2.16 -> -1.67
  - mass +1.62 -> +1.72
- **So the extractor artefact was 0.5-0.7 points on EI/GJ.** T6 -> T7d is the real effect of the 0.25 m bays
  (EI_edge +0.35, EI_flap +1.08, GJ -0.30 points). It still mixes the bay-mean laminate and the finer spanwise mesh;
  T7e separates them.

**Finding: the bay limit couples to the spanwise mesh.** NuMAD meshes each bay with `ceil(bay/size)` rows, so 0.25 m
bays halve the spanwise element size:

| | native | bays <= 0.25 m |
|---|---|---|
| nodes / elements | 32,515 / 33,673 | 47,894 / 49,610 (+47 %) |
| span size r/R 0.1 / 0.3 / 0.6 / 0.9 | 0.239 / 0.217 / 0.171 / 0.133 m | 0.119 everywhere |

So "es 0.25" no longer means the Rule 9 mesh. Part of T6 -> B0 in S-2 can be discretisation, and the Rule 9
convergence (0.25 vs 0.125 m) was done with native bays. User L12: decouple fully (T7e).
- T7e | S5 | worker | done: 7ef630e + f28d0ae on fix/input-T7e (off T7d); reviewed and measured (session 3) | spanwise rows from element_size and its grading only (native discretisation); laminate evaluated per element row, replacing the max_bay_length coupling

**T7e review (session 3, 2026-10-23).** The worker had committed before the close; the Handoff 2 "uncommitted" was
stale.
- **7ef630e, finding: inserted stations lost the root TE band.** `KeyPoints.generate` gives every station before the
  first nonzero band the first nonzero width; the linear interpolation of the yaml 0 -> 750 mm step put a tiny width at
  the first inserted station, so the root TE reinforcement came out 0.08 m wide instead of 0.67 m. Any dense station set
  triggers it, so **T7 (0.25 m bays) and T8w (built on T7) carry it**; B0/T7d numbers include it (the worker's
  commit message; not re-verified in session 3). Fix:
  `Blade._interpolate_band` keeps the leading zero run.
- **Tests.** 11/11 in the worktree, imported through the T7e shim (`test_blade_span_decoupling`, `_bay_length`,
  `_te_band_width`, `_material_constant_order`). Counts equal native: 32,515 / 33,673; row z to < 1e-6 m.
- **Geometry (not pinned by the tests; `logs/check/input-fix/t7e-geom.py`, `t7e-surf.py`, es 0.25).**

  | decoupled vs native | p50 | p95 | p99 | max |
  |---|---|---|---|---|
  | nearest node | 0.65 mm | 2.5 cm | 4.1 cm | 15.4 cm (z 25-28) |
  | distance to the other surface | 0.08 mm | 3.2 mm | 8.0 mm | 2.4 cm (z 3.3-3.6, root circle); 2.9 cm reverse (z 25-28, x -3.8) |

  The nodes mostly slide along the same surface (chordwise partition changes). The residual shape change, up to
  ~3 cm (~0.5 % of the local chord), comes from the inserted stations being `Blade`-interpolated sections where the
  native mesh uses the mesher's in-bay patch. The yaml defines nothing between stations, so neither is the reference:
  it is an interpolation choice, recorded as a finding. The f28d0ae message ("the mesh equals the native one") holds
  for counts and row z only.
- **T9 merge map (`git merge-tree`, no writes).** One textual conflict: T7e x T8w in `mesh_gen.py`. T7e x T10,
  T8 x T8w, T8w x T11, T8w x T12 and T7e x `l3/twist-p1` merge cleanly. The conflict is also semantic: T8w loops
  cap `pieces` inside region `j`, and the piece count varies per bay at the root (no B piece for eta <= 0.082), so the
  `SpanPlan` key `("shell", j)` must include the piece, and every inserted row must carry the pieces of its native
  bay, or the counts stop being native. T8w's 49,890 / 51,631 were on 0.25 m bays; after the merge expect native
  counts plus the piece splits.
- **Measured (SLURM 11623208, `submit_wait.sh T7e`, es 0.25, `compare.py capture-T6 capture-T7d capture-T7e`).**
  Same 6 failing tests as T7d and B0 (4 S-1 bands, S-2 LC1/LC4 official), none new. Signed %:

  | metric | T6 (native mesh, yaml-bay laminate) | T7d (0.119 m span mesh, 0.25 m bays) | T7e (native mesh, row laminate) |
  |---|---|---|---|
  | S-1 EA | +8.03 | +8.21 | +8.03 |
  | S-1 EI_edge | -3.77 | -3.42 | -4.18 |
  | S-1 EI_flap | +3.77 | +4.85 | +4.61 |
  | S-1 GJ | -1.37 | -1.67 | -1.26 |
  | S-1 mass | +1.74 | +1.72 | +1.58 |
  | S-2 LC1 gravity | -6.44 | -6.45 | -6.56 |
  | S-2 LC2 tip flap | -3.49 | -2.82 | -3.42 |
  | S-2 LC3 tip edge | +12.00 | +11.75 | +12.11 |
  | S-2 LC4 distributed | -6.28 | -6.25 | -6.33 |
  | V-03 mass [kg] | 66,903.4 | 66,904.1 | 66,907.9 |

  - **T6 -> T7e is the laminate resolution alone** (same discretisation; T6 has no inserted stations, so the TE band
    defect does not apply to it; the inserted-section geometry, <= 3 cm, is in it). EA 0.00, EI_edge -0.42,
    EI_flap +0.84, GJ +0.11, mass -0.16 points; S-2 moves <= 0.12 points on all four LCs.
  - **L3 prediction "A and B remove material, so EA worsens": not for B.** EA is unchanged (8.029 -> 8.027) and mass
    drops only 0.16 points.
  - **Finding: the T7d S-2 LC2 shift (-3.49 -> -2.82) was discretisation, not laminate.** It reverts to -3.42 on the
    native mesh. A 0.6-point tip-flap move from 0.25 -> 0.119 m spanwise is the order of the Rule 9 convergence
    (0.84 % between 0.25 and 0.125), so the T7 0.25 m bays default silently changed the mesh class of every L3 number;
    T7e restores the Rule 9 mesh.
  - Status: T7e done; enters T9 first after T7d. No band, anchor or contract.yaml touched.

**G3 (c157c62, fix/input-T8w).**
- **Geometry.** Webs attach at the yaml `start/end_nd_arc`. The caps are cut into HP/LP_SPAR_A/B/C pieces that keep the
  parent laminate and share nodes with the web edges.
- **Audit F4.** Residuals are <= 4 mm (before: up to 0.157 m); F4 passes 21/21 at 0.30, 0.60 and 0.90.
- **Mesh.** 49,890 nodes / 51,631 elements against B0's 47,894 / 49,610. The 132 triangles include collapsed wedges
  at stations 40 and 479, where two webs meet the cap.
- **Two findings along the way:**
  - The old pairing of webs to yaml webs was by name. It was wrong, and harmless only because both layups are identical.
  - The mesh HP side is the yaml pressure side: the HP arc is 1 - end_nd_arc.
- **Measured vs B0 (job via submit_wait T8w), signed %:**
  - GJ -2.16 -> **-11.04**, distributed over r/R 0.2-0.85: median -8 %, range -17 to -4 %. A singular zone at r/R
    0.085-0.105 swings -92 % / +171 %.
  - EI_flap -0.67 points; S-2 flap about 0.5 % softer; LC3 +0.4 points.
  - EA and mass are unchanged.
- **Not accepted as physics:** a 0.1 m web shift should move GJ by 1-2 %.
  - Topology is intact: 12 interior free edges, all at web tip ends, against 38 in B0; no element area under 1e-5 m^2.
  - In the root bay the two webs share one line.
  - A verify subagent is testing four hypotheses: H1 physical (independent multi-cell Bredt), H2 S-1 twist extraction
    (which nodes), H3 property and fibre-angle assignment of the cap pieces, H4 element quality. G3 does not enter T9
    until this is resolved.
- **Open risk.** Thin cap pieces (5-11 mm wide) near r/R 0.1. The Rust triangle path for the wedges is exercised
  only by the SLURM run.

**T17 inventory (explore, 2026-10-23, web read only, not verified where marked):**
- **All on windIO 2.x now.** Every current master yaml in IEAWindSystems is windIO 2.x; our `tests/IEA-15-240-RWT.yaml` is
  v1. v1 copies exist only in the windIO repo `test/turbine/v1p0/` (15, 22) or in pre-update commits (3.4 `824f8d6`,
  10 `350cf06`). The official converter goes v1 -> v2.1 only, so T18's canonical form should be v2.
- **v2.1 yamls embed the 6x6** (`elastic_properties`), and every v2 file has a spanwise `rthick` (relevant to G1).

| RWT | 6x6 source | webs | caps | flatback | use |
|---|---|---|---|---|---|
| 3.4-130 | SONATA (same chain as 15); ANBA mass structural only (issue #5) | 2 | glass | root (issue #2) | 2nd blade for S-1 |
| 10-198 | converted from HAWC2 st (BECAS, unconfirmed) | 3 (TE web 0.05-0.475) | glass | no | mesher (3 webs, ~30 DP layers); cheapest v1 ingest |
| 22-280 | BECAS/AESOpt, independent | 3 (TE web 0.05-0.5) | carbon, 1.0 m | FB90-60 | 3rd blade: mesher stress test, Campbell in repo; issue #176: TE laminate overlap past r/R 0.69 (input defect) |
| LB 6.3/6.9/7.2 | likely SONATA | 2 | carbon | root | v2 only, no report |

- **Link to the deck-mass hypothesis (T6):** SONATA/ANBA mass is structural only (3.4 #5), and the IEA-15 deck shares
  that chain. Check whether the IEA-15 deck mass excludes some yaml layers before blaming the TE band.
- **Mesher gaps exposed:** 3 webs, a partial-span TE web, and v2 `anchors`/`plane_intersection` placement are not
  supported by the current mesher or reader. These are T18/T17 scope.

## Handoff 2 (2026-10-23, session close, user: "cerramos aqui")

**Resume**:
1. `mem_search "blade-input-audit"` and read this file up to `## Log`.
2. Check the two agent results appended below under "Agent results at close". If they are absent, the agents did not
   finish: rerun them from the T7e and G3 paragraphs.

**Branch state.** Everything is local; nothing is pushed.
- Main branch: `l3/twist-p1` = base aa62f91 (T5 D, T6 A, T7 B, T7a/b/c), plus docs.
- Candidate branches, each one commit unless noted, worktrees in `$SCRATCH/wt/input-fix-<TAG>`:

| branch | commit | content | measured vs B0 |
|---|---|---|---|
| fix/input-T7d | 133edc5 | S-1 fits on native bay centres | yes: extractor artefact 0.5-0.7 points |
| fix/input-T7e | 7ef630e, f28d0ae (session 3: was committed, not uncommitted) | root TE band at inserted stations; decouple: native rows, laminate per row (user L12) | yes (session 3): vs T6 laminate-only <= 0.84 pt S-1, <= 0.12 pt S-2 |
| fix/input-T10 | 8db6c19 | F1 EA metric | yes: EA -1.8 points |
| fix/input-T8 | c5ad881 | G2 chord grid | yes |
| fix/input-T8w | c157c62 | G3 webs, caps split | yes: GJ -8.9 points, **physical** (verify, below); enters T9 per the rules |
| fix/input-T11 | 7242e92 | G5 LE band | yes |
| fix/input-T12 | ccc0987 | F2 ply order | yes |
| fix/input-T14 | 61b0d52 | F4 span_direction default | no (after T9) |
| fix/input-T15 | 63760fb..4ada981 (4 commits) | G8, G9, ValueErrors, d from the corner | no (after T9) |

**Next, in order:**
1. **T7e.** Review and commit the worker diff in its worktree. Check its geometry (node displacement vs native) and its
   counts (must equal native: 32,515 / 33,673). Measure with
   `bash logs/check/input-fix/submit_wait.sh T7e`; the comparison is against T7d.
2. **G3.** The verify verdict is physical (see "Agent results at close"). G3 enters T9 under rule "a moved row is a
   finding". Open question for the next session: why the yaml-faithful shell is about 11 % below the BeamDyn GJ.
3. **T9.**
   - Merge into `l3/twist-p1` in this order: T7d, T7e, T10, T8, T8w, T11, T12. Expect conflicts between T7e and T8w in
     `mesh_gen.py`/`generators.py`.
   - Run `check.sh quick` and the blade tests.
   - Re-audit: `tools/audit_blade_input.py`.
   - Measure the combination against 9f72c6b (same extractor as T7d).
   - Restamp the L3 rows, and register the moved anchors as findings.
4. **Merge T14 and T15**, then measure.
5. **Second iteration (user L8).** T16 (multi-blade fibre direction), G6 (prebend tilt), the untraced
   `element_utils.py:78` divide warning, `definition.chord` on the coarse grid for the BEM.
6. **Future features:**
   - T18: generic windIO importer, canonical form v2. It is a prerequisite for T17.
   - T17: multi-turbine cross-check. Robustness of the mesher first (IEA-22, IEA-10); then S-1 on IEA-3.4.
7. **Before any merge to main**, run the full suite on SLURM; crates are untouched.

**Rules that still hold:**
- One candidate per commit.
- RED first.
- No band, anchor or contract.yaml change.
- A moved row is a finding.
- Measure on SLURM with `logs/check/input-fix/measure.sbatch` via `submit_wait.sh <TAG>`:
  - it needs worktree `$SCRATCH/wt/input-fix-<TAG>` plus the shim `<TAG>.pp/aeroelast` symlink;
  - it picks sequana_cpu or sequana_cpu_bigmem, 50 min.

**Cleanup when done:** `git worktree remove` for input-fix-T5/T6/T7/B0 and the merged branches;
`$SCRATCH/tmp/input-fix-*`.

**Agent results at close:** (appended below when they arrive)

**T7e at close:** the worker (task mvlioo6b) was still running when the user closed the session. Its diff is
uncommitted in `$SCRATCH/wt/input-fix-T7e`. Next session:
- check `git -C $SCRATCH/wt/input-fix-T7e status` and `git diff --stat`;
- run its tests (`tests/test_blade_span_decoupling.py`, `tests/test_blade_bay_length.py`);
- check the node/element counts against native 32,515 / 33,673 and the geometry displacement;
- if incomplete, relaunch from the T7e task line (design: span_plan from a native planning pass, one row per bay,
  native chordwise counts, `BladeMesh(laminate_per_row=True)`).

**G3 GJ drop, verify verdict: H1 physical** (`logs/check/input-fix/t8w-diag/`, captures reproduced to 0.0).
- **An independent multi-cell Bredt-Batho reproduces the drop.** Wall `1/(A^-1)66` from each layup, 3 cells:

  | r/R | Bredt T8w/B0 | FE T8w/B0 |
  |---|---|---|
  | 0.3 | 0.934 | 0.9365 |
  | 0.6 | 0.8975 | 0.8999 |

  The FE-Bredt absolute offset (+3 % / -3.3 %) is the same on both meshes.
- **The whole drop is the webs' share.** Skin-only J is unchanged. The webs add +23 % (0.3R) and +39 % (0.6R) to GJ in
  B0; in T8w that share drops about 35 %, because the webs move toward each other (endpoint extent 0.897 -> 0.773 m at
  0.3R, 0.896 -> 0.723 m at 0.6R). **The parent's premise ("0.1 m moves GJ 1-2 %") was wrong.**
- **H2 extraction, ruled out.** Twist fitted on skin-only nodes gives -11.02 % against -11.04 %, and the FE twist itself
  differs: theta(0.85R) 0.690 vs 0.775 rad.
  - The r/R 0.085-0.105 swing is extractor noise. The twist there is ~2e-4 rad against the tolerance of the
    `UnivariateSpline(s=1e-3)` at `tools/run_s1_sectional.py:187`; Bredt gives 0.98-1.00 there. It lies outside the
    0.2-0.85 test window.
  - Optional fix: theta-scaled `s` or a mask below r/R 0.13.
- **H3 properties and H4 element quality, ruled out.**
  - All 2,840 SPAR_A/B/C layups equal their parent's.
  - Normals are 100 % outward.
  - Cap aspect ratio: median 1.93 vs 1.84, maximum 14.3.
- **Root zone.** The yaml webs coincide for eta <= 0.082, so 50 bays have no B piece and the LP C piece narrows to
  8.5 mm at bay 50. This has no GJ effect.
- **Consequence (finding).** With the yaml web positions the shell GJ is about 11 % below BeamDyn (B0: -2 %). B0's
  agreement was a compensation: the webs sat too wide.
  - The remaining gap needs its own hypotheses (inference, untested): web-to-cap overlap at the mid-surface, cap
    thickness idealisation, or the SONATA web modelling.
  - No band was touched.

## Handoff (2026-10-23, first close at ~190 k context; superseded by the second batch above)

Resume: `mem_search "blade-input-audit"`, then read this file up to `## Log`. Next actions, in order:
1. T7a: fix the extractor and remeasure T7. Record S-1 after B, and check the EA prediction for B.
2. T7b: the user decides the default.
3. The user decides which audit candidates join the plan. The parent's recommendation: F1 first (it is a
   metric fix and changes how the EA row reads), then G2+G3 as T8 (they are C and E), then G5, then F2.
   G1 needs its reference established first. F4 is a separate workstream.
4. T9.

Tools: `logs/check/input-fix/{measure.sbatch, compare.py, submit_wait.sh}`. Worktrees live in
`$SCRATCH/wt/input-fix-T5/T6/T7`; remove them with `git worktree remove` when done.

Delivery: the branch is past the ~400-line slice since origin/main (3198 lines, 11 commits). Each
candidate is its own commit, so the review slices are fb273f2, 0c110d3 and f5f860f. Push/PR are out of scope
for this session.

## Log

- L1 (2026-10-23) User: "en lo que esperamos vamos a hacer la auditoria campo por campo de la malla contra el yaml". This follows the four-way split agreed in the same session: implementation, input translation, model class, reference quality. Input translation goes first because it is the cheapest discriminator.
- L2 (2026-10-23) User: "y para corregir lo que encontraste en la auditoria de malla? yo creo que so seria lo principal para comenzar, porque todo lo demas depende de la malla". Decision: correct A-E before T4, C1 or P3, one at a time with before/after measurement (T5-T9).
- L3 (2026-10-23) User, resuming: "Trabajo, en este orden y con estas reglas (ya fijadas en el documento): T5 D → T6 A → T7 B → T8 C/E → T9 re-auditoría. Un candidato por commit. Primero diagnosticá: si es una elección deliberada de NuMAD, registrá la justificación en vez de 'corregir'. […] Predicción a verificar: A y B quitan material, así que EA debería empeorar. Si no empeora, es un hallazgo. No toques bandas, anclas ni contract.yaml."
- L4 (2026-10-23) User: "que tal si en lo que esperamos haces una auditoria profunda teorica del mallador de numad, porque es un mayador que puede estar dise~nado para otros propositos y tal vez tenga errores que podamos arreglar o mejorar" / "la verdad numad esta pensado para ansys y abaqus donde utiliza otro tipo de elementos de pegar capas (no recuerdo bien) yo hice la adaptacion a shell de nuestro aeroelast"
- L5 (2026-10-23) User: "vamos a poner el por defecto en 0.25 e implementar todo en esta misma sesion en el orden de prioridades recomendado"
- L6 (2026-10-23) User: "haz tareas en paralelo para no sobrecargar al orquestador" / "que no este en el orden recomendado, no quiere decir que no hay que corregirlo, asi que no lo olvides dejalo guardado en memoria en alun lugar"
- L7 (2026-10-23) User chose G3 option 3: "Opción 3: partir el cap en la posición del alma".
- L8 (2026-10-23) User: "creo que todavia quedan cosas por corregir de la malla como el F4, entonces vamos a notar ese T16 para cuando terminemos lo que estamos haciendo ahora, en la segunda iteraccion de correcion del mallado"
- L9 (2026-10-23) User: "en la cuenta de github de nrel hay varios repositorios que presentan el modelo de una tubina diferente con su fichero de windio y sus tests de openfast [...] tal vez tener mas geometrias y ensallos podria respadar o refutar algunas hipotesis comparando comportamientos [...] https://github.com/IEAWindSystems ten en cuenta que no todos estan tan completos como el de la iea15mw" / on the inventory: "si por favor"
- L10 (2026-10-23) User: "la idea es cerrar todo lo que fuera posibles herrores de mallado"
- L11 (2026-10-23) User: "no sabia que habian varias versiones de windio, por lo que hay que planificar para un futuro un importador generico de windio para cualquier version"
- L12 (2026-10-23) User chose "Desacoplar del todo (recomendado)": element rows from element_size and its grading as before; laminate evaluated at each element row's z, not per bay.
- L13 (2026-10-23, session 3) User, on the plan "revisar el diff de 7ef630e y f28d0ae; correr tests y chequear conteos/geometría; mandar submit_wait.sh T7e; registrar 7ef630e": "si".
