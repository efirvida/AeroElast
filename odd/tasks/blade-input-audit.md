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
- T7 | S5 | parent | open f5f860f | B: bays limited to 0.5 m; S-2 measured; S-1 blocked by the extractor (T7a)
- T7a | S5 | - | pending | `postprocess/sectional.py:325-333`: the tangent-frame local fit uses a window of 2x the median station spacing, which assumes uniform stations; with B's stations one station is alone (`polyfit` degree 0, IndexError). Fix: at least the 3 nearest stations (identical when the window already holds 3+). Then remeasure T7 (`submit_wait.sh T7`, worktree pinned to the fix)
- T7b | S5 | - | pending, needs user | default of `max_bay_length`: global in `BladeMesh` (now) or production `Blade` only until the full suite has run
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

## Handoff (2026-10-23, session closed at ~190 k context)

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
