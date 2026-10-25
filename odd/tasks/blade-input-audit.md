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
- T14 | S5 | worker | merged ad5f3ac (quick OK, 17/17); measured with T15 | F4: `ElementConfig.span_direction` default None (`core/config.py:661`), fibre along edge 0->1 in FSI configs that omit it; the key is misplaced in `tests/test_kg_prestress_coupling.py:36` and `tests/validation/element/test_prestress_centrifugal_analytic.py:104,214`
- T15 | S5 | worker | merged d5d0fd2 (web blocks: T8w corners + G8 thirds; 62 passed). Measured, T15m against T14m (G10b is null): EI_edge -3.71 -> -1.89, EA +6.24 -> +6.56, GJ -10.34 -> -9.78, mass +1.65 -> +2.26 % (V-03 +211 kg), LC1 -5.55 -> -5.26, LC3 +11.26 -> +10.26, same 5 failing tests. Audit F5 (skin rings) is unchanged at 0.30/0.60, so the mass sits outside the skin rings or comes from moved stations; attributed (`logs/check/t15-mass/REPORT.md`): all +211 kg from 4ada981 (keypoint d from the corner), all skin, TE_REINF +269 kg and TE_PANEL -58 kg at z 20-70 m; G8, G9 and the ValueErrors are mass-null. **A fidelity gain:** the TE reinforcement goes from -5.4/-7.8 % to -1.0/-3.4 % against the yaml (HP/LP), and the total from -0.79 % to -0.47 % against the yaml (67,470 kg). The BeamDyn mass row worsens only because the short band had cancelled a surplus; was: committed fix/input-T15 (63760fb G8, b487763 G9, 7453c0c ValueErrors, 4ada981 d from corner); G6/G7 reported; integrates after T9 |
- T16 | S5 | - | deferred to the second mesh-correction iteration (user L8), after T9/T14/T15 | multi-blade `RotorMesh` fibre direction: blades i>0 are rotated by i*2pi/n about y (`generators.py:1569`) but `span_direction` is one global vector (`assembler.py:457`), so on the 3-blade rotor (e.g. `tests/IEA15MW/frontiersin_2025_yaw/solid_corotational_sharpy_yaw_0.yaml`, explicit [0,0,1]) blades 2 and 3 get the fibre ~60 deg off their span (estimate). Pre-existing; found while reviewing T14. Likely fix: per-element radial span direction about the rotor axis (assembler or Rust). G6-G9: planar z sections ignore the prebend tilt; flat-TE half-face cubic with repeated control points (fold at >= 4 elements); web control weights 0.6666/0.3333; `L = ceil(z_tip)`; also the open `d` measured from the flat-back face midpoint and the never-raised `ValueError`s in `yaml_to_blade.py`
- T17 | - | explore | inventory running; separate feature after the second iteration (user L9) | multi-turbine cross-check: inventory the IEAWindSystems RWT repos (windIO version, BeamDyn/6x6 and the tool that produced it, OpenFAST tests, docs), primary goal (user L10): close every possible meshing error, so run the input audit and the mesh checks on each ingestible blade first (clamps and guards the IEA 15 MW never binds); then pick a second blade for S-1, pre-register which hypotheses it discriminates (EA/EI_edge bias systematic or turbine-specific; deck missing TE band mass; G1 lofting; mesher clamps not exercised by IEA 15 MW)
- T18 | - | - | future, separate feature (user L11) | generic windIO importer for any version: detect the version, normalise to one canonical internal form (prefer the official windIO converters/schema if they cover v1->v2; to check), validate against the windIO schema, one reader downstream. Acceptance idea: the upstream IEA-15 windIO v2 file and our v1 file give the same mesh and the same audit. Today `yaml_to_blade.py` reads v1 only. Prerequisite for T17 on v2 turbines
- T19 | S5 | parent | done 554fefd (fix/input-G10), merged df4c210; measured vs T9: S-1 <= 0.014 pt, S-2 LC2 -0.10 pt, LC1 +0.04, LC3 -0.05, mass -2 kg (SLURM, `capture-G10.pkl`) | uniform linear span interpolation of layer thickness (`yaml_to_blade.py` imethod and the `update_internal_structure` resample); user L16
- T20 | - | - | planned, not implemented, on demand (user L16) | `span_scheme="pchip"` option, one scheme for every layer, for T17 turbines whose reference comes from WISDEM >= PR #419. Also open: `stackdb.py:196` pchip on nd_arc_start along the span; `excel_to_blade.py` pchip (UTD xlsx path, disowned input)
- T23 | S5 | worker | done 028d717, merged a680e56 (RED/GREEN; 10/728 stacks change, ply count <= 0.56); measured vs G10: null, every S-1/S-2 row <= 0.003 pt, mass +0.4 kg, same 5 failing tests | G10b: layer `start/end_nd_arc` (and any other layer/web bound) linear along the span; `stackdb.py:196-198` was pchip (Shell_skin, Shell_skin_inner). Outer-shape geometry (`geometry.py`, windIO PCHIP for airfoils) and the `mesh_gen.py` surface loft stay pchip; user L17
- T24 | S5 | worker | done: (a) defect candidate, fix deferred to T28; (b), (c) yaml-faithful | Baseline reproduced at d8a1c3a (67,150.0 / skin 64,648.8 / webs 2,500.6 kg). (a) TE_FLAT glass_uni 366.7 kg is a **translation defect** inherited from NuMAD: `yaml_to_blade.py:449-460` gives the TE band the extent ["d","te"], `component.py:98` maps it to te-e (TE_FLAT) plus e-d, `stackdb.py:91,107` adds the ply to both, while the yaml band is 0.75 m from the TE corner and excludes the flat face. Impact +0.55 % mass (glass_uni +4.3 -> -2.2 % without it), EI_edge +3-7 % at z 15-60 (ring estimate), EA +0.5-1.1 %. Candidate fix: extent `e` (or no glass_uni on `*_TE_FLAT`); it lowers EI_edge (now -1.9 % vs BD, S-4 1E -5.6 %), so it must be judged on S-4, V-02 and S-2, not on mass. (b) triax -215 / carbon -214 kg and (c) webs -4.83 % are the **mid-surface inset** (`skin_reference="mid"`, `skin_reference.py:8`) against a yaml integral on outer arcs: on the OML mesh webs +0.30 %, carbon -0.13 %, triax +0.81 %; density/units, arc truncation and span interpolation ruled out. Open: an OML root surplus +265 kg at z 0-10 (+1.1..+2.4 % at z 3-10), check per-bay thickness sampling against the linear taper there. Note: the audit's F4/F5 passes describe the OML mesh, not the mid mesh the anchors use. `logs/check/t24/FINDING.md`
- T21 | S5 | worker | finding, not merged (user L20: "si no, registrarlo como hallazgo") | CalculiX `tip_axial_tension` parity with the load-work metric W = 1/2 sum f.u (SLURM 11623281-83, read by worker): T7e 9.56 %, T8w 10.63 %, T14m 20.53 % vs bound 5 % (old tip node-mean gap 4.64 / 4.63 / 7.07 %). FAIL on every tree; the work metric does not rescue the band, so the gap is not a node-count artefact: it is distributed along the span (W gap > tip gap). The gap doubles at T14 (span_direction default, fibre along the span), so the next lead is a frame/fibre-orientation mismatch between the AE and CCX models of the same mesh (AGENTS rule 4: frames before physics), not the tolerance. `fix/ccx-parity-work` (worktree input-fix-CCXW) keeps its uncommitted metric diff (unit tests 2 passed, ruff clean); 11623284 full parity timed out, not resubmitted (the gating case already fails). `$SCRATCH/wt/input-fix-CCXW/logs/check/t21/FINDING.md` **Superseded cause (T31 Part 3):** the gap is tip load introduction, not a frame or stiffness mismatch; the `tip_*` W metric is replaced, not widened.
- T22 | S5 | parent | done 091e17d | `test_innermost_yaml_layer_is_on_the_inner_face` was red on T12 alone (it read the audit's outer-first listing as raw laminate order); now: every skin stack outer to inner is a subsequence of the yaml layer materials
- T25 | S5 | parent + worker | rows restamped (commit below), full suite 11623361 pending | Handoff 4 item 1. Suite 11623302 timed out (4 h, ~65 %, OOM, no summary); re-measured: SLURM 11623360 (`regression --write` for G5/11/20/23/27/32/34/43/58 in the pinned worktree a35bba3, S-4 and V-05 captures), summary in `logs/check/restamp-d5d0fd2/SUMMARY.md`. **Moved anchors (findings):** S-7 torsion ratio (G34) 1.00095 (0.09 %) -> 1.0725 (7.25 %), bound 5 %: **in -> OUT**, consistent with the S-1 GJ move 5.05 -> 10.18 % and the T27 null result; V-02 1F 6.07 -> 7.09 %, 1E 8.50 -> 8.84 % (still OUT). S-4 vs OpenFAST v5 BeamDyn p20, all in: parked 1F +3.42 -> +2.33, 1E +0.20 -> +0.56, 2F +1.80 -> +0.20 %; rated 1F +2.83 -> +1.85, 1E +0.12 -> +0.44, 2F +1.62 -> +0.06 %. V-05 mass -0.47 -> +0.23 %. G23 berna.3 newly measured 4.46 % (in). G58 13 same, G11 same. Rows 5 (mesh convergence) and 20 (twist anchor beam) stay stale: their tests pass but print no `signed error ... bound` residual, so the mapper reads 0 residuals (a test/mapper fix, out of scope). S-0 is out of scope, S-5/S-6 are level evidence only, the CCX parity is T21: none has a store row. `check`: 313 rows, 0 errors, 1 warning (group 27 on no level); `check.sh quick` 75 passed. **Full suite 11623361 (session 5):** 28 failed, 1552 passed, 44 skipped, 18 xfailed, 3 xpassed (4 h 07 min). `scripts/known-red-full.txt` never existed (not on the branch nor origin/main), so the gate counted all 28 as outside. Fast 10 re-run with messages (11623488) and at the merge-base 9f72c6b (11623491): pre-existing at origin/main: `test_blade_mesh` x2 (9357 nodes at base vs 9277 expected; 10637 on the branch) and `test_thin_walled_tube_projection` (torque sign +1.24e5 vs -1.24e5, a sign finding older than the branch); branch-caused: `test_blade_no_webs` x2 (T32, fixed), `test_wall_flow_activation` x2 (T33, fixed), `test_force_projection_load_frame` x3 (T34, open). Slow 18 (blade anchors) re-run with messages: SLURM 11623489, `logs/check/t25/`. Known-red not recorded yet: record only after the 18 are classified. **Slow 18 classified (11623489, a35bba3, 1 h 17 min, 18 failed):** (a) known findings, 9: CCX `tip_axial_tension`/`compression` 7.3 % on the test's tip metric (T21/T31: ill-posed case), S-1 GJ 10.18 %, S-7 7.25 %, V-02 1F 7.09 % / 1E 8.84 %, S-2 `LC1_gravity-utd` 7.60 % (documented red), mesh convergence (mode 6: 0.35 m gap 1.01 % not below the 0.5 m gap 0.89 %; row 5); (b) anchors moved by the branch geometry (findings, not widened), 4: S-1 EI_edge mid-span mean |rel err| 8.01 %, S-1 EA(CLT) 6.68 % (tol 5 %), S-2 official LC1_gravity uy 5.26 % and LC4_distributed uy 5.48 % vs BeamDyn p60 (tol 5 %); (c) disowned UTD input, 2: S-2 `LC2_tip_flap-utd` 17.63 %, `LC4_distributed-utd` 18.33 %; (d) tripwires that fire because the rated-twist findings moved, 4 (T35): `test_rated_twist_under_production_loads` (tip rotation now 1.0438 x Zhou, inside 5 %: promote the row), `test_mesh_free_edges_document_the_open_finding` (no mid-span free edges: tear fixed), `test_section_moment_as_a_couple_sets_the_physical_sign` (nodal-rotation over-report 0.244 no longer > 2 x 0.312), `test_blade_twist_is_load_path_dominated` (4.67 not < 0.1 x 18.26). Every one of the 28 now has a cause. Known-red still not recorded: `--record` re-runs the suite (~4 h), do it after the T32/T33 merge decision so the record matches the tree. **Session 6:** T32/T33 merged (8a81f71); `scripts/check_sbatch.sh full --record` submitted from the pinned worktree `$SCRATCH/wt/full-8a81f71` (shim `.pp`), SLURM 11623853 (bigmem, 8 h), log `$SCRATCH/wt/full-8a81f71/logs/check/full-slurm-11623853.log`; the recorded `scripts/known-red-full.txt` is written in that worktree and has to be copied back and committed. Not merged on purpose: `fix/input-T30` (T24(a) decision open, GJ 100 %), `fix/s1-gj-metric` (null result), `fix/ccx-parity-work` (uncommitted W metric, superseded by T31).
- T26 | S5 | worker | done (finding, no change) | Handoff 4 item 3, DECKY vs T14m (both ad5f3ac, only the yaml differs; pickles read, nothing rerun): LC3 ux +11.27 -> +11.74 % (+0.47 pp, more flexible: the wrong sign for added root TE glass), LC1 -0.30, LC2 +0.94, LC4 -0.01 pp. The deck-era yaml adds TE_reinforcement_SS/PS glass_uni of only 2.98 mm (r/R 0.061) and 8.39 mm (0.082), active only in the width ramp 0.082-0.102. Confounded with the tip chord grid (deck-era 1.16 m vs official 1.708 m at r/R 0.99), the likelier cause of the LC2 move (inferred). Verdict: the root TE reinforcement does not explain the LC3 gap (size <= 0.5 pp vs ~10 pp). Separating root TE and tip chord would need two single-change yamls; not run, the magnitude already answers the question. Note: the S-2 row's +10.26 % is T15m; T14m is +11.27 %. `logs/check/decky/FINDING.md`
- T27 | S5 | worker | done, null result (no change) | Handoff 4 item 4, GJ gap at r/R 0.70 (S-1 FE/BD -19.4 %; checks at 0.45/0.50/0.60/0.80; mesh f354972 = d5d0fd2 blade code, es 0.25 mid). Per-region A66 mesher vs the yaml stack SONATA was fed: <= 0.7 % (residuals are region averaging), so laminate shear contributes +0.2 pt. Thin Bredt -14.8 % (yaml laminate -14.6 %); thick-wall open-section term +2.0 pt (-12.8 %); wall path at the shear centroid -0.2 pt; FE below its own Bredt -4.6 pt here, scattered -7.2..+9.0 pt by station. Unexplained -12.8 pt, growing smoothly with span (-2.9 at 0.45 to -14.1 at 0.80); closing it needs an enclosed area ~8 % larger (wall 17.5 mm outside the mid-surface, past the OML) or wall G x1.17. Lead (inferred): reference-side section geometry, which joins T13 (reference lofting of relative thickness not established). Limit: no per-section SONATA/ANBA output, only BeamDyn K66 and the input yaml; ANBA's 3D lamina law not quantified. `logs/check/gj/FINDING.md`
- T28 | S5 | - | deferred to its own session (would move the geometry under the T25 restamp) | Handoff 4 item 6: second mesh iteration, T16, G6, `element_utils.py:78` divide warning; plus T24(a) TE_FLAT glass extent, measured on S-4/V-02/S-2
- T29 | S5 | worker | row 20 done 88b001e; row 5: cause found, CCX mode 7 non-convergence finding (open) | No test change needed: both tests already print through `assert_residual_below` (`tests/support/assertions.py`) in the declared pattern, and at HEAD `regression --group 20` maps all four comparisons (exit 0, no RED). Row 20 restamped: EA 0.0023 %, GKt 0.0049 % (bound 1 %), FlpStff 1.9239 %, EdgStff 1.8011 % (bound 2 %), all `same`. The T25 `unmapped` for G5/G20 is therefore an artefact of the restamp harness, not the tests; suspect (unverified): the git shim in PATH answers every git command with the hash. G5 rerun from the main repo, its stored high-mode gap 13.98 % already exceeds 5 %. **Correction:** SLURM 11623416 ran `regression --group 5 --write` from the main repo (no git shim) and still returned `unmapped` (0 printed residuals for 2 asserted comparisons, wrote 0, `logs/check/t29/g5.json`), so for G5 the shim hypothesis is not the explanation and the worker's "already prints" was not verified for G5. Next: run the G5 test once with `-s` and read what it prints (it is ~25 min, SLURM) or read the test for the branch where `assert_residual_below` is not reached; row 5 stays stale. **Session 6 (fe47cf7, SLURM 11623984, 22 min):** cause of `unmapped` found by reading the test: the convergence loop asserts before `assert_residual_below`, so a failing loop prints no residual. Also fixed: `_matched_gaps` returned the gaps sorted by size, so the loop compared size ranks, not modes; now indexed by CCX mode. Gaps by CCX mode 1-8 at 0.5 / 0.35 / 0.25 m: 0.04/0.03/0.08, 0.15/0.03/0.05, 0.20/0.06/0.00, 0.43/0.18/0.05, 0.67/0.31/0.14, 1.84/1.34/0.98, **0.89/1.01/1.01**, 13.19/7.97/3.83 %. Prediction partly refuted: **CCX mode 7 does not converge** (gap plateau ~1.0 %, both frequencies fall: AE 4.216/4.199/4.189 Hz; older tree 0.95/0.73/0.52 %, so the branch geometry moved it; the rank bug was not the cause, rank 6 was mode 7 on both meshes). Mode 8 now 3.83 % at 0.25 m, inside 5 % (predicted to fail: refuted). Open decision (user): the code applies NOISE_FLOOR 0.2 % to the coarse gap while its comment calls it the resolvable gap *difference*; mode 7's rise is 0.12 pt, below it; changing that after seeing the data is not done here. Row 5 stays unmapped while the loop fails first.
- T30 | S5 | worker | recorded, NOT merged (fix/input-T30 5edd82c, worktree input-fix-T30) | T24(a) fix: TE band extent [d,e] (`yaml_to_blade.py:449-460`); RED/GREEN on tests/test_blade_te_flat_glass.py, ruff clean; neighbouring tests/test_blade_*.py not completed (login timeout; one F in test_blade_ccx_parity, not compared with d5d0fd2). SLURM 11623402 (S-1/S-2/V-03), 11623403 (S-4/V-02) vs d5d0fd2: as predicted V-03 mass -0.55 % (66,783 kg), EI_edge signed -1.89 -> -5.08 %, S-4 1E +0.56 -> -3.09 % (parked) / +0.44 -> -3.00 % (rated), V-02 1E 8.84 -> 4.90 % (now in band), EI_flap and 1F unchanged, LC3 +10.26 -> +15.67 %; **not predicted:** S-1 GJ mean |rel err| 10.18 -> 100.36 % and S-4 2F +0.20 -> -7.53 % (parked) / +0.06 -> -7.39 % (rated), 2F newly OUT. Parent check (inline, Blade stackdb on both trees): the 104 `*_TE_FLAT` stacks drop from a median summed ply thickness of 30.9 to 1.0 (same units), only Shell_skin, Shell_skin_inner and UV_protection remain: the flat-back face, which closes the TE cell, becomes almost empty, which plausibly explains GJ and the torsion-coupled 2F (inferred, not measured). Open question that decides T24(a): what closes the flat-back in the reference (the yaml arcs, and how SONATA/ANBA built the BeamDyn K66 there); a 1-unit flat-back wall is not physically credible for this blade. `logs/check/t30-FINDING.md` (worker table).
- T31 | S5 | parent | done: tip-load cases ill-posed, replace (not widen); commit below | T21 lead (fibre/material frame) ruled out: the parity test passes span_direction=(0,0,1) explicitly (`tests/test_blade_ccx_parity.py:135,272`), so the T14 default never enters; our fibre (`assembler.rs:257-325`, `geometry.rs:37-118`) vs the CCX *ORIENTATION (`writers.py:1480-1640`) on the whole es 0.25 mesh (35,012 elements): median 0.0004 deg, p95 0.011, max 0.050 (half the 0.1 deg rounding); all plies 0 deg. Sign also excludes it: W_AE < W_CCX, our model is the stiffer one. The 9.6 -> 20.5 % growth comes from another change in the T14m combination tree. Next suspects: elements dropped from the CCX deck (`writers.py:120`), section content per set (thickness, EA), S8R vs MITC4 under the lumped tip load, root BC/load split. Probe caveat: our side is a NumPy port of the Rust formula. `logs/check/t31/FINDING.md` Part 2 (deck read on T7e and T14m, CCX not run): same structure on both sides, 0 dropped elements, 0 missing sections, ply thickness and A11 deck vs ours <= 4.6e-7, sum(A11*area) ratio 1.000000; root BC and tip load equivalent (corner vs corner+midside nodes, 1e6 N). EA per 20 m band moves < 0.6 % T7e -> T14m, while W rises 24 % (ours) and 41 % (CCX): `tip_axial_tension` is an eccentric-load flexure/torsion case (tip response 6.2 m vs ~5 mm pure extension on the prebent, twisted blade), so the growth is in flap/edge/torsion stiffness, matching the S-1 GJ move. Likely (not measured): MITC4 vs S8R on warped, high-aspect TE strips (462 elements aspect > 50 at z > 60 m, warp <= 0.48 edge). Next: CCX `tip_flap`, `tip_edge`, `tip_torsion`, `tip_axial_tension` on both trees, components and W per case; if flap/edge/torsion < 5 % and axial 20 %, replace the case (an eccentricity amplifier), never widen it; es 0.125 as the convergence check. Submitted: SLURM 11623421-24 (T7e flap/edge/torsion/axial, es 0.25), 11623425-28 (T14m same), 11623429 (T14m axial es 0.125, 8 h); prediction in `logs/check/t31/job.sbatch`; torsion is a 5e6 N m force couple about z as tip-node forces (no couple path in the writer); read with `scripts/aeroenv.sh python logs/check/t31/read_percase.py` (untested until outputs exist). **Per-case result (session 5, 11623421-28, `read_percase.py`):** W gap AE/CCX-1, T7e / T14m: flap -2.12 / -2.11 %, edge -0.55 / -0.61 %, torsion -23.01 / -26.29 %, axial -9.56 / -20.53 %. Against the rule fixed before the run: flap and edge in, torsion OUT (predicted 5-12 %), so the pure eccentricity-amplifier reading is refuted; the axial gap tracks the torsional response. Confounder before calling it MITC4 vs S8R: the torsion load is z x r nodal forces on 16 tip nodes with no tip diaphragm (not a wall-tangent shear flow), so W includes loaded-end section distortion; the sign is also odd (ours 23-26 % stiffer than CCX in torsion while the S-1 GJ FE is 19 % below BeamDyn). Probe: `logs/check/t31/station_twist.py` (interior twist rate from the preserved FRDs), SLURM 11623455/56, prediction in `station.sbatch`. 11623429 (T14m axial es 0.125) failed: ccx rc 16 about 1 min into the solve, stdout discarded by percase.py; diagnostic rerun on the kept deck 11623454. **Verdict (station probe 11623455/56/58/59, `logs/check/t31/FINDING.md` Part 3):** the two FE models agree on the blade, GJ twist rate within 0.7 % and axial-case translations within 1 % up to z 82-94 m on both trees; the gap is born in the outer ~10 m and concentrates in the last 1-5 m at the 16 point-loaded tip nodes (no diaphragm; fitted tip rotation 430 -> 2000 deg between 111 and 116 m), so W on the tip nodes measures load introduction, not stiffness. Replace the tip-load parity cases with a well-posed metric (interior station response, or a distributing coupling / diaphragm on both sides); never widen. Corollary (inferred): S-1 GJ -19.4 % vs BeamDyn is not element formulation, supports T27/T13. Open: z 82-111 m torsion rate -4.5 %, edge translation up to +6.9 % (T28). es 0.125: CCX out of memory (`u_realloc`, 21.7 GB at --ntasks=8), not resubmitted.
- T32 | S5 | parent | done, merged into l3/twist-p1 73d5b3f (session 6, user L23) | T25 triage: `test_blade_no_webs` (2) fails on the branch, passes at the merge-base 9f72c6b. Cause: c157c62 splits the spar caps at the web attachment arcs only `if include_webs` (`mesh_gen.py:344`), so the no-webs surface mesh (the CFD/coupling surface, `generators.py:1054`) has a different skin discretisation (2810 vs 3516 outer-shell elements at es 1.0) and no longer matches the structural skin. Fix: plan the cap split whatever `include_webs`; `webHP/webLP` are only read in the web loop, which `include_webs=False` skips. RED 2 failed (2810 == 3516, 4065-2810 == 549); GREEN `test_blade_no_webs.py` + `test_blade_web_position.py` 11 passed, ruff clean. Not re-measured: the with-webs mesh is unchanged by construction (only the no-webs branch of the condition moves).
- T33 | S5 | worker | wall-flow done, merged into l3/twist-p1 8a81f71 (session 6, user L23); load-frame split to T34 | T25 triage: `test_force_projection_load_frame` (3) and `test_wall_flow_activation` (2) fail on the branch, pass at 9f72c6b. Probe `logs/check/t25/ring_probe.py` (SLURM 11623493): every one of the 50 None rings (42 no closed cell, 7 unusable area, 1 group<3) sits exactly on a BEM strip edge (z 1.194, 3.582, 5.969, ... step 2.388 m) and is split between the two neighbour strips (e.g. 35 + 31 nodes at z 1.194); the 19.87 deg chord angle is the split ring at z 3.582. Cause: `ForceProjector.__init__` assigns each node by the half-open `[r_low, r_high)` of its own span coordinate (`force_projection.py:375-378`), so a mesh row on a strip edge with round-off straddles it; the branch's mesher now puts rows on those edges. Fix direction: assign whole physical rings (or snap span coordinates within `_ring_gap_tolerance` of an edge) so a ring never straddles a strip. Note (unverified): the probe also printed `props for 11758 ids, coupling elements without props 28483`. **Result (worker):** assignment-only snap of span coordinates within `_ring_gap_tolerance` (1e-4 x span) of a strip edge, interior edge to the higher band; new `tests/validation/bem/test_force_projection_strip_edge_ring.py` RED 2 failed (`strip-edge ring split across strips: [4, 4]`) -> GREEN 4 passed; `tests/validation/bem/` 3 failed, 133 passed, 12 skipped, 1 xfailed; `test_wall_flow_activation.py` 7/7 passed; the 3 load-frame failures keep `z=3.582 m: c_hat is 19.874 deg` (T34); thin-walled tube message unchanged (pre-existing); ruff clean.
- T34 | S5 | parent | done: guard ill-posed on the near-round root ring, made well-posed (commit below) | Probe `logs/check/t34/probe.py`/`probe2.py` (SLURM 11623859/60, 11623861/62): base 9f72c6b and 8a81f71 pick the **same** frac 0.03 station z 3.5816 (68 nodes) with the **same** outline, min/max extent 0.9672 / 0.9678 (prediction that the base picks another station: refuted; root outline unchanged, so not T24). angle(c_hat, max-extent) z 3.10/3.34/3.58/3.82/4.06: base 9.3/8.1/6.9/5.9/5.0, head 9.9/21.4/19.9/7.6/5.7 deg. Arc-length-weighted PCA vs max-extent: base 0.2-1.2 deg, head 12.4/12.8 deg at 3.34/3.58: the max-extent argmax jumps on a flat maximum, and the node-weighted c_hat is itself 5-8 deg off the arc-length axis on both trees. So at aspect ~0.97 both estimators follow node placement. Fix (test only, production projector untouched): the angle check is asserted only where the outline resolves its major axis at the guard's own 10 deg bound (extent drop at +-10 deg > circular n-gon error 1-cos(pi/n)), the model applied only to near-round rings (`NEAR_ROUND_ASPECT` 0.90, a domain, not a tolerance: any value in (0.51, 0.96) gives the same verdicts; first version without it flagged the 16-node tip ring, caught by the new confinement test); the chord-extent check stays on every ring; unresolved rings leave the twist witnesses; new `test_unresolved_major_axis_only_at_the_round_root` (frac >= 0.15 must resolve). SLURM 11623871: 6 passed at 8a81f71 and 6 passed with the same file at 9f72c6b. Marginal by construction: z 3.58 drop 1.28e-3 (base, resolved, 6.9 deg asserted) vs 1.02e-3 (head, skipped) against 1.07e-3. Derived after seeing the data (stated). Finding kept: the node-weighted c_hat is not a discretisation-invariant chord estimator near the root.
- T35 | S5 | - | open | Four rated-twist tripwires fired on the branch (T25 (d)): each failure message says what to do (promote the rated-twist row to an asserted magnitude vs Zhou, update the free-edge test and the WU7 section of `odd/tasks/rated-twist-sign-convention.md`, re-derive the couple-vs-nodal-rotation and load-path-dominance claims). Measure first which branch commit moved them (mesh tear fixed: likely the T8-T15 mesher fixes); then rewrite the guards with the new numbers and references, in their own session.
- T36 | S5 | worker + parent | done, merged (fix/ccx-parity-station dcaa07f) | Handoff 6 item 3, T31 verdict: replace the tip-node metric of `tests/test_blade_ccx_parity.py::test_blade_parity_ccx` (mean |u| over the 16 point-loaded tip nodes = load introduction, not stiffness) with an interior-station metric, on branch `fix/ccx-parity-station` (worktree `$SCRATCH/wt/ccx-station`, from a9ca498; supersedes `fix/ccx-parity-work`'s W metric). **Design fixed before any run (session 6):** per case, the mean translation vector over a node slab |z - zc| < 0.6 x ELEMENT_SIZE (CCX corner nodes matched by coordinates, cKDTree, 1e-4 m) at zc = z0 + f (z1 - z0); metric |u_AE - u_CCX| / |u_CCX|; **asserted** at f = 0.5 and 0.7 with the unchanged 5 % bound; f = 0.8, 0.9 and the old tip metric **printed, not asserted** (0.8-0.9 R is the open T28 band: edge translation up to +6.9 %); all four cases (flap, edge, axial +/-). **Prediction:** axial +/- <= ~1 % at 0.5/0.7 (T31 Part 3: <= 0.85 % up to z 82 m); flap and edge < 3 %. If a case is > 5 % at 0.5 or 0.7 R it is a finding (report, never widen, never move the stations). **Result (session 6, SLURM 11623955-58, es 0.25, 7.7 min each):** all 4 pass. Slab-mean vector gap AE/CCX at f 0.5 / 0.7 (asserted) / 0.8 / 0.9 (printed): flap 0.16 / 0.11 / 0.08 / 0.04 %, edge 0.17 / 0.08 / 0.11 / 0.16 %, axial +/- 0.15 / 0.09 / 0.06 / 0.09 %; every slab fully matched (104/78/72/60 nodes); old tip metric (diagnostic) flap 2.14, edge 0.75, axial 7.25 %. Prediction (axial <= ~1 %, flap/edge < 3 %) held. f 0.9 shows no excess in the vector metric; T31's +6.9 % at 82-111 m was the minor x component of the axial case (here 43 vs 828 mm at 0.7 R), not contradicted. **Matching, four iterations (findings):** (1) nearest node within a fixed 1e-4 m misses everything beyond z 100 m (FRD E12.5 rounding 5e-4 m there); (2) the FRD carries 1.83 M nodes (CCX expands composite shells through the thickness) vs 102,667 in the deck; (3) even on deck corners a coordinate match is ambiguous for 579 nodes: our mesh has distinct nodes 0.25 mm apart; (4) fix: the quadratic writer labels our node i as i + 1 (`_build_quadratic_mesh_data`), so pairing is by identity and coordinates only verify it within the rounding (corner set must be exactly 1..N). The run-2 nearest-node numbers coincided with these at 0.5/0.7, but only the label-paired run is evidence. Commits 5282da8, 125c604, 12ec1ef, 4a69c4c, 57a7ff9, dcaa07f; logs `$SCRATCH/wt/ccx-station/logs/check/t36/`.
- T37 | S5 | parent + worker | in progress (session 7): reference confirmed, fill is glass not adhesive, see section 'T37 flat-back closure' | Flat-back (TE) closure to match the reference. Findings (session 6, two read-only explores): (1) the yaml places **no layer on the flat-back face** outboard of r/R ~0.245: the FFA-W3-360/330blend/301/270blend polylines are open at the TE, so the face has zero nd-arc; Shell_skin starts at nd 0.01 / Shell_skin_inner at 0.08 outboard of 0.3/0.25, TE_reinforcement_SS/PS (glass_uni 30 mm, 0.75 m) run from each corner forward; the face laminate is a tool assumption. (2) Ours: `*_TE_FLAT` key segment has zero arc length (`keypoints.py:237-239,285-287`), `stackdb._map_continuous_layers` `if b <= a:` gives fraction 0/1 by the yaml start/end; HEAD (`yaml_to_blade.py:449-460` extent [d,te]) puts glass_uni 30 + gelcoat 1 = 31 mm on the face, T30 ([d,e]) leaves gelcoat 1 mm only (G t ~1.3 MN/m: cell effectively open, hence GJ 100 % and 2F -7.5 %). Units: mm. (3) Reference = SONATA/ANBA4 (IEA-15-240-RWT ReleaseNotes v1.1.4: 'Blade elastic properties consistently originate from SONATA/ANBA4'). SONATA master: closes the profile with a straight line between the corners (`classAirfoil.gen_OCCtopo` :125, `BSplineLst_utils.py:481-484`), s=0 at the face midpoint (`classAirfoil.py:86-92`, `classBlade.py:402-411`); [0,1] layers wrap across the face (`utl/converter_WT.py:540-557` snaps near-0/1); **outboard of 20 % span** (`span_adhesive = 0.2`, :131) TE reinforcement is trimmed back to `adhesive_extent = min(0.04, 0.04/chord)` (:323-331), an **adhesive web** is inserted at [adhesive_extent, 1-adhesive_extent] (:521-524) and the TE wedge behind it is **solid adhesive** (:422-432); inboard of 20 %: TE glass_uni + skins + gelcoat on the face. Confidence medium on the adhesive details (master, not the v1.1.4 commit). **So neither HEAD (31 mm glass) nor T30 (1 mm) is the reference:** the reference torsion cell is closed by an adhesive bond line ~4 % of arc forward of the TE plus a solid adhesive wedge. Next: confirm on `IEA-15-240-RWT/CAD/SONATA/blade_section_0.3.pdf` / `_0.5.pdf`; then model it in the shell mesh (a TE shear web of adhesive at adhesive_extent, or an equivalent face laminate with the wedge's G t and E A), keep T30's [d,e] for glass_uni outboard of 0.2, measure S-1 GJ/EI_edge/EA, S-7, V-02 1E/1F, S-4 1E/2F, V-03 mass against the current numbers; prediction to write before the run. The yaml's `adhesive` material exists unused (explore 1).
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

## T37 flat-back closure (session 7, 2026-10-25; user L25)

**Reference geometry confirmed on the SONATA PDFs** (`IEA-15-240-RWT/CAD/SONATA/blade_section_<r>.pdf`, vector triangles
read with PyMuPDF, scale from the axis ticks; scripts `logs/check/t37/wedge.py`, `mass.py`, `stiff.py`):

| r/R | 0.15 | 0.2 | 0.25 | 0.3 | 0.4 | 0.5 | 0.65 | 0.75 | 0.9 |
|---|---|---|---|---|---|---|---|---|---|
| solid TE region behind a straight bond line [mm chordwise] | 42 | 50 | 166 | 150 | 162 | 184 | 198 | 208 | 206 |
| its area [cm2] (glass_uni-coloured / total) | 19/45 | 18/30 | 166/181 | 170/177 | 147/155 | 112/120 | 52/61 | 22/31 | 18/25 |

Inboard (0.15, 0.2) the face is closed by the laminate only (triax + glass_uni + gelcoat, 42-50 mm = the corner);
from 0.25 outward a straight bond line closes the cell 150-210 mm ahead of the TE and the region behind it is solid,
wrapped by gelcoat (4 cm2), the skins trimmed at the bond line (triax 2-11 cm2).
**Finding 1: the fill is glass, not adhesive.** Its triangles carry exactly the glass_uni colour; the colour bar lists 9
materials (Gelcoat ... medium_density_foam) and no adhesive. SONATA master (`converter_WT.py`) fills with `Adhesive`, and
the yaml `adhesive` is E 4.56e6 Pa, G 1.52e6 Pa (1000x below an epoxy, structurally void; likely a unit slip).
**Finding 2: the PDF sections are the 6x6.** Integrating the PDF triangles (E1 x area per material, fill = glass_uni):

| r/R | 0.25 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 |
|---|---|---|---|---|---|---|
| mass vs BeamDyn | +0.00 % | +0.00 | +0.00 | +0.00 | +0.00 | +0.00 |
| EA | -0.1 % | -0.1 | -0.1 | -0.1 | -0.1 | -0.1 |
| EI_edge (principal) | -0.2 % | -0.1 | -0.1 | -0.1 | -0.2 | -0.3 |
| EI_flap (principal) | -0.1 % | -0.1 | -0.0 | -0.1 | -0.1 | -0.1 |
| solid TE region share of EA / of chordwise EI | 3.7 / 26.8 % | 3.6 / 28.5 | 3.0 / 27.5 | 2.5 / 24.0 | 1.8 / 22.4 | 1.2 / 26.5 |

With an adhesive fill (rho 1100) the mass would be ~3 % low at 0.3: the reference fill has glass density and stiffness.
So the 6x6 is reproduced to <= 0.3 % by these sections, which makes them a per-region arbiter (also for T27/T13), and
the TE solid region carries a quarter of EI_edge: the TE construction is a first-order EI_edge input, not a GJ one.
**Correction to Handoff 7:** by Bredt the closure is ~1 % of the cell's integral of ds/(G t) at HEAD (31 mm glass_uni
face, G t ~ 1e8 N/m against ~1e8 N/m of skin over ~11 m); a rigid closure can raise GJ by at most ~1 %. T37 cannot
close the -10 % GJ gap or S-7.

**Modelling decision (no web).** The solid region is 150-210 mm long, shorter than one 0.25 m element, and it is glass:
a bond-line web would create a sub-element cell and a third web in the mesher topology to carry a material that acts
as a lumped axial member. Represent it as an equivalent glass_uni laminate on the `*_TE_FLAT` elements outboard of
r/R 0.2 (SONATA `span_adhesive`), with area = the cavity between the inner skin surfaces over the last L_b = 0.19 m
(median of the PDF row above, a construction input, not a fit to the 6x6), so skins + fill = the solid OML region;
the TE_reinforcement glass_uni comes off the face (T30's [d,e]); inboard of 0.2 the HEAD face stack (skins + gelcoat
+ glass_uni) is already the reference.

**Section-integral probe (SLURM 11624106, worktree `wt/t37` a2dc3b0, es 0.25 mid mesh, tools `tools/sonata_pdf_section.py`,
`tools/t37_section_compare.py`, log `wt/t37/logs/check/t37/`).** E1 x area integrals vs BeamDyn, signed, current -> T37
(face-equivalent rule above, L_b 0.19 m):

| r/R | EA | EI_edge | EI_flap | mass | TE solid region area, ours / PDF [cm2] | face glass now -> T37 [mm] |
|---|---|---|---|---|---|---|
| 0.25 | +5.8 -> +7.1 | -0.8 -> +9.0 | +4.4 -> +4.8 | -1.5 -> +0.7 | 142 / 181 | 30 -> 105 |
| 0.3 | +6.3 -> +7.9 | -1.6 -> +11.8 | +4.7 -> +5.1 | -0.2 -> +2.9 | 138 / 176 | 30 -> 95 |
| 0.4 | +7.9 -> +8.7 | +3.1 -> +10.6 | +6.5 -> +6.6 | +3.1 -> +4.9 | 139 / 155 | 30 -> 69 |
| 0.5 | +8.6 -> +8.3 | +11.0 -> +7.9 | +7.0 -> +7.0 | +5.8 -> +5.1 | 143 / 120 | 30 -> 8 |
| 0.6 | +8.0 -> +7.9 | +11.2 -> +9.1 | +6.6 -> +6.6 | +5.7 -> +5.3 | 97 / 79 | 20 -> 2 |
| 0.7 | +5.7 -> +6.0 | -5.6 -> +1.9 | +4.9 -> +4.9 | +2.4 -> +3.3 | - | 4 -> 82 |

Findings: (a) inboard (0.25-0.4) our TE region lacks 10-22 % of the reference solid; outboard (0.5-0.6) it holds
19-24 % MORE than the solid, because the two mid-surface skins (~34 mm each) overlap in a 35-61 mm thick TE: the face
rule cannot remove that excess and drops the face glass to 2-8 mm there. (b) The whole section is too stiff/heavy
beyond the TE: EA +6-9 %, EI_flap +4-7 %, mass up to +6 % at 0.5-0.6 (T27/T13 territory). (c) The section integral
does not reproduce the FE S-1 EI_edge sign (integral -1.6..+11 %, FE mean signed about -7.5 %); S-4 1E is +0.56 % vs
BeamDyn, consistent with the integral, so the S-1 FE EI_edge extraction looks biased low by ~8 pt (open).

**Prediction (written before any anchor run), T37 face-equivalent vs the T25 numbers:**
- S-1 EI_edge: FE rises ~+5 % mid-span (integral deltas +10/+14/+7/-3/-2/+8 %): mean |rel| 8.0 -> 5-6 %.
- S-1 EA(CLT): 6.7 -> ~7.0 % (+0.3-0.5 pt, worse). Mass per length +0.5..+3 % inboard; V-03/V-05 blade mass
  +0.23 -> +0.7..+1.5 %.
- S-1 GJ: closure stiffer at 0.25-0.4 (<= +1 %), weaker at 0.5-0.6 (face glass 2-8 mm: -1..-4 %): mean |rel| 10.2 ->
  10-12 %; S-7 ratio 1.0725 -> 1.07-1.09 (no improvement).
- 1E: inboard EI_edge up with TE mass up: S-4 1E +0.56 -> +2..+3.5 % (worse vs BeamDyn), V-02 1E +8.8 -> +10..+12 %.
- Flap: EI_flap <= +0.4 %, inboard mass up: V-02 1F -7.1 -> about -7.5 %, S-4 2F +0.20 -> -0.3..-1 %.
So the face-equivalent T37 is predicted NOT to improve the anchors: the reference wedge is real, but our section
already compensates it with excess elsewhere. Run it anyway to test the prediction and the integral -> FE translation.

## Handoff 7 (2026-10-25, session 6; user L23-L24)

**Resume.** Read this section and the Tasks lines T25, T29, T36, T37. Engram still down; this file is the record.

**Done.** T32+T33 merged (73d5b3f, 8a81f71). T34 done (a9ca498, test only: the root-ring guard measured node placement). T36 done and merged (7eaba0b): CCX parity on interior slabs 0.04-0.17 %, so AE and CCX agree as solvers and the remaining anchor gaps vs BeamDyn are input/section physics. T29 row 5 (fe47cf7, d9309d5): gaps now per CCX mode; **CCX mode 7 does not converge** (0.89/1.01/1.01 %); mode 8 3.83 % at 0.25 m. User L24: the priority is to **correct the physics**, not test plumbing.

**Running.** SLURM 11623882: `check.sh full --record` on the pinned worktree `$SCRATCH/wt/full-a9ca498` (a9ca498, before T36/T29). When done: copy `$SCRATCH/wt/full-a9ca498/scripts/known-red-full.txt` to the main repo, commit it (stamp a9ca498); compare its failures with the 28 classified in T25 (expect 26: T32/T33 green; T34 green); T36 will show `now-passes` for the CCX `tip_axial_*` on the next run. Log `$SCRATCH/wt/full-a9ca498/logs/check/full-slurm-11623882.log`.

**Next, in order (physics first).** 1. T37 flat-back closure as SONATA builds it (adhesive bond line + wedge outboard of 20 % span): expected to move S-1 GJ (-19 %), S-7 (7.25 %), V-02 1E (8.8 %) and S-4 2F together. 2. T27/T13 the rest of the GJ gap (reference section geometry). 3. T35 rated-twist tripwires (own session). 4. T28 mesh iteration (mode 7 non-convergence belongs here).

**Cleanup.** Worktrees `full-8a81f71`, `ccx-station`, `base-9f72c6b` (+ `.pp`) and those of Handoffs 5-6; `$SCRATCH/tmp/t34-*`, `t36-*`.

## Handoff 6 (2026-10-24, session 5; user L22)

**Resume.** Read this section and the Tasks lines T25, T31-T34. Engram not used this session; this file is the record.

**Done this session.**
- T31 closed (b85665c): the CCX parity gap is tip load introduction (16 point-loaded tip nodes, no diaphragm); the two FE models
  agree on interior GJ within 0.7 % and on axial-case translations within 1 % to z 82-94 m. Replace the tip-load cases, never widen.
  Corollary (inferred): S-1 GJ -19.4 % vs BeamDyn is not element formulation (supports T27/T13).
- T25 suite 11623361: 28 failed / 1552 passed. Fast 10 classified (pre-existing 3, branch-caused 7). `known-red-full.txt` does not exist.
- T32 (fix/input-T32 e6327ec): the no-webs surface mesh now keeps the cap split, so its skin matches the structural one. Not merged.
- T33 (fix/bem-T33 a7a3e59): BEM strip assignment no longer splits a ring lying on a strip edge; wall-flow guard green. Not merged.

**Running.** SLURM 11623489: the 18 slow T25 failures (blade anchors) with messages, `logs/check/t25/rerun-t25-slow-11623489.log`
(`grep -aE "^(_{3,}|E  )"`). Expected from earlier findings: S-7, V-02 1F/1E, S-1 GJ, CCX `tip_axial_*`, S-2 `LC1_gravity-utd`, mesh convergence
(row 5). To classify: S-1 edgewise EI and EA CLT, S-2 official LC1/LC4 and UTD LC2/LC4, rated twist x2, twist mechanism x2; run any that look
branch-caused at the merge-base worktree `$SCRATCH/wt/base-9f72c6b` with `logs/check/t25/rerun-any.sbatch` (`--export=ALL,W=...`).

**Next, in order.** 1. The 18 are classified (T25 line): decide the T32/T33 merges, then `scripts/check.sh full --record` (~4 h, SLURM) on that tree. T35 (rated-twist tripwires) in its own session.
2. T34 (load-frame station at z 3.582). 3. Replace the CCX tip-load parity cases (T31 verdict) on `fix/ccx-parity-work`.
4. T30/T24(a) flat-back decision, T29 row 5, T27, T28 as in Handoff 5. 5. Merges of fix/input-T32, fix/bem-T33 into l3/twist-p1 are the
user's call (T32 moves the no-webs mesh only; T33 moves BEM strip ownership of edge rings only).

**Cleanup.** Worktrees `base-9f72c6b`, `input-fix-T32`, `bem-fix-T33` (+ `.pp`) and those listed in Handoff 5; `$SCRATCH/tmp/t31-*`, `t25-*`, `t33`.

## Handoff 5 (2026-10-24, session 4 close at ~170 k context; user L20/L21)

**Resume.** Read this section and the Tasks lines T21, T24-T31. Engram was unavailable all session (server ownership
mismatch at 127.0.0.1:7437), so this file is the only record; fix Engram or keep resuming from here.

**State.** `l3/twist-p1`, local only, not pushed, about 8.2 k authored lines in 90 commits since origin/main; every
commit of this session is docs/validation or docs/odd. No blade code changed since d5d0fd2. Candidates not merged:
`fix/ccx-parity-work` (T21, uncommitted metric diff in worktree input-fix-CCXW) and `fix/input-T30` (5edd82c).

**Done this session.**
- T25: blade rows restamped at d5d0fd2 (c4b548d). S-7 torsion ratio 0.09 -> 7.25 %, crosses 5 % (finding). V-02 1F 7.09, 1E 8.84 % (still out).
  S-4 all six in band vs BeamDyn p20. Row 20 restamped (88b001e). Row 5 still stale (T29).
- T21: load-work CCX parity fails 5 % on every tree (9.6 / 10.6 / 20.5 %); recorded, not merged.
- T26 DECKY: root TE reinforcement does not move LC3 (+0.47 pp, wrong sign).
- T27 GJ: null for laminate shear (+0.2 pt) and thick wall (+2 pt); -12.8 pt unexplained, reference geometry suspected (joins T13).
- T24: (a) TE_FLAT glass is a translation defect; (b), (c) are the mid-surface inset (yaml-faithful); OML root surplus open.
- T31: no fibre-frame or deck-content mismatch in the CCX parity; `tip_axial_tension` is eccentric flexure/torsion.
- T30: the T24(a) fix empties the flat-back wall (median 30.9 -> 1.0); GJ 100 %, 2F -7.5 %; recorded, not merged.

**Running on SLURM; read the results, do not rerun.**

| job | what | read with |
|---|---|---|
| 11623361 | full suite at a35bba3 (pinned worktree `$SCRATCH/wt/full-a35bba3`), 8 h bigmem; was at ~66 % after 2 h 47 min | `logs/check/full-a35bba3-11623361.log`; check log copied to `logs/check/full-*.log` |
| 11623421-28 | T31 CCX per case (flap, edge, torsion, axial) on T7e and T14m, es 0.25 | `scripts/aeroenv.sh python logs/check/t31/read_percase.py` (untested) |
| 11623429 | T31 T14m axial at es 0.125 (convergence) | same |

**Next session, in order.**
1. **Close T25** from 11623361: compare the failures with `scripts/known-red-full.txt`; re-derive known-red only with
   `scripts/check.sh full --record` on that run; any new failure outside the blade anchors is a finding.
2. **T31 verdict** from the per-case table, with the rule fixed before the run: flap/edge/torsion < 5 % with axial ~20 %
   means the case is an eccentricity amplifier, to be replaced by a well-posed one (never widened); otherwise a real
   MITC4 vs S8R stiffness gap on the warped TE strips (462 elements aspect > 50 at z > 60 m).
3. **T30/T24(a) decision:** what closes the flat-back in the reference? Read the yaml arcs at the flat-back and how
   SONATA/ANBA meshed it for the BeamDyn K66. Then: keep the glass (if the reference closes the cell with it), or move
   the closing laminate (if it does not), and re-measure GJ and 2F. Do not merge 5edd82c as it stands.
4. **T29 row 5:** run the G5 test once with `-s` (~25 min, SLURM) and find why `assert_residual_below` prints nothing
   for the mapper.
5. **T27 lead:** the reference section geometry (T13, relative-thickness lofting), the -12.8 pt GJ gap growing with span.
6. **T28** (own session): second mesh iteration, T16, G6, `element_utils.py:78`.
7. Before any merge to main: the full suite and maintainer review; slicing the ~8.2 k-line branch is the maintainer's call.

**Cleanup when done.** Worktrees: `restamp-a35bba3` (+ `.pp`, `$SCRATCH/gitshim-d5d0fd2`), `input-fix-T30` (+ `.pp`),
`full-a35bba3` after 11623361; decks under `$SCRATCH/tmp/t31-*`.

## Handoff 4 (2026-10-24, session 3 close; user: "haz el resellado de L3 y luego cierra la sesion")

**Resume.** `mem_search "blade-input-audit"`, read this section, then "GJ gap and reference provenance" and the
Tasks list (T21, T24). Use a new, bounded session; this one carried the whole T9 integration.

**State.** `l3/twist-p1` at the Handoff 4 commit, local only, not pushed. The blade candidates of the first mesh
iteration are merged:
- T7d, T7e, T10, T8, T8w, T11, T12;
- G10 and G10b (linear span interpolation of the layers);
- T14 and T15;
- T22.

The final-tree measurement is d5d0fd2, `capture-T15m.pkl`, SLURM 11623286.

**Restamped (a35bba3).** The ten S-1/S-2/V-03 store rows are at d5d0fd2. `contract.yaml` is not touched (audit
rule). Every other store row whose test builds the IEA-15 blade is stale: S-0, S-3, S-4, S-5/S-6, S-7, V-01/V-02/V-05,
the CCX parity, and anything else using `BladeMesh`/`Blade` on the official yaml.

**Running on SLURM; read the results, do not rerun.**

| job | what | output |
|---|---|---|
| 11623302 | full suite at a35bba3 (pinned worktree `$SCRATCH/wt/full-a35bba3` + shim; `logs/check/full-a35bba3.sbatch`) | `logs/check/full-a35bba3-11623302.log`, check log copied to `logs/check/full-*.log` |
| 11623281-83 | T21 work metric `tip_axial_tension` on T7e, T8w, T14m | `$SCRATCH/wt/input-fix-CCXW/logs/check/t21/work-w-<tree>-tip_axial_tension.log` |
| 11623284 | T21 full CCX parity, moved to bigmem | `$SCRATCH/wt/input-fix-CCXW/logs/check/t21/ccxw-full.log` |
| DECKY | ad5f3ac with the deck-era yaml | `logs/check/input-fix/capture-DECKY.pkl`; compare with `capture-T14m.pkl` |

The T21 diff (`tests/test_blade_ccx_parity.py`) is uncommitted in `$SCRATCH/wt/input-fix-CCXW` (branch
`fix/ccx-parity-work`). Its unit test showed RED (2 failed), then GREEN (2 passed); ruff is clean.

**Next session, in order.**
1. **Restamp every blade-dependent store group from 11623302**, the L39 procedure:
   - compare failures with `scripts/known-red-full.txt`;
   - run `python tools/validation_matrix.py regression --group <N> --scope <test files> --write` per group;
   - record moved anchors as findings;
   - re-derive the known-red list only with `scripts/check.sh full --record` on that run.
2. **T21.**
   - Read the work-metric logs.
   - If the combination passes 5 %: commit in CCXW, merge, and restamp the CCX parity rows.
   - If not: record it as a finding and change nothing.
3. **DECKY.** Does the root TE reinforcement of the deck-era yaml move S-2 (LC3) at all?
4. **GJ gap** (open; web placement ruled out by the SONATA PDFs). Next cheap step:
   - compare one section's laminate shear (A66 per region) against SONATA/ANBA, or run a thick-wall/warping
     Bredt variant;
   - the root reference interpolation (log-interp `gj_bd` at 0.19-0.25) is a separate metric check;
   - the FE-side wobble at 0.817, where the TE reinforcement, the fillers and the LE reinforcement end together.
5. **T24 mass findings.**
   - TE_FLAT carries 367 kg of glass_uni on the flat-back face (the band carried onto it?).
   - Skin triax -215 kg and carbon -214 kg against the yaml.
   - Webs -4.8 %.
6. **Second mesh iteration** (user L8): T16 (multi-blade fibre direction), G6 (prebend tilt), the
   `element_utils.py:78` divide warning, and the T20 pchip option only on demand.
7. Before any merge to main: the full suite (item 1) and maintainer review. The branch is about 5.9 k authored
   lines in 74 commits since origin/main; each candidate is its own commit, so the slices exist.

**Cleanup when done.** `git worktree list | grep -E 'input-fix|full-a35bba3|DECKY'`, then `git worktree remove` for
the merged or finished ones. `$SCRATCH/tmp/{input-fix-*,t9-*,g10*,t15*,gj-*,sonata-*,deck-yaml,ccx*}`.

## Handoff 3 (2026-10-24, session 3; written while jobs run)

**Resume.** `mem_search "blade-input-audit"`, then read this section, then "T9 integration", "GJ gap and reference
provenance" and the Tasks list.

**Branch.** `l3/twist-p1` at bd18fae, local only. Merged this session:
- the T9 candidates (T7d through T12);
- G10 (linear layer thickness) and G10b (linear nd_arc bounds);
- T14 (span_direction default) and T15 (G8, G9, ValueErrors, d from corner);
- T22 (ply-order test fix).

Not merged:
- `fix/s1-gj-metric` (a null result);
- `fix/ccx-parity-work` (T21, worker in flight).

**Results that land on disk; read them, do not rerun.**

| what | output | compare against |
|---|---|---|
| G10b | `logs/check/input-fix/capture-G10b.pkl`, `pytest-G10b.log` | `capture-G10.pkl` |
| T15m (d5d0fd2: G10b + T15) | `capture-T15m.pkl` | `capture-T14m.pkl` (ad5f3ac) |
| DECKY (ad5f3ac with the deck-era yaml `deef925e`) | `capture-DECKY.pkl` | `capture-T14m.pkl`; separates the root TE reinforcement and the tip chord grid |
| T21 work metric | branch `fix/ccx-parity-work`, worktree `$SCRATCH/wt/input-fix-CCXW`, `logs/check/t21/` | |
| GJ dip 0.66-0.74 | `logs/check/gj-dip/REPORT.md` | |

Run each comparison with `python logs/check/input-fix/compare.py <a.pkl> <b.pkl>`.

**Next, in order:**
1. Read the five results.
2. T21: if the work metric passes on the combination, merge `fix/ccx-parity-work` and run the blade suite with CCX
   on SLURM (`logs/check/input-fix/t9-blade.sbatch` pattern).
3. Act on the GJ dip report.
4. T9.7: restamp the L3 rows on the final tree, against 9f72c6b and T9. Register every moved anchor as a finding.
5. Full suite on SLURM before any merge to main.

**Rules (user L18).** Do any unblocked work toward a physically realistic model without asking. No band, anchor or
contract.yaml change. One candidate per commit, RED first, measured on SLURM.

## T9 integration (session 3, 2026-10-23, user L14)

**Merges into `l3/twist-p1`.** ceca579 T7e (contains T7d), b7dafdd T10, 93df15e T8, 4c9372d T8w, ea28216 T11,
61cf5b0 T12.
- **Conflict.** One textual conflict, T7e x T8w in `mesh_gen.py`. The chordwise counts are keyed
  `("shell", j, k)`, with `k` the cap-piece letter, taken before the zero-width filter.
- **Finding: sliver cap pieces in the native T8w mesh.** T8w decides piece presence per bay and tapers the width
  linearly, so the native mesh (`laminate_per_row=False`) carries a 0-5 mm wide cap piece in 14 rows. Native minimum
  node gap 0.3-4.4 mm against 75-77 mm decoupled.
  - Root bay where the yaml webs coincide: z 7.40-9.31.
  - Tip: z 114.7-115.1.
  - The decoupled default uses the exact web arcs at every inserted station and drops them: 28 nodes fewer, same
    elements (34,026 vs 34,054 / 35,211). Logs: `logs/check/input-fix/t9-nodes.log`, `t9-sliver.log`.
  - `test_blade_span_decoupling` now states this contract (elements and row z equal; every node difference is a
    native sliver row). Fragile point: its decoupled minimum-gap floor of 1 mm against 1.2 mm measured at the tip.
- **Checks.** `check.sh quick` OK (75 passed, 1 xfailed).
- **Finding: the CalculiX `tip_axial` parity breaks in the combination.** `test_blade_ccx_parity[tip_axial_tension]`
  and `[tip_axial_compression]` pass on T7e alone and on T8w alone (SLURM 11623249: 2 passed each,
  `logs/check/input-fix/t9-blade-T7e.log`, `t9-blade-T8w.log`). They failed on the T7e x T8w merge (login run before
  T11/T12). Full-combination result: `t9-blade-T9.log`.
- **Measured, combination at 61cf5b0** (SLURM 11623247, `submit_wait.sh T9`; `compare.py
  ../restamp-9f72c6b/capture-s1s2v3.pkl capture-T6 capture-T7e capture-T8w capture-T9`). Signed %:

  | metric | 9f72c6b | T6 | T7e | T8w | T9 |
  |---|---|---|---|---|---|
  | S-1 EA | +7.57 | +8.03 | +8.03 | +8.11 | +6.25 |
  | S-1 EI_edge | -7.51 | -3.77 | -4.18 | -2.96 | -3.71 |
  | S-1 EI_flap | +3.79 | +3.77 | +4.61 | +4.88 | +3.89 |
  | S-1 GJ | -2.76 | -1.37 | -1.26 | -11.04 | -10.36 |
  | S-1 mass | +0.73 | +1.74 | +1.58 | +1.65 | +1.65 |
  | S-2 LC1 gravity | -7.44 | -6.44 | -6.56 | -6.02 | -5.59 |
  | S-2 LC2 tip flap | -3.47 | -3.49 | -3.42 | -2.31 | -2.43 |
  | S-2 LC3 tip edge | +15.57 | +12.00 | +12.11 | +12.13 | +11.31 |
  | S-2 LC4 distributed | -6.27 | -6.28 | -6.33 | -5.81 | -5.47 |
  | V-03 mass [kg] | 66,682 | 66,903 | 66,908 | 66,921 | 66,941 |

  - Failing anchor tests: 5, against 6 on T7e. `test_flapwise_ei_within_band` now passes. Still failing: S-1 EA,
    EI_edge and GJ bands; S-2 LC1 and LC4 official.
  - **Moved rows, findings.** EA -1.3 points, from the T10 metric (1/(A^-1)11). EI_edge improves 3.8 points over
    9f72c6b (mostly T6). GJ -7.6 points, from the G3 webs (physical, see G3). All four S-2 LCs move toward BeamDyn;
    LC3 15.6 -> 11.3. Mass +259 kg.
  - EA in T9 is not T8w's +8.11 minus the T10 shift alone: the T10 metric change applies to every row (inference;
    the T10 standalone measurement is in its section).

**Re-audit at 61cf5b0** (`logs/check/input-audit-T9/`, mesh 33,935 / 35,120). Fail counts, original -> T9:

| r/R | F1 | F3 | F4 | F5 |
|---|---|---|---|---|
| 0.30 | 0 -> 0 | 72 -> 22 | 1 -> 0 | 0 -> 0 |
| 0.60 | 0 -> 0 | 17 -> 2 | 0 -> 0 | 2 -> 0 |
| 0.90 | 0 -> 0 | 44 -> 8 | 4 -> 0 | 2 -> 0 |
| 0.97 | 1 -> 0 | 8 -> 28 | 0 -> 0 | 2 -> 1 (axial -1.58 %; mass passes) |

F2, F6 and F7 are unchanged (pass or not gated).

F3/F5 diagnosis (explore, from the logs on disk plus hand calculation; per-branch audit runs not done):
- **Justified: 40 fails, the T6 mass-conserving partial-layer smear.** TE_PANEL/TE_REINF `Shell_skin(_inner)` at
  every station; the layer edge falls inside the region (`stackdb.py:207-224`).
  - The ply miss accounts for the whole total-thickness miss (0.90: 0.2297/3.073 = 7.47 %).
  - These are not "boundary" by S4's literal es/2 rule; they are justified only through the T6 "deliberate NuMAD
    choice" rule.
  - The 0.97 TE_PANEL fails were already in the original audit (-29.6/-30.8 %).
- **Candidate G10: spar-cap thickness is interpolated in span with pchip, every other layer linear**
  (`yaml_to_blade.py:380-383`, applied in `component.py:94-96`).
  - At 0.97 the cap row (z 113.481-113.607) carries 2.91815 mm. Hand pchip at the row centre gives 2.9181 mm;
    linear gives 3.0107 at the point and 2.9733 as the row mean.
  - Of the -3.07 % cap ply residual, about -1.2 % is the row mean against the point value, and about -1.9 % is the
    scheme. The F5 axial -1.58 % splits the same way (-0.85 % row mean, -0.74 % scheme).
  - **The original pass was a cancellation.** The bay mean sat +1.6 % above the point value and the pchip deficit
    was -1.2 %.
  - The fail first appears with the <= 0.25 m bays (T7 audit `bay0.25`: -2.93 %) and reaches T9 through T7e.
    T8/T8w/T11/T12 do not move it.
  - The yaml does not define the span scheme, so this is the same reference question as G1 (T13). It is not "fixed"
    toward linear without an established reference.
  - Not checked: other stations. Side note, not followed up: `yaml_to_blade.py:619` pchip-resamples every layer key
    onto the outer-shape z grid.
- **G10 reference (session 3, user L15).** Raw WISDEM sources, `wisdem/glue_code/gc_WT_InitModel.py`, fetched into
  `$SCRATCH/tmp/g10-ref/`:
  - **WISDEM never singles out the spar caps.** Every layer thickness uses one scheme. v3.0.0 and v3.2.0 use
    `np.interp(nd_span, grid, values, left=0, right=0)` (v3.2.0 `:277-283`). master uses
    `PchipInterpolator(grid, values, extrapolate=False)` (`:579-585`), since PR #419, merged 2023-03-16.
  - **The NuMAD rule matches neither version** (name contains "spar" -> pchip, else linear).
  - **The reference deck predates #419.** `tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat` is the 2020
    deck; the S-1 docstring says "PreComp/VABS/BECAS". A WISDEM-based chain of that time was linear
    (inference: the producing tool and its scheme are not established).
  - **The anchor cannot resolve G10.** BeamDyn has 26 stations, at eta 0.90, 0.95 and 1.00 near the tip, so the
    0.97 difference is below the anchor's resolution.
  - **Magnitude (worker, measured; `logs/check/g10/g10.py`, `g10.log`).**
    - Production equals pchip-of-yaml for `Spar_Cap_SS/PS` (deviation 0) and linear-of-yaml for every other layer
      (< 1e-6 relative).
    - The `update_internal_structure` pchip resample (`yaml_to_blade.py:~614-617`) changes no layer on this yaml,
      because every layer grid equals the outer-shape grid. It would act on a yaml with other grids (inference;
      relevant to T17/T18).
    - Mean cap thickness, pchip against linear: +0.033 % over r/R 0.2-0.85 (77.340 vs 77.314 mm, the S-1 window)
      and -0.37 % over 0.85-1.0.
    - Locally the difference reaches 34 % near the root (eta < 0.074, where the cap is thin) and a few percent at
      0.95-0.99.
    - So G10 cannot move the S-1 window numbers. It is a consistency and robustness question (the name heuristic),
      plus the 0.97 audit row.
- **Re-audit after G10** (df4c210, `logs/check/input-audit-G10/`). F5 passes at every station; the 0.97 axial goes
  from -1.58 % to a pass. F3 fails per station: 0.30: 22, 0.60: 2, 0.90: 8, 0.97: 28 -> 18.
  - The 0.97 cap totals now pass.
  - The 10 remaining cap fails are `Spar_Cap_SS/PS` ply -1.24 % (tolerance 1 %), on the A/B/C pieces of row 635.
    This is the row-mean share the diagnosis predicted. The laminate is the mean over the 0.126 m element row
    (T7e), and the audit cut sits 9 mm inside that row's inboard end, where the cap drops from 6 to 1 mm.
  - **Justified as the Rule 9 discretisation** (the laminate resolution is the element row), not a translation
    defect.
  - With the 40 T6 smear fails, every remaining F3 fail is justified. `check.sh quick` OK on df4c210.
- **Ruled out:** T12 ply order (per-ply checks run only when the order matches, `audit_blade_input.py:617-629`); the
  T6/T11 clamps (they move the edges, not the mechanism); audit-tool artefacts.

**Regression: CalculiX parity `tip_axial`.** `tests/test_blade_ccx_parity.py -k tip_axial` passes on T7e alone and on
T8w alone. On the combination it fails: tip 7.2 % off CCX (max 5 %), same for tension and compression. A worker is
bisecting ceca579 / 93df15e / 4c9372d.

**Measured, PROVISIONAL until the CCX regression is explained** (SLURM `submit_wait.sh T9`, worktree
`input-fix-T9` at 61cf5b0, es 0.25). Signed %:

| metric | 9f72c6b (`capture-s1s2v3`) | T7e | T9 |
|---|---|---|---|
| S-1 EA | +7.57 | +8.03 | +6.25 (T10 metric, about -1.8 points) |
| S-1 EI_edge | -7.51 | -4.18 | -3.71 |
| S-1 EI_flap | +3.79 | +4.61 | +3.89 |
| S-1 GJ | -2.76 | -1.26 | -10.36 (G3, physical) |
| S-1 mass | +0.73 | +1.58 | +1.65 |
| S-2 LC1 gravity | -7.44 | -6.56 | -5.59 |
| S-2 LC2 tip flap | -3.47 | -3.42 | -2.43 |
| S-2 LC3 tip edge | +15.57 | +12.11 | +11.31 |
| S-2 LC4 distributed | -6.27 | -6.33 | -5.47 |
| V-03 mass [kg] | 66,682 | 66,908 | 66,941 |

Failing tests go from 6 to 5: `test_flapwise_ei_within_band` now passes. The other five are the same as T7e (S-1 EA,
EI_edge, GJ bands; S-2 LC1/LC4 official).

## GJ gap and reference provenance (session 3, 2026-10-24, user L18)

With the yaml webs (G3) the S-1 GJ sits at a median of -10.2 % against BeamDyn over r/R 0.2-0.85. Diagnosis by a
worker; section-level calculations, no FE solve (`logs/check/gj-gap/`, `results.log`):
- **H1 metric: ruled out.** The deck's column order is x, y, z, rx, ry, rz; `K[5,5]` is GJ, as the loader reads it.
  K66 against 1/(K^-1)66 differs by -0.69 % (0.3R) and -0.99 % (0.6R). That is about 1 point, and in the wrong
  direction.
- **H2 provenance: mostly ruled out.**
  - The deck K data equals IEA-15-240-RWT commit `deef925e` (2021-11-29), not the 2020 release.
  - That commit's yaml (`$SCRATCH/tmp/deck-yaml/`) equals ours inside 0.2-0.85, webs included. It differs in three
    places only:
    - TE_reinforcement_SS/PS at the root, eta 0.04-0.10: the deck ramps from 0.04 (8.49 mm at 0.082 against our
      0.25 mm; upstream `e347d5eb`, 2022-11-06).
    - The tip chord grid (50 against 53 points, our G2 knots).
    - Gelcoat strength units.
  - Not tested: whether SONATA placed the webs where the yaml says (the cross-section images, upstream `827a69ee`).
- **H3 shell modelling: about 1 point.** The webs are uncoupled (A16 = A26 = 0, so 1/(A^-1)66 = A66). A web ending
  at the cap inner face gives -0.9 % (0.3R) and -1.8 % (0.6R).
- **Bredt-Batho against BeamDyn K66.**

  | webs | r/R 0.3 | r/R 0.6 |
  |---|---|---|
  | yaml | -3.6 % | -8.2 % |
  | spar-cap edge (B0) | +3.3 % | +2.2 % |
  | FE (T8w) | -0.5 % | -11.2 % |

- **Finding: the S-1 GJ per station is extraction noise on top of a real offset.** Shell/BeamDyn - 1 by r/R (T9):

  | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.85 |
  |---|---|---|---|---|---|---|---|
  | -41 % | +2.6 % | -9.5 % | -10.1 % | -8.9 % | -19.6 % | -3.3 % | -11.8 % |

  No smooth stiffness ratio does that. The real offset is the Bredt -4 to -8 %, plus SONATA's web placement (open).
- **Correction, measured (SLURM, `capture-GJM.pkl` against `capture-T14m.pkl`, both at ad5f3ac): the per-station
  scatter is mostly real profile, not extraction noise.** The window-integral estimator (`fix/s1-gj-metric` 2e79750,
  synthetic evidence only, no true RED) leaves median and spread about unchanged.

  | estimator | median | sd |
  |---|---|---|
  | spline | -10.2 % | 7.3 |
  | window | -9.8 % | 7.2 |

  Per station:

  | r/R | 0.21 | 0.26 | 0.30 | 0.34-0.58 | 0.62 | 0.66-0.74 | 0.79 | 0.83 |
  |---|---|---|---|---|---|---|---|---|
  | spline | -32 | -0 | +3 | -2 to -12 | -12 | -17 to -20 | -10 | -3 |
  | window | -33 | -17 | -9 | -7 to -10 | -14 | -18 to -20 | -9 | -17 |

  - The two disagree only at the window edges (0.26-0.30, 0.83). At 0.30 the spline is closer to Bredt (-3.6 %).
  - **Not merged: a null result.** The S-1 signed mean would move -10.3 -> -13.2 % from edge stations alone.
  - **Finding:** a GJ dip of -17 to -20 % at r/R 0.66-0.74 and -32 % at 0.21, present in both estimators. It needs its
    own diagnosis: what changes in the layup, the webs or the reference there.
- **GJ dip diagnosis (worker; `logs/check/gj-dip/REPORT.md`; S-1 torque replicate SLURM 11623290).**
  - **About 75 % of the 0.66-0.74 dip is section against reference, not FE.** Bredt (yaml webs) against BD drifts
    smoothly: -5 % at 0.45, -15 % from 0.70 to 0.85. No layup, airfoil or element-set event lies in 0.55-0.80, and
    the BD K66/K44/K55 are smooth.
  - **Bredt with the webs at the cap edges (B0 placement) stays within +4/-4 % along the span.** That closes about 10
    of the 15 points (indicative: B0 is a different mesh).
    - Inference: the BD deck's section model placed the webs wider than the yaml.
    - Pending proof: SONATA cross-section images, upstream `827a69ee`.
  - **About 25 % is FE below its own Bredt** (-5 to -7 % at 0.66-0.74). The jump to -3 % at 0.81-0.83 is FE-side:
    the twist is non-monotone at 0.817, where the TE reinforcement, the fillers and the LE reinforcement all end.
  - **Root 0.19-0.24: a reference-interpolation artefact.** BD K66 drops 2.7x between deck stations 0.20 and 0.25,
    and the linear interpolation inflates the reference +8-12 %. Against a log-interpolated reference, Bredt sits at
    0 to -5 %. The rest of the -32 % is FE below Bredt, unresolved.
- **Next (cheapest first).**
  1. Log-interpolate `gj_bd` near the root (a metric check).
  2. Full S-1 on a cap-edge-web mesh: is FE/BD flat?
  3. Inspect the rings at 0.81-0.835 offline (`logs/check/gj-dip/mesh_ad5f3ac.pkl`, `solve-ad5f3ac.pkl`).
  4. The root FE deficit (es 0.125 or a thick-wall Bredt).
  5. The user decides the reference policy for GJ (see L19).
- **SONATA web placement, measured: the premise of L19 is refuted** (worker; `logs/check/sonata-webs/REPORT.md`; vector
  extraction of the upstream 827a69ee section PDFs at r/R 0.2-0.8).
  - SONATA put the webs at the yaml `start/end_nd_arc` (mean +0.003, max 0.010 nd_arc), 0.69-0.71 m apart
    (2 x `offset_y_pa`).
  - The CarbonUD caps (0.89-0.91 m chordwise) match the yaml cap arcs. Each web sits 0.09-0.10 m of arc inside the
    nearest cap end.
  - The parent's visual reading of the 0.6 plot ("webs at the cap ends") was wrong: the web meets the cap taper.
  - The 827a69ee yaml and the deef925e yaml differ by <= 0.003 nd_arc on 0.2-0.8. deef925e equals ours.
  - **Consequences.**
    - Our yaml-faithful web geometry is the reference's geometry; G3 stands.
    - The B0 (cap-edge webs) agreement with BeamDyn GJ was a cancellation.
    - **The S-1 GJ gap (-5 % at 0.45 to -15 % from 0.70, Bredt and FE alike) is real and unexplained by web placement.**
      Open hypotheses (inference): web/cap/skin laminate in-plane shear (A66, the biax/triax G12), and how the
      reference section tool (SONATA/ANBA4: warping, thick walls) computes torsion against the thin-wall shell.
  - The GJ anchor stays red as an open finding. No band touched; L19's "reference built with another web placement"
    is withdrawn.
  - Measure the deck-era yaml as input (worktree DECKY, `submit_wait.sh DECKY`) against the same tree with the
    official yaml (T14m, ad5f3ac). The root TE reinforcement can move S-2 even though S-1 cannot see it.

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
- L14 (2026-10-23, session 3) User, on "next step: T9 merge, better in a new session": "METELE".
- L15 (2026-10-23, session 3) User asked "que me recomiendas?" on G10. The parent recommended: do not change it now;
  measure pchip vs linear along the span, read the WISDEM scheme, then decide (a uniform reference scheme ->
  adopt it for every layer and drop the name heuristic; none established -> record as G1, replace the heuristic by an
  explicit parameter keeping today's default). User: "si hazlo".
- L16 (2026-10-23, session 3) User, on why WISDEM moved to pchip (answer: PR #419, optimizer smoothness on resampled control nodes; windIO prescribes no span scheme for layers): "me va asi, vamos a hacer esto, ponemos lineal, y dejamos el pchip como no implementado pero planificado, asi podemos avanzar mas rapido, y lo implementamos a demanda cuando haga falta si hiciera falta"
- L17 (2026-10-24, session 3) User asked "por fin corregiste todo a lineal?" The parent answered: thickness only; still pchip are layer nd_arc bounds in span (`stackdb.py:196-198`, should be linear for consistency), outer geometry (`geometry.py`, pchip correct per windIO/WISDEM) and the in-bay surface loft (`mesh_gen.py`); recommended G10b for the first. User: "si".
- L18 (2026-10-24, session 3) User: "siempre que tengas trabajo que hacer que no este bloqueado por algun job tienes permiso para hacerlo. el objetivo es tener un modelo fisicamente realista. entonces cualquier cosa que nos lleve ahi hayque hacerla"
- L19 (2026-10-24, session 3) GJ reference policy. Parent recommendation: keep the webs where the yaml puts them (the design); record the S-1 GJ anchor as a reference built with another web placement, no band touched; confirm with the SONATA cross-section images before closing. User: "continuamos segun tus recomendaciones".
- L20 (2026-10-24, session 4) User: "Abrí una sesión nueva y arrancá con mem_search \"blade-input-audit\" y la sección Handoff 4 de odd/tasks/blade-input-audit.md (commit d8a1c3a). El orden: 1. Re-sellar todas las filas del store que dependen de la pala (S-0, S-3, S-4, S-5/S-6, S-7, V-01/V-02/V-05, paridad con CalculiX) con la suite 11623302. Hoy están desactualizadas, porque la geometría cambió. 2. T21: si la métrica de trabajo pasa el 5 %, commitear y mergear; si no, registrarlo como hallazgo. 3. DECKY: ver si el refuerzo de TE en la raíz del yaml del deck mueve LC3. 4. La brecha de GJ: comparar la rigidez al corte del laminado de una sección contra SONATA/ANBA, y probar una variante de pared gruesa. 5. Las masas abiertas de T24: vidrio en la cara plana del TE, faltantes de triax y carbono, y almas 4.8 % por debajo del yaml. 6. La segunda iteración de mallado: T16, G6 y el warning de element_utils.py:78." Later: "utiliza workers en background para optimizar el contexto del orquestador". Evidence: Engram unavailable (server ownership mismatch at 127.0.0.1:7437); suite 11623302 TIMEOUT 4 h with an OOM kill, no summary, so item 1 is re-measured (T25).
- L21 (2026-10-24, session 4) User: "no hay algo que podamos hacer en paralelo mientras esperamos?" Parent: T29-T31, none moves the geometry of l3/twist-p1 (T30 lives in its own worktree), so the T25 restamp stays valid; standing permission L18.
- L22 (2026-10-24, session 5) User: "leé el Handoff 5 y las líneas T21 y T24–T31", then "si" to start with the T31 verdict and the 11623429 failure.
- L23 (2026-10-25, session 6) User: "leé el Handoff 6 y las líneas T25 y T31–T35 de odd/tasks/blade-input-audit.md. Ese archivo es el único registro: Engram sigue caído por el mismo conflicto de servidor en 127.0.0.1:7437, y falló el intento de guardar el resumen de sesión." then "no se que dices que esta corriendo paero ya no hay nada corriendo. mezcla todo lo que entiendas que sea necesario y continua la validacion". Parent: merged T32 and T33 only (T30 waits for the T24(a) decision), full --record 11623853, T34 station probe.
- L24 (2026-10-25, session 6) User: "no me digas esas cosas y dime si mejoro la fisica del problema q es lo que me importa", then "si quiero corregir la fisica". Parent: this session improved the evidence, not the model physics; next is T37 (flat-back closure per SONATA).
- L25 (2026-10-25, session 7) User: copy and commit the 11623882 known-red record, compare with T25; then T37: close the flat-back as SONATA/ANBA4 does (bond line ~4 % of arc ahead of the TE + solid adhesive wedge outboard of 20 %, TE glass_uni off the face; inboard skins + gelcoat + glass_uni on the face); confirm on the 0.3/0.5 PDFs first; propose web vs equivalent face laminate and go on if clear; write the prediction before running; measure at es 0.25 S-1 GJ/EI_edge/EA, S-7, V-02 1F/1E, S-4 1E/2F, V-03 mass; never widen. Parent: 11623882 done, 21 failed = the 28 of T25 minus T32 (2), T33 (2), T34 (3), committed 9b5de65.
