# L3 closure: blade structure, sensitivity first

Step 2 of `odd/tasks/thesis-validation-plan.md` (L3, IEA 15 MW official blade). Store `docs/validation/`, tool
`tools/validation_matrix.py`. Model history and diagnosis of the reds: `odd/tasks/blade-input-audit.md` Handoff 10.

## Specs

- S1 (user, verbatim): "vamos con L3"
- S2 (user decision, option chosen verbatim): "1. Sensibilidad primero (recomendado)" — "Re-medir las stale y auditar las pendientes. Para cada red, medir su efecto en torque/potencia con BEM rígido y twist perturbado (barato, no depende de L5). Si el efecto es chico, se declara como limitación. Si es grande, recién ahí se toca el modelo (cierre TE fiel, T37)."
- S3 (thesis plan step 2 acceptance): 5 %; every remaining red gets its torque/power sensitivity measured and is declared as a limitation (EA: no torque path, declared).
- S4 (repository physics discipline): no pass band widened; a moved anchor is a finding; `regression --write` never over a moved residual without reporting it.
- S6 (user L4, verbatim): "si hay opcion de mejorar el modelo hay que mejorarlo. Las formas modales son muy importantes porque luego voy a necesitar hacer analisis de frecuencias en rotacion para ver como esas viraciones modifican o no la se~nal de torque. una pala 10% mas blanda (no se contra que referencia) hace que sus oscilaciones sean diferente 10% diferentes a la respuesta fisica"
- S7 (user L5/L6): second blade from the official repositories to tell model defects from solver-core defects; option chosen verbatim: "1. IEA-3.4 antes de T6 (recomendado)".
- S8 (user choice for core-torsion C5, verbatim): "Cajón para validar L3; el BEM se decide en L5 (recomendado)" — S-7, S-1 GJ and twist vs Zhou are measured with the box rigid fit; the BEM feedback twist is decided at L5 with the C5 numbers.
- S9 (user choice, verbatim): "Importar ambos; Campbell contra BeamDyn p20 (recomendado)" — `sx_campbell` re-pointed to BeamDyn p20 per rpm, `sx_mesh_convergence` with its artifact regenerated at HEAD.
- S10 (user, verbatim): "recuerda que puedes usar las dos palas ademas de ejemplos teoricos para comprobar los analisi hechos" — every L3 analysis is checked on IEA-15 and IEA-3.4 plus a closed-form case.
- S11 (user choice, verbatim): "Arreglar la introducción de carga primero (recomendado)" — "Tracción de borde equivalente en trabajo, ponderada por longitud de arco, y media del anillo también por arco. Primero se valida en el tubo que converja al refinar; después se re-corren las dos palas y recién ahí van los veredictos."
- S5 Done means: `contract --level L3` closed, or every open entry listed with its cause and, for reds on a torque path, its measured torque/power sensitivity.

## Starting state (2026-10-26, HEAD 332c69f, `contract --level L3 --rows`, `logs/check/l1l2/contract-L3-start.txt`)

L3 [open]: undecided 10, stale 15, unaudited 9, absent 2; gaps 1 blocking (`rotating_modal_shift`).
Stale rows carry `audited_rev 9f72c6b` (before the cap-arc / nd_arc geometry commits e6327ec, 028d717 and T40
d81931b/a4b3ca2). Known reds (blade-input-audit Handoff 10 / T25): S-7 1.0717, S-1 GJ 10.74 %, S-1 EI_edge 8.01 %,
S-1 EA 6.68 %, S-2 LC1 uy 5.26 %, LC4 uy 5.48 %, V-02 1F/1E reference conflict. In: S-4 (1F/1E/2F within 2.4 %),
masses, rated tip twist 1.044 x Zhou.

Torque paths: twist (S-7, S-1 GJ, rated twist vs Zhou); edgewise (S-1 EI_edge, S-4 1E); flapwise (S-1 EI_flap,
S-2, S-4 1F); EA none (declare); mass (V-03, V-05).

## Tasks

- T1 | S2 S4 | parent, SLURM | done (store pass 1, 2026-10-27; see "Store pass 1 result") | Re-measure stale groups at HEAD (`regression --group G --scope <row nodes>`, no `--write`): 38 S-1, 34 S-7, 39 S-2, 40 S-4, 41 V-03, 42 V-05, 32 V-02, 20 anchor beam; restamp unchanged rows, report moved ones.
- T2 | S2 S4 | parent | pending | Audit group 23 (CCX parity, Bernardi, Escalera) and group 5 (mesh convergence) rows; decide S-2 LC2/LC3 undecided.
- T3 | S5 | parent | done: gap bounded with both S-4 rows; shift delta-f 7-13 % below BeamDyn noted | Gap `rotating_modal_shift`: its cause was "no row"; group 40 rows exist, link them and set the status.
- T4 | S2 S3 | parent | done: twist +0.34 %, flap -0.34 %, edge 0.005 % torque | Sensitivity tool: rigid BEM at rated with the twist distribution perturbed by each red's effect; torque/power deltas.
- T9 | S7 S4 | parent/worker | in progress: (a) d7d000e, (b) 188c997, (c) S-1 + modes measured 94237f7; modes refuted; H6 held (secondary), audit fix b9192df, yaml Bredt 0630e71; mesh Bredt 863ff38 (mesh = yaml), H7 refuted, FE torsionally soft 0.65-0.85 (twist integral); element vs distortion open (Handoff 4) | Second blade IEA-3.4-130 before T6 (T17 of blade-input-audit; v1 yaml at IEAWindSystems commit 824f8d6, 2 webs, glass caps, SONATA 6x6 = same chain as the 15; ANBA mass structural only, issue #5): (a) pre-register the expected signed biases (EA, EI_flap, EI_edge, GJ outboard, 1F/1E/2F) from the IEA-15 diagnosis; (b) ingest + input/mesh audit (clamps the 15 never binds); (c) S-1 section vs its SONATA 6x6 and S-4 modes vs its OpenFAST deck; (d) T6/T7 must then improve both blades. BECAS blades (10, 22: 3 webs) out of this task.
- T6 | S6 S4 | parent/worker | pending | Cap chamfer (SONATA convention, blade-input-audit Handoff 8 item 1): write the prediction (EA, EI_flap, S-2 LC1/LC4, 1F) first, section check vs the SONATA PDF, then the FE anchors. Target: flap -5.5 % and EA +6-9 % reds.
- T7 | S6 S4 | parent/worker | pending | Faithful T37 TE closure (bond-line keypoint + solid closure; Handoff 9 item 3, predicted +5..7 pt GJ at 0.7-0.9): prediction, section check, then S-7 / S-1 GJ / torsional mode. One rule at a time, after T6.
- T8 | S6 | parent | pending | Mode-shape rows: MAC per mode (1F 1E 2F 2E 1T, parked and rated) of our modes vs BeamDyn p20 mode shapes; frequencies alone do not pin the shapes the rotating torque analysis needs.
- T10 | S10 S4 L12(4) | parent, SLURM | in progress (branch l3/pass2-s2-ring) | S-2 tip-ring-mean metric: (a) closed-form tube, ring mean vs Timoshenko and node-vs-ring spread under the S-2 load introduction; (b) node / ring mean on both blades at es 0.25 / 0.125, IEA-3.4 BeamDyn static LC1-LC4 as its reference; (c) test asserts the ring mean LC1-LC4 vs BeamDyn p60, LC2/LC3 verdicts, store rows.
- T11 | S9 | parent | pending (own session, after the full suite 11625665) | Import `sx_mesh_convergence` (artifact regenerated at HEAD) and `sx_campbell` re-pointed to BeamDyn p20 per rpm into the L3 store.
- T5 | S3 S5 | parent | pending | Declare the remaining reds with cause and sensitivity; `contract --level L3`; decide model work (T37) only if a sensitivity is large.

## Session 10 notes

- T1 jobs: `regression` for groups 23+5 (11625012), 34 (in 11625009), 32+20 (in 11625011); scopes in
  `logs/check/l3/scope-<g>.txt`, outputs `logs/check/l3/reg-<g>.json|.err`. Groups 38-42 declare no residual
  pattern (`regression` exits 2: "declares no residual patterns"), so they run as plain pytest on the row nodes
  (`logs/check/l3/pyt.sbatch`): 40 S-4 11625023, 38 S-1 11625024, 39+41+42 11625025; outputs `pyt-<g>.txt`,
  to be compared by hand with the values stored in `docs/validation/rows/<g>-*.yaml`.
- T4 design change: the production two-way coupling feeds back to the BEM only the deformed strip radius and the
  elastic twist (`fsi_participant.py` `_compute_deformed_geometry`; no precurve/presweep), so the sensitivity
  perturbs those two channels. Tool `tools/l3_torque_sensitivity.py` (one rated shell solve under the rigid
  production loads, as `test_bem_fsi_deformed_geometry.py::fsi`; then rigid-geometry BEM per correction).
  Job 11625022 (dev) with provisional inputs: twist ratio 1.0717 (S-7), flap error +5.26 % (S-2 LC1, sign to be
  confirmed by T1), edge error +8.01 % (S-1 EI_edge mid-span). Known from the existing test (one-way): twist alone
  -0.81 % power, deformed radii +0.85 % power (prebend straightening).

- T1 partial (332c69f): S-1 38: EI_flap, mass pass; EI_edge, GJ, EA fail (known reds). S-2 39 official, signed vs
  BeamDyn p60: LC1 uy -5.26 % (fail), LC2 uy -2.5 % (pass), LC3 ux +10.26 % (edgewise ratio test passes), LC4 uy
  -5.48 % (fail): our blade is STIFFER in flap and softer edgewise. V-03 41, V-05 42 pass. S-7 34: regression
  `changed` 7.2496 -> 7.1712 % (the documented T40 move, 1.0725 -> 1.0717, not new).
- T4b first run (11625022, provisional flap sign +5.26 %, wrong): rated one-way, tip elastic twist +3.758 deg
  (BEM sense, toward feather), dr -0.814..0 m; deformed vs rigid torque -4.40 %, thrust -9.09 %; twist alone
  -3.18 % torque, radius alone -1.36 %. Removing the S-7 twist excess (/1.0717): torque +0.34 % = 7.4 % of the
  flexible-minus-rigid torque delta. Edge channel (+8 %): +0.004 % (negligible). Flap re-run with the measured
  sign (-5.48 %, LC4) and edge +10.26 %: 11625034. Observation for L5 (not chased here): the production radius
  only shortens (dr <= 0) on a blade prebent upwind; the deformed-geometry test docstring says the radii grow.

- T4b signed run (11625034): flap error -5.48 % (LC4, ours stiffer), edge +10.26 % (LC3). Removing each red from its
  channel, rated one-way, torque (= power at fixed rpm) vs the deformed blade / as a fraction of the
  flexible-minus-rigid torque delta (-4.40 % of rigid): twist (S-7 /1.0717) +0.340 % / -7.4 %; flap (radius with
  u_y x 1/(1-0.0548)) -0.343 % / +7.5 %; edge +0.005 % / -0.1 %. The two torque-path reds cancel to ~0.00 % (sum of
  the separate cases, the combined case not run); the cancellation is coincidental, so each is declared on its own.
  Thrust: twist +0.59 %, flap -0.33 %. JSON `logs/check/l3/sensitivity-signed.json`.
- T1 done: S-4 40 all pass vs BeamDyn p20, values identical to T25 (parked 1F +2.33, 1E +0.56, 2F +0.20 %; rated 1F
  +1.85, 1E +0.44 %); the ElastoDyn MBC3 comparison (not asserted) reads -4.4..-7.8 %. Group 23: 11 same, 5 unmapped
  (escalera static deflection, bernardi 4-7: no residual printed). Group 5: unmapped, test exit 1 (the mode-6 gap
  finding of T25). Groups 32/20: see `reg-32.json`, `reg-20.json`.

## Handoff 1 (2026-10-26, session 10 close, ~175 k context)

**Resume.** Read this file only (Specs, Session 10 notes); logs in `logs/check/l3/`. Engram is down.
**Decision taken (user L4): improve the model (T6, T7) and add mode shapes (T8); declaration only for what stays red after them.** Superseded: threshold for "small" sensitivity. Parent recommendation: declare both reds (S-7 torsion,
S-2 flap) with no model work: each moves rated torque by 0.34 % (one order below the 5 % acceptance) and they
cancel; the 7.4 % share of the flexible-minus-rigid delta is stated in the declaration so L5 can carry it.
**Next, in order.** 0. Session 11: T9 (IEA-3.4: pre-registered biases first), then T6, T7 on both blades, then T8; items 1-3 below can run alongside. 1. Restamp the stale rows whose values did not move (S-1 mass/EI_flap, S-4 x2, V-03, V-05,
anchor beam 20, V-02 32 if `same`), record moved ones as findings (S-7 7.25 -> 7.17 %, T40). 2. T3: group 40 rows
exist -> `rotating_modal_shift` gap to `bounded` with the rows listed. 3. T2: group 23 unmapped rows and group 5
need a printed residual or a manual verdict; S-2 LC2 (-2.5 %, in) and LC3 (+10.26 %, the edgewise ratio test
passes with a different band: read it before deciding). 4. T5: declarations with the sensitivity numbers above
(EA: no torque path), then `contract --level L3`.

## T9 IEA-3.4-130 (session 11, 2026-10-26)

**Inputs and provenance (checked before any run).** Clone `$SCRATCH/IEA-3.4-130-RWT` (HEAD aacff39, windIO 2.1).
- yaml v1 = `824f8d6:yaml/IEA-3.4-130-RWT.yaml` (2025-09-26); BeamDyn deck `openfast/IEA-3.4-130-RWT_BeamDyn_Blade.dat`,
  21 stations, last touched in a42ea70 (2020-11-09, "BD input files, webs thick min 20 mm"); that commit also adds
  `yaml/figures/blade_section_*.png` (SONATA-style section plots, like the IEA-15 PDFs) and edits the yaml.
- The v2.1 yaml at HEAD embeds `structure/elastic_properties` = the same deck (K33, K66, mass 0.0 % at 0.1-0.9), not
  a recomputation; the `six_x_six` in a42ea70 is a 2-point placeholder. The deck is the only 6x6 of this blade.
- **The pair is consistent:** layer-by-name diff a42ea70 -> 824f8d6 finds no physical change in the blade (layup,
  webs, materials, outer shape): only `TE_reinforcement` (one ud layer, width 0.8 m, midpoint TE) split into
  `_SS`/`_PS` (0.4 m each, start/end at TE; equivalent), cap names capitalised, airfoil digits, paint fvf/fwf
  (cost only). An index-based diff first showed x5 thickness changes: an artefact of the inserted layer.
- Issue #5 ("ANBA mass structural only") concerns the TUM spreadsheet (ANBA) vs the ElastoDyn file, not the deck;
  the windIO model was a NREL conversion at the blade-structure level. ElastoDyn (PreComp, WISDEM) vs the deck:
  mass -3..+5 % mid-span, EI_flap +2.5..+26 %, EI_edge +6..+16 %: PreComp is not an arbiter here.
- Trap found: `tools/beam_reference.load_beamdyn_blade`/`load_elastodyn_blade` default `blade_length=117.0`
  (IEA-15); every 3.4 call has to pass its axis length (63.0 m in the yaml reference axis; deck eta is arc length).

**T9a pre-registration (written before any IEA-3.4 shell run).** Signed %, ours vs the deck, S-1 metric (mid-span
0.2-0.85, mean signed; mean |rel| in brackets for the 15), es 0.25, mid-surface skin, same code as the 15 at 1902eed.
Scaling rule: the cap-chamfer convention (blade-input-audit Handoff 8 item 1) removes ~t^2 per cap end pair, so
`dEA/EA ~ 2 E_cap t^2 / EA_deck`. Over eta 0.2-0.8 it gives **IEA-15 +7.8 %** (measured EA +6.56) and **IEA-3.4
+4.8 %** (glass ud E1 42 GPa, t 30-63 mm, w 0.8 m, cap EA share 68-94 %): ratio r = 0.62 for chamfer-driven metrics.

| metric | IEA-15 at HEAD | mechanism (15 diagnosis) | IEA-3.4 prediction | refuted if |
|---|---|---|---|---|
| EA(CLT) | +6.56 (6.68) | cap chamfer convention | **+4.8, band +3.3..+6.3** | >= +6.5 (systematic, not chamfer) or <= +2 |
| EI_flap | +3.8 (4.53; signed from 11622464, HEAD prints only mean abs) | chamfer + cap position | **+2.4 (r x 3.8), band 0..+4.5** | >= +5 or < -1 |
| EI_edge | -1.89 (8.01) | TE solid closure absent (T37 sonata_wedge +6.9 pt) | **-2, band -7..+3; mean abs scatter 5-9** | > +4 or < -10 |
| GJ 0.3-0.6 / 0.6-0.85 | -5.6 / -13.0 | TE solid region outboard (Bredt per wall) | **negative in both bands, outboard worse; -2..-10 / -5..-18** | either band >= -1 |
| mass | +2.26 (3.54) | foam/TE bands | **-3..+3** | outside +-4 |
| S-4 parked 1F / 1E / 2F vs BeamDyn modal | +2.33 / +0.56 / +0.20 | follows EI_flap, mass | **+1.4 / +0.5 / 0, each +-1.5** | 1F >= +3 or <= -1; 1E outside +-2 |

What each outcome discriminates. (i) EA and EI_flap scale with r (both within band): the flap/EA reds are the
chamfer convention, turbine-specific, and T6 applies to both blades. (ii) EA at ~ +6.5 % or more on a glass-cap blade:
a systematic bias of our section (mesher/CLT), T6 alone will not close it. (iii) GJ negative with the same outboard
shape: the TE closure (T7) is a modelling convention shared by both SONATA decks; GJ near 0 on the 3.4: the 15's
deficit is turbine-specific (its TE band geometry). (iv) Any metric off by the same % on both blades regardless of
r: suspect the solver core / extraction (but L1/L2 closed and S-4 within 2.4 % bound that). (v) Mesher or reader
clamps, guards or warnings that fire on the 3.4 and never on the 15: T9b findings, listed before any anchor.
Mesh: the 3.4 chord is 0.68-0.76 of the 15 at equal eta, so es 0.25 is coarser per chord; run 0.25 and 0.15 and
trust 0.15 if any S-1 metric moves > 1 pt.

**T9b ingest (188c997).** Vendored `tests/IEA-3.4-130-RWT.yaml` (824f8d6) and `tests/IEA34/reference/` (BeamDyn
blade + driver deck, 824f8d6 = a42ea70 bytes; Apache 2.0). Reader clamps the IEA-15 never binds, each fixed with a
RED test (`tests/software/mesh/test_yaml_materials_optional.py`, 4 passed; `check.sh quick` 75 passed):
(1) `Xt/Xc/S` required on every material (3.4 `resin`, `steel` lack them; unused by the blade) -> NaN;
(2) isotropic `G` required (`resin`, `steel`) -> E/(2(1+nu)); (3) `fiber_orientation` required on every layer
(absent on most 3.4 layers; the yaml path never consumes `fabricangle`, stackdb sets ply angle 0) -> 0 deg;
(4) exactly two `shell` components (15: outer + inner skin; 3.4: one) -> one or two.
Checked harmless: missing `ply_t` (triax, biax, ud) falls back to 1 mm, but `nPlies` is a float
(`continuous_layers`), so thickness is exact. Still name-driven (substring matching in `_add_components`): caps
`spar`+`ss`/`ps`, webs by `fore`/`rear` in the layer's `web` key (3.4: `Web_aft_*` belong to `fore_web` -> group 1,
keypoint b; `Web_rear_*` -> group 2, c); the input audit (mesh vs yaml) is the check that the mapping is right.
Input audit es 0.25 (`logs/check/l3/t9/audit-es025/`, 13592 nodes): 0.6 clean; 0.3: 42 fails (Shell_skin ply
1.41 vs 1.44 mm, chord -1.5 %, one cap-edge element where the yaml has balsa); 0.9 / 0.97: the ud layer missing in
the LE/TE panels, axial sum(Ex t ds) -26 % / -72 %, mass -54 %, 0 webs vs 2 at 0.97. **Cause: partial-span layers,
which the IEA-15 does not have (all its layers span 0-1).** 3.4: caps and webs on grid 0.1-0.94, LE/TE reinforcement
and all fillers on 0.1-0.8, each ending at a nonzero thickness (fillers 7-13 mm at 0.8). The mesh drops a layer
outside its grid (WISDEM convention, thickness 0); the audit's yaml side uses `np.interp`, which holds the end value:
the outboard fails are an audit artefact (fix: zero outside the layer grid), to confirm against the mesher rule at
the grid ends (abrupt or over one bay). Root 0-0.1: skin only (triax 65 -> 20 mm, paint), deck EA 11.3 -> 4.7 GN.

**T9c-2 result (worker, not committed yet; `tools/l3_s4_blade.py`, `tools/bd_rotating_freevib.py` beam `iea34`).**
Parked, Hz: ours es 0.25 0.7369 / 0.9485 / 2.1929, es 0.15 0.7364 / 0.9470 / 2.1936; BeamDyn iea34 p20 0.7565 /
0.9869 / 2.2703 (p15-p30 within 0.5 %; load and dt halved move <= 0.001 Hz; fit residual 0.25 vs 0.038 on the 15,
from unfitted modes above 8 Hz, a 14 Hz band leaves the three unchanged). **Signed vs BeamDyn (es 0.15): 1F -2.66,
1E -4.04, 2F -3.38 %. Pre-registration refuted** (+1.4 / +0.5 / 0): every mode is low, opposite to the 15. Same
code on the IEA-15 reproduces +2.33 / +0.56 / +0.20 % (job 11625095). Mode labels: shell beam-fraction 1.000
(1F, 1E, 2F, 2E, 3F), BeamDyn by dominant tip channel; no swap. Jobs: BD 11625092/113/147/173/174/214, shell
11625093/94/95.

**T9c-1 S-1 result (`tools/l3_s1_blade.py`, c2db0be; jobs IEA-15 11625089, 3.4 es 0.25 11625230, es 0.15 11625231).**
Same code on the IEA-15 reproduces the stored S-1 exactly (EI_edge -1.89 (8.01), GJ -10.74, EA +6.56 (6.68), mass
+2.26 (3.54)) and gives the signed EI_flap at HEAD: **+3.87 (4.53)**. IEA-3.4, signed mid-span (mean |rel|):

| metric | predicted | es 0.25 | es 0.15 | verdict |
|---|---|---|---|---|
| EA(CLT) | +4.8 (+3.3..+6.3) | +4.31 (4.33) | +4.38 (4.39) | **held** -> chamfer scaling r confirmed, outcome (i) |
| EI_flap | +2.4 (0..+4.5) | +0.76 (3.59) | +1.45 (3.48) | **held** |
| EI_edge | -2 (-7..+3), scatter 5-9 | -1.37 (4.70) | -1.82 (4.72) | **held** (scatter 4.7, just under) |
| mass | -3..+3 | +1.43 (2.34) | +1.49 (2.38) | **held** |
| GJ 0.3-0.6 / 0.6-0.85 | -2..-10 / -5..-18 | -13.9 / -31.8 | -13.4 / -32.2 | sign held, **magnitude refuted** (0.2-0.3: -60.6) |

The mid-span sections cannot explain the modes (EI_flap +1.5 and mass +1.5 give 1F ~ 0, measured -2.7). Outside
0.2-0.85 (es 0.15, signed % per station): EI_flap -19..-31 at 0.04-0.11 and EI_edge -18..-26 at 0.08-0.18 with EA
and mass within -5..+3; EI_edge -21 at 0.81-0.83; EA +48..+60 and mass +34..+43 at 0.93-0.95; EI -38..-100 at
0.96-1.0 (free-end curvature, not a section value). Same code on the IEA-15 root: EI_flap -2..-6, EI_edge -4..-14:
the 3.4 root deficit is not a generic extraction artefact.

**Post-hoc hypothesis H6 (not pre-registered; to be tested with a written prediction).** The 3.4 layers start and
end abruptly at their grid ends (caps 9.9 mm at 0.1 and 2 mm at 0.94; LE/TE reinforcement and fillers 1-17 mm at 0.1
and 7-13 mm at 0.8; webs at 0.1/0.94). The deck has 21 stations every 0.05 eta and BeamDyn interpolates linearly, so
the beam smears every start/end over 3.15 m (0.05-0.10, 0.80-0.85, 0.90-0.95), where it is stiffer and heavier than
the shell step. The shell sees the step at the grid end. That matches the station pattern above and the sign of
all three modes; the IEA-15 has no partial-span layer, so it never showed it. Not yet separated from: shear lag at
the abrupt cap start (real shell physics a beam lacks), and the GJ deficit (-29 % mean; torsion couples little
into 1F/1E).

**H6 pre-registration (session 12, written before any ramped run).** Correction to the premise: the shell does
not see a step. `update_internal_structure` (`yaml_to_blade.py:602`) resamples every layer onto the 50-point
reference-axis grid with zero outside its own grid, and the laminate is linear between grid stations, so the shell
already ramps over one grid interval: 0.0833-0.10, 0.80-0.8175, 0.94-0.955 (~1 m). The deck ramps over 0.05-0.10,
0.80-0.85, 0.90-0.95 (3.15 m). Prior from a decoupled Euler-Bernoulli beam on the deck scaled by the S-1 station
ratios (`tools/l3_beam_station_ratio.py`, es 0.15 json): all ratios <= 0.955 give 1F/1E/2F -3.4 / -5.4 /
-3.8 % (measured -2.66 / -4.04 / -3.38): the station deficits do explain the modes, mostly the root 0-0.20
(-2.8 / -3.6 / -2.3). Compared at the deck stations themselves (ours interpolated, no deck interpolation),
EI_flap is -20.9 at 0.05 and EI_edge -20.6 at 0.15, where no layer ends in either model: H6 cannot close those.
H6 bands set to ratio 1 (upper bound): +2.4 / +1.7 / +2.6; bands bridged between their neighbours (realistic):
+0.8 / +0.6 / +1.6. Ramped yaml (`tools/l3_h6_ramp_yaml.py` -> `logs/check/l3/t9/h6/`, never in tests/): each partial-span layer and web
becomes the deck's linear interpolation across each deck interval that contains its grid end (zero at 0.05,
yaml value at 0.10; yaml value at 0.80, zero at 0.85; t(0.90) to zero at 0.95), other keys held at their end.

| quantity (es 0.15, signed vs BeamDyn) | now | H6 predicts | refuted if |
|---|---|---|---|
| 1F / 1E / 2F shift | | **+0.8 / +0.6 / +1.6** | any shift < 0, or above +2.9 / +2.2 / +3.1 (bound + 0.5) |
| 1F / 1E / 2F after | -2.66 / -4.04 / -3.38 | -1.9 / -3.5 / -1.8 | (follows) |
| S-1 EI_flap, EI_edge at 0.0583-0.0917 | -19..-31 | rise toward the 0.0417 / 0.11 neighbours | no rise |
| S-1 at 0.0417 and EI_edge at 0.11-0.18 | -19, -18..-26 | unchanged within 2 pt (no layer end) | moves > 2 pt |
| S-1 EI_edge 0.81-0.83 | -21 | closes to the 0.75-0.78 level | stays below -10 |
| S-1 mass / EA 0.93-0.95 | +34..+43 / +48..+60 | shrink (less cap and web past 0.90) | grow |

If H6 holds as predicted it is real but secondary, and the modes stay 1.5-3.5 % low: the lever is the root
section 0.04-0.18, where mass and EA agree (+0.4 / +0.8 at 0.05) and EI_flap does not (-21). Same EA and mass,
less EI points at the axis about which the deck states EI (elastic centre vs reference axis, a parallel-axis
term EA d^2) before the material: candidate H7, to check on the deck's 6x6 off-diagonals at 0.0-0.15.

**H6 result (session 12; jobs modes 11625285, S-1 11625286, es 0.15; yaml `logs/check/l3/t9/h6/`).** Modes,
signed shift vs the unramped run: **1F +0.83, 1E +0.52, 2F +0.98 %** (predicted +0.8 / +0.6 / +1.6): held. Vs
BeamDyn now -1.84 / -3.55 / -2.43 %: H6 is real and secondary. S-1 rows: 0.0594-0.0927 EI_flap -23/-31/-29 ->
-16/-22/-22 and EI_edge +0.7..+1.3 pt (rise: held); EI_edge 0.81-0.83 -21.5/-21.7 -> +1.3/-0.7 (held); mass / EA
0.93-0.95 +43/+34 / +60/+48 -> +7/+4 / +11/+5 (held). **Refuted row: 0.0417 EI_flap -19.2 -> -10.9** (8.3 pt;
also 0.026 +2.3, 0.11 +6.7, 0.13 +3.6) while mass and EA there stay within 0.2 pt: the S-1 EI_flap of a root
station is not a local section value (the response extraction carries the neighbouring cap start), so the
"-21 % at 0.05 with EA and mass agreeing" that motivated H7 is partly the cap-start step. EI_edge 0.11-0.18
unchanged (-26..-20: held, this deficit is real and not H6). Side effects: mid-span EI_edge -1.82 -> +0.22 (the
0.80-0.85 stations), GJ 0.6-0.85 -32.2 -> -35.1, EI_flap 0.91-0.93 +8/-5 -> -10/-26 (free-end extraction zone).

**Audit fix (b9192df, 4 RED tests then 28 passed; rerun `logs/check/l3/t9/audit-es025-zero/`, job 11625287).**
Fails per station 0.3 / 0.6 / 0.9 / 0.97: 42 / 0 / 26 / 7 -> 42 / 0 / 14 / 0. The 0.97 artefact is gone. Left
at 0.9: caps 3.09 vs 3.01 mm (+2.7 %, tol 1 %) on every cap element, and F5 mass +1.2 % / axial +1.8 %: the mesh
samples each layer on the 50-point reference-axis grid (0.87, 0.905) and loses the yaml's own taper knots in
between. Real, small, mesher-side; the 0.3 fails are unchanged and still to read.

**Bredt of the yaml section (0630e71, `tools/l3_bredt_yaml.py`, worker; `logs/check/l3/t9/bredt/`).** Thin-wall
multi-cell, G12 t walls, mid-line offset by t/2 (outer mould line in brackets), vs BeamDyn GJ_sc (K66 condensed
over the shear strains). IEA-3.4: 0.2 -0.1 (+6.4), 0.3 -5.6, 0.4 -9.2, 0.5 -9.0, 0.6 -8.6, 0.7 -6.5, 0.8 -6.6,
0.85 -3.2 %: **the deck GJ is reachable from the yaml geometry** (it lies between the offset and the outer-line
Bredt); our shell FE is -13 / -32 % in the bands, well below. TE cell exterior wall = 61 % (0.5) / 55 % (0.8) of
the flexibility, LE cell 28-30 %, box 1-5 %, each web 4-6 %. IEA-15 control: -1.4 / -4.7 / -12.0 / -19.4 % at
0.3 / 0.5 / 0.7 / 0.8 (outboard a thin wall misses SONATA's solid TE fill), TE share 51 / 68 %.

**Mesh Bredt (session 13, Handoff 3 item 1; `tools/l3_bredt_mesh.py`, selftest: thin circle -0.0003 %, two cells =
`multicell_bredt`).** Same graph Bredt as `tools/bredt_gj.py model_bredt` on the S-1 element row (es 0.15, mid surface,
wall G t = A66), at the yaml-Bredt stations; FE = S-1 GJ (`s1-iea34-es0.15.json`) interpolated at the row z. FE vs
deck there: 0.3 -53, 0.4 -4.7, 0.5 -6.9, 0.6 -15, 0.7 -29, 0.8 -43, 0.85 -41 %.
**Prediction (written before the run).** (1) 3.4 mesh / yaml Bredt within +-5 % at 0.3-0.85 (same mid line; on the 15
the mesher A66 is within 1 % of the yaml stack, T27, and cell areas cost 2-3 pt). Then FE / mesh Bredt is -10..-35 %
at 0.6-0.85 and the outboard deficit is element or extraction, not section. Refuted if mesh / yaml < -5 % at
0.6-0.85 (mesher side: TE closure, offset; T7). (2) At 0.3 the FE -53 % is extraction (mesh Bredt near the yaml).
(3) IEA-15 control (es 0.25): mesh / yaml within +-5 % at 0.3-0.8, FE / mesh within +-10 pt.
**Result (jobs 3.4 11625341, 15 11625338; `logs/check/l3/t9/bredt/mesh-iea*.json`, `bm34c-*.log`).** Section = every
element the plane z crosses (a row picked by centroid z fails on the 3.4: its node rings are not at constant z, the
plane cuts the neighbouring rows, first runs left 0.28 / 0.78-0.80 open; the 15's rows are planar, same numbers both
ways). All stations close as 3 cells, no open wall. % vs deck GJ_sc at the row z:

| eta | 0.2 | 0.25 | 0.3 | 0.35 | 0.4 | 0.5 | 0.6 | 0.7 | 0.75 | 0.8 | 0.82 | 0.85 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3.4 mesh Bredt | -2.8 | -9.8 | -10.9 | -9.0 | -8.2 | -8.2 | -8.2 | -6.7 | -6.7 | -8.2 | -13.1 | -4.0 |
| 3.4 mesh / yaml | -3.1 | | -5.3 | | +0.7 | +0.9 | +0.8 | +0.3 | 0.0 | -2.3 | | -1.2 |
| 3.4 FE / mesh | | -67 | -44 | -14 | +4.8 | +4.7 | -2.7 | -17.7 | -21.0 | -32.7 | -32.0 | -31.1 |
| 15 mesh / yaml | | | +0.2 | | | 0.0 | | +0.4 | | +11.4 | | |
| 15 FE / mesh | | | -1.5 | | | -0.6 | | -4.2 | | +12.7 | | |

TE-cell share of the flexibility, mesh / yaml: 0.54-0.62 / 0.55-0.62 on the 3.4 (equal); 15 at 0.8 0.58 / 0.68.
Verdicts: (1) **held**: our 3.4 section is the yaml section in torsion (within 1 % at 0.4-0.85, -3..-5 % at
0.2-0.3), so the mesher (TE closure, offset) is not the 3.4 outboard GJ deficit; the deficit of the S-1 FE GJ
against its own section grows smoothly from -3 % (0.6) to -18 / -21 / -33 % (0.7 / 0.75 / 0.8) and stays there:
**element or extraction, and 3.4-specific** (the 15's FE tracks its section within 4 pt to 0.7). (2) **held**: at
0.25-0.35 the mesh Bredt is -9..-11 % and the FE -22..-71 %, the inboard extraction swing seen on the 15 (T39).
(3) **held to 0.7, refuted at 0.8**: the 15 mesh is +11 % over the yaml Bredt there, with a stiffer TE cell (share
0.58 vs 0.68); the 15's FE follows its mesh (+12.7 vs deck +1.5), so at 0.8 the IEA-15 S-1 GJ is near the deck by a
mesh-vs-yaml TE difference, not by agreement (the yaml Bredt is -19 % there; SONATA's solid TE fill). Not chased.
The T40 extractor was proven on prismatic tubes only (T39/T40, residual 3-14 % on long or oblique tubes from section
distortion); a short, thin tip section of the 3.4 is outside what it was proven on.
**Next (element vs extraction, no derivative):** the S-7-style integral on the 3.4: FE twist increment over
0.6-0.85 under the same shear-flow tip torque vs the integral of T / GJ of the mesh Bredt. Ratio ~1: the FE is
torsionally right and the S-1 station derivative is wrong there (extraction); ratio ~1.2-1.3: the FE is soft
(element, then a tapered closed-form tube through MITC4 at es 0.15).
**Result (worker, `tools/l3_twist_integral.py`; jobs 3.4 11625348, 15 11625349; `logs/check/l3/t9/twist/`).** Same
T40 shear-flow tip torque and linear rigid-rotation fit (`rigid_rotation_fit`, lstsq, so no large-angle bias) at every
mesh station; integrand T / (GJ t_z) (rate per arc length); mesh Bredt every 0.5 m (all closed); delivered torque
about the tangent 0.991-1.000e6 (3.4) / 0.962-1.000e6 (15) over 0.3-0.9. FE delta-theta / integral:

| band | 3.4 / mesh Bredt | 3.4 / deck | 15 / mesh Bredt | 15 / deck |
|---|---|---|---|---|
| 0.3-0.6 | 1.014 | 1.116 | 1.002 | 1.054 |
| 0.4-0.6 | 0.982 | 1.078 | 1.000 | 1.056 |
| 0.6-0.85 | **1.314** | 1.434 | 0.952 | 1.069 |
| 0.7-0.85 | **1.380** | 1.506 | 0.932 | 1.050 |

3.4 slices (0.05): 0.55 1.025, 0.60 1.043, 0.65 1.176, 0.70 1.237, 0.75 1.281, 0.80 1.530 (0.85-0.95 tip-affected:
1.09, 1.79; 0.25-0.35 inboard 2.67, 1.36). 15 slices 0.98-1.03 at 0.3-0.75, 0.93 / 0.88 at 0.75-0.85.
**Verdict: pre-registered "FE soft" reading holds** (1.31 over 0.6-0.85); the S-1 station derivative is right there
(its -18..-33 % is the same softness), and right at 0.4-0.6 for both blades. Caveat: the twist is a least-squares
rigid fit over all section nodes, so in-plane section distortion (no diaphragms; thin, tapered outboard 3.4 section)
reads as twist: "soft" = element stiffness **or** section distortion, not yet split. The 15 does not show it
(0.93-1.0), so it is specific to the 3.4 outboard geometry (chord ~2 m, TE cell area 0.09-0.14 m2 at 0.7-0.8).

**H7 on the deck side (session 13, Handoff 3 item 2; read of the 3.4 deck 6x6, no run).** The S-1 reference is
already the elastic-centre principal value (`beam_reference.load_beamdyn_blade` `ei_*_ec`: Schur complement of K33,
diagonalised). At 0.05 / 0.10 / 0.15 / 0.20 the elastic centre sits 0.14 / 0.42 / 0.51 / 0.40 m from the reference
axis along the edge direction (K23/EA; flap offset <= 0.04 m), so raw K44 is 1-15 % above the principal edge value
(6.66 vs 5.50e9 at 0.15) and that term is already removed; principal axes turn <= 11 deg (hi/lo 1.2-2.3), worth
<= 2 % on EI_edge. **H7 cannot explain the root EI_edge -20..-26 %: the deficit is ours (section or extraction).**
Next for the root: a yaml-side EI (line integral of E t about the elastic centre, as `t37_section_compare` does on
the mesh) at 0.05-0.20, to split section from extraction the way the Bredt pair does for GJ.

## Handoff 2 (2026-10-26, session 11 close, ~180 k context)

**Resume.** Read section "T9 IEA-3.4-130" only. Engram still down. Logs `logs/check/l3/t9/`, sbatch
`logs/check/l3/t9/job.sbatch` (runs "$@"); 3.4 jobs take 3-5 min (S-1) and ~10 min (modes) on `sequana_cpu`.
**State.** T9a/b done; T9c measured: S-1 sections as predicted (EA chamfer scaling holds on a glass-cap blade, so
T6 applies to both); GJ much worse than predicted; parked modes -2.7 / -4.0 / -3.4 % vs BeamDyn (refuted).
**Next, in order.** 1. Test H6 with a prediction written first: a "deck-equivalent" copy of the 3.4 yaml (under
`logs/check/l3/t9/`, never in tests/) where every partial-span layer ramps linearly to zero over the adjacent deck
interval (0.05-0.10, 0.80-0.85, 0.94-0.95) instead of stepping; re-run `tools/l3_s4_blade.py` es 0.15 and
`tools/l3_s1_blade.py`. H6 predicts 1F/1E/2F move toward 0..+1.5 % and the 0.04-0.18 / 0.81-0.83 EI deficits close.
If the modes stay low: shear lag at the cap start, then the GJ deficit. 2. Fix the input audit's yaml side: a layer
is zero outside its grid (`np.interp` holds the end value today); RED test with a partial-span layer, then re-run the
3.4 audit (the 0.9/0.97 fails should vanish; the 0.3 ones stay to read). 3. GJ: per-band Bredt of the 3.4 section
(`tools/bredt_gj.py` lives on `fix/cap-arcs`) vs K66 before blaming the TE closure; the 3.4 has a 0.8 m ud TE band.
4. Then T6/T7 on both blades (S6), T8. Uncommitted: `tests/IEA34/bd_freevib/` (generated inputs, rebuilt by
`tools/bd_rotating_freevib.py build --beam iea34`; kept untracked as for the IEA-15).

## Handoff 3 (2026-10-26, session 12)

**State.** Handoff 2 items 1-3 done (H6 held and secondary; audit fixed; yaml Bredt). The two blades agree on
the mid-span section signs (EA, EI_flap, EI_edge, mass, GJ outboard-worse), so T6/T7 apply to both; they differ
in the modes, and the beam check (`tools/l3_beam_station_ratio.py`) puts that in the 3.4 root sections, not in
the eigen solve. **Next, in order.** 1. Mesh-side Bredt on our 3.4 section (generalise `tools/bredt_gj.py`
`model_bredt` from `fix/cap-arcs`): mesh Bredt ~ yaml Bredt with FE low -> element or GJ extraction (then a
closed-form thin tube through the FE extraction); mesh Bredt ~ FE -> mesher TE closure / offset = T7. 2. Root
0.04-0.18: EI_edge -20..-26 is real (H6-proof); check the deck 6x6 axis (elastic centre vs reference axis,
off-diagonals at 0.0-0.20) before our section (H7, weakened for EI_flap by the refuted row). 3. Then T6/T7 on
both blades, T8.

## Handoff 4 (2026-10-26, session 13 close, ~165 k context)

**Resume.** Read section "T9 IEA-3.4-130" from "Mesh Bredt (session 13" to this handoff. Engram: `l3/t9/mesh-bredt`.
**State.** Handoff 3 items 1-2 done: the 3.4 mesh is the yaml section in torsion (Bredt within 1 %); H7 refuted on
the deck side (root EI_edge deficit is ours); the 3.4 FE is torsionally soft at 0.65-0.85 by 1.2-1.5x against its
own section with no derivative (twist integral), the 15 is not. Commits 863ff38 (mesh Bredt) and the twist integral.
**Next, in order.** 1. Split element vs section distortion on the 3.4 at 0.65-0.85 (prediction first): (a) the
least-squares fit residual per station (distortion makes it grow outboard) and the twist fitted on the box nodes
only (caps + webs) vs all nodes; (b) a tapered closed-form tube through MITC4 at es 0.15 with the 3.4 outboard
perimeter, wall G t and taper (extend `tools/diagnose_s1_gj_tube.py` / `tests/test_s1_gj_tube.py`); (c) element
quality of the 3.4 outboard rows (aspect ratio, warp of the quads at the TE / LE, es 0.15 vs 0.10). 2. Root
0.04-0.18 EI_edge: yaml-side EI (line integral E t about the elastic centre) vs mesh vs FE. 3. Then T6/T7 on both
blades, T8. The 15 at 0.8: mesh Bredt +11 % over yaml (stiffer TE cell), noted, not chased.

## Store action map (session 13, explore; built from docs/validation/*.yaml and validation_matrix.py logic)

All 34 L3 rows are blade structure only (no BEM/FSI file in L1-L3 `depends_on`; to confirm with
`git log 4809427..HEAD -- <depends_on>`). Row files were re-measured at d5d0fd2 but `contract.yaml` verdict notes cite
older values, so **14 of 15 stale rows changed value: new verdict, never restamp** (`validation-l3-b3-audit.md:23`);
only group 20 (identity) is a pure restamp. The full suite (`check.sh -q`) prints no residuals: values come from
`regression --group G --scope S` (groups 5, 20, 23, 32, 34) and `pytest -s` / `tools/run_*` (38-42).
Next action per non-trusted row (K3: > 5 % by >= 1 pp rejected; within 1 pp of 5 % without es 0.125 undecided):
- S-1 (38): EI_edge 8.01 rejected; EI_flap 4.53 in band -> es 0.125 then verdict; mass 3.54 trusted; GJ 10.18 rejected; EA 6.68 rejected.
- S-2 (39): LC1 5.26, LC4 5.48 -> es 0.125 then verdict; LC2 (2.54 row) and LC3 (+17 % edgewise) undecided: USER (metric: tip node vs ring spread 9-11 %).
- S-4 (40): parked 2.33/0.56/0.20, rated 1.85/0.44/0.06 -> new verdicts trusted; then gap `rotating_modal_shift` -> bounded with both row ids.
- V-03 (41) 1.14 trusted; V-05 (42) +0.23 trusted (sign flipped); S-7 (34) 7.25 rejected.
- V-02 (32) 7.09 / 8.84 vs TP-75698 (v1.0 blade): USER (v1.1 source or drop).
- Anchor beam (20): pure restamp after `regression --group 20`.
- Group 23: mass vs Escalera = UTD conversion (disowned): USER keep/drop; CCX modes 0-4 (0.01-0.12 %) and static (0.37 %) audit then second arbiter (BeamDyn p20 modes, `tools/bd_static_reference.py`); Escalera DLC 1.4 (xfail 13.2 %): USER drop/pin; Bernardi modes undecided: USER drop / re-point to p20.
- Group 5 mesh convergence: fails its own 5 % at the finest mesh (high modes): finding, USER before audit.
- Absent `sx_mesh_convergence`, `sx_campbell` (ElastoDyn targets): import or out of scope; Campbell re-point to p20 per rpm.
Single/no arbiter (28): S-1 x5, S-7, V-02 x2, V-03, group 5, group 23 x16, 2 absent. CalculiX can be a second arbiter
only for S-7 (`tools/ccx_blade_twist_arbitration.py` at es 0.25, not its 1.0 default; `test_blade_ccx_parity.py`
torsion); S-1 are section integrals (no FE). USER: does the official BeamDyn deck count as `reference_model`
(it does in 39, 40, 42 but is `published numerical` in 34, 38)? That alone settles S-1/S-7 arbitration.
`test_blade_ccx_parity.py` (6 static cases vs CCX) has no row nor group. L3 cannot print CLOSED with rejected rows:
it closes by S5 (every open entry listed with cause and torque sensitivity).
Missing capability rows (tolerances to fix before measuring): MAC per mode parked/rated vs BeamDyn p20 eigenvectors
(to generate) and CCX (`postprocess_ccx_modal.py:201`); 1T and 2E frequencies; centrifugal shift delta-f; Campbell 3-11
rpm vs p20; blade-only transient vs BeamDyn dynamic driver; damping per mode vs deck mu; ply stress if used; rated
twist vs Zhou (group 27, on no level): USER.
Procedure: verdicts are hand entries in `contract.yaml` (target, class, verdict, audited_rev, files, arbiters, note);
values by `regression --group G [--scope S] --write`; new rows: test, `extract --group G --write`, `regression --write`,
plus `groups.yaml`, `residual-patterns.json`, `levels.yaml`; gate `contract --level L3 --rows`.

## Store pass 1 result (2026-10-27, session 14, logs/check/l3/store-pass1/)

Measured at d50d64a / f1c9e91 (f1c9e91 touches only FSI/BEM, no L3 `depends_on`). Verdicts in `contract.yaml`:
- S-1 (es 0.25 / 0.125, mean |rel|): EI_edge 8.01 / 8.34 rejected; EI_flap 4.53 / 4.86 **trusted** (0.14 pp
  margin, refinement moves toward the ceiling); mass 3.54 / 3.64 trusted; GJ 10.74 / 10.70 rejected (signed -10.7,
  shell soft everywhere; was 10.18 before d81931b/a4b3ca2); EA 6.68 / 6.83 rejected. Arbiter -> reference_model.
- S-2 vs BeamDyn p60: LC1 -5.26 / -5.28 rejected; LC4 -5.48 / -5.90 rejected; LC2 node -2.54 / -3.56 and LC3 node
  +10.26 / +11.44 (ring mean +9.01 / +9.78) undecided until the test asserts the ring mean (pass 2). **Finding:** refinement
  moves every S-2 case away from BeamDyn, so the residual is not discretisation.
- S-4 parked +2.33 / +0.56 / +0.20 %, rated +1.85 / +0.44 / +0.06 % (1F/1E/2F) trusted. Centrifugal shift parked -> rated
  7-13 % below BeamDyn (1F -8.1 %, half of it the MBC3 HubRad 3.97 m): gap `rotating_modal_shift` bounded, delta-f row is pass 2.
- S-7 1.07171 (7.17 %) rejected, arbiter reference_model. V-03 -1.14 % trusted (xlsx, still single arbiter). V-05 +0.23 % trusted.
- Group 20 restamped (identity, margins same). V-02: group 32 removed from L3 (`levels.yaml`). Group 23: the two Escalera/UTD
  tests declared `non_validation_tests` and their rows removed by `extract --group 23 --write` (14 rows left).
- Not touched (pass 2 / user): group 23 CCX audit + second arbiter, Bernardi -> p20, group 5 (prints no residual; finding),
  the two absent sx files, S-2 ring-mean metric, missing capability rows (MAC, 1T/2E, delta-f, Campbell).

## Handoff 5 (2026-10-27, session 14 close)

**Resume.** Read "Store pass 1 result" and this handoff; Engram `l3/store-pass1`. Branch l3/core-torsion merged to main
(this session's PR); the full suite never completed on it (SLURM 11625408 / 11625463 TIMEOUT at 4 h, 73 % / 71 %, no
failures beyond the three known reds up to there): a full suite on main with `--time` 8 h is pending and is the first
check of the next session (crates/ changed: rotor_fsi.rs, rotor_physics.rs).
**Next, pass 2, in order.** 1. S-2 ring-mean metric (L12 (4)): the test asserts the tip-ring mean for LC1-LC4, then
LC2/LC3 verdicts (LC3 ring +9.78 % at es 0.125 -> rejected). 2. Bernardi rows -> BeamDyn p20 modes (L12 (3)). 3. Group
23 CCX rows: audit, second arbiter p20. 4. Group 5 prints no residual: make it print, then audit (finding: fails its own
5 % at the finest mesh). 5. Absent sx_mesh_convergence / sx_campbell: import or declare out of scope (USER). 6. Capability
rows: delta-f centrifugal shift (ours 7-13 % below BeamDyn), MAC per mode vs p20, 1T/2E, Campbell. Then T5 (declare reds
with torque sensitivity), T6/T7 model work per core-torsion.md.

## T10a prediction (written before the run, 2026-10-27, session 15, at 4b29bc2)

Closed-form check of the S-2 metric on prismatic thin-wall tubes (`tools/s2_ring_metric_tube.py`, the S-2 load
introduction of `tools/run_s2_static_cases.py`: tip force split equally over the tip-ring nodes, root ring clamped, es
0.25, L 117 m, isotropic wall). Reference: Timoshenko cantilever `P L^3 / (3 EI) + P L / (G A_s)` with the thin-wall
mid-line EI (exact segment integral) and A_s the walls parallel to the load (rect) or A/2 (circle). Cases: rect 4 x 2 m,
t 0.02, load along the 2 m (flap-like) and along the 4 m (edge-like); circle of the same perimeter, one direction.
Prediction: tip-ring mean within 1 % of the closed form in every case; the first-ring node (`argmax z`, the metric the
test used) within 1 % of the ring mean, and the ring max-min spread of the load-direction displacement < 2 %.
Refuted if the ring mean misses 1 %: then the ring mean is not the beam-equivalent translation and the L12(4) metric
needs a different estimator (finding, back to the user). A node-vs-ring spread > 1 % on a prismatic tube means local
load introduction alone produces the blade's 9-11 % node-vs-ring spread order, not the blade section.

## T10a RESULT (2026-10-27, login, `tools/s2_ring_metric_tube.py`, logs `logs/check/l3/t10a/L20.log`, `L117.log`, `L117-sv.log`)

| tube (L 117 m, 1 MN) | tip ring mean vs closed form | node `argmax z` vs ring | ring 1 depth inboard | Saint-Venant tip |
|---|---|---|---|---|
| rect flap (load normal to the 4 m walls) | **+1.241 %** | -1.654 % | -0.181 % | -0.424 % |
| rect edge (load normal to the 2 m walls) | +0.080 % | -0.390 % | -0.312 % | -0.308 % |
| circle | -0.274 % | -0.005 % | -0.279 % | -0.279 % |

L 20 m (0.1 MN): tip ring mean +322 % (flap), +74 % (edge), +0.40 % (circle); 2 depths inboard +2.8 / -0.01 / -0.25 %.
**Prediction refuted for the flat-walled tube.** A tip force split equally over the ring pushes flat walls out of
their plane (no end diaphragm); the tip ring mean carries that plate bending (~0.55 m on the 4 m flange, independent
of L), which decays within about one section depth (Saint-Venant). The circle (no flat wall normal to the load) holds
the prediction. The `argmax z` node is a corner of the rect tube: close to the beam by luck, not by construction.
Estimator `saint_venant_tip` (`tools/run_s2_static_cases.py`): ring means at 1, 1.5, 2, 2.5, 3 x the tip-ring extent
inboard, cubic in z, evaluated at the tip: within 0.43 % of the closed form in all three; its residual equals the
inboard rings' uniform -0.2..-0.4 % (element / thin-wall), not the estimator. Consequence for L12(4): the tip-ring mean
is the beam-equivalent for distributed loads (LC1, LC4) but not for the tip loads (LC2, LC3) on a flat-walled section;
T10b measures the size of that term on both blades before the metric goes back to the user.

## T10b prediction (written before any blade number, 2026-10-27, session 15; jobs 11625682-85 + BeamDyn 3.4 chained)

`tools/s2_ring_metric_blade.py` on IEA-15 and IEA-3.4, es 0.25 / 0.125: node, tip-ring mean, `sv` (rings 1-3 tip chords
inboard) and `sv_far` (3-6 chords); IEA-3.4 reference: BeamDyn v5 p60 static, scales 1e-4 / 2e-4
(`tools/bd_static_reference.py --blade iea34`). Tip chord 0.50 m (15, tip ring 16 nodes), 0.20 m (3.4, 12 nodes).
P1 (estimator stable): `sv_far / sv - 1` within 1 % in every case, blade and mesh; else the decay is longer than the
tip chord and `sv` is not admissible on that blade (finding).
P2 (distributed loads): LC1 and LC4 tip-ring mean within 0.5 % of `sv` (no load introduction at the tip).
P3 (tip loads): LC2 and LC3 tip-ring mean above `sv` (same sign as the tube), by the local term; node vs `sv` moves with
the node's position on the ring, so the known node-vs-ring 9-11 % (LC2) is mostly that local term.
P4 (both blades): the 3.4 shows the same signs as the 15 for P2-P3. No verdict vs BeamDyn is predicted here; the S-2
verdicts are written after the metric decision (L12(4) revisited with these numbers).
BeamDyn linearity gate (as K9): the two scales agree within 0.1 % per case, else that case has no 3.4 reference.

## T10b RESULT (2026-10-27; SLURM 11625682-85 shell, 11625690 BeamDyn 3.4; `logs/check/l3/t10b/*.json|log`,
`logs/check/bd_static/iea34/*/result.json`; tube refinement `logs/check/l3/t10a/L20-es*.log`)

Signed error vs BeamDyn p60, node / tip-ring mean / sv (`saint_venant_tip`):

| blade es | LC1 gravity | LC2 tip flap | LC3 tip edge | LC4 distributed |
|---|---|---|---|---|
| 15, 0.25 | -5.26 / -5.04 / -5.33 | -2.54 / +1.83 / +0.41 | +10.26 / +9.01 / +8.49 | -5.48 / -5.06 / -5.76 |
| 15, 0.125 | -5.28 / -5.07 / -5.18 | -3.56 / +3.09 / +0.94 | +11.44 / +9.78 / +9.02 | -5.90 / -5.46 / -5.74 |
| 3.4, 0.25 | +1.96 / +1.96 / +1.96 | +18.25 / +20.78 / +20.64 | +5.94 / +4.91 / +4.78 | -0.02 / +0.23 / +0.23 |
| 3.4, 0.125 | +2.04 / +2.05 / +2.05 | +15.37 / +19.18 / +17.59 | +6.92 / +6.30 / +5.97 | -0.50 / -0.27 / -0.41 |

BeamDyn 3.4 linearity gate: scales 1e-4 / 2e-4 agree within 0.04 % (LC2), < 0.001 % otherwise: PASS.
- P1 **refuted**: `sv_far / sv` reaches -2.5 % (15 LC4 0.125), -1.9 % (15 LC3), +1.6 % (3.4 LC2 0.125). The arithmetic
  ring mean steps where the ring node count changes (20 -> 26 nodes on the 15; same mechanism as core-torsion C6), and
  the cubic extrapolation amplifies the steps. `sv` is not admissible on the blades as built.
- P2 partly held: LC1 / LC4 ring vs sv within 0.3 % except 15 LC4 es 0.25 (+0.74 %).
- P3 held on both blades: LC2 / LC3 ring above sv in all 8 (+0.1..+2.1 %), and the node-vs-ring spread grows with
  refinement (15 LC2 node -2.5 -> -3.6 %, ring +1.8 -> +3.1 %). P4 held (same signs on the 3.4).
- Theory (tube L 20, es 0.25 -> 0.125): circle tip-ring mean +0.40 -> +1.00 % vs closed form, node vs ring -0.89 ->
  -1.26 %, ring 2 depths inboard -0.25 -> -0.04 %. Equal nodal forces are point loads on a Reissner-Mindlin shell
  (local deflection log-singular): the tip-ring term grows under refinement while the interior converges. This is a
  candidate cause of the pass-1 finding "refinement moves every S-2 case away from BeamDyn" for the tip-load cases;
  LC1 (consistent body force) moves 0.02 pp only.
- **Verdicts do not depend on the metric** (15, es 0.125, K3): LC1 rejected (5.07-5.28), LC2 in band (0.94-3.56),
  LC3 rejected (9.0-11.4), LC4 rejected (5.46-5.90) under node, ring and sv. The metric moves LC2 / LC3 by up to 4 pp.
- **Second-blade finding**: IEA-3.4 LC2 (tip flap) +15..+21 % softer than its BeamDyn, while its LC1 / LC4 flap agree
  within 2 %; the 15 shows the opposite pattern (distributed flap stiff 5-6 %, tip flap in band). Tip loads weight the
  outboard stiffness, distributed loads the inboard: the 3.4 outboard flap is soft (cf. T9c modes refuted). Not chased.
- Slice note: commit 9814a37 is 322 authored lines; the other 1852 are generated BeamDyn inputs (precedent: the 15's
  s1e-4 decks are tracked, other scales are not).

## T10c prediction (written before the run, 2026-10-27, session 15)

Load introduction (`tools/run_s2_static_cases.py`): tip load `distribution: arc` = uniform line traction on the free tip
edge (nodal forces by arc-length weight, `ring_arc_weights`); LC4 `distribution: area` = uniform surface traction per
slice (tributary-area nodal forces), included because its equal per-node forces are the same point-load defect.
Metric: arc-weighted tip-ring mean. Tube (L 20, es 0.25 -> 0.125): circle arc mean within 0.3 % of the closed form at
both sizes and moving < 0.2 pp between them (converges; equal-split moved +0.6 pp); rect flap keeps a local flange
term (physical for a line load on a free flat wall) but moves < 0.5 % of itself between sizes. Blades (T10c run):
arc tip mean moves < 0.5 pp between es 0.25 and 0.125 in LC2 / LC3 (equal-split ring moved 1.3 pp on the 15 LC2).
Refuted if the circle still drifts > 0.2 pp: then the line traction is not the cure and the metric goes back to the user.

## T10c RESULT, tube (2026-10-27, login; logs `logs/check/l3/t10a/L20-arc-es*.log`, `L117-sv.log`, `L117-es0125.log`)

**Prediction refuted, and the T10b mechanism with it.** On a uniformly spaced ring the arc-weighted line traction gives
exactly the equal nodal split (identical numbers to 6 digits), so the circle's +0.40 -> +1.00 % at L 20 is not a
point-load singularity. Reading: any edge load not tangent to the wall has an out-of-plane component on the free edge
(uniform y-traction on a circle: normal part ~ sin theta); its transverse-shear edge layer is of order t (0.02 m) and
es 0.25 / 0.125 do not resolve it, so the tip-ring value is mesh dependent under any non-tangent edge load.
L 117, es 0.25 -> 0.125, vs closed form:

| tube | tip ring | ring 1 depth inboard | ring 2 depths inboard |
|---|---|---|---|
| rect flap | +1.24 -> **+1.54** | -0.18 -> +0.10 | -0.43 -> -0.15 |
| rect edge | +0.08 -> **+0.29** | -0.31 -> -0.11 | -0.31 -> -0.12 |
| circle | -0.27 -> -0.09 | -0.28 -> -0.10 | -0.28 -> -0.10 |

Interior rings at >= 2 section depths converge to the closed form (error / 3 per halving); the tip ring on flat walls
grows. Without MPC (no rigid tip coupling) and without a shear-flow (wall-tangent) tip traction, the converging
comparison is at an interior station against BeamDyn's own displacement at that station. LC4 `area` weighting and
`ring_arc_weights` stay (they remove the node-count steps of T10b on non-uniform blade rings); not yet run on blades.

## Log

- L1 (2026-10-26, session 10) User: "vamos con L3"
- L2 (2026-10-26, session 10) User chose option 1 "Sensibilidad primero (recomendado)" over "Arreglar el modelo primero" and "Solo declarar".
- L3 (2026-10-26, session 10) User: "revisa si se esta ejecutando bien para no perder el tiempo". Parent: jobs healthy; flap sign corrected (ours stiffer), signed sensitivity re-run.
- L4 (2026-10-26, session 10) User: S6 verbatim. Parent: the "10 % softer" is S-2 LC3 static tip-edge (point load, undecided row) and S-1 EI_edge 8 % (mid-span extraction), both vs the BeamDyn/SONATA 6x6; the edgewise dynamics agree (S-4 1E +0.56 / +0.44 % vs BeamDyn, a 10 % EI error would move 1E ~5 %). Model reds with diagnosed causes: flap/EA (cap chamfer, T6), torsion (TE solid closure, T7). Mode shapes are not in the store: T8.
- L5 (2026-10-26, session 10) User: "en algun momento te sugeri utilizar otra pala del repositorio oficial de NREL y comparar tambien para ver si son defectos de modelo o de implementacion del core de solucion". Parent: suggested 2026-10-23 (blade-input-audit L9), inventoried as T17, parked behind T18 (windIO v2 importer), never run. A SONATA-chain blade (3.4) tells systematic-vs-turbine-specific and doubles the evidence for T6/T7; core is bounded by L1/L2 closed, S-4 within 2.4 %, and the section-only PDF check; SONATA convention vs physics needs a BECAS blade (3 webs, not supported).
- L6 (2026-10-26, session 10) User chose "1. IEA-3.4 antes de T6 (recomendado)" over "3.4 y además una pala BECAS" and "T6/T7 primero".
- L7 (2026-10-26, session 11) User: "seguimos con L3, T9". Parent: deck/yaml pairing checked consistent (section T9), T9a pre-registered before any 3.4 run; T9b reader fixes; T9c sections held, modes refuted, H6 recorded (Handoff 2).
- L8 (2026-10-26, session 12) User: Handoff 2 items 1-4 in order, verbatim: "1. Probar H6 con la predicción escrita antes: una copia del yaml con las capas en rampa en vez de escalón, y volver a correr modos y S-1. 2. Corregir la auditoría para que una capa valga cero fuera de su grilla, test primero. 3. Bredt de la 3.4 contra el GJ del deck antes de culpar al cierre del borde de fuga. 4. Recién después, T6 y T7 sobre las dos palas." Then asked whether the two blades behave consistently and whether the inconsistencies are model, solver or element; parent: model conventions and reference sampling, solver core unlikely, torsion still open between mesher, extraction and element.
- L9 (2026-10-26, session 13) User: "continuamos". Parent: Handoff 3 items 1-2 (mesh Bredt, H7 deck check), then the twist integral (worker); Handoff 4.
- L10 (2026-10-26, session 13) User, verbatim: "recuerda que no es que yo quiera arreglar el modelo de una u otra pala, yo quero que aeroelast sea cientificamente confiable para el proposito que te di". Parent: the blades are instruments to find code defects; the 3.4 outboard torsional softness (1.3x its own section) is a core-reliability question (element or section-distortion measure), which is why Handoff 4 item 1 comes before any blade-model work. Also asked for one PR of l3/twist-p1 to main.
- L11 (2026-10-26, session 13) User: "si, eso que me estas proponiendo es lo que yo crei que estavamos haciendo desde que desarrolle el store de pruebas de validacion" and "ok puedes empezar con eso en paralelo?". Parent: L3 is closed from the store contract (restamp stale, audit unaudited, decide undecided, second arbiter CalculiX where single, missing capability rows such as mode shapes); store snapshot 2026-10-26: L1 closed (157 trusted), L2 closed (79), L3 open (undecided 10, stale 15, unaudited 9, absent 2, 28 single/no arbiter, 1 blocking gap), L4 open (20 stale), L5 open and blocked (unaudited 11, absent 8; one row touches coupling code). Core-torsion findings in `odd/tasks/core-torsion.md`. Read-only action map launched (explore).
- L12 (2026-10-26, session 13) User decisions on the store action map, verbatim choices: (1) BeamDyn deck class: "reference_model en todas (recomendado)"; (2) V-02: "Darla de baja (recomendado)"; (3) group 23: "Bajarlas; Bernardi re-apuntado a p20 (recomendado)"; (4) S-2 metric: "Media del anillo de punta (recomendado)". Parent: pass 1 (measure at HEAD without --write, es 0.125 runs in the 1 pp band, new K3 verdicts, deck reclassified, V-02 and UTD/Escalera rows dropped, group 20 restamp, gap rotating_modal_shift bounded) by a worker; pass 2 (S-2 ring-mean metric, Bernardi -> p20, missing capability rows) after.
- L13 (2026-10-27, session 14) User: "ya terminaron las tareas del hpc vamos a continuar cerrando l3". Parent: read the 12 pass-1 jobs, ran the S-1 metric dump at es 0.25 / 0.125 and `regression --write` for groups 20 and 34, wrote the pass-1 verdicts, bounded `rotating_modal_shift`, applied the L12 drops (section "Store pass 1 result").
- L14 (2026-10-27, session 14) User: "si mezcla todo y vamos a otra sesion"; chose "Mergear ya (recomendado: dejo la falta escrita)" over running the full suite first. Parent: merged l3/core-torsion with the incomplete full suite recorded (Handoff 5).
- L15 (2026-10-27, session 15) User: "dame un orden de prioridad a resolver de los reds de L3". Parent: order 0 full suite, 1 S-2 ring mean, 2 delta-f centrifugal (HubRad), 3 twist definition then S-7/GJ re-measure, 4 modes (group 5, MAC, 1T/2E, Campbell), 5 T6 then T7, 6 EI_edge (extractor vs W), 7 store hygiene, 8 T5.
- L16 (2026-10-27, session 15) User: "ok vamos a correr la suite completa, para dejar el estado actual, pero como son 8h parado podemos hacer otras cosas en paralelo, dime cuales?". Parent: full suite already running (SLURM 11625665, 8 h, 3185da6 = origin/main tree); parallel work only in a separate worktree, no maturin rebuild (shared venv `_aeroelast`).
- L17 (2026-10-27, session 15) User: "si" + S10; decisions S8, S9. Parent: worktree `../wt/l3-pass2` on branch l3/pass2-s2-ring from origin/main 4b29bc2; T10, T11 added. Finding: the S-2 "shell" value is `argmax(z)`, the first tip-ring node in mesh order (`tools/run_s2_static_cases.py:119`), an arbitrary material point, while the tip load is split equally over the ring.
- L18 (2026-10-27, session 15) Parent: T10b numbers and the metric question; user chose S11. The solver has no MPC / rigid tip coupling, so BeamDyn's axis point load cannot be reproduced exactly; the line traction is the non-singular substitute.
