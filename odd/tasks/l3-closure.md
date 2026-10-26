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

- T1 | S2 S4 | parent, SLURM | measured, restamp pending | Re-measure stale groups at HEAD (`regression --group G --scope <row nodes>`, no `--write`): 38 S-1, 34 S-7, 39 S-2, 40 S-4, 41 V-03, 42 V-05, 32 V-02, 20 anchor beam; restamp unchanged rows, report moved ones.
- T2 | S2 S4 | parent | pending | Audit group 23 (CCX parity, Bernardi, Escalera) and group 5 (mesh convergence) rows; decide S-2 LC2/LC3 undecided.
- T3 | S5 | parent | pending | Gap `rotating_modal_shift`: its cause was "no row"; group 40 rows exist, link them and set the status.
- T4 | S2 S3 | parent | done: twist +0.34 %, flap -0.34 %, edge 0.005 % torque | Sensitivity tool: rigid BEM at rated with the twist distribution perturbed by each red's effect; torque/power deltas.
- T9 | S7 S4 | parent/worker | in progress: (a) d7d000e, (b) 188c997, (c) S-1 + modes measured 94237f7; modes refuted, diagnosis open (Handoff 2) | Second blade IEA-3.4-130 before T6 (T17 of blade-input-audit; v1 yaml at IEAWindSystems commit 824f8d6, 2 webs, glass caps, SONATA 6x6 = same chain as the 15; ANBA mass structural only, issue #5): (a) pre-register the expected signed biases (EA, EI_flap, EI_edge, GJ outboard, 1F/1E/2F) from the IEA-15 diagnosis; (b) ingest + input/mesh audit (clamps the 15 never binds); (c) S-1 section vs its SONATA 6x6 and S-4 modes vs its OpenFAST deck; (d) T6/T7 must then improve both blades. BECAS blades (10, 22: 3 webs) out of this task.
- T6 | S6 S4 | parent/worker | pending | Cap chamfer (SONATA convention, blade-input-audit Handoff 8 item 1): write the prediction (EA, EI_flap, S-2 LC1/LC4, 1F) first, section check vs the SONATA PDF, then the FE anchors. Target: flap -5.5 % and EA +6-9 % reds.
- T7 | S6 S4 | parent/worker | pending | Faithful T37 TE closure (bond-line keypoint + solid closure; Handoff 9 item 3, predicted +5..7 pt GJ at 0.7-0.9): prediction, section check, then S-7 / S-1 GJ / torsional mode. One rule at a time, after T6.
- T8 | S6 | parent | pending | Mode-shape rows: MAC per mode (1F 1E 2F 2E 1T, parked and rated) of our modes vs BeamDyn p20 mode shapes; frequencies alone do not pin the shapes the rotating torque analysis needs.
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

## Log

- L1 (2026-10-26, session 10) User: "vamos con L3"
- L2 (2026-10-26, session 10) User chose option 1 "Sensibilidad primero (recomendado)" over "Arreglar el modelo primero" and "Solo declarar".
- L3 (2026-10-26, session 10) User: "revisa si se esta ejecutando bien para no perder el tiempo". Parent: jobs healthy; flap sign corrected (ours stiffer), signed sensitivity re-run.
- L4 (2026-10-26, session 10) User: S6 verbatim. Parent: the "10 % softer" is S-2 LC3 static tip-edge (point load, undecided row) and S-1 EI_edge 8 % (mid-span extraction), both vs the BeamDyn/SONATA 6x6; the edgewise dynamics agree (S-4 1E +0.56 / +0.44 % vs BeamDyn, a 10 % EI error would move 1E ~5 %). Model reds with diagnosed causes: flap/EA (cap chamfer, T6), torsion (TE solid closure, T7). Mode shapes are not in the store: T8.
- L5 (2026-10-26, session 10) User: "en algun momento te sugeri utilizar otra pala del repositorio oficial de NREL y comparar tambien para ver si son defectos de modelo o de implementacion del core de solucion". Parent: suggested 2026-10-23 (blade-input-audit L9), inventoried as T17, parked behind T18 (windIO v2 importer), never run. A SONATA-chain blade (3.4) tells systematic-vs-turbine-specific and doubles the evidence for T6/T7; core is bounded by L1/L2 closed, S-4 within 2.4 %, and the section-only PDF check; SONATA convention vs physics needs a BECAS blade (3 webs, not supported).
- L6 (2026-10-26, session 10) User chose "1. IEA-3.4 antes de T6 (recomendado)" over "3.4 y además una pala BECAS" and "T6/T7 primero".
- L7 (2026-10-26, session 11) User: "seguimos con L3, T9". Parent: deck/yaml pairing checked consistent (section T9), T9a pre-registered before any 3.4 run; T9b reader fixes; T9c sections held, modes refuted, H6 recorded (Handoff 2).
