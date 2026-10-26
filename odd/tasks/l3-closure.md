# L3 closure: blade structure, sensitivity first

Step 2 of `odd/tasks/thesis-validation-plan.md` (L3, IEA 15 MW official blade). Store `docs/validation/`, tool
`tools/validation_matrix.py`. Model history and diagnosis of the reds: `odd/tasks/blade-input-audit.md` Handoff 10.

## Specs

- S1 (user, verbatim): "vamos con L3"
- S2 (user decision, option chosen verbatim): "1. Sensibilidad primero (recomendado)" — "Re-medir las stale y auditar las pendientes. Para cada red, medir su efecto en torque/potencia con BEM rígido y twist perturbado (barato, no depende de L5). Si el efecto es chico, se declara como limitación. Si es grande, recién ahí se toca el modelo (cierre TE fiel, T37)."
- S3 (thesis plan step 2 acceptance): 5 %; every remaining red gets its torque/power sensitivity measured and is declared as a limitation (EA: no torque path, declared).
- S4 (repository physics discipline): no pass band widened; a moved anchor is a finding; `regression --write` never over a moved residual without reporting it.
- S6 (user L4, verbatim): "si hay opcion de mejorar el modelo hay que mejorarlo. Las formas modales son muy importantes porque luego voy a necesitar hacer analisis de frecuencias en rotacion para ver como esas viraciones modifican o no la se~nal de torque. una pala 10% mas blanda (no se contra que referencia) hace que sus oscilaciones sean diferente 10% diferentes a la respuesta fisica"
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
**Next, in order.** 0. Session 11: T6 (prediction before the run), then T7, then T8; items 1-3 below can run alongside. 1. Restamp the stale rows whose values did not move (S-1 mass/EI_flap, S-4 x2, V-03, V-05,
anchor beam 20, V-02 32 if `same`), record moved ones as findings (S-7 7.25 -> 7.17 %, T40). 2. T3: group 40 rows
exist -> `rotating_modal_shift` gap to `bounded` with the rows listed. 3. T2: group 23 unmapped rows and group 5
need a printed residual or a manual verdict; S-2 LC2 (-2.5 %, in) and LC3 (+10.26 %, the edgewise ratio test
passes with a different band: read it before deciding). 4. T5: declarations with the sensitivity numbers above
(EA: no torque path), then `contract --level L3`.

## Log

- L1 (2026-10-26, session 10) User: "vamos con L3"
- L2 (2026-10-26, session 10) User chose option 1 "Sensibilidad primero (recomendado)" over "Arreglar el modelo primero" and "Solo declarar".
- L3 (2026-10-26, session 10) User: "revisa si se esta ejecutando bien para no perder el tiempo". Parent: jobs healthy; flap sign corrected (ours stiffer), signed sensitivity re-run.
- L4 (2026-10-26, session 10) User: S6 verbatim. Parent: the "10 % softer" is S-2 LC3 static tip-edge (point load, undecided row) and S-1 EI_edge 8 % (mid-span extraction), both vs the BeamDyn/SONATA 6x6; the edgewise dynamics agree (S-4 1E +0.56 / +0.44 % vs BeamDyn, a 10 % EI error would move 1E ~5 %). Model reds with diagnosed causes: flap/EA (cap chamfer, T6), torsion (TE solid closure, T7). Mode shapes are not in the store: T8.
