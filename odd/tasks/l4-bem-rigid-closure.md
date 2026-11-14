# L4: built-in BEM, rigid blade, vs OpenFAST AeroDyn (thesis plan step 3)
Step 3 of `odd/tasks/thesis-validation-plan.md`, opened after L3 closed (plan L4). Builds on the B4 audit
`odd/tasks/validation-l4-b4-audit.md` (verdicts at f1a1a18, 2026-10-20; L4 left open). Runner: `scripts/aeroenv.sh`, `PYTHONPATH=$PWD`; SLURM `--ntasks` only.

## Specs
- S1 (user, verbatim): "cerrar los pendientes antes de L3 Luego Cerrar L3, y comienza con L4"
- S2 (thesis plan, step 3): built-in BEM, rigid blade, against OpenFAST AeroDyn with the same corrections, tower and wake models off;
  Zhou rigid means (16.11 MW, 2.53 MN); acceptance 5 % on power, thrust, torque, spanwise loads.
- S3 (thesis plan S3): every compared number has an independent reference and a written tolerance (5 % relative ceiling); a red that
  stays red is reported with its cause and its torque/power sensitivity, never hidden by widening.
- S4 (user, B4 audit S5, verbatim): "para la turbina IEA 15mw aunque no tenemos arbitro analitico o pulicaciones tenemos el repositorio
  oficial de OpenFAST por lo que ese arbitro es valido" — the official OpenFAST model on the official deck settles alone (`reference_model`).

## Tasks
- T1 | S2, S3 | parent | done (L2; baseline write awaits user) | Entry measurement at HEAD: `tools/validation_matrix.py contract --level L4 --rows`; list what changed in
  `solvers/bem/` and force projection since f1a1a18; decide which B4 verdicts need re-measurement (group regression, not the suite).
- T2 | S2, S3 | worker+SLURM | pending | Spanwise loads (B4 follow-up 1, F7): criterion written first, then AeroDyn dt 0.25 s on A1; decide
  whether dt = 1 s is a valid spanwise reference; verdict for spanwise .0-.2.
- T3 | S2, S4 | worker | pending | V-01 (B4 follow-up 2): re-point at OpenFAST on the official deck at rated, or record it unarbitrated.
- T4 | S2 | worker | pending | Zhou et al. 2025 rigid means (16.11 MW, 2.53 MN): matched operating point, our BEM vs Zhou vs AeroDyn.
- T5 | S2, S3 | worker | pending | Matched-options audit: our BEM corrections (tip/hub loss, skew, induction, polars, Viterna) vs the AeroDyn
  deck used as arbiter; tower and wake off on both sides; record any mismatch and its torque/thrust effect.
- T6 | S3 | parent | pending | Residual patterns for groups 59/60 and the spanwise two-quantity print (B4 follow-up 3); L4 closure or declarations.

## Log
- L1 (2026-11-14, user, verbatim): "cerrar los pendientes antes de L3 Luego Cerrar L3, y comienza con L4". Parent: L3 closed
  (thesis plan L4); this document opened from the B4 audit follow-ups 1-3 and the thesis plan step 3 acceptance.
- L2 (2026-11-14, parent, T1 RESULT; logs `logs/check/l4/contract-l4-rows-entry.log`, `regress-g{11,59,60}.log`, SLURM 11738995,
  OpenFAST v5.0.0 `$SCRATCH/conda-envs/openfast/bin/openfast`). `contract --level L4 --rows`: L4 open, all 24 rows stale by
  5acbbd8 5cccc88 2d6f448 5c15e8d 75ba976 (only 75ba976 touches `solvers/bem/`: skip an unrealisable ring). `regression --group 11`
  (no --write): 15 changed, 3 same, 3 unmapped (spanwise, known), 2 declared; source digest unchanged. Every move is <= 0.011 pp
  (largest: yaw/shear .0 torque 0.3567 -> 0.3674 %; A1 .0 torque 0.6646 -> 0.6716 %; A3 .2 torque 2.1526 -> 2.1531 %); thrust moves
  -0.002 pp, torque +0.005..+0.011 pp; all rows stay far below 5 %. Groups 59 and 60: `regression` refuses, no residual patterns
  (B4 follow-up 3 = T6). Not written: the rule forbids `--write` over moved residuals; the moved baseline is reported for the user.
