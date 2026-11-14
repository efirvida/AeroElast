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
- T1 | S2, S3 | parent | pending | Entry measurement at HEAD: `tools/validation_matrix.py contract --level L4 --rows`; list what changed in
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
