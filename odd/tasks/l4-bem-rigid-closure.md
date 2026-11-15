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
- T1 | S2, S3 | parent | done (L2-L3) | Entry measurement at HEAD: `tools/validation_matrix.py contract --level L4 --rows`; list what changed in
  `solvers/bem/` and force projection since f1a1a18; decide which B4 verdicts need re-measurement (group regression, not the suite).
- T2 | S2, S3 | worker+SLURM | pending | Spanwise loads: B4 F7 was resolved on 2026-10-28 (32bd208, contract.yaml note: azimuth sampling,
  reader now takes last-revolution means; trusted, e_alpha 4.30 -> 4.15 %, K5 passes). Re-measure e_alpha / mean|dCn| at HEAD for .0-.2
  with the 10-28 recipe (rows are UNMAPPED, so regression cannot read them), restamp if within 0.1 pp, else report.
- T3 | S2, S4 | worker | pending | V-01 (B4 follow-up 2): re-point at OpenFAST on the official deck at rated, or record it unarbitrated.
- T4 | S2 | worker | done (L5) | Zhou et al. 2025 rigid means (16.11 MW, 2.53 MN): matched operating point, our BEM vs Zhou vs AeroDyn.
- T5 | S2, S3 | worker | done (L5) | Matched-options audit: our BEM corrections (tip/hub loss, skew, induction, polars, Viterna) vs the AeroDyn
  deck used as arbiter; tower and wake off on both sides; record any mismatch and its torque/thrust effect.
- T6 | S3 | parent | patterns done (L4); closure pending | Residual patterns for groups 59/60 and the spanwise two-quantity print (B4 follow-up 3); L4 closure or declarations.

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
- L3 (2026-11-14, user choice verbatim): "Aceptar el corrimiento y escribir la nueva base (Recomendado)". Parent: `regression --group 11
  --write` (SLURM 11739002): wrote 18 measurements and the source digest; `check` 307 rows / 426 comparisons, 0 errors, 2 warnings
  (groups 27 and 32 on no level, pre-existing). contract.yaml: 9 group-11 verdicts restamped at a6faa90 with a note (A1 .0-.2, yaw/shear
  .0-.1, A3 polars .0-.2 regression_pin, Viterna like-for-like); verdicts unchanged. Spanwise .0-.2 NOT restamped: UNMAPPED, the
  regression did not read their residuals -> T2. Contract L4 after: trusted 6, regression_pin 3, stale 14 (spanwise 3, V-01 10,
  viterna_manual 1); over ceiling 0.
- L4 (2026-11-14, parent + worker mwh2dzlq-4-pn8i (patterns) + worker mwh4tymc-6-byy5 (test prints; died on quota, edits
  reviewed and finished inline), T6 RESULT, commit 1f0adaf). residual-patterns.json: group 11 reads `%` and ` deg` (spanwise alpha
  and Cn are two printed lines, no print change); groups 59/60 declared. V-01 and Viterna-manual tests now print through
  `assert_residual_below` (same bounds, same strict `<`; Viterna loop: worst residual over 20/30/45/60/90 deg). Local runs: 31 passed;
  `regression --group 59 --scope tests/test_iea15mw_v01_rotor_performance.py` 16 new_baseline then same (rated thrust 2.6708 %, torque
  3.3442 %, equal to B4); `--group 60` 4 new_baseline then same; both written; check 307/426, 0 errors. contract.yaml: 11 verdicts
  restamped at 1f0adaf. Contract L4: trusted 7, regression_pin 3, stale 3 (spanwise, T2), unarbitrated 10 (V-01, T3); over ceiling 0.
  T2 spanwise job 11739029 died in prolog requeues (cancelled), resubmitted as 11739057. Workers for T3-T5 failed on quota.
- L5 (2026-11-14, worker mwh59too-7-ztxn, T4+T5 RESULT; tools `tools/l4_t5_options.py`, `tools/l4_t4_zhou.py`, logs `logs/l4_t4t5/`;
  local runs, OpenFAST 5.0.0, AeroDyn dt 0.25 s / tmax 40 s, no K5 dt-halving control for these single points).
  T5 matched options (rated 10.59 m/s, 7.56 rpm, pitch 0, tilt -6, precone -4): tip/hub loss, TanInd/drag terms, Buhl-Glauert high
  induction, tower off, DBEMT/UA off, density, viscosity (single-Re polars, 0.000 %), nSector (0.000 %) MATCH. Mismatches: skewed-wake
  correction (CCBlade none, AeroDyn Skew_Mod 1): +0.012 % thrust, +0.020 % torque; polar interpolation (CCBlade smoothed spline vs
  AeroDyn table): bounded by a linear-interpolation probe (+0.18 % thrust, +0.71 % torque, moves away from AeroDyn, so not the cause).
  Ours (deck polars) vs AeroDyn at rated: thrust +0.28 %, torque/power +0.39 %; ours (YAML blade, production polars) +1.76 % / +1.40 %.
  FINDING for L5: production configs run tilt 0 / precone 0 (e.g. tests/smoke_fix/base_fix/fluid_yaw_0.yaml:44-45) while the parity
  tests and the deck use -6 / -4: rated thrust +1.38 %, torque/power +2.44 %. Not a BEM defect; the two-way comparison must match it.
  T4 Zhou 2025 Table 6 (10.59 m/s, 7.55 rpm, cone and tilt zero; source odd/tasks/composite-bend-twist-verdict.md:688, 2600-2610;
  pitch/density/shear not stated, assumed 0 / 1.225 / uniform): P/T Zhou 16.11 MW / 2.53 MN; AeroDyn 16.416 (+1.90 %) / 2.535
  (+0.21 %); ours deck polars 16.327 (+1.34 %) / 2.522 (-0.31 %); ours YAML 16.498 (+2.41 %) / 2.560 (+1.18 %). Ours vs AeroDyn at
  Zhou's point: deck polars -0.52 % thrust / -0.54 % power; YAML +0.97 % / +0.50 %. Zhou uses LL-FVW + Beddoes-Leishman: secondary
  reference, not an arbiter; agreement within 2.5 %. Dominant uncertainty: Zhou's unstated pitch (+-0.5 deg -> thrust +-3.2 %).
- L6 (2026-11-14, explorer mwh59tpl-8-vsr4, T3 evidence). V-01 points (V/rpm/pitch): 5.006/5.0/2.893, 7.159/5.10/0, 9.027/6.43/0,
  10.659/7.518/0 (rated rows), 12.259/7.518/6.76, 15.471/7.518/12.19 (tests/test_iea15mw_v01_rotor_performance.py:55-125); none shared
  with group 11 (8/9/12 m/s at 7.56 rpm). `write_aerodyn_dvr` takes (V, rpm, pitch) tuples (tests/support/openfast_bem.py:377-415), so the
  points can be run against AeroDyn as-is. V-01 like-for-like issues vs its WISDEM reference: BEMSolver defaults shear 0.2, tilt 0,
  precone 0 (engine.py:132-141); `_ROTOR_RADIUS` 120.675 vs YAML 121.119 vs deck 120.97 (Cp/Ct normalisation); mu 1.81206e-5.
  Decision (parent, per S4 and explorer): both — add the V-01 points to group 11 against official AeroDyn (arbiter), keep V-01 vs WISDEM
  as a declared unarbitrated secondary record.
