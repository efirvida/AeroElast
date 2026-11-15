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
- T2 | S2, S3 | parent | done (L9) | Spanwise loads: B4 F7 was resolved on 2026-10-28 (32bd208, contract.yaml note: azimuth sampling,
  reader now takes last-revolution means; trusted, e_alpha 4.30 -> 4.15 %, K5 passes). Re-measure e_alpha / mean|dCn| at HEAD for .0-.2
  with the 10-28 recipe (rows are UNMAPPED, so regression cannot read them), restamp if within 0.1 pp, else report.
- T3 | S2, S4 | worker | points added (L7); 2 pitched points red, cause partly located | V-01 (B4 follow-up 2): re-point at OpenFAST on the official deck at rated, or record it unarbitrated.
- T4 | S2 | worker | done (L5) | Zhou et al. 2025 rigid means (16.11 MW, 2.53 MN): matched operating point, our BEM vs Zhou vs AeroDyn.
- T5 | S2, S3 | worker | done (L5) | Matched-options audit: our BEM corrections (tip/hub loss, skew, induction, polars, Viterna) vs the AeroDyn
  deck used as arbiter; tower and wake off on both sides; record any mismatch and its torque/thrust effect.
- T6 | S3 | parent | patterns done (L4); closure pending | Residual patterns for groups 59/60 and the spanwise two-quantity print (B4 follow-up 3); L4 closure or declarations.
- T7 | S2, S3 | worker | implemented (L11), awaiting tests + store decision | (L10) BEM model: blade prebend/cant (CCBlade precurve from the blade reference axis) and exact table
  interpolation of polars (no spline smoothing), test-first; re-measure groups 11 and 59 and the V-01-vs-AeroDyn points; production default.
- T8 | S3 | parent | pending | Store rows/verdicts for the V-01-vs-AeroDyn ids; L4 closure or declarations.

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
- L7 (2026-11-14, workers mwh5fgte-9-swsa (test) and mwh5j1ad-a-r4c8 (diagnosis), T3 RESULT; test `test_v01_points_match_aerodyn`
  in tests/validation/bem/test_bem_openfast_parity.py, logs `logs/check/l4/t3/`, `logs/l4_t3diag/`, tools `tools/l4_t3_{span,polar,pitchshift}.py`).
  V-01 points vs official AeroDyn, A1 setup and A1 source bound 1.5 % (ceiling 5 %), thrust / torque: 5.006/5.0/2.893 +1.370 / +0.426 %;
  7.159/5.10/0 +0.283 / +0.403 %; 9.027/6.43/0 +0.283 / +0.403 %; 10.659/7.518/0 +0.278 / +0.422 %; 12.259/7.518/6.76 +1.864 / +1.007 %
  (RED, thrust); 15.471/7.518/12.19 +3.708 / +3.022 % (RED). K5 dt 1 -> 0.5 s <= 0.005 pp. The two pitched reds are committed as
  strict xfail citing this entry (no bound moved). Diagnosis: (1) pitch sign/axis identical on both sides (alpha agrees <= 0.15 deg);
  (2) CCBlade CCAirfoil smoothed spline (s 0.1) vs AeroDyn table: Cl +0.007..+0.013 outboard in pitched cases; exact/linear
  interpolation moves 15.471 m/s to +2.66 / +1.78 % and rated pitch 0 to +0.47 / +1.09 % (the pitch-0 pass is partly a cancellation);
  (3) residual ~+0.08..+0.15 deg effective alpha offset outboard of r ~80 m (pitch shift that zeroes thrust: 0.076/0.075/0.106/0.149 deg
  at pitch 0/2.9/6.8/12.2), amplified by dT/T per deg falling from -6 %/deg (pitch 0) to -17.8 %/deg (pitch 12.2); origin unlocated
  (candidates: twist/tip interpolation, small outboard a/a' difference); (4) tip node r = Rtip carries load in ours, AeroDyn a = 1 there:
  ~0.1-0.3 pp. Open: locate (3) (per-station phi and twist, tip node removed), then the user decides whether production polars use
  exact table interpolation (moves rated pitch-0 torque by about +0.67 pp vs AeroDyn).
- L8 (2026-11-14, worker mwh5oq87-b-de6e, L7 item (3) LOCATED; tools `tools/l4_t3_locate_{ad,cmp,var,detail,final}.py`, logs
  `logs/l4_t3locate/` (final.log)). Twist, chord, node r, Rhub/Rtip, tilt, precone identical on both sides. The AeroDyn side carries the
  deck blade prebend (BlCrvAC) and its cant angle (BlCrvAng; local cone = precone + cant, -1.5 deg at r 61 m to -5.6 deg at the tip);
  our BEM carries neither (only precone). Effects: AeroDyn Vx falls outboard (10.56 -> 10.46 m/s, pitch 0) vs ours constant 10.575;
  phi +0.023 deg at 61 m to ~+0.06 deg outboard; pitch-dependent twist term pitch*(1/cos(cant) - 1), 0.060 deg at the tip at 12.19 deg
  (form inferred from outputs, not from the Fortran). Mean d_alpha 80-116 m, ours - AeroDyn: pitch 0 +0.056 -> -0.008 deg with
  CCBlade precurve = BlCrvAC (sign -1, empirical); pitch 12.19 +0.111 -> +0.046 (prebend) -> +0.018 deg (+ twist term).
  Thrust/torque, exact linear polars, ours - AeroDyn: 15.471/12.19 baseline +2.66 / +1.77 %, + prebend +1.23 / +0.23 %, + twist term
  +0.81 / -0.22 %; 10.659/0 baseline +0.46 / +1.09 %, + prebend -0.53 / -0.33 %. Tip node (AeroDyn a = 1, negative Fx there) ~0.9 pp
  thrust at 15.471 (tip-removed numbers +-0.5 pp: reference rebuilt from nodal loads). Residual dphi 0.01-0.03 deg unlocated. Sweep
  (BlSwpAC) not consistent, left out. Classification: (a) model-input difference — the parity setup gives AeroDyn the prebent blade
  and ours a straight one; not a defect in the BEM equations. Using prebend in production is a solver change -> user decision.
- L9 (2026-11-14, parent, T2 RESULT + store; SLURM 11739029/11739057 never started (prolog requeues on dev and sequana_cpu,
  cancelled); run locally in 12 s, `logs/check/l4/t2-e-alpha.log`). Spanwise at HEAD, e_alpha / mean|dCn|: dt 1 s 4.30/3.34/4.78 % and
  1.58/1.54/1.71 % (8/9/12 m/s); dt 0.5 s 4.23/3.29/4.74 %; dt 0.25 s 4.15/3.19/4.60 % — identical to 2026-10-28, K5 dt-halving
  <= 0.15 pp on e_alpha. `regression --group 11 --write`: 18 same, 6 new_baseline (spanwise alpha and Cn now read), 6 unclaimed (the new
  V-01-vs-AeroDyn ids have no store row yet); wrote 24 measurements. 12 group-11 verdicts restamped. check 307/426, 0 errors.
  CONTRACT L4: trusted 10, regression_pin 3 (A3 production polars), unarbitrated 10 (V-01 vs WISDEM, declared secondary); stale 0;
  over ceiling 0. Open for L4 closure: (i) store rows/verdicts for the 6 V-01-vs-AeroDyn ids (4 pass, 2 documented red);
  (ii) user decision on prebend/cant and polar interpolation in our BEM (L7-L8); (iii) V-01 unarbitrated rows declared.
- L10 (2026-11-14, user choices verbatim): "Agregar prebend + cant al BEM (Recomendado)"; "Interpolación exacta de tabla (Recomendado)".
  Parent: T7 (solver change, production default, it moves the L5 torque baseline) and T8 (store + closure).
- L11 (2026-11-15, worker mwh64iai-c-arfd / mwh7ncr2-d-fm0i, T7 RESULT; logs `logs/check/l4/t7/`). engine.py: `BEMSolver(...,
  use_prebend=True, polar_interp="spline")`; prebend from `blade_aero.prebend` as CCBlade `precurve = -prebend` (+ precurveTip), only
  when the blade carries prebend; polar modes `spline` (interpolating, s = 0; AeroDyn InterpOrd default 3 = cubic), `linear`,
  `smoothed` (old). Deck blade builder reads BlCrvAC. Sign confirmed from CCBlade Fortran (x_az = -r sin(precone) + precurve
  cos(precone)) and numerically (Vx at r 99.5 m, 15.471/12.19: AeroDyn 15.196, straight 15.349, shipped 15.222, flipped 15.386 m/s).
  Cant twist term NOT added: not documented in the AeroDyn theory manual available (.sources/papers/36881.pdf); residual outboard alpha
  +0.046 deg at pitch 12.19. RED 14 -> GREEN; new tests `tests/validation/bem/test_bem_prebend_polar_interp.py` (8) and
  `test_prebend_outboard_alpha_matches_aerodyn` (bound 0.07 deg on the 80-116 m mean offset: a presence/sign pin between the measured
  straight 0.056/0.111 and prebent -0.008/+0.046 deg, not a physics acceptance). Strict xfail on the pitched V-01 points removed.
  Parity after (thrust / torque vs AeroDyn): A1 8/9/12 -0.684/-0.732, -0.608/-0.488, -0.479/-0.279 %; A2 yaw -0.664/-0.658, shear
  -0.790/-0.254 %; A3 +0.52/-0.72, +0.50/-0.29, +0.31/+0.03 % (was up to +2.15 %); V-01 points -0.378/-0.915, -0.521/-0.322 (x3 at pitch 0),
  12.259/6.76 +0.157/-0.300 %, 15.471/12.19 +1.234/+0.234 % (was +3.708/+3.022 %); spanwise max|dalpha| 0.317/0.445/0.954 deg, mean|dCn|
  1.46/1.59/1.73 %. All within the 1.5 % A1 bound. Production-like rated (YAML blade, tilt -6, precone -4, 10.59 m/s, 7.56 rpm, pitch 0):
  thrust 2.5274 -> 2.4940 MN (-1.32 %), torque 20.334 -> 20.047 MN m (-1.41 %), power 16.10 -> 15.87 MW (prebend -1.02 / -1.50 %,
  polars -0.32 / +0.04 %). V-01 vs WISDEM: rated thrust 2.671 -> 2.104 %, torque 3.344 -> 2.875 %. Regression (no --write): group 11
  22 changed / 2 same / 12 unclaimed; group 59 16 changed; group 60 same. check.sh quick OK.
