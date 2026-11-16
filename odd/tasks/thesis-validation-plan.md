# Thesis validation plan: structural model + simple BEM, one-way and two-way coupling

What the PhD needs from aeroelast: a defensible, component-partitioned validation of the shell structural model
and of its one-way and two-way coupling with the built-in BEM. Tower and wake effects are out of scope here
(they come later from the CFD coupling, another project).

## Specs

- S1 (user, verbatim): "el problema es que nadie va a creer en mis resultados si no presento una validacion cientifica del mismo que  pueda ser aveptada. Lo que quiero estudiar de forma general es medir como influye la deformacion de las palas cuando el aerogenerador produce energia sobre la se~nal de torque del mismo, incluyendo cuando la pala pasa frente a la torre. Pero de momento para validacion no hay muchos datos por lo que el proceso de validacion tiene que ser particionado por componentes, lo unico que he encontrado para validar el es articulo de Zhou."
- S2 (user, verbatim): "el efecto de la torre lo medire cuando tenga los dos codigos acoplados, nuestro trabajo es validar sobretodo el modelo estructural contra un BEM simple y que se acople bien unidireccionalmente y bidireccionalmente, luego el efecto de la estela y la torre se calcula cuando tenga un CFD"
- S3 (repository rule): every compared number has an independent reference and a written tolerance (5 % relative ceiling); a red that stays red is reported with its cause and, where it can move the thesis quantity, with its measured effect on rotor torque and power (sensitivity), never hidden by widening.

## Plan (V&V hierarchy: code verification -> components -> coupled system)

**Strict order (user L3).** Steps 1 -> 2 -> 3 -> 4 -> 5 are validated in this order and no step starts before
the previous one is closed (its rows trusted, or every remaining red declared with cause and torque/power
sensitivity). Work outside the current step needs the user's explicit go-ahead.

| step | level | what | reference | acceptance |
|---|---|---|---|---|
| 1 | L1, L2 | elements and linear solver | analytic and published benchmarks, CalculiX | `contract` closed (odd/tasks/l1-l2-closure.md) |
| 2 | L3 | blade structure, focused on what moves torque: elastic twist (GJ, S-7, twist under load), edgewise (EI_edge, 1E), flapwise (EI_flap, static deflection, 1F); mass | BeamDyn 6x6, OpenFAST v5 modes, Zhou et al. 2025 (Energy 336:138488) twist | 5 %; every remaining red gets its torque/power sensitivity measured and is declared as a limitation (EA: no torque path, declared) |
| 3 | L4 | built-in BEM, rigid blade | OpenFAST AeroDyn with the same corrections, tower and wake models off; Zhou rigid means (16.11 MW, 2.53 MN) | 5 % on power, thrust, torque, spanwise loads |
| 4 | L5 one-way | BEM loads -> structure, no feedback | OpenFAST BeamDyn under the same loads (S-2 type), S-5 rated run | 5 % on deflection and twist; S-6 NaN closed |
| 5 | L5 two-way | deformation feeds back into BEM | OpenFAST BeamDyn + AeroDyn, steady rated wind, tower off; Zhou flexible-blade twist, power, thrust | 5 % on mean torque/power/thrust and tip twist; **the flexible-minus-rigid torque and power deltas** compared between the two codes (the thesis quantity) |
| - | L6 | preCICE + OpenFOAM, tower, wake | CFD project | out of scope here |

## Tasks

- T1 | S2 S3 | worker | done (L1 157/2, L2 79/4 CLOSED; l1-l2-closure.md L3) | step 1 (odd/tasks/l1-l2-closure.md)
- T2 | S1-S3 | parent+workers | done (L4; odd/tasks/l3-closure.md, odd/tasks/l3-outboard-section-stiffness.md L30-L35) | step 2: list the L3 rows per torque path, decide close vs declare, sensitivity runs
- T3 | S2 S3 | parent+workers | done (L5; odd/tasks/l4-bem-rigid-closure.md L14) | step 3: BEM vs AeroDyn (rigid), matched options
- T4 | S2 S3 | - | pending | step 4: one-way coupling
- T5 | S2 S3 | - | pending | step 5: two-way coupling vs OpenFAST and Zhou, flexible-minus-rigid deltas

## Log

- L1 (2026-10-26, session 9) User: S1 (full message in odd/tasks/l1-l2-closure.md Log L2).
- L2 (2026-10-26, session 9) User: S2. Parent: tower and wake removed from this plan; OpenFAST comparisons run with tower influence off.

- L3 (2026-10-26, session 9) User: "recuerda esto como un orden estricto que vamos a validar. y cierra la sesion". Parent: strict-order rule above; pointer added to AGENTS.md.
- L4 (2026-11-14) User: "cerrar los pendientes antes de L3 Luego Cerrar L3, y comienza con L4". Parent: step 2 (L3) CLOSED with
  declared limitations, each with its torque/power sensitivity (all below ~1 %): D1 outboard GJ, deck TE adhesive closure (a), torque/power
  +0.54 %, thrust +0.90 %; D2 outboard EI_edge, same deck provenance (a), torque -0.006 %; D3 tip end-region compliance (b), torque/power
  <= +0.27 %, thrust <= +0.37 %; z 90-92 distortion (b), torque/power <= +0.154 %; D4 S-7 metric fixed (c), real-mesh test committed;
  EI_edge CLT mean +0.4 % (hoop-free), the old -0.40 % superseded. Yaml used as written (user decision, l3-outboard L27). The
  production-twist estimator (S8) stays deferred to L5. Evidence: odd/tasks/l3-outboard-section-stiffness.md L30-L35. Step 3 (L4) opened:
  odd/tasks/l4-bem-rigid-closure.md.
- L5 (2026-11-15) User: "1". Parent: step 3 (L4) CLOSED: our rigid BEM matches official OpenFAST 5.0.0 AeroDyn (tower and wake off,
  matched options) within 1.5 % on thrust and torque at every arbitrated point, including pitched V-01 points after adding blade
  prebend and exact polar interpolation (production rated torque -1.41 %, power 16.10 -> 15.87 MW); spanwise alpha/Cn within bounds;
  Zhou rigid means within 2.5 % (secondary). Declared: cant twist term absent, V-01 vs WISDEM unarbitrated (contract shows L4 open for
  them). For step 4/5: production tilt/precone must match the deck (-6 / -4) or the difference be declared (+2.44 % torque). Pending:
  T7 blade-test check (SLURM 11739091). Next: step 4 (L5 one-way) only with the user's go-ahead.
- L7 (2026-11-16) Parent: step-5 OpenFAST reference decks run on the login node (60 s, $SCRATCH/l5_twoway/openfast/{flex,rigid}/run.out).
  Rigid twin: steady, P 15.73 MW, thrust 2.476 MN, torque 19.87 MN m, tip TDxr 0.020 m (x1e3 stiffness is rigid enough;
  consistent with the L4 closure, production 15.87 MW). Flex: NOT a usable reference. Power 7.4 +/- 6.0 MW (5 s window std),
  thrust -23 %, power -53 % vs rigid, tip flap 11.7 +/- 0.9 m, tip WM rotation z -0.05 +/- 0.30 (about +/-17 deg) sustained
  from t = 10 s to 60 s, 1.6 %/rev drift: a sustained torsional/limit-cycle oscillation (numerical instability at DT 0.01 with
  BeamDyn after the v4.2->v5 conversion, or a deck defect) -> diagnose before any step-5 number (DT/DTBeam, rhoinf, damping,
  the converted BeamDyn blocks). Login limit: 1800 CPU-s per process; the flex 60 s took 1436 CPU-s single-threaded, so a
  smaller DT needs SLURM or a restart split.
- L8 (2026-11-16) User: "si" (diagnose the flex deck on the login node). Parent, 12 s variants ($SCRATCH/l5_twoway/openfast/diag/):
  oscillation 3.62 Hz (power, tip flap, tip torsion; blades in phase; 1P 0.126 Hz) independent of DT 0.01/0.005/0.0025 and of
  NumCrctn 1 (not numerical). No aero (CompAero 0): steady (tip flap 1.134 m, WM rz -0.0086, std ~1e-16). UA_Mod 3 + DBEMT 2:
  P 14.58 +/- 0.21 MW. Structural damping mu x10: P 15.34 +/- 0.04 MW. Verdict: an aeroelastic torsional instability of the
  converted deck with quasi-steady aero (UA 0, DBEMT 0) and the deck's light stiffness-proportional damping (mu 0.003/0.002/
  0.0008) -- quasi-steady lift has no lag and pumps the ~3.6 Hz torsion-flap mode. Next: 60 s with UA 3 / DBEMT 0 (steady
  state converges to the static polars, as our quasi-steady BEM) and with mu x10, to check that the reference MEANS agree.
- L9 (2026-11-16) Parent: 60 s robustness check: flex UA_Mod 3 / DBEMT 0 vs UA 0 with BeamDyn mu x10, last revolution:
  P 15.351 vs 15.344 MW (0.05 %), thrust 2.1857 vs 2.1860 MN, torque 19.391 vs 19.380 MN m, tip flap 14.582 vs 14.583 m,
  tip WM rz -0.06478 vs -0.06479 -> the reference means do not depend on the cure. Generator switched to UA_Mod 3 for both
  decks (tools/l5_twoway_openfast_deck.py AD_OVERRIDES, comment cites L7-L8); decks regenerated and rerun (login, 60 s, flex
  ~100 CPU-s once stable): flex P 15.351 MW, thrust 2.18573 MN, torque 19.3908 MN m, tip flap 14.58 m, edge -1.324 m, tip
  rotation -3.711 deg (WM -0.06478), drift last/prev rev <= 0.007 %; rigid P 15.74 MW (UA 0: 15.73; 4-digit output),
  thrust 2.476 MN, torque 19.88 MN m, tip 0.020 m. OpenFAST flexible-minus-rigid: power -2.47 %, torque -2.46 %, thrust
  -11.72 % (the step-5 thesis quantity, OpenFAST side; reference data only, step 4 not closed).
  logs/check/l5_twoway/openfast_means.json.
