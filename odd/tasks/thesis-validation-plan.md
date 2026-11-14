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
- T3 | S2 S3 | parent+workers | in progress (odd/tasks/l4-bem-rigid-closure.md) | step 3: BEM vs AeroDyn (rigid), matched options
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
