# L3 v2 closure: diagnose, fix or declare every blade residual above 5 % (gmsh N20, official windIO v2)

L3 was closed on 2026-11-14 with D1-D4 declared (thesis-validation-plan.md:42-47). That closure was measured on the
earlier meshes and revision. The v2 numbers of L-S86 (2026-11-23) postdate it and nobody re-declared them. This plan
re-closes L3 on the current pipeline in one working day, so L5 (FSI) can start. It does not reopen settled physics.

## Specs

- S1 (user, verbatim, 2026-11-24): "haz un plan detallado para cerrar estos gams de una vez que llevamos un mes y no logramos llegar a la validacion FSI entonces escribe un plan para cerrar mejorar o justificar estas diferencias moyores de 5% diagnosticarlas resolverlas si tienen solucion o justificarlas(documentarlas) si fueran diferencias fisicas de modelo. y poder avanzar, HOY tenemros que resolver esto para ir a FSI"
- S2 (user, verbatim, 2026-11-24): "nota que te hablo de todos los que en la tabla anterior supera el 5% de diferencia"
- S3 (scope, derived from S2: the rows of the v1-vs-v2 table above 5 % on v2, L-S86). Each one ends either *fixed* (below 5 %, with
  an independent reason) or *declared* (cause, independent evidence, torque/power/thrust sensitivity, verdict):

  | id | residual (v2 gmsh N20) | band | cause already diagnosed | where |
  |---|---|---|---|---|
  | R1 | S-1 GJ 12.58 % | 5 % | D1: deck carries a TE adhesive wedge the yaml and the shell do not (verdict a). Plus an undiagnosed jump 10.85 -> 12.58 % when ring recovery moved from z strips to gmsh rows (L-S72) | l3-outboard-section-stiffness.md L23, L25, L29, L30; gmsh-official-mesher.md L-S72 |
  | R2 | S-7 torsion ratio 1.086 (8.60 %) | 5 % | same D1 deficit seen globally; metric fixed (D4) | l3-outboard-section-stiffness.md L20, L25, L30, L34 |
  | R3 | S-1 EI_edge 8.57 % | 5 % | D2: deck TE filler stiffness; shell and three section solvers agree outboard (verdict a) | l3-outboard-section-stiffness.md L21, L25, L30, L32 |
  | R4 | S-1 EA 7.00 % | 5 % | deck under-represents the yaml spar cap ~9 % at eta 0.50 (SONATA layer-end rounding), verdict a; plan says "EA: no torque path, declared" | thesis-validation-plan.md:22; l3-outboard-section-stiffness.md L19, L25 |
  | R5 | V-02 1F 7.17 %, 1E 6.88 % | 5 % | reference problem: NREL/TP-5000-75698 is the v1.0 blade, the yaml/decks are v1.1; dropped from L3 2026-10-26 ("Darla de baja (recomendado)"); mesh-converged (L-S80) | l3-closure.md L12 (2); levels.yaml:98-100 |
  | R6 | S-2 edgewise tip ratio: v1 10.26 %, **v2 value not read** | ? | unknown on v2 | to measure in T0 |

- S4 (physics discipline, AGENTS.md 1-8): no band is widened; a residual is *declared* only with a named cause, an
  independent cross-check and a measured sensitivity on the current revision; negative results are reported. The 5 %
  ceiling (levels.yaml:13-14) and the closure rule (thesis-validation-plan.md:22: "every remaining red gets its
  torque/power sensitivity measured and is declared as a limitation") are the criterion, written before measuring.
- S5 (standing decision, l3-outboard-section-stiffness.md L27, verbatim): "Quedarnos con el yaml tal cual (recomendado)"
  -- the shell follows the yaml; the deck's TE closure is not added. This plan does not re-open it unless the user does
  (decision D-A in T6).
- S6 (time box): one day. Every task has its budget; a task that overruns stops and goes to the user as a decision,
  it does not grow.

## Tasks

| id | S | route | budget | status | commit |
|---|---|---|---|---|---|
| T0 | S3 | inline | 30 min | pending | |
| T1 | S3 R1, S4 | worker, SLURM | 2 h | pending | |
| T2 | S3 R1-R4, S4 | SLURM (3 jobs, 2-3 min each) | 1 h | pending | |
| T3 | S3 R4 | inline | 30 min | pending | |
| T4 | S3 R5 | inline | 20 min | pending | |
| T5 | S3 R6 | inline/SLURM | 30 min | pending | |
| T6 | S3, S4, S5 | user decision | 15 min | pending | |
| T7 | S1 | inline | 1 h | pending | |
| T8 | S1 | inline | 30 min | pending | |

### T0 -- baseline on HEAD (30 min)

- `tools/validation_matrix.py contract --level L3 --rows` and the same for L1, L2, L4: which verdicts are stale and why.
- Read the v2 S-2 rows (group 39) for the edgewise ratio (R6) and the S-1 rows (group 38) at HEAD.
- Record revision, environment, SLURM ids in the Log. Output: the R1-R6 table with HEAD values, nothing else.

### T1 -- is the GJ jump 10.85 -> 12.58 % harness or physics? (2 h, the only open diagnosis)

L-S72 changed how rings are recovered (z strips -> gmsh rows) and GJ moved +1.7 pp with no physics change. Also open:
"row z-spread up to 0.15 m exceeds ~0.1 m tip ring pitch".
- Same mesh, same revision, three GJ extractions per station: (a) gmsh rows (current harness), (b) z strips (old
  harness), (c) the independent Bredt/St-Venant evaluation on the same yaml section (l3-outboard-section-stiffness.md).
  Mid-span (eta 0.3-0.7) and outboard (0.8-0.95) separately.
- Acceptance written now: if (a) and (b) differ by > 1 pp where (c) sides with one, the other is a harness defect ->
  fix the harness, re-measure S-1 (SLURM, --ntasks only). If (a) ~ (b) within 1 pp at stations and the 12.58 % mean
  comes from which stations enter the mid-span mean, document the mean definition and keep 12.58 %.
- Either way R1's *cause* stays D1; T1 only decides which number is the true residual.

### T2 -- re-measure the D1/D2/EA sensitivities on HEAD (1 h)

The L22 sensitivities (+0.54 % torque, +0.90 % thrust for GJ; -0.006 % for EI_edge) were taken at revision 1f6ab2b,
before the centrifugal K_G static-solve fix (T17, +10 % flap rotating), the BEM prebend/polar change (L4 T7), the
pretwist sign fix and the gmsh node weld. They are not citable for v2.
- Re-run `tools/l3_w6_scale.py` + `tools/l3_w6_run.py` (L22 recipe: rated one-way, 10.59 m/s, 7.56 rpm, pitch 0, gmsh
  chord-scaled N20 AR 1) at HEAD: base, GJ scaled to the deck (V1), EI_edge scaled (V2), EA scaled (V3 -- new, to give
  R4 a number instead of "no torque path"). SLURM `--ntasks=1`, `sequana_cpu_dev`, 2-3 min each.
- Output: torque, power, thrust, tip flap/edge/twist deltas per variant vs base. Acceptance (from the closure rule):
  a residual is declarable when its sensitivity on torque/power is below 1 %; above 1 % it goes to T6 as a decision.

### T3 -- EA (R4): independent evidence (30 min)

- Already measured (L-S83): layup mass 67 740 kg vs yaml 67 015 kg (+1.08 %); thin-wall EA vs deck +2.9..+10.4 % over
  eta 0.25-0.8. Pair it with T2's V3 sensitivity and with S-4 rated (centrifugal stiffening path, 1F 1.77 %, passes):
  EA's only route into the rotor response is the centrifugal stretch, and S-4 measures it.
- Output: the R4 declaration paragraph.

### T4 -- V-02 (R5): close as a reference limitation (20 min)

- No physics work: the reference describes a different blade revision (v1.0 vs v1.1) and the rows are already off L3.
  The substitute evidence for the same quantity is S-4 parked vs the official v1.1 deck: 1F 2.24 %, 1E 1.25 %, 2F 0.04 %.
- Output: R5 paragraph in the closure; group 32 keeps its rows on no level (the `check` warning stays, explained).

### T5 -- S-2 edgewise ratio on v2 (R6) (30 min)

- If T0 finds it below 5 %: nothing to do. If above: it is the static face of R3 (EI_edge, D2) -- show it with the
  T2 V2 variant (does scaling EI_edge to the deck bring the edge tip to the beam?). If V2 closes it, it is declared
  with D2; if not, it is a new finding reported to the user, not absorbed.

### T6 -- user decisions (15 min, the only stop)

Present the R1-R6 table with: residual at HEAD, cause, cross-check, sensitivity (T2), proposed verdict. Two decisions:
- D-A: accept D1/D2 as deck-model differences again on v2 (S5 holds) **or** reopen S5 and add the TE adhesive to the
  shell (l3-closure.md T7, predicted +5..7 pt GJ at eta 0.7-0.9: closes R1/R2 and moves R3; cost 1-2 days, outside today).
- D-B: accept the L3 v2 closure and open L5.

### T7 -- write the closure (1 h)

- `docs/validation_closures.md`: L3 v2 section, one row per R1-R6 (number, cause, evidence path, sensitivity, verdict).
- `odd/tasks/thesis-validation-plan.md`: L3 re-closed on v2 with the date, revision and the user's verbatim acceptance.
- Store: `gaps.yaml` / `levels.yaml` limitations for D1, D2, EA, V-02 with the T2 numbers; `contract` verdicts restamped
  only where the residual is `same`; `tools/validation_matrix.py check` 0 errors.
- Thesis text: one paragraph per declared limitation (what differs, why, how much it costs on torque/power).

### T8 -- open L5 one-way (30 min)

- Confirm L1, L2, L4 contract status (T0) so L5's implicit `requires` (all lower levels, levels.yaml:18-19) holds.
- Next session starts at l5-one-way.md T1 (BeamDyn under the S-5 AeroDyn loads) and the S8 twist-estimator decision.

### Out of today (recorded, not done)

- Store debt (no number moves): groups 68-70/72 residual patterns, group 7 id, 15 stale re-derive groups + group 60
  target (L-S86, L-S88). After L5 starts.
- l3-closure.md T6 (cap chamfer) and T7 (TE closure): only if D-A reopens S5.
- `tools/l5_s1_ea_gj_diag.py` sweeps `GMSH_STRIP_M`, which no longer has an effect (L-S72/L-S75).

## Log

- L1 (2026-11-24) User request, verbatim: see S1 and S2.
- L2 (2026-11-24) Evidence map compiled from l3-outboard-section-stiffness.md, l3-closure.md, thesis-validation-plan.md,
  levels.yaml, gmsh-official-mesher.md L-S72/L-S80/L-S83/L-S86 (explorer mwuldbd2-2-fb4w). Gap found: the 2026-11-14
  L3 closure predates the v2 numbers; no user acceptance of the v2 values exists; the L22 sensitivities predate T17
  (centrifugal K_G), L4 T7 (BEM prebend/polars), the pretwist fix and the node weld.
