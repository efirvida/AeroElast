# Core torsion reliability (shell FE torsional stiffness)

Gate before continuing L3 (T6/T7/T8) and L4/L5: decide whether the shell FE (MITC4) or the twist measure is soft
in torsion on thin, twisted or tapered sections, and fix it if it is the core. Origin: `odd/tasks/l3-closure.md` T9,
Handoff 4. Branch `l3/core-torsion` from main 3a2d222.

## Specs

- S1 (user, verbatim): "recuerda que no es que yo quiera arreglar el modelo de una u otra pala, yo quero que aeroelast sea cientificamente confiable para el proposito que te di"
- S2 (user, verbatim, after the reliability summary): "entonces sabes ien lo que tienes que corregir antes de poder seguir"
- S3 (repository physics discipline): every number with a metric, an independent reference and a written tolerance;
  prediction written before each run; no pass band widened; a moved anchor is a finding.
- S4 Done means: the IEA-3.4 outboard softness (FE twist 1.2-1.5x the integral of T/GJ of its own mesh section at
  eta 0.65-0.85) is attributed to element or twist measure with closed-form evidence, and, if element, fixed with
  the element benchmarks (D-Tube, KO2017 3.5 twisted beam, thin tubes) and S-7 / S-4 re-measured.

## Tasks

- C1 | S3 | parent | pending | Baseline: full suite on main 3a2d222 (SLURM), compare with known-red-full at a9ca498.
- C2 | S3 S4 | parent | pending | Collect the three symptoms with their current numbers: D-Tube torsion vs Bredt, KO2017 3.5 twisted beam (MITC4 thin, 8.5 % vs 1 %), 3.4 outboard 1.3x (`tools/l3_twist_integral.py`); IEA-15 S-7 1.072. Prediction of what an element defect vs a distortion measure would show in each.
- C3 | S3 S4 | parent/worker | pending | Split element vs measure on the 3.4 at 0.65-0.85: rigid-fit residual per station; twist fitted on the box nodes (caps + webs) vs all nodes; element quality of the outboard rows (aspect, warp) and es 0.15 vs 0.10.
- C4 | S3 S4 | parent/worker | pending | Tapered closed-form tube through MITC4 at es 0.15 with the 3.4 outboard perimeter, wall G t and taper (extend `tools/diagnose_s1_gj_tube.py` / `tests/test_s1_gj_tube.py`), RED first if it fails.
- C5 | S4 | parent | pending | If element: fix (test first), then element benchmarks, S-7, S-4 on both blades; if measure: fix the twist extraction and re-measure S-1 GJ / S-7.

## C2 result (main 3a2d222, login, 54 s; `logs/check/core-torsion/c2-benchmarks.log`, 14 passed)

| benchmark | measured | bound |
|---|---|---|
| D-Tube GJ vs Bredt-Batho closed cell (`test_dtube_torsion_matches_bredt`) | 1.33336e11 vs 1.34554e11, -0.91 % | 2 % |
| KO2017 3.5 twisted beam MITC4, t 0.02667 / 0.0002667, in / out of plane | 0.10 / 0.25 / 0.06 / 0.01 % | 1 % |
| Thin twisted beam vs Simo, Fox & Rifai 1989 | 0.24 % | 5 % |
| S-1 tube truth tests (rectangle, circle, oblique, S-7 slices, torque helper) | 7 passed | per test |

**Correction to L2:** the "element reds" (D-Tube, KO2017 twisted beam 8.5 %) are not red on main; the AGENTS.md
note dates from 890d4e6. The element is right in torsion on prismatic thin closed sections and on a twisted
thin strip. What none of them has: **taper** (perimeter and wall changing along the span), a multi-cell section
with webs, and the distorted outboard quads of a real blade. The 3.4 softness has to come from one of those or
from the twist measure (C3, C4).

## C3 prediction (written before the run)

On the IEA-3.4 at eta 0.65-0.85, same solve as `tools/l3_twist_integral.py` (T40 shear-flow tip torque):
- **Measure (section distortion)**: the twist fitted on the box nodes only (elements of the spar caps `*SPAR*` and
  webs `*_SW`) gives an increment within 0.95-1.05 of the mesh-Bredt integral while the all-node fit stays at
  ~1.3; the rigid-fit residual (RMS of u - (u0 + theta x d) over |theta| RMS |d|) grows outboard above 5 %.
- **Element**: the box-only ratio is also >= 1.2 and the residual stays below 2 %.
- Mixed (box 1.05-1.2): both; then C4 decides the element share.
Control IEA-15 at 0.3-0.85: box and all-node ratios within 0.05 of each other, residual < 2 %.

## C3 result (worker; jobs 11625399 3.4, 11625400 15; `logs/check/core-torsion/c3-*.log|json`; `tools/l3_twist_integral.py --c3`)

Box = element sets `*SPAR_A/B/C` (caps) and `*_SW` (webs). FE delta-theta / mesh-Bredt integral, all nodes -> box:

| blade | 0.3-0.6 | 0.6-0.85 | 0.7-0.85 |
|---|---|---|---|
| IEA-3.4 | 1.014 -> 1.083 | 1.314 -> **0.926** | 1.380 -> **0.911** |
| IEA-15 | 1.002 -> 0.971 | 0.952 -> 1.018 | 0.932 -> 1.032 |

Rigid-fit residual (all / box), 3.4: 0.7 0.13 / 0.045, 0.75 0.14 / 0.034, 0.8 0.24 / 0.020, 0.85 0.36 / 0.051
(inboard 0.46-0.69 / 0.04-0.08, small |theta| in the denominator); 15: 0.04-0.10 / 0.009-0.04 at 0.7-0.85.
Elements benign (median aspect 1.2-2.1, warp < 1 deg) except the thin TE reinf/flat strips (max aspect 31-84).
Verdict: **measure** (prediction held in kind; box 0.91-0.93 is below the 0.95-1.05 band, and the 15 control missed
its 0.05 / 2 % limits: distortion is present on both blades, larger on the 3.4). The FE box is not soft; the
section distorts in its plane, and an all-node rigid fit reads that as twist. Open: is the distortion physical
(non-prismatic: taper, pretwist, prebend; no diaphragms) or an artefact; and which twist the BEM feedback uses
(`fsi_participant.py _compute_deformed_geometry`).

## C4 (revised) prediction (written before the run)

Prismatic extrusion of the real IEA-3.4 section at eta 0.8 (a copy of the yaml, under `logs/check/core-torsion/`,
with chord, airfoil/relative thickness, pitch axis and every layer and web value held at their eta-0.8 value over
the whole span; twist, prebend and sweep zero), same mesher, es 0.15, same solve. Saint-Venant/Bredt is exact for
a prismatic tube away from the ends. Predicted, at eta 0.3-0.7 (away from root and tip):
- all-node and box ratios both within 0.97-1.03 of the mesh Bredt, residual all-node < 3 %: then the 3.4
  outboard distortion is a **non-prismatic physical effect**, element and measure are right, and the open item is
  which twist definition the aerodynamics and the beam comparison must use;
- all-node ratio >= 1.1 or residual >= 5 % on the prismatic: the distortion is produced by the section itself
  (shell model of the TE strips, load path), an artefact to fix in the model or element.

## C4 result (worker; job 11625404; `logs/check/core-torsion/prismatic/`; `tools/core_torsion_prismatic_yaml.py`)

Prismatic copy of the 3.4 at eta 0.8 (chord 1.9, pitch axis 0.308, twist 0, no prebend/sweep, DU91/DU08 blend
w 0.5576 -> t/c 0.2277, all 17 layers and both webs held at their eta-0.8 values): mesh Bredt GJ 9.561e6 at all 100
stations; FE / mesh Bredt **0.982 (all nodes) = 0.982 (box)** at every slice 0.25-0.85 (0.989 at 0.9, tip);
residual 0.15-0.45 % (all and box). **First reading holds:** element and twist measure are right on the real
multi-cell section; the 3.4 outboard distortion (C3) is a non-prismatic physical effect (taper, pretwist, prebend)
that a rigid-section beam cannot represent. The uniform 1.8 % stiffer-than-Bredt is the thick-wall / open-wall term
(T27 measured +2.0 pt on the IEA-15), not distortion. TE strips aspect 49-51 do no harm.
**Consequences.** (1) The core is reliable in torsion: no element fix (C5 as written does not apply). (2) Against
the beam references (S-1 GJ, S-7) the beam-equivalent twist is the box twist; the all-node twist carries physical
section distortion. (3) Production BEM feedback (`fsi_participant._section_rotation`: antisymmetric part of an
affine fit over the whole ring) includes the non-affine distortion; the aerodynamic twist is the chord-line
rotation. Which one the BEM must take is a modelling decision, to be made on its measured torque effect.

## C5 (revised): twist definition sensitivity; prediction (written before the run)

IEA-15, rated, one-way, same solve as `tools/l3_torque_sensitivity.py` (rigid production loads, es 0.25): per BEM
strip, elastic twist by (a) the production ring-affine rotation, (b) the chord line (LE and TE nodes of the ring),
(c) the box rigid fit; rigid-geometry BEM torque/thrust with each. Prediction: (b) and (c) within 5 % of (a) in
tip elastic twist (the 15's all-node residual is 4-10 % at 0.7-0.85), so torque within 0.15 % of (a) (twist alone
moves torque -3.2 %). Refuted if any definition moves torque > 0.3 % or tip twist > 10 %: then the definition is a
thesis-level choice and goes to the user with the numbers.

## Log

- L1 (2026-10-26) User: S1 verbatim (during the PR of l3/twist-p1, merged as #68 / 3a2d222).
- L2 (2026-10-26) Parent reliability summary: bending and frequencies within 2.4 % (IEA-15) and 4 % (IEA-3.4) of the
  BeamDyn 6x6; torsion not reliable (IEA-15 S-7 7 % soft; 3.4 outboard FE 1.3x its own section; KO2017 twisted beam
  and D-Tube element reds); mode shapes unmeasured (T8). Recommendation: core torsion first, then T8.
- L3 (2026-10-26) User: S2 verbatim.
- L4 (2026-10-26) User: "continua en esta sesion". Parent: C2 (benchmarks green), C3 (measure), C4 (prismatic: core right).
