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

## Log

- L1 (2026-10-26) User: S1 verbatim (during the PR of l3/twist-p1, merged as #68 / 3a2d222).
- L2 (2026-10-26) Parent reliability summary: bending and frequencies within 2.4 % (IEA-15) and 4 % (IEA-3.4) of the
  BeamDyn 6x6; torsion not reliable (IEA-15 S-7 7 % soft; 3.4 outboard FE 1.3x its own section; KO2017 twisted beam
  and D-Tube element reds); mode shapes unmeasured (T8). Recommendation: core torsion first, then T8.
- L3 (2026-10-26) User: S2 verbatim.
