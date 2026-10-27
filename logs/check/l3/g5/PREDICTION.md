# T12 group 5 — predictions (written before any 0.5 / 0.35 / 0.25 m run, 2026-10-27, session 17)

Worktree `wt/l3-g5`, branch l3/g5-mode7 on 4b29bc2, test diff uncommitted at submission (residuals print before
the monotonicity assertion; no bound moved). Plumbing only, NOT evidence: `diagnose_g5_mode7.py --meshes 1.0`
on the login node paired modes 1-7 identically by frequency and by MAC (MAC 0.85-1.00); CCX 7 edgewise (x 0.90).

Context: T29 (SLURM 11623984) read CCX mode 7 as 0.89 / 1.01 / 1.01 % (unsigned) and mode 8 as 13.19 / 7.97 /
3.83 %; the test docstring table (SLURM 11622478, before Hungarian pairing) had mode 7 +0.95 / +0.73 / +0.52 %.

- P1 (H1, pairing): frequency and MAC pairing agree for CCX modes 1-7 on all three meshes, MAC >= 0.8. If they
  disagree for mode 7 on any mesh, the "non-convergence" is a pairing artefact and the test must pair by shape.
- P2 (H2, character): CCX mode 7 is the second edgewise mode (x share >= 0.85 of the squared translations), not
  the torsion-dominated mode, so the core-torsion softness does not explain it.
- P3 (H3, plateau): if P1 holds, the signed mode 7 gap is positive and flat within 0.2 pp from 0.35 to 0.25 m;
  a positive plateau of about +1 % is an element-order / modelling difference to explain, not discretisation.
- P4: CCX mode 8 keeps MAC < 0.8 with every AeroElast mode at 0.5 m and its MAC rises with refinement: the
  3.83 % is a gap between partially-paired shapes, so it is not a clean residual either.
- P5 (test run): both residual lines print; low worst < 0.2 %; high worst about 3.8 % (< 5 %); the test still
  FAILS on the monotonicity assertion for CCX mode 7 (unchanged physics, only the order of the checks moved).

## Jobs
- 11625834 diag (sequana_cpu, 8 tasks): logs/check/l3/g5/diag-11625834.log, mac.json
- 11625835 test -s (sequana_cpu_bigmem, 8 tasks): logs/check/l3/g5/test-11625835.log

## Note added after submission, before any result (2026-10-27)
The cluster `ccx` resolved by `ccx_bin_or_skip()` is `/scratch/leahk/eduardo.donestevez/venv/bin/ccx`, **Version 2.20** (T13 audit). The test docstring and store label say 2.23, and issue #8 measured the high-mode worst gap 4.03 % under 2.20 vs 2.12 % under 2.23. Both jobs above, and T29 11623984, run 2.20: a mode 7/8 reading must be checked against the CCX version before it is called element physics.
