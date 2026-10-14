# Feature: UL elastica rejected rows - register converged evidence, resolve the alpha=2 defect

Worktree `fem-shell-worktrees/rotor-frame`, branch `fix/bem-participant-rotor-frame`. Follows PRE-2 of `odd/tasks/validation-l3-b3-audit.md` (3 L1 rows `rejected`: G33 UL alpha=1 and 2, G22). Diagnosis in Engram `fem-shell/validation/ul-elastica-rejected-rows-diagnosis` (obs 544). No P1 execution, no public GitHub mutation, no push.

## Specs

- S1: "haz las dos cosas" - (a) register the converged-mesh evidence for G33 alpha=1 and G22 and correct the stale docstrings, (b) investigate the alpha=2 defect.
- S2: "investiga el defecto de α=2 hasta resolverlo" - find the cause of the mesh- and step-converged +5.29 % (w) / +13 % (u) alpha=2 gap and fix it, or conclude with evidence that it is not a defect.
- S3: No pass band is widened to make a row green (AGENTS.md physics discipline); a moved anchor is reported.
- S4: Mandatory environment: `PYTHONPATH=$SCRATCH/validenv:$PWD scripts/aeroenv.sh python ...`; verify the import first. SLURM only `-p sequana_cpu --ntasks=...`.

## Tasks

- [ ] T1 (S1, S3) Register converged-mesh rows for G33 alpha=1 and G22; fix docstrings (6.0 %, "converges ~4 % below"). Commit.
- [ ] T2 (S2) Independent arbiter for alpha=2: CalculiX NLGEOM S8R dead tip load, same plate, paired convergence. Commit.
- [x] T3 (S2) Cause established with a 2-D corotational beam (`$SCRATCH/studies/beam_scheme.py`): tangent flow (K_mat on current geometry, no K_sigma, no equilibrium) converges AWAY from the elastica (alpha=1 w/L 0.31576 vs 0.30172; alpha=2 0.54916 vs 0.49346, u/L -0.2054 vs -0.1606) while Newton on the consistent tangent gives 0.30173 / 0.49349. Commit: see T4.
- [x] T4 (S2, S3) Opt-in `scheme: newton` (total-Lagrangian Newton on fint/kt, planar-valid), UL stays the default (blade twist is multi-axis). Narrow strip b=0.025 nx=64: alpha=1 -0.36 %, alpha=2 -0.39 % vs elastica; b=0.1 converges to -4.5 / -4.1 % (plate physics: -5.0/-1.2/-0.36 % for b=0.1/0.05/0.025). Tests: `test_static_nonlinear_newton_elastica.py` (5), UL numbers unchanged (3.63 / 4.27 / 5.29 %).
- [ ] T5 (S1, S2) Re-register rows/verdicts, update audit doc, run proportionate checks. Commit.

## Log

- L1 (user, verbatim): "haz las dos cosas e investiga el defecto de α=2 hasta resolverlo"
