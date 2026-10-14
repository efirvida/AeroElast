# Feature: UL elastica rejected rows - register converged evidence, resolve the alpha=2 defect

Worktree `fem-shell-worktrees/rotor-frame`, branch `fix/bem-participant-rotor-frame`. Follows PRE-2 of `odd/tasks/validation-l3-b3-audit.md` (3 L1 rows `rejected`: G33 UL alpha=1 and 2, G22). Diagnosis in Engram `fem-shell/validation/ul-elastica-rejected-rows-diagnosis` (obs 544). No P1 execution, no public GitHub mutation, no push.

## Specs

- S1: "haz las dos cosas" - (a) register the converged-mesh evidence for G33 alpha=1 and G22 and correct the stale docstrings, (b) investigate the alpha=2 defect.
- S2: "investiga el defecto de α=2 hasta resolverlo" - find the cause of the mesh- and step-converged +5.29 % (w) / +13 % (u) alpha=2 gap and fix it, or conclude with evidence that it is not a defect.
- S3: No pass band is widened to make a row green (AGENTS.md physics discipline); a moved anchor is reported.
- S4: Mandatory environment: `PYTHONPATH=$SCRATCH/validenv:$PWD scripts/aeroenv.sh python ...`; verify the import first. SLURM only `-p sequana_cpu --ntasks=...`.

- S5: "el objetivo de los tests es corregir el codigo de produccion. entonces claro que debes corregir el kernel del solver que valla a dar mejores resultados aun si es mas trabajo."
- S6: "arreglar el StaticNonlinearSolver tal vez no corrija el solver FSI linear o el de rotor. pero vamos paso por pase porque es muy probable que ese defecto tambien este alla. entonces corrije el/los solver/s y continuamos"

## Tasks

- [x] T1 (S1, S3) G22 nonlinear solve moved to 64 x 4 (1.3863 % vs 5 %, trusted); G33 alpha=1/2 stay rejected with the corrected cause, bounds untouched; docstrings and group notes corrected (2462e1f, 71eebaf). Not registered as converged-mesh G33 rows on purpose: the alpha=1 agreement with the elastica at nx>=64 is compensation (plate -4.3 % vs scheme +4.4 %), so a green row would be misleading.
- [x] T2 (S2) CalculiX 2.20 S8R NLGEOM arbiter on the b=0.1 plate (`test_strip_elastica_calculix.py`, `nlgeom=` in `ccx_surface.py`): w/L alpha=1 0.28811/0.28880 (nu=32/64), alpha=2 0.47329/0.47424; u/L -0.05154/-0.05175, -0.14712/-0.14761; insensitive to increments (10..40). Newton agrees to 0.16-0.26 % in w and 0.56-0.65 % in u; the default UL sits +5.6 % (alpha=2) from the elastica and ~9 points from CalculiX. Paired bounds 1 % / 2 % (both codes still drift 0.1-0.3 % per doubling).
- [x] T3 (S2) Cause established with a 2-D corotational beam (`$SCRATCH/studies/beam_scheme.py`): tangent flow (K_mat on current geometry, no K_sigma, no equilibrium) converges AWAY from the elastica (alpha=1 w/L 0.31576 vs 0.30172; alpha=2 0.54916 vs 0.49346, u/L -0.2054 vs -0.1606) while Newton on the consistent tangent gives 0.30173 / 0.49349. Commit: see T4.
- [x] T4 (S2, S3) Opt-in `scheme: newton` (total-Lagrangian Newton on fint/kt, planar-valid), UL stays the default (blade twist is multi-axis). Narrow strip b=0.025 nx=64: alpha=1 -0.36 %, alpha=2 -0.39 % vs elastica; b=0.1 converges to -4.5 / -4.1 % (plate physics: -5.0/-1.2/-0.36 % for b=0.1/0.05/0.025). Tests: `test_static_nonlinear_newton_elastica.py` (5), UL numbers unchanged (3.63 / 4.27 / 5.29 %).
- [x] T5 (S1, S2) Groups 52 (newton vs elastica, narrow strip, 2 rows) and 53 (newton vs CalculiX NLGEOM, 6 rows) registered, residual patterns, 8 stale + 8 new verdicts. L1: trusted 129, regression_pin 2, rejected 2 (was 120/2/3), single-arbiter 0, over ceiling 0. Verification: see Log.

- [ ] T6 (S5) Probe the kernel for general (multi-axis) rotations: rigid-rotation objectivity of `assemble_fint` about oblique axes, and `assemble_kt` against finite differences of `assemble_fint` at random multi-axis states (the 2026-09-18 note says the tangent is inconsistent at large rotations).
- [ ] T7 (S5) Fix the kernel (Rust element/assembler) so TL Newton is objective and consistent for general rotations; rebuild with maturin; element-level regression tests.
- [ ] T8 (S5) Make `newton` robust (drilling-DOF spike, step control), validate multi-axis against CalculiX and against the UL on a small-rotation blade case, make it the default, keep `ul` selectable.
- [ ] T9 (S6) Assess the FSI linear and rotor solvers for the same defect (stress-stiffening K_G, K_SP, rotating-frame tangent); fix what is found, one solver at a time, with an independent reference.
- [ ] T10 (S3, S5, S6) Re-register rows and verdicts (G33 rows), update the audit doc, proportionate checks. Commit.

## Log

- L1 (user, verbatim): "haz las dos cosas e investiga el defecto de α=2 hasta resolverlo"
- L2 (agent, 2026-10-14): CORRECTION to the first diagnosis. The plate effect is not a constant -1.1 %: CalculiX reads the b=0.1 strip -3.3 / -4.3 / -3.9 % from the elastica at alpha 4/7, 1, 2. The UL scheme is therefore +2.0 / +4.4 / +9.6 % above CalculiX, and its agreement with the elastica at alpha <= 1 is compensation.
- L3 (agent): adding `static_nonlinear.py` to L1 `depends_on` staled 125 verdicts (the dependency is per level, not per group); reverted. Consequence: a future change to the solver does not stale the G22/G33/G52/G53 verdicts automatically.
- L4 (agent): open decision for the maintainer - make `scheme: newton` the default for planar problems, or add the initial-stress stiffness to the UL loop. The two G33 rows stay rejected until then. `scheme` is already reachable from YAML (`solver:` params; runner.py:1305).
- L5 (agent): a login-node process that builds several solver instances in sequence got SIGKILLed (memory accumulation across instances, not scheme-specific); each case alone peaks at 0.34 GB.
- L6 (user, verbatim): "OK, entonces asumo que hiciste un solver independiente para diagnosticar el que tenemos en produccion y resulto en que tien errores, y eso esta bien, pero el objetivo de los tests es corregir el codigo de produccion. entonces claro que debes corregir el kernel del solver que valla a dar mejores resultados aun si es mas trabajo.  Ten en cuenta que arreglar el StaticNonlinearSolver tal vez no corrija el solver FSI linear o el de rotor. pero vamos paso por pase porque es muy probable que ese defecto tambien este alla. entonces corrije el/los solver/s y continuamos"
