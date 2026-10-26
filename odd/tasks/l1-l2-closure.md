# L1/L2 closure: re-run and re-audit the element and linear-solver validation levels

Close the validation levels that only need tests re-run and rows re-audited (L1 MITC shell elements, L2 linear
structural solver), before deciding the PhD-use plan. Store: `docs/validation/`, tool `tools/validation_matrix.py`.

## Specs

- S1 (user L1, verbatim): "yo te propongo una cosa vamos a cerrar lo que depende de correr tests, como por ejemplo L1 y L2, despues pesamos en esa propuesta"
- S2 (repository physics discipline, AGENTS.md): every row keeps its independent reference and tolerance; a pass band is never widened; a moved anchor is a finding to report, not a baseline to overwrite silently; a row with no independent reference is reported, not promoted.
- S3 Done means: `tools/validation_matrix.py contract --level L1` and `--level L2` report the level closed, or every entry that keeps it open is listed with its cause (finding, missing reference, or a decision for the user).

## Starting state (2026-10-26, HEAD 617ed3e, `contract`)

- L1 [open]: trusted 127, regression_pin 2, stale 18 (stale by 7dfe78d 7d20c81 / 61b0d52), unaudited 12 (all `material_suite`); over ceiling 0; gaps 0. Rows: `logs/check/contract-L1.txt`.
- L2 [open, blocked by L1]: stale 83; over ceiling 0; single or no arbiter 4; gaps 0 blocking / 1. (`contract --level L2 --rows` takes > 60 s.)

## Tasks

- T1 | S1-S3 | worker (background) | in progress | L1: re-run the stale groups (`regression --group`), refresh verdicts whose residual is unchanged within its tolerance, audit the 12 unaudited `material_suite` rows, report the 2 regression pins.
- T2 | S1-S3 | worker (background) | in progress | L2: same for the 83 stale rows; report the 4 single/no-arbiter rows and the gap.
- T3 | S3 | parent | pending | review, commit, report to the user; then the PhD-use proposal.

## Handoff 1 (2026-10-26, session 9 close)

**State at close.** Worker `mvp8tox3-1-qu6u` (T1+T2) was still running: no change in `docs/validation` yet;
its outputs are in `logs/check/l1l2/` (`contract-L2-start.txt`, `reg-15..18.txt`, `j-a..e.log`) and SLURM jobs
`run.sbatch` 11624500 / 11624504 were running. Its final report may not survive the session close.
**Resume.** Read this file and `odd/tasks/thesis-validation-plan.md` (strict order: this is step 1).
1. `git status docs/validation` and `ls logs/check/l1l2/`: if the worker wrote rows, review them (no `--write`
   over a moved residual, no widened tolerance) and commit; if not, re-run T1/T2 from the row lists.
2. `scripts/aeroenv.sh python tools/validation_matrix.py contract --level L1` (and L2, > 60 s): report closed
   rows and what stays open with its cause. 3. Only then step 2 of the thesis plan.

## Log

- L1 (2026-10-26, session 9) User: "yo te propongo una cosa vamos a cerrar lo que depende de correr tests, como por ejemplo L1 y L2, despues pesamos en esa propuesta"
- L2 (2026-10-26, session 9) User (thesis goal, verbatim): "volviendo a lo que necesito para mi tesis, el problema es que nadie va a creer en mis resultados si no presento una validacion cientifica del mismo que  pueda ser aveptada. Lo que quiero estudiar de forma general es medir como influye la deformacion de las palas cuando el aerogenerador produce energia sobre la se~nal de torque del mismo, incluyendo cuando la pala pasa frente a la torre. Pero de momento para validacion no hay muchos datos por lo que el proceso de validacion tiene que ser particionado por componentes, lo unico que he encontrado para validar el es articulo de Zhou. En otro proyecto estoy haciendo la validacion numerica de OpenFoam con metodos ALM/ASM pero como es numericamente costozo aqui tenemos un BEM simple para validar la parte estructural y el acoplamiento". Evidence: the built-in BEM (`src/aeroelast/solvers/bem/`) has no tower influence model (no file mentions "tower"); Zhou et al. 2025 (Energy 336:138488) is used by `test_blade_rated_twist.py` (tip twist -3.60 deg, rigid 2.53 MN / 16.11 MW), `test_bem_load_frame.py`, V-02.
