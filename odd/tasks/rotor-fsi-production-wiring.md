# Rotor FSI test wiring to the production solver

Date: 2026-10-06 · Branch: `integrate/origin-main-2026-09-30`

Trigger: the coupling regression (17-30 sub-iterations per window against 7-12 in the campaign
era) went undetected because **no FSI simulation runs in the test suite** — every one is a manual
SLURM script (`run_step1b_smoke.srm`, `run_step1c_mp_off.srm`, `run_zeta_sweep.srm`,
`run_step1_revalidate.srm`, `submit_exp_feedback_off.sh`).

## Scope

Target solver: `LinearDynamicFSIRotorCorotationalSolver` (production, `rotor.py:232`).
Out of scope by maintainer decision: `LinearDynamicFSIRotorInertialSolver` (WIP, `rotor_inertial.py`).

## Findings that shape the plan (read-only exploration, 2026-10-06)

- **No in-process bidirectional coupling exists.** Every FSI coupling is a preCICE participant over
  sockets: `run_bem_fsi.py:378` (`BEMFSIParticipant.run()`), `run_aero_fsi.py:257`
  (`AeroFSIParticipant.run()`). `BEMStandaloneSolver.solve()` (`bem/standalone.py:51`) is a
  single-pass projection with no loop. A pytest FSI test must therefore launch **two processes**
  and is slow by construction.
- **Cheapest bidirectional case on disk**: `tests/IEA15MW/convergence_b1_b2/{solid,fluid}_h_coarse.yaml`
  (`element_size: 0.5`) against the campaigns' 0.25 m baseline.
- **Start-up cost**: a 0.25 m case spends ~10-20 min (mesh ~50 s for the fluid, ~193.5k free DOF
  assembly, KSP/MUMPS factorization, preCICE init) before window 1. 0.5 m is roughly 10x fewer DOF.
- **66 YAML files** declare `type: LinearDynamicFSIRotor`, the legacy alias. `rotor.py:2536` defines
  `LinearDynamicFSIRotorSolver = LinearDynamicFSIRotorCorotationalSolver` — the same object — so
  they already dispatch to production via `runner.py:1326-1328`. Renaming all 66 is churn with no
  behaviour change.
- **All 18 scipy "replica" solves in tests are DELIBERATE**: each compares against an external
  reference (CCX, OpenFAST, a beam, a benchmark, a parity check) and needs to control the
  formulation. Replacing them would destroy what they validate. The gap is the missing production
  coverage, not the presence of replicas.
- **The replicas are natural arbiters of the production solver (maintainer's insight, verified).**
  All 23 test files that call `spsolve` build `K`/`M` with the PRODUCTION Rust assembler
  (`PyMeshAssembler` / `MeshAssembler` / `blade_validation._to_rust_mesh`; 23/23 with strict
  evidence) and NONE of them ever invokes a production solver class. Same problem, same assembly,
  different solve: so the replica and production can be compared directly, which makes the replica
  an independent arbiter rather than a piece of test-only code. Today no test asserts anything
  about the production solve path.
- Production step call sites: `_aeroelast.run_rotor_fsi_solver` (`rotor.py:2377`) and
  `LinearDynamicFSIRotorCorotationalSolver.solve()` via `FSIRunner.run()` (`runner.py:126`) both
  require preCICE; `LinearStaticSolver.solve()` (`elasticity/static_linear.py:27`) does not.

## Tasks

1. [ ] **WU-1 Wiring guard (fast, no mesh, no preCICE).** Assert the alias IS the corotational
   class, and that every FSI case YAML on disk declares either the alias or the explicit
   corotational name — never the inertial variant, with the known inertial cases whitelisted.
   Fails if a case silently switches to the half-built solver.
2. [ ] **WU-2 Measure the cheapest viable bidirectional case** (0.5 m, 3-5 windows): start-up and
   per-window wall time. Prerequisite for WU-3's shape.
3. [ ] **WU-3 Arbiter parity: the production solver must reproduce the replica.** Keep every
   replica test as it is (it validates against its external reference) and ADD a production-side
   assertion on the same mesh, loads and BCs. Staged to protect the review:
   - **WU-3a** Feasibility prototype on the smallest case (`test_box_torsion_bending_benchmark.py`
     or `test_thin_walled_tube_torsion.py`): drive the in-process production static/dynamic solver
     on the existing fixture and compare. Resolves the `model_config` plumbing, since production
     solvers take `(mesh, model_config)`.
   - **WU-3b** One shared helper (e.g. `tests/_production_arbiter.py`) so each test adds ~3 lines,
     then extend in small batches. Acceptance tolerance is a maintainer decision: state whether the
     arbiter uses the same tolerance as the external reference or a separate one.
   - **WU-3c** The coupled path has no in-process equivalent: a pytest FSI test, marked slow, that
     launches the two participants, asserts both finalize, every window converges (`conv = true`)
     and the de-loading signature holds (flexible thrust below the rigid anchor). This is Gate 1's
     automatable form. If WU-2's measured runtime makes it impractical in the suite, split it: a
     scripted runner target plus a fast test over the stored artifacts.
   Scope note: the production *rotor* solver requires preCICE, so only the non-FSI replicas
   (~20 files) can be compared in-process; the rotor-FSI path is covered by WU-3c.
4. [ ] **WU-4 Campaign cases only**: rename the alias to the explicit corotational name in the
   cases that will be re-run (yaw 5, convergence_b1_b2 4, ch6 corotational 3). Leave the other
   ~54 files untouched to keep the review focused.
5. [ ] **WU-5 Gate 1 update**: the FSI gate becomes the WU-3 test instead of manual job `11609164`,
   and the store/docs record it.

## Constraints

- Do not modify the 9 inertial case YAMLs.
- Do not replace the deliberate scipy replicas: add production coverage where it is missing.
- Every work unit closes with a work-unit commit; `tools/validation_matrix.py check` must stay at
  0 errors.
