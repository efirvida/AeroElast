# AGENTS.md

Guidance for coding agents in this repository. The installed Python package is `aeroelast`
(`src/aeroelast/`, `import aeroelast`); the name `fem-shell` survives in the README, the
long-standing CLI names (`fem-shell-fsi`, `fem-shell-monitor`, `fem-shell-bem-fsi`) and the
SLURM scripts. Grep and import `aeroelast`, never `fem_shell`.

## Read this first

- **Bootstrap every command that imports the package.** The extension is built against GCC 14;
  a bare invocation dies at collection with `CXXABI_1.3.15 not found`, and shell state does not
  survive between commands:

  ```bash
  scripts/aeroenv.sh python -c "import aeroelast"     # or: source scripts/aeroenv.sh
  scripts/check.sh quick      # ~13 s smoke gate; exit 1 only for failures outside known-red
  scripts/check.sh full       # the real suite (~25+ min); log kept under logs/check/
  ```

- `check.sh quick` is a **smoke gate**: it never says anything about physics. Escalate by what
  changed (table in §Validation anchors).
- `.pi/skills/fem-shell/SKILL.md` carries the same workflow in detail — bootstrap, canonical
  commands, landmark map, memory protocol, validity rules. Load it instead of re-deriving.
- Work in **bounded sessions**. Resume from `odd/tasks/<feature>.md` plus one narrow
  `mem_search`; do not continue a session that already carries hundreds of turns, because
  context cost grows with the square of session length. `scripts/token_audit.py` measures it.

## Session economy

Measured here: 94 % of the billed tokens are accumulated conversation history, tool payload is
0.15 %, and cost grows with the square of the session length. Behaviour that follows from it:

- **One thing per session.** If a request mixes unrelated workstreams, say so and propose the
  split instead of doing both in one context.
- **Propose the close before the context forces it.** `.pi/extensions/session-guard.ts` shows the
  budget in the status bar; when it passes ~200 k, finish the work unit, write the handoff into
  `odd/tasks/<feature>.md` and hand back the resume recipe. `/handoff` does exactly that, and
  costs the user one keystroke.
- **Never re-read** what is already in context or what a subagent already summarised. Spot-check
  at most one detail.
- **Bound every tool output**: `--stat`, `tail`, `grep -c`, `-q`, or redirect to `logs/check/`
  and read the summary line. Never paste a log you can summarise.
- **Price the verification before running it.** `scripts/check.sh quick` is ~13 s; the full suite
  is 25+ min and gets killed on a login node. State which one you are choosing and why.
- **Do not edit `AGENTS.md` or a skill mid-session** unless asked: the prompt is the front of the
  request prefix, and Pi documents a prompt change as a checkpoint that can invalidate the
  cached prefix for the rest of the session.

## Physics discipline

The product of this repository is numbers that survive scrutiny, not code that runs.

1. **Every number needs a metric, an independent reference, and an accepted tolerance.** A
   deflection with no reference is not evidence. If the tolerance is not written down, the task
   is to find or justify it — not to pick one that passes.
2. **Never widen a pass band to make a test green.** A moved anchor is a finding to report.
   Three of the four documented-red tests here are accuracy questions, not bugs.
3. **Convergence, not single values.** Prefer a mesh-refinement trend or a global scalar such as
   the work done by the load (`W = ½ Σ f·u`) over one node's displacement, which moves with node
   count (`docs/origin_main_integration_2026-09-30.md` §B).
4. **Suspect units, signs, and frames before suspecting the physics.** Nearly every "wrong"
   result in this repo turned out to be one of the three; the `tools/diagnose_*` scripts exist
   for exactly those hunts.
5. **State the baseline.** A before/after claim needs the revision, the environment, and the
   same command on both sides. Read `docs/origin_main_integration_2026-09-30.md` before
   measuring: it holds the suite before/after counts and the anchors that moved.
6. **Rebuild before believing a Rust-side number**: `cd crates/aeroelast-py && maturin develop
   --release`. `pip install -e .` does not recompile Rust, and a stale `.so` silently answers.
7. **Watch the degenerate inputs.** New validation work uses the official
   `tests/IEA-15-240-RWT.yaml`. The UTD xlsx blade has a 30–75 % thinner TE reinforcement
   outboard (≈30 % less edgewise stiffness; +40–55 % static tip deflection) and survives only
   as an input-sensitivity record in `convergence_b1_b2`
   (`docs/blade_input_divergence_utd_vs_official.md`).
8. **Report negative and null results.** A reference that disagrees is the deliverable.

## Conventions that bite

| trap | the rule |
|---|---|
| stress units | `element.Cm()` returns the **integrated** `D = C·h` (N/m). For Pa divide by thickness: `C_mat = element.Cm()/h`. |
| geometric stiffness input | `assemble_geometric_k` treats its σ input as **Pa** and multiplies by `h` internally — do not pre-multiply by thickness. |
| composite stress recovery | resolve the **ply**, not the section mean; the MIDDLE tie-break is deliberate (upstream #27). |
| Voigt layout | solid: 6 components `[σxx,σyy,σzz,τxy,τyz,τzx]`; shell: 3 in-plane, but result arrays keep the 6-slot layout with out-of-plane terms zeroed. |
| DOF layout | `ElementFamily` is a **Rust** enum re-exported by `src/aeroelast/elements/__init__.py`. Single source of truth: `MeshAssembler._FAMILY_PROPERTIES` (`core/assembler.py`): `SHELL → (6, [u,v,w,θx,θy,θz])`, `PLANE → (2, [u,v])`. **There is no `SOLID` family** — 3D solid support was removed on 2026-09-30. `dofs_count` prefers the Rust authoritative value; reduced-DOF scatter must use the family's `dofs_per_node`, and mixed meshes infer from `dofs_per_node == 6`. |
| frames and signs | rotating frame with explicit Coriolis; `CoordinateTransforms` (Rodrigues), `InertialForcesCalculator` and the Omega providers live in `fsi/corotational.py`. Sign conventions are recorded in `odd/tasks/rated-twist-sign-convention.md`. |
| centrifugal `K_G` | solved **statically** (`K u = f_cf` on free DOFs), then membrane stresses recovered from `u`; `assemble_geometric_stiffness` **requires `free_dofs`** (clamped-root elimination). The old local `σ = ρω²r·l_char` was 2–3 orders of magnitude too small, so rotating modes looked ≈ parked. |
| solver naming | `LinearDynamicFSIRotorSolver` is an alias of `LinearDynamicFSIRotorCorotationalSolver` (`rotor.py`); the ~66 case YAMLs that still use the alias are correct, not stale. |
| hot-loop boundary | inside a converged preCICE window the loop never re-enters Python: the KSP is factorized once and reused (`solve_with_cached_ksp`, `dynamic_newmark.rs`) because `K_eff` is constant for fixed Δt. Refactorization happens only when `K_G`/`K_SP` change. |

## Where things live

| need | path |
|---|---|
| YAML → solver dispatch | `src/aeroelast/solvers/fsi/runner.py` (keys on `SolverType`) |
| production rotor solver | `src/aeroelast/solvers/fsi/rotor.py` (`:259`) |
| inertial variant (WIP) | `src/aeroelast/solvers/fsi/rotor_inertial.py` |
| stress-stiffened base (`_solve_via_rust`) | `src/aeroelast/solvers/fsi/stress_stiffened_dynamic.py` |
| solver ABC and re-exports | `src/aeroelast/solvers/solver.py`, `solvers/__init__.py` (aliases live here) |
| config dataclasses, `SolverType` | `src/aeroelast/core/config.py` (loads preCICE XML via `auto_complete_from_precice()`) |
| assembler, DOF conventions, `K_G` helpers | `src/aeroelast/core/assembler.py` |
| stress / ply recovery, sectional post | `src/aeroelast/postprocess/` |
| Rust per-sub-iteration FSI loop | `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` |
| Rust Newmark stepper, cached KSP | `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs` |
| PyO3 crate (maturin root) | `crates/aeroelast-py/Cargo.toml` |
| theory: FSI rotor formulation | `docs/formulations/teoria_formulacion_fsi_rotor.md` |
| theory: frames, inertial forces, solvers | `docs/formulations/solvers.md` |
| theory: shell elements, MITC4+/D | `docs/formulations/shell-elements.md`, `mitc4plusd-2025-extract.md` |
| merge state, anchors, open decisions | `docs/origin_main_integration_2026-09-30.md` |
| closures and re-run triggers | `docs/validation_closures.md` |
| consolidated rotor validation report | `docs/validation_corotational_rotor_solver.md` |
| validation matrix (`tools/validation_matrix.py check`) | `odd/tasks/validation-matrix-store.md` |
| feature/ODD documents | `odd/tasks/*.md` (up to 156 kB — read the section, not the file) |

`solver.py` defines the ABC; `StaticLinearSolver`, `StaticNonlinearSolver`, `DynamicNewmarkSolver`
(aliased `LinearDynamicSolver`), and `ModalSolver` build on it; `LinearDynamicFSISolver`
(`fsi/linear_dynamic.py`) is the preCICE base for the two rotor solvers.

### Rotor FSI hot loop

```text
FSIRunner.run() → rotor.solve()
  → _assemble_system_matrices()            # PETSc: K, M, K_G, K_SP
  → _solve_via_rust(..., step_callback=_step_cb)
    → _aeroelast.run_rotor_fsi_solver()    # Rust owns the loop from here
      ├─ per preCICE sub-iteration: read_data → transform forces →
      │  centrifugal/coriolis/euler → Newmark step → write_data
      └─ per converged window: calls back into Python `_step_cb` (rotor.py:2021)
           logging, metrics, stress recovery, checkpoints
           (it re-computes inertial forces for output — known wart)
```

## Validation anchors

Anchors are the physics contracts. They are not touched to make a run pass, and their accepted
deltas live in the documents listed here — do not duplicate numbers into this file.

| anchor | what it closes | where its deltas are recorded |
|---|---|---|
| S-0, S-1, S-3 | omega-zero consistency, sectional, Newmark transient | the `tests/test_iea15mw_s*.py` docstrings |
| S-2 | static prescribed (gravity and LC sets) vs beam references | `docs/origin_main_integration_2026-09-30.md`; official blade closes at +8.1 % |
| S-4 | rotating modal vs OpenFAST v5 MBC3 at 7.56 rpm | inline: 1st flap +0.6 %, 1st edge −4.5 %, 2nd flap −1.5 % |
| S-5, S-6 | one-way 100 s rated run; one-way dynamic (S-6 NaN open) | `docs/validation_closures.md` |
| S-7 | torsion ratio | `docs/validation_closures.md` (superseded 2026-10-07: 0.9396 → 0.9336) |
| V-01/V-02/V-03/V-05 | rotor performance, natural frequencies (2F open), gravity, structural properties | `docs/validation_results_v01_v02_v04.md` |
| box EI, D-Tube, UL elastica | element-level benchmarks | test docstrings + the integration doc |

Minimum evidence by change type: docs/refactor → `check.sh quick`; element, DOF, convention or
BC → quick **plus** the matching reference test; solver formulation → the anchor number with its
reference, never "it runs"; anything under `crates/` → `maturin develop --release`, then a
full-suite measurement before a merge.

## Build, environment, dependencies

- Python 3.12 + PETSc/SLEPc (`petsc4py`, `slepc4py`), `mpi4py`, Gmsh, meshio are required;
  preCICE is optional (guarded import — without it the package imports and the FSI solvers are
  disabled); OpenFOAM only for coupled runs. Top-level imports tolerate a missing PETSc so
  `--preview` / `--validate` work in a light environment.
- Extras: `pip install -e .[dev]` (pytest), `.[mesh]` (Triangle), `.[bem]` (CCBlade, manual
  build — see `pyproject.toml`).
- HPC dependency installs use `Makefile.fedora`, `Makefile.sdumont`, `Makefile.wsl`; never run
  `make all` casually, they compile PETSc/preCICE/OpenFOAM from source for hours.
- Cluster build details (module load, `maturin` flags, stale flat `_aeroelast.so`) are in
  `.github/skills/build-aeroelast/SKILL.md`.
- OpenFAST S-4/S-5 reference decks, the v4.2→v5 conversions and the ROSCO build are recorded in
  Engram (`config/openfast-5-0-0-corriendo-iea-15-mw-deck-v4-2-v5-0-convertido`); the `/tmp`
  paths in older notes do not survive a reboot, so rebuild rather than reuse them.

## Suite state

- The full suite takes ~25 min on a good node (a 2026-10-09 run reached 53 % in 40 min) and
  emits 190 k–270 k warnings; the failure set moves between runs on the same tree. Treat it as
  a measurement, not a smoke test. On a login node it can be killed part-way (two attempts died
  around 17 % on 2026-10-09), so submit it as a SLURM job when the result matters.
- `scripts/known-red-quick.txt` / `known-red-full.txt` hold the documented-red tests per mode,
  stamped with the revision they were derived at. Re-derive with
  `scripts/check.sh <mode> --record`; never hand-edit to make a run green.
- The four documented-red tests: `test_corotational_is_frame_objective_tl_is_not` (deliberate —
  the MITC4 path of `assemble_kt_corotational` now uses upstream's total-Lagrangian tangent, so
  the asserted property no longer differs), the D-Tube and UL-elastica accuracy cases, and
  `test_iea15mw_s2[LC1_gravity-utd]` (the disowned UTD input).
- `tests/validation/benchmarks/test_ko2017_performance.py` has 8 pre-existing failures on coarse
  distributed meshes — not regressions.
- Removed in the 2026-09-30 merge and therefore absent: 3D solid support and its tests
  (`test_vol_mesh.py`, `test_solid_elements.py`, `test_beam_4cases_parity.py`), plus
  `test_blade_mesh.py` and `test_rotor_inertial.py`. CalculiX helpers read
  `tests/support/ccx_io.py`.

## Adjacent workflows

- **SDD**: `.github/agents/`, `.github/skills/` and `openspec/` carry the Spec-Driven
  Development flow; Engram is the active artifact store. Prefer it for substantial changes.
- **OpenFOAM + preCICE cases** live in a `simulations/` workspace outside this tree
  (`precice-config.xml`, `solid/simulation.yaml`, `fluid/`, `runAll.srm`). Pitfalls:
  `.github/instructions/openfoam-fsi.instructions.md`; SLURM lifecycle:
  `.github/skills/run-fsi-simulation/SKILL.md`.
