# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Identity Quirk

The repo is named `fem-shell` and the README/copilot docs still reference that name, but **the installed Python package is `aeroelast`** (`src/aeroelast/`, `import aeroelast`). The CLI scripts keep both naming conventions:

- `fem-shell-fsi`, `fem-shell-monitor`, `fem-shell-reconstruct-csv`, `fem-shell-bem-fsi` (kept for backwards compatibility)
- `aeroelast` (newer entry point)

When grepping or writing imports, use `aeroelast`, not `fem_shell`.

The `LinearDynamicFSIRotorSolver` was renamed to `LinearDynamicFSIRotorCorotationalSolver`; an alias is kept under the legacy name (`LINEAR_DYNAMIC_FSI_ROTOR` enum value). A second variant `LinearDynamicFSIRotorInertialSolver` is in development (see `docs/rotor_inertial_solver_design.md`).

## Build, Test, Lint

The project is a hybrid **Python + Rust** workspace. Rust crates compile to `_aeroelast` (a PyO3 extension module).

```bash
# Editable Python install (also triggers Rust build via setuptools shim)
pip install -e .

# Optional extras
pip install -e .[mesh]   # adds Triangle for 2D meshing
pip install -e .[dev]    # pytest stack
pip install -e .[bem]    # BEM aero (requires CCBlade — see pyproject.toml for manual build)

# Lint (Ruff, line-length=100, target py38 for ruff but project requires py3.12+)
ruff check

# Full test suite (excludes stale imports and the WIP 3D volumetric-mesh pipeline)
python -m pytest tests/ -q --tb=short \
    --ignore=tests/test_blade_mesh.py \
    --ignore=tests/test_rotor_inertial.py \
    --ignore=tests/test_vol_mesh.py

# Single test
pytest tests/test_rotor_rust_parity.py::test_name -v

# Skip slow / benchmark tests
pytest tests/ -m "not slow and not benchmark"
```

For HPC/system-level installs (PETSc, SLEPc, preCICE, OpenFOAM, OpenFOAM-preCICE adapter) use the platform-specific Makefiles in the repo root: `Makefile.fedora`, `Makefile.sdumont`, `Makefile.wsl`. Never run `make all` casually — these compile dependencies from source and take hours.

## Architecture — Python + Rust Hybrid

The hot path of FSI rotor simulations runs in Rust; Python orchestrates configuration, mesh handling, and per-converged-window callbacks. Understanding this split is essential before editing any solver.

```
src/aeroelast/                       # Python package
  cli/                               # YAML-driven entry points
    run_fsi.py                       # fem-shell-fsi (FSI orchestrator)
    aeroelast.py                     # aeroelast (newer multi-solver dispatcher)
    run_bem_fsi.py                   # BEM-only coupling variant
    fsi_monitor.py                   # Live TUI monitor for running cases
    reconstruct_rotor_csv.py         # Recover rotor_performance.csv from checkpoints
  core/
    config.py                        # FSISimulationConfig + 30+ dataclasses, parses preCICE XML
    mesh/, material.py, bc.py, laminate.py, assembler.py
  elements/                          # ElementFamily enum; element classes live in Rust now
  solvers/
    linear.py, modal.py              # Static + modal (PETSc/SLEPc)
    elasticity/                      # Static linear/nonlinear, dynamic Newmark
    fsi/
      base.py                        # preCICE Adapter wrapper
      linear_dynamic.py              # LinearDynamicFSISolver (base FSI class)
      stress_stiffened_dynamic.py    # Linear dynamic with K_G updates
      rotor.py                       # LinearDynamicFSIRotorCorotationalSolver (rotating frame)
      rotor_inertial.py              # LinearDynamicFSIRotorInertialSolver (alternative formulation, in development)
      corotational.py                # CoordinateTransforms (Rodrigues), InertialForcesCalculator, OmegaProvider hierarchy
      runner.py                      # FSIRunner — dispatches YAML → solver
      time_integration.py, force_clipper.py
    bem/                             # CCBlade-based BEM aero engine
  postprocess/                       # StressRecovery, watchpoint plots

crates/                              # Rust workspace
  aeroelast-core/                    # Element kernels (MITC3/4, HEXA, TETRA, ...)
  aeroelast-mesh/                    # Mesh data structures
  aeroelast-solvers/                 # PETSc-backed assembly + Newmark stepper + FSI rotor loop
    src/petsc/elasticity/dynamic_newmark.rs  # NewmarkDynamicStepper with cached KSP factorization
    src/petsc/fsi/rotor_fsi.rs               # Main per-sub-iteration FSI loop
  aeroelast-py/                      # PyO3 bindings → builds _aeroelast.so
```

### Hybrid call graph for the rotor FSI hot loop

```
Python: FSIRunner.run()
  → LinearDynamicFSIRotorCorotationalSolver.solve()
    → _assemble_system_matrices()              # PETSc: K, M, K_G assembly
    → _solve_via_rust(..., step_callback=_step_cb)
      → _aeroelast.run_rotor_fsi_solver(...)   # Rust takes over
        ├── per preCICE sub-iteration (Rust):
        │     read_data → transform forces → centrifugal/coriolis/euler → Newmark step → write_data
        └── per converged window: invokes Python _step_cb
              ↓
              Python _step_cb (rotor.py:1827)
                logging, metrics, stress recovery, checkpoint handling
                (currently re-computes inertial forces for output — known wart)
```

Key takeaway: **inside a converged window the loop never re-enters Python**. The KSP is factorized once and reused (`dynamic_newmark.rs:492-496`); refactorizations only happen when K_G or K_SP change.

## Solver Hierarchy

```
Solver (ABC)
├── LinearStaticSolver
├── LinearDynamicSolver (Newmark-β trapezoidal)
├── ModalSolver (SLEPc)
└── LinearDynamicFSISolver (preCICE)
    ├── LinearDynamicFSIRotorCorotationalSolver  # rotating frame, K_G + K_SP, explicit Coriolis
    └── LinearDynamicFSIRotorInertialSolver       # in development — inertial frame variant
```

The corotational solver is the production rotor solver. See `docs/teoria_formulacion_fsi_rotor.md` for the full theoretical derivation, `docs/validez_teorica_fsi_rotor_corotational.md` for the validity assessment, and `docs/mejoras_rendimiento_fsi_rotor_corotational.md` for the performance audit.

## Configuration Conventions

- Top-level config is `FSISimulationConfig` (`core/config.py`), loaded via `FSISimulationConfig.from_yaml(path)`.
- `config.auto_complete_from_precice()` reads `total_time`, `time_step`, and watchpoint files from the preCICE XML — don't duplicate them in the YAML.
- Solver dispatch in `runner.py` keys on `SolverType` enum (`config.py`). Both `LINEAR_DYNAMIC_FSI_ROTOR_COROTATIONAL` and the legacy alias `LINEAR_DYNAMIC_FSI_ROTOR` route to the corotational solver.
- `RotorConfig.to_dict()` serializes to the dict format the rotor solver constructor expects.
- `moment_of_inertia: "auto"` triggers automatic I computation from mesh (lumped mass × r_perp²).
- `send_omega_to_precice: true` causes the solver to write `AngularVelocity` on a `GlobalSolidMesh` (single vertex at rotation center). This is handled internally — do NOT add it to `coupling.write_data`.

### Canonical blade input (since 2026-09-08)

The official IEA 15 MW blade in WindIO format is `tests/IEA-15-240-RWT.yaml`
(`BladeMesh(yaml_file=...)` / `mesh.generator.params.yaml_file` in campaign
YAMLs). The previous input `NuMAD_utd_iea15mw.xlsx` (UTD re-modeling,
AIAA 2023-2093) has a TE reinforcement 30–75% thinner than the official
layup, which reduced edgewise stiffness by ~30% outboard (S-2 static tip
deflection +40–55% vs beam references; official blade closes at +8.1%).
Full evidence: `docs/blade_input_divergence_utd_vs_official.md`. The
convergence_b1_b2 campaign runs intentionally remain on the UTD blade as an
input-sensitivity record. Do not re-introduce `excel_file:
NuMAD_utd_iea15mw.xlsx` in new validation work.

### Centrifugal geometric stiffness (since 2026-09-09)

`MeshAssembler.assemble_geometric_stiffness(omega=...)` solves the centrifugal
pre-stress problem **statically** (`K u = f_cf` on the free DOFs, then
`assemble_geometric_k_from_disp` recovers the membrane stresses from `u`).
The old approximation `σ = ρω²r·l_char` (local, per-element) was 2–3 orders of
magnitude too small — rotating modes were ≈ parked. The method now **requires
`free_dofs`** (clamped-root elimination); callers in `rotor.py` and the
diagnose tools pass `bc_manager.free_dofs`. `assemble_geometric_k` interprets
its σ input as **Pa** and multiplies by `h` internally (do not pre-multiply by
thickness). S-4 (rotating modal vs OpenFAST v5 MBC3, `tests/test_iea15mw_s4_rotating_modal.py`)
closes: 1st flap +0.6%, 1st edge −4.5%, 2nd flap −1.5% at 7.56 rpm.

### OpenFAST environment (for S-4/S-5 targets)

OpenFAST v5.0.0 (conda-forge) runs the IEA 15 MW monopile deck; the repo's
v4.2 inputs were converted to the v5 schema (deltas documented in Engram,
`config/openfast-5-0-0-corriendo-iea-15-mw-deck-v4-2-v5-0-convertido`).
Key setup: micromamba env `openfast`, ROSCO `libdiscon.so` compiled with
gfortran in `/tmp/opencode/ROSCO/rosco/controller/build`, run decks in
`/tmp/opencode/ofrun*` (S-4 = linearization at 7.56 rpm with only blade
DOFs; S-5 = 100 s rated run with 313 OutList channels). MBC3 post-processing
via `openfast_toolbox` (`pip install git+https://github.com/OpenFAST/openfast_toolbox`,
`openfast_toolbox.linearization.mbc.fx_mbc3`).

### Omega Provider modes

| `moment_of_inertia` | `omega_ramp_time` | Provider class       |
|---------------------|-------------------|----------------------|
| `null`              | `0`               | ConstantOmega        |
| `null`              | `> 0`             | RampedOmega          |
| `"auto"` / float    | `0`               | ComputedOmega        |
| `"auto"` / float    | `> 0`             | RampedComputedOmega  |

Other types: `TableOmega` (tabulated time series), `FunctionOmega` (user callable).

## Element & DOF Conventions

| Family    | DOFs/node | Layout                       | Elements                                              |
|-----------|-----------|------------------------------|-------------------------------------------------------|
| `PLANE`   | 2         | `[u, v]`                     | QUAD4, QUAD8                                          |
| `SHELL`   | 6         | `[u, v, w, θx, θy, θz]`      | MITC3, MITC4, MITC3Composite, MITC4Composite          |
| `SOLID`   | 3         | `[u, v, w]`                  | HEXA8/20/27, TETRA4/10, WEDGE6/15, PYRAMID5/13        |

Mixed-element meshes use the **max stride (6)** for DOF indexing — be aware when scattering reduced-DOF vectors back to nodal fields.

### Voigt notation

- Solid: 6 components `[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_zx]`
- Shell: 3 components `[σ_xx, σ_yy, τ_xy]` (plane stress); the result arrays still use the 6-component layout with the out-of-plane terms zeroed.

### Stress recovery gotcha

`element.Cm()` returns the **integrated** membrane stiffness `D = C·h` (units N/m). To get actual stress (Pa), divide by thickness: `C_mat = element.Cm() / h`. This trips up most people the first time.

## External dependencies

| Package           | Purpose                                  | Required?                  |
|-------------------|------------------------------------------|----------------------------|
| PETSc (`petsc4py`)| Sparse assembly, KSP solvers             | Yes                        |
| SLEPc (`slepc4py`)| Eigenvalue problems (modal)              | For modal analysis         |
| preCICE (`precice`)| FSI coupling                            | Optional (guarded import)  |
| MPI (`mpi4py`)    | Parallel computing                       | Yes                        |
| Gmsh (`gmsh`)     | Mesh generation                          | Yes                        |
| meshio            | VTU/PVD checkpoint export                | Yes                        |
| OpenFOAM          | Fluid solver (CFD side)                  | For full FSI runs          |

If preCICE is missing, the package still imports — FSI solvers are simply disabled. Top-level imports also tolerate missing PETSc so `--preview` and `--validate` work in lightweight environments.

## Known test-suite quirks

- `tests/test_blade_mesh.py` and `tests/test_rotor_inertial.py` have stale imports — exclude them by default (see test command above).
- `tests/test_vol_mesh.py` no longer exists: the 2026-09-30 merge accepted upstream's removal of 3D solid support (Rust `7b295f3`, Python `9230ea2`), which also deleted `tests/test_solid_elements.py` and `tests/test_beam_4cases_parity.py`. The `--ignore=tests/test_vol_mesh.py` in the test command is now a harmless no-op. The surviving CalculiX checks read `tests/_ccx_io.py`.
- `tests/test_ko2017_performance.py` has 8 pre-existing failures with tight tolerances on coarse distributed meshes — these are NOT regressions.
- Some benchmarks under `tests/` are intentionally heavy and unsuitable for quick smoke tests.

## Element line and validation state (since 2026-09-30)

The element kernels and the CalculiX writer now come from `origin/main` (the reviewed
MITC4+/D element), not from the older local line: `crates/aeroelast-core/src/assembly/`
was taken from upstream in the merge, and only the K_G-static helpers
(`centrifugal_load`, `assemble_geometric_k_from_disp`, `update_node_coordinates`,
`assemble_kt_corotational`) and the `body_load_rho` fix are ours on top.

Before measuring anything, read `docs/origin_main_integration_2026-09-30.md`: it holds the
before/after suite counts, the 28 CalculiX tests blocked by one degenerate blade-mesh
element (`*ERROR in e_c3d: nonpositive jacobian`, element 2790), the anchors that moved
(V-02 2F, S-4 rotating 1F, S-7 torsion, S-6 NaN, box EI, D-Tube, UL elastica), and the
open decisions. `docs/validation_closures.md` carries the same re-run trigger at the top.

## SDD workflow

The `.github/agents/`, `.github/skills/`, and `openspec/` directories indicate this project uses Spec-Driven Development (SDD). When proposing substantial changes, prefer the SDD flow: `/sdd-new <change>` → propose → spec → design → tasks → apply → verify → archive. Engram is the active artifact store.

## OpenFOAM + preCICE simulations

Sample cases live under a `simulations/` workspace (not in this repo's tree directly). Each case has:
- `precice-config.xml` — coupling scheme, data mappings, convergence
- `solid/simulation.yaml` — YAML config consumed by `fem-shell-fsi`
- `fluid/` — OpenFOAM case (`system/`, `constant/`, `0/`)
- `runAll.srm` — SLURM submission script

See `.github/instructions/openfoam-fsi.instructions.md` for OpenFOAM coupling pitfalls and `.github/skills/run-fsi-simulation/SKILL.md` for the SLURM job lifecycle.
