# AeroElast

AeroElast is a high-performance finite element toolkit for structural and
aeroelastic simulation of shell and plane structures, with a focus on
wind turbine blade FSI (fluid–structure interaction).

The project combines a **Rust computation core** with **Python ergonomics**:

- `crates/aeroelast-core` — pure Rust FEM kernels (MITC3/MITC4 shell elements,
  materials, assembly, quadrature), no PETSc or Python dependencies
- `crates/aeroelast-mesh` — Rust mesh data types (Node, Element, MeshModel)
- `crates/aeroelast-solvers` — PETSc/SLEPc-backed assemblers and modal/static/
  dynamic/FSI solvers exposed through FFI
- `crates/aeroelast-py` — PyO3 bindings that call the Rust kernels from Python
- `src/aeroelast` — the Python package: YAML-driven CLI, mesh generators,
  preCICE coupling, post-processing, and monitoring tools

AeroElast is designed for research and engineering evaluation. Interfaces and
workflows are evolving rather than frozen.

## Current Scope

- finite element support for plane and shell elements
- isotropic, orthotropic, and laminated composite material models
- static, dynamic, modal, and FSI structural solvers
- rotor-oriented FSI workflows with angular-velocity feedback to CFD
- BEM-based aerodynamic coupling for blade performance studies
- mesh generation, mesh import/export, and node-set selection tools
- checkpointing, restart, and post-processing helpers
- command-line tools for running, monitoring, and reconstructing simulations

## Project Status

This repository is a research and engineering prototype. It is a good fit for:

- evaluating FEM formulations and solver behavior
- building custom structural or FSI workflows
- running rotor and blade-oriented studies
- testing integration patterns around OpenFOAM, preCICE, and PETSc

It is not yet presented as a finished end-user product with a frozen API,
turnkey installers, or broad industrial validation.

## Main Capabilities

### Structural analysis

- `LinearStaticSolver` for linear static problems
- `LinearDynamicSolver` for transient structural dynamics (Newmark-β)
- `ModalSolver` for modal analysis (SLEPc)
- geometric (stress-stiffened) nonlinearity and co-rotational formulations

### FSI workflows

- `LinearDynamicFSISolver` for structural coupling through preCICE
- `LinearDynamicFSIRotorSolver` for rotating systems with inertial effects,
  torque-driven omega updates, and CFD feedback
- BEM-FSI coupling for aerodynamic blade analysis
- automatic checkpoint/restart support for long-running coupled simulations

### Rust computation core

- MITC3 and MITC4 shell element kernels
- composite laminate materials with failure criteria
- zero-allocation hot paths in the dynamic Newmark stepper
- direct LU solvers for robust transient dynamics
- PyO3 bindings exposed as the `aeroelast` Python extension

### Mesh and model utilities

- built-in mesh generators such as `SquareShapeMesh`, `BoxSurfaceMesh`,
  `MultiFlapMesh`, `BladeMesh`, and `RotorMesh`; the blade mesher is gmsh
  (`GmshBladeMesh`; the YAML key `BladeMesh` is kept as an alias) and reads the
  windIO v2 blade definition (`tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml`).
  The reference mesh is chord-scaled N20, aspect ratio 1; element-size meshes
  use 0.15-0.25 m. windIO v1 yaml and Excel (xlsx) blade input are rejected.
- turbine component meshers (`TowerMesh`, `NacelleMesh`) and a
  full `TurbineMesh` assembly of blades + nacelle + tower
- mesh import/export utilities for common engineering formats
- geometric node-set creation from coordinate, box, distance, and direction
  criteria

### Monitoring and recovery tools

- `aeroelast-monitor` for live TUI-based monitoring of FSI runs
- `aeroelast-reconstruct-csv` to recover rotor performance histories from
  checkpoint data

## Installation

The Python package targets Python 3.12+. The Rust extension (`_aeroelast`) is
required and is built automatically during `pip install` via maturin.

Recommended dev install using a dedicated conda environment with the native
solver dependencies (PETSc/SLEPc and their Python bindings, MPI, preCICE).
`hdf5` is pinned to 1.14 because the Rust HDF5 bindings do not support the
HDF5 2.x series.

```bash
conda create -n aeroelast-dev -c conda-forge python=3.12 petsc slepc petsc4py slepc4py openmpi precice pyprecice hdf5=1.14
conda activate aeroelast-dev
cd AeroElast
export HDF5_DIR="$CONDA_PREFIX"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig"

# On SDumont instead: `module load glu gcc/14.2.0_sequana` and see Makefile.sdumont.
pip install -e .
```

`pip install -e .` builds the Rust extension (`_aeroelast`) with maturin and
registers the available CLI commands. The two `export` lines point the Rust
build scripts at the conda-provided HDF5, PETSc/SLEPc and preCICE. Without the
compiled extension the package cannot import, so the conda env (or an
equivalent environment providing PETSc/SLEPc/preCICE) is mandatory.

On SDumont, load GCC before any `aeroelast`/`fem-shell-*` execution and keep
its `libstdc++` first in `LD_LIBRARY_PATH` to satisfy `GLIBCXX_3.4.32+`
runtime requirements.

If you need Triangle-based meshing helpers, install the optional extra:

```bash
pip install -e .[mesh]
```

For BEM-coupled aerodynamic analysis (CCBlade + NeuralFoil):

```bash
pip install -e .[bem]
```

> Note: CCBlade ships a Fortran extension that pip cannot build reliably on its
> own. See the `bem` extra comments in `pyproject.toml` for the manual
> two-step build, or use `make -f Makefile.sdumont ccblade` on SDumont.

The SHARPy and VLM aerodynamic backends are legacy, outside the validation
path (BEM is the validated aerodynamic model), and their removal is planned
(`odd/tasks/aero-vlm-sharpy-removal.md`); the code and the `sharpy` extra still
exist. Install the extra only if you need it:

```bash
pip install -e .[sharpy]
```

### Optional external dependencies

Some workflows require software that is not bundled with this repository:

- preCICE for FSI coupling
- OpenFOAM for CFD-side coupled simulations
- PETSc/SLEPc (C libraries) for the native solver crates, plus `petsc4py`
  and `slepc4py` for the Python solver layer and test suite
- an MPI runtime for distributed runs where applicable

If preCICE is not available in the environment, the Rust extension cannot be
linked (the `fsi` feature is on by default) and the build will fail.

## Quick Start

Generate a simulation template:

```bash
aeroelast-fsi --template > simulation.yaml
```

Validate the configuration:

```bash
aeroelast-fsi simulation.yaml --validate
```

Preview the resolved configuration without running the solver:

```bash
aeroelast-fsi simulation.yaml --preview
```

Run the simulation:

```bash
aeroelast-fsi simulation.yaml
```

Monitor a running case:

```bash
aeroelast-monitor /path/to/workdir
```

Reconstruct rotor performance CSV data from checkpoints:

```bash
aeroelast-reconstruct-csv results/
```

Generate meshes with the `mesh` subcommand (see
[docs/mesh-cli.md](docs/mesh-cli.md)):

```bash
# Blade shell mesh (gmsh reference mesh N20, AR 1), windIO v2 input
aeroelast mesh blade tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml --out blade.stl
# --no-webs gives the CFD outer surface (the webbed mesh with the web elements and orphan nodes dropped)

# Tapered tower from the WindIO definition (`nacelle` is a separate subcommand)
aeroelast mesh tower tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml --out tower.obj

# Whole turbine, one file per component
aeroelast mesh turbine tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml --out-dir meshes/ --format stl
```

## CLI Commands

The package provides these command-line entry points:

- `aeroelast` - main CLI entry point (`aeroelast mesh ...` generates meshes)
- `aeroelast-fsi` - run or inspect YAML-defined simulations
- `aeroelast-bem-fsi` - run BEM-coupled FSI simulations
- `aeroelast-monitor` - monitor rotor/FSI runs from CSV output
- `aeroelast-reconstruct-csv` - rebuild missing rotor performance histories
- `fem-shell-fsi`, `fem-shell-aero-fsi`, `fem-shell-aero-report-info`,
  `fem-shell-bem-fsi`, `fem-shell-monitor`, `fem-shell-reconstruct-csv` -
  deprecated aliases, kept for backwards compatibility

Detailed CLI documentation is available in
[docs/cli-reference.md](docs/cli-reference.md) for the simulation runner and
[docs/mesh-cli.md](docs/mesh-cli.md) for mesh generation.

Generic Aero-FSI report consumers should also read
[docs/aero-fsi-report-schema.md](docs/aero-fsi-report-schema.md), which explains
the generated CSV files and the accompanying schema manifest.

A historical IEA-15 MW VLM-versus-BEM record (not part of the current
validation order) is kept in
[docs/vlm-bem-iea15mw-validation.md](docs/vlm-bem-iea15mw-validation.md).

## Repository Layout

```text
AeroElast/
├── crates/                # Rust workspace
│   ├── aeroelast-core/    # Pure FEM kernels (elements, materials, assembly)
│   ├── aeroelast-mesh/    # Mesh data types
│   ├── aeroelast-solvers/ # PETSc/SLEPc-backed solvers
│   └── aeroelast-py/      # PyO3 bindings
├── docs/                  # User-facing documentation
├── scripts/               # aeroenv.sh bootstrap, check.sh test gate
├── tools/                 # Validation-store and diagnostic tools
├── src/aeroelast/         # Python package source
│   ├── cli/               # CLI entry points
│   ├── core/              # Mesh, BCs, materials, config
│   ├── constitutive/      # Failure criteria and constitutive helpers
│   ├── elements/          # Element implementations
│   ├── models/            # Blade model (windIO v2 reader, vendored pyNuMAD objects)
│   ├── postprocess/       # Output and analysis utilities
│   └── solvers/           # Static, dynamic, modal, FSI, BEM and aero solvers
└── tests/                 # Unit, regression and validation tests (tests/validation/)
```

## Tests and Validation

- Automated tests live under `tests/`; the benchmark and reference cases are
  under `tests/validation/` (e.g. `tests/validation/benchmarks/`)
- Every command that imports the package needs the environment bootstrap
  (the extension is built against GCC 14):

```bash
scripts/aeroenv.sh python -c "import aeroelast"
scripts/check.sh quick   # ~13 s smoke gate, says nothing about physics
scripts/check.sh full    # the full suite (long; submit it as a SLURM job)
```

Some benchmark cases are intentionally heavier and are not suitable for a
quick local smoke test. A few tests are documented as known-red (listed in
`scripts/known-red-quick.txt` and `scripts/known-red-full.txt`); several of
them are element-accuracy questions rather than solver bugs. The validation
rules are in [docs/validation-policy.md](docs/validation-policy.md).

## Documentation

- [docs/cli-reference.md](docs/cli-reference.md) — CLI reference for the
  simulation runner
- [docs/mesh-cli.md](docs/mesh-cli.md) — mesh generation CLI
  (`aeroelast mesh blade|rotor|nacelle|tower|turbine`)
- [docs/formulations/](docs/formulations/) — theory documents: the shell element,
  the material models, the solvers, the two paper extracts, and the FSI rotor
  formulation (Spanish)
- [docs/validation-policy.md](docs/validation-policy.md) — validation policy
- [docs/windio2-input.md](docs/windio2-input.md) — windIO v2 blade input
- [docs/blade_mesh_size_rule.md](docs/blade_mesh_size_rule.md) — blade mesh
  size rule and its evidence

Generic Aero-FSI runs now emit an `aero_report_schema.json` manifest alongside
their CSV outputs. It records the report files, schema versions, and stable
field ordering so downstream tools do not have to infer column meaning from the
first emitted row.

See [docs/aero-fsi-report-schema.md](docs/aero-fsi-report-schema.md) for the
quick path and field-level guidance.

The `StressRecovery` class (`src/aeroelast/postprocess/stress_recovery.py`)
computes element-centroidal and node-averaged stress and strain fields from the
displacement solution.


## Third-Party Code and References

AeroElast builds on external codes and scientific literature. As a research
project, we give credit where credit is due and respect every license
involved.

### Embedded / vendored code

- **pyNuMAD** (Sandia National Laboratories) — a modified copy of its blade
  definition objects, I/O and gmsh skin-meshing helpers is vendored
  (attribution kept; its own mesher was retired in favour of `GmshBladeMesh`) under `src/aeroelast/models/blade/numad/`, under the
  **BSD 3-Clause License**. See the [LICENSE](src/aeroelast/models/blade/numad/LICENSE)
  and [NOTICE](src/aeroelast/models/blade/numad/NOTICE) files shipped with it.
  Upstream: <https://github.com/sandialabs/pyNuMAD>

### Coupling and solver ecosystem

| Code | Purpose | License | Upstream |
|------|---------|---------|----------|
| preCICE | Partitioned FSI coupling | LGPL-3.0 | <https://precice.org> |
| OpenFOAM | CFD-side coupled simulations | GPL-3.0 | <https://openfoam.org> |
| PETSc | Sparse linear algebra, KSP solvers | BSD-2-Clause | <https://petsc.org> |
| SLEPc | Eigenvalue problems (modal analysis) | BSD-2-Clause | <https://slepc.upv.es> |
| CCBlade | Blade element momentum (BEM) aerodynamics | Apache-2.0 | <https://github.com/WISDEM/CCBlade> |
| NeuralFoil | Neural-network airfoil aerodynamics | Apache-2.0 | <https://github.com/peterdsharples/NeuralFoil> |

### Core Python / Rust dependencies

The Python package depends on NumPy, SciPy, Matplotlib, meshio, Gmsh,
Polars, PyVista, Shapely, Trimesh, NetworkX, Rich, Textual, MPI4Py, HDF5,
and related libraries; the Rust crates depend on nalgebra, rayon, PyO3, and
numpy. Each library is distributed under its own license, which is reproduced
in the corresponding package metadata.

### Scientific references

The canonical bibliography for this repository is
[docs/validation/references.yaml](docs/validation/references.yaml). It lists the shell element
formulations, the constitutive and failure models, the time-integration and
finite element references, the partitioned FSI coupling and BEM aerodynamics
references, and the third-party software this project builds on, with full
author lists and a DOI wherever one could be verified.

The formulation implemented in this repository is documented in
[docs/formulations/teoria_formulacion_fsi_rotor.md](docs/formulations/teoria_formulacion_fsi_rotor.md). Publications are
cited by bibliographic reference (authors, venue, DOI) and by their key in
`docs/validation/references.yaml`, never by a file path: the repository does not
distribute articles. Copies used during development are local only, in the
gitignored `.sources/papers/`; `references.yaml` records which entries have a held
copy (`held_files`) and marks every field that still needs verification.

## Validation Summary

Validation follows a strict ladder set by the thesis plan
(`odd/tasks/thesis-validation-plan.md`); the evidence lives in the validation
store `docs/validation/` (levels in
[docs/validation/levels.yaml](docs/validation/levels.yaml)), driven by
`tools/validation_matrix.py` and governed by
[docs/validation-policy.md](docs/validation-policy.md). The acceptance ceiling
is 5 % relative error against an independent reference. No numbers are
repeated here; read them in the store and in
[docs/validation_closures.md](docs/validation_closures.md).

| level | scope | status |
|---|---|---|
| L1 | MITC shell elements against published results | closed |
| L2 | linear structural solver against published results and CalculiX | closed |
| L3 | blade structure against OpenFAST/BeamDyn decks | closed 2026-11-14 with declared deck differences D1-D4; re-closure on the windIO v2 pipeline planned (`odd/tasks/l3-v2-closure.md`) |
| L4 | built-in BEM (rigid blade) against OpenFAST AeroDyn | closed |
| L5 | coupled FSI, one-way and two-way, against OpenFAST and Zhou | open |
| L6 | preCICE + OpenFOAM, tower, wake | out of scope (CFD project) |

> Note: the benchmark tests under `tests/validation/benchmarks/` use tight
> tolerances; the documented failures there are element-accuracy limits, not
> solver regressions (see `scripts/known-red-full.txt` and
> `docs/origin_main_integration_2026-09-30.md`).

## Limitations

At this stage, users should expect some rough edges:

- documentation is still being expanded
- workflows are stronger for research and custom engineering use than for
  packaged end-user operation
- some advanced FSI setups depend heavily on external solver configuration,
  mesh quality, and restart discipline
- APIs and configuration details may evolve as the project matures

## License and Usage

This repository includes a restrictive research license in [LICENSE](LICENSE).

In short:

- use is allowed for non-commercial research, academic, educational, and
  evaluation purposes
- commercial use is not allowed without prior written permission
- contributions are welcome and governed by [CONTRIBUTING.md](CONTRIBUTING.md)

This is not an OSI-approved open-source license. If you later want broader
adoption, external packaging, or community growth, you should expect to revisit
this choice.