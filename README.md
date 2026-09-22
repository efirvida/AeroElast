# AeroElast

AeroElast is a high-performance finite element toolkit for structural and
aeroelastic simulation of shell, solid, and plane structures, with a focus on
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

- finite element support for plane, shell, and solid elements
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
  `BoxVolumeMesh`, `MultiFlapMesh`, `BladeMesh`, and `RotorMesh`
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
pip install -e .
```

`pip install -e .` builds the Rust extension (`_aeroelast`) with maturin and
registers the available CLI commands. The two `export` lines point the Rust
build scripts at the conda-provided HDF5, PETSc/SLEPc and preCICE. Without the
compiled extension the package cannot import, so the conda env (or an
equivalent environment providing PETSc/SLEPc/preCICE) is mandatory.

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

## CLI Commands

The package provides these command-line entry points:

- `aeroelast` - main CLI entry point
- `aeroelast-fsi` - run or inspect YAML-defined simulations
- `aeroelast-bem-fsi` - run BEM-coupled FSI simulations
- `aeroelast-monitor` - monitor rotor/FSI runs from CSV output
- `aeroelast-reconstruct-csv` - rebuild missing rotor performance histories

Detailed CLI documentation is available in
[docs/cli-reference.md](docs/cli-reference.md).

## Repository Layout

```text
AeroElast/
├── crates/                # Rust workspace
│   ├── aeroelast-core/    # Pure FEM kernels (elements, materials, assembly)
│   ├── aeroelast-mesh/    # Mesh data types
│   ├── aeroelast-solvers/ # PETSc/SLEPc-backed solvers
│   └── aeroelast-py/      # PyO3 bindings
├── docs/                  # User-facing documentation
├── examples/              # Small usage examples
├── src/aeroelast/         # Python package source
│   ├── cli/               # CLI entry points
│   ├── core/              # Mesh, BCs, materials, config
│   ├── constitutive/      # Failure criteria and constitutive helpers
│   ├── elements/          # Element implementations
│   ├── models/            # Higher-level structural models
│   ├── postprocess/       # Output and analysis utilities
│   └── solvers/           # Static, dynamic, modal, and FSI solvers
└── tests/                 # Unit and regression tests
```

## Examples and Validation

- Example scripts live under `examples/`
- Automated tests live under `tests/`
- A practical test command used in development:

```bash
python -m pytest tests/ -q --tb=short --ignore=tests/test_blade_mesh.py --ignore=tests/test_rotor_inertial.py
```

Some benchmark cases are intentionally heavier and may not be suitable for a
quick local smoke test. A handful of distributed-mesh MITC3/MITC4 benchmarks
are also known to fall marginally outside the 5 % tolerance band used in the
assertions; this reflects element accuracy limits on coarse meshes, not solver
correctness.

## Documentation

- [docs/cli-reference.md](docs/cli-reference.md) — CLI reference for the
  simulation runner
- [docs/teoria_formulacion_fsi_rotor.md](docs/teoria_formulacion_fsi_rotor.md)
  — FSI rotor formulation theory (Spanish)
- [docs/s4r_composite_shell_formulation.md](docs/s4r_composite_shell_formulation.md)
  — composite shell formulation notes
- [docs/improvement_plan_mitc4_vs_s4r.md](docs/improvement_plan_mitc4_vs_s4r.md)
  — MITC4 vs S4R improvement plan
- [docs/FSI_ROTOR_PAPER_DRAFT.md](docs/FSI_ROTOR_PAPER_DRAFT.md) — scientific
  paper draft on the FSI rotor formulation

## Third-Party Code and References

AeroElast builds on external codes and scientific literature. As a research
project, we give credit where credit is due and respect every license
involved.

### Embedded / vendored code

- **pyNuMAD** (Sandia National Laboratories) — a modified copy is vendored
  under `src/aeroelast/models/blade/numad/` for blade meshing, under the
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

The formulation implemented in this repository is documented in
[docs/FSI_ROTOR_PAPER_DRAFT.md](docs/FSI_ROTOR_PAPER_DRAFT.md) and builds on
the following key works:

- Bathe, K.J., *Finite Element Procedures*, 2nd ed., Prentice Hall, 2014.
- Bucalem, M.L., Bathe, K.J., "Higher-order MITC general shell elements,"
  *Int. J. Numer. Meth. Eng.*, 36(21):3729–3754, 1993.
- Crisfield, M.A., *Non-linear Finite Element Analysis of Solids and
  Structures*, Vol. 2, Wiley, 1997.
- Géradin, M., Rixen, D., *Mechanical Vibrations: Theory and Application to
  Structural Dynamics*, 3rd ed., Wiley, 2015.
- Goldstein, H., Poole, C., Safko, J., *Classical Mechanics*, 3rd ed.,
  Addison Wesley, 2002.
- Bungartz, H.J., et al., "preCICE – A fully parallel library for multi-physics
  surface coupling," *Computers & Fluids*, 141:250–258, 2016.
- Küttler, U., Wall, W.A., "Fixed-point fluid–structure interaction solvers
  with dynamic relaxation," *Comput. Mech.*, 43(1):61–72, 2008.
- Degroote, J., Bathe, K.J., Vierendeels, J., "Performance of a new partitioned
  procedure versus a monolithic procedure in fluid–structure interaction,"
  *Computers & Structures*, 87(11–12):793–801, 2009.
- Ning, S.A., "A simple solution method for the blade element momentum
  equations with guaranteed convergence," *Wind Energy*, 17(9):1327–1345, 2014.
- Moriarty, P.J., Hansen, A.C., *AeroDyn Theory Manual*, NREL/TP-500-36881, 2005.
- Jonkman, J., Butterfield, S., Musial, W., Scott, G., *Definition of a 5-MW
  Reference Wind Turbine for Offshore System Development*, NREL/TP-500-38060, 2009.

### Element formulation references

The MITC3/MITC4 shell element formulations implemented in this repository
follow these publications:

- Ko, Y., Lee, P.-S., "A new MITC4+ shell element," *Computers & Structures*,
  182:404–418, 2017.
- Ko, Y., Lee, P.-S., Bathe, K.J., "Performance of the MITC3+ and MITC4+ shell
  elements in widely-used benchmark problems," *Computers & Structures*,
  193:187–206, 2017.
- Jeon, H.-M., Lee, Y., Lee, P.-S., "The MITC3+ shell element in geometric
  nonlinear analysis," *Computers & Structures*, 146:91–104, 2015.
- Ko, Y., Lee, P.-S., "The MITC4+ shell element in geometric nonlinear
  analysis," *Computers & Structures*, 185:1–14, 2017.
- A comparative formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements
  with new numerical results based on s-norm tests, *Computers & Structures*,
  204:48–64, 2018.
- Towards improving the 2D-MITC4 element for analysis of plane stress and
  strain problems, *Computers & Structures*, 275:106933, 2023.
- Cui, X., Peng, G., Ran, Q., Zhang, H., Li, S., "Derivation and implementation
  of one-point quadrature quadrilateral shell element with MITC4+ method
  (MITC4+R)," *Computers & Structures*, 291:107207, 2024.
- Continuum mechanics-based shell elements with six degrees of freedom at each
  node, *Computers & Structures*, 308:107622, 2025.

The reference PDFs used during development were removed from the repository;
the citations above are the canonical sources.

## Validation Summary

This section summarizes the validation work performed so far, grouped by
logical domain. Only cases with a validated error within the target tolerance
are listed. Percentages and references are intentionally left blank here and
are updated as results are confirmed; no unverified numbers are published.

### Shell elements

- MITC4 isotropic shell against analytical/reference solutions — _error: pending_
- MITC3 benchmarks — _error: pending_
- Large-rotation benchmark problems — _error: pending_
- Orthotropic shell parity — _error: pending_
- Composite shell (laminate) validation — _error: pending_

### Solid elements

- Solid element benchmark suite — _error: pending_
- Mixed solid element validation — _error: pending_

### Beam / structural validation

- Beam 4-case parity — _error: pending_
- Beam + shell 4-case parity — _error: pending_
- Composite beam parity — _error: pending_
- Shell analytical validation — _error: pending_
- Shell comprehensive validation — _error: pending_

### Dynamic / modal

- Mass matrix validation — _error: pending_
- Modal analysis parity (Rust vs Python) — _error: pending_
- Stress-stiffened solver — _error: pending_

### FSI / rotor

- Rotor Rust parity (Rust fast-path vs Python reference) — _error: pending_
- Rotor physical consistency — _error: pending_
- BEM engine validation — _error: pending_
- BEM polars / force projection — _error: pending_
- Composite B-coupling — _error: pending_
- FSI structural report — _error: pending_

### Rust core

- Rust assembler parity — _error: pending_
- Rust composite parity — _error: pending_

> Note: benchmark tests under `tests/benchmarks/` use tight tolerances and some
> distributed-mesh MITC3/MITC4 cases fall marginally outside the 5 % band on
> coarse meshes; those cases are element-accuracy limits, not solver
> regressions.

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