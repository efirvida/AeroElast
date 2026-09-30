# Feature: Rust PyO3 binding modularization

## Objective

Split `crates/aeroelast-py/src/lib.rs` (4,319 lines, 47 `#[pyfunction]` + 6
`#[pyclass]`, zero tests) into cohesive modules without changing a single
behavior. This is the top change hotspot in the repository (33 commits) and
the largest single Rust file.

Pure structural move: no logic change, no API change, no dependency change.

## Scope

### In scope
- Extract implementation items from `lib.rs` into 6 new modules.
- Keep a single `register_module` in `lib.rs` so the Python-facing public API
  stays visible in exactly one place.
- Keep `#[pymodule] fn _aeroelast` at the crate root (required for the cdylib
  init symbol).

### Out of scope
- Splitting `mitc4.rs` (2,837), `solid.rs` (1,715), `assembler.rs` (1,610),
  `dynamic_newmark.rs` (1,456), `mitc3.rs` (1,285) — separate work units.
- Removing dead code, fixing warnings, changing any algorithm.
- Touching any other crate.

## Target module map

```
crates/aeroelast-py/src/
├── lib.rs          mod decls + register_module + #[pymodule]      (~120 lines)
├── elements.rs     12 mitc3/mitc4 batch fns + 22 solid/quad batch
│                   fns + PyElementFamily                          (~1,050 lines)
├── assembler.rs    PyMeshAssembler (struct + impl + pymethods)    (~815 lines)
├── solvers.rs      coo_assembly, compute_nnz, linear_static_solve_coo,
│                   newmark_beta_solve_coo, PETSC_MAT_CAPSULE_NAME,
│                   petsc_assemble_matrix, petsc_modal_solve,
│                   modal_solve_coo, nonlinear_static_solve_coo,
│                   compute_rayleigh_auto                          (~700 lines)
├── mesh.rs         PyMeshModel                                    (~197 lines)
├── materials.rs    parse_material, PyOrthotropicMaterial, PyPly,
│                   PyLaminate                                     (~320 lines)
└── fsi.rs          run_linear_elastic_fsi, run_fsi_solver,
                    run_stress_stiffened_fsi_solver,
                    run_rotor_fsi_solver   (feature-gated)         (~920 lines)
```

## Item assignment (source line refs are pre-refactor)

| Item | Source line | Target module |
|---|---|---|
| `batch_ke_mitc3`, `batch_me_mitc3`, `batch_kt_mitc3`, `batch_fint_mitc3` | 27, 72, 118, 171 | `elements.rs` |
| `batch_ke_mitc4`, `batch_me_mitc4`, `batch_kt_mitc4`, `batch_fint_mitc4` | 227, 271, 316, 368 | `elements.rs` |
| `batch_ke_mitc3_composite`, `batch_me_mitc3_composite` | 422, 486 | `elements.rs` |
| `batch_ke_mitc4_composite`, `batch_me_mitc4_composite` | 537, 601 | `elements.rs` |
| `batch_ke_hexa8`, `batch_me_hexa8`, `batch_ke_tetra4`, `batch_me_tetra4` | 1933-2049 | `elements.rs` |
| `batch_ke_wedge6`, `batch_me_wedge6`, `batch_ke_tetra10`, `batch_me_tetra10` | 2087-2201 | `elements.rs` |
| `batch_ke_wedge15`, `batch_me_wedge15`, `batch_ke_hexa20`, `batch_me_hexa20` | 2233-2331 | `elements.rs` |
| `batch_ke_pyramid5`, `batch_me_pyramid5`, `batch_ke_pyramid13`, `batch_me_pyramid13` | 2363-2461 | `elements.rs` |
| `batch_ke_quad4`, `batch_me_quad4`, `batch_ke_quad8`, `batch_me_quad8` | 2497-2590 | `elements.rs` |
| `batch_ke_quad9`, `batch_me_quad9` | 2623, 2658 | `elements.rs` |
| `PyElementFamily` (`#[pyclass(name = "ElementFamily")]`) | 2897 | `elements.rs` |
| `PyMeshAssembler` (struct 1118, impl 1206, pymethods 1213) | 1117-1931 | `assembler.rs` |
| `coo_assembly`, `compute_nnz` | 652, 681 | `solvers.rs` |
| `linear_static_solve_coo`, `newmark_beta_solve_coo` | 710, 801 | `solvers.rs` |
| `PETSC_MAT_CAPSULE_NAME`, `petsc_assemble_matrix`, `petsc_modal_solve` | 924, 943, 974 | `solvers.rs` |
| `modal_solve_coo` | 1027 | `solvers.rs` |
| `nonlinear_static_solve_coo` | 3144 | `solvers.rs` |
| `compute_rayleigh_auto` | 4154 | `solvers.rs` |
| `PyMeshModel` (struct 2701, pymethods 2705) | 2700-2896 | `mesh.rs` |
| `parse_material` | 1129 | `materials.rs` |
| `PyOrthotropicMaterial`, `PyPly`, `PyLaminate` | 2911, 2952, 2979 | `materials.rs` |
| `run_linear_elastic_fsi`, `run_fsi_solver` | 3228, 3357 | `fsi.rs` |
| `run_stress_stiffened_fsi_solver`, `run_rotor_fsi_solver` | 3553, 3740 | `fsi.rs` |
| `register_module`, `#[pymodule] fn _aeroelast` | 4250, 4316 | **stays in `lib.rs`** |

## Constraints

1. `#[pymodule] fn _aeroelast` MUST remain at the crate root (`lib.rs`).
2. Registration stays centralized: one `register_module` in `lib.rs` with all
   47 `add_function!` + 6 `add_class!` calls. Submodules expose items as
   `pub(crate)`.
3. `parse_material` lives in `materials.rs` and is consumed by
   `assembler.rs` — it must be `pub(crate)`.
4. `fsi.rs` keeps the same feature gating the current code has.
5. Zero behavior change. Zero API change. No new dependency.

## Acceptance criteria

- [ ] `cargo check -p aeroelast-py --manifest-path crates/Cargo.toml` exits 0.
- [ ] `cargo check -p aeroelast-py --no-default-features --manifest-path crates/Cargo.toml` exits 0.
- [ ] Warning count is not worse than the 34 baseline.
- [ ] `_aeroelast` still exposes exactly the same 53 public symbols.
- [ ] `lib.rs` is under 150 lines.
- [ ] No file in the crate exceeds ~1,100 lines.

## Verification command

```bash
export CONDA_PREFIX=/home/efirvida/miniconda3/envs/aeroelast-dev
export PETSC_DIR=$CONDA_PREFIX SLEPC_DIR=$CONDA_PREFIX HDF5_DIR=$CONDA_PREFIX
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig
cargo check -p aeroelast-py --manifest-path crates/Cargo.toml
```

Symbol baseline (must match exactly):
`ElementFamily Laminate MeshModel OrthotropicMaterial Ply PyMeshAssembler
batch_fint_mitc3 batch_fint_mitc4 batch_ke_hexa20 batch_ke_hexa8
batch_ke_mitc3 batch_ke_mitc3_composite batch_ke_mitc4
batch_ke_mitc4_composite batch_ke_pyramid13 batch_ke_pyramid5 batch_ke_quad4
batch_ke_quad8 batch_ke_quad9 batch_ke_tetra10 batch_ke_tetra4
batch_ke_wedge15 batch_ke_wedge6 batch_kt_mitc3 batch_kt_mitc4
batch_me_hexa20 batch_me_hexa8 batch_me_mitc3 batch_me_mitc3_composite
batch_me_mitc4 batch_me_mitc4_composite batch_me_pyramid13
batch_me_pyramid5 batch_me_quad4 batch_me_quad8 batch_me_quad9
batch_me_tetra10 batch_me_tetra4 batch_me_wedge15 batch_me_wedge6
compute_nnz compute_rayleigh_auto coo_assembly linear_static_solve_coo
modal_solve_coo newmark_beta_solve_coo nonlinear_static_solve_coo
petsc_assemble_matrix petsc_modal_solve run_fsi_solver
run_linear_elastic_fsi run_rotor_fsi_solver run_stress_stiffened_fsi_solver`

## Tasks

- [x] T1 Extract `materials.rs` (`parse_material` + 3 pyclasses).
- [x] T2 Extract `mesh.rs` (`PyMeshModel`).
- [x] T3 Extract `solvers.rs` (9 fns + capsule const).
- [x] T4 Extract `assembler.rs` (`PyMeshAssembler`).
- [x] T5 Extract `elements.rs` (34 batch fns + `PyElementFamily`).
- [x] T6 Extract `fsi.rs` (4 runner fns, feature-gated).
- [x] T7 Verify both feature configurations and the 53-symbol baseline.

## Result

Completed. Final sizes: `elements.rs` 1,427 · `fsi.rs` 946 ·
`assembler.rs` 744 · `solvers.rs` 661 · `materials.rs` 313 · `mesh.rs` 198 ·
`lib.rs` 95 (was 4,319).

### Verification evidence

| Check | Baseline | After |
|---|---|---|
| `cargo check` (default) | 0 errors, py 8 warnings | 0 errors, py 7 |
| `cargo check --no-default-features` | 0 errors, py 7 warnings | 0 errors, py 6 |
| `#[pyfunction]` / `#[pyclass]` | 47 / 6 | 47 / 6 |
| `_aeroelast` public symbols (after maturin rebuild) | 53 | 53, set identical |
| Reference test subset | 72 passed, 1 skipped | 72 passed, 1 skipped |
| `register_module` block | — | byte-identical to pre-refactor |

The single failure in the broader test subset
(`test_mass_matrix_validation.py::test_mass_consistency[tri3]`, mass error
12.5%) predates this refactor and is unrelated to it.

### Deviations from plan

1. **`elements.rs` is 1,427 lines, above the ~1,100 soft cap.** The prescribed
   item set is ~1,412 moved lines, so the plan's ~1,050 estimate was wrong.
   Splitting it further needs an 8th module (shell vs solid/quad), which is
   outside this work unit's scope. Deferred to a follow-up task.
2. **Two lines added:** `#[cfg(feature = "fsi")]` on `mod fsi;` and on
   `use fsi::*;` in `lib.rs`. Without them an empty module plus an unguarded
   glob would raise `unused_imports`. Item-level attributes are unchanged.
3. **Visibility:** `pub(crate)` added to the 47 registered functions,
   `parse_material`, `PETSC_MAT_CAPSULE_NAME`, 5 pyclasses, and two `inner`
   fields that `assembler.rs` reads across the module boundary
   (`mesh::PyMeshModel.inner`, `materials::PyLaminate.inner`).
   `PyMeshAssembler` and its `pub` methods were left untouched so the baseline
   `private_interfaces` warnings are unchanged.
4. **`use nalgebra::{Matrix2, Matrix3};` dropped** from the crate root as dead:
   it was already emitting an unused-import warning, and `elements.rs` uses
   `nalgebra::Matrix3` fully qualified. This is the single cause of the
   py-warning improvement from 8 to 7.

### Deferred

- Split `elements.rs` into `elements_shell.rs` and `elements_solid.rs`
  (~1,412 lines is still too large for one file).
- Next targets: `mitc4.rs` (2,837), `solid.rs` (1,715),
  `assembly/assembler.rs` (1,610), `dynamic_newmark.rs` (1,456),
  `mitc3.rs` (1,285).
