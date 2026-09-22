# Feature: Remove 3D solid element support

## Objective

Delete all 3D solid element support and keep only plane (Quad4/Quad8/Quad9) and
shell (MITC3/MITC4 + composite) elements. The blade/rotor application path is
shell-only and never uses solids, so solids only add maintenance surface.

## Decisions taken (user, 2026-09-22)

1. **Delete everything solid, including `BoxVolumeMesh`** and its CLI branch.
   Accepted consequence: a simulation YAML using `type: "BoxVolumeMesh"` stops
   working, and `family: SOLID` is no longer a valid configuration.
2. **Keep the vendored NuMAD tree intact** (`src/aeroelast/models/blade/numad/`,
   10,199 lines). Its volumetric pipeline (`get_vol_mesh`,
   `create_offset_layers`, `get_vertex_normals`) stays even though it emits
   hex8, has no consumers, and loses its tests.
3. **Delete `tests/test_solid_elements.py` and `tests/test_vol_mesh.py`.**

## Accepted debt (explicitly acknowledged)

- NuMAD's volumetric pipeline survives with **no consumers, no tests, and
  output that the solver can no longer assemble** (hex8). It is dead code kept
  only to stay close to NuMAD upstream. Revisit if NuMAD is ever updated.
- The plane elements (Quad4/Quad8/Quad9) are kept by decision, but they are not
  exercised by the blade/rotor FSI path either — they are candidates for the
  same treatment later.

## Honest error accounting

| | Before | After |
|---|---|---|
| Rust test failures | 2 | **1** |
| Python non-green | 48 | **11** |
| **Total non-green** | **50** | **12** |

Removed: the solid `wedge6_me_trace_matches_volume` Rust failure plus the 37
`test_vol_mesh.py` Python failures.

**NOT fixed by this work — remaining real debt:**

- `elements::mitc4::tests::test_drill_membrane_operator_detects_drill_gradient`
  (`mitc4.rs:2480`) — a shell MITC4 drilling failure, real physics debt.
- 10 failures in `tests/test_rotor_physical_consistency.py` — a mathematically
  wrong assertion (`rel_error` is `δ/(R+δ)` but the test compares it to `δ/R`).
- 1 failure in `tests/test_mass_matrix_validation.py[tri3]` — 12.5% mass error
  vs a 5% tolerance. **`tri3` is a shell/surface element, not a solid**, so it
  is untouched by this work.

## Scope

### In scope
- Rust: solid element kernels, enums, reference-element trait, material spec,
  PyO3 batch functions and code arms, HDF5 mappings.
- Python: solid helpers, volume mesh generators, config values, solver and
  writer branches, CLI branches, exports.
- Tests: delete the two solid test files; remove the solid branch from
  `test_beam_4cases_parity.py`.
- Docs: README, `docs/cli-reference.md`, `.github/**`, `tools/`.

### Out of scope
- The vendored NuMAD tree (decision 2).
- Plane and shell elements.
- The remaining 12 non-green items listed above.
- Fixing the 34 pre-existing Rust warnings.

## Inventory

### Rust

| Item | Location |
|---|---|
| `solid.rs` (whole file, 1,715 lines, 7 tests) | `crates/aeroelast-core/src/elements/solid.rs` |
| `mod solid` | `elements/mod.rs:5` |
| `ReferenceElement3D` + 8 impls, `integrate_ke_3d_flat`, `integrate_me_3d_flat`, 9 solid tests | `elements/reference.rs:44-61,261,342,467-660` |
| 8 enum variants + `dofs_per_node` solid arm | `assembly/topology.rs:25-41,50-51` |
| `PrecomputedElem` solid variants, `MaterialSpec::Solid3D`, dispatch arms, `solid3d_en`/`solid3d_rho` | `assembly/assembler.rs:72-73,96-103,381-409,472-500,628-656,810-837` |
| `Material::Solid = 4` | `materials/mod.rs:14` |
| 8 enum variants + `node_count`/`assembler_code` solid arms | `crates/aeroelast-mesh/src/entities.rs:25-32,37-44,58-65` |
| HDF5 code mappings | `crates/aeroelast-mesh/src/io/hdf5.rs:143-150` |
| `MaterialSpec::Solid3D` match arm | `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs:474` |
| 16 `batch_ke_*`/`batch_me_*` solid functions | `crates/aeroelast-py/src/elements.rs` |
| 16 registrations | `crates/aeroelast-py/src/lib.rs:39-54` |
| Solid element code arms | `crates/aeroelast-py/src/mesh.rs:172-179`, `assembler.rs:89-96,153-160,465-472` |

### Python

| Item | Location |
|---|---|
| 11 solid helpers + `SOLID_CHECKS` | `core/mesh/utils.py:270,314,885,907,971,1004,1036,1079,1127,1179,1223-1341` |
| `SOLID_ELEMENT_NODES_MAP` + solid `ElementType` variants | `core/mesh/entities.py:33-41,58-68` |
| `BoxVolumeMesh`, `CylinderVolumeMesh`, `MixedElementBeamMesh`, `PyramidTransitionMesh` | `core/mesh/generators.py:2394,2673,2859,3066` |
| `*SOLID SECTION`, C3D maps | `core/mesh/io/writers.py:376-384,395-403,1550-1560` |
| Solid type lists | `core/mesh/model.py:371-379,413-447` |
| VTK solid maps | `core/viewer.py:47-48,906,966` |
| `has_solid` branches | `postprocess/stress_recovery.py:429,523,710,769` |
| `ElementFamily.SOLID`, validation, warning, `MeshGeneratorType.BOX_VOLUME` | `core/config.py:211,233,595-597,1374-1375` |
| `_FAMILY_PROPERTIES`/`_FAMILY_VECTOR_FORM` solid entries | `core/assembler.py:156,162,192` |
| Solid quality check branch | `solvers/fsi/runner.py:451-455` |
| `has_solid` branches | `solvers/fsi/rotor.py:2232`, `linear_dynamic.py:686` |
| CLI branches + templates | `cli/run_fsi.py:286-344`, `cli/run_bem_fsi.py:124,207` |
| Exports | `aeroelast/__init__.py:32,126`, `core/mesh/__init__.py:53,86,103,129` |
| CCX modal tool reference | `tools/postprocess_ccx_modal.py:28` |

### Tests and docs

| Item | Location |
|---|---|
| Delete (9 solid tests) | `tests/test_solid_elements.py` |
| Delete | `tests/test_vol_mesh.py` |
| Remove HEXA8 branch | `tests/test_beam_4cases_parity.py:247-480` |
| Update | `README.md:23,346`, `docs/cli-reference.md:290-851`, `.github/copilot-instructions.md:112-160`, `.github/agents/fem-reviewer.agent.md:13` |

## Execution slices

**Order is Rust first, then Python.** The Rust removal is compiler-driven (the
compiler enumerates every non-exhaustive match), which makes it the safest step
and the best way to discover sites. No Python production code calls the solid
batch functions, so the Python side is not blocked by it.

**Why the Rust work is sliced along crate boundaries.** A first attempt gave one
writer all four Rust crates at once (~16 files). It ran out of capacity partway:
it left `aeroelast-core/src/assembly/assembler.rs` (145 solid references, the bulk
of the work) untouched, never reached `aeroelast-mesh/src/entities.rs`,
`io/hdf5.rs`, `rotor_fsi.rs` or `linear_elastic.rs`, and reformatted unrelated
files with a formatter it had been told not to run. Each slice below is small
enough for one writer run.

**Why core and py cannot be committed separately.** Removing the `ElemType` solid
variants from `aeroelast-core` immediately breaks `aeroelast-py`, whose code arms
map integer codes onto those variants. So `cargo check -p aeroelast-core` passing
does NOT mean the workspace compiles: slice A is verified with core's own gate and
slice B closes the workspace. **Do not commit between A and B.**

| Slice | Content | Gate |
|---|---|---|
| A | `aeroelast-core` only (6 files): `elements/solid.rs`, `elements/mod.rs`, `elements/reference.rs`, `assembly/topology.rs`, `assembly/assembler.rs`, `materials/mod.rs` | `cargo check -p aeroelast-core` → 0 errors; `cargo test -p aeroelast-core` → exactly 1 failure (the known mitc4 drill). `-p aeroelast-py` is EXPECTED to fail here. |
| B | `aeroelast-mesh` + `aeroelast-solvers` + `aeroelast-py`: `mesh/src/entities.rs`, `mesh/src/io/hdf5.rs`, `solvers/.../rotor_fsi.rs`, `solvers/.../linear_elastic.rs`, `py/src/elements.rs`, `py/src/lib.rs`, `py/src/mesh.rs`, `py/src/assembler.rs` | `cargo check -p aeroelast-py` in BOTH feature configurations → 0 errors; exactly 31 `wrap_pyfunction!` + 6 `add_class` (37 registrations) |
| C | Python mesh layer, config, assembler maps, solvers, CLI, exports, tools | `python -c "import aeroelast"`; no dangling reference to a removed symbol |
| D | Delete `tests/test_vol_mesh.py`; fix `tests/test_beam_4cases_parity.py:277`; clean the stale solid prose in `postprocess/stress_recovery.py`; update README, `docs/cli-reference.md` and `.github/**` | full `pytest -m "not slow"`; expected **37 failed / 330 passed / 27 skipped / 11 errors**, the pre-existing baseline |

## Progress

| Slice | State |
|---|---|
| A | done, committed with B |
| B | done, committed with A as `7b295f3` |
| C | done, split in two: `dd2e05a` (callers) and `9230ea2` (definitions) |
| D | **not started** |

### Two things the slices learned that the plan did not predict

1. **`tests/test_solid_elements.py` was deleted in slice C, not D.** All 59 of its cases
exercise code removed by C — the element codes no longer resolve — so it belongs with
the removal. Committing 59 red tests would have left a commit that `git bisect` flags
and no reviewer can attribute. With it gone the suite returns to its exact baseline.
2. **`core/mesh/io/readers.py` had to be added to slice C's surfaces.** Its module-level
`MESHIO_TYPE_MAP` references all 9 solid `ElementType` variants and it is imported at
`import aeroelast` time, so deleting the variants without it breaks the import. Two
further removals were also authorised: `volumetric_remesh` and `check_mesh_quality`,
which read as generic helpers but are solid-only with zero consumers.

### Slice D remaining work

- Delete `tests/test_vol_mesh.py` (the accepted debt recorded above: it covers NuMAD's
  volumetric pipeline, which is kept, but its 37 failures are fixture-driven and it has
  no consumer).
- `tests/test_beam_4cases_parity.py:277` still references `ElementType.hexahedron` inside
  a helper used only by slow-marked tests. It does not affect the non-slow counts.
- `src/aeroelast/postprocess/stress_recovery.py` still has docstrings naming solid
  element families (lines 16, 299, 316-323, 385). No code identifier dangles; this is
  prose residue.
- README, `docs/cli-reference.md` and `.github/**` still advertise the removed volume
  generators and the `SOLID` family.

Every writer task for these slices must state explicitly: **never run `cargo fmt`,
`rustfmt`, or any formatter**, touch only the named files, and report
`git diff --stat` so scope violations are visible.

## Verification commands

```bash
export CONDA_PREFIX=/home/efirvida/miniconda3/envs/aeroelast-dev
export PETSC_DIR=$CONDA_PREFIX SLEPC_DIR=$CONDA_PREFIX HDF5_DIR=$CONDA_PREFIX
export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig
cargo test  -p aeroelast-core --manifest-path crates/Cargo.toml
cargo check -p aeroelast-py   --manifest-path crates/Cargo.toml
cargo check -p aeroelast-py   --no-default-features --manifest-path crates/Cargo.toml
~/miniconda3/envs/aeroelast-dev/bin/python -m pytest -m "not slow" -q -p no:cacheprovider
```

## Acceptance criteria

- [ ] No `solid`, `hexa`, `tetra`, `wedge`, `pyramid` symbol remains in `crates/` or `src/` outside the vendored NuMAD tree.
- [ ] `cargo test -p aeroelast-core` reports exactly 1 failure (the mitc4 drill).
- [ ] `cargo check` passes in both feature configurations.
- [ ] `_aeroelast` exposes exactly 37 public symbols (53 minus the 16 solid batches).
- [ ] `pytest -m "not slow"` has no new failures; total non-green is 11.
- [ ] No dangling import or reference to a removed symbol.

## Tasks

- [x] A Remove solid support from `aeroelast-core` (6 files). Verified by core's own gate; workspace expected to be broken until B.
- [x] B Remove the remaining solid referrers from `aeroelast-mesh`, `aeroelast-solvers` and `aeroelast-py`; close the workspace. Committed with A as `7b295f3`.
- [x] C Remove the Python callers (`dd2e05a`) and then the Python definitions (`9230ea2`).
- [ ] D Delete `test_vol_mesh.py`, fix the HEXA8 helper in `test_beam_4cases_parity.py`, clean the stale solid prose in `stress_recovery.py`, update the docs.
- [ ] V Final verification and commit.
