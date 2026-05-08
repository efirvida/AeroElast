# Tasks: Rust Migration of LinearDynamicFSIRotorInertialSolver

## Phase 1: Inertial Physics Utilities (6 tasks) ✅ COMPLETE

- [x] 1.1 Add `compute_rigid_body_acceleration_inertial()` to `rotor_physics.rs` — computes `a_ref = α×r + ω×(ω×r)` in global frame
- [x] 1.2 Add `compute_reference_load_vector()` to `rotor_physics.rs` — computes `F_ref = -M·a_ref` element-wise 
- [x] 1.3 Add unit tests for `compute_rigid_body_acceleration_inertial()` in `rotor_physics.rs` — pure rotation, pure angular accel, combined, zero at rest cases
- [x] 1.4 Add `update_elastic_stiffness(&mut self, k_vals: &[f64])` to `NewmarkStepper` in `dynamic_newmark.rs` — updates K and refactorizes
- [x] 1.5 Verify `MeshAssembler::update_node_coordinates()` exists (no work needed — already implemented)
- [x] 1.6 Unit test `NewmarkStepper::update_elastic_stiffness()` — 1-DOF spring-mass with K change `k₁ → k₂`

## Phase 2: InertialRotorFsiSolver Structure (5 tasks) ✅ COMPLETE

- [x] 2.1 Create `InertialRotorFsiConfig` struct in `rotor_inertial.rs` — fields: `fsi: FsiConfig`, rotation geometry, gravity, `k_update_interval`, `omega_rebuild_threshold`
- [x] 2.2 Create `DisplacementMode` enum in `rotor_inertial.rs` — variants: `Elastic`, `Total`
- [x] 2.3 Create `InertialRotorFsiSolver` struct in `rotor_inertial.rs` — fields: stepper, omega_provider, coords tracking, masses, config
- [x] 2.4 Implement `InertialRotorFsiSolver::new()` constructor — validates rotation_axis, clones coords_ref, normalizes axis
- [x] 2.5 Implement builder methods `with_initial_state()` and `with_step_callback()` in `InertialRotorFsiSolver`

## Phase 3: Setup Utilities (3 tasks) ✅ COMPLETE

- [x] 3.1 Add `extract_masses(lumped_diag: &[f64], dofs_per_node: usize) -> Vec<f64>` to `setup.rs` — extracts first translational DOF per node
- [x] 3.2 Add `rotate_mesh_coords(coords_ref, transforms, theta) -> Vec<f64>` to `setup.rs` — applies Rodrigues rotation `R(θ)` 
- [x] 3.3 Add `reassemble_k(assembler, coords_rotated) -> (rows, cols, vals)` to `setup.rs` — calls `update_node_coordinates()` + `assemble_k()`

## Phase 4: FSI Coupling Loop (8 tasks) ✅ COMPLETE

- [x] 4.1 Implement `expand_to_full()` helper in `InertialRotorFsiSolver` — maps reduced DOF array to full structural mesh
- [x] 4.2 Implement `scatter_node_forces()` helper — converts interface force vector to full DOF array
- [x] 4.3 Implement `reassemble_k_if_needed()` — checks ω² change + step interval, rotates mesh, reassembles K(θ), calls `update_elastic_stiffness()`
- [x] 4.4 Implement `run()` method skeleton — preCICE initialization, main coupling loop, checkpoint/restore pattern
- [x] 4.5 Add F_ref computation in `run()` — calls `compute_rigid_body_acceleration_inertial()` + `compute_reference_load_vector()`
- [x] 4.6 Add gravity handling in `run()` — constant global vector scattered to reduced DOFs
- [x] 4.7 Add DisplacementMode handling in `run()` — write elastic-only vs total displacement to preCICE based on enum
- [x] 4.8 Add K(θ) rebuild logic in `run()` — call `reassemble_k_if_needed()` at window convergence, update stepper

## Phase 5: PyO3 Bindings (5 tasks) ✅ COMPLETE

- [x] 5.1 Add `run_inertial_rotor_fsi_solver()` function to `aeroelast-py/src/lib.rs` — PyO3 wrapper mirroring `run_rotor_fsi_solver` pattern
- [x] 5.2 Export `run_inertial_rotor_fsi_solver` in `#[pymodule] fn aeroelast_py()` — add `m.add_function(wrap_pyfunction!(...))`
- [x] 5.3 Add `displacement_mode: &str` parameter — parse "elastic" | "total" string to `DisplacementMode` enum
- [x] 5.4 Add `all_node_masses: PyReadonlyArray1<f64>` parameter — nodal mass array from Python
- [x] 5.5 Call `extract_masses()` in wrapper to get masses from lumped diagonal (already integrated in 5.1)

## Phase 6: Integration Tests (10 tasks) — IN PROGRESS

- [x] 6.1 Create `tests/test_inertial_rotor_rust_parity.py` — smoke test file structure
- [x] 6.2 Add `TestInertialRotorRustBinding::test_signature_matches_expected()` — verify function signature
- [x] 6.3 Add `TestInertialRotorRustBinding::test_displacement_mode_elastic_accepted()` — verify "elastic" mode accepted
- [x] 6.4 Add `TestInertialRotorRustBinding::test_displacement_mode_total_accepted()` — verify "total" mode accepted
- [ ] 6.5 Add unit test: pure rigid rotation (ω const, no aero, no gravity) → verify `max(|u_e|) < 1e-9`
- [ ] 6.6 Add unit test: gravity-only rotor (ω const) → verify 1P torque modulation in frequency domain
- [ ] 6.7 Add analytical test: F_ref for single-node rotor → verify against hand calculation `m·ω²·r_perp + m·α×r`
- [ ] 6.8 Add integration test: 3-blade IEA-15 stub → compare Rust vs Python u_e, tau_aero, omega(t) within 1e-6 relative tolerance
- [ ] 6.9 Add checkpoint/restart test → restart from mid-simulation, verify continuation matches non-restarted run
- [ ] 6.10 Add preCICE contract test → verify inertial solver writes elastic-only displacement (not total rigid)

## Phase 7: Documentation & Polish (14 tasks)

- [ ] 7.1 Add module-level docstring to `crates/aeroelast-solvers/src/petsc/fsi/rotor_inertial.rs` — explain inertial frame formulation
- [ ] 7.2 Add docstring to `DisplacementMode` enum — explain "elastic" vs "total" and CFD coupling implications
- [ ] 7.3 Add docstring to `InertialRotorFsiConfig` — document all fields, units, and valid ranges
- [ ] 7.4 Add docstring to `InertialRotorFsiSolver::new()` — constructor contract and validation rules
- [ ] 7.5 Add docstring to `InertialRotorFsiSolver::run()` — full preCICE coupling loop behavior
- [ ] 7.6 Add docstring to `compute_rigid_body_acceleration_inertial()` in `rotor_physics.rs` — formula and coordinate frame
- [ ] 7.7 Add docstring to `compute_reference_load_vector()` in `rotor_physics.rs` — sign convention and DOF layout
- [ ] 7.8 Update `docs/teoria_formulacion_fsi_rotor.md` — add section on inertial solver formulation vs corotational
- [ ] 7.9 Update `docs/cli-reference.md` — document `LinearDynamicFSIRotorInertial` solver type and parameters
- [ ] 7.10 Add example YAML config for inertial solver in `examples/rotor-inertial-example.yaml`
- [ ] 7.11 Add Python wrapper docstring to `run_inertial_rotor_fsi_solver` in `aeroelast-py/src/lib.rs`
- [ ] 7.12 Remove temporary debug prints / comments from implementation files
- [ ] 7.13 Run `cargo fmt` and `cargo clippy` on all modified Rust files
- [ ] 7.14 Update `CHANGELOG.md` — add entry for inertial solver Rust migration

## Dependency Order

```
Phase 1 (physics utils)
  ↓
Phase 2 (solver struct)
  ↓
Phase 3 (setup helpers)
  ↓
Phase 4 (FSI loop)
  ↓
Phase 5 (Python bindings)
  ↓
Phase 6 (tests)
  ↓
Phase 7 (docs & polish)
```

**Critical paths**:
- 1.1-1.2 → 4.5 (F_ref physics)
- 1.4 → 4.8 (K(θ) update in stepper)
- 3.2-3.3 → 4.3 (geometry rotation + K reassembly)
- 5.1-5.2 → 6.* (Python binding → tests require it)
- 6.8 (parity test) → determines if migration is successful

## Exit Criteria

- [x] Rust implementation complete (Phases 1-5)
- [ ] All unit tests pass (6.5-6.7)
- [ ] Integration parity test shows <1e-6 relative error vs Python baseline (6.8)
- [ ] Documentation complete (Phase 7)
- [ ] Python wrapper functional and validated end-to-end
