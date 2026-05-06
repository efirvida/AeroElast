# Tasks: LinearDynamicFSIRotor Inertial Solver

## Phase 0: Rename Current Solver and Preserve Compatibility

- [x] 0.1 **Update** `src/aeroelast/core/config.py` — add explicit solver types for `LinearDynamicFSIRotorCorotational` and `LinearDynamicFSIRotorInertial`. Keep `LinearDynamicFSIRotor` available as the legacy/default entry point. Verify: YAML validation accepts the new names and old cases remain valid.
- [x] 0.2 **Rename** the main class in `src/aeroelast/solvers/fsi/rotor.py` from `LinearDynamicFSIRotorSolver` to `LinearDynamicFSIRotorCorotationalSolver`. Verify: import path still works through a temporary alias.
- [x] 0.3 **Add** a compatibility alias in `src/aeroelast/solvers/fsi/rotor.py` so `LinearDynamicFSIRotorSolver = LinearDynamicFSIRotorCorotationalSolver`. Verify: existing tests and runtime imports do not change behavior.
- [x] 0.4 **Update** lazy exports in `src/aeroelast/solvers/fsi/__init__.py` — expose `LinearDynamicFSIRotorCorotationalSolver`, keep `LinearDynamicFSIRotorSolver`, and reserve/export `LinearDynamicFSIRotorInertialSolver`. Verify: all three names are importable from `aeroelast.solvers.fsi`.
- [x] 0.5 **Update** `src/aeroelast/solvers/fsi/runner.py` — dispatch `LinearDynamicFSIRotor` to the corotational solver for now, and add explicit dispatch branches for the new solver names. Verify: runner creates the expected class instance for all three solver type strings.
- [x] 0.6 **Update** docs in `docs/cli-reference.md` and `docs/teoria_formulacion_fsi_rotor.md` — mark the current implementation as corotational and note that the inertial solver is a separate path under development. Verify: terminology is consistent across both documents.

## Phase 1: Shared Utilities for Internal Rigid Rotation

- [x] 1.1 **Extract or add** a reusable helper in `src/aeroelast/solvers/fsi/corotational.py` or a new shared utility module to compute Rodrigues rotation matrices and rotate 3D point clouds about an arbitrary axis/center. Verify: unit test round-trip rotation on coordinates and vectors passes at 1e-12 tolerance.
- [x] 1.2 **Add** a helper that computes rigid-body nodal acceleration

  $$
  \ddot{\hat{\mathbf{x}}}_i = \boldsymbol{\alpha} \times \mathbf{r}_i + \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}_i)
  $$

  for all nodes in global coordinates. Verify: analytical test for a node on a circular trajectory matches closed-form centrifugal/Euler acceleration.
- [x] 1.3 **Add** a helper that converts rigid-body nodal accelerations into a reduced DOF load vector `F_ref = -M * a_ref`, respecting the repo DOF stride rules. Verify: translational DOFs receive the expected values and constrained DOFs are skipped correctly.
- [x] 1.4 **Decide and document** where the rotated structural geometry lives during a time window: cloned `MeshModel`, temporary assembler input arrays, or a dedicated helper object. Verify: the design decision is reflected in code comments/docstrings where the geometry update occurs.

## Phase 2: New Inertial Solver Skeleton

- [x] 2.1 **Create** `src/aeroelast/solvers/fsi/rotor_inertial.py` — add `LinearDynamicFSIRotorInertialSolver` subclassing the same FSI base path used by the rotor solver family. Verify: module imports without triggering preCICE or PETSc runtime errors at import time.
- [x] 2.2 **Implement** `_init_rotor_config()` in `src/aeroelast/solvers/fsi/rotor_inertial.py` by reusing the existing `OmegaProvider` parsing rules from the corotational solver. Verify: `ConstantOmega`, `RampedOmega`, `ComputedOmega`, and `RampedComputedOmega` map exactly as in the current solver.
- [x] 2.3 **Implement** `_init_solver_config()` and `_init_state_tracking()` for the inertial solver, mirroring only the parts that remain valid without `K_G`, `K_SP`, and `G_cor`. Verify: object construction succeeds from a minimal rotor YAML dict.
- [x] 2.4 **Add** clear solver-level docstrings in `src/aeroelast/solvers/fsi/rotor_inertial.py` stating that the unknown is elastic displacement in global coordinates over a rigidly rotated reference. Verify: the docstring does not reuse corotational terminology incorrectly.

## Phase 3: Structural Assembly on Internally Rotated Geometry

- [ ] 3.1 **Implement** internal rotation of all structural node coordinates at the start of each FSI window using the representative window kinematics `theta_target`, `omega_window`, `alpha_window`. Verify: rotated coordinates preserve pairwise distances and the rotor radius remains constant under pure rigid rotation.
- [ ] 3.2 **Build** a safe path to assemble `K(theta)` on the rotated geometry for the prototype implementation, even if that requires reconstructing the Rust-backed assembler per time window. Verify: a small integration test shows `assemble_stiffness_matrix()` changes with orientation for anisotropic/composite cases or remains numerically identical for isotropic symmetry cases where expected.
- [ ] 3.3 **Reuse** the original lumped or consistent mass matrix `M` without rebuilding it. Verify: `M` before and after geometry rotation is identical within floating-point tolerance.
- [ ] 3.4 **Rebuild** Rayleigh damping as `C(theta) = eta_m M + eta_k K(theta)` whenever `K(theta)` changes. Verify: with `eta_k = 0`, the damping matrix remains identical across orientations; with `eta_k != 0`, the orientation-dependent branch is exercised.
- [ ] 3.5 **Ensure** root boundary conditions remain elastic constraints only, not time-dependent rigid-body displacement constraints. Verify: constrained DOFs at the root stay homogeneous in the elastic unknown space.

## Phase 4: Newmark Step for the Inertial Formulation

- [ ] 4.1 **Implement** assembly of the inertial solver effective system

  $$
  [\mathbf{K}_{eff}] = [\mathbf{K}(\theta)] + a_0[\mathbf{M}] + a_1[\mathbf{C}(\theta)]
  $$

  without adding `K_G`, `K_SP`, or `G_cor`. Verify: the effective operator matches the standalone `LinearDynamicFSI` structure when rotation is disabled.
- [ ] 4.2 **Assemble** the RHS using aerodynamic forces in global coordinates, gravity in global coordinates, the rigid-body reference load `-M a_ref`, and the standard Newmark history terms. Verify: sign convention of `-M a_ref` matches analytical rigid rotation tests.
- [ ] 4.3 **Keep** the structural state variable as elastic displacement `u_e`, not total displacement. Verify: logging and checkpoint code do not silently reinterpret `u_e` as absolute displacement.
- [ ] 4.4 **Implement** the state update and checkpoint/rollback protocol for `u_e`, `v_e`, `a_e`, `theta`, `omega`, and `alpha`. Verify: rollback restores the same elastic state and rigid-body kinematics used at the start of the window.

## Phase 5: preCICE Contract with Fixed Interface Mesh

- [ ] 5.1 **Register** interface vertices in preCICE using the original reference coordinates only; do not rotate `SolidMesh` coordinates on the preCICE side. Verify: repeated windows do not attempt to redefine or move coupling vertices.
- [ ] 5.2 **Read** aerodynamic loads from preCICE directly in global coordinates, without applying `R^T(\theta)` as in the corotational solver. Verify: the inertial solver force path contains no force-frame transform.
- [ ] 5.3 **Write** only elastic displacement in global coordinates back to `SolidMesh`:

  $$
  \mathbf{u}_{fsi} = \mathbf{x} - \hat{\mathbf{x}} = \mathbf{u}_e^{glob}
  $$

  Verify: an internal assertion or test fails if total rigid rotation is accidentally added to the outgoing field.
- [ ] 5.4 **Retain** `GlobalSolidMesh` output for representative angular velocity exactly as in the current solver. Verify: a fixed-omega case writes the same scalar history as the corotational solver.
- [ ] 5.5 **Add** an explicit compatibility note or guard for participants that expect total displacement instead of elastic displacement. Verify: unsupported coupling modes fail with a clear error message rather than silent wrong physics.

## Phase 6: Omega Dynamics and Torque Accounting

- [ ] 6.1 **Reuse** the existing `OmegaProvider` update workflow so that `omega` continues to evolve only after a converged FSI window. Verify: `ComputedOmega` and `RampedComputedOmega` preserve the same state machine semantics as the current solver.
- [ ] 6.2 **Compute** driving torque using only external forces: aerodynamic torque, gravity torque, and shaft torque. Verify: no rigid-body reference load contribution is included in `tau_driving`.
- [ ] 6.3 **Log** both aerodynamic torque and total structural response torque for comparison with the corotational solver. Verify: output files/report columns clearly distinguish `tau_aero`, `tau_gravity`, `tau_total`, and any non-aero residual.
- [ ] 6.4 **Use** the representative window angular velocity consistently for interface output and power coefficient postprocessing. Verify: `cp`, `cq`, `ct`, and `tsr` remain coherent with the current postprocess conventions.

## Phase 7: Prototype Validation

- [ ] 7.1 **Add** a unit/integration test for pure rigid rotation with zero aero and zero gravity: elastic displacement should remain numerically near zero for all time steps. Verify: `max(|u_e|)` stays below a tight tolerance.
- [ ] 7.2 **Add** a gravity-only rotor test with prescribed `omega`: the response should show the expected 1P modulation in torque or displacement. Verify: the dominant frequency in the signal is 1P.
- [ ] 7.3 **Add** a test for the rigid-body reference load `-M a_ref` against a simple analytical rotor with one or a few lumped nodes. Verify: nodal loads match hand calculations.
- [ ] 7.4 **Add** a preCICE-contract regression test or mock verifying that the inertial solver writes elastic displacement only, while the corotational solver still writes transformed displacement from its own formulation. Verify: the two outgoing interface fields differ exactly by the rigid-body component in a controlled test.
- [ ] 7.5 **Add** checkpoint/restart coverage for the inertial solver. Verify: restarting from a converged state reproduces the same subsequent step history within tolerance.

## Phase 8: Comparative Benchmark Against Corotational Solver

- [ ] 8.1 **Create** a reproducible benchmark case shared by both solvers with identical mesh, materials, aero forcing, damping, and omega provider configuration. Verify: input YAML differs only in the solver type and any formulation-specific flags.
- [ ] 8.2 **Measure** and report per-window timings: geometry rotation, assembler reconstruction/update, matrix assembly, factorization, solve, and total wall time. Verify: benchmark output includes a structured timing breakdown for both solvers.
- [ ] 8.3 **Compare** physical outputs between both formulations: tip displacement, aerodynamic torque, total torque, `omega(t)`, and key frequency content. Verify: results are saved in a side-by-side report or CSV suitable for paper figures.
- [ ] 8.4 **Document** the expected mismatch caused by omitting `K_G`/`K_SP` in the inertial prototype, so that the first comparison is interpreted as framework comparison, not full physics equivalence. Verify: benchmark notes and docs state this limitation explicitly.

## Phase 9: Performance Follow-Up in Rust

- [ ] 9.1 **Profile** whether assembler reconstruction dominates the inertial solver cost. Verify: profiling data identifies whether the bottleneck is geometry transfer, assembler creation, COO assembly, or PETSc factorization.
- [ ] 9.2 **If needed, add** a Rust-side API in `crates/aeroelast-core/src/assembly/assembler.rs` and `crates/aeroelast-py/src/lib.rs` to refresh node coordinates without rebuilding topology/material state. Verify: new API reproduces the same `K(theta)` as full reconstruction.
- [ ] 9.3 **Evaluate** a matrix-transform shortcut

  $$
  \mathbf{K}(\theta) = \mathbf{T}(\theta)^T \mathbf{K}_0 \mathbf{T}(\theta)
  $$

  for rigid-body rotation as an optimization candidate. Verify: transformed matrices match direct reassembly within tolerance on representative cases.
- [ ] 9.4 **Only after parity**, consider porting the inertial FSI loop to Rust, similar to the current corotational fast path. Verify: Rust and Python inertial paths produce the same results on the shared benchmark.

## Phase 10: Product Decision and Cleanup

- [ ] 10.1 **Decide** whether `LinearDynamicFSIRotor` should remain an alias to the corotational solver, switch to the inertial solver, or be deprecated in favor of explicit naming. Verify: decision is documented in `docs/cli-reference.md` and release notes/changelog if applicable.
- [ ] 10.2 **Update** user-facing docs, examples, and any template YAMLs to expose both formulations and their tradeoffs. Verify: users can discover both solver modes from the CLI docs without reading internal design notes.
- [ ] 10.3 **Remove** temporary compatibility glue only if a migration decision is made and downstream cases are updated. Verify: no example, simulation template, or test still depends on obsolete solver names.

## Dependency Order

```text
Phase 0
  -> Phase 1
  -> Phase 2
  -> Phase 3
  -> Phase 4
  -> Phase 5
  -> Phase 6
  -> Phase 7
  -> Phase 8
  -> Phase 9
  -> Phase 10
```

Critical dependencies inside the implementation:

- `0.1` to `0.5` must land before any runtime comparison work, otherwise both solver families cannot coexist cleanly.
- `1.2` and `1.3` must exist before `4.2`, because the rigid-body reference load is the core physical term of the inertial formulation.
- `3.2` and `3.4` must exist before `4.1`, because `K(theta)` and `C(theta)` define the effective system.
- `5.3` must be correct before any coupled benchmark, or the comparison against the corotational solver will be contaminated by an interface-contract mismatch.
- `8.4` must be documented before making product decisions from benchmark results.

## Exit Criteria

- [ ] The current rotor solver remains available explicitly as corotational.
- [ ] The new inertial solver runs a coupled FSI window with fixed `SolidMesh` and active `GlobalSolidMesh`.
- [ ] The inertial solver writes elastic displacement only, never total rigid-body motion, to preCICE.
- [ ] The rigid-body reference load `-M a_ref` is validated analytically.
- [ ] A reproducible A/B benchmark exists with timing and physics output for both formulations.