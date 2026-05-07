# Inertial vs Corotational Solver Comparison Guide

## Current Implementation Status (May 6, 2026)

### ✅ Implemented (Phases 0-6, 33/44 tasks)

**Infrastructure & Building Blocks:**
- Phase 0: Config, enums, dispatch, backward compatibility
- Phase 1: Shared utilities (rotation, inertial forces, reference load)
- Phase 2: Solver skeleton (OmegaProvider, config, state tracking)
- Phase 3: Structural assembly on rotated geometry (K(θ), C(θ), BCs)
- Phase 4: Newmark integration (K_eff, F_eff, checkpoint/rollback)
- Phase 5: preCICE contract (fixed interface, elastic displacement)
- Phase 6: Omega dynamics and torque accounting

**Total:** 1918 lines of Python code, 24 methods in `LinearDynamicFSIRotorInertialSolver`

### ❌ NOT Implemented (Phases 7-10, 11/44 tasks)

**Critical for execution:**
1. `solve()` main method - **STUB ONLY** (raises `NotImplementedError`)
   - preCICE adapter integration
   - Time loop with sub-iterations
   - State management (u_e, v_e, a_e)
   - Checkpoint/restart flow
   
2. `_solve_fsi_step()` - per-timestep assembly and solve
   - Assemble K_eff, F_eff using Phase 3-4 methods
   - Solve linear system
   - Update Newmark velocity/acceleration
   
3. Integration with parent class
   - Checkpoint file I/O (inherited from `LinearDynamicFSI`)
   - VTK output writers
   - Performance metric logging
   
4. Validation & testing (Phase 7-10)
   - Unit tests for building blocks
   - Integration tests with mock preCICE
   - Benchmark comparison with corotational

---

## How to Compare with Corotational Solver

### Current Limitation

**You CANNOT execute an inertial solver case yet.** The `solve()` method raises:
```python
NotImplementedError:
    LinearDynamicFSIRotorInertialSolver.solve() is not yet implemented.
    Use LinearDynamicFSIRotorCorotationalSolver for production work.
```

### What You CAN Do Now

1. **Code review comparison:**
   - Compare method signatures between `rotor.py` (corotational) and `rotor_inertial.py`
   - Verify physical formulation differences (see design doc)
   - Review preCICE contract differences (Phase 5 methods)

2. **Unit test individual methods** (requires writing tests):
   ```python
   # Example: Test coordinate rotation
   solver = LinearDynamicFSIRotorInertialSolver(...)
   coords_rotated = solver._coord_transforms.rotate_point_cloud(coords, theta)
   
   # Example: Test reference load computation
   a_ref = solver._inertial_calculator.compute_rigid_body_acceleration_inertial(
       coords, omega, alpha
   )
   F_ref = solver._inertial_calculator.compute_reference_load_vector(
       M_diag, coords, omega, alpha
   )
   ```

3. **Static analysis:**
   - Verify imports and dependencies
   - Check method call graph
   - Confirm enum values and config structure

---

## Implementation Roadmap (To Enable Comparison)

### Minimal Path to Execution (~500-1000 lines)

#### 1. Implement `solve()` main loop
**File:** `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Lines:** ~200-300  
**Dependencies:**
- Inherit checkpoint methods from `LinearDynamicFSISolver`
- Use existing `Adapter` class from `base.py`
- Reuse output infrastructure from parent class

**Key sections:**
```python
def solve(self):
    # 1. Checkpoint/restart
    checkpoint_state = self._try_restore_checkpoint()
    
    # 2. Assembly
    (K, M), bc_manager = self._assemble_system_matrices()
    
    # 3. Extract interface
    interface_coords, interface_dofs = self._extract_interface_nodes()
    
    # 4. Initialize preCICE
    adapter = Adapter(participant, config_file, 0, 1)
    self._register_interface_at_reference_coords(adapter, "SolidMesh", interface_coords)
    dt = adapter.initialize()
    
    # 5. Time loop
    while adapter.is_coupling_ongoing():
        # 5a. Checkpoint
        if adapter.requires_writing_checkpoint():
            self._checkpoint_elastic_state()
        
        # 5b. Read forces
        F_aero = self._read_forces_from_precice_global(adapter, "SolidMesh", "Force")
        
        # 5c. Solve step
        u_e_new, v_e_new, a_e_new = self._solve_fsi_step(F_aero, dt, ...)
        
        # 5d. Write displacement
        self._write_elastic_displacement_to_precice(adapter, "SolidMesh", "Displacement", u_e_interface, theta)
        
        # 5e. Write omega
        self._write_omega_to_global_mesh(adapter, "GlobalSolidMesh", "AngularVelocity", omega)
        
        # 5f. Advance preCICE
        adapter.advance(dt)
        
        # 5g. Rollback or converge
        if adapter.requires_reading_checkpoint():
            self._rollback_elastic_state()
        else:
            # Update omega, rotate geometry, rebuild matrices
            tau_driving = self._compute_driving_torque(...)
            self._update_omega_after_converged_window(tau_driving, dt)
            theta += omega * dt
            self._rotate_structural_geometry_internal(theta)
            self._rebuild_assembler_with_rotated_geometry()
            
    adapter.finalize()
    return u_e, v_e, a_e
```

#### 2. Implement `_solve_fsi_step()`
**Lines:** ~100-150  
**Purpose:** Assemble and solve K_eff·u_e = F_eff for one sub-iteration

```python
def _solve_fsi_step(self, F_aero, dt, theta, omega, alpha, u_e_prev, v_e_prev, a_e_prev, bc_manager):
    # 1. Compute reference load
    F_ref = self._inertial_calculator.compute_reference_load_vector(
        M_diag, coords_rotated, omega, alpha
    )
    
    # 2. Compute gravity load
    F_gravity = self._compute_gravity_load_vector()
    
    # 3. Assemble effective system
    K_eff = self._assemble_inertial_effective_system()
    
    # 4. Assemble RHS
    F_eff = self._assemble_inertial_rhs(F_aero, F_gravity, F_ref, dt, u_e_prev, v_e_prev, a_e_prev)
    
    # 5. Solve linear system
    u_e_new = self._solve_linear_system(K_eff, F_eff, bc_manager)
    
    # 6. Update velocity and acceleration (Newmark formulas)
    v_e_new = self._newmark_velocity_update(u_e_new, u_e_prev, v_e_prev, a_e_prev, dt)
    a_e_new = self._newmark_acceleration_update(u_e_new, u_e_prev, v_e_prev, a_e_prev, dt)
    
    return u_e_new, v_e_new, a_e_new
```

#### 3. Integrate with parent class methods
**Lines:** ~50-100  
**Methods to override/extend:**
- `_try_restore_checkpoint()` - inherit from `LinearDynamicFSISolver`
- `_save_checkpoint()` - extend with `theta`, `omega`, `alpha`
- `_write_vtk_output()` - inherit from parent
- `_solve_linear_system()` - inherit from parent (PETSc KSP)

#### 4. Add Newmark update helpers
**Lines:** ~50-100  
**Purpose:** Compute velocity and acceleration from displacement

```python
def _newmark_velocity_update(self, u_new, u_prev, v_prev, a_prev, dt):
    # v_new = γ/(β·Δt) · (u_new - u_prev) + (1 - γ/β)·v_prev + Δt·(1 - γ/(2β))·a_prev
    beta = 0.25
    gamma = 0.5
    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    
    v_new = a1 * (u_new - u_prev) + (1.0 - gamma/beta) * v_prev + dt * (1.0 - gamma/(2.0*beta)) * a_prev
    return v_new

def _newmark_acceleration_update(self, u_new, u_prev, v_prev, a_prev, dt):
    # a_new = 1/(β·Δt²) · (u_new - u_prev - Δt·v_prev) - (1/(2β) - 1)·a_prev
    beta = 0.25
    a0 = 1.0 / (beta * dt * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    
    a_new = a0 * (u_new - u_prev) - a2 * v_prev - a3 * a_prev
    return a_new
```

---

## Testing Strategy

### Phase 7: Prototype Validation (Before Comparison)

**Test 1: Pure rigid rotation (no aero, no gravity)**
```yaml
solver:
  type: LinearDynamicFSIRotorInertial
  rotor:
    omega: 10.0  # rad/s, constant
    include_gravity: false
```
**Expected:** `||u_e|| < 1e-10` (elastic displacement remains zero)

**Test 2: Gravity-only rotor (prescribed omega)**
```yaml
solver:
  rotor:
    omega: 5.0
    gravity: [0, 0, -9.81]
```
**Expected:** 1P torque modulation, `τ_gravity = τ_gravity(t + 2π/ω)`

**Test 3: Reference load analytical verification**
- Single lumped mass at known radius
- Hand-calculate `a_ref = α×r + ω×(ω×r)`
- Verify `F_ref = -M·a_ref` matches

### Phase 8: Comparative Benchmarking

**Benchmark 1: Identical external loads**
- Run same case with corotational and inertial solvers
- Compare:
  - Displacement magnitude: `||u_total_corot|| ≈ ||u_e_inertial||` ?
  - Torque balance: `τ_aero + τ_gravity ≈ τ_total` ?
  - Omega evolution: `ω(t)` histories should match for `ComputedOmega`

**Benchmark 2: Fixed-omega convergence**
- Constant `ω = 10 rad/s`, no torque dynamics
- preCICE coupling with identical CFD participant
- Verify:
  - Steady-state aerodynamic loads match
  - Elastic displacement fields differ by rigid-body component only
  - Performance coefficients (Cp, Cq, Ct) match

---

## Expected Differences

### Physical Behavior

| Aspect | Corotational | Inertial |
|--------|--------------|----------|
| **Reference frame** | Rotating (with rotor) | Inertial (global) |
| **Unknown variable** | `u_local` (displacement in rotating frame) | `u_e` (elastic displacement in global frame) |
| **Fictitious forces** | K_SP, K_G, Coriolis, Euler terms | None (replaced by -M·a_ref) |
| **preCICE mesh** | Rotates with rotor (or transforms displacements) | Fixed at reference position |
| **Displacement output** | `u_global = R(θ)·u_local` | `u_e` (elastic only, CFD handles rigid rotation) |

### Numerical Behavior

- **Convergence rates:** Should be similar for implicit coupling
- **Matrix sparsity:** Inertial K(θ), C(θ) rebuilt each window (more assembly cost), corotational K_SP, K_G fixed (less assembly cost)
- **Stability:** Both should be unconditionally stable for β=0.25, γ=0.5
- **Accuracy:** Should match to machine precision for small rotations, diverge for large rotations if preCICE contract mismatch

---

## Decision Tree for Implementation

```
Do you want to...

├─ Execute a working inertial solver case NOW?
│  └─> Implement solve() + _solve_fsi_step() (~400 lines)
│      Estimated effort: 1-2 days for experienced developer
│
├─ Just validate the formulation is correct?
│  └─> Write unit tests for Phases 0-6 methods (~200 lines test code)
│      Estimated effort: 4-8 hours
│
├─ Compare performance with corotational?
│  └─> Need full solve() implementation first (see above)
│      Then add Phase 8 benchmarking scripts
│
└─ Understand the architecture only?
   └─> Read design doc and compare method signatures (0 code)
       Files: docs/rotor_inertial_solver_design.md
              src/aeroelast/solvers/fsi/rotor_inertial.py
              src/aeroelast/solvers/fsi/rotor.py
```

---

## Quick Start for Implementation

### Step 1: Create a feature branch
```bash
git checkout -b feature/inertial-solver-integration
```

### Step 2: Implement `solve()` main loop
- Copy structure from `rotor.py` lines 1141-1300
- Replace corotational-specific code with inertial equivalents
- Use Phase 0-6 methods already implemented

### Step 3: Implement `_solve_fsi_step()`
- Assemble K_eff using `_assemble_inertial_effective_system()`
- Assemble F_eff using `_assemble_inertial_rhs()`
- Solve with inherited `_solve_linear_system()` from parent

### Step 4: Test with simple case
- Square plate, constant omega, no coupling
- Verify it runs without crashing
- Check outputs are reasonable

### Step 5: Compare with corotational
- Same geometry, same omega, same loads
- Compare displacement fields, torque histories, performance metrics

---

## Current Git Status

```bash
Branch: rust_implementation
Latest commit: c452105 feat(fsi): add solve() stub and helper methods to inertial solver
File: src/aeroelast/solvers/fsi/rotor_inertial.py (1918 lines)
Status: Building blocks complete, orchestration stub added, NOT FUNCTIONAL END-TO-END
```

**Next commit should be:**
```
feat(fsi): implement solve() and _solve_fsi_step() for inertial solver
- Full preCICE time loop integration
- Per-step assembly and linear solve
- Checkpoint/restart flow
- VTK output integration
```

---

## Estimated Effort

| Task | Lines of Code | Estimated Time |
|------|---------------|----------------|
| Implement `solve()` | 200-300 | 4-6 hours |
| Implement `_solve_fsi_step()` | 100-150 | 2-3 hours |
| Newmark update helpers | 50-100 | 1-2 hours |
| Integration with parent class | 50-100 | 2-3 hours |
| Unit tests for solve() | 200-300 | 4-6 hours |
| **Total** | **600-950** | **13-20 hours** |

For an experienced developer familiar with the codebase.

---

## Contact & Documentation

- **Design document:** `docs/rotor_inertial_solver_design.md`
- **Task breakdown:** `docs/rotor_inertial_solver_tasks.md`
- **Implementation:** `src/aeroelast/solvers/fsi/rotor_inertial.py`
- **Corotational reference:** `src/aeroelast/solvers/fsi/rotor.py`

**Questions?** Check the implementation notes in each method's docstring. All Phase 0-6 methods include verification criteria and physical interpretation.
