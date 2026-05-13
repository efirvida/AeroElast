/// Co-rotational FSI solver for rotating structures (rotor blades).
///
/// Implements the full preCICE coupling loop for a rotor in a rotating reference
/// frame, including:
/// - Rodrigues coordinate transforms (forces: global → rotating, disps: rotating → global)
/// - Centrifugal, Coriolis, and Euler inertial body forces on the RHS
/// - Spin-softening K_SP diagonal update on significant Δω
/// - Geometric stiffness K_G update at configurable intervals
/// - Four angular-velocity provider modes: Constant, Ramped, Computed, RampedComputed
/// - Checkpoint/restore for preCICE implicit IQN-ILS coupling
///
/// # Feature gate
/// Compiled only with `--features fsi`.

use std::time::Instant;

use aeroelast_core::assembly::assembler::MeshAssembler;

use crate::petsc::elasticity::dynamic_newmark::{NewmarkCheckpoint, NewmarkStepper};
use crate::petsc::fsi::force_utils::{apply_cap, apply_ramp};
use crate::petsc::fsi::linear_elastic::{FsiConfig, FsiError, FsiInitialState, FsiResult};
use crate::petsc::fsi::profiling::WindowProfiler;
use crate::petsc::fsi::rotor_physics::{
    OmegaCheckpoint, OmegaProvider, PerformanceCoefficients, RotorTransforms,
    build_coriolis_matrix, build_ksp_vals,
    compose_deformed_coords, compute_centrifugal_force, compute_euler_force,
    compute_gravity_force, compute_performance_coefficients, compute_thrust,
    compute_torque, max_radial_deflection_ratio,
};

// Re-export `omega_changed_significantly` so the test block (which uses `use super::*`)
// can continue to call it by its original name without modification.
pub(crate) use crate::petsc::fsi::rotor_physics::omega_changed_significantly;

// ── Step callback type ─────────────────────────────────────────────────────────

/// Per-step callback for the rotor FSI solver.
///
/// Invoked once per **converged** time window with the full rotor state.
///
/// Arguments:
/// `(t, time_step, dt, u_red, v_red, a_red, force_mag,
///   forces_iface_global, omega_state, alpha_state, theta, tau_aero, perf)`
pub type RotorStepCallback = Box<
    dyn Fn(
            f64,                     // t: converged time
            usize,                   // time_step: 1-based step index
            f64,                     // dt: time window width
            &[f64],                  // u_red: displacement (reduced DOFs)
            &[f64],                  // v_red: velocity (reduced DOFs)
            &[f64],                  // a_red: acceleration (reduced DOFs)
            f64,                     // force_mag: ||F_aero|| Frobenius
            &[f64],                  // forces_iface_global: aero forces, global frame
            f64,                     // omega_state: dynamic state after convergence
            f64,                     // alpha_state: dynamic state after convergence
            f64,                     // theta: cumulative rotation angle
            f64,                     // tau_aero: aerodynamic torque (scalar, about axis)
            PerformanceCoefficients, // ct, cp, cq, tsr
            &[f64],                  // applied_inertial_forces: centrifugal+Euler per-node (rotating frame)
        ) -> Result<(), FsiError>
        + Send,
>;

// ── Configuration ──────────────────────────────────────────────────────────────

/// Configuration for the rotor FSI solver.
#[derive(Debug, Clone)]
pub struct RotorFsiConfig {
    /// Base preCICE FSI configuration (participant, mesh, read/write data names).
    pub fsi: FsiConfig,

    // ── Rotor geometry ────────────────────────────────────────────────────────
    /// Rotation axis unit vector (will be normalized if not already).
    pub rotation_axis: [f64; 3],
    /// Rotation center in the global frame.
    pub rotation_center: [f64; 3],

    // ── Inertial forces ───────────────────────────────────────────────────────
    /// Gravity vector in the GLOBAL (inertial) frame [m/s²]. Use `[0,0,0]` to disable.
    pub gravity: [f64; 3],
    /// Enable centrifugal forces on the RHS.
    pub include_centrifugal: bool,
    /// Enable Coriolis forces on the RHS.
    pub include_coriolis: bool,
    /// Enable Euler forces on the RHS (angular-acceleration term).
    pub include_euler: bool,

    // ── Stiffness updates ─────────────────────────────────────────────────────
    /// Enable geometric stiffness K_G updates (requires `assembler` to be `Some`).
    pub include_kg: bool,
    /// Rebuild K_G every N converged steps. `0` and `1` both mean every step.
    pub kg_update_interval: usize,
    /// Enable corotational tangent stiffness rebuilds K_T(u) in the FSI loop.
    ///
    /// Default: `false` for strict backward compatibility (no K_T(u) updates).
    pub use_corotational_kt: bool,
    /// Rebuild corotational K_T every N converged steps.
    /// `0` and `1` both mean every step.
    pub kt_coro_update_freq: u32,
    /// Enable spin-softening K_SP diagonal update.
    pub include_ksp: bool,
    /// Rebuild K_SP when |Δω| exceeds this threshold [rad/s].
    /// Kept for schema backward-compatibility; no longer used as the rebuild
    /// criterion — the relative predicate `omega_changed_significantly` is used
    /// instead, with `omega_rebuild_rel_high` / `omega_rebuild_rel_low`.
    pub ksp_omega_threshold: f64,
    /// Relative |Δ(ω²)|/ω² threshold to trigger a K_G or K_SP rebuild.
    /// Applied to both matrices (both ∝ ω²). Default: `0.005` (0.5%).
    pub omega_rebuild_rel_high: f64,
    /// Relative |Δ(ω²)|/ω² threshold below which a rebuild is skipped once
    /// one was just performed (hysteresis lower band). Default: `0.003` (0.3%).
    pub omega_rebuild_rel_low: f64,

    // ── K_G(u) — Nivel 2 reassembly with deformed geometry ────────────────────
    /// Master toggle for the foreshortening-aware K_G assembly. When `true`,
    /// `K_G` is built at deformed coordinates `X_ref + u` (rotating frame: the
    /// reference geometry doesn't rotate, so `X_ref` is the canonical mesh).
    /// Default `false` for parity with the explicit-in-geometry baseline.
    pub kg_use_deformed_coords: bool,
    /// High-band relative threshold for the K_G deflection-driven rebuild
    /// trigger. When the radial deflection ratio drift since the last K_G
    /// rebuild exceeds this value, the next converged window will reassemble
    /// `K_G` with deformed coords. Default `0.01` (1%).
    pub kg_deflection_rebuild_rel_high: f64,
    /// Low-band relative threshold for the K_G deflection-driven rebuild
    /// trigger. Reserved for the hysteresis predicate; not consulted by the
    /// simple drift check. Default `0.005` (0.5%).
    pub kg_deflection_rebuild_rel_low: f64,

    // ── DOF layout ────────────────────────────────────────────────────────────
    /// DOFs per FEM node (typically 6 for shells).
    pub dofs_per_node: usize,

    // ── Performance coefficients ──────────────────────────────────────────────
    /// Fluid (air) density [kg/m³].
    pub fluid_density: f64,
    /// Freestream flow velocity [m/s].
    pub flow_velocity: f64,
    /// Rotor outer radius [m].
    pub rotor_radius: f64,

    // ── Optional GlobalSolidMesh for ω communication to the CFD solver ─────────
    /// preCICE mesh name for the GlobalSolidMesh that communicates ω.
    /// `None` disables this feature.
    pub omega_mesh_name: Option<String>,
    /// Data name to write ω to on the GlobalSolidMesh.
    pub omega_write_data: Option<String>,
    /// Single vertex coordinate for the GlobalSolidMesh (typically the rotation center).
    pub omega_vertex_coord: Option<[f64; 3]>,
}

impl RotorFsiConfig {
    /// Return the gravity magnitude squared.
    fn gravity_active(&self) -> bool {
        self.gravity.iter().map(|x| x * x).sum::<f64>() > 1e-24
    }

    /// Normalize the K_G cadence so legacy `0` means "every converged step".
    fn normalized_kg_update_interval(&self) -> usize {
        self.kg_update_interval.max(1)
    }

    /// Whether K_G should be rebuilt on this 1-based converged step index.
    fn should_update_kg_on_step(&self, time_step: usize) -> bool {
        time_step % self.normalized_kg_update_interval() == 0
    }

    /// Normalize corotational K_T cadence so legacy `0` means every step.
    fn normalized_kt_coro_update_freq(&self) -> u32 {
        self.kt_coro_update_freq.max(1)
    }

    /// Whether the corotational element frame should be refreshed this step.
    fn should_update_coro_frame_on_step(&self, coro_step_counter: u32) -> bool {
        coro_step_counter % self.normalized_kt_coro_update_freq() == 0
    }

    /// Whether K_T corotational update should run on this converged step.
    fn should_update_corotational_kt(&self, coro_step_counter: u32) -> bool {
        self.use_corotational_kt && self.should_update_coro_frame_on_step(coro_step_counter)
    }
}

// ── Solver struct ──────────────────────────────────────────────────────────────

/// Co-rotational FSI solver for rotating blades.
///
pub struct RotorFsiSolver {
    // ── Newmark structural integrator ─────────────────────────────────────────
    stepper: NewmarkStepper,

    // ── Configuration ─────────────────────────────────────────────────────────
    config: RotorFsiConfig,

    // ── preCICE interface ─────────────────────────────────────────────────────
    /// Flat interface node reference coordinates in the rotating frame, length `n_iface * 3`.
    ///
    /// Numerically these match the initial global coordinates at `t = 0`, but the
    /// values are treated as the co-rotational reference geometry `X0` when
    /// evaluating torques and deformed positions.
    interface_coords_global: Vec<f64>,
    /// Indices into the REDUCED DOF vector for each interface component
    /// (length `n_iface * dofs_dim`, typically 3 per node).
    interface_dofs: Vec<usize>,
    /// Number of spatial dimensions on the coupling mesh (3 for 3D).
    mesh_dims: usize,
    // ── Coordinate transforms ─────────────────────────────────────────────────
    transforms: RotorTransforms,

    // ── Angular velocity ──────────────────────────────────────────────────────
    omega_provider: OmegaProvider,
    /// Cumulative rotation angle θ = ∫ω dt [rad].
    theta: f64,

    // ── Inertial forces: full-mesh nodal data ─────────────────────────────────
    /// Reference (undeformed) coordinates of ALL structural nodes in the
    /// ROTATING frame, flat `[x0,y0,z0,…]`, length `n_all_nodes * 3`.
    all_node_coords: Vec<f64>,
    /// Per-node scalar lumped masses, length `n_all_nodes`.
    all_node_masses: Vec<f64>,

    // ── Matrix updates ────────────────────────────────────────────────────────
    /// Optional mesh assembler for K_G computation.
    assembler: Option<MeshAssembler>,
    /// Lumped mass diagonal in the FULL (unreduced) DOF space, length `n_full_dofs`.
    m_lumped_full: Vec<f64>,
    /// Free global DOF indices (ascending), length `n_reduced_dofs`.
    free_dofs_i32: Vec<i32>,
    /// Lookup: `full_to_red[global_dof]` = reduced index, or -1 if constrained.
    full_to_red: Vec<i32>,
    /// Total DOF count in the full (unreduced) system.
    n_full_dofs: usize,
    /// COO row indices of the REDUCED stiffness matrix (for K_SP and K_G expansion).
    k_rows: Vec<i32>,
    /// COO column indices of the REDUCED stiffness matrix.
    k_cols: Vec<i32>,
    /// Precomputed mapping: full K_G COO index → output index in reduced K sparsity
    /// (or -1 for entries outside the free-DOF set). Computed once in `new()`,
    /// eliminates all per-timestep HashMap allocations in `update_geometric_stiffness`.
    kg_coo_map: Vec<i32>,
    /// Precomputed mapping: full K/K_T COO index → output index in reduced K sparsity
    /// (or -1 for entries outside free DOFs). Computed once in `new()`.
    k_coo_map: Vec<i32>,
    /// ω² at the last K_SP rebuild. Used by `omega_changed_significantly`.
    /// Initialized to `NEG_INFINITY` to guarantee the first rebuild always fires.
    omega_sq_at_last_ksp: f64,
    /// Step number of last K_SP rebuild (for hysteresis logic).
    last_ksp_rebuild_step: usize,
    /// ω² at the last K_G rebuild. Used to skip rebuilds when ω is stable.
    /// Initialized to `NEG_INFINITY` to guarantee the first rebuild always fires.
    omega_sq_at_last_kg: f64,
    /// Step number of last K_G rebuild (for hysteresis logic).
    last_kg_rebuild_step: usize,
    /// Radial deflection ratio at the last K_G rebuild. Used by the K_G(u)
    /// trigger when `kg_use_deformed_coords` is enabled. Initialized to `0.0`
    /// so the first non-zero deflection above the high-band threshold triggers
    /// the first deformed-coords K_G assembly.
    deflection_ratio_at_last_kg_rebuild: f64,

    // ── Optional restart state ────────────────────────────────────────────────
    initial_state: Option<FsiInitialState>,
    /// Initial cumulative angle for restarts.
    initial_theta: f64,

    // ── Per-window inertial force accumulator (Fix #5) ───────────────────────
    /// Pre-allocated scratch buffer for the combined per-node inertial forces
    /// (centrifugal + Euler) in the rotating frame, flat `[fx0,fy0,fz0,…]`,
    /// length `n_all_nodes * 3`. Zeroed at the start of each sub-iteration;
    /// passed to the step callback at each converged window so Python
    /// `_step_cb` can consume it directly without recomputing.
    /// Coriolis is on the LHS (G_cor matrix) and is NOT accumulated here.
    inertial_scratch: Vec<f64>,
    /// Independent converged-window counter for corotational K_T cadence.
    coro_step_counter: u32,

    // ── Optional per-step callback ────────────────────────────────────────────
    step_callback: Option<RotorStepCallback>,
}

impl RotorFsiSolver {
    /// Create a new `RotorFsiSolver`.
    ///
    /// # Arguments
    ///
    /// * `stepper`                 — pre-assembled Newmark integrator (assembled with K, M, C)
    /// * `config`                  — rotor FSI configuration
    /// * `interface_coords_global` — flat interface node reference coords (`X0`) in the
    ///   rotating frame, `n_iface * 3`
    /// * `interface_dofs`          — reduced DOF indices for each interface component
    /// * `mesh_dims`               — coupling mesh spatial dimension (3)
    /// * `omega_provider`          — angular-velocity provider (Constant/Ramped/Computed/…)
    /// * `all_node_coords`         — ALL structural node reference coords (rotating frame), `n * 3`
    /// * `all_node_masses`         — per-node scalar masses, length `n_all_nodes`
    /// * `m_lumped_full`           — lumped mass diagonal in the full DOF space
    /// * `k_rows`, `k_cols`        — COO sparsity of the reduced K (for K_SP expansion)
    /// * `free_dofs_i32`           — sorted free global DOF indices
    /// * `n_full_dofs`             — total DOF count (unreduced)
    /// * `assembler`               — mesh assembler for K_G; `None` to disable K_G updates
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        stepper: NewmarkStepper,
        config: RotorFsiConfig,
        interface_coords_global: Vec<f64>,
        interface_dofs: Vec<usize>,
        mesh_dims: usize,
        omega_provider: OmegaProvider,
        all_node_coords: Vec<f64>,
        all_node_masses: Vec<f64>,
        m_lumped_full: Vec<f64>,
        k_rows: Vec<i32>,
        k_cols: Vec<i32>,
        free_dofs_i32: Vec<i32>,
        n_full_dofs: usize,
        assembler: Option<MeshAssembler>,
    ) -> Self {
        // Build full→reduced lookup.
        let mut full_to_red = vec![-1i32; n_full_dofs];
        for (i, &dof) in free_dofs_i32.iter().enumerate() {
            full_to_red[dof as usize] = i as i32;
        }

        // Precompute K_G COO → output index map (once, avoids per-timestep HashMap).
        let (kg_coo_map, k_coo_map) = if let Some(ref asm) = assembler {
            let dummy = vec![[0.0f64; 3]; asm.topology.n_elems];
            let (rows, cols, _) = asm.assemble_geometric_k(&dummy);
            let kg_map = crate::petsc::fsi::setup::build_kg_coo_map(
                &rows,
                &cols,
                &free_dofs_i32,
                &k_rows,
                &k_cols,
            );
            let (k_rows_full, k_cols_full, _) = asm.assemble_k();
            let k_map = crate::petsc::fsi::setup::build_kg_coo_map(
                &k_rows_full,
                &k_cols_full,
                &free_dofs_i32,
                &k_rows,
                &k_cols,
            );
            (kg_map, k_map)
        } else {
            (Vec::new(), Vec::new())
        };

        let axis = config.rotation_axis;
        let center = config.rotation_center;
        let transforms = RotorTransforms::new(axis, center);
        let initial_omega = omega_provider.initial_omega();
        let n_inertial = all_node_masses.len() * 3;

        Self {
            stepper,
            config,
            interface_coords_global,
            interface_dofs,
            mesh_dims,
            transforms,
            omega_provider,
            theta: 0.0,
            all_node_coords,
            all_node_masses,
            assembler,
            m_lumped_full,
            free_dofs_i32,
            full_to_red,
            n_full_dofs,
            k_rows,
            k_cols,
            kg_coo_map,
            k_coo_map,
            omega_sq_at_last_ksp: f64::NEG_INFINITY,
            last_ksp_rebuild_step: 0,
            omega_sq_at_last_kg: f64::NEG_INFINITY,
            last_kg_rebuild_step: 0,
            deflection_ratio_at_last_kg_rebuild: 0.0,
            initial_state: None,
            initial_theta: 0.0,
            inertial_scratch: vec![0.0f64; n_inertial],
            coro_step_counter: 0,
            step_callback: None,
        }
    }

    /// Set the initial structural state for a restart.
    pub fn with_initial_state(mut self, state: FsiInitialState, theta: f64) -> Self {
        self.initial_state = Some(state);
        self.initial_theta = theta;
        self
    }

    /// Register a per-step callback invoked after each converged time window.
    pub fn with_step_callback<F>(mut self, cb: F) -> Self
    where
        F: Fn(
            f64,
            usize,
            f64,
            &[f64],
            &[f64],
            &[f64],
            f64,
            &[f64],
            f64,
            f64,
            f64,
            f64,
            PerformanceCoefficients,
            &[f64],         // applied_inertial_forces: centrifugal+Euler per-node (rotating frame)
        ) -> Result<(), FsiError>
            + Send
            + 'static,
    {
        self.step_callback = Some(Box::new(cb));
        self
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    /// Expand a reduced-DOF vector to the full DOF space (zeros at constrained DOFs).
    fn expand_to_full(&self, v_red: &[f64]) -> Vec<f64> {
        let mut v_full = vec![0.0f64; self.n_full_dofs];
        for (i, &dof) in self.free_dofs_i32.iter().enumerate() {
            v_full[dof as usize] = v_red[i];
        }
        v_full
    }

    /// Extract per-node translational components from a full-DOF vector.
    ///
    /// Returns a flat array `[vx0,vy0,vz0,…]` of length `n_all_nodes * 3`.
    fn extract_node_translations(&self, v_full: &[f64]) -> Vec<f64> {
        let n_nodes = self.all_node_masses.len();
        let dofs = self.config.dofs_per_node;
        let mut out = vec![0.0f64; n_nodes * 3];
        for i in 0..n_nodes {
            for j in 0..3 {
                let gdof = i * dofs + j;
                if gdof < v_full.len() {
                    out[i * 3 + j] = v_full[gdof];
                }
            }
        }
        out
    }

    /// Scatter flat node forces (n_nodes × 3) into the reduced-DOF force vector.
    ///
    /// Only free DOFs receive contributions; constrained DOFs are silently skipped.
    fn scatter_node_forces(&self, f_node: &[f64], f_global: &mut [f64]) {
        let n_nodes = self.all_node_masses.len();
        let dofs = self.config.dofs_per_node;
        for i in 0..n_nodes {
            let nb = i * 3;
            for j in 0..3 {
                let gdof = i * dofs + j;
                if gdof < self.full_to_red.len() {
                    let red = self.full_to_red[gdof];
                    if red >= 0 {
                        f_global[red as usize] += f_node[nb + j];
                    }
                }
            }
        }
    }

    /// Build spin-softening K_SP and Coriolis gyroscopic matrix, apply to stepper.
    ///
    /// K_SP: spin-softening translational block operator
    /// `-ω² M_lump (I - n̂⊗n̂)` aligned to K sparsity.
    /// G_cor: antisymmetric Coriolis matrix (implicit treatment for stability).
    fn apply_ksp(&mut self, omega: f64) -> Result<(), FsiError> {
        if !self.config.include_ksp {
            return Ok(());
        }
        let ksp = build_ksp_vals(
            &self.m_lumped_full,
            &self.transforms.axis,
            omega,
            self.config.dofs_per_node,
            &self.k_rows,
            &self.k_cols,
            &self.free_dofs_i32,
            self.n_full_dofs,
        );
        
        // Build Coriolis gyroscopic matrix for implicit treatment.
        // This provides unconditional stability at high rotation rates.
        let (g_rows, g_cols, g_vals) = if self.config.include_coriolis {
            build_coriolis_matrix(
                &self.all_node_masses,
                &self.transforms.axis,
                omega,
                &self.free_dofs_i32,
                self.config.dofs_per_node,
            )
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };
        
        self.stepper
            .update_spin_softening_and_gyroscopic(&ksp, &g_rows, &g_cols, &g_vals)
            .map_err(FsiError::StepperError)?;

        self.omega_sq_at_last_ksp = omega * omega;

        if !g_vals.is_empty() {
            log::info!(
                "RotorFsi: K_SP + G_cor updated, ω={:.4} rad/s, |G_cor| entries={}",
                omega, g_vals.len()
            );
        }
        Ok(())
    }

    /// Rebuild K_G from centrifugal prestress when the configured cadence says it is due.
    ///
    /// Calls `stepper.set_initial_geometric_stiffness` (writes to `kg_base_vals`) so the
    /// centrifugal prestress is always included in K_eff regardless of runtime deformation.
    ///
    /// When `kg_use_deformed_coords` is enabled, K_G is reassembled at deformed
    /// coords `X_ref + u` (Nivel 2 foreshortening contribution). `u_red` is the
    /// current reduced displacement, used both for the deflection-driven rebuild
    /// trigger and for the K_G assembly itself. After the deformed-coords
    /// assembly the original `all_node_coords` are pushed back to the assembler
    /// so that any later assembly path sees the canonical reference geometry
    /// (Option A discipline: K stays linear, only K_G uses X_ref + u).
    fn update_kg_if_needed(
        &mut self,
        omega: f64,
        u_red: &[f64],
        time_step: usize,
    ) -> Result<(), FsiError> {
        if !self.config.include_kg {
            return Ok(());
        }
        if !self.config.should_update_kg_on_step(time_step) {
            return Ok(());
        }

        // Trigger #1: ω-change predicate (existing behavior).
        // ksp_omega_threshold² seeds the near-zero eps guard.
        let steps_since_rebuild = time_step.saturating_sub(self.last_kg_rebuild_step);
        let recently_rebuilt = steps_since_rebuild < 10;
        let eps_sq = self.config.ksp_omega_threshold * self.config.ksp_omega_threshold;
        let omega_significant = omega_changed_significantly(
            omega,
            self.omega_sq_at_last_kg,
            self.config.omega_rebuild_rel_high,
            self.config.omega_rebuild_rel_low,
            recently_rebuilt,
            eps_sq,
        );

        // Trigger #2: K_G(u) deflection-change predicate (Nivel 2 only).
        // When `kg_use_deformed_coords` is off the ratio is forced to 0 and
        // the predicate stays inert.
        let current_deflection_ratio = if self.config.kg_use_deformed_coords {
            let u_full = self.expand_to_full(u_red);
            max_radial_deflection_ratio(
                &u_full,
                &self.all_node_coords,
                &self.transforms.axis,
                &self.transforms.center,
                self.config.dofs_per_node,
            )
        } else {
            0.0
        };
        let deflection_significant = self.config.kg_use_deformed_coords && {
            let drift = (current_deflection_ratio
                - self.deflection_ratio_at_last_kg_rebuild)
                .abs();
            drift > self.config.kg_deflection_rebuild_rel_high
        };

        if !omega_significant && !deflection_significant {
            log::trace!(
                "RotorFsi: K_G rebuild skipped at step {}, ω={:.4} rad/s, def_ratio={:.4e} (no trigger)",
                time_step, omega, current_deflection_ratio,
            );
            return Ok(());
        }

        // Bail out if there's no assembler (K_G updates require one).
        if self.assembler.is_none() {
            return Ok(());
        }

        // ── Pre-extract values needed under the &mut self.assembler borrow ───
        let axis = self.config.rotation_axis;
        let center = self.config.rotation_center;
        let dofs_per_node = self.config.dofs_per_node;
        let n_red_nnz = self.k_rows.len();
        let kg_use_deformed = self.config.kg_use_deformed_coords;

        // Compose deformed coords / snapshot original coords up-front so the
        // mut-borrow scope below stays contained.
        let (coords_def, coords_orig) = if kg_use_deformed {
            let u_full = self.expand_to_full(u_red);
            let coords_def =
                compose_deformed_coords(&self.all_node_coords, &u_full, dofs_per_node);
            let coords_orig = self.all_node_coords.clone();
            (coords_def, coords_orig)
        } else {
            (Vec::new(), Vec::new())
        };

        // ── K_G assembly under a scoped &mut assembler borrow ─────────────────
        let kg_vals_full = {
            use aeroelast_core::assembly::assembler::MaterialSpec;
            let asm = self
                .assembler
                .as_mut()
                .expect("assembler presence checked above");

            let rho_per_elem: Vec<f64> = asm
                .materials
                .iter()
                .map(|m| match m {
                    MaterialSpec::Isotropic { rho, .. } => *rho,
                    MaterialSpec::Composite { mass_per_area, .. } => *mass_per_area,
                    MaterialSpec::PlaneStress { rho, .. } => *rho,
                    MaterialSpec::Solid3D { rho, .. } => *rho,
                })
                .collect();

            if kg_use_deformed {
                asm.update_node_coordinates(&coords_def);
                let (_, _, vals) =
                    asm.assemble_centrifugal_k(omega, axis, center, &rho_per_elem);
                // Restore X_ref so any later assembly path sees the canonical
                // reference geometry.
                asm.update_node_coordinates(&coords_orig);
                vals
            } else {
                let (_, _, vals) =
                    asm.assemble_centrifugal_k(omega, axis, center, &rho_per_elem);
                vals
            }
        };

        let kg_red = crate::petsc::fsi::setup::apply_kg_coo_map(
            &self.kg_coo_map,
            &kg_vals_full,
            n_red_nnz,
        );

        self.stepper
            .set_initial_geometric_stiffness(&kg_red)
            .map_err(FsiError::StepperError)?;

        self.omega_sq_at_last_kg = omega * omega;
        self.last_kg_rebuild_step = time_step;
        if kg_use_deformed {
            self.deflection_ratio_at_last_kg_rebuild = current_deflection_ratio;
        }

        let kg_norm: f64 = kg_red.iter().map(|x| x * x).sum::<f64>().sqrt();
        log::info!(
            "RotorFsi: K_G rebuilt at step {time_step}, ω={omega:.4} rad/s, K_G(u)={kg_use_deformed}, def_ratio={current_deflection_ratio:.4e}, ||K_G||_F = {kg_norm:.3e}"
        );
        Ok(())
    }

    /// Rebuild K_T (linear or corotational) from current displacement when enabled.
    fn update_kt_if_needed(&mut self, u_red: &[f64], coro_step_counter: u32) -> Result<(), FsiError> {
        if !self.config.should_update_corotational_kt(coro_step_counter) {
            return Ok(());
        }
        if self.assembler.is_none() {
            return Ok(());
        }

        let u_full = self.expand_to_full(u_red);
        let n_red_nnz = self.k_rows.len();
        let k_vals_full = {
            let asm = self
                .assembler
                .as_mut()
                .expect("assembler presence checked above");

            // Frame consistency: u_full is in the rotating frame. Keep
            // assembler reference coords in the same frame at current θ.
            let coords_rotated = self
                .transforms
                .disps_to_inertial(&self.all_node_coords, self.theta);
            asm.update_node_coordinates(&coords_rotated);

            let (_, _, vals) = asm.assemble_kt_corotational(&u_full);
            vals
        };

        let k_red = crate::petsc::fsi::setup::apply_kg_coo_map(&self.k_coo_map, &k_vals_full, n_red_nnz);
        self.stepper
            .update_elastic_stiffness(&k_red)
            .map_err(FsiError::StepperError)?;

        log::info!(
            "RotorFsi: K_T rebuilt at corotational step {} (corotational={})",
            coro_step_counter,
            self.config.use_corotational_kt
        );
        Ok(())
    }

    // ── Main coupling loop ────────────────────────────────────────────────────

    /// Run the preCICE co-rotational FSI coupling loop.
    ///
    /// # Returns
    /// [`FsiResult`] containing the time history of displacement, velocity, and
    /// acceleration in the ROTATING frame (reduced DOF space).
    pub fn run(&mut self) -> Result<FsiResult, FsiError> {
        // ── Initialize preCICE ────────────────────────────────────────────────
        let mut participant = precice::Participant::new(
            &self.config.fsi.participant_name,
            &self.config.fsi.config_file,
            0,
            1,
        )?;

        let n_vertices = self.interface_coords_global.len() / self.mesh_dims.max(1);
        let mut vertex_ids = vec![0i32; n_vertices];
        participant.set_mesh_vertices(
            &self.config.fsi.coupling_mesh,
            &self.interface_coords_global,
            &mut vertex_ids,
        )?;

        // Optional GlobalSolidMesh for ω communication.
        let omega_vertex_ids: Option<Vec<i32>> =
            if let (Some(mesh), Some(coord)) = (
                &self.config.omega_mesh_name,
                &self.config.omega_vertex_coord,
            ) {
                let mut ids = vec![0i32; 1];
                participant.set_mesh_vertices(mesh, coord, &mut ids)?;
                Some(ids)
            } else {
                None
            };

        let restart_time = self.initial_state.as_ref().map(|state| state.t).unwrap_or(0.0);

        // Write initial omega before initialize() when preCICE requires it.
        // This satisfies the initialize="true" exchange for AngularVelocity.
        if participant.requires_initial_data()? {
            if let (Some(mesh), Some(wdata), Some(ids)) = (
                &self.config.omega_mesh_name,
                &self.config.omega_write_data,
                &omega_vertex_ids,
            ) {
                let omega_init = self.omega_provider.get(restart_time).0;
                participant.write_data(mesh, wdata, ids, &[omega_init])?;
            }
        }

        participant.initialize()?;
        let mut dt = participant.get_max_time_step_size()?;

        // ── Optional restart ──────────────────────────────────────────────────
        if let Some(ref state) = self.initial_state {
            self.stepper.set_state(&state.u, &state.v, &state.a, state.t);
            self.theta = self.initial_theta;
        }

        // ── Initial K_SP build ────────────────────────────────────────────────
        let omega_init = self.omega_provider.get(restart_time).0;
        self.apply_ksp(omega_init)?;

        let n_dofs = self.stepper.n_dofs();
        let n_data = n_vertices * self.mesh_dims;

        // Checkpoint storage.
        let mut newmark_cp: Option<NewmarkCheckpoint> = None;
        let mut omega_cp: Option<OmegaCheckpoint> = None;
        let mut theta_cp: f64 = self.theta;

        let mut result = FsiResult::default();
        let mut prof = WindowProfiler::new("RotorFsi");

        // ── Coupling loop ─────────────────────────────────────────────────────
        while participant.is_coupling_ongoing()? {
            // ── Save checkpoint ───────────────────────────────────────────────
            if participant.requires_writing_checkpoint()? {
                prof.start_window();
                let t_cp = Instant::now();
                newmark_cp = Some(self.stepper.checkpoint());
                omega_cp = Some(self.omega_provider.checkpoint());
                theta_cp = self.theta;
                prof.record_since("checkpoint_save", t_cp);
            }

            let t = self.stepper.current_time();
            let window_kin = self.omega_provider.kinematics_over_step(t, dt);
            let omega_step = window_kin.omega_step;
            let alpha_step = window_kin.alpha_step;
            let theta_target = self.theta + window_kin.theta_increment;

            // ── Read forces from preCICE (global frame) ───────────────────────
            let mut forces_global = vec![0.0f64; n_data];
            let t_io = Instant::now();
            participant.read_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.read_data,
                &vertex_ids,
                dt,
                &mut forces_global,
            )?;
            prof.record_since("precice_read", t_io);

            // ── Transform aero forces to rotating frame ───────────────────────
            let t_tr = Instant::now();
            let mut forces_local = self
                .transforms
                .forces_to_rotating(&forces_global, theta_target);
            prof.record_since("force_transform", t_tr);

            // ── Force pre-processing (ramp + cap in rotating frame) ───────────
            let t_fp = Instant::now();
            apply_ramp(&mut forces_local, t, self.config.fsi.ramp_time);
            if let Some(max_f) = self.config.fsi.force_max {
                apply_cap(&mut forces_local, max_f, self.mesh_dims);
            }
            prof.record_since("force_preproc", t_fp);

            // ── Assemble global reduced force vector ──────────────────────────
            let t_inertial = Instant::now();
            let mut f_red = vec![0.0f64; n_dofs];

            // Scatter aero forces at interface DOFs.
            for (li, &rdof) in self.interface_dofs.iter().enumerate() {
                if li < forces_local.len() && rdof < n_dofs {
                    f_red[rdof] += forces_local[li];
                }
            }

            // Inertial forces (all nodes, rotating frame).
            // Centrifugal force coordinate selection depends on include_ksp:
            //   include_ksp=true  → evaluate at X₀ (undeformed).
            //     K_SP·u already accounts for the displacement-dependent
            //     linearisation on the LHS; using X₀+u here would double-count
            //     the O(ω²|u|) term.
            //   include_ksp=false → evaluate at X₀+u (deformed).
            //     No LHS counterpart exists; the full nonlinear correction must
            //     appear in the force vector.
            //
            // Fix #5: zero the per-iteration inertial scratch buffer before
            // accumulating contributions from this sub-iteration.
            for v in self.inertial_scratch.iter_mut() { *v = 0.0; }

            if self.config.include_centrifugal {
                let f_cf = if self.config.include_ksp {
                    // K_SP active: evaluate centrifugal at undeformed coords X₀.
                    compute_centrifugal_force(
                        &self.all_node_coords,
                        &self.all_node_masses,
                        &self.transforms.axis,
                        &self.transforms.center,
                        omega_step,
                    )
                } else {
                    // K_SP inactive: evaluate centrifugal at deformed coords X₀+u.
                    let u_full = self.expand_to_full(self.stepper.current_u());
                    let deformed_coords: Vec<f64> = self
                        .all_node_coords
                        .iter()
                        .enumerate()
                        .map(|(k, &x0)| {
                            let node = k / 3;
                            let comp = k % 3;
                            let gdof = node * self.config.dofs_per_node + comp;
                            if gdof < u_full.len() {
                                x0 + u_full[gdof]
                            } else {
                                x0
                            }
                        })
                        .collect();
                    compute_centrifugal_force(
                        &deformed_coords,
                        &self.all_node_masses,
                        &self.transforms.axis,
                        &self.transforms.center,
                        omega_step,
                    )
                };
                // Fix #5: accumulate centrifugal into scratch.
                for (dst, &src) in self.inertial_scratch.iter_mut().zip(f_cf.iter()) {
                    *dst += src;
                }
                self.scatter_node_forces(&f_cf, &mut f_red);
            }

            if self.config.include_euler && alpha_step.abs() > 1e-14 {
                // Euler evaluated at deformed coords X₀ + u.
                let u_full = self.expand_to_full(self.stepper.current_u());
                let deformed: Vec<f64> = self
                    .all_node_coords
                    .iter()
                    .enumerate()
                    .map(|(k, &x0)| {
                        let node = k / 3;
                        let comp = k % 3;
                        let gdof = node * self.config.dofs_per_node + comp;
                        if gdof < u_full.len() {
                            x0 + u_full[gdof]
                        } else {
                            x0
                        }
                    })
                    .collect();
                let f_euler = compute_euler_force(
                    &deformed,
                    &self.all_node_masses,
                    &self.transforms.axis,
                    &self.transforms.center,
                    alpha_step,
                );
                // Fix #5: accumulate Euler into scratch.
                for (dst, &src) in self.inertial_scratch.iter_mut().zip(f_euler.iter()) {
                    *dst += src;
                }
                self.scatter_node_forces(&f_euler, &mut f_red);
            }

            if self.config.gravity_active() {
                let g_rot =
                    self.transforms.gravity_to_rotating(self.config.gravity, theta_target);
                let f_g = compute_gravity_force(&self.all_node_masses, &g_rot);
                self.scatter_node_forces(&f_g, &mut f_red);
            }
            prof.record_since("inertial_forces", t_inertial);

            // ── Advance structural state ──────────────────────────────────────
            // step() returns only the updated time; the structural state (u,v,a)
            // is stored in self.stepper and accessed via current_u/v/a() below.
            let t_step = Instant::now();
            let step_t = self.stepper.step(&f_red, dt)?.t;
            prof.record_since("newmark_step", t_step);

            // ── Gather interface displacements (rotating frame) ───────────────
            let t_gather = Instant::now();
            let disp_iface_local: Vec<f64> = self
                .interface_dofs
                .iter()
                .map(|&rdof| {
                    if rdof < self.stepper.n_dofs() {
                        self.stepper.current_u()[rdof]
                    } else {
                        0.0
                    }
                })
                .collect();

            // ── Transform displacements to global frame ───────────────────────
            let disp_iface_global = self
                .transforms
                .disps_to_inertial(&disp_iface_local, theta_target);
            prof.record_since("disp_gather_transform", t_gather);

            // ── Write displacements to preCICE ────────────────────────────────
            let t_wr = Instant::now();
            participant.write_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.write_data,
                &vertex_ids,
                &disp_iface_global,
            )?;

            // ── Write ω to GlobalSolidMesh (if used) ──────────────────────────
            if let (Some(ref ids), Some(ref mesh), Some(ref wdata)) = (
                &omega_vertex_ids,
                &self.config.omega_mesh_name,
                &self.config.omega_write_data,
            ) {
                participant.write_data(mesh, wdata, ids, &[omega_step])?;
            }
            prof.record_since("precice_write", t_wr);

            let t_adv = Instant::now();
            participant.advance(dt)?;
            prof.record_since("precice_advance", t_adv);

            // ── Implicit coupling: restore or commit ──────────────────────────
            if participant.requires_reading_checkpoint()? {
                let t_rs = Instant::now();
                match (newmark_cp.as_ref(), omega_cp.as_ref()) {
                    (Some(ncp), Some(ocp)) => {
                        self.stepper.restore(ncp);
                        self.omega_provider.restore(ocp);
                        self.theta = theta_cp;
                    }
                    _ => {
                        return Err(FsiError::PreciceError(
                            "preCICE requires checkpoint restore but no checkpoint was saved"
                                .to_string(),
                        ))
                    }
                }
                prof.record_since("checkpoint_restore", t_rs);
            } else {
                // ── Converged time window ─────────────────────────────────────
                // 1. Commit the converged target angle for this time window.
                self.theta = theta_target;

                // 2. Compute aerodynamic torque (for ω dynamics and logging).
                //    Use interface reference coords + displacements in the ROTATING frame.
                let n_iface = self.interface_dofs.len() / self.mesh_dims.max(1);
                let iface_coords_rot: Vec<f64> = self
                    .interface_coords_global
                    .iter()
                    .take(n_iface * 3)
                    .copied()
                    .collect();
                let iface_disps_rot: Vec<f64> = {
                    let mut buf = vec![0.0f64; n_iface * 3];
                    for k in 0..n_iface * 3 {
                        let rdof = self.interface_dofs.get(k).copied().unwrap_or(0);
                        buf[k] = if rdof < self.stepper.n_dofs() {
                            self.stepper.current_u()[rdof]
                        } else {
                            0.0
                        };
                    }
                    buf
                };
                let (_, tau_aero) = compute_torque(
                    &iface_coords_rot,
                    &iface_disps_rot,
                    &forces_local[..n_iface * 3],
                    &self.transforms.axis,
                    &self.transforms.center,
                );

                // 3. Gravity torque in the rotating frame.
                //    τ_gravity = torque from gravity body forces about the rotation axis.
                //    Needed so that ω dynamics use τ_driving = τ_aero + τ_gravity
                //    (gravitational sag changes the driving torque; centrifugal/Coriolis/Euler
                //    are fictitious forces and must NOT be included).
                let tau_gravity = if self.config.gravity_active() {
                    let g_rot = self.transforms.gravity_to_rotating(self.config.gravity, self.theta);
                    let f_grav = compute_gravity_force(&self.all_node_masses, &g_rot);
                    let u_full = self.expand_to_full(self.stepper.current_u());
                    let u_nodes = self.extract_node_translations(&u_full);
                    // Torque about axis at the deformed nodal positions X0 + u.
                    let (_, tau_g) = compute_torque(
                        &self.all_node_coords,
                        &u_nodes,
                        &f_grav,
                        &self.transforms.axis,
                        &self.transforms.center,
                    );
                    tau_g
                } else {
                    0.0
                };

                // 4. Update ω from driving torque = τ_aero + τ_gravity.
                //    (excludes centrifugal, Coriolis, Euler — those are fictitious forces)
                self.omega_provider.update_from_torque(tau_aero + tau_gravity, dt, t + dt);

                let (omega_new, _) = self.omega_provider.get(t + dt);
                let time_step = result.times.len() + 1; // 1-based, before push

                // 5. Rebuild K_SP when ω changes significantly (unified relative predicate).
                // ksp_omega_threshold² seeds the near-zero eps guard.
                let steps_since_ksp = time_step.saturating_sub(self.last_ksp_rebuild_step);
                let ksp_recently_rebuilt = steps_since_ksp < 10;
                let ksp_eps_sq = self.config.ksp_omega_threshold * self.config.ksp_omega_threshold;
                if omega_changed_significantly(
                    omega_new,
                    self.omega_sq_at_last_ksp,
                    self.config.omega_rebuild_rel_high,
                    self.config.omega_rebuild_rel_low,
                    ksp_recently_rebuilt,
                    ksp_eps_sq,
                ) {
                    self.last_ksp_rebuild_step = time_step;
                    let t_ksp = Instant::now();
                    self.apply_ksp(omega_new)?;
                    prof.record_since("ksp_rebuild", t_ksp);
                }

                // 6. Rebuild K_G from centrifugal prestress on the configured cadence.
                //    Snapshot the current reduced displacement so the K_G(u) trigger and
                //    deformed-coords assembly path see a consistent state and we don't
                //    conflict with the &mut self borrow inside update_kg_if_needed.
                let u_red_snapshot = self.stepper.current_u().to_vec();
                let t_kg = Instant::now();
                self.update_kg_if_needed(omega_new, &u_red_snapshot, time_step)?;
                prof.record_since("kg_update_if_needed", t_kg);

                // 7. Optional K_T(u) rebuild in the structural frame.
                let t_kt = Instant::now();
                self.coro_step_counter = self.coro_step_counter.saturating_add(1);
                self.update_kt_if_needed(&u_red_snapshot, self.coro_step_counter)?;
                prof.record_since("kt_update_if_needed", t_kt);

                // 8. Overwrite final state (no history accumulation — O(n_dofs) RAM).
                result.u_final = self.stepper.current_u().to_vec();
                result.v_final = self.stepper.current_v().to_vec();
                result.a_final = self.stepper.current_a().to_vec();
                result.times.push(step_t);

                // 9. Performance coefficients.
                let thrust = compute_thrust(&forces_global, &self.transforms.axis);
                let power_aero = tau_aero * omega_step;
                let perf = compute_performance_coefficients(
                    thrust,
                    power_aero,
                    tau_aero,
                    omega_step,
                    self.config.rotor_radius,
                    self.config.fluid_density,
                    self.config.flow_velocity,
                );

                let (omega_state, alpha_state) = self.omega_provider.get(t + dt);

                log::info!(
                    "RotorFsi step {time_step}: t={:.4} θ={:.4}rad ω_window={omega_step:.4}rad/s ω_state={omega_state:.4}rad/s \
                     τ_aero={tau_aero:.3e}N·m τ_grav={tau_gravity:.3e}N·m  Ct={:.4} Cp={:.4} TSR={:.3}",
                    step_t, self.theta, perf.ct, perf.cp, perf.tsr
                );

                // 10. Per-step callback.
                if let Some(ref cb) = self.step_callback {
                    let t_cb = Instant::now();
                    let force_mag = forces_global.iter().map(|x| x * x).sum::<f64>().sqrt();
                    // Fix #5: pass the last-applied inertial forces (centrifugal+Euler)
                    // to the callback so Python _step_cb can consume them directly
                    // without recomputing.
                    let inertial_snap = &self.inertial_scratch;
                    cb(
                        step_t,
                        time_step,
                        dt,
                        self.stepper.current_u(),
                        self.stepper.current_v(),
                        self.stepper.current_a(),
                        force_mag,
                        &forces_global,
                        omega_state,
                        alpha_state,
                        self.theta,
                        tau_aero,
                        perf,
                        inertial_snap,
                    )?;
                    prof.record_since("callback", t_cb);
                }

                dt = participant.get_max_time_step_size()?;

                // ── Per-window profiling summary (debug log) ──────────────────
                prof.log_summary(step_t, time_step);
            }
        }

        participant.finalize()?;
        Ok(result)
    }

    /// Number of free DOFs in the structural system.
    pub fn n_dofs(&self) -> usize {
        self.stepper.n_dofs()
    }

    /// Number of interface nodes registered with the coupling mesh.
    pub fn n_interface_nodes(&self) -> usize {
        self.interface_coords_global.len() / self.mesh_dims.max(1)
    }

    /// Current cumulative rotation angle θ [rad].
    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Current angular velocity ω [rad/s] from the provider.
    pub fn omega(&self) -> f64 {
        self.omega_provider.get(self.stepper.current_time()).0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_config(kg_update_interval: usize) -> RotorFsiConfig {
        RotorFsiConfig {
            fsi: FsiConfig {
                participant_name: "Solid".to_string(),
                config_file: "precice-config.xml".to_string(),
                coupling_mesh: "Solid-Mesh".to_string(),
                write_data: "Displacement".to_string(),
                read_data: "Force".to_string(),
                ramp_time: 0.0,
                force_max: None,
            },
            rotation_axis: [0.0, 0.0, 1.0],
            rotation_center: [0.0, 0.0, 0.0],
            gravity: [0.0, 0.0, 0.0],
            include_centrifugal: true,
            include_coriolis: false,
            include_euler: false,
            include_kg: true,
            kg_update_interval,
            use_corotational_kt: false,
            kt_coro_update_freq: 1,
            include_ksp: true,
            ksp_omega_threshold: 1e-4,
            omega_rebuild_rel_high: 0.005,
            omega_rebuild_rel_low: 0.003,
            kg_use_deformed_coords: false,
            kg_deflection_rebuild_rel_high: 0.01,
            kg_deflection_rebuild_rel_low: 0.005,
            dofs_per_node: 6,
            fluid_density: 1.225,
            flow_velocity: 10.0,
            rotor_radius: 1.0,
            omega_mesh_name: None,
            omega_write_data: None,
            omega_vertex_coord: None,
        }
    }

    #[test]
    fn kg_update_interval_zero_means_every_step() {
        let cfg = dummy_config(0);
        assert_eq!(cfg.normalized_kg_update_interval(), 1);
        assert!(cfg.should_update_kg_on_step(1));
        assert!(cfg.should_update_kg_on_step(2));
        assert!(cfg.should_update_kg_on_step(3));
    }

    #[test]
    fn kg_update_interval_updates_on_multiples_only() {
        let cfg = dummy_config(3);
        assert_eq!(cfg.normalized_kg_update_interval(), 3);
        assert!(!cfg.should_update_kg_on_step(1));
        assert!(!cfg.should_update_kg_on_step(2));
        assert!(cfg.should_update_kg_on_step(3));
        assert!(!cfg.should_update_kg_on_step(4));
        assert!(!cfg.should_update_kg_on_step(5));
        assert!(cfg.should_update_kg_on_step(6));
    }

    #[test]
    fn kt_coro_update_freq_zero_means_every_step() {
        let mut cfg = dummy_config(1);
        cfg.kt_coro_update_freq = 0;
        assert_eq!(cfg.normalized_kt_coro_update_freq(), 1);
        assert!(cfg.should_update_coro_frame_on_step(1));
        assert!(cfg.should_update_coro_frame_on_step(2));
        assert!(cfg.should_update_coro_frame_on_step(3));
    }

    #[test]
    fn kt_coro_update_freq_updates_on_multiples_only() {
        let mut cfg = dummy_config(1);
        cfg.kt_coro_update_freq = 5;
        assert_eq!(cfg.normalized_kt_coro_update_freq(), 5);
        assert!(!cfg.should_update_coro_frame_on_step(1));
        assert!(!cfg.should_update_coro_frame_on_step(4));
        assert!(cfg.should_update_coro_frame_on_step(5));
        assert!(!cfg.should_update_coro_frame_on_step(9));
        assert!(cfg.should_update_coro_frame_on_step(10));
    }

    #[test]
    fn corotational_kt_disabled_keeps_baseline_path() {
        let mut cfg = dummy_config(1);
        cfg.use_corotational_kt = false;
        for step in 1..=10 {
            assert!(!cfg.should_update_corotational_kt(step));
        }
    }

    // ── Tests for omega_changed_significantly predicate (task 3.1 / 3.2) ──────
    //
    // eps = (1e-4)² = 1e-8 mirrors the default ksp_omega_threshold² used in
    // production call sites, per task 3.2: ksp_omega_threshold is rerouted as
    // the eps guard for the near-zero ω² branch.

    /// Stable ω (change well below both thresholds) → skip.
    #[test]
    fn omega_pred_stable_omega_skips() {
        let eps = 1e-8_f64; // (1e-4)^2
        let omega_base = 100.0f64;
        let omega_sq_last = omega_base * omega_base;
        // 0.1% change — below both 0.5% and 0.3%
        let omega_new = omega_base * (1.0 + 0.001);
        // Not recently rebuilt
        assert!(!omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, false, eps));
        // Recently rebuilt
        assert!(!omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, true, eps));
    }

    /// Large jump above high threshold → rebuild regardless of hysteresis state.
    #[test]
    fn omega_pred_large_jump_rebuilds() {
        let eps = 1e-8_f64;
        let omega_base = 100.0f64;
        let omega_sq_last = omega_base * omega_base;
        // 0.8% change — above both 0.5% and 0.3%
        let omega_new = omega_base * (1.0 + 0.004);
        // Not recently rebuilt: 0.8% > 0.3% → rebuild
        assert!(omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, false, eps));
        // Recently rebuilt: 0.8% > 0.5% → rebuild
        assert!(omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, true, eps));
    }

    /// ω → 0: both ω² values below eps → always rebuild (avoid divide-by-zero).
    #[test]
    fn omega_pred_omega_to_zero_always_rebuilds() {
        // Use eps = 1e-8; values whose squares are both 0 → below eps → rebuild.
        let eps = 1e-8_f64;
        // max(0, 0) = 0 < eps → rebuild
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, false, eps));
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, true, eps));
    }

    /// First call (omega_sq_at_last = NEG_INFINITY) → always rebuild.
    #[test]
    fn omega_pred_first_call_always_rebuilds() {
        let eps = 1e-8_f64;
        let omega_new = 50.0;
        assert!(omega_changed_significantly(
            omega_new,
            f64::NEG_INFINITY,
            0.005,
            0.003,
            false,
            eps,
        ));
        assert!(omega_changed_significantly(
            omega_new,
            f64::NEG_INFINITY,
            0.005,
            0.003,
            true,
            eps,
        ));
        // Also for zero omega on first call
        assert!(omega_changed_significantly(
            0.0,
            f64::NEG_INFINITY,
            0.005,
            0.003,
            false,
            eps,
        ));
    }

    /// Hysteresis: change between skip and rebuild thresholds behaves correctly.
    ///
    /// Change of 0.4% (between 0.3% skip and 0.5% rebuild thresholds):
    /// - `currently_rebuilt = true`  (high bar = 0.5%) → 0.4% < 0.5% → skip
    /// - `currently_rebuilt = false` (low  bar = 0.3%) → 0.4% > 0.3% → rebuild
    #[test]
    fn omega_pred_hysteresis_between_thresholds() {
        let eps = 1e-8_f64;
        let omega_base = 100.0f64;
        let omega_sq_last = omega_base * omega_base;
        // 0.4% ω change → 0.8% ω² change (approximation: Δω²/ω² ≈ 2Δω/ω = 0.8%)
        // Use exact calculation to be safe.
        let omega_new = omega_base * (1.0 + 0.002);
        let omega_sq_new = omega_new * omega_new;
        let rel = (omega_sq_new - omega_sq_last).abs() / omega_sq_new.max(omega_sq_last);
        assert!(rel > 0.003 && rel < 0.005, "rel={rel} must be in (0.003, 0.005)");

        // recently rebuilt → use high threshold (0.005) → change < 0.5% → skip
        assert!(!omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, true, eps));
        // not recently rebuilt → use low threshold (0.003) → change > 0.3% → rebuild
        assert!(omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, false, eps));
    }

    /// Both omega² values exactly at the eps guard → rebuild.
    #[test]
    fn omega_pred_eps_boundary_rebuilds() {
        // max(ω², ω_last²) = 0.0 < eps → below guard → rebuild
        let eps = 1e-8_f64;
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, false, eps));
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, true, eps));
    }

    #[test]
    fn corotational_kt_coords_follow_theta_for_zero_displacement() {
        use std::f64::consts::FRAC_PI_2;

        let transforms = RotorTransforms::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let coords_ref = vec![1.0, 0.0, 0.0];

        // This is exactly the transform used in update_kt_if_needed.
        let coords_rotated = transforms.disps_to_inertial(&coords_ref, FRAC_PI_2);

        // (1,0,0) rotated 90° around +Z => (0,1,0)
        assert!(coords_rotated[0].abs() < 1e-12, "x={}", coords_rotated[0]);
        assert!((coords_rotated[1] - 1.0).abs() < 1e-12, "y={}", coords_rotated[1]);
        assert!(coords_rotated[2].abs() < 1e-12, "z={}", coords_rotated[2]);
    }
}
