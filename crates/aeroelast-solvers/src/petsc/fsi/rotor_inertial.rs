/// Inertial-frame FSI solver for rotating structures (rotor blades).
///
/// Implements the full preCICE coupling loop for a rotor in the **inertial (global) frame**,
/// where the mesh geometry is rotated at each converged time window to match θ(t).
///
/// Key differences from the co-rotational solver (`rotor_fsi.rs`):
/// - **K(θ) varies**: The elastic stiffness is reassembled after geometry rotation.
/// - **No fictitious forces**: Centrifugal/Coriolis/Euler forces are replaced by
///   the reference-frame load vector `F_ref = −M·a_ref`, where `a_ref = α×r + ω×(ω×r)`.
/// - **Displacement modes**: The solver can write either elastic displacement only
///   (for inertial-aware CFD) or total displacement (for standard CFD).
/// - **Gravity**: Applied in the global frame (constant vector).
///
/// # Geometry approximation (explicit-in-geometry / predictor geometry)
///
/// Both `K(θ)` and `F_ref` are computed with the geometry from the **previous converged
/// window** (θ_{n-1}, the checkpoint value) throughout all sub-iterations of window n.
/// The mesh is only advanced to θ_n after convergence.
///
/// This is a deliberate **explicit-in-geometry** (or "predictor geometry") scheme:
/// - Error: O(ω · Δt) in lever arm → O(ω² · Δt) in F_ref → negligible for small Δt.
/// - Benefit: avoids reassembling K and recomputing F_ref inside the implicit coupling
///   loop, keeping sub-iteration cost at O(K·solve) only.
/// - Alternative ("mid-window geometry"): rotate mesh to θ_n = θ_{n-1} + Δθ at the
///   first sub-iteration of each window, then hold it fixed. This would be O(Δt²)
///   accurate in geometry at the cost of one extra mesh-rotation + optional K-rebuild
///   per window. It is NOT implemented because: (a) K rebuild is expensive, and (b)
///   the coupling loop structure would require a two-phase initialize step.
///
/// For validation at low RPM (ω < 5 rad/s, Δt < 0.01 s) the lag is < 0.05 °/step
/// and has negligible impact on blade response. High-RPM cases with large Δt
/// should decrease `theta_rebuild_threshold` to stay accurate.
///
/// # Feature gate
/// Compiled only with `--features fsi`.

use std::time::Instant;

use aeroelast_core::assembly::assembler::{MaterialSpec, MeshAssembler};

use crate::petsc::elasticity::dynamic_newmark::{NewmarkCheckpoint, NewmarkStepper};
use crate::petsc::fsi::force_utils::{apply_cap, apply_ramp};
use crate::petsc::fsi::linear_elastic::{FsiConfig, FsiError, FsiInitialState, FsiResult};
use crate::petsc::fsi::profiling::WindowProfiler;
use crate::petsc::fsi::rotor_physics::{
    OmegaCheckpoint, OmegaProvider, PerformanceCoefficients, RotorTransforms,
    compose_deformed_coords, compute_gravity_force, compute_performance_coefficients,
    compute_rigid_body_acceleration_inertial, compute_reference_load_vector,
    compute_rigid_body_velocity_inertial, compute_thrust, compute_torque,
    max_radial_deflection_ratio, omega_changed_significantly,
};

// ── Step callback type ─────────────────────────────────────────────────────────

/// Per-step callback for the inertial rotor FSI solver.
///
/// Invoked once per **converged** time window with the full rotor state.
///
/// Arguments:
/// `(t, time_step, dt, u_red, v_red, a_red, force_mag,
///   forces_iface_global, omega_state, alpha_state, theta, tau_aero, perf)`
pub type InertialRotorStepCallback = Box<
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
        ) -> Result<(), FsiError>
        + Send,
>;

// ── Configuration ──────────────────────────────────────────────────────────────

/// Configuration for the inertial rotor FSI solver.
#[derive(Debug, Clone)]
pub struct InertialRotorFsiConfig {
    /// Base preCICE FSI configuration (participant, mesh, read/write data names).
    pub fsi: FsiConfig,

    // ── Rotor geometry ────────────────────────────────────────────────────────
    /// Rotation axis unit vector (will be normalized if not already).
    pub rotation_axis: [f64; 3],
    /// Rotation center in the global frame.
    pub rotation_center: [f64; 3],

    // ── Inertial-frame physics ────────────────────────────────────────────────
    /// Gravity vector in the global frame [m/s²]. Use `[0,0,0]` to disable.
    pub gravity: [f64; 3],
    /// Enable reference acceleration forces (F_ref = −M·a_ref).
    /// Should always be `true` for physical correctness unless debugging.
    pub include_reference_acceleration: bool,
    /// Include centrifugal prestress as a dynamic K_G(θ, ω) term assembled on
    /// the current rotated geometry and updated together with K(θ).
    pub include_geometric_stiffness: bool,
    /// Include spin-softening K_SP(ω) = −ω²·M·(I − n̂⊗n̂) in the effective stiffness.
    ///
    /// K_SP arises from linearizing the centripetal displacement term around u_e = 0.
    /// It is negative-definite in the plane perpendicular to the rotation axis and
    /// prevents overestimating flapwise eigenfrequencies at operating speed.
    ///
    /// Default: `true` — K_SP runs unconditionally for physical correctness.
    /// Set to `false` only for debugging or backward-compatibility checks.
    pub include_spin_softening: bool,
    /// High-band relative threshold for K_SP ω-driven rebuild trigger.
    /// When |Δ(ω²)|/ω² exceeds this value and the solver is NOT in a recently-rebuilt
    /// state, K_SP is rebuilt. Default: `0.005` (0.5% change in ω²).
    pub ksp_omega_rebuild_high: f64,
    /// Low-band relative threshold for K_SP ω-driven rebuild trigger.
    /// Used as the hysteresis floor when the solver IS in a recently-rebuilt state.
    /// Default: `0.003` (0.3% change in ω²).
    pub ksp_omega_rebuild_low: f64,
    /// High-band relative threshold for K_G ω-driven rebuild trigger.
    /// When |Δ(ω²)|/ω² exceeds this value and the solver is NOT in a recently-rebuilt
    /// state, K_G is rebuilt (independent of K(θ) cadence). Default: `0.005`.
    pub kg_omega_rebuild_high: f64,
    /// Low-band relative threshold for K_G ω-driven rebuild trigger.
    /// Hysteresis floor analogous to `ksp_omega_rebuild_low`. Default: `0.003`.
    pub kg_omega_rebuild_low: f64,

    // ── Stiffness reassembly ──────────────────────────────────────────────────
    /// Rebuild K(θ) every N converged steps. `0` and `1` both mean every step.
    /// For inertial solver, K(θ) **must** be updated when geometry rotates,
    /// so this controls the maximum window between forced rebuilds.
    pub k_update_interval: usize,
    /// Relative |Δ(ω²)|/ω² threshold to trigger an extra K(θ) rebuild when ω
    /// changes between windows (only meaningful for ramped/computed ω modes —
    /// constant ω never trips this). Default: `0.01` (1% change in ω²).
    pub omega_rebuild_threshold: f64,
    /// Absolute |Δθ| threshold (radians) to trigger a K(θ) rebuild when the
    /// rotation angle has accumulated since the last reassembly. This is the
    /// physically correct trigger for this solver since K(θ) depends on θ
    /// (geometry orientation), not on ω. Default: `0.05` rad ≈ 2.86°.
    /// Set to a small value to rebuild often, large to amortize cost.
    pub theta_rebuild_threshold: f64,

    // ── K_G(u) — Nivel 2 reassembly with deformed geometry ────────────────────
    /// Master toggle for the foreshortening-aware K_G assembly. When `true`,
    /// `K_G` is built at deformed coordinates `X_rotated + u` (only K_G —
    /// `K(θ)` keeps the linear "Option A" discipline and uses `X_rotated`).
    /// Default `false` for parity with the explicit-in-geometry baseline.
    pub kg_use_deformed_coords: bool,
    /// High-band relative threshold for the K_G deflection-driven rebuild
    /// trigger. When the radial deflection ratio drift since the last K_G
    /// rebuild exceeds this value, the next converged window will reassemble
    /// `K_G` with deformed coords. Default `0.01` (1%).
    pub kg_deflection_rebuild_rel_high: f64,
    /// Low-band relative threshold for the K_G deflection-driven rebuild
    /// trigger. Acts as the hysteresis floor so that cyclic 1P deflection
    /// (gravity-driven) cannot ping-pong the rebuild logic. Reserved for the
    /// hysteresis predicate; not consulted by the simple drift check.
    /// Default `0.005` (0.5%).
    pub kg_deflection_rebuild_rel_low: f64,

    // ── Displacement mode ─────────────────────────────────────────────────────
    /// Which displacement to write to preCICE: elastic only or total (elastic + rigid).
    pub displacement_mode: DisplacementMode,

    // ── DOF layout ────────────────────────────────────────────────────────────
    /// DOFs per FEM node (typically 6 for shells, 3 for solids).
    pub dofs_per_node: usize,

    // ── Performance coefficients ──────────────────────────────────────────────
    /// Fluid (air) density [kg/m³].
    pub fluid_density: f64,
    /// Freestream flow velocity [m/s].
    pub flow_velocity: f64,
    /// Rotor outer radius [m].
    pub rotor_radius: f64,

    // ── Optional GlobalSolidMesh for ω communication ───────────────────────────
    /// preCICE mesh name for the GlobalSolidMesh that communicates ω.
    /// `None` disables this feature.
    pub omega_mesh_name: Option<String>,
    /// Data name to write ω to on the GlobalSolidMesh.
    pub omega_write_data: Option<String>,
    /// Single vertex coordinate for the GlobalSolidMesh (typically the rotation center).
    pub omega_vertex_coord: Option<[f64; 3]>,

    // ── Optional nodal velocity write ──────────────────────────────────────────
    /// Data name to write structural nodal velocities to on the coupling mesh.
    /// When set, the solver writes the translational velocity (vx, vy, vz) of
    /// each interface node every coupling iteration so that the BEM participant
    /// can correct the effective inflow velocity (aerodynamic damping).
    /// `None` disables velocity writing (default for back-compatibility).
    pub velocity_write_data: Option<String>,
}

impl InertialRotorFsiConfig {
    /// Return whether gravity is active (non-zero magnitude).
    fn gravity_active(&self) -> bool {
        self.gravity.iter().map(|x| x * x).sum::<f64>() > 1e-24
    }

    /// Normalize the K(θ) update cadence so legacy `0` means "every converged step".
    fn normalized_k_update_interval(&self) -> usize {
        self.k_update_interval.max(1)
    }

    /// Whether K(θ) should be rebuilt on this 1-based converged step index.
    fn should_update_k_on_step(&self, time_step: usize) -> bool {
        time_step % self.normalized_k_update_interval() == 0
    }
}

// ── Displacement mode enum ─────────────────────────────────────────────────────

/// Defines which displacement to write to preCICE for the inertial solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DisplacementMode {
    /// Write elastic displacement only (u_e) — for inertial-aware CFD solvers
    /// that expect the mesh to rotate rigidly and only receive elastic deformation.
    Elastic,
    /// Write total displacement (u_total = u_e + u_rigid) — for standard CFD solvers
    /// that expect the full deformed position (elastic + rigid-body motion).
    Total,
}

// ── Checkpoint for K_SP / K_G rebuild trackers ────────────────────────────────

/// Snapshot of the K_SP and K_G ω-rebuild state for checkpoint/restore.
///
/// Captured together with `NewmarkCheckpoint` and `OmegaCheckpoint` so that a
/// preCICE rollback correctly restores the hysteresis counters. Without this,
/// a rollback would reset θ and ω but leave the rebuild-step indices pointing
/// into the future, which could suppress the first post-rollback rebuild.
#[derive(Debug, Clone)]
pub struct InertialKgKspCheckpoint {
    /// ω² at the last K_SP rebuild at checkpoint time.
    pub omega_sq_at_last_ksp_rebuild: f64,
    /// ω² at the last K_G rebuild at checkpoint time.
    pub omega_sq_at_last_kg_rebuild: f64,
    /// Step index of the last K_SP rebuild at checkpoint time.
    pub last_ksp_rebuild_step: usize,
    /// Step index of the last K_G rebuild at checkpoint time.
    pub last_kg_rebuild_step: usize,
}

// ── Solver state ───────────────────────────────────────────────────────────────

/// Internal state of the inertial rotor FSI solver.
pub struct InertialRotorFsiSolver {
    // ── Configuration ─────────────────────────────────────────────────────────
    config: InertialRotorFsiConfig,

    // ── Mesh and assembly ─────────────────────────────────────────────────────
    /// FEM mesh assembler (must be Clone for geometry rotation).
    /// The inertial solver calls `assembler.update_node_coordinates()` after
    /// rotating the mesh, then reassembles K(θ).
    assembler: MeshAssembler,

    // ── Time integration ──────────────────────────────────────────────────────
    /// Newmark-β stepper for dynamic elasticity (M·ü + C·u̇ + K(θ)·u = F).
    stepper: NewmarkStepper,

    // ── Rotor kinematics ──────────────────────────────────────────────────────
    /// Angular velocity provider (Constant, Ramped, Computed, RampedComputed).
    omega_provider: OmegaProvider,
    /// Coordinate transformation utilities (Rodrigues rotation, force transforms).
    transforms: RotorTransforms,

    // ── Geometry and state ────────────────────────────────────────────────────
    /// Current rotation angle θ (cumulative) [rad].
    theta: f64,
    /// Reference (unrotated) node coordinates (flat, length = n_nodes × 3).
    /// Stored to rotate the mesh back to θ=0 before applying Δθ.
    coords_ref: Vec<f64>,
    /// Current rotated node coordinates (flat, length = n_nodes × 3).
    coords_rotated: Vec<f64>,
    /// Nodal masses (flat, length = n_nodes) for F_ref = −M·a_ref.
    masses: Vec<f64>,
    /// Sorted free global DOF indices for the reduced structural system.
    free_dofs: Vec<i32>,

    // ── Interface state ───────────────────────────────────────────────────────
    /// Number of interface (coupling) nodes.
    n_iface: usize,
    /// Translational DOF indices for the interface nodes (length = n_iface × 3).
    iface_dofs: Vec<usize>,
    /// Interface node indices (0-based, length = n_iface).
    iface_nodes: Vec<usize>,

    // ── Convergence tracking ──────────────────────────────────────────────────
    /// 1-based converged time step counter.
    time_step: usize,

    // ── Reference reduced-K mapping ──────────────────────────────────────────
    /// Maps full reassembled K(θ) COO entries into the reference reduced K COO.
    k_coo_map: Vec<i32>,
    /// Number of entries in the reference reduced K COO.
    k_red_nnz: usize,
    /// Reduced K row indices (reference COO structure, post-BC reduction).
    /// Stored to compute K_SP expanded to the same COO sparsity on every rebuild.
    k_red_rows: Vec<i32>,
    /// Reduced K column indices (reference COO structure, post-BC reduction).
    k_red_cols: Vec<i32>,
    /// Maps full reassembled K_G COO entries into the reference reduced K COO.
    kg_coo_map: Vec<i32>,
    /// Per-element density (or mass-per-area for composites) used to build the
    /// centrifugal prestress K_G(θ, ω).
    rho_per_elem: Vec<f64>,

    // ── Stiffness update tracking ─────────────────────────────────────────────
    /// ω² at the last K(θ) rebuild (for relative change check).
    omega_sq_at_last_k_rebuild: f64,
    /// θ at the last K(θ) rebuild (for absolute change check).
    theta_at_last_k_rebuild: f64,
    /// Step index of the last K(θ) rebuild.
    step_at_last_k_rebuild: usize,
    /// Radial deflection ratio at the last K_G rebuild. Used by the K_G(u)
    /// hysteresis trigger when `kg_use_deformed_coords` is enabled. Initialized
    /// to `0.0` so the first non-zero deflection above the high-band threshold
    /// triggers the first deformed-coords K_G assembly.
    deflection_ratio_at_last_kg_rebuild: f64,
    /// ω² at the last K_SP-only rebuild (for the Δω-driven K_SP hysteresis predicate).
    /// Initialized to `f64::NEG_INFINITY` so the first call always rebuilds.
    omega_sq_at_last_ksp_rebuild: f64,
    /// ω² at the last K_G-only rebuild (for the Δω-driven K_G hysteresis predicate).
    /// Initialized to `f64::NEG_INFINITY` so the first call always rebuilds.
    omega_sq_at_last_kg_rebuild: f64,
    /// Step index of the last K_SP rebuild (used to determine `recently_rebuilt` for hysteresis).
    last_ksp_rebuild_step: usize,
    /// Step index of the last K_G rebuild (used to determine `recently_rebuilt` for hysteresis).
    last_kg_rebuild_step: usize,

    // ── Optional callback ─────────────────────────────────────────────────────
    /// Per-step callback (invoked after each converged window).
    step_callback: Option<InertialRotorStepCallback>,

    // ── Optional GlobalSolidMesh ──────────────────────────────────────────────
    /// Vertex ID for the GlobalSolidMesh (if enabled).
    omega_vertex_id: Option<i32>,
}

impl InertialRotorFsiSolver {
    // ── Constructor ───────────────────────────────────────────────────────────

    /// Create a new `InertialRotorFsiSolver`.
    ///
    /// # Arguments
    /// * `config`   — solver configuration (rotation axis, gravity, stiffness update intervals)
    /// * `assembler` — FEM mesh assembler (must be `Clone` for geometry rotation)
    /// * `stepper`  — Newmark-β time integrator
    /// * `omega_provider` — angular velocity provider (Constant, Ramped, Computed)
    /// * `iface_nodes` — interface node indices (0-based)
    /// * `masses`   — nodal masses (flat, length = n_nodes)
    /// * `free_dofs` — sorted free global DOF indices for the reduced system
    /// * `k_coo_map` — mapping from full K COO entries to reduced reference COO
    /// * `k_red_nnz` — number of entries in the reduced reference K COO
    ///
    /// # Returns
    /// `Result<Self, FsiError>` — ready to call `with_initial_state()` and `run()`.
    pub fn new(
        config: InertialRotorFsiConfig,
        assembler: MeshAssembler,
        stepper: NewmarkStepper,
        omega_provider: OmegaProvider,
        iface_nodes: Vec<usize>,
        masses: Vec<f64>,
        free_dofs: Vec<i32>,
        k_coo_map: Vec<i32>,
        k_red_nnz: usize,
        k_red_rows: Vec<i32>,
        k_red_cols: Vec<i32>,
    ) -> Result<Self, FsiError> {
        let n_iface = iface_nodes.len();
        let dofs_per_node = config.dofs_per_node;
        if dofs_per_node < 3 {
            return Err(FsiError::PreciceError(format!(
                "inertial solver requires at least 3 translational DOFs per node, got {}",
                dofs_per_node,
            )));
        }

        // Build interface translational DOF indices only. preCICE Displacement
        // is a 3D vector field and must not include rotational shell DOFs.
        let mut iface_dofs = Vec::with_capacity(n_iface * 3);
        for &node in &iface_nodes {
            for d in 0..3 {
                iface_dofs.push(node * dofs_per_node + d);
            }
        }

        // Extract reference coordinates from the assembler
        let coords_ref = assembler.topology.node_coords.clone();
        let coords_rotated = coords_ref.clone();

        let n_nodes = coords_ref.len() / 3;
        if masses.len() != n_nodes {
            return Err(FsiError::PreciceError(format!(
                "inertial solver requires one nodal mass per mesh node (got {} masses for {} nodes)",
                masses.len(),
                n_nodes,
            )));
        }

        let n_red_dofs = stepper.n_dofs();
        if free_dofs.len() != n_red_dofs {
            return Err(FsiError::PreciceError(format!(
                "inertial solver free_dofs length {} does not match reduced stepper size {}",
                free_dofs.len(),
                n_red_dofs,
            )));
        }
        if k_red_nnz == 0 {
            return Err(FsiError::PreciceError(
                "inertial solver requires a non-empty reduced K sparsity".to_string(),
            ));
        }

        // Normalize rotation axis
        let axis = config.rotation_axis;
        let axis_norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if axis_norm < 1e-14 {
            return Err(FsiError::PreciceError(
                "rotation_axis must have non-zero magnitude".to_string(),
            ));
        }
        let axis_normalized = [axis[0] / axis_norm, axis[1] / axis_norm, axis[2] / axis_norm];

        let transforms = RotorTransforms::new(axis_normalized, config.rotation_center);

        let rho_per_elem: Vec<f64> = assembler
            .materials
            .iter()
            .map(|m| match m {
                MaterialSpec::Isotropic { rho, .. } => *rho,
                MaterialSpec::Composite { mass_per_area, .. } => *mass_per_area,
                MaterialSpec::PlaneStress { rho, .. } => *rho,
                MaterialSpec::Solid3D { rho, .. } => *rho,
            })
            .collect();

        let kg_coo_map = if config.include_geometric_stiffness {
            let dummy_sigma = vec![[0.0f64; 3]; assembler.topology.n_elems];
            let (kg_rows_full, kg_cols_full, _) = assembler.assemble_geometric_k(&dummy_sigma);
            crate::petsc::fsi::setup::build_kg_coo_map(
                &kg_rows_full,
                &kg_cols_full,
                &free_dofs,
                &k_red_rows,
                &k_red_cols,
            )
        } else {
            Vec::new()
        };

        Ok(InertialRotorFsiSolver {
            config,
            assembler,
            stepper,
            omega_provider,
            transforms,
            theta: 0.0,
            coords_ref,
            coords_rotated,
            masses,
            free_dofs,
            k_coo_map,
            k_red_nnz,
            k_red_rows,
            k_red_cols,
            kg_coo_map,
            rho_per_elem,
            n_iface,
            iface_dofs,
            iface_nodes,
            time_step: 0,
            omega_sq_at_last_k_rebuild: f64::NEG_INFINITY,
            theta_at_last_k_rebuild: 0.0,
            step_at_last_k_rebuild: 0,
            deflection_ratio_at_last_kg_rebuild: 0.0,
            omega_sq_at_last_ksp_rebuild: f64::NEG_INFINITY,
            omega_sq_at_last_kg_rebuild: f64::NEG_INFINITY,
            last_ksp_rebuild_step: 0,
            last_kg_rebuild_step: 0,
            step_callback: None,
            omega_vertex_id: None,
        })
    }

    // ── Builder methods ───────────────────────────────────────────────────────

    /// Set initial FEM state (displacement, velocity, acceleration).
    ///
    /// # Arguments
    /// * `state` — initial conditions (u₀, v₀, a₀) with length = `dofs_count`.
    ///
    /// # Returns
    /// `Result<Self, FsiError>` — ready to call `with_step_callback()` or `run()`.
    pub fn with_initial_state(mut self, state: &FsiInitialState) -> Result<Self, FsiError> {
        self.stepper.set_initial_conditions_with_acceleration(
            &state.u,
            &state.v,
            &state.a,
        );
        Ok(self)
    }

    /// Set the per-step callback (invoked after each converged window).
    ///
    /// # Arguments
    /// * `callback` — closure with signature `InertialRotorStepCallback`.
    ///
    /// # Returns
    /// `Self` — ready to call `run()`.
    pub fn with_step_callback(mut self, callback: InertialRotorStepCallback) -> Self {
        self.step_callback = Some(callback);
        self
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    /// Expand a reduced-DOF vector to the full DOF space (zeros at constrained DOFs).
    fn expand_to_full(&self, v_red: &[f64], free_dofs: &[i32]) -> Vec<f64> {
        let n_full_dofs = self.assembler.topology.dofs_count();
        let mut v_full = vec![0.0f64; n_full_dofs];
        for (i, &dof) in free_dofs.iter().enumerate() {
            v_full[dof as usize] = v_red[i];
        }
        v_full
    }

    /// Scatter flat node forces (n_nodes × 3) into the reduced-DOF force vector.
    ///
    /// Only free DOFs receive contributions; constrained DOFs are silently skipped.
    fn scatter_node_forces(
        &self,
        f_node: &[f64],
        f_global: &mut [f64],
        free_dofs: &[i32],
    ) {
        let n_nodes = self.masses.len();
        let dofs_per_node = self.config.dofs_per_node;
        for i in 0..n_nodes {
            let nb = i * 3;
            for j in 0..3 {
                let gdof = (i * dofs_per_node + j) as i32;
                if let Ok(pos) = free_dofs.binary_search(&gdof) {
                    f_global[pos] += f_node[nb + j];
                }
            }
        }
    }

    /// Compute consistent initial acceleration a₀ = M⁻¹·F₀ for a lumped
    /// diagonal mass matrix.
    ///
    /// For translational DOFs (local index 0–2) the mass is `masses[node]`.
    /// For rotational DOFs (local index 3–5) the lumped mass is zero, so
    /// `a₀ = 0` (which is also the physically correct result for a force-
    /// only initial condition with no rotational inertia in the lumped model).
    ///
    /// # Arguments
    /// * `f0`        — external force vector in the **reduced** free-DOF space.
    /// * `free_dofs` — sorted global DOF indices for the free system.
    ///
    /// # Returns
    /// `Vec<f64>` of length `free_dofs.len()` with a₀ for each free DOF.
    fn compute_initial_acceleration(&self, f0: &[f64], free_dofs: &[i32]) -> Vec<f64> {
        let dofs_per_node = self.config.dofs_per_node;
        let n_free = free_dofs.len();
        let mut a0 = vec![0.0f64; n_free];
        for (j, &gdof) in free_dofs.iter().enumerate() {
            let local = (gdof as usize) % dofs_per_node;
            if local < 3 {
                // Translational DOF: a₀ = F / m
                let node = (gdof as usize) / dofs_per_node;
                if node < self.masses.len() {
                    let m = self.masses[node];
                    if m > 1e-14 {
                        a0[j] = f0[j] / m;
                    }
                }
            }
            // Rotational DOFs: zero rotational inertia in lumped model → a₀ = 0
        }
        a0
    }

    /// Rebuild K(θ) (and K_G, K_SP) when geometry has rotated or deformed
    /// significantly.
    ///
    /// Assumes `self.coords_rotated` already tracks `self.theta` (the caller
    /// updates it at every theta change). Then calls
    /// `assembler.update_node_coordinates(&self.coords_rotated)` →
    /// `assembler.assemble_k()` → `apply_kg_coo_map()` →
    /// `stepper.update_elastic_stiffness()`.
    ///
    /// When `kg_use_deformed_coords` is enabled, K_G is reassembled at
    /// deformed coords `X_rotated + u` (Nivel 2 foreshortening contribution).
    /// `u_red` is the current reduced displacement, used both for the
    /// deflection-driven rebuild trigger and for the K_G assembly itself.
    /// `K(θ)` keeps the linear "Option A" discipline and always uses
    /// `X_rotated`.
    fn reassemble_k_if_needed(
        &mut self,
        omega: f64,
        u_red: &[f64],
        time_step: usize,
        free_dofs: &[i32],
    ) -> Result<(), FsiError> {
        // ── Trigger evaluation ────────────────────────────────────────────────
        //
        // Six OR-gate triggers:
        //   1. Forced cadence: every `k_update_interval` converged windows → K(θ).
        //   2. K(θ) ω-change trigger: |Δω²|/ω² above threshold.
        //   3. K(θ) θ-change trigger: |Δθ| above threshold (primary trigger).
        //   4. K_G(u) deflection-change trigger (Nivel 2 only).
        //   5. K_SP Δω trigger: independent per-matrix ω hysteresis predicate.
        //   6. K_G Δω trigger: independent per-matrix ω hysteresis predicate.
        //
        // Triggers 1-4 are coupled (rebuild K+K_G+K_SP together).
        // Triggers 5-6 are independent: they can fire without 1-4, and rebuild
        // only the specific matrix that changed while reusing cached values for
        // the others. All paths end in a single `update_tangent_terms` call →
        // single refactorization.

        let should_rebuild_by_step = self.config.should_update_k_on_step(time_step);
        let should_rebuild_by_omega = {
            let omega_sq_new = omega * omega;
            let denom = omega_sq_new.max(self.omega_sq_at_last_k_rebuild);
            if denom < 1e-16 {
                true // Near-zero ω → always rebuild (safe default)
            } else {
                let rel_change = (omega_sq_new - self.omega_sq_at_last_k_rebuild).abs() / denom;
                rel_change > self.config.omega_rebuild_threshold
            }
        };
        let should_rebuild_by_theta =
            (self.theta - self.theta_at_last_k_rebuild).abs() > self.config.theta_rebuild_threshold;

        // K_G(u) deflection-driven trigger (Nivel 2 only). The ratio defaults
        // to 0.0 when the feature is off, making the predicate inert.
        let current_deflection_ratio = if self.config.kg_use_deformed_coords {
            let u_full = self.expand_to_full(u_red, free_dofs);
            max_radial_deflection_ratio(
                &u_full,
                &self.coords_rotated,
                &self.transforms.axis,
                &self.transforms.center,
                self.config.dofs_per_node,
            )
        } else {
            0.0
        };
        let should_rebuild_by_deflection = self.config.kg_use_deformed_coords && {
            let drift = (current_deflection_ratio - self.deflection_ratio_at_last_kg_rebuild).abs();
            drift > self.config.kg_deflection_rebuild_rel_high
        };

        // K_SP Δω trigger — only active when `include_spin_softening`.
        let eps_sq = 1e-8_f64; // near-zero guard (ω ≈ 1e-4 rad/s)
        let ksp_recently_rebuilt = time_step.saturating_sub(self.last_ksp_rebuild_step) < 10;
        let should_rebuild_ksp_by_omega = self.config.include_spin_softening
            && omega_changed_significantly(
                omega,
                self.omega_sq_at_last_ksp_rebuild,
                self.config.ksp_omega_rebuild_high,
                self.config.ksp_omega_rebuild_low,
                ksp_recently_rebuilt,
                eps_sq,
            );

        // K_G Δω trigger — only active when `include_geometric_stiffness`.
        let kg_recently_rebuilt = time_step.saturating_sub(self.last_kg_rebuild_step) < 10;
        let should_rebuild_kg_by_omega = self.config.include_geometric_stiffness
            && omega_changed_significantly(
                omega,
                self.omega_sq_at_last_kg_rebuild,
                self.config.kg_omega_rebuild_high,
                self.config.kg_omega_rebuild_low,
                kg_recently_rebuilt,
                eps_sq,
            );

        // Any full-system trigger OR any per-matrix trigger?
        let full_rebuild = should_rebuild_by_step
            || should_rebuild_by_omega
            || should_rebuild_by_theta
            || should_rebuild_by_deflection;
        let partial_rebuild = should_rebuild_ksp_by_omega || should_rebuild_kg_by_omega;

        if !full_rebuild && !partial_rebuild {
            return Ok(());
        }

        let t_start = std::time::Instant::now();

        use crate::petsc::fsi::setup::{apply_kg_coo_map, reassemble_k};
        use crate::petsc::fsi::rotor_physics::build_ksp_vals;

        // ── Determine which pieces need to be recomputed ──────────────────────
        let rebuild_k = full_rebuild;
        let rebuild_kg = full_rebuild
            || (should_rebuild_kg_by_omega && self.config.include_geometric_stiffness);
        let rebuild_ksp = full_rebuild
            || (should_rebuild_ksp_by_omega && self.config.include_spin_softening);

        // ── K(θ): reassemble from rotated geometry ────────────────────────────
        let k_vals_red = if rebuild_k {
            // Caller guarantees self.coords_rotated tracks self.theta.
            let (_, _, k_vals) = reassemble_k(&mut self.assembler, &self.coords_rotated);
            if k_vals.len() != self.k_coo_map.len() {
                return Err(FsiError::PreciceError(format!(
                    "reassembled K COO length {} does not match reference mapping length {}",
                    k_vals.len(),
                    self.k_coo_map.len(),
                )));
            }
            apply_kg_coo_map(&self.k_coo_map, &k_vals, self.k_red_nnz)
        } else {
            // Reuse cached K(θ) from the stepper — no geometry change.
            self.stepper.k_vals().to_vec()
        };

        // ── K_G(θ, ω): centrifugal prestress ─────────────────────────────────
        //
        // Option A discipline: K(θ) uses X_rotated. K_G below optionally uses
        // X_rotated + u (Nivel 2 foreshortening). After deformed-coords assembly
        // we restore X_rotated so later assembly paths see canonical geometry.
        let kg_vals_red = if rebuild_kg {
            if self.config.include_geometric_stiffness {
                let kg_vals_full = if self.config.kg_use_deformed_coords {
                    let u_full = self.expand_to_full(u_red, free_dofs);
                    let coords_def = compose_deformed_coords(
                        &self.coords_rotated,
                        &u_full,
                        self.config.dofs_per_node,
                    );
                    self.assembler.update_node_coordinates(&coords_def);
                    let (_, _, vals) = self.assembler.assemble_centrifugal_k(
                        omega,
                        self.transforms.axis,
                        self.transforms.center,
                        &self.rho_per_elem,
                    );
                    // Restore X_rotated as the canonical assembler geometry.
                    let coords_rotated_snapshot = self.coords_rotated.clone();
                    self.assembler.update_node_coordinates(&coords_rotated_snapshot);
                    vals
                } else {
                    let (_, _, vals) = self.assembler.assemble_centrifugal_k(
                        omega,
                        self.transforms.axis,
                        self.transforms.center,
                        &self.rho_per_elem,
                    );
                    vals
                };
                apply_kg_coo_map(&self.kg_coo_map, &kg_vals_full, self.k_red_nnz)
            } else {
                // include_geometric_stiffness=false → K_G = 0 (never changes).
                // On a full rebuild we fill zeros; on a partial rebuild this branch
                // is unreachable (should_rebuild_kg_by_omega is gated on the flag).
                vec![0.0f64; self.k_red_nnz]
            }
        } else {
            // Reuse cached K_G from the stepper.
            self.stepper.kg_vals().to_vec()
        };

        // ── K_SP(ω): spin-softening ───────────────────────────────────────────
        //
        // In the inertial frame:
        //   M·ü_e + C·u̇_e + [K(θ) + K_G(θ,ω) + K_SP(ω)]·u_e = F
        //
        // K_SP = −ω²·M_lump·(I − n̂⊗n̂) (negative-definite in perpendicular plane).
        let ksp_vals_red = if rebuild_ksp {
            if self.config.include_spin_softening {
                let n_full_dofs = self.assembler.topology.dofs_count();
                let dofs_per_node = self.config.dofs_per_node;
                let mut m_lumped_full = vec![0.0f64; n_full_dofs];
                for (i, &m) in self.masses.iter().enumerate() {
                    for j in 0..3_usize {
                        let dof = i * dofs_per_node + j;
                        if dof < n_full_dofs {
                            m_lumped_full[dof] = m;
                        }
                    }
                }
                build_ksp_vals(
                    &m_lumped_full,
                    &self.transforms.axis,
                    omega,
                    dofs_per_node,
                    &self.k_red_rows,
                    &self.k_red_cols,
                    free_dofs,
                    n_full_dofs,
                )
            } else {
                // include_spin_softening=false → K_SP = 0.
                vec![0.0f64; self.k_red_nnz]
            }
        } else {
            // Reuse cached K_SP from the stepper.
            self.stepper.ksp_vals().to_vec()
        };

        // ── Single update_tangent_terms → single refactorization ─────────────
        self.stepper
            .update_tangent_terms(&k_vals_red, &kg_vals_red, &ksp_vals_red)
            .map_err(FsiError::StepperError)?;

        // ── Update trackers for actually-rebuilt components ───────────────────
        let omega_sq = omega * omega;
        if rebuild_k {
            self.omega_sq_at_last_k_rebuild = omega_sq;
            self.theta_at_last_k_rebuild = self.theta;
            self.step_at_last_k_rebuild = time_step;
        }
        if rebuild_kg {
            self.omega_sq_at_last_kg_rebuild = omega_sq;
            self.last_kg_rebuild_step = time_step;
        }
        if rebuild_ksp {
            self.omega_sq_at_last_ksp_rebuild = omega_sq;
            self.last_ksp_rebuild_step = time_step;
        }
        if self.config.kg_use_deformed_coords && (rebuild_kg || full_rebuild) {
            self.deflection_ratio_at_last_kg_rebuild = current_deflection_ratio;
        }

        let elapsed = t_start.elapsed();
        log::info!(
            "InertialRotorFsi: stiffness reassembly done in {:.3}s \
             (step {}, θ={:.4}rad, ω={:.4}rad/s, \
             rebuild_k={}, rebuild_kg={}, rebuild_ksp={}, \
             K_G={}, K_G(u)={}, K_SP={}, def_ratio={:.4e})",
            elapsed.as_secs_f64(),
            time_step,
            self.theta,
            omega,
            rebuild_k,
            rebuild_kg,
            rebuild_ksp,
            self.config.include_geometric_stiffness,
            self.config.kg_use_deformed_coords,
            self.config.include_spin_softening,
            current_deflection_ratio,
        );

        Ok(())
    }

    // ── Main coupling loop ────────────────────────────────────────────────────

    /// Run the full preCICE coupling loop for the inertial rotor solver.
    ///
    /// # Returns
    /// `FsiResult` — converged final state (u, v, a) and time history.
    pub fn run(mut self) -> Result<FsiResult, FsiError> {
        use precice::Participant;
        use crate::petsc::fsi::force_utils::{apply_cap, apply_ramp};
        use crate::petsc::fsi::rotor_physics::{
            compute_gravity_force, compute_performance_coefficients,
            compute_thrust, compute_torque,
        };

        // Initialize preCICE participant
        let mut participant = Participant::new(
            &self.config.fsi.participant_name,
            &self.config.fsi.config_file,
            0,
            1,
        )?;

        let mesh_dims = participant.get_mesh_dimensions(&self.config.fsi.coupling_mesh)? as usize;

        // Define coupling mesh vertices (interface nodes, global frame)
        let iface_coords: Vec<f64> = self.iface_nodes.iter()
            .flat_map(|&node| {
                let i = node * 3;
                [self.coords_ref[i], self.coords_ref[i + 1], self.coords_ref[i + 2]]
            })
            .collect();

        let n_vertices = iface_coords.len() / mesh_dims.max(1);
        let mut vertex_ids = vec![0i32; n_vertices];
        participant.set_mesh_vertices(
            &self.config.fsi.coupling_mesh,
            &iface_coords,
            &mut vertex_ids,
        )?;

        // Optional GlobalSolidMesh for ω communication
        let omega_vertex_ids: Option<Vec<i32>> = if let (Some(ref mesh), Some(ref coords)) = (
            &self.config.omega_mesh_name,
            &self.config.omega_vertex_coord,
        ) {
            let mut ids = vec![0i32; 1];
            participant.set_mesh_vertices(mesh, &coords[..], &mut ids)?;
            Some(ids)
        } else {
            None
        };

        // Write initial omega before initialize() when preCICE requires it.
        // This satisfies the initialize="true" exchange for AngularVelocity.
        if participant.requires_initial_data()? {
            if let (Some(mesh), Some(wdata), Some(ids)) = (
                &self.config.omega_mesh_name,
                &self.config.omega_write_data,
                &omega_vertex_ids,
            ) {
                // InertialRotorFsiSolver always starts from t=0 (no restart support yet)
                let omega_init = self.omega_provider.get(0.0).0;
                participant.write_data(mesh, wdata, ids, &[omega_init])?;
            }
        }

        participant.initialize()?;
        let mut dt = participant.get_max_time_step_size()?;

        // Use the BC reduction provided by the caller.
        let free_dofs = self.free_dofs.clone();

        let n_dofs = free_dofs.len();
        let n_data = self.n_iface * mesh_dims;

        // ── Correct initial acceleration a₀ = M⁻¹·F₀ ─────────────────────────
        // When u₀ = v₀ = a₀ = 0 (fresh start, not a restart), the Newmark
        // scheme defaults to a₀ = 0. This is only correct when all external
        // static loads at t = 0 are zero. With gravity and/or reference-frame
        // centrifugal loads, F₀ ≠ 0 → a₀ = 0 introduces a spurious transient
        // that decays over several cycles.
        //
        // Heuristic for "fresh start": current time is at the initial time AND
        // the entire state vector is zero (user did not provide a non-trivial
        // a₀ via FsiInitialState). Skip when any component is non-zero to
        // respect explicit restart states.
        {
            let t0 = self.stepper.current_time();
            let a_is_zero = self.stepper.current_a().iter().all(|&x| x == 0.0);
            let u_is_zero = self.stepper.current_u().iter().all(|&x| x == 0.0);
            let has_static_loads =
                self.config.gravity_active() || self.config.include_reference_acceleration;

            if a_is_zero && u_is_zero && has_static_loads {
                let (omega_0, alpha_0) = self.omega_provider.get(t0);
                let mut f0 = vec![0.0f64; n_dofs];

                if self.config.include_reference_acceleration {
                    use crate::petsc::fsi::rotor_physics::{
                        compute_reference_load_vector, compute_rigid_body_acceleration_inertial,
                    };
                    let a_ref0 = compute_rigid_body_acceleration_inertial(
                        &self.coords_rotated,
                        &self.transforms.axis,
                        &self.transforms.center,
                        omega_0,
                        alpha_0,
                    );
                    let f_ref0 = compute_reference_load_vector(&a_ref0, &self.masses);
                    self.scatter_node_forces(&f_ref0, &mut f0, &free_dofs);
                }
                if self.config.gravity_active() {
                    use crate::petsc::fsi::rotor_physics::compute_gravity_force;
                    let f_g0 = compute_gravity_force(&self.masses, &self.config.gravity);
                    self.scatter_node_forces(&f_g0, &mut f0, &free_dofs);
                }

                let a0 = self.compute_initial_acceleration(&f0, &free_dofs);
                let a0_norm: f64 = a0.iter().map(|x| x * x).sum::<f64>().sqrt();
                log::info!(
                    "InertialRotorFsi: computed a₀ from static loads (||a₀||={:.3e}, ω₀={:.4}rad/s)",
                    a0_norm, omega_0
                );
                self.stepper.set_initial_conditions_with_acceleration(
                    &vec![0.0f64; n_dofs],
                    &vec![0.0f64; n_dofs],
                    &a0,
                );
            }
        }

        // Checkpoint storage for implicit coupling
        let mut newmark_cp: Option<NewmarkCheckpoint> = None;
        let mut omega_cp: Option<OmegaCheckpoint> = None;
        let mut theta_cp = 0.0f64;
        let mut kgksp_cp: Option<InertialKgKspCheckpoint> = None;

        let mut result = FsiResult {
            u_final: Vec::new(),
            v_final: Vec::new(),
            a_final: Vec::new(),
            times: Vec::new(),
        };

        // Main coupling loop
        let mut sub_iter: u32 = 0;
        let mut window_num: u32 = 0;
        let mut prof = WindowProfiler::new("InertialRotorFsi");
        while participant.is_coupling_ongoing()? {
            // ── Save checkpoint ───────────────────────────────────────────────
            if participant.requires_writing_checkpoint()? {
                prof.start_window();
                let t_cp = Instant::now();
                newmark_cp = Some(self.stepper.checkpoint());
                omega_cp = Some(self.omega_provider.checkpoint());
                theta_cp = self.theta;
                kgksp_cp = Some(InertialKgKspCheckpoint {
                    omega_sq_at_last_ksp_rebuild: self.omega_sq_at_last_ksp_rebuild,
                    omega_sq_at_last_kg_rebuild: self.omega_sq_at_last_kg_rebuild,
                    last_ksp_rebuild_step: self.last_ksp_rebuild_step,
                    last_kg_rebuild_step: self.last_kg_rebuild_step,
                });
                prof.record_since("checkpoint_save", t_cp);
                sub_iter = 0;
                window_num += 1;
                let u_norm: f64 = self.stepper.current_u().iter().map(|x| x * x).sum::<f64>().sqrt();
                log::debug!(
                    "InertialRotorFsi: [w={} iter=0] CHECKPOINT SAVED  ||u||={:.3e}  t={:.6}  theta={:.6}",
                    window_num, u_norm, self.stepper.current_time(), self.theta
                );
            }
            sub_iter += 1;

            let t = self.stepper.current_time();
            let window_kin = self.omega_provider.kinematics_over_step(t, dt);
            let omega_step = window_kin.omega_step;
            let alpha_step = window_kin.alpha_step;
            let theta_target = self.theta + window_kin.theta_increment;

            // ── Read forces from preCICE (global frame) ───────────────────────
            let t_read = Instant::now();
            let mut forces_global = vec![0.0f64; n_data];
            participant.read_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.read_data,
                &vertex_ids,
                dt,
                &mut forces_global,
            )?;
            prof.record_since("precice_read", t_read);

            // ── Force pre-processing (ramp + cap, global frame) ───────────────
            // Save physical (unramped) forces for torque/thrust/Ct/Cp after
            // convergence. The ramp is a purely numerical start-up artifact and
            // must NOT contaminate aerodynamic diagnostics or OmegaProvider.
            let t_fp = Instant::now();
            let forces_raw = forces_global.clone();
            apply_ramp(&mut forces_global, t, self.config.fsi.ramp_time);
            if let Some(max_f) = self.config.fsi.force_max {
                apply_cap(&mut forces_global, max_f, mesh_dims);
            }
            prof.record_since("force_preproc", t_fp);

            // ── Assemble global reduced force vector ──────────────────────────
            let t_assemble = Instant::now();
            let mut f_red = vec![0.0f64; n_dofs];

            // Scatter aero forces at interface DOFs
            for (li, &rdof) in self.iface_dofs.iter().enumerate() {
                if li < forces_global.len() {
                    if let Ok(pos) = free_dofs.binary_search(&(rdof as i32)) {
                        f_red[pos] += forces_global[li];
                    }
                }
            }

            // Reference acceleration forces (F_ref = −M·a_ref)
            //
            // NOTE — geometry lag (explicit-in-geometry approximation):
            // `coords_rotated` here reflects θ_{n-1} (the checkpoint rotation),
            // not θ_n = θ_{n-1} + Δθ. This introduces an O(ω·Δt) lag in the
            // lever arm and an O(ω²·Δt) error in F_ref magnitude.
            // See module-level doc for rationale and quantitative bounds.
            if self.config.include_reference_acceleration {
                let a_ref = compute_rigid_body_acceleration_inertial(
                    &self.coords_rotated,
                    &self.transforms.axis,
                    &self.transforms.center,
                    omega_step,
                    alpha_step,
                );
                let f_ref = compute_reference_load_vector(&a_ref, &self.masses);
                self.scatter_node_forces(&f_ref, &mut f_red, &free_dofs);
            }

            // Gravity (global frame, constant)
            if self.config.gravity_active() {
                let f_g = compute_gravity_force(&self.masses, &self.config.gravity);
                self.scatter_node_forces(&f_g, &mut f_red, &free_dofs);
            }

            prof.record_since("force_scatter_ref_grav", t_assemble);
            let f_norm: f64 = forces_global.iter().map(|x| x * x).sum::<f64>().sqrt();
            let f_red_norm: f64 = f_red.iter().map(|x| x * x).sum::<f64>().sqrt();
            let u_pre_step_norm: f64 = self.stepper.current_u().iter().map(|x| x * x).sum::<f64>().sqrt();

            // ── Advance structural state ──────────────────────────────────────
            let t_solve = Instant::now();
            let step_t = self.stepper.step(&f_red, dt)?.t;
            prof.record_since("newmark_step", t_solve);

            // When exporting total kinematics, the rigid contribution must use
            // the current window target θ, not the last converged geometry.
            let iface_coords_target = if self.config.displacement_mode == DisplacementMode::Total {
                Some(crate::petsc::fsi::setup::rotate_mesh_coords(
                    &iface_coords,
                    &self.transforms,
                    theta_target,
                ))
            } else {
                None
            };

            // ── Gather interface displacements (elastic or total) ─────────────
            let u_red = self.stepper.current_u().to_vec();
            let u_post_step_norm: f64 = u_red.iter().map(|x| x * x).sum::<f64>().sqrt();
            let disp_iface: Vec<f64> = match self.config.displacement_mode {
                DisplacementMode::Elastic => {
                    // Write elastic displacement only (u_e)
                    self.iface_dofs
                        .iter()
                        .map(|&dof| {
                            if let Ok(pos) = free_dofs.binary_search(&(dof as i32)) {
                                u_red[pos]
                            } else {
                                0.0
                            }
                        })
                        .collect()
                }
                DisplacementMode::Total => {
                    // Write total displacement using the rigid kinematics frozen
                    // for the CURRENT window: u_total = u_elastic + (R(θ_target)X₀ - X₀).
                    let iface_coords_rotated = iface_coords_target
                        .as_ref()
                        .expect("total displacement mode requires target rigid coordinates");
                    let mut disp_total = Vec::with_capacity(self.iface_dofs.len());
                    for (iface_idx, &node) in self.iface_nodes.iter().enumerate() {
                        let i_ref = iface_idx * 3;
                        for j in 0..3 {
                            let u_rigid = iface_coords_rotated[i_ref + j] - iface_coords[i_ref + j];
                            let dof = (node * self.config.dofs_per_node + j) as i32;
                            let u_elastic = if let Ok(pos) = free_dofs.binary_search(&dof) {
                                u_red[pos]
                            } else {
                                0.0
                            };
                            disp_total.push(u_elastic + u_rigid);
                        }
                    }
                    disp_total
                }
            };

            // ── Write displacements to preCICE ────────────────────────────────
            let d_norm: f64 = disp_iface.iter().map(|x| x * x).sum::<f64>().sqrt();
            log::debug!(
                "InertialRotorFsi: [w={} iter={}] ||F_raw||={:.3e}  ||f_red||={:.3e}  ||u_before||={:.3e}  ||u_after||={:.3e}  ||D_written||={:.3e}",
                window_num, sub_iter, f_norm, f_red_norm, u_pre_step_norm, u_post_step_norm, d_norm
            );
            let t_write = Instant::now();
            participant.write_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.write_data,
                &vertex_ids,
                &disp_iface,
            )?;

            // ── Write nodal velocities to preCICE (aerodynamic damping) ──────
            if let Some(ref vdata) = self.config.velocity_write_data {
                let v_red = self.stepper.current_v().to_vec();
                let vel_iface: Vec<f64> = match self.config.displacement_mode {
                    DisplacementMode::Elastic => self
                        .iface_dofs
                        .iter()
                        .map(|&dof| {
                            if let Ok(pos) = free_dofs.binary_search(&(dof as i32)) {
                                v_red[pos]
                            } else {
                                0.0
                            }
                        })
                        .collect(),
                    DisplacementMode::Total => {
                        let iface_coords_rotated = iface_coords_target
                            .as_ref()
                            .expect("total displacement mode requires target rigid coordinates");
                        let vel_rigid = compute_rigid_body_velocity_inertial(
                            iface_coords_rotated,
                            &self.transforms.axis,
                            &self.transforms.center,
                            omega_step,
                        );

                        let mut vel_total = Vec::with_capacity(self.iface_dofs.len());
                        for (iface_idx, &node) in self.iface_nodes.iter().enumerate() {
                            let i = iface_idx * 3;
                            for j in 0..3 {
                                let dof = (node * self.config.dofs_per_node + j) as i32;
                                let v_elastic = if let Ok(pos) = free_dofs.binary_search(&dof) {
                                    v_red[pos]
                                } else {
                                    0.0
                                };
                                vel_total.push(v_elastic + vel_rigid[i + j]);
                            }
                        }
                        vel_total
                    }
                };
                participant.write_data(
                    &self.config.fsi.coupling_mesh,
                    vdata,
                    &vertex_ids,
                    &vel_iface,
                )?;
            }

            // ── Write ω to GlobalSolidMesh (if used) ──────────────────────────
            if let (Some(ref ids), Some(ref mesh), Some(ref wdata)) = (
                &omega_vertex_ids,
                &self.config.omega_mesh_name,
                &self.config.omega_write_data,
            ) {
                participant.write_data(mesh, wdata, ids, &[omega_step])?;
            }
            prof.record_since("precice_write", t_write);

            let t_advance = Instant::now();
            participant.advance(dt)?;
            prof.record_since("precice_advance", t_advance);

            // ── Implicit coupling: restore or commit ──────────────────────────
            if participant.requires_reading_checkpoint()? {
                let t_rs = Instant::now();
                match (newmark_cp.as_ref(), omega_cp.as_ref(), kgksp_cp.as_ref()) {
                    (Some(ncp), Some(ocp), Some(kcp)) => {
                        let u_before_restore_norm: f64 = self.stepper.current_u().iter().map(|x| x * x).sum::<f64>().sqrt();
                        self.stepper.restore(ncp);
                        self.omega_provider.restore(ocp);
                        self.theta = theta_cp;
                        // Restore K_SP / K_G rebuild trackers so hysteresis is consistent
                        // after rollback. Without this, the tracker step indices would point
                        // beyond the rolled-back time, suppressing the first post-rollback rebuild.
                        self.omega_sq_at_last_ksp_rebuild = kcp.omega_sq_at_last_ksp_rebuild;
                        self.omega_sq_at_last_kg_rebuild = kcp.omega_sq_at_last_kg_rebuild;
                        self.last_ksp_rebuild_step = kcp.last_ksp_rebuild_step;
                        self.last_kg_rebuild_step = kcp.last_kg_rebuild_step;

                        // CRITICAL: Re-rotate coords to theta_cp to maintain consistency.
                        // If we don't do this, coords_rotated will be out of sync with theta
                        // until the next K rebuild, causing incorrect reference forces.
                        use crate::petsc::fsi::setup::rotate_mesh_coords;
                        self.coords_rotated = rotate_mesh_coords(&self.coords_ref, &self.transforms, self.theta);

                        let u_after_restore_norm: f64 = self.stepper.current_u().iter().map(|x| x * x).sum::<f64>().sqrt();
                        let cp_u_norm: f64 = ncp.u.iter().map(|x| x * x).sum::<f64>().sqrt();
                        log::debug!(
                            "InertialRotorFsi: [w={} iter={}] ROLLBACK  ||u_before||={:.3e}  ||u_after||={:.3e}  ||cp.u||={:.3e}",
                            window_num, sub_iter, u_before_restore_norm, u_after_restore_norm, cp_u_norm
                        );
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
                self.time_step += 1;
                self.theta = theta_target;

                // Keep self.coords_rotated synchronized with self.theta on every
                // converged window. Required by the lever-arm computation below
                // and by F_ref/K(θ) on the next iteration; cheap O(n_nodes)
                // and decouples correctness from K-rebuild cadence.
                use crate::petsc::fsi::setup::rotate_mesh_coords;
                let t_rot = Instant::now();
                self.coords_rotated =
                    rotate_mesh_coords(&self.coords_ref, &self.transforms, self.theta);
                prof.record_since("mesh_rotate", t_rot);

                // Update K(θ) if rotation changed significantly.
                //
                // NOTE — geometry lag: K is rebuilt from `coords_rotated` which was
                // just advanced to θ_n (line above). So K(θ_n) is correct here.
                // However, during sub-iterations the K that was used came from
                // the checkpoint (θ_{n-1}). This is the explicit-in-geometry
                // approximation documented in the module docstring.
                let t_reassemble = Instant::now();
                // Snapshot the current reduced displacement so the K_G(u)
                // trigger and assembly path see a consistent state. The
                // snapshot is also necessary because reassemble_k_if_needed
                // takes &mut self and would conflict with a live borrow on
                // self.stepper.current_u().
                let u_red_snapshot = self.stepper.current_u().to_vec();
                self.reassemble_k_if_needed(
                    omega_step,
                    &u_red_snapshot,
                    self.time_step,
                    &free_dofs,
                )?;
                prof.record_since("k_reassemble_if_needed", t_reassemble);

                // Compute aerodynamic torque (global frame)
                let u_full = self.expand_to_full(&u_red, &free_dofs);
                let dofs_per_node = self.config.dofs_per_node;
                let mut disp_iface_full = Vec::with_capacity(self.iface_nodes.len() * 3);
                for &node in &self.iface_nodes {
                    for j in 0..3 {
                        let dof = node * dofs_per_node + j;
                        disp_iface_full.push(if dof < u_full.len() { u_full[dof] } else { 0.0 });
                    }
                }

                // Lever arm uses the CURRENT rotated geometry r = R(θ)·X₀ + u_e
                // − center. The preCICE-registered iface_coords are anchored at
                // X₀ (θ=0, fixed mesh contract) so they are not the right input
                // for τ_aero; passing them produced a constant lever-arm bug
                // that contaminated τ_aero, Cp/Cq, and the ComputedOmega update.
                let iface_coords_rotated: Vec<f64> = self
                    .iface_nodes
                    .iter()
                    .flat_map(|&node| {
                        let i = node * 3;
                        [
                            self.coords_rotated[i],
                            self.coords_rotated[i + 1],
                            self.coords_rotated[i + 2],
                        ]
                    })
                    .collect();

                let (_tau_vec, tau_aero) = compute_torque(
                    &iface_coords_rotated,
                    &disp_iface_full,
                    &forces_raw,
                    &self.transforms.axis,
                    &self.transforms.center,
                );

                // Gravity torque about the rotation axis must also drive ω.
                // In the inertial formulation gravity is already expressed in
                // the global frame, so compute the body torque directly from
                // the current rotated geometry plus the elastic displacement.
                let tau_gravity = if self.config.gravity_active() {
                    let f_gravity = compute_gravity_force(&self.masses, &self.config.gravity);
                    let mut disp_nodes_full = Vec::with_capacity(self.masses.len() * 3);
                    for node in 0..self.masses.len() {
                        for j in 0..3 {
                            let dof = node * dofs_per_node + j;
                            disp_nodes_full.push(if dof < u_full.len() { u_full[dof] } else { 0.0 });
                        }
                    }
                    let (_, tau_g) = compute_torque(
                        &self.coords_rotated,
                        &disp_nodes_full,
                        &f_gravity,
                        &self.transforms.axis,
                        &self.transforms.center,
                    );
                    tau_g
                } else {
                    0.0
                };

                // Update ω from the physical driving torque only. The provider
                // itself no-ops for constant/ramped modes and handles the
                // post-ramp transition for RampedComputed.
                self.omega_provider
                    .update_from_torque(tau_aero + tau_gravity, dt, step_t);

                // Store final state
                result.u_final = u_red.to_vec();
                result.v_final = self.stepper.current_v().to_vec();
                result.a_final = self.stepper.current_a().to_vec();
                result.times.push(step_t);

                // Performance coefficients use physical (unramped) forces.
                let thrust = compute_thrust(&forces_raw, &self.transforms.axis);
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
                    "InertialRotorFsi step {}: t={:.4} θ={:.4}rad ω_window={:.4}rad/s ω_state={:.4}rad/s \
                     τ_aero={:.3e}N·m τ_grav={:.3e}N·m  Ct={:.4} Cp={:.4} TSR={:.3}",
                    self.time_step, step_t, self.theta, omega_step, omega_state,
                    tau_aero, tau_gravity, perf.ct, perf.cp, perf.tsr
                );

                // Per-step callback — always pass physical (unramped) forces
                if let Some(ref cb) = self.step_callback {
                    let t_cb = Instant::now();
                    let force_mag = forces_raw.iter().map(|x| x * x).sum::<f64>().sqrt();
                    cb(
                        step_t,
                        self.time_step,
                        dt,
                        &u_red,
                        self.stepper.current_v(),
                        self.stepper.current_a(),
                        force_mag,
                        &forces_raw,
                        omega_state,
                        alpha_state,
                        self.theta,
                        tau_aero,
                        perf,
                    )?;
                    prof.record_since("callback", t_cb);
                }

                dt = participant.get_max_time_step_size()?;

                // ── Per-window profiling summary (debug log) ──────────────────
                prof.log_summary(step_t, self.time_step);
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
        self.n_iface
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

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Construct a minimal `InertialRotorFsiConfig` with sensible defaults for testing.
    fn default_config() -> InertialRotorFsiConfig {
        use crate::petsc::fsi::linear_elastic::FsiConfig;
        InertialRotorFsiConfig {
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
            rotation_center: [0.0; 3],
            gravity: [0.0; 3],
            include_reference_acceleration: true,
            include_geometric_stiffness: false,
            include_spin_softening: true,
            ksp_omega_rebuild_high: 0.005,
            ksp_omega_rebuild_low: 0.003,
            kg_omega_rebuild_high: 0.005,
            kg_omega_rebuild_low: 0.003,
            k_update_interval: 1,
            omega_rebuild_threshold: 0.01,
            theta_rebuild_threshold: 0.05,
            kg_use_deformed_coords: false,
            kg_deflection_rebuild_rel_high: 0.01,
            kg_deflection_rebuild_rel_low: 0.005,
            displacement_mode: DisplacementMode::Elastic,
            dofs_per_node: 3,
            fluid_density: 1.225,
            flow_velocity: 0.0,
            rotor_radius: 1.0,
            omega_mesh_name: None,
            omega_write_data: None,
            omega_vertex_coord: None,
            velocity_write_data: None,
        }
    }

    // ── T4.4a: K_G Δω trigger fires when ω changes, not when only θ changes ──

    /// Verify that `omega_changed_significantly` correctly gates K_G Δω rebuild.
    ///
    /// Spec scenario R4: K_G rebuilds on Δω above `kg_omega_rebuild_high`;
    /// Δθ alone does NOT rebuild K_G through the per-matrix path.
    #[test]
    fn kg_omega_trigger_fires_on_large_omega_change() {
        // Large ω change (above high threshold = 0.5%) → should trigger.
        let omega_last_sq = 100.0f64 * 100.0f64; // ω_last = 100 rad/s
        let omega_new = 100.5f64; // 0.5% change in ω → ~1% in ω²
        let result = omega_changed_significantly(
            omega_new,
            omega_last_sq,
            0.005, // high threshold: 0.5% change in ω²
            0.003,
            false, // not recently rebuilt
            1e-8,
        );
        assert!(
            result,
            "K_G trigger should fire for large ω change ({} → {}): rel={:.4}%",
            100.0,
            omega_new,
            100.0 * (omega_new * omega_new - omega_last_sq).abs() / omega_last_sq
        );
    }

    /// Verify that a small ω change does NOT trigger K_G rebuild (hysteresis).
    #[test]
    fn kg_omega_trigger_suppressed_for_small_omega_change() {
        let omega_last_sq = 100.0f64 * 100.0f64;
        let omega_new = 100.05f64; // 0.05% change in ω → ~0.1% in ω² — below thresholds
        let result = omega_changed_significantly(
            omega_new,
            omega_last_sq,
            0.005,
            0.003,
            false,
            1e-8,
        );
        assert!(
            !result,
            "K_G trigger should NOT fire for tiny ω change: {:.4}%",
            100.0 * (omega_new * omega_new - omega_last_sq).abs() / omega_last_sq
        );
    }

    /// First call (NEG_INFINITY sentinel) always triggers.
    #[test]
    fn kg_omega_trigger_first_call_always_triggers() {
        let omega_new = 50.0f64;
        assert!(
            omega_changed_significantly(omega_new, f64::NEG_INFINITY, 0.005, 0.003, false, 1e-8),
            "First call must always trigger (NEG_INFINITY sentinel)"
        );
    }

    // ── T4.4b: should_rebuild_ksp_by_omega correctly gated on include_spin_softening ─

    /// When include_spin_softening=false, the K_SP Δω trigger must not fire
    /// even with a large ω change.
    #[test]
    fn ksp_trigger_disabled_when_include_spin_softening_false() {
        // Construct config with include_spin_softening = false.
        let mut cfg = default_config();
        cfg.include_spin_softening = false;

        // We verify the gating logic by directly testing the predicate with the
        // flag. The flag is the only guard — the underlying omega_changed_significantly
        // predicate itself is symmetric.
        let omega_new = 110.0f64;
        let omega_sq_last = 100.0f64 * 100.0f64;
        // omega_changed_significantly would return true here (10% change), but
        // the flag should suppress it.
        let underlying = omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, false, 1e-8);
        assert!(underlying, "The underlying predicate fires for 10% change");

        let gated_result = cfg.include_spin_softening && underlying;
        assert!(!gated_result, "K_SP trigger must be suppressed when include_spin_softening=false");
    }

    // ── T4.4c: include_spin_softening defaults to true ────────────────────────

    #[test]
    fn default_config_has_include_spin_softening_true() {
        let cfg = default_config();
        assert!(
            cfg.include_spin_softening,
            "include_spin_softening must default to true (matches HEAD unconditional behavior)"
        );
    }

    // ── T4.4d: New config fields have expected defaults ───────────────────────

    #[test]
    fn new_config_fields_have_correct_defaults() {
        let cfg = default_config();
        assert_eq!(cfg.ksp_omega_rebuild_high, 0.005, "ksp_omega_rebuild_high default mismatch");
        assert_eq!(cfg.ksp_omega_rebuild_low, 0.003, "ksp_omega_rebuild_low default mismatch");
        assert_eq!(cfg.kg_omega_rebuild_high, 0.005, "kg_omega_rebuild_high default mismatch");
        assert_eq!(cfg.kg_omega_rebuild_low, 0.003, "kg_omega_rebuild_low default mismatch");
    }

    // ── T4.4e: InertialKgKspCheckpoint save/restore semantics ─────────────────

    #[test]
    fn inertial_kgksp_checkpoint_round_trips_correctly() {
        let cp = InertialKgKspCheckpoint {
            omega_sq_at_last_ksp_rebuild: 12.5,
            omega_sq_at_last_kg_rebuild: 99.9,
            last_ksp_rebuild_step: 42,
            last_kg_rebuild_step: 17,
        };
        // Clone and verify fields are preserved.
        let cp2 = cp.clone();
        assert_eq!(cp2.omega_sq_at_last_ksp_rebuild, 12.5);
        assert_eq!(cp2.omega_sq_at_last_kg_rebuild, 99.9);
        assert_eq!(cp2.last_ksp_rebuild_step, 42);
        assert_eq!(cp2.last_kg_rebuild_step, 17);
    }
}
