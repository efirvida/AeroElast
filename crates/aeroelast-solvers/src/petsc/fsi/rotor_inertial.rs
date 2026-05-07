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
/// # Feature gate
/// Compiled only with `--features fsi`.

use aeroelast_core::assembly::assembler::MeshAssembler;

use crate::petsc::elasticity::dynamic_newmark::{NewmarkCheckpoint, NewmarkStepper};
use crate::petsc::fsi::force_utils::{apply_cap, apply_ramp};
use crate::petsc::fsi::linear_elastic::{FsiConfig, FsiError, FsiInitialState, FsiResult};
use crate::petsc::fsi::rotor_physics::{
    OmegaCheckpoint, OmegaProvider, PerformanceCoefficients, RotorTransforms,
    compute_gravity_force, compute_performance_coefficients,
    compute_rigid_body_acceleration_inertial, compute_reference_load_vector,
    compute_thrust, compute_torque,
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

    // ── Stiffness reassembly ──────────────────────────────────────────────────
    /// Rebuild K(θ) every N converged steps. `0` and `1` both mean every step.
    /// For inertial solver, K(θ) **must** be updated when geometry rotates,
    /// so this controls the frequency (typically `1` for every window).
    pub k_update_interval: usize,
    /// Relative |Δ(ω²)|/ω² threshold to trigger a K(θ) rebuild when θ changes
    /// significantly between windows. Default: `0.01` (1% change in ω²).
    pub omega_rebuild_threshold: f64,

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

// ── Solver state ───────────────────────────────────────────────────────────────

/// Internal state of the inertial rotor FSI solver.
struct InertialRotorFsiSolver {
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

    // ── Interface state ───────────────────────────────────────────────────────
    /// Number of interface (coupling) nodes.
    n_iface: usize,
    /// DOF indices for the interface nodes (length = n_iface × dofs_per_node).
    iface_dofs: Vec<usize>,
    /// Interface node indices (0-based, length = n_iface).
    iface_nodes: Vec<usize>,

    // ── Convergence tracking ──────────────────────────────────────────────────
    /// 1-based converged time step counter.
    time_step: usize,

    // ── Stiffness update tracking ─────────────────────────────────────────────
    /// ω² at the last K(θ) rebuild (for relative change check).
    omega_sq_at_last_k_rebuild: f64,
    /// Step index of the last K(θ) rebuild.
    step_at_last_k_rebuild: usize,

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
    ) -> Result<Self, FsiError> {
        let n_iface = iface_nodes.len();
        let dofs_per_node = config.dofs_per_node;

        // Build interface DOF indices
        let mut iface_dofs = Vec::with_capacity(n_iface * dofs_per_node);
        for &node in &iface_nodes {
            for d in 0..dofs_per_node {
                iface_dofs.push(node * dofs_per_node + d);
            }
        }

        // Extract reference coordinates from the assembler
        let coords_ref = assembler.topology.node_coords.clone();
        let coords_rotated = coords_ref.clone();

        // Normalize rotation axis
        let axis = config.rotation_axis;
        let axis_norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if axis_norm < 1e-14 {
            return Err(FsiError::InvalidConfig(
                "rotation_axis must have non-zero magnitude".to_string(),
            ));
        }
        let axis_normalized = [axis[0] / axis_norm, axis[1] / axis_norm, axis[2] / axis_norm];

        let transforms = RotorTransforms::new(axis_normalized, config.rotation_center);

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
            n_iface,
            iface_dofs,
            iface_nodes,
            time_step: 0,
            omega_sq_at_last_k_rebuild: f64::NEG_INFINITY,
            step_at_last_k_rebuild: 0,
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
            &state.displacement,
            &state.velocity,
            &state.acceleration,
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

    /// Rebuild K(θ) when geometry has rotated significantly.
    ///
    /// Calls `rotate_mesh_coords()` → `assembler.update_node_coordinates()` →
    /// `assembler.assemble_k()` → `reduce_coo()` → `stepper.update_elastic_stiffness()`.
    fn reassemble_k_if_needed(
        &mut self,
        omega: f64,
        time_step: usize,
        free_dofs: &[i32],
    ) -> Result<(), FsiError> {
        // Check if rebuild is due (by step interval or ω² change)
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

        if !should_rebuild_by_step && !should_rebuild_by_omega {
            return Ok(());
        }

        // Rotate mesh coordinates to current θ
        use crate::petsc::fsi::setup::rotate_mesh_coords;
        self.coords_rotated = rotate_mesh_coords(&self.coords_ref, &self.transforms, self.theta);

        // Reassemble K(θ) at rotated geometry
        use crate::petsc::fsi::setup::reassemble_k;
        let (k_rows, k_cols, k_vals) = reassemble_k(&mut self.assembler, &self.coords_rotated);

        // Reduce to free-DOF system
        use crate::petsc::fsi::setup::reduce_coo;
        let (_, _, k_vals_red) = reduce_coo(&k_rows, &k_cols, &k_vals, free_dofs);

        // Update stepper K
        self.stepper.update_elastic_stiffness(&k_vals_red)
            .map_err(FsiError::StepperError)?;

        self.omega_sq_at_last_k_rebuild = omega * omega;
        self.step_at_last_k_rebuild = time_step;

        log::info!(
            "InertialRotorFsi: K(θ) rebuilt at step {}, θ={:.4}rad, ω={:.4}rad/s",
            time_step, self.theta, omega
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

        let mesh_dims = participant.get_mesh_dimensions(&self.config.fsi.coupling_mesh)?;
        let dt = participant.initialize()?;

        // Define coupling mesh vertices (interface nodes, global frame)
        let iface_coords: Vec<f64> = self.iface_nodes.iter()
            .flat_map(|&node| {
                let i = node * 3;
                [self.coords_ref[i], self.coords_ref[i + 1], self.coords_ref[i + 2]]
            })
            .collect();

        let vertex_ids = participant.set_mesh_vertices(
            &self.config.fsi.coupling_mesh,
            &iface_coords,
        )?;

        // Optional GlobalSolidMesh for ω communication
        let omega_vertex_ids = if let (Some(ref mesh), Some(ref coords)) = (
            &self.config.omega_mesh_name,
            &self.config.omega_vertex_coord,
        ) {
            Some(participant.set_mesh_vertices(mesh, &coords[..])?)
        } else {
            None
        };

        // Determine free DOFs (assuming all DOFs are free for now — caller should provide)
        // For a real implementation, this would come from boundary conditions.
        let n_full_dofs = self.assembler.topology.dofs_count();
        let free_dofs: Vec<i32> = (0..n_full_dofs as i32).collect();

        let n_dofs = free_dofs.len();
        let n_data = self.n_iface * mesh_dims;

        // Checkpoint storage for implicit coupling
        let mut newmark_cp: Option<NewmarkCheckpoint> = None;
        let mut omega_cp: Option<OmegaCheckpoint> = None;
        let mut theta_cp = 0.0f64;

        let mut result = FsiResult {
            u_final: Vec::new(),
            v_final: Vec::new(),
            a_final: Vec::new(),
            times: Vec::new(),
        };

        // Main coupling loop
        let mut dt = dt;
        while participant.is_coupling_ongoing()? {
            // ── Save checkpoint ───────────────────────────────────────────────
            if participant.requires_writing_checkpoint()? {
                newmark_cp = Some(self.stepper.checkpoint());
                omega_cp = Some(self.omega_provider.checkpoint());
                theta_cp = self.theta;
            }

            let t = self.stepper.current_time();
            let window_kin = self.omega_provider.kinematics_over_step(t, dt);
            let omega_step = window_kin.omega_step;
            let alpha_step = window_kin.alpha_step;
            let theta_target = self.theta + window_kin.theta_increment;

            // ── Read forces from preCICE (global frame) ───────────────────────
            let mut forces_global = vec![0.0f64; n_data];
            participant.read_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.read_data,
                &vertex_ids,
                dt,
                &mut forces_global,
            )?;

            // ── Force pre-processing (ramp + cap, global frame) ───────────────
            apply_ramp(&mut forces_global, t, self.config.fsi.ramp_time);
            if let Some(max_f) = self.config.fsi.force_max {
                apply_cap(&mut forces_global, max_f, mesh_dims);
            }

            // ── Assemble global reduced force vector ──────────────────────────
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

            // ── Advance structural state ──────────────────────────────────────
            let step_t = self.stepper.step(&f_red, dt)?.t;

            // ── Gather interface displacements (elastic or total) ─────────────
            let u_red = self.stepper.current_u();
            let disp_iface = match self.config.displacement_mode {
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
                    // Write total displacement (u_elastic + u_rigid)
                    // Compute rigid-body displacement: u_rigid = rotated_coords - ref_coords
                    let mut disp_total = Vec::with_capacity(self.iface_dofs.len());
                    for &node in &self.iface_nodes {
                        let i_ref = node * 3;
                        let i_rot = node * 3;
                        for j in 0..3 {
                            let u_rigid = self.coords_rotated[i_rot + j] - self.coords_ref[i_ref + j];
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
            participant.write_data(
                &self.config.fsi.coupling_mesh,
                &self.config.fsi.write_data,
                &vertex_ids,
                &disp_iface,
            )?;

            // ── Write ω to GlobalSolidMesh (if used) ──────────────────────────
            if let (Some(ref ids), Some(ref mesh), Some(ref wdata)) = (
                &omega_vertex_ids,
                &self.config.omega_mesh_name,
                &self.config.omega_write_data,
            ) {
                participant.write_data(mesh, wdata, ids, &[omega_step])?;
            }

            participant.advance(dt)?;

            // ── Implicit coupling: restore or commit ──────────────────────────
            if participant.requires_reading_checkpoint()? {
                match (newmark_cp.as_ref(), omega_cp.as_ref()) {
                    (Some(ncp), Some(ocp)) => {
                        self.stepper.restore(ncp);
                        self.omega_provider.restore(ocp);
                        self.theta = theta_cp;
                        // Note: coords_rotated will be recomputed in next K rebuild
                    }
                    _ => {
                        return Err(FsiError::PreciceError(
                            "preCICE requires checkpoint restore but no checkpoint was saved"
                                .to_string(),
                        ))
                    }
                }
            } else {
                // ── Converged time window ─────────────────────────────────────
                self.time_step += 1;
                self.theta = theta_target;

                // Update K(θ) if rotation changed significantly
                self.reassemble_k_if_needed(omega_step, self.time_step, &free_dofs)?;

                // Compute aerodynamic torque (global frame)
                let u_full = self.expand_to_full(u_red, &free_dofs);
                let disp_iface_full: Vec<f64> = self.iface_nodes.iter()
                    .flat_map(|&node| {
                        (0..3).map(move |j| {
                            let dof = node * self.config.dofs_per_node + j;
                            if dof < u_full.len() {
                                u_full[dof]
                            } else {
                                0.0
                            }
                        })
                    })
                    .collect();

                let (tau_vec, tau_aero) = compute_torque(
                    &iface_coords,
                    &disp_iface_full,
                    &forces_global,
                    &self.transforms.axis,
                    &self.transforms.center,
                );

                // Update omega provider if using Computed mode
                if let OmegaProvider::Computed { .. } = &mut self.omega_provider {
                    self.omega_provider.update_from_torque(tau_aero, dt, step_t);
                }

                // Store final state
                result.u_final = u_red.to_vec();
                result.v_final = self.stepper.current_v().to_vec();
                result.a_final = self.stepper.current_a().to_vec();
                result.times.push(step_t);

                // Performance coefficients
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
                    "InertialRotorFsi step {}: t={:.4} θ={:.4}rad ω_window={:.4}rad/s ω_state={:.4}rad/s \
                     τ_aero={:.3e}N·m  Ct={:.4} Cp={:.4} TSR={:.3}",
                    self.time_step, step_t, self.theta, omega_step, omega_state,
                    tau_aero, perf.ct, perf.cp, perf.tsr
                );

                // Per-step callback
                if let Some(ref cb) = self.step_callback {
                    let force_mag = forces_global.iter().map(|x| x * x).sum::<f64>().sqrt();
                    cb(
                        step_t,
                        self.time_step,
                        dt,
                        u_red,
                        self.stepper.current_v(),
                        self.stepper.current_a(),
                        force_mag,
                        &forces_global,
                        omega_state,
                        alpha_state,
                        self.theta,
                        tau_aero,
                        perf,
                    )?;
                }

                dt = participant.get_max_time_step_size()?;
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
