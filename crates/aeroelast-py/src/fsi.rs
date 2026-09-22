//! FSI runner bindings — only compiled with the `fsi` feature.

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::prelude::*;

use crate::assembler::PyMeshAssembler;

// ============================================================================
// FSI binding — only compiled with feature "fsi"
// ============================================================================

/// Run a linear-elastic FSI simulation using preCICE for coupling.
///
/// This function drives the full partitioned FSI loop:
/// registers the structural mesh with preCICE, reads forces,
/// advances the Newmark integrator, writes displacements, and
/// iterates until the coupling is complete.
///
/// The GIL is released for the duration of the coupling loop so that
/// Python threads can run concurrently (though in practice preCICE is
/// single-threaded on the structural side).
///
/// Args:
///   k_rows, k_cols, k_vals   – COO stiffness matrix (i32 indices, n_dofs × n_dofs)
///   m_rows, m_cols, m_vals   – COO mass matrix (same sparsity as K)
///   c_rows, c_cols, c_vals   – COO damping matrix (same sparsity as K; can be all zeros)
///   n_dofs                   – number of free DOFs
///   beta, gamma, dt          – Newmark-β parameters and initial time step
///   interface_coords         – flat array [x0,y0,z0, x1,y1,z1, …] of interface node coords
///   interface_dofs           – global DOF indices of interface nodes in the reduced system
///   mesh_dims                – spatial dimension of the interface (2 or 3)
///   participant_name         – preCICE participant name
///   config_file              – path to precice-config.xml
///   coupling_mesh            – mesh name in preCICE
///   write_data_name          – displacement data field name written to preCICE
///   read_data_name           – force data field name read from preCICE
///   ramp_time                – force ramp-up duration (0.0 = no ramp)
///   force_max                – optional per-component force cap (None = no cap)
///
/// Returns:
///   (u_final, v_final, a_final, times)
///   u_final/v_final/a_final are flat arrays of the last converged step; times is a flat list.
#[pyfunction]
#[cfg(feature = "fsi")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_linear_elastic_fsi(
    py: Python<'_>,
    // Stiffness matrix K (COO)
    k_rows: PyReadonlyArray1<i32>,
    k_cols: PyReadonlyArray1<i32>,
    k_vals: PyReadonlyArray1<f64>,
    // Mass matrix M (COO)
    m_rows: PyReadonlyArray1<i32>,
    m_cols: PyReadonlyArray1<i32>,
    m_vals: PyReadonlyArray1<f64>,
    // Damping matrix C (COO)
    c_rows: PyReadonlyArray1<i32>,
    c_cols: PyReadonlyArray1<i32>,
    c_vals: PyReadonlyArray1<f64>,
    // System size
    n_dofs: usize,
    // Newmark parameters
    beta: f64,
    gamma: f64,
    dt: f64,
    // Interface definition
    interface_coords: PyReadonlyArray1<f64>,
    interface_dofs: PyReadonlyArray1<usize>,
    mesh_dims: usize,
    // preCICE configuration
    participant_name: &str,
    config_file: &str,
    coupling_mesh: &str,
    write_data_name: &str,
    read_data_name: &str,
    ramp_time: f64,
    force_max: Option<f64>,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    use aeroelast_solvers::petsc::elasticity::dynamic_newmark::NewmarkStepper;
    use aeroelast_solvers::petsc::fsi::linear_elastic::{FsiConfig, LinearElasticFsiSolver};
    use pyo3::exceptions::PyRuntimeError;

    // Validate and extract array slices
    let kr = k_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kc = k_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kv = k_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mr = m_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mc = m_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mv = m_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let cr = c_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let cc = c_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let cv = c_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let icoords: Vec<f64> = interface_coords.as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
    let idofs: Vec<usize> = interface_dofs.as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();

    // Owned String copies for the closure (cannot borrow &str across GIL release)
    let participant_name = participant_name.to_string();
    let config_file = config_file.to_string();
    let coupling_mesh = coupling_mesh.to_string();
    let write_data_name = write_data_name.to_string();
    let read_data_name = read_data_name.to_string();

    // Build the Newmark stepper — must happen before GIL release (PETSc init)
    let stepper = NewmarkStepper::new(
        kr, kc, kv,
        mr, mc, mv,
        cr, cc, cv,
        n_dofs, beta, gamma, dt,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let config = FsiConfig {
        participant_name,
        config_file,
        coupling_mesh,
        write_data: write_data_name,
        read_data: read_data_name,
        ramp_time,
        force_max,
    };

    let mut solver = LinearElasticFsiSolver::new(stepper, config, icoords, idofs, mesh_dims);

    let result = solver
        .run()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok((
        result.u_final,
        result.v_final,
        result.a_final,
        result.times,
    ))
}

/// High-level FSI entry point — Rust internalises setup.
///
/// Unlike `run_linear_elastic_fsi`, this binding accepts the **global** (not
/// already-reduced) stiffness and mass matrices, plus a `free_dofs` array.
/// Rust performs:
///   1. BC reduction (extract free-DOF sub-system from global COO)
///   2. Row-sum mass lumping
///   3. Rayleigh damping  C = η_k · K_red + η_m · M_lump
///   4. Interface DOF remapping from global to reduced indices
///   5. Full preCICE coupling loop via [`LinearElasticFsiSolver`]
///
/// # Arguments
///
///   k_rows / k_cols / k_vals  – global stiffness K (COO)
///   m_rows / m_cols / m_vals  – global consistent mass M (COO)
///   free_dofs                 – sorted free (unconstrained) global DOF indices
///   eta_k                     – stiffness-proportional Rayleigh coefficient
///   eta_m                     – mass-proportional Rayleigh coefficient
///   beta / gamma / dt         – Newmark-β parameters
///   interface_coords          – flat coupling-mesh node coordinates
///   interface_dofs_global     – DOF indices in the **global** system
///   mesh_dims                 – spatial dimension of coupling mesh (2 or 3)
///   participant_name / config_file / coupling_mesh / write_data / read_data
///                             – preCICE configuration
///   ramp_time                 – force ramp-up duration (0.0 = no ramp)
///   force_max                 – optional per-node force cap (None = no cap)
///   u0 / v0 / a0              – optional initial displacement / velocity /
///                               acceleration in the **reduced** DOF space
///                               (pass None for a cold start)
///   t0                        – simulation start time (0.0 for cold start)
///
/// # Returns
///   (u_final, v_final, a_final, times)
///   Flat reduced-DOF arrays for the last converged step; times accumulates all step times.
#[pyfunction]
#[cfg(feature = "fsi")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_fsi_solver(
    py: Python<'_>,
    // Global stiffness K (COO)
    k_rows: PyReadonlyArray1<i32>,
    k_cols: PyReadonlyArray1<i32>,
    k_vals: PyReadonlyArray1<f64>,
    // Global consistent mass M (COO)
    m_rows: PyReadonlyArray1<i32>,
    m_cols: PyReadonlyArray1<i32>,
    m_vals: PyReadonlyArray1<f64>,
    // BC reduction
    free_dofs: PyReadonlyArray1<i32>,
    // Rayleigh damping
    eta_k: f64,
    eta_m: f64,
    // Newmark parameters
    beta: f64,
    gamma: f64,
    dt: f64,
    // Interface definition (global DOF numbering)
    interface_coords: PyReadonlyArray1<f64>,
    interface_dofs_global: PyReadonlyArray1<usize>,
    mesh_dims: usize,
    // preCICE configuration
    participant_name: &str,
    config_file: &str,
    coupling_mesh: &str,
    write_data_name: &str,
    read_data_name: &str,
    ramp_time: f64,
    force_max: Option<f64>,
    // Optional restart state (reduced DOF space)
    u0: Option<PyReadonlyArray1<f64>>,
    v0: Option<PyReadonlyArray1<f64>>,
    a0: Option<PyReadonlyArray1<f64>>,
    t0: f64,
    // Optional per-step callback: callable(t, time_step, dt, u_red, v_red, a_red)
    step_callback: Option<Py<PyAny>>,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    use aeroelast_solvers::petsc::elasticity::dynamic_newmark::NewmarkStepper;
    use aeroelast_solvers::petsc::fsi::linear_elastic::{
        FsiConfig, FsiInitialState, LinearElasticFsiSolver,
    };
    use aeroelast_solvers::petsc::fsi::setup;
    use pyo3::exceptions::PyRuntimeError;

    // Extract input slices.
    let kr = k_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kc = k_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kv = k_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mr = m_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mc = m_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mv = m_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let fd = free_dofs.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let icoords: Vec<f64> = interface_coords
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();
    let idofs_global: Vec<usize> = interface_dofs_global
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();

    // ── 1. BC reduction ───────────────────────────────────────────────────────
    let (kr_red, kc_red, kv_red) = setup::reduce_coo(kr, kc, kv, fd);
    let (_, _, mv_red) = setup::reduce_coo(mr, mc, mv, fd);
    let n_free = fd.len();

    // ── 2. Extract lumped diagonal from the reduced M ─────────────────────────
    // Python always passes M already row-sum lumped (diagonal COO).
    // After reduce_coo the diagonal entries are preserved.  We collect
    // them into a plain f64 vector indexed by reduced DOF.
    let mut m_diag = vec![0.0f64; n_free];
    for (row, &v) in (0..).zip(mv_red.iter()) {
        // The reduced M diagonal has exactly n_free entries in sorted order.
        if row < n_free {
            m_diag[row] += v.abs();
        }
    }

    // ── 3. Expand diagonal M into K's sparsity (NewmarkStepper requirement) ───
    // NewmarkStepper requires K, M, C to share the SAME sparsity pattern.
    // Expanding M_lump (diagonal) to K's full sparsity inserts zeros at
    // off-diagonal positions, which is semantically correct for lumped M.
    let mv_expanded = setup::expand_diag_to_sparsity(&m_diag, &kr_red, &kc_red);

    // ── 4. Rayleigh damping C = η_k·K + η_m·M (in K's sparsity) ─────────────
    let cv_red: Vec<f64> = kv_red
        .iter()
        .zip(mv_expanded.iter())
        .map(|(&k, &m)| eta_k * k + eta_m * m)
        .collect();

    // ── 5. Interface DOF remapping (global → reduced) ─────────────────────────
    let idofs_red: Vec<usize> = setup::remap_interface_dofs(&idofs_global, fd);

    // ── 6. Owned String copies (cannot borrow &str across potential GIL ops) ──
    let participant_name = participant_name.to_string();
    let config_file = config_file.to_string();
    let coupling_mesh = coupling_mesh.to_string();
    let write_data_name = write_data_name.to_string();
    let read_data_name = read_data_name.to_string();

    // ── 7. Build Newmark stepper ──────────────────────────────────────────────
    // K, M_expanded, and C all share the same sparsity (kr_red/kc_red).
    let stepper = NewmarkStepper::new(
        &kr_red, &kc_red, &kv_red,
        &kr_red, &kc_red, &mv_expanded,
        &kr_red, &kc_red, &cv_red,
        n_free, beta, gamma, dt,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let config = FsiConfig {
        participant_name,
        config_file,
        coupling_mesh,
        write_data: write_data_name,
        read_data: read_data_name,
        ramp_time,
        force_max,
    };

    // ── 8. Optionally extract initial state for restart ───────────────────────
    let initial_state = if let (Some(u_arr), Some(v_arr), Some(a_arr)) = (u0, v0, a0) {
        let u = u_arr
            .as_slice()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .to_vec();
        let v = v_arr
            .as_slice()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .to_vec();
        let a = a_arr
            .as_slice()
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
            .to_vec();
        Some(FsiInitialState { u, v, a, t: t0 })
    } else {
        None
    };

    // ── 9. Run solver ─────────────────────────────────────────────────────────
    let mut solver = {
        let s = LinearElasticFsiSolver::new(stepper, config, icoords, idofs_red, mesh_dims);
        let s = if let Some(state) = initial_state {
            s.with_initial_state(state)
        } else {
            s
        };
        if let Some(py_cb) = step_callback {
            s.with_step_callback(move |t, time_step, dt, u, v, a, force_mag, forces_iface| {
                Python::attach(|py| {
                    let u_arr = Array1::from(u.to_vec()).into_pyarray(py);
                    let v_arr = Array1::from(v.to_vec()).into_pyarray(py);
                    let a_arr = Array1::from(a.to_vec()).into_pyarray(py);
                    let fi_arr = Array1::from(forces_iface.to_vec()).into_pyarray(py);
                    py_cb
                        .call1(py, (t, time_step as u64, dt, u_arr, v_arr, a_arr, force_mag, fi_arr))
                        .map(|_| ())
                        .map_err(|e| {
                            aeroelast_solvers::petsc::fsi::linear_elastic::FsiError::CallbackError(
                                e.to_string(),
                            )
                        })
                })
            })
        } else {
            s
        }
    };

    let _ = py; // GIL is held; preCICE may release internally via C
    let result = solver
        .run()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok((
        result.u_final,
        result.v_final,
        result.a_final,
        result.times,
    ))
}

// ── run_stress_stiffened_fsi_solver ───────────────────────────────────────────
/// Identical to `run_fsi_solver` but uses `StressStiffenedFsiSolver` which rebuilds
/// the geometric stiffness K_G every `kg_update_interval` converged steps.
///
/// Extra arguments compared to `run_fsi_solver`:
///   assembler          – `PyMeshAssembler` instance (for stress recovery and K_G assembly)
///   n_full_dofs        – total DOF count in the unreduced system
///   kg_update_interval – rebuild K_G every N converged steps (1 = every step)
#[pyfunction]
#[cfg(feature = "fsi")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_stress_stiffened_fsi_solver(
    py: Python<'_>,
    // Stress-stiffening extras
    assembler: &PyMeshAssembler,
    n_full_dofs: usize,
    kg_update_interval: usize,
    // Global stiffness K (COO)
    k_rows: PyReadonlyArray1<i32>,
    k_cols: PyReadonlyArray1<i32>,
    k_vals: PyReadonlyArray1<f64>,
    // Global consistent mass M (COO)
    m_rows: PyReadonlyArray1<i32>,
    m_cols: PyReadonlyArray1<i32>,
    m_vals: PyReadonlyArray1<f64>,
    // BC reduction
    free_dofs: PyReadonlyArray1<i32>,
    // Rayleigh damping
    eta_k: f64,
    eta_m: f64,
    // Newmark parameters
    beta: f64,
    gamma: f64,
    dt: f64,
    // Interface definition (global DOF numbering)
    interface_coords: PyReadonlyArray1<f64>,
    interface_dofs_global: PyReadonlyArray1<usize>,
    mesh_dims: usize,
    // preCICE configuration
    participant_name: &str,
    config_file: &str,
    coupling_mesh: &str,
    write_data_name: &str,
    read_data_name: &str,
    ramp_time: f64,
    force_max: Option<f64>,
    // Optional restart state (reduced DOF space)
    u0: Option<PyReadonlyArray1<f64>>,
    v0: Option<PyReadonlyArray1<f64>>,
    a0: Option<PyReadonlyArray1<f64>>,
    t0: f64,
    // Optional per-step callback
    step_callback: Option<Py<PyAny>>,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    use aeroelast_solvers::petsc::elasticity::dynamic_newmark::NewmarkStepper;
    use aeroelast_solvers::petsc::fsi::linear_elastic::{FsiConfig, FsiInitialState};
    use aeroelast_solvers::petsc::fsi::stress_stiffened::StressStiffenedFsiSolver;
    use aeroelast_solvers::petsc::fsi::setup;
    use pyo3::exceptions::PyRuntimeError;

    let kr = k_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kc = k_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kv = k_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mr = m_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mc = m_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mv = m_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let fd = free_dofs.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let icoords: Vec<f64> = interface_coords
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();
    let idofs_global: Vec<usize> = interface_dofs_global
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();

    // ── BC reduction ──────────────────────────────────────────────────────────
    let (kr_red, kc_red, kv_red) = setup::reduce_coo(kr, kc, kv, fd);
    let (_, _, mv_red) = setup::reduce_coo(mr, mc, mv, fd);
    let n_free = fd.len();

    let mut m_diag = vec![0.0f64; n_free];
    for (row, &v) in (0..).zip(mv_red.iter()) {
        if row < n_free {
            m_diag[row] += v.abs();
        }
    }
    let mv_expanded = setup::expand_diag_to_sparsity(&m_diag, &kr_red, &kc_red);

    let cv_red: Vec<f64> = kv_red
        .iter()
        .zip(mv_expanded.iter())
        .map(|(&k, &m)| eta_k * k + eta_m * m)
        .collect();

    let idofs_red: Vec<usize> = setup::remap_interface_dofs(&idofs_global, fd);

    let participant_name = participant_name.to_string();
    let config_file = config_file.to_string();
    let coupling_mesh = coupling_mesh.to_string();
    let write_data_name = write_data_name.to_string();
    let read_data_name = read_data_name.to_string();

    let stepper = NewmarkStepper::new(
        &kr_red, &kc_red, &kv_red,
        &kr_red, &kc_red, &mv_expanded,
        &kr_red, &kc_red, &cv_red,
        n_free, beta, gamma, dt,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let config = FsiConfig {
        participant_name,
        config_file,
        coupling_mesh,
        write_data: write_data_name,
        read_data: read_data_name,
        ramp_time,
        force_max,
    };

    let initial_state = if let (Some(u_arr), Some(v_arr), Some(a_arr)) = (u0, v0, a0) {
        let u = u_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let v = v_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let a = a_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        Some(FsiInitialState { u, v, a, t: t0 })
    } else {
        None
    };

    // Clone the inner assembler (MeshAssembler: Clone) for the Rust solver.
    let rust_assembler = assembler.inner().clone();

    let mut solver = {
        let s = StressStiffenedFsiSolver::new(
            stepper,
            config,
            icoords,
            idofs_red,
            mesh_dims,
            rust_assembler,
            fd.to_vec(),
            n_full_dofs,
            kg_update_interval,
        );
        let s = if let Some(state) = initial_state { s.with_initial_state(state) } else { s };
        if let Some(py_cb) = step_callback {
            s.with_step_callback(move |t, time_step, dt, u, v, a, force_mag, forces_iface| {
                Python::attach(|py| {
                    let u_arr = Array1::from(u.to_vec()).into_pyarray(py);
                    let v_arr = Array1::from(v.to_vec()).into_pyarray(py);
                    let a_arr = Array1::from(a.to_vec()).into_pyarray(py);
                    let fi_arr = Array1::from(forces_iface.to_vec()).into_pyarray(py);
                    py_cb
                        .call1(py, (t, time_step as u64, dt, u_arr, v_arr, a_arr, force_mag, fi_arr))
                        .map(|_| ())
                        .map_err(|e| {
                            aeroelast_solvers::petsc::fsi::linear_elastic::FsiError::CallbackError(
                                e.to_string(),
                            )
                        })
                })
            })
        } else {
            s
        }
    };

    let _ = py;
    let result = solver
        .run()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok((
        result.u_final,
        result.v_final,
        result.a_final,
        result.times,
    ))
}

// ── run_rotor_fsi_solver ──────────────────────────────────────────────────────
/// Co-rotational FSI solver for rotating structures (rotor blades).
///
/// Wraps `RotorFsiSolver` with the full preCICE coupling loop including:
/// - Rodrigues coordinate transforms (forces global→rotating, disps rotating→global)
/// - Centrifugal, Coriolis, Euler, and gravity inertial body forces
/// - Spin-softening K_SP diagonal update when |Δω| > `ksp_omega_threshold`
/// - Geometric stiffness K_G rebuild every `kg_update_interval` converged steps
///   (`0` and `1` both mean every converged step for backward compatibility)
/// - Four angular-velocity modes: "constant", "ramped", "computed", "ramped_computed"
///
/// # Returns
/// `(u_final, v_final, a_final, times)`
/// Flat reduced-DOF arrays for the last converged step; times accumulates all step times.
#[pyfunction]
#[cfg(feature = "fsi")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_rotor_fsi_solver(
    py: Python<'_>,
    // ── Optional assembler (K_G) ──────────────────────────────────────────────
    assembler: Option<&PyMeshAssembler>,
    n_full_dofs: usize,
    kg_update_interval: usize,
    // ── Rotor geometry ────────────────────────────────────────────────────────
    rotation_axis: Vec<f64>,
    rotation_center: Vec<f64>,
    // ── All-node data (full structural mesh, rotating frame) ──────────────────
    all_node_coords: PyReadonlyArray1<f64>,
    all_node_masses: PyReadonlyArray1<f64>,
    // ── Angular-velocity provider ─────────────────────────────────────────────
    omega_mode: &str,
    omega: f64,
    omega_target: Option<f64>,
    t_ramp: Option<f64>,
    moment_of_inertia: Option<f64>,
    shaft_torque: Option<f64>,
    // ── Inertial forces ───────────────────────────────────────────────────────
    gravity: Vec<f64>,
    include_centrifugal: bool,
    include_coriolis: bool,
    include_euler: bool,
    // ── Stiffness updates ─────────────────────────────────────────────────────
    include_kg: bool,
    include_ksp: bool,
    ksp_omega_threshold: f64,
    // ── DOF layout ────────────────────────────────────────────────────────────
    dofs_per_node: usize,
    // ── Performance coefficients ──────────────────────────────────────────────
    fluid_density: f64,
    flow_velocity: f64,
    rotor_radius: f64,
    // ── Global stiffness K (COO) ──────────────────────────────────────────────
    k_rows: PyReadonlyArray1<i32>,
    k_cols: PyReadonlyArray1<i32>,
    k_vals: PyReadonlyArray1<f64>,
    // ── Global consistent mass M (COO) ────────────────────────────────────────
    m_rows: PyReadonlyArray1<i32>,
    m_cols: PyReadonlyArray1<i32>,
    m_vals: PyReadonlyArray1<f64>,
    // ── BC reduction ──────────────────────────────────────────────────────────
    free_dofs: PyReadonlyArray1<i32>,
    // ── Rayleigh damping ──────────────────────────────────────────────────────
    eta_k: f64,
    eta_m: f64,
    // ── Newmark parameters ────────────────────────────────────────────────────
    beta: f64,
    gamma: f64,
    dt: f64,
    // ── Interface (global DOF numbering) ──────────────────────────────────────
    interface_coords: PyReadonlyArray1<f64>,
    interface_dofs_global: PyReadonlyArray1<usize>,
    mesh_dims: usize,
    // ── preCICE configuration ─────────────────────────────────────────────────
    participant_name: &str,
    config_file: &str,
    coupling_mesh: &str,
    write_data_name: &str,
    read_data_name: &str,
    ramp_time: f64,
    force_max: Option<f64>,
    // ── Optional GlobalSolidMesh for ω communication ──────────────────────────
    omega_mesh_name: Option<String>,
    omega_write_data: Option<String>,
    omega_vertex_coord: Option<Vec<f64>>,
    // ── Optional restart state (reduced DOF space) ────────────────────────────
    u0: Option<PyReadonlyArray1<f64>>,
    v0: Option<PyReadonlyArray1<f64>>,
    a0: Option<PyReadonlyArray1<f64>>,
    t0: f64,
    theta0: f64,
    restart_omega: Option<f64>,
    restart_alpha: Option<f64>,
    restart_ramp_completed: Option<bool>,
    restart_current_time: Option<f64>,
    // ── Initial geometric stiffness K_G (centrifugal prestress, full DOF space) ─
    kg0_rows: Option<PyReadonlyArray1<i64>>,
    kg0_cols: Option<PyReadonlyArray1<i64>>,
    kg0_vals: Option<PyReadonlyArray1<f64>>,
    // ── Optional per-step callback ────────────────────────────────────────────
    step_callback: Option<Py<PyAny>>,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    use aeroelast_solvers::petsc::elasticity::dynamic_newmark::NewmarkStepper;
    use aeroelast_solvers::petsc::fsi::linear_elastic::{FsiConfig, FsiInitialState};
    use aeroelast_solvers::petsc::fsi::rotor_fsi::{RotorFsiConfig, RotorFsiSolver};
    use aeroelast_solvers::petsc::fsi::rotor_physics::OmegaProvider;
    use aeroelast_solvers::petsc::fsi::setup;
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::exceptions::PyValueError;

    let kr = k_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kc = k_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kv = k_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mr = m_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mc = m_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mv = m_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let fd = free_dofs.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let icoords: Vec<f64> = interface_coords
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();
    let idofs_global: Vec<usize> = interface_dofs_global
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();
    let node_coords: Vec<f64> = all_node_coords
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();
    let node_masses: Vec<f64> = all_node_masses
        .as_slice()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
        .to_vec();

    // ── Validate geometry vectors ─────────────────────────────────────────────
    if rotation_axis.len() != 3 || rotation_center.len() != 3 {
        return Err(PyValueError::new_err(
            "rotation_axis and rotation_center must have length 3",
        ));
    }
    if gravity.len() != 3 {
        return Err(PyValueError::new_err("gravity must have length 3"));
    }

    // ── BC reduction ──────────────────────────────────────────────────────────
    let (kr_red, kc_red, kv_red) = setup::reduce_coo(kr, kc, kv, fd);
    let (_, _, mv_red) = setup::reduce_coo(mr, mc, mv, fd);
    let n_free = fd.len();

    let mut m_diag = vec![0.0f64; n_free];
    for (row, &v) in (0..).zip(mv_red.iter()) {
        if row < n_free {
            m_diag[row] += v.abs();
        }
    }
    let mv_expanded = setup::expand_diag_to_sparsity(&m_diag, &kr_red, &kc_red);

    let cv_red: Vec<f64> = kv_red
        .iter()
        .zip(mv_expanded.iter())
        .map(|(&k, &m)| eta_k * k + eta_m * m)
        .collect();

    let idofs_red: Vec<usize> = setup::remap_interface_dofs(&idofs_global, fd);

    // ── m_lumped_full: lump M COO over the full (unreduced) DOF space ─────────
    let mut m_lumped_full = vec![0.0f64; n_full_dofs];
    for (&row, &val) in mr.iter().zip(mv.iter()) {
        if row >= 0 && (row as usize) < n_full_dofs {
            m_lumped_full[row as usize] += val.abs();
        }
    }

    // ── Build OmegaProvider ───────────────────────────────────────────────────
    let omega_provider = match omega_mode {
        "constant" => OmegaProvider::Constant { omega },
        "ramped" => OmegaProvider::Ramped {
            omega_target: omega_target.ok_or_else(|| {
                PyValueError::new_err("omega_mode='ramped' requires omega_target")
            })?,
            t_ramp: t_ramp.ok_or_else(|| {
                PyValueError::new_err("omega_mode='ramped' requires t_ramp")
            })?,
        },
        "computed" => OmegaProvider::Computed {
            moment_of_inertia: moment_of_inertia.ok_or_else(|| {
                PyValueError::new_err("omega_mode='computed' requires moment_of_inertia")
            })?,
            shaft_torque: shaft_torque.unwrap_or(0.0),
            omega: restart_omega.unwrap_or(omega),
            alpha: restart_alpha.unwrap_or(0.0),
        },
        "ramped_computed" => {
            let omega_target = omega_target.ok_or_else(|| {
                PyValueError::new_err("omega_mode='ramped_computed' requires omega_target")
            })?;
            let t_ramp = t_ramp.ok_or_else(|| {
                PyValueError::new_err("omega_mode='ramped_computed' requires t_ramp")
            })?;
            let current_time = restart_current_time.unwrap_or(t0);
            let ramp_completed =
                restart_ramp_completed.unwrap_or(current_time >= t_ramp - 1.0e-12);
            let omega_state = if ramp_completed {
                restart_omega.unwrap_or_else(|| if omega.abs() > 0.0 { omega } else { omega_target })
            } else if let Some(omega_restart) = restart_omega {
                omega_restart
            } else if current_time < t_ramp {
                omega_target * current_time / t_ramp
            } else {
                omega_target
            };
            let alpha_state = if ramp_completed {
                restart_alpha.unwrap_or(0.0)
            } else if current_time < t_ramp {
                omega_target / t_ramp
            } else {
                restart_alpha.unwrap_or(0.0)
            };

            OmegaProvider::RampedComputed {
                omega_target,
                t_ramp,
                moment_of_inertia: moment_of_inertia.ok_or_else(|| {
                PyValueError::new_err("omega_mode='ramped_computed' requires moment_of_inertia")
                })?,
                shaft_torque: shaft_torque.unwrap_or(0.0),
                omega: omega_state,
                alpha: alpha_state,
                ramp_completed,
                current_time,
            }
        }
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown omega_mode '{}'; expected 'constant', 'ramped', 'computed', or 'ramped_computed'",
                other
            )))
        }
    };

    // ── Geometry arrays ───────────────────────────────────────────────────────
    let axis: [f64; 3] = rotation_axis[..3]
        .try_into()
        .map_err(|_| PyValueError::new_err("rotation_axis must have exactly 3 elements"))?;
    let center: [f64; 3] = rotation_center[..3]
        .try_into()
        .map_err(|_| PyValueError::new_err("rotation_center must have exactly 3 elements"))?;
    let grav: [f64; 3] = gravity[..3]
        .try_into()
        .map_err(|_| PyValueError::new_err("gravity must have exactly 3 elements"))?;

    let omega_vertex_coord_arr: Option<[f64; 3]> = omega_vertex_coord
        .as_ref()
        .map(|v| {
            if v.len() != 3 {
                return Err(PyValueError::new_err("omega_vertex_coord must have length 3"));
            }
            Ok([v[0], v[1], v[2]])
        })
        .transpose()?;

    // ── preCICE strings ───────────────────────────────────────────────────────
    let participant_name = participant_name.to_string();
    let config_file = config_file.to_string();
    let coupling_mesh = coupling_mesh.to_string();
    let write_data_name = write_data_name.to_string();
    let read_data_name = read_data_name.to_string();

    // ── Newmark stepper ───────────────────────────────────────────────────────
    let mut stepper = NewmarkStepper::new(
        &kr_red, &kc_red, &kv_red,
        &kr_red, &kc_red, &mv_expanded,
        &kr_red, &kc_red, &cv_red,
        n_free, beta, gamma, dt,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // ── Initial geometric stiffness (centrifugal prestress) ───────────────────
    if let (Some(kg0r), Some(kg0c), Some(kg0v)) = (kg0_rows, kg0_cols, kg0_vals) {
        let kg0r_s = kg0r.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let kg0c_s = kg0c.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let kg0v_s = kg0v.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let kg_coo_map = setup::build_kg_coo_map(kg0r_s, kg0c_s, fd, &kr_red, &kc_red);
        let kg0_red = setup::apply_kg_coo_map(&kg_coo_map, kg0v_s, kr_red.len());
        stepper
            .set_initial_geometric_stiffness(&kg0_red)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    }

    // ── FSI base config ───────────────────────────────────────────────────────
    let fsi_config = FsiConfig {
        participant_name,
        config_file,
        coupling_mesh,
        write_data: write_data_name,
        read_data: read_data_name,
        ramp_time,
        force_max,
    };

    // ── Rotor FSI config ──────────────────────────────────────────────────────
    let config = RotorFsiConfig {
        fsi: fsi_config,
        rotation_axis: axis,
        rotation_center: center,
        gravity: grav,
        include_centrifugal,
        include_coriolis,
        include_euler,
        include_kg,
        kg_update_interval,
        include_ksp,
        ksp_omega_threshold,
        dofs_per_node,
        fluid_density,
        flow_velocity,
        rotor_radius,
        omega_mesh_name,
        omega_write_data,
        omega_vertex_coord: omega_vertex_coord_arr,
    };

    // ── Initial state ─────────────────────────────────────────────────────────
    let initial_state = if let (Some(u_arr), Some(v_arr), Some(a_arr)) = (u0, v0, a0) {
        let u = u_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let v = v_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let a = a_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        Some(FsiInitialState { u, v, a, t: t0 })
    } else {
        None
    };

    // ── Build and configure solver ────────────────────────────────────────────
    let rust_assembler = assembler.map(|a| a.inner().clone());
    let mut solver = {
        let s = RotorFsiSolver::new(
            stepper,
            config,
            icoords,
            idofs_red,
            mesh_dims,
            omega_provider,
            node_coords,
            node_masses,
            m_lumped_full,
            kr_red,
            kc_red,
            fd.to_vec(),
            n_full_dofs,
            rust_assembler,
        );
        let s = if let Some(state) = initial_state {
            s.with_initial_state(state, theta0)
        } else {
            s
        };
        if let Some(py_cb) = step_callback {
            s.with_step_callback(move |t, step, dt, u, v, a, force_mag, forces_iface,
                                        omega, alpha, theta, tau_aero, perf| {
                Python::attach(|py| {
                    let u_arr = Array1::from(u.to_vec()).into_pyarray(py);
                    let v_arr = Array1::from(v.to_vec()).into_pyarray(py);
                    let a_arr = Array1::from(a.to_vec()).into_pyarray(py);
                    let fi_arr = Array1::from(forces_iface.to_vec()).into_pyarray(py);
                    // PyO3 tuple conversion is limited to 12 elements.
                    // Pack performance coefficients as a sub-tuple to stay within the limit.
                    // Python callback receives:
                    //   (t, step, dt, u, v, a, force_mag, fi,
                    //    omega, alpha, theta, (tau_aero, ct, cp, cq, tsr))
                    let perf_tuple = (perf.ct, perf.cp, perf.cq, perf.tsr);
                    py_cb
                        .call1(
                            py,
                            (
                                t, step as u64, dt,
                                u_arr, v_arr, a_arr,
                                force_mag, fi_arr,
                                omega, alpha, theta,
                                (tau_aero, perf_tuple.0, perf_tuple.1, perf_tuple.2, perf_tuple.3),
                            ),
                        )
                        .map(|_| ())
                        .map_err(|e| {
                            aeroelast_solvers::petsc::fsi::linear_elastic::FsiError::CallbackError(
                                e.to_string(),
                            )
                        })
                })
            })
        } else {
            s
        }
    };

    let _ = py;
    let result = solver
        .run()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    Ok((
        result.u_final,
        result.v_final,
        result.a_final,
        result.times,
    ))
}
