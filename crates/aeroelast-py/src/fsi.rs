//! FSI runner bindings — only compiled with the `fsi` feature.

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::prelude::*;

use crate::assembler::PyMeshAssembler;

// ============================================================================
// FSI binding — only compiled with feature "fsi"
// ============================================================================


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
    .and_then(|stepper| stepper.with_rayleigh_damping(eta_k, eta_m))
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
    velocity_write_data: Option<String>,
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
    .and_then(|stepper| stepper.with_rayleigh_damping(eta_k, eta_m))
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
            velocity_write_data,
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


#[pyfunction]
#[pyo3(signature = (
    assembler, n_full_dofs, kg_update_interval,
    rotation_axis, rotation_center,
    all_node_coords, all_node_masses,
    omega_mode, omega,
    omega_target, t_ramp, moment_of_inertia, shaft_torque,
    gravity, include_centrifugal, include_coriolis, include_euler,
    include_kg, include_ksp, ksp_omega_threshold,
    omega_rebuild_rel_high, omega_rebuild_rel_low,
    kg_use_deformed_coords, kg_deflection_rebuild_rel_high, kg_deflection_rebuild_rel_low,
    dofs_per_node, fluid_density, flow_velocity, rotor_radius,
    k_rows, k_cols, k_vals,
    m_rows, m_cols, m_vals,
    free_dofs,
    eta_k, eta_m,
    beta, gamma, dt,
    interface_coords, interface_dofs_global, mesh_dims,
    participant_name, config_file, coupling_mesh,
    write_data_name, read_data_name, ramp_time, force_max,
    omega_mesh_name, omega_write_data, omega_vertex_coord,
    u0, v0, a0, t0, theta0,
    restart_omega, restart_alpha, restart_ramp_completed, restart_current_time,
    kg0_rows, kg0_cols, kg0_vals,
    step_callback,
    use_corotational_kt=false,
    kt_coro_update_freq=1
))]
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
    omega_rebuild_rel_high: f64,
    omega_rebuild_rel_low: f64,
    // ── K_G(u) — Nivel 2 deformed-coords K_G ──────────────────────────────────
    kg_use_deformed_coords: bool,
    kg_deflection_rebuild_rel_high: f64,
    kg_deflection_rebuild_rel_low: f64,
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
    // ── Optional corotational K_T wiring (backward-compatible defaults) ──────
    use_corotational_kt: bool,
    kt_coro_update_freq: u32,
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
            alpha_prev: None,
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
                alpha_prev: None,
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
    .and_then(|stepper| stepper.with_rayleigh_damping(eta_k, eta_m))
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
        use_corotational_kt,
        kt_coro_update_freq,
        include_ksp,
        ksp_omega_threshold,
        omega_rebuild_rel_high,
        omega_rebuild_rel_low,
        kg_use_deformed_coords,
        kg_deflection_rebuild_rel_high,
        kg_deflection_rebuild_rel_low,
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
                                        omega, alpha, theta, tau_aero, perf, inertial| {
                Python::attach(|py| {
                    let u_arr = Array1::from(u.to_vec()).into_pyarray(py);
                    let v_arr = Array1::from(v.to_vec()).into_pyarray(py);
                    let a_arr = Array1::from(a.to_vec()).into_pyarray(py);
                    let fi_arr = Array1::from(forces_iface.to_vec()).into_pyarray(py);
                    // PyO3 tuple conversion is limited to 12 elements.
                    // Pack performance coefficients + inertial forces as a sub-tuple.
                    // Fix #5: inertial forces (centrifugal+Euler, rotating frame, per-node)
                    // are appended as the 6th element of the perf sub-tuple so that the
                    // outer 12-element call1 limit is not exceeded.
                    // Python callback receives:
                    //   (t, step, dt, u, v, a, force_mag, fi,
                    //    omega, alpha, theta, (tau_aero, ct, cp, cq, tsr, f_inertial))
                    let inertial_arr = Array1::from(inertial.to_vec()).into_pyarray(py);
                    let perf_tuple = (perf.ct, perf.cp, perf.cq, perf.tsr);
                    py_cb
                        .call1(
                            py,
                            (
                                t, step as u64, dt,
                                u_arr, v_arr, a_arr,
                                force_mag, fi_arr,
                                omega, alpha, theta,
                                (tau_aero, perf_tuple.0, perf_tuple.1, perf_tuple.2, perf_tuple.3,
                                 inertial_arr),
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


#[pyfunction]
#[pyo3(signature = (
    assembler, rotation_axis, rotation_center,
    all_node_masses, omega_mode, omega,
    omega_target, t_ramp, moment_of_inertia, shaft_torque,
    gravity, include_reference_acceleration,
    k_update_interval, omega_rebuild_threshold, theta_rebuild_threshold,
    kg_use_deformed_coords, kg_deflection_rebuild_rel_high, kg_deflection_rebuild_rel_low,
    displacement_mode, dofs_per_node,
    fluid_density, flow_velocity, rotor_radius,
    k_rows, k_cols, k_vals,
    m_rows, m_cols, m_vals,
    free_dofs, eta_k, eta_m,
    beta, gamma, dt,
    interface_nodes, mesh_dims,
    participant_name, config_file, coupling_mesh,
    write_data_name, read_data_name, ramp_time, force_max,
    omega_mesh_name, omega_write_data, omega_vertex_coord,
    velocity_write_data,
    u0, v0, a0, t0, theta0,
    restart_omega, restart_alpha, restart_ramp_completed, restart_current_time,
    kg0_rows, kg0_cols, kg0_vals,
    step_callback,
    include_geometric_stiffness=None,
    include_spin_softening=true,
    ksp_omega_rebuild_high=0.005, ksp_omega_rebuild_low=0.003,
    kg_omega_rebuild_high=0.005, kg_omega_rebuild_low=0.003,
    use_corotational_kt=false, kt_coro_update_freq=1
))]
#[cfg(feature = "fsi")]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_inertial_rotor_fsi_solver(
    py: Python<'_>,
    // ── Assembler (required for K(θ) reassembly) ──────────────────────────────
    assembler: &PyMeshAssembler,
    // ── Rotor geometry ────────────────────────────────────────────────────────
    rotation_axis: Vec<f64>,
    rotation_center: Vec<f64>,
    // ── Nodal data (full structural mesh) ─────────────────────────────────────
    all_node_masses: PyReadonlyArray1<f64>,
    // ── Angular-velocity provider ─────────────────────────────────────────────
    omega_mode: &str,
    omega: f64,
    omega_target: Option<f64>,
    t_ramp: Option<f64>,
    moment_of_inertia: Option<f64>,
    shaft_torque: Option<f64>,
    // ── Inertial-frame physics ────────────────────────────────────────────────
    gravity: Vec<f64>,
    include_reference_acceleration: bool,
    // ── Stiffness reassembly ──────────────────────────────────────────────────
    k_update_interval: usize,
    omega_rebuild_threshold: f64,
    theta_rebuild_threshold: f64,
    // ── K_G(u) — Nivel 2 deformed-coords K_G ──────────────────────────────────
    kg_use_deformed_coords: bool,
    kg_deflection_rebuild_rel_high: f64,
    kg_deflection_rebuild_rel_low: f64,
    // ── Displacement mode ─────────────────────────────────────────────────────
    displacement_mode: &str,
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
    interface_nodes: PyReadonlyArray1<usize>,
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
    // ── Optional nodal velocity write (aerodynamic damping) ───────────────────
    velocity_write_data: Option<String>,
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
    // ── Optional K_G enable sentinel (deprecated — use include_geometric_stiffness) ──
    kg0_rows: Option<PyReadonlyArray1<i64>>,
    kg0_cols: Option<PyReadonlyArray1<i64>>,
    kg0_vals: Option<PyReadonlyArray1<f64>>,
    // ── Optional per-step callback ────────────────────────────────────────────
    step_callback: Option<Py<PyAny>>,
    // ── New keyword-only args with defaults (backward-compatible additions) ────
    include_geometric_stiffness: Option<bool>,
    include_spin_softening: bool,
    ksp_omega_rebuild_high: f64,
    ksp_omega_rebuild_low: f64,
    kg_omega_rebuild_high: f64,
    kg_omega_rebuild_low: f64,
    use_corotational_kt: bool,
    kt_coro_update_freq: u32,
) -> PyResult<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    use aeroelast_solvers::petsc::elasticity::dynamic_newmark::NewmarkStepper;
    use aeroelast_solvers::petsc::fsi::linear_elastic::{FsiConfig, FsiInitialState};
    use aeroelast_solvers::petsc::fsi::rotor_inertial::{
        DisplacementMode, InertialRotorFsiConfig, InertialRotorFsiSolver,
    };
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
    let iface_nodes: Vec<usize> = interface_nodes
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

    // ── Parse displacement mode ───────────────────────────────────────────────
    let disp_mode = match displacement_mode {
        "elastic" => DisplacementMode::Elastic,
        "total" => DisplacementMode::Total,
        other => {
            return Err(PyValueError::new_err(format!(
                "unknown displacement_mode '{}'; expected 'elastic' or 'total'",
                other
            )))
        }
    };

    // ── BC reduction ──────────────────────────────────────────────────────────
    let (kr_red, kc_red, kv_red) = setup::reduce_coo(kr, kc, kv, fd);
    let (_, _, mv_red) = setup::reduce_coo(mr, mc, mv, fd);
    let n_free = fd.len();

    // NOTE: k_coo_map is built later from the Rust assembler's own K sparsity
    // (after rust_assembler.assemble_k()), which is the authoritative source.
    // The Python-side COO (kr/kc as i32) was previously used here to build a
    // first map, but it was immediately shadowed → dead code, now removed.

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

    // ── Use full-mesh nodal masses provided by Python ─────────────────────────
    // The inertial solver computes F_ref on the full structural mesh, so it
    // needs one scalar mass per original node, including constrained nodes.
    // Re-deriving masses from the reduced diagonal would silently drop fixed
    // root nodes and break the n_nodes × 3 reference-acceleration layout.
    let masses = node_masses;

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
            alpha_prev: None,
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
                alpha_prev: None,
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
    .and_then(|stepper| stepper.with_rayleigh_damping(eta_k, eta_m))
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // ── Resolve include_geometric_stiffness (explicit arg overrides sentinel) ──
    let kg_sentinel_active = kg0_rows.is_some() || kg0_cols.is_some() || kg0_vals.is_some();
    let include_geometric_stiffness_resolved = match include_geometric_stiffness {
        Some(b) => b,
        None => {
            // Backward-compat shim: if any kg0_* sentinel was passed, treat as true.
            if kg_sentinel_active {
                log::warn!(
                    "run_inertial_rotor_fsi_solver: kg0_rows/cols/vals sentinel is deprecated. \
                     Pass `include_geometric_stiffness=True` explicitly instead."
                );
                true
            } else {
                false
            }
        }
    };

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

    // ── Inertial Rotor FSI config ─────────────────────────────────────────────
    let config = InertialRotorFsiConfig {
        fsi: fsi_config,
        rotation_axis: axis,
        rotation_center: center,
        gravity: grav,
        include_reference_acceleration,
        include_geometric_stiffness: include_geometric_stiffness_resolved,
        include_spin_softening,
        ksp_omega_rebuild_high,
        ksp_omega_rebuild_low,
        kg_omega_rebuild_high,
        kg_omega_rebuild_low,
        k_update_interval,
        omega_rebuild_threshold,
        use_corotational_kt,
        kt_coro_update_freq,
        theta_rebuild_threshold,
        kg_use_deformed_coords,
        kg_deflection_rebuild_rel_high,
        kg_deflection_rebuild_rel_low,
        displacement_mode: disp_mode,
        dofs_per_node,
        fluid_density,
        flow_velocity,
        rotor_radius,
        omega_mesh_name,
        omega_write_data,
        omega_vertex_coord: omega_vertex_coord_arr,
        velocity_write_data,
    };

    // ── Initial state ─────────────────────────────────────────────────────────
    let initial_state = if let (Some(u_arr), Some(v_arr), Some(a_arr)) = (u0, v0, a0) {
        let u = u_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let v = v_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        let a = a_arr.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?.to_vec();
        FsiInitialState { u, v, a, t: t0 }
    } else {
        FsiInitialState {
            u: vec![0.0; n_free],
            v: vec![0.0; n_free],
            a: vec![0.0; n_free],
            t: t0,
        }
    };

    // ── Build and configure solver ────────────────────────────────────────────
    let rust_assembler = assembler.inner().clone();
    let (k_full_rows, k_full_cols, _) = rust_assembler.assemble_k();
    let k_coo_map = setup::build_kg_coo_map(&k_full_rows, &k_full_cols, fd, &kr_red, &kc_red);
    let k_red_nnz = kr_red.len();
    let mut solver = InertialRotorFsiSolver::new(
        config,
        rust_assembler,
        stepper,
        omega_provider,
        iface_nodes,
        masses,
        fd.to_vec(),
        k_coo_map,
        k_red_nnz,
        kr_red.to_vec(),
        kc_red.to_vec(),
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?
    .with_initial_state(&initial_state)
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // ── Attach callback if provided ───────────────────────────────────────────
    if let Some(py_cb) = step_callback {
        solver = solver.with_step_callback(Box::new(
            move |t, step, dt, u, v, a, force_mag, forces_iface,
                  omega, alpha, theta, tau_aero, perf| {
                Python::attach(|py| {
                    let u_arr = Array1::from(u.to_vec()).into_pyarray(py);
                    let v_arr = Array1::from(v.to_vec()).into_pyarray(py);
                    let a_arr = Array1::from(a.to_vec()).into_pyarray(py);
                    let fi_arr = Array1::from(forces_iface.to_vec()).into_pyarray(py);
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
            },
        ));
    }

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
