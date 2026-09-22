//! Sparse (COO) and PETSc-backed solver bindings exposed to Python.

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;

use aeroelast_core::assembly;

use crate::assembler::PyMeshAssembler;

/// Generate COO triplets for sparse assembly from element matrices.
#[pyfunction]
pub(crate) fn coo_assembly<'py>(
    py: Python<'py>,
    dofs: PyReadonlyArray2<'py, i64>,
    ke_flat: PyReadonlyArray1<'py, f64>,
    n_dof_per_elem: usize,
) -> (
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
) {
    let dofs_arr = dofs.as_array();
    let ke_arr = ke_flat.as_array();
    let n_elem = dofs_arr.nrows();

    let dofs_flat: Vec<i64> = dofs_arr.iter().copied().collect();
    let ke_flat_slice: Vec<f64> = ke_arr.iter().copied().collect();

    let (rows, cols, vals) =
        assembly::coo_from_batch(&dofs_flat, &ke_flat_slice, n_elem, n_dof_per_elem);

    (
        Array1::from(rows).into_pyarray(py),
        Array1::from(cols).into_pyarray(py),
        Array1::from(vals).into_pyarray(py),
    )
}

/// Compute NNZ per row for PETSc preallocation.
#[pyfunction]
pub(crate) fn compute_nnz<'py>(
    py: Python<'py>,
    dofs: PyReadonlyArray2<'py, i64>,
    n_dof_total: usize,
) -> Bound<'py, PyArray1<i64>> {
    let dofs_arr = dofs.as_array();
    let n_elem = dofs_arr.nrows();
    let dof_per_elem = dofs_arr.ncols();

    let dofs_flat: Vec<i64> = dofs_arr.iter().copied().collect();
    let nnz = assembly::compute_nnz_per_row(&dofs_flat, n_elem, dof_per_elem, n_dof_total);

    Array1::from(nnz).into_pyarray(py)
}

/// Linear static solve from full-system COO triplets + free DOF list.
///
/// Assembles the reduced stiffness matrix and load vector, solves K_red·u_red = F_red
/// with PETSc KSP (CG + ICC), then expands the solution back to the full DOF space.
///
/// Args:
///   k_rows, k_cols, k_vals  – COO triplets for the full stiffness matrix (i64 indices)
///   f_full                  – full load vector (length n_dof_total)
///   n_dof_total             – total DOF count (full system size)
///   free_dofs               – 1-D array of free DOF indices (i64, sorted)
///
/// Returns:
///   u_full: np.ndarray shape (n_dof_total,) — displacement at all DOFs (0 at constrained)
#[pyfunction]
pub(crate) fn linear_static_solve_coo<'py>(
    py: Python<'py>,
    k_rows: PyReadonlyArray1<'py, i64>,
    k_cols: PyReadonlyArray1<'py, i64>,
    k_vals: PyReadonlyArray1<'py, f64>,
    f_full: PyReadonlyArray1<'py, f64>,
    n_dof_total: usize,
    free_dofs: PyReadonlyArray1<'py, i64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let kr = k_rows.as_slice()?;
    let kc = k_cols.as_slice()?;
    let kv = k_vals.as_slice()?;
    let f_arr = f_full.as_slice()?;
    let free = free_dofs.as_slice()?;

    // Build free DOF set and remapping: global → reduced index
    let mut free_sorted: Vec<i64> = free.to_vec();
    free_sorted.sort_unstable();
    let n_free = free_sorted.len();
    let free_map: std::collections::HashMap<i64, i32> = free_sorted
        .iter()
        .enumerate()
        .map(|(new_idx, &old_dof)| (old_dof, new_idx as i32))
        .collect();

    // Restrict K COO to free×free submatrix and remap indices to [0, n_free)
    let mut rr = Vec::new();
    let mut cc = Vec::new();
    let mut vv = Vec::new();
    for ((r, c), v) in kr.iter().zip(kc.iter()).zip(kv.iter()) {
        if let (Some(&ri), Some(&ci)) = (free_map.get(r), free_map.get(c)) {
            rr.push(ri);
            cc.push(ci);
            vv.push(*v);
        }
    }

    // Restrict F to free DOFs
    let f_red: Vec<f64> = free_sorted.iter().map(|&d| f_arr[d as usize]).collect();

    // Assemble PETSc K_red
    let mat_k = aeroelast_solvers::petsc::assembler::assemble_seq_aij(&rr, &cc, &vv, n_free)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Build PETSc F_red vector
    let f_vec = aeroelast_solvers::petsc::elasticity::static_linear::build_vec_from_slice(&f_red)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Solve K_red · u_red = F_red
    let result = aeroelast_solvers::petsc::elasticity::static_linear::linear_static_solve(&mat_k, &f_vec, n_free)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    if result.converged_reason <= 0 {
        return Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
            "KSP did not converge: reason={}, iters={}, rnorm={:.3e}",
            result.converged_reason, result.iterations, result.residual_norm
        )));
    }

    // Expand u_red → u_full (0 at constrained DOFs)
    let mut u_full = vec![0.0f64; n_dof_total];
    for (i, &global_dof) in free_sorted.iter().enumerate() {
        u_full[global_dof as usize] = result.displacements[i];
    }

    Ok(Array1::from(u_full).into_pyarray(py))
}

/// Newmark-β dynamic solve from full-system COO triplets + free DOF list.
///
/// Integrates M·ü + C·u̇ + K·u = F(t) using the Newmark-β scheme.
/// Rayleigh damping: C = η_k·K + η_m·M.
///
/// Args:
///   k_rows, k_cols, k_vals  – COO triplets for the full stiffness matrix (i64 indices)
///   m_rows, m_cols, m_vals  – COO triplets for the full mass matrix (i64 indices)
///   f_history_flat          – force history, shape (n_steps+1, n_dof_total), row-major
///   n_dof_total             – total DOF count (full system)
///   free_dofs               – 1-D array of free DOF indices (i64, sorted)
///   dt                      – time step (seconds)
///   n_steps                 – number of time steps
///   eta_k                   – Rayleigh stiffness proportional damping (default 0.0)
///   eta_m                   – Rayleigh mass proportional damping (default 0.0)
///   beta                    – Newmark-β parameter (default 0.25)
///   gamma                   – Newmark-γ parameter (default 0.5)
///
/// Returns:
///   (u_hist, v_hist, a_hist) each shape (n_steps+1, n_dof_total), row-major flat
#[pyfunction]
#[pyo3(signature = (k_rows, k_cols, k_vals, m_rows, m_cols, m_vals, f_history_flat, n_dof_total, free_dofs, dt, n_steps, eta_k=0.0, eta_m=0.0, beta=0.25, gamma=0.5))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn newmark_beta_solve_coo<'py>(
    py: Python<'py>,
    k_rows: PyReadonlyArray1<'py, i64>,
    k_cols: PyReadonlyArray1<'py, i64>,
    k_vals: PyReadonlyArray1<'py, f64>,
    m_rows: PyReadonlyArray1<'py, i64>,
    m_cols: PyReadonlyArray1<'py, i64>,
    m_vals: PyReadonlyArray1<'py, f64>,
    f_history_flat: PyReadonlyArray2<'py, f64>,
    n_dof_total: usize,
    free_dofs: PyReadonlyArray1<'py, i64>,
    dt: f64,
    n_steps: usize,
    eta_k: f64,
    eta_m: f64,
    beta: f64,
    gamma: f64,
) -> PyResult<(
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let kr = k_rows.as_slice()?;
    let kc = k_cols.as_slice()?;
    let kv = k_vals.as_slice()?;
    let mr = m_rows.as_slice()?;
    let mc = m_cols.as_slice()?;
    let mv = m_vals.as_slice()?;
    let free = free_dofs.as_slice()?;
    let f_arr = f_history_flat.as_array();

    // ── Build free DOF remapping ─────────────────────────────────────────────
    let mut free_sorted: Vec<i64> = free.to_vec();
    free_sorted.sort_unstable();
    let n_free = free_sorted.len();
    let free_map: std::collections::HashMap<i64, i32> = free_sorted
        .iter()
        .enumerate()
        .map(|(new_idx, &old_dof)| (old_dof, new_idx as i32))
        .collect();

    // ── Restrict K and M COO to free×free ────────────────────────────────────
    let restrict = |rows: &[i64], cols: &[i64], vals: &[f64]| -> (Vec<i32>, Vec<i32>, Vec<f64>) {
        let mut rr = Vec::new();
        let mut cc = Vec::new();
        let mut vv = Vec::new();
        for ((r, c), v) in rows.iter().zip(cols.iter()).zip(vals.iter()) {
            if let (Some(&ri), Some(&ci)) = (free_map.get(r), free_map.get(c)) {
                rr.push(ri);
                cc.push(ci);
                vv.push(*v);
            }
        }
        (rr, cc, vv)
    };

    let (rk, ck, vk) = restrict(kr, kc, kv);
    let (rm, cm, vm) = restrict(mr, mc, mv);

    // ── Restrict F history ───────────────────────────────────────────────────
    let n_time = n_steps + 1;
    if f_arr.nrows() != n_time {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "f_history_flat must have {} rows (n_steps+1), got {}",
            n_time, f_arr.nrows()
        )));
    }
    if f_arr.ncols() != n_dof_total {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "f_history_flat must have {} cols (n_dof_total), got {}",
            n_dof_total, f_arr.ncols()
        )));
    }

    let f_history_red: Vec<Vec<f64>> = (0..n_time)
        .map(|t| free_sorted.iter().map(|&d| f_arr[[t, d as usize]]).collect())
        .collect();

    // ── Solve ────────────────────────────────────────────────────────────────
    let result = aeroelast_solvers::petsc::elasticity::dynamic_newmark::newmark_beta_solve(
        &rk, &ck, &vk,
        &rm, &cm, &vm,
        eta_k, eta_m,
        &f_history_red,
        dt,
        n_steps,
        n_free,
        beta,
        gamma,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // ── Expand histories to full DOF space ───────────────────────────────────
    let expand = |hist: Vec<Vec<f64>>| -> Vec<f64> {
        let mut out = vec![0.0f64; n_time * n_dof_total];
        for (t, step) in hist.iter().enumerate() {
            for (i, &global_dof) in free_sorted.iter().enumerate() {
                out[t * n_dof_total + global_dof as usize] = step[i];
            }
        }
        out
    };

    let u_flat = expand(result.displacements);
    let v_flat = expand(result.velocities);
    let a_flat = expand(result.accelerations);

    Ok((
        Array1::from(u_flat).into_pyarray(py),
        Array1::from(v_flat).into_pyarray(py),
        Array1::from(a_flat).into_pyarray(py),
    ))
}

// ============================================================================
// PETSc pipeline: assemble + modal solve via SLEPc
// ============================================================================

/// Stable name for capsules wrapping a `PetscMat`.
///
/// `petsc_modal_solve` verifies this name via `PyCapsule::pointer_checked`
/// before casting the capsule pointer to `PetscMat`, preventing UB from
/// arbitrary capsules being passed in.
pub(crate) const PETSC_MAT_CAPSULE_NAME: &std::ffi::CStr =
    unsafe { std::ffi::CStr::from_bytes_with_nul_unchecked(b"aeroelast.PetscMat\0") };

/// Assemble a global stiffness or mass matrix into a PETSc Mat (AIJ sequential).
///
/// Accepts COO triplets produced by `coo_assembly` and returns an opaque
/// Python capsule wrapping the assembled `PetscMat`. The capsule is consumed
/// by `petsc_modal_solve`.
///
/// Args:
///   rows:        COO row indices, dtype=int32
///   cols:        COO col indices, dtype=int32
///   vals:        COO values, dtype=float64
///   n_dof:       total number of DOFs (matrix is n_dof × n_dof)
///
/// Returns:
///   A Python capsule wrapping the PETSc Mat handle.
///   The handle is freed automatically when the capsule is garbage-collected.
#[pyfunction]
pub(crate) fn petsc_assemble_matrix<'py>(
    py: Python<'py>,
    rows: PyReadonlyArray1<'py, i32>,
    cols: PyReadonlyArray1<'py, i32>,
    vals: PyReadonlyArray1<'py, f64>,
    n_dof: usize,
) -> PyResult<Bound<'py, PyCapsule>> {
    let rows_s = rows.as_slice()?;
    let cols_s = cols.as_slice()?;
    let vals_s = vals.as_slice()?;

    let mat = aeroelast_solvers::petsc::assembler::assemble_seq_aij(rows_s, cols_s, vals_s, n_dof)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Move `mat` into the capsule. PyO3 boxes it internally; Drop runs MatDestroy.
    // Name the capsule so `petsc_modal_solve` can validate it before casting.
    PyCapsule::new(py, mat, Some(std::ffi::CString::from(PETSC_MAT_CAPSULE_NAME)))
}

/// Modal solve using SLEPc EPS on two PETSc matrices (K and M).
///
/// Args:
///   k_capsule:  Python capsule wrapping the stiffness PetscMat
///   m_capsule:  Python capsule wrapping the mass PetscMat
///   n_modes:    number of modes to compute
///
/// Returns:
///   (eigenvalues, eigenvectors_flat)
///   eigenvalues:      np.ndarray shape (n_conv,)  — ω² (rad²/s²)
///   eigenvectors_flat: np.ndarray shape (n_conv * n_dof,) — row-major mode shapes
#[pyfunction]
pub(crate) fn petsc_modal_solve<'py>(
    py: Python<'py>,
    k_capsule: &Bound<'py, PyCapsule>,
    m_capsule: &Bound<'py, PyCapsule>,
    n_modes: usize,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    // SAFETY: the capsule name must match the one set by `petsc_assemble_matrix`
    // (`aeroelast.PetscMat`). `pointer_checked(Some(name))` validates the capsule
    // name AND non-null pointer before we cast and reborrow as a shared ref, so a
    // mismatched/wrong capsule returns a Python error instead of UB.
    let k = unsafe {
        k_capsule
            .pointer_checked(Some(PETSC_MAT_CAPSULE_NAME))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?
            .cast::<aeroelast_solvers::PetscMat>()
            .as_ref()
    };
    let m = unsafe {
        m_capsule
            .pointer_checked(Some(PETSC_MAT_CAPSULE_NAME))
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?
            .cast::<aeroelast_solvers::PetscMat>()
            .as_ref()
    };

    let result = aeroelast_solvers::petsc::elasticity::modal::modal_solve(k, m, n_modes)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let eigenvalues = Array1::from(result.eigenvalues).into_pyarray(py);
    let eigenvectors_flat: Vec<f64> = result.eigenvectors.into_iter().flatten().collect();
    let eigenvectors = Array1::from(eigenvectors_flat).into_pyarray(py);

    Ok((eigenvalues, eigenvectors))
}

/// All-in-one modal solve from full-system COO triplets + free DOF list.
///
/// Restricts the COO system to the free DOFs submatrix, assembles the
/// reduced PETSc Mats, solves the generalized eigenvalue problem K·x = λ·M·x,
/// and returns natural frequencies in Hz (sorted ascending).
///
/// Args:
///   k_rows, k_cols, k_vals  – COO triplets for the full stiffness matrix (i64 indices)
///   m_rows, m_cols, m_vals  – COO triplets for the full mass matrix (i64 indices)
///   n_dof_total             – total DOF count (full system size)
///   free_dofs               – 1-D array of free DOF indices (i64, sorted)
///   n_modes                 – number of modes to compute
///
/// Returns:
///   (frequencies_hz, modes_flat)
///   frequencies_hz:  np.ndarray shape (n_conv,)  — natural frequencies in Hz
///   modes_flat:      np.ndarray shape (n_conv * n_free,) — row-major, reduced space
#[pyfunction]
pub(crate) fn modal_solve_coo<'py>(
    py: Python<'py>,
    k_rows: PyReadonlyArray1<'py, i64>,
    k_cols: PyReadonlyArray1<'py, i64>,
    k_vals: PyReadonlyArray1<'py, f64>,
    m_rows: PyReadonlyArray1<'py, i64>,
    m_cols: PyReadonlyArray1<'py, i64>,
    m_vals: PyReadonlyArray1<'py, f64>,
    n_dof_total: usize,
    free_dofs: PyReadonlyArray1<'py, i64>,
    n_modes: usize,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let kr = k_rows.as_slice()?;
    let kc = k_cols.as_slice()?;
    let kv = k_vals.as_slice()?;
    let mr = m_rows.as_slice()?;
    let mc = m_cols.as_slice()?;
    let mv = m_vals.as_slice()?;
    let free = free_dofs.as_slice()?;

    // Build free DOF set and remapping: old global → new reduced index
    let mut free_sorted: Vec<i64> = free.to_vec();
    free_sorted.sort_unstable();
    let n_free = free_sorted.len();
    // HashMap: global dof → reduced index
    let free_map: std::collections::HashMap<i64, i32> = free_sorted
        .iter()
        .enumerate()
        .map(|(new_idx, &old_dof)| (old_dof, new_idx as i32))
        .collect();

    // Restrict COO to free×free submatrix and remap indices to [0, n_free)
    let restrict = |rows: &[i64], cols: &[i64], vals: &[f64]| -> (Vec<i32>, Vec<i32>, Vec<f64>) {
        let mut rr = Vec::new();
        let mut cc = Vec::new();
        let mut vv = Vec::new();
        for ((r, c), v) in rows.iter().zip(cols.iter()).zip(vals.iter()) {
            if let (Some(&ri), Some(&ci)) = (free_map.get(r), free_map.get(c)) {
                rr.push(ri);
                cc.push(ci);
                vv.push(*v);
            }
        }
        (rr, cc, vv)
    };

    let (rk, ck, vk) = restrict(kr, kc, kv);
    let (rm, cm, vm) = restrict(mr, mc, mv);

    // Assemble reduced PETSc matrices
    let mat_k = aeroelast_solvers::petsc::assembler::assemble_seq_aij(&rk, &ck, &vk, n_free)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
    let mat_m = aeroelast_solvers::petsc::assembler::assemble_seq_aij(&rm, &cm, &vm, n_free)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Solve eigenvalue problem
    let result = aeroelast_solvers::petsc::elasticity::modal::modal_solve(&mat_k, &mat_m, n_modes)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    // Convert eigenvalues ω² → frequencies in Hz, filter positive, sort
    let mut pairs: Vec<(f64, Vec<f64>)> = result
        .eigenvalues
        .into_iter()
        .zip(result.eigenvectors.into_iter())
        .filter(|(lam, _)| *lam > 1e-8)
        .map(|(lam, vec)| (lam.sqrt() / (2.0 * std::f64::consts::PI), vec))
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    pairs.truncate(n_modes);

    let frequencies: Vec<f64> = pairs.iter().map(|(f, _)| *f).collect();
    let modes_flat: Vec<f64> = pairs.into_iter().flat_map(|(_, v)| v).collect();

    let _ = n_dof_total; // available if caller needs expansion; not used here

    Ok((
        Array1::from(frequencies).into_pyarray(py),
        Array1::from(modes_flat).into_pyarray(py),
    ))
}

// ============================================================================
// Nonlinear static solver via PETSc SNES
// ============================================================================

/// Solve the nonlinear static FEM system R(u) = F_int(u) - F_ext = 0
/// using PETSc SNES (Newton-Raphson with line search by default).
///
/// The solver operates on the full DOF space — Dirichlet BCs are enforced
/// internally via identity rows in the Jacobian and zeroing the residual.
///
/// To switch to arc-length continuation pass ``petsc_options="-snes_type newtonal"``.
///
/// Args:
///   assembler:       PyMeshAssembler wrapping the assembled mesh
///   f_ext:           External load vector (length = n_dof), dtype=float64
///   dirichlet_dofs:  Indices of constrained DOFs, dtype=int64
///   atol:            Absolute residual tolerance (default 1e-10)
///   rtol:            Relative residual tolerance (default 1e-8)
///   stol:            Step-length tolerance (default 1e-8)
///   max_it:          Maximum Newton iterations (default 50)
///   x0:              Optional initial guess displacement vector (length = n_dof)
///
/// Returns:
///   Tuple ``(displacements, iterations, residual_norm, converged_reason)``
#[pyfunction]
#[pyo3(signature = (assembler, f_ext, dirichlet_dofs, atol=1e-10, rtol=1e-8, stol=1e-8, max_it=50, x0=None, diagnostics=false, diagnostics_every=1))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn nonlinear_static_solve_coo<'py>(
    py: Python<'py>,
    assembler: &PyMeshAssembler,
    f_ext: PyReadonlyArray1<'py, f64>,
    dirichlet_dofs: PyReadonlyArray1<'py, i64>,
    atol: f64,
    rtol: f64,
    stol: f64,
    max_it: i32,
    x0: Option<PyReadonlyArray1<'py, f64>>,
    diagnostics: bool,
    diagnostics_every: i32,
) -> PyResult<(Bound<'py, PyArray1<f64>>, i32, f64, i32)> {
    let f_ext_slice = f_ext.as_slice()?;
    let dofs_i64 = dirichlet_dofs.as_slice()?;
    let dofs_usize: Vec<usize> = dofs_i64.iter().map(|&d| d as usize).collect();
    // Own optional x0 to keep memory alive across the Rust solver call.
    let x0_owned: Option<Vec<f64>> = match x0 {
        Some(arr) => Some(arr.as_slice()?.to_vec()),
        None => None,
    };
    let x0_slice: Option<&[f64]> = x0_owned.as_deref();

    let config = aeroelast_solvers::petsc::elasticity::static_nonlinear::NonlinearConfig {
        atol,
        rtol,
        stol,
        max_it,
        max_funcs: -1,
        diagnostics,
        diagnostics_every,
    };

    let result = aeroelast_solvers::petsc::elasticity::static_nonlinear::nonlinear_static_solve_with_guess(
        assembler.inner(),
        f_ext_slice,
        &dofs_usize,
        x0_slice,
        &config,
    )
    .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;

    let u_arr = Array1::from(result.displacements).into_pyarray(py);
    Ok((u_arr, result.iterations, result.residual_norm, result.converged_reason))
}

// ── compute_rayleigh_auto ─────────────────────────────────────────────────────
/// Compute Rayleigh damping coefficients (η_k, η_m) from the global stiffness
/// and mass COO matrices using the two-point formula.
///
/// Reduces K and M by the `free_dofs` mask, assembles PETSc matrices, runs a
/// SLEPc modal analysis via the Rust `modal_solve` function, then applies:
///
///   η_k = 2 (ζ_i·ω_i − ζ_j·ω_j) / (ω_i² − ω_j²)  
///   η_m = 2·ω_i·ω_j·(ζ_j·ω_i − ζ_i·ω_j) / (ω_i² − ω_j²)
///
/// This replaces the Python SLEPc boilerplate in `_compute_rayleigh_auto` and
/// re-uses the same `modal.rs` solver that is used by the Rust modal analysis.
///
/// # Arguments
/// * `k_rows`, `k_cols`, `k_vals` — global stiffness COO (full system)
/// * `m_rows`, `m_cols`, `m_vals` — global consistent mass COO (full system)
/// * `free_dofs`                   — sorted free (unconstrained) DOF indices
/// * `num_modes`                   — number of modes to request from SLEPc
/// * `mode_i`, `mode_j`            — 1-based mode indices for the two-point formula
/// * `zeta_i`, `zeta_j`            — target damping ratios for modes i and j
///
/// # Returns
/// `(eta_k, eta_m)` — stiffness-proportional and mass-proportional Rayleigh coefficients.
#[pyfunction]
#[cfg(feature = "fsi")]
pub(crate) fn compute_rayleigh_auto(
    k_rows: PyReadonlyArray1<i32>,
    k_cols: PyReadonlyArray1<i32>,
    k_vals: PyReadonlyArray1<f64>,
    m_rows: PyReadonlyArray1<i32>,
    m_cols: PyReadonlyArray1<i32>,
    m_vals: PyReadonlyArray1<f64>,
    free_dofs: PyReadonlyArray1<i32>,
    num_modes: usize,
    mode_i: usize,
    mode_j: usize,
    zeta_i: f64,
    zeta_j: f64,
) -> PyResult<(f64, f64)> {
    use aeroelast_solvers::petsc::fsi::setup;
    use aeroelast_solvers::petsc::assembler::assemble_seq_aij;
    use aeroelast_solvers::petsc::elasticity::modal::modal_solve;
    use pyo3::exceptions::PyRuntimeError;
    use pyo3::exceptions::PyValueError;

    let kr = k_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kc = k_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let kv = k_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mr = m_rows.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mc = m_cols.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let mv = m_vals.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let fd = free_dofs.as_slice().map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    if mode_i == 0 || mode_j == 0 || mode_i == mode_j {
        return Err(PyValueError::new_err(
            "mode_i and mode_j must be distinct 1-based mode indices",
        ));
    }

    let n_free = fd.len();

    // ── Reduce K (and M sharing K's sparsity pattern) ─────────────────────
    let (kr_red, kc_red, kv_red) = setup::reduce_coo(kr, kc, kv, fd);
    let (_, _, mv_red)            = setup::reduce_coo(mr, mc, mv, fd);

    // Lump mass diagonal from the reduced consistent mass
    let mut m_diag = vec![0.0f64; n_free];
    for (&row, &val) in kr_red.iter().zip(mv_red.iter()) {
        if row >= 0 && (row as usize) < n_free {
            m_diag[row as usize] += val.abs();
        }
    }
    // Distribute diagonal values onto K's sparsity for PETSc assembly
    let mv_expanded = setup::expand_diag_to_sparsity(&m_diag, &kr_red, &kc_red);

    // ── Assemble PETSc matrices ────────────────────────────────────────────
    let k_mat = assemble_seq_aij(&kr_red, &kc_red, &kv_red, n_free)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    let m_mat = assemble_seq_aij(&kr_red, &kc_red, &mv_expanded, n_free)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // ── Modal analysis ─────────────────────────────────────────────────────
    let result = modal_solve(&k_mat, &m_mat, num_modes)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    // Keep only positive eigenvalues (λ = ω²) — rigid body modes are near zero.
    let pos_eigs: Vec<f64> = result
        .eigenvalues
        .iter()
        .copied()
        .filter(|&v| v > 1e-8)
        .collect();

    let max_mode = mode_i.max(mode_j);
    if pos_eigs.len() < max_mode {
        return Err(PyValueError::new_err(format!(
            "Only {} positive modes converged, but mode {} was requested. \
             Increase num_modes or check the model.",
            pos_eigs.len(),
            max_mode
        )));
    }

    // ── Two-point Rayleigh formula ─────────────────────────────────────────
    let omega_i = pos_eigs[mode_i - 1].sqrt();
    let omega_j = pos_eigs[mode_j - 1].sqrt();
    let denom = omega_i * omega_i - omega_j * omega_j;

    if denom.abs() < 1e-16 {
        return Err(PyValueError::new_err(
            "modes i and j have identical frequencies; cannot compute Rayleigh coefficients",
        ));
    }

    let eta_k = 2.0 * (zeta_i * omega_i - zeta_j * omega_j) / denom;
    let eta_m = 2.0 * omega_i * omega_j * (zeta_j * omega_i - zeta_i * omega_j) / denom;

    Ok((eta_k, eta_m))
}
