//! PyO3 bindings for the aeroelast finite-element stack.
//!
//! Item bodies live in the submodules below; the Python-facing API is declared
//! once in `register_module`, and the `#[pymodule]` entry point stays at the
//! crate root so the cdylib init symbol is unchanged.

mod assembler;
mod elements;
#[cfg(feature = "fsi")]
mod fsi;
mod materials;
mod mesh;
mod solvers;

use assembler::*;
use elements::*;
#[cfg(feature = "fsi")]
use fsi::*;
use materials::*;
use mesh::*;
use solvers::*;

use pyo3::prelude::*;

/// Register all aeroelast functions into a PyModule.
pub fn register_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Initialize Rust logger once (env_logger ignores subsequent calls)
    // Logs appear on stderr; Python CLI can redirect or capture as needed.
    // Default level: INFO. Override with RUST_LOG env var (e.g., RUST_LOG=debug).
    let _ = env_logger::try_init();

    m.add_function(wrap_pyfunction!(batch_ke_mitc3, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_mitc3, m)?)?;
    m.add_function(wrap_pyfunction!(batch_kt_mitc3, m)?)?;
    m.add_function(wrap_pyfunction!(batch_fint_mitc3, m)?)?;
    m.add_function(wrap_pyfunction!(assemble_smoothed_mitc3, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_mitc4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_mitc4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_kt_mitc4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_fint_mitc4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_mitc3_composite, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_mitc3_composite, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_mitc4_composite, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_mitc4_composite, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_quad4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_quad4, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_quad8, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_quad8, m)?)?;
    m.add_function(wrap_pyfunction!(batch_ke_quad9, m)?)?;
    m.add_function(wrap_pyfunction!(batch_me_quad9, m)?)?;
    m.add_function(wrap_pyfunction!(coo_assembly, m)?)?;
    m.add_function(wrap_pyfunction!(compute_nnz, m)?)?;
    m.add_function(wrap_pyfunction!(petsc_assemble_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(petsc_modal_solve, m)?)?;
    m.add_function(wrap_pyfunction!(modal_solve_coo, m)?)?;
    m.add_function(wrap_pyfunction!(linear_static_solve_coo, m)?)?;
    m.add_function(wrap_pyfunction!(newmark_beta_solve_coo, m)?)?;
    m.add_function(wrap_pyfunction!(nonlinear_static_solve_coo, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(run_linear_elastic_fsi, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(run_fsi_solver, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(run_stress_stiffened_fsi_solver, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(run_rotor_fsi_solver, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(run_inertial_rotor_fsi_solver, m)?)?;
    #[cfg(feature = "fsi")]
    m.add_function(wrap_pyfunction!(compute_rayleigh_auto, m)?)?;
    m.add_class::<PyMeshModel>()?;

    m.add_class::<PyMeshAssembler>()?;
    m.add_class::<PyOrthotropicMaterial>()?;
    m.add_class::<PyPly>()?;
    m.add_class::<PyLaminate>()?;
    m.add_class::<PyElementFamily>()?;
    Ok(())
}

// ============================================================================
// Module entry point (named `aeroelast`)
// ============================================================================

#[pymodule]
fn _aeroelast(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_module(m)
}
