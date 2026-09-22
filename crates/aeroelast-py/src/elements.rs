//! Batch element kernels exposed to Python (MITC3/MITC4 shells, solid and quad elements).

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::prelude::*;

use aeroelast_core::elements::mitc3::{self, Mitc3Precomputed};
use aeroelast_core::elements::mitc4::{self, Mitc4Precomputed};
use aeroelast_core::elements::quad::{Quad4Precomputed, Quad8Precomputed, Quad9Precomputed};
use aeroelast_core::elements::solid;
use aeroelast_core::materials::composite::composite_constitutive;
use aeroelast_core::materials::isotropic::IsotropicMaterial;
use aeroelast_core::Material;

// ============================================================================
// PyO3 exposed functions
// ============================================================================

/// Batch-compute MITC3+ element stiffness matrices (global coords).
#[pyfunction]
#[pyo3(signature = (coords, e_mod, nu, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_ke_mitc3<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc3Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);
            let ke = mitc3::compute_ke_global(&pre);
            let mut flat = [0.0f64; 324];
            for i in 0..18 {
                for j in 0..18 {
                    flat[i * 18 + j] = ke[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 324;
    let mut output = Vec::with_capacity(total);
    for ke_flat in &results {
        output.extend_from_slice(ke_flat);
    }

    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC3+ element mass matrices (global coords).
#[pyfunction]
#[pyo3(signature = (coords, e_mod, nu, rho, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_me_mitc3<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    rho: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, rho);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc3Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);
            let me = mitc3::compute_me_global(&pre, rho);
            let mut flat = [0.0f64; 324];
            for i in 0..18 {
                for j in 0..18 {
                    flat[i * 18 + j] = me[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 324;
    let mut output = Vec::with_capacity(total);
    for me_flat in &results {
        output.extend_from_slice(me_flat);
    }

    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC3+ tangent stiffness matrices (global coords, nonlinear).
#[pyfunction]
#[pyo3(signature = (coords, displacements, e_mod, nu, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_kt_mitc3<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    displacements: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let disp_arr = displacements.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc3Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);

            let mut u = mitc3::Vec18::zeros();
            for i in 0..18 {
                u[i] = disp_arr[[e, i]];
            }

            let kt = mitc3::compute_kt_global(&pre, &u);
            let mut flat = [0.0f64; 324];
            for i in 0..18 {
                for j in 0..18 {
                    flat[i * 18 + j] = kt[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 324;
    let mut output = Vec::with_capacity(total);
    for kt_flat in &results {
        output.extend_from_slice(kt_flat);
    }

    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC3+ internal force vectors (global coords).
#[pyfunction]
#[pyo3(signature = (coords, displacements, e_mod, nu, thickness, shear_correction=5.0/6.0, nonlinear=true))]
pub(crate) fn batch_fint_mitc3<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    displacements: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
    nonlinear: bool,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let disp_arr = displacements.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 18]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc3Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);

            let mut u = mitc3::Vec18::zeros();
            for i in 0..18 {
                u[i] = disp_arr[[e, i]];
            }

            let f = mitc3::compute_fint_global(&pre, &u, nonlinear);
            let mut flat = [0.0f64; 18];
            for i in 0..18 {
                flat[i] = f[i];
            }
            flat
        })
        .collect();

    let total = n_elem * 18;
    let mut output = Vec::with_capacity(total);
    for f_flat in &results {
        output.extend_from_slice(f_flat);
    }

    Array1::from(output).into_pyarray(py)
}

// ============================================================================
// MITC4 batch functions
// ============================================================================

/// Batch-compute MITC4+ element stiffness matrices (global coords).
#[pyfunction]
#[pyo3(signature = (coords, e_mod, nu, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_ke_mitc4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);
            let ke = mitc4::compute_ke_global(&pre);
            let mut flat = [0.0f64; 576];
            for i in 0..24 {
                for j in 0..24 {
                    flat[i * 24 + j] = ke[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 576;
    let mut output = Vec::with_capacity(total);
    for ke_flat in &results {
        output.extend_from_slice(ke_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC4+ element mass matrices (global coords).
#[pyfunction]
#[pyo3(signature = (coords, e_mod, nu, rho, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_me_mitc4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    rho: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, rho);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);
            let me = mitc4::compute_me_global(&pre, rho);
            let mut flat = [0.0f64; 576];
            for i in 0..24 {
                for j in 0..24 {
                    flat[i * 24 + j] = me[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 576;
    let mut output = Vec::with_capacity(total);
    for me_flat in &results {
        output.extend_from_slice(me_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC4+ tangent stiffness matrices (global coords, nonlinear).
#[pyfunction]
#[pyo3(signature = (coords, displacements, e_mod, nu, thickness, shear_correction=5.0/6.0))]
pub(crate) fn batch_kt_mitc4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    displacements: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let disp_arr = displacements.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);

            let mut u = mitc4::Vec24::zeros();
            for i in 0..24 {
                u[i] = disp_arr[[e, i]];
            }

            let kt = mitc4::compute_kt_global(&pre, &u);
            let mut flat = [0.0f64; 576];
            for i in 0..24 {
                for j in 0..24 {
                    flat[i * 24 + j] = kt[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 576;
    let mut output = Vec::with_capacity(total);
    for kt_flat in &results {
        output.extend_from_slice(kt_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC4+ internal force vectors (global coords).
#[pyfunction]
#[pyo3(signature = (coords, displacements, e_mod, nu, thickness, shear_correction=5.0/6.0, nonlinear=true))]
pub(crate) fn batch_fint_mitc4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    displacements: PyReadonlyArray2<'py, f64>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
    nonlinear: bool,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let disp_arr = displacements.as_array();
    let n_elem = coords_arr.nrows();

    let mat = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = mat.constitutive(thickness, shear_correction);

    let results: Vec<[f64; 24]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, e_mod, 1.0);

            let mut u = mitc4::Vec24::zeros();
            for i in 0..24 {
                u[i] = disp_arr[[e, i]];
            }

            let f = mitc4::compute_fint_global(&pre, &u, nonlinear);
            let mut flat = [0.0f64; 24];
            for i in 0..24 {
                flat[i] = f[i];
            }
            flat
        })
        .collect();

    let total = n_elem * 24;
    let mut output = Vec::with_capacity(total);
    for f_flat in &results {
        output.extend_from_slice(f_flat);
    }
    Array1::from(output).into_pyarray(py)
}

// ============================================================================
// Composite batch functions
// ============================================================================

/// Batch-compute MITC3+ composite element stiffness matrices (global coords).
#[pyfunction]
pub(crate) fn batch_ke_mitc3_composite<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    cm_flat: PyReadonlyArray2<'py, f64>,
    b_coupling_flat: PyReadonlyArray2<'py, f64>,  // ADDED
    cb_flat: PyReadonlyArray2<'py, f64>,
    cs_flat: PyReadonlyArray2<'py, f64>,
    thickness: PyReadonlyArray1<'py, f64>,
    e_equiv: PyReadonlyArray1<'py, f64>,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let cm_arr = cm_flat.as_array();
    let b_arr = b_coupling_flat.as_array();  // ADDED
    let cb_arr = cb_flat.as_array();
    let cs_arr = cs_flat.as_array();
    let h_arr = thickness.as_array();
    let e_arr = e_equiv.as_array();
    let n_elem = coords_arr.nrows();

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let mut a = [0.0f64; 9];
            let mut b = [0.0f64; 9];
            let mut d = [0.0f64; 9];
            let mut cs = [0.0f64; 4];
            for i in 0..9 {
                a[i] = cm_arr[[e, i]];
                b[i] = b_arr[[e, i]];
                d[i] = cb_arr[[e, i]];
            }
            for i in 0..4 {
                cs[i] = cs_arr[[e, i]];
            }
            let h = h_arr[e];
            let e_eq = e_arr[e];

            let constitutive = composite_constitutive(&a, &b, &d, &cs, h);
            let pre = Mitc3Precomputed::new(&node_coords, constitutive, h, e_eq, 1.0);
            let ke = mitc3::compute_ke_global(&pre);
            let mut flat = [0.0f64; 324];
            for i in 0..18 {
                for j in 0..18 {
                    flat[i * 18 + j] = ke[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 324;
    let mut output = Vec::with_capacity(total);
    for ke_flat in &results {
        output.extend_from_slice(ke_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC3+ composite element mass matrices (global coords).
#[pyfunction]
pub(crate) fn batch_me_mitc3_composite<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    mass_per_area: PyReadonlyArray1<'py, f64>,
    rotational_inertia: PyReadonlyArray1<'py, f64>,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let mpa_arr = mass_per_area.as_array();
    let ri_arr = rotational_inertia.as_array();
    let n_elem = coords_arr.nrows();

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 9];
            for i in 0..9 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let m_trans = mpa_arr[e];
            let m_rot = ri_arr[e];

            // Dummy constitutive for mass matrix (not used for stress computation)
            let dummy = aeroelast_core::materials::ShellConstitutive {
                cm: nalgebra::Matrix3::identity(),
                cb_coupling: nalgebra::Matrix3::zeros(),
                cb: nalgebra::Matrix3::identity(),
                cs: nalgebra::Matrix2::identity(),
                cm_raw: nalgebra::Matrix3::identity(),
            };
            let pre = Mitc3Precomputed::new(&node_coords, dummy, 1.0, 1.0, 1.0);
            let me = mitc3::compute_me_composite_global(&pre, m_trans, m_rot);
            let mut flat = [0.0f64; 324];
            for i in 0..18 {
                for j in 0..18 {
                    flat[i * 18 + j] = me[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 324;
    let mut output = Vec::with_capacity(total);
    for me_flat in &results {
        output.extend_from_slice(me_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC4+ composite element stiffness matrices (global coords).
#[pyfunction]
pub(crate) fn batch_ke_mitc4_composite<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    cm_flat: PyReadonlyArray2<'py, f64>,
    b_coupling_flat: PyReadonlyArray2<'py, f64>,  // ADDED
    cb_flat: PyReadonlyArray2<'py, f64>,
    cs_flat: PyReadonlyArray2<'py, f64>,
    thickness: PyReadonlyArray1<'py, f64>,
    e_equiv: PyReadonlyArray1<'py, f64>,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let cm_arr = cm_flat.as_array();
    let b_arr = b_coupling_flat.as_array();  // ADDED
    let cb_arr = cb_flat.as_array();
    let cs_arr = cs_flat.as_array();
    let h_arr = thickness.as_array();
    let e_arr = e_equiv.as_array();
    let n_elem = coords_arr.nrows();

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let mut a = [0.0f64; 9];
            let mut b = [0.0f64; 9];  // ADDED
            let mut d = [0.0f64; 9];
            let mut cs = [0.0f64; 4];
            for i in 0..9 {
                a[i] = cm_arr[[e, i]];
                b[i] = b_arr[[e, i]];  // ADDED
                d[i] = cb_arr[[e, i]];
            }
            for i in 0..4 {
                cs[i] = cs_arr[[e, i]];
            }
            let h = h_arr[e];
            let e_eq = e_arr[e];

            let constitutive = composite_constitutive(&a, &b, &d, &cs, h);
            let pre = Mitc4Precomputed::new(&node_coords, constitutive, h, e_eq, 1.0);
            let ke = mitc4::compute_ke_global(&pre);
            let mut flat = [0.0f64; 576];
            for i in 0..24 {
                for j in 0..24 {
                    flat[i * 24 + j] = ke[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 576;
    let mut output = Vec::with_capacity(total);
    for ke_flat in &results {
        output.extend_from_slice(ke_flat);
    }
    Array1::from(output).into_pyarray(py)
}

/// Batch-compute MITC4+ composite element mass matrices (global coords).
#[pyfunction]
pub(crate) fn batch_me_mitc4_composite<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    mass_per_area: PyReadonlyArray1<'py, f64>,
    rotational_inertia: PyReadonlyArray1<'py, f64>,
) -> Bound<'py, PyArray1<f64>> {
    let coords_arr = coords.as_array();
    let mpa_arr = mass_per_area.as_array();
    let ri_arr = rotational_inertia.as_array();
    let n_elem = coords_arr.nrows();

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|e| {
            let mut node_coords = [0.0f64; 12];
            for i in 0..12 {
                node_coords[i] = coords_arr[[e, i]];
            }
            let m_trans = mpa_arr[e];
            let m_rot = ri_arr[e];

            // Dummy constitutive for mass matrix (not used for stress computation)
            let dummy = aeroelast_core::materials::ShellConstitutive {
                cm: nalgebra::Matrix3::identity(),
                cb_coupling: nalgebra::Matrix3::zeros(),
                cb: nalgebra::Matrix3::identity(),
                cs: nalgebra::Matrix2::identity(),
                cm_raw: nalgebra::Matrix3::identity(),
            };
            let pre = Mitc4Precomputed::new(&node_coords, dummy, 1.0, 1.0, 1.0);
            let me = mitc4::compute_me_composite_global(&pre, m_trans, m_rot);
            let mut flat = [0.0f64; 576];
            for i in 0..24 {
                for j in 0..24 {
                    flat[i * 24 + j] = me[(i, j)];
                }
            }
            flat
        })
        .collect();

    let total = n_elem * 576;
    let mut output = Vec::with_capacity(total);
    for me_flat in &results {
        output.extend_from_slice(me_flat);
    }
    Array1::from(output).into_pyarray(py)
}

// ============================================================================
// Solid element batch functions: HEXA8, TETRA4, WEDGE6
// ============================================================================

/// Batch-compute HEXA8 element stiffness matrices.
///
/// coords: shape (n_elem, 24) — each row: 8 nodes × [x,y,z]
/// Returns flat array (n_elem × 576).
#[pyfunction]
pub(crate) fn batch_ke_hexa8<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 8];
            for n in 0..8 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::hexa8_ke(&c, e, nu);
            let mut flat = [0.0f64; 576];
            for r in 0..24 {
                for col in 0..24 {
                    flat[r * 24 + col] = ke[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 576);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute HEXA8 element mass matrices.
///
/// coords: shape (n_elem, 24)
/// Returns flat array (n_elem × 576).
#[pyfunction]
pub(crate) fn batch_me_hexa8<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 576]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 8];
            for n in 0..8 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::hexa8_me(&c, rho);
            let mut flat = [0.0f64; 576];
            for r in 0..24 {
                for col in 0..24 {
                    flat[r * 24 + col] = me[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 576);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute TETRA4 element stiffness matrices.
///
/// coords: shape (n_elem, 12) — each row: 4 nodes × [x,y,z]
/// Returns flat array (n_elem × 144).
#[pyfunction]
pub(crate) fn batch_ke_tetra4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 144]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 4];
            for n in 0..4 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::tetra4_ke(&c, e, nu);
            let mut flat = [0.0f64; 144];
            for r in 0..12 {
                for col in 0..12 {
                    flat[r * 12 + col] = ke[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 144);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute TETRA4 element mass matrices.
///
/// coords: shape (n_elem, 12)
/// Returns flat array (n_elem × 144).
#[pyfunction]
pub(crate) fn batch_me_tetra4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 144]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 4];
            for n in 0..4 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::tetra4_me(&c, rho);
            let mut flat = [0.0f64; 144];
            for r in 0..12 {
                for col in 0..12 {
                    flat[r * 12 + col] = me[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 144);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute WEDGE6 element stiffness matrices.
///
/// coords: shape (n_elem, 18) — each row: 6 nodes × [x,y,z]
/// Returns flat array (n_elem × 324).
#[pyfunction]
pub(crate) fn batch_ke_wedge6<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 6];
            for n in 0..6 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::wedge6_ke(&c, e, nu);
            let mut flat = [0.0f64; 324];
            for r in 0..18 {
                for col in 0..18 {
                    flat[r * 18 + col] = ke[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 324);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute WEDGE6 element mass matrices.
///
/// coords: shape (n_elem, 18)
/// Returns flat array (n_elem × 324).
#[pyfunction]
pub(crate) fn batch_me_wedge6<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();

    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 6];
            for n in 0..6 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::wedge6_me(&c, rho);
            let mut flat = [0.0f64; 324];
            for r in 0..18 {
                for col in 0..18 {
                    flat[r * 18 + col] = me[(r, col)];
                }
            }
            flat
        })
        .collect();

    let mut out = Vec::with_capacity(n_elem * 324);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

// ============================================================================
// Higher-order solid element batch functions
// ============================================================================

/// Batch-compute TETRA10 element stiffness matrices.
///
/// coords: shape (n_elem, 30) — each row: 10 nodes × [x, y, z]
/// Returns flat array (n_elem × 900).
#[pyfunction]
pub(crate) fn batch_ke_tetra10<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 900]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 10];
            for n in 0..10 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::tetra10_ke(&c, e, nu);
            let mut flat = [0.0f64; 900];
            for r in 0..30 { for col in 0..30 { flat[r * 30 + col] = ke[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 900);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute TETRA10 consistent mass matrices.
///
/// coords: shape (n_elem, 30), rho: density
/// Returns flat array (n_elem × 900).
#[pyfunction]
pub(crate) fn batch_me_tetra10<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 900]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 10];
            for n in 0..10 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::tetra10_me(&c, rho);
            let mut flat = [0.0f64; 900];
            for r in 0..30 { for col in 0..30 { flat[r * 30 + col] = me[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 900);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute WEDGE15 element stiffness matrices.
///
/// coords: shape (n_elem, 45) — each row: 15 nodes × [x, y, z]
/// Returns flat array (n_elem × 2025).
#[pyfunction]
pub(crate) fn batch_ke_wedge15<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 2025]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 15];
            for n in 0..15 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::wedge15_ke(&c, e, nu);
            let mut flat = [0.0f64; 2025];
            for r in 0..45 { for col in 0..45 { flat[r * 45 + col] = ke[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 2025);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute WEDGE15 consistent mass matrices.
///
/// coords: shape (n_elem, 45), rho: density
/// Returns flat array (n_elem × 2025).
#[pyfunction]
pub(crate) fn batch_me_wedge15<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 2025]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 15];
            for n in 0..15 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::wedge15_me(&c, rho);
            let mut flat = [0.0f64; 2025];
            for r in 0..45 { for col in 0..45 { flat[r * 45 + col] = me[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 2025);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute HEXA20 element stiffness matrices.
///
/// coords: shape (n_elem, 60) — each row: 20 nodes × [x, y, z]
/// Returns flat array (n_elem × 3600).
#[pyfunction]
pub(crate) fn batch_ke_hexa20<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 3600]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 20];
            for n in 0..20 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::hexa20_ke(&c, e, nu);
            let mut flat = [0.0f64; 3600];
            for r in 0..60 { for col in 0..60 { flat[r * 60 + col] = ke[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 3600);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute HEXA20 consistent mass matrices.
///
/// coords: shape (n_elem, 60), rho: density
/// Returns flat array (n_elem × 3600).
#[pyfunction]
pub(crate) fn batch_me_hexa20<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 3600]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 20];
            for n in 0..20 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::hexa20_me(&c, rho);
            let mut flat = [0.0f64; 3600];
            for r in 0..60 { for col in 0..60 { flat[r * 60 + col] = me[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 3600);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute PYRAMID5 element stiffness matrices.
///
/// coords: shape (n_elem, 15) — each row: 5 nodes × [x, y, z]
/// Returns flat array (n_elem × 225).
#[pyfunction]
pub(crate) fn batch_ke_pyramid5<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 225]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 5];
            for n in 0..5 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::pyramid5_ke(&c, e, nu);
            let mut flat = [0.0f64; 225];
            for r in 0..15 { for col in 0..15 { flat[r * 15 + col] = ke[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 225);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute PYRAMID5 consistent mass matrices.
///
/// coords: shape (n_elem, 15), rho: density
/// Returns flat array (n_elem × 225).
#[pyfunction]
pub(crate) fn batch_me_pyramid5<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 225]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 5];
            for n in 0..5 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::pyramid5_me(&c, rho);
            let mut flat = [0.0f64; 225];
            for r in 0..15 { for col in 0..15 { flat[r * 15 + col] = me[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 225);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute PYRAMID13 element stiffness matrices.
///
/// coords: shape (n_elem, 39) — each row: 13 nodes × [x, y, z]
/// Returns flat array (n_elem × 1521).
#[pyfunction]
pub(crate) fn batch_ke_pyramid13<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 1521]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 13];
            for n in 0..13 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let ke = solid::pyramid13_ke(&c, e, nu);
            let mut flat = [0.0f64; 1521];
            for r in 0..39 { for col in 0..39 { flat[r * 39 + col] = ke[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 1521);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute PYRAMID13 consistent mass matrices.
///
/// coords: shape (n_elem, 39), rho: density
/// Returns flat array (n_elem × 1521).
#[pyfunction]
pub(crate) fn batch_me_pyramid13<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 1521]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let mut c = [[0.0f64; 3]; 13];
            for n in 0..13 {
                c[n][0] = arr[[i, 3 * n]];
                c[n][1] = arr[[i, 3 * n + 1]];
                c[n][2] = arr[[i, 3 * n + 2]];
            }
            let me = solid::pyramid13_me(&c, rho);
            let mut flat = [0.0f64; 1521];
            for r in 0..39 { for col in 0..39 { flat[r * 39 + col] = me[(r, col)]; } }
            flat
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 1521);
    for f in &results { out.extend_from_slice(f); }
    Array1::from(out).into_pyarray(py)
}

// ============================================================================
// QUAD plane element batch functions
// ============================================================================

/// Batch-compute QUAD4 element stiffness matrices (plane strain).
///
/// coords: shape (n_elem, 8) — each row: 4 nodes × [x, y]
/// Returns flat array (n_elem × 64).
#[pyfunction]
pub(crate) fn batch_ke_quad4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 64]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
            ];
            Quad4Precomputed::new(&c).compute_ke_global(e, nu)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 64);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute QUAD4 consistent mass matrices.
///
/// coords: shape (n_elem, 8), rho: density
/// Returns flat array (n_elem × 64).
#[pyfunction]
pub(crate) fn batch_me_quad4<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 64]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
            ];
            Quad4Precomputed::new(&c).compute_me_global(rho)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 64);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute QUAD8 element stiffness matrices (plane strain).
///
/// coords: shape (n_elem, 16) — 8 nodes × [x, y]
/// Returns flat array (n_elem × 256).
#[pyfunction]
pub(crate) fn batch_ke_quad8<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 256]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
                [arr[[i, 8]], arr[[i, 9]]],
                [arr[[i, 10]], arr[[i, 11]]],
                [arr[[i, 12]], arr[[i, 13]]],
                [arr[[i, 14]], arr[[i, 15]]],
            ];
            Quad8Precomputed::new(&c).compute_ke_global(e, nu)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 256);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute QUAD8 consistent mass matrices.
///
/// coords: shape (n_elem, 16), rho: density
/// Returns flat array (n_elem × 256).
#[pyfunction]
pub(crate) fn batch_me_quad8<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 256]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
                [arr[[i, 8]], arr[[i, 9]]],
                [arr[[i, 10]], arr[[i, 11]]],
                [arr[[i, 12]], arr[[i, 13]]],
                [arr[[i, 14]], arr[[i, 15]]],
            ];
            Quad8Precomputed::new(&c).compute_me_global(rho)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 256);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute QUAD9 element stiffness matrices (plane strain).
///
/// coords: shape (n_elem, 18) — 9 nodes × [x, y]
/// Returns flat array (n_elem × 324).
#[pyfunction]
pub(crate) fn batch_ke_quad9<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    e: f64,
    nu: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
                [arr[[i, 8]], arr[[i, 9]]],
                [arr[[i, 10]], arr[[i, 11]]],
                [arr[[i, 12]], arr[[i, 13]]],
                [arr[[i, 14]], arr[[i, 15]]],
                [arr[[i, 16]], arr[[i, 17]]],
            ];
            Quad9Precomputed::new(&c).compute_ke_global(e, nu)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 324);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

/// Batch-compute QUAD9 consistent mass matrices.
///
/// coords: shape (n_elem, 18), rho: density
/// Returns flat array (n_elem × 324).
#[pyfunction]
pub(crate) fn batch_me_quad9<'py>(
    py: Python<'py>,
    coords: PyReadonlyArray2<'py, f64>,
    rho: f64,
) -> Bound<'py, PyArray1<f64>> {
    let arr = coords.as_array();
    let n_elem = arr.nrows();
    let results: Vec<[f64; 324]> = (0..n_elem)
        .into_par_iter()
        .map(|i| {
            let c = [
                [arr[[i, 0]], arr[[i, 1]]],
                [arr[[i, 2]], arr[[i, 3]]],
                [arr[[i, 4]], arr[[i, 5]]],
                [arr[[i, 6]], arr[[i, 7]]],
                [arr[[i, 8]], arr[[i, 9]]],
                [arr[[i, 10]], arr[[i, 11]]],
                [arr[[i, 12]], arr[[i, 13]]],
                [arr[[i, 14]], arr[[i, 15]]],
                [arr[[i, 16]], arr[[i, 17]]],
            ];
            Quad9Precomputed::new(&c).compute_me_global(rho)
        })
        .collect();
    let mut out = Vec::with_capacity(n_elem * 324);
    for r in &results { out.extend_from_slice(r); }
    Array1::from(out).into_pyarray(py)
}

// ============================================================================
// ElementFamily
// ============================================================================

/// Element family classification exposed to Python.
///
/// Mirrors `aeroelast.elements.ElementFamily` (IntEnum).
/// Values are kept identical for drop-in compatibility:
///   SHELL = 2, PLANE = 3, SOLID = 4
#[pyclass(name = "ElementFamily", eq, eq_int, hash, frozen)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum PyElementFamily {
    #[pyo3(name = "SHELL")]
    Shell = 2,
    #[pyo3(name = "PLANE")]
    Plane = 3,
    #[pyo3(name = "SOLID")]
    Solid = 4,
}
