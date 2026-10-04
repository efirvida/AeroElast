//! Batch element kernels exposed to Python (MITC3/MITC4 shells and quad plane elements).

use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use rayon::prelude::*;

use aeroelast_core::elements::mitc3::{self, Mitc3Precomputed};
use aeroelast_core::elements::mitc4::{self, Mitc4Precomputed};
use aeroelast_core::elements::smoothing::{self, TriangleFrame};
use aeroelast_core::elements::quad::{Quad4Precomputed, Quad8Precomputed, Quad9Precomputed};
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
// Strain-smoothed MITC3+ (Lee & Lee 2019)
// ============================================================================

/// Assemble the global stiffness of a triangular mesh of strain-smoothed MITC3+
/// shell elements (Lee & Lee 2019, *Computers and Structures* 223:106096).
///
/// Each element's membrane strain is smoothed with its three edge neighbours
/// (Eqs. 15-18), so the element stiffness lives over the six-node union layout of
/// [`mitc3::compute_ke_union_smoothed`] and is rotated to global coordinates with
/// each owning element's own frame.  Bending, transverse shear and drilling stay
/// on the element's own three nodes, exactly as the paper prescribes.
///
/// Parameters
/// ----------
/// node_coords : np.ndarray (n_nodes, 3)
/// connectivity : list[list[int]]
///     One 3-node triangle per element, 0-based, consistently oriented.
/// e_mod, nu, thickness : float
/// shear_correction : float, default 5/6
///
/// Returns (rows, cols, vals) as numpy int64/int64/float64 COO triplets.
#[pyfunction]
#[pyo3(signature = (node_coords, connectivity, e_mod, nu, thickness, shear_correction=5.0/6.0))]
pub(crate) fn assemble_smoothed_mitc3<'py>(
    py: Python<'py>,
    node_coords: PyReadonlyArray2<'py, f64>,
    connectivity: Vec<Vec<usize>>,
    e_mod: f64,
    nu: f64,
    thickness: f64,
    shear_correction: f64,
) -> PyResult<(
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let coords = node_coords.as_array();
    let n_elem = connectivity.len();

    let material = IsotropicMaterial::new(e_mod, nu, 0.0);
    let constitutive = material.constitutive(thickness, shear_correction);

    // Precompute every triangle, then its smoothing frame.
    let mut elements = Vec::with_capacity(n_elem);
    let mut frames = Vec::with_capacity(n_elem);
    for nodes in &connectivity {
        if nodes.len() != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "strain-smoothed MITC3 needs 3-node triangle connectivity",
            ));
        }
        let mut c = [0.0f64; 9];
        for (local, &node) in nodes.iter().enumerate() {
            for axis in 0..3 {
                c[3 * local + axis] = coords[[node, axis]];
            }
        }
        let pre = Mitc3Precomputed::new(&c, constitutive.clone(), thickness, e_mod, 1.0);
        frames.push(TriangleFrame {
            j_mat: pre.j_mat,
            area: pre.area,
            normal: nalgebra::Vector3::new(pre.t3[(2, 0)], pre.t3[(2, 1)], pre.t3[(2, 2)]),
            frame: pre.t3,
        });
        elements.push(pre);
    }

    let neighbours = smoothing::edge_neighbours(&connectivity);
    let operators = smoothing::smoothed_membrane_strain(&frames, &neighbours);

    let mut rows: Vec<i64> = Vec::new();
    let mut cols: Vec<i64> = Vec::new();
    let mut vals: Vec<f64> = Vec::new();

    for element in 0..n_elem {
        let own = &elements[element];

        // Union slots: the target's three nodes, then each neighbour's unique node.
        let mut slots = [[0usize, 1, 2]; 4];
        let mut union_node_ids = [
            connectivity[element][0],
            connectivity[element][1],
            connectivity[element][2],
            0,
            0,
            0,
        ];
        let mut entries: [Option<&Mitc3Precomputed>; 4] = [Some(own), None, None, None];
        let mut slot_frames = [own.t3; mitc3::SMOOTHED_UNION_NODES];

        for edge in 0..3 {
            let Some(neighbour) = neighbours[element][edge] else {
                continue;
            };
            entries[1 + edge] = Some(&elements[neighbour]);
            slot_frames[3 + edge] = elements[neighbour].t3;

            // Shared nodes map to the target's local position; the neighbour's
            // remaining node becomes the union's unique node for this edge.
            let a = connectivity[element][edge];
            let b = connectivity[element][(edge + 1) % 3];
            let mut unique = None;
            let mut mapping = [0usize; 3];
            for (local, &node) in connectivity[neighbour].iter().enumerate() {
                if node == a {
                    mapping[local] = edge;
                } else if node == b {
                    mapping[local] = (edge + 1) % 3;
                } else {
                    mapping[local] = 3 + edge;
                    unique = Some(local);
                }
            }
            slots[1 + edge] = mapping;
            if let Some(unique) = unique {
                union_node_ids[3 + edge] = connectivity[neighbour][unique];
            }
        }

        let bm_union =
            mitc3::smoothed_membrane_b(own, &entries, &slots, &operators[element].weights);
        let k_local = mitc3::compute_ke_union_smoothed(own, &bm_union);
        let k_global = mitc3::transform_union_to_global(&k_local, &slot_frames);

        for row in 0..mitc3::SMOOTHED_UNION_DOFS {
            let global_row = 6 * union_node_ids[row / 6] + (row % 6);
            for col in 0..mitc3::SMOOTHED_UNION_DOFS {
                let value = k_global[(row, col)];
                if value == 0.0 {
                    continue;
                }
                let global_col = 6 * union_node_ids[col / 6] + (col % 6);
                rows.push(global_row as i64);
                cols.push(global_col as i64);
                vals.push(value);
            }
        }
    }

    Ok((
        Array1::from(rows).into_pyarray(py),
        Array1::from(cols).into_pyarray(py),
        Array1::from(vals).into_pyarray(py),
    ))
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
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, shear_correction);
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
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, shear_correction);
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
///
/// Uses the faithful total-Lagrangian tangent of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eq. (24a) (`n_gamma_kt_global`).
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
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, shear_correction);

            let mut u = mitc4::Vec24::zeros();
            for i in 0..24 {
                u[i] = disp_arr[[e, i]];
            }

            // Faithful total-Lagrangian tangent, Eq. (24a). KNOWN LIMITATION:
            // `n_gamma_kt_global` builds its `GlCurrentState` from the INITIAL
            // element geometry, so this is the first-step
            // (reference-configuration) linearisation; threading the
            // last-converged state belongs to the solver and is a follow-up.
            let kt = mitc4::n_gamma_kt_global(&pre, &u);
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
///
/// `nonlinear = true` uses the faithful total-Lagrangian internal force of
/// Ko, Lee & Bathe (2017), C&S 185:1-14, Eq. (24b) (`n_gamma_fint_global`);
/// `nonlinear = false` keeps the linear path `K_0 u`
/// (`compute_fint_global(.., false)`).
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
            let pre = Mitc4Precomputed::new(&node_coords, constitutive.clone(), thickness, shear_correction);

            let mut u = mitc4::Vec24::zeros();
            for i in 0..24 {
                u[i] = disp_arr[[e, i]];
            }

            // KNOWN LIMITATION: `n_gamma_fint_global` builds its
            // `GlCurrentState` from the INITIAL element geometry, so the
            // nonlinear branch is the first-step (reference-configuration)
            // linearisation; threading the last-converged state belongs to the
            // solver and is a follow-up.
            let f = if nonlinear {
                mitc4::n_gamma_fint_global(&pre, &u)
            } else {
                mitc4::compute_fint_global(&pre, &u, false)
            };
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
    let _e_arr = e_equiv.as_array();
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

            let constitutive = composite_constitutive(&a, &b, &d, &cs, h);
            // `e_equiv` has no MITC4+/D constructor argument (out-of-scope PyO3 change).
            let pre = Mitc4Precomputed::new(&node_coords, constitutive, h, 1.0);
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
            let pre = Mitc4Precomputed::new(&node_coords, dummy, 1.0, 1.0);
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
///   SHELL = 2, PLANE = 3
#[pyclass(name = "ElementFamily", eq, eq_int, hash, frozen)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum PyElementFamily {
    #[pyo3(name = "SHELL")]
    Shell = 2,
    #[pyo3(name = "PLANE")]
    Plane = 3,
}
