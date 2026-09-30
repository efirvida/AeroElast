//! `PyMeshAssembler` — Python-facing wrapper around the Rust mesh assembler.

use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;

use aeroelast_core::assembly::{ElemType, MaterialSpec, MeshAssembler, MeshTopology};
use aeroelast_core::Material;

use crate::materials::{parse_material, PyLaminate};
use crate::mesh::PyMeshModel;

// ============================================================================
// PyMeshAssembler: PyO3 class wrapping aeroelast_core::assembly::MeshAssembler
// ============================================================================

/// Python-accessible mesh assembler backed by the Rust `MeshAssembler`.
///
/// Construct from flat numpy arrays (matching the existing Python patterns),
/// then call assembly methods that return COO triplets as numpy arrays.
/// The Python shim (`assembler.py`) calls `_coo_to_petsc()` to convert.
#[pyclass]
pub struct PyMeshAssembler {
    inner: MeshAssembler,
}


impl PyMeshAssembler {
    /// Returns a reference to the inner `MeshAssembler` (Rust-only, not exposed to Python).
    pub fn inner(&self) -> &MeshAssembler {
        &self.inner
    }
}

#[pymethods]
impl PyMeshAssembler {
    /// Construct a `PyMeshAssembler` from Python mesh data.
    ///
    /// Parameters
    /// ----------
    /// node_coords : np.ndarray shape (n_nodes, 3)
    ///     Node coordinates in global frame.
    /// connectivity : list[list[int]]
    ///     Per-element node index lists (0-based).
    /// elem_types : list[int]
    ///     Element type code per element:
    ///       3   = MITC3,       4   = MITC4
    ///       33  = MITC3Composite, 44 = MITC4Composite
    ///       104 = QUAD4,       108 = QUAD8,       109 = QUAD9
    /// materials : list[dict]
    ///     Per-element material dicts.  Supported types:
    ///     isotropic   — {type, e, nu, rho, thickness, shear_correction}
    ///     composite   — {type, cm, cb, cs, thickness, e_equiv, mass_per_area, rotational_inertia}
    ///     plane_stress — {type, e, nu, rho, thickness}
    #[new]
    pub fn new(
        node_coords: PyReadonlyArray2<f64>,
        connectivity: Vec<Vec<usize>>,
        elem_types: Vec<u16>,
        materials: Vec<Py<PyAny>>,
        py: Python,
    ) -> PyResult<Self> {
        let coords_arr = node_coords.as_array();
        let n_nodes = coords_arr.nrows();

        // Build flat node_coords vec
        let mut flat_coords = Vec::with_capacity(n_nodes * 3);
        for i in 0..n_nodes {
            flat_coords.push(coords_arr[[i, 0]]);
            flat_coords.push(coords_arr[[i, 1]]);
            flat_coords.push(coords_arr[[i, 2]]);
        }

        // Build elem_types vec
        let rust_elem_types: Vec<ElemType> = elem_types
            .iter()
            .map(|&t| match t {
                3   => Ok(ElemType::Mitc3),
                4   => Ok(ElemType::Mitc4),
                33  => Ok(ElemType::Mitc3Composite),
                44  => Ok(ElemType::Mitc4Composite),
                104 => Ok(ElemType::Quad4),
                108 => Ok(ElemType::Quad8),
                109 => Ok(ElemType::Quad9),
                other => Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown elem_type code {}: see PyMeshAssembler docstring for valid codes", other
                ))),
            })
            .collect::<PyResult<Vec<_>>>()?;

        // Build topology
        let topology = MeshTopology::new(flat_coords, connectivity, rust_elem_types);

        // Parse material specs
        let rust_materials: Vec<MaterialSpec> = materials
            .iter()
            .map(|m| parse_material(py, m))
            .collect::<PyResult<Vec<_>>>()?;

        let inner = MeshAssembler::new(topology, rust_materials);
        Ok(PyMeshAssembler { inner })
    }

    /// Construct a `PyMeshAssembler` from a `MeshModel` already loaded in Rust.
    ///
    /// This is the fast path: topology (node_coords, connectivity, elem_type codes)
    /// is extracted from Rust-native `MeshModel` without any Python loop.
    /// Only the `materials` list still comes from Python.
    ///
    /// Parameters
    /// ----------
    /// mesh : MeshModel
    ///     A `MeshModel` loaded via `MeshModel.from_hdf5(...)`.
    /// materials : list[dict]
    ///     Per-element material dicts in the same element order as `mesh.element_ids()`.
    ///     Same format as `PyMeshAssembler.__init__` (isotropic / composite / ...).
    #[staticmethod]
    pub fn from_mesh_model(
        mesh: &mut PyMeshModel,
        materials: Vec<Py<PyAny>>,
        py: Python,
    ) -> PyResult<Self> {
        // Extract flat coords from Rust MeshModel — no Python loop needed
        let flat_coords = mesh.inner.node_coords_flat().to_vec();

        // Extract connectivity and elem_type codes — fully in Rust
        let (connectivity_usize, code_i32) = mesh.inner.build_connectivity_arrays();

        // Convert i32 codes to ElemType
        let rust_elem_types: Vec<ElemType> = code_i32
            .iter()
            .map(|&t| match t as u16 {
                3   => Ok(ElemType::Mitc3),
                4   => Ok(ElemType::Mitc4),
                33  => Ok(ElemType::Mitc3Composite),
                44  => Ok(ElemType::Mitc4Composite),
                104 => Ok(ElemType::Quad4),
                108 => Ok(ElemType::Quad8),
                109 => Ok(ElemType::Quad9),
                other => Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown elem_type code {other} from MeshModel::build_connectivity_arrays"
                ))),
            })
            .collect::<PyResult<Vec<_>>>()?;

        let topology = MeshTopology::new(flat_coords, connectivity_usize, rust_elem_types);

        // Parse material specs (still from Python)
        let rust_materials: Vec<MaterialSpec> = materials
            .iter()
            .map(|m| parse_material(py, m))
            .collect::<PyResult<Vec<_>>>()?;

        let inner = MeshAssembler::new(topology, rust_materials);
        Ok(PyMeshAssembler { inner })
    }

    /// Build a `PyMeshAssembler` directly from a `PyMeshModel` + Python properties map.
    ///
    /// This replaces `_build_py_mesh_assembler` entirely — no Python loops for
    /// node_coords, connectivity, or angle offsets. Only the properties dict is
    /// still passed from Python.
    ///
    /// Parameters
    /// ----------
    /// mesh : MeshModel
    ///     Rust-native mesh loaded via `MeshModel.from_hdf5(...)`.
    /// properties : dict[str, PyLaminate | dict]
    ///     Mapping from element-set name to property. Supported values:
    ///     - `PyLaminate` instance → composite MITC3/4 element
    ///     - dict with keys `type="isotropic"`, `e`, `nu`, `rho`, `thickness`[, `shear_correction`]
    /// span_direction : list[float] | None
    ///     3-component span direction for fibre angle correction. None to skip.
    /// fallback_material : dict | None
    ///     Isotropic material dict used for elements not in any property set.
    ///     If None, elements without properties raise an error.
    #[staticmethod]
    #[pyo3(signature = (mesh, properties, span_direction=None, fallback_material=None))]
    pub fn from_model(
        mesh: &mut PyMeshModel,
        properties: &pyo3::Bound<'_, pyo3::types::PyDict>,
        span_direction: Option<Vec<f64>>,
        fallback_material: Option<Py<PyAny>>,
        py: Python,
    ) -> PyResult<Self> {
        use aeroelast_mesh::geometry::batch_angle_offsets;
        use aeroelast_core::materials::laminate::Laminate;
        use std::collections::HashMap;

        // ── topology (all in Rust, no Python loops) ──────────────────────────
        let flat_coords = mesh.inner.node_coords_flat().to_vec();
        let (connectivity, code_i32) = mesh.inner.build_connectivity_arrays();

        // ── parse Python properties dict: set_name → MaterialSpec template ───
        // We pre-compute ABD for each (laminate, angle_bucket) pair (the cache)
        // Properties can be PyLaminate or dict
        let mut set_to_prop: HashMap<String, Py<PyAny>> = HashMap::new();
        for (key, val) in properties.iter() {
            let set_name: String = key.extract()?;
            set_to_prop.insert(set_name, val.unbind());
        }

        // ── per-element property lookup: elem_id → (set_name, property) ──────
        // We build this from element_sets in MeshModel (already in Rust).
        // An element may appear in multiple sets (e.g. "allOuterShellEls" and
        // "05_00_HP_LE"). Only sets that have a corresponding property matter;
        // sets without a property are skipped so they don't overwrite a valid
        // mapping with a last-write-wins collision.
        let mut elem_to_set: HashMap<u64, String> = HashMap::new();
        for (set_name, eset) in &mesh.inner.element_sets {
            if !set_to_prop.contains_key(set_name) {
                continue; // skip sets that have no property — avoids last-write-wins clobber
            }
            for &eid in &eset.element_ids {
                elem_to_set.insert(eid, set_name.clone());
            }
        }

        // ── Phase 1: batch angle offsets (Rayon parallel, pure Rust) ─────────
        let span_dir: Option<[f64; 3]> = span_direction.as_ref().and_then(|v| {
            if v.len() >= 3 { Some([v[0], v[1], v[2]]) } else { None }
        });

        // Only compute for composite quad4 elements that have a laminate property
        let mut elem_angle_offsets: HashMap<u64, f64> = HashMap::new();
        if let Some(sd) = span_dir {
            // Identify composite quad4 elements
            let quad_indices: Vec<usize> = mesh.inner.elements.iter().enumerate()
                .filter_map(|(i, elem)| {
                    let set = elem_to_set.get(&elem.id)?;
                    let prop = set_to_prop.get(set)?;
                    // Check if it's a PyLaminate (composite) and has 4 nodes
                    let is_composite = prop.bind(py).is_instance_of::<PyLaminate>();
                    if is_composite && elem.node_ids.len() == 4 { Some(i) } else { None }
                })
                .collect();

            if !quad_indices.is_empty() {
                let quad_conns: Vec<Vec<usize>> = quad_indices.iter()
                    .map(|&i| connectivity[i].clone())
                    .collect();
                let offsets = batch_angle_offsets(&flat_coords, &quad_conns, sd);
                for (i_local, &i_elem) in quad_indices.iter().enumerate() {
                    elem_angle_offsets.insert(mesh.inner.elements[i_elem].id, offsets[i_local]);
                }
            }
        }

        // ── Phase 2a: pre-compute MaterialSpec per (set, angle_bucket) ──────────
        // This avoids PyO3 downcast + GIL ops inside the hot 86k-element loop.
        // We extract all Python objects here (with GIL), compute ABD matrices,
        // and cache the results in a pure-Rust HashMap indexed by (set_name, bucket).
        //
        // Key insight: there are only O(N_sets × N_angle_buckets) unique specs,
        // not O(N_elements). For a blade with 722 sets and typical angle spread,
        // this is ~722 cache misses total instead of ~86k PyO3 ops per run.

        // Helper closure: compute MaterialSpec from a Laminate + angle_bucket_tenths
        let compute_composite_spec = |lam: &Laminate, bucket_tenths: i32| -> PyResult<MaterialSpec> {
            let corrected_lam = if bucket_tenths == 0 {
                lam.clone()
            } else {
                let angle_deg = bucket_tenths as f64 / 10.0;
                let corrected_plies: Vec<_> = lam.plies.iter().map(|p| {
                    aeroelast_core::materials::laminate::Ply {
                        material: p.material.clone(),
                        thickness: p.thickness,
                        angle: p.angle + angle_deg,
                        z_bottom: p.z_bottom,
                        z_top: p.z_top,
                    }
                }).collect();
                Laminate::new(corrected_plies, lam.shear_correction_factor)
                    .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?
            };

            let h = corrected_lam.total_thickness;
            let a_trace = corrected_lam.a[(0,0)] + corrected_lam.a[(1,1)] + corrected_lam.a[(2,2)];
            let e_equiv = if h > 0.0 { a_trace / (3.0 * h) } else { 0.0 };
            let mass_per_area: f64 = corrected_lam.plies.iter()
                .map(|p| p.material.density() * p.thickness)
                .sum();
            let rotational_inertia: f64 = corrected_lam.plies.iter()
                .map(|p| p.material.density() * (p.z_top.powi(3) - p.z_bottom.powi(3)) / 3.0)
                .sum();

            let mut cm = [0.0f64; 9];
            let mut cb_coupling = [0.0f64; 9];  // ADDED
            let mut cb = [0.0f64; 9];
            let mut cs = [0.0f64; 4];
            for ii in 0..3 { for jj in 0..3 {
                cm[ii*3+jj] = corrected_lam.a[(ii, jj)];
                cb_coupling[ii*3+jj] = corrected_lam.b[(ii, jj)];  // ADDED
                cb[ii*3+jj] = corrected_lam.d[(ii, jj)];
            }}
            cs[0] = corrected_lam.cs[(0,0)];
            cs[1] = corrected_lam.cs[(0,1)];
            cs[2] = corrected_lam.cs[(1,0)];
            cs[3] = corrected_lam.cs[(1,1)];

            // ADR-1 (amended): the uncorrected section shear the element
            // consumes. The energy-equivalent `cs` above carries the section's
            // own correction -- the 5/6 of a homogeneous stack -- so it must not
            // be what the element's shear block sees.
            let cs_plain = corrected_lam.shear_stiffness_uncorrected();
            let cs_uncorrected = [
                cs_plain[(0, 0)],
                cs_plain[(0, 1)],
                cs_plain[(1, 0)],
                cs_plain[(1, 1)],
            ];

            Ok(MaterialSpec::Composite {
                cm,
                cb_coupling,
                cb,
                cs,
                thickness: h,
                e_equiv,
                mass_per_area,
                rotational_inertia,
                applied_shear_correction: corrected_lam.applied_shear_correction_factor(),
                cs_uncorrected,
            })
        };

        // Enum to distinguish composite vs isotropic without PyO3 objects in hot path
        #[derive(Clone)]
        enum SetSpec {
            Composite(Laminate),  // raw Laminate for angle-corrected lookup
            Isotropic(MaterialSpec), // pre-parsed, no angle correction needed
        }

        // Pre-compute SetSpec per set (O(N_sets) PyO3 ops, not O(N_elements))
        let mut set_to_setspec: HashMap<String, SetSpec> = HashMap::with_capacity(set_to_prop.len());
        for (set_name, prop) in &set_to_prop {
            let bound = prop.bind(py);
            if let Ok(lam_ref) = bound.downcast::<PyLaminate>() {
                let lam = lam_ref.borrow().inner.clone();
                set_to_setspec.insert(set_name.clone(), SetSpec::Composite(lam));
            } else {
                let spec = parse_material(py, prop)?;
                set_to_setspec.insert(set_name.clone(), SetSpec::Isotropic(spec));
            }
        }

        // Pre-compute fallback spec if provided (avoids per-element parse)
        let fallback_spec: Option<MaterialSpec> = fallback_material.as_ref()
            .map(|fb| parse_material(py, fb))
            .transpose()?;

        // ABD cache: (set_name, angle_bucket_tenths) → MaterialSpec (composite only)
        // bucket_tenths=0 covers the vast majority of elements when no span_direction.
        let mut abd_cache: HashMap<(String, i32), MaterialSpec> = HashMap::new();

        // ── Phase 2b: build MaterialSpec per element (pure Rust, no PyO3) ──────
        let n_elems = mesh.inner.elements.len();
        let mut rust_materials: Vec<MaterialSpec> = Vec::with_capacity(n_elems);
        let mut rust_elem_types: Vec<ElemType> = Vec::with_capacity(n_elems);

        for (i, elem) in mesh.inner.elements.iter().enumerate() {
            let n_nodes_elem = elem.node_ids.len();

            let set_name_opt: Option<&String> = elem_to_set.get(&elem.id);

            let (etype, mat_spec) = if let Some(set_name) = set_name_opt {
                match set_to_setspec.get(set_name) {
                    Some(SetSpec::Composite(lam)) => {
                        // ── Composite shell ───────────────────────────────────
                        let etype = if n_nodes_elem == 3 { ElemType::Mitc3Composite } else { ElemType::Mitc4Composite };

                        let angle_offset = elem_angle_offsets.get(&elem.id).copied().unwrap_or(0.0);
                        let bucket_tenths = (angle_offset * 10.0).round() as i32;

                        // Use (set_name, bucket) as cache key — pure Rust, no PyO3
                        let spec = if let Some(cached) = abd_cache.get(&(set_name.clone(), bucket_tenths)) {
                            cached.clone()
                        } else {
                            let spec = compute_composite_spec(lam, bucket_tenths)?;
                            abd_cache.insert((set_name.clone(), bucket_tenths), spec.clone());
                            spec
                        };
                        (etype, spec)
                    }
                    Some(SetSpec::Isotropic(spec)) => {
                        // ── Isotropic dict ────────────────────────────────────
                        let etype = if n_nodes_elem == 3 { ElemType::Mitc3 } else { ElemType::Mitc4 };
                        (etype, spec.clone())
                    }
                    None => {
                        return Err(pyo3::exceptions::PyValueError::new_err(format!(
                            "element {} belongs to set '{}' which has no property",
                            elem.id, set_name
                        )));
                    }
                }
            } else if let Some(ref fb_spec) = fallback_spec {
                // ── Fallback isotropic ─────────────────────────────────────────
                let etype = if n_nodes_elem == 3 { ElemType::Mitc3 } else { ElemType::Mitc4 };
                (etype, fb_spec.clone())
            } else {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "element {} has no property and no fallback_material was provided",
                    elem.id
                )));
            };

            // Override elem_type from mesh entities (e.g. plane element codes)
            let final_etype = {
                let mesh_code = code_i32[i];
                match mesh_code as u16 {
                    3   => ElemType::Mitc3,
                    4   => ElemType::Mitc4,
                    33  => ElemType::Mitc3Composite,
                    44  => ElemType::Mitc4Composite,
                    104 => ElemType::Quad4,
                    108 => ElemType::Quad8,
                    109 => ElemType::Quad9,
                    // For shell elements loaded from HDF5 without composite info,
                    // use the composite/isotropic determination from property
                    _ => etype,
                }
            };

            rust_elem_types.push(final_etype);
            rust_materials.push(mat_spec);
        }
        let topology = MeshTopology::new(flat_coords, connectivity, rust_elem_types);
        let inner = MeshAssembler::new(topology, rust_materials);
        Ok(PyMeshAssembler { inner })
    }

    /// Total number of DOFs in the system.
    #[getter]
    pub fn dofs_count(&self) -> usize {
        self.inner.dofs_count
    }

    #[getter]
    pub fn n_elems(&self) -> usize {
        self.inner.topology.n_elems
    }

    #[getter]
    pub fn n_nodes(&self) -> usize {
        self.inner.topology.n_nodes
    }

    /// Assemble the global elastic stiffness matrix K.
    ///
    /// Returns (rows, cols, vals) as numpy int64/float64 arrays.
    pub fn assemble_k<'py>(
        &self,
        py: Python<'py>,
    ) -> (
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    ) {
        let (rows, cols, vals) = self.inner.assemble_k();
        (
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        )
    }

    /// Assemble the global consistent mass matrix M.
    ///
    /// Returns (rows, cols, vals) as numpy int64/float64 arrays.
    pub fn assemble_m<'py>(
        &self,
        py: Python<'py>,
    ) -> (
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    ) {
        let (rows, cols, vals) = self.inner.assemble_m();
        (
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        )
    }

    /// Assemble the row-sum lumped mass diagonal.
    ///
    /// Returns a 1-D float64 array of length `dofs_count` containing the
    /// diagonal entries after row-sum lumping and the standard mass floor.
    pub fn assemble_m_lumped<'py>(&self, py: Python<'py>) -> pyo3::Bound<'py, PyArray1<f64>> {
        Array1::from(self.inner.assemble_m_lumped()).into_pyarray(py)
    }

    /// Compute the total model mass from material properties and geometry.
    ///
    /// For shell elements (MITC3/MITC4) uses the precomputed element area.
    /// For plane/solid elements uses the partition-of-unity property of the
    /// consistent mass matrix: ``Σ_ij (Me_x)_ij = m_elem``.
    ///
    /// Returns
    /// -------
    /// float
    ///     Total mass in kg (or consistent units of the model).
    pub fn total_elemental_mass(&self) -> f64 {
        self.inner.total_elemental_mass()
    }

    /// Assemble the global body load vector.
    ///
    /// Parameters
    /// ----------
    /// gravity : list[f64] of length 3 — body acceleration [gx, gy, gz] m/s²
    ///
    /// Returns np.ndarray of length dofs_count.
    pub fn assemble_f_body<'py>(
        &self,
        py: Python<'py>,
        gravity: [f64; 3],
    ) -> pyo3::Bound<'py, PyArray1<f64>> {
        let f = self.inner.assemble_f_body(gravity);
        Array1::from(f).into_pyarray(py)
    }

    /// Assemble the global geometric stiffness matrix K_σ.
    ///
    /// Parameters
    /// ----------
    /// sigma : np.ndarray shape (n_elems, 3) — membrane stress [σxx, σyy, σxy] per element
    ///
    /// Returns (rows, cols, vals) as numpy int64/float64 arrays.
    pub fn assemble_geometric_k<'py>(
        &self,
        py: Python<'py>,
        sigma: PyReadonlyArray2<f64>,
    ) -> PyResult<(
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    )> {
        let sigma_arr = sigma.as_array();
        let n_elems = sigma_arr.nrows();
        if n_elems != self.inner.topology.n_elems {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "sigma must have shape (n_elems, 3); got ({}, {})",
                n_elems, sigma_arr.ncols()
            )));
        }
        let sigma_vecs: Vec<[f64; 3]> = (0..n_elems)
            .map(|e| [sigma_arr[[e, 0]], sigma_arr[[e, 1]], sigma_arr[[e, 2]]])
            .collect();
        let (rows, cols, vals) = self.inner.assemble_geometric_k(&sigma_vecs);
        Ok((
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        ))
    }

    /// Assemble the nonlinear tangent stiffness matrix K_T(u).
    ///
    /// Parameters
    /// ----------
    /// u : np.ndarray shape (dofs_count,) — global displacement vector
    ///
    /// Returns (rows, cols, vals) as numpy int64/float64 arrays.
    pub fn assemble_kt<'py>(
        &self,
        py: Python<'py>,
        u: PyReadonlyArray1<f64>,
    ) -> PyResult<(
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    )> {
        let u_slice = u.as_slice()?;
        let (rows, cols, vals) = self.inner.assemble_kt(u_slice);
        Ok((
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        ))
    }

    /// Assemble the global internal force vector f_int(u).
    ///
    /// Parameters
    /// ----------
    /// u : np.ndarray shape (dofs_count,) — global displacement vector
    /// nonlinear : bool — include geometric-nonlinear contributions
    ///
    /// Returns np.ndarray of length dofs_count.
    pub fn assemble_fint<'py>(
        &self,
        py: Python<'py>,
        u: PyReadonlyArray1<f64>,
        nonlinear: bool,
    ) -> PyResult<pyo3::Bound<'py, PyArray1<f64>>> {
        let u_slice = u.as_slice()?;
        let f = self.inner.assemble_fint(u_slice, nonlinear);
        Ok(Array1::from(f).into_pyarray(py))
    }

    /// Update the reference (undeformed) configuration for Updated-Lagrangian
    /// incremental nonlinear analysis.
    ///
    /// Call after each converged load step. `u_inc` is the incremental
    /// displacement vector (length = dofs_count). Only the translational DOFs
    /// (indices 0,1,2 of each 6-DOF node block) update the node coordinates;
    /// the precomputed element data is rebuilt from the new geometry.
    pub fn update_reference(
        &mut self,
        u_inc: PyReadonlyArray1<f64>,
    ) -> PyResult<()> {
        let u_slice = u_inc.as_slice()?;
        self.inner.update_reference(u_slice);
        Ok(())
    }

    /// NNZ per row for PETSc preallocation.
    ///
    /// Returns np.ndarray of length dofs_count (int64).
    pub fn nnz_per_row<'py>(&self, py: Python<'py>) -> pyo3::Bound<'py, PyArray1<i64>> {
        Array1::from(self.inner.nnz_per_row().to_vec()).into_pyarray(py)
    }

    /// Return mass-per-area (kg/m²) for each element.
    ///
    /// For composite elements this is Σ(ρ_k · t_k) over all plies.
    /// For isotropic shell elements this is ρ (kg/m³) — the caller is
    /// responsible for multiplying by thickness when needed.
    ///
    /// Returns
    /// -------
    /// np.ndarray shape (n_elems,), dtype float64
    pub fn rho_per_elem<'py>(&self, py: Python<'py>) -> pyo3::Bound<'py, PyArray1<f64>> {
        let rho: Vec<f64> = self.inner.materials.iter().map(|m| match m {
            MaterialSpec::Composite { mass_per_area, .. } => *mass_per_area,
            MaterialSpec::Isotropic { rho, .. } => *rho,
            _ => 0.0,
        }).collect();
        Array1::from(rho).into_pyarray(py)
    }

    /// Assemble the geometric stiffness matrix K_σ from centrifugal loading.
    ///
    /// Computes centrifugal prestress per element from first principles and
    /// assembles the geometric stiffness matrix in a single Rust pass.
    ///
    /// Parameters
    /// ----------
    /// omega : float
    ///     Angular velocity magnitude (rad/s).
    /// rotation_axis : list[float] of length 3
    ///     Unit vector of the rotation axis [ax, ay, az].
    /// rotation_center : list[float] of length 3
    ///     A point on the rotation axis [cx, cy, cz] (m).
    /// rho_per_elem : np.ndarray shape (n_elems,)
    ///     Density (kg/m³) or mass-per-area (kg/m²) per element.
    ///     For isotropic shell elements use ρ (kg/m³);
    ///     for composite laminates use mass_per_area (kg/m²).
    ///
    /// Returns (rows, cols, vals) as numpy int64/float64 arrays.
    pub fn assemble_centrifugal_k<'py>(
        &self,
        py: Python<'py>,
        omega: f64,
        rotation_axis: [f64; 3],
        rotation_center: [f64; 3],
        rho_per_elem: PyReadonlyArray1<f64>,
    ) -> PyResult<(
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    )> {
        let rho_slice = rho_per_elem.as_slice()?;
        if rho_slice.len() != self.inner.topology.n_elems {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "rho_per_elem must have length n_elems={}; got {}",
                self.inner.topology.n_elems,
                rho_slice.len()
            )));
        }
        let (rows, cols, vals) = self.inner.assemble_centrifugal_k(
            omega,
            rotation_axis,
            rotation_center,
            rho_slice,
        );
        Ok((
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        ))
    }

    /// Compute element-centroid stress and strain for every element.
    ///
    /// Parameters
    /// ----------
    /// u : np.ndarray shape (dofs_count,) — global displacement vector
    /// z_factor : float — non-dimensional through-thickness location for shell
    ///     elements: z = z_factor × h.  Typical values: -0.5 (bottom),
    ///     0.0 (mid-surface, default), +0.5 (top).  Ignored for plane and
    ///     solid elements.
    /// stress_type : int — shell stress contribution (ignored for non-shell):
    ///     0 = membrane only, 1 = bending only, 2 = total (default).
    ///
    /// Returns
    /// -------
    /// (sigma, epsilon) : tuple of np.ndarray, each shape (n_elems, 6)
    ///     Voigt [σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_zx] per element.
    ///     Out-of-plane entries are zero for shell / plane elements.
    pub fn compute_stress_field<'py>(
        &self,
        py: Python<'py>,
        u: PyReadonlyArray1<f64>,
        z_factor: f64,
        stress_type: u8,
    ) -> PyResult<(
        pyo3::Bound<'py, PyArray2<f64>>,
        pyo3::Bound<'py, PyArray2<f64>>,
    )> {
        let u_slice = u.as_slice()?;
        let (sigma, epsilon) = self.inner.compute_stress_field(u_slice, z_factor, stress_type);
        let n = sigma.len();
        // Flatten into (n_elems, 6) arrays
        let sigma_flat: Vec<f64> = sigma.into_iter().flat_map(|s| s.into_iter()).collect();
        let eps_flat: Vec<f64> = epsilon.into_iter().flat_map(|e| e.into_iter()).collect();
        let sigma_arr = Array2::from_shape_vec((n, 6), sigma_flat)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?
            .into_pyarray(py);
        let eps_arr = Array2::from_shape_vec((n, 6), eps_flat)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?
            .into_pyarray(py);
        Ok((sigma_arr, eps_arr))
    }

    /// Replace the mesh node coordinates in place.
    ///
    /// Used by the S-4/S-5 rotating-frame checks that re-assemble stiffness and
    /// mass matrices on a rotated geometry without rebuilding the assembler.
    ///
    /// Parameters
    /// ----------
    /// coords : np.ndarray shape (n_nodes, 3)
    ///     New nodal coordinates, in the same node order as the original mesh.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the number of coordinates does not match the mesh node count.
    #[pyo3(name = "update_node_coordinates")]
    pub fn update_node_coordinates(
        &mut self,
        coords: PyReadonlyArray2<f64>,
    ) -> PyResult<()> {
        let n_nodes = self.inner.topology.n_nodes;
        let arr = coords.as_array();
        let rows = arr.nrows();
        let cols = arr.ncols();
        if rows != n_nodes || cols != 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "update_node_coordinates: expected shape ({n_nodes}, 3), got ({rows}, {cols})",
            )));
        }
        let flat: Vec<f64> = arr.iter().copied().collect();
        self.inner.update_node_coordinates(&flat);
        Ok(())
    }

    /// Nodal centrifugal load vector for the pre-stress static solve.
    ///
    /// Returns a np.ndarray of shape (dofs_count,) with the lumped nodal
    /// centrifugal forces `rho_A * omega^2 * r` (radial direction).
    ///
    /// This is the right-hand side of `K u = f_cf`, the static pre-stress
    /// problem whose displacement solution feeds
    /// `assemble_geometric_k_from_disp`.  There is no point in
    /// `assemble_centrifugal_k` alone: on a rotating blade the old local
    /// approximation `sigma = rho*omega^2*r*l_char` was two to three orders of
    /// magnitude too small, so rotating modes came out parked.
    pub fn centrifugal_load<'py>(
        &self,
        py: Python<'py>,
        omega: f64,
        rotation_axis: [f64; 3],
        rotation_center: [f64; 3],
        rho_per_elem: PyReadonlyArray1<f64>,
    ) -> PyResult<pyo3::Bound<'py, PyArray1<f64>>> {
        let rho_slice = rho_per_elem.as_slice()?;
        if rho_slice.len() != self.inner.topology.n_elems {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "rho_per_elem must have length n_elems={}; got {}",
                self.inner.topology.n_elems,
                rho_slice.len()
            )));
        }
        let f = self
            .inner
            .centrifugal_load(omega, rotation_axis, rotation_center, rho_slice);
        Ok(Array1::from(f).into_pyarray(py))
    }

    /// Geometric stiffness K_sigma from a full displacement field.
    ///
    /// Recovers membrane force resultants at each shell element centroid and
    /// assembles the geometric stiffness.  Use with the displacement solution
    /// of the centrifugal pre-stress static solve `K u = f_cf`.
    pub fn assemble_geometric_k_from_disp<'py>(
        &self,
        py: Python<'py>,
        u: PyReadonlyArray1<f64>,
    ) -> PyResult<(
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<i64>>,
        pyo3::Bound<'py, PyArray1<f64>>,
    )> {
        let u_slice = u.as_slice()?;
        if u_slice.len() != self.inner.dofs_count {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "u must have length dofs_count={}; got {}",
                self.inner.dofs_count,
                u_slice.len()
            )));
        }
        let (rows, cols, vals) = self.inner.assemble_geometric_k_from_disp(u_slice);
        Ok((
            Array1::from(rows).into_pyarray(py),
            Array1::from(cols).into_pyarray(py),
            Array1::from(vals).into_pyarray(py),
        ))
    }
}
