//! Material parsing helpers and Python-facing material types.

use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use aeroelast_core::assembly::MaterialSpec;
use aeroelast_core::materials::laminate::{Laminate, Ply};
use aeroelast_core::materials::orthotropic::OrthotropicMaterial;
use aeroelast_core::Material;

/// Parse a Python material dict into a `MaterialSpec`.
///
/// Expected keys (isotropic): type="isotropic", e, nu, rho, thickness, shear_correction, drilling_scale(optional)
/// Expected keys (composite): type="composite", cm=[9], cb=[9], cs=[4],
///                            thickness, e_equiv, mass_per_area, rotational_inertia
/// Expected keys (plane_stress): type="plane_stress", e, nu, rho, thickness
/// Expected keys (solid_3d): type="solid_3d", e, nu, rho
pub(crate) fn parse_material(py: Python, obj: &Py<PyAny>) -> PyResult<MaterialSpec> {
    let dict = obj.bind(py);
    let mat_type: String = dict.get_item("type")
        .map_err(|_| pyo3::exceptions::PyKeyError::new_err("material dict missing 'type'"))?
        .extract()?;

    match mat_type.as_str() {
        "isotropic" => {
            let e: f64 = dict.get_item("e")?.extract()?;
            let nu: f64 = dict.get_item("nu")?.extract()?;
            let rho: f64 = dict.get_item("rho")?.extract()?;
            let thickness: f64 = dict.get_item("thickness")?.extract()?;
            let shear_correction: f64 = dict.get_item("shear_correction")
                .ok()
                .and_then(|v| v.extract::<f64>().ok())
                .unwrap_or(5.0 / 6.0);
            let drilling_scale: f64 = dict.get_item("drilling_scale")
                .ok()
                .and_then(|v| v.extract::<f64>().ok())
                .unwrap_or(1.0);
            Ok(MaterialSpec::Isotropic {
                e,
                nu,
                rho,
                thickness,
                shear_correction,
                drilling_scale,
            })
        }
        "composite" => {
            let cm_list: Vec<f64> = dict.get_item("cm")?.extract()?;
            let cb_list: Vec<f64> = dict.get_item("cb")?.extract()?;
            let cs_list: Vec<f64> = dict.get_item("cs")?.extract()?;
            let b_list: Vec<f64> = dict.get_item("b_coupling")?.extract()?;  // ADDED
            let thickness: f64 = dict.get_item("thickness")?.extract()?;
            let e_equiv: f64 = dict.get_item("e_equiv")?.extract()?;
            let mass_per_area: f64 = dict.get_item("mass_per_area")?.extract()?;
            let rotational_inertia: f64 = dict.get_item("rotational_inertia")?.extract()?;

            if cm_list.len() != 9 || cb_list.len() != 9 || cs_list.len() != 4 || b_list.len() != 9 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "composite material: cm must have 9 elements, b_coupling 9, cb 9, cs 4"
                ));
            }
            let mut cm = [0.0f64; 9];
            let mut cb_coupling = [0.0f64; 9];  // ADDED
            let mut cb = [0.0f64; 9];
            let mut cs = [0.0f64; 4];
            cm.copy_from_slice(&cm_list);
            cb_coupling.copy_from_slice(&b_list);  // ADDED
            cb.copy_from_slice(&cb_list);
            cs.copy_from_slice(&cs_list);

            Ok(MaterialSpec::Composite { cm, cb_coupling, cb, cs, thickness, e_equiv, mass_per_area, rotational_inertia })
        }
        "plane_stress" => {
            let e: f64 = dict.get_item("e")?.extract()?;
            let nu: f64 = dict.get_item("nu")?.extract()?;
            let rho: f64 = dict.get_item("rho")?.extract()?;
            let thickness: f64 = dict.get_item("thickness")
                .ok()
                .and_then(|v| v.extract::<f64>().ok())
                .unwrap_or(1.0);
            Ok(MaterialSpec::PlaneStress { e, nu, rho, thickness })
        }
        "solid_3d" => {
            let e: f64 = dict.get_item("e")?.extract()?;
            let nu: f64 = dict.get_item("nu")?.extract()?;
            let rho: f64 = dict.get_item("rho")?.extract()?;
            Ok(MaterialSpec::Solid3D { e, nu, rho })
        }
        other => Err(pyo3::exceptions::PyValueError::new_err(format!(
            "unknown material type '{}', expected 'isotropic', 'composite', 'plane_stress', or 'solid_3d'", other
        ))),
    }
}

// ============================================================================
// Python-exposed material types: OrthotropicMaterial, Ply, Laminate
// ============================================================================

// ============================================================================
// PyOrthotropicMaterial, PyPly, PyLaminate
// ============================================================================

/// Python wrapper for `OrthotropicMaterial`.
///
/// Mirrors the Python `OrthotropicMaterial` dataclass from `material.py`.
#[pyclass(name = "OrthotropicMaterial")]
#[derive(Clone)]
pub(crate) struct PyOrthotropicMaterial {
    inner: OrthotropicMaterial,
}

#[pymethods]
impl PyOrthotropicMaterial {
    #[new]
    #[pyo3(signature = (e1, e2, e3, g12, g23, g13, nu12, nu23, nu31, rho))]
    fn new(
        e1: f64, e2: f64, e3: f64,
        g12: f64, g23: f64, g13: f64,
        nu12: f64, nu23: f64, nu31: f64,
        rho: f64,
    ) -> Self {
        Self {
            inner: OrthotropicMaterial::new(e1, e2, e3, g12, g23, g13, nu12, nu23, nu31, rho),
        }
    }

    #[getter] fn e1(&self) -> f64 { self.inner.e1 }
    #[getter] fn e2(&self) -> f64 { self.inner.e2 }
    #[getter] fn e3(&self) -> f64 { self.inner.e3 }
    #[getter] fn g12(&self) -> f64 { self.inner.g12 }
    #[getter] fn g23(&self) -> f64 { self.inner.g23 }
    #[getter] fn g13(&self) -> f64 { self.inner.g13 }
    #[getter] fn nu12(&self) -> f64 { self.inner.nu12 }
    #[getter] fn nu23(&self) -> f64 { self.inner.nu23 }
    #[getter] fn nu31(&self) -> f64 { self.inner.nu31 }
    #[getter] fn rho(&self) -> f64 { self.inner.rho }

    fn __repr__(&self) -> String {
        format!(
            "OrthotropicMaterial(E1={:.3e}, E2={:.3e}, G12={:.3e}, nu12={:.3}, rho={:.1})",
            self.inner.e1, self.inner.e2, self.inner.g12, self.inner.nu12, self.inner.rho
        )
    }
}


/// Python wrapper for a single `Ply`.
#[pyclass(name = "Ply")]
#[derive(Clone)]
pub(crate) struct PyPly {
    inner: Ply,
}

#[pymethods]
impl PyPly {
    #[new]
    #[pyo3(signature = (material, thickness, angle))]
    fn new(material: &PyOrthotropicMaterial, thickness: f64, angle: f64) -> Self {
        Self { inner: Ply::new(material.inner, thickness, angle) }
    }

    #[getter] fn thickness(&self) -> f64 { self.inner.thickness }
    #[getter] fn angle(&self) -> f64 { self.inner.angle }
    #[getter] fn z_bottom(&self) -> f64 { self.inner.z_bottom }
    #[getter] fn z_top(&self) -> f64 { self.inner.z_top }

    fn __repr__(&self) -> String {
        format!("Ply(angle={:.1}°, t={:.4e}m)", self.inner.angle, self.inner.thickness)
    }
}


/// Python wrapper for `Laminate` (CLT — Classical Lamination Theory).
///
/// Computes A, B, D and Cs matrices from the ply stack.
#[pyclass(name = "Laminate")]
pub(crate) struct PyLaminate {
    pub(crate) inner: Laminate,
}

#[pymethods]
impl PyLaminate {
    #[new]
    #[pyo3(signature = (plies, shear_correction_factor=0.75))]
    fn new(plies: Vec<PyPly>, shear_correction_factor: f64) -> PyResult<Self> {
        let rust_plies: Vec<Ply> = plies.into_iter().map(|p| p.inner).collect();
        let lam = Laminate::new(rust_plies, shear_correction_factor)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;
        Ok(Self { inner: lam })
    }

    /// 3×3 extensional stiffness matrix A [N/m] — row-major flat array (9 values).
    fn a_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut arr = Array2::<f64>::zeros((3, 3));
        for i in 0..3 {
            for j in 0..3 {
                arr[[i, j]] = self.inner.a[(i, j)];
            }
        }
        arr.into_pyarray(py)
    }

    /// 3×3 coupling stiffness matrix B [N] — row-major flat array (9 values).
    fn b_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut arr = Array2::<f64>::zeros((3, 3));
        for i in 0..3 {
            for j in 0..3 {
                arr[[i, j]] = self.inner.b[(i, j)];
            }
        }
        arr.into_pyarray(py)
    }

    /// 3×3 bending stiffness matrix D [N·m] — row-major flat array (9 values).
    fn d_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut arr = Array2::<f64>::zeros((3, 3));
        for i in 0..3 {
            for j in 0..3 {
                arr[[i, j]] = self.inner.d[(i, j)];
            }
        }
        arr.into_pyarray(py)
    }

    /// 2×2 transverse shear stiffness matrix Cs [N/m].
    fn cs_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let mut arr = Array2::<f64>::zeros((2, 2));
        for i in 0..2 {
            for j in 0..2 {
                arr[[i, j]] = self.inner.cs[(i, j)];
            }
        }
        arr.into_pyarray(py)
    }

    /// Full 6×6 ABD matrix [[A, B], [B, D]] as a numpy array.
    fn abd_matrix<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let flat = self.inner.abd_matrix_flat();
        let mut arr = Array2::<f64>::zeros((6, 6));
        for i in 0..6 {
            for j in 0..6 {
                arr[[i, j]] = flat[i * 6 + j];
            }
        }
        arr.into_pyarray(py)
    }

    #[getter] fn total_thickness(&self) -> f64 { self.inner.total_thickness }
    #[getter] fn n_plies(&self) -> usize { self.inner.plies.len() }
    #[getter] fn shear_correction_factor(&self) -> f64 { self.inner.shear_correction_factor }
    #[getter] fn is_symmetric(&self) -> bool { self.inner.is_symmetric() }
    #[getter] fn is_balanced(&self) -> bool { self.inner.is_balanced() }

    /// Areal density of the laminate [kg/m²] — sum of (ply density × ply thickness).
    ///
    /// Divide by ``total_thickness`` to get the equivalent volumetric density [kg/m³]
    /// needed for CalculiX *DENSITY cards.
    #[getter]
    fn areal_density(&self) -> f64 {
        self.inner.plies.iter()
            .map(|p| p.material.density() * p.thickness)
            .sum()
    }

    /// List of plies (bottom → top) as dicts.
    ///
    /// Each dict has keys: ``thickness``, ``angle``, ``z_bottom``, ``z_top``
    /// and ``material`` (sub-dict with e1, e2, e3, g12, g23, g13, nu12, nu23, nu31, rho).
    #[getter]
    fn plies(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let list = pyo3::types::PyList::empty(py);
        for ply in &self.inner.plies {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("thickness", ply.thickness)?;
            dict.set_item("angle",     ply.angle)?;
            dict.set_item("z_bottom",  ply.z_bottom)?;
            dict.set_item("z_top",     ply.z_top)?;
            let mdict = pyo3::types::PyDict::new(py);
            mdict.set_item("e1",   ply.material.e1)?;
            mdict.set_item("e2",   ply.material.e2)?;
            mdict.set_item("e3",   ply.material.e3)?;
            mdict.set_item("g12",  ply.material.g12)?;
            mdict.set_item("g23",  ply.material.g23)?;
            mdict.set_item("g13",  ply.material.g13)?;
            mdict.set_item("nu12", ply.material.nu12)?;
            mdict.set_item("nu23", ply.material.nu23)?;
            mdict.set_item("nu31", ply.material.nu31)?;
            mdict.set_item("rho",  ply.material.rho)?;
            dict.set_item("material", mdict)?;
            list.append(dict)?;
        }
        Ok(list.unbind().into_any())
    }

    /// Backward/interop helper: explicit method form of ``plies``.
    ///
    /// Some call sites can prefer a method instead of a property lookup.
    fn plies_data(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.plies(py)
    }

    fn __repr__(&self) -> String {
        let angles: Vec<String> = self.inner.plies.iter()
            .map(|p| format!("{:.0}", p.angle))
            .collect();
        format!(
            "Laminate([{}], h={:.4}mm)",
            angles.join("/"),
            self.inner.total_thickness * 1000.0
        )
    }
}
