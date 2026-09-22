//! `PyMeshModel` — Python wrapper for the `aeroelast_mesh` mesh model.

use numpy::ndarray::{Array1, Array2};
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::prelude::*;

// ============================================================================
// PyMeshModel — wraps aeroelast_mesh::MeshModel
// ============================================================================

/// Python wrapper for `MeshModel` (aeroelast-mesh crate).
///
/// Holds the full mesh topology (nodes, elements, sets) loaded from HDF5.
/// Use `PyMeshModel.from_hdf5(path)` to load, then pass to
/// `PyMeshAssembler.from_mesh_model(mesh, properties)`.
#[pyclass(name = "MeshModel")]
pub(crate) struct PyMeshModel {
    pub(crate) inner: aeroelast_mesh::model::MeshModel,
}

#[pymethods]
impl PyMeshModel {
    /// Load a MeshModel from an HDF5 file written by the Python writer.
    #[staticmethod]
    fn from_hdf5(path: &str) -> PyResult<Self> {
        let inner = aeroelast_mesh::io::load_hdf5(path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
        Ok(PyMeshModel { inner })
    }

    /// Number of nodes in the mesh.
    #[getter]
    fn node_count(&self) -> usize {
        self.inner.node_count()
    }

    /// Number of elements in the mesh.
    #[getter]
    fn element_count(&self) -> usize {
        self.inner.element_count()
    }

    /// Node names of sets as a list of strings.
    fn node_set_names(&self) -> Vec<String> {
        self.inner.node_sets.keys().cloned().collect()
    }

    /// Element set names as a list of strings.
    fn element_set_names(&self) -> Vec<String> {
        self.inner.element_sets.keys().cloned().collect()
    }

    /// Node coordinates as a flat (n_nodes × 3) numpy array (row-major).
    fn node_coords_flat<'py>(&mut self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        let n = self.inner.node_count();
        let flat = self.inner.node_coords_flat().to_vec();
        Array2::from_shape_vec((n, 3), flat)
            .expect("node_coords_flat shape mismatch")
            .into_pyarray(py)
    }

    /// Node IDs in index order (matches row order of node_coords_flat).
    fn node_ids<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        let ids: Vec<u64> = self.inner.nodes.iter().map(|n| n.id).collect();
        Array1::from(ids).into_pyarray(py)
    }

    /// Element IDs in index order.
    fn element_ids<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u64>> {
        let ids: Vec<u64> = self.inner.elements.iter().map(|e| e.id).collect();
        Array1::from(ids).into_pyarray(py)
    }

    /// Returns a list of node-id lists, one per element, in element order.
    fn element_node_ids(&self) -> Vec<Vec<u64>> {
        self.inner.elements.iter().map(|e| e.node_ids.clone()).collect()
    }

    /// Node IDs for elements in a named element set. Returns flat list.
    fn element_set_node_ids(&self, set_name: &str) -> PyResult<Vec<u64>> {
        let set = self.inner.element_sets.get(set_name).ok_or_else(|| {
            pyo3::exceptions::PyKeyError::new_err(format!("element set '{set_name}' not found"))
        })?;
        Ok(set.element_ids.clone())
    }

    /// Node IDs for a named node set.
    fn node_set_node_ids(&self, set_name: &str) -> PyResult<Vec<u64>> {
        let set = self.inner.node_sets.get(set_name).ok_or_else(|| {
            pyo3::exceptions::PyKeyError::new_err(format!("node set '{set_name}' not found"))
        })?;
        Ok(set.node_ids.clone())
    }

    fn __repr__(&self) -> String {
        format!(
            "MeshModel(nodes={}, elements={}, node_sets={}, element_sets={})",
            self.inner.node_count(),
            self.inner.element_count(),
            self.inner.node_sets.len(),
            self.inner.element_sets.len(),
        )
    }

    /// Construct a MeshModel from raw Python data.
    ///
    /// Parameters
    /// ----------
    /// node_ids : list[int]
    ///     Node IDs in index order.
    /// node_coords_flat : list[float]
    ///     Flat (n_nodes × 3) coordinate array — row-major.
    /// element_ids : list[int]
    ///     Element IDs in index order.
    /// element_node_ids : list[list[int]]
    ///     Per-element node-ID lists (variable length).
    /// element_type_codes : list[int]
    ///     Assembler type code per element (3, 4, 33, 44, 208, …).
    /// element_sets : dict[str, list[int]]
    ///     Named element sets — values are lists of element IDs.
    /// node_sets : dict[str, list[int]]
    ///     Named node sets — values are lists of node IDs.
    #[staticmethod]
    #[pyo3(signature = (node_ids, node_coords_flat, element_ids, element_node_ids, element_type_codes, element_sets, node_sets))]
    fn from_raw_data(
        node_ids: Vec<u64>,
        node_coords_flat: Vec<f64>,
        element_ids: Vec<u64>,
        element_node_ids: Vec<Vec<u64>>,
        element_type_codes: Vec<i32>,
        element_sets: std::collections::HashMap<String, Vec<u64>>,
        node_sets: std::collections::HashMap<String, Vec<u64>>,
    ) -> PyResult<Self> {
        use aeroelast_mesh::entities::{Element, ElementSet, ElementType, Node, NodeSet};
        use aeroelast_mesh::model::MeshModel;

        let n_nodes = node_ids.len();
        if node_coords_flat.len() != n_nodes * 3 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "node_coords_flat length {} != node_ids.len() * 3 = {}",
                node_coords_flat.len(), n_nodes * 3
            )));
        }

        let mut model = MeshModel::new();

        // Build nodes
        for (i, &nid) in node_ids.iter().enumerate() {
            let x = node_coords_flat[i * 3];
            let y = node_coords_flat[i * 3 + 1];
            let z = node_coords_flat[i * 3 + 2];
            model.add_node(Node::new(nid, x, y, z))
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        }

        // Build elements
        let elem_type_from_code = |code: i32, n_nodes_e: usize| -> ElementType {
            match code {
                3  => ElementType::Triangle3,
                33 => ElementType::CompTri3,
                4  => ElementType::Quad4,
                44 => ElementType::CompQuad4,
                104 => ElementType::Quad4,
                108 => ElementType::Quad8,
                109 => ElementType::Quad9,
                _ => if n_nodes_e == 3 { ElementType::Triangle3 } else { ElementType::Quad4 },
            }
        };

        for (i, &eid) in element_ids.iter().enumerate() {
            let node_ids_elem = element_node_ids[i].clone();
            let code = element_type_codes.get(i).copied().unwrap_or(0);
            let etype = elem_type_from_code(code, node_ids_elem.len());
            model.add_element(Element::new(eid, node_ids_elem, etype))
                .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        }

        // Build element sets
        for (set_name, eids) in element_sets {
            model.element_sets.insert(set_name.clone(), ElementSet::with_ids(set_name, eids));
        }

        // Build node sets
        for (set_name, nids) in node_sets {
            model.node_sets.insert(set_name.clone(), NodeSet::with_ids(set_name, nids));
        }

        Ok(PyMeshModel { inner: model })
    }
}
