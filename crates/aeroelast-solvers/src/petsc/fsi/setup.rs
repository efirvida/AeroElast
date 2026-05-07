/// FSI setup utilities: BC reduction, mass lumping, and interface DOF remapping.
///
/// These routines bridge the gap between the global assembled COO matrices
/// (produced by Python via the Rust assembler) and the reduced system that
/// [`crate::petsc::elasticity::dynamic_newmark::NewmarkStepper`] expects.
///
/// # Workflow
/// 1. [`reduce_coo`]         — extract the free-DOF sub-system from global COO.
/// 2. [`lump_mass_coo`]      — convert consistent mass to lumped diagonal form.
/// 3. [`rayleigh_c_coo`]     — build Rayleigh C = η_k·K + η_m·M from free-DOF COO.
/// 4. [`remap_interface_dofs`] — translate global DOF indices to reduced indices.
/// 5. [`extract_masses`]     — extract nodal masses from lumped mass diagonal.
/// 6. [`rotate_mesh_coords`] — rotate mesh coordinates using Rodrigues formula.
/// 7. [`reassemble_k`]       — reassemble K after coordinate update (inertial solver).

use std::collections::HashMap;

use aeroelast_core::assembly::assembler::MeshAssembler;
use crate::petsc::fsi::rotor_physics::RotorTransforms;

// ── BC reduction ──────────────────────────────────────────────────────────────

/// Extract the free-DOF sub-system from a global COO matrix.
///
/// Given global COO triplets `(rows, cols, vals)` and a sorted list of free
/// DOF indices `free_dofs`, returns new COO triplets whose row/column indices
/// run from `0` to `free_dofs.len() - 1`.
///
/// Entries whose row **or** column refers to a constrained (non-free) DOF are
/// discarded — this is the standard condensation for homogeneous Dirichlet BCs.
///
/// # Arguments
/// * `rows`      — global row indices (i32)
/// * `cols`      — global column indices (i32)
/// * `vals`      — matrix values
/// * `free_dofs` — sorted array of free (unconstrained) global DOF indices
///
/// # Returns
/// `(rows_red, cols_red, vals_red)` — reduced COO triplets.
pub fn reduce_coo(
    rows: &[i32],
    cols: &[i32],
    vals: &[f64],
    free_dofs: &[i32],
) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    // Build global_dof → reduced_index lookup.
    let dof_map: HashMap<i32, i32> = free_dofs
        .iter()
        .enumerate()
        .map(|(i, &d)| (d, i as i32))
        .collect();

    let mut r_out = Vec::with_capacity(vals.len());
    let mut c_out = Vec::with_capacity(vals.len());
    let mut v_out = Vec::with_capacity(vals.len());

    for ((&row, &col), &val) in rows.iter().zip(cols.iter()).zip(vals.iter()) {
        if let (Some(&ri), Some(&ci)) = (dof_map.get(&row), dof_map.get(&col)) {
            r_out.push(ri);
            c_out.push(ci);
            v_out.push(val);
        }
    }

    (r_out, c_out, v_out)
}

/// Precompute a per-entry mapping from the **full** (element-assembled) K_G COO to
/// output indices in the reduced reference sparsity pattern.
///
/// Call this ONCE in the solver constructor (e.g. with a zero-stress dummy assembly)
/// to build a `Vec<i32>` of the same length as the full COO.  Each element is either:
/// * `>= 0`  — the index into the reduced output vector where this value must be summed, or
/// * `-1`    — the entry is outside the free-DOF set and must be skipped.
///
/// At runtime use [`apply_kg_coo_map`] to accumulate values with no heap allocation.
///
/// # Arguments
/// * `full_rows`, `full_cols` — raw element-assembled COO row/col indices (i64).
/// * `free_dofs` — sorted list of free global DOF indices (i32).
/// * `ref_rows`, `ref_cols` — the REDUCED K sparsity (unique, sorted, 0-based).
///
/// # Assumptions
/// `ref_rows` is non-decreasing and, within each row, `ref_cols` is sorted.
/// This is always the case when the reference COO comes from PETSc CSR → COO.
pub fn build_kg_coo_map(
    full_rows: &[i64],
    full_cols: &[i64],
    free_dofs: &[i32],
    ref_rows: &[i32],
    ref_cols: &[i32],
) -> Vec<i32> {
    // global_dof (i64) → reduced index (i32)
    let dof_map: HashMap<i64, i32> = free_dofs
        .iter()
        .enumerate()
        .map(|(i, &d)| (d as i64, i as i32))
        .collect();

    full_rows
        .iter()
        .zip(full_cols.iter())
        .map(|(&r, &c)| {
            let (Some(&ri), Some(&ci)) = (dof_map.get(&r), dof_map.get(&c)) else {
                return -1i32;
            };
            // ref_rows is sorted non-decreasing; find the row block for ri.
            let lo = ref_rows.partition_point(|&x| x < ri);
            let hi = ref_rows.partition_point(|&x| x <= ri);
            // Within the row block ref_cols is sorted; binary search for ci.
            match ref_cols[lo..hi].binary_search(&ci) {
                Ok(pos) => (lo + pos) as i32,
                Err(_) => -1,
            }
        })
        .collect()
}

/// Apply a precomputed K_G COO map (from [`build_kg_coo_map`]) to a value slice.
///
/// Accumulates `full_vals[i]` into `out[map[i]]` for every entry where `map[i] >= 0`.
/// This is an O(N) scan with no heap allocation — suitable for hot per-timestep paths.
///
/// # Arguments
/// * `kg_coo_map` — precomputed index map, same length as `full_vals`.
/// * `full_vals`  — raw element-assembled K_G values.
/// * `n_ref`      — length of the reduced output vector (= `ref_rows.len()`).
///
/// # Returns
/// A `Vec<f64>` of length `n_ref` with accumulated values in reduced sparsity order.
pub fn apply_kg_coo_map(kg_coo_map: &[i32], full_vals: &[f64], n_ref: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n_ref];
    for (&idx, &val) in kg_coo_map.iter().zip(full_vals.iter()) {
        if idx >= 0 {
            out[idx as usize] += val;
        }
    }
    out
}

// ── Mass lumping ──────────────────────────────────────────────────────────────

/// Convert a consistent mass matrix (COO) to a diagonal (lumped) form.
///
/// Uses the **row-sum** lumping method: the lumped diagonal coefficient for
/// DOF `i` is the sum of the absolute values of all entries in row `i`.
///
/// A mass floor of `1e-4 × global_max` is applied after lumping. This is
/// required to handle rotational DOFs and tip-closure elements whose thickness
/// → 0: without it those DOFs get M ≈ 0, which makes K_eff = K + a0·0 = K
/// for those rows and creates extreme ill-conditioning in the Newmark effective
/// stiffness (ratio of K_eff diagonal entries between translational and
/// rotational DOFs can reach 1e7+). The Python implementation (`lump_mass_matrix`)
/// applied the same floor before this logic was migrated to Rust.
///
/// Returns diagonal COO triplets: `rows = cols = [0, 1, …, n_dofs-1]`.
///
/// # Arguments
/// * `rows`   — COO row indices (in the **reduced** system, 0-based)
/// * `cols`   — COO column indices
/// * `vals`   — matrix values
/// * `n_dofs` — number of rows/columns in the reduced system
///
/// # Returns
/// `(rows_diag, cols_diag, vals_diag)` — diagonal COO.
pub fn lump_mass_coo(
    rows: &[i32],
    cols: &[i32],
    vals: &[f64],
    n_dofs: usize,
) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    let _ = cols; // not needed for row-sum lumping
    let mut diag = vec![0.0f64; n_dofs];

    for (&r, &v) in rows.iter().zip(vals.iter()) {
        let idx = r as usize;
        if idx < n_dofs {
            diag[idx] += v.abs();
        }
    }

    // Apply the same mass floor that the original Python lump_mass_matrix used.
    // Prevents zero-mass DOFs (rotational DOFs of shell elements, tip-closure
    // elements) from making K_eff singular/ill-conditioned in the Newmark solve.
    let m_max = diag.iter().cloned().fold(0.0f64, f64::max);
    if m_max > 0.0 {
        let floor = m_max * 1e-4;
        for v in &mut diag {
            if *v < floor {
                *v = floor;
            }
        }
    }

    let indices: Vec<i32> = (0..n_dofs as i32).collect();
    (indices.clone(), indices, diag)
}

// ── Rayleigh damping ──────────────────────────────────────────────────────────

/// Build Rayleigh damping C = η_k·K + η_m·M from reduced COO triplets.
///
/// K and M must share the same sparsity pattern (same `rows` and `cols`
/// arrays), which is always true for matrices assembled from the same mesh
/// topology under Rayleigh damping.
///
/// Returns COO triplets for C with the same sparsity pattern.
///
/// # Arguments
/// * `k_vals` — stiffness matrix values (reduced)
/// * `m_vals` — mass matrix values (reduced, same sparsity as K)
/// * `rows`   — shared row indices
/// * `cols`   — shared column indices
/// * `eta_k`  — stiffness-proportional coefficient (α)
/// * `eta_m`  — mass-proportional coefficient (β)
///
/// # Returns
/// `(rows, cols, c_vals)` — the `rows` and `cols` are cloned from input.
pub fn rayleigh_c_coo(
    k_vals: &[f64],
    m_vals: &[f64],
    rows: &[i32],
    cols: &[i32],
    eta_k: f64,
    eta_m: f64,
) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    assert_eq!(
        k_vals.len(),
        m_vals.len(),
        "K and M must have the same number of non-zeros for Rayleigh damping"
    );
    let c_vals: Vec<f64> = k_vals
        .iter()
        .zip(m_vals.iter())
        .map(|(&k, &m)| eta_k * k + eta_m * m)
        .collect();
    (rows.to_vec(), cols.to_vec(), c_vals)
}

// ── Interface DOF remapping ───────────────────────────────────────────────────

/// Translate global DOF indices to their positions in the reduced system.
///
/// Given a sorted `free_dofs` array (the output of BC reduction), each global
/// DOF `g` in `interface_dofs_global` is mapped to its **reduced** index via
/// binary search.
///
/// # Constrained interface DOFs
/// If a DOF is not found in `free_dofs` (i.e. the interface node shares a
/// Dirichlet BC — e.g. root nodes on an aerodynamic coupling surface), the
/// entry is mapped to `usize::MAX`.  All callers guard reads/writes with a
/// bounds check, so:
/// - Force scatter: `usize::MAX < n_dofs` → force silently dropped (goes to
///   reaction at the support, which is physically correct).
/// - Displacement gather: `usize::MAX < u.len()` → returns 0.0 (prescribed
///   displacement, also correct).
///
/// # Arguments
/// * `interface_dofs_global` — DOF indices in the global (unreduced) system
/// * `free_dofs`             — sorted array of free global DOF indices
///
/// # Returns
/// DOF indices in the reduced (0-based) system, or `usize::MAX` for fixed DOFs.
pub fn remap_interface_dofs(
    interface_dofs_global: &[usize],
    free_dofs: &[i32],
) -> Vec<usize> {
    interface_dofs_global
        .iter()
        .map(|&d| {
            let d_i32 = d as i32;
            let pos = free_dofs.partition_point(|&f| f < d_i32);
            if pos < free_dofs.len() && free_dofs[pos] == d_i32 {
                pos
            } else {
                usize::MAX
            }
        })
        .collect()
}

// ── Sparsity alignment ────────────────────────────────────────────────────────

/// Expand a lumped-mass diagonal into the sparsity pattern of matrix K.
///
/// [`NewmarkStepper`] requires K, M, and C to share the same COO sparsity
/// pattern (same `rows` and `cols` arrays, aligned entry-by-entry).  When
/// the mass matrix has been row-sum lumped to a diagonal, this function
/// "pads" it back to K's full sparsity by inserting zeros at off-diagonal
/// positions.
///
/// The resulting `m_vals` slice is aligned with `k_rows/k_cols` so that it
/// can be passed directly to `NewmarkStepper::new` alongside K's triplets.
///
/// # Arguments
/// * `lumped_diag` — diagonal mass values, length = number of free DOFs
/// * `k_rows`      — COO row indices for K (reduced system)
/// * `k_cols`      — COO column indices for K (reduced system)
///
/// # Returns
/// `m_vals` aligned with `k_rows/k_cols`.
pub fn expand_diag_to_sparsity(lumped_diag: &[f64], k_rows: &[i32], k_cols: &[i32]) -> Vec<f64> {
    k_rows
        .iter()
        .zip(k_cols.iter())
        .map(|(&r, &c)| {
            if r == c {
                let idx = r as usize;
                if idx < lumped_diag.len() {
                    lumped_diag[idx]
                } else {
                    0.0
                }
            } else {
                0.0
            }
        })
        .collect()
}

// ── Inertial solver helpers ───────────────────────────────────────────────────

/// Extract nodal masses from a lumped mass diagonal.
///
/// The inertial rotor solver needs nodal masses (one scalar per node) to compute
/// the reference load vector `F_ref = −M·a_ref`. This function extracts the first
/// translational DOF mass for each node from the lumped diagonal.
///
/// # Arguments
/// * `lumped_diag` — lumped mass diagonal (length = n_nodes × dofs_per_node)
/// * `dofs_per_node` — DOF stride (typically 6 for shells, 3 for solids)
///
/// # Returns
/// `Vec<f64>` — nodal masses (length = n_nodes), one scalar per node.
///
/// # Panics
/// Panics if `lumped_diag.len()` is not divisible by `dofs_per_node`.
pub fn extract_masses(lumped_diag: &[f64], dofs_per_node: usize) -> Vec<f64> {
    assert_eq!(
        lumped_diag.len() % dofs_per_node,
        0,
        "lumped_diag length must be a multiple of dofs_per_node"
    );
    let n_nodes = lumped_diag.len() / dofs_per_node;
    let mut masses = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        // Extract the first translational DOF mass (DOF 0 of the node).
        // For lumped mass, all translational DOFs (0, 1, 2) have the same value.
        masses.push(lumped_diag[i * dofs_per_node]);
    }
    masses
}

/// Rotate mesh coordinates using Rodrigues' formula.
///
/// Given reference (unrotated) coordinates and a rotation angle θ about the
/// axis defined in `transforms`, returns the rotated coordinates.
///
/// # Arguments
/// * `coords_ref` — reference node coordinates (flat, length = n_nodes × 3)
/// * `transforms` — `RotorTransforms` with rotation axis and center
/// * `theta` — rotation angle [rad]
///
/// # Returns
/// `Vec<f64>` — rotated coordinates (flat, length = n_nodes × 3)
pub fn rotate_mesh_coords(
    coords_ref: &[f64],
    transforms: &RotorTransforms,
    theta: f64,
) -> Vec<f64> {
    let n_nodes = coords_ref.len() / 3;
    let mut coords_rotated = Vec::with_capacity(coords_ref.len());

    let r_mat = transforms.rotation_matrix(theta);

    for i in 0..n_nodes {
        let x0 = coords_ref[3 * i];
        let y0 = coords_ref[3 * i + 1];
        let z0 = coords_ref[3 * i + 2];

        // Translate to rotation center
        let dx = x0 - transforms.center[0];
        let dy = y0 - transforms.center[1];
        let dz = z0 - transforms.center[2];

        // Apply rotation matrix
        let rx = r_mat[0][0] * dx + r_mat[0][1] * dy + r_mat[0][2] * dz;
        let ry = r_mat[1][0] * dx + r_mat[1][1] * dy + r_mat[1][2] * dz;
        let rz = r_mat[2][0] * dx + r_mat[2][1] * dy + r_mat[2][2] * dz;

        // Translate back
        coords_rotated.push(rx + transforms.center[0]);
        coords_rotated.push(ry + transforms.center[1]);
        coords_rotated.push(rz + transforms.center[2]);
    }

    coords_rotated
}

/// Reassemble the elastic stiffness matrix K after updating mesh coordinates.
///
/// Used by the **inertial rotor solver** to rebuild K(θ) after the mesh has
/// been rotated to a new angle θ. The function:
/// 1. Calls `assembler.update_node_coordinates(coords_rotated)` to rebuild
///    element precomputed data (Jacobians, normals, etc.).
/// 2. Calls `assembler.assemble_k()` to get the new K in COO format.
///
/// The caller is responsible for reducing the COO to the free-DOF system
/// (via `reduce_coo`) before passing to `NewmarkStepper::update_elastic_stiffness`.
///
/// # Arguments
/// * `assembler` — mutable reference to the `MeshAssembler`
/// * `coords_rotated` — new (rotated) node coordinates (flat, length = n_nodes × 3)
///
/// # Returns
/// `(rows, cols, vals)` — COO triplets for the **global** (unreduced) K matrix.
///
/// # Notes
/// - This is an O(n_elems) operation (recomputes all element Jacobians).
/// - The co-rotational solver does **not** use this (its K is constant).
/// - For the inertial solver, call this every time the mesh rotates significantly.
pub fn reassemble_k(
    assembler: &mut MeshAssembler,
    coords_rotated: &[f64],
) -> (Vec<i64>, Vec<i64>, Vec<f64>) {
    // Step 1: Update mesh coordinates (rebuilds element precomputed data)
    assembler.update_node_coordinates(coords_rotated);

    // Step 2: Reassemble K
    assembler.assemble_k()
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reduce_coo_2dof_system() {
        // 4×4 system, DOFs 0-3; constrain DOF 1 → free = [0, 2, 3]
        // Matrix: identity
        let rows = vec![0i32, 1, 2, 3];
        let cols = vec![0i32, 1, 2, 3];
        let vals = vec![1.0, 2.0, 3.0, 4.0];
        let free_dofs = vec![0i32, 2, 3];

        let (r, c, v) = reduce_coo(&rows, &cols, &vals, &free_dofs);

        // DOF 1 is constrained → only entries for DOFs 0, 2, 3 survive.
        // Remapped: 0→0, 2→1, 3→2
        assert_eq!(r, vec![0i32, 1, 2]);
        assert_eq!(c, vec![0i32, 1, 2]);
        assert_eq!(v, vec![1.0, 3.0, 4.0]);
    }

    #[test]
    fn test_reduce_coo_off_diagonal() {
        // 3×3 system, free = [0, 2]; off-diagonal (0,2) should survive remapped to (0,1)
        let rows = vec![0i32, 0, 1, 2, 2];
        let cols = vec![0i32, 2, 1, 0, 2];
        let vals = vec![10.0, 5.0, 99.0, 5.0, 20.0];
        let free_dofs = vec![0i32, 2];

        let (r, c, v) = reduce_coo(&rows, &cols, &vals, &free_dofs);

        // Entry (1,1)=99 is dropped (DOF 1 is constrained)
        assert_eq!(r.len(), 4);
        assert!(r.contains(&0) && r.contains(&1));
        let _ = (c, v); // just check lengths
    }

    #[test]
    fn test_lump_mass_coo_diagonal() {
        // 2×2 diagonal mass: [[3, 0], [0, 5]] → lumped = [3, 5]
        let rows = vec![0i32, 1];
        let cols = vec![0i32, 1];
        let vals = vec![3.0, 5.0];
        let (_, _, v) = lump_mass_coo(&rows, &cols, &vals, 2);
        assert_eq!(v, vec![3.0, 5.0]);
    }

    #[test]
    fn test_lump_mass_coo_consistent() {
        // 2×2 consistent: [[2, 1], [1, 2]] → row sums = [3, 3]
        let rows = vec![0i32, 0, 1, 1];
        let cols = vec![0i32, 1, 0, 1];
        let vals = vec![2.0, 1.0, 1.0, 2.0];
        let (_, _, v) = lump_mass_coo(&rows, &cols, &vals, 2);
        assert!((v[0] - 3.0).abs() < 1e-12);
        assert!((v[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_c_coo() {
        let k_vals = vec![1.0, 2.0];
        let m_vals = vec![3.0, 4.0];
        let rows = vec![0i32, 1];
        let cols = vec![0i32, 1];
        let (_, _, c) = rayleigh_c_coo(&k_vals, &m_vals, &rows, &cols, 0.1, 0.2);
        // c[0] = 0.1*1 + 0.2*3 = 0.7,  c[1] = 0.1*2 + 0.2*4 = 1.0
        assert!((c[0] - 0.7).abs() < 1e-12);
        assert!((c[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_remap_interface_dofs() {
        // free_dofs = [2, 5, 7, 10]; interface global DOFs: [5, 10]
        // Expected: 5→1, 10→3
        let free_dofs = vec![2i32, 5, 7, 10];
        let interface = vec![5usize, 10];
        let remapped = remap_interface_dofs(&interface, &free_dofs);
        assert_eq!(remapped, vec![1, 3]);
    }

    #[test]
    fn test_remap_interface_dofs_constrained() {
        // DOF 3 is NOT in free_dofs → should map to usize::MAX (not panic)
        let free_dofs = vec![2i32, 5, 7, 10];
        let interface = vec![2usize, 3, 7];
        let remapped = remap_interface_dofs(&interface, &free_dofs);
        assert_eq!(remapped[0], 0);           // 2 → index 0
        assert_eq!(remapped[1], usize::MAX);  // 3 → constrained
        assert_eq!(remapped[2], 2);           // 7 → index 2
    }

    #[test]
    fn test_expand_diag_to_sparsity() {
        // 3×3 K: diag + one off-diagonal. lumped_diag = [10, 20, 30]
        // K triplets: (0,0), (0,1), (1,1), (2,2)
        let k_rows = vec![0i32, 0, 1, 2];
        let k_cols = vec![0i32, 1, 1, 2];
        let diag = vec![10.0, 20.0, 30.0];
        let m = expand_diag_to_sparsity(&diag, &k_rows, &k_cols);
        // (0,0)→10, (0,1)→0, (1,1)→20, (2,2)→30
        assert_eq!(m, vec![10.0, 0.0, 20.0, 30.0]);
    }

    // ── Inertial solver helper tests ─────────────────────────────────────────

    #[test]
    fn test_extract_masses_6dof() {
        // 2 nodes × 6 DOFs = 12 entries. Lumped mass: [m, m, m, 0, 0, 0] per node.
        let lumped_diag = vec![
            5.0, 5.0, 5.0, 0.0, 0.0, 0.0,  // node 0: m=5
            3.0, 3.0, 3.0, 0.0, 0.0, 0.0,  // node 1: m=3
        ];
        let masses = extract_masses(&lumped_diag, 6);
        assert_eq!(masses.len(), 2);
        assert_eq!(masses[0], 5.0);
        assert_eq!(masses[1], 3.0);
    }

    #[test]
    fn test_extract_masses_3dof() {
        // 3 nodes × 3 DOFs = 9 entries. Lumped mass: [m, m, m] per node.
        let lumped_diag = vec![
            2.0, 2.0, 2.0,  // node 0: m=2
            4.0, 4.0, 4.0,  // node 1: m=4
            1.0, 1.0, 1.0,  // node 2: m=1
        ];
        let masses = extract_masses(&lumped_diag, 3);
        assert_eq!(masses, vec![2.0, 4.0, 1.0]);
    }

    #[test]
    fn test_rotate_mesh_coords_90deg_z_axis() {
        // Single node at [1, 0, 0], rotate 90° about Z axis → [0, 1, 0]
        use std::f64::consts::FRAC_PI_2;
        let coords_ref = vec![1.0, 0.0, 0.0];
        let transforms = RotorTransforms::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let coords_rotated = rotate_mesh_coords(&coords_ref, &transforms, FRAC_PI_2);

        assert_eq!(coords_rotated.len(), 3);
        assert!(coords_rotated[0].abs() < 1e-14, "x={}", coords_rotated[0]);
        assert!((coords_rotated[1] - 1.0).abs() < 1e-14, "y={}", coords_rotated[1]);
        assert!(coords_rotated[2].abs() < 1e-14, "z={}", coords_rotated[2]);
    }

    #[test]
    fn test_rotate_mesh_coords_identity() {
        // Rotate by θ=0 → identity
        let coords_ref = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let transforms = RotorTransforms::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let coords_rotated = rotate_mesh_coords(&coords_ref, &transforms, 0.0);

        assert_eq!(coords_rotated.len(), coords_ref.len());
        for (a, b) in coords_ref.iter().zip(coords_rotated.iter()) {
            assert!((a - b).abs() < 1e-14, "a={a}, b={b}");
        }
    }

    #[test]
    fn test_rotate_mesh_coords_with_center() {
        // Node at [2, 0, 0], rotate 90° about Z with center=[1, 0, 0] → [1, 1, 0]
        use std::f64::consts::FRAC_PI_2;
        let coords_ref = vec![2.0, 0.0, 0.0];
        let transforms = RotorTransforms::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let coords_rotated = rotate_mesh_coords(&coords_ref, &transforms, FRAC_PI_2);

        // r = [2,0,0] - [1,0,0] = [1,0,0]
        // R(90°)·[1,0,0] = [0,1,0]
        // r' + center = [0,1,0] + [1,0,0] = [1,1,0]
        assert!((coords_rotated[0] - 1.0).abs() < 1e-14, "x={}", coords_rotated[0]);
        assert!((coords_rotated[1] - 1.0).abs() < 1e-14, "y={}", coords_rotated[1]);
        assert!(coords_rotated[2].abs() < 1e-14, "z={}", coords_rotated[2]);
    }
}
