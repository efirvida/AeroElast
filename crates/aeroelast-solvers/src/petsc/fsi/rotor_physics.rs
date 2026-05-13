/// Rotor physics utilities for FSI solvers (co-rotational and inertial).
///
/// Provides:
/// - [`RotorTransforms`]  — Rodrigues rotation, force/displacement coordinate transforms.
/// - Rotating-frame inertial forces — centrifugal, Coriolis, Euler (co-rotational solver).
/// - Inertial-frame physics — [`compute_rigid_body_acceleration_inertial`],
///   [`compute_rigid_body_velocity_inertial`], [`compute_reference_load_vector`]
///   (inertial solver).
/// - [`build_ksp_vals`]   — spin-softening K_SP values aligned to the K sparsity pattern.
/// - [`compute_torque`]   — scalar torque about the rotation axis from nodal forces.
/// - [`compute_performance_coefficients`] — Ct, Cp, Cq, TSR.
/// - [`OmegaProvider`]    — angular-velocity provider (constant, ramped, computed, mixed).
///
/// All coordinate arrays use the flat layout `[x0, y0, z0, x1, y1, z1, …]` with stride 3.
/// Force arrays returned by the inertial-force functions follow the same layout (n_nodes × 3).
///
/// # Feature gate
/// Compiled only with `--features fsi`.

// ── Skew-symmetric helper ──────────────────────────────────────────────────────

/// Build the 3×3 skew-symmetric (cross-product) matrix K such that `K·x = axis × x`.
fn skew(axis: &[f64; 3]) -> [[f64; 3]; 3] {
    let [nx, ny, nz] = *axis;
    [
        [0.0, -nz, ny],
        [nz, 0.0, -nx],
        [-ny, nx, 0.0],
    ]
}

/// 3×3 matrix product A·B.
fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

// ── RotorTransforms ────────────────────────────────────────────────────────────

/// Coordinate transformation utilities for a rotating reference frame.
///
/// Uses Rodrigues' formula to compute rotation matrices for an arbitrary
/// rotation axis:
///
/// ```text
/// R(θ) = I + sin(θ)·K + (1 − cos(θ))·K²
/// ```
///
/// where K is the skew-symmetric matrix of the rotation axis.
///
/// # Conventions
/// * `to_rotating`   — `v_local  = R^T · v_global`
/// * `to_inertial`   — `v_global = R   · v_local`
#[derive(Debug, Clone)]
pub struct RotorTransforms {
    /// Normalized rotation axis unit vector.
    pub axis: [f64; 3],
    /// Rotation center (default: origin).
    pub center: [f64; 3],
    /// Precomputed skew-symmetric matrix K of the axis.
    k: [[f64; 3]; 3],
    /// Precomputed K² = K·K.
    k2: [[f64; 3]; 3],
}

impl RotorTransforms {
    /// Create a new `RotorTransforms` with the given rotation axis.
    ///
    /// # Panics
    /// Panics if `axis` has zero norm (cannot normalize).
    pub fn new(axis: [f64; 3], center: [f64; 3]) -> Self {
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        assert!(norm > 1e-12, "rotation axis must be non-zero");
        let n = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        let k = skew(&n);
        let k2 = mat3_mul(&k, &k);
        Self { axis: n, center, k, k2 }
    }

    /// Compute the 3×3 rotation matrix for angle `theta` (radians).
    ///
    /// `R(θ) = I + sin(θ)·K + (1 − cos(θ))·K²`
    pub fn rotation_matrix(&self, theta: f64) -> [[f64; 3]; 3] {
        let s = theta.sin();
        let c = 1.0 - theta.cos();
        let mut r = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let ident = if i == j { 1.0 } else { 0.0 };
                r[i][j] = ident + s * self.k[i][j] + c * self.k2[i][j];
            }
        }
        r
    }

    /// Transform a flat force array from the global (inertial) frame to the
    /// rotating (local) frame: `F_local = R^T · F_global`.
    ///
    /// `forces_global` — flat `[fx0, fy0, fz0, fx1, fy1, fz1, …]`.
    /// Returns a new flat array with the same layout.
    pub fn forces_to_rotating(&self, forces_global: &[f64], theta: f64) -> Vec<f64> {
        let r = self.rotation_matrix(theta);
        let n = forces_global.len() / 3;
        let mut out = vec![0.0f64; n * 3];
        for i in 0..n {
            let b = i * 3;
            let f = [forces_global[b], forces_global[b + 1], forces_global[b + 2]];
            // v_local = R^T · v_global  ↔  dot with rows of R
            for j in 0..3 {
                out[b + j] = r[0][j] * f[0] + r[1][j] * f[1] + r[2][j] * f[2];
            }
        }
        out
    }

    /// Transform a flat displacement array from the rotating (local) frame to
    /// the global (inertial) frame: `u_global = R · u_local`.
    ///
    /// `disps_local` — flat `[ux0, uy0, uz0, ux1, uy1, uz1, …]`.
    /// Returns a new flat array with the same layout.
    pub fn disps_to_inertial(&self, disps_local: &[f64], theta: f64) -> Vec<f64> {
        let r = self.rotation_matrix(theta);
        let n = disps_local.len() / 3;
        let mut out = vec![0.0f64; n * 3];
        for i in 0..n {
            let b = i * 3;
            let u = [disps_local[b], disps_local[b + 1], disps_local[b + 2]];
            for j in 0..3 {
                out[b + j] = r[j][0] * u[0] + r[j][1] * u[1] + r[j][2] * u[2];
            }
        }
        out
    }

    /// Transform gravity vector from the inertial frame to the rotating frame.
    ///
    /// `g_global` — gravity vector `[gx, gy, gz]` in the inertial frame.
    /// Returns `R^T · g_global`.
    pub fn gravity_to_rotating(&self, g_global: [f64; 3], theta: f64) -> [f64; 3] {
        let r = self.rotation_matrix(theta);
        let mut out = [0.0f64; 3];
        for j in 0..3 {
            out[j] = r[0][j] * g_global[0] + r[1][j] * g_global[1] + r[2][j] * g_global[2];
        }
        out
    }
}

// ── Inertial forces ────────────────────────────────────────────────────────────

/// Compute centrifugal forces for all nodes.
///
/// `F_cf,i = m_i · ω² · r_⊥,i`  where  `r_⊥ = r − (r·n̂)·n̂`
///
/// # Arguments
/// * `coords`  — flat node coordinates `[x0,y0,z0,…]` in the rotating frame
/// * `masses`  — per-node scalar masses, length `n_nodes`
/// * `axis`    — rotation axis unit vector
/// * `center`  — rotation center coordinates
/// * `omega`   — angular velocity (rad/s)
///
/// # Returns
/// Flat force array `[fx0,fy0,fz0,…]` length `n_nodes * 3`.
pub fn compute_centrifugal_force(
    coords: &[f64],
    masses: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
    omega: f64,
) -> Vec<f64> {
    let n = masses.len();
    let mut out = vec![0.0f64; n * 3];
    if omega == 0.0 {
        return out;
    }
    let omega_sq = omega * omega;
    for i in 0..n {
        let b = i * 3;
        // r = coords - center
        let rx = coords[b] - center[0];
        let ry = coords[b + 1] - center[1];
        let rz = coords[b + 2] - center[2];
        // r_dot_axis = r · n̂
        let rda = rx * axis[0] + ry * axis[1] + rz * axis[2];
        // r_perp = r - (r·n̂) * n̂
        let rpx = rx - rda * axis[0];
        let rpy = ry - rda * axis[1];
        let rpz = rz - rda * axis[2];
        let m = masses[i];
        out[b] = m * omega_sq * rpx;
        out[b + 1] = m * omega_sq * rpy;
        out[b + 2] = m * omega_sq * rpz;
    }
    out
}

/// Compute Coriolis forces for all nodes.
///
/// `F_cor,i = −2 · m_i · (ω_vec × v_i)`
///
/// # Arguments
/// * `velocities` — flat node velocities `[vx0,vy0,vz0,…]` in the rotating frame
/// * `masses`     — per-node scalar masses
/// * `axis`       — rotation axis unit vector
/// * `omega`      — angular velocity (rad/s)
///
/// # Returns
/// Flat force array, length `n_nodes * 3`.
pub fn compute_coriolis_force(
    velocities: &[f64],
    masses: &[f64],
    axis: &[f64; 3],
    omega: f64,
) -> Vec<f64> {
    let n = masses.len();
    let mut out = vec![0.0f64; n * 3];
    if omega == 0.0 {
        return out;
    }
    // ω_vec = ω · n̂
    let wx = omega * axis[0];
    let wy = omega * axis[1];
    let wz = omega * axis[2];
    for i in 0..n {
        let b = i * 3;
        let vx = velocities[b];
        let vy = velocities[b + 1];
        let vz = velocities[b + 2];
        // ω × v
        let cx = wy * vz - wz * vy;
        let cy = wz * vx - wx * vz;
        let cz = wx * vy - wy * vx;
        let m = masses[i];
        out[b] = -2.0 * m * cx;
        out[b + 1] = -2.0 * m * cy;
        out[b + 2] = -2.0 * m * cz;
    }
    out
}

/// Compute Euler forces for all nodes (due to angular acceleration α = dω/dt).
///
/// `F_euler,i = −m_i · (α_vec × r_i)`
///
/// Evaluate at **deformed** coordinates `coords = X₀ + u`.
///
/// # Arguments
/// * `coords` — flat node coordinates in the rotating frame (may be deformed)
/// * `masses` — per-node scalar masses
/// * `axis`   — rotation axis unit vector
/// * `center` — rotation center
/// * `alpha`  — angular acceleration (rad/s²)
///
/// # Returns
/// Flat force array, length `n_nodes * 3`.
pub fn compute_euler_force(
    coords: &[f64],
    masses: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
    alpha: f64,
) -> Vec<f64> {
    let n = masses.len();
    let mut out = vec![0.0f64; n * 3];
    if alpha == 0.0 {
        return out;
    }
    // α_vec = α · n̂
    let ax = alpha * axis[0];
    let ay = alpha * axis[1];
    let az = alpha * axis[2];
    for i in 0..n {
        let b = i * 3;
        let rx = coords[b] - center[0];
        let ry = coords[b + 1] - center[1];
        let rz = coords[b + 2] - center[2];
        // α_vec × r
        let cx = ay * rz - az * ry;
        let cy = az * rx - ax * rz;
        let cz = ax * ry - ay * rx;
        let m = masses[i];
        out[b] = -m * cx;
        out[b + 1] = -m * cy;
        out[b + 2] = -m * cz;
    }
    out
}

/// Compute gravity forces for all nodes in the rotating frame.
///
/// `F_g,i = m_i · g_rot`
///
/// # Arguments
/// * `masses` — per-node scalar masses
/// * `g_rot`  — gravity vector already transformed to the rotating frame
///
/// # Returns
/// Flat force array, length `n_nodes * 3`.
pub fn compute_gravity_force(masses: &[f64], g_rot: &[f64; 3]) -> Vec<f64> {
    let n = masses.len();
    let mut out = vec![0.0f64; n * 3];
    for i in 0..n {
        let b = i * 3;
        let m = masses[i];
        out[b] = m * g_rot[0];
        out[b + 1] = m * g_rot[1];
        out[b + 2] = m * g_rot[2];
    }
    out
}

// ── Inertial-frame rigid-body acceleration ────────────────────────────────────

/// Compute rigid-body reference acceleration in the INERTIAL (global) frame.
///
/// For a node at position `r_i = X_rot_i − center` (where `X_rot_i` is the
/// node coordinate already rotated by R(θ)), the inertial acceleration is:
///
///   `a_ref_i = α × r_i + ω × (ω × r_i)`
///
/// This is the reference acceleration field used by the inertial rotor solver
/// to compute the load vector `F_ref = −M·a_ref` that accounts for the rigid-body
/// inertial effects without introducing fictitious forces.
///
/// # Arguments
/// * `coords_rotated` — flat node coordinates `[x0,y0,z0,…]` already rotated
///   to the current angle θ (i.e., `R(θ)·X₀`), in the global frame.
/// * `axis` — rotation axis unit vector `n̂`
/// * `center` — rotation center coordinates
/// * `omega` — angular velocity ω [rad/s]
/// * `alpha` — angular acceleration α [rad/s²]
///
/// # Returns
/// Flat acceleration array `[ax0,ay0,az0,…]` in the global frame, length `n_nodes * 3`.
///
/// # Notes
/// The two terms:
/// - `α × r`: tangential acceleration (Euler term)
/// - `ω × (ω × r)`: centripetal acceleration pointing toward the axis
///
/// This is the ONLY inertial reference load needed in the inertial formulation;
/// centrifugal, Coriolis, and Euler forces are NOT applied separately (they
/// would be fictitious forces in the rotating frame).
pub fn compute_rigid_body_acceleration_inertial(
    coords_rotated: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
    omega: f64,
    alpha: f64,
) -> Vec<f64> {
    let n = coords_rotated.len() / 3;
    let mut out = vec![0.0f64; n * 3];

    // Early return if both ω and α are zero
    if omega.abs() < 1e-14 && alpha.abs() < 1e-14 {
        return out;
    }

    let omega_sq = omega * omega;

    // α_vec = α · n̂
    let ax = alpha * axis[0];
    let ay = alpha * axis[1];
    let az = alpha * axis[2];

    for i in 0..n {
        let b = i * 3;
        // r = coords_rotated - center
        let rx = coords_rotated[b] - center[0];
        let ry = coords_rotated[b + 1] - center[1];
        let rz = coords_rotated[b + 2] - center[2];

        // Tangential term: a_tang = α × r
        let tang_x = ay * rz - az * ry;
        let tang_y = az * rx - ax * rz;
        let tang_z = ax * ry - ay * rx;

        // Centripetal term: a_cent = ω×(ω×r) = -ω² · r_perp (radially INWARD)
        // where r_perp = r - (r·n̂)·n̂ is the component of r perpendicular to the axis.
        // Note the negative sign: ω×(ω×r) points toward the axis (centripetal),
        // not away from it (centrifugal).
        let r_dot_axis = rx * axis[0] + ry * axis[1] + rz * axis[2];
        let rpx = rx - r_dot_axis * axis[0];
        let rpy = ry - r_dot_axis * axis[1];
        let rpz = rz - r_dot_axis * axis[2];

        let cent_x = -omega_sq * rpx;
        let cent_y = -omega_sq * rpy;
        let cent_z = -omega_sq * rpz;

        // Total reference acceleration: a_ref = a_tang + a_cent
        out[b] = tang_x + cent_x;
        out[b + 1] = tang_y + cent_y;
        out[b + 2] = tang_z + cent_z;
    }
    out
}

/// Compute rigid-body velocity in the INERTIAL (global) frame.
///
/// For a node at position `r_i = X_rot_i − center` (where `X_rot_i` is the
/// current rigidly rotated coordinate), the rigid-body velocity is:
///
///   `v_rigid_i = ω × r_i`
///
/// This is needed when the inertial solver exports **total** kinematics to a
/// CFD participant. In `Elastic` mode the rigid part is intentionally omitted.
///
/// # Arguments
/// * `coords_rotated` — flat node coordinates `[x0,y0,z0,…]` already rotated
///   to the target angle θ, in the global frame.
/// * `axis` — rotation axis unit vector `n̂`
/// * `center` — rotation center coordinates
/// * `omega` — angular velocity ω [rad/s]
///
/// # Returns
/// Flat velocity array `[vx0,vy0,vz0,…]` in the global frame, length `n_nodes * 3`.
pub fn compute_rigid_body_velocity_inertial(
    coords_rotated: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
    omega: f64,
) -> Vec<f64> {
    let n = coords_rotated.len() / 3;
    let mut out = vec![0.0f64; n * 3];

    if omega.abs() < 1e-14 {
        return out;
    }

    let wx = omega * axis[0];
    let wy = omega * axis[1];
    let wz = omega * axis[2];

    for i in 0..n {
        let b = i * 3;
        let rx = coords_rotated[b] - center[0];
        let ry = coords_rotated[b + 1] - center[1];
        let rz = coords_rotated[b + 2] - center[2];

        out[b] = wy * rz - wz * ry;
        out[b + 1] = wz * rx - wx * rz;
        out[b + 2] = wx * ry - wy * rx;
    }

    out
}

/// Compute the reference load vector `F_ref = −M_lumped ⊙ a_ref`.
///
/// This converts the rigid-body reference acceleration field (computed by
/// `compute_rigid_body_acceleration_inertial`) into a nodal force vector
/// ready to be scattered into the reduced-DOF global force vector.
///
/// # Arguments
/// * `a_ref` — flat acceleration array `[ax0,ay0,az0,…]` length `n_nodes * 3`
/// * `masses` — per-node scalar masses, length `n_nodes`
///
/// # Returns
/// Flat force array `[fx0,fy0,fz0,…]` length `n_nodes * 3`, where each
/// component is `f[i*3+j] = -masses[i] * a_ref[i*3+j]`.
///
/// # Notes
/// The negative sign is critical: the load vector `F_ref` appears on the RHS
/// as `−M·a_ref` to account for the d'Alembert inertial force in the global frame.
pub fn compute_reference_load_vector(
    a_ref: &[f64],
    masses: &[f64],
) -> Vec<f64> {
    let n = masses.len();
    assert_eq!(
        a_ref.len(),
        n * 3,
        "a_ref must have length n_nodes * 3 (got {} for {} nodes)",
        a_ref.len(),
        n
    );

    let mut out = vec![0.0f64; n * 3];
    for i in 0..n {
        let b = i * 3;
        let m = masses[i];
        out[b] = -m * a_ref[b];
        out[b + 1] = -m * a_ref[b + 1];
        out[b + 2] = -m * a_ref[b + 2];
    }
    out
}

// ── Coriolis gyroscopic matrix ─────────────────────────────────────────────────

/// Build the Coriolis gyroscopic matrix G_cor (antisymmetric, 3×3 blocks).
///
/// For a rotating frame with angular velocity ω, Coriolis acceleration is:
///   a_cor = -2·ω × v
///
/// In matrix form: F_cor = G_cor · v, where G_cor is antisymmetric.
///
/// For each node with mass m_i, the 3×3 block is:
///   G_i = -2·m_i·Ω, where Ω is the skew-symmetric cross-product matrix:
///
///   Ω = [  0   -ωz   ωy ]
///       [  ωz   0   -ωx ]
///       [ -ωy   ωx   0  ]
///
/// This function returns the REDUCED matrix in COO format, ready to add to
/// the Newmark effective damping matrix for implicit treatment.
///
/// # Arguments
/// * `masses` — per-node scalar masses (length n_nodes)
/// * `axis` — rotation axis unit vector
/// * `omega` — angular velocity (rad/s)
/// * `free_dofs_i32` — global DOF indices of free DOFs (reduced system)
/// * `dofs_per_node` — typically 6 for shells (u,v,w,rx,ry,rz)
///
/// # Returns
/// * `(rows, cols, vals)` — COO format for the REDUCED Coriolis matrix.
///   Only includes entries where both row and col are in free_dofs.
pub fn build_coriolis_matrix(
    masses: &[f64],
    axis: &[f64; 3],
    omega: f64,
    free_dofs_i32: &[i32],
    dofs_per_node: usize,
) -> (Vec<i32>, Vec<i32>, Vec<f64>) {
    if omega.abs() < 1e-14 {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    // Build free DOF lookup: full_dof -> reduced_index
    let max_dof = free_dofs_i32.iter().max().copied().unwrap_or(0) as usize;
    let mut full_to_red = vec![-1i32; max_dof + 1];
    for (red_idx, &gdof) in free_dofs_i32.iter().enumerate() {
        if (gdof as usize) < full_to_red.len() {
            full_to_red[gdof as usize] = red_idx as i32;
        }
    }

    // Coriolis only couples translational DOFs (u,v,w), not rotational
    let n_nodes = masses.len();
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    // Skew-symmetric matrix Ω for ω × v
    let wx = omega * axis[0];
    let wy = omega * axis[1];
    let wz = omega * axis[2];

    // 3×3 block pattern (antisymmetric):
    // [ 0    -wz   wy ]
    // [ wz    0   -wx ]
    // [-wy   wx    0  ]
    
    for node in 0..n_nodes {
        let m = masses[node];
        let coeff = -2.0 * m;

        // Global DOF base for this node's translations
        let base = node * dofs_per_node;

        // 6 off-diagonal entries per node (diagonal is zero in antisymmetric matrix)
        let entries = [
            // Row 0 (u): couples with v, w
            (0, 1, -wz * coeff),  // G[u, v] = +2m·ωz
            (0, 2,  wy * coeff),  // G[u, w] = -2m·ωy
            // Row 1 (v): couples with u, w
            (1, 0,  wz * coeff),  // G[v, u] = -2m·ωz
            (1, 2, -wx * coeff),  // G[v, w] = +2m·ωx
            // Row 2 (w): couples with u, v
            (2, 0, -wy * coeff),  // G[w, u] = +2m·ωy
            (2, 1,  wx * coeff),  // G[w, v] = -2m·ωx
        ];

        for (local_row, local_col, value) in &entries {
            if value.abs() < 1e-20 {
                continue;
            }
            let gdof_row = base + local_row;
            let gdof_col = base + local_col;

            // Check if both are free DOFs
            if gdof_row < full_to_red.len() && gdof_col < full_to_red.len() {
                let red_row = full_to_red[gdof_row];
                let red_col = full_to_red[gdof_col];
                if red_row >= 0 && red_col >= 0 {
                    rows.push(red_row);
                    cols.push(red_col);
                    vals.push(*value);
                }
            }
        }
    }

    (rows, cols, vals)
}

// ── Spin-softening K_SP ────────────────────────────────────────────────────────

/// Build spin-softening `K_SP` aligned with the reduced K COO sparsity.
///
/// ANSYS Eq. 3-74 / 14-55 for lumped mass:
/// ```text
/// K_SP,node = −ω² · m_node · (I − n̂⊗n̂)   on the 3×3 translational block
/// K_SP,node = 0                            on rotational rows/cols and inter-node couplings
/// ```
///
/// For axes aligned with a global basis vector this reduces to a diagonal block.
/// For arbitrary axes the translational block contains off-diagonal couplings,
/// so treating `K_SP` as a pure diagonal is physically wrong.
///
/// # Arguments
/// * `m_lumped_full` — lumped mass diagonal of the **full** (unreduced) system,
///   length `n_full_dofs`. Only translational DOFs carry mass (rotational = 0).
/// * `axis`          — rotation axis unit vector
/// * `omega`         — angular velocity (rad/s)
/// * `dofs_per_node` — DOFs per FEM node (typically 6 for shells)
/// * `k_rows`, `k_cols` — COO sparsity of the **reduced** stiffness matrix
/// * `free_dofs`     — sorted free-DOF global indices (length = n_reduced_dofs)
///
/// # Returns
/// `Vec<f64>` of length `k_rows.len()` ready to pass to `update_spin_softening`.
#[allow(clippy::too_many_arguments)]
pub fn build_ksp_vals(
    m_lumped_full: &[f64],
    axis: &[f64; 3],
    omega: f64,
    dofs_per_node: usize,
    k_rows: &[i32],
    k_cols: &[i32],
    free_dofs: &[i32],
    n_full_dofs: usize,
) -> Vec<f64> {
    if omega == 0.0 {
        return vec![0.0; k_rows.len()];
    }

    let omega_sq = omega * omega;
    let projector = [
        [1.0 - axis[0] * axis[0], -axis[0] * axis[1], -axis[0] * axis[2]],
        [-axis[1] * axis[0], 1.0 - axis[1] * axis[1], -axis[1] * axis[2]],
        [-axis[2] * axis[0], -axis[2] * axis[1], 1.0 - axis[2] * axis[2]],
    ];

    k_rows
        .iter()
        .zip(k_cols.iter())
        .map(|(&r_red, &c_red)| {
            let (Some(&r_global_i32), Some(&c_global_i32)) =
                (free_dofs.get(r_red as usize), free_dofs.get(c_red as usize))
            else {
                return 0.0;
            };

            let r_global = r_global_i32 as usize;
            let c_global = c_global_i32 as usize;
            if r_global >= n_full_dofs || c_global >= n_full_dofs {
                return 0.0;
            }

            let r_local = r_global % dofs_per_node;
            let c_local = c_global % dofs_per_node;
            if r_local >= 3 || c_local >= 3 {
                return 0.0;
            }

            let r_node = r_global / dofs_per_node;
            let c_node = c_global / dofs_per_node;
            if r_node != c_node {
                return 0.0;
            }

            let m_row = m_lumped_full[r_global];
            -omega_sq * m_row * projector[r_local][c_local]
        })
        .collect()
}

// ── Torque computation ─────────────────────────────────────────────────────────

/// Compute the scalar driving torque about the rotation axis from interface nodal forces.
///
/// `τ_scalar = (Σᵢ (X₀,ᵢ + uᵢ − center) × Fᵢ) · n̂`
///
/// # Arguments
/// * `iface_coords` — flat interface node coordinates X₀ in the rotating frame,
///   length `n_nodes * 3`
/// * `iface_disps`  — flat interface displacements u, length `n_nodes * 3`
/// * `iface_forces` — flat interface forces F, length `n_nodes * 3`
/// * `axis`         — rotation axis unit vector
/// * `center`       — rotation center coordinates
///
/// # Returns
/// `(torque_vec [f64; 3], torque_scalar f64)`
pub fn compute_torque(
    iface_coords: &[f64],
    iface_disps: &[f64],
    iface_forces: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
) -> ([f64; 3], f64) {
    let n = iface_coords.len() / 3;
    let mut tau = [0.0f64; 3];
    for i in 0..n {
        let b = i * 3;
        // r = (X₀ + u) − center
        let rx = iface_coords[b] + iface_disps[b] - center[0];
        let ry = iface_coords[b + 1] + iface_disps[b + 1] - center[1];
        let rz = iface_coords[b + 2] + iface_disps[b + 2] - center[2];
        let fx = iface_forces[b];
        let fy = iface_forces[b + 1];
        let fz = iface_forces[b + 2];
        // r × F
        tau[0] += ry * fz - rz * fy;
        tau[1] += rz * fx - rx * fz;
        tau[2] += rx * fy - ry * fx;
    }
    let scalar = tau[0] * axis[0] + tau[1] * axis[1] + tau[2] * axis[2];
    (tau, scalar)
}

// ── Performance coefficients ───────────────────────────────────────────────────

/// Non-dimensional rotor performance coefficients.
#[derive(Debug, Clone, Copy)]
pub struct PerformanceCoefficients {
    pub ct: f64,
    pub cp: f64,
    pub cq: f64,
    pub tsr: f64,
}

/// Kinematic quantities used to keep the rotor angle update consistent over one time window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowKinematics {
    /// Representative angular velocity held constant during the converged window.
    pub omega_step: f64,
    /// Angular acceleration used in the Euler body-force term during the window.
    pub alpha_step: f64,
    /// Exact or approximated angular increment over the full window.
    pub theta_increment: f64,
}

fn ramp_window_kinematics(omega_target: f64, t_ramp: f64, t: f64, dt: f64) -> WindowKinematics {
    if dt <= 0.0 {
        let omega = if t >= t_ramp {
            omega_target
        } else {
            omega_target * (t / t_ramp)
        };
        let alpha = if t < t_ramp { omega_target / t_ramp } else { 0.0 };
        return WindowKinematics {
            omega_step: omega,
            alpha_step: alpha,
            theta_increment: 0.0,
        };
    }

    if t_ramp <= 0.0 || t >= t_ramp {
        return WindowKinematics {
            omega_step: omega_target,
            alpha_step: 0.0,
            theta_increment: omega_target * dt,
        };
    }

    let alpha = omega_target / t_ramp;
    let omega_start = omega_target * (t / t_ramp);
    let t_end = t + dt;

    let theta_increment = if t_end <= t_ramp {
        omega_start * dt + 0.5 * alpha * dt * dt
    } else {
        let dt_ramp = t_ramp - t;
        let dtheta_ramp = omega_start * dt_ramp + 0.5 * alpha * dt_ramp * dt_ramp;
        let dtheta_const = omega_target * (dt - dt_ramp);
        dtheta_ramp + dtheta_const
    };

    WindowKinematics {
        omega_step: theta_increment / dt,
        alpha_step: alpha,
        theta_increment,
    }
}

/// Compute rotor performance coefficients.
///
/// Standard aerodynamic definitions:
/// ```text
/// Ct  = Thrust      / (½·ρ·V∞²·π·R²)
/// Cp  = P_aero      / (½·ρ·V∞³·π·R²)
/// Cq  = τ_aero      / (½·ρ·V∞²·π·R²·R)
/// TSR = |ω|·R / V∞
/// ```
///
/// All coefficients use AERODYNAMIC forces only (not total).
///
/// # Arguments
/// * `thrust`      — axial thrust [N] (aero force projected on rotation axis)
/// * `power_aero`  — aerodynamic power [W] = τ_aero · ω
/// * `torque_aero` — aerodynamic torque [N·m]
/// * `omega`       — angular velocity [rad/s]
/// * `radius`      — rotor radius [m]
/// * `fluid_density`   — fluid density [kg/m³]
/// * `flow_velocity`   — freestream velocity [m/s]
pub fn compute_performance_coefficients(
    thrust: f64,
    power_aero: f64,
    torque_aero: f64,
    omega: f64,
    radius: f64,
    fluid_density: f64,
    flow_velocity: f64,
) -> PerformanceCoefficients {
    let min_denom = 1e-6_f64;
    let area = std::f64::consts::PI * radius * radius;
    let q_dyn = 0.5 * fluid_density * flow_velocity * flow_velocity;
    let denom_force = q_dyn * area;
    let denom_power = q_dyn * area * flow_velocity;
    let denom_torque = q_dyn * area * radius;

    // Sign convention: the structural solver uses RHR around +Y (omega > 0,
    // clockwise from -Y). The driving tangential force projects onto -X, so
    // tau_aero < 0 and power_aero < 0 for a wind-turbine extracting energy.
    // Ct, Cp, Cq are reported as positive quantities (standard HAWT convention)
    // by negating tau and power before normalisation.
    // Thrust projects onto +Y (axial) and is already positive — no sign flip.
    let ct = if denom_force.abs() > min_denom { thrust / denom_force } else { 0.0 };
    let cp = if denom_power.abs() > min_denom { -power_aero / denom_power } else { 0.0 };
    let cq = if denom_torque.abs() > min_denom { -torque_aero / denom_torque } else { 0.0 };
    let tsr = if flow_velocity.abs() > min_denom {
        omega.abs() * radius / flow_velocity
    } else {
        0.0
    };

    PerformanceCoefficients { ct, cp, cq, tsr }
}

// ── OmegaProvider ──────────────────────────────────────────────────────────────

/// Saved state for checkpoint/restore of an `OmegaProvider`.
#[derive(Debug, Clone)]
pub struct OmegaCheckpoint {
    pub omega: f64,
    pub alpha: f64,
    /// Previous alpha (t_{n-1}) for Adams-Bashforth 2 on restore. `None` on first step.
    pub alpha_prev: Option<f64>,
    /// For `RampedComputed`: whether the ramp phase has been completed.
    pub ramp_completed: bool,
    /// Internal time tracked by the provider (for `RampedComputed` phase detection).
    pub current_time: f64,
}

/// Angular-velocity provider for the rotor FSI solver.
///
/// The provider is updated **once per converged time window** (not per sub-iteration).
/// During implicit FSI sub-iterations, ω is held constant.
///
/// # Variants
/// * `Constant`      — ω = const, α = 0 always.
/// * `Ramped`        — linear ramp ω(t) = ω_target · min(t/t_ramp, 1).
/// * `Computed`      — dynamic ω from `I·dω/dt = τ_driving + τ_shaft`.
///   Uses **Adams-Bashforth 2** (AB2) once a previous α is available, falling
///   back to Euler on the first step.  This gives O(Δt²) accuracy in ω at
///   the cost of storing one extra scalar `alpha_prev`.
/// * `RampedComputed`— ramp until `t ≥ t_ramp`, then `Computed`.
#[derive(Debug, Clone)]
pub enum OmegaProvider {
    Constant {
        omega: f64,
    },
    Ramped {
        omega_target: f64,
        t_ramp: f64,
    },
    Computed {
        moment_of_inertia: f64,
        shaft_torque: f64,
        /// Current dynamic state.
        omega: f64,
        alpha: f64,
        /// Previous α (step n−1) for Adams-Bashforth 2. `None` until the
        /// second converged window (first step uses Euler fallback).
        alpha_prev: Option<f64>,
    },
    RampedComputed {
        omega_target: f64,
        t_ramp: f64,
        moment_of_inertia: f64,
        shaft_torque: f64,
        /// Current dynamic state (used after ramp completion).
        omega: f64,
        alpha: f64,
        /// Previous α for AB2 (set after ramp completion). `None` until the
        /// second post-ramp converged window.
        alpha_prev: Option<f64>,
        ramp_completed: bool,
        current_time: f64,
    },
}

impl OmegaProvider {
    /// Get the current (ω, α) at time `t`.
    ///
    /// For `Computed` and `RampedComputed` (post-ramp), `t` is ignored —
    /// the state is maintained internally and updated via [`Self::update_from_torque`].
    pub fn get(&self, t: f64) -> (f64, f64) {
        match self {
            OmegaProvider::Constant { omega } => (*omega, 0.0),
            OmegaProvider::Ramped { omega_target, t_ramp } => {
                if t >= *t_ramp {
                    (*omega_target, 0.0)
                } else {
                    let ratio = t / t_ramp;
                    (*omega_target * ratio, *omega_target / t_ramp)
                }
            }
            OmegaProvider::Computed { omega, alpha, .. } => (*omega, *alpha),
            OmegaProvider::RampedComputed {
                omega_target,
                t_ramp,
                omega,
                alpha,
                ramp_completed,
                ..
            } => {
                if !ramp_completed {
                    if t < *t_ramp {
                        let ratio = t / t_ramp;
                        (*omega_target * ratio, *omega_target / t_ramp)
                    } else {
                        (*omega_target, 0.0)
                    }
                } else {
                    (*omega, *alpha)
                }
            }
        }
    }

    /// Return the kinematics used during a single converged time window.
    ///
    /// The rotor solver holds a representative angular velocity constant during
    /// the window, while the cumulative angle is advanced with a second-order
    /// accurate increment whenever the provider supplies a meaningful `α`.
    pub fn kinematics_over_step(&self, t: f64, dt: f64) -> WindowKinematics {
        match self {
            OmegaProvider::Constant { omega } => WindowKinematics {
                omega_step: *omega,
                alpha_step: 0.0,
                theta_increment: *omega * dt,
            },
            OmegaProvider::Ramped { omega_target, t_ramp } => {
                ramp_window_kinematics(*omega_target, *t_ramp, t, dt)
            }
            OmegaProvider::Computed { omega, alpha, .. } => {
                let theta_increment = *omega * dt + 0.5 * *alpha * dt * dt;
                WindowKinematics {
                    omega_step: if dt > 0.0 { theta_increment / dt } else { *omega },
                    alpha_step: *alpha,
                    theta_increment,
                }
            }
            OmegaProvider::RampedComputed {
                omega_target,
                t_ramp,
                omega,
                alpha,
                ramp_completed,
                ..
            } => {
                if !ramp_completed {
                    ramp_window_kinematics(*omega_target, *t_ramp, t, dt)
                } else {
                    let theta_increment = *omega * dt + 0.5 * *alpha * dt * dt;
                    WindowKinematics {
                        omega_step: if dt > 0.0 { theta_increment / dt } else { *omega },
                        alpha_step: *alpha,
                        theta_increment,
                    }
                }
            }
        }
    }

    /// Integrate the equation of motion for one converged time window.
    ///
    /// `I · dω/dt = τ_driving + τ_shaft`  →  Euler step: `ω += α · dt`
    ///
    /// **`tau_driving` must equal τ_aero + τ_gravity.**
    /// Gravitational sag creates a non-zero torque about the rotation axis that
    /// contributes to the angular acceleration and must be included here.
    /// Fictitious forces (centrifugal, Coriolis, Euler) appear on the RHS of the
    /// structural equation but are NOT physical driving torques and must be excluded.
    ///
    /// For `Constant` and `Ramped` this is a no-op.
    /// For `RampedComputed`, ignored while still in the ramp phase.
    pub fn update_from_torque(&mut self, tau_driving: f64, dt: f64, t: f64) {
        match self {
            OmegaProvider::Computed {
                moment_of_inertia,
                shaft_torque,
                omega,
                alpha,
                alpha_prev,
            } => {
                let alpha_new = (tau_driving + *shaft_torque) / *moment_of_inertia;
                // Adams-Bashforth 2 when a previous step is available,
                // otherwise fall back to Euler on the first step.
                if let Some(a_prev) = *alpha_prev {
                    *omega += (1.5 * alpha_new - 0.5 * a_prev) * dt;
                } else {
                    *omega += alpha_new * dt;
                }
                *alpha_prev = Some(*alpha);
                *alpha = alpha_new;
            }
            OmegaProvider::RampedComputed {
                omega_target,
                t_ramp,
                moment_of_inertia,
                shaft_torque,
                omega,
                alpha,
                alpha_prev,
                ramp_completed,
                current_time,
            } => {
                *current_time = t;
                if !*ramp_completed {
                    if t >= *t_ramp {
                        // Transition to dynamic phase.
                        *omega = *omega_target;
                        *alpha = 0.0;
                        *alpha_prev = None; // Reset AB2 history at ramp-to-dynamic transition.
                        *ramp_completed = true;
                    }
                    // Still ramping — torque update ignored.
                    return;
                }
                let alpha_new = (tau_driving + *shaft_torque) / *moment_of_inertia;
                if let Some(a_prev) = *alpha_prev {
                    *omega += (1.5 * alpha_new - 0.5 * a_prev) * dt;
                } else {
                    *omega += alpha_new * dt;
                }
                *alpha_prev = Some(*alpha);
                *alpha = alpha_new;
            }
            _ => {}
        }
    }

    /// Save the current state for implicit coupling checkpoint/restore.
    pub fn checkpoint(&self) -> OmegaCheckpoint {
        match self {
            OmegaProvider::Constant { omega } => OmegaCheckpoint {
                omega: *omega,
                alpha: 0.0,
                alpha_prev: None,
                ramp_completed: true,
                current_time: 0.0,
            },
            OmegaProvider::Ramped { .. } => OmegaCheckpoint {
                omega: 0.0,
                alpha: 0.0,
                alpha_prev: None,
                ramp_completed: false,
                current_time: 0.0,
            },
            OmegaProvider::Computed { omega, alpha, alpha_prev, .. } => OmegaCheckpoint {
                omega: *omega,
                alpha: *alpha,
                alpha_prev: *alpha_prev,
                ramp_completed: true,
                current_time: 0.0,
            },
            OmegaProvider::RampedComputed {
                omega,
                alpha,
                alpha_prev,
                ramp_completed,
                current_time,
                ..
            } => OmegaCheckpoint {
                omega: *omega,
                alpha: *alpha,
                alpha_prev: *alpha_prev,
                ramp_completed: *ramp_completed,
                current_time: *current_time,
            },
        }
    }

    /// Restore the state from a checkpoint.
    pub fn restore(&mut self, cp: &OmegaCheckpoint) {
        match self {
            OmegaProvider::Computed { omega, alpha, alpha_prev, .. } => {
                *omega = cp.omega;
                *alpha = cp.alpha;
                *alpha_prev = cp.alpha_prev;
            }
            OmegaProvider::RampedComputed {
                omega,
                alpha,
                alpha_prev,
                ramp_completed,
                current_time,
                ..
            } => {
                *omega = cp.omega;
                *alpha = cp.alpha;
                *alpha_prev = cp.alpha_prev;
                *ramp_completed = cp.ramp_completed;
                *current_time = cp.current_time;
            }
            _ => {}
        }
    }

    /// Return the angular velocity at t=0 (used for initial K_SP / K_G assembly).
    pub fn initial_omega(&self) -> f64 {
        match self {
            OmegaProvider::Constant { omega } => *omega,
            OmegaProvider::Ramped { .. } => 0.0,
            OmegaProvider::Computed { omega, .. } => *omega,
            OmegaProvider::RampedComputed { .. } => 0.0,
        }
    }
}

// ── Thrust calculation ──────────────────────────────────────────────────────────

/// Compute the axial thrust (aerodynamic force projected on the rotation axis).
///
/// `Thrust = (Σᵢ F_aero,i) · n̂`
///
/// # Arguments
/// * `aero_forces` — flat aerodynamic force array in the **inertial** frame,
///   length `n_nodes * 3`
/// * `axis`        — rotation axis unit vector
pub fn compute_thrust(aero_forces: &[f64], axis: &[f64; 3]) -> f64 {
    let n = aero_forces.len() / 3;
    let mut fsum = [0.0f64; 3];
    for i in 0..n {
        let b = i * 3;
        fsum[0] += aero_forces[b];
        fsum[1] += aero_forces[b + 1];
        fsum[2] += aero_forces[b + 2];
    }
    fsum[0] * axis[0] + fsum[1] * axis[1] + fsum[2] * axis[2]
}

// ── Moment of inertia ────────────────────────────────────────────────────────────

/// Estimate the moment of inertia about the rotation axis from nodal masses
/// and coordinates (lumped mass approximation).
///
/// `I = Σᵢ mᵢ · |r_⊥,i|²`
///
/// # Arguments
/// * `coords` — flat node coordinates in the rotating frame, length `n_nodes * 3`
/// * `masses` — per-node scalar masses, length `n_nodes`
/// * `axis`   — rotation axis unit vector
/// * `center` — rotation center
pub fn estimate_moment_of_inertia(
    coords: &[f64],
    masses: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
) -> f64 {
    let n = masses.len();
    let mut inertia = 0.0f64;
    for i in 0..n {
        let b = i * 3;
        let rx = coords[b] - center[0];
        let ry = coords[b + 1] - center[1];
        let rz = coords[b + 2] - center[2];
        let rda = rx * axis[0] + ry * axis[1] + rz * axis[2];
        let rpx = rx - rda * axis[0];
        let rpy = ry - rda * axis[1];
        let rpz = rz - rda * axis[2];
        inertia += masses[i] * (rpx * rpx + rpy * rpy + rpz * rpz);
    }
    inertia
}

// ── Geometric stiffness K_G(u) helpers ─────────────────────────────────────────

/// Compose deformed nodal coordinates from rotated reference coords + nodal
/// translation DOFs.
///
/// `coords_def[3·i + k] = coords[3·i + k] + u_full[i·dofs_per_node + k]`
/// for `k ∈ {0,1,2}` (the translational DOFs). Rotational DOFs in `u_full`
/// are ignored.
///
/// Used to assemble `K_G(θ, u)` capturing the foreshortening contribution from
/// the elastic deformation. Call sites:
/// - `crates/aeroelast-solvers/src/petsc/fsi/rotor_inertial.rs` (inertial solver)
/// - `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs`     (corotational solver)
///
/// # Arguments
/// * `coords`         — flat node coordinates `[x0,y0,z0,…]`, length `n_nodes * 3`
///                      (typically the rotated reference geometry `X_ref · R(θ)`)
/// * `u_full`         — flat full DOF displacement vector,
///                      length `n_nodes * dofs_per_node`
/// * `dofs_per_node`  — DOFs per FEM node (must be ≥ 3)
///
/// # Returns
/// Flat deformed coords `[x0+ux0, y0+uy0, z0+uz0, …]`, length `n_nodes * 3`.
///
/// # Panics
/// If `dofs_per_node < 3` or if the nodal counts derived from `coords` and
/// `u_full` disagree.
pub fn compose_deformed_coords(
    coords: &[f64],
    u_full: &[f64],
    dofs_per_node: usize,
) -> Vec<f64> {
    assert!(
        dofs_per_node >= 3,
        "dofs_per_node must be >= 3 for compose_deformed_coords, got {}",
        dofs_per_node
    );
    let n_nodes = coords.len() / 3;
    assert_eq!(
        n_nodes,
        u_full.len() / dofs_per_node,
        "compose_deformed_coords: coords ({} nodes) and u_full ({} nodes) disagree",
        n_nodes,
        u_full.len() / dofs_per_node,
    );
    let mut out = vec![0.0f64; n_nodes * 3];
    for i in 0..n_nodes {
        let cb = i * 3;
        let ub = i * dofs_per_node;
        out[cb] = coords[cb] + u_full[ub];
        out[cb + 1] = coords[cb + 1] + u_full[ub + 1];
        out[cb + 2] = coords[cb + 2] + u_full[ub + 2];
    }
    out
}

/// Compute the maximum perpendicular-to-axis deflection ratio.
///
/// `ratio = max_i ||u_⊥,i|| / max_j ||r_⊥,j||`
///
/// where the perpendicular projections are taken with respect to the rotation
/// axis `n̂`:
/// - `u_⊥,i = u_trans,i − (u_trans,i · n̂)·n̂`  (perpendicular component of the
///   nodal translation)
/// - `r_⊥,j = (X_j − center) − ((X_j − center) · n̂)·n̂`  (perpendicular
///   radial position of node `j` relative to the rotation axis)
///
/// Acts as the trigger metric for `K_G(θ, u)` reassembly: comparing the change
/// in this dimensionless ratio across windows captures foreshortening-dominant
/// deformation evolution while remaining invariant to absolute blade size.
///
/// # Arguments
/// * `u_full`         — flat full DOF displacement, length `n_nodes * dofs_per_node`
/// * `coords`         — flat node coordinates `[x0,y0,z0,…]`, length `n_nodes * 3`
///                      (rotated reference geometry, NOT deformed)
/// * `axis`           — rotation axis unit vector
/// * `center`         — rotation center
/// * `dofs_per_node`  — DOFs per FEM node (must be ≥ 3)
///
/// # Returns
/// Dimensionless ratio in `[0, +∞)`. Returns `0.0` when:
/// - all nodes are on the rotation axis (`max_r_⊥ ≈ 0`), or
/// - `u_full` is the zero vector.
///
/// # Panics
/// If `dofs_per_node < 3` or if the nodal counts derived from `coords` and
/// `u_full` disagree.
pub fn max_radial_deflection_ratio(
    u_full: &[f64],
    coords: &[f64],
    axis: &[f64; 3],
    center: &[f64; 3],
    dofs_per_node: usize,
) -> f64 {
    assert!(
        dofs_per_node >= 3,
        "dofs_per_node must be >= 3 for max_radial_deflection_ratio, got {}",
        dofs_per_node
    );
    let n_nodes = coords.len() / 3;
    if n_nodes == 0 {
        return 0.0;
    }
    assert_eq!(
        n_nodes,
        u_full.len() / dofs_per_node,
        "max_radial_deflection_ratio: coords ({} nodes) and u_full ({} nodes) disagree",
        n_nodes,
        u_full.len() / dofs_per_node,
    );

    let mut max_u_perp_sq = 0.0f64;
    let mut max_r_perp_sq = 0.0f64;

    for i in 0..n_nodes {
        let cb = i * 3;
        let ub = i * dofs_per_node;

        // r_perp = (coords - center) - ((coords - center) · n̂) · n̂
        let rx = coords[cb] - center[0];
        let ry = coords[cb + 1] - center[1];
        let rz = coords[cb + 2] - center[2];
        let rda = rx * axis[0] + ry * axis[1] + rz * axis[2];
        let rpx = rx - rda * axis[0];
        let rpy = ry - rda * axis[1];
        let rpz = rz - rda * axis[2];
        let r_perp_sq = rpx * rpx + rpy * rpy + rpz * rpz;
        if r_perp_sq > max_r_perp_sq {
            max_r_perp_sq = r_perp_sq;
        }

        // u_perp = u_trans - (u_trans · n̂) · n̂
        let ux = u_full[ub];
        let uy = u_full[ub + 1];
        let uz = u_full[ub + 2];
        let uda = ux * axis[0] + uy * axis[1] + uz * axis[2];
        let upx = ux - uda * axis[0];
        let upy = uy - uda * axis[1];
        let upz = uz - uda * axis[2];
        let u_perp_sq = upx * upx + upy * upy + upz * upz;
        if u_perp_sq > max_u_perp_sq {
            max_u_perp_sq = u_perp_sq;
        }
    }

    if max_r_perp_sq < 1e-30 {
        return 0.0;
    }

    (max_u_perp_sq / max_r_perp_sq).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotation_matrix_identity_at_zero() {
        let t = RotorTransforms::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let r = t.rotation_matrix(0.0);
        let expected = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for i in 0..3 {
            for j in 0..3 {
                assert!((r[i][j] - expected[i][j]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn rotation_90deg_z_axis() {
        // R(90°) about Z: X→Y, Y→-X, Z→Z
        use std::f64::consts::FRAC_PI_2;
        let t = RotorTransforms::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let r = t.rotation_matrix(FRAC_PI_2);
        // R·[1,0,0] should be [0,1,0]
        let rx0 = r[0][0]; // should be ≈0
        let ry0 = r[1][0]; // should be ≈1
        assert!(rx0.abs() < 1e-14, "R[0][0] = {rx0}");
        assert!((ry0 - 1.0).abs() < 1e-14, "R[1][0] = {ry0}");
    }

    #[test]
    fn centrifugal_zero_at_rest() {
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let coords = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let masses = vec![1.0, 1.0];
        let f = compute_centrifugal_force(&coords, &masses, &axis, &center, 0.0);
        for v in &f {
            assert!(*v == 0.0);
        }
    }

    #[test]
    fn centrifugal_radially_outward() {
        // Node at [1, 0, 0], axis Z, omega=2.0 → r_perp=[1,0,0]
        // F = m * omega^2 * r_perp = 1 * 4 * [1,0,0]
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let coords = vec![1.0, 0.0, 0.0];
        let masses = vec![1.0];
        let f = compute_centrifugal_force(&coords, &masses, &axis, &center, 2.0);
        assert!((f[0] - 4.0).abs() < 1e-13, "fx={}", f[0]);
        assert!(f[1].abs() < 1e-14, "fy={}", f[1]);
        assert!(f[2].abs() < 1e-14, "fz={}", f[2]);
    }

    #[test]
    fn coriolis_force_z_axis() {
        // Node velocity [1,0,0], axis Z, omega=1 → ω×v=[0,0,1]×[1,0,0]=[0*0-1*0, 1*1-0*0, 0*0-0*1]
        // ω=[0,0,1], v=[1,0,0] → ω×v = [0*0-1*0, 1*1-0*0, 0*0-0*1] = [0,1,0]
        // F_cor = -2*1*[0,1,0] = [0,-2,0]
        let axis = [0.0, 0.0, 1.0];
        let velocities = vec![1.0, 0.0, 0.0];
        let masses = vec![1.0];
        let f = compute_coriolis_force(&velocities, &masses, &axis, 1.0);
        assert!(f[0].abs() < 1e-14, "fx={}", f[0]);
        assert!((f[1] + 2.0).abs() < 1e-14, "fy={}", f[1]);
        assert!(f[2].abs() < 1e-14, "fz={}", f[2]);
    }

    #[test]
    fn torque_single_node() {
        // Node at [1,0,0], force=[0,0,1] → r×F = [1,0,0]×[0,0,1] = [0*1-0*0, 0*0-1*1, 1*0-0*0] = [0,-1,0]
        // axis=[0,0,1] → τ_scalar = [0,-1,0]·[0,0,1] = 0
        let coords = vec![1.0, 0.0, 0.0];
        let disps = vec![0.0, 0.0, 0.0];
        let forces = vec![0.0, 0.0, 1.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let (tau, scalar) = compute_torque(&coords, &disps, &forces, &axis, &center);
        assert!((tau[1] + 1.0).abs() < 1e-14, "tau_y={}", tau[1]);
        assert!(scalar.abs() < 1e-14, "scalar={scalar}");
    }

    #[test]
    fn torque_uses_deformed_coords_and_rotation_center() {
        // center = [1,0,0], X0 = [2,0,0], u = [0,1,0], F = [10,0,0]
        // r = (X0 + u - center) = [1,1,0]
        // r × F = [1,1,0] × [10,0,0] = [0,0,-10]
        // τ_scalar about +Z = -10
        let coords = vec![2.0, 0.0, 0.0];
        let disps = vec![0.0, 1.0, 0.0];
        let forces = vec![10.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [1.0, 0.0, 0.0];

        let (tau, scalar) = compute_torque(&coords, &disps, &forces, &axis, &center);

        assert!(tau[0].abs() < 1e-14, "tau_x={}", tau[0]);
        assert!(tau[1].abs() < 1e-14, "tau_y={}", tau[1]);
        assert!((tau[2] + 10.0).abs() < 1e-14, "tau_z={}", tau[2]);
        assert!((scalar + 10.0).abs() < 1e-14, "scalar={scalar}");
    }

    #[test]
    fn omega_provider_constant() {
        let p = OmegaProvider::Constant { omega: 5.0 };
        assert_eq!(p.get(0.0), (5.0, 0.0));
        assert_eq!(p.get(100.0), (5.0, 0.0));
    }

    #[test]
    fn omega_provider_ramped() {
        let p = OmegaProvider::Ramped { omega_target: 10.0, t_ramp: 5.0 };
        let (w, a) = p.get(2.5);
        assert!((w - 5.0).abs() < 1e-12, "w={w}");
        assert!((a - 2.0).abs() < 1e-12, "a={a}");
        let (w2, a2) = p.get(10.0);
        assert_eq!((w2, a2), (10.0, 0.0));
    }

    #[test]
    fn omega_provider_ramped_window_kinematics_inside_ramp() {
        let p = OmegaProvider::Ramped { omega_target: 10.0, t_ramp: 5.0 };
        let kin = p.kinematics_over_step(2.0, 1.0);
        assert!((kin.omega_step - 5.0).abs() < 1e-12, "omega_step={}", kin.omega_step);
        assert!((kin.alpha_step - 2.0).abs() < 1e-12, "alpha_step={}", kin.alpha_step);
        assert!((kin.theta_increment - 5.0).abs() < 1e-12, "dtheta={}", kin.theta_increment);
    }

    #[test]
    fn omega_provider_ramped_window_kinematics_crossing_ramp_end() {
        let p = OmegaProvider::Ramped { omega_target: 10.0, t_ramp: 5.0 };
        let kin = p.kinematics_over_step(4.5, 1.0);
        assert!((kin.omega_step - 9.75).abs() < 1e-12, "omega_step={}", kin.omega_step);
        assert!((kin.alpha_step - 2.0).abs() < 1e-12, "alpha_step={}", kin.alpha_step);
        assert!((kin.theta_increment - 9.75).abs() < 1e-12, "dtheta={}", kin.theta_increment);
    }

    #[test]
    fn omega_provider_computed_update() {
        let mut p = OmegaProvider::Computed {
            moment_of_inertia: 100.0,
            shaft_torque: 0.0,
            omega: 0.0,
            alpha: 0.0,
            alpha_prev: None,
        };
        // tau_driving = 100 N·m, dt = 1.0 → alpha = 1.0, omega = 1.0 (Euler first step)
        p.update_from_torque(100.0, 1.0, 0.0);
        let (w, a) = p.get(0.0);
        assert!((w - 1.0).abs() < 1e-12, "w={w}");
        assert!((a - 1.0).abs() < 1e-12, "a={a}");
    }

    #[test]
    fn omega_provider_ramped_computed_ignores_torque_during_ramp() {
        let mut p = OmegaProvider::RampedComputed {
            omega_target: 10.0,
            t_ramp: 5.0,
            moment_of_inertia: 100.0,
            shaft_torque: 0.0,
            omega: 0.0,
            alpha: 0.0,
            alpha_prev: None,
            ramp_completed: false,
            current_time: 0.0,
        };

        p.update_from_torque(200.0, 1.0, 2.0);

        let cp = p.checkpoint();
        assert!(!cp.ramp_completed, "torque must be ignored while still ramping");
        assert!((cp.omega - 0.0).abs() < 1e-12, "stored omega changed during ramp");
        assert!((cp.alpha - 0.0).abs() < 1e-12, "stored alpha changed during ramp");
        assert!((cp.current_time - 2.0).abs() < 1e-12, "current_time not tracked during ramp");
    }

    #[test]
    fn omega_provider_ramped_computed_updates_after_ramp_completion() {
        let mut p = OmegaProvider::RampedComputed {
            omega_target: 10.0,
            t_ramp: 5.0,
            moment_of_inertia: 100.0,
            shaft_torque: 0.0,
            omega: 0.0,
            alpha: 0.0,
            alpha_prev: None,
            ramp_completed: false,
            current_time: 0.0,
        };

        // First call completes the ramp and seeds the dynamic phase.
        p.update_from_torque(200.0, 1.0, 5.0);
        let cp_after_ramp = p.checkpoint();
        assert!(cp_after_ramp.ramp_completed, "ramp should complete at t=t_ramp");
        assert!((cp_after_ramp.omega - 10.0).abs() < 1e-12, "omega should jump to target at ramp completion");

        // Next call must use the dynamic torque update (Euler on first dynamic step).
        p.update_from_torque(200.0, 1.0, 6.0);
        let (w, a) = p.get(6.0);
        assert!((w - 12.0).abs() < 1e-12, "w={w}");
        assert!((a - 2.0).abs() < 1e-12, "a={a}");
    }

    #[test]
    fn omega_provider_computed_window_kinematics_is_second_order() {
        let p = OmegaProvider::Computed {
            moment_of_inertia: 100.0,
            shaft_torque: 0.0,
            omega: 5.0,
            alpha: 2.0,
            alpha_prev: None,
        };
        let kin = p.kinematics_over_step(0.0, 0.5);
        assert!((kin.omega_step - 5.5).abs() < 1e-12, "omega_step={}", kin.omega_step);
        assert!((kin.alpha_step - 2.0).abs() < 1e-12, "alpha_step={}", kin.alpha_step);
        assert!((kin.theta_increment - 2.75).abs() < 1e-12, "dtheta={}", kin.theta_increment);
    }

    #[test]
    fn omega_provider_checkpoint_restore() {
        let mut p = OmegaProvider::Computed {
            moment_of_inertia: 100.0,
            shaft_torque: 0.0,
            omega: 5.0,
            alpha: 1.0,
            alpha_prev: None,
        };
        let cp = p.checkpoint();
        p.update_from_torque(200.0, 1.0, 0.0); // changes state
        p.restore(&cp);
        let (w, a) = p.get(0.0);
        assert!((w - 5.0).abs() < 1e-12, "w={w}");
        assert!((a - 1.0).abs() < 1e-12, "a={a}");
    }

    // ── Inertial-frame physics tests ───────────────────────────────────────────

    #[test]
    fn rigid_body_acceleration_pure_rotation() {
        // Pure rotation (ω≠0, α=0) → centripetal only: a = ω×(ω×r) = -ω²·r_perp (inward)
        // Node at [1,0,0], axis Z, center origin, ω=2.0:
        //   ω×r = (0,0,2)×(1,0,0) = (0,2,0)
        //   ω×(ω×r) = (0,0,2)×(0,2,0) = (-4,0,0)  ← inward (-x direction)
        let coords = vec![1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let omega = 2.0;
        let alpha = 0.0;

        let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

        assert_eq!(a_ref.len(), 3);
        assert!((a_ref[0] + 4.0).abs() < 1e-13, "ax={}", a_ref[0]);
        assert!(a_ref[1].abs() < 1e-14, "ay={}", a_ref[1]);
        assert!(a_ref[2].abs() < 1e-14, "az={}", a_ref[2]);
    }

    #[test]
    fn rigid_body_acceleration_pure_angular_acceleration() {
        // Pure angular acceleration (ω=0, α≠0) → tangential only: a = α×r
        // Node at [1,0,0], axis Z, center origin, α=3.0 → a = [0,0,3]×[1,0,0] = [0,3,0]
        let coords = vec![1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let omega = 0.0;
        let alpha = 3.0;

        let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

        assert_eq!(a_ref.len(), 3);
        assert!(a_ref[0].abs() < 1e-14, "ax={}", a_ref[0]);
        assert!((a_ref[1] - 3.0).abs() < 1e-13, "ay={}", a_ref[1]);
        assert!(a_ref[2].abs() < 1e-14, "az={}", a_ref[2]);
    }

    #[test]
    fn rigid_body_acceleration_both_components() {
        // Combined ω and α → tangential + centripetal
        // Node at [1,0,0], axis Z, ω=1.0, α=2.0
        // Tangential: α×r = [0,0,2]×[1,0,0] = [0,2,0]
        // Centripetal: ω×(ω×r) = -ω²·r_perp = -1·[1,0,0] = [-1,0,0]  (inward)
        // Total: [-1,2,0]
        let coords = vec![1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let omega = 1.0;
        let alpha = 2.0;

        let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, omega, alpha);

        assert_eq!(a_ref.len(), 3);
        assert!((a_ref[0] + 1.0).abs() < 1e-13, "ax={}", a_ref[0]);
        assert!((a_ref[1] - 2.0).abs() < 1e-13, "ay={}", a_ref[1]);
        assert!(a_ref[2].abs() < 1e-14, "az={}", a_ref[2]);
    }

    #[test]
    fn rigid_body_acceleration_zero_at_rest() {
        // ω=0, α=0 → early return with zeros
        let coords = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];

        let a_ref = compute_rigid_body_acceleration_inertial(&coords, &axis, &center, 0.0, 0.0);

        assert_eq!(a_ref.len(), 6);
        for v in &a_ref {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn rigid_body_velocity_pure_rotation() {
        // Node at [1,0,0], axis Z, ω=2 → v = [0,0,2]×[1,0,0] = [0,2,0]
        let coords = vec![1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];

        let v_ref = compute_rigid_body_velocity_inertial(&coords, &axis, &center, 2.0);

        assert_eq!(v_ref.len(), 3);
        assert!(v_ref[0].abs() < 1e-14, "vx={}", v_ref[0]);
        assert!((v_ref[1] - 2.0).abs() < 1e-13, "vy={}", v_ref[1]);
        assert!(v_ref[2].abs() < 1e-14, "vz={}", v_ref[2]);
    }

    #[test]
    fn rigid_body_velocity_zero_at_rest() {
        let coords = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];

        let v_ref = compute_rigid_body_velocity_inertial(&coords, &axis, &center, 0.0);

        assert_eq!(v_ref.len(), 6);
        for v in &v_ref {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn reference_load_vector_single_node() {
        // F_ref = -M·a_ref for a_ref=[2,3,4], m=5 → F=[-10,-15,-20]
        let a_ref = vec![2.0, 3.0, 4.0];
        let masses = vec![5.0];

        let f_ref = compute_reference_load_vector(&a_ref, &masses);

        assert_eq!(f_ref.len(), 3);
        assert!((f_ref[0] + 10.0).abs() < 1e-13, "fx={}", f_ref[0]);
        assert!((f_ref[1] + 15.0).abs() < 1e-13, "fy={}", f_ref[1]);
        assert!((f_ref[2] + 20.0).abs() < 1e-13, "fz={}", f_ref[2]);
    }

    #[test]
    fn reference_load_vector_multi_node() {
        // Two nodes: a_ref=[1,0,0, 0,2,0], masses=[3,4]
        // F_ref = [-3,0,0, 0,-8,0]
        let a_ref = vec![1.0, 0.0, 0.0, 0.0, 2.0, 0.0];
        let masses = vec![3.0, 4.0];

        let f_ref = compute_reference_load_vector(&a_ref, &masses);

        assert_eq!(f_ref.len(), 6);
        assert!((f_ref[0] + 3.0).abs() < 1e-13, "fx0={}", f_ref[0]);
        assert!(f_ref[1].abs() < 1e-14, "fy0={}", f_ref[1]);
        assert!(f_ref[2].abs() < 1e-14, "fz0={}", f_ref[2]);
        assert!(f_ref[3].abs() < 1e-14, "fx1={}", f_ref[3]);
        assert!((f_ref[4] + 8.0).abs() < 1e-13, "fy1={}", f_ref[4]);
        assert!(f_ref[5].abs() < 1e-14, "fz1={}", f_ref[5]);
    }

    #[test]
    fn build_ksp_vals_includes_off_diagonal_terms_for_tilted_axis() {
        let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
        let axis = [inv_sqrt2, inv_sqrt2, 0.0];
        let m_lumped_full = vec![3.0, 3.0, 3.0, 0.0, 0.0, 0.0];
        let k_rows = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let k_cols = vec![0, 1, 2, 0, 1, 2, 0, 1, 2];
        let free_dofs = vec![0, 1, 2];

        let ksp = build_ksp_vals(
            &m_lumped_full,
            &axis,
            2.0,
            6,
            &k_rows,
            &k_cols,
            &free_dofs,
            6,
        );

        let expected = [-6.0, 6.0, 0.0, 6.0, -6.0, 0.0, 0.0, 0.0, -12.0];
        for (actual, expected) in ksp.iter().zip(expected.iter()) {
            assert!((actual - expected).abs() < 1e-12, "actual={actual}, expected={expected}");
        }
    }

    // ── compose_deformed_coords ────────────────────────────────────────────────

    #[test]
    fn compose_deformed_coords_zero_u_returns_reference() {
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let u_full = vec![0.0; 12]; // 2 nodes × 6 DOFs (SHELL)
        let out = compose_deformed_coords(&coords, &u_full, 6);
        assert_eq!(out, coords);
    }

    #[test]
    fn compose_deformed_coords_dofs_per_node_3_uses_full_u() {
        // dofs_per_node = 3 (e.g., SOLID): all three components of u contribute
        let coords = vec![1.0, 0.0, 0.0];
        let u_full = vec![0.5, -0.2, 0.1];
        let out = compose_deformed_coords(&coords, &u_full, 3);
        assert!((out[0] - 1.5).abs() < 1e-14, "x={}", out[0]);
        assert!((out[1] + 0.2).abs() < 1e-14, "y={}", out[1]);
        assert!((out[2] - 0.1).abs() < 1e-14, "z={}", out[2]);
    }

    #[test]
    fn compose_deformed_coords_dofs_per_node_6_ignores_rotational_dofs() {
        // SHELL layout: u_full = [u, v, w, θx, θy, θz] per node
        // Only the translational DOFs (first 3) must update coords; rotations are ignored.
        let coords = vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let u_full = vec![
            0.1, 0.2, 0.3, 1.0, 2.0, 3.0, // node 0: trans + rot
            0.4, 0.5, 0.6, 4.0, 5.0, 6.0, // node 1: trans + rot
        ];
        let out = compose_deformed_coords(&coords, &u_full, 6);
        let expected = vec![1.1, 0.2, 0.3, 2.4, 0.5, 0.6];
        for (a, e) in out.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 1e-14, "got {a}, expected {e}");
        }
    }

    #[test]
    #[should_panic(expected = "dofs_per_node must be >= 3 for compose_deformed_coords")]
    fn compose_deformed_coords_panics_on_low_dofs_per_node() {
        let coords = vec![0.0, 0.0, 0.0];
        let u_full = vec![0.0, 0.0];
        let _ = compose_deformed_coords(&coords, &u_full, 2);
    }

    #[test]
    #[should_panic(expected = "compose_deformed_coords: coords")]
    fn compose_deformed_coords_panics_on_node_count_mismatch() {
        // 2 nodes in coords, 1 node in u_full
        let coords = vec![1.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let u_full = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // only 1 node × 6 DOFs
        let _ = compose_deformed_coords(&coords, &u_full, 6);
    }

    // ── max_radial_deflection_ratio ────────────────────────────────────────────

    #[test]
    fn max_radial_deflection_ratio_zero_u_returns_zero() {
        let coords = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let u_full = vec![0.0; 12];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!(r.abs() < 1e-14, "got {r}");
    }

    #[test]
    fn max_radial_deflection_ratio_axial_u_returns_zero() {
        // Axis Z, u purely along Z → u_perp = 0 (axial sliding doesn't affect K_G)
        let coords = vec![1.0, 0.0, 0.0];
        let u_full = vec![0.0, 0.0, 0.5, 0.0, 0.0, 0.0]; // dofs_per_node=6
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!(r.abs() < 1e-14, "got {r}");
    }

    #[test]
    fn max_radial_deflection_ratio_radial_u_single_node() {
        // Node at (1,0,0), axis Z. u = (0.1, 0, 0):
        //   u_perp = (0.1, 0, 0), |u_perp| = 0.1
        //   r_perp = (1, 0, 0),   |r_perp| = 1.0
        //   ratio = 0.1
        let coords = vec![1.0, 0.0, 0.0];
        let u_full = vec![0.1, 0.0, 0.0, 0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!((r - 0.1).abs() < 1e-13, "got {r}");
    }

    #[test]
    fn max_radial_deflection_ratio_uses_independent_max_u_and_max_r() {
        // Node 0: r_perp=(2,0,0), u_perp=(0.2,0,0)
        // Node 1: r_perp=(1,0,0), u_perp=(0.5,0,0)
        // max|u_perp|² = 0.25 (node 1), max|r_perp|² = 4.0 (node 0)
        // ratio = sqrt(0.25 / 4.0) = 0.25
        let coords = vec![2.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let u_full = vec![
            0.2, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.5, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!((r - 0.25).abs() < 1e-13, "got {r}");
    }

    #[test]
    fn max_radial_deflection_ratio_all_nodes_on_axis_returns_zero() {
        // All nodes on the rotation axis → max_r_perp ≈ 0 → guarded return 0
        let coords = vec![0.0, 0.0, 1.0, 0.0, 0.0, 2.0];
        let u_full = vec![1.0; 12]; // arbitrary nonzero u
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert_eq!(r, 0.0);
    }

    #[test]
    fn max_radial_deflection_ratio_respects_center_offset() {
        // Center at (10,0,0), node at (11,0,0), axis Z
        //   r_perp = (11-10, 0, 0) = (1, 0, 0), |r_perp| = 1
        //   u_perp = (0.3, 0, 0)
        //   ratio = 0.3
        let coords = vec![11.0, 0.0, 0.0];
        let u_full = vec![0.3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [10.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!((r - 0.3).abs() < 1e-13, "got {r}");
    }

    #[test]
    fn max_radial_deflection_ratio_tilted_axis() {
        // Axis = (1,0,0)/√1 = X. Node at (0, 1, 0), u_trans = (0.4, 0.3, 0).
        //   r_perp from X axis = (0, 1, 0) − 0·(1,0,0) = (0, 1, 0), |r_perp| = 1
        //   u_perp = (0.4, 0.3, 0) − 0.4·(1,0,0) = (0, 0.3, 0), |u_perp| = 0.3
        //   ratio = 0.3
        let coords = vec![0.0, 1.0, 0.0];
        let u_full = vec![0.4, 0.3, 0.0, 0.0, 0.0, 0.0];
        let axis = [1.0, 0.0, 0.0];
        let center = [0.0, 0.0, 0.0];
        let r = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 6);
        assert!((r - 0.3).abs() < 1e-13, "got {r}");
    }

    #[test]
    #[should_panic(expected = "dofs_per_node must be >= 3 for max_radial_deflection_ratio")]
    fn max_radial_deflection_ratio_panics_on_low_dofs_per_node() {
        let coords = vec![1.0, 0.0, 0.0];
        let u_full = vec![0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let center = [0.0, 0.0, 0.0];
        let _ = max_radial_deflection_ratio(&u_full, &coords, &axis, &center, 2);
    }
}

// ── omega_changed_significantly predicate ─────────────────────────────────────

/// Hysteresis predicate for ω-driven matrix rebuilds (K_SP or K_G).
///
/// Returns `true` when a rebuild should fire, `false` when it can be skipped.
///
/// # Arguments
/// * `omega_new`       — current angular velocity [rad/s]
/// * `omega_sq_at_last` — ω² at the last rebuild (`f64::NEG_INFINITY` on first call)
/// * `threshold_rebuild` — high-band relative threshold (triggers rebuild)
/// * `threshold_skip`    — low-band relative threshold (suppresses rebuild after a recent one)
/// * `currently_rebuilt` — `true` if a rebuild was performed recently (use `threshold_rebuild`)
/// * `eps`               — near-zero guard for ω² (use `ksp_omega_threshold²`)
///
/// # Returns
/// `true`  → rebuild the matrix
/// `false` → skip the rebuild
pub fn omega_changed_significantly(
    omega_new: f64,
    omega_sq_at_last: f64,
    threshold_rebuild: f64,
    threshold_skip: f64,
    currently_rebuilt: bool,
    eps: f64,
) -> bool {
    // First call: no prior rebuild — always rebuild.
    if omega_sq_at_last == f64::NEG_INFINITY {
        return true;
    }

    let omega_sq_new = omega_new * omega_new;
    let denom = omega_sq_new.max(omega_sq_at_last);

    // ω → 0 guard: both values near zero — rebuild to stay in a safe state.
    // eps is seeded from ksp_omega_threshold² (legacy absolute threshold squared).
    if denom < eps {
        return true;
    }

    let rel_change = (omega_sq_new - omega_sq_at_last).abs() / denom;

    // Hysteresis: use a higher bar when we recently rebuilt (prevents chattering).
    let threshold = if currently_rebuilt {
        threshold_rebuild
    } else {
        threshold_skip
    };

    rel_change >= threshold
}

// ── Tests for omega_changed_significantly and build_ksp_vals (T4.5) ───────────

#[cfg(test)]
mod tests_spin_softening {
    use super::*;

    // ── omega_changed_significantly ───────────────────────────────────────────

    /// First call (NEG_INFINITY) always triggers regardless of omega value.
    #[test]
    fn omega_pred_first_call_always_rebuilds() {
        let eps = 1e-8_f64;
        // Non-zero omega on first call.
        assert!(omega_changed_significantly(50.0, f64::NEG_INFINITY, 0.005, 0.003, false, eps));
        assert!(omega_changed_significantly(50.0, f64::NEG_INFINITY, 0.005, 0.003, true, eps));
        // Zero omega on first call.
        assert!(omega_changed_significantly(0.0, f64::NEG_INFINITY, 0.005, 0.003, false, eps));
    }

    /// Both ω values below the near-zero guard → always rebuild.
    #[test]
    fn omega_pred_near_zero_guard_always_rebuilds() {
        let eps = 1e-8_f64;
        // max(0², 0²) = 0 < eps → rebuild.
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, false, eps));
        assert!(omega_changed_significantly(0.0, 0.0, 0.005, 0.003, true, eps));
    }

    /// Change well below both thresholds → skip rebuild.
    #[test]
    fn omega_pred_stable_omega_skips() {
        let eps = 1e-8_f64;
        let omega_last_sq = 100.0f64 * 100.0f64;
        // 0.1% change in ω → ~0.2% in ω², well below both 0.3% and 0.5%.
        let omega_new = 100.0 * (1.0 + 0.001);
        assert!(!omega_changed_significantly(omega_new, omega_last_sq, 0.005, 0.003, false, eps));
        assert!(!omega_changed_significantly(omega_new, omega_last_sq, 0.005, 0.003, true, eps));
    }

    /// Large change (above both thresholds) always rebuilds.
    #[test]
    fn omega_pred_large_change_always_rebuilds() {
        let eps = 1e-8_f64;
        let omega_last_sq = 100.0f64 * 100.0f64;
        // 1% change in ω → ~2% in ω², well above both thresholds.
        let omega_new = 100.0 * (1.0 + 0.005);
        assert!(omega_changed_significantly(omega_new, omega_last_sq, 0.005, 0.003, false, eps));
        assert!(omega_changed_significantly(omega_new, omega_last_sq, 0.005, 0.003, true, eps));
    }

    /// Hysteresis: change between thresholds behaves correctly.
    #[test]
    fn omega_pred_hysteresis_between_thresholds() {
        let eps = 1e-8_f64;
        let omega_base = 100.0f64;
        let omega_sq_last = omega_base * omega_base;
        // 0.2% change in ω → ~0.4% in ω² (between 0.3% skip and 0.5% rebuild).
        let omega_new = omega_base * (1.0 + 0.002);
        let omega_sq_new = omega_new * omega_new;
        let rel = (omega_sq_new - omega_sq_last).abs() / omega_sq_new.max(omega_sq_last);
        assert!(
            rel > 0.003 && rel < 0.005,
            "Test setup: rel {rel:.4} must be in (0.003, 0.005)"
        );
        // recently_rebuilt=true → use high threshold (0.005) → rel < 0.005 → skip
        assert!(!omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, true, eps));
        // recently_rebuilt=false → use low threshold (0.003) → rel > 0.003 → rebuild
        assert!(omega_changed_significantly(omega_new, omega_sq_last, 0.005, 0.003, false, eps));
    }

    // ── build_ksp_vals ────────────────────────────────────────────────────────

    /// Zero omega → all zeros (no spin-softening at rest).
    #[test]
    fn ksp_zero_for_zero_omega() {
        let m_lumped = vec![1.0, 0.0, 0.0]; // 1 node × 3 DOFs
        let axis = [0.0f64, 0.0, 1.0];
        let omega = 0.0f64;
        let rows = vec![0i32, 1, 2];
        let cols = vec![0i32, 1, 2];
        let free_dofs = vec![0i32, 1, 2];
        let vals = build_ksp_vals(&m_lumped, &axis, omega, 3, &rows, &cols, &free_dofs, 3);
        assert!(vals.iter().all(|&v| v == 0.0), "Expected all zeros for omega=0: {vals:?}");
    }

    /// For Z-axis rotation, DOF along Z is not spin-softened (axis direction).
    #[test]
    fn ksp_z_axis_rotation_no_softening_along_axis() {
        // 1 node × 3 DOFs; mass only on DOF 0 (x-direction).
        let m_lumped = vec![1.0, 0.0, 0.0];
        let axis = [0.0f64, 0.0, 1.0]; // Z rotation axis
        let omega = 10.0f64;
        let rows = vec![0i32, 1, 2]; // diagonal only
        let cols = vec![0i32, 1, 2];
        let free_dofs = vec![0i32, 1, 2];
        let vals = build_ksp_vals(&m_lumped, &axis, omega, 3, &rows, &cols, &free_dofs, 3);
        // Z DOF (index 2) is along rotation axis → K_SP = 0.
        assert_eq!(
            vals[2], 0.0,
            "K_SP for Z DOF (along rotation axis) must be zero"
        );
        // X DOF (index 0) is perpendicular to Z → K_SP = -ω²·m·(1 - nz²) = -100·1·1 = -100.
        let expected_x = -omega * omega * m_lumped[0] * (1.0 - axis[2] * axis[2]);
        let expected_x_actual = -omega * omega * 1.0 * (1.0 - 0.0 * 0.0); // = -100
        // Note: m_lumped[0] = 1.0 applies to DOF 0. For Z axis, projector[0][0] = 1 - nz²·0 = 1.
        // Actually projector = I - n̂⊗n̂ = I - [0,0,1]⊗[0,0,1]; P[0][0] = 1 - 0 = 1.
        assert!(
            (vals[0] - expected_x_actual).abs() < 1e-10,
            "K_SP for X DOF: expected {expected_x_actual}, got {}",
            vals[0]
        );
    }

    /// K_SP scales as ω² (spec scenario R3 verification).
    #[test]
    fn ksp_scales_with_omega_squared() {
        let m_lumped = vec![2.0, 0.0, 0.0]; // mass on DOF 0
        let axis = [0.0f64, 1.0, 0.0]; // Y rotation axis
        let rows = vec![0i32];
        let cols = vec![0i32];
        let free_dofs = vec![0i32];

        let ksp_10 = build_ksp_vals(&m_lumped, &axis, 10.0, 3, &rows, &cols, &free_dofs, 3);
        let ksp_20 = build_ksp_vals(&m_lumped, &axis, 20.0, 3, &rows, &cols, &free_dofs, 3);

        if ksp_10[0] != 0.0 {
            let ratio = ksp_20[0] / ksp_10[0];
            assert!(
                (ratio - 4.0).abs() < 1e-10,
                "K_SP should scale as ω²: ratio {ratio} != 4"
            );
        }
    }

    /// K_SP matrix is symmetric (projector I - n̂⊗n̂ is symmetric).
    #[test]
    fn ksp_matrix_is_symmetric() {
        let n_nodes = 2;
        let dofs_per_node = 3;
        let n_full = n_nodes * dofs_per_node;
        let m_lumped: Vec<f64> = (0..n_full).map(|i| if i < 3 { 1.0 } else { 0.0 }).collect();
        let axis = [0.0f64, 0.0, 1.0];
        let omega = 5.0f64;
        let free_dofs: Vec<i32> = (0..n_full as i32).collect();
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        for r in 0..n_full as i32 {
            for c in 0..n_full as i32 {
                rows.push(r);
                cols.push(c);
            }
        }
        let ksp = build_ksp_vals(&m_lumped, &axis, omega, dofs_per_node, &rows, &cols, &free_dofs, n_full);
        // Verify K[i,j] == K[j,i]
        let n = n_full;
        for i in 0..n {
            for j in 0..n {
                let kij = ksp[i * n + j];
                let kji = ksp[j * n + i];
                assert!(
                    (kij - kji).abs() < 1e-14,
                    "K_SP[{i},{j}]={kij} != K_SP[{j},{i}]={kji}: not symmetric"
                );
            }
        }
    }
}
