// MITC3+ Shell Element Kernel
//
// High-performance implementation of the MITC3+ triangular shell element
// (Lee, Lee & Bathe, 2014).
//
// All geometry-constant quantities are precomputed once per element.
// No heap allocations in the hot path — fixed-size nalgebra matrices throughout.

use nalgebra::{Matrix2, Matrix3, SMatrix, SVector, Vector2, Vector3};

use crate::materials::ShellConstitutive;

// ============================================================================
// Type aliases for fixed-size element matrices
// ============================================================================

/// 18×18 element stiffness/mass matrix (3 nodes × 6 DOFs)
pub type Mat18 = SMatrix<f64, 18, 18>;
/// 20×20 extended matrix (18 nodal + 2 bubble DOFs)
pub type Mat20 = SMatrix<f64, 20, 20>;
/// 18-element force/displacement vector
pub type Vec18 = SVector<f64, 18>;
/// 20-element extended vector
pub type Vec20 = SVector<f64, 20>;

// ============================================================================
// Gauss quadrature (Hammer points, degree 2 on triangle)
// ============================================================================
// N_GAUSS, GAUSS_R, GAUSS_S, GAUSS_W are the Hammer degree-2 triangle rule
// (points 1/6,2/3 mixed, weight 1/3 each); kept local here for zero-overhead
// reference in the hot path.

pub const N_GAUSS: usize = 3;

const GAUSS_R: [f64; N_GAUSS] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
const GAUSS_S: [f64; N_GAUSS] = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
const GAUSS_W: [f64; N_GAUSS] = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

// ============================================================================
// Tying points for MITC3+ assumed shear strain
// ============================================================================

const D_PARAM: f64 = 1.0e-4;

// Constant part: A, B, C
const TP_A: (f64, f64) = (1.0 / 6.0, 2.0 / 3.0);
const TP_B: (f64, f64) = (2.0 / 3.0, 1.0 / 6.0);
const TP_C: (f64, f64) = (1.0 / 6.0, 1.0 / 6.0);

// Linear part: D, E, F
const TP_D: (f64, f64) = (1.0 / 3.0 + D_PARAM, 1.0 / 3.0 - 2.0 * D_PARAM);
const TP_E: (f64, f64) = (1.0 / 3.0 - 2.0 * D_PARAM, 1.0 / 3.0 + D_PARAM);
const TP_F: (f64, f64) = (1.0 / 3.0 + D_PARAM, 1.0 / 3.0 + D_PARAM);

// ============================================================================
// Precomputed element data (geometry-constant, computed once in `new`)
// ============================================================================

/// All data that is constant for a given element geometry.
/// Computed once, reused for every assembly call.
#[derive(Clone)]
pub struct Mitc3Precomputed {
    /// Element area
    pub area: f64,
    /// Jacobian matrix (constant for linear triangle)
    pub j_mat: Matrix2<f64>,
    /// Jacobian inverse
    pub j_inv: Matrix2<f64>,
    /// Jacobian determinant
    pub det_j: f64,
    /// Shape function derivatives in Cartesian coords: dH[0,i] = dhi/dx, dH[1,i] = dhi/dy
    pub dh: Matrix2x3,
    /// Local coordinates (x1,y1, x2,y2, x3,y3)
    pub local_coords: [f64; 6],
    /// Transformation matrix local<->global (3×3 rotation part)
    pub t3: Matrix3<f64>,
    /// Covariant base vectors (constant for linear triangle)
    pub g_r: Vector2<f64>,
    pub g_s: Vector2<f64>,
    /// Constitutive matrices (material + thickness dependent, but constant per element)
    pub constitutive: ShellConstitutive,
    /// Drilling stiffness factor
    pub k_drill: f64,
    pub drilling_scale: f64,
    /// Thickness
    pub thickness: f64,
    /// Precomputed covariant shear at the 6 tying points (extended 20-DOF space)
    /// Each is (B_ert, B_est) as Vec20
    pub tying_ext_a: (Vec20, Vec20),
    pub tying_ext_b: (Vec20, Vec20),
    pub tying_ext_c: (Vec20, Vec20),
    pub tying_ext_d: (Vec20, Vec20),
    pub tying_ext_e: (Vec20, Vec20),
    pub tying_ext_f: (Vec20, Vec20),
}

/// 2×3 matrix type for shape function derivatives
type Matrix2x3 = SMatrix<f64, 2, 3>;

// ============================================================================
// Shape functions (pure, no state)
// ============================================================================

/// Linear shape functions: h1 = 1-r-s, h2 = r, h3 = s
#[inline(always)]
fn shape_functions(r: f64, s: f64) -> Vector3<f64> {
    Vector3::new(1.0 - r - s, r, s)
}

/// Bubble function: f4 = 27·r·s·(1-r-s)
#[inline(always)]
fn bubble(r: f64, s: f64) -> f64 {
    27.0 * r * s * (1.0 - r - s)
}

/// Bubble derivatives: (df4/dr, df4/ds)
#[inline(always)]
fn bubble_deriv(r: f64, s: f64) -> (f64, f64) {
    (
        27.0 * s * (1.0 - 2.0 * r - s),
        27.0 * r * (1.0 - r - 2.0 * s),
    )
}

/// Enriched shape functions: [f1, f2, f3, f4]
/// fi = hi - f4/3 for i=1,2,3
#[inline(always)]
fn enriched_shape(r: f64, s: f64) -> [f64; 4] {
    let h = shape_functions(r, s);
    let f4 = bubble(r, s);
    let f4_3 = f4 / 3.0;
    [h[0] - f4_3, h[1] - f4_3, h[2] - f4_3, f4]
}

/// Enriched shape function derivatives in natural coords
/// Returns (df/dr[4], df/ds[4])
#[inline(always)]
fn enriched_shape_deriv(r: f64, s: f64) -> ([f64; 4], [f64; 4]) {
    let (df4_dr, df4_ds) = bubble_deriv(r, s);
    let df4_dr_3 = df4_dr / 3.0;
    let df4_ds_3 = df4_ds / 3.0;
    (
        [-1.0 - df4_dr_3, 1.0 - df4_dr_3, -df4_dr_3, df4_dr],
        [-1.0 - df4_ds_3, -df4_ds_3, 1.0 - df4_ds_3, df4_ds],
    )
}

// ============================================================================
// Material frame rotation
// ============================================================================

/// Rotate the shell constitutive matrices from the global material reference
/// frame to the element local frame.
///
/// The global material reference direction is taken as global X = (1,0,0)
/// projected onto the shell plane. When the element's local e1 differs from
/// this reference (e.g. diagonal triangles formed by splitting quads), the
/// ABD matrices must be rotated to avoid a frame mismatch.
///
/// Rotation by angle β (angle of e1 from global X in the shell plane):
///   - cos_b = e1 · x_ref,  sin_b = e1 · y_ref  (y_ref = e3 × x_ref)
///   - Reuter T matrix (angle = -β): c = cos_b, s = -sin_b
///   - A/B/D_local = T(-β) · A/B/D_global · T(-β)^T
///   - Cs_local = R · Cs_global · R^T  (R is 2×2 in-plane rotation)
fn rotate_constitutive_to_local(
    constitutive: ShellConstitutive,
    e1: &Vector3<f64>,
    e3: &Vector3<f64>,
) -> ShellConstitutive {
    // Project global X onto the shell plane
    let x_global = Vector3::new(1.0, 0.0, 0.0);
    let x_proj = x_global - x_global.dot(e3) * e3;
    let x_ref_len = x_proj.norm();

    // If global X is nearly normal to shell, use global Z as reference instead
    let x_ref = if x_ref_len > 1e-10 {
        x_proj / x_ref_len
    } else {
        let z_global = Vector3::new(0.0, 0.0, 1.0);
        let z_proj = z_global - z_global.dot(e3) * e3;
        let z_len = z_proj.norm();
        if z_len > 1e-10 { z_proj / z_len } else { return constitutive; }
    };

    // cos and sin of angle β (angle of e1 from x_ref, measured in shell plane)
    let cos_b = e1.dot(&x_ref).clamp(-1.0, 1.0);
    let y_ref = e3.cross(&x_ref); // second in-plane reference direction
    let sin_b = e1.dot(&y_ref).clamp(-1.0, 1.0);

    // If already aligned (cos_b ≈ 1, sin_b ≈ 0), skip rotation
    if (cos_b - 1.0).abs() < 1e-12 {
        return constitutive;
    }

    // Reuter transformation T for rotation by (+β):
    // C_local = T(β) · C_global · T^T(β)
    // where β = angle of element e1 from global reference x_ref
    let c = cos_b;
    let s = sin_b;

    // T = [[c²,    s²,    2cs  ],
    //      [s²,    c²,   -2cs  ],
    //      [-cs,   cs,  c²-s² ]]
    let t = Matrix3::new(
        c * c,           s * s,           2.0 * c * s,
        s * s,           c * c,          -2.0 * c * s,
        -c * s,          c * s,           c * c - s * s,
    );

    // Rotate in-plane constitutive matrices: M_local = T · M_global · T^T
    let cm_rot           = t * constitutive.cm           * t.transpose();
    let cb_rot           = t * constitutive.cb           * t.transpose();
    let cb_coupling_rot  = t * constitutive.cb_coupling  * t.transpose();
    let cm_raw_rot       = t * constitutive.cm_raw       * t.transpose();

    // Rotate transverse shear (2×2): Cs_local = R · Cs_global · R^T
    // where R = [[cos_b, sin_b], [-sin_b, cos_b]] rotates (x_ref,y_ref) → (e1,e2)
    let r = Matrix2::new(cos_b, sin_b, -sin_b, cos_b);
    let cs_rot = r * constitutive.cs * r.transpose();

    ShellConstitutive {
        cm:          cm_rot,
        cb:          cb_rot,
        cb_coupling: cb_coupling_rot,
        cs:          cs_rot,
        cm_raw:      cm_raw_rot,
    }
}

// ============================================================================
// Precomputation
// ============================================================================

impl Mitc3Precomputed {
    /// Create precomputed data for a MITC3+ element.
    ///
    /// `coords` is [x1,y1,z1, x2,y2,z2, x3,y3,z3] (9 values, row-major 3×3)
    /// `material_constitutive` is the precomputed constitutive matrices.
    /// `thickness` is the shell thickness.
    /// `shear_correction` is the shear correction factor (default 5/6).
    pub fn new(
        coords: &[f64; 9],
        constitutive: ShellConstitutive,
        thickness: f64,
        e_modulus: f64,
        drilling_scale: f64,
    ) -> Self {
        let p1 = Vector3::new(coords[0], coords[1], coords[2]);
        let p2 = Vector3::new(coords[3], coords[4], coords[5]);
        let p3 = Vector3::new(coords[6], coords[7], coords[8]);

        // Local coordinate system
        let v12 = p2 - p1;
        let v13 = p3 - p1;
        let l12 = v12.norm();
        let e1 = v12 / l12;
        let n = v12.cross(&v13);
        let n_len = n.norm();
        let e3 = n / n_len;
        let e2 = e3.cross(&e1);

        // Rotation matrix (rows = local basis vectors)
        let t3 = Matrix3::new(
            e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2],
        );

        // Local 2D coordinates
        let x1 = 0.0;
        let y1 = 0.0;
        let x2 = l12;
        let y2 = 0.0;
        let x3 = v13.dot(&e1);
        let y3 = v13.dot(&e2);
        let local_coords = [x1, y1, x2, y2, x3, y3];

        // Jacobian (constant for linear triangle)
        let dx_dr = x2 - x1;
        let dy_dr = y2 - y1;
        let dx_ds = x3 - x1;
        let dy_ds = y3 - y1;
        let j_mat = Matrix2::new(dx_dr, dy_dr, dx_ds, dy_ds);
        let det_j = dx_dr * dy_ds - dy_dr * dx_ds;
        let inv_det = 1.0 / det_j;
        let j_inv = Matrix2::new(
            dy_ds * inv_det,
            -dy_dr * inv_det,
            -dx_ds * inv_det,
            dx_dr * inv_det,
        );

        // Area = |detJ| / 2
        let area = det_j.abs() / 2.0;

        // Shape function derivatives in Cartesian (constant)
        // dh/dr = [-1, 1, 0], dh/ds = [-1, 0, 1]
        let dh_dr = Vector3::new(-1.0, 1.0, 0.0);
        let dh_ds = Vector3::new(-1.0, 0.0, 1.0);
        let dh_dx = j_inv[(0, 0)] * dh_dr + j_inv[(0, 1)] * dh_ds;
        let dh_dy = j_inv[(1, 0)] * dh_dr + j_inv[(1, 1)] * dh_ds;
        let dh = Matrix2x3::from_rows(&[dh_dx.transpose(), dh_dy.transpose()]);

        // Covariant base vectors (2D, in local plane)
        let g_r = Vector2::new(x2 - x1, y2 - y1);
        let g_s = Vector2::new(x3 - x1, y3 - y1);

        // Rotate constitutive matrices from the global material reference frame
        // (direction 1 = global X projected onto shell plane) to the element
        // local frame (direction 1 = e1).
        //
        // This is necessary because the laminate ABD matrices are computed with
        // a fixed global reference direction (global X), but the MITC3 element
        // computes strains in its own local frame (e1 = edge 0→1 direction).
        // For elements where e1 ≠ global-X (e.g., diagonal triangles from quad
        // splitting), the material must be rotated to avoid a frame mismatch.
        let constitutive = rotate_constitutive_to_local(constitutive, &e1, &e3);

        // Drilling stiffness. No source is established for the 0.15 factor;
        // it is the same uncited penalty as in mitc4.rs.
        let k_drill = e_modulus * thickness * thickness * 0.15 * drilling_scale;

        // Precompute tying point shear evaluations (extended space)
        let tying_ext_a =
            eval_covariant_shear_ext(TP_A.0, TP_A.1, &g_r, &g_s);
        let tying_ext_b =
            eval_covariant_shear_ext(TP_B.0, TP_B.1, &g_r, &g_s);
        let tying_ext_c =
            eval_covariant_shear_ext(TP_C.0, TP_C.1, &g_r, &g_s);
        let tying_ext_d =
            eval_covariant_shear_ext(TP_D.0, TP_D.1, &g_r, &g_s);
        let tying_ext_e =
            eval_covariant_shear_ext(TP_E.0, TP_E.1, &g_r, &g_s);
        let tying_ext_f =
            eval_covariant_shear_ext(TP_F.0, TP_F.1, &g_r, &g_s);

        Self {
            area,
            j_mat,
            j_inv,
            det_j,
            dh,
            local_coords,
            t3,
            g_r,
            g_s,
            constitutive,
            k_drill,
            drilling_scale,
            thickness,
            tying_ext_a,
            tying_ext_b,
            tying_ext_c,
            tying_ext_d,
            tying_ext_e,
            tying_ext_f,
        }
    }
}

// ============================================================================
// Covariant shear evaluation at a tying point (extended 20-DOF space)
// ============================================================================

/// Evaluate covariant transverse shear B-vectors (e_rt, e_st) at a parametric
/// point in the extended 20-DOF space (18 nodal + 2 bubble rotations).
///
/// Returns (B_ert, B_est) each as Vec20.
fn eval_covariant_shear_ext(
    r: f64,
    s: f64,
    g_r: &Vector2<f64>,
    g_s: &Vector2<f64>,
) -> (Vec20, Vec20) {
    let f = enriched_shape(r, s);

    let mut b_ert = Vec20::zeros();
    let mut b_est = Vec20::zeros();

    // Corner nodes
    for i in 0..3 {
        let w_idx = 6 * i + 2;
        let thx_idx = 6 * i + 3;
        let thy_idx = 6 * i + 4;

        // dh/dr = [-1, 1, 0], dh/ds = [-1, 0, 1]
        let dhi_dr = [-1.0, 1.0, 0.0][i];
        let dhi_ds = [-1.0, 0.0, 1.0][i];

        // w uses linear interpolation
        b_ert[w_idx] = dhi_dr;
        b_est[w_idx] = dhi_ds;

        // Rotation contribution: physical Reissner–Mindlin director.
        //
        // From Ko, Bathe & Zhang 2025 (MITC4/D), Eq. (3a), the offset
        // displacement interpolates as
        //     u(r,s,t) = Σ h_i u_i + (t/2) Σ a_i h_i (θ_i × V_in),
        // where θ is the physical rotation vector and V_in the shell director.
        // Linearizing the director about the flat reference gives
        //     V3 = e3 + θ × e3 = e3 + θy·e1 − θx·e2,
        // hence
        //     e_rt = dw/dr + V3·g_r = dw/dr + θy·g_r[0] − θx·g_r[1].
        // This is the Cartesian convention γ_xz = w_,x + θy and
        // γ_yz = w_,y − θx, the same one b_gamma_mitc4 uses.
        b_ert[thy_idx] = f[i] * g_r[0];
        b_ert[thx_idx] = -f[i] * g_r[1];

        b_est[thy_idx] = f[i] * g_s[0];
        b_est[thx_idx] = -f[i] * g_s[1];
    }

    // Bubble rotations (indices 18, 19)
    b_ert[19] = f[3] * g_r[0]; // thy4
    b_ert[18] = -f[3] * g_r[1]; // thx4

    b_est[19] = f[3] * g_s[0];
    b_est[18] = -f[3] * g_s[1];

    (b_ert, b_est)
}

// ============================================================================
// Extended bending B-matrix (3×20)
// ============================================================================

type Mat3x20 = SMatrix<f64, 3, 20>;

fn b_kappa_ext(r: f64, s: f64, j_inv: &Matrix2<f64>) -> Mat3x20 {
    let (df_dr, df_ds) = enriched_shape_deriv(r, s);

    // Convert to Cartesian
    let mut df_dx = [0.0; 4];
    let mut df_dy = [0.0; 4];
    for i in 0..4 {
        df_dx[i] = j_inv[(0, 0)] * df_dr[i] + j_inv[(0, 1)] * df_ds[i];
        df_dy[i] = j_inv[(1, 0)] * df_dr[i] + j_inv[(1, 1)] * df_ds[i];
    }

    let mut bk = Mat3x20::zeros();

    // Corner nodes
    for i in 0..3 {
        let thx = 6 * i + 3;
        let thy = 6 * i + 4;

        // κxx = ∂θy/∂x
        bk[(0, thy)] = df_dx[i];
        // κyy = -∂θx/∂y
        bk[(1, thx)] = -df_dy[i];
        // κxy = ∂θy/∂y - ∂θx/∂x
        bk[(2, thy)] = df_dy[i];
        bk[(2, thx)] = -df_dx[i];
    }

    // Bubble (indices 18=thx4, 19=thy4)
    bk[(0, 19)] = df_dx[3];
    bk[(1, 18)] = -df_dy[3];
    bk[(2, 19)] = df_dy[3];
    bk[(2, 18)] = -df_dx[3];

    bk
}

// ============================================================================
// Extended shear B-matrix (2×20) using precomputed tying point data
// ============================================================================

type Mat2x20 = SMatrix<f64, 2, 20>;

fn b_gamma_ext(
    r: f64,
    s: f64,
    pre: &Mitc3Precomputed,
) -> Mat2x20 {
    let (ref b_ert_a, ref b_est_a) = pre.tying_ext_a;
    let (ref b_ert_b, ref b_est_b) = pre.tying_ext_b;
    let (ref b_ert_c, ref b_est_c) = pre.tying_ext_c;
    let (ref b_ert_d, ref _b_est_d) = pre.tying_ext_d;
    let (ref _b_ert_e, ref b_est_e) = pre.tying_ext_e;
    let (ref b_ert_f, ref b_est_f) = pre.tying_ext_f;

    // Constant part (Eq. 15)
    let b_ert_const =
        (2.0 / 3.0) * (b_ert_b - 0.5 * b_est_b) + (1.0 / 3.0) * (b_ert_c + b_est_c);
    let b_est_const =
        (2.0 / 3.0) * (b_est_a - 0.5 * b_ert_a) + (1.0 / 3.0) * (b_ert_c + b_est_c);

    // Linear part (Eq. 16)
    let b_c_hat = (b_ert_f - b_ert_d) - (b_est_f - b_est_e);
    let b_ert_linear = (1.0 / 3.0) * &b_c_hat * (3.0 * s - 1.0);
    let b_est_linear = (1.0 / 3.0) * &b_c_hat * (1.0 - 3.0 * r);

    // Total assumed
    let b_ert = b_ert_const + b_ert_linear;
    let b_est = b_est_const + b_est_linear;

    // Transform covariant → Cartesian shear: γ = J^{-1} · e_cov
    // Derivation: e_cov = J · γ_cart (e_rt = g_r·γ, e_st = g_s·γ)
    // so γ_xz = J^{-1}[0,·]·e_cov,  γ_yz = J^{-1}[1,·]·e_cov
    let j_inv = &pre.j_inv;
    let b_xz = j_inv[(0, 0)] * &b_ert + j_inv[(0, 1)] * &b_est;
    let b_yz = j_inv[(1, 0)] * &b_ert + j_inv[(1, 1)] * &b_est;

    let mut bg = Mat2x20::zeros();
    for j in 0..20 {
        bg[(0, j)] = b_xz[j];
        bg[(1, j)] = b_yz[j];
    }
    bg
}

// ============================================================================
// Membrane B-matrix (3×18) — standard CST
// ============================================================================

pub type Mat3x18 = SMatrix<f64, 3, 18>;

/// The union DOF count of a strain-smoothed element: six nodes times six DOFs.
pub const SMOOTHED_UNION_DOFS: usize = 6 * SMOOTHED_UNION_NODES;

/// The membrane strain operator of a smoothed element, over the union layout.
pub type Mat3Union = SMatrix<f64, 3, { SMOOTHED_UNION_DOFS }>;

/// The element stiffness of a smoothed element, over the union layout.
pub type MatUnion = SMatrix<f64, { SMOOTHED_UNION_DOFS }, { SMOOTHED_UNION_DOFS }>;

/// The union displacement vector of a smoothed element.
pub type VecUnion = SMatrix<f64, { SMOOTHED_UNION_DOFS }, 1>;

fn b_membrane(dh: &Matrix2x3) -> Mat3x18 {
    let mut bm = Mat3x18::zeros();

    for i in 0..3 {
        let u_idx = 6 * i;
        let v_idx = 6 * i + 1;
        let dhi_dx = dh[(0, i)];
        let dhi_dy = dh[(1, i)];

        // εxx = ∂u/∂x
        bm[(0, u_idx)] = dhi_dx;
        // εyy = ∂v/∂y
        bm[(1, v_idx)] = dhi_dy;
        // γxy = ∂u/∂y + ∂v/∂x
        bm[(2, u_idx)] = dhi_dy;
        bm[(2, v_idx)] = dhi_dx;
    }

    bm
}

// ============================================================================
// Drilling B-matrix (row for θz only)
// ============================================================================

type RowVec18 = SMatrix<f64, 1, 18>;

/// Drilling B-row (1×18): gamma = (dv/dx - du/dy)/2 - theta_z.
/// No source is cited here; `mitc4.rs::b_drill` applies the same operator with an
/// unverified "Hughes & Brezzi" attribution (see docs/references.md §1).
fn b_drill(r: f64, s: f64, dh: &Matrix2x3) -> RowVec18 {
    let h = shape_functions(r, s);
    let mut bd = RowVec18::zeros();

    for i in 0..3 {
        let u_idx = 6 * i;
        let v_idx = 6 * i + 1;
        let thz_idx = 6 * i + 5;
        let dhi_dx = dh[(0, i)];
        let dhi_dy = dh[(1, i)];

        // (∂v/∂x - ∂u/∂y)/2 - θz
        bd[(0, u_idx)] = -0.5 * dhi_dy;
        bd[(0, v_idx)] = 0.5 * dhi_dx;
        bd[(0, thz_idx)] = -h[i];
    }

    bd
}

// ============================================================================
// Stiffness matrix computation
// ============================================================================

/// Compute the 18×18 element stiffness matrix in LOCAL coordinates.
///
/// This computes K = K_membrane + K_bending_shear_condensed + K_coupling.
/// The bending-shear part uses the MITC3+ bubble enrichment and is statically
/// condensed from a 20-DOF to 18-DOF space.
pub fn compute_ke_local(pre: &Mitc3Precomputed) -> Mat18 {
    let bm = b_membrane(&pre.dh);
    compute_ke_local_with_membrane(pre, &[bm; N_GAUSS])
}

/// Compute the 18×18 element stiffness in LOCAL coordinates using a prescribed
/// membrane strain operator at each Gauss point.
///
/// `bm_gp[gp]` is the local Cartesian membrane B matrix used at Gauss point
/// `gp`.  Passing the element's own `b_membrane` at all three points reproduces
/// [`compute_ke_local`] exactly, which is what keeps the MITC3+ of Lee, Lee &
/// Bathe 2014 available.  Passing the strain-smoothed operators built by
/// [`smoothed_membrane_b`] gives the strain-smoothed MITC3+ of Lee & Lee 2019.
///
/// Only the membrane field is parameterised.  The bending and transverse shear
/// fields below are the MITC3+ ones either way, which is what the 2019 paper
/// prescribes: *"We use the originally defined b1 eij and b2 eij for the
/// covariant bending strains. For the covariant transverse shear strains, we
/// adopt the assumed strains of the MITC3+ shell element."*
pub fn compute_ke_local_with_membrane(pre: &Mitc3Precomputed, bm_gp: &[Mat3x18; N_GAUSS]) -> Mat18 {
    let area = pre.area;
    let cm = &pre.constitutive.cm;
    let cb = &pre.constitutive.cb;
    let cb_coupling = &pre.constitutive.cb_coupling;
    let cs = &pre.constitutive.cs;
    let dh = &pre.dh;
    let j_inv = &pre.j_inv;

    // --- Drilling stiffness ---
    let mut k_drill_total = Mat18::zeros();
    for gp in 0..N_GAUSS {
        let bd = b_drill(GAUSS_R[gp], GAUSS_S[gp], dh);
        k_drill_total += (GAUSS_W[gp] * area * pre.k_drill) * (bd.transpose() * &bd);
    }

    // --- Membrane, bending, shear and B-coupling over the Gauss points ---
    //
    // The membrane term is accumulated inside this loop rather than before it,
    // because a smoothed membrane operator differs from Gauss point to Gauss
    // point.  The weights sum to one, so a constant operator reproduces the
    // previous `area * bm^T cm bm` exactly.
    let mut km = Mat18::zeros();
    let mut k_ext = Mat20::zeros();
    // B-coupling: bm^T · B · bk_ext  (18×20) — full extended space
    let mut k_mb_ext = SMatrix::<f64, 18, 20>::zeros();
    for gp in 0..N_GAUSS {
        let r = GAUSS_R[gp];
        let s = GAUSS_S[gp];
        let w = GAUSS_W[gp];
        let bm = &bm_gp[gp];

        let bk = b_kappa_ext(r, s, j_inv);
        let bg = b_gamma_ext(r, s, pre);

        km += (w * area) * (bm.transpose() * cm * bm);
        k_ext += w * area * (bk.transpose() * cb * &bk + bg.transpose() * cs * &bg);

        // Membrane-bending coupling
        k_mb_ext += (w * area) * (bm.transpose() * cb_coupling * &bk);
    }

    // Extract nodal and bubble blocks
    let k_uu = k_ext.fixed_view::<18, 18>(0, 0).into_owned();
    let k_uq = k_ext.fixed_view::<18, 2>(0, 18).into_owned();
    let k_qu = k_ext.fixed_view::<2, 18>(18, 0).into_owned();
    let k_qq = k_ext.fixed_view::<2, 2>(18, 18).into_owned();

    // B-coupling blocks
    let k_mb_uu: SMatrix<f64, 18, 18> = k_mb_ext.fixed_view::<18, 18>(0, 0).into_owned();
    let k_mb_uq: SMatrix<f64, 18, 2>  = k_mb_ext.fixed_view::<18, 2>(0, 18).into_owned();

    // 2×2 inverse (explicit for performance)
    let det_qq = k_qq[(0, 0)] * k_qq[(1, 1)] - k_qq[(0, 1)] * k_qq[(1, 0)];
    let inv_qq = SMatrix::<f64, 2, 2>::new(
        k_qq[(1, 1)] / det_qq,
        -k_qq[(0, 1)] / det_qq,
        -k_qq[(1, 0)] / det_qq,
        k_qq[(0, 0)] / det_qq,
    );

    // Modified [uq] block: include B-coupling to bubble DOFs
    let k_uq_full = k_uq + &k_mb_uq;
    let k_qu_full = k_qu + k_mb_uq.transpose();

    // Condensed bending+shear+B-coupling stiffness
    let k_mb_sym = &k_mb_uu + k_mb_uu.transpose();
    let k_bs_cond = k_uu + k_mb_sym - &k_uq_full * &inv_qq * &k_qu_full;

    // Total local stiffness
    let k_local = km + k_drill_total + k_bs_cond;

    // Symmetrize
    0.5 * (&k_local + k_local.transpose())
}

/// The element stiffness over the union DOF layout of a strain-smoothed element.
///
/// Rows and columns `0..18` are the target's own DOFs; `18..36` are the unique
/// node of each of the three edge neighbours, in edge order.  The bending,
/// transverse shear and drilling terms act only on the target's own DOFs - the
/// 2019 paper keeps them as they are - while the membrane and the
/// membrane-bending coupling use the supplied union membrane operators, which is
/// what couples the target to its neighbours.
///
/// The result is in the union's **local** DOF space; rotate it with
/// [`transform_union_to_global`] before assembling.  With a membrane operator
/// whose first 18 columns are the element's own `b_membrane` and whose remaining
/// columns are zero, the top-left 18x18 block reproduces [`compute_ke_local`]
/// exactly.  That is the property the unit tests pin down - for an isotropic
/// material and again with the membrane-bending coupling switched on, because
/// the isotropic case alone cannot see the coupling - and it is what makes this a
/// strict extension of the MITC3+ of Lee, Lee & Bathe 2014 rather than a
/// replacement.
///
/// **Known defect (issue #2, not localized).** With the rotational block of
/// [`union_rotation`] fixed, the union path still over-stiffens a curved shell:
/// Scordelis-Lo on a triangular mesh gives 0.0145 (N=8) against the un-smoothed
/// MITC3+ 0.8561, where Lee & Lee 2019 (CAS 223:106096, Table 6) reports the
/// smoothing *improving* the un-smoothed element.  The cause is in the
/// curved-shell covariant handling or the union assembly, and the kernel stays
/// unwired (nothing consumes it) until it is found.
pub fn compute_ke_union_smoothed(pre: &Mitc3Precomputed, bm_union: &[Mat3Union; N_GAUSS]) -> MatUnion {
    const U: usize = SMOOTHED_UNION_DOFS;
    let area = pre.area;
    let cm = &pre.constitutive.cm;
    let cb = &pre.constitutive.cb;
    let cb_coupling = &pre.constitutive.cb_coupling;
    let cs = &pre.constitutive.cs;
    let dh = &pre.dh;
    let j_inv = &pre.j_inv;

    // --- Drilling stiffness: the target's own DOFs only ---
    let mut k_drill = Mat18::zeros();
    for gp in 0..N_GAUSS {
        let bd = b_drill(GAUSS_R[gp], GAUSS_S[gp], dh);
        k_drill += (GAUSS_W[gp] * area * pre.k_drill) * (bd.transpose() * &bd);
    }
    let mut k_drill_union = MatUnion::zeros();
    k_drill_union
        .fixed_view_mut::<18, 18>(0, 0)
        .copy_from(&k_drill);

    // --- Membrane and coupling over the union; bending and shear on the element ---
    let mut km = MatUnion::zeros();
    let mut k_ext = Mat20::zeros();
    let mut k_mb_ext = SMatrix::<f64, U, 20>::zeros();
    for gp in 0..N_GAUSS {
        let r = GAUSS_R[gp];
        let s = GAUSS_S[gp];
        let w = GAUSS_W[gp];
        let bm = &bm_union[gp];

        let bk = b_kappa_ext(r, s, j_inv);
        let bg = b_gamma_ext(r, s, pre);

        km += (w * area) * (bm.transpose() * cm * bm);
        k_ext += w * area * (bk.transpose() * cb * &bk + bg.transpose() * cs * &bg);
        k_mb_ext += (w * area) * (bm.transpose() * cb_coupling * &bk);
    }

    // The element-only bending and shear block lives in the target's own DOFs.
    let mut k_uu = MatUnion::zeros();
    k_uu.fixed_view_mut::<18, 18>(0, 0)
        .copy_from(&k_ext.fixed_view::<18, 18>(0, 0));
    let k_uq = k_ext.fixed_view::<18, 2>(0, 18).into_owned();
    let k_qu = k_ext.fixed_view::<2, 18>(18, 0).into_owned();
    let k_qq = k_ext.fixed_view::<2, 2>(18, 18).into_owned();

    let k_mb_uu: SMatrix<f64, U, 18> = k_mb_ext.fixed_view::<U, 18>(0, 0).into_owned();
    let k_mb_uq: SMatrix<f64, U, 2> = k_mb_ext.fixed_view::<U, 2>(0, 18).into_owned();

    // 2x2 inverse of the bubble block (explicit, as in the un-smoothed path).
    let det_qq = k_qq[(0, 0)] * k_qq[(1, 1)] - k_qq[(0, 1)] * k_qq[(1, 0)];
    let inv_qq = SMatrix::<f64, 2, 2>::new(
        k_qq[(1, 1)] / det_qq,
        -k_qq[(0, 1)] / det_qq,
        -k_qq[(1, 0)] / det_qq,
        k_qq[(0, 0)] / det_qq,
    );

    // The target's own u-q block, then the coupling's.
    let mut k_uq_union = SMatrix::<f64, U, 2>::zeros();
    k_uq_union.fixed_view_mut::<18, 2>(0, 0).copy_from(&k_uq);
    let k_uq_full = k_uq_union + &k_mb_uq;
    let k_qu_full = k_uq_full.transpose();

    // The coupling's u-u contribution.  In the element's own space it is
    // `k_mb_uu + k_mb_uu^T`, the two cross terms of the coupling's virtual work.
    // In the union space those two terms are transposes of each other and land in
    // different sub-blocks of the same 36x36 matrix: the first spans all union
    // rows against the target's own columns, the second the target's own rows
    // against all union columns.
    let mut coupling_rows = MatUnion::zeros();
    coupling_rows
        .fixed_view_mut::<U, 18>(0, 0)
        .copy_from(&k_mb_uu);
    let mut coupling_columns = MatUnion::zeros();
    coupling_columns
        .fixed_view_mut::<18, U>(0, 0)
        .copy_from(&k_mb_uu.transpose());
    let k_mb_sym = coupling_rows + coupling_columns;

    let k_bs_cond = k_uu + k_mb_sym - &k_uq_full * &inv_qq * &k_qu_full;
    let k_total = km + k_drill_union + k_bs_cond;

    0.5 * (&k_total + k_total.transpose())
}


/// The internal force of a strain-smoothed element, over the union layout.
///
/// `u_union` holds the union's **local** DOF components.  The linear part is
/// `K_union · u_union`, with `K_union` the stiffness of
/// [`compute_ke_union_smoothed`], so at zero displacement this agrees with it
/// exactly - which is what keeps `K` and `K_T` consistent.
///
/// The geometric (Green-Lagrange) correction is evaluated on the target's own
/// DOFs and added to the target's own rows: the 2019 paper presents the smoothed
/// element for nonlinear analysis too, but how the smoothing interacts with the
/// geometric terms is not extracted from it yet, so those stay on the element's
/// own field.  That is a documented limitation, not an accident.
pub fn compute_fint_union_smoothed(
    pre: &Mitc3Precomputed,
    u_union: &VecUnion,
    bm_union: &[Mat3Union; N_GAUSS],
    nonlinear: bool,
) -> VecUnion {
    let mut f = compute_ke_union_smoothed(pre, bm_union) * u_union;

    if nonlinear {
        // The union vector is already in local components, so the target's own
        // rows are its local DOFs and no rotation is needed here.
        let u_local = u_union.fixed_rows::<18>(0).into_owned();
        let h_mat = displacement_gradient(pre, &u_local);
        let eps_gl = gl_strain_voigt(&h_mat);
        let sigma_m = &pre.constitutive.cm_raw * &eps_gl;

        let b_l = compute_b_l(&pre.dh);
        let b_nl = compute_b_nl(&pre.dh, &h_mat);
        let b_total = b_l + b_nl;

        let f_nonlinear = (pre.area * pre.thickness) * (b_total.transpose() * &sigma_m);
        let mut target_rows = f.fixed_rows::<18>(0).into_owned();
        target_rows += f_nonlinear;
        f.fixed_rows_mut::<18>(0).copy_from(&target_rows);
    }

    f
}

/// The tangent stiffness of a strain-smoothed element, over the union layout.
///
/// `u_union` holds the union's **local** DOF components.  The linear part is
/// exactly [`compute_ke_union_smoothed`], so at zero displacement this equals the
/// assembled `K` - the consistency the assembler's `K`/`K_T` test checks.
///
/// The geometric (K_sigma) and initial-displacement (K_L) terms are evaluated on
/// the target's own DOFs and added to its own rows, for the reason recorded on
/// [`compute_fint_union_smoothed`].
pub fn compute_kt_union_smoothed(
    pre: &Mitc3Precomputed,
    u_union: &VecUnion,
    bm_union: &[Mat3Union; N_GAUSS],
) -> MatUnion {
    let k0 = compute_ke_union_smoothed(pre, bm_union);

    // The union vector is already in local components.
    let u_local = u_union.fixed_rows::<18>(0).into_owned();

    let h_mat = displacement_gradient(pre, &u_local);
    let eps_gl = gl_strain_voigt(&h_mat);
    let sigma_m = &pre.constitutive.cm_raw * &eps_gl;
    let k_sigma_local = compute_k_sigma_local(pre, &sigma_m);

    let dh = &pre.dh;
    let b_l = compute_b_l(dh);
    let b_nl = compute_b_nl(dh, &h_mat);
    let cm_raw = &pre.constitutive.cm_raw;
    let k_l_local = pre.area
        * pre.thickness
        * (b_l.transpose() * cm_raw * &b_nl
            + b_nl.transpose() * cm_raw * &b_l
            + b_nl.transpose() * cm_raw * &b_nl);

    let nonlinear_local = k_l_local + k_sigma_local;

    let mut out = k0;
    let mut target_block = out.fixed_view::<18, 18>(0, 0).into_owned();
    target_block += nonlinear_local;
    out.fixed_view_mut::<18, 18>(0, 0).copy_from(&target_block);
    0.5 * (&out + out.transpose())
}


/// The union's local-to-global rotation: the target's frame for its own nodes
/// (slots 0..2) and each neighbour's own frame for its unique node (slot 3 + k).
///
/// This is what makes the union assembly frame-consistent.  The smoothed
/// operator maps each element's own local DOFs, so the finished stiffness must be
/// rotated block by block with the frame of the element that owns each slot -
/// using the target's frame for all six slots would be wrong for the neighbours.
///
/// `frames[slot]` is the owning element's `t3`; a padded boundary slot never
/// carries weight, so its frame is irrelevant.
pub fn union_rotation(frames: &[Matrix3<f64>; SMOOTHED_UNION_NODES]) -> MatUnion {
    let mut t = MatUnion::zeros();
    for (slot, frame) in frames.iter().enumerate() {
        let base = 6 * slot;
        for a in 0..3 {
            for b in 0..3 {
                // Both the translational and the rotational 3x3 slots carry the
                // same frame.  Leaving the rotational block zero (issue #2) made
                // T singular and annihilated every rotational DOF in
                // `transform_union_to_global` = T^T K T.
                t[(base + a, base + b)] = frame[(a, b)];
                t[(base + 3 + a, base + 3 + b)] = frame[(a, b)];
            }
        }
    }
    t
}

/// Rotate a union stiffness from the union's local DOFs to global ones.
pub fn transform_union_to_global(
    k_local: &MatUnion,
    frames: &[Matrix3<f64>; SMOOTHED_UNION_NODES],
) -> MatUnion {
    let t = union_rotation(frames);
    t.transpose() * k_local * t
}

/// Rotate a union force from the union's local DOFs to global ones.
pub fn union_force_to_global(
    f_local: &VecUnion,
    frames: &[Matrix3<f64>; SMOOTHED_UNION_NODES],
) -> VecUnion {
    union_rotation(frames).transpose() * f_local
}

/// Compute the 18×18 element stiffness in GLOBAL coordinates.
pub fn compute_ke_global(pre: &Mitc3Precomputed) -> Mat18 {
    let k_local = compute_ke_local(pre);
    transform_to_global(&k_local, &pre.t3)
}

/// Number of nodes in the union layout a smoothed element spans: the target's
/// three nodes plus the unique node of each of its three edge neighbours.
pub const SMOOTHED_UNION_NODES: usize = 6;

/// The strain-smoothed membrane strain operators of a target triangle, one per
/// Gauss point, expressed over the union DOF layout.
///
/// Reference: Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite
/// element", *Computers and Structures* 223:106096, 2019.  The weights come from
/// [`crate::elements::smoothing::smoothed_membrane_strain`], which also holds
/// the derivation: Eq. (15) carries each neighbour's covariant strain into the
/// target's convected coordinates, Eqs. (16) and (17) average it with the target
/// weighted by projected area, the boundary rule falls back to the target's own
/// strain when an edge has no neighbour, and Eq. (18) assigns the three pairwise
/// strains to the three Gauss points cyclically.
///
/// `entries[0]` is the target and `entries[1 + k]` is the element across local
/// edge `k`; `slots[entry][local_node]` is the union slot (0..6) that the
/// entry's local node occupies.  A `None` entry contributes nothing because the
/// boundary rule has already folded its weight onto the target.
///
/// The returned matrices act on the union's **local** DOF components - the
/// target's own for its nodes, each neighbour's own for its unique node - and
/// produce the target's local Cartesian membrane strain.  The finished stiffness
/// is rotated to global by [`transform_union_to_global`], which uses each
/// element's own frame, so no frame is mixed anywhere inside the operator.
pub fn smoothed_membrane_b(
    target: &Mitc3Precomputed,
    entries: &[Option<&Mitc3Precomputed>; 4],
    slots: &[[usize; 3]; 4],
    weights: &[[Matrix3<f64>; 4]; N_GAUSS],
) -> [Mat3Union; N_GAUSS] {
    // Each entry's membrane operator, converted to the target's convected
    // coordinates.  The target's own operator needs no convected transform.
    let mut covariant: [Option<SMatrix<f64, 3, 18>>; 4] = [None, None, None, None];
    for entry in 0..4 {
        let Some(pre) = entries[entry] else {
            continue;
        };
        // Everything here stays in the union's LOCAL DOFs: the target's own
        // components for its nodes and each neighbour's own components for its
        // unique node.  Mixing frames inside the operator is what made the
        // coupling wrong - `bm_union` would act on global components while `bk`,
        // `bg` and the drilling operator act on the target's local ones.  The
        // caller rotates the finished stiffness with `transform_union_to_global`.
        let own = crate::elements::smoothing::tensor_operator(&pre.j_mat) * b_membrane(&pre.dh);
        covariant[entry] = Some(if entry == 0 {
            own
        } else {
            crate::elements::smoothing::convected_operator(&target.j_mat, &pre.j_mat)
                .map_or(own.clone(), |transform| transform * own)
        });
    }

    // Back to the target's local Cartesian frame.  A singular target Jacobian is
    // degenerate geometry; the identity keeps the operator finite rather than
    // producing NaN.
    let back = crate::elements::smoothing::tensor_operator(
        &target.j_mat.try_inverse().unwrap_or_else(Matrix2::identity),
    );

    let mut out = [Mat3Union::zeros(); N_GAUSS];
    for gp in 0..N_GAUSS {
        let mut cov = Mat3Union::zeros();
        for entry in 0..4 {
            let Some(b_entry) = &covariant[entry] else {
                continue;
            };
            let w = &weights[gp][entry];
            if w.norm() < f64::MIN_POSITIVE {
                continue;
            }
            for local_node in 0..3 {
                let column = 6 * local_node;
                let target_column = 6 * slots[entry][local_node];
                for dof in 0..6 {
                    for row in 0..3 {
                        cov[(row, target_column + dof)] += w[(row, 0)] * b_entry[(0, column + dof)]
                            + w[(row, 1)] * b_entry[(1, column + dof)]
                            + w[(row, 2)] * b_entry[(2, column + dof)];
                    }
                }
            }
        }
        out[gp] = back * cov;
    }
    out
}

// ============================================================================
// Mass matrix computation
// ============================================================================

/// Compute the 18×18 consistent mass matrix in GLOBAL coordinates.
pub fn compute_me_global(pre: &Mitc3Precomputed, rho: f64) -> Mat18 {
    let area = pre.area;
    let h = pre.thickness;
    let m_trans = rho * h;
    let m_rot = rho * h * h * h / 12.0;

    let mut m_local = Mat18::zeros();

    for gp in 0..N_GAUSS {
        let hf = shape_functions(GAUSS_R[gp], GAUSS_S[gp]);
        let w = GAUSS_W[gp];

        for i in 0..3 {
            for j in 0..3 {
                let val_t = hf[i] * hf[j] * m_trans * w * area;
                let val_r = hf[i] * hf[j] * m_rot * w * area;

                for k in 0..3 {
                    m_local[(6 * i + k, 6 * j + k)] += val_t;
                }
                for k in 3..6 {
                    m_local[(6 * i + k, 6 * j + k)] += val_r;
                }
            }
        }
    }

    transform_to_global(&m_local, &pre.t3)
}

/// Compute the 18×18 consistent mass matrix for composite elements (global coords).
///
/// Uses pre-computed ply-integrated mass parameters instead of ρ·h.
pub fn compute_me_composite_global(
    pre: &Mitc3Precomputed,
    mass_per_area: f64,
    rotational_inertia: f64,
) -> Mat18 {
    let area = pre.area;
    let m_trans = mass_per_area;
    let m_rot = rotational_inertia;

    let mut m_local = Mat18::zeros();

    for gp in 0..N_GAUSS {
        let hf = shape_functions(GAUSS_R[gp], GAUSS_S[gp]);
        let w = GAUSS_W[gp];

        for i in 0..3 {
            for j in 0..3 {
                let val_t = hf[i] * hf[j] * m_trans * w * area;
                let val_r = hf[i] * hf[j] * m_rot * w * area;

                for k in 0..3 {
                    m_local[(6 * i + k, 6 * j + k)] += val_t;
                }
                for k in 3..6 {
                    m_local[(6 * i + k, 6 * j + k)] += val_r;
                }
            }
        }
    }

    transform_to_global(&m_local, &pre.t3)
}

// ============================================================================
// Transformation local <-> global
// ============================================================================

/// Transform an 18×18 matrix from local to global: T^T · K · T
fn transform_to_global(k_local: &Mat18, t3: &Matrix3<f64>) -> Mat18 {
    let t = t18(t3);
    t.transpose() * k_local * &t
}

/// The 18×18 block-diagonal rotation taking an element's global DOF components to
/// its local ones.  `t3`'s rows are the element's local basis vectors, so a
/// global vector's local components are `t3 * global`.
fn t18(t3: &Matrix3<f64>) -> Mat18 {
    let mut out = Mat18::zeros();
    for block in 0..6 {
        let row = 3 * block;
        for a in 0..3 {
            for b in 0..3 {
                out[(row + a, row + b)] = t3[(a, b)];
            }
        }
    }
    out
}

// ============================================================================
// Nonlinear: Displacement gradient, Green-Lagrange strain
// ============================================================================
// Uncited: Cartesian total-Lagrangian formulation with the CST field. Jeon, Lee,
// Lee & Bathe (2015), C&S 146:91-104, is the MITC3+ nonlinear paper, but it uses
// MITC-interpolated covariant strains, so no source is confirmed for this path.

/// Compute displacement gradient H = ∂u/∂X at centroid.
///
/// `u_local` is the 18-DOF displacement in local coordinates.
/// Returns 3×3 H matrix (only first 2 columns are non-zero for flat shell).
pub fn displacement_gradient(pre: &Mitc3Precomputed, u_local: &Vec18) -> Matrix3<f64> {
    let dh = &pre.dh;
    let mut h_mat = Matrix3::zeros();

    for i in 0..3 {
        let ux = u_local[6 * i];
        let uy = u_local[6 * i + 1];
        let uz = u_local[6 * i + 2];
        let dhi_dx = dh[(0, i)];
        let dhi_dy = dh[(1, i)];

        // H[comp, deriv_dir]
        h_mat[(0, 0)] += ux * dhi_dx;
        h_mat[(0, 1)] += ux * dhi_dy;
        h_mat[(1, 0)] += uy * dhi_dx;
        h_mat[(1, 1)] += uy * dhi_dy;
        h_mat[(2, 0)] += uz * dhi_dx;
        h_mat[(2, 1)] += uz * dhi_dy;
    }

    h_mat
}

/// Compute Green-Lagrange strain tensor: E = 0.5·(H + H^T + H^T·H)
pub fn green_lagrange_strain(h_mat: &Matrix3<f64>) -> Matrix3<f64> {
    0.5 * (h_mat + h_mat.transpose() + h_mat.transpose() * h_mat)
}

/// Green-Lagrange strain in membrane Voigt form: [Exx, Eyy, 2·Exy]
pub fn gl_strain_voigt(h_mat: &Matrix3<f64>) -> Vector3<f64> {
    let e = green_lagrange_strain(h_mat);
    Vector3::new(e[(0, 0)], e[(1, 1)], 2.0 * e[(0, 1)])
}

// ============================================================================
// Nonlinear: B_L, B_NL, B_geometric matrices
// ============================================================================

type Mat6x18 = SMatrix<f64, 6, 18>;

/// Linear strain-displacement matrix B_L (6×18) for TL formulation.
pub fn compute_b_l(dh: &Matrix2x3) -> Mat3x18 {
    // For membrane only (rows: εxx, εyy, γxy)
    b_membrane(dh)
}

/// Nonlinear strain-displacement matrix B_NL (3×18) depending on current displacement.
pub fn compute_b_nl(dh: &Matrix2x3, h_mat: &Matrix3<f64>) -> Mat3x18 {
    let mut b_nl = Mat3x18::zeros();

    for i in 0..3 {
        let col = 6 * i;
        let dni_dx = dh[(0, i)];
        let dni_dy = dh[(1, i)];

        // Exx nonlinear: (∂u/∂x)·δ(∂u/∂x) + (∂v/∂x)·δ(∂v/∂x) + (∂w/∂x)·δ(∂w/∂x)
        b_nl[(0, col)] = h_mat[(0, 0)] * dni_dx;
        b_nl[(0, col + 1)] = h_mat[(1, 0)] * dni_dx;
        b_nl[(0, col + 2)] = h_mat[(2, 0)] * dni_dx;

        // Eyy nonlinear
        b_nl[(1, col)] = h_mat[(0, 1)] * dni_dy;
        b_nl[(1, col + 1)] = h_mat[(1, 1)] * dni_dy;
        b_nl[(1, col + 2)] = h_mat[(2, 1)] * dni_dy;

        // 2*Exy nonlinear
        b_nl[(2, col)] = h_mat[(0, 0)] * dni_dy + h_mat[(0, 1)] * dni_dx;
        b_nl[(2, col + 1)] = h_mat[(1, 0)] * dni_dy + h_mat[(1, 1)] * dni_dx;
        b_nl[(2, col + 2)] = h_mat[(2, 0)] * dni_dy + h_mat[(2, 1)] * dni_dx;
    }

    b_nl
}

/// Geometric B-matrix (6×18) for K_σ computation.
pub fn compute_b_geometric(dh: &Matrix2x3) -> Mat6x18 {
    let mut bg = Mat6x18::zeros();

    for i in 0..3 {
        let col = 6 * i;
        let dni_dx = dh[(0, i)];
        let dni_dy = dh[(1, i)];

        bg[(0, col)] = dni_dx; // ∂u/∂x
        bg[(1, col)] = dni_dy; // ∂u/∂y
        bg[(2, col + 1)] = dni_dx; // ∂v/∂x
        bg[(3, col + 1)] = dni_dy; // ∂v/∂y
        bg[(4, col + 2)] = dni_dx; // ∂w/∂x
        bg[(5, col + 2)] = dni_dy; // ∂w/∂y
    }

    bg
}

// ============================================================================
// Nonlinear: Geometric stiffness K_σ
// ============================================================================

type Mat6 = SMatrix<f64, 6, 6>;

/// Compute geometric stiffness K_σ (18×18) in local coordinates.
///
/// `sigma_membrane` = [σxx, σyy, σxy]
pub fn compute_k_sigma_local(
    pre: &Mitc3Precomputed,
    sigma_membrane: &Vector3<f64>,
) -> Mat18 {
    let h = pre.thickness;

    // 2×2 stress block
    let s_m = Matrix2::new(
        sigma_membrane[0] * h,
        sigma_membrane[2] * h,
        sigma_membrane[2] * h,
        sigma_membrane[1] * h,
    );

    // 6×6 block diagonal
    let mut s_tilde = Mat6::zeros();
    s_tilde.fixed_view_mut::<2, 2>(0, 0).copy_from(&s_m);
    s_tilde.fixed_view_mut::<2, 2>(2, 2).copy_from(&s_m);
    s_tilde.fixed_view_mut::<2, 2>(4, 4).copy_from(&s_m);

    // For CST, B_G is constant → K_σ = area * B_G^T S̃ B_G
    let bg = compute_b_geometric(&pre.dh);
    let k_sigma = pre.area * (bg.transpose() * &s_tilde * &bg);

    0.5 * (&k_sigma + k_sigma.transpose())
}

// ============================================================================
// Nonlinear: Tangent stiffness K_T
// ============================================================================

/// Compute tangent stiffness K_T (18×18) in GLOBAL coordinates.
///
/// K_T = K_0 + K_L(u) + K_σ(u) for Total Lagrangian formulation.
///
/// `u_global` is the 18-DOF displacement vector in global coordinates.
pub fn compute_kt_global(pre: &Mitc3Precomputed, u_global: &Vec18) -> Mat18 {
    // Transform displacement to local
    let u_local = global_to_local_disp(u_global, &pre.t3);

    // Reference stiffness K_0
    let k0 = compute_ke_local(pre);

    // Displacement gradient
    let h_mat = displacement_gradient(pre, &u_local);

    // Current membrane stress
    let eps_gl = gl_strain_voigt(&h_mat);
    let sigma_m = &pre.constitutive.cm_raw * &eps_gl;

    // Geometric stiffness
    let k_sigma = compute_k_sigma_local(pre, &sigma_m);

    // Initial displacement stiffness K_L
    let dh = &pre.dh;
    let b_l = compute_b_l(dh);
    let b_nl = compute_b_nl(dh, &h_mat);
    let cm_raw = &pre.constitutive.cm_raw;
    let area = pre.area;
    let h = pre.thickness;

    // K_L = area·h·(B_L^T·C·B_NL + B_NL^T·C·B_L + B_NL^T·C·B_NL)
    let k_l = area
        * h
        * (b_l.transpose() * cm_raw * &b_nl
            + b_nl.transpose() * cm_raw * &b_l
            + b_nl.transpose() * cm_raw * &b_nl);

    let k_t_local = k0 + k_l + k_sigma;
    let k_t_local = 0.5 * (&k_t_local + k_t_local.transpose());

    transform_to_global(&k_t_local, &pre.t3)
}

// ============================================================================
// Nonlinear: Internal forces
// ============================================================================

/// Compute internal force vector f_int (18) in GLOBAL coordinates.
pub fn compute_fint_global(pre: &Mitc3Precomputed, u_global: &Vec18, nonlinear: bool) -> Vec18 {
    let u_local = global_to_local_disp(u_global, &pre.t3);
    let area = pre.area;
    let cm_raw = &pre.constitutive.cm_raw;
    let h = pre.thickness;
    let dh = &pre.dh;

    let mut f_int = Vec18::zeros();

    // Bending + shear (linear, condensed)
    let k_bs = compute_ke_local(pre) - {
        // Recompute just membrane + drilling to subtract
        let bm = b_membrane(dh);
        let km = area * (bm.transpose() * &pre.constitutive.cm * &bm);
        let mut kd = Mat18::zeros();
        for gp in 0..N_GAUSS {
            let bd = b_drill(GAUSS_R[gp], GAUSS_S[gp], dh);
            kd += (GAUSS_W[gp] * area * pre.k_drill) * (bd.transpose() * &bd);
        }
        km + kd
    };
    f_int += k_bs * &u_local;

    // Drilling
    for gp in 0..N_GAUSS {
        let bd = b_drill(GAUSS_R[gp], GAUSS_S[gp], dh);
        let drill_strain = &bd * &u_local;
        f_int += (pre.k_drill * GAUSS_W[gp] * area) * (bd.transpose() * &drill_strain);
    }

    // Membrane
    if nonlinear {
        let h_mat = displacement_gradient(pre, &u_local);
        let eps_gl = gl_strain_voigt(&h_mat);
        let sigma_m = cm_raw * &eps_gl;

        let b_l = compute_b_l(dh);
        let b_nl = compute_b_nl(dh, &h_mat);
        let b_total = b_l + b_nl;

        f_int += (area * h) * (b_total.transpose() * &sigma_m);
    } else {
        let bm = b_membrane(dh);
        let eps_m = &bm * &u_local;
        // Use cm_raw * eps → σ, then Bm^T * σ * h * area
        let sigma_m = cm_raw * &eps_m;
        f_int += (area * h) * (bm.transpose() * &sigma_m);
    }

    // Transform to global
    local_to_global_force(&f_int, &pre.t3)
}

// ============================================================================
// Displacement / force transformations
// ============================================================================

fn global_to_local_disp(u_global: &Vec18, t3: &Matrix3<f64>) -> Vec18 {
    let mut u_local = Vec18::zeros();
    for i in 0..6 {
        let g = u_global.fixed_rows::<3>(3 * i);
        let l = t3 * g;
        u_local.fixed_rows_mut::<3>(3 * i).copy_from(&l);
    }
    u_local
}

fn local_to_global_force(f_local: &Vec18, t3: &Matrix3<f64>) -> Vec18 {
    let mut f_global = Vec18::zeros();
    let t3t = t3.transpose();
    for i in 0..6 {
        let l = f_local.fixed_rows::<3>(3 * i);
        let g = &t3t * &l;
        f_global.fixed_rows_mut::<3>(3 * i).copy_from(&g);
    }
    f_global
}

// ============================================================================
// Body load
// ============================================================================

/// Compute body load vector f_body (18) in GLOBAL coordinates.
///
/// Integrates `f = ∫ Nᵀ · (ρ · h · g) dA` over the element area.
/// Only translational DOFs (indices 0,1,2 of each node block) receive
/// contributions; rotational DOFs (3,4,5) are zero.
///
/// `gravity` is the body-force acceleration vector in global coordinates [gx, gy, gz].
pub fn compute_body_load_global(
    pre: &Mitc3Precomputed,
    rho: f64,
    gravity: &Vector3<f64>,
) -> Vec18 {
    let area = pre.area;
    let h = pre.thickness;
    let rho_h = rho * h;

    let mut f_local = Vec18::zeros();

    for gp in 0..N_GAUSS {
        let hf = shape_functions(GAUSS_R[gp], GAUSS_S[gp]);
        let w = GAUSS_W[gp];

        // Transform gravity to local coordinates
        let g_local = pre.t3 * gravity;

        for i in 0..3 {
            let contrib = hf[i] * rho_h * w * area;
            // Translational DOFs only (0,1,2 of each 6-DOF node block)
            for k in 0..3 {
                f_local[6 * i + k] += contrib * g_local[k];
            }
        }
    }

    // Transform to global
    local_to_global_force(&f_local, &pre.t3)
}

// ============================================================================
// Geometric stiffness (global)
// ============================================================================

/// Compute geometric stiffness K_σ (18×18) in GLOBAL coordinates.
///
/// `sigma_membrane` = [σxx, σyy, σxy] in local coordinates.
pub fn compute_k_sigma_global(
    pre: &Mitc3Precomputed,
    sigma_membrane: &Vector3<f64>,
) -> Mat18 {
    let k_local = compute_k_sigma_local(pre, sigma_membrane);
    transform_to_global(&k_local, &pre.t3)
}

// ============================================================================
// Centrifugal prestress
// ============================================================================

/// Compute centrifugal prestress [σxx, σyy, σxy] in local coordinates.
///
/// Returns the membrane stress state due to centrifugal loading:
///     σ_cf ≈ ρ · ω² · r_radial · L_char
///
/// The stress is then projected into the local membrane plane using
/// the angle of the radial direction.
///
/// `omega`           : angular velocity (rad/s)
/// `rotation_axis`   : unit vector of rotation axis in global coords
/// `rotation_center` : a point on the rotation axis in global coords
/// `centroid`        : element centroid in global coords
/// `rho`             : material density (kg/m³)
pub fn compute_centrifugal_prestress(
    pre: &Mitc3Precomputed,
    omega: f64,
    rotation_axis: &Vector3<f64>,
    rotation_center: &Vector3<f64>,
    centroid: &Vector3<f64>,
    rho: f64,
) -> Vector3<f64> {
    // Normalize rotation axis
    let axis = rotation_axis.normalize();

    // Vector from rotation center to centroid
    let r_vec = centroid - rotation_center;

    // Project out the component along the axis → radial vector
    let r_parallel = r_vec.dot(&axis) * axis;
    let r_radial_vec = r_vec - r_parallel;
    let r_radial = r_radial_vec.norm();

    if r_radial < 1.0e-10 {
        // Element is on the rotation axis — no centrifugal stress
        return Vector3::zeros();
    }

    let radial_dir = r_radial_vec / r_radial;

    // Characteristic element length and stress magnitude
    let l_char = pre.area.sqrt();
    let sigma_cf = rho * omega * omega * r_radial * l_char;

    // Transform radial direction to local coordinates using the 3×3 rotation part
    let radial_local = pre.t3 * radial_dir;
    let cos_theta = radial_local[0];
    let sin_theta = radial_local[1];

    Vector3::new(
        sigma_cf * cos_theta * cos_theta,
        sigma_cf * sin_theta * sin_theta,
        sigma_cf * cos_theta * sin_theta,
    )
}

// ============================================================================
// Stress / Strain Recovery
// ============================================================================

/// Compute element stress and strain at the element centroid (r=1/3, s=1/3).
///
/// Returns `([sx, sy, sxy, 0, 0, 0], [exx, eyy, exy, 0, 0, 0])`.
///
/// # Arguments
/// * `u_global` - 18-DOF global displacement vector
/// * `z_factor` - Normalized through-thickness position: 0.0=mid, ±0.5=top/bottom
/// * `stress_type` - 0=membrane only, 1=bending only, 2=total
pub fn compute_element_stress(
    pre: &Mitc3Precomputed,
    u_global: &Vec18,
    z_factor: f64,
    stress_type: u8,
) -> ([f64; 6], [f64; 6]) {
    let u_local = global_to_local_disp(u_global, &pre.t3);
    let cm_raw = &pre.constitutive.cm_raw;
    let h = pre.thickness;
    let dh = &pre.dh;

    // Membrane (constant for CST)
    let bm = b_membrane(dh);
    let eps_m = bm * &u_local;
    let sig_m = cm_raw * eps_m;

    // Bending at centroid (r=1/3, s=1/3)
    let r_c = 1.0_f64 / 3.0;
    let s_c = 1.0_f64 / 3.0;
    let bk_ext = b_kappa_ext(r_c, s_c, &pre.j_inv);
    // Extract the first 18 columns (condensed, ignoring bubble DOFs)
    let bk18 = bk_ext.fixed_columns::<18>(0).into_owned();
    let kappa = bk18 * &u_local;
    let z = z_factor * h;
    let sig_b = cm_raw * kappa * z;

    // Apply stress_type flag
    let sig: nalgebra::Vector3<f64> = match stress_type {
        0 => sig_m,
        1 => sig_b,
        _ => sig_m + sig_b,
    };
    let eps_total: nalgebra::Vector3<f64> = match stress_type {
        0 => eps_m,
        1 => kappa * z,
        _ => eps_m + kappa * z,
    };

    let sigma6 = [sig[0], sig[1], 0.0, sig[2], 0.0, 0.0];
    let eps6 = [eps_total[0], eps_total[1], 0.0, eps_total[2], 0.0, 0.0];
    (sigma6, eps6)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::isotropic::IsotropicMaterial;
    use crate::materials::Material;

    /// Build a simple flat right-triangle in the XY plane.
    /// Nodes: (0,0,0), (1,0,0), (0,1,0)
    fn make_pre() -> Mitc3Precomputed {
        let thickness = 0.01_f64;
        let mat = IsotropicMaterial::new(2.0e11, 0.3, 7800.0);
        let shell = mat.constitutive(thickness, 5.0 / 6.0);
        let node_coords: [f64; 9] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        Mitc3Precomputed::new(&node_coords, shell, thickness, 2.0e11, 1.0)
    }

    #[test]
    fn test_body_load_global_zero_gravity() {
        let pre = make_pre();
        let g = Vector3::zeros();
        let f = compute_body_load_global(&pre, 7800.0, &g);
        assert!(f.norm() < 1e-12, "zero gravity → zero body load");
    }

    #[test]
    fn test_body_load_global_z_gravity() {
        let pre = make_pre();
        let g = Vector3::new(0.0, 0.0, -9.81);
        let f = compute_body_load_global(&pre, 7800.0, &g);

        // Only translational z-DOFs should be non-zero
        for i in 0..3 {
            // x, y translational and all rotational DOFs must be ~0
            assert!(f[6 * i].abs() < 1e-10, "node {i} fx should be ~0");
            assert!(f[6 * i + 1].abs() < 1e-10, "node {i} fy should be ~0");
            assert!(f[6 * i + 2].abs() > 1e-6, "node {i} fz should be nonzero");
            for k in 3..6 {
                assert!(f[6 * i + k].abs() < 1e-12, "node {i} rotational dof {k} should be 0");
            }
        }

        // Total z-force = ρ·h·g·area (area = 0.5 for this triangle)
        let area = 0.5_f64;
        let h = 0.01_f64;
        let rho = 7800.0_f64;
        let expected_total_fz = rho * h * (-9.81) * area;
        let total_fz: f64 = (0..3).map(|i| f[6 * i + 2]).sum();
        assert!(
            (total_fz - expected_total_fz).abs() < 1e-6,
            "total fz: got {total_fz}, expected {expected_total_fz}"
        );
    }

    #[test]
    fn test_k_sigma_global_zero_stress() {
        let pre = make_pre();
        let sigma = Vector3::zeros();
        let k = compute_k_sigma_global(&pre, &sigma);
        assert!(k.norm() < 1e-12, "zero stress → zero K_sigma");
    }

    #[test]
    fn test_k_sigma_global_symmetric() {
        let pre = make_pre();
        let sigma = Vector3::new(1.0e6, 0.5e6, 0.2e6);
        let k = compute_k_sigma_global(&pre, &sigma);
        let diff = k - k.transpose();
        assert!(
            diff.norm() < 1e-6 * k.norm().max(1.0),
            "K_sigma_global must be symmetric"
        );
    }

    #[test]
    fn test_k_sigma_global_matches_local_transformed() {
        let pre = make_pre();
        let sigma = Vector3::new(1.0e6, 0.5e6, 0.2e6);
        let k_local = compute_k_sigma_local(&pre, &sigma);
        let k_global_direct = compute_k_sigma_global(&pre, &sigma);
        let k_global_manual = transform_to_global(&k_local, &pre.t3);
        let diff = k_global_direct - k_global_manual;
        assert!(diff.norm() < 1e-6, "compute_k_sigma_global must equal transform(k_local)");
    }

    #[test]
    fn test_centrifugal_prestress_on_axis() {
        let pre = make_pre();
        // Centroid at (1/3, 1/3, 0); rotation axis = Z; center = origin
        // Element is at r_radial = sqrt((1/3)^2+(1/3)^2)
        // But if we put centroid ON the axis, expect zeros
        let axis = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::new(1.0 / 3.0, 1.0 / 3.0, 0.0); // same as centroid XY
        let centroid = Vector3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
        let sigma = compute_centrifugal_prestress(&pre, 100.0, &axis, &center, &centroid, 7800.0);
        assert!(sigma.norm() < 1e-6, "element on axis → zero centrifugal stress");
    }

    #[test]
    fn test_centrifugal_prestress_nonzero() {
        let pre = make_pre();
        let axis = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::zeros();
        let centroid = Vector3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
        let sigma = compute_centrifugal_prestress(&pre, 100.0, &axis, &center, &centroid, 7800.0);
        // All components should be finite and the trace (σxx + σyy) ≥ 0
        assert!(sigma[0].is_finite());
        assert!(sigma[1].is_finite());
        assert!(sigma[2].is_finite());
        assert!(sigma[0] + sigma[1] >= 0.0, "centrifugal stress trace must be non-negative");
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Global stiffness invariants: symmetry, PSD, rigid body, patch tests
    // ─────────────────────────────────────────────────────────────────────────

    /// Global coordinates of node `i` reconstructed from the element frame.
    ///
    /// `t3` maps global -> local (``local = t3 * global``), so the inverse is
    /// ``global = t3^T * local`` with the local z coordinate zero. This keeps
    /// the rigid-body helper generic even though MITC3 exposes no
    /// ``initial_coords_3d`` field.
    fn node_global_coords(pre: &Mitc3Precomputed, i: usize) -> Vector3<f64> {
        let x_local = pre.local_coords[2 * i];
        let y_local = pre.local_coords[2 * i + 1];
        pre.t3.transpose() * Vector3::new(x_local, y_local, 0.0)
    }

    /// Element centroid from the reconstructed global coordinates.
    fn element_centroid(pre: &Mitc3Precomputed) -> Vector3<f64> {
        let mut c = Vector3::zeros();
        for i in 0..3 {
            c += node_global_coords(pre, i);
        }
        c / 3.0
    }

    /// Physical rigid-body field: ``u = t + omega x (x - x_c)``, ``theta = omega``.
    ///
    /// This is the definition of a rigid body motion and does not depend on the
    /// element's own kinematic conventions, so it is the right way to test
    /// rigid-body invariance.
    fn rigid_body_mode(
        pre: &Mitc3Precomputed,
        translation: Vector3<f64>,
        omega: Vector3<f64>,
    ) -> Vec18 {
        let centroid = element_centroid(pre);
        let mut u = Vec18::zeros();
        for i in 0..3 {
            let x = node_global_coords(pre, i);
            let disp = translation + omega.cross(&(x - centroid));
            for k in 0..3 {
                u[6 * i + k] = disp[k];
                u[6 * i + 3 + k] = omega[k];
            }
        }
        u
    }

    fn rigid_body_modes(pre: &Mitc3Precomputed) -> [(&'static str, Vec18); 6] {
        let z = Vector3::zeros();
        let ex = Vector3::new(1.0, 0.0, 0.0);
        let ey = Vector3::new(0.0, 1.0, 0.0);
        let ez = Vector3::new(0.0, 0.0, 1.0);
        [
            ("translation x", rigid_body_mode(pre, ex, z)),
            ("translation y", rigid_body_mode(pre, ey, z)),
            ("translation z", rigid_body_mode(pre, ez, z)),
            ("rotation x", rigid_body_mode(pre, z, ex)),
            ("rotation y", rigid_body_mode(pre, z, ey)),
            ("rotation z", rigid_body_mode(pre, z, ez)),
        ]
    }

    /// ``|K u| / (|K| |u|)``: scale-free measure of an invariant residual.
    fn scaled_residual(k: &Mat18, u: &Vec18) -> f64 {
        let denom = k.norm() * u.norm();
        if denom > 0.0 { (k * u).norm() / denom } else { 0.0 }
    }

    #[test]
    fn test_ke_global_is_symmetric() {
        let pre = make_pre();
        let k = compute_ke_global(&pre);
        let asymmetry = (&k - k.transpose()).norm() / k.norm();
        assert!(
            asymmetry < 1e-12,
            "K_global must be symmetric: relative asymmetry = {asymmetry:.3e}"
        );
    }

    #[test]
    fn test_ke_global_is_positive_semidefinite() {
        let pre = make_pre();
        let k = compute_ke_global(&pre);
        let symmetric_part = (&k + k.transpose()) * 0.5;
        let eigenvalues = nalgebra::SymmetricEigen::new(symmetric_part).eigenvalues;
        let lambda_min = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        let lambda_max = eigenvalues.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            lambda_min > -1e-9 * lambda_max,
            "K_global must be positive semi-definite: \
             lambda_min = {lambda_min:.6e}, lambda_max = {lambda_max:.6e}"
        );
    }

    #[test]
    fn test_ke_global_leaves_all_six_rigid_body_modes_free() {
        let pre = make_pre();
        let k = compute_ke_global(&pre);

        let mut worst_label = "";
        let mut worst_residual = 0.0_f64;
        for (label, u) in rigid_body_modes(&pre) {
            let residual = scaled_residual(&k, &u);
            if residual > worst_residual {
                worst_label = label;
                worst_residual = residual;
            }
        }

        assert!(
            worst_residual < 1e-10,
            "all six rigid-body modes must be in the null space of K_global; \
             worst is '{worst_label}' with |K u| / (|K| |u|) = {worst_residual:.3e}"
        );
    }

    #[test]
    fn test_membrane_patch_reproduces_constant_strain_at_every_gauss_point() {
        let pre = make_pre();

        // u = a x + b y, v = c x + d y  =>  eps = [a, d, b + c], constant.
        let (a, b, c, d) = (1.0e-3, -4.0e-4, 2.0e-4, 7.0e-4);
        let mut u = Vec18::zeros();
        for i in 0..3 {
            let x = pre.local_coords[2 * i];
            let y = pre.local_coords[2 * i + 1];
            u[6 * i] = a * x + b * y;
            u[6 * i + 1] = c * x + d * y;
        }

        let expected = Vector3::new(a, d, b + c);
        let bm = b_membrane(&pre.dh);
        for g in 0..N_GAUSS {
            let eps = &bm * &u;
            let error = (eps - expected).norm() / expected.norm();
            // 1e-10 relative is far above the measured round-off of this
            // operator (~4e-12 for a 1e-3 field) and far below any real
            // formulation error, which would be O(1) relative.
            assert!(
                error < 1e-10,
                "membrane patch test at Gauss point {g}: eps = {eps:?}, \
                 expected = {expected:?}, relative error = {error:.3e}"
            );
        }
    }

    #[test]
    fn test_bending_patch_reproduces_constant_curvature_at_every_gauss_point() {
        let pre = make_pre();

        // w = 0.5 kxx x^2 with theta_y = dw/dx = kxx x  =>  kappa = [kxx, 0, 0].
        let kxx = 1.0e-3;
        let mut u20 = Vec20::zeros();
        for i in 0..3 {
            let x = pre.local_coords[2 * i];
            u20[6 * i + 2] = 0.5 * kxx * x * x;
            u20[6 * i + 4] = kxx * x;
        }

        // Entries 18 (theta_x4) and 19 (theta_y4) are the two internal bubble
        // rotations. `b_kappa_ext` is built from the ENRICHED shape derivatives
        // fi = hi - f4/3, so a constant-curvature state is represented only when
        // the bubble rotations carry the field value at the bubble node (the
        // element centroid in local coordinates); zeroing them instead feeds the
        // operator a NON-constant-curvature field and it correctly returns a
        // non-constant curvature.
        //
        // This is the API-forced deviation from the MITC4 reference: MITC4's
        // nodal `b_kappa` is unenriched, so its patch test needs no bubble DOF.
        // MITC3 exposes only the enriched 3x20 operator.
        let x_centroid = (pre.local_coords[0] + pre.local_coords[2] + pre.local_coords[4]) / 3.0;
        u20[19] = kxx * x_centroid; // theta_y4 = theta_y(centroid)
        u20[18] = 0.0; // theta_x4 = theta_x(centroid) = 0

        let expected = Vector3::new(kxx, 0.0, 0.0);
        for g in 0..N_GAUSS {
            let kappa = b_kappa_ext(GAUSS_R[g], GAUSS_S[g], &pre.j_inv) * u20;
            let error = (kappa - expected).norm() / expected.norm();
            // Same rationale as the membrane patch: ~4e-12 round-off for a
            // 1e-3 field vs. O(1) for a real formulation error.
            assert!(
                error < 1e-10,
                "bending patch test at Gauss point {g}: kappa = {kappa:?}, \
                 expected = {expected:?}, relative error = {error:.3e}"
            );
        }
    }
    // ─────────────────────────────────────────────────────────────────────────
    // Mass matrix invariants
    // ─────────────────────────────────────────────────────────────────────────

    /// Total translational mass per direction: the sum of the 3x3 sub-block of
    /// one direction.  For a consistent mass matrix built from a partition of
    /// unity this equals `rho * h * A` exactly, for every element type.
    fn translational_mass_per_direction(m: &Mat18) -> [f64; 3] {
        let mut totals = [0.0; 3];
        for (d, total) in totals.iter_mut().enumerate() {
            let mut sum = 0.0;
            for i in (d..18).step_by(6) {
                for j in (d..18).step_by(6) {
                    sum += m[(i, j)];
                }
            }
            *total = sum;
        }
        totals
    }

    /// Total rotary-inertia mass per rotation direction: `rho * h^3 / 12 * A`.
    fn rotary_mass_per_direction(m: &Mat18) -> [f64; 3] {
        let mut totals = [0.0; 3];
        for (k, total) in totals.iter_mut().enumerate() {
            let mut sum = 0.0;
            for i in (3 + k..18).step_by(6) {
                for j in (3 + k..18).step_by(6) {
                    sum += m[(i, j)];
                }
            }
            *total = sum;
        }
        totals
    }

    #[test]
    fn test_me_global_is_symmetric_and_positive_semidefinite() {
        let pre = make_pre();
        let m = compute_me_global(&pre, 7800.0);

        let asymmetry = (&m - m.transpose()).norm() / m.norm();
        assert!(
            asymmetry < 1e-14,
            "M_global must be symmetric: relative asymmetry = {asymmetry:.3e}"
        );

        let symmetric_part = (&m + m.transpose()) * 0.5;
        let eigenvalues = nalgebra::SymmetricEigen::new(symmetric_part).eigenvalues;
        let lambda_min = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        let lambda_max = eigenvalues
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            lambda_min > -1e-12 * lambda_max,
            "M_global must be positive semi-definite: \
             lambda_min = {lambda_min:.6e}, lambda_max = {lambda_max:.6e}"
        );
    }

    #[test]
    fn test_me_global_total_translational_mass_is_rho_h_a() {
        let pre = make_pre();
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let expected = rho * pre.thickness * pre.area;

        for (d, total) in translational_mass_per_direction(&m).iter().enumerate() {
            let error = (total - expected).abs() / expected;
            assert!(
                error < 1e-14,
                "direction {d}: total mass {total:.6e} != rho*h*A {expected:.6e} \
                 (relative error {error:.3e})"
            );
        }
    }

    #[test]
    fn test_me_global_matches_the_exact_linear_triangle_coefficients() {
        let pre = make_pre();
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let mass = rho * pre.thickness * pre.area;

        // Linear-triangle consistent mass: per translational direction
        // M_ii = m/6 and M_ij = m/12.  The 3-point rule is exact for these
        // quadratic products, so the tolerance is round-off.
        for d in 0..3 {
            for i in 0..3 {
                for j in 0..3 {
                    let expected = if i == j { mass / 6.0 } else { mass / 12.0 };
                    let actual = m[(6 * i + d, 6 * j + d)];
                    let error = (actual - expected).abs() / mass;
                    assert!(
                        error < 1e-14,
                        "M[{i},{j}] direction {d}: {actual:.6e} != {expected:.6e} \
                         (relative error {error:.3e})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_me_global_rotary_inertia_is_rho_h3_a_over_12() {
        let pre = make_pre();
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let expected = rho * pre.thickness.powi(3) / 12.0 * pre.area;

        for (k, total) in rotary_mass_per_direction(&m).iter().enumerate() {
            let error = (total - expected).abs() / expected;
            assert!(
                error < 1e-14,
                "rotation direction {k}: rotary mass {total:.6e} != \
                 rho*h^3/12*A {expected:.6e} (relative error {error:.3e})"
            );
        }
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Strain-smoothed membrane field (Lee & Lee 2019, Eqs. 15-18)
    // ─────────────────────────────────────────────────────────────────────────

    /// The triangle frame the smoothing needs, taken from the precomputed data.
    fn triangle_frame(pre: &Mitc3Precomputed) -> crate::elements::smoothing::TriangleFrame {
        crate::elements::smoothing::TriangleFrame {
            j_mat: pre.j_mat,
            area: pre.area,
            normal: Vector3::new(pre.t3[(2, 0)], pre.t3[(2, 1)], pre.t3[(2, 2)]),
        }
    }

    #[test]
    fn union_rotation_rotates_the_full_six_dof_block() {
        // Issue #2: `union_rotation` filled only the translational 3x3 block, so
        // `transform_union_to_global` = T^T K T annihilated every rotational DOF.
        // The 6-DOF node block must carry the frame in BOTH the translational and
        // the rotational slots, which also makes T orthogonal and non-singular.
        let pre = make_pre();
        let frames = [pre.t3; SMOOTHED_UNION_NODES];
        let t = union_rotation(&frames);

        assert!(
            t.determinant().abs() > 1e-12,
            "union_rotation must be non-singular; a zero rotational block makes it singular"
        );
        for slot in 0..SMOOTHED_UNION_NODES {
            let base = 6 * slot;
            for a in 0..3 {
                for b in 0..3 {
                    assert!((t[(base + a, base + b)] - pre.t3[(a, b)]).abs() < 1e-12);
                    assert!(
                        (t[(base + 3 + a, base + 3 + b)] - pre.t3[(a, b)]).abs() < 1e-12,
                        "the rotational block must equal the frame (issue #2)"
                    );
                }
            }
        }

        // An orthogonal transform preserves the Frobenius norm of any stiffness.
        let mut k = MatUnion::zeros();
        for i in 0..k.nrows() {
            k[(i, i)] = 1.0 + i as f64;
        }
        let rotated = transform_union_to_global(&k, &frames);
        assert!(
            (rotated.norm() - k.norm()).abs() < 1e-10 * k.norm(),
            "an orthogonal union rotation must preserve the stiffness norm"
        );
    }

    /// The 18 local nodal displacements of the linear field
    /// `u = [a x + b y, c x + d y]`, whose constant strain is `[a, d, b + c]`.
    fn linear_field(pre: &Mitc3Precomputed, a: f64, b: f64, c: f64, d: f64) -> Vec18 {
        let mut u = Vec18::zeros();
        for i in 0..3 {
            let x = pre.local_coords[2 * i];
            let y = pre.local_coords[2 * i + 1];
            u[6 * i] = a * x + b * y;
            u[6 * i + 1] = c * x + d * y;
        }
        u
    }

    /// Scatter an entry's local displacements into the union layout.
    fn place_in_union(
        slots: [usize; 3],
        u_local: &Vec18,
        u_union: &mut SMatrix<f64, { 6 * SMOOTHED_UNION_NODES }, 1>,
    ) {
        for local_node in 0..3 {
            for dof in 0..6 {
                u_union[6 * slots[local_node] + dof] = u_local[6 * local_node + dof];
            }
        }
    }

    #[test]
    fn the_union_stiffness_reduces_to_the_element_with_the_coupling_active() {
        // The isotropic reduction test cannot see the membrane-bending coupling,
        // because cb_coupling is identically zero for an isotropic material.
        // This one switches the coupling on, so the union treatment of the
        // coupling's two cross terms is actually exercised.
        let mut shell =
            IsotropicMaterial::new(2.0e11, 0.3, 7800.0).constitutive(0.01, 5.0 / 6.0);
        shell.cb_coupling = SMatrix::<f64, 3, 3>::identity() * 1.0e6;
        let coords: [f64; 9] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let pre = Mitc3Precomputed::new(&coords, shell, 0.01, 2.0e11, 1.0);

        let frames = [triangle_frame(&pre)];
        let neighbours = [[None, None, None]];
        let operator = crate::elements::smoothing::smoothed_membrane_strain(&frames, &neighbours);
        let entries = [Some(&pre), None, None, None];
        let slots = [[0usize, 1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]];
        let bm_union = smoothed_membrane_b(&pre, &entries, &slots, &operator[0].weights);

        let union = compute_ke_union_smoothed(&pre, &bm_union);
        let own = compute_ke_local(&pre);

        let top_left = union.fixed_view::<18, 18>(0, 0).into_owned();
        let difference = (top_left - own).norm() / own.norm();
        assert!(
            difference < 1e-12,
            "with the coupling active the union stiffness must still reduce to the element, \
             relative difference = {difference:.3e}"
        );
    }

    #[test]
    fn the_union_stiffness_reduces_to_the_element_without_neighbours() {
        // With a membrane operator whose first 18 columns are the element's own
        // b_membrane and whose remaining columns are zero - exactly what the
        // boundary rule produces for a triangle with no neighbours - the union
        // stiffness must reproduce compute_ke_local in its top-left block.  This
        // is the property that makes the smoothed path a strict extension of the
        // 2014 MITC3+ rather than a replacement, and it covers the drilling term,
        // the bubble condensation and the coupling's two cross terms at once.
        let pre = make_pre();
        let frames = [triangle_frame(&pre)];
        let neighbours = [[None, None, None]];
        let operator = crate::elements::smoothing::smoothed_membrane_strain(&frames, &neighbours);
        let entries = [Some(&pre), None, None, None];
        let slots = [[0usize, 1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]];
        let bm_union = smoothed_membrane_b(&pre, &entries, &slots, &operator[0].weights);

        let union = compute_ke_union_smoothed(&pre, &bm_union);
        // Both sides in the target's local DOFs: the union assembly is local, and
        // transform_union_to_global does the rotation for the assembled system.
        let own = compute_ke_local(&pre);

        let top_left = union.fixed_view::<18, 18>(0, 0).into_owned();
        let difference = (top_left - own).norm() / own.norm();
        // A relative tolerance: the union path performs the same arithmetic in a
        // different order and on wider matrices, so the two agree to round-off
        // rather than bit for bit.  The stiffness norm is of the order of 1e9, so
        // an absolute threshold here would be measuring nothing.
        assert!(
            difference < 1e-12,
            "the union stiffness must reduce to the element, relative difference = {difference:.3e}"
        );

        // The unused neighbour slots must stay uncoupled.
        let coupling = union.fixed_view::<18, 18>(0, 18).into_owned();
        assert!(
            coupling.norm() < 1e-12 * own.norm(),
            "unexpected neighbour coupling = {:.3e}",
            coupling.norm()
        );
    }

    #[test]
    fn default_membrane_path_is_the_un_smoothed_one() {
        // Passing the element's own membrane operator at every Gauss point must
        // reproduce compute_ke_local exactly.  This is what keeps the MITC3+ of
        // Lee, Lee & Bathe 2014 available while the smoothed field is added.
        let pre = make_pre();
        let bm = b_membrane(&pre.dh);
        let with = compute_ke_local_with_membrane(&pre, &[bm; N_GAUSS]);
        let without = compute_ke_local(&pre);
        let difference = (with - without).norm();
        assert!(
            difference < 1e-15,
            "the un-smoothed path must be unchanged, difference = {difference:.3e}"
        );
    }

    #[test]
    fn a_boundary_triangle_reproduces_its_own_membrane_operator() {
        // With no neighbours the boundary rule makes every pairwise strain the
        // element's own, so the smoothed operator must equal b_membrane in the
        // target's columns and put nothing on the neighbour slots.
        let pre = make_pre();
        let frames = [triangle_frame(&pre)];
        let neighbours = [[None, None, None]];
        let operator = crate::elements::smoothing::smoothed_membrane_strain(&frames, &neighbours);

        let entries = [Some(&pre), None, None, None];
        let slots = [[0usize, 1, 2], [3, 4, 5], [3, 4, 5], [3, 4, 5]];
        let smoothed = smoothed_membrane_b(&pre, &entries, &slots, &operator[0].weights);

        let bm = b_membrane(&pre.dh);
        for gp in 0..N_GAUSS {
            let own = smoothed[gp].fixed_columns::<18>(0).into_owned();
            assert!(
                (own - bm).norm() < 1e-15,
                "Gauss point {gp}: own columns differ by {:.3e}",
                (own - bm).norm()
            );
            let others = smoothed[gp].fixed_columns::<18>(18).into_owned();
            assert!(
                others.norm() < 1e-15,
                "Gauss point {gp}: neighbour block should be zero, got {:.3e}",
                others.norm()
            );
        }
    }

    #[test]
    fn the_smoothed_field_passes_the_membrane_patch_test() {
        // A constant strain field must survive the smoothing: each element's
        // smoothed operator, applied to that field, must return the same strain.
        // This exercises the covariant transforms of Eq. (15) and the assignment
        // of Eq. (18), which are the parts that can silently scale or rotate the
        // field.
        let target = make_pre();
        let other = make_pre();
        let frames = [triangle_frame(&target), triangle_frame(&other)];
        let neighbours = [[Some(1), Some(1), Some(1)], [Some(0), None, None]];
        let operator = crate::elements::smoothing::smoothed_membrane_strain(&frames, &neighbours);

        let entries = [Some(&target), Some(&other), Some(&other), Some(&other)];
        let slots = [[0usize, 1, 2], [3, 1, 2], [3, 1, 2], [3, 1, 2]];
        let smoothed = smoothed_membrane_b(&target, &entries, &slots, &operator[0].weights);

        let (a, b, c, d) = (1.0e-3, -4.0e-4, 2.0e-4, 7.0e-4);
        let expected = Vector3::new(a, d, b + c);
        let u_local = linear_field(&target, a, b, c, d);

        let mut u_union = SMatrix::<f64, { 6 * SMOOTHED_UNION_NODES }, 1>::zeros();
        place_in_union(slots[0], &u_local, &mut u_union);
        place_in_union(slots[1], &u_local, &mut u_union);
        place_in_union(slots[2], &u_local, &mut u_union);
        place_in_union(slots[3], &u_local, &mut u_union);

        for gp in 0..N_GAUSS {
            let strain = smoothed[gp] * u_union;
            let error = (strain - expected).norm() / expected.norm();
            assert!(
                error < 1e-12,
                "Gauss point {gp}: smoothed strain {strain:?} != {expected:?} (rel {error:.3e})"
            );
        }
    }
}
