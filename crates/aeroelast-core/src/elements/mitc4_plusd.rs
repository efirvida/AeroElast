// MITC4+/D Shell Element Kernel — WU1: patch-test fixtures
//
// The element is the MITC4+/D: the 2017 MITC4+ of Ko, Lee & Bathe (2017),
// Computers and Structures (C&S) 182:404–418, plus the 2025 penalty-free
// drill-membrane strain of Ko, Bathe & Zhang (2025), C&S 308:107622.
//
// SHORTHAND NOTE: the extract documents under `docs/formulations/` label these
// two papers with a letter shorthand. `mitc4plusd-2025-extract.md` calls the
// 2017 formulation paper "Paper A" and the 2025 MITC4+/D paper "Paper B";
// `mitc4plus-2017-extract.md` reuses the same letters for a different pair (its
// "Paper B" is Ko, Lee, Lee & Bathe (2017), C&S 193:187–206). This module never
// relies on that shorthand: every citation below is self-contained as
// author-year plus journal, volume and pages.
//
// Ko, Lee & Bathe (2017), C&S 182:404–418, gives the 3D continuum kinematics
// (Eqs. 1–3), the assumed membrane field (Eqs. 17–27), the displacement-based
// bending terms (Eqs. 7c/7d) and the MITC4 / Dvorkin–Bathe transverse-shear
// field (DB84 Eq. 3, reproduced in Ko, Lee & Bathe (2017), C&S 182:404–418,
// p. 405). Ko, Bathe & Zhang (2025), C&S 308:107622, gives the sixth (drilling)
// DOF through the drill-membrane strain of Eq. (26) with the operator of
// Eq. (18). The element contains no numerical factor: no penalty, no shear
// correction factor, no selective reduced integration.
//
// SCOPE (WU1): only the patch-test fixtures and their self-tests. The
// production element lands in WU2–WU5. Every coordinate is a figure read
// recorded in `docs/formulations/mitc4plusd-2025-extract.md` §Fig. 7(a);
// Ko, Lee & Bathe (2017), C&S 182:404–418, Fig. 5 is the same mesh.
//
// CORRECTIONS APPLIED AND RECORDED
//   1. TASK 2.4 IS DEFERRED: it assembles a 48-DOF dense patch matrix and needs
//      the `Mitc4PlusDPrecomputed` type, which does not exist until WU2–WU4.
//      Writing it now would not produce a red test but a broken build, so it
//      lands with WU4 and stays unchecked in `tasks.md`.
//   2. MODULE-LAYOUT DEVIATION: `tasks.md` names
//      `elements/mitc4_plusd/tests/fixtures.rs`, but the design fixes
//      `elements/mitc4_plusd.rs` and every existing element (`mitc3.rs`,
//      `mitc4.rs`, `quad.rs`) is a single file with an inline `#[cfg(test)]
//      mod tests`. Per `openspec/config.yaml`'s apply guideline (follow the
//      existing Rust element patterns), the fixtures live inline here and this
//      file is registered in `elements/mod.rs`.
//   3. BC_2017_PATCH DERIVATION: the 2017 Tier-1a minimum boundary conditions
//      are derived, not figure-read. Design §4.1's set does not remove all six
//      rigid-body modes (its `u_x = 0` at C and D is the same condition because
//      C = (0,0) and D = (10,0) share `y = 0`, and it constrains neither θ_x nor
//      θ_y, so the two tilt modes survive). The corrected set is on
//      `BC_2017_PATCH`.

// ============================================================================
// WU2 — the element core: geometry, 2017 membrane coefficients, nodal
//        directors and kinematics
// ============================================================================
//
// Citations are self-contained (author-year plus journal, volume and pages),
// anchored to `docs/references.md`; the extract documents' local "paper A/B"
// shorthand is never used here.
//
//   Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
//   Structures 182:404-418.
//   Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements with
//   six degrees of freedom at each node - the MITC4/D and MITC4+/D elements",
//   Computers and Structures 308:107622.

use nalgebra::{Matrix2, Matrix3, Vector3};

use crate::materials::ShellConstitutive;

// ============================================================================
// Shape functions
// ============================================================================

/// Bilinear shape functions `h_i(r,s)`, `i = 0..3`, of Ko, Bathe & Zhang
/// (2025), C&S 308:107622, Eq. (2), p. 2 — written in the repository's node
/// order (node 0 at `(r,s) = (-1,-1)`, counter-clockwise).
#[inline(always)]
fn shape_functions(r: f64, s: f64) -> [f64; 4] {
    [
        0.25 * (1.0 - r) * (1.0 - s),
        0.25 * (1.0 + r) * (1.0 - s),
        0.25 * (1.0 + r) * (1.0 + s),
        0.25 * (1.0 - r) * (1.0 + s),
    ]
}

/// Derivatives `(dh_i/dr, dh_i/ds)` of the bilinear shape functions of
/// Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (2), p. 2.
/// Returns `(dh/dr[4], dh/ds[4])`.
#[inline(always)]
fn shape_function_derivatives(r: f64, s: f64) -> ([f64; 4], [f64; 4]) {
    let dn_dr = [
        -0.25 * (1.0 - s),
        0.25 * (1.0 - s),
        0.25 * (1.0 + s),
        -0.25 * (1.0 + s),
    ];
    let dn_ds = [
        -0.25 * (1.0 - r),
        -0.25 * (1.0 + r),
        0.25 * (1.0 + r),
        0.25 * (1.0 - r),
    ];
    (dn_dr, dn_ds)
}

// ============================================================================
// Local orthonormal frame
// ============================================================================

/// Local orthonormal frame `(e1, e2, e3)` and the projected 2D node
/// coordinates of the four nodes.
///
/// `e3` is the average of the two diagonal triangle normals of the quad, as
/// used by the repository's existing shell elements; `e1` is edge `0 -> 1`
/// orthogonalized against `e3`; `e2 = e3 x e1`.
///
/// The frame is the repository's local reference system, not a paper equation;
/// it is the same construction the assembly layers already use. The paper's own
/// plane normal is `n_vec` (Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (10))
/// and is computed separately by [`compute_characteristic_vectors`].
///
/// Returns `(local_coords[4][2], e1, e2, e3)`.
fn compute_local_coordinate_system(
    coords_3d: &[[f64; 3]; 4],
) -> ([[f64; 2]; 4], Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let nodes: [Vector3<f64>; 4] = [
        Vector3::new(coords_3d[0][0], coords_3d[0][1], coords_3d[0][2]),
        Vector3::new(coords_3d[1][0], coords_3d[1][1], coords_3d[1][2]),
        Vector3::new(coords_3d[2][0], coords_3d[2][1], coords_3d[2][2]),
        Vector3::new(coords_3d[3][0], coords_3d[3][1], coords_3d[3][2]),
    ];

    // Average normal of the two diagonal triangles (0,1,2) and (0,2,3).
    let n1 = (nodes[1] - nodes[0]).cross(&(nodes[2] - nodes[0]));
    let n2 = (nodes[2] - nodes[0]).cross(&(nodes[3] - nodes[0]));
    let mut e3 = Vector3::zeros();
    let mut count = 0;
    if n1.norm() > 1e-12 {
        e3 += n1.normalize();
        count += 1;
    }
    if n2.norm() > 1e-12 {
        e3 += n2.normalize();
        count += 1;
    }
    if count > 0 {
        e3 = e3.normalize();
    } else {
        e3 = Vector3::new(0.0, 0.0, 1.0);
    }

    // e1 from edge 0 -> 1, orthogonalized against e3.
    let mut e1 = nodes[1] - nodes[0];
    e1 -= e1.dot(&e3) * e3;
    if e1.norm() < 1e-12 {
        e1 = nodes[2] - nodes[0];
        e1 -= e1.dot(&e3) * e3;
    }
    e1 = e1.normalize();
    let e2 = e3.cross(&e1).normalize();

    let mut local_coords = [[0.0f64; 2]; 4];
    for i in 0..4 {
        local_coords[i][0] = nodes[i].dot(&e1);
        local_coords[i][1] = nodes[i].dot(&e2);
    }
    (local_coords, e1, e2, e3)
}

// ============================================================================
// Characteristic vectors and the dual basis
// ============================================================================

/// Characteristic vectors and the dual basis of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (9)-(11), p. 406.
///
/// ```text
/// x_r = 1/4 sum xi_i x_i,  x_s = 1/4 sum eta_i x_i,
/// x_d = 1/4 sum xi_i eta_i x_i                                    (Eq. 9)
/// n   = (x_r x x_s)/||x_r x x_s||                                 (Eq. 10)
/// m^r.x_r = m^s.x_s = 1,  m^r.x_s = m^s.x_r = 0,  m^r.n = m^s.n = 0 (Eq. 11)
/// ```
///
/// The dual basis is obtained by solving the 3x3 Gram system
/// `[x_r x_s n]^T [m^r m^s] = [e_r e_s]`; the two right-hand sides share the
/// same factorization. Returns `(x_r, x_s, x_d, n_vec, m_r, m_s)`.
fn compute_characteristic_vectors(
    coords_3d: &[[f64; 3]; 4],
) -> (
    Vector3<f64>,
    Vector3<f64>,
    Vector3<f64>,
    Vector3<f64>,
    Vector3<f64>,
    Vector3<f64>,
) {
    let nodes: [Vector3<f64>; 4] = [
        Vector3::new(coords_3d[0][0], coords_3d[0][1], coords_3d[0][2]),
        Vector3::new(coords_3d[1][0], coords_3d[1][1], coords_3d[1][2]),
        Vector3::new(coords_3d[2][0], coords_3d[2][1], coords_3d[2][2]),
        Vector3::new(coords_3d[3][0], coords_3d[3][1], coords_3d[3][2]),
    ];

    let x_r = 0.25 * (-nodes[0] + nodes[1] + nodes[2] - nodes[3]);
    let x_s = 0.25 * (-nodes[0] - nodes[1] + nodes[2] + nodes[3]);
    let x_d = 0.25 * (nodes[0] - nodes[1] + nodes[2] - nodes[3]);

    let mut n_vec = x_r.cross(&x_s);
    if n_vec.norm() < 1e-12 {
        let v = nodes[1] - nodes[0];
        let w = nodes[2] - nodes[0];
        n_vec = v.cross(&w);
    }
    n_vec = n_vec.normalize();

    let a_mat = Matrix3::new(
        x_r.dot(&x_r),
        x_r.dot(&x_s),
        x_r.dot(&n_vec),
        x_s.dot(&x_r),
        x_s.dot(&x_s),
        x_s.dot(&n_vec),
        n_vec.dot(&x_r),
        n_vec.dot(&x_s),
        n_vec.dot(&n_vec),
    );

    let (m_r, m_s) = if let Some(a_inv) = a_mat.try_inverse() {
        let cr = a_inv * Vector3::new(1.0, 0.0, 0.0);
        let cs = a_inv * Vector3::new(0.0, 1.0, 0.0);
        (
            cr[0] * x_r + cr[1] * x_s + cr[2] * n_vec,
            cs[0] * x_r + cs[1] * x_s + cs[2] * n_vec,
        )
    } else {
        // Degenerate geometry: fall back to the orthogonal projections.
        (x_r / x_r.dot(&x_r), x_s / x_s.dot(&x_s))
    };

    (x_r, x_s, x_d, n_vec, m_r, m_s)
}

/// Regularized inverse of a 2x2 matrix, used by the covariant-to-local maps.
///
/// A tiny diagonal shift (relative to the matrix scale) keeps the inverse
/// finite on degenerate geometry. It is a division guard, not a numerical
/// factor in the formulation.
#[inline(always)]
fn regularized_inverse_2x2(m: &Matrix2<f64>) -> Matrix2<f64> {
    let scale = m[(0, 0)]
        .abs()
        .max(m[(0, 1)].abs())
        .max(m[(1, 0)].abs())
        .max(m[(1, 1)].abs())
        .max(1.0);
    let mut reg = *m;
    let eps = scale * 1.0e-12;
    reg[(0, 0)] += eps;
    reg[(1, 1)] += eps;
    reg.try_inverse().unwrap_or_else(Matrix2::identity)
}

/// Covariant-to-local membrane strain mapping of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (15)-(19), p. 408.
///
/// Maps the covariant strain triple `(e_rr, e_ss, 2 e_rs)` to the local
/// orthonormal-frame triple `(e_11, e_22, 2 e_12)` through the point-wise
/// Jacobian `j_loc`. The mapping is `J0` in the paper's notation.
fn covariant_to_local_mapping(j_loc: &Matrix2<f64>) -> Matrix3<f64> {
    let j_inv = regularized_inverse_2x2(j_loc);
    let j11 = j_inv[(0, 0)];
    let j12 = j_inv[(0, 1)];
    let j21 = j_inv[(1, 0)];
    let j22 = j_inv[(1, 1)];
    Matrix3::new(
        j11 * j11,
        j21 * j21,
        j11 * j21,
        j12 * j12,
        j22 * j22,
        j12 * j22,
        2.0 * j11 * j12,
        2.0 * j21 * j22,
        j11 * j22 + j12 * j21,
    )
}

/// Covariant-to-local mapping for the transverse shear.
///
/// Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (4), p. 405, defines the
/// covariant transverse shear as `e_rt`, `e_st`; the paper does not print the
/// metric normalization to the local orthonormal frame (recorded as WU0's
/// task 1.4 result). Following the design's exact definition, the contraction
/// with the 3D dual basis is
///
/// ```text
/// gamma_a3 = 2 e_i3 (g^i . e_a)(g^t . e_3),   i = r, s
/// g^r = (g_s x g_t)/j,  g^s = (g_t x g_r)/j,  g^t = (g_r x g_s)/j,
/// j = det[g_r g_s g_t]
/// ```
///
/// with no invented factor. Returns the 2x2 matrix `T` such that
/// `[gamma_13, gamma_23] = T [e_rt, e_st]`.
fn shear_covariant_to_local(
    g_r: &Vector3<f64>,
    g_s: &Vector3<f64>,
    g_t: &Vector3<f64>,
    e1: &Vector3<f64>,
    e2: &Vector3<f64>,
    e3: &Vector3<f64>,
) -> Matrix2<f64> {
    let j = g_r.cross(g_s).dot(g_t);
    if j.abs() < 1e-30 {
        return Matrix2::zeros();
    }
    let g_r_dual = g_s.cross(g_t) / j;
    let g_s_dual = g_t.cross(g_r) / j;
    let g_t_dual = g_r.cross(g_s) / j;
    let ft = g_t_dual.dot(e3);
    Matrix2::new(
        2.0 * g_r_dual.dot(e1) * ft,
        2.0 * g_s_dual.dot(e1) * ft,
        2.0 * g_r_dual.dot(e2) * ft,
        2.0 * g_s_dual.dot(e2) * ft,
    )
}

// ============================================================================
// Nodal directors (ADR-4 option B)
// ============================================================================

/// Per-node directors `V_n^i`, `V_1^i`, `V_2^i` of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (1), (3), (8a), pp. 405-406.
///
/// There is no per-node director input anywhere in the repository, so the
/// element builds it from its own four coordinates (ADR-4 option B): split the
/// quad at its centre `x_c = 1/4 sum x_i` into the four triangles
/// `(x_c, x_k, x_{k+1})`; node `i` is shared by the triangles `k = i-1` and
/// `k = i`, so
///
/// ```text
/// V_n^i = normalize( A_{i-1} n_{i-1} + A_i n_i )
/// ```
///
/// with `n_k` the outward (sign-aligned to `n_vec`) unit normal and `A_k` the
/// area of triangle `k`. This is the local analogue of a mesh-averaged nodal
/// normal and is exact (`V_n^i = n_vec`) for flat geometry.
///
/// The in-plane pair is derived from the stored director and the local frame's
/// `e1`: `V_1^i = normalize(e1 - (e1.V_n^i) V_n^i)`, `V_2^i = V_n^i x V_1^i`,
/// so `(V_1^i, V_2^i, V_n^i)` is right-handed and Eq. (3d) holds.
///
/// Returns `(vn, v1, v2)`.
fn compute_node_directors(
    coords_3d: &[[f64; 3]; 4],
    n_vec: &Vector3<f64>,
    e1: &Vector3<f64>,
) -> ([Vector3<f64>; 4], [Vector3<f64>; 4], [Vector3<f64>; 4]) {
    let x: [Vector3<f64>; 4] = [
        Vector3::new(coords_3d[0][0], coords_3d[0][1], coords_3d[0][2]),
        Vector3::new(coords_3d[1][0], coords_3d[1][1], coords_3d[1][2]),
        Vector3::new(coords_3d[2][0], coords_3d[2][1], coords_3d[2][2]),
        Vector3::new(coords_3d[3][0], coords_3d[3][1], coords_3d[3][2]),
    ];
    let c = 0.25 * (x[0] + x[1] + x[2] + x[3]);

    // Sub-triangle area vectors, sign-aligned to the element normal.
    let mut tri = [Vector3::zeros(); 4];
    for k in 0..4 {
        let mut nk = (x[k] - c).cross(&(x[(k + 1) % 4] - c));
        if nk.dot(n_vec) < 0.0 {
            nk = -nk;
        }
        tri[k] = nk;
    }

    let mut vn = [Vector3::zeros(); 4];
    let mut v1 = [Vector3::zeros(); 4];
    let mut v2 = [Vector3::zeros(); 4];
    for i in 0..4 {
        let k_prev = (i + 3) % 4;
        let v = tri[k_prev] + tri[i];
        vn[i] = if v.norm() > 1e-30 {
            v.normalize()
        } else {
            *n_vec
        };

        let mut t1 = e1 - e1.dot(&vn[i]) * vn[i];
        if t1.norm() < 1e-12 {
            let alt = if vn[i][0].abs() < 0.9 {
                Vector3::new(1.0, 0.0, 0.0)
            } else {
                Vector3::new(0.0, 1.0, 0.0)
            };
            t1 = alt - alt.dot(&vn[i]) * vn[i];
        }
        v1[i] = t1.normalize();
        v2[i] = vn[i].cross(&v1[i]);
    }
    (vn, v1, v2)
}

// ============================================================================
// 2017 membrane coefficients
// ============================================================================

/// The five MITC4+ membrane coefficients of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (23)-(25) and (27a-c), pp. 409-410.
///
/// ```text
/// c_r = x_d . m^r,   c_s = x_d . m^s,   d = c_r^2 + c_s^2 - 1        (Eq. 25)
/// a_A = c_r(c_r-1)/(2d),  a_B = c_r(c_r+1)/(2d),
/// a_C = c_s(c_s-1)/(2d),  a_D = c_s(c_s+1)/(2d),
/// a_E = +2 c_r c_s / d                                               (Eq. 27c)
/// ```
///
/// `a_E` is **positive**, as printed on p. 410; the deleted implementation used
/// the negated form. Returns `(c_r, c_s, d, [a_A, a_B, a_C, a_D, a_E])`.
fn compute_membrane_coefficients_2017(
    x_d: &Vector3<f64>,
    m_r: &Vector3<f64>,
    m_s: &Vector3<f64>,
) -> (f64, f64, f64, [f64; 5]) {
    let c_r = x_d.dot(m_r);
    let c_s = x_d.dot(m_s);
    let d = c_r * c_r + c_s * c_s - 1.0;
    if d.abs() < 1e-14 {
        // Degenerate geometry (c_r^2 + c_s^2 = 1): the coefficients are
        // singular. Guard the division and return the flat limit.
        return (c_r, c_s, d, [0.0; 5]);
    }
    let a_a = c_r * (c_r - 1.0) / (2.0 * d);
    let a_b = c_r * (c_r + 1.0) / (2.0 * d);
    let a_c = c_s * (c_s - 1.0) / (2.0 * d);
    let a_d = c_s * (c_s + 1.0) / (2.0 * d);
    let a_e = 2.0 * c_r * c_s / d;
    (c_r, c_s, d, [a_a, a_b, a_c, a_d, a_e])
}

// ============================================================================
// Enriched 3D Jacobian
// ============================================================================

/// Enriched 3D base vectors `(g_r, g_s, g_t)` of Ko, Bathe & Zhang (2025),
/// C&S 308:107622, Eqs. (1), (4a-b), (5), pp. 2-3, and Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (1), (8a), pp. 405-406.
///
/// With `x(r,s,t) = sum h_i x_i + (t/2) sum a_i h_i V_n^i` (Eq. (1)),
///
/// ```text
/// g_r = sum (dh_i/dr) (x_i + (t/2) a_i V_n^i) = x_r + s x_d + (t/2)(...)
/// g_s = sum (dh_i/ds) (x_i + (t/2) a_i V_n^i) = x_s + r x_d + (t/2)(...)
/// g_t = (1/2) sum a_i h_i V_n^i = x_b
/// ```
///
/// so `j = det[g_r g_s g_t]|_{(r,s,0)} = (g_r x g_s) . g_t` and
/// `j0 = j(0,0,0)`, the ratio `j0/j` of Eq. (17b). The `t = 0` special case
/// recovers the paper's `g_r = x_r + s x_d`, `g_s = x_s + r x_d`.
fn compute_j3d_enriched(
    coords_3d: &[[f64; 3]; 4],
    vn: &[Vector3<f64>; 4],
    a_i: &[f64; 4],
    r: f64,
    s: f64,
    t: f64,
) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
    let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
    let h = shape_functions(r, s);
    let mut g_r = Vector3::zeros();
    let mut g_s = Vector3::zeros();
    let mut g_t = Vector3::zeros();
    for i in 0..4 {
        let x = Vector3::new(coords_3d[i][0], coords_3d[i][1], coords_3d[i][2]);
        let director = 0.5 * a_i[i] * vn[i];
        g_r += dn_dr[i] * (x + t * director);
        g_s += dn_ds[i] * (x + t * director);
        g_t += h[i] * director;
    }
    (g_r, g_s, g_t)
}

// ============================================================================
// Precomputed element data
// ============================================================================

/// All data that is constant for a given element geometry.
///
/// Every stored quantity names its source below. The 2017 core (Ko, Lee &
/// Bathe (2017), C&S 182:404-418) supplies the geometry and the membrane
/// coefficients; the 2025 drill element (Ko, Bathe & Zhang (2025),
/// C&S 308:107622) supplies the single drill normal `V^D` and the Jacobian
/// ratio `j0/j`.
#[derive(Clone)]
pub struct Mitc4PlusDPrecomputed {
    /// Local 2D coordinates of the four nodes, repository frame (4x2).
    pub local_coords: [[f64; 2]; 4],
    /// Local-to-global rotation, rows `(e1, e2, e3)`.
    pub t3: Matrix3<f64>,
    /// Local orthonormal frame, `e2 = e3 x e1`.
    pub e1: Vector3<f64>,
    /// Local orthonormal frame.
    pub e2: Vector3<f64>,
    /// Local orthonormal frame.
    pub e3: Vector3<f64>,
    /// Initial 3D node coordinates, used by the covariant operators.
    pub initial_coords_3d: [[f64; 3]; 4],
    /// Mid-surface constitutive (ADR-6).
    pub constitutive: ShellConstitutive,
    /// Shell thickness `h`; the 2017 `a_i` is the per-node form of this.
    pub thickness: f64,
    /// Per-node thickness `a_i` of Ko, Lee & Bathe (2017), C&S 182:404-418,
    /// Eq. (1), p. 405. No per-node thickness input exists in this change, so
    /// all four entries equal `thickness` (recorded limitation).
    pub a_i: [f64; 4],
    /// Per-node directors `V_n^i` of Ko, Lee & Bathe (2017), C&S 182:404-418,
    /// Eqs. (1)/(8a), pp. 405-406, built locally by `compute_node_directors`.
    pub vn: [Vector3<f64>; 4],
    /// Per-node in-plane directors `V_1^i` of Eq. (3d), p. 405.
    pub v1: [Vector3<f64>; 4],
    /// Per-node in-plane directors `V_2^i = V_n^i x V_1^i` of Eq. (3d).
    pub v2: [Vector3<f64>; 4],
    /// Characteristic vector `x_r`, Ko, Lee & Bathe (2017), C&S 182:404-418,
    /// Eq. (9), p. 406.
    pub x_r: Vector3<f64>,
    /// Characteristic vector `x_s`, Eq. (9).
    pub x_s: Vector3<f64>,
    /// Distortion vector `x_d`, Eq. (9).
    pub x_d: Vector3<f64>,
    /// Element plane normal `n = (x_r x x_s)/||x_r x x_s||`, Eq. (10), p. 406.
    pub n_vec: Vector3<f64>,
    /// Dual basis vector `m^r`, Eq. (11), p. 406.
    pub m_r: Vector3<f64>,
    /// Dual basis vector `m^s`, Eq. (11).
    pub m_s: Vector3<f64>,
    /// Paper A `c_r = x_d . m^r`, Ko, Lee & Bathe (2017), C&S 182:404-418,
    /// Eq. (25), p. 409. NOT the 2025 drill `c_r` of Eq. (18).
    pub c_r_mem: f64,
    /// Paper A `c_s = x_d . m^s`, Eq. (25). NOT the 2025 drill `c_s`.
    pub c_s_mem: f64,
    /// Paper A `d = c_r^2 + c_s^2 - 1`, Eq. (25).
    pub d_mem: f64,
    /// Paper A `[a_A, a_B, a_C, a_D, a_E]`, Eqs. (27a-c), p. 410. `a_E` is
    /// positive.
    pub a_coeffs: [f64; 5],
    /// Drill normal `V^D = (x_r x x_s)/||x_r x x_s||` at the centre, Ko, Bathe
    /// & Zhang (2025), C&S 308:107622, Eq. (5), p. 3. One vector per element.
    pub v_d: Vector3<f64>,
    /// `j0 = j(0,0,0)` with `j = det[g_r g_s g_t]|_{(r,s,0)}`, Ko, Bathe &
    /// Zhang (2025), C&S 308:107622, Eq. (17b) context, p. 8.
    pub j0: f64,
}

impl Mitc4PlusDPrecomputed {
    /// Build the precomputed element data from the 12 nodal coordinates
    /// `[x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4]`, the mid-surface constitutive
    /// and the shell thickness.
    ///
    /// The nodal directors are built locally (ADR-4 option B); no director or
    /// connectivity input is taken. The shear-correction channel of ADR-1 is
    /// added in WU4.
    pub fn new(node_coords: &[f64; 12], constitutive: ShellConstitutive, thickness: f64) -> Self {
        let mut coords_3d = [[0.0f64; 3]; 4];
        for i in 0..4 {
            coords_3d[i][0] = node_coords[3 * i];
            coords_3d[i][1] = node_coords[3 * i + 1];
            coords_3d[i][2] = node_coords[3 * i + 2];
        }

        let (local_coords, e1, e2, e3) = compute_local_coordinate_system(&coords_3d);
        let t3 = Matrix3::new(
            e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], e3[0], e3[1], e3[2],
        );

        let (x_r, x_s, x_d, n_vec, m_r, m_s) = compute_characteristic_vectors(&coords_3d);
        let a_i = [thickness; 4];
        let (vn, v1, v2) = compute_node_directors(&coords_3d, &n_vec, &e1);
        let (c_r_mem, c_s_mem, d_mem, a_coeffs) =
            compute_membrane_coefficients_2017(&x_d, &m_r, &m_s);

        // Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (5), p. 3: V^D is the
        // normal to the plane P at the element centre. At t = 0 the enriched
        // base vectors are g_r(0,0,0) = x_r and g_s(0,0,0) = x_s, so
        // V^D = normalize(x_r x x_s) = n_vec; it is stored once per element.
        let (g_r0, g_s0, g_t0) = compute_j3d_enriched(&coords_3d, &vn, &a_i, 0.0, 0.0, 0.0);
        let centre_normal = g_r0.cross(&g_s0);
        let v_d = if centre_normal.norm() > 1e-14 {
            centre_normal.normalize()
        } else {
            n_vec
        };
        let j0 = centre_normal.dot(&g_t0);

        Mitc4PlusDPrecomputed {
            local_coords,
            t3,
            e1,
            e2,
            e3,
            initial_coords_3d: coords_3d,
            constitutive,
            thickness,
            a_i,
            vn,
            v1,
            v2,
            x_r,
            x_s,
            x_d,
            n_vec,
            m_r,
            m_s,
            c_r_mem,
            c_s_mem,
            d_mem,
            a_coeffs,
            v_d,
            j0,
        }
    }
}

// ============================================================================
// Kinematics interpolation
// ============================================================================

/// Position interpolation of Ko, Bathe & Zhang (2025), C&S 308:107622,
/// Eq. (1), p. 2, equivalently Ko, Lee & Bathe (2017), C&S 182:404-418,
/// Eq. (1), p. 405:
///
/// ```text
/// x(r,s,t) = sum h_i x_i + (t/2) sum a_i h_i V_n^i,   t in [-1, 1]
/// ```
fn interpolate_position(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64, t: f64) -> Vector3<f64> {
    let h = shape_functions(r, s);
    let mut x = Vector3::zeros();
    for i in 0..4 {
        let xi = Vector3::new(
            pre.initial_coords_3d[i][0],
            pre.initial_coords_3d[i][1],
            pre.initial_coords_3d[i][2],
        );
        x += h[i] * xi;
        x += (t * 0.5) * pre.a_i[i] * h[i] * pre.vn[i];
    }
    x
}

/// Displacement interpolation of Ko, Bathe & Zhang (2025), C&S 308:107622,
/// Eq. (3a), p. 2, equivalently Ko, Lee & Bathe (2017), C&S 182:404-418,
/// Eq. (3), p. 405:
///
/// ```text
/// u(r,s,t) = sum h_i u_i + (t/2) sum a_i h_i (theta_i x V_n^i)
/// ```
///
/// `dofs` is the repository's 24-vector, six per node in the local frame:
/// `(u, v, w, theta_x, theta_y, theta_z)` at slots `6i .. 6i+5`. The rotation
/// enters only through `theta_i x V_n^i`, so the component of `theta_i` along
/// `V_n^i` is annihilated: `theta_z` is carried in the DOF vector but is
/// **unconsumed by the 2017 core**. The drill rotation enters through the 2025
/// operator (WU3), where it is `theta_i^D = theta_i . V^D`.
fn interpolate_displacement(
    pre: &Mitc4PlusDPrecomputed,
    r: f64,
    s: f64,
    t: f64,
    dofs: &[f64; 24],
) -> Vector3<f64> {
    let h = shape_functions(r, s);
    let mut u = Vector3::zeros();
    for i in 0..4 {
        let ui = Vector3::new(dofs[6 * i], dofs[6 * i + 1], dofs[6 * i + 2]);
        let theta = Vector3::new(dofs[6 * i + 3], dofs[6 * i + 4], dofs[6 * i + 5]);
        u += h[i] * ui;
        u += (t * 0.5) * pre.a_i[i] * h[i] * theta.cross(&pre.vn[i]);
    }
    u
}

#[cfg(test)]
mod tests {
    use super::{
        compute_membrane_coefficients_2017, interpolate_displacement, interpolate_position,
        shape_functions, Mitc4PlusDPrecomputed,
    };
    use crate::materials::{isotropic::IsotropicMaterial, Material, ShellConstitutive};
    use nalgebra::{DMatrix, Vector3};

    // ========================================================================
    // The 2017 and 2025 boundary-condition fixtures — the star patch
    // (Ko, Lee & Bathe (2017), C&S 182:404–418, Fig. 5 = Ko, Bathe & Zhang
    //  (2025), C&S 308:107622, Fig. 7(a))
    // ========================================================================
    //
    // Source: `docs/formulations/mitc4plusd-2025-extract.md` §Fig. 7(a), a
    // vision read of Ko, Bathe & Zhang (2025), C&S 308:107622, Fig. 7(a), p. 5.
    // The mesh is NOT a 3x3 grid: it is a five-element star with eight nodes on
    // the 10x10 square.
    //
    //   corners  B(0,10)  A(10,10)  D(10,0)  C(0,0)
    //   interior (4,7) (8,7) (8,3) (2,2)
    //
    // Connectivity is written in the repository's counter-clockwise (r,s)
    // convention (node 0 at (-1,-1), 1 at (+1,-1), 2 at (+1,+1), 3 at (-1,+1));
    // the extract's prose order is clockwise, so the central element is listed
    // here as (2,2) -> (8,3) -> (8,7) -> (4,7).

    /// The eight star-patch nodes `[x, y, z]`. Nodes 0..3 are the corners
    /// C, D, A, B; nodes 4..7 are the interior nodes (2,2), (8,3), (8,7),
    /// (4,7). Every coordinate is a figure read (2025 extract §Fig. 7(a)).
    pub const STAR_NODES: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],   // 0  C
        [10.0, 0.0, 0.0],  // 1  D
        [10.0, 10.0, 0.0], // 2  A
        [0.0, 10.0, 0.0],  // 3  B
        [2.0, 2.0, 0.0],   // 4
        [8.0, 3.0, 0.0],   // 5
        [8.0, 7.0, 0.0],   // 6
        [4.0, 7.0, 0.0],   // 7
    ];

    /// The five star-patch elements, counter-clockwise.
    pub const STAR_ELEMS: [[usize; 4]; 5] = [
        [4, 5, 6, 7], // central (2,2) (8,3) (8,7) (4,7)
        [7, 6, 2, 3], // top     (4,7) (8,7) A     B
        [6, 5, 1, 2], // right   (8,7) (8,3) D     A
        [0, 1, 5, 4], // bottom  C     D     (8,3) (2,2)
        [4, 7, 3, 0], // left    (2,2) (4,7) B     C
    ];

    /// Signed area (shoelace) of a quadrilateral in the xy-plane; positive
    /// for a counter-clockwise node order.
    fn signed_area(nodes: &[[f64; 3]; 8], elem: &[usize; 4]) -> f64 {
        let mut a = 0.0;
        for k in 0..4 {
            let p = nodes[elem[k]];
            let q = nodes[elem[(k + 1) % 4]];
            a += p[0] * q[1] - q[0] * p[1];
        }
        0.5 * a
    }

    #[test]
    fn test_star_patch_has_five_elements_and_eight_nodes() {
        assert_eq!(STAR_NODES.len(), 8, "star patch must have 8 nodes");
        assert_eq!(STAR_ELEMS.len(), 5, "star patch must have 5 elements");
        for e in &STAR_ELEMS {
            for &n in e {
                assert!(
                    n < STAR_NODES.len(),
                    "element references node {n} out of range"
                );
            }
        }
    }

    #[test]
    fn test_star_patch_all_signed_areas_positive_and_sum_to_100() {
        let areas: Vec<f64> = STAR_ELEMS
            .iter()
            .map(|e| signed_area(&STAR_NODES, e))
            .collect();
        for (i, a) in areas.iter().enumerate() {
            assert!(
                *a > 0.0,
                "element {i} is not counter-clockwise: signed area {a}"
            );
        }
        let total: f64 = areas.iter().sum();
        assert!(
            (total - 100.0).abs() < 1e-12,
            "areas sum to {total}, expected 100"
        );
    }

    #[test]
    fn test_star_patch_has_no_duplicate_coordinates() {
        for (i, a) in STAR_NODES.iter().enumerate() {
            for (j, b) in STAR_NODES.iter().enumerate().skip(i + 1) {
                let d: f64 = (0..3).map(|k| (a[k] - b[k]).powi(2)).sum();
                assert!(d > 1e-24, "nodes {i} and {j} share coordinates");
            }
        }
    }

    #[test]
    fn test_star_patch_central_element_matches_extract() {
        // The extract's central element is (2,2) -> (4,7) -> (8,7) -> (8,3),
        // i.e. clockwise; the fixture stores the counter-clockwise rotation
        // (2,2) -> (8,3) -> (8,7) -> (4,7) required by the element's (r,s)
        // convention. The vertex set is the extract's.
        let expected = [[2.0, 2.0], [8.0, 3.0], [8.0, 7.0], [4.0, 7.0]];
        for (k, e) in STAR_ELEMS[0].iter().enumerate() {
            let got = STAR_NODES[*e];
            assert_eq!([got[0], got[1]], expected[k], "central element vertex {k}");
        }
    }

    // ========================================================================
    // Boundary-condition sets
    // ========================================================================

    /// One scalar boundary condition: DOF `dof` (0..2 = u_x,u_y,u_z;
    /// 3..5 = θ_x,θ_y,θ_z) at patch node `node` is prescribed to `value`.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Bc {
        pub node: usize,
        pub dof: usize,
        pub value: f64,
    }

    const fn bc(node: usize, dof: usize) -> Bc {
        Bc {
            node,
            dof,
            value: 0.0,
        }
    }

    /// BC_2017_PATCH — the T1.3 minimum boundary conditions for the 2017 patch
    /// tests of Ko, Lee & Bathe (2017), C&S 182:404–418.
    ///
    /// PROVENANCE: **derived**, not figure-read. Ko, Lee & Bathe (2017),
    /// C&S 182:404–418, Fig. 5 publishes only the mesh; §4 (p. 410) says the
    /// patch is subjected to "the minimum number of constraints to prevent
    /// rigid body motions" and does not print them.
    ///
    /// DEVIATION FROM DESIGN §4.1 (recorded): the design lists `u_x=u_y=u_z=0`
    /// at C, `u_x=0` at D, `u_y=0` at A and `θ_z=0` at C. That set does not
    /// remove the six rigid-body modes: C=(0,0) and D=(10,0) share `y=0`, so
    /// their `u_x=0` rows project onto the *same* rigid-body equation
    /// (`u_x = t_x - ω_z y`), and neither `θ_x` nor `θ_y` is constrained, so the
    /// two tilt modes survive (rank 4 on the six rigid-body modes). The set
    /// below removes exactly six independent rigid-body constraints:
    ///   C(0): u_x=u_y=u_z=0   (the three translations)
    ///   A(2): u_y=0           (with C.u_y, removes the drill ω_z)
    ///   C(0): θ_x=θ_y=0       (the two tilt modes)
    /// `θ_z` is left free: the 2017 core is blind to it. The prescribed values
    /// of the constant-straining case are supplied at solve time by the test
    /// that owns the mode; here they are stored as zero.
    pub const BC_2017_PATCH: &[Bc] = &[
        bc(0, 0),
        bc(0, 1),
        bc(0, 2), // C: u_x, u_y, u_z
        bc(2, 1), // A: u_y
        bc(0, 3),
        bc(0, 4), // C: θ_x, θ_y
    ];

    /// BC_2025_PATCH — the T1.6 boundary conditions for the strong-form patch
    /// tests of Ko, Bathe & Zhang (2025), C&S 308:107622 (Fig. 7(b)(c)(d),
    /// p. 5).
    ///
    /// PROVENANCE: **figure-read**, recorded verbatim in
    /// `docs/formulations/mitc4plusd-2025-extract.md` §Fig. 7(b)(c)(d). Every
    /// prescribed value is zero. In all three cases `θ_z` is free at every node
    /// except corner B.
    pub struct Patch2025Bc {
        pub extension: &'static [Bc],
        pub bending: &'static [Bc],
        pub shearing: &'static [Bc],
    }

    pub const BC_2025_PATCH: Patch2025Bc = Patch2025Bc {
        // (d) extension: B fully clamped; C: u_x=u_z=0; load at A in +x.
        extension: &[
            bc(3, 0),
            bc(3, 1),
            bc(3, 2),
            bc(3, 3),
            bc(3, 4),
            bc(3, 5),
            bc(0, 0),
            bc(0, 2),
        ],
        // (b) bending: B fully clamped; C: u_x=u_z=0, θ_y=0; load at A in +y.
        bending: &[
            bc(3, 0),
            bc(3, 1),
            bc(3, 2),
            bc(3, 3),
            bc(3, 4),
            bc(3, 5),
            bc(0, 0),
            bc(0, 2),
            bc(0, 4),
        ],
        // (c) shearing: B fully clamped; C: u_x=u_z=0, θ_x=θ_y=0;
        // interior {4,5,6,7}: u_x=0, θ_x=θ_y=0; load at A in +y.
        shearing: &[
            bc(3, 0),
            bc(3, 1),
            bc(3, 2),
            bc(3, 3),
            bc(3, 4),
            bc(3, 5),
            bc(0, 0),
            bc(0, 2),
            bc(0, 3),
            bc(0, 4),
            bc(4, 0),
            bc(4, 3),
            bc(4, 4),
            bc(5, 0),
            bc(5, 3),
            bc(5, 4),
            bc(6, 0),
            bc(6, 3),
            bc(6, 4),
            bc(7, 0),
            bc(7, 3),
            bc(7, 4),
        ],
    };

    /// The six rigid-body fields of the flat star patch in the 6-DOF/node
    /// layout: `u_i = t + ω × (x_i − x_c)`, `θ_i = ω`, columns ordered
    /// (t_x, t_y, t_z, ω_x, ω_y, ω_z).
    fn rigid_body_modes() -> DMatrix<f64> {
        let n = STAR_NODES.len();
        let mut c = Vector3::new(0.0, 0.0, 0.0);
        for x in STAR_NODES.iter() {
            c += Vector3::new(x[0], x[1], x[2]);
        }
        c /= n as f64;

        let unit = [
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut m = DMatrix::zeros(6 * n, 6);
        for (i, x) in STAR_NODES.iter().enumerate() {
            let r = Vector3::new(x[0], x[1], x[2]) - c;
            for k in 0..3 {
                for d in 0..3 {
                    m[(6 * i + d, k)] = unit[k][d]; // translation t_k
                    m[(6 * i + 3 + d, 3 + k)] = unit[k][d]; // rotation θ_k
                }
                let u = unit[k].cross(&r);
                for d in 0..3 {
                    m[(6 * i + d, 3 + k)] = u[d];
                }
            }
        }
        m
    }

    /// Selection matrix of a boundary-condition set: one row per constraint.
    fn constraint_matrix(bcs: &[Bc]) -> DMatrix<f64> {
        let mut c = DMatrix::zeros(bcs.len(), 6 * STAR_NODES.len());
        for (row, b) in bcs.iter().enumerate() {
            c[(row, 6 * b.node + b.dof)] = 1.0;
        }
        c
    }

    fn rank(m: &DMatrix<f64>) -> usize {
        let svd = m.clone().svd(false, false);
        let tol = 1e-9 * svd.singular_values.max().max(1.0);
        svd.singular_values.iter().filter(|&&s| s > tol).count()
    }

    #[test]
    fn test_bc_2017_patch_constrains_exactly_six_rigid_body_modes() {
        let crb = constraint_matrix(BC_2017_PATCH) * rigid_body_modes();
        assert_eq!(crb.ncols(), 6);
        assert_eq!(
            rank(&crb),
            6,
            "BC_2017_PATCH must remove exactly the six rigid-body modes"
        );
    }

    #[test]
    fn test_bc_2025_patch_constrains_exactly_six_rigid_body_modes() {
        for (name, set) in [
            ("extension", BC_2025_PATCH.extension),
            ("bending", BC_2025_PATCH.bending),
            ("shearing", BC_2025_PATCH.shearing),
        ] {
            let crb = constraint_matrix(set) * rigid_body_modes();
            assert_eq!(
                rank(&crb),
                6,
                "BC_2025_PATCH {name} must remove the six rigid-body modes"
            );
        }
    }

    #[test]
    fn test_bc_2025_patch_matches_fig7_bcd() {
        // Expected (node, dof) triples from Ko, Bathe & Zhang (2025),
        // C&S 308:107622, Fig. 7(b)(c)(d) (2025 extract §Fig. 7(b)(c)(d)).
        // dof: 0=u_x, 1=u_y, 2=u_z, 3=θ_x, 4=θ_y, 5=θ_z.
        let clamped_b = [(3, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5)];
        let mut expected_extension = clamped_b.to_vec();
        expected_extension.extend_from_slice(&[(0, 0), (0, 2)]);
        let mut expected_bending = clamped_b.to_vec();
        expected_bending.extend_from_slice(&[(0, 0), (0, 2), (0, 4)]);
        let mut expected_shearing = clamped_b.to_vec();
        expected_shearing.extend_from_slice(&[(0, 0), (0, 2), (0, 3), (0, 4)]);
        for n in [4, 5, 6, 7] {
            expected_shearing.extend_from_slice(&[(n, 0), (n, 3), (n, 4)]);
        }

        let actual = |set: &[Bc]| set.iter().map(|b| (b.node, b.dof)).collect::<Vec<_>>();
        assert_eq!(
            actual(BC_2025_PATCH.extension),
            expected_extension,
            "extension (Fig. 7d)"
        );
        assert_eq!(
            actual(BC_2025_PATCH.bending),
            expected_bending,
            "bending (Fig. 7b)"
        );
        assert_eq!(
            actual(BC_2025_PATCH.shearing),
            expected_shearing,
            "shearing (Fig. 7c)"
        );

        // θ_z (dof 5) is constrained only at corner B, in every case, and every
        // figure-read value is zero.
        for (name, set) in [
            ("extension", BC_2025_PATCH.extension),
            ("bending", BC_2025_PATCH.bending),
            ("shearing", BC_2025_PATCH.shearing),
        ] {
            let drill: Vec<usize> = set.iter().filter(|b| b.dof == 5).map(|b| b.node).collect();
            assert_eq!(
                drill,
                vec![3],
                "{name}: θ_z must be free except at corner B"
            );
            assert!(
                set.iter().all(|b| b.value == 0.0),
                "{name}: every Fig. 7 value is zero"
            );
        }
    }

    // ========================================================================
    // F-W — the derived warped star patch
    // ========================================================================

    /// F-W — the warped star patch.
    ///
    /// PROVENANCE: **derived**. The same eight nodes as the flat 2017/2025 star
    /// patch with `z` offsets on the four interior nodes only, `±0.5`
    /// alternating. It is used only for the warped-variant separation of T1.6d
    /// (`θ_z` constrained at every node vs free) and for the 2×2×2 vs
    /// surface-only integration discrimination. A warped patch is NOT a valid
    /// constant-stress patch test and F-W is never used for a constant-stress
    /// assertion.
    pub const STAR_NODES_WARPED: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],   // 0  C
        [10.0, 0.0, 0.0],  // 1  D
        [10.0, 10.0, 0.0], // 2  A
        [0.0, 10.0, 0.0],  // 3  B
        [2.0, 2.0, 0.5],   // 4
        [8.0, 3.0, -0.5],  // 5
        [8.0, 7.0, 0.5],   // 6
        [4.0, 7.0, -0.5],  // 7
    ];

    /// The four sub-quad normals of an element, split by its centre and
    /// sign-aligned to `e3`. This is the element-local normal construction the
    /// design's `compute_node_directors` averages per node.
    fn sub_quad_normals(nodes: &[[f64; 3]; 8], elem: &[usize; 4]) -> [Vector3<f64>; 4] {
        let v = |n: usize| {
            let x = nodes[n];
            Vector3::new(x[0], x[1], x[2])
        };
        let c = (v(elem[0]) + v(elem[1]) + v(elem[2]) + v(elem[3])) / 4.0;

        let mut out = [Vector3::zeros(); 4];
        for k in 0..4 {
            let mut n = (v(elem[(k + 1) % 4]) - c).cross(&(v(elem[(k + 2) % 4]) - v(elem[k])));
            if n.z < 0.0 {
                n = -n;
            }
            out[k] = n.normalize();
        }
        out
    }

    #[test]
    fn test_f_w_z_offsets_on_interior_nodes_only() {
        for (i, n) in STAR_NODES_WARPED.iter().take(4).enumerate() {
            assert_eq!(n[2], 0.0, "corner node {i} must stay flat");
        }
        for (i, n) in STAR_NODES_WARPED.iter().enumerate().skip(4) {
            assert!(n[2].abs() > 0.0, "interior node {i} needs a z offset");
        }
        for (i, (fw, flat)) in STAR_NODES_WARPED.iter().zip(STAR_NODES.iter()).enumerate() {
            assert_eq!(
                [fw[0], fw[1]],
                [flat[0], flat[1]],
                "F-W node {i} xy changed"
            );
        }
    }

    #[test]
    fn test_f_w_element_normals_differ_from_flat_n_vec() {
        let e3 = Vector3::new(0.0, 0.0, 1.0);
        // Flat patch: every sub-quad normal is exactly the flat n_vec.
        for e in &STAR_ELEMS {
            for n in sub_quad_normals(&STAR_NODES, e) {
                assert!(
                    (n - e3).norm() < 1e-14,
                    "flat sub-quad normal {n} differs from n_vec"
                );
            }
        }
        // Warped patch: every element has a sub-quad normal that deviates from
        // the flat n_vec by a measurable angle.
        for (i, e) in STAR_ELEMS.iter().enumerate() {
            let max_dev = sub_quad_normals(&STAR_NODES_WARPED, e)
                .iter()
                .map(|n| n.dot(&e3).clamp(-1.0, 1.0).acos())
                .fold(0.0_f64, f64::max);
            assert!(
                max_dev > 1e-3,
                "warped element {i}: max normal deviation {max_dev}"
            );
        }
    }

    #[test]
    fn test_f_w_is_not_a_constant_stress_fixture() {
        // A constant-stress patch requires the mesh to carry the constant
        // state; F-W is deliberately warped, so it must never be used for a
        // constant-stress assertion. Assert the property that makes that true:
        // F-W is genuinely non-planar, i.e. no element is flat.
        assert!(
            STAR_NODES_WARPED.iter().any(|n| n[2] != 0.0),
            "F-W must be warped"
        );
        for (i, e) in STAR_ELEMS.iter().enumerate() {
            let z0 = STAR_NODES_WARPED[e[0]][2];
            assert!(
                e.iter().any(|&n| STAR_NODES_WARPED[n][2] != z0),
                "F-W element {i} is flat"
            );
        }
    }

    // ========================================================================
    // WU2 — element core: geometry, coefficients, directors and kinematics
    // (tasks 3.1-3.4).  Tests are written before the production code; each is
    // shown RED by a deliberate perturbation of the implementation.
    // ========================================================================

    /// Flatten a `[[f64;3];4]` node list into the constructor's `[f64;12]`.
    fn coords12(c: &[[f64; 3]; 4]) -> [f64; 12] {
        let mut n = [0.0f64; 12];
        for i in 0..4 {
            for j in 0..3 {
                n[3 * i + j] = c[i][j];
            }
        }
        n
    }

    /// A real isotropic shell constitutive, so the geometry tests exercise the
    /// same constructor path the assembly layers will.
    fn shell_iso() -> ShellConstitutive {
        IsotropicMaterial::new(2.0e11, 0.3, 7800.0).constitutive(1.0, 5.0 / 6.0)
    }

    fn pre_from(c: &[[f64; 3]; 4]) -> Mitc4PlusDPrecomputed {
        Mitc4PlusDPrecomputed::new(&coords12(c), shell_iso(), 1.0)
    }

    /// A flat rectangle in the repository's `(r,s)` node order.
    const RECT: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];

    /// A flat, distorted (non-parallelogram) quad: `x_d` is in-plane and non-zero.
    const FLAT_DISTORTED: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [2.5, 0.4, 0.0],
        [2.1, 1.8, 0.0],
        [0.2, 1.2, 0.0],
    ];

    /// A ruled warped quad.
    const RULED_WARPED: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.4],
        [2.0, 1.0, -0.4],
        [0.0, 1.0, 0.0],
    ];

    /// A doubly warped quad (all four nodes off the plane).
    const DOUBLY_WARPED: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.1],
        [2.0, 0.0, 0.4],
        [2.0, 1.0, -0.3],
        [0.0, 1.0, -0.2],
    ];

    #[test]
    fn test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18() {
        let pre = pre_from(&RECT);

        // Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (9): x_d vanishes for a
        // flat rectangle, so Eq. (23)-(25)'s c_r, c_s vanish and d = -1.
        let scale = pre.x_r.norm().max(pre.x_s.norm());
        assert!(
            pre.x_d.norm() <= 1e-14 * scale,
            "x_d = {} (scale {scale})",
            pre.x_d
        );
        assert!(pre.c_r_mem.abs() <= 1e-14, "c_r = {}", pre.c_r_mem);
        assert!(pre.c_s_mem.abs() <= 1e-14, "c_s = {}", pre.c_s_mem);
        assert!((pre.d_mem + 1.0).abs() <= 1e-14, "d = {}", pre.d_mem);
        for (k, a) in pre.a_coeffs.iter().enumerate() {
            assert!(a.abs() <= 1e-14, "a[{k}] = {a}");
        }

        // With all five coefficients zero, Eq. (27a-c) must collapse term by
        // term onto the linear assumed field of Eqs. (18)/(19).  The two sides
        // are evaluated by different formulas, so a wrong coefficient or a
        // wrong polynomial fails even though the coefficients are zero.
        let (r, s) = (0.37, -0.62);
        let (e_a, e_b, e_c, e_d, e_e) = (1.3f64, -0.7f64, 2.1f64, 0.4f64, -1.9f64);
        let a = pre.a_coeffs;

        let rr27 = 0.5 * (1.0 - 2.0 * a[0] + s + 2.0 * a[0] * s * s) * e_a
            + 0.5 * (1.0 - 2.0 * a[1] - s + 2.0 * a[1] * s * s) * e_b
            + a[2] * (-1.0 + s * s) * e_c
            + a[3] * (-1.0 + s * s) * e_d
            + a[4] * (-1.0 + s * s) * e_e;
        let rr18 = 0.5 * (1.0 + s) * e_a + 0.5 * (1.0 - s) * e_b;
        assert!(
            (rr27 - rr18).abs() <= 1e-13,
            "Eq. (27a) {rr27} vs Eq. (18a) {rr18}"
        );

        let ss27 = a[0] * (-1.0 + r * r) * e_a
            + a[1] * (-1.0 + r * r) * e_b
            + 0.5 * (1.0 - 2.0 * a[2] + r + 2.0 * a[2] * r * r) * e_c
            + 0.5 * (1.0 - 2.0 * a[3] - r + 2.0 * a[3] * r * r) * e_d
            + a[4] * (-1.0 + r * r) * e_e;
        let ss18 = 0.5 * (1.0 + r) * e_c + 0.5 * (1.0 - r) * e_d;
        assert!(
            (ss27 - ss18).abs() <= 1e-13,
            "Eq. (27b) {ss27} vs Eq. (18b) {ss18}"
        );

        // Eq. (27c) collapses onto Eq. (19c): the bilinear shear plus the two
        // linear shear terms that the patch test requires.
        let rs27 = 0.25 * (r + 4.0 * a[0] * r * s) * e_a
            + 0.25 * (-r + 4.0 * a[1] * r * s) * e_b
            + 0.25 * (s + 4.0 * a[2] * r * s) * e_c
            + 0.25 * (-s + 4.0 * a[3] * r * s) * e_d
            + (1.0 + a[4] * r * s) * e_e;
        let e_rr_lin = 0.5 * (e_a - e_b); // Eq. (17): e(A) - e(B) = 2 e_lin
        let e_ss_lin = 0.5 * (e_c - e_d);
        let rs19 = e_e + 0.5 * e_rr_lin * r + 0.5 * e_ss_lin * s;
        assert!(
            (rs27 - rs19).abs() <= 1e-13,
            "Eq. (27c) {rs27} vs Eq. (19c) {rs19}"
        );
    }

    #[test]
    fn test_geometry_dual_basis_identities_eq11() {
        // Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (11): the dual basis
        // m^r, m^s on the plane P satisfies m^r.x_r = m^s.x_s = 1,
        // m^r.x_s = m^s.x_r = 0 and m^r.n = m^s.n = 0.
        for (name, c) in [
            ("flat-distorted", FLAT_DISTORTED),
            ("ruled-warped", RULED_WARPED),
            ("doubly-warped", DOUBLY_WARPED),
        ] {
            let pre = pre_from(&c);
            assert!(
                (pre.m_r.dot(&pre.x_r) - 1.0).abs() <= 1e-12,
                "{name}: m^r.x_r = {}",
                pre.m_r.dot(&pre.x_r)
            );
            assert!(
                (pre.m_s.dot(&pre.x_s) - 1.0).abs() <= 1e-12,
                "{name}: m^s.x_s = {}",
                pre.m_s.dot(&pre.x_s)
            );
            assert!(
                pre.m_r.dot(&pre.x_s).abs() <= 1e-12 * pre.m_r.norm() * pre.x_s.norm(),
                "{name}: m^r.x_s = {}",
                pre.m_r.dot(&pre.x_s)
            );
            assert!(
                pre.m_s.dot(&pre.x_r).abs() <= 1e-12 * pre.m_s.norm() * pre.x_r.norm(),
                "{name}: m^s.x_r = {}",
                pre.m_s.dot(&pre.x_r)
            );
            assert!(
                pre.m_r.dot(&pre.n_vec).abs() <= 1e-12 * pre.m_r.norm(),
                "{name}: m^r.n = {}",
                pre.m_r.dot(&pre.n_vec)
            );
            assert!(
                pre.m_s.dot(&pre.n_vec).abs() <= 1e-12 * pre.m_s.norm(),
                "{name}: m^s.n = {}",
                pre.m_s.dot(&pre.n_vec)
            );
        }
    }

    #[test]
    fn test_geometry_a_E_is_positive_eq27c() {
        // Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (27c), p. 410: the paper
        // prints a_E = +2 c_r c_s / d, POSITIVE.  The deleted implementation
        // used the negated form; that sign error is what this test rejects.
        //
        // RECORDED DEVIATION: for a convex element c_r^2 + c_s^2 > 1 (i.e.
        // d > 0), which is what makes +2 c_r c_s / d positive for same-sign
        // c_r, c_s, requires an extreme warped distortion; no flat quad in the
        // repository's fixtures reaches it.  The formula is therefore exercised
        // directly on its own inputs (c_r, c_s > 0, d > 0), and the same closed
        // form is re-checked on a real flat distorted element.
        let x_d = Vector3::new(1.0, 1.0, 0.0);
        let m_r = Vector3::new(0.8, 0.0, 0.0);
        let m_s = Vector3::new(0.0, 0.8, 0.0);
        let (c_r, c_s, d, a) = compute_membrane_coefficients_2017(&x_d, &m_r, &m_s);
        assert!(
            c_r * c_s > 0.0 && d > 0.0,
            "c_r = {c_r}, c_s = {c_s}, d = {d}"
        );
        let a_e = a[4];
        let expected = 2.0 * c_r * c_s / d;
        assert!(a_e > 0.0, "a_E must be positive, got {a_e}");
        assert!(
            (a_e - expected).abs() <= 1e-14 * expected.abs(),
            "a_E = {a_e}, +2 c_r c_s / d = {expected}"
        );
        assert!(
            (a_e - (-expected)).abs() > 1e-6 * expected.abs(),
            "the negated (deleted) form must be rejected"
        );

        // On a real flat distorted element the stored a_E matches the same
        // closed form and the negated form is still rejected.
        let pre = pre_from(&FLAT_DISTORTED);
        let a_e = pre.a_coeffs[4];
        let expected = 2.0 * pre.c_r_mem * pre.c_s_mem / pre.d_mem;
        assert!(
            (a_e - expected).abs() <= 1e-14 * expected.abs(),
            "element a_E = {a_e}, +2 c_r c_s / d = {expected}"
        );
        assert!(
            (a_e - (-expected)).abs() > 1e-6 * expected.abs(),
            "element: the negated (deleted) form must be rejected"
        );
    }

    #[test]
    fn test_geometry_node_directors_reduce_to_n_vec_when_flat() {
        // ADR-4 option B: the per-node director V_n^i is the area-weighted
        // average of the two centre sub-quad normals adjacent to node i,
        // sign-aligned to n_vec.  On flat geometry every sub-quad normal is
        // n_vec, so every director is n_vec.
        let flat_square = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        for (name, c) in [
            ("flat square", flat_square),
            ("flat distorted", FLAT_DISTORTED),
        ] {
            let pre = pre_from(&c);
            for i in 0..4 {
                assert!(
                    (pre.vn[i] - pre.n_vec).norm() <= 1e-14,
                    "{name}: V_n^{i} = {} != n_vec {}",
                    pre.vn[i],
                    pre.n_vec
                );
                assert!(
                    pre.v1[i].dot(&pre.vn[i]).abs() <= 1e-14,
                    "{name}: V_1^{i} not normal to V_n"
                );
                assert!(
                    pre.v2[i].dot(&pre.vn[i]).abs() <= 1e-14,
                    "{name}: V_2^{i} not normal to V_n"
                );
                assert!(
                    (pre.v1[i].cross(&pre.v2[i]) - pre.vn[i]).norm() <= 1e-12,
                    "{name}: (V_1 x V_2)^{i} != V_n^{i} (not right-handed)"
                );
            }
        }

        // Non-vacuity: on the warped fixture the four directors must be
        // measurably different, otherwise the flat assertion above is empty.
        let e = STAR_ELEMS[0];
        let mut wc = [[0.0f64; 3]; 4];
        for k in 0..4 {
            wc[k] = STAR_NODES_WARPED[e[k]];
        }
        let pre = pre_from(&wc);
        let mut max_diff = 0.0f64;
        for i in 0..4 {
            for j in 0..4 {
                max_diff = max_diff.max((pre.vn[i] - pre.vn[j]).norm());
            }
        }
        assert!(
            max_diff > 1e-6,
            "directors are indistinguishable on F-W (max pairwise diff {max_diff})"
        );
    }

    #[test]
    fn test_kinematics_displacement_field_matches_eq1_to_eq3() {
        // Ko, Bathe & Zhang (2025), C&S 308:107622, Eqs. (1)/(3a), p. 2:
        //   x(r,s,t) = sum h_i x_i + (t/2) sum a_i h_i V_n^i
        //   u(r,s,t) = sum h_i u_i + (t/2) sum a_i h_i (theta_i x V_n^i)
        let pre = pre_from(&DOUBLY_WARPED);

        // A generic 24-vector: (u,v,w, theta_x,theta_y,theta_z) per node.
        let mut dofs = [0.0f64; 24];
        for i in 0..4 {
            for k in 0..6 {
                dofs[6 * i + k] = 0.3 + 0.17 * (6 * i + k) as f64;
            }
        }

        // Eq. (1): the position interpolation, evaluated independently from the
        // stored directors and thicknesses.
        for &(r, s, t) in &[
            (-0.5, 0.3, 0.0),
            (0.4, -0.7, 1.0),
            (0.2, 0.5, -1.0),
            (0.0, 0.0, 1.0),
        ] {
            let h = shape_functions(r, s);
            let mut x_ref = Vector3::zeros();
            for i in 0..4 {
                let xi = Vector3::new(
                    pre.initial_coords_3d[i][0],
                    pre.initial_coords_3d[i][1],
                    pre.initial_coords_3d[i][2],
                );
                x_ref += h[i] * xi + (t * 0.5) * pre.a_i[i] * h[i] * pre.vn[i];
            }
            let x = interpolate_position(&pre, r, s, t);
            assert!(
                (x - x_ref).norm() <= 1e-12 * x_ref.norm(),
                "x({r},{s},{t}) = {x}, Eq. (1) = {x_ref}"
            );
        }

        // Eq. (3d) + the rotation identity of Eq. (3a): for a right-handed
        // (V_1, V_2, V_n), theta x V_n = -V_2 alpha + V_1 beta with
        // alpha = theta.V_1 and beta = theta.V_2.
        for i in 0..4 {
            let theta = Vector3::new(dofs[6 * i + 3], dofs[6 * i + 4], dofs[6 * i + 5]);
            let alpha = theta.dot(&pre.v1[i]);
            let beta = theta.dot(&pre.v2[i]);
            let lhs = theta.cross(&pre.vn[i]);
            let rhs = -pre.v2[i] * alpha + pre.v1[i] * beta;
            assert!(
                (lhs - rhs).norm() <= 1e-12 * theta.norm(),
                "node {i}: theta x V_n {lhs} != -V_2 a + V_1 b {rhs}"
            );
        }

        // Eq. (3a): u at t = 0 is the membrane interpolation, and the t-linear
        // part is exactly (t/2) sum a_i h_i (theta_i x V_n^i).
        for &(r, s) in &[(0.3, -0.2), (-0.6, 0.4)] {
            let h = shape_functions(r, s);
            let mut um = Vector3::zeros();
            let mut ub = Vector3::zeros();
            for i in 0..4 {
                let ui = Vector3::new(dofs[6 * i], dofs[6 * i + 1], dofs[6 * i + 2]);
                let theta = Vector3::new(dofs[6 * i + 3], dofs[6 * i + 4], dofs[6 * i + 5]);
                um += h[i] * ui;
                ub += pre.a_i[i] * h[i] * theta.cross(&pre.vn[i]);
            }
            let u0 = interpolate_displacement(&pre, r, s, 0.0, &dofs);
            assert!(
                (u0 - um).norm() <= 1e-12 * um.norm(),
                "u({r},{s},0) = {u0}, membrane {um}"
            );
            let u1 = interpolate_displacement(&pre, r, s, 1.0, &dofs);
            let um1 = interpolate_displacement(&pre, r, s, -1.0, &dofs);
            assert!(
                (u1 - u0 - 0.5 * ub).norm() <= 1e-12 * ub.norm(),
                "u(t=1)-u(t=0) != half the director rotation"
            );
            assert!(
                (u1 - um1 - ub).norm() <= 1e-12 * ub.norm(),
                "u(t=1)-u(t=-1) != the full director rotation"
            );
        }

        // theta_z is carried in the DOF vector but only its components in the
        // plane of V_n matter: theta parallel to V_n produces no director
        // rotation, which is why the 2017 core is blind to the drill DOF.
        let mut drill = [0.0f64; 24];
        for i in 0..4 {
            let tz = 0.7 + i as f64;
            drill[6 * i + 3] = tz * pre.vn[i][0];
            drill[6 * i + 4] = tz * pre.vn[i][1];
            drill[6 * i + 5] = tz * pre.vn[i][2];
        }
        for &(r, s, t) in &[(0.3, -0.2, 1.0), (-0.6, 0.4, -1.0)] {
            let u = interpolate_displacement(&pre, r, s, t, &drill);
            assert!(
                u.norm() <= 1e-12,
                "theta parallel to V_n must give no director rotation, got {u}"
            );
        }
    }
}
