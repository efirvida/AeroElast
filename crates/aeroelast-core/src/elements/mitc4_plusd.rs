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

use nalgebra::{Matrix2, Matrix3, SMatrix, Vector3};

use crate::materials::ShellConstitutive;

/// 24×24 element matrix (4 nodes × 6 DOF/node).
pub type Mat24 = SMatrix<f64, 24, 24>;

/// The 2×2 Gauss rule on the element surface, `xi, eta = ±1/sqrt(3)` with unit
/// weights. Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410: "2 × 2 × 2 Gauss
/// integration over the element domain"; the surface part of that rule is the
/// 2×2 rule used here, and Ko, Bathe & Zhang (2025), C&S 308:107622, §3.1,
/// p. 13 requires "only the 2 × 2 Gauss integration over the element surfaces"
/// for the 2025 six-DOF contribution.
const N_GAUSS: usize = 4;
const GAUSS_XI: [f64; N_GAUSS] = [
    -0.577_350_269_189_625_8,
    0.577_350_269_189_625_8,
    0.577_350_269_189_625_8,
    -0.577_350_269_189_625_8,
];
const GAUSS_ETA: [f64; N_GAUSS] = [
    -0.577_350_269_189_625_8,
    -0.577_350_269_189_625_8,
    0.577_350_269_189_625_8,
    0.577_350_269_189_625_8,
];
const GAUSS_W: [f64; N_GAUSS] = [1.0, 1.0, 1.0, 1.0];

/// The bilinear node coordinates `(xi_i, eta_i)` of Ko, Bathe & Zhang (2025),
/// C&S 308:107622, Eq. (2), p. 2, in the repository's node order (node 0 at
/// `(-1,-1)`, counter-clockwise):
///
/// ```text
/// [xi_1 xi_2 xi_3 xi_4]  = [ 1 -1 -1  1 ]   (paper numbering)
/// [eta_1 eta_2 eta_3 eta_4] = [ 1  1 -1 -1 ]
/// ```
const NODE_XI: [f64; 4] = [-1.0, 1.0, 1.0, -1.0];
const NODE_ETA: [f64; 4] = [-1.0, -1.0, 1.0, 1.0];

/// Read node `i` of a four-node element as a 3D vector.
#[inline(always)]
fn node_vec(coords_3d: &[[f64; 3]; 4], i: usize) -> Vector3<f64> {
    Vector3::new(coords_3d[i][0], coords_3d[i][1], coords_3d[i][2])
}

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
    /// The scalar shear-correction factor the material model applied when it
    /// built `constitutive.cs` (ADR-1). `1.0` when it applied none — e.g. a
    /// multi-ply laminate's energy-equivalent section stiffness.
    pub applied_shear_correction: f64,
    /// The **uncorrected** transverse-shear stiffness `G·h` the shear block
    /// consumes (ADR-1). Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410: the
    /// element formulation "does not include any numerical factor". Equal to
    /// `constitutive.cs / applied_shear_correction` when that factor is positive
    /// (see `ShellConstitutive::transverse_shear_uncorrected`).
    pub cs_uncorrected: Matrix2<f64>,
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
    /// The five covariant membrane tying rows of Ko, Lee & Bathe (2017),
    /// C&S 182:404-418, Eqs. (15)-(17), sampled at the Fig. 4 tying points
    /// A(0,+1), B(0,-1), C(+1,0), D(-1,0), E(0,0). Translational slots only;
    /// each row is `e_rr`, `e_ss` or `e_rs` (not doubled).
    pub b_rr_a: [f64; 24],
    pub b_rr_b: [f64; 24],
    pub b_ss_c: [f64; 24],
    pub b_ss_d: [f64; 24],
    pub b_rs_e: [f64; 24],
    /// The four covariant transverse-shear tying operators of Dvorkin & Bathe
    /// (1984), Engineering Computations 1:77-88, Eq. (3), as reproduced in
    /// Ko, Lee & Bathe (2017), C&S 182:404-418, p. 405, ordered
    /// `[A(top,s=+1), B(bottom,s=-1), C(right,r=+1), D(left,r=-1)]`. Row 0 is
    /// the `e_rt` component, row 1 the `e_st` component.
    pub b_shear_tie: [SMatrix<f64, 2, 24>; 4],
    /// The four drill-membrane edge terms of Ko, Bathe & Zhang (2025),
    /// C&S 308:107622, Eqs. (13c)/(18), in the paper's edge order
    /// `l = 5..8` = [right, top, left, bottom].
    pub drill_edges: [DrillEdgeTerm; 4],
}

impl Mitc4PlusDPrecomputed {
    /// Build the precomputed element data from the 12 nodal coordinates
    /// `[x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4]`, the mid-surface constitutive
    /// and the shell thickness.
    ///
    /// The nodal directors are built locally (ADR-4 option B); no director or
    /// connectivity input is taken. `applied_shear_correction` is the scalar
    /// factor the material model applied to `constitutive.cs` (ADR-1); the
    /// constructor stores the uncorrected `G·h` the shear block consumes.
    pub fn new(
        node_coords: &[f64; 12],
        constitutive: ShellConstitutive,
        thickness: f64,
        applied_shear_correction: f64,
    ) -> Self {
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

        // The five covariant membrane tying rows of Ko, Lee & Bathe (2017),
        // C&S 182:404-418, Eqs. (15)-(17), at Fig. 4's A, B, C, D, E.
        let b_rr_a = covariant_membrane_b_row(&x_r, &x_s, &x_d, &e1, &e2, &e3, 0.0, 1.0, 0);
        let b_rr_b = covariant_membrane_b_row(&x_r, &x_s, &x_d, &e1, &e2, &e3, 0.0, -1.0, 0);
        let b_ss_c = covariant_membrane_b_row(&x_r, &x_s, &x_d, &e1, &e2, &e3, 1.0, 0.0, 1);
        let b_ss_d = covariant_membrane_b_row(&x_r, &x_s, &x_d, &e1, &e2, &e3, -1.0, 0.0, 1);
        let b_rs_e = covariant_membrane_b_row(&x_r, &x_s, &x_d, &e1, &e2, &e3, 0.0, 0.0, 2);
        // The four transverse-shear tying operators (DB84 Eq. 3): A(top),
        // B(bottom), C(right), D(left).
        let b_shear_tie = [
            compute_shear_tie(&coords_3d, &vn, &a_i, &e1, &e2, &e3, 0.0, 1.0),
            compute_shear_tie(&coords_3d, &vn, &a_i, &e1, &e2, &e3, 0.0, -1.0),
            compute_shear_tie(&coords_3d, &vn, &a_i, &e1, &e2, &e3, 1.0, 0.0),
            compute_shear_tie(&coords_3d, &vn, &a_i, &e1, &e2, &e3, -1.0, 0.0),
        ];
        // The four drill-membrane edge terms of Ko, Bathe & Zhang (2025),
        // C&S 308:107622, Eq. (13c)/(18).
        let drill_edges = compute_drill_edges(&coords_3d, &v_d);

        // ADR-1: remove the scalar the material model applied, so the shear
        // block is the uncorrected `G·h` of Ko, Lee & Bathe (2017),
        // C&S 182:404-418, p. 410.
        let cs_uncorrected = constitutive.transverse_shear_uncorrected(applied_shear_correction);

        Mitc4PlusDPrecomputed {
            local_coords,
            t3,
            e1,
            e2,
            e3,
            initial_coords_3d: coords_3d,
            constitutive,
            applied_shear_correction,
            cs_uncorrected,
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
            b_rr_a,
            b_rr_b,
            b_ss_c,
            b_ss_d,
            b_rs_e,
            b_shear_tie,
            drill_edges,
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
        let vn_local = local_components(pre, &pre.vn[i]);
        u += (t * 0.5) * pre.a_i[i] * h[i] * theta.cross(&vn_local);
    }
    u
}

// ============================================================================
// WU3 — the B-operators
// ============================================================================
//
// Citations are self-contained (author-year plus journal, volume and pages):
//
//   Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
//   Structures 182:404-418.
//   Dvorkin & Bathe (1984), "A continuum mechanics based four-node shell
//   element for general nonlinear analysis", Engineering Computations 1:77-88
//   (reproduced in Ko, Lee & Bathe (2017), C&S 182:404-418, p. 405).
//   Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements with
//   six degrees of freedom at each node - the MITC4/D and MITC4+/D elements",
//   Computers and Structures 308:107622.

/// The mid-surface Jacobian `j_loc[alpha][a] = g_a . e_alpha` (alpha = 1, 2;
/// a = r, s) at `(r, s)`.
///
/// With `g_r = x_r + s x_d` and `g_s = x_s + r x_d` (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eq. (9), p. 406), this is the matrix inverted by Eq. (11)
/// and used by [`covariant_to_local_mapping`].
fn j_loc_at(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> Matrix2<f64> {
    let (g_r, g_s, _) = compute_j3d_enriched(&pre.initial_coords_3d, &pre.vn, &pre.a_i, r, s, 0.0);
    Matrix2::new(
        g_r.dot(&pre.e1),
        g_s.dot(&pre.e1),
        g_r.dot(&pre.e2),
        g_s.dot(&pre.e2),
    )
}

/// Components of a 3D vector in the element's local orthonormal frame.
#[inline(always)]
fn local_components(pre: &Mitc4PlusDPrecomputed, v: &Vector3<f64>) -> Vector3<f64> {
    Vector3::new(v.dot(&pre.e1), v.dot(&pre.e2), v.dot(&pre.e3))
}

/// One displacement-based covariant membrane strain row (1x24) of
/// Ko, Lee & Bathe (2017), C&S 182:404-418, Eqs. (15)-(17), p. 408, sampled at
/// `(r, s)`. `component`: 0 = `e_rr`, 1 = `e_ss`, 2 = `e_rs` (not doubled).
///
/// ```text
/// e_rr^m(r,s) = (x_r + s x_d) . (u_r + s u_d)
/// e_ss^m(r,s) = (x_s + r x_d) . (u_s + r u_d)
/// e_rs^m(r,s) = 1/2[(x_r + s x_d).(u_s + r u_d) + (x_s + r x_d).(u_r + s u_d)]
/// ```
///
/// with `u_r = 1/4 sum xi_i u_i`, `u_s = 1/4 sum eta_i u_i`,
/// `u_d = 1/4 sum xi_i eta_i u_i` (Eq. 9). The rows carry translational slots
/// only. The five rows are sampled at the Fig. 4 tying points A(0,+1),
/// B(0,-1), C(+1,0), D(-1,0), E(0,0).
fn covariant_membrane_b_row(
    x_r: &Vector3<f64>,
    x_s: &Vector3<f64>,
    x_d: &Vector3<f64>,
    e1: &Vector3<f64>,
    e2: &Vector3<f64>,
    e3: &Vector3<f64>,
    r: f64,
    s: f64,
    component: usize,
) -> [f64; 24] {
    let g_r = x_r + s * x_d;
    let g_s = x_s + r * x_d;
    let gr = [g_r.dot(e1), g_r.dot(e2), g_r.dot(e3)];
    let gs = [g_s.dot(e1), g_s.dot(e2), g_s.dot(e3)];
    let mut b = [0.0f64; 24];
    for i in 0..4 {
        let a_r = 0.25 * NODE_XI[i];
        let a_s = 0.25 * NODE_ETA[i];
        let a_d = 0.25 * NODE_XI[i] * NODE_ETA[i];
        for k in 0..3 {
            b[6 * i + k] = match component {
                0 => (a_r + s * a_d) * gr[k],
                1 => (a_s + r * a_d) * gs[k],
                _ => 0.5 * ((a_s + r * a_d) * gr[k] + (a_r + s * a_d) * gs[k]),
            };
        }
    }
    b
}

/// The assumed covariant membrane strain rows `[e~_rr, e~_ss, e~_rs]` of
/// Ko, Lee & Bathe (2017), C&S 182:404-418, Eqs. (27a-c), p. 410.
///
/// Note F1 of `docs/formulations/mitc4plus-2017-extract.md`: Eqs. (27a-c) is
/// the closed form of Eqs. (21)+(26) and its Eq. (27c) coefficient
/// `(1 + a_E r s)` carries the leading `e_rs^m|bil` term the printed Eq. (21)
/// omits (note F2). The engineering shear `2 e_rs` is applied by
/// [`b_membrane_2017`] at the covariant-to-local mapping, not here.
fn b_membrane_covariant_2017(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 3, 24> {
    let a = pre.a_coeffs;
    let mut b = SMatrix::<f64, 3, 24>::zeros();
    for j in 0..24 {
        b[(0, j)] = 0.5 * (1.0 - 2.0 * a[0] + s + 2.0 * a[0] * s * s) * pre.b_rr_a[j]
            + 0.5 * (1.0 - 2.0 * a[1] - s + 2.0 * a[1] * s * s) * pre.b_rr_b[j]
            + a[2] * (-1.0 + s * s) * pre.b_ss_c[j]
            + a[3] * (-1.0 + s * s) * pre.b_ss_d[j]
            + a[4] * (-1.0 + s * s) * pre.b_rs_e[j];
        b[(1, j)] = a[0] * (-1.0 + r * r) * pre.b_rr_a[j]
            + a[1] * (-1.0 + r * r) * pre.b_rr_b[j]
            + 0.5 * (1.0 - 2.0 * a[2] + r + 2.0 * a[2] * r * r) * pre.b_ss_c[j]
            + 0.5 * (1.0 - 2.0 * a[3] - r + 2.0 * a[3] * r * r) * pre.b_ss_d[j]
            + a[4] * (-1.0 + r * r) * pre.b_rs_e[j];
        b[(2, j)] = 0.25 * (r + 4.0 * a[0] * r * s) * pre.b_rr_a[j]
            + 0.25 * (-r + 4.0 * a[1] * r * s) * pre.b_rr_b[j]
            + 0.25 * (s + 4.0 * a[2] * r * s) * pre.b_ss_c[j]
            + 0.25 * (-s + 4.0 * a[3] * r * s) * pre.b_ss_d[j]
            + (1.0 + a[4] * r * s) * pre.b_rs_e[j];
    }
    b
}

/// The MITC4+ assumed membrane strain of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (27a-c) mapped point-wise to the local orthonormal
/// frame. Returns the local `[e_11, e_22, 2 e_12]` rows.
fn b_membrane_2017(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 3, 24> {
    let mut cov = b_membrane_covariant_2017(pre, r, s);
    for j in 0..24 {
        cov[(2, j)] *= 2.0;
    }
    let t = covariant_to_local_mapping(&j_loc_at(pre, r, s));
    t * cov
}

/// The displacement-based bending strain rows of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (7c)/(7d), p. 406, including the `dx_b . du_m`
/// contribution of Eq. (8a), p. 406.
///
/// Returns the covariant `(e^b1, e^b2)` triples, each `[e_rr, e_ss, 2 e_rs]`,
/// so the caller maps them with [`covariant_to_local_mapping`].
fn b_bending_covariant_2017(
    pre: &Mitc4PlusDPrecomputed,
    r: f64,
    s: f64,
) -> (SMatrix<f64, 3, 24>, SMatrix<f64, 3, 24>) {
    let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
    // x_m tangents (Eq. 9), in the local frame.
    let xm_r = local_components(pre, &(pre.x_r + s * pre.x_d));
    let xm_s = local_components(pre, &(pre.x_s + r * pre.x_d));
    // x_b tangents (Eq. 8a); V_n^i are stored in the local frame.
    let mut xb_r = Vector3::zeros();
    let mut xb_s = Vector3::zeros();
    for i in 0..4 {
        let vn_local = local_components(pre, &pre.vn[i]);
        xb_r += 0.5 * pre.a_i[i] * dn_dr[i] * vn_local;
        xb_s += 0.5 * pre.a_i[i] * dn_ds[i] * vn_local;
    }

    let mut b1 = SMatrix::<f64, 3, 24>::zeros();
    let mut b2 = SMatrix::<f64, 3, 24>::zeros();
    for i in 0..4 {
        // Coefficients of u_i in du_m/dr and du_m/ds.
        let um_r = 0.25 * NODE_XI[i] * (1.0 + s * NODE_ETA[i]);
        let um_s = 0.25 * NODE_ETA[i] * (1.0 + r * NODE_XI[i]);
        // Coefficients of (theta_i x V_n^i) in du_b/dr and du_b/ds.
        let ub_r = 0.5 * pre.a_i[i] * dn_dr[i];
        let ub_s = 0.5 * pre.a_i[i] * dn_ds[i];
        // The rotation DOFs are local-frame components, so the director must be
        // projected to the local frame before `theta x V_n` is formed (ADR-2:
        // the stored `pre.vn` is a physical, global-frame vector).
        let vn = local_components(pre, &pre.vn[i]);
        // (theta x V)[k] = sum_bb L[k][bb] theta[bb].
        let cross = [
            [0.0, vn[2], -vn[1]],
            [-vn[2], 0.0, vn[0]],
            [vn[1], -vn[0], 0.0],
        ];
        for k in 0..3 {
            // Eq. (7c), the dx_b . du_m part of Eq. (8a).
            b1[(0, 6 * i + k)] += xb_r[k] * um_r;
            b1[(1, 6 * i + k)] += xb_s[k] * um_s;
            b1[(2, 6 * i + k)] += 0.5 * (xb_r[k] * um_s + xb_s[k] * um_r);
            for bb in 0..3 {
                let ubk = ub_r * cross[k][bb];
                let ubs = ub_s * cross[k][bb];
                // Eq. (7c), the dx_m . du_b part.
                b1[(0, 6 * i + 3 + bb)] += xm_r[k] * ubk;
                b1[(1, 6 * i + 3 + bb)] += xm_s[k] * ubs;
                b1[(2, 6 * i + 3 + bb)] += 0.5 * (xm_r[k] * ubs + xm_s[k] * ubk);
                // Eq. (7d).
                b2[(0, 6 * i + 3 + bb)] += xb_r[k] * ubk;
                b2[(1, 6 * i + 3 + bb)] += xb_s[k] * ubs;
                b2[(2, 6 * i + 3 + bb)] += 0.5 * (xb_r[k] * ubs + xb_s[k] * ubk);
            }
        }
    }
    // Engineering shear: double the third covariant row.
    for j in 0..24 {
        b1[(2, j)] *= 2.0;
        b2[(2, j)] *= 2.0;
    }
    (b1, b2)
}

/// The MITC4+ bending strain operators of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eqs. (7c)/(7d), mapped to the local frame. Returns
/// `(B_b1, B_b2)`, each 3x24. The matrix is 24x24 in the assembly: there is no
/// condensed internal DOF.
fn b_bending_2017(
    pre: &Mitc4PlusDPrecomputed,
    r: f64,
    s: f64,
) -> (SMatrix<f64, 3, 24>, SMatrix<f64, 3, 24>) {
    let (c1, c2) = b_bending_covariant_2017(pre, r, s);
    let t = covariant_to_local_mapping(&j_loc_at(pre, r, s));
    (t * c1, t * c2)
}

// ============================================================================
// MITC4 assumed transverse shear (Dvorkin & Bathe 1984, Eq. 3, reproduced in
// Ko, Lee & Bathe (2017), C&S 182:404-418, p. 405)
// ============================================================================

/// Build the covariant transverse-shear tying operator (2x24) at one tying
/// point. Row 0 is `e_rt`, row 1 is `e_st` (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eq. (4), p. 405):
///
/// ```text
/// e_rt = 1/2 (g_r . u_t + g_t . u_r)     e_st = 1/2 (g_s . u_t + g_t . u_s)
/// ```
///
/// with `u_t = u_b = 1/2 sum a_i h_i (theta_i x V_n^i)` (Eq. 3, Eq. 8b) and
/// `u_r = du_m/dr`, `u_s = du_m/ds` at `t = 0` (Eq. 9). The metric
/// normalization to the local frame is applied point-wise afterwards by
/// [`shear_covariant_to_local`]; no factor is introduced here.
#[allow(clippy::too_many_arguments)]
fn compute_shear_tie(
    coords_3d: &[[f64; 3]; 4],
    vn: &[Vector3<f64>; 4],
    a_i: &[f64; 4],
    e1: &Vector3<f64>,
    e2: &Vector3<f64>,
    e3: &Vector3<f64>,
    r: f64,
    s: f64,
) -> SMatrix<f64, 2, 24> {
    let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
    let h = shape_functions(r, s);
    let (g_r_g, g_s_g, g_t_g) = compute_j3d_enriched(coords_3d, vn, a_i, r, s, 0.0);
    let g_r = Vector3::new(g_r_g.dot(e1), g_r_g.dot(e2), g_r_g.dot(e3));
    let g_s = Vector3::new(g_s_g.dot(e1), g_s_g.dot(e2), g_s_g.dot(e3));
    let g_t = Vector3::new(g_t_g.dot(e1), g_t_g.dot(e2), g_t_g.dot(e3));
    let mut b = SMatrix::<f64, 2, 24>::zeros();
    for i in 0..4 {
        // g_r . u_t = 1/2 a_i h_i theta_i . (V_n^i x g_r).
        // Local-frame director, so `V_n x g_a` pairs with the local DOFs.
        let vn_i = Vector3::new(vn[i].dot(e1), vn[i].dot(e2), vn[i].dot(e3));
        let vxg_r = vn_i.cross(&g_r);
        let vxg_s = vn_i.cross(&g_s);
        for k in 0..3 {
            b[(0, 6 * i + k)] = 0.5 * g_t[k] * dn_dr[i];
            b[(1, 6 * i + k)] = 0.5 * g_t[k] * dn_ds[i];
            b[(0, 6 * i + 3 + k)] = 0.25 * a_i[i] * h[i] * vxg_r[k];
            b[(1, 6 * i + 3 + k)] = 0.25 * a_i[i] * h[i] * vxg_s[k];
        }
    }
    b
}

/// The MITC4 assumed transverse shear of Dvorkin & Bathe (1984),
/// Engineering Computations 1:77-88, Eq. (3), as reproduced in Ko, Lee &
/// Bathe (2017), C&S 182:404-418, p. 405:
///
/// ```text
/// e~_rt = 1/2(1+s) e_rt^A + 1/2(1-s) e_rt^B      tying A (top, s=+1), B (bottom, s=-1)
/// e~_st = 1/2(1+r) e_st^C + 1/2(1-r) e_st^D      tying C (right, r=+1), D (left, r=-1)
/// ```
///
/// The covariant field is mapped to the local frame by
/// [`shear_covariant_to_local`] (the 3D dual-basis form of the design, since
/// Ko, Lee & Bathe (2017), C&S 182:404-418 prints the covariant definition
/// (Eq. 4) but defers the tying-point construction to Dvorkin & Bathe (1984)).
/// Returns the local `[gamma_13, gamma_23]` rows (2x24).
fn b_shear_mitc4(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 2, 24> {
    let mut cov = SMatrix::<f64, 2, 24>::zeros();
    for j in 0..24 {
        cov[(0, j)] = 0.5 * (1.0 + s) * pre.b_shear_tie[0][(0, j)]
            + 0.5 * (1.0 - s) * pre.b_shear_tie[1][(0, j)];
        cov[(1, j)] = 0.5 * (1.0 + r) * pre.b_shear_tie[2][(1, j)]
            + 0.5 * (1.0 - r) * pre.b_shear_tie[3][(1, j)];
    }
    let (g_r, g_s, g_t) =
        compute_j3d_enriched(&pre.initial_coords_3d, &pre.vn, &pre.a_i, r, s, 0.0);
    let t = shear_covariant_to_local(&g_r, &g_s, &g_t, &pre.e1, &pre.e2, &pre.e3);
    t * cov
}

// ============================================================================
// Drill-membrane strain (Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. 18)
// ============================================================================

/// The four edge terms `(start, end, c_r^l, c_s^l)` of the Eq. (18)
/// drill-membrane operator, in the paper's edge order `l = 5..8` =
/// [right, top, left, bottom].
#[derive(Clone, Copy, Debug)]
pub struct DrillEdgeTerm {
    /// Code node `i` of the paper's edge node pair `(i, i+1)` (subtracted).
    pub start: usize,
    /// Code node `i+1` (added).
    pub end: usize,
    /// `c_r^l = x_m^l . (-x_r^l x V^D)`, Eq. (18)/(19c).
    pub c_r: f64,
    /// `c_s^l = x_m^l . ( x_s^l x V^D)`, Eq. (18)/(19c).
    pub c_s: f64,
}

/// The code node pair `(i, i+1)` of each drill edge, in the paper's order
/// `l = 5..8` = [right, top, left, bottom]. Derived from Eq. (2)'s paper node
/// numbering `1=(+1,+1), 2=(-1,+1), 3=(-1,-1), 4=(+1,-1)` and Eq. (13c)'s
/// edge vectors `x_m^5 = 1/8(x_4 - x_1)`, `x_m^6 = 1/8(x_1 - x_2)`,
/// `x_m^7 = 1/8(x_2 - x_3)`, `x_m^8 = 1/8(x_3 - x_4)`.
const DRILL_EDGE_NODES: [(usize, usize); 4] = [(1, 2), (2, 3), (3, 0), (0, 1)];

/// The mid-point `(r, s)` of each drill edge, in the same order.
const DRILL_EDGE_MID: [(f64, f64); 4] = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];

/// Compute the four edge terms of Eq. (13c)/(18): `x_m^l = 1/8 (x_i - x_{i+1})`
/// and `c_r^l`, `c_s^l` from the mid-surface tangents at the edge mid-point.
fn compute_drill_edges(coords_3d: &[[f64; 3]; 4], v_d: &Vector3<f64>) -> [DrillEdgeTerm; 4] {
    let mut out = [DrillEdgeTerm {
        start: 0,
        end: 0,
        c_r: 0.0,
        c_s: 0.0,
    }; 4];
    for e in 0..4 {
        let (start, end) = DRILL_EDGE_NODES[e];
        let (r_m, s_m) = DRILL_EDGE_MID[e];
        // Eq. (13c): x_m^l = 1/8 (x_i - x_{i+1}), so ||x_m^l|| = L_l / 8.
        let x_m = 0.125 * (node_vec(coords_3d, start) - node_vec(coords_3d, end));
        // Eq. (13b): the mid-surface tangents at the edge mid-point.
        let (dn_dr, dn_ds) = shape_function_derivatives(r_m, s_m);
        let mut x_r = Vector3::zeros();
        let mut x_s = Vector3::zeros();
        for j in 0..4 {
            x_r += dn_dr[j] * node_vec(coords_3d, j);
            x_s += dn_ds[j] * node_vec(coords_3d, j);
        }
        out[e] = DrillEdgeTerm {
            start,
            end,
            c_r: x_m.dot(&(-x_r.cross(v_d))),
            c_s: x_m.dot(&(x_s.cross(v_d))),
        };
    }
    out
}

/// The simplified ("curl") derivatives of the four drill mid-side functions,
/// Ko, Bathe & Zhang (2025), C&S 308:107622, Eqs. (10), (11a)/(11b), p. 5:
///
/// ```text
/// [h~_5,r h~_6,r h~_7,r h~_8,r] = [ 0, -r(1+s), 0, -r(1-s) ]      (11a)
/// [h~_5,s h~_6,s h~_7,s h~_8,s] = [ -s(1+r), 0, -s(1-r), 0 ]      (11b)
/// ```
///
/// The zeros are the paper's deliberate simplification (text below Eq. 11b),
/// which keeps the required integration order without losing the patch tests.
/// The return order is the paper's [right, top, left, bottom] = [5, 6, 7, 8].
#[inline(always)]
fn drill_midside_shape_derivatives(r: f64, s: f64) -> [(f64, f64); 4] {
    [
        (0.0, -s * (1.0 + r)),
        (-r * (1.0 + s), 0.0),
        (0.0, -s * (1.0 - r)),
        (-r * (1.0 - s), 0.0),
    ]
}

/// The ratio `j0 / j` of Eq. (17b), Ko, Bathe & Zhang (2025),
/// C&S 308:107622, p. 8, with `j = det[g_r g_s g_t]|_{(r,s,0)}` and
/// `j0 = j(0,0,0)`. Exact scalar triple product, not a surface-measure ratio.
fn drill_jacobian_ratio(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> f64 {
    let (g_r, g_s, g_t) =
        compute_j3d_enriched(&pre.initial_coords_3d, &pre.vn, &pre.a_i, r, s, 0.0);
    let j = g_r.cross(&g_s).dot(&g_t);
    if j.abs() > 1e-30 {
        pre.j0 / j
    } else {
        0.0
    }
}

/// The drill-membrane strain-displacement operator `B_md` (3x24) of
/// Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (18), p. 8:
///
/// ```text
/// e~_rr^md = (j0/j) h~_m,r^l (theta_{i+1}^D - theta_i^D) x_m^l . (-x_r^l x V^D)
/// e~_ss^md = -(j0/j) h~_m,s^l (theta_{i+1}^D - theta_i^D) x_m^l . ( x_s^l x V^D)
/// e~_rs^md = 1/2 (j0/j) [ h~_m,s^l x_m^l . (-x_r^l x V^D)
///                        - h~_m,r^l x_m^l . ( x_s^l x V^D) ] (theta_{i+1}^D - theta_i^D)
/// ```
///
/// with the drill rotation `theta_i^D = theta_i . V^D` of Eq. (6a)/(6f), one
/// element-normal `V^D` of Eq. (5), and the edge `l` joining nodes `i` and
/// `i+1` (Eq. 16a/16b, confirmed by Eq. 12d). The covariant operator is mapped
/// with the constant element-centre base vectors of Eq. (21) (ADR-3).
fn b_drill_membrane_2025(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 3, 24> {
    let ratio = drill_jacobian_ratio(pre, r, s);
    let dh = drill_midside_shape_derivatives(r, s);
    // `theta^D = theta . V^D` with `theta` a local-frame DOF triple, so V^D is
    // projected to the local frame (the `c_r`/`c_s` edge scalars use the global
    // V^D and are frame-invariant).
    let vd_local = local_components(pre, &pre.v_d);
    let mut b_cov = SMatrix::<f64, 3, 24>::zeros();
    for e in 0..4 {
        let edge = pre.drill_edges[e];
        let (h_r, h_s) = dh[e];
        let coeff_rr = ratio * h_r * edge.c_r;
        let coeff_ss = -ratio * h_s * edge.c_s;
        let coeff_rs = 0.5 * ratio * (h_s * edge.c_r - h_r * edge.c_s);
        for beta in 0..3 {
            let vd = vd_local[beta];
            b_cov[(0, 6 * edge.start + 3 + beta)] -= coeff_rr * vd;
            b_cov[(0, 6 * edge.end + 3 + beta)] += coeff_rr * vd;
            b_cov[(1, 6 * edge.start + 3 + beta)] -= coeff_ss * vd;
            b_cov[(1, 6 * edge.end + 3 + beta)] += coeff_ss * vd;
            b_cov[(2, 6 * edge.start + 3 + beta)] -= 2.0 * coeff_rs * vd;
            b_cov[(2, 6 * edge.end + 3 + beta)] += 2.0 * coeff_rs * vd;
        }
    }
    let t = covariant_to_local_mapping(&j_loc_at(pre, 0.0, 0.0));
    t * b_cov
}

// ============================================================================
// WU4 — the stiffness assembly
// ============================================================================
//
// Citations are self-contained (author-year plus journal, volume and pages):
//
//   Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
//   Structures 182:404-418 (Eqs. 1-27; the 2x2x2 Gauss rule, p. 410).
//   Dvorkin & Bathe (1984), "A continuum mechanics based four-node shell
//   element for general nonlinear analysis", Engineering Computations 1:77-88
//   (Eq. 3, reproduced in Ko, Lee & Bathe (2017), C&S 182:404-418, p. 405).
//   Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements with
//   six degrees of freedom at each node - the MITC4/D and MITC4+/D elements",
//   Computers and Structures 308:107622 (Eqs. 5, 18, 21, 26; the 2x2 surface
//   rule of SS3.1, p. 13).

/// The through-thickness resultant/moment block matrix `W` (9x9) of the
/// design's Section 2.3.
///
/// Writing the strain decomposition of Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eq. (7a), as the **resultant** triple `[eps_m, kappa, E2]`
/// with `kappa = (2/h) e^b1` and `E2 = (4/h^2) e^b2`, the through-thickness
/// energy `int e(z)^T C e(z) dz` with
/// `e(z) = e^m + (2z/h) e^b1 + (4z^2/h^2) e^b2` expands to the blocks
///
/// ```text
/// W_00 = int C dz      = cm            (membrane)
/// W_01 = int z C dz    = cb_coupling   (membrane-bending coupling)
/// W_02 = int z^2 C dz  = cb            (membrane-E2 coupling)
/// W_11 = int z^2 C dz  = cb            (bending)
/// W_12 = int z^3 C dz  = 0             (odd moment)
/// W_22 = int z^4 C dz  = cm/9          (E2)
/// ```
///
/// so the repository's ABD blocks act verbatim and no scaling factor appears in
/// the stiffness assembly. The `cm/9` of `W_22` is **the paper's own 2x2 rule in
/// the through-thickness coordinate `t`** (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, p. 410: "2 x 2 x 2 Gauss integration over the element
/// domain"): two-point Gauss in `t` gives `int t^4 dt -> 2/9` where exact
/// integration gives `2/5`, so
/// `int z^4 C dz -> (16/h^4)(h/2)(h/2)^4(2/9) C = h C/9 = cm/9`. It is the
/// paper's own under-integration, not a factor invented here.
///
/// RECORDED APPROXIMATION (design Section 9, risk 7). For a homogeneous section
/// `cm/9` is exact under the paper's rule. For a **multi-ply laminate** it uses
/// the section's thickness-average membrane stiffness (the same smearing
/// `cm_raw` already documents) in the `E2-E2` block only, instead of the true
/// `int z^4 C(z) dz`; the term it multiplies is second order in the warping, and
/// the surface-only discrimination of
/// `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only`
/// makes the approximation visible to a test rather than hidden.
pub fn resultant_moment_matrix(constitutive: &ShellConstitutive) -> SMatrix<f64, 9, 9> {
    let cm = &constitutive.cm;
    let cb = &constitutive.cb;
    let cbc = &constitutive.cb_coupling;
    let mut w = SMatrix::<f64, 9, 9>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            w[(i, j)] = cm[(i, j)];
            w[(i, 3 + j)] = cbc[(i, j)];
            w[(3 + i, j)] = cbc[(j, i)];
            w[(i, 6 + j)] = cb[(i, j)];
            w[(6 + i, j)] = cb[(j, i)];
            w[(3 + i, 3 + j)] = cb[(i, j)];
            w[(6 + i, 6 + j)] = cm[(i, j)] / 9.0;
        }
    }
    w
}

/// Mid-surface area measure `sqrt_g = ||g_r x g_s||` at `(r, s, t = 0)`, with
/// `g_r = x_r + s x_d` and `g_s = x_s + r x_d` (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eq. (9), p. 406). Every surface contribution integrates with
/// it: the 2x2x2 rule of paper A p. 410 and the 2x2 surface rule of paper B
/// §3.1 p. 13, including the drill term, whose `j0/j` factor of Eq. (17b) is
/// already inside [`b_drill_membrane_2025`].
fn surface_measure(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> f64 {
    let (g_r, g_s, _) = compute_j3d_enriched(&pre.initial_coords_3d, &pre.vn, &pre.a_i, r, s, 0.0);
    g_r.cross(&g_s).norm()
}

/// The membrane block `sum_g B_m^T cm B_m w sqrt_g` (local, 24x24). Exposed to
/// the identity lock, which compares it with the test-local reference block.
fn membrane_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
    let cm = &pre.constitutive.cm;
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let bm = b_membrane_2017(pre, r, s);
        k += (bm.transpose() * cm * bm) * (GAUSS_W[g] * sqrt_g);
    }
    k
}

/// The transverse-shear block
/// `sum_g B_gamma^T cs_uncorrected B_gamma w sqrt_g` (local, 24x24), with the
/// **uncorrected** `G h` of ADR-1 (Ko, Lee & Bathe (2017), C&S 182:404-418,
/// p. 410: the element "does not include any numerical factor") and the MITC4
/// assumed field of Dvorkin & Bathe (1984), Engineering Computations 1:77-88,
/// Eq. (3).
fn shear_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
    let cs = pre.cs_uncorrected;
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let bg = b_shear_mitc4(pre, r, s);
        k += (bg.transpose() * cs * bg) * (GAUSS_W[g] * sqrt_g);
    }
    k
}

/// The drill-membrane block `sum_g B_md^T cm B_md w sqrt_g` (local, 24x24) of
/// Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (26) with the operator of
/// Eq. (18) and the 2x2 surface rule of §3.1, p. 13. Symmetrised exactly so the
/// provenance assertions of
/// `test_identity_drill_stiffness_comes_only_from_eq26` can be exact; the
/// symmetrisation is a round-off guard, not a formulation factor.
fn drill_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
    let cm = &pre.constitutive.cm;
    let mut m = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let bmd = b_drill_membrane_2025(pre, r, s);
        m += (bmd.transpose() * cm * bmd) * (GAUSS_W[g] * sqrt_g);
    }
    let mut out = m;
    for i in 0..24 {
        for j in 0..24 {
            out[(i, j)] = 0.5 * (m[(i, j)] + m[(j, i)]);
        }
    }
    out
}

/// The full local stiffness with the B Eq. (26) drill contribution switched on
/// or off, so the provenance test can form the exact difference
/// `K(operator) - K(operator := 0)` the spec fixes.
fn compute_ke_local_with_drill(pre: &Mitc4PlusDPrecomputed, use_drill: bool) -> Mat24 {
    let w = resultant_moment_matrix(&pre.constitutive);
    let s1 = 2.0 / pre.thickness;
    let s2 = 4.0 / (pre.thickness * pre.thickness);
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let bm = b_membrane_2017(pre, r, s);
        let (bb1, bb2) = b_bending_2017(pre, r, s);
        let mut b = SMatrix::<f64, 9, 24>::zeros();
        for i in 0..3 {
            for j in 0..24 {
                b[(i, j)] = bm[(i, j)];
                b[(3 + i, j)] = s1 * bb1[(i, j)];
                b[(6 + i, j)] = s2 * bb2[(i, j)];
            }
        }
        k += (b.transpose() * w * b) * (GAUSS_W[g] * sqrt_g);
    }
    k += shear_ke_local(pre);
    if use_drill {
        k += drill_ke_local(pre);
    }
    // Round-off guard only (no formulation factor): the assembled sum is
    // exactly symmetric in exact arithmetic, and the exact-symmetry and
    // exact-zero assertions of the drill provenance test need it to be so in
    // floating point too. The base and the drill block are each symmetric, so
    // this changes no value beyond the last bit.
    let mut out = k;
    for i in 0..24 {
        for j in 0..24 {
            out[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
        }
    }
    out
}

/// Compute the element stiffness in LOCAL coordinates (24x24):
///
/// ```text
/// K = sum_g [B_m; (2/h) B_b1; (4/h^2) B_b2]^T W [..] w sqrt_g
///   + sum_g B_gamma^T cs_uncorrected B_gamma w sqrt_g
///   + sum_g B_md^T cm B_md w sqrt_g
/// ```
///
/// with the moment matrix `W` of [`resultant_moment_matrix`], the 2x2 surface
/// rule (Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410; Ko, Bathe & Zhang
/// (2025), C&S 308:107622, §3.1, p. 13) and **no numerical factor**: the drill
/// term is the penalty-free B Eq. (26) strain and the shear block is the
/// uncorrected `G h` of ADR-1.
pub fn compute_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
    compute_ke_local_with_drill(pre, true)
}

/// 24x24 global-to-local transformation matrix `T = diag(t3, ..., t3)`, one
/// 3x3 block per translational and per rotational nodal triple.
fn build_t24(pre: &Mitc4PlusDPrecomputed) -> SMatrix<f64, 24, 24> {
    let mut t24 = SMatrix::<f64, 24, 24>::zeros();
    for i in 0..8 {
        let r = 3 * i;
        for a in 0..3 {
            for b in 0..3 {
                t24[(r + a, r + b)] = pre.t3[(a, b)];
            }
        }
    }
    t24
}

/// Transform a 24x24 local matrix to global coordinates: `T^T M T`.
fn transform_to_global(pre: &Mitc4PlusDPrecomputed, m_local: &Mat24) -> Mat24 {
    let t24 = build_t24(pre);
    let mut t24t = t24;
    for i in 0..24 {
        for j in 0..24 {
            t24t[(i, j)] = t24[(j, i)];
        }
    }
    let k = t24t * m_local * t24;
    let mut out = k;
    for i in 0..24 {
        for j in 0..24 {
            out[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
        }
    }
    out
}

/// Compute the element stiffness in GLOBAL coordinates (24x24).
pub fn compute_ke_global(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
    transform_to_global(pre, &compute_ke_local(pre))
}

#[cfg(test)]
mod tests {
    use super::{
        b_bending_2017, b_bending_covariant_2017, b_drill_membrane_2025, b_membrane_2017,
        b_membrane_covariant_2017, b_shear_mitc4, compute_ke_global, compute_ke_local,
        compute_ke_local_with_drill, compute_membrane_coefficients_2017,
        covariant_to_local_mapping, drill_ke_local, interpolate_displacement, interpolate_position,
        j_loc_at, membrane_ke_local, node_vec, resultant_moment_matrix, shape_function_derivatives,
        shape_functions, shear_ke_local, surface_measure, Mat24, Mitc4PlusDPrecomputed,
        DRILL_EDGE_MID, GAUSS_ETA, GAUSS_W, GAUSS_XI, NODE_ETA, NODE_XI,
    };
    use crate::materials::laminate::{Laminate, Ply};
    use crate::materials::orthotropic::OrthotropicMaterial;
    use crate::materials::{isotropic::IsotropicMaterial, Material, ShellConstitutive};
    use nalgebra::{DMatrix, Matrix2, SMatrix, Vector3};

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
        Mitc4PlusDPrecomputed::new(&coords12(c), shell_iso(), 1.0, 5.0 / 6.0)
    }

    /// A flat unit square in the repository's `(r,s)` node order.
    const FLAT_SQUARE: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];

    /// A strongly warped quad: the node `z` offsets are of the order of the
    /// element's in-plane size, so the `t^2` (`E2`) term of the through-thickness
    /// rule is not a borderline effect. Used only for the 2x2x2 vs surface-only
    /// discrimination (design Section 2.3); it is not a constant-stress fixture.
    const STRONGLY_WARPED: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.9],
        [2.0, 1.0, -0.9],
        [0.0, 1.0, 0.0],
    ];

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
                let vn_local = local_components(&pre, &pre.vn[i]);
                ub += pre.a_i[i] * h[i] * theta.cross(&vn_local);
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
            let vn_local = local_components(&pre, &pre.vn[i]);
            drill[6 * i + 3] = tz * vn_local[0];
            drill[6 * i + 4] = tz * vn_local[1];
            drill[6 * i + 5] = tz * vn_local[2];
        }
        for &(r, s, t) in &[(0.3, -0.2, 1.0), (-0.6, 0.4, -1.0)] {
            let u = interpolate_displacement(&pre, r, s, t, &drill);
            assert!(
                u.norm() <= 1e-12,
                "theta parallel to V_n must give no director rotation, got {u}"
            );
        }
    }

    // ========================================================================
    // WU3 — the B-operators (tasks 4.1-4.5)
    // ========================================================================

    /// The 2x2 Gauss coordinate `±1/sqrt(3)`.
    const GP: f64 = 0.577_350_269_189_625_8;

    /// A generic 24-vector with distinct, non-trivial entries.
    fn generic_dofs() -> [f64; 24] {
        let mut d = [0.0f64; 24];
        let mut v = 0.037_f64;
        for (k, dk) in d.iter_mut().enumerate() {
            v = v * 1.7 + 0.019 * (k as f64 + 1.0);
            *dk = v.sin();
        }
        d
    }

    fn row_dot_3(m: &SMatrix<f64, 3, 24>, row: usize, dofs: &[f64; 24]) -> f64 {
        (0..24).map(|j| m[(row, j)] * dofs[j]).sum()
    }

    fn row_dot_2(m: &SMatrix<f64, 2, 24>, row: usize, dofs: &[f64; 24]) -> f64 {
        (0..24).map(|j| m[(row, j)] * dofs[j]).sum()
    }

    fn scale_3(m: &SMatrix<f64, 3, 24>) -> f64 {
        let mut s = 0.0f64;
        for i in 0..3 {
            for j in 0..24 {
                s = s.max(m[(i, j)].abs());
            }
        }
        s
    }

    /// Components of a 3D vector in the element's local orthonormal frame.
    fn local_components(pre: &Mitc4PlusDPrecomputed, v: &Vector3<f64>) -> Vector3<f64> {
        Vector3::new(v.dot(&pre.e1), v.dot(&pre.e2), v.dot(&pre.e3))
    }

    #[test]
    fn test_t1a_membrane_eq22_flat_tying_condition() {
        // Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (22): on a flat element
        // (x_d . n = 0) the assumed bilinear shear coefficient e~_rs^m|bil
        // equals its displacement-based counterpart e_rs^m|bil = x_d . u_d
        // (Eq. 16).
        let dofs = generic_dofs();

        // The assumed bilinear coefficient, extracted exactly from the
        // production covariant field
        //   e~_rs(r,s) = e_rs(E) + 1/2 e_rr|lin r + 1/2 e_ss|lin s + bil r s
        // (Eq. 27c) by the four-corner second difference.
        let assumed_bil = |pre: &Mitc4PlusDPrecomputed| {
            let f = |r: f64, s: f64| row_dot_3(&b_membrane_covariant_2017(pre, r, s), 2, &dofs);
            0.25 * (f(1.0, 1.0) - f(1.0, -1.0) - f(-1.0, 1.0) + f(-1.0, -1.0))
        };

        // The independent Eq. (16)/(17) strain parts, recomputed from the DOFs
        // and the stored geometry: (e_rr|con, e_ss|con, e_rs|con, e_rr|lin,
        // e_ss|lin, e_rs|bil = x_d . u_d).
        let parts = |pre: &Mitc4PlusDPrecomputed| {
            let mut u_r = Vector3::zeros();
            let mut u_s = Vector3::zeros();
            let mut u_d = Vector3::zeros();
            for i in 0..4 {
                let ui = Vector3::new(dofs[6 * i], dofs[6 * i + 1], dofs[6 * i + 2]);
                u_r += 0.25 * NODE_XI[i] * ui;
                u_s += 0.25 * NODE_ETA[i] * ui;
                u_d += 0.25 * NODE_XI[i] * NODE_ETA[i] * ui;
            }
            let x_r = local_components(pre, &pre.x_r);
            let x_s = local_components(pre, &pre.x_s);
            let x_d = local_components(pre, &pre.x_d);
            (
                x_r.dot(&u_r),
                x_s.dot(&u_s),
                0.5 * (x_r.dot(&u_s) + x_s.dot(&u_r)),
                x_r.dot(&u_d) + x_d.dot(&u_r),
                x_s.dot(&u_d) + x_d.dot(&u_s),
                x_d.dot(&u_d),
            )
        };
        let true_bil = |pre: &Mitc4PlusDPrecomputed| parts(pre).5;

        // Flat distorted element: x_d is in-plane (x_d . n = 0), so Eq. (22)
        // must hold exactly.
        let flat = pre_from(&FLAT_DISTORTED);
        assert!(
            flat.x_d.dot(&flat.n_vec).abs() <= 1e-14 * flat.x_d.norm().max(1.0),
            "flat fixture is not flat: x_d . n = {}",
            flat.x_d.dot(&flat.n_vec)
        );
        let a_flat = assumed_bil(&flat);
        let t_flat = true_bil(&flat);
        assert!(
            t_flat.abs() > 1e-6,
            "flat fixture gives a zero bilinear coefficient; the check would be vacuous"
        );
        assert!(
            (a_flat - t_flat).abs() <= 1e-14,
            "flat element: Eq. (22) violated: e~_rs^m|bil = {a_flat}, e_rs^m|bil = {t_flat}"
        );

        // The engineering-shear factor `2 e_rs` is applied at the
        // covariant-to-local mapping, not inside Eq. (27c): the mapped operator
        // equals the covariant field with its third row doubled, followed by
        // the point-wise mapping.
        let (rr, ss) = (0.37, -0.62);
        let mapped = b_membrane_2017(&flat, rr, ss);
        let mut cov2 = b_membrane_covariant_2017(&flat, rr, ss);
        for j in 0..24 {
            cov2[(2, j)] *= 2.0;
        }
        let expect = covariant_to_local_mapping(&j_loc_at(&flat, rr, ss)) * cov2;
        for i in 0..3 {
            for j in 0..24 {
                assert!(
                    (mapped[(i, j)] - expect[(i, j)]).abs() <= 1e-14,
                    "engineering-shear placement [{i}][{j}]: {} vs {}",
                    mapped[(i, j)],
                    expect[(i, j)]
                );
            }
        }

        // Warped element: x_d . n != 0, so Eq. (22) must fail measurably.
        let warped = pre_from(&DOUBLY_WARPED);
        assert!(
            warped.x_d.dot(&warped.n_vec).abs() > 1e-3,
            "warped fixture is flat"
        );
        let a_w = assumed_bil(&warped);
        let t_w = true_bil(&warped);
        let scale = a_w.abs().max(t_w.abs()).max(1e-30);
        assert!(
            (a_w - t_w).abs() > 1e-6 * scale,
            "warped element: Eq. (22) not violated, no separation ({a_w} vs {t_w})"
        );

        // Note F2 of `docs/formulations/mitc4plus-2017-extract.md`: on a flat
        // RECTANGLE (c_r = c_s = 0) Eq. (27c) must reduce to Eq. (18c)'s
        //   e_rs(E) + 1/2 e_rr|lin r + 1/2 e_ss|lin s,
        // which requires the leading `1` of `(1 + a_E r s)`. Dropping it is
        // caught here (the four-corner extraction above cannot see a constant
        // offset).
        let rect = pre_from(&RECT);
        let (_, _, e_rs_con, e_rr_lin, e_ss_lin, _) = parts(&rect);
        for &(r, s) in &[(0.3, -0.4), (-0.7, 0.2), (0.0, 0.0)] {
            let got = row_dot_3(&b_membrane_covariant_2017(&rect, r, s), 2, &dofs);
            let expected = e_rs_con + 0.5 * e_rr_lin * r + 0.5 * e_ss_lin * s;
            assert!(
                (got - expected).abs() <= 1e-14,
                "flat rectangle Eq. (18c) reduction at ({r},{s}): {got} vs {expected}"
            );
        }
    }

    /// Test-local evaluation of Ko, Lee & Bathe (2017), C&S 182:404-418,
    /// Eqs. (7c)/(7d) with the `x_b` enrichment of Eq. (8a). Returns the
    /// covariant `[e_rr, e_ss, 2 e_rs]` rows of `e^b1` and `e^b2`.
    fn bending_reference(
        pre: &Mitc4PlusDPrecomputed,
        r: f64,
        s: f64,
        include_xb_um: bool,
    ) -> (SMatrix<f64, 3, 24>, SMatrix<f64, 3, 24>) {
        let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
        // x_m tangents (Eq. 9) in local components.
        let mut xm_r = Vector3::zeros();
        let mut xm_s = Vector3::zeros();
        for i in 0..4 {
            xm_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
            xm_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
        }
        let xm_r = local_components(pre, &xm_r);
        let xm_s = local_components(pre, &xm_s);
        // x_b tangents (Eq. 8a) are already in the local frame.
        let mut xb_r = Vector3::zeros();
        let mut xb_s = Vector3::zeros();
        for i in 0..4 {
            let vn_local = local_components(pre, &pre.vn[i]);
            xb_r += 0.5 * pre.a_i[i] * dn_dr[i] * vn_local;
            xb_s += 0.5 * pre.a_i[i] * dn_ds[i] * vn_local;
        }

        let mut b1 = SMatrix::<f64, 3, 24>::zeros();
        let mut b2 = SMatrix::<f64, 3, 24>::zeros();
        for i in 0..4 {
            let um_r = 0.25 * NODE_XI[i] * (1.0 + s * NODE_ETA[i]);
            let um_s = 0.25 * NODE_ETA[i] * (1.0 + r * NODE_XI[i]);
            let ub_r = 0.5 * pre.a_i[i] * dn_dr[i];
            let ub_s = 0.5 * pre.a_i[i] * dn_ds[i];
            let vn = local_components(pre, &pre.vn[i]);
            let cross = [
                [0.0, vn[2], -vn[1]],
                [-vn[2], 0.0, vn[0]],
                [vn[1], -vn[0], 0.0],
            ];
            for k in 0..3 {
                if include_xb_um {
                    b1[(0, 6 * i + k)] += xb_r[k] * um_r;
                    b1[(1, 6 * i + k)] += xb_s[k] * um_s;
                    b1[(2, 6 * i + k)] += 0.5 * (xb_r[k] * um_s + xb_s[k] * um_r);
                }
                for bb in 0..3 {
                    let ubk = ub_r * cross[k][bb];
                    let ubs = ub_s * cross[k][bb];
                    b1[(0, 6 * i + 3 + bb)] += xm_r[k] * ubk;
                    b1[(1, 6 * i + 3 + bb)] += xm_s[k] * ubs;
                    b1[(2, 6 * i + 3 + bb)] += 0.5 * (xm_r[k] * ubs + xm_s[k] * ubk);
                    b2[(0, 6 * i + 3 + bb)] += xb_r[k] * ubk;
                    b2[(1, 6 * i + 3 + bb)] += xb_s[k] * ubs;
                    b2[(2, 6 * i + 3 + bb)] += 0.5 * (xb_r[k] * ubs + xb_s[k] * ubk);
                }
            }
        }
        for j in 0..24 {
            b1[(2, j)] *= 2.0;
            b2[(2, j)] *= 2.0;
        }
        (b1, b2)
    }

    #[test]
    fn test_identity_bending_operator_matches_eq7c_eq7d() {
        let points = [(GP, GP), (GP, -GP), (-GP, GP), (-GP, -GP)];
        for (name, geo, expect_xb_um) in [
            ("flat-rect", RECT, false),
            ("doubly-warped", DOUBLY_WARPED, true),
        ] {
            let pre = pre_from(&geo);

            // The local bending matrices are 3x24 each: no condensed DOF.
            let (bb1, bb2) = b_bending_2017(&pre, 0.0, 0.0);
            assert_eq!(bb1.shape(), (3, 24), "{name}: B_b1 shape");
            assert_eq!(bb2.shape(), (3, 24), "{name}: B_b2 shape");

            let mut max_sep = 0.0f64;
            for &(r, s) in &points {
                let (c1, c2) = b_bending_covariant_2017(&pre, r, s);
                let (r1, r2) = bending_reference(&pre, r, s, true);
                let s1 = scale_3(&r1).max(1e-30);
                let s2 = scale_3(&r2).max(1e-30);
                for k in 0..3 {
                    for j in 0..24 {
                        assert!(
                            (c1[(k, j)] - r1[(k, j)]).abs() <= 1e-10 * s1 + 1e-13,
                            "{name} B_b1[{k}][{j}] at ({r},{s}): {} vs Eq. (7c) {}",
                            c1[(k, j)],
                            r1[(k, j)]
                        );
                        assert!(
                            (c2[(k, j)] - r2[(k, j)]).abs() <= 1e-10 * s2 + 1e-13,
                            "{name} B_b2[{k}][{j}] at ({r},{s}): {} vs Eq. (7d) {}",
                            c2[(k, j)],
                            r2[(k, j)]
                        );
                    }
                }
                // Non-vacuity of the Eq. (8a) `dx_b . du_m` term.
                let (r1_without, _) = bending_reference(&pre, r, s, false);
                for k in 0..3 {
                    for j in 0..24 {
                        max_sep = max_sep.max((r1[(k, j)] - r1_without[(k, j)]).abs());
                    }
                }
            }
            if expect_xb_um {
                assert!(
                    max_sep > 1e-6,
                    "{name}: dx_b . du_m term of Eq. (8a) is absent (separation {max_sep})"
                );
            } else {
                assert!(
                    max_sep <= 1e-14,
                    "{name}: dx_b . du_m must vanish on flat geometry (got {max_sep})"
                );
            }
        }
    }

    /// The standard Mindlin transverse shear of a flat plate in the element's
    /// local frame: `gamma_13 = w_,x + theta_y`, `gamma_23 = w_,y - theta_x`.
    fn mindlin_shear(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64, dofs: &[f64; 24]) -> [f64; 2] {
        let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
        let h = shape_functions(r, s);
        let mut j = [[0.0f64; 2]; 2];
        for i in 0..4 {
            let lx = pre.local_coords[i][0];
            let ly = pre.local_coords[i][1];
            j[0][0] += dn_dr[i] * lx;
            j[0][1] += dn_ds[i] * lx;
            j[1][0] += dn_dr[i] * ly;
            j[1][1] += dn_ds[i] * ly;
        }
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        let inv = [
            [j[1][1] / det, -j[0][1] / det],
            [-j[1][0] / det, j[0][0] / det],
        ];
        let (mut w_x, mut w_y, mut thx, mut thy) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..4 {
            let dn_dx = inv[0][0] * dn_dr[i] + inv[0][1] * dn_ds[i];
            let dn_dy = inv[1][0] * dn_dr[i] + inv[1][1] * dn_ds[i];
            w_x += dn_dx * dofs[6 * i + 2];
            w_y += dn_dy * dofs[6 * i + 2];
            thx += h[i] * dofs[6 * i + 3];
            thy += h[i] * dofs[6 * i + 4];
        }
        [w_x + thy, w_y - thx]
    }

    #[test]
    fn test_shear_mitc4_flat_reduces_to_mindlin_assumed_field() {
        // DB84 Eq. (3) as reproduced in Ko, Lee & Bathe (2017), C&S 182:404-418,
        // p. 405. On a flat element the tying-point values are the standard
        // Mindlin shears gamma_13 = w_,x + theta_y, gamma_23 = w_,y - theta_x.
        let pre = pre_from(&RECT);
        let dofs = generic_dofs();

        for &(r, s, comp) in &[
            (0.0f64, 1.0f64, 0usize),
            (0.0, -1.0, 0),
            (1.0, 0.0, 1),
            (-1.0, 0.0, 1),
        ] {
            let got = row_dot_2(&b_shear_mitc4(&pre, r, s), comp, &dofs);
            let exp = mindlin_shear(&pre, r, s, &dofs)[comp];
            assert!(
                (got - exp).abs() <= 1e-12 * exp.abs().max(1.0),
                "tying point ({r},{s}) component {comp}: MITC4 {got} vs Mindlin {exp}"
            );
        }

        // Non-vacuity: away from the tying points the assumed field is NOT the
        // point-wise Mindlin field, so the check above is not the naive operator.
        let b = b_shear_mitc4(&pre, GP, GP);
        let got = [row_dot_2(&b, 0, &dofs), row_dot_2(&b, 1, &dofs)];
        let exp = mindlin_shear(&pre, GP, GP, &dofs);
        let sep = (got[0] - exp[0]).abs().max((got[1] - exp[1]).abs());
        let scale = got[0]
            .abs()
            .max(got[1].abs())
            .max(exp[0].abs())
            .max(exp[1].abs())
            .max(1e-30);
        assert!(
            sep > 1e-6 * scale,
            "assumed and point-wise shear coincide at a Gauss point ({sep} vs {scale})"
        );
    }

    /// Test-local reference for the drill-membrane operator, written from the
    /// printed Eq. (18) alone (Ko, Bathe & Zhang (2025), C&S 308:107622, p. 8).
    /// It shares no code with `b_drill_membrane_2025`: it recomputes `x_m^l`,
    /// `x_r^l`, `x_s^l`, `j`, `j0`, `theta^D` and the paper's edge order
    /// independently. If it ever disagrees with the production operator, the
    /// resolution is a recorded vision re-read of the paper, never a silent
    /// edit of this reference.
    mod drill {
        use super::*;

        /// Paper node `p` (1..4) -> code node, Ko, Bathe & Zhang (2025),
        /// C&S 308:107622, Eq. (2): paper 1=(+1,+1), 2=(-1,+1), 3=(-1,-1),
        /// 4=(+1,-1).
        const PAPER_TO_CODE: [usize; 4] = [2, 3, 0, 1];

        /// Eq. (13c): the four edge vectors, paper edge order `l = 5..8` =
        /// [right, top, left, bottom]. Entry `(a, b)` is the paper node pair of
        /// edge `l` and `x_m^l = 1/8 (x_a - x_b)`.
        const EDGE_PAPER: [(usize, usize); 4] = [(4, 1), (1, 2), (2, 3), (3, 4)];

        /// Mid-point `(r, s)` of each edge `l = 5..8`.
        const EDGE_MID: [(f64, f64); 4] = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)];

        fn shape(r: f64, s: f64) -> [f64; 4] {
            [
                0.25 * (1.0 - r) * (1.0 - s),
                0.25 * (1.0 + r) * (1.0 - s),
                0.25 * (1.0 + r) * (1.0 + s),
                0.25 * (1.0 - r) * (1.0 + s),
            ]
        }

        fn dshape(r: f64, s: f64) -> ([f64; 4], [f64; 4]) {
            (
                [
                    -0.25 * (1.0 - s),
                    0.25 * (1.0 - s),
                    0.25 * (1.0 + s),
                    -0.25 * (1.0 + s),
                ],
                [
                    -0.25 * (1.0 - r),
                    -0.25 * (1.0 + r),
                    0.25 * (1.0 + r),
                    0.25 * (1.0 - r),
                ],
            )
        }

        /// The exact `j(r,s) = det[g_r g_s g_t]|_{(r,s,0)}` of Eq. (17b).
        fn jacobian(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> f64 {
            let (dn_dr, dn_ds) = dshape(r, s);
            let h = shape(r, s);
            let mut g_r = Vector3::zeros();
            let mut g_s = Vector3::zeros();
            let mut g_t = Vector3::zeros();
            for i in 0..4 {
                g_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
                g_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
                g_t += 0.5 * pre.a_i[i] * h[i] * pre.vn[i];
            }
            g_r.cross(&g_s).dot(&g_t)
        }

        /// The `V^D` of Eq. (5), recomputed from the geometry.
        pub fn drill_normal(pre: &Mitc4PlusDPrecomputed) -> Vector3<f64> {
            let (dn_dr, dn_ds) = dshape(0.0, 0.0);
            let mut x_r = Vector3::zeros();
            let mut x_s = Vector3::zeros();
            for i in 0..4 {
                x_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
                x_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
            }
            x_r.cross(&x_s).normalize()
        }

        /// The exact 3x3 covariant-to-local map of Eq. (21): the constant
        /// element-centre base vectors.
        fn centre_map(
            pre: &Mitc4PlusDPrecomputed,
            b_cov: &SMatrix<f64, 3, 24>,
        ) -> SMatrix<f64, 3, 24> {
            let (dn_dr, dn_ds) = dshape(0.0, 0.0);
            let mut g_r = Vector3::zeros();
            let mut g_s = Vector3::zeros();
            for i in 0..4 {
                g_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
                g_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
            }
            let (j00, j01) = (g_r.dot(&pre.e1), g_s.dot(&pre.e1));
            let (j10, j11) = (g_r.dot(&pre.e2), g_s.dot(&pre.e2));
            // The same de-singularising shift the production mapping uses; a
            // numerical guard, not a formulation term.
            let scale = j00
                .abs()
                .max(j01.abs())
                .max(j10.abs())
                .max(j11.abs())
                .max(1.0);
            let eps = scale * 1.0e-12;
            let (a, b, c, d) = (j00 + eps, j01, j10, j11 + eps);
            let det = a * d - b * c;
            let (i11, i12, i21, i22) = (d / det, -b / det, -c / det, a / det);
            let t = [
                [i11 * i11, i21 * i21, i11 * i21],
                [i12 * i12, i22 * i22, i12 * i22],
                [2.0 * i11 * i12, 2.0 * i21 * i22, i11 * i22 + i12 * i21],
            ];
            let mut out = SMatrix::<f64, 3, 24>::zeros();
            for i in 0..3 {
                for j in 0..24 {
                    out[(i, j)] =
                        t[i][0] * b_cov[(0, j)] + t[i][1] * b_cov[(1, j)] + t[i][2] * b_cov[(2, j)];
                }
            }
            out
        }

        /// The clean reference: Eq. (18) verbatim.
        pub fn b_md_reference(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 3, 24> {
            let vd = drill_normal(pre);
            let vd_local = local_components(pre, &vd);
            b_md_parameterised(pre, r, s, &vd, &vd_local, 0.125, &[0, 1, 2, 3], 1.0)
        }

        /// The reference with one ingredient deliberately wrong, so the five
        /// known-risk deviations are shown to be rejected: `xm_scale` (the 1/8
        /// of Eq. 13c), `edge_perm` (the `h~` pairing order), the
        /// edge-difference `sign`, `v_d` (the vector in `c_r`/`c_s`) and
        /// `v_theta` (the vector of `theta^D = theta . v_theta`).
        #[allow(clippy::too_many_arguments)]
        pub fn b_md_parameterised(
            pre: &Mitc4PlusDPrecomputed,
            r: f64,
            s: f64,
            v_d: &Vector3<f64>,
            v_theta: &Vector3<f64>,
            xm_scale: f64,
            edge_perm: &[usize; 4],
            sign: f64,
        ) -> SMatrix<f64, 3, 24> {
            let j = jacobian(pre, r, s);
            let j0 = jacobian(pre, 0.0, 0.0);
            let ratio = if j.abs() > 1e-30 { j0 / j } else { 0.0 };

            // Eq. (11a)/(11b): the four "curl" derivatives, paper order
            // [right, top, left, bottom].
            let dh = [
                (0.0, -s * (1.0 + r)),
                (-r * (1.0 + s), 0.0),
                (0.0, -s * (1.0 - r)),
                (-r * (1.0 - s), 0.0),
            ];

            let mut b_cov = SMatrix::<f64, 3, 24>::zeros();
            for l in 0..4 {
                let (a, b) = EDGE_PAPER[l];
                let icode = PAPER_TO_CODE[a - 1];
                let kcode = PAPER_TO_CODE[b - 1];
                // Eq. (13c): x_m^l = 1/8 (x_i - x_{i+1}).
                let x_m = xm_scale
                    * (node_vec(&pre.initial_coords_3d, icode)
                        - node_vec(&pre.initial_coords_3d, kcode));
                // Eq. (13b): the mid-surface tangents at the edge mid-point.
                let (rm, sm) = EDGE_MID[l];
                let (dn_dr, dn_ds) = dshape(rm, sm);
                let mut x_r = Vector3::zeros();
                let mut x_s = Vector3::zeros();
                for i in 0..4 {
                    x_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
                    x_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
                }
                let c_r = x_m.dot(&(-x_r.cross(v_d)));
                let c_s = x_m.dot(&(x_s.cross(v_d)));
                let (h_r, h_s) = dh[edge_perm[l]];
                let coeff_rr = ratio * h_r * c_r;
                let coeff_ss = -ratio * h_s * c_s;
                let coeff_rs = 0.5 * ratio * (h_s * c_r - h_r * c_s);
                for beta in 0..3 {
                    let vt = v_theta[beta];
                    b_cov[(0, 6 * icode + 3 + beta)] -= sign * coeff_rr * vt;
                    b_cov[(0, 6 * kcode + 3 + beta)] += sign * coeff_rr * vt;
                    b_cov[(1, 6 * icode + 3 + beta)] -= sign * coeff_ss * vt;
                    b_cov[(1, 6 * kcode + 3 + beta)] += sign * coeff_ss * vt;
                    b_cov[(2, 6 * icode + 3 + beta)] -= sign * 2.0 * coeff_rs * vt;
                    b_cov[(2, 6 * kcode + 3 + beta)] += sign * 2.0 * coeff_rs * vt;
                }
            }
            centre_map(pre, &b_cov)
        }
    }

    #[test]
    fn test_identity_drill_operator_matches_eq18_term_by_term() {
        // Oracle 1: the production Eq. (18) operator vs the independent
        // reference, entry by entry at 9 sample points on four quads.
        let points = [
            (GP, GP),
            (GP, -GP),
            (-GP, GP),
            (-GP, -GP),
            (0.0, 1.0),
            (1.0, 0.0),
            (0.0, -1.0),
            (-1.0, 0.0),
            (0.0, 0.0),
        ];
        let flat_square = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        for (name, geo) in [
            ("flat-square", flat_square),
            ("flat-distorted", FLAT_DISTORTED),
            ("ruled-warped", RULED_WARPED),
            ("doubly-warped", DOUBLY_WARPED),
        ] {
            let pre = pre_from(&geo);
            for &(r, s) in &points {
                let prod = b_drill_membrane_2025(&pre, r, s);
                let refr = drill::b_md_reference(&pre, r, s);
                for i in 0..3 {
                    for j in 0..24 {
                        assert!(
                            (prod[(i, j)] - refr[(i, j)]).abs() <= 1e-12,
                            "{name} ({r},{s}) entry [{i}][{j}]: production {} vs Eq. (18) {} ",
                            prod[(i, j)],
                            refr[(i, j)]
                        );
                    }
                }
            }
        }

        // Oracle 2: the five known-risk deviations are asserted to be rejected.
        let pre = pre_from(&DOUBLY_WARPED);
        let vd = drill::drill_normal(&pre);
        let e3 = Vector3::new(0.0, 0.0, 1.0);
        let refr = drill::b_md_reference(&pre, GP, GP);
        let vd_local = local_components(&pre, &vd);
        let ref_scale = scale_3(&refr).max(1e-30);
        let rel = |m: &SMatrix<f64, 3, 24>| {
            let mut d = 0.0f64;
            for i in 0..3 {
                for j in 0..24 {
                    d = d.max((m[(i, j)] - refr[(i, j)]).abs());
                }
            }
            d / ref_scale
        };
        let variants: [(&str, SMatrix<f64, 3, 24>, f64); 5] = [
            (
                "V^D -> e3",
                drill::b_md_parameterised(&pre, GP, GP, &e3, &e3, 0.125, &[0, 1, 2, 3], 1.0),
                1e-6,
            ),
            (
                "missing 1/8",
                drill::b_md_parameterised(&pre, GP, GP, &vd, &vd_local, 1.0, &[0, 1, 2, 3], 1.0),
                1e-3,
            ),
            (
                "edge order [bottom,right,top,left]",
                drill::b_md_parameterised(&pre, GP, GP, &vd, &vd_local, 0.125, &[1, 2, 3, 0], 1.0),
                1e-6,
            ),
            (
                "flipped edge difference",
                drill::b_md_parameterised(&pre, GP, GP, &vd, &vd_local, 0.125, &[0, 1, 2, 3], -1.0),
                1e-6,
            ),
            (
                "theta_z alone",
                drill::b_md_parameterised(&pre, GP, GP, &vd, &e3, 0.125, &[0, 1, 2, 3], 1.0),
                1e-6,
            ),
        ];
        for (name, variant, margin) in variants {
            assert!(
                rel(&variant) > margin,
                "rejection '{name}': not distinguished from Eq. (18) (relative {})",
                rel(&variant)
            );
        }
    }

    #[test]
    fn test_identity_cr_cs_2017_and_2025_are_different_quantities() {
        // The two papers reuse the symbols c_r, c_s for different quantities:
        //  - Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (25): c_r = x_d . m^r.
        //  - Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (18)/(19c):
        //    c_r^l = x_m^l . (-x_r^l x V^D).
        // Conflating them would silently corrupt the drill operator.
        let pre = pre_from(&DOUBLY_WARPED);

        // The 2017 quantity is the stored x_d . m^r / x_d . m^s.
        assert!((pre.c_r_mem - pre.x_d.dot(&pre.m_r)).abs() <= 1e-14);
        assert!((pre.c_s_mem - pre.x_d.dot(&pre.m_s)).abs() <= 1e-14);

        // The 2025 quantity, recomputed independently for each edge.
        let mut min_rel = f64::INFINITY;
        for e in 0..4 {
            let term = pre.drill_edges[e];
            let x_m = (node_vec(&pre.initial_coords_3d, term.start)
                - node_vec(&pre.initial_coords_3d, term.end))
                / 8.0;
            let (rm, sm) = DRILL_EDGE_MID[e];
            let (dn_dr, dn_ds) = shape_function_derivatives(rm, sm);
            let mut x_r = Vector3::zeros();
            let mut x_s = Vector3::zeros();
            for i in 0..4 {
                x_r += dn_dr[i] * node_vec(&pre.initial_coords_3d, i);
                x_s += dn_ds[i] * node_vec(&pre.initial_coords_3d, i);
            }
            let cr_md = x_m.dot(&(-x_r.cross(&pre.v_d)));
            let cs_md = x_m.dot(&(x_s.cross(&pre.v_d)));
            assert!(
                (term.c_r - cr_md).abs() <= 1e-14 * cr_md.abs().max(1.0),
                "edge {e} stored c_r {} != Eq. (18) {cr_md}",
                term.c_r
            );
            assert!(
                (term.c_s - cs_md).abs() <= 1e-14 * cs_md.abs().max(1.0),
                "edge {e} stored c_s {} != Eq. (18) {cs_md}",
                term.c_s
            );

            let scale_r = cr_md.abs().max(pre.c_r_mem.abs()).max(1e-30);
            let scale_s = cs_md.abs().max(pre.c_s_mem.abs()).max(1e-30);
            min_rel = min_rel
                .min((cr_md - pre.c_r_mem).abs() / scale_r)
                .min((cs_md - pre.c_s_mem).abs() / scale_s);
        }
        assert!(
            min_rel > 1e-6,
            "2017 and 2025 c_r/c_s are not distinguished (min relative diff {min_rel})"
        );
    }

    // ========================================================================
    // WU4 — the stiffness assembly (tasks 5.1-5.8).
    //
    // The identity lock (5.2, 5.5) compares the production assembly with an
    // INDEPENDENT test-local reference (`mod ke_ref`) that does not call
    // `compute_ke_local`/`compute_ke_global`: it recomputes the characteristic
    // vectors, the dual basis, the directors, `V^D`, `j0`, the membrane and
    // shear operators, the moment matrix `W` and the Gauss summation itself. It
    // reuses only the already-independent test-local `bending_reference` (task
    // 4.2) and `drill::b_md_reference` (task 4.4), which are themselves locked
    // to the printed equations, and the repository's local-frame convention
    // (`e1`/`e2`/`e3`, which cancels in the global matrix).
    // ========================================================================

    /// Independent test-local reference for the stiffness assembly.
    mod ke_ref {
        use super::*;
        use nalgebra::Matrix2;

        /// The bilinear node coordinates in the repository order.
        const XI: [f64; 4] = [-1.0, 1.0, 1.0, -1.0];
        const ETA: [f64; 4] = [-1.0, -1.0, 1.0, 1.0];

        fn n(r: f64, s: f64) -> [f64; 4] {
            [
                0.25 * (1.0 - r) * (1.0 - s),
                0.25 * (1.0 + r) * (1.0 - s),
                0.25 * (1.0 + r) * (1.0 + s),
                0.25 * (1.0 - r) * (1.0 + s),
            ]
        }

        fn dn(r: f64, s: f64) -> ([f64; 4], [f64; 4]) {
            (
                [
                    -0.25 * (1.0 - s),
                    0.25 * (1.0 - s),
                    0.25 * (1.0 + s),
                    -0.25 * (1.0 + s),
                ],
                [
                    -0.25 * (1.0 - r),
                    -0.25 * (1.0 + r),
                    0.25 * (1.0 + r),
                    0.25 * (1.0 - r),
                ],
            )
        }

        /// Independently recomputed geometry: `x_r`, `x_s`, `x_d`, the plane
        /// normal, the dual basis (cross-product form, not the production Gram
        /// solve), the directors (ADR-4 option B, re-derived), `V^D`, `j0` and
        /// the five 2017 coefficients.
        pub(super) struct Geom {
            x_r: Vector3<f64>,
            x_s: Vector3<f64>,
            x_d: Vector3<f64>,
            vn: [Vector3<f64>; 4],
            j0: f64,
            a: [f64; 5],
        }

        pub(super) fn geom(pre: &Mitc4PlusDPrecomputed) -> Geom {
            let x: [Vector3<f64>; 4] = std::array::from_fn(|i| node_vec(&pre.initial_coords_3d, i));
            let x_r = 0.25 * (-x[0] + x[1] + x[2] - x[3]);
            let x_s = 0.25 * (-x[0] - x[1] + x[2] + x[3]);
            let x_d = 0.25 * (x[0] - x[1] + x[2] - x[3]);
            let nv = x_r.cross(&x_s).normalize();
            // Dual basis by cross products (independent of the production solve).
            let m_r = x_s.cross(&nv) / x_r.dot(&x_s.cross(&nv));
            let m_s = nv.cross(&x_r) / x_s.dot(&nv.cross(&x_r));
            // Directors: the same ADR-4 option B definition, re-derived here.
            let c = 0.25 * (x[0] + x[1] + x[2] + x[3]);
            let mut tri = [Vector3::zeros(); 4];
            for k in 0..4 {
                let mut t = (x[k] - c).cross(&(x[(k + 1) % 4] - c));
                if t.dot(&nv) < 0.0 {
                    t = -t;
                }
                tri[k] = t;
            }
            let vn: [Vector3<f64>; 4] =
                std::array::from_fn(|i| (tri[(i + 3) % 4] + tri[i]).normalize());
            let mut g_t0 = Vector3::zeros();
            for (a, v) in pre.a_i.iter().zip(vn.iter()) {
                g_t0 += 0.5 * a * v;
            }
            let j0 = x_r.cross(&x_s).dot(&g_t0);
            let c_r = x_d.dot(&m_r);
            let c_s = x_d.dot(&m_s);
            let d = c_r * c_r + c_s * c_s - 1.0;
            let a = [
                c_r * (c_r - 1.0) / (2.0 * d),
                c_r * (c_r + 1.0) / (2.0 * d),
                c_s * (c_s - 1.0) / (2.0 * d),
                c_s * (c_s + 1.0) / (2.0 * d),
                2.0 * c_r * c_s / d,
            ];
            Geom {
                x_r,
                x_s,
                x_d,
                vn,
                j0,
                a,
            }
        }

        /// Point-wise covariant-to-local map, recomputed independently.
        fn map_local(
            pre: &Mitc4PlusDPrecomputed,
            g: &Geom,
            r: f64,
            s: f64,
            cov: &SMatrix<f64, 3, 24>,
        ) -> SMatrix<f64, 3, 24> {
            let gr = g.x_r + s * g.x_d;
            let gs = g.x_s + r * g.x_d;
            let j00 = gr.dot(&pre.e1);
            let j01 = gs.dot(&pre.e1);
            let j10 = gr.dot(&pre.e2);
            let j11 = gs.dot(&pre.e2);
            let det = j00 * j11 - j01 * j10;
            let (i11, i12, i21, i22) = (j11 / det, -j01 / det, -j10 / det, j00 / det);
            let t = [
                [i11 * i11, i21 * i21, i11 * i21],
                [i12 * i12, i22 * i22, i12 * i22],
                [2.0 * i11 * i12, 2.0 * i21 * i22, i11 * i22 + i12 * i21],
            ];
            let mut out = SMatrix::<f64, 3, 24>::zeros();
            for i in 0..3 {
                for j in 0..24 {
                    out[(i, j)] =
                        t[i][0] * cov[(0, j)] + t[i][1] * cov[(1, j)] + t[i][2] * cov[(2, j)];
                }
            }
            out
        }

        /// The MITC4+ assumed membrane operator of Ko, Lee & Bathe (2017),
        /// C&S 182:404-418, Eqs. (15)-(27), re-derived from the printed field.
        pub(super) fn membrane_local(
            pre: &Mitc4PlusDPrecomputed,
            g: &Geom,
            r: f64,
            s: f64,
        ) -> SMatrix<f64, 3, 24> {
            let lc =
                |v: &Vector3<f64>| Vector3::new(v.dot(&pre.e1), v.dot(&pre.e2), v.dot(&pre.e3));
            let row_rr = |ss: f64| -> [f64; 24] {
                let gr = lc(&(g.x_r + ss * g.x_d));
                let mut b = [0.0f64; 24];
                for i in 0..4 {
                    let coef = 0.25 * XI[i] * (1.0 + ss * ETA[i]);
                    for k in 0..3 {
                        b[6 * i + k] = coef * gr[k];
                    }
                }
                b
            };
            let row_ss = |rr: f64| -> [f64; 24] {
                let gs = lc(&(g.x_s + rr * g.x_d));
                let mut b = [0.0f64; 24];
                for i in 0..4 {
                    let coef = 0.25 * ETA[i] * (1.0 + rr * XI[i]);
                    for k in 0..3 {
                        b[6 * i + k] = coef * gs[k];
                    }
                }
                b
            };
            let row_rs = |rr: f64, ss: f64| -> [f64; 24] {
                let gr = lc(&(g.x_r + ss * g.x_d));
                let gs = lc(&(g.x_s + rr * g.x_d));
                let mut b = [0.0f64; 24];
                for i in 0..4 {
                    let cr = 0.25 * XI[i] * (1.0 + ss * ETA[i]);
                    let cs = 0.25 * ETA[i] * (1.0 + rr * XI[i]);
                    for k in 0..3 {
                        b[6 * i + k] = 0.5 * (cs * gr[k] + cr * gs[k]);
                    }
                }
                b
            };
            let e_a = row_rr(1.0);
            let e_b = row_rr(-1.0);
            let e_c = row_ss(1.0);
            let e_d = row_ss(-1.0);
            let e_e = row_rs(0.0, 0.0);
            let a = g.a;
            let mut cov = SMatrix::<f64, 3, 24>::zeros();
            for j in 0..24 {
                cov[(0, j)] = 0.5 * (1.0 - 2.0 * a[0] + s + 2.0 * a[0] * s * s) * e_a[j]
                    + 0.5 * (1.0 - 2.0 * a[1] - s + 2.0 * a[1] * s * s) * e_b[j]
                    + a[2] * (-1.0 + s * s) * e_c[j]
                    + a[3] * (-1.0 + s * s) * e_d[j]
                    + a[4] * (-1.0 + s * s) * e_e[j];
                cov[(1, j)] = a[0] * (-1.0 + r * r) * e_a[j]
                    + a[1] * (-1.0 + r * r) * e_b[j]
                    + 0.5 * (1.0 - 2.0 * a[2] + r + 2.0 * a[2] * r * r) * e_c[j]
                    + 0.5 * (1.0 - 2.0 * a[3] - r + 2.0 * a[3] * r * r) * e_d[j]
                    + a[4] * (-1.0 + r * r) * e_e[j];
                cov[(2, j)] = 0.25 * (r + 4.0 * a[0] * r * s) * e_a[j]
                    + 0.25 * (-r + 4.0 * a[1] * r * s) * e_b[j]
                    + 0.25 * (s + 4.0 * a[2] * r * s) * e_c[j]
                    + 0.25 * (-s + 4.0 * a[3] * r * s) * e_d[j]
                    + (1.0 + a[4] * r * s) * e_e[j];
            }
            for j in 0..24 {
                cov[(2, j)] *= 2.0;
            }
            map_local(pre, g, r, s, &cov)
        }

        /// The MITC4 assumed transverse-shear operator of Dvorkin & Bathe
        /// (1984), Engineering Computations 1:77-88, Eq. (3), re-derived.
        pub(super) fn shear_local(
            pre: &Mitc4PlusDPrecomputed,
            g: &Geom,
            r: f64,
            s: f64,
        ) -> SMatrix<f64, 2, 24> {
            let tie = |rr: f64, ss: f64| -> SMatrix<f64, 2, 24> {
                let (dnr, dns) = dn(rr, ss);
                let h = n(rr, ss);
                let mut gr = Vector3::zeros();
                let mut gs = Vector3::zeros();
                let mut gt = Vector3::zeros();
                for i in 0..4 {
                    gr += dnr[i] * node_vec(&pre.initial_coords_3d, i);
                    gs += dns[i] * node_vec(&pre.initial_coords_3d, i);
                    gt += 0.5 * pre.a_i[i] * h[i] * g.vn[i];
                }
                let gr = Vector3::new(gr.dot(&pre.e1), gr.dot(&pre.e2), gr.dot(&pre.e3));
                let gs = Vector3::new(gs.dot(&pre.e1), gs.dot(&pre.e2), gs.dot(&pre.e3));
                let gt = Vector3::new(gt.dot(&pre.e1), gt.dot(&pre.e2), gt.dot(&pre.e3));
                let mut b = SMatrix::<f64, 2, 24>::zeros();
                for i in 0..4 {
                    let vn_local = Vector3::new(
                        g.vn[i].dot(&pre.e1),
                        g.vn[i].dot(&pre.e2),
                        g.vn[i].dot(&pre.e3),
                    );
                    let vxg_r = vn_local.cross(&gr);
                    let vxg_s = vn_local.cross(&gs);
                    for k in 0..3 {
                        b[(0, 6 * i + k)] = 0.5 * gt[k] * dnr[i];
                        b[(1, 6 * i + k)] = 0.5 * gt[k] * dns[i];
                        b[(0, 6 * i + 3 + k)] = 0.25 * pre.a_i[i] * h[i] * vxg_r[k];
                        b[(1, 6 * i + 3 + k)] = 0.25 * pre.a_i[i] * h[i] * vxg_s[k];
                    }
                }
                b
            };
            let ba = tie(0.0, 1.0);
            let bb = tie(0.0, -1.0);
            let bc = tie(1.0, 0.0);
            let bd = tie(-1.0, 0.0);
            let mut cov = SMatrix::<f64, 2, 24>::zeros();
            for j in 0..24 {
                cov[(0, j)] = 0.5 * (1.0 + s) * ba[(0, j)] + 0.5 * (1.0 - s) * bb[(0, j)];
                cov[(1, j)] = 0.5 * (1.0 + r) * bc[(1, j)] + 0.5 * (1.0 - r) * bd[(1, j)];
            }
            // 3D dual-basis metric normalization.
            let gr = g.x_r + s * g.x_d;
            let gs = g.x_s + r * g.x_d;
            let h = n(r, s);
            let mut gt = Vector3::zeros();
            for ((a, hi), vn) in pre.a_i.iter().zip(h.iter()).zip(g.vn.iter()) {
                gt += 0.5 * a * hi * vn;
            }
            let jdet = gr.cross(&gs).dot(&gt);
            let grd = gs.cross(&gt) / jdet;
            let gsd = gt.cross(&gr) / jdet;
            let gtd = gr.cross(&gs) / jdet;
            let ft = gtd.dot(&pre.e3);
            let t = Matrix2::new(
                2.0 * grd.dot(&pre.e1) * ft,
                2.0 * gsd.dot(&pre.e1) * ft,
                2.0 * grd.dot(&pre.e2) * ft,
                2.0 * gsd.dot(&pre.e2) * ft,
            );
            t * cov
        }

        /// The resultant/moment matrix `W`, recomputed independently.
        fn moment_matrix(c: &ShellConstitutive, include_w22: bool) -> SMatrix<f64, 9, 9> {
            let mut w = SMatrix::<f64, 9, 9>::zeros();
            for i in 0..3 {
                for j in 0..3 {
                    w[(i, j)] = c.cm[(i, j)];
                    w[(i, 3 + j)] = c.cb_coupling[(i, j)];
                    w[(3 + i, j)] = c.cb_coupling[(j, i)];
                    w[(i, 6 + j)] = c.cb[(i, j)];
                    w[(6 + i, j)] = c.cb[(j, i)];
                    w[(3 + i, 3 + j)] = c.cb[(i, j)];
                    if include_w22 {
                        w[(6 + i, 6 + j)] = c.cm[(i, j)] / 9.0;
                    }
                }
            }
            w
        }

        fn assemble(pre: &Mitc4PlusDPrecomputed, g: &Geom, include_w22: bool) -> Mat24 {
            let w = moment_matrix(&pre.constitutive, include_w22);
            let s1 = 2.0 / pre.thickness;
            let s2 = 4.0 / (pre.thickness * pre.thickness);
            let mut k = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let gr = g.x_r + s * g.x_d;
                let gs = g.x_s + r * g.x_d;
                let sqrt_g = gr.cross(&gs).norm();
                let bm = membrane_local(pre, g, r, s);
                let (c1, c2) = bending_reference(pre, r, s, true);
                let b1 = map_local(pre, g, r, s, &c1);
                let b2 = map_local(pre, g, r, s, &c2);
                let mut b = SMatrix::<f64, 9, 24>::zeros();
                for i in 0..3 {
                    for j in 0..24 {
                        b[(i, j)] = bm[(i, j)];
                        b[(3 + i, j)] = s1 * b1[(i, j)];
                        b[(6 + i, j)] = s2 * b2[(i, j)];
                    }
                }
                k += (b.transpose() * w * b) * (GAUSS_W[gi] * sqrt_g);
            }
            k += shear_local_block(pre, g);
            k += drill_block(pre, g);
            k
        }

        fn shear_local_block(pre: &Mitc4PlusDPrecomputed, g: &Geom) -> Mat24 {
            let mut k = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let gr = g.x_r + s * g.x_d;
                let gs = g.x_s + r * g.x_d;
                let sqrt_g = gr.cross(&gs).norm();
                let bg = shear_local(pre, g, r, s);
                k += (bg.transpose() * pre.cs_uncorrected * bg) * (GAUSS_W[gi] * sqrt_g);
            }
            k
        }

        fn drill_block(pre: &Mitc4PlusDPrecomputed, g: &Geom) -> Mat24 {
            let cm = &pre.constitutive.cm;
            let mut m = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let gr = g.x_r + s * g.x_d;
                let gs = g.x_s + r * g.x_d;
                let sqrt_g = gr.cross(&gs).norm();
                let bmd = drill::b_md_reference(pre, r, s);
                m += (bmd.transpose() * cm * bmd) * (GAUSS_W[gi] * sqrt_g);
            }
            m
        }

        /// The full independent reference stiffness (local, 24x24).
        pub fn ke_local(pre: &Mitc4PlusDPrecomputed, include_w22: bool) -> Mat24 {
            let g = geom(pre);
            let _ = g.j0;
            assemble(pre, &g, include_w22)
        }

        /// The independent reference membrane block.
        pub fn membrane_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            let g = geom(pre);
            let cm = &pre.constitutive.cm;
            let mut k = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let gr = g.x_r + s * g.x_d;
                let gs = g.x_s + r * g.x_d;
                let sqrt_g = gr.cross(&gs).norm();
                let bm = membrane_local(pre, &g, r, s);
                k += (bm.transpose() * cm * bm) * (GAUSS_W[gi] * sqrt_g);
            }
            k
        }

        /// The independent reference transverse-shear block.
        pub fn shear_ke_local(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            let g = geom(pre);
            shear_local_block(pre, &g)
        }
    }

    /// The largest-magnitude eigenvalue of a symmetric matrix (the `λ_max`
    /// scale the spec's tolerances use).
    fn lambda_max(k: &Mat24) -> f64 {
        let mut sym = *k;
        for i in 0..24 {
            for j in 0..24 {
                sym[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
            }
        }
        sym.symmetric_eigenvalues()
            .iter()
            .fold(0.0f64, |m, &v| m.max(v.abs()))
    }

    /// The quadratic form `u^T K u`.
    fn quadratic(u: &[f64; 24], k: &Mat24) -> f64 {
        let mut s = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                s += u[i] * k[(i, j)] * u[j];
            }
        }
        s
    }

    fn norm2(u: &[f64; 24]) -> f64 {
        u.iter().map(|v| v * v).sum()
    }

    /// The six rigid-body fields in the element's local 6-DOF layout:
    /// translations `t_k` and rotations `omega_k` about the local axes through
    /// the origin, with `u_i = omega x x_i` and `theta_i = omega`.
    fn rigid_body_fields(pre: &Mitc4PlusDPrecomputed) -> [[f64; 24]; 6] {
        let mut out = [[0.0f64; 24]; 6];
        for d in 0..3 {
            for i in 0..4 {
                out[d][6 * i + d] = 1.0;
            }
        }
        for a in 0..3 {
            let mut omega = [0.0f64; 3];
            omega[a] = 1.0;
            for i in 0..4 {
                let x = [
                    pre.local_coords[i][0],
                    pre.local_coords[i][1],
                    node_vec(&pre.initial_coords_3d, i).dot(&pre.e3),
                ];
                let u = [
                    omega[1] * x[2] - omega[2] * x[1],
                    omega[2] * x[0] - omega[0] * x[2],
                    omega[0] * x[1] - omega[1] * x[0],
                ];
                for k in 0..3 {
                    out[3 + a][6 * i + k] = u[k];
                    out[3 + a][6 * i + 3 + k] = omega[k];
                }
            }
        }
        out
    }

    fn single_ply_lam(h: f64, k: f64) -> Laminate {
        let e = 2.0e11;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let mat = OrthotropicMaterial::new(e, e, e, g, g, g, nu, nu, nu, 2700.0);
        Laminate::new(vec![Ply::new(mat, h, 0.0)], k).unwrap()
    }

    fn max_abs_diff(a: &Mat24, b: &Mat24) -> f64 {
        let mut m = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                m = m.max((a[(i, j)] - b[(i, j)]).abs());
            }
        }
        m
    }

    fn max_abs(k: &Mat24) -> f64 {
        let mut m = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                m = m.max(k[(i, j)].abs());
            }
        }
        m
    }

    // ------------------------------------------------------------------
    // 5.1 — the resultant/moment matrix
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_resultant_moment_matrix_blocks_match_closed_forms() {
        let h = 0.02;
        let e = 2.0e11;
        let nu = 0.3;
        let c = IsotropicMaterial::new(e, nu, 7800.0).constitutive(h, 5.0 / 6.0);
        let w = resultant_moment_matrix(&c);
        assert_eq!(w.shape(), (9, 9));
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(w[(i, j)], c.cm[(i, j)], "W_00 = cm");
                assert_eq!(w[(3 + i, 3 + j)], c.cb[(i, j)], "W_11 = cb");
                assert_eq!(w[(i, 3 + j)], c.cb_coupling[(i, j)], "W_01 = cb_coupling");
                assert_eq!(w[(3 + i, j)], c.cb_coupling[(j, i)], "W_10 = cb_coupling^T");
                assert_eq!(w[(i, 6 + j)], c.cb[(i, j)], "W_02 = cb");
                assert_eq!(w[(6 + i, j)], c.cb[(j, i)], "W_20 = cb^T");
                assert_eq!(w[(6 + i, 6 + j)], c.cm[(i, j)] / 9.0, "W_22 = cm/9");
                assert_eq!(w[(3 + i, 6 + j)], 0.0, "W_12 = 0 (odd moment)");
                assert_eq!(w[(6 + i, 3 + j)], 0.0, "W_21 = 0 (odd moment)");
            }
        }
        for i in 0..9 {
            for j in 0..9 {
                assert_eq!(w[(i, j)], w[(j, i)], "W must be exactly symmetric");
            }
        }
        // Isotropic closed forms (Ko, Lee & Bathe (2017), C&S 182:404-418,
        // Eq. (7a)'s moments).
        let cm00 = e * h / (1.0 - nu * nu);
        let cb00 = e * h * h * h / (12.0 * (1.0 - nu * nu));
        assert!((w[(0, 0)] - cm00).abs() <= 1e-12 * cm00, "W_00 closed form");
        assert!((w[(3, 3)] - cb00).abs() <= 1e-12 * cb00, "W_11 closed form");
        assert!(
            (w[(6, 6)] - cm00 / 9.0).abs() <= 1e-12 * cm00,
            "W_22 = cm/9"
        );
        // The paper's 2x2 t-rule gives int t^4 dt = 2/9 (not the exact 2/5).
        let t = 1.0 / 3.0f64.sqrt();
        let t4 = 2.0 * t.powi(4);
        assert!((t4 - 2.0 / 9.0).abs() <= 1e-15, "2-point t-rule gives 2/9");
        assert!(
            (w[(6, 6)] - cm00 / 5.0).abs() > 1e-3 * (cm00 / 5.0),
            "W_22 must be the paper's cm/9, not the exact cm/5"
        );
    }

    // ------------------------------------------------------------------
    // 5.2 — the identity lock
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_ke_lock_matches_2017_core_plus_2025_drill() {
        for (name, geo) in [
            ("flat-square", FLAT_SQUARE),
            ("flat-distorted", FLAT_DISTORTED),
            ("ruled-warped", RULED_WARPED),
            ("doubly-warped", DOUBLY_WARPED),
        ] {
            let pre = pre_from(&geo);
            let k_prod = compute_ke_local(&pre);
            let k_ref = ke_ref::ke_local(&pre, true);
            assert_eq!(k_prod.shape(), (24, 24), "{name}: K shape");
            let scale = max_abs(&k_ref);
            assert!(scale > 0.0, "{name}: reference is zero");
            // Membrane block.
            let km_prod = membrane_ke_local(&pre);
            let km_ref = ke_ref::membrane_ke_local(&pre);
            let sm = max_abs(&km_ref).max(1e-30);
            let dm = max_abs_diff(&km_prod, &km_ref);
            assert!(
                dm <= 1e-10 * sm,
                "{name}: membrane block diff {dm}, scale {sm}"
            );
            // Transverse-shear block.
            let ks_prod = shear_ke_local(&pre);
            let ks_ref = ke_ref::shear_ke_local(&pre);
            let ss = max_abs(&ks_ref).max(1e-30);
            let ds = max_abs_diff(&ks_prod, &ks_ref);
            assert!(
                ds <= 1e-10 * ss,
                "{name}: shear block diff {ds}, scale {ss}"
            );
            let d = max_abs_diff(&k_prod, &k_ref);
            assert!(
                d <= 1e-10 * scale,
                "{name}: max|K_prod - K_ref| = {d}, scale {scale}"
            );
        }
    }

    // ------------------------------------------------------------------
    // 5.3 — drill stiffness provenance
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_drill_stiffness_comes_only_from_eq26() {
        // (a) Warped geometry: the drill operator is live.
        let warped = pre_from(&DOUBLY_WARPED);
        let k_full = compute_ke_local_with_drill(&warped, true);
        let k_zero = compute_ke_local_with_drill(&warped, false);
        let diff = k_full - k_zero;
        let scale = max_abs(&k_full).max(1e-30);
        let dmax = max_abs(&diff);
        assert!(
            dmax > 1e-10 * scale,
            "drill operator is inert: max|K(op) - K(op:=0)| = {dmax}"
        );
        // (b) Exactly symmetric.
        for i in 0..24 {
            for j in 0..24 {
                assert_eq!(
                    diff[(i, j)],
                    diff[(j, i)],
                    "drill difference not exactly symmetric"
                );
            }
        }
        // (c) Every translational row/column block is exactly zero.
        for i in 0..4 {
            for j in 0..4 {
                for a in 0..3 {
                    for b in 0..3 {
                        assert_eq!(diff[(6 * i + a, 6 * j + b)], 0.0, "t-t block");
                        assert_eq!(diff[(6 * i + a, 6 * j + 3 + b)], 0.0, "t-r block");
                        assert_eq!(diff[(6 * i + 3 + a, 6 * j + b)], 0.0, "r-t block");
                    }
                }
            }
        }
        // (d) On flat geometry V^D = e3, so every rotation block other than the
        // drill's own is exactly zero.
        let flat = pre_from(&RECT);
        let df =
            compute_ke_local_with_drill(&flat, true) - compute_ke_local_with_drill(&flat, false);
        for i in 0..4 {
            for j in 0..4 {
                for a in 0..2 {
                    for b in 0..2 {
                        assert_eq!(df[(6 * i + 3 + a, 6 * j + 3 + b)], 0.0, "flat θx/θy block");
                    }
                    assert_eq!(df[(6 * i + 3 + a, 6 * j + 5)], 0.0, "flat θx-θz block");
                    assert_eq!(df[(6 * i + 5, 6 * j + 3 + a)], 0.0, "flat θz-θx block");
                }
            }
        }
        // Non-vacuity: the drill's own flat block is non-zero.
        let mut nz = 0.0f64;
        for i in 0..4 {
            for j in 0..4 {
                nz = nz.max(df[(6 * i + 5, 6 * j + 5)].abs());
            }
        }
        assert!(
            nz > 1e-10 * max_abs(&k_full).max(1e-30),
            "flat drill block is inert"
        );
        // (e) The six rigid-body fields carry zero energy, with the drill DOF
        // present.
        let lmax = lambda_max(&k_full);
        for (r, u_rb) in rigid_body_fields(&warped).iter().enumerate() {
            let e = quadratic(u_rb, &k_full).abs();
            let n2 = norm2(u_rb);
            assert!(
                e <= 1e-12 * lmax * n2,
                "rigid-body field {r}: |u^T K u| = {e}, bound {}",
                1e-12 * lmax * n2
            );
        }
    }

    // ------------------------------------------------------------------
    // 5.4 — uncorrected transverse shear
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_transverse_shear_uses_uncorrected_shear_modulus() {
        let h = 1.0;
        let e = 2.0e11;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let c = IsotropicMaterial::new(e, nu, 7800.0).constitutive(h, 5.0 / 6.0);
        let pre = Mitc4PlusDPrecomputed::new(&coords12(&RECT), c, h, 5.0 / 6.0);
        let k_shear = shear_ke_local(&pre);

        let closed_form = |cs: Matrix2<f64>| -> Mat24 {
            let mut k = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let sqrt_g = surface_measure(&pre, r, s);
                let bg = b_shear_mitc4(&pre, r, s);
                k += (bg.transpose() * cs * bg) * (GAUSS_W[gi] * sqrt_g);
            }
            k
        };
        // C_s = G·h·I, the uncorrected plane-stress shear stiffness.
        let k_ref = closed_form(Matrix2::new(g * h, 0.0, 0.0, g * h));
        let scale = max_abs(&k_ref).max(1e-30);
        let d = max_abs_diff(&k_shear, &k_ref);
        assert!(
            d <= 1e-10 * scale,
            "element shear block must be the uncorrected G·h closed form: diff {d}, scale {scale}"
        );
        // The 5/6 value is rejected by > 1e-3 relative.
        let k_wrong = closed_form(Matrix2::new(
            (5.0 / 6.0) * g * h,
            0.0,
            0.0,
            (5.0 / 6.0) * g * h,
        ));
        let rel = max_abs_diff(&k_wrong, &k_shear) / scale;
        assert!(
            rel > 1e-3,
            "the 5/6 value must be rejected (relative {rel})"
        );
    }

    #[test]
    fn test_identity_transverse_shear_invariant_to_shear_correction_factor() {
        let h = 1.0;
        let e = 2.0e11;
        let nu = 0.3;
        let g = e / (2.0 * (1.0 + nu));
        let coords = coords12(&RECT);

        // Single-ply laminate family: two `pre` values differing only in k.
        let lam_a = single_ply_lam(h, 5.0 / 6.0);
        let lam_b = single_ply_lam(h, 0.5);
        let pre_a = Mitc4PlusDPrecomputed::new(
            &coords,
            lam_a.to_shell_constitutive(),
            h,
            lam_a.applied_shear_correction_factor(),
        );
        let pre_b = Mitc4PlusDPrecomputed::new(
            &coords,
            lam_b.to_shell_constitutive(),
            h,
            lam_b.applied_shear_correction_factor(),
        );
        let ka = shear_ke_local(&pre_a);
        let kb = shear_ke_local(&pre_b);
        let scale = max_abs(&ka).max(1e-30);
        let d = max_abs_diff(&ka, &kb);
        assert!(
            d <= 1e-10 * scale,
            "single-ply shear block must be invariant to k: diff {d}, scale {scale}"
        );
        // Both equal the uncorrected closed form.
        let mut k_ref = Mat24::zeros();
        for gi in 0..4 {
            let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
            let sqrt_g = surface_measure(&pre_a, r, s);
            let bg = b_shear_mitc4(&pre_a, r, s);
            let cs = Matrix2::new(g * h, 0.0, 0.0, g * h);
            k_ref += (bg.transpose() * cs * bg) * (GAUSS_W[gi] * sqrt_g);
        }
        assert!(max_abs_diff(&ka, &k_ref) <= 1e-10 * max_abs(&k_ref).max(1e-30));
        // Non-vacuity control: a pre that carries the k factor in cs (applied_k
        // = 1.0) must differ by > 1e-3 relative.
        let pre_ctrl = Mitc4PlusDPrecomputed::new(&coords, lam_a.to_shell_constitutive(), h, 1.0);
        let k_ctrl = shear_ke_local(&pre_ctrl);
        let rel = max_abs_diff(&k_ctrl, &ka) / scale;
        assert!(
            rel > 1e-3,
            "the k-carrying control must be separated (relative {rel})"
        );

        // Isotropic constitutive: same construction and same controls.
        let iso_a = IsotropicMaterial::new(e, nu, 7800.0).constitutive(h, 5.0 / 6.0);
        let iso_b = IsotropicMaterial::new(e, nu, 7800.0).constitutive(h, 0.5);
        let ipre_a = Mitc4PlusDPrecomputed::new(&coords, iso_a, h, 5.0 / 6.0);
        let ipre_b = Mitc4PlusDPrecomputed::new(&coords, iso_b, h, 0.5);
        let ia = shear_ke_local(&ipre_a);
        let ib = shear_ke_local(&ipre_b);
        let iscale = max_abs(&ia).max(1e-30);
        assert!(
            max_abs_diff(&ia, &ib) <= 1e-10 * iscale,
            "isotropic k-invariance"
        );
        let ipre_ctrl = Mitc4PlusDPrecomputed::new(
            &coords,
            IsotropicMaterial::new(e, nu, 7800.0).constitutive(h, 5.0 / 6.0),
            h,
            1.0,
        );
        let ic = shear_ke_local(&ipre_ctrl);
        let irel = max_abs_diff(&ic, &ia) / iscale;
        assert!(
            irel > 1e-3,
            "isotropic k-carrying control must be separated (relative {irel})"
        );
    }

    // ------------------------------------------------------------------
    // 5.5 — the integration rule
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only() {
        let pre = pre_from(&STRONGLY_WARPED);
        let k_prod = compute_ke_local(&pre);
        let k_three = ke_ref::ke_local(&pre, true);
        let k_surface = ke_ref::ke_local(&pre, false);
        let scale = max_abs(&k_three).max(1e-30);
        let d = max_abs_diff(&k_prod, &k_three);
        assert!(
            d <= 1e-10 * scale,
            "production must match the three-term (2x2x2) reference: diff {d}, scale {scale}"
        );
        let rel = max_abs_diff(&k_three, &k_surface) / scale;
        assert!(
            rel > 1e-4,
            "the surface-only (W_22 = 0) reference must be separated (relative {rel})"
        );
    }

    // ------------------------------------------------------------------
    // 5.6 — local matrix shapes
    // ------------------------------------------------------------------
    #[test]
    fn test_kinematics_local_matrices_are_24x24() {
        let pre = pre_from(&DOUBLY_WARPED);
        let k_local = compute_ke_local(&pre);
        let k_global = compute_ke_global(&pre);
        assert_eq!(k_local.shape(), (24, 24), "K local must be 24x24");
        assert_eq!(k_global.shape(), (24, 24), "K global must be 24x24");
        assert_eq!(k_local.iter().count(), 576);
        assert_eq!(k_global.iter().count(), 576);
        // The DOF layout: 6 DOF/node with the drilling DOF at slot 6i+5. The
        // 2017-only operators are blind to that slot (they act through
        // `theta x V_n`), while the 2025 drill operator is live there.
        // The blindness of the 2017-only operators to slot 6i+5 is exactly
        // testable only on flat geometry, where the director V_n^i = e3 and
        // `theta x V_n` annihilates the local theta_z. (On warped geometry the
        // 2017 core does couple to the local theta_z, because e3 != V_n^i; the
        // drill rotation is theta.V^D, not theta_z -- see the WU4 findings.)
        let flat = pre_from(&RECT);
        for gi in 0..4 {
            let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
            let bm = b_membrane_2017(&flat, r, s);
            let (b1, b2) = b_bending_2017(&flat, r, s);
            let bg = b_shear_mitc4(&flat, r, s);
            for i in 0..4 {
                for row in 0..3 {
                    assert_eq!(bm[(row, 6 * i + 5)], 0.0, "membrane touches the drill slot");
                    assert_eq!(b1[(row, 6 * i + 5)], 0.0, "B_b1 touches the drill slot");
                    assert_eq!(b2[(row, 6 * i + 5)], 0.0, "B_b2 touches the drill slot");
                }
                assert_eq!(bg[(0, 6 * i + 5)], 0.0, "shear touches the drill slot");
                assert_eq!(bg[(1, 6 * i + 5)], 0.0, "shear touches the drill slot");
            }
        }
        // The drill operator is live at slot 6i+5 on flat geometry (V^D = e3):
        // the drill stiffness block has non-zero entries in that slot.
        let kd = drill_ke_local(&flat);
        let mut nz = 0.0f64;
        for i in 0..4 {
            nz = nz.max(kd[(6 * i + 5, 6 * i + 5)].abs());
        }
        assert!(nz > 0.0, "the drill slot 6i+5 is inert on flat geometry");
    }

    // ------------------------------------------------------------------
    // 5.7 — drill-DOF energy behaviour
    // ------------------------------------------------------------------
    #[test]
    fn test_kinematics_drill_dof_is_theta_z_through_eq26_operator() {
        // (a) A pure rigid rotation about V_n on a flat element: theta_i = γ e3
        // and u_i = γ e3 x x_i. No penalty, so the energy is zero.
        let flat = pre_from(&RECT);
        let kf = compute_ke_local(&flat);
        let lmax_f = lambda_max(&kf);
        let gamma = 0.7;
        let mut u_rb = [0.0f64; 24];
        for i in 0..4 {
            let x = flat.local_coords[i][0];
            let y = flat.local_coords[i][1];
            u_rb[6 * i] = -gamma * y;
            u_rb[6 * i + 1] = gamma * x;
            u_rb[6 * i + 5] = gamma;
        }
        let e_rb = quadratic(&u_rb, &kf).abs();
        assert!(
            e_rb <= 1e-12 * lmax_f * norm2(&u_rb),
            "rigid rotation about V_n carries energy {e_rb} (no penalty must couple to θ_z)"
        );

        // (b) A warped drill-rotation pattern theta_i = γ_i V_n^i violates the
        // constant-strain state. The 2017 core is blind to it (theta_i x V_n^i
        // = 0), so the energy is non-zero only through the Eq. (26) operator.
        let warped = pre_from(&DOUBLY_WARPED);
        let kd = compute_ke_local_with_drill(&warped, true);
        let kn = compute_ke_local_with_drill(&warped, false);
        let lmax_d = lambda_max(&kd);
        let mut u_z = [0.0f64; 24];
        for i in 0..4 {
            let gamma_i = if i % 2 == 0 { 1.0 } else { -1.0 };
            let vn = local_components(&warped, &warped.vn[i]);
            for k in 0..3 {
                u_z[6 * i + 3 + k] = gamma_i * vn[k];
            }
        }
        let n2 = norm2(&u_z);
        let e_without = quadratic(&u_z, &kn).abs();
        let e_with = quadratic(&u_z, &kd).abs();
        assert!(
            e_without <= 1e-12 * lmax_d * n2,
            "the 2017-only blocks must be blind to the drill pattern (energy {e_without})"
        );
        assert!(
            e_with > 1e-6 * lmax_d * n2,
            "the drill pattern must carry energy through Eq. (26) (energy {e_with})"
        );
    }

    // ------------------------------------------------------------------
    // 5.8 — mid-surface restriction (ADR-6 / G7)
    // ------------------------------------------------------------------
    #[test]
    fn test_identity_element_uses_midsurface_constitutive() {
        let h = 0.01;
        let lam = single_ply_lam(h, 5.0 / 6.0);
        let shell_mid = lam.to_shell_constitutive();
        let applied = lam.applied_shear_correction_factor();
        let pre =
            Mitc4PlusDPrecomputed::new(&coords12(&FLAT_SQUARE), shell_mid.clone(), h, applied);
        // The element's constitutive is exactly the mid-surface one.
        assert_eq!(pre.constitutive.cm, shell_mid.cm);
        assert_eq!(pre.constitutive.cb_coupling, shell_mid.cb_coupling);
        assert_eq!(pre.constitutive.cb, shell_mid.cb);
        assert_eq!(pre.constitutive.cs, shell_mid.cs);
        assert_eq!(pre.constitutive.cm_raw, shell_mid.cm_raw);
        assert_eq!(pre.applied_shear_correction, applied);
        // Non-vacuity: an offset reference surface would change the coupling
        // block to `B - z0 A` (ADR-6 / evidence gap G7), so the equality above
        // is only meaningful if the test can tell that block apart. The offset
        // block is recomputed here from its definition, so the element path
        // itself references no offset variant (the supporting static check).
        let z0 = 0.5 * h;
        let b_offset = lam.b - lam.a * z0;
        let mut diff = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                diff = diff.max((b_offset[(i, j)] - shell_mid.cb_coupling[(i, j)]).abs());
            }
        }
        assert!(
            diff > 0.0,
            "the offset coupling block must differ, otherwise the mid-surface assertion is vacuous"
        );
    }
}
