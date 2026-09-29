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

use nalgebra::{Matrix2, Matrix3, SMatrix, SVector, Vector2, Vector3, Vector4};

use crate::materials::ShellConstitutive;

/// 24×24 element matrix (4 nodes × 6 DOF/node).
pub type Mat24 = SMatrix<f64, 24, 24>;

/// 24-long element vector (4 nodes × 6 DOF/node), the layout the assembly
/// layers use for `f_int` and the element displacement.
pub type Vec24 = SVector<f64, 24>;

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
/// `e3` is the **area-weighted** normal of the two diagonal triangles
/// `(0,1,2)` and `(0,2,3)`, i.e. the normalized sum of their raw
/// (unnormalized) cross products. For a bilinear quad this is exactly the
/// paper's plane normal `n = (x_r x x_s)/||x_r x x_s||` of Ko, Lee & Bathe
/// (2017), C&S 182:404-418, Eq. (10), p. 406 — the identity
/// `x_r x x_s = (n1 + n2)/8` holds for every quad — and is therefore a
/// function of the geometry alone, independent of the node-numbering sequence.
///
/// `e1` is edge `0 -> 1` orthogonalized against `e3`; `e2 = e3 x e1`. `e1` is
/// not unique under a cyclic renumbering (it rotates in the `(e1, e2)` plane),
/// but the operators below are frame-covariant under an in-plane rotation about
/// `e3`, so the global stiffness does not depend on it.
///
/// WU6d finding: averaging the two **unit** triangle normals instead of the raw
/// ones gives a different `e3` on warped geometry with unequal triangle areas,
/// one that swaps with the diagonal under a cyclic node renumbering; the
/// covariant-to-local mapping is then no longer frame-covariant and `K_global`
/// acquires a spurious node-order dependence. See
/// `test_t1a_isotropy_element_orientation_and_node_sequence_invariant`.
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

    // Area-weighted normal of the two diagonal triangles (0,1,2) and (0,2,3):
    // sum the RAW cross products (not their unit normals), so that
    // `normalize(n1 + n2)` is the paper's Eq. (10) normal `normalize(x_r x x_s)`
    // for every geometry. Normalizing each triangle first would make `e3`
    // depend on which diagonal the node sequence selects (WU6d).
    let n1 = (nodes[1] - nodes[0]).cross(&(nodes[2] - nodes[0]));
    let n2 = (nodes[2] - nodes[0]).cross(&(nodes[3] - nodes[0]));
    let mut e3 = Vector3::zeros();
    let mut count = 0;
    if n1.norm() > 1e-12 {
        e3 += n1;
        count += 1;
    }
    if n2.norm() > 1e-12 {
        e3 += n2;
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
#[cfg(test)]
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
#[cfg(test)]
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
/// Eqs. (27a-c) is the closed form of Eqs. (17)+(18)+(19)+(21)/(25)+(26):
/// Appendix A Eq. (A.7), p. 416, is Eq. (21) with the printed `B_1..B_5`
/// substituted and retains the `e_rs^m|bil` terms, and substituting
/// Eqs. (25)/(26) into the assumed field reproduces (27a-c) term for term
/// (audit record: `openspec/changes/mitc4plusd-faithful/fidelity-audit.md`,
/// §B-audit record, 2026-09-24). The `(1 + a_E r s)` coefficient of Eq. (27c) is
/// the coefficient of the sampled `e_rs^m(E)`, not of `e_rs^m|bil`.
/// The engineering shear `2 e_rs` is applied by
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
/// W_00 = int C dz      = cm               (membrane)
/// W_01 = int z C dz    = cb_coupling      (membrane-bending coupling)
/// W_02 = int z^2 C dz  = cb               (membrane-E2 coupling)
/// W_11 = int z^2 C dz  = cb               (bending)
/// W_12 = int z^3 C dz  = 0                (odd moment)
/// W_22 = int z^4 C dz  = cm h^4 / 144     (E2)
/// ```
///
/// The `cm h^4 / 144` of `W_22` is the **4th moment** under the paper's own
/// 2x2 rule in the through-thickness coordinate `t` (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, p. 410: "2 x 2 x 2 Gauss integration over the element
/// domain"): two-point Gauss in `t = 2z/h` gives `int t^4 dt -> 2/9` where exact
/// integration gives `2/5`, so with `cm = int C dz = C h`,
/// `int z^4 C dz -> (h/2)(h/2)^4(2/9) C = h^5 C / 144 = cm h^4 / 144`.
///
/// The `(4/h^2)` of `E2` is carried by the B-operator, so `W_22` is the raw
/// 4th moment and **not** its `(4/h^2)^2`-scaled value: the effective `E2-E2`
/// coefficient the stiffness sees is
/// `(4/h^2)^2 W_22 = (16/h^4)(cm h^4/144) = cm/9`, the paper's
/// under-integrated value. Applying the `(4/h^2)^2` scaling to both `B_b2` and
/// `W_22` (the pre-WU9b form, `W_22 = cm/9`) multiplies the term by `16/h^4`,
/// which over-stiffens thin warped elements by orders of magnitude (the
/// MacNeal-Harder twisted beam, `tests/test_ko2017_performance.py::test_3_5`).
///
/// RECORDED APPROXIMATION (design Section 9, risk 7). For a homogeneous section
/// `cm h^4/144` is exact under the paper's rule. For a **multi-ply laminate** it
/// uses the section's thickness-average membrane stiffness (the same smearing
/// `cm_raw` already documents) in the `E2-E2` block only, instead of the true
/// `int z^4 C(z) dz`; the term it multiplies is second order in the warping, and
/// the surface-only discrimination of
/// `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only`
/// makes the approximation visible to a test rather than hidden.
pub fn resultant_moment_matrix(
    constitutive: &ShellConstitutive,
    thickness: f64,
) -> SMatrix<f64, 9, 9> {
    let cm = &constitutive.cm;
    let cb = &constitutive.cb;
    let cbc = &constitutive.cb_coupling;
    // `int z^4 C dz` under the paper's 2x2 through-thickness rule, whose value
    // is `cm h^4 / 144` -- see this function's docstring for the derivation and
    // for why the exact `h^5/80` is deliberately NOT used. The approximation that
    // remains is the multi-ply `cm`: the section's thickness-average membrane
    // stiffness stands in for `sum_k int_{z_k}^{z_k+1} z^4 C_k dz`, and the two
    // differ whenever the plies differ (design open item 6).
    let w22 = cm * (thickness.powi(4) / 144.0);
    let mut w = SMatrix::<f64, 9, 9>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            w[(i, j)] = cm[(i, j)];
            w[(i, 3 + j)] = cbc[(i, j)];
            w[(3 + i, j)] = cbc[(j, i)];
            w[(i, 6 + j)] = cb[(i, j)];
            w[(6 + i, j)] = cb[(j, i)];
            w[(3 + i, 3 + j)] = cb[(i, j)];
            w[(6 + i, 6 + j)] = w22[(i, j)];
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
#[cfg(test)]
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
#[cfg(test)]
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
    let w = resultant_moment_matrix(&pre.constitutive, pre.thickness);
    let s1 = 2.0 / pre.thickness;
    let s2 = 4.0 / (pre.thickness * pre.thickness);
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        // Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (22a), p. 10:
        //   e_ij = e_ij^m + e_ij^md + t e_ij^b1 + t^2 e_ij^b2.
        // The drill-membrane strain is ADDED to the assumed membrane strain, so
        // the t^0 row is `B_m + B_md` and the single `W` block supplies the
        // paper's cross terms. Adding the drill as a separate energy block would
        // drop them (and make the drill inert at the solution).
        let mut bm = b_membrane_2017(pre, r, s);
        if use_drill {
            bm += b_drill_membrane_2025(pre, r, s);
        }
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
/// K = sum_g [B_m + B_md; (2/h) B_b1; (4/h^2) B_b2]^T W [..] w sqrt_g
///   + sum_g B_gamma^T cs_uncorrected B_gamma w sqrt_g
/// ```
///
/// The drill-membrane strain is **added to the assumed membrane strain**, per
/// Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (22a), p. 10:
/// `e_ij = e_ij^m + e_ij^md + t e_ij^b1 + t^2 e_ij^b2`. The single `W` block
/// therefore carries the paper's cross terms as well: `(e^m)^T C e^md` at the
/// `t^0` level and `(e^m + e^md)^T C e^b2` at the `t^2` level (the `W`
/// `(m, b2)` entry `cb`). Adding the drill as an independent energy block would
/// drop both cross terms and would make the drill inert at the solution (a free
/// `theta_z` drives its own strain to zero).
///
/// with the moment matrix `W` of [`resultant_moment_matrix`], the 2x2 surface
/// rule (Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410; Ko, Bathe & Zhang
/// (2025), C&S 308:107622, §3.1, p. 13) and **no numerical factor**: the drill
/// term is the penalty-free 2025 drill-membrane strain (Eqs. (18)/(19d)/(22a))
/// and the shear block is the uncorrected `G h` of ADR-1.
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

// ============================================================================
// WU5 — the assembly-facing API
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

/// Shape-function derivatives with respect to the local orthonormal axes
/// `(e1, e2)` at `(r, s)`, i.e. `dh[beta, i] = dN_i/d(e_beta)`.
///
/// `j_loc[alpha, beta] = g_alpha . e_beta` (Ko, Lee & Bathe (2017),
/// C&S 182:404-418, Eq. (9)/(11), p. 406), so the chain rule is
/// `dN/d(e_beta) = sum_alpha (dN/dr_alpha) (j_loc^{-1})[alpha, beta]`.
fn local_shape_derivatives(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 2, 4> {
    let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
    let j_inv = regularized_inverse_2x2(&j_loc_at(pre, r, s));
    let mut dh = SMatrix::<f64, 2, 4>::zeros();
    for i in 0..4 {
        dh[(0, i)] = j_inv[(0, 0)] * dn_dr[i] + j_inv[(1, 0)] * dn_ds[i];
        dh[(1, i)] = j_inv[(0, 1)] * dn_dr[i] + j_inv[(1, 1)] * dn_ds[i];
    }
    dh
}

/// Displacement gradient `H = du/dX` (3x3, the third column zero) from the
/// local shape-function derivatives and the local 24-vector. Used only by the
/// superseded bounded path, hence `#[cfg(test)]`.
#[cfg(test)]
fn displacement_gradient(dh: &SMatrix<f64, 2, 4>, u: &Vec24) -> Matrix3<f64> {
    let mut u_nodes = [[0.0f64; 3]; 4];
    for i in 0..4 {
        u_nodes[i][0] = u[6 * i];
        u_nodes[i][1] = u[6 * i + 1];
        u_nodes[i][2] = u[6 * i + 2];
    }
    let mut h = Matrix3::zeros();
    for j in 0..2 {
        for comp in 0..3 {
            let mut val = 0.0;
            for nd in 0..4 {
                val += u_nodes[nd][comp] * dh[(j, nd)];
            }
            h[(comp, j)] = val;
        }
    }
    h
}

/// Green-Lagrange strain correction `1/2 (H^T H)` in Voigt form
/// `[E_xx, E_yy, 2 E_xy]` (the through-thickness terms vanish at `t = 0`).
///
/// This is the repository's total-Lagrangian covariant correction applied to
/// the new operators; it is the bounded nonlinear extension of the design's
/// open item 5 (risk 9). Neither Ko, Lee & Bathe (2017), C&S 182:404-418 nor
/// Ko, Bathe & Zhang (2025), C&S 308:107622 prints a nonlinear MITC4+/D
/// formulation for this repository's updated-Lagrangian form, so the linear
/// path is paper-faithful and this correction is bounded. It is superseded by
/// the faithful total-Lagrangian pair and is used only by the `#[cfg(test)]`
/// historical reference, hence `#[cfg(test)]`.
#[cfg(test)]
fn membrane_strain_nl(h: &Matrix3<f64>) -> Vector3<f64> {
    Vector3::new(
        0.5 * (h[(0, 0)].powi(2) + h[(1, 0)].powi(2) + h[(2, 0)].powi(2)),
        0.5 * (h[(0, 1)].powi(2) + h[(1, 1)].powi(2) + h[(2, 1)].powi(2)),
        h[(0, 0)] * h[(0, 1)] + h[(1, 0)] * h[(1, 1)] + h[(2, 0)] * h[(2, 1)],
    )
}

/// Derivative of [`membrane_strain_nl`] with respect to the local 24-vector,
/// `B_nl` (6x24; rows 0, 1, 3 carry `[E_xx, E_yy, 2 E_xy]`). Used only by the
/// superseded bounded path, hence `#[cfg(test)]`.
#[cfg(test)]
fn compute_b_nl(dh: &SMatrix<f64, 2, 4>, h: &Matrix3<f64>) -> SMatrix<f64, 6, 24> {
    let mut bnl = SMatrix::<f64, 6, 24>::zeros();
    for i in 0..4 {
        let col = 6 * i;
        let dni_dx = dh[(0, i)];
        let dni_dy = dh[(1, i)];

        bnl[(0, col)] = h[(0, 0)] * dni_dx;
        bnl[(0, col + 1)] = h[(1, 0)] * dni_dx;
        bnl[(0, col + 2)] = h[(2, 0)] * dni_dx;

        bnl[(1, col)] = h[(0, 1)] * dni_dy;
        bnl[(1, col + 1)] = h[(1, 1)] * dni_dy;
        bnl[(1, col + 2)] = h[(2, 1)] * dni_dy;

        bnl[(3, col)] = h[(0, 0)] * dni_dy + h[(0, 1)] * dni_dx;
        bnl[(3, col + 1)] = h[(1, 0)] * dni_dy + h[(1, 1)] * dni_dx;
        bnl[(3, col + 2)] = h[(2, 0)] * dni_dy + h[(2, 1)] * dni_dx;
    }
    bnl
}

/// Extract the membrane rows `[0, 1, 3]` from a 6x24 matrix into a 3x24. Used
/// only by the superseded bounded path, hence `#[cfg(test)]`.
#[cfg(test)]
fn extract_membrane_rows(b6: &SMatrix<f64, 6, 24>) -> SMatrix<f64, 3, 24> {
    let mut b3 = SMatrix::<f64, 3, 24>::zeros();
    for j in 0..24 {
        b3[(0, j)] = b6[(0, j)];
        b3[(1, j)] = b6[(1, j)];
        b3[(2, j)] = b6[(3, j)];
    }
    b3
}

/// Geometric B-matrix `B_geo` (6x24) for the initial-stress stiffness: rows
/// `[du/dx, du/dy, dv/dx, dv/dy, dw/dx, dw/dy]` per node.
fn compute_b_geometric(dh: &SMatrix<f64, 2, 4>) -> SMatrix<f64, 6, 24> {
    let mut bg = SMatrix::<f64, 6, 24>::zeros();
    for i in 0..4 {
        let col = 6 * i;
        let dni_dx = dh[(0, i)];
        let dni_dy = dh[(1, i)];

        bg[(0, col)] = dni_dx;
        bg[(1, col)] = dni_dy;
        bg[(2, col + 1)] = dni_dx;
        bg[(3, col + 1)] = dni_dy;
        bg[(4, col + 2)] = dni_dx;
        bg[(5, col + 2)] = dni_dy;
    }
    bg
}

/// Initial-stress (geometric) stiffness contribution at Gauss point `g` for the
/// membrane resultant state `sigma = [N_xx, N_yy, N_xy]`:
/// `B_geo^T blockdiag(sigma) B_geo w sqrt_g`.
fn geometric_stiffness_contribution(
    pre: &Mitc4PlusDPrecomputed,
    g: usize,
    sigma: &Vector3<f64>,
) -> Mat24 {
    let s_m = Matrix2::new(sigma[0], sigma[2], sigma[2], sigma[1]);
    let mut s_tilde = SMatrix::<f64, 6, 6>::zeros();
    for i in 0..2 {
        for j in 0..2 {
            s_tilde[(i, j)] = s_m[(i, j)];
            s_tilde[(i + 2, j + 2)] = s_m[(i, j)];
            s_tilde[(i + 4, j + 4)] = s_m[(i, j)];
        }
    }
    let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
    let sqrt_g = surface_measure(pre, r, s);
    let bg = compute_b_geometric(&local_shape_derivatives(pre, r, s));
    (bg.transpose() * s_tilde * bg) * (GAUSS_W[g] * sqrt_g)
}

/// Initial-stress stiffness in LOCAL coordinates from a pre-computed membrane
/// resultant state `sigma = [N_xx, N_yy, N_xy]`.
fn geometric_stiffness_from_stress(pre: &Mitc4PlusDPrecomputed, sigma: &Vector3<f64>) -> Mat24 {
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        k += geometric_stiffness_contribution(pre, g, sigma);
    }
    0.5 * (k + k.transpose())
}

/// SUPERSEDED bounded nonlinear membrane correction — retained ONLY as the
/// historical "before" reference for the `#[ignore]`d N-alpha and N-gamma
/// instruments.
///
/// ```text
/// f_corr = sum_g [ B_total^T cm eps - B_m^T cm (B_m u) ] w sqrt_g
/// B_total = B_m + B_nl,   eps = B_m u + 1/2 (H^T H)|Voigt
/// ```
///
/// so that `f_int(nonlinear) = K u + f_corr` and the linear part is exactly the
/// paper-faithful stiffness. `f_corr` is `O(||u||^2)` and vanishes at `u = 0`.
///
/// BOUNDED NONLINEAR PATH (design open item 5, risk 9): neither
/// Ko, Lee & Bathe (2017), C&S 182:404-418 nor Ko, Bathe & Zhang (2025),
/// C&S 308:107622 provides a nonlinear MITC4+/D formulation for this
/// repository's updated-Lagrangian form, so only this total-Lagrangian covariant
/// correction is carried. It was superseded by the faithful total-Lagrangian
/// pair `n_gamma_fint_global` / `n_gamma_kt_global` (N-delta, option B); this
/// body is `#[cfg(test)]`-only and is called solely by the two `#[ignore]`d
/// instruments as the recorded "before".
#[cfg(test)]
fn membrane_nonlinear_correction_bounded_superseded(
    pre: &Mitc4PlusDPrecomputed,
    u_local: &Vec24,
) -> Vec24 {
    let cm = &pre.constitutive.cm;
    let mut f = Vec24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let dh = local_shape_derivatives(pre, r, s);
        let h_mat = displacement_gradient(&dh, u_local);
        let bm = b_membrane_2017(pre, r, s);
        let bm_nl = extract_membrane_rows(&compute_b_nl(&dh, &h_mat));
        let mut b_total = bm;
        for i in 0..3 {
            for j in 0..24 {
                b_total[(i, j)] += bm_nl[(i, j)];
            }
        }
        let eps_lin = bm * u_local;
        let mut eps = eps_lin;
        let eps_nl = membrane_strain_nl(&h_mat);
        for i in 0..3 {
            eps[i] += eps_nl[i];
        }
        f += (b_total.transpose() * (cm * eps) - bm.transpose() * (cm * eps_lin))
            * (GAUSS_W[g] * sqrt_g);
    }
    f
}

/// SUPERSEDED bounded internal force — the historical "before" reference for
/// the `#[ignore]`d N-alpha and N-gamma instruments. Not the production path:
/// the production nonlinear path is [`n_gamma_fint_global`].
#[cfg(test)]
fn compute_fint_global_bounded_superseded(
    pre: &Mitc4PlusDPrecomputed,
    u_global: &Vec24,
) -> Vec24 {
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;
    let f_local = compute_ke_local(pre) * u_local
        + membrane_nonlinear_correction_bounded_superseded(pre, &u_local);
    t24.transpose() * f_local
}

/// Compute the internal force vector in GLOBAL coordinates (24 long).
///
/// `nonlinear = false` returns `K u` exactly (the linear path, unchanged).
///
/// `nonlinear = true` delegates to the faithful total-Lagrangian formulation
/// [`n_gamma_fint_global`] (Ko, Lee & Bathe (2017), C&S 185:1-14, Eq. (24b) on
/// the §2.4 incremental strains of Eqs. (19)-(23)), which is now the production
/// nonlinear behaviour of this element (N-delta, option B). The superseded
/// bounded correction survives only as the `#[cfg(test)]` historical reference
/// [`compute_fint_global_bounded_superseded`] for the `#[ignore]`d
/// instruments.
pub fn compute_fint_global(
    pre: &Mitc4PlusDPrecomputed,
    u_global: &Vec24,
    nonlinear: bool,
) -> Vec24 {
    if nonlinear {
        return n_gamma_fint_global(pre, u_global);
    }
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;
    t24.transpose() * (compute_ke_local(pre) * u_local)
}

/// Compute the tangent stiffness `K_T` in GLOBAL coordinates (24x24).
///
/// Delegates to the faithful total-Lagrangian tangent [`n_gamma_kt_global`]
/// (Ko, Lee & Bathe (2017), C&S 185:1-14, Eq. (24a)). That tangent is the
/// derivative of [`compute_fint_global`]`(.., true)` to round-off and satisfies
/// `K_T(0) = K_0` bit for bit, so the linear and nonlinear paths share one
/// operator.
pub fn compute_kt_global(pre: &Mitc4PlusDPrecomputed, u_global: &Vec24) -> Mat24 {
    n_gamma_kt_global(pre, u_global)
}

// ============================================================================
// N-alpha — faithful incremental (total-Lagrangian) kinematics
// ============================================================================
//
// Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear
// analysis", Computers and Structures 185:1-14, §2.1-2.2, pp. 2-3.
//
// This block adds the paper's incremental kinematic building blocks as
// self-contained functions over a plain-data current state. They are now the
// production nonlinear path: the public `compute_fint_global(.., true)` and
// `compute_kt_global` delegate to the N-gamma pair built on them (N-delta,
// option B), so the faithful formulation is what every solver calls.
//
// Notation (`^t` = the current configuration at time `t`, `^0` = the initial
// configuration, `zeta` = the thickness coordinate, `r_1 = r`, `r_2 = s`,
// `r_3 = zeta`):
//
//   Eq. (4c) `^{t+dt}V_n^i - ^t V_n^i = theta_i x ^t V_n^i
//            + 1/2 theta_i x (theta_i x ^t V_n^i)`,
//            `theta_i = ^t V_1^i alpha_i + ^t V_2^i beta_i`.
//   Eq. (5b) `u_m = sum_{i=1..4} h_i u_i`.
//   Eq. (5c) `u_b1 = 1/2 sum_i a_i h_i (-^t V_2^i alpha_i + ^t V_1^i beta_i)`,
//            `u_b2 = -1/4 sum_i a_i h_i (alpha_i^2 + beta_i^2) ^t V_n^i`.
//   Eq. (6)  `u_1 = u_m + zeta u_b1`, `u_2 = zeta u_b2`.
//   Eq. (7)  `^t_0 e_ij = 1/2(^t g_i . ^t g_j - ^0 g_i . ^0 g_j)` with
//            `^t g_i = d ^t x / d r_i`.
//   Eq. (8)  `_0 e_ij = 1/2(^t g_i . u_j + u_i . ^t g_j + u_i . u_j)` with
//            `u_i = d u / d r_i`.
//   Eq. (9)  `_0 e_ij = _0 e_ij + _0 eta_ij`, with
//            `_0 e_ij  = 1/2(^t g_i . u_{1,j} + u_{1,i} . ^t g_j)` (linear) and
//            `_0 eta_ij = 1/2(u_{1,i} . u_{1,j} + ^t g_i . u_{2,j}
//                        + u_{2,i} . ^t g_j)` (nonlinear, three terms).

/// The current configuration (`^t`, the state at time `t`) the faithful
/// incremental kinematics of Ko, Lee & Bathe (2017), C&S 185:1-14, Eqs. (1),
/// (4c), (7)-(9), needs.
///
/// It is plain data, decoupled from [`Mitc4PlusDPrecomputed`], so the same
/// functions can later be driven by a re-precomputed (updated) state or by the
/// frozen initial state. Eq. (1) makes the directors `^t V_n^i` current
/// quantities and Eq. (7) makes the covariant base vectors `^t g_i` current.
#[derive(Clone, Copy, Debug)]
pub struct GlCurrentState {
    /// Current node coordinates `^t x_i` (global).
    pub coords: [[f64; 3]; 4],
    /// Current nodal directors `^t V_n^i`, Eq. (1).
    pub vn: [Vector3<f64>; 4],
    /// Current in-plane directors `^t V_1^i` of Eq. (4c).
    pub v1: [Vector3<f64>; 4],
    /// Current in-plane directors `^t V_2^i = ^t V_n^i x ^t V_1^i` of Eq. (4c).
    pub v2: [Vector3<f64>; 4],
    /// Per-node thickness `a_i` of Eq. (1).
    pub a_i: [f64; 4],
}

/// The incremental unknowns of one load step: the nodal displacement
/// increments `u_i` of Eq. (4b) and the two director-rotation components
/// `alpha_i`, `beta_i` of Eqs. (4c)/(5c).
///
/// `alpha_i` and `beta_i` are the rotations of `^t V_n^i` about `^t V_1^i` and
/// `^t V_2^i`, i.e. `alpha_i = theta_i . ^t V_1^i` and
/// `beta_i = theta_i . ^t V_2^i` for a nodal rotation vector `theta_i`.
#[derive(Clone, Copy, Debug)]
pub struct GlIncrement {
    /// Nodal displacement increments `u_i` (global), Eq. (4b).
    pub u: [[f64; 3]; 4],
    /// Director rotation about `^t V_1^i`, Eq. (4c).
    pub alpha: [f64; 4],
    /// Director rotation about `^t V_2^i`, Eq. (4c).
    pub beta: [f64; 4],
}

/// The director increment `^{t+dt}V_n - ^t V_n` of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eq. (4c), to quadratic order:
///
/// ```text
/// theta x V_n + 1/2 theta x (theta x V_n)
/// ```
///
/// with `theta = ^t V_1 alpha + ^t V_2 beta` the nodal rotation vector. Calling
/// it with `theta` already formed from the `(alpha, beta)` pair is Eq. (4c)
/// verbatim.
pub fn director_increment_quadratic(vn: &Vector3<f64>, theta: &Vector3<f64>) -> Vector3<f64> {
    theta.cross(vn) + 0.5 * theta.cross(&theta.cross(vn))
}

/// The per-node director contributions of Eqs. (5c) at `h_i = 1`:
///
/// ```text
/// t1_i = 1/2 a_i (-^t V_2^i alpha_i + ^t V_1^i beta_i)
/// t2_i = -1/4 a_i (alpha_i^2 + beta_i^2) ^t V_n^i
/// ```
///
/// so that `u_b1 = sum_i h_i t1_i` and `u_b2 = sum_i h_i t2_i`. Shared by
/// [`incremental_disp_split`] and [`incremental_disp_gradients`].
#[inline]
fn node_director_terms(
    state: &GlCurrentState,
    inc: &GlIncrement,
    i: usize,
) -> (Vector3<f64>, Vector3<f64>) {
    let alpha = inc.alpha[i];
    let beta = inc.beta[i];
    let a2 = alpha * alpha + beta * beta;
    (
        0.5 * state.a_i[i] * (-state.v2[i] * alpha + state.v1[i] * beta),
        -0.25 * state.a_i[i] * a2 * state.vn[i],
    )
}

/// The Eq. (6) split `u = u_1 + u_2` at `(r, s, zeta)`, built from Eq. (5b)
/// (`u_m`) and Eq. (5c) (`u_b1`, `u_b2`, the **quadratic** director term):
///
/// ```text
/// u_1 = u_m + zeta u_b1,   u_2 = zeta u_b2
/// ```
///
/// Returns `(u_1, u_2)`.
pub fn incremental_disp_split(
    state: &GlCurrentState,
    inc: &GlIncrement,
    r: f64,
    s: f64,
    zeta: f64,
) -> (Vector3<f64>, Vector3<f64>) {
    let h = shape_functions(r, s);
    let mut u_m = Vector3::zeros();
    let mut u_b1 = Vector3::zeros();
    let mut u_b2 = Vector3::zeros();
    for i in 0..4 {
        let ui = Vector3::new(inc.u[i][0], inc.u[i][1], inc.u[i][2]);
        let (t1, t2) = node_director_terms(state, inc, i);
        u_m += h[i] * ui;
        u_b1 += h[i] * t1;
        u_b2 += h[i] * t2;
    }
    (u_m + zeta * u_b1, zeta * u_b2)
}

/// Derivatives of the Eq. (6) split with respect to `(r, s, zeta)`, i.e. the
/// `u_{1,i}` and `u_{2,i}` of Eq. (9).
///
/// `u_1 = u_m(r,s) + zeta u_b1(r,s)` and `u_2 = zeta u_b2(r,s)`, so the `r`/`s`
/// derivatives act on the shape functions and the `zeta` derivative is `u_b1` /
/// `u_b2`. Returns `(u_{1,r}, u_{1,s}, u_{1,t})` and `(u_{2,r}, u_{2,s}, u_{2,t})`
/// in that order (`r_3 = zeta`).
fn incremental_disp_gradients(
    state: &GlCurrentState,
    inc: &GlIncrement,
    r: f64,
    s: f64,
    zeta: f64,
) -> ([Vector3<f64>; 3], [Vector3<f64>; 3]) {
    let h = shape_functions(r, s);
    let (dn_dr, dn_ds) = shape_function_derivatives(r, s);
    let mut u_m = Vector3::zeros();
    let mut u_b1 = Vector3::zeros();
    let mut u_b2 = Vector3::zeros();
    let mut u_m_r = Vector3::zeros();
    let mut u_m_s = Vector3::zeros();
    let mut u_b1_r = Vector3::zeros();
    let mut u_b1_s = Vector3::zeros();
    let mut u_b2_r = Vector3::zeros();
    let mut u_b2_s = Vector3::zeros();
    for i in 0..4 {
        let ui = Vector3::new(inc.u[i][0], inc.u[i][1], inc.u[i][2]);
        let (t1, t2) = node_director_terms(state, inc, i);
        u_m += h[i] * ui;
        u_b1 += h[i] * t1;
        u_b2 += h[i] * t2;
        u_m_r += dn_dr[i] * ui;
        u_m_s += dn_ds[i] * ui;
        u_b1_r += dn_dr[i] * t1;
        u_b1_s += dn_ds[i] * t1;
        u_b2_r += dn_dr[i] * t2;
        u_b2_s += dn_ds[i] * t2;
    }
    (
        [
            u_m_r + zeta * u_b1_r,
            u_m_s + zeta * u_b1_s,
            u_b1,
        ],
        [zeta * u_b2_r, zeta * u_b2_s, u_b2],
    )
}

/// The current covariant base vectors `[^t g_r, ^t g_s, ^t g_zeta]` of
/// Eq. (7) at `(r, s, zeta)`, formed from the current node coordinates
/// `^t x_i` and the current directors `^t V_n^i`:
///
/// ```text
/// ^t x(r,s,zeta) = sum_i h_i ^t x_i + (zeta/2) sum_i a_i h_i ^t V_n^i
/// ^t g_a = d ^t x / d r_a,  ^t g_zeta = d ^t x / d zeta
/// ```
///
/// which is the enriched base-vector construction of [`compute_j3d_enriched`]
/// evaluated on the current state.
fn gl_current_base_vectors(
    state: &GlCurrentState,
    r: f64,
    s: f64,
    zeta: f64,
) -> [Vector3<f64>; 3] {
    let (g_r, g_s, g_t) =
        compute_j3d_enriched(&state.coords, &state.vn, &state.a_i, r, s, zeta);
    [g_r, g_s, g_t]
}

/// Assemble the six independent covariant components of Eq. (9) from the
/// current base vectors `^t g_i` and the split gradients `u_{1,i}`, `u_{2,i}`.
///
/// Component order is `[rr, ss, tt, rs, rt, st]` and the shear components are
/// the tensor components (not doubled). The linear and nonlinear parts are
/// returned separately:
///
/// ```text
/// _0 e_ij   = 1/2(^t g_i . u_{1,j} + u_{1,i} . ^t g_j)
/// _0 eta_ij = 1/2(u_{1,i} . u_{1,j} + ^t g_i . u_{2,j} + u_{2,i} . ^t g_j)
/// ```
fn gl_strain_increment_components(
    g: &[Vector3<f64>; 3],
    u1: &[Vector3<f64>; 3],
    u2: &[Vector3<f64>; 3],
) -> ([f64; 6], [f64; 6]) {
    // The six independent index pairs (r,s,zeta) = (0,1,2).
    const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
    let mut e_lin = [0.0f64; 6];
    let mut eta_nl = [0.0f64; 6];
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        e_lin[k] = 0.5 * (g[i].dot(&u1[j]) + u1[i].dot(&g[j]));
        eta_nl[k] = 0.5 * (u1[i].dot(&u1[j]) + g[i].dot(&u2[j]) + u2[i].dot(&g[j]));
    }
    (e_lin, eta_nl)
}

/// The incremental covariant Green-Lagrange strains of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eq. (9), at `(r, s, zeta)`: the linear part `_0 e_ij` and the
/// nonlinear part `_0 eta_ij` (with its three printed terms), using the CURRENT
/// base vectors `^t g_i` of Eq. (7).
///
/// Returns two six-component vectors in the order `[rr, ss, tt, rs, rt, st]`
/// with tensor (undoubled) shear components. A finite rigid-body rotation must
/// make `_0 e_ij + _0 eta_ij` vanish to the paper's retained order; the
/// companion test `n_alpha_rigid_rotation_gives_zero_gl_strain_increment`
/// measures that identity and the necessity of the quadratic `u_2` term.
pub fn gl_strain_increment(
    state: &GlCurrentState,
    inc: &GlIncrement,
    r: f64,
    s: f64,
    zeta: f64,
) -> ([f64; 6], [f64; 6]) {
    let g = gl_current_base_vectors(state, r, s, zeta);
    let (u1, u2) = incremental_disp_gradients(state, inc, r, s, zeta);
    gl_strain_increment_components(&g, &u1, &u2)
}

// ============================================================================
// N-beta — the paper's assumed strain fields on the total-Lagrangian strains
// ============================================================================
//
// Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear
// analysis", Computers and Structures 185:1-14, §2.3-2.4, pp. 3-4.
//
// This block adds the paper's assumed strain fields on the total-Lagrangian
// (Green-Lagrange) strains of the N-alpha block. They feed the N-gamma pair,
// which the public `compute_fint_global(.., true)` / `compute_kt_global` now
// delegate to (N-delta, option B). The `^0`
// (initial-configuration) input is the precomputed element `pre`; the `^t`
// (current-configuration) input is the plain-data [`GlCurrentState`] the
// N-alpha block already uses.
//
// The assumed fields:
//
//   Eq. (10)   §2.2, assumed covariant transverse shear,
//              `^t_0 e~_rzeta = 1/2(1+s) ^t_0 e_rzeta^(A)
//                              + 1/2(1-s) ^t_0 e_rzeta^(B)`,
//              `^t_0 e~_szeta = 1/2(1+r) ^t_0 e_szeta^(C)
//                              + 1/2(1-r) ^t_0 e_szeta^(D)`.
//   Eq. (11a)  the through-thickness split `^t_0 e_ij = ^t_0 e_ij^m
//              + zeta ^t_0 e_ij^b1 + zeta^2 ^t_0 e_ij^b2` (this block supplies
//              the `m` term only; `b1`/`b2` are N-gamma).
//   Eqs. (12b)/(12c): `^t g_ij^m = ^t x_{m,i} . ^t x_{m,j}
//              + ^t x_{m,j} . ^t x_{m,i}`.
//   Eq. (13)   `^t x_{m,r} = ^t x_r + s ^t x_d`,
//              `^t x_{m,s} = ^t x_s + r ^t x_d`.
//   Eq. (14)   `^t n = (^t x_r x ^t x_s)/||^t x_r x ^t x_s||`, dual basis
//              `^t m^r`, `^t m^s` on the plane.
//   Eq. (15a)  `^t_0 e~_ij^m = 1/2 ^t g~_ij^m - 1/2 ^0 g~_ij^m`,
//   Eqs. (15b)/(15c)/(15d) give `^t g~_ij^m` from the five tying metrics with
//              the CURRENT coefficients `^t a_A..^t a_E` (same formulas as the
//              linear element, evaluated on the current configuration).

/// The current characteristic vectors and dual basis of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eqs. (13)/(14) / C&S 182:404-418, Eq. (9)-(11):
///
/// ```text
/// ^t x_r = 1/4 sum_i xi_i  ^t x_i
/// ^t x_s = 1/4 sum_i eta_i ^t x_i
/// ^t x_d = 1/4 sum_i xi_i eta_i ^t x_i
/// ^t n = (^t x_r x ^t x_s)/||^t x_r x ^t x_s||
/// ^t m^r, ^t m^s : ^t m^{r_i} . ^t x_{r_j} = delta_ij, ^t m^{r_i} . ^t n = 0
/// ```
///
/// Returns `([^t x_r, ^t x_s, ^t x_d], ^t n, ^t m^r, ^t m^s)`.
pub fn gl_current_characteristic_vectors(
    state: &GlCurrentState,
) -> (
    [Vector3<f64>; 3],
    Vector3<f64>,
    Vector3<f64>,
    Vector3<f64>,
) {
    let (x_r, x_s, x_d, n_vec, m_r, m_s) = compute_characteristic_vectors(&state.coords);
    ([x_r, x_s, x_d], n_vec, m_r, m_s)
}

/// The mid-surface metric components of Ko, Lee & Bathe (2017), C&S 185:1-14,
/// Eqs. (12b)/(12c), at `(r, s)`, built from the Eq. (13) mid-surface tangents
/// `^t x_{m,r} = ^t x_r + s ^t x_d` and `^t x_{m,s} = ^t x_s + r ^t x_d`:
///
/// ```text
/// ^t g_rr^m = 2 ^t x_{m,r} . ^t x_{m,r}
/// ^t g_ss^m = 2 ^t x_{m,s} . ^t x_{m,s}
/// ^t g_rs^m = ^t x_{m,r} . ^t x_{m,s} + ^t x_{m,s} . ^t x_{m,r}
/// ```
///
/// Returns `[^t g_rr^m, ^t g_ss^m, ^t g_rs^m]`. The same formula gives
/// `^0 g_ij^m` when fed the initial characteristic vectors.
pub fn gl_mid_metrics(
    x_r: &Vector3<f64>,
    x_s: &Vector3<f64>,
    x_d: &Vector3<f64>,
    r: f64,
    s: f64,
) -> [f64; 3] {
    let x_mr = x_r + s * x_d;
    let x_ms = x_s + r * x_d;
    [
        2.0 * x_mr.dot(&x_mr),
        2.0 * x_ms.dot(&x_ms),
        x_mr.dot(&x_ms) + x_ms.dot(&x_mr),
    ]
}

/// The five tying metrics of the MITC4+ assumed membrane field,
/// Ko, Lee & Bathe (2017), C&S 185:1-14, Eqs. (15b)/(15c)/(15d): evaluated at
/// the same five tying points as the linear element (Fig. 4 of
/// C&S 182:404-418), A(0,+1) for `g_rr^m(A)`, B(0,-1) for `g_rr^m(B)`,
/// C(+1,0) for `g_ss^m(C)`, D(-1,0) for `g_ss^m(D)`, E(0,0) for `g_rs^m(E)`.
///
/// Returns `[g_rr^m(A), g_rr^m(B), g_ss^m(C), g_ss^m(D), g_rs^m(E)]`.
pub fn gl_tying_metrics(x_r: &Vector3<f64>, x_s: &Vector3<f64>, x_d: &Vector3<f64>) -> [f64; 5] {
    let a = gl_mid_metrics(x_r, x_s, x_d, 0.0, 1.0)[0];
    let b = gl_mid_metrics(x_r, x_s, x_d, 0.0, -1.0)[0];
    let c = gl_mid_metrics(x_r, x_s, x_d, 1.0, 0.0)[1];
    let d = gl_mid_metrics(x_r, x_s, x_d, -1.0, 0.0)[1];
    let e = gl_mid_metrics(x_r, x_s, x_d, 0.0, 0.0)[2];
    [a, b, c, d, e]
}

/// The assumed mid-surface metric `[^t g~_rr^m, ^t g~_ss^m, ^t g~_rs^m]` of
/// Ko, Lee & Bathe (2017), C&S 185:1-14, Eqs. (15b)/(15c)/(15d), at `(r, s)`.
///
/// `tie = [g_rr^m(A), g_rr^m(B), g_ss^m(C), g_ss^m(D), g_rs^m(E)]` are the five
/// tying metrics of [`gl_tying_metrics`]; `coeffs = [a_A, a_B, a_C, a_D, a_E]`
/// are the [`compute_membrane_coefficients_2017`] coefficients evaluated on the
/// configuration the tying metrics come from. The coefficient structure is
/// exactly the linear element's (Eqs. 27a-c), only the coefficients and the
/// tying metrics are current quantities.
pub fn gl_assumed_mid_metric(coeffs: &[f64; 5], tie: &[f64; 5], r: f64, s: f64) -> [f64; 3] {
    let a = coeffs;
    let g = tie;
    [
        0.5 * (1.0 - 2.0 * a[0] + s + 2.0 * a[0] * s * s) * g[0]
            + 0.5 * (1.0 - 2.0 * a[1] - s + 2.0 * a[1] * s * s) * g[1]
            + a[2] * (-1.0 + s * s) * g[2]
            + a[3] * (-1.0 + s * s) * g[3]
            + a[4] * (-1.0 + s * s) * g[4],
        a[0] * (-1.0 + r * r) * g[0]
            + a[1] * (-1.0 + r * r) * g[1]
            + 0.5 * (1.0 - 2.0 * a[2] + r + 2.0 * a[2] * r * r) * g[2]
            + 0.5 * (1.0 - 2.0 * a[3] - r + 2.0 * a[3] * r * r) * g[3]
            + a[4] * (-1.0 + r * r) * g[4],
        0.25 * (r + 4.0 * a[0] * r * s) * g[0]
            + 0.25 * (-r + 4.0 * a[1] * r * s) * g[1]
            + 0.25 * (s + 4.0 * a[2] * r * s) * g[2]
            + 0.25 * (-s + 4.0 * a[3] * r * s) * g[3]
            + (1.0 + a[4] * r * s) * g[4],
    ]
}

/// The total-Lagrangian assumed membrane strain of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eq. (15a):
///
/// ```text
/// ^t_0 e~_ij^m = 1/2 ^t g~_ij^m - 1/2 ^0 g~_ij^m,   i, j = 1, 2
/// ```
///
/// CONVENTION. Eq. (12b) prints the assumed metric as the symmetric Gram sum
/// `^t g_ij^m = ^t x_{m,i} . ^t x_{m,j} + ^t x_{m,j} . ^t x_{m,i}`, i.e. TWICE
/// the metric `^t m_ij^m = ^t x_{m,i} . ^t x_{m,j}`. [`gl_tying_metrics`] and
/// [`gl_assumed_mid_metric`] therefore return the paper's doubled `g~`, so the
/// `1/2` of Eq. (15a) already converts it to the metric (`1/2 g~ = m~`). The
/// Green-Lagrange strain carried by the validated linear operator
/// [`b_membrane_covariant_2017`] is the further `1/2` of the metric difference,
///
/// ```text
/// ^t_0 e~_ij^m = 1/2(^t m~_ij^m - ^0 m~_ij^m) = 1/4(^t g~_ij^m - ^0 g~_ij^m),
/// ```
///
/// so the total conversion from the doubled tying metric is `1/4` -- exactly the
/// `0.25 * (next_tie - cur_tie)` that [`n_gamma_local_strain`] uses. An earlier
/// implementation applied the literal `1/2` of Eq. (15a) to the already-doubled
/// `g~`, returning twice the linear path's strain; the error cancelled against
/// the doubling of `gl_mid_metrics` in its own measurement, which is why that
/// check was vacuous.
///
/// The current assumed metric `^t g~_ij^m` is Eqs. (15b-d) with the current
/// tying metrics [`gl_tying_metrics`] of `state` and the current coefficients
/// [`compute_membrane_coefficients_2017`] on `state`; the initial
/// `^0 g~_ij^m` uses the precomputed `pre` geometry and `pre.a_coeffs`. Returns
/// `[^t_0 e~_rr^m, ^t_0 e~_ss^m, ^t_0 e~_rs^m]` (tensor, undoubled shear).
pub fn gl_assumed_membrane_strain(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    r: f64,
    s: f64,
) -> [f64; 3] {
    let cur_assumed = gl_assumed_metric_of_state(state, r, s);

    let init_tie = gl_tying_metrics(&pre.x_r, &pre.x_s, &pre.x_d);
    let init_assumed = gl_assumed_mid_metric(&pre.a_coeffs, &init_tie, r, s);

    [
        0.25 * cur_assumed[0] - 0.25 * init_assumed[0],
        0.25 * cur_assumed[1] - 0.25 * init_assumed[1],
        0.25 * cur_assumed[2] - 0.25 * init_assumed[2],
    ]
}

/// The next configuration `^{t+dt}` of one increment, using the paper's own
/// solution-procedure updates:
///
///   Eq. (4b)  `^{t+dt}x_i = ^t x_i + u_i` (nodal geometry),
///   Eq. (26)  `^{t+dt}V_n^i = Q ^t V_n^i`, `^{t+dt}V_1^i = Q ^t V_1^i`,
///             `^{t+dt}V_2^i = Q ^t V_2^i`, with `Q` the quaternion rotation of
///             `theta_i = alpha_i ^t V_1^i + beta_i ^t V_2^i`.
///
/// Eq. (26) prints the exact finite rotation (`q_0 = cos(theta_i/2)`,
/// `[q_1 q_2 q_3]^T = theta_i/theta_i sin(theta_i/2)`), and its use here is what
/// makes a finite rigid-body field leave the total metric — and therefore every
/// Green-Lagrange strain — invariant to round-off. Shared by
/// [`gl_assumed_membrane_increment`] and the N-gamma block.
pub fn advance_gl_state(state: &GlCurrentState, inc: &GlIncrement) -> GlCurrentState {
    let mut next = *state;
    for i in 0..4 {
        for k in 0..3 {
            next.coords[i][k] += inc.u[i][k];
        }
        let theta = inc.alpha[i] * state.v1[i] + inc.beta[i] * state.v2[i];
        let q = Mitc4PlusDPrecomputed::quaternion_from_vector(&theta);
        next.vn[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.vn[i], &q);
        next.v1[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.v1[i], &q);
        next.v2[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.v2[i], &q);
    }
    next
}

/// The assumed mid-surface metric `[^t g~_rr^m, ^t g~_ss^m, ^t g~_rs^m]` of one
/// configuration `^t` alone, Ko, Lee & Bathe (2017), C&S 185:1-14,
/// Eqs. (15b)/(15c)/(15d): the tying metrics [`gl_tying_metrics`] of `state`
/// combined with the coefficients [`compute_membrane_coefficients_2017`] of the
/// **same** configuration (Eqs. 15b-d with Eqs. 12b/12c). Returns the paper's
/// **doubled** `g~` (Eq. 12b), not the metric; [`gl_assumed_membrane_strain`]
/// applies the `1/4` conversion of Eq. (15a) plus the Green-Lagrange `1/2`.
fn gl_assumed_metric_of_state(state: &GlCurrentState, r: f64, s: f64) -> [f64; 3] {
    let (vecs, _n, m_r, m_s) = gl_current_characteristic_vectors(state);
    let (_c_r, _c_s, _d, coeff) = compute_membrane_coefficients_2017(&vecs[2], &m_r, &m_s);
    let tie = gl_tying_metrics(&vecs[0], &vecs[1], &vecs[2]);
    gl_assumed_mid_metric(&coeff, &tie, r, s)
}

/// The assumed **increment** of the membrane strain over one step,
/// `^{t+dt}_0 e~_ij^m - ^t_0 e~_ij^m`.
///
/// CONVENTION (the same as [`gl_assumed_membrane_strain`]). Eq. (12b) prints
/// `^t g_ij^m = ^t x_{m,i} . ^t x_{m,j} + ^t x_{m,j} . ^t x_{m,i}`, twice the
/// metric, and Eq. (15a) is `^t_0 e~ = 1/2 ^t g~ - 1/2 ^0 g~`; the `1/2` of
/// Eq. (15a) turns `g~` into the metric and the Green-Lagrange `1/2` of the
/// metric difference gives the strain, so the total conversion from the doubled
/// tying metric is `1/4`: `0.25 * (^{t+dt}tie - ^t tie)`.
///
/// Ko, Lee & Bathe (2017), C&S 185:1-14, §2.4 (Eqs. 21a-c and 22) states the
/// assumed INCREMENTAL in-plane field with the CURRENT coefficients
/// `^t a_A..^t a_E` (Eq. 15e), evaluated on `state` and **frozen**, applied to
/// the one-step increment of the tying metrics. That is the construction used
/// here. Differencing the two total assumed metrics `1/4(^{t+dt}g~ - ^t g~)`,
/// each with its own coefficients, would instead add the Eq. (15e)
/// coefficient-change term (which is `O(1)` in the initial metric), and the
/// `u = 0` operator would then not be the element's linear membrane operator;
/// the total form belongs to the state ([`gl_assumed_membrane_strain`]), not to
/// the tangent. This is the same frozen-coefficient form [`n_gamma_local_strain`]
/// uses, so the `u = 0` derivative of the returned field is
/// [`b_membrane_covariant_2017`].
///
/// Returns `[^{t+dt}e~_rr^m - ^t e~_rr^m, ...]`, tensor shear undoubled.
pub fn gl_assumed_membrane_increment(
    state: &GlCurrentState,
    inc: &GlIncrement,
    r: f64,
    s: f64,
) -> [f64; 3] {
    let next = advance_gl_state(state, inc);
    let (cur_vecs, _n, m_r, m_s) = gl_current_characteristic_vectors(state);
    let (_c_r, _c_s, _d, cur_coeff) =
        compute_membrane_coefficients_2017(&cur_vecs[2], &m_r, &m_s);
    let (next_vecs, _n_next, _m_r_next, _m_s_next) = gl_current_characteristic_vectors(&next);
    let cur_tie = gl_tying_metrics(&cur_vecs[0], &cur_vecs[1], &cur_vecs[2]);
    let next_tie = gl_tying_metrics(&next_vecs[0], &next_vecs[1], &next_vecs[2]);
    let mut dtie = [0.0f64; 5];
    for k in 0..5 {
        dtie[k] = 0.25 * (next_tie[k] - cur_tie[k]);
    }
    gl_assumed_mid_metric(&cur_coeff, &dtie, r, s)
}

/// The assumed covariant transverse shear of Ko, Lee & Bathe (2017),
/// C&S 185:1-14, Eq. (10):
///
/// ```text
/// ^t_0 e~_rzeta = 1/2(1+s) ^t_0 e_rzeta^(A) + 1/2(1-s) ^t_0 e_rzeta^(B)
/// ^t_0 e~_szeta = 1/2(1+r) ^t_0 e_szeta^(C) + 1/2(1-r) ^t_0 e_szeta^(D)
/// ```
///
/// The tying values are the Eq. (9) TL shear increments at the MITC4 tying
/// points A(0,+1), B(0,-1), C(+1,0), D(-1,0), i.e. the `rt = 4` and `st = 5`
/// components of [`gl_strain_increment`] (linear + nonlinear), sampled at
/// `zeta`. Returns `(^t_0 e~_rzeta, ^t_0 e~_szeta)`, tensor shear undoubled.
pub fn gl_assumed_transverse_shear(
    state: &GlCurrentState,
    inc: &GlIncrement,
    r: f64,
    s: f64,
    zeta: f64,
) -> (f64, f64) {
    let rt = |r_t: f64, s_t: f64| {
        let (e_lin, eta_nl) = gl_strain_increment(state, inc, r_t, s_t, zeta);
        e_lin[4] + eta_nl[4]
    };
    let st = |r_t: f64, s_t: f64| {
        let (e_lin, eta_nl) = gl_strain_increment(state, inc, r_t, s_t, zeta);
        e_lin[5] + eta_nl[5]
    };
    (
        0.5 * (1.0 + s) * rt(0.0, 1.0) + 0.5 * (1.0 - s) * rt(0.0, -1.0),
        0.5 * (1.0 + r) * st(1.0, 0.0) + 0.5 * (1.0 - r) * st(-1.0, 0.0),
    )
}

// ============================================================================
// N-gamma — faithful total-Lagrangian internal force and tangent stiffness
// ============================================================================
//
// Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear
// analysis", Computers and Structures 185:1-14, §2.4-2.5, pp. 6-7.
//
// What the paper prints (all C&S 185:1-14):
//
//   Eq. (8)/(9)   `_0 e_ij = 1/2(^t g_i . u_j + u_i . ^t g_j + u_i . u_j)`
//                 with `u = u_1 + u_2` of Eqs. (5b)/(5c)/(6) — the INCREMENTAL
//                 Green-Lagrange strain about the current configuration `^t`;
//   Eq. (20a-g)   its point-wise `zeta^0/zeta^1/zeta^2` coefficients,
//                 `_0 e_ij` (Eqs. 20b-d) and `_0 eta_ij` (Eqs. 20e-g);
//   Eq. (21a-c)   the assumed incremental in-plane field, applied by Eq. (22)
//                 to the incremental tying strains of Eqs. (20b)/(20e);
//   Eq. (19)      the assumed incremental transverse shear (Eq. 10 for `^t_0`);
//   Eq. (23)      `_0 e~_ij = _0 e~_kl (^0 L_i . ^0 g^k)(^0 L_j . ^0 g^l)`;
//   Eq. (24b)     `^t_0 F_e = int_{0V} B_ij^T ^t_0 S_ij d0V`  — internal force;
//   Eq. (24a)     `^t K_e = int_{0V} B_ij^T C_ijkl B_kl d0V
//                            + int_{0V} ^t_0 S_ij N_ij d0V`   — tangent stiffness;
//   Eq. (25)      `_0 e_ij = B_ij U_e`, `delta _0 eta_ij = delta U_e^T N_ij U_e`.
//
// The paper does not print `B_ij` or `N_ij` in closed form. This block therefore
// implements Eq. (24a)/(24b) with the consistent linearization of exactly the
// printed INCREMENTAL strain expressions of §2.4, evaluated about the current
// configuration `^t = state` supplied by the caller. §2.4's incremental object
// is not §2.3's total TL strain: Eqs. (15a-d) apply the assumed field to the
// metric of EACH configuration with that configuration's own coefficients
// `^t a_A..^t a_E` (Eq. 15e), so differencing them injects the coefficient-change
// term into `u = 0`; Eqs. (21)/(22) apply the assumed field to the incremental
// tying STRAINS with the CURRENT coefficients only, which is the object whose
// `u = 0` derivative is the element's linear strain operator.
//
// The strain vector is the nine-component through-thickness triple
// `[_0 e~^m + _0 eta~^m; (2/h)(_0 e~^b1 + _0 eta~^b1); (4/h^2)(_0 e~^b2 +
// _0 eta~^b2)]` of the linear element (same `W` block and the same `2 x 2`
// surface rule) plus the two local transverse shears. The internal force is
// `f = int B^T S dV`, and the tangent keeps Eq. (24a)'s two terms:
// `int B^T C B dV` (material) and the `N`-term `int (dB/du)^T S dV`
// (geometric), with `B = d(e~)/du` and `S = W e~`. Because `B` is not printed,
// it is obtained by central differences of the printed strain expressions; the
// same construction is used for the internal force, so the pair is consistent
// by construction and Eq. (24a) holds for the implemented strains.

/// The N-gamma local strain of Eqs. (20a)-(23), SPLIT into its linear and its
/// quadratic (nonlinear) part.
///
/// WHY THE SPLIT, AND WHY IT IS EXACT. The paper's Eq. (25) defines the tangent's
/// two operators by their property -- `_0 e~_ij = B_ij U_e` for the linear part,
/// `delta _0 eta~_ij = delta U_e^T N_ij U_e` for the nonlinear one -- and prints
/// neither matrix. Both are extracted exactly from this split by
/// [`n_gamma_b_linear`] and [`n_gamma_n_matrix`], with no finite difference and no
/// hand derivation.
///
/// The split is an algebraic identity, not a difference quotient. Every quantity
/// below is a difference of products of two vectors, `a . b - a0 . b0`, and with
/// `d = a - a0` that is exactly
///
/// ```text
/// a . b - a0 . b0 = (a0 . d + d . b0)  +  d . d ,
/// ```
///
/// a term linear in the increment plus a term quadratic in it. The same identity
/// covers the metric form `2 a . a - 2 a0 . a0 = (4 a0 . d) + (2 d . d)` and the
/// cross form `2 a . b - 2 a0 . b0 = 2(a0 . d + d . b0) + 2 d . d`, so the five
/// tying metrics of Eqs. (15b-d) split the same way. The assumed field of Eq. (21)
/// and the Eq. (23) map are linear, so they carry the split through unchanged, and
/// the drill row of Eq. (22a) is linear and contributes nothing to the quadratic
/// part.
///
/// [`n_gamma_local_strain`] is the elementwise sum of the two halves, which is what
/// makes the split checkable: the whole suite pins the total.
///
/// Because every slice is an increment about `state` and the assumed field's
/// coefficients and tying values are both taken from `state` and frozen, the
/// `u = 0` limit of the linear part is the element's LINEAR strain operator -- the
/// Eq. (15a) total form differenced one assumed total per configuration and its
/// `u = 0` derivative therefore carried the Eq. (15e) coefficient-change term
/// (which is `O(1)` in the initial metric, not `O(u)`).
///
/// A finite rigid-body field leaves every metric (and therefore every slice)
/// invariant to round-off, so the assembled force of [`n_gamma_fint_local`] is
/// exactly zero for translations and finite rotations.
fn n_gamma_local_strain_parts(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
    r: f64,
    s: f64,
) -> ([f64; 11], [f64; 11]) {
    let next = gl_state_after_local(pre, state, u_local);
    // Eq. (23): the covariant-to-local map uses the `^0` metric and `^0` frame.
    let tmap = covariant_to_local_mapping(&j_loc_at(pre, r, s));
    let s1 = 2.0 / pre.thickness;
    let s2 = 4.0 / (pre.thickness * pre.thickness);

    // Eqs. (21a-c) + (22): the assumed in-plane field on the INCREMENTAL tying
    // strains of Eq. (20b)/(20e), with the CURRENT coefficients `^t a_A..^t a_E`
    // of Eq. (15e) (evaluated on `state` and frozen, so no coefficient increment
    // enters). [`gl_tying_metrics`] reports the paper's symmetric metric form
    // `^t g_ij^m = ^t x_{m,i} . ^t x_{m,j} + ^t x_{m,j} . ^t x_{m,i}`
    // (Eqs. 12b/12c), i.e. twice the Gram metric, and the tying strain of
    // Eq. (20b) is the `1/2` of the `1/2 ^t g~ - 1/2 ^0 g~` of Eq. (15a): the
    // `1/4` below is exactly that pair of factors, and it makes the `u = 0`
    // derivative of this slice equal `pre.b_rr_a[j]` etc. -- the operator
    // [`b_membrane_2017`] of the validated linear path uses.
    let (cur_vecs, _n, m_r, m_s) = gl_current_characteristic_vectors(state);
    let (_c_r, _c_s, _d, cur_coeff) =
        compute_membrane_coefficients_2017(&cur_vecs[2], &m_r, &m_s);
    let (next_vecs, _n_next, _m_r_next, _m_s_next) = gl_current_characteristic_vectors(&next);

    // Eqs. (15b-d): the five tying metrics are `2 a . a` (rr at A and B), `2 a . a`
    // (ss at C and D) and `2 a . b` (rs at E), with the Eq. (13) mid-surface
    // tangents `a = x_r + s x_d` and `b = x_s + r x_d`. Split exactly.
    let mid = |v: &[Vector3<f64>; 3], r_t: f64, s_t: f64| (v[0] + s_t * v[2], v[1] + r_t * v[2]);
    let tie_split = |r_t: f64, s_t: f64, kind: usize| -> (f64, f64) {
        let (a0, b0) = mid(&cur_vecs, r_t, s_t);
        let (a1, b1) = mid(&next_vecs, r_t, s_t);
        let da = a1 - a0;
        let db = b1 - b0;
        match kind {
            0 => (4.0 * a0.dot(&da), 2.0 * da.dot(&da)),
            1 => (4.0 * b0.dot(&db), 2.0 * db.dot(&db)),
            _ => (2.0 * (a0.dot(&db) + da.dot(&b0)), 2.0 * da.dot(&db)),
        }
    };
    // Fig. 4's tying points A(0,1), B(0,-1), C(1,0), D(-1,0), E(0,0) with the
    // metric component each one supplies: rr, rr, ss, ss, rs.
    let tying_points: [(f64, f64, usize); 5] = [
        (0.0, 1.0, 0),
        (0.0, -1.0, 0),
        (1.0, 0.0, 1),
        (-1.0, 0.0, 1),
        (0.0, 0.0, 2),
    ];
    let mut dtie_lin = [0.0f64; 5];
    let mut dtie_nl = [0.0f64; 5];
    for (k, &(r_t, s_t, kind)) in tying_points.iter().enumerate() {
        let (lin, nl) = tie_split(r_t, s_t, kind);
        dtie_lin[k] = 0.25 * lin;
        dtie_nl[k] = 0.25 * nl;
    }
    let am_lin = gl_assumed_mid_metric(&cur_coeff, &dtie_lin, r, s);
    let am_nl = gl_assumed_mid_metric(&cur_coeff, &dtie_nl, r, s);
    // Engineering shear `2 _0 e~_12` for the Eq. (23) mapping, matching
    // [`b_membrane_2017`].
    let m_loc_lin = tmap * Vector3::new(am_lin[0], am_lin[1], 2.0 * am_lin[2]);
    let m_loc_nl = tmap * Vector3::new(am_nl[0], am_nl[1], 2.0 * am_nl[2]);
    // Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (22a): the drill-membrane
    // strain is added to the membrane row (linear, constant operator), so it lands
    // entirely in the linear part.
    let md = b_drill_membrane_2025(pre, r, s) * u_local;

    // Eq. (20a): the exact covariant through-thickness slices, engineering shear
    // triple `[e_rr, e_ss, 2 e_rs]` for the Eq. (23) mapping. The reference
    // configuration is `state` (`^t`), so this is the increment of Eqs. (7)/(8)
    // and not the total TL strain about `^0`.
    let cov_parts = |zeta: f64| -> (Vector3<f64>, Vector3<f64>) {
        let (g_r, g_s, _) =
            compute_j3d_enriched(&next.coords, &next.vn, &state.a_i, r, s, zeta);
        let (g0_r, g0_s, _) =
            compute_j3d_enriched(&state.coords, &state.vn, &state.a_i, r, s, zeta);
        let dr = g_r - g0_r;
        let ds = g_s - g0_s;
        (
            Vector3::new(
                g0_r.dot(&dr),
                g0_s.dot(&ds),
                g0_r.dot(&ds) + dr.dot(&g0_s),
            ),
            Vector3::new(0.5 * dr.dot(&dr), 0.5 * ds.dot(&ds), dr.dot(&ds)),
        )
    };
    let (e0_lin, e0_nl) = cov_parts(0.0);
    let (ep_lin, ep_nl) = cov_parts(1.0);
    let (em_lin, em_nl) = cov_parts(-1.0);
    let b1_lin = tmap * (0.5 * (ep_lin - em_lin));
    let b1_nl = tmap * (0.5 * (ep_nl - em_nl));
    let b2_lin = tmap * (0.5 * (ep_lin + em_lin) - e0_lin);
    let b2_nl = tmap * (0.5 * (ep_nl + em_nl) - e0_nl);

    // Eq. (19): assumed transverse shear on the incremental tying shears.
    let tie_rz_parts = |rt: f64, st: f64| -> (f64, f64) {
        let (g_r, _, g_t) =
            compute_j3d_enriched(&next.coords, &next.vn, &state.a_i, rt, st, 0.0);
        let (g0_r, _, g0_t) =
            compute_j3d_enriched(&state.coords, &state.vn, &state.a_i, rt, st, 0.0);
        let dr = g_r - g0_r;
        let dt = g_t - g0_t;
        (0.5 * (g0_r.dot(&dt) + dr.dot(&g0_t)), 0.5 * dr.dot(&dt))
    };
    let tie_sz_parts = |rt: f64, st: f64| -> (f64, f64) {
        let (_, g_s, g_t) =
            compute_j3d_enriched(&next.coords, &next.vn, &state.a_i, rt, st, 0.0);
        let (_, g0_s, g0_t) =
            compute_j3d_enriched(&state.coords, &state.vn, &state.a_i, rt, st, 0.0);
        let ds = g_s - g0_s;
        let dt = g_t - g0_t;
        (0.5 * (g0_s.dot(&dt) + ds.dot(&g0_t)), 0.5 * ds.dot(&dt))
    };
    let rz_top = tie_rz_parts(0.0, 1.0);
    let rz_bot = tie_rz_parts(0.0, -1.0);
    let sz_right = tie_sz_parts(1.0, 0.0);
    let sz_left = tie_sz_parts(-1.0, 0.0);
    let (w_top, w_bot) = (0.5 * (1.0 + s), 0.5 * (1.0 - s));
    let (w_right, w_left) = (0.5 * (1.0 + r), 0.5 * (1.0 - r));
    let e_rz_lin = w_top * rz_top.0 + w_bot * rz_bot.0;
    let e_rz_nl = w_top * rz_top.1 + w_bot * rz_bot.1;
    let e_sz_lin = w_right * sz_right.0 + w_left * sz_left.0;
    let e_sz_nl = w_right * sz_right.1 + w_left * sz_left.1;
    // Eq. (23): mapped with the `^0` metric and `^0` frame, as the linear path's
    // [`b_shear_mitc4`] does.
    let (g0_r0, g0_s0, g0_t0) =
        compute_j3d_enriched(&pre.initial_coords_3d, &pre.vn, &pre.a_i, r, s, 0.0);
    let tshear = shear_covariant_to_local(&g0_r0, &g0_s0, &g0_t0, &pre.e1, &pre.e2, &pre.e3);
    let gamma_lin = tshear * Vector2::new(e_rz_lin, e_sz_lin);
    let gamma_nl = tshear * Vector2::new(e_rz_nl, e_sz_nl);

    (
        [
            m_loc_lin[0] + md[0],
            m_loc_lin[1] + md[1],
            m_loc_lin[2] + md[2],
            s1 * b1_lin[0],
            s1 * b1_lin[1],
            s1 * b1_lin[2],
            s2 * b2_lin[0],
            s2 * b2_lin[1],
            s2 * b2_lin[2],
            gamma_lin[0],
            gamma_lin[1],
        ],
        [
            m_loc_nl[0],
            m_loc_nl[1],
            m_loc_nl[2],
            s1 * b1_nl[0],
            s1 * b1_nl[1],
            s1 * b1_nl[2],
            s2 * b2_nl[0],
            s2 * b2_nl[1],
            s2 * b2_nl[2],
            gamma_nl[0],
            gamma_nl[1],
        ],
    )
}

/// [`n_gamma_local_strain_parts`]' two halves summed: the strain vector
/// `[m(3), b1(3), b2(3), gamma(2)]` consumed by the internal force and by the
/// geometric term's stress.
fn n_gamma_local_strain(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
    r: f64,
    s: f64,
) -> [f64; 11] {
    let (lin, nl) = n_gamma_local_strain_parts(pre, state, u_local, r, s);
    core::array::from_fn(|a| lin[a] + nl[a])
}

/// The next configuration of one local 24-DOF increment `u_local`: Eq. (4b) for
/// the nodal geometry and Eq. (26) for the directors.
///
/// FRAME. `u_local` is the LOCAL-frame increment (`[`build_t24`]` is
/// `u_local = T u_global`), while `state.coords`, `state.vn`, `state.v1`,
/// `state.v2` are physical GLOBAL-frame vectors (`pre.initial_coords_3d` and the
/// Eq. (1)/(3d) nodal directors; see ADR-2). Eq. (4b) prints
/// `^{t+dt}x_i = ^t x_i + u_i` and Eq. (5b) `u_m = sum_i h_i u_i` with `u_i` a
/// *vector*, so its representation must be carried into the frame of `^t x_i`
/// before it is added: `u_i = u_x e1 + u_y e2 + u_z e3`. Both the translational
/// and the rotational triples are therefore mapped through the local frame, and
/// `d(next)/d(u_local)` is then the element's LINEAR operator [`b_membrane_2017`]
/// etc. Adding the local components to the global coordinates directly (the
/// earlier form) instead made `B(0)` the operator of a fictitious element whose
/// translational DOFs act along the global axes: measured at `u = 0`, relatives
/// `1.2e-2` (membrane), `2.0e-5` (b1) and `1.2e-2` (shear) on a warped element
/// and `4.5e-1` on a 60-degree-rotated flat element, while an axis-aligned flat
/// element was unaffected.
fn gl_state_after_local(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
) -> GlCurrentState {
    let mut next = *state;
    for i in 0..4 {
        let du = pre.e1 * u_local[6 * i]
            + pre.e2 * u_local[6 * i + 1]
            + pre.e3 * u_local[6 * i + 2];
        for k in 0..3 {
            next.coords[i][k] += du[k];
        }
        let theta = pre.e1 * u_local[6 * i + 3]
            + pre.e2 * u_local[6 * i + 4]
            + pre.e3 * u_local[6 * i + 5];
        let q = Mitc4PlusDPrecomputed::quaternion_from_vector(&theta);
        next.vn[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.vn[i], &q);
        next.v1[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.v1[i], &q);
        next.v2[i] = Mitc4PlusDPrecomputed::rotate_vector_by_quaternion(&state.v2[i], &q);
    }
    next
}

/// The 11x11 constitutive-resultant block of the N-gamma element: the nine-row
/// through-thickness moment matrix [`resultant_moment_matrix`] (Eq. 16's
/// `[1, zeta, zeta^2]` moments) in the first nine rows, and the uncorrected
/// transverse-shear stiffness `G h` of ADR-1 in the last two.
fn n_gamma_w11(pre: &Mitc4PlusDPrecomputed) -> SMatrix<f64, 11, 11> {
    let w9 = resultant_moment_matrix(&pre.constitutive, pre.thickness);
    let mut w = SMatrix::<f64, 11, 11>::zeros();
    for i in 0..9 {
        for j in 0..9 {
            w[(i, j)] = w9[(i, j)];
        }
    }
    for i in 0..2 {
        for j in 0..2 {
            w[(9 + i, 9 + j)] = pre.cs_uncorrected[(i, j)];
        }
    }
    w
}

/// The strain-displacement matrix `B = d(_0 e~)/du` (11x24) of the N-gamma
/// strain vector, obtained by central differences of [`n_gamma_local_strain`]
/// with a step of `H = 2.0e-5`. This is the consistent linearization of the
/// printed incremental strain expressions (the paper prints Eq. (25)
/// `_0 e_ij = B_ij U_e` but not `B_ij`, p. 7 right).
///
/// The step is the round-off/truncation balance of that difference: the printed
/// strain is a difference of `O(1)` metric quantities, so its absolute round-off
/// is `~eps` independently of the strain value, and the resulting `B` carries
/// `~eps/H` noise against `O(H^2)` truncation. Measured on the `u = 0` assembly
/// residual of `n_gamma_rigid_body_zero_force_and_consistent_tangent` (relative):
/// `3.2e-10` at `H = 1e-6` (noise) and `5.5e-10` at `H = 1e-4` (truncation), so
/// the two terms balance near `H = 2e-5`, where the residual is `7e-11` -- under
/// the `1e-10` the measurement asks for, and six orders below its `1e-6`
/// consistency gate.
///
/// At `u = 0` this matrix is the element's LINEAR strain operator: the strain of
/// [`n_gamma_local_strain`] is an increment about `state` with the assumed
/// field's coefficients and tying values frozen at `state`, so its central
/// difference at `u = 0` reproduces the linear rows of the element
/// (Eqs. 27a-c / 7c / 7d of C&S 182:404-418 plus the Eq. 26 drill row of
/// C&S 308:107622), to round-off.
fn n_gamma_b_matrix(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
    r: f64,
    s: f64,
) -> SMatrix<f64, 11, 24> {
    n_gamma_b_matrix_at(pre, state, u_local, r, s, N_GAMMA_B_H)
}

/// The step of [`n_gamma_b_matrix`], i.e. of the `B` used by the internal force
/// and by the material term of the tangent. It is pinned from BELOW by the
/// `u = 0` identity the element has to satisfy (`K_T(0) == K_0` to `<= 1e-10`
/// relative), so it cannot be relaxed to buy accuracy for the geometric term;
/// that term uses its own step, [`N_GAMMA_GEO_H`].
const N_GAMMA_B_H: f64 = 2.0e-5;

/// [`n_gamma_b_matrix`] with an explicit central-difference step. Only the
/// geometric term of [`n_gamma_kt_local`] passes a step other than
/// [`N_GAMMA_B_H`].
fn n_gamma_b_matrix_at(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
    r: f64,
    s: f64,
    h: f64,
) -> SMatrix<f64, 11, 24> {
    let mut b = SMatrix::<f64, 11, 24>::zeros();
    for j in 0..24 {
        let mut up = *u_local;
        up[j] += h;
        let mut um = *u_local;
        um[j] -= h;
        let ep = n_gamma_local_strain(pre, state, &up, r, s);
        let em = n_gamma_local_strain(pre, state, &um, r, s);
        for a in 0..11 {
            b[(a, j)] = 0.5 * (ep[a] - em[a]) / h;
        }
    }
    b
}

/// The central-difference step of the geometric term of [`n_gamma_kt_local`],
/// used by BOTH of its nested differences.
///
/// THE STEP CANNOT BE SPLIT, and this is a measurement, not a preference. The
/// geometric term is the Hessian `sum_a (d2 _0 e~_a / du_i du_j) S_a`, evaluated
/// as `FD_{H_o}(FD_{H_i}(_0 e~))`. Written out, that composition IS the
/// four-point cross difference
///
/// ```text
/// [e(u + H_i e_i + H_o e_j) - e(u - H_i e_i + H_o e_j)
///  - e(u + H_i e_i - H_o e_j) + e(u - H_i e_i - H_o e_j)] / (4 H_i H_o),
/// ```
///
/// whose round-off bound is `eps |e_us| / (H_i H_o)` -- the four evaluations'
/// round-off over `4 H_i H_o` -- against `(H_i^2/12) |e''''|` truncation. The
/// form this replaced used `(H_i, H_o) = (2e-5, 1e-6)`, whose bound is
/// `eps/(H_i H_o) = 1.1e-5`; the measured asymmetry was `3.4e-6 .. 8.1e-6`, i.e.
/// `0.2 .. 0.8` of that bound, on all five fixtures of
/// `n_gamma_geo_symmetry_diagnostic`. With both steps at `h` the bound becomes
/// `eps/h^2` against `O(h^2)`, whose optimum is the classic
/// `h* = (12 eps)^(1/4) ~ 2.3e-4` with total error `~eps^(1/2) ~ 1e-8`.
///
/// WHY THIS IS NOT [`N_GAMMA_B_H`], and why the two optima differ. `B` is a
/// FIRST central difference: its balance is `eps|e_us|/h` against `(h^2/6)|e'''|`,
/// whose optimum is the classic `h* = (6 eps)^(1/3) ~ 1.1e-5` -- the order of
/// [`N_GAMMA_B_H`], and also what the `u = 0` identity requires (the measured
/// minimum near `H = 2e-5` recorded in [`n_gamma_b_matrix`]). This term is a
/// SECOND (cross) difference, one order worse in both directions, so its optimum
/// is about twenty times larger. Sharing one step with `B` is exactly what made
/// the composition's floor `eps/(H_i H_o)`.
///
/// Instrument `n_gamma_geo_stencil_probe` (local frame): the single-step
/// estimates at `h = 3e-5, 1e-4, 3e-4` agree with each other to `3e-8` while
/// differing from the split-step form by `7e-6`, FLAT in `h` (`6.655e-6`,
/// `6.680e-6`, `6.684e-6`). A truncation-limited stencil would move by four
/// orders of magnitude over that range, so the split-step form is the inaccurate
/// one and this step gives the Hessian to about `3e-8` -- two orders inside the
/// `1e-6` tangent-consistency gate.
///
/// A consequence worth keeping: both index orders now combine the SAME four
/// points, so `geo` is exactly symmetric in `(i, j)` by construction instead of
/// to round-off. Instrument `n_gamma_geo_symmetry_diagnostic`.
const N_GAMMA_GEO_H: f64 = 1.0e-4;

/// The faithful total-Lagrangian internal force of Eq. (24b),
/// `^t_0 F_e = int_{0V} B^T ^t_0 S d0V`, on the local 24-DOF increment from the
/// current configuration `state`. The 2x2 surface rule carries
/// `^t_0 S = W _0 e~` with the through-thickness moments of [`n_gamma_w11`], and
/// `B` is the Eq. (25) matrix of [`n_gamma_b_matrix`] built from the
/// [`n_gamma_local_strain`] increment of Eqs. (19)-(22).
///
/// At `u = 0` the force is exactly zero, and for any rigid-body field it is
/// zero to round-off (see `n_gamma_rigid_body_zero_force_and_consistent_tangent`).
pub fn n_gamma_fint_local(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
) -> Vec24 {
    let w = n_gamma_w11(pre);
    let mut f = Vec24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sg = surface_measure(pre, r, s);
        let mut e = SVector::<f64, 11>::zeros();
        for (a, v) in n_gamma_local_strain(pre, state, u_local, r, s).iter().enumerate() {
            e[a] = *v;
        }
        let b = n_gamma_b_matrix(pre, state, u_local, r, s);
        f += (b.transpose() * (w * e)) * (GAUSS_W[g] * sg);
    }
    f
}

/// The faithful total-Lagrangian tangent of Eq. (24a) on the local 24-DOF
/// increment from `state`:
///
/// ```text
/// ^t K_e = int B^T C B d0V + int (dB/du)^T S d0V
/// ```
///
/// The material term `int B^T C B d0V` is [`n_gamma_w11`] sandwiched by
/// [`n_gamma_b_matrix`]; the geometric `N`-term is the central difference of `B`
/// contracted with `S = W _0 e~`, i.e. exactly `int (dB/du)^T S d0V`, taken with
/// the single step [`N_GAMMA_GEO_H`] in both of its nested differences. The
/// split-step form this replaced carried an `eps/(H_i H_o) = 1.1e-5` round-off
/// bound (measured `3.4e-6 .. 8.1e-6`) and made
/// `K_t != dF/du` on every element whose local frame is not axis-aligned -- the
/// producer of the `SNES diverged` large-rotation failures.
///
/// The returned matrix is the derivative of [`n_gamma_fint_local`]:
/// `F(u) = int B(u)^T W e~(u) d0V`, so
/// `dF/du = int (dB/du)^T W e~ d0V + int B^T W (d e~/du) d0V`. The two terms are
/// the contraction above and `mat = int B^T W B d0V`; the paper prints Eq. (25)
/// but not `B`, and [`n_gamma_b_matrix`]'s finite-difference `B` is what both
/// terms use, so the pair is consistent by construction.
///
/// There is deliberately NO additive `compute_ke_local(pre)` reference. An
/// earlier construction returned `compute_ke_local(pre) + (mat - mat0)` with
/// `mat0 = int B(0)^T W11 B(0)`, which forces `K_t(0) == K_0` bit for bit, but
/// `mat0` equals `compute_ke_local(pre)` only for a FLAT element; on a warped
/// (non-planar) element the two differ by a state-constant matrix (measured
/// relative `2.95e-3` on a warped fixture), and adding that constant back makes
/// `K_t` differ from `dF/du` by the same constant -- which broke the incremental
/// solver's Newton line search. `K_t(0)` is now the exact elementwise `u = 0`
/// limit of `dF/du`, namely `mat0`; on a flat element it reproduces the linear
/// stiffness to round-off (relative `2.3e-11` on `RECT`, see
/// `n_gamma_rigid_body_zero_force_and_consistent_tangent`).
pub fn n_gamma_kt_local(
    pre: &Mitc4PlusDPrecomputed,
    state: &GlCurrentState,
    u_local: &Vec24,
) -> Mat24 {
    // Both nested differences use the SAME step: the floor of the composition is
    // `eps / (H_i H_o)`, so splitting the steps cannot buy accuracy. See
    // `N_GAMMA_GEO_H` for the measurement behind the value.
    const H: f64 = N_GAMMA_GEO_H;
    let w = n_gamma_w11(pre);
    let mut mat = Mat24::zeros();
    let mut geo = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sg = surface_measure(pre, r, s);
        let wq = GAUSS_W[g] * sg;
        let mut e = SVector::<f64, 11>::zeros();
        for (a, v) in n_gamma_local_strain(pre, state, u_local, r, s).iter().enumerate() {
            e[a] = *v;
        }
        let b = n_gamma_b_matrix(pre, state, u_local, r, s);
        mat += (b.transpose() * w * b) * wq;
        // Geometric term: `geo[i][j] = sum_a (dB[a][i]/du[j]) S[a]`. The
        // quadrature weight `wq` is the same one the material term uses; without
        // it the term is the unweighted Gauss-point sum, which is not
        // `int (dB/du)^T S d0V` for a general element and is off by the
        // (constant) measure factor on a flat one.
        let svec = w * e;
        for j in 0..24 {
            let mut up = *u_local;
            up[j] += H;
            let mut um = *u_local;
            um[j] -= H;
            let bp = n_gamma_b_matrix_at(pre, state, &up, r, s, H);
            let bm = n_gamma_b_matrix_at(pre, state, &um, r, s, H);
            for a in 0..11 {
                for i in 0..24 {
                    geo[(i, j)] += 0.5 * (bp[(a, i)] - bm[(a, i)]) / H * svec[a] * wq;
                }
            }
        }
    }
    mat + geo
}

/// Global-coordinate wrapper of [`n_gamma_fint_local`] for the initial element
/// geometry: `f_global = T^T f_local`, `u_local = T u_global`.
///
/// **Known limitation (not fixed here).** The `GlCurrentState` is rebuilt from
/// the INITIAL element geometry (`pre.initial_coords_3d`, `pre.vn`, `pre.v1`,
/// `pre.v2`, `pre.a_i`), so this is the first-step
/// (reference-configuration) linearisation of the formulation, not the
/// last-converged state of a load step. Threading the last-converged state is
/// the solver's job and is a follow-up.
pub fn n_gamma_fint_global(pre: &Mitc4PlusDPrecomputed, u_global: &Vec24) -> Vec24 {
    let t24 = build_t24(pre);
    let state = GlCurrentState {
        coords: pre.initial_coords_3d,
        vn: pre.vn,
        v1: pre.v1,
        v2: pre.v2,
        a_i: pre.a_i,
    };
    t24.transpose() * n_gamma_fint_local(pre, &state, &(t24 * u_global))
}

/// Global-coordinate wrapper of [`n_gamma_kt_local`] for the initial element
/// geometry: `K_global = T^T K_local T`.
///
/// **Known limitation (not fixed here).** The `GlCurrentState` is rebuilt from
/// the INITIAL element geometry (`pre.initial_coords_3d`, `pre.vn`, `pre.v1`,
/// `pre.v2`, `pre.a_i`), so this is the first-step
/// (reference-configuration) linearisation of the formulation, not the
/// last-converged state of a load step. Threading the last-converged state is
/// the solver's job and is a follow-up.
pub fn n_gamma_kt_global(pre: &Mitc4PlusDPrecomputed, u_global: &Vec24) -> Mat24 {
    let t24 = build_t24(pre);
    let state = GlCurrentState {
        coords: pre.initial_coords_3d,
        vn: pre.vn,
        v1: pre.v1,
        v2: pre.v2,
        a_i: pre.a_i,
    };
    t24.transpose() * n_gamma_kt_local(pre, &state, &(t24 * u_global)) * t24
}

/// Mid-surface area of the element: `sum_g w_g ||g_r x g_s||` under the 2x2
/// rule of Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410.
fn element_area(pre: &Mitc4PlusDPrecomputed) -> f64 {
    (0..N_GAUSS)
        .map(|g| surface_measure(pre, GAUSS_XI[g], GAUSS_ETA[g]) * GAUSS_W[g])
        .sum()
}

/// Consistent translational/rotary mass matrix in GLOBAL coordinates (24x24),
/// with `m_trans = rho h` and `m_rot = rho h^3 / 12`.
///
/// Integrates `sum_ij N_i N_j m w sqrt_g` under the 2x2 surface rule; only the
/// diagonal 3x3 blocks per node pair are filled (translations and rotations).
pub fn compute_me_global(pre: &Mitc4PlusDPrecomputed, rho: f64) -> Mat24 {
    compute_me_with_inertias(pre, rho * pre.thickness, rho * pre.thickness.powi(3) / 12.0)
}

/// Consistent mass matrix for composite elements (global coordinates, 24x24),
/// using the pre-integrated ply mass parameters instead of `rho h`.
pub fn compute_me_composite_global(
    pre: &Mitc4PlusDPrecomputed,
    mass_per_area: f64,
    rotational_inertia: f64,
) -> Mat24 {
    compute_me_with_inertias(pre, mass_per_area, rotational_inertia)
}

/// Shared body of the two consistent mass matrices.
fn compute_me_with_inertias(pre: &Mitc4PlusDPrecomputed, m_trans: f64, m_rot: f64) -> Mat24 {
    let mut m_local = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let w = GAUSS_W[g];
        let n = shape_functions(r, s);
        for i in 0..4 {
            for j in 0..4 {
                let val_t = n[i] * n[j] * m_trans * w * sqrt_g;
                for k in 0..3 {
                    m_local[(6 * i + k, 6 * j + k)] += val_t;
                }
                let val_r = n[i] * n[j] * m_rot * w * sqrt_g;
                for k in 3..6 {
                    m_local[(6 * i + k, 6 * j + k)] += val_r;
                }
            }
        }
    }
    transform_to_global(pre, &m_local)
}

/// Body-load vector (24 long) in GLOBAL coordinates from a constant body-force
/// acceleration `gravity` and density `rho`:
/// `f = sum_g N^T (rho h) g_local w sqrt_g`, transformed with `T^T`. Only the
/// translational slots receive contributions.
pub fn compute_body_load_global(
    pre: &Mitc4PlusDPrecomputed,
    rho: f64,
    gravity: &Vector3<f64>,
) -> Vec24 {
    let rho_h = rho * pre.thickness;
    let g_local = pre.t3 * gravity;

    let mut f_local = Vec24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let n = shape_functions(r, s);
        for i in 0..4 {
            let contrib = n[i] * rho_h * GAUSS_W[g] * sqrt_g;
            for k in 0..3 {
                f_local[6 * i + k] += contrib * g_local[k];
            }
        }
    }
    build_t24(pre).transpose() * f_local
}

/// Initial-stress stiffness `K_sigma` in GLOBAL coordinates from a membrane
/// resultant state `sigma = [N_xx, N_yy, N_xy]` in LOCAL coordinates.
pub fn compute_k_sigma_global(pre: &Mitc4PlusDPrecomputed, sigma_membrane: &Vector3<f64>) -> Mat24 {
    let k_local = geometric_stiffness_from_stress(pre, sigma_membrane);
    transform_to_global(pre, &k_local)
}

/// Centrifugal prestress `[N_xx, N_yy, N_xy]` in LOCAL coordinates:
/// `sigma_cf ~ rho omega^2 r_radial L_char` with `L_char = sqrt(area)`. The
/// repository's centrifugal model, not a paper equation.
#[allow(clippy::too_many_arguments)]
pub fn compute_centrifugal_prestress(
    pre: &Mitc4PlusDPrecomputed,
    omega: f64,
    rotation_axis: &Vector3<f64>,
    rotation_center: &Vector3<f64>,
    centroid: &Vector3<f64>,
    rho: f64,
) -> Vector3<f64> {
    let axis = rotation_axis.normalize();
    let r_vec = centroid - rotation_center;
    let r_parallel = r_vec.dot(&axis) * axis;
    let r_radial_vec = r_vec - r_parallel;
    let r_radial = r_radial_vec.norm();
    if r_radial < 1.0e-10 {
        return Vector3::zeros();
    }
    let radial_dir = r_radial_vec / r_radial;
    let l_char = element_area(pre).sqrt();
    let sigma_cf = rho * omega * omega * r_radial * l_char;

    let radial_local = pre.t3 * radial_dir;
    let cos_theta = radial_local[0];
    let sin_theta = radial_local[1];
    Vector3::new(
        sigma_cf * cos_theta * cos_theta,
        sigma_cf * sin_theta * sin_theta,
        sigma_cf * cos_theta * sin_theta,
    )
}

/// Element stress and strain at the element centre `(r, s) = (0, 0)` in the
/// element's LOCAL frame, `([sx, sy, sxy, 0, 0, 0], [exx, eyy, exy, 0, 0, 0])`.
///
/// The membrane strain is `B_m u` from Ko, Lee & Bathe (2017), C&S 182:404-418,
/// Eqs. (27a-c); the through-thickness strain is the paper's Eq. (7a)
/// decomposition `t e^b1 + t^2 e^b2` with `t = 2 z_factor`, so `z_factor = 0` is
/// the mid-surface and `z_factor = +/- 1/2` the faces. `stress_type`: 0 =
/// membrane, 1 = bending, 2 = total.
pub fn compute_element_stress(
    pre: &Mitc4PlusDPrecomputed,
    u_global: &Vec24,
    z_factor: f64,
    stress_type: u8,
) -> ([f64; 6], [f64; 6]) {
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;

    let cm_raw = &pre.constitutive.cm_raw;
    let bm = b_membrane_2017(pre, 0.0, 0.0);
    let eps_m = bm * u_local;
    let sig_m = cm_raw * eps_m;

    let (bb1, bb2) = b_bending_2017(pre, 0.0, 0.0);
    let t = 2.0 * z_factor;
    let eps_b = t * (bb1 * u_local) + (t * t) * (bb2 * u_local);
    let sig_b = cm_raw * eps_b;

    let sig = match stress_type {
        0 => sig_m,
        1 => sig_b,
        _ => sig_m + sig_b,
    };
    let eps = match stress_type {
        0 => eps_m,
        1 => eps_b,
        _ => eps_m + eps_b,
    };

    (
        [sig[0], sig[1], 0.0, sig[2], 0.0, 0.0],
        [eps[0], eps[1], 0.0, eps[2], 0.0, 0.0],
    )
}

/// Extract a 24-DOF element displacement from a global vector and a DOF map.
pub fn extract_elem_disp_24(u: &[f64], dofs: &[usize]) -> Vec24 {
    let mut ue = Vec24::zeros();
    for (local, &global) in dofs.iter().enumerate() {
        ue[local] = u[global];
    }
    ue
}

// ============================================================================
// Finite-rotation helpers — the paper's Eq. (26)
// ============================================================================
//
// Ko, Lee & Bathe (2017), C&S 185:1-14, Eq. (26): the director frame of the
// current configuration is the exact finite rotation of the reference frame,
// applied to `^t V_n^i`, `^t V_1^i`, `^t V_2^i`. `gl_state_after_local` is the
// only consumer of these helpers.
//
// The corotational scaffolding that used to live here (`GpLocalFrame`, polar
// decomposition, logarithmic strain, `update_corotational_frame`,
// `frame_incremental_rotation`, `update_normals_with_displacements`) was
// unsourced — the paper is total Lagrangian, not corotational — and had zero
// call sites; N-delta deleted it.

impl Mitc4PlusDPrecomputed {
    /// Rotation matrix from a unit quaternion `q = [q0, qx, qy, qz]`:
    /// `R = I + 2 q0 [q^] + 2 [q^]^2`.
    pub fn quaternion_to_matrix(q: &Vector4<f64>) -> Matrix3<f64> {
        let q0 = q[0];
        let qx = q[1];
        let qy = q[2];
        let qz = q[3];
        let q_hat = Matrix3::new(0.0, -qz, qy, qz, 0.0, -qx, -qy, qx, 0.0);
        Matrix3::identity() + (2.0 * q0 * q_hat) + (2.0 * q_hat * q_hat)
    }

    /// Quaternion from a rotation vector `theta = [theta_x, theta_y, theta_z]`.
    pub fn quaternion_from_vector(theta: &Vector3<f64>) -> Vector4<f64> {
        let theta_norm = theta.norm();
        if theta_norm < 1e-15 {
            return Vector4::new(1.0, 0.0, 0.0, 0.0);
        }
        let half_angle = 0.5 * theta_norm;
        let sin_half = half_angle.sin();
        let cos_half = half_angle.cos();
        let inv_norm = 1.0 / theta_norm;
        Vector4::new(
            cos_half,
            sin_half * theta[0] * inv_norm,
            sin_half * theta[1] * inv_norm,
            sin_half * theta[2] * inv_norm,
        )
    }

    /// Rotate a vector by a quaternion: `v_new = R(q) v`.
    pub fn rotate_vector_by_quaternion(v: &Vector3<f64>, q: &Vector4<f64>) -> Vector3<f64> {
        Self::quaternion_to_matrix(q) * v
    }
}

#[cfg(test)]
mod tests {
    use super::{
        b_bending_2017, b_bending_covariant_2017, b_drill_membrane_2025, b_membrane_2017,
        b_membrane_covariant_2017, b_shear_mitc4, build_t24, compute_body_load_global,
        compute_centrifugal_prestress, compute_fint_global,
        compute_fint_global_bounded_superseded, compute_k_sigma_global,
        compute_ke_global, compute_ke_local, compute_ke_local_with_drill, compute_kt_global,
        compute_me_composite_global, compute_me_global, compute_membrane_coefficients_2017,
        covariant_to_local_mapping, director_increment_quadratic, drill_jacobian_ratio,
        drill_ke_local, drill_midside_shape_derivatives, element_area,
        geometric_stiffness_from_stress, gl_assumed_membrane_increment,
        gl_assumed_membrane_strain, gl_assumed_mid_metric, gl_assumed_transverse_shear,
        gl_current_base_vectors, gl_current_characteristic_vectors,
        gl_strain_increment, gl_strain_increment_components, gl_tying_metrics,
        incremental_disp_gradients, incremental_disp_split,
        interpolate_displacement, interpolate_position, j_loc_at, membrane_ke_local,
        n_gamma_b_matrix, n_gamma_b_matrix_at, n_gamma_fint_global, n_gamma_fint_local,
        n_gamma_kt_global,
        n_gamma_kt_local, n_gamma_w11, node_vec,
        resultant_moment_matrix, shape_function_derivatives, shape_functions, shear_ke_local,
        surface_measure, transform_to_global, GlCurrentState, GlIncrement, Mat24,
        Mitc4PlusDPrecomputed, Vec24, DRILL_EDGE_MID, GAUSS_ETA, GAUSS_W, GAUSS_XI, NODE_ETA,
        NODE_XI, N_GAUSS, N_GAMMA_GEO_H,
    };
    use crate::materials::laminate::{Laminate, Ply};
    use crate::materials::orthotropic::OrthotropicMaterial;
    use crate::materials::{isotropic::IsotropicMaterial, Material, ShellConstitutive};
    use nalgebra::{DMatrix, DVector, Matrix2, Matrix3, SMatrix, Vector3};

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
        //
        // WU6d note on `theta_z alone`: after the local frame's `e3` was corrected
        // to the paper's Eq. (10) normal, `V^D = e3` identically
        // (Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (5)/(10); Ko, Lee &
        // Bathe (2017), C&S 182:404-418, Eq. (10)), so `theta^D = theta . V^D`
        // IS the local drill slot. The implementation risk the row guards is
        // therefore "use the GLOBAL vertical instead of `V^D`", which is what
        // `local_components(&pre, &e3)` with `e3 = (0,0,1)_global` expresses.
        // Passing the raw global vector as local components would be vacuous now.
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
                "theta_z alone (global z instead of V^D)",
                drill::b_md_parameterised(
                    &pre,
                    GP,
                    GP,
                    &vd,
                    &local_components(&pre, &e3),
                    0.125,
                    &[0, 1, 2, 3],
                    1.0,
                ),
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

        /// The resultant/moment matrix `W`, recomputed independently. `W_22` is
        /// the raw 4th moment `int z^4 C dz` under the paper's 2x2 t-rule
        /// (`cm h^4/144`); the `(4/h^2)` of `E2` is carried by `b` below.
        fn moment_matrix(
            c: &ShellConstitutive,
            thickness: f64,
            include_w22: bool,
        ) -> SMatrix<f64, 9, 9> {
            let w22 = c.cm * (thickness.powi(4) / 144.0);
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
                        w[(6 + i, 6 + j)] = w22[(i, j)];
                    }
                }
            }
            w
        }

        fn assemble(pre: &Mitc4PlusDPrecomputed, g: &Geom, include_w22: bool) -> Mat24 {
            let w = moment_matrix(&pre.constitutive, pre.thickness, include_w22);
            let s1 = 2.0 / pre.thickness;
            let s2 = 4.0 / (pre.thickness * pre.thickness);
            let mut k = Mat24::zeros();
            for gi in 0..4 {
                let (r, s) = (GAUSS_XI[gi], GAUSS_ETA[gi]);
                let gr = g.x_r + s * g.x_d;
                let gs = g.x_s + r * g.x_d;
                let sqrt_g = gr.cross(&gs).norm();
                // Ko, Bathe & Zhang (2025), C&S 308:107622, Eq. (22a), p. 10:
                // the drill-membrane strain is SUMMED into the membrane row, so
                // the single `W` block carries the paper's cross terms too.
                let mut bm = membrane_local(pre, g, r, s);
                bm += drill::b_md_reference(pre, r, s);
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
        let w = resultant_moment_matrix(&c, h);
        assert_eq!(w.shape(), (9, 9));
        // W_22 is the raw 4th moment `int z^4 C dz` under the paper's 2x2
        // t-rule; the `(4/h^2)` of `E2` is carried by the B-operator, so the
        // **effective** E2-E2 coefficient is `(4/h^2)^2 W_22 = cm/9`.
        let w22 = c.cm * (h.powi(4) / 144.0);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(w[(i, j)], c.cm[(i, j)], "W_00 = cm");
                assert_eq!(w[(3 + i, 3 + j)], c.cb[(i, j)], "W_11 = cb");
                assert_eq!(w[(i, 3 + j)], c.cb_coupling[(i, j)], "W_01 = cb_coupling");
                assert_eq!(w[(3 + i, j)], c.cb_coupling[(j, i)], "W_10 = cb_coupling^T");
                assert_eq!(w[(i, 6 + j)], c.cb[(i, j)], "W_02 = cb");
                assert_eq!(w[(6 + i, j)], c.cb[(j, i)], "W_20 = cb^T");
                assert_eq!(w[(6 + i, 6 + j)], w22[(i, j)], "W_22 = int z^4 C dz");
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
            (w[(6, 6)] - cm00 * h.powi(4) / 144.0).abs() <= 1e-12 * (cm00 * h.powi(4) / 144.0),
            "W_22 = cm h^4/144 (int z^4 C dz under the 2x2 t-rule)"
        );
        // The paper's 2x2 t-rule gives int t^4 dt = 2/9 (not the exact 2/5).
        let t = 1.0 / 3.0f64.sqrt();
        let t4 = 2.0 * t.powi(4);
        assert!((t4 - 2.0 / 9.0).abs() <= 1e-15, "2-point t-rule gives 2/9");
        // The effective E2-E2 coefficient the stiffness sees is cm/9, not the
        // exact cm/5: negligible at h=0.02, so assert it on the scaled value.
        let s2 = 4.0 / (h * h);
        let eff = s2 * s2 * w[(6, 6)];
        assert!(
            (eff - cm00 / 9.0).abs() <= 1e-9 * cm00,
            "effective E2-E2 coefficient must be the paper's cm/9, got {eff}"
        );
        assert!(
            (eff - cm00 / 5.0).abs() > 1e-3 * (cm00 / 5.0),
            "effective E2-E2 must be the paper's cm/9, not the exact cm/5"
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
        // (c) Eq. (22a) structure: the drill-membrane strain is SUMMED into the
        // membrane row, so the drill's stiffness difference is
        //   `B_m^T W B_md + B_md^T W B_m + B_md^T W B_md`
        // and therefore couples the translational and rotational slots. Two
        // exact consequences are asserted (the old "only the rotational blocks
        // are touched" form belonged to the superseded separate-block split):
        //   (i) every TRANSLATIONAL-translational block is still exactly zero —
        //       no term of the folded difference has both factors on
        //       translational slots; and
        //   (ii) the translational-rotational coupling is NON-VACUOUS: the
        //       paper's membrane-drill cross term is present. Its absence was
        //       exactly the deviation this identity lock used to miss, because
        //       its reference shared the split.
        let mut coupling = 0.0f64;
        for i in 0..4 {
            for j in 0..4 {
                for a in 0..3 {
                    for b in 0..3 {
                        assert_eq!(diff[(6 * i + a, 6 * j + b)], 0.0, "t-t block");
                        coupling = coupling
                            .max(diff[(6 * i + a, 6 * j + 3 + b)].abs())
                            .max(diff[(6 * i + 3 + a, 6 * j + b)].abs());
                    }
                }
            }
        }
        assert!(
            coupling > 1e-10 * scale,
            "the Eq. (22a) membrane-drill coupling is missing: max |t-r| = {coupling}"
        );
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

    // ========================================================================
    // WU5 — the assembly-facing API: consistency and mass invariants
    // ========================================================================
    //
    // The T2B consistency set (Requirement 2 / Requirement 12) and the mass
    // invariants (Requirement 12), retargeted from `Mitc4Precomputed` to
    // `Mitc4PlusDPrecomputed`. The consistency guards are load-bearing: a
    // partial wiring of the nonlinear path or the tangent makes them fail.

    /// Total translational mass per direction: the sum of the 4x4 sub-block of
    /// one direction. For a consistent mass matrix built from a partition of
    /// unity this equals `rho h A` exactly, for every element type.
    fn translational_mass_per_direction(m: &Mat24) -> [f64; 3] {
        let mut totals = [0.0; 3];
        for (d, total) in totals.iter_mut().enumerate() {
            let mut sum = 0.0;
            for i in (d..24).step_by(6) {
                for j in (d..24).step_by(6) {
                    sum += m[(i, j)];
                }
            }
            *total = sum;
        }
        totals
    }

    /// Total rotary-inertia mass per rotation direction: `rho h^3 / 12 * A`.
    fn rotary_mass_per_direction(m: &Mat24) -> [f64; 3] {
        let mut totals = [0.0; 3];
        for (k, total) in totals.iter_mut().enumerate() {
            let mut sum = 0.0;
            for i in (3 + k..24).step_by(6) {
                for j in (3 + k..24).step_by(6) {
                    sum += m[(i, j)];
                }
            }
            *total = sum;
        }
        totals
    }

    // ------------------------------------------------------------------
    // 6.1 — K_T / f_int consistency (the T2B consistency oracle)
    // ------------------------------------------------------------------
    #[test]
    fn test_kt_zero_matches_ke() {
        let pre = pre_from(&FLAT_SQUARE);
        let u_zero = Vec24::zeros();

        let k_linear = compute_ke_local(&pre);
        let k_t_zero = compute_kt_global(&pre, &u_zero);

        let t24 = build_t24(&pre);
        let k_linear_global = t24.transpose() * k_linear * t24;

        let diff = k_t_zero - k_linear_global;
        // `K_T(0)` is the exact elementwise `u = 0` limit of `dF/du`,
        // `int B(0)^T W11 B(0)`, with NO additive `compute_ke_local` reference
        // (see `n_gamma_kt_local`). The earlier additive reference made the two
        // bit-identical but broke `K_t = dF/du` on warped elements; the
        // finite-difference limit reproduces the closed-form linear stiffness to
        // a tight relative bound instead.
        let rel = diff.norm() / k_linear_global.norm().max(1e-30);
        assert!(
            rel < 1e-9,
            "K_T(u=0) must match K_linear_global to round-off, rel = {rel:.3e} (diff norm = {})",
            diff.norm()
        );
    }

    #[test]
    fn test_fint_linear_nonlinear_parity() {
        let pre = pre_from(&FLAT_SQUARE);
        let k_local = compute_ke_local(&pre);
        let t24 = build_t24(&pre);

        let mut u_global = Vec24::zeros();
        u_global[0] = 0.0001;
        u_global[1] = 0.00005;
        u_global[2] = 0.00002;
        u_global[6] = 0.0001;
        u_global[7] = -0.00005;
        u_global[8] = 0.00002;
        u_global[12] = 0.0001;
        u_global[13] = 0.00005;
        u_global[14] = 0.00002;
        u_global[18] = 0.0001;
        u_global[19] = -0.00005;
        u_global[20] = 0.00002;

        let u_local = t24 * u_global;
        let f_linear_local = k_local * u_local;

        let f_nonlinear_global = compute_fint_global(&pre, &u_global, true);
        let f_nonlinear_local = t24.transpose() * f_nonlinear_global;

        let diff = f_nonlinear_local - f_linear_local;
        let u_norm = u_local.norm();
        let rel_err = diff.norm() / f_linear_local.norm();
        assert!(
            rel_err < 1e-1,
            "f_int(nonlinear) - K u should be O(u^2): rel_err = {rel_err:.2e}, u_norm = {u_norm:.2e}"
        );
    }

    #[test]
    fn test_kt_fint_directional_derivative() {
        let pre = pre_from(&FLAT_SQUARE);

        let mut u_base = Vec24::zeros();
        for i in 0..4 {
            u_base[6 * i] = 0.0005;
            u_base[6 * i + 1] = 0.0002;
        }

        let k_t = compute_kt_global(&pre, &u_base);

        let delta = 1e-6_f64;
        let mut du = Vec24::zeros();
        du[0] = delta;
        du[1] = delta;
        du[6] = -delta;
        du[7] = delta;
        du[12] = delta;
        du[13] = -delta;
        du[18] = -delta;
        du[19] = -delta;

        let k_t_du = k_t * du;

        let f_plus = compute_fint_global(&pre, &(u_base + du), true);
        let f_base = compute_fint_global(&pre, &u_base, true);
        let f_diff = f_plus - f_base;

        let num = (k_t_du - f_diff).norm();
        let denom = f_diff.norm().max(1.0);
        let rel_err = if denom > 1e-30 { num / denom } else { 0.0 };

        assert!(
            rel_err < 0.05,
            "K_T du ~ f_int(u+du) - f_int(u): rel_err = {rel_err:.2e} (want < 5e-2)"
        );
    }

    #[test]
    fn test_kt_fint_directional_derivative_rotations() {
        let pre = pre_from(&FLAT_SQUARE);

        let mut u_base = Vec24::zeros();
        for i in 0..4 {
            u_base[6 * i] = 2.0e-4 * (i as f64 + 1.0);
            u_base[6 * i + 2] = -1.0e-4 * (i as f64 + 1.0);
            u_base[6 * i + 4] = 3.0e-4;
        }

        let k_t = compute_kt_global(&pre, &u_base);

        let delta = 1.0e-6_f64;
        let mut du = Vec24::zeros();
        du[3] = delta;
        du[4] = delta;
        du[9] = -delta;
        du[10] = delta;
        du[15] = delta;
        du[16] = -delta;
        du[21] = -delta;
        du[22] = -delta;

        let k_t_du = k_t * du;
        let f_plus = compute_fint_global(&pre, &(u_base + du), true);
        let f_base = compute_fint_global(&pre, &u_base, true);
        let f_diff = f_plus - f_base;

        let num = (k_t_du - f_diff).norm();
        let denom = f_diff.norm().max(1.0);
        let rel_err = if denom > 1e-30 { num / denom } else { 0.0 };

        assert!(
            rel_err < 0.05,
            "K_T du ~ f_int(u+du) - f_int(u) for rotational DOFs: rel_err = {rel_err:.2e} (want < 5e-2)"
        );
    }

    #[test]
    fn test_kt_fint_directional_derivative_with_drill_dofs() {
        // The translational and rotational guards above do not excite the
        // drill DOF, so none can see the drill-membrane path. This one excites
        // slot 6i+5 in both the base state and the perturbation.
        let pre = pre_from(&FLAT_SQUARE);

        let mut u_base = Vec24::zeros();
        for i in 0..4 {
            u_base[6 * i] = 2.0e-4 * (i as f64 + 1.0);
            u_base[6 * i + 1] = -1.0e-4 * (i as f64 + 1.0);
            u_base[6 * i + 5] = 3.0e-4 * (i as f64 + 1.0);
        }

        let k_t = compute_kt_global(&pre, &u_base);

        let delta = 1.0e-6_f64;
        let mut du = Vec24::zeros();
        for i in 0..4 {
            du[6 * i] = delta;
            du[6 * i + 1] = -delta;
            du[6 * i + 5] = delta;
        }

        let k_t_du = k_t * du;
        let f_plus = compute_fint_global(&pre, &(u_base + du), true);
        let f_base = compute_fint_global(&pre, &u_base, true);
        let f_diff = f_plus - f_base;

        let num = (k_t_du - f_diff).norm();
        let denom = f_diff.norm().max(1.0);
        let rel_err = if denom > 1e-30 { num / denom } else { 0.0 };

        assert!(
            rel_err < 0.05,
            "K_T du ~ f_int(u+du) - f_int(u) with the drill DOF excited: rel_err = {rel_err:.2e} (want < 5e-2)"
        );
    }

    // ------------------------------------------------------------------
    // 6.2 — mass invariants
    // ------------------------------------------------------------------
    #[test]
    fn test_me_global_is_symmetric_and_positive_semidefinite() {
        let pre = pre_from(&FLAT_SQUARE);
        let m = compute_me_global(&pre, 7800.0);

        let asymmetry = (m - m.transpose()).norm() / m.norm();
        assert!(
            asymmetry < 1e-14,
            "M_global must be symmetric: relative asymmetry = {asymmetry:.3e}"
        );

        let symmetric_part = (m + m.transpose()) * 0.5;
        let eigenvalues = nalgebra::SymmetricEigen::new(symmetric_part).eigenvalues;
        let lambda_min = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        let lambda_max = eigenvalues
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            lambda_min > -1e-12 * lambda_max,
            "M_global must be positive semi-definite: lambda_min = {lambda_min:.6e}, lambda_max = {lambda_max:.6e}"
        );
    }

    #[test]
    fn test_me_global_total_translational_mass_is_rho_h_a() {
        let pre = pre_from(&FLAT_SQUARE);
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let expected = rho * pre.thickness * element_area(&pre);

        for (d, total) in translational_mass_per_direction(&m).iter().enumerate() {
            let error = (total - expected).abs() / expected;
            assert!(
                error < 1e-14,
                "direction {d}: total mass {total:.6e} != rho h A {expected:.6e} (relative error {error:.3e})"
            );
        }
    }

    #[test]
    fn test_me_global_matches_the_exact_bilinear_coefficients() {
        let pre = pre_from(&FLAT_SQUARE);
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let mass = rho * pre.thickness * element_area(&pre);

        // Bilinear consistent mass on a rectangle: M_ii = m/9, adjacent
        // M_ij = m/18 and opposite M_ij = m/36. The 2x2 quadrature is exact for
        // these products, so the tolerance is round-off, not discretisation.
        let coefficient = |i: i64, j: i64| match (i - j).abs() {
            0 => 4.0 / 36.0,
            1 | 3 => 2.0 / 36.0,
            _ => 1.0 / 36.0,
        };

        for d in 0..3 {
            for i in 0..4_i64 {
                for j in 0..4_i64 {
                    let expected = coefficient(i, j) * mass;
                    let actual = m[(6 * i as usize + d, 6 * j as usize + d)];
                    let error = (actual - expected).abs() / mass;
                    assert!(
                        error < 1e-14,
                        "M[{i},{j}] direction {d}: {actual:.6e} != {expected:.6e} (relative error {error:.3e})"
                    );
                }
            }
        }
    }

    #[test]
    fn test_me_global_rotary_inertia_is_rho_h3_a_over_12() {
        let pre = pre_from(&FLAT_SQUARE);
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let expected = rho * pre.thickness.powi(3) / 12.0 * element_area(&pre);

        for (k, total) in rotary_mass_per_direction(&m).iter().enumerate() {
            let error = (total - expected).abs() / expected;
            assert!(
                error < 1e-14,
                "rotation direction {k}: rotary mass {total:.6e} != rho h^3/12 A {expected:.6e} (relative error {error:.3e})"
            );
        }
    }

    #[test]
    fn test_me_composite_global_matches_the_rho_h_construction() {
        let pre = pre_from(&FLAT_SQUARE);
        let rho = 7800.0;
        let m = compute_me_global(&pre, rho);
        let m_composite = compute_me_composite_global(
            &pre,
            rho * pre.thickness,
            rho * pre.thickness.powi(3) / 12.0,
        );
        assert!(
            max_abs_diff(&m, &m_composite) <= 1e-12 * max_abs(&m),
            "the composite mass path must reproduce the rho h / rho h^3/12 construction"
        );
    }

    // ========================================================================
    // WU6 - Tier 1a: the 2017 paper's own basic tests
    //
    // Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
    // Structures 182:404-418, Section 4 ("Basic numerical tests"), pp. 410-411:
    //
    //   Isotropy    - "the element passes the test of spatial isotropy": the
    //                 stiffness is invariant under the element's orientation and
    //                 the node-numbering sequence.
    //   Zero energy - "the number of zero eigenvalues of the stiffness matrix of
    //                 a single unsupported element are counted ... For the new
    //                 element only the six zero eigenvalues corresponding to the
    //                 six rigid body modes are obtained."
    //   Patch tests - "We perform three patch tests: the membrane, bending and
    //                 shearing patch tests ... The patch of elements is subjected
    //                 to the minimum number of constraints to prevent rigid body
    //                 motions and the nodal point forces on the boundary
    //                 corresponding to the constant stress states are applied.
    //                 The patch tests are passed if the correct values of
    //                 constant stress fields are calculated at any location
    //                 within the mesh."
    //
    // The Fig. 5 mesh is `STAR_NODES`/`STAR_ELEMS` (WU1: Ko, Lee & Bathe (2017),
    // C&S 182:404-418, Fig. 5, p. 407, = Ko, Bathe & Zhang (2025),
    // C&S 308:107622, Fig. 7(a)); `BC_2017_PATCH` is the minimum constraint set.
    // The load magnitudes are NOT published: design Section 4.2 derives them
    // element-independently from the constant state, which is what these
    // fixtures do. No load below is built from the element stiffness, so a
    // wrong element cannot pass by construction.
    // ========================================================================

    /// Eigenvalues of the symmetrised 24x24 stiffness, sorted ascending.
    fn sorted_eigenvalues(k: &Mat24) -> Vec<f64> {
        let mut sym = *k;
        for i in 0..24 {
            for j in 0..24 {
                sym[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
            }
        }
        let mut eig: Vec<f64> = sym.symmetric_eigenvalues().iter().cloned().collect();
        eig.sort_by(|a, b| a.total_cmp(b));
        eig
    }

    /// Rotation matrix by `angle` about `axis` (Rodrigues).
    fn rotation_matrix(axis: [f64; 3], angle: f64) -> Matrix3<f64> {
        let a = Vector3::new(axis[0], axis[1], axis[2]).normalize();
        let (s, c) = angle.sin_cos();
        let k = Matrix3::new(0.0, -a[2], a[1], a[2], 0.0, -a[0], -a[1], a[0], 0.0);
        Matrix3::identity() + s * k + (1.0 - c) * (k * k)
    }

    /// Rigidly rotate a four-node element geometry.
    fn rotate_geom(c: &[[f64; 3]; 4], r: &Matrix3<f64>) -> [[f64; 3]; 4] {
        let mut out = [[0.0f64; 3]; 4];
        for i in 0..4 {
            let v = r * Vector3::new(c[i][0], c[i][1], c[i][2]);
            out[i] = [v[0], v[1], v[2]];
        }
        out
    }

    /// Co-rotate a 24-DOF field with `R`: both the translational and the
    /// rotational nodal triples are rotated (the DOFs co-rotate with the
    /// element, Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (3)).
    fn rotate_dofs(u: &[f64; 24], r: &Matrix3<f64>) -> [f64; 24] {
        let mut out = [0.0f64; 24];
        for i in 0..4 {
            for k in 0..3 {
                for b in 0..3 {
                    out[6 * i + k] += r[(k, b)] * u[6 * i + b];
                    out[6 * i + 3 + k] += r[(k, b)] * u[6 * i + 3 + b];
                }
            }
        }
        out
    }

    /// Rotate an engineering 3-vector `[v11, v22, 2 v12]` from the element's
    /// local frame to the global frame (`t3` maps global -> local).
    fn rotate_eng3_to_global(t3: &Matrix3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
        let t = Matrix3::new(v[0], 0.5 * v[2], 0.0, 0.5 * v[2], v[1], 0.0, 0.0, 0.0, 0.0);
        let r = t3.transpose();
        let tg = r * t * r.transpose();
        Vector3::new(tg[(0, 0)], tg[(1, 1)], 2.0 * tg[(0, 1)])
    }

    /// Rotate the stress-like 3-vector `[s11, s22, s12]` from the element's local
    /// frame to the global frame. Unlike an engineering strain, the third
    /// component is the tensor shear `s12` (not `2 s12`).
    fn rotate_stress_to_global(t3: &Matrix3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
        let t = Matrix3::new(v[0], v[2], 0.0, v[2], v[1], 0.0, 0.0, 0.0, 0.0);
        let r = t3.transpose();
        let tg = r * t * r.transpose();
        Vector3::new(tg[(0, 0)], tg[(1, 1)], tg[(0, 1)])
    }

    /// Rotate the transverse-shear 2-vector `[gamma_13, gamma_23]` to the
    /// global frame (`t3` maps global -> local).
    fn rotate_shear_to_global(t3: &Matrix3<f64>, g: &[f64; 2]) -> [f64; 2] {
        let t = Matrix3::new(
            0.0,
            0.0,
            0.5 * g[0],
            0.0,
            0.0,
            0.5 * g[1],
            0.5 * g[0],
            0.5 * g[1],
            0.0,
        );
        let r = t3.transpose();
        let tg = r * t * r.transpose();
        [2.0 * tg[(0, 2)], 2.0 * tg[(1, 2)]]
    }

    // ------------------------------------------------------------------
    // 7.1 - the rigid-body fixture and the zero-energy test (T1.2)
    // ------------------------------------------------------------------

    /// The drill-rotation operator `B~` (3x4) of Ko, Bathe & Zhang (2025),
    /// C&S 308:107622, Eq. (19a)/(19b), pp. 10 and 12, recovered from the
    /// production Eq. (18) operator. The production operator stores
    /// `B_md[:, 6 i + 3 + beta] = B~[:, i] vd_local[beta]` with `vd_local` the
    /// unit drill vector `V^D` in the element's local frame, so projecting a
    /// node's rotation triple onto `V^D` recovers `B~[:, i]` exactly.
    fn drill_operator_eq19(pre: &Mitc4PlusDPrecomputed, r: f64, s: f64) -> SMatrix<f64, 3, 4> {
        let bmd = b_drill_membrane_2025(pre, r, s);
        let vd = Vector3::new(
            pre.v_d.dot(&pre.e1),
            pre.v_d.dot(&pre.e2),
            pre.v_d.dot(&pre.e3),
        );
        let mut bt = SMatrix::<f64, 3, 4>::zeros();
        for i in 0..4 {
            for beta in 0..3 {
                for row in 0..3 {
                    bt[(row, i)] += bmd[(row, 6 * i + 3 + beta)] * vd[beta];
                }
            }
        }
        bt
    }

    /// Rebuild `B~` (3x4) directly from the four edge terms of Eq. (19b), used
    /// to show the recovered operator is a genuine Eq. (19b) object and not an
    /// artefact of the projection: it must equal `drill_operator_eq19`.
    fn drill_operator_eq19_edges(
        pre: &Mitc4PlusDPrecomputed,
        r: f64,
        s: f64,
    ) -> SMatrix<f64, 3, 4> {
        let ratio = drill_jacobian_ratio(pre, r, s);
        let dh = drill_midside_shape_derivatives(r, s);
        let mut bt = SMatrix::<f64, 3, 4>::zeros();
        for e in 0..4 {
            let edge = pre.drill_edges[e];
            let (h_r, h_s) = dh[e];
            let coeff_rr = ratio * h_r * edge.c_r;
            let coeff_ss = -ratio * h_s * edge.c_s;
            let coeff_rs = 0.5 * ratio * (h_s * edge.c_r - h_r * edge.c_s);
            bt[(0, edge.start)] -= coeff_rr;
            bt[(0, edge.end)] += coeff_rr;
            bt[(1, edge.start)] -= coeff_ss;
            bt[(1, edge.end)] += coeff_ss;
            bt[(2, edge.start)] -= 2.0 * coeff_rs;
            bt[(2, edge.end)] += 2.0 * coeff_rs;
        }
        covariant_to_local_mapping(&j_loc_at(pre, 0.0, 0.0)) * bt
    }

    /// The drill block `B~^T C B~` (4x4) of Eq. (19a)/(19b) with the 2x2
    /// surface rule of B Section 3.1, p. 13, over the four corner drill
    /// rotations. `mode = 0` is the production operator; `mode = 1` makes the
    /// block inert (a zero operator, null dimension 4); `mode = 2` is a
    /// deliberate synthetic rank gain (a positive multiple of the identity
    /// added, null dimension 0). Modes 1 and 2 are perturbations used only to
    /// show the null-space assertions are load-bearing, not vacuous; they are
    /// not paper variants and no production path is affected.
    fn drill_block_eq19(pre: &Mitc4PlusDPrecomputed, mode: u8) -> SMatrix<f64, 4, 4> {
        let cm = &pre.constitutive.cm;
        let mut m = SMatrix::<f64, 4, 4>::zeros();
        for g in 0..N_GAUSS {
            let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
            let sqrt_g = surface_measure(pre, r, s);
            let bt = match mode {
                0 | 2 => drill_operator_eq19(pre, r, s),
                1 => SMatrix::<f64, 3, 4>::zeros(),
                _ => unreachable!(),
            };
            m += (bt.transpose() * cm * bt) * (GAUSS_W[g] * sqrt_g);
        }
        if mode == 2 {
            let scale = m.iter().fold(0.0f64, |a, &v| a.max(v.abs())).max(1e-30);
            for i in 0..4 {
                m[(i, i)] += scale;
            }
        }
        let mut sym = m;
        for i in 0..4 {
            for j in 0..4 {
                sym[(i, j)] = 0.5 * (m[(i, j)] + m[(j, i)]);
            }
        }
        sym
    }

    /// Sorted eigenvalues and null space of the 4x4 drill block, with the
    /// spec's threshold `|lambda| <= 1e-10 lambda_max`.
    fn drill_block_null_space(pre: &Mitc4PlusDPrecomputed, mode: u8) -> (Vec<f64>, Vec<[f64; 4]>) {
        let m = drill_block_eq19(pre, mode);
        let eig = m.symmetric_eigen();
        let mut idx: Vec<usize> = (0..4).collect();
        idx.sort_by(|&a, &b| eig.eigenvalues[a].total_cmp(&eig.eigenvalues[b]));
        let lam: Vec<f64> = idx.iter().map(|&i| eig.eigenvalues[i]).collect();
        let lam_max = lam.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        let mut null = Vec::new();
        for &i in idx.iter() {
            if eig.eigenvalues[i].abs() <= 1e-10 * lam_max.max(1e-30) {
                let mut v = [0.0f64; 4];
                for row in 0..4 {
                    v[row] = eig.eigenvectors[(row, i)];
                }
                null.push(v);
            }
        }
        (lam, null)
    }

    /// Eigenvalues (ascending) of the symmetrised stiffness with the `fixed`
    /// DOFs removed. The only constraint used here is the paper's own drill
    /// device: `theta_z` free at every node except corner B (code node 3) --
    /// Ko, Bathe & Zhang (2025), C&S 308:107622, Fig. 7(b)(c)(d), B Section 3.1
    /// pp. 13-14.
    fn reduced_eigenvalues(k: &Mat24, fixed: &[usize]) -> Vec<f64> {
        let mut is_fixed = [false; 24];
        for &d in fixed {
            is_fixed[d] = true;
        }
        let keep: Vec<usize> = (0..24).filter(|&i| !is_fixed[i]).collect();
        let n = keep.len();
        let mut red = DMatrix::zeros(n, n);
        for (a, &i) in keep.iter().enumerate() {
            for (b, &j) in keep.iter().enumerate() {
                red[(a, b)] = k[(i, j)];
            }
        }
        let sym = (&red + red.transpose()) * 0.5;
        let mut eig: Vec<f64> = sym.symmetric_eigenvalues().iter().cloned().collect();
        eig.sort_by(|a, b| a.total_cmp(b));
        eig
    }

    fn max_abs_34(a: &SMatrix<f64, 3, 4>, b: &SMatrix<f64, 3, 4>) -> f64 {
        let mut m = 0.0f64;
        for i in 0..3 {
            for j in 0..4 {
                m = m.max((a[(i, j)] - b[(i, j)]).abs());
            }
        }
        m
    }

    #[test]
    fn test_t1a_zero_energy_modes_single_unsupported_element_six_or_seven() {
        // Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
        // Structures 182:404-418, Section 4 ("Basic numerical tests"), p. 410:
        // "In the zero energy mode test, the number of zero eigenvalues of the
        //  stiffness matrix of a single unsupported element are counted ...
        //  For the new element only the six zero eigenvalues corresponding to
        //  the six rigid body modes are obtained.  That is, the element passes
        //  the zero energy mode test."
        //
        // SPEC REV 6 (Requirement 5). The count is taken with the drill DOF
        // constrained the way paper B's own patch tests constrain it --
        // `theta_z` free at every node except the corner node B, i.e. code
        // node 3 (Ko, Bathe & Zhang (2025), C&S 308:107622, Fig. 7(b)(c)(d),
        // B Section 3.1 pp. 13-14) -- and the drill operator's own null space
        // is characterized separately. Eq. (19b)'s columns are DIFFERENCES of
        // edge terms, so a constant drill rotation `theta^D` is a zero-energy
        // direction of the paper's own equations and is not a rigid-body mode:
        // it is counted here as the drill null space, not as a rigid-body mode.
        // No penalty, constraint or numerical factor is added.
        //
        // The count is exactly 7 on the flat rectangle (the six rigid-body
        // modes plus the single-corner `theta_z` condition leaving one of the
        // flat rectangle's two curl-induced drill directions) and exactly 6 on
        // the flat distorted, ruled-warped and doubly-warped quads. The flat
        // rectangle's surplus is the Eq. (19d) curl-induced `theta_z` hourglass,
        // not a spurious 2017-core mode.
        //
        // The rigid-body representation is verified in TWO parts. With `theta_z`
        // free, all six rigid-body fields `u = t + omega x (x - x_c)`,
        // `theta = omega` (constant), satisfy the `1e-12 lambda_max` residual
        // bound. The paper's single-corner `theta_z` condition is one linear
        // condition; besides removing one of the two flat-rectangle drill
        // directions it also removes the `omega_z` rigid-body rotation, whose
        // `theta_z` is constant (`omega_z` . `V^D`), so with the constraint
        // applied only the FIVE fields satisfying `theta_z(B) = 0` can satisfy
        // the bound. The element is linear, so the constrained residual of a
        // field that already satisfies `theta_z(B) = 0` is its free residual
        // restricted to the free rows.
        //
        // The six rigid-body fields' drill images ARE the constant drill
        // rotation -- itself a null direction of Eq. (19b) by the telescoping
        // columns -- so the test MUST NOT assert that they are not all
        // annihilated by the drill block. The satisfiable non-vacuity
        // distinction is the drill block's own null space: every null vector is
        // a pure drill-rotation field (zero translations, zero alpha/beta) at a
        // measurable distance from the rigid-body space, and the block is live.
        //
        // The six rigid-body fields are the design 4.3 construction
        // `u_i = t + omega x (x_i - x_c)`, `theta_i = omega` (constant), in the
        // element's 6-DOF layout (Ko, Lee & Bathe (2017), C&S 182:404-418,
        // Eq. (3)); the element-local form is `rigid_body_fields`.
        let geoms = [
            // (name, geometry, expected constrained zero-count, expected drill
            //  null-space dimension)
            ("flat-rectangle", RECT, 7usize, 2usize),
            ("flat-distorted", FLAT_DISTORTED, 6, 1),
            ("ruled-warped", RULED_WARPED, 6, 1),
            ("doubly-warped", DOUBLY_WARPED, 6, 1),
        ];

        // The paper's own constraint: `theta_z` free except at corner node B
        // (code node 3), i.e. slot 6*3+5 = 23.
        let drill_fixed = [6 * 3 + 5];
        // Collect all measured failures so every geometry is reported, not just
        // the first one to fail.
        let mut failures: Vec<String> = Vec::new();

        for (name, c, expected_count, expected_null) in geoms {
            let pre = pre_from(&c);
            let k = compute_ke_local(&pre);
            let lam_max = lambda_max(&k);

            // (1) The zero-eigenvalue count with the paper's drill constraint.
            let lam_c = reduced_eigenvalues(&k, &drill_fixed);
            let lam_c_max = lam_c.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
            let zeros = lam_c
                .iter()
                .filter(|v| v.abs() <= 1e-10 * lam_c_max)
                .count();
            let sep = lam_c
                .iter()
                .map(|v| v.abs())
                .nth(6)
                .expect("the constrained 23x23 system has at least seven eigenvalues");

            // (2) The six physical rigid-body fields must be annihilated with
            // `theta_z` free; with the paper's constraint applied, the five
            // fields satisfying `theta_z(B) = 0` must satisfy the same bound.
            let mut worst_rb = 0.0f64;
            let mut worst_field = 0usize;
            let mut worst_rb_constrained = 0.0f64;
            let mut worst_field_constrained = 0usize;
            let mut rb_surviving = 0usize;
            for (i, u) in rigid_body_fields(&pre).iter().enumerate() {
                let ku = k * Vec24::from_column_slice(u);
                let inf = ku.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                let un = u.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                let ratio = inf / (lam_max * un);
                if ratio > worst_rb {
                    worst_rb = ratio;
                    worst_field = i;
                }
                // The fields satisfying the paper's own drill constraint
                // `theta_z(node B) = 0`: apply the constraint by zeroing slot
                // 6*3+5 and restricting the residual to the free rows.
                if u[6 * 3 + 5].abs() <= 1e-12 {
                    rb_surviving += 1;
                    let mut uc = *u;
                    uc[6 * 3 + 5] = 0.0;
                    let kuc = k * Vec24::from_column_slice(&uc);
                    let inf_free = kuc
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| *j != 6 * 3 + 5)
                        .fold(0.0f64, |m, (_, &v)| m.max(v.abs()));
                    let unc = uc.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                    let ratio_c = inf_free / (lam_c_max * unc);
                    if ratio_c > worst_rb_constrained {
                        worst_rb_constrained = ratio_c;
                        worst_field_constrained = i;
                    }
                }
            }

            // (3) The drill block's own null space, Eq. (19a)/(19b).
            let (lam4, null) = drill_block_null_space(&pre, 0);
            let kd = drill_ke_local(&pre);
            let vd = Vector3::new(
                pre.v_d.dot(&pre.e1),
                pre.v_d.dot(&pre.e2),
                pre.v_d.dot(&pre.e3),
            );
            // Every null vector must be a PURE drill-rotation field: lift it to
            // the 24-DOF layout as `theta_i = n_i V^D` (zero translations, zero
            // alpha/beta), it must be annihilated by the drill block, and it
            // must not be a rigid-body field (a rigid rotation carries
            // `u_i = omega x x_i != 0`), so the drill null space is the drill
            // block's own and not the rigid-body space.
            let rb = rigid_body_fields(&pre);
            let mut worst_null_energy = 0.0f64;
            let mut worst_null_rb_distance = 0.0f64;
            for n in null.iter() {
                let mut u = [0.0f64; 24];
                for i in 0..4 {
                    for beta in 0..3 {
                        u[6 * i + 3 + beta] = n[i] * vd[beta];
                    }
                }
                let q = quadratic(&u, &kd).abs();
                worst_null_energy = worst_null_energy.max(q / (lam_max * norm2(&u)));
                // Distance from the six-dimensional rigid-body space.
                let mut g = DMatrix::<f64>::zeros(6, 6);
                let mut rhs = DVector::<f64>::zeros(6);
                for a in 0..6 {
                    for b in 0..6 {
                        g[(a, b)] = rb[a].iter().zip(rb[b].iter()).map(|(x, y)| x * y).sum();
                    }
                    rhs[a] = rb[a].iter().zip(u.iter()).map(|(x, y)| x * y).sum();
                }
                let coef: DVector<f64> = g
                    .lu()
                    .solve(&rhs)
                    .expect("the six rigid-body fields are linearly independent");
                let mut resid = 0.0f64;
                for a in 0..6 {
                    for j in 0..24 {
                        resid += (u[j] - coef[a] * rb[a][j]).powi(2);
                    }
                }
                worst_null_rb_distance =
                    worst_null_rb_distance.max(resid.sqrt() / norm2(&u).sqrt());
            }

            // (4) Non-vacuity. The drill block must be live: its null-space
            // dimension must be the paper's own value, an inert block must
            // give 4, and a rank gain must drop it below the paper's value.
            // The six rigid-body fields' drill images are the CONSTANT drill
            // rotation, which Eq. (19b)'s telescoping columns make a null
            // direction of the paper's own operator; what distinguishes the
            // rigid-body space from the drill null space is that every drill
            // null vector is a pure drill rotation (zero translations, zero
            // alpha/beta) and is not a rigid-body field.
            let (_lam4_inert, null_inert) = drill_block_null_space(&pre, 1);
            let (_lam4_rank, null_rank) = drill_block_null_space(&pre, 2);
            let mut worst_rb_drill = 0.0f64;
            for u in rigid_body_fields(&pre).iter() {
                let ku = kd * Vec24::from_column_slice(u);
                let inf = ku.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                let un = u.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
                worst_rb_drill = worst_rb_drill.max(inf / (lam_max * un));
            }

            println!(
                "T1.2 rev6 {name}: lambda_max={lam_max:.6e}; constrained zero-count={zeros} \
                 (expected {expected_count}); (|lambda_1..7|={:?}); |lambda_7|={sep:.6e} \
                 ({lam_c_max:.3e} lambda_max); \
                 worst ||K u_rb||_inf/(lambda_max ||u_rb||_inf)={worst_rb:.3e} (field {worst_field}); \
                 with constraint applied: {worst_rb_constrained:.3e} (field {worst_field_constrained}); \
                 rigid-body fields satisfying theta_z(B)=0: {rb_surviving}/6; \
                 drill block |lambda_1..4|={lam4:?} null-dim={} (expected {expected_null}); \
                 inert null-dim={} rank-gain null-dim={}; \
                 worst drill-null energy/(lambda_max ||u||^2)={worst_null_energy:.3e}; \
                 worst drill-null distance-from-rigid-body-space={worst_null_rb_distance:.3e}; \
                 worst ||K_drill u_rb||_inf/(lambda_max ||u_rb||_inf)={worst_rb_drill:.3e}",
                &lam_c.iter().map(|v| v.abs()).take(7).collect::<Vec<_>>(),
                null.len(),
                null_inert.len(),
                null_rank.len(),
            );

            // (1) the per-geometry zero-eigenvalue count under the paper's own
            // drill constraint: exactly 7 on the flat rectangle, 6 elsewhere.
            if zeros != expected_count {
                failures.push(format!(
                    "{name}: expected exactly {expected_count} zero eigenvalues with the drill \
                     DOF constrained as Ko, Bathe & Zhang (2025), C&S 308:107622, Fig. 7 \
                     constrains it, got {zeros} (lambda_max {lam_c_max:.3e})"
                ));
            }
            // The count-6 geometries must show an unambiguous gap: the seventh
            // eigenvalue at least 1e-9 lambda_max, strictly above the 1e-10
            // count threshold. The flat rectangle's seventh eigenvalue is ~0
            // (it has seven zero modes), so no separation is asserted there.
            if expected_count == 6 && sep < 1e-9 * lam_c_max {
                failures.push(format!(
                    "{name}: the seventh eigenvalue {sep:.3e} is not separated from the six \
                     zero modes by 1e-9 lambda_max ({:.3e})",
                    1e-9 * lam_c_max
                ));
            }
            // (2) rigid-body verification, two parts: the six fields with
            // `theta_z` free, and the five fields with `theta_z(B) = 0` under
            // the paper's constraint.
            if worst_rb > 1e-12 {
                failures.push(format!(
                    "{name}: ||K u_rb||_inf/(lambda_max ||u_rb||_inf) = {worst_rb:.3e} \
                     (field {worst_field})"
                ));
            }
            if rb_surviving != 5 {
                failures.push(format!(
                    "{name}: expected the paper's single-corner `theta_z(B) = 0` condition to \
                     leave exactly five rigid-body fields, got {rb_surviving}"
                ));
            }
            if worst_rb_constrained > 1e-12 {
                failures.push(format!(
                    "{name}: with the constraint applied, \
                     ||K u_rb||_inf/(lambda_max ||u_rb||_inf) = {worst_rb_constrained:.3e} \
                     (field {worst_field_constrained})"
                ));
            }
            // (3) the drill block's own null space, dimension and purity.
            if null.len() != expected_null {
                failures.push(format!(
                    "{name}: the drill block B~^T C B~ null-space dimension is {}, expected \
                     {expected_null} (constant drill rotation{}); lambda_4 = {lam4:?}",
                    null.len(),
                    if expected_null == 2 {
                        " + Eq. (19d) curl-induced theta_z hourglass"
                    } else {
                        ""
                    }
                ));
            }
            if worst_null_energy > 1e-12 {
                failures.push(format!(
                    "{name}: a drill null vector carries energy {worst_null_energy:.3e} \
                     (it must be a pure drill-rotation field)"
                ));
            }
            if worst_null_rb_distance <= 2.4 {
                failures.push(format!(
                    "{name}: a drill null vector lies in the rigid-body space \
                     (distance {worst_null_rb_distance:.3e}); the drill null space must be its own"
                ));
            }
            // (4) non-vacuity: the block is live (null dim < 4), an inert block
            // would give 4, and a rank gain drops the dimension below the
            // paper's value.
            if null.len() >= 4 {
                failures.push(format!(
                    "{name}: the drill block is inert (null-space dimension {})",
                    null.len()
                ));
            }
            if null_inert.len() != 4 {
                failures.push(format!(
                    "{name}: an inert drill block must have null-space dimension 4, got {} \
                     (the null-space assertion must catch an inert block)",
                    null_inert.len()
                ));
            }
            if null_rank.len() != 0 {
                failures.push(format!(
                    "{name}: a rank gain in the drill block must drop its null-space dimension \
                     to 0, got {} (the assertion must catch a rank gain)",
                    null_rank.len()
                ));
            }
        }

        assert!(
            failures.is_empty(),
            "zero-energy-mode test failures:\n{}",
            failures.join("\n")
        );

        // The rank-change perturbation is a genuine reconstruction of Eq. (19b):
        // with the curl zeros restored it equals the production operator.
        for c in [RECT, FLAT_DISTORTED, RULED_WARPED, DOUBLY_WARPED] {
            let pre = pre_from(&c);
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let a = drill_operator_eq19(&pre, r, s);
                let b = drill_operator_eq19_edges(&pre, r, s);
                let scale = a.iter().fold(0.0f64, |m, &v| m.max(v.abs())).max(1e-30);
                assert!(
                    max_abs_34(&a, &b) <= 1e-12 * scale,
                    "the Eq. (19b) edge reconstruction must equal the production operator \
                     (max diff {:.3e}, scale {scale:.3e})",
                    max_abs_34(&a, &b)
                );
            }
        }
    }

    // ------------------------------------------------------------------
    // 7.2 - isotropy of the 2017 core (T1.1)
    // ------------------------------------------------------------------

    #[test]
    fn test_t1a_isotropy_element_orientation_and_node_sequence_invariant() {
        // Ko, Lee & Bathe (2017), C&S 182:404-418, Section 4, p. 410:
        // "The element behavior should not depend on the ... sequence of node
        //  numbering, i.e. on the element orientation ... The element passes
        //  the test of spatial isotropy."
        //
        // The rotated element's stiffness is formed in global coordinates so
        // that the rotations act on the co-rotated DOF, and the symmetrised
        // eigenvalues are compared to the reference; the node-sequence variants
        // rebuild the element with a permuted connectivity and must reproduce
        // the exactly-permuted entries `P K P^T` of the reference matrix.
        let axis = [1.0, 2.0, 3.0];
        let angles = [
            0.0,
            0.7,
            core::f64::consts::FRAC_PI_4,
            core::f64::consts::FRAC_PI_2,
        ];
        let geoms = [
            ("flat-square", FLAT_SQUARE),
            ("flat-distorted", FLAT_DISTORTED),
            ("ruled-warped", RULED_WARPED),
        ];
        let u_ref = generic_dofs();

        for (name, c) in geoms {
            let k0 = compute_ke_global(&pre_from(&c));
            let lam0 = sorted_eigenvalues(&k0);
            let lam_max0 = lam0.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
            let q0 = quadratic(&u_ref, &k0);

            let mut worst_eig = 0.0f64;
            let mut worst_energy = 0.0f64;
            for &theta in angles.iter() {
                let r = rotation_matrix(axis, theta);
                let k = compute_ke_global(&pre_from(&rotate_geom(&c, &r)));
                let lam = sorted_eigenvalues(&k);
                for i in 0..24 {
                    worst_eig = worst_eig.max((lam[i] - lam0[i]).abs());
                }
                let u_rot = rotate_dofs(&u_ref, &r);
                let q = quadratic(&u_rot, &k);
                worst_energy = worst_energy.max((q - q0).abs() / q0.abs().max(1e-30));
            }
            println!(
                "T1.1 {name}: worst |dlambda|={worst_eig:.3e} (1e-10 lambda_max={:.3e}); \
                 worst |duKu|/|uKu|={worst_energy:.3e}",
                1e-10 * lam_max0
            );
            assert!(
                worst_eig <= 1e-10 * lam_max0,
                "{name}: orientation changes an eigenvalue by {worst_eig:.3e} > 1e-10 lambda_max ({:.3e})",
                1e-10 * lam_max0
            );
            assert!(
                worst_energy <= 1e-10,
                "{name}: the co-rotated strain energy is not invariant ({worst_energy:.3e} relative)"
            );

            // Node-numbering sequences (Ko, Lee & Bathe (2017), C&S
            // 182:404-418, p. 410). `seq[k]` is the reference vertex placed at
            // the new element's local slot `k`.
            let seqs = [
                [0usize, 1, 2, 3],
                [1, 2, 3, 0],
                [2, 3, 0, 1],
                [3, 0, 1, 2],
                [0, 3, 2, 1],
            ];
            let kmax = max_abs(&k0);
            for seq in seqs {
                let mut cp = [[0.0f64; 3]; 4];
                for kk in 0..4 {
                    cp[kk] = c[seq[kk]];
                }
                let kp = compute_ke_global(&pre_from(&cp));
                let mut expected = Mat24::zeros();
                for a in 0..4 {
                    for b in 0..4 {
                        for ka in 0..6 {
                            for kb in 0..6 {
                                expected[(6 * a + ka, 6 * b + kb)] =
                                    k0[(6 * seq[a] + ka, 6 * seq[b] + kb)];
                            }
                        }
                    }
                }
                let d = max_abs_diff(&kp, &expected);
                let lam_new = sorted_eigenvalues(&kp);
                let dlam = (0..24)
                    .map(|i| (lam_new[i] - lam0[i]).abs())
                    .fold(0.0f64, f64::max);
                // Physical (mapping-independent) check: the same physical field
                // expressed in the two node orders must carry the same energy.
                let mut u_perm = [0.0f64; 24];
                for a in 0..4 {
                    for kk in 0..6 {
                        u_perm[6 * a + kk] = u_ref[6 * seq[a] + kk];
                    }
                }
                let q_perm = quadratic(&u_perm, &kp);
                let q_rel = (q_perm - q0).abs() / q0.abs().max(1e-30);
                println!(
                    "T1.1 {name}: node-sequence {seq:?}: max|K_new - P K_ref P^T|={d:.3e} ({:.3e} max|K|); \
                     eigenvalue delta={dlam:.3e} ({:.3e} lambda_max); u^T K u rel={q_rel:.3e}",
                    d / kmax,
                    dlam / lam_max0
                );
                assert!(
                    d <= 1e-12 * kmax,
                    "{name}: node sequence {seq:?} gives {d:.3e} > 1e-12 max|K| ({:.3e})",
                    1e-12 * kmax
                );
            }
        }
    }

    // ------------------------------------------------------------------
    // 7.3/7.4 - the patch tests: the Fig. 5 mesh, the minimum constraints and
    // the boundary nodal loads derived element-independently from the constant
    // state (design 4.2).
    // ------------------------------------------------------------------

    /// The star patch's outer boundary, counter-clockwise, as `(start, end)`
    /// patch-node pairs; the outward normal of `(p -> q)` is `(dy, -dx)/L`.
    const STAR_BOUNDARY: [(usize, usize); 4] = [(0, 1), (1, 2), (2, 3), (3, 0)];

    /// Node coordinates of star-patch element `e` as a `[[f64; 3]; 4]` list.
    fn star_elem_coords(nodes: &[[f64; 3]; 8], e: usize) -> [[f64; 3]; 4] {
        let el = STAR_ELEMS[e];
        [nodes[el[0]], nodes[el[1]], nodes[el[2]], nodes[el[3]]]
    }

    /// The five element `pre` values of the star patch.
    fn star_patch_pres(nodes: &[[f64; 3]; 8]) -> [Mitc4PlusDPrecomputed; 5] {
        [
            pre_from(&star_elem_coords(nodes, 0)),
            pre_from(&star_elem_coords(nodes, 1)),
            pre_from(&star_elem_coords(nodes, 2)),
            pre_from(&star_elem_coords(nodes, 3)),
            pre_from(&star_elem_coords(nodes, 4)),
        ]
    }

    /// Assemble the 48-DOF star-patch stiffness (task 2.4, design 4.2 step 3).
    fn assemble_star_patch(pres: &[Mitc4PlusDPrecomputed; 5]) -> DMatrix<f64> {
        let mut k = DMatrix::zeros(48, 48);
        for (e, pre) in pres.iter().enumerate() {
            let ke = compute_ke_global(pre);
            let el = STAR_ELEMS[e];
            for a in 0..4 {
                for b in 0..4 {
                    for ka in 0..6 {
                        for kb in 0..6 {
                            k[(6 * el[a] + ka, 6 * el[b] + kb)] += ke[(6 * a + ka, 6 * b + kb)];
                        }
                    }
                }
            }
        }
        k
    }

    /// Integrate a distributed load `per_length(point, outward_normal) ->
    /// [f64; 6]` over the patch's outer boundary with 2-point Gauss per edge
    /// (design 4.2 step 2).
    fn boundary_integrate<F>(nodes: &[[f64; 3]; 8], mut per_length: F) -> DVector<f64>
    where
        F: FnMut([f64; 3], [f64; 2]) -> [f64; 6],
    {
        let mut f = DVector::zeros(48);
        for &(s, e) in STAR_BOUNDARY.iter() {
            let p = nodes[s];
            let q = nodes[e];
            let d = [q[0] - p[0], q[1] - p[1]];
            let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let n = [d[1] / l, -d[0] / l];
            for xi in [(-GP, 1.0f64), (GP, 1.0)] {
                let x = [
                    p[0] + 0.5 * (1.0 + xi.0) * d[0],
                    p[1] + 0.5 * (1.0 + xi.0) * d[1],
                    0.0,
                ];
                let load = per_length(x, n);
                let ns = 0.5 * (1.0 - xi.0);
                let ne = 0.5 * (1.0 + xi.0);
                for kk in 0..6 {
                    f[6 * s + kk] += xi.1 * (l * 0.5) * ns * load[kk];
                    f[6 * e + kk] += xi.1 * (l * 0.5) * ne * load[kk];
                }
            }
        }
        f
    }

    /// Solve the constrained patch system `K_ff u_f = f_f - K_fc u_c` with a
    /// dense LU (design 4.2 step 3).
    fn solve_constrained(
        k: &DMatrix<f64>,
        f: &DVector<f64>,
        prescribed: &[(usize, f64)],
    ) -> DVector<f64> {
        let n = k.nrows();
        let mut fixed = vec![false; n];
        let mut u = DVector::zeros(n);
        for &(d, v) in prescribed {
            assert!(!fixed[d], "duplicate constrained dof {d}");
            fixed[d] = true;
            u[d] = v;
        }
        let free: Vec<usize> = (0..n).filter(|&i| !fixed[i]).collect();
        let m = free.len();
        let mut kff = DMatrix::zeros(m, m);
        let mut rhs = DVector::zeros(m);
        for (a, &i) in free.iter().enumerate() {
            rhs[a] = f[i];
            for (b, &j) in free.iter().enumerate() {
                kff[(a, b)] = k[(i, j)];
            }
            for &(d, v) in prescribed.iter() {
                rhs[a] -= k[(i, d)] * v;
            }
        }
        let sol = kff
            .lu()
            .solve(&rhs)
            .expect("star-patch constrained system must be non-singular");
        for (a, &i) in free.iter().enumerate() {
            u[i] = sol[a];
        }
        u
    }

    /// The element displacement (global, 24) of star-patch element `e` from a
    /// 48-DOF patch vector.
    fn star_elem_disp(u: &DVector<f64>, e: usize) -> Vec24 {
        let el = STAR_ELEMS[e];
        let mut ue = Vec24::zeros();
        for a in 0..4 {
            for kk in 0..6 {
                ue[6 * a + kk] = u[6 * el[a] + kk];
            }
        }
        ue
    }

    /// Plane-stress membrane strain `eps = cm_raw^{-1} sigma` (engineering
    /// shear) for the isotropic shell constitutive `shell_iso`.
    fn membrane_strain(sigma: &Vector3<f64>) -> Vector3<f64> {
        shell_iso()
            .cm_raw
            .try_inverse()
            .expect("cm_raw must be invertible")
            * sigma
    }

    /// Exact in-plane displacement of a constant membrane strain at `(x, y)`.
    fn membrane_uv(eps: &Vector3<f64>, x: f64, y: f64) -> [f64; 2] {
        [eps[0] * x + 0.5 * eps[2] * y, 0.5 * eps[2] * x + eps[1] * y]
    }

    #[test]
    fn test_assemble_star_patch_is_symmetric_and_rigid_body_free() {
        // Task 2.4 self-test (landed with WU6, which needs the assembler): the
        // assembled 48x48 star-patch matrix is symmetric and its six
        // rigid-body fields carry zero energy.
        let pres = star_patch_pres(&STAR_NODES);
        let k = assemble_star_patch(&pres);
        assert_eq!(k.shape(), (48, 48));
        let mut asym = 0.0f64;
        let mut kmax = 0.0f64;
        for i in 0..48 {
            for j in 0..48 {
                asym = asym.max((k[(i, j)] - k[(j, i)]).abs());
                kmax = kmax.max(k[(i, j)].abs());
            }
        }
        assert!(
            asym <= 1e-12 * kmax,
            "assembled star-patch K is not symmetric: {asym:.3e} vs 1e-12 max|K| ({:.3e})",
            1e-12 * kmax
        );

        let mut sym = k.clone();
        for i in 0..48 {
            for j in 0..48 {
                sym[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
            }
        }
        let eig = sym.symmetric_eigen();
        let lam_max = eig.eigenvalues.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
        let rb = rigid_body_modes();
        for j in 0..6 {
            let u = rb.column(j);
            let e = (u.transpose() * &k * u)[(0, 0)];
            let un = u.norm_squared();
            assert!(
                e.abs() <= 1e-12 * lam_max * un,
                "rigid-body column {j} carries energy {e:.3e} > 1e-12 lambda_max ||u||^2 ({:.3e})",
                1e-12 * lam_max * un
            );
        }
    }

    #[test]
    fn test_t1a_membrane_patch_constant_stress_fig5_mesh() {
        // Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
        // Structures 182:404-418, Section 4, pp. 410-411: the membrane patch
        // test on the Fig. 5 mesh (Ko, Lee & Bathe (2017), C&S 182:404-418,
        // Fig. 5, p. 407), with the minimum constraints against rigid-body
        // motion and the boundary nodal forces of the constant stress state.
        //
        // The paper publishes no load magnitudes. The load is derived here
        // element-independently (design 4.2 step 2) as the boundary traction of
        // the constant state integrated with the bilinear shape functions,
        //   f_i = contour_integral N_i (N . n) dGamma,   N = sigma h,
        // 2-point Gauss per boundary edge. Nothing uses the element stiffness,
        // so a wrong element cannot pass by construction; zeroing the load
        // makes the test fail (the non-vacuity control recorded in WU6).
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let pres = star_patch_pres(&nodes);
        let k = assemble_star_patch(&pres);

        let states = [
            ("sigma_xx", [1.0f64, 0.0, 0.0]),
            ("sigma_yy", [0.0, 1.0, 0.0]),
            ("tau_xy", [0.0, 0.0, 1.0]),
        ];
        for (name, sig) in states {
            let sigma = Vector3::new(sig[0], sig[1], sig[2]);
            let eps = membrane_strain(&sigma);

            let f = boundary_integrate(&nodes, |_x, n| {
                [
                    (sigma[0] * n[0] + sigma[2] * n[1]) * h,
                    (sigma[2] * n[0] + sigma[1] * n[1]) * h,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            });

            let mut prescribed = Vec::new();
            for bc in BC_2017_PATCH {
                let uv = membrane_uv(&eps, nodes[bc.node][0], nodes[bc.node][1]);
                let val = match bc.dof {
                    0 => uv[0],
                    1 => uv[1],
                    _ => 0.0,
                };
                prescribed.push((6 * bc.node + bc.dof, val));
            }

            let u = solve_constrained(&k, &f, &prescribed);
            let snorm = sigma.norm();
            let floor = 1e-10 * snorm;
            let mut err_max = 0.0f64;
            let mut vals: Vec<Vector3<f64>> = Vec::new();
            for (e, pe) in pres.iter().enumerate() {
                let ue = star_elem_disp(&u, e);
                let ul = build_t24(pe) * ue;
                for g in 0..GAUSS_XI.len() {
                    let bm = b_membrane_2017(pe, GAUSS_XI[g], GAUSS_ETA[g]);
                    let sig_l = pe.constitutive.cm_raw * (bm * ul);
                    let sig_g = rotate_stress_to_global(&pe.t3, &sig_l);
                    err_max = err_max.max((sig_g - sigma).norm());
                    vals.push(sig_g);
                }
            }
            let spread = (0..3)
                .map(|c| {
                    let mut lo = f64::INFINITY;
                    let mut hi = f64::NEG_INFINITY;
                    for v in &vals {
                        lo = lo.min(v[c]);
                        hi = hi.max(v[c]);
                    }
                    hi - lo
                })
                .fold(0.0f64, f64::max);
            println!(
                "T1.3a {name}: max|sigma_gp - sigma|={err_max:.3e} (rel {:.3e}); \
                 spread={spread:.3e} (rel {:.3e}); floor {floor:.3e}",
                err_max / snorm,
                spread / snorm
            );
            assert!(
                err_max <= 1e-8 * snorm + floor,
                "{name}: recovered membrane stress error {err_max:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
                1e-8 * snorm + floor
            );
            assert!(
                spread <= 1e-8 * snorm + floor,
                "{name}: recovered membrane stress spread {spread:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
                1e-8 * snorm + floor
            );
        }
    }

    /// The exact Mindlin displacement field of a constant curvature state
    /// `kappa = [kxx, kyy, 2 kxy]` in the patch frame, in the sign convention of
    /// the element's own flat bending operator (Ko, Lee & Bathe (2017),
    /// C&S 182:404-418, Eqs. (7c)/(7d), p. 406, verified by the flat reduction
    /// `gamma_13 = w,x + theta_y`, `gamma_23 = w,y - theta_x`):
    /// `kappa_11 = theta_y,x`, `kappa_22 = -theta_x,y`,
    /// `2 kappa_12 = theta_y,y - theta_x,x`. Hence
    /// `w = -1/2 (kxx x^2 + 2 kxy x y + kyy y^2)`,
    /// `theta_x = -(kxy x + kyy y)`, `theta_y = kxx x + kxy y`.
    fn bending_field(kappa: &[f64; 3], x: f64, y: f64) -> [f64; 6] {
        let (kxx, kyy, kxy) = (kappa[0], kappa[1], 0.5 * kappa[2]);
        let w = -0.5 * (kxx * x * x + 2.0 * kxy * x * y + kyy * y * y);
        [0.0, 0.0, w, -(kxy * x + kyy * y), kxx * x + kxy * y, 0.0]
    }

    #[test]
    fn test_t1a_bending_patch_constant_curvature_fig5_mesh() {
        // Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
        // Structures 182:404-418, Section 4, pp. 410-411: the bending patch
        // test on the Fig. 5 mesh, with the minimum constraints and the
        // boundary nodal moments of a constant curvature (constant moment)
        // state (Ko, Lee & Bathe (2017), C&S 182:404-418, Eqs. (7c)/(7d)).
        //
        // The paper publishes no load magnitudes. The load is derived
        // element-independently (design 4.2 step 2) from the constant moment
        // resultant `M = cb kappa`: for the virtual work
        // `integral M : delta kappa dA` the divergence theorem gives the
        // boundary moment tractions
        //   m_theta_x = -(Mxy nx + Myy ny),  m_theta_y = Mxx nx + Mxy ny,
        // with no transverse force (a constant moment field is in equilibrium
        // with `Q_n = 0`). Nothing uses the element stiffness.
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let pres = star_patch_pres(&nodes);
        let k = assemble_star_patch(&pres);
        let cb = shell_iso().cb;

        let states = [
            ("kappa_xx", [1.0e-3f64, 0.0, 0.0]),
            ("kappa_yy", [0.0, 1.0e-3, 0.0]),
            ("kappa_xy", [0.0, 0.0, 1.0e-3]),
        ];
        for (name, kap) in states {
            let kappa = Vector3::new(kap[0], kap[1], kap[2]);
            let m = cb * kappa; // [Mxx, Myy, Mxy]

            let f = boundary_integrate(&nodes, |_x, n| {
                [
                    0.0,
                    0.0,
                    0.0,
                    -(m[2] * n[0] + m[1] * n[1]),
                    m[0] * n[0] + m[2] * n[1],
                    0.0,
                ]
            });

            // The prescribed values are the exact field's values at the
            // constrained DOFs of BC_2017_PATCH; `bending_field` makes them
            // explicit (all happen to vanish here: C is the origin and A lies
            // on the u_y = 0 line).
            let prescribed: Vec<(usize, f64)> = BC_2017_PATCH
                .iter()
                .map(|bc| {
                    let fld = bending_field(&kap, nodes[bc.node][0], nodes[bc.node][1]);
                    (6 * bc.node + bc.dof, fld[bc.dof])
                })
                .collect();

            let u = solve_constrained(&k, &f, &prescribed);
            let knorm = kappa.norm();
            let floor = 1e-10 * knorm;
            let mut err_max = 0.0f64;
            let mut vals: Vec<Vector3<f64>> = Vec::new();
            for (e, pe) in pres.iter().enumerate() {
                let ue = star_elem_disp(&u, e);
                let ul = build_t24(pe) * ue;
                for g in 0..GAUSS_XI.len() {
                    let (b1, _b2) = b_bending_2017(pe, GAUSS_XI[g], GAUSS_ETA[g]);
                    let kap_l = (2.0 / h) * (b1 * ul);
                    let kap_g = rotate_eng3_to_global(&pe.t3, &kap_l);
                    err_max = err_max.max((kap_g - kappa).norm());
                    vals.push(kap_g);
                }
            }
            let spread = (0..3)
                .map(|c| {
                    let mut lo = f64::INFINITY;
                    let mut hi = f64::NEG_INFINITY;
                    for v in &vals {
                        lo = lo.min(v[c]);
                        hi = hi.max(v[c]);
                    }
                    hi - lo
                })
                .fold(0.0f64, f64::max);
            println!(
                "T1.3b {name}: max|kappa_gp - kappa|={err_max:.3e} (rel {:.3e}); \
                 spread={spread:.3e} (rel {:.3e})",
                err_max / knorm,
                spread / knorm
            );
            assert!(
                err_max <= 1e-8 * knorm + floor,
                "{name}: recovered curvature error {err_max:.3e} > 1e-8 ||kappa|| + floor ({:.3e})",
                1e-8 * knorm + floor
            );
            assert!(
                spread <= 1e-8 * knorm + floor,
                "{name}: recovered curvature spread {spread:.3e} > 1e-8 ||kappa|| + floor ({:.3e})",
                1e-8 * knorm + floor
            );
        }
    }

    #[test]
    fn test_t1a_shearing_patch_constant_stress_fig5_mesh() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Fig. 7(c) (p. 5) and
        // Section 3.1 (p. 14); Ko, Lee & Bathe (2017), "A new MITC4+ shell
        // element", Computers and Structures 182:404-418, Section 4
        // (pp. 410-411) and Eqs. (17)-(27). Spec rev 7, Requirement "Tier 1a -
        // shearing patch test".
        //
        // STATE (the 2025 strong form). A constant IN-PLANE shear
        // `tau_xy = tau` with `sigma_xx = sigma_yy = 0` and every
        // transverse-shear and moment resultant zero. Its exact displacement
        // field is the simple shear
        //   u_x = 0,  u_y = (tau / G_xy) x,
        // a constant in-plane shear strain plus the rigid rotation that makes
        // it satisfy the Fig. 7(c) constraint set. This is an in-plane state,
        // NOT a transverse one. (Dvorkin & Bathe (1984), Engineering
        // Computations 1:77-88, Eq. (3), reproduced in Ko, Lee & Bathe (2017),
        // C&S 182:404-418, p. 405, is the element's transverse-shear field and
        // is not exercised by this state.)
        //
        // BOUNDARY CONDITIONS. The figure-read Fig. 7(c) set
        // (`BC_2025_PATCH.shearing`): B(0,10) fully clamped; C(0,0):
        // u_x=u_z=0, theta_x=theta_y=0; the four interior nodes
        // (4,7),(8,7),(8,3),(2,2): u_x=theta_x=theta_y=0; theta_z free except
        // at B; the load at A(10,10) in +y. Every prescribed value is zero
        // because the exact field takes the value zero at every constrained
        // DOF; that is asserted below rather than assumed. (The Tier-1a
        // membrane/bending tests use the derived `BC_2017_PATCH`; this test uses
        // the figure-read Fig. 7(c) set the requirement fixes.)
        //
        // LOAD (well posed). The constant in-plane state's consistent boundary
        // tractions
        //   f_i += integral N_i (sigma . n) dGamma,  sigma = [[0, tau],[tau,0]],
        // integrated with the bilinear boundary shape functions (2-point Gauss
        // per edge), independently of the element stiffness. `sigma_ij,j = 0`
        // for the constant state, so the tractions balance: this derivation IS
        // well posed where the withdrawn transverse one was not. The 2025
        // figure shows a single +y arrow at A and does not publish the load's
        // magnitude or distribution; the constant state's complete
        // boundary-traction vector also loads C and D, and that figure-schematic
        // deviation is recorded (proposal 2.1, spec Evidence gap G9).
        //
        // TOLERANCES (unchanged). tau_xy to 1e-8 relative (absolute floor
        // 1e-10 ||sigma||); the spread of tau_xy across all Gauss points
        // <= 1e-8 ||sigma||; every analytically-zero component (sigma_xx,
        // sigma_yy, the transverse-shear and moment resultants) <= 1e-10 ||sigma||.
        //
        // NON-VACUITY. The test MUST fail if (a) the load is zeroed, or (b) the
        // load is derived from a constant transverse shear resultant (the
        // withdrawn design 4.2 derivation). Both controls are run and asserted
        // below.
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let tau = 1.0f64;
        let sigma = Vector3::new(0.0, 0.0, tau);
        let snorm = sigma.norm();
        let floor = 1e-10 * snorm;
        let pres = star_patch_pres(&nodes);
        let k = assemble_star_patch(&pres);
        let gamma = membrane_strain(&sigma)[2]; // tau / G_xy

        // (1) The state satisfies every Fig. 7(c) constraint: the exact field
        // equals the fixture's prescribed value at each constrained DOF.
        let exact = |x: f64, _y: f64| -> [f64; 6] { [0.0, gamma * x, 0.0, 0.0, 0.0, 0.0] };
        for bc in BC_2025_PATCH.shearing {
            let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
            assert!(
                (e[bc.dof] - bc.value).abs() <= 1e-14,
                "Fig. 7(c) constraint (node {}, dof {}): exact field {} != prescribed {}",
                bc.node,
                bc.dof,
                e[bc.dof],
                bc.value
            );
        }
        let prescribed: Vec<(usize, f64)> = BC_2025_PATCH
            .shearing
            .iter()
            .map(|bc| (6 * bc.node + bc.dof, bc.value))
            .collect();

        // Non-vacuity of the constraint check: the four interior nodes are
        // where Fig. 7(c) pins u_x and leaves u_y free, and the exact field's
        // u_y = gamma x is non-zero there, so the field is non-trivial and the
        // constraint set is not over-constrained.
        for n in [4usize, 5, 6, 7] {
            let uy = exact(nodes[n][0], nodes[n][1])[1];
            assert!(
                uy.abs() > 1e-15,
                "interior node {n}: exact u_y = gamma x must be non-zero, got {uy}"
            );
        }

        // The constant in-plane state's consistent boundary tractions,
        // element-independent (spec Requirement 8, "the load derivation").
        let inplane_load = boundary_integrate(&nodes, |_x, n| {
            [
                (sigma[0] * n[0] + sigma[2] * n[1]) * h,
                (sigma[2] * n[0] + sigma[1] * n[1]) * h,
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        });

        // Recovery at every Gauss point of every element. Returns
        // (tau_xy error, tau_xy spread, max analytically-zero membrane
        // component, max bending (moment) component, max transverse-shear
        // resultant component).
        let measure = |u: &DVector<f64>| -> (f64, f64, f64, f64, f64) {
            let mut err = 0.0f64;
            let mut zero_mem = 0.0f64;
            let mut bending = 0.0f64;
            let mut shear = 0.0f64;
            let mut taus: Vec<f64> = Vec::new();
            for (e, pe) in pres.iter().enumerate() {
                let ue = star_elem_disp(u, e);
                let ul = build_t24(pe) * ue;
                for g in 0..GAUSS_XI.len() {
                    let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                    let sig_m = pe.constitutive.cm_raw * (b_membrane_2017(pe, r, s) * ul);
                    let sig_g = rotate_stress_to_global(&pe.t3, &sig_m);
                    err = err.max((sig_g[2] - tau).abs());
                    zero_mem = zero_mem.max(sig_g[0].abs()).max(sig_g[1].abs());
                    taus.push(sig_g[2]);
                    // The moment resultant, in stress units: the through-thickness
                    // bending stress at the surface t = 1, `cm_raw (e_b1 + e_b2)`.
                    let (b1, b2) = b_bending_2017(pe, r, s);
                    let sig_b = pe.constitutive.cm_raw * ((b1 + b2) * ul);
                    let sig_bg = rotate_stress_to_global(&pe.t3, &sig_b);
                    for c in 0..3 {
                        bending = bending.max(sig_bg[c].abs());
                    }
                    // The transverse-shear resultant `Q = G h gamma`
                    // (uncorrected, Ko, Lee & Bathe (2017), C&S 182:404-418,
                    // p. 410), in force-per-length; with h = 1 it is the
                    // recovered transverse shear stress.
                    let q = pe.cs_uncorrected * (b_shear_mitc4(pe, r, s) * ul);
                    let qg = rotate_shear_to_global(&pe.t3, &[q[0], q[1]]);
                    shear = shear.max(qg[0].abs()).max(qg[1].abs());
                }
            }
            let spread = {
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for v in &taus {
                    lo = lo.min(*v);
                    hi = hi.max(*v);
                }
                hi - lo
            };
            (err, spread, zero_mem, bending, shear)
        };

        let u = solve_constrained(&k, &inplane_load, &prescribed);
        let (err, spread, zero_mem, bending, shear) = measure(&u);
        println!(
            "T1.3c in-plane tau_xy={tau}: max|tau_gp - tau|={:.3e} (rel {:.3e}); \
             spread={:.3e} (rel {:.3e}); zero membrane={:.3e}; \
             moment/bending={:.3e}; transverse shear={:.3e}; floor={floor:.3e}",
            err,
            err / snorm,
            spread,
            spread / snorm,
            zero_mem,
            bending,
            shear
        );
        assert!(
            err <= 1e-8 * snorm + floor,
            "tau_xy: recovered error {err:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
            1e-8 * snorm + floor
        );
        assert!(
            spread <= 1e-8 * snorm + floor,
            "tau_xy: recovered spread {spread:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
            1e-8 * snorm + floor
        );
        assert!(
            zero_mem <= 1e-10 * snorm,
            "sigma_xx/sigma_yy: {zero_mem:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            bending <= 1e-10 * snorm,
            "moment resultant (surface bending stress): {bending:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            shear <= 1e-10 * snorm,
            "transverse-shear resultant: {shear:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );

        // (2) Non-vacuity (a): zeroing the load recovers a zero state, not the
        // prescribed tau_xy, so the tolerance above fails.
        let u_zero = solve_constrained(&k, &DVector::zeros(48), &prescribed);
        let (err_zero, _, _, _, _) = measure(&u_zero);
        println!(
            "T1.3c non-vacuity (a) zeroed load: max|tau_gp - tau|={:.3e} (rel {:.3e})",
            err_zero,
            err_zero / snorm
        );
        assert!(
            err_zero / snorm > 1e-3,
            "zeroed-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err_zero / snorm
        );

        // (3) Non-vacuity (b): the withdrawn design 4.2 derivation, a load from
        // a constant TRANSVERSE shear resultant q = G h gamma, does not
        // reproduce the in-plane state (spec Evidence gap G9). The same Fig.
        // 7(c) BC set is used, so the only difference from the main solve is
        // the load derivation.
        let gamma_t = 1.0e-3f64;
        let q = [2.0e11 / (2.0 * (1.0 + 0.3)) * gamma_t, 0.0];
        let transverse_load = boundary_integrate(&nodes, |_x, n| {
            [0.0, 0.0, q[0] * n[0] + q[1] * n[1], 0.0, 0.0, 0.0]
        });
        let u_trans = solve_constrained(&k, &transverse_load, &prescribed);
        let (err_trans, _, _, _, shear_trans) = measure(&u_trans);
        println!(
            "T1.3c non-vacuity (b) withdrawn transverse load: max|tau_gp - tau|={:.3e} \
             (rel {:.3e}); transverse shear={:.3e}",
            err_trans,
            err_trans / snorm,
            shear_trans
        );
        assert!(
            err_trans / snorm > 1e-3,
            "withdrawn transverse-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err_trans / snorm
        );
    }

    // ========================================================================
    // WU7 (tasks 8.1-8.4) - Tier 1b: the 2025 six-DOF element's own basic
    // tests.  Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell
    // elements with six degrees of freedom at each node - the MITC4/D and
    // MITC4+/D elements", Computers and Structures 308:107622, Section 3.1
    // (pp. 13-14) and Fig. 7 (p. 5).  The citations are self-contained and
    // anchored to `docs/references.md`; the "paper A/B" shorthand is not used.
    // ========================================================================

    /// One Gauss point's strong-form state, in the patch's global frame.
    struct StrongFormGp {
        /// Mid-surface membrane stress `cm_raw (e_m)`.
        membrane: Vector3<f64>,
        /// Surface (`t = 1`) bending stress `cm_raw (e_b1 + e_b2)`.
        bending: Vector3<f64>,
        /// Uncorrected transverse-shear resultant `Q = G h gamma`
        /// (Dvorkin & Bathe (1984), Engineering Computations 1:77-88, Eq. (3),
        /// reproduced in Ko, Lee & Bathe (2017), C&S 182:404-418, p. 405).
        shear: [f64; 2],
        /// Curvature `(2/h) e_b1`.
        curvature: Vector3<f64>,
    }

    /// The strong-form state at every Gauss point of every element of the star
    /// patch, recovered from a 48-DOF patch solution `u`.
    fn strong_form_gauss_points(nodes: &[[f64; 3]; 8], u: &DVector<f64>) -> Vec<StrongFormGp> {
        let pres = star_patch_pres(nodes);
        let h = 1.0f64;
        let mut out = Vec::new();
        for (e, pe) in pres.iter().enumerate() {
            let ue = star_elem_disp(u, e);
            let ul = build_t24(pe) * ue;
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let sig_m = pe.constitutive.cm_raw * (b_membrane_2017(pe, r, s) * ul);
                let (b1, b2) = b_bending_2017(pe, r, s);
                let sig_b = pe.constitutive.cm_raw * ((b1 + b2) * ul);
                let kap_l = (2.0 / h) * (b1 * ul);
                let q = pe.cs_uncorrected * (b_shear_mitc4(pe, r, s) * ul);
                out.push(StrongFormGp {
                    membrane: rotate_stress_to_global(&pe.t3, &sig_m),
                    bending: rotate_stress_to_global(&pe.t3, &sig_b),
                    shear: rotate_shear_to_global(&pe.t3, &[q[0], q[1]]),
                    curvature: rotate_eng3_to_global(&pe.t3, &kap_l),
                });
            }
        }
        out
    }

    /// Solve the star-patch strong-form problem on `nodes`: assemble the 48-DOF
    /// patch stiffness and solve with the prescribed DOF values `prescribed`.
    fn solve_strong_patch(
        nodes: &[[f64; 3]; 8],
        f: &DVector<f64>,
        prescribed: &[(usize, f64)],
    ) -> DVector<f64> {
        let k = assemble_star_patch(&star_patch_pres(nodes));
        solve_constrained(&k, f, prescribed)
    }

    /// Maximum absolute value of a 48-DOF vector.
    fn max_abs_vec(u: &DVector<f64>) -> f64 {
        u.iter().fold(0.0f64, |m, &v| m.max(v.abs()))
    }

    /// Maximum absolute entrywise difference of two 48-DOF vectors.
    fn max_abs_diff_vec(a: &DVector<f64>, b: &DVector<f64>) -> f64 {
        a.iter()
            .zip(b.iter())
            .fold(0.0f64, |m, (x, y)| m.max((x - y).abs()))
    }

    /// The constant-`sigma_xx` extension state's consistent boundary-traction
    /// load, derived element-independently from the constant stress state
    /// (Ko, Lee & Bathe (2017), C&S 182:404-418, Section 4, p. 410).
    fn extension_load(nodes: &[[f64; 3]; 8], sigma0: f64, h: f64) -> DVector<f64> {
        boundary_integrate(nodes, |_x, n| [sigma0 * n[0] * h, 0.0, 0.0, 0.0, 0.0, 0.0])
    }

    /// The constant-moment `Mxx` bending state's consistent boundary-moment
    /// load `m_theta_x = -(Mxy nx + Myy ny)`, `m_theta_y = Mxx nx + Mxy ny`
    /// (Ko, Lee & Bathe (2017), C&S 182:404-418, Eqs. (7c)/(7d), p. 406).
    fn bending_load(nodes: &[[f64; 3]; 8], kappa0: f64) -> DVector<f64> {
        let m = shell_iso().cb * Vector3::new(kappa0, 0.0, 0.0);
        boundary_integrate(nodes, |_x, n| {
            [
                0.0,
                0.0,
                0.0,
                -(m[2] * n[0] + m[1] * n[1]),
                m[0] * n[0] + m[2] * n[1],
                0.0,
            ]
        })
    }

    /// The constant in-plane `tau_xy` shearing state's consistent
    /// boundary-traction load (Ko, Bathe & Zhang (2025), C&S 308:107622,
    /// Fig. 7(c), p. 5; Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410).
    fn shearing_load(nodes: &[[f64; 3]; 8], tau: f64, h: f64) -> DVector<f64> {
        boundary_integrate(nodes, |_x, n| {
            [tau * n[1] * h, tau * n[0] * h, 0.0, 0.0, 0.0, 0.0]
        })
    }

    /// The distance of a 24-DOF field from the six-dimensional rigid-body
    /// space, relative to the field's L2 norm.
    fn distance_from_rigid_body_space(pre: &Mitc4PlusDPrecomputed, u: &[f64; 24]) -> f64 {
        let rb = rigid_body_fields(pre);
        let mut g = DMatrix::<f64>::zeros(6, 6);
        let mut rhs = DVector::<f64>::zeros(6);
        for a in 0..6 {
            for b in 0..6 {
                g[(a, b)] = rb[a].iter().zip(rb[b].iter()).map(|(x, y)| x * y).sum();
            }
            rhs[a] = rb[a].iter().zip(u.iter()).map(|(x, y)| x * y).sum();
        }
        let coef: DVector<f64> = g
            .lu()
            .solve(&rhs)
            .expect("the six rigid-body fields are linearly independent");
        let mut resid = 0.0f64;
        for a in 0..6 {
            for j in 0..24 {
                resid += (u[j] - coef[a] * rb[a][j]).powi(2);
            }
        }
        resid.sqrt() / norm2(u).sqrt().max(1e-30)
    }

    // ------------------------------------------------------------------
    // 8.1 - spatial isotropy of the 2025 six-DOF element (T1.4)
    // ------------------------------------------------------------------

    #[test]
    fn test_t1b_spatial_isotropy() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1, p. 13:
        // "MITC/D and MITC4+/D elements pass the spatial isotropy test."
        //
        // Spec rev 6, Requirement 9. The 24x24 stiffness is formed in global
        // coordinates with the drilling DOF free at all nodes and compared
        // across at least four global orientations; for a co-rotated field that
        // includes a non-zero drilling component the strain energy `u^T K u`
        // must be invariant to 1e-10 relative. The drilling component and the
        // field's drill-block energy are asserted non-zero, so the drill DOF is
        // exercised rather than bypassed.
        let axis = [1.0, 2.0, 3.0];
        let angles = [
            0.0,
            0.7,
            core::f64::consts::FRAC_PI_4,
            core::f64::consts::FRAC_PI_2,
        ];
        let geoms = [
            ("flat-square", FLAT_SQUARE),
            ("flat-distorted", FLAT_DISTORTED),
            ("ruled-warped", RULED_WARPED),
        ];

        // A field that explicitly includes a non-zero drilling component: a
        // generic field plus a constant drill rotation at every node.
        let mut u_ref = generic_dofs();
        for i in 0..4 {
            u_ref[6 * i + 5] += 0.13;
        }

        for (name, c) in geoms {
            let pre0 = pre_from(&c);
            let k0 = compute_ke_global(&pre0);
            let lam0 = sorted_eigenvalues(&k0);
            let lam_max0 = lam0.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
            let q0 = quadratic(&u_ref, &k0);

            // Non-vacuity: the field's drilling component `theta_i . V^D` and
            // its drill-block energy are non-zero.
            let kd0 = drill_ke_local(&pre0);
            let mut worst_drill = 0.0f64;
            for i in 0..4 {
                let theta = Vector3::new(u_ref[6 * i + 3], u_ref[6 * i + 4], u_ref[6 * i + 5]);
                worst_drill = worst_drill.max(pre0.v_d.dot(&theta).abs());
            }
            let u_loc = build_t24(&pre0) * Vec24::from_column_slice(&u_ref);
            let mut u_loc_arr = [0.0f64; 24];
            u_loc_arr.copy_from_slice(u_loc.as_slice());
            let drill_energy = quadratic(&u_loc_arr, &kd0).abs();
            assert!(
                worst_drill > 1e-3,
                "{name}: the field's drilling component is zero ({worst_drill:.3e})"
            );
            assert!(
                drill_energy > 0.0,
                "{name}: the field carries no drill-block energy"
            );

            let mut worst_eig = 0.0f64;
            let mut worst_energy = 0.0f64;
            for &theta in angles.iter() {
                let r = rotation_matrix(axis, theta);
                let pre = pre_from(&rotate_geom(&c, &r));
                let k = compute_ke_global(&pre);
                let lam = sorted_eigenvalues(&k);
                for i in 0..24 {
                    worst_eig = worst_eig.max((lam[i] - lam0[i]).abs());
                }
                let u_rot = rotate_dofs(&u_ref, &r);
                let q = quadratic(&u_rot, &k);
                worst_energy = worst_energy.max((q - q0).abs() / q0.abs().max(1e-30));
            }
            println!(
                "T1.4 {name}: worst |dlambda|={worst_eig:.3e} (1e-10 lambda_max={:.3e}); \
                 worst |duKu|/|uKu|={worst_energy:.3e}; worst drill component={worst_drill:.3e}; \
                 drill-block energy={drill_energy:.3e}",
                1e-10 * lam_max0
            );
            assert!(
                worst_eig <= 1e-10 * lam_max0,
                "{name}: orientation changes an eigenvalue by {worst_eig:.3e} > 1e-10 lambda_max ({:.3e})",
                1e-10 * lam_max0
            );
            assert!(
                worst_energy <= 1e-10,
                "{name}: the co-rotated strain energy with a drilling component is not invariant ({worst_energy:.3e} relative)"
            );
        }
    }

    // ------------------------------------------------------------------
    // 8.2 - six-or-seven zero-energy modes with the drill DOF (T1.5)
    // ------------------------------------------------------------------

    #[test]
    fn test_t1b_zero_energy_modes_six_or_seven_with_drill_dof() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1,
        // pp. 13-14: "The elements pass the zero energy mode tests, and rigid
        // body modes are properly represented."
        //
        // Spec rev 6, Requirement 10. With the drill DOF constrained the way
        // the paper's own patch tests constrain it (`theta_z` free at every
        // node except corner B; Ko, Bathe & Zhang (2025), C&S 308:107622,
        // Fig. 7(b)(c)(d), B Section 3.1 pp. 13-14), the count of eigenvalues
        // with `|lambda| <= 1e-10 lambda_max` is exactly 7 on the flat
        // rectangle and exactly 6 on each of the flat distorted, ruled-warped
        // and doubly-warped quads. The rigid-body fields include the `theta_z`
        // component of the rigid rotation about `V_n`; each of the six has
        // energy `<= 1e-12 lambda_max` at `||u|| = 1` (with the constraint
        // applied, the five fields with `theta_z(B) = 0`). The flat rectangle's
        // surplus is the second, curl-induced direction of the drill operator's
        // two-dimensional null space (Eq. (19d)), which the paper's single-
        // corner condition -- one linear condition -- cannot remove. Eq. (19b)'s
        // columns are DIFFERENCES of edge terms, so a constant drill rotation is
        // a zero-energy direction of the paper's own equations. No penalty,
        // constraint or numerical factor is added.
        let geoms = [
            // (name, geometry, expected constrained zero-count, expected drill
            //  null-space dimension)
            ("flat-rectangle", RECT, 7usize, 2usize),
            ("flat-distorted", FLAT_DISTORTED, 6, 1),
            ("ruled-warped", RULED_WARPED, 6, 1),
            ("doubly-warped", DOUBLY_WARPED, 6, 1),
        ];
        let drill_fixed = [6 * 3 + 5];
        let mut failures: Vec<String> = Vec::new();

        for (name, c, expected_count, expected_null) in geoms {
            let pre = pre_from(&c);
            let k = compute_ke_local(&pre);
            let lam_max = lambda_max(&k);

            // (1) The zero-eigenvalue count under the paper's drill constraint.
            let lam_c = reduced_eigenvalues(&k, &drill_fixed);
            let lam_c_max = lam_c.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
            let zeros = lam_c
                .iter()
                .filter(|v| v.abs() <= 1e-10 * lam_c_max)
                .count();
            let sep = lam_c
                .iter()
                .map(|v| v.abs())
                .nth(6)
                .expect("the constrained 23x23 system has at least seven eigenvalues");

            // (2) The rigid-body representation: energy `<= 1e-12 lambda_max`
            // at `||u|| = 1`, for the six fields with `theta_z` free and the
            // five satisfying `theta_z(B) = 0` under the paper's constraint.
            let mut worst_rb = 0.0f64;
            let mut worst_rb_c = 0.0f64;
            let mut surviving = 0usize;
            for u in rigid_body_fields(&pre).iter() {
                let nrm = norm2(u).sqrt();
                let mut un = *u;
                for v in un.iter_mut() {
                    *v /= nrm;
                }
                worst_rb = worst_rb.max(quadratic(&un, &k).abs() / lam_max);
                if u[6 * 3 + 5].abs() <= 1e-12 {
                    surviving += 1;
                    let mut uc = un;
                    uc[6 * 3 + 5] = 0.0;
                    let nc = norm2(&uc).sqrt();
                    for v in uc.iter_mut() {
                        *v /= nc;
                    }
                    worst_rb_c = worst_rb_c.max(quadratic(&uc, &k).abs() / lam_max);
                }
            }

            // (3) The drill block's own null space, Eq. (19a)/(19b).
            let (lam4, null) = drill_block_null_space(&pre, 0);
            let kd = drill_ke_local(&pre);
            let vd = Vector3::new(
                pre.v_d.dot(&pre.e1),
                pre.v_d.dot(&pre.e2),
                pre.v_d.dot(&pre.e3),
            );
            let mut worst_null_energy = 0.0f64;
            let mut worst_dist = 0.0f64;
            for n in null.iter() {
                let mut u = [0.0f64; 24];
                for i in 0..4 {
                    for b in 0..3 {
                        u[6 * i + 3 + b] = n[i] * vd[b];
                    }
                }
                worst_null_energy =
                    worst_null_energy.max(quadratic(&u, &kd).abs() / (lam_max * norm2(&u)));
                worst_dist = worst_dist.max(distance_from_rigid_body_space(&pre, &u));
            }

            // (4) The drilling DOF carries non-zero stiffness in at least one
            // non-rigid mode: the drill block's largest eigenvector, lifted to
            // a pure drill-rotation field, is not in the rigid-body space.
            let blk = drill_block_eq19(&pre, 0);
            let eig = blk.symmetric_eigen();
            let imax = (0..4)
                .max_by(|&a, &b| {
                    eig.eigenvalues[a]
                        .abs()
                        .total_cmp(&eig.eigenvalues[b].abs())
                })
                .expect("the 4x4 drill block has four eigenvalues");
            let mut u_max = [0.0f64; 24];
            for i in 0..4 {
                let n_i = eig.eigenvectors[(i, imax)];
                for b in 0..3 {
                    u_max[6 * i + 3 + b] = n_i * vd[b];
                }
            }
            let drill_energy = quadratic(&u_max, &kd).abs();
            let dist_max = distance_from_rigid_body_space(&pre, &u_max);

            // (5) Non-vacuity: an inert drill block would give null dim 4 and
            // a rank gain would drop it below the paper's value.
            let (_li, null_inert) = drill_block_null_space(&pre, 1);
            let (_lr, null_rank) = drill_block_null_space(&pre, 2);

            println!(
                "T1.5 {name}: lambda_max={lam_max:.6e}; constrained zero-count={zeros} \
                 (expected {expected_count}); |lambda_7|={sep:.6e} ({lam_c_max:.3e} lambda_max); \
                 six-field rb energy/lambda_max={worst_rb:.3e}; constrained {worst_rb_c:.3e} \
                 (fields with theta_z(B)=0: {surviving}/6); drill block |lambda_1..4|={lam4:?} \
                 null-dim={} (expected {expected_null}); inert={} rank-gain={}; \
                 drill-null energy/(lambda_max ||u||^2)={worst_null_energy:.3e}; \
                 drill-null distance={worst_dist:.3e}; largest drill mode energy={drill_energy:.3e} \
                 distance={dist_max:.3e}",
                null.len(),
                null_inert.len(),
                null_rank.len(),
            );

            if zeros != expected_count {
                failures.push(format!(
                    "{name}: expected exactly {expected_count} zero eigenvalues with the drill \
                     DOF constrained as Ko, Bathe & Zhang (2025), C&S 308:107622, Fig. 7 \
                     constrains it, got {zeros} (lambda_max {lam_c_max:.3e})"
                ));
            }
            if expected_count == 6 && sep < 1e-9 * lam_c_max {
                failures.push(format!(
                    "{name}: the seventh eigenvalue {sep:.3e} is not separated from the six \
                     zero modes by 1e-9 lambda_max ({:.3e})",
                    1e-9 * lam_c_max
                ));
            }
            if worst_rb > 1e-12 {
                failures.push(format!(
                    "{name}: rigid-body energy/lambda_max = {worst_rb:.3e} at ||u|| = 1"
                ));
            }
            if surviving != 5 {
                failures.push(format!(
                    "{name}: expected the single-corner theta_z(B) = 0 condition to leave \
                     exactly five rigid-body fields, got {surviving}"
                ));
            }
            if worst_rb_c > 1e-12 {
                failures.push(format!(
                    "{name}: with the constraint applied, rigid-body energy/lambda_max = \
                     {worst_rb_c:.3e}"
                ));
            }
            if null.len() != expected_null {
                failures.push(format!(
                    "{name}: the drill block null-space dimension is {}, expected {expected_null} \
                     (constant drill rotation{}); lambda_4 = {lam4:?}",
                    null.len(),
                    if expected_null == 2 {
                        " + Eq. (19d) curl-induced theta_z hourglass"
                    } else {
                        ""
                    }
                ));
            }
            if worst_null_energy > 1e-12 {
                failures.push(format!(
                    "{name}: a drill null vector carries energy {worst_null_energy:.3e}"
                ));
            }
            if worst_dist <= 2.4 {
                failures.push(format!(
                    "{name}: a drill null vector lies in the rigid-body space (distance \
                     {worst_dist:.3e})"
                ));
            }
            if eig.eigenvalues[imax].abs() <= 1e-10 * lam_max {
                failures.push(format!(
                    "{name}: the drill block has no non-zero eigenvalue (largest {:.3e})",
                    eig.eigenvalues[imax].abs()
                ));
            }
            if drill_energy <= 0.0 {
                failures.push(format!(
                    "{name}: the drilling DOF carries no stiffness in its largest mode"
                ));
            }
            if dist_max <= 1e-6 {
                failures.push(format!(
                    "{name}: the drill block's largest mode is a rigid-body field (distance \
                     {dist_max:.3e})"
                ));
            }
            if null.len() >= 4 {
                failures.push(format!(
                    "{name}: the drill block is inert (null-space dimension {})",
                    null.len()
                ));
            }
            if null_inert.len() != 4 {
                failures.push(format!(
                    "{name}: an inert drill block must have null-space dimension 4, got {}",
                    null_inert.len()
                ));
            }
            if null_rank.len() != 0 {
                failures.push(format!(
                    "{name}: a rank gain in the drill block must drop its null-space dimension \
                     to 0, got {}",
                    null_rank.len()
                ));
            }
        }

        assert!(
            failures.is_empty(),
            "Tier-1b zero-energy-mode test failures:\n{}",
            failures.join("\n")
        );
    }

    // ------------------------------------------------------------------
    // 8.3 - the strong-form patch tests: extension, bending, shearing (T1.6a-c)
    // ------------------------------------------------------------------

    #[test]
    fn test_t1b_strong_patch_extension_constant_and_zero_stress() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1, p. 14:
        // "we consider the 'strong form' of the patch tests ... we require that
        // the calculations give the analytical solutions of constant and zero
        // stresses throughout the patch due to the applied loading."  The
        // extension BC set is Fig. 7(d) (p. 5): B fully clamped; C: u_x = u_z =
        // 0; the load at A in +x; theta_z free except at B (B Section 3.1
        // pp. 13-14).
        //
        // Spec rev 8, Requirement 11. The constant axial state `sigma_xx` is
        // reproduced to 1e-8 relative in every element, and every
        // analytically-zero component (sigma_yy, tau_xy, the moment resultants
        // and the transverse-shear resultants) is <= 1e-10 of |sigma_axial|
        // throughout the patch. The load is the constant state's consistent
        // boundary tractions (Ko, Lee & Bathe (2017), C&S 182:404-418,
        // Section 4, p. 410), derived element-independently. Zeroing the load
        // is the non-vacuity control.
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let sigma0 = 1.0f64;
        let sigma = Vector3::new(sigma0, 0.0, 0.0);
        let snorm = sigma.norm();
        let floor = 1e-10 * snorm;
        let eps = membrane_strain(&sigma);

        // The exact constant-strain field, made admissible by the rigid
        // translation `u_y = eps_yy (y - 10)` that satisfies B(0,10) fully
        // clamped (`eps_xy = 0` for `sigma_xx` alone).
        let exact = |x: f64, y: f64| -> [f64; 6] {
            [
                eps[0] * x + 0.5 * eps[2] * y,
                0.5 * eps[2] * x + eps[1] * (y - 10.0),
                0.0,
                0.0,
                0.0,
                0.0,
            ]
        };
        for bc in BC_2025_PATCH.extension {
            let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
            assert!(
                (e[bc.dof] - bc.value).abs() <= 1e-14,
                "Fig. 7(d) constraint (node {}, dof {}): exact field {} != prescribed {}",
                bc.node,
                bc.dof,
                e[bc.dof],
                bc.value
            );
        }
        let prescribed: Vec<(usize, f64)> = BC_2025_PATCH
            .extension
            .iter()
            .map(|bc| {
                let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
                (6 * bc.node + bc.dof, e[bc.dof])
            })
            .collect();

        let f = extension_load(&nodes, sigma0, h);
        let u = solve_strong_patch(&nodes, &f, &prescribed);
        let gps = strong_form_gauss_points(&nodes, &u);

        let mut err = 0.0f64;
        let mut zero_mem = 0.0f64;
        let mut bend = 0.0f64;
        let mut shear = 0.0f64;
        for gp in &gps {
            err = err.max((gp.membrane[0] - sigma0).abs());
            zero_mem = zero_mem.max(gp.membrane[1].abs()).max(gp.membrane[2].abs());
            for c in 0..3 {
                bend = bend.max(gp.bending[c].abs());
            }
            shear = shear.max(gp.shear[0].abs()).max(gp.shear[1].abs());
        }
        println!(
            "T1.6a extension sigma_xx={sigma0}: max|sigma_xx,gp - sigma|={err:.3e} (rel {:.3e}); \
             zero membrane={zero_mem:.3e}; moment={bend:.3e}; transverse shear={shear:.3e}; \
             floor={floor:.3e}",
            err / snorm
        );
        assert!(
            err <= 1e-8 * snorm + floor,
            "sigma_xx: recovered error {err:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
            1e-8 * snorm + floor
        );
        assert!(
            zero_mem <= 1e-10 * snorm,
            "sigma_yy/tau_xy: {zero_mem:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            bend <= 1e-10 * snorm,
            "moment resultant: {bend:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            shear <= 1e-10 * snorm,
            "transverse-shear resultant: {shear:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );

        // Non-vacuity: zeroing the load recovers a zero state, not sigma_xx.
        let u0 = solve_strong_patch(&nodes, &DVector::zeros(48), &prescribed);
        let gps0 = strong_form_gauss_points(&nodes, &u0);
        let err0 = gps0
            .iter()
            .map(|g| (g.membrane[0] - sigma0).abs())
            .fold(0.0f64, f64::max);
        println!(
            "T1.6a non-vacuity zeroed load: max|sigma_xx,gp - sigma|={err0:.3e} (rel {:.3e})",
            err0 / snorm
        );
        assert!(
            err0 / snorm > 1e-3,
            "zeroed-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err0 / snorm
        );
    }

    #[test]
    fn test_t1b_strong_patch_bending_constant_and_zero_stress() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1, p. 14
        // (the strong form: "constant and zero stresses throughout the patch"),
        // with the Fig. 7(b) (p. 5) bending BC set: B fully clamped; C:
        // u_x = u_z = 0, theta_y = 0; theta_z free except at B (B Section 3.1
        // pp. 13-14).
        //
        // Spec rev 8, Requirement 11. The bending-produced constant-moment
        // state is reproduced to 1e-8 relative in every element, and every
        // analytically-zero component is <= 1e-10 of the state's magnitude.
        //
        // ADMISSIBILITY (recorded). The strong form requires the solution to
        // BE the constant straining mode, so the exact field must satisfy the
        // Fig. 7(b) zero BCs. Of the three independent constant-moment states,
        // only `Mxx` (kappa = [kappa_xx, 0, 0]) does: `w = -kappa_xx x^2/2`,
        // `theta_y = kappa_xx x` vanish at B(0,10) and at C(0,0). For
        // `kappa_yy` alone the exact field needs `theta_x = -kappa_yy y`, and
        // the rigid motion that makes `theta_x(B) = 0` and `w(C) = 0` leaves
        // `w(B) = 50 kappa_yy != 0`; for `kappa_xy` alone it needs
        // `theta_y = kappa_xy y`, and the motion that makes `theta_y(B) = 0`
        // leaves `theta_y(C) = -10 kappa_xy != 0`. Those two states are
        // therefore not admissible under Fig. 7(b)'s zero-valued constraints,
        // and the strong-form bending state realized is `Mxx`. The load is that
        // constant state's consistent boundary moments (Ko, Lee & Bathe (2017),
        // C&S 182:404-418, Eqs. (7c)/(7d), p. 406), derived element-
        // independently; zeroing it is the non-vacuity control.
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let kappa0 = 1.0e-3f64;
        let kappa = Vector3::new(kappa0, 0.0, 0.0);
        let knorm = kappa.norm();
        let floor = 1e-10 * knorm;
        // The state's stress magnitude: the surface bending stress of the
        // prescribed constant curvature.
        let sig_scale = (shell_iso().cm_raw * (kappa * (0.5 * h))).norm();

        let exact = |x: f64, y: f64| -> [f64; 6] { bending_field(&[kappa0, 0.0, 0.0], x, y) };
        for bc in BC_2025_PATCH.bending {
            let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
            assert!(
                (e[bc.dof] - bc.value).abs() <= 1e-14,
                "Fig. 7(b) constraint (node {}, dof {}): exact field {} != prescribed {}",
                bc.node,
                bc.dof,
                e[bc.dof],
                bc.value
            );
        }
        let prescribed: Vec<(usize, f64)> = BC_2025_PATCH
            .bending
            .iter()
            .map(|bc| {
                let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
                (6 * bc.node + bc.dof, e[bc.dof])
            })
            .collect();

        let f = bending_load(&nodes, kappa0);
        let u = solve_strong_patch(&nodes, &f, &prescribed);
        let gps = strong_form_gauss_points(&nodes, &u);

        let mut err = 0.0f64;
        let mut zero_kap = 0.0f64;
        let mut zero_mem = 0.0f64;
        let mut shear = 0.0f64;
        for gp in &gps {
            err = err.max((gp.curvature[0] - kappa0).abs());
            zero_kap = zero_kap
                .max(gp.curvature[1].abs())
                .max(gp.curvature[2].abs());
            zero_mem = zero_mem
                .max(gp.membrane[0].abs())
                .max(gp.membrane[1].abs())
                .max(gp.membrane[2].abs());
            shear = shear.max(gp.shear[0].abs()).max(gp.shear[1].abs());
        }
        println!(
            "T1.6b bending Mxx={kappa0:.3e}: max|kappa_xx,gp - kappa|={err:.3e} (rel {:.3e}); \
             zero curvature={zero_kap:.3e}; zero membrane={zero_mem:.3e}; transverse shear={shear:.3e}; \
             stress scale={sig_scale:.3e}; floor={floor:.3e}",
            err / knorm
        );
        assert!(
            err <= 1e-8 * knorm + floor,
            "kappa_xx: recovered error {err:.3e} > 1e-8 ||kappa|| + floor ({:.3e})",
            1e-8 * knorm + floor
        );
        assert!(
            zero_kap <= 1e-10 * knorm,
            "kappa_yy/2kappa_xy: {zero_kap:.3e} > 1e-10 ||kappa|| ({:.3e})",
            1e-10 * knorm
        );
        assert!(
            zero_mem <= 1e-10 * sig_scale,
            "membrane stress: {zero_mem:.3e} > 1e-10 of the state ({:.3e})",
            1e-10 * sig_scale
        );
        assert!(
            shear <= 1e-10 * sig_scale,
            "transverse-shear resultant: {shear:.3e} > 1e-10 of the state ({:.3e})",
            1e-10 * sig_scale
        );

        // Non-vacuity: zeroing the load recovers a zero state, not Mxx.
        let u0 = solve_strong_patch(&nodes, &DVector::zeros(48), &prescribed);
        let gps0 = strong_form_gauss_points(&nodes, &u0);
        let err0 = gps0
            .iter()
            .map(|g| (g.curvature[0] - kappa0).abs())
            .fold(0.0f64, f64::max);
        println!(
            "T1.6b non-vacuity zeroed load: max|kappa_xx,gp - kappa|={err0:.3e} (rel {:.3e})",
            err0 / knorm
        );
        assert!(
            err0 / knorm > 1e-3,
            "zeroed-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err0 / knorm
        );
    }

    #[test]
    fn test_t1b_strong_patch_shearing_constant_and_zero_stress() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1, p. 14
        // (the strong form: "constant and zero stresses throughout the patch"),
        // Fig. 7(c) (p. 5).
        //
        // Spec rev 8, Requirement 11. The 2025 paper's own shearing patch is
        // the IN-PLANE patch: a constant in-plane `tau_xy` with
        // `sigma_xx = sigma_yy = 0` and every transverse-shear and moment
        // resultant zero -- NOT a transverse one (the transverse reading was
        // withdrawn; Evidence gap G9). The exact field is the simple shear
        // `u_x = 0`, `u_y = (tau / G_xy) x`. The BC set is the figure-read
        // Fig. 7(c) set (`BC_2025_PATCH.shearing`: B fully clamped; C:
        // u_x = u_z = 0, theta_x = theta_y = 0; the four interior nodes
        // u_x = theta_x = theta_y = 0; theta_z free except at B; the load at A
        // in +y). The load is the constant in-plane state's consistent
        // boundary tractions, an equilibrium state (`sigma_ij,j = 0`).
        //
        // Non-vacuity: the test fails if the load is zeroed, and fails if the
        // load is derived from a constant transverse shear resultant (the
        // withdrawn design 4.2 derivation).
        let nodes = STAR_NODES;
        let h = 1.0f64;
        let tau = 1.0f64;
        let sigma = Vector3::new(0.0, 0.0, tau);
        let snorm = sigma.norm();
        let floor = 1e-10 * snorm;
        let gamma = membrane_strain(&sigma)[2];

        let exact = |x: f64, _y: f64| -> [f64; 6] { [0.0, gamma * x, 0.0, 0.0, 0.0, 0.0] };
        for bc in BC_2025_PATCH.shearing {
            let e = exact(nodes[bc.node][0], nodes[bc.node][1]);
            assert!(
                (e[bc.dof] - bc.value).abs() <= 1e-14,
                "Fig. 7(c) constraint (node {}, dof {}): exact field {} != prescribed {}",
                bc.node,
                bc.dof,
                e[bc.dof],
                bc.value
            );
        }
        let prescribed: Vec<(usize, f64)> = BC_2025_PATCH
            .shearing
            .iter()
            .map(|bc| (6 * bc.node + bc.dof, bc.value))
            .collect();

        let f = shearing_load(&nodes, tau, h);
        let u = solve_strong_patch(&nodes, &f, &prescribed);
        let gps = strong_form_gauss_points(&nodes, &u);

        let measure = |gps: &[StrongFormGp]| -> (f64, f64, f64, f64, f64) {
            let mut err = 0.0f64;
            let mut zero_mem = 0.0f64;
            let mut bend = 0.0f64;
            let mut shear = 0.0f64;
            let mut taus: Vec<f64> = Vec::new();
            for gp in gps {
                err = err.max((gp.membrane[2] - tau).abs());
                zero_mem = zero_mem.max(gp.membrane[0].abs()).max(gp.membrane[1].abs());
                for c in 0..3 {
                    bend = bend.max(gp.bending[c].abs());
                }
                shear = shear.max(gp.shear[0].abs()).max(gp.shear[1].abs());
                taus.push(gp.membrane[2]);
            }
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for v in &taus {
                lo = lo.min(*v);
                hi = hi.max(*v);
            }
            (err, hi - lo, zero_mem, bend, shear)
        };

        let (err, spread, zero_mem, bend, shear) = measure(&gps);
        println!(
            "T1.6c in-plane tau_xy={tau}: max|tau_gp - tau|={err:.3e} (rel {:.3e}); \
             spread={spread:.3e} (rel {:.3e}); zero membrane={zero_mem:.3e}; moment={bend:.3e}; \
             transverse shear={shear:.3e}; floor={floor:.3e}",
            err / snorm,
            spread / snorm
        );
        assert!(
            err <= 1e-8 * snorm + floor,
            "tau_xy: recovered error {err:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
            1e-8 * snorm + floor
        );
        assert!(
            spread <= 1e-8 * snorm + floor,
            "tau_xy: recovered spread {spread:.3e} > 1e-8 ||sigma|| + floor ({:.3e})",
            1e-8 * snorm + floor
        );
        assert!(
            zero_mem <= 1e-10 * snorm,
            "sigma_xx/sigma_yy: {zero_mem:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            bend <= 1e-10 * snorm,
            "moment resultant: {bend:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );
        assert!(
            shear <= 1e-10 * snorm,
            "transverse-shear resultant: {shear:.3e} > 1e-10 ||sigma|| ({:.3e})",
            1e-10 * snorm
        );

        // Non-vacuity (a): zeroing the load recovers a zero state.
        let u_zero = solve_strong_patch(&nodes, &DVector::zeros(48), &prescribed);
        let gps_zero = strong_form_gauss_points(&nodes, &u_zero);
        let (err_zero, _, _, _, _) = measure(&gps_zero);
        println!(
            "T1.6c non-vacuity (a) zeroed load: max|tau_gp - tau|={:.3e} (rel {:.3e})",
            err_zero,
            err_zero / snorm
        );
        assert!(
            err_zero / snorm > 1e-3,
            "zeroed-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err_zero / snorm
        );

        // Non-vacuity (b): the withdrawn design 4.2 derivation, a load from a
        // constant TRANSVERSE shear resultant `q = G h gamma`, does not
        // reproduce the in-plane state (Evidence gap G9). The same Fig. 7(c)
        // BC set is used, so only the load derivation differs.
        let gamma_t = 1.0e-3f64;
        let g = 2.0e11 / (2.0 * (1.0 + 0.3));
        let q = [g * gamma_t, 0.0];
        let transverse_load = boundary_integrate(&nodes, |_x, n| {
            [0.0, 0.0, q[0] * n[0] + q[1] * n[1], 0.0, 0.0, 0.0]
        });
        let u_trans = solve_strong_patch(&nodes, &transverse_load, &prescribed);
        let gps_trans = strong_form_gauss_points(&nodes, &u_trans);
        let (err_trans, _, _, _, shear_trans) = measure(&gps_trans);
        println!(
            "T1.6c non-vacuity (b) withdrawn transverse load: max|tau_gp - tau|={:.3e} \
             (rel {:.3e}); transverse shear={:.3e}",
            err_trans,
            err_trans / snorm,
            shear_trans
        );
        assert!(
            err_trans / snorm > 1e-3,
            "withdrawn transverse-load control must fail the 1e-8 tolerance, measured {:.3e}",
            err_trans / snorm
        );
    }

    // ------------------------------------------------------------------
    // 8.4 - the drill DOF's `theta_z`-free-except-corner-B device (T1.6d)
    // ------------------------------------------------------------------

    #[test]
    fn test_t1b_drill_theta_z_free_except_corner_b() {
        // Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements
        // with six degrees of freedom at each node - the MITC4/D and MITC4+/D
        // elements", Computers and Structures 308:107622, Section 3.1, p. 14:
        // "in all cases theta_z is left free at the element nodes except at the
        // corner node B ... the use of theta_z = 0 at the corner node C does
        // not affect the results."  The BC sets are Fig. 7(b)(c)(d) (p. 5).
        //
        // Spec rev 8, Requirement 11. Three variants are solved for the
        // extension, bending and shearing states, on the flat Fig. 7 star patch
        // and on the derived warped patch `STAR_NODES_WARPED`: (a) theta_z free
        // at every node except corner B; (b) theta_z = 0 additionally imposed
        // at corner C; (c) theta_z constrained at every node. Variants (a) and
        // (b) must agree to 1e-10 relative (the corner-C fixing is immaterial),
        // and variant (c) must differ from (a) by more than 1e-8 relative on at
        // least one warped patch (so the theta_z-free condition is asserted to
        // matter where it should, and the pass cannot come from an
        // over-constrained model).  The warped patch is NOT used for any
        // constant-stress assertion (it is not a constant-stress fixture).
        let h = 1.0f64;
        let sigma0 = 1.0f64;
        let kappa0 = 1.0e-3f64;
        let tau = 1.0f64;

        // The outer boundary is the flat 10x10 square for both meshes, so the
        // constant-state loads are the same.
        let ext_f = extension_load(&STAR_NODES, sigma0, h);
        let ben_f = bending_load(&STAR_NODES, kappa0);
        let she_f = shearing_load(&STAR_NODES, tau, h);
        let cases: [(&str, &DVector<f64>, &[Bc]); 3] = [
            ("extension", &ext_f, BC_2025_PATCH.extension),
            ("bending", &ben_f, BC_2025_PATCH.bending),
            ("shearing", &she_f, BC_2025_PATCH.shearing),
        ];

        // Variant (c) constrains theta_z at every node except B (B is already
        // clamped in the base set); variant (b) adds it at corner C only.
        let c_extra: Vec<(usize, f64)> = (0..8)
            .filter(|&n| n != 3)
            .map(|n| (6 * n + 5, 0.0))
            .collect();
        let b_extra: Vec<(usize, f64)> = vec![(6 * 0 + 5, 0.0)];

        // The paper's `theta_z = 0` at C immateriality is stated for the
        // strong-form patch tests (the constant-state patches); the warped
        // patch is used only for the variant-(c) separation, so the (a) == (b)
        // bound is asserted on the flat patch and the warped (a) vs (b)
        // measurement is reported as a finding.
        let mut worst_warped_c = 0.0f64;
        let mut worst_warped_ab = 0.0f64;
        for (name, f, bcs) in cases {
            let base: Vec<(usize, f64)> = bcs
                .iter()
                .map(|bc| (6 * bc.node + bc.dof, bc.value))
                .collect();
            let mut b = base.clone();
            b.extend_from_slice(&b_extra);
            let mut c = base.clone();
            c.extend_from_slice(&c_extra);

            for (mesh_name, nodes) in [("flat", &STAR_NODES), ("warped", &STAR_NODES_WARPED)] {
                let ua = solve_strong_patch(nodes, f, &base);
                let ub = solve_strong_patch(nodes, f, &b);
                let uc = solve_strong_patch(nodes, f, &c);
                let scale = max_abs_vec(&ua).max(1e-30);
                let ab = max_abs_diff_vec(&ub, &ua) / scale;
                let ca = max_abs_diff_vec(&uc, &ua) / scale;
                println!(
                    "T1.6d {name}/{mesh_name}: variant (b) vs (a) rel={ab:.3e}; \
                     variant (c) vs (a) rel={ca:.3e}"
                );
                if mesh_name == "flat" {
                    assert!(
                        ab <= 1e-10,
                        "{name}/flat: variant (b) (theta_z = 0 at C) differs from (a) by {ab:.3e} > 1e-10"
                    );
                } else {
                    worst_warped_ab = worst_warped_ab.max(ab);
                    worst_warped_c = worst_warped_c.max(ca);
                }
            }
        }
        println!(
            "T1.6d worst warped: variant (b) vs (a) rel={worst_warped_ab:.3e}; \
             variant (c) vs (a) rel={worst_warped_c:.3e}"
        );
        assert!(
            worst_warped_c > 1e-8,
            "variant (c) (theta_z constrained at every node) differs from (a) by only \
             {worst_warped_c:.3e} on the warped patch (must exceed 1e-8)"
        );
    }

    // ========================================================================
    // WU8 (S2) - the layout-bound Tier-2 tests, moved onto the new element
    //
    // These are the T2A/T2B and T2I tests of design Section 5.1. They were moved
    // out of the hybrid `mitc4.rs` test module (S2, before the flip) and
    // retargeted from `Mitc4Precomputed` to `Mitc4PlusDPrecomputed`. The
    // assertions, fields and tolerances are the originals: no bound was widened,
    // no assertion dropped and no test renamed. `make_pre()` reproduces the
    // hybrid fixture (a flat unit square, thickness 0.01, isotropic
    // E = 2.0e11, nu = 0.3, rho = 7800, k = 5/6) so the moved tests' hardcoded
    // fixture constants still describe the element they exercise.
    // ========================================================================

    /// The hybrid `mitc4.rs` fixture, rebuilt for `Mitc4PlusDPrecomputed`.
    fn make_pre() -> Mitc4PlusDPrecomputed {
        let thickness = 0.01_f64;
        let mat = IsotropicMaterial::new(2.0e11, 0.3, 7800.0);
        let shell = mat.constitutive(thickness, 5.0 / 6.0);
        let node_coords: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0];
        let pre = Mitc4PlusDPrecomputed::new(&node_coords, shell, thickness, 5.0 / 6.0);
        // Explicit retarget evidence (not vacuous): the MITC4+/D element stores
        // the single drill normal `V^D` of Ko, Bathe & Zhang (2025),
        // C&S 308:107622, Eq. (5), and the ADR-1 uncorrected transverse-shear
        // stiffness `cs / applied_k` - neither exists on `Mitc4Precomputed`.
        assert!(
            (pre.v_d.norm() - 1.0).abs() < 1e-14,
            "V^D must be a unit vector, got {}",
            pre.v_d.norm()
        );
        let cs_diff = (pre.cs_uncorrected - pre.constitutive.cs / (5.0 / 6.0)).norm();
        assert!(
            cs_diff < 1e-12 * pre.constitutive.cs.norm(),
            "cs_uncorrected must be cs / applied_k, diff = {cs_diff}"
        );
        pre
    }

    #[test]
    fn test_ke_local_flat_plate_parity() {
        let pre = make_pre();
        let ke = compute_ke_local(&pre);

        // We expect Ke to be symmetric
        let diff = &ke - ke.transpose();
        assert!(diff.norm() < 1e-10, "Ke must be symmetric");

        // For a unit square, the membrane part of Ke should be non-zero
        assert!(ke.norm() > 1e-6, "Ke should not be zero");
    }

    /// Element centroid from the initial global coordinates.
    fn element_centroid(pre: &Mitc4PlusDPrecomputed) -> Vector3<f64> {
        let mut c = Vector3::zeros();
        for i in 0..4 {
            c += Vector3::new(
                pre.initial_coords_3d[i][0],
                pre.initial_coords_3d[i][1],
                pre.initial_coords_3d[i][2],
            );
        }
        c / 4.0
    }

    /// Physical rigid-body field: ``u = t + omega x (x - x_c)``, ``theta = omega``.
    fn rigid_body_mode(
        pre: &Mitc4PlusDPrecomputed,
        translation: Vector3<f64>,
        omega: Vector3<f64>,
    ) -> Vec24 {
        let centroid = element_centroid(pre);
        let mut u = Vec24::zeros();
        for i in 0..4 {
            let x = Vector3::new(
                pre.initial_coords_3d[i][0],
                pre.initial_coords_3d[i][1],
                pre.initial_coords_3d[i][2],
            );
            let disp = translation + omega.cross(&(x - centroid));
            for k in 0..3 {
                u[6 * i + k] = disp[k];
                u[6 * i + 3 + k] = omega[k];
            }
        }
        u
    }

    /// The six *labelled* physical rigid-body fields (helper name differs from
    /// the hybrid's `rigid_body_modes` because the Tier-1 BC fixture already
    /// owns that name in this module; the fields and the assertions are the
    /// hybrid's).
    fn physical_rigid_body_modes(pre: &Mitc4PlusDPrecomputed) -> [(&'static str, Vec24); 6] {
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
    fn scaled_residual(k: &Mat24, u: &Vec24) -> f64 {
        let denom = k.norm() * u.norm();
        if denom > 0.0 {
            (k * u).norm() / denom
        } else {
            0.0
        }
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
        let lambda_max = eigenvalues
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
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
        for (label, u) in physical_rigid_body_modes(&pre) {
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
        let mut u = Vec24::zeros();
        for i in 0..4 {
            let x = pre.initial_coords_3d[i][0];
            let y = pre.initial_coords_3d[i][1];
            u[6 * i] = a * x + b * y;
            u[6 * i + 1] = c * x + d * y;
        }

        let expected = Vector3::new(a, d, b + c);
        for g in 0..N_GAUSS {
            let eps = b_membrane_2017(&pre, GAUSS_XI[g], GAUSS_ETA[g]) * u;
            let error = (eps - expected).norm() / expected.norm();
            // 1e-10 relative is far above the measured round-off of the
            // Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (27) operator
            // (~4e-12 for this field, i.e. ~4e-15 absolute on a 1e-3 strain)
            // and far below any real formulation error, which would be O(1).
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
        let mut u = Vec24::zeros();
        for i in 0..4 {
            let x = pre.initial_coords_3d[i][0];
            u[6 * i + 2] = 0.5 * kxx * x * x;
            u[6 * i + 4] = kxx * x;
        }

        let expected = Vector3::new(kxx, 0.0, 0.0);
        for g in 0..N_GAUSS {
            // `b_bending_2017` returns the paper's t-linear bending strain
            // measure e^b1 (Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (7a),
            // with the through-thickness coordinate t = 2z/h), so the physical
            // curvature is kappa = (2/h) e^b1. The invariant and its 1e-10
            // relative bound are the hybrid's; only the operator's
            // normalization is expressed.
            let (b1, _b2) = b_bending_2017(&pre, GAUSS_XI[g], GAUSS_ETA[g]);
            let kappa = (2.0 / pre.thickness) * (b1 * u);
            let error = (kappa - expected).norm() / expected.norm();
            assert!(
                error < 1e-10,
                "bending patch test at Gauss point {g}: kappa = {kappa:?}, \
                 expected = {expected:?}, relative error = {error:.3e}"
            );
        }
    }

    // ------------------------------------------------------------------
    // T2I - the body-load, initial-stress and stress-recovery entry points
    // ------------------------------------------------------------------

    #[test]
    fn test_body_load_global_zero_gravity() {
        let pre = make_pre();
        let g = Vector3::zeros();
        let f = compute_body_load_global(&pre, 7800.0, &g);
        assert!(f.norm() < 1e-12, "zero gravity -> zero body load");
    }

    #[test]
    fn test_body_load_global_z_gravity() {
        let pre = make_pre();
        let g = Vector3::new(0.0, 0.0, -9.81);
        let f = compute_body_load_global(&pre, 7800.0, &g);

        // Only translational z-DOFs should be non-zero
        for i in 0..4 {
            assert!(f[6 * i].abs() < 1e-10, "node {i} fx should be ~0");
            assert!(f[6 * i + 1].abs() < 1e-10, "node {i} fy should be ~0");
            assert!(f[6 * i + 2].abs() > 1e-6, "node {i} fz should be nonzero");
            for k in 3..6 {
                assert!(
                    f[6 * i + k].abs() < 1e-12,
                    "node {i} rotational dof {k} should be 0"
                );
            }
        }

        // Total z-force = rho h |g| area (area = 1.0 for unit square)
        let area = 1.0_f64;
        let h = 0.01_f64;
        let rho = 7800.0_f64;
        let expected_total_fz = rho * h * (-9.81) * area;
        let total_fz: f64 = (0..4).map(|i| f[6 * i + 2]).sum();
        assert!(
            (total_fz - expected_total_fz).abs() < 1e-4,
            "total fz: got {total_fz}, expected {expected_total_fz}"
        );
    }

    #[test]
    fn test_k_sigma_global_zero_stress() {
        let pre = make_pre();
        let sigma = Vector3::zeros();
        let k = compute_k_sigma_global(&pre, &sigma);
        assert!(k.norm() < 1e-12, "zero stress -> zero K_sigma");
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
        let k_local = geometric_stiffness_from_stress(&pre, &sigma);
        let k_global_direct = compute_k_sigma_global(&pre, &sigma);
        let k_global_manual = transform_to_global(&pre, &k_local);
        let diff = k_global_direct - k_global_manual;
        assert!(
            diff.norm() < 1e-6,
            "compute_k_sigma_global must equal transform(k_local)"
        );
    }

    #[test]
    fn test_centrifugal_prestress_on_axis() {
        let pre = make_pre();
        let axis = Vector3::new(0.0, 0.0, 1.0);
        // Place center at the centroid (0.5, 0.5, 0) -> r_radial ~ 0
        let center = Vector3::new(0.5, 0.5, 0.0);
        let centroid = Vector3::new(0.5, 0.5, 0.0);
        let sigma = compute_centrifugal_prestress(&pre, 100.0, &axis, &center, &centroid, 7800.0);
        assert!(
            sigma.norm() < 1e-6,
            "element on axis -> zero centrifugal stress"
        );
    }

    #[test]
    fn test_centrifugal_prestress_nonzero() {
        let pre = make_pre();
        let axis = Vector3::new(0.0, 0.0, 1.0);
        let center = Vector3::zeros();
        let centroid = Vector3::new(0.5, 0.5, 0.0);
        let sigma = compute_centrifugal_prestress(&pre, 100.0, &axis, &center, &centroid, 7800.0);
        assert!(sigma[0].is_finite());
        assert!(sigma[1].is_finite());
        assert!(sigma[2].is_finite());
        assert!(
            sigma[0] + sigma[1] >= 0.0,
            "centrifugal stress trace must be non-negative"
        );
    }

    // ========================================================================
    // WU9i — block isolation of defect #4
    // (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// WU9i — block isolation of defect #4 (MEASUREMENT INSTRUMENT, NOT A GATE).
    ///
    /// Defect #4 (`openspec/changes/mitc4plusd-faithful/fidelity-audit.md`): the
    /// flat in-plane bending of the 8x4 cantilever strip is ~37% too stiff. The
    /// WU9f flip measured `uY/uX = 251.55` against the beam-theory 400.00 while
    /// the hybrid is 0.04% off. On a flat element the transverse-shear block and
    /// the bending block have non-zero entries only on {u_z, theta_x, theta_y},
    /// so for the symmetric in-plane solution the flat in-plane stiffness is
    /// exactly `membrane + drill`. `test_fx` (axial) passes while `test_fy` and
    /// `test_ratio_physical` fail, so the defect is in the GRADIENT response.
    ///
    /// This instrument SOLVES the exact failing case densely (mirrors
    /// `tests/test_shell_validation_fixed.py::_build_cantilever_mesh`, without
    /// importing it) and decomposes the FULL solution's energy by block. It
    /// MEASURES; it does not confirm a hypothesis and it fixes nothing. Run it
    /// with `cargo test -p aeroelast-core wu9i -- --ignored --nocapture`.
    ///
    /// RESOLVED 2026-09-24: defect #4 was the missing Eq. (22a) coupling — the
    /// drill was added as an independent energy block instead of being summed
    /// into the membrane strain. With the fold the FULL ratio is 400.3877
    /// (+0.10% vs beam theory) and the cross terms appear in the shares; the
    /// pre-fold value 251.5515 is kept in the printed output as the recorded
    /// baseline.
    #[test]
    #[ignore = "WU9i measurement instrument (defect #4 block isolation); run with --ignored --nocapture"]
    fn wu9i_block_isolation_flat_inplane_strip() {
        // ---- the exact failing case: 8x4 cantilever strip, L=1, b=0.1, h=1 mm ----
        const L: f64 = 1.0;
        const WIDTH: f64 = 0.1;
        const H: f64 = 0.001;
        const E: f64 = 2.1e11;
        const NU: f64 = 0.3;
        const NX: usize = 8;
        const NY: usize = 4;
        const NNODES: usize = (NX + 1) * (NY + 1); // 45
        const NDOF: usize = 6 * NNODES; // 270
        const NELEM: usize = NX * NY; // 32
        const SHEAR_CORRECTION: f64 = 5.0 / 6.0;
        const P_TOTAL: f64 = 600.0;
        const P_PER_NODE: f64 = P_TOTAL / (NY as f64 + 1.0); // 120 N on each free-edge node

        // ---- element-local stiffness builders (nested fn: no captures) ----
        fn symmetrise(k: &Mat24) -> Mat24 {
            let mut out = *k;
            for i in 0..24 {
                for j in 0..24 {
                    out[(i, j)] = 0.5 * (k[(i, j)] + k[(j, i)]);
                }
            }
            out
        }

        /// The production 9-row `[bm; s1 bb1; s2 bb2]^T W [..]` contribution of
        /// `compute_ke_local_with_drill` (identical recipe, symmetrised).
        fn ke_bend9(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            let w = resultant_moment_matrix(&pre.constitutive, pre.thickness);
            let s1 = 2.0 / pre.thickness;
            let s2 = 4.0 / (pre.thickness * pre.thickness);
            let mut k = Mat24::zeros();
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let sqrt_g = surface_measure(pre, r, s);
                // Eq. (22a): the same folded `t^0` row the production uses.
                let mut bm = b_membrane_2017(pre, r, s);
                bm += b_drill_membrane_2025(pre, r, s);
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
            symmetrise(&k)
        }

        /// Variant 3/4: ONLY the membrane field changes. The covariant rows are
        /// the COMPATIBLE (displacement-based) ones of
        /// `covariant_membrane_b_row` (`comp = 0,1,2`, row 2 doubled as the
        /// engineering shear `2 e_rs`), mapped by the production point-wise
        /// covariant-to-local map (`covariant_to_local_mapping(j_loc_at(..))`,
        /// numerically identical to `ke_ref::map_local`). Bending, shear and
        /// drill are the production blocks, untouched.
        fn ke_compat(pre: &Mitc4PlusDPrecomputed, use_drill: bool) -> Mat24 {
            let w = resultant_moment_matrix(&pre.constitutive, pre.thickness);
            let s1 = 2.0 / pre.thickness;
            let s2 = 4.0 / (pre.thickness * pre.thickness);
            let mut k = Mat24::zeros();
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let sqrt_g = surface_measure(pre, r, s);
                let map = covariant_to_local_mapping(&j_loc_at(pre, r, s));
                let mut cov = SMatrix::<f64, 3, 24>::zeros();
                for comp in 0..3 {
                    let row = super::covariant_membrane_b_row(
                        &pre.x_r, &pre.x_s, &pre.x_d, &pre.e1, &pre.e2, &pre.e3, r, s, comp,
                    );
                    let factor = if comp == 2 { 2.0 } else { 1.0 };
                    for j in 0..24 {
                        cov[(comp, j)] = factor * row[j];
                    }
                }
                let mut bm = map * cov;
                if use_drill {
                    // Eq. (22a): the drill-membrane strain is summed into the
                    // membrane row, exactly as the production now does.
                    bm += b_drill_membrane_2025(pre, r, s);
                }
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
            let mut out = symmetrise(&k);
            out += shear_ke_local(pre);
            symmetrise(&out)
        }

        fn ke_full(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            compute_ke_local_with_drill(pre, true)
        }
        fn ke_no_drill(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            compute_ke_local_with_drill(pre, false)
        }
        fn ke_compat_full(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            ke_compat(pre, true)
        }
        fn ke_compat_no_drill(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            ke_compat(pre, false)
        }
        fn ke_memb(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            membrane_ke_local(pre)
        }
        fn ke_shear(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            shear_ke_local(pre)
        }
        fn ke_drill(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            drill_ke_local(pre)
        }

        /// The membrane block alone, but with the COMPATIBLE covariant rows in
        /// place of the assumed MITC4+ field (used to show that on this flat,
        /// undistorted mesh the two fields coincide).
        fn compat_membrane_ke(pre: &Mitc4PlusDPrecomputed) -> Mat24 {
            let cm = &pre.constitutive.cm;
            let mut k = Mat24::zeros();
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let sqrt_g = surface_measure(pre, r, s);
                let map = covariant_to_local_mapping(&j_loc_at(pre, r, s));
                let mut cov = SMatrix::<f64, 3, 24>::zeros();
                for comp in 0..3 {
                    let row = super::covariant_membrane_b_row(
                        &pre.x_r, &pre.x_s, &pre.x_d, &pre.e1, &pre.e2, &pre.e3, r, s, comp,
                    );
                    let factor = if comp == 2 { 2.0 } else { 1.0 };
                    for j in 0..24 {
                        cov[(comp, j)] = factor * row[j];
                    }
                }
                let bm = map * cov;
                k += (bm.transpose() * cm * bm) * (GAUSS_W[g] * sqrt_g);
            }
            k
        }

        // ---- mesh (node index = j*(NX+1) + i, exactly as the reference mesh) ----
        let nidx = |i: usize, j: usize| j * (NX + 1) + i;
        let gdof = |n: usize, k: usize| 6 * n + k;
        let centre = nidx(NX, NY / 2); // free-edge centre node (x = L, y = b/2)
        let node_xyz: Vec<[f64; 3]> = (0..NNODES)
            .map(|n| {
                let (i, j) = (n % (NX + 1), n / (NX + 1));
                [
                    L * (i as f64) / (NX as f64),
                    WIDTH * (j as f64) / (NY as f64),
                    0.0,
                ]
            })
            .collect();
        let mut elem_coords: Vec<[[f64; 3]; 4]> = Vec::with_capacity(NELEM);
        let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(NELEM);
        for j in 0..NY {
            for i in 0..NX {
                let n = [nidx(i, j), nidx(i + 1, j), nidx(i + 1, j + 1), nidx(i, j + 1)];
                elem_coords.push([node_xyz[n[0]], node_xyz[n[1]], node_xyz[n[2]], node_xyz[n[3]]]);
                elem_nodes.push(n);
            }
        }

        let build_pre = |c: &[[f64; 3]; 4]| {
            let constitutive =
                IsotropicMaterial::new(E, NU, 7800.0).constitutive(H, SHEAR_CORRECTION);
            Mitc4PlusDPrecomputed::new(&coords12(c), constitutive, H, SHEAR_CORRECTION)
        };

        // ---- global assembly `K = sum_e T24^T K_e T24` ----
        let assemble_global = |ke_of: fn(&Mitc4PlusDPrecomputed) -> Mat24| -> DMatrix<f64> {
            let mut kg = DMatrix::zeros(NDOF, NDOF);
            for e in 0..NELEM {
                let pre = build_pre(&elem_coords[e]);
                let ke_local = ke_of(&pre);
                let t24 = build_t24(&pre);
                let ke = t24.transpose() * ke_local * t24;
                for a in 0..4 {
                    for b in 0..4 {
                        for ka in 0..6 {
                            for kb in 0..6 {
                                kg[(gdof(elem_nodes[e][a], ka), gdof(elem_nodes[e][b], kb))] +=
                                    ke[(6 * a + ka, 6 * b + kb)];
                            }
                        }
                    }
                }
            }
            kg
        };

        // ---- dense reduced solve (tiny: <= 240x240) ----
        let solve_reduced = |kg: &DMatrix<f64>, f: &DVector<f64>, fixed: &[usize]| -> Option<DVector<f64>> {
            let n = kg.nrows();
            let mut is_fixed = vec![false; n];
            for &d in fixed {
                is_fixed[d] = true;
            }
            let free: Vec<usize> = (0..n).filter(|&i| !is_fixed[i]).collect();
            let m = free.len();
            let mut kff = DMatrix::zeros(m, m);
            for (a, &i) in free.iter().enumerate() {
                for (b, &jj) in free.iter().enumerate() {
                    kff[(a, b)] = kg[(i, jj)];
                }
            }
            let mut rhs = DVector::zeros(m);
            for (a, &i) in free.iter().enumerate() {
                rhs[a] = f[i];
            }
            let sol = kff.lu().solve(&rhs)?;
            let mut u = DVector::zeros(n);
            for (a, &i) in free.iter().enumerate() {
                u[i] = sol[a];
            }
            Some(u)
        };

        // ---- singularity evidence: spectrum of the constrained K ----
        let reduced_spectrum = |kg: &DMatrix<f64>, fixed: &[usize]| -> (Vec<f64>, f64) {
            let n = kg.nrows();
            let mut is_fixed = vec![false; n];
            for &d in fixed {
                is_fixed[d] = true;
            }
            let free: Vec<usize> = (0..n).filter(|&i| !is_fixed[i]).collect();
            let m = free.len();
            let mut red = DMatrix::zeros(m, m);
            for (a, &i) in free.iter().enumerate() {
                for (b, &j) in free.iter().enumerate() {
                    red[(a, b)] = kg[(i, j)];
                }
            }
            let sym = (&red + red.transpose()) * 0.5;
            let mut eig: Vec<f64> = sym.symmetric_eigenvalues().iter().cloned().collect();
            eig.sort_by(|a, b| a.total_cmp(b));
            let lam_max = eig.iter().fold(0.0f64, |mx, &v| mx.max(v.abs()));
            (eig, lam_max)
        };

        // ---- boundary conditions and loads ----
        let mut clamped: Vec<usize> = Vec::new();
        for j in 0..=NY {
            for k in 0..6 {
                clamped.push(gdof(nidx(0, j), k)); // every DOF of the x = 0 nodes
            }
        }
        let mut theta_z: Vec<usize> = Vec::new();
        for n in 0..NNODES {
            theta_z.push(gdof(n, 5));
        }
        // Every DOF other than the in-plane translations (u_x, u_y) at every
        // node: the membrane block has no stiffness for the remaining DOFs on a
        // flat element, so this isolates the purely in-plane subsystem.
        let mut non_inplane: Vec<usize> = Vec::new();
        for n in 0..NNODES {
            for k in 2..6 {
                non_inplane.push(gdof(n, k));
            }
        }
        let load_vector = |dir: usize| -> DVector<f64> {
            let mut f = DVector::zeros(NDOF);
            for j in 0..=NY {
                f[gdof(nidx(NX, j), dir)] += P_PER_NODE;
            }
            f
        };
        let f_x = load_vector(0);
        let f_y = load_vector(1);

        // ---- one variant: assemble, count near-zero modes, solve both loads ----
        let run = |label: &str,
                   ke_of: fn(&Mitc4PlusDPrecomputed) -> Mat24,
                   extra_fixed: &[usize]|
         -> (DMatrix<f64>, Option<DVector<f64>>, Option<DVector<f64>>) {
            let kg = assemble_global(ke_of);
            let mut fixed = clamped.clone();
            fixed.extend_from_slice(extra_fixed);
            fixed.sort_unstable();
            fixed.dedup();
            let reduced = NDOF - fixed.len();
            let (eig, lam_max) = reduced_spectrum(&kg, &fixed);
            let n_zero_8 = eig.iter().filter(|&&v| v.abs() < 1e-8 * lam_max).count();
            let n_zero_12 = eig.iter().filter(|&&v| v.abs() < 1e-12 * lam_max).count();
            let smallest: Vec<String> = eig.iter().take(6).map(|v| format!("{v:.4e}")).collect();
            let sol_x = solve_reduced(&kg, &f_x, &fixed);
            let sol_y = solve_reduced(&kg, &f_y, &fixed);
            match (&sol_x, &sol_y) {
                (Some(ux), Some(uy)) => {
                    let (rx, ry) = (ux[gdof(centre, 0)], uy[gdof(centre, 1)]);
                    println!(
                        "WU9i {label:<31} ux = {rx:.12e} m   uy = {ry:.12e} m   ratio uy/ux = {:>9.4}   [reduced {reduced} dof, lam_max {lam_max:.4e}, near-zero(1e-8) {n_zero_8}, near-zero(1e-12) {n_zero_12}, smallest eps = {}]",
                        ry / rx,
                        smallest.join(", ")
                    );
                }
                _ => println!(
                    "WU9i {label:<31} SINGULAR: dense LU reports no solution for the reduced {reduced}-dof constrained K; near-zero(1e-8) {n_zero_8}, near-zero(1e-12) {n_zero_12} of lam_max = {lam_max:.4e}; smallest eps = {}",
                    smallest.join(", ")
                ),
            }
            (kg, sol_x, sol_y)
        };

        let ux_ref = P_TOTAL * L / (E * WIDTH * H);
        println!(
            "\n=== WU9i: flat in-plane block isolation (8x4 cantilever, L={L}, b={WIDTH}, h={H}, E={E:.1e}, nu={NU}) ==="
        );
        println!(
            "nodes = {NNODES}, dofs = {NDOF}, elements = {NELEM}, clamped dofs = {}, analytical ux = 600 L/(E b h) = {ux_ref:.9e} m",
            clamped.len()
        );

        let (k_full, sol_x_full, sol_y_full) = run("1. FULL", ke_full, &[]);
        let (_k_nd, _sx_nd, _sy_nd) = run("2. NO_DRILL", ke_no_drill, &[]);
        let (_k_ndt, _sx_ndt, _sy_ndt) = run("2b. NO_DRILL + theta_z=0", ke_no_drill, &theta_z);
        let (_k_c, _sx_c, _sy_c) = run("3. COMPAT_MEMB", ke_compat_full, &[]);
        let (_k_cn, _sx_cn, _sy_cn) = run("4. COMPAT_MEMB_NO_DRILL", ke_compat_no_drill, &[]);
        let (_k_cnt, _sx_cnt, _sy_cnt) =
            run("4b. COMPAT_MEMB_NO_DRILL + theta_z=0", ke_compat_no_drill, &theta_z);
        let (_k_mo, _sx_mo, _sy_mo) = run("5. MEMBRANE_ONLY (in-plane only)", ke_memb, &non_inplane);

        // Is the assumed MITC4+ membrane field different from the compatible
        // one on this flat, undistorted mesh? (element-wise, first element)
        let pre0 = build_pre(&elem_coords[0]);
        let km_assumed = membrane_ke_local(&pre0);
        let km_compat = compat_membrane_ke(&pre0);
        println!(
            "WU9i membrane block, assumed MITC4+ vs compatible (element 0): max|K_assumed - K_compat| = {:.3e}, max|K_assumed| = {:.4e}, relative {:.3e}",
            max_abs_diff(&km_assumed, &km_compat),
            max_abs(&km_assumed),
            max_abs_diff(&km_assumed, &km_compat) / max_abs(&km_assumed).max(f64::MIN_POSITIVE)
        );

        let (ux_full, uy_full) = match (&sol_x_full, &sol_y_full) {
            (Some(ux), Some(uy)) => (ux.clone(), uy.clone()),
            _ => panic!("WU9i: the FULL variant must be non-singular (the drill block regularises it)"),
        };
        let ratio_full = uy_full[gdof(centre, 1)] / ux_full[gdof(centre, 0)];
        let ux_full_centre = ux_full[gdof(centre, 0)];
        println!(
            "WU9i instrument validation: FULL ratio = {ratio_full:.4} vs beam theory 400.0 -> {:.2}% (the pre-Eq.(22a) separate-block value was 251.5515, -37.1%); FULL ux = {ux_full_centre:.9e} vs analytical {ux_ref:.9e} -> {:.2}%",
            100.0 * (ratio_full - 400.0).abs() / 400.0,
            100.0 * (ux_full_centre - ux_ref).abs() / ux_ref
        );

        // ---- block energy decomposition on the FULL solution ----
        let k_memb = assemble_global(ke_memb);
        let k_shear = assemble_global(ke_shear);
        let k_drill = assemble_global(ke_drill);
        // The bending+coupling share is taken BY DIFFERENCE: with Eq. (22a)
        // `K_full = K_9row(folded) + K_shear`, and the folded `t^0` row is
        // `B_m + B_md`, so subtracting the pure membrane and pure drill blocks
        // leaves the paper's cross terms plus the bending rows.
        let k_bend = &k_full - &k_memb - &k_shear - &k_drill;
        // Cross-check against the explicit folded 9-row recipe minus the pure
        // membrane and the pure drill blocks: that residual is
        // `k_cross + k_bending`, exactly what `k_bend` is by difference.
        let k9 = assemble_global(ke_bend9);
        let k_bend_recipe = &k9 - &k_memb - &k_drill;
        let mut max_diff = 0.0f64;
        let mut max_ref = 0.0f64;
        for i in 0..NDOF {
            for j in 0..NDOF {
                max_diff = max_diff.max((k_bend[(i, j)] - k_bend_recipe[(i, j)]).abs());
                max_ref = max_ref.max(k_bend_recipe[(i, j)].abs());
            }
        }
        println!(
            "WU9i bending share route: by difference == ([bm; s1 bb1; s2 bb2]^T W [..] - K_membrane) to {max_diff:.3e} abs / {:.3e} rel (max|K_bend| = {max_ref:.4e})",
            max_diff / max_ref.max(f64::MIN_POSITIVE)
        );

        let energy = |k: &DMatrix<f64>, u: &DVector<f64>| -> f64 {
            let mut s = 0.0;
            for i in 0..k.nrows() {
                let mut row = 0.0;
                for j in 0..k.ncols() {
                    row += k[(i, j)] * u[j];
                }
                s += u[i] * row;
            }
            s
        };

        let mut shares_all: Vec<[f64; 4]> = Vec::new();
        for (load_label, u) in [("+x (axial)", &ux_full), ("+y (in-plane bending)", &uy_full)] {
            let e_total = energy(&k_full, u);
            let shares = [
                energy(&k_memb, u) / e_total,
                energy(&k_bend, u) / e_total,
                energy(&k_shear, u) / e_total,
                energy(&k_drill, u) / e_total,
            ];
            println!(
                "WU9i block energy shares (FULL solution, load {load_label}): membrane {:.12}, bending {:.12}, shear {:.12}, drill {:.12} | sum {:.16} (u^T K u = {e_total:.6e})",
                shares[0],
                shares[1],
                shares[2],
                shares[3],
                shares.iter().sum::<f64>()
            );
            shares_all.push(shares);
        }

        // (a) the four shares must be finite and partition u^T K u.
        for shares in &shares_all {
            let total: f64 = shares.iter().sum();
            assert!(
                shares.iter().all(|s| s.is_finite()),
                "block energy shares must be finite: {shares:?}"
            );
            assert!(
                (total - 1.0).abs() <= 1e-10,
                "block shares must sum to 1 within 1e-10: sum = {total:.17}"
            );
        }

        // (b) The faithful (Eq. (22a)) element must match the beam-theory
        // `4 L^2 / b^2 = 400` within the repo test's 2% window. Before the fold
        // the same instrument measured 251.5515 (-37.1%); that value is
        // recorded in the output above, not asserted.
        assert!(
            (ratio_full - 400.0).abs() <= 0.02 * 400.0,
            "WU9i instrument validation: FULL ratio {ratio_full:.4} is not within 2% of the beam-theory 400.0"
        );
    }

    // ========================================================================
    // t2025 — the paper's own Table 1 slender plane-stress cantilever
    // (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// Ko, Bathe & Zhang (2025), C&S 308:107622, §3.2 (Fig. 8, mesh type 1) and
    /// Table 1, the "vertical deflection (`-u_y`) at the tip (average of
    /// displacements at nodes near point A) for the slender cantilever
    /// problems" (MEASUREMENT INSTRUMENT, NOT A GATE).
    ///
    /// The cell reproduced here is the flat plane-stress **MITC4+/D** column for
    /// **mesh type 1** (the regular mesh): the paper publishes `0.0976755`,
    /// against its own reference solution `0.1081`. The other mesh-type-1 cells
    /// are MITC4 `0.010088` and MITC4-IC `0.107328`.
    ///
    /// Geometry (Fig. 8a, self-consistent with the reference value):
    /// `L = 6.0` (x), in-plane height `h = 0.2` (y), out-of-plane thickness
    /// `t = 0.1`, `E = 1.0e7`, `nu = 0.3`, tip load `P = 1.0` in `-y`. The
    /// Euler-Bernoulli tip deflection `P L^3/(3 E I)` with
    /// `I = t h^3/12 = 6.6667e-5` is `0.108`, i.e. the paper's reference. Mesh
    /// type 1 is the regular `6 x 1` quad mesh: 14 nodes, 84 DOF, 6 elements.
    ///
    /// Why this exists: `defect #4` (the missing Eq. (22a) coupling, fixed) was
    /// measured with `wu9i_block_isolation_flat_inplane_strip`, whose strip has
    /// a different aspect ratio and no published counterpart. This instrument is
    /// the harder sibling: slender (element aspect ratio 5), dominated by the
    /// in-plane bending response, and it compares our element against the
    /// paper's **own published numbers** rather than against beam theory alone.
    ///
    /// It MEASURES; it confirms nothing and it changes no production code. The
    /// assembly reuses the repository's own global element matrix
    /// (`compute_ke_global`, i.e. `T^T K_e T`) and the module's dense
    /// constrained solver (`solve_constrained`); only the mesh scatter is local.
    ///
    /// Run with:
    /// `cargo test -p aeroelast-core t2025_table1 -- --ignored --nocapture`.
    #[test]
    #[ignore = "2025 paper Table 1 measurement instrument; run with --ignored --nocapture"]
    fn t2025_table1_slender_plane_stress_cantilever() {
        // ---- published cells (Ko, Bathe & Zhang (2025), C&S 308:107622, Table 1,
        //      mesh type 1 = regular mesh; see the fidelity audit's citation). ----
        const PUB_MITC4: f64 = 0.010088;
        const PUB_MITC4_IC: f64 = 0.107328;
        const PUB_MITC4_D: f64 = 0.0976755; // the column our element must match (or not)
        const PUB_REFERENCE: f64 = 0.1081;
        const BEND_THEORY: f64 = 0.108; // P L^3 / (3 E I)
        // Deliberately broad physical bracket for the instrument: from the
        // paper's MITC4 cell up to a little above its reference solution.
        const BRACKET_LO: f64 = PUB_MITC4;
        const BRACKET_HI: f64 = 0.115;

        // ---- Fig. 8a geometry, material and load ----
        const L: f64 = 6.0;
        const HEIGHT: f64 = 0.2;
        const THICKNESS: f64 = 0.1;
        const E: f64 = 1.0e7;
        const NU: f64 = 0.3;
        const SHEAR_CORRECTION: f64 = 5.0 / 6.0;
        const P_TOTAL: f64 = 1.0; // downward (-y) at the free edge
        // ---- mesh type 1: the regular `NX x NY` quad mesh ----
        const NX: usize = 6;
        const NY: usize = 1;
        const NNODES: usize = (NX + 1) * (NY + 1); // 14
        const NDOF: usize = 6 * NNODES; // 84
        const NELEM: usize = NX * NY; // 6

        let nidx = |i: usize, j: usize| j * (NX + 1) + i;
        let gdof = |n: usize, k: usize| 6 * n + k;

        let node_xyz: Vec<[f64; 3]> = (0..NNODES)
            .map(|n| {
                let (i, j) = (n % (NX + 1), n / (NX + 1));
                [
                    L * (i as f64) / (NX as f64),
                    HEIGHT * (j as f64) / (NY as f64),
                    0.0,
                ]
            })
            .collect();
        let mut elem_coords: Vec<[[f64; 3]; 4]> = Vec::with_capacity(NELEM);
        let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(NELEM);
        for j in 0..NY {
            for i in 0..NX {
                let n = [nidx(i, j), nidx(i + 1, j), nidx(i + 1, j + 1), nidx(i, j + 1)];
                elem_coords.push([
                    node_xyz[n[0]],
                    node_xyz[n[1]],
                    node_xyz[n[2]],
                    node_xyz[n[3]],
                ]);
                elem_nodes.push(n);
            }
        }

        // ---- global assembly `K = sum_e T^T K_e T` (production global matrix) ----
        let mut kg = DMatrix::zeros(NDOF, NDOF);
        for e in 0..NELEM {
            let constitutive =
                IsotropicMaterial::new(E, NU, 7800.0).constitutive(THICKNESS, SHEAR_CORRECTION);
            let pre = Mitc4PlusDPrecomputed::new(
                &coords12(&elem_coords[e]),
                constitutive,
                THICKNESS,
                SHEAR_CORRECTION,
            );
            let ke = compute_ke_global(&pre);
            for a in 0..4 {
                for b in 0..4 {
                    for ka in 0..6 {
                        for kb in 0..6 {
                            kg[(gdof(elem_nodes[e][a], ka), gdof(elem_nodes[e][b], kb))] +=
                                ke[(6 * a + ka, 6 * b + kb)];
                        }
                    }
                }
            }
        }

        // ---- clamp all 6 DOF at `x = 0` (the paper's clamped edge) ----
        let mut prescribed: Vec<(usize, f64)> = Vec::new();
        for j in 0..=NY {
            for k in 0..6 {
                prescribed.push((gdof(nidx(0, j), k), 0.0));
            }
        }

        // ---- tip load: `P` in `-y`, one nodal load per free-edge node ----
        let p_per_node = P_TOTAL / ((NY + 1) as f64);
        let mut f: DVector<f64> = DVector::zeros(NDOF);
        for j in 0..=NY {
            f[gdof(nidx(NX, j), 1)] -= p_per_node;
        }
        let applied_total: f64 = (0..NDOF).map(|i| -f[i]).sum();

        // ---- solve with the module's dense constrained solver ----
        let u = solve_constrained(&kg, &f, &prescribed);

        // ---- measure: `-u_y` at A = average of `u_y` over the free-edge nodes ----
        let uy_avg: f64 =
            (0..=NY).map(|j| u[gdof(nidx(NX, j), 1)]).sum::<f64>() / ((NY + 1) as f64);
        let tip = -uy_avg;

        let r_pub = tip / PUB_MITC4_D;
        let r_ref = tip / PUB_REFERENCE;
        let r_eb = tip / BEND_THEORY;
        let pos_in_bracket = (tip - BRACKET_LO) / (BRACKET_HI - BRACKET_LO);

        println!("\n=== 2025 Table 1 instrument: slender plane-stress cantilever, mesh type 1 ===");
        println!(
            "geometry: L = {L}, h = {HEIGHT}, t = {THICKNESS}, E = {E:.1e}, nu = {NU}, shear correction {SHEAR_CORRECTION:.6}"
        );
        println!(
            "mesh: {NX}x{NY} regular quads, {NNODES} nodes, {NDOF} dof, {NELEM} elements, {} clamped dof (all 6 at x = 0); load P = {P_TOTAL} in -y -> {p_per_node} per free-edge node (applied total {applied_total})",
            prescribed.len()
        );
        println!(
            "measured -u_y at A (average over the {}-node free edge) = {tip:.12e}",
            NY + 1
        );
        println!("published 2025 Table 1, mesh type 1 (ours / published):");
        println!(
            "  MITC4     = {PUB_MITC4:.6}  -> {:.6}",
            tip / PUB_MITC4
        );
        println!(
            "  MITC4-IC  = {PUB_MITC4_IC:.6}  -> {:.6}",
            tip / PUB_MITC4_IC
        );
        println!("  MITC4+/D  = {PUB_MITC4_D:.7}  -> {r_pub:.9}  (difference {:+.3e})", tip - PUB_MITC4_D);
        println!("  reference = {PUB_REFERENCE:.4}  -> {r_ref:.6}");
        println!("  Euler-Bernoulli P L^3/(3 E I) = {BEND_THEORY:.4}  -> {r_eb:.6}");
        println!(
            "bracket [{BRACKET_LO:.6}, {BRACKET_HI}]: measured value sits {:.3}% into it",
            100.0 * pos_in_bracket
        );
        println!(
            "reading: vs published MITC4+/D -> {:+.2}%; vs reference -> {:+.2}%; vs beam theory -> {:+.2}%",
            100.0 * (r_pub - 1.0),
            100.0 * (r_ref - 1.0),
            100.0 * (r_eb - 1.0)
        );

        // ---- loose instrument assertions (the parent draws the verdict) ----
        assert!(
            tip.is_finite(),
            "the measured tip deflection must be finite, got {tip}"
        );
        assert!(
            tip >= BRACKET_LO && tip <= BRACKET_HI,
            "the measured tip deflection {tip:.9e} must lie inside the broad physical bracket [{BRACKET_LO}, {BRACKET_HI}]"
        );
        // Non-vacuity: the load must be the one that produces the measured value.
        assert!(
            (applied_total - P_TOTAL).abs() < 1e-15,
            "the applied nodal load must total P = {P_TOTAL}, got {applied_total}"
        );
        assert!(
            tip > 1e-3,
            "the measured tip deflection {tip:.9e} must be non-vacuous (> 1e-3)"
        );
    }

    // ========================================================================
    // t2025 Table 2 — the paper's own CURVED (annular-sector) plane-stress
    // cantilever (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// Ko, Bathe & Zhang (2025), C&S 308:107622, Table 2, "Vertical
    /// displacements at the tip (point A) for the curved beam problem"
    /// (Fig. 9): plane stress, `E = 1.0e3`, `nu = 0.0`, unit out-of-plane
    /// thickness, `1 x N` meshes (N straight-sided quads along the arc, one
    /// radially). (MEASUREMENT INSTRUMENT, NOT A GATE.)
    ///
    /// Why this exists (task T3): the Table 1 instrument reproduces the paper's
    /// regular-mesh MITC4+/D cell `0.0976755` to machine precision, but Table
    /// 1's REGULAR mesh has `x_d = 0`, so the 2017 assumed-membrane coefficients
    /// `a_A..a_E` are INACTIVE there. The one open fidelity question is whether
    /// our assumed membrane equals the 2025 one where those coefficients are
    /// live, i.e. on a curvature-capturing mesh with `x_d != 0`. This curved
    /// beam is exactly such a case: every flat trapezoidal element has
    /// `x_d != 0`, so `a_A..a_E` load the solution.
    ///
    /// Geometry (Fig. 9; cross-checked against beam theory): a 90-degree
    /// annular sector, inner radius `Ri = 10`, outer radius `Ro = 15` (radial
    /// wall thickness 5), unit thickness, `E = 1.0e3`, `nu = 0.0`. The sector
    /// is placed so the CLAMPED edge is the radial segment at `theta = 0`
    /// (from `(10, 0)` to `(15, 0)`) and the FREE edge is at `theta = 90`
    /// (from `(0, 10)` to `(0, 15)`). Point A is the free-end INNER corner,
    /// `(0, 10)` = node `inner(N)`. The tip load `P = 600` acts at A in `+y`,
    /// which at `theta = 90` is the OUTWARD RADIAL direction — a transverse
    /// (bending) load. Beam theory `delta = pi P R^3 / (4 E I)` with
    /// `R = (Ri + Ro)/2 = 12.5`, `I = 1 * 5^3 / 12` gives `P ~ 611` for
    /// `delta = 90`, matching the paper's `P = 600` and reference `90.1`.
    ///
    /// Published cells (Table 2, Fig. 9):
    /// ```text
    /// element     1x2        1x4        1x8
    /// MITC4       22.5988    57.9325    79.9218
    /// MITC4-IC    52.2291    84.6070    89.3023
    /// MITC4/D     51.2489    84.1086    89.1698
    /// MITC4+/D    51.2692    84.1444    89.2077   <- ours must match, or not
    /// reference   90.1 (all meshes)
    /// ```
    ///
    /// LOAD PLACEMENT / MEASUREMENT NODE (the reconstruction choices). A `1 x N`
    /// mesh has NO mid-surface node, so `P = 600` has to be reconstructed from
    /// the two free-edge nodes. The PRIMARY case places ALL 600 at point A
    /// (the inner tip node) and measures `u_y` at point A. Because beam theory
    /// needs the resultant to act at the mid-surface radius `R = 12.5`, this
    /// instrument ALSO sweeps the three natural placements (inner / outer /
    /// split 300-300) crossed with the three measurement nodes (inner / outer /
    /// average) and both element orientations, and prints all of them. Only a
    /// placement that puts the resultant at `R = 12.5` (i.e. the split) can
    /// reproduce the reference `90.1`; a single-corner placement at `R = 10` or
    /// `R = 15` shifts the effective lever arm by `R^3`.
    ///
    /// CURVED vs STRAIGHT-SIDED: this kernel is a 4-node element, so the arc is
    /// approximated by straight-sided chords (there is no mid-side-node element
    /// available). The mesh-refinement cases 1x2 / 1x4 / 1x8 in the paper are
    /// the counterpart of that approximation; the sweep prints all three for
    /// the reconstruction that best matches the published `1x4` cell.
    ///
    /// The assembly reuses the repository's own global element matrix
    /// (`compute_ke_global`, i.e. `T^T K_e T`) and the module's dense
    /// constrained solver (`solve_constrained`); only the curved mesh build and
    /// the scatter are local to this instrument. It MEASURES; it confirms
    /// nothing and changes no production line.
    ///
    /// Run with:
    /// `cargo test -p aeroelast-core t2025_table2 -- --ignored --nocapture`.
    #[test]
    #[ignore = "2025 paper Table 2 curved-beam measurement instrument; run with --ignored --nocapture"]
    fn t2025_table2_curved_plane_stress_beam() {
        // ---- published Table 2 cells (Ko, Bathe & Zhang (2025), C&S 308:107622) ----
        const PUB_MITC4_1X2: f64 = 22.5988;
        const PUB_MITC4_1X4: f64 = 57.9325;
        const PUB_MITC4_1X8: f64 = 79.9218;
        const PUB_MITC4_IC_1X4: f64 = 84.6070;
        const PUB_MITC4_D_1X4: f64 = 84.1086;
        const PUB_MITC4PLUS_D_1X2: f64 = 51.2692;
        const PUB_MITC4PLUS_D_1X4: f64 = 84.1444; // the column ours must match (or not)
        const PUB_MITC4PLUS_D_1X8: f64 = 89.2077;
        const PUB_REFERENCE: f64 = 90.1;
        // Deliberately broad physical bracket for the instrument: from the
        // paper's MITC4 1x2 cell up to just above its reference solution.
        const BRACKET_LO: f64 = PUB_MITC4_1X2;
        const BRACKET_HI: f64 = 92.0;
        // "Tight margin" for deciding whether the 1x4 reconstruction matched.
        const MATCH_MARGIN: f64 = 0.01; // 1%

        // ---- the instrument's dense solve over one curved `1 x N` mesh ----
        //
        // Returns `(u_y at the inner tip node, u_y at the outer tip node)` for
        // the given along-arc element count, nodal tip loads (in `+y`) and
        // element orientation. It reuses `compute_ke_global` + `solve_constrained`.
        fn solve_curved_beam(
            n_along: usize,
            inner_load: f64,
            outer_load: f64,
            natural_order: bool,
        ) -> (f64, f64) {
            // Fig. 9 geometry, material and load.
            const RI: f64 = 10.0;
            const RO: f64 = 15.0;
            const THICKNESS: f64 = 1.0;
            const E: f64 = 1.0e3;
            const NU: f64 = 0.0;
            const SHEAR_CORRECTION: f64 = 5.0 / 6.0;

            let n_per_arc = n_along + 1;
            let nnodes = 2 * n_per_arc;
            let ndof = 6 * nnodes;
            let nelem = n_along;

            let gdof = |n: usize, k: usize| 6 * n + k;
            let inner = |i: usize| i;
            let outer = |i: usize| n_per_arc + i;
            let theta =
                |i: usize| (i as f64) * std::f64::consts::FRAC_PI_2 / (n_along as f64);

            // Nodes on the two arcs at theta = 0 .. 90 degrees; clamped edge at
            // theta = 0, free edge at theta = 90, tip load pushing in +y.
            let mut node_xyz: Vec<[f64; 3]> = Vec::with_capacity(nnodes);
            for i in 0..n_per_arc {
                let t = theta(i);
                node_xyz.push([RI * t.cos(), RI * t.sin(), 0.0]);
            }
            for i in 0..n_per_arc {
                let t = theta(i);
                node_xyz.push([RO * t.cos(), RO * t.sin(), 0.0]);
            }

            // Straight-sided quads. `natural_order` is the inner-arc then
            // outer-arc listing [inner_e, inner_e+1, outer_e+1, outer_e]; the
            // reversed cycle flips the element normal (orientation variant).
            let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(nelem);
            let mut elem_coords: Vec<[[f64; 3]; 4]> = Vec::with_capacity(nelem);
            for e in 0..nelem {
                let n = if natural_order {
                    [inner(e), inner(e + 1), outer(e + 1), outer(e)]
                } else {
                    [inner(e), outer(e), outer(e + 1), inner(e + 1)]
                };
                elem_nodes.push(n);
                elem_coords.push([
                    node_xyz[n[0]],
                    node_xyz[n[1]],
                    node_xyz[n[2]],
                    node_xyz[n[3]],
                ]);
            }

            // Global assembly `K = sum_e T^T K_e T` (production global matrix).
            let mut kg = DMatrix::zeros(ndof, ndof);
            for e in 0..nelem {
                let constitutive = IsotropicMaterial::new(E, NU, 7800.0)
                    .constitutive(THICKNESS, SHEAR_CORRECTION);
                let pre = Mitc4PlusDPrecomputed::new(
                    &coords12(&elem_coords[e]),
                    constitutive,
                    THICKNESS,
                    SHEAR_CORRECTION,
                );
                let ke = compute_ke_global(&pre);
                for a in 0..4 {
                    for b in 0..4 {
                        for ka in 0..6 {
                            for kb in 0..6 {
                                kg[(gdof(elem_nodes[e][a], ka), gdof(elem_nodes[e][b], kb))] +=
                                    ke[(6 * a + ka, 6 * b + kb)];
                            }
                        }
                    }
                }
            }

            // Clamp all 6 DOF at the theta = 0 edge (inner + outer node 0).
            let mut prescribed: Vec<(usize, f64)> = Vec::new();
            for k in 0..6 {
                prescribed.push((gdof(inner(0), k), 0.0));
                prescribed.push((gdof(outer(0), k), 0.0));
            }

            // Tip load in `+y` at the free edge.
            let mut f: DVector<f64> = DVector::zeros(ndof);
            f[gdof(inner(n_along), 1)] += inner_load;
            f[gdof(outer(n_along), 1)] += outer_load;

            let u = solve_constrained(&kg, &f, &prescribed);
            (u[gdof(inner(n_along), 1)], u[gdof(outer(n_along), 1)])
        }

        // ---- PRIMARY: all P at point A (inner tip node), measured at A ----
        let (primary_inner, _primary_outer) = solve_curved_beam(4, 600.0, 0.0, true);
        let primary = primary_inner;

        println!("\n=== 2025 Table 2 instrument: CURVED plane-stress cantilever (Fig. 9) ===");
        println!(
            "geometry: 90-deg annular sector, Ri = 10, Ro = 15, t = 1, E = 1.0e3, nu = 0.0, shear correction {:.6}",
            5.0 / 6.0
        );
        println!(
            "primary case: 1x4 mesh (4 along arc x 1 radial = 10 nodes / 60 dof); P = 600 in +y at point A = inner tip node; measure u_y at point A"
        );
        println!("  our 1x4 value at A = {primary:.12e}");
        println!("  ratio to published MITC4+/D (84.1444) = {:.9}", primary / PUB_MITC4PLUS_D_1X4);
        println!("  ratio to published MITC4/D  (84.1086) = {:.9}", primary / PUB_MITC4_D_1X4);
        println!("  ratio to published reference(90.1  ) = {:.9}", primary / PUB_REFERENCE);
        println!(
            "  difference vs MITC4+/D = {:+.6e}; vs reference = {:+.6e}",
            primary - PUB_MITC4PLUS_D_1X4,
            primary - PUB_REFERENCE
        );

        // ---- RECONSTRUCTION SWEEP (1x4): load placement x measurement node x orientation ----
        // A 1xN mesh has no mid-surface node: the load resultant must sit at
        // R = (Ri+Ro)/2 = 12.5 to reproduce the reference 90.1, which is the
        // 300/300 split. This sweep is the evidence the parent needs to decide
        // whether any mismatch is a formulation difference or a reconstruction
        // artefact.
        let loads: [(&str, f64, f64); 3] = [
            ("all 600 at inner(A)", 600.0, 0.0),
            ("all 600 at outer", 0.0, 600.0),
            ("split 300/300", 300.0, 300.0),
        ];
        let measures = ["inner(A)", "outer", "average"];

        println!("\nreconstruction sweep on the 1x4 mesh (value at each measurement node):");
        println!(
            "  {:<20} {:<10} {:>16} {:>16} {:>16} {:>12}",
            "load placement", "orientation", "u_y inner(A)", "u_y outer", "u_y average", "dA/84.1444"
        );

        // Best reconstruction of the published 1x4 cell, tracked by index for
        // the refinement cases below: (relative error, load index, measure
        // index, natural-order flag).
        let mut best: (f64, usize, usize, bool) = (f64::INFINITY, 0, 0, true);
        for (load_idx, &(label, li, lo)) in loads.iter().enumerate() {
            for &natural in &[true, false] {
                let (ui, uo) = solve_curved_beam(4, li, lo, natural);
                let vals = [ui, uo, 0.5 * (ui + uo)];
                let orient = if natural { "natural" } else { "reversed" };
                println!(
                    "  {:<20} {:<10} {:>16.9} {:>16.9} {:>16.9} {:>12.6}",
                    label,
                    orient,
                    vals[0],
                    vals[1],
                    vals[2],
                    vals[0] / PUB_MITC4PLUS_D_1X4
                );
                for (m, &v) in vals.iter().enumerate() {
                    let err = (v - PUB_MITC4PLUS_D_1X4).abs() / PUB_MITC4PLUS_D_1X4;
                    if err < best.0 {
                        best = (err, load_idx, m, natural);
                    }
                }
            }
        }

        let (best_err, best_load, best_measure, best_natural) = best;
        let best_label = loads[best_load].0;
        println!(
            "\nclosest 1x4 reconstruction to 84.1444: load = {best_label}, measure = {}, orientation = {} -> relative error {:.4}%",
            measures[best_measure],
            if best_natural { "natural" } else { "reversed" },
            100.0 * best_err
        );
        let matched = best_err <= MATCH_MARGIN;
        println!(
            "1x4 match within the {:.1}% tight margin: {}",
            100.0 * MATCH_MARGIN,
            matched
        );

        // ---- REFINEMENT CASES: add 1x2 / 1x8 only when the 1x4 matched ----
        if matched {
            // Recover the best reconstruction's load pair.
            let (li, lo) = (loads[best_load].1, loads[best_load].2);
            println!(
                "\n1x2 / 1x8 for the matched reconstruction (load = {best_label}, measure = {}, orientation = {}):",
                measures[best_measure],
                if best_natural { "natural" } else { "reversed" }
            );
            for (n, published) in
                [(2usize, PUB_MITC4PLUS_D_1X2), (8usize, PUB_MITC4PLUS_D_1X8)]
            {
                let (ui, uo) = solve_curved_beam(n, li, lo, best_natural);
                let v = match best_measure {
                    0 => ui,
                    1 => uo,
                    _ => 0.5 * (ui + uo),
                };
                println!(
                    "  1x{n}: our value = {v:.9}  published MITC4+/D = {published:.4}  ratio = {:.6}  difference = {:+.6e}",
                    v / published,
                    v - published
                );
            }
            println!("  reference (all meshes) = {PUB_REFERENCE:.4}");
        } else {
            println!(
                "\nno 1x4 reconstruction is within {:.1}% of 84.1444; 1x2 and 1x8 are NOT added (per the task instruction).",
                100.0 * MATCH_MARGIN
            );
        }

        // ---- loose instrument assertions on the PRIMARY value (the parent draws the verdict) ----
        assert!(
            primary.is_finite(),
            "the measured tip displacement must be finite, got {primary}"
        );
        assert!(
            primary >= BRACKET_LO && primary <= BRACKET_HI,
            "the measured tip displacement {primary:.9e} must lie inside the broad physical bracket [{BRACKET_LO}, {BRACKET_HI}]"
        );
        assert!(
            primary > 1e-3,
            "the measured tip displacement {primary:.9e} must be non-vacuous (> 1e-3)"
        );
    }

    // ========================================================================
    // N-alpha — faithful incremental kinematics: the rigid-rotation identity
    // (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// N-alpha — the rigid-body-rotation identity of Ko, Lee & Bathe (2017),
    /// C&S 185:1-14, Eq. (9) (MEASUREMENT INSTRUMENT, NOT A GATE).
    ///
    /// A finite rigid-body rotation about an in-plane axis (here `e1`, so
    /// `theta = phi e1` is perpendicular to `V_n = e3`) must produce a zero
    /// incremental Green-Lagrange strain `_0 e_ij + _0 eta_ij` to the paper's
    /// retained order. The current state is the flat `RECT` element (time `t`);
    /// the increment is `u_i = (R - I) x_i`, `alpha_i = theta . V_1^i`,
    /// `beta_i = theta . V_2^i`.
    ///
    /// The measurement prints, for angles 10/30/60/90 degrees:
    ///   * the max `|_0 e + _0 eta|` over several `(r,s,zeta)` points, `zeta` in
    ///     {-1, 0, 1}, and its `phi^3`-normalised value (the retained order);
    ///   * the same residual with the `u_2` term forced to zero (the
    ///     non-vacuity control: without the quadratic `u_b2` of Eq. (5c) the
    ///     residual is `O(phi^2)` and much larger);
    ///   * the production nonlinear path's residual on the same rigid
    ///     rotation, `||compute_fint_global(.., true)|| / (||K_0|| ||u||)`,
    ///     which after N-delta (option B) IS the faithful pair and must be at
    ///     round-off;
    ///   * the SUPERSEDED bounded path's residual,
    ///     `||compute_fint_global_bounded_superseded(..)|| / (||K_0|| ||u||)`,
    ///     kept only as the recorded historical "before".
    ///
    /// Run it with `cargo test -p aeroelast-core n_alpha -- --ignored --nocapture`.
    #[test]
    #[ignore = "N-alpha rigid-rotation identity measurement; run with --ignored --nocapture"]
    fn n_alpha_rigid_rotation_gives_zero_gl_strain_increment() {
        let pre = pre_from(&RECT);
        let state = GlCurrentState {
            coords: pre.initial_coords_3d,
            vn: pre.vn,
            v1: pre.v1,
            v2: pre.v2,
            a_i: pre.a_i,
        };

        // Several (r, s, zeta) points, including the surfaces zeta = +-1.
        let points: [(f64, f64, f64); 9] = [
            (0.0, 0.0, -1.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.577_350_269_189_625_8, -0.577_350_269_189_625_8, -1.0),
            (0.577_350_269_189_625_8, 0.577_350_269_189_625_8, 1.0),
            (-0.4, 0.7, 1.0),
            (0.8, 0.3, 0.5),
            (1.0, -1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ];

        // Frobenius scale of the linear stiffness, for the current-path residual.
        let k0 = compute_ke_local(&pre);
        let mut k_fro = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                k_fro += k0[(i, j)] * k0[(i, j)];
            }
        }
        let k_fro = k_fro.sqrt();

        const AXIS: [f64; 3] = [1.0, 0.0, 0.0];
        let angles_deg = [10.0f64, 30.0, 60.0, 90.0];

        println!("\nN-alpha rigid rotation about e1 (flat RECT, current state = flat):");
        println!(
            "{:>8}  {:>14}  {:>16}  {:>14}  {:>16}  {:>16}  {:>14}  {:>14}",
            "phi[deg]",
            "max|e+eta|",
            "max|e+eta|/phi^3",
            "max|e|(u2=0)",
            "max|e|/phi^2",
            "ratio(no/yes)",
            "prod fint/Ku",
            "bounded fint/Ku"
        );

        let mut worst_phi3 = 0.0f64;
        for &deg in &angles_deg {
            let phi = deg.to_radians();
            let r_mat = rotation_matrix(AXIS, phi);
            let theta = phi * Vector3::new(AXIS[0], AXIS[1], AXIS[2]);

            // Increment: u_i = (R - I) x_i; alpha/beta from theta and the
            // current directors (Eq. 4c).
            let mut inc = GlIncrement {
                u: [[0.0f64; 3]; 4],
                alpha: [0.0f64; 4],
                beta: [0.0f64; 4],
            };
            let mut director_diff = 0.0f64;
            for i in 0..4 {
                let x = Vector3::new(
                    state.coords[i][0],
                    state.coords[i][1],
                    state.coords[i][2],
                );
                let urel = r_mat * x - x;
                inc.u[i] = [urel[0], urel[1], urel[2]];
                inc.alpha[i] = theta.dot(&state.v1[i]);
                inc.beta[i] = theta.dot(&state.v2[i]);

                // Eq. (4c) consistency: the quadratic director increment must
                // reproduce the exact rotated director to O(phi^3).
                let d_inc = director_increment_quadratic(&state.vn[i], &theta);
                let d_exact = r_mat * state.vn[i] - state.vn[i];
                director_diff = director_diff.max((d_inc - d_exact).norm());
            }

            let mut with_u2 = 0.0f64;
            let mut without_u2 = 0.0f64;
            for &(r, s, zeta) in &points {
                let (e_lin, eta_nl) = gl_strain_increment(&state, &inc, r, s, zeta);
                for k in 0..6 {
                    with_u2 = with_u2.max((e_lin[k] + eta_nl[k]).abs());
                }
                // Non-vacuity control: force the Eq. (6) quadratic part u_2 to
                // zero, keeping only the paper's three Eq. (9) terms.
                let g = gl_current_base_vectors(&state, r, s, zeta);
                let (u1, u2) = incremental_disp_gradients(&state, &inc, r, s, zeta);
                let u2_zero = [Vector3::zeros(); 3];
                let (e_lin0, eta_nl0) = gl_strain_increment_components(&g, &u1, &u2_zero);
                for k in 0..6 {
                    without_u2 = without_u2.max((e_lin0[k] + eta_nl0[k]).abs());
                }
                // The point split of Eq. (6) and the gradient set must agree:
                // d/dzeta of the split is `u_b1` / `u_b2` (the third gradient
                // entries), which is exactly what the strain assembly uses.
                let (u1_at_0, u2_at_0) = incremental_disp_split(&state, &inc, r, s, 0.0);
                let (u1_at_1, u2_at_1) = incremental_disp_split(&state, &inc, r, s, 1.0);
                let dn1 = ((u1_at_1 - u1_at_0) - u1[2]).norm();
                let dn2 = ((u2_at_1 - u2_at_0) - u2[2]).norm();
                assert!(
                    dn1 < 1.0e-12 && dn2 < 1.0e-12,
                    "Eq. (6) split gradient mismatch at ({r},{s},{zeta}): {dn1:.3e}, {dn2:.3e}"
                );
            }

            // Production nonlinear path (N-delta, option B): the faithful
            // pair. It must satisfy the rigid-rotation identity to round-off.
            let mut u_rot = Vec24::zeros();
            for i in 0..4 {
                u_rot[6 * i] = inc.u[i][0];
                u_rot[6 * i + 1] = inc.u[i][1];
                u_rot[6 * i + 2] = inc.u[i][2];
                u_rot[6 * i + 3] = theta[0];
                u_rot[6 * i + 4] = theta[1];
                u_rot[6 * i + 5] = theta[2];
            }
            let u_norm = u_rot.iter().map(|v| v * v).sum::<f64>().sqrt();
            let f_cur = compute_fint_global(&pre, &u_rot, true);
            let f_norm = f_cur.iter().map(|v| v * v).sum::<f64>().sqrt();
            let f_rel = f_norm / (k_fro * u_norm);
            // Superseded bounded path: the historical "before" number, printed
            // for the record only.
            let f_bounded = compute_fint_global_bounded_superseded(&pre, &u_rot);
            let fb_norm = f_bounded.iter().map(|v| v * v).sum::<f64>().sqrt();
            let f_bounded_rel = fb_norm / (k_fro * u_norm);

            let phi3 = phi.powi(3);
            let phi2 = phi * phi;
            worst_phi3 = worst_phi3.max(with_u2 / phi3);
            println!(
                "{:>8.1}  {:>14.6e}  {:>16.6e}  {:>14.6e}  {:>16.6e}  {:>16.3e}  {:>14.6e}  {:>14.6e}  (Eq.4c dn={:.2e})",
                deg,
                with_u2,
                with_u2 / phi3,
                without_u2,
                without_u2 / phi2,
                without_u2 / with_u2,
                f_rel,
                f_bounded_rel,
                director_diff
            );

            // The paper's identity, to the retained order: the residual is
            // bounded by O(phi^3).
            assert!(
                with_u2 <= 1.0e-12 + 1.0 * phi3,
                "phi = {deg} deg: max |_0 e + _0 eta| = {with_u2:.3e} exceeds the retained order phi^3 = {phi3:.3e}"
            );
            // Non-vacuity: with u_2 zeroed the residual is O(phi^2), so the
            // quadratic director term of Eq. (5c) is what removes the leading
            // error and this measurement discriminates.
            assert!(
                without_u2 > 1.5 * with_u2,
                "phi = {deg} deg: dropping u_2 must break the identity; with = {with_u2:.3e}, without = {without_u2:.3e}"
            );
            assert!(
                without_u2 > 1.0e-3,
                "phi = {deg} deg: the u_2 control must be non-vacuous, got {without_u2:.3e}"
            );
            // Eq. (4c) is reproduced to quadratic order.
            assert!(
                director_diff <= 1.0 * phi3,
                "phi = {deg} deg: Eq. (4c) director error {director_diff:.3e} exceeds phi^3 = {phi3:.3e}"
            );
            // The production nonlinear path (N-delta, option B) is now the
            // faithful pair, so it satisfies the rigid-rotation identity to
            // round-off; this is the OPPOSITE of the superseded bounded path's
            // recorded failure below.
            assert!(
                f_rel < 1.0e-10,
                "phi = {deg} deg: the production nonlinear path must satisfy the rigid rotation to round-off, got {f_rel:.3e}"
            );
            // The superseded bounded path does not satisfy the identity: the
            // historical "before" that motivated the faithful implementation.
            assert!(
                f_bounded_rel > 1.0e-5,
                "phi = {deg} deg: the superseded bounded path was expected to fail the rigid rotation, got {f_bounded_rel:.3e}"
            );
        }
        println!(
            "worst max|e+eta|/phi^3 over the sweep: {worst_phi3:.6e} (roughly angle independent = O(phi^3))"
        );
    }

    // ========================================================================
    // N-beta — the paper's assumed fields on the TL strains
    // (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// N-beta — the assumed strain fields of Ko, Lee & Bathe (2017),
    /// C&S 185:1-14, Eqs. (10) and (15a-d) (MEASUREMENT INSTRUMENT, NOT A GATE).
    ///
    /// (1) **DECISIVE acceptance.** The N-beta assumed INCREMENTAL membrane
    ///     operator at `u = 0` must equal the LINEAR assumed membrane operator
    ///     (the same class of check that exposed the N-gamma factor-2 defect):
    ///     the covariant rows are compared to `b_membrane_covariant_2017` and,
    ///     with the engineering-shear/`^0`-frame mapping applied, to
    ///     `b_membrane_2017`, at the sample points AND the five tying points.
    ///     `max|difference| / max|linear| <= 1e-10`. The N-beta operator is
    ///     built by central differences of [`gl_assumed_membrane_increment`],
    ///     so it cannot pass by sharing a factor with the reference.
    /// (2) **Flat reduction.** On a flat parallelogram (`^t x_d = 0`) the current
    ///     coefficients `^t a_A..^t a_E` vanish, Eqs. (15b-d) collapse onto the
    ///     linear interpolation, and the assumed membrane must reproduce the
    ///     displacement-based TL membrane. Both sides now use the SAME,
    ///     explicitly-stated convention `1/2(m_cur - m_init)` with
    ///     `m = x_m . x_m` -- NOT `gl_mid_metrics`, whose `2 m` doubling
    ///     cancelled the old defect. Two cases: the literal `^t x = ^0 x` and a
    ///     stretched/sheared flat parallelogram with a nonzero membrane.
    /// (3) **Non-vacuity.** Perturbing one coefficient away from its natural
    ///     zero (`^t a_E = 0.5`) breaks the reduction; setting `^t a_E = 0` on a
    ///     distorted flat quad where it is naturally nonzero changes the field.
    /// (4) **Rigid rotation.** With the assumed fields (Eqs. 15b-d and Eq. 10)
    ///     applied to the N-alpha Eq. (9) increment of a finite rigid rotation
    ///     about the in-plane axis `e1`, the residual must stay at `O(phi^3)`.
    ///
    /// A `1e-12`-scale magnitude report accompanies (1): a common-factor
    /// regression in the assumed field makes `max|difference| / max|linear|`
    /// jump from round-off to `O(1)`, which the report makes visible.
    ///
    /// Run it with `cargo test -p aeroelast-core n_beta -- --ignored --nocapture`.
    #[test]
    #[ignore = "N-beta assumed-field measurement; run with --ignored --nocapture"]
    fn n_beta_assumed_fields_reduce_to_displacement_based_on_flat() {
        let pre = pre_from(&RECT);
        let state0 = GlCurrentState {
            coords: pre.initial_coords_3d,
            vn: pre.vn,
            v1: pre.v1,
            v2: pre.v2,
            a_i: pre.a_i,
        };

        let init_vecs = [pre.x_r, pre.x_s, pre.x_d];
        let init_tie = gl_tying_metrics(&init_vecs[0], &init_vecs[1], &init_vecs[2]);

        let points: [(f64, f64); 5] = [
            (0.0, 0.0),
            (0.37, -0.62),
            (0.8, 0.3),
            (-1.0, 1.0),
            (0.577_350_269_189_625_8, -0.577_350_269_189_625_8),
        ];

        // Displacement-based TL membrane with the SAME convention as Eq. (15a)
        // after undoing the Eq. (12b) doubling: `1/2(m_cur - m_init)` with
        // `m = x_m . x_m`, built directly from the geometry and NOT from
        // `gl_mid_metrics` (which returns the doubled `2 m`).
        let direct_membrane = |state: &GlCurrentState, r: f64, s: f64| -> [f64; 3] {
            let metric_of = |v: &[Vector3<f64>; 3]| -> [f64; 3] {
                let g_r = v[0] + s * v[2];
                let g_s = v[1] + r * v[2];
                [g_r.dot(&g_r), g_s.dot(&g_s), g_r.dot(&g_s)]
            };
            let (cur, _n, _mr, _ms) = gl_current_characteristic_vectors(state);
            let mc = metric_of(&cur);
            let mi = metric_of(&init_vecs);
            [
                0.5 * (mc[0] - mi[0]),
                0.5 * (mc[1] - mi[1]),
                0.5 * (mc[2] - mi[2]),
            ]
        };

        // ---- (1) DECISIVE: the N-beta assumed INCREMENTAL membrane operator at
        //      u = 0 equals the LINEAR assumed membrane operator ----
        //
        // Central differences of `gl_assumed_membrane_increment` in the
        // increment's nodal DOFs. The membrane metric depends only on the nodal
        // positions (not the directors), so the rotational columns are
        // identically zero; the flat RECT's local frame is the global one.
        let nb_operator = |state: &GlCurrentState, r: f64, s: f64| -> [[f64; 24]; 3] {
            let h = 1.0e-4;
            let mut op = [[0.0f64; 24]; 3];
            for j in 0..12 {
                let node = j / 3;
                let k = j % 3;
                let col = 6 * node + k;
                let mut inc_p = GlIncrement {
                    u: [[0.0; 3]; 4],
                    alpha: [0.0; 4],
                    beta: [0.0; 4],
                };
                let mut inc_m = inc_p;
                inc_p.u[node][k] = h;
                inc_m.u[node][k] = -h;
                let ep = gl_assumed_membrane_increment(state, &inc_p, r, s);
                let em = gl_assumed_membrane_increment(state, &inc_m, r, s);
                for row in 0..3 {
                    op[row][col] = (ep[row] - em[row]) / (2.0 * h);
                }
            }
            op
        };

        let op_points: [(f64, f64); 10] = [
            (0.0, 0.0),
            (0.37, -0.62),
            (0.8, 0.3),
            (-1.0, 1.0),
            (0.577_350_269_189_625_8, -0.577_350_269_189_625_8),
            (0.0, 1.0),  // A
            (0.0, -1.0), // B
            (1.0, 0.0),  // C
            (-1.0, 0.0), // D
            (0.5, -0.25),
        ];
        let mut op_diff_cov = 0.0f64;
        let mut op_scale_cov = 0.0f64;
        let mut op_diff_loc = 0.0f64;
        let mut op_scale_loc = 0.0f64;
        for &(r, s) in &op_points {
            let nb = nb_operator(&state0, r, s);
            let bcov = b_membrane_covariant_2017(&pre, r, s);
            let bloc = b_membrane_2017(&pre, r, s);
            let tmap = covariant_to_local_mapping(&j_loc_at(&pre, r, s));
            for j in 0..24 {
                // Covariant rows (tensor shear).
                let nb_cov = Vector3::new(nb[0][j], nb[1][j], nb[2][j]);
                let lin_cov = Vector3::new(bcov[(0, j)], bcov[(1, j)], bcov[(2, j)]);
                op_diff_cov = op_diff_cov.max((nb_cov - lin_cov).norm());
                op_scale_cov = op_scale_cov.max(bcov[(0, j)].abs());
                op_scale_cov = op_scale_cov.max(bcov[(1, j)].abs());
                op_scale_cov = op_scale_cov.max(bcov[(2, j)].abs());
                // Local rows: engineering shear (`2 e_12`) and the Eq. (23)
                // `^0`-metric mapping, matching `b_membrane_2017`.
                let nb_loc = tmap * Vector3::new(nb[0][j], nb[1][j], 2.0 * nb[2][j]);
                let lin_loc = Vector3::new(bloc[(0, j)], bloc[(1, j)], bloc[(2, j)]);
                op_diff_loc = op_diff_loc.max((nb_loc - lin_loc).norm());
                op_scale_loc = op_scale_loc.max(bloc[(0, j)].abs());
                op_scale_loc = op_scale_loc.max(bloc[(1, j)].abs());
                op_scale_loc = op_scale_loc.max(bloc[(2, j)].abs());
            }
        }
        let op_rel_cov = op_diff_cov / op_scale_cov;
        let op_rel_loc = op_diff_loc / op_scale_loc;
        println!("\nN-beta (1) DECISIVE: incremental assumed membrane operator at u = 0");
        println!(
            "  covariant rows vs b_membrane_covariant_2017: max|diff| = {op_diff_cov:.6e}, max|linear| = {op_scale_cov:.6e}, rel = {op_rel_cov:.6e}"
        );
        println!(
            "  local rows (2 e_12) vs b_membrane_2017:      max|diff| = {op_diff_loc:.6e}, max|linear| = {op_scale_loc:.6e}, rel = {op_rel_loc:.6e}"
        );
        println!(
            "  a common-factor regression would make either rel O(1) instead of ~1e-12"
        );
        assert!(
            op_rel_cov <= 1.0e-10,
            "N-beta incremental covariant operator must equal b_membrane_covariant_2017, rel = {op_rel_cov:.3e}"
        );
        assert!(
            op_rel_loc <= 1.0e-10,
            "N-beta incremental local operator must equal b_membrane_2017, rel = {op_rel_loc:.3e}"
        );

        // ---- (1a) literal: ^t x = ^0 x ----
        let (vecs0, _n0, mr0, ms0) = gl_current_characteristic_vectors(&state0);
        assert!(
            vecs0[2].norm() < 1.0e-14,
            "^t x_d = {} must vanish on the flat RECT",
            vecs0[2]
        );
        let (_cr0, _cs0, _d0, coeff0) = compute_membrane_coefficients_2017(&vecs0[2], &mr0, &ms0);
        for (k, a) in coeff0.iter().enumerate() {
            assert!(
                a.abs() < 1.0e-14,
                "^t a[{k}] = {a} must vanish on the flat RECT"
            );
        }
        let mut res_a = 0.0f64;
        for &(r, s) in &points {
            let assumed = gl_assumed_membrane_strain(&pre, &state0, r, s);
            let direct = direct_membrane(&state0, r, s);
            for k in 0..3 {
                res_a = res_a.max((assumed[k] - direct[k]).abs());
            }
        }
        println!("\nN-beta (1a) flat reduction, literal ^t x = ^0 x: max|e~^m - e^m| = {res_a:.6e}");
        assert!(res_a < 1.0e-14, "residual {res_a:.3e}");

        // ---- (1b) flat parallelogram, genuinely deformed (nonzero membrane) ----
        // An affine image of the flat RECT: still flat, still a parallelogram,
        // so ^t x_d = 0 and every current coefficient vanishes.
        let deformed: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.4, 1.1, 0.0],
            [0.4, 1.1, 0.0],
        ];
        let state = GlCurrentState {
            coords: deformed,
            vn: pre.vn,
            v1: pre.v1,
            v2: pre.v2,
            a_i: pre.a_i,
        };
        let (vecs, _n, mr, ms) = gl_current_characteristic_vectors(&state);
        let scale = vecs[0].norm().max(vecs[1].norm());
        assert!(
            vecs[2].norm() <= 1.0e-14 * scale,
            "^t x_d = {} must vanish on a parallelogram",
            vecs[2]
        );
        let (_cr, _cs, _d, coeff) = compute_membrane_coefficients_2017(&vecs[2], &mr, &ms);
        for (k, a) in coeff.iter().enumerate() {
            assert!(
                a.abs() <= 1.0e-12,
                "^t a[{k}] = {a} must vanish on a parallelogram"
            );
        }
        let mut res_b = 0.0f64;
        let mut mem_norm = 0.0f64;
        for &(r, s) in &points {
            let assumed = gl_assumed_membrane_strain(&pre, &state, r, s);
            let direct = direct_membrane(&state, r, s);
            for k in 0..3 {
                res_b = res_b.max((assumed[k] - direct[k]).abs());
                mem_norm = mem_norm.max(direct[k].abs());
            }
        }
        println!(
            "N-beta (1b) flat reduction, deformed parallelogram: max|e~^m - e^m| = {res_b:.6e}, max|e^m| = {mem_norm:.6e}"
        );
        assert!(
            mem_norm > 1.0e-3,
            "the deformed-parallelogram membrane must be non-vacuous, got {mem_norm:.3e}"
        );
        assert!(res_b < 1.0e-12, "residual {res_b:.3e}");

        // ---- (2a) non-vacuity: perturb one coefficient away from zero ----
        let cur_tie = gl_tying_metrics(&vecs[0], &vecs[1], &vecs[2]);
        let mut coef_pert = coeff;
        coef_pert[4] = 0.5;
        let mut res_pert = 0.0f64;
        for &(r, s) in &points {
            let assumed = gl_assumed_mid_metric(&coef_pert, &cur_tie, r, s);
            let init_assumed = gl_assumed_mid_metric(&pre.a_coeffs, &init_tie, r, s);
            let direct = direct_membrane(&state, r, s);
            for k in 0..3 {
                let e = 0.25 * assumed[k] - 0.25 * init_assumed[k];
                res_pert = res_pert.max((e - direct[k]).abs());
            }
        }
        println!(
            "N-beta (2a) non-vacuity, ^t a_E perturbed 0.0 -> 0.5 on the flat parallelogram: max|e~^m - e^m| = {res_pert:.6e}"
        );
        assert!(
            res_pert > 1.0e-3,
            "perturbing ^t a_E must break the flat reduction, got {res_pert:.3e}"
        );

        // ---- (2b) the literal control: set ^t a_E = 0 where it is nonzero ----
        let pre_d = pre_from(&FLAT_DISTORTED);
        let state_d = GlCurrentState {
            coords: FLAT_DISTORTED,
            vn: pre_d.vn,
            v1: pre_d.v1,
            v2: pre_d.v2,
            a_i: pre_d.a_i,
        };
        let (vecs_d, _n_d, mr_d, ms_d) = gl_current_characteristic_vectors(&state_d);
        let (_cr_d, _cs_d, _d_d, coeff_d) =
            compute_membrane_coefficients_2017(&vecs_d[2], &mr_d, &ms_d);
        assert!(
            coeff_d[4].abs() > 1.0e-6,
            "FLAT_DISTORTED must have ^t a_E != 0, got {}",
            coeff_d[4]
        );
        let tie_d = gl_tying_metrics(&vecs_d[0], &vecs_d[1], &vecs_d[2]);
        let mut coef_zero = coeff_d;
        coef_zero[4] = 0.0;
        let mut d_zero = 0.0f64;
        for &(r, s) in &points {
            let natural = gl_assumed_mid_metric(&coeff_d, &tie_d, r, s);
            let zeroed = gl_assumed_mid_metric(&coef_zero, &tie_d, r, s);
            for k in 0..3 {
                d_zero = d_zero.max((natural[k] - zeroed[k]).abs());
            }
        }
        println!(
            "N-beta (2b) non-vacuity, ^t a_E set to 0 on FLAT_DISTORTED (natural a_E = {:.6e}): max|e~^m(a_E) - e~^m(0)| = {d_zero:.6e}",
            coeff_d[4]
        );
        assert!(
            d_zero > 1.0e-6,
            "setting ^t a_E = 0 must change the assumed field, got {d_zero:.3e}"
        );

        // ---- (3) rigid rotation with the assumed fields applied ----
        //
        // Two in-plane axes: `e1` (the N-alpha instrument's axis) and the
        // in-plane diagonal `(e1 + e2)/sqrt(2)`. About `e1` the membrane
        // residual is identically zero by symmetry and only the assumed shear
        // is exercised; the diagonal axis makes the assumed membrane residual
        // nonzero, so the membrane field is measured too.
        let diag = (Vector3::new(1.0, 0.0, 0.0) + Vector3::new(0.0, 1.0, 0.0)).normalize();
        let axes: [(&str, [f64; 3]); 2] = [
            ("e1", [1.0, 0.0, 0.0]),
            ("(e1+e2)/sqrt2", [diag[0], diag[1], diag[2]]),
        ];
        let angles_deg = [10.0f64, 30.0, 60.0, 90.0];
        let rot_points: [(f64, f64, f64); 9] = [
            (0.0, 0.0, -1.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.577_350_269_189_625_8, -0.577_350_269_189_625_8, -1.0),
            (0.577_350_269_189_625_8, 0.577_350_269_189_625_8, 1.0),
            (-0.4, 0.7, 1.0),
            (0.8, 0.3, 0.5),
            (1.0, -1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ];

        for &(axis_name, axis) in &axes {
            println!(
                "\nN-beta (3) rigid rotation about the in-plane axis {axis_name} with the assumed fields applied (flat RECT, current state = flat):"
            );
            println!(
                "{:>8}  {:>16}  {:>16}  {:>16}  {:>16}  {:>14}  {:>14}",
                "phi[deg]",
                "max|e~^m|",
                "/phi^3",
                "max|e~_sh|",
                "/phi^3",
                "ratio m/d",
                "ratio s/d"
            );

            let mut worst_mem_phi3 = 0.0f64;
            let mut worst_sh_phi3 = 0.0f64;
            for &deg in &angles_deg {
                let phi = deg.to_radians();
                let r_mat = rotation_matrix(axis, phi);
                let theta = phi * Vector3::new(axis[0], axis[1], axis[2]);

                let mut inc = GlIncrement {
                    u: [[0.0f64; 3]; 4],
                    alpha: [0.0f64; 4],
                    beta: [0.0f64; 4],
                };
                for i in 0..4 {
                    let x = Vector3::new(
                        state0.coords[i][0],
                        state0.coords[i][1],
                        state0.coords[i][2],
                    );
                    let urel = r_mat * x - x;
                    inc.u[i] = [urel[0], urel[1], urel[2]];
                    inc.alpha[i] = theta.dot(&state0.v1[i]);
                    inc.beta[i] = theta.dot(&state0.v2[i]);
                }

                let mut mem_assumed = 0.0f64;
                let mut mem_direct = 0.0f64;
                let mut sh_assumed = 0.0f64;
                let mut sh_direct = 0.0f64;
                for &(r, s, zeta) in &rot_points {
                    let a = gl_assumed_membrane_increment(&state0, &inc, r, s);
                    let (e_lin, eta_nl) = gl_strain_increment(&state0, &inc, r, s, 0.0);
                    for k in 0..3 {
                        mem_assumed = mem_assumed.max(a[k].abs());
                        mem_direct = mem_direct.max((e_lin[k] + eta_nl[k]).abs());
                    }
                    let (sh_r, sh_s) = gl_assumed_transverse_shear(&state0, &inc, r, s, zeta);
                    let (e_lin6, eta_nl6) = gl_strain_increment(&state0, &inc, r, s, zeta);
                    for (idx, v) in [(4usize, sh_r), (5usize, sh_s)] {
                        sh_assumed = sh_assumed.max(v.abs());
                        sh_direct = sh_direct.max((e_lin6[idx] + eta_nl6[idx]).abs());
                    }
                }

                let phi3 = phi.powi(3);
                worst_mem_phi3 = worst_mem_phi3.max(mem_assumed / phi3);
                worst_sh_phi3 = worst_sh_phi3.max(sh_assumed / phi3);
                println!(
                    "{:>8.1}  {:>16.6e}  {:>16.6e}  {:>16.6e}  {:>16.6e}  {:>14.3}  {:>14.3}",
                    deg,
                    mem_assumed,
                    mem_assumed / phi3,
                    sh_assumed,
                    sh_assumed / phi3,
                    mem_assumed / mem_direct.max(1.0e-30),
                    sh_assumed / sh_direct.max(1.0e-30)
                );

                // The assumed fields are bounded linear combinations of the
                // Eq. (9) tying values; both residuals must stay at the retained
                // order `O(phi^3)` relative to their own displacement-based value.
                assert!(
                    mem_assumed <= 1.0e-12 + 20.0 * mem_direct,
                    "phi = {deg} deg: assumed membrane {mem_assumed:.3e} exceeds 20x the direct {mem_direct:.3e}"
                );
                assert!(
                    sh_assumed <= 1.0e-12 + 20.0 * sh_direct,
                    "phi = {deg} deg: assumed shear {sh_assumed:.3e} exceeds 20x the direct {sh_direct:.3e}"
                );
                assert!(
                    mem_direct <= 1.0 * phi3,
                    "phi = {deg} deg: direct membrane {mem_direct:.3e} exceeds phi^3 = {phi3:.3e}"
                );
                assert!(
                    sh_direct <= 1.0 * phi3,
                    "phi = {deg} deg: direct shear {sh_direct:.3e} exceeds phi^3 = {phi3:.3e}"
                );
            }
            println!(
                "worst assumed-field max/phi^3 on {axis_name}: membrane {worst_mem_phi3:.6e}, shear {worst_sh_phi3:.6e} (roughly angle independent = O(phi^3))"
            );
        }
    }

    // ========================================================================
    // N-gamma — the faithful TL internal force and tangent
    // (MEASUREMENT INSTRUMENT, NOT A GATE; `#[ignore]`d on purpose)
    // ========================================================================

    /// N-gamma — the implementation-consistency identities of the faithful
    /// total-Lagrangian pair `n_gamma_fint_local` (Eq. (24b)) and
    /// `n_gamma_kt_local` (Eq. (24a)):
    ///
    /// 1. `^t_0 F_e = 0` at `u = 0`, exactly (it is `B^T W e~` with `e~ = 0`).
    /// 2. `^t K_e = d ^t_0 F_e / d u` by central differences — Eq. (24a) must be
    ///    the consistent linearisation of Eq. (24b). This is measured for the 12
    ///    rigid motions below AND for general random increments at several `||u||`
    ///    (check 4), because the geometric term vanishes on a rigid motion.
    /// 3. `^t K_e(0) = K_0` to round-off (the exact `u = 0` limit of `dF/du`; the
    ///    earlier bit-for-bit equality was a construction artefact of a removed
    ///    additive reference — see `n_gamma_kt_local`).
    ///
    /// The flat `RECT` element is used because its local frame coincides with
    /// the global one, so `build_t24` is the identity and the local and global
    /// 24-vectors agree; the measurement is then about the formulation, not the
    /// frame.
    ///
    /// The rigid-body residual `||F|| / (||K_0|| ||u||)` is PRINTED for the
    /// faithful local pair, for the PRODUCTION global path
    /// (`compute_fint_global(.., true)`), which after N-delta (option B) IS the
    /// faithful pair, and for the SUPERSEDED bounded path
    /// (`compute_fint_global_bounded_superseded(..)`) kept as the recorded
    /// historical "before". The production value is ASSERTED at round-off; the
    /// bounded value is printed only. The exact strain-level discriminator lives
    /// in `n_alpha_rigid_rotation_gives_zero_gl_strain_increment`.
    #[test]
    #[ignore = "N-gamma rigid-body force + consistent-tangent measurement; run with --ignored --nocapture"]
    fn n_gamma_rigid_body_zero_force_and_consistent_tangent() {
        let pre = pre_from(&RECT);
        let state = GlCurrentState {
            coords: pre.initial_coords_3d,
            vn: pre.vn,
            v1: pre.v1,
            v2: pre.v2,
            a_i: pre.a_i,
        };
        let k0 = compute_ke_local(&pre);
        let k0_scale = max_abs(&k0).max(1e-30);
        let zero = Vec24::zeros();
        let vec_norm = |v: &[f64; 24]| -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() };

        // (1) exactly zero at u = 0.
        let f0 = n_gamma_fint_local(&pre, &state, &zero);
        let mut f0_buf = [0.0f64; 24];
        for i in 0..24 {
            f0_buf[i] = f0[i];
        }
        let f0_norm = vec_norm(&f0_buf);
        println!("\nN-gamma (1) |F_int(0)| = {f0_norm:.3e} (exactly zero)");
        assert_eq!(f0_norm, 0.0, "F_int must be exactly zero at u = 0");

        // (3) K_t(0) == K_0 to round-off.
        let kt0 = n_gamma_kt_local(&pre, &state, &zero);
        let mut kt0_dev = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                if kt0[(i, j)] != k0[(i, j)] {
                    kt0_dev = kt0_dev.max((kt0[(i, j)] - k0[(i, j)]).abs());
                }
            }
        }
        println!(
            "N-gamma (3) K_t(0) vs K_0: max|difference| = {kt0_dev:.3e} (relative {:.3e})",
            kt0_dev / k0_scale
        );
        // `K_t(0)` is the exact u = 0 limit of `dF/du`, `int B(0)^T W11 B(0)`,
        // NOT an additive `compute_ke_local` reference (see `n_gamma_kt_local`).
        // On the flat RECT that limit reproduces the linear stiffness to
        // round-off; a bit-for-bit equality was a construction artefact of the
        // removed additive reference and is not a paper requirement.
        assert!(
            kt0_dev / k0_scale < 1.0e-9,
            "K_t(0) must equal K_0 to a tight relative bound, got {:.3e}",
            kt0_dev / k0_scale
        );

        // (0) Diagnostic: does the 11-component ZERO-STATE assembly,
        // `int B(0)^T W11 B(0)`, reproduce the nine-row linear stiffness
        // `compute_ke_local(pre)`? If it does not, `n_gamma_kt_local`'s additive
        // `compute_ke_local(pre)` reference injects the spurious difference
        // `compute_ke_local - mat0` and the returned matrix cannot be the
        // derivative of `n_gamma_fint_local` — which is what the (2) check shows.
        let w11 = n_gamma_w11(&pre);
        let mut mat0 = Mat24::zeros();
        for g in 0..N_GAUSS {
            let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
            let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
            let b0 = n_gamma_b_matrix(&pre, &state, &zero, r, s);
            mat0 += (b0.transpose() * w11 * b0) * wq;
        }
        let mut dm = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                dm = dm.max((mat0[(i, j)] - k0[(i, j)]).abs());
            }
        }
        println!(
            "N-gamma (0) max|B(0)^T W11 B(0) - compute_ke_local| = {dm:.3e} (scale {k0_scale:.3e}, relative {:.3e})",
            dm / k0_scale
        );
        // C&S 185 is the 2017 MITC4+ formulation, which has NO drilling DOF, while
        // the production `compute_ke_local` is the 2025 MITC4+/D (the drill folded
        // into the membrane, Eq. (22a)). The correct u = 0 reference for the
        // N-gamma block is therefore the NO-DRILL linear stiffness; comparing
        // against the drilled one conflates "the drill is missing" with "the
        // 11-component assembly is wrong".
        let k0_nodrill = compute_ke_local_with_drill(&pre, false);
        let k0n_scale = max_abs(&k0_nodrill).max(1e-30);
        let mut dmn = 0.0f64;
        for i in 0..24 {
            for j in 0..24 {
                dmn = dmn.max((mat0[(i, j)] - k0_nodrill[(i, j)]).abs());
            }
        }
        println!(
            "N-gamma (0b) max|B(0)^T W11 B(0) - compute_ke_local_with_drill(false)| = {dmn:.3e} (scale {k0n_scale:.3e}, relative {:.3e})",
            dmn / k0n_scale
        );
        println!(
            "N-gamma (0c) max|compute_ke_local - compute_ke_local_with_drill(false)| = {:.3e} (the drill block's size)",
            {
                let mut d = 0.0f64;
                for i in 0..24 {
                    for j in 0..24 {
                        d = d.max((k0[(i, j)] - k0_nodrill[(i, j)]).abs());
                    }
                }
                d
            }
        );

        let axes: [([f64; 3], &str); 3] = [
            ([1.0, 0.0, 0.0], "e1 (in-plane)"),
            ([0.0, 0.0, 1.0], "e3 (element normal)"),
            (
                [
                    0.577_350_269_189_625_8,
                    0.577_350_269_189_625_8,
                    0.577_350_269_189_625_8,
                ],
                "skew (1,1,1)/sqrt3",
            ),
        ];
        let angles_deg = [1.0f64, 10.0, 45.0, 90.0];

        println!(
            "{:>20}  {:>8}  {:>16}  {:>16}  {:>16}  {:>14}  {:>14}",
            "axis",
            "phi[deg]",
            "faithful |F|/Ku",
            "prod |F|/Ku",
            "bounded |F|/Ku",
            "max|Kt-FD|",
            "rel"
        );

        for &(axis, label) in &axes {
            for &deg in &angles_deg {
                let phi = deg.to_radians();
                let r_mat = rotation_matrix(axis, phi);
                let theta = phi * Vector3::new(axis[0], axis[1], axis[2]);

                let mut u = Vec24::zeros();
                for i in 0..4 {
                    let x = Vector3::new(
                        state.coords[i][0],
                        state.coords[i][1],
                        state.coords[i][2],
                    );
                    let urel = r_mat * x - x;
                    u[6 * i] = urel[0];
                    u[6 * i + 1] = urel[1];
                    u[6 * i + 2] = urel[2];
                    u[6 * i + 3] = theta[0];
                    u[6 * i + 4] = theta[1];
                    u[6 * i + 5] = theta[2];
                }
                let mut u_buf = [0.0f64; 24];
                for i in 0..24 {
                    u_buf[i] = u[i];
                }
                let u_norm = vec_norm(&u_buf).max(1e-30);

                let f_faithful = n_gamma_fint_local(&pre, &state, &u);
                let mut ff_buf = [0.0f64; 24];
                for i in 0..24 {
                    ff_buf[i] = f_faithful[i];
                }
                let f_production = compute_fint_global(&pre, &u, true);
                let mut fp_buf = [0.0f64; 24];
                for i in 0..24 {
                    fp_buf[i] = f_production[i];
                }
                let f_bounded = compute_fint_global_bounded_superseded(&pre, &u);
                let mut fb_buf = [0.0f64; 24];
                for i in 0..24 {
                    fb_buf[i] = f_bounded[i];
                }
                let rel_faithful = vec_norm(&ff_buf) / (k0_scale * u_norm);
                let rel_production = vec_norm(&fp_buf) / (k0_scale * u_norm);
                let rel_bounded = vec_norm(&fb_buf) / (k0_scale * u_norm);

                // (2) K_t = dF/du by central differences.
                const H: f64 = 1.0e-6;
                let kt = n_gamma_kt_local(&pre, &state, &u);
                let mut fd = Mat24::zeros();
                for j in 0..24 {
                    let mut up = u;
                    up[j] += H;
                    let mut um = u;
                    um[j] -= H;
                    let fp = n_gamma_fint_local(&pre, &state, &up);
                    let fm = n_gamma_fint_local(&pre, &state, &um);
                    for i in 0..24 {
                        fd[(i, j)] = (fp[i] - fm[i]) / (2.0 * H);
                    }
                }
                let mut dmax = 0.0f64;
                let mut smax = 0.0f64;
                for i in 0..24 {
                    for j in 0..24 {
                        dmax = dmax.max((kt[(i, j)] - fd[(i, j)]).abs());
                        smax = smax.max(fd[(i, j)].abs());
                    }
                }
                let rel_fd = dmax / smax.max(1e-30);

                println!(
                    "{label:>20}  {deg:>8.1}  {rel_faithful:>16.3e}  {rel_production:>16.3e}  {rel_bounded:>16.3e}  {dmax:>14.3e}  {rel_fd:>14.3e}"
                );

                // (2) is the gate: the printed tangent must be the derivative of
                // the printed force (central differences are O(H^2) exact).
                assert!(
                    rel_fd < 1.0e-6,
                    "{label} {deg} deg: K_t differs from dF/du by {rel_fd:.3e}"
                );
                // The PRODUCTION nonlinear path (N-delta, option B) is the
                // faithful pair, so its rigid-body force residual must be at
                // round-off (the superseded bounded path is only printed).
                assert!(
                    rel_production < 1.0e-10,
                    "{label} {deg} deg: the production nonlinear path rigid-body residual {rel_production:.3e} is not round-off"
                );
            }
        }
        // ------------------------------------------------------------------
        // (4) GENERAL-increment consistency. On a rigid motion `S = W e~` is
        // zero, so the geometric term `int (dB/du)^T S d0V` is invisible and
        // the table above cannot discriminate it. `K_t = dF/du` must also hold
        // for a general random increment, and the relative residual must NOT
        // grow with `||u||` because `S` is proportional to the strain. This is
        // the identity the production Newton line search exercises.
        // ------------------------------------------------------------------
        let fd_rel = |k: &Mat24, fd: &Mat24| -> f64 {
            let mut d = 0.0f64;
            for i in 0..24 {
                for j in 0..24 {
                    d = d.max((k[(i, j)] - fd[(i, j)]).abs());
                }
            }
            d / max_abs(fd).max(1e-30)
        };
        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (((seed >> 11) as f64) / ((1u64 << 53) as f64)) * 2.0 - 1.0
        };
        println!("\nN-gamma (4) general-increment consistency (initial state, RECT):");
        println!("{:>8}  {:>16}  {:>16}", "|u|", "local rel", "global rel");
        for &scale in &[1.0e-4f64, 1.0e-3, 1.0e-2, 1.0e-1] {
            let mut u = Vec24::zeros();
            for x in u.iter_mut() {
                *x = scale * rnd();
            }
            const HG: f64 = 1.0e-6;
            let kt = n_gamma_kt_local(&pre, &state, &u);
            let fd = fd_matrix(|x| n_gamma_fint_local(&pre, &state, x), &u, HG);
            let ktg = n_gamma_kt_global(&pre, &u);
            let fdg = fd_matrix(|x| n_gamma_fint_global(&pre, x), &u, HG);
            let rl = fd_rel(&kt, &fd);
            let rg = fd_rel(&ktg, &fdg);
            println!("{scale:>8.0e}  {rl:>16.3e}  {rg:>16.3e}");
            assert!(
                rl < 1.0e-6,
                "general increment |u|={scale:.0e}: local K_t != dF/du (rel {rl:.3e})"
            );
            assert!(
                rg < 1.0e-6,
                "general increment |u|={scale:.0e}: global K_t != dF/du (rel {rg:.3e})"
            );
        }
        println!("(\"prod |F|/Ku\" is `compute_fint_global(.., true)`, the production nonlinear path; \"bounded |F|/Ku\" is the superseded `compute_fint_global_bounded_superseded(..)`, the recorded historical before.)");
    }

    /// Central-difference Jacobian of a force function. Diagnostic only.
    fn fd_matrix<F: Fn(&Vec24) -> Vec24>(f: F, u: &Vec24, h: f64) -> Mat24 {
        let mut fd = Mat24::zeros();
        for j in 0..24 {
            let mut up = *u;
            up[j] += h;
            let mut um = *u;
            um[j] -= h;
            let fp = f(&up);
            let fm = f(&um);
            for i in 0..24 {
                fd[(i, j)] = (fp[i] - fm[i]) / (2.0 * h);
            }
        }
        fd
    }




    /// Largest absolute entry of an 11x24 operator (diagnostic helper).
    fn bl_matrix_max(b: &SMatrix<f64, 11, 24>) -> f64 {
        let mut m = 0.0f64;
        for a in 0..11 {
            for j in 0..24 {
                m = m.max(b[(a, j)].abs());
            }
        }
        m
    }

    /// DIAGNOSTIC (temporary): decompose the `K_t - dF/du` residual of
    /// `n_gamma_kt_local` into its `mat`/`geo` parts and test whether the
    /// geometric term's nested finite difference is the source.
    #[test]
    #[ignore = "diagnostic: K_t - dF/du decomposition"]
    fn n_gamma_kt_residual_diagnostic() {
        const BENT: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.01],
            [2.0, 1.0, 0.025],
            [1.0, 1.0, 0.011],
        ];
        for (name, c) in [("RECT", &RECT), ("BENT", &BENT)] {
            let pre = pre_from(c);
            let state = GlCurrentState {
                coords: pre.initial_coords_3d,
                vn: pre.vn,
                v1: pre.v1,
                v2: pre.v2,
                a_i: pre.a_i,
            };
            let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
            let mut rnd = move || {
                seed = seed
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (((seed >> 11) as f64) / ((1u64 << 53) as f64)) * 2.0 - 1.0
            };
            for &scale in &[1.0e-3f64, 0.1] {
                let mut u = Vec24::zeros();
                for x in u.iter_mut() {
                    *x = scale * rnd();
                }
                let w = n_gamma_w11(&pre);
                let t24 = build_t24(&pre);
                let u_local = t24 * u;
                // Replicate n_gamma_kt_local with a configurable outer step
                // (local coordinates; the caller maps to global).
                let kt_with = |h_outer: f64| -> (Mat24, Mat24, Mat24) {
                    let mut mat = Mat24::zeros();
                    let mut geo = Mat24::zeros();
                    for g in 0..N_GAUSS {
                        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                        let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
                        let mut e = nalgebra::SVector::<f64, 11>::zeros();
                        for (a, v) in
                            super::n_gamma_local_strain(&pre, &state, &u_local, r, s).iter().enumerate()
                        {
                            e[a] = *v;
                        }
                        let b = n_gamma_b_matrix(&pre, &state, &u_local, r, s);
                        mat += (b.transpose() * w * b) * wq;
                        let svec = w * e;
                        for j in 0..24 {
                            let mut up = u_local;
                            up[j] += h_outer;
                            let mut um = u_local;
                            um[j] -= h_outer;
                            // NOTE: both nested differences use the SAME step, which
                            // is what production does (`N_GAMMA_GEO_H`); see the
                            // constant's docstring for the measurement.
                            let bp = n_gamma_b_matrix_at(&pre, &state, &up, r, s, h_outer);
                            let bm = n_gamma_b_matrix_at(&pre, &state, &um, r, s, h_outer);
                            for a in 0..11 {
                                for i in 0..24 {
                                    geo[(i, j)] +=
                                        0.5 * (bp[(a, i)] - bm[(a, i)]) / h_outer * svec[a] * wq;
                                }
                            }
                        }
                    }
                    (mat + geo, mat, geo)
                };
                let fd = fd_matrix(|x| n_gamma_fint_global(&pre, x), &u, 1.0e-6);
                let fds = max_abs(&fd).max(1e-30);
                print!("  {name} scale={scale:.0e}: rel(kt(h), fd) ");
                for &ho in &[1.0e-6f64, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4] {
                    let (kt, _m, _g) = kt_with(ho);
                    let kt = t24.transpose() * kt * t24;
                    let mut d = 0.0f64;
                    for i in 0..24 {
                        for j in 0..24 {
                            d = d.max((kt[(i, j)] - fd[(i, j)]).abs());
                        }
                    }
                    print!("h={ho:.0e}:{:.2e} ", d / fds);
                }
                println!();
                // PRODUCTION: the acceptance metric, `|n_gamma_kt_local - dF/du|`
                // in LOCAL coordinates with a clean reference step. The replica
                // above is the same construction written out; this line is the
                // real element.
                let kt_prod = n_gamma_kt_local(&pre, &state, &u_local);
                let fd_loc =
                    fd_matrix(|x| n_gamma_fint_local(&pre, &state, x), &u_local, 1.0e-6);
                let mut dprot = 0.0f64;
                let mut atprot = (0usize, 0usize);
                for i in 0..24 {
                    for j in 0..24 {
                        let v = (kt_prod[(i, j)] - fd_loc[(i, j)]).abs();
                        if v > dprot {
                            dprot = v;
                            atprot = (i, j);
                        }
                    }
                }
                println!(
                    "    PRODUCTION rel(n_gamma_kt_local, dF/du local) = {:.3e} at {atprot:?}",
                    dprot / max_abs(&fd_loc).max(1e-30)
                );
                let (_, mat, geo) = kt_with(1.0e-6);
                let mat = t24.transpose() * mat * t24;
                let geo = t24.transpose() * geo * t24;
                println!(
                    "    max|mat|={:.3e} max|geo|={:.3e} max|fd|={:.3e}",
                    max_abs(&mat),
                    max_abs(&geo),
                    max_abs(&fd)
                );
                // Direct: is the B operator the derivative of the strain?
                let (r, s) = (GAUSS_XI[0], GAUSS_ETA[0]);
                let b2 = n_gamma_b_matrix(&pre, &state, &u_local, r, s);
                let mut dfdu = SMatrix::<f64, 11, 24>::zeros();
                for j in 0..24 {
                    let mut up = u_local;
                    up[j] += 1.0e-6;
                    let mut um = u_local;
                    um[j] -= 1.0e-6;
                    let ep = super::n_gamma_local_strain(&pre, &state, &up, r, s);
                    let em = super::n_gamma_local_strain(&pre, &state, &um, r, s);
                    for a in 0..11 {
                        dfdu[(a, j)] = 0.5 * (ep[a] - em[a]) / 1.0e-6;
                    }
                }
                let mut db = 0.0f64;
                for a in 0..11 {
                    for j in 0..24 {
                        db = db.max((b2[(a, j)] - dfdu[(a, j)]).abs());
                    }
                }
                println!(
                    "    g0: max|B(2e-5) - de/du(1e-6)|/max|B| = {:.3e}",
                    db / bl_matrix_max(&b2).max(1e-30)
                );
                // Predicted residual `sum_g B^T W (B - de/du) wq` at every
                // Gauss point, plus the actual `kt - fd` and its argmax.
                let mut pred = Mat24::zeros();
                for g in 0..N_GAUSS {
                    let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                    let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
                    let bg = n_gamma_b_matrix(&pre, &state, &u_local, r, s);
                    let mut dg = SMatrix::<f64, 11, 24>::zeros();
                    for j in 0..24 {
                        let mut up = u_local;
                        up[j] += 1.0e-6;
                        let mut um = u_local;
                        um[j] -= 1.0e-6;
                        let ep = super::n_gamma_local_strain(&pre, &state, &up, r, s);
                        let em = super::n_gamma_local_strain(&pre, &state, &um, r, s);
                        for a in 0..11 {
                            dg[(a, j)] = 0.5 * (ep[a] - em[a]) / 1.0e-6;
                        }
                    }
                    let mut dmax = 0.0f64;
                    for a in 0..11 {
                        for j in 0..24 {
                            dmax = dmax.max((bg[(a, j)] - dg[(a, j)]).abs());
                        }
                    }
                    print!("    g{g} max|B-de/du|={dmax:.2e} ");
                    pred += (bg.transpose() * w * (bg - dg)) * wq;
                }
                println!();
                let (kt, _m, _gg) = kt_with(1.0e-6);
                let kt = t24.transpose() * kt * t24;
                let mut actual = 0.0f64;
                let mut ai = (0usize, 0usize);
                let mut pi = (0usize, 0usize);
                let mut pmax = 0.0f64;
                for i in 0..24 {
                    for j in 0..24 {
                        let d = (kt[(i, j)] - fd[(i, j)]).abs();
                        if d > actual {
                            actual = d;
                            ai = (i, j);
                        }
                        let dp = pred[(i, j)].abs();
                        if dp > pmax {
                            pmax = dp;
                            pi = (i, j);
                        }
                    }
                }
                println!(
                    "    max|kt-fd|={actual:.2e} at {ai:?}; max|pred|={pmax:.2e} at {pi:?}; rel(actual)={:.2e}",
                    actual / fds
                );
            }
        }
    }

    /// DIAGNOSTIC (temporary): the REFERENCE-FREE error of the N-gamma tangent,
    /// and the geometry class that error belongs to.
    ///
    /// The tangent is a Hessian, so `geo` must be symmetric in `(i, j)` in every
    /// configuration. The antisymmetric part `X - X^T` is therefore pure error,
    /// and it separates the two candidate causes with no finite-difference
    /// reference at all:
    ///
    /// - ROUND-OFF of the nested difference: the two index orders evaluate
    ///   different four-point stencils (`H_i` and `H_o` swap roles), so their
    ///   round-offs are uncorrelated and `|geo - geo^T|` sits at the noise floor
    ///   `~ eps/(H_i H_o) max|W e|`.
    /// - A MISSING OR WRONG SYMMETRIC TERM: it contributes equally to both index
    ///   orders and leaves `|geo - geo^T|` at the floor while `|kt - dF/du|`
    ///   stays large.
    ///
    /// The fixtures separate the other two candidate causes, non-planarity
    /// (`x_d` out of the mid-surface, the `zeta`/director path) and the assumed
    /// membrane coefficients `a_A..a_E` of Eq. (15e) (identically zero only when
    /// `x_d = 0`):
    ///
    /// ```text
    /// fixture          planar  x_d   a_A..a_E
    /// RECT             yes     0     0
    /// ROT_RECT(60deg)  yes     0     0
    /// FLAT_DISTORTED   yes     != 0  != 0
    /// STRONGLY_WARPED  no      != 0  != 0
    /// BENT             no      != 0  != 0
    /// ```
    #[test]
    #[ignore = "diagnostic: reference-free tangent symmetry by geometry"]
    fn n_gamma_geo_symmetry_diagnostic() {
        fn asym(m: &Mat24) -> (f64, f64, (usize, usize)) {
            let mut best = 0.0f64;
            let mut at = (0usize, 0usize);
            let mut scale = 0.0f64;
            for i in 0..24 {
                for j in 0..24 {
                    scale = scale.max(m[(i, j)].abs());
                    let d = (m[(i, j)] - m[(j, i)]).abs();
                    if d > best {
                        best = d;
                        at = (i, j);
                    }
                }
            }
            (best, scale, at)
        }

        const ROT_RECT: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.866_025_403_784_438_6, 1.0],
            [0.0, 0.866_025_403_784_438_6, 1.0],
        ];
        const BENT: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.01],
            [2.0, 1.0, 0.025],
            [1.0, 1.0, 0.011],
        ];

        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (((seed >> 11) as f64) / ((1u64 << 53) as f64)) * 2.0 - 1.0
        };
        let mut cross_checked = false;

        for (name, c) in [
            ("RECT", &RECT),
            ("ROT_RECT(60)", &ROT_RECT),
            ("FLAT_DISTORTED", &FLAT_DISTORTED),
            ("STRONGLY_WARPED", &STRONGLY_WARPED),
            ("BENT", &BENT),
        ] {
            let pre = pre_from(c);
            let state = GlCurrentState {
                coords: pre.initial_coords_3d,
                vn: pre.vn,
                v1: pre.v1,
                v2: pre.v2,
                a_i: pre.a_i,
            };
            let (cur_vecs, _n, m_r, m_s) = gl_current_characteristic_vectors(&state);
            let (_cr, _cs, d, coeff) =
                compute_membrane_coefficients_2017(&cur_vecs[2], &m_r, &m_s);
            println!(
                "\n=== {name}: |x_d|/|x_r|={:.3e} d={:.3e} max|a|={:.3e}",
                cur_vecs[2].norm() / cur_vecs[0].norm(),
                d,
                coeff.iter().fold(0.0f64, |m, v| m.max(v.abs()))
            );
            for &scale in &[1.0e-3f64, 1.0e-1] {
                // Local-frame increment, as the other N-gamma diagnostics use.
                let mut u = Vec24::zeros();
                for x in u.iter_mut() {
                    *x = scale * rnd();
                }
                // Production steps: both nested differences use the SAME step
                // (`N_GAMMA_GEO_H`), so the replica must too or it stops being a
                // replica -- the cross-check below catches it if it drifts.
                const HO: f64 = N_GAMMA_GEO_H;
                let w = n_gamma_w11(&pre);
                let mut mat = Mat24::zeros();
                let mut geo = Mat24::zeros();
                let mut we_max = 0.0f64;
                for g in 0..N_GAUSS {
                    let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                    let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
                    let mut e = nalgebra::SVector::<f64, 11>::zeros();
                    for (a, v) in
                        super::n_gamma_local_strain(&pre, &state, &u, r, s).iter().enumerate()
                    {
                        e[a] = *v;
                    }
                    let b = n_gamma_b_matrix(&pre, &state, &u, r, s);
                    mat += (b.transpose() * w * b) * wq;
                    let svec = w * e;
                    for a in 0..11 {
                        we_max = we_max.max(svec[a].abs());
                    }
                    for j in 0..24 {
                        let mut up = u;
                        up[j] += HO;
                        let mut um = u;
                        um[j] -= HO;
                        let bp = n_gamma_b_matrix_at(&pre, &state, &up, r, s, HO);
                        let bm = n_gamma_b_matrix_at(&pre, &state, &um, r, s, HO);
                        for a in 0..11 {
                            for i in 0..24 {
                                geo[(i, j)] += 0.5 * (bp[(a, i)] - bm[(a, i)]) / HO * svec[a] * wq;
                            }
                        }
                    }
                }
                let kt = mat + geo;
                if !cross_checked {
                    let kt_prod = n_gamma_kt_local(&pre, &state, &u);
                    let mut dd = 0.0f64;
                    for i in 0..24 {
                        for j in 0..24 {
                            dd = dd.max((kt_prod[(i, j)] - kt[(i, j)]).abs());
                        }
                    }
                    println!("    replication cross-check vs n_gamma_kt_local: {dd:.3e}");
                    cross_checked = true;
                }
                let (ak, sk, atk) = asym(&kt);
                let (am, sm, _) = asym(&mat);
                let (ag, sg, atg) = asym(&geo);
                // Round-off scale of the four-point stencil: eps / h^2, times the
                // stress-weighted strain magnitude.
                let floor = f64::EPSILON / (HO * HO) * we_max;
                println!(
                    "  scale={scale:.0e}: max|We|={we_max:.3e} max|kt|={sk:.3e} max|mat|={sm:.3e} max|geo|={sg:.3e}"
                );
                println!(
                    "    |kt-kt^T|/|kt|={:.3e} at {atk:?} | |mat-mat^T|/|mat|={:.3e} | |geo-geo^T|/|geo|={:.3e} at {atg:?}",
                    ak / sk.max(1e-300),
                    am / sm.max(1e-300),
                    ag / sg.max(1e-300)
                );
                println!(
                    "    |geo-geo^T| / [eps/h^2*max|We|] = {:.3e}",
                    ag / floor.max(1e-300)
                );
            }
        }
    }

    /// DIAGNOSTIC (temporary): does a SINGLE-STEP four-point stencil for the
    /// geometric term reach the accuracy the nested difference cannot?
    ///
    /// `geo` is the Hessian `sum_a (d2 _0 e~_a / du_i du_j) S_a`, which the
    /// production code estimates as `FD(FD(e))` with an inner step `H_i = 2e-5`
    /// (set by the `B(0)` accuracy) and an outer step `H_o = 1e-6`. That
    /// composition has the round-off floor `eps / (H_i H_o)` and the two steps
    /// cannot both be relaxed, so the floor is what the previous instrument
    /// measured. The same Hessian can be estimated by the single-step four-point
    /// cross difference
    ///
    /// ```text
    /// d2 e / du_i du_j ~ [e(u + h e_i + h e_j) - e(u + h e_i - h e_j)
    ///                     - e(u - h e_i + h e_j) + e(u - h e_i - h e_j)] / (4 h^2),
    /// ```
    ///
    /// whose round-off is `eps / h^2` and whose truncation is `O(h^2)`, balanced
    /// at the classic `h ~ eps^(1/4) ~ 1e-4` with total error `~ eps^(1/2) ~ 1e-8`
    /// -- two orders below the `1e-6` the consistency gate asks for. The stencil
    /// is also EXACTLY symmetric in `(i, j)` (both index orders combine the same
    /// four points), so it is measurable with no reference at all through
    /// `|X - X^T| / |X|`.
    #[test]
    #[ignore = "diagnostic: geometric-term stencil accuracy probe"]
    fn n_gamma_geo_stencil_probe() {
        fn asym(m: &Mat24) -> (f64, f64) {
            let mut best = 0.0f64;
            let mut scale = 0.0f64;
            for i in 0..24 {
                for j in 0..24 {
                    scale = scale.max(m[(i, j)].abs());
                    best = best.max((m[(i, j)] - m[(j, i)]).abs());
                }
            }
            (best, scale)
        }
        fn rel(a: &Mat24, b: &Mat24) -> (f64, f64) {
            let mut d = 0.0f64;
            let mut s = 0.0f64;
            for i in 0..24 {
                for j in 0..24 {
                    d = d.max((a[(i, j)] - b[(i, j)]).abs());
                    s = s.max(b[(i, j)].abs());
                }
            }
            (d, s)
        }

        const ROT_RECT: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 0.866_025_403_784_438_6, 1.0],
            [0.0, 0.866_025_403_784_438_6, 1.0],
        ];
        const BENT: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.01],
            [2.0, 1.0, 0.025],
            [1.0, 1.0, 0.011],
        ];

        let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rnd = move || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (((seed >> 11) as f64) / ((1u64 << 53) as f64)) * 2.0 - 1.0
        };

        for (name, c) in [
            ("RECT", &RECT),
            ("ROT_RECT(60)", &ROT_RECT),
            ("FLAT_DISTORTED", &FLAT_DISTORTED),
            ("BENT", &BENT),
        ] {
            let pre = pre_from(c);
            let state = GlCurrentState {
                coords: pre.initial_coords_3d,
                vn: pre.vn,
                v1: pre.v1,
                v2: pre.v2,
                a_i: pre.a_i,
            };
            let mut u = Vec24::zeros();
            for x in u.iter_mut() {
                *x = 1.0e-1 * rnd();
            }
            let w = n_gamma_w11(&pre);
            let hs: &[f64] = if name == "BENT" {
                &[3.0e-5, 1.0e-4, 3.0e-4]
            } else {
                &[1.0e-4]
            };
            let mut mat = Mat24::zeros();
            let mut geo_nested = Mat24::zeros();
            for g in 0..N_GAUSS {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
                let mut e = nalgebra::SVector::<f64, 11>::zeros();
                for (a, v) in
                    super::n_gamma_local_strain(&pre, &state, &u, r, s).iter().enumerate()
                {
                    e[a] = *v;
                }
                let b = n_gamma_b_matrix(&pre, &state, &u, r, s);
                mat += (b.transpose() * w * b) * wq;
                let svec = w * e;
                const HO: f64 = 1.0e-6;
                for j in 0..24 {
                    let mut up = u;
                    up[j] += HO;
                    let mut um = u;
                    um[j] -= HO;
                    let bp = n_gamma_b_matrix(&pre, &state, &up, r, s);
                    let bm = n_gamma_b_matrix(&pre, &state, &um, r, s);
                    for a in 0..11 {
                        for i in 0..24 {
                            geo_nested[(i, j)] +=
                                0.5 * (bp[(a, i)] - bm[(a, i)]) / HO * svec[a] * wq;
                        }
                    }
                }
            }
            let (an, sn) = asym(&geo_nested);
            println!(
                "\n=== {name}: scale=1e-1  max|mat|={:.3e} max|geo_nest|={sn:.3e}",
                mat.abs().max()
            );
            println!(
                "    superseded split-step (H_i=2e-5, H_o=1e-6): |geo-geo^T|/|geo|={:.3e}",
                an / sn.max(1e-300)
            );
            for &h in hs {
                let mut g4 = Mat24::zeros();
                for g in 0..N_GAUSS {
                    let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                    let wq = GAUSS_W[g] * surface_measure(&pre, r, s);
                    let mut e = nalgebra::SVector::<f64, 11>::zeros();
                    for (a, v) in
                        super::n_gamma_local_strain(&pre, &state, &u, r, s).iter().enumerate()
                    {
                        e[a] = *v;
                    }
                    let svec = w * e;
                    for i in 0..24 {
                        for j in i..24 {
                            let mut pp = u;
                            pp[i] += h;
                            pp[j] += h;
                            let mut pm = u;
                            pm[i] += h;
                            pm[j] -= h;
                            let mut mp = u;
                            mp[i] -= h;
                            mp[j] += h;
                            let mut mm = u;
                            mm[i] -= h;
                            mm[j] -= h;
                            let epp = super::n_gamma_local_strain(&pre, &state, &pp, r, s);
                            let epm = super::n_gamma_local_strain(&pre, &state, &pm, r, s);
                            let emp = super::n_gamma_local_strain(&pre, &state, &mp, r, s);
                            let emm = super::n_gamma_local_strain(&pre, &state, &mm, r, s);
                            let mut acc = 0.0f64;
                            for a in 0..11 {
                                acc += ((epp[a] - epm[a] - emp[a] + emm[a]) / (4.0 * h * h))
                                    * svec[a];
                            }
                            g4[(i, j)] += acc * wq;
                            if i != j {
                                g4[(j, i)] += acc * wq;
                            }
                        }
                    }
                }
                let (a4, s4) = asym(&g4);
                let (d, sc) = rel(&g4, &geo_nested);
                println!(
                    "    4-pt h={h:.0e}          : |geo-geo^T|/|geo|={:.3e} | |g4-gn|/|gn|={:.3e} (|g4|={s4:.3e}, |gn|={sc:.3e})",
                    a4 / s4.max(1e-300),
                    d / sc.max(1e-300)
                );
            }
        }
    }
}
