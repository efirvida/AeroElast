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

use nalgebra::{Matrix2, Matrix3, SMatrix, SVector, Vector3, Vector4};

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
/// local shape-function derivatives and the local 24-vector.
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
/// path is paper-faithful and this correction is bounded and oracle-checked by
/// the T2B consistency set.
fn membrane_strain_nl(h: &Matrix3<f64>) -> Vector3<f64> {
    Vector3::new(
        0.5 * (h[(0, 0)].powi(2) + h[(1, 0)].powi(2) + h[(2, 0)].powi(2)),
        0.5 * (h[(0, 1)].powi(2) + h[(1, 1)].powi(2) + h[(2, 1)].powi(2)),
        h[(0, 0)] * h[(0, 1)] + h[(1, 0)] * h[(1, 1)] + h[(2, 0)] * h[(2, 1)],
    )
}

/// Derivative of [`membrane_strain_nl`] with respect to the local 24-vector,
/// `B_nl` (6x24; rows 0, 1, 3 carry `[E_xx, E_yy, 2 E_xy]`).
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

/// Extract the membrane rows `[0, 1, 3]` from a 6x24 matrix into a 3x24.
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

/// Initial-stress stiffness `K_sigma` in LOCAL coordinates from the current
/// displacement state: the membrane resultant `sigma = cm (B_m u)` is evaluated
/// at every Gauss point and contracted with `B_geo`.
fn geometric_stiffness_local(pre: &Mitc4PlusDPrecomputed, u_local: &Vec24) -> Mat24 {
    let cm = &pre.constitutive.cm;
    let mut k = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let bm = b_membrane_2017(pre, r, s);
        let sigma_g = cm * (bm * u_local);
        k += geometric_stiffness_contribution(pre, g, &sigma_g);
    }
    0.5 * (k + k.transpose())
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

/// The bounded nonlinear membrane correction to the internal force:
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
/// correction is carried. The T2B consistency set
/// (`test_kt_zero_matches_ke`, `test_fint_linear_nonlinear_parity`, the
/// directional-derivative triple) is its oracle; a paper-faithful nonlinear
/// derivation is explicitly deferred, not invented.
fn membrane_nonlinear_correction(pre: &Mitc4PlusDPrecomputed, u_local: &Vec24) -> Vec24 {
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

/// Compute the internal force vector in GLOBAL coordinates (24 long).
///
/// `nonlinear = false` returns `K u` exactly. `nonlinear = true` adds the
/// bounded total-Lagrangian membrane correction of
/// [`membrane_nonlinear_correction`]; the bending, transverse-shear and drill
/// contributions stay linear (their geometric parts are carried by `K_sigma` of
/// [`compute_kt_global`]).
///
/// The nonlinear path is bounded: see the note on
/// [`membrane_nonlinear_correction`].
pub fn compute_fint_global(
    pre: &Mitc4PlusDPrecomputed,
    u_global: &Vec24,
    nonlinear: bool,
) -> Vec24 {
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;

    let f_local = if !nonlinear {
        compute_ke_local(pre) * u_local
    } else {
        compute_ke_local(pre) * u_local + membrane_nonlinear_correction(pre, &u_local)
    };

    t24.transpose() * f_local
}

/// Compute the tangent stiffness `K_T` in GLOBAL coordinates (24x24).
///
/// Total Lagrangian: `K_T = K_0 + K_L + K_sigma`, where `K_0` is the linear
/// stiffness, `K_L` the initial-displacement stiffness of the bounded membrane
/// correction and `K_sigma` the initial-stress (geometric) stiffness. At
/// `u = 0` both corrections vanish, so `K_T(0) = K_0` exactly.
///
/// The tangent is symmetrised in the local frame before the transformation, but
/// the transformation itself is **not** post-symmetrised, so `K_T(0)` equals
/// `T^T K_0 T` bit for bit (the `test_kt_zero_matches_ke` guard).
pub fn compute_kt_global(pre: &Mitc4PlusDPrecomputed, u_global: &Vec24) -> Mat24 {
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;

    let k0 = compute_ke_local(pre);
    let cm = &pre.constitutive.cm;

    // K_L: initial-displacement stiffness of the nonlinear membrane strain.
    let mut k_l = Mat24::zeros();
    for g in 0..N_GAUSS {
        let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
        let sqrt_g = surface_measure(pre, r, s);
        let dh = local_shape_derivatives(pre, r, s);
        let h_mat = displacement_gradient(&dh, &u_local);
        let bnl = compute_b_nl(&dh, &h_mat);
        let bm = b_membrane_2017(pre, r, s);
        let bm_nl = extract_membrane_rows(&bnl);
        k_l += (bm.transpose() * cm * bm_nl
            + bm_nl.transpose() * cm * bm
            + bm_nl.transpose() * cm * bm_nl)
            * (GAUSS_W[g] * sqrt_g);
    }

    let k_sigma = geometric_stiffness_local(pre, &u_local);

    let k_t = k0 + k_l + k_sigma;
    let mut k_t_sym = k_t;
    for i in 0..24 {
        for j in 0..24 {
            k_t_sym[(i, j)] = 0.5 * (k_t[(i, j)] + k_t[(j, i)]);
        }
    }
    t24.transpose() * k_t_sym * t24
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
// Corotational machinery (retargeted, formulation-independent)
// ============================================================================
//
// Ported from the repository's existing S4R-style machinery. These helpers are
// independent of the element formulation: they rotate the nodal directors and
// extract the logarithmic membrane strain for the corotational frame update.
// The large-rotation Python benchmarks (`tests/test_large_rotation_benchmarks.py`)
// are their oracle, judged at the S3 gate (task 10.7).

/// Local orthonormal frame at a Gauss point.
#[derive(Clone, Copy)]
pub struct GpLocalFrame {
    /// First in-plane tangent (e1).
    pub e1: Vector3<f64>,
    /// Second in-plane tangent (e2).
    pub e2: Vector3<f64>,
    /// Normal (e3).
    pub e3: Vector3<f64>,
}

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

    /// Compose quaternions: `q = q1 (x) q2`.
    pub fn quaternion_multiply(q1: &Vector4<f64>, q2: &Vector4<f64>) -> Vector4<f64> {
        let (q1_0, q1x, q1y, q1z) = (q1[0], q1[1], q1[2], q1[3]);
        let (q2_0, q2x, q2y, q2z) = (q2[0], q2[1], q2[2], q2[3]);
        Vector4::new(
            q1_0 * q2_0 - q1x * q2x - q1y * q2y - q1z * q2z,
            q1_0 * q2x + q1x * q2_0 + q1y * q2z - q1z * q2y,
            q1_0 * q2y - q1x * q2z + q1y * q2_0 + q1z * q2x,
            q1_0 * q2z + q1x * q2y - q1y * q2x + q1z * q2_0,
        )
    }

    /// Updated nodal directors from a local displacement increment. The initial
    /// directors are the element's own `V_n^i` (ADR-4 option B), not a
    /// placeholder field.
    pub fn update_normals_with_displacements(
        pre: &Mitc4PlusDPrecomputed,
        delta_u_local: &Vec24,
    ) -> [Vector3<f64>; 4] {
        let mut updated_normals: [Vector3<f64>; 4] = [Vector3::zeros(); 4];
        for i in 0..4 {
            let theta = Vector3::new(
                delta_u_local[6 * i + 3],
                delta_u_local[6 * i + 4],
                delta_u_local[6 * i + 5],
            );
            let q_inc = Self::quaternion_from_vector(&theta);
            updated_normals[i] = Self::rotate_vector_by_quaternion(&pre.vn[i], &q_inc);
        }
        updated_normals
    }

    /// Polar decomposition `F = R U` by Newton iteration on `U = sqrt(F^T F)`;
    /// returns `(R, U)`.
    pub fn polar_decomposition(h: &Matrix3<f64>) -> (Matrix3<f64>, Matrix3<f64>) {
        let f = Matrix3::identity() + h;
        let ct = f.transpose() * f;
        let mut u = 0.5 * (ct + Matrix3::identity());
        for _ in 0..3 {
            if let Some(u_inv) = u.try_inverse() {
                u = 0.5 * (u + ct * u_inv);
            } else {
                break;
            }
        }
        let r_inc = if let Some(u_inv) = u.try_inverse() {
            f * u_inv
        } else {
            Matrix3::identity()
        };
        (r_inc, u)
    }

    /// Logarithmic strain from the right stretch tensor `U` (small-strain
    /// approximation `ln(U) ~ U - I`).
    pub fn log_strain_from_polar(u: &Matrix3<f64>) -> Matrix3<f64> {
        u - Matrix3::identity()
    }

    /// Membrane strain in Voigt form from the logarithmic strain of the
    /// displacement gradient.
    pub fn compute_membrane_strain_log(h: &Matrix3<f64>) -> Vector3<f64> {
        let (_r_inc, u_inc) = Self::polar_decomposition(h);
        let eps_log = Self::log_strain_from_polar(&u_inc);
        Vector3::new(eps_log[(0, 0)], eps_log[(1, 1)], 2.0 * eps_log[(0, 1)])
    }

    /// Updated corotational frame at every Gauss point from the current nodal
    /// coordinates. Falls back to the reference frame at that point when the
    /// deformed tangents degenerate.
    pub fn update_corotational_frame(
        &self,
        current_coords: &[[f64; 3]; 4],
    ) -> [GpLocalFrame; N_GAUSS] {
        let mut updated_frames: [GpLocalFrame; N_GAUSS] = [GpLocalFrame {
            e1: Vector3::zeros(),
            e2: Vector3::zeros(),
            e3: Vector3::zeros(),
        }; N_GAUSS];

        for g in 0..N_GAUSS {
            let xi = GAUSS_XI[g];
            let eta = GAUSS_ETA[g];

            // Reference frame at this Gauss point (fallback).
            let (gr0, gs0, _) =
                compute_j3d_enriched(&self.initial_coords_3d, &self.vn, &self.a_i, xi, eta, 0.0);
            let n0 = gr0.cross(&gs0);
            let e3_0 = if n0.norm() > 1e-12 {
                n0.normalize()
            } else {
                self.e3
            };
            let mut e1_0 = gr0 - gr0.dot(&e3_0) * e3_0;
            if e1_0.norm() < 1e-12 {
                e1_0 = self.e1 - self.e1.dot(&e3_0) * e3_0;
            }
            let e1_0 = e1_0.normalize();
            let e2_0 = e3_0.cross(&e1_0).normalize();

            // Deformed tangents.
            let (g_r_def, g_s_def, _) =
                compute_j3d_enriched(current_coords, &self.vn, &self.a_i, xi, eta, 0.0);
            let n_cross = g_r_def.cross(&g_s_def);
            let e3_new = if n_cross.norm() > 1e-12 {
                n_cross.normalize()
            } else {
                e3_0
            };
            let g_r_proj = g_r_def - g_r_def.dot(&e3_new) * e3_new;
            let e1_new = if g_r_proj.norm() > 1e-12 {
                g_r_proj.normalize()
            } else {
                e1_0
            };
            let e2_new = if e3_new.norm() > 1e-12 && e1_new.norm() > 1e-12 {
                e3_new.cross(&e1_new).normalize()
            } else {
                e2_0
            };

            updated_frames[g] = GpLocalFrame {
                e1: e1_new,
                e2: e2_new,
                e3: e3_new,
            };
        }

        updated_frames
    }

    /// Incremental rotation from the old to the new corotational frame:
    /// `R_inc = R_new R_old^T`.
    pub fn frame_incremental_rotation(
        old_frame: &GpLocalFrame,
        new_frame: &GpLocalFrame,
    ) -> Matrix3<f64> {
        let r_old = Matrix3::from_columns(&[old_frame.e1, old_frame.e2, old_frame.e3]);
        let r_new = Matrix3::from_columns(&[new_frame.e1, new_frame.e2, new_frame.e3]);
        r_new * r_old.transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::{
        b_bending_2017, b_bending_covariant_2017, b_drill_membrane_2025, b_membrane_2017,
        b_membrane_covariant_2017, b_shear_mitc4, build_t24, compute_fint_global,
        compute_ke_global, compute_ke_local, compute_ke_local_with_drill, compute_kt_global,
        compute_me_composite_global, compute_me_global, compute_membrane_coefficients_2017,
        covariant_to_local_mapping, drill_jacobian_ratio, drill_ke_local,
        drill_midside_shape_derivatives, element_area, interpolate_displacement,
        interpolate_position, j_loc_at, membrane_ke_local, node_vec, resultant_moment_matrix,
        shape_function_derivatives, shape_functions, shear_ke_local, surface_measure, Mat24,
        Mitc4PlusDPrecomputed, Vec24, DRILL_EDGE_MID, GAUSS_ETA, GAUSS_W, GAUSS_XI, NODE_ETA,
        NODE_XI, N_GAUSS,
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
        assert!(
            diff.norm() < 1e-10,
            "K_T(u=0) must be identical to K_linear_global, diff norm = {}",
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

    /// The interior moment conjugate to the rotation DOFs for a constant
    /// transverse-shear resultant `q = [q13, q23]`:
    /// `f_theta_x,i = -integral N_i q23 dA`, `f_theta_y,i = +integral N_i q13 dA`
    /// (design 4.2 completed; see the shearing test's derivation).
    fn interior_shear_moment(pres: &[Mitc4PlusDPrecomputed; 5], q: [f64; 2]) -> DVector<f64> {
        let mut f = DVector::zeros(48);
        for (e, pe) in pres.iter().enumerate() {
            let el = STAR_ELEMS[e];
            for g in 0..GAUSS_XI.len() {
                let (r, s) = (GAUSS_XI[g], GAUSS_ETA[g]);
                let w = GAUSS_W[g] * surface_measure(pe, r, s);
                let n = shape_functions(r, s);
                for a in 0..4 {
                    f[6 * el[a] + 3] += w * n[a] * (-q[1]);
                    f[6 * el[a] + 4] += w * n[a] * q[0];
                }
            }
        }
        f
    }

    #[test]
    fn test_t1a_shearing_patch_constant_stress_fig5_mesh() {
        // Ko, Lee & Bathe (2017), "A new MITC4+ shell element", Computers and
        // Structures 182:404-418, Section 4, pp. 410-411, and the MITC4 assumed
        // transverse shear of Dvorkin & Bathe (1984), "A continuum mechanics
        // based four-node shell element for general nonlinear analysis",
        // Engineering Computations 1:77-88, Eq. (3), reproduced in Ko, Lee &
        // Bathe (2017), C&S 182:404-418, p. 405. No shear correction factor is
        // applied: the element consumes the uncorrected `G h` (Ko, Lee &
        // Bathe (2017), C&S 182:404-418, p. 410).
        //
        // DERIVATION OF THE LOAD (design 4.2, as written). The paper publishes
        // no load magnitudes. The design derives the boundary nodal forces of
        // the constant state as `f_i = contour_integral N_i (q . n) dGamma` with
        // `q = (q13, q23) = G h (gamma_13, gamma_23)`, integrated with 2-point
        // Gauss per boundary edge and computed element-independently from `q`
        // and the mesh (nothing uses the element stiffness).
        //
        // FINDING (recorded in `apply-progress.md`, WU6): that load does not
        // produce a constant transverse shear. A constant shear resultant is
        // not an equilibrium state of a Mindlin plate: the rotation rows of the
        // element's internal force `integral B_gamma^T q dA` are not balanced by
        // boundary tractions alone (the pointwise-mindlin completion of the
        // design's derivation adds `f_theta_x,i = -integral N_i q23 dA`,
        // `f_theta_y,i = +integral N_i q13 dA`, printed below as the
        // `complete` diagnostic, and it too does not recover the constant
        // state, because the assumed MITC4 operator's rotation rows are not the
        // pointwise ones). Both measured values are printed before the
        // assertion; the assertion is on the design's boundary-only load. The
        // design 4.2 derivation is therefore wrong for the shearing patch, and
        // the paper's own shearing patch (Ko, Bathe & Zhang (2025),
        // C&S 308:107622, Fig. 7(c): u_x constrained at the interior nodes, load
        // in +y at A) is an in-plane shear state, not a transverse one.
        let nodes = STAR_NODES;
        let e_mod = 2.0e11f64;
        let nu = 0.3f64;
        let gh = e_mod / (2.0 * (1.0 + nu)); // the uncorrected G h, h = 1
        let pres = star_patch_pres(&nodes);
        let k = assemble_star_patch(&pres);

        let states = [("gamma_13", [1.0e-3f64, 0.0]), ("gamma_23", [0.0, 1.0e-3])];
        for (name, gam) in states {
            let gamma = [gam[0], gam[1]];
            let q = [gh * gamma[0], gh * gamma[1]];

            let f_boundary = boundary_integrate(&nodes, |_x, n| {
                [0.0, 0.0, q[0] * n[0] + q[1] * n[1], 0.0, 0.0, 0.0]
            });
            let f_full = &f_boundary + interior_shear_moment(&pres, q);

            // The exact field: constant rotations `theta_x = -gamma_23`,
            // `theta_y = gamma_13`, `w = 0` (zero curvature, constant shear).
            let prescribed: Vec<(usize, f64)> = BC_2017_PATCH
                .iter()
                .map(|bc| {
                    let v = match bc.dof {
                        3 => -gamma[1],
                        4 => gamma[0],
                        _ => 0.0,
                    };
                    (6 * bc.node + bc.dof, v)
                })
                .collect();

            let u_full = solve_constrained(&k, &f_full, &prescribed);
            let u_boundary = solve_constrained(&k, &f_boundary, &prescribed);

            let gnorm = (gamma[0] * gamma[0] + gamma[1] * gamma[1]).sqrt();
            let floor = 1e-10 * gnorm;
            let measure = |u: &DVector<f64>| -> (f64, f64) {
                let mut err = 0.0f64;
                let mut vals: Vec<[f64; 2]> = Vec::new();
                for (e, pe) in pres.iter().enumerate() {
                    let ue = star_elem_disp(u, e);
                    let ul = build_t24(pe) * ue;
                    for g in 0..GAUSS_XI.len() {
                        let gl = b_shear_mitc4(pe, GAUSS_XI[g], GAUSS_ETA[g]) * ul;
                        let gg = rotate_shear_to_global(&pe.t3, &[gl[0], gl[1]]);
                        let d = ((gg[0] - gamma[0]).powi(2) + (gg[1] - gamma[1]).powi(2)).sqrt();
                        err = err.max(d);
                        vals.push(gg);
                    }
                }
                let spread = (0..2)
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
                (err, spread)
            };

            let (err_full, spread_full) = measure(&u_full);
            let (err_b, spread_b) = measure(&u_boundary);
            println!(
                "T1.3c {name}: complete load (boundary + interior moment): \
                 max|gamma_gp - gamma|={err_full:.3e} (rel {:.3e}); spread={spread_full:.3e} (rel {:.3e})",
                err_full / gnorm,
                spread_full / gnorm
            );
            println!(
                "T1.3c {name}: boundary-only (design 4.2 as written): \
                 max|gamma_gp - gamma|={err_b:.3e} (rel {:.3e}); spread={spread_b:.3e} (rel {:.3e})",
                err_b / gnorm,
                spread_b / gnorm
            );
            assert!(
                err_b <= 1e-8 * gnorm + floor,
                "{name}: recovered shear error {err_b:.3e} > 1e-8 ||gamma|| + floor ({:.3e})",
                1e-8 * gnorm + floor
            );
            assert!(
                spread_b <= 1e-8 * gnorm + floor,
                "{name}: recovered shear spread {spread_b:.3e} > 1e-8 ||gamma|| + floor ({:.3e})",
                1e-8 * gnorm + floor
            );
        }
    }
}
