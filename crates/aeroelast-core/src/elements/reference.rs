// Reference Element Traits
//
// Defines `ReferenceElement2D` plus generic Gauss integration helpers for
// stiffness (Ke) and consistent mass (Me).
//
// Design decisions:
// - Return `&'static [[f64; D]]` slices for Gauss points/weights — these are
//   static references to const arrays, zero heap allocation in the hot path.
// - Avoid GATs / const-generic arithmetic in return types — requires unstable
//   `generic_const_exprs`. Callers pass pre-allocated `&mut [f64]` slices,
//   which they stack-allocate using fixed-size arrays.
// - MITC3/MITC4 and composite elements are intentionally excluded — they have
//   special shear-correction terms that do not fit the plain Bᵀ·C·B pattern.

// ============================================================================
// 2D trait (plane stress / plane strain quad elements)
// ============================================================================

/// A 2D reference element defined on [−1,1]² (quads).
///
/// `NODES`: number of nodes (e.g. 4 for QUAD4, 8 for QUAD8, 9 for QUAD9).
pub trait ReferenceElement2D<const NODES: usize> {
    /// Shape function values at reference point (xi, eta).
    fn shape_functions(&self, xi: f64, eta: f64) -> [f64; NODES];

    /// Shape function derivatives w.r.t. (xi, eta).
    /// Returns (dN/dxi, dN/deta) each of length NODES.
    fn shape_derivs(&self, xi: f64, eta: f64) -> ([f64; NODES], [f64; NODES]);

    /// Gauss integration points as (xi, eta) pairs.
    fn gauss_points(&self) -> &'static [[f64; 2]];

    /// Gauss weights (same length as `gauss_points`).
    fn gauss_weights(&self) -> &'static [f64];
}

// ============================================================================
// Internal helpers
// ============================================================================

#[inline(always)]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn dot3x3(c: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        c[0][0] * v[0] + c[0][1] * v[1] + c[0][2] * v[2],
        c[1][0] * v[0] + c[1][1] * v[1] + c[1][2] * v[2],
        c[2][0] * v[0] + c[2][1] * v[1] + c[2][2] * v[2],
    ]
}

fn symmetrize_flat(ke: &mut [f64], n: usize) {
    for i in 0..n {
        for j in i + 1..n {
            let avg = 0.5 * (ke[i * n + j] + ke[j * n + i]);
            ke[i * n + j] = avg;
            ke[j * n + i] = avg;
        }
    }
}

// ============================================================================
// Generic 2D integrators
// ============================================================================

/// Compute a plane-strain 2D stiffness matrix using the generic Gauss loop.
///
/// Accumulates into `out` (flat row-major, length `(NODES*2)²`), then symmetrises.
/// Caller must pass a zeroed slice.
pub fn integrate_ke_2d<const NODES: usize, E: ReferenceElement2D<NODES>>(
    elem: &E,
    x: &[f64; NODES],
    y: &[f64; NODES],
    c: &[[f64; 3]; 3],
    out: &mut [f64],
) {
    let n_dof = 2 * NODES;
    debug_assert_eq!(out.len(), n_dof * n_dof);

    for (gp, &w) in elem.gauss_points().iter().zip(elem.gauss_weights()) {
        let xi  = gp[0];
        let eta = gp[1];

        let (dn_xi, dn_eta) = elem.shape_derivs(xi, eta);

        let j00: f64 = (0..NODES).map(|k| dn_xi[k]  * x[k]).sum();
        let j01: f64 = (0..NODES).map(|k| dn_xi[k]  * y[k]).sum();
        let j10: f64 = (0..NODES).map(|k| dn_eta[k] * x[k]).sum();
        let j11: f64 = (0..NODES).map(|k| dn_eta[k] * y[k]).sum();

        let det_j   = j00 * j11 - j01 * j10;
        let inv_det = 1.0 / det_j;
        let ji00 =  j11 * inv_det;
        let ji01 = -j01 * inv_det;
        let ji10 = -j10 * inv_det;
        let ji11 =  j00 * inv_det;

        let mut dn_dx = [0.0f64; NODES];
        let mut dn_dy = [0.0f64; NODES];
        for k in 0..NODES {
            dn_dx[k] = ji00 * dn_xi[k] + ji01 * dn_eta[k];
            dn_dy[k] = ji10 * dn_xi[k] + ji11 * dn_eta[k];
        }

        let scale = det_j * w;

        for i in 0..NODES {
            let bi  = [dn_dx[i], 0.0,      dn_dy[i]];
            let bi1 = [0.0,      dn_dy[i], dn_dx[i]];

            for j in 0..NODES {
                let bj  = [dn_dx[j], 0.0,      dn_dy[j]];
                let bj1 = [0.0,      dn_dy[j], dn_dx[j]];

                let cbj  = dot3x3(c, &bj);
                let cbj1 = dot3x3(c, &bj1);

                out[(2 * i) * n_dof + 2 * j]             += scale * dot3(&bi,  &cbj);
                out[(2 * i) * n_dof + 2 * j + 1]         += scale * dot3(&bi,  &cbj1);
                out[(2 * i + 1) * n_dof + 2 * j]         += scale * dot3(&bi1, &cbj);
                out[(2 * i + 1) * n_dof + 2 * j + 1]     += scale * dot3(&bi1, &cbj1);
            }
        }
    }

    symmetrize_flat(out, n_dof);
}

/// Compute a 2D consistent mass matrix using the generic Gauss loop.
///
/// Caller must pass a zeroed slice of length `(NODES*2)²`.
pub fn integrate_me_2d<const NODES: usize, E: ReferenceElement2D<NODES>>(
    elem: &E,
    x: &[f64; NODES],
    y: &[f64; NODES],
    rho: f64,
    out: &mut [f64],
) {
    let n_dof = 2 * NODES;
    debug_assert_eq!(out.len(), n_dof * n_dof);

    for (gp, &w) in elem.gauss_points().iter().zip(elem.gauss_weights()) {
        let xi  = gp[0];
        let eta = gp[1];

        let n_sf            = elem.shape_functions(xi, eta);
        let (dn_xi, dn_eta) = elem.shape_derivs(xi, eta);

        let j00: f64 = (0..NODES).map(|k| dn_xi[k]  * x[k]).sum();
        let j01: f64 = (0..NODES).map(|k| dn_xi[k]  * y[k]).sum();
        let j10: f64 = (0..NODES).map(|k| dn_eta[k] * x[k]).sum();
        let j11: f64 = (0..NODES).map(|k| dn_eta[k] * y[k]).sum();

        let det_j = j00 * j11 - j01 * j10;
        let scale = rho * det_j * w;

        for p in 0..NODES {
            for q in 0..NODES {
                let v = scale * n_sf[p] * n_sf[q];
                out[(2 * p) * n_dof + 2 * q]         += v;
                out[(2 * p + 1) * n_dof + 2 * q + 1] += v;
            }
        }
    }

    symmetrize_flat(out, n_dof);
}

/// Compute a 2D body load vector using the generic Gauss loop.
///
/// Caller must pass a zeroed slice of length `NODES*2`.
pub fn integrate_fext_2d<const NODES: usize, E: ReferenceElement2D<NODES>>(
    elem: &E,
    x: &[f64; NODES],
    y: &[f64; NODES],
    b: [f64; 2],
    out: &mut [f64],
) {
    debug_assert_eq!(out.len(), 2 * NODES);

    for (gp, &w) in elem.gauss_points().iter().zip(elem.gauss_weights()) {
        let xi  = gp[0];
        let eta = gp[1];

        let n_sf            = elem.shape_functions(xi, eta);
        let (dn_xi, dn_eta) = elem.shape_derivs(xi, eta);

        let j00: f64 = (0..NODES).map(|k| dn_xi[k]  * x[k]).sum();
        let j01: f64 = (0..NODES).map(|k| dn_xi[k]  * y[k]).sum();
        let j10: f64 = (0..NODES).map(|k| dn_eta[k] * x[k]).sum();
        let j11: f64 = (0..NODES).map(|k| dn_eta[k] * y[k]).sum();

        let det_j = j00 * j11 - j01 * j10;
        let scale = det_j * w;

        for p in 0..NODES {
            out[2 * p]     += scale * n_sf[p] * b[0];
            out[2 * p + 1] += scale * n_sf[p] * b[1];
        }
    }
}

// ============================================================================
// Tests: verify generic integrators match specialized element implementations
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elements::quad::{Quad4Precomputed, Quad8Precomputed, Quad9Precomputed};
    use crate::elements::quad::{Quad4Ref, Quad8Ref, Quad9Ref};

    fn constitutive_plane_strain(e: f64, nu: f64) -> [[f64; 3]; 3] {
        let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        [
            [lambda + 2.0 * mu, lambda, 0.0],
            [lambda, lambda + 2.0 * mu, 0.0],
            [0.0, 0.0, mu],
        ]
    }

    #[test]
    fn quad4_ke_generic_matches_original() {
        let coords: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        let e = 2.0e11;
        let nu = 0.3;

        // Original
        let pre = Quad4Precomputed::new(&coords);
        let orig = pre.compute_ke_global(e, nu);

        // Generic
        let mut x = [0.0f64; 4];
        let mut y = [0.0f64; 4];
        for i in 0..4 { x[i] = coords[i][0]; y[i] = coords[i][1]; }
        let c = constitutive_plane_strain(e, nu);
        let mut ke_generic = [0.0f64; 64];
        integrate_ke_2d(&Quad4Ref, &x, &y, &c, &mut ke_generic);

        for i in 0..64 {
            let diff = (orig[i] - ke_generic[i]).abs();
            let scale = orig[i].abs().max(1.0);
            assert!(diff / scale < 1e-10, "QUAD4 Ke mismatch at [{},{}]: orig={}, generic={}", i/8, i%8, orig[i], ke_generic[i]);
        }
    }

    #[test]
    fn quad4_me_generic_matches_original() {
        let coords: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
        let rho = 7800.0;

        let pre = Quad4Precomputed::new(&coords);
        let orig = pre.compute_me_global(rho);

        let mut x = [0.0f64; 4];
        let mut y = [0.0f64; 4];
        for i in 0..4 { x[i] = coords[i][0]; y[i] = coords[i][1]; }
        let mut me_generic = [0.0f64; 64];
        integrate_me_2d(&Quad4Ref, &x, &y, rho, &mut me_generic);

        for i in 0..64 {
            let diff = (orig[i] - me_generic[i]).abs();
            assert!(diff < 1e-8, "QUAD4 Me mismatch at [{},{}]: {} vs {}", i/8, i%8, orig[i], me_generic[i]);
        }
    }
}
