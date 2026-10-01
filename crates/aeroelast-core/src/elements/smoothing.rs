//! Strain-smoothed element (SSE) membrane field for the MITC3+ triangle.
//!
//! Reference: Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite
//! element", *Computers and Structures* 223:106096, 2019.
//!
//! The membrane strain of a target triangle is smoothed with its three edge
//! neighbours.  The paper's recipe, in the order it applies it:
//!
//! 1. the covariant membrane strain is evaluated at each element centre
//!    (`r = s = 1/3`, `t = 0`);
//! 2. the neighbour's strain is carried into the target's convected
//!    coordinates, Eq. (15):
//!    `e_ij^(k) = e_ln^(k) (g_i^(e)·g^l^(k)) (g_j^(e)·g^n^(k))`;
//! 3. the target's and the neighbour's strains are averaged weighted by area,
//!    Eq. (16), where the neighbour's area is first projected onto the target's
//!    mid-surface plane, Eq. (17): `A_bar^(k) = (n^(e)·n^(k)) A^(k)`.  The
//!    smoothing therefore fades to nothing as the angle between the two
//!    elements approaches 90 degrees;
//! 4. if an edge has no neighbour, the pairwise strain falls back to the
//!    target's own strain - the boundary rule stated in the text after
//!    Eq. (17);
//! 5. the three pairwise results are assigned to the three Gauss points in a
//!    cyclic pattern, Eq. (18):
//!    `e^(A) = (e^(3) + e^(1))/2`, `e^(B) = (e^(1) + e^(2))/2`,
//!    `e^(C) = (e^(2) + e^(3))/2`.
//!
//! Bending and transverse shear are untouched: the paper keeps the original
//! `b1 eij` and `b2 eij` for the bending strains and adopts the MITC3+ assumed
//! strains for the transverse shear, so this module only ever touches the
//! membrane field.
//!
//! This module holds the geometry and linear algebra only.  Composing it with
//! the element's membrane B-matrices, and the assembly, belong to the caller.
//!
//! ## Strain convention
//!
//! Every operator here acts on the **engineering** strain vector
//! `[e11, e22, gamma12]` with `gamma12 = 2*e12`, which is what the element
//! kernels carry, in the element's local Cartesian frame.  Covariant components
//! use the same convention, so `gamma_cov = 2*e_cov12`.

use nalgebra::{Matrix2, Matrix3, Vector2, Vector3};

/// The 3x3 operator of `e -> m e m^T` on the engineering strain vector.
///
/// With `e` the symmetric tensor and `[e11, e22, 2 e12]` its vector form, the
/// transformed components are
///
/// ```text
/// (m e m^T)11 = m11^2 e11 + m12^2 e22 + m11 m12 (2 e12)
/// (m e m^T)22 = m21^2 e11 + m22^2 e22 + m21 m22 (2 e12)
/// 2 (m e m^T)12 = 2 m11 m21 e11 + 2 m12 m22 e22 + (m11 m22 + m12 m21)(2 e12)
/// ```
pub fn tensor_operator(m: &Matrix2<f64>) -> Matrix3<f64> {
    let (m11, m12) = (m[(0, 0)], m[(0, 1)]);
    let (m21, m22) = (m[(1, 0)], m[(1, 1)]);
    Matrix3::new(
        m11 * m11,
        m12 * m12,
        m11 * m12,
        m21 * m21,
        m22 * m22,
        m21 * m22,
        2.0 * m11 * m21,
        2.0 * m12 * m22,
        m11 * m22 + m12 * m21,
    )
}

/// Local Cartesian strain to covariant strain: `e_cov = J e J^T`.
///
/// `j_mat` follows the element kernel's convention, its **rows** are the
/// covariant base vectors at the element centre:
/// `j_mat = [[g_r.x, g_r.y], [g_s.x, g_s.y]]`, so `(J e J^T)_ij = g_i·e·g_j`
/// is the covariant strain.
pub fn to_covariant(strain: &Vector3<f64>, j_mat: &Matrix2<f64>) -> Vector3<f64> {
    tensor_operator(j_mat) * strain
}

/// Covariant strain to local Cartesian strain, the inverse of [`to_covariant`].
/// Returns `None` for a singular Jacobian.
pub fn from_covariant(covariant: &Vector3<f64>, j_mat: &Matrix2<f64>) -> Option<Vector3<f64>> {
    Some(tensor_operator(&j_mat.try_inverse()?) * covariant)
}

/// The operator of Eq. (15): the neighbour's covariant strain carried into the
/// target's convected coordinates.
///
/// Eq. (15) is `e~_ij = e_ln (g_i^(e)·g^l^(k))(g_j^(e)·g^n^(k))`.  The target's
/// covariant base vectors and the neighbour's contravariant base vectors live in
/// different tangent planes, so the dot products need the relative rotation
/// between the two element frames.  With `R = t3` mapping global coordinates to
/// an element's local ones, and `J` the 2x2 covariant Jacobian, the operator is
///
/// ```text
/// M = J_target · Q · J_neighbour^-1,   Q = (R_target R_neighbour^T)[0..2, 0..2]
/// ```
///
/// Dropping `Q` made the two transforms cancel and left the smoothing as an
/// unrotated average of local Cartesian strains, which over-stiffens a curved
/// shell (issue #2).  Out-of-plane strain components are neglected, as stated
/// after Eq. (15).  Returns `None` for a singular neighbour Jacobian.
pub fn convected_operator(
    j_target: &Matrix2<f64>,
    frame_target: &Matrix3<f64>,
    j_neighbour: &Matrix2<f64>,
    frame_neighbour: &Matrix3<f64>,
) -> Option<Matrix3<f64>> {
    let relative = frame_target * frame_neighbour.transpose();
    let q = Matrix2::new(
        relative[(0, 0)],
        relative[(0, 1)],
        relative[(1, 0)],
        relative[(1, 1)],
    );
    Some(tensor_operator(&(j_target * q * j_neighbour.try_inverse()?)))
}

/// Eq. (17): the neighbour's area projected onto the target's mid-surface plane.
///
/// `n_target` and `n_neighbour` are unit centre normals, so the factor is
/// `cos(theta)` and the projected area vanishes at 90 degrees.
pub fn projected_area(
    n_target: &Vector3<f64>,
    n_neighbour: &Vector3<f64>,
    area_neighbour: f64,
) -> f64 {
    n_target.dot(n_neighbour) * area_neighbour
}

/// Eq. (16): the area-weighted pairwise smoothed covariant strain.
///
/// Both strains must already be covariant and expressed in the target's
/// convected coordinates; `area_neighbour` must be the projected area from
/// [`projected_area`].
pub fn pairwise_smoothed(
    target: &Vector3<f64>,
    neighbour: &Vector3<f64>,
    area_target: f64,
    area_neighbour: f64,
) -> Vector3<f64> {
    let denominator = area_target + area_neighbour;
    if denominator.abs() < f64::MIN_POSITIVE {
        return *target;
    }
    (target * area_target + neighbour * area_neighbour) / denominator
}

/// Eq. (18): assign the three pairwise strains to the three Gauss points.
///
/// The pairing is cyclic - `A` pairs edges 3 and 1, `B` pairs 1 and 2, `C` pairs
/// 2 and 3 - so every edge appears in exactly two of the three assigned strains.
pub fn assign_to_gauss_points(pairwise: &[Vector3<f64>; 3]) -> [Vector3<f64>; 3] {
    let [e1, e2, e3] = pairwise;
    [(e3 + e1) * 0.5, (e1 + e2) * 0.5, (e2 + e3) * 0.5]
}

/// The three edge neighbours of each element, keyed by local edge.
///
/// Edge `k` of a triangle is the edge from local node `k` to local node
/// `(k + 1) % 3`, so the edges follow the element's own node ordering.  Two
/// elements are neighbours when they share exactly that pair of nodes.  `None`
/// means the edge is on a boundary (or the edge is shared by more than two
/// elements, which is not a manifold mesh and is left unsmoothed).
pub fn edge_neighbours(connectivity: &[Vec<usize>]) -> Vec<[Option<usize>; 3]> {
    use std::collections::HashMap;

    let mut edge_map: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
    for (element, nodes) in connectivity.iter().enumerate() {
        if nodes.len() != 3 {
            continue;
        }
        for k in 0..3 {
            let a = nodes[k];
            let b = nodes[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_map.entry(key).or_default().push((element, k));
        }
    }

    let mut neighbours = vec![[None; 3]; connectivity.len()];
    for entries in edge_map.values() {
        if entries.len() != 2 {
            continue;
        }
        let (element_a, edge_a) = entries[0];
        let (element_b, edge_b) = entries[1];
        neighbours[element_a][edge_a] = Some(element_b);
        neighbours[element_b][edge_b] = Some(element_a);
    }
    neighbours
}

/// Covariant geometry of one triangle, as the smoothing needs it.
#[derive(Debug, Clone, Copy)]
pub struct TriangleFrame {
    /// Rows are the covariant base vectors `g_r`, `g_s` at the element centre.
    pub j_mat: Matrix2<f64>,
    /// Mid-surface area.
    pub area: f64,
    /// Unit normal at the element centre, in global coordinates.
    pub normal: Vector3<f64>,
    /// The element's local-to-global rotation `t3` (rows are the local basis
    /// vectors in global coordinates).  Eq. (15) needs it to relate the target's
    /// tangent plane to the neighbour's.
    pub frame: Matrix3<f64>,
}

/// The smoothed membrane strain operator of one triangle.
///
/// Each Gauss point's covariant strain is a weighted sum of four covariant
/// strains: the target's own and the three neighbours', the latter already
/// carried into the target's convected coordinates by Eq. (15).  Entry 0 is the
/// target; entry `1 + k` is the element across local edge `k`.
#[derive(Debug, Clone)]
pub struct SmoothedMembraneStrain {
    /// Per Gauss point, the 3x3 weights on `[target, edge 0, edge 1, edge 2]`.
    pub weights: [[Matrix3<f64>; 4]; 3],
    /// Element index of each entry; `None` for a boundary edge.
    pub sources: [Option<usize>; 4],
}

/// Build the smoothed membrane strain operator of every triangle.
///
/// The weights are purely geometric, so they are built once; the caller applies
/// them to whatever covariant membrane strains it has.
///
/// A degenerate frame (non-invertible Jacobian) gets the identity smoothing -
/// weight 1 on its own strain - so a broken element cannot poison its neighbours.
pub fn smoothed_membrane_strain(
    frames: &[TriangleFrame],
    neighbours: &[[Option<usize>; 3]],
) -> Vec<SmoothedMembraneStrain> {
    let mut out = Vec::with_capacity(frames.len());

    for (element, frame) in frames.iter().enumerate() {
        let mut pairwise_weights = [[Matrix3::zeros(); 4]; 3];
        let mut sources = [Some(element), None, None, None];

        for edge in 0..3 {
            let neighbour = neighbours
                .get(element)
                .and_then(|n| n[edge])
                .filter(|k| *k < frames.len());

            let Some(neighbour) = neighbour else {
                // Boundary rule: the pairwise strain is the element's own.
                pairwise_weights[edge][0] = Matrix3::identity();
                continue;
            };

            let Some(convected) = convected_operator(
                &frame.j_mat,
                &frame.frame,
                &frames[neighbour].j_mat,
                &frames[neighbour].frame,
            ) else {
                // Degenerate neighbour: fall back to the element's own strain.
                pairwise_weights[edge][0] = Matrix3::identity();
                continue;
            };

            let projected = projected_area(
                &frame.normal,
                &frames[neighbour].normal,
                frames[neighbour].area,
            );
            let denominator = frame.area + projected;
            if denominator.abs() < f64::MIN_POSITIVE {
                pairwise_weights[edge][0] = Matrix3::identity();
                continue;
            }

            pairwise_weights[edge][0] = Matrix3::identity() * (frame.area / denominator);
            pairwise_weights[edge][1 + edge] = convected * (projected / denominator);
            sources[1 + edge] = Some(neighbour);
        }

        // Eq. (18): A pairs edges 3 and 1, B pairs 1 and 2, C pairs 2 and 3.
        const PAIRS: [[usize; 2]; 3] = [[2, 0], [0, 1], [1, 2]];
        let mut weights = [[Matrix3::zeros(); 4]; 3];
        for (gauss_point, pair) in PAIRS.iter().enumerate() {
            for edge in pair {
                for entry in 0..4 {
                    weights[gauss_point][entry] += pairwise_weights[*edge][entry] * 0.5;
                }
            }
        }

        // Drop sources that ended up with no weight at all.
        for entry in 1..4 {
            if weights.iter().all(|w| w[entry].norm() < f64::MIN_POSITIVE) {
                sources[entry] = None;
            }
        }

        out.push(SmoothedMembraneStrain { weights, sources });
    }

    out
}

/// Apply a smoothed operator to a set of covariant strains.
pub fn apply_smoothed(
    operator: &SmoothedMembraneStrain,
    covariant_strains: &[Vector3<f64>],
) -> [Vector3<f64>; 3] {
    let mut out = [Vector3::zeros(); 3];
    for (gauss_point, weights) in operator.weights.iter().enumerate() {
        for (entry, weight) in weights.iter().enumerate() {
            if let Some(source) = operator.sources[entry] {
                out[gauss_point] += weight * covariant_strains[source];
            }
        }
    }
    out
}

/// The covariant base vectors of a frame.
pub fn base_vectors(j_mat: &Matrix2<f64>) -> (Vector2<f64>, Vector2<f64>) {
    (
        Vector2::new(j_mat[(0, 0)], j_mat[(0, 1)]),
        Vector2::new(j_mat[(1, 0)], j_mat[(1, 1)]),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Matrix2;

    /// A right triangle (0,0), (1,0), (0,1): J = identity, area 1/2.
    fn unit_triangle() -> TriangleFrame {
        TriangleFrame {
            j_mat: Matrix2::identity(),
            area: 0.5,
            normal: Vector3::new(0.0, 0.0, 1.0),
            frame: Matrix3::identity(),
        }
    }

    fn strain(e11: f64, e22: f64, gamma12: f64) -> Vector3<f64> {
        Vector3::new(e11, e22, gamma12)
    }

    #[test]
    fn tensor_operator_reproduces_the_metric_transform() {
        // A strain along x, in a frame rotated 45 degrees: the covariant
        // components must come out of e_cov = J e J^T.
        let j = Matrix2::new(1.0, 0.0, 1.0, 1.0);
        let e = strain(1.0, 0.0, 0.0);
        let covariant = to_covariant(&e, &j);
        // g_r = (1,0), g_s = (1,1): e_rr = 1, e_ss = 1, e_rs = 1 -> gamma = 2.
        assert!((covariant[0] - 1.0).abs() < 1e-14);
        assert!((covariant[1] - 1.0).abs() < 1e-14);
        assert!((covariant[2] - 2.0).abs() < 1e-14);

        let back = from_covariant(&covariant, &j).expect("invertible");
        assert!((back - e).norm() < 1e-14);
    }

    #[test]
    fn convected_operator_is_identity_for_equal_frames() {
        let j = Matrix2::new(0.3, 0.1, -0.2, 0.4);
        let identity = Matrix3::identity();
        let operator = convected_operator(&j, &identity, &j, &identity).expect("invertible");
        assert!((operator - Matrix3::identity()).norm() < 1e-14);
    }

    #[test]
    fn boundary_rule_keeps_the_elements_own_strain() {
        let frames = vec![unit_triangle()];
        let neighbours = vec![[None, None, None]];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let own = strain(1.0e-3, -2.0e-3, 5.0e-4);
        let assigned = apply_smoothed(&operators[0], &[own]);
        for value in assigned {
            assert!(
                (value - own).norm() < 1e-15,
                "a boundary triangle must reproduce its own strain, got {value:?}"
            );
        }
        assert_eq!(operators[0].sources, [Some(0), None, None, None]);
    }

    #[test]
    fn each_gauss_point_weight_sums_to_one() {
        // Two triangles sharing an edge, so one neighbour exists and two edges are
        // on the boundary: every Gauss point must still be a partition of unity.
        let frames = vec![unit_triangle(), unit_triangle()];
        let neighbours = vec![[Some(1), None, None], [Some(0), None, None]];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        for operator in &operators {
            for weights in &operator.weights {
                let sum: Matrix3<f64> = weights.iter().sum();
                assert!(
                    (sum - Matrix3::identity()).norm() < 1e-14,
                    "weights must sum to the identity, got {sum:?}"
                );
            }
        }
    }

    #[test]
    fn a_constant_strain_field_is_reproduced_exactly() {
        // The patch test: if every element carries the same covariant strain, the
        // smoothing cannot change it, whatever the weights are.
        let frames = vec![unit_triangle(), unit_triangle(), unit_triangle()];
        let neighbours = vec![
            [Some(1), Some(2), None],
            [Some(0), None, Some(2)],
            [None, Some(0), Some(1)],
        ];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let constant = strain(2.0e-4, -3.0e-4, 1.0e-4);
        let all = vec![constant, constant, constant];
        for operator in &operators {
            for value in apply_smoothed(operator, &all) {
                assert!(
                    (value - constant).norm() < 1e-15,
                    "constant strain must survive smoothing, got {value:?}"
                );
            }
        }
    }

    #[test]
    fn equal_areas_average_the_two_elements_equally() {
        let frames = vec![unit_triangle(), unit_triangle()];
        let neighbours = vec![[Some(1), None, None], [Some(0), None, None]];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let own = strain(1.0, 0.0, 0.0);
        let other = strain(0.0, 1.0, 0.0);
        let assigned = apply_smoothed(&operators[0], &[own, other]);

        // Only edge 0 has a neighbour, so pairwise_0 = (own + other)/2 while
        // pairwise_1 and pairwise_2 fall back to the element's own strain
        // (the boundary rule).  Eq. (18) then gives
        //   A = (pairwise_2 + pairwise_0)/2 = 3/4 own + 1/4 other
        //   B = (pairwise_0 + pairwise_1)/2 = 3/4 own + 1/4 other
        //   C = (pairwise_1 + pairwise_2)/2 = own
        let smoothed = own * 0.75 + other * 0.25;
        assert!(
            (assigned[0] - smoothed).norm() < 1e-15,
            "got {:?}",
            assigned[0]
        );
        assert!(
            (assigned[1] - smoothed).norm() < 1e-15,
            "got {:?}",
            assigned[1]
        );
        assert!((assigned[2] - own).norm() < 1e-15, "got {:?}", assigned[2]);
    }

    #[test]
    fn all_three_neighbours_exercise_the_cyclic_pairing() {
        // A central triangle with a neighbour across every edge, all flat and of
        // equal area, so every pairwise strain is the mean of the two elements and
        // Eq. (18) is exercised without any boundary fallback.
        let frames = vec![unit_triangle(); 4];
        let neighbours = vec![
            [Some(1), Some(2), Some(3)],
            [Some(0), None, None],
            [Some(0), None, None],
            [Some(0), None, None],
        ];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let own = strain(1.0, 0.0, 0.0);
        let e1 = strain(0.0, 1.0, 0.0);
        let e2 = strain(0.0, 0.0, 1.0);
        let e3 = strain(1.0, 1.0, 1.0);
        let assigned = apply_smoothed(&operators[0], &[own, e1, e2, e3]);

        // pairwise_k = (own + e_k)/2, so
        //   A = (pairwise_3 + pairwise_1)/2 = 1/2 own + 1/4 e3 + 1/4 e1
        //   B = (pairwise_1 + pairwise_2)/2 = 1/2 own + 1/4 e1 + 1/4 e2
        //   C = (pairwise_2 + pairwise_3)/2 = 1/2 own + 1/4 e2 + 1/4 e3
        let expected = [
            own * 0.5 + e3 * 0.25 + e1 * 0.25,
            own * 0.5 + e1 * 0.25 + e2 * 0.25,
            own * 0.5 + e2 * 0.25 + e3 * 0.25,
        ];
        for (gauss_point, (got, want)) in assigned.iter().zip(expected.iter()).enumerate() {
            assert!(
                (got - want).norm() < 1e-15,
                "Gauss point {gauss_point}: got {got:?}, expected {want:?}"
            );
        }
        // Each edge appears in exactly two Gauss points, so the total is preserved.
        let total: Vector3<f64> = assigned.iter().sum();
        let expected_total: Vector3<f64> = expected.iter().sum();
        assert!((total - expected_total).norm() < 1e-15);
    }

    #[test]
    fn the_projected_area_switches_the_smoothing_off_at_ninety_degrees() {
        let mut neighbour = unit_triangle();
        neighbour.normal = Vector3::new(1.0, 0.0, 0.0); // perpendicular to the target
        // A frame whose third row is that normal: rotate 90 degrees about y.
        neighbour.frame = Matrix3::new(0.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0);
        let frames = vec![unit_triangle(), neighbour];
        let neighbours = vec![[Some(1), None, None], [None, None, None]];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let own = strain(1.0, 0.0, 0.0);
        let other = strain(0.0, 1.0, 0.0);
        // The neighbour contributes zero, so the pairwise strain is the element's
        // own and every Gauss point reproduces it.
        for value in apply_smoothed(&operators[0], &[own, other]) {
            assert!((value - own).norm() < 1e-15, "got {value:?}");
        }
    }

    #[test]
    fn unequal_areas_weight_the_average_by_projected_area() {
        let mut neighbour = unit_triangle();
        neighbour.area = 1.5; // three times the target
        let frames = vec![unit_triangle(), neighbour];
        let neighbours = vec![[Some(1), None, None], [None, None, None]];
        let operators = smoothed_membrane_strain(&frames, &neighbours);

        let own = strain(1.0, 0.0, 0.0);
        let other = strain(0.0, 1.0, 0.0);
        let assigned = apply_smoothed(&operators[0], &[own, other]);

        // pairwise_0 = (0.5 own + 1.5 other) / 2 = 1/4 own + 3/4 other, and edges 1
        // and 2 are on the boundary so pairwise_1 = pairwise_2 = own.  Hence
        //   B = (pairwise_0 + pairwise_1)/2 = 1/2 (1/4 own + 3/4 other + own)
        //     = 5/8 own + 3/8 other
        let expected = own * 0.625 + other * 0.375;
        assert!(
            (assigned[1] - expected).norm() < 1e-15,
            "got {:?}",
            assigned[1]
        );
        // C pairs the two boundary edges, so it is the element's own strain.
        assert!((assigned[2] - own).norm() < 1e-15, "got {:?}", assigned[2]);
    }

    #[test]
    fn assign_to_gauss_points_is_cyclic() {
        let pairwise = [
            strain(1.0, 0.0, 0.0),
            strain(0.0, 1.0, 0.0),
            strain(0.0, 0.0, 1.0),
        ];
        let assigned = assign_to_gauss_points(&pairwise);
        assert!((assigned[0] - (pairwise[2] + pairwise[0]) * 0.5).norm() < 1e-15);
        assert!((assigned[1] - (pairwise[0] + pairwise[1]) * 0.5).norm() < 1e-15);
        assert!((assigned[2] - (pairwise[1] + pairwise[2]) * 0.5).norm() < 1e-15);
        // Every edge appears exactly twice across the three Gauss points.
        let total = assigned[0] + assigned[1] + assigned[2];
        let expected_total = pairwise[0] + pairwise[1] + pairwise[2];
        assert!((total - expected_total).norm() < 1e-15);
    }

    #[test]
    fn edge_neighbours_finds_the_shared_edge() {
        // Two triangles sharing the node pair (0,1).  That pair is local edge 0 of
        // BOTH elements: edge k runs from local node k to local node (k+1) % 3, so
        // element 0 has edges (0,1), (1,2), (2,0) and element 1 has (1,0), (0,3),
        // (3,1).
        let connectivity = vec![vec![0, 1, 2], vec![1, 0, 3]];
        let neighbours = edge_neighbours(&connectivity);
        assert_eq!(neighbours[0][0], Some(1));
        assert_eq!(neighbours[0][1], None);
        assert_eq!(neighbours[0][2], None);
        assert_eq!(neighbours[1][0], Some(0));
        assert_eq!(neighbours[1][1], None);
        assert_eq!(neighbours[1][2], None);
    }

    #[test]
    fn edge_neighbours_leaves_a_non_manifold_edge_unsmoothed() {
        // Three elements sharing one edge: not a manifold, so no neighbour.
        let connectivity = vec![vec![0, 1, 2], vec![1, 0, 3], vec![0, 1, 4]];
        let neighbours = edge_neighbours(&connectivity);
        for element in 0..3 {
            for edge in 0..3 {
                assert_eq!(
                    neighbours[element][edge], None,
                    "a non-manifold edge must not be smoothed"
                );
            }
        }
    }

    #[test]
    fn edge_neighbours_ignores_non_triangles() {
        let connectivity = vec![vec![0, 1, 2, 3], vec![0, 1, 2]];
        let neighbours = edge_neighbours(&connectivity);
        assert_eq!(neighbours[0], [None, None, None]);
        assert_eq!(neighbours[1], [None, None, None]);
    }
}
