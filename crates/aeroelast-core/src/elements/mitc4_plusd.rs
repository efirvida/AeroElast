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

#[cfg(test)]
mod tests {
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
}
