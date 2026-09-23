# Explore — `mitc4plus-2017-faithful`

Read-only exploration for the proposal that scopes a faithful MITC4+ (Ko, Lee &
Bathe 2017, paper A) as a **new** element, 5 DOF/node, separate from the current
`mitc4.rs` 6-DOF flat-projection element.

Scope of this document: inventory of the theoretical tests (two tiers), the
integration surface the new element must fit, and the gap/risks. **No code was
written or modified.** Evidence is a static read of the tree; the tool set for
this phase had **no shell and no PDF tooling**, so paper A was not re-read (the
transcription in `docs/formulations/mitc4plus-2017-extract.md` was used) and no
suite was run (the session baselines were taken as given).

Sources used:

- `docs/formulations/mitc4plus-2017-extract.md` (paper A transcription, with page
  citations).
- `odd/tasks/mitc4plus-faithful.md` (plan sketch).
- `docs/formulations/shell-elements.md` (current element, code↔equation).
- `docs/validation-matrix.md` (test-by-test references; snapshot of commit
  `b2c62ff`, **stale** relative to the session baseline — see §5).
- Source read directly: `crates/aeroelast-core/src/elements/mitc4.rs`,
  `.../elements/mod.rs`, `crates/aeroelast-core/src/assembly/{assembler.rs,topology.rs}`,
  `crates/aeroelast-py/src/{elements.rs,lib.rs,assembler.rs}`,
  `src/aeroelast/core/assembler.py`, `src/aeroelast/elements/__init__.py`, and the
  test modules named below.

---

## 0. What is being scoped

Paper A defines a **5-DOF/node** element (3 translations + 2 director rotations
`α_i`, `β_i`; **no drilling DOF**), 3D continuum kinematics with the director
`V_n^i` and nodal thickness `a_i`, only the membrane term modified (Eqs. 17–27),
MITC4 shear (Dvorkin & Bathe 1984), full 2×2×2 Gauss, **no numerical factor**
(`docs/formulations/mitc4plus-2017-extract.md`, re-read section; §"Assumed
membrane strain field"; §"Transverse shear"). The current `mitc4.rs` is a
**different element**: 6 DOF/node, flat projection, a 2-DOF rotation bubble,
selective reduced integration, the Choi–Paik starting membrane field, a shear
correction factor, a Winkler & Plakomytis ERC drilling treatment and a `β_w`
warping penalty (`odd/tasks/mitc4plus-faithful.md`, deviations table). The
proposal is for a new element, not a patch.

### Tying-point / figure conventions I rely on

- **Fig. 3(d) (paper A).** The red arrows are the *edge vectors* `x_e^k` of
  Eq. (13) — `x_e^1=(x_2−x_1)/2`, etc. They are **not** the thickness `a_i` and
  **not** the director `V_n^i`. Extract, §"Eq. (13) — the edge vectors", and the
  warning in §"Element geometry and kinematics" (`a_i`/`V_n^i` are stated in words
  after Eq. 1). *Do not confuse them when reading the figure.*
- **Fig. 4 (paper A).** Five membrane tying points: `A(0,+1)`, `B(0,−1)`,
  `C(+1,0)`, `D(−1,0)`, `E(0,0)`, sampling `e_rr(A)`, `e_rr(B)`, `e_ss(C)`,
  `e_ss(D)`, `e_rs(E)`. Extract, §"Assumed membrane strain field (Section 3.2,
  p. 408)". The extract notes the coordinates are the code's parameterisation;
  the paper presents them graphically.
- **Fig. 2(a) of paper B**, shared MITC4 transverse-shear tying points: `A`
  (top, `s=+1`), `B` (bottom, `s=−1`), `C` (right, `r=+1`), `D` (left, `r=−1`),
  giving `ẽ_rt=½(1+s)e_rt^A+½(1−s)e_rt^B`, `ẽ_st=½(1+r)e_st^C+½(1−r)e_st^D`.
  Extract, §"Transverse shear".

---

## 1. Tier 1 — the paper's own basic tests

Paper A's Section 4 is titled **"Basic numerical tests"** and states, verbatim in
the extract (p. 410): *"the element formulation does not include any numerical
factor, and consider next the **isotropy, zero energy mode and patch tests**."*
That sentence is the **entire** Section-4 content transcribed in the extract —
the quantitative tables of Section 4 were **not** transcribed (gap, §5). The
plan sketch `odd/tasks/mitc4plus-faithful.md` F2 defines the intended checks.
What follows is therefore "what the paper requires + what the check would be",
with the exact published numbers flagged as unread.

| # | Paper A test | What the paper requires | The check the new element must pass |
| --- | --- | --- | --- |
| T1.1 | **Isotropy** | The element stiffness must not depend on the element's spatial orientation (rotating the same element and its DOF must give the same physical stiffness, i.e. the strain energy is invariant under rigid rotation of the element). | Build one element in several global orientations, push each through the local-frame transform, and assert the eigenvalues (or strain energy `uᵀKu` for a physically co-rotated `u`) are invariant to round-off. |
| T1.2 | **Zero-energy modes** | The element stiffness has **exactly the six rigid-body modes** (3 translations + 3 rotations) and **no spurious zero modes**, on flat *and* warped geometry. | Symmetrise `K`, count eigenvalues with `|λ| < tol·λ_max`; assert `count == 6`, and assert `‖K·u_rb‖ = 0` for the six physical rigid-body fields. Note the 5-DOF element represents the drilling rigid rotation through the in-plane translations with `α=β=0` (since `ω×V_n = 0`), so the six RB modes remain expressible — the count must still be 6, but the *mode construction* differs from the 6-DOF test. |
| T1.3 | **Patch tests** | Constant-strain membrane and constant-curvature bending states must be reproduced exactly at every quadrature point, including on **distorted** geometry; the assumed membrane field's Eq. (22) condition exists precisely to satisfy the patch test. | Assemble a patch of distorted elements with a linear in-plane field (`u=ax+by, v=cx+dy` → `ε=[a,d,b+c]`) and a constant-curvature field, and assert the sampled strain equals the exact constant at every Gauss point. Add the paper's own patch-test condition of Eq. (22) as a unit test (`ẽ_rs|bil = e_rs|bil` for flat geometry). |

The paper's **no-numerical-factor** statement is a constraint on all three: no
SRI, no shear correction factor, no penalty. `a_E` is **positive**
(`+2 c_r c_s / d`) per the extract's sign correction.

> **Acceptance gate (plan sketch F1–F3):** T1.1–T1.3 are the gate *before* any
> benchmark number is consulted. Only after these pass is Tier 2 read.

---

## 2. Tier 2 — the repository's existing theoretical tests

Session baseline (given, not re-run): Rust `cargo test -p aeroelast-core`
**120 passed / 0 failed** (workspace root is `crates/`); Python
`python -m pytest -m "not slow" -q` **345 passed / 2 failed / 2 skipped**. The
two failures are pre-existing: `test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple`
and `test_shell_convergence.py::test_in_bending_convergence`. The two skips are
not element-related (CalculiX/BEM/blade per the validation matrix). No active
`xfail` marker exists anywhere (`_TWISTED_BEAM_CASES` has `xfail_reason=None`
for all four).

Below: each group, exact node ids, file, reference, and current status.

### 2A. Patch and rigid-body invariants (Rust, in `crates/aeroelast-core/src/elements/mitc4.rs`)

Module `tests` starts at `mitc4.rs:2188`. These are the invariants the new
element must preserve (a 20-DOF equivalent, see §2H).

| Test (node id) | What it checks | Reference | Status |
| --- | --- | --- | --- |
| `tests::test_ke_global_leaves_all_six_rigid_body_modes_free` | exact null space for the six physical RB fields | physical RB field `u=t+ω×(x−x_c)`, `θ=ω` | pass |
| `tests::test_ke_global_has_exactly_six_zero_modes` | the null space is *exactly* 6 — no hourglass, no spurious drill-shear mode | eigenvalue guard | pass |
| `tests::test_membrane_patch_reproduces_constant_strain_at_every_gauss_point` | constant-strain patch, `u=ax+by,v=cx+dy` | analytical constant strain | pass |
| `tests::test_bending_patch_reproduces_constant_curvature_at_every_gauss_point` | constant-curvature patch, `w=½kxx x²` | analytical constant curvature | pass |
| `tests::test_b_erc_membrane_rows_annihilate_rigid_body_modes` | MITC4+ assumed membrane rows annihilate RB on flat *and* warped | physical RB field | pass |
| `tests::test_b_erc_drill_row_is_rigid_body_invariant_on_flat_geometry` | ERC drill row RB-invariant on flat | physical RB field | pass |
| `tests::test_b_erc_drill_row_warping_residual_is_characterized` | **pins a known defect**: warped drill-row RB residual in `[1e-4,1e-2]` | measured `~1.1e-3` | pass (it pins the defect) |

### 2B. Eigenvalue / zero-mode guards

| Test | File | Reference | Status |
| --- | --- | --- | --- |
| `tests::test_ke_global_is_symmetric` | `mitc4.rs` | symmetry | pass |
| `tests::test_ke_global_is_positive_semidefinite` | `mitc4.rs` | PSD | pass |
| `tests::test_ke_local_eigenvalues_nonsymmetric` | `mitc4.rs` | PSD with `cb_coupling` | pass |
| `tests::test_ke_local_flat_plate_parity` | `mitc4.rs` | symmetry/nonzero | pass |
| `test_rust_composite.py::TestMITC4BatchSanity::test_ke_shape_and_symmetry` | Python | 24×24 symmetry | pass |
| `...::TestMITC4BatchSanity::test_ke_positive_semidefinite` | Python | PSD | pass |
| `...::TestMITC4BatchSanity::test_me_shape_and_symmetry`, `test_me_positive_semidefinite` | Python | mass symmetry/PSD | pass |
| `...::TestBatchComposite::test_batch_ke_mitc4_multiple` | Python | 3-element batch symmetry + PSD incl. one warped quad | **FAILS (pre-existing)** |
| `test_material_suite.py::TestStiffnessProperties::test_k_positive_definite` / `test_k_symmetric` | Python | assembled K PSD/symmetry | pass |
| `test_rust_assembler.py::TestTangentStiffness::test_kt_at_zero_equals_k`, `test_kt_symmetric` | Python | `K_T(0)=K` | pass |
| `test_rust_assembler.py::TestInternalForces::test_fint_zero_at_zero`, `test_fint_linear_equals_ku` | Python | linear identity | pass |
| `test_rust_assembler.py::TestNewtonRaphsonConsistency::test_nr_consistency` | Python | `K_T(u)δu ≈ Δf_int` | pass |
| `tests::test_kt_fint_directional_derivative`, `..._rotations`, `..._with_drill_dofs` | `mitc4.rs` | `K_T`↔`f_int` consistency | pass |
| `tests::test_fint_linear_nonlinear_parity`, `test_kt_zero_matches_ke` | `mitc4.rs` | linear/nonlinear parity | pass |

### 2C. Analytical shell validation

| Test | File | Reference | Status |
| --- | --- | --- | --- |
| `TestIsotropicAnalytical::test_mitc4_out_of_plane`, `test_mitc4_in_plane_lateral`, `test_mitc4_axial_stiffness` | `test_material_suite.py:215` | Euler–Bernoulli with CLT `A`/`D` | pass |
| `TestOrthotropicSinglePly` (6), `TestSymmetricLaminates` (7), `TestAsymmetricLaminates` (6) | `test_material_suite.py` | E–B + Reddy CLT B-coupling | pass |
| `TestCantileverBeam::test_tip_deflection[4x2/8x4/16x8]`, `TestSimplySupportedBeam::test_center_deflection`, `TestBeamBending::test_bending_convergence`, `TestMembraneStretching::test_axial_extension`, `TestShearLocking::test_thin_plate_convergence` | `test_shell_analytical_validation.py` | `P L³/(3EI)`, `P L/(EA)`, etc. (matrix flags name-vs-body defects) | pass |
| `test_shell_comprehensive.py` (`fx/fy/fz`, ratio, modal) | Python | beam theory + Timoshenko | pass |
| `test_shell_validation_fixed.py` (`fx/fy/fz/ratio/axial convergence/modal`) | Python | beam theory | pass |
| `test_large_rotation_benchmarks.py` (`test_linear_tip_deflection_euler_bernoulli`, `test_cantilever_large_rotation_half_circle`, `test_equilibrium_path`, `test_simo_vu_quoc_rollup_360`) | Python, **MITC4 element type code 4** | elastica analytical + Simo & Vu-Quoc 1986 (weak: `REFERENCE_TABLE` is the rounded own formula) | pass |
| `test_rust_modal.py::TestSimplySupportedPlate::test_analytical_convergence` | Python | Kirchhoff `f11` | pass |
| `test_rust_modal.py` MITC4 frequency/match groups | Python | SLEPc + analytical | pass |
| `test_mass_matrix_validation.py` exact coefficients (`TestExactConsistentMassCoefficients`, `TestConsistentMassTotal`, `TestElementMassVsTotalMass`) | Python | exact `∫ρh N_iN_j`, `ρhA`, `ρh³A/12` | pass |
| CCX parity: `test_isotropic_shell_parity.py::TestIsotropicShellParity::test_transverse_tip_displacement`, `test_beam_shell_4cases_parity.py` (9), `test_composite_beam_parity.py` (5), `test_orthotropic_shell_parity.py` (3) | Python | **CalculiX 2.23 S4 / S8R** + `*SHELL SECTION, COMPOSITE` | pass (when `ccx` on PATH; otherwise skip) |

### 2D. Convergence order

| Test | File | Reference | Status |
| --- | --- | --- | --- |
| `test_shell_convergence.py::test_in_plane_bending_convergence` | Python | analytical `P L³/(3EI)`, Richardson self-convergence, `MIN_ORDER=1.5`, `EXTRAPOLATED_TOL=0.02` | **FAILS (pre-existing)** |
| `test_shell_convergence.py::test_composite_laminate_gap_mesh_study` | Python | CalculiX S8R across 4 meshes; **no value asserted** | pass (cannot fail) |

### 2E. Published benchmark columns — Ko, Lee, Lee & Bathe 2017 (paper **B**)

`tests/test_ko2017_performance.py`, 31 cases, `rtol=0.05` (3.5 per-case `tol`).
Reference: `.sources/papers/1-s2.0-S0045794917309550-main.pdf`, Tables 2–19,
MITC4 column at the mesh the test builds. All currently pass (none is in the
2-failure list).

| Test function | Tables | Benchmark | Status |
| --- | --- | --- | --- |
| `test_3_1_square_plate_tables_2_to_5` | 2–5 | clamped/SS square plate, regular+distorted | pass |
| `test_3_2_circular_plate_tables_6_to_7` | 6–7 | clamped/SS circular plate | pass (matrix flags SS rows use Table 6 cells) |
| `test_3_3_pinched_cylinder_tables_8_to_9` | 8–9 | pinched cylinder, regular+distorted | pass (matrix flags distorted cell mismatch) |
| `test_3_4_scordelis_lo_tables_10_to_11` | 10–11 | Scordelis-Lo roof, regular+distorted | pass |
| `test_3_5_twisted_beam_tables_12_to_13` | 12–13 | MacNeal–Harder twisted beam, 4 cases, N=16 | pass (input mesh `n_width=16, n_length=96`) |
| `test_3_6_hook_table_14_minimal_fix` | 14 | Raasch hook | pass |
| `test_3_7_hemisphere_cutout_tables_15_to_16` | 15–16 | hemisphere w/ cutout | pass (matrix flags two cell mismatches) |
| `test_3_8_full_hemisphere_table_17` | 17 | full hemisphere | pass |
| `test_3_9_hyperbolic_paraboloid_tables_18_to_19` | 18–19 | hyperbolic paraboloid | pass |

Documented as *targets only* (not correctness evidence): the published MITC4+ N=16
cells for Scordelis-Lo regular `0.9973`, distorted `0.9942` (`validation-matrix.md`
§10.2). The Scordelis-Lo/pinched/square-plate rows above are the "Scordelis-Lo,
pinched shell, square plates" group the request asks for.

### 2F. Composite / orthotropic parity

Covered by §2C: `test_material_suite.py` laminate groups, `test_composite_b_coupling.py`
(reference mostly commented out — sign + floor only), `test_rust_composite.py`
(invariants only), `test_orthotropic_shell_parity.py`, `test_composite_beam_parity.py`.
All pass. Composite uses element codes `33`/`44` → `Mitc4Composite` (`assembler.py:347`).

### 2G. Rust unit tests **about the current 6-DOF architecture** that need a 5-DOF equivalent

Everything that indexes with `6*i`, `Vec24`, `Mat24`, `build_t24`, or the
`θ_z`/drill DOF `6*i+5` is written against the 6-DOF layout. For a 5-DOF
(20-DOF) element, the *invariant* must be re-expressed; the following need an
equivalent:

- `test_ke_global_leaves_all_six_rigid_body_modes_free` — `rigid_body_mode()` writes
  `θ=ω` into DOF slots `3..6`; the 5-DOF version must map `ω` onto `(α,β)` in the
  per-node `(V_1,V_2)` basis.
- `test_ke_global_has_exactly_six_zero_modes` — `Mat24`/`compute_ke_global` → 20×20.
- `test_membrane_patch...`, `test_bending_patch...` — `b_m_mitc4_plus`/`b_kappa` are
  3×24/3×24; need 5-DOF operators.
- `test_ke_global_is_symmetric`, `test_ke_global_is_positive_semidefinite`.
- `test_fint_linear_nonlinear_parity`, `test_kt_zero_matches_ke`,
  `test_kt_fint_directional_derivative`, `..._rotations` — `Vec24`, `build_t24`.
- All `me_global*` mass tests (`Mat24`, `6*i` strides) — a 5-DOF mass matrix has no
  drilling rotary-inertia block and is 20×20.
- `test_body_load_global_*`, `test_k_sigma_global_*`,
  `test_centrifugal_prestress_*`, `compute_element_stress` tests — assumed `6*i`.

**Do NOT port — delete with the old element** (they are about the 6-DOF
*deviation*, not the paper): `test_drill_midside_derivatives_match_paper_eq11`,
`test_drill_membrane_operator_*` (5 tests), `test_b_erc_*` (5 tests),
`test_erc_covariant_to_local_*` (2 tests), `test_ke_local_erc_flat_rigid_body_and_symmetry`,
`test_enhanced_drill_stiffness`, `test_drill_warping_moment_flat_element`,
`test_ke_local_eigenvalues_nonsymmetric` (uses `cb_coupling` intentionally non-symmetric).
Quaternion/polar/corotational tests (`test_quaternion_*`, `test_polar_decomposition`,
`test_log_strain_small_deformation`, `test_corotational_frame_update`,
`test_frame_incremental_rotation`) are architecture-neutral but currently unused by
the linear MITC4+ path — decide per F4/F6.

### 2H. Assembly / Python invariants that constrain the new element

- `test_rust_assembler.py` (18 tests) — COO vs PETSc, symmetry, `K_T(0)=K`,
  `fint(0)=0`, NR consistency. Must hold for the new element.
- `test_mass_matrix_validation.py` (15) — exact mass coefficients; must hold.
- `test_quad_elements.py` (21) — **plane** QUAD4/8/9, 2 DOF/node; not applicable to
  the shell element.

---

## 3. Integration surface the new element must fit

### 3.1 How an element is declared and reached

| Layer | Exact location | Current assumption |
| --- | --- | --- |
| Rust module | `crates/aeroelast-core/src/elements/mod.rs` (`pub mod mitc4;`) | add `pub mod mitc4_plus;` |
| Element type enum | `crates/aeroelast-core/src/assembly/topology.rs` `enum ElemType`, `impl ElemType::dofs_per_node()` | `Mitc3|Mitc4|Mitc3Composite|Mitc4Composite => 6`; new variant needs a `5` arm |
| Global stride | `MeshTopology` (`topology.rs:58`), inferred as `max(dofs_per_node())` (`:84`); `global_dof_indices` uses `elem_dofs = elem_type.dofs_per_node()` (`:107`) | all nodes share one stride; mixed meshes assume the max |
| Precompute | `assembly/assembler.rs` `enum PrecomputedElem::Quad(Mitc4Precomputed)`, `build_constitutive_mitc4`, `MeshAssembler::new`/`update_geometry`/`update_reference` | `Mitc4Precomputed::new` + `Quad`; `update_reference` hardcodes `let dofs_per_node = 6;` (`assembler.rs:224`) |
| Element consumption | `assembler.rs`: `assemble_k` (K), `assemble_m`/`assemble_m_lumped`, `assemble_kt` (`extract_elem_disp_24`, `:1079`), `assemble_fint`, `assemble_geometric_k`, `compute_element_stress` | `mitc4::compute_ke_global`, `compute_kt_global`, `compute_fint_global` all consume `mitc4::Vec24` = `SVector<f64,24>` |
| PyO3 element kernels | `crates/aeroelast-py/src/elements.rs` `batch_ke_mitc4`/`batch_me_mitc4`/`batch_kt_mitc4` return `[f64;576]` (24×24), `batch_fint_mitc4` returns `[f64;24]`; `batch_ke_mitc4_composite` etc. | hardcoded 12 coords, `24`, `576` |
| PyO3 registration | `crates/aeroelast-py/src/lib.rs` `register_module` (lines 30–37) | one symbol per element function |
| PyO3 class | `crates/aeroelast-py/src/assembler.rs` `PyMeshAssembler`; `materials.rs` `PyElementFamily` (SHELL=2, PLANE=3) | element type codes 3/4/33/44 mapped in `new` and `from_model` |
| Python family | `src/aeroelast/core/assembler.py` `_FAMILY_PROPERTIES = {SHELL:(6,3), PLANE:(2,2)}` (line 164), `_FAMILY_VECTOR_FORM = {SHELL:{U:(Ux,Uy,Uz), θ:(θx,θy,θz)}}` (170), `_node_dofs_map` (`:258`, stride `dpn`), `_build_py_mesh_assembler` legacy path (`:271`) and composite path (`:288`) | `dofs_per_node=6` for shells throughout |
| Python package | `src/aeroelast/elements/__init__.py` | only re-exports `ElementFamily` — no per-element Python classes exist |

### 3.2 Where local stiffness/force/tangent are consumed

- Local → global: `mitc4::compute_ke_global = t3·K_local·t3ᵀ` via `build_t24`
  (`mitc4.rs:1971`) and `transform_to_global` (`:1985`). `t3` is a 3×3 rotation
  built from the element's mean normal `e3` and edge `e1`
  (`compute_local_coordinate_system`). `build_t24` block-diagonalises `t3` over
  the 6-DOF-per-node layout (3 translations + 3 rotations).
- `compute_ke_local` (`mitc4.rs:1405`) is **6-DOF-specific**: `b_m_mitc4_plus`
  (3×24), `b_gamma_mitc4`, `b_drill` (1×24), the 26-DOF bubble, the SRI centre
  point, `k_drill`.
- Stress recovery: `assembler.rs` `PrecomputedElem::Quad(pre) => extract_elem_disp_24`
  → `mitc4::compute_element_stress`.

### 3.3 What must change for 3 translations + 2 director rotations (5 DOF)

Exact hardcoded assumptions:

1. **`ElemType::dofs_per_node()` → 6** (`topology.rs:31`). Needs a 5 arm; the
   global stride (max over elements) then becomes non-uniform if shells of both
   kinds coexist.
2. **`update_reference`: `let dofs_per_node = 6;`** (`assembler.rs:224`) — updates
   node coords from `u_inc[6i..6i+2]`; must become stride-aware (and the rotation
   update convention must be defined).
3. **`extract_elem_disp_24` / `Vec24` / `Mat24` / `[f64;576]` / `[f64;24]`** — the
   whole PyO3 + `MeshAssembler` kernel path assumes 4×6. A 5-DOF element is 20×20
   (400) with a 20-vector.
4. **`build_t24` / `transform_to_global`** — assumes 3 global rotations per node.
   For 5 DOF the local rotation space is the per-node `(V_1,V_2)` pair; the
   transform is no longer a single block-diagonal `t3` but a per-node 2×3 (or
   3×2) director frame map, and the frame uses the *nodal director* `V_n^i`, not
   a single element `e3`.
5. **Rotation-vector convention.** Current element (`shell-elements.md` §1.3) uses
   a single physical rotation vector `θ` with `γ_xz = w_,x+θ_y`,
   `γ_yz = w_,y−θ_x`. Paper A uses `α_i` (about `V_1^i`) and `β_i` (about
   `V_2^i`) in `u_b = ½Σ a_i h_i (−V_2^i α_i + V_1^i β_i)` (extract Eq. 8b).
   The two are not the same parameterisation and the sign/frame mapping must be
   stated explicitly.
6. **Python `_FAMILY_PROPERTIES` / `_FAMILY_VECTOR_FORM`** — `(6,3)` and the
   `θ=(θx,θy,θz)` output form (VTK) assume 6. A 5-DOF element has no `θ_z`
   output slot.
7. **Per-node geometry (`V_n^i`, `a_i`)** does not exist in the current element;
   the current code uses a single element normal and a scalar thickness. The new
   element must ingest nodal directors and nodal thicknesses through
   `Mitc4Precomputed::new`-equivalent, and through `PyMeshAssembler` (which today
   passes only 12 floats: 4×3 coords).

---

## 4. Gap and risks

### 4.1 What exists today that the new element must NOT inherit

The deviations list (`odd/tasks/mitc4plus-faithful.md`; `docs/formulations/shell-elements.md` §2):

1. 6 DOF/node with a drilling rotation (`θ_z`).
2. Flat projection onto a local 2D frame + ABD resultants, instead of the 3D
   continuum kinematics with director `V_n^i` and thickness coordinate `t`.
3. A separate 2-DOF `(1−ξ²)(1−η²)` rotation bubble, instead of the per-node
   director enrichment `u_b` (Eq. 8b) and `x_b` (Eq. 8a).
4. Selective reduced integration of the in-plane shear (`compute_ke_local`'s
   `cm_normal`/centre-point split).
5. Membrane field of Eqs. (18)–(19) (`b_m_mitc4_plus`, installed by `dc7593e`)
   instead of the new field Eqs. (21)–(27).
6. The shear correction factor `5/6` in the constitutive path.
7. Winkler & Plakomytis ERC (`compute_ke_local_erc`, `b_erc`,
   `erc_covariant_to_local`) + the `β_w` warping penalty.
8. Uncited `k_drill = 0.15·E·h²·drilling_scale` penalty (`shell-elements.md` §2.6).
9. Dead code confirmed by the compiler (`b_md_mitc4_plus`,
   `drill_midside_shape_derivatives`, `b_m_standard`, `green_lagrange_strain`,
   `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`) plus the unused
   hourglass scaffolding (`hg_factor`, `h_vec`, `hg_stiffness_factor`, `h_orth`).
   Prefer a new module; leave the old one intact as reference until F4/F5.

### 4.2 What must be built new

- `elements/mitc4_plus.rs`: geometry/kinematics (`x_r,x_s,x_d,n,m^r,m^s`,
  `c_r,c_s,d`, `a_A..a_E`), the 20-DOF layout, membrane (Eqs. 15–16, then 17–27),
  displacement bending (7c/7d, incl. `∂x_b`), MITC4 shear, 2×2×2 Gauss,
  plane-stress material — **no numerical factor**.
- A 20-DOF local→global transform with nodal directors.
- PyO3 `batch_*_mitc4_plus` (400/20) + registration.
- Rust `ElemType` variant with `dofs_per_node()=5`; `PrecomputedElem` arm;
  stride-aware `update_reference`; `[f64;400]`/`[f64;20]` kernel path.
- Python family/vector-form entry for the 5-DOF shell (or a documented 6-DOF
  embedding, §4.3).
- Tests: 5-DOF equivalents of §2G, plus the Tier-1 T1.1–T1.3.

### 4.3 Where 5-vs-6 bites hardest, and the options (no choice made here)

The whole repository is built on a **uniform 6-DOF/node global stride** for
shells: `ElemType::dofs_per_node()`, `MeshTopology.dofs_per_node`,
`_FAMILY_PROPERTIES[SHELL]=(6,3)`, the `θ=(θx,θy,θz)` VTK form, all
`test_ko2017_performance.py` (`DOF=6`), `test_shell_convergence.py`
(`_DOFS_PER_NODE=6`), and the 24-wide Rust kernels. A paper-faithful element has
**20 local DOF and no drilling**; a global system with an unconstrained drilling
DOF is singular if that DOF gets zero stiffness. Options:

- **A. Keep the 6-DOF global stride; contribute zero drilling stiffness.** Smallest
  change to assembly/Python. Risk: the global `K` has a zero row/column at every
  `θ_n` (a mechanism) unless the user constrains drilling or the neighbours
  couple it. The current element avoided this with the `k_drill` penalty — which
  is exactly the deviation.
- **B. Global stride 5 (true 20 DOF).** Matches the paper exactly. Cost: touches
  `ElemType`, `MeshTopology`, `update_reference`, every PyO3 batch signature,
  `assembler.py`, VTK output, and the DOF indexing of every pure-shell test;
  breaks mixed shell/plane meshes (max-stride assumption) unless reworked.
- **C. 6-DOF global stride; express the 5 local DOF as 2 rotations in the per-node
  `(V_1,V_2)` frame, and let `θ_n` be a decoupled zero-stiffness DOF handled by BC
  or a static condensation / explicit constraint.** A middle path: kernels gain a
  20×24 (or 20×6·n) operator, assembly unchanged. Cost: the mapping is
  element-local and per-node, and the `θ_n` DOF still needs a policy.
- **D. Port `mitc4.rs` in place to 5 DOF and absorb the test/DOF churn.** Not a
  "new element", but avoids two element families; heavier risk to a large passing
  suite.
- **E. Keep 6 DOF and add the drill DOF penalty-free later (the 2025 MITC4/D path,
  plan sketch F6).** Out of scope for this change; noted as the intended successor.

Each has a different blast radius across Tier 2 (§2) — in particular the
`DOF=6` indexing in `test_ko2017_performance.py` and the 6-DOF clamps in
`test_shell_convergence.py`/`test_large_rotation_benchmarks.py`.

### 4.4 Existing tests whose expectation would have to change (and why it is legitimate)

- **`test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple`**
  (currently **failing**): a correct 5-DOF element should make this pass; the fix
  is the formulation, not a loosened tolerance. No weakening.
- **`test_shell_convergence.py::test_in_plane_bending_convergence`** (currently
  **failing**): the tolerance is derived in its own docstring (Timoshenko shear
  floor ≈ 0.8% on 1230 µm, window 2%). If the new element changes the strip
  stiffness, `MIN_ORDER`/`EXTRAPOLATED_TOL` may need re-derivation against the
  *same physical reference*. That is legitimate only if the new value is justified
  by the reference, not by the measured output.
- **DOF-indexing constants** in `test_ko2017_performance.py` (`DOF = 6`),
  `test_shell_convergence.py` (`_DOFS_PER_NODE = 6`), and
  `test_large_rotation_benchmarks.py` (`elem_types=[4]`, 6-DOF `update_reference`):
  these are mechanical consequences of the 5-DOF choice, not tolerance changes.
- **Published benchmark cells already flagged as mis-sourced** in
  `validation-matrix.md` §3 (`test_3_2` SS rows use Table 6; `test_3_3[dist]`
  hardcodes a Table-12 value; `test_3_7[reg]` uses N=8/S4 cells; `test_3_5` docstring
  vs no xfail). Correcting these to the true paper cells **strengthens** the tests.
- **`test_b_erc_drill_row_warping_residual_is_characterized`**: pins the warped
  drill defect in `[1e-4,1e-2]`. The new element must not inherit it; the test
  should not be ported (its guard is meaningless without the ERC).

### 4.5 What could NOT be determined, and why

1. **Paper A Section 4's quantitative content** (isotropy error, zero-mode
   eigenvalues, patch-test meshes/values). Only the one sentence on p. 410 is in
   the extract; Figures 5+/Tables of Section 4 were never transcribed. This phase
   had **no PDF tooling** (`pdftoppm` unavailable), so paper A was not re-read.
   The T1.1–T1.3 checks in §1 are the plan sketch's design, not quoted paper
   numbers.
2. **Exact current pass/fail per test.** No shell. The session baseline was taken
   as given (Rust 120/0; Python 345/2/2) and `docs/validation-matrix.md` is a
   snapshot of `b2c62ff` (112/31/…) that disagrees with the baseline — treat the
   matrix as reference provenance, not as current status.
3. **The current thin twisted-beam value at N=8 / N=16.** The extract's handoff
   records `0.7313` at N=8 (quad) while the test comment records `0.9982` at N=16
   with the ERC; `_TWISTED_BEAM_CASES` currently expects the paper's
   `0.9975/0.9980` at N=16 with `tol=0.01` and has no active xfail. Which number
   the working tree produces today could not be measured here.
4. **Whether the "2 xfailed" figure (older baseline) still applies.** No active
   `xfail` marker was found in `tests/`; only the empty `_TWISTED_BEAM_CASES`
   machinery. Recorded as a discrepancy.
5. **The exact 5-DOF global-assembly design** (options A–E in §4.3) — a proposal
   decision, deliberately not made here.
6. **How `update_reference` / corotational geometry should map the 5 DOF** — not
   derivable from source alone; requires the formulation decision.

---

## 5. Session baseline vs `validation-matrix.md` (staleness note)

The matrix (§2) records `test_ko2017_performance.py` 31/31 pass, all sheets
green, at commit `b2c62ff`. The session baseline (given) reports
`345 passed / 2 failed / 2 skipped` with `test_batch_ke_mitc4_multiple` and
`test_in_bending_convergence` failing, and Rust 120 passed. The matrix's suite
snapshot (364 passed) is therefore **stale**. Provenance in the matrix (which
table cell each expectation came from, and the flagged mis-sourcings) remains
valid; its pass/fail columns do not.

---

## 6. Result contract

- **status:** ok
- **artifact:** `openspec/changes/mitc4plus-2017-faithful/explore.md`
- **tier_1:** 3 tests (isotropy, zero-energy modes, patch tests) — paper A §4;
  quantitative cells not transcribed.
- **tier_2:** 9 benchmark cases (paper B Tables 2–19), 4 convergence/analytical
  clusters, CCX parity (4 files), Rust `mitc4.rs` invariants (~25 tests) + 18
  assembly tests; 2 currently failing, 0 xfail.
- **integration:** 6-DOF stride hardcoded in `topology.rs`,
  `assembly/assembler.rs:224`, `aeroelast-py/src/elements.rs` (`[f64;576]`),
  `aeroelast-py/src/assembler.rs`, `src/aeroelast/core/assembler.py:164`.
- **gap:** 9 deviations the new element must not inherit; 7 build items; 5 options
  for the 5-vs-6 DOF mismatch (unchosen).
- **could_not_determine:** paper A §4 numbers (no PDF tooling), exact per-test
  status (no shell), current twisted-beam value, 5-DOF assembly design.
