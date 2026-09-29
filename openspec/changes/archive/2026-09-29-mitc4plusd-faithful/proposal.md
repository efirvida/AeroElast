# Proposal — the MITC4+/D element, faithfully implemented

Change: `mitc4plusd-faithful` · store: `openspec` + Engram mirror · phase: **proposal**

## 0. Intent

The repository does not currently know which element it implements. `crates/aeroelast-core/src/elements/mitc4.rs`
is a 6-DOF flat-projection element that mixes ingredients from at least four
different papers (Choi–Paik membrane, Winkler & Plakomytis ERC, an uncited
`0.15·E·h²` drilling penalty, a shear correction factor, a selective reduced
integration fix) with pieces of the MITC4+ formulation. The goal, in the user's
own words:

> "el goal es hasta que pase todos los test teoricos, primero los test de su mismo
> articulo, luego los que ya tenemos de otras referencias."

This change delivers **one element, built from its papers, with nothing the papers
do not have**, and makes its identity and acceptance criteria explicit and
testable. Nothing is inherited from the current hybrid. Acceptance is
test-driven in two tiers: (1) the papers' own basic theoretical tests, then
(2) the repository's existing theoretical tests from the other references.

**Non-negotiable outcome:** after this change, the repository can state exactly
which element it runs, against which paper section and equation, and can prove it
with tests that fail if an ingredient drifts.

## 1. The intended behavior, with the reference that defines it

**The element is the MITC4+/D of Ko, Bathe & Zhang (2025):** the 2017 MITC4+
formulation **plus** the 2025 penalty-free drilling DOF. There is no other element
in scope. Both papers are held locally:

- **Paper A (formulation):** Ko, Y., Lee, P.-S., Bathe, K.-J., *"A new MITC4+
  shell element"*, **Computers and Structures 182:404–418, 2017** —
  `.sources/papers/A_new_MITC4+_shell_element.pdf`. Transcription with page
  citations: `docs/formulations/mitc4plus-2017-extract.md`.
- **Paper B (2025 drilling DOF + tests):** Ko, Y., Bathe, K.-J., Zhang, X.,
  *"Continuum mechanics-based shell elements with six degrees of freedom at each
  node — the MITC4/D and MITC4+/D elements"*, **Computers and Structures
  308:107622, 2025** — `.sources/papers/1-s2.0-S0045794924003511-main.pdf`.
- **Shear field:** Dvorkin, E. N., Bathe, K.-J. (1984), *Eng. Comput.* 1(1):77–88,
  as reproduced in paper A (p. 405).
- **Benchmarks (Tier 2 only, not the formulation):** Ko, Y., Lee, P.-S., Lee, K.,
  Bathe, K.-J., *"Performance of the MITC3+ and MITC4+ shell elements in
  widely-used benchmark problems"*, **C&S 193:187–206, 2017**.

### 1.1 Ingredient → reference, equation by equation

| # | Ingredient | Reference and equation |
| --- | --- | --- |
| 1 | 3D continuum kinematics, nodal thickness `a_i`, nodal director `V_n^i` | Paper A, **Eqs. (1)–(3)**, p. 405 |
| 2 | DOF: 3 translations + 2 director rotations `α_i` (about `V_1^i`), `β_i` (about `V_2^i`). **No drilling rotation in the 2017 element** | Paper A, **Eq. (3)**, p. 405 (the rotation about `V_n^i` does not appear) |
| 3 | Covariant 3D strain `e_ij = ½(g_i·u_j + g_j·u_i)`, `r_1=r, r_2=s, r_3=t` | Paper A, **Eqs. (4)–(5)** |
| 4 | In-plane decomposition `e_ij = e_ij^m + t·e_ij^b1 + t²·e_ij^b2` | Paper A, **Eq. (7a)**, p. 406 |
| 5 | **Membrane** assumed field — the only modified term: sampled at the five tying points A(0,+1), B(0,−1), C(+1,0), D(−1,0), E(0,0) of Fig. 4; final field with `a_A..a_E` and `c_r = x_d·m^r`, `c_s = x_d·m^s`, `d = c_r²+c_s²−1` | Paper A, **Eqs. (17)–(27)** (17, 18, 19, 20, 21, 22, 23–25, 26, 27a–c), Section 3.2, p. 408–410 |
| 6 | **Bending** terms stay displacement-based (not assumed): `e_ij^b1`, `e_ij^b2` including their `∂x_b` parts | Paper A, **Eqs. (7c)/(7d)**, p. 406 ("we leave the other terms in Eq. (7a) as they are and evaluate them using the displacement formulation") |
| 7 | **Transverse shear** — MITC4 assumed field, tying points A(top, s=+1), B(bottom, s=−1), C(right, r=+1), D(left, r=−1): `e_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, `e_st = ½(1+r)e_st^C + ½(1−r)e_st^D` | Dvorkin & Bathe (1984) **Eq. (3)**, as reproduced in Paper A, p. 405 |
| 8 | **Drill-membrane strain** `e_ij^md` — the 2025 sixth-DOF contribution, so that `e_ij = ẽ_ij^m + e_ij^md + t·e_ij^b1 + t²·e_ij^b2` | Paper B, **Eq. (26)**, plus the equations the existing operator cites: **Eqs. (5), (10), (11a–b), (13b–c), (17a–b), (18), (19a–c), (21)** (`b_md_mitc4_plus`, `drill_midside_shape_derivatives` in `mitc4.rs`) |
| 9 | **Integration:** 2×2×2 Gauss over the element domain for the 2017 core | Paper A, p. 410: *"we use 2 × 2 × 2 Gauss integration over the element domain"* |
| 10 | **Integration:** the 2025 six-DOF elements require only **2×2 Gauss over the element surfaces** | Paper B, Section 3.1, pp. 13–14 |
| 11 | **No numerical factor anywhere** | Paper A, p. 410: *"the element formulation does not include any numerical factor"* |
| 12 | Plane-stress material in the 3D covariant frame (no shear correction factor) | Paper A constitutive path; consequence of #11 |

### 1.2 Identity statement (the one sentence this change exists to make true)

> The element is the **MITC4+/D**: the 2017 MITC4+ (paper A) — 3D continuum
> kinematics, assumed membrane Eqs. (17)–(27), displacement-based bending Eqs.
> (7c)/(7d), MITC4/Dvorkin–Bathe transverse shear, 2×2×2 integration, plane-stress
> material, **no numerical factor** — **plus** the 2025 penalty-free drill-membrane
> strain of paper B (Eq. 26 and the operator equations above), giving 6 DOF/node:
> 3 translations + 2 director rotations `(α_i, β_i)` + 1 drilling rotation.

The drilling DOF is a **paper-B ingredient, not a penalty**. The current
`k_drill = 0.15·E·h²·drilling_scale` and the Winkler & Plakomytis ERC treatment are
explicitly *not* part of this element.

### 1.3 Preserved invariants — the laminate/composite material models (hard constraint)

The user's constraint, verbatim: *"si vas a sustituir todo el codigo viejo, recuerda
que hay que mantener los modelos de material laminado"*. It is a **preserved invariant**
of this change, not an optional regression watch:

- `crates/aeroelast-core/src/materials/{laminate,orthotropic,composite,failure,isotropic}.rs`
  and the Python `src/aeroelast/core/laminate.py` + `src/aeroelast/constitutive/failure.py`
  SHALL keep their public behaviour and their numerical results. No laminate capability
  may be dropped, weakened or re-scoped by this change.
- The element consumes `ShellConstitutive { cm, cb_coupling, cb, cs, cm_raw }` — the A, B
  and D matrices plus the transverse shear — which is exactly how the composite path
  already feeds the element. That interface SHALL be preserved.
- Any existing test that exercises laminate/composite behaviour SHALL keep passing; those
  tests are part of Tier 2 (T2C and T2G, §2).

**A laminate capability regression is a change failure.** The retirement of the old
element (§4) may not be used to justify touching, simplifying or removing any of these
surfaces. Concretely, the invariant protects: (a) the Rust material modules named above,
(b) the Python laminate and failure modules named above, (c) the `ShellConstitutive`
interface and its five fields, and (d) the laminate/composite tests in Tier 2 — the
`TestOrthotropicSinglePly` / `TestSymmetricLaminates` / `TestAsymmetricLaminates` groups
of `test_material_suite.py`, `test_composite_b_coupling.py`, `test_rust_composite.py`,
`test_orthotropic_shell_parity.py` and `test_composite_beam_parity.py`.

### 1.4 Interface constraint — uncorrected transverse shear (design decision, not a proposal decision)

The paper's element has **no numerical factor** (paper A, p. 410), so it needs the
**uncorrected** transverse shear stiffness `G·h`. The repository's `ShellConstitutive`
(`crates/aeroelast-core/src/materials/mod.rs`) is in `k·G·h` form — the struct's own doc
comment says `Transverse shear stiffness (2×2): force-shear strain [k·G·h] form` — and
`ShellConstitutive` does **not** carry the correction factor `k`; `Laminate` does
(`shear_correction_factor`). This is a real interface mismatch between the faithful
element and the existing constitutive contract.

This amendment records the constraint and the candidate mechanisms; **the mechanism is a
design decision, not a proposal decision**:

- (a) extend `ShellConstitutive` with the correction factor;
- (b) add an uncorrected transverse-shear accessor; or
- (c) divide by `k` at the element boundary.

The choice is left to the design phase, which must also keep §1.3's preserved invariant
intact (whatever the mechanism, `Laminate`'s public behaviour and numerical results must
not change).

Related **open point, also a design decision**: `Laminate::to_shell_constitutive_with_offset`
exists and applies a reference-surface offset, while the paper's kinematics assumes the
reference surface **is the mid-surface**. Whether the element path may use the offset
variant, or must be restricted to the mid-surface form, is a design-phase decision.

## 2. Acceptance criteria

Acceptance is the user's two tiers, stated as testable statements. **Tier 1 is the
gate:** Tier 2 is only consulted once Tier 1 passes. No benchmark number may be
used to justify a Tier-1 failure.

### 2.1 Fixture policy for the Tier 1 patch tests (resolves the spec's G1/G3)

The two papers do **not** publish the patch-test boundary-node coordinates numerically:
paper A's **Fig. 5** gives the 10×10 square with interior nodes (2,2), (4,7), (8,7),
(8,3); paper B's **Fig. 7** gives the extension/bending/shearing patch meshes and the
minimum boundary conditions, including *"θ_z free at the element nodes except at the
corner node B"*. The policy for this change SHALL be:

1. Read both figures with **vision** (`pdftoppm`, never `pdftotext`) and record in the
   change what was actually read, **including any coordinate that cannot be determined**.
2. Where a boundary-node coordinate cannot be read reliably from the figure, use a
   **documented equivalent distorted patch** and record the deviation explicitly. The
   requirement is then verified against the documented equivalent, not against a silently
   invented mesh.
3. **Never invent a mesh, and never edit the reference to match the code.**

Which Tier 1 requirements are verified against which fixture:

| Requirement | Fixture |
| --- | --- |
| **T1.3** membrane/bending/shearing patch (paper A §4, pp. 410–411) | **Figure-read mesh** of paper A **Fig. 5** (10×10, interior nodes (2,2), (4,7), (8,7), (8,3)); where a boundary coordinate is unreadable, the **documented equivalent distorted patch** is used and the deviation recorded. |
| **T1.6** strong-form extension/bending/shearing patch (paper B §3.1, pp. 13–14) | **Figure-read mesh** of paper B **Fig. 7** and its minimum BCs (`θ_z` free except at corner B); where a boundary coordinate is unreadable, the **documented equivalent distorted patch** is used and the deviation recorded. |
| T1.1/T1.2 and T1.4/T1.5 (isotropy, exactly-six zero modes) | No patch fixture: single-element tests. |
| T1.3's paper condition Eq. (22) | Analytic flat-element unit test; no fixture. |

### Tier 1 — the papers' own basic theoretical tests

#### Tier 1a — Paper A (2017), Section 4 "Basic numerical tests", pp. 410–411

- **T1.1 Isotropy.** The element stiffness is invariant under the element's
  orientation and the node-numbering sequence: building the same element in
  several global orientations and co-rotating the DOF yields the same physical
  stiffness (eigenvalues / `uᵀKu` invariant to round-off). Paper A: *"The element
  passes the test of spatial isotropy."*
- **T1.2 Zero-energy modes — exactly six.** On a **single unsupported element**,
  the count of zero eigenvalues of `K` is **exactly 6**, corresponding to the six
  rigid-body modes; no spurious zero mode. Paper A: *"For the new element only the
  six zero eigenvalues corresponding to the six rigid body modes are obtained."*
  Checked on flat **and** warped geometry.
- **T1.3 Three patch tests — membrane, bending, shearing.** On the **distorted**
  mesh of paper A **Fig. 5** (10×10 square, interior nodes at (2,2), (4,7), (8,7),
  (8,3)), with the **minimum number of constraints** to prevent rigid-body motion
  and boundary nodal forces from the constant stress states. Passed if the correct
  **constant stress fields** are recovered anywhere in the mesh. Plus the paper's
  own patch-test condition **Eq. (22)** (`ẽ_rs^m|bil = e_rs^m|bil` for flat
  geometry) as a unit test, and the reduction test: for a flat rectangle
  `x_d = 0` ⇒ `c_r = c_s = 0` ⇒ all five coefficients vanish and Eq. (27) reduces
  to Eq. (18).

#### Tier 1b — Paper B (2025), Section 3.1 "Basic tests including the patch tests", pp. 13–14

- **T1.4 Spatial isotropy (2025).** The MITC4+/D passes the spatial isotropy test
  as stated in paper B: *"MITC/D and MITC4+/D elements pass the spatial isotropy
  test."*
- **T1.5 Zero-energy modes (2025).** The 2025 element passes the zero-energy-mode
  test and the rigid-body modes are properly represented. As in T1.2, the count
  must remain **exactly 6** — the drilling DOF must not add a seventh zero mode.
- **T1.6 Strong-form patch tests (2025).** Paper B's own patch tests use the
  **'strong form'**: the calculations must give the analytical solutions of
  **constant AND zero stresses throughout the patch** due to the applied loading,
  with the **'minimum boundary conditions'** imposed to prevent rigid-body motion;
  for each test the solution is the constant straining mode. Figure 7 of paper B
  gives the extension, bending and shearing patch tests **including the drill
  DOF**: **`θ_z` is left free at the element nodes except at the corner node B**,
  and `θ_z = 0` at corner node C does not affect the results. Three tests
  (extension, bending, shearing), each asserted for constant *and* zero stress
  states.

### Tier 2 — the repository's existing theoretical tests that must keep passing

Tier 2 means: **the existing theoretical suite keeps passing after the element
change**, with the two pre-existing failures reported as pre-existing and **not**
introduced or tolerated as new regressions. Names below are exact, as enumerated
in `openspec/changes/mitc4plusd-faithful/explore.md`.

**Suite baselines (must not regress):** `cd crates && cargo test -p aeroelast-core`
→ **120 passed / 0 failed**; `python -m pytest -m "not slow" -q` → **345 passed /
2 failed / 2 skipped**, the 2 failures being pre-existing and out of scope (§5).

**T2A — patch and rigid-body invariants** (`crates/aeroelast-core/src/elements/mitc4.rs`,
module `tests` from line 2188; the new element must reproduce the invariant in its
6-DOF layout):

- `tests::test_ke_global_leaves_all_six_rigid_body_modes_free`
- `tests::test_ke_global_has_exactly_six_zero_modes`
- `tests::test_membrane_patch_reproduces_constant_strain_at_every_gauss_point`
- `tests::test_bending_patch_reproduces_constant_curvature_at_every_gauss_point`

  *Note:* explore §2A also lists `tests::test_b_erc_membrane_rows_annihilate_rigid_body_modes`,
  `tests::test_b_erc_drill_row_is_rigid_body_invariant_on_flat_geometry` and
  `tests::test_b_erc_drill_row_warping_residual_is_characterized`. Those three test
  the **ERC deviation**, and explore §2G/§4.4 explicitly says they must **not** be
  ported — they are deleted with the old element. They are therefore *not* Tier 2
  acceptance; the underlying invariants (RB annihilated, exactly six zero modes)
  are covered by T1.2/T1.5 and the four tests above.

**T2B — eigenvalue / symmetry / tangent consistency:**

- `tests::test_ke_global_is_symmetric`, `tests::test_ke_global_is_positive_semidefinite`,
  `tests::test_ke_local_flat_plate_parity`
- `tests::test_kt_fint_directional_derivative`, `tests::test_kt_fint_directional_derivative_rotations`,
  `tests::test_kt_fint_directional_derivative_with_drill_dofs`
- `tests::test_fint_linear_nonlinear_parity`, `tests::test_kt_zero_matches_ke`
- `test_rust_composite.py::TestMITC4BatchSanity::test_ke_shape_and_symmetry`,
  `...::test_ke_positive_semidefinite`, `...::test_me_shape_and_symmetry`,
  `...::test_me_positive_semidefinite`
- `test_material_suite.py::TestStiffnessProperties::test_k_positive_definite`,
  `test_material_suite.py::TestStiffnessProperties::test_k_symmetric`
- `test_rust_assembler.py::TestTangentStiffness::test_kt_at_zero_equals_k`,
  `...::test_kt_symmetric`
- `test_rust_assembler.py::TestInternalForces::test_fint_zero_at_zero`,
  `...::test_fint_linear_equals_ku`
- `test_rust_assembler.py::TestNewtonRaphsonConsistency::test_nr_consistency`

  *Excluded:* `tests::test_ke_local_eigenvalues_nonsymmetric` — explore §2G says it
  uses `cb_coupling` intentionally non-symmetric and is a 6-DOF deviation test, not
  a paper invariant.

**T2C — analytical shell validation:**

- `test_material_suite.py::TestIsotropicAnalytical::test_mitc4_out_of_plane`,
  `...::test_mitc4_in_plane_lateral`, `...::test_mitc4_axial_stiffness`
- `test_material_suite.py::TestOrthotropicSinglePly` (6 tests),
  `TestSymmetricLaminates` (7 tests), `TestAsymmetricLaminates` (6 tests)
- `test_shell_analytical_validation.py::TestCantileverBeam::test_tip_deflection[4x2|8x4|16x8]`,
  `TestSimplySupportedBeam::test_center_deflection`,
  `TestBeamBending::test_bending_convergence`,
  `TestMembraneStretching::test_axial_extension`,
  `TestShearLocking::test_thin_plate_convergence`
- `test_shell_comprehensive.py` (fx/fy/fz, ratio, modal);
  `test_shell_validation_fixed.py` (fx/fy/fz/ratio/axial convergence/modal)
- `test_large_rotation_benchmarks.py::test_linear_tip_deflection_euler_bernoulli`,
  `...::test_cantilever_large_rotation_half_circle`, `...::test_equilibrium_path`,
  `...::test_simo_vu_quoc_rollup_360`
- `test_rust_modal.py::TestSimplySupportedPlate::test_analytical_convergence`, plus
  the `test_rust_modal.py` MITC4 frequency/match groups
- `test_mass_matrix_validation.py::TestExactConsistentMassCoefficients`,
  `...::TestConsistentMassTotal`, `...::TestElementMassVsTotalMass`
- CalculiX parity: `test_isotropic_shell_parity.py::TestIsotropicShellParity::test_transverse_tip_displacement`,
  `test_beam_shell_4cases_parity.py` (9), `test_composite_beam_parity.py` (5),
  `test_orthotropic_shell_parity.py` (3) — reference CalculiX 2.23 S4/S8R with
  `*SHELL SECTION, COMPOSITE`; skipped when `ccx` is not on `PATH`

**T2D — convergence:**

- `test_shell_convergence.py::test_in_plane_bending_convergence`
- `test_shell_convergence.py::test_composite_laminate_gap_mesh_study` (no value
  asserted — cannot fail)

**T2E — published benchmark columns, paper B (Ko, Lee, Lee & Bathe 2017),
`tests/test_ko2017_performance.py`, 31 cases, `rtol=0.05`:**

- `test_3_1_square_plate_tables_2_to_5`
- `test_3_2_circular_plate_tables_6_to_7`
- `test_3_3_pinched_cylinder_tables_8_to_9`
- `test_3_4_scordelis_lo_tables_10_to_11`
- `test_3_5_twisted_beam_tables_12_to_13`
- `test_3_6_hook_table_14_minimal_fix`
- `test_3_7_hemisphere_cutout_tables_15_to_16`
- `test_3_8_full_hemisphere_table_17`
- `test_3_9_hyperbolic_paraboloid_tables_18_to_19`

**T2F — assembly / mass invariants:**

- `test_rust_assembler.py` (18 tests): COO vs PETSc, symmetry, `K_T(0)=K`,
  `fint(0)=0`, NR consistency
- `test_mass_matrix_validation.py` (15 tests)

**T2G — composite / orthotropic parity** (as enumerated in explore §2F):
`test_material_suite.py` laminate groups, `test_composite_b_coupling.py`,
`test_rust_composite.py` invariants, `test_orthotropic_shell_parity.py`,
`test_composite_beam_parity.py`.

**T2H — explicitly not applicable:** `test_quad_elements.py` (21 tests) covers the
2-DOF/node **plane** QUAD4/8/9, not the shell element.

**T2I — Rust unit tests written against the 6-DOF architecture** that need a
6-DOF-layout equivalent for the new element (explore §2G): the rigid-body-mode
construction, the patch tests, `ke_global` symmetry/PSD, `fint`/`kt` parity,
`me_global*`, `test_body_load_global_*`, `test_k_sigma_global_*`,
`test_centrifugal_prestress_*`, and the `compute_element_stress` tests. These are
**mechanical 6-DOF-layout equivalents**, not tolerance changes.

**T2J — retired with the old element** (explore §2G, "Do NOT port — delete with
the old element"): `test_drill_midside_derivatives_match_paper_eq11`,
`test_drill_membrane_operator_*` (5), `test_b_erc_*` (5),
`test_erc_covariant_to_local_*` (2),
`test_ke_local_erc_flat_rigid_body_and_symmetry`, `test_enhanced_drill_stiffness`,
`test_drill_warping_moment_flat_element`, `test_ke_local_eigenvalues_nonsymmetric`.
They encode the deviations this change removes.

## 3. The gap

### 3.1 What exists today

- `crates/aeroelast-core/src/elements/mitc4.rs` — a **6-DOF flat-projection
  element**: 2-DOF rotation bubble, selective reduced integration, assumed
  membrane from **Eqs. (18)–(19)** (the Choi–Paik starting field, **not** the
  paper's Eqs. 21–27), a Winkler & Plakomytis ERC drilling treatment plus a `β_w`
  warping penalty (a *different* element, *different* paper), a shear correction
  factor, and an uncited `0.15·E·h²` drilling penalty.
- `b_md_mitc4_plus` + `drill_midside_shape_derivatives` — the **2025 paper-B
  drill-membrane operator**, transcribed correctly (including the `1/8` factor of
  Eq. 13c and the `[right, top, left, bottom]` edge ordering) but **currently dead
  code**, never wired into any stiffness/force path. This change revives it as the
  element's sixth DOF, which is only possible once the SRI split is removed (the
  SRI is exactly what made the drill term inert — measured bit-identical even when
  amplified ×1000).
- Other compiler-confirmed dead code in the same file: `b_m_standard`,
  `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`.

### 3.2 What must be built new

- The **2017 core** as the papers write it: geometry/kinematics (`x_r, x_s, x_d, n,
  m^r, m^s, c_r, c_s, d, a_A..a_E`), the membrane Eqs. (15)–(16) and the assumed
  field **Eqs. (17)–(27)** with `a_E = +2c_r c_s/d` (**positive** — the extract's
  §2.3 and the deleted `compute_membrane_coefficients` had it negative),
  displacement-based bending Eqs. (7c)/(7d) including `∂x_b`, MITC4 transverse
  shear (Dvorkin & Bathe 1984), **2×2×2** integration, plane-stress material, **no
  numerical factor**.
- The **2025 drill-membrane** contribution from paper B **Eq. (26)** with the
  operator equations the existing code already cites, wired into the stiffness,
  internal-force and tangent paths (and thus into the 2×2 element-surface
  integration of paper B §3.1).
- A 6-DOF local↔global transform consistent with the **nodal director** `V_n^i`
  and the per-node `(V_1^i, V_2^i)` basis (not a single element `e3`).
- Tier-1a and Tier-1b tests (T1.1–T1.6) and the Tier-2 equivalents listed in T2I.
- Documentation that states the identity: rewrite `docs/formulations/shell-elements.md`
  §2 and the `mitc4.rs` module header to describe **what actually runs**, and
  update `docs/validation-matrix.md`.

### 3.3 What must NOT be inherited (the deviations list)

Not inherited, in full:

1. 6-DOF **drilling penalty** (`k_drill = 0.15·E·h²·drilling_scale`) — uncited.
2. Winkler & Plakomytis **ERC** (`compute_ke_local_erc`, `b_erc`,
   `erc_covariant_to_local`) and the **`β_w` warping penalty**.
3. **Selective reduced integration** of the in-plane shear (the `cm_normal` /
   centre-point split).
4. The 2-DOF **`(1−ξ²)(1−η²)` rotation bubble** instead of the paper's director
   enrichment `u_b` (Eq. 8b) / `x_b` (Eq. 8a).
5. The membrane field of **Eqs. (18)–(19)** instead of **Eqs. (21)–(27)**.
6. The **shear correction factor** `5/6`.
7. **Flat projection + ABD resultants** instead of 3D continuum kinematics with
   `V_n^i` and the thickness coordinate `t`.
8. Any **numerical factor** — paper A, p. 410 states there is none.
9. The 2-DOF **hourglass scaffolding** (`hg_factor`, `h_vec`,
   `hg_stiffness_factor`, `h_orth`).

## 4. Scope boundaries and non-goals

**In scope:** the MITC4+/D element as §1 defines it, its Tier-1a/Tier-1b tests,
the Tier-2 suite staying green, the wiring into the existing 6-DOF kernel/assembly
path at **24×24**, the documentation that names the element, and — **gated on the new
element passing Tier 1 and Tier 2** — the retirement of the superseded hybrid from
`crates/aeroelast-core/src/elements/mitc4.rs`.

The retirement SHALL remove: the Winkler & Plakomytis ERC drilling treatment, the
`beta_w` warping penalty, the selective reduced integration, the 2-DOF
`(1−ξ²)(1−η²)` rotation bubble, the Eqs. (18)–(19) membrane, the shear correction
factor usage, the uncited `0.15·E·h²` `k_drill`, and the compiler-confirmed dead code
(`b_m_standard`, `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`,
`Mat26`/`Vec26`). It SHALL **not** remove `b_md_mitc4_plus` or
`drill_midside_shape_derivatives`: those are **revived** as the 2025 drill ingredient
(§3.1). The old element remains the reference until the gate passes, so the retirement
expands the change; the tasks must slice it into work units within the accepted
700-line budget.

**Explicit non-goals (out of scope, do not do):**

1. **No shell-element trait.** `crates/aeroelast-core/src/elements/mod.rs` is 5
   lines of `pub mod`; the assembler dispatches by a hardcoded `match` on
   `elem_types[e]`, and the PyO3 batch functions are per-element with the size
   baked into the type. Adding a shell-element trait is real debt but is **not on
   this change's critical path**.
2. **No migration of the assembly to a 5-DOF stride.** The element stays 24×24, so
   `ElemType::dofs_per_node()`, `MeshTopology`'s global stride, `update_reference`'s
   `let dofs_per_node = 6;`, `[f64;576]`/`[f64;24]` in `crates/aeroelast-py/src/elements.rs`,
   `_FAMILY_PROPERTIES[SHELL] = (6, 3)` and the VTK `θ=(θx,θy,θz)` form are
   **untouched**. The 5-DOF-only options A–D of explore §4.3 are closed by the
   user's decision; explore §4.3 option **E** is what this change implements.
3. **No touching `crates/aeroelast-solvers`.** It cannot build in this environment
   (the preCICE crate's build script needs the preCICE library).
4. **No fixing the two pre-existing Python failures.**
   `test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple` and
   `test_shell_convergence.py::test_in_bending_convergence` are pre-existing and
   must be reported as such, not as regressions of this change. If the new element
   incidentally makes them pass, that is recorded as a bonus, not as a requirement
   or an acceptance criterion.
5. **No parameter fitting and no tolerance loosening.** No tolerance may be
   widened to accommodate the new element; a Tier-1 failure is a formulation bug,
   never a tolerance problem.
6. **No new numerical factor of any kind**, and no re-introduction of any item in
   §3.3.
7. **No reading equations with `pdftotext`.** Papers are read with vision
   (`pdftoppm -png -r 300`); `pdftotext` has already produced wrong readings in
   this project.
8. **No running the twisted-beam benchmarks at N ≥ 32.**

## 5. Affected crates and build impact

| Layer | Crate / path | Impact |
| --- | --- | --- |
| Element kernels | **`crates/aeroelast-core`** — `src/elements/mitc4.rs` (and/or a new module), `src/elements/mod.rs` | the formulation itself; Tier-1 and Tier-2 Rust tests; **and, after the Tier 1 + Tier 2 gate, removal of the superseded hybrid code from `mitc4.rs` (§4)** |
| Assembly | `crates/aeroelast-core` — `src/assembly/assembler.rs`, `src/assembly/topology.rs` | `PrecomputedElem` arm + dispatch; stride **unchanged** (24) |
| PyO3 bindings | **`crates/aeroelast-py`** — `src/elements.rs`, `src/lib.rs`, `src/assembler.rs` | batch kernels stay `[f64;576]`/`[f64;24]`; symbol registration |
| Python family table | `src/aeroelast/core/assembler.py` | `_FAMILY_PROPERTIES[SHELL] = (6, 3)` **unchanged** |
| Python tests | `tests/` | Tier-2 suite; no tolerance changes |

**Retirement and the PyO3 surface:** the change now also **removes code from
`aeroelast-core`**, so a maturin rebuild is still required (below) and the retirement
SHALL NOT change the PyO3 surface — `crates/aeroelast-py/src/elements.rs` keeps
`[f64;576]`/`[f64;24]`, no `#[pyfunction]`/`#[pyclass]` symbol is removed or renamed,
and `src/aeroelast/core/assembler.py`'s `_FAMILY_PROPERTIES[SHELL] = (6, 3)` is
unchanged.

**Build impact: a maturin rebuild is required before any Python test.** The Python
suite drives the Rust extension through PyO3, so after any Rust change:

```bash
cd /home/efirvida/Desktop/dev/fem-shell && python -m maturin develop --release
python -m pytest -m "not slow" -q
```

Rust tests do not need the rebuild:

```bash
cd crates && cargo test -p aeroelast-core
```

## 6. Risks and tradeoffs

1. **`strict_tdd: false`** (`openspec/config.yaml`, decision-gate reason: no
   workspace-level test command; `crates/aeroelast-solvers` does not build here).
   Consequence: the apply phase will **not** be forced test-first automatically.
   The tasks phase must therefore carry an **explicit test-first instruction** —
   Tier-1a/Tier-1b tests written and failing **before** the formulation is
   implemented — or this change loses its only real safety net.
2. **The element is 24×24, so the assembly layers are untouched.** This is the
   change's main risk-reducer: the PyO3 signatures, the global stride, the Python
   family table and the VTK output form do not move, which is why the blast radius
   stays small. The corollary risk is that **no test currently proves the assembly
   layer is element-agnostic** — it is proven only by the kernel size not changing.
3. **The 2017 Tier-1 tests must still hold *with* the drilling DOF present.** This
   is the sharpest risk in the change:
   - **the zero-energy count must stay exactly six.** The drilling DOF must not
     introduce a seventh zero mode (a drilling mechanism) nor a spurious
     non-zero mode. T1.2 and T1.5 are both stated as exactly six for this reason.
   - the 2017 membrane/bending/shearing patch tests (T1.3) must pass **with**
     `e_ij^md` present, and the 2025 strong-form tests (T1.6, constant **and** zero
     stresses, `θ_z` free except at corner B) must pass simultaneously. Passing
     one and failing the other means the two ingredients were not merged correctly.
   - the 2017 element's six rigid-body modes are expressible with `α=β=0` for the
     drilling rotation (`ω × V_n = 0`); the 6-DOF mode construction must reproduce
     that mapping or T1.2 will be reported as a false failure.
4. **The drill-membrane strain as currently transcribed needs re-verification
   against the 2025 paper.** `b_md_mitc4_plus` was transcribed from paper B and its
   doc comment corrects two known errors of an earlier deleted version (the `1/8`
   factor of Eq. 13c and the `[right, top, left, bottom]` edge ordering), but it
   has **never been exercised by any stiffness path** — it was inert by
   construction (SRI). Reviving it means trusting a transcription whose only
   evidence so far is that it did not change a result *because it could not*.
   **Mitigation:** re-read paper B Eqs. (5), (10), (11), (13b–c), (17a–b), (18),
   (19a–c), (21) and Eq. (26) with vision and verify the operator term by term
   before wiring it in; the existing `b_md` unit tests are known to be weak
   (explore §4.4 / `shell-elements.md` §4.2 records that the earlier ones were
   vacuous because they evaluated at `(0,0)`, where every mid-side derivative
   vanishes by construction).
5. **Removing the SRI changes behaviour.** The drill term is only non-inert once
   the SRI split is gone, so this change is not purely additive: the 2017 core
   itself changes the element's numbers. Every Tier-2 benchmark that currently
   passes with the hybrid must be re-checked; the expectation is that a
   paper-faithful element still passes, but a value that moves and stays inside
   its published tolerance is a legitimate outcome, while a value that moves
   *outside* tolerance is a formulation bug.
6. **Documentation drift is part of the defect.** `docs/validation-matrix.md` is a
   snapshot of commit `b2c62ff` and its pass/fail columns are stale relative to the
   session baseline; it also carries flagged mis-sourcings in
   `test_ko2017_performance.py` (`test_3_2` SS rows use Table 6; `test_3_3[dist]`
   hardcodes a Table-12 value; `test_3_7[reg]` uses N=8/S4 cells; `test_3_5`
   docstring vs no xfail). Correcting those to the true paper cells **strengthens**
   the tests and is in scope for the documentation work; leaving them mis-sourced
   would let this change's success criteria rest on the wrong numbers.
7. **Explore-artifact naming drift (known, benign).** `explore.md` is titled
   `# Explore — mitc4plus-2017-faithful` and its §0/§4.3 scope a **new 5-DOF**
   element, with the 5-vs-6 DOF decision left open (options A–E). The change was
   subsequently renamed `mitc4plusd-faithful` and the user closed that decision in
   favour of option E (6 DOF, 2025 penalty-free drill). The explore artifact's
   inventory and evidence remain valid; its **scope framing and §4.3 options are
   superseded by this proposal**.
8. **Fixture gap replaced by an explicit policy (§2.1).** The two papers do not
   publish the patch-test boundary-node coordinates numerically, and the explore
   phase had no PDF tooling, so paper A §4's quantitative content (isotropy error,
   zero-mode eigenvalues, patch-test meshes/values) was never transcribed and no
   suite was run. This is no longer an open-ended gap: §2.1 fixes the policy — read
   the figures with vision and record what was read, fall back to a **documented
   equivalent distorted patch** where a coordinate cannot be read, and never invent a
   mesh or edit the reference to match the code. Paper B's Section 3.1 numbers were
   read with vision by the requester but are not yet in a persistent extract, so the
   tasks phase must still include a **paper-B extraction step** before the drill term
   is wired, or the Tier-1b criteria rest on a summary rather than on the paper.
9. **The retirement enlarges the change beyond the original estimate.** The original
   scope was additive (a new element plus its tests); retiring the old hybrid adds the
   deletion of the ERC / `beta_w` / SRI / bubble / `k_drill` / dead-code surfaces and
   the retirement of the T2J tests. It is gated on the new element passing Tier 1 and
   Tier 2 (§4), so the old element remains the reference until then, but the change now
   carries a removal step the original 700-line budget did not plan for.
   **Mitigation:** the tasks phase slices the retirement into its own work units and
   keeps the deviation-code deletion separate from the formulation work, so the
   retirement can be reviewed as a distinct unit and reverted on its own (§7).
10. **Tier 2 — especially the laminate/composite tests — is the safety net that proves
   the retirement did not cost a capability.** The retirement touches the same file the
   composite path feeds through `ShellConstitutive`; §1.3 makes the laminate/composite
   material models a preserved invariant and any laminate capability regression a change
   failure. The evidence that the retirement is safe is therefore Tier 2's
   laminate/composite groups (T2C, T2G: `test_material_suite.py`'s
   `TestOrthotropicSinglePly` / `TestSymmetricLaminates` / `TestAsymmetricLaminates`,
   `test_composite_b_coupling.py`, `test_rust_composite.py`,
   `test_orthotropic_shell_parity.py`, `test_composite_beam_parity.py`), not a code
   review of the deletion. If those groups do not stay green, the retirement is a
   failure regardless of the new element passing.

## 7. Rollback

- The new element is **additive**: a new module and a new `PrecomputedElem` arm.
  The old `mitc4.rs` stays intact as a reference until the new element passes
  Tier 1 and Tier 2 (plan-sketch constraint: *"nothing is deleted before then"*).
- Rollback is therefore: revert the element-type registration (one `match` arm and
  one `ElemType` variant), keeping the old element path and all Tier-2 tests as
  they were. No data migration, no API break, no assembly-format change.
- The retired-deviation tests (T2J) and the deletion of the ERC / `β_w` / SRI /
  bubble / `k_drill` / dead-code surfaces happen **only after the new element passes
  Tier 1 and Tier 2** (§4); until then they remain in place and green.
- Once the retirement has happened, rollback is a revert of the retirement commit:
  the deleted old-element code returns and the old element path is restored. Until the
  gate passes, rollback remains the registration revert of the bullet above. Either
  way: no data migration, no API break, no assembly-format change.

## 8. Success criteria

1. **Identity is provable.** The repository states, in one place, that the element
   is the MITC4+/D = paper A (2017) core + paper B (2025) drill, with the
   equation-level citations of §1.1, and a test fails if an ingredient drifts.
2. **Tier 1a passes:** isotropy, exactly six zero eigenvalues, and the membrane /
   bending / shearing patch tests on the Fig. 5 distorted mesh (T1.1–T1.3).
3. **Tier 1b passes:** 2025 spatial isotropy, zero-energy modes, and strong-form
   patch tests with constant **and** zero stresses, minimum BCs, `θ_z` free except
   at corner B (T1.4–T1.6).
4. **Tier 2 passes:** `cd crates && cargo test -p aeroelast-core` → 120 passed /
   0 failed; `python -m pytest -m "not slow" -q` → 345 passed / 2 failed / 2
   skipped with the same two pre-existing failures and no new ones.
5. **No deviation survives:** no ERC, no `β_w`, no SRI, no rotation bubble, no
   shear correction factor, no `k_drill` penalty, no numerical factor of any kind,
   and no live `Eqs. (18)–(19)` membrane path in the element that runs.
6. **No scope creep:** no shell-element trait, no 5-DOF stride migration, no
   `crates/aeroelast-solvers` work, no pre-existing-failure fixes.

## 9. Proposal question round

This executor could not put questions to the user directly (parent owns
interaction), so the proposed product/business round is recorded here for review.
These are the product-shaped unknowns that would make the proposal ambiguous or
easy to overbuild. The **current assumptions** are stated so the proposal is
actionable if the user skips the round.

| # | Question | Current assumption if unanswered |
| --- | --- | --- |
| Q1 | **Business problem / cost of not doing it.** Is the pain "we cannot trust our own shell results" (correctness for blade design), or "we cannot tell which element we run" (auditability), or both? Which one justifies the effort if only one is true? | Both: correctness first, auditability as the mechanism. |
| Q2 | **Who is affected and when.** Does this element need to become the **default** shell element for wind-turbine blade runs, or does it start as an opt-in element selectable next to the current one? | Opt-in: a new element that coexists; default switch is a later, separate decision. |
| Q3 | **Business rules / invariants the element must respect.** Are there blade-specific constraints (warped/curved geometry, thickness taper, composite layups, large rotations) that must be part of acceptance beyond the papers' patch tests? | No extra invariants in this change; blade fidelity is delegated to the existing Tier-2 analytical and CalculiX-parity tests. |
| Q4 | **Product outcome.** After this change, what should become possible that is not possible today — published-paper parity claims, certification evidence, a defensible element-identity statement? Which artifact is the deliverable? | A defensible identity statement plus Tier-1 evidence; published-benchmark parity is claimed only where Tier-2 already passes. |
| Q5 | **Edge cases and failure policy.** If the faithful element passes Tier 1 but moves a currently-passing Tier-2 benchmark outside tolerance, is the correct response (a) treat it as a formulation bug and block, or (b) investigate the benchmark's sourcing first (several are known mis-sourced)? | (b) first — correct the mis-sourced cells, then re-judge; (a) only after sourcing is verified. |
| Q6 | **Scope boundary on the old element.** Should the old hybrid element be **deleted** in this change, or kept selectable indefinitely as a fallback? | Kept as reference until Tier 1 passes; deletion of the deviation code is in scope, removal of the element family is not. |
| Q7 | **Business risk / tradeoff.** What downside matters most if this change chooses the wrong direction — a slower path to blade results, or an element that is faithful to the papers but worse on the blade use case? | Faithfulness first, by explicit user decision; performance regressions on benchmarks are a finding to report, not a reason to deviate. |
| Q8 | **Non-goal confirmation.** Confirm that `strict_tdd: false` is accepted, and that the tasks phase should carry an explicit test-first instruction for Tier 1 rather than relying on the harness (see risk §6.1). | Confirmed as stated in the preflight; test-first is instructed in tasks. |

## 10. Provenance

Inputs read for this proposal:

- `openspec/changes/mitc4plusd-faithful/explore.md` (446 lines) — the test
  inventory (Tier 1/Tier 2), integration surface, deviations list, and options;
  note §6.7 above on its superseded scope framing.
- `docs/formulations/mitc4plus-2017-extract.md` — paper A transcription with page
  citations (Eqs. 1–27, the sign correction `a_E = +2c_r c_s/d`, integration and
  no-numerical-factor statements).
- `crates/aeroelast-core/src/elements/mitc4.rs` — `b_m_mitc4_plus` (Eqs. 18–19,
  the Choi–Paik field), `b_md_mitc4_plus` + `drill_midside_shape_derivatives`
  (paper B Eqs. 5, 10, 11, 13b–c, 17a–b, 18, 19a–c, 21), and the recorded
  measurement that the drill term is inert under the current SRI.
- `openspec/config.yaml` — `strict_tdd: false` and its decision-gate reason,
  project test commands, the maturin requirement.
- `odd/tasks/mitc4plus-faithful.md` — plan sketch and the user's scope decision.
- `docs/formulations/shell-elements.md`, `docs/formulations/solvers.md` — the
  deviation and inertness records.

Paper B (2025) itself was **not** read in this phase; §1.1's paper-B citations
come from the requester's vision read and from the operator's own doc comments.
Per risk §6.4, re-verification against `.sources/papers/1-s2.0-S0045794924003511-main.pdf`
is a required tasks-phase step before the drill term is wired.
