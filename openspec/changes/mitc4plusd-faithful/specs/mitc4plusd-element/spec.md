# MITC4+/D Shell Element Specification

Change: `mitc4plusd-faithful` · domain: `mitc4plusd-element` · phase: **spec**
Status: new domain (no `openspec/specs/mitc4plusd-element/spec.md` existed before this change).

## Purpose

Define, testably, the element this repository runs: the **MITC4+/D** of Ko, Bathe &
Zhang (2025) — the 2017 MITC4+ of Ko, Lee & Bathe (2017) plus the 2025 penalty-free
drilling DOF. The spec fixes (a) the ingredient set and the exact absence of everything
the papers do not contain, (b) the DOF layout, (c) the five geometric membrane
coefficients, (d) the papers' own basic tests (Tier 1a for the 2017 core, Tier 1b for the
2025 six-DOF element), (e) the repository's existing theoretical tests that must keep
passing (Tier 2), and (f) the provenance requirement that every live ingredient names its
paper/section/equation. The amendment adds (g) the gated retirement of the superseded
hybrid, (h) the preserved-invariant requirement for the laminate/composite material
models, and (i) the uncorrected transverse-shear requirement with its offset deferral.

## Documents referenced

| Ref | Short name | Locator |
| --- | --- | --- |
| **A** | Ko, Lee & Bathe 2017, *A new MITC4+ shell element*, C&S 182:404–418 | `.sources/papers/A_new_MITC4+_shell_element.pdf`; transcription with page citations `docs/formulations/mitc4plus-2017-extract.md` |
| **B** | Ko, Bathe & Zhang 2025, *Continuum mechanics-based shell elements with six degrees of freedom at each node — the MITC4/D and MITC4+/D elements*, C&S 308:107622 | `.sources/papers/1-s2.0-S0045794924003511-main.pdf` |
| **DB84** | Dvorkin & Bathe 1984, *Eng. Comput.* 1(1):77–88 | transverse shear assumed field, as reproduced in A p. 405 |
| **BENCH** | Ko, Lee, Lee & Bathe 2017, C&S 193:187–206 | Tier 2 benchmark columns only; **not** a formulation source |

Papers MUST be read with vision (`pdftoppm -png -r 300`), never `pdftotext`.

## Verification commands

| Scope | Command |
| --- | --- |
| Tier 1a + Tier 1b + Rust Tier 2 | `cd crates && cargo test -p aeroelast-core` |
| One Rust test (substring match on the test path) | `cd crates && cargo test -p aeroelast-core <test_name>` |
| Tier 1a gate only | `cd crates && cargo test -p aeroelast-core test_t1a_` |
| Tier 1b gate only | `cd crates && cargo test -p aeroelast-core test_t1b_` |
| Retirement gate only | `cd crates && cargo test -p aeroelast-core test_retirement_` |
| Uncorrected-shear discriminator | `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_invariant_to_shear_correction_factor` |
| Laminate/composite preserved invariant | `python -m pytest "tests/test_material_suite.py::TestOrthotropicSinglePly" "tests/test_material_suite.py::TestSymmetricLaminates" "tests/test_material_suite.py::TestAsymmetricLaminates" "tests/test_composite_b_coupling.py" "tests/test_rust_composite.py" "tests/test_orthotropic_shell_parity.py" "tests/test_composite_beam_parity.py" -q` |
| Python Tier 2 (requires the extension rebuild first) | `cd /home/efirvida/Desktop/dev/fem-shell && python -m maturin develop --release && python -m pytest -m "not slow" -q` |
| One Python test | `python -m pytest "tests/<file>.py::<Class>::<test>" -q` |

Test names marked **NEW** do not exist yet and are created by this change; the enclosing
module path is a design decision, so the `cargo test` invocation above uses the
substring form that is insensitive to it. Names marked **EXISTING** are in the
repository today and MUST NOT be renamed.

Any numeric tolerance stated in this spec is a fixed acceptance threshold. Loosening a
tolerance to make a Tier-1 requirement pass is a spec violation (a Tier-1 failure is a
formulation bug, never a tolerance problem).

## Requirements

### Requirement: Element identity and formulation

The element SHALL be the **MITC4+/D**: the 2017 MITC4+ core of paper A plus the 2025
penalty-free drill-membrane strain of paper B. Its live ingredient set SHALL be exactly:

1. 3D continuum kinematics with nodal thickness `a_i` and nodal director `V_n^i` — A **Eqs. (1)–(3)**, p. 405.
2. Director enrichment `u_b = ½ Σ a_i h_i (−V_2^i α_i + V_1^i β_i)` and geometry enrichment `x_b = ½ Σ a_i h_i V_n^i` — A **Eqs. (8a)/(8b)**, p. 406.
3. Covariant 3D strain `e_ij = ½(g_i·u_j + g_j·u_i)`, `r_1=r, r_2=s, r_3=t` — A **Eqs. (4)–(5)**.
4. Decomposition `e_ij = e_ij^m + t·e_ij^b1 + t²·e_ij^b2` (A **Eq. (7a)**) with the membrane term `e_ij^m` (A **Eq. (7b)**) and the **displacement-based** bending terms (A **Eqs. (7c)/(7d)**), including their `∂x_b` parts.
5. Assumed membrane field — the only modified term — sampled at the five tying points A(0,+1), B(0,−1), C(+1,0), D(−1,0), E(0,0) (A Fig. 4), with the five geometric coefficients of A **Eqs. (23)–(25)** and the field of A **Eqs. (17)–(27)**, including **Eq. (26)** and the efficient form **Eqs. (27a–c)**.
6. Transverse shear as the MITC4 assumed field of DB84 **Eq. (3)** as reproduced in A p. 405: `ẽ_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, `ẽ_st = ½(1+r)e_st^C + ½(1−r)e_st^D`.
7. Drill-membrane strain `e_ij^md`: B **Eq. (26)** (`e_ij = ẽ_ij^m + e_ij^md + t·e_ij^b1 + t²·e_ij^b2`) with the operator equations the 2025 element is built from: B **Eqs. (5), (10), (11a–b), (13b–c), (17a–b), (18), (19a–c), (21)**.
8. Integration: **2×2×2** Gauss over the element domain — A p. 410; and **2×2 over the element surfaces** for the 2025 six-DOF contribution — B §3.1, p. 13.

The element SHALL contain **no numerical factor**: no penalty term of any kind, no shear
correction factor, no selective reduced integration, no incompatible-mode (bubble)
rotation DOF, and none of the Winkler & Plakomytis ERC or `beta_w` ingredients. The
forbidden, non-inherited ingredients are exactly: (1) the uncited drilling penalty
`k_drill = 0.15·E·h²·drilling_scale`; (2) the ERC drilling treatment and `beta_w` warping
penalty; (3) selective reduced integration of the in-plane shear (the
`cm_normal`/centre-point split); (4) the 2-DOF `(1−ξ²)(1−η²)` rotation bubble; (5) the
shear correction factor `5/6`; (6) the membrane field of A **Eqs. (18)–(19)** *as the final
field* (it is admissible only as the exact flat-rectangle reduction of Eq. 27, see
Requirement 3); (7) flat projection with ABD resultants in place of the 3D continuum
kinematics; (8) hourglass scaffolding; (9) any other numerical factor.

A states (p. 410): *"the element formulation does not include any numerical factor"*.

#### Scenario: All live ingredients are locked to a paper-faithful reference

- GIVEN a test-local reference implementation, written from the ingredient list above (A Eqs. 1–27 with Eqs. 7c/7d displacement-based, DB84 Eq. 3 as in A p. 405, B Eq. 26, 2×2×2 Gauss, plane-stress material, every factor at its printed value)
- WHEN the production element stiffness `K` is computed for (i) a flat unit square, (ii) a flat distorted quad, (iii) a ruled warped quad, (iv) a doubly warped quad
- THEN `max|K_production − K_reference| ≤ 1e-10 · max|K_reference|` for every geometry, and the element's membrane/transverse-shear bending operator blocks match the reference blocks with the same bound.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_ke_lock_matches_2017_core_plus_2025_drill` — **NEW** `test_identity_ke_lock_matches_2017_core_plus_2025_drill` — encodes A Eqs. (1)–(27), (7c)/(7d), DB84 Eq. (3), B Eq. (26).

#### Scenario: The drilling stiffness comes only from the 2025 drill-membrane strain, not from a penalty

- GIVEN the production element and a variant reference in which the B Eq. (26) drill-membrane operator is replaced by the zero operator
- WHEN the two local stiffness matrices are computed on flat and warped geometry
- THEN the difference `K(operator) − K(operator := 0)` is non-zero on warped geometry (the operator is live) and vanishes identically when the operator is zeroed (so a test that passes by making the operator inert is caught), is exactly symmetric, is exactly zero on every **translational** row/column block, and is exactly zero on every rotation row/column block **other than the drill's own** whenever the element is flat (so that `V^D = e3` and the operator touches only slot `6·i + 5`; the flat case is what makes this clause exactly testable); on warped geometry the drill DOF's coupling to the nodal rotation components along `V^D` is required by paper B p. 8 and is therefore permitted, while every translational block remains exactly zero and the six rigid-body fields still have strain energy exactly zero (`|u_rbᵀ K u_rb| ≤ 1e-12 · λ_max(K) · ‖u_rb‖²`), so no `θ_z`-only stiffness term can exist.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_drill_stiffness_comes_only_from_eq26` — **NEW** `test_identity_drill_stiffness_comes_only_from_eq26` — encodes B Eq. (26) with the `θ^D = θ·V^D` coupling of B p. 8 (permitted on warped geometry, exactly zero on flat rotation blocks other than the drill's own); the rigid-body condition encodes A Eq. (3).
- Supporting static check (no penalty, no ERC, no `beta_w`, no SRI, no bubble in the code path that runs): `grep -nE 'k_drill|drilling_scale|compute_ke_local_erc|beta_w|hg_stiffness_factor|cm_normal' <elements that the MITC4+/D type dispatches to>` MUST return no match. The module paths are fixed by design; this check is supporting, not the normative evidence.

#### Scenario: Transverse shear uses the uncorrected plane-stress shear modulus

- GIVEN a flat rectangular element with isotropic material `(E, ν)` and thickness `h`
- WHEN the transverse-shear block of the local stiffness is compared with the closed-form `∫ B_γᵀ C_s B_γ dA` over the 2×2 rule with `C_s = G·I`, `G = E/(2(1+ν))`
- THEN they agree to `1e-10` relative, and the value obtained with the `5/6` shear correction factor is rejected by a margin of at least `1e-3` relative (the test asserts the two candidate values are separated, so it cannot pass vacuously).

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_uses_uncorrected_shear_modulus` — **NEW** — encodes A p. 410 (no numerical factor) and DB84 Eq. (3).

#### Scenario: The integration rule is 2×2×2 and is discriminated from a surface-only rule

- GIVEN a warped quad on which the 2×2×2 rule and the 2×2 surface-only rule produce measurably different stiffness contributions
- WHEN the production `K` is compared with both references
- THEN `K` matches the 2×2×2 reference to `1e-10` relative, and differs from the 2×2 surface-only reference by more than `1e-4` relative (asserted, so the discrimination is real), and the surface-only summation is permitted to apply only to the B Eq. (26) drill contribution per B §3.1 p. 13.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` — **NEW** — encodes A p. 410 and B §3.1 p. 13.

#### Scenario: Bending is displacement-based

- GIVEN a flat and a warped element and a nodal field with non-zero director rotations `α_i, β_i`
- WHEN the bending strain `(e_ij^b1, e_ij^b2)` produced by the element is compared with a test-local evaluation of A Eqs. (7c)/(7d) including the `∂x_b·∂u_m` contribution of A Eq. (8a)
- THEN they agree to `1e-10` relative at every Gauss point, and the rotation contribution is carried by the per-node enrichment of A Eq. (8b) rather than by any internal bubble DOF (the local matrix is 24×24 with no condensed internal DOF, see Requirement 2).

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_bending_operator_matches_eq7c_eq7d` — **NEW** — encodes A Eqs. (7c)/(7d) and (8a)/(8b).

### Requirement: Kinematics and DOF

The element's core kinematics SHALL be the five-DOF continuum formulation of A Eqs.
(1)–(3): three translations `u_i` plus two director rotations `α_i` (about `V_1^i`) and
`β_i` (about `V_2^i`), with the through-thickness coordinate `t ∈ [−1, 1]` and the
per-node thickness `a_i` and director `V_n^i`. The **sixth (drilling) DOF** SHALL be
supplied by the B Eq. (26) drill-membrane strain `e_ij^md`, not by a penalty. The
element's local matrices SHALL remain **24×24** (6 DOF/node, drilling at slot `6·i + 5`),
so the existing assembly stride, PyO3 batch shapes and Python family table are unchanged.

#### Scenario: The displacement field matches the paper's kinematics

- GIVEN a single element with nodal coordinates `x_i`, nodal thicknesses `a_i`, nodal directors `V_n^i` and nodal DOF `(u_i, α_i, β_i)`
- WHEN the element interpolates its displacement at sampled `(r, s, t)` points
- THEN it equals A Eq. (3) (`Σ h_i u_i + (t/2) Σ a_i h_i (−V_2^i α_i + V_1^i β_i)`) with `h_i` bilinear and `t` scaled by `t/2`, and the geometry equals A Eq. (1), to `1e-12` relative, including at `t = ±1`.

  Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_displacement_field_matches_eq1_to_eq3` — **NEW** — encodes A Eqs. (1)–(3).

#### Scenario: The sixth DOF is the 2025 drill-membrane strain

- GIVEN a nodal field that is a pure rigid rotation of the element about its own `V_n` (constant `θ_z ≠ 0`, with the accompanying rigid translations and `α_i = β_i = 0`)
- WHEN the element's strain energy is evaluated
- THEN the energy is exactly zero (`≤ 1e-12 · λ_max · ‖u‖²`), proving no penalty couples to `θ_z`; and on a warped element a nodal `θ_z` pattern that violates the constant-strain state produces non-zero energy through `e_ij^md` only (the B Eq. (26) operator), so the DOF is neither inert nor penalized.

  Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_drill_dof_is_theta_z_through_eq26_operator` — **NEW** — encodes B Eq. (26) and A Eq. (3).

#### Scenario: Local matrices stay 24×24

- GIVEN the element type registered for the MITC4+/D
- WHEN its local and global stiffness, mass and tangent matrices are produced
- THEN each is exactly 24×24 (576 entries) and the internal force vector is exactly 24 long, with the drilling DOF at global slot `6·i + 5` for node `i`.

  Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_local_matrices_are_24x24` — **NEW**. Cross-check at the binding level: `python -m pytest "tests/test_rust_composite.py::TestMITC4BatchSanity::test_ke_shape_and_symmetry" -q` — **EXISTING**.

### Requirement: The five geometric coefficients

The element SHALL compute the characteristic vectors of A Eqs. (9)–(13) — `x_r = ¼Σ ξ_i x_i`, `x_s = ¼Σ η_i x_i`, `x_d = ¼Σ ξ_i η_i x_i`, the normal `n` of A **Eq. (10)**, and the dual basis `m^r, m^s` of A **Eq. (11)** — and the five coefficients `a_A = c_r(c_r−1)/(2d)`, `a_B = c_r(c_r+1)/(2d)`, `a_C = c_s(c_s−1)/(2d)`, `a_D = c_s(c_s+1)/(2d)`, `a_E = +2c_r c_s/d` (A **Eqs. 27a–c**) with `c_r = x_d·m^r`, `c_s = x_d·m^s`, `d = c_r² + c_s² − 1` (A **Eq. 23–25**). `a_E` SHALL be **positive**, as printed (the negative form implemented by the deleted code is wrong).

#### Scenario: A flat rectangle produces zero distortion and zero coefficients

- GIVEN a flat rectangular element
- WHEN `x_d`, `c_r`, `c_s`, `d` and `a_A..a_E` are computed
- THEN `|x_d| ≤ 1e-14 · max(|x_r|, |x_s|)`, `|c_r| ≤ 1e-14`, `|c_s| ≤ 1e-14`, `d = −1` to `1e-14`, all five coefficients are zero to `1e-14`, and the assumed membrane field of A Eq. (27) reduces term-by-term to A Eq. (18) to `1e-17` absolute.

  Verification: `cd crates && cargo test -p aeroelast-core test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18` — **NEW** — encodes A Eqs. (18), (27a–c).

#### Scenario: The dual-basis identities hold

- GIVEN a flat distorted quad, a ruled warped quad and a doubly warped quad
- WHEN `x_r, x_s, n, m^r, m^s` are computed
- THEN `m^r·x_r = 1`, `m^s·x_s = 1`, `m^r·x_s = 0`, `m^s·x_r = 0` (each to `1e-12` absolute) and `m^r·n = 0`, `m^s·n = 0` (each to `1e-12 · |m^r|·|n|`), for every geometry.

  Verification: `cd crates && cargo test -p aeroelast-core test_geometry_dual_basis_identities_eq11` — **NEW** — encodes A Eq. (11) with A Eqs. (9)–(10).

#### Scenario: a_E has the printed positive sign

- GIVEN a distorted geometry with `c_r ≠ 0` and `c_s ≠ 0` of the same sign
- WHEN `a_E = 2c_r c_s/d` is evaluated
- THEN `a_E > 0` and the value equals `2c_r c_s/(c_r² + c_s² − 1)` to `1e-14`; the test also asserts that the negated value differs from the computed one by more than `1e-6` relative, so the sign is discriminated rather than merely tolerated.

  Verification: `cd crates && cargo test -p aeroelast-core test_geometry_a_E_is_positive_eq27c` — **NEW** — encodes A Eq. (27c) with the p. 410 sign.

### Requirement: Tier 1a — isotropy of the 2017 core

The element SHALL pass the isotropy test of A Section 4 (p. 410): the single-element stiffness SHALL be invariant under the element's spatial orientation and under the node-numbering sequence.

#### Scenario: Orientation and node-numbering invariance

- GIVEN one flat square, one flat distorted quad and one warped quad, each built (i) unrotated and (ii) rotated by a rigid rotation `R` about an axis not aligned with the element normal (at least four distinct orientations, including one with a `π/2` rotation), and each rebuilt with node sequences `(1,2,3,4)`, `(2,3,4,1)`, `(3,4,1,2)`, `(4,1,2,3)` and `(1,4,3,2)`
- WHEN the 24×24 stiffness is formed in global coordinates, so that the rotated element's rotations act on the co-rotated DOF, and the symmetrised eigenvalues are sorted
- THEN `|λ_i(rotated) − λ_i(reference)| ≤ 1e-10 · λ_max(reference)` for every `i`, the node-sequence variants reproduce the exactly-permuted entries of the reference matrix to `1e-12 · max|K|`, and for a physically co-rotated field `u` the strain energy `uᵀKu` is invariant to `1e-10` relative.
- GIVEN the paper publishes no numeric tolerance for this test, the `1e-10`-relative thresholds above are the round-off acceptance criteria fixed by this spec; a real orientation dependence is a formulation defect and MUST NOT be absorbed by relaxing them.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_isotropy_element_orientation_and_node_sequence_invariant` — **NEW** — encodes A Section 4 p. 410 (isotropy), A Eq. (3) (DOF co-rotation).

### Requirement: Tier 1a — exactly six zero-energy modes on a single unsupported element

The stiffness of a **single unsupported element** SHALL have exactly six zero eigenvalues, corresponding to the six rigid-body modes, with no spurious (hourglass or drilling-mechanism) zero mode.

#### Scenario: Exactly six zero eigenvalues, flat and warped

- GIVEN a single element with no boundary conditions — (i) flat rectangle, (ii) flat distorted quad, (iii) warped quad — and the six physical rigid-body fields `u = t + ω×(x − x_c)` with `θ = ω`, expressed in the element's 6-DOF layout (for the 2017 core with `α = β = 0`, since `ω × V_n = 0` there)
- WHEN the symmetrised stiffness is eigendecomposed with `λ_max = max|λ_i|` as the scale
- THEN the count of eigenvalues with `|λ_i| ≤ 1e-10 · λ_max` is **exactly 6**, and `‖K u_rb‖∞ ≤ 1e-10 · λ_max · ‖u_rb‖∞` for each of the six fields, and no eigenvalue lies in the ambiguous band between the six zero modes and the first structural mode (the test asserts the six smallest `|λ|` are separated from the seventh by at least `1e-6 · λ_max`, so the count cannot be an artefact of a loose threshold).

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_zero_energy_modes_single_unsupported_element_exactly_six` — **NEW** — encodes A Section 4 p. 410 (zero energy mode test) with A Eq. (3).

### Requirement: Tier 1a — membrane patch test

On the distorted mesh of A **Fig. 5** (10×10 square, interior nodes at (2,2), (4,7), (8,7), (8,3)) — or, where a boundary-node coordinate cannot be read reliably from the figure, on the **documented equivalent distorted patch** of the fixture policy (proposal §2.1) with the deviation recorded — the element SHALL reproduce the correct **constant** membrane stress field anywhere in the mesh under boundary nodal forces from that constant stress state and the minimum number of constraints needed to prevent rigid-body motion.

#### Scenario: Constant membrane stress on the Fig. 5 distorted mesh

- GIVEN the Fig. 5 distorted patch (fixture: the **figure-read mesh** of A Fig. 5, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — never an invented mesh, proposal §2.1), with the minimum constraints against rigid-body motion and boundary nodal forces equivalent to a prescribed constant membrane stress state `σ = [σ_xx, σ_yy, τ_xy]`
- WHEN the patch is solved and the element stress is recovered at every Gauss point of every element
- THEN the recovered membrane stress equals the prescribed constant state to `1e-8` relative (with an absolute floor of `1e-10 · ‖σ‖`), for each of the three independent constant states (`σ_xx` alone, `σ_yy` alone, `τ_xy` alone), and the recovered field is constant, i.e. the spread across Gauss points is `≤ 1e-8·‖σ‖`.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_membrane_patch_constant_stress_fig5_mesh` — **NEW** — encodes A Eqs. (17)–(27) with A Section 4 p. 410 (patch tests).

#### Scenario: The paper's own patch-test condition Eq. (22)

- GIVEN a flat element (any shape, `x_d·n = 0`) and arbitrary nodal displacements
- WHEN the assumed bilinear membrane shear strain `ẽ_rs^m|bil` of A Eq. (22)/(26) is compared with the displacement-based `e_rs^m|bil`
- THEN they are equal to `1e-14` absolute for all `(r, s)` in the element, and on a warped element the same comparison is non-zero beyond `1e-6` relative (the condition is asserted to be non-vacuous).

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_membrane_eq22_flat_tying_condition` — **NEW** — encodes A Eq. (22).

### Requirement: Tier 1a — bending patch test

On the same Fig. 5 distorted mesh — or the documented equivalent distorted patch under the fixture policy (proposal §2.1) — the element SHALL reproduce the correct **constant curvature** (constant bending moment) state anywhere in the mesh, with the minimum constraints against rigid-body motion and boundary nodal forces from that constant bending state.

#### Scenario: Constant curvature on the Fig. 5 distorted mesh

- GIVEN the Fig. 5 patch (fixture: the **figure-read mesh** of A Fig. 5, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — proposal §2.1) with the minimum constraints and boundary loading equivalent to a prescribed constant curvature state
- WHEN the patch is solved and the curvature / bending moment is recovered at every Gauss point of every element
- THEN the recovered curvature equals the prescribed constant value to `1e-8` relative, per independent constant state, and the spread across Gauss points is `≤ 1e-8` relative to the state's magnitude.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_bending_patch_constant_curvature_fig5_mesh` — **NEW** — encodes A Eqs. (7c)/(7d) with A Section 4 p. 410.

### Requirement: Tier 1a — shearing patch test

On the same Fig. 5 distorted mesh — or the documented equivalent distorted patch under the fixture policy (proposal §2.1) — the element SHALL reproduce the correct **constant transverse shear stress** state anywhere in the mesh, with the minimum constraints against rigid-body motion and boundary nodal forces from that constant shear state.

#### Scenario: Constant shear stress on the Fig. 5 distorted mesh

- GIVEN the Fig. 5 patch (fixture: the **figure-read mesh** of A Fig. 5, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — proposal §2.1) with the minimum constraints and boundary loading equivalent to a prescribed constant transverse shear stress state
- WHEN the patch is solved and the transverse shear strain/stress is recovered at every Gauss point of every element
- THEN the recovered shear equals the prescribed constant value to `1e-8` relative, and the spread across Gauss points is `≤ 1e-8` relative, without any shear correction factor (Requirement 1, third scenario).

  Verification: `cd crates && cargo test -p aeroelast-core test_t1a_shearing_patch_constant_stress_fig5_mesh` — **NEW** — encodes DB84 Eq. (3) as reproduced in A p. 405, with A Section 4 p. 410.

### Requirement: Tier 1b — spatial isotropy of the 2025 six-DOF element

The element **with the drilling DOF active** SHALL pass the spatial isotropy test of B Section 3.1 (pp. 13–14), as paper B states: *"MITC/D and MITC4+/D elements pass the spatial isotropy test."*

#### Scenario: Isotropy with the drill DOF free

- GIVEN a flat and a warped single element, built in at least four global orientations, with the drilling DOF free at all nodes (no `θ_z` constraint)
- WHEN the 24×24 stiffness is formed in global coordinates and compared across orientations
- THEN the sorted eigenvalues agree to `1e-10 · λ_max`, and for a co-rotated field `u` that includes a non-zero drilling component the strain energy `uᵀKu` is invariant to `1e-10` relative, so the drill contribution does not break isotropy.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_spatial_isotropy` — **NEW** — encodes B Section 3.1 pp. 13–14 with B Eq. (26).

### Requirement: Tier 1b — zero-energy modes and rigid-body representation with the drill DOF

The 2025 element SHALL pass the zero-energy-mode test with the rigid-body modes properly represented: the count of zero eigenvalues of a single unsupported element SHALL remain **exactly 6**; the drilling DOF SHALL NOT add a seventh zero mode (a drilling mechanism) and SHALL NOT remove a rigid-body mode.

#### Scenario: Exactly six zero eigenvalues with the drill DOF present

- GIVEN a single unsupported element with all six DOF/node free — (i) flat, (ii) warped, (iii) flat distorted — and the six physical rigid-body fields including the `θ_z` component of the rigid rotation about `V_n`
- WHEN the symmetrised stiffness is eigendecomposed with `λ_max = max|λ_i|`
- THEN the count with `|λ_i| ≤ 1e-10 · λ_max` is **exactly 6**, each of the six rigid-body fields has energy `≤ 1e-12 · λ_max` under the normalisation `‖u‖ = 1`, the six smallest `|λ_i|` are separated from the seventh by at least `1e-6 · λ_max`, and the drilling DOF carries non-zero stiffness in at least one non-rigid mode.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_zero_energy_modes_exactly_six_with_drill_dof` — **NEW** — encodes B Section 3.1 pp. 13–14 with B Eq. (26).

### Requirement: Tier 1b — strong-form patch tests

The element SHALL pass the **strong-form** patch tests of B Section 3.1: extension, bending and shearing patches, with the **minimum boundary conditions** to prevent rigid-body motion, such that the calculations give the analytical solutions of **constant AND zero stresses** throughout the patch. The drill DOF `θ_z` SHALL be free at the element nodes **except at the corner node B**, and setting `θ_z = 0` at corner node C SHALL NOT affect the results. Each patch uses the **figure-read mesh** of B Fig. 7 or, where a boundary coordinate cannot be read reliably, the **documented equivalent distorted patch** of the fixture policy (proposal §2.1) with the deviation recorded; an extension/bending/shearing case that can be established as neither remains untestable (Evidence gaps G3).

Interpretation used here (recorded as an assumption, see Evidence gaps): for each test the loaded direction carries the constant analytical stress and the other stress components are zero.

#### Scenario: Extension patch — constant and zero stresses

- GIVEN the B Fig. 7 extension patch (fixture: the **figure-read mesh** of B Fig. 7, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — proposal §2.1) with the minimum boundary conditions and `θ_z` free except at corner node B
- WHEN the patch is solved
- THEN the axial stress equals the analytical constant value to `1e-8` relative in every element, and every stress component that the analytical solution gives as zero has magnitude `≤ 1e-10 · |σ_axial|` throughout the patch.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_strong_patch_extension_constant_and_zero_stress` — **NEW** — encodes B Section 3.1 pp. 13–14, Fig. 7, with B Eq. (26).

#### Scenario: Bending patch — constant and zero stresses

- GIVEN the B Fig. 7 bending patch (fixture: the **figure-read mesh** of B Fig. 7, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — proposal §2.1) with the minimum boundary conditions and `θ_z` free except at corner node B
- WHEN the patch is solved
- THEN the bending-produced stress distribution equals the analytical (constant-moment) solution to `1e-8` relative, and the analytically-zero components are `≤ 1e-10` of the state's magnitude throughout the patch.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_strong_patch_bending_constant_and_zero_stress` — **NEW** — encodes B Section 3.1 pp. 13–14, Fig. 7, with B Eq. (26).

#### Scenario: Shearing patch — constant and zero stresses

- GIVEN the B Fig. 7 shearing patch (fixture: the **figure-read mesh** of B Fig. 7, or the **documented equivalent distorted patch** where a boundary coordinate cannot be read reliably — proposal §2.1) with the minimum boundary conditions and `θ_z` free except at corner node B
- WHEN the patch is solved
- THEN the transverse shear stress equals the analytical constant value to `1e-8` relative, and the analytically-zero components are `≤ 1e-10` of that magnitude throughout the patch.

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_strong_patch_shearing_constant_and_zero_stress` — **NEW** — encodes B Section 3.1 pp. 13–14, Fig. 7, with B Eq. (26).

#### Scenario: The θ_z boundary condition is free except at corner B, and corner-C fixing is immaterial

- GIVEN the extension patch (and repeated for bending and shearing)
- WHEN three variants are solved: (`a`) `θ_z` free at every node except corner node B; (`b`) `θ_z = 0` additionally imposed at corner node C; (`c`) `θ_z` constrained at every node
- THEN variants `(a)` and `(b)` produce nodal displacements and stresses that agree to `1e-10` relative (so the `θ_z = 0` condition at C does not affect the results), and variant `(c)` is asserted to differ from `(a)` beyond `1e-8` relative on at least one warped patch (so the `θ_z`-free condition is asserted to matter where it should, and the pass cannot come from an over-constrained model).

  Verification: `cd crates && cargo test -p aeroelast-core test_t1b_drill_theta_z_free_except_corner_b` — **NEW** — encodes B Fig. 7 and B Section 3.1 pp. 13–14.

### Requirement: Tier 2 — existing theoretical tests keep passing

After the element change the repository's existing theoretical suite SHALL keep passing, with **no new failures**, and the two pre-existing failures (`test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple`, `test_shell_convergence.py::test_in_bending_convergence`) SHALL be reported as pre-existing and SHALL be **excluded from the pass condition**. Suite baselines that MUST NOT regress: `cd crates && cargo test -p aeroelast-core` → **120 passed / 0 failed**; `python -m pytest -m "not slow" -q` → **345 passed / 2 failed / 2 skipped**.

The tests retired with the old element (T2J of the proposal §2: `test_drill_midside_derivatives_match_paper_eq11`, `test_drill_membrane_operator_*` (5), `test_b_erc_*` (5), `test_erc_covariant_to_local_*` (2), `test_ke_local_erc_flat_rigid_body_and_symmetry`, `test_enhanced_drill_stiffness`, `test_drill_warping_moment_flat_element`, `test_ke_local_eigenvalues_nonsymmetric`) and `test_quad_elements.py` (plane QUAD, not a shell element) are explicitly **not** part of the pass condition.

#### Scenario: Patch, rigid-body and zero-mode invariants in the 6-DOF layout

- GIVEN the element change applied
- WHEN `cd crates && cargo test -p aeroelast-core` runs
- THEN these names (T2A/T2B, `mitc4.rs` module `tests`) pass: `test_ke_global_leaves_all_six_rigid_body_modes_free`, `test_ke_global_has_exactly_six_zero_modes`, `test_membrane_patch_reproduces_constant_strain_at_every_gauss_point`, `test_bending_patch_reproduces_constant_curvature_at_every_gauss_point`, `test_ke_global_is_symmetric`, `test_ke_global_is_positive_semidefinite`, `test_ke_local_flat_plate_parity`, `test_kt_fint_directional_derivative`, `test_kt_fint_directional_derivative_rotations`, `test_kt_fint_directional_derivative_with_drill_dofs`, `test_fint_linear_nonlinear_parity`, `test_kt_zero_matches_ke`; where a name is layout-bound it is satisfied by its mechanical 6-DOF-layout equivalent (proposal T2I), not by a tolerance change.
- Verification: `cd crates && cargo test -p aeroelast-core` (**EXISTING** names).

#### Scenario: Assembly, tangent and mass invariants

- GIVEN the rebuilt extension (`python -m maturin develop --release`)
- WHEN `python -m pytest -m "not slow" -q` runs
- THEN `tests/test_rust_assembler.py` (all `TestTangentStiffness`, `TestInternalForces`, `TestNewtonRaphsonConsistency` tests, including `test_kt_at_zero_equals_k`, `test_kt_symmetric`, `test_fint_zero_at_zero`, `test_fint_linear_equals_ku`, `test_nr_consistency`), `tests/test_mass_matrix_validation.py` (all 15, including `TestExactConsistentMassCoefficients`, `TestConsistentMassTotal`, `TestElementMassVsTotalMass`) and `tests/test_rust_composite.py::TestMITC4BatchSanity` (`test_ke_shape_and_symmetry`, `test_ke_positive_semidefinite`, `test_me_shape_and_symmetry`, `test_me_positive_semidefinite`) pass.
- Verification: `python -m pytest "tests/test_rust_assembler.py" "tests/test_mass_matrix_validation.py" "tests/test_rust_composite.py::TestMITC4BatchSanity" -q` (**EXISTING**).

#### Scenario: Analytical shell validation stays green

- GIVEN the rebuilt extension
- WHEN the analytical groups run
- THEN these pass: `tests/test_material_suite.py::TestIsotropicAnalytical` (`test_mitc4_out_of_plane`, `test_mitc4_in_plane_lateral`, `test_mitc4_axial_stiffness`), `tests/test_material_suite.py::TestStiffnessProperties` (`test_k_positive_definite`, `test_k_symmetric`), `tests/test_material_suite.py::TestOrthotropicSinglePly`, `::TestSymmetricLaminates`, `::TestAsymmetricLaminates`, `tests/test_shell_analytical_validation.py::TestCantileverBeam` (incl. `test_tip_deflection[4x2]`, `[8x4]`, `[16x8]`), `::TestSimplySupportedBeam::test_center_deflection`, `::TestBeamBending::test_bending_convergence`, `::TestMembraneStretching::test_axial_extension`, `::TestShearLocking::test_thin_plate_convergence`, `tests/test_shell_comprehensive.py` (fx/fy/fz, ratio, modal), `tests/test_shell_validation_fixed.py` (fx/fy/fz/ratio/axial convergence/modal), `tests/test_large_rotation_benchmarks.py::test_linear_tip_deflection_euler_bernoulli`, `::test_cantilever_large_rotation_half_circle`, `::test_equilibrium_path`, `::test_simo_vu_quoc_rollup_360`, `tests/test_rust_modal.py::TestSimplySupportedPlate::test_analytical_convergence`, `tests/test_composite_b_coupling.py`, and the CCX-parity files `tests/test_isotropic_shell_parity.py::TestIsotropicShellParity::test_transverse_tip_displacement`, `tests/test_beam_shell_4cases_parity.py`, `tests/test_composite_beam_parity.py`, `tests/test_orthotropic_shell_parity.py` (the CCX files skip when `ccx` is not on `PATH`).
- Verification: `python -m pytest "tests/test_material_suite.py" "tests/test_shell_analytical_validation.py" "tests/test_shell_comprehensive.py" "tests/test_shell_validation_fixed.py" "tests/test_large_rotation_benchmarks.py" "tests/test_rust_modal.py" "tests/test_composite_b_coupling.py" "tests/test_isotropic_shell_parity.py" "tests/test_beam_shell_4cases_parity.py" "tests/test_composite_beam_parity.py" "tests/test_orthotropic_shell_parity.py" -q` (**EXISTING**).

#### Scenario: Published benchmark columns stay inside their published tolerance

- GIVEN the rebuilt extension
- WHEN `tests/test_ko2017_performance.py` runs (31 cases, `rtol=0.05`)
- THEN all nine test functions pass: `test_3_1_square_plate_tables_2_to_5`, `test_3_2_circular_plate_tables_6_to_7`, `test_3_3_pinched_cylinder_tables_8_to_9`, `test_3_4_scordelis_lo_tables_10_to_11`, `test_3_5_twisted_beam_tables_12_to_13`, `test_3_6_hook_table_14_minimal_fix`, `test_3_7_hemisphere_cutout_tables_15_to_16`, `test_3_8_full_hemisphere_table_17`, `test_3_9_hyperbolic_paraboloid_tables_18_to_19`; a value that moves but stays inside the published tolerance is a legitimate outcome (proposal risk §6.5) and MUST NOT be reported as a regression; a value outside tolerance is a formulation bug.
- Verification: `python -m pytest "tests/test_ko2017_performance.py" -q` (**EXISTING**).

#### Scenario: Convergence suite, with the pre-existing failure excluded

- GIVEN the rebuilt extension
- WHEN the convergence module runs
- THEN `tests/test_shell_convergence.py::test_composite_laminate_gap_mesh_study` passes and no test other than the two excluded pre-existing failures fails.
- Verification: `python -m pytest "tests/test_shell_convergence.py" -q -k "not test_in_bending_convergence"`; the excluded node ids are `tests/test_shell_convergence.py::test_in_bending_convergence` and `tests/test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple`. If either starts passing, it is recorded as a bonus, never as an acceptance criterion (proposal §4 non-goal 4).

### Requirement: Traceability of every live ingredient

The element's documentation SHALL name, for each live ingredient in Requirement 1, the paper, section and equation that defines it (A's 2017 core, DB84's shear field as reproduced in A, and B's 2025 drill-membrane strain), and SHALL NOT document any ingredient as sourced when it is not implemented, and SHALL NOT omit an ingredient that is implemented. The documented ingredient table is the human-readable form of Requirement 1 and MUST stay in agreement with the code.

#### Scenario: Every live ingredient has a resolvable citation

- GIVEN the element's documentation table (in `docs/formulations/shell-elements.md` §2 and the element's module header, per the proposal §3.2)
- WHEN the traceability test parses the table
- THEN every ingredient key of Requirement 1 is present with a non-empty paper/section/equation citation, and every cited equation number is present in the corresponding local artifact (`docs/formulations/mitc4plus-2017-extract.md` for A, `.sources/papers/1-s2.0-S0045794924003511-main.pdf` for B) or is explicitly recorded as pending a vision re-read (see Evidence gaps G4).

  Verification: `python -m pytest "tests/test_mitc4plusd_traceability.py::test_traceability_identity_table_names_source_for_every_ingredient" -q` — **NEW**.

#### Scenario: Nothing is documented as sourced that is not implemented

- GIVEN the documentation table and the element source that runs
- WHEN the test cross-checks the forbidden ingredient list of Requirement 1 against both
- THEN no forbidden ingredient (`k_drill`, ERC, `beta_w`, selective reduced integration, `(1−ξ²)(1−η²)` bubble, shear correction factor `5/6`, hourglass scaffolding) appears as an implemented ingredient in the table, and the source check of Requirement 1 (second scenario) returns no match; any remaining mention MUST be inside an explicit "not part of this element" list.

  Verification: `python -m pytest "tests/test_mitc4plusd_traceability.py::test_traceability_no_forbidden_ingredient_claimed_or_present" -q` — **NEW**.

#### Scenario: Equation citations resolve rather than being decorative

- GIVEN the same table
- WHEN each citation is resolved against the local artifacts and, for each, the tests that encode it are looked up by the naming convention in this spec (`test_identity_*`, `test_kinematics_*`, `test_geometry_*`, `test_t1a_*`, `test_t1b_*`)
- THEN every cited ingredient is referenced by at least one existing test name, so no row is unsupported by evidence and no test cites an unsourced ingredient.

  Verification: `python -m pytest "tests/test_mitc4plusd_traceability.py::test_traceability_equation_citations_resolve_to_extract_or_paper" -q` — **NEW**.

### Requirement: Retirement of the superseded hybrid, gated on Tier 1 and Tier 2

After — and only after — the MITC4+/D passes Tier 1 (T1.1–T1.6) **and** Tier 2 (the pass condition of the Requirement "Tier 2 — existing theoretical tests keep passing"), the superseded hybrid SHALL be removed from `crates/aeroelast-core/src/elements/mitc4.rs` (proposal §4, §7). Until that gate passes, the superseded code SHALL remain in place and the old element SHALL remain the reference. The removal SHALL delete: the Winkler & Plakomytis ERC drilling treatment, the `beta_w` warping penalty, the selective reduced integration, the 2-DOF `(1−ξ²)(1−η²)` rotation bubble, the A **Eqs. (18)–(19)** membrane as the live field, the shear-correction-factor usage in the element, the uncited `k_drill = 0.15·E·h²·drilling_scale`, and the dead code `b_m_standard`, `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`. The removal SHALL **NOT** delete `b_md_mitc4_plus` or `drill_midside_shape_derivatives`, which are revived as the 2025 drill ingredient (Requirement "Element identity and formulation", item 7). The removal SHALL **NOT** touch the laminate/composite material modules (Requirement "Preserved invariant — laminate/composite material models").

#### Scenario: The removal happens only after Tier 1 and Tier 2 pass

- GIVEN the superseded hybrid still in the tree and the new element under construction
- WHEN the retirement is scheduled
- THEN the superseded code is removed only after `cd crates && cargo test -p aeroelast-core test_t1a_`, `cd crates && cargo test -p aeroelast-core test_t1b_` and the Tier-2 pass condition (the Requirement "Tier 2 — existing theoretical tests keep passing") are all green; if any is not green the removal MUST NOT be applied and the superseded code MUST remain as the reference.

  Verification: process gate — the `test_retirement_*` tests MUST NOT be present-and-green while any Tier-1 or Tier-2 acceptance test is failing; asserted by running the gate commands (`cd crates && cargo test -p aeroelast-core` and the Python Tier-2 command) before the retirement tests are accepted.

#### Scenario: Every named deviation surface is gone from the live element path

- GIVEN the retirement applied
- WHEN the live element path is inspected and its tests run
- THEN none of `compute_ke_local_erc`, `beta_w`, the selective reduced integration (`cm_normal`/centre-point split), the `(1−ξ²)(1−η²)` bubble, the Eqs. (18)–(19) live membrane path, the element-side shear-correction-factor usage, `k_drill`/`drilling_scale`, `b_m_standard`, `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`, `Mat26` or `Vec26` remains in the code the element dispatches to, and the supporting static check `grep -nE 'k_drill|drilling_scale|compute_ke_local_erc|beta_w|hg_stiffness_factor|cm_normal|b_m_standard|green_lagrange_strain|compute_b_l|compute_membrane_stress|Mat26|Vec26' <elements that the MITC4+/D type dispatches to>` returns no match.

  Verification: `cd crates && cargo test -p aeroelast-core test_retirement_removes_hybrid_deviation_surfaces` — **NEW** — plus the supporting `grep` above (supporting, not the normative evidence).

#### Scenario: The revived 2025 drill operator is retained

- GIVEN the retirement applied
- WHEN the element source and the drill-membrane stiffness path are inspected
- THEN `b_md_mitc4_plus` and `drill_midside_shape_derivatives` are present and exercised by the stiffness, internal-force and tangent paths (Requirement "Element identity and formulation", second scenario), and the B **Eq. (26)** drill contribution is non-inert on warped geometry.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_drill_stiffness_comes_only_from_eq26` — **NEW** — plus the supporting static check that both symbols still exist.

#### Scenario: The PyO3 surface is unchanged by the removal

- GIVEN the retirement applied and the extension rebuilt (`python -m maturin develop --release`)
- WHEN the PyO3 batch kernels and the Python family table are inspected
- THEN `crates/aeroelast-py/src/elements.rs` still declares `[f64; 576]`/`[f64; 24]`, no `#[pyfunction]`/`#[pyclass]` symbol is removed or renamed, and `src/aeroelast/core/assembler.py`'s `_FAMILY_PROPERTIES[SHELL] = (6, 3)` is unchanged.

  Verification: `python -m pytest "tests/test_rust_composite.py::TestMITC4BatchSanity" -q` (**EXISTING**, after the rebuild) plus the static checks on the two files above (supporting).

### Requirement: Preserved invariant — laminate/composite material models

The laminate and composite material models SHALL be preserved by this change (hard user constraint, verbatim: *"si vas a sustituir todo el codigo viejo, recuerda que hay que mantener los modelos de material laminado"*). The Rust modules `crates/aeroelast-core/src/materials/{laminate,orthotropic,composite,failure,isotropic}.rs` and the Python modules `src/aeroelast/core/laminate.py` and `src/aeroelast/constitutive/failure.py` SHALL keep their public behaviour and their numerical results; no laminate capability may be dropped, weakened or re-scoped. The element consumes `ShellConstitutive { cm, cb_coupling, cb, cs, cm_raw }`, and that interface (the five fields and their meanings) SHALL be preserved. The retirement (Requirement "Retirement of the superseded hybrid, gated on Tier 1 and Tier 2") MUST NOT be used to justify touching, simplifying or removing any of these surfaces. **A laminate capability regression SHALL be treated as a change failure.**

#### Scenario: Laminate/composite tests keep passing with unchanged expectations

- GIVEN the element change and the retirement applied, and the extension rebuilt
- WHEN the laminate/composite groups run
- THEN `tests/test_material_suite.py::TestOrthotropicSinglePly`, `::TestSymmetricLaminates`, `::TestAsymmetricLaminates`, `tests/test_composite_b_coupling.py`, `tests/test_rust_composite.py` (including `TestMITC4BatchSanity`), `tests/test_orthotropic_shell_parity.py` and `tests/test_composite_beam_parity.py` all pass with their existing assertions unmodified (no expectation, tolerance or reference value changed by this change).

  Verification: `python -m pytest "tests/test_material_suite.py::TestOrthotropicSinglePly" "tests/test_material_suite.py::TestSymmetricLaminates" "tests/test_material_suite.py::TestAsymmetricLaminates" "tests/test_composite_b_coupling.py" "tests/test_rust_composite.py" "tests/test_orthotropic_shell_parity.py" "tests/test_composite_beam_parity.py" -q` (**EXISTING**), after `python -m maturin develop --release`.

#### Scenario: The ShellConstitutive interface is preserved and a capability regression fails the change

- GIVEN the element change and the retirement applied
- WHEN the Rust `ShellConstitutive` type and the laminate/composite public surfaces are inspected
- THEN `ShellConstitutive` still exposes exactly the five fields `cm`, `cb_coupling`, `cb`, `cs`, `cm_raw`; `Laminate` still exposes `to_shell_constitutive`, `to_shell_constitutive_with_offset`, `shear_correction_factor` and `abd_matrix_flat`; and the Rust materials unit tests pass — any removed or renamed public item, or any changed numerical result, is a change failure (not a regression to be tolerated).

  Verification: `cd crates && cargo test -p aeroelast-core materials::` (**EXISTING** unit tests) plus `python -m pytest "tests/test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged" -q` — **NEW** guard test.

### Requirement: Uncorrected transverse shear stiffness

The element's transverse shear stiffness SHALL be the **uncorrected** `G·h` — paper A p. 410: *"the element formulation does not include any numerical factor"* — with no shear correction factor `k`. The repository's `ShellConstitutive.cs` is in `k·G·h` form (its doc comment: `force-shear strain [k·G·h] form`) and `ShellConstitutive` does not carry `k` (while `Laminate::shear_correction_factor` does), so the element must obtain the uncorrected stiffness by one of the candidate mechanisms — (a) extend `ShellConstitutive` with the correction factor, (b) add an uncorrected transverse-shear accessor, or (c) divide by `k` at the element boundary. **Which mechanism is used is a design decision, not a spec decision**; this requirement fixes only the observable result and MUST be verifiable without knowing the mechanism. The chosen mechanism SHALL keep the laminate/composite invariant of the Requirement "Preserved invariant — laminate/composite material models" intact.

Related open point: `Laminate::to_shell_constitutive_with_offset` applies a reference-surface offset while the paper's kinematics assumes the reference surface is the mid-surface. The element path SHALL use the mid-surface form (`to_shell_constitutive`, i.e. `z_offset = 0`) unless and until a test proves the offset variant preserves the mid-surface kinematics; whether the offset variant may ever be used by the element is explicitly **deferred** (Evidence gap G7), and the element path SHALL NOT use the offset variant in this change.

#### Scenario: A residual k factor in the element's transverse shear stiffness is detected

- GIVEN a flat element with an isotropic or laminate constitutive whose transverse-shear block is `cs = k·G·h`, and the same constitutive with the laminate's `shear_correction_factor` changed to a different positive value `k' ≠ k` (so `cs` changes by the factor `k'/k` while `G·h` is unchanged)
- WHEN the element's transverse-shear contribution to the local stiffness is compared between the two constitutions
- THEN the two contributions agree to `1e-10` relative, proving the element carries no residual `k` (a residual factor would change the contribution by `k'/k` and fail the assertion by a margin `> 1e-3` relative, which the test asserts so it cannot pass vacuously); and the same contribution matches the closed-form uncorrected `∫ B_γᵀ (G·h) B_γ dA` to `1e-10` relative.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_invariant_to_shear_correction_factor` — **NEW** — encodes A p. 410 (no numerical factor) and DB84 Eq. (3).

#### Scenario: The element path uses the mid-surface constitutive, and the offset variant is deferred

- GIVEN the element's constitutive construction path
- WHEN the code path that feeds the element is inspected
- THEN it uses the mid-surface form `Laminate::to_shell_constitutive` (`z_offset = 0`) and does not call `to_shell_constitutive_with_offset`; whether the offset variant may ever be used by the element remains an explicit deferral recorded as Evidence gap G7, and the mid-surface assumption is asserted by a test that the element's constitutive equals `to_shell_constitutive()`.

  Verification: `cd crates && cargo test -p aeroelast-core test_identity_element_uses_midsurface_constitutive` — **NEW** — plus the supporting static check that the element path does not reference `to_shell_constitutive_with_offset` (supporting, not normative).

## Evidence gaps and blocking preconditions

These are stated so that no requirement above is read as vague. Pass conditions are fully
specified; the following fixtures or citations are one-time reads that must happen before
the corresponding test can be written.

**Fixture policy (replaces the former open-ended G1/G3 gaps, proposal §2.1).** The two
papers do not publish the patch-test boundary-node coordinates numerically. For every
patch fixture the change SHALL: (1) read paper A **Fig. 5** and paper B **Fig. 7** with
vision (`pdftoppm -png -r 300`, never `pdftotext`) and record in the change what was
actually read, **including any coordinate that cannot be determined**; (2) where a
boundary-node coordinate cannot be read reliably, use a **documented equivalent distorted
patch** and record the deviation explicitly — the requirement is then verified against that
documented equivalent, not against a silently invented mesh; (3) never invent a mesh and
never edit the reference to match the code.

What this makes **testable now**: the Tier 1a membrane/bending/shearing patch tests and the
Tier 1b extension/bending/shearing strong-form patch tests, each against its named fixture
(figure-read mesh, documented equivalent patch, or — for A Eq. (22) — analytic, no fixture).
What **remains genuinely untestable**: the paper's exact boundary-node coordinates where the
figure read cannot recover them (the papers do not print them); the accepted response is the
documented equivalent patch with the deviation recorded, so this is a **documented
deviation, not a failure**. The B Fig. 7 mesh/BC sets are the same case (G3).

| # | Affects | Gap | Precondition (and honest fallback) |
| --- | --- | --- | --- |
| G1 | Tier 1a membrane/bending/shearing patch requirements | Paper A **Fig. 5**'s boundary-node coordinates are not published numerically; the interior nodes (2,2), (4,7), (8,7), (8,3) of the 10×10 square are transcribed but the connectivity (inferred here: 3×3 elements, 16 nodes) is not. | **Resolved by policy** (§2.1 of the proposal): one vision read of A pp. 410–411 (`pdftoppm -png -r 300 .sources/papers/A_new_MITC4+_shell_element.pdf`), recording what is actually read **including any coordinate that cannot be determined**. Where a boundary-node coordinate cannot be read reliably, the tests use a **documented equivalent distorted patch** and record the deviation explicitly, and the requirement is verified against that documented equivalent. The exact paper boundary-node coordinates are the only genuinely untestable part and are recorded as such; inventing a mesh or editing the reference to match the code is forbidden. |
| G2 | Tier 1a isotropy, zero-energy | Paper A §4 publishes no numeric tolerance (only the pass/fail statement). | None: the `1e-10`-relative / `1e-6`-separation thresholds are fixed by this spec as round-off criteria and MUST NOT be relaxed (Requirement 1 of the proposal's non-goals). |
| G3 | Tier 1b strong-form patch requirements | Paper B Fig. 7's patch geometry, the exact minimum-BC sets and the "constant AND zero stresses" phrasing are held only as a summary (B pp. 13–14). The interpretation "loaded direction constant, other components zero" is an assumption of this spec. | **Resolved by policy** (§2.1 of the proposal): one vision read of B pp. 13–14 + Fig. 7, recording what is read **including any unreadable coordinate**; where a boundary coordinate cannot be read reliably, use a **documented equivalent distorted patch** and record the deviation. The interpretation "loaded direction constant, other components zero" remains a recorded assumption of this spec. An extension/bending/shearing case that can be established as neither a figure-read mesh nor a documented equivalent remains **untestable** and MUST be reported rather than implemented on an invented mesh. |
| G4 | Requirement 1, first scenario (ingredient lock) | The B **Eq. (26)** drill-membrane operator exists in the tree as dead code (`b_md_mitc4_plus`) that has never been exercised by a stiffness path, so its transcription is unverified (proposal risk §6.4). | A vision re-read of B Eqs. (5), (10), (11a–b), (13b–c), (17a–b), (18), (19a–c), (21), (26) to build the test-local reference. If the transcription is wrong, `test_identity_ke_lock_matches_2017_core_plus_2025_drill` fails — that is the intended signal, and it MUST NOT be silenced by changing the reference to match the code. |
| G5 | Tier 2 | Per-test current pass/fail values were not re-measured in this phase, and `docs/validation-matrix.md` is a stale `b2c62ff` snapshot. | None for testability: the pass condition is expressed as node-id selections plus the suite baselines, so it is verifiable by command (`cargo test` / `pytest -m "not slow" -q`) even though the snapshot is stale. |
| G6 | Requirement 7 | The 2025 paper's Eq. (26)-related citations in the element's doc comments were not verified in this phase. | Covered by G4; until then, the table MUST mark those rows "pending paper-B re-read" rather than claim a verified source. |
| G7 | Requirement "Uncorrected transverse shear stiffness", second scenario | `Laminate::to_shell_constitutive_with_offset` applies a reference-surface offset while the paper's kinematics assumes the reference surface is the mid-surface. | Explicitly **deferred** to design: this change restricts the element path to the mid-surface form `to_shell_constitutive` (`z_offset = 0`) and asserts it (`test_identity_element_uses_midsurface_constitutive`); whether the offset variant may ever be used by the element is not decided here. The deferral is testable as "the element path does not use the offset variant". |

## Out of scope (do not spec, do not implement under this change)

No shell-element trait; no 5-DOF stride migration and no assembly/PyO3 signature change
(the element stays 24×24); no work in `crates/aeroelast-solvers` (unbuildable here —
preCICE); no fixing the two pre-existing Python failures; no parameter fitting and no
tolerance loosening; no new numerical factor and no re-introduction of any forbidden
ingredient; no `pdftotext` reading of the papers; no twisted-beam benchmark runs at N ≥ 32.
