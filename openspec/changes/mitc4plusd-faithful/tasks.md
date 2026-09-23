# Tasks — the MITC4+/D element, faithfully implemented

Change: `mitc4plusd-faithful` · phase: **tasks** · store: openspec (Engram mirror:
`sdd/mitc4plusd-faithful/tasks`).

Backbone: the design's work-unit slicing (WU0…WU11) is preserved verbatim as the section
order; each section names the WU it implements, so every task is traceable to the design
and every WU keeps its own start, finish, verification and rollback boundary.

`strict_tdd: false` (see `openspec/config.yaml`), so the harness will not enforce
test-first. **Tests live in the same work unit as the behaviour they verify**, and every
task from WU2 onward adds its named test and observes it failing before the production
code is written.

## Review Workload Forecast

| Field | Value |
|-------|-------|
| Estimated changed lines | **≈ 3,760 total** (≈ 2,740 added + ≈ 1,020 removed); ≈ 1,720 net |
| 400-line budget risk | **High** (the canonical 400 threshold; the session budget is 700 and is *also* exceeded) |
| Chained PRs recommended | **No** — the session resolved delivery as **single-pr with an explicitly accepted `size:exception`** |
| Suggested split | Single PR. Fallback split points if the reviewer objects, recorded and **not chosen**: S1–S2 (additive element + Tier 1) → S3 (the flip) → S4 (the retirement) → WU11 (docs) |
| Delivery strategy | `single-pr` with accepted `size:exception` (per the session preflight) |
| Chain strategy | `size-exception` |

```text
Decision needed before apply: No
Chained PRs recommended: No
Chain strategy: size-exception
400-line budget risk: High
```

**Stated plainly.** The design's ledger forecasts ≈ 2,740 added / ≈ 1,020 removed
(≈ 1,720 net) changed lines. That exceeds the 700-line session review budget several
times over. This task list does **not** shrink the work to fit: no test, doc, comment or
citation is dropped, no WU is merged to hide lines, and no chain strategy is invented.
The accepted `single-pr (size:exception)` is what covers the total. What *is* held to the
budget is the **per-unit** bound (≤ ~350 changed lines per WU), so each WU is a
self-contained review slice with a green tree, its own verification, and its own
rollback.

## Read-only inputs (do not edit)

| Input | Why read-only |
| --- | --- |
| `.sources/papers/A_new_MITC4+_shell_element.pdf`, `.sources/papers/1-s2.0-S0045794924003511-main.pdf` | Reference papers. Vision only (`pdftoppm -png -r 300`), never `pdftotext`. |
| `openspec/changes/mitc4plusd-faithful/proposal.md` (rev 2), `.../specs/mitc4plusd-element/spec.md` (rev 4), `.../design.md` (rev 2) | Frozen upstream artifacts of this change. |
| `crates/aeroelast-core/src/assembly/topology.rs` | **Untouched by design** (ADR-4 option B: the element computes its own per-node directors, so `MeshTopology` gains no field). |
| `crates/aeroelast-core/src/elements/mitc4.rs` | Read-only until WU8 (tests only) and WU10 (deletions only); the live hybrid must stay byte-identical while the new element is built and tested. |
| `src/aeroelast/core/assembler.py` | Frozen: `_FAMILY_PROPERTIES[SHELL] = (6, 3)`. |
| `crates/aeroelast-core/src/materials/{orthotropic,composite,failure}.rs`, `src/aeroelast/core/laminate.py`, `src/aeroelast/constitutive/failure.py` | Preserved invariant (Requirement 15); no task may touch them. |
| `docs/formulations/mitc4plus-2017-extract.md`, `docs/formulations/mitc4plusd-2025-extract.md` | Read-only **except in WU0**, which is the only unit authorized to write them. |
| `tests/test_quad_elements.py`, `tests/test_shell_convergence.py`, `tests/test_rust_composite.py` | Tier-2 sources; expectations and tolerances unmodified. Only the sourcing corrections of task 10.6 touch any test file. |

Authorized edit roots for the whole change: `crates/aeroelast-core/src/elements/mitc4_plusd.rs`
(+ `…/mitc4_plusd/tests.rs`, `…/tests/fixtures.rs`, `…/tests/drill.rs`),
`crates/aeroelast-core/src/elements/mod.rs`, `crates/aeroelast-core/src/assembly/assembler.rs`,
`crates/aeroelast-core/src/materials/{mod,laminate}.rs`, `crates/aeroelast-py/src/{elements,assembler}.rs`,
`crates/aeroelast-core/src/elements/mitc4.rs` (tests at WU8, deletions at WU10), the two
`docs/formulations/*-extract.md` (WU0 only), `docs/formulations/shell-elements.md`,
`docs/validation-matrix.md`, `tests/test_mitc4plusd_traceability.py`,
`tests/test_laminate_invariant_guard.py`.

---

## 1. Precondition — paper vision re-read (WU0)

**WU0** · Touches: `docs/formulations/mitc4plusd-2025-extract.md`, `docs/formulations/mitc4plus-2017-extract.md` · Verification: the extract self-consistency checks below; downstream, the traceability test (WU11) · Rollback: revert the two Markdown files · **~140 changed lines.**

This unit is a **hard precondition** for the drill operator: the design's term-by-term
verification of B Eq. (18) (task 4.4) is only as good as the transcription this unit
produces. It settles open items 2, 3 and 4 of the design's table.

- [x] **1.1 Transcribe paper B Eqs. (1)–(16), the `h̃_m` mid-side functions and the `θ_z` interpolation.** Touches `docs/formulations/mitc4plusd-2025-extract.md`. Verification: the equations are present with page citations from a `pdftoppm -png -r 300` vision read; the downstream consumer is `test_traceability_equation_citations_resolve_to_extract_or_paper` (task 12.1). No automated oracle exists for Markdown math — the check is a recorded vision re-read, stated here rather than invented. Satisfies: Requirement 1 item 7; Requirement 13. <!-- sdd-owner: implementation -->
- [x] **1.2 Settle the B Eq. (16b)/(18) edge convention** (`θ_{i-1}^D − θ_i^D` as printed in the extract vs the dead code's `θ_{i+1} − θ_i`, `mitc4.rs:866-869`) and record the resolved convention plus the node-pair mapping with a page citation. Touches `docs/formulations/mitc4plusd-2025-extract.md`. Verification: the extract's Eq. (18) block states one unambiguous convention; the automated fallback detector is Oracle 2's sign-flip row in `test_identity_drill_operator_matches_eq18_term_by_term` (task 4.4). Satisfies: Requirement 1 item 7. <!-- sdd-owner: implementation -->
- [x] **1.3 Record F2 (paper A Eq. (21)'s missing leading term) as a note against the printed equation — do not "fix" the quote — and record the F1 Eq. (27a–c) → Eq. (18) reduction derivation.** Touches `docs/formulations/mitc4plus-2017-extract.md`. Verification: the printed Eq. (21) quote (`:347-349`) and the Eq. (27c) block (`:388-390`, carrying `(1 + a_E·rs) e_rs^m(E)`) are unchanged; the F1/F2 note is present; the reduction is re-checked by hand against `test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18` (task 3.2). Satisfies: Requirement 1 item 5; Requirement 3. <!-- sdd-owner: implementation -->
- [x] **1.4 Search paper A p. 405 for the paper's own covariant transverse-shear metric normalization.** If it is printed, replace the design's §2.2 definition with it; if not, record "not printed" and keep §2.2. Touches `docs/formulations/mitc4plus-2017-extract.md`. Verification: the extract records either the printed quote with a page citation or an explicit "not printed" line; the implementation in task 4.3 follows whichever is recorded. Satisfies: Requirement 1 item 6; Requirement 16. <!-- sdd-owner: implementation -->
- [x] **1.5 Add the `transcription-verified: <date> <page>` lines** to both extracts for every row the traceability test resolves, at minimum the Eq. (18) drill row. Touches both extract files. Verification: `grep -n "transcription-verified" docs/formulations/mitc4plusd-2025-extract.md docs/formulations/mitc4plus-2017-extract.md` returns the lines, and task 12.1's `test_traceability_equation_citations_resolve_to_extract_or_paper` requires the Eq. (18) row's line to be present. Satisfies: Requirement 13. <!-- sdd-owner: implementation -->

## 2. Fixtures (WU1)

**WU1** · Touches: `crates/aeroelast-core/src/elements/mitc4_plusd.rs` (inline `#[cfg(test)] mod tests`) · Verification: the fixture self-tests below · Rollback: delete the fixtures module · **~200 changed lines.** Fixture self-test names are not spec-fixed; the design calls them "fixture self-tests". **Recorded deviation:** the tasks name `…/mitc4_plusd/tests/fixtures.rs`, but the design fixes `…/mitc4_plusd.rs` and every existing element is a single file with an inline test module; per `openspec/config.yaml`'s apply guideline the fixtures live inline in `mitc4_plusd.rs` and the deviation is recorded in that file's header comment. The verification command is therefore `cd crates && cargo test -p aeroelast-core mitc4_plusd::tests::`.

- [x] **2.1 Star-patch fixture: 8 nodes / 5 elements, CCW connectivity in the element's `(r,s)` convention.** Touches `crates/aeroelast-core/src/elements/mitc4_plusd.rs` (`STAR_NODES`, `STAR_ELEMS`). Verification: `cd crates && cargo test -p aeroelast-core mitc4_plusd::tests::` — the self-test asserts 5 elements / 8 nodes, every signed area positive, the five areas summing to 100, no duplicate coordinates, and the central element's connectivity equal to the extract's. Satisfies: Requirement 6/7/8 (fixtures); Requirement 1 (figure-read mesh policy). <!-- sdd-owner: implementation -->
- [x] **2.2 The T1.3 minimum BC set (F-A-BC, derived) and the T1.6 BC sets (F-B-BC, figure-read) as code, each with its provenance recorded in a doc comment.** Touches `crates/aeroelast-core/src/elements/mitc4_plusd.rs`. Verification: self-test asserting each BC set constrains exactly the six rigid-body modes, and that F-B-BC reproduces `docs/formulations/mitc4plusd-2025-extract.md` §Fig. 7(b)(c)(d) including `θ_z` free except at corner B. F-A-BC is derived (paper A prints no BCs); design §4.1's literal set has rank 4 on the six rigid-body modes, so the fixture uses a corrected minimal set, recorded on `F_A_BC` and in the module header. Satisfies: Requirement 6/7/8; Requirement 11. <!-- sdd-owner: implementation -->
- [x] **2.3 Warped star patch F-W (documented as derived) plus the geometry self-tests.** Touches `crates/aeroelast-core/src/elements/mitc4_plusd.rs`. Verification: self-test asserting the `z` offsets sit on the four interior nodes only and that the warped patch's element normals differ from `n_vec`; the fixture is asserted **not** to be used for any constant-stress assertion. Satisfies: Requirement 1 item 8 (warped discrimination); Requirement 11 (variant (c)). <!-- sdd-owner: implementation -->
- [ ] **2.4 Boundary-traction loader and the 48-DOF dense patch assembler** `assemble_star_patch(&[Mitc4PlusDPrecomputed; 5]) -> DMatrix<f64>`, with Gauss-rule boundary integration and a dense LU solve. Touches `…/tests/fixtures.rs`. Verification: the self-test asserts the assembled 48×48 matrix is symmetric and that its six rigid-body fields carry zero energy; it is added and observed **failing** (RED) until WU4 lands. Satisfies: Requirement 6/7/8; Requirement 12. <!-- sdd-owner: implementation -->
  - **Deferred to WU4 (recorded):** this task needs `Mitc4PlusDPrecomputed`, which does not exist until WU2–WU4; writing the assembler now would produce a broken build rather than a red test. It lands with WU4 once the element type exists and stays unchecked here.

## 3. Element core — geometry, coefficients, directors, kinematics (WU2)

**WU2** · Touches: `crates/aeroelast-core/src/elements/mitc4_plusd.rs`, `…/mitc4_plusd/tests.rs`, `crates/aeroelast-core/src/elements/mod.rs` · Verification: the geometry and kinematics tests below · Rollback: delete the module + the `mod.rs` line · **~300 changed lines.**

- [ ] **3.1 Module skeleton.** Create `crates/aeroelast-core/src/elements/mitc4_plusd.rs` with `#[cfg(test)] mod tests;`, and add the single line `pub mod mitc4_plusd;` to `crates/aeroelast-core/src/elements/mod.rs`. Verification: `cd crates && cargo build -p aeroelast-core` compiles, and `git diff --stat crates/aeroelast-core/src/elements/mitc4.rs` is empty (the hybrid is byte-identical). Satisfies: Requirement 1; Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **3.2 `Mitc4PlusDPrecomputed` struct and constructor**, with `compute_characteristic_vectors`, `compute_membrane_coefficients_2017`, `compute_j3d_enriched`, `regularized_inverse_2x2`, `covariant_to_local_mapping`, `shear_covariant_to_local`, `shape_functions`, `shape_function_derivatives`; every stored quantity carries its paper/equation doc comment. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_geometry_` — `test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18`, `test_geometry_dual_basis_identities_eq11`, `test_geometry_a_E_is_positive_eq27c`. Satisfies: Requirement 3; Requirement 13. <!-- sdd-owner: implementation -->
- [ ] **3.3 `compute_node_directors` (ADR-4 option B)**: area-weighted average of the element's own four sub-quad normals, sign-aligned to `n_vec`, storing `vn`, `v1`, `v2` and `a_i: [f64; 4]`. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_geometry_node_directors_reduce_to_n_vec_when_flat` (design-added name) — `V_n^i = n_vec` to `1e-14` on a flat square and flat distorted quad, with a non-vacuity assertion that the four directors differ by `> 1e-6` relative on F-W. Satisfies: Requirement 1 item 2 (A Eq. 8a); Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **3.4 Kinematics interpolation** (A Eqs. (1)–(3)), including the `u_b = ½ Σ a_i h_i (θ_i × V_n^i)` rotation convention and `θ_z` carried but unconsumed by the 2017 core. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_displacement_field_matches_eq1_to_eq3` — including `t = ±1` and the identity `θ × V_n = −V_2 α + V_1 β`. Satisfies: Requirement 2. <!-- sdd-owner: implementation -->

## 4. Element core — B-operators (WU3)

**WU3** · Touches: `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`, `…/tests/drill.rs` · Verification: the operator tests below · Rollback: revert the operator block (the module does not yet dispatch) · **~350 changed lines.** Depends on WU0 (task 1.2) and WU2.

- [ ] **4.1 The five covariant membrane tying rows** (`b_rr_a` … `b_rs_e`, A Eqs. (15)–(17) at A/B/C/D/E) **and `b_membrane_2017`** (A Eq. 27, with the five `a_*` and the point-wise covariant→local mapping). Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1a_membrane_eq22_flat_tying_condition` — `ẽ_rs^m|bil == e_rs^m|bil` to `1e-14` absolute on the flat element, and the warped comparison non-zero beyond `1e-6` relative (non-vacuity). Satisfies: Requirement 1 item 5; Requirement 6 (Eq. 22 scenario). <!-- sdd-owner: implementation -->
- [ ] **4.2 `b_bending_2017` → `(B_b1, B_b2)`** from A Eqs. (7c)/(7d), displacement-based, including the `∂x_b·∂u_m` contribution of A Eq. (8a). Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_bending_operator_matches_eq7c_eq7d` — agree to `1e-10` relative at every Gauss point, and the local matrix is 24×24 with no condensed internal DOF. Satisfies: Requirement 1 item 4; Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **4.3 `b_shear_mitc4` (DB84 Eq. 3) + the four stored tying operators + `shear_covariant_to_local`**, using the metric normalization recorded in task 1.4 (the design's §2.2 3D-dual-basis form if the paper does not print one) — no invented factor. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: a flat-reduction-to-Mindlin unit test (no spec-fixed name; the design names it oracle (a)); the load-driven oracles are task 7.4 (`test_t1a_shearing_patch_constant_stress_fig5_mesh`), task 8.3 (`test_t1b_strong_patch_shearing_constant_and_zero_stress`) and task 5.4. Satisfies: Requirement 1 item 6. <!-- sdd-owner: implementation -->
- [ ] **4.4 `b_drill_membrane_2025` (B Eq. 18) + `drill_midside_shape_derivatives` (kept name) + the test-local reference `tests::drill::b_md_reference`** written from the printed Eq. (18) alone (independently recomputed `x_m^I`, `x_r^I`, `x_s^I`, `j`, `θ^D`, paper edge order). Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`, `…/tests/drill.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_drill_operator_matches_eq18_term_by_term` — entry-by-entry agreement at 9 sample points (4 Gauss + 4 edge mid-points + centre) to `1e-12` absolute on flat square, flat distorted, ruled and doubly warped quads, **plus Oracle 2's five asserted rejections** (`V^D → e3`; missing `1/8`; edge order `[bottom,right,top,left]`; flipped edge-difference sign; `θ_z` alone instead of `θ·V^D`). The reference must catch a wrong transcription and must **not** be silenceable by editing the reference to match the code: if the two disagree, the resolution is a recorded vision re-read in the extract with a page citation. Satisfies: Requirement 1 item 7; Requirement 14 (revived operator retained). <!-- sdd-owner: implementation -->
- [ ] **4.5 The `c_r`/`c_s` collision test**: assert `pre.c_r_mem`/`c_s_mem` (paper A, `x_d·m^r`, `x_d·m^s`) and the locals `cr_md`/`cs_md` (paper B, `x_m^I·(−x_r^I × V^D)`, `x_m^I·(x_s^I × V^D)`) are different quantities on a warped element. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core` with the new test (no spec-fixed name; the design fixes it as "a test asserts they are different quantities on a warped element"). Satisfies: Requirement 3; Requirement 1 items 5 and 7. <!-- sdd-owner: implementation -->

## 5. Element core — stiffness assembly (WU4)

**WU4** · Touches: `mitc4_plusd.rs`, `mitc4_plusd/tests.rs` · Verification: the identity/kinematics tests below · Rollback: revert the assembly block · **~340 changed lines.** Depends on WU3.

- [ ] **5.1 `resultant_moment_matrix`** (the `W_00 … W_22` block matrix of design §2.3, with `W_22 = cm/9` from paper A's 2×2 `t`-rule), with the multi-ply approximation documented in the doc comment. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: a unit test of the moment entries against their closed forms; no spec-fixed name exists, so this is stated as a design-derived unit test rather than a spec test. Satisfies: Requirement 1 item 8. <!-- sdd-owner: implementation -->
- [ ] **5.2 `compute_ke_local` / `compute_ke_global` + the drill contribution + the `cs_uncorrected` wiring + the test-local reference implementation.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_ke_lock_matches_2017_core_plus_2025_drill` — `max|K_prod − K_ref| ≤ 1e-10·max|K_ref|` on flat square, flat distorted, ruled warped and doubly warped quads, with the same bound on the membrane and transverse-shear blocks. Satisfies: Requirement 1 scenario 1. <!-- sdd-owner: implementation -->
- [ ] **5.3 Drill stiffness provenance.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_drill_stiffness_comes_only_from_eq26` — non-zero on warped geometry; exactly symmetric; **exactly zero on every translational row/column block**; exactly zero on every rotation block other than the drill's own on flat geometry; and `|u_rbᵀ K u_rb| ≤ 1e-12·λ_max·‖u_rb‖²` for the six rigid-body fields. Satisfies: Requirement 1 scenario 2 (design §9.1 interpretation, already reflected in spec rev 4). <!-- sdd-owner: implementation -->
- [ ] **5.4 Uncorrected transverse shear — the discriminating test.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_uses_uncorrected_shear_modulus` (block vs closed-form `∫ B_γᵀ (G·I) B_γ dA` to `1e-10`; the `5/6` value rejected by `> 1e-3` relative, asserted) and `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_invariant_to_shear_correction_factor` (ADR-1's two-`pre` construction on isotropic and single-ply laminate, plus the two non-vacuity controls: the `k`-carrying control differs by `> 1e-3` relative). Satisfies: Requirement 16. <!-- sdd-owner: implementation -->
- [ ] **5.5 Integration rule.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` on the strongly warped quad — matches the three-term reference to `1e-10` and differs from the surface-only reference by `> 1e-4` relative (asserted). Satisfies: Requirement 1 item 8. <!-- sdd-owner: implementation -->
- [ ] **5.6 Local matrix shapes.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_local_matrices_are_24x24` — `K`, `M`, `K_T` exactly 24×24 and `f_int` exactly 24 long, drilling at slot `6i+5`. Satisfies: Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **5.7 Drill-DOF energy behaviour.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_drill_dof_is_theta_z_through_eq26_operator` — energy `≤ 1e-12·λ_max·‖u‖²` for a pure rigid rotation about `V_n`, and non-zero energy on the warped `θ_z` pattern only through the Eq. (26) operator. Satisfies: Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **5.8 Mid-surface restriction (ADR-6 / G7).** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_element_uses_midsurface_constitutive` — the element's constitutive equals `Laminate::to_shell_constitutive()` — plus the supporting static check `grep -n "to_shell_constitutive_with_offset" crates/aeroelast-core/src/elements/mitc4_plusd.rs` returns no match (supporting, not normative). Satisfies: Requirement 16 scenario 2. <!-- sdd-owner: implementation -->

## 6. Element core — assembly-facing API (WU5)

**WU5** · Touches: `mitc4_plusd.rs` · Verification: the consistency and mass tests below · Rollback: revert the API block · **~300 changed lines.** Depends on WU4.

- [ ] **6.1 `compute_fint_global` (linear + the bounded nonlinear path) and `compute_kt_global`.** Touches `mitc4_plusd.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kt_zero_matches_ke`, `… test_fint_linear_nonlinear_parity`, and the directional-derivative triple including `test_kt_fint_directional_derivative_with_drill_dofs`. Satisfies: Requirement 2; Requirement 12 (T2B consistency oracle). <!-- sdd-owner: implementation -->
- [ ] **6.2 `compute_me_global` / `compute_me_composite_global`.** Touches `mitc4_plusd.rs`. Verification: `cd crates && cargo test -p aeroelast-core me_global` — the mass invariants including the design's `test_me_global_total_translational_mass_is_rho_h_a`. Satisfies: Requirement 12. <!-- sdd-owner: implementation -->
- [ ] **6.3 `compute_body_load_global`, `compute_k_sigma_global`, `compute_centrifugal_prestress`, `compute_element_stress`.** Touches `mitc4_plusd.rs`. Verification: the spec fixes no new name for these; the verification is the retargeted T2I names once moved (task 9.2). Stated rather than invented. Satisfies: Requirement 12. <!-- sdd-owner: implementation -->
- [ ] **6.4 `build_t24` / `transform_to_global` / `extract_elem_disp_24` and the corotational machinery** (`quaternion_to_matrix`, `quaternion_from_vector`, `rotate_vector_by_quaternion`, `quaternion_multiply`, `update_normals_with_displacements`, `polar_decomposition`, `log_strain_from_polar`, `compute_membrane_strain_log`, `update_corotational_frame`, `frame_incremental_rotation`). Touches `mitc4_plusd.rs`. Verification: `test_fint_linear_nonlinear_parity` and the large-rotation Python benchmarks (`tests/test_large_rotation_benchmarks.py`, judged at the S3 gate, task 10.7). Satisfies: Requirement 12. <!-- sdd-owner: implementation -->

## 7. Tier 1a tests — the 2017 core (WU6)

**WU6** · Touches: `mitc4_plusd/tests.rs` · Verification: T1.1–T1.3 below · Rollback: revert the test block · **~300 changed lines.** Depends on WU1; the tests are written and observed RED before WU4/WU5 complete.

- [ ] **7.1 Rigid-body fixture + T1.2.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1a_zero_energy_modes_single_unsupported_element_exactly_six` — exactly 6 eigenvalues with `|λ| ≤ 1e-10·λ_max` on flat, flat-distorted and warped elements, `‖K u_rb‖∞ ≤ 1e-10·λ_max·‖u_rb‖∞` for each of the six fields, separation `≥ 1e-6·λ_max`. Satisfies: Requirement 5. <!-- sdd-owner: implementation -->
- [ ] **7.2 T1.1.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1a_isotropy_element_orientation_and_node_sequence_invariant` — eigenvalues to `1e-10·λ_max` across ≥ 4 orientations including a `π/2` rotation, permuted node sequences to `1e-12·max|K|`, `uᵀKu` invariant to `1e-10` relative. Satisfies: Requirement 4. <!-- sdd-owner: implementation -->
- [ ] **7.3 T1.3a membrane patch.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1a_membrane_patch_constant_stress_fig5_mesh` — the constant state to `1e-8` relative (absolute floor `1e-10·‖σ‖`) for each of `σ_xx`, `σ_yy`, `τ_xy` alone, spread `≤ 1e-8·‖σ‖`. Satisfies: Requirement 6. <!-- sdd-owner: implementation -->
- [ ] **7.4 T1.3b bending and T1.3c shearing patches.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1a_bending_patch_constant_curvature_fig5_mesh` (constant curvature to `1e-8` relative, spread `≤ 1e-8`) and `… test_t1a_shearing_patch_constant_stress_fig5_mesh` (constant shear to `1e-8` relative, spread `≤ 1e-8`, no shear correction factor). Satisfies: Requirement 7; Requirement 8. <!-- sdd-owner: implementation -->

## 8. Tier 1b tests — the 2025 six-DOF element (WU7)

**WU7** · Touches: `mitc4_plusd/tests.rs` · Verification: T1.4–T1.6 below · Rollback: revert the test block · **~330 changed lines.** Depends on WU1 and WU5.

- [ ] **8.1 T1.4.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1b_spatial_isotropy` — eigenvalues to `1e-10·λ_max` across ≥ 4 orientations with the drill free, and `uᵀKu` invariant to `1e-10` relative for a co-rotated field with a non-zero drilling component. Satisfies: Requirement 9. <!-- sdd-owner: implementation -->
- [ ] **8.2 T1.5.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1b_zero_energy_modes_exactly_six_with_drill_dof` — exactly 6 zero eigenvalues on flat, warped and flat-distorted elements, each rigid-body field `≤ 1e-12·λ_max` at `‖u‖ = 1`, separation `≥ 1e-6·λ_max`, and non-zero drill stiffness in at least one non-rigid mode. Satisfies: Requirement 10. <!-- sdd-owner: implementation -->
- [ ] **8.3 T1.6a–c strong-form patches.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1b_strong_patch_extension_constant_and_zero_stress`, `… test_t1b_strong_patch_bending_constant_and_zero_stress`, `… test_t1b_strong_patch_shearing_constant_and_zero_stress` — the constant state to `1e-8` relative in every element and every analytically-zero component `≤ 1e-10` of the state. Satisfies: Requirement 11. <!-- sdd-owner: implementation -->
- [ ] **8.4 T1.6d.** Touches `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_t1b_drill_theta_z_free_except_corner_b` — variants (a) and (b) agree to `1e-10` relative, and variant (c) differs from (a) by `> 1e-8` relative on the warped patch F-W. Satisfies: Requirement 11. <!-- sdd-owner: implementation -->

## 9. Move the layout-bound Tier-2 Rust tests onto the new element (S2, WU8)

**WU8** · Touches: `crates/aeroelast-core/src/elements/mitc4.rs` (tests only), `mitc4_plusd/tests.rs` · Verification: `cd crates && cargo test -p aeroelast-core` at the baseline · Rollback: revert the move · **~300 changed lines (mostly moves; ~0 net).**

- [ ] **9.1 Move T2A/T2B** (the twelve names of design §5.1) from `mitc4.rs`'s test module into `mitc4_plusd/tests.rs`, retargeted from `Mitc4Precomputed` to `Mitc4PlusDPrecomputed`; **exclude** `test_ke_local_eigenvalues_nonsymmetric`. Touches both files. Verification: the moved names pass against the new element; `cd crates && cargo test -p aeroelast-core` reports 120+ passed / 0 failed; each moved name exists exactly once (`grep -rn "fn test_ke_global_has_exactly_six_zero_modes" crates/` returns one hit). Satisfies: Requirement 12. <!-- sdd-owner: implementation -->
- [ ] **9.2 Move T2I** (`me_global*`, `body_load`, `k_sigma`, `centrifugal_prestress`, `compute_element_stress`), retargeted the same way. Touches both files. Verification: same command and the same "exactly once" grep; assertions are mechanical 6-DOF-layout equivalents, with no tolerance changed. Satisfies: Requirement 12. <!-- sdd-owner: implementation -->
- [ ] **9.3 Record the S2 gate.** Touches the change record. Verification: `cd crates && cargo test -p aeroelast-core` at the baseline (120 passed / 0 failed) with the hybrid's copies removed and the moved tests passing against `Mitc4PlusDPrecomputed` **before** the flip. Satisfies: Requirement 14 scenario 1 (ordering evidence). <!-- sdd-owner: implementation -->

## 10. The flip (S3, WU9)

**WU9** · Touches: `assembly/assembler.rs`, `aeroelast-py/src/{elements,assembler}.rs`, `materials/{mod,laminate}.rs`, `tests/test_ko2017_performance.py` · Verification: the S3 gate (task 10.7) · Rollback: one line per dispatch site — the hybrid path returns · **~220 changed lines.** Depends on WU8.

- [ ] **10.1 ADR-1 material channel: `ShellConstitutive::transverse_shear_uncorrected(applied_k)` and `Laminate::applied_shear_correction_factor()`**, both additive; `ShellConstitutive` keeps exactly its five fields and `Laminate`'s public behaviour and numbers are untouched. Touches `materials/mod.rs`, `materials/laminate.rs`. Verification: `cd crates && cargo test -p aeroelast-core materials::` plus new unit tests for both accessors (`k > 0` and `k = 0`); no spec-fixed name exists for these accessor unit tests, stated rather than invented. Satisfies: Requirement 16; Requirement 15. <!-- sdd-owner: implementation -->
- [ ] **10.2 `MaterialSpec::Composite` gains the internal `applied_shear_correction` field**, filled in `crates/aeroelast-py/src/assembler.rs` from `corrected_lam.applied_shear_correction_factor()`. Touches `assembly/assembler.rs`, `aeroelast-py/src/assembler.rs`. Verification: `cd crates && cargo build -p aeroelast-py` (via the maturin build in task 10.7) and the T2G laminate groups at the gate; a single-ply laminate routed through the production path removes its `0.75`, a multi-ply laminate's factor is `1.0` and its `cs` passes through verbatim. Satisfies: Requirement 16; Requirement 15. <!-- sdd-owner: implementation -->
- [ ] **10.3 Dispatch the element.** `PrecomputedElem::Quad` payload → `Mitc4PlusDPrecomputed`; both `Mitc4Precomputed::new` sites (`assembler.rs:161-167`, `:252-254`); `mitc4::…` → `mitc4_plusd::…`; `build_constitutive_mitc4_plusd`; `extract_elem_disp_24`'s return type; `update_reference`; `ElemType::{Mitc4, Mitc4Composite}`. Touches `assembly/assembler.rs`. Verification: `cd crates && cargo test -p aeroelast-core` (all Tier 1 + Tier 2 Rust green at the gate, task 10.7) and `git diff --stat crates/aeroelast-core/src/assembly/topology.rs` is empty. Satisfies: Requirement 1; Requirement 14 scenario 1. <!-- sdd-owner: implementation -->
- [ ] **10.4 PyO3 entry points.** Swap the MITC4 kernels' internals to `mitc4_plusd`; the Python-facing signatures, `[f64; 576]` / `[f64; 24]` shapes and every `#[pyfunction]`/`#[pyclass]` name stay frozen; the now-unused `e_mod` argument is kept and ignored with a comment naming the out-of-scope signature change. Touches `aeroelast-py/src/elements.rs`. Verification: the maturin rebuild plus `python -m pytest "tests/test_rust_composite.py::TestMITC4BatchSanity" -q` (after the rebuild), and the supporting static checks that `[f64; 576]`/`[f64; 24]` and the function names are still declared. Satisfies: Requirement 2; Requirement 14 scenario 4. <!-- sdd-owner: implementation -->
- [ ] **10.5 Record the composite batch entry point's open factor channel.** `batch_ke_mitc4_composite` (`aeroelast-py/src/elements.rs:532-540`) takes no factor argument, so `applied_shear_correction = 1.0` is passed there and the caller's `cs` is consumed verbatim; document this in the code and in the change record. Touches `aeroelast-py/src/elements.rs`. Verification: the code comment plus the change record; **no test can cover the missing channel** (a factor argument would be an out-of-scope PyO3 signature change) — stated, not invented. Satisfies: Requirement 16 (open item). <!-- sdd-owner: implementation -->
- [ ] **10.6 Correct the four flagged benchmark mis-sourcings** in `tests/test_ko2017_performance.py` (`test_3_2` SS rows vs Table 6; `test_3_3[dist]`'s hardcoded Table-12 value; `test_3_7[reg]`'s N=8/S4 cells; the `test_3_5` docstring/xfail mismatch) to the true paper cells. Touches `tests/test_ko2017_performance.py`. Verification: `python -m pytest "tests/test_ko2017_performance.py" -q` — correcting a *source* is not loosening a *tolerance*, and no `rtol` changes. Must land **before** the S3 gate is judged, because a mis-sourced cell could otherwise let a formulation error pass. Satisfies: Requirement 12. <!-- sdd-owner: implementation -->
- [ ] **10.7 The S3 gate.** Touches nothing (a recorded run). Verification, all of: `cd crates && cargo test -p aeroelast-core test_t1a_`, `… test_t1b_`, `cd crates && cargo test -p aeroelast-core` (Tier 2 Rust, baseline 120/0), the maturin rebuild, `python -m pytest -m "not slow" -q` (baseline 345 passed / 2 failed / 2 skipped, the 2 failures excluded as pre-existing), the laminate/composite preserved-invariant command of the spec, `python -m pytest "tests/test_ko2017_performance.py" -q`, and the PyO3 surface check. Record the run. Satisfies: Requirement 12; Requirement 14 scenario 1; Requirement 15. <!-- sdd-owner: implementation -->

## 11. The retirement (S4, WU10)

**WU10** · Touches: `crates/aeroelast-core/src/elements/mitc4.rs` · Verification: the S4 gate (task 11.5) · Rollback: revert the retirement commit — the hybrid returns · **net-negative lines (≈ 940 removed).** Depends on task 10.7 being green.

- [ ] **11.1 Delete the hybrid's deviation surfaces and dead code** from `mitc4.rs`: the Winkler & Plakomytis ERC block, `beta_w`/`BETA_W`, `k_drill`/`drilling_scale`, the SRI split, the 2-DOF `(1−ξ²)(1−η²)` rotation bubble with `GpBubble`/`b_kappa_bubble`/`b_gamma_mitc4_plus`, the A Eqs. (18)–(19) membrane as the live field, the element-side shear-correction usage, and `b_m_standard`, `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`. Touches `mitc4.rs`. Verification: the S4 gate (task 11.5) plus the grep of task 11.4. `b_md_mitc4_plus` / `drill_midside_shape_derivatives` are **not** deleted — they are revived in `mitc4_plusd.rs` (task 4.4). Satisfies: Requirement 14 scenario 2. <!-- sdd-owner: implementation -->
- [ ] **11.2 Delete the T2J tests** (`test_drill_midside_derivatives_match_paper_eq11`, `test_drill_membrane_operator_*` ×5, `test_b_erc_*` ×5, `test_erc_covariant_to_local_*` ×2, `test_ke_local_erc_flat_rigid_body_and_symmetry`, `test_enhanced_drill_stiffness`, `test_drill_warping_moment_flat_element`, `test_ke_local_eigenvalues_nonsymmetric`) from `mitc4.rs`. Touches `mitc4.rs`. Verification: `cd crates && cargo test -p aeroelast-core` at the baseline; the names are explicitly **not** part of the pass condition. Satisfies: Requirement 12; Requirement 14. <!-- sdd-owner: implementation -->
- [ ] **11.3 Rewrite the `mitc4.rs` module header** so it no longer documents the deviation text or the retired element. Touches `mitc4.rs`. Verification: `grep -nE 'k_drill|drilling_scale|beta_w|ERC|shear correction' crates/aeroelast-core/src/elements/mitc4.rs` returns no match; the header states what remains. Satisfies: Requirement 13; Requirement 14 scenario 2. <!-- sdd-owner: implementation -->
- [ ] **11.4 The retirement tests and the static check.** Add `test_retirement_removes_hybrid_deviation_surfaces` and the `test_retirement_*` gate tests. Touches `mitc4.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_retirement_` plus the supporting grep `grep -nE 'k_drill|drilling_scale|compute_ke_local_erc|beta_w|hg_stiffness_factor|cm_normal|b_m_standard|green_lagrange_strain|compute_b_l|compute_membrane_stress|Mat26|Vec26'` over the module the MITC4+/D type dispatches to returns no match (supporting, not normative). **Process obligation:** these tests must not be present-and-green while any Tier-1 or Tier-2 acceptance test is failing — the S3 gate run of task 10.7 is the *before* evidence and is recorded in the change. Satisfies: Requirement 14 scenarios 1–2. <!-- sdd-owner: implementation -->
- [ ] **11.5 The S4 gate, re-run.** Touches nothing (a recorded run). Verification: task 10.7's full command set re-run after the deletion (Tier 1 + Tier 2 + the static checks + the PyO3 surface check + `python -m pytest "tests/test_rust_composite.py::TestMITC4BatchSanity" -q` after the rebuild), with both baselines unchanged. Satisfies: Requirement 12; Requirement 14. <!-- sdd-owner: implementation -->

## 12. Docs and guard tests (WU11)

**WU11** · Touches: `docs/formulations/shell-elements.md`, `docs/validation-matrix.md`, `tests/test_mitc4plusd_traceability.py`, `tests/test_laminate_invariant_guard.py` · Verification: the three traceability scenarios, the laminate-guard scenario and the Python baseline · Rollback: revert the docs and the two test files · **~300 changed lines.**

- [ ] **12.1 `tests/test_mitc4plusd_traceability.py`** with its three tests, and the ingredient table it parses: `docs/formulations/shell-elements.md` §2 and the `mitc4_plusd.rs` module header. Touches both docs and the new test file. Verification: `python -m pytest "tests/test_mitc4plusd_traceability.py::test_traceability_identity_table_names_source_for_every_ingredient" "tests/test_mitc4plusd_traceability.py::test_traceability_no_forbidden_ingredient_claimed_or_present" "tests/test_mitc4plusd_traceability.py::test_traceability_equation_citations_resolve_to_extract_or_paper" -q` — every Requirement-1 key carries a citation that resolves against an extract, no forbidden ingredient is claimed or present, and every cited ingredient is referenced by a `test_identity_* | test_kinematics_* | test_geometry_* | test_t1a_* | test_t1b_*` name. Source scope is the fixed live list `["mitc4_plusd.rs"]`, unchanged before and after the retirement. Satisfies: Requirement 13. <!-- sdd-owner: implementation -->
- [ ] **12.2 `tests/test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged`.** Touches the new test file. Verification: `python -m pytest "tests/test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged" -q` — `ShellConstitutive` still has exactly `cm`, `cb_coupling`, `cb`, `cs`, `cm_raw`, and `Laminate` still exposes `to_shell_constitutive`, `to_shell_constitutive_with_offset`, `shear_correction_factor`, `abd_matrix_flat`. Satisfies: Requirement 15. <!-- sdd-owner: implementation -->
- [ ] **12.3 Rewrite `docs/formulations/shell-elements.md` §2** so it describes the element that actually runs (the MITC4+/D) with the parseable ingredient table, and resolve §4.2's "pending MITC4/D drill". Touches `shell-elements.md`. Verification: task 12.1's three scenarios pass against the rewritten §2; §2 no longer documents the SRI split, the `k_drill` penalty or the ERC treatment as live. Satisfies: Requirement 13; Requirement 14. <!-- sdd-owner: implementation -->
- [ ] **12.4 Refresh `docs/validation-matrix.md`** (the pass/fail columns, the §4.2 pending entry, and the four mis-sourcings already corrected in code by task 10.6). Touches `validation-matrix.md`. Verification: the refreshed matrix agrees with the recorded S4 gate run; the four corrected cells match `tests/test_ko2017_performance.py`'s true paper sources; `python -m pytest -m "not slow" -q` returns the baseline. Satisfies: Requirement 12. <!-- sdd-owner: implementation -->

## 13. Change close-out

- [ ] **13.1 Run the change close-out gate and confirm the three end states.** Touches nothing (a recorded run). Verification, all of: (a) the **Tier 1 + Tier 2 gate** — `cd crates && cargo test -p aeroelast-core test_t1a_`, `… test_t1b_`, `cd crates && cargo test -p aeroelast-core` (120 passed / 0 failed), the maturin rebuild, `python -m pytest -m "not slow" -q` (345 passed / 2 failed / 2 skipped with the two pre-existing failures excluded and named), the spec's laminate/composite preserved-invariant command, and `python -m pytest "tests/test_ko2017_performance.py" -q`; (b) the **retirement** — `cd crates && cargo test -p aeroelast-core test_retirement_` green *and* the task 11.4 grep returning no match, with the S3 gate run recorded as the before-evidence; (c) the **documentation rewrite** — the `mitc4.rs` module header (task 11.3), `docs/formulations/shell-elements.md` §2 (task 12.3) and `docs/validation-matrix.md` (task 12.4) no longer describe the retired element, and the three traceability scenarios plus the laminate guard pass. Record the whole run in the change. Satisfies: Requirement 12; Requirement 13; Requirement 14; Requirement 15. <!-- sdd-owner: implementation -->
- [ ] **13.2 Start or reuse the bounded post-apply review and record the lifecycle-gate decision**, including the accepted `size:exception` against the 700-line session budget and the forecast total of ≈ 3,760 changed lines. Verification: the review's own report against the task list, and the recorded decision; no test covers this (it is a process gate). Satisfies: the change's delivery decision. <!-- sdd-owner: parent -->

---

## Open items carried, and where each is settled

No open item of the design disappears. Each is either a task, an explicit deferral with
the WU that carries it, or a stated assumption.

| Design open item | Carried as |
| --- | --- |
| 1 — the composite per-element PyO3 batch entry point has no factor channel | Task **10.5**: `applied_shear_correction = 1.0` there, documented; the production path is fully fixed in 10.1/10.2. Explicit deferral until PyO3 signature changes come into scope. |
| 2 — the Eq. (18) edge convention | Task **1.2** (settles it) with task **4.4**'s Oracle 2 sign-flip row as the automated detector. |
| 3 — paper B Eqs. (1)–(16), the `h̃_m` functions and the `θ_z` interpolation are untranscribed | Task **1.1** — the explicit first task, with its own verification (recorded vision re-read + the traceability line of task 1.5); it is the precondition for task 4.4's term-by-term drill verification. |
| 4 — the covariant transverse-shear metric normalization is not printed in paper A | Task **1.4** (search + record) / task **4.3** (implement + the three oracles). |
| 5 — the nonlinear path's fidelity | Task **6.1** implements and bounds it; T2B's consistency set is the oracle. Explicit deferral of a paper-faithful nonlinear derivation. |
| 6 — the `t⁴` moment for multi-ply laminates uses `cm/9` | Task **5.1** names the approximation in the doc comment; task **5.5** makes it visible to a test. |
| 7 — a tapered per-node thickness | Task **3.3** stores `a_i: [f64; 4]` in the per-node form; **no taper case is tested** (no input exists). |
| 8 — the reading of the spec's "row/column block" clause | Resolved in spec rev 4 (§9.1's wording was adopted); task **5.3** asserts the interpretation and task **4.4**'s Oracle 2 asserts the `θ_z`-only variant is rejected. |
| 9 — the offset deferral (G7) | Task **5.8** asserts the mid-surface path and the supporting static check; explicit deferral recorded. |
| 10 — the element-local director field (ADR-4 option B) | Task **3.3** pins the flat reduction and the warped non-vacuity; option (A) is the named follow-up. Stated limitation, not a task. |

## What could not be turned into a task

- **A paper-faithful nonlinear MITC4+/D formulation.** Neither paper provides one for this
  repository's UL formulation, so task 6.1 bounds the nonlinear path to the repository's
  existing total-Lagrangian covariant correction and relies on T2B's consistency set.
  Deferred, not invented.
- **A test for the missing factor channel of `batch_ke_mitc4_composite`.** Closing it needs
  an out-of-scope PyO3 signature change; task 10.5 records the gap instead.
- **A tapered-thickness case.** No per-node thickness input exists in this change; task 3.3
  keeps the door open but tests no taper.
- **Automated verification of WU0's Markdown derivations.** No oracle exists for Markdown
  math; WU0's checks are recorded vision re-reads, with the downstream traceability test
  (task 12.1) as the only automated hook. Stated rather than faked with an invented test.
- **The exact paper boundary-node coordinates** where the figure read cannot recover them
  (G1/G3). The documented-equivalent-patch fallback of proposal §2.1 applies; the deviation
  is recorded, and no mesh is invented.
- **A chained-PR split.** The session accepted `single-pr (size:exception)`; the S1–S2 / S3
  / S4 / WU11 split is recorded as a fallback in the forecast and is deliberately not chosen.
