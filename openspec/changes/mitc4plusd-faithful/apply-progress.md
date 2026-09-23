# Apply progress — the MITC4+/D element, faithfully implemented

Change: `mitc4plusd-faithful` · phase: **apply** · store: openspec (Engram mirror:
`sdd/mitc4plusd-faithful/apply-progress`).

Work units: **WU1 — fixtures (tasks 2.1, 2.2, 2.3)** and **WU2 — element core:
geometry, coefficients, directors and kinematics (tasks 3.1, 3.2, 3.3, 3.4)**. This
file is cumulative; the WU1 record is preserved below and the WU2 record follows it.

## Structured status consumed

- Source: native SDD status engine (authoritative, `artifactStore: openspec`).
- `changeName`: `mitc4plusd-faithful`; `applyState`: `ready`; `nextRecommended`: `apply`.
- `actionContext.mode`: `repo-local`; `workspaceRoot`:
  `/home/efirvida/Desktop/dev/fem-shell`; `allowedEditRoots`:
  `/home/efirvida/Desktop/dev/fem-shell`. All edits stayed inside the workspace and
  inside the change's authorized edit roots.
- `taskProgress` at entry: 59 total / 5 completed / 54 pending.
- Review workload gate: `Decision needed before apply: No`, `Chained PRs recommended:
  No`, `Chain strategy: size-exception`, `400-line budget risk: High`. The session
  resolved delivery as **single-pr with an explicitly accepted `size:exception`** and a
  700-line review budget, so WU1 proceeded. WU1 landed at ~494 changed lines, above the
  design's ~200 forecast and the ~350 per-unit guide; this is covered by the accepted
  session `size:exception` and recorded here.

## Completed tasks (persisted checkboxes updated in `tasks.md`)

- [x] **2.1** Star-patch fixture (`STAR_NODES`, `STAR_ELEMS`).
- [x] **2.2** F-A-BC (derived) and F-B-BC (figure-read) boundary-condition sets.
- [x] **2.3** F-W warped star patch plus geometry self-tests.

Task 2.4 remains `- [ ]` and its text now carries the recorded deferral to WU4.

## Files changed

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | **new** — module header + inline `#[cfg(test)] mod tests` with the fixtures and 10 self-tests (~494 lines) |
| `crates/aeroelast-core/src/elements/mod.rs` | +1 line: `pub mod mitc4_plusd;` |
| `openspec/changes/mitc4plusd-faithful/tasks.md` | 2.1–2.3 checked; 2.4 deferral note; WU1 layout-deviation note |
| `openspec/changes/mitc4plusd-faithful/apply-progress.md` | this file |

Diff stat (intent-to-add for the new file):

```text
crates/aeroelast-core/src/elements/mitc4_plusd.rs | 494 ++++++++++++++++++++++
crates/aeroelast-core/src/elements/mod.rs         |   1 +
2 files changed, 495 insertions(+)
```

No file outside the two authorized Rust files and the change's openspec artifacts was
touched. `mitc4.rs`, `mitc3.rs`, `quad.rs`, the materials and the assembler are
byte-identical (`git diff --stat` empty).

## Verification

Command (workspace root is `crates/`):

```text
cd crates && cargo test -p aeroelast-core
```

Result: **130 passed / 0 failed** (baseline 120 + the 10 new fixture self-tests).
`cargo clippy -p aeroelast-core --all-targets` reports **no warning attributable to
`mitc4_plusd.rs`**. `rustfmt --check` is clean for the new file.

## Fixture self-tests and what each asserts

| Test | Asserts |
| --- | --- |
| `test_star_patch_has_five_elements_and_eight_nodes` | `STAR_ELEMS` has 5 elements, `STAR_NODES` has 8 nodes, every connectivity index is in range |
| `test_star_patch_all_signed_areas_positive_and_sum_to_100` | every element's signed area is `> 0` (CCW) and the five areas sum to 100 (the 10×10 square) |
| `test_star_patch_has_no_duplicate_coordinates` | no two of the eight nodes share coordinates |
| `test_star_patch_central_element_matches_extract` | the central element's vertices are exactly `(2,2) → (8,3) → (8,7) → (4,7)`, the CCW rotation of the extract's `(2,2) → (4,7) → (8,7) → (8,3)` |
| `test_f_a_bc_constrains_exactly_six_rigid_body_modes` | `rank(constraint_matrix(F_A_BC) · rigid_body_modes) == 6` — the derived set removes exactly the six rigid-body modes |
| `test_f_b_bc_constrains_exactly_six_rigid_body_modes` | the same rank-6 assertion for the extension, bending and shearing F-B-BC sets |
| `test_f_b_bc_matches_fig7_bcd` | the exact `(node, dof)` triples of the 2025 extract §Fig. 7(b)(c)(d); every value is zero; `θ_z` (dof 5) is constrained only at corner B |
| `test_f_w_z_offsets_on_interior_nodes_only` | F-W keeps the four corners at `z = 0` and carries `±0.5` on the four interior nodes only, with every `xy` unchanged from F-A/F-B |
| `test_f_w_element_normals_differ_from_flat_n_vec` | the flat patch's four sub-quad normals equal `n_vec = e3` to `1e-14`, and every F-W element has a sub-quad normal deviating by `> 1e-3` rad (non-vacuity) |
| `test_f_w_is_not_a_constant_stress_fixture` | F-W is genuinely non-planar (no flat element), i.e. it cannot serve as a constant-stress fixture |

## TDD evidence (explicit test-first; `strict_tdd: false`)

The fixtures and their self-tests were written together; each self-test was then shown
to fail for a wrong fixture and pass for the right one.

| RED perturbation | Observed failure | Restored |
| --- | --- | --- |
| central element set clockwise `[4,7,6,5]` | `test_star_patch_all_signed_areas_positive_and_sum_to_100` → "element 0 … signed area -22"; `test_star_patch_central_element_matches_extract` → vertex 1 mismatch | yes |
| F-A-BC replaced by design §4.1's literal set (`D u_x`, `C θ_z`) | `test_f_a_bc_constrains_exactly_six_rigid_body_modes` → rank ≠ 6 | yes |
| F-W flattened (all interior `z = 0`) | `test_f_w_z_offsets_on_interior_nodes_only`, `test_f_w_element_normals_differ_from_flat_n_vec`, `test_f_w_is_not_a_constant_stress_fixture` all fail | yes |

After restoration: **130 passed / 0 failed**.

## Deviations from design

1. **Module layout (recorded in the file header).** Tasks name
   `…/mitc4_plusd/tests/fixtures.rs`; the design fixes `…/mitc4_plusd.rs` and every
   existing element is a single file with an inline `#[cfg(test)] mod tests`. Per
   `openspec/config.yaml`'s apply guideline, the fixtures live inline in
   `mitc4_plusd.rs`. The task's verification substring
   `mitc4_plusd::tests::fixtures` therefore becomes `mitc4_plusd::tests::`.

2. **F-A-BC correction (recorded on `F_A_BC` and in the file header).** Design §4.1's
   derived set (`u_x=u_y=u_z=0` at C; `u_x=0` at D; `u_y=0` at A; `θ_z=0` at C) has
   rank 4 on the six rigid-body modes: C=(0,0) and D=(10,0) share `y=0`, so their
   `u_x=0` rows are the same rigid-body equation, and neither `θ_x` nor `θ_y` is
   constrained. The fixture uses the corrected minimal set `C: u_x=u_y=u_z=0`,
   `A: u_y=0`, `C: θ_x=θ_y=0` (rank 6). The RED perturbation above is the proof that
   the design's literal set fails the required assertion.

3. **WU1 size.** ~494 changed lines vs the ~200 forecast (still above the ~350 per-unit
   guide). Covered by the accepted session `size:exception`; no test, doc or citation
   was dropped.

## Remaining tasks

All WU1 tasks 2.1–2.3 are complete. The exact remaining unchecked lines in section 2:

```text
- [ ] **2.4 Boundary-traction loader and the 48-DOF dense patch assembler** `assemble_star_patch(&[Mitc4PlusDPrecomputed; 5]) -> DMatrix<f64>`, with Gauss-rule boundary integration and a dense LU solve. Touches `…/tests/fixtures.rs`. Verification: the self-test asserts the assembled 48×48 matrix is symmetric and that its six rigid-body fields carry zero energy; it is added and observed **failing** (RED) until WU4 lands. Satisfies: Requirement 6/7/8; Requirement 12. <!-- sdd-owner: implementation -->
  - **Deferred to WU4 (recorded):** this task needs `Mitc4PlusDPrecomputed`, which does not exist until WU2–WU4; writing the assembler now would produce a broken build rather than a red test. It lands with WU4 once the element type exists and stays unchecked here.
```

WU2–WU11 remain pending.

## Workload / PR boundary

Single PR, accepted `size:exception`. WU1 is one self-contained review slice: new file
plus one `mod.rs` line, green tree, its own verification and rollback (delete the
module and the `mod.rs` line).

## Next

`next_recommended: parent-lifecycle` for WU1's own boundary is not applicable — the
change has pending tasks. The next implementable unit is **WU2** (task 3.1 module
skeleton, then 3.2–3.4), which builds on this fixtures module.

---

## WU2 — element core: geometry, coefficients, directors, kinematics

Work unit: **WU2 (tasks 3.1, 3.2, 3.3, 3.4)**. This section is appended to the WU1
record above; WU1's bytes are preserved.

### Structured status consumed (WU2)

- Source: native SDD status engine (authoritative, `artifactStore: openspec`).
- `changeName`: `mitc4plusd-faithful`; `applyState`: `ready`; `nextRecommended`:
  `apply`.
- `actionContext.mode`: `repo-local`; `workspaceRoot`:
  `/home/efirvida/Desktop/dev/fem-shell`; `allowedEditRoots`:
  `/home/efirvida/Desktop/dev/fem-shell`. All edits stayed inside the workspace and
  inside WU2's authorized edit roots (`mitc4_plusd.rs`, `tasks.md`, `apply-progress.md`).
- `taskProgress` at entry: 59 total / 8 completed / 51 pending.
- Review workload gate: `Decision needed before apply: No`, `Chained PRs recommended:
  No`, `Chain strategy: size-exception`, `400-line budget risk: High`. The session
  resolved delivery as **single-pr with an explicitly accepted `size:exception`** and a
  700-line review budget, so WU2 proceeded.

### Completed tasks (WU2, persisted checkboxes updated in `tasks.md`)

- [x] **3.1** Module skeleton — already satisfied by WU1; the recorded deviation (inline
  `#[cfg(test)] mod tests`, no separate `mitc4_plusd/tests.rs`) is noted on the task.
- [x] **3.2** `Mitc4PlusDPrecomputed` struct + constructor + the eight geometry helpers.
- [x] **3.3** `compute_node_directors` (ADR-4 option B), `vn`/`v1`/`v2`/`a_i`.
- [x] **3.4** Kinematics interpolation (`interpolate_position`, `interpolate_displacement`).

### Files changed (WU2)

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | extended: the WU2 production core (geometry, coefficients, directors, kinematics) plus the 5 new tests in the existing inline test module |
| `openspec/changes/mitc4plusd-faithful/tasks.md` | 3.1–3.4 checked; the 3.1 deviation note added |
| `openspec/changes/mitc4plusd-faithful/apply-progress.md` | this WU2 section |

Diff stat (`crates/aeroelast-core/src/elements/mitc4_plusd.rs`, tracked):

```text
crates/aeroelast-core/src/elements/mitc4_plusd.rs | 1013 +++++++++++++++++++++
1 file changed, 1013 insertions(+)
```

`git diff --stat crates/aeroelast-core/src/elements/mitc4.rs` is **empty**: the hybrid
is byte-identical. No file outside the authorized set was touched. `elements/mod.rs`
was not edited (the module was already registered in WU1).

### Verification (WU2)

Command (workspace root is `crates/`):

```text
cd crates && cargo test -p aeroelast-core
```

Result: **135 passed / 0 failed** (WU1's 130 + the 5 new WU2 tests).
`rustfmt --edition 2021 --check` is clean for the file.

### New tests and what each asserts

| Test | Asserts |
| --- | --- |
| `test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18` | on a flat rectangle `x_d = 0`, `c_r = c_s = 0`, `d = -1`, all five `a_*` zero; Eq. (27a)/(27b) collapse term-by-term onto Eq. (18a)/(18b) and Eq. (27c) onto Eq. (19c) for arbitrary sampled tying strains |
| `test_geometry_dual_basis_identities_eq11` | Eq. (11) on flat-distorted, ruled-warped and doubly-warped quads: `m^r.x_r = m^s.x_s = 1`, `m^r.x_s = m^s.x_r = 0`, `m^r.n = m^s.n = 0` |
| `test_geometry_a_E_is_positive_eq27c` | `a_E = +2 c_r c_s / d` (positive sign as printed, p. 410) on the coefficients' own inputs with `c_r c_s > 0`, `d > 0`; the same closed form on a real flat distorted element; the negated (deleted) form rejected by `> 1e-6` relative |
| `test_geometry_node_directors_reduce_to_n_vec_when_flat` | `V_n^i = n_vec` to `1e-14` for all four nodes on a flat square and a flat distorted quad; `V_1^i`, `V_2^i` orthonormal with a right-handed `(V_1, V_2, V_n)`; non-vacuity: the four directors differ by `> 1e-6` on the warped fixture F-W |
| `test_kinematics_displacement_field_matches_eq1_to_eq3` | Eq. (1) position at `t = 0, ±1`; the identity `θ x V_n = -V_2 α + V_1 β`; Eq. (3a) `u(t=0)` is the membrane interpolation and `u(t=1)-u(t=0) = ½ Σ a_i h_i (θ_i x V_n^i)`; `θ` parallel to `V_n` produces no director rotation (the 2017 core is blind to the drill component) |

### TDD evidence (WU2, explicit test-first; `strict_tdd: false`)

The 5 tests were written first and observed **RED** (unresolved API: `E0425` on
`interpolate_displacement`/`interpolate_position`/`Mitc4PlusDPrecomputed::new`). The
production code was then added and the suite went **GREEN** (135/0). Each test was
then shown to be able to fail by a deliberate perturbation of the implementation:

| # | Perturbation | Test shown RED | Observed failure |
| --- | --- | --- | --- |
| 1 | `x_d` wrong sign pattern in `compute_characteristic_vectors` | `test_geometry_flat_rectangle_...` | `x_d = [-1,0,0]`, expected `0` |
| 2 | `m_r`/`m_s` swapped in the dual-basis return | `test_geometry_dual_basis_identities_eq11` | `flat-distorted: m^r.x_r = 0` |
| 3 | `a_e = -2 c_r c_s / d` (the deleted sign) | `test_geometry_a_E_is_positive_eq27c` | `a_E must be positive, got -4.5714...` |
| 4 | director sub-normal sign alignment flipped | `test_geometry_node_directors_...` | director returned `-n_vec` |
| 5 | `theta.cross(&vn)` replaced by `theta` in `interpolate_displacement` | `test_kinematics_displacement_field_...` | `u(t=1)-u(t=0) != half the director rotation` |
| 6 | director term dropped from `interpolate_position` | `test_kinematics_displacement_field_...` | Eq. (1) position mismatch at `t = ±1` |

After restoring every perturbation the suite is **135 passed / 0 failed**.

### Deviations from design (WU2)

1. **3.1 module layout (recorded in the WU1 file header and on the task).** No separate
   `mitc4_plusd/tests.rs`; the tests live inline, per the repository convention adopted
   in WU1. The task's verification substring `mitc4_plusd::tests::fixtures` is therefore
   `mitc4_plusd::tests::`.
2. **`test_geometry_a_E_is_positive_eq27c` fixture.** The task and design ask for
   `a_E > 0` "for an element with `c_r c_s > 0`". Since `a_E = +2 c_r c_s / d`, that
   requires `d > 0`, i.e. `c_r^2 + c_s^2 > 1`, which no convex quad of the repository's
   fixtures reaches (a brute-force search over convex 2D integer quads found none; it
   needs an extreme 3D warped distortion). The test therefore exercises the closed form
   directly on its own inputs (`x_d`, `m_r`, `m_s` with `c_r c_s > 0`, `d > 0`) and
   re-checks the same closed form on a real flat distorted element, rejecting the
   negated form both times. The sign error the test exists to catch is fully covered.
3. **Constructor signature.** `Mitc4PlusDPrecomputed::new(node_coords, constitutive,
   thickness)` for WU2. The ADR-1 `applied_shear_correction` argument and the stored
   `cs_uncorrected` are **not** added here; they land in WU4 (task 5.2). This is the
   design's own slicing, and WU4 owns the call-site updates in this same file.
4. **WU2 size.** 1013 changed lines vs the design's ~300 forecast, above the ~350
   per-unit guide and the 700 session budget. Covered by the accepted session
   `size:exception`; the overrun is the paper/equation doc comments required by the task
   (every stored quantity and helper) and the 6-perturbation RED evidence. No test, doc
   or citation was dropped.
5. **Staged `dead_code` warnings.** Three helpers used only by WU3
   (`regularized_inverse_2x2`, `covariant_to_local_mapping`, `shear_covariant_to_local`)
   and the two kinematics helpers (used only by the tests until WU5) warn as "never
   used" in the lib build. This matches the existing `mitc4.rs`/`quad.rs` state, which
   already carries such staged-helper warnings; no `#[allow(dead_code)]` was added.
   `cargo test` is unaffected.

### Remaining tasks (WU2)

All WU2 tasks 3.1–3.4 are complete. WU3 (task 4.1 onwards) is next. The section-2
unchecked line is unchanged from the WU1 record (task 2.4, deferred to WU4).

### Workload / PR boundary (WU2)

Single PR, accepted `size:exception`. WU2 is a self-contained review slice: it extends
`mitc4_plusd.rs` with the geometry/coefficient/director/kinematics core and its five
tests, with a green tree, its own verification, and rollback = delete the WU2 block
(the WU1 fixtures and `mod.rs` line are untouched).

---

## WU3 — the B-operators (tasks 4.1, 4.2, 4.3, 4.4, 4.5)

Work unit: **WU3**. Appended to the WU1/WU2 record above; earlier bytes preserved.

### Structured status consumed (WU3)

- Source: native SDD status engine (authoritative, `artifactStore: openspec`).
- `changeName`: `mitc4plusd-faithful`; `applyState`: `ready`; `nextRecommended`: `apply`.
- `actionContext.mode`: `repo-local`; `workspaceRoot`:
  `/home/efirvida/Desktop/dev/fem-shell`; `allowedEditRoots`:
  `/home/efirvida/Desktop/dev/fem-shell`. All edits stayed inside the workspace and
  inside WU3's authorized edit roots (`mitc4_plusd.rs`, `tasks.md`,
  `apply-progress.md`).
- `taskProgress` at entry: 59 total / 12 completed / 47 pending.
- Review workload gate: `Decision needed before apply: No`, `Chained PRs recommended:
  No`, `Chain strategy: size-exception`, `400-line budget risk: High`. The session
  resolved delivery as **single-pr with an explicitly accepted `size:exception`** and a
  700-line review budget, so WU3 proceeded.

### Completed tasks (persisted checkboxes updated in `tasks.md`)

- [x] **4.1** the five covariant membrane tying rows + `b_membrane_2017` (Eq. 27).
- [x] **4.2** `b_bending_2017` → `(B_b1, B_b2)` (Eqs. 7c/7d + Eq. 8a).
- [x] **4.3** `b_shear_mitc4` + the four stored tying operators +
  `shear_covariant_to_local`.
- [x] **4.4** `b_drill_membrane_2025` (Eq. 18) + `drill_midside_shape_derivatives` +
  `tests::drill::b_md_reference` + the five asserted rejections.
- [x] **4.5** the `c_r`/`c_s` collision test.

### Files changed (WU3)

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | extended: the WU3 production operators (`j_loc_at`, `local_components`, `covariant_membrane_b_row`, `b_membrane_covariant_2017` / `b_membrane_2017`, `b_bending_covariant_2017` / `b_bending_2017`, `compute_shear_tie` / `b_shear_mitc4`, `DrillEdgeTerm` / `compute_drill_edges`, `drill_midside_shape_derivatives`, `drill_jacobian_ratio`, `b_drill_membrane_2025`), the new `Mitc4PlusDPrecomputed` fields (`b_rr_a`…`b_rs_e`, `b_shear_tie`, `drill_edges`) and the five WU3 tests in the inline test module (with `mod drill`) |
| `openspec/changes/mitc4plusd-faithful/tasks.md` | 4.1–4.5 checked; the 4.4 module-layout deviation recorded |
| `openspec/changes/mitc4plusd-faithful/apply-progress.md` | this WU3 section |

Diff stat (tracked file, vs the WU2 commit `2f6cdf0`):

```text
crates/aeroelast-core/src/elements/mitc4_plusd.rs | 1180 ++++++++++++++++++++-
1 file changed, 1176 insertions(+), 4 deletions(-)
```

`git diff --stat crates/aeroelast-core/src/elements/mitc4.rs` is **empty**: the hybrid
is byte-identical. No file outside the authorized set was touched. `elements/mod.rs`
was not edited.

### Verification (WU3)

Command (workspace root is `crates/`):

```text
cd crates && cargo test -p aeroelast-core
```

Result: **140 passed / 0 failed** (the WU2 baseline 135 + the 5 new WU3 tests).
`rustfmt --edition 2021 --check` is clean for the file. `cargo clippy -p
aeroelast-core --all-targets` reports no lint on the WU3 code other than the same
staged `dead_code` "never used" warnings WU2 already documented (the operators are
exercised by the tests and are wired into the stiffness in WU4/WU5). No
`#[allow(dead_code)]` was added.

### New tests and what each asserts

| Test | Asserts |
| --- | --- |
| `test_t1a_membrane_eq22_flat_tying_condition` | Eq. (22) `ẽ_rs^m\|bil = e_rs^m\|bil = x_d·u_d` to `1e-14` absolute on the flat distorted quad (non-vacuous, `\|bil\| > 1e-6`); the same comparison separates by `> 1e-6` relative on the doubly warped quad; the mapped operator equals the covariant field with its third row doubled followed by the point-wise mapping (proves the `2 e_rs` factor sits at the mapping, not in Eq. 27c); and Eq. (27c) reduces to Eq. (18c)'s `e_rs(E) + ½e_rr\|lin r + ½e_ss\|lin s` on a flat rectangle (note F2 / the leading `1`). |
| `test_identity_bending_operator_matches_eq7c_eq7d` | `b_bending_covariant_2017 == bending_reference` (Eqs. 7c/7d including `∂x_b·∂u_m` of Eq. 8a) to `1e-10` relative at the four Gauss points on the flat rectangle and the doubly warped quad; both operators are 3×24 with no condensed internal DOF; the `∂x_b·∂u_m` term is present (separates by `> 1e-6` on the warped quad) and vanishes `≤ 1e-14` on the flat one. |
| `test_shear_mitc4_flat_reduces_to_mindlin_assumed_field` | at the four DB84 tying points the local operator equals the standard Mindlin shears `γ_13 = w_,x + θ_y`, `γ_23 = w_,y − θ_x` to `1e-12`; non-vacuity: at a Gauss point the assumed field differs from the point-wise Mindlin field by `> 1e-6` relative. |
| `test_identity_drill_operator_matches_eq18_term_by_term` | **Oracle 1**: `b_drill_membrane_2025` vs `tests::drill::b_md_reference` entry by entry (72 entries) to `1e-12` absolute at 9 points (4 Gauss + 4 edge mid-points + centre) on flat square, flat distorted, ruled warped and doubly warped quads. **Oracle 2**: the five wrong variants each separate from Eq. (18) above their margin (`V^D → e3` `>1e-6`; missing `1/8` `>1e-3`; edge order `[bottom,right,top,left]` `>1e-6`; flipped edge difference `>1e-6`; `θ_z` alone `>1e-6`). |
| `test_identity_cr_cs_2017_and_2025_are_different_quantities` | `pre.c_r_mem`/`c_s_mem` equal `x_d·m^r`/`x_d·m^s` (Eq. 25); the stored `drill_edges[e].c_r/c_s` equal the independently recomputed Eq. (18)/(19c) values; and the 2017 and 2025 quantities separate by `> 1e-6` relative on the warped quad. |

### TDD evidence (WU3, explicit test-first; `strict_tdd: false`)

The five tests were written first and observed **RED**: `cargo test` failed to compile
with `error[E0432]: unresolved imports super::b_bending_2017, …, super::drill_edges …`
and `error[E0609]: no field drill_edges on type Mitc4PlusDPrecomputed`. The production
operators were then added and the suite went **GREEN** (140/0). Each test was
afterwards shown to fail for a wrong implementation; restoring the file returns 140/0.

| # | Test shown RED | Perturbation | Observed failure |
| --- | --- | --- | --- |
| 1 | `test_t1a_membrane_eq22_flat_tying_condition` | drop the Eq. (27c) leading `1` (`(a_E rs)` instead of `(1 + a_E rs)`) | flat-rectangle Eq. (18c) reduction mismatch |
| 2 | same | negate `a_E` in Eq. (27c) | flat bilinear-coefficient mismatch |
| 3 | `test_identity_bending_operator_matches_eq7c_eq7d` | zero the `∂x_b·∂u_m` term of Eq. (8a) | B_b1 mismatch vs Eq. (7c) at a Gauss point |
| 4 | `test_shear_mitc4_flat_reduces_to_mindlin_assumed_field` | swap the A/B tying rows | tying-point Mindlin mismatch |
| 5 | `test_identity_drill_operator_matches_eq18_term_by_term` | `V^D → e3` (both `c_r/c_s` and `θ^D`) | `ruled-warped [0][4]`: `0` vs `−1.88e-4` |
| 6 | same | drop the `1/8` of Eq. (13c) | `flat-square [2][11]`: `−2.309` vs `−0.2887` |
| 7 | same | `h̃` edge order `[bottom, right, top, left]` | `flat-square [0][5]`: `−6.10e-2` vs `0` |
| 8 | same | flip the edge-difference sign | `flat-square [2][11]`: `+0.2887` vs `−0.2887` |
| 9 | same | `θ_z` alone in place of `θ·V^D` | `ruled-warped [0][4]`: `0` vs `−1.88e-4` |
| 10 | `test_identity_cr_cs_2017_and_2025_are_different_quantities` | set the stored `drill_edges` c-values to `0.0` | stored `c_r` mismatch |
| 11 | `test_identity_drill_operator_matches_eq18_term_by_term` (the 4.5 conflation) | substitute the 2017 `c_r_mem`/`c_s_mem` into the drill operator | `flat-square [2][11]`: `0` vs `−0.2887` |

Rows 5–9 are the design's five asserted rejections; the in-test variant assertions and
these production perturbations are the two sides of the same evidence.

### Findings (WU3)

1. **The Eq. (22) test needed one extra assertion to cover note F2.** The task's
   specified comparison extracts the bilinear (`r·s`) coefficient; a four-corner second
   difference is blind to a constant offset, so dropping Eq. (27c)'s leading `1` does
   not move that coefficient. The test now additionally asserts the flat-rectangle
   reduction to Eq. (18c), which the leading `1` *is* required for — that assertion is
   what row 1 above shows RED. No extract defect is implied; the extract's Eq. (18) and
   note F2 agree with the implementation.
2. **No error found in either extract.** The production Eq. (18) operator agrees with
   the independently written `b_md_reference` to `1e-12` at every sample point on all
   four quads, and Eq. (22) holds exactly on flat geometry (also re-derived
   algebraically: `a_A e_rr(A) + … + a_E e_rs(E) = x_d·u_d` when `x_d·n = 0`). No vision
   re-read was needed and neither extract was edited.
3. **A test-side frame bug was found and fixed during test-first.** The first `true_bil`
   applied the local-frame projection twice (the DOFs are already local components); the
   failure exposed it and the test now projects `x_d` only. The production operator was
   correct throughout.

### Deviations from design (WU3)

1. **4.4 module layout (recorded on the task).** Tasks name `…/tests/drill.rs`; per the
   WU1 repository-convention decision the independent reference lives in `mod drill`
   inside the inline `#[cfg(test)] mod tests`. It still shares no code with
   `b_drill_membrane_2025`.
2. **WU3 size.** ~1180 changed lines vs the design's ~350 forecast and the 700-line
   session budget. Covered by the accepted session `size:exception`. The overrun is the
   paper/equation doc comments on every operator plus the independent drill reference,
   its one-parameter wrong variants, the 4-quad × 9-point × 72-entry comparison and the
   11-row RED evidence. No test, doc or citation was dropped.
3. **`b_membrane_covariant_2017` / `b_bending_covariant_2017` helpers.** The design lists
   only `b_membrane_2017` and `b_bending_2017 (→ B_b1, B_b2)`. The covariant halves are
   split out so the Eq. (27) / Eqs. (7c)-(7d) assembly can be tested independently of
   the covariant→local mapping, and so the `2 e_rs` factor's placement is explicit. The
   designed entry points are unchanged and delegate to them.

### Remaining tasks (WU3)

The next implementable unit is **WU4** (tasks 5.1–5.8); the exact remaining unchecked
lines of section 5 are:

```text
- [ ] **5.1 `resultant_moment_matrix`** (the `W_00 … W_22` block matrix of design §2.3, with `W_22 = cm/9` from paper A's 2×2 `t`-rule), with the multi-ply approximation documented in the doc comment. Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: a unit test of the moment entries against their closed forms; no spec-fixed name exists, so this is stated as a design-derived unit test rather than a spec test. Satisfies: Requirement 1 item 8. <!-- sdd-owner: implementation -->
- [ ] **5.2 `compute_ke_local` / `compute_ke_global` + the drill contribution + the `cs_uncorrected` wiring + the test-local reference implementation.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_ke_lock_matches_2017_core_plus_2025_drill` — `max|K_prod − K_ref| ≤ 1e-10·max|K_ref|` on flat square, flat distorted, ruled warped and doubly warped quads, with the same bound on the membrane and transverse-shear blocks. Satisfies: Requirement 1 scenario 1. <!-- sdd-owner: implementation -->
- [ ] **5.3 Drill stiffness provenance.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_drill_stiffness_comes_only_from_eq26` — non-zero on warped geometry; exactly symmetric; **exactly zero on every translational row/column block**; exactly zero on every rotation block other than the drill's own on flat geometry; and `|u_rbᵀ K u_rb| ≤ 1e-12·λ_max·‖u_rb‖²` for the six rigid-body fields. Satisfies: Requirement 1 scenario 2 (design §9.1 interpretation, already reflected in spec rev 4). <!-- sdd-owner: implementation -->
- [ ] **5.4 Uncorrected transverse shear — the discriminating test.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_uses_uncorrected_shear_modulus` (block vs closed-form `∫ B_γᵀ (G·I) B_γ dA` to `1e-10`; the `5/6` value rejected by `> 1e-3` relative, asserted) and `cd crates && cargo test -p aeroelast-core test_identity_transverse_shear_invariant_to_shear_correction_factor` (ADR-1's two-`pre` construction on isotropic and single-ply laminate, plus the two non-vacuity controls: the `k`-carrying control differs by `> 1e-3` relative). Satisfies: Requirement 16. <!-- sdd-owner: implementation -->
- [ ] **5.5 Integration rule.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` on the strongly warped quad — matches the three-term reference to `1e-10` and differs from the surface-only reference by `> 1e-4` relative (asserted). Satisfies: Requirement 1 item 8. <!-- sdd-owner: implementation -->
- [ ] **5.6 Local matrix shapes.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_local_matrices_are_24x24` — `K`, `M`, `K_T` exactly 24×24 and `f_int` exactly 24 long, drilling at slot `6i+5`. Satisfies: Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **5.7 Drill-DOF energy behaviour.** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_kinematics_drill_dof_is_theta_z_through_eq26_operator` — energy `≤ 1e-12·λ_max·‖u‖²` for a pure rigid rotation about `V_n`, and non-zero energy on the warped `θ_z` pattern only through the Eq. (26) operator. Satisfies: Requirement 2. <!-- sdd-owner: implementation -->
- [ ] **5.8 Mid-surface restriction (ADR-6 / G7).** Touches `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`. Verification: `cd crates && cargo test -p aeroelast-core test_identity_element_uses_midsurface_constitutive` — the element's constitutive equals `Laminate::to_shell_constitutive()` — plus the supporting static check `grep -n "to_shell_constitutive_with_offset" crates/aeroelast-core/src/elements/mitc4_plusd.rs` returns no match (supporting, not normative). Satisfies: Requirement 16 scenario 2. <!-- sdd-owner: implementation -->
```

Sections 6–13 remain pending in `tasks.md`; task 2.4 remains the recorded WU4 deferral.

### Workload / PR boundary (WU3)

Single PR, accepted `size:exception`. WU3 is one review slice: the four B-operators and
the drill verification, extending `mitc4_plusd.rs` with a green tree, its own
verification `command` and rollback = revert the WU3 block (WU1/WU2 and `mod.rs` are
untouched).

---

## WU4 (constitutive half) — the ADR-1 uncorrected transverse-shear mechanism (task 5.4's material half)

Work unit: the **constitutive half of WU4**, deliberately split out. This unit lands only
the ADR-1 material-side mechanism (`materials/mod.rs`, `materials/laminate.rs`) and its
own unit tests. The stiffness assembly — `resultant_moment_matrix`, `compute_ke_local` /
`compute_ke_global`, the drill contribution, the `cs_uncorrected` element wiring and the
element-level discriminating tests (tasks 5.1–5.3, 5.5–5.8, and task 5.4's own
verification commands) — is **not** in scope here; it is the next unit. Appended to the
WU1/WU2/WU3 record above; earlier bytes are preserved.

### Structured status consumed (WU4-constitutive)

- Source: native SDD status engine (authoritative, `artifactStore: openspec`).
- `changeName`: `mitc4plusd-faithful`; `applyState`: `ready`; `nextRecommended`: `apply`.
- `actionContext.mode`: `repo-local`; `workspaceRoot`:
  `/home/efirvida/Desktop/dev/fem-shell`; `allowedEditRoots`:
  `/home/efirvida/Desktop/dev/fem-shell`. All edits stayed inside the workspace and
  inside this unit's authorized edit roots (`materials/mod.rs`, `materials/laminate.rs`,
  `tasks.md`, `apply-progress.md`).
- `taskProgress` at entry: 59 total / 17 completed / 42 pending. It is **unchanged** on
  exit: task 5.4's own verification is element-level and stays unchecked; the mechanism
  is recorded in a note on that task rather than as a completed checkbox.
- Review workload gate: `Decision needed before apply: No`, `Chained PRs recommended:
  No`, `Chain strategy: size-exception`, `400-line budget risk: High`. The session
  resolved delivery as **single-pr with an explicitly accepted `size:exception`** and a
  700-line review budget, so this unit proceeded.

### What landed (additive only)

- `ShellConstitutive::transverse_shear_uncorrected(&self, applied_k: f64) -> Matrix2<f64>`
  — returns `cs / applied_k` when `applied_k > 0`, else `cs` unchanged. `ShellConstitutive`
  keeps **exactly its five fields**; no sixth field was added.
- `Laminate::applied_shear_correction_factor(&self) -> f64` — returns
  `shear_correction_factor` for a single-ply laminate, `1.0` for a multi-ply one (the
  factor `compute_shear_stiffness` actually applied).

This is the same material channel that task 10.1 names; it landed early here because the
uncorrected-shear mechanism is what task 5.4 is about. A cross-reference note was added
on task 10.1; its wiring (`MaterialSpec::Composite` field, PyO3 call sites) remains for
WU9.

### Files changed (WU4-constitutive)

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/materials/mod.rs` | **+126 lines**: the `transverse_shear_uncorrected` accessor (+ doc comment) and a new `#[cfg(test)] mod tests` with 5 tests |
| `crates/aeroelast-core/src/materials/laminate.rs` | **+42 lines**: the `applied_shear_correction_factor` accessor (+ doc comment) and 1 test in the existing test module |
| `openspec/changes/mitc4plusd-faithful/tasks.md` | the task 5.4 mechanism note and the task 10.1 cross-reference note |
| `openspec/changes/mitc4plusd-faithful/apply-progress.md` | this WU4-constitutive section |

Diff stat (the two Rust files):

```text
crates/aeroelast-core/src/materials/laminate.rs |  42 ++++++++
crates/aeroelast-core/src/materials/mod.rs      | 126 ++++++++++++++++++++++++
2 files changed, 168 insertions(+)
```

`git diff --numstat` reports **168 added / 0 deleted** for both files: every added line is
new, no pre-existing line was modified or removed. `git diff --stat
crates/aeroelast-core/src/elements/mitc4.rs` is **empty** (the hybrid is byte-identical),
and no file outside the authorized set was touched (`git status --short` shows only the
two materials files modified plus the pre-existing untracked `.pi/`).

### Verification (WU4-constitutive)

Command (workspace root is `crates/`):

```text
cd crates && cargo test -p aeroelast-core
```

Result: **146 passed / 0 failed** (baseline 140 + the 6 new tests). Focused run
`cargo test -p aeroelast-core materials::` → **18 passed / 0 failed**. `cargo clippy -p
aeroelast-core --all-targets` reports **no lint on the new lines** (the pre-existing
`op_ref` lints in the older laminate tests were left untouched). `rustfmt --check` on
these files is not clean, but **that is pre-existing** (the materials files are not
rustfmt-formatted at HEAD, e.g. `composite.rs`/`failure.rs`); the new code follows the
file's existing compact style and the diff stays purely additive — reformatting the whole
file would have violated the additive-only constraint.

### New tests and what each asserts

| Test | Asserts |
| --- | --- |
| `materials::tests::test_transverse_shear_uncorrected_single_isotropic_ply_equals_g_h` | a single isotropic ply's `cs` is exactly `k·G·h` (`≤ 1e-12` rel), the reported applied factor is `k`, and `transverse_shear_uncorrected(k)` recovers the uncorrected `G·h` (`≤ 1e-12` rel) |
| `materials::tests::test_transverse_shear_uncorrected_isotropic_constitutive_removes_k` | the isotropic `ShellConstitutive` (`cs = k·G·h`) yields `G·h` after removing `k` (`≤ 1e-12` rel) |
| `materials::tests::test_transverse_shear_uncorrected_multi_ply_is_cs_unchanged` | a multi-ply laminate reports applied factor `1.0` and `transverse_shear_uncorrected(1.0)` equals `cs` unchanged (`≤ 1e-12` rel) — i.e. it does **not** divide by `k` there |
| `materials::tests::test_transverse_shear_uncorrected_discriminates_naive_multi_ply_division` | the rejected mechanism (c) — divide the multi-ply `cs` by the laminate's scalar `shear_correction_factor` — differs from the correct value by **`> 1e-3` relative** (asserted, so the test cannot pass vacuously) |
| `materials::tests::test_transverse_shear_uncorrected_nonpositive_factor_returns_cs` | `applied_k = 0` and `applied_k < 0` both return `cs` unchanged, with no division by zero |
| `materials::laminate::tests::test_applied_shear_correction_factor_single_and_multi_ply` | the accessor returns `shear_correction_factor` (`5/6`) for a single-ply laminate and exactly `1.0` for a three-ply one |

### TDD evidence (explicit test-first; `strict_tdd: false`)

**RED (all six).** The tests were written first, against accessors that did not exist.
`cd crates && cargo test -p aeroelast-core materials::` failed to compile with **12
errors**, every one of the form:

```text
error[E0599]: no method named `applied_shear_correction_factor` found for struct `laminate::Laminate`
error[E0599]: no method named `transverse_shear_uncorrected` found for struct `ShellConstitutive`
```

The two accessors were then added and the suite went **GREEN** (146/0). Each test was
afterwards shown to fail for a deliberately wrong implementation; restoring the file
returns 146/0.

| # | Perturbation | Tests shown RED | Observed failure |
| --- | --- | --- | --- |
| 1 | `applied_shear_correction_factor` returns `self.shear_correction_factor` always (the unsound mechanism) | `test_applied_shear_correction_factor_single_and_multi_ply`; `test_transverse_shear_uncorrected_multi_ply_is_cs_unchanged`; `test_transverse_shear_uncorrected_discriminates_naive_multi_ply_division` | multi-ply returned `0.8333…` vs `1.0` (×2); the discriminator's `rel = 0` (`> 1e-3` failed), proving the discriminator catches the unsound mechanism |
| 2 | `transverse_shear_uncorrected` divides unconditionally (drops the `applied_k > 0` guard) | `test_transverse_shear_uncorrected_nonpositive_factor_returns_cs` | `left: [[inf, NaN], [NaN, inf]]` vs `cs` — the divide-by-zero is caught |
| 3 | `transverse_shear_uncorrected` multiplies by `applied_k` instead of dividing | `test_transverse_shear_uncorrected_single_isotropic_ply_equals_g_h`; `test_transverse_shear_uncorrected_isotropic_constitutive_removes_k` | uncorrected returned `k²·G·h` instead of `G·h` |

### Preserved-laminate invariant (hard user constraint)

- **Additive proof.** `git diff --numstat` for `materials/mod.rs` + `materials/laminate.rs`
is `168  0` (and `168  0` per file): no pre-existing line was changed or removed. The two
accessors are new `pub` methods; `ShellConstitutive`'s five fields, `Laminate`'s fields and
`compute_shear_stiffness` are untouched.
- **Existing tests unchanged and passing.** The 7 pre-existing laminate tests
(`test_laminate_symmetric`, `test_laminate_to_shell_constitutive`, `test_z_offset_zero`,
`test_z_offset_transforms_b_d`, `test_asymmetric_laminate_has_nonzero_b`, plus the
`materials::failure` set) are part of the 146/0 run with their original assertions and
no edit to any expectation or tolerance. `git diff` contains no modification to them.
- **No other material file touched.** `orthotropic.rs`, `composite.rs`, `failure.rs`,
`isotropic.rs` are byte-identical.

### Deviations from design (WU4-constitutive)

1. **Scope split (user-directed, recorded).** The material-side mechanism is nominally
task 10.1's code; it was requested here as "task 5.4's mechanism" and landed early. This
is not a design deviation — the mechanism is exactly ADR-1's chosen accessors — only a
sequencing note, and it is cross-referenced on both tasks 5.4 and 10.1.
2. **No design deviation in the accessor bodies.** `transverse_shear_uncorrected` is
byte-for-byte the ADR-1 snippet (`if applied_k > 0.0 { self.cs / applied_k } else {
self.cs }`), and `applied_shear_correction_factor` is the ADR-1 snippet verbatim.
3. **`rustfmt` not applied to the files** (see Verification): pre-existing non-clean
state; applying it would break the additive-only diff. No line of new code depends on it.

### Remaining tasks (WU4-constitutive)

The next implementable unit is the **stiffness assembly** (tasks 5.1–5.3, 5.5–5.8 and the
element-level half of 5.4). The exact remaining unchecked lines of section 5:

```text
- [ ] **5.1 `resultant_moment_matrix`** … <!-- sdd-owner: implementation -->
- [ ] **5.2 `compute_ke_local` / `compute_ke_global` + the drill contribution + the `cs_uncorrected` wiring + the test-local reference implementation.** … <!-- sdd-owner: implementation -->
- [ ] **5.3 Drill stiffness provenance.** … <!-- sdd-owner: implementation -->
- [ ] **5.4 Uncorrected transverse shear — the discriminating test.** … <!-- sdd-owner: implementation -->
  - **Mechanism done (constitutive half, this unit); element-level tests deferred to the assembly unit.** …
- [ ] **5.5 Integration rule.** … <!-- sdd-owner: implementation -->
- [ ] **5.6 Local matrix shapes.** … <!-- sdd-owner: implementation -->
- [ ] **5.7 Drill-DOF energy behaviour.** … <!-- sdd-owner: implementation -->
- [ ] **5.8 Mid-surface restriction (ADR-6 / G7).** … <!-- sdd-owner: implementation -->
```

(The lines are elided with `…` for readability; the full text is in `tasks.md` at section
5, and each line's bytes are unchanged apart from the note appended under 5.4.) Task 2.4
remains the recorded WU4 deferral. Sections 6–13 remain pending.

### Workload / PR boundary (WU4-constitutive)

Single PR, accepted `size:exception`. This unit is one small, self-contained review slice:
168 additive lines across the two material files, a green tree, its own verification and
rollback = revert the two accessors and their tests (the WU1–WU3 bytes are untouched).

---

## WU4 (stiffness assembly) — tasks 5.1–5.8, the element-level half of 5.4, and the frame-convention fix

Work unit: **the WU4 stiffness assembly** — `resultant_moment_matrix`, `compute_ke_local` /
`compute_ke_global`, the drill contribution, the `cs_uncorrected` wiring, the test-local
independent reference, and the element-level tests 5.1–5.8. Appended to the WU1–WU3 and
WU4-constitutive record above; earlier bytes are preserved.

### Structured status consumed (WU4-assembly)

- Source: native SDD status engine (authoritative, `artifactStore: openspec`).
- `changeName`: `mitc4plusd-faithful`; `applyState`: `ready`; `nextRecommended`: `apply`.
- `actionContext.mode`: `repo-local`; `workspaceRoot`:
  `/home/efirvida/Desktop/dev/fem-shell`; `allowedEditRoots`:
  `/home/efirvida/Desktop/dev/fem-shell`. All edits stayed inside the workspace and
  inside this unit's authorized edit roots (`crates/aeroelast-core/src/elements/mitc4_plusd.rs`,
  `openspec/changes/mitc4plusd-faithful/tasks.md`,
  `openspec/changes/mitc4plusd-faithful/apply-progress.md`).
- `taskProgress` at entry: 59 total / 17 completed / 42 pending.
- Review workload gate: `Decision needed before apply: No`, `Chained PRs recommended:
  No`, `Chain strategy: size-exception`, `400-line budget risk: High`. The session
  resolved delivery as **single-pr with an explicitly accepted `size:exception`** and a
  700-line review budget, so this unit proceeded.

### Completed tasks (persisted checkboxes updated in `tasks.md`)

- [x] **5.1** `resultant_moment_matrix` (the `W_00 … W_22` block matrix; `W_22 = cm/9`).
- [x] **5.2** `compute_ke_local` / `compute_ke_global`, the drill contribution, the
  `cs_uncorrected` wiring, and the independent test-local reference.
- [x] **5.3** drill-stiffness provenance.
- [x] **5.4** the element-level half: `test_identity_transverse_shear_uses_uncorrected_shear_modulus`
  and `test_identity_transverse_shear_invariant_to_shear_correction_factor`.
- [x] **5.5** integration rule (2×2×2 vs surface-only).
- [x] **5.6** local matrix shapes (`K` 24×24 and the drill-slot layout; `M`/`K_T`/`f_int`
  deferred to WU5, recorded on the task).
- [x] **5.7** drill-DOF energy behaviour.
- [x] **5.8** mid-surface restriction (ADR-6 / G7).

### Files changed (WU4-assembly)

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | extended: the WU4 production assembly (`resultant_moment_matrix`, `surface_measure`, `membrane_ke_local`, `shear_ke_local`, `drill_ke_local`, `compute_ke_local_with_drill`, `compute_ke_local`, `compute_ke_global`, `build_t24`, `transform_to_global`), the ADR-1 fields `applied_shear_correction` / `cs_uncorrected` and the constructor argument, the independent `ke_ref` reference, the 9 new tests, and the WU2/WU3 frame-convention fix |
| `openspec/changes/mitc4plusd-faithful/tasks.md` | 5.1–5.8 checked; notes under 5.2, 5.3 (the frame fix), 5.4, 5.6, 5.8 |
| `openspec/changes/mitc4plusd-faithful/apply-progress.md` | this WU4-assembly section |

Diff stat (tracked file, vs the WU4-constitutive commit `ee13216`):

```text
crates/aeroelast-core/src/elements/mitc4_plusd.rs | 1278 ++++++++++++++++++++-
1 file changed, 1278 insertions(+), 27 deletions(-)
```

`git diff --stat crates/aeroelast-core/src/elements/mitc4.rs` is **empty** and
`git diff --numstat` reports nothing for it: the hybrid is byte-identical. No file outside
the authorized set was touched (`git status --short` shows only `mitc4_plusd.rs` modified
plus the pre-existing untracked `.pi/`). `elements/mod.rs` was not edited.

### Verification (WU4-assembly)

Command (workspace root is `crates/`):

```text
cd crates && cargo test -p aeroelast-core
```

Result: **155 passed / 0 failed** (the WU4-constitutive baseline 146 + the 9 new tests).
`rustfmt --edition 2021 --check` is clean for the file. `cargo clippy -p aeroelast-core
--all-targets` reports no lint on the new lines other than the staged
`membrane_ke_local` "never used" in the lib build (the same staged-helper situation WU2
and WU3 documented: it is exercised by the identity lock in the test build and is wired
into `compute_ke_local`; no `#[allow(dead_code)]` was added).

### New tests and what each asserts

| Test | Asserts |
| --- | --- |
| `test_identity_resultant_moment_matrix_blocks_match_closed_forms` | `W` is 9×9, exactly symmetric, with `W_00 = cm`, `W_01 = cb_coupling` (`W_10 = cb_couplingᵀ`), `W_02 = cb` (`W_20 = cbᵀ`), `W_11 = cb`, `W_12 = 0`, `W_22 = cm/9`; the isotropic closed forms `cm00 = E h/(1−ν²)`, `cb00 = E h³/(12(1−ν²))`, `W_22 = cm00/9`; the paper's 2-point `t`-rule `Σ w_i t_i⁴ = 2/9`; and `W_22` is **not** the exact `cm/5` (non-vacuity) |
| `test_identity_ke_lock_matches_2017_core_plus_2025_drill` | production `K` vs the independent `ke_ref::ke_local` to `≤ 1e-10·max|K_ref|` on flat square, flat distorted, ruled warped and doubly warped quads, with the same bound on the membrane block and the transverse-shear block |
| `test_identity_drill_stiffness_comes_only_from_eq26` | `K(op) − K(op := 0)` is non-zero on warped geometry (the operator is live), **exactly symmetric**, **exactly zero on every translational row/column block**, exactly zero on every rotation block other than the drill's own on flat geometry (where `V^D = e3`), the flat drill block is non-zero (non-vacuity), and the six rigid-body fields carry `|u_rbᵀ K u_rb| ≤ 1e-12·λ_max·‖u_rb‖²` |
| `test_identity_transverse_shear_uses_uncorrected_shear_modulus` | the element's transverse-shear block equals the closed form `Σ_g B_γᵀ (G·h·I) B_γ w √g` to `1e-10` relative, and the `5/6` value is rejected by `> 1e-3` relative (asserted) |
| `test_identity_transverse_shear_invariant_to_shear_correction_factor` | two single-ply-laminate `pre` values with `k = 5/6` and `k = 0.5` give the same shear block to `1e-10`, both match the uncorrected closed form, and the `k`-carrying non-vacuity control (`applied_k = 1.0`) differs by `> 1e-3` relative; the same construction and controls on an isotropic constitutive |
| `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` | on the strongly warped quad, production matches the three-term (`W_22 = cm/9`) reference to `1e-10` and differs from the surface-only (`W_22 = 0`) reference by `> 1e-4` relative (asserted) |
| `test_kinematics_local_matrices_are_24x24` | `K` local and global are exactly 24×24 (576 entries); the 2017-only operators are exactly blind to slot `6i+5` on flat geometry; the drill block is live at `6i+5` (the `M`/`K_T`/`f_int` shapes are WU5, recorded on the task) |
| `test_kinematics_drill_dof_is_theta_z_through_eq26_operator` | a pure rigid rotation about `V_n` carries energy `≤ 1e-12·λ_max·‖u‖²`; a warped drill-rotation pattern `θ_i = γ_i V_n^i` carries zero energy through the 2017-only blocks and non-zero energy through the Eq. (26) operator |
| `test_identity_element_uses_midsurface_constitutive` | the element's constitutive equals `Laminate::to_shell_constitutive()` field by field and `applied_shear_correction` equals `Laminate::applied_shear_correction_factor()`; the offset coupling block `B − z₀A` differs (non-vacuity) |

### TDD evidence (explicit test-first; `strict_tdd: false`)

**RED (compile).** The 9 tests were written first, against production functions that did
not exist. `cd crates && cargo test -p aeroelast-core` failed to compile with:

```text
error[E0432]: unresolved imports `super::compute_ke_global`, `super::compute_ke_local`,
`super::compute_ke_local_with_drill`, `super::drill_ke_local`, `super::membrane_ke_local`,
`super::resultant_moment_matrix`, `super::shear_ke_local`, `super::surface_measure`
```

The production assembly was then added and the suite went **GREEN** (155/0). Each test was
afterwards shown to fail for a deliberately wrong implementation; restoring the file
returns 155/0.

| # | Test shown RED | Perturbation | Observed failure |
| --- | --- | --- | --- |
| 1 | `test_identity_resultant_moment_matrix_blocks_match_closed_forms` | `W_22 = cm/5` (exact `t`-integration) | `W_22 = cm/9` `left: 879120879.12` `right: 488400488.40` |
| 2 | `test_identity_ke_lock_matches_2017_core_plus_2025_drill` | `s1 = 2/h → 4/h` | full-K mismatch vs the reference |
| 3 | `test_identity_drill_stiffness_comes_only_from_eq26` | drill operator zeroed (`b_md := 0`) | the "operator is inert" assertion fails — the non-vacuity control |
| 4 | `test_identity_transverse_shear_uses_uncorrected_shear_modulus` | `cs_uncorrected → constitutive.cs` (residual `k`) | shear block ≠ the `G·h` closed form |
| 5 | `test_identity_transverse_shear_invariant_to_shear_correction_factor` | `cs_uncorrected → constitutive.cs` | the two `pre` shear blocks no longer agree |
| 6 | `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` | `W_22 = 0` (surface-only) | production no longer matches the three-term reference |
| 7 | `test_kinematics_local_matrices_are_24x24` | membrane operator writes slot `6i+5` | `membrane touches the drill slot` `left: 0.99999…` `right: 0.0` |
| 8 | `test_kinematics_drill_dof_is_theta_z_through_eq26_operator` | a `θ_z` diagonal penalty added to the drill block | `rigid rotation about V_n carries energy 1960000.000038147` |
| 9 | `test_identity_transverse_shear_uses_uncorrected_shear_modulus` (5.4a control) | the `5/6` candidate made equal to the correct value | `the 5/6 value must be rejected (relative 0)` |
| 10 | `test_identity_transverse_shear_invariant_to_shear_correction_factor` (5.4b control) | the control `pre` given the correct `applied_k` | `the k-carrying control must be separated (relative 0)` |
| 11 | `test_identity_element_uses_midsurface_constitutive` | the element path given the offset constitutive | `cb_coupling` mismatch (`left: [[-10989010.98, …]]`, `right: [[0.0, …]]`) |

Rows 3, 9 and 10 are the non-vacuity controls the task names for 5.3 and 5.4: each was
shown capable of failing.

### Finding: a WU2/WU3 frame-convention defect, exposed by 5.3 and fixed

**What was wrong.** The stored `pre.vn` / `pre.v1` / `pre.v2` and `pre.v_d` are
**global-frame** vectors (`compute_node_directors` builds them from the global node
coordinates and the global `e1`; the geometry test compares `pre.vn[i]` with the global
`n_vec`). Four operators paired them with the **local-frame** DOF triples:

- `b_bending_covariant_2017` — `xm_r`/`xb_r` local, `cross` built from `pre.vn` global;
- `compute_shear_tie` — `g_r`/`g_s`/`g_t` projected to local, `vxg_r = pre.vn × g_r`;
- `interpolate_displacement` — `theta` local, `theta.cross(&pre.vn[i])`;
- `b_drill_membrane_2025` — `theta` local, `theta^D = theta · pre.v_d`.

**How it showed.** `test_identity_drill_stiffness_comes_only_from_eq26` measured the six
rigid-body fields on warped geometry before the fix:

```text
flat-square     |u_rbᵀ K u_rb| / (λ_max ‖u_rb‖²) = 2.7e-17
flat-distorted                                     3.7e-17
ruled-warped                                       2.5e-3
doubly-warped                                      5.4e-3
```

The spec requires `≤ 1e-12`. The 2017 core computes `u_b = ½ Σ a_i h_i (θ_i × V_n^i)`,
which is only a rigid rotation when `θ_i` and `V_n^i` are expressed in the same basis;
mixing them leaves a spurious strain on warped geometry. (The identity lock of 5.2 did
**not** catch this, because its independent reference shared the same convention — which
is exactly why the rigid-body assertion exists as a second, independent oracle.)

**The fix.** The four operators now project the director to the local frame
(`local_components(pre, &pre.vn[i])`, `local_components(pre, &pre.v_d)`), and the three
matching test-local references (`bending_reference`, `drill::b_md_reference` /
`b_md_parameterised`, and the kinematics test's `u_b` reconstruction) were updated the same
way. `compute_j3d_enriched`, `interpolate_position`, `j_loc_at`, `surface_measure`,
`drill_jacobian_ratio` and the shear metric continue to use the global director (they work
in global coordinates), so the change is exactly a same-basis pairing at the local-DOF use
sites. After the fix the rigid-body ratio is within the spec's `1e-12` bound on all four
geometries and the full 155-test suite is green. This is a WU2/WU3 correctness fix
discovered by WU4 and recorded on task 5.3.

### Deviations from design (WU4-assembly)

1. **`build_t24` / `transform_to_global` landed here, not in WU5.** Task 5.2 requires
   `compute_ke_global`, which needs them; the design's WU5 list also names them. They are
   the retargeted repository transformation (`Tᵀ M T`), so WU5's remaining work is
   unaffected.
2. **`membrane_ke_local` / `shear_ke_local` / `drill_ke_local` helpers.** The design names
   only `resultant_moment_matrix` and `compute_ke_local`/`compute_ke_global`. The three
   blocks are split out so the identity lock can compare the membrane and transverse-shear
   blocks separately (the spec's "with the same bound on the membrane and transverse-shear
   blocks") and so the drill provenance test can form the exact
   `K(op) − K(op := 0)` difference. The designed entry points are unchanged and delegate
   to them.
3. **`compute_ke_local_with_drill(pre, use_drill)` is private.** It implements the spec's
   "variant reference in which the B Eq. (26) operator is replaced by the zero operator"
   without adding a second public API.
4. **End-of-assembly symmetrisation.** `compute_ke_local_with_drill` and
   `transform_to_global` symmetrise their result (`0.5 (M + Mᵀ)`), and `drill_ke_local`
   symmetrises the drill block. This is a round-off guard, not a formulation factor: the
   spec's drill-provenance scenario demands **exact** symmetry and exact zeros, and the
   assembled `Bᵀ W B` is only exactly symmetric in exact arithmetic. No value changes
   beyond the last bit.
5. **5.6's `M`/`K_T`/`f_int` shapes are deferred to WU5.** Those functions are tasks
   6.1/6.2 and do not exist in this unit; the WU4-owned part of the scenario (`K` local and
   global 24×24 and the drill-slot layout) is asserted. Recorded on task 5.6.
6. **5.8's static check.** The non-vacuity control recomputes the offset coupling block
   `B − z₀A` from its definition instead of calling `to_shell_constitutive_with_offset`, so
   the task's supporting grep over the module returns no match.
7. **WU4 size.** 1278 added / 27 removed lines vs the design's ~340 forecast, above the
   ~350 per-unit guide and the 700 session budget. Covered by the accepted session
   `size:exception`. The overrun is the paper/equation doc comments on every new function,
   the fully independent `ke_ref` reference (~250 lines), and the frame fix with its
   matching reference updates. No test, doc or citation was dropped.

### Findings on the design and the extracts

1. **Design §2.2's "the 2017 core is structurally blind to `θ_z`" is only true on flat
   geometry.** The 2017 core is blind to the component of `θ_i` **along `V_n^i`**; the
   local DOF `θ_z` is along `e3`, which equals `V_n^i` only when the element is flat. On
   warped geometry the 2017 bending/shear operators do couple to `θ_z` (WU4 measured it:
   `test_kinematics_local_matrices_are_24x24` asserted blindness and failed on the doubly
   warped quad before the assertion was scoped to flat geometry). The spec's own qualifier
   ("the flat case is what makes this clause exactly testable", rev 4) already anticipates
   this; the design's §2.2 sentence should be read with that qualifier.
2. **Design §2.1's `vn: // V_n^i, local frame` is wrong as implemented in WU2** — the
   stored vectors are global-frame (see the finding above). The doc comment and the field
   should either be corrected to "global frame" (with the operators projecting at the
   local-DOF use sites, as now) or `compute_node_directors` should be changed to store
   local vectors (which would also require `compute_j3d_enriched` / `interpolate_position`
   to convert back). WU4 took the first route; the design's intent (operators in the local
   frame) is preserved either way.
3. **Design §2.3's "`sqrt_g = j(r,s)`" is loose.** The integration measure is the
   mid-surface area measure `‖g_r × g_s‖`, not the drill Jacobian `j = det[g_r g_s g_t]`;
   using the triple product would scale every block by `h/2` and contradict the ABD
   resultant formulation. WU4 implemented `‖g_r × g_s‖` (the repository's existing
   `sqrt_g`), which is what makes the closed forms and the Tier-2 parity meaningful.
4. **No error found in either extract.** The independent `ke_ref` reference reproduces the
   production stiffness to `≤ 1e-10` on all four geometries, and the rigid-body /
   uncorrected-shear / integration-rule assertions hold after the frame fix. No vision
   re-read was needed and neither extract was edited.
5. **The `M`/`K_T`/`f_int` scenario is mis-placed in WU4.** Requirement 2's scenario names
   the local and global mass, tangent and internal-force matrices, which the design assigns
   to WU5 (tasks 6.1/6.2); task 5.6 inherits the mismatch. WU4 asserts the part it owns and
   records the deferral.

### Remaining tasks (WU4-assembly)

All WU4 tasks 5.1–5.8 are complete. Task 2.4 (the boundary-traction loader and 48-DOF
dense patch assembler) remains `- [ ]` and was **not** in this unit's scope; its recorded
deferral text still says "it lands with WU4". The exact remaining unchecked lines in
section 5 are none. WU5 (tasks 6.1–6.4) and sections 7–13 remain pending.

### Workload / PR boundary (WU4-assembly)

Single PR, accepted `size:exception`. WU4 is one review slice: the stiffness assembly, its
independent reference and its nine tests, plus the frame-convention fix, extending
`mitc4_plusd.rs` with a green tree, its own verification and rollback = revert the WU4
block (the WU1–WU3 and WU4-constitutive bytes are otherwise untouched; note that the frame
fix does modify four WU2/WU3 operator lines, which is called out on task 5.3).
