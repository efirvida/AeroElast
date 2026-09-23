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
