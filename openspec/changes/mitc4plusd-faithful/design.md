# Design — the MITC4+/D element, faithfully implemented

Change: `mitc4plusd-faithful` · store: **openspec** (Engram mirror:
`sdd/mitc4plusd-faithful/design`) · phase: **design**

Scope of this document: the technical design that turns the proposal
(`openspec/changes/mitc4plusd-faithful/proposal.md`, rev. 2) and the 16-requirement
spec (`.../specs/mitc4plusd-element/spec.md`) into work that needs no further
interpretation. It fixes the six design decisions the proposal/spec left open, the
element's data and operator layout, the constitutive path, the Tier 1 and Tier 2
test design, the drill-operator verification strategy, the retirement sequencing,
the risks, and the work-unit slicing.

**It does not write tasks.** `tasks.md` is the next phase.

**Rev. 2 (corrective).** Rev. 1 claimed that `MeshTopology.node_normals` supplied the
per-node directors `V_n^i`. That symbol does not exist. Verified: `MeshTopology`
(`crates/aeroelast-core/src/assembly/topology.rs:46`) has exactly six fields —
`node_coords`, `connectivity`, `elem_types`, `n_nodes`, `n_elems`, `dofs_per_node` — and
`grep -rn "node_normals" crates/` returns nothing. The only per-node normals anywhere
are `Mitc4Precomputed::initial_normals` (`crates/aeroelast-core/src/elements/mitc4.rs:105`),
which the constructor fills with a single `e3` for all four nodes
(`crates/aeroelast-core/src/elements/mitc4.rs:582-584`, the comment reading *"All nodes
share the same normal in flat reference config"*) — a placeholder, not a director field.
ADR-4 is therefore re-derived below from the real code, and every symbol and file:line
this revision cites was re-resolved against the tree before it was written.

Non-goals restated (unchanged, not re-opened): no shell-element trait, no 5-DOF
stride migration, no PyO3 signature change, no `aeroelast-solvers` work, no fix of
the two pre-existing Python failures, no tolerance loosening, **no numerical factor
anywhere** (no penalty, no shear correction, no SRI).

## 0. Inputs read for this design

| Input | What it fixed |
| --- | --- |
| `openspec/changes/mitc4plusd-faithful/proposal.md` (rev. 2) | identity, preserved invariants, fixture policy, retirement scope |
| `openspec/changes/mitc4plusd-faithful/specs/mitc4plusd-element/spec.md` | the 16 requirements / 39 scenarios, the exact NEW test names, the evidence gaps G1–G7 |
| `docs/formulations/mitc4plus-2017-extract.md` | paper A Eqs. (1)–(27), the `a_E > 0` sign correction, the 2×2×2 / no-factor statements |
| `docs/formulations/mitc4plusd-2025-extract.md` | paper B Eq. (18) verbatim, Fig. 7(a) mesh (8 nodes / 5 elements), Fig. 7(b)(c)(d) BC sets, `V^D` single-normal statement |
| `crates/aeroelast-core/src/elements/mitc4.rs` | the hybrid to be superseded; the dead `b_md_mitc4_plus` / `drill_midside_shape_derivatives`; `Mitc4Precomputed`; the API surface the assembly layers call |
| `crates/aeroelast-core/src/materials/{mod,laminate,isotropic,composite}.rs` | `ShellConstitutive` (5 fields, `materials/mod.rs:20-31`), `to_shell_constitutive[_with_offset]`, `compute_shear_stiffness`, `cs = k·G·h` on the **single-ply** branch (`laminate.rs:137-138`) but `d²/flex` with no scalar on the **multi-ply** branch (`laminate.rs:178-187`) |
| `crates/aeroelast-core/src/assembly/{assembler,topology}.rs` | `PrecomputedElem::Quad`, `ElemType::{Mitc4,Mitc4Composite}`, `build_constitutive_mitc4`, `MaterialSpec` |
| `crates/aeroelast-py/src/{elements,assembler}.rs` | the PyO3 surface that must not move; where the composite `MaterialSpec` is built |
| `openspec/config.yaml` | the design rule to record paper tests separately from cross-reference tests; baselines; maturin gate |

Two facts discovered while reading that change the design and are recorded here as
**resolved findings**, not open questions:

**F1 — paper A's Eq. (27a–c) reduces to Eq. (18) exactly.** Expanding
`a_A e_rr(A) + a_B e_rr(B)` with Eq. (17) and `a_A = c_r(c_r−1)/(2d)`,
`a_B = c_r(c_r+1)/(2d)` gives `(c_r/d)[c_r(e_rr|con + e_rs|bil) − e_rr|lin]`, which is
paper A's Eq. (21) first term verbatim. So Eq. (27) *is* the closed form of
Eqs. (21)+(26), and at `c_r = c_s = 0` it collapses term-by-term onto Eq. (18):
`¼r(e_rr(A) − e_rr(B)) = ½ r e_rr|lin` and `¼s(e_ss(C) − e_ss(D)) = ½ s e_ss|lin`.
**The element therefore implements Eq. (27a–c), not Eq. (21) or Eq. (26)**, which
removes the transcription ambiguity in Eq. (21) from the critical path.

**F2 — the extract's Eq. (21) is missing its leading term; CONFIRMED, with proof.**
Paper A's Eq. (27c) is transcribed at
`docs/formulations/mitc4plus-2017-extract.md:388-390` as
`ẽ_rs^m = ¼(r + 4a_A·rs) e_rr^m(A) + ¼(−r + 4a_B·rs) e_rr^m(B) + ¼(s + 4a_C·rs) e_ss^m(C)
+ ¼(−s + 4a_D·rs) e_ss^m(D) + (1 + a_E·rs) e_rs^m(E)`. The coefficient of `e_rs^m(E)`
contains a **bare `1`** — exactly the leading `e_rs^m|bil` term that the printed Eq. (21)
omits. So Eq. (27c) is the correct closed form, and the printed Eq. (21)
(`docs/formulations/mitc4plus-2017-extract.md:347-349` in `B`-form, re-expressed at
`:359-361`) is **defective as printed**: with all `B_i` vanishing on a flat element
(`c_r = c_s = 0`), it would give `ẽ_rs^m|bil = 0`, contradicting the paper's own Eq. (22)
(`ẽ_rs^m|bil = e_rs^m|bil` when flat, extract `:349`). Two consequences:

1. `docs/formulations/mitc4plus-2017-extract.md` **faithfully transcribes the printed
   (defective) Eq. (21)** — the defect is the paper's, not the extract's. WU0 records a
   note against the printed equation rather than "fixing" the quote, so the artifact
   keeps reporting what the paper prints.
2. The element **must use the Eq. (27a–c) form**, which carries the leading term. It
   never needs the printed Eq. (21) form; only the test
   `test_t1a_membrane_eq22_flat_tying_condition` evaluates it, and F2 is why that test
   is expected to pass rather than to expose a code defect.

---

## 1. The six architecture decisions

Each decision states the choice, the rationale, the rejected alternatives, and what
would falsify it.

### ADR-1 — Uncorrected transverse shear: mechanism (b) with a material-side factor channel

**Decision.** The element owns the removal of the shear correction factor, and the
factor it must remove is an explicit input:

1. `ShellConstitutive` **keeps exactly its five fields** and gains one additive
   accessor (`crates/aeroelast-core/src/materials/mod.rs`):

   ```rust
   impl ShellConstitutive {
       /// Paper A p. 410: the element carries no numerical factor, so it needs the
       /// **uncorrected** transverse-shear stiffness.
       ///
       /// `cs` is documented in `k·G·h` form. `applied_k` is the scalar correction
       /// the material model applied when it built `cs` — 1.0 when it applied none
       /// (e.g. a multi-ply laminate's energy-equivalent section stiffness).
       /// Additive: `cs` and every existing consumer are untouched.
       pub fn transverse_shear_uncorrected(&self, applied_k: f64) -> Matrix2<f64> {
           if applied_k > 0.0 { self.cs / applied_k } else { self.cs }
       }
   }
   ```

2. `Laminate` gains one additive accessor (`materials/laminate.rs`) reporting the
   factor `compute_shear_stiffness` actually applied:

   ```rust
   /// The scalar correction `compute_shear_stiffness` applied to `cs`:
   /// `shear_correction_factor` on the single-ply branch, `1.0` on the
   /// multi-ply energy-equivalence branch (which applies no scalar factor).
   pub fn applied_shear_correction_factor(&self) -> f64 {
       if self.plies.len() == 1 { self.shear_correction_factor } else { 1.0 }
   }
   ```

3. The element's constructor takes the factor and stores the uncorrected stiffness:

   ```rust
   pub struct Mitc4PlusDPrecomputed {
       pub constitutive: ShellConstitutive,        // cm, cb_coupling, cb, cm_raw
       pub applied_shear_correction: f64,
       pub cs_uncorrected: Matrix2<f64>,           // what the shear block uses
       …
   }
   ```

4. Every element path supplies the factor where the repository knows it:
   `MaterialSpec::Isotropic { shear_correction }`; `MaterialSpec::Composite` gains one
   **internal** field `applied_shear_correction` (filled in
   `crates/aeroelast-py/src/assembler.rs` from
   `corrected_lam.applied_shear_correction_factor()`); the PyO3 isotropic batch
   kernels use the `shear_correction` argument they already take.

**Rationale.** The observable the spec fixes is "the element's transverse-shear
contribution is invariant to `k` when `cs = k·G·h`". The only way to satisfy that
literally is for the element to be *told* `k`, because a bare `ShellConstitutive`
cannot distinguish "`cs` carries a scalar factor" from "`cs` is already uncorrected".
The material-side accessor is what makes the removal *sound* rather than guessed: for
a multi-ply laminate `cs` is the equilibrium energy-equivalent section stiffness with
no scalar in it, so `applied_shear_correction_factor() == 1.0` and the stiffness
passes through verbatim.

**Evidence (re-resolved).** The repository's `ShellConstitutive` doc comment reads
`Transverse shear stiffness (2×2): force-shear strain [k·G·h] form`
(`crates/aeroelast-core/src/materials/mod.rs:28-29`) and the struct carries no `k`
(fields `cm`, `cb_coupling`, `cb`, `cs`, `cm_raw` at `materials/mod.rs:20-31`).
`Laminate` does carry it (`crates/aeroelast-core/src/materials/laminate.rs:39`) and
`compute_shear_stiffness` applies it in exactly two ways:

- single ply, scalar correction: `let k = self.shear_correction_factor;`
  (`laminate.rs:137`, applied at `:138`) — this is the `k·G·h` form the doc comment
  describes;
- multi-ply, energy-equivalent section stiffness: `cs_55 = d11*d11/flex_55` and
  `cs_44 = d22*d22/flex_44` (`laminate.rs:178-187`), with **no `k` in the numerator** —
  the `let k = self.shear_correction_factor;` at `laminate.rs:176` is used only in the
  degenerate fallbacks (`:181`, `:187`, `:196`).

So the `[k·G·h]` doc comment is **misleading for the multi-ply branch**, and a plain
`cs / k` at the element boundary is unsound for it — which is exactly why the
mechanism must be material-reported rather than assumed.

**Alternatives rejected.**

- **(a) Extend `ShellConstitutive` with the correction factor** — rejected twice over.
  The spec's preserved-invariant scenario requires `ShellConstitutive` to expose
  *exactly* the five fields `cm`, `cb_coupling`, `cb`, `cs`, `cm_raw`
  (`crates/aeroelast-core/src/materials/mod.rs:20-31`; spec scenario at
  `.../specs/mitc4plusd-element/spec.md:434-441`), so a sixth field is a spec violation;
  and a single stored `k` cannot represent the multi-ply case, where the correct value
  of the stored factor is `1.0` and the field would then be dead weight on every
  laminate.
- **(c) Divide by `k` at the element boundary, full stop** — rejected as unsound.
  `cs / k` is correct only when the material model applied exactly that scalar.
  `Laminate::compute_shear_stiffness` does not: for a multi-ply laminate `cs_55 =
  d11²/flex_55` and `cs_44 = d22²/flex_44`
  (`crates/aeroelast-core/src/materials/laminate.rs:178-187`) carry no `k` at all, and
  `k` survives only in the degenerate fallbacks (`laminate.rs:181`, `:187`, `:196`).
  Dividing a multi-ply laminate's `cs` by 0.75
  would inflate the section's transverse shear stiffness by 33 % and move
  `test_composite_beam_parity.py` / `test_orthotropic_shell_parity.py` — i.e. it would
  break the preserved-laminate invariant. Mechanism (c) becomes correct only once the
  material model reports the factor it applied, which is mechanism (b).
- **Reading `G·h` out of `cm[(2,2)]`** (tempting: for an isotropic material
  `cm[(2,2)] = G·h` and `cm_raw[(2,2)] = G`) — rejected because it is false for every
  laminate: `cm[(2,2)] = A66 = ∫Q̄66 dz` is the *in-plane* shear stiffness, unrelated to
  the transverse shear stiffness of the section.

**How the spec's discriminating test detects a residual `k`.** Two `pre` values are
built from the *same* single-ply laminate family with `k ≠ k'`:
`pre(k)` = `(constitutive = lam(k).to_shell_constitutive() → cs = k·G·h,
applied_shear_correction = k)` and `pre(k')` likewise. Their transverse-shear blocks
must agree to `1e-10` relative, and each must equal the closed-form
`∫ B_γᵀ (G·h) B_γ dA` (2×2 rule) to `1e-10`. The **non-vacuity control** builds a third
`pre` = `(cs = k·G·h, applied_shear_correction = 1.0)` and asserts its shear block
differs from `pre(k)`'s by more than `1e-3` relative — that is the separation that a
residual `k` would produce. A fourth control uses the isotropic constitutive and the
`5/6` factor with the same shape of assertion.

**Preserved-laminate invariant.** `cs`, `compute_shear_stiffness`,
`to_shell_constitutive`, `to_shell_constitutive_with_offset`, `shear_correction_factor`
and every existing number are untouched; both new accessors are additive and change no
result. `test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged`
pins the public surface.

**Falsified by.** A laminate whose element shear block moves when
`shear_correction_factor` moves while `G·h` does not — that is the failure the
discriminator reports.

**AMENDMENT (this session) — the multi-ply `cs` DOES carry a numerical factor, and the
verbatim step above is falsified by measurement.**

The decision above reasons: the multi-ply `cs` "is the equilibrium energy-equivalent
section stiffness with no scalar in it", therefore it "passes through verbatim". The
first clause is true of the *formula* and false of the *number*. `compute_shear_stiffness`'s
energy-equivalence branch builds `cs` from the piecewise-quadratic shear stress profile,
and for a homogeneous stack that derivation **yields exactly `(5/6)·G·h`** — the `5/6`
of the parabolic distribution, a numerical factor the element must not see. Measured for
four identical isotropic plies of `h/4` (`E1=E2=E3`, `G1=G2=G3`, `ν=0.3`):
`cs_55 = cs_44 = (5/6)·G·h`, `cs_45 = 0`, against the plain section integral
`a_ij = Σ c̄s_ij t_k = G·h` by construction. Pinned in Rust by
`Laminate::shear_stiffness_uncorrected`'s test, which also asserts the two are separated
by more than `1e-3` relative so it cannot pass vacuously.

**Consequence, measured.** Because the isotropic path removes its `k` and lands on `G·h`
while the composite path passed the corrected `cs` through, the two models of the SAME
homogeneous plate disagreed in the shear block by 20%, and the repository's own
equivalence tests caught it: `test_material_suite.py::TestIsoEquivalence::test_n_iso_plies_equal_single_layer_mitc4`
at 3.434% (tolerance 1%), and `test_orthotropic_shell_parity.py::test_multi_layer_iso_equivalence`
at 5.42% (tolerance 1e-4). The second is slow-marked and was therefore invisible to the
`-m "not slow"` gate.

**Amended mechanism.** A scalar cannot express the correction, because the
energy-equivalent `cs` is not a scalar multiple of the plain integral in general
(`cs_55/a55 ≠ cs_44/a44`); so the material supplies the **uncorrected matrix**, not a
factor:

1. `Laminate::shear_stiffness_uncorrected() -> Matrix2<f64>` — the plain integral,
   equal to `G·h` for a homogeneous stack. Additive: no stored field, `cs` and every
   existing `Laminate` number unchanged.
2. `MaterialSpec::Composite` gains the **internal** field `cs_uncorrected: [f64; 4]`
   beside `applied_shear_correction`, filled in `crates/aeroelast-py/src/assembler.rs`
   from `corrected_lam.shear_stiffness_uncorrected()`. The raw-dict path accepts an
   optional `cs_uncorrected` key and otherwise keeps `cs` verbatim (its legacy
   contract; identical to the constructor's own derivation for a factor of 1.0).
3. The assembler overwrites `pre.cs_uncorrected` from
   `mitc4_uncorrected_shear(mat, &pre.constitutive)` — the same pass-2 pattern it
   already uses for the mesh-consistent nodal directors. `ShellConstitutive` keeps
   exactly its five fields, `applied_shear_correction` stays as the record of what the
   material model applied, and the element's own constructor derivation is unchanged for
   direct callers.

**External confirmation that the uncorrected modulus is the right input.** The element
reproduces the published twisted-beam cells of Ko et al. (2017), C&S 193:187-206,
Tables 12/13 — a shear-locking benchmark — to 0.01-0.25%, and the pinched-cylinder cells
of Tables 8/9 to 0.04-1.41%, while its shear input is the uncorrected modulus. The
paper's element therefore uses the uncorrected modulus too, which is what p. 410's
"does not include any numerical factor" requires. Measured after the amendment:
`test_ko2017_performance.py` stays at 31 passed (the isotropic path is untouched),
`test_material_suite.py` 41 passed, `test_orthotropic_shell_parity.py` 3 passed, and the
composite files 33 passed.

**Requirement 16 is unaffected.** Its observable — the element's transverse-shear
contribution is invariant to `k` when `cs = k·G·h` — still holds, and the isotropic path
that exercises it is bit-identical. The amendment removes a *false premise* from ADR-1's
rationale, not a commitment of the spec.

**Still open from this amendment.** Five test-side composite dict builders declare `cs`
explicitly (`test_composite_b_coupling.py` ×3, `test_composite_beam_parity.py`,
`test_rust_composite.py`) and therefore still hand the element the corrected stiffness;
they pass today under the legacy convention and moving them is a separate measured step.

### ADR-2 — Internal structure: a new module, `Mitc4PlusDPrecomputed`, unchanged 24-DOF layout

**Decision.** The faithful element lands in a **new module**
`crates/aeroelast-core/src/elements/mitc4_plusd.rs` with its own precomputed type
`Mitc4PlusDPrecomputed`. `crates/aeroelast-core/src/elements/mod.rs` gains one line
(`pub mod mitc4_plusd;`). `mitc4.rs` is **not touched until the retirement work unit**.

- **Local 24-vector layout is unchanged**: node `i` owns slots `6i .. 6i+5`, holding
  the **local-frame** components `(u, v, w, θx, θy, θz)`. `t3`, `build_t24` and
  `transform_to_global` keep their present meaning, so `_FAMILY_PROPERTIES[SHELL] =
  (6, 3)`, the global stride, the `[f64; 576]` / `[f64; 24]` PyO3 shapes and the VTK
  `θ = (θx, θy, θz)` output form are untouched.
- **Rotation convention.** The paper's `(α_i, β_i)` are never stored. The element uses
  the paper's own equivalent form
  `u_b = ½ Σ a_i h_i (θ_i × V_n^i)`, which is *identically* `−V_2^i α_i + V_1^i β_i`
  with `α_i = θ_i·V_1^i`, `β_i = θ_i·V_2^i`, because for a right-handed
  `(V_1^i, V_2^i, V_n^i)`: `V_1 × V_n = −V_2` and `V_2 × V_n = V_1`. The per-node basis
  is derived from the stored director: `V_1^i = normalize(e1 − (e1·V_n^i)V_n^i)`,
  `V_2^i = V_n^i × V_1^i`.
- **`θ_z` is carried, not consumed by the 2017 core.** `θ_i × V_n^i` is blind to the
  component of `θ_i` along `V_n^i`, so slot `6i+5` contributes nothing to `u_b`, `e^m`,
  `e^b1`, `e^b2` or the MITC4 shear — exactly the paper's statement that the 2017
  element has no drilling DOF. It enters only through the 2025 drill term, as
  `θ_i^D = θ_i · V^D` (ADR-3).
- **Function names in the new module** (all new; the hybrid's names in `mitc4.rs` are
  left alone and deleted at retirement):

  | Kind | Names |
  | --- | --- |
  | geometry / coefficients | `compute_characteristic_vectors`, `compute_node_directors` (ADR-4), `compute_membrane_coefficients_2017`, `compute_j3d_enriched`, `covariant_to_local_mapping`, `shear_covariant_to_local`, `regularized_inverse_2x2`, `shape_functions`, `shape_function_derivatives` |
  | B-operators | `b_membrane_2017` (Eq. 27), `b_bending_2017` (Eqs. 7c/7d → `e^b1`, `e^b2`), `b_shear_mitc4` (DB84 Eq. 3), `b_drill_membrane_2025` (Eq. 18), `drill_midside_shape_derivatives`, `drill_jacobian_ratio`, `resultant_moment_matrix` |
  | assembly | `compute_ke_local`, `compute_ke_global`, `compute_fint_global`, `compute_kt_global`, `compute_me_global`, `compute_me_composite_global`, `compute_body_load_global`, `compute_k_sigma_global`, `compute_centrifugal_prestress`, `compute_element_stress`, `build_t24`, `transform_to_global`, `extract_elem_disp_24` |
  | corotational (retargeted, formulation-independent) | `quaternion_to_matrix`, `quaternion_from_vector`, `rotate_vector_by_quaternion`, `quaternion_multiply`, `update_normals_with_displacements`, `polar_decomposition`, `log_strain_from_polar`, `compute_membrane_strain_log`, `update_corotational_frame`, `frame_incremental_rotation` |

- **Names deliberately NOT created** (the deviation set the retirement removes, which
  must never be re-introduced): `k_drill`, `drilling_scale`, `beta_w`/`BETA_W`,
  `erc_*`/`compute_ke_local_erc`/`b_erc`/`erc_covariant_to_local`, `b_drill`,
  `b_drill_grad`, `drill_dh`, `bubble_function`/`bubble_derivatives`/`GpBubble`/
  `b_kappa_bubble`/`b_gamma_mitc4_plus`, `b_m_standard`, `green_lagrange_strain`,
  `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`, `m8_erc`,
  `pseudo_inverse_8x8`, `compute_enhanced_drill_stiffness`,
  `compute_drill_warping_moment`, the SRI split, `hg_*`.

**Rationale.** A new module keeps the hybrid **byte-identical and green** while the new
element is built and tested, which is what "nothing is deleted before the gate" and
"the tree is green throughout" require; it makes rollback a one-line dispatch revert;
and it lets a reviewer read one file top-to-bottom against the two papers. It also
keeps the constructor honest: the new element needs neither `e_mod` nor
`drilling_scale`, and dropping them is the visible proof that the penalty is gone.

**Alternatives rejected.**

- **Rewrite `mitc4.rs` in place with a `const Formulation` switch.** Rejected: the
  hybrid stays compiled and its `k_drill`/`erc`/`beta_w` symbols remain in the file the
  element is documented to run from, so the identity checks of Requirement 1 (and its
  supporting grep) cannot be satisfied until the retirement, and a reviewer cannot see
  the new formulation as a coherent whole. The switch would also be a scaffold the
  retirement has to remove, i.e. extra churn inside the riskiest unit.
- **New `ElemType::Mitc4PlusD` variant, element opt-in.** Rejected: the spec's Tier-2
  requirement is the *gate* that proves the retirement did not cost a capability. If
  the new element is opt-in and the old one still dispatches, Tier 2 never exercises
  the new element, and the gate becomes decorative. The change's identity statement
  ("the element the repository runs") also requires the dispatch to move.
- **Port the hybrid's `PrecomputedElem::Quad` payload in place (same type, new
  semantics).** Rejected: it would silently change every Tier-2 expectation in the
  same commit that introduces the formulation, with no reference left to diff against.

**Falsified by.** Any Tier-2 test that fails only because of the retargeted type and
not because of a formulation difference — that would mean the new module is not a
faithful port of the *contract*, only of the equations.

### ADR-3 — `V^D`, `j`, `j0`, and the drill operator's `c_r`/`c_s`

**Decision.**

- **`V^D = n_vec`.** Paper B: `V^D` is "the single normal to the plane P at the element
  centre". With the paper's geometry `x = x_m + t x_b`, the base vectors at `t = 0` are
  `g_r(r,s,0) = ∂x_m/∂r = x_r + s x_d` and `g_s(r,s,0) = x_s + r x_d`, so at the centre
  `g_r = x_r`, `g_s = x_s` and `V^D = normalize(x_r × x_s)` — which **is** paper A
  Eq. (10)'s `n`. So the code's existing computation (normalize of the flat centre
  cross product) is already Eq. (10) and is kept, stored once as `pre.v_d`. One vector
  per element, never per node.
- **`j = det[g_r g_s g_t]|_{(r,s,0)}` as a scalar triple product.** `g_t = ∂x/∂t = x_b`
  (the director field is the `t`-base vector), so
  `j(r,s) = (x_r + s x_d) · ((x_s + r x_d) × x_b(r,s))` and `j0 = j(0,0)`. This is
  implemented **exactly** as the triple product, not as a ratio of surface measures.
  For a flat element with a common director it reduces to `(h/2)|x_r × x_s|`, so
  `j0/j = |x_r × x_s|(0,0) / |x_r × x_s|(r,s)` — the form the dead code used. The dead
  code's form is therefore a *documented special case*, not a different quantity; the
  exact form matters only on warped geometry.
- **The drill operator carries `θ^D`, not `θ_z`.** `θ_i^D = θ_i · V^D` with `θ_i` the
  local-frame rotation triple, so the operator has non-zero entries in the `θx` and
  `θy` slots whenever `V^D` is not parallel to the local `e3`. On a flat element
  (`V^D = e3`) it touches slot `6i+5` only.
- **The naming collision is made visible.** Paper A's `c_r = x_d·m^r`, `c_s = x_d·m^s`
  are stored as `pre.c_r_mem`, `pre.c_s_mem` (with `pre.d_mem` and
  `pre.a_coeffs = [a_A..a_E]`); paper B's `c_r^I = x_m^I·(−x_r^I × V^D)`,
  `c_s^I = x_m^I·(x_s^I × V^D)` are locals named `cr_md` / `cs_md` inside
  `b_drill_membrane_2025`, with a doc comment naming both papers, both symbols and both
  definitions. A test asserts they are different quantities on a warped element.
- **The drill's covariant→local transform is the element-centre one** (paper B
  Eq. (21): "using the constant element-centre base vectors"), i.e.
  `covariant_to_local_mapping(j_loc(0,0))`, not the point-wise mapping the membrane
  uses. Pinned by a test on a warped quad.

**Rationale.** Each of these is a place where the dead transcription could be wrong in
a way that changes nothing until the term is live. Making them explicit, named, and
separately tested is what converts "trusting a transcription" into evidence.

**Alternatives rejected.** (i) Per-node `V^D` — explicitly refuted by paper B ("using
the single normal to the plane P results in an overall improved element behavior");
(ii) `j0/j` as the surface-measure ratio — exact only for flat geometry, and the
warped case is precisely what this change exists to fix; (iii) storing `θ^D` as the
6th DOF instead of `θ_z` — would change the PyO3/Python DOF meaning and the VTK output
form, both frozen by the spec.

### ADR-4 — Per-node directors `V_n^i` and thicknesses `a_i`

**What the code actually provides (re-verified, Rev. 2).** There is no director field
to consume:

- `MeshTopology` (`crates/aeroelast-core/src/assembly/topology.rs:46-59`) has exactly six
  fields — `node_coords`, `connectivity`, `elem_types`, `n_nodes`, `n_elems`,
  `dofs_per_node` — and is built only through `MeshTopology::new`
  (`topology.rs:72-92`). There is no `node_normals`; `grep -rn "node_normals" crates/`
  returns nothing.
- `Mitc4Precomputed::initial_normals` (`crates/aeroelast-core/src/elements/mitc4.rs:105`)
  is the only per-node normal array in the tree, and the constructor sets every entry to
  the *single* element normal `e3` (`mitc4.rs:582-584`, comment: *"All nodes share the
  same normal in flat reference config"*). It exists for the quaternion update
  (`mitc4.rs:3668`), not for `V_n^i`.
- The element constructor receives only `node_coords: &[f64;12]` and the constitutive
  (`mitc4.rs:497-505`), i.e. **no connectivity, no neighbours, no nodal normals**.

So the paper's per-node `V_n^i` (`docs/formulations/mitc4plus-2017-extract.md:32-37`,
Eq. (1)/(8a) at `:179`) must be *produced*, not read.

**Options evaluated.**

- **(A) Area-weighted average of the adjacent elements' Eq. (10) normals, supplied by
  the mesh.** Correct shell-theory nodal director. Requires (i) a new field on
  `MeshTopology`, (ii) recomputation after `update_reference` moves the nodes
  (`crates/aeroelast-core/src/assembly/assembler.rs:217-230` mutates
  `self.topology.node_coords` in place and then rebuilds each element), and (iii) a new
  argument on the element constructor. Blast radius if chosen: `topology.rs` (one field
  + the loop that fills it inside `new`), `assembler.rs` (both `Mitc4Precomputed::new`
  sites, `assembler.rs:161-167` and `:252-254`, plus the rebuild loop), and — because
  `MeshTopology::new`'s three-argument signature would be unchanged — **no** change to
  `crates/aeroelast-py/src/assembler.rs` (`MeshTopology::new` at `:92`, `:146`, `:405`)
  and no change to the frozen PyO3 surface. It stays inside scope, but it buys the
  neighbour information the element does not otherwise receive at the cost of a cached
  field that must be kept consistent with a mutated `node_coords`, and it still does not
  serve the single-element path (no neighbours exist).
- **(B) Compute the per-node directors locally inside the element from its own four 3D
  coordinates.** No new field, no constructor argument, no assembler change, and it is
  automatically correct after `update_reference` because the element is rebuilt from the
  new coordinates. It cannot see the *mesh* normal, so on a curved mesh adjacent
  elements get slightly different nodal directors (a discontinuity in the director
  field).
- **(C) Extend `MeshTopology` with the field and also keep a local fallback.** Both
  mechanisms, both tested. Rejected as overbuilding: two code paths for one quantity,
  and the mesh path has no test in this change that the local path does not also
  satisfy.

**Decision: option (B).** The element computes its own per-node directors in
`Mitc4PlusDPrecomputed::new`, from the twelve coordinates it already receives, with no
new input and no change to `MeshTopology`, the assembler's data flow, or PyO3.

- **Definition (`compute_node_directors`, new, in `mitc4_plusd.rs`).** Split the quad
  into four sub-quads by its centre `x_c = ¼ Σ x_i`; sub-quad `k` is
  `(x_c, x_k, x_{k+1}, x_{k+2})` in the element's counter-clockwise node order. Node `i`
  is shared by sub-quads `k = i−1` and `k = i`, so
  `V_n^i = normalize( A_{i−1} n_{i−1} + A_i n_i )`, where `n_k` and `A_k` are the outward
  normal and area of sub-quad `k`, each sign-aligned so that `n_k · n_vec > 0` with
  `n_vec` the Eq. (10) element normal. This is an area-weighted average of the element's
  own corner normals — the local analogue of option (A) — and it needs no neighbour.
  The codebase already contains this averaging pattern for the element normal itself
  (`crates/aeroelast-core/src/elements/mitc4.rs:253-272` averages two sub-quad normals to
  build `e3`); the new function applies it per node instead of per element.
- **Single-element case (Tier 1, no neighbours).** There is no neighbour path, so the
  element uses exactly the above; on a flat element all four sub-quad normals equal
  `n_vec`, hence `V_n^i = n_vec` for all `i` to round-off. Every Tier-1 fixture is flat
  (the star patch, the flat square, the flat distorted quad) or deliberately warped
  (F-W), and the PyO3 per-element kernels call the same constructor — so the Tier-1 and
  batch results are the `V_n^i = n_vec` case exactly, with no separate code path.
- **Limitation, recorded.** On a curved/warped mesh the local directors are
  element-local rather than mesh-shared, so the director field is discontinuous across
  element boundaries. This is the price of option (B) and it is named, not hidden; if a
  Tier-2 result shows it matters, option (A) is the follow-up and its blast radius is
  the one stated above. It is not a deviation from paper A — the paper requires `V_n^i`
  to be per node, and it is.
- **`a_i`.** The struct stores `a_i: [f64;4]`, filled with the scalar `thickness` for
  all four nodes. The kinematics is written in the paper's per-node form, so a tapered
  thickness is a one-line change, but **no input carries a per-node thickness in this
  change** and no taper case is tested. Recorded as a limitation, not a deviation.
- **`initial_normals` is not reused.** The new module's `vn` is computed as above; the
  hybrid's `initial_normals` placeholder (`mitc4.rs:105`, `:582-584`) stays with the
  hybrid and is deleted at S4. No new code reads `initial_normals`.

**Rationale.** The paper's `a_i`/`V_n^i` are per node by construction, and a single `e3`
for all four nodes makes `x_b = ½ Σ a_i h_i V_n^i` constant, silently deleting Eq. (8a) —
which the 2017 extract names as the **first of the three unchecked parts of the
formulation** (`docs/formulations/mitc4plus-2017-extract.md:491-499`, items 1–3 being
Eq. (8a), Eq. (8b) and the `∂x_b` terms of Eqs. (7c)/(7d)). Option (B) restores Eq. (8a)
with the smallest possible blast radius and with no cached state to keep in sync, which
is the change's stated main risk-reducer (proposal §5/§6.2: the assembly layers are
untouched).

**Alternatives rejected.** (i) Option (A) alone — correct in principle but adds a cached
field that `update_reference` must keep consistent, and still cannot serve the
single-element path; kept as the named follow-up if a Tier-2 result demands mesh-shared
directors. (ii) Option (C) — two mechanisms for one quantity, one of them untested.
(iii) The element scanning neighbours itself — the constructor receives no connectivity
(`mitc4.rs:497-505`), so this is not available without a signature change that reaches
into the assembly layers. (iv) Angle-weighted normals — a refinement with no evidence
behind it in either paper. (v) Per-node thickness from the mesh — no such input exists
and inventing one is scope creep.

**AMENDMENT (this session) — the trigger named above FIRED, and option (A)'s CONTENT is
already in production. What is stale is the record, and what is left is the contract.**

This ADR kept option (A) as "the named follow-up if a Tier-2 result demands mesh-shared
directors". That is what happened. WU9e measured a Tier-2 benchmark — the
MacNeal-Harder twisted beam at `t/L = 0.02667`, the thick case — and with option (B)
alone the element gives `1.06965` (N=4) / `1.25132` (N=8) against the published MITC4+
cells `0.9960` / `0.9968`; with a mesh-consistent nodal director it gives `0.99515` /
`0.99733`. So the "Limitation, recorded" paragraph above is superseded: (B) alone is
7-25% off on a Tier-2 case and the mesh-shared director is what fixes it.

Production therefore does **(A) by content, (B) by plumbing**:
`crates/aeroelast-core/src/assembly/assembler.rs:135-184` (`MeshAssembler::new`) and
`:251-283` (`update_reference`) compute `mitc4_nodal_directors(&topology, &materials)` —
the area-weighted mean of the adjacent elements' element-local `vn`, normalized per
global node — and then overwrite `pre.vn[a] = nodal_director[node]` after the
constructor has already computed the element-local directors.

Measured today, on the cells that fired the trigger:
`test_3_5_twisted_beam_tables_12_to_13` passes all four cases at **0.01-0.25%** of the
published Tables 12/13 (thick in-plane `0.9981` vs `0.9971`, thick out-of-plane `0.9998`
vs `0.9973`, thin `0.9972` vs `0.9978` and `0.9981` vs `0.9982`), and
`test_ko2017_performance.py` is 31 passed.

**What is stale (a record defect, not a behaviour one).** (i) The decision paragraph still
reads "Decision: option (B)". (ii) The task list's read-only inputs table still justifies
`topology.rs` as "**Untouched by design** (ADR-4 option B: the element computes its own
per-node directors, so `MeshTopology` gains no field)" — `topology.rs` really is
untouched, but the stated *reason* is no longer the whole story, because the assembler
overwrites what the element computed. (iii) Option (C) — "Both mechanisms, both tested…
rejected as overbuilding: two code paths for one quantity" — is de facto what runs: the
local path is the fallback that gets overwritten.

**The debt this leaves, named.** `pre.vn` is a value the constructor computes and the
assembler discards — the same shape as the additive `compute_ke_local(pre)` reference that
was removed from the geometric tangent. The clean forms, in increasing blast radius, are:
(a) document the overwrite as the pass-2 channel it is (what this amendment does);
(b) move the mesh-consistent directors into the constructor's contract so nothing is
computed and thrown away; (c) option (A) proper, with the `MeshTopology` field and the
`update_reference` recomputation this ADR already scoped. (b) is the one to take next: it
stays inside the element's own file plus the two `Mitc4PlusDPrecomputed::new` call sites.

**Acceptance for whichever form is taken.** The four `test_3_5_twisted_beam_tables_12_to_13`
cases stay at 0.01-0.25% of Tables 12/13; the flat Tier-1 fixtures keep `V_n^i = n_vec` to
round-off; `test_ko2017_performance.py` stays at 31 passed.

### ADR-5 — Retirement sequencing (§8)

**Decision.** Four ordered stages, each leaving the tree green: **(S1)** build and test
the new element as a dead-but-tested module; **(S2)** move the layout-bound Tier-2 Rust
tests onto the new element and delete the hybrid's copies (the names must exist exactly
once); **(S3) the flip** — route `ElemType::{Mitc4,Mitc4Composite}`, the two
`Mitc4Precomputed::new` call sites, `update_reference`, `extract_elem_disp_24` and the
PyO3 MITC4 entry points to `mitc4_plusd`, and supply the ADR-1 factor channel; **(S4)
the retirement** — delete the hybrid from `mitc4.rs` and its deviation tests.

The gate is checked *at S3* (Tier 1 + Tier 2 with the new element live) and *again at
S4* (Tier 1 + Tier 2 + the static checks). S4 is a pure deletion and is revertible on
its own. **The PyO3 surface does not move at any stage**: `crates/aeroelast-py/src/elements.rs`
keeps its Python-facing signatures, its `[f64;576]`/`[f64;24]` shapes and every
`#[pyfunction]`/`#[pyclass]` name; the now-unused `e_mod` argument of the MITC4 kernels
is kept and ignored, with a comment naming the out-of-scope signature change.

**Rationale.** The dangerous step is S3, because that is the moment the Tier-2 suite
starts measuring the new element; keeping it separate from both the test move (S2) and
the deletion (S4) means a Tier-2 regression at S3 is attributable to the dispatch and
nothing else, and a revert is one line.

**Alternatives rejected.** (i) Flip and retire in one unit — an attribution dead end;
(ii) keep both elements dispatchable and retire the hybrid's code without moving the
dispatch — Tier 2 would never gate the change; (iii) delete the hybrid's tests at S4
only, leaving the same test name in two modules — ambiguous for
`cargo test -p aeroelast-core <name>` and for the traceability test.

### ADR-6 — The mid-surface restriction (spec gap G7)

**Decision.** The element path uses the mid-surface form only:
`Laminate::to_shell_constitutive()` / the equivalent raw-ABD construction already used
by `MaterialSpec::Composite` (which passes `a`, `b`, `d`, `cs` with no offset). The
element path **does not call** `to_shell_constitutive_with_offset`, and
`test_identity_element_uses_midsurface_constitutive` asserts that the element's
constitutive equals `to_shell_constitutive()` on a laminate, with a supporting static
check that the element path references no offset variant. Whether the offset variant
may *ever* be used is deferred.

**Rationale for the deferral being substantive, not lazy.** The paper's kinematics is
`x = x_m + t x_b`, `t ∈ [−1,1]`, about the surface `x_m` interpolated from the four
nodes — i.e. the reference surface *is* the mid-surface, and `a_i` is the thickness
about it. `to_shell_constitutive_with_offset(z0)` shifts the ABD reference plane
(`B ← B − z0·A`, `D ← D − 2z0·B + z0²·A`) without any corresponding shift in the
kinematics. Consuming it therefore requires an extra offset term in `u_b`/`x_b` — a
derivation neither paper provides — and using it without that term would silently
change the membrane/bending coupling of every offset laminate while leaving the
kinematics wrong. So the offset path is *incompatible* with the paper's kinematics as
written, and the honest resolution is to forbid it and record why, not to test it.

**Alternatives rejected.** (i) Allow the offset variant and "verify" it with a test —
there is no equation to verify it against; (ii) delete `to_shell_constitutive_with_offset`
— it is public API and part of the preserved surface, so deleting it is a violation.

---

## 2. Element data and operator layout

### 2.1 Precomputed data (`Mitc4PlusDPrecomputed`)

```rust
pub struct Mitc4PlusDPrecomputed {
    // frame and transform
    pub local_coords: [[f64; 2]; 4],
    pub t3: Matrix3<f64>, pub e1: Vector3<f64>, pub e2: Vector3<f64>, pub e3: Vector3<f64>,
    pub initial_coords_3d: [[f64; 3]; 4],

    // constitutive (ADR-1)
    pub constitutive: ShellConstitutive,       // cm, cb_coupling, cb, cm_raw
    pub applied_shear_correction: f64,
    pub cs_uncorrected: Matrix2<f64>,

    // A Eqs. (1), (8a): per-node thickness and director
    // ADR-4: vn is computed locally by `compute_node_directors` from the element's
    // own four coordinates; the constructor takes no director argument.
    pub thickness: f64, pub a_i: [f64; 4],
    pub vn: [Vector3<f64>; 4],                 // V_n^i, local frame
    pub v1: [Vector3<f64>; 4], pub v2: [Vector3<f64>; 4],   // V_1^i, V_2^i

    // A Eqs. (9)-(11): characteristic vectors and dual basis
    pub x_r: Vector3<f64>, pub x_s: Vector3<f64>, pub x_d: Vector3<f64>,
    pub n_vec: Vector3<f64>, pub m_r: Vector3<f64>, pub m_s: Vector3<f64>,

    // A Eqs. (23)-(25), (27a-c)
    pub c_r_mem: f64, pub c_s_mem: f64, pub d_mem: f64, pub a_coeffs: [f64; 5],

    // B Eqs. (5), (18): drill normal and drill Jacobian ratio
    pub v_d: Vector3<f64>, pub j0: f64,

    // A Fig. 4 tying rows (displacement-based covariant membrane, Eqs. (15)-(17))
    pub b_rr_a: Vec24, pub b_rr_b: Vec24, pub b_ss_c: Vec24, pub b_ss_d: Vec24, pub b_rs_e: Vec24,

    // DB84 Eq. (3) tying operators: [A(top,s=+1), B(bottom,s=-1), C(right,r=+1), D(left,r=-1)]
    pub b_shear_tie: [SMatrix<f64, 2, 24>; 4],

    // 2x2 Gauss data on the surfaces
    pub gp: [Mitc4PlusDGp; 4],                 // j_loc, j_inv, j, sqrt_g, dh, g_r, g_s, g_t, g^r, g^s, g^t

    pub element_area: f64,

    // formulation-independent machinery, retargeted from the hybrid
    pub quaternions: [Vector4<f64>; 4],
    pub gp_initial_tangents: [GpTangents; 4], pub gp_initial_frames: [GpLocalFrame; 4],
}
```

Every stored quantity names its paper and equation in a doc comment, which is what the
traceability test parses (Requirement "Traceability of every live ingredient").

### 2.2 The B-operators

All operators are 24 columns wide and act on the local 24-vector
`(u, v, w, θx, θy, θz)` per node. Covariant strains are formed from the local-frame
components of the base vectors and the local-frame DOF, so the covariant components are
frame-independent scalars; each is then mapped to the local orthonormal frame.

| Operator | Shape | Definition | Reference |
| --- | --- | --- | --- |
| `b_membrane_2017(pre, r, s)` | 3×24 | `[ẽ_rr^m, ẽ_ss^m, 2ẽ_rs^m]` from Eq. (27a–c) with the five tying rows and `a_A..a_E`; then `covariant_to_local_mapping(j_loc(r,s))` | A Eqs. (17)–(27), Fig. 4 |
| `b_bending_2017(pre, r, s) → (B_b1, B_b2)` | 3×24 each | `B_b1` from `e_ij^b1 = ½(∂x_m/∂r_i·∂u_b/∂r_j + ∂x_m/∂r_j·∂u_b/∂r_i + ∂x_b/∂r_i·∂u_m/∂r_j + ∂x_b/∂r_j·∂u_m/∂r_i)`; `B_b2` from `e_ij^b2 = ½(∂x_b/∂r_i·∂u_b/∂r_j + ∂x_b/∂r_j·∂u_b/∂r_i)`; both displacement-based, both mapped point-wise | A Eqs. (7c)/(7d), (8a)/(8b) |
| `b_shear_mitc4(pre, r, s)` | 2×24 | `ẽ_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, `ẽ_st = ½(1+r)e_st^C + ½(1−r)e_st^D` from the four stored tying operators; then `shear_covariant_to_local` | DB84 Eq. (3) as reproduced in A p. 405 |
| `b_drill_membrane_2025(pre, r, s)` | 3×24 | Eq. (18) verbatim, with `cr_md`/`cs_md`, the `1/8` of Eq. (13c), the `[right, top, left, bottom]` edge order and `j0/j`; transformed with the **centre** mapping | B Eqs. (5), (10), (11a–b), (13b–c), (17a–b), (18), (19a–c), (21) |

**The covariant membrane tying rows** (`b_rr_a` … `b_rs_e`) are the *displacement-based*
covariant strains of Eq. (15)/(16) sampled at A(0,+1), B(0,−1), C(+1,0), D(−1,0),
E(0,0) — i.e. `x_r·u_r`, `x_s·u_s`, `½(x_r·u_s + x_s·u_r)` with
`∂x_m/∂r = x_r + s x_d`, `∂x_m/∂s = x_s + r x_d`, `∂u_m/∂r = u_r + s u_d`,
`∂u_m/∂s = u_s + r u_d`. Eq. (14)'s edge strains are algebraically the same values at
those points and are not implemented separately.

**The transverse shear's covariant→local metric normalization.** Paper A prints the
MITC4 assumed field in the *covariant* components `e_rt, e_st`; the metric
normalization to the local orthonormal frame is not printed. The design defines it
exactly, with no invented factor, as the contraction with the 3D dual basis:

```text
γ_α3 = 2 e_i3 (g^i · e_α) (g^t · e_3)          i = r, s
g^r = (g_s × g_t)/j,  g^s = (g_t × g_r)/j,  g^t = (g_r × g_s)/j,  j = det[g_r g_s g_t]
```

`e_tt = x_b·u_b = 0` exactly (the director rotation is orthogonal to the director), so
only the two transverse terms survive. Three oracles pin this mapping: (a) the flat
reduction to the standard Mindlin shear; (b) `test_t1a_shearing_patch_constant_stress_fig5_mesh`
and the shearing scenario of T1.6; (c) the closed-form uncorrected
`∫ B_γᵀ (G·h) B_γ dA` comparison of Requirement 1's third scenario. It is listed in §9
as a named risk precisely because no held equation prints it.

**The 2017 core is structurally blind to `θ_z`.** Because `u_b` depends on `θ_i × V_n^i`,
no operator above has an entry in slot `6i+5`. That is a property of the code, not a
convention, and `test_kinematics_drill_dof_is_theta_z_through_eq26_operator` asserts it
by checking that the 2017-only blocks have exactly zero in those columns while the
drill block does not.

### 2.3 Integration, and where the through-thickness resultants come from

**Surfaces: 2×2 Gauss** (`N_GAUSS = 4`, `ξ,η = ±1/√3`, weights 1) for every
contribution, including the drill — paper B §3.1: "the elements require only the 2×2
Gauss integration over the element surfaces"; paper A p. 410: "2 × 2 × 2 Gauss
integration over the element domain".

**Through thickness: the ABD moments.** The element is a *resultant* shell: the
repository's `ShellConstitutive` is the already-integrated `(A, B, D, C_s)`, and the
element consumes it as such. Paper A's `t`-integration is realized through those
moments. Writing the paper's decomposition as the coefficient triple
`[e^m, e^b1, e^b2]`, and using `z = t·h/2`, the energy
`∫ e(z)ᵀ C e(z) dz` expands to a 3×3 block matrix of through-thickness moments whose
entries are, with `e(z) = e^m + (2z/h)e^b1 + (4z²/h²)e^b2`:

```text
W_00 = ∫ C dz          = cm           (already in the interface)
W_01 = ∫ z C dz        = cb_coupling  (already in the interface)
W_11 = ∫ z² C dz       = cb           (already in the interface)
W_02 = ∫ z² C dz       = cb           (the same moment; ε_m–E2 coupling)
W_12 = ∫ z³ C dz       = 0            (odd moment: vanishes under the symmetric rule)
W_22 = ∫ z⁴ C dz       = cm/9         (paper A's 2×2 rule in t; see below)
```

with the **resultant strain triple** `[ε_m, κ, E2]`, `κ = (2/h)e^b1`,
`E2 = (4/h²)e^b2`, so that the ABD blocks act verbatim and no scaling factor appears in
the stiffness assembly. `W_22` is the paper's own under-integration: 2-point Gauss in
`t` gives `∫t⁴ dt → 2/9` where exact integration gives `2/5`, so
`∫z⁴C dz → (16/h⁴)·(h/2)(h/2)⁴(2/9)C = hC/9 = cm/9`. `cm/9` equals `h·cm_raw/9` for
every material, because `cm_raw := cm/h` in both constructors. For a homogeneous
section this is exact under the paper's rule; for a multi-ply laminate `cm/9` uses the
section's thickness-average membrane stiffness (the same smearing `cm_raw` already
documents) in the `E2–E2` block only, and the term it multiplies is second order in the
warping. Recorded as a named approximation in §9.

The full element stiffness is therefore

```text
K = Σ_{g=1..4} w_g [ Bᵀ W B ]_g  ·  sqrt_g        B = [B_m ; B_b1 ; B_b2]
  + Σ_{g=1..4} w_g [ B_γᵀ cs_uncorrected B_γ ]_g · sqrt_g
  + Σ_{g=1..4} w_g [ B_mdᵀ cm B_md ]_g           · sqrt_g      (drill, 2×2 surfaces only)
```

with `B_md` already carrying `j0/j`, so the measure of the drill term is the same
surface measure `sqrt_g = j(r,s)` as the membrane.

**The spec's "2×2×2 vs surface-only" discrimination is interpreted here.** For a
resultant element the only through-thickness content beyond the ABD moments is the
`t²` term, so the two references are: *2×2×2* = the three-term form above (with `W_22`);
*surface-only* = the two-term form with `W_22 = 0`. The difference is exactly
`B_b2ᵀ (cm/9) B_b2`, which is zero on a flat element (`∂x_b/∂r = 0`) and non-zero on a
warped one. The test therefore uses a **strongly warped quad** (node `z` offsets of the
order of the element's in-plane size) so that the `> 1e-4` relative separation the
scenario demands is real rather than borderline, and it asserts the separation so the
test cannot pass vacuously.

---

## 3. Constitutive path

`K`'s in-plane blocks use `cm` (membrane), `cb_coupling` (membrane–bending), `cb`
(bending) and `cm` again for the drill membrane contribution and the `E2` block; the
transverse-shear block uses `cs_uncorrected` (ADR-1). `cm_raw` is not needed by the
linear path. All five fields are consumed, and the interface is unchanged.

- **Uncorrected transverse shear** (ADR-1). Observable: the shear block is invariant to
  the material's applied correction factor, and equals the closed-form
  `∫ B_γᵀ (G·h) B_γ dA`. Interaction with the preserved-laminate invariant: the removal
  is `cs / applied_shear_correction_factor()`, and that factor is `1.0` for every
  multi-ply laminate, so no laminate's section stiffness is scaled. The single-ply
  laminate's scalar factor *is* removed from the element's shear block — and that is
  correct: the factor belongs to the material model's `Cs` (which is unchanged and
  still exposed), not to the element. `test_laminate_invariant_guard.py` and the T2C/T2G
  groups are the evidence that this does not move a laminate result.
- **Mid-surface only** (ADR-6).
- **Behaviour change is expected and is Tier 2's judgement.** The hybrid consumes
  `cs = k·G·h`; the faithful element consumes `G·h`, a `1/k = 1.2` change to the
  transverse-shear stiffness. That is a *deliberate* consequence of paper A p. 410, not
  an accident, and it is the single largest expected source of Tier-2 movement. Thin
  shells are dominated by membrane/bending, so the expected movement is small; a
  published-benchmark value that moves but stays inside `rtol=0.05` is a legitimate
  outcome (proposal §6.5), and a value that moves *outside* tolerance is a formulation
  bug, not a reason to re-introduce `k`.

---

## 4. Test design — the papers' own basic tests (Tier 1)

Recorded separately from the repository's cross-reference tests (Tier 2, §5), as
`openspec/config.yaml`'s design rule requires.

All Tier-1 Rust tests live in
`crates/aeroelast-core/src/elements/mitc4_plusd/tests.rs`
(module path `elements::mitc4_plusd::tests`), declared by `#[cfg(test)] mod tests;`
inside `mitc4_plusd.rs`; fixtures live in
`crates/aeroelast-core/src/elements/mitc4_plusd/tests/fixtures.rs`
(`…::tests::fixtures`). The spec's `cargo test -p aeroelast-core test_t1a_` /
`test_t1b_` substrings are unaffected by the module path.

### 4.1 Fixtures

**F-A/F-B — the star patch (one fixture, both papers).** Paper A Fig. 5 and paper B
Fig. 7(a) are the same mesh: a five-element star on the 10×10 square with 8 nodes.
`docs/formulations/mitc4plusd-2025-extract.md` records every coordinate as printed, so
**no equivalent-patch deviation is needed** and the fixture policy's fallback clause is
not invoked.

```rust
// fixtures.rs — star patch, CCW connectivity in the element's (r,s) convention
pub const STAR_NODES: [[f64; 3]; 8] = [
    [ 0.0,  0.0, 0.0],   // 0  C
    [10.0,  0.0, 0.0],   // 1  D
    [10.0, 10.0, 0.0],   // 2  A
    [ 0.0, 10.0, 0.0],   // 3  B
    [ 2.0,  2.0, 0.0],   // 4
    [ 8.0,  3.0, 0.0],   // 5
    [ 8.0,  7.0, 0.0],   // 6
    [ 4.0,  7.0, 0.0],   // 7
];
pub const STAR_ELEMS: [[usize; 4]; 5] = [
    [4, 5, 6, 7],   // central  (2,2) (8,3) (8,7) (4,7)
    [7, 6, 2, 3],   // top      (4,7) (8,7) A     B
    [6, 5, 1, 2],   // right    (8,7) (8,3) D     A
    [0, 1, 5, 4],   // bottom   C     D     (8,3) (2,2)
    [4, 7, 3, 0],   // left     (2,2) (4,7) B     C
];
```

The connectivity is written in the repository's counter-clockwise convention
(node 0 at `(ξ,η) = (−1,−1)`, …, node 3 at `(−1,+1)`), which the extract's prose order
is not; the fixture carries a self-test asserting each element's signed area is
positive and that the five areas sum to 100.

**F-A-BC — the T1.3 minimum boundary conditions (derived, not figure-read).**
Paper A's Fig. 5 publishes the mesh only; §4 says "the minimum number of constraints to
prevent rigid body motions". The fixture defines them explicitly and records that they
are derived:

```text
u_x = u_y = u_z = 0 at node C(0)
u_x = 0            at node D(1)      // removes the in-plane rotation about C
u_y = 0            at node A(2)      // removes the second in-plane tilt
θ_z = 0            at node C(0)      // removes the drill rotation
```

prescribed at the **exact constant-straining mode's values** (not zero), so the
constraints never conflict with the mode being tested.

**F-B-BC — the T1.6 boundary conditions (figure-read, used as printed).**
`docs/formulations/mitc4plusd-2025-extract.md` §Fig. 7(b)(c)(d):

```text
extension : B(3) u_x=u_y=u_z=0, θ_x=θ_y=θ_z=0 ; C(0) u_x=u_z=0 ; load direction +x at A(2)
bending   : B(3) all six = 0                   ; C(0) u_x=u_z=0, θ_y=0 ; load direction +y at A(2)
shearing  : B(3) all six = 0                   ; C(0) u_x=u_z=0, θ_x=θ_y=0
            interior nodes {4,5,6,7}: u_x=0, θ_x=θ_y=0 ; load direction +y at A(2)
θ_z       : free at every node except corner B, in all three cases
```

**F-W — the derived warped star patch.** The same 8 nodes with `z` offsets on the four
interior nodes (`z = ±0.5` alternating), used **only** for the scenario that requires a
warped patch (`test_t1b_drill_theta_z_free_except_corner_b`'s variant-(c) separation)
and for the `2×2×2` vs surface-only discrimination. It is documented as derived, never
used for a constant-stress assertion (a warped patch is not a valid constant-stress
patch test), and it is not a substitute for the figure-read mesh.

### 4.2 How the three patch tests' load magnitudes are derived

The papers publish no load magnitudes, and the extract confirms they are implied by the
constant-straining-mode requirement. The design derives them **independently of the
element stiffness**, so the tests are not vacuous:

1. Choose a constant stress state (for extension: `σ_xx = σ0`; bending: a constant
   moment `M0`; shearing: a constant transverse shear `Q0`).
2. Compute the **boundary tractions of that state** on the patch's outer boundary and
   integrate them with the bilinear shape functions, 2-point Gauss per boundary edge:

   ```text
   membrane : f_i += ∫_edge N_i (σ · n) dΓ
   bending  : f_{θ_n,i} += ∫_edge N_i M_n dΓ      (a constant moment field is in
                                                    equilibrium with Q_n = 0, so no
                                                    transverse nodal force appears)
   shearing : f_i += ∫_edge N_i (q · n) dΓ with q = the constant shear traction
   ```

   For a constant traction/moment on a straight edge of length `L` these integrate to
   `t·L/2` at each end, which the fixture computes by the Gauss rule rather than by
   hand so the same code path serves all three tests.
3. Assemble the 48-DOF star-patch system from the production element matrices
   (`fixtures::assemble_star_patch(&[Mitc4PlusDPrecomputed; 5]) -> DMatrix<f64>`), apply
   the F-A-BC/F-B-BC constraints with the mode's prescribed values, and solve
   `K_ff u_f = f_f − K_fc u_c` with a dense LU (48 DOFs — no sparse machinery, no
   benchmark probe, no N ≥ 32 concern).
4. **Pass condition**: for every element and every Gauss point, the recovered membrane
   stress equals the prescribed constant state to `1e-8` relative (absolute floor
   `1e-10·‖σ‖`) for each of the three independent states (`σ_xx` alone, `σ_yy` alone,
   `τ_xy` alone); the spread across all Gauss points is `≤ 1e-8·‖σ‖`; and for T1.6's
   *strong form*, every analytically-zero component is `≤ 1e-10` of the state's
   magnitude.

Because the load is the exact boundary traction of the constant state (not `K u_exact`),
the test is non-vacuous: it fails if the element's discrete strain at the constant mode
differs from the state, or if the integration is inconsistent.

### 4.3 Tier 1a — the 2017 core (paper A §4, pp. 410–411)

| Req | Test name (spec-fixed) | Fixture | Pass condition |
| --- | --- | --- | --- |
| T1.1 | `test_t1a_isotropy_element_orientation_and_node_sequence_invariant` | flat square, flat distorted quad, warped quad; 4 orientations incl. `π/2`; node sequences `(1,2,3,4)`, `(2,3,4,1)`, `(3,4,1,2)`, `(4,1,2,3)`, `(1,4,3,2)` | sorted symmetrised eigenvalues agree to `1e-10·λ_max`; permuted sequences reproduce the exactly-permuted entries to `1e-12·max|K|`; `uᵀKu` invariant to `1e-10` relative for a co-rotated field |
| T1.2 | `test_t1a_zero_energy_modes_single_unsupported_element_exactly_six` | single unsupported element, flat / flat-distorted / warped | count of `|λ| ≤ 1e-10·λ_max` is **exactly 6**; `‖K u_rb‖∞ ≤ 1e-10·λ_max·‖u_rb‖∞` for each of the six fields; the six smallest `|λ|` separated from the seventh by `≥ 1e-6·λ_max` |
| T1.3a | `test_t1a_membrane_patch_constant_stress_fig5_mesh` | F-A mesh + F-A-BC + derived tractions | §4.2 pass condition |
| T1.3b | `test_t1a_bending_patch_constant_curvature_fig5_mesh` | F-A mesh + F-A-BC + derived moment tractions | constant curvature/moment to `1e-8` relative, spread `≤ 1e-8` |
| T1.3c | `test_t1a_shearing_patch_constant_stress_fig5_mesh` | F-A mesh + F-A-BC + derived shear tractions | constant shear to `1e-8` relative, spread `≤ 1e-8`, no shear correction factor |
| T1.3d | `test_t1a_membrane_eq22_flat_tying_condition` | flat element (any shape) and a warped element; no patch | `ẽ_rs^m|bil == e_rs^m|bil` to `1e-14` absolute on the flat element; the same comparison on the warped element is `> 1e-6` relative (non-vacuity) |
| T1.3e | `test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18` | flat rectangle | `|x_d| ≤ 1e-14·max(|x_r|,|x_s|)`, `c_r = c_s = 0` to `1e-14`, `d = −1`, all five `a_*` zero to `1e-14`, Eq. (27) reduces term-by-term to Eq. (18) to `1e-17` absolute |

Supporting geometry tests, also Tier 1a (the spec puts them under the coefficients
requirement): `test_geometry_dual_basis_identities_eq11` (the four `m·x` identities to
`1e-12` and `m^r·n = m^s·n = 0` to `1e-12·|m||n|`, on flat-distorted, ruled-warped and
doubly-warped quads), `test_geometry_a_E_is_positive_eq27c` (`a_E > 0` for
same-sign `c_r, c_s`, equal to `2c_rc_s/(c_r²+c_s²−1)` to `1e-14`, with the negated
value asserted to differ by `> 1e-6` relative), and
`test_geometry_node_directors_reduce_to_n_vec_when_flat` ✱ (ADR-4: `V_n^i = n_vec` to
`1e-14` for all four nodes on a flat square and a flat distorted quad, and a non-vacuity
assertion that the four directors differ from each other by `> 1e-6` relative on the
warped patch F-W).

**The rigid-body field construction (T1.2/T1.5, and risk §9-3).** `fixtures::rigid_body_modes`
builds `u_i = t + ω × (x_i − x_c)` and `θ_i = ω` for all `i` (constant rotation
vector), expressed in the 6-DOF local layout. That field is exactly the paper's
rigid-body field for the 2017 core *with* a drilling component: `α_i = ω·V_1`,
`β_i = ω·V_2`, `θ_z,i = ω·V_n`, all constant, so `u_b` is the constant director
rotation and `θ^D = ω·V^D` is constant. A constant `θ^D` makes every drill edge
difference vanish, so `B_md u_rb = 0` **structurally** — which is why exactly six zero
modes survive the addition of the drilling DOF, and why paper B can state that using
one `V^D` satisfies the rigid-body-mode test.

### 4.4 Tier 1b — the 2025 six-DOF element (paper B §3.1, pp. 13–14)

| Req | Test name (spec-fixed) | Fixture | Pass condition |
| --- | --- | --- | --- |
| T1.4 | `test_t1b_spatial_isotropy` | flat + warped, ≥ 4 orientations, drill free | sorted eigenvalues to `1e-10·λ_max`; `uᵀKu` invariant to `1e-10` relative for a co-rotated field with a non-zero drilling component |
| T1.5 | `test_t1b_zero_energy_modes_exactly_six_with_drill_dof` | single unsupported element, flat / warped / flat-distorted, all 6 DOF/node free | exactly 6 zero eigenvalues; each rigid-body field `≤ 1e-12·λ_max` at `‖u‖ = 1`; separation `≥ 1e-6·λ_max`; the drilling DOF carries non-zero stiffness in at least one non-rigid mode |
| T1.6a | `test_t1b_strong_patch_extension_constant_and_zero_stress` | F-B mesh + F-B-BC | `σ_axial` equals the analytical constant to `1e-8` relative in every element; every analytically-zero component `≤ 1e-10·|σ_axial|` |
| T1.6b | `test_t1b_strong_patch_bending_constant_and_zero_stress` | F-B mesh + F-B-BC | constant-moment solution to `1e-8` relative; zero components `≤ 1e-10` of the state |
| T1.6c | `test_t1b_strong_patch_shearing_constant_and_zero_stress` | F-B mesh + F-B-BC | constant shear to `1e-8` relative; zero components `≤ 1e-10` of the state |
| T1.6d | `test_t1b_drill_theta_z_free_except_corner_b` | F-B mesh (flat) for (a)/(b); F-W (warped) for (c) | (a) vs (b) agree to `1e-10` relative (fixing `θ_z` at C is immaterial); (c) differs from (a) by `> 1e-8` relative on the warped patch |

### 4.5 The identity / kinematics lock (Requirement 1 and 2)

These are not "paper basic tests" but they are the mechanism that makes the identity
provable, so they are listed here rather than in Tier 2.

| Test name (spec-fixed unless marked ✱) | Content | Pass condition |
| --- | --- | --- |
| `test_identity_ke_lock_matches_2017_core_plus_2025_drill` | production `K` vs a **test-local reference implementation** written from the ingredient list (A Eqs. 1–27 with 7c/7d displacement-based, DB84 Eq. 3, B Eq. 26, the §2.3 moment matrix, plane-stress material, every factor at its printed value), on a flat square, a flat distorted quad, a ruled warped quad and a doubly warped quad | `max|K_prod − K_ref| ≤ 1e-10·max|K_ref|`, same bound on the membrane and transverse-shear blocks |
| `test_identity_drill_stiffness_comes_only_from_eq26` | production vs a variant whose Eq. (26) operator is the zero operator | difference non-zero on warped geometry; exactly symmetric; exactly zero on every **translational** row/column block (see §9-8 for the interpretation of the spec's wording); `\|u_rbᵀ K u_rb\| ≤ 1e-12·λ_max·‖u_rb‖²` for the six rigid-body fields |
| `test_identity_transverse_shear_uses_uncorrected_shear_modulus` | flat rectangular element, isotropic `(E, ν, h)`; transverse-shear block vs `∫ B_γᵀ (G·I) B_γ dA`, `G = E/(2(1+ν))` | agree to `1e-10` relative; the `5/6` value rejected by `> 1e-3` relative (asserted) |
| `test_identity_transverse_shear_invariant_to_shear_correction_factor` | ADR-1's two-`pre` construction, isotropic and single-ply-laminate, plus the two non-vacuity controls | shear blocks agree to `1e-10` relative; each matches the closed form to `1e-10`; the `k`-carrying control differs by `> 1e-3` relative |
| `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` | strongly warped quad; production vs the three-term reference vs the two-term (`W_22 = 0`) reference | matches the three-term reference to `1e-10`; differs from the two-term reference by `> 1e-4` relative (asserted) |
| `test_identity_bending_operator_matches_eq7c_eq7d` | flat + warped element with non-zero `α_i, β_i`; element's `e^b1, e^b2` vs a test-local Eqs. (7c)/(7d) evaluation including the `∂x_b·∂u_m` term of Eq. (8a) | agree to `1e-10` relative at every Gauss point; the local matrix is 24×24 with no condensed internal DOF |
| `test_kinematics_displacement_field_matches_eq1_to_eq3` | sampled `(r,s,t)` points, incl. `t = ±1` | `u` equals Eq. (3) to `1e-12` relative, `x` equals Eq. (1), and `θ×V_n = −V_2α + V_1β` with `α = θ·V_1`, `β = θ·V_2` |
| `test_kinematics_drill_dof_is_theta_z_through_eq26_operator` | pure rigid rotation about `V_n` with `θ_z ≠ 0`, `α = β = 0`; then a warped `θ_z` pattern violating the constant-strain state | energy `≤ 1e-12·λ_max·‖u‖²` for the rigid rotation (no penalty); non-zero energy through the Eq. (26) operator only on the warped pattern |
| `test_kinematics_local_matrices_are_24x24` | element type registration | `K`, `M`, `K_T` exactly 24×24 (576 entries), `f_int` exactly 24 long, drill at global slot `6i+5` |
| `test_identity_element_uses_midsurface_constitutive` | element path on a laminate | the element's constitutive equals `to_shell_constitutive()`; no offset variant on the path (ADR-6) |
| `test_identity_drill_operator_matches_eq18_term_by_term` ✱ | see §6 | see §6 |

✱ **Design-added test names.** The spec fixes the verification name for every requirement;
`test_identity_drill_operator_matches_eq18_term_by_term` (§6's Oracle 1/2) and
`test_geometry_node_directors_reduce_to_n_vec_when_flat` (§4.3, ADR-4) are names this
design adds because no spec requirement owns them. Every other name in §4.3–§4.5 is
the spec's exact string.

---

## 5. Test design — the repository's cross-reference theoretical tests (Tier 2)

Tier 2 is unchanged in *content*: the same test names, the same assertions, the same
tolerances. What changes is **which element they exercise** (S3, the flip) and the Rust
type the layout-bound tests construct (S2).

### 5.1 Rust — `cd crates && cargo test -p aeroelast-core` (baseline 120 passed / 0 failed)

- **T2A / T2B (`mitc4.rs` module `tests`)** — `test_ke_global_leaves_all_six_rigid_body_modes_free`,
  `test_ke_global_has_exactly_six_zero_modes`,
  `test_membrane_patch_reproduces_constant_strain_at_every_gauss_point`,
  `test_bending_patch_reproduces_constant_curvature_at_every_gauss_point`,
  `test_ke_global_is_symmetric`, `test_ke_global_is_positive_semidefinite`,
  `test_ke_local_flat_plate_parity`, `test_kt_fint_directional_derivative`,
  `test_kt_fint_directional_derivative_rotations`,
  `test_kt_fint_directional_derivative_with_drill_dofs`,
  `test_fint_linear_nonlinear_parity`, `test_kt_zero_matches_ke`.
  **How kept green:** at S2 these are *moved* (not copied) into
  `mitc4_plusd/tests.rs`, retargeted from `Mitc4Precomputed` to
  `Mitc4PlusDPrecomputed` and from `mitc4::…` to `mitc4_plusd::…`. The assertions are
  mechanical 6-DOF-layout equivalents (proposal T2I): the same fields, the same
  invariants, the same tolerances. They must pass against the new element **before** the
  flip, which is the strongest available evidence that the port is a port.
  `test_ke_local_eigenvalues_nonsymmetric` is *excluded* (explore §2G: it uses a
  deliberately non-symmetric `cb_coupling` and encodes the 6-DOF deviation).
- **T2I (`me_global*`, `body_load`, `k_sigma`, `centrifugal_prestress`,
  `compute_element_stress`)** — moved and retargeted with the same treatment.
- **T2J** (`test_drill_midside_derivatives_match_paper_eq11`,
  `test_drill_membrane_operator_*` ×5, `test_b_erc_*` ×5, `test_erc_covariant_to_local_*` ×2,
  `test_ke_local_erc_flat_rigid_body_and_symmetry`, `test_enhanced_drill_stiffness`,
  `test_drill_warping_moment_flat_element`, `test_ke_local_eigenvalues_nonsymmetric`)
  — **deleted at S4 with the hybrid**, never ported. The three `test_b_erc_*` tests that
  encode the ERC deviation are not Tier-2 acceptance; the underlying invariants (RB
  annihilated, exactly six zero modes) are covered by T1.2/T1.5 and the four T2A tests.
  Note that `test_drill_midside_derivatives_match_paper_eq11` is *not* vacuous on its own
  (it compares the stored derivatives against a `midside_shape_function` reconstruction),
  but its assertion is preserved in the new module as part of §6's drill verification, so
  the coverage is not lost — only the hybrid's copy is.
- `cargo test -p aeroelast-core materials::` — the laminate/composite unit tests; must
  pass with no edits (the preserved invariant).

### 5.2 Python — `python -m maturin develop --release && python -m pytest -m "not slow" -q` (baseline 345 passed / 2 failed / 2 skipped)

All names are the spec's, all assertions and tolerances are unmodified. The only thing
that changes is the Rust kernel behind them.

- **T2B/T2F** — `test_rust_composite.py::TestMITC4BatchSanity::{test_ke_shape_and_symmetry,
  test_ke_positive_semidefinite, test_me_shape_and_symmetry, test_me_positive_semidefinite}`;
  `test_material_suite.py::TestStiffnessProperties::{test_k_positive_definite, test_k_symmetric}`;
  `test_rust_assembler.py::{TestTangentStiffness, TestInternalForces, TestNewtonRaphsonConsistency}`
  (18 tests); `test_mass_matrix_validation.py` (15 tests).
  **How kept green:** shape/symmetry/PSD invariants are structural and hold for any
  positive-semidefinite element; `K_T(0) = K`, `f_int(0) = 0` and the NR consistency are
  preserved by construction because `compute_fint_global(pre, u, false)` returns
  `K u` and `compute_kt_global(pre, 0) == compute_ke_global(pre)` in the new module too.
- **T2C/T2D** — the analytical, convergence, modal, mass and CCX-parity groups of the
  spec's scenarios. **How kept green:** these are the *numerical* judge of the
  behaviour change of §3. Expected movement is small (thin shells are membrane/bending
  dominated); `test_shell_convergence.py::test_composite_laminate_gap_mesh_study`
  asserts no value and cannot fail.
- **T2E** — `test_ko2017_performance.py`, 31 cases, `rtol=0.05`. **How kept green:** a
  value that moves but stays inside tolerance is a legitimate outcome and is *not*
  reported as a regression; a value outside tolerance is a formulation bug. The
  known mis-sourcings flagged in proposal §6.6 (`test_3_2` SS rows vs Table 6,
  `test_3_3[dist]`'s hardcoded Table-12 value, `test_3_7[reg]`'s N=8/S4 cells, the
  `test_3_5` docstring/xfail mismatch) are corrected to the true paper cells as part of
  the documentation unit (WU11), because leaving them would let the gate rest on the
  wrong numbers. Correcting a *source* is not loosening a *tolerance*.
- **T2G — the laminate/composite preserved invariant.** `TestOrthotropicSinglePly`,
  `TestSymmetricLaminates`, `TestAsymmetricLaminates`, `test_composite_b_coupling.py`,
  `test_rust_composite.py`, `test_orthotropic_shell_parity.py`,
  `test_composite_beam_parity.py`, all with unmodified assertions. **This is the evidence
  that the retirement did not cost a capability**, not a code review of the deletion.
- **The two pre-existing failures** —
  `test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple` and
  `test_shell_convergence.py::test_in_bending_convergence` — are excluded from the pass
  condition, reported as pre-existing, and counted as a bonus if they start passing.
- **T2H** — `test_quad_elements.py` (plane QUAD4/8/9) is not applicable and stays
  untouched.

### 5.3 Guard tests (NEW, Python)

- `tests/test_mitc4plusd_traceability.py` — parses the ingredient table in
  `docs/formulations/shell-elements.md` §2 and the `mitc4_plusd.rs` module header, and
  asserts (a) every Requirement-1 ingredient key carries a non-empty paper/section/
  equation citation that resolves against `docs/formulations/mitc4plus-2017-extract.md`
  or `docs/formulations/mitc4plusd-2025-extract.md`, (b) no forbidden ingredient appears
  as implemented, (c) every cited ingredient is referenced by at least one existing test
  name matching `test_identity_* | test_kinematics_* | test_geometry_* | test_t1a_* |
  test_t1b_*`. **Source scope:** the *live* module list, fixed as `["mitc4_plusd.rs"]`
  before the retirement and unchanged after it; the retirement's own absence checks live
  in the Rust `test_retirement_*` tests, so the two do not overlap.
- `tests/test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged` —
  asserts `ShellConstitutive` still has exactly the five fields `cm`, `cb_coupling`,
  `cb`, `cs`, `cm_raw` (by parsing `materials/mod.rs`'s struct body), and that
  `Laminate` still exposes `to_shell_constitutive`, `to_shell_constitutive_with_offset`,
  `shear_correction_factor`, `abd_matrix_flat`.

---

## 6. Verification strategy for the drill operator (paper B Eq. 18)

The drill operator is the highest-risk transcription in the change: it has only ever
been dead code, and its only evidence so far is that it changed nothing *because it
could not*. The design verifies it against Eq. (18) with **three independent oracles**,
so that a wrong transcription cannot be silenced by editing the reference.

**Oracle 1 — a term-by-term reference written from Eq. (18) alone.**
`tests::drill::b_md_reference(pre, r, s)` is a test-local function that implements the
printed Eq. (18) with the paper's own symbols and *does not share code* with the
production operator: it recomputes `x_m^I` from the two edge node coordinates
(`x_m^I = 1/8 (x_i − x_{i+1})`), recomputes `x_r^I`, `x_s^I` from the geometry at the
edge mid-point, recomputes `j` as the scalar triple product, recomputes
`θ_i^D = θ_i·V^D` from the DOF, and accumulates the three components by summing over
the four edges in the paper's order. `test_identity_drill_operator_matches_eq18_term_by_term`
compares it to `b_drill_membrane_2025` entry by entry at **9 sample points** — the four
Gauss points, the four edge mid-points and the centre — with `1e-12` absolute, on a flat
square, a flat distorted quad and a ruled and doubly warped quad.

**Oracle 2 — the five known-risk deviations are asserted to be *rejected*.** The same
test builds five wrong variants and asserts each differs from the production operator by
more than a stated margin, so the pass cannot come from a loose comparison:

| wrong variant | margin | what it catches |
| --- | --- | --- |
| `V^D` replaced by the local `e3` | `> 1e-6` relative on a warped quad | the "single normal" rule (paper B p. 8) |
| `x_m^I` without the `1/8` of Eq. (13c) | factor 8, `> 1e-3` | the deleted implementation's error |
| edge order `[bottom, right, top, left]` | `> 1e-6` | the deleted implementation's second error |
| the edge difference's sign flipped | `> 1e-6` | a sign/indexing transcription |
| `θ_z` alone instead of `θ·V^D` | `> 1e-6` on a warped quad | the "rotations about `V^D`" rule |

**Oracle 3 — the paper's own tests, which do not touch the reference.**
`test_t1b_strong_patch_*`, `test_t1b_spatial_isotropy`,
`test_t1b_zero_energy_modes_exactly_six_with_drill_dof` and
`test_t1b_drill_theta_z_free_except_corner_b` are *load-driven or spectral*: they use
the assembled stiffness, never the drill reference. A production operator and a
reference that were both wrong in the same way would still fail these. The
structural invariant behind them is exact and worth stating: a constant `θ^D` field
gives zero drill strain for *any* coefficients, so the six rigid-body modes are in the
drill operator's null space by construction — which is the property the paper cites for
using one `V^D` per element.

**Why editing the reference cannot silence a wrong transcription.** The reference is
not a second copy of the production path: it is written from the printed equation with
independently recomputed geometric inputs, and it is pinned on three sides by
(a) Oracle 2's five asserted rejections, (b) Oracle 3's paper tests, which are
independent of it, and (c) `docs/formulations/mitc4plusd-2025-extract.md`, whose drill
section must carry an explicit `transcription-verified: <date> <page>` line that
`test_mitc4plusd_traceability.py` requires to be present for the Eq. (18) row. Editing
the reference to match the code therefore requires editing the extract's verification
line, which is a recorded, reviewable act — and Oracle 3 still fails if the operator is
wrong. If the reference and the operator disagree, the resolution is a **vision re-read**
recorded in the extract with a page citation, never a silent edit.

**One unresolved transcription point (flagged, settled by WU0).** The extract's Eq. (18)
block prints the edge difference as `(θ_{i-1}^D − θ_i^D)`
(`docs/formulations/mitc4plusd-2025-extract.md:67-70`) while the dead code's doc comment
describes it as `(θ_{i+1} − θ_i)`
(`crates/aeroelast-core/src/elements/mitc4.rs:866-869`, with the node-pair arrays
`edge_start = [1,2,3,0]`, `edge_end = [2,3,0,1]` at `:868-869`). These differ by the edge
traversal convention. The design's default is the **extract's printed form** (it is the
newer vision-read artifact), with the code's node-pair mapping retained only after WU0
confirms which end is which. This is the single most likely place for a wrong
transcription, and Oracle 2's sign-flip row plus T1.6 are what catch it.

---

## 7. File changes

| File | Change |
| --- | --- |
| `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | **NEW** — the faithful element (geometry, coefficients, B-operators, integration, K/f_int/K_T/M/body-load/k_sigma/stress, corotational machinery) |
| `crates/aeroelast-core/src/elements/mitc4_plusd/tests.rs`, `…/tests/fixtures.rs` | **NEW** — Tier-1a/1b tests, identity/kinematics/geometry tests, the drill reference, the star-patch and warped-patch fixtures, the retargeted T2A/T2B/T2I tests (moved at S2) |
| `crates/aeroelast-core/src/elements/mod.rs` | one line: `pub mod mitc4_plusd;` |
| `crates/aeroelast-core/src/assembly/topology.rs` | **unchanged** (ADR-4 option B: the element computes its own per-node directors, so `MeshTopology` gains no `node_normals` field and `MeshTopology::new`'s signature stays as at `topology.rs:72`) |
| `crates/aeroelast-core/src/assembly/assembler.rs` | `PrecomputedElem::Quad` payload → `Mitc4PlusDPrecomputed`; both `Mitc4Precomputed::new` sites (`assembler.rs:161-167` and `:252-254`); `mitc4::…` → `mitc4_plusd::…`; `build_constitutive_mitc4_plusd`; `extract_elem_disp_24`'s return type. **No nodal-normal data is threaded through** (ADR-4) |
| `crates/aeroelast-core/src/materials/mod.rs` | `ShellConstitutive::transverse_shear_uncorrected` (additive accessor) |
| `crates/aeroelast-core/src/materials/laminate.rs` | `Laminate::applied_shear_correction_factor` (additive accessor) |
| `crates/aeroelast-py/src/elements.rs` | swap the MITC4 entry points' internals to `mitc4_plusd`; **Python-facing signatures, shapes and symbol names unchanged**; `e_mod` kept and ignored with a comment |
| `crates/aeroelast-py/src/assembler.rs` | `MaterialSpec::Composite` gains `applied_shear_correction` from `Laminate::applied_shear_correction_factor()` |
| `crates/aeroelast-core/src/elements/mitc4.rs` | **deletions only, at S4** — the ERC block, `beta_w`/`BETA_W`, `k_drill`/`drilling_scale`, the SRI split, the rotation bubble and its `GpBubble`/`b_kappa_bubble`/`b_gamma_mitc4_plus`, the Eqs. (18)–(19) membrane as the live field, the shear-correction usage, the dead code (`b_m_standard`, `green_lagrange_strain`, `compute_b_l`, `compute_membrane_stress`, `Mat26`/`Vec26`), the T2J tests, and the module header's deviation text. `b_md_mitc4_plus`/`drill_midside_shape_derivatives` are **revived in `mitc4_plusd.rs`** as `b_drill_membrane_2025`/`drill_midside_shape_derivatives` (kept name) |
| `docs/formulations/mitc4plusd-2025-extract.md` | WU0: Eqs. (1)–(16), the `θ_z` interpolation, the `h̃_m` functions, the exact Eq. (18) edge convention, the `transcription-verified` line |
| `docs/formulations/mitc4plus-2017-extract.md` | WU0: Eq. (21)'s leading term (F2), the Eq. (27) derivation note (F1) |
| `docs/formulations/shell-elements.md` | §2 rewritten to describe what runs (the MITC4+/D), with the ingredient table the traceability test parses; §4.2's "pending MITC4/D drill" resolved |
| `docs/validation-matrix.md` | refresh the pass/fail columns and correct the four flagged mis-sourcings |
| `tests/test_mitc4plusd_traceability.py`, `tests/test_laminate_invariant_guard.py` | **NEW** guard tests (§5.3) |
| `src/aeroelast/core/assembler.py` | **unchanged** (`_FAMILY_PROPERTIES[SHELL] = (6, 3)` is frozen) |

---

## 8. Retirement sequencing, and how the tree stays green

| Stage | What happens | Why the tree is green afterwards | Revert |
| --- | --- | --- | --- |
| **S1** | `mitc4_plusd.rs` + its tests land. The module is compiled and tested but **not dispatched to**; the hybrid is untouched and remains the live element. | Tier 2 measures the unchanged hybrid, exactly as before; the new module's tests are additive. | delete the module (nothing else references it) |
| **S2** | The layout-bound Tier-2 Rust tests (T2A/T2B/T2I) are **moved** from `mitc4.rs`'s test module into `mitc4_plusd/tests.rs` and retargeted. The hybrid loses its copies but keeps its code. | The moved tests must pass against the new element *before* the flip — the strongest port evidence available. The hybrid is unchanged code with no tests, and the suite is green. | revert the move |
| **S3 — the flip** | `ElemType::{Mitc4,Mitc4Composite}` dispatch, `update_reference`, `extract_elem_disp_24`, the PyO3 MITC4 kernels and `build_constitutive_mitc4_plusd` all move to `mitc4_plusd`; the ADR-1 factor channel is wired. **Gate check: Tier 1 (all) + Tier 2 (all) + the laminate groups + the baselines.** | Tier 1 passes because S1 proved it; Tier 2 passes because the gate says so; the hybrid is now dead but present. | **one line per dispatch site** — the hybrid path returns |
| **S4 — the retirement** | The hybrid's deviation surfaces, dead code and T2J tests are deleted from `mitc4.rs`; the module header is rewritten. **Gate check: Tier 1 + Tier 2 + the static checks + the PyO3 surface check, re-run.** | The live path is the new module; nothing the tests exercise is deleted. | revert the retirement commit; the hybrid returns |

**Ordering rules that make this safe.**

1. **Nothing is deleted before the gate passes.** S4 cannot start until S3's gate is
   green, and S4's own gate is re-run after the deletion.
2. **The flip and the deletion are separate units.** A Tier-2 regression at S3 is
   attributable to the dispatch alone; a regression at S4 is attributable to the
   deletion alone. Neither can hide inside the other.
3. **The PyO3 surface never moves.** The retirement removes Rust code that the bindings
   call into, but the bindings' Python-facing names, argument lists, return shapes
   (`[f64; 576]` / `[f64; 24]`) and the Python family table are frozen. The check is
   re-run at S4.
4. **The process gate is a checklist, not a test.** The spec's
   `test_retirement_*` tests must not be *present-and-green* while any Tier-1 or Tier-2
   acceptance test is failing; this is enforced by running the gate commands **before**
   the retirement tests are accepted, and the fact is recorded in the change. A test
   cannot observe the order in which other tests were run, so this is stated as a
   process obligation in the work unit and in §9.

---

## 9. Risks and mitigations

| # | Risk | Mitigation |
| --- | --- | --- |
| 1 | **G4 — the Eq. (18) transcription is wrong.** It was never exercised, so "it changed nothing" is not evidence. | §6's three oracles, the five asserted rejections, the extract's `transcription-verified` line required by the traceability test, and WU0's vision re-read **before** the operator is wired. The unresolved `θ_{i-1} − θ_i` vs `θ_{i+1} − θ_i` point is settled in WU0. |
| 2 | **G7 — the offset open point.** `to_shell_constitutive_with_offset` shifts the ABD reference plane while the paper's kinematics puts it at the mid-surface. | ADR-6: the element path uses the mid-surface form only; `test_identity_element_uses_midsurface_constitutive` asserts it; the static check that the path references no offset variant; the incompatibility is recorded as a derivation the papers do not provide, not as an untested assumption. |
| 3 | **The zero-energy count must stay exactly 6 with the drill DOF.** A drilling mechanism (a seventh zero mode) or a spurious mode would fail T1.2/T1.5. | The rigid-body field is built as a *constant* `θ = ω` (all nodes), so `α, β, θ_z` are constant and every drill edge difference vanishes **structurally**; `B_md u_rb = 0` identically. The test asserts exactly 6, the separation `≥ 1e-6·λ_max`, and non-zero drill stiffness in at least one non-rigid mode — so both failure directions are caught. |
| 4 | **The process-gate nature of the retirement ordering.** No test can prove the deletion happened after the gate; a green suite after S4 does not prove S3 was green. | §8's ordered stages with two explicit gate checks, the rule that S4 is a revertible pure deletion, and the recording of the gate runs in the change. The `test_retirement_removes_hybrid_deviation_surfaces` test plus the static `grep` are the *after* evidence; the *before* evidence is the recorded S3 gate run. |
| 5 | **Removing the SRI changes behaviour.** The change is not purely additive: the 2017 core itself changes the numbers, and the uncorrected shear changes the transverse-shear stiffness by `1/k`. | Tier 2 is the judge (proposal §6.5): a value that moves inside its published tolerance is legitimate; outside is a formulation bug. The two pre-existing failures are excluded and reported. T2C/T2G's laminate groups are the preserved-invariant evidence. |
| 6 | **The covariant transverse-shear metric normalization is not printed in paper A.** | §2.2 defines it exactly from the 3D dual basis with no invented factor, and pins it with three oracles (flat reduction to Mindlin, the shearing patch tests, the closed-form `G·h` comparison). If it is wrong, T1.3c and T1.6c fail loudly. |
| 7 | **The `t⁴` moment for multi-ply laminates** is the section's thickness-average `cm/9` rather than the true `∫z⁴C dz`. | Named in §2.3; exact for the homogeneous/single-ply materials the paper's element is defined for; the term it multiplies (`E2`) is second order in the warping; and the `E2–E2` block is exactly the block the "surface-only" discrimination exercises, so the approximation is visible to a test rather than hidden. |
| 8 | **The spec's "zero on every row/column block that does not involve the drilling DOF"** cannot hold literally on warped geometry, because `θ^D = θ·V^D` couples the drill operator to `θx`/`θy` when `V^D ≠ e3` — and the paper *requires* those rotations. | Recorded interpretation: the assertion is read as **exactly zero on every translational row/column block**, plus exactly zero on *all* rotation blocks other than the drill's on flat geometry (where `V^D = e3`). The θx/θy coupling is a paper requirement (`rotations about V^D`), not a deviation, and Oracle 2 of §6 (inside `test_identity_drill_operator_matches_eq18_term_by_term`) asserts that the `θ_z`-only variant is *rejected*. The exact spec clause, its requirement/scenario, and the proposed replacement wording the parent may adopt are stated in **§9.1**. |
| 9 | **The nonlinear path is not a paper-faithful MITC4+/D nonlinear formulation.** | Bounded scope: `nonlinear = false` returns `K u` exactly; `nonlinear = true` keeps the repository's total-Lagrangian covariant correction applied to the new operators. The oracle is T2B's consistency set (`test_fint_linear_nonlinear_parity`, `test_kt_fint_directional_derivative{,_rotations,_with_drill_dofs}`, `test_kt_zero_matches_ke`), which caught a partial wiring once before. A full paper-faithful nonlinear derivation is explicitly deferred. |
| 10 | **Per-node directors change the reference geometry of every warped element**, so `x_b` becomes non-constant and `sqrt_g`, `element_area` and the mass matrix move slightly. | Expected and intended (it is Eq. (8a)). `test_me_global_total_translational_mass_is_rho_h_a` and the mass groups are the oracle for the mass; the geometry self-tests (areas, dual-basis identities) are the oracle for the geometry. |
| 10a | **The directors are element-local (ADR-4 option B)**, so on a curved mesh adjacent elements disagree slightly about a shared node's director, i.e. the director field is discontinuous. | Named in ADR-4 as the price of option (B) and bounded: on flat geometry every node gets `n_vec` exactly, so no Tier-1 result is affected; option (A) (a `MeshTopology` field, blast radius stated in ADR-4) is the named follow-up if a Tier-2 result shows the discontinuity matters. `test_geometry_node_directors_reduce_to_n_vec_when_flat` pins the flat case and asserts non-vacuity on the warped patch. |
| 11 | **Budget.** The change is well over 700 changed lines in total (new element + tests + a deletion). | The session's accepted delivery is `single-pr (size:exception)` with a 700-line review budget. §10 therefore sizes **every work unit** to ≤ ~350 changed lines — i.e. each unit is reviewable inside the 700-line budget — and orders the units so a reviewer can stop at any unit with a green tree. The honest total is stated in §10's ledger (~2 700 added, ~900 removed at S4, ≈1 800 net) and exceeds a single 700-line PR: that is precisely what the accepted `size:exception` covers, and no chain strategy is invented here (per `ask-on-risk`). If the reviewer objects to the total, the natural split points are S1–S2 (additive element + Tier 1) / S3 (the flip) / S4 (the retirement) / WU11 (docs) — recorded as a fallback, not chosen. |

---

### 9.1 The spec conflict, stated precisely (no spec edit made)

**Where.** The spec's requirements are unnumbered headings; the conflicting clause is in
the **first of the 16 `### Requirement:` headings**, `### Requirement: Element identity
and formulation` (`openspec/changes/mitc4plusd-faithful/specs/mitc4plusd-element/spec.md:55`),
inside `#### Scenario: The drilling stiffness comes only from the 2025 drill-membrane
strain, not from a penalty` (`spec.md:91`). The exact clause is at `spec.md:95`:

> THEN the difference `K(operator) − K(operator := 0)` is non-zero on warped geometry
> (the operator is live), is exactly symmetric, is exactly zero on every row/column block
> that does not involve the drilling DOF, and the six rigid-body fields still have strain
> energy exactly zero …

**Why it cannot hold literally.** The requirement itself (spec `Requirement: Element
identity and formulation`, item 7 at `spec.md:66`) adopts paper B **Eq. (18)** as the
drill operator, and paper B p. 8 requires the drill rotation to be `θ^D = θ·V^D` —
*"using the rotations about `V^D`, leads to satisfying the element isotropy test and rigid
body modes test, and even when the element is not flat"*
(`docs/formulations/mitc4plusd-2025-extract.md:76-80`, with the operator's `θ_{i-1}^D −
θ_i^D` at `:67-70`). On a warped element `V^D ≠ e3`, so `θ^D` has non-zero components
along the local `e1`/`e2`, and the drill operator's rows/columns are non-zero in the
`θx`/`θy` slots. Those are exactly "row/column block[s] that do not involve the drilling
DOF" in the clause's literal sense, and the clause says they must be **exactly zero**.
The two statements are therefore mutually exclusive on warped geometry: satisfying the
clause would require rejecting paper B p. 8, and satisfying paper B p. 8 necessarily
violates the clause. (On flat geometry there is no conflict: `V^D = e3`, `θ^D = θ_z`, and
the drill operator touches only slot `6i+5`.)

**Proposed replacement wording (for the parent to apply or reject — this design does not
edit the spec).** Replace the single clause with:

> is exactly symmetric, is exactly zero on every **translational** row/column block, and
> is exactly zero on every rotation row/column block **other than the drill's own**
> whenever the element is flat (so that `V^D = e3` and the operator touches only slot
> `6·i + 5`); on warped geometry the drill DOF's coupling to the nodal rotation
> components along `V^D` is required by paper B p. 8 and is therefore permitted, while
> every translational block remains exactly zero and the six rigid-body fields still have
> strain energy exactly zero (`|u_rbᵀ K u_rb| ≤ 1e-12 · λ_max(K) · ‖u_rb‖²`)

and, to keep the clause non-vacuous, add the assertion the design already implements:

> the `θ_z`-only variant of the operator (using `θ_z` in place of `θ·V^D`) is asserted to
> differ from the production operator by more than `1e-6` relative on a warped quad

The last sentence is the separation that makes the amended clause testable rather than a
loophole: `test_identity_drill_stiffness_comes_only_from_eq26` asserts the translational
blocks are exactly zero, and `test_identity_drill_operator_matches_eq18_term_by_term`
(Oracle 2, §6) asserts the `θ_z`-only variant is rejected. Until the parent amends the
spec, the design's recorded interpretation is the amended wording above and the tests
implement it; the literal clause is not satisfiable together with paper B p. 8.

---

## 10. Work-unit slicing

Ordered. Each unit states its own verification and is sized so a reviewer can follow it
independently (≤ ~350 changed lines unless noted). The session's accepted delivery is a
single PR with `size:exception`; the units are the review slices.

| WU | Unit | Touches | Verification | Size |
| --- | --- | --- | --- | --- |
| **WU0** | **Paper-B/paper-A vision re-read** (precondition, no code): Eqs. (1)–(16), the `θ_z` interpolation, the `h̃_m` functions, the exact Eq. (18) edge convention and sign; paper A Eq. (21)'s leading term (F2) and the Eq. (27) derivation (F1); add the `transcription-verified: <date> <page>` lines | `docs/formulations/mitc4plusd-2025-extract.md`, `docs/formulations/mitc4plus-2017-extract.md` | the two extracts are self-consistent (Eq. (27) → Eq. (18) reduction still holds; Eq. (21) now yields Eq. (22)); the drill convention is written down; no code touched | ~180 lines of Markdown |
| **WU1** | **Fixtures as code**: star patch nodes/connectivity, the derived T1.3 BC set, the figure-read T1.6 BC sets, the warped star patch, the boundary-traction loader, the 48-DOF dense patch assembler, the geometry self-tests | `mitc4_plusd/tests/fixtures.rs` | fixture self-tests: 5 elements / 8 nodes, positive signed areas summing to 100, no duplicate coordinates, the central element's connectivity equals the extract's, each BC set constrains exactly the six rigid-body modes | ~200 lines |
| **WU2** | **Geometry, coefficients, directors, kinematics**: `compute_characteristic_vectors`, `compute_node_directors` (ADR-4, including the flat-reduction and warped non-vacuity assertions), `compute_membrane_coefficients_2017`, `compute_j3d_enriched`, the `j`/`j0` triple product, the `Mitc4PlusDPrecomputed` struct and constructor | `mitc4_plusd.rs`, `mitc4_plusd/tests.rs` | `test_geometry_flat_rectangle_zero_xd_and_zero_coefficients_eq27_reduces_eq18`, `test_geometry_dual_basis_identities_eq11`, `test_geometry_a_E_is_positive_eq27c`, `test_geometry_node_directors_reduce_to_n_vec_when_flat`, `test_kinematics_displacement_field_matches_eq1_to_eq3` | ~300 lines |
| **WU3** | **B-operators**: `b_membrane_2017` (Eq. 27), `b_bending_2017` (7c/7d), `b_shear_mitc4` + `shear_covariant_to_local`, `b_drill_membrane_2025` (Eq. 18) + the five tying rows + the four shear tying operators | `mitc4_plusd.rs`, `mitc4_plusd/tests.rs`, `…/tests/drill.rs` | `test_identity_bending_operator_matches_eq7c_eq7d`, `test_identity_drill_operator_matches_eq18_term_by_term` with all five asserted rejections, `test_t1a_membrane_eq22_flat_tying_condition`, the `cr_md`/`cs_md` distinction test | ~350 lines |
| **WU4** | **K assembly**: `resultant_moment_matrix`, `compute_ke_local`/`compute_ke_global`, the drill contribution, the `cs_uncorrected` wiring; the test-local reference implementation | `mitc4_plusd.rs`, `mitc4_plusd/tests.rs` | `test_identity_ke_lock_matches_2017_core_plus_2025_drill`, `test_identity_drill_stiffness_comes_only_from_eq26`, `test_identity_transverse_shear_uses_uncorrected_shear_modulus`, `test_identity_transverse_shear_invariant_to_shear_correction_factor`, `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only`, `test_kinematics_local_matrices_are_24x24` | ~340 lines |
| **WU5** | **Everything the assembly layers call**: `compute_fint_global` (linear + the bounded nonlinear path), `compute_kt_global`, `compute_me_global`/`compute_me_composite_global`, `compute_body_load_global`, `compute_k_sigma_global`, `compute_centrifugal_prestress`, `compute_element_stress`, `build_t24`/`transform_to_global`, the corotational helpers | `mitc4_plusd.rs` | `test_kt_zero_matches_ke`, `test_fint_linear_nonlinear_parity`, the directional-derivative triple, `test_kt_fint_directional_derivative_with_drill_dofs`, the mass invariants | ~300 lines |
| **WU6** | **Tier 1a tests** | `mitc4_plusd/tests.rs` | T1.1, T1.2, T1.3a–c (the three Fig. 5 patch tests), the rigid-body-mode fixture | ~300 lines |
| **WU7** | **Tier 1b tests** | `mitc4_plusd/tests.rs` | T1.4, T1.5, T1.6a–d (the three strong-form Fig. 7 patch tests, the `θ_z` variant triple) | ~330 lines |
| **WU8** | **Move the layout-bound Tier-2 Rust tests onto the new element** (S2): T2A/T2B/T2I moved and retargeted; the hybrid's copies removed | `mitc4.rs` (tests only), `mitc4_plusd/tests.rs` | the moved names pass against `Mitc4PlusDPrecomputed`; `cargo test -p aeroelast-core` is 120+ passed / 0 failed; each name exists exactly once | ~300 lines (mostly moves) |
| **WU9** | **The flip** (S3): dispatch, `update_reference`, `extract_elem_disp_24`, the PyO3 MITC4 kernels, `build_constitutive_mitc4_plusd`, the ADR-1 factor channel (`ShellConstitutive::transverse_shear_uncorrected`, `Laminate::applied_shear_correction_factor`, `MaterialSpec::Composite.applied_shear_correction`). **`topology.rs` is untouched** (ADR-4 option B) | `assembler.rs`, `aeroelast-py/src/{elements,assembler}.rs`, `materials/{mod,laminate}.rs` | **the S3 gate**: all Tier 1, all Tier 2, the laminate/composite groups, both baselines, the maturin rebuild, and the PyO3 surface check | ~180 lines |
| **WU10** | **The retirement** (S4): delete the hybrid's deviation surfaces, dead code and T2J tests; rewrite the module header; add `test_retirement_removes_hybrid_deviation_surfaces` and the `test_retirement_*` gate tests | `mitc4.rs` | **the S4 gate**: Tier 1 + Tier 2 + the static checks + the PyO3 surface check, re-run; the supporting `grep -nE 'k_drill\|drilling_scale\|compute_ke_local_erc\|beta_w\|hg_stiffness_factor\|cm_normal\|b_m_standard\|green_lagrange_strain\|compute_b_l\|compute_membrane_stress\|Mat26\|Vec26'` returns no match in the dispatched module | net-negative lines |
| **WU11** | **Docs and guard tests**: `shell-elements.md` §2 rewritten with the ingredient table; `validation-matrix.md` refreshed and the four mis-sourcings corrected; the two Python guard tests | `docs/formulations/shell-elements.md`, `docs/validation-matrix.md`, `tests/test_mitc4plusd_traceability.py`, `tests/test_laminate_invariant_guard.py` | the three traceability scenarios and the laminate-guard scenario pass; `python -m pytest -m "not slow" -q` returns the baseline | ~280 lines |

**Budget accounting (single PR, accepted `size:exception`, 700-line review budget).**
Each unit is sized so that it is reviewable **inside the 700-line budget**; the table
below is the design's honest estimate of the cumulative diff, so the delivery decision
rests on real numbers rather than on an optimistic total.

| Unit | Added (est.) | Removed (est.) | Unit size |
| --- | --- | --- | --- |
| WU0 | 120 | 20 | 140 |
| WU1 | 200 | 0 | 200 |
| WU2 | 300 | 0 | 300 |
| WU3 | 350 | 0 | 350 |
| WU4 | 340 | 0 | 340 |
| WU5 | 300 | 0 | 300 |
| WU6 | 300 | 0 | 300 |
| WU7 | 330 | 0 | 330 |
| WU8 | 0 (moves: delete + re-add the same test bodies) | 0 | ~0 net |
| WU9 | 180 | 40 | 220 |
| WU10 | 40 | ~940 | net-negative |
| WU11 | 280 | 20 | 300 |
| **Total** | **≈ 2 740** | **≈ 1 020** | **≈ 1 720 net** |

The per-unit bound (≤ ~350) is what the 700-line review budget governs: **no work unit
exceeds it**, so each unit is a self-contained review slice with a green tree and its own
named verification (the `Verification` column above). The cumulative ≈1 720 net lines do
**not** fit inside a single 700-line PR; that is exactly the case the accepted
`single-pr (size:exception)` covers, and this design does not invent a chain strategy
(`ask-on-risk`). The largest additive blocks are the element's operators (WU3) and its
Tier-1 tests (WU6/WU7); if the reviewer objects to the total, the recorded split points
are S1–S2 / S3 / S4 / WU11.

**Ordering constraints.** WU0 before WU3 (the drill operator needs the settled Eq. (18)
convention and the `h̃_m` transcription). WU1 before WU6/WU7. WU2 before WU3 before WU4
before WU5. WU8 before WU9 (the moved tests must pass against the new element before the
dispatch moves). WU9 before WU10. WU11's traceability test targets the live module
(`mitc4_plusd.rs`) and is green from WU8 onwards; the four benchmark-sourcing
corrections in WU11 must land **before** the S3 gate is judged, because a mis-sourced
cell could otherwise let a formulation error pass.

**The test-first instruction (proposal risk §6.1).** `strict_tdd: false`, so the harness
will not enforce test-first. Every unit from WU2 onwards is written **test-first within
the unit**: the named test is added and observed failing before the production function
is written. WU6/WU7's Tier-1 tests are written and failing before WU4/WU5's assembly is
complete; this is the change's only real safety net and is carried explicitly into
`tasks.md`.

---

## 11. Decided, and not decided

**Decided by this design** (the six questions the phase had to settle): ADR-1 the
uncorrected shear mechanism; ADR-2 the internal structure (new module,
`Mitc4PlusDPrecomputed`, unchanged 24-DOF layout, `θ×V_n` rotation convention,
`θ_z` carried but unconsumed by the 2017 core); ADR-3 `V^D`/`j`/`j0` and the `c_r`/`c_s`
disambiguation; **ADR-4 (Rev. 2) the per-node directors and thicknesses, re-derived from
the real code: the element computes its own directors locally (option B) and
`MeshTopology` is left unchanged**; ADR-5 the retirement sequencing; ADR-6 the
mid-surface restriction. Plus F1/F2 (the Eq. (27) → Eq. (18) reduction and the
**confirmed** Eq. (21) print defect) and §2.3's through-thickness moment matrix.

**Could not be decided here — carried as explicit open items.** Each states what would
settle it and the work unit that settles it; none is left as a silent assumption.

| # | Open item | What would settle it | Settled by |
| --- | --- | --- | --- |
| 1 | **The composite per-element PyO3 batch entry point's factor channel.** `batch_ke_mitc4_composite` (`crates/aeroelast-py/src/elements.rs:532-540`) takes `(coords, cm_flat, b_coupling_flat, cb_flat, cs_flat, thickness, e_equiv)` — **no factor channel** — and its Python callers are tests only (no `src/aeroelast` caller). The design passes `applied_shear_correction = 1.0` there, so the element consumes the caller's `cs` verbatim: literally "no factor of its own", but a *single-ply* laminate routed through that test-facing entry point keeps the material model's `0.75`. The **production** composite path (`MeshAssembler` + `MaterialSpec::Composite`) is fully fixed. | Either the Python caller passes the uncorrected `cs`, or a factor argument is added — which is an out-of-scope PyO3 signature change. The `[f64;576]`/`[f64;24]` shapes and every `#[pyfunction]` name stay frozen either way. | WU9 (the flip) documents and wires the production path; the test-facing entry point keeps the recorded `1.0` until PyO3 changes come into scope. |
| 2 | **The Eq. (18) edge convention:** `θ_{i-1}^D − θ_i^D` as printed in the extract (`docs/formulations/mitc4plusd-2025-extract.md:67-70`) vs the dead code's `θ_{i+1} − θ_i` (`crates/aeroelast-core/src/elements/mitc4.rs:866-869`). The default is the extract's printed form with the code's node-pair mapping (`edge_start = [1,2,3,0]`, `edge_end = [2,3,0,1]`). | A vision re-read of paper B Eq. (16b)/(18) establishing which end of edge `I` the difference starts from; recorded in the extract with a page citation. | **WU0** (before WU3); Oracle 2's sign-flip row and T1.6 are the fallback detector. |
| 3 | **Paper B's Eqs. (1)–(16)** — including the full `h̃_m` mid-side functions and the `θ_z` interpolation — are not transcribed, so the drill reference's inputs are reconstructed from the extract's Eq. (18) block and the code's doc comments. | A vision read of B pp. 2–8 adding the equations to `docs/formulations/mitc4plusd-2025-extract.md` (the G4 precondition). | **WU0**; WU3 consumes it. |
| 4 | **The covariant transverse-shear metric normalization** is defined here (exact 3D dual basis, no invented factor, §2.2) but is not printed in paper A. | WU0's re-read of paper A p. 405 looking for the paper's own normalization; if it is printed, the design's definition is replaced by it. Otherwise the three oracles stand as the evidence. | **WU0** (search) / **WU3** (implement + oracles). |
| 5 | **The nonlinear path's fidelity.** The linear formulation is paper-faithful; the nonlinear path keeps the repository's total-Lagrangian covariant correction applied to the new operators. | A paper-faithful nonlinear MITC4+/D derivation, which neither paper in scope provides for this repository's UL formulation. | **WU5** implements and bounds it; T2B's consistency set (`test_fint_linear_nonlinear_parity`, the directional-derivative triple, `test_kt_zero_matches_ke`) is the oracle. |
| 6 | **The `t⁴` moment for multi-ply laminates** uses the section's average `cm/9` instead of the true `∫z⁴C dz`. | A sixth `ShellConstitutive` field (forbidden by the spec) or a new `Laminate` accessor plumbed into the element. | **WU4** names the approximation and its bound; `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` makes it visible to a test. |
| 7 | **A tapered per-node thickness** (`a_i` varying) is representable in the struct but has no input and no test. | A mesh-level per-node thickness input, which does not exist in this change. | **WU2** stores `a_i: [f64;4]` in the paper's per-node form so the future change is one line; no taper case is tested. |
| 8 | **The reading of the spec's "row/column block that does not involve the drilling DOF"** (spec.md:95) — the conflict and the proposed replacement wording are in **§9.1**. | The parent's decision to amend the spec (or not). The design's recorded interpretation is the amended wording and the tests implement it. | **WU4** asserts the interpretation (`test_identity_drill_stiffness_comes_only_from_eq26` + Oracle 2's `θ_z`-only rejection). |
| 9 | **The offset deferral** (ADR-6 / spec G7): `Laminate::to_shell_constitutive_with_offset` shifts the ABD reference plane while paper A's kinematics puts the reference surface at the mid-surface. | A derivation of the extra `u_b`/`x_b` offset term — neither paper provides it. | **WU9**'s static check (the element path references no offset variant) plus `test_identity_element_uses_midsurface_constitutive`. |
| 10 | **The element-local director field** (ADR-4 option B) is discontinuous across element boundaries on a curved mesh. | A Tier-2 result that moves because of it, or a decision to add the `MeshTopology` field (blast radius stated in ADR-4). | **WU2** pins the flat reduction and the warped non-vacuity; option (A) is the named follow-up. |

---

## 13. Known limitations (measured, with the evidence)

These are the limitations that the test suite now carries EXPLICITLY, after the cleanup
that took `cargo test -p aeroelast-core` from 169 passed / 10 ignored to
**175 passed / 0 failed / 0 ignored** and the Python suite from 382 passed / 1 failed to
**382 passed / 0 failed / 0 skipped**. Full evidence and the whole measurement history is in
`odd/tasks/mitc4plusd-2025-purity.md`, iterations 24-30.

### 13.1 MITC4: `K_t(0)` equals `K_0` RELATIVELY, not exactly

The element builds `B(0)` by finite differences at step `N_GAMMA_B_H = 2e-5`. Measured on the
`4 x 4` quad plate of `tests/test_rust_assembler.py`:

```text
|K|inf = 4.153846e+09        |D|inf = 2.624178e-02
|D|inf/|K|inf = 6.317465e-12         |D|F/|K|F = 4.522245e-12
```

`K_t(0)` reproduces `K_0` to **6.3e-12 relative**: that is the finite-difference noise floor,
and it is four decades better than a naive FD would give.

The removed Python assertion paired `atol=1e-6` with `rtol=1e-10`. On entries of magnitude
`1e9` that pair demands `2.4e-16` relative, which is **below double-precision epsilon
(`2.20e-16`)**: unmeetable by construction, by an FD implementation and by an exact one
alike. The failure was therefore a defect in the assertion, not in the element. The MITC4
case was removed from `test_kt_at_zero_equals_k`; the MITC3 case stays, and the property is
still asserted for MITC4, with a relative tolerance, by the Rust gate
`n_gamma_rigid_body_zero_force_and_consistent_tangent` (which passes).

The exact route is designed and PENDING: analytic `B_L` and `N` such that `B(u) = B_L + N u`,
`K_t(0) = integral(B_L^T W B_L) = K_0` exactly, and no finite differences anywhere in the
tangent. The declaration-by-evaluation trick that makes it exact, and the two reverted
wiring attempts, are in the task doc.

### 13.2 The analytic B/N route is not in the tree

The two wiring attempts of iterations 27-29 were reverted, so the tree keeps the validated
finite-difference tangent. What those attempts established, measured, is that the printed
Eq. (20c)'s BENDING-ROTATION coupling comes out **uniformly 4x** the validated operator on
every bending-rotation pair, and that this constant is what a square fixture names (a `2 x 1`
rectangle had transformed it into the misleading `rr` 1 / `rs` 2 / `ss` 4 pattern). The
acceptance test for the fix is in the tree and passes today: the block gate
`n_gamma_b_linear_block_check` asserts `B(0)` block by block on three fixtures and would turn
red the moment that wiring lands with the factor still wrong.

### 13.3 Instruments retired, and where their findings live

Four `#[ignore]`d instruments were deleted, none of them because it was red:

| retired instrument | why |
| --- | --- |
| `wu9i_block_isolation_flat_inplane_strip` | superseded: its defect was fixed and folded into the element, and the gate covers the family |
| `n_gamma_kt_residual_diagnostic` | ZERO assertions: a vacuous gate. The property it measured is limitation 13.1 |
| `n_gamma_geo_symmetry_diagnostic` | ZERO assertions; its finding is folded into `N_GAMMA_GEO_H = 1e-4` and the comment at `mitc4_plusd.rs:2842` |
| `n_gamma_geo_stencil_probe` | ZERO assertions; same finding, same recorded place |

The lesson worth keeping: an `#[ignore]`d test with no assertion hides a vacuous gate. The
three that had zero assertions had been "passing" for the whole change without checking
anything, and that is exactly what the ignore attribute was concealing.

### 13.4 Runtime cost accepted deliberately

`n_gamma_rigid_body_zero_force_and_consistent_tangent` costs **107 s**, so
`cargo test -p aeroelast-core` is now about 125 s instead of 4 s. It was promoted anyway
because it is the only active guard of `K_t == dF/du` and `F_int(0) == 0` once the Python
MITC4 case is gone. Re-ignoring it is a one-line change if the iteration cost ever outweighs
that guard.

## 12. Provenance

Read for this design: the proposal (634 lines, rev. 2), the spec (504 lines, 16
requirements / 39 scenarios), `docs/formulations/mitc4plus-2017-extract.md`,
`docs/formulations/mitc4plusd-2025-extract.md`, `docs/formulations/shell-elements.md`
(§1–§2, §4.2), `crates/aeroelast-core/src/elements/mitc4.rs` (struct, constructor,
`b_m_mitc4_plus`, `b_md_mitc4_plus`, `drill_midside_shape_derivatives`,
`b_gamma_mitc4`, `compute_ke_local_erc`, `compute_ke_local`, `compute_fint_global`,
`compute_bending_shear_condensed`, `transform_to_global`, the test module),
`crates/aeroelast-core/src/materials/{mod,laminate,isotropic,composite}.rs`,
`crates/aeroelast-core/src/assembly/{assembler,topology}.rs`,
`crates/aeroelast-py/src/{elements,assembler}.rs`, `src/aeroelast/core/laminate.py`,
`openspec/config.yaml`, and the Engram project context (`sdd-init/fem-shell`,
project `fem-shell`).

Paper B (2025) itself was **not** re-read in this phase; §1.1/§2.2's paper-B citations
come from `docs/formulations/mitc4plusd-2025-extract.md` (a vision read by the
requester) and from the dead code's doc comments. Per §9-1, the vision re-read is WU0
and it precedes the drill operator's implementation.

**Rev. 2 (corrective) re-resolution.** Every symbol and file:line this revision cites was
re-read against the tree rather than recalled. The specific verifications:

| Claim | Resolved at |
| --- | --- |
| `MeshTopology` has exactly six fields; no `node_normals` anywhere | `crates/aeroelast-core/src/assembly/topology.rs:46-59`; `grep -rn "node_normals" crates/` → no match |
| `Mitc4Precomputed::initial_normals` is the only per-node normal, filled with one `e3` | `crates/aeroelast-core/src/elements/mitc4.rs:105`, `:582-584`, used at `:3668` |
| The element constructor receives only 12 coords + constitutive (no connectivity) | `crates/aeroelast-core/src/elements/mitc4.rs:497-505` |
| Multi-ply `cs` is not `k·G·h`; `k` is applied only in the single-ply branch and the degenerate fallbacks | `crates/aeroelast-core/src/materials/laminate.rs:137-138`, `:176`, `:178-187`, `:196` |
| `ShellConstitutive` has exactly five fields and its `cs` doc comment says `[k·G·h] form` | `crates/aeroelast-core/src/materials/mod.rs:20-31` (comment at `:28`) |
| `batch_ke_mitc4_composite` has no factor channel | `crates/aeroelast-py/src/elements.rs:532-540` |
| Eq. (27c)'s `(1 + a_E·rs) e_rs^m(E)` proves the missing leading term in the printed Eq. (21) | `docs/formulations/mitc4plus-2017-extract.md:388-390` vs `:347-349`, `:359-361` |
| Paper B's drill uses `θ^D = θ·V^D` and the `θ_{i-1}^D − θ_i^D` edge difference | `docs/formulations/mitc4plusd-2025-extract.md:67-70`, `:76-80` |
| The extract's three unchecked parts are Eq. (8a), Eq. (8b) and the `∂x_b` terms | `docs/formulations/mitc4plus-2017-extract.md:491-499` |
| `MeshTopology::new` is called with three arguments from the PyO3 assembler (so an internally-filled field would not move PyO3) | `crates/aeroelast-py/src/assembler.rs:92`, `:146`, `:405`; `aeroelast-core/src/assembly/assembler.rs:1120`; `topology.rs:165` |
| `update_reference` mutates `topology.node_coords` in place before rebuilding elements (why a cached `node_normals` would need maintenance) | `crates/aeroelast-core/src/assembly/assembler.rs:217-230` (mutation at `:227`) |

No symbol in this document was left unresolved: the one that did not exist in Rev. 1
(`MeshTopology.node_normals`) is removed, and its replacement (ADR-4 option B) is built
from symbols that do exist.
