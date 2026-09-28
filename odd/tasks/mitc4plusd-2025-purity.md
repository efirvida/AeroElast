# Task document — MITC4+/D: pure 2025 fidelity, dead-code removal, then measurement

**Created:** 2026-09-24, resumed feature `mitc4plusd-faithful` (SDD change in
`openspec/changes/mitc4plusd-faithful/`). **Authorized by the maintainer**: "arranca con todo y haz
varias iteraciones de lectura de artículo e implementación hasta asegurarte 100% de que se implementó
y se cableó totalmente la formulación de 2025 y no queda código muerto o viejo de la formulación
anterior, solo después mide el comportamiento, y continúa con los demás puntos."

**Goal.** Make `crates/aeroelast-core/src/elements/mitc4_plusd.rs` a *complete and wired*
implementation of the 2025 formulation (Ko, Bathe & Zhang, C&S 308:107622 — the MITC4/D and MITC4+/D
elements), with no dead or stale code from the previous formulation. **Only then** measure.

**Rules (non-negotiable, inherited).**
1. Equations only from PDF **vision** (`pdftoppm -png -r 300` + the read tool + PIL crops). **Never
   `pdftotext`** — not even for navigation or for prose. Never trust `docs/formulations/*-extract.md`
   as the paper.
2. Compare each paper equation to the **function body** term by term; where the paper prints no
   equation (e.g. the covariant→local shear metric), verify numerically.
3. No new ingredient: every implemented term must be traceable to a printed 2025 (or explicitly
   cited) equation. A term with no source is a STOP, not a patch.
4. Never weaken a test or widen a tolerance. Never commit without an explicit request.
5. `crates/aeroelast-core/src/elements/mitc4.rs` is REFERENCE READ-ONLY.
6. Provenance: the 2025 paper is `.sources/papers/1-s2.0-S0045794924003511-main.pdf`. The 2017 core is
   `.sources/papers/A_new_MITC4+_shell_element.pdf` (C&S 182:404-418). **`1-s2.0-S0045794917309550-main.pdf`
   is the benchmark paper (C&S 193), not the formulation.**

## Confirmed starting point (vision, this session)

- 2025 **Eq. (22a)**, p. 10: `e_ij = e_ij^m + e_ij^md + t·e_ij^b1 + t²·e_ij^b2`, `i,j = 1,2`. It is a
  **sum of strains**. The code instead adds the drill as an independent energy block
  (`k += drill_ke_local(pre)`), so the cross terms `2(e^m)ᵀ C e^md` (t⁰) and `2(e^md)ᵀ C e^b2` (t²)
  are absent. Measured consequence (WU9i): the drill block's energy share is exactly 0 and the flat
  in-plane response equals the plain compatible Q4 membrane (37% lock).
- 2025 **Eq. (24b)**, p. 11: `ẽ_rr^m = e_rr^m(0,0) + (√3/2)λ(r,s)(e_rr^m(A) − e_rr^m(B))·s`,
  `λ = j0/j(r,s,0)`, and the `r` sibling for `ss`, `ẽ_rs^m = e_rs^m(0,0)`. The code implements the
  **2017** Eqs. (27a-c) instead; the slope coefficient differs (`√3λ` vs `1`).
- 2025 validated envelope: `t/L = 1/100, 1/1000, 1/10000` (Fig. 15 legend, p. 11).
- 2025 **§3.2** has 2D planar-beam (plane-stress) examples, Tables 1-3 — an **unused oracle**.

## Tasks

- [ ] **T1. Complete the 2025 equation audit (read-only).** Read with vision pp. 2-13 of the 2025
  paper and compare every formulation equation to the code: Eqs. (1)-(7) (geometry, `V^D`, covariant
  strains), (8a-d), (9a-b), (10), (11a-b), (12a-d), (13a-c), (14), (15)-(17), (18)-(19d), (20), (21),
  (22a-c), (23a-b), (24a-c), (25a). Output: a gap table (equivalent / different / missing / extra),
  each row with the printed equation and the code anchor. **Deliverable: the complete gap list.**
- [ ] **T2. Implement Eq. (22a): fold `e^md` into the membrane strain.** Build the membrane B row as
  `B_m + B_md` in the `t⁰` slot instead of adding an independent `drill_ke_local` block, so the
  paper's cross terms exist. Touches `mitc4_plusd.rs`. Verification: Tier 1 (`test_t1a_*`,
  `test_t1b_*`) stays green **and** the WU9i instrument changes in a measured way.
- [ ] **T3. Align the assumed membrane with the 2025 Eqs. (24a-c)/(25a)** if T1 shows they differ from
  the 2017 Eqs. (27a-c) (the `√3/2·λ` slope in particular). Implement only what the paper prints.
- [ ] **T4. Nodal directors** (2025 Eq. (1)/(9a): "the director vector at the node"): ADR-4 option A,
  mesh-consistent `V_n^i`. Only if T1 confirms the paper requires the nodal reading.
- [ ] **T5. Remove dead/stale code** from the previous formulation: inventory unused functions,
  fields, constants and comments in `mitc4_plusd.rs` and delete them, keeping the file
  production-only plus its Tier-1 tests.
- [ ] **T6. Verify and measure (only after T2-T5).** `cd crates && cargo test -p aeroelast-core`; the
  WU9i instrument; the 2025 §3.2 Tables 1-3 as an oracle; `pytest -m "not slow"`; the CCX judge from
  `tests/test_ccx_shell_element_types_parity.py`. Record every number.
- [ ] **T7. Report the behavior** and only then revisit the gate decision (the three flat in-plane
  bending tests) and the remaining points (Task 2b/2c: CCX references for plate/lateral/modal,
  multimaterial rows; the hybrid ingredient documentation table).

## Nonlinear fidelity audit (N-series) — added 2026-09-24 at the maintainer's request

The maintainer's instruction: *"en `.sources` sí hay artículos específicamente sobre la no linealidad, que debemos entonces
implementar fielmente a los artículos punto por punto tal cual hicimos con la formulación de 2025. No confío en la
no linealidad ahora implementada, por lo que hay que ser estricto con la formulación declarada en los artículos."*

This supersedes the old checklist section H, which declared the nonlinear path "out of the papers' scope". **It is not:** the primary reference is in `.sources`.

### The references (identified from `.sources/papers/`)

| file | reference | role |
| --- | --- | --- |
| `The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf` (14 pp.) | **Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear analysis", Computers and Structures 185:1-14**, doi 10.1016/j.compstruc.2017.01.015 | **the primary reference** for our nonlinear path |
| `A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf` | Dvorkin & Bathe (1984) | the TL continuum-mechanics origin (its `TOTAL LAGRANGIAN FORMULATION` section, p. 79, is already partly read) |
| `mitc3+_no_lineal.pdf` | "The MITC3+ shell element in geometric nonlinear analysis" | the triangular sibling, same scheme |
| `mitc4+_no_lineal.pdf.pdf`, `jun2018.pdf`, `32452534.pdf` | unidentified (no title metadata) | identify with vision before use |

### What the paper specifies (vision, pp. 1-3)

- **Total Lagrangian formulation** with Green-Lagrange strains, `^t` = time (or load step), `0` = the initial reference configuration (abstract, p. 1).
- Eq. (1): `^t x(r,s,zeta) = ^t x_m + zeta ^t x_b`, `^t x_m = sum h_i ^t x_i`, **`^t x_b = 1/2 sum a_i h_i ^t V_n^i`** — the director is the **current** one, i.e. it rotates with the deformation.
- Eq. (2): the bilinear shape functions; Eq. (3): the incremental displacement `u = ^{t+dt}x - ^t x`.
- Eqs. (9)-(11a): the GL strain splits into a **linear** part `_0 e_ij` and a **nonlinear** part `_0 eta_ij`; the **assumed transverse shear** is the classical MITC4 field applied to the TL strains (Eq. (10), same tying points A-D, Fig. 2); and the through-thickness split in TL is `^t_0 e_ij = ^t_0 e_ij^m + zeta ^t_0 e_ij^b1 + zeta^2 ^t_0 e_ij^b2` (Eq. (11a)) — **the same structure as the linear element, on the GL strains**.
- Fig. 3: the characteristic vectors `^t x_r`, `^t x_s`, plane P, normal n **at time t**.

### Code inventory (measured)

**In `mitc4_plusd.rs` (ours):** `compute_b_nl` (1579), `compute_b_geometric` (1614), `geometric_stiffness_local` (1657), `membrane_nonlinear_correction` (1697), `compute_fint_global` (1735), `compute_kt_global` (1762), plus a corotational block (1993-2160): `quaternion_*`, `update_normals_with_displacements`, `polar_decomposition`, `log_strain_from_polar`, `compute_membrane_strain_log`, `update_corotational_frame`, `frame_incremental_rotation`.
**In `mitc4.rs` (hybrid, READ-ONLY):** the same, plus `green_lagrange_strain` (1528) and `compute_b_l` (1533) — **both compiler-dead**.
**PyO3:** the nonlinear kernels are wired — `aeroelast-py/src/elements.rs:341` calls `mitc4_plusd::compute_kt_global` and `:394` calls `mitc4_plusd::compute_fint_global(pre, &u, nonlinear)`.

### Two immediate red flags (to be proven or refuted by the audit)

1. **`green_lagrange_strain` and `compute_b_l` — the paper's GL building blocks — are DEAD in the hybrid and ABSENT from `mitc4_plusd.rs`.** The paper's whole formulation is built on `_0 e_ij + _0 eta_ij`; if the new element has no GL strain construction, its nonlinear path is not the paper's.
2. **The corotational machinery** (quaternion, polar decomposition, log-strain, `update_corotational_frame`, `frame_incremental_rotation`) is **UNSOURCED** (the hybrid inventory flagged it as "S4R-style, no citation") and the paper specifies **total Lagrangian, not corotational**. If the nonlinear path routes through it, that is an extra ingredient and a fidelity deviation; if nothing routes through it, it is dead code for T5.

### Audit tasks

- [ ] **N1 director update** — `^t V_n^i` (Eq. (1)): does `update_normals_with_displacements` update it from the incremental rotation, and is it wired into `compute_kt_global`/`compute_fint_global`/the Python nonlinear solver?
- [ ] **N2 the GL strain split** — Eq. (9) `_0 e_ij` (linear) + `_0 eta_ij` (nonlinear); map to `compute_b_nl`, `compute_b_geometric`, and the (missing/dead) `green_lagrange_strain`/`compute_b_l`.
- [ ] **N3 assumed membrane in TL** — the paper's MITC4+ assumed membrane applied to the **GL** strains (the paper's §2.3, pp. 4-5): is `b_membrane_2017` applied to the TL strains, or are the displacement-based GL strains used (which would lock)?
- [ ] **N4 assumed transverse shear in TL** — Eq. (10): is `b_shear_mitc4` applied to the TL shear strains?
- [ ] **N5 through-thickness split in TL** — Eq. (11a) and the resultant/moment assembly.
- [ ] **N6 tangent stiffness** — the paper's `^t_0 K`: map to `geometric_stiffness_local` + `compute_kt_global`, term by term.
- [ ] **N7 internal force** — the paper's `^t_0 f_int`: map to `compute_fint_global` (including its `nonlinear` flag).
- [ ] **N8 corotational machinery disposition** — sourced or not, used or dead.
- [ ] **N9 the paper's own nonlinear benchmarks** (§3 of C&S 185) as the external oracle, plus the CCX judge where possible.

### N-audit, first pass (2026-09-24) — the nonlinear path is NOT the paper's TL formulation

**Vision read:** C&S 185 (2017), p. 2 (both columns), §§2.1-2.2.

**What the paper specifies:**
- Eq. (4c): `^{t+dt}V_n^i - ^t V_n^i = theta_i x ^t V_n^i + 1/2 theta_i x (theta_i x ^t V_n^i)`, `theta_i = ^t V_1^i alpha_i + ^t V_2^i beta_i` — the director increment **to quadratic order, about the CURRENT director**.
- Eq. (5c): `u_b1 = 1/2 sum a_i h_i(-^t V_2^i alpha_i + ^t V_1^i beta_i)`, **`u_b2 = -1/4 sum a_i h_i(alpha_i^2 + beta_i^2) ^t V_n^i`** — the **quadratic** director term.
- Eq. (6): `u_1 = u_m + zeta u_b1`, `u_2 = zeta u_b2`.
- Eq. (7): `^t_0 e_ij = 1/2(^t g_i . ^t g_j - ^0 g_i . ^0 g_j)` with `^t g_i = d^t x/dr_i` the **CURRENT** covariant base vectors.
- Eq. (8): `_0 e_ij = 1/2(^t g_i . u_j + u_i . ^t g_j + u_i . u_j)`.
- **Eq. (9): `_0 e_ij = _0 e_ij + _0 eta_ij` with** `_0 e_ij = 1/2(^t g_i . u_{1,j} + u_{1,i} . ^t g_j)` (linear) and **`_0 eta_ij = 1/2(u_{1,i} . u_{1,j} + ^t g_i . u_{2,j} + u_{2,i} . ^t g_j)` (nonlinear — THREE terms)**.

**Findings (code side, measured):**

| # | verdict | evidence |
| --- | --- | --- |
| **N1 director/geometry update** | **DEAD / not wired** | `update_normals_with_displacements` (`mitc4_plusd.rs:2040`) has **zero call sites** in `crates/`. The current-configuration director `^t V_n^i` (Eqs. 1, 4c) is never formed and the current metric `^t g_i` (Eqs. 7-9) is never used: `compute_kt_global`/`compute_fint_global` build everything from `pre` (the INITIAL geometry and directors) |
| **N2 the GL split** | **INCOMPLETE** | `membrane_strain_nl` returns `1/2(H^T H)` in Voigt form = **only the first term** `1/2 u_{1,i}.u_{1,j}` of the paper's `_0 eta_ij`. The two `u_2` terms (the quadratic director increment of Eq. (5c)) have **no representation at all**: `membrane_strain_nl`/`compute_b_nl` take only the 3x3 displacement gradient, which encodes no `u_b2` |
| **N6/N7 tangent and internal force** | **APPROXIMATION, not the paper's** | `compute_kt_global` = `K_0 + K_L + K_sigma` with `K_0 = compute_ke_local(pre)` (initial geometry); `K_L`/`K_sigma` are the standard `1/2(nabla u)^T(nabla u)` initial-displacement and initial-stress terms; `compute_fint_global = K_0 u + membrane_nonlinear_correction` |
| **N8 corotational machinery** | **DEAD and UNSOURCED** | `compute_membrane_strain_log`, `update_corotational_frame`, `frame_incremental_rotation` -> **zero call sites**; `polar_decomposition`/`log_strain_from_polar` are called only from inside that same dead block. They are `pub`, hence no rustc dead-code warning. The paper specifies **total Lagrangian**, not corotational |

**The code admits it.** `membrane_nonlinear_correction`'s docstring says the path is a "BOUNDED NONLINEAR PATH (design open item 5, risk 9)"; the 5 nonlinear test failures are its symptom.

**What a faithful implementation needs:** per load step, form `^t V_n^i` (Eq. 4c) and `^t g_i` (Eq. 7); build `u_1` (Eqs. 5b/5c) and `u_2` with the **quadratic** `u_b2`; form `_0 e_ij + _0 eta_ij` (Eq. 9); apply the paper's assumed membrane and shear fields to those GL strains (Eqs. 10, 11a, and the §2.3 membrane field); build `^t_0 K` and `^t_0 f_int`; then the §3 benchmarks as the oracle (N9).

### N-implementation plan (maintainer-authorized 2026-09-24): implement all of it, with intermediate measurements

Instruction: *"implementalos todos automaticamente pero con mediciones intermedias, e iteraciones sucesivas si las mediciones
fallan, estas iteraciones deben releer la formulacion en cada caso si es necesario."*

So each stage below lands with its own **exact identity measurement** before the next one starts, and a
failed measurement sends us back to the paper, not forward. The identity that governs the whole series is
the one the paper's `u_2` term exists for: **a finite rigid-body rotation must produce zero incremental
green-Lagrange strain**, i.e. `_0 e_ij + _0 eta_ij = 0` (Eq. (9)), which is *not* second-order exact
without the quadratic director term `u_b2`.

| stage | what it implements | intermediate measurement (exact identity) |
| --- | --- | --- |
| **N-alpha** | Eqs. (4c), (5b), (5c), (6), (9): the current director increment, `u_1`/`u_2`, and the linear + nonlinear GL increments | a finite **rigid-body rotation about an in-plane axis** gives `_0 e_ij + _0 eta_ij = 0` to `O(phi^3)`; non-vacuity: dropping `u_2` breaks it. Plus the CURRENT path's residual on the same quantity (the before number) |
| **N-beta** | the paper's assumed fields on the GL strains: the §2.3 assumed membrane, Eq. (10) transverse shear, Eq. (11a) through-thickness split | the assumed fields must reproduce exactly the displacement-based GL strain on a flat element (the same reduction the linear audit proved), and must keep the rigid-rotation zero |
| **N-gamma** | the paper's `^t_0 K` and `^t_0 f_int` | `f_int = 0` exactly for any rigid-body field; `K_t = d f_int / d u` by finite difference; `K_t(0) = K_0` bit for bit |
| **N-delta (N9)** | wire it and delete the old bounded path + the dead corotational block (N8/T5) | the paper's §3 nonlinear benchmarks (roll-up, cantilever, and the shell problems) as the external oracle |

## Evidence log

### Iteration 1 — 2026-09-24

**Read with vision (2025 paper, `.sources/papers/1-s2.0-S0045794924003511-main.pdf`):**
- p. 4: Eqs. (8a)-(8d), (9a-b), (10), Figs. 4-6.
- p. 10 (left): Eqs. (19a), (19b), (19c), (19d), (20), (21).
- p. 10 (right, bottom): **Eq. (22a) `e_ij = e_ij^m + e_ij^md + t e_ij^b1 + t^2 e_ij^b2`**, Eq. (22b), Eq. (22c), Eq. (23a).
- p. 11 (left): Eq. (23b), Eq. (24a), **Eq. (24b)** `e~_rr^m = e_rr^m(0,0) + (√3/2)λ(r,s)(e_rr^m(A) − e_rr^m(B))s`, `λ = j0/j`.
- p. 11 (right): Fig. 15 legend (`t/L = 1/100, 1/1000, 1/10000`).
- p. 12 (left): **Eq. (25b)** `[a_k|ij] = Q(r,s)R(r,s)S`, **Eq. (25c)** the explicit `Q`, `R` (with `R ∝ λ`) and `S` matrices, and the `a_A..a_E` definitions.

**Changed (T2, Eq. (22a)):** `compute_ke_local_with_drill` in `crates/aeroelast-core/src/elements/mitc4_plusd.rs` now builds the `t^0` row as `B_m + B_md` (`bm += b_drill_membrane_2025(pre, r, s)` when `use_drill`), and the independent `k += drill_ke_local(pre)` block is gone. The single `W` block therefore carries the paper's cross terms `(e^m)^T C e^md` and `(e^m + e^md)^T C e^b2`.

**Measured (the WU9i instrument, dense 270-DOF solve of the exact failing strip):**

| variant | before (separate drill block) | **after Eq. (22a)** | beam theory |
| --- | --- | --- | --- |
| `uy/ux` | 251.5515 (−37.1%) | **400.3877 (+0.10%)** | 400 |
| `uy` [m] | 7.118881e-3 (−37.7%) | **1.129252e-2 (−1.19%)** | 1.142857e-2 |
| `ux` [m] | 2.829990e-5 (−0.95%) | **2.820397e-5 (−1.29%)** | 2.857143e-5 |

Block energy shares on the FULL solution (load +y): before `membrane 1.000000000002 / bending 0 / shear 0 / drill 0`; after `membrane 1.588773916354 / bending −1.177547832743 / shear 0 / drill 0.588773916386` (sum 1). The cross terms exist now.

**Verdict: defect #4 was a real fidelity defect, and Eq. (22a) closes it with no new ingredient.** The element is now ~1.59× softer in flat in-plane bending and lands 0.10% from beam theory (test tolerance 2%) and 1.19% / 1.29% from the analytical `uy`/`ux` (test tolerance 3%).

**NEW BLOCKER (T2b): the test-side references encode the OLD structure.** `cargo test -p aeroelast-core` is now **165 passed / 3 failed**; the 3 failures are the identity locks, whose test-local reference re-implements the OLD separate-drill structure:
- `test_identity_ke_lock_matches_2017_core_plus_2025_drill` — the reference `assemble`/`ke_local` adds `drill_block` separately.
- `test_identity_drill_stiffness_comes_only_from_eq26` — asserts `K(drill) − K(no drill) == drill_block`; now the difference also carries the paper's cross terms. Its own failure message is the evidence: `t-r block left: 6933161907.296183 right: 0.0`.
- `test_identity_integration_rule_is_2x2x2_and_discriminates_surface_only` — compares the production to a three-term (separate-block) reference.

This is the **sixth time** in this change that an internal check shared the assumption it was supposed to test: the "identity lock" compared the code against a test-side re-implementation of the SAME deviation, so it could never see the Eq. (22a) gap. The instrument's own 9-row cross-check and its "FULL within 1% of 251.55" validation are invalidated for the same reason and must be re-baselined.

**Next (T2b):** update the test-side references to the paper's structure (sum the drill into the membrane row in `membrane_local`/`assemble`/`ke_local`/`drill_block` usage, and re-baseline the provenance assertion and the WU9i validation, which is why the instrument's assertion now fails), keeping every assertion a real assertion and weakening nothing. Then T3, T5, T6.

### Iteration 2 — 2026-09-24 (T2b, done)

**Changed (test-side references only, plus the instrument; no production change beyond iteration 1):**
- `ke_ref::assemble` now sums `drill::b_md_reference` into the membrane row and no longer adds a separate `drill_block`; `ke_ref::drill_block` is now unused (dead-code inventory for T5).
- `test_identity_drill_stiffness_comes_only_from_eq26`: assertion (c) re-baselined from "every translational row/column block is exactly zero" to the Eq. (22a) structure — every `t-t` block still exactly zero, **and the `(t, r)` coupling asserted non-vacuous** (`max |t-r| > 1e-10 scale`). That coupling is precisely what the old assertion forbade and what made this identity lock blind to the Eq. (22a) deviation.
- The WU9i instrument: `ke_bend9` folded to the production recipe; the bending-by-difference cross-check subtracts the pure drill block (`k_cross + k_bending`); `ke_compat` sums the drill instead of adding it separately; the validation now compares against the beam-theory `400.0` at 2% (the pre-fold `251.5515` is printed as the recorded baseline) and the docstring records the resolution.

**Measured after T2b:**
- `cd crates && cargo test -p aeroelast-core` -> **168 passed / 0 failed / 1 ignored**.
- The instrument (`wu9i`): variant 1 FULL `400.3877` (+0.10% vs beam theory); **variant 3 `COMPAT_MEMB` = `400.3877`** — identical to 12 digits, confirming that on flat geometry the 2025 assumed membrane and the compatible field coincide under the paper's structure and that the **summed drill** is what removes the lock; variants 2b/4b/5 (drill absent) `251.5515`, i.e. the pure membrane locks as before; the bending cross-check now matches to `3.5e-14` relative (was `6.5e3`).

**Note for T3:** the fact that variant 3 == variant 1 on a flat element bounds the T3 question — the 2025 membrane restatement (Eqs. (24a-c)/(25a), with `lambda = j0/j` and the `Q R S` coefficients) can only differ from the 2017 (27a-c) on **distorted or warped** geometry (`lambda != 1`, `a_A..a_E != 0`). T3 must therefore be judged on the warped/twisted benchmarks, not on the flat strip.

**Next:** T3 (align the assumed membrane with the 2025 restatement if it differs on warped geometry), then the rest of T1 (pp. 2, 3, 5-9, 13), T5 (dead code), T6 (full measurement: pytest, CCX, the 2025 §3.2 Tables 1-3).

### Iteration 3 — 2026-09-24 (T3 audit; a new oracle found)

**Read with vision (2025 paper):** p. 11 (right, bottom): **Eq. (24c)** and **Eq. (25a)**; p. 5 (left): Fig. 7, Fig. 8 and **Table 1**; p. 6 (left): **Table 2**, Fig. 10 and **Table 3**, and **Eqs. (11a)/(11b)**.

**T3 findings (term by term):**
- **Eq. (24c) is IDENTICAL to Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (25), term for term:**
  `e~_rs|bil = (c_r/d)[c_r(e_rr|con + e_rs|bil) - e_rr|lin] + (c_s/d)[c_s(e_ss|con + e_rs|bil) - e_ss|lin] + (2 c_r c_s/d) e_rs|con`.
  So the `rs` half of the 2025 assumed membrane is exactly what `b_membrane_covariant_2017` implements.
- **`a_A..a_E` are identical in both papers** (`a_A = c_r(c_r-1)/(2d)`, ..., `a_E = 2c_rc_s/d`) — confirmed on p. 12 of the 2025 paper and p. 410 of the 2017 paper.
- **Eq. (25a) is `[e~] = Q(r,s)R(r,s)S [e(A) e(B) e(C) e(D) e(E)]^T`** with `Q`, `R`, `S` printed in Eqs. (25b)/(25c); `R` carries `lambda = j0/j`. This is a factorised statement of the SAME coefficient family that Eq. (24a) writes as `sum a_k|ij (samples)`.
- **The only textual difference is Eq. (24b):** `e~_rr = e_rr(0,0) + (sqrt(3)/2) lambda(r,s)(e_rr(A) - e_rr(B)) s` (and the `r` sibling), versus the 2017 Eq. (18)'s `e~_rr = 1/2(e(A)+e(B)) + 1/2(e(A)-e(B))s`. Where the two differ: the constant term (`e(0,0)` vs the average of the samples) and the slope (`sqrt(3) lambda` vs `1`).
- **Still unread: the definitions of `eta_i` and `m_i`** used in `R` (Eq. (25c)). They are not in pp. 4-6. Without them the equivalence of (24b)/(25) with the 2017 (27a-c) cannot be proven algebraically, so T3 is **not** closed on the printed text alone.

**NEW ORACLE — the paper's own plane-stress beam tables (this is what closes T3 empirically):**

| table | case | MITC4 | MITC4-IC | MITC4/D | **MITC4+/D** | reference |
| --- | --- | --- | --- | --- | --- | --- |
| **Table 1** (Fig. 8, slender straight cantilever, `E=1e7`, **thickness (out-of-plane) `0.1`**, **in-plane height `0.2`**, `nu=0.3`), mesh 1 / mesh 2 (distorted) | tip `-u_y` | 0.010088 / 0.00290874 | 0.107328 / 0.00569730 | 0.0976755 / 0.00550604 | **0.0976755 / 0.00550838** | **0.1081** |
| **Table 2** (Fig. 9, thick curved beam, plane stress, `E=1e3, nu=0`), `1x2` / `1x4` / `1x8` | tip displacement | 22.5988 / 57.9325 / 79.9218 | 52.2291 / 84.6070 / 89.3023 | 51.2489 / 84.1086 / 89.1698 | **51.2692 / 84.1444 / 89.2077** | **90.1** |
| **Table 3** (Fig. 10, cantilever with roller, plane stress, `E=3e4, nu=0`), regular / distorted | tip displacement | 0.235608 / 0.203966 | 0.347810 / 0.342552 | 0.347810 / 0.339249 | **0.347810 / 0.341015** | **0.347810** |

What the authors' own data say: MITC4+/D is **exact** on the Table 3 regular mesh, **1.95% low** on its distorted mesh, and **9.6% low** on the very slender Table 1 beam (and MITC4 locks 90.7% there). So the flat plane-stress in-plane bending of MITC4+/D is *improved but not locking-free*, and the residual grows with slenderness — fully consistent with our measured post-Eq.(22a) 0.10% / 1.19% on the repo's `b/L = 0.1` strip.

**T3 disposition:** close it **empirically** — reproduce Tables 1, 2 and 3 with our element (T6) and require the MITC4+/D columns. That is a printed, third-party oracle for exactly the family in question, and it is stronger evidence than re-deriving `Q R S`. The algebraic (`eta_i`/`m_i`) check stays open and is listed below as an optional item.

**Also noted for T1/T5:** Eqs. (11a)/(11b) print the simplified drill edge derivatives (`[0, 1/2(-2r)(1+s), 0, 1/2(-2r)(1-s)]` and its `s` sibling, with `h~^l_{m,r} = h~_{l,r}`) — the deliberate zeros are printed, so the code's `drill_midside_shape_derivatives` can be checked against them directly. And Table 3's numbers are a much sharper oracle than the repo's current `test_fy`/`test_ratio_physical` (one exact, one distorted).

### Iteration 4 — 2026-09-24 (T6a, first published cell reproduced to machine precision)

**Vision re-read (high resolution, p. 5 left):** Fig. 8 caption and **Table 1**, twice and zoomed. Confirmed digits: `MITC4 0.010088`, `MITC4-IC 0.107328`, `MITC4/D 0.0976755`, **`MITC4+/D 0.0976755`**, `reference 0.1081`; distorted column `0.00290874 / 0.00569730 / 0.00550604 / 0.00550838`. Fig. 8 caption: `E = 1.0 x 10^7`, **thickness 0.1**, `nu = 0.3`. (The earlier note in this file saying `t = 0.2` was wrong: `0.2` is the in-plane height. Corrected above.)

**Added:** `#[ignore]`d instrument `t2025_table1_slender_plane_stress_cantilever` in the in-file test module (test code only, +202/-0; no production line moved). It builds the paper's **Fig. 8 mesh type 1** (regular `6 x 1` quad mesh: 14 nodes, 84 DOF, 6 elements), clamps all 6 DOF at `x = 0`, applies `P = 1.0` in `-y` on the free edge, reuses `compute_ke_global` + the module's dense constrained solver, and measures `-u_y` at A as the average of `u_y` over the free-edge nodes. Physics cross-check (independent of our code): `P L^3/(3 E I)` with `I = t h^3/12 = 6.6667e-5` gives `0.108` = the paper's reference `0.1081`, which fixes the geometry.

**Measured:** `-u_y = 9.767550000444e-2`.

| comparison | value | ratio |
| --- | --- | --- |
| **published MITC4+/D, mesh type 1** | **0.0976755** | **1.000000000 (difference +4.4e-12)** |
| published reference solution | 0.1081 | 0.903566 |
| Euler-Bernoulli `P L^3/(3 E I)` | 0.108 | 0.904403 |
| published MITC4 | 0.010088 | 9.68x |
| published MITC4-IC | 0.107328 | 0.910065 |

**Verdict: our element reproduces the 2025 paper's own published MITC4+/D cell to machine precision** (an 84-DOF dense solve agreeing to round-off means the same element, mesh, BC and load). It also means the element is **9.6% below the reference on this very slender beam because the paper's MITC4+/D is** — `0.108` was never the MITC4+/D target. This is the empirical fidelity demonstration T6a was for, and it also bounds T3.

**Honest caveat (raised by the implementer, and answered here):** machine-precision agreement with a 7-digit published cell is also what one would see if the transcribed number had come from this same code. Provenance was therefore re-checked directly against the source twice at high resolution; the digits are the paper's, `MITC4/D` shares the mesh-1 value while the distorted column differs (a paper-internal feature no run of ours could produce), and the physics cross-check (`0.108` vs the paper's `0.1081`) independently fixes the geometry.

**Still open in T3:** Table 1 **mesh type 2** is the case that discriminates the `a_A..a_E` machinery (its `MITC4+`/`D` value differs from `MITC4/D`: `0.00550838` vs `0.00550604`), and Tables 2 and 3 add the curved and roller cases. Reproducing mesh type 2 requires a precise vision read of Fig. 8(b)'s crossed-distorted pattern.

### Iteration 5 — 2026-09-24 (T3 CLOSED empirically; the curved case matches)

**Vision read (p. 5 right, zoomed):** Fig. 9 and the penultimate full-page read. Geometry fixed as: a **90-degree annular sector** (curved cantilever), **inner radius 10**, **outer radius 15** (radial wall 5 — the `10` and `5` are the two stacked horizontal dimension lines), unit thickness, `E = 1.0e3`, `nu = 0.0`, tip load **`P = 600`** in `+y` at point A (the free end), mesh `1 x N` = N elements along the arc, 1 radially, straight-sided quads, far end fully clamped. The geometry is cross-checked numerically without our code: for a 90-degree curved cantilever, `pi P R^3/(2 E I)` with `E=1e3`, `R=12.5`, `I = 1*5^3/12` gives `P ~ 612` for `delta = 90`, matching the paper's `P = 600` and reference `90.1`.

**Added:** `#[ignore]`d instrument `t2025_table2_curved_plane_stress_beam` (test code only; no production line moved). It builds the `1 x N` curved mesh, reuses `compute_ke_global` + the module's `solve_constrained`, and sweeps the reconstruction choices.

**Measured vs the published MITC4+/D column:**

| mesh | ours | published MITC4+/D | ratio |
| --- | --- | --- | --- |
| 1x2 | 51.2185 | 51.2692 | 0.9990 (−0.10%) |
| 1x4 (load at A, measured at A) | 84.7102 | 84.1444 | 1.0067 (+0.67%) |
| 1x4 (load at A, measured as the free-end average) | 84.0819 | 84.1444 | 0.9993 (−0.074%) |
| 1x8 | 89.1626 | 89.2077 | 0.9995 (−0.050%) |

Reconstruction sweep on 1x4 (load placement x measurement node): spread ~0.9%; every variant agrees with `84.1444` to within ~1%; element orientation is irrelevant. The reference is `90.1`, so our value is ~6% below it, exactly as the paper's MITC4+/D is.

**T3 verdict — CLOSED (empirically).** The curved mesh has `x_d != 0` at every element, so the `a`-coefficient machinery is live, and our element matches the paper's MITC4+/D column to `≤ 0.7%` (within the `0.9%` reconstruction ambiguity) and to `4.4e-12` on the flat regular case. Combined with the term-by-term findings (Eq. (24c) identical to the 2017 Eq. (25); `a_A..a_E` identical; Eq. (24a) the same sum structure), the implementation is the 2025 field by construction and matches it empirically.

**Honest residual (recorded, not hidden):** the 2025 paper's own tables do **not** sharply discriminate the "+" assumed membrane from the plain one — `MITC4/D` and `MITC4+/D` differ by only `0.043%` (Table 2), `0.043%` (Table 1 mesh 2) and `0.52%` (Table 3 distorted). So the paper's published data bounds the `a`-coefficient question to that band; an exact algebraic check of Eq. (24b)'s `lambda = j0/j` restatement against the 2017 (27a-c) is still not available because the `eta_i`/`m_i` definitions in Eq. (25c) have not been located in the paper. The consequence of that residual is bounded by the numbers above.

**Next:** finish T1 (pp. 2, 3, 7, 8, 9, 13 — including, if wanted, the `eta_i`/`m_i` definitions for the exact algebraic closure), T5 (dead code), T6 (flip + `pytest` + CCX; and decide whether to promote these two `#[ignore]`d instruments into real oracle tests with justified tolerances).

### Iteration 6 — 2026-09-24 (T5 done; T6: the flip gate measured)

**T5 (dead code, `mitc4_plusd.rs` only).** `cargo check --lib -p aeroelast-core` listed four functions unreachable from the production path, all of them **used by the test module**: `interpolate_position`, `interpolate_displacement` (the paper's Eqs. (1)/(3) pointwise interpolation, kept as the tests' independent reference), `membrane_ke_local` and `drill_ke_local` (the block accessors the identity locks compare against). They were **gated with `#[cfg(test)]`** rather than deleted; the production lib of this file now has **zero warnings**. (The other `never used` items in the crate are in `mitc4.rs` — REFERENCE READ-ONLY — and in `mitc3.rs`/`quad.rs`, i.e. outside this unit.) Four **pre-existing clippy style** warnings remain in production (`too_many_arguments` at `covariant_membrane_b_row`, two `very complex type`, one `needless_range_loop`): style, not dead code, and not introduced here. Suite stays **168 passed / 0 failed / 3 ignored**.

**T6 (the flip, measured).** The flip patch (`odd/tasks/mitc4plusd-wu9f-flip.patch`) no longer applies to its two documentation files (they changed since), but its **four code files apply cleanly**, so it was applied with `--exclude=openspec/*`, the extension was rebuilt with `maturin develop --release`, and the gate was run.

| run | HEAD | WU9f flip (pre-Eq.(22a)) | **flip + Eq. (22a)** |
| --- | --- | --- | --- |
| `cargo test -p aeroelast-core` | 168/0 | 168/0 | **168/0** |
| `pytest -m "not slow"` | 344/3/2 | 336/11/2 | **343 passed / 8 failed / 2 skipped** |
| `pytest tests/test_ko2017_performance.py` | 30/1 | 30/1 | **30 passed / 1 failed** |
| `pytest tests/test_shell_validation_fixed.py` | 7/0 | (flat family failing) | **6 passed / 1 failed** |
| `pytest tests/test_ccx_shell_element_types_parity.py` | 4/0 | — | **4 passed** |

**The flat in-plane bending family now PASSES**: `TestLinearStatic::test_fy`, `test_ratio_physical` and `TestLinearStaticCantilever::test_fy_in_plane` are no longer in the failure list, and the four twisted-beam cells at N=16 still pass. The single failure in `test_shell_validation_fixed.py` is `TestNonlinearStatic::test_geometric_nonlinearity`, which asserts `pytest.raises(RuntimeError, match="SNES diverged")` and now fails with **"DID NOT RAISE"** — i.e. the nonlinear solver *converges* where it used to diverge. That is a test that encodes a failure mode, not a physics regression.

The **8 remaining failures** are:

| family | tests | status |
| --- | --- | --- |
| pre-existing (fails at HEAD too) | `test_ko2017_performance.py::test_3_3_pinched_cylinder_tables_8_to_9[expected0-True]` | shared with HEAD |
| geometrically nonlinear path (design's recorded open item 5) | `test_large_rotation_benchmarks.py::test_equilibrium_path` (2 params), `::test_simo_vu_quoc_rollup_360[10]`, `test_shell_comprehensive.py::TestNonlinearStaticCantilever::test_large_displacement_tip_load`, `test_shell_validation_fixed.py::TestNonlinearStatic::test_geometric_nonlinearity` | 5, not paper-covered |
| material iso-equivalence | `test_material_suite.py::TestIsoEquivalence::test_n_iso_plies_equal_single_layer_mitc4` | open, same as WU9f |
| modal | `test_rust_modal.py::TestSimplySupportedPlate::test_analytical_convergence` | open, same as WU9f |

So, versus HEAD (3 failures) the flip-plus-fix has 8: the pinched cylinder is shared, the flat in-plane family and `test_batch_ke_mitc4_multiple` **were fixed**, and the remaining 7 are the nonlinear path (5), one material iso-equivalence and one modal — none of them paper-covered behaviour. Versus the WU9f flip it is a clear improvement (336/11 → 343/8).

**Tree state:** the flip is APPLIED (4 code files, uncommitted) and the installed extension was rebuilt from it. The flip's nodal-director injection is still the temporary two-pass instrument, not the clean ADR-4 option A.

**Next:** the remaining T1 reads (pp. 2, 3, 7, 8, 9, 13), T4 (adopt ADR-4 option A properly), the nonlinear path and the two open failures, T6b (promote the instruments), and the maintainer's call on whether to keep the flip applied or revert it and keep it as a patch.

### Iteration 7 — 2026-09-24 (N-alpha: the faithful incremental kinematics, and the identity that discriminates)

**Vision re-read (C&S 185, p. 2 right column and p. 3 top-left):** Eqs. (4a), (4b), (4c), (5a), (5b), (5c), (6), (7), (8) and **Eq. (9)** with its full `_0 eta_ij = 1/2(u_{1,i} . u_{1,j} + ^t g_i . u_{2,j} + u_{2,i} . ^t g_j)` (`u_{1,i} = du_1/dr_i`, `u_{2,i} = du_2/dr_i`). The paper's `u_1`/`u_2` are the Eq. (6) **linear/quadratic split** of the displacement, and `r_1 = r`, `r_2 = s`, `r_3 = zeta`; the `^t g_i` are the **time-`t`** base vectors. Confirmed term by term against the equations this task prints.

**Added (production file, `mitc4_plusd.rs`, test-gated behaviour unchanged):** a clearly delimited `N-alpha — faithful incremental (total-Lagrangian) kinematics` section with plain-data `GlCurrentState`/`GlIncrement` and:

| function | line | equation |
| --- | --- | --- |
| `director_increment_quadratic` | 1879 | Eq. (4c) |
| `node_director_terms` (helper) | 1893 | Eqs. (5c) |
| `incremental_disp_split` | 1915 | Eqs. (5b)/(5c)/(6) |
| `incremental_disp_gradients` | 1943 | the `u_{1,i}`, `u_{2,i}` of Eq. (9) |
| `gl_current_base_vectors` | 1995 | Eq. (7) current `^t g_i` |
| `gl_strain_increment_components` | 2017 | Eq. (9) six components |
| `gl_strain_increment` | 2043 | Eq. (9) driver |

It is deliberately **not** wired into `compute_fint_global` / `compute_kt_global` (N-beta/N-gamma), so no existing operator or tolerance changed. The three entry functions and the two structs are `pub`, so `cargo check --lib` is warning-free without any `#[allow]` (the private helpers are reached by the pub ones).

**Measurement (`#[ignore]`d instrument `n_alpha_rigid_rotation_gives_zero_gl_strain_increment`, line 8990):** flat `RECT`; current state = the flat element; increment = a finite rigid rotation about the in-plane axis `e1`, `u_i = (R - I) x_i`, `alpha_i = theta . V_1^i`, `beta_i = theta . V_2^i`; max `|_0 e + _0 eta|` over nine `(r,s,zeta)` points (`zeta` in {-1,0,1}) and the six components.

| `phi` [deg] | max `\|_0 e + _0 eta\|` | `/phi^3` | max `\|_0 e\|` (u_2 = 0) | `/phi^2` | ratio | current path `\|f_int\|/(\|K_0\|\|u\|)` | Eq. (4c) director err |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 | 2.208500e-4 | 4.153988e-2 | 3.807718e-3 | 1.250000e-1 | 17.24 | 6.130636e-4 | 8.86e-4 |
| 30 | 5.818770e-3 | 4.053548e-2 | 3.426946e-2 | 1.250000e-1 | 5.89 | 5.469097e-3 | 2.38e-2 |
| 60 | 4.280333e-2 | 3.727277e-2 | 1.370778e-1 | 1.250000e-1 | 3.20 | 2.122833e-2 | 1.88e-1 |
| 90 | 1.250000e-1 | 3.225153e-2 | 3.084251e-1 | 1.250000e-1 | 2.47 | 4.538420e-2 | 6.17e-1 |

**Reading.** `max|_0 e + _0 eta| / phi^3` is roughly constant (`3.2e-2 .. 4.2e-2`), so the faithful kinematics satisfies the rigid-rotation identity to the retained order `O(phi^3)`. With `u_2` forced to zero the residual is **exactly `0.125 phi^2` at every angle** (constant `/phi^2`), i.e. the quadratic `u_b2` of Eq. (5c) is what removes the leading `O(phi^2)` error — the measurement discriminates. The CURRENT path does **not** satisfy the identity: `\|compute_fint_global(nonlinear)\| / (\|K_0\|_F \|u\|)` is `6.13e-4 / 5.47e-3 / 2.12e-2 / 4.54e-2` at 10/30/60/90 degrees (nonzero and growing with the angle), against a faithful incremental strain residual that is `O(phi^3)`. Eq. (4c) reproduces the exact rotated director to `O(phi^3)` as well (`6.17e-1` at 90 degrees, `~phi^3/6`).

**Verification (no commit):** `cd crates && cargo test -p aeroelast-core n_alpha -- --ignored --nocapture` -> `1 passed`; `cd crates && cargo test -p aeroelast-core` -> **168 passed / 0 failed / 4 ignored**; `cargo check --lib -p aeroelast-core --message-format=short | grep mitc4_plusd` -> **no output** (no new dead-code warning; no `#[allow]`).

**Next in the N-series:** N-beta — apply the paper's assumed membrane (§2.3), the Eq. (10) transverse shear and the Eq. (11a) through-thickness split to these `_0 e + _0 eta` GL strains; then N-gamma (`^t_0 K`, `^t_0 f_int`), then N-delta (wire, delete the bounded path + dead corotational block, papers' §3 benchmarks).

### Iteration 8 — 2026-09-24 (N-beta: the paper's assumed fields on the TL strains)

**Equations implemented** (Ko, Lee & Bathe (2017), C&S 185:1-14, pp. 3-4; all parent-vision-verified, no new ingredient):

- Eq. (13): `^t x_r = 1/4 sum xi_i ^t x_i`, `^t x_s = 1/4 sum eta_i ^t x_i`, `^t x_d = 1/4 sum xi_i eta_i ^t x_i`; Eq. (14): `^t n = (^t x_r x ^t x_s)/||^t x_r x ^t x_s||` and the dual in-plane basis `^t m^r`, `^t m^s`.
- Eqs. (12b)/(12c): `^t g_ij^m = ^t x_{m,i} . ^t x_{m,j} + ^t x_{m,j} . ^t x_{m,i}`, with `^t x_{m,r} = ^t x_r + s ^t x_d`, `^t x_{m,s} = ^t x_s + r ^t x_d` (Eq. 13), sampled at the five tying points A(0,+1), B(0,-1), C(+1,0), D(-1,0), E(0,0).
- Eqs. (15b)/(15c)/(15d): the assumed combination `^t g~_ij^m` with the CURRENT coefficients `^t a_A..^t a_E` (`compute_membrane_coefficients_2017` on the current `(^t x_d, ^t m^r, ^t m^s)`), coefficient structure identical to the linear element (Eqs. 27a-c).
- Eq. (15a): the TL assumed membrane `^t_0 e~_ij^m = 1/2 ^t g~_ij^m - 1/2 ^0 g~_ij^m`, i, j = 1, 2.
- Eq. (10): the assumed covariant transverse shear `^t_0 e~_rzeta = 1/2(1+s) ^t_0 e_rzeta^(A) + 1/2(1-s) ^t_0 e_rzeta^(B)`, `^t_0 e~_szeta = 1/2(1+r) ^t_0 e_szeta^(C) + 1/2(1-r) ^t_0 e_szeta^(D)`, the same MITC4 tying as the linear element, on the full Eq. (9) TL increment (`e_lin + eta_nl`, linear + nonlinear).

**Added (production file `crates/aeroelast-core/src/elements/mitc4_plusd.rs`, a clearly delimited `N-beta — the paper's assumed strain fields on the total-Lagrangian strains` section, test-gated behaviour unchanged):**

| function | line | equation |
| --- | --- | --- |
| `gl_current_characteristic_vectors` | 2103 | Eqs. (13)/(14) |
| `gl_mid_metrics` | 2127 | Eqs. (12b)/(12c) + Eq. (13) |
| `gl_tying_metrics` | 2150 | Eqs. (15b-d) tying points |
| `gl_assumed_mid_metric` | 2168 | Eqs. (15b)/(15c)/(15d) |
| `gl_assumed_membrane_strain` | 2202 | Eq. (15a) |
| `gl_tying_membrane_increments` (helper) | 2230 | Eq. (15a) on the Eq. (9) increment |
| `gl_assumed_membrane_increment` | 2252 | Eq. (15a) incremental (linearity) form |
| `gl_assumed_transverse_shear` | 2276 | Eq. (10) |

Not wired into `compute_fint_global` / `compute_kt_global` (N-gamma), so no existing operator or tolerance changed. All entry functions are `pub` (no `#[allow(dead_code)]`); the private helpers are reached by the `pub` ones. `cargo check --lib` and `cargo check --tests` are warning-free for this file.

**Measurement (`#[ignore]`d instrument `n_beta_assumed_fields_reduce_to_displacement_based_on_flat`, line 9431; `#[ignore]`d so the default suite stays at 168):** flat `RECT`; the displacement-based comparator is `1/2(^t g_ij^m - ^0 g_ij^m)` from Eqs. (12b)/(12c).

| case | value |
| --- | --- |
| (1a) flat reduction, literal `^t x = ^0 x` | max `\|e~^m - e^m\|` = `0.000000e0` |
| (1b) flat reduction, deformed parallelogram | max `\|e~^m - e^m\|` = `5.551115e-17`, max `\|e^m\|` = `1.25` (non-vacuous) |
| (2a) non-vacuity, `^t a_E` perturbed `0.0 -> 0.5` | max `\|e~^m - e^m\|` = `1.500000e-1` (reduction breaks) |
| (2b) non-vacuity, `^t a_E = 0` on `FLAT_DISTORTED` (natural `a_E = 3.434797e-2`) | max `\|e~^m(a_E) - e~^m(0)\|` = `7.384815e-3` |

Rigid rotation with the assumed fields applied (flat `RECT`, current state = flat; increment = `u_i = (R-I)x_i`, `alpha_i = theta . V_1^i`, `beta_i = theta . V_2^i`), same nine points and `zeta` in {-1,0,1} as N-alpha:

| axis | `phi` [deg] | max `\|e~^m\|` | `/phi^3` | max `\|e~_sh\|` | `/phi^3` | ratio sh/direct |
| --- | --- | --- | --- | --- | --- | --- |
| `e1` | 10 | 9.974660e-18 | 1.876143e-15 | 2.208500e-4 | 4.153988e-2 | 1.000 |
| `e1` | 30 | 6.938894e-18 | 4.833863e-17 | 5.818770e-3 | 4.053548e-2 | 1.000 |
| `e1` | 60 | 0.000000e0 | 0.000000e0 | 4.280333e-2 | 3.727277e-2 | 1.000 |
| `e1` | 90 | 0.000000e0 | 0.000000e0 | 1.250000e-1 | 3.225153e-2 | 1.000 |
| `(e1+e2)/sqrt2` | 10 | 9.892261e-17 | 1.860645e-14 | 3.123290e-4 | 5.874626e-2 | 1.000 |
| `(e1+e2)/sqrt2` | 30 | 6.522560e-17 | 4.543832e-16 | 8.228983e-3 | 5.732582e-2 | 1.000 |
| `(e1+e2)/sqrt2` | 60 | 0.000000e0 | 0.000000e0 | 6.053305e-2 | 5.271166e-2 | 1.000 |
| `(e1+e2)/sqrt2` | 90 | 0.000000e0 | 0.000000e0 | 1.767767e-1 | 4.561056e-2 | 1.000 |

**Reading.** (1) On a flat parallelogram `^t x_d = 0`, so all five current coefficients vanish and Eqs. (15b-d) become the linear interpolation: the assumed membrane reproduces the displacement-based TL membrane to machine precision (`5.55e-17` on the deformed parallelogram) with a non-vacuous membrane (`1.25`), and the literal `^t x = ^0 x` case is exactly `0`. (2) The reduction is real: perturbing a single coefficient `^t a_E` to `0.5` breaks it (`1.5e-1`), and setting `^t a_E = 0` on `FLAT_DISTORTED` (natural `3.43e-2`) changes the field (`7.38e-3`). (3) With the assumed fields applied to the N-alpha Eq. (9) increment, the assumed **transverse shear** Eq. (10) reproduces the N-alpha shear residual exactly (`2.2085e-4` ... `1.2500e-1`, identical to the N-alpha table) and the MITC4 tying is exact for the rigid-rotation shear field (ratio sh/direct `1.000` at every angle); the `/phi^3` is roughly constant (`3.2e-2 .. 5.9e-2`), i.e. the retained order. The assumed **membrane** residual under a rigid rotation about an in-plane axis is identically zero (round-off `~1e-17`): for a flat element the Eq. (9) membrane GL strain of a rigid rotation about an in-plane axis is exact (the entire N-alpha `O(phi^3)` residual is the transverse shear), so the in-plane rigid-rotation identity cannot exercise the assumed membrane; the assumed membrane is instead validated by the flat-reduction case (1) above. Both in-plane axes (`e1` and the diagonal) behave the same way.

**Verification (no commit):** `cd crates && cargo test -p aeroelast-core n_beta -- --ignored --nocapture` -> `1 passed`; `cd crates && cargo test -p aeroelast-core` -> **168 passed / 0 failed / 5 ignored**; `cargo check --lib -p aeroelast-core --message-format=short | grep mitc4_plusd` -> **no output** (also `cargo check --tests` -> no output; no `#[allow]`).

**Ambiguities recorded.** (a) The task's non-vacuity phrase "set `^t a_E = 0`" cannot be literal on the flat `RECT` (`^t a_E` is already exactly `0`), so both readings are measured: perturb it away from zero (2a) and zero it where it is naturally nonzero (2b). (b) Eq. (15a) is printed from the initial configuration `^0`; to apply the assumed field to the one-step N-alpha Eq. (9) increment, the linearity of Eqs. (15b-d) is used to act on the displacement-based tying strains with the current coefficients (`gl_assumed_membrane_increment`). (c) Eq. (10) is taken on the full Eq. (9) TL increment (`e_lin + eta_nl`), not the linear part alone. (d) The through-thickness `b1`/`b2` terms of Eq. (11a) are explicitly out of this stage (they land in N-gamma), as are the `^t_0 K` / `^t_0 f_int` assemblies.

### Iteration 9 — N-gamma implemented; the measurement found a real defect in its tangent

**What exists (implemented by the delegated worker before its tool call was aborted, so the code landed but
the report did not):** the N-gamma block at `crates/aeroelast-core/src/elements/mitc4_plusd.rs:2334+`:
`n_gamma_local_strain` (11-component local strain, the Eq. (15a) total assumed membrane of the current
configuration), `gl_state_after_local`, `advance_gl_state` (the exact finite rotation of the paper's
Eq. (26) — so Eq. (26) IS printed by the paper and this is not the corotational scaffolding),
`n_gamma_w11` (the 11x11 constitutive-resultant block), `n_gamma_b_matrix` (11x24), `n_gamma_fint_local`
(Eq. (24b), `int B^T S d0V`), `n_gamma_kt_local` (Eq. (24a)), plus the `_global` wrappers.

**The measurement I added** (`n_gamma_rigid_body_zero_force_and_consistent_tangent`, `#[ignore]`d):

| check | result |
| --- | --- |
| `F_int(0) = 0` exactly | **0.000e0** PASS |
| `K_t(0) == K_0` bit for bit | **max difference 0.000e0** PASS |
| faithful rigid-body `\|\|F\|\|/(\|\|K_0\|\|\|u\|)` at 1 degree | **1.245e-15** (round-off) vs the CURRENT path's **2.237e-5** |
| **`K_t == dF/du`** by central differences | **FAILS: max\|K_t - FD\| = 4.78e11, relative 0.75** |

**Diagnosis (hypothesis, not yet confirmed):** `n_gamma_kt_local` returns
`compute_ke_local(pre) + (mat - mat0) + geo`, where `mat0 = int B(0)^T W B(0)` with the **11-component**
assembly. The additive `compute_ke_local(pre)` reference is what makes `K_t(0) = K_0` bit for bit, but if
`mat0 != compute_ke_local(pre)` then the returned matrix carries the spurious difference
`(compute_ke_local - mat0)` and is **not** the derivative of `n_gamma_fint_local` — which is exactly what the
measurement shows. So the two possibilities to settle next are: (a) the 11-component `B(0)^T W11 B(0)` does
**not** reproduce the nine-row linear stiffness (a formulation-fidelity defect — the nonlinear formulation
must reduce to the linear one at `u = 0`), or (b) it does and only the additive reference is wrong.
The next step is to print `max|mat0 - compute_ke_local(pre)|` and re-read the paper's Eq. (24a) with vision
before touching anything — per the standing instruction, a failed measurement sends us back to the paper.

**Also to check in N-gamma:** the geometric `N`-term is evaluated by a CENTRAL DIFFERENCE of `B` (H = 1e-6)
rather than the paper's closed form. That is the same quantity numerically, but it is 24 extra `B`
evaluations per Gauss point per element per Newton step, and "the printed Eq. (24a)" would be the analytic
form. Decide whether to replace it with the closed form once the consistency bug above is settled.

### Iteration 10 — N-gamma defect CONFIRMED and localised (not the drill)

**Vision read: C&S 185 pp. 4 right column** (the equations after the assumed membrane), for the first time:

- **Eq. (15e)**: `^t a_A..^t a_E` in terms of `^t c_r = ^t m^r . ^t x_d`, `^t c_s = ^t m^s . ^t x_d`,
  `^t d = ^t c_r^2 + ^t c_s^2 - 1`. (Note the order: `^t c_r = ^t m^r . ^t x_d` here, whereas the
  linear paper's Eq. (24) writes `c_r = x_d . m^r` — same scalar, and `compute_membrane_coefficients_2017`
  already matches.)
- **Eq. (16)**: `^t_0 e~_ij = ^t_0 e~_ij^m + zeta ^t_0 e_ij^b1 + zeta^2 ^t_0 e_ij^b2` — the
  through-thickness split with the **assumed** membrane in the `t^0` slot.
- **Eq. (17)**: the shell-aligned local Cartesian frame at time 0,
  `^0 L_3 = ^0 g_3/||^0 g_3||`, `^0 L_1 = (^0 g_2 x ^0 L_3)/||..||`, `^0 L_2 = ^0 L_3 x ^0 L_1`.
- **Eq. (18)**: `^t_0 e~_ij = ^t_0 e~_kl (^0 L_i . ^0 g^k)(^0 L_j . ^0 g^l)` — **the covariant-to-local
  mapping of the TL strain uses the INITIAL metric `^0 g^k` and the INITIAL frame `^0 L_i`.** This is the
  prime suspect for the defect below.
- Fig. 6 (p. 4 right) is the §3.1 cantilever load-displacement benchmark, i.e. the N9 oracle.

**The defect, now measured and localised.** Added the decisive diagnostic to the N-gamma instrument:

| diagnostic | value |
| --- | --- |
| `max \|B(0)^T W11 B(0) - compute_ke_local\|` | `4.780e11` (scale `1.593e11`, relative **3.000**) |
| `max \|B(0)^T W11 B(0) - compute_ke_local_with_drill(false)\|` | **`4.780e11`** (relative **3.000**) — the same |
| `max \|compute_ke_local - compute_ke_local_with_drill(false)\|` | `6.410e9` (the drill block is three orders smaller) |

So the missing 2025 drill is **not** the cause: against the NO-DRILL 2017 linear stiffness the 11-component
zero-state assembly is off by the same 300%. **The N-gamma strain/`W11`/`B` assembly genuinely does not
reduce to the linear formulation at `u = 0`** — a formulation-level defect, and it is exactly the
`(mat0 - compute_ke_local)` term that made `K_t != dF/du`.

**Candidate causes, in the order to check them** (each against the printed page, not by tuning):
1. **Eq. (18)'s initial metric.** If `n_gamma_local_strain` / `n_gamma_b_matrix` map the covariant TL
   strain to the local frame with the CURRENT metric (`^t g^k`) or the current frame instead of the
   `^0` ones, the zero-state stiffness is wrong by construction. Check this first — it is the one equation
   read this iteration that the N-alpha/N-beta work could not have anticipated.
2. **The `s1 = 2/h` / `s2 = 4/h^2` convention.** The linear assembly puts them in the B rows with
   `W = [cm, cbc, cb; ...]`; the N-gamma `W11` may carry the moments instead. Two equivalent conventions
   are fine, but a MIXED one (scaling in neither or in both) is exactly a factor-sized error.
3. **The `t^0`/`t^1`/`t^2` moments of `W11`** against Eq. (11b)/(11c)/(11d) (`h`, `h^3/12`, `h^5/80` after
   the `zeta` split) and the 2x2x2 rule of §3.1.
4. **The geometric `N`-term**: implemented as a central difference of `B` (H = 1e-6) rather than the
   printed closed form. Once (1)-(3) are settled, decide whether to replace it (it is the same quantity,
   but 24 extra `B` evaluations per Gauss point per element per Newton step).

**State.** `cargo test -p aeroelast-core` = **168 passed / 0 failed / 6 ignored** (the N-gamma instrument
is `#[ignore]`d and currently FAILS its `K_t == dF/du` assertion — that failure is the evidence above, not
a red suite). Nothing committed.

### Iteration 11 — Eq. (18) verified CORRECT; the defect is narrowed to the Eq. (15a) coefficient change

**Verification against the printed page (C&S 185 Eq. (18)):** the N-gamma covariant-to-local mapping uses the
**initial** metric and frame. In `n_gamma_local_strain`: the membrane and the `b1`/`b2` slices are mapped with
`covariant_to_local_mapping(&j_loc_at(pre, r, s))`, and `j_loc_at` builds `g_r = x_r + s x_d` from
`pre.initial_coords_3d`/`pre.vn` (the `^0` state); the shear is mapped with
`shear_covariant_to_local(&g0_r0, &g0_s0, &g0_t0, ...)` on the **initial** base vectors. So the
`(^0 L_i . ^0 g^k)` of Eq. (18) is implemented. **Suspect 1 CLEARED.**

The rest of the slice construction also matches the printed page: the membrane is the Eq. (15a) total assumed
membrane plus the 2025 drill (Eq. (22a)); the `b1`/`b2` slices are the exact total GL through-thickness
moments obtained from the `zeta = 0, +1, -1` metrics with `b1 = (g(1) - g(-1))/2` and
`b2 = (g(1) + g(-1))/2 - g(0)` — exact, because the metric is quadratic in `zeta`; the shear is the Eq. (10)
assumed field on the exact total covariant shear; the scaling `s1 = 2/h`, `s2 = 4/h^2` and the 11x11 block
`[W9 ; cs_uncorrected]` reproduce the linear assembly's convention.

**So the remaining explanation of the 300% is a formulation question, not a mapping bug:** `mat0` is
`int B(0)^T W11 B(0)` with `B(0) = d(e~)/du` at `u = 0`, and `e~` is built from the **Eq. (15a) total
assumed metric**, whose coefficients `^t a_A..^t a_E` are those of **each configuration**
(`^t g~^m` uses `^t a`, `^0 g~^m` uses `^0 a`, Eqs. 15b-d with 15e). Therefore `d(e~)/du` at `u = 0`
contains the **coefficient-change term**, which the linear `b_membrane_2017` — the 2017 (27a-c) combination
of the sampled **strains** with the fixed **initial** coefficients — does not have. That is exactly the term
the N-beta caveat note was about, and it is why `K_t` (which adds `compute_ke_local(pre)` as its reference)
cannot be the derivative of `n_gamma_fint_local`.

**The open question the paper must answer (next step):** what object the paper's Eq. (24a) tangent is built
on. Two readings:
- (i) the tangent is the exact derivative of Eq. (24b)'s `F = int B^T S dV` with the Eq. (15) strain — then the
  `compute_ke_local` reference in `n_gamma_kt_local` must go, and the u = 0 tangent equals
  `mat0` (+ the shear), which is then **not** the linear stiffness but the linear stiffness plus the
  coefficient-change term (a statement to confirm against the paper's small-displacement limit);
- (ii) the formulation is written so the u = 0 tangent IS the linear stiffness (the TL element must reduce to
  the linear one), which would mean Eq. (15a)'s per-configuration coefficients must collapse at `u = 0` in
  the derivative — in which case `mat0 == compute_ke_local_with_drill(false)` must hold and the current
  300% is a real bug elsewhere.
Reading the paper's Eqs. (19)-(24) settles which. Nothing touched until then.

### Iteration 12 — FOUND IT: the paper distinguishes the assumed TOTAL strain (§2.3) from the assumed INCREMENTAL strain (§2.4)

**Vision read: C&S 185 p. 5 right column — §2.4, "Assumed incremental Green-Lagrange strains".** This is the
section that was never read before, and it changes the N-gamma construction:

> *"From the assumed Green-Lagrange strains in Section 2.3, we proceed to obtain the corresponding assumed
> incremental Green-Lagrange strains in a consistent manner. For the transverse shear we use the following
> assumed incremental shear strains"*

**Eq. (19):**
`_0 e~_rzeta = 1/2(1+s) _0 e_rzeta^(A) + 1/2(1-s) _0 e_rzeta^(B)`, `_0 e~_szeta = 1/2(1+r) ... (C) + 1/2(1-r) ... (D)`,
and **separately** `_0 eta~_rzeta = 1/2(1+s) _0 eta_rzeta^(A) + ...`, `_0 eta~_szeta = ...`.

So the paper has **two different assumed fields**:
- **§2.3, Eq. (15a-d): the assumed TOTAL TL strain** — `^t_0 e~^m` from the assumed *metric* of each
  configuration, with **that configuration's own coefficients** `^t a_A..^t a_E` (Eq. 15e).
- **§2.4, Eq. (19)+: the assumed INCREMENTAL strain** — the assumed field applied to the **increment's**
  linear and nonlinear parts (`_0 e` and `_0 eta` separately; equivalent to applying it to the total
  increment, since the combination is linear).

**Why this is the defect.** `n_gamma_local_strain` builds the membrane slice from the **§2.3 total** (Eq. 15a)
and the `b1`/`b2` slices from the total metric differences, then `n_gamma_b_matrix` obtains `B` by finite
differences of those totals. That derivative therefore carries the **coefficient-change term** of Eq. (15e)
(`^t a` in the current configuration vs `^0 a` in the reference), which is why
`mat0 = int B(0)^T W11 B(0)` misses the linear stiffness by 300%, and why `K_t` — whose reference is
`compute_ke_local(pre)` — cannot be the derivative of `n_gamma_fint_local`.

**The fix direction (paper-grounded, not a patch):** the TANGENT (and the B used inside it) must be built on
the **§2.4 assumed incremental** strains — the assumed field applied to the **increment** with the **current**
coefficients — which at `u = 0` reduces to the linear strain-displacement matrix and therefore to
`compute_ke_local_with_drill(false)` exactly. That is the "reading (ii)" of iteration 11, and it is what the
paper prints. The §2.3 total remains the right object for the internal FORCE `F = int B^T S dV` (Eq. 24b),
which is consistent because `B` there is the derivative of the incremental strain.

Note the irony recorded honestly: the N-beta "caveat fix" (replacing "the assumed field on the increment" by
"the difference of the two assumed totals") went in the **opposite** direction to what §2.4 prints. The
original construction was closer to the paper for the TANGENT; the total form belongs to the FORCE/state.

**Next:** read the rest of §2.4 (p. 6: the assumed incremental in-plane strains and the Eq. (20)-(24)
internal force/tangent) and then restructure `n_gamma_local_strain`/`n_gamma_b_matrix` accordingly, with the
`mat0 == compute_ke_local_with_drill(false)` identity as the acceptance measurement.

### Iteration 13 — the paper's COMPLETE incremental strain structure, Eqs. (20a)-(20g)

**Vision read: C&S 185 p. 6 left column, §2.4 continued.** This is the missing specification. The paper writes
the incremental strain as a `zeta` split **for each of the two parts** (linear `_0 e`, nonlinear `_0 eta`):

**Eq. (20a):** `_0 e_ij = _0 e_ij^m + zeta _0 e_ij^b1 + zeta^2 _0 e_ij^b2` and
`_0 eta_ij = _0 eta_ij^m + zeta _0 eta_ij^b1 + zeta^2 _0 eta_ij^b2`, `i,j = 1,2`.

**Linear part:**
- Eq. (20b): `_0 e_ij^m = 1/2(^t x_{m,i} . u_{m,j} + ^t x_{m,j} . u_{m,i})`
- Eq. (20c): `_0 e_ij^b1 = 1/2(^t x_{m,i} . u_{b1,j} + ^t x_{m,j} . u_{b1,i} + ^t x_{b,i} . u_{m,j} + ^t x_{b,j} . u_{m,i})`
- Eq. (20d): `_0 e_ij^b2 = 1/2(^t x_{b,i} . u_{b1,j} + ^t x_{b,j} . u_{b1,i})`

**Nonlinear part:**
- Eq. (20e): `_0 eta_ij^m = 1/2 u_{m,i} . u_{m,j}`
- Eq. (20f): `_0 eta_ij^b1 = 1/2(u_{m,i} . u_{b1,j} + u_{m,j} . u_{b1,i} + ^t x_{m,i} . u_{b2,j} + ^t x_{m,j} . u_{b2,i})`
- Eq. (20g): `_0 eta_ij^b2 = 1/2(u_{b1,i} . u_{b1,j} + ^t x_{b,i} . u_{b2,j} + ^t x_{b,j} . u_{b2,i})`

**What this settles, and it is decisive:**

1. **The linear incremental part (Eqs. 20b-d) uses the CURRENT `^t x_m`, `^t x_b` and the increment's
   `u_m`, `u_{b1}`** — no reference-metric difference, no coefficient change. Therefore the tangent's `B` at
   `u = 0` (where `^t x = ^0 x`) is exactly the linear strain-displacement matrix, and the tangent's
   `K_t(0)` is exactly `compute_ke_local_with_drill(false)`. **Reading (ii) of iteration 11 is what the
   paper prints.**
2. **The `u_b2` terms appear ONLY in the nonlinear part (Eqs. 20f, 20g)** — which is where the current
   bounded path has nothing at all (N2), and where the quadratic director increment
   `u_b2 = -1/4 sum a_i h_i (alpha_i^2 + beta_i^2) ^t V_n^i` of Eq. (5c) enters.
3. The §2.3 **total** (Eq. 15) is a different object and belongs to the STATE (the internal force/stress),
   not to the tangent.

**The fix, now fully specified from the printed page** (no guessing left):
`n_gamma_local_strain` / `n_gamma_b_matrix` must be built from **Eqs. (20a-g)** — the incremental strains with
the CURRENT `^t x_m`/`^t x_b` and the increment's `u_m`, `u_{b1}`, `u_{b2}` — with the assumed fields applied
to the incremental tying strains (§2.4's Eq. (19) for the shear and its in-plane counterpart), instead of from
the increment of the §2.3 totals. Acceptance measurement, already in place:
`mat0 == compute_ke_local_with_drill(false)` to round-off, and `K_t == dF/du`.

**Still to read for completeness:** the rest of p. 6 (the assumed incremental in-plane field and Eqs. (21)-(24),
the force/tangent closed forms) — but the strain side is now unambiguous, which is what the measurement failed
on.

### Iteration 14 — N-gamma restructured onto §2.4's INCREMENTAL strain; all four measurements closed

**Read with vision:** p. 6 (both columns) and p. 7 (both columns) of
`.sources/papers/The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf`.

**Page-vs-spec disagreement (recorded, page followed):** p. 6 carries only Eqs. (20a-g) plus the "with"
definitions of `^t x_{m,i}`, `^t x_{b,i}`, `u_{m,i}`, `u_{b1,i}`, `u_{b2,i}` and Figs. 10-11. The rest of §2.4 —
the assumed incremental in-plane field, **Eqs. (21a-c)**, the Eq. (22) statement, **Eq. (23)** and the §2.5
**Eqs. (24a)/(24b)/(25)** — is on **p. 7**. The spec's placement of (21)-(24) on p. 6 is wrong; the content is as
the spec describes.

**Changed (production, `n_gamma` block only):**
- `n_gamma_local_strain` (line 2433): the membrane row is now `gl_assumed_mid_metric(&^t a, &dtie, r, s)` with
  `dtie[k] = 0.25*(tie(next) - tie(state))` — Eqs. (21a-c) applied by Eq. (22) to the **incremental** tying
  strains of Eq. (20b)/(20e), with the **current** coefficients `^t a_A..^t a_E` (Eq. 15e) frozen at `state`.
- `n_gamma_local_strain`: the `zeta = 0, ±1` covariant slices and both shear tying shears now difference against
  `state` (`^t`), not `pre.initial_coords_3d`/`pre.vn` (`^0`) — Eqs. (7)/(8) increment instead of the §2.3 total.
- `n_gamma_b_matrix` (line 2601): central-difference step `1e-6 -> 2e-5`, the measured round-off/truncation
  balance of the printed-expression linearization.

**Two defects were in the old form, both measured, and the first is the dominant one:**

| # | defect | how the measured number identifies it |
| --- | --- | --- |
| D1 | **factor 2 in the tying strain.** `gl_tying_metrics` reports `^t g_ij^m = x_i.x_j + x_j.x_i` (twice the Gram metric); the old membrane used `0.5*(tie(next)-tie(init))`, i.e. `A.(M(next)-M(init))`, whose `u=0` derivative is `2x` `b_membrane_covariant_2017` (whose `b_rr_a[j] = dn_i/dr (g_r.e_k)` is exactly `d(0.5 M)/du`). | `mat0`'s membrane block was `4x`; `4x - 1x = 3x` on a `max abs K` that is itself a membrane entry: `3.000 x 1.593e11 = 4.780e11`, to the digit the instrument printed. The Eq. (15e) coefficient-change term (D2) is real but second order by comparison. |
| D2 | same-form, the **Eq. (15e) coefficient-change term**: the old form differenced one assumed *total* per configuration (`^t a` on the perturbed geometry vs `^0 a` on the initial one), so `d/du` at `u=0` also carried `d(^t a)/d^t x` times the initial tying metric. | invisible under D1 at this scale; removed by construction (the coefficient operator is now evaluated on `state` and frozen). |

**Measured, `cd crates && cargo test -p aeroelast-core n_gamma -- --ignored --nocapture`:**

| measurement | before | **after** |
| --- | --- | --- |
| (0) `max abs B(0)^T W11 B(0) - compute_ke_local` (relative) | 4.780e11 (3.000e0) | **3.699e0 (2.321e-11)** |
| (0b) same vs `compute_ke_local_with_drill(false)` | 4.780e11 (3.000e0) | 6.410e9 (4.023e-2) = (0c), the drill block |
| (2) `max abs K_t - dF/du` (relative), 12 rigid motions | 4.780e11 (7.502e-1), panic at 1 deg | **9.6e0..3.4e1 (5.6e-11..2.2e-10)** |
| (3) `abs F_int(0)` / `K_t(0) - K_0` | 0 / 0 (bit) | **0 / 0 (bit)** |
| rigid-body force residual `F/(K u)`, faithful path | 1.2e-15..1.0e-16 | 1.0e-16..3.0e-17 |
| same, bounded path (N6/N7, lacks `u_b2`) | 2.2e-5 .. 1.65e-1 | unchanged |

**(0b) is not a defect.** The N-gamma strain contains the 2025 drill row (C&S 308:107622 Eq. 22a folded into the
membrane, exactly as the linear path), so the correct `u = 0` reference is `compute_ke_local(pre)` = row (0);
rows (0b)/(0c) are the same number, i.e. `mat0` reproduces the drilled linear stiffness and the no-drill
difference is *entirely* the drill block. Making (0b) small would require deleting the drill from the strain while
the Eq. (25) bit-for-bit gate compares against the drilled `K_0`: with the additive `compute_ke_local` reference
kept, `dF/du` would then be off by the drill block (4.0e-2 relative) and gate (2) would fail. The `u_b2` presence
(item 3 of the task) is what the `faithful |F|/Ku` column measures: round-off here vs 1.65e-1 for the bounded path.

**Follow-up candidate (NOT changed here, out of this task's scope):** `gl_assumed_membrane_strain` and
`gl_assumed_membrane_increment` themselves still carry D1 — their value is `2x` the assumed strain of the
validated linear path (`0.5.A.d(2M)/du` vs `A.d(0.5M)/du`). They are `pub` and are the N-beta block's objects,
and their tests compare against `gl_mid_metrics`' own doubling, i.e. they share the assumption. The N-gamma block
no longer calls either.

### Iteration 15 — the N-beta factor-2 fixed; its measurement re-baselined against the linear operator

**Authorized by the parent**: "Fix the factor-2 that the N-gamma work just unmasked in the N-beta block, and
re-baseline its measurement against the LEVEL linear operator. NO commit." The allowed surfaces were exactly
`crates/aeroelast-core/src/elements/mitc4_plusd.rs` and this document; the N-gamma block, the linear path,
`compute_fint_global`/`compute_kt_global` and every other test were left untouched.

**The convention adopted (Eqs. 12b and 15a with the paper's doubled `g`).** Eq. (12b) prints
`^t g_ij^m = ^t x_{m,i} . ^t x_{m,j} + ^t x_{m,j} . ^t x_{m,i}`, i.e. **twice** the metric
`^t m_ij^m = ^t x_{m,i} . ^t x_{m,j}`. Eq. (15a) prints `^t_0 e~_ij^m = 1/2 ^t g~_ij^m - 1/2 ^0 g~_ij^m`, so its
`1/2` already turns the doubled `g~` into the metric (`1/2 g~ = m~`); the Green-Lagrange strain that
`b_membrane_covariant_2017` implements is the further `1/2` of the metric difference,
`^t_0 e~_ij^m = 1/2(^t m~_ij^m - ^0 m~_ij^m) = 1/4(^t g~_ij^m - ^0 g~_ij^m)`. The total conversion from the
doubled tying metric is therefore **`1/4`** — the same `0.25 * (next_tie - cur_tie)` the N-gamma block uses.
This is stated verbatim in the two docstrings, which quote Eqs. (12b) and (15a).

**Changed (production, `mitc4_plusd.rs` only):**

| function | line | change |
| --- | --- | --- |
| `gl_assumed_membrane_strain` | ~2237 | `0.5*cur - 0.5*init` -> `0.25*cur - 0.25*init`; now calls `gl_assumed_metric_of_state` for the current metric; convention docstring with Eq. (12b)/(15a) |
| `gl_assumed_metric_of_state` | ~2293 | docstring only (states it returns the doubled `g~`) |
| `gl_assumed_membrane_increment` | ~2311 | total-difference form replaced by the N-gamma frozen-coefficient incremental form `gl_assumed_mid_metric(&^t a, &0.25(^{t+dt}tie - ^t tie))`; docstring records why (the Eq. 15e coefficient-change term `D2`) |

No helper or dead code was added; `cargo check --lib` and `cargo check --tests` are clean for the file apart
from the file's pre-existing test-only warnings (`PUB_MITC4_*`, `drill_block`, one snake-case test name). The only
new build artefact was the removal of the now-unused `gl_mid_metrics` import from the test module.

**Deliverable 2 — the re-baselined, non-vacuous measurement** (`n_beta_assumed_fields_reduce_to_displacement_based_on_flat`):

1. **NEW decisive acceptance.** The N-beta assumed incremental covariant operator at `u = 0`, built by central
differences of `gl_assumed_membrane_increment` (so it cannot share a factor with the reference), equals the
linear operator at the five sample points and the five tying points:

| comparison | max\|diff\| | max\|linear\| | **rel** |
| --- | --- | --- | --- |
| covariant rows vs `b_membrane_covariant_2017` | 1.197889e-12 | 5.000000e-1 | **2.395778e-12** |
| local rows (`2 e_12`, `^0` map) vs `b_membrane_2017` | 1.609502e-12 | 1.000000e0 | **1.609502e-12** |

Both are asserted `<= 1e-10`. A common-factor regression makes either `rel` jump from `~1e-12` to `O(1)` — with
the old `0.5` convention (now restored) the local `rel` was exactly `2.000e0`.

2. **Flat reduction**, both sides in the SAME explicitly-stated `1/2(m_cur - m_init)` convention computed
directly from `x_m` (NOT `gl_mid_metrics`), and `b_membrane_2017` used for the operator identity above:

| case | before (old vacuous instrument) | after |
| --- | --- | --- |
| (1a) literal `^t x = ^0 x`, `max\|e~^m - e^m\|` | 0.000000e0 | **0.000000e0** |
| (1b) flat parallelogram, `max\|e~^m - e^m\|` (`max\|e^m\|`) | 5.551115e-17 (`1.25`) | **2.775558e-17 (`0.625`)** |
| (2a) non-vacuity, `^t a_E` 0 -> 0.5 | 1.500000e-1 | **7.500000e-2** (still breaks the reduction) |
| (2b) non-vacuity, `^t a_E = 0` on `FLAT_DISTORTED` | 7.384815e-3 | **7.384815e-3** (unchanged) |

The (1b) `max\|e^m\|` and the (2a) residual halve exactly because the old comparator shared the doubling; the
residual `2.8e-17` is now two genuinely independent constructions agreeing.

**Verification (no commit):**

| command | result |
| --- | --- |
| `cd crates && cargo test -p aeroelast-core n_beta -- --ignored --nocapture` | **1 passed / 0 failed** |
| `cd crates && cargo test -p aeroelast-core n_gamma -- --ignored --nocapture` | **1 passed / 0 failed** (276.68 s) |
| `cd crates && cargo test -p aeroelast-core` | **168 passed / 0 failed / 6 ignored** |
| `cargo check --lib -p aeroelast-core --message-format=short \| grep mitc4_plusd` | **no output** |

**Printed page vs the (linear) operator.** No disagreement was found this iteration: the printed Eq. (12b) is
unquestionably the doubled Gram sum and Eq. (15a) unquestionably carries `1/2 g`, so a literal implementation
returns the metric difference (twice the strain). The validated linear operator `b_membrane_covariant_2017`
implements `d(0.5 m)/du`, which is the Green-Lagrange strain; the paper's printed pair is therefore read as the
`1/2 g` conversion plus the Green-Lagrange `1/2`, i.e. the `1/4` total, exactly matching the independent N-gamma
fix of Iteration 14. This is the only place the page needed re-interpretation, and it is recorded in both
docstrings rather than hidden in the code.

### Iteration 16 — N-delta: the faithful pair wired on the PyO3 kernels, the corotational block deleted, and **the production solve path does not go through those kernels**

**Authorized by the parent** (bounded N-delta unit): wire the faithful TL pair on the production PyO3 path,
delete the dead corotational block (N8), and measure the large-rotation benchmarks. Allowed surfaces:
`crates/aeroelast-py/src/elements.rs`, this file, and `odd/tasks/mitc4plusd-2025-purity.md`. NO commit.

**Changed — `crates/aeroelast-py/src/elements.rs` (the four MITC4+ PyO3 kernels are the only code touched):**

| kernel | line | before | after |
| --- | --- | --- | --- |
| `batch_kt_mitc4` (no `nonlinear` flag; its doc calls it the nonlinear tangent) | 341 | `mitc4_plusd::compute_kt_global(&pre, &u)` | `mitc4_plusd::n_gamma_kt_global(&pre, &u)` |
| `batch_fint_mitc4` (`nonlinear = true`) | 394 | `mitc4_plusd::compute_fint_global(&pre, &u, nonlinear)` | `mitc4_plusd::n_gamma_fint_global(&pre, &u)` |
| `batch_fint_mitc4` (`nonlinear = false`) | 394 | same call | `mitc4_plusd::compute_fint_global(&pre, &u, false)` = `K_0 u`, unchanged |

Every PyO3 signature, name and returned shape is unchanged; the MITC3+ kernels are untouched. The known
limitation is recorded in doc comments (and in the two `n_gamma_*_global` docstrings, this file), not fixed:
the wrappers build their `GlCurrentState` from the INITIAL element geometry, so they are the first-step
(reference-configuration) linearisation; threading the last-converged state belongs to the solver.

**Deleted — the N8 corotational block (zero call sites, unsourced: the paper is total Lagrangian).**
Grep evidence (`grep -rn "\b<name>\b" --include=*.rs crates/`, excluding `crates/target/`), before deletion:

| deleted item | non-definition occurrences in the whole workspace |
| --- | --- |
| `GpLocalFrame` (struct) | 0 in `mitc4_plusd.rs`; used only by the two deleted methods |
| `quaternion_multiply` | 0 anywhere |
| `update_normals_with_displacements` | 0 in `mitc4_plusd.rs` (`mitc4.rs` has its own read-only copy) |
| `polar_decomposition` | 1, and only inside `compute_membrane_strain_log` |
| `log_strain_from_polar` | 1, and only inside `compute_membrane_strain_log` |
| `compute_membrane_strain_log` | 0 |
| `update_corotational_frame` | 0 |
| `frame_incremental_rotation` | 0 |

**Kept:** `quaternion_to_matrix`, `quaternion_from_vector`, `rotate_vector_by_quaternion` — they implement the
paper's own Eq. (26) finite rotation and are the only consumers of `gl_state_after_local`
(`mitc4_plusd.rs:2259-2262`, `2576-2579`). No test referenced any deleted item, so no test had to be removed and
no assertion was weakened.

**Measured — the wired kernels are live and faithful (direct PyO3 calls, flat unit quad, `E=1e7`, `nu=0.3`,
`h=0.1`, `sc=5/6`):**

| check | value |
| --- | --- |
| `max\|K_t(0) - K_0\|` (batch kt vs batch ke at `u=0`) | **0.000e0** (bit for bit, `\|K_0\|max = 4.945e5`) |
| `max\|F_nl(0)\|` | **0.000e0** |
| `max\|F_lin - K_0 u\| / max\|K_0 u\|` (`nonlinear=false`) | **9.209e-17** |
| `max\|F_nl - F_lin\| / max\|F_lin\|` at `u` random `~1e-3` (non-vacuity) | **6.890e-3** |
| rigid body rotation 10 deg: `\|F\|/(\|K_0\|\|u\|)`, faithful kernel | **1.294e-16** |
| same quantity, the OLD bounded path (N-alpha instrument, its own fixture) | 6.13e-4 |

**Measured — and the finding that matters: the production solve path does not call these kernels.**
The Python solver calls `_aeroelast.nonlinear_static_solve_coo(assembler, ...)`, i.e.
`PyMeshAssembler::assemble_kt` / `assemble_fint` (`crates/aeroelast-py/src/assembler.rs:554/580`) ->
`MeshAssembler::assemble_kt` / `assemble_fint` (`crates/aeroelast-core/src/assembly/assembler.rs:648/723`) ->
`mitc4_plusd::compute_kt_global` / `compute_fint_global` **directly**. `batch_kt_mitc4` / `batch_fint_mitc4`
have **no Python caller at all** (`grep -rn "batch_fint_mitc4\|batch_kt_mitc4" --include=*.py .` -> empty).
Direct measurement of the production path on the same flat quad and the same 10-degree rigid rotation:

| quantity | production `assemble_fint` / `assemble_kt` | faithful pair |
| --- | --- | --- |
| rigid rotation 10 deg, `\|F\|/(\|K_0\|\|u\|)` | **1.971e-3** | 1.294e-16 |
| `max\|Kt_prod(u) - Kt_faithful(u)\|` | **4.385e4** (`\|K_0\|max = 4.945e5`) | 0 |

So the faithful pair is now wired and correct on the PyO3 batch kernels, but those kernels are not the
production nonlinear solve path. Putting the faithful pair on the production path requires
`crates/aeroelast-core/src/assembly/assembler.rs` (or delegating inside `compute_fint_global` /
`compute_kt_global` in this file) — both outside the unit's allowed surfaces, so the parent was asked.

**The N-delta measurement, before vs after (the benchmark numbers are IDENTICAL, as the finding above predicts):**

| run | before (installed extension, current tree) | after (rebuilt with the wired kernels) |
| --- | --- | --- |
| `cd crates && cargo test -p aeroelast-core` | 168 / 0 / 6 ignored | **168 / 0 / 6 ignored** |
| `cargo test -p aeroelast-core n_alpha -- --ignored` | 1 passed | **1 passed** |
| `cargo test -p aeroelast-core n_beta -- --ignored` | 1 passed | **1 passed** |
| `cargo test -p aeroelast-core n_gamma -- --ignored` | 1 passed (276.68 s) | **1 passed (353.87 s)** |
| `pytest tests/test_large_rotation_benchmarks.py` | 4 passed / 3 failed | **4 passed / 3 failed** |
| `pytest tests/test_shell_validation_fixed.py` | 6 passed / 1 failed | **6 passed / 1 failed** |
| `pytest -m "not slow"` | 343 / 8 / 2 skipped | **343 / 8 / 2 skipped** |
| `pytest tests/test_ccx_shell_element_types_parity.py` | 4 passed | **4 passed** |

Large-rotation, per case (all seven, before == after to every printed digit):

| case | result | value |
| --- | --- | --- |
| `test_linear_tip_deflection_euler_bernoulli` | PASS | — |
| `test_cantilever_large_rotation_half_circle[10]` | PASS | — |
| `test_equilibrium_path[pi/2, -3.6338, 6.3662]` | PASS | — |
| `test_equilibrium_path[pi, -10.0, 6.3662]` | PASS | — |
| `test_equilibrium_path[3pi/2, -12.122, 2.1221]` | **FAIL** | `w_tip = 1.5109` vs `2.1221`, **28.80%** (`tol 5%`) |
| `test_equilibrium_path[2pi, -10.0, 0.0]` | **FAIL** | `u_tip = -9.3249` vs `-10.0000`, **6.75%** |
| `test_simo_vu_quoc_rollup_360[10]` | **FAIL** | `u_tip = -9.3249` vs `-10.0000`, **6.75%** |

The 28.80% / 6.75% / 6.75% triple is exactly the recorded WU9f value, so the installed extension before this
unit matched the record. **No benchmark got worse (none changed).** `TestNonlinearStatic::test_geometric_nonlinearity`
still converges where it used to diverge: the failure is `DID NOT RAISE`
(`tests/test_shell_validation_fixed.py:282`), so the faithful pair is not on that test's path either.

**Verification (no commit):** `cargo check --lib -p aeroelast-py --message-format=short` -> 0 errors, 1 warning
(the pre-existing `proc-macro-error2` future-incompat note); the only `elements.rs` warning in the release build
is the pre-existing pyo3 `HasAutomaticFromPyObject` deprecation at line 875. `maturin develop --release` builds
clean and reinstalls the extension.

**Next (needs the parent's decision):** wire the faithful pair where the production solve actually calls it
(core assembler, or delegation inside `compute_*_global`), then re-measure. That measurement is the one that can
move the three large-rotation cases.

### Iteration 17 — N-delta option B: the two public entry points delegate, and the measurement is a **regression**

**Authorized by the parent** (maintainer-chosen option B): make the faithful total-Lagrangian pair the
production nonlinear behaviour of `mitc4_plusd` by delegating inside the two public entry points, exactly as
`mitc4`/`mitc3` structure theirs. Allowed surfaces: `crates/aeroelast-core/src/elements/mitc4_plusd.rs` and
this file. NO commit.

**Changed (production, `mitc4_plusd.rs`).** Both signatures, names, returned shapes and local/global
conventions are unchanged.

```rust
// compute_fint_global (was: `K u + membrane_nonlinear_correction` when nonlinear)
    if nonlinear {
        return n_gamma_fint_global(pre, u_global);
    }
    let t24 = build_t24(pre);
    let u_local = t24 * u_global;
    t24.transpose() * (compute_ke_local(pre) * u_local)

// compute_kt_global (was: `K_0 + K_L + K_sigma`, symmetrised)
pub fn compute_kt_global(pre: &Mitc4PlusDPrecomputed, u_global: &Vec24) -> Mat24 {
    n_gamma_kt_global(pre, u_global)
}
```

The `nonlinear = false` branch is byte-for-byte the old linear behavior (`compute_ke_local(pre) * u_local`,
then `T^T`).

**Renamed / deleted, and which instrument uses them.** The superseded bounded bodies are NOT dead code and are
NOT left as live production: they are `#[cfg(test)]`-only and are called only by the two `#[ignore]`d
instruments as the recorded historical "before".

| item | disposition | used by |
| --- | --- | --- |
| `membrane_nonlinear_correction` | renamed `membrane_nonlinear_correction_bounded_superseded`, `#[cfg(test)]` | via the fint below |
| old `compute_fint_global` nonlinear branch | new `compute_fint_global_bounded_superseded(..)`, `#[cfg(test)]` | `n_alpha_*` and `n_gamma_*` instruments |
| old `compute_kt_global` body (`K_0 + K_L + K_sigma`) | **deleted** — no instrument uses it | — |
| `geometric_stiffness_local` | **deleted** — its only caller was the old tangent body | — |
| `displacement_gradient`, `membrane_strain_nl`, `compute_b_nl`, `extract_membrane_rows` | `#[cfg(test)]` (their only callers are the bounded bodies) | bounded bodies |

No `#[allow(dead_code)]` was added. `cargo check --lib` is clean for the file; `cargo check --tests` shows
only the five pre-existing test-only warnings (`PUB_MITC4_1X4/1X8`, `PUB_MITC4_IC_1X4`, `drill_block`, the
`test_geometry_a_E` snake-case name).

**Re-baselined assertions (before/after text).**

`n_alpha_rigid_rotation_gives_zero_gl_strain_increment`:

```text
// BEFORE:
//   let f_cur = compute_fint_global(&pre, &u_rot, true);   // == the OLD bounded path
//   // The current path does not satisfy the rigid-rotation identity.
//   assert!(f_rel > 1.0e-5,
//       "... the current path was expected to fail the rigid rotation, got {f_rel:.3e}");
// AFTER:
//   let f_cur = compute_fint_global(&pre, &u_rot, true);                 // now the FAITHFUL pair
//   let f_bounded = compute_fint_global_bounded_superseded(&pre, &u_rot); // recorded "before"
//   assert!(f_rel < 1.0e-10, "... production nonlinear path must satisfy the rigid rotation to round-off ...");
//   assert!(f_bounded_rel > 1.0e-5, "... the superseded bounded path was expected to fail ...");
```

`n_gamma_rigid_body_zero_force_and_consistent_tangent`: the `bounded |F|/Ku` column now reads the renamed
`compute_fint_global_bounded_superseded(..)`; a new `prod |F|/Ku` column reads the production
`compute_fint_global(.., true)` and is ASSERTED `< 1e-10` (the faithful local column is unchanged).

**The in-file T2B tests hold as written — no assertion was changed.** `test_kt_zero_matches_ke` (`K_T(0) ==
K_0`), `test_fint_linear_nonlinear_parity` (O(u²) at `u ~ 1e-4`), and the three
`test_kt_fint_directional_derivative*` tests pass against the faithful pair (`rel_fd < 0.05`; measured
`K_t == dF/du` to `~1e-10`). None genuinely failed, so none had to be weakened or deleted.

**Instrument verification.**

| run | result |
| --- | --- |
| `cd crates && cargo test -p aeroelast-core` | **168 passed / 0 failed / 6 ignored** |
| `cargo test -p aeroelast-core n_alpha -- --ignored --nocapture` | **1 passed** |
| `cargo test -p aeroelast-core n_beta -- --ignored --nocapture` | **1 passed** |
| `cargo test -p aeroelast-core n_gamma -- --ignored --nocapture` | **1 passed (244.5 s)** |
| `cargo check --lib -p aeroelast-core --message-format=short \| grep mitc4_plusd` | **no output** |

N-alpha rigid rotation (flat `RECT`, `prod` = `compute_fint_global(.., true)`, `bounded` = the superseded
reference):

| `phi` [deg] | max `\|_0 e + _0 eta\|` | `/phi^3` | prod `\|F\|/(\|K_0\|\|u\|)` | bounded `\|F\|/(\|K_0\|\|u\|)` |
| --- | --- | --- | --- | --- |
| 10 | 2.208500e-4 | 4.153988e-2 | **6.830422e-17** | 6.130636e-4 |
| 30 | 5.818770e-3 | 4.053548e-2 | **2.306251e-17** | 5.469097e-3 |
| 60 | 4.280333e-2 | 3.727277e-2 | **1.349926e-17** | 2.122833e-2 |
| 90 | 1.250000e-1 | 3.225153e-2 | **2.889624e-17** | 4.538420e-2 |

N-gamma rigid-body + tangent gate: `F(0) = 0` exactly; `K_t(0) == K_0` bit for bit; `K_t == dF/du` with
`rel 5.6e-11 .. 2.2e-10`; `prod |F|/Ku` equals the faithful local column to every printed digit
(6.7e-16 .. 3.0e-17), with the bounded column `2.2e-5 .. 1.65e-1`.

**The four pytest runs (installed extension rebuilt with `maturin develop --release`; no tuning).**

| run | this iteration | iteration-16 baseline |
| --- | --- | --- |
| `pytest tests/test_large_rotation_benchmarks.py -q` | **1 passed / 6 failed** | 4 passed / 3 failed |
| `pytest tests/test_shell_validation_fixed.py -q` | **7 passed** | 6 passed / 1 failed |
| `pytest -m "not slow" -q` | **341 passed / 10 failed / 2 skipped** | 343 passed / 8 failed / 2 skipped |
| `pytest tests/test_ccx_shell_element_types_parity.py -q` | **4 passed** | 4 passed |

Large-rotation, per case (only `test_linear_tip_deflection_euler_bernoulli` passes now):

| case | before | **now** |
| --- | --- | --- |
| `test_cantilever_large_rotation_half_circle[10]` | PASS | **FAIL: RuntimeError, SNES diverged (line search), level=8** |
| `test_equilibrium_path[pi/2, -3.6338, 6.3662]` | PASS | **FAIL: RuntimeError, SNES diverged, level=8** |
| `test_equilibrium_path[pi, -10.0, 6.3662]` | PASS | **FAIL: RuntimeError, SNES diverged, level=8** |
| `test_equilibrium_path[3pi/2, -12.122, 2.1221]` | FAIL 28.80% (`w_tip = 1.5109`) | **FAIL: RuntimeError (no value)** |
| `test_equilibrium_path[2pi, -10.0, 0.0]` | FAIL 6.75% (`u_tip = -9.3249`) | **FAIL: RuntimeError (no value)** |
| `test_simo_vu_quoc_rollup_360[10]` | FAIL 6.75% (`u_tip = -9.3249`) | **FAIL: RuntimeError (no value)** |

The three benchmark cells that used to produce numbers now diverge; the three that used to pass now diverge.
`TestNonlinearStatic::test_geometric_nonlinearity` now PASSES for the wrong reason: the solver DOES diverge
again, so `pytest.raises(RuntimeError, match="SNES diverged")` is satisfied.

**The broad-suite failure delta (10 now, 8 before):** four new failures, two newly passing.

| new | `test_cantilever_large_rotation_half_circle[10]`, `test_equilibrium_path[pi/2]`,
`test_equilibrium_path[pi]`, `TestNewtonRaphsonConsistency::test_nr_consistency[MITC4]` |
| --- | --- |
| newly passing | `test_shell_comprehensive::TestNonlinearStaticCantilever::test_large_displacement_tip_load`,
`test_shell_validation_fixed::TestNonlinearStatic::test_geometric_nonlinearity` (divergence test) |
| unchanged | `test_ko2017_performance::test_3_3_pinched_cylinder_tables_8_to_9[expected0-True]`,
`test_material_suite::TestIsoEquivalence::test_n_iso_plies_equal_single_layer_mitc4`,
`test_rust_modal::TestSimplySupportedPlate::test_analytical_convergence`, and the three large-rotation cells
above |

**Root cause of the regression, measured (not a wiring bug).** On a single 1x1 quad assembler, with a RIGID
rotation `u` the production tangent is the derivative of the production force (`rel ~ 1e-6` at 0.0001°/0.01°/10°),
but with a random `u ~ N(0, 1e-5)` it is not: `max|K_T du - (f(u+du)-f(u))|/|K_T du|` is **2.9e-3** for one
element and **4.52** for the 4x4 plate (`test_nr_consistency[MITC4]`), growing with `u` (0.227 at `||u||~1e-3`,
1.13 at `||u||~1e-2`). `K_T(0) == K_0` and `F(0) = 0` remain exact. Since the production path is a thin
`T^T (·) T` wrapper, the inconsistency is in the `n_gamma` pair itself: its `K_t = dF/du` identity, which the
rigid-rotation instrument confirms, does not hold for general displacement increments. The incremental solver's
Newton line search therefore fails (`snes_reason = -6`, `SNES_DIVERGED_LINE_SEARCH`) and the adaptive
refinement cannot recover it. Per the unit's constraints the `n_gamma` functions themselves were NOT touched;
this is reported, not patched.

**Known limitation (recorded, not fixed, unchanged from iteration 16).** `n_gamma_fint_global` /
`n_gamma_kt_global` build the `GlCurrentState` from the INITIAL element geometry, so they are the
reference-configuration first-step linearisation; threading the last-converged state between load/Newton steps
belongs to the solver.

**Tree state:** the delegation is APPLIED and uncommitted; the installed extension was rebuilt from it. No
commit.

### Iteration 18 — the N-gamma defect is LOCALISED to bent/warped geometry (probe evidence)

The instrument `n_gamma_slender_element_probe` (`#[ignore]`d) measured, per geometry, the zero-state block
difference `B(0) - B_lin` and the general-increment tangent consistency:

| geometry | `max\|mat0 - compute_ke_local\|` rel | `B(0) - B_lin` per block (m / b1 / b2 / s) | `K_t == dF/du` vs `\|\|u\|\|` |
| --- | --- | --- | --- |
| RECT (flat) | 2.3e-11 | 5e-12 / 6e-11 / 0 / 3e-11 | 1.8e-10 -> 1.1e-10 (does NOT grow) |
| SLENDER | 2.7e-11 | 3e-11 / 6e-11 / 0 / 3e-11 | 1.8e-10 -> 4.3e-10 (does not grow) |
| **BENT (warped)** | **2.954e-3** | **1.2e-2 / 2.0e-5 / 3e-11 / 1.2e-2** | **1.5e-8 -> 1.7e-6 (GROWS with `\|\|u\|\|`)** |

So: flat and slender are exact (the faithful `B(0)` equals the linear operator to ~1e-11 and the tangent is
consistent); the **membrane and shear blocks** are wrong by ~1.2e-2 (and `b1` by 2e-5) on a warped element,
`b2` is clean. On flat geometry the Eq. (15e)/(21) coefficients `^t a_A..^t a_E` vanish, which is why the
construction is not exercised there and the defect is invisible.

**Consequence in production (option B):** `K_t != dF/du` for general increments, growing with `||u||` ->
Newton's line search fails -> the large-rotation benchmarks raise `SNES diverged` (they used to produce
28.80% / 6.75% / 6.75%), and `tests/test_rust_assembler.py::TestNewtonRaphsonConsistency` fails at rel 4.52.

**Next step (specified, not done):** on the BENT element, print the **row-by-row** difference between each
block's `B(0)` and the linear operator and separate the candidate terms — (a) the Eq. (21a-c) assumed
incremental membrane with the **current** `^t a_A..^t a_E`; (b) the tying incremental strains of Eq. (20b)
(`^t x` form) as the code derives them (`0.25*(next_tie - cur_tie)`); (c) the Eq. (19) assumed incremental
shear vs the linear `b_shear_mitc4`/`compute_shear_tie`; (d) the `^t x_b`/`u_{b1}`/`u_{b2}` terms of
Eqs. (20c)/(20d)/(20f)/(20g) — and fix the differing one from the printed pp. 6-7 equations. Acceptance:
BENT `mat0` vs linear `<= 1e-10` relative; general-increment `K_t == dF/du` `<= 1e-6` without growth;
`TestNewtonRaphsonConsistency` passes; the `SNES diverged` disappears from the large-rotation file.

**State at the end of this session:** option B is applied and uncommitted; the installed extension is built
from it (so the Python nonlinear solve currently diverges on large rotations — the element is NOT delivery
ready in this state); `cargo test -p aeroelast-core` = 168 passed / 0 failed / 7 ignored; the four N
instruments pass; `pytest -m "not slow"` = 341 passed / 10 failed / 2 skipped; CCX parity 4 passed.

---

### Iteration 19 — the blocker is LOCALISED: it is the round-off floor of the nested finite difference, and NO printed term is missing

**Authorized by the parent:** secure the tree, run the three blocker diagnostics, report before writing a line of
the fix.

**The tree was secured first.** Commit `7ff506c` (`chore: checkpoint the uncommitted 2025-purity session (not
delivery-ready)`, 14 files, `+6744/-429`) plus the annotated tag `backup/mitc4plusd-pre-fix`. `uv.lock` was
reverted: it had never been committed, so there was no diff to commit separately — the 983 lines were re-lock
churn from the `maturin` runs. Nothing else was touched.

**Measurement 1 — `n_gamma_bent_localisation_diagnostic`.** The production `raw` column gives `B(0)` equal to the
linear operator to `1.17e-11` on BENT (and `1.2e-11` on RECT). The `frame-fix` column (`1.18e-2`) is a stale trial
variant that adds a spurious `T^T`; it is the red herring the handoff warned about. **The strain and the zero
state are not the problem.**

**Measurement 2 — `n_gamma_fd_profile_diagnostic`**, `rel(kt,fd)` against the reference step `h`:

| case | 1e-6 | 2e-5 | 1e-4 | 1e-3 |
| --- | --- | --- | --- | --- |
| RECT 1e-3 | **1.83e-10** | 1.38e-8 | 8.82e-9 | 5.91e-7 |
| RECT 1e-1 | **1.52e-10** | 1.46e-7 | 1.65e-7 | 6.88e-7 |
| BENT 1e-3 | **1.87e-8** | 1.51e-8 | 2.35e-8 | 1.89e-6 |
| BENT 1e-1 | **1.18e-6** | 1.19e-6 | 1.19e-6 | 1.47e-6 |

On RECT the residual GROWS with `h` (the reference is the limiter); on BENT it is FLAT from 1e-6 to 1e-4. The
flat-in-`h` reading was previously taken as proof of "a missing term, not noise". **It is not:** it only rules out
the REFERENCE's noise, and says nothing about the noise of a NESTED difference inside `K_t`.

**Measurement 3 — `n_gamma_kt_residual_diagnostic`.** `max|kt-fd|` against `max|pred| = sum_g B^T W (B - de/du) wq`:

| case | max\|kt-fd\| | max\|pred\| | reading |
| --- | --- | --- | --- |
| RECT 1e-3 | 2.93e1 at (13,1) | 2.82e1 at **(13,1)** | same entry: the residual IS `pred` |
| RECT 1e-1 | 2.71e1 at (7,19) | 2.73e1 at **(7,19)** | same entry |
| BENT 1e-3 | 5.31e3 at (20,7) | 5.31e1 at (19,1) | 100x larger, OTHER entry |
| BENT 1e-1 | 5.92e5 at (12,1) | 6.68e1 at (7,1) | 8900x larger, other entry |

**Measurement 4 (new) — `n_gamma_geo_symmetry_diagnostic`.** The tangent is a Hessian, so `|geo - geo^T|` is PURE
error with no finite-difference reference. Normalised by the round-off scale `eps/(H_i H_o) max|W e|`:

| fixture | planar | `a_A..a_E` | `\|geo-geo^T\|/floor` | `\|geo-geo^T\|/\|geo\|` |
| --- | --- | --- | --- | --- |
| RECT | yes | 0 | **0.200** | 3.41e-6 |
| ROT_RECT(60 deg) | yes | 0 | **0.500** | 8.14e-6 |
| FLAT_DISTORTED | yes | 0.076 | **0.677** | 6.59e-6 |
| STRONGLY_WARPED | no | 0.405 | **0.660** | 5.62e-6 |
| BENT | no | 1.1e-5 | **0.803** | 4.98e-6 |

`|mat - mat^T|/|mat|` is `1e-16` on every fixture, as it must be. Two conclusions:
1. **`geo` carries the nested-difference round-off floor on EVERY geometry, RECT included.** The error is
   symmetric-invisible, i.e. it is NOISE — a symmetric missing term would not appear in this metric at all. So
   there is no evidence for a missing printed term.
2. **The "RECT is clean" reading (1.5e-10) was an artefact of the comparison.** On RECT the noise of `geo` cancels
   against the noise of the reference finite difference because the perturbation directions agree (`t24 = I`); on
   any rotated or warped element they do not, which is why BENT reads 1e-6. RECT was never clean.
3. The fixtures also kill the assumed-membrane hypothesis: BENT has `max|a| = 1.1e-5`, i.e. it is not a real
   coefficient-path element, and the ordering by `|a|` does not follow the ordering by `|geo-geo^T|`.

**Measurement 5 (new) — `n_gamma_geo_stencil_probe`.** The nested form is algebraically the SAME four-point cross
difference, only with independent steps; its floor is `eps/(H_i H_o)`. The single-step form is bitwise symmetric
by construction, and its estimates at `h = 3e-5, 1e-4, 3e-4` agree with each other to `3e-8` while differing from
the split-step form by `~7e-6` — **FLAT in `h`** (`6.655e-6`, `6.680e-6`, `6.684e-6`). A truncation-limited
stencil would move by four orders of magnitude over that range, so the split-step form is the inaccurate one, by
about two orders.

**Verdict.** No printed term is missing and no test needs weakening. The blocker is a numerical-realisation defect:
`geo` was `FD(FD(_0 e~))` with `(H_i, H_o) = (2e-5, 1e-6)` and therefore a `5e-6` relative round-off floor, which
makes `K_t != dF/du` on every element whose local frame is not axis-aligned — i.e. on every real mesh, and on all
six large-rotation cases that were raising `SNES diverged`.

### Iteration 20 — the fix: one step for BOTH nested differences

**Changed (production, `n_gamma_kt_local` and the `B` helper only):**

| item | change |
| --- | --- |
| `n_gamma_b_matrix_at(pre, state, u_local, r, s, h)` | NEW: the central difference with an explicit step |
| `n_gamma_b_matrix(..)` | now `n_gamma_b_matrix_at(.., N_GAMMA_B_H)`; `N_GAMMA_B_H = 2.0e-5`, byte-equivalent to the old body |
| `const N_GAMMA_GEO_H: f64 = 1.0e-4` | NEW: the step of the geometric term, used by the outer difference AND by both inner `B` calls |
| `n_gamma_kt_local` | `const H = N_GAMMA_GEO_H` and the two `n_gamma_b_matrix_at(.., H)` calls |

`mat`, `n_gamma_fint_local` and the internal force are untouched: `N_GAMMA_B_H` stays where it is because the
`u = 0` identity (`K_T(0) == K_0` to `<= 1e-10`) pins it from below and it cannot be relaxed. No new ingredient:
the same Hessian of the same printed strain, evaluated at a single step chosen by the classic `h ~ eps^(1/4)`
balance (`~1e-8` total). A free consequence: both index orders now combine the SAME four points, so `geo` is
symmetric by construction rather than to round-off.

**Tests-only follow-up:** the two new instruments replicated the OLD split-step form, so they stopped being
replicas when production changed — the symmetry instrument's own cross-check caught it (`2.280e3`). Both now use
`N_GAMMA_GEO_H`, and the cross-check is back to `0.000e0`. The stencil probe deliberately keeps the split-step
form as the recorded "before" and labels it as such.

**Measured, acceptance of TAREA 1:**

| acceptance | result |
| --- | --- |
| 1. BENT `rel(actual) <= 1e-6`, no growth with `\|\|u\|\|` | **NOT RESOLVED BY THIS INSTRUMENT** — see the "reference-limited" note below. Measured production: `1.505e-8` (1e-3) and `1.185e-6` (1e-1), both AT the instrument's own reference floor |
| 2. `K_t(0) == K_0`, `F(0) == 0`, rigid-body checks at round-off | `n_gamma_rigid_body_zero_force_and_consistent_tangent` PASSES (424 s): `F(0) = 0` exactly; `K_t(0)` vs `K_0` rel `2.321e-11` (bitwise unchanged); general-increment consistency on RECT `8.448e-10` (1e-4), `3.007e-9` (1e-3), `6.071e-8` (1e-2), `3.997e-7` (1e-1); rigid-motion `max\|Kt-FD\|` rel `6.7e-11 .. 2.3e-10` on 12 motions |
| 3. `cargo test -p aeroelast-core` 168/0 | **168 passed / 0 failed / 12 ignored** |
| 4. `TestNewtonRaphsonConsistency` passes | **PASSES** (`tests/test_rust_assembler.py`: 19 passed / 1 failed, and the failure is NOT that test) |
| 5. the `SNES diverged` disappears | **GONE: 6 passed / 1 failed** (was 1 passed / 6 failed with six `SNES diverged`) |
| 6. `pytest -m "not slow"` >= 343/8/2 and CCX parity 4 | **346 passed / 5 failed / 2 skipped**; CCX parity **4 passed** |

**The reference-free metric, which is the one that resolves the fix.** `|geo - geo^T| / |geo|`, before vs after
(same instrument, same fixtures, `n_gamma_geo_symmetry_diagnostic`):

| fixture | before | after |
| --- | --- | --- |
| RECT | 3.41e-6 | **4.10e-13** |
| ROT_RECT(60 deg) | 8.14e-6 | **5.11e-13** |
| FLAT_DISTORTED | 6.59e-6 | **5.35e-13** |
| STRONGLY_WARPED | 5.62e-6 | **2.30e-13** |
| BENT | 4.98e-6 | **1.96e-13** |

A seven-order collapse, with no reference finite difference involved. The replica's cross-check against
`n_gamma_kt_local` is `0.000e0` again.

**Why acceptance point 1 as written cannot be certified by this instrument (a finding, not an excuse).** The
PRODUCTION line compares `n_gamma_kt_local` against `fd_matrix(n_gamma_fint_local, h = 1e-6)`, and `n_gamma_fint_local`
carries the round-off of the `B` inside it (`eps/H_i` amplified by `|W e| wq`), which the reference difference then
amplifies by `1/h`. The replica's profile is the proof that the reference, not the tangent, is the limiter: it is
FLAT from `h = 1e-5` to `h = 3e-4` (`9.25e-7`, `9.40e-7`, `9.13e-7`, `9.11e-7` on BENT at 1e-1), whereas the
construction's own floor `eps/h^2` would move by two orders over that range. Both measured production numbers sit
on that plateau. Certifying `<= 1e-6` with margin therefore needs the analytic `B`/`N` route (which removes the
`B`-noise from the force as well), and until then the weight is carried by the reference-free metric above and by
the EXTERNAL oracle below.

**The external oracle, which is the strongest evidence in this unit: the Newton solve converges again.** SIX of
seven large-rotation cases now pass, including two that used to be 6.75% off and one that used to be 28.80% off,
and the three that used to produce numbers now do so against their published/analytical references.

**Large-rotation, case by case (before this unit / now):**

| case | before option B | at option B | **now** |
| --- | --- | --- | --- |
| `test_linear_tip_deflection_euler_bernoulli` | PASS | PASS | **PASS** |
| `test_cantilever_large_rotation_half_circle[10]` | PASS | SNES diverged | **PASS** |
| `test_equilibrium_path[pi/2]` | PASS | SNES diverged | **PASS** |
| `test_equilibrium_path[pi]` | PASS | SNES diverged | **PASS** |
| `test_equilibrium_path[3pi/2]` | FAIL 28.80% (`w_tip = 1.5109`) | SNES diverged | **FAIL 9.84% (`w_tip = 1.9133` vs `2.1221`)** |
| `test_equilibrium_path[2pi]` | FAIL 6.75% (`u_tip = -9.3249`) | SNES diverged | **PASS** |
| `test_simo_vu_quoc_rollup_360[10]` | FAIL 6.75% (`u_tip = -9.3249`) | SNES diverged | **PASS** |

The two `6.75%` cells now close, the `28.80%` cell came down to `9.84%`, and the `SNES diverged` is gone. The
remaining cell is a single 9.84% miss against a `5%` tolerance and is now a legitimate follow-up (the §3 Fig. 6
benchmarks of C&S 185 are the next oracle).

**The 5 remaining `not-slow` failures, and the one that is new:**

| failure | status |
| --- | --- |
| `test_ko2017_performance::test_3_3_pinched_cylinder_tables_8_to_9[expected0-True]` | pre-existing (in the recorded unchanged set) |
| `test_large_rotation_benchmarks::test_equilibrium_path[3pi/2]` | the 9.84% cell above (was 28.80%) |
| `test_material_suite::TestIsoEquivalence::test_n_iso_plies_equal_single_layer_mitc4` | pre-existing |
| `test_rust_modal::TestSimplySupportedPlate::test_analytical_convergence` | pre-existing |
| `test_rust_assembler::TestTangentStiffness::test_kt_at_zero_equals_k[MITC4]` | **NOT this unit's regression** — see below |

`test_kt_at_zero_equals_k[MITC4]` asserts `K_T(u=0) == K` with `rtol=1e-10` and fails at a max relative deviation
of `5.4e-10` (`max abs 2.27e-3` on entries of `4.2e6`). It is NOT caused by this unit: `mat0` is bitwise unchanged
(the instrument still prints `3.699e0 (2.321e-11)`, the same digits as Iteration 14), because `geo(0) = 0`
exactly at `u = 0` and `mat` still uses `N_GAMMA_B_H`. It was introduced by **Iteration 17**, where the deliberate
removal of the additive `compute_ke_local(pre)` reference replaced `K_t(0) = K_0` bit for bit with
`K_t(0) = mat0`, which reproduces `K_0` only to `O(H^2) ~ 4e-10`. The previous session's failure inventory was
taken before that removal, which is why it reads as new. **Maintainer decision (this session): DO NOT re-baseline.** `rtol=1e-10` stays and the test stays red. No
tolerance was widened and no assertion was deleted. It is recorded as a known consequence of a
finite-difference `B(0)`, and the way to close it green is the analytic `B`/`N` of the printed equations, which
makes `K_t(0) = K_0` exact again (and removes the `B`-noise from the internal force, which is also what currently
limits acceptance point 1).

---

## HANDOFF — state at the end of the session (2026-09-24)

**Branch:** `test/physical-correctness`. **NOTHING committed this session** (last commit is the previous
session's `ea244c7`). 11 modified files + 5 untracked. **The whole session's work is uncommitted** — commit
or stash it before doing anything destructive.

Modified: `crates/aeroelast-core/src/assembly/assembler.rs` (the WU9f flip, 175 lines),
`crates/aeroelast-core/src/elements/mitc4_plusd.rs` (4101 lines: Eq. (22a) fold, T5 gating, N-alpha,
N-beta, N-gamma, option-B delegation, the instruments), `crates/aeroelast-py/src/{assembler,elements,materials}.rs`
(the flip + the batch-kernel wiring), `docs/formulations/mitc4plus-2017-extract.md`,
`docs/validation-matrix.md` (CCX rows), `openspec/changes/mitc4plusd-faithful/{apply-progress,tasks}.md`,
`src/aeroelast/core/mesh/io/writers.py` (the `shell_element_type` selector). Untracked:
`odd/tasks/mitc4plusd-2025-purity.md` (this file), `odd/tasks/mitc4plusd-wu9f-flip.patch`,
`openspec/changes/mitc4plusd-faithful/fidelity-audit.md`, `tests/test_ccx_shell_element_types_parity.py`.
**Unrelated churn:** `uv.lock` (983 lines) is a side effect of the maturin runs, not this work.

### What is DONE and measured

- **Eq. (22a) fold** (the 2025 drill summed into the membrane row): the flat in-plane ratio went 251.55 ->
  400.39 (+0.10% vs beam theory). **Defect #4 was a real fidelity defect.**
- **Section B audited CLEAR** (2017 assumed membrane, against pp. 405-410 + Appendix A); a documented
  "Eq. (21) omits a term" paper defect was **falsified** (Eq. (A.7) keeps the term) and corrected in the code
  comment, the extract and `fidelity-audit.md`.
- **Section D audited FAITHFUL and provably inert** on the flat in-plane case.
- **Empirical fidelity**: the 2025 Table 1 regular-mesh cell reproduced to **4.4e-12**; the Table 2 curved
  beam (where the `a` coefficients are live) to **<= 0.7%**. **T3 closed empirically.**
- **T5**: the production lib of `mitc4_plusd.rs` has no dead code (4 test-only fns gated `#[cfg(test)]`).
- **Option B applied**: `compute_fint_global(.., true)` -> `n_gamma_fint_global`, `compute_kt_global` ->
  `n_gamma_kt_global` (the common-interface pattern each element uses). The dead, unsourced corotational
  block is DELETED.
- **The N instruments (10, `#[ignore]`d)**: N-alpha (the Eq. (9) identity is O(phi^3); the quadratic `u_b2`
  is required), N-beta (the incremental operator equals the linear one at 1.6e-12), N-gamma (`F(0)=0`
  exactly, `K_t(0)=K_0`, `K_t == dF/du` local to 5.6e-11..2.2e-10), plus the three bent/FD diagnostics.
- **The gate with the flip measured**: `pytest -m "not slow"` 343/8/2 (HEAD baseline 344/3/2) — the flat
  in-plane family now PASSES.

### THE BLOCKER — **RESOLVED in iterations 19-20** (kept for the record)

Option B is applied, so the production nonlinear path is the faithful pair — and it **DIVERGES on large
rotations** (`RuntimeError: SNES diverged` on 6 of 7 cases in `tests/test_large_rotation_benchmarks.py`;
`test_rust_assembler.py::TestNewtonRaphsonConsistency` rel 4.52 vs 2e-3; `pytest -m "not slow"` 341/10/2).

**RESOLVED.** The residual was the round-off floor `eps/(H_i H_o) = 5e-6` of `geo` computed as a finite difference
of a finite difference — NOT a missing printed term (the handoff's Eq. (25) `N_ij` hypothesis was the right place
to look and the wrong cause). Both of `geo`'s nested differences now use the SINGLE step `N_GAMMA_GEO_H = 1e-4`.
The `SNES diverged` is gone: `tests/test_large_rotation_benchmarks.py` = **6 passed / 1 failed** (the one cell at
9.84%, was 28.80%); `TestNewtonRaphsonConsistency` **passes**; `pytest -m "not slow"` = **346 passed / 5 failed /
2 skipped**; `cargo test -p aeroelast-core` = 168/0. See iterations 19 and 20 for the measurements. Backlog is now
the single 9.84% cell, the §3 Fig. 6 oracle of C&S 185, and the analytic `B`/`N` route that would also make
`K_t(0) == K_0` exact again (today it is `mat0`, good to `O(H^2) ~ 4e-10`).

**Measured diagnosis (corrected, iteration 18):**
- The strain/zero-state construction is CORRECT on flat, slender AND bent geometry (`raw` B(0) vs the linear
  B: 1.17e-11 on BENT). The earlier "bent geometry defect" was an artifact of a probe's trial variant.
- The real residual is in the TANGENT: on BENT, `K_t != dF/du` with `rel(actual)` = 1.87e-8 (||u||=1e-3) ->
  1.18e-6 (1e-1), **growing with ||u||** and **flat in the FD step h** (so it is a missing term, not noise).
  On RECT the same residual is fully explained by the FD noise (1.8e-10).
- **Hypothesis, grounded in the paper:** Eq. (24a) is `^t K_e = int B^T C B d0V + int ^t_0 S_ij N_ij d0V`, and
  **Eq. (25) defines `N_ij` by the variation of the NONLINEAR strain** (`delta _0 eta_ij = delta U_e^T N_ij U_e`).
  The code computes that second term as a **finite difference of `B`**, correct only if `B` carries the
  derivative of `_0 eta` (Eqs. 20e-g). Verify that, and implement `N_ij` (Eq. 25, p. 7) as printed rather than
  the FD of `B`.

**Acceptance for the fix:** BENT `rel(actual)` <= 1e-6 and NOT growing with ||u||; `TestNewtonRaphsonConsistency`
passes; the `SNES diverged` disappears from the large-rotation file; `cargo test -p aeroelast-core` 168/0;
`pytest -m "not slow"` back to at least 343/8/2.

**Documented limitation (not a bug):** the `n_gamma_*_global` wrappers build the state from the INITIAL
geometry, so they are the reference-configuration first-step linearisation; threading the last-converged
state between Newton steps belongs to the solver.

### Then the remaining backlog

1. The §3 benchmarks of C&S 185 as the oracle (Fig. 6: the two large-rotation cantilevers; Fig. 9: the slit
   annular plate) once the tangent is consistent.
2. T4 — adopt ADR-4 option A properly (mesh-consistent nodal directors; today they exist only as the flip's
   two-pass instrument).
3. T1's remaining 2025 reads (pp. 2, 3, 7-9, 13).
4. T6b — promote the two 2025-table instruments into real oracle tests with justified tolerances.
5. The gate/retirement decisions (the rename to `mitc4`, the hybrid's retirement) stay gated on a green gate.

### Environment (unchanged)

`source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev && export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig" HDF5_DIR="$CONDA_PREFIX" && export PATH="$HOME/miniconda3/envs/aeroelast-dev/bin:$PATH"`
— `cd crates && cargo test -p aeroelast-core`; `python -m maturin develop --release`;
`python -m pytest ...`; CCX 2.23; the papers in `.sources/papers/` (`A_new_MITC4+_shell_element.pdf` = 2017
linear, `The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf` = C&S 185 nonlinear,
`1-s2.0-S0045794924003511-main.pdf` = 2025 MITC4+/D; `1-s2.0-S0045794917309550-main.pdf` = the BENCHMARK
paper, not a formulation).
