# Shell element formulations: code-to-equation reference

This is the production reference for the shell and plane elements in
`crates/aeroelast-core/src/elements/`: what each element implements, which published
element it is, and which equation of which paper each part comes from. Equations
were read from the PDFs held in `.sources/papers/`; where an equation could not be
read, this document says so instead of reconstructing it. The full bibliography,
including DOIs and the provenance of every held copy, is `docs/references.md`.

## 1. Scope and conventions

### 1.1 Implemented elements

| Element | File | Nodes / DOF | Paper |
| --- | --- | --- | --- |
| MITC4+/D | `elements/mitc4.rs` | 4 × 6 = 24 | Ko, Lee & Bathe (2017), C&S 182:404–418, with the six-DOF drill-membrane strain of Ko, Bathe & Zhang (2025), C&S 308:107622 |
| MITC3+ | `elements/mitc3.rs` | 3 × 6 = 18 (+2 bubble) | Lee, Lee & Bathe (2014), C&S 171:21–34 |
| Strain-smoothed MITC3+ | `elements/mitc3.rs` + `elements/smoothing.rs` | union layout, six nodes | Lee & Lee (2019), C&S 206:181–191 |
| Plane quadrilaterals and triangles | `elements/quad.rs` (integration via `elements/reference.rs`) | 4/8/9 nodes | **no published formulation cited** — see §4.4 |

Documented in detail here: §2 (MITC4+/D) and §3 (MITC3+). The plane elements are
implemented and tested but carry no literature citation, so §4.4 states that gap
rather than inventing a source for them.

### 1.2 DOF ordering

Both elements use the same per-node ordering, six global Cartesian components per
node, packed node-major:

```text
(u, v, w, theta_x, theta_y, theta_z)     6 per node, global Cartesian components
```

- MITC4+: `24` DOF, node `i` occupies indices `6i .. 6i+5`
  (`mitc4.rs`: `u_idx = 6*i`, `v_idx = 6*i+1`, `thz_idx = 6*i+5`).
- MITC3+: `18` nodal DOF in the same layout, plus two internal bubble rotation DOF
  at indices `18`, `19` in the extended `20`-DOF space, condensed out of the
  returned `18×18` matrix.

### 1.3 Rotation convention

`theta` is the **physical rotation vector**, not a director-rotation pair. The
offset displacement interpolates as

```text
u(r,s,t) = sum_i h_i u_i + (t/2) sum_i a_i h_i (theta_i x V_in)                    (3a)
```

Ko, Y., Bathe, K.-J., Zhang, X., "Continuum mechanics-based shell elements with six
degrees of freedom at each node — the MITC4/D and MITC4+/D elements", *Computers
and Structures* 308:107622, 2025 — `1-s2.0-S0045794924003511-main.pdf`, journal
p. 2. Verified verbatim in the extraction; the same paper defines
`theta_i = theta_ix i_x + theta_iy i_y + theta_iz i_z` as Eq. (3c).

Linearizing `theta x V_in` with `V_in = e3` gives `V3 = e3 + theta_y e1 - theta_x e2`,
hence the Cartesian transverse-shear convention

```text
gamma_xz = w_,x + theta_y
gamma_yz = w_,y - theta_x
```

consequences of Eq. (3a). MITC3 was brought onto this convention by commit
`d6f37fb` ("fix(mitc3): restore the physical rotation sign in the transverse
shear"), which restored the pre-`b136ce5` form; commit `b136ce5` had inverted it.
The MITC4 shear operator (`b_gamma_mitc4`) already used this convention before that
fix. The derivation comment now lives in `mitc3.rs` at the covariant shear assembly
("From Ko, Bathe & Zhang 2025 (MITC4/D), Eq. (3a) ...").

This reference document does **not** claim that the MITC3+ 2014 paper's own
rotation parameters `alpha`, `beta` carry the same sign; see §4.5.

### 1.4 Local frame, and what "local" means

Both elements precompute, at construction time, an orthonormal triad
`(e1, e2, e3)` for the element (`t3` in `Mitc4Precomputed` / `Mitc3Precomputed`),
together with the Cartesian node coordinates projected onto it
(`local_coords`). All strain operators referenced below are evaluated **in that
local frame**:

- `compute_ke_local` assembles the local `24×24` (resp. `18×18`) stiffness; the
  constitutive matrices are rotated into the local frame once
  (`rotate_constitutive_to_local` in `mitc3.rs`), and
  `compute_ke_global = t3 * k_local * t3^T`.
- The covariant quantities (`g_r`, `g_s`, `x_d`, `m_r`, `m_s`) used by the assumed
  strain fields are 3D vectors, but the B-matrix rows that multiply the DOF vector
  are expressed against the **local Cartesian** basis. The covariant-to-local
  mappings are `covariant_to_local_mapping(j_loc)` in `mitc4.rs` and
  `j_inv` in `mitc3.rs` (`gamma_xz = J^{-1}[0,.] . e_cov`, `gamma_yz = J^{-1}[1,.] . e_cov`).
- Therefore "local" for a strain operator means the element-local Cartesian frame;
  the covariant sampling is an intermediate purely for the assumed-field
  construction.

---

## 2. MITC4+/D (quadrilateral, 24 DOF) — the production shell element

The element is the MITC4+ of Ko, Lee & Bathe (2017), *Computers and Structures*
182:404–418, doi:10.1016/j.compstruc.2016.11.004 — with the sixth, drilling degree
of freedom supplied by the drill-membrane strain of Ko, Bathe & Zhang (2025),
*Computers and Structures* 308:107622, doi:10.1016/j.compstruc.2024.107622. That
combination is what the 2025 paper calls **MITC4+/D**. It is the only shell
quadrilateral in production: the earlier MITC4+ hybrid that occupied this file was
retired in full, so nothing below describes a penalty, a reduced-integration split
or a bubble that the element does not actually contain.

The element introduces **no numerical factor**: no penalty, no factor on the
transverse shear stiffness, no selective reduced integration. Transverse shear is
handled by the assumed MITC4 field alone, and the drilling degree of freedom by the
2025 assumed drill-membrane strain rather than by a penalty.

### 2.1 As implemented

- **Nodes and integration.** 4 nodes × 6 DOF. `N_GAUSS = 4`; `GAUSS_XI/ETA = ±1/√3`
  with weights `1.0`, i.e. the 2×2 Gauss–Legendre rule on `[-1,1]²`, and `ζ`
  integrated with the same rule: the paper's "we use 2 × 2 × 2 Gauss integration"
  (journal p. 410). The through-thickness rule is exact for the bending moments the
  element forms, because `ζ²` is integrated exactly by the two-point rule.
- **Membrane field.** The assumed covariant membrane field of Ko, Lee & Bathe
  (2017): the mid-surface metric assumed through the five tying points (A)–(E) of
  their Fig. 4, blended by the coefficients of their Eq. (24) and the block printed
  after Eq. (27c). §2.3 quotes the equations; §2.4 covers the transverse shear.
- **Transverse shear.** The rotation-based MITC4 (Dvorkin & Bathe 1984)
  interpolation, sampled on the element edges — §2.4.
- **Drilling degree of freedom.** The drill-membrane strain of Ko, Bathe & Zhang
  (2025), Eq. (26), built with the operator of their Eq. (18). It is a *strain*
  contribution, not a penalty: there is no stiffness coefficient to tune and no
  scale factor in the code.
- **Finite rotations.** The kinematics are the Green–Lagrange strain of the
  offset displacement field of Eq. (3a), linearised in the increment around the
  current director state (`gl_strain_increment` and the `gl_*` family in
  `mitc4.rs`). This is what makes the geometric stiffness
  (`compute_k_sigma_global`) consistent with the internal forces rather than a
  separate approximation.
- **No hourglass control, and none needed.** `compute_ke_local` assembles the
  membrane, bending and transverse-shear blocks of the assumed fields above; there
  is no hourglass term, and no hourglass stiffness is precomputed.

### 2.2 Code → equation table

| Code (`mitc4.rs`) | Paper and equation | Verified |
| --- | --- | --- |
| `gl_strain_increment`, `gl_assumed_mid_metric` | Ko, Lee & Bathe (2017), Eqs. (9) and (21a)–(21c) | yes, term by term |
| `gl_tying_metrics` — the five tying points (A)–(E) | Ko et al. (2017), Fig. 4 | structurally: the five parametric positions in §2.3 |
| `compute_membrane_coefficients_2017` | Ko et al. (2017), Eq. (24) and the block after Eq. (27c) | yes, term by term (§2.3) |
| `b_membrane_2017` rows `rr`, `ss`, `rs` | Ko et al. (2017), Eqs. (27a)–(27c) | yes, term by term (§2.3) |
| `b_shear_mitc4` | Dvorkin & Bathe (1984), Eq. (3), reproduced unnumbered in Ko et al. (2017) §2 | partially — §2.4 states exactly what could not be read |
| `b_drill_membrane_2025`, `drill_midside_shape_derivatives` | Ko, Bathe & Zhang (2025), Eqs. (18) and (26) | yes: the operator and the strain it multiplies |
| `compute_ke_local`, `compute_ke_global` | Ko, Lee & Bathe (2017), p. 410 — 2 × 2 × 2 Gauss, no numerical factor | yes |
| `compute_k_sigma_global` | Ko, Lee & Bathe (2017) — the finite-rotation formulation whose geometric stiffness is consistent with the Green–Lagrange kinematics above | yes, by the Newton tangent tests |

Where the table says "partially", the missing part is named in the section it points
to and is never presented as verified.

### 2.3 Membrane field — Eqs. (27a), (27b), (27c)

Quoted from journal p. 409 (Eq. (27a)) and p. 410 (Eqs. (27b), (27c)), with
`r ≡ xi`, `s ≡ eta`:

```text
e~rr = 1/2 (1 - 2aA + s + 2aA s^2) err(A)
     + 1/2 (1 - 2aB - s + 2aB s^2) err(B)
     + aC (-1 + s^2) ess(C) + aD (-1 + s^2) ess(D) + aE (-1 + s^2) ers(E)      (27a)

e~ss = aA (-1 + r^2) err(A) + aB (-1 + r^2) err(B)
     + 1/2 (1 - 2aC + r + 2aC r^2) ess(C)
     + 1/2 (1 - 2aD - r + 2aD r^2) ess(D)
     + aE (-1 + r^2) ers(E)                                                    (27b)

e~rs = 1/4 (r + 4aA r s) err(A) + 1/4 (-r + 4aB r s) err(B)
     + 1/4 (s + 4aC r s) ess(C) + 1/4 (-s + 4aD r s) ess(D)
     + (1 + aE r s) ers(E)                                                     (27c)
```

with (printed immediately after Eq. (27c); the `c_r`, `c_s`, `d` definitions are
Eq. (24)):

```text
aA = cr(cr - 1)/(2d)   aB = cr(cr + 1)/(2d)   aC = cs(cs - 1)/(2d)
aD = cs(cs + 1)/(2d)   aE = -2 cr cs / d                                   (after 27c)

cr = xd . mr      cs = xd . ms      d = cr^2 + cs^2 - 1                          (24)
```

`aA..aE` are **not** Eq. (24): Eq. (24) defines `B1..B5` (`cr²/d`, `cs²/d`,
`2 cr cs/d`, `-cr/d`, `-cs/d`) and `cr, cs, d`, used in the intermediate Eq. (25).
The `aA..aE` block is printed separately, after Eq. (27c), in the "with" clause.
Appendix A is titled "Derivation of the constants in Eq. (24)".

**Where this lives in the code.** The implemented rows are in `b_membrane_2017` and
the assumed mid-metric in `gl_assumed_mid_metric` (`elements/mitc4.rs`); the
equations above are the 2017 paper's and were not changed when the element was
rewritten to the printed formulation.

One structural detail: the code stacks the covariant rows as
`B_cov = [B_rr; B_ss; 2*B_rs]` (`b_cov[(2,j)] = 2.0*b_rs[j]`). The `2*B_rs` is the
engineering-shear factor applied at the covariant→local mapping, not a change to
Eq. (27c).

### 2.4 Transverse shear — Dvorkin & Bathe 1984

Reference: Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node
shell element for general non-linear analysis", *Engineering Computations*
1:77–88, 1984 — the page footer of the held scan prints
`Eng. Comput., 1984, Vol. 1, March 77`, and the last page prints 88; no DOI is
printed — `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`.

**Status of the verification: partial. Report honestly.**

- The **displayed equations of the original 1984 paper could not be read.** The
  scan has no text layer over its equations; `pdftotext -layout` and
  `pdftotext -raw` both omit them. Page 2 reads, verbatim, "The choice of the
  interpolation for the transverse shear strain components is the key assumption
  in our element formulation ... For our element we use (see Figure 2):" and then
  the displayed equation is absent from the extraction. The equation number of
  that equation therefore cannot be reported, and is not invented here.
- What **is** readable in the 1984 text (p. 3): "the interpolation employed for
  the transverse shear strains shows that [gamma_rt] is constant with r1 and in
  general discontinuous at [element boundaries] ... similarly [gamma_st] is
  constant with r2"; the interpolations are performed "in a convected coordinate
  system"; and the tying points are the element edges (Figure 2).
- The same field is reproduced in **Ko et al. 2017, §2** ("The transverse shear
  strain field is based on assuming constant covariant transverse shear strain
  conditions along the edges, see Ref. [2]", i.e. Dvorkin & Bathe), with the
  equations **unnumbered** in the extraction:

  ```text
  e~rt = 1/2 (1 + s) ert(A) + 1/2 (1 - s) ert(B)
  e~st = 1/2 (1 + r) est(C) + 1/2 (1 - r) est(D)
  ```

  Tying points A–D lie on the element edges (Ko et al. 2017, Fig. 2: "Tying
  positions (A)–(D) for the assumed transverse shear strain field of the MITC4
  shell element. The constant transverse shear strain conditions are imposed
  along its edges.").

**Correspondence to `b_gamma_mitc4`.** The code matches this field structurally:
it samples four edges (ordered `4-1`, `1-2`, `2-3`, `3-4`) into a `4×12` `G`
matrix and multiplies by an `Ms` matrix carrying exactly the `(1±xi)`, `(1±eta)`
linear weights, then scales each row by `r1/(8 |g_r x g_s|)` and `r2/(8 |g_r x g_s|)`
with `r1 = |C + xi B|`, `r2 = |A + eta B|` built from the element's edge-difference
vectors. **What could not be verified:** whether the `(1/8)·|g_r x g_s|^{-1}`
scaling and the angle pair `alpha = atan2(Ay, Ax)`, `beta = pi/2 - atan2(Cx, Cy)`
reproduce the original paper's normalisation and its alpha/beta coordinate system,
because the source equations are unreadable and the reproduction in Ko et al. 2017
carries neither the angles nor the scaling. This is an open item, not a claim.

### 2.5 The drilling degree of freedom (2025)

The sixth DOF is not a penalty. Ko, Bathe & Zhang (2025) introduce a
drill-membrane strain, Eq. (26) of that paper, built with the midside operator of
their Eq. (18), and add it to the in-plane strain field. The element consumes it as
a strain contribution, so its stiffness follows from the constitutive law like any
other membrane component.

The practical consequences, and the reason the section exists at all: there is **no
stiffness coefficient, no scale factor and no user-tunable drilling parameter** in
the element. A reader arriving from the older documentation — this file used to
describe a `k_drill = 0.15 · E · h² · drilling_scale` penalty attributed to "Hughes
& Brezzi" — will not find any of it, because the element no longer contains it. The
Hughes & Brezzi reference itself does exist, as `docs/references.md` §1 records; it
is simply no longer the source of this element's drilling stiffness.

### 2.6 What this element does not contain

Stated explicitly because each of these was part of the documentation before the
element was rewritten, and none of them is part of the element now:

- no selective reduced integration: the in-plane shear is integrated by the same
  2 × 2 rule as everything else, and there is no centre-point substitution. The
  retired scheme was Hughes, Taylor & Kanoknukulchai 1977 (References item 2, kept
  as a retired source for exactly this reason);
- no factor on the transverse shear stiffness: the assumed MITC4 field is the only
  mechanism, so the classical `5/6` reference value appears in this document only
  where a benchmark is compared against a classical solution;
- no drilling penalty and no scale factor (§2.5);
- no rotation bubble and no static condensation;
- no hourglass stiffness.

The retirement test `test_retirement_removes_hybrid_deviation_surfaces` in
`elements/mitc4.rs` asserts that these names are absent from the module, so this
section cannot silently become false again.

## 3. MITC3+ (triangle, 18 DOF)

Reference: Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element and its
performance", *Computers and Structures* 138:12–23, 2014,
doi:10.1016/j.compstruc.2014.02.005 —
`.sources/papers/The_MITC3+_shell_element_and_its_performance.pdf`.

### 3.1 As implemented

- **Topology.** 3 nodes × 6 DOF, plus **one bubble node carrying 2 rotation DOF
  only** (`a4`, `b4`), positioned on the flat surface of the corner nodes.
  Extended space `20`, statically condensed to `18`. Integration: Hammer degree-2
  triangle rule, `N_GAUSS = 3`, weights `1/3`.
- **Enriched shape functions.** Code `shape_functions` / `bubble`:

  ```rust
  h1 = 1 - r - s;  h2 = r;  h3 = s;
  f4 = 27 * r * s * (1 - r - s);
  ```

  Paper Eqs. (1) and (7): `h1 = 1 - r - s; h2 = r; h3 = s` (1) and
  `f1 = h1 - f4; f2 = h2 - f4; f3 = h3 - f4; f4 = 27 r s (1 - r - s)` (7). The code
  uses the same `f4` magnitude (`27 * r * s * (1 - r - s)`) and the same `h_i`, but
  a **different enrichment coefficient**: `enriched_shape` returns
  `fi = hi - f4/3` for `i = 1, 2, 3`, not `hi - f4`.

  **Extraction caveat, quoted as it stands.** The `-layout` extraction renders
  Eq. (7) as

  ```text
  f1 ¼ h1  f4 ; f 2 ¼ h2  f4 ; f 3 ¼ h3  f4 ; f 4 ¼ 27rsð1  r  sÞ:            ð7Þ
  ```

  The operator glyph between `h_i` and `f4` is lost, exactly as the two minus signs
  in `(1 - r - s)` are lost in the same line, so the extraction cannot distinguish
  `hi - f4` from `hi - f4/3` with certainty; the visible character count is
  consistent with a single operator. This is recorded as an unresolved transcription
  question: either the paper prints `hi - f4` and the code deviates, or it prints
  `hi - f4/3` and the extraction is misleading. The code's `f4/3` is the coefficient
  that preserves partition of unity (`f1 + f2 + f3 + f4 = h1 + h2 + h3 = 1`), which
  favours the latter, but that is an inference and not evidence.
- **Membrane field.** Constant strain triangle, §3.3. The paper's displacement
  interpolation is Eq. (2); the covariant strain definition is Eq. (3)
  (`e_ij = 1/2 (g_i . u_,j + g_j . u_,i)`).
- **Assumed transverse shear.** Tying points A–F:

  ```rust
  A = (1/6, 2/3)   B = (2/3, 1/6)   C = (1/6, 1/6)
  D = (1/3 + d, 1/3 - 2d)   E = (1/3 - 2d, 1/3 + d)   F = (1/3 + d, 1/3 + d)
  d = D_PARAM = 1.0e-4
  ```

  The paper chooses `d = 1/10,000` and states "we choose and always use the
  distance d = 1/10,000" (journal p. 15), after a strain-energy study in which the
  element becomes "rapidly more flexible in the in-plane twisting mode" as `d -> 0`
  (Fig. 6; Tables 2, 3). The code's `1.0e-4` matches. Tying points (D), (E), (F)
  "are positioned on the three internal lines from the barycenter to the centers of
  the edges" (journal p. 15).
- **Static condensation.** `compute_ke_local` extracts the nodal and bubble blocks
  and returns `k_bs_cond = k_uu + k_mb_sym - k_uq_full * inv_qq * k_qu_full` (18×18);
  the full local matrix is `km + k_drill_total + k_bs_cond`.
- **Drilling.** A single-row penalty on `theta_z`, assembled over the 3-point rule
  with `k_drill = E * h^2 * 0.15 * drilling_scale` — an expression that carries no
  source comment here and is not attributed to a work. Note that this is MITC3+'s own
  construction and is unrelated to the MITC4+/D drilling of §2.5, which is a strain
  contribution with no penalty coefficient at all.

### 3.2 Code → equation table

| Code | Paper | Verified |
| --- | --- | --- |
| `b_gamma_ext`, constant part `b_ert_const` | Eq. (15), second line | yes, term by term |
| `b_gamma_ext`, constant part `b_est_const` | Eq. (15), third line | yes, term by term |
| `b_gamma_ext`, linear part via `b_c_hat` | Eq. (16) | yes, term by term |
| total assumed shear `b_ert + b_ert_linear`, `b_est + b_est_linear` | Eq. (17) | yes, term by term |
| `b_membrane` | CST operator; the paper prints the shape functions as Eq. (1) and the covariant strain definition as Eq. (3), but no numbered Cartesian CST operator was located in the extracted pages | definition only |
| shear sign convention in `b_gamma_ext`/`eval_covariant_shear_ext` | Ko, Bathe & Zhang 2025 Eq. (3a) | yes |

**Correction to a commonly repeated claim.** The constant part is **Eq. (15)**, not
Eq. (16). Eq. (16) is the *linear* part. The code's own comments
(`// Constant part (Eq. 15)`, `// Linear part (Eq. 16)`) are correct. Verified
against the paper directly.

### 3.3 Constant part — Eqs. (15) and (17)

Quoted from journal p. 15 (Eq. (15) is printed as two lines under one number):

```text
e^const_rt = e^rt = 2/3 ( ert(B) - 1/2 est(B) ) + 1/3 ( ert(C) + est(C) )         (15)
e^const_st = e^st = 2/3 ( est(A) - 1/2 ert(A) ) + 1/3 ( ert(C) + est(C) )         (15)
```

and the total assumed field, Eq. (17):

```text
e^rt = e^const_rt + e^linear_rt = 2/3 ( ert(B) - 1/2 est(B) ) + 1/3 ( ert(C) + est(C) )
       + (1/3) c^ (3s - 1)                                                       (17)
e^st = e^const_st + e^linear_st = 2/3 ( est(A) - 1/2 ert(A) ) + 1/3 ( ert(C) + est(C) )
       + (1/3) c^ (1 - 3r)                                                       (17)
```

Note the asymmetry, which is easy to get wrong: `e^const_rt` samples at **(B)** and
`e^const_st` samples at **(A)**. Code:

```rust
// mitc3.rs, b_gamma_ext — comment: "Constant part (Eq. 15)"
let b_ert_const = (2.0/3.0) * (b_ert_b - 0.5*b_est_b) + (1.0/3.0) * (b_ert_c + b_est_c);
let b_est_const = (2.0/3.0) * (b_est_a - 0.5*b_ert_a) + (1.0/3.0) * (b_ert_c + b_est_c);
```

matches Eq. (15) exactly, with `a <-> (A)`, `b <-> (B)`, `c <-> (C)`. Eq. (15) is in
turn derived from Eq. (14), which assumes the constant covariant edge shear strains
along the three edge directions.

### 3.4 Linear part — Eq. (16), and tying points (D), (E), (F)

Quoted from journal p. 15:

```text
e^linear_rt = (1/3) c^ (3s - 1);   e^linear_st = (1/3) c^ (1 - 3r)   with

c^ = ( ert(F) - ert(D) ) - ( est(F) - est(E) )                                    (16)
```

The paper's stated purpose for (D), (E), (F): "In order to render the in-plane
twisting stiffness more flexible, the linear part is modified by using three new
tying points (D), (E) and (F) instead of the tying points (1)–(3) when we evaluate
c in Eq. (11)" (journal p. 15, immediately before Eq. (16)). The earlier
(d = 1/6, i.e. edge midpoints) version of this field defines the element the paper
labels **MITC3i**; MITC3+ is the version with D/E/F and `d = 1/10,000`.

Code:

```rust
// mitc3.rs, b_gamma_ext — comment: "Linear part (Eq. 16)"
let b_c_hat = (b_ert_f - b_ert_d) - (b_est_f - b_est_e);
let b_ert_linear = (1.0/3.0) * &b_c_hat * (3.0*s - 1.0);
let b_est_linear = (1.0/3.0) * &b_c_hat * (1.0 - 3.0*r);
```

matches Eq. (16) term by term, including the single shared scalar `c^` used for
both components.

### 3.5 Membrane field — CST

`b_membrane` is the constant-strain-triangle membrane operator over the three
corner nodes:

```rust
bm[(0, u_idx)] = dhi_dx;                     // eps_xx = du/dx
bm[(1, v_idx)] = dhi_dy;                     // eps_yy = dv/dy
bm[(2, u_idx)] = dhi_dy; bm[(2, v_idx)] = dhi_dx;   // gamma_xy = du/dy + dv/dx
```

with `dh` the constant Cartesian derivatives of the linear shape functions
`h1 = 1 - r - s`, `h2 = r`, `h3 = s` — paper Eq. (1). **Note on citations:** the
paper prints the shape functions as Eq. (1) and the covariant strain definition as
Eq. (3), but no numbered *Cartesian* CST operator was located on the extracted
pages (journal pp. 14–15); the correspondence here is therefore to Eq. (1) +
Eq. (3), not to a dedicated equation number.

This constant membrane field is exactly the field that the **strain-smoothed
MITC3+** replaces. That element is a forward reference only and is not described in
this revision (§4.3).

### 3.6 Transverse shear sign convention

The code's covariant shear assembly carries the derivation comment:

```text
// From Ko, Bathe & Zhang 2025 (MITC4/D), Eq. (3a), the offset
// displacement interpolates as
//     u(r,s,t) = sum h_i u_i + (t/2) sum a_i h_i (theta_i x V_in),
// where theta is the physical rotation vector and V_in the shell director.
// Linearizing the director about the flat reference gives
//     V3 = e3 + theta x e3 = e3 + theta_y e1 - theta_x e2,
// hence
//     e_rt = dw/dr + V3.g_r = dw/dr + theta_y g_r[0] - theta_x g_r[1].
// This is the Cartesian convention gamma_xz = w_,x + theta_y and
// gamma_yz = w_,y - theta_x, the same one b_gamma_mitc4 uses.
```

History: commit `b136ce5` ("fix(mitc3): correct covariant shear sign convention
...") **inverted** the director to `V3 ~ e3 - theta_y e1 + theta_x e2`, the inverse
of the Reissner–Mindlin director; commit `d6f37fb` **restored** it, verified
byte-identical to `b136ce5^` for all eight rotation assignments. The consequences
are quantified in §4.1.

---

## 4. Limitations, history and open questions

### 4.1 MITC3 rotation-sign defect and its correction

- Introduced by `b136ce5`; corrected by `d6f37fb`.
- Measured consequences, from the `d6f37fb` commit message (before → after):

  | symptom | before | after |
  | --- | --- | --- |
  | physical rigid-body mode, MITC3 `compute_ke_global` | penalized (`\|Ku\|/\|K_bs\| = 0.786`) | free |
  | mixed mesh, triangles on the last row only | 27.637% error | 0.682% |
  | mixed mesh, triangles on alternating rows | 99.813% error | 1.009% |
  | MITC3Comp B-coupling, [0/90] under axial load | +2.627691e-03 (inverted sign) | -2.627691e-03 |

- Pure meshes were unchanged (all-quad 0.682%, all-triangle 1.015%), which was
  used as the control that this was a convention fix and not a stiffness change.
- Tests were corrected alongside, because they had encoded the inverted convention
  rather than the physics; an `xfail(strict=True)` marker on
  `test_material_suite.py::test_axial_produces_bending_mitc3comp` was removed.
- **This restored fidelity to the 2014 paper, not just internal consistency.** §4.5
  reads Eq. (2) from the rendered page: the paper's own convention is the physical
  one, so `b136ce5` had deviated from MITC3+'s source and `d6f37fb` put the code back
  on it.

### 4.2 Resolved: the MITC4/D drill-membrane operator

This section used to describe the drill-membrane operator as an uncommitted
working-tree experiment that was "being measured". That is no longer true and is
kept here only as history: the operator is committed, it is the source of the
element's sixth DOF, and the element it belongs to is the only shell quadrilateral
in production. Its formulation is §2.5 and its validation is
`docs/validation-matrix.md`.

### 4.3 The strain-smoothed MITC3+

Reference: Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite element",
*Computers and Structures* 223:106096, 2019, doi:10.1016/j.compstruc.2019.07.005 —
`.sources/papers/lee2019.pdf`. The first author is Chaemin Lee, **not** Youngyu Lee
of item 5.

**What it changes, and what it does not.** Only the membrane field. The paper is
explicit: *"We use the originally defined b1 eij and b2 eij in Eqs. (11) and (12)
for the covariant bending strains. For the covariant transverse shear strains, we
adopt the assumed strains of the MITC3+ shell element, in Eqs. (7) and (8)."* So the
bending field (§3.3, §3.4), the transverse shear field (§3.1) and the rotation
convention (§3.6) are the MITC3+ ones unchanged. The smoothing is therefore
orthogonal to the rotation-sign correction of §4.1: it could have been applied
before that fix and would have been equally wrong, and applying it changes nothing
about the convention.

The membrane strain of a target triangle `e` is evaluated at the element centre
(`r = s = 1/3`, `t = 0`) and smoothed with the strains of the three elements across
its edges.

**Eq. (15) — the neighbour's strain in the target's convected coordinates.**

```text
e_ij^(k) = e_ln^(k) (g_i^(e) . g^l^(k)) (g_j^(e) . g^n^(k)),    i, j = 1, 2
```

Because a contravariant base vector is a row of `J^-T` (from `g_i . g^j =
delta_i^j`), `g_i^(e) . g^l^(k) = (J_e J_k^-1)_il`, so the whole transform is
`e~ = M e M^T` with `M = J_e J_k^-1`. `smoothing.rs::convected_operator` is exactly
that, and returns `None` for a singular neighbour Jacobian.

**Eq. (17) — the neighbour's area projected on the target's mid-surface.**

```text
A_bar^(k) = (n^(e) . n^(k)) A^(k)
```

`n` are the unit centre normals, so the factor is `cos(theta)` and the projected area
vanishes at 90 degrees: the smoothing fades to nothing as the two elements become
perpendicular. `smoothing.rs::projected_area`.

**Eq. (16) — the area-weighted pairwise smoothed strain.**

```text
ê_ij^(k) = ( e_ij^(e) A^(e) + e~_ij^(k) A_bar^(k) ) / ( A^(e) + A_bar^(k) )
```

`smoothing.rs::pairwise_smoothed`. **Boundary rule**, stated in the text right after
Eq. (17): *"we use ê_ij = e_ij if the kth edge of the target element is located along
boundary"*. An edge with no neighbour, or with a degenerate one, falls back to the
element's own strain; `smoothed_membrane_strain` implements that by putting the
identity on the target entry.

**Eq. (18) — assignment to the three Gauss points, cyclic.**

```text
e^(A) = (ê^(3) + ê^(1))/2     e^(B) = (ê^(1) + ê^(2))/2     e^(C) = (ê^(2) + ê^(3))/2
```

`smoothing.rs::assign_to_gauss_points`, and the `PAIRS = [[2,0],[0,1],[1,2]]` table in
`smoothed_membrane_strain`. Every edge appears in exactly two of the three assigned
strains.

Eq. (19) gives the equivalent explicit interpolation (`p = 1/6`, `q = 2/3`), but the
paper states it *"is not utilized in actual computation of the stiffness matrix. We
use the assigned strains in Eq. (18) directly in the 3-point Gauss integration"*, so
the code uses Eq. (18) as well.

| Code | Paper | Verified |
| --- | --- | --- |
| `smoothing.rs::tensor_operator` | the `M e M^T` form of Eq. (15) on the engineering strain vector | yes, algebraically |
| `smoothing.rs::convected_operator` | Eq. (15), with `M = J_e J_k^-1` | yes, from the `g_i . g^j = delta_i^j` identity |
| `smoothing.rs::projected_area` | Eq. (17) | yes |
| `smoothing.rs::pairwise_smoothed` | Eq. (16) | yes |
| `smoothing.rs::smoothed_membrane_strain` boundary branch | the boundary rule after Eq. (17) | yes |
| `smoothing.rs::assign_to_gauss_points`, `PAIRS` | Eq. (18) | yes |
| `mitc3.rs::compute_ke_local_with_membrane` | the smoothing changes the membrane term only; `bm_gp[gp]` is the per-Gauss-point membrane operator | yes |
| `mitc3.rs` union layout, `SMOOTHED_UNION_NODES = 6`, `SMOOTHED_UNION_DOFS = 36` | the six-node union of the target and its edge-neighbour node sets | yes |

**Implementation status.** The smoothing primitives are `elements/smoothing.rs`: the
covariant tensor operator of Eq. (15), the edge-neighbour connectivity, the projected
area of Eq. (17), the pairwise average of Eq. (16), the boundary rule and the
assignment of Eq. (18), each with its own unit test. `elements/mitc3.rs` carries the
union layout and `compute_ke_local_with_membrane`, which takes the membrane operator
per Gauss point: passing the element's own `b_membrane` reproduces MITC3+ exactly,
which is what keeps the 2014 element available, and passing the smoothed operators
gives the smoothed element.

**It is not wired into the production path.** Nothing under `crates/aeroelast-py` or
`crates/aeroelast-core/src/assembly` references the smoothed entry points, so it is a
Rust kernel with unit tests rather than an element a user can select. The production
shell quadrilateral is the MITC4+/D of §2 and the production triangle is the MITC3+
of §3.

### 4.4 `quad.rs` has no literature citation

`crates/aeroelast-core/src/elements/quad.rs` carries the comment "elasticity (plane
strain formulation matching the Python reference)" and "Constitutive matrix: plane
strain with Lamé coefficients (matches Python). No thickness multiplier (matching
Python reference)." There is no paper, no author, no DOI, and the justification is
circular: the Rust implementation is validated against the Python implementation, and
nothing establishes which formulation either one is.

**What the code is, read rather than guessed.** The standard displacement-based
plane-strain finite element of any textbook, with isoparametric 4-, 8- and 9-node
quadrilaterals, the Lamé plane-strain constitutive matrix and Gauss integration via
`elements/reference.rs`. That is a formulation one can *name* but not one this
repository *took from* a paper, which is why no citation is invented for it here.

**The gap is declared, not papered over.** Two tempting wrong closures are ruled out
explicitly. The held Choi & Lee 2023 paper
(`.sources/papers/1-s2.0-S0045794922001936-main.pdf`, "Towards improving the 2D-MITC4
element for analysis of plane stress and strain problems") is a **different** element
— an assumed-strain MITC4 for the in-plane problem — and must not be attributed to
this displacement-based kernel. And "matches the Python reference" is not a source
at all. Closing the gap means either citing the textbook the implementation follows,
which is a documentation change, or leaving it explicit; either way this section is
the record that the gap is known and not an oversight.

### 4.5 Resolved: the MITC3+ 2014 rotation convention is the physical one

**Question.** Do the MITC3+ 2014 rotation parameters `alpha`, `beta` in its Eq. (2)
carry the same sign convention as the physical rotation vector `theta` in Ko, Bathe
& Zhang 2025 Eq. (3a), which the code uses?

**Answer: yes, the same convention.** Read from the **rendered page** rather than
from text extraction, the 2014 paper's Eq. (2) is

```text
u(r,s,t) = sum h_i u_i + (t/2) sum a_i h_i ( - V_2^i alpha_i + V_1^i beta_i )
```

with the prose: *"V_1^i and V_2^i are unit vectors orthogonal to V_n^i and to each
other, and alpha_i and beta_i are the rotations of the director vector V_n^i about
V_1^i and V_2^i, respectively, at node i"*
(`.sources/papers/The_MITC3+_shell_element_and_its_performance.pdf`, PDF page 2,
journal p. 13).

A rotation of the director by `alpha` about `V_1` moves a point at offset `t/2` along
the director by `(t/2) alpha (V_1 x V_n)`, and one by `beta` about `V_2` by
`(t/2) beta (V_2 x V_n)`. For a right-handed triad `(V_1, V_2, V_n)`,
`V_1 x V_n = -V_2` and `V_2 x V_n = V_1`, so the offset term is
`(t/2)(-alpha V_2 + beta V_1)`, which is Eq. (2) exactly as printed. That same
expression is `(t/2)(theta x V_n)` for `theta = alpha V_1 + beta V_2`, i.e. Ko, Bathe
& Zhang 2025 Eq. (3a) term for term. **The two papers use the same convention**, and
the code's `V3 = e3 + theta_y e1 - theta_x e2` is the 2014 paper's own form for a
right-handed triad, not a choice imposed from outside it.

**Why this was recorded as open, and what was wrong.** The earlier entry quoted
Eq. (2) as `(V_i2 alpha_i + V_i1 beta_i)`, without the minus sign on the `V_2` term.
That minus is printed in the paper; the text extraction dropped it, and its loss is
what made the two conventions look opposite. The extraction also read the prose as
giving no Cartesian dictionary, which is true but irrelevant: the prose pins the
rotations to `V_1` and `V_2`, which is enough once Eq. (2)'s signs are read
correctly. The conclusion: there was never a conflict to resolve, and the sign
correction of §4.1 restored **fidelity to the 2014 paper**, not merely internal
consistency.

**Caveat that remains, and it is narrow.** The mapping above assumes `(V_1, V_2, V_n)`
is right-handed. The paper leaves the pair `V_1`, `V_2` free up to that choice, so a
reader who takes the opposite handedness gets `(alpha, beta) = (-theta_x, -theta_y)`.
Nothing in the paper fixes the handedness, so "the same convention" is asserted for
the right-handed triad the code uses, which is the standard one for a shell element.

### References

1. Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node shell
   element for general non-linear analysis", *Engineering Computations* 1:77–88,
   1984 (footer prints `Eng. Comput., 1984, Vol. 1, March 77`; last page 88; no DOI
   printed) —
   `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`.
   (Displayed equations not text-extractable from the held scan — §2.4.)
2. Hughes, T.J.R., Taylor, R.L., Kanoknukulchai, W., "A simple and efficient finite
   element for plate bending", *International Journal for Numerical Methods in
   Engineering* 11(10):1529–1543, 1977, doi:10.1002/nme.1620111005. No held copy.
   **Retired source**: this is where the selective reduced integration of the
   superseded hybrid came from. The element described here has none, so the entry is
   kept for the history in §2.6, not because §2 uses it. The DOI was verified against
   the Wiley record and the MaRDI portal entry.
3. Ko, Y., Lee, P.-S., Bathe, K.-J., "A new MITC4+ shell element", *Computers and
   Structures* 182:404–418, 2017, doi:10.1016/j.compstruc.2016.11.004.
   `.sources/papers/1-s2.0-S0045794916309464-main.pdf`.
4. Ko, Y., Bathe, K.-J., Zhang, X., "Continuum mechanics-based shell elements with
   six degrees of freedom at each node — the MITC4/D and MITC4+/D elements",
   *Computers and Structures* 308:107622, 2025,
   doi:10.1016/j.compstruc.2024.107622.
   `.sources/papers/1-s2.0-S0045794924003511-main.pdf`. (`.sources/papers/MITC_D_elements_published.pdf`
   is a scanned duplicate of the same paper with no text layer.)
5. Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element and its
   performance", *Computers and Structures* 138:12–23, 2014,
   doi:10.1016/j.compstruc.2014.02.005.
   `.sources/papers/The_MITC3+_shell_element_and_its_performance.pdf`.
6. Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite element",
   *Computers and Structures* 223:106096, 2019,
   doi:10.1016/j.compstruc.2019.07.005.
   `.sources/papers/lee2019.pdf` (authors: Chaemin Lee, Phill-Seung Lee; the first
   author is **not** Youngyu Lee of item 5). **Implemented** — §4.3 describes it
   equation by equation; verified against the held copy; entry present in
   `docs/references.md`.

`docs/references.md` is the canonical bibliography for the repository. Where this
document and `docs/references.md` disagree, `docs/references.md` should be
corrected — this document deliberately does not duplicate its per-entry verification
annotations. The Hughes–Brezzi drilling attribution that §2.6 used to depend on is
recorded there in §1, together with the note that it is not the source of the
current element.
