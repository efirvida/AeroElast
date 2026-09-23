# Shell element formulations: code-to-equation reference

Scope of this revision: the two quadrilateral/triangular continuum-mechanics-based
elements that are currently in production in `crates/aeroelast-core/src/elements/`
— `mitc4.rs` (MITC4+) and `mitc3.rs` (MITC3+). Every equation below was read from
the PDF in `.sources/papers/` named with it; equations that could not be read are
marked as such rather than reconstructed.

---

## 1. Scope and conventions

### 1.1 Elements covered in this version

| Element | File | Nodes / DOF | Status |
| --- | --- | --- | --- |
| MITC4+ | `crates/aeroelast-core/src/elements/mitc4.rs` | 4 × 6 = 24 | documented here |
| MITC3+ | `crates/aeroelast-core/src/elements/mitc3.rs` | 3 × 6 + 2 = 20 → 18 | documented here |

**Pending, not covered by this revision:**

- MITC4 (original) and MITC4/D — the six-DOF (drill-including) variant of Ko, Bathe
  & Zhang 2025. The drill-membrane operator exists only as an uncommitted
  working-tree experiment (§4.2).
- The strain-smoothed MITC3+ of Lee & Lee 2019 (§4.3).
- The plane-stress/plane-strain quadrilaterals and triangles in
  `crates/aeroelast-core/src/elements/` (`quad.rs` and relatives), which carry no
  literature citation (§4.4).

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

## 2. MITC4+ (quadrilateral, 24 DOF)

Reference: Ko, Y., Lee, P.-S., Bathe, K.-J., "A new MITC4+ shell element",
*Computers and Structures* 182:404–418, 2017, doi:10.1016/j.compstruc.2016.11.004 —
`.sources/papers/1-s2.0-S0045794916309464-main.pdf` (duplicate:
`A_new_MITC4+_shell_element.pdf`).

### 2.1 As implemented

- **Nodes / integration.** 4 nodes × 6 DOF. `N_GAUSS = 4`; `GAUSS_XI/ETA = ±1/√3`,
  weights `1.0` — i.e. the 2×2 Gauss-Legendre rule on `[-1,1]^2`, matching the
  paper's "we use 2×2×2 Gauss integration over the element domain" (journal p. 410).
- **Membrane field.** The assumed covariant membrane strain of Eq. (27) below,
  sampled at five tying points. Code tying points (`mitc4.rs`, precomputed in
  `Mitc4Precomputed::new`):

  ```text
  A(0,+1)  -> pre.b_rr_a   covariant row  rr
  B(0,-1)  -> pre.b_rr_b   covariant row  rr
  C(+1,0)  -> pre.b_ss_c   covariant row  ss
  D(-1,0)  -> pre.b_ss_d   covariant row  ss
  E( 0,0)  -> pre.b_rs_e   covariant row  rs
  ```

  The paper names the same five tying points (A)–(E) (Fig. 4, "Tying positions
  (A)–(E) for the assumed membrane shear strain field") but presents them
  graphically: the parametric coordinates quoted above are the code's own
  parameterisation, not values read off the paper.
- **Blending coefficients.** `compute_membrane_coefficients(&x_d, &m_r, &m_s)`
  returns `(a_a..a_e)` from `c_r = x_d·m_r`, `c_s = x_d·m_s`,
  `d = c_r² + c_s² - 1` — the paper's Eq. (24) plus the `a_A..a_E` block printed
  immediately after Eq. (27c), both reproduced in §2.3.
- **Bubble.** `N_b = (1 - xi²)(1 - eta²)` (`bubble_function`), carrying **2 rotation
  DOF only**, statically condensed:
  `K = K_nn - K_nb K_bb^{-1} K_nb^T` with a regularized 2×2 inverse
  (`regularized_inverse_2x2`), symmetrized as
  `K = 0.5 (K_raw + K_raw^T)`. The 26-DOF extended space appears only inside
  `compute_ke_local`; the bubble block is never exported.
- **Transverse shear.** `b_gamma_mitc4` is the rotation-based MITC4 (Dvorkin &
  Bathe) interpolation — see §2.4. `b_gamma_mitc4_plus` adds the bubble coupling
  `Bs_bubble = [[0, N_b], [-N_b, 0]]`.
- **Drilling.** A separate rank-1 penalty, §2.6.
- **Not implemented.** No hourglass control is applied in the production
  stiffness path despite the "S4R-style" scaffolding (`hg_factor`, `h_vec`,
  `hg_stiffness_factor`, `h_orth`) precomputed in the struct and stored in
  `Mitc4Precomputed`. `compute_ke_local` returns
  `k_m + k_mb_coup + k_mb + k_bs + k_drill` with no hourglass term, and
  `compute_hourglass_stiffness` / `compute_hourglass_forces` are referenced only
  from `#[cfg(test)]` tests in the same file. The scaffolding is dead weight as
  far as the assembled stiffness is concerned.

### 2.2 Code → equation table

| Code | Paper | Verified |
| --- | --- | --- |
| `b_m_mitc4_plus` row `b_rr` | Eq. (27a) | yes, term by term |
| `b_m_mitc4_plus` row `b_ss` | Eq. (27b) | yes, term by term |
| `b_m_mitc4_plus` row `b_rs`, stacked as `2*b_rs` | Eq. (27c) | yes, term by term |
| `pre.a_a .. pre.a_e` | `a_A .. a_E` after Eq. (27c), from Eq. (24) | yes, term by term |
| `b_gamma_mitc4` | Dvorkin & Bathe 1984 transverse shear, reproduced unnumbered in Ko et al. 2017 §2 | partially — see §2.4 |
| `cm_normal` / centre-point shear term in `compute_ke_local` (SRI) | Hughes, Taylor & Kanoknukulchai 1977 | not verifiable from a held copy; see §2.5 |
| `b_drill` + `k_drill` | code comment names "Hughes & Brezzi" with no title, year or venue | unattributed in practice; see §2.6 |

The code's own comments agree with these numbers for the membrane rows
(`// Blended covariant B-rows (Ko et al. 2017, Eqs. 27a-c)`).

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

**Code correspondence (verified term by term).**

```rust
// mitc4.rs, b_m_mitc4_plus — rows b_rr, b_ss, b_rs
b_rr = 0.5*(1 - 2*a_a + s + 2*a_a*s*s) * pre.b_rr_a
     + 0.5*(1 - 2*a_b - s + 2*a_b*s*s) * pre.b_rr_b
     + a_c*(-1 + s*s) * pre.b_ss_c
     + a_d*(-1 + s*s) * pre.b_ss_d
     + a_e*(-1 + s*s) * pre.b_rs_e;                                   // = (27a)

b_ss = a_a*(-1 + r*r) * pre.b_rr_a
     + a_b*(-1 + r*r) * pre.b_rr_b
     + 0.5*(1 - 2*a_c + r + 2*a_c*r*r) * pre.b_ss_c
     + 0.5*(1 - 2*a_d - r + 2*a_d*r*r) * pre.b_ss_d
     + a_e*(-1 + r*r) * pre.b_rs_e;                                   // = (27b)

b_rs = 0.25*(r + 4*a_a*r*s) * pre.b_rr_a
     + 0.25*(-r + 4*a_b*r*s) * pre.b_rr_b
     + 0.25*(s + 4*a_c*r*s) * pre.b_ss_c
     + 0.25*(-s + 4*a_d*r*s) * pre.b_ss_d
     + (1 + a_e*r*s) * pre.b_rs_e;                                    // = (27c)
```

and in `compute_membrane_coefficients`:

```rust
a_a = c_r*(c_r - 1.0)/(2.0*d);   a_b = c_r*(c_r + 1.0)/(2.0*d);
a_c = c_s*(c_s - 1.0)/(2.0*d);   a_d = c_s*(c_s + 1.0)/(2.0*d);
a_e = -2.0*c_r*c_s/d;
```

One structural detail: the code stacks the covariant rows as
`B_cov = [B_rr; B_ss; 2*B_rs]` (`b_cov[(2,j)] = 2.0*b_rs[j]`). The `2*B_rs` is the
engineering-shear factor applied at the covariant→local mapping, not a change to
Eq. (27c).

### 2.4 Transverse shear — Dvorkin & Bathe 1984

Reference: Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node
shell element for general non-linear analysis", *Engineering Computations*,
1984 — `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`.

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

### 2.5 Selective reduced integration (SRI) in `compute_ke_local`

Code, `compute_ke_local` (comment block beginning "Despite MITC4+ blending,
residual in-plane shear locking persists ..."):

```rust
let mut cm_normal = *cm;
cm_normal[(0,2)] = 0.0; cm_normal[(1,2)] = 0.0;
cm_normal[(2,0)] = 0.0; cm_normal[(2,1)] = 0.0;
cm_normal[(2,2)] = 0.0;
let c_shear = cm[(2,2)];
// 4-GP loop:  k_m += (bm^T * cm_normal * bm) * (w * sqrt_g)
// centre point: k_m += b_shear^T * c_shear * b_shear * (4 * sqrt_g_c)
//               with b_shear = row 2 of b_m_mitc4_plus(pre, 0.0, 0.0)
```

That is: **full 2×2 for the normal membrane components** (`eps_xx`, `eps_yy`, and
the coupling terms struck out) and **one centre point at `xi = eta = 0` for the
in-plane shear** (`2*eps_xy`), with weight `4 * sqrt_g` (the 2×2 rule's total
weight) at the centre. The same split is repeated in `compute_fint_global` for
Newton tangent consistency (comments at `mitc4.rs` on the "SRI-consistent virtual
work").

**Citation.** The comment names "Hughes, Taylor & Kanoknukulchai (1977) — SRI for Q4
membrane element", and commit `929db32` ("fix(mitc4): apply selective reduced
integration to eliminate in-plane bending locking") ends with `Ref: Hughes, Taylor
& Kanoknukulchai (1977)`.

**What was and was not confirmed:**

- Confirmed: the code and the commit both attribute the split to that paper;
  `docs/references.md` lists Hughes, Taylor & Kanoknukulchai, "A simple and
  efficient finite element for plate bending", *IJNME* 11(10):1529–1543, 1977.
- **Not confirmed:** no copy of that paper is in `.sources/papers/`, so the
  specific claim — that this paper prescribes exactly a 2×2 rule for the normal
  membrane components and a one-point rule for the in-plane shear — could not be
  verified. `docs/references.md` itself marks the DOI "to verify" and states that
  the journal, volume and pages come from the published record and are **not**
  re-verified against a held copy. The 1977 attribution should be treated as
  repository assertion, not as verified evidence.

### 2.6 Drilling penalty

`b_drill` (`mitc4.rs`), assembled in `compute_ke_local` over the same 2×2 rule:

```rust
bd[u_idx]   = -0.5 * dni_dy;
bd[v_idx]   =  0.5 * dni_dx;
bd[thz_idx] = -n_vals[i];
// k_drill += (bd * bd^T) * (pre.k_drill * w * sqrt_g)
// pre.k_drill = e_mod * thickness^2 * 0.15 * drilling_scale
```

i.e. the sum of `k_drill = 0.15 * E * h^2 * drilling_scale` weighted against the
drilling strain `gamma_xz_drill = (v_,x - u_,y)/2 - theta_z` (the `mitc3.rs`
`b_drill` carries the same comment `// (dv/dx - du/dy)/2 - thetaz`).

**Sourcing, plainly stated.** The doc comment on `b_drill` reads
`/// Drilling B-vector (Hughes & Brezzi, 1×24)`. That is the whole attribution: two
author surnames and no title, year, journal or volume. `docs/references.md` has **no**
Hughes–Brezzi entry, and no such paper is held in `.sources/papers/`. The `0.15`
penalty factor and the `drilling_scale` multiplier are **not** attributed to
anything at all in the code. So: the operator cites two names without a work, and
the specific penalty expression — including the `0.15` coefficient — is uncited.
This should be resolved before the MITC4/D work (§4.2) lands.

---

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
  with `k_drill = E * h^2 * 0.15 * drilling_scale` — the same uncited expression as
  §2.6, and here with **no** source comment at all.

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

## 4. Known deviations, pending work and references

### 4.1 MITC3 rotation-sign defect and its correction

- Introduced by `b136ce5`; corrected by `d6f37fb`.
- Measured consequences, from the `d6f37fb` commit message (before → after):

  | symptom | before | after |
  | --- | --- | --- |
  | physical rigid-body mode, MITC3 `compute_ke_global` | penalized (`|Ku|/|K_bs| = 0.786`) | free |
  | mixed mesh, triangles on the last row only | 27.637% error | 0.682% |
  | mixed mesh, triangles on alternating rows | 99.813% error | 1.009% |
  | MITC3Comp B-coupling, [0/90] under axial load | +2.627691e-03 (inverted sign) | -2.627691e-03 |

- Pure meshes were unchanged (all-quad 0.682%, all-triangle 1.015%), which was
  used as the control that this was a convention fix and not a stiffness change.
- Tests were corrected alongside, because they had encoded the inverted convention
  rather than the physics; an `xfail(strict=True)` marker on
  `test_material_suite.py::test_axial_produces_bending_mitc3comp` was removed.

### 4.2 Pending: MITC4/D drill-membrane operator

- `b_md_mitc4_plus` and `drill_midside_shape_derivatives` in `mitc4.rs`. The
  comment cites "Ko et al. 2025, Eq (11)" for the simplified midside shape
  derivatives, with deliberately zeroed entries "to avoid higher-order integration
  than the base element".
- History: introduced in `929db32`; deleted in `4a46af2` ("delete the unwired
  drill-membrane operator and its tests" — its two unit tests were vacuous, because
  they evaluated the operator at `(0,0)` where all four entries are identically
  zero by construction). It is **present again in the current working tree only as
  an uncommitted change** (`git diff --stat` shows `mitc4.rs` +95 lines), where it
  is wired into the SRI loop as
  `let bm = b_m_mitc4_plus(pre, xi, eta) + b_md_mitc4_plus(pre, xi, eta);`.
- It is being measured. No values are recorded in this document because no
  validated measurement was read.

### 4.3 Pending: strain-smoothed MITC3+

Lee, Y., Lee, P.-S., "The strain-smoothed MITC3+ shell element", planned. The
intent is to replace the constant-strain triangle membrane field of §3.5 with a
smoothed field. `docs/references.md` currently has **no** entry for this work and
there is no held copy, so the volume/pages/DOI must be verified before it is cited
anywhere in the repository. Not implemented in this revision.

### 4.4 `quad.rs` has no literature citation

`crates/aeroelast-core/src/elements/quad.rs` carries the comment "elasticity (plane
strain formulation matching the Python reference)" and "Constitutive matrix: plane
strain with Lamé coefficients (matches Python). No thickness multiplier (matching
Python reference)." There is no paper, no author, no DOI. The justification is
circular: the Rust implementation is validated against the Python implementation,
and nothing establishes which formulation either one is. This is a documentation
defect, not a numerical one, but it blocks any future claim that the plane
quadrilaterals are a published formulation.

### 4.5 Open: sign of `alpha`, `beta` in the MITC3+ 2014 paper

**Question.** Do the MITC3+ 2014 rotation parameters `alpha`, `beta` in its Eq. (2)
carry the same sign convention as the physical rotation vector `theta` in Ko, Bathe
& Zhang 2025 Eq. (3a), which the code now uses?

**What was tried.** Reading the MITC3+ 2014 paper directly:

- Eq. (2) reads `u(r,s,t) = sum h_i u_i + (t/2) sum a_i h_i (V_i2 alpha_i + V_i1 beta_i)`.
- The only definition given is prose: "`alpha_i` and `beta_i` are the rotations of
  the director vector `V_in` about `V_i1` and `V_i2`, respectively" (journal p. 14).
  There is no Cartesian dictionary and no relation to a physical rotation vector.
- Eqs. (3), (8) and (14)–(17) are written entirely in covariant components, so they
  do not disambiguate either.

**Result: unresolved.** Against Ko et al. 2025 Eq. (3a), `theta x V_in` with
`V_in = e3` expands to `theta_y V_i2 - theta_x V_i1` (given `V_i1 = e1`,
`V_i2 = e2`), which would map `alpha -> theta_y` and `beta -> -theta_x`. Nothing in
the 2014 paper confirms or refutes that mapping. It is stated here as open.

### References

1. Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node shell
   element for general non-linear analysis", *Engineering Computations*, 1984.
   `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`.
   (Displayed equations not text-extractable from the held scan — §2.4.)
2. Hughes, T.J.R., Taylor, R.L., Kanoknukulchai, W., "A simple and efficient finite
   element for plate bending", *International Journal for Numerical Methods in
   Engineering* 11(10):1529–1543, 1977. DOI to verify. No held copy — §2.5.
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
6. Lee, Y., Lee, P.-S., "The strain-smoothed MITC3+ shell element", *Computers and
   Structures* 223:106096, 2019 — **planned; not held, not verified, not in
   `docs/references.md`**. §4.3.

`docs/references.md` is the canonical bibliography for the repository. Where this
document and `docs/references.md` disagree, or where `docs/references.md` is
missing an entry (items 6 and the Hughes–Brezzi attribution of §2.6),
`docs/references.md` should be corrected — this document deliberately does not
duplicate its per-entry verification annotations.
