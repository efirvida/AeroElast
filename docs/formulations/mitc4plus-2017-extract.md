# MITC4+ formulation extract — Ko, Lee & Bathe (2017), C&S 182:404–418

> **Internal working notes, not production documentation.** This file is a page-by-page
> transcription of the paper with the observations made while reading it: work-unit ids,
> gate outcomes, commit ids and open questions from the development process. It exists so
> that a later reader can check a specific equation against the PDF. For the element as it
> is in production, read `docs/formulations/shell-elements.md`; for the bibliography,
> `docs/references.md`. Nothing here should be cited as a specification.

A persistent transcription of the equations and figures this repository's MITC4+
implementation is checked against, so that verifying the code does not require
re-reading the PDF each time.

**This is an extract, not a substitute for the paper.** Every equation below was
read from the PDF with vision. Nothing here is reconstructed from memory, and the
gaps are marked as gaps rather than filled in.

## The two 2017 papers, and why the distinction matters

| | Paper | Role | File |
|---|---|---|---|
| **A** | Ko, Lee & Bathe (2017), *"A new MITC4+ shell element"*, **Computers and Structures 182:404–418** | **The formulation.** Defines the assumed membrane and shear fields. This document. | `.sources/papers/A_new_MITC4+_shell_element.pdf` |
| **B** | Ko, Lee, Lee & Bathe (2017), *"Performance of the MITC3+ and MITC4+ shell elements in widely-used benchmark problems"*, **Computers and Structures 193:187–206** | **The benchmarks.** Publishes the values the tests compare against (Tables 3–19). | `.sources/papers/1-s2.0-S0045794917309550-main.pdf` |

Paper **B** is the judge, paper **A** is the law. Citing one for the other is a
real error: an earlier comment in `mitc4.rs` cited *"Ko et al. 2017, Eqs. 27a-c"*
for the membrane blending, and neither paper has such equations.

## The element as paper A actually defines it (full re-read, 2026-09-23)

This section was written after reading PDF pages 2–7 (journal pp. 405–410) with
vision. It exists because the code had drifted far enough that the element it
implements is **not** this one, and the difference is architectural, not a
coefficient.

**The kinematics is continuum-mechanics-based and 3D (Eqs. 1–3, p. 405).**

```text
x(r,s,t) = Σ_{i=1..4} h_i(r,s) x_i + (t/2) Σ_{i=1..4} a_i h_i(r,s) V_n^i      (1)
u(r,s,t) = Σ_{i=1..4} h_i(r,s) u_i + (t/2) Σ_{i=1..4} a_i h_i(r,s) (−V_2^i α_i + V_1^i β_i)   (3)
```

- `a_i` is the **shell thickness at node i** and `V_n^i` the **director vector**
  at node i (the paper says so in words, right after Eq. 1). Fig. 3(d)'s red
  arrows are the *edge* vectors `x_e^k` of Eq. (13), which are a different thing.
- The through-thickness coordinate is `t ∈ [−1, 1]`, scaled by `t/2`.
- **The DOF per node are 5, not 6**: three translations `u_i` plus two director
  rotations `α_i` (about `V_1^i`) and `β_i` (about `V_2^i`). **There is no
  drilling rotation in this element at all.** The rotation of the director about
  `V_n^i` does not appear in Eq. (3), so it is not a DOF.

**The strain measure is the covariant 3D strain (Eqs. 4–5):**

```text
e_ij = ½ (g_i · u_j + g_j · u_i),   g_i = ∂x/∂r_i,  u_i = ∂u/∂r_i,
with r_1 = r, r_2 = s, r_3 = t.                                              (4), (5)
```

**Only the membrane term is modified (p. 406, in words):**

> "The first term `e_ij^m` in Eq. (7a) is the covariant in-plane membrane strain
> at the shell mid-surface (t = 0), and the remaining terms are the covariant
> in-plane strains due to bending.  The in-plane membrane strain, see Eq. (7b),
> can in general induce locking and it is this term that we modify as described
> below; **we leave the other terms in Eq. (7a) as they are and evaluate them
> using the displacement formulation**."

So the bending terms (7c)/(7d) — including their `∂x_b` warping parts — are **not**
assumed; they are the plain displacement-based ones.  The MITC4+ contribution is
the assumed membrane field alone.

**Integration and the absence of any numerical factor (p. 410, in words):**

> "In the numerical solutions, we use **2 × 2 × 2 Gauss integration** over the
> element domain for all shell elements considered."
>
> "We note that **the element formulation does not include any numerical
> factor**, and consider next the isotropy, zero energy mode and patch tests."

**The paper's own basic tests are: isotropy, zero-energy modes, patch tests**
(Section 4, "Basic numerical tests", pp. 410–411).  That is the acceptance list
for *"pasa los tests teóricos de su formulación"*, and paper A states the three
precisely:

> **Isotropy.** "... sequence of node numbering, i.e. on the element orientation
> [1,3,15–18].  The element **passes the test of spatial isotropy**."
> → the single-element stiffness must be invariant under the element's
> orientation / the node-numbering sequence.
>
> **Zero energy mode.** "In the zero energy mode test, the number of zero
> eigenvalues of the stiffness matrix of a **single unsupported element** are
> counted [1–3,9–18].  For the new element only the **six zero eigenvalues
> corresponding to the six rigid body modes** are obtained.  That is, the element
> passes the zero energy mode test."
> → exactly six, no more, on an unsupported element.
>
> **Patch tests.** "We perform **three patch tests: the membrane, bending and
> shearing patch tests**, see Refs. [1–3,9–18].  The mesh geometry is shown in
> **Fig. 5**.  The patch of elements is subjected to the minimum number of
> constraints to prevent rigid body motions and the nodal point forces on the
> boundary corresponding to the constant stress states are applied.  The patch
> tests are passed if the correct values of **constant stress fields** are
> calculated at any location within the mesh.  The element passes the membrane,
> bending and shearing patch tests."
> → three separate patch tests, on the **distorted** mesh of Fig. 5 (the 10 × 10
> square whose interior nodes are at (2,2), (4,7), (8,7), (8,3)), with the minimum
> constraints against rigid-body motion and boundary nodal forces from the
> constant stress states.
The convergence studies (Section 5) measure the error in the **s-norm of Hiller
& Bathe**, Eq. (28), with the exact solution replaced by a very fine-mesh
reference `u_ref` when no analytical solution exists.

**The transverse shear (p. 405):** the MITC4 assumed field of Dvorkin & Bathe
(1984), tying points A (top, s=+1), B (bottom, s=−1), C (right, r=+1),
D (left, r=−1):

```text
e_rt = ½(1+s) e_rt^(A) + ½(1−s) e_rt^(B)     e_st = ½(1+r) e_st^(C) + ½(1−r) e_st^(D)
```

**Sign correction to §2.3 below.** The paper prints (p. 410, after Eq. 27c)
`a_E = 2 c_r c_s / d` — **positive**.  §2.3 of this extract wrote
`aE = -2 cr cs / d`, and the deleted `compute_membrane_coefficients` implemented
the negative form.  Both are wrong.

**The five architectural deviations of this repository's element.**  For the
audit, the things that do NOT belong to paper A:

| # | paper A | this repository |
|---|---|---|
| 1 | 5 DOF/node (3 translations + 2 director rotations), no drilling | 6 DOF/node with a drilling rotation |
| 2 | 3D continuum kinematics, director `V_n^i`, thickness coordinate `t` | flat projection onto a local 2D frame, ABD resultants |
| 3 | director enrichment (Eq. 8b) is part of the kinematics | a separate 2-DOF `(1−ξ²)(1−η²)` rotation bubble |
| 4 | full 2×2×2 Gauss integration, no numerical factor | selective reduced integration of the in-plane shear |
| 5 | membrane Eqs. (21)–(27) | membrane Eqs. (18)–(19) (`compute_membrane_coefficients` deleted) |

Plus a sixth that is material rather than kinematic: the repository's transverse
shear uses a shear correction factor, while paper A states there is no numerical
factor.  And, on top of the element itself, the repository layers a Winkler &
Plakomytis ERC drilling treatment and a `beta_w` warping penalty — both from a
different element in a different paper.

## How to read the PDF

```bash
pdftoppm -png -f <page> -l <page> -r 320 <pdf> /tmp/out
# then crop with PIL and read the image
```

**Never use `pdftotext` for equations.** Its default extraction interleaves the
two columns and loses superscripts; with `-layout` it still mangles two-column
equation blocks. Both failure modes have already produced wrong readings in this
project.

## Element geometry and kinematics

**Eq. (7a)** — the covariant in-plane strain decomposition (p. 406):

```text
e_ij = e_ij^m + t·e_ij^b1 + t²·e_ij^b2        with i,j = 1,2
```

**Eq. (7b)** — the membrane part. *This is the term the MITC4+ modifies*:

```text
e_ij^m = ½( ∂x_m/∂r_i · ∂u_m/∂r_j + ∂x_m/∂r_j · ∂u_m/∂r_i )
```

**Eq. (7c)** — the first bending part:

```text
e_ij^b1 = ½( ∂x_m/∂r_i · ∂u_b/∂r_j + ∂x_m/∂r_j · ∂u_b/∂r_i
           + ∂x_b/∂r_i · ∂u_m/∂r_j + ∂x_b/∂r_j · ∂u_m/∂r_i )
```

**Eq. (7d)**:

```text
e_ij^b2 = ½( ∂x_b/∂r_i · ∂u_b/∂r_j + ∂x_b/∂r_j · ∂u_b/∂r_i )
```

**Eq. (8a)** — characteristic geometry and its enrichment (p. 406):

```text
x_m = Σ_{i=1..4} h_i(r,s) x_i
x_b = ½ Σ_{i=1..4} a_i h_i(r,s) V_n^i
```

**Eq. (8b)** — the displacement field and its enrichment:

```text
u_m = Σ_{i=1..4} h_i(r,s) u_i
u_b = ½ Σ_{i=1..4} a_i h_i(r,s) ( −V_2^i α_i + V_1^i β_i )
```

`a_i` is the **shell thickness at node i** (paper A, p. 405, right after Eq. 1:
*"a_i and V_n^i denote the shell thickness and the director vector at the node"*),
`V_n^i` the nodal normal, `V_1^i`, `V_2^i` the nodal in-plane vectors and `α_i`,
`β_i` the two director rotations.  Fig. 3(d)'s red arrows are the *edge* vectors
`x_e^k` of Eq. (13), a different quantity.  Note that **`x_b` vanishes for a flat
element** (the edge vectors have no normal component), so this enrichment is the
element's *warping* treatment, and it is part of the kinematics rather than a
separate bubble.

**Relations following from Eq. (2) in Eqs. (8a) and (8b)** (p. 406):

```text
∂x_m/∂r = x_r + s·x_d        ∂x_m/∂s = x_s + r·x_d
∂u_m/∂r = u_r + s·u_d        ∂u_m/∂s = u_s + r·u_d
```

**Eq. (9)** — the characteristic vectors, with ξ_i, η_i = ±1:

```text
x_r = ¼ Σ ξ_i x_i     x_s = ¼ Σ η_i x_i     x_d = ¼ Σ ξ_i η_i x_i
u_r = ¼ Σ ξ_i u_i     u_s = ¼ Σ η_i u_i     u_d = ¼ Σ ξ_i η_i u_i
```

`x_d` connects the centres of the two diagonals. It is what makes the
displacement-based `e_rr^m = (x_r + s·x_d)·(u_r + s·u_d)` **quadratic** in `s`.

**Eq. (10)** — the element plane normal:

```text
n = (x_r × x_s) / ‖x_r × x_s‖
```

**Eq. (11)** — the dual basis vectors on the plane P:

```text
m^ri · x_rj = δ_ij,   m^ri · n = 0     with r_1 = r, r_2 = s
```

**Eq. (12)** — the distortion vector decomposed into in-plane and out-of-plane
parts:

```text
x_d = (x_d · m^r) x_r + (x_d · m^s) x_s + (x_d · n) n
```

`m^r · x_d` and `m^s · x_d` are the in-plane distortions, `x_d · n` the
out-of-plane distortion. The paper notes the length of `x_d` becomes non-zero for
both in-plane and out-of-plane distortions, and that membrane locking occurs
because of out-of-plane distortions of the geometry.

**Eq. (13)** — the edge vectors:

```text
x_e^1 = (x_2 − x_1)/2 = −x_r − x_d = −∂x_m/∂r (0,  1)
x_e^2 = (x_3 − x_2)/2 = −x_s + x_d = −∂x_m/∂s (−1, 0)
x_e^3 = (x_4 − x_3)/2 =  x_r − x_d =  ∂x_m/∂r (0, −1)
x_e^4 = (x_1 − x_4)/2 =  x_s + x_d =  ∂x_m/∂s (1,  0)
```

**Eq. (14)** — the edge strains, each containing only two nodal displacements:

```text
e_rr^m(0,  1) = x_e^1 · u_e^1        e_rr^m(0, −1) = x_e^3 · u_e^3
e_ss^m(1,  0) = x_e^4 · u_e^4        e_ss^m(−1, 0) = x_e^2 · u_e^2
```

The paper is explicit about why this matters: *"The use of the edge strains in
Eq. (14) is important to establish an improved behavior in bending-dominated
problems."* Checked against Eq. (15) at the tying points: the two agree, because
Eq. (17) says the sampled value at A *is* `e_rr|con + e_rr|lin + e_rs|bil`, which
is exactly Eq. (15) evaluated at (0,1). So Eq. (14) is an efficient way to compute
the same quantity, not a different one.

## Assumed membrane strain field (Section 3.2, p. 408)

**Fig. 4** — the five tying points for the assumed membrane field:

| point | (r, s) | sampled component |
|---|---|---|
| A | (0, +1) | `e_rr^m(A)` |
| B | (0, −1) | `e_rr^m(B)` |
| C | (+1, 0) | `e_ss^m(C)` |
| D | (−1, 0) | `e_ss^m(D)` |
| E | (0, 0) | `e_rs^m(E)` |

**Eq. (17)** — the five sampled strains, decomposed:

```text
e_rr^m(A) = e_rr^m|con + e_rr^m|lin + e_rs^m|bil
e_rr^m(B) = e_rr^m|con − e_rr^m|lin + e_rs^m|bil
e_ss^m(C) = e_ss^m|con + e_ss^m|lin + e_rs^m|bil
e_ss^m(D) = e_ss^m|con − e_ss^m|lin + e_rs^m|bil
e_rs^m(E) = e_rs^m|con
```

**Eq. (18)** — the assumed field, linear in the respective coordinate:

```text
ẽ_rr^m = ½(e_rr^m(A) + e_rr^m(B)) + ½(e_rr^m(A) − e_rr^m(B))·s
ẽ_ss^m = ½(e_ss^m(C) + e_ss^m(D)) + ½(e_ss^m(C) − e_ss^m(D))·r
ẽ_rs^m = e_rs^m(E) = e_rs^m|con
```

**Eq. (19)** — with the linear shear terms the patch test requires:

```text
ẽ_rr^m = ẽ_rr^m
ẽ_ss^m = ẽ_ss^m
ẽ_rs^m = ẽ_rs^m + ½·e_rr^m|lin·r + ½·e_ss^m|lin·s
```

**Eq. (15)** — the displacement-based membrane strains in terms of the
characteristic vectors:

```text
e_rr^m = e_rr^m|con + e_rr^m|lin·s + e_rs^m|bil·s²
e_ss^m = e_ss^m|con + e_ss^m|lin·r + e_rs^m|bil·r²
e_rs^m = e_rs^m|con + ½e_rr^m|lin·r + ½e_ss^m|lin·s + e_rs^m|bil·rs
```

**Eq. (16)** — the parts:

```text
e_rr^m|con = x_r · u_r        e_ss^m|con = x_s · u_s
                            e_rs^m|con = ½(x_r · u_s + x_s · u_r)
e_rr^m|lin = x_r · u_d + x_d · u_r
e_ss^m|lin = x_s · u_d + x_d · u_s
e_rs^m|bil = x_d · u_d
```

This is the field our `compute_covariant_membrane_b_row` builds: its
`g_r·(∂u/∂r)` with `g_r = x_r + s·x_d` and `∂u/∂r = u_r + s·u_d` expands to
`x_r·u_r + s(x_r·u_d + x_d·u_r) + s²(x_d·u_d)`, term for term Eq. (15).

**Eq. (20)** — the inverse relations, comparing Eq. (19) with the
displacement-based field of Eq. (15):

```text
e_rr^m = ẽ_rr^m − e_rs^m|bil + e_rs^m|bil·s²
e_ss^m = ẽ_ss^m − e_rs^m|bil + e_rs^m|bil·r²
e_rs^m = ẽ_rs^m + e_rs^m|bil·r·s
```

The paper's own summary: the assumed field is *"one order lower than implicitly
given in the original displacement-based element"* — Eq. (18) deliberately
discards the quadratic part that Eq. (20) shows the displacement-based field has.

**Correction (this repository, 2026-09-23).** The paragraph that used to sit here
claimed *"there are no geometry-dependent coefficients in Eqs. (18)–(19)"* and
that the five `a_a`..`a_e` coefficients were not in the paper. **That reading was
wrong.** Eq. (18) is only the Choi–Paik starting field; the paper's *new MITC4+*
field is Eqs. (21)–(27), which do carry five geometry-dependent coefficients.
Commit `dc7593e` removed them as "unsourced" on the basis of that wrong reading.
They are reproduced below.

**Eq. (21)–(22)** — the assumed bilinear term and the patch-test condition:

```text
ẽ_rs^m|bil = B1·(e_rr^m|con + e_rs^m|bil) + B2·(e_ss^m|con + e_rs^m|bil)
          + B3·e_rs^m|con + B4·e_rr^m|lin + B5·e_ss^m|lin
ẽ_rs^m|bil = e_rs^m|bil   when the element geometry is flat (x_d · n = 0)
```

`transcription-verified: 2026-09-23 p.409`

**Note F2 — CORRECTED 2026-09-24: the printed Eq. (21) is NOT defective.** The
earlier note here (2026-09-23, WU0) claimed that "as printed, Eq. (21) has no
leading term", "proved" by Eq. (27c)'s `(1 + a_E·rs)` coefficient. A vision
re-read of p. 409 and of **Appendix A, p. 416** refutes it: Eq. (21) printed does
contain `e_rs^m|bil`, inside the `B_1` and `B_2` parentheses, with coefficient
`B_1 + B_2 = (c_r²+c_s²)/d = 1 + 1/d`; and the paper's own **Eq. (A.7)** — the
substitution of Eq. (A.6) into Eq. (21) — is
`(c_r²/d)(e_rr|con + e_rs|bil) + (c_s²/d)(e_ss|con + e_rs|bil) + (2c_rc_s/d)e_rs|con − (c_r/d)e_rr|lin − (c_s/d)e_ss|lin`,
i.e. Eq. (25) with the `e_rs|bil` terms retained. Eq. (21) and Eq. (25) are the
same equation with `B_1..B_5` substituted. The `(1 + a_E·rs)` coefficient of
Eq. (27c) multiplies `e_rs^m(E) = e_rs^m|con` (the bare `1` comes from Eq. (19)'s
added linear terms with `ē_rs = e_rs^m|con`, the `a_E·rs` from Eq. (26)'s
`ẽ_rs^m|bil·rs`), so it proves nothing about Eq. (21). The old "flat rectangle"
argument is void too: on a flat rectangle `x_d = 0`, so `e_rs^m|bil = 0` and
Eq. (22) requires `ẽ_rs^m|bil = 0` — no contradiction. The implementation
implements Eq. (27a-c), which is the printed closed form and is correct; only the
motivating claim was wrong. Full record:
`openspec/changes/mitc4plusd-faithful/fidelity-audit.md` §B-audit record.

**Note F1 (2026-09-23, WU0) — Eq. (27a–c) is the closed form of Eqs. (21)+(26).**
Substituting Eq. (25) into Eq. (26) and collecting terms gives Eq. (27a–c)
exactly, and at `c_r = c_s = 0` (flat rectangle) Eq. (27) reduces term by term to
Eq. (18).  So the element implements **Eq. (27)**, not the intermediate Eq. (21).

**Eq. (23)–(25)** — the constants, with `c_r = x_d·m^r`, `c_s = x_d·m^s`,
`d = c_r² + c_s² − 1`:

```text
B1 = c_r²/d      B2 = c_s²/d      B3 = 2c_r c_s/d
B4 = −c_r/d      B5 = −c_s/d

ẽ_rs^m|bil = (c_r/d)[c_r(e_rr^m|con + e_rs^m|bil) − e_rr^m|lin]
          + (c_s/d)[c_s(e_ss^m|con + e_rs^m|bil) − e_ss^m|lin]
          + (2c_r c_s/d)·e_rs^m|con
```

**Eq. (26)** — the final assumed field:

```text
ẽ_rr^m = ê_rr^m − ẽ_rs^m|bil + ẽ_rs^m|bil·s²
ẽ_ss^m = ê_ss^m − ẽ_rs^m|bil + ẽ_rs^m|bil·r²
ẽ_rs^m = ê_rs^m + ẽ_rs^m|bil·r·s
```

**Eq. (27a–c)** — the efficient form, with
`a_A = c_r(c_r−1)/(2d)`, `a_B = c_r(c_r+1)/(2d)`, `a_C = c_s(c_s−1)/(2d)`,
`a_D = c_s(c_s+1)/(2d)`, `a_E = 2c_r c_s/d`  ← **positive**; see the sign
correction in the re-read section above.

```text
ẽ_rr^m = ½(1 − 2a_A + s + 2a_A s²) e_rr^m(A)
       + ½(1 − 2a_B − s + 2a_B s²) e_rr^m(B)
       + a_C(−1 + s²) e_ss^m(C) + a_D(−1 + s²) e_ss^m(D)
       + a_E(−1 + s²) e_rs^m(E)

ẽ_ss^m = a_A(−1 + r²) e_rr^m(A) + a_B(−1 + r²) e_rr^m(B)
       + ½(1 − 2a_C + r + 2a_C r²) e_ss^m(C)
       + ½(1 − 2a_D − r + 2a_D r²) e_ss^m(D)
       + a_E(−1 + r²) e_rs^m(E)

ẽ_rs^m = ¼(r + 4a_A·rs) e_rr^m(A) + ¼(−r + 4a_B·rs) e_rr^m(B)
       + ¼(s + 4a_C·rs) e_ss^m(C) + ¼(−s + 4a_D·rs) e_ss^m(D)
       + (1 + a_E·rs) e_rs^m(E)
```

For a flat **rectangle** `x_d = 0`, so `c_r = c_s = 0`, `d = −1` and all five
coefficients vanish: Eq. (27) reduces to Eq. (18). That is why removing them was
bit-identical on every rectangular benchmark. For a flat **distorted** element
`x_d ≠ 0` in-plane, and for a **warped** quad `x_d·n ≠ 0`; the coefficients are
then non-zero. For the twisted-beam mesh the quads are ruled, `x_d` is purely
out-of-plane, `c_r = c_s = 0`, and they vanish again — so they are **not** the
thin twisted-beam fix either.

## Transverse shear

The paper states (p. 405): *"The MITC4+ shell element uses the same assumed
transverse shear strain fields as the MITC4 shell element, but also assumed
membrane strains to also alleviate membrane locking."*

**Fig. 2(a)** of paper **B** draws that shared field, with tying points A (top,
s=1), B (bottom, s=−1), C (right, r=1), D (left, r=−1):

```text
ẽ_rt = ½(1+s)·e_rt^(A) + ½(1−s)·e_rt^(B)
ẽ_st = ½(1+r)·e_st^(C) + ½(1−r)·e_st^(D)
```

The MITC4 original is Dvorkin & Bathe (1984), Engineering Computations 1:77–88.

## Gaps in this extract

Not yet transcribed, and deliberately not guessed:

- **Eqs. (1)–(6)** — the shell kinematics and the geometry definitions feeding
  Eq. (7a). On PDF pages 2–3.
- **The assumed transverse shear construction itself** (the MITC4 tying that
  Eq. 8 of paper B shares). The page-finding from `pdftotext` put it on page 6,
  which turned out to be the benchmark section — the extraction is unreliable for
  this, so the next reader should page through 4–5 visually instead of trusting a
  text search.
- **Eqs. (21)–(27)** — transcribed above (the new MITC4+ assumed field with its
  five coefficients).
- **Tables 1–2 and the benchmark sections** — these live in paper **B**, and
  `docs/validation-matrix.md` already records the values the tests use.

Note on the two reads that failed: `pdftotext` reported Eq. (16) on PDF page 6
and Eqs. (17)–(19) on page 5. The first was wrong — page 6 is Fig. 9 and the
cylindrical shell benchmark — and the second was right. A text search is a hint
about where to look, never a substitute for looking.

## How this maps to the code

| paper | code | status |
|---|---|---|
| Fig. 4 tying points A–E | `compute_covariant_membrane_b_row` call sites in `Mitc4Precomputed::new` | **match** (A(0,1), B(0,−1), C(1,0), D(−1,0), E(0,0)) |
| Eq. (17) sampling | `compute_covariant_membrane_b_row` | **match** |
| Eqs. (18)–(19) assumed membrane | `b_m_mitc4_plus` | **match** since commit `dc7593e` |
| Eq. (9) characteristic vectors, `∂x_m/∂r = x_r + s·x_d` | `compute_j3d` via the bilinear shape derivatives | **match** (verified algebraically) |
| Eq. (7b) membrane | `compute_covariant_membrane_b_row` | **match** |
| Eq. (7c)/(7d) bending, `∂x_b` terms | `b_kappa`, `b_kappa_bubble` | **not verified** — the `∂x_b·∂u_m` warping term is not obviously present |
| Eq. (8a) `x_b` | none | **missing** — no warping enrichment in the geometry |
| Eq. (8b) `u_b` | `bubble_function` = `(1−ξ²)(1−η²)`, 2 DOFs | **diverges** — the paper's enrichment is per-node and carries the rotations `α_i`, `β_i` |
| Eq. (10) normal | `compute_local_coordinate_system` | **match** (e3 = mean of the two diagonal normals) |

The open question this extract exists to answer: our MITC4 gives **0.7313** on
the thin twisted beam at N=8 where paper **B**'s MITC4 gives **0.9959**, while our
MITC3 gives **0.9932** where paper **B**'s MITC3+ gives **0.9932** exactly. The
defect is in the quad path, and the last two rows above are the only parts of the
formulation not yet checked against the paper.

---

# Handoff: the open defect and how to attack it

**Goal.** Make our MITC4 reproduce the paper's values on the thin twisted beam, so
the two `xfail(strict)` cases in `tests/test_ko2017_performance.py` flip on their
own.

**The defect, measured.** Thin twisted beam (`t/L = 0.0002667`), in-plane load,
90° twist, **N = 8**:

| element | ours | paper B |
|---|---|---|
| MITC3 (triangles) | **0.9932** | **0.9932** (MITC3+, Table 12) — exact |
| MITC4 (quads) | **0.7313** | **0.9959** (MITC4, Table 12) |

At N=16 the quad gives 0.9131 where the paper gives 0.9975; N=4 gives 0.4207. It
converges, but far too slowly — the locking signature. The quad is right on every
flat and curved **non-warped** case (Scordelis-Lo +0.15% against the published
column, the square plates), and wrong only on the warped twisted beam. So the
defect is the **warped-quad treatment**.

**Already verified against paper A** (so do not re-check these): the Fig. 4 tying
points, Eqs. (7b), (9)–(19), (10), and the shared transverse shear field of
paper B's Fig. 2(a).

**Already refuted by measurement** (so do not re-try these): the membrane-bending
coupling in the union assembly; the frame mixing there (a real defect, fixed in
`3c7230e`, but not this); the SRI centre-shear patch (0.7307 vs 0.7313); the
assumed membrane blending (rewritten per Eqs. 18–19 in `dc7593e`, bit-identical);
the shear bubble (0.7312 vs 0.7313); and the edge strains of Eq. (14), which are
algebraically the same as Eq. (15) at the tying points.

**What is left — the only unchecked part of the formulation:**

1. **Eq. (8a), `x_b = ½ Σ a_i h_i V_n^i`** — no implementation. `x_b` is the
   element's warping enrichment and it vanishes for a flat element.
2. **Eq. (8b), `u_b = ½ Σ a_i h_i (−V_2^i α_i + V_1^i β_i)`** — our
   `bubble_function` is `(1−ξ²)(1−η²)` with **2 DOFs**; the paper's enrichment is
   **per node and carries the rotations**. This is a real divergence.
3. **Eqs. (7c)/(7d), the `∂x_b` terms** — the bending strain's warping
   contribution. `b_kappa` / `b_kappa_bubble` do not obviously have it.

Changing (2) makes `kbb` stop being 2×2, so `compute_ke_local`, `compute_fint` and
`compute_kt` all change together. The consistency guards
(`test_fint_linear_nonlinear_parity`, `test_kt_fint_directional_derivative`,
`test_assemble_kt_linear_matches_k_at_zero`) are load-bearing and will catch a
partial wiring — they already did once, during the drill work.

**Method that works.** Read equations with vision, never `pdftotext`. Validate
every change with the A/B against the paper's triangular column: the triangular
element matching to four digits is what turned "we are 8.5% off" into "the quad is
27% off and the triangle is perfect". A change that produces an identical number
is still a result — it eliminates a hypothesis.

**Machine limits.** Never run these benchmarks at N ≥ 32: the assembly is sparse
but `spsolve` does a SuperLU LU whose fill-in needs tens of GB, and a probe at N=32
was killed for memory. N ≤ 16 is ~47 s per case.

**Commands.**

```bash
# build the extension (conda env, preCICE pkg-config for the solvers crate)
source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig" HDF5_DIR="$CONDA_PREFIX"
python -m maturin develop --release

# the A/B probe: quad vs triangle on the thin twisted beam at N=8
export PATH="$HOME/miniconda3/envs/aeroelast-dev/bin:$PATH"
python /tmp/probe_tri.py          # rebuild this if /tmp was cleared

# suites
cd crates && cargo test -p aeroelast-core -q          # 112 passed
cargo test -p aeroelast-solvers --lib -- --test-threads=1   # 41 passed; parallel CRASHES (PETSc/preCICE are not thread-safe)
python -m pytest -m "not slow" -q                     # 345 passed, 2 xfailed, 2 skipped

# the papers
.sources/papers/A_new_MITC4+_shell_element.pdf        # paper A: formulation
.sources/papers/1-s2.0-S0045794917309550-main.pdf     # paper B: benchmarks
```

**Acceptance.** The twisted beam at N=8 and N=16 reaches ~0.9959 / 0.9975, the
`xfail(strict)` cases start failing as XPASS, and removing the markers leaves the
suite green.
