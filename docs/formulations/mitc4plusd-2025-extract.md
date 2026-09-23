# MITC4/D and MITC4+/D extract — Ko, Bathe & Zhang (2025), C&S 308:107622

A persistent transcription of the equations, figures and tests this repository's
MITC4+/D implementation is checked against, so that verifying the code does not
require re-reading the PDF each time.

**This is an extract, not a substitute for the paper.** Everything below was read
from the PDF with vision (`pdftoppm -png -r 300`, never `pdftotext` for
equations). Gaps are marked as gaps.

- Paper: Ko, Y., Bathe, K.-J., Zhang, L., *"Continuum mechanics-based shell
  elements with six degrees of freedom at each node − the MITC4/D and MITC4+/D
  elements"*, **Computers and Structures 308 (2025) 107622**.
- File: `.sources/papers/1-s2.0-S0045794924003511-main.pdf` (21 pages).
- Companion/duplicate: `.sources/papers/MITC_D_elements_published.pdf`.

## What this element is, and its relation to the 2017 MITC4+

The 2017 element (Ko, Lee & Bathe, C&S 182:404-418 — see
`mitc4plus-2017-extract.md`) has **5 DOF/node** and no drilling rotation. This
2025 element adds the **sixth DOF (the drill rotation `θ_z`)** on top of the same
formulation, penalty-free. Eq. (26) states the total strain:

```text
e_ij = ẽ_ij^m + e_ij^md + t·e_ij^b1 + t²·e_ij^b2      with i, j = 1, 2        (26)
```

i.e. **the assumed membrane (2017) + the drill-membrane strain (new) + the
displacement-based bending terms (2017)**. So MITC4+/D = 2017 MITC4+ plus
`e_ij^md`, and the 2017 paper's own basic tests remain the core's acceptance
criteria.

The paper's own words on the interpolation cost (p. 13): *"the above
interpolation used for the MITC4 + shell element might appear rather complicated
when in fact it is reasonably straightforward in application."*

## The drill-membrane strain — Eq. (18), p. 8

Context, verbatim: *"with the indices i and j = 1,2 corresponding to r_1 = r and
r_2 = s. We also use `j = det[g_r g_s g_t]|_{(r,s,0)}` and `j0 = j(0,0,0)`. The
strain in Eq. (17a) is the usual displacement-based strain and the strain in
Eq. (17b) is an approximation that is simple, relative to the value used in
Ref. [20], and when used leads to satisfying the patch tests."*

```text
ẽ_ij(θ) = ½ (ū_i,j + ū_j,i)                                          (17a)
ẽ_ij(θ) = (j0/j) ½ (ũ_i,j + ũ_j,i)                                   (17b)
```

Substituting Eq. (16b) into Eq. (17b) gives the individual "drill-membrane"
components (superscript `md`):

```text
ẽ_rr^md = (j0/j) ũ_rr(θ) = (j0/j) h̃_m,r^i (θ_{i-1}^D − θ_i^D) x_m^i · (−x_r^i × V^D)
ẽ_ss^md = (j0/j) ũ_ss(θ) = −(j0/j) h̃_m,s^i (θ_{i-1}^D − θ_i^D) x_m^i · ( x_s^i × V^D)
ẽ_rs^md = ½ (j0/j) (ũ_r,s(θ) + ũ_s,r(θ))
        = ½ (j0/j) [ h̃_m,s^i x_m^i · (−x_r^i × V^D) − h̃_m,r^i x_m^i · (x_s^i × V^D) ] (θ_{i-1}^D − θ_i^D)
                                                                              (18)
```

**Two facts that matter for the implementation** (p. 8, in words):

> *"using the rotations about `V^D`, leads to satisfying the element isotropy test
> and rigid body modes test, and even when the element is not flat (Fig. 1(c) and
> 2(b)). In addition, it is valuable to note that **using the single normal to the
> plane P results in an overall improved element behavior** when compared to using
> a changing normal over the element geometry."*

So `V^D` is **one vector per element** — the normal to the plane `P` at the
element centre — not a per-node normal.  And the strain carries the ratio `j0/j`.

### Correspondence with the code

`b_md_mitc4_plus` and `drill_midside_shape_derivatives` in
`crates/aeroelast-core/src/elements/mitc4.rs` cite 2025 Eqs. (5), (10), (11a-b),
(13b-c), (17a-b), (18), (19a-c) and (21), and their doc comment summarises the
operator as

```text
coeff_rr = (j0/j) h̃_r c_r
coeff_ss = −(j0/j) h̃_s c_s
coeff_rs = ½ (j0/j) (h̃_s c_r − h̃_r c_s)
```

which matches the structure of Eq. (18) with `c_r ≡ x_m^i · (−x_r^i × V^D)` and
`c_s ≡ x_m^i · (x_s^i × V^D)`.  **Note that these are NOT the 2017 membrane
coefficients of the same name** (`c_r = x_d·m^r`, `c_s = x_d·m^s`): same symbol,
different quantity, different paper.  Term-by-term verification against Eq. (18)
is still owed — the operator has only ever been dead code.

## Section 3.1 — "Basic tests including the patch tests" (pp. 13-14)

Verbatim, p. 13:

> *"MITC/D and MITC4+/D elements pass the spatial isotropy test. The elements pass
> the zero energy mode tests, and rigid body modes are properly represented."*

And p. 13 on integration: *"the elements require only the **2 × 2 Gauss
integration over the element surfaces**."*

Verbatim, p. 14 — the patch tests:

> *"The patch tests are shown in **Fig. 7**. We consider the **'strong form' of the
> patch tests** [5,6,22,23], i.e. we require that the calculations give the
> **analytical solutions of constant and zero stresses throughout the patch** due
> to the applied loading. When considering the patch tests, it is important that
> the **'minimum boundary conditions'** are imposed to constrain the patch from
> undergoing a rigid body motion. For each test, the solution should then be the
> constant straining mode."*
>
> *"Fig. 7 shows in detail what boundary conditions are imposed when including the
> drill degree of freedom `θ_z` at the nodes for the extension, bending, and
> shearing tests. Note that in all cases **`θ_z` is left free at the element nodes
> except at the corner node B**."*
>
> *"It is important that this strict form of the patch tests is performed, and for
> the elements presented all tests are strictly passed. Also, the use of `θ_z = 0`
> at the corner node C does not affect the results."*

### Fig. 7(a) — the mesh (read from the figure; p. 5)

**It is not a 3 × 3 grid.**  It is a **five-element "star" patch with 8 nodes**,
and every node coordinate is labelled on the figure:

```text
corners  B(0, 10)   A(10, 10)   D(10, 0)   C(0, 0)
interior          (4, 7)   (8, 7)   (8, 3)   (2, 2)

central element : (2,2) → (4,7) → (8,7) → (8,3)
surrounding     : B-A-(8,7)-(4,7)     A-D-(8,3)-(8,7)
                  D-C-(2,2)-(8,3)     C-B-(4,7)-(2,2)
```

This is the **same mesh as paper A's Fig. 5** (the 2017 extract).  The mesh
generation is therefore fully determined and needs no equivalent-patch deviation.

### Fig. 7(b)(c)(d) — the boundary conditions (read from the figure; p. 5)

Each figure is a schematic of the same square; the constrained nodes and their
constrained components are printed in boxes.

```text
(b) Bending
    B : u_x = u_y = u_z = 0 ,  θ_x = θ_y = θ_z = 0      (fully clamped)
    C : u_x = u_z = 0 ,        θ_y = 0
    load at A in the +y direction ;  θ_z free

(c) Shearing
    B : u_x = u_y = u_z = 0 ,  θ_x = θ_y = θ_z = 0
    C : u_x = u_z = 0 ,        θ_x = θ_y = 0
    interior : u_x = 0 ,       θ_x = 0 , θ_y = 0
    load at A in the +y direction ;  θ_z free

(d) Extension
    B : u_x = u_y = u_z = 0 ,  θ_x = θ_y = θ_z = 0
    C : u_x = u_z = 0
    load at A in the +x direction ;  θ_z free
```

Not published in the figure: the **magnitude** of the applied nodal loads.  They
are implied by the requirement that the solution be the constant straining mode,
so a test must either scale them from the constant stress state or assert the
constant-stress recovery rather than a displacement value.

## 2017 vs 2025: which tests apply to what

| test | 2017 (paper A) | 2025 (paper B) |
| --- | --- | --- |
| isotropy | yes (spatial) | yes (spatial) |
| zero-energy modes | exactly six, single unsupported element | yes, rigid body modes properly represented |
| patch tests | membrane, bending, shearing (Fig. 5 mesh, minimum constraints, boundary forces from constant stress states) | extension, bending, shearing, **strong form** (constant AND zero stresses), **minimum BCs**, `θ_z` free except corner B (Fig. 7 mesh and BCs) |

Because MITC4+/D = 2017 core + the 2025 drill (Eq. 26), **both** sets apply to
this repository's target element, and the 2017 set must still hold *with* the
drilling DOF present — the zero-energy count in particular must remain exactly
six.

## Gaps in this extract

- **Eqs. (1)-(16)** of the 2025 paper — the kinematics, the geometry and the
  displacement interpolation feeding Eqs. (17)/(18) — are **not yet transcribed**.
  They live on PDF pages 2-8 and are needed for the term-by-term verification of
  the drill operator and for the `h̃_m` mid-side functions.
- **Eqs. (19)-(25)** and the full `θ_z` interpolation are not yet transcribed.
- The **load magnitudes** of the three patch tests (see above).
- The **convergence-study norm** used in Section 3 (the paper reuses the
  s-norm of Hiller & Bathe from the 2017 paper, Eq. 28 there).
