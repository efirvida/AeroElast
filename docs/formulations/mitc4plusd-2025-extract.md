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

**No letter shorthand in this document.**  `mitc4plus-2017-extract.md` labels its
two 2017 papers "A" (the formulation, C&S 182) and "B" (the benchmarks, C&S 193).
That shorthand is **local to that document and is NOT reused here**: in this file
every reference names its paper, and the 2025 MITC4+/D paper is never called
"paper B".  Where a cross-reference is needed it is spelled out as
"Ko, Lee & Bathe (2017), C&S 182:404–418".

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
ẽ_rr^md = (j0/j) ũ_rr(θ) = (j0/j) h̃_m,r^l (θ_{i+1}^D − θ_i^D) x_m^l · (−x_r^l × V^D)
ẽ_ss^md = (j0/j) ũ_ss(θ) = −(j0/j) h̃_m,s^l (θ_{i+1}^D − θ_i^D) x_m^l · ( x_s^l × V^D)
ẽ_rs^md = ½ (j0/j) (ũ_r,s(θ) + ũ_s,r(θ))
        = ½ (j0/j) [ h̃_m,s^l x_m^l · (−x_r^l × V^D) − h̃_m,r^l x_m^l · (x_s^l × V^D) ] (θ_{i+1}^D − θ_i^D)
                                                                              (18)
```

**Correction (2026-09-23, WU0).** The first transcription of this equation read
`(θ_{i-1}^D − θ_i^D)`. That was a **misreading**: Eq. (16a)/(16b) on p. 7 print
`(θ_{i+1}^D − θ_i^D)` with *"the edge `l` corresponds to nodes `i` and `i + 1`"*, and
Eq. (12d) confirms it independently (`θ_5^D = (L_5/8)(θ_4^D − θ_1^D)`, edge 5 joining
nodes 1 and 4).  The dead code's `θ_{i+1} − θ_i` was therefore **right** and the
extract was wrong.  The superscript is `l` (the edge / fictitious mid-side node),
not `i`.

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

This is the **same mesh as Ko, Lee & Bathe (2017), C&S 182:404–418, Fig. 5**
(recorded in `mitc4plus-2017-extract.md`).  The mesh
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

| test | Ko, Lee & Bathe (2017), C&S 182 | Ko, Bathe & Zhang (2025), C&S 308 |
| --- | --- | --- |
| isotropy | yes (spatial) | yes (spatial) |
| zero-energy modes | exactly six, single unsupported element | yes, rigid body modes properly represented |
| patch tests | membrane, bending, shearing (Fig. 5 mesh, minimum constraints, boundary forces from constant stress states) | extension, bending, shearing, **strong form** (constant AND zero stresses), **minimum BCs**, `θ_z` free except corner B (Fig. 7 mesh and BCs) |

Because MITC4+/D = 2017 core + the 2025 drill (Eq. 26), **both** sets apply to
this repository's target element, and the 2017 set must still hold *with* the
drilling DOF present — the zero-energy count in particular must remain exactly
six.

## The full formulation, Eqs. (1)–(16) (transcribed 2026-09-23, WU0)

Read from PDF pp. 2–7 with vision. This closes the gap that the drill operator's
term-by-term verification needed.

### 2.1 Geometry and displacement interpolation (pp. 2, 4)

```text
x(r,s,t) = Σ_{i=1..4} h_i(r,s) x_i + (t/2) Σ_{i=1..4} a_i h_i(r,s) V_n^i                  (1)

h_i(r,s) = ¼ (1 + ξ_i r)(1 + η_i s),   i = 1,2,3,4
[ξ_1 ξ_2 ξ_3 ξ_4] = [ 1  −1  −1   1 ]
[η_1 η_2 η_3 η_4] = [ 1   1  −1  −1 ]                                                   (2)

u(r,s,t) = Σ h_i(r,s) u_i + (t/2) Σ a_i h_i(r,s) (θ_i × V_n^i)                        (3a)
u_i = u_x^i i_x + u_y^i i_y + u_z^i i_z                                                 (3b)
θ_i = θ_x^i i_x + θ_y^i i_y + θ_z^i i_z                                                 (3c)
V_1^i = (i_y × V_n^i)/‖i_y × V_n^i‖ ;   V_2^i = V_n^i × V_1^i                          (3d)
```

**The 2025 element's rotation DOF is the global rotation VECTOR `θ_i` with three
components `(θ_x, θ_y, θ_z)`**, not the 2017 paper's two director angles
`(α_i, β_i)`.  This is the **sixth DOF**: the third component is the drill.  It
also means the paper's DOF convention is the same as this repository's —
`(θ_x, θ_y, θ_z)` per node — so no `(α, β)` conversion is needed at the boundary.

### Strain decomposition (p. 4)

```text
e_ij = e_ij^m + t·e_ij^b1 + t²·e_ij^b2        with i, j = 1, 2                      (8a)

e_ij^m  = ½ (∂x_m/∂r_i · ∂u_m/∂r_j + ∂x_m/∂r_j · ∂u_m/∂r_i)                        (8b)
e_ij^b1 = ½ (∂x_m/∂r_i·∂u_b/∂r_j + ∂x_m/∂r_j·∂u_b/∂r_i
            + ∂x_b/∂r_i·∂u_m/∂r_j + ∂x_b/∂r_j·∂u_m/∂r_i)                             (8c)
e_ij^b2 = ½ (∂x_b/∂r_i · ∂u_b/∂r_j + ∂x_b/∂r_j · ∂u_b/∂r_i)                        (8d)

x_m = Σ h_i(r,s) x_i ,   x_b = ½ Σ a_i h_i(r,s) V_n^i                              (9a)
u_m = Σ h_i(r,s) u_i ,   u_b = ½ Σ a_i h_i(r,s) (θ_i × V_n^i)                       (9b)
```

*"The first term `e_ij^m` in Eq. (8a) is the covariant in-plane membrane strain at
the shell mid-surface (t = 0), and the remaining terms are the covariant in-plane
strains due to bending."* — i.e. the same split as the 2017 paper.

### Characteristic vectors and the plane P (p. 3)

```text
g_i = ∂x/∂r_i ,  u_i = ∂u/∂r_i      with r_1 = r, r_2 = s, r_3 = t              (4a)
g_r = x_r + s·x_d ,  u_,r = u_r + s·u_d
g_s = x_s + r·x_d ,  u_,s = u_s + r·u_d                                            (4b)
x_r = ¼ Σ ξ_i x_i ,  x_s = ¼ Σ η_i x_i ,  x_d = ¼ Σ ξ_i η_i x_i
u_r = ¼ Σ ξ_i u_i ,  u_s = ¼ Σ η_i u_i ,  u_d = ¼ Σ ξ_i η_i u_i                     (4c)

V^D = (x_r × x_s)/‖x_r × x_s‖ ,   x_r = g_r(0,0,0) ,   x_s = g_s(0,0,0)              (5)
```

`V^D` is the normal to the plane `P` at the element **centre** — one vector per
element.

### The rotation decomposition (p. 3)

```text
θ_i^D = θ_i^D V^D                              (drill rotation about V^D)           (6a)
θ_i^S = α_i V_1^i + β_i V_2^i                  (standard rotation)                  (6b)
θ_i    = θ_i^S + θ_i^D                                                               (6c)
θ_i    = α_i V_1^i + β_i V_2^i + γ_i V_n^i                                           (6d)
δγ_i   = δθ_i^D (V_n^i · V^D)                                                        (6e)
θ_i^D  = α_i(V_1^i·V^D) + β_i(V_2^i·V^D) + γ_i(V_n^i·V^D)                            (6f)

e_ij = ½ (g_i · u_j + g_j · u_i)                                                     (7)
```

The paper states that *"the components in Eq. (6d) give directly the global
components in Eq. (3c), and vice versa"*, and that `(α_i, β_i, γ_i)` are used to
update the directors `(V_1^i, V_2^i, V_n^i)`.

### 2.2 New interpolation functions (p. 4)

*"In the formulation to include the drill rotation we create just a **single
additional interpolation function (along the edges)** to have stability of the
element for in-plane rotations, pass the patch test when the rotation is free
(without use of a special procedure) and improve the element behavior."*

```text
[h_5 h_6 h_7 h_8] = [ ½(1−s²)(1+r) , ½(1−r²)(1+s) , ½(1−s²)(1−r) , ½(1−r²)(1−s) ]    (10)
```

**Eq. (10) is byte-identical to the repository's `midside_shape_function`.**

```text
[h̃_5,r h̃_6,r h̃_7,r h̃_8,r] = [ 0 , ½(−2r)(1+s) , 0 , ½(−2r)(1−s) ] ,  h̃'_m,r = h̃_l,r   (11a)
[h̃_5,s h̃_6,s h̃_7,s h̃_8,s] = [ ½(−2s)(1+r) , 0 , ½(−2s)(1−r) , 0 ] ,  h̃'_m,s = h̃_l,s   (11b)
```

**Eqs. (11a)/(11b) are byte-identical to the repository's
`drill_midside_shape_derivatives`** — including the deliberate zeros.  The paper's
justification is explicit: *"the zeros in Eq. (11b) are introduced to avoid a
higher order of numerical integration than used with the original elements.  The
important point is that this simplification can be used with the basic and patch
tests still satisfied, and we also see good convergence behavior."*  This is the
"curl" simplification the code's comment calls it.

### The fictitious mid-side nodes (p. 6)

```text
u_t^5(l) = (1 − l/L_5) u_t^4 + (l/L_5) u_t^1                                       (12a)
u_n^5(l) = (1 − l/L_5) u_n^4 + (l/L_5) u_n^1 + (4l/L_5)(1 − l/L_5) θ_5^D            (12b)

∂u_n^5/∂l |_{l=L_5} − ∂u_n^5/∂l |_{l=0} = −θ_4 + θ_1                               (12c)

θ_5^D = (L_5/8)(θ_4 − θ_1) = (L_5/8)(θ_4^D − θ_1^D)                                 (12d)
```

**Edge 5 joins nodes 1 and 4** (the paper says so in words: *"the normal
displacement at node 5 (along the edge of node 5 between nodes 1 and 4)"*).
Eq. (12d) is therefore the authoritative statement of the edge convention:
**`θ_m^D = (L_m/8)(θ_b^D − θ_a^D)` for the edge joining corner nodes `a` and
`b`**, which for edge 5 is `θ_4^D − θ_1^D`.  That is the same ordering as
Eq. (16a)/(16b)'s `(θ_{i+1}^D − θ_i^D)` and it is what settles the WU0
convention question: the **`θ_{i+1} − θ_i` order is correct**.

### Drill kinematics (pp. 6–7)

```text
u_r^l m^r + u_s^l m^s = u_n^l ĩ_n + u_t^l ĩ_t                                    (13a)
ĩ_t = x_m^l/‖x_m^l‖ ,   ĩ_n = ĩ_t × V^D                                            (13b)
[x_m^5 x_m^6 x_m^7 x_m^8] = [ ⅛(x_4 − x_1) , ⅛(x_1 − x_2) , ⅛(x_2 − x_3) , ⅛(x_3 − x_4) ]   (13c)

u_r^l = u_n^l ĩ_n·x_r^l + u_t^l ĩ_t·x_r^l ,   u_s^l = u_n^l ĩ_n·x_s^l + u_t^l ĩ_t·x_s^l   (14a)
u_t^l = u_θ^l(l) ĩ_n·x_r^l ,   u_s^l = u_θ^l(l) ĩ_n·x_s^l                         (14b)
u_θ^5(l) = (4l/L_5)(1 − l/L_5) θ_5^D ,   similarly for edges 6, 7, 8               (14c)
u_r^l(l) = (1/‖x_m^l‖)( u_θ^l(l) x_m^l × V^D · x_r^l )                              (14d)
u_s^l(l) = (1/‖x_m^l‖)( u_θ^l(l) x_m^l × V^D · x_s^l )                              (14e)

u_r = h̃_m^l u_r^l(L_l/2) = (1/‖x_m^l‖) [ −h̃_m^l u_θ^l(L_l/2) x_m^l × V^D ] · x_?^l   (15a)
u_s = h̃_m^l u_s^l(L_l/2) = (1/‖x_m^l‖) [ −h̃_m^l u_θ^l(L_l/2) x_m^l × V^D ] · x_?^l   (15b)
```

*"in which we use `l = L_l/2` for each edge because the tying is performed at the
fictitious nodes and also `‖x_m^l‖ = L_l/8`."*  (Consistent with Eq. (13c):
`‖x_m^5‖ = ⅛‖x_4 − x_1‖ = ⅛ L_5`.)

**Open transcription detail:** the trailing dot product in Eqs. (15a)/(15b) prints
as `x_m^l` in both, which cannot be right — by analogy with (14d)/(14e) it should
be `x_r^l` in (15a) and `x_s^l` in (15b).  Recorded as an ambiguity rather than
silently corrected; the implementation must resolve it against the strain
symmetry.

### The displacement fields the strain is built from (p. 7)

```text
ū_r(θ) =  h_l (θ_{i+1}^D − θ_i^D) x_m^l · (−x_r^l × V^D)
ū_s(θ) = −h_l (θ_{i+1}^D − θ_i^D) x_m^l · ( x_s^l × V^D)                            (16a)

ũ_r(θ) =  h̃_m^l (θ_{i+1}^D − θ_i^D) x_m^l · (−x_r^l × V^D)
ũ_s(θ) = −h̃_m^l (θ_{i+1}^D − θ_i^D) x_m^l · ( x_s^l × V^D)                          (16b)
```

*"where the edge `l` corresponds to nodes `i` and `i + 1`.  Eq. (16a) is the
displacement-based interpolation of the drill effect with the standard mid-side
interpolation `h_m^l = h_l` given in Eq. (10).  Next use the 'assumed
interpolation' at the edges.  By applying Eq. (11) to Eq. (16a), we obtain"*
Eq. (16b).  Substituting (16b) into (17b) then gives the strain, Eq. (18).

**A key statement, p. 6, verbatim:** *"we utilize the vector `V^D` at the element
center about which to measure all rotations at the nodal points (that is, we use
`θ_i = θ_i^D` for all nodes `i = 1, 2, ..., 8`, hence all nodal drill rotations are
measured as drill rotations about the same vector `V^D`).  It then follows that
continuity between edge displacements is always obtained, and the rigid body mode
test is always satisfied even when the shell element is curved."*

### What this settles

| question | answer, with the source |
| --- | --- |
| DOF convention | the global rotation vector `θ_i = (θ_x, θ_y, θ_z)` per node — Eq. (3c). 6 DOF/node, same as this repository |
| the drill component | `θ_i^D = θ_i · V^D` — Eq. (6a)/(6f) |
| the mid-side functions | Eq. (10), byte-identical to the repo's `midside_shape_function` |
| the "curl" derivatives | Eqs. (11a)/(11b), byte-identical to the repo's `drill_midside_shape_derivatives`, zeros included and justified by the paper |
| the edge ordering | `(θ_{i+1}^D − θ_i^D)`, edge `l` joining nodes `i` and `i+1` — Eq. (16a)/(16b), confirmed by Eq. (12d) |
| `V^D` | one vector per element, the plane-P normal at the centre — Eq. (5) |
| `‖x_m^l‖` | `L_l/8` — Eq. (13c) and the p. 7 text |

## The drill-membrane strain as a matrix, and the zero-energy-mode question
### Eqs. (19)–(25), pp. 10 and 12 (transcribed 2026-09-23, from the WU6 failure)

This section was written to answer a question WU6 raised: the element has EIGHT
zero modes on a flat element and SEVEN on the distorted and warped ones, where
the paper says exactly six.  Reading Eqs. (19)–(25) answers it, and the answer is
not the one the design expected.

**Eq. (19a)** — the drill-membrane strain as a matrix times the four corner drill
rotations:

```text
[ẽ_rr^md  ẽ_ss^md  ẽ_rs^md]^T  =  B̃ [θ_1^D  θ_2^D  θ_3^D  θ_4^D]^T ,   B̃ = (j0/j) [B̃_rr ; B̃_ss ; B̃_rs]     (19a)
```

**Eq. (19b)** — the full form, with the edge index `l` and the edge coefficients
`c_r^l`, `c_s^l`:

```text
B̃_rr = [ h̃_m,r^5 c_r^5 − h̃_m,r^6 c_r^6 ,  h̃_m,r^6 c_r^6 − h̃_m,r^7 c_r^7 ,
          h̃_m,r^7 c_r^7 − h̃_m,r^8 c_r^8 ,  h̃_m,r^8 c_r^8 − h̃_m,r^5 c_r^5 ]

B̃_ss = [ −h̃_m,s^5 c_s^5 + h̃_m,s^6 c_s^6 ,  −h̃_m,s^6 c_s^6 + h̃_m,s^7 c_s^7 ,
          −h̃_m,s^7 c_s^7 + h̃_m,s^8 c_s^8 ,  −h̃_m,s^8 c_s^8 + h̃_m,s^5 c_s^5 ]

B̃_rs = ½ ( [ −h̃_m,r^5 c_s^5 + h̃_m,r^6 c_s^6 ,  −h̃_m,r^6 c_s^6 + h̃_m,r^7 c_s^7 ,
              −h̃_m,r^7 c_s^7 + h̃_m,r^8 c_s^8 ,  −h̃_m,r^8 c_s^8 + h̃_m,r^5 c_s^5 ]
          + [  h̃_m,s^5 c_r^5 − h̃_m,s^6 c_r^6 ,   h̃_m,s^6 c_r^6 − h̃_m,s^7 c_r^7 ,
               h̃_m,s^7 c_r^7 − h̃_m,s^8 c_r^8 ,   h̃_m,s^8 c_r^8 − h̃_m,s^5 c_r^5 ] )   (19b)
```

**Eq. (19c)** — the edge coefficients, which are the 2025 quantity (NOT the 2017
`c_r = x_d·m^r`):

```text
c_r^l = x_m^l · (−x_r^l × V^D) ,      c_s^l = x_m^l · (x_s^l × V^D)                    (19c)
```

**Eq. (19d)** is the same `B̃` in a "simpler form" using Eq. (11)'s curl
derivatives, in which `h̃_m,r^5 = h̃_m,r^7 = 0` and `h̃_m,s^6 = h̃_m,s^8 = 0`.

**Eq. (20)** — the MITC4 assumed transverse shear, as in the 2017 paper:
`ẽ_rt = ½(1+s)e_rt^(A) + ½(1−s)e_rt^(B)`, `ẽ_st = ½(1+r)e_st^(C) + ½(1−r)e_st^(D)`.

**Eq. (21)** — the drill-membrane strain transformed to the natural coordinate
system: `e_ij^md = (g_i·g^k)(g_j·g^l) ẽ_kl^md`.  The paper notes that `g_i` and
`g^i` are evaluated at the element **centre** and are constant over the element,
and that *"the strain in Eq. (21) is calculated only once through the thickness t,
and the calculation of the coefficients `c_r^l`, `c_s^l` and `(g_i·g^k)` is only
performed once per each element."*

**Eq. (22a)** — the total strain of the element the paper calls **MITC4/D**:

```text
e_ij = e_ij^m + e_ij^md + t·e_ij^b1 + t²·e_ij^b2      with i, j = 1, 2                 (22a)
```

which is the same structure as Eq. (26), the **MITC4+/D**.  Eqs. (22b)/(22c)
relate the "standard" displacement-based strain to the natural-coordinate one
through the centre base vectors.  Eq. (23a) repeats the 2017 decomposition
(`e_rr^m = e_rr|con + e_rr|lin·s + e_rs|bil·s²`, etc.).

**Eqs. (24)–(25)** are **not** a drill term.  They are the derivation of the 2017
membrane coefficients by the transformation chain
`[a_0|rr … a_4|rs] = Q(r,s) R(r,s) S`, with `Q`, `R`, `S` printed in (25b)–(25c)
and the `n_i`, `m_i` in (25d).  The paper says the coefficients "require only few
bases to be calculated, and the computation is performed only once per element",
and that Eq. (24) "include the improvements proposed in Ref. [8] for the MITC4+
element".  **There is no additional term that penalizes or constrains a constant
drill rotation.**

### An inconsistency in the paper's own reduction, Eqs. (14)–(18) (found 2026-09-24)

Read from p. 7 with vision, to settle why the two THICK twisted-beam cells come out
2–3x too soft while the three thin ones match.

**The paper's derivation carries a per-edge `1/‖x_m^I‖`:**

```text
u_r^I(l) = (1/‖x_m^I‖)(u_θ^I(l) x_m^I × V^D · x_r^I)                        (14d)
u_s^I(l) = (1/‖x_m^I‖)(u_θ^I(l) x_m^I × V^D · x_s^I)                        (14e)
u_r = h_m^I u_r^I(L_I/2) = (1/‖x_m^I‖)[ −h_m^I u_θ^I(L_I/2) x_r^I × V^D ]·x_m^I     (15a)
u_s = h_m^I u_s^I(L_I/2) = (1/‖x_m^I‖)[ −h_m^I u_θ^I(L_I/2) x_s^I × V^D ]·x_m^I     (15b)
```

and it states in words, on the same page: *"in which we use `l = L_I/2` for each edge
because the tying is performed at the fictitious nodes and also `‖x_m^I‖ = L_I/8`."*
So the factor is `8/L_I`.

**But Eq. (16a) — which the paper presents as the reduction of Eq. (15), *"Using the
geometric relations (Eq. (12d) and (13c)), the displacement fields assumed in Eq. (15)
reduce to"* — drops it:**

```text
ū_r(θ) =  h_I (θ_{i+1}^D − θ_i^D) x_m^I · (−x_r^I × V^D)
ū_s(θ) = −h_I (θ_{i+1}^D − θ_i^D) x_m^I · ( x_s^I × V^D)                     (16a)
```

Eq. (16b) keeps the same shape with `h̃_m^I` in place of `h_I`, and Eq. (18) — the
operative strain — inherits it.  **Eq. (16a) is dimensionally inconsistent without
that factor**: `x_m^I·(x_r^I × V^D)` already carries `L²`, and the result is a
displacement, so the `θ·L²` needs the `1/L` that `1/‖x_m^I‖` supplies (`θ·L²/L = θ·L`).

**Consequence for this repository, measured:** the implementation follows Eq. (18), so
its drill block `B̃ᵀ C B̃` is `(8/L_I)²` too small.  The earlier edge-convention note
(above) and this one are the paper's **second** internal inconsistency; the first is
the missing leading term in the 2017 paper's Eq. (21) (note F2 of that extract).  Both
were found by measuring, not by trusting the equations.

**But adopting the Eq. (15) form is not a fix by itself — it is a trade.** Measured on
the five twisted-beam cells (this repository's element, `h = 0.32` for the thick cases):

| case | published | with Eq. (18) | with Eq. (15)'s factor |
| --- | --- | --- | --- |
| thin N=8 in | 0.9959 | 0.9958 | 0.9732 |
| thin N=16 in | 0.9975 | 0.9978 | 0.9728 |
| thin N=16 out | 0.9980 | 0.9987 | 0.9684 |
| thick N=16 in | 0.9972 | **1.9623** | **1.0141** |
| thick N=16 out | 0.9972 | **2.9933** | **0.9028** |

So the factor repairs the thick in-plane cell but overshoots the thick out-of-plane by
9.5% and pushes all three thin cells 2.3–3.1% below the published columns — the same
thin-wants-soft / thick-wants-stiff conflict, reversed.  **At least one further error
remains**, and it is not a drill *scale*: neither paper has one, and none may be added.
Settling it needs a deeper pass over the drill's cross terms, the Eq. (21) centre
metric and the through-thickness integration than this session had.

### What this settles, and what it breaks

**Every column of `B̃` is a DIFFERENCE of edge terms** — `B̃_rr`'s columns are
`(5→6)`, `(6→7)`, `(7→8)`, `(8→5)`, so they sum to zero, and the same telescoping
holds for `B̃_ss` and `B̃_rs`.  Therefore:

> **A constant drill rotation `θ^D` produces zero drill-membrane strain by
> construction of the paper's own Eq. (19).**

And the 2017 core is blind to the drill direction, so a field with `u = 0` and
`θ_z` constant carries zero energy in **every** term.  That is exactly the surplus
zero mode WU6 measured, and it is **the paper's behaviour, not an implementation
defect**.  It also explains why the paper's own patch tests (Fig. 7(b)(c)(d))
leave `θ_z` free at every node **except the corner node B**: constraining it at one
node is precisely what removes the constant mode.

Consequences that must be carried forward:

1. **The constant-drill zero mode is by design.**  Requirement 5's "exactly six
   zero eigenvalues" is therefore **not reconcilable with Eq. (19) as printed**,
   unless the drill is constrained (as the patch tests do) or the requirement is
   restated as "exactly six rigid-body modes plus the drill operator's own null
   space".  This is a **spec defect**, not an element defect.
2. **The second surplus mode on the flat rectangle** follows from Eq. (19d)'s
   curl zeros: with `h̃_m,r^5 = h̃_m,r^7 = 0`, `B̃_rr` has rank 2 instead of 4, so
   the drill block loses two directions on a flat element and one on a warped one
   — which matches the measured 8 and 7.
3. The paper's Section 3.1 statement that the elements pass the zero-energy-mode
   test with rigid body modes "properly represented" is therefore either
   qualified (the drill constrained) or loose.  **It cannot be reproduced from
   Eqs. (19)–(25) as printed**, and this extract records that rather than
   pretending the equations agree.

## Gaps in this extract

- ~~**Eqs. (1)-(16)**~~ — **transcribed** in the section above (WU0, 2026-09-23).
  `transcription-verified: 2026-09-23 pp.2-7`
- **Eqs. (19)-(25)** and the full `θ_z` interpolation are not yet transcribed.
  `transcription-verified: 2026-09-23 p.8` for Eq. (18) only.
- The **load magnitudes** of the three patch tests (see above).
- The **convergence-study norm** used in Section 3 (the paper reuses the
  s-norm of Hiller & Bathe from the 2017 paper, Eq. 28 there).

## WU0 result: the transverse-shear metric question (task 1.4)

Task 1.4 asked whether Ko, Lee & Bathe (2017) prints its own covariant
transverse-shear metric normalization.  **It does, as Eq. (4)**:
`e_ij = ½(g_i·u_j + g_j·u_i)` with `r_1 = r, r_2 = s, r_3 = t`, so the covariant
transverse shear is `e_rt = ½(g_r·u_t + g_t·u_r)` and likewise for `e_st`.  What
that paper does **not** print is the construction of the tying-point values
`e_rt^(A)` etc. beyond the interpolation itself (`ẽ_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, p. 405), which it
defers to **Dvorkin & Bathe (1984), Engineering Computations 1:77-88** — not held
in `.sources/papers/`.  So the design's §2.2 definition stands for the tying-point
metric, and this is the recorded answer: printed for the covariant definition,
not printed for the tying-point construction.
