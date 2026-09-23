# Feature: a faithful MITC4+ (Ko, Lee & Bathe 2017), as the base for MITC4/D 2025

## Objective

Implement paper A's element as it is written, with nothing added, so that it
passes **the paper's own basic tests** — isotropy, zero-energy modes, patch tests
— and then the published benchmark columns. This element is the base on which the
2025 MITC4/D (Ko, Bathe & Zhang, *continuum mechanics-based shell elements with
six degrees of freedom at each node*) will later be built. Application target:
wind-turbine blades — flexible, curved, warped.

## Why a new element rather than a patch

Full re-read of paper A (PDF pp. 2–7 / journal pp. 405–410, with vision) on
2026-09-23. The element the repository has is **not** the paper's element, and the
difference is architectural:

| # | paper A | this repository |
| --- | --- | --- |
| 1 | **5 DOF/node**: 3 translations + 2 director rotations `α_i`, `β_i`. **No drilling rotation exists in the element** | 6 DOF/node with a drilling rotation |
| 2 | 3D continuum kinematics: `x(r,s,t) = Σh_i x_i + (t/2)Σa_i h_i V_n^i`, director `V_n^i`, thickness `a_i` | flat projection onto a local 2D frame, ABD resultants |
| 3 | the director enrichment (Eq. 8b) **is part of the kinematics** | a separate 2-DOF `(1−ξ²)(1−η²)` rotation bubble |
| 4 | full **2×2×2** Gauss integration; *"the element formulation does not include any numerical factor"* | selective reduced integration of the in-plane shear |
| 5 | membrane field **Eqs. (21)–(27)** | Eqs. (18)–(19) (`compute_membrane_coefficients` deleted) |
| 6 | no numerical factor anywhere | shear correction factor, plus a Winkler & Plakomytis ERC drilling treatment and a `β_w` warping penalty from a *different* element in a *different* paper |

Consequence, and the reason this is a new element: the repository's drilling
"problem" — the `k_drill` sensitivity, the ERC, the spurious `theta_z` modes —
**does not exist in paper A**, because paper A has no drilling DOF. It is an
artefact of our 6-DOF flat-projection architecture. Adding the drilling DOF is
exactly what the 2025 MITC4/D does, and it does it penalty-free.

Also corrected from the re-read: `a_E = +2 c_r c_s / d` (positive). The extract's
§2.3 and the deleted `compute_membrane_coefficients` both had it negative.
And `a_i` is the **thickness**, not an edge vector.

## The formulation, in one place

See `docs/formulations/mitc4plus-2017-extract.md` for the full transcription with
page citations. The essentials:

- **Geometry / kinematics (Eqs. 1–3, 5 DOF/node).**
- **Strain (Eqs. 4–5)**: covariant 3D `e_ij = ½(g_i·u_j + g_j·u_i)`,
  `r_1 = r, r_2 = s, r_3 = t`; in-plane decomposition
  `e_ij = e_ij^m + t e_ij^b1 + t² e_ij^b2` (Eq. 7a).
- **Only `e_ij^m` is modified.** The bending terms (7c)/(7d), including their
  `∂x_b` parts, stay displacement-based (paper A p. 406, in words).
- **Membrane**: sampled at the five tying points A(0,1) B(0,−1) C(1,0) D(−1,0)
  E(0,0) (Eq. 17), then the assumed field Eqs. (21)–(27) with
  `c_r = x_d·m^r`, `c_s = x_d·m^s`, `d = c_r² + c_s² − 1`,
  `a_A = c_r(c_r−1)/(2d)`, `a_B = c_r(c_r+1)/(2d)`, `a_C = c_s(c_s−1)/(2d)`,
  `a_D = c_s(c_s+1)/(2d)`, `a_E = +2c_r c_s/d`.
  Flat rectangle → all five vanish → Eq. (18) is recovered.
- **Transverse shear**: the MITC4 assumed field of Dvorkin & Bathe (1984),
  `e_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, `e_st = ½(1+r)e_st^C + ½(1−r)e_st^D`.
- **Integration**: 2×2×2, no numerical factor.

## Work units

- **F1** Module + geometry + kinematics. New `elements/mitc4_plus.rs`: the
  precomputed element data (nodal coords, directors `V_n^i`, thicknesses `a_i`,
  `x_r`, `x_s`, `x_d`, `n`, `m^r`, `m^s`, `c_r`, `c_s`, `d`, `a_A..a_E`), the
  20-DOF layout, and the shape/derivative machinery. Unit tests: flat element
  gives `x_d = 0` and all coefficients zero; `m^r·x_r = 1`, `m^s·x_s = 1`,
  `m^r·x_s = m^s·x_r = 0`, `m^{r,s}·n = 0`.
- **F2** Displacement-based operators: the membrane (Eq. 15/16), the bending
  (Eqs. 7c/7d, displacement form), the transverse shear (Eq. 4 with the MITC4
  tying), and the 2×2×2 integration with the plane-stress material. Unit tests:
  **patch test** (constant strain reproduced), **isotropy** (the element stiffness
  does not depend on the element's orientation), **zero-energy modes** (exactly
  the rigid-body modes, no more).
- **F3** The assumed membrane field (Eqs. 17–27) replacing the displacement-based
  one, with the patch-test condition of Eq. (22) as a unit test, plus a test that
  the coefficients vanish for a flat rectangle (so Eq. 27 reduces to Eq. 18).
- **F4** Wire into the repository: decide the 6-DOF assembly treatment of the
  absent drilling DOF, and document it as the single, explicit, non-paper
  ingredient. Measure the thin twisted beam and the published columns.
- **F5** Retire the superseded path once F4 passes: remove the ERC, the `β_w`
  penalty, the SRI, the 2-DOF bubble and their dead code, and rewrite
  `shell-elements.md` §2, the module header and the validation matrix to describe
  what actually runs.
- **F6** (after) MITC4/D 2025: add the drilling DOF with the penalty-free
  drill-membrane strain.

## Constraints

- Vision reads only (`pdftoppm`), never `pdftotext` for equations.
- No parameter fitting, no numerical factor that the paper does not have.
- Never run the twisted-beam benchmarks at N ≥ 32.
- The paper's basic tests (isotropy, zero-energy modes, patch tests) are the
  acceptance gate for F1–F3, before any benchmark number is consulted.
- The old element stays as a reference until the new one passes; nothing is
  deleted before then.
