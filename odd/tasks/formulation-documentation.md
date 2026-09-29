# Feature: Formulation documentation and the MITC4/D + MITC3 corrections

## Objective

Make every formulation in this repository traceable to its source: a definitive,
well-referenced document of the formulations we implement, a canonical
bibliography, and tests whose reference values and tolerances are documented
against the literature. Along the way, correct the two formulation defects this
work has already exposed (the MITC3 rotation sign, and the unwired MITC4/D drill
operator), always against the papers.

## Decisions taken (user, 2026-09-23)

1. Authorised to find and download the referenced papers into a **gitignored**
   folder.
2. **MITC4/D drill operator**: restore and wire it **only if it significantly
   improves the element's behaviour**, always checking it against the paper's
   formulation first.
3. **MITC3 sign defect: fix it.** Also search for the most recent MITC3 papers
   and implement the most current version, keeping the references.
4. Create **one definitive, well-referenced document of the implemented
   formulations, in English**. The existing docs are the user's to correct later.
5. Replace the ANSYS citations with the **original source, verified before
   replacing**.
6. Priority order: **(i) element kernels, (ii) materials, (iv) tests**.
7. The validation matrix is wanted, but the better shape is: locate the reference
   papers, document the validation properly, improve the tests, and have the
   single reference table documented and referenced to the literature.

## Recovered source material

The reference PDFs were committed and later deleted; the README even said so
("The reference PDFs used during development were removed from the repository").
They are still in git history and were recovered into **`.sources/papers/`**
(already covered by the `.sources` entry in `.gitignore`) with
`git -c core.quotepath=false rev-list --all --objects | grep '\.pdf$'`.

| File in `.sources/papers/` | Reference |
| --- | --- |
| `1-s2.0-S0045794924003511-main.pdf` | Ko, Y., Bathe, K.-J., Zhang, X., "Continuum mechanics-based shell elements with six degrees of freedom at each node - the MITC4/D and MITC4+/D elements", *Computers & Structures*, 308:107622, 2025. doi:10.1016/j.compstruc.2024.107622 |
| `MITC_D_elements_published.pdf` | Same paper, scanned copy without a text layer (duplicate; use the one above) |
| `1-s2.0-S0045794916309464-main.pdf`, `A_new_MITC4+_shell_element.pdf` | Ko, Y., Lee, P.-S., Bathe, K.-J., "A new MITC4+ shell element", *Computers & Structures*, 182:404-418, 2017 |
| `The_MITC3+_shell_element_and_its_performance.pdf` | Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element and its performance", *Computers & Structures*, 138:12-23, 2014 |
| `mitc3+_no_lineal.pdf` | Jeon, H.-M., Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element in geometric nonlinear analysis", *Computers & Structures*, 146:91-104, 2015 |
| `mitc4+_no_lineal.pdf.pdf`, `The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf` | Ko, Y., Lee, P.-S., Bathe, K.-J., "The MITC4+ shell element in geometric nonlinear analysis", *Computers & Structures*, 185:1-14, 2017 |
| `1-s2.0-S0045794917309550-main.pdf`, `Performance_of_the_MITC3+_and_MITC4+_shell_elements_in_widely_used_benchmark_problems.pdf` | Ko, Y., Lee, Y., Lee, P.-S., Bathe, K.-J., "Performance of the MITC3+ and MITC4+ shell elements in widely-used benchmark problems", *Computers & Structures*, 193:187-206, 2017 |
| `A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf` | Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node shell element for general nonlinear analysis", *Engineering Computations*, 1(1):77-88, 1984 - the original MITC4 |
| `1-s2.0-S0045794922001936-main.pdf` | Choi, H.-G., Lee, P.-S., "Towards improving the 2D-MITC4 element for analysis of plane stress and strain problems", *Computers & Structures*, 275:106933, 2023 |
| `1-s2.0-S0045794923002377-main.pdf` | Cui, X., Peng, G., Ran, Q., Zhang, H., Li, S., "Derivation and implementation of one-point quadrature quadrilateral shell element with MITC4+ method (MITC4+R)", *Computers & Structures*, 291:107207, 2024 |
| `1-s2.0-S0045794917317078-main.pdf` | Katili, I., Batoz, J.-L., Maknun, I.J., Lardeur, P., "A comparative formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements with new numerical results based on s-norm tests", *Computers & Structures*, 204:48-64, 2018 |
| `32452534.pdf` | Wong, F.T., "Testing of Shell Elements using Challenging Benchmark Problems", 2nd Indonesian Structural Engineering and Materials Symposium |
| `U0707194204.pdf` | American Journal of Engineering Research, 7(7):194-204, 2018 - identify on read |
| `Implementación del MITC4_MITC4+ en MITC4.py vs formulaciones teóricas.pdf` | Internal Spanish comparison note (not a paper) |

Reading them is `pdftotext <file> -` (poppler is installed; no Python PDF library
is available in the conda env).

Still to obtain: Reddy (2004) for CLT, Tsai & Wu (1971) and Hashin (1980) for the
failure criteria, Newmark (1959) for the time integration, and the most recent
MITC3 work (directive 3).

## Work units

### U1 - Canonical bibliography

- New `docs/references.md`: every reference with full author list, journal,
  volume, pages, year and DOI, grouped by subject (element formulations,
  constitutive/CLT, failure criteria, time integration, FSI coupling, BEM,
  software).
- Fix the README's "Element formulation references": three entries currently
  carry a title only (the 2018 Katili, the 2023 Choi & Lee and the 2025 Ko/Bathe/
  Zhang papers); there are no DOIs anywhere; the linear MITC3+ paper that
  `mitc3.rs` cites is missing; and two entries omit Bathe while the code cites him.
- Make the README link to `docs/references.md` instead of duplicating the list.

### U2 - Shell element formulations (i)

- New `docs/formulations/shell-elements.md`, English, definitive: MITC3+, MITC4+,
  MITC4/D and the plane quads, each with the formulation as implemented, and a
  **code -> paper -> equation map** (file, function, line, paper, equation).
- Covers the current gap that `quad.rs` has no external reference at all (it says
  "matching the Python reference", which is circular) and that `mitc3.rs:470,476`
  cites "Eq. 15"/"Eq. 16" without naming the paper in the comment.
- Must state the rotation/director convention explicitly, since that is what the
  MITC3 sign defect turns on.

### U3 - Constitutive and material formulations (ii)

- Same document or a sibling: CLT/ABD, orthotropic, composite shell constitutive,
  and the failure criteria. `materials/failure.rs` already has a good `# References`
  block (Tsai & Wu 1971, Hashin 1980); the others need the same treatment.

### U4 - Solver formulations

- Newmark-beta (`solvers/elasticity/dynamic_newmark.py` has equations and no
  citation), the co-rotational/rotating-frame formulation
  (`solvers/fsi/corotational.py`, `rotor.py`), the updated-Lagrangian path, the
  SRI scheme (`mitc4.rs:946-953`, already cited to Hughes, Taylor &
  Kanoknukulchai 1977) and the assumed-strain/MITC machinery.

### U5 - MITC3 sign fix (DONE, commit `d6f37fb`) and modernisation (open)

**Done.** `eval_covariant_shear_ext` restored to the physical director, verified
byte-identical to `b136ce5^` for all eight rotation assignments, with the
derivation from Ko, Bathe & Zhang 2025 Eq. (3a) in the comment. The five MITC3
tests that encoded the inverted convention were corrected (`_analytical_tip`
negates its transverse component, `REFERENCE_TABLE`'s four w_ref become
-6.3662, -6.3662, -2.1221, 0.0000, `test_linear_tip_moment_sign` asserts
`w_tip < 0`), and the now-obsolete `xfail(strict=True)` on
`test_material_suite.py::test_axial_produces_bending_mitc3comp` was removed.

Measured, three independent symptoms resolved:

| symptom | before | after |
| --- | --- | --- |
| physical rigid-body mode, MITC3 `compute_ke_global` | penalised, `norm(Ku)/norm(K_bs) = 0.786` | free, worst residual 1.252576e-17 |
| mixed mesh, triangles on the last row only | 27.637% | **0.682%** (= all-quad) |
| mixed mesh, alternating triangle rows | 99.813% | **1.009%** (= all-triangle) |
| MITC3Comp B-coupling, [0/90] axial | +2.627691e-03 | **-2.627691e-03**, 0.45% from CLT |

Pure meshes unchanged (0.682% / 1.015%), which is the control that this is a
convention fix and not a stiffness change. `cargo test -p aeroelast-core` 89
passed / 0 failed; `pytest -m "not slow"` 345 passed / 4 skipped / 0 failed.

**Open, and it changes how this is documented, not what was done.** The two
reference papers write the rotation differently:

| Paper | Offset displacement |
| --- | --- |
| MITC3+ 2014, Eq. (2) | `(t/2) Σ a_i h_i (V_i2·α_i + V_i1·β_i)` - rotation parameters α, β |
| MITC4/D 2025, Eq. (3a) | `(t/2) Σ a_i h_i (θ_i × V_in)` - the rotation vector θ |

With `V_in = e3` and `V_i1 = e1`, Eq. (3a) gives `u_y = -z·θx` while Eq. (2) gives
`u_y = +z·α` if `V_i2 = V_in × V_i1` and `u_y = -z·α` if `V_i2 = V_i1 × V_in`. So
MITC3+'s α is either `+θx` or `-θx` depending on how the paper defines `V_i2`,
and that definition could not be extracted reliably from the recovered PDF (the
subscripts come out mangled). Either way the fix is right for a codebase that
mixes elements, but the document must say which of these it is:

- if MITC3+ defines the basis as MITC4/D does, then `b136ce5` **deviated from
  MITC3's own paper** and the fix restores fidelity;
- if MITC3+ uses the opposite sign, then the code was faithful to each paper and
  the **two papers use different DOF conventions**, so mixing them requires
  standardising on one, with the 2025 paper as the tiebreaker.

Also still to do, per directive 3: search for the most recent MITC3 work (a
MITC3/D equivalent of the 2025 paper would align the two families by
construction) and record the reference either way.

### U6 - MITC4/D drill operator, evaluated before it is wired

- Restore `b_md_mitc4_plus` and `drill_midside_shape_derivatives` from `929db32`
  (deleted in `4a46af2`), re-implement them faithfully to the 2025 paper, and
  validate.
- Wire it **only if** it significantly improves the element: measure with the
  benchmark suite (the paper's tables), the mixed-mesh cases and the patch/rigid
  tests. Record the measured before/after either way.
- Note the deleted version's test evaluated the operator at the element centre,
  where its midside derivatives are identically zero, so it proved nothing.

### U7 - Replace the ANSYS citations with the original source

- `src/aeroelast/solvers/fsi/rotor.py` ("cf. ANSYS MAPDL Theory Reference,
  Eq. 14-57, §14.4.1") and
  `crates/aeroelast-solvers/src/petsc/fsi/rotor_physics.rs:432` ("ANSYS Eq. 3-74 /
  14-55").
- Verify each equation against the original source before replacing; the likely
  sources are Geradin & Rixen (already in the README) and Bathe's FEP.

### U8 - Tests: one documented reference table and a validation matrix

- Single, documented reference table for the benchmark values, referenced to the
  paper's table and page (today `REFERENCE_VALUES` in
  `tests/test_ko2017_performance.py` is dead and every test hardcodes its own
  expectation, which is how the Hook's 1.12 got in).
- `docs/validation-matrix.md`: one row per test with what it validates, the exact
  reference (paper + table/equation), the tolerance and its justification, and the
  measured margin.

## Verified against the recovered papers

### The MITC4/D drill operator was a faithful implementation (correction)

Checked line by line against Ko, Bathe & Zhang 2025, and the earlier claim that
it carried an inverted director was **wrong** - that inverted-director comment
belongs to MITC3's transverse shear (`mitc3.rs:391-392`, commit `b136ce5`), a
different code path. What the paper actually says:

- Eq. (11a) defines the midside drill derivatives as
  `[h~5,r h~6,r h~7,r h~8,r] = [0, (1/2)(-2r)(1+s), 0, (1/2)(-2r)(1-s)]` and
  `[h~5,s h~6,s h~7,s h~8,s] = [(1/2)(-2s)(1+r), 0, (1/2)(-2s)(1-r), 0]`, i.e.
  `[0, -r(1+s), 0, -r(1-s)]` and `[-s(1+r), 0, -s(1-r), 0]`. The deleted
  `drill_midside_shape_derivatives` is literally those expressions.
- Eq. (5) defines `V_D = (x_r x x_s)/||x_r x x_s||` with `x_r = g_r(0,0,0)`,
  `x_s = g_s(0,0,0)` - the element-centre normal. The deleted code's
  `vd = e3/|e3|` is that vector.
- Eqs. (17a)/(17b) give the drill-membrane strains with the `(j0/j)` factor and
  the `x_m . (+/- x_r x V_D)` terms; the deleted code's `jac_ratio = sqrt_g0/sqrt_g`
  and its `coeff_rr/coeff_ss/coeff_rs` match them, including the factor 2 on the
  shear row for the engineering strain.

So the operator was faithful, and the user's recollection is right. What was
actually wrong with it: it was never wired into `compute_ke_local` (only its two
tests referenced it), and its test evaluated at `(xi, eta) = (0, 0)`, where
Eq. (11a) makes all four midside derivatives vanish identically, so the test
proved nothing. One comment also misattributed "Ko et al. 2025 Eq (21)" to the
centre-metric statement; Eq. (21) is the transverse-shear tying interpolation,
while the `j0/j` factor and the centre metric come from Eqs. (5) and (17).

### The rotation convention is settled by the paper (U5b's basis)

Ko, Bathe & Zhang 2025 Eq. (3a) interpolates the displacement as

```text
u(r,s,t) = sum_i h_i(r,s) u_i + (t/2) sum_i a_i h_i(r,s) (theta_i x V_in)
```

so the offset displacement is `(t/2)(theta x V_in)`: **theta is the physical
rotation vector**. Eq. (6a)-(6f) then split it into a standard part
`theta_S = alpha V_i1 + beta V_i2` and a drill part `theta_D = theta_D V_D`
measured about the element-centre normal, with
`delta gamma_i = delta theta_Di (V_in . V_D)`.

That is exactly MITC4's convention in this repository, and the opposite of what
MITC3's transverse shear implements after `b136ce5` (which makes `theta_x, theta_y`
the negation of the physical rotation while its curvature and drilling operators
keep the physical one). The most modern formulation in this element family
therefore confirms Finding 1 independently of the kinematic derivation: MITC3's
current shear sign is the odd one out, and the fix authorised as directive 3 is
justified by the source, not only by my own reading of the kinematics.

Still to extract for U2/U5: the paper's explicit covariant transverse-shear
definition in terms of `w` and the rotations (section 2.4 and the interpolation
Eq. (3)), to quote it directly in the code -> equation map.

## Verified code -> equation mappings

These are the entries of U2's map that have been checked term by term against the
recovered papers. Each is a verified identity, not a resemblance.

| Code | Paper | Status |
| --- | --- | --- |
| `mitc4.rs` `b_m_mitc4_plus` rows `b_rr`, `b_ss`, `b_rs` | Ko, Lee & Bathe 2017 ("A new MITC4+ shell element", CAS 182:404-418) Eqs. (27a), (27b), (27c) | **verified, exact** |
| `mitc4.rs` `drill_midside_shape_derivatives` (deleted in `4a46af2`) | Ko, Bathe & Zhang 2025 Eq. (11a) | **verified, exact** |
| `mitc4.rs` `b_md_mitc4_plus` (deleted in `4a46af2`) | Ko, Bathe & Zhang 2025 Eqs. (5), (17a), (17b) | **verified, exact** |
| Rotation convention of `compute_ke_global` | Ko, Bathe & Zhang 2025 Eq. (3a) | MITC4 matches; **MITC3 does not** (U5b) |
| SRI split in `compute_ke_local` (`mitc4.rs:946-953`) | Hughes, Taylor & Kanoknukulchai 1977 | cited, not yet re-verified |
| `mitc3.rs` `b_gamma_ext` constant part | Lee, Lee & Bathe 2014 **Eq. (15)** | **verified, exact** (Eq. (16) is the linear part; Eq. (17) the total) |
| `mitc3.rs` `b_gamma_ext` linear part | Lee, Lee & Bathe 2014 Eq. (16), tying points D/E/F | **verified, exact** |
| `mitc4.rs` `b_gamma_mitc4` | Dvorkin & Bathe 1984 (the original MITC4) | to verify |

### MITC4+ membrane operator, verified term by term

Paper Eqs. (27a-c), transcribed:

```text
(27a) e~rr = 1/2(1 - 2aA + s + 2aA s^2) e_rr(A)
           + 1/2(1 - 2aB - s + 2aB s^2) e_rr(B)
           + aC(-1 + s^2) e_ss(C) + aD(-1 + s^2) e_ss(D) + aE(-1 + s^2) e_rs(E)

(27b) e~ss = aA(-1 + r^2) e_rr(A) + aB(-1 + r^2) e_rr(B)
           + 1/2(1 - 2aC + r + 2aC r^2) e_ss(C)
           + 1/2(1 - 2aD - r + 2aD r^2) e_ss(D) + aE(-1 + r^2) e_rs(E)

(27c) e~rs = 1/4(r + 4aA rs) e_rr(A) + 1/4(-r + 4aB rs) e_rr(B)
           + 1/4(s + 4aC rs) e_ss(C) + 1/4(-s + 4aD rs) e_ss(D)
           + (1 + aE rs) e_rs(E)
```

The code's three rows are the same expressions with `a_a..a_e` and
`pre.b_rr_a`, `pre.b_rr_b`, `pre.b_ss_c`, `pre.b_ss_d`, `pre.b_rs_e` in the same
order. The tying points A-E and the blending coefficients are therefore in the
right places, which is worth having on record: this operator is the production
membrane path for MITC4, and it had no reference in the file beyond
"(Ko et al. 2017, Eqs. 27a-c)".

### U5c - The most recent MITC3 work (searched)

**There is no MITC3/D.** The 2025 `/D` paper covers the 4-node elements only
(its own title: "the MITC4/D and MITC4+/D elements"), so the two families cannot
be aligned by adopting a triangular `/D` element. The modern MITC3 line is
different work, and it targets a different weakness: the MITC3+ membrane field is
that of the displacement-based constant-strain triangle, and the improvements
below enrich the membrane behaviour.

| Option | Reference | What it changes | Cost |
| --- | --- | --- | --- |
| **Interpolation covers** | Jun, H., Yoon, K., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element enriched in membrane displacements by interpolation covers", *Comput. Methods Appl. Mech. Engrg.* 337:458-480, 2018, doi:10.1016/j.cma.2018.04.007 | Only the membrane displacement field, via cover-based enrichment `u = standard + Σ H_i û_i`. **Does not change the rotation convention and does not change the assumed transverse shear field.** Passes the isotropy, patch and zero-energy tests. | 9 DOFs per corner node instead of 5 (element = 27 nodal + 2 internal, the 2 condensed) |
| Strain-smoothed MITC3+ | *Computers and Structures*, 2019, doi:10.1016/j.compstruc.2019.07.005 | Improves the membrane strain with the strain-smoothed element method | No special smoothing domains, no additional DOFs |
| CC-MITC3+ | *Archive of Applied Mechanics*, 2025, doi:10.1007/s00419-025-02922-4 | Constant-curvature correction in a Hu-Washizu three-field framework, plus an Allman-like membrane, for rib-stiffened shells | not extracted |
| CS-MITC3+ / CS-MITC18+ | *Thin-Walled Structures*, 2024 (S0263823124006955); HSDT-type variant, 2026 | Cell-based smoothed, Allman-type with **real drilling DOFs** | not extracted |
| Benchmark suite | "Benchmark tests of MITC triangular shell elements", *Structural Engineering and Mechanics* 68 | Not a formulation: a benchmark reference for the MITC3 family, useful for the validation matrix | - |

**The reported gain from interpolation covers is large**, which matters because it
says the deficiency is real rather than cosmetic. From the paper's own tables,
normalized values against the reference solution:

- Cook's skew beam, mesh I: **0.95 / 0.99 / 1.00 / 1.00** versus MITC3+
  0.50 / 0.76 / 0.92 / 0.98; mesh II: **0.84 / 0.96 / 0.99 / 1.00** versus
  MITC3+ 0.28 / 0.47 / 0.72 / 0.90;
- MacNeal's cantilever, tip shear and tip moment: **~0.95-1.00** versus MITC3+
  **~0.011-0.037**, i.e. the standard MITC3+ is roughly a factor of 30 off on
  that problem;
- Scordelis-Lo roof: **0.9610 / 0.9931 / 0.9983** versus MITC3+
  0.7312 / 0.8743 / 0.9593.

**Recommendation.** The current MITC3+ is now correct (rigid-body invariant,
mixed meshes at the pure-mesh accuracy, B-coupling sign right) and it is what the
recovered 2014 paper specifies. Moving to an enriched MITC3 is a **new feature**,
not a defect fix, and the interpolation-cover option changes the number of DOFs
per node, so it is an architectural decision rather than a local edit. It is also
orthogonal to the convention work already done, since it leaves the rotation
convention and the transverse shear field untouched. That decision is the user's;
the references and the measured baseline are recorded here either way.

To make the decision concrete, one cheap measurement is still missing: whether
**this repository's** MITC3+ reproduces the reported membrane deficiency (the
MacNeal cantilever and Cook's skew beam numbers above), because if our element
already behaves better than the published MITC3+, the upgrade is less urgent.

**The 2019 strain-smoothed element is the best upgrade for this codebase**, and
the user supplied both papers locally (`.sources/papers/jun2018.pdf`,
`.sources/papers/lee2019.pdf`).

`lee2019.pdf` = Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite
element", *Computers and Structures* 223:106096, 2019. Its own abstract states the
key property: *"The major advantage of the SSE method is that no additional degree
of freedom is required for solution improvement."*

What it does, from the formulation section:

- computes the covariant membrane strain of the target element and of its three
  edge neighbours **at the element centres** (`r = s = 1/3`, `t = 0`);
- transforms each neighbour's strain into the target's convected coordinates
  (Eq. 15) and averages it with the target's, weighted by area (Eq. 16):
  `e_ij^(k) = (e_ij^(e) A^(e) + e~_ij^(k) Ā^(k)) / (A^(e) + Ā^(k))`, with
  `Ā^(k) = (n^(e) · n^(k)) A^(k)` (Eq. 17), so the smoothing fades to nothing as
  the angle between the two elements approaches 90 degrees;
- the smoothed membrane strain **replaces** the covariant membrane strain, used
  directly in the 3-point Gauss integration;
- and explicitly: *"We use the originally defined b1 eij and b2 eij for the
  covariant bending strains. For the covariant transverse shear strains, we adopt
  the assumed strains of the MITC3+ shell element"*. So **the bending field, the
  transverse shear field and the rotation convention are untouched** - it is
  orthogonal to the sign fix already done, exactly like the covers.

It passes the patch, isotropy and zero-energy-mode tests.

### Measured comparison, from the papers' own tables

Normalized displacement against the reference solution (1.0 = exact).
Scordelis-Lo roof, `t/L = 1/100`, reference `w_ref = 0.3024`, Table 6 of the 2019
paper - note our repository already runs this exact benchmark and uses
`3.0240e-1` as its reference in `tests/test_ko2017_performance.py`:

| element | DOFs/element | 4x4 | 8x8 | 16x16 |
| --- | --- | --- | --- | --- |
| MITC3+ (what we implement) | 15 | 0.7409 | 0.8793 | 0.9618 |
| Enriched MITC3+ (covers, 2018) | 27 | 0.9610 | 0.9931 | 0.9983 |
| **Strain-smoothed MITC3+ (2019)** | **15** | **1.1017** | **1.0323** | **1.0075** |
| MITC4+ | 20 | 1.0476 | 1.0053 | 0.9977 |

Cook's skew beam, mesh II (Table 3 of the 2019 paper):

| element | DOFs/element | 2x2 | 4x4 | 8x8 |
| --- | --- | --- | --- | --- |
| MITC3+ | 15 | 0.2815 | 0.4698 | 0.7236 |
| Enriched MITC3+ | 27 | 0.8393 | 0.9611 | 0.9916 |
| **Strain-smoothed MITC3+** | **15** | 0.5154 | 0.8873 | 0.9830 |

So the smoothed element gets most of the enrichment's benefit at **zero DOF cost**,
and on Scordelis-Lo it is even better than the enriched one. The covers are better
on Cook's coarse mesh II (0.8393 versus 0.5154) but cost +80% DOFs per element.

### Implementation shape for us

Not a local element-kernel edit, but bounded and with no system-size change: the
smoothed membrane B for each triangle is a weighted combination of the target's
centre-B and its three edge neighbours' transformed centre-B's, so it needs
neighbour lookup (the topology already knows element connectivity) and a
pre-pass before assembly. The assembled matrix stays symmetric because each
element's contribution remains a quadratic form `B^T C B`. The coupling term
`k_mb_ext` would use the smoothed membrane B as well. Table 7 of the paper gives
timings, which I did not extract.

### Recommendation

The current MITC3+ is now correct and matches its 2014 paper. Between the two
modern options, the **strain-smoothed MITC3+ (2019) is the one to implement** if
the user wants the upgrade: same DOFs, membrane-only, orthogonal to the
convention work, passes the basic tests, and it closes most of a deficiency that
the papers measure as a factor of two to thirty depending on the benchmark. The
interpolation covers (2018) buy a little more on some problems at a real
architectural cost. Either way this is a **new feature, not a defect fix**, so it
is the user's call.

To justify it with our own numbers rather than the papers', one measurement is
still missing and is now cheap: `_assemble_global` in
`tests/test_ko2017_performance.py` hardcodes element code 4, so the benchmark
suite can only run MITC4 even though the source paper publishes both the MITC3+
and MITC4+ columns. Deriving the code from the element's node count (3 -> MITC3,
4 -> MITC4) unlocks running the published MITC3+ values against our element, and
feeds the validation matrix of U8 at the same time.

### S1-S4 - Implement the strain-smoothed MITC3+ (AUTHORISED by the user)

Reference: Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite
element", *Computers and Structures* 223:106096, 2019 (`.sources/papers/lee2019.pdf`).
The formulation below was extracted from the paper with `pdftotext -layout` on
page 5, which recovers the superscripts that the default extraction mangles.

#### The formulation, faithfully

Per target triangle `e`, with mid-surface area `A^(e)` and unit centre normal
`n^(e)`, and its three edge neighbours `k`:

1. **Neighbour strain in the target's convected coordinates**, Eq. (15):
   `e_ij^(k) = e_ln^(k) (g_i^(e)·g^l^(k)) (g_j^(e)·g^n^(k))`, `i,j = 1,2`, using the
   covariant base vectors of the target and the contravariant base vectors of the
   neighbour (`g_i^(k)·g^j^(k) = delta_i^j`). Out-of-plane strains are neglected.
2. **Pairwise smoothing**, Eq. (16), weighted by area:
   `ê_ij^(k) = ( e_ij^(e) A^(e) + e~_ij^(k) Ā^(k) ) / ( A^(e) + Ā^(k) )`, with the
   neighbour's area **projected onto the target's mid-surface plane**, Eq. (17):
   `Ā^(k) = (n^(e)·n^(k)) A^(k)`. So the smoothing fades to nothing as the angle
   between the two elements approaches 90 degrees.
3. **Boundary rule**, stated in the text right after Eq. (17): *"we use
   m ê_ij = m e_ij if the kth edge of the target element is located along
   boundary"*. A boundary edge has no neighbour, so the pairwise strain falls
   back to the target's own strain.
4. **Assignment to the three Gauss points**, Eq. (18), cyclic pairing:
   `e^(A) = (ê^(3) + ê^(1))/2`, `e^(B) = (ê^(1) + ê^(2))/2`,
   `e^(C) = (ê^(2) + ê^(3))/2`.
5. Eq. (19) gives the equivalent explicit interpolation with `p = 1/6`,
   `q = 2/3`, but the paper states it *"is not utilized in actual computation of
   the stiffness matrix. We use the assigned strains in Eq. (18) directly in the
   3-point Gauss integration"*.
6. **Everything else stays MITC3+**: *"We use the originally defined b1 eij and
   b2 eij in Eqs. (11) and (12) for the covariant bending strains. For the
   covariant transverse shear strains, we adopt the assumed strains of the MITC3+
   shell element, in Eqs. (7) and (8)."* So the bending field, the transverse
   shear field and the rotation convention are untouched - the change is
   orthogonal to the sign fix already made.

#### Work units

- **S1 - the smoothing operator.** Edge-neighbour topology for triangles plus the
  smoothed covariant membrane B per element (Eqs. 15-18), as a pure function with
  unit tests: a constant-strain patch must return that same strain (averaging
  identical values), the boundary rule must fall back to the element's own strain,
  and a flat two-element case must reproduce the expected area-weighted average.
- **S2 - wire it into the stiffness.** The smoothed membrane B replaces the
  covariant membrane strain in the MITC3 stiffness, including the
  membrane-bending coupling term; bending and transverse shear untouched. The
  element kernel must gain an entry point that accepts the precomputed smoothed B,
  since the smoothing is not a per-element-local quantity.
- **S3 - validate.** Patch, isotropy and zero-energy-mode tests (the paper says it
  passes all three); the published columns versus our element (Scordelis-Lo,
  Cook's skew beam, hyperbolic paraboloid); mixed MITC3/MITC4 meshes; the rigid
  body invariant; and the full suite. Note the measured target: on the
  Scordelis-Lo roof the paper reports 1.1017/1.0323/1.0075 for the smoothed
  element against 0.7409/0.8793/0.9618 for MITC3+, and on the von Mises stress
  error at point B, mesh I 24.76/13.30/6.99 against 45.56/22.52/10.66.
- **S4 - document.** Add the element to `docs/formulations/shell-elements.md` with
  its code -> paper -> equation map, and add the 2018 and 2019 papers to
  `docs/references.md` (the 2019 one is not there yet).

#### Prerequisite for S3

`_assemble_global` in `tests/test_ko2017_performance.py` hardcodes element code
4, so the benchmark suite can only run MITC4 even though the source paper
publishes both the MITC3+ and the MITC4+ columns. Deriving the code from the
element's node count (3 -> MITC3, 4 -> MITC4) is required to compare against the
published MITC3+ column, and it feeds U8's validation matrix at the same time.

#### What must not change

The 2014 MITC3+ behaviour has to remain available and correct: the smoothed field
is an addition, and if the smoothing is not applied (or the topology is absent)
the element must fall back to the current MITC3+ membrane field. The existing
tests, including the corrected moment-sign expectations, must keep passing.

### U2 - Shell element formulations (i) - FIRST VERSION DELIVERED

`docs/formulations/shell-elements.md`, 645 lines: scope and conventions, MITC4+,
MITC3+, and known deviations plus references. Written against the recovered PDFs
page by page, with every equation read from the named page and an explicit
verification status per mapping. MITC4/D, the strain-smoothed MITC3+ and the
plane quads are declared pending in later revisions.

**It corrected three of my own briefed claims**, each verified afterwards:

1. The MITC3+ **constant** shear part is **Eq. (15)**, not Eq. (16); Eq. (16) is
the **linear** part and Eq. (17) is the total field. I had it wrong in this
document and in memory. Re-extracting page 4 of the 2014 paper with
`pdftotext -layout` (which preserves the superscripts the default extraction
mangles) confirms it, and the code's own comments - "Constant part (Eq. 15)",
"Linear part (Eq. 16)" - were right all along. `b_gamma_ext` matches Eq. (15) for
the constant part and Eq. (16) for the linear part **exactly**, including
`c_hat = (e_rt(F) - e_rt(D)) - (e_st(F) - e_st(E))` and the `1/3`, `(3s-1)`,
`(1-3r)` factors.
2. `b_drill` **does** carry an attribution: `/// Drilling B-vector (Hughes &
   Brezzi, 1x24)`. It is incomplete - no title, year or venue, and
   `docs/references.md` has no Hughes-Brezzi entry - but it is not "no citation".
3. The `1/8` scaling and the rotation of `b_gamma_mitc4` **cannot** be verified
   against Dvorkin & Bathe 1984: the recovered scan has no text layer over its
   displayed equations. The document reproduces instead the same field as it
   appears unnumbered in Ko et al. 2017 section 2 and flags the rest as open.

### New findings from that work

- **No hourglass control reaches the assembled stiffness.**
  `compute_hourglass_stiffness` and `compute_hourglass_forces` are referenced only
  from `#[cfg(test)]` tests and defined at `mitc4.rs:2787,2794`;
  `compute_ke_local` returns `k_m + k_mb_coup + k_mb + k_bs + k_drill` with no
  hourglass term, and nothing in `assembly/assembler.rs` calls them. The
  "S4R-style" scaffolding (`hg_factor`, `h_vec`, `h_orth`, `hg_stiffness_factor`)
  is precomputed and unused. Either it is dead code to delete or a missing
  stabilisation to wire - a decision for the user.
- The MITC3+ 2014 enrichment coefficient could not be confirmed from the
  extraction: the operator glyph between `h_i` and `f4` is lost, so `hi - f4` and
  `hi - f4/3` are indistinguishable in the text. The code uses `hi - f4/3`.
  Recorded as unresolved, not as a defect.
- The sign of `alpha`, `beta` in the MITC3+ 2014 Eq. (2) versus the rotation
  vector in the MITC4/D 2025 Eq. (3a) remains open; the document records what was
  tried.

## Open questions

- Which of the recovered papers actually correspond to the implemented code, and
  where the code deviates from them. That is the subject of U2's code -> paper ->
  equation map. The MITC3 and MITC4 shear/curvature operators are the first
  entries, since their conventions are now known to disagree.
- Whether the MITC3+ 2014 paper (or a newer MITC3 paper, per directive 3) defines
  the shear with the same convention as MITC4/D. This is what decides whether the
  fix is a one-line sign flip or a re-derivation.

## Out of scope

- Rewriting the existing Spanish formulation docs
  (`docs/teoria_formulacion_fsi_rotor.md`, `docs/s4r_composite_shell_formulation.md`);
  the user will correct those separately.
- Committing the papers: they stay in the gitignored `.sources/papers/`.

## Status reconciliation (2026-09-29, tree at `e879eba`)

This file was written as a working plan and stopped at the S1-S4 plan; the work then continued
past it, so the units below were still listed as open here long after they closed. Re-read
against the tree rather than against this file:

| unit | state in this file | actual state at `e879eba` | evidence |
| --- | --- | --- | --- |
| U1 canonical bibliography | open | **closed** | `docs/references.md` (290 lines); commits `8aa65fd`, `6b13278`, `0e5f0a4` |
| U2 shell element reference | "first version delivered" | **closed**, including §2.3 and §3 | `docs/formulations/shell-elements.md` §2.1-2.6, §3.1-3.6, §4.1-4.5; commits `80f50c3`, `b097549`, `4d296d6` |
| U3 constitutive/materials | open | **closed** | `docs/formulations/materials.md` (508 lines); commit `95cfd16` |
| U4 solver formulations | open | **closed** | `docs/formulations/solvers.md` (532 lines); commit `95cfd16` |
| U5 MITC3 sign fix | done (and modernisation open) | **closed** | commit `d6f37fb`; the measured before/after is §4.1 |
| U6 MITC4/D drill operator | open | **closed by the archived change** | the faithful MITC4+/D carries the 2025 drill strain (§2.5); the dead `b_md_mitc4_plus` was deleted, not wired |
| U7 ANSYS citations | open | **closed** | no `ANSYS` under `crates/`; `src/aeroelast/solvers/fsi/rotor.py` records what the old equation numbers were |
| U8 validation matrix | open | **closed, then refreshed** | `docs/validation-matrix.md` re-measured at `e879eba` (386 passed / 0 failed / 0 skipped) |
| S1-S3 strain-smoothed MITC3+ | plan | **implemented and Rust-unit-tested, not wired** | `elements/smoothing.rs`; `elements/mitc3.rs` union layout + `smoothed_membrane_b`; `shell-elements.md` §4.3 |
| S4 document the smoothed element | plan | **OPEN** | §4.3 says it "is not described equation by equation in this document yet" |

The SDD change `mitc4plusd-faithful` that owned most of this is archived with **0 unchecked
tasks** (`openspec/changes/archive/2026-09-29-mitc4plusd-faithful/tasks.md`, 59 checked).

### What is genuinely left (documentation completeness, not defects)

1. **S4 - §4.3 equation-by-equation.** The strain-smoothed MITC3+ is implemented and has Rust
   unit tests, but no equation-level section and no Python surface.
2. **Decide the smoothed element's delivery.** It is not reachable from PyO3 or the assembler:
   `grep -rn "union\|smoothed" crates/aeroelast-py/src crates/aeroelast-core/src/assembly`
   returns nothing. Today it is a Rust kernel, not a usable element, so either it gets wired
   and a Python test, or the gap is stated as deliberate.
3. **§4.4 `quad.rs` citation.** Stated as a gap; needs a decision (cite a textbook source or
   keep the gap explicit).
4. **§4.5 alpha/beta sign** in the MITC3+ 2014 paper. Unresolved from the held copy; recorded
   as open rather than guessed.
5. **Bibliography.** `references.md`'s "Hughes, Taylor & Kanoknukulchai 1977" still carries
   "DOI to verify. No held copy" in `shell-elements.md` §References item 2.
6. **Minor:** `shell-elements.md` §References' closing paragraph still parenthetically lists
   the Hughes-Brezzi attribution as a missing bibliography entry, but `references.md:104` now
   has it; the paragraph is stale.

None of 1-6 is a defect in the element or the suite. They are the honest residue of the
documentation work, and the previous iteration's self-description ("left half-done") reads as
this residue, not as unfinished core work.
