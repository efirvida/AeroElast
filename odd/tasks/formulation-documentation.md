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

### U5 - MITC3 sign fix and modernisation

- Establish from the MITC4/D paper (and the MITC3+ papers) which rotation
  convention is authoritative, then fix `eval_covariant_shear_ext` and the five
  tests whose expectations encode the inverted convention
  (`test_linear_tip_moment_sign`, `test_cantilever_large_rotation_half_circle`,
  three `test_equilibrium_path` cases).
- Search for the most recent MITC3 papers and decide whether to move to a newer
  version; record the reference either way.
- Unblocks: the mixed-mesh correctness, the composite B-coupling sign, and the
  `xfail(strict=True)` marker added in `test_axial_produces_bending_mitc3comp`.

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
| `mitc3.rs` `b_gamma_ext` ("Eq. 15"/"Eq. 16") | Lee, Lee & Bathe 2014 (MITC3+, CAS 138:12-23) | **to verify**: the comment names the equations but not the paper |
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
