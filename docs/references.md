# Canonical bibliography

This file is the canonical bibliography for AeroElast. It replaces the two
reference lists that used to be duplicated in the README ("Scientific
references" and "Element formulation references"), and the README now points
here. Entries are grouped by subject and carry the full author list, title,
journal or publisher, volume, pages and year, plus a DOI **where one could be
verified**.

## About the recovered PDFs and the "to verify" markers

The reference PDFs used during development were committed, later deleted from
the working tree, and have been recovered into `.sources/papers/`, which is
gitignored (the `.gitignore` entry is `.sources`). They were recovered from git
history with
`git -c core.quotepath=false rev-list --all --objects | grep '\.pdf$'` and are
still reachable there (deletions in `10b860c` and `8baa4af`).

The README previously stated that "the reference PDFs used during development
were removed from the repository". That statement was misleading: the files were
removed from the working tree, but they remained recoverable from git history
and have now been recovered. The README no longer makes that claim.

Verification convention used below:

- `Verified against the recovered PDF` — authors, title, journal, volume, pages
  and year were read off the first page (and the printed footer or the PDF
  metadata where a DOI is given) of the file named in the entry.
- `Source: repository citation` — the entry was carried over from an existing
  citation in the code, README or `docs/`; it has **not** been re-verified
  against a held copy.
- `DOI: to verify` — no DOI could be read from the recovered material or from
  any source readable in this repository. Nothing was inferred or invented. This
  marker also covers items that may carry no DOI at all, such as conference
  proceedings and journals that do not register DOIs.

`.sources/papers/Implementación del MITC4_MITC4+ en MITC4.py vs formulaciones
teóricas.pdf` is a local working note (internal Spanish comparison of the
implementation against the formulations), **not** a publication, and therefore
has no entry in this bibliography.

Each entry below carries a `Cited by the code` line naming the code location
that grounds it, or a `Source of an implemented feature` line when the code
implements the formulation but does not cite it. Entries that the code neither
cites nor implements were removed.

## 1. Shell element formulations

- Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node shell
  element for general nonlinear analysis", *Engineering Computations*, 1(1):77–88,
  1984. DOI: to verify.
  *Verified against `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf` for authors, title, journal, volume, year, first page and length (12 pages, `© 1984 Pineridge Press Ltd`, printed code `0264-4401/84/010077-12`, footer "Eng. Comput., 1984, Vol. 1, March"). The issue number is not printed on the copy. This is the original MITC4 element.*
  *Cited by the code: `tests/test_ko2017_performance.py:1360`, `tests/test_shell_comprehensive.py:13`; it is also the source of the implemented MITC4 assumed transverse shear (`b_gamma_mitc4`) in `crates/aeroelast-core/src/elements/mitc4.rs`.*

- Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element and its
  performance", *Computers and Structures*, 138:12–23, 2014.
  DOI: 10.1016/j.compstruc.2014.02.005.
  *Verified against `.sources/papers/The_MITC3+_shell_element_and_its_performance.pdf` (authors: Youngyu Lee, Phill-Seung Lee, Klaus-Jürgen Bathe; DOI printed as `http://dx.doi.org/10.1016/j.compstruc.2014.02.005`). Cited by the MITC3+ kernel, `crates/aeroelast-core/src/elements/mitc3.rs:1-4`. It was missing from the README list.*
  *Cited by the code: `crates/aeroelast-core/src/elements/mitc3.rs:4`.*

- Jeon, H.-M., Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element in
  geometric nonlinear analysis", *Computers and Structures*, 146:91–104, 2015.
  DOI: 10.1016/j.compstruc.2014.09.004.
  *Verified against `.sources/papers/mitc3+_no_lineal.pdf` (authors: Hyeong-Min Jeon, Youngyu Lee, Phill-Seung Lee, Klaus-Jürgen Bathe; DOI printed as `http://dx.doi.org/10.1016/j.compstruc.2014.09.004`).*
  *Source of an implemented feature; the code does not cite it (see the audit).*

- Ko, Y., Lee, P.-S., Bathe, K.-J., "A new MITC4+ shell element",
  *Computers and Structures*, 182:404–418, 2017.
  DOI: 10.1016/j.compstruc.2016.11.004.
  *Verified against `.sources/papers/1-s2.0-S0045794916309464-main.pdf` and the duplicate `.sources/papers/A_new_MITC4+_shell_element.pdf` (authors: Yeongbin Ko, Phill-Seung Lee, Klaus-Jürgen Bathe; DOI printed as `http://dx.doi.org/10.1016/j.compstruc.2016.11.004`). The README entry omitted Bathe.*
  *Cited by the code: `crates/aeroelast-core/src/elements/mitc4.rs:4`, `crates/aeroelast-core/src/elements/mitc4.rs:954`.*

- Ko, Y., Lee, P.-S., Bathe, K.-J., "The MITC4+ shell element in geometric
  nonlinear analysis", *Computers and Structures*, 185:1–14, 2017.
  DOI: 10.1016/j.compstruc.2017.01.015.
  *Verified against `.sources/papers/mitc4+_no_lineal.pdf.pdf` and the duplicate `.sources/papers/The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf` (authors: Yeongbin Ko, Phill-Seung Lee, Klaus-Jürgen Bathe; DOI printed as `http://dx.doi.org/10.1016/j.compstruc.2017.01.015`). The README entry omitted Bathe.*
  *Cited by the code: `src/aeroelast/core/assembler.py:681`, `src/aeroelast/solvers/fsi/stress_stiffened_dynamic.py:52`.*

- Ko, Y., Lee, Y., Lee, P.-S., Bathe, K.-J., "Performance of the MITC3+ and
  MITC4+ shell elements in widely-used benchmark problems", *Computers and
  Structures*, 193:187–206, 2017. DOI: 10.1016/j.compstruc.2017.08.003.
  *Verified against `.sources/papers/1-s2.0-S0045794917309550-main.pdf` and the duplicate `.sources/papers/Performance_of_the_MITC3+_and_MITC4+_shell_elements_in_widely_used_benchmark_problems.pdf` (authors: Yeongbin Ko, Youngyu Lee, Phill-Seung Lee, Klaus-Jürgen Bathe; DOI printed as `http://dx.doi.org/10.1016/j.compstruc.2017.08.003`).*
  *Cited by the code: `tests/test_ko2017_performance.py:4-7`, `src/aeroelast/core/mesh/generators.py:65` (as "Ko2017 ratio-based mesh distortion").*

- Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite element",
  *Computers and Structures*, 223:106096, 2019.
  DOI: 10.1016/j.compstruc.2019.07.005.
  *Verified against `.sources/papers/lee2019.pdf` (authors: Chaemin Lee,
  Phill-Seung Lee — note the first author is **not** Youngyu Lee of the 2014 and
  2015 papers; DOI printed on the first page). It smooths the membrane strain
  with the target element and its three edge neighbours, at the same DOF count as
  MITC3+, and leaves the bending field, the transverse shear field and the
  rotation convention untouched. This is the formulation being implemented in
  this repository; see `odd/tasks/formulation-documentation.md` units S1-S4.*
  *Cited by the code: `crates/aeroelast-core/src/elements/smoothing.rs:3-4`.*

- Ko, Y., Bathe, K.-J., Zhang, X., "Continuum mechanics-based shell elements
  with six degrees of freedom at each node — the MITC4/D and MITC4+/D elements",
  *Computers and Structures*, 308:107622, 2025.
  DOI: 10.1016/j.compstruc.2024.107622.
  *Verified against `.sources/papers/1-s2.0-S0045794924003511-main.pdf` (authors: Yeongbin Ko, Klaus-Jürgen Bathe, Xinwei Zhang; DOI printed as `https://doi.org/10.1016/j.compstruc.2024.107622`). `.sources/papers/MITC_D_elements_published.pdf` is a scanned copy of the same paper with no text layer, so nothing can be read from it; it is a duplicate of the file above. The README entry carried a title only.*
  *Cited by the code: `crates/aeroelast-core/src/elements/mitc3.rs:393`.*

- Hughes, T.J.R., Brezzi, F., "On drilling degrees of freedom", *Computer
  Methods in Applied Mechanics and Engineering*, 72(1):105–121, 1989.
  DOI: to verify.
  *Identified from the reference list of the held copy of Ko, Bathe & Zhang 2025
  (`.sources/papers/1-s2.0-S0045794924003511-main.pdf`, item [13]), which cites it
  for the drilling degrees of freedom. The Hughes–Brezzi paper itself is **not**
  held, so this identification is second-hand and the DOI could not be read. The
  code names the pair with no work at all: `/// Drilling B-vector (Hughes &
  Brezzi, 1×24)` (`mitc4.rs:1085`) and `// Priority 3: Enhanced Drill Rotation
  (Hughes-Brezzi with variable penalty)` (`mitc4.rs:3335`). The `0.15` penalty
  factor used with that operator is attributed to nothing, in either element.*
  *Cited by the code: `crates/aeroelast-core/src/elements/mitc4.rs:1085`,
  `crates/aeroelast-core/src/elements/mitc4.rs:3335`.*

## 2. Constitutive and failure models

- Reddy, J.N., *Mechanics of Laminated Composite Plates and Shells: Theory and
  Analysis*, 2nd ed., CRC Press, 2004. DOI: to verify.
  *Source: repository citation — the classical lamination theory implementation cites this work at `crates/aeroelast-core/src/materials/laminate.rs:35` and `src/aeroelast/core/laminate.py:15`.*
  *Cited by the code: `crates/aeroelast-core/src/materials/laminate.rs:35`, `src/aeroelast/core/laminate.py:15`.*

- Jones, R.M., *Mechanics of Composite Materials*, 2nd ed., Taylor & Francis,
  1999. DOI: to verify.
  *Source: repository citation — cited alongside Reddy in the same CLT code paths (`crates/aeroelast-core/src/materials/laminate.rs:34`, `src/aeroelast/core/laminate.py:14`).*
  *Cited by the code: `crates/aeroelast-core/src/materials/laminate.rs:34`, `src/aeroelast/core/laminate.py:14`, `crates/aeroelast-core/src/materials/orthotropic.rs:69`.*

- Tsai, S.W., Wu, E.M., "A general theory of strength for anisotropic
  materials", *Journal of Composite Materials*, 5(1):58–80, 1971.
  DOI: to verify.
  *Source: repository citation — `crates/aeroelast-core/src/materials/failure.rs:10` gives authors, journal, volume and pages; the article title comes from the published record and is not re-verified against a held copy.*
  *Cited by the code: `crates/aeroelast-core/src/materials/failure.rs:9`, `src/aeroelast/constitutive/failure.py:12`.*

- Hashin, Z., "Failure criteria for unidirectional fiber composites",
  *Journal of Applied Mechanics*, 47(2):329–334, 1980. DOI: to verify.
  *Source: repository citation — `crates/aeroelast-core/src/materials/failure.rs:11` gives authors, journal, volume and pages; the article title comes from the published record and is not re-verified against a held copy.*
  *Cited by the code: `crates/aeroelast-core/src/materials/failure.rs:10`, `src/aeroelast/constitutive/failure.py:14`.*

## 3. Time integration and the finite element method

- Newmark, N.M., "A method of computation for structural dynamics", *Journal of
  the Engineering Mechanics Division, ASCE*, 85(EM3):67–94, 1959.
  DOI: to verify.
  *Source: published record; there is no PDF in `.sources/papers/` and the code cites nothing. The Newmark-β time integration is implemented in `src/aeroelast/solvers/elasticity/dynamic_newmark.py` and `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs`.*
  *Source of an implemented feature; the code does not cite it (see the audit).*

- Hughes, T.J.R., Taylor, R.L., Kanoknukulchai, W., "A simple and efficient
  finite element for plate bending", *International Journal for Numerical
  Methods in Engineering*, 11(10):1529–1543, 1977,
  doi:10.1002/nme.1620111005.
  *DOI verified against the Wiley record (`onlinelibrary.wiley.com/doi/10.1002/nme.1620111005`)
  and the MaRDI portal entry, both of which give volume 11, issue 10, pages
  1529–1543; no PDF is held. The previous entry said "DOI: to verify" and named
  `mitc4.rs:953` as the citing location, which is no longer true.*
  *Historical source of a RETIRED ingredient. It is the selective reduced
  integration of the in-plane shear term that the superseded hybrid MITC4+ carried.
  The faithful MITC4+/D that replaced it contains no selective reduced integration
  (`docs/formulations/shell-elements.md` §2.6, guarded by
  `tests/test_mitc4plusd_traceability.py`), and no line of the code cites this work
  any more, so it is kept as the origin of an ingredient whose retirement is
  documented rather than as a reference the code uses. The `Hughes, T.J.R., Brezzi,
  F.` entry below is a different work and is still the attributed origin of the
  retired drilling penalty.*

- Bathe, K.J., *Finite Element Procedures*, 2nd ed., Prentice Hall, 2014.
  DOI: to verify. *Source: repository citation (README, `docs/FSI_ROTOR_PAPER_DRAFT.md:647`, `docs/teoria_formulacion_fsi_rotor.md:862`).*
  *Cited by the code: `src/aeroelast/postprocess/stress_recovery.py:53` (2nd ed.); `src/aeroelast/solvers/fsi/stress_stiffened_dynamic.py:54` (1996 ed., §6.3 Updated Lagrangian).*

- Géradin, M., Rixen, D., *Mechanical Vibrations: Theory and Application to
  Structural Dynamics*, 3rd ed., Wiley, 2015. DOI: to verify.
  *Source: repository citation (README, `docs/FSI_ROTOR_PAPER_DRAFT.md:639`). Used for the gyroscopic matrices of rotating systems.*
  *Source of an implemented feature; the code does not cite it (see the audit).*

- Goldstein, H., Poole, C., Safko, J., *Classical Mechanics*, 3rd ed., Addison
  Wesley, 2002, §4.9–4.10. DOI: to verify.
  *Source: repository citation (README, `docs/FSI_ROTOR_PAPER_DRAFT.md:637`, `docs/teoria_formulacion_fsi_rotor.md:855`). A classical mechanics text rather than a finite element reference; it is listed here because it is the reference for the non-inertial (rotating) frame treatment used in the rotor FSI formulation, and it fits no other group.*
  *Source of an implemented feature; the code does not cite it (see the audit).*

- Cook, R.D., Malkus, D.S., Plesha, M.E., Witt, R.J., *Concepts and Applications
  of Finite Element Analysis*, 4th ed., 2002. DOI: to verify.
  *Source: repository citation — `src/aeroelast/postprocess/stress_recovery.py:53-55`
  gives the author list, the title and the edition; the publisher is not given and
  the work is not held. Also cited as "Cook et al. 'Concepts and Applications of
  Finite Element Analysis'" at `tests/test_shell_comprehensive.py:15`.*
  *Cited by the code: `src/aeroelast/postprocess/stress_recovery.py:53`,
  `tests/test_shell_comprehensive.py:15`.*

- Hinton, E., Campbell, J.S., "Local and Global Smoothing of Discontinuous Finite
  Element Functions Using a Least Squares Method", *International Journal for
  Numerical Methods in Engineering*, 8:461–480, 1974. DOI: to verify.
  *Source: repository citation — `src/aeroelast/postprocess/stress_recovery.py:57-59`
  gives authors, title, journal, volume and pages (no issue number); it is not
  re-verified against a held copy. It grounds the Gauss-to-node extrapolation in
  the stress recovery.*
  *Cited by the code: `src/aeroelast/postprocess/stress_recovery.py:57`.*

- Zienkiewicz, O.C., Zhu, J.Z., "The Superconvergent Patch Recovery and a
  posteriori error estimates", *International Journal for Numerical Methods in
  Engineering*, 33:1331–1364, 1992. DOI: to verify.
  *Source: repository citation — `src/aeroelast/postprocess/stress_recovery.py:60-63`
  gives authors, title, journal, volume and pages; it is not re-verified against a
  held copy. It is named as the theoretical basis for the nodal stress smoothing.*
  *Cited by the code: `src/aeroelast/postprocess/stress_recovery.py:60`.*

- Knight, N.F., "Raasch challenge for shell elements", *AIAA Journal*,
  35(2):375–381, 1997. DOI: to verify.
  *Read from the reference list of two held copies, `.sources/papers/1-s2.0-S0045794917309550-main.pdf`
  item [30] and `.sources/papers/mitc4+_no_lineal.pdf.pdf` item [36]; the paper
  itself is not held, so the DOI could not be read. The code cites it as
  "Knight (1997)" alone (`src/aeroelast/core/mesh/generators.py:2010,2023`,
  `tests/test_ko2017_performance.py:1455`), for the Raasch hook geometry and its
  benchmark parameters.*
  *Cited by the code: `src/aeroelast/core/mesh/generators.py:2010`,
  `tests/test_ko2017_performance.py:1455`.*

## 4. Software and vendored code

Software is identified by upstream URL and licence; DOIs do not apply.

- **pyNuMAD** (Sandia National Laboratories), BSD-3-Clause. A modified copy is
  vendored under `src/aeroelast/models/blade/numad/` for blade meshing; see the
  `LICENSE` and `NOTICE` files shipped with it. Upstream:
  <https://github.com/sandialabs/pyNuMAD>
  *Used by the code: `src/aeroelast/models/blade/aerodynamics.py:394`.*

- **preCICE**, LGPL-3.0 — partitioned FSI coupling.
  Upstream: <https://precice.org>
  *Used by the code: `crates/aeroelast-solvers/Cargo.toml:19`, `src/aeroelast/solvers/fsi/base.py:17`.*

- **OpenFOAM**, GPL-3.0 — CFD-side coupled simulations.
  Upstream: <https://openfoam.org>
  *Used by the code: `src/aeroelast/solvers/checkpoint.py:5` (restart and time-directory conventions).*

- **PETSc**, BSD-2-Clause — sparse linear algebra and KSP solvers.
  Upstream: <https://petsc.org>
  *Used by the code: `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs:1`, `src/aeroelast/solvers/elasticity/dynamic_newmark.py:9`.*

- **SLEPc**, BSD-2-Clause — eigenvalue problems for modal analysis.
  Upstream: <https://slepc.upv.es>
  *Used by the code: `src/aeroelast/solvers/elasticity/dynamic_newmark.py:11`.*

- **CCBlade**, Apache-2.0 — blade element momentum (BEM) aerodynamics.
  Upstream: <https://github.com/WISDEM/CCBlade>
  *Used by the code: `src/aeroelast/solvers/bem/engine.py:143`.*

- **NeuralFoil**, Apache-2.0 — neural-network airfoil aerodynamics.
  Upstream: <https://github.com/peterdsharples/NeuralFoil>
  *Used by the code: `src/aeroelast/models/blade/aerodynamics.py:316`.*

## 5. Verification benchmarks and manuals

These are not journal articles; DOIs do not apply. Neither manual or report is
held, so the section and report numbers below are repository assertions.

- **NAFEMS** benchmark reports — cited by report number only: "NAFEMS R7191,
  R7301 (linear benchmarks)" and "NAFEMS R0024 (nonlinear benchmarks)". NAFEMS
  is the organisation that publishes the benchmark reports.
  *Cited by the code: `tests/test_shell_comprehensive.py:11` ("NAFEMS benchmarks
  (linear and nonlinear)"), `tests/test_shell_comprehensive.py:82-83`. The code
  gives no report titles and none is held, so no title is asserted here.*

- **Abaqus** (Dassault Systèmes Simulia) — *Abaqus Theory Guide* 2016, §3.6.x
  (S4R composite shell), and the *ABAQUS Verification Manual*. Proprietary; DOIs
  do not apply.
  *Cited by the code: `tests/test_shell_comprehensive.py:12,84` ("ABAQUS
  verification manual"), `tests/test_shell_analytical_validation.py:9` ("Abaqus
  Theory Guide 2016, Section 3.6.x"), `crates/aeroelast-core/src/elements/mitc4.rs:95`
  ("S4R-style enhancements (Abaqus formulation)") and
  `docs/s4r_composite_shell_formulation.md:3` (Theory Guide 2016 §3.6.1, §3.6.5,
  §3.6.8). Neither manual is held.*

## 6. Aerodynamics and post-stall extrapolation

- Viterna, Corrigan (1981) — the post-stall polar extrapolation. **Attribution
  incomplete:** the code gives the author surnames, the year and the equation
  only — `# Viterna & Corrigan (1981) eq. 8: Cd_max = 1.11 + 0.018·AR` — with no
  title, venue, volume or pages, and no copy is held, so those are not asserted
  here. DOI: to verify.
  *Cited by the code: `src/aeroelast/models/blade/aerodynamics.py:209` (the brief
  named `src/aeroelast/constitutive/aerodynamics.py:209`; that path does not exist
  in this tree).*

## 7. Gap list

Entries here are **not** bibliography entries. They are sources the code cites
that must not be added, or citations that must be resolved before an entry can be
written.

- **ANSYS MAPDL Theory Reference** — cited at `src/aeroelast/solvers/fsi/rotor.py`
  ("cf. ANSYS MAPDL Theory Reference, Eq. 14-57, §14.4.1"; also lines 24, 77, 84,
  221, 234, 263, 287) and
  `crates/aeroelast-solvers/src/petsc/fsi/rotor_physics.rs:432` ("ANSYS Eq. 3-74 /
  14-55"), for the gyroscopic matrix, spin softening and the rotating-frame
  equation of motion. **Status: to be replaced.** The manual is proprietary and is
  not a verifiable source, and the user's decision is to replace it with the
  original source rather than add it here. The likely original for the gyroscopic
  and rotating-frame terms is Géradin & Rixen (already an entry in §3), with
  Goldstein, Poole & Safko for the non-inertial frame treatment. **The replacement
  must be verified against the original before it is made**, and none of the ANSYS
  equation numbers may be carried over unverified.
