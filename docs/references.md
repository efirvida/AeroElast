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
  Methods in Engineering*, 11(10):1529–1543, 1977. DOI: to verify.
  *Source: repository citation (`crates/aeroelast-core/src/elements/mitc4.rs:953`) for authors and year; the journal, volume and pages come from the published record and are not re-verified against a held copy. It is the reference for the selective reduced integration split of the in-plane shear term in `compute_ke_local`.*
  *Cited by the code: `crates/aeroelast-core/src/elements/mitc4.rs:953`.*

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
