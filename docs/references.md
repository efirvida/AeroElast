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
that grounds it, a `Source of an implemented feature` line when the code
implements the formulation but does not cite it, or a `Planned reference` /
`not yet cited by the code` line when the entry is a held source for the
validation campaign that the code does not use yet. Entries that the code neither
cites, implements, nor plans to validate against were removed.

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

- **OpenFAST** (NREL), Apache-2.0 — reference wind-turbine aeroelastic solver
  (AeroDyn BEM aerodynamics + ElastoDyn structural dynamics), planned as the
  code-to-code reference for the IEA 15 MW aerodynamic and aeroelastic cases.
  Upstream: <https://github.com/OpenFAST/openfast>
  *Installed for the validation campaign: conda-forge `openfast=4.2.1` in a
  dedicated `openfast` environment (`~/miniconda3/envs/openfast/bin/openfast`),
  which also ships the standalone drivers (`aerodyn_driver`, `inflowwind_driver`,
  …). The version is pinned to match the official IEA 15 MW model, whose CI
  requires `openfast<5.0` (its `tests/environment.yml`) because the reference
  deck uses AeroDyn input version 15.03 (`Wake_Mod`, `UA_Mod`); OpenFAST 3.5.4
  expects the older 15.00 naming (`WakeMod`) and rejects it. The tests will skip
  when the binary is absent, following `conftest.ccx_bin_or_skip`. Not yet used
  by the code.*

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

## 5. Verification benchmarks, reference models and validation literature

This section mixes two kinds of entry. The **benchmarks and manuals** first
(NAFEMS, Abaqus) are not journal articles and carry no DOI; neither manual is
held, so their report numbers are repository assertions. The **reference wind
turbine models and validation studies** after them are journal articles,
conference papers and reports: DOIs apply where one is printed or verifiable, and
each entry states whether the PDF is held. **Numbers that will be used as
validation references** (mode tables, material tables, stress values, equations)
must be read from the rendered page with vision before they are quoted, never by
`pdftotext`; text extraction silently drops or reorders the mathematics of the
recovered scans.

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

**Reference wind turbine models and validation studies.** These are the held
sources for the IEA 15 MW blade/blade-structure and aeroelastic validation
campaign. They are not yet cited by the code unless an entry says so.

- Gaertner, E., Rinker, J., Sethuraman, L., Zahle, F., Anderson, B., Barter, G.,
  Abbas, N., Meng, F., Bortolotti, P., Skrzypinski, W., Scott, G., Feil, R.,
  Bredmose, H., Dykes, K., Shields, M., Allen, C., Viselli, A., *Definition of
  the IEA Wind 15-Megawatt Offshore Reference Wind Turbine*, IEA Wind TCP Task
  37 / NREL/TP-5000-75698, March 2020, doi:10.2172/1603478.
  *Verified against `.sources/papers/75698.pdf` (title, full author list and date
  read from the first page; the DOI was confirmed against the Crossref record,
  which carries the alternate title "IEA Wind TCP Task 37: Definition of the IEA
  15-Megawatt Offshore Reference Wind Turbine"). This is the NREL reference
  definition of the turbine; the repository ships its WindIO rendering as
  `tests/IEA-15-240-RWT.yaml`. The report states a blade mass about 65 t, which
  the test suite compares against.*
  *Cited by the code: `tests/test_blade_iea15mw_validation.py` (mass vs the
  report's ~65 t).*

- Escalera Mendoza, A.S., Mishra, I., Griffith, D.T., "An Open-Source NuMAD
  Model for the IEA 15 MW Blade with Baseline Structural Analysis", AIAA SciTech
  2023 Forum, AIAA 2023-2093, doi:10.2514/6.2023-2093.
  *Verified against `.sources/papers/AIAA2023Paper_EscaleraMendoza.pdf` (title,
  authors and DOI printed on the first page). The model itself is distributed as
  **UTD-IEA15MWBlade** v1.0.0 (Apache-2.0), archived at
  doi:10.5281/zenodo.7392283, upstream
  <https://github.com/UTDGriffithLab/UTD-IEA15MWBlade>; the held copy ships the
  NuMAD `.nmd`/`.xlsx`, the ElastoDyn blade file and the BModes mode output.*
  *Cited by the code: `tests/test_blade_iea15mw_validation.py:70-81` (68,077 kg
  mass; Table 3 first two modes; section V DLC 1.4 root moment and tip
  deflection).*

- Luo, Y., Gao, Y., "Modeling and FEM Analysis of the IEA15MW Wind Turbine
  Blade", in Wan, D., Ning, D., Tian, H. (eds.), *Advanced and Emerging Marine
  Engineering Technologies*, Lecture Notes in Civil Engineering, Springer
  (Singapore), pp. 933–941. Book DOI: 10.1007/978-981-95-1487-8.
  *Verified against `.sources/papers/978-981-95-1487-8.pdf` (chapter at PDF
  pages 947–955; authors, title and printed page range read there; book DOI and
  ISBNs read from the front matter). The chapter's own DOI was not printed in the
  held copy. **Its tables and equations are a lead, not evidence**: read them
  with vision before use. Planned structural validation source (stress
  distribution and the most loaded regions under different load directions).*
  *Not yet cited by the code.*

- Toledo de Almeida, I., Lapa, G.V.P., Gay Neto, A., Müller de Almeida, S.F.,
  "Design and extreme structural analysis of wind turbine blades: Beam and shell
  model comparison and discussion for a 10-MW reference turbine", *Engineering
  Structures*, 334:120155, 2025, doi:10.1016/j.engstruct.2025.120155.
  *Verified against `.sources/papers/1-s2.0-S0141029625005462-main.pdf` (title,
  authors, journal, volume, article number and DOI printed on the first page).*
  *Planned methodology reference for the beam-vs-shell comparison; not yet cited
  by the code. It is a 10 MW turbine, not the IEA 15 MW.*

- Höning, L., Herráez, I., Stoevesandt, B., Peinke, J., "Aeroelastic
  instabilities of the IEA 15 MW rotor during extreme yaw maneuvers", *Wind
  Energy Science*, 11:1531–1552, 2026, doi:10.5194/wes-11-1531-2026.
  *Verified against `.sources/papers/wes-11-1531-2026.pdf` (title, authors,
  journal, volume, pages, year and DOI printed on the first page).*
  *Planned aeroelastic-stability validation reference; not yet cited by the code.*

- Bernardi, C., Cherubini, S., Manganelli, F., Della Posta, G., Leonardi, S.,
  De Palma, P., "Large Eddy Simulation of the IEA 15-MW Wind Turbine Using a
  Two-Way Coupled Fluid-Structure Interaction Model", *Wind Energy Science*
  discussion preprint, `wes-2025-120`, doi:10.5194/wes-2025-120.
  *Verified against the held revised manuscript `.sources/papers/we_fsi_nrel15mw.pdf`
  (title, full author list and affiliations) and the author reply
  `.sources/papers/wes-2025-120-AR2.pdf` (paper number `wes-2025-120`). The DOI
  was confirmed against the Crossref record; the held copy is the discussion
  preprint, so the final journal DOI, volume and pages are not yet assigned. The
  study itself compares its LES FSI against **OpenFAST**, which is what makes it
  the high-fidelity reference for this campaign.*
  *Planned FSI validation reference; not yet cited by the code.*

- Ramos-García, N., Kontos, S., Pegalajar-Jurado, A., González Horcas, S.,
  Bredmose, H., "Investigation of the floating IEA Wind 15 MW RWT using vortex
  methods Part I: Flow regimes and wake recovery", *Wind Energy*, 2021,
  doi:10.1002/we.2682.
  *Verified against `.sources/papers/Wind Energy - 2021 - Ramos‐García -
  Investigation of the floating IEA Wind 15 MW RWT using vortex methods Part I
  Flow.pdf` (title, authors and DOI printed on the first page).*
  *Planned aerodynamic validation reference (free-vortex method); not yet cited by
  the code.*

## 6. Aerodynamics, actuator-line/vortex modelling and post-stall extrapolation

### 6.1 BEM theory and post-stall extrapolation

- Viterna, Corrigan (1981) — the post-stall polar extrapolation. **Attribution
  incomplete:** the code gives the author surnames, the year and the equation
  only — `# Viterna & Corrigan (1981) eq. 8: Cd_max = 1.11 + 0.018·AR` — with no
  title, venue, volume or pages, and no copy is held, so those are not asserted
  here. DOI: to verify.
  *Cited by the code: `src/aeroelast/models/blade/aerodynamics.py:209` (the brief
  named `src/aeroelast/constitutive/aerodynamics.py:209`; that path does not exist
  in this tree).*

- Ning, S.A., "A simple solution method for the blade element momentum equations
  with guaranteed convergence", *Wind Energy*, 17(9):1327–1345, 2014,
  doi:10.1002/we.1636.
  *Verified against the held preprint `.sources/papers/ning2013.pdf` (author
  S. Andrew Ning, NREL; title, journal and DOI `10.1002/we.1636` printed on the
  first page). This is the BEM solution method implemented by CCBlade, the solver
  the repository wraps, so it is the theory behind the aero side.*
  *Used by the code: `src/aeroelast/solvers/bem/engine.py` (through CCBlade).*

- Moriarty, P.J., Hansen, A.C., *AeroDyn Theory Manual*, NREL/TP-500-36881,
  January 2005, 47 pp. DOI: to verify (NREL technical report; no DOI printed).
  *Verified against the held copy `.sources/papers/36881.pdf` (title, authors,
  report number, month/year and page count read off the first page and PDF
  metadata). It is the theory source for the BEM aerodynamic model and its
  tip/root losses, and the planned reference for the BEM-vs-OpenFAST comparison.
  Its equations must be read with vision before use.*
  *Cited by the docs: `docs/teoria_formulacion_fsi_rotor.md:867`,
  `docs/FSI_ROTOR_PAPER_DRAFT.md:655`.*

- Shen, W.Z., Mikkelsen, R., Sørensen, J.N., "Tip loss corrections for wind
  turbine computations", *Wind Energy*, 8(4):457–475, 2005, doi:10.1002/we.153.
  *Verified against `.sources/papers/Wind Energy - 2005 - Shen - Tip loss
  corrections for wind turbine computations.pdf`; DOI, volume and pages confirmed
  against Crossref. Held. Reference for the BEM tip-loss correction in the
  BEM-vs-OpenFAST comparison.*
  *Not yet cited by the code.*

- Dağ, K.O., Sørensen, J.N., "A new tip correction for actuator line
  computations", *Wind Energy*, 23(2):148–160, 2020, doi:10.1002/we.2419.
  *Verified against `.sources/papers/Wind Energy - 2019 - Dağ - A new tip
  correction for actuator line computations.pdf`; DOI, volume and pages confirmed
  against Crossref. Held.*
  *Not yet cited by the code.*

### 6.2 Actuator-line, actuator-surface and vortex methods (CFD side)

These are the reference methods for the OpenFOAM track, not the BEM track; they
are catalogued here because aero validation has to name the method it compares
against. All held, and all DOIs cross-checked against Crossref except where a
venue has no DOI.

- Martínez-Tossas, L.A., Churchfield, M.J., Meneveau, C., "Optimal smoothing
  length scale for actuator line models of wind turbine blades based on Gaussian
  body force distribution", *Wind Energy*, 20(6):1083–1096, 2017,
  doi:10.1002/we.2081.
  *Verified against `.sources/papers/Wind Energy - 2017 - Martínez‐Tossas - Optimal
  smoothing length scale ... pdf`. Held.*

- Jost, E., Klein, L., Leipprand, H., Lutz, T., Krämer, E., "Extracting the
  angle of attack on rotor blades from CFD simulations", *Wind Energy*,
  21(10):807–822, 2018, doi:10.1002/we.2196.
  *Verified against the held copy; DOI and pages confirmed against Crossref. It is
  the reference for extracting the section angle of attack from CFD, which is what
  the A3 comparison needs.*

- Zormpa, M., Zilic de Arcos, F., Chen, X., "The Effect of Flow Sampling on the
  Robustness of the Actuator Line Method", *Wind Energy*, 28:e2965, 2025,
  doi:10.1002/we.2965.
  *Verified against `.sources/papers/Wind Energy - 2024 - Zormpa - ... pdf`; DOI
  confirmed against Crossref. Held.*

- Churchfield, M., Schreck, S., Martínez-Tossas, L.A., Meneveau, C., "An
  Advanced Actuator Line Method for Wind Energy Applications and Beyond:
  Preprint", AIAA SciTech 2017, NREL/CP-5000-67611. DOI: to verify (conference
  preprint).
  *Verified against the held copy `.sources/papers/67611.pdf` (title, authors,
  report number and venue read off the first page).*

- Jha, P.K., Churchfield, M.J., Moriarty, P.J., Schmitz, S., "Guidelines for
  Volume Force Distributions Within Actuator Line Modeling of Wind Turbines on
  Large-Eddy Simulation-Type Grids", *Journal of Solar Energy Engineering*,
  136(3):031003, 2014, doi:10.1115/1.4026252.
  *Verified against the held copy `.sources/papers/sol_136_03_031003.pdf`; DOI,
  journal, volume and article number confirmed against Crossref. Held.*

- Meyer Forsting, A.R., Pirrung, G.R., Ramos-García, N., "Brief communication: A
  fast vortex-based smearing correction for the actuator line", *Wind Energy
  Science*, 5(1):349–353, 2020, doi:10.5194/wes-5-349-2020.
  *Verified against the held copy `.sources/papers/MeyerForsting2020_fast-vortex-smearing-correction.pdf`
  (title, authors, journal, volume, pages and DOI printed on the first page).*

- Mohammadi, M.M., Olivares-Espinosa, H., Navarro Diaz, G.P., Ivanell, S., "An
  actuator sector model for wind power applications: a parametric study", *Wind
  Energy Science*, 9(6):1305–1321, 2024, doi:10.5194/wes-9-1305-2024.
  *Verified against the held copy `.sources/papers/Mohammadi2024_actuator-sector-model-parametric.pdf`
  (title, authors, journal, volume, pages and DOI printed on the first page).*

- Kim, T., Oh, S., Yee, K., "Improved actuator surface method for wind turbine
  application", *Renewable Energy*, 76:16–26, 2015,
  doi:10.1016/j.renene.2014.11.002.
  *Verified against the held copy `.sources/papers/1-s2.0-S0960148114007162-main.pdf`;
  DOI, journal, volume and pages confirmed against Crossref.*

- Shen, W.Z., Sørensen, J.N., Zhang, J., "Actuator surface model for wind
  turbine flow computations", *Proceedings of the European Wind Energy
  Conference (EWEC) 2007*. DOI: to verify (conference proceedings).
  *Verified against the held copy `.sources/papers/proceeding_shen.pdf` (DTU
  Orbit record: title, authors, venue, year).*

- Yang, X., Sotiropoulos, F., "A new class of actuator surface models for wind
  turbines", *Wind Energy*, 21(5):285–302, 2018, doi:10.1002/we.2162.
  *Verified against the held preprint `.sources/papers/1702.02108v4.pdf`
  (arXiv:1702.02108); DOI, journal, volume and pages confirmed against Crossref.*

- Santoni, C., Sotiropoulos, F., Khosronejad, A., "A Comparative Analysis of
  Actuator-Based Turbine Structure Parametrizations for High-Fidelity Modeling of
  Utility-Scale Wind Turbines under Neutral Atmospheric Conditions", *Energies*,
  17(3):753, 2024, doi:10.3390/en17030753.
  *Verified against the held copy `.sources/papers/energies-17-00753.pdf`; DOI,
  journal, volume and article number confirmed against Crossref.*

- Martínez-Tossas, L.A., Sakievich, P., Churchfield, M.J., Meneveau, C.,
  "Generalized filtered lifting line theory for arbitrary chord lengths and
  application to wind turbine blades", *Wind Energy*, 27(1):101–106, 2024,
  doi:10.1002/we.2872.
  *Verified against `.sources/papers/Wind Energy - 2023 - Martínez‐Tossas -
  Generalized filtered lifting line theory ... pdf`; DOI, volume and pages
  confirmed against Crossref. Held.*

- Troldborg, N., *Actuator Line Modeling of Wind Turbine Wakes*, PhD thesis,
  Technical University of Denmark (DTU Wind Energy / Risø), 2009. DOI: does not
  apply.
  *Verified against the held copy `.sources/papers/niels_troldborg.pdf` (DTU
  Orbit record: author, title, year).*

- Frontera Pericàs, P., *An FSI approach based on the actuator line and
  relaxation zone methods*, master thesis, 126 pp. DOI: does not apply.
  *Verified against the held copy `.sources/papers/Memòria.pdf` (title and author
  from the report cover, 126 pages). The awarding institution is not printed on
  the held first page; to verify.*

- *Aerodynamic Performance Simulation of the NREL Phase VI Wind Turbine*, held
  presentation (`.sources/papers/Aerodynamic Performance Simulation of the NREL
  Phase VI Wind Turbine.pptx`). No author, year or venue is printed on the first
  slide, so this is a **lead for the Phase VI reference, not a citation**. The
  Phase VI experimental data itself (Hand et al., NREL/TP-500-29955) would be the
  proper reference and is **not** held.

### 6.3 Rotorcraft and immersed-boundary CFD coupling (OpenFOAM side)

These are held references for the CFD coupling machinery (immersed boundary,
actuator surface) that the CFD-side workflow uses; they are rotorcraft or
wind-farm rather than blade-resolved wind turbine references.

- Linton, D., Barakos, G., Widjaja, R., Thornber, B., "A New Actuator Surface
  Model with Improved Wake Model for CFD Simulations of Rotorcraft",
  *Proceedings of the Vertical Flight Society 73rd Annual Forum and Technology
  Display*, pp. 1–10, 2017, doi:10.4050/f-0073-2017-12010.
  *Verified against the held copy `.sources/papers/138994.pdf`; DOI and venue
  confirmed against Crossref.*

- Linton, D., Barakos, G., Widjaja, R., Thornber, B., "Coupling of an Unsteady
  Aerodynamics Model with a Computational Fluid Dynamics Solver", *AIAA Journal*,
  56(8):3153–3166, 2018, doi:10.2514/1.J056784.
  *Verified against the held accepted manuscript `.sources/papers/161926.pdf`
  (title, authors, journal, volume, pages and DOI printed on the cover sheet).*

- Park, H.S., Linton, D., Thornber, B., "Towards DES of Flow Around a Rotorcraft
  Fuselage Using an Immersed Boundary Method", AIAA SciTech 2020. DOI: to verify.
  *Verified against the held manuscript
  `.sources/papers/SciTech2020_Manuscript_HSPark_DLinton_BThornberf.pdf`.*

- Park, J., Linton, D., Thornber, B., "Wind Farm Detached-Eddy Simulations Using
  an Immersed Boundary Method Actuator Surface Model Solver", presentation, DLR
  Institute of Software Methods for Product Virtualization, 04.06.2024. DOI: does
  not apply.
  *Verified against the held slides `.sources/papers/JackPark_04Jun24_elib.pdf`.*

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

### To obtain for the validation campaign

The aerodynamic campaign is the active one. **Free and public** unless marked
otherwise. The reachability note is from this sandbox: GitHub and Crossref
resolve, `nrel.gov` and `docs.openfast.org` do not.

| source | status | needed for |
| --- | --- | --- |
| OpenFAST source/release (`OpenFAST/openfast`) | to fetch and pin | the reference solver itself |
| OpenFAST `r-test` IEA 15 MW cases (`OpenFAST/r-test`) | to fetch (GitHub reachable) | the code-to-code case files (not a paper) |
| OpenFAST user documentation (`openfast.readthedocs.io`) | to fetch (reachable) | setting the OpenFAST cases and the pinned version |

**Now held, no action:** the **AeroDyn Theory Manual** (`.sources/papers/36881.pdf`)
and **Ning 2014** (`.sources/papers/ning2013.pdf`) were downloaded by the user and
are verified in §6 above.

### Deferred: structural and time-integration sources

Decision (user, 2026-09-30): postpone these books and paywalled structural
references until they are strictly necessary. None blocks the aerodynamic
campaign; they close the citation gaps listed in §1–§3.

| source | needed for |
| --- | --- |
| Newmark 1959, *J. Eng. Mech. Div. ASCE* 85(EM3):67–94 | the time-integration citation the code lacks |
| Géradin & Rixen, *Mechanical Vibrations*, 3rd ed. | replacing the unverifiable ANSYS gyroscopic/rotating-frame citation |
| Goldstein, Poole & Safko, *Classical Mechanics*, 3rd ed. | the non-inertial frame (centrifugal/Coriolis/Euler) treatment |
| Bathe, *Finite Element Procedures*, 2nd ed. | geometric stiffness; the UL question in `solvers.md` §4 |
| Hughes & Brezzi 1989, *CMAME* 72(1):105–121 | the drilling degree-of-freedom penalty |
| Hughes, Taylor & Kanoknukulchai 1977, *IJNME* 11(10):1529–1543 | the retired selective reduced integration (historical) |
| Reddy, *Mechanics of Laminated Composite Plates and Shells*, 2nd ed. | CLT and the B-coupling formula |
| Jones, *Mechanics of Composite Materials*, 2nd ed. | the constitutive/CLT citation |
| Tsai & Wu 1971; Hashin 1980 | the failure criteria |
| Viterna & Corrigan 1981 | the post-stall polar extrapolation |
| Hinton & Campbell 1974; Zienkiewicz & Zhu 1992 | the stress-recovery smoothing |

Priority: the aero track (OpenFAST source + docs + `r-test`) is next, with the
AeroDyn manual and Ning 2014 already held. The deferred structural books do not
block the IEA 15 MW aerodynamic campaign.
