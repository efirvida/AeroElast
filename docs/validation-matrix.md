# Validation matrix

One auditable row per test (or per coherent parametrised group) of the suite under
`tests/`, recording **what each test validates**, **the exact reference it validates
against**, **the tolerance as the code states it**, and **the margin the run actually
achieved**.

This document exists because the suite scatters reference values and tolerances
across 36 files: each test hardcodes its own expectation, so a wrong benchmark value
can survive inside an assertion that cannot fail. A single matrix makes every
reference and every margin visible.

## Guide: how an LLM must use this file to produce a validation paper

**Purpose.** This file is the single source of truth for AeroElast's numerical
validation evidence. It is written so that an LLM (or a careful human), without
reading the 36 test files, can (a) know exactly what each test checks, (b) know which
rows are strong evidence and which are flagged, and (c) regenerate every number a
validation paper will report. Treat it as the **input contract** for the paper's
experimental section.

**Which rows are usable as evidence.** Only rows whose *Notes* column carries no flag
from §9. In particular:

- A row flagged in §9.1 (tautology), §9.2 (assertion that cannot fail), §9.3 (name
  promises more than the body) or §9.4 (unjustified tolerance) may be cited as a
  **limitation or a negative finding, never as a validation result**.
- §9.5 lists values that are printed and never asserted: do not quote them as results.
- A reference that is the same formula re-implemented in the test is not a reference.
- §10.1-§10.2 are the strongest element-versus-published-columns evidence; §10.1 is not
  covered by a live test, so cite it explicitly as out-of-band evidence.

**Row schema.** Every row gives: the test (grouped by parametrisation), what it
validates, the independent reference (paper cell, `CCX <version>, <element>`, or an
analytical formula), the tolerance as the code states it, the measured margin, and
notes (flags per §9). `not printed` means the test asserts without printing a residual;
`not measured` means the residual was not obtained.

**Producing the data (recipe).**

1. `python -m pytest -o addopts="" -q -s` runs the whole suite and prints every
   residual; the *measured margin* column is exactly that output. Run it with `CCX_BIN`
   set, OpenFAST reachable and `neuralfoil` installed (§12), or the relevant rows skip.
2. `python -m pytest -o addopts="" --collect-only -q` reproduces the file-by-file
   inventory in §2 (441 tests / 36 files) that proves nothing drifted.
3. §2.1 is the headline table; §3-§8 are the per-row evidence behind it; §12 has the
   per-tool commands and the skip conditions.
4. For a paper table, copy the row's *reference*, *tolerance* and *measured margin*
   verbatim; never re-derive a margin that a test does not print.

**Mapping matrix sections to a validation paper.**

| paper section | matrix | what it supports |
| --- | --- | --- |
| Element formulation verification | §3, §10.1-§10.2, §6.1, §6.5 | MITC4/MITC3+ against the published benchmark columns and the K/mass invariants |
| Code-to-code verification | §4, §6.4, §8.4 | AeroElast vs CalculiX (S4/S8/S8R), vs the Rust/PETSc paths, and vs OpenFAST AeroDyn |
| Analytical benchmarks | §5 | closed-form references (Euler-Bernoulli, Timoshenko, Kirchhoff, elastica) |
| Constitutive and composite verification | §6.2, §6.3, §6.8 | CLT/ABD, failure criteria, B-coupling |
| Rotor and FSI verification | §7 | rotating-frame inertia, gyroscopic/implicit terms, checkpointing |
| Aerodynamic and aerodynamic-load validation | §8.1-§8.4 | polars, BEM, force projection, BEM-vs-AeroDyn, Viterna |
| Reproducibility | §2, §12 | the counts, the commands and the skip gates |

**Reading paper values (non-negotiable).** Numbers taken from the reference PDFs in
`.sources/papers/` must be read from the rendered page with vision
(`pdftoppm -png`), never with `pdftotext`, because the recovered scans have no text
layer over their mathematics. Any matrix cell sourced from a PDF is a **lead, not
evidence**, until re-read visually; the cells re-read this way are named where they
appear (§3).

**Honesty contract.** A green suite does not unflag a row (§9). A row that cannot fail
is not a result. When a row's reference is a preprint or a simplified model (e.g. the
Luo & Gao blade, the Bernardi modes), say so in the paper; the caveat travels with the
number.

## Index

Six groups, 36 files, 441 tests. Every collected test node is accounted for below;
§2 carries the file-by-file inventory table that proves the sum.

| group | section | files | tests | passed | failed | measured margin range |
| --- | --- | --- | --- | --- | --- | --- |
| Ko, Lee, Lee & Bathe 2017 benchmarks | [§3](#3-teststest_ko2017_performancepy-ko-lee-lee--bathe-2017) | 1 | 31 | 31 | 0 | 0.01% – 2.73% |
| CCX parity | [§4](#4-ccx-parity-group) | 10 | 66 | 66 | 0 | 0.37% – 12.57% |
| Analytical | [§5](#5-analytical-group) | 5 | 40 | 40 | 0 | 0.10% – 2.42% |
| Element and assembly invariants | [§6](#6-element-and-assembly-invariants) | 10 | 151 | 151 | 0 | 0.45% – 1.68% |
| Rotor and FSI | [§7](#7-rotor-and-fsi-group) | 5 | 95 | 95 | 0 | algebraic / invariant (`not printed`) |
| BEM, aero and mesh | [§8](#8-bem-aero-and-mesh-group) | 5 | 58 | 58 | 0 | 0.33% – 40.00% |
| **Total** | | **36** | **441** | **441** | **0** | |

The *measured margin range* covers the rows that print a numeric residual; the per-row detail
is in §13. The suite is green: 441 passed, 0 failed.

By file and subsection:

- [Guide: how an LLM must use this file](#guide-how-an-llm-must-use-this-file-to-produce-a-validation-paper)
- [1. How this matrix was produced](#1-how-this-matrix-was-produced)
- [2. Suite snapshot](#2-suite-snapshot)
  - [2.1 Results at a glance](#21-results-at-a-glance)
- [3. `test_ko2017_performance.py` (31)](#3-teststest_ko2017_performancepy-ko-lee-lee--bathe-2017)
- [4. CCX parity group](#4-ccx-parity-group)
  - [4.1 `test_beam_shell_4cases_parity.py` (9)](#41-test_beam_shell_4cases_paritypy-9)
  - [4.2 `test_isotropic_shell_parity.py` (1)](#42-test_isotropic_shell_paritypy-1)
  - [4.3 `test_composite_beam_parity.py` (5)](#43-test_composite_beam_paritypy-5)
  - [4.4 `test_orthotropic_shell_parity.py` (3)](#44-test_orthotropic_shell_paritypy-3)
  - [4.5 `test_shell_convergence.py` (2)](#45-test_shell_convergencepy-2)
  - [4.6 `test_ccx_shell_element_types_parity.py` (4)](#46-test_ccx_shell_element_types_paritypy-4)
  - [4.7 `test_composite_layup_parity.py` (18)](#47-test_composite_layup_paritypy-18)
  - [4.8 `test_blade_iea15mw_validation.py` (18)](#48-test_blade_iea15mw_validationpy-18)
  - [4.9 `test_ccx_writer_ids.py` (3)](#49-test_ccx_writer_idspy-3)
  - [4.10 `test_shell_stress_ccx_parity.py` (3)](#410-test_shell_stress_ccx_paritypy-3)
- [5. Analytical group](#5-analytical-group)
  - [5.1 `test_shell_analytical_validation.py` (11)](#51-test_shell_analytical_validationpy-11)
  - [5.2 `test_shell_comprehensive.py` (7)](#52-test_shell_comprehensivepy-7)
  - [5.3 `test_shell_validation_fixed.py` (7)](#53-test_shell_validation_fixedpy-7)
  - [5.4 `test_large_rotation_benchmarks.py` (7) and `test_mitc3_benchmarks.py` (8)](#54-test_large_rotation_benchmarkspy-7-and-test_mitc3_benchmarkspy-8)
- [6. Element and assembly invariants](#6-element-and-assembly-invariants)
  - [6.1 `test_quad_elements.py` (21)](#61-test_quad_elementspy-21)
  - [6.2 `test_material_suite.py` (41)](#62-test_material_suitepy-41)
  - [6.3 `test_mass_matrix_validation.py` (15)](#63-test_mass_matrix_validationpy-15)
  - [6.4 `test_rust_assembler.py` (19)](#64-test_rust_assemblerpy-19)
  - [6.5 `test_rust_composite.py` (24)](#65-test_rust_compositepy-24)
  - [6.6 `test_rust_modal.py` (10)](#66-test_rust_modalpy-10)
  - [6.7 `test_stress_stiffened_solver.py` (13)](#67-test_stress_stiffened_solverpy-13)
  - [6.8 `test_composite_b_coupling.py` (4)](#68-test_composite_b_couplingpy-4)
  - [6.9 Documentation and contract guards](#69-documentation-and-contract-guards)
- [7. Rotor and FSI group](#7-rotor-and-fsi-group)
  - [7.1 `test_rotor_inertial.py` (34)](#71-test_rotor_inertialpy-34)
  - [7.2 `test_rotor_physical_consistency.py` (22)](#72-test_rotor_physical_consistencypy-22)
  - [7.3 `test_rotor_rust_parity.py` (37)](#73-test_rotor_rust_paritypy-37)
  - [7.4 `test_rotor_performance_report.py` (1) and `test_fsi_structural_report.py` (1)](#74-test_rotor_performance_reportpy-1-and-test_fsi_structural_reportpy-1)
- [8. BEM, aero and mesh group](#8-bem-aero-and-mesh-group)
  - [8.1 `test_bem_polars.py` (20)](#81-test_bem_polarspy-20)
  - [8.2 `test_bem_engine.py` (14) and `test_blade_mesh.py` (1)](#82-test_bem_enginepy-14-and-test_blade_meshpy-1)
  - [8.3 `test_force_projection.py` (10)](#83-test_force_projectionpy-10)
  - [8.4 `test_bem_openfast_parity.py` (13)](#84-test_bem_openfast_paritypy-13)
- [9. Flag summary](#9-flag-summary)
  - [9.1 Tautological references](#91-tautological-references-the-arithmetic-under-test-re-implemented-in-the-test)
  - [9.2 Assertions that cannot fail](#92-assertions-that-cannot-fail)
  - [9.3 Names that promise more than the body delivers](#93-names-that-promise-more-than-the-body-delivers)
  - [9.4 Tolerances with no stated justification](#94-tolerances-with-no-stated-justification-tolerancecomment-mismatches-dead-conditionals)
  - [9.5 Informational-only comparisons](#95-informational-only-comparisons-that-are-printed-and-never-asserted)
  - [9.6 Configuration drift](#96-configuration-drift)
- [10. Validation evidence that exists outside `tests/`](#10-validation-evidence-that-exists-outside-tests)
  - [10.1 Element-versus-published-columns comparison](#101-the-element-versus-published-columns-comparison)
  - [10.2 The published columns, for the record](#102-the-published-columns-for-that-benchmark-for-the-record)
  - [10.3 The target for the strain-smoothed MITC3+](#103-the-target-for-the-strain-smoothed-mitc3)
- [11. Cross-references](#11-cross-references)
- [12. Reproducing this matrix](#12-reproducing-this-matrix)
- [13. Consolidated results matrix](#13-consolidated-results-matrix)

## 1. How this matrix was produced

**Current tree.** The **full `-s` run at this tree is `441 passed, 0 failed, 0 skipped` in
999.82s (16:39)**, with CalculiX 2.23, OpenFAST 4.2.1 and `neuralfoil` present so no row
skipped; the collected suite is **441 tests** (`python -m pytest -o addopts=""
--collect-only -q`), matching the file-by-file inventory in §2. The rows of §3-§7 also carry
the margins measured at the `e879eba` refresh (417 tests); the sections added since — the
BEM-vs-OpenFAST parity A1/A2/A3 of §8.4, the blade's eight Bernardi modes in §4.8 and the
isotropic stress recovery of §4.10 — carry their own measured margins from this tree, and
their tests are named where they appear. §2.1 collects the headline numbers and §13 the
per-row matrix, so every result the suite produces is visible in one place.

- Working tree: `e879eba`, clean; branch `test/physical-correctness`. This refresh re-ran
  every numeric row rather than restating the first pass, which was captured at `b2c62ff`
  (hundreds of commits behind) and mixed numbers from the retired MITC4+ hybrid with numbers
  from the faithful MITC4+/D that replaced it. That mix is what made the first version read as
  a defect where there was none (section 3, the thin twisted-beam rows).
- Interpreter: `~/miniconda3/envs/aeroelast-dev/bin/python`, after
  `source activate aeroelast-dev`. The repository's own `.venv` is an empty virtualenv; the
  dependencies and the built `_aeroelast` extension live in the conda environment.
- CalculiX: **2.23** from `~/miniconda3/envs/aeroelast-dev/bin/ccx`, on `PATH`, so the
  CCX parity tests ran instead of skipping.
- Suite run (baseline production of the `e879eba` rows):
  `python -m pytest -q` -> `386 passed, 0 failed, 0 skipped` in 633.45s.
- Printed margins: `python -m pytest -o addopts="" -q -s` -> every *measured margin* below is
  the value the test itself printed in that run. Rows whose test does not print are marked
  `not printed`, and that is a statement about the test, not about the element.
- Per-file test counts in the section headings were re-collected at this tree with
  `python -m pytest -o addopts="" --collect-only -q`; they sum to the 441 total in §2 (the
  `e879eba` subset of that collection summed to 386, and the 417-refresh subset to 417).
- Paper values were read from the recovered PDFs in `.sources/papers/`, not from
  second-hand notes. The Ko et al. 2017 benchmark tables cited below are from
  `.sources/papers/1-s2.0-S0045794917309550-main.pdf`.
- **Correction to how those values were read, and what follows from it.** This matrix was
  first assembled with `pdftotext -layout`, and the project's rule is the opposite: equations
  are read from the rendered page with vision (`pdftoppm -png`), never by text extraction,
  because the recovered scans have no text layer over their mathematics and text extraction
  silently drops or reorders it. Numbers inside plain tables usually survive that, equations
  do not. Treat therefore any single cell of this matrix as a **lead, not as evidence**: re-read
  it from the PDF visually before relying on it. The cells this change re-read with vision are
  named where they appear (the twisted-beam and cylindrical-patch columns of section 3).

Column conventions:

- **Reference** — the independent quantity. `Analytical P L^3/(3 E I)` is a reference;
  `CCX 2.23, S4` is a reference (element and version named); `Ko, Lee, Lee & Bathe 2017,
  Table 10, N=16` is a reference; *the same formula re-implemented in the test* is **not**
  a reference and is flagged as a tautology.
- **Tolerance** — quoted as the code has it (relative unless stated).
- **Measured margin** — the printed value. `not printed` means the test asserts without
  printing its residual; `not measured` means it was not obtainable in this pass.
- **Notes** — flags per section 9.

Shorthand: S4 = CCX 4-node linear shell, S8 = CCX 8-node quadratic full-integration
shell, S8R = CCX 8-node quadratic reduced-integration shell (used for
`*SHELL SECTION, COMPOSITE`), MITC4+ / MITC3+ = the elements documented in
`docs/formulations/shell-elements.md` §2 and §3. `write_ccx_mesh(quadratic=False)` maps
quads to `S4` and `quadratic=True` maps them to `S8R`. The explicit
`write_ccx_mesh(shell_element_type=...)` selector overrides that mapping with `"S4"`,
`"S8"` or `"S8R"` (and rejects `"S4"`/`"S8"` for a composite section, which CalculiX
only accepts as S8R/S6; `None` keeps the `quadratic` behaviour)
(`src/aeroelast/core/mesh/io/writers.py`).

## 2. Suite snapshot

The collected suite is **441 tests / 36 files** (`python -m pytest -o addopts=""
--collect-only -q`), and the file-by-file inventory below sums to the same 441. The path to
that number, so a reader can tell real drift from a stale cell:

| step | tests | what it added |
| --- | --- | --- |
| `e879eba` baseline | 386 | the first full refresh |
| + `test_composite_layup_parity.py` (18), `test_blade_iea15mw_validation.py` (10), `test_ccx_writer_ids.py` (3) | 417 | composite-layup and blade CCX parity, and the writer id-scheme test |
| + `test_bem_openfast_parity.py` (13) | 430 | A1/A2/A3 BEM vs OpenFAST AeroDyn (§8.4) |
| + `test_blade_iea15mw_validation.py` +8 (Bernardi modes) | 438 | the eight Bernardi et al. blade modes (§4.8) |
| + `test_shell_stress_ccx_parity.py` (3) | **441** | outer-fibre stress recovery vs CalculiX `OUTPUT=3D` (§4.10) |

The **full `-s` run at this tree is `441 passed, 0 failed, 0 skipped` in 999.82s (16:39)**,
with CalculiX 2.23, OpenFAST 4.2.1 and `neuralfoil` present so no row skipped. Earlier, for
reference: `e879eba` was `386 passed` in 633.45s and the 417-refresh was `417 passed` in
1025.15s. The Rust side is green too: `cargo test --manifest-path crates/Cargo.toml -p
aeroelast-core` -> **155 passed, 0 failed, 0 ignored** (the Cargo workspace root is `crates/`,
not the repository root).

### 2.1 Results at a glance

Every headline result the suite measures, in one table. The sections below carry the exact
reference, the tolerance as the code states it, and the flag for each row.

| group | headline measured result | section |
| --- | --- | --- |
| Ko et al. 2017, 31 benchmark cases | worst **2.73%** (Scordelis-Lo regular), tolerance 5%; most < 1% | §3 |
| Beam/shell vs CCX S4 | analytical ≤ 0.87%, CCX ≤ 0.70%, modal ≤ 0.73% | §4.1 |
| Isotropic shell vs CCX S4 | 2.8% (the CCX side is a max over all FRD nodes) | §4.2 |
| Composite layup vs CCX S8R | axial ≤ 1.83%, bending ≤ 0.51%, modal ≤ 1.26% | §4.7 |
| IEA 15 MW blade | mass **+3.7%** vs Escalera; CCX modes ≤ 1.65%; Escalera modes ≤ 8.1%; Bernardi modes ≤ 12.3%; static 1.9% vs CCX, -7.7% vs article | §4.8 |
| Isotropic stress vs CCX `OUTPUT=3D` | outer fibre **1.58%**, mid-surface 0, TOP = BOT | §4.10 |
| Rotor/foundation CCX element families | family spread 1.11%, each ≤ 0.61% vs analytical | §4.6 |
| CCX writer id-scheme | byte-identical decks, 0.0 displacement difference | §4.9 |
| **BEM vs OpenFAST AeroDyn, identical polars** | thrust ≤ **0.45%**, torque ≤ **0.80%**, interior span Δα ≤ 1.36° | §8.4 |
| BEM vs AeroDyn, yaw / shear | ≤ 0.33% (both), azimuth-averaged | §8.4 |
| BEM with the repo's own polars | sensitivity ≤ 2.12% (a sensitivity, not a parity) | §8.4 |
| NeuralFoil + Viterna vs official post-stall | attached Cl ≤ 4%; post-stall Cl ≤ 40%, Cd ≤ 17% | §8.4 |
| Large-rotation elastica | 5% relative per component | §5.4 |
| Shell K/mass invariants | symmetry, PSD, rigid-body and exact mass coefficients to 1e-12 | §6.1-§6.3 |

Two skips that the first version of this matrix recorded as verified are **resolved**, and
the two tests they hid now run:

| file | then | now |
| --- | --- | --- |
| `test_bem_engine.py` | 0 collected, module skipped (`ccblade` missing) | **14 passed** — CCBlade 1.3.1 is installed in the environment |
| `test_blade_mesh.py` | 0 collected, empty parameter set | **1 passed** — the module scanned a directory that does not exist and never ran; it now scans its own directory and raises if the parameter set is empty |

Two files were added by the change that produced the previous refresh and are now described
here: `test_mitc4plusd_traceability.py` (3, the documentation-to-code gate) and
`test_laminate_invariant_guard.py` (1, the composite surface guard).

**Per-file drift, corrected rather than left implicit.** The first version's per-file counts
were captured at `b2c62ff` and several no longer held. Corrected in the section headings below:
`test_rust_assembler.py` 19 (was 18), `test_material_suite.py` 41 (was 44),
`test_rotor_rust_parity.py` 37 (was 49), `test_rotor_inertial.py` 34 (was 38),
`test_rotor_physical_consistency.py` 22 (was 20), `test_bem_polars.py` 20 (was 18),
`test_rust_composite.py` 24 (was 22), `test_stress_stiffened_solver.py` 13 (was 14).

Four files the first version did not cover at all are added: `test_large_rotation_benchmarks.py`
(7) and `test_mitc3_benchmarks.py` (8) in section 5.4, `test_mitc4plusd_traceability.py` (3) and
`test_laminate_invariant_guard.py` (1) in section 6.9.

Files added after the first refresh: `test_composite_layup_parity.py` (18, composite layups
and their modal frequencies against CCX S8R, §4.7), `test_blade_iea15mw_validation.py` (18,
the IEA 15 MW blade against CCX, the Escalera Mendoza 2023 article **and the eight Bernardi
et al. modes**, §4.8), `test_ccx_writer_ids.py` (3, the CCX writer's id-scheme invariance,
§4.9), `test_shell_stress_ccx_parity.py` (3, outer-fibre stress recovery vs CalculiX
`OUTPUT=3D`, §4.10) and `test_bem_openfast_parity.py` (13, the A1/A2/A3 BEM-vs-OpenFAST
campaign, §8.4).

**File-by-file inventory.** The collected set at this tree, file by file, and the section
that documents it. This table is the index's audit: its column sums to `441`, and it is the
one place to check whether a file has drifted out of the matrix. Reproduce with
`python -m pytest -o addopts="" --collect-only -q | grep -c '::'`.

| file | tests | section |
| --- | --- | --- |
| `test_ko2017_performance.py` | 31 | §3 |
| `test_beam_shell_4cases_parity.py` | 9 | §4.1 |
| `test_isotropic_shell_parity.py` | 1 | §4.2 |
| `test_composite_beam_parity.py` | 5 | §4.3 |
| `test_orthotropic_shell_parity.py` | 3 | §4.4 |
| `test_shell_convergence.py` | 2 | §4.5 |
| `test_ccx_shell_element_types_parity.py` | 4 | §4.6 |
| `test_composite_layup_parity.py` | 18 | §4.7 |
| `test_blade_iea15mw_validation.py` | 18 | §4.8 |
| `test_ccx_writer_ids.py` | 3 | §4.9 |
| `test_shell_analytical_validation.py` | 11 | §5.1 |
| `test_shell_comprehensive.py` | 7 | §5.2 |
| `test_shell_validation_fixed.py` | 7 | §5.3 |
| `test_large_rotation_benchmarks.py` | 7 | §5.4 |
| `test_mitc3_benchmarks.py` | 8 | §5.4 |
| `test_quad_elements.py` | 21 | §6.1 |
| `test_material_suite.py` | 41 | §6.2 |
| `test_mass_matrix_validation.py` | 15 | §6.3 |
| `test_rust_assembler.py` | 19 | §6.4 |
| `test_rust_composite.py` | 24 | §6.5 |
| `test_rust_modal.py` | 10 | §6.6 |
| `test_stress_stiffened_solver.py` | 13 | §6.7 |
| `test_composite_b_coupling.py` | 4 | §6.8 |
| `test_mitc4plusd_traceability.py` | 3 | §6.9 |
| `test_laminate_invariant_guard.py` | 1 | §6.9 |
| `test_rotor_inertial.py` | 34 | §7.1 |
| `test_rotor_physical_consistency.py` | 22 | §7.2 |
| `test_rotor_rust_parity.py` | 37 | §7.3 |
| `test_rotor_performance_report.py` | 1 | §7.4 |
| `test_fsi_structural_report.py` | 1 | §7.4 |
| `test_bem_polars.py` | 20 | §8.1 |
| `test_bem_engine.py` | 14 | §8.2 |
| `test_blade_mesh.py` | 1 | §8.2 |
| `test_force_projection.py` | 10 | §8.3 |
| `test_bem_openfast_parity.py` | 13 | §8.4 |
| `test_shell_stress_ccx_parity.py` | 3 | §4.10 |
| **36 files** | **441** | |

**Row-level inventory corrections made with this refresh.** Four headings carried a group
count that did not sum to the file's collected total; the rows below were the cause and are
now fixed in place:

| file | group | was | now | the tests that were uncounted |
| --- | --- | --- | --- | --- |
| `test_bem_engine.py` | `TestBEMSolverParked` | 4 named | 9 named | `test_result_is_bem_result`, `test_r_matches_blade`, `test_Np_shape`, `test_Tp_shape`, `test_alpha_shape` (shape/type smoke tests, §8.2) |
| `test_bem_polars.py` | `TestBladeAeroYAML` | 11 | 13 | `test_polar_alpha_covers_range` and the three property-shape tests were counted as one; the enumeration is now explicit in §8.1 |
| `test_material_suite.py` | `TestABDMatrices` | 8 | 9 | the row listed nine checks under an 8 label; `test_isotropic_ab_d_ratios` was the ninth (§6.2) |
| `test_rust_composite.py` | `TestCompositeSanity` | 4 named | 6 named | `test_ke_positive_semidefinite`, `test_me_positive_semidefinite` (§6.5) |

One section-3 label caveat, so the count and the rows can be reconciled: `test_3_1`
parametrises only `distorted × t/L` (6 nodes) but runs the clamped and the simply-supported
case inside each node, so the four `test_3_1[reg/dist × clamped/SS]` rows describe 12 logical
benchmarks on 6 collected nodes. Every other section-3 row is one collected node per row.

## 3. `tests/test_ko2017_performance.py` (Ko, Lee, Lee & Bathe 2017)

All 31 cases print `Norm vs Kirchhoff: <value> (expected: <cell>, error: <x>%)`. The
*measured margin* column below is the value printed at `e879eba` and reproduced digit for
digit at `34e2328` (the current refresh). Every hardcoded
`expected_normalized` is the paper's MITC4/MITC4+ N=16 cell for the named table; the cells
that were mis-sourced (the distorted pinched cylinder, the four twisted-beam cells, the two
regular hemisphere-cutout cells) were corrected by the change archived at `2fd847d`, and the
*notes* column records what each row used to compare against. **A stale row here is how this
document can manufacture a false defect: the first version recorded an 8.5% deviation on the
thin twisted beam that no longer exists (the element now measures 0.9972 / 0.9981 against
0.9978 / 0.9982).**

Common tolerance: `rtol=0.05` via `assert_relative_error` (3.1-3.4, 3.6-3.9) and `tol=0.01`
per case in 3.5.

| test | what it validates | reference (paper cell) | tolerance | measured margin at `e879eba` (= `34e2328`) | notes |
| --- | --- | --- | --- | --- | --- |
| `test_3_1_square_plate_tables_2_to_5[reg clamped, t/L=1/100..1/10000]` (3 cases) | clamped square plate centre deflection, N=16 quarter model, uniform pressure | **Table 2**, MITC4 / MITC4+ (identical columns) N=16: 0.9984 / 0.9980 / 0.9979; `w_ref` alpha = 1.267e-3 in `p L^4/D` | `rtol=0.05` | 0.9996 (0.12%), 0.9980 (0.00%), 0.9979 (0.00%) | cell match |
| `...[dist clamped]` (3 cases) | same, distorted mesh | **Table 3**, N=16: 1.002 / 1.001 / 1.001 | `rtol=0.05` | 0.9991 (0.29%), 0.9975 (0.35%), 0.9974 (0.35%) | cell match |
| `...[reg SS]` (3 cases) | simply supported square plate centre deflection | **Table 4**, N=16: 1.000 / 0.9998 / 0.9998; alpha = 4.062e-3 | `rtol=0.05` | 1.0026 (0.26%), 0.9998 (0.00%), 0.9998 (0.00%) | cell match |
| `...[dist SS]` (3 cases) | same, distorted | **Table 5**, N=16: 1.003 / 1.003 / 1.003 | `rtol=0.05` | 1.0074 (0.44%), 0.9998 (0.32%), 0.9996 (0.34%) | cell match |
| `test_3_2_circular_plate_tables_6_to_7[clamped, 3 t/L]` | clamped circular plate centre deflection, N=16 quarter disk | **Table 6**, MITC4 N=16: 1.001 / 0.9997 / 0.9997; `w_ref = alpha p R^4/D`, alpha = 1/64 | `rtol=0.05` | 0.9971 (0.39%), 0.9967 (0.30%), 0.9967 (0.30%) | cell match |
| `...[ss, 3 t/L]` | simply supported circular plate | **Table 7**, MITC4 N=16: 0.9991 / 0.9988 / 0.9988; alpha = (5+nu)/(64(1+nu)) | `rtol=0.05` | 0.9975 (0.16%), 0.9974 (0.14%), 0.9974 (0.14%) | **RESOLVED**: the SS case now carries its own Table 7 tuple (`expected_clamped` / `expected_ss` are separate parametrisations). It used to share the Table 6 row, which the first version flagged. |
| `test_3_3_pinched_cylinder_tables_8_to_9[reg]` | pinched cylinder load-point displacement, N=16, R=300, t=3 | **Table 8**, MITC4+ N=16: 0.9313; `w_ref = 1.8248e-5` | `rtol=0.05` | 0.9182 (1.41%) | cell match; one of the two widest live margins in the module (Scordelis-Lo regular is the other). |
| `test_3_3_pinched_cylinder_tables_8_to_9[dist]` | same, distorted mesh | **Table 9**, MITC4+ N=16: 0.9321 | `rtol=0.05` | 0.9317 (0.04%) | **RESOLVED**: the expectation used to be 0.9892, which is a Table 12 N=2 MITC4 cell, not a Table 9 value. The pre-fix element measured 0.9943 (6.7% off); the faithful element measures 0.04% off. |
| `test_3_4_scordelis_lo_tables_10_to_11[reg]` | Scordelis-Lo roof free-edge centre, N=16, R=25, L=50, th=40deg, t=0.25, rho=360 | **Table 10**, MITC4+ N=16: 0.9973; `w_ref = 3.0240e-1` | `rtol=0.05` | 0.9701 (2.73%) | cell match; **the widest live margin in the module**. |
| `test_3_4_scordelis_lo_tables_10_to_11[dist]` | same, distorted | **Table 11**, MITC4+ N=16: 0.9942 | `rtol=0.05` | 0.9809 (1.34%) | cell match |
| `test_3_5_twisted_beam_tables_12_to_13[0.02667, In-plane]` | MacNeal-Harder twisted beam, thick, in-plane, N=16x96 | **Table 12**, MITC4+ N=16: 0.9971; `w_ref = 5.4240e-3` | `tol=0.01` | 0.9981 (0.10%) | **RESOLVED**: the hardcoded 1.02 was not a paper cell and is now 0.9971. |
| `...[0.02667, Out-of-plane]` | same, out-of-plane | **Table 13**, MITC4+ N=16: 0.9973; `w_ref = 1.7540e-3` | `tol=0.01` | 0.9998 (0.25%) | **RESOLVED**: hardcoded 0.99 replaced by 0.9973. |
| `...[0.0002667, In-plane]` | thin twisted beam, in-plane | **Table 12**, MITC4+ N=16: 0.9978; `w_ref = 5.2560e-3` | `tol=0.01` | 0.9972 (0.06%) | **RESOLVED**: hardcoded 0.92 ("what the element achieves") replaced by the paper's 0.9978; the element now reaches 0.9972, so the 8.5% deviation the first version recorded is gone. |
| `...[0.0002667, Out-of-plane]` | thin twisted beam, out-of-plane | **Table 13**, MITC4+ N=16: 0.9982; `w_ref = 1.2940e-3` | `tol=0.01` | 0.9981 (0.01%) | **RESOLVED**, same story. |
| `test_3_6_hook_table_14_minimal_fix[0.9782]` | Raasch hook tip deflection, N=8x48, t=2, E=3.3e3 | **Table 14**, MITC4+ N=8: 0.9782; `w_ref = 4.82482` (MITC9, N=64) | `rel_err < 0.03`, justified and measured in the comment | 0.9814 (0.33%) | cell match; the first version's 0.9927 and the historical 1.12 / `assert norm > 0.1` are recorded in the test comment. |
| `test_3_7_hemisphere_cutout_tables_15_to_16[reg, 4/1000]` | hemisphere with cut-out, regular mesh, N=16, R=10, th0=18deg | **Table 15**, MITC4+ N=16: 1.003 | `rtol=0.05` | 0.9963 (0.67%) | **RESOLVED**: hardcoded 1.009 was the Table 15 MITC4+ *N=8* cell. |
| `...[reg, 4/10000]` | same, thin | **Table 15**, MITC4+ N=16: 0.9834 | `rtol=0.05` | 0.9816 (0.18%) | **RESOLVED**: hardcoded 0.9811 was the Table 15 *S4* N=16 cell. |
| `...[dist, 4/1000]` | same, distorted mesh | **Table 16**, MITC4+ N=16: 0.9958 | `rtol=0.05` | 1.0011 (0.53%) | cell match |
| `...[dist, 4/10000]` | same, distorted, thin | **Table 16**, MITC4+ N=16: 0.9736 | `rtol=0.05` | 0.9871 (1.39%) | cell match |
| `test_3_8_full_hemisphere_table_17[4/1000]` | full hemisphere, N=16, th0=2deg | **Table 17**, MITC4+ N=16: 0.9960; `w_ref = 9.24e-2` | `rtol=0.05` | 0.9912 (0.49%) | cell match; the mesh skips the pole (th0 = 2 deg) whereas the paper's mesh is not stated to do so |
| `...[4/10000]` | same, thin | **Table 17**, MITC4+ N=16: 0.9798 | `rtol=0.05` | 0.9772 (0.26%) | cell match |
| `test_3_9_hyperbolic_paraboloid_tables_18_to_19[reg, ...]` (2 cases) | hyperbolic paraboloid free-edge centre, 32x32 on the full saddle | **Table 18**, MITC4+ N=16: 0.9762 / 0.9777 | `rtol=0.05` | 0.9681 (0.83%), 0.9681 (0.99%) | cell match |
| `...[dist, ...]` (2 cases) | same, distorted | **Table 19**, MITC4+ N=16: 0.9904 / 0.9936 | `rtol=0.05` | 0.9955 (0.52%), 1.0138 (2.03%) | cell match |

Structural notes on this module:

- `_assemble_global` derives the element code from the node count
  (`3 if len(node_ids_elem) == 3 else 4`), which is what allows the same benchmark to be
  run with quads or triangles; every live call site currently passes
  `use_triangular = False`, so the triangular path is unreachable from the suite (section 4
  exercises it out of band).
- The module-level `REFERENCE_VALUES` table named in the audit no longer exists; only
  `PAPER_REFS` exists (line 56) and it *is* read, by `test_3_1` only.
- Every case builds a `_Case` with `expected_paper`, and `_run_case` prints
  `[x] Norm vs Paper 3D` for it. That line is informational and never asserted, and it is
  wrong by construction: it divides `|disp|` by `PAPER_REFS[...] * pressure` and compares
  the quotient against `expected_normalized` (a *normalized* value). Printed errors run
  from 523% to 6.2e8%. Nothing fails on it (section 9.5).

## 4. CCX parity group

### 4.1 `test_beam_shell_4cases_parity.py` (9)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_linear_static_with_analytical[tension_axial]` | free-face mean axial extension of a 4x20 MITC4 cantilever, point load 1000 N | analytical `F L/(E A)`, stated as a formula in the docstring | `TOL_ANALYTICAL = 0.02`, justified in a comment against measured errors | 0.612% | the test never skips: it needs no CalculiX |
| `...[compression_axial]` | same, reversed | analytical `-F L/(E A)` | 0.02 | 0.612% | – |
| `...[bending_fx]` | in-plane bending | analytical `F L^3/(3 E I_y)` | 0.02 | 0.873% | – |
| `...[transverse_fy]` | out-of-plane bending | analytical `F L^3/(3 E I_x)` | 0.02 | 0.749% | – |
| `test_linear_static_vs_ccx[4 cases]` | MITC4 vs CCX at the loaded node | **CCX 2.23, S4** (quadratic=False), same 4x20 mesh, point load at free-face centre | `rel_err <= 0.05` | 0.70% (tension), 0.70% (compression), 0.62% (bending_fx), 0.28% (transverse_fy) | the loaded node over-reads the axial cases vs analytical (AE 4.37%, CCX 5.10%), for both codes; that is why the analytical test uses the face mean. The first version's 0.01% axial figure was the hybrid element's. |
| `test_modal_first_five_modes` | first 5 modes from 12 requested, matched by Hungarian assignment | **CCX 2.23, S4** eigenfrequencies, same mesh | `tol=0.05` | max rel err 7.27e-03 (0.73%) | the matching step (not the physical modes) is what keeps this test robust; the printed table shows modes 1-12 on both sides |

### 4.2 `test_isotropic_shell_parity.py` (1)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_transverse_tip_displacement` | MITC4 vs CCX out-of-plane tip displacement, 2x10 mesh | **CCX 2.23, S4**; the module also defines an analytical Mindlin tip formula that the test never uses | `tol = 0.05` | AE 0.149886 m, CCX 0.145747 m, ratio 0.972, difference 2.8% | (a) **RESOLVED**: the comment used to say "Allow 10% tolerance" over the enforced 5%; it now records that the comment was what was wrong. (b) The CCX side is `max abs(V)` over *all* FRD nodes (`:326-356`) while the AeroElast side is the loaded centre node, so the comparison is not like-for-like. The unused `analytical_tip_displacement` helper is dead reference code. |

### 4.3 `test_composite_beam_parity.py` (5)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestCompositeMaterial::test_laminate_abd_matrices` | ABD dict keys, lengths 9/9/9, thickness passthrough | none (schema test) | exact keys/`==` | not applicable | no numerical reference |
| `TestCompositeMaterial::test_mesh_connectivity` | node sets exist and sit at y=0 / y=L | mesh construction invariant | `atol=1e-12` | not printed | – |
| `test_composite_axial_tension` | [0/90/45/-45]s 8-ply laminate axial, 4x10 mesh | **CCX 2.23, S8R** + `*SHELL SECTION, COMPOSITE` (quadratic=True) | `rel_error < 0.1` | AE 45.87 um, CCX 48.09 um, 4.62% | the 10% window has no stated justification. The AeroElast block (mesh, material, assembly) is duplicated verbatim in the body. |
| `test_composite_isotropic_equiv` | same mesh with an isotropic equivalent mapped through a single-ply laminate | **CCX 2.23, S8R**, single isotropic ply (so only the element formulation differs) | `rel_error < 0.1` | AE 33.45 um, CCX 34.73 um, 3.67% | same unjustified 10% window |
| `test_composite_bending` | laminate transverse bending, 100 N at the free centre | **CCX 2.23, S8R** | `rel_error < 0.1` | AE 1671.71 um, CCX 1700.86 um, 1.71% | same unjustified 10% window |

### 4.4 `test_orthotropic_shell_parity.py` (3)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_orthotropic_axial` | single 0-ply orthotropic laminate, transverse Fy, 4x10 mesh, `elem_type=44` | **CCX 2.23, S8R** + composite section | `rel_error < 0.05` | CCX 3 172 580.00 um, 0.54% | name says "axial" but the body and docstring apply a transverse load |
| `test_orthotropic_bending` | [0/90/90/0] in-plane lateral bending, load spread over the tip | analytical `F L^3/(3 EI_in)` with `EI_in = A11 B^3/12`, written out in the docstring | `rel_error < 0.05` | AE 1266.10 um, analytical 1221.54 um, 3.6% | the module docstring claims the lateral-bending tolerance was widened for a ~25-30% MITC4-vs-S8R gap, yet this test compares against the analytical value at 5% and passes at 3.6%; no CalculiX in this test |
| `test_multi_layer_iso_equivalence` | 4 laminar plies vs 1 layer of the same total thickness, MITC4 vs MITC4Composite | internal identity (A, B, D algebra), documented in the docstring | `rel_diff < 1e-4` | 9.181e-17 | a consistency identity rather than an external reference; the docstring states the algebra, so it is a legitimate duplicate-path check |

### 4.5 `test_shell_convergence.py` (2)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_in_plane_bending_convergence` | observed order of the MITC4 in-plane tip displacement over 4 meshes, by Richardson self-convergence, plus the extrapolated limit | analytical `P L^3/(3 E I)` with `I = t B^3/12`; the Timoshenko-vs-Euler-Bernoulli shear floor is derived in the comment (9.6 um on 1230 um ≈ 0.8%) | `MIN_ORDER = 1.5` and `EXTRAPOLATED_TOL = 0.02`, both justified in the module docstring | orders 1.7314 and 1.7561 (agree, delta = 0.0247); Richardson limit 1238.1367 um vs 1230.7692 um = 0.5986% | the strongest tolerance justification in the suite; the raw pairwise orders (2.865, 0.164, -0.703) are printed and explicitly not asserted |
| `test_composite_laminate_gap_mesh_study` | whether the AeroElast-vs-CCX laminate gap shrinks (mesh artifact) or plateaus (formulation/ABD) | **CCX 2.23, S8R** across 4 meshes; no gap value asserted | none: the only assertions are `np.all(np.isfinite(...))` and `... > 0` | gaps -4.1765 / -1.7137 / -1.3512 / -1.7290%; verdict printed: `PLATEAUS -> formulation / ABD` | an analysis script, not a test: it cannot fail for any formulation. The docstring says "No gap value is asserted yet", so this is deliberate, but it should be read as a measurement, not as coverage |

### 4.6 `test_ccx_shell_element_types_parity.py` (4)

In-plane bending strip (defect #4): `L x b x h = 1.0 x 0.1 x 0.001 m`, isotropic
`E = 2.1e11`, `nu = 0.3`, `rho = 7800`, regular **8x4** quad mesh, every node at `x = 0`
clamped in 6 DOF, total `600 N` in `+y` on the free edge `x = L`, measured `uy` at the
free-edge centre node `(L, b/2)`. Analytical reference `F L^3/(3 E I)`, `I = h b^3/12` ->
`1.142857E-02 m`. Decks are written through the new
`write_ccx_mesh(shell_element_type=...)` selector; the strip is statically determinate, so
`S4` uses a consistent 4-node edge traction and `S8`/`S8R` a consistent 3-node edge
traction (corner `1/6`, midside `2/3`), each with resultant `600 N`.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_ccx_element_types_agree_with_each_other` | S4 vs S8 vs S8R centre deflection, 8x4 | **CCX 2.23, S4 / S8 / S8R** (three independent runs) | `TOL_AGREEMENT = 0.02`, justified in the module docstring against the mesh study | S4 1.135840E-02, S8 1.145430E-02, S8R 1.148490E-02 m; spread (max-min)/min = 1.1137% | the same deck at 4x2 / 8x4 / 16x8 gives spreads 2.67% / 1.11% / 0.42%, so the family difference falls with refinement |
| `test_ccx_element_type_matches_analytical[S4]` | S4 centre `uy` vs beam theory, 8x4 | analytical `F L^3/(3 E I) = 1.142857E-02 m` | `TOL_ANALYTICAL = 0.02`, justified in the module docstring | 1.135840E-02 m, 0.6140% | consistent 4-node edge traction, resultant 600 N |
| `test_ccx_element_type_matches_analytical[S8]` | S8 centre `uy` vs beam theory, 8x4 | same analytical | 0.02 | 1.145430E-02 m, 0.2251% | quadratic mesh; S8 is CCX full integration and is not in `ELEMENTS_TO_CALCULIX` -- the type is threaded through `_build_quadratic_mesh_data` |
| `test_ccx_element_type_matches_analytical[S8R]` | S8R centre `uy` vs beam theory, 8x4 | same analytical | 0.02 | 1.148490E-02 m, 0.4929% | same load; S8R is the `*SHELL SECTION, COMPOSITE` element |

The three independent CalculiX versions bracket the analytical value and agree with each
other within 1.11%, so the CCX reference model for the strip is not the source of any
AeroElast-vs-CCX gap on this problem. The first version of this matrix closed the paragraph
with a "37.7% AeroElast gap on the new shell element"; that number was carried over from the
pre-flip element and has no measured provenance at `e879eba`, so it was removed rather than
repeated. The module skips cleanly when CalculiX is absent (`conftest.ccx_bin_or_skip`, path
overridable with `CCX_BIN`).

### 4.7 `test_composite_layup_parity.py` (18)

Five layups of the same clamped-free composite strip (CFRP, L=1.0 m, B=0.1 m, t=5 mm),
compared against **CCX 2.23 S8R** with `*SHELL SECTION, COMPOSITE` on the same mesh. Load is
the total resultant distributed over the free edge and the measured quantity is the edge mean,
so the comparison is a structural response rather than a single-node local effect. This is the
CCX reference the composite material suite was missing; the layups are the ones that had only
CLT-analytical or invariant checks.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_axial_extension_matches_ccx[5 layups]` | mean axial extension under a 1000 N edge resultant, 4x10 mesh | **CCX 2.23, S8R** composite section | 2.5% | 0.10% (`uni_0`) to 1.83% (`uni_90`) | the `uni_90` limit is the soft axial direction |
| `test_transverse_bending_matches_ccx[5 layups]` | mean out-of-plane deflection under a 100 N edge resultant | **CCX 2.23, S8R** | 1% | 0.11% (`asym_0_90`) to 0.51% (`uni_0`) | – |
| `test_asymmetric_b_coupling_matches_ccx` | the asymmetric `[0/90]` strip bends out of plane under axial load | **CCX 2.23, S8R** | 2% | 0.37% (`aero 1.6816e-2` vs `ccx 1.6754e-2`) | **the case section 6.8 used to flag as sign-and-floor only**; the test also asserts the CCX reference is macroscopic (`> 1e-3`), so it cannot pass on round-off |
| `test_symmetric_laminates_have_no_b_coupling[sym_0_90s, quasi_iso]` | B = 0 control: symmetric laminates stay flat | **CCX 2.23, S8R** plus an absolute bound | 1e-11 absolute | aero ~1e-19, ccx ~1e-13 | an absolute bound because both sides are round-off; a relative test here would divide by zero |
| `test_modal_frequencies_match_ccx[5 layups]` | first five matched eigenfrequencies of the clamped-free strip | **CCX 2.23, S8R** `*FREQUENCY` | 3% | worst 1.26% (`uni_0`, 8x20 mesh) | the modal case runs on 8x20, not 4x10: at 4x10 the highest matched `uni_0` mode is 6.72% off, refining to 1.26% at 8x20 and 0.53% at 16x40. Matching is Hungarian over 10 requested modes |

### 4.8 `test_blade_iea15mw_validation.py` (18)

The IEA 15 MW reference blade, meshed from `tests/IEA-15-240-RWT.yaml` with this repository's
own `Blade` model at `element_size = 1.0` and given the composite shell properties the model
derives. Four references: CCX S8R on the same mesh (tight), the published mass and modes
(Escalera Mendoza et al. 2023, Table 3), the **first eight modes of Bernardi et al. 2025**
(`wes-2025-120`, Table 2, a beam-based CSD used in an LES FSI solver), and the published
DLC 1.4 response. Runtime is ~4.3 min, dominated by the two CCX runs.

**Both solvers must be told the span direction** (`(0, 0, 1)` here). The blade's ply angles are
defined relative to the span, and the assembler uses it to compute a per-element angle offset;
without it every ply sits in the element-local frame, the blade comes out about three times too
soft, and the first modes are 0.188 / 0.440 Hz with no 1st edgewise. That was a setup error in
the first revision of this test, not a property of the shell model, and the corrected setup is
what the numbers below measure.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_blade_mass_matches_published_models` | total elemental mass of the meshed blade | **Escalera Mendoza et al. 2023 (AIAA 2023-2093)**: 68,077 kg for the UTD NuMAD conversion; **Gaertner et al. 2020 (NREL/TP-5000-75698)**: about 65 t | 10% vs the article, and above the report but within 20% | **70,623 kg** = +3.7% over the article, +8.7% over the report | the article's own conversion is +4.33% over the report, so being above it is expected and asserted as a sign, not parity |
| `test_blade_modal_frequencies_match_ccx[0..4]` | first five matched eigenfrequencies of the clamped-root blade | **CCX 2.23, S8R** modal on the identical mesh, properties and span direction | 10% | worst 1.65%; all five: 0.44%, 0.71%, 0.78%, 0.84%, 1.65% | pairing is Hungarian over 10 requested modes |
| `test_blade_first_modes_match_article[flapwise, edgewise]` | the first two computed frequencies | **Escalera Mendoza et al. 2023, Table 3**: 0.57 Hz (1st flapwise), 0.65 Hz (1st edgewise) | 15% | 0.526 Hz (-7.6%) and 0.702 Hz (+8.1%) | with the span direction supplied the ordering maps directly onto the article's: the shell model is 7.6% softer than the BModes beam on flapwise and 8.1% stiffer on edgewise, which is the expected direction for a shell that restrains cross-section warping. The 2.0 m mesh gave the same ordering (0.528 / 0.708) |
| `test_blade_modal_frequencies_match_bernardi[0..7]` | the first eight matched frequencies | **Bernardi et al., Wind Energy Science preprint `wes-2025-120`, Table 2**: 0.5369 / 0.7267 / 1.577 / 2.267 / 3.113 / 3.642 / 4.571 / 5.385 Hz (1F, 1E, 2F, 2E, 3F, 1T, 3E, 4F) | 15% | worst 12.3% (highest pair); lower modes 1.9% / 3.4% / 4.4% / 5.3% / 5.8% / 6.8% / 9.8% | an independent beam-based CSD reference from the FSI literature; pairing is Hungarian over 10 computed modes. The shell sits progressively above the beam as the modes go up, the expected direction for a shell that restrains cross-section warping |
| `test_blade_static_tip_deflection_matches_ccx` | static flapwise tip deflection under a uniform load scaled to the article's DLC 1.4 root moment | **CCX 2.23, S8R** static on the identical mesh, properties and span direction | 15% | aero 21.69 m vs ccx 22.10 m = 1.9% | – |
| `test_blade_static_deflection_matches_article_dlc` | the same tip deflection against the article's reported value | **Escalera Mendoza et al. 2023, section V**: DLC 1.4 max root moment 90.4 MNm, max out-of-plane tip deflection 23.49 m | 15% | 21.69 m = -7.7% | the load distribution is a **proxy** (DLC 1.4 is aero-elastic); with the root moment matched, the tip deflection is the comparable quantity. A uniform-cantilever beam estimate from the article's own 1st flapwise frequency is **not usable** here: 0.46 m against the shell's 7.5 m for a tip load, a factor of 14, because the blade tapers hard and the tip-load compliance is dominated by the soft outboard section |

### 4.9 `test_ccx_writer_ids.py` (3)

The CCX writer's id handling. The entity id counters are process-global, so a mesh built after
another one has ids that are neither 0-based nor contiguous; the writer used to label the
`.msh` nodes and element blocks by 1-based index but write the connectivity and the `.nam`
set members from the raw ids, producing a deck that referenced entities that did not exist.
The fix is in `src/aeroelast/core/mesh/io/writers.py` (commit `930d055`).

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_deck_is_id_offset_invariant` | the exported `.msh`/`.nam`/`.inp` are byte-identical for 0-based and offset ids | the 0-based deck | exact equality | identical | also pins the deterministic ELSET order, which comes from a set |
| `test_deck_labels_are_internally_consistent` | every connectivity label and every set member points at an entity that exists | the `.msh` node/element blocks | set containment | consistent | the element `EPLATE` must equal the full element block |
| `test_offset_ids_give_the_same_ccx_result` | the displacement does not depend on the id scheme | **CCX 2.23, S4** on both decks | 1e-9 relative | 0.0 | the end-to-end proof that the deck is right, not only self-consistent |

### 4.10 `test_shell_stress_ccx_parity.py` (3)

The stress-recovery path, validated on an **isotropic** cantilever plate
(`L x B x H = 1.0 x 0.1 x 0.01 m`, steel, 8x2 quads, out-of-plane tip load 100 N).
Two controlled probes established that CalculiX `*EL FILE, OUTPUT=2D` + `S`
reports the **mid-surface** stress (0 for this pure-bending case) and that
`OUTPUT=3D` exposes the outer fibre, which is what AeroElast `StressRecovery`
returns at `TOP`/`BOTTOM`. The composite
`*SHELL SECTION, COMPOSITE` blade does **not** respond to `OUTPUT=3D` (its `S`
stayed at 899.67 MPa in both modes), so it is out of scope here.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_bending_has_no_membrane_stress` | mid-surface von Mises under pure out-of-plane bending | beam theory (0) | `< 2%` of the outer fibre | 0.00 MPa | the membrane/bending split is correct |
| `test_outer_fibre_stress_is_symmetric` | TOP and BOTTOM carry equal magnitude | symmetric section | exact | 52.46 / 52.46 MPa | - |
| `test_outer_fibre_stress_matches_ccx_and_analytical` | outer-fibre von Mises | **CCX 2.23, S8R `OUTPUT=3D`** and analytical `M c / I` = 60 MPa | 15% / 20% | 52.46 vs 51.64 MPa (1.58%); vs analytical 12.57% | the analytical gap is the coarse 8x2 linear mesh |

## 5. Analytical group

### 5.1 `test_shell_analytical_validation.py` (11)

No test in this module prints its residual.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestCantileverBeam::test_tip_deflection[4x2, 8x4, 16x8]` | cantilever tip load | analytical `P L^3/(3 E I)`, `I = b h^3/12` | `rel_error < 0.05` | not printed | the measured quantity is `np.linalg.norm([u_x, u_y, u_z])` of a corner node, not the transverse component the formula predicts; the 5% has no stated justification |
| `TestSimplySupportedBeam::test_center_deflection` | "simply supported" centre point load | analytical `P L^3(1-nu^2)/(192 E I)` — a documented plate-bending modification of the beam formula | `rel_error < 0.05` | not printed | **the boundary condition is not simply supported**: the "roller" loop (`:275`) applies a Dirichlet condition for *every* DOF of the right edge, so both ends are clamped. The name and docstring promise a support the model does not have. 5% unjustified. |
| `TestBeamBending::test_bending_convergence[4x2, 8x4, 16x8]` | "convergence" of tip deflection | analytical `P L^3/(3 E I)` | `rel_error < 0.05` | not printed | the name promises a convergence study; the body is three independent single-mesh assertions with no cross-mesh comparison |
| `TestMembraneStretching::test_axial_extension` | membrane tension | analytical `P L/(E A)`, `A = b h` | `rel_error < 0.05` | not printed | 5% unjustified |
| `TestCompositeMatrices::test_laminate_solve[[0], [0,90]]` | "solve" of a laminate | none: only `total_thickness > 0`, ply count, and A/B/D shapes | exact shapes; `> 0` for thickness | not applicable | the name promises a solve; the body checks object shapes |
| `TestShearLocking::test_thin_plate_convergence[0.1]` | thin-plate tip deflection vs beam theory | analytical `P L^3/(3 E I)` | `tol = 0.05 if thickness_ratio < 0.01 else 0.05` — both branches identical | not printed | dead conditional: the branch can never change the tolerance. The assertion is one-sided (`ratio > 1-tol`), so an over-flexible element is undetected. The module also defines `plate_bending_stiffness_D`, `torsion_angle`, `beam_shear_angle` and `simply_supported_plate_deflection`, none of which is called. |

### 5.2 `test_shell_comprehensive.py` (7)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestLinearStaticCantilever::test_fx_in_plane` | axial membrane load, 8x4 mesh | analytical `P L/(E A)` | `error < 5.0` (percent) | FEM 2.820e-05, ana 2.857e-05 -> 1.3% | the comment used to say "Should match within 10%" (`:400`) while the assertion is 5%; the test now records that the comment was the wrong part, not the tolerance |
| `test_fy_in_plane` | in-plane bending, ny=8 | analytical `P L^3/(3 E I_z)`, `I_z = h b^3/12` | `error < 5.0` | FEM 1.134e-02 vs ana 1.143e-02 (raw values only; no error line printed) | 5% unjustified |
| `test_fz_out_of_plane` | out-of-plane flexure | analytical Timoshenko `P L^3/(3 E I) + P L/(k G A)`, `k=5/6` | `error < 5.0` | FEM 1.120e+02 vs ana(bending only) 1.143e+02 (raw values only) | the printed reference is the bending-only term; 5% unjustified |
| `test_in_plane_ratio_constraint` | uY/uX ratio | beam theory `4 (L/b)^2 = 400` | `0.5 ratio_ref <= ratio <= 1.2 ratio_ref` | 400.39 vs 400.00 (0.10%) | the window is asymmetric (-50% / +20%) and unjustified; `test_shell_validation_fixed.py` asserts the same ratio at ±2% |
| `TestNonlinearStaticCantilever::test_large_displacement_tip_load` | "large displacement" response | none: the body asserts the linear estimate exceeds L, then that the solve raises `RuntimeError` | `abs(dz_lin) > L` | linear estimate 1.121e+02 | the name promises a numerical large-displacement validation; the body validates the divergence *error path* only |
| `TestModalAnalysis::test_first_mode_frequency` | first cantilever mode | analytical `(1.875104^2 / 2pi) sqrt(E I/(rho A L^4))`, `I = b h^3/12` | `error < 5.0` | 0.848 Hz vs 0.838 Hz -> 1.2% | 5% unjustified |
| `test_higher_modes` | "higher modes" | **none** | only strict monotonicity `f2 > f1`, `f3 > f2` | [0.848, 5.501, 16.542] Hz | the name promises validation; no reference value is compared. `simply_supported_plate_central_load`, `simply_supported_plate_uniform_pressure`, `nonlinear_large_displacement_cantilever` and `build_simply_supported_mesh` are all dead. |

### 5.3 `test_shell_validation_fixed.py` (7)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestLinearStatic::test_fx` | axial, 8x4 | analytical `600 L/(E b h)` | `TOL_STATIC = 3.0` (percent), justified in a comment with the measured errors | 1.29% | – |
| `test_fy` | in-plane bending | analytical `600 L^3/(3 E (h b^3/12))` | 3.0% | 1.19% | – |
| `test_fz` | out-of-plane bending | analytical `600 L^3/(3 E (b h^3/12))` | 3.0% | 1.97% | – |
| `test_ratio_physical` | uY/uX ratio | beam theory 400 (documented derivation in the docstring) | ±2%, justified against the measured 400.39 | 400.39 (0.10%) | the docstring records that the previous window was -50%/+20% |
| `test_axial_load_converges_to_the_analytical_solution` | monotone convergence, 4 meshes | analytical `P L/(E A)` | monotone decrease and finest `error < 0.01`; the exclusion of (2,1) is justified | 2.422%, 1.286%, 0.768%, 0.489% | – |
| `TestNonlinearStatic::test_geometric_nonlinearity` | "geometric nonlinearity" | none: asserts `dz_lin > L`, then `RuntimeError` matching "SNES diverged" | `dz_lin > L` | linear estimate 1.121e+02 | same shape as its `test_shell_comprehensive` twin; no numerical nonlinear reference |
| `TestModal::test_first_mode` | first cantilever mode | analytical `(1.875104^2 / 2pi) sqrt(E I/(rho A L^4))` | `error < 2.0` (percent) | 0.848 Hz vs 0.838 Hz -> 1.2% | the tolerance is the only one in this module without a justification comment, in a module whose other tolerances all carry one |

### 5.4 `test_large_rotation_benchmarks.py` (7) and `test_mitc3_benchmarks.py` (8)

Both files drive the Updated-Lagrangian incremental solve (`nonlinear_static_solve_coo` +
`PyMeshAssembler.update_reference`) on a cantilever under pure end moment, with the exact
elastica as the reference. No test in either file prints its residual, so the measured margin
column is `not printed` throughout; the assertion is what the code states.

| file / test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `large_rotation` `test_linear_tip_deflection_euler_bernoulli` | linear tip deflection and rotation, small moment | Euler-Bernoulli `P L^3/(3 E I)` | `rtol=1e-8`, `atol=1e-10` | not printed | an independent analytical reference |
| `large_rotation` `test_cantilever_large_rotation_half_circle[n_elem]` | λ = π half-circle (u_tip = -10.0, w_tip = 6.3662) | exact elastica | `tol = 0.05` relative, each component | not printed | – |
| `large_rotation` `test_equilibrium_path[lam]` | the four λ control points of `REFERENCE_TABLE` | exact elastica, tabulated in the same file | `tol_rel` for the non-zero component, `tol_abs` for the near-zero one | not printed | **weak reference** (section 9.1): `REFERENCE_TABLE` holds the rounded outputs of `_analytical_tip`, defined in the same file; the formula is correct, but a wrong formula would not be caught. |
| `large_rotation` `test_simo_vu_quoc_rollup_360[n_elem]` | λ = 2π full roll-up (tip returns to (0,0,0)) | Simo & Vu-Quoc 1986 / Bathe & Bolourchi 1979 (cited in the docstring), compared against the same `_analytical_tip` | relative for u_tip, absolute for w_tip | not printed | same weak-reference note |
| `mitc3_benchmarks` `test_linear_tip_deflection_euler_bernoulli` | same as the MITC4 case, MITC3+ mesh | Euler-Bernoulli | `rel_err < 0.02` | not printed | – |
| `mitc3_benchmarks` `test_linear_tip_moment_sign` | that a positive `M_y` produces `w_tip < 0` | sign of the physical rotation | `w_tip < 0` plus `rel_err < 0.02` | not printed | this is the MITC3 sign fix's regression guard |
| `mitc3_benchmarks` `test_cantilever_large_rotation_half_circle[n_elem]`, `test_equilibrium_path[lam]`, `test_simo_vu_quoc_rollup_360[n_elem]` | same three large-rotation cases for MITC3+ | exact elastica / `REFERENCE_TABLE` | same shape as the MITC4 file | not printed | same weak-reference note |

## 6. Element and assembly invariants

### 6.1 `test_quad_elements.py` (21)

| test group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestStiffnessMatrix::test_dimensions[Quad4/8/9]` | K shape 8x8 / 16x16 / 18x18 | plane-stress DOF count | exact | n/a | – |
| `TestStiffnessMatrix::test_symmetry` | `K == K^T` | matrix invariant | `atol=1e-6` (absolute) | not printed | absolute atol on a matrix whose entries scale with `E=210e9`; tolerance not justified |
| `TestStiffnessMatrix::test_rigid_body_modes` | exactly 3 near-zero eigenvalues | plane-stress rigid modes (tx, ty, rz) | `max_eig * 1e-10` | not printed | threshold defined in code, not justified |
| `TestStiffnessMatrix::test_positive_semi_definite` | no negative eigenvalue | PSD | `-max abs(lambda) * 1e-10` | not printed | – |
| `TestMassMatrix::test_symmetry` / `test_positive_semi_definite` | M symmetric, PSD | matrix invariants | `atol=1e-10` / `-1e-10` | not printed | – |
| `TestRigidBodyModes::test_rigid_modes_kx_ky_rz` | `K @ mode = 0` for tx, ty, rz built from node coordinates | analytical rigid-body modes `ux=-y, uy=x` | relative `max abs(residual) / norm(K)_F < 1e-10` | not printed | the most meaningful test in the file; reference is the analytical rigid-body field |

### 6.2 `test_material_suite.py` (41)

Reference model: Euler-Bernoulli with `EI` taken from the CLT `A` or `D` block, plus the
Reddy CLT B-coupling formula. No test prints its margin.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestIsotropicAnalytical` (4) | MITC4/MITC3 out-of-plane, MITC4 in-plane and axial | analytical Euler-Bernoulli with `I_out = B t^3/12`, `I_in = t B^3/12`, `P L/(E A)` | `TOL = 0.05` | not printed | 5% class tolerance, no per-case justification |
| `TestOrthotropicSinglePly` (6) | [0], [90] out-of-plane (D22), [0] and [90] axial (A22, A11), [45] ordering | analytical `P L^3/(3 D22 B)` / `P L/(A B)`, with the index choice explained per docstring | `TOL = 0.05`, [45] uses ordering plus `> 1.5x` E-B | not printed | `test_mitc4comp_ply45_out_of_plane` is a qualitative ordering check, not a numerical one; the docstring explains why E-B is invalid at 45 deg |
| `TestSymmetricLaminates` (7) | B=0 for [0/90/90/0] and quasi-iso, out-of-plane, in-plane, axial | analytical with A11/D22; quasi-iso uses `err < 0.20` for the documented ~16% D16/D26 inflation | `TOL = 0.05` (0.20 for quasi-iso) | not printed | the 0.20 is justified in a comment; `test_symmetric_b_is_zero` thresholds `max abs(B) < 1e-6` against A entries of order 1e7 — an absolute threshold with no scale argument |
| `TestAsymmetricLaminates` (6) | [0/90] B != 0, B11 < 0, MITC3Comp and MITC4Comp B-coupling vs CLT, bending-to-extension, symmetric control | Reddy CLT `w_tip = B11 F L^2/(2 A11 D11_eff b)` | `TOL = 0.10`, justified as shear-correction uncertainty | not printed; the docstring records the measured MITC3Comp `-2.627691e-03` vs CLT `-2.616014e-03` (0.45%) | the module records that the previous assertion was `abs(uy) > 1e-8`, sign-blind; the current test asserts the sign and 10% against CLT. `test_asymmetric_b_nonzero` uses `max abs(B) > 1.0` — an absolute threshold |
| `TestIsoEquivalence` (2) | 4 isotropic plies vs 1 layer | internal identity (K trace, tip displacement) | `rel < 0.01` | not printed | duplicate-path consistency check |
| `TestABDMatrices` (9) | D11 = A11 h^2/12; B exactly zero; B11 hand formula; A11 = Q11 h; D11 = Q11 h^3/12; Qbar 45/0/90 identities; Cs positive definite | closed-form CLT algebra | `1e-6` .. `1e-10` | not printed | `test_asymmetric_b11_formula` recomputes `B11_hand` from the same ply z-integrals as the implementation: **tautological** reference. `test_qbar_45_symmetry` uses `max(abs(x), 1.0)` as denominator (absolute below 1.0) and its comment promises "both positive at +45", which is never asserted. |
| `TestStiffnessProperties` (7) | global K with BC is positive definite / symmetric, 4 element types | matrix invariants | `min_eig > 0`; `max abs(K-K^T) / max abs(K) < 1e-10` | not printed | – |

### 6.3 `test_mass_matrix_validation.py` (15)

The strongest file in the suite on the reference side: all references are exact
closed-form integrals, and the module docstring explains why the translational block sums
to `rho*h*A` while the trace does not.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestElementMassVsTotalMass::test_mass_consistency[tri3/quad4]` | translational block sums to the physical mass | partition-of-unity identity `sum M = rho h A` | `rtol=1e-12` | not printed | – |
| `TestLumpedMassMatrix::test_rust_lumped_binding_matches_domain_matrix` | Rust `assemble_m_lumped` vs the domain matrix diagonal | cross-path identity | exact (`assert_allclose`, default rtol 1e-7) | not printed | – |
| `TestLumpedMassMatrix::test_lumped_vs_analytical` | row-sum lumped mass totals | `rho h A` per direction | `rtol=1e-10` | not printed | – |
| `TestModalMassConvergence::test_first_mode_mass[3 meshes]` | "first mode effective mass should converge" | analytical `f1 = 3.516/(2pi) sqrt(EI/(rho A L^4))` | `tol = 0.05 if nx <= 2 else 0.05` | not printed | three issues: the name promises a *mass*; the body asserts a *frequency*; the tolerance branch is dead (both arms 0.05); `except Exception: pytest.skip` can silently turn a solver failure into a skip |
| `TestConsistentMassTotal::test_consistent_mass_total_equals_analytical` | consistent mass total per direction | `rho h A` | `rtol=1e-12` | not printed | the docstring records the rejected `4/3 m` trace artefact |
| `TestMassMatrixSymmetry::test_symmetry` | `M - M^T` | symmetry | `max_diff < 1e-10` | not printed | absolute on an assembled matrix; justified by the fact that entries are O(1) here |
| `TestExactConsistentMassCoefficients` (4) | single-element coefficients, tri3 (1/6, 1/12) and quad4 (4/36, 2/36, 1/36), and rotary `rho h^3 A/12` | exact `int rho h N_i N_j dA` | `rtol=1e-12, atol=1e-18` | not printed | the only place in the suite that pins a *distribution* rather than a total; the docstring states why a global sum cannot catch a wrong quadrature |
| `TestConsistentMassDistribution::test_row_sums_match_tributary_areas[tri3/quad4]` | row sums equal `rho h * integral(N_i)` | shape-function tributary areas, computed by shoelace | `rel=1e-12` | not printed | – |

### 6.4 `test_rust_assembler.py` (19)

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestRustCOOAssembly::test_stiffness_rust_vs_petsc_loop[MITC3/MITC4]`, `test_mass_...` | Rust COO path vs the PETSc `setValuesLocal` fallback | the same code's independent Python path | `atol=1e-6, rtol=1e-10` | not printed | a genuine duplicate-implementation comparison; the reference is the other implementation, not a physical law |
| `...test_stiffness_symmetric`, `test_mass_symmetric` | K, M symmetric | invariants | `atol=1e-6` | not printed | – |
| `TestTangentStiffness::test_kt_at_zero_equals_k` | `KT(u=0) = K` | linear identity | `atol=1e-6, rtol=1e-10` | not printed | – |
| `TestTangentStiffness::test_kt_symmetric` | `KT(0)` symmetric | invariant | `atol=1e-6` | not printed | evaluated at `u=0` only, so the geometric part is inert |
| `TestInternalForces::test_fint_zero_at_zero` | `fint(0) = 0` | invariant | `atol=1e-10` | not printed | – |
| `TestInternalForces::test_fint_linear_equals_ku` | `fint(u, linear) = K u` for random `u ~ N(0, 1e-6)` | linear identity | `rtol=1e-6` on components above `1e-12 max` | not printed | – |
| `TestNewtonRaphsonConsistency::test_nr_consistency` | `KT(u) du = fint(u+du) - fint(u)` | Taylor consistency | `rel_err < 2e-3`, justified with the measured ~1e-3 | not printed | the only tolerance in the file with a justification |
| `TestRustGroupCoverage::test_py_mesh_assembler_built` | `assembler._rust is not None` | wiring | `is not None` | n/a | a smoke test; cannot fail given the module-level import guard |

### 6.5 `test_rust_composite.py` (24)

All references are element-matrix invariants or ordering statements; no test prints.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestMITC3BatchSanity` (8), `TestMITC4BatchSanity` (8) | batch `ke`/`me` shapes (18x18, 24x24), symmetry, PSD, per-layup symmetry (0, 45, quasi-iso, glass/epoxy, 3D tilted) | invariants | `atol=1e-6` (K), `1e-10` (M), PSD `-1e-6 max(eig)` | not printed | no external reference; the docstring states "physical sanity checks only" |
| `TestBatchComposite` (2) | 3-element batches: symmetry + PSD per element | invariants | `atol=1e-6` | not printed | – |
| `TestCompositeSanity::test_stiffer_fiber_direction` | `ke_0[0,0] > ke_90[0,0]` | ordering | strict `>` | not printed | sign/order only, magnitude unconstrained |
| `TestCompositeSanity::test_ke_positive_semidefinite` / `test_me_positive_semidefinite` (2) | `ke`/`me` eigenvalues non-negative for the `[0/90/90/0]` MITC3 element | PSD | `>= -1e-6 max(eig)` (K), `>= -1e-10` (M) | not printed | the PSD guards the section-6.5 count was missing; same class as `TestMITC3BatchSanity` |
| `test_thicker_laminate_stiffer` | `norm(ke_thick) > norm(ke_thin)` | ordering | strict `>` | not printed | ordering only |
| `test_assembler_composite_mitc3/mitc4_symmetry` | `PyMeshAssembler` composite K/M symmetry | invariants | `atol=1e-6` / `1e-10` | not printed | – |

### 6.6 `test_rust_modal.py` (10)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestMITC4ModalCantilever::test_frequencies_match` | Rust modal vs SLEPc, 6x4 clamped plate, 6 modes | `_python_modal_solve` (PETSc + SLEPc, `_aeroelast.petsc_modal_solve` vs `modal_solve_coo`) | `rtol=1e-4` | worst rel err 4.68e-12 (6 modes) | duplicate-implementation comparison |
| `TestMITC4ModalCantilever::test_mode_shapes_orthogonal` | "mode shapes are M-orthogonal" | none: asserts `norm(mode) > 1e-10` | cannot fail (an eigenvector of a GHEP is nonzero by construction) | not printed | name promises M-orthogonality; body checks a non-zero norm |
| `TestMITC3ModalCantilever::test_frequencies_match` | same for MITC3 | SLEPc | `rtol=1e-4` | max 9.00e-15 | – |
| `TestSimplySupportedPlate::test_frequencies_match_python` | pinned-edge plate, 8x8 | SLEPc | `rtol=1e-4` | max 1.32e-08 | – |
| `TestSimplySupportedPlate::test_analytical_convergence` | first mode of a pinned plate | Kirchhoff `f11 = (pi/L^2) sqrt(D/(rho h))`, `D = E h^3/(12(1-nu^2))` | `rel_err < 0.05`, "expect ~1-5% for 8x8" | 50.1553 Hz vs 49.3288 Hz = 1.68% | – |
| `TestCOOModalSolve::test_coo_matches_element_path` | `modal_solve_coo` vs the `PyMeshAssembler` path | cross-path identity | `rtol=1e-10` | not printed | – |
| `TestModalBenchmark::test_benchmark_mitc4` | 8x8 and 12x12 timing | SLEPc frequencies as the correctness guard | `rtol=1e-3` | 0.3x "speedup" printed for both sizes (Rust slower) | a timing test scored as a test; the print shows Rust losing by ~3x, which is the opposite of what the class name implies |
| `TestCompositeModal::test_composite_frequencies_match` | quasi-iso composite cantilever | SLEPc | `rtol=1e-4` | max 1.61e-10 | – |
| `TestCompositeModal::test_composite_stiffer_than_isotropic` | name says "stiffer" | none: asserts the two frequency vectors are *not* `allclose` | `rtol=0.01` inside `allclose` | first mode 0.70 Hz vs 8.72 Hz | the body asserts "different", not "stiffer" |
| `TestCompositeModalMITC3::test_composite_mitc3_frequencies_match` | same, triangles | SLEPc | `rtol=1e-4` | max 1.36e-14 | – |

### 6.7 `test_stress_stiffened_solver.py` (13)

No test prints a margin. References are PETSc matrix invariants plus the physics of a
tensile geometric stiffness.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestKGAssemblyPipeline::test_assemble_geometric_stiffness_from_stress_field` | `K_G` is a `PETSc.Mat` of the right size | type/shape | exact | n/a | – |
| `...test_K_G_symmetry` | `K_G - K_G^T ~ 0` | symmetry | FROBENIUS `1e-8` relative | not printed | relative — the same class of fix as the assembled-K check in `test_ko2017_performance` |
| `...test_K_G_is_positive_semidefinite` | tensile `K_G` PSD | PSD | `-1e-6 max abs(lambda)` | not printed | – |
| `...test_stress_recovery_element_stresses_returns_arrays` | recovery returns one value per element | array length | exact | n/a | smoke test |
| `...test_stress_field_dict_from_recovery` | non-empty `stress_field`, `K_G` builds | `len > 0` | `> 0` | not printed | cannot fail once stresses are non-zero by construction (u = 5e-3 x) |
| `...test_keff_with_KG_larger_than_without` | `sum(diag(K+KG+a0 M)) > sum(diag(K+a0 M))` | sign of stress stiffening under tension | strict `>` on a matrix-trace sum | not printed | a sign-only check; a 1e-9 relative change passes |
| `TestStressStiffenedHook::test_hook_returns_none_for_zero_displacement` | hook returns `None` at `u=0` | behavioural | `is None` | n/a | – |
| `...test_hook_returns_new_keff_under_membrane_load` | hook returns a *new* `PETSc.Mat` | behavioural | `is not` + `isinstance` | not printed | has a legitimate `pytest.skip` escape when all stresses fall below threshold |
| `...test_update_interval_skips_rebuild` | no rebuild when `step % interval != 0` | behavioural | `is None` | not printed | the second assertion is `result_5 is None or isinstance(result_5, PETSc.Mat)` — a tautology over the return type |
| `TestStressStiffenedConfig` (4) | enum value, inheritance, default and custom `update_interval` | configuration | exact | n/a | – |

### 6.8 `test_composite_b_coupling.py` (4)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_b_matrix_nonzero_for_asymmetric_laminate` | `max abs(B) > 1.0` for [0/90] | none beyond the CLT structure | absolute `> 1.0` | not printed | absolute threshold on a stiffness-scale quantity; `B` here is O(1e4) |
| `test_b_coupling_produces_bending_under_axial_load` | `w_tip < -1e-6` under axial tension, MITC4 strip | the docstring quotes `w_tip = -2.161433e-05` (CLT) and the measured `-2.1660e-05` (0.21%) | `w_tip < -1e-6` | not printed | the quoted margin (0.21%) is not asserted: the reference is commented out, only the sign and a magnitude floor are checked. The CLT formula is re-derived in the module, so the docstring's reference is the same arithmetic. Section 4.7 now adds the independent CCX S8R comparison of the same B-coupling signature (0.37%). |
| `test_symmetric_laminate_no_bending_under_axial_load` | `abs(w_tip) < 1e-9` for [0/90/90/0] | symmetry (B = 0) | `1e-9` absolute | not printed | – |
| `test_b_coupling_sign` | `B11 < 0` and `w_tip < -1e-6` | as above | `1e-6` floor | not printed | duplicates the previous test's load case; the two share ~80 lines of identical body |

### 6.9 Documentation and contract guards

These files are gates over text, not over numerics: they read the repository's own sources and
fail when a guard they pin changes. They were added by the archived change and are why section 2
can no longer silently disagree with the code.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_mitc4plusd_traceability.py` scenario 1 (1) | every ingredient row of `shell-elements.md` §2.2 names a code symbol that exists in `mitc4.rs` | the code itself | exact membership | not printed | reads files only; no Rust build |
| `...` scenario 2 (1) | the production path (before `#[cfg(test)]`) contains no drilling penalty, ERC/`beta_w`, SRI `cm_normal`, rotation bubble, `5/6` shear correction or hourglass scaffolding | the spec's forbidden set | exact absence | not printed | the non-vacuity control runs the same detectors against `mitc3.rs` and asserts they fire (where `k_drill` legitimately lives) |
| `...` scenario 3 (1) | each author-year citation resolves in `docs/references.md` and each quoted equation in the extracts | the canonical bibliography and the two extract files | exact match | not printed | – |
| `test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged` (1) | the Rust `ShellConstitutive` field set, the Rust `Laminate` method set and the Python `Laminate`/`Ply` surfaces keep their exact shape | the three declared frozensets | exact set equality | not printed | a change here is a contract change in the composite path, to be made on purpose |

## 7. Rotor and FSI group

### 7.1 `test_rotor_inertial.py` (34)

No test prints. The module loads `corotational.py` by file path to dodge the package
`__init__` (PETSc). References are closed-form rigid-body inertia and rotation algebra.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestCoordinateTransforms` (11) | identity at theta=0, orthogonality, `det R = 1`, 90deg axis mapping, round trips, magnitude preservation, axis normalisation, zero-axis error | rotation algebra | `decimal=10..12` / `atol=1e-10` | not printed | direct checks of the code under test against textbook algebra |
| `TestInertialForcesCalculator` (11) | centrifugal `m w^2 r` (direction, magnitude, zero-omega), Coriolis `-2 m (w x v)` (value, perpendicularity, zero-omega), Euler `-m (alpha x r)`, combined sum, selective inclusion, off-axis node | closed-form rigid-body inertia | `rtol=1e-10` / `atol=1e-10` | not printed | genuine formula tests; the reference is independent of the implementation |
| `TestOmegaProviders` (9) | constant/table/function providers: interpolation, extrapolation, alpha from gradient, validation errors, initial omega | provider semantics (piecewise-linear `omega`) | `rtol=0.01` for numerical alpha, exact elsewhere | not printed | `test_table_omega_interpolation` asserts `omega == 0.0` / `== 10.0` exactly at the nodes |
| `TestIntegration` (3) | `test_force_transform_and_inertial_combination` | none: `total_mag > 0` | cannot fail | not printed | the pipeline is exercised but nothing is asserted about values |
| | `test_theta_accumulation_simulation` | none: the test integrates `theta += omega*dt` itself | `rtol=0.001` vs 10 rad | not printed | **tautological**: the accumulation loop is in the test, not in the solver. Only `ConstantOmega.get_omega` is exercised. |
| | `test_displacement_consistency_over_rotation` | rotation algebra | `decimal=12` | not printed | round-trip identity of the class under test |

### 7.2 `test_rotor_physical_consistency.py` (22)

This file is the suite's worst case for tautological references. No test prints, and **no
test in this file calls the production code path it claims to protect** — the only import
is a `_aeroelast` existence probe used to decide whether to skip.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_centrifugal_deformed_geometry[16 combos]` | that a cached reference geometry produces `delta/(R+delta)` error | analytical ratio, derived in the docstring | `rel_error == expected_rel_error, abs=1e-10`; plus `> 0.09` for the 10% case | not printed (pure algebra) | **tautology**: `F_exact` and `F_cached` are both computed in the test from the same formula; `_aeroelast` is imported only to decide a skip. The physics is right; the code under test is never touched. |
| `test_kg_hysteresis_prevents_chattering` | that the 0.5%/0.3% omega^2 hysteresis gate avoids chattering | none: the gate is re-implemented in Python | `0 in rebuild_steps`, `4 not in`, `len <= 3` | not printed | **tautology**: the docstring says "The gate mirrors the Rust implementation (`rotor_fsi.rs::update_kg_if_needed`)"; the mirrored copy is what is asserted on. Assertions are structural, not numerical. |
| `test_coriolis_matrix_antisymmetry` | that `G^T = -G` | none: the matrix is hand-built in the test (`:149` "this would call `build_coriolis_matrix` in Rust") | `atol=1e-14` | not printed | **tautology**: a hand-written antisymmetric assignment is asserted to be antisymmetric. |
| `test_coriolis_implicit_stability` | that implicit Coriolis treatment keeps `K_eff` positive definite | none: `K_eff = K + a1 G + a0 M` is built in the test | `all(eigvals > 0)`, `cond < 1e6` | not printed | **tautology** for the same reason; the solver's Newmark assembly is never called. The stated stability claim is not a property of the code. |
| `test_stress_gate_checkpoint_consistency[1, 5, 10]` | that the stress gate writes stress at checkpoint steps | none: the gate (`stress_interval <= 1 or step % interval == 0 or is_checkpoint`) is re-implemented in the test with a `# <-- CRITICAL` marker | membership of checkpoint steps | not printed | **tautology**: the predicate the test asserts on is written in the test, not read from `rotor.py`. |

### 7.3 `test_rotor_rust_parity.py` (37)

No test prints. Two reference classes: (a) the Python implementation of the same helper,
compared across a stub, and (b) marshalling smoke tests whose expected outcome is "raises
something other than a TypeError/AttributeError".

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestMapOmegaProvider` (8) | that each `OmegaProvider` subclass maps to the right Rust tuple | **`_RotorStub._map_omega_provider` — a hand-copied mirror of the production method** | `assert_allclose` default (1e-7) | not printed | **tautology**: the stub's docstring says "The logic mirrors `LinearDynamicFSIRotorSolver._map_omega_provider` exactly". The test therefore validates the copy. `TestUseRustFlag::test_map_omega_round_trip_*` does call the real method, but only for 2 of 4 provider modes. |
| `TestUseRustFlag::test_use_rust_default_is_true` / `_true` / `_false_explicit` / `_truthy_int` (4) | name promises `use_rust` flag parsing | none: each asserts `not hasattr(solver, "_use_rust_fsi")` | `not hasattr` | n/a | four near-identical tests for a flag that no longer exists; the assertions cannot fail and say nothing about the four input values they are named for |
| `TestUseRustFlag::test_omega_provider_type_*` (4) | provider selection from rotor config | `_init_rotor_config` (real code) | `assert_allclose` default | not printed | genuine |
| `TestRotorAutoInertia` (11) | `_compute_estimated_inertia`, auto-inertia/auto-radius solve wiring, restart-omega resolution, `_compute_rotor_radius`, `_compute_axis_torque`, performance coefficients | hand-computed expected scalars (e.g. inertia 14.0 from masses 2 and 3 at radii 1 and 2; radius 2.0 from the coordinate set) | `assert_allclose` default | not printed | genuine, and the only part of this file where an independent hand value is compared |
| `TestRotorRustBinding::test_symbol_exists` | the symbol is callable | wiring | `callable` | n/a | – |
| `...test_call_without_precice_raises_runtime_or_os_error` | the Rust side raises rather than a marshalling error | none | any `Exception` except `TypeError`/`AttributeError` | not printed | the name promises `RuntimeError`/`OSError`; the body accepts any exception. The helper `_PRECICE_ERRORS` is defined and never used. |
| `...test_call_{ramped,computed,ramped_computed,with_callback}_marshalling` (4) | argument marshalling through the binding | none | same broad-exception filter | not printed | also verifies the callback kwarg is accepted |
| `...test_wrong_omega_mode_raises` | an unknown mode string reaches Rust | none | same | not printed | – |
| `...test_mismatched_dofs_raises` | out-of-range free DOFs | none | `BaseException` | observed: the run prints `thread '<unnamed>' panicked at aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs:241:24: index out of bounds: the len is 12 but the index is 999` | the test passes because a Rust panic surfaces as a `BaseException`; the production failure mode for a bad free-DOF array is a panic, not a clean error. Worth its own issue. |

### 7.4 `test_rotor_performance_report.py` (1) and `test_fsi_structural_report.py` (1)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_rotor_performance_report_creates_output_folder_and_csv` | that the CSV writer creates nested folders and writes the fields it is given | the injected values (`torque_aero=20`, `max_displacement=0.02`, `deformed_radius=7.5`, total-Z 6.0) | `pytest.approx` default (rel 1e-6) | not printed | a round-trip of values the test itself passed in: an I/O contract test, not a validation |
| `test_structural_report_logs_max_displacement_components` | that `max disp` and its components/owner node/position are logged | `np.linalg.norm([0.4, -0.5, 0.6])` of the array the test constructed; node `"20"` | `pytest.approx` default; exact strings for node/position | not printed | same: reference is the arithmetic re-implemented on the injected input |

## 8. BEM, aero and mesh group

### 8.1 `test_bem_polars.py` (20)

No test prints.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestPolarData` (4) | interpolation at 0 and at 5deg, periodic wrapping, vectorisation | the fixture's own `cl = 2 pi sin(alpha)` definition | `atol=0.05` (Cl), `atol=1e-6` (wrap), `< 1e-3` (zero) | not printed | **tautological reference**: `cl` is generated as `2 pi sin(alpha)` and then compared to `2 pi sin(alpha)`. The test validates the interpolator round-trip, not aerodynamics. `atol=0.05` is not justified. |
| `TestAirfoilAero` (3) | exact-Re selection, nearest-Re selection, single-polar passthrough | library semantics | exact `==` on `polar.re` | not printed | – |
| `TestBladeAeroYAML` (13) | blade length, hub radius, rotor radius, n_blades, station ordering, chord bounds, twist bounds, polars present, alpha range, `std(cl)`, and the `r`/`chord`/`twist` property shapes | the WindIO YAML that the code under test also loads | assorted absolute thresholds (`> 100`, `> 0`, `< 10`, `< 30deg`, `std > 0.1`) | not printed | the assertions are thresholds chosen by inspection, not from a source: they detect a broken parse, not a wrong value. `rotor_radius > hub_radius` and `alpha_range >= pi` are the meaningful ones. |

`tests/IEA-15-240-RWT.yaml` is the fixture both this module and `test_bem_engine.py` load;
it is a copy of the public IEA-15-240-RWT definition, not an independent reference.

### 8.2 `test_bem_engine.py` (14) and `test_blade_mesh.py` (1)

Both modules now run. The first version of this matrix recorded them as a skipped module and an
empty parameter set; that is no longer true (`CCBlade 1.3.1` is installed and the blade module
scans its own directory). Their tests do not print a residual, so the measured margin column
below reads `not printed`, not `not measured`.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestBEMSolverParked` shape/type smoke tests: `test_result_is_bem_result`, `test_r_matches_blade`, `test_Np_shape`, `test_Tp_shape`, `test_alpha_shape` (5) | result type; the BEM radial stations are the blade's; `Np`/`Tp`/`alpha` array shapes | the blade fixture's own station count and array shapes | exact type/`len`/`shape` | not printed | shape/type guards only, so they cannot detect a wrong value; they are what makes `TestBEMSolverParked` 9 tests, not 4 |
| `TestBEMSolverParked::test_parked_alpha_close_to_twist` | parked blade AoA = 90deg - twist (geometric identity) | geometry | `residual < 15deg` for `> 80%` of stations, justified as "generous" | not printed | thresholds with no derivation (`> 0.8` of stations within 15deg) |
| `TestBEMSolverParked::test_parked_thrust_positive`, `..._torque_near_zero` (`abs(power) < 1e6`), `..._Np_mostly_positive` (`> 0.7`) | sign/order sanity | none | loose thresholds | not printed | cannot detect a magnitude error |
| `TestBEMSolverRotating` (5) | rated power/thrust/torque positive, induction in range, Cl non-trivial | none | `frac_ok > 0.7`, `np.any(abs(cl) > 0.1)` | not printed | same class |
| `test_blade_mesh_generation` | that meshing a reference turbine produces nodes, elements and root/outer-shell/shear-web node sets | none | `node_count > 0`, `elements_count > 0` | not printed (the run prints `Blade mesh generated: 9277 nodes, 9867 elements`) | structural smoke test; the assertion cannot fail for a mesh with any content |

### 8.3 `test_force_projection.py` (10)

No test prints.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestForceProjectorConstruction::test_creates_strips` | one strip per BEM station | count identity | exact `== 10` | n/a | – |
| `...test_all_nodes_assigned` | every mesh node belongs to a strip | partition | exact | n/a | – |
| `TestForceConservation::test_uniform_Np_conservation` | `sum(f_nodes) = int Np dr`; y and z components zero | the trapezoidal integral of the same `Np` array, computed by `projector.verify` (production code) | `force_error < 1e-6`; `atol=1e-6` | not printed | the conservation check calls production `verify`, so a sign or factor error in both would cancel; the independent part is that the y/z sums are zero |
| `...test_uniform_Tp_conservation` | same for tangential load | as above | `< 1e-6` | not printed | – |
| `...test_combined_Np_Tp_conservation` | both components | as above | `< 1e-6` | not printed | – |
| `...test_varying_Np_conservation` | linearly varying `Np` (2000 -> 0) | trapezoidal integral | `< 50 N` on a 20 m blade, justified in a comment as strip-discretisation mismatch | not printed | the relaxed window is documented |
| `TestForceProjectionOutput::test_output_shape` | `(n_nodes, 3)` | shape | exact | n/a | – |
| `...test_zero_load_gives_zero_forces` | zero load -> zero force | linearity | `atol=1e-12` | not printed | – |
| `...test_forces_only_in_load_direction` | Np produces force only along the normal | direction | `atol=1e-8` on y/z; `sum(fx) > 0` | not printed | the comment records the measured nodal value (+250 N each, sum +6000 N) and that the sign check replaced a sign-blind one |
| `TestSingleNodeStrip::test_single_node_per_strip` | one chordwise node per strip gets the whole strip force | `F = Np dr` | `force_error < 1.0` ("relaxed for coarse discretisation") | not printed | `np.all(np.abs(forces[:,0]) > 0)` is sign-blind; the documented warning about a dropped strip moment is never asserted |

### 8.4 `test_bem_openfast_parity.py` (13)

A1 of the wind-turbine validation campaign: the AeroElast BEM (CCBlade) against
**OpenFAST 4.2.1 AeroDyn** on the IEA 15 MW. Apples-to-apples by construction: both codes get
the **same** geometry and the **same** official AeroDyn polar tables, taken from the vendored
deck `tests/reference/iea15mw_openfast/` (Apache-2.0; see its NOTICE). The helper is
`tests/_openfast_bem.py`.

OpenFAST runs through `aerodyn_driver` in combined-case mode (one invocation per session); the
AeroElast side is `BEMSolver` fed a `BladeAero` built from the same deck. The IEA primary file
ships with unsteady aerodynamics (`UA_Mod=3`) and tower influence/shadow enabled, which would
compare a quasi-steady BEM against a Beddoes-Leishman model; `write_bem_primary` turns those
off (`UA_Mod=0`, `TwrPotent=TwrShadow=0`, `TwrAero=False`, `Wake_Mod=1`, `BEM_Mod=1`) and that
is the record of the choice.

**Hub-radius note.** The AeroElast WindIO loader sets `rotor_radius = rotor_diameter/2`
(121.119 m) while its stations use `hub_radius + eta*L` (120.97 m). Both sides here are built
from the deck (`r = HubRad + BlSpn`, `HubRad = 3.97`, `TipRad = 120.97`), so that
inconsistency is not in the loop.

| test | what it validates | reference | tolerance | measured margin (V = 8 / 9 / 12 m/s) | notes |
| --- | --- | --- | --- | --- | --- |
| `test_rotor_performance_matches_aerodyn[0..2]` | rotor thrust (`RtFldFxh`) and torque (`RtFldMxh`) at 8/9/12 m/s, 7.56 rpm | **OpenFAST 4.2.1 AeroDyn**, same deck | 1.5% each | thrust +0.445% / +0.393% / +0.329%; torque -0.653% / -0.045% / +0.794% | 0 shear, 0 pitch; both margins printed |
| `test_spanwise_loads_match_aerodyn[0..2]` | interior-span angle of attack and normal coefficient (0.15R < r < 0.985R) | **OpenFAST 4.2.1 AeroDyn** | 2 deg (alpha), 2.5% mean (Cn) | max abs d alpha 0.676 / 0.707 / 1.361 deg; mean abs rel d Cn 0.85% / 0.75% / 1.69% | the exact hub and tip nodes are excluded |
| `test_yaw_and_shear_match_aerodyn[0..1]` | rotor thrust and torque with `yaw = 10 deg` and with `shear_exp = 0.2`, both azimuth-averaged | **OpenFAST 4.2.1 AeroDyn**, same deck | 1.5% each | yaw: thrust +0.334%, torque -0.274%; shear: thrust +0.186%, torque +0.327% | CCBlade integrates across azimuth, so the AeroDyn side is averaged over the last revolution (dt 0.25 s); A2 of the campaign |
| `test_bem_with_repo_default_polars_matches_aerodyn[0..2]` | integrated thrust and torque with the repository's own WindIO polars | **OpenFAST 4.2.1 AeroDyn** (official polars) | 3% band | thrust +1.178% / +1.394% / +1.968%; torque -1.364% / +0.089% / +2.117% | a **sensitivity**, not parity: A1 already shows <0.5% with identical polars |
| `test_viterna_post_stall_matches_aerodyn` | the repo's NeuralFoil + Viterna polar vs the official AeroDyn `FFA-W3-211` table (deck node at 0.9R) | **OpenFAST 4.2.1 AeroDyn** | 10% attached Cl, 45% post-stall | attached Cl within 4% (alpha 10: 1.509 vs 1.513); post-stall worst ~40% Cl at 30 deg, Cd within 17%; both reach cd_max at 90 deg | validates the repository's post-stall model; skipped when `neuralfoil` is absent |
| `test_reference_polars_are_the_official_aerodyn_tables` | the AeroElast side really parses the 50 AirfoilInfo files with 200-point tables | the deck | exact counts | 50 / 200 | non-vacuity guard against a silently empty airfoil list |

The exact hub and tip nodes are excluded from the span comparison: AeroDyn drives the axial
induction to 1 at the singular tip node and CCBlade does not, a local end-node loss convention
rather than a bulk difference. Tests skip when the OpenFAST binary is absent.

## 9. Flag summary

Grouped by defect class, with the file and the evidence. Re-read against the code at
`e879eba`: these flags are about the **tests**, not about the element's physics, and the
entries that no longer hold are marked RESOLVED in place. A flag here does not make the suite
red, and a green suite does not remove a flag.

### 9.1 Tautological references (the arithmetic under test re-implemented in the test)

| file | test | evidence |
| --- | --- | --- |
| `test_rotor_physical_consistency.py` | `test_coriolis_matrix_antisymmetry` | `:149` "this would call `build_coriolis_matrix` in Rust"; the matrix is assigned by hand and then asserted antisymmetric |
| `test_rotor_physical_consistency.py` | `test_coriolis_implicit_stability` | `K_eff` is built in the test; the solver is never called |
| `test_rotor_physical_consistency.py` | `test_kg_hysteresis_prevents_chattering` | `:71` "The gate mirrors the Rust implementation (`rotor_fsi.rs::update_kg_if_needed`)" |
| `test_rotor_physical_consistency.py` | `test_stress_gate_checkpoint_consistency` | the gate predicate is re-implemented with a `# <-- CRITICAL` comment |
| `test_rotor_physical_consistency.py` | `test_centrifugal_deformed_geometry` | `F_exact` and `F_cached` are both computed in the test; `_aeroelast` is only a skip probe |
| `test_rotor_inertial.py` | `TestIntegration::test_theta_accumulation_simulation` | `theta += omega*dt` runs in the test body |
| `test_rotor_rust_parity.py` | `TestMapOmegaProvider` (8) | `_RotorStub._map_omega_provider` is a hand-copied mirror of the production method |
| `test_bem_polars.py` | `TestPolarData::test_evaluate_at_known_alpha` | `cl` is generated as `2 pi sin(alpha)` and compared against `2 pi sin(alpha)` |
| `test_material_suite.py` | `TestABDMatrices::test_asymmetric_b11_formula` | `B11_hand` is re-derived from the same ply z-integrals the implementation uses |
| `test_fsi_structural_report.py`, `test_rotor_performance_report.py` | both tests | the reference is the norm/values the test itself injected |
| `test_large_rotation_benchmarks.py`, `test_mitc3_benchmarks.py` | `test_equilibrium_path[*]`, `test_cantilever_large_rotation_half_circle`, `test_simo_vu_quoc_rollup_360` | `REFERENCE_TABLE` holds the rounded outputs of `_analytical_tip`, the function defined in the same file; the Simo & Vu-Quoc 1986 / Bathe & Bolourchi 1979 citation in the docstring is not what the assertion compares against. The formulas are correct and the tabulated digits match them to 4 decimals, so this is weak, not wrong — but it cannot detect a wrong formula. |

### 9.2 Assertions that cannot fail

| file:line | test | assertion |
| --- | --- | --- |
| `test_shell_convergence.py:392-393` | `test_composite_laminate_gap_mesh_study` | `np.all(np.isfinite(aero)) and np.all(aero > 0)` is the entire assertion; the measured gap trend is printed and discarded |
| `test_shell_comprehensive.py:609` | `TestNonlinearStaticCantilever::test_large_displacement_tip_load` | `assert abs(dz_lin) > L` on a linear estimate of 1.121e+02 m |
| `test_shell_validation_fixed.py:280` | `TestNonlinearStatic::test_geometric_nonlinearity` | `assert dz_lin > L`, same value |
| `test_stress_stiffened_solver.py:256` | `test_stress_field_dict_from_recovery` | `len(stress_field) > 0` where the field is built from a prescribed `u = 5e-3 x` |
| `test_stress_stiffened_solver.py` | `test_update_interval_skips_rebuild` | `result_5 is None or isinstance(result_5, PETSc.Mat)` — the second disjunct exhausts the return type |
| `test_stress_stiffened_solver.py` | `test_keff_with_KG_larger_than_without` | strict `>` between two diagonal sums |
| `test_rust_modal.py:305` | `test_mode_shapes_orthogonal` | `assert norm > 1e-10` |
| `test_rotor_inertial.py:481` | `test_force_transform_and_inertial_combination` | `assert total_mag > 0` |
| `test_force_projection.py:333` | `test_single_node_per_strip` | `np.all(np.abs(forces[:, 0]) > 0)` — sign-blind |
| `test_blade_mesh.py:39-40` | `test_blade_mesh_generation` | `node_count > 0`, `elements_count > 0` |
| `test_bem_engine.py:105-123` | `TestBEMSolverRotating` | `power > 0`, `frac_ok > 0.7`, `np.any(abs(cl) > 0.1)` |
| `test_bem_polars.py:189` | `test_polar_cl_not_constant` | `std(cl) > 0.1` |
| `test_rotor_rust_parity.py` | `TestUseRustFlag` (4 tests) | `not hasattr(solver, "_use_rust_fsi")` for four different inputs |
| test_ko2017 (historical) | `test_3_6_hook_table_14_minimal_fix` | `assert norm > 0.1`; the current file records this in a comment and now asserts a 3% window |

### 9.3 Names that promise more than the body delivers

| file | test | the gap |
| --- | --- | --- |
| `test_shell_analytical_validation.py` | `TestSimplySupportedBeam::test_center_deflection` | both ends get full Dirichlet conditions (`:275`); the model is clamped-clamped |
| `test_shell_analytical_validation.py` | `TestBeamBending::test_bending_convergence` | three independent single-mesh assertions; no convergence |
| `test_shell_analytical_validation.py` | `TestCompositeMatrices::test_laminate_solve` | shapes only; no solve |
| `test_shell_comprehensive.py` | `TestModalAnalysis::test_higher_modes` | monotonicity only; no reference mode values |
| `test_shell_comprehensive.py` | `TestNonlinearStaticCantilever::test_large_displacement_tip_load` | asserts a divergence exception, not a large-displacement result |
| `test_rust_modal.py` | `test_mode_shapes_orthogonal` | M-orthogonality is never computed |
| `test_rust_modal.py` | `test_composite_stiffer_than_isotropic` | asserts "different", not "stiffer" |
| `test_rust_modal.py` | `TestModalBenchmark::test_benchmark_mitc4` | a benchmark named as a test; it prints Rust at 0.3x of Python |
| `test_mass_matrix_validation.py:615` | `TestModalMassConvergence::test_first_mode_mass` | asserts a frequency, not a mass |
| `test_rotor_rust_parity.py` | `test_call_without_precice_raises_runtime_or_os_error` | accepts any exception except `TypeError`/`AttributeError`; `_PRECICE_ERRORS` is defined and unused |
| `test_orthotropic_shell_parity.py` | `test_orthotropic_axial` | applies a transverse load |
| `test_ko2017_performance.py` | `test_3_2_circular_plate_tables_6_to_7` (SS rows) | **RESOLVED** — the SS case now carries its own Table 7 tuple (`expected_clamped` / `expected_ss` are separate parametrisations) |
| `test_ko2017_performance.py` | `test_3_5_twisted_beam_tables_12_to_13` | **RESOLVED** — the docstring and the code now agree: the four cases run with no xfail marker and `_TWISTED_BEAM_CASES` carries `xfail_reason=None` throughout, which is what this row used to flag as a mismatch |

### 9.4 Tolerances with no stated justification, tolerance/comment mismatches, dead conditionals

| file:line | item | evidence |
| --- | --- | --- |
| `test_isotropic_shell_parity.py:365-366` | comment "Allow 10% tolerance" vs `tol = 0.05` | **RESOLVED**: the comment now records that it was the wrong part; the code was 0.05 |
| `test_shell_comprehensive.py:400` | comment "Should match within 10%" vs `assert error < 5.0` | **RESOLVED**: the comment now records the same, and the internal message no longer says FY inside the FX test |
| `test_shell_comprehensive.py:544` | `0.5*ratio_ref <= ratio <= 1.2*ratio_ref` | asymmetric and unexplained; `test_shell_validation_fixed` uses ±2% for the same ratio |
| `test_shell_comprehensive.py:446,491,655` | `error < 5.0` for FY, FZ, modal | no justification comments, unlike `test_shell_validation_fixed.py` where the same physics carries a documented 3% |
| `test_shell_validation_fixed.py:313` | modal `error < 2.0` | the only unjustified tolerance in a module that justifies all others |
| `test_shell_analytical_validation.py:379,606` | `tol = 0.05` (multiple) and `0.05 if thickness_ratio < 0.01 else 0.05` | both branches identical |
| `test_mass_matrix_validation.py:615` | `tol = 0.05 if nx <= 2 else 0.05` | both branches identical |
| `test_composite_beam_parity.py` (3 tests) | `rel_error < 0.1` for 8-ply composite vs S8R | no justification for 10% |
| `test_shell_analytical_validation.py` (`TestCantileverBeam`, `TestMembraneStretching`) | `rel_error < 0.05` | no justification |
| `test_quad_elements.py` | `atol=1e-6` on K, `max_eig*1e-10` | absolute tolerances on stiffness-scale matrices, no derivation |
| `test_force_projection.py` | `< 1.0` and `< 50 N` | the 50 N is justified in a comment; the 1.0 is "relaxed for coarse discretisation" without a bound |
| `test_bem_engine.py`, `test_bem_polars.py` | assorted percentage thresholds | chosen by inspection |
| `test_material_suite.py:452,459,575` | `max abs(B) < 1e-6` / `> 1.0` | absolute thresholds on `B`, which is O(1e4) for this layup; the scale argument is missing |
| `test_material_suite.py:819,823` | `... / max(abs(x), 1.0) < 1e-8` | the `max(..., 1.0)` denominator makes the tolerance absolute for small entries and the comment promises a positivity check that is never made |

### 9.5 Informational-only comparisons that are printed and never asserted

| file | printed line | evidence |
| --- | --- | --- |
| `test_ko2017_performance.py` | `[x] Norm vs Paper 3D: ...` via `wref_paper` / `expected_paper` | read only by `_run_case`'s print; the printed errors are 523% to 6.2e8%. The comparison divides by a 3D reference and compares against a *normalized* expectation, and only `test_3_1` sets the fields. It is a broken diagnostic, not a validation. |
| `test_shell_convergence.py` | `raw-error log-log fit`, pairwise orders, `gap trend`, `DECISIVE ANSWER` | deliberately not asserted (documented in both cases) |
| `test_shell_comprehensive.py:542` | `Beam-theory ratio: 400.00` (measured 400.39) | followed by a ±(-50%/+20%) assertion, so it does carry an assertion, just a very wide one |
| `test_rust_modal.py` | per-mode tables, benchmark speedups | the frequency tables back a real `rtol=1e-4` assertion; the speedup column backs nothing |

### 9.6 Configuration drift

- `tests/test_blade_mesh.py` **RESOLVED**: it used to read only the non-existent
  `examples/reference_turbines/yamls`, so its single test never ran. It now scans its own
  directory first and raises on an empty parameter set instead of vanishing; it runs and
  passes (section 8.2).
- `tests/test_orthotropic_shell_parity.py` resolves CalculiX with its own
  `shutil.which("ccx")` instead of `conftest.ccx_bin_or_skip()`, so it does not honour
  `CCX_BIN` or the documented fallbacks.
- `test_ko2017_performance.py` writes VTK meshes on every run to
  `Path(__file__).parent.parent.parent / "output" / "test_meshes"` — i.e. *outside* the
  repository (`/home/efirvida/Desktop/dev/output/test_meshes`). It is a test side effect,
  not a fixture.
- `test_ko2017_performance.py` builds elements with the code derived from the node count,
  but every call site passes `use_triangular = False`, so the triangular path is dead in
  the suite.

## 10. Validation evidence that exists outside `tests/`

### 10.1 The element-versus-published-columns comparison

The strongest check available in this repository, because it compares the element against
the paper's own element column rather than against a rescaling of an analytical value.

Reference: **Ko, Lee, Lee & Bathe 2017, Table 10** — Scordelis-Lo roof, regular mesh,
`w_ref = 3.0240e-1` (`.sources/papers/1-s2.0-S0045794917309550-main.pdf`, journal p. 197).

Reproduced on this working tree by driving the module's own machinery with the triangular
mesh enabled, which is what the node-count-derived element code allows:

```text
PYTHONPATH=tests python -c "<build the Scordelis-Lo _Case with triangular=... and call _run_case>"
```

| element | mesh | measured (this tree) | published Table 10 | deviation |
| --- | --- | --- | --- | --- |
| MITC4+ | N=16, quads | **0.9988** | **0.9973** | **+0.15%** |
| MITC3+ | N=16, triangles | **0.9545** | **0.9550** | **-0.05%** |
| MITC3+ | N=8, triangles | 0.8561 | 0.8577 | -0.19% |
| MITC3+ | N=32, triangles | 0.9851 | 0.9851 | 0.00% |

The N=32 triangular value reproduces the published digit for digit, which is the strongest
single piece of numerical evidence in the repository for the MITC3+ implementation: the
same benchmark, the same normalization, the same reference.

Note that the *live* test (`test_3_4_scordelis_lo_tables_10_to_11`) only exercises the quad
path, so this comparison is not covered by the suite as it stands.

### 10.2 The published columns for that benchmark, for the record

Ko, Lee, Lee & Bathe 2017, Table 10 (Scordelis-Lo roof, **regular** mesh, N = 4/8/16/32/64,
reference solution `3.0240 x 10^-1`), read from the recovered PDF:

| element | N=4 | N=8 | N=16 | N=32 | N=64 |
| --- | --- | --- | --- | --- | --- |
| MITC3+ | 0.6695 | 0.8577 | 0.9550 | 0.9851 | 0.9932 |
| MITC4+ | 1.048 | 1.005 | 0.9973 | 0.9958 | 0.9960 |
| MITC4 | 0.9432 | 0.9726 | 0.9886 | 0.9936 | 0.9955 |

The same paper's Table 11 (distorted mesh) gives MITC3+ 0.7982 / 0.9223 / 0.9757 / 0.9909 /
0.9948 and MITC4+ 1.040 / 0.9973 / 0.9942 / 0.9948 / 0.9957; the live test's distorted
expectation (0.9942) is the Table 11 N=16 MITC4+ cell, so the distorted row is correctly
sourced while the pinched-cylinder and twisted-beam distorted rows are not.

### 10.3 The target for the strain-smoothed MITC3+

The element is implemented but not enabled (`docs/formulations/shell-elements.md` §4.3;
`crates/aeroelast-core/src/elements/smoothing.rs`).

- Lee, C., Lee, P.-S., "The strain-smoothed MITC3+ shell finite element", *Computers and
  Structures* 223:106096, 2019 — `lee2019.pdf`. **Table 6** ("Normalized vertical
  displacements (w/wref) at point B in the Scordelis-Lo roof shell problem when
  t/L = 1/100", Mesh I): **Smoothed MITC3+ 1.1017 / 1.0323 / 1.0075** at 4x4 / 8x8 /
  **16x16**. In the same table: MITC3+ 0.7409 / 0.8793 / 0.9618, Enriched MITC3+ 0.9610 /
  0.9931 / 0.9983, MITC4+ (Mesh I) 1.0476 / 1.0053 / 0.9977; reference solution
  `w_ref = 0.3024`.
- The target for an N=16 implementation is therefore **1.0075**, against the **0.9545**
  the unsmoothed element reaches (10.1). The gap is ~5.3% of the reference, and it is the
  quantity the smoothed variant's refinement should be measured against (the variant IS
  implemented -- `elements/smoothing.rs` plus the union layout in `mitc3.rs`; see
  `docs/formulations/shell-elements.md` section 4.3 -- so this is a comparison target, not
  pending work).
- Caveat, stated rather than assumed: Lee & Lee measure at "point B" on Mesh I, while 10.1
  and the live test measure Ko's point A on the regular mesh. The two papers agree on the
  problem definition (`L = 25`, `R = 25`, `t/L = 1/100`, self-weight 90 per unit area,
  `E = 4.32e8`, `nu = 0`, `w_ref = 0.3024`), but the measurement point and mesh pattern
  labels differ, and the published MITC3+ columns themselves differ between them
  (0.9618 Mesh I vs 0.9566 Mesh II vs Ko's 0.9550). Treat 1.0075 as the target for the
  same variant, not as a value already shown to be directly comparable.

## 11. Cross-references

- Citations and their verification status: `docs/references.md` (Ko et al. 2017 Table-10
  workbench paper, Lee & Lee 2019, Dvorkin & Bathe 1984, Hughes et al. 1977, Ko/Bathe/Zhang
  2025).
- Element formulation, DOF ordering, rotation convention and the drilling degree of freedom:
  `docs/formulations/shell-elements.md` §1 (conventions), §2 (MITC4+/D, the production shell
  element, including §2.5 for the drilling strain and §2.6 for what the element does NOT
  contain), §3 (MITC3+), §4 (limitations, history and open questions).
- Constitutive and CLT: `docs/formulations/materials.md` §2 (A, B, D, `Cs`), §3 (Q, Qbar),
  §4 (ABD to shell mapping).
- Solvers and time integration: `docs/formulations/solvers.md` §2 (Newmark-beta).
- The audit that motivated this document: `odd/tasks/test-suite-physical-correctness.md`
  (findings 1-3) and `odd/tasks/formulation-documentation.md` (unit U8).
- The wind-turbine validation references (OpenFAST, AeroDyn manual, the IEA 15 MW model,
  Bernardi, Escalera Mendoza, Luo & Gao): `docs/references.md` §4-§6.

## 12. Reproducing this matrix

```bash
source activate aeroelast-dev        # or: export PATH="$HOME/miniconda3/envs/aeroelast-dev/bin:$PATH"
cd /home/efirvida/Desktop/dev/fem-shell

# Python counts and failures (addopts= defeats the -v in pyproject.toml)
python -m pytest -q

# every printed margin in one run
python -m pytest -o addopts="" -q -s

# per-file counts used in the section headings
python -m pytest -o addopts="" --collect-only -q

# the Rust side (workspace root is crates/, not the repository root)
cargo test --manifest-path crates/Cargo.toml -p aeroelast-core
```

CCX must be on `PATH` (or `CCX_BIN` set) or the parity rows in §4.1-§4.6, §4.8 and §4.10
turn into skips silently: `conftest.ccx_bin_or_skip()` calls `pytest.skip`, so a machine
without CalculiX reports a green suite with those tests absent. §8.4 needs the OpenFAST
AeroDyn driver (`OPENFAST_BIN`, `PATH`, or the `openfast` conda env) and skips without it;
the Viterna row additionally skips when `neuralfoil` is absent. §8.4's reference deck is
vendored under `tests/reference/iea15mw_openfast/`, so it needs no network.

## 13. Consolidated results matrix

One row per matrix row of §3-§8 (a test group), with the tolerance as the code states it, the
measured margin and the **slack to fail** = tolerance − margin. `n/a` means the row prints no
numeric residual (the test asserts without printing). A **near** flag marks `slack < 1`; a
**>5%** flag marks a tolerance above the suite rule below; `not printed` is the same statement
§1 makes. Use this table to find directly which row is closest to failing.

**Suite rule: no test tolerance may exceed 5%.** Rows above it are flagged `>5%` and
justified in §13.1. Tightening one is a test change, not a document change.

| § | test | tolerance | measured margin | slack to fail | flags |
| --- | --- | --- | --- | --- | --- |
| 3. | `test_3_1_square_plate_tables_2_to_5[reg clamped, t/L=1/100..1/10000]` (3 cases) | `rtol=0.05` | 0.9996 (0.12%), 0.9980 (0.00%), 0.9979 (0.00%) | +4.88% | - |
| 3. | `...[dist clamped]` (3 cases) | `rtol=0.05` | 0.9991 (0.29%), 0.9975 (0.35%), 0.9974 (0.35%) | +4.65% | - |
| 3. | `...[reg SS]` (3 cases) | `rtol=0.05` | 1.0026 (0.26%), 0.9998 (0.00%), 0.9998 (0.00%) | +4.74% | - |
| 3. | `...[dist SS]` (3 cases) | `rtol=0.05` | 1.0074 (0.44%), 0.9998 (0.32%), 0.9996 (0.34%) | +4.56% | - |
| 3. | `test_3_2_circular_plate_tables_6_to_7[clamped, 3 t/L]` | `rtol=0.05` | 0.9971 (0.39%), 0.9967 (0.30%), 0.9967 (0.30%) | +4.61% | - |
| 3. | `...[ss, 3 t/L]` | `rtol=0.05` | 0.9975 (0.16%), 0.9974 (0.14%), 0.9974 (0.14%) | +4.84% | - |
| 3. | `test_3_3_pinched_cylinder_tables_8_to_9[reg]` | `rtol=0.05` | 0.9182 (1.41%) | +3.59% | - |
| 3. | `test_3_3_pinched_cylinder_tables_8_to_9[dist]` | `rtol=0.05` | 0.9317 (0.04%) | +4.96% | - |
| 3. | `test_3_4_scordelis_lo_tables_10_to_11[reg]` | `rtol=0.05` | 0.9701 (2.73%) | +2.27% | - |
| 3. | `test_3_4_scordelis_lo_tables_10_to_11[dist]` | `rtol=0.05` | 0.9809 (1.34%) | +3.66% | - |
| 3. | `test_3_5_twisted_beam_tables_12_to_13[0.02667, In-plane]` | `tol=0.01` | 0.9981 (0.10%) | +0.90% | near |
| 3. | `...[0.02667, Out-of-plane]` | `tol=0.01` | 0.9998 (0.25%) | +0.75% | near |
| 3. | `...[0.0002667, In-plane]` | `tol=0.01` | 0.9972 (0.06%) | +0.94% | near |
| 3. | `...[0.0002667, Out-of-plane]` | `tol=0.01` | 0.9981 (0.01%) | +0.99% | near |
| 3. | `test_3_6_hook_table_14_minimal_fix[0.9782]` | `rel_err < 0.03`, justified and measured in the comment | 0.9814 (0.33%) | +2.67% | - |
| 3. | `test_3_7_hemisphere_cutout_tables_15_to_16[reg, 4/1000]` | `rtol=0.05` | 0.9963 (0.67%) | +4.33% | - |
| 3. | `...[reg, 4/10000]` | `rtol=0.05` | 0.9816 (0.18%) | +4.82% | - |
| 3. | `...[dist, 4/1000]` | `rtol=0.05` | 1.0011 (0.53%) | +4.47% | - |
| 3. | `...[dist, 4/10000]` | `rtol=0.05` | 0.9871 (1.39%) | +3.61% | - |
| 3. | `test_3_8_full_hemisphere_table_17[4/1000]` | `rtol=0.05` | 0.9912 (0.49%) | +4.51% | - |
| 3. | `...[4/10000]` | `rtol=0.05` | 0.9772 (0.26%) | +4.74% | - |
| 3. | `test_3_9_hyperbolic_paraboloid_tables_18_to_19[reg, ...]` (2 cases) | `rtol=0.05` | 0.9681 (0.83%), 0.9681 (0.99%) | +4.01% | - |
| 3. | `...[dist, ...]` (2 cases) | `rtol=0.05` | 0.9955 (0.52%), 1.0138 (2.03%) | +2.97% | - |
| 4.1 | `test_linear_static_with_analytical[tension_axial]` | `TOL_ANALYTICAL = 0.02`, justified in a comment against measured errors | 0.612% | n/a | - |
| 4.1 | `...[compression_axial]` | 0.02 | 0.612% | n/a | - |
| 4.1 | `...[bending_fx]` | 0.02 | 0.873% | n/a | - |
| 4.1 | `...[transverse_fy]` | 0.02 | 0.749% | n/a | - |
| 4.1 | `test_linear_static_vs_ccx[4 cases]` | `rel_err <= 0.05` | 0.70% (tension), 0.70% (compression), 0.62% (bending_fx), 0.28% (transverse_fy) | +4.30% | - |
| 4.1 | `test_modal_first_five_modes` | `tol=0.05` | max rel err 7.27e-03 (0.73%) | +4.27% | - |
| 4.2 | `test_transverse_tip_displacement` | `tol = 0.05` | AE 0.149886 m, CCX 0.145747 m, ratio 0.972, difference 2.8% | +2.20% | - |
| 4.3 | `TestCompositeMaterial::test_laminate_abd_matrices` | exact keys/`==` | not applicable | n/a | - |
| 4.3 | `TestCompositeMaterial::test_mesh_connectivity` | `atol=1e-12` | not printed | n/a | not printed |
| 4.3 | `test_composite_axial_tension` | `rel_error < 0.1` | AE 45.87 um, CCX 48.09 um, 4.62% | +5.38% | >5% |
| 4.3 | `test_composite_isotropic_equiv` | `rel_error < 0.1` | AE 33.45 um, CCX 34.73 um, 3.67% | +6.33% | >5% |
| 4.3 | `test_composite_bending` | `rel_error < 0.1` | AE 1671.71 um, CCX 1700.86 um, 1.71% | +8.29% | >5% |
| 4.4 | `test_orthotropic_axial` | `rel_error < 0.05` | CCX 3 172 580.00 um, 0.54% | +4.46% | - |
| 4.4 | `test_orthotropic_bending` | `rel_error < 0.05` | AE 1266.10 um, analytical 1221.54 um, 3.6% | +1.40% | - |
| 4.4 | `test_multi_layer_iso_equivalence` | `rel_diff < 1e-4` | 9.181e-17 | n/a | - |
| 4.5 | `test_in_plane_bending_convergence` | `MIN_ORDER = 1.5` and `EXTRAPOLATED_TOL = 0.02`, both justified in the module docstring | orders 1.7314 and 1.7561 (agree, delta = 0.0247); Richardson limit 1238.1367 um vs 1230.7692 um = 0.5986% | +1.40% | - |
| 4.5 | `test_composite_laminate_gap_mesh_study` | none: the only assertions are `np.all(np.isfinite(...))` and `... > 0` | gaps -4.1765 / -1.7137 / -1.3512 / -1.7290%; verdict printed: `PLATEAUS -> formulation / ABD` | n/a | - |
| 4.6 | `test_ccx_element_types_agree_with_each_other` | `TOL_AGREEMENT = 0.02`, justified in the module docstring against the mesh study | S4 1.135840E-02, S8 1.145430E-02, S8R 1.148490E-02 m; spread (max-min)/min = 1.1137% | n/a | - |
| 4.6 | `test_ccx_element_type_matches_analytical[S4]` | `TOL_ANALYTICAL = 0.02`, justified in the module docstring | 1.135840E-02 m, 0.6140% | n/a | - |
| 4.6 | `test_ccx_element_type_matches_analytical[S8]` | 0.02 | 1.145430E-02 m, 0.2251% | n/a | - |
| 4.6 | `test_ccx_element_type_matches_analytical[S8R]` | 0.02 | 1.148490E-02 m, 0.4929% | n/a | - |
| 4.7 | `test_axial_extension_matches_ccx[5 layups]` | 2.5% | 0.10% (`uni_0`) to 1.83% (`uni_90`) | +0.67% | near |
| 4.7 | `test_transverse_bending_matches_ccx[5 layups]` | 1% | 0.11% (`asym_0_90`) to 0.51% (`uni_0`) | +0.49% | near |
| 4.7 | `test_asymmetric_b_coupling_matches_ccx` | 2% | 0.37% (`aero 1.6816e-2` vs `ccx 1.6754e-2`) | +1.63% | - |
| 4.7 | `test_symmetric_laminates_have_no_b_coupling[sym_0_90s, quasi_iso]` | 1e-11 absolute | aero ~1e-19, ccx ~1e-13 | n/a | - |
| 4.7 | `test_modal_frequencies_match_ccx[5 layups]` | 3% | worst 1.26% (`uni_0`, 8x20 mesh) | +1.74% | - |
| 4.8 | `test_blade_mass_matches_published_models` | 10% vs the article, and above the report but within 20% | **70,623 kg** = +3.7% over the article, +8.7% over the report | +11.30% | >5% |
| 4.8 | `test_blade_modal_frequencies_match_ccx[0..4]` | 10% | worst 1.65%; all five: 0.44%, 0.71%, 0.78%, 0.84%, 1.65% | +8.35% | >5% |
| 4.8 | `test_blade_first_modes_match_article[flapwise, edgewise]` | 15% | 0.526 Hz (-7.6%) and 0.702 Hz (+8.1%) | +6.90% | >5% |
| 4.8 | `test_blade_modal_frequencies_match_bernardi[0..7]` | 15% | worst 12.3% (highest pair); lower modes 1.9% / 3.4% / 4.4% / 5.3% / 5.8% / 6.8% / 9.8% | +2.70% | >5% |
| 4.8 | `test_blade_static_tip_deflection_matches_ccx` | 15% | aero 21.69 m vs ccx 22.10 m = 1.9% | +13.10% | >5% |
| 4.8 | `test_blade_static_deflection_matches_article_dlc` | 15% | 21.69 m = -7.7% | +7.30% | >5% |
| 4.9 | `test_deck_is_id_offset_invariant` | exact equality | identical | n/a | - |
| 4.9 | `test_deck_labels_are_internally_consistent` | set containment | consistent | n/a | - |
| 4.9 | `test_offset_ids_give_the_same_ccx_result` | 1e-9 relative | 0.0 | n/a | - |
| 4.10 | `test_bending_has_no_membrane_stress` | `< 2%` of the outer fibre | 0.00 MPa | n/a | - |
| 4.10 | `test_outer_fibre_stress_is_symmetric` | exact | 52.46 / 52.46 MPa | n/a | - |
| 4.10 | `test_outer_fibre_stress_matches_ccx_and_analytical` | 15% / 20% | 52.46 vs 51.64 MPa (1.58%); vs analytical 12.57% | +7.43% | >5% |
| 5.1 | `TestCantileverBeam::test_tip_deflection[4x2, 8x4, 16x8]` | `rel_error < 0.05` | not printed | n/a | not printed |
| 5.1 | `TestSimplySupportedBeam::test_center_deflection` | `rel_error < 0.05` | not printed | n/a | not printed |
| 5.1 | `TestBeamBending::test_bending_convergence[4x2, 8x4, 16x8]` | `rel_error < 0.05` | not printed | n/a | not printed |
| 5.1 | `TestMembraneStretching::test_axial_extension` | `rel_error < 0.05` | not printed | n/a | not printed |
| 5.1 | `TestCompositeMatrices::test_laminate_solve[[0], [0,90]]` | exact shapes; `> 0` for thickness | not applicable | n/a | - |
| 5.1 | `TestShearLocking::test_thin_plate_convergence[0.1]` | `tol = 0.05 if thickness_ratio < 0.01 else 0.05` — both branches identical | not printed | n/a | not printed |
| 5.2 | `TestLinearStaticCantilever::test_fx_in_plane` | `error < 5.0` (percent) | FEM 2.820e-05, ana 2.857e-05 -> 1.3% | +3.70% | - |
| 5.2 | `test_fy_in_plane` | `error < 5.0` | FEM 1.134e-02 vs ana 1.143e-02 (raw values only; no error line printed) | n/a | - |
| 5.2 | `test_fz_out_of_plane` | `error < 5.0` | FEM 1.120e+02 vs ana(bending only) 1.143e+02 (raw values only) | n/a | - |
| 5.2 | `test_in_plane_ratio_constraint` | `0.5 ratio_ref <= ratio <= 1.2 ratio_ref` | 400.39 vs 400.00 (0.10%) | n/a | - |
| 5.2 | `TestNonlinearStaticCantilever::test_large_displacement_tip_load` | `abs(dz_lin) > L` | linear estimate 1.121e+02 | n/a | - |
| 5.2 | `TestModalAnalysis::test_first_mode_frequency` | `error < 5.0` | 0.848 Hz vs 0.838 Hz -> 1.2% | +3.80% | - |
| 5.2 | `test_higher_modes` | only strict monotonicity `f2 > f1`, `f3 > f2` | [0.848, 5.501, 16.542] Hz | n/a | - |
| 5.3 | `TestLinearStatic::test_fx` | `TOL_STATIC = 3.0` (percent), justified in a comment with the measured errors | 1.29% | n/a | - |
| 5.3 | `test_fy` | 3.0% | 1.19% | +1.81% | - |
| 5.3 | `test_fz` | 3.0% | 1.97% | +1.03% | - |
| 5.3 | `test_ratio_physical` | ±2%, justified against the measured 400.39 | 400.39 (0.10%) | +1.90% | - |
| 5.3 | `test_axial_load_converges_to_the_analytical_solution` | monotone decrease and finest `error < 0.01`; the exclusion of (2,1) is justified | 2.422%, 1.286%, 0.768%, 0.489% | -1.42% | near |
| 5.3 | `TestNonlinearStatic::test_geometric_nonlinearity` | `dz_lin > L` | linear estimate 1.121e+02 | n/a | - |
| 5.3 | `TestModal::test_first_mode` | `error < 2.0` (percent) | 0.848 Hz vs 0.838 Hz -> 1.2% | +0.80% | near |
| 5.4 | `large_rotation` `test_linear_tip_deflection_euler_bernoulli` | `rtol=1e-8`, `atol=1e-10` | not printed | n/a | not printed |
| 5.4 | `large_rotation` `test_cantilever_large_rotation_half_circle[n_elem]` | `tol = 0.05` relative, each component | not printed | n/a | not printed |
| 5.4 | `large_rotation` `test_equilibrium_path[lam]` | `tol_rel` for the non-zero component, `tol_abs` for the near-zero one | not printed | n/a | not printed |
| 5.4 | `large_rotation` `test_simo_vu_quoc_rollup_360[n_elem]` | relative for u_tip, absolute for w_tip | not printed | n/a | not printed |
| 5.4 | `mitc3_benchmarks` `test_linear_tip_deflection_euler_bernoulli` | `rel_err < 0.02` | not printed | n/a | not printed |
| 5.4 | `mitc3_benchmarks` `test_linear_tip_moment_sign` | `w_tip < 0` plus `rel_err < 0.02` | not printed | n/a | not printed |
| 5.4 | `mitc3_benchmarks` `test_cantilever_large_rotation_half_circle[n_elem]`, `test_equilibrium_path[lam]`, `test_simo_vu_quoc_rollup_360[n_elem]` | same shape as the MITC4 file | not printed | n/a | not printed |
| 6.1 | `TestStiffnessMatrix::test_dimensions[Quad4/8/9]` | exact | n/a | n/a | - |
| 6.1 | `TestStiffnessMatrix::test_symmetry` | `atol=1e-6` (absolute) | not printed | n/a | not printed |
| 6.1 | `TestStiffnessMatrix::test_rigid_body_modes` | `max_eig * 1e-10` | not printed | n/a | not printed |
| 6.1 | `TestStiffnessMatrix::test_positive_semi_definite` | `-max abs(lambda) * 1e-10` | not printed | n/a | not printed |
| 6.1 | `TestMassMatrix::test_symmetry` / `test_positive_semi_definite` | `atol=1e-10` / `-1e-10` | not printed | n/a | not printed |
| 6.1 | `TestRigidBodyModes::test_rigid_modes_kx_ky_rz` | relative `max abs(residual) / norm(K)_F < 1e-10` | not printed | n/a | not printed |
| 6.2 | `TestIsotropicAnalytical` (4) | `TOL = 0.05` | not printed | n/a | not printed |
| 6.2 | `TestOrthotropicSinglePly` (6) | `TOL = 0.05`, [45] uses ordering plus `> 1.5x` E-B | not printed | n/a | not printed |
| 6.2 | `TestSymmetricLaminates` (7) | `TOL = 0.05` (0.20 for quasi-iso) | not printed | n/a | not printed |
| 6.2 | `TestAsymmetricLaminates` (6) | `TOL = 0.10`, justified as shear-correction uncertainty | not printed; the docstring records the measured MITC3Comp `-2.627691e-03` vs CLT `-2.616014e-03` (0.45%) | +9.55% | >5%, not printed |
| 6.2 | `TestIsoEquivalence` (2) | `rel < 0.01` | not printed | n/a | not printed |
| 6.2 | `TestABDMatrices` (9) | `1e-6` .. `1e-10` | not printed | n/a | not printed |
| 6.2 | `TestStiffnessProperties` (7) | `min_eig > 0`; `max abs(K-K^T) / max abs(K) < 1e-10` | not printed | n/a | not printed |
| 6.3 | `TestElementMassVsTotalMass::test_mass_consistency[tri3/quad4]` | `rtol=1e-12` | not printed | n/a | not printed |
| 6.3 | `TestLumpedMassMatrix::test_rust_lumped_binding_matches_domain_matrix` | exact (`assert_allclose`, default rtol 1e-7) | not printed | n/a | not printed |
| 6.3 | `TestLumpedMassMatrix::test_lumped_vs_analytical` | `rtol=1e-10` | not printed | n/a | not printed |
| 6.3 | `TestModalMassConvergence::test_first_mode_mass[3 meshes]` | `tol = 0.05 if nx <= 2 else 0.05` | not printed | n/a | not printed |
| 6.3 | `TestConsistentMassTotal::test_consistent_mass_total_equals_analytical` | `rtol=1e-12` | not printed | n/a | not printed |
| 6.3 | `TestMassMatrixSymmetry::test_symmetry` | `max_diff < 1e-10` | not printed | n/a | not printed |
| 6.3 | `TestExactConsistentMassCoefficients` (4) | `rtol=1e-12, atol=1e-18` | not printed | n/a | not printed |
| 6.3 | `TestConsistentMassDistribution::test_row_sums_match_tributary_areas[tri3/quad4]` | `rel=1e-12` | not printed | n/a | not printed |
| 6.4 | `TestRustCOOAssembly::test_stiffness_rust_vs_petsc_loop[MITC3/MITC4]`, `test_mass_...` | `atol=1e-6, rtol=1e-10` | not printed | n/a | not printed |
| 6.4 | `...test_stiffness_symmetric`, `test_mass_symmetric` | `atol=1e-6` | not printed | n/a | not printed |
| 6.4 | `TestTangentStiffness::test_kt_at_zero_equals_k` | `atol=1e-6, rtol=1e-10` | not printed | n/a | not printed |
| 6.4 | `TestTangentStiffness::test_kt_symmetric` | `atol=1e-6` | not printed | n/a | not printed |
| 6.4 | `TestInternalForces::test_fint_zero_at_zero` | `atol=1e-10` | not printed | n/a | not printed |
| 6.4 | `TestInternalForces::test_fint_linear_equals_ku` | `rtol=1e-6` on components above `1e-12 max` | not printed | n/a | not printed |
| 6.4 | `TestNewtonRaphsonConsistency::test_nr_consistency` | `rel_err < 2e-3`, justified with the measured ~1e-3 | not printed | n/a | not printed |
| 6.4 | `TestRustGroupCoverage::test_py_mesh_assembler_built` | `is not None` | n/a | n/a | - |
| 6.5 | `TestMITC3BatchSanity` (8), `TestMITC4BatchSanity` (8) | `atol=1e-6` (K), `1e-10` (M), PSD `-1e-6 max(eig)` | not printed | n/a | not printed |
| 6.5 | `TestBatchComposite` (2) | `atol=1e-6` | not printed | n/a | not printed |
| 6.5 | `TestCompositeSanity::test_stiffer_fiber_direction` | strict `>` | not printed | n/a | not printed |
| 6.5 | `TestCompositeSanity::test_ke_positive_semidefinite` / `test_me_positive_semidefinite` (2) | `>= -1e-6 max(eig)` (K), `>= -1e-10` (M) | not printed | n/a | not printed |
| 6.5 | `test_thicker_laminate_stiffer` | strict `>` | not printed | n/a | not printed |
| 6.5 | `test_assembler_composite_mitc3/mitc4_symmetry` | `atol=1e-6` / `1e-10` | not printed | n/a | not printed |
| 6.6 | `TestMITC4ModalCantilever::test_frequencies_match` | `rtol=1e-4` | worst rel err 4.68e-12 (6 modes) | n/a | - |
| 6.6 | `TestMITC4ModalCantilever::test_mode_shapes_orthogonal` | cannot fail (an eigenvector of a GHEP is nonzero by construction) | not printed | n/a | not printed |
| 6.6 | `TestMITC3ModalCantilever::test_frequencies_match` | `rtol=1e-4` | max 9.00e-15 | n/a | - |
| 6.6 | `TestSimplySupportedPlate::test_frequencies_match_python` | `rtol=1e-4` | max 1.32e-08 | n/a | - |
| 6.6 | `TestSimplySupportedPlate::test_analytical_convergence` | `rel_err < 0.05`, "expect ~1-5% for 8x8" | 50.1553 Hz vs 49.3288 Hz = 1.68% | +3.32% | - |
| 6.6 | `TestCOOModalSolve::test_coo_matches_element_path` | `rtol=1e-10` | not printed | n/a | not printed |
| 6.6 | `TestModalBenchmark::test_benchmark_mitc4` | `rtol=1e-3` | 0.3x "speedup" printed for both sizes (Rust slower) | n/a | - |
| 6.6 | `TestCompositeModal::test_composite_frequencies_match` | `rtol=1e-4` | max 1.61e-10 | n/a | - |
| 6.6 | `TestCompositeModal::test_composite_stiffer_than_isotropic` | `rtol=0.01` inside `allclose` | first mode 0.70 Hz vs 8.72 Hz | n/a | - |
| 6.6 | `TestCompositeModalMITC3::test_composite_mitc3_frequencies_match` | `rtol=1e-4` | max 1.36e-14 | n/a | - |
| 6.7 | `TestKGAssemblyPipeline::test_assemble_geometric_stiffness_from_stress_field` | exact | n/a | n/a | - |
| 6.7 | `...test_K_G_symmetry` | FROBENIUS `1e-8` relative | not printed | n/a | not printed |
| 6.7 | `...test_K_G_is_positive_semidefinite` | `-1e-6 max abs(lambda)` | not printed | n/a | not printed |
| 6.7 | `...test_stress_recovery_element_stresses_returns_arrays` | exact | n/a | n/a | - |
| 6.7 | `...test_stress_field_dict_from_recovery` | `> 0` | not printed | n/a | not printed |
| 6.7 | `...test_keff_with_KG_larger_than_without` | strict `>` on a matrix-trace sum | not printed | n/a | not printed |
| 6.7 | `TestStressStiffenedHook::test_hook_returns_none_for_zero_displacement` | `is None` | n/a | n/a | - |
| 6.7 | `...test_hook_returns_new_keff_under_membrane_load` | `is not` + `isinstance` | not printed | n/a | not printed |
| 6.7 | `...test_update_interval_skips_rebuild` | `is None` | not printed | n/a | not printed |
| 6.7 | `TestStressStiffenedConfig` (4) | exact | n/a | n/a | - |
| 6.8 | `test_b_matrix_nonzero_for_asymmetric_laminate` | absolute `> 1.0` | not printed | n/a | not printed |
| 6.8 | `test_b_coupling_produces_bending_under_axial_load` | `w_tip < -1e-6` | not printed | n/a | not printed |
| 6.8 | `test_symmetric_laminate_no_bending_under_axial_load` | `1e-9` absolute | not printed | n/a | not printed |
| 6.8 | `test_b_coupling_sign` | `1e-6` floor | not printed | n/a | not printed |
| 6.9 | `test_mitc4plusd_traceability.py` scenario 1 (1) | exact membership | not printed | n/a | not printed |
| 6.9 | `...` scenario 2 (1) | exact absence | not printed | n/a | not printed |
| 6.9 | `...` scenario 3 (1) | exact match | not printed | n/a | not printed |
| 6.9 | `test_laminate_invariant_guard.py::test_laminate_public_surface_unchanged` (1) | exact set equality | not printed | n/a | not printed |
| 7.1 | `TestCoordinateTransforms` (11) | `decimal=10..12` / `atol=1e-10` | not printed | n/a | not printed |
| 7.1 | `TestInertialForcesCalculator` (11) | `rtol=1e-10` / `atol=1e-10` | not printed | n/a | not printed |
| 7.1 | `TestOmegaProviders` (9) | `rtol=0.01` for numerical alpha, exact elsewhere | not printed | n/a | not printed |
| 7.1 | `TestIntegration` (3) | cannot fail | not printed | n/a | not printed |
| 7.1 | - | `rtol=0.001` vs 10 rad | not printed | n/a | not printed |
| 7.1 | - | `decimal=12` | not printed | n/a | not printed |
| 7.2 | `test_centrifugal_deformed_geometry[16 combos]` | `rel_error == expected_rel_error, abs=1e-10`; plus `> 0.09` for the 10% case | not printed (pure algebra) | n/a | >5%, not printed |
| 7.2 | `test_kg_hysteresis_prevents_chattering` | `0 in rebuild_steps`, `4 not in`, `len <= 3` | not printed | n/a | not printed |
| 7.2 | `test_coriolis_matrix_antisymmetry` | `atol=1e-14` | not printed | n/a | not printed |
| 7.2 | `test_coriolis_implicit_stability` | `all(eigvals > 0)`, `cond < 1e6` | not printed | n/a | not printed |
| 7.2 | `test_stress_gate_checkpoint_consistency[1, 5, 10]` | membership of checkpoint steps | not printed | n/a | not printed |
| 7.3 | `TestMapOmegaProvider` (8) | `assert_allclose` default (1e-7) | not printed | n/a | not printed |
| 7.3 | `TestUseRustFlag::test_use_rust_default_is_true` / `_true` / `_false_explicit` / `_truthy_int` (4) | `not hasattr` | n/a | n/a | - |
| 7.3 | `TestUseRustFlag::test_omega_provider_type_*` (4) | `assert_allclose` default | not printed | n/a | not printed |
| 7.3 | `TestRotorAutoInertia` (11) | `assert_allclose` default | not printed | n/a | not printed |
| 7.3 | `TestRotorRustBinding::test_symbol_exists` | `callable` | n/a | n/a | - |
| 7.3 | `...test_call_without_precice_raises_runtime_or_os_error` | any `Exception` except `TypeError`/`AttributeError` | not printed | n/a | not printed |
| 7.3 | `...test_call_{ramped,computed,ramped_computed,with_callback}_marshalling` (4) | same broad-exception filter | not printed | n/a | not printed |
| 7.3 | `...test_wrong_omega_mode_raises` | same | not printed | n/a | not printed |
| 7.3 | `...test_mismatched_dofs_raises` | `BaseException` | observed: the run prints `thread '<unnamed>' panicked at aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs:241:24: index out of bounds: the len is 12 but the index is 999` | n/a | - |
| 7.4 | `test_rotor_performance_report_creates_output_folder_and_csv` | `pytest.approx` default (rel 1e-6) | not printed | n/a | not printed |
| 7.4 | `test_structural_report_logs_max_displacement_components` | `pytest.approx` default; exact strings for node/position | not printed | n/a | not printed |
| 8.1 | `TestPolarData` (4) | `atol=0.05` (Cl), `atol=1e-6` (wrap), `< 1e-3` (zero) | not printed | n/a | not printed |
| 8.1 | `TestAirfoilAero` (3) | exact `==` on `polar.re` | not printed | n/a | not printed |
| 8.1 | `TestBladeAeroYAML` (13) | assorted absolute thresholds (`> 100`, `> 0`, `< 10`, `< 30deg`, `std > 0.1`) | not printed | n/a | not printed |
| 8.2 | `TestBEMSolverParked` shape/type smoke tests: `test_result_is_bem_result`, `test_r_matches_blade`, `test_Np_shape`, `test_Tp_shape`, `test_alpha_shape` (5) | exact type/`len`/`shape` | not printed | n/a | not printed |
| 8.2 | `TestBEMSolverParked::test_parked_alpha_close_to_twist` | `residual < 15deg` for `> 80%` of stations, justified as "generous" | not printed | n/a | >5%, not printed |
| 8.2 | `TestBEMSolverParked::test_parked_thrust_positive`, `..._torque_near_zero` (`abs(power) < 1e6`), `..._Np_mostly_positive` (`> 0.7`) | loose thresholds | not printed | n/a | not printed |
| 8.2 | `TestBEMSolverRotating` (5) | `frac_ok > 0.7`, `np.any(abs(cl) > 0.1)` | not printed | n/a | not printed |
| 8.2 | `test_blade_mesh_generation` | `node_count > 0`, `elements_count > 0` | not printed (the run prints `Blade mesh generated: 9277 nodes, 9867 elements`) | n/a | not printed |
| 8.3 | `TestForceProjectorConstruction::test_creates_strips` | exact `== 10` | n/a | n/a | - |
| 8.3 | `...test_all_nodes_assigned` | exact | n/a | n/a | - |
| 8.3 | `TestForceConservation::test_uniform_Np_conservation` | `force_error < 1e-6`; `atol=1e-6` | not printed | n/a | not printed |
| 8.3 | `...test_uniform_Tp_conservation` | `< 1e-6` | not printed | n/a | not printed |
| 8.3 | `...test_combined_Np_Tp_conservation` | `< 1e-6` | not printed | n/a | not printed |
| 8.3 | `...test_varying_Np_conservation` | `< 50 N` on a 20 m blade, justified in a comment as strip-discretisation mismatch | not printed | n/a | not printed |
| 8.3 | `TestForceProjectionOutput::test_output_shape` | exact | n/a | n/a | - |
| 8.3 | `...test_zero_load_gives_zero_forces` | `atol=1e-12` | not printed | n/a | not printed |
| 8.3 | `...test_forces_only_in_load_direction` | `atol=1e-8` on y/z; `sum(fx) > 0` | not printed | n/a | not printed |
| 8.3 | `TestSingleNodeStrip::test_single_node_per_strip` | `force_error < 1.0` ("relaxed for coarse discretisation") | not printed | n/a | not printed |
| 8.4 | `test_rotor_performance_matches_aerodyn[0..2]` | 1.5% each | thrust +0.445% / +0.393% / +0.329%; torque -0.653% / -0.045% / +0.794% | +0.71% | near |
| 8.4 | `test_spanwise_loads_match_aerodyn[0..2]` | 2 deg (alpha), 2.5% mean (Cn) | max abs d alpha 0.676 / 0.707 / 1.361 deg; mean abs rel d Cn 0.85% / 0.75% / 1.69% | +0.81% | near |
| 8.4 | `test_yaw_and_shear_match_aerodyn[0..1]` | 1.5% each | yaw: thrust +0.334%, torque -0.274%; shear: thrust +0.186%, torque +0.327% | +1.17% | - |
| 8.4 | `test_bem_with_repo_default_polars_matches_aerodyn[0..2]` | 3% band | thrust +1.178% / +1.394% / +1.968%; torque -1.364% / +0.089% / +2.117% | +0.88% | near |
| 8.4 | `test_viterna_post_stall_matches_aerodyn` | 10% attached Cl, 45% post-stall | attached Cl within 4% (alpha 10: 1.509 vs 1.513); post-stall worst ~40% Cl at 30 deg, Cd within 17%; both reach cd_max at 90 deg | +5.00% | >5% |
| 8.4 | `test_reference_polars_are_the_official_aerodyn_tables` | exact counts | 50 / 200 | n/a | - |

### 13.1 Tolerance audit (> 5%)

The rows whose stated tolerance is above the 5% rule, with the measured margin and the reason.
Every other row is at or below 5%.

| § | test | tolerance | measured margin | why above 5% |
| --- | --- | --- | --- | --- |
| 4.3 | `test_composite_axial_tension` | `rel_error < 0.1` | AE 45.87 um, CCX 48.09 um, 4.62% | no stated justification (§9.4) |
| 4.3 | `test_composite_isotropic_equiv` | `rel_error < 0.1` | AE 33.45 um, CCX 34.73 um, 3.67% | no stated justification (§9.4) |
| 4.3 | `test_composite_bending` | `rel_error < 0.1` | AE 1671.71 um, CCX 1700.86 um, 1.71% | no stated justification (§9.4) |
| 4.8 | `test_blade_mass_matches_published_models` | 10% vs the article, and above the report but within 20% | **70,623 kg** = +3.7% over the article, +8.7% over the report | mass vs two published models; asserted as a sign (+3.7%) |
| 4.8 | `test_blade_modal_frequencies_match_ccx[0..4]` | 10% | worst 1.65%; all five: 0.44%, 0.71%, 0.78%, 0.84%, 1.65% | same shell mesh in CCX; measured 1.65% |
| 4.8 | `test_blade_first_modes_match_article[flapwise, edgewise]` | 15% | 0.526 Hz (-7.6%) and 0.702 Hz (+8.1%) | beam (BModes) vs shell; measured 8.1% |
| 4.8 | `test_blade_modal_frequencies_match_bernardi[0..7]` | 15% | worst 12.3% (highest pair); lower modes 1.9% / 3.4% / 4.4% / 5.3% / 5.8% / 6.8% / 9.8% | beam-based CSD vs shell; measured 12.3% |
| 4.8 | `test_blade_static_tip_deflection_matches_ccx` | 15% | aero 21.69 m vs ccx 22.10 m = 1.9% | measured 1.9%, so the window is loose |
| 4.8 | `test_blade_static_deflection_matches_article_dlc` | 15% | 21.69 m = -7.7% | proxy load (DLC 1.4 is aero-elastic); measured 7.7% |
| 4.10 | `test_outer_fibre_stress_matches_ccx_and_analytical` | 15% / 20% | 52.46 vs 51.64 MPa (1.58%); vs analytical 12.57% | coarse 8x2 linear mesh; measured 12.57% vs analytical |
| 6.2 | `TestAsymmetricLaminates` (6) | `TOL = 0.10`, justified as shear-correction uncertainty | not printed; the docstring records the measured MITC3Comp `-2.627691e-03` vs CLT `-2.616014e-03` (0.45%) | justified as shear-correction uncertainty (§6.2) |
| 8.4 | `test_viterna_post_stall_matches_aerodyn` | 10% attached Cl, 45% post-stall | attached Cl within 4% (alpha 10: 1.509 vs 1.513); post-stall worst ~40% Cl at 30 deg, Cd within 17%; both reach cd_max at 90 deg | post-stall model difference; preprint reference |

Two rows were **excluded by inspection**: `test_centrifugal_deformed_geometry[16 combos]`
(§7.2, its `abs=1e-10` is an absolute algebra tolerance, not a percentage) and
`TestBEMSolverParked::test_parked_alpha_close_to_twist` (§8.2, `residual < 15 deg` on `> 80%`
of stations — the 80% is a station fraction, not a tolerance). Both are flagged in §9 instead.
So the rule is exceeded by **12 rows**, all in the CCX-parity and BEM-parity families where the
reference itself is a different model (beam vs shell, or a proxy load).

