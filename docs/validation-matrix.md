# Validation matrix

One auditable row per test (or per coherent parametrised group) of the suite under
`tests/`, recording **what each test validates**, **the exact reference it validates
against**, **the tolerance as the code states it**, and **the margin the run actually
achieved**.

This document exists because the suite scatters reference values and tolerances
across 29 files: each test hardcodes its own expectation, so a wrong benchmark value
can survive inside an assertion that cannot fail. A single matrix makes every
reference and every margin visible.

## 1. How this matrix was produced

- Working tree: `c7f2bbc`, clean (this refresh). The matrix was first produced at `b2c62ff`,
  which is hundreds of commits behind; the per-row margins below were captured then and have
  not all been re-measured since, which is stated again in section 2.
- Interpreter: `~/miniconda3/envs/aeroelast-dev/bin/python`.
- CalculiX: **2.23** from `~/miniconda3/envs/aeroelast-dev/bin/ccx`, on `PATH`, so the
  CCX parity tests ran instead of skipping.
- Suite run:
  `python -c "import pytest,sys; sys.exit(pytest.main(['tests','-o','addopts=','-q','--tb=no']))"`
  -> `364 passed, 2 skipped, 1 warning in 138.09s`.
  (The `-o addopts=` is needed because `pyproject.toml` adds `-v`.)
- Printed margins were read from a second pass per file using `-s`; the number recorded
  in the *measured margin* column is the one the test itself printed.
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
- **Notes** — flags per section 5.

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

Measured at `c7f2bbc` with `python -m pytest -q` -> **386 passed, 0 failed, 0 skipped**,
and reproduced by `python -m pytest tests --collect-only -q` -> `386 tests collected`.

Two skips that the first version of this matrix recorded as verified are **resolved**, and
the two tests they hid now run:

| file | then | now |
| --- | --- | --- |
| `test_bem_engine.py` | 0 collected, module skipped (`ccblade` missing) | **14 passed** — CCBlade 1.3.1 is installed in the environment |
| `test_blade_mesh.py` | 0 collected, empty parameter set | **1 passed** — the module scanned a directory that does not exist and never ran; it now scans its own directory and raises if the parameter set is empty |

Three files were added by the change that produced this refresh: `test_mitc4plusd_traceability.py`
(3, the documentation-to-code gate) and `test_laminate_invariant_guard.py` (1, the composite
surface guard).

**Known drift, stated rather than papered over.** The per-file table that used to stand here was
captured at `b2c62ff` and has drifted as the suite changed: `test_rust_assembler.py` is 19 now,
not the 18 it recorded, because the MITC4 case of `test_kt_at_zero_equals_k` was removed while
the file gained tests. Rather than restate a table that would drift again, the authoritative
number is the total above, which `pytest --collect-only -q` reproduces in seconds. The per-file
sections from section 3 onward remain the audit trail of the run that first built them.

## 3. `tests/test_ko2017_performance.py` (Ko, Lee, Lee & Bathe 2017)

All 31 cases print `Norm vs Kirchhoff: <value> (expected: <cell>, error: <x>%)`. Every
hardcoded `expected_normalized` was checked against the cited table cell of the recovered
PDF; the *reference* column names the cell, and "cell match" records whether the hardcoded
number is that cell.

Common tolerance: `rtol=0.05` (tests 3.1, 3.2, 3.3, 3.4, 3.6, 3.7, 3.8, 3.9 via
`assert_relative_error`) and per-case `tol` in 3.5.

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_3_1_square_plate_tables_2_to_5[reg, t/L=1/100..1/10000]` (6 cases: 3 t/L x distorted False) | clamped square plate centre deflection, N=16 quarter model, uniform pressure | Ko 2017 **Table 2** (regular), MITC4 column, N=16: 0.9984 / 0.9980 / 0.9979; `w_ref` alpha = 1.267e-3 in `p L^4/D` | `np.isclose(..., rtol=0.05)` | 1.0029 (0.45%), 1.0004 (0.24%), 0.9981 (0.02%) | cell match |
| `test_3_1_square_plate_tables_2_to_5[dist, ...]` (3 cases) | same, distorted mesh | Ko 2017 **Table 3**, MITC4 N=16: 1.002 / 1.001 / 1.001 | `rtol=0.05` | 1.0028 (0.08%), 0.9998 (0.12%), 0.9978 (0.32%) | cell match |
| `test_3_1_square_plate_tables_2_to_5[reg SS]` (3 cases) | simply supported square plate centre deflection | **Table 4**, MITC4 N=16: 1.000 / 0.9998 / 0.9998; alpha = 4.062e-3 | `rtol=0.05` | 1.0057 (0.57%), 1.0004 (0.06%), 0.9998 (0.00%) | cell match |
| `test_3_1_square_plate_tables_2_to_5[dist SS]` (3 cases) | same, distorted | **Table 5**, MITC4 N=16: 1.003 / 1.003 / 1.003 | `rtol=0.05` | 1.0094 (0.64%), 1.0009 (0.21%), 0.9997 (0.33%) | cell match |
| `test_3_2_circular_plate_tables_6_to_7[clamped, 3 t/L]` | clamped circular plate centre deflection, N=16 quarter disk | **Table 6** MITC4 N=16: 1.001 / 0.9997 / 0.9997; `w_ref = alpha p R^4/D`, alpha = 1/64 | `rtol=0.05` | 0.9998 (0.12%), 0.9980 (0.17%), 0.9968 (0.29%) | cell match |
| `test_3_2_circular_plate_tables_6_to_7[ss, 3 t/L]` | simply supported circular plate | **Table 7** MITC4 N=16: 0.9991 / 0.9988 / 0.9988; alpha = (5+nu)/(64(1+nu)) | `rtol=0.05` | 0.9982 (0.28%), 0.9976 (0.21%), 0.9974 (0.23%) | **cell mismatch**: the same `expected_mitc4` tuple (Table 6 values 1.001 / 0.9997 / 0.9997) is passed to both the clamped and the simply supported case, so the SS rows are compared against the clamped column. The name promises Tables 6 *and* 7. |
| `test_3_3_pinched_cylinder_tables_8_to_9[reg]` | pinched cylinder load-point displacement, N=16, R=300, t=3 | **Table 8** MITC4+ N=16: 0.9313; `w_ref = 1.8248e-5` | `assert_relative_error tol=0.05` | 0.9674 (3.88%) | cell match; only 1.12 pp inside a 5% window |
| `test_3_3_pinched_cylinder_tables_8_to_9[dist]` | same, distorted mesh | **Table 9** MITC4+ N=16: **0.9321** | `tol=0.05` | 0.9932 (0.41% from the hardcoded value) | **cell mismatch**: the hardcoded 0.9892 is not a Table 9 value. `0.9892` occurs in this paper only in Table 12 (twisted beam, in-plane, t/L=0.02667, N=2, MITC4). Against the true Table 9 N=16 cell the measured value is 6.6% off, i.e. outside the 5% window. |
| `test_3_4_scordelis_lo_tables_10_to_11[reg]` | Scordelis-Lo roof free-edge centre, N=16, R=25, L=50, th=40deg, t=0.25, rho=360 | **Table 10** MITC4+ N=16: 0.9973; `w_ref = 3.0240e-1` | `tol=0.05` | 0.9988 (0.15%) | cell match |
| `test_3_4_scordelis_lo_tables_10_to_11[dist]` | same, distorted | **Table 11** MITC4+ N=16: 0.9942 | `tol=0.05` | 1.0006 (0.64%) | cell match |
| `test_3_5_twisted_beam_tables_12_to_13[0.02667, In-plane]` | MacNeal-Harder twisted beam, N=16x96, tip centre, in-plane | **Table 12** MITC4+ N=16: **0.9971**; `w_ref = 5.4240e-3` | `tol=0.10` | 1.0001 (1.95% from the hardcoded 1.02) | **cell mismatch**: hardcoded 1.02 appears nowhere in Tables 12/13. The 10% window is what absorbs the difference. |
| `test_3_5_twisted_beam_tables_12_to_13[0.02667, Out-of-plane]` | same, out-of-plane | **Table 13** MITC4+ N=16: 0.9973; `w_ref = 1.7540e-3` | `tol=0.05` | 0.9999 (1.00%) | **cell mismatch** (hardcoded 0.99; 0.99 is not a Table 13 cell either) |
| `test_3_5_twisted_beam_tables_12_to_13[0.0002667, In-plane]` | thin twisted beam, in-plane | **Table 12** MITC4+ N=16: **0.9978**; `w_ref = 5.2560e-3` | `tol=0.05` | 0.9131 (0.75% from the hardcoded 0.92) | **cell mismatch, expectation tuned to the implementation**: the comment says "MITC4+ achieves ~91% of reference", and 0.92 is that measurement, not the paper's 0.9978. Against the paper cell this case fails by 8.5%, i.e. the window would have to be ~19% wide. Same defect class as the historical Hook 1.12. |
| `test_3_5_twisted_beam_tables_12_to_13[0.0002667, Out-of-plane]` | thin twisted beam, out-of-plane | **Table 13** MITC4+ N=16: **0.9982**; `w_ref = 1.2940e-3` | `tol=0.05` | 0.9112 (0.95%) | **cell mismatch, expectation tuned to the implementation** (hardcoded 0.92) |
| `test_3_6_hook_table_14_minimal_fix[0.9782]` | Raasch hook tip deflection, N=8x48, t=2, E=3.3e3 | **Table 14** MITC4+ N=8: 0.9782; `w_ref = 4.82482` (MITC9, N=64) | `rel_err < 0.03`, justified and measured in the comment | 0.9927 (1.48%) | cell match; the reference and the 3% window are documented, and the previous 1.12 / `assert norm > 0.1` are recorded in the comment |
| `test_3_7_hemisphere_cutout_tables_15_to_16[reg, 4/1000]` | hemisphere with cut-out, regular mesh, N=16, R=10, th0=18deg | **Table 15** MITC4+ N=16: **1.003** | `tol=0.05` | 1.0047 (0.43% from the hardcoded 1.009) | **cell mismatch**: 1.009 is the Table 15 MITC4+ **N=8** cell. Against the correct N=16 cell the margin is 0.17%. |
| `test_3_7_hemisphere_cutout_tables_15_to_16[reg, 4/10000]` | same, thin | **Table 15** MITC4+ N=16: **0.9834** | `tol=0.05` | 0.9786 (0.25% from the hardcoded 0.9811) | **cell mismatch**: 0.9811 is the Table 15 **S4** N=16 cell. Against MITC4+ N=16 the margin is 0.49%. |
| `test_3_7_hemisphere_cutout_tables_15_to_16[dist, ...]` (2 cases) | same, distorted mesh | **Table 16** MITC4+ N=16: 0.9958 / 0.9736 | `tol=0.05` | 1.0042 (0.84%), 0.9641 (0.98%) | cell match |
| `test_3_8_full_hemisphere_table_17[4/1000]` | full hemisphere, N=16, th0=2deg | **Table 17** MITC4+ N=16: 0.9960; `w_ref = 9.24e-2` | `tol=0.05` | 0.9982 (0.22%) | cell match; the mesh skips the pole (th0 = 2 deg) whereas the paper's mesh is not stated to do so |
| `test_3_8_full_hemisphere_table_17[4/10000]` | same, thin | **Table 17** MITC4+ N=16: 0.9798 | `tol=0.05` | 0.9720 (0.80%) | cell match |
| `test_3_9_hyperbolic_paraboloid_tables_18_to_19[reg, ...]` (2 cases) | hyperbolic paraboloid free-edge centre, 32x32 on the full saddle | **Table 18** MITC4+ N=16: 0.9762 / 0.9777 | `tol=0.05` | 0.9761 (0.01%), 0.9770 (0.07%) | cell match |
| `test_3_9_hyperbolic_paraboloid_tables_18_to_19[dist, ...]` (2 cases) | same, distorted | **Table 19** MITC4+ N=16: 0.9904 / 0.9936 | `tol=0.05` | 0.9975 (0.72%), 1.0168 (2.34%) | cell match |

Structural notes on this module:

- `_assemble_global` derives the element code from the node count
  (`3 if len(node_ids_elem) == 3 else 4`), which is what allows the same benchmark to be
  run with quads or triangles; every live call site currently passes
  `use_triangular = False`, so the triangular path is unreachable from the suite (section 6
  exercises it out of band).
- The module-level `REFERENCE_VALUES` table named in the audit no longer exists; only
  `PAPER_REFS` exists (line 56) and it *is* read, by `test_3_1` only.
- Every case builds a `_Case` with `expected_paper`, and `_run_case` prints
  `[x] Norm vs Paper 3D` for it. That line is informational and never asserted, and it is
  wrong by construction: it divides `|disp|` by `PAPER_REFS[...] * pressure` and compares
  the quotient against `expected_normalized` (a *normalized* value). Printed errors run
  from 527% to 6.2e8%. Nothing fails on it (section 5.5).

## 4. CCX parity group

### 4.1 `test_beam_shell_4cases_parity.py` (9)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_linear_static_with_analytical[tension_axial]` | free-face mean axial extension of a 4x20 MITC4 cantilever, point load 1000 N | analytical `F L/(E A)`, stated as a formula in the docstring | `TOL_ANALYTICAL = 0.02`, justified in a comment against measured errors | 0.560% | the test never skips: it needs no CalculiX |
| `...[compression_axial]` | same, reversed | analytical `-F L/(E A)` | 0.02 | 0.560% | – |
| `...[bending_fx]` | in-plane bending | analytical `F L^3/(3 E I_y)` | 0.02 | 0.890% | – |
| `...[transverse_fy]` | out-of-plane bending | analytical `F L^3/(3 E I_x)` | 0.02 | 0.682% | – |
| `test_linear_static_vs_ccx[4 cases]` | MITC4 vs CCX at the loaded node | **CCX 2.23, S4** (quadratic=False), same 4x20 mesh, point load at free-face centre | `rel_err <= 0.05` | 0.01% (tension), 0.01% (compression), 0.63% (bending_fx), 0.35% (transverse_fy) | the loaded node over-reads the axial cases by 5.09% vs analytical, for both codes; that is why the analytical test uses the face mean |
| `test_modal_first_five_modes` | first 5 modes from 12 requested, matched by Hungarian assignment | **CCX 2.23, S4** eigenfrequencies, same mesh | `tol=0.05` | max rel err 5.44e-03 (0.54%) | the matching step (not the physical modes) is what keeps this test robust; the printed table shows modes 1-12 on both sides |

### 4.2 `test_isotropic_shell_parity.py` (1)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_transverse_tip_displacement` | MITC4 vs CCX out-of-plane tip displacement, 2x10 mesh | **CCX 2.23, S4**; the module also defines an analytical Mindlin tip formula that the test never uses | `tol = 0.05` | AE 0.150095 m, CCX 0.145747 m, ratio 0.971, difference 2.9% | two defects: (a) the comment says "Allow 10% tolerance" while the code is 0.05 (`:365-366`); (b) the CCX side is `max abs(V)` over *all* FRD nodes (`:326-356`) while the AeroElast side is the loaded centre node, so the comparison is not like-for-like. The unused `analytical_tip_displacement` helper is dead reference code. |

### 4.3 `test_composite_beam_parity.py` (5)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestCompositeMaterial::test_laminate_abd_matrices` | ABD dict keys, lengths 9/9/9, thickness passthrough | none (schema test) | exact keys/`==` | not applicable | no numerical reference |
| `TestCompositeMaterial::test_mesh_connectivity` | node sets exist and sit at y=0 / y=L | mesh construction invariant | `atol=1e-12` | not printed | – |
| `test_composite_axial_tension` | [0/90/45/-45]s 8-ply laminate axial, 4x10 mesh | **CCX 2.23, S8R** + `*SHELL SECTION, COMPOSITE` (quadratic=True) | `rel_error < 0.1` | AE 46.49 um, CCX 48.09 um, 3.33% | the 10% window has no stated justification. The AeroElast block (mesh, material, assembly) is duplicated verbatim in the body. |
| `test_composite_isotropic_equiv` | same mesh with an isotropic equivalent mapped through a single-ply laminate | **CCX 2.23, S8R**, single isotropic ply (so only the element formulation differs) | `rel_error < 0.1` | AE 33.90 um, CCX 34.73 um, 2.37% | same unjustified 10% window |
| `test_composite_bending` | laminate transverse bending, 100 N at the free centre | **CCX 2.23, S8R** | `rel_error < 0.1` | AE 1671.42 um, CCX 1700.86 um, 1.73% | same unjustified 10% window |

### 4.4 `test_orthotropic_shell_parity.py` (3)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_orthotropic_axial` | single 0-ply orthotropic laminate, transverse Fy, 4x10 mesh, `elem_type=44` | **CCX 2.23, S8R** + composite section | `rel_error < 0.05` | AE 3 190 709.55 um, CCX 3 172 580.00 um, 0.57% | name says "axial" but the body and docstring apply a transverse load |
| `test_orthotropic_bending` | [0/90/90/0] in-plane lateral bending, load spread over the tip | analytical `F L^3/(3 EI_in)` with `EI_in = A11 B^3/12`, written out in the docstring | `rel_error < 0.05` | AE 1265.97 um, analytical 1221.54 um, 3.6% | the module docstring claims the lateral-bending tolerance was widened for a ~25-30% MITC4-vs-S8R gap, yet this test compares against the analytical value at 5% and passes at 3.6%; no CalculiX in this test |
| `test_multi_layer_iso_equivalence` | 4 laminar plies vs 1 layer of the same total thickness, MITC4 vs MITC4Composite | internal identity (A, B, D algebra), documented in the docstring | `rel_diff < 1e-4` | 5.437e-05 | a consistency identity rather than an external reference; the docstring states the algebra, so it is a legitimate duplicate-path check |

### 4.5 `test_shell_convergence.py` (2)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `test_in_plane_bending_convergence` | observed order of the MITC4 in-plane tip displacement over 4 meshes, by Richardson self-convergence, plus the extrapolated limit | analytical `P L^3/(3 E I)` with `I = t B^3/12`; the Timoshenko-vs-Euler-Bernoulli shear floor is derived in the comment (9.6 um on 1230 um ≈ 0.8%) | `MIN_ORDER = 1.5` and `EXTRAPOLATED_TOL = 0.02`, both justified in the module docstring | orders 1.7416 and 1.7638 (agree, delta = 0.022); Richardson limit 1238.1350 um vs 1230.7692 um = 0.5985% | the strongest tolerance justification in the suite; the raw pairwise orders (2.878, 0.153, -0.698) are printed and explicitly not asserted |
| `test_composite_laminate_gap_mesh_study` | whether the AeroElast-vs-CCX laminate gap shrinks (mesh artifact) or plateaus (formulation/ABD) | **CCX 2.23, S8R** across 4 meshes; no gap value asserted | none: the only assertions are `np.all(np.isfinite(...))` and `... > 0` | gaps -4.2527 / -1.7312 / -1.3546 / -1.7294%; verdict printed: `PLATEAUS -> formulation / ABD` | an analysis script, not a test: it cannot fail for any formulation. The docstring says "No gap value is asserted yet", so this is deliberate, but it should be read as a measurement, not as coverage |

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

The three independent CalculiX formulations bracket the analytical value and agree with
each other within 1.11%, so the CCX reference model for the strip is not the source of the
37.7% AeroElast gap on the new shell element. The module skips cleanly when CalculiX is
absent (`conftest.ccx_bin_or_skip`, path overridable with `CCX_BIN`).

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
| `TestLinearStaticCantilever::test_fx_in_plane` | axial membrane load, 8x4 mesh | analytical `P L/(E A)` | `error < 5.0` (percent) | FEM 2.824e-05, ana 2.857e-05 -> 1.2% | the comment says "Should match within 10%" (`:400`) while the assertion is 5% |
| `test_fy_in_plane` | in-plane bending, ny=8 | analytical `P L^3/(3 E I_z)`, `I_z = h b^3/12` | `error < 5.0` | FEM 1.134e-02 vs ana 1.143e-02 -> 0.8% (only raw values printed) | 5% unjustified |
| `test_fz_out_of_plane` | out-of-plane flexure | analytical Timoshenko `P L^3/(3 E I) + P L/(k G A)`, `k=5/6` | `error < 5.0` | FEM 1.121e+02 vs ana 1.143e+02 -> 1.92% | the failure message says "FX error" for the FZ case (`:491`) |
| `test_in_plane_ratio_constraint` | uY/uX ratio | beam theory `4 (L/b)^2 = 400` | `0.5 ratio_ref <= ratio <= 1.2 ratio_ref` | 399.84 vs 400.00 (0.04%) | the window is asymmetric (-50% / +20%) and unjustified; `test_shell_validation_fixed.py` asserts the same ratio at ±2% |
| `TestNonlinearStaticCantilever::test_large_displacement_tip_load` | "large displacement" response | none: the body asserts the linear estimate exceeds L, then that the solve raises `RuntimeError` | `abs(dz_lin) > L` | linear estimate 1.121e+02 | the name promises a numerical large-displacement validation; the body validates the divergence *error path* only |
| `TestModalAnalysis::test_first_mode_frequency` | first cantilever mode | analytical `(1.875104^2 / 2pi) sqrt(E I/(rho A L^4))`, `I = b h^3/12` | `error < 5.0` | 0.848 Hz vs 0.838 Hz -> 1.2% | 5% unjustified |
| `test_higher_modes` | "higher modes" | **none** | only strict monotonicity `f2 > f1`, `f3 > f2` | [0.848, 5.501, 16.542] Hz | the name promises validation; no reference value is compared. `simply_supported_plate_central_load`, `simply_supported_plate_uniform_pressure`, `nonlinear_large_displacement_cantilever` and `build_simply_supported_mesh` are all dead. |

### 5.3 `test_shell_validation_fixed.py` (7)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestLinearStatic::test_fx` | axial, 8x4 | analytical `600 L/(E b h)` | `TOL_STATIC = 3.0` (percent), justified in a comment with the measured errors | 1.15% | – |
| `test_fy` | in-plane bending | analytical `600 L^3/(3 E (h b^3/12))` | 3.0% | 1.19% | – |
| `test_fz` | out-of-plane bending | analytical `600 L^3/(3 E (b h^3/12))` | 3.0% | 1.92% | – |
| `test_ratio_physical` | uY/uX ratio | beam theory 400 (documented derivation in the docstring) | ±2%, justified against the measured 399.84 | 399.84 (0.04%) | the docstring records that the previous window was -50%/+20% |
| `test_axial_load_converges_to_the_analytical_solution` | monotone convergence, 4 meshes | analytical `P L/(E A)` | monotone decrease and finest `error < 0.01`; the exclusion of (2,1) is justified | 2.702%, 1.154%, 0.748%, 0.487% | – |
| `TestNonlinearStatic::test_geometric_nonlinearity` | "geometric nonlinearity" | none: asserts `dz_lin > L`, then `RuntimeError` matching "SNES diverged" | `dz_lin > L` | linear estimate 1.121e+02 | same shape as its `test_shell_comprehensive` twin; no numerical nonlinear reference |
| `TestModal::test_first_mode` | first cantilever mode | analytical `(1.875104^2 / 2pi) sqrt(E I/(rho A L^4))` | `error < 2.0` (percent) | 0.848 Hz vs 0.838 Hz -> 1.2% | the tolerance is the only one in this module without a justification comment, in a module whose other tolerances all carry one |

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

### 6.2 `test_material_suite.py` (44)

Reference model: Euler-Bernoulli with `EI` taken from the CLT `A` or `D` block, plus the
Reddy CLT B-coupling formula. No test prints its margin.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestIsotropicAnalytical` (4) | MITC4/MITC3 out-of-plane, MITC4 in-plane and axial | analytical Euler-Bernoulli with `I_out = B t^3/12`, `I_in = t B^3/12`, `P L/(E A)` | `TOL = 0.05` | not printed | 5% class tolerance, no per-case justification |
| `TestOrthotropicSinglePly` (6) | [0], [90] out-of-plane (D22), [0] and [90] axial (A22, A11), [45] ordering | analytical `P L^3/(3 D22 B)` / `P L/(A B)`, with the index choice explained per docstring | `TOL = 0.05`, [45] uses ordering plus `> 1.5x` E-B | not printed | `test_mitc4comp_ply45_out_of_plane` is a qualitative ordering check, not a numerical one; the docstring explains why E-B is invalid at 45 deg |
| `TestSymmetricLaminates` (7) | B=0 for [0/90/90/0] and quasi-iso, out-of-plane, in-plane, axial | analytical with A11/D22; quasi-iso uses `err < 0.20` for the documented ~16% D16/D26 inflation | `TOL = 0.05` (0.20 for quasi-iso) | not printed | the 0.20 is justified in a comment; `test_symmetric_b_is_zero` thresholds `max abs(B) < 1e-6` against A entries of order 1e7 — an absolute threshold with no scale argument |
| `TestAsymmetricLaminates` (6) | [0/90] B != 0, B11 < 0, MITC3Comp and MITC4Comp B-coupling vs CLT, bending-to-extension, symmetric control | Reddy CLT `w_tip = B11 F L^2/(2 A11 D11_eff b)` | `TOL = 0.10`, justified as shear-correction uncertainty | not printed; the docstring records the measured MITC3Comp `-2.627691e-03` vs CLT `-2.616014e-03` (0.45%) | the module records that the previous assertion was `abs(uy) > 1e-8`, sign-blind; the current test asserts the sign and 10% against CLT. `test_asymmetric_b_nonzero` uses `max abs(B) > 1.0` — an absolute threshold |
| `TestIsoEquivalence` (2) | 4 isotropic plies vs 1 layer | internal identity (K trace, tip displacement) | `rel < 0.01` | not printed | duplicate-path consistency check |
| `TestABDMatrices` (8) | D11 = A11 h^2/12; B exactly zero; B11 hand formula; A11 = Q11 h; D11 = Q11 h^3/12; Qbar 45/0/90 identities; Cs positive definite | closed-form CLT algebra | `1e-6` .. `1e-10` | not printed | `test_asymmetric_b11_formula` recomputes `B11_hand` from the same ply z-integrals as the implementation: **tautological** reference. `test_qbar_45_symmetry` uses `max(abs(x), 1.0)` as denominator (absolute below 1.0) and its comment promises "both positive at +45", which is never asserted. |
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

### 6.4 `test_rust_assembler.py` (18)

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

### 6.5 `test_rust_composite.py` (22)

All references are element-matrix invariants or ordering statements; no test prints.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestMITC3BatchSanity` (8), `TestMITC4BatchSanity` (8) | batch `ke`/`me` shapes (18x18, 24x24), symmetry, PSD, per-layup symmetry (0, 45, quasi-iso, glass/epoxy, 3D tilted) | invariants | `atol=1e-6` (K), `1e-10` (M), PSD `-1e-6 max(eig)` | not printed | no external reference; the docstring states "physical sanity checks only" |
| `TestBatchComposite` (2) | 3-element batches: symmetry + PSD per element | invariants | `atol=1e-6` | not printed | – |
| `TestCompositeSanity::test_stiffer_fiber_direction` | `ke_0[0,0] > ke_90[0,0]` | ordering | strict `>` | not printed | sign/order only, magnitude unconstrained |
| `test_thicker_laminate_stiffer` | `norm(ke_thick) > norm(ke_thin)` | ordering | strict `>` | not printed | ordering only |
| `test_assembler_composite_mitc3/mitc4_symmetry` | `PyMeshAssembler` composite K/M symmetry | invariants | `atol=1e-6` / `1e-10` | not printed | – |

### 6.6 `test_rust_modal.py` (10)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestMITC4ModalCantilever::test_frequencies_match` | Rust modal vs SLEPc, 6x4 clamped plate, 6 modes | `_python_modal_solve` (PETSc + SLEPc, `_aeroelast.petsc_modal_solve` vs `modal_solve_coo`) | `rtol=1e-4` | rel err 1.76e-13 .. 3.87e-16 (6 modes) | duplicate-implementation comparison |
| `TestMITC4ModalCantilever::test_mode_shapes_orthogonal` | "mode shapes are M-orthogonal" | none: asserts `norm(mode) > 1e-10` | cannot fail (an eigenvector of a GHEP is nonzero by construction) | not printed | name promises M-orthogonality; body checks a non-zero norm |
| `TestMITC3ModalCantilever::test_frequencies_match` | same for MITC3 | SLEPc | `rtol=1e-4` | max 9.00e-15 | – |
| `TestSimplySupportedPlate::test_frequencies_match_python` | pinned-edge plate, 8x8 | SLEPc | `rtol=1e-4` | max 9.81e-15 | – |
| `TestSimplySupportedPlate::test_analytical_convergence` | first mode of a pinned plate | Kirchhoff `f11 = (pi/L^2) sqrt(D/(rho h))`, `D = E h^3/(12(1-nu^2))` | `rel_err < 0.05`, "expect ~1-5% for 8x8" | 49.9985 Hz vs 49.3288 Hz = 1.36% | – |
| `TestCOOModalSolve::test_coo_matches_element_path` | `modal_solve_coo` vs the `PyMeshAssembler` path | cross-path identity | `rtol=1e-10` | not printed | – |
| `TestModalBenchmark::test_benchmark_mitc4` | 8x8 and 12x12 timing | SLEPc frequencies as the correctness guard | `rtol=1e-3` | 0.3x "speedup" printed for both sizes (Rust slower) | a timing test scored as a test; the print shows Rust losing by ~3x, which is the opposite of what the class name implies |
| `TestCompositeModal::test_composite_frequencies_match` | quasi-iso composite cantilever | SLEPc | `rtol=1e-4` | max 3.70e-13 | – |
| `TestCompositeModal::test_composite_stiffer_than_isotropic` | name says "stiffer" | none: asserts the two frequency vectors are *not* `allclose` | `rtol=0.01` inside `allclose` | first mode 0.69 Hz vs 8.68 Hz | the body asserts "different", not "stiffer" |
| `TestCompositeModalMITC3::test_composite_mitc3_frequencies_match` | same, triangles | SLEPc | `rtol=1e-4` | max 1.36e-14 | – |

### 6.7 `test_stress_stiffened_solver.py` (14)

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
| `test_b_coupling_produces_bending_under_axial_load` | `w_tip < -1e-6` under axial tension, MITC4 strip | the docstring quotes `w_tip = -2.161433e-05` (CLT) and the measured `-2.1660e-05` (0.21%) | `w_tip < -1e-6` | not printed | the quoted margin (0.21%) is not asserted: the reference is commented out, only the sign and a magnitude floor are checked. The CLT formula is re-derived in the module, so the docstring's reference is the same arithmetic. |
| `test_symmetric_laminate_no_bending_under_axial_load` | `abs(w_tip) < 1e-9` for [0/90/90/0] | symmetry (B = 0) | `1e-9` absolute | not printed | – |
| `test_b_coupling_sign` | `B11 < 0` and `w_tip < -1e-6` | as above | `1e-6` floor | not printed | duplicates the previous test's load case; the two share ~80 lines of identical body |

## 7. Rotor and FSI group

### 7.1 `test_rotor_inertial.py` (38)

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

### 7.2 `test_rotor_physical_consistency.py` (20)

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

### 7.3 `test_rotor_rust_parity.py` (49)

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

### 8.1 `test_bem_polars.py` (18)

No test prints.

| group | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestPolarData` (4) | interpolation at 0 and at 5deg, periodic wrapping, vectorisation | the fixture's own `cl = 2 pi sin(alpha)` definition | `atol=0.05` (Cl), `atol=1e-6` (wrap), `< 1e-3` (zero) | not printed | **tautological reference**: `cl` is generated as `2 pi sin(alpha)` and then compared to `2 pi sin(alpha)`. The test validates the interpolator round-trip, not aerodynamics. `atol=0.05` is not justified. |
| `TestAirfoilAero` (3) | exact-Re selection, nearest-Re selection, single-polar passthrough | library semantics | exact `==` on `polar.re` | not printed | – |
| `TestBladeAeroYAML` (11) | blade length, hub radius, rotor radius, n_blades, station ordering, chord bounds, twist bounds, polars present, alpha range, `std(cl)`, property shapes | the WindIO YAML that the code under test also loads | assorted absolute thresholds (`> 100`, `> 0`, `< 10`, `< 30deg`, `std > 0.1`) | not printed | the assertions are thresholds chosen by inspection, not from a source: they detect a broken parse, not a wrong value. `rotor_radius > hub_radius` and `alpha_range >= pi` are the meaningful ones. |

`tests/IEA-15-240-RWT.yaml` is the fixture both this module and `test_bem_engine.py` load;
it is a copy of the public IEA-15-240-RWT definition, not an independent reference.

### 8.2 `test_bem_engine.py` (module skipped) and `test_blade_mesh.py` (0 collected)

| test | what it validates | reference | tolerance | measured margin | notes |
| --- | --- | --- | --- | --- | --- |
| `TestBEMSolverParked::test_parked_alpha_close_to_twist` | parked blade AoA = 90deg - twist (geometric identity) | geometry | `residual < 15deg` for `> 80%` of stations, justified as "generous" | **not measured** (module skipped: `ccblade` absent) | thresholds with no derivation (`> 0.8` of stations within 15deg) |
| `TestBEMSolverParked::test_parked_thrust_positive`, `..._torque_near_zero` (`abs(power) < 1e6`), `..._Np_mostly_positive` (`> 0.7`) | sign/order sanity | none | loose thresholds | not measured | cannot detect a magnitude error |
| `TestBEMSolverRotating` (5) | rated power/thrust/torque positive, induction in range, Cl non-trivial | none | `frac_ok > 0.7`, `np.any(abs(cl) > 0.1)` | not measured | same class |
| `test_blade_mesh_generation[NOTSET]` | that meshing a reference turbine produces nodes, elements and root/outer-shell/shear-web node sets | none | `node_count > 0`, `elements_count > 0` | not measured (empty parameter set) | structural smoke test; the assertion cannot fail for a mesh with any content |

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

## 9. Flag summary

Grouped by defect class, with the file and the evidence.

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
| `test_ko2017_performance.py` | `test_3_2_circular_plate_tables_6_to_7` (SS rows) | validates against Table 6, not Table 7 |
| `test_ko2017_performance.py` | `test_3_5_twisted_beam_tables_12_to_13` | **RESOLVED** — the docstring and the code now agree: the four cases run with no xfail marker and `_TWISTED_BEAM_CASES` carries `xfail_reason=None` throughout, which is what this row used to flag as a mismatch |

### 9.4 Tolerances with no stated justification, tolerance/comment mismatches, dead conditionals

| file:line | item | evidence |
| --- | --- | --- |
| `test_isotropic_shell_parity.py:365-366` | comment "Allow 10% tolerance" vs `tol = 0.05` | quoted code |
| `test_shell_comprehensive.py:400` | comment "Should match within 10%" vs `assert error < 5.0` | quoted code |
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
| `test_ko2017_performance.py` | `[x] Norm vs Paper 3D: ...` via `wref_paper` / `expected_paper` | read only by `_run_case`'s print; the printed errors are 527% to 6.2e8%. The comparison divides by a 3D reference and compares against a *normalized* expectation, and only `test_3_1` sets the fields. It is a broken diagnostic, not a validation. |
| `test_shell_convergence.py` | `raw-error log-log fit`, pairwise orders, `gap trend`, `DECISIVE ANSWER` | deliberately not asserted (documented in both cases) |
| `test_shell_comprehensive.py:542` | `Beam-theory ratio: 400.00` | followed by a ±(-50%/+20%) assertion, so it does carry an assertion, just a very wide one |
| `test_rust_modal.py` | per-mode tables, benchmark speedups | the frequency tables back a real `rtol=1e-4` assertion; the speedup column backs nothing |

### 9.6 Configuration drift

- `tests/test_blade_mesh.py` reads only `examples/reference_turbines/yamls`, which does not
  exist, so its single test never runs; `tests/IEA-15-240-RWT.yaml` exists and is not in
  `conftest._BLADE_YAML_CANDIDATES` or in the module's list.
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

## 12. Reproducing this matrix

```bash
export PATH="$HOME/miniconda3/envs/aeroelast-dev/bin:$PATH"
cd /home/efirvida/Desktop/dev/fem-shell

# counts only (addopts= defeats the -v in pyproject.toml)
python -c "import pytest,sys; sys.exit(pytest.main(['tests','-o','addopts=','-q','--tb=no']))"

# printed margins, per file
python -c "import pytest,sys; sys.exit(pytest.main(['tests/test_ko2017_performance.py','-o','addopts=','-q','-s','--tb=line']))"
```

CCX must be on `PATH` (or `CCX_BIN` set) or the parity rows in 4.1-4.4 turn into skips
silently: `conftest.ccx_bin_or_skip()` calls `pytest.skip`, so a machine without CalculiX
reports a green suite with those tests absent.
