# Fidelity audit — `mitc4_plusd` against its own papers, with CCX as an external judge

**Status:** plan + checklist, opened at the end of the WU9e/WU9f session. Nothing in this
file authorizes a code change; each item is a *work unit* with its own measurement.

**Why this exists.** Four defects have been found so far in the MITC4+/D implementation,
and **every one was caught by an external oracle, never by internal reasoning**:

| # | defect | oracle that caught it | status |
| --- | --- | --- | --- |
| 1 | `e3` built from unit diagonal normals → a 2.4° frame error | Ko, Lee & Bathe (2017), C&S 182:404-418, Eq. (10) / the isotropy test | **fixed** (WU6d, `426249d`) |
| 2 | `W_22 = cm/9` **and** the `s2 = 4/h²` B-scaling → the `t`-rule entered twice, `16/h⁴` too stiff | the twisted-beam benchmark | **fixed** (WU9b, `5e90b4b`) |
| 3 | the nodal director `V_n^i` built **element-locally** instead of as a node property | the benchmark + Eq. (1) ("the director vector **at the node**"), p. 405 | **found, measured, fix pending** (WU9e/WU9f; option A) |
| 4 | in-plane bending of a flat strip **locks ~1.59×** (uY/uX = 251.55 vs beam theory 400.00; the hybrid is 0.04% off) | the repo's analytical beam reference; CCX 2.23 S4/S8/S8R; and **the paper's own words, p. 410** | **NOT A DEFECT — it is the paper's declared behaviour.** The element's flat plane-stress membrane *is* the displacement-based Q4 membrane (measured: assumed vs compatible `1.9e-17`; the paper says so verbatim), so it carries the classical in-plane shear locking. Not in section B (faithful), not in D (inert), not in E (inert) — the carrier is the membrane block itself, and it is faithful. The tests encode the hybrid's non-paper ingredients (SRI/ERC/bubble). **Do not fix by adding an ingredient the papers do not have.** |
| 5 | a **false provenance claim**: "the printed Eq. (21) omits the leading `e_rs^m\|bil` term" (checklist B1, the `b_membrane_covariant_2017` comment, `mitc4plus-2017-extract.md` note F2) | Appendix A, Eq. (A.7), read with vision | **falsified** — the term is printed; see the B-audit record |

So the hypothesis "there are implementation errors left" is not speculation: it has a
track record. A targeted hunt is the wrong response to it; a **point-by-point audit of the
code against the papers** is the right one. The goal is a line-by-line traceability
argument: every formulation line attributable to a paper equation, plus external numerical
validation. That is what makes the element defensible.

## The three oracles (and what each one can see)

| oracle | source | what it sees | what it cannot see |
| --- | --- | --- | --- |
| Tier 1 | the authors' own basic tests (isotropy, zero-energy, patch tests; both papers) | constant-stress/constant-curvature states, frame and node-order invariance, rigid-body modes | **locking** (non-constant fields), the drill's behaviour at thickness, anything the papers do not test |
| published cells | Ko, Lee, Lee & Bathe (2017), C&S 193:187-206, Tables 8-19 (the MITC4+ column, vision-read) | discretisation behaviour on the standard benchmarks, against 3-D/analytical reference solutions | configurations the benchmark does not include (e.g. a flat in-plane bending strip) |
| **CCX 2.23** | `~/miniconda3/envs/aeroelast-dev/bin/ccx`, on `PATH` when the `aeroelast-dev` env is active | an **independent numerical** answer for any case we can model with S4/SC6, including cases no paper and no analytical formula covers | it is not truth: it needs a converged mesh of its own, and its number is only as good as its model |

**CCX validated as a judge (measured this session, hybrid in production, `pytest tests/test_beam_shell_4cases_parity.py -s`, 9 passed in 7.15 s):**

| case | analytical | CCX 2.23 | AeroElast (hybrid) | AE vs CCX |
| --- | --- | --- | --- | --- |
| tension_axial | 1.904762e-5 | 2.001930e-5 (5.10%) | 2.004412e-5 (5.23%) | 0.12% |
| compression_axial | −1.904762e-5 | −2.001930e-5 (5.10%) | −2.004412e-5 (5.23%) | 0.12% |
| **bending_fx** (in-plane, membrane-dominated) | 3.047619e-2 | **3.039750e-2 (0.26%)** | 3.040036e-2 (0.25%) | **0.01%** |
| transverse_fy | 3.047619e0 | 3.016300e0 (1.03%) | 3.026857e0 (0.68%) | 0.35% |

So **CCX independently backs the analytical in-plane bending reference** (0.26%) and the hybrid (0.01%). Defect #4 is therefore confirmed as an implementation defect of the new element, not a test artefact: the element locks ~1.59× in in-plane bending. Its first audit target is section B (the 2017 assumed membrane). The loaded node of the axial cases over-reads by ~5.2% for **both** codes (why the analytical test uses the free-face mean — `docs/validation-matrix.md:149`). The same run also shows the modal gap (CCX 12 modes vs AeroElast 5) that `test_rust_modal` reports.

### CCX element types and the composite section (user requirement: S4 + S8 + S8R, and multi-material)

`src/aeroelast/core/mesh/io/writers.py` already carries most of this:

| piece | where | state |
| --- | --- | --- |
| element map | `_CCX_ELEMENT_MAP` (line ~379): `"quad" -> "S4"`, `"quad8" -> "S8R"` | **present**; **`S8` (full integration) is absent** |
| quadratic export | `write_ccx_mesh(..., quadratic=True)`: inserts midside nodes and writes `*SHELL SECTION, COMPOSITE` with per-ply data | **present** |
| linear export | `quadratic=False`: S4/S3 + `*SHELL SECTION, MATERIAL=` (equivalent orthotropic) | **present** |
| the hard CCX constraint | writer docstring lines 603-608: **"CalculiX requires S8R/S6 for `*SHELL SECTION, COMPOSITE`"** | **a constraint, not a choice**: any multi-material CCX comparison is *inherently* quadratic; S4 + composite does not exist in CCX |
| a worked S8R + COMPOSITE case | `tests/test_composite_beam_parity.py` (`_write_composite_ccx_inp`, "S8R + COMPOSITE section", lines 209/233-244) + `_parse_frd_disp` | **present** |
| the judge for the failing family | no S8/S8R reference exists for the in-plane bending strip, the plate, the in-plane lateral or the modal cases | **to add** |

**Why second order matters for the audit.** S8R/S8 (and S8R with `*SHELL SECTION, COMPOSITE`) converge better than S4, so a three-way comparison (S4 / S8 / S8R) separates **discretisation** from **formulation**: if all three CCX element types agree with the analytical reference while our 4-node element is 37.7% off, the error is in the formulation; if they disagree among themselves, the reference model (mesh, section, BCs) is the first thing to fix. Note the honest asymmetry: the repo's FEM stays 4-node (MITC4 / MITC4+/D), so the second-order CCX model is a *reference*, and its mesh must be converged (or stated) for the comparison to be a judge rather than a mesh artefact.

**To add (a separate unit, tests + docs only, no element change):** (1) an explicit shell element-type selector `S4 | S8 | S8R` in `write_ccx_mesh`, keeping the `quadratic` boolean backward-compatible; (2) S4/S8/S8R references for every case in the failing family; (3) multi-material CCX rows (composite section, hence S8R by CCX's own restriction) for the laminate cases against the new element; (4) one `docs/validation-matrix.md` row per (case, element type) naming the element type, the mesh, the reference and the tolerance.

CCX is not new infrastructure: `tests/conftest.py::ccx_bin_or_skip()` resolves the binary,
`tests/test_isotropic_shell_parity.py` writes the `.inp`, runs `ccx` and compares
(`write_ccx_mesh(quadratic=False)`, S4). `docs/validation-matrix.md` is the canonical
record of what has been compared. **Known gap:** `tests/test_orthotropic_shell_parity.py`
resolves CalculiX with its own `shutil.which("ccx")` instead of the conftest helper
(`docs/validation-matrix.md:520`).

## Audit rule (non-negotiable)

1. Read the equation **from the PDF with vision** (`pdftoppm -png -r 300` + the read tool),
   never `pdftotext`, and **never trust `docs/formulations/*-extract.md` as the paper** —
   both extracts contain at least one wrong note each, and both were corrected by vision
   re-reads in this change.
2. Compare the code function **term by term**, with units, against that read.
3. Record one verdict per item: `faithful` / `defective` / `paper-printed-form-differs`
   (with the proof) / `not-applicable-to-this-element`.
4. A `defective` verdict is only accepted with a **measurement** that shows it, and the fix
   must not add an ingredient the papers do not have, must not change a test or a tolerance.
5. Every closed item gets a `docs/validation-matrix.md` row (test, case, reference source,
   tolerance, measured value, skip behaviour).

## Checklist — code site × paper equation

`M17` = Ko, Lee & Bathe (2017), "A new MITC4+ shell element", C&S 182:404-418.
`M25` = Ko, Bathe & Zhang (2025), "Continuum mechanics-based shell elements with six
degrees of freedom at each node — the MITC4/D and MITC4+/D elements", C&S 308:107622.
`DB84` = Dvorkin & Bathe (1984), Engineering Computations 1:77-88.

### A. Geometry and kinematics

| # | code site (`mitc4_plusd.rs`) | paper | status |
| --- | --- | --- | --- |
| A1 | `compute_local_coordinate_system` | M17 Eq. (10) (the plane normal `n`) | **defective once** → fixed (2.4° frame, WU6d) |
| A2 | `compute_characteristic_vectors` (`x_r`, `x_s`, `x_d`, `n`, `m_r`, `m_s`) | M17 Eqs. (4c)/(9)/(11) | audited (WU0-WU3); re-check `m^r`/`m^s` dual-basis normalisation |
| A3 | `regularized_inverse_2x2`, `j_loc_at`, `covariant_to_local_mapping` | M17 Eqs. (17b)/(21) context; the extract flags the covariant→local metric as "defined from the dual basis, not printed" | **partially audited** — the design's own open item 4 |
| A4 | `compute_node_directors` (`V_n^i`, `V_1^i`, `V_2^i`) | M17 Eqs. (1)/(8a): "the director vector **at the node**" | **DEFECTIVE — element-local, not nodal** (WU9e/WU9f). Fix = ADR-4 option A. Fix pending |
| A5 | `compute_j3d_enriched` (`g_r`, `g_s`, `g_t`) | M17 Eqs. (4b)/(8a) | audited (WU2) |
| A6 | `interpolate_position`, `interpolate_displacement` | M17 Eqs. (1)/(3)/(8b) | audited (WU2), and Eq. (3) re-read at p. 405 in this session |
| A7 | `local_components`, `build_t24`, `transform_to_global` | not a paper item — a framing convention | audited (WU4b) |

### B. The 2017 assumed membrane — **primary suspect for the in-plane locking**

| # | code site | paper | status |
| --- | --- | --- | --- |
| B1 | `compute_membrane_coefficients_2017` (`c_r`, `c_s`, `d`, `a_A…a_E`) | M17 Eq. (24), p. 409 + the printed `a_A…a_E`, p. 410 (= Appendix A Eq. (A.6), p. 416) | **faithful** (audited 2026-09-24, vision). `c_r = x_d·m^r`, `c_s = x_d·m^s`, `d = c_r²+c_s²−1`, `a_A = c_r(c_r−1)/(2d)`, `a_B = c_r(c_r+1)/(2d)`, `a_C = c_s(c_s−1)/(2d)`, `a_D = c_s(c_s+1)/(2d)`, `a_E = 2c_rc_s/d` (positive). The claim of a printed Eq. (21) defect is **falsified** — see the B-audit record |
| B2 | `covariant_membrane_b_row` (the 5 tying rows at Fig. 4's A-E) | M17 Eqs. (7b)/(9)/(15)-(17), pp. 406-408 + Fig. 4, p. 406 | **faithful** (audited, vision + algebra). Tying points A(0,+1)/B(0,−1) sample `e_rr`; C(+1,0)/D(−1,0) sample `e_ss`; E(0,0) samples `e_rs`; `a_r = ¼ξ_i`, `a_s = ¼η_i`, `a_d = ¼ξ_iη_i` match Eq. (9) |
| B3 | `b_membrane_covariant_2017` (Eqs. 27a-c) | M17 Eqs. (27a-c), p. 409-410 | **faithful** (audited, vision). All 15 coefficients match the printed page term for term, and the printed chain (17)+(18)+(19)+(21)/(25)+(26) reproduces (27a-c) exactly |
| B4 | `b_membrane_2017` (covariant→local, `2 e_rs` doubling) | M17 Eq. (7a) context; the covariant→local metric itself is **not printed** (Appendix A / §2 context) | **faithful**: `covariant_to_local_mapping` verified to `8.9e-16` against the analytic `T = diag(1,1,2)·M⁻¹·diag(1,1,½)` with `M = [[J11², J21², 2J11J21],[…]]`; the `2e_rs` doubling before the mapping is correct for the engineering-shear triple |

#### B-audit record (2026-09-24) — section B is **clear**, and the "Eq. (21) defect" is a misreading

Read with vision, `pdftoppm -png -r 300`, from `A_new_MITC4+_shell_element.pdf` (C&S 182:404-418) for
pp. 405-410 and **Appendix A pp. 416-417**. Every equation below was cropped and read at full
resolution, never with `pdftotext`.

**1. The assumed membrane field (B1/B2/B3) is faithful, term by term.**

- Eq. (24), p. 409: `B_1 = c_r²/d`, `B_2 = c_s²/d`, `B_3 = 2c_rc_s/d`, `B_4 = −c_r/d`, `B_5 = −c_s/d`,
  with `c_r = x_d·m^r`, `c_s = x_d·m^s` and `d = c_r²+c_s²−1`. The printed second form of `d`
  `= (x_e²·m^s)(x_e⁴·m^s) + (x_e³·m^r)(x_e¹·m^r) + 1` is **not** a sign error: with Eq. (13)'s edge
  vectors and Eq. (11)'s dual-basis identities it expands to `(c_s²−1) + (c_r²−1) + 1 = c_r²+c_s²−1`.
- The printed `a_A…a_E`, p. 410 (`a_A = c_r(c_r−1)/(2d)`, `a_B = c_r(c_r+1)/(2d)`, `a_C = c_s(c_s−1)/(2d)`,
  `a_D = c_s(c_s+1)/(2d)`, `a_E = 2c_rc_s/d`) are reproduced in `compute_membrane_coefficients_2017`
  with the array order `[a_A, a_B, a_C, a_D, a_E] = a_coeffs`. Appendix A Eq. (A.6), p. 416, repeats
  the five constants identically.
- Eq. (27a-c), pp. 409-410, match `b_membrane_covariant_2017`'s three rows coefficient for coefficient:
  row 0 `½(1−2a_A+s+2a_As²)`, `½(1−2a_B−s+2a_Bs²)`, `a_C(−1+s²)`, `a_D(−1+s²)`, `a_E(−1+s²)`;
  row 1 `a_A(−1+r²)`, `a_B(−1+r²)`, `½(1−2a_C+r+2a_Cr²)`, `½(1−2a_D−r+2a_Dr²)`, `a_E(−1+r²)`;
  row 2 `¼(r+4a_Ars)`, `¼(−r+4a_Brs)`, `¼(s+4a_Crs)`, `¼(−s+4a_Drs)`, `(1+a_Ers)`.
- Eq. (15), p. 407 (`e_rr^m = e_rr^m|con + e_rr^m|lin·s + e_rs^m|bil·s²`, and the two siblings) and
  Eq. (16) are exactly what `covariant_membrane_b_row` expands to, because
  `g_r·(u_r + s u_d) = (x_r+s x_d)·(u_r+s u_d)`. The tying points are Fig. 4, p. 406, as the code
  claims, and Eq. (17) assigns the sampled components (e_rr at A/B, e_ss at C/D, e_rs at E).

**2. The claim that the printed Eq. (21) omits the leading `e_rs^m|bil` term is FALSE.**

The claim appears in three places (checklist B1, the `b_membrane_covariant_2017` code comment,
`docs/formulations/mitc4plus-2017-extract.md` Note F2) with the stated proof "Eq. (27c)'s
`(1 + a_E·rs)` coefficient". The vision read contradicts it on its own:

- Eq. (21), p. 409, printed, reads
  `ẽ_rs^m|bil = B_1·(e_rr^m|con + e_rs^m|bil) + B_2·(e_ss^m|con + e_rs^m|bil) + B_3·e_rs^m|con + B_4·e_rr^m|lin + B_5·e_ss^m|lin`.
  The `e_rs^m|bil` term **is** there, inside the `B_1` and `B_2` parentheses, with coefficient `B_1+B_2 = (c_r²+c_s²)/d = 1 + 1/d` — the `1` the note (and Eq. (27c)'s `1`) refers to.
- **Decisive**: Appendix A, p. 416, performs the substitution the note claims was missing. Eq. (A.3)
  plus the condition Eq. (A.2) (`ẽ_rs^m|bil = 0` for the `Eq. (A.1)` patch-test mode) are substituted
  into Eq. (21) to give Eq. (A.4)/(A.5), hence Eq. (A.6), and the paper's own **Eq. (A.7)** is
  `ẽ_rs^m|bil = (c_r²/d)(e_rr|con + e_rs|bil) + (c_s²/d)(e_ss|con + e_rs|bil) + (2c_rc_s/d)e_rs|con − (c_r/d)e_rr|lin − (c_s/d)e_ss|lin` —
  i.e. Eq. (25) **with the `e_rs|bil` terms retained**. Eq. (21) and Eq. (25) are therefore the same
  equation with `B_1..B_5` substituted, and neither is defective.
- The `(1 + a_E·rs)` coefficient of Eq. (27c) multiplies `e_rs^m(E) = e_rs^m|con`, not `e_rs^m|bil`.
  The bare `1` comes from Eq. (19)'s added linear terms (`ẽ_rs = ē_rs + ½e_rr|lin·r + ½e_ss|lin·s`)
  with `ē_rs = e_rs^m|con`, and the `a_E·rs` from Eq. (26)'s `ẽ_rs^m|bil·rs`. It proves nothing about Eq. (21).
- Note F2's supporting argument ("for a flat rectangle `c_r=c_s=0` makes `B_1..B_5` zero, so Eq. (21)
  would give `ẽ_rs^m|bil = 0`, contradicting Eq. (22)") is void: on a flat **rectangle** `x_d = 0`, so
  `e_rs^m|bil = x_d·u_d = 0` and Eq. (22) requires `ẽ_rs^m|bil = 0` too — no contradiction. On a flat
  **distorted** element the coefficients are non-zero and Eq. (A.7) supplies the correct value.

**3. Independent re-derivation.** With `p,q,w,u,v,z` the five strain coefficients and `e_rs^m|bil`,
Eqs. (17)+(18)+(19)+(21)/(25)+(26) give `z̃ = a_A A + a_B B + a_C C + a_D D + a_E E` and then
`ẽ_rr^m = ½(1+s)A + ½(1−s)B + (s²−1)z̃`, reproducing (27a); `ẽ_ss^m` reproduces (27b); and
`ẽ_rs^m = E + (r/4)(A−B) + (s/4)(C−D) + rs·z̃` reproduces (27c) exactly. So the printed chain is
internally consistent and (27a-c) is its closed form — which is what the code implements.

**4. Consequence for defect #4.** Section B is **cleared**: the 2017 assumed membrane is a faithful
implementation of the printed formulation (and of Appendix A), and the in-plane bending lock is
**not** there. The audit moves on down the checklist: A3's metric is already verified above (`8.9e-16`);
D (assumed transverse shear) and the constitutive/integration path (F) are the next candidates, plus
the CCX judge of work unit 2/5. The record above also removes a documented "paper defect" that never
existed — a traceability defect in its own right, since the code comment and the extract asserted it.

### C. The 2017 bending

| # | code site | paper | status |
| --- | --- | --- | --- |
| C1 | `b_bending_covariant_2017` (`dx_b·du_m` of Eq. (8a), `dx_m·du_b`, `dx_b·du_b`) | M17 Eqs. (7c)/(7d)/(8a), p. 406 | audited (WU3) + the `∂x_b` refinement refuted by measurement in an earlier session |
| C2 | `b_bending_2017` | M17 Eq. (7a) `t`-linear measure | audited (WU8, retargeted normalisation) |

### D. The assumed transverse shear

| # | code site | paper | status |
| --- | --- | --- | --- |
| D1 | `compute_shear_tie` (tying A, B, C, D) | M17 Eq. (4)/(5), p. 405 | **faithful** (audited 2026-09-24, vision + algebra). Row 0 = `½(g_r·u_b + g_t·u_r)`, row 1 = `½(g_s·u_b + g_t·u_s)`, i.e. Eq. (4)'s covariant `e_rt`/`e_st`; `g_t = ½Σ a_i h_i V_n^i` and `u_b = ½Σ a_i h_i (θ_i × V_n^i)` (Eq. (3)/(8b)); `g_r·(θ×V_n) = θ·(V_n×g_r)` gives the code's `¼a_i h_i (V_n^i×g_r)`. Called **only** at the four tying points |
| D2 | `shear_covariant_to_local`, `b_shear_mitc4` | M17 p. 405 assumed field + DB84 Eq. (3), p. 78 | **faithful**. `b_shear_mitc4`: `ẽ_rt = ½(1+s)e_rt^A + ½(1−s)e_rt^B`, `ẽ_st = ½(1+r)e_st^C + ½(1−r)e_st^D` exactly as M17 p. 405 prints. **Note on the labels**: DB84 Eq. (3) writes `ẽ_13 = ½(1+r_2)e_13^A + ½(1−r_2)e_13^C` and `ẽ_23 = ½(1+r_1)e_23^D + ½(1−r_1)e_23^B` (its points A top, B left, C bottom, D right); M17 **relabels** them to A(top)/B(bottom)/C(right)/D(left). The code follows **M17's** labels, and the precompute supplies exactly A(0,+1), B(0,−1), C(+1,0), D(−1,0) — consistent. The local-frame map `shear_covariant_to_local` is **printed by neither paper** (M17 prints the covariant Eq. (4) and defers the tying construction to DB84; DB84 p. 78 keeps the covariant form and imposes the constraint via Lagrange multipliers, printing no metric map) — it is the design's 3D dual-basis definition, verified below |

### D-audit record (2026-09-24) — section D is **faithful**, and it is **provably inert** on the flat in-plane case

Vision reads: M17 p. 405 (Eqs. (4)-(5) and the assumed `ẽ_rt`/`ẽ_st` with Fig. 2) and
Dvorkin & Bathe (1984), *Eng. Comput.* 1:77-88, p. 78, Eq. (3) with Figs. 2-3
(`.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`).

**1. The sampled covariant shear (D1).** Eq. (4) is `e_ij = ½(g_i·u_j + g_j·u_i)` with
`g_i = ∂x/∂r_i`, `u_i = ∂u/∂r_i` (Eq. (5)). `compute_shear_tie` builds, for each node, exactly
`b[(0,6i+k)] = ½ g_t[k]·∂h_i/∂r` (from `u_r = Σ(∂h_i/∂r)u_i`) and
`b[(0,6i+3+k)] = ¼ a_i h_i (V_n^i × g_r)[k]` (from `u_b = ½Σ a_j h_j (θ_j×V_n^j)` and
`g_r·(θ×V_n) = θ·(V_n×g_r)`), i.e. row 0 = `½(g_r·u_b + g_t·u_r) = e_rt`; row 1 the `e_st`
sibling. Called only at A(0,+1), B(0,−1), C(+1,0), D(−1,0) (`grep` confirms the only call sites).

**2. The assumed field and the tying labels (D2).** `b_shear_mitc4` reproduces the printed
M17 interpolation term for term. The one subtlety found: **M17 and DB84 use different letters
for the same four points** — DB84 Eq. (3) has `(A,C)` on the `r_2` pair and `(D,B)` on the `r_1`
pair, M17 renames them `A(top)/B(bottom)/C(right)/D(left)`. The code follows M17 and the
precompute supplies the M17 order, so the pairing is right; a reader checking against DB84
alone would see a mismatch. Recorded.

**3. The local-frame metric is not in the papers; it checks out numerically.** The design's
`γ_a3 = 2 e_i3 (g^i·e_a)(g^t·e_3)` was verified against the analytic scaling on the flat
aligned element: it reduces **exactly** to `γ_13 = 8/(L h)·e_rt`, `γ_23 = 8/(b h)·e_st`
(computed: `T = diag(8000, 80000)` for `L=1, b=0.1, h=0.001`), and the dual-basis terms it
drops (`g^r·e3`, `g^s·e3`, `g^t·e1`, `g^t·e2`) are `~1e-19` for a 0.1% warp and change `T` by
`< 0.006%` for a 1% warp. So the approximation is immaterial at shell warps, and exact on the
flat elements of the failing strip. (D2's WU9b check — equality with the hybrid's `b_gamma_mitc4`
— is an internal cross-check and is **not** offered as evidence here.)

**4. THE FINDING: the transverse-shear block cannot cause defect #4.** On a flat element with
`V_n^i = e3` and the local frame aligned with the edges, `g_t = ½h e3`, `g_r = (L/2)e1`,
`g_s = (b/2)e2`, so in `compute_shear_tie`:
- `g_t[k] ≠ 0` only for `k = 2` → the translational rows touch **only `u_z`**;
- `V_n × g_r = e3 × e1 = e2` → `k = 1` only → **`θ_y`**;
- `V_n × g_s = e3 × e2 = −e1` → `k = 0` only → **`θ_x`**.

So the shear block has non-zero entries only on the set `{u_z, θ_x, θ_y}`. The exact symmetric
in-plane solution of the failing strip has `u_z = θ_x = θ_y = 0` (the load is in the plane and
the geometry is symmetric through the mid-surface), so the shear block contributes **zero
energy**, and the same argument kills the bending block (`∂x_b/∂r = 0` on a flat element ⇒
`e_b1` is rotation-only and `e_b2 ≡ 0`). **Therefore the flat in-plane stiffness is the sum of
exactly two blocks: the assumed membrane (B, cleared) and the 2025 drill (E).**

**5. New, sharper signature of defect #4.** WU9f's failure list for this family is `test_fy`,
`test_ratio_physical`, `test_fy_in_plane` — **`test_fx` (axial) is not on it**. So the element
is *correct in uniform/axial in-plane* and *37% too stiff in in-plane bending*: the defect is in
the **gradient** response of a block that is inert for a constant strain. With B cleared, D
inert, A3 verified and F already audited, the load falls on the **2025 drill block (E)** — the
only 2025 ingredient, the only remaining active block, and a gradient-only operator. The
decisive measurement is a block-isolation harness on the flat strip (drill on/off; assumed vs
displacement-based membrane), the same instrument method as WU9e.

### Block-isolation record (WU9i, 2026-09-24) — defect #4 is the **paper's declared behaviour**, not an implementation defect

Instrument: the `#[ignore]`d in-crate test `wu9i_block_isolation_flat_inplane_strip`
(`crates/aeroelast-core/src/elements/mitc4_plusd.rs`, test module only; reproducible with
`cd crates && cargo test -p aeroelast-core wu9i -- --ignored --nocapture`). Dense 270-DOF solve of
the exact failing case (8x4 strip, `L=1, b=0.1, h=0.001`, `E=2.1e11`, `nu=0.3`, root clamped, 600 N
in the measured direction on the free edge).

**Validation.** FULL reproduces the WU9f flip exactly — `ratio = 251.5515` vs `251.55` (0.00%),
`ux = 2.829989765e-5` vs the analytical `2.857142857e-5` (0.95%). The instrument is trustworthy.

**Block energy shares on the FULL solution** (membrane / bending / shear / drill):
`1.000000000000 / 0 / 0 / 0` under +x (axial) and `1.000000000002 / 0 / 0 / 0` under +y (in-plane
bending). **The flat in-plane stiffness is the membrane block, 100%, and nothing else.**

**The refutation (my hypothesis was wrong).** I predicted the 2025 drill. Measurement says the
drill is **inert**: `NO_DRILL` is singular only because the 40 `theta_z` DOFs lose their stiffness
(40 null modes at `1e-12 lambda_max`), and pinning those 40 to zero returns the FULL numbers to 12
digits; a membrane-only system (80 DOF) also returns `251.5515`. Shear and bending carry zero share
(as the D record proves algebraically). Confirmed by the worker and re-run by the parent.

**The mechanism.** On a flat element `x_d = 0`, so `a_A..a_E` all vanish and Eqs. (27a-c) collapse:
`e~_rr = ½(1+s)e_rr(A) + ½(1-s)e_rr(B)`, etc. — which is **exactly** the compatible
displacement-based field of Eqs. (15)/(16), because that field is linear in `s`/`r` on flat geometry.
Element-level `max|K_assumed - K_compat|/max|K_assumed| = 1.9e-17`. The 2017 assumption therefore buys
nothing on flat, undistorted geometry: its relief lives in the `a_E·rs` term, i.e. in `e_rs^m|bil =
x_d·u_d`, which is zero when `x_d = 0`.

**And the paper says so, verbatim** (Ko, Lee & Bathe (2017), C&S 182:404-418, p. 410, read with vision):

> "We note that the membrane part of the new MITC4+ shell element is identical to that of the
displacement-based element when the element geometry is flat. That is, in two-dimensional plane
stress problems, both shell elements always yield the identical solutions."

So the element is **faithful**, and its ~37% in-plane bending excess on the repo's flat regular mesh
is the published MITC4+'s behaviour: the classical 4-node plane-stress shear locking of the
compatible membrane. The paper's own numerical section never tests a flat in-plane bending strip
(its flat case is a plate in out-of-plane bending), so this behaviour is declared but not
benchmarked by the authors.

**Consequence for the S3 flip gate — a maintainer decision, and a STOP.**
`TestLinearStatic::test_fy`, `test_ratio_physical` (`uY/uX = 400 ± 2%`) and
`TestLinearStaticCantilever::test_fy_in_plane` encode the **hybrid's** in-plane-bending accuracy.
The hybrid reaches it with ingredients no paper has (SRI, the Winkler-Plakomytis ERC, the 2-DOF
rotation bubble, the `beta_w` warping penalty, `k_drill`). A faithful MITC4+/D cannot pass those
three tests, and passing them would require adding exactly those ingredients — which the audit's
rule 4 forbids. Per the standing instruction ("if the only way to pass is an ingredient the paper
does not have: STOP and report"), this unit **stops and reports** instead of adding one.

The faithful options are maintainer decisions, none of them a formulation change:
1. re-scope/retire the three flat in-plane bending tests with the p. 410 quote as the recorded
   rationale, and record the paper's own declared behaviour as the new expectation;
2. keep them as a documented, non-blocking deviation for the new element (they are satisfied by the
   hybrid, which is not paper-faithful);
3. restrict the fidelity claim to the papers' own validated envelope (`t/L <= 1/100`, curved/
   warped and out-of-plane cases), where the flip's remaining failures are only the flat family and
   the geometrically nonlinear path (the design's recorded open item 5).

### CORRECTION to the WU9i record (2026-09-24, same session) — the "NOT A DEFECT" verdict above is **WITHDRAWN**

The WU9i verdict compared the element's flat behaviour with the **2017** paper's p. 410 statement
("the membrane part ... is identical to that of the displacement-based element when the element
geometry is flat"). That statement is about the **2017 MITC4+ (five DOF, no drill)**. This element is
**MITC4+/D**, and the **2025** paper adds the drill *to the membrane strain itself*:

> Ko, Bathe & Zhang (2025), C&S 308:107622, **Eq. (22a), p. 10** (read with vision):
> `e_ij = e_ij^m + e_ij^md + t e_ij^b1 + t^2 e_ij^b2` with `i,j = 1,2`. "This element is denoted as 'MITC4/D'."

It is a **sum of strains**, so the strain energy carries the cross term `2 (e^m)^T C e^md` (and
`2 (e^md)^T C e^b2` at the `t^2` level). `compute_ke_local_with_drill` instead adds the drill as an
**independent block** (`k += drill_ke_local(pre)`, `sum_g B_md^T cm B_md`), so **no cross term exists**.
With a free `theta_z` and a positive-semidefinite independent block, the minimiser drives the drill
strain to zero and the block contributes nothing to the displacement response — which is exactly
what WU9i measured (drill share 0; identical numbers with the drill off and `theta_z` pinned).
**The measurement was the symptom of the missing coupling, not a clearance of the drill.**

So: section B is faithful, section D is faithful and inert, and the flat in-plane response equals the
plain compatible Q4 membrane — because the 2025 coupling that would enrich it is **absent from the
code**. `defect #4` is therefore, with high confidence, a **fidelity defect with a paper-faithful fix
(build the membrane row as `B_m + B_md`), not the paper's declared behaviour, and not a new
ingredient**. It is not yet measured: the size of the effect is unknown until the fix is implemented.

Two further purity items confirmed/flagged in the same read:
1. **The 2025 paper restates the assumed membrane** — Eq. (24a) `e~_ij^m = a_0|ij e_rr^m(A) + ...`,
   Eq. (24b) (vision-confirmed) `e~_rr^m = e_rr^m(0,0) + (sqrt(3)/2) lambda(r,s) (e_rr^m(A) - e_rr^m(B)) s`,
   `e~_ss^m` likewise in `r`, `e~_rs^m = e_rs^m(0,0)`, with `lambda(r,s) = j0/j(r,s,0)`. Our code
   implements the **2017** Eqs. (27a-c) instead. On flat geometry they coincide to the leading term
   but the slope coefficient differs (`sqrt(3) lambda` vs `1`), so this needs its own audit before
   claiming "pure 2025". (Eq. (24c)/(25a) still to read.)
2. **The validated envelope** of the 2025 study is `t/L = 1/100, 1/1000, 1/10000` (Fig. 15 legend,
   vision-confirmed, p. 11), and its own 2D planar-beam examples are plane-stress cases on regular
   and distorted meshes (§3.2, Tables 1-3) — **an oracle we have not yet used**.

### E. The 2025 drill

| # | code site | paper | status |
| --- | --- | --- | --- |
| E1 | `compute_drill_edges` (`c_r`, `c_s`, `x_m^l`) | M25 Eqs. (13b)/(13c)/(19c), pp. 6/10 | audited: Eq. (18)'s printed form (no `1/‖x_m^l‖`, WU9d) |
| E2 | `drill_midside_shape_derivatives` | M25 Eqs. (11a)/(11b), p. 5 | audited (byte-identical to the printed rules, delib. zeros) |
| E3 | `drill_jacobian_ratio` (`j0/j`) | M25 Eq. (17b), p. 8 | audited: numerically inert on the benchmark (WU9d) |
| E4 | `b_drill_membrane_2025` | M25 Eq. (18), p. 8 | audited term by term to `1e-12` (Tier-1 identity test) |
| E5 | the edge pairing / `θ_{i+1}^D − θ_i^D` telescoping | M25 Eqs. (12d)/(16a)/(19b), pp. 6/7/10 | audited (WU9d) |

### F. Through-thickness integration and the constitutive

| # | code site | paper | status |
| --- | --- | --- | --- |
| F1 | `resultant_moment_matrix` (9×9 `W`) | M17 Eq. (7a), p. 406 + the 2×2×2 rule p. 410 | **defective once** → fixed (`W_22`, WU9b) |
| F2 | the integration rules (2×2 surface, the `2×2 t` rule) | M17 p. 410; M25 §3.1 p. 13 | audited |
| F3 | `ShellConstitutive::transverse_shear_uncorrected` (ADR-1) | M17 p. 410 (no numerical factor) | audited (WU4a) |
| F4 | the `e^md` `t⁰` slot of Eq. (26) | M25 Eq. (26), p. 4: `e_ij = ẽ^m + e^md + t·e^b1 + t²·e^b2`, `i,j = 1,2` | **DEVIATION IDENTIFIED, not fixed**: the code adds `∫B_mdᵀ cm B_md` as a separate block, so the cross terms of `e^md` with the `t²` block (`W_02`) are absent. Measured: folding changes the N=8 cells by ≤3% and fixes neither thick case (WU9d) — but the paper's form is the folded one. **Tracked, low priority** |

### G. Assembly, boundary conditions and the surrounding kernels

| # | code site | paper | status |
| --- | --- | --- | --- |
| G1 | `compute_ke_local_with_drill`, `compute_ke_global` | assembly of the above | audited (WU4/WU5) |
| G2 | the BC/load conventions used by the benchmarks | the benchmark's own definitions | audited for the twisted beam; the flat in-plane armature of the repo's tests is now **audited against CCX**: `tests/test_ccx_shell_element_types_parity.py` reproduces `test_shell_validation_fixed.py`'s exact strip (8x4, clamped root, 600 N in +y) with CCX S4/S8/S8R and all three match the analytical value (0.61/0.23/0.49% at 8x4) — so those BCs/loads are the physical ones |
| G3 | `compute_me_*`, `compute_body_load_global`, `compute_k_sigma_global`, `compute_element_stress`, `compute_centrifugal_prestress` | mostly repo-level, not paper items | not audited (not part of the papers' claims) |

### H. The geometrically nonlinear path — **CORRECTED 2026-09-24: it is NOT out of scope**

The original text here declared the nonlinear path "out of the papers' scope". **That was wrong and is
withdrawn.** `.sources/papers/` holds the papers that specify it, and the primary one is:

> **Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear analysis", Computers and
> Structures 185:1-14** (`The_MITC4+_shell_element_in_geometric_nonlinear_analysis.pdf`), plus Dvorkin &
> Bathe (1984) for the total-Lagrangian origin and `mitc3+_no_lineal.pdf` for the triangular sibling.

The paper specifies the large-displacement and large-rotation formulation of **this same MITC4+ element**
using the **total Lagrangian formulation with Green-Lagrange strains** and, per its abstract, the
"assumed shear and membrane fields of the element for the total Lagrangian formulation".

**First-pass verdict (vision-read paper equations vs the code, recorded in full in
`odd/tasks/mitc4plusd-2025-purity.md`, §"N-audit, first pass"):** the implemented path is **not** the
paper's formulation. It is a bounded approximation — the code's own docstrings call it "BOUNDED NONLINEAR
PATH (design open item 5, risk 9)" — and the audit found the specific gaps:

| # | code site | paper | status |
| --- | --- | --- | --- |
| N1 | `update_normals_with_displacements`; the initial-geometry `pre` used by `compute_kt_global`/`compute_fint_global` | Eqs. (1)/(4c)/(7): the CURRENT director `^t V_n^i` and CURRENT covariant base vectors `^t g_i` | **defective — the director update has zero call sites, so the current-configuration geometry is never formed** |
| N2 | `membrane_strain_nl`, `compute_b_nl` | Eq. (9)'s `_0 eta_ij = 1/2(u_{1,i}.u_{1,j} + ^t g_i.u_{2,j} + u_{2,i}.^t g_j`) | **defective — only the first of the three terms is implemented; the quadratic director increment `u_b2 = -1/4 sum a_i h_i(alpha_i^2+beta_i^2)^t V_n^i` (Eq. (5c)) has no representation** |
| N3/N4/N5 | `b_membrane_2017`, `b_shear_mitc4` | Eqs. (10)/(11a) and the §2.3 assumed membrane, applied to the **GL** strains | **not yet audited** |
| N6/N7 | `geometric_stiffness_local`, `compute_kt_global`, `compute_fint_global` | the paper's `^t_0 K` and `^t_0 f_int` | **approximation: `K_0 + K_L + K_sigma` with `K_0` from the initial geometry** |
| N8 | the corotational block (`compute_membrane_strain_log`, `update_corotational_frame`, `frame_incremental_rotation`, `polar_decomposition`, `log_strain_from_polar`, `quaternion_*`) | no paper — the formulation is **total Lagrangian**, not corotational | **dead (zero call sites) and UNSOURCED → delete (T5)** |

The 5 nonlinear test failures at the gate are the symptom. The modal path (`test_rust_modal`) remains
partly hybrid-calibrated and is likewise not a paper item.

## CCX judge record (2026-09-24) — the in-plane bending reference is confirmed by three independent CCX element types

`tests/test_ccx_shell_element_types_parity.py` (new) runs the **exact** case of
`tests/test_shell_validation_fixed.py` (`L x b x h = 1.0 x 0.1 x 0.001` m, `E = 2.1e11`,
`nu = 0.3`, regular `8x4` quad mesh, every node at `x = 0` clamped in 6 DOF, total 600 N in
`+y` on the free edge `x = L`, measuring `uy` at `(L, b/2)`) through CCX 2.23 with three
element formulations, using the new `write_ccx_mesh(shell_element_type=...)` selector and a
consistent edge traction of resultant 600 N:

| element | 8x4 `uy` [m] | vs analytical `F L^3/(3 E I) = 1.142857e-2` |
| --- | --- | --- |
| CCX 2.23 **S4** (linear, Dvorkin-Bathe shear) | `1.135840e-2` | 0.614% |
| CCX 2.23 **S8** (quadratic, full integration) | `1.145430e-2` | 0.225% |
| CCX 2.23 **S8R** (quadratic, reduced integration) | `1.148490e-2` | 0.493% |

Mutual spread `(max−min)/min = 1.11%`, and the mesh study `4x2 / 8x4 / 16x8` gives spreads
`2.67% / 1.11% / 0.42%` — the family difference shrinks under refinement, so the reference
model is converged, not a mesh artifact. Note the physical floor: **S4/S8R agree in bending,
but S8R is the element CCX requires for `*SHELL SECTION, COMPOSITE`**, so the composite rows
(open) are inherently S8R.

**Consequence.** The analytical beam reference for the flat in-plane bending strip is
third-party confirmed (three element types, one of them full integration). The new shell
element's `uY/uX = 251.55` (37.1% below the beam-theory 400) is therefore a **formulation
defect of our element, not a reference-model artifact** — and it is not in section B. The
next audit targets are D (assumed transverse shear), F (through-thickness/constitutive) and
the flat-plane plate/modal family of the flip's failures.

**Files.** `src/aeroelast/core/mesh/io/writers.py` (`shell_element_type` selector, S8 added
to the quadratic builder, composite guard), `tests/test_ccx_shell_element_types_parity.py`
(new, 4 tests, skips without CCX), `docs/validation-matrix.md` §4.6. `pytest -m "not slow"`
moves from `344 passed / 3 failed / 2 skipped` to `348 passed / 3 failed / 2 skipped` (the
same 3 pre-existing failures). No element line changed, no tolerance widened.

## First work units (in order)

1. **CCX tie-breaker on the in-plane bending strip.** The new element passes the existing
   CCX-parity tests (`test_beam_shell_4cases_parity.py::test_linear_static_vs_ccx`,
   0.01%-0.63%) yet fails the repo's beam-theory reference by **37.7%** on
   `test_fy` / `test_ratio_physical`. **Steps 1 and 2 are DONE** (2026-09-24): CCX 2.23 is
   installed and runs; and the *failing test's own* geometry (`L = 1.0, b = 0.1, h = 0.001`,
   the `nx = 8, ny = 4` mesh, clamped root, the 600 N in-plane edge load) now has a
   third-party number: `tests/test_ccx_shell_element_types_parity.py` runs **S4, S8 and
   S8R** on that exact strip and all three agree with the analytical `F L^3/(3 E I)` value
   (S4 0.614%, S8 0.225%, S8R 0.493% at 8x4; spreads 2.67% → 1.11% → 0.42% over 4x2/8x4/16x8,
   i.e. converged). So the reference model is **not** the cause and the 37.7% is ours. Rows
   recorded in `docs/validation-matrix.md` §4.6. (Element-type selector: `write_ccx_mesh`
   gained `shell_element_type=None|"S4"|"S8"|"S8R"`, backward compatible with `quadratic`.)
2. **Audit B1-B3 against pp. 409-410 with vision** (the membrane tying rows and the
   `a_A…a_E` coefficients). This is the prime suspect for defect #4 if CCX agrees with beam
   theory. Check in particular the five tying points of Fig. 4, the `1 + a·rs` forms of
   Eqs. (27a-c), and the sign/structure of Eq. (25)'s `d = c_r² + c_s² − 1`.
3. **Audit A3** (the covariant→local metric and the dual basis) — the design's own
   unresolved open item 4, and it feeds every membrane/bending row.
4. **CCX judge for the rest of the failing family** (plate convergence, in-plane lateral,
   modal) — each needs its own `.inp`; the modal one needs CCX's `*FREQUENCY` step.
5. **CCX element-type matrix and multi-material (user requirement).** Extend
   `write_ccx_mesh` with an explicit `S4 | S8 | S8R` selector (the `quadratic` boolean stays
   backward-compatible; CCX's own restriction makes a composite section inherently S8R) —
   **DONE 2026-09-24** for the selector plus the S4/S8/S8R reference of the in-plane bending
   strip. **Still open:** the S4/S8/S8R (+ `*FREQUENCY`) references for the plate, the
   in-plane lateral and the modal cases, and the multi-material (composite-section) rows
   against the element. Tests and docs only; no element change.
6. **Adopt ADR-4 option A** (mesh-consistent nodal directors) once B1-B3 and A3 have been
   audited, so the director change lands on an audited core.
7. **CCX as the closing criterion for the flip**: the S3 gate's Tier 2 plus a CCX row for
   every case in the failing family, all of them recorded in `docs/validation-matrix.md`.
