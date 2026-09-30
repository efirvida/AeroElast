# Feature: MITC4 warped-quad enrichment (the ERC drill fix)

## Objective

Make the MITC4 quad reproduce the published MITC4 value on the thin twisted beam
(MacNeal-Harder), so the two `xfail(strict)` cases in
`tests/test_ko2017_performance.py` flip.

## Status: the ERC formulation WORKS (measured)

Winkler & Plakomytis' ERC (enhanced rotation constraint) plus the element's
existing SRI takes the thin twisted beam from **0.7313 to 1.0016** (paper
0.9959) at N=8, and frees the drilling penalty.

| case | baseline | ERC beta=1 | ERC beta=0.1 | ERC beta=0.01 | paper |
| --- | --- | --- | --- | --- | --- |
| thin N=8 in-plane | 0.7313 | **1.0016** | 0.9977 | 0.9977 | 0.9959 |
| thick N=16 in-plane | 1.0001 | **0.9984** | 1.0014 | 1.0317 | 0.9972 |

Every number above was reproduced by the parent, not only by the writer.

**Reading.** The ERC EAS frees the drill constraint: thin is flat for
`beta <= 0.1` and within 0.18% of the paper.  SRI is load-bearing: without it the
ERC block uses the full constitutive and the thin beam collapses to 0.0062
(classic Q4 shear locking).  MITC4/D and SRI were "mutually exclusive" only
because the earlier probe kept the drill in the *local* constitutive slot; the
ERC's split constitutive `C4 = blkdiag(A, beta*t*mu/4)` composes with SRI.  The
thick beam is *not* flat as `beta -> 0` (it softens), so `beta` must stay
`>= 0.1`; the paper uses `beta = 1`, where both cases land within 0.6%.

### E1 (done)

`b_erc` (4x24, LOCAL frame) and `erc_covariant_to_local` (4x4) in
`crates/aeroelast-core/src/elements/mitc4.rs`, with 7 tests.  `b_erc` rows 0..2
are the MITC4+ assumed membrane; row 3 is `-2 * b_drill`, so
`row3 * u = 2 e[12] = 2 c` (Winkler Eq. 104/115).

### E2 (probe done, not wired yet)

`m8_erc` (Eq. 103), `regularized_inverse_8x8`,
`compute_ke_local_erc(pre, beta, use_sri)` and one rigid-body/symmetry test.  No
production path calls it yet; `compute_ke_local` is byte-identical to before.

### Source ambiguity resolved (the part the paper gets wrong)

Eq. (102) as printed gives rows 3-4 equal to the plain `e12` and `e21` transform
rows, so it CANNOT produce the Eq. (107) left-hand side.  The consistent reading:
`J0^(4x4)` is the plain nonsymmetric transform on `[e11,e22,e12,e21]^(A)`, and
the ERC vector `[e11,e22,2e(12),2e[12]]` is formed afterwards by
`S = [[1,0,0,0],[0,1,0,0],[0,0,1,1],[0,0,1,-1]]`.  Hence
`erc_covariant_to_local = S * J0^(4x4)`, and `M^(8)` is applied to the plain
covariant components.  Row 3 of the transform, `det(j_inv) * (e12 - e21)`,
reproduces the correct antisymmetric transform `e[12]^(T) = det(J) e[12]^(A)`.

### The previous probe's five errors (all fixed here)

1. drill modulus `0.15 E h^2` -> `beta*t*mu/4` (`= beta*cm[2][2]/4` in the ERC
   vector);
2. the missing `/4` bookkeeping of the `2e[12]` component;
3. the constitutive must act in the SPLIT basis (`blkdiag(A, beta*t*mu/4)`), not
   `blkdiag(cm, k_drill)` on the nonsymmetric vector -- this was the real reason
   the first probe measured no change;
4. `b_erc` row 3 sign/factor (`-2*b_drill`, giving `2e[12]`, not `b_drill`);
5. (frame note) the two readings of Eq. (107) -- covariant-then-transform vs
   local -- are EQUIVALENT for the compatible strain, so "frame mixing" was a red
   herring.

### New pre-existing defect found (independent of the ERC)

The drill constraint is NOT rigid-body invariant on warped geometry: a rigid
rotation carries the nodal out-of-plane offset `z_I` into the local in-plane
displacement, so `1/2*(du/dy - dv/dx) = -1/2*omega*dz/dy != 0` while `thz = 0`.
Measured: the ERC drill row leaves ~1.1e-3 (relative) on `make_pre_warped()`
against ~1e-16 on the flat element, and `compute_ke_global` on the warped element
leaves a rigid-body rotation residual of ~6e-7 against ~1e-18 flat.  Winkler &
Plakomytis' own remedy is the nodal warping correction of Eqs. (130)-(132)
(`u_l = u_G + skew(w_I t3) theta_G`, `w_I` = distance to the middle surface).
Pinned by `test_b_erc_drill_row_warping_residual_is_characterized`.

## Status: E3 DONE — the ERC is the production path, Rust suite green

The ERC is wired into `compute_ke_local` / `compute_fint_global` /
`compute_element_stress`, and `cargo test -p aeroelast-core` is **120 passed / 0
failed** (112 before this feature + 8 new tests).

| case | baseline | ERC | paper |
| --- | --- | --- | --- |
| thin N=8 in-plane | 0.7313 | **0.9976** | 0.9959 |
| thin N=16 in-plane | 0.9131 | **0.9982** | 0.9975 |
| thin N=16 out-of-plane | 0.9112 | **0.9986** | 0.9980 |
| thick N=16 in-plane | 1.0001 | **0.9984** | 0.9972 |

Every number above was reproduced by the parent.  The paper's own LFS4-ERC
converges to `5.2463 / 5.256 = 0.9982` on the same thin pre-twisted beam, so this
matches *their* element, not only MITC4.

### The two fixes that unblocked E3

1. **Symmetric rank-revealing pseudo-inverse** for the 8x8 enhanced-strain block
   (`pseudo_inverse_8x8`, `SymmetricEigen`, keep `lam > 1e-10 * max|lam|`).  The
   previous `regularized_inverse_8x8` injected a rank-1 NEGATIVE term (the block is
   singular under SRI because `m8_erc(0,0) = 0` leaves the enhanced in-plane shear
   unpenalized at the centre) and broke symmetry.  With the pseudo-inverse the
   `lmin/lmax = -7.4e-10` indefiniteness and the rigid-body residual disappear.
2. **The Eq. (110) warping penalty**, `beta_w * (t^3/12) * mu * K_3a * dK_3a` with
   `K_3a = d(theta_z)/dx_a` (the drill-rotation gradient).  This is the paper's own
   term, and `beta_w = 0.01` is their Table 1 value for LFS4-MP.  It removes exactly
   the spurious non-constant `theta_z` zero modes: flat 9 -> **exactly 6** (all six
   rigid-body), warped 10 -> 7.  It is genuinely non-critical: the four benchmarks
   move by <= 1e-4 across `beta_w` in `{0.001, 0.01, 0.1}`, matching their
   Figs. 3/7/9/21.  LFS4-ERC itself sets `beta_w = 0`; the nonzero value is a
   well-posedness stabiliser taken from the paper's sibling element, not a fit.

The drill row and the warping operator now share one Jacobian (`drill_dh`), so the
constraint and the stabiliser differentiate the same `theta_z` field.

### Confirmed from the paper (kills two earlier hypotheses)

Table 1: **LFS4-ERC uses `beta_t = beta_w = 0` and `beta = 1`.**  So the
drill-curvature/twist terms are not the missing piece in the paper's reference
configuration.  And the mode-set screen showed the enhanced-mode set is not the
lever either: `row3 = row2` (no drill enhancement) gives exactly 6 zero modes but
collapses the thin beam to 0.0050 — the thin pre-twisted beam genuinely needs the
released drill.

## Remaining work

- **E4** two real Python failures survive (`pytest -m "not slow"` = 343 passed /
  4 failed, of which 2 are the desired XPASS):
  1. `test_rust_composite.py::TestBatchComposite::test_batch_ke_mitc4_multiple` —
     a negative eigenvalue (`~ -4.7e2`) on the WARPED element of the batch, with a
     composite layup.  Candidate cause: the ERC's softer released drill lets the
     indefinite membrane-bending `B` coupling dominate.  Measure before fixing.
  2. `test_shell_convergence.py::test_in_plane_bending_convergence` — measured
     order 1.16 against the required 1.5.

  Also open: the warped element still has ONE spurious `theta_z` mode, the
  pre-existing warped-geometry drill non-invariance; the paper's remedy is the
  nodal warping correction of Eqs. (130)-(132) (`u_l = u_G + skew(w_I t3) theta_G`).
- **E5** flip the two `xfail(strict)` markers (they now XPASS), update the extract
  and the validation matrix, and re-baseline the Python expectations the ERC
  legitimately moved.

---

<details><summary>Historical record: the refuted handoff hypothesis</summary>

## Status (superseded): open, and the handoff hypothesis is refuted

The three unchecked pieces of Ko, Lee & Bathe (2017) were implemented and
measured. The measurement refutes the handoff's hypothesis: the `∂x_b` warping
term of Eq. (7c) is *not* the fix, and the dominant defect is the drilling
treatment.

### What was implemented, faithfully

- `nodal_normals: [Vector3<f64>; 4]`, the nodal directors `V_n^i` of Eq. (8a)
  (`normalize(g_r × g_s)` at each node; `a_i` is the shell thickness, confirmed on
  paper A p. 405 Eq. (1): *"a_i and V_n^i denote the shell thickness and the
  director vector at the node"* — the extract's "edge vectors" note is wrong).
- `b_kappa_covariant(pre, xi, eta)`: Eq. (7c) in full, both the rotation pair
  `½(x_m,i·u_b,j + x_m,j·u_b,i)` and the warping pair
  `½(x_b,i·u_m,j + x_b,j·u_m,i)`, with `x_b = ½Σa_i h_i V_n^i` (Eq. 8a) and
  `u_b = ½Σa_i h_i(−V_2^i α_i + V_1^i β_i)` (Eq. 8b). Verified to reduce exactly
  to `b_kappa` on a flat element (unit test) and to cancel identically under a
  rigid-body motion (algebraic check).

### Measured result (thin twisted beam, `t/L = 0.0002667`, in-plane)

| configuration | N=8 | paper |
| --- | --- | --- |
| baseline `b_kappa` | 0.7313 | 0.9959 |
| `b_kappa_covariant` (∂x_b warping on) | **0.4609** | 0.9959 |

The warping term makes the element much *stiffer*, the wrong direction. The
rotation pair alone (`b_kappa_covariant` with the translation block disabled)
reproduces `b_kappa` to the digit (0.7313), so the whole difference is the
warping pair. **Refuted.**

### The real defect: the drilling treatment

Scaling the drilling penalty `k_drill` in `compute_ke_local` shows it dominates
the warped-quad response:

| `drilling_scale` | thin N=8 | thick N=16 |
| --- | --- | --- |
| 1.0 (current) | 0.7313 | 1.0001 |
| 0.1 | 0.9619 | 1.0265 |
| 0.03 | 0.9874 | 1.0926 |
| 0.01 | 0.9966 | 1.2786 |
| 0.003 | 1.0054 | 1.9021 |

Thin wants a soft penalty, thick wants a stiff one. The required stiffness scales
as `E·h³` for both (the current code uses `0.15·E·h²`, off by one factor of `h`):
thin wants `k ≈ 0.44`, thick wants `k ≈ 4.5e5`, and `E·h³` is `0.95` and `9.5e5`
respectively, so a single `k ≈ 0.46·E·h³` fits both. **This is a measured
hypothesis, not a derivation, and the user has (correctly) asked for the
derivation from the papers rather than a fitted constant.**

### The paper-faithful 6-DOF drill operator does not compose cleanly

`b_md_mitc4_plus` was verified against Ko, Bathe & Zhang (2025), Eq. (18), term by
term (`coeff_rr = (j0/j) h̃_r c_r`, `coeff_ss = −(j0/j) h̃_s c_s`,
`coeff_rs = ½(j0/j)(h̃_s c_r − h̃_r c_s)`), and the paper has **no drilling
penalty**: the drill DOF is carried by this assumed drill-membrane strain,
integrated with the full constitutive matrix. Wiring it that way:

| configuration | thin N=8 | thick N=16 |
| --- | --- | --- |
| `b_m + b_md`, full `cm`, no SRI, no penalty | 0.9951 | **3.77** |
| `b_m` only, full `cm`, no SRI | 1.2420 | **18.26** |
| `b_m` + SRI + `b_md` shear at 4 GP, tiny penalty | 0.9973 | **3.73** |

The full-integration membrane is required by the paper and is exactly what our
element cannot use: it needs the SRI split (a deviation from the paper, commit
`929db32`) to keep the thick case, and the SRI zeroes the shear row at the 4
Gauss points, which is where the drill shear lives. The MITC4/D operator and our
SRI are mutually exclusive. This is the architectural conflict the independent
review flagged: a flat-projection Mindlin/ABD shell cannot adopt a
continuum-mechanics drill operator piecemeal.

### Corrected paper reading (independent of the fix)

Ko, Lee & Bathe (2017), Eqs. (21)–(27), *do* carry five geometry-dependent
coefficients `a_A..a_E`; commit `dc7593e` removed them as "unsourced" and that
reading was wrong. For the twisted-beam mesh the quads are ruled, `x_d` is purely
out-of-plane, so `c_r = c_s = 0` and all five vanish — which is why removing them
was bit-identical. They are not this fix, but the extract and the code comment
must stop claiming the paper has no such coefficients.

## Next step (principled, not fitted)

### What the literature actually says (derived, not tuned)

- Hughes & Brezzi (1989) is a **mixed variational principle with an independent
  rotation field**, not a material stiffness. Their own descendants report it is
  *too stiff* numerically.
- The drill constraint is `c = ω + ½(u_{1,2} − u_{2,1}) = 0` — identical to our
  `b_drill` — imposed by a **penalty** `F = β/2·c²` (Winkler & Plakomytis,
  ECCOMAS 2016, Eqs. 115–117).
- LS-DYNA's production drilling constraint (Erhart & Borrvall, 9th European
  LS-DYNA Conf. 2013) uses a **fitted** `K = 0.005·k·E·V·BᵀB` and calls it
  *"a non-physical additional stiffness"*.
- Winkler & Plakomytis: *"Formulations involving a rotation constraint imply the
  application of a problem-dependent penalty or regularization parameter"*, and
  the fix is **not** the value: with properly designed **enhanced strain fields**
  *"the results are independent of the penalty parameter ... this threshold is
  independent of the problem"*. Also: *"the stiffness related to the nonsymmetric
  membrane strain plays the role of a penalty parameter rather than that of a
  material parameter"*.

**Conclusion.** There is no unique `k_drill` to derive: it is a penalty, and the
measured thin/thick conflict (soft for thin, stiff for thick) is the symptom of a
**missing enhanced strain field**, not of a wrong constant. The `E·h³` scaling
above is a diagnostic, not a fix.

### Concrete next step

The principled fix is Winkler & Plakomytis' **ERC** (enhanced rotation
constraint). Their LFS4-ERC element matches the thin pre-twisted beam
(`h = 0.0032`) to 0.1% at 4x24 and is flat against the penalty parameter over
`β ∈ [1e-8, 1e8]` (their Fig. 17, Tables 9–11) — exactly our failing case.

**The formulation, verified from the papers (Winkler & Plakomytis, ECCOMAS 2016):**

- Modified Hu–Washizu functional (Eq. 108):
  `∫ [ t Q^abcd (ε_ab + ε̃_ab)(ε_cd + ε̃_cd) + t³/12 Q^abcd κ_ab κ_cd
     + κ_s t γ_a γ_a + κ t c c ] √A dξ¹dξ² = W_ext`,
  with `c = E_[12]` the drill constraint and `ε_ab` the full **4-component**
  non-symmetric strain.
- ERC strain vector (Eq. 104): `ε̃ = [ε11, ε22, 2ε(12), 2ε[12]]ᵀ`.
- Transformation (Eq. 107): `ε̃^(T) = J0^(4×4) ε̃^(A)`, with `J0^(4×4)` (Eq. 102):
  `[(J1¹)², (J1²)², J1¹J1², J1²J1¹; (J2¹)², (J2²)², J2¹J2², J2²J2¹;
    J1¹J2¹, J1²J2², J1¹J2², J1²J2¹; J2¹J1¹, J2²J1², J2¹J1², J2²J1¹]`.
- Enhanced strain (Eq. 96): `ε̃ = (j0/j) J0 M α_e`.
- **M^(8)** (Eq. 103), 4×8, the mode set that performs best:

  ```text
  M^(8) = [ ξ¹  0   0   0   ξ¹ξ² 0    0    0    ]
          [ 0   ξ²  0   0   0    ξ¹ξ² 0    0    ]
          [ 0   0   ξ¹  0   0    0    ξ¹ξ² 0    ]
          [ 0   0   0   ξ²  0    0    0    ξ¹ξ² ]
  ```

  (The symmetric-only `M^(7)`, Eq. 100, is the Andelfinger–Ramm 7-mode EAS; the
  ERC splits the `2ε12` slot into `2ε(12)` and `2ε[12]` and gives each its own
  mode. `M^(8)` is from Chróścielewski & Witkowski 2006, whose EAS(14)m1 uses 14
  parameters, `H` is 6×14, and enhances `{e11, e22, e12, e21}`; bending is **not**
  enhanced — "semi-EAS".)
- Static condensation (Eqs. 5.38–5.45 of C&W):
  `K = K_M + K_G − K_βqᵀ K_ββ⁻¹ K_βq`,
  with `K_ββ = ∫PᵀCP`, `K_βq = ∫PᵀCB`, `P = (j0/j) J0 M`.
- The drill constraint `c = ω + ½(u_{1,2} − u_{2,1}) = 0` (Winkler Eq. 115) is
  the same operator as our `b_drill`; it is enforced by `κ t c c`, and the EAS
  enhancement of `2ε[12]` is what removes the penalty sensitivity.

**Work units**

- **E1** Add the 4-component ERC strain (`b_erc`): the compatible 4-component
  strain `[ε11, ε22, 2ε(12), 2ε[12]]` from the existing MITC4+ membrane plus the
  drill constraint, transformed with `J0^(4×4)`.
- **E2** Add the `P = (j0/j) J0 M^(8)` enhanced operator and the `K_ββ`, `K_βq`
  condensation; unit-test that the penalty becomes insensitive (the Winkler
  Fig. 2/17 flatness) and that the flat rigid-body/nullspace invariants hold.
- **E3** Wire into `compute_ke_local` / `compute_fint_global` / `compute_kt_global`
  with the existing consistency guards.
- **E4** Measure the thin and thick twisted beams; expect ~0.996 for both, flat
  against `k_drill`.
- **E5** Flip the `xfail` markers, update the extract and the validation matrix.

This is an element rewrite (the membrane becomes a mixed EAS field, the static
condensation grows from 2 to 8 incompatible modes, and the drill penalty stops
being the load path). It should be done as its own tracked feature, not as a
patch.

### Measured attempt at E1+E2 (this session, reverted)

A first probe was implemented in `compute_ke_local`: the 4-component compatible
strain `[e11, e22, 2e(12), 2e[12]]` from `b_m_mitc4_plus` + `b_drill`, the
`J0^(4×4)` transform built from `j_inv` (Eq. 102), the `M^(8)` modes (Eq. 103),
and the 8-mode static condensation with `C = blkdiag(cm, k_drill)`. Result:

| `drilling_scale` | thin N=8 (probe) | thin N=8 (baseline) |
| --- | --- | --- |
| 1.0 | 0.7314 | 0.7313 |
| 0.1 | 0.9619 | 0.9619 |
| 0.01 | 0.9967 | 0.9966 |
| 0.001 | 1.0241 | 1.0209 |

The penalty sensitivity is **unchanged**, so the EAS did not free the drill
constraint. What is missing, read with vision from Winkler's Eq. (110):

- the drilling term is `β t μ (E_[12] + ε̃_[12])(δE_[12] + δε̃_[12])` — the
drilling modulus is **`β·t·μ` (β·G·h)**, not our `0.15·E·h²`;
- the ERC strain vector's 4th entry is `2ε[12]` while the penalty is written on
`E_[12]`, so the 4×4 constitutive entry is `β t μ / 4` (bookkeeping the probe got
wrong);
- the compatible 4-component strain must be the **covariant** strain transformed
by `J0^(4×4)` (Eq. 107), not the already-local `b_m`/`b_drill`; the probe mixed
the frames.

Also: `M^(8₁)`/`M^(8₂)` (Eqs. 105–106) define ERC1/ERC2, and `M^(8)` (Eq. 103)
is the one that performs best (Winkler Fig. 2). The reference implementation uses
`β = 1` (their Table 1: `LFS4-ERC ... β = 1`).

## Constraints

- Vision reads only (`pdftoppm`), never `pdftotext` for equations.
- Never run the benchmarks at N ≥ 32.
- Never weaken a guard test to make the target pass.
- No parameter fitting against a benchmark.

</details>
