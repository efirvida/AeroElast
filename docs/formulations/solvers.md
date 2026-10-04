# Solver formulations: code-to-equation reference

Scope of this revision: the time-integration and formulation paths that the
solvers use — Newmark-β (`src/aeroelast/solvers/elasticity/dynamic_newmark.py`,
`crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs`), the
co-rotational and rotating-frame formulation
(`src/aeroelast/solvers/fsi/corotational.py`,
`src/aeroelast/solvers/fsi/rotor.py`,
`crates/aeroelast-solvers/src/petsc/fsi/rotor_physics.rs`), the stress-stiffened
incremental path (`src/aeroelast/solvers/fsi/stress_stiffened_dynamic.py`), and
the element-level mechanics that the solvers consume — selective reduced
integration and the MITC assumed-strain machinery in
`crates/aeroelast-core/src/elements/`. Element-level detail is *pointed at*
`docs/formulations/shell-elements.md`, not duplicated.

Every equation below was read from the source named with it. Where the code names
a work that is not held, the attribution is reported as a repository assertion,
not as verified evidence. Where the code implements a formula with no citation at
all, it is stated.

---

## 1. Scope and conventions

| Code | Role | Status |
| --- | --- | --- |
| `src/aeroelast/solvers/elasticity/dynamic_newmark.py` | Newmark-β, PETSc KSP | documented here; **not cited** |
| `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs` | Newmark-β, PETSc KSP | documented here; **not cited** |
| `src/aeroelast/solvers/fsi/corotational.py` | Rodrigues rotation, frame transforms, fictitious forces | documented here; **not cited** |
| `src/aeroelast/solvers/fsi/rotor.py` | rotating-frame governing equation, K_G, K_SP, torque, ω dynamics | documented here; **ANSYS citation unverifiable** |
| `crates/aeroelast-solvers/src/petsc/fsi/rotor_physics.rs` | gyroscopic `G`, spin softening K_SP, torque | documented here; **ANSYS citation unverifiable** |
| `src/aeroelast/solvers/fsi/stress_stiffened_dynamic.py` | frozen-tangent K_G incremental loop | documented here; **Bathe §6.3 mismatch** |
| `crates/aeroelast-core/src/elements/mitc4.rs` | SRI split, MITC assumed-strain operators | documented here; details in `shell-elements.md` |

Conventions used in this document: `dt` is the time step, `u`, `v`, `a` are the
displacement, velocity and acceleration, `M`, `C`, `K` the mass, damping
(Rayleigh) and elastic stiffness, `K_eff` the effective (dynamic) stiffness and
`F_eff` the effective load. `n` indexes the converged step; `n+1` the unknown.

### 1.1 What "verified"/"unverifiable" means here

- `verified` — read in a held source, or read in the code and cross-checked
  between the Python and Rust implementations.
- `unverifiable` — the code names a work that is **not** in `.sources/papers/`
  and whose equations cannot be reached from any held copy.
- `repository assertion` — the code (or a commit message) names a work; the work
  is not held.

---

## 2. Newmark-β time integration

### 2.1 Python implementation

`DynamicNewmarkSolver` (`dynamic_newmark.py:22-...`). Defaults set in
`_validate_params`: `beta = 0.25`, `gamma = 0.5`, `eta_m = eta_k = 1e-4`. The
class docstring names the scheme — "the implicit constant-average-acceleration
method (β=0.25, γ=0.5 by default)" — and cites nothing. The file even carries a
`# FIXME: This class is not well implemented` on the class.

```text
a0   = 1 / (beta dt^2)
a1_v = 1 / (beta dt)
a1_c = gamma / (beta dt)
a3   = 1 / (2 beta) - 1

K_eff = K + a0 M + a1_c C
F_eff = F(t) + M (a0 u + a1_v v + a3 a)
a_new = a0 (u_new - u) - a1_v v - a3 a
v_new = v + dt [ (1 - gamma) a + gamma a_new ]                          (code)
```

Note, recorded as found: `_compute_effective_force` takes `C` in its signature
but **never uses it**. The Python RHS omits the `C (a1 u + a4 v + a5 a)` term
that the standard scheme requires when `C != 0`, even though `C` is assembled and
`eta_k = eta_m = 1e-4` by default. The Rust path includes it (§2.2); the two
implementations therefore disagree whenever damping is non-zero.

### 2.2 Rust implementation

`newmark_beta_solve` (`dynamic_newmark.rs:...`). The module doc and the function
doc give the scheme and defaults (`beta = 0.25`, `gamma = 0.5`), and no work is
cited.

```text
a0 = 1/(beta dt^2)   a1 = gamma/(beta dt)   a2 = 1/(beta dt)
a3 = 1/(2 beta) - 1  a4 = gamma/beta - 1    a5 = (dt/2)(gamma/beta - 2)

K_eff = K + a0 M + a1 C
rhs   = F_{n+1} + M (a0 u + a2 v + a3 a) + C (a1 u + a4 v + a5 a)
K_eff u_{n+1} = rhs
a_{n+1} = a0 (u_{n+1} - u) - a2 v - a3 a
v_{n+1} = v + dt [ (1 - gamma) a + gamma a_{n+1} ]                      (code)
```

Initial acceleration: with `u0 = v0 = 0`, `a0` is obtained from `M a0 = F0`. The
function doc also prints a predictor/corrector form (`u*`, `v*` and the
`a0..a5` block); the time loop does not use an explicit predictor, but the RHS
above is the equivalent standard effective-load form. `C = eta_k K + eta_m M`
is Rayleigh damping.

**Gap.** Neither implementation cites Newmark, N.M. (1959).
`docs/validation/references.yaml`
lists "A method of computation for structural dynamics", *Journal of the
Engineering Mechanics Division, ASCE*, 85(EM3):67–94, 1959, with DOI **"to
verify"**, notes that there is **no PDF in `.sources/papers/`**, and records the
scheme as "source of an implemented feature; the code does not cite it". The
citation should be added to both source files, and the DOI verified against a
held copy or the published record before it is written down as verified.

### 2.3 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `dynamic_newmark.py` `_assemble_effective_matrix` / `_update_acceleration` / `_update_velocity` | §2.1 | Newmark 1959 — **not cited** | **yes** (code); citation **missing** |
| `dynamic_newmark.rs` `assemble_keff` / `matvec_add` / update block | §2.2 | Newmark 1959 — **not cited** | **yes** (code); citation **missing** |
| Python `_compute_effective_force` | omits the `C` term | - | **yes**, and inconsistent with Rust |

---

## 3. Co-rotational and rotating-frame formulation

### 3.1 `corotational.py` — transforms and fictitious forces

The module has a `Mathematical Formulation` section and **cites nothing**.

```text
x_global = R(theta) (x_ref + u_local)                                   (code)
R(theta) = I + sin(theta) K + (1 - cos(theta)) K^2,   K = skew(axis)    (Rodrigues)
dR/dtheta = cos(theta) K + sin(theta) K^2
F_local = R^T F_global          u_global = R u_local                    (code)
```

`InertialForcesCalculator` implements the fictitious forces of a rotating frame
(rotating frame ⇒ `F_inertial = -m (omega x (omega x r) + 2 omega x v + alpha x r)`):

```text
F_cf    = m omega^2 r_perp,   r_perp = r - (r . n) n                    (code)
F_cor   = -2 m (omega x v)
F_euler = -m (alpha x r)                                                (code)
```

`r = nodal_coords - center`, `omega_vec = omega * axis`, `alpha_vec = alpha *
axis`. The sign convention of the fictitious-force sum is stated in the class
docstring. `OmegaProvider` subclasses (`ConstantOmega`, `RampedOmega`,
`ComputedOmega`, `RampedComputedOmega`, `TableOmega`, `FunctionOmega`) only
prescribe `omega(t)` and `alpha(t)`; they are not formulations.

### 3.2 `rotor.py` — the rotating-frame equation of motion

The module docstring gives the equation "cf. ANSYS MAPDL Theory Reference,
Eq. 14-57, §14.4.1":

```text
[M]{u''} + [C]{u'} + ([K] + [K_G] + [K_SP]){u}
   = {F_aero} + {F_cf} + {F_cor} + {F_euler} + {F_g}                    (rotor.py:14-17)
```

with, as documented in the module:

```text
[C]    = eta_m [M] + eta_k [K]                                          (rotor.py)
[K_G]  = integral B_G^T S~ B_G dA     (from centrifugal prestress
                                        sigma_cf ~ rho omega^2 r L_char) (rotor.py:20-23)
[K_SP] = -omega^2 [M] (I - n tensor n)                                  (rotor.py:24-25)
         (ANSYS Eq. 3-74 / 14-55)
F_cor  = -2 m (omega x v)      F_euler = -m (alpha x r)                 (rotor.py)
```

The module also states the LHS/RHS split of the physical effects — `K_G` and
`K_SP` on the LHS, centrifugal `F_cf` at **undeformed** coordinates `X_0`,
Coriolis explicit and lagged, Euler at deformed coordinates — and the relationship
`[K_total] = [K] + [S] + [S~2]` "(ANSYS §3.4–3.5, Eq. 3-88)", where `[S]` is
stress stiffening and `[S~2]` spin softening.

**Spin softening (`rotor_physics.rs` `build_ksp_vals`)** is implemented per
translational DOF of each node, with `n` the unit rotation axis:

```text
K_SP[dof_base + j] = -omega^2 m_node (1 - n_j^2),  j = 0, 1, 2
K_SP[dof_base + j] = 0,                             j >= 3  (rotational)  (code)
```

The comment heading this function is exactly `ANSYS Eq. 3-74 / 14-55 for lumped
mass` (`rotor_physics.rs:432`).

**Coriolis gyroscopic matrix (`rotor_physics.rs` `build_coriolis_matrix`)** is
built per node as `G_i = -2 m_i Omega`, with `Omega` the skew matrix of `omega`
(`F_cor = G_cor v`); the comment writes the block as

```text
Omega = [  0   -wz   wy ]   w = omega * axis
        [  wz   0   -wx ]
        [ -wy   wx    0  ]                                              (code)
```

**Time integration in the rotor path** reuses the Newmark effective form printed
in the docstring: `K_eff = [K] + [K_G] + [K_SP] + a0 [M] + a1 [C]`,
`F_eff = {F} + [M](a0 u + a2 v + a3 a) + [C](a1 u + a4 v + a5 a)` — note this
matches the **Rust** Newmark form of §2.2, not the Python one.

**Rotor angular dynamics** (also in the module docstring):

```text
I = sum_i m_i r_perp,i^2                       (or user-prescribed)
tau_aero = n . sum_i (r_i x F_cfd,i)
alpha = (tau_driving + tau_shaft) / I
omega^{n+1} = omega^n + alpha dt               (explicit Euler)
theta^{n+1} = theta^n + omega_bar^n dt,  omega_bar^n = omega^n + 0.5 alpha^n dt
                                                                        (code)
```

**Performance coefficients** `Ct`, `Cp`, `Cq`, `TSR` are defined in the module
docstring with the standard normalisation `1/2 rho V_inf^2 pi R^2`; no work is
cited.

### 3.3 Verification status of the ANSYS attribution

`rotor.py` cites the **ANSYS MAPDL Theory Reference** at lines 11, 24, 77, 84,
221, 234, 263 and 287; `rotor_physics.rs` cites it at line 432.

**The ANSYS manual is not a verifiable source and is not held.** It is
proprietary documentation; it cannot be reproduced here, and its equation numbers
(`14-57`, `3-74`, `14-55`, `3-88`) cannot be checked against a held copy. The
user has asked for it to be replaced by the original source.

**Likely original sources** (to be verified before any replacement is made):

- **Géradin, M., Rixen, D.**, *Mechanical Vibrations: Theory and Application to
  Structural Dynamics*, 3rd ed., Wiley, 2015 — listed in `docs/validation/references.yaml` §3
  and used there for "the gyroscopic matrices of rotating systems". This is the
  best candidate for the gyroscopic/Coriolis term `G_cor` and the rotating-frame
  terms.
- **Goldstein, H., Poole, C., Safko, J.**, *Classical Mechanics*, 3rd ed., Addison
  Wesley, 2002, §4.9–4.10 — listed in `docs/validation/references.yaml` for "the non-inertial
  (rotating) frame treatment". This is the candidate for the centrifugal,
  Coriolis and Euler accelerations of §3.1.
- **Bathe, K.J.**, *Finite Element Procedures* — the candidate for the geometric
  stiffness `K_G` (see §4).

**This replacement is not made in this revision.** Each equation must be read in
the original work first; the ANSYS numbers are not evidence of the original
equation numbering, and no equation number here should be transferred to Géradin
& Rixen without reading it. `docs/validation/references.yaml` marks both Géradin/Rixen and
Goldstein with DOI "to verify" and records them as "source of an implemented
feature; the code does not cite it".

### 3.4 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `corotational.py` Rodrigues `R`, `dR/dtheta` | §3.1 | **uncited** | **yes** |
| `corotational.py` frame transforms | §3.1 | **uncited** | **yes** |
| `corotational.py` `compute_centrifugal_force` / `_coriolis_` / `_euler_` | §3.1 | **uncited** | **yes** |
| `rotor.py` governing equation | §3.2 | ANSYS 14-57 | **unverifiable** |
| `rotor_physics.rs` `build_ksp_vals` | §3.2 | ANSYS 3-74 / 14-55 | **unverifiable** |
| `rotor_physics.rs` `build_coriolis_matrix` | §3.2 | docstring only, no external work | **yes** (code) |
| `rotor.py` θ/ω integrations and torque | §3.2 | **uncited** | **yes** |
| `rotor.py` Ct / Cp / Cq / TSR | §3.2 | **uncited** | **yes** |

---

## 4. The stress-stiffened incremental path

### 4.1 What the code does

`stress_stiffened_dynamic.py` (`StressStiffenedFSISolver`) documents its own
method precisely, in its own words:

- "This is a *linearized incremental* (also called "frozen-tangent") approach:
  ONE linear solve per time step — no Newton-Raphson inner iterations; ONE `K_eff`
  factorization per `update_interval` converged steps; converges to the correct
  nonlinear solution as `dt -> 0`."

```text
K_eff^{n+1} = K + K_G(sigma^n) + a0 M + a1 C
sigma^n = membrane stress recovered from u^n                                 (code)
```

The geometric stiffness is assembled (both the Rust fast-path and the Python
`_post_convergence_hook`) as `K_G = integral B_G^T S~ B_G dA`
(`assembler.py` `assemble_geometric_stiffness` docstring). The element kernel is
`compute_geometric_stiffness_from_stress` (`mitc4.rs:1516-1524`):

```text
S_m     = [[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]]
S_tilde = block-diag(S_m, S_m, S_m)                 (copies for u, v, w)
B_G     = 6x24, rows (d/dx, d/dy) of (u, v, w)
K_sigma = sum_{4 Gauss pts} B_G^T S_tilde B_G w sqrt(g),  symmetrised       (code)
```

`sigma` is the **membrane resultant force per unit length** (`cm * eps_m`,
`cm = A`), evaluated with the MITC4+ covariant membrane B-matrix at the element
centre (`compute_membrane_stress`) or per Gauss point
(`compute_geometric_stiffness_local`).

### 4.2 The Bathe §6.3 "Updated Lagrangian" citation does not match the code

The module `References` block (`stress_stiffened_dynamic.py:52-54`) names both

- Ko, Lee & Bathe (2017), "The MITC4+ shell element in geometric nonlinear
  analysis", *Computers & Structures* 185:1-14, and
- **Bathe (1996), *Finite Element Procedures*, §6.3 Updated Lagrangian.**

§6.3 of Bathe is the updated-Lagrangian (UL) incremental procedure: strains and
stresses are referred to the **last computed configuration**, and the tangent
stiffness is reassembled and **iterated** (Newton-Raphson) to equilibrium at each
increment. The code does **not** do that:

- it performs **one linear solve per step, no inner iterations** (its own
  docstring, §4.1);
- the stiffness it adds is a total-Lagrangian geometric stiffness:
  `compute_kt_global` in `mitc4.rs:1408-1462` is documented and implemented as
  `K_T = K_0 + K_L + K_sigma (Total Lagrangian)`, with the linearization taken
  about the **undeformed** configuration.

So the code is a frozen-tangent incremental scheme using a total-Lagrangian
`K_sigma`, and the "Updated Lagrangian" attribution is a code-comment claim that
does not describe the implementation.

**Verification status:** Bathe (1996) is not held, so §6.3 cannot be read here.
What is verifiable is the mismatch between the citation and the code's own
described method. `docs/validation/references.yaml` lists Bathe, *Finite Element Procedures*,
**2nd ed., 2014** at `stress_recovery.py:53`, while this file cites the **1996**
edition; the edition discrepancy is another item to resolve. The likely original
source for the geometric stiffness itself is Ko, Lee & Bathe (2017) — a held PDF
(`.sources/papers/mitc4+_no_lineal.pdf.pdf`) — but that equation was **not**
extracted for this revision, so it is not claimed here.

### 4.3 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `stress_stiffened_dynamic.py` + Rust loop | §4.1 | Ko et al. 2017 + Bathe 1996 §6.3 | frozen-tangent construction **yes** (code); Bathe §6.3 **unverifiable and mismatched** |
| `compute_geometric_stiffness_from_stress` / `_local` | §4.1 | none | **yes** |
| `compute_kt_global` | `K_0 + K_L + K_sigma (Total Lagrangian)` | none (docstring says TL) | **yes** (code) |

---

## 5. Selective reduced integration (SRI)

### 5.1 The split exactly as implemented

`compute_ke_local` (`mitc4.rs:1107`, comment block starting "Despite MITC4+
blending, residual in-plane shear locking persists ...", reference line "Hughes,
Taylor & Kanoknukulchai (1977) — SRI for Q4 membrane element", `mitc4.rs:1123`):

```text
cm_normal   = cm, then
              cm_normal[(0,2)] = cm_normal[(1,2)] = 0
              cm_normal[(2,0)] = cm_normal[(2,1)] = cm_normal[(2,2)] = 0
c_shear     = cm[(2,2)]

k_m  = sum_{4 Gauss pts} Bm(xi,eta)^T cm_normal Bm(xi,eta) (w sqrt(g))
k_m += b_shear^T c_shear b_shear (4 sqrt(g_c)),
       b_shear = row 2 of b_m_mitc4_plus at (xi, eta) = (0, 0),
       sqrt(g_c) = |g_r x g_s| at the centre
k_m  = 0.5 (k_m + k_m^T)                                                (code)
```

That is: **full 2×2 for the normal membrane components** (`eps_xx`, `eps_yy` and
their coupling, the shear row/column struck out) and **one centre point at
`xi = eta = 0` for the in-plane shear** (`2 eps_xy`), with weight `4 sqrt(g_c)`
(the total weight of the 2×2 rule). The same split is repeated in
`compute_fint_global` for Newton consistency. This is a *data* description of the
scheme; the element treatment lives in `docs/formulations/shell-elements.md` §2.5.

### 5.2 Verification status

The comment cites **Hughes, T.J.R., Taylor, R.L., Kanoknukulchai, W.**, "A simple
and efficient finite element for plate bending", *IJNME* 11(10):1529–1543, 1977,
and commit `929db32` repeats `Ref: Hughes, Taylor & Kanoknukulchai (1977)`.

**No copy of that paper is held in `.sources/papers/`.** The claim that this paper
prescribes exactly a 2×2 rule for the normal membrane components and a one-point
rule for the in-plane shear **cannot be verified here**, and must be treated as a
**repository assertion**, not as evidence. `docs/validation/references.yaml` also marks its
DOI "to verify" and states that the journal, volume and pages come from the
published record and are not re-verified against a held copy.

### 5.3 The measured consequence for the drill-membrane term

The split has a consequence already measured elsewhere (recorded in
`mitc4.rs:802-825` and `docs/formulations/shell-elements.md` §2.5/§4.2): the
MITC4/D drill-membrane strain of Ko, Bathe & Zhang (2025) is a **pure in-plane
engineering shear** — its `rr` and `ss` rows vanish identically and only the `rs`
row survives. Because SRI zeroes the in-plane shear in the 2×2 loop and integrates
it at the centre, where all mid-side derivatives of the drill interpolation
vanish (`b_md(0,0) = 0`), the drill term contributes nothing to the assembled
stiffness. It was measured as bit-identical (sha256 of the COO values) with the
operator wired in and amplified by 1000. This is why the MITC4/D drill-membrane
term is inert in this element; element-level detail is in
`docs/formulations/shell-elements.md`.

### 5.4 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `compute_ke_local` SRI block | §5.1 | Hughes, Taylor & Kanoknukulchai 1977 (`mitc4.rs:1123`) | split **yes** (code); attribution **repository assertion** |
| `compute_fint_global` SRI block | §5.1 | same | **yes** (code) |
| drill-membrane inertness | §5.3 | measured, `mitc4.rs:802-825` | **yes** (design note) |

---

## 6. The MITC assumed-strain machinery

This section records only *which paper ties each operator to which equation*. The
element-level correspondence is in `docs/formulations/shell-elements.md`; it is
not repeated here.

| Operator | File | Paper | Equation | Status |
| --- | --- | --- | --- | --- |
| `b_m_mitc4_plus` | `mitc4.rs:657` | Ko, Lee & Bathe 2017, "A new MITC4+ shell element", *C&S* 182:404–418 | Eqs. (27a), (27b), (27c), blending `aA..aE` printed after (27c) from Eq. (24) | **verified term by term** in `shell-elements.md` §2.3 |
| `b_gamma_mitc4` | `mitc4.rs:963` | Dvorkin & Bathe 1984, *Eng. Comput.* 1(1):77–88 | Eq. (3), journal p. 78 | **equation now read from the scan image**; correspondence detail in `shell-elements.md` §2.4 |
| `b_gamma_ext` | `mitc3.rs:468` | Lee, Lee & Bathe 2014, *C&S* 138:12–23 | constant part Eq. (15), linear part Eq. (16), total Eq. (17) | **verified term by term** in `shell-elements.md` §3.3–3.4 |
| `b_gamma_ext` sign convention | `mitc3.rs:393` | Ko, Bathe & Zhang 2025, *C&S* 308:107622 | Eq. (3a) | **verified** in `shell-elements.md` §3.6 |
| `b_md_mitc4_plus` `drill_midside_shape_derivatives` | `mitc4.rs:738`, `mitc4.rs:829` | Ko, Bathe & Zhang 2025 | Eqs. (5), (10), (11), (13c), (17a-b), (18) | **verified** in `shell-elements.md` §4.2 |

### 6.1 Note on the Dvorkin & Bathe 1984 equation

`docs/formulations/shell-elements.md` §2.4 recorded that the displayed equations
of the 1984 paper "could not be read" because the held scan has no text layer over
them. The scan *does* carry them as an image, and rendering journal page 78 at
300 dpi and reading it recovers the equation that `shell-elements.md` reported as
unreadable and therefore unnumbered. It is **Eq. (3)**, journal p. 78:

```text
e~13 = 1/2 (1 + r2) e~A13 + 1/2 (1 - r2) e~C13
e~23 = 1/2 (1 + r1) e~D23 + 1/2 (1 - r1) e~B23                            (3)

with the tying points A, C on the r1-constant pair and B, D on the r2-constant
pair (Fig. 2, "Interpolation functions for the transverse shear strains"),
    source: Dvorkin & Bathe 1984, journal p. 78 (PDF p. 2), read as an image.
```

The same page also prints the displacement interpolation as Eq. (2),
`u_i = sum_k h_k u_i^k + (r3/2) sum_k a_k h_k ( -^0V_k^2 alpha_k + ^0V_k^1 beta_k )`,
whose rotation parameters `alpha_k`, `beta_k` are the subject of the open sign
question in `shell-elements.md` §4.5. This document does not resolve that
question; it only fixes the transverse-shear equation number and form, which
should be folded back into `shell-elements.md` §2.4.

---

## 7. Gaps: uncited implementations and unverifiable citations

1. **Newmark-β is implemented with no citation** in
   `dynamic_newmark.py` and `dynamic_newmark.rs`. Newmark (1959) is the source;
   `docs/validation/references.yaml` lists it with DOI "to verify" and no held copy. **Fix
   needed in the code.**
2. **The Python Newmark RHS omits the damping term** while the Rust RHS includes
   it. With the default `eta_k = eta_m = 1e-4` this is an active inconsistency.
3. **`corotational.py` cites nothing**, despite carrying a "Mathematical
   Formulation" section: Rodrigues' formula, the frame transforms and the
   centrifugal/Coriolis/Euler forces have no source.
4. **`rotor.py` and `rotor_physics.rs` cite the ANSYS MAPDL Theory Reference**,
   which is not verifiable. The user wants it replaced by the original source;
   the likely original is Géradin & Rixen (gyroscopic and rotating-frame terms),
   with Goldstein, Poole & Safko for the non-inertial frame treatment. **The
   replacement must be verified against the original before it is made**; none of
   the ANSYS equation numbers may be carried over unverified.
5. **`stress_stiffened_dynamic.py` cites Bathe (1996) §6.3 "Updated Lagrangian"
   but implements a frozen-tangent incremental scheme on a total-Lagrangian
   geometric stiffness** (`compute_kt_global` says "Total Lagrangian"). The
   citation and the code disagree; Bathe is not held, so §6.3 itself is
   unverifiable. The edition also disagrees with `docs/validation/references.yaml` (2014 vs
   1996).
6. **The SRI attribution is a repository assertion.** Hughes, Taylor &
   Kanoknukulchai (1977) is named in the code and in commit `929db32`, but no
   copy is held and the specific claim (2×2 for normal membrane, one point for
   in-plane shear) could not be checked.
7. **`b_gamma_mitc4`'s normalisation** — the `(1/8) |g_r x g_s|^{-1}` scaling and
   the `alpha`/`beta` angle pair — still cannot be tied to the original paper's
   normalisation, because Eq. (3) has no scaling; the correspondence remains open
   in `shell-elements.md` §2.4.
8. **The drilling penalty** (`b_drill`, `k_drill = 0.15 E h^2 drilling_scale`) is
   element code, not solver code, and is documented as uncited in
   `docs/formulations/shell-elements.md` §2.6. It is re-listed here only because
   the SRI split (§5) is what makes the connected drill-membrane term inert.

No other solver-side formula was found to be implemented without attribution.

---

## References

`docs/validation/references.yaml` is the canonical bibliography. Where this document and that
file disagree, `docs/validation/references.yaml` should be corrected; this document does not
duplicate its per-entry verification annotations.

1. Newmark, N.M., "A method of computation for structural dynamics", *Journal of
   the Engineering Mechanics Division, ASCE*, 85(EM3):67–94, 1959. DOI: to
   verify. **Not held; the code cites nothing.** §2.
2. Dvorkin, E.N., Bathe, K.-J., "A continuum mechanics based four-node shell
   element for general non-linear analysis", *Engineering Computations*
   1(1):77–88, 1984. DOI: to verify.
   `.sources/papers/A_Continuum_Mechanics_Based_Four-Node_Shell_Element_for_General_Nonlinear_Analysis.pdf`.
   Eq. (3) read from the scan image (journal p. 78). §6.
3. Ko, Y., Lee, P.-S., Bathe, K.-J., "A new MITC4+ shell element", *Computers
   and Structures* 182:404–418, 2017, doi:10.1016/j.compstruc.2016.11.004.
   `.sources/papers/1-s2.0-S0045794916309464-main.pdf`. §6.
4. Lee, Y., Lee, P.-S., Bathe, K.-J., "The MITC3+ shell element and its
   performance", *Computers and Structures* 138:12–23, 2014,
   doi:10.1016/j.compstruc.2014.02.005.
   `.sources/papers/The_MITC3+_shell_element_and_its_performance.pdf`. §6.
5. Ko, Y., Bathe, K.-J., Zhang, X., "Continuum mechanics-based shell elements
   with six degrees of freedom at each node — the MITC4/D and MITC4+/D elements",
   *Computers and Structures* 308:107622, 2025,
   doi:10.1016/j.compstruc.2024.107622.
   `.sources/papers/1-s2.0-S0045794924003511-main.pdf`. §6.
6. Ko, Y., Lee, P.-S., Bathe, K.-J., "The MITC4+ shell element in geometric
   nonlinear analysis", *Computers and Structures* 185:1–14, 2017,
   doi:10.1016/j.compstruc.2017.01.015.
   `.sources/papers/mitc4+_no_lineal.pdf.pdf`. Equation **not extracted** in this
   revision. §4, §6.
7. Hughes, T.J.R., Taylor, R.L., Kanoknukulchai, W., "A simple and efficient
   finite element for plate bending", *International Journal for Numerical
   Methods in Engineering* 11(10):1529–1543, 1977. DOI to verify. **Not held —
   repository assertion.** §5.
8. Bathe, K.J., *Finite Element Procedures*, Prentice Hall, 1996 (§6.3 Updated
   Lagrangian) / 2nd ed. 2014. DOI: to verify. **Not held.** §4.
9. Géradin, M., Rixen, D., *Mechanical Vibrations: Theory and Application to
   Structural Dynamics*, 3rd ed., Wiley, 2015. DOI: to verify. **Not held**;
   likely original for the gyroscopic and rotating-frame terms. §3.
10. Goldstein, H., Poole, C., Safko, J., *Classical Mechanics*, 3rd ed., Addison
    Wesley, 2002, §4.9–4.10. DOI: to verify. **Not held**; likely original for
    the non-inertial frame treatment. §3.
11. ANSYS MAPDL Theory Reference — proprietary, **not a verifiable source, not
    held**, to be replaced by items 9 and 10 and/or 8. §3.

`docs/formulations/shell-elements.md` is the element-level companion to this
document and carries the paper-page-level derivation of the MITC operators
referenced in §6.
