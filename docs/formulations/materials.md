# Constitutive and failure models: code-to-equation reference

Scope of this revision: the material constitutive models and the failure criteria
in `crates/aeroelast-core/src/materials/`, plus their Python mirrors
`src/aeroelast/core/laminate.py` and `src/aeroelast/constitutive/failure.py`.
Every equation below was read from the source named with it — the code itself, or
a PDF held in `.sources/papers/`. Where the code attributes a formula to a work
that is **not** held, the attribution is reported as a repository citation and
not as verified evidence. Nothing is reconstructed.

---

## 1. Scope and conventions

### 1.1 Files covered in this version

| File | Contents | Status |
| --- | --- | --- |
| `crates/aeroelast-core/src/materials/laminate.rs` | CLT: plies, ABD, `Cs`, `to_shell_constitutive*` | documented here |
| `crates/aeroelast-core/src/materials/orthotropic.rs` | plane-stress `Q`, `Qbar`, `Cbar_shear` | documented here |
| `crates/aeroelast-core/src/materials/composite.rs` | ABD → `ShellConstitutive` constructor | documented here |
| `crates/aeroelast-core/src/materials/isotropic.rs` | isotropic plane-stress shell law | documented here |
| `crates/aeroelast-core/src/materials/failure.rs` | Tsai-Wu, Hashin, max-stress | documented here |
| `crates/aeroelast-core/src/materials/mod.rs` | `ShellConstitutive`, `Material` trait | documented here |
| `src/aeroelast/core/laminate.py` | Python mirror of the CLT and `Cs` code | documented here |
| `src/aeroelast/constitutive/failure.py` | Python mirror of the failure criteria | documented here |
| `src/aeroelast/core/material.py` | material *definition* dataclasses (no constitutive math) | definition only |

The Python and Rust CLT paths are line-for-line mirrors: `compute_Q` /
`compute_q`, `compute_Qbar` / `compute_qbar`, `compute_shear_Cbar` /
`compute_cbar_shear`, `_compute_ABD_matrices` / `compute_abd_matrices`,
`_compute_shear_stiffness` / `compute_shear_stiffness`. The Python side does not
implement its own isotropic or composite shell constitutive; those go through the
Rust assembler.

### 1.2 Conventions

- **Ply Voigt order** is `[1, 2, 12]`; the laminate/element order is `[x, y, xy]`.
  `tau12`/`gamma12` always denote the **engineering** shear.
- `ShellConstitutive` (`mod.rs`) carries five matrices, all in the element-local
  frame:

  ```text
  cm         3x3   membrane stiffness        N = A . eps        [N/m]
  cb_coupling 3x3  membrane-bending coupling  N = B . eps        [N]
  cb         3x3   bending stiffness         M = D . kappa      [N.m]
  cs         2x2   transverse shear stiffness Q = Cs . gamma     [N/m]
  cm_raw     3x3   stress-strain (no thickness) sigma = C_raw . eps
  ```

- The ABD system is stored as `[[A, B], [B, D]]` (`abd_matrix_flat`,
  `get_ABD_matrix`).
- Angles are in **degrees** at every entry point and converted with
  `to_radians()` / `np.radians()` immediately.

### 1.3 What "verified" means in this document

- `verified` means the equation was read in the named **held** source, or (for
  code-implemented formulas) read in the code and cross-checked between the Rust
  and Python mirrors.
- `repository citation` means the code names a work that is **not** in
  `.sources/papers/`, so the attribution cannot be checked here.
- `uncited` means no work is named anywhere for that formula.

---

## 2. Classical lamination theory and the ABD matrices

### 2.1 Code path and citations

`Laminate::new` (`laminate.rs:56`) computes ply positions, the ABD matrices and
`Cs`, in that order. The struct doc comment carries a `# References` block naming
**Jones, R.M. (1999)**, *Mechanics of Composite Materials*, 2nd ed., at
`laminate.rs:34`, and **Reddy, J.N. (2004)**, *Mechanics of Laminated Composite
Plates and Shells*, at `laminate.rs:35`. The Python class docstring names the same
two works at `laminate.py:14-15`. Neither work is held in `.sources/papers/`;
`docs/references.md` records both as `Source: repository citation`.

### 2.2 Layer positions

`compute_layer_positions` / `_compute_layer_positions`, bottom → top:

```text
z_0 = -h/2,   z_k = z_{k-1} + t_k,   h = sum_k t_k                      (code)
```

`h` is `total_thickness`. The mid-plane is `z = 0`.

### 2.3 The A, B, D definitions as implemented

`compute_abd_matrices` (`laminate.rs:95`) and `_compute_ABD_matrices`
(`laminate.py`), with `Qbar^(k)` the ply-angle-transformed reduced stiffness of
ply `k` (see §3.2), over the ply interval `[z_k, z_{k+1}]`:

```text
A_ij = sum_k Qbar_ij^(k) (z_{k+1}   - z_k  )                (laminate.rs:106)
B_ij = (1/2) sum_k Qbar_ij^(k) (z_{k+1}^2 - z_k^2)          (laminate.rs:107)
D_ij = (1/3) sum_k Qbar_ij^(k) (z_{k+1}^3 - z_k^3)          (laminate.rs:108)
```

These are the classical CLT through-thickness integrals. They are **not** printed
as numbered equations in the code; the block above is the code's own comment,
transcribed verbatim in structure and confirmed against the loop body. The
Jones/Reddy attribution sits on the struct, not on these lines, and is a
repository citation.

### 2.4 Reference-surface offset

`to_shell_constitutive_with_offset` (`laminate.rs:220-249`) shifts the reference
surface by `z0`:

```text
A_offset = A
B_offset = B - z0 A
D_offset = D - 2 z0 B + z0^2 A                                          (code)
```

This is stated in the code comment only; **no work is cited** for it. It is the
standard parallel-axis transform of CLT but the document does not claim a source
number for it.

### 2.5 Transverse shear stiffness `Cs` — the energy-equivalence method

`compute_shear_stiffness` (`laminate.rs:116-202`) and `_compute_shear_stiffness`
(`laminate.py`) first form the uncorrected shear stiffnesses

```text
A55 = sum_k Cbar55^(k) t_k,  A45 = sum_k Cbar45^(k) t_k,  A44 = sum_k Cbar44^(k) t_k
                                                                        (code)
```

For a **single ply** the result is `Cs = k * [[A55, A45], [A45, A44]]` with `k`
the user's `shear_correction_factor` (default `5/6` in the shell path, `0.75` in
the Python `Laminate` dataclass).

For **multi-ply** laminates the code instead builds an equilibrium shear-stress
shape function and energy-equivalences it:

```text
Phi_aa(z) = integral_{-h/2}^{z} Q_aa(z') z' dz'
   within ply k:  Phi_aa(z) = A_k + B_k z^2,  B_k = Q_aa^(k)/2,
                  A_k = cum_aa - (Q_aa^(k)/2) z_bot^2
   cum_aa += Q_aa^(k) (z_top^2 - z_bot^2)/2          (carried ply to ply)
I_aa^(k) = A_k^2 t_k + 2 A_k B_k dz3_k + B_k^2 dz5_k
           dz3_k = (z_top^3 - z_bot^3)/3,  dz5_k = (z_top^5 - z_bot^5)/5
flex_aa  = sum_k I_aa^(k) / C_aa^(k)
Cs_55 = D11^2 / flex_55      Cs_44 = D22^2 / flex_44
Cs_45 = sqrt(k55 k44) A45,   k_aa = Cs_aa / A_aa
                                                                        (code)
```

Fallbacks: if `flex_aa = 0` or `D_aa = 0`, `Cs_aa = k A_aa`; if `A55 = 0` or
`A44 = 0`, `Cs_45 = k A45`. The single-ply fallback and the `D_aa^2 / flex_aa`
correction are the **only** place in the material code where the shear correction
`k` is computed rather than supplied.

**Sourcing, plainly stated.** The docstring describes the intent ("energy
equivalence", "yields `k = 5/6` for homogeneous plates") but names **no paper,
author or equation**. `D_aa^2 / I_aa` and the geometric-mean rule for the
off-diagonal term are **uncited**.

### 2.6 Code → equation table

| Code | Equation | Source of the equation | Verified |
| --- | --- | --- | --- |
| `compute_abd_matrices` (A, B, D) | §2.3 | code comment + loop body | **yes**, Rust↔Python identical |
| `compute_layer_positions` | §2.2 | code | **yes** |
| `to_shell_constitutive_with_offset` | §2.4 | code | **yes** (no citation exists) |
| `compute_shear_stiffness` (multi-ply) | §2.5 | code | **yes** (uncited formula) |
| `is_symmetric` / `is_balanced` | `B ≈ 0` / `A16 = A26 ≈ 0` | code | **yes** |
| `get_equivalent_properties` (Python) | `Ex_m = 1/(h a11)`, `Ex_b = 12/(h^3 d11)` | code | **yes** (uncited) |

`is_symmetric` uses a relative tolerance `1e-10 * |A|_max`; `is_balanced` the
same. `get_equivalent_properties` uses `a = A^{-1}`, `d = D^{-1}` and the
`12/h^3` bending factor, with no citation.

---

## 3. The orthotropic plane-stress constitutive law

### 3.1 `Q` in the principal material axes

`OrthotropicMaterial::compute_q` (`orthotropic.rs:43-64`) and `compute_Q`
(`laminate.py`), plane stress (`sigma3 = 0`):

```text
nu21 = nu12 E2 / E1
denom = 1 - nu12 nu21
Q11 = E1 / denom,  Q22 = E2 / denom,  Q12 = nu12 E2 / denom,  Q66 = G12   (code)
```

`E3`, `G23`, `G13`, `nu23`, `nu31` do not enter `Q`; they enter only the shear
law (§3.3). A warning is printed if `|denom| < 1e-6`.

### 3.2 Rotation to the lamina / element frame — `Qbar`

`compute_qbar` (`orthotropic.rs:70-100`) and `compute_Qbar` (`laminate.py`), with
`c = cos(theta)`, `s = sin(theta)`:

```text
Qbar11 = Q11 c^4 + 2 (Q12 + 2 Q66) s^2 c^2 + Q22 s^4
Qbar22 = Q11 s^4 + 2 (Q12 + 2 Q66) s^2 c^2 + Q22 c^4
Qbar12 = (Q11 + Q22 - 4 Q66) s^2 c^2 + Q12 (s^4 + c^4)
Qbar16 = (Q11 - Q12 - 2 Q66) s c^3 + (Q12 - Q22 + 2 Q66) s^3 c
Qbar26 = (Q11 - Q12 - 2 Q66) s^3 c + (Q12 - Q22 + 2 Q66) s c^3
Qbar66 = (Q11 + Q22 - 2 Q12 - 2 Q66) s^2 c^2 + Q66 (s^4 + c^4)          (code)
```

The only attribution is the inline comment `etc. (Jones, 1999 — eq. 2.78)` at
`orthotropic.rs:69`. `docs/references.md` lists Jones as a repository citation
and **no held copy exists**, so the equation number `2.78` is a repository
assertion and could not be checked. The formulas themselves were read from the
code and are identical in Rust and Python.

**What `Qbar` rotates between.** `theta_deg` is the **ply** angle from the
laminate x-y axes. The rotation of a laminate into the **element-local** frame is
*not* done here: it is done in the element kernels by
`rotate_constitutive_to_local` (`mitc3.rs:165`), which builds the Reuter
transformation for the angle `beta = angle(e1, x_ref)` in the shell plane and
applies `A_local = T(-beta) A_global T(-beta)^T` (and likewise `B`, `D`) plus a
2×2 rotation of `Cs`. **MITC4 does not perform this rotation**:
`compute_ke_local` uses `pre.constitutive` directly and `mitc4.rs` contains no
equivalent call. For an anisotropic material on a surface whose element `e1`
differs from the global x-y axes this is an unreconciled frame mismatch; it is
recorded here as a finding, not fixed. Element-level detail belongs in
`docs/formulations/shell-elements.md`.

### 3.3 Transverse shear transformation — `Cbar_shear`

`compute_cbar_shear` (`orthotropic.rs:109-125`) and `compute_shear_Cbar`
(`laminate.py`):

```text
C44 = G23 c^2 + G13 s^2
C55 = G13 c^2 + G23 s^2
C45 = (G13 - G23) s c
returns [[C55, C45], [C45, C44]]                                        (code)
```

**Uncited** in both files.

### 3.4 The single-ply shell constitutive

`impl Material for OrthotropicMaterial::constitutive` (`orthotropic.rs:129-153`),
with `h` the thickness, `h3_12 = h^3/12`, and `shear_correction` the shell
correction factor:

```text
qbar   = compute_qbar(0.0)              # principal axes, angle 0
cbar_s = compute_cbar_shear(0.0)
cm          = qbar * h
cb          = qbar * h3_12
cb_coupling = 0                          # no membrane-bending coupling, single ply
cs          = cbar_s * (shear_correction * h)
cm_raw      = qbar
                                                                        (code)
```

### 3.5 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `compute_q` | §3.1 | repository (implied by the CLT refs) | **yes**, Rust↔Python |
| `compute_qbar` | §3.2 | `Jones 1999 eq. 2.78` (`orthotropic.rs:69`) | equation **yes**; attribution **repository citation** |
| `compute_cbar_shear` | §3.3 | **uncited** | **yes** |
| `OrthotropicMaterial::constitutive` | §3.4 | none | **yes** |

---

## 4. The composite shell constitutive

### 4.1 `composite_constitutive` — the ABD → shell mapping

`composite.rs:23-58` reshapes flat arrays into `ShellConstitutive`:

```text
cm          = A   (a_flat, row-major 3x3)
cb_coupling = B   (b_flat)
cb          = D   (d_flat)
cs          = Cs  (cs_flat, row-major 2x2)
cm_raw      = A / h          (h = thickness argument, 0 if |h| < 1e-30)  (code)
```

Note the naming: **`cm` is `A`, `cb` is `D`** — the field names do not follow the
letter each one stores. The docstring states the mapping explicitly
(`composite.rs:9-16`), and it is the only place the ABD system becomes what the
element kernels consume.

### 4.2 `Laminate::to_shell_constitutive*` — the same mapping with the offset

`to_shell_constitutive` / `to_shell_constitutive_with_offset` (`laminate.rs:204-249`)
do the same with the offset ABD of §2.4 and the energy-equivalent `Cs` of §2.5:

```text
cm          = A
cb_coupling = B - z0 A
cb          = D - 2 z0 B + z0^2 A
cs          = Cs                (from compute_shear_stiffness)
cm_raw      = A / h             (0 if |h| < 1e-30)                     (code)
```

### 4.3 Code → equation table

| Code | Output field | Input | Verified |
| --- | --- | --- | --- |
| `composite_constitutive` | `cm, cb_coupling, cb, cs, cm_raw` | `A, B, D, Cs, h` | **yes**, matches docstring and caller layout |
| `Laminate::to_shell_constitutive_with_offset` | same | offset ABD + `Cs` | **yes** |

No work is cited for the mapping itself; it is a data-layout convention, not a
formulation.

---

## 5. Failure criteria

Both `failure.rs` and `failure.py` carry a `# References` / `References` block
naming **Tsai & Wu (1971)** and **Hashin (1980)**. The entries in
`docs/references.md` mark both as repository citations — the authors, journal,
volume and pages come from the published record, **no copy is held**, and the
article titles are not re-verified. The criteria below were read from the code;
their attribution is therefore self-attested, not independently verified.

### 5.1 Tsai-Wu

`tsai_wu` (`failure.rs:91`) and `tsai_wu_failure_index` (`failure.py:76-137`):

```text
F1 = 1/Xt - 1/Xc          F2 = 1/Yt - 1/Yc
F11 = 1/(Xt Xc)           F22 = 1/(Yt Yc)          F66 = 1/S12^2
F12 = f12_coefficient / sqrt(Xt Xc Yt Yc)          (default f12_coefficient = -0.5)

FI = F1 s1 + F2 s2 + F11 s1^2 + F22 s2^2 + F66 t12^2 + 2 F12 s1 s2       (code)
```

Failure when `FI >= 1`; the reported mode is `Combined` (`COMBINED`) when it
fails, `NoFailure` otherwise. `reserve_factor = 1/sqrt(FI)`,
`margin_of_safety = reserve_factor - 1` (from `FailureResult`). The default
`f12_coefficient = -0.5` makes `F12 = -1/(2 sqrt(Xt Xc Yt Yc))`, the classical
`F12* = -1/2 sqrt(F11 F22)` normalisation; the code does not print that identity
or cite it.

### 5.2 Hashin (mode-specific)

`hashin` (`failure.rs:122`) and `hashin_failure_indices` (`failure.py:150-243`).
`ST` is taken as `s23`, and `s23` defaults to `s12/2` when not supplied
(`failure.rs:34`, `laminate.py` `StrengthProperties.__post_init__`):

```text
fiber tension      (s1 >= 0):  FI = (s1/Xt)^2 + (t12/S12)^2
fiber compression  (s1 <  0):  FI = (|s1|/Xc)^2
matrix tension     (s2 >= 0):  FI = (s2/Yt)^2 + (t12/S12)^2
matrix compression (s2 <  0):  FI = (s2/(2 ST))^2
                                   + ((Yc/(2 ST))^2 - 1) (s2/Yc)
                                   + (t12/S12)^2                        (code)
```

**A docstring/code discrepancy on fibre compression, recorded as found.**
`failure.py` describes the mode in its `Notes` as "`|sigma1/Xc| = 1`" (linear),
but the implementation computes `(abs(s1)/Xc)**2` in both Python and Rust. The
classical Hashin fibre-compression index is of the form
`(sigma1/Xc)^2` *with* a shear term, or a largely linear form depending on the
variant; without a copy of Hashin (1980) the code's simplification cannot be
attributed or checked. What is verifiable is the mismatch: the docstring's stated
formula and the code's formula differ.

### 5.3 Maximum stress

`max_stress` (`failure.rs:173`) and `max_stress_failure_index`
(`failure.py:245-310`):

```text
R1  = s1/Xt if s1 >= 0 else |s1|/Xc
R2  = s2/Yt if s2 >= 0 else |s2|/Yc
R12 = |t12|/S12
FI  = max(R1, R2, R12);   mode = argmax                                   (code)
```

**Uncited.** Both `# References` blocks name only Tsai-Wu and Hashin, so the
maximum-stress criterion is implemented with no source at all.

### 5.4 Stress and strain transformations

`stress_transformation` / `stress_transformation_matrix` and
`strain_transformation` / `strain_transformation_matrix`, with `c = cos(theta)`,
`s = sin(theta)`, `sigma = [sigma_xx, sigma_yy, tau_xy]`,
`eps = [eps_xx, eps_yy, gamma_xy]`:

```text
T   = [ c^2   s^2   2 s c ]        such that sigma_12  = T   sigma_xy
      [ s^2   c^2  -2 s c ]
      [-s c   s c   c^2-s^2 ]

T_e = [ c^2   s^2    s c ]         such that eps_12    = T_e eps_xy
      [ s^2   c^2   -s c ]
      [-2 s c 2 s c  c^2-s^2 ]                                          (code)
```

**Uncited** in both files. The two matrices were confirmed identical between Rust
and Python.

### 5.5 Code → equation table

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `tsai_wu` / `tsai_wu_failure_index` | §5.1 | Tsai & Wu 1971 (self-attested) | equation **yes** (code); attribution **repository citation** |
| `hashin` / `hashin_failure_indices` | §5.2 | Hashin 1980 (self-attested) | equation **yes** (code); attribution **repository citation**; fibre-compression docstring mismatch |
| `max_stress` / `max_stress_failure_index` | §5.3 | **uncited** | **yes** |
| `stress_transformation` / `strain_transformation` | §5.4 | **uncited** | **yes** |
| `s23 = s12/2` default | - | **uncited** | **yes** |
| `reserve_factor = 1/sqrt(FI)` | - | **uncited** | **yes** |

---

## 6. The isotropic law

`impl Material for IsotropicMaterial::constitutive` (`isotropic.rs:20-68`), with
`h` the thickness and `shear_correction` the shell correction:

```text
factor_m   = E h / (1 - nu^2)                # plane-stress equivalent modulus * h
factor_b   = E h^3 / (12 (1 - nu^2))
G          = E / (2 (1 + nu))
factor_raw = E / (1 - nu^2)                  # plane-stress equivalent modulus
half_1_minus_nu = (1 - nu)/2

cm = factor_m   * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
cb = factor_b   * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
cs = shear_correction * G * h * I2
cm_raw = factor_raw * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]       (code)
cb_coupling = 0
```

The **plane-stress equivalent modulus the shell path uses** is
`E/(1 - nu^2)`: `cm = factor_m` integrates it through the thickness and
`cm_raw = factor_raw` exposes it without thickness. The in-plane shear entry
`factor_m (1-nu)/2` equals `G h`, since `E (1-nu) / (2 (1-nu^2)) = E/(2(1+nu))`.
No work is cited; this is textbook plane-stress elasticity and the code names no
source.

| Code | Equation | Attribution | Verified |
| --- | --- | --- | --- |
| `IsotropicMaterial::constitutive` | §6 | **uncited** | **yes** |

---

## 7. Uncited formulas, gaps and mismatches

This is the list the revision deliberately records rather than fills in.

1. **Transverse shear correction (`Cs`) has no source.** The energy-equivalence
   construction of §2.5 — including `Cs_aa = D_aa^2 / I_aa` and the geometric-mean
   off-diagonal `sqrt(k55 k44) A45` — names no paper, author or equation. This is
   the material code's equivalent of the `0.15` drilling factor in the element
   code.
2. **`Cbar_shear` (§3.3) is uncited.**
3. **The reference-surface offset transform (§2.4) is uncited.**
4. **The maximum-stress criterion (§5.3) is uncited.** The `# References` blocks
   of both `failure.rs` and `failure.py` name only Tsai-Wu and Hashin.
5. **The stress/strain transformations and the `S23 = S12/2` default (§5.4) are
   uncited.**
6. **Tsai-Wu and Hashin are repository citations.** Both are self-attested in the
   code, neither paper is held, and `docs/references.md` marks both DOIs "to
   verify". The criteria as implemented are therefore verified only against the
   code, not against the sources the code names.
7. **Hashin fibre compression is quadratic in the code, linear in the docstring.**
   `failure.py` states `|sigma1/Xc| = 1`; the code computes `(|sigma1|/Xc)^2`.
   Which one is Hashin (1980) cannot be decided without the paper.
8. **MITC4 does not rotate the constitutive into the element frame.** Only MITC3
   calls `rotate_constitutive_to_local`. For anisotropic materials this is a
   frame mismatch in `mitc4.rs`; recorded, not fixed.
9. **The CLT citations carry no equation-level reference.** Jones (1999) and
   Reddy (2004) are named on the struct/class but not on `compute_abd_matrices`,
   `compute_q` or `compute_qbar` (the sole exception is the `eq. 2.78` comment at
   `orthotropic.rs:69`).

No `0.15`-style magic constant exists in the material code; the uncited formulas
above are the material-side gaps.

---

## References

`docs/references.md` is the canonical bibliography. Where this document and that
file disagree, `docs/references.md` should be corrected; this document does not
duplicate its per-entry verification annotations.

1. Jones, R.M., *Mechanics of Composite Materials*, 2nd ed., Taylor & Francis,
   1999. DOI: to verify. **Not held.** Cited by the code at
   `crates/aeroelast-core/src/materials/laminate.rs:34`,
   `src/aeroelast/core/laminate.py:14` and (for `Qbar`)
   `crates/aeroelast-core/src/materials/orthotropic.rs:69` (`eq. 2.78`).
2. Reddy, J.N., *Mechanics of Laminated Composite Plates and Shells: Theory and
   Analysis*, 2nd ed., CRC Press, 2004. DOI: to verify. **Not held.** Cited by
   the code at `crates/aeroelast-core/src/materials/laminate.rs:35` and
   `src/aeroelast/core/laminate.py:15`.
3. Tsai, S.W., Wu, E.M., "A general theory of strength for anisotropic
   materials", *Journal of Composite Materials*, 5(1):58–80, 1971. DOI: to
   verify. **Not held.** Self-attested by the code at
   `crates/aeroelast-core/src/materials/failure.rs:9` and
   `src/aeroelast/constitutive/failure.py:12`.
4. Hashin, Z., "Failure criteria for unidirectional fiber composites",
   *Journal of Applied Mechanics*, 47(2):329–334, 1980. DOI: to verify. **Not
   held.** Self-attested by the code at
   `crates/aeroelast-core/src/materials/failure.rs:10` and
   `src/aeroelast/constitutive/failure.py:14`.

No PDF in `.sources/papers/` corresponds to any of these four works; all four
attributions are repository citations.
