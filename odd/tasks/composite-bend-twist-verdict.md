# Task document — composite shell bend-twist: settle the reference, then give the element a verdict (issue #9)

- **Opened**: 2026-10-01
- **Phase**: OBSERVE / RESEARCH (no production code touched yet)
- **Trigger**: GitHub issue #9 (`[bug] bend-twist: MITC4 matches CalculiX but is 3.4-4.8x the
  CLT coupon — the coupon reference needs settling`) plus its two comments of 2026-10-01.
- **Base**: `main @ c465aac`, working tree clean, 9 commits ahead of `origin/main`.
- **Environment**: conda `aeroelast-dev` (`~/miniconda3/envs/aeroelast-dev/bin/python`),
  `ccx` on `PATH`, Rust extension `_aeroelast` built. Rust workspace root is `crates/`.

---

## 1. Objective

Decide, with reproducible evidence, whether the composite MITC4/MITC3 shell path
(`crates/aeroelast-core/src/elements/`, constitutive from `materials/composite.rs`) mispredicts
**bend-twist coupling** (`D16`/`D26`), and then:

- if the element is wrong → localise, fix, and add the guard test that would have caught it;
- if the element is right → land the residual as a **documented validity limit** and replace the
  wrong coupon reference, so the FSI/aeroelastic numbers stop inheriting an unproved bias.

The deliverable is a **verdict with evidence**, not a green test. The acceptance bar is that a
third party can re-run the coupon, the reference, and the convergence study and reach the same
conclusion.

## 1b. VERDICT (2026-10-01, after WU1-WU4)

The composite MITC4 path is **not** the defect. The 3.4-4.8x is a **units/convention error in the
reference formula**: the issue's `theta = (D^-1)[2,0] * P L^2 / 2` omits

* the moment **per unit width**, `Mx = P(L-x)/b`, and
* the **geometric** twist rate, `theta' = kappa_xy / 2` (the CLT `kappa_xy` is the *engineering*
  twist curvature, twice the geometric one),

which together are a factor `1/(2b) = 5` for the issue's own coupon (`b = 0.1 m`). With both
factors in place the element, a converged Rayleigh-Ritz solution of the same plate BVP, and a
hand-authored CalculiX S8R composite deck agree within ~1-3%:

| route | tip twist, `[45,0,0,45]s` coupon [deg] | vs element |
| --- | ---: | ---: |
| MITC4 composite, mesh 16x96 | -0.282051 | 1.0000 |
| Rayleigh-Ritz, converged (independent discretisation, same BVP) | -0.291250 | 0.9709 |
| CalculiX S8R, hand-authored deck (independent of our writer) | -0.281413 | 1.0023 |
| free-edge CLT reference, correctly normalised | -0.289944 | 0.9728 |
| **the issue's formula** | -0.057990 | **4.864** |

Full evidence in §12. Hypotheses H2 (metric), H3 (mesh) and H5 (element defect) are **refuted**;
H1 (reference) is **confirmed**; H4 (shared laminate) is guarded by a test; H6 (3D section
distortion) is still **open** and belongs to WU5/WU7, not to this verdict.

**This verdict does not clear the blade-level claim.** A coupon cannot explain a 9x blade twist
gap (+8.9 deg vs the BeamDyn anchor +0.98 deg); the coupon only removes the *element's D16/D26
coupling* as a candidate cause, and it fixes the reference that was making a correct element look
wrong. WU5/WU7 remain.

## 2. What the issue claims (and why it matters)

| Claim | Source |
| --- | --- |
| Coupon `[45,0,0,45]s`, L=1 m, b=0.1 m, t=8 mm, tip load, mesh `NX=2, NY=12`: shell twist = **-0.44132 deg**, CLT = -0.09257 deg, ratio **4.767** | issue body §1 |
| Ratio is layup dependent (3.36-4.77) and exactly 0 when `D16 = 0` | issue body §2 |
| The CLT reference is internally ambiguous (`(D^-1)[2,0]` vs `[2,1]`), and may not be the right closed form | issue body §3 |
| The same strip matches CalculiX S8R to 0.4%, but that deck is written by **our** writer | issue body §4 |
| Blade level: shell tip twist +8.9 deg vs BeamDyn +0.98 deg; FSI de-loads -25%/-35% vs literature -8%/-13% | issue body §Impact, comment 2 |
| The excess scales with aerodynamic load (`CT~0.5` -> -35%; feathered `CT~0` -> none) | issue comment 2 |

The blade/FSI claim is the *motivation*: every normal-operation aeroelastic number inherits the
elastic twist fed back into the BEM. It is **not** evidence about the element, and it is not
decidable at coupon resolution.

## 3. Census: what the issue cites vs what this tree contains (measured 2026-10-01)

### 3.1 Cited artifacts

| Artifact the issue cites | In tree | In git history | Method |
| --- | --- | --- | --- |
| `tests/test_laminate_bend_twist.py` (and its `sweep_laminate_bend_twist.py`) | **no** | **never** | `ls tests/`, `git log --all -S"_analytic_tip_twist_deg"`, `git log --all --name-only --diff-filter=A \| grep -i twist`, `git ls-tree` over every local+remote branch, `git fsck --unreachable` blob scan |
| `tests/test_composite_ccx_parity.py` | **no** | **never** | same |
| `tests/test_laminate_shear_coupling.py` | **no** | **never** | same |
| `tests/test_box_torsion_bending_benchmark.py` | **no** | **never** | same |
| `docs/twist_distortion_diagnosis.md` ("open since 2026-09-16") | **no** | **never** | `ls docs/` |
| FSI campaign numbers (`16.51 MW`, `2.541 MN`, `12.34 MW`, `-25%/-35%`, `Zhou 2025`) | **no** | **no** | repo-wide grep over `*.md/*.py/*.json/*.yaml`; `examples/fsi/` holds only the preCICE perpendicular-flap tutorial |

**Consequence.** The issue's reproduction cannot be executed against this tree. Either the
reproduction lives in a checkout that is not this one, or it does not exist. WU1 therefore
*rebuilds* the coupon from the issue's description; it does not "fix" an existing test.

### 3.2 The issue's CLT column does not reproduce with this repo's CLT

`[45,0,0,45]s`-family, symmetric expansion, uniformly 1 mm plies, h = 8 mm:

| layup (as named) | issue `D16` | repo `Laminate` `D16` (CFRP 120/10/5, nu=0.3) | independent textbook Qbar CLT |
| --- | ---: | ---: | ---: |
| `[45,0,0,45]s` | 264.52 | 701.93 | 701.93 |
| `[45,-45,-45,45]s` | 226.73 | 221.66 | 221.66 |
| `[0,45,45,0]s` | 37.79 | 480.27 | 480.27 |
| `[45,45,45,45]s` | 302.31 | 1182.20 | 1182.20 |
| `[30,0,0,30]s` | 394.89 | 895.67 | 895.67 |
| `[60,0,0,60]s` | 63.27 | 320.11 | 320.11 |
| `[0,0,0,0]` | 0.00 | 0.00 | 0.00 |

The mismatch is **not** a material or thickness rescale (issue/repo ratios: 0.38, 1.02, 0.08,
0.26, 0.44, 0.20). No single scale factor maps one column onto the other.

## 4. Findings of the first analysis pass

### F1 — The repo's CLT is correct (this is not where the bug is)

`aeroelast.core.laminate.Laminate.get_ABD_matrix()` was cross-checked against a from-scratch
Qbar/T implementation written for this pass (standard explicit `Qbar16 = (Q11-Q12-2Q66)c^3 s - (Q22-Q12-2Q66)c s^3`).
Agreement on `D16`, `D26`, `D11`, `D22`, `D12` is **exact to the last printed digit** for all
seven layups above. The ABD that feeds the element is trustworthy; any coupon verdict must be
settled on the *reference formula*, the *metric*, or the *assembled* strip — not on our ABD.

### F2 — The issue's `D16` pattern is not a physical laminate with the stacking the issue names

`D16` of a uniform-ply stack is `(1/3) * sum_j Qbar16(theta_j) * (z_j^3 - z_{j-1}^3)`; the
*z^3 weights* are pure geometry, and the *material* enters only through three combinations
(`Q11-Q22` for the +/-45 family, `X = Q11-Q12-2Q66`, `Y = Q22-Q12-2Q66` for the 30/60 family).

Fitting the issue's six `D16` values with a single orthotropic material inside a **physical**
box (E1 in [30,400] GPa, E2 in [1,40] GPa, G12 in [1,25] GPa, nu12 in [0,0.5]):

| stacking hypothesis | best-fit max |error| | verdict |
| --- | ---: | --- |
| **8-ply symmetric** (`[45,0,0,45]s`, what the issue names) | **88.0%** | **refuted** — no physical material |
| 16-ply symmetric | 88.0% | refuted (same weights as 8-ply up to scale) |
| **4-ply literal** (`[45,0,0,45]`, h = 4 mm) | **12.7%** | closest, but needs E1=257 GPa, E2=40 GPa, G12=1 GPa, nu12=0 — not a real CFRP |

The four +/-45/0/90 layups are *exactly proportional* in the 4-ply-literal reading (constant
ratio 2.046), which is what a single-material stack must do; the 8-ply reading cannot be made
proportional at all. **So the issue's reference column was very probably computed for a
different laminate than the one its prose names** (a 4-ply stack, or a stack built by a
different convention), with an anisotropic material that is not the repo's CFRP.

### F3 — A self-consistent sweep must share **one** laminate object between the solve and the reference

The cited sweep takes the reference `D` from `t._laminate()` and the shell displacement from
`t._solve_strip(MESH, P)`. If those two paths do not build *exactly* the same ABD (different
ply count, different total thickness, different angle expansion, different ply thickness), the
ratio is measuring the *difference between two laminates*, not an element error. F2 shows the
reference side is already suspect, so WU1 must make the laminate a single immutable object
consumed by both sides — and assert that the ABD the element received equals the ABD the
reference used.

### F4 — The element's `D16` chain is already exercised locally

| Surface | Where | What it proves |
| --- | --- | --- |
| `kappa = [kxx, kyy, 2 kxy]` convention, third covariant row doubled for engineering shear | `mitc4.rs:995-1001`, `mitc4.rs:7042-7052` | the twist row is the *engineering* twist, not the tensor twist |
| `cb` = the CLT `D` matrix, all 9 entries (incl. `D16`,`D26`) | `materials/composite.rs:23-70`, `materials/mod.rs:20-35` | the constitutive carries the full `D` |
| through-thickness blocks `W_01=B`, `W_11=D`, `W_02=D`, `W_22=cm h^4/144` | `mitc4.rs:1253-1315` (`resultant_moment_matrix`) | correct pairing of B and D with the `epsilon_m + z kappa + z^2 E2` decomposition |
| patch test T1a applies constant curvature states, including `("kappa_xy", [0,0,1e-3])`, and checks `M = cb kappa` | `mitc4.rs:7058-7090` | a rectangular patch maps each `D` entry to the right moment |

Composite elements differ from isotropic ones **only** by the `ShellConstitutive` they are
handed (`assembler.rs:1019-1071`, `aeroelast-py/src/elements.rs:628,742`); they share the
kernel. A 4.7x bend-twist error is therefore **not** consistent with a local constitutive
defect — the surviving loci are the *assembled strip* (BC/load path, membrane-bending
interaction on a coarse mesh), the *metric*, or the *reference*.

### F5 — A one-line "CLT closed form" for the strip is not available (analytic)

The cited reference evaluates `theta_tip = (D^-1)[2,0] * P L^2 / 2`. Three independent problems:

1. **Dimensional**: `PL^2/2` is a moment; the plate relation `kappa = D^-1 M` needs the moment
   **per unit width**, `M(x) = P(L-x)/b`. A missing `1/b` is a factor 10 at `b = 0.1 m`. The
   observed 4.77 is *not* 10, so `1/b` alone does not explain it — but the dimensional check
   must be done explicitly, because a total-force load and a per-width load differ by exactly
   that factor in both the reference and the FEM.
2. **Index ambiguity**: `(D^-1)[2,0]` is `kxy` per `Mx`; `(D^-1)[2,1]` is `kxy` per `My`. Taking
   `My = 0` is a *model choice* (free edges in a beam-like reduction). The two columns disagree
   by layup — and the issue itself shows `[30,0,0,30]` flipping sign between them. Settling the
   index is settling the *boundary-value problem*, not a typo.
3. **The model is over-constrained**. Assume a rigid cross-section, `w(x,y) = W(x) + y*phi(x)`.
   Then `kxy = 2 phi'` is uniform in `y`, and the free-edge conditions `My(y=±b/2) = 0`,
   `Mxy(y=±b/2) = 0` reduce to two *homogeneous* equations in `(W'', 2 phi')`:
   `D12 W'' + D26 2phi' = 0`, `D16 W'' + D66 2phi' = 0`. A non-trivial solution requires
   `det = D12 D66 - D26 D16 = 0`, which is false in general (it is false for every layup in the
   issue's table). The same failure appears in the membrane: with a rigid cross-section
   `eps_yy = 0`, so `N_yy = A12 eps_xx != 0` at a free edge. **A free-edged strip must distort
   (anticlastic and warping); a rigid-cross-section beam formula cannot be its reference.**
   This is precisely why the issue says the formula "may not be the plate-bending relation at
   all", and it means the correct reference has to be *computed*, not quoted.

### F6 — The CCX counter-evidence is not independent

The issue's own caveat is correct and this repo's memory confirms it: the CCX deck is written
by our writer, and a previously localised defect was in that writer's `*SHELL SECTION,
COMPOSITE` block (per-ply orientation realisation: `[0/90]2s` -34.5%, mixed -16.6% against an
independent CLT, while a smeared `MATERIAL=` section matched CLT to 0.6%). "MITC4 matches CCX"
is therefore not evidence while both sides route through the same writer.

## 5. Hypotheses (each with its falsification test)

| id | hypothesis | falsified by |
| --- | --- | --- |
| **H1** | the *reference* is wrong (index, dimensional factor, or a rigid-section beam formula applied to a free-edged strip) | WU4: a converged independent solution of the *same* plate BVP; if the shell converges to it, H1 holds and the element is exonerated |
| **H2** | the *metric* is wrong (LSQ tip twist aliases local tip deformation, clamped-edge artefacts) | WU3: metric sweep (LSQ vs edge slope vs end rotation vs Saint-Venant) on the *same* solution |
| **H3** | the reported ratio is a **coarse-mesh artifact** (`NX=2`, 2 elements along the length) | WU2: mesh convergence NX x NY = (2,12) -> (4,24) -> (8,48) -> (16,96); the ratio must move |
| **H4** | the *laminate identity* is not shared between solve and reference (F2/F3) | WU1: single immutable ABD object; assert element ABD == reference ABD |
| **H5** | a real defect in the assembled composite path (load/BC realisation, mixed membrane-bending assembly, distorted-element behaviour) | WU2/WU3 stable at 4.7x + WU4 reproduced by an independent FEM + the patch tests of F4 still passing |
| **H6** | CLT plate theory itself is insufficient for this strip (3D section distortion, Brazier-type ovalization) | WU5: 3D solid layer-wise model. Bounds what a coupon can say about the blade; does **not** excuse a 4.7x CLT-vs-CLT gap |

**Status after WU1-WU4**: H1 confirmed (the reference is wrong by `1/(2b)`); H2 refuted (the three
twist metrics agree to 0.0% on a refined mesh); H3 refuted (the ratio is flat under refinement in
both directions); H4 converted into a test that asserts the element's ABD equals the reference's;
H5 refuted (two independent routes reproduce the element to 1-3%); H6 **refuted at coupon scale**
by WU5 (3D layer-wise elasticity reproduces the shell to 0.2% - see §13), which leaves the
blade-level gap to WU7 as a blade-model question, not a theory-of-the-coupon question.

## 6. Work units

### WU0 — Census and conventions (DONE, this document)
Artifact census (3.1), CLT cross-check (F1), reference forensics (F2/F3), element audit (F4),
analytic objection to the closed form (F5), independence audit (F6).

### WU1 — Rebuild the coupon deterministically in-repo
New `tests/test_laminate_bend_twist.py` (the file the issue names, authored from its
description), containing: strip mesh builder, **one** immutable laminate object, a tip-load
application whose resultant is mesh independent, an explicit twist metric with its definition
in the docstring, and assertions that (a) the element received the same ABD as the reference,
(b) the `D16 = 0` control gives ~0 twist. No assertion on the absolute twist yet — WU1 produces
the *measurement*.

### WU2 — Mesh convergence of the coupon
Sweep NX x NY and (for the winning metric) report shell twist vs (NX,NY). If the ratio drifts
with refinement, the issue's number is a coarse-mesh artifact (H3) and the fix is the test.

### WU3 — Metric sensitivity
Compute the tip twist four ways on the same displacement field: LSQ plane fit of tip nodes,
edge slope at the mid-line, rotation of the end cross-section, and the Saint-Venant torsional
rotation. Report the spread. A metric-dependent 4.7x is a test defect, not an element defect.

### WU4 — Settle the reference (the core of the issue)
Establish the correct reference for the *free-edged anisotropic strip* by **two independent
routes**, and require them to agree:

1. **Converged Rayleigh-Ritz** on the same CLT/Mindlin plate BVP: expand
   `w(x,y) = sum_ij a_ij X_i(x) Y_j(y)` with `X_i` clamped-free beam functions and `Y_j`
   polynomial/Chebyshev, minimise the CLT energy, and show convergence by increasing
   `(i,j)` until the tip twist changes by < 0.1%. This is an *upper bound on the compliance*
   and is independent of the FE discretisation under test.
2. **Independent FE / independent deck** (choose one in §9): a CalculiX deck authored by hand
   (not by `write_ccx_mesh`), validated by a same-model A/B (composite section card vs
   isotropic-equivalent `MATERIAL=` card) on a **quasi-isotropic** coupon where both must
   agree; or a 3D solid layer-wise model.

Acceptance: the element's converged coupon twist agrees with the settled reference within the
usual suite tolerance (5%) **or** the element is declared wrong with the discrepancy localised.

### WU5 — Independent judge for the strip physics (3D)
3D solid layer-wise (CalculiX `C3D20R`, one element per ply, independently authored geometry)
of the coupon and of an outboard blade segment. Answers H6 and bounds the validity limit.

### WU6 — Verdict and landing
Either (a) the Rust fix + a guard test that fails on the old code, with a before/after coupon
table, or (b) a documented validity limit + the corrected reference test replacing the wrong
one. Either way: `docs/validation-matrix.md` row, and the FSI/aeroelastic numbers either
re-derived or explicitly marked as inheriting an unproved twist bias.

### WU7 — Blade-level follow-up (only if WU1-WU6 conclude the element is right)
The blade twist gap (+8.9 vs +0.98 deg) is a *section/theory* question (shear deformation,
warping, 3D effects), not a coupon question. Scope it separately rather than letting it
contaminate the element verdict.

## 7. Acceptance criteria

- [ ] The coupon, its laminate, its load and its metric are defined by one executable artifact
      in this tree, and the element's ABD is asserted equal to the reference's ABD.
- [ ] Mesh convergence and metric sensitivity are reported as tables, not asserted away.
- [ ] A reference independent of both the element and our writer reproduces (or refutes) the
      coupon twist, with its own convergence evidence.
- [ ] The verdict is one of: **element wrong (fix + guard)**, or **element right (documented
      validity limit + corrected reference)** — never "the test is xfail".
- [ ] No citation of the missing artifacts remains; the issue is updated with the measured
      tables and the residual is stated in the validity envelope.

## 8. Verification commands

```bash
export PATH=$HOME/miniconda3/envs/aeroelast-dev/bin:$PATH
cd /home/efirvida/Desktop/dev/fem-shell

# coupon + metric + convergence (WU1-WU3)
python -m pytest -q -o addopts="" tests/test_laminate_bend_twist.py

# the element's own invariants must stay green (F4)
cargo test --manifest-path crates/Cargo.toml -p aeroelast-core -q

# existing composite coverage must not regress
python -m pytest -q -o addopts="" tests/test_material_suite.py tests/test_composite_b_coupling.py \
    tests/test_composite_layup_parity.py tests/test_rust_composite.py

# independent judge (WU4/WU5), written by hand, run through CCX
ccx -i <hand-authored-deck>.inp
```

## 9. Decisions needed (user)

| id | decision | options |
| --- | --- | --- |
| **D1** | Does issue #9's reproduction exist in some checkout that is not this one? | (a) no — rebuild it here from the issue text; (b) yes — give the path/repo and I start from it; (c) the target is not #9 — say what `elementos` means |
| **D2** | Which independent judge to invest in | (a) Rayleigh-Ritz in-repo (cheap, same theory); (b) hand-authored CCX deck + quasi-iso A/B (medium); (c) 3D layer-wise solid (expensive, decisive for H6); (d) chain a -> b -> c, stopping as soon as two agree |
| **D3** | Does the blade-level FSI claim block this work? | (a) coupon verdict first, blade scoped separately (recommended); (b) blade twist is the real deliverable and the coupon is a means |

## 10. Risks, cost, out of scope

- **Risk R1 — wrong target.** If the issue's reproduction lives elsewhere, WU1 rebuilds a
  coupon that may not match its mesh/load/metric. Mitigation: WU1 states its definitions
  explicitly and D1 resolves the ambiguity before any verdict.
- **Risk R2 — reference chasing.** F5 says no short formula exists; a "reference" that is not
  converged will produce another 4.7x-style artefact. Mitigation: two independent routes plus
  convergence evidence, not a single quoted number.
- **Risk R3 — theory confound.** H6 can make the *blade* claim true while the *coupon* verdict
  is "element right". Mitigation: keep WU6 and WU7 separate.
- **Cost.** WU1-WU3: minutes. WU4 route (a): hours. WU4 route (b) / WU5: a session each
  (meshing + CCX runs). WU5 blade segment: a session.
- **Out of scope**: isotropic shells (validated 1-2% vs CCX), the BEM/aero side (rigid numbers
  match the literature), and the FSI coupling numerics (the parked V50 divergence is a separate
  instability symptom).

## 11. Progress log

- **2026-10-01, pass 1 (OBSERVE)**: artifact census; CLT cross-check (F1); reference forensics
  (F2/F3: the named 8-ply stacking is refuted at 88%, the 4-ply literal reading fits to 13%
  with an unphysical material); element audit (F4); analytic objection to the closed form (F5);
  independence audit (F6). Baseline suites re-run and green (mesh 51 passed; material composite
  subset 12 passed; b-coupling 4 passed). Body of the issue's second comment recorded (blade
  level corroboration, not verifiable in-tree).
- **2026-10-01, pass 2 (WU1-WU4, DONE)**: `tests/test_laminate_bend_twist.py` authored (the file
  the issue named, which never existed here) with a single-laminate design rule, three twist
  metrics, a mesh-convergence study, a constant-moment constitutive probe, a converged
  Rayleigh-Ritz reference and a hand-authored CCX S8R deck. 11 tests, all green. Verdict in §1b;
  evidence in §12. Committed as `939e49a` after RDD review `review-4af7f1c0f12b1d65` (4/4 lenses,
  approved, authority burned) with one bounded correction for the CRITICAL finding
  `R4-ccx-timeout` (the CCX subprocess now runs under a timeout and fails typed).
- **2026-10-01, pass 3 (WU5, DONE)**: `tests/test_laminate_bend_twist_3d.py` - hand-authored
  C3D20R layer-wise judge with per-layer orientation and a consistent traction. Controls pass,
  the load-discretisation artifact is pinned, and the verdict is 3D/shell = 1.0021 at matched
  refinement: **H6 refuted**. Two methodological traps recorded (§13.2, §13.4), both of which
  first produced a wrong answer - the same failure mode as the Ritz load-vector bug in §12.8.
  Committed as `c2528cf` after RDD review `review-73c6a2cfcf1a75ae` (4/4, approved, no
  correction).
- **2026-10-01, pass 4 (WU7, DONE)**: `tests/test_blade_twist_mechanism.py` - the blade twist is
  measured, its mechanism is identified, and the issue's causal attribution is excluded (§14).
- **2026-10-01, pass 5 (E1, DONE)**: the section stiffness against the BeamDyn anchor (§15). The
  anchor deck was in the tree all along; the modal route puts every section stiffness inside the
  reference scatter, which **retracts** the ~400x torsional claim of §14.3 (it was a measurement
  artifact). Three estimators were tried for a per-station number and each carries a bias, so the
  per-station comparison is left as an open methodological item that needs the section frame.

---

## 14. WU7 - the blade-level twist (DONE, with one open lead)

### 14.0 What cannot be reproduced here

The issue's blade numbers (+8.9 deg shell vs the BeamDyn anchor +0.98 deg, and the FSI de-loading)
come from a campaign whose artifacts are absent from this tree (§3.1), and the BeamDyn anchor is
not in this tree either. So WU7 does not re-run that comparison. It measures the blade's twist
mechanism on the IEA-15-240-RWT model that *does* exist here, at the load level the validation
suite already uses (root moment 90.4 MNm, the article's DLC 1.4 maximum).

### 14.1 Finding 1 - the blade has no laminate bend-twist coupling at all

All **696 sections** are balanced: `max |D16| + |D26| = 0.000e+00` against a `|D|max` scale of
1.0. The issue's causal story ("the same bend-twist is the mechanism behind the blade's elastic
over-twist") is therefore **impossible for this model**: there is no coupling to over-predict.
The coupon work (§12-§13) and the blade mechanism are *not* the same question.

### 14.2 Finding 2 - the twist is a load-path response

One root moment, three lines of action for the resultant, same model, same mesh. The twist is the
mean nodal rotation about the blade axis at the tip ring:

| line of action | tip flapwise [m] | tip twist [deg] |
| --- | ---: | ---: |
| spread over all surface nodes | 21.691 | **-0.570** |
| leading edge | 20.553 | **-46.920** |
| mid-chord | 20.567 | **+23.947** |

The deflection barely moves (21.7 -> 20.6 m, 5%) while the twist swings across **70 degrees**.
That is the right order of magnitude to explain a 9x disagreement between two models that are each
internally consistent about a *different* load line of action - a shell spreading the aerodynamic
load over the surface and a beam applying it along the aerodynamic centre are not modelling the
same load. It also means the blade twist gap is a **convention** question before it is a physics
question.

### 14.3 Finding 3 - the mesh has free edges at two mid-span stations (bounded effect)

**44 free edges: 8 at z ~ 11.94 m and 2 at z ~ 112.22 m** - mid-span, at the shear-web junctions,
where the webs are not connected across those stations. The node dedup step reports 0% reduction,
so it does not merge them. That is a real connectivity defect.

**Its effect on the section stiffness is bounded** (see §15): the modal route puts the shell's 1st
torsion at 4.000 Hz against the reference 4.290 Hz (GJ ratio 0.87, reference scatter 15.1%), so the
tear is a defect to fix but it is **not** the blade-twist explanation. Two repair attempts are still
recorded because they are informative: merging coincident nodes naively degenerates the thin
transition elements (the matrix becomes singular), and a careful merge that skips pairs sharing an
element finds only *one* pair, so the tear is not a simple duplicate-node problem.

**RETRACTED (2026-10-01, E1).** This section previously reported the shell's torsion as ~400x too
soft, from a pure tip torque giving -4989.59 deg/MNm against a closed-section estimate of order
10 deg/MNm. **That was an artifact of the measurement, not a property of the model**: the metric was
the *mean nodal rotation at the loaded tip ring*, where the local shell deformation dominates - the
same contamination class as §12.8 and §13.2. The modal cross-check refutes it and a per-station
rigid-body fit gives GJ ratios of order 0.8-1.5. The claim was removed from
`tests/test_blade_twist_mechanism.py`, whose docstring now carries the retraction.

### 14.4 A metric lesson

The first metric tried was the LSQ slope of the flapwise displacement against the chordwise
coordinate at the tip. With a ~21 m tip deflection that measures the *bending* rotation, not the
twist: it read **-18.3 deg** on a case whose actual twist is **-0.570 deg**. The twist about the
blade axis is the mean nodal rotation about z. Same class of trap as the coupon's load and the 3D
model's traction.

### 14.5 What WU7 settles for the issue

- The element is not the cause (WU4), coupon plate theory is not the cause (WU5), and the blade's
  laminate coupling does not exist (14.1). Three of the issue's implied causes are excluded.
- The remaining candidates are the **load line of action** each model assumes (14.2) and the
  **section's torsional path**, whose measured softness (14.3) is an open defect in the mesh.
- The +8.9 deg figure is plausible for *some* load line of action and implausible for others; the
  number alone, without the load's line of action, cannot settle the disagreement.

---

## 12. Second pass — measured evidence (2026-10-01)

### 12.1 The coupon reproduces the reported ratio against the reported formula

`L = 1.0 m`, `b = 0.1 m`, 8 plies of 1 mm (`[45,0,0,45]s`), CFRP 120/10/5, `P = 10 N`, mesh
2x12. `issue20` is the issue's own formula; `free-edge` is the same free-edge state with the two
unit factors restored (`theta(L) = (D^-1)[2,0] P L^2 / (4 b)`).

| layup | D16 | D26 | shell [deg] | issue20 | free-edge | shell/issue | shell/free-edge |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `[45,0,0,45]s` | 701.93 | 701.93 | -0.287813 | -0.05799 | -0.28994 | **4.963** | 0.9926 |
| `[45,-45,-45,45]s` | 221.66 | 221.66 | -0.080064 | -0.01672 | -0.08361 | 4.788 | 0.9576 |
| `[0,45,45,0]s` | 480.27 | 480.27 | -0.181582 | -0.03719 | -0.18593 | 4.883 | 0.9766 |
| `[45,45,45,45]s` | 1182.20 | 1182.20 | -1.379655 | -0.30774 | -1.53871 | 4.483 | 0.8966 |
| `[30,0,0,30]s` | 895.67 | 320.11 | -0.481091 | -0.11046 | -0.55231 | 4.355 | 0.8710 |
| `[60,0,0,60]s` | 320.11 | 895.67 | -0.112080 | -0.01693 | -0.08466 | 6.619 | 1.3239 |
| `[0,0,0,0]` | 0.00 | 0.00 | 0.000000 | 0.00000 | 0.00000 | -- | -- |

The ratio against the issue's formula is 4.4-6.6, i.e. the issue's 3.4-4.8 is reproduced. The
same sweep with T300/5208 plies gives the same picture (shell/free-edge 0.88-1.13).

### 12.2 Not the mesh

`[45,0,0,45]s`, `NX x NY` refined in both directions (`shell/free-edge`):

| mesh | shell [deg] | shell/issue | shell/free-edge |
| --- | ---: | ---: | ---: |
| 2x12 | -0.287813 | 4.963 | 0.9926 |
| 4x24 | -0.283103 | 4.882 | 0.9764 |
| 8x48 | -0.282779 | 4.876 | 0.9753 |
| 16x96 | -0.282051 | 4.864 | 0.9728 |

Width refinement alone is flat from `NY = 12` (0.9759 -> 0.9752 for `[45,0,0,45]s`; 1.0732 ->
1.0709 for `[60,0,0,60]s`). The residual is a converged property, not a discretisation artefact.

### 12.3 Not the metric

`[45,0,0,45]s`, three metrics on the same displacement field (LSQ slope of the free edge, exact
corner slope, mean nodal rotation about the length axis): spread 0.2% at 2x12 and **0.0%** at
8x48.

### 12.4 The coupling itself is the free-edge state (constant-moment probe)

A constant tip moment is the one load case whose exact plate solution is the free-edge state
(`M = (Mx, 0, 0)`, constant curvatures), so it isolates the constitutive chain from the BVP. With
`Mx = -10 N` per unit width (repo sign convention, locked by an assertion on `w_tip < 0`):

| quantity | measured | CLT `D^-1 (Mx,0,0)` | ratio |
| --- | ---: | ---: | ---: |
| interior twist rate `theta'` | +1.016533e-03 rad/m | `kappa_xy/2` = +1.012095e-03 | **1.0044** |
| interior twist rate | +1.016533e-03 rad/m | `kappa_xy` = +2.024191e-03 | 0.5022 |
| tip twist | +9.998324e-04 rad | `kappa_xy/2 * L` = +1.012095e-03 | 0.9879 |

The `kappa_xy/2` column is the one that matches, which is the unit factor the issue's formula
drops.

### 12.5 Independent judge 1 — converged Rayleigh-Ritz (same BVP, other discretisation)

`w(x,y) = sum_{p>=2,q} a_pq (x/L)^p P_q(2y/b)`, minimising the CLT bending energy; clamped edge
imposed exactly by `p >= 2`; load = uniform free-edge line load. Machinery control: with the
**constant moment** the Ritz tip twist is `1.0046x` the free-edge state (the exact interior
solution), and for a decoupled laminate its tip deflection (6.4995e-03 m) sits between the beam
value (6.4616e-03) and the free-edge value (6.5104e-03), as the b/h = 12.5 theory requires.

| layup | shell (8x48) [deg] | Ritz (P=10, Q=24) | shell/Ritz |
| --- | ---: | ---: | ---: |
| `[45,0,0,45]s` | -0.282779 | -0.291250 | 0.9709 |
| `[45,-45,-45,45]s` | -0.084758 | -0.088373 | 0.9591 |
| `[0,45,45,0]s` | -0.178899 | -0.185380 | 0.9650 |
| `[45,45,45,45]s` | -1.479731 | -1.488781 | 0.9939 |
| `[30,0,0,30]s` | -0.510887 | -0.513121 | 0.9956 |
| `[60,0,0,60]s` | -0.090694 | -0.103072 | 0.8799 |

### 12.6 Independent judge 2 — hand-authored CalculiX S8R composite deck

Written from the CalculiX manual (not through `write_ccx_mesh`, whose composite-section block a
previous audit found defective, so a deck built by it could not judge the element). Controls:
the decoupled `[0,0,0,0]` layup gives **exactly** 0 twist in CCX too (deck, section, load path and
parser validated), and the `+/-(45,0,0,45)s` pair gives equal and opposite twists (per-ply
orientations actually realised).

| layup | CCX (8x48) [deg] | element (16x96) | element/CCX |
| --- | ---: | ---: | ---: |
| `[45,0,0,45]s` | -0.281413 | -0.282051 | 1.0023 |
| `[60,0,0,60]s` | -0.086234 | -0.087915 | 1.0195 |
| `[45,-45,-45,45]s` | -0.086028 | -0.085524 | 0.9941 |

Matched-refinement sweep for the two layups that look furthest apart at coarse meshes:
`[60,0,0,60]s` element/CCX 0.8745 -> 0.9440 -> 0.9694 -> 0.9809, and `[45,-45,-45,45]s`
1.0659 -> 1.0202 -> 1.0114 -> 1.0059. Both converge to 1 from below/above as expected of two
different element families on different meshes.

### 12.7 Then why is the free-edge formula not exact either?

Compatibility ties the curvatures of a Kirchhoff plate: `kappa_xy,x = 2 kappa_x,y`. The free-edge
state has `kappa_x = kappa_x(x)` and `kappa_xy = kappa_xy(x)`, both proportional to `M(x)`, so
`kappa_xy,x = r_xy M'(x) != 0` while `kappa_x,y = 0`. The free-edge state is therefore admissible
**only for a constant moment**, and for a varying moment the true plate solution must redistribute
`kappa_x` across the width. That is exactly the 1-13% residual that `shell/free-edge` shows, and
it is why the reference has to be a converged computation rather than a one-line formula. Both
converged routes (12.5, 12.6) put the element within ~3% of the truth.

### 12.8 Honest notes and limits

- **One of my own tools was wrong first.** The first Rayleigh-Ritz implementation evaluated the
  Legendre basis at a single quadrature point instead of over the full rule, which produced a
  nonzero twist on a decoupled laminate (physically impossible). The CCX route exposed it; the
  load vector was fixed and the control now passes at 1.0046. This is the concrete reason the
  chain "two independent routes" was the right decision: a single reference, however converged,
  would have carried the bug into the verdict.
- **The blade-level claim is untouched.** The coupon verdict removes the D16/D26 coupling as the
  cause of a 4.7x inflation; it does not explain a 9x blade gap, which lives in the
  section/theory question (shear deformation, warping, 3D effects, or the beam anchor's own
  assumptions). WU5/WU7 stay open, and the FSI de-loading numbers in the issue inherit whatever
  the blade answer turns out to be.
- **`pyright` cannot be run from the shell here** (`reportMissingImports` for numpy/pytest/scipy);
  the same errors appear on pre-existing test files, so it is the environment, not this file.
- **The issue's own artifacts remain absent** from this tree and from its history (§3.1): the
  coupon, its reference and its judge now exist because they were authored here, not recovered.

---

## 13. WU5 — the 3D layer-wise judge (DONE): H6 refuted

`tests/test_laminate_bend_twist_3d.py` (4 tests). Hand-authored CalculiX deck, never
`write_ccx_mesh`: one **C3D20R** element per ply, each layer carrying its own `*ORIENTATION` ply
angle, clamped `x = 0` face, loaded free tip, and the **same metric** as the shell (LSQ slope of
`w` against `y` over the tip edge) read on the mid-plane nodes.

### 13.1 Controls

| control | result |
| --- | --- |
| decoupled `[0,0,0,0]` | `-0.000000` deg (exact zero) |
| `[45,0,0,45]s` vs `[-45,0,0,-45]s` | equal magnitude, opposite sign (this is what validates the per-layer orientation; an all-zero layup would pass under any convention) |

### 13.2 The trap that first looked like a theory verdict

A **uniform force per node over the brick face is not a uniform traction**: it over-weights the
perimeter nodes and injects a spurious through-thickness moment. The resulting gap is
thickness- and mesh-independent, which is exactly what a "CLT is wrong by 21%" conclusion looks
like:

| load discretisation | mesh 10x6 | mesh 16x10 |
| --- | ---: | ---: |
| uniform nodal force over the whole face | 0.7876 | 0.7781 |
| mid-plane line only | 1.0058 | 1.0025 |
| **consistent traction** (8-node serendipity face integrals) | **1.0050** | **1.0021** |

Ratios are `3D / shell`. `consistent_traction_loads` integrates the face shape functions against a
constant traction, which is the honest 3D counterpart of a shell's mid-surface line load. The
difference is pinned by `test_a_uniform_nodal_force_on_the_face_is_not_a_uniform_traction` so it
cannot come back unnoticed.

### 13.3 Convergence and verdict

In-plane refinement with the consistent load (`3D / shell`): 6x4 = 1.0145, 10x6 = 1.0050,
16x10 = **1.0021**; the 3D answer settles **toward** the shell as it refines.

**H6 refuted at coupon scale**: real 3D elasticity, with full transverse constants the CLT cannot
see, reproduces the composite MITC4 shell twist to **0.2%**. The plate theory is not the coupon's
problem, and neither is the element.

### 13.4 Second trap: orphan nodes in a structured C3D20R grid

Nodes at an odd index in *both* in-plane directions are the would-be centre of a 20-node brick's
face: no element owns them, CalculiX reports no result for them, and a naive sampler dies on the
missing key. Sampling element-corner stations only avoids it. Related: the repo helper
`tests/_ccx_io.py::parse_frd_disp` keeps only the **last** `-4 DISP` block, which is fine for small
models and silently incomplete for large ones; the 3D module carries its own `parse_all_disp`.

### 13.5 What this does to the blade claim

Neither the element (WU4) nor coupon-scale plate theory (WU5) can explain a 9x blade twist gap.
WU7 therefore has to look at the **blade-level model itself** - mesh, load path, boundary
conditions, section properties, or the beam anchor's own assumptions - not at the shell theory.
That is a sharper starting point than the one the issue provided.

---

## 15. E1 - section stiffness against the BeamDyn anchor (DONE, with an open methodological item)

**The anchor was in the tree all along.** `.sources/openfast/iea15mw/` holds the BeamDyn blade file
(26 stations x 6x6 stiffness + mass matrices), the BeamDyn primary deck, the ElastoDyn blade file
and the AeroDyn decks. §14.0 said the anchor was absent - **that was wrong** and is corrected here.

The 6x6 ordering was taken from `openfast_toolbox.converters.beam.K66toPropsDecoupled(convention='BeamDyn')`
rather than from memory: `K[0,0] K[1,1]` shear, `K[2,2]` axial, `K[3,3] K[4,4] K[3,4]` bending,
`K[5,5]` torsion. At the root: EA 4.605e10 N, GKt 8.749e10 N.m^2, EI 1.4963e11 / 1.4973e11 N.m^2.
Those two bendings cross-check against ElastoDyn's `FlpStff`/`EdgStff` (1.5253e11 / 1.5248e11) to
2%, and that cross-check is what fixes the reading: it is what rules out the competing ordering in
which the torsion would sit at index 3.

### 15.1 The modal route (the reliable measure)

The validation suite already compares the shell's parked modes with the published references:

| mode | shell | reference | gap | reference scatter |
| --- | ---: | ---: | ---: | ---: |
| 1st flapwise | 0.526 Hz | 0.570 Hz (article) | 7.6% | 5.8% |
| 1st edgewise | 0.702 Hz | 0.650 Hz (article) | 8.1% | 11.8% |
| 1st torsion | 4.000 Hz | 4.290 Hz (NuMAD) | 6.8% | 15.1% |

The torsional frequency gives `GJ_shell / GJ_ref ~ (4.000/4.290)^2 = 0.87`.

**Verdict**: every section stiffness is inside the reference scatter. Neither the mesh tear (§14.3)
nor the section stiffness explains a 9x blade twist gap. With §14.1 (no laminate coupling) and
§14.2 (the twist is load-path dominated), the **load's line of action remains the prime suspect** -
which is what the issue's checklist items A1-A3 have to supply.

### 15.2 Open methodological item: a per-station number needs the section frame

Three estimators were tried for a per-station stiffness and each carries a bias of its own:

| estimator | failure |
| --- | --- |
| mean nodal rotation at the loaded ring | local shell deformation dominates: under-reports GJ by ~400x (this is the retracted §14.3 claim) |
| mean ring displacement + quadratic fit | the mean carries the section's twist times the ring's asymmetry (a real airfoil is not symmetric), giving sign-flipping EI |
| 6-parameter rigid-body fit on *global* axes | the structural twist (15.6 deg at the root) and the prebend mean the global-z rotation is not the section's torsion; it read a negative GJ at mid-span |

The reliable route is therefore the modal one above. A precise per-station comparison needs each
station's own frame (rotate by the structural twist, and use the section's principal axes), which is
a well-defined next step rather than a mystery.

### 15.3 What E1 changes in the plan

- The `~400x too soft` claim is retracted in the module, in §14.3 and on issue #9.
- The mesh tear stays a real defect with a *bounded* effect: worth fixing, not the explanation.
- The next decisive item is the load's line of action (checklist A1), because it is the only
  candidate left standing for the blade twist gap.

---

## 16. Verified references (extracted from the PDFs, not from the issue text)

The maintainer asked for the anchors to be extracted from the papers rather than trusted as
quoted. What follows is my own extraction, with the discrepancies found.

| source | what I read | status |
| --- | --- | --- |
| **Zhou et al. 2025**, *Energy* **336:138488** | Table 4, mean tip deflections at rated: ALM-GEBT flap 14.10 m / edge -1.27 m / **torsion -3.77 deg**; Present (LL-FVW&GEBT) 13.86 / -1.22 / **-3.60**. Table 6: rigid 16.11 MW / 2.53 MN, flexible 14.76 / 2.20, reductions 8.38% / 13.04%. Sign: flap positive downstream, edge positive toward the trailing edge; and "the aerodynamic moments at the blade sections make the airfoil sections **twist towards feather and reduce the angle of attack**" | **matches the issue's numbers exactly**; the issue's article number (136488) is wrong, it is **138488** |
| **Ma et al. 2025**, *Front. Energy Res.* **13:1571567** (DOI verified in the PDF) | torsion is **nose-down**; the yaw-sweep values are read off Figs. 15-17 - the paper has only Table 1 (parameters) and Table 2 (frequencies) | the issue's sweep numbers are **figure-read**, i.e. lower confidence than a table |
| **Gaertner et al. 2020**, NREL/TP-5000-75698 | blade mass **65,250 kg** (twice in the report); "worst-case out-of-plane tip deflection is **22.8 m**" | matches the issue |
| **Escalera Mendoza et al. 2023**, AIAA 2023-2093 | mass **68,077 kg**, "4.33% greater compared to blade model reported in Ref. [1]"; Table 3 modes **0.57 / 0.65 / 1.72 / 2.08 / 3.41 / 4.29** (1F, 1E, 2F, 2E, 3F, **1T**) | matches the issue and the repo's own test constants |
| **Bernardi et al.** | the published version is *Wind Energ. Sci.* **11:2345-2367 (2026)**, DOI 10.5194/wes-11-2345-2026; `wes-2025-120-*` are the discussion preprints | the issue cites the **preprint** DOI; the flap ~16 m figure is not yet verified from the published text |
| **Zhou, "stiffness Table 6"** | Zhou's Table 6 is **power/thrust**, not stiffness | the issue's table mapping is wrong |

### 16.1 The finding that reorders the problem: the anchors disagree in sign

Three numbers, two conventions, one physical question:

| source | tip torsion at rated | physical sense (per Zhou's own text) |
| --- | ---: | --- |
| Zhou 2025 (LL-FVW + GEBT) | **-3.60 deg** | toward feather / nose-down, reducing the angle of attack |
| Ma 2025 (LL-FVW + GEBT) | about **-3.9 deg** (figure-read) | nose-down |
| OpenFAST/BeamDyn anchor (per the issue) | **+0.98 deg** | sign not yet established |
| our shell, one-way S-8c (per the issue) | **+8.9 deg** | sign not yet established |

Zhou's own physics statement is unambiguous: at rated the aerodynamic pitching moment is nose-down,
so the *physically expected* sense is Zhou's negative one. That makes the **sign** - not the
magnitude - the first thing to settle, and it changes the size of the disagreement: if our +8.9 deg
is the same physical sense as Zhou's -3.60 deg, the gap is **2.5x**, not 9x; and the BeamDyn
anchor's +0.98 deg would be the outlier, in the opposite sense to both the literature and the
aerodynamics.

Corroboration by magnitude: our FSI de-loading is -25%/-35% and Zhou's is **-8.38%/-13.04%**
(verified from Table 6). The ratio of the excess (~2-3x) matches the ratio of the twist (~2.5x),
which is what a coherent single-cause story looks like - and it is a different story from "9x".

## 17. Plan after these findings

### 17.1 The load cases are not the same load case (from the maintainer's own A1)

The maintainer's A1 answer settles a structural point: the **one-way S-8c path** - the source of the
+8.9 deg - applies the AeroDyn normal and tangential distributed loads to every section node with
spanwise tributary weights and **does not apply the pitching moment** (`apply_ad_loads` loads only
the two translational DOFs per node). The **FSI path** does apply the aerodynamic moment about the
aerodynamic centre (`force_projection` adds `Mp[k]*dr`). And the OpenFAST/BeamDyn anchor runs the
official coupled deck, which applies the aerodynamic pitching moment too.

So the reported comparison is between a shell that omits the pitching moment and an anchor that
includes it. That is a concrete, checkable difference, and it is the cheapest thing to close.

### 17.2 Work units

**WU-A - settle the sign (no external input needed).** Map all four numbers into one physical
criterion: "twist that reduces the angle of attack (toward feather) is positive". For the shell,
rotate the section rotation into the *section's own frame* (using the station's chord direction from
the undeformed mesh) instead of the global z, and state the sense. For the BeamDyn channel `RDxr`,
read what the rotation is measured about and with which sign, and map it the same way. Deliverable: a
sign table with the four entries and a test that pins it.
*Acceptance*: the shell, the anchor and Zhou/Ma are expressed in one sense, and the sense is
justified physically (nose-down at rated), not by convention.

**WU-B - load-case parity.** Re-run the one-way case **with** the AeroDyn pitching moment about the
aerodynamic centre added, everything else identical, and report the twist before and after.
*Acceptance*: the twist with the pitching moment is compared against the FSI path's own number
(+8.8 deg per the issue) and against Zhou, in the settled sense.

**WU-C - anchor reconciliation.** With WU-A and WU-B done, compare shell vs BeamDyn on the *same*
load case, or state precisely which load-case difference remains. The anchor's own load line of
action is the aerodynamic centre (official deck), so this is where the two should line up.
*Acceptance*: either the two agree within the reference scatter, or the residual is attributed to a
named difference in load or boundary condition - not to "the element".

**WU-D - magnitude against the literature and the FSI claim.** Compare the settled twist against
Zhou (-3.60 deg, verified) and Ma (about -3.9 deg, figure-read) and the de-loading against
-8.38%/-13.04% (verified). *Acceptance*: a single table in one sense, one load case, with the
figure-read source marked as such.

**WU-E - cleanup, in this order.** (1) The mesh tear at z ~ 11.94 / 112.22 m: real, bounded, worth
fixing (the dedup does not merge the coincident nodes). (2) The per-station section stiffness needs
the section frame - three estimators failed for three different reasons (§15.2). (3) State the
`D16 = D26 = 0` fact prominently: the blade has no laminate bend-twist coupling, so the coupon's
mechanism is not the blade's, and issue #9's original causal story should be corrected there.

### 17.3 Standing rule for this investigation

Four times now the first measurement was an artifact of my own scaffolding: the Ritz load vector,
the 3D model's load discretisation, the blade twist metric, and the tip-torque stiffness metric.
Each was caught by a *second independent route*, never by re-checking the same path. So: before
comparing two models, fix (a) the metric's physical meaning and (b) the load's line of action; and
when a number looks dramatic (400x, 9x), look in the scaffolding first.

---

## 18. WU-B second pass - the beam reference with the full load set (OPEN finding)

### 18.1 What was built

A beam-level twist reference from the **anchor's own section data** and **our own loads**, so the
comparison is a self-consistency check between two of our models rather than a literature quote:

- `GJ(z)` and the **shear centre** `(xS, yS)` per station from the BeamDyn 6x6 via
  `openfast_toolbox.converters.beam.K66toPropsDecoupled(convention='BeamDyn')` (the library's own
  convention handling, not a hand-rolled reading of the matrix);
- the **pitch axis** per station from the ElastoDyn blade file (`BlastFract`, `PitchAxis`);
- the aerodynamic loads `Np`, `Tp`, `Mp` from our own BEM at the rated point;
- the effective distributed torque about the shear centre,
  `m_eff = Mp + (x_AC - xS) * Np - (y_AC - yS) * Tp`, with `x_AC = (0.25 - pitch_axis) * chord`;
- `theta(z) = int_0^z T_eff(s)/GJ(s) ds` with `T_eff(z) = int_z^L m_eff ds`.

### 18.2 The tip station is degenerate and must not be integrated through

The BeamDyn deck's last station sits at the very tip where the section is a point: `GKt` runs from
**8.749e10** at the root down to **5.9e4** at the tip, a factor of 1.5e6. Integrating `T/GJ` through
that last interval alone contributes about **-9.4 deg** and swamps the result (the naive
tip twist comes out at -10.5 deg). The comparison therefore has to be made **station by station,
stopping the integral at the station**, which is what the table below does.

### 18.3 The measured discrepancy (same loads, same section data)

| station | beam | shell | shell / beam |
| --- | ---: | ---: | ---: |
| 80.1% span | -6.28 deg | -1.25 deg | **0.20** |
| 89.8% span | -8.60 deg | -6.46 deg | **0.75** |
| 95.2% span | -9.86 deg | -7.23 deg | **0.73** |

### 18.4 The internal inconsistency that has to be resolved first

The same comparison with the **pitching moment alone** points the other way:

| load set | shell | beam | shell / beam |
| --- | ---: | ---: | ---: |
| pitching moment alone | -7.35 deg | -5.07 deg | **1.45** |
| full load set | -2.72 deg (tip) | -9.86 deg (95%) | **0.73** |

Adding the forces flips the shell from *softer* than the beam to *stiffer* than it. Since E1 (§15)
independently measures the shell's torsional stiffness as **0.87x** the reference (i.e. softer, so a
`shell/beam > 1` is expected), the flip is not a property of the structure - it points at **the
force application in the shell**, which is the part that differs between the two rows: the forces
are spread over the ring and then moved to the aerodynamic centre with a corrective couple.

So this is an **open finding**, and per `CONTRIBUTING.md` ("Test rules", rule 3) it is recorded as a
residual rather than asserted: the test asserts what is physical (the sign, and that the line of
action dominates) and reports the magnitude.

### 18.5 Localisation plan (in this order)

1. **Controlled load-path test.** Apply a *single* known eccentric force (a pure `Np` with a known
   lever arm about the section's own centroid) and compare the shell's twist against the analytic
   `theta = int T/GJ` with the *same* lever arm. That isolates the corrective-couple implementation
   from everything else.
2. **Frame reconciliation.** The anchor's `xS` is measured from the pitch axis; the mesh's `x_AC` is
   built from the chord extremes. Verify that both reduce to the same physical point at two or
   three stations (a geometry check, no solve).
3. **Only then** re-run the full comparison. If the flip survives, the residual is a genuine
   beam-vs-shell difference (warping restraint) and belongs in the validity envelope.

### 18.6 The application-defect chain, and the withdrawal of the magnitude

The controlled eccentric-force test (step 1 of §18.5) did its job: it did not close the gap, it
**localised the cause**, and it produced a chain of load-application defects in my own harness.
Each was measured:

| # | defect | evidence | fix |
| --- | --- | --- | --- |
| 1 | a section moment applied as a **force pair** at the LE/TE nodes | erratic twist ratios against the analytic `int T/GJ`: 1.94 / 2.04 / 2.59 / 3.39 / 2.30 | apply it as a **shear flow** `q = M/2A` (how a closed thin-walled tube carries torsion); the ratios become consistent: 1.79 / 1.82 / 1.87 / 1.95 / 2.02 |
| 2 | the analytic reference integrated the torque only up to the measurement station (`m(z-s)`) instead of out to the tip (`m(L-s)`) | it shifted every ratio by 1.1-1.5x | corrected; the shear-flow ratios then settle at 1.8-2.0 |
| 3 | a distributed force applied as an **equal force per node** | the node spacing around a real airfoil is not uniform, so that is not a uniform traction | apply it as a **consistent traction** (length-weighted per edge) |
| 4 | the **sign** of the corrective couple | moving a resultant in y from `x_c` to `x_ac` needs `(x_ac - x_c) * F`, not `(x_c - x_ac) * F` | corrected |
| 5 | the resultant position under a consistent traction | it sits at the **edge-length-weighted centroid**, not at the mean of the node positions | corrected |

**The magnitude is therefore withdrawn as a result.** For the *same* load set:

| application | tip twist |
| --- | ---: |
| pitching moment alone (shear flow) | -9.79 deg |
| Np+Tp equal force per node, no moment (the previous one-way convention) | -14.25 deg |
| Np+Tp equal force per node + moment | -4.35 deg |
| Np+Tp consistent traction at the aerodynamic centre + moment | -34.13 deg |
| **Zhou 2025 Table 4** | **-3.60 deg** |

A factor of **7.8** between the smallest and the largest, from the application alone. The
`0.755x Zhou` reported in the previous revision of the rated test and in the matrix was computed
with defects 1-5 and is **withdrawn**. **The 7.8x spread itself was also computed with a
load-magnitude defect**: the per-metre aero loads were applied once per mesh station with no
spanwise tributary weight on a 266-bucket, 117 m, tip-refined mesh, so the applied resultant was
1.72-2.19x the BEM's own integral. With the corrected loads the spread is **4.29x** by `theta_z`
and **5.46x** by the section rotation `omega` - see §20.

What survives is what does not depend on the application:

- the **sign**: nose-down at rated, the same physical sense as Zhou's -3.60 deg and Ma's about
  -3.9 deg, opposite to the +0.98 deg previously attributed to the BeamDyn anchor;
- the **rigid aero side**: our BEM 2.525 MN / 16.327 MW against Zhou's 2.53 MN / 16.11 MW;
- the **measured shear centre**: 0.477 of the chord (the twist under a uniform flapwise line
  load vanishes there), consistent with the deck's pitch axis;
- the **application rules**: a section moment is a shear flow, not a force pair.

**Next step (bounded, and the only defensible one).** Validate the application on a case whose
answer is exact *before* touching the blade again: a **rectangular closed thin-walled tube** of
known dimensions and laminate under a known torque, where `theta = T L / GJ` with
`GJ = 4 A^2 / oint(ds / (G t))` is exact. The application is validated when the shell reproduces
that inside the 5% rule; only then can the blade comparison be promoted from a reported residual
to an asserted row.

---

## 19. WU-B2 step 2 - the load application and the metric on an exact case

**What was built.** `tests/test_thin_walled_tube_torsion.py` (4 tests). A rectangular **closed**
thin-walled tube, mid-surface mesh, axis along z, loaded only at the rings. The reference is
hand-written Bredt-Batho and independent of every path the repository computes. All numbers were
measured on this tree with the module's own `-s` output; the module reports `4 passed`.

| item | value |
| --- | --- |
| mid-line section | `b = 1.0 m` (x) x `h = 0.6 m` (y), corners at ( +-0.5, +-0.3 ), centroid at the origin |
| mesh | `n_seg = 4` per side (16-node closed ring, corners shared), 0.2 m elements along z; `L = 6, 12, 24 m` |
| isotropic wall | MITC4 (code 4): `E = 2.1e11`, `nu = 0.3`, `t = 0.01 m` (and 0.5 mm for the thin-wall comparison) |
| laminate wall | MITC4Composite (code 44): `[45,-45]s` CFRP, `E1 = 120e9`, `E2 = 10e9`, `G12 = 5e9`, `nu12 = 0.3`, ply `0.125e-3 m`, `t_total = 0.5e-3 m` |
| load | `T = 1.0e4 N.m`; a ring shear flow `q = T/(2A)`, `A = b h = 0.6 m^2`; or a single-wall traction |
| reference | `theta' = T/GJ`, `GJ = 4 A^2 A66 / perimeter`, `perimeter = 2(b+h) = 3.2 m`; `A66 = G t` (isotropic) or `t_total * Qbar66(45 deg)` (laminate) |

**A correction to the reference as the plan wrote it.** §18.6 wrote the reference as
`GJ = 4 A^2 / (perimeter / t)`. That expression drops `G`: it has units of `m^4`, not `N.m^2`, and
would make every ratio meaningless. The dimensionally consistent form - the one the plan's own
opening line and its laminate row already had - is `GJ = 4 A^2 A66 / perimeter`, i.e.
`4 A^2 / oint(ds / (G t))`. It is the form used and asserted here; nothing was tuned around it.

The two twist metrics are computed and reported in every case: `theta_fit` (in-plane rigid-ring
fit) and `theta_z` (mean DOF-5 rotation about z). In the exact Saint-Venant solution the section
does not distort in plane, so both must equal the section rotation.

### 19.1 Self-equilibrated torsion - the boundary-layer-free primary validation

`+T` at the tip ring and `-T` at the root ring (both as `q = T/(2A)`) is fully self-equilibrated:
all six rigid-mode resultants vanish (measured `<= 3.7e-18` normalised) and the realised end torques
are `+-T` to `1e-9`. With no clamped ring there is no boundary layer and the solution is uniform
torsion. The six rigid modes are removed by the exact constraint `Phi^T u = 0` (see the penalty
finding below); the removed solution satisfies `|Phi^T u| / |u| < 1e-10`.

| material | `theta_fit` ratio | `theta_z` ratio | metric difference | profile linearity (0.1-0.5L vs 0.5-0.9L) |
| --- | ---: | ---: | ---: | ---: |
| isotropic `t = 10 mm` | 1.00309 | 1.00674 | 0.36% | 0.0000% |
| laminate `[45,-45]s` | 0.99860 | 1.00703 | 0.84% | 0.0000% |

Both metrics reproduce `T/GJ` inside 5%, they agree with **each other** inside 2%, and the profile
is exactly linear - the two half-windows have the same slope to the printed precision. **This is the
metric validation:** with no end boundary layer, `theta_fit` and `theta_z` are the same section
rotation, and the disagreement seen in the clamped case is a boundary-layer artefact, not a property
of the element.

**The rigid-mode penalty is not exactly solution-neutral (measured).** The plan prescribes a penalty
`k Phi Phi^T` with `k = 1e-6 * median(diag K)` and an invariance check at `1000k`. Measured: the
assembled `K` annihilates the analytic rigid modes only to round-off, and the twist rate drifts by
**1.93e-2** (isotropic, `k = 5.60e2`) and **1.06e-1** (laminate, `k = 1.78e0`) between `k` and
`1000k`. The physical spectrum floor of `K` was measured at `lambda_min ~ 5.06e3`, so the prescribed
`k` sits only ~9x below it and `1000k` above it. A penalty at `1e-14 * median(diag K)` reproduces the
exact constrained solution to `~1e-10`, but the exact saddle-point constraint `Phi^T u = 0` is
k-free and is what the tests use. No tolerance was widened; the prescribed penalty simply does not
meet the `1e-9` invariance it was supposed to.

### 19.2 Clamped root - which metric is the section rotation, and how long is the layer

Clamped root ring, tip shear flow, `theta_fit` and `theta_z` over `0.4L-0.9L` and over the fixed
window `2.4-5.4 m`:

| t | L | 0.4L-0.9L `theta_fit` | 0.4L-0.9L `theta_z` | 2.4-5.4 m `theta_fit` | 2.4-5.4 m `theta_z` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10 mm | 6 m | 0.9916 | 1.0063 | 0.9916 | 1.0063 |
| 10 mm | 12 m | 1.0315 | 1.0067 | 1.0162 | 1.0062 |
| 10 mm | 24 m | 1.0082 | 1.0068 | 1.0169 | 1.0062 |
| 0.5 mm | 6 m | 0.8825 | 1.0070 | 0.8825 | 1.0070 |
| 0.5 mm | 12 m | 0.8933 | 1.0070 | 0.8898 | 1.0070 |
| 0.5 mm | 24 m | 0.9750 | 1.0070 | 0.9207 | 1.0070 |

Local secant slope of `theta_fit` / `T/GJ` along z (station z in m):

- `t = 10 mm, L = 6 m`: `0.1:0.771 0.7:0.914 1.3:0.941 1.9:0.961 2.5:0.975 3.1:0.985 3.7:0.992`
  `4.3:0.996 4.9:0.998 5.5:0.999`
- `t = 0.5 mm, L = 6 m`: `0.1:0.765 0.7:0.879 1.3:0.882 1.9:0.882 2.5:0.882 3.1:0.882 3.7:0.883`
  `4.3:0.883 4.9:0.883 5.5:0.885`

**Verdict.** `theta_z` converges to `~1.007` of `T/GJ` for every L, both thicknesses and both
windows - it is the section rotation. `theta_fit` does **not** converge: it is polluted by the
clamped-root boundary layer (which raises it from `0.77` near the root and can overshoot) and, for
the thin wall, by a thickness-dependent section in-plane shear that keeps it below `1` even at
`L = 24 m`. The test asserts only the asymptotic metric (`theta_z` inside 5%) and reports the other.

### 19.3 The one-wall traction - a statically equivalent load that is not equivalent

At the tip, a uniform `+y` traction on the wall at `x = +b/2` has the same `T_eff = 1.0e4 N.m` and
zero net axial force as the shear flow. It is not equivalent:

| load | `theta_fit` ratio | `theta_z` ratio | `alpha = du_x/dy` | `beta = du_y/dx` | shear `(alpha+beta)/2` | departure from rigid rotation | energy | tip in-plane |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| shear flow (tip) | 0.9916 | 1.0063 | -9.31e-5 | +7.14e-5 | -1.09e-5 | 0.125 | 0.827 J | 7.65e-5 m |
| wall traction (`x=+b/2`) | 20.1512 | 0.9366 | +1.28e-3 | +1.45e-3 | +1.37e-3 | 0.858 | 57.22 J | 5.72e-3 m |

The distortion amplitudes are at `z = 3.0 m`; the one-wall shear is **125.7x** the shear flow's, its
elastic energy is **69.2x** and its tip in-plane displacement **74.8x**. The mirrored wall (same `+y`
traction at `x = -b/2`) is exact: `|slope_left| = |slope_right|` to the last digit, opposite sign.
The two metrics disagree by `2050%` here because `theta_fit` reads the section parallelogram as
rotation while `theta_z` does not.

**Verdict.** The closed shear flow is the application that reproduces the exact answer. A localized
end load on one wall does not: it excites a soft section-distortion mode (a parallelogram field)
that the exact Saint-Venant traction does not. Saint-Venant's principle is **not** a usable
justification for redistributing a torque over the ring of this flat-faceted tube.

**What the blade investigation takes from it.** A thin-walled closed section loaded at one end with a
clamped root can have its measured twist dominated by a section-distortion mode rather than the
section rotation. The blade's 7.8x application spread (§18.6) is a candidate for the same mechanism,
and the next work unit must **measure the blade's distortion** (the section parallelogram and the
departure from rigid rotation), not just its rotation. Promoting the blade row from a reported
residual to an asserted one remains the NEXT work unit; nothing here changes a blade number.


---

## 20. WU-B3 - the blade's section distortion and the rated twist re-derived

**What was built.** Two tests added to `tests/test_blade_rated_twist.py` (4 -> 6): the
section-distortion test and the load-resultant invariant. Both use the four application vectors
`_rated_load_cases` provides (`mp_only`, `uniform`, `uniform_plus_mp`, `at_ac`). The four original
tests were not touched.

**A load-magnitude defect that invalidated this section's first numbers (measured).**
`_rated_load_cases` multiplied each station's per-metre load by nothing: the blade mesh has 266
raw z-buckets over 117 m (tip-refined, non-uniform), so the applied resultant was inflated and
concentrated. The BEM's own per-blade integrals are `I_Np = 8.4748e+05 N`, `I_Tp = 1.0371e+05 N`,
`I_Mp = -2.7659e+05 N.m` (trapezoid over the 50 radial stations; `bem.thrust/3 = 8.4169e+05 N`
agrees with `I_Np`). Measured on the assembled vectors before the fix:

| vector | force / moment | ratio to the BEM integral |
| --- | ---: | ---: |
| `at_ac` | `F_y = +1.4896e+06 N` | **1.76x** `I_Np` |
| `uniform` | `F_y = +1.8536e+06 N` | **2.19x** `I_Np` |
| `mp_only` | `M_z = -4.7659e+05 N.m` | **1.72x** `I_Mp` |

The per-metre intensities were always correct (at the BEM stations `Np/(q c) = 1.12..1.48`, a sane
`Cn`, and `Mp/(q c^2) = -0.099..-0.148`, a sane `Cm`), so this was a summation defect. Two
corrections:

1. **Spanwise tributary weight.** Every per-station contribution is multiplied by its tributary
   length `dz[k]` (interior stations `(z[k+1] - z[k-1])/2`; the two end stations carry half a cell,
   so `sum(dz)` is exactly the 117 m blade length - the plan's literal end cells over-counted by
   0.47%).
2. **Physical stations.** The blade is **prebent**: a physical ring's nodes span ~1e-3 m in z, so
   the raw unique-z set split each ring into near-duplicate buckets and the 1e-6 in-plane edge test
   found no edges at **62 of the 266 buckets**. The edge-based paths (consistent traction, shear
   flow) therefore applied nothing at 14% of the normal force; the node-based `uniform` path never
   had the problem. Merging z values closer than `STATION_GAP_TOLERANCE = 0.02 m` recovers 186
   physical stations whose rings are closed contours.

`test_rated_aero_loads_reproduce_the_bem_resultants` asserts the corrected resultants against the
BEM's integrals to 0.5%. Measured:

| vector | quantity | applied | ratio | bound |
| --- | --- | ---: | ---: | ---: |
| `at_ac` | `sum(f_y)` | +8.481267e+05 N | **1.000767** | 0.5% |
| `at_ac` | `sum(f_x)` | +1.033349e+05 N | **0.996377** | 0.5% |
| `uniform` | `sum(f_y)` | +8.469990e+05 N | **0.999436** | 0.5% |
| `mp_only` | `sum(x F_y - y F_x)` | -2.765723e+05 N.m | **0.999924** | 0.5% |

The corrected `at_ac` maximum nodal displacement (one-way) is **17.156 m** against Zhou's verified
coupled tip flapwise deflection **13.86 m** (ratio 1.24): the same order, which is the independent
confirmation that the load magnitude is now right. That comparison is **reported, not asserted** -
a one-way load against a coupled aeroelastic solution has no defensible numeric bound.

**The first-pass §20.2-§20.4 numbers were computed with the defect and are replaced below.**

### 20.1 The measurement: affine and rigid ring kinematics

A ring is the set of nodes at one constant z (`shell["rings"][zz]`), so it is planar. A helper
`_ring_kinematics(coords, u, ring)` fits the ring's in-plane field `(u_x, u_y)` at `(x, y)` by least
squares (`np.linalg.lstsq`) to the six-parameter affine map

    u_x = a11 x + a12 y + tx
    u_y = a21 x + a22 y + ty

using `(x, y)` as they are (the fitted translations absorb the origin), and reports:

| key | definition | meaning |
| --- | --- | --- |
| `omega` | `0.5 (a21 - a12)` | the **section rotation** [rad], the antisymmetric part, invariant to the origin |
| `shear` | `0.5 (a12 + a21)` | the **parallelogram** distortion (`exy`), the tube's measure |
| `dilatation` | `a11 + a22` | the relative area change (breathing) |
| `distortion` | `sqrt(a11^2 + a22^2 + 2 shear^2)` | the symmetric (strain) part magnitude |
| `residual` | `||u_ip - (A x + t)|| / ||u_ip||` | the non-affine part (0 if exactly affine) |
| `rigid_rotation` | fitted rotation of the **strain-free** model (2 translations + 1 rotation about the ring node mean) | to compare against the mean DOF-5 `theta_z` |
| `rigid_residual` | the same non-affine departure for that rigid-only model | how much of the field is not even a rigid motion |

The ring node ordering is **not used** for anything (the blade rings do not follow the contour
order), so there is no shoelace area and no ordering dependence.

### 20.2 Measured numbers (tip ring)

IEA-15-240-RWT at rated (V = 10.59 m/s, 7.56 rpm, pitch 0); dimensionless except the angles, the
energy and the displacement. Reproduced with
`python -m pytest tests/test_blade_rated_twist.py -v -s`.

| application | `theta_z` [deg] | `omega` [deg] | `rigid_rotation` [deg] | `shear` | `dilatation` | `distortion` | `residual` | energy [J] | max in-plane [m] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mp_only` (pitching moment, shear flow) | -4.8337 | -4.7851 | -4.7468 | +7.850e-04 | -4.713e-03 | 4.988e-03 | 0.0013 | 3.9006e+03 | 8.561e-02 |
| `uniform` (Np+Tp equal force per node) | -8.4522 | -10.4193 | -11.6487 | -2.100e-02 | -3.893e-02 | 5.068e-02 | 0.0003 | 2.8267e+06 | 1.699e+01 |
| `uniform_plus_mp` (per node + Mp) | -5.7845 | -11.2589 | -14.7356 | -6.120e-02 | -2.824e-02 | 9.202e-02 | 0.0007 | 2.8346e+06 | 1.701e+01 |
| `at_ac` (consistent traction at the AC + Mp) | -20.7348 | -26.1215 | -29.5274 | -5.975e-02 | -3.734e-02 | 9.375e-02 | 0.0007 | 2.9392e+06 | 1.708e+01 |

`omega` along the span [deg] at 50 / 80 / 100 % of the 117 m span (merged physical rings):

| application | 50 % | 80 % | 100 % |
| --- | ---: | ---: | ---: |
| `mp_only` | -0.9223 | -3.0287 | -4.7851 |
| `uniform` | -1.8881 | -7.0825 | -10.4193 |
| `uniform_plus_mp` | -2.1965 | -7.9352 | -11.2589 |
| `at_ac` | -5.0706 | -17.6041 | -26.1215 |

(The test prints exactly these three probes, so this table is reproducible from its `-s` output.)

The affine residual is small everywhere (0.0003-0.0013), so the ring field **is** well described by
the affine map: the split into rotation and strain is meaningful, not a fit failure. The merged
physical rings (12-26 nodes each) are what the fits run on.

### 20.3 Spreads and the verdict on the hypothesis

| quantity | value |
| --- | ---: |
| `spread_omega = max\|omega\| / min\|omega\|` | **5.459** |
| `spread_theta_z = max\|theta_z\| / min\|theta_z\|` | **4.290** |
| `ratio_omega = omega_at_ac / Zhou (-3.60 deg)` | **7.2560** |
| `ratio_theta_z = theta_z_at_ac / Zhou (-3.60 deg)` | **5.7597** |
| `max(distortion[uniform], distortion[uniform_plus_mp]) / distortion[at_ac]` | **0.982** |

**Status of the four assertions, on the corrected loads.**

| assertion | measured | status |
| --- | --- | --- |
| every application gives `omega < 0` (nose-down) | -4.79 / -10.42 / -11.26 / -26.12 deg | **green** |
| `spread_omega > 2.0` (the spread survives the distortion-free metric) | 5.459 (`spread_theta_z` 4.290) | **green** |
| `distortion[at_ac] > 5.0 * distortion[mp_only]` (the tube's contrast, stated on the correct pair) | 9.375e-02 vs 5 x 4.988e-03 = 2.494e-02, i.e. **18.8x** | **green** |
| promotion guard `abs(1 - ratio_omega) > 0.05` | 7.256 | **green** |

**Two predictions were tested and falsified, and each was replaced by the true restrictive claim.**
A known-false assertion must not stay in the suite as a false claim, so the two refuted predictions
were rewritten - not softened - into the claims the measurement supports; the refutations themselves
are recorded here and in the test docstring, with their raw numbers.

1. **"The application spread lives in the section distortion, not in the rotation" - refuted.**
   The distortion-free section-rotation spread (**5.459**) is *larger* than the mean-`theta_z` spread
   (**4.290**), so removing the symmetric part does not shrink the spread: the application moves the
   **section rotation itself**. The section-rotation metric also gives a *worse* residual against Zhou
   (`ratio_omega = 7.256` against `ratio_theta_z = 5.760`). The asserted claim is therefore
   `spread_omega > 2.0` - the spread survives the distortion-free metric, so the rated twist magnitude
   is **not** promotable - with the promotion guard against Zhou beside it.
2. **The tube's "off-path excites more distortion" contrast - mis-specified, then refuted.** Measured
   **0.982x**, not `> 1.5x`. The premise was wrong: the tube's contrast was a **localized one-wall
   traction against a closed shear flow**, while on the blade both compared cases are
   ring-distributed, and the `at_ac` path (a distributed consistent traction plus the Mp shear flow)
   excites as much section distortion as the off-path ones. The correct analogue of the tube's
   contrast is the **validated shear-flow moment application against the force-loaded cases**: an
   in-plane traction around the ring loads the section's in-plane flexibility, a pure torsion shear
   flow does not. That holds strongly - `distortion[at_ac] / distortion[mp_only] = 18.8x` - and is
   what the test now asserts `> 5x`.

The off-path numbers are still printed (the record keeps them); they are no longer asserted.

### 20.4 Residual against Zhou, and the blade row's status

The sign is green (every application gives `omega < 0`, nose-down at rated, the sense of Zhou's
-3.60 deg). The magnitude residual is **not** inside 5 %: `ratio_omega = 7.2560`, `ratio_theta_z =
5.7597`. The distortion-free metric makes the `at_ac` residual *worse* against Zhou, not better.

**The load magnitudes are now correct and invariant-checked.** The four resultants reproduce the
BEM's own integrals inside 0.5 % (`at_ac` `F_y` 1.000767, `at_ac` `F_x` 0.996377, `uniform` `F_y`
0.999436, `mp_only` `M_z` 0.999924), and the corrected `at_ac` tip displacement (17.156 m) is the
same order as Zhou's coupled 13.86 m. The application is validated against the aero side; the
remaining 5.76x residual is a structure/reference difference, not a load-magnitude artefact.

**The blade row stays a reported residual.** It does not become promotable: the promotion guard
`abs(1 - ratio_omega) > 0.05` holds (7.2560), the section-rotation spread (5.459) survives and is
now *larger* than the `theta_z` spread, and the sign is the only settled claim. Promotion still
requires removing the application spread itself - a load case whose line of action is settled
against the reference - which is a separate work unit.

### 20.5 What the corrected loads reveal: the twist is largely a sectional deformation

With the loads now invariant-checked, the physical application `at_ac` gives a section rotation
`omega = -26.12 deg` and a section strain `distortion = 9.375e-2`, i.e. `distortion / |omega| =
0.359`: the section strains are the **same order as the rotation**. The measured "twist" of this
shell blade at rated is therefore substantially a **sectional deformation**, not a rigid section
rotation - which is exactly the mechanism the tube test isolated (a thin-walled closed section
whose in-plane flexibility the load path excites), now quantified on the blade itself.

The dominant contribution is the **load's line of action**, not the pitching moment:
`mp_only` (the distributed pitching moment as the validated shear flow) gives `omega = -4.79 deg`,
while `at_ac` (the same aero side with the normal force at the 25 %-chord aerodynamic centre) gives
`-26.12 deg`. The eccentricity torque of the normal force about the section's own shear centre is
about `Np * (x_ac - xS)` with the mesh's measured shear centre at **0.477 of the chord** (§15.2) and
the aerodynamic centre at 0.25 c, so the eccentricity term is several times the pitching moment
itself. The `mp_only` result is also the one that independently agrees with the beam built from the
anchor's own 6x6 sections (§18.4: shell/beam = 1.45 with `Mp` alone, consistent with the measured
GJ ratio 0.87), which is what makes the eccentricity the prime suspect rather than the structure.

**Next unit (WU-C, now with a sharp entry point).** Compare the mesh's per-station shear centre
against the BeamDyn anchor's `xS` from its 6x6 sections (via
`openfast_toolbox.converters.beam.K66toPropsDecoupled`, the same route §14.2 used for GJ). If the
anchor places the shear centre near the pitching axis and the mesh places it at 0.477 c, the
eccentricity - and therefore the twist - is a **structural/geometry difference to settle at the
anchor**, not a load-application defect; if the two agree, the residual moves back onto the load
case and the aero model. That comparison is a geometry/section-data check, needs no large solve, and
is the cheapest decisive step left.

### 20.6 Third defect, found in review: the ring kinematics must use the merged physical rings

The RDD review of this unit raised a CRITICAL on the new kinematics: the fit ran on
`shell["rings"][zz]`, the fixture's **raw 1e-6 z buckets**, while the load path had already moved to
the merged physical stations. On a **prebent** blade a raw bucket is only a *slice* of a physical
ring, so the affine split could silently run on a partial section. Measured at the three probes:

| probe | nearest raw 1e-6 bucket | merged physical ring |
| --- | ---: | ---: |
| 50 % | 14 nodes | **16 nodes** |
| 80 % | 14 nodes | 14 nodes |
| tip | 12 nodes | 12 nodes |

At the 50 % ring the raw bucket was two nodes short, and the fitted `omega` moved from `-0.898` to
`-0.9223` deg when the complete ring was used (a 2.7 % change; the tip and 80 % rings were already
complete, which is why the tip table and the spreads are unchanged). The fixture now exposes
`phys_stations` / `phys_rings` / `phys_tip`, the new test's kinematics and tip mean use them, and a
guard asserts that each merged ring is a **superset** of the nearest raw bucket (the exact defect
class). The four original tests keep the raw sets they were written against; moving them to the
merged rings is a separate question for whoever next touches their assertions.

### 20.7 Review state of this unit

The RDD review could not be closed by **transport**, not by findings. The candidate was frozen twice
(lineages `review-4e1db2062df50a6e` then `review-95f93ea2e95257a0`; tier medium, one consolidated
`review-reliability` lens) and the lens relay returned an **empty reviewer output** on four attempts
(`stopReason: length` twice, `stop` twice, ~150 s and ~2 s respectively); the first attempt's payload
was also cut mid-JSON at byte 283 with "3 objects opened, 1 closed". No slot was ever consumed and no
acknowledgement was burned in either lineage; both remain open and unconsumed.

That first, partial payload did surface the **CRITICAL recorded in §20.6** (the ring kinematics fitted
on the raw 1e-6 z buckets), which is fixed and now guarded. Local evidence at this commit:
`python -m pytest tests/test_blade_rated_twist.py -v` -> **6 passed**; `ruff check` clean; exactly the
three files changed. The maintainer authorized committing this unit **with the partial review on the
record**, which is what this section documents.

---

## 21. Production defect P1: the ForceProjector places the aerodynamic centre at 0.75 c

This section records a **production** defect in the load path of issue #9, found while re-deriving the
blade twist, and its fix. It is the same family as the harness datum defect #6 (a datum chosen from
the wrong end of the chord) but it lives in `src/`, not in a test, so it affects the numbers of the
whole campaign: both `BEMStandaloneSolver` (`standalone.py`) and the BEM FSI participant
(`fsi_participant.py`) construct this projector on the real blade mesh.

### 21.1 The code site

`src/aeroelast/solvers/bem/force_projection.py`, `ForceProjector.__init__`. The leading edge of every
chordwise strip was inferred as the node with the **minimum** chordwise projection:

```python
chord_proj = strip_pts @ chord_dir
le_proj = float(np.min(chord_proj))
ac_proj = le_proj + station.airfoil.aerodynamic_center * station.chord
```

`chord_dir` is the PCA chord axis of the strip, oriented by an **unrelated** reference axis
(`normal_direction`, default `[1, 0, 0]`). The strip's `+x` therefore had nothing to do with the
leading edge; `min` merely selected whichever end that reference axis pointed away from. The AC is
then placed 0.25 c from that end, and the offset `(centroid_proj - ac_proj) * chord_dir` is fed to
`_distribute()` as the moment arm. Because the arm term usually dominates `M_AC`, the strip moment
`M_strip = M_ac + cross(self._strip_ac_offsets[k], F_strip)` carried a ~0.5 c force-couple error.

### 21.2 Proof A - the pitch-axis mapping

On the repo's own IEA-15-240-RWT blade mesh (generated by `Blade(str(tests/IEA-15-240-RWT.yaml),
element_size=1.0)` + `generate_mesh()`), the ring's chordwise extent is
`[-(1 - pitch_axis)·chord, +pitch_axis·chord]` exactly at `z = 68.05 m` and at the tip, with
`pitch_axis` read from the ElastoDyn blade deck. The cell spanning the negative end is the
`1 - pitch_axis` side, and the **mesh `+x` points toward the leading edge**.

### 21.3 Proof B - the blunt/sharp kink angle

At `z = 68.05 m` the interior angle at the `min(x)` node is **9.35 deg** (a sharp corner) and at the
`max(x)` node it is **97.3 deg** (a blunt corner). The aerodynamic leading edge of an airfoil is the
blunt one, so `min(x)` is the **trailing** edge. Independent check with the geometry-only rule now in
the code (in-plane node spread in the outer quarter of the chord): on the pre-fix grid the blunt end
is the `+chord_dir` end at every non-empty strip.

### 21.4 The measurement

On a 266-strip real-class measurement (task brief) the AC landed at **0.75-0.94 c from the true
leading edge in 241 of 266 strips**. The independent 50-station AeroDyn set measured here (`_openfast_bem.
build_blade_aero_from_aerodyn`, all airfoils at `aerodynamic_center = 0.25`) gives a uniform
**0.750 c from the true LE** on all 48 non-empty strips (the other two strips were empty, see P4), i.e.
a datum error of `0.50 c` per strip.

### 21.5 The moment arithmetic

For a strip of chord `c` and force `F`, displacing the AC by `~0.5 c` along the chord changes the
applied moment by `~0.5 c × F`. At mid-span (`c ≈ 4.7 m`, `Np·dr ≈ 1400 N`) the tangential component
gives `~48 kN·m` per strip against the physical pitching moment `|Mp|·dr ≈ 10 kN·m`. The datum term,
not `M_ac`, dominated the strip balance. `ForceProjector.verify()` is blind to this because it
conserves **force only** (`force_error = ||F_mesh - F_bem||`); it never inspects `M_strip`, and no
existing test measured a moment or an AC position (the 10 tests in `test_force_projection.py` use a
synthetic flat-plate mesh and a synthetic `BladeAero`).

### 21.6 The test added

`tests/test_force_projection_ac_datum.py` (3 tests) uses the **real** blade mesh and the **real**
AeroDyn `BladeAero` and pins three independent properties:

1. every mesh node belongs to a strip (no hub-offset datum gap, no empty strips);
2. the applied AC sits at `airfoil.aerodynamic_center` of the chord from the **geometrically**
   identified LE (the blunt end), within `0.02 c`;
3. the projected nodal forces reproduce the analytic moment about the origin
   `sum_k [r_ac_k x F_k + Mp_k·dr_k·span_dir]` with the **3-D** geometric AC point, within `1 %`.

Assertion 3 is exactly the check `verify()` lacks.

### 21.7 The fix

In `force_projection.py`, the `min`-end inference was replaced by the geometry-only blunt/sharp rule:
the LE is the end whose in-plane thickness (the spread perpendicular to `chord_dir`, in the plane
normal to the span) is larger in the outer 25 % of the chord; the TE is the sharp end. The AC is
interpolated `ac_proj = le_proj + ac_frac * (te_proj - le_proj)`, which is sign-free, and ties
(`<= 10 %` thickness difference, e.g. the circular root) fall back to the `+chord_dir` end. Two
further corrections were needed to meet assertion 3's 1 % bound and are recorded here:

* the moment arm is now the **full 3-D** `centroid - r_ac` instead of its chordwise projection. The
  out-of-chord (span/thickness) component is not negligible: the chordwise-only arm left a `1.14 %`
  moment residual (mostly a span lever arm of `~0.3-0.8 m`) on the real blade;
* the arm's **sign** was inverted. `_distribute()` solves `sum_j d_j x f_j = M_strip` about the
  centroid, so the correct transfer is `M_strip = M_ac + (r_ac - r_centroid) x F = M_ac -
  cross(ac_offset, F)` with `ac_offset = centroid - r_ac`; the shipped code added `+cross(ac_offset,
  F)`, i.e. the opposite couple. The comment's formula `M_centroid = M_AC + r_{AC->centroid} x F` was
  itself wrong.

### 21.8 P4 - the hub offset in the strip assignment (resolved here)

Assertion 1 failed independently: the strip boundaries were built from `blade_aero.r` (hub-referenced,
`3.97..120.97 m`) while the single-blade mesh's span coordinate is blade-root-referenced
(`0..117 m`). Measured before the fix: **88 unassigned root nodes** (`z in [0.000, 2.388] m`) and
**2 empty strips** (`[48, 49]`), so the tip loads were dropped. The fix anchors the station grid on
the mesh's own span origin (`r_stations += span_coords.min() - r_stations[0]`): a shift by the hub
radius for a blade-root mesh, a no-op for a hub-referenced rotor mesh (and for the synthetic tests'
hub-referenced mesh). The unused `hub_radius` constructor argument now names that datum when a caller
supplies it. After the fix: 0 unassigned nodes, 0 empty strips.

The FSI participant's deformed geometry remains out of scope and is **not** touched: `_ref_r` is
hub-referenced while `r_def` (mean mesh span) is blade-root-referenced. That datum question is a
follow-up for whoever next touches `fsi_participant.py`.

### 21.9 Before / after

Datum, measured on the 50-station AeroDyn set (fraction of the chord from the true LE toward the TE):

| station | chord [m] | thickness lo / hi before | AC before | thickness lo / hi after | AC after |
| ---: | ---: | --- | ---: | --- | ---: |
| 0 | 5.29 / 5.22 | 3.513 / 4.303 | **0.750** | 4.032 / 3.806 (tie -> fallback) | **0.250** |
| 10 | 5.74 / 5.80 | 0.700 / 1.989 | **0.750** | 0.910 / 2.249 | **0.250** |
| 20 | 4.50 / 4.69 | 0.301 / 1.302 | **0.750** | 0.304 / 1.404 | **0.250** |
| 30 | 3.50 / 3.66 | 0.247 / 0.825 | **0.750** | 0.256 / 0.892 | **0.250** |
| 40 | 2.57 / 2.74 | 0.227 / 0.584 | **0.750** | 0.251 / 0.629 | **0.250** |
| 48 | (empty) / 1.88 | - | - | 0.183 / 0.506 | **0.250** |
| 49 | (empty) / 1.37 | - | - | 0.067 / 0.288 | **0.250** |

Worst datum error: `0.504 c` before (`0.0000 c` after).

Moment about the origin at the repo's rated point (`V = 10.59 m/s`, `7.56 rpm`, BEM built by the
repo's own engine over the same AeroDyn polars):

| quantity | before | after |
| --- | ---: | ---: |
| `abs(M_applied)` [N·m] | 6.473352e7 | 6.422547e7 |
| `abs(M_expected)` [N·m] | 6.400697e7 | 6.422547e7 |
| `abs(M_applied - M_expected)` [N·m] | **7.269812e5** | **2.987159e-8** |
| relative to `abs(M_expected)` | **1.1358e-2** | **4.6511e-16** |
| relative to `sum_k abs(Mp_k·dr_k)` | **2.671** | **~1.1e-13** |

The after-fix `|M_applied| = |M_expected| = 6.422547e7 N·m` is the geometric analytic moment to machine
precision: the force now acts at the 3-D AC point and the strip balance reproduces it exactly.

### 21.10 Guard rails

The 10 pre-existing tests in `test_force_projection.py` are unchanged and stay green (`python -m
pytest tests/test_force_projection.py -v` -> **10 passed**; the task text said 11, the file collects
10). `ruff check` is clean on the modified source and the new test. The new file is the regression
guard for P1 and P4.
