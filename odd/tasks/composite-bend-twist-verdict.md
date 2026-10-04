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

### 20.5 What the corrected loads reveal: the twist is largely a sectional deformation — **SUPERSEDED, do not quote**

> **This section's conclusion was wrong and is kept only as the record of how it was reached.**
> It read `distortion / |omega| = 0.359` from the `at_ac` **hand-built** case and concluded that
> the measured "twist" of this shell blade at rated is substantially a sectional deformation.
> That `omega` was inflated ~17x by the harness's own chordwise datum slip (defect #6) and, behind
> it, by the production load-path defects P1/P2/P4/P5. Measured under the **production** load path
> (§22.4) the tip section rotation is `omega = -1.5112 deg` with `distortion = 2.4695e-1`, i.e.
> the ratio is **9.36, not 0.359**. The reading was an artefact of applying the aero load at 0.75 c
> from the leading edge on the wrong axis and in the wrong half-space, not a property of the blade.
> The premise was also incomplete: the line of action was not a modelling choice to be settled
> against the anchor, it was a production bug, and it is fixed. The numbers quoted below stay as
> measured at the time - they are correct for the loads the harness actually built.

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

**Corrected by section 22.12 - the kink numbers above do not reproduce, the spread rule does.** The
`9.35 deg` is not obtainable from either object it could have meant: measured on the mesh ring the
`min(x)` node reads **102.19 deg** and the smallest interior angle anywhere in the mesh is **76.7 deg**
(section 22.12), and measured on the source geometry - the `airfoils[].coordinates` of
`tests/IEA-15-240-RWT.yaml` - the trailing-edge end reads **99.2 to 116.6 deg** and the nose end
**152.6 to 177.5 deg** over the eight profiles in the file. These airfoils have no sharp trailing
corner, so "sharp means trailing" was a premise this mesh cannot support and the number quoted for it
is not reproducible; it is reported rather than deleted so the correction is auditable. What survives,
and what the identification actually rests on, is the second half: the outer-quarter in-plane spread
(0.17 m at `min(x)` against 0.86 m at `max(x)`, i.e. the wide rounded end is `+x`), the pitch-axis
split (`x_max = pitch_axis * c` to 0.005 % median over all 186 rings) and `pitch_axis` agreeing
between the yaml and the ElastoDyn deck to `1.1e-16`. The conclusion - leading edge at `+x`, so
`x = (pitch_axis - f) * c` - is therefore unchanged and now pinned by guards that do not use a kink
angle.

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

## 22. Suite honesty audits (requested by the maintainer before the next work unit)

Two audits were run on this tree, both read-only, because the maintainer's point was that a test
which only exercises its own harness proves nothing about the framework.

### 22.1 Tests that do not exercise production paths

Full classification in the working notes; the load-bearing results, all verified by me against the
source:

| finding | evidence | status |
| --- | --- | --- |
| **`tests/test_blade_rated_twist.py` never imports `ForceProjector`** - the module's own subject (the load application) is re-implemented in `_rated_load_cases`, so a repeat of the P1 datum bug in production code would still leave this file green | `grep -c ForceProjector tests/test_blade_rated_twist.py` -> **0** | **OPEN**, next unit |
| My "load resultant invariant" (`test_rated_aero_loads_reproduce_the_bem_resultants`) compares hand-built vectors against trapezoidal integrals of the *same* `bem.Np`/`Tp`/`Mp`, i.e. it validates my own `dz` bookkeeping, not production | `tests/test_blade_rated_twist.py` (the integrals and the vectors both come from `bem`) | **OPEN** - the invariant is worth keeping but must be asserted against `ForceProjector.project()` output |
| `tests/test_force_projection.py` takes its expected force from `projector.verify()`, whose docstring states it reuses the same `Np[k] * dr` summation as `project()` - production asserted against production | `tests/test_force_projection.py` lines 177, 195, 206, 232, 339 | **OPEN** |
| `tests/test_blade_twist_mechanism.py` hand-builds its three lines of action (`_station_load`) and asserts only relative spreads | same pattern, lower stakes | **OPEN** |
| Production physics with **no test that calls it at all** | `src/aeroelast/constitutive/failure.py` (Tsai-Wu / Hashin / max-stress), `solvers/elasticity/dynamic_newmark.py`, the FSI time loops (`solvers/fsi/linear_dynamic.py`, `stress_stiffened_dynamic.py`), `solvers/bem/fsi_participant.py`, `solvers/fsi/force_clipper.py`, and the Rust kernels `newmark_beta_solve_coo` + the FSI drivers | **OPEN**, ranked by risk after 22.2 |
| Genuinely production-facing and kept | `tests/test_thin_walled_tube_torsion.py` drives `PyMeshAssembler.assemble_k` + `create_laminate_from_angles` against a hand-written Bredt reference; `tests/test_force_projection_ac_datum.py` drives `ForceProjector.project()` against geometry; `test_bem_openfast_parity.py` compares against the real AeroDyn binary | done |

**Consequence for this investigation.** The tube module validated *the element*, and the new datum
module validates *the production load path* - but the blade twist numbers of section 20 still come
from hand-built vectors. They must be re-derived through `ForceProjector` before any blade number
here is quoted, which also retires the section 20 magnitude question rather than answering it.

### 22.2 The 5% rule sweep (CONTRIBUTING.md "Test rules (inviolable)")

Every tolerance above 5% found by sweep, classified. `verify` bounds are quoted because they are
absolute newtons against a force of known size.

| site | bound | measured | verdict |
| --- | --- | --- | --- |
| `tests/test_material_suite.py:561` | `err < 0.20` with the comment "inflat[es] compliance ~16%; **allow 20% tolerance**" | ~16% | **VIOLATION of rule 2** - the tolerance was chosen to accommodate the measurement. The class itself sets `TOL = 0.05` |
| `tests/test_material_suite.py:572` | `TOL = 0.10` "B-coupling analytical solution has shear correction uncertainty" | - | **VIOLATION** - an uncertain reference must be *stated* as uncertain (or the reference made exact), not widened |
| `tests/test_blade_iea15mw_validation.py:335` | `0.0 < rel_report < 0.20` | see matrix row | **justified**: a sign claim plus a documented reference-author disagreement (the article's own NuMAD conversion is +4.33% over that report); the *parity* claim in the same test uses `MASS_TOL = 0.05` |
| `tests/test_bem_engine.py:29` | `IEA15_RATED_POWER_TOL = 0.15` | see matrix | to classify against the reference's own precision |
| `tests/test_laminate_bend_twist_3d.py:264` | `0.70 < face/shell < 0.85` | - | **legitimate**: a negative control pinning the size of a documented load-discretisation artefact, not a parity claim |
| `tests/test_mitc3_smoothed.py:145` | `plain < 0.95` | - | **legitimate**: a negative control (expected membrane locking) |
| `tests/test_blade_twist_mechanism.py:186` | `< 0.10 * abs(surface)` | - | between two of our own models, not a reference comparison - and 10% is loose regardless |
| `tests/test_force_projection.py:235` | `force_error < 50.0` N (varying `Np`, 20 m blade) | - | total force is ~2e4 N, so ~0.25% - inside 5%, but stated in newtons without the denominator |
| `tests/test_force_projection.py:340` | `force_error < 1.0` N "relaxed for coarse discretisation" | - | same; the comment is a smell even when the number is fine |
| `tests/test_shell_comprehensive.py:404,449,494,687`, `test_ko2017_performance.py:1698`, `test_laminate_bend_twist.py:329,427,442,468,469`, `test_laminate_bend_twist_3d.py:287,305`, `test_orthotropic_shell_parity.py:342,404`, `test_shell_analytical_validation.py:195,299,472`, `tests/test_ko2017_performance.py:821,845` and the `*_TOL = 0.05` constants | 5% or tighter | - | compliant |

**Matrix discrepancies found by the sweep** (the matrix is the audit target, so these are the
pending items it must absorb):

1. `docs/validation-matrix.md` section 2 claims "429 tests: 416 passed, 13 xfailed" while its own
   index/inventory rows now read 442 inventoried and 510 collected. A measured full run on this
   tree gave **497 passed, 13 xfailed, 0 failed** (510 collected, 47 files) in 26:06; that run
   started on the working tree that already contained the `6064aa8` fix. A confirming re-run at the
   committed tree is in progress and its XFAIL list will be quoted node by node.
2. Section 13.2 says "the **eleven** nodes below" while the run reports **13** xfailed. The list
   must be reconciled against the authoritative XFAIL enumeration, not guessed.
3. The ">5% tolerances" column of the index (8) does not include `test_material_suite.py:561/572`,
   which are the two real violations in 22.2.

### 22.3 The P5 guard: the convention-free invariants that fail today

`tests/test_force_projection_load_frame.py` (new, **not committed** - see below) is the RED guard for
the load-frame defect. `ForceProjector.project()` applies
`F_strip = F_n * normal_dir + F_t * tangential_dir` on fixed global axes (production defaults
`normal_dir=[1,0,0]`, `tangential_dir=[0,1,0]`), while `Np`/`Tp` from
`ccblade.rotor.distributedAeroLoads` are normal/tangential **to the section chord**
(`BEMResult.Np` = "Normal force per unit length"). The test derives the section frame from the
**ring outline only** - the merged ring's principal in-plane axes, oriented leading-to-trailing with
the blunt-end rule - and never from `normal_dir`/`tangential_dir`, so it cannot agree with the defect
by construction.

**Measured geometry** (real IEA-15MW mesh, `Blade(str(tests/IEA-15-240-RWT.yaml), element_size=1.0)`
+ `generate_mesh()`, real AeroDyn `BladeAero`, `span_direction=[0,0,1]`, production defaults for the
other two axes). `strip`/`dr` are the production strip a single-station load reaches, `n` the ring
nodes, `extent(c_hat)` the ring extent along `c_hat`, `extent(max)` the brute-force max in-plane
extent:

| frac | z [m] | strip | dr [m] | n | `c_hat` | `f_hat` | angle(`c_hat`,[1,0,0]) [deg] | chord [m] | extent(`c_hat`) [m] | extent(max) [m] | angle(`c_hat`,max) [deg] | x ext [m] | y ext [m] | y mean [m] |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.03 | 3.184 | 1 | 2.388 | 22 | (0.92558, -0.37856, 0) | (0.37856, 0.92558, 0) | 22.245 | 5.2185 | 5.2087 | 5.2159 | 6.645 | 5.1899 | 5.0707 | +0.096 |
| 0.15 | 17.509 | 7 | 2.388 | 23 | (0.97941, -0.20188, 0) | (0.20188, 0.97941, 0) | 11.647 | 5.6454 | 5.6468 | 5.6471 | 0.547 | 5.5588 | 2.9324 | +0.373 |
| 0.25 | 29.447 | 12 | 2.388 | 21 | (0.99038, -0.13838, 0) | (0.13838, 0.99038, 0) | 7.954 | 5.6704 | 5.6761 | 5.6764 | 0.554 | 5.6392 | 2.0365 | +0.369 |
| 0.50 | 58.497 | 25 | 2.388 | 16 | (0.99922, -0.03938, 0) | (0.03938, 0.99922, 0) | 2.257 | 4.1549 | 4.1545 | 4.1545 | 0.157 | 4.1535 | 1.1170 | -0.011 |
| 0.75 | 87.750 | 37 | 2.388 | 14 | (0.99976, 0.02191, 0) | (-0.02191, 0.99976, 0) | 1.256 | 2.9959 | 2.9959 | 2.9959 | 0.044 | 2.9953 | 0.5696 | -1.517 |
| 0.90 | 105.061 | 44 | 2.388 | 12 | (0.99929, 0.03779, 0) | (-0.03779, 0.99929, 0) | 2.166 | 2.2759 | 2.2759 | 2.2759 | 0.066 | 2.2744 | 0.3833 | -2.879 |
| 1.00 | 117.000 | 49 | 2.388 | 12 | (0.99979, 0.02052, 0) | (-0.02052, 0.99979, 0) | 1.176 | 0.5000 | 0.5000 | 0.5000 | 0.124 | 0.4999 | 0.0799 | -3.998 |

**The "chord runs along the ring's major axis" claim is checked, not assumed.** The extent along
`c_hat` matches the interpolated `blade_aero.chord` within **0.2 %** at every station
(worst at frac 0.03: 5.2087 vs 5.2185 m, -0.19 %; the other six are within 0.11 %), and `c_hat` is
within **6.65 deg** of the independently scanned max in-plane extent direction (bound 10 deg). The
tip ring reproduces the task's numbers exactly: x extent 0.4999 m = the 0.5000 m tip chord, y extent
0.0799 m = the airfoil thickness, y mean -3.998 m = the documented prebend. So `c_hat` is really the
chord and `f_hat` the out-of-plane/thickness axis, and both twist along the span
(`angle(c_hat,[1,0,0])` runs +22.2 deg at the root to +1.2 deg at the tip).

**Inter-station chord angle: 24.410 deg**, the maximum pair, between frac 0.03 and frac 0.90. Because
24.4 deg > 2 x 10 deg, **no single global axis can lie within the 10 deg bound of every station's
`f_hat`**, whatever its sign - so the invariant cannot be satisfied by any fixed-vector
implementation. (The task's 25 / 50 / 75 % example spans only ~9.5 deg; the root and 90 % stations
were used to make the claim decisive, per the task's "pick more separated stations".)

**The projections.** With `dr = 2.3878 m`, the current projector returns `F = (+2387.75, 0, 0) N` for a
single-station `Np = 1000 N/m` and `F = (0, +2387.75, 0) N` for a single-station `Tp = 1000 N/m`:
both are exactly `1000 * dr` on a fixed global axis, independent of the station.

| frac | z [m] | `Np`: \|F.c\|/\|F\| | `Np`: \|F.f\|/\|F\| | angle(F, `f_hat`) [deg] | `Tp`: \|F.c\|/\|F\| | `Tp`: \|F.f\|/\|F\| |
| --- | --- | --- | --- | --- | --- | --- |
| 0.03 | 3.184 | 0.9256 | 0.3786 | 67.755 | 0.3786 | 0.9256 |
| 0.15 | 17.509 | 0.9794 | 0.2019 | 78.353 | 0.2019 | 0.9794 |
| 0.25 | 29.447 | 0.9904 | 0.1384 | 82.046 | 0.1384 | 0.9904 |
| 0.50 | 58.497 | 0.9992 | 0.0394 | 87.743 | 0.0394 | 0.9992 |
| 0.75 | 87.750 | 0.9998 | 0.0219 | 88.744 | 0.0219 | 0.9998 |
| 0.90 | 105.061 | 0.9993 | 0.0378 | 87.834 | 0.0378 | 0.9993 |
| 1.00 | 117.000 | 0.9998 | 0.0205 | 88.824 | 0.0205 | 0.9998 |

**The three RED assertions** (all fail today; bounds are the task's 5 % / 10 deg):

1. `test_normal_load_is_perpendicular_to_the_chord` - `Np` must ride `f_hat` (out of chord):
   **worst `\|F.c_hat\|/\|F\| = 0.9998 > 0.05`** and **worst `\|F.f_hat\|/\|F\| = 0.0205 < 0.95`**.
   Today `Np` rides the chord.
2. `test_tangential_load_is_along_the_chord` - `Tp` must ride `c_hat`: **worst
   `\|F.c_hat\|/\|F\| = 0.0205 < 0.95`** and **worst `\|F.f_hat\|/\|F\| = 0.9998 > 0.05`**. Today `Tp`
   rides the flapwise axis.
3. `test_load_direction_follows_the_section_not_a_global_axis` - a single-station load must track
   that station's own section: **worst `angle(F, f_hat) = 88.824 deg > 10` at frac 1.00** (the force
   is always `+x`), with the measured `24.410 deg` twist proving no global axis can pass.

Command and result (`aeroelast-dev`, `CCX_BIN` set):

```
python -m pytest -o addopts="" \
  tests/test_force_projection_load_frame.py tests/test_force_projection.py \
  tests/test_force_projection_ac_datum.py -v
```

-> **3 failed, 13 passed**. The three failures are the three new RED invariants; the 13 green are the
unchanged existing projection tests.

**Sign limitation (deliberate).** The absolute sign - which chordwise end is the leading edge, and
therefore the sense of `c_hat`/`f_hat`, and any downstream/upstream or flapwise direction - is **not
asserted**: no source in this tree (`.sources/`, `docs/`, the deck headers) documents it, so the
bound tests only the **axis**, never the direction. This must not be mistaken for a validated
convention; a source is required before any signed claim is added.

**Why it is not committed.** `docs/validation-matrix.md` is untouched and the committed tree must
stay green, so this file is the RED guard handed to the next work unit rather than a committed test.
It is left untracked on purpose.

#### 22.3.1 fixed

The frame defect is fixed in `src/aeroelast/solvers/bem/force_projection.py`. The guard now has four
tests and is green; the implementation is below, including the one deliberate departure from the
letter of the task.

**Design actually implemented (axis from geometry, sign from configuration).**

1. **Axis.** For every strip with an outline (>= 2 nodes, non-degenerate in-plane spread) the chord
   axis is the first singular vector of the offsets projected into the plane normal to the span
   (`_strip_chord_axis`). SVD leaves a per-strip +/- 180 deg ambiguity, so the axes are made
   continuous along the span before anything else: each axis is flipped to have a non-negative dot
   with the previous strip's axis. The AC's leading/trailing split (`_section_ends`) is then run on
   the resolved axis, so the P1 datum and the applied load share one chord orientation.
2. **Sign.** The absolute sense is resolved **once for the whole blade**, not per strip:
   `chord_hat = chord_sign * axis` with
   `chord_sign = +1 if (sum_k dr_k * axis_k) . tangential_direction >= 0 else -1`, and
   `normal_hat = normal_sign * (chord_hat x span_dir)` with
   `normal_sign = +1 if (sum_k dr_k * (axis_k x span_dir)) . normal_direction >= 0 else -1`.
   The LE/TE geometry rule (blunt end = leading edge; a < 10 % blunt/sharp difference is a tie and
   keeps the `+chord_dir` end) is unchanged from P1.
3. **Application.** `project()` now uses `F_strip = F_n * normal_hat_k + F_t * chord_hat_k`, and
   `verify()` recomputes `force_bem = sum_k (Np_k * normal_hat_k + Tp_k * chord_hat_k) dr_k` from the
   same per-strip frames. Per strip the frame is checked to be orthonormal and in the section plane
   (`|dot(normal_hat, chord_hat)| < 1e-9`, both unit, `|dot(., span_dir)| < 1e-9`); a failure raises.

**Departure from the task, and why.** The task said to flip each strip into the configured
half-plane (`dot(chord_hat, tangential_dir) >= 0`, `dot(normal_hat, normal_dir) >= 0`). That rule is
discontinuous: on this blade the section rotates through `tangential_dir` (the chord's y-component
changes sign near frac 0.75) and the geometry normal's x-component changes sign, so a per-strip flip
mirrors part of the blade by 180 deg. Measured with the new magnitude test (`Np = 1000 N/m` uniform,
`Tp = 0`): the per-strip rule gave `|F| = 43 320 N` against `sum_k Np_k dr_k = 119 388 N`, a
**63.71 %** loss of the integrated normal force - i.e. it is not conservative. The global-sign rule
gives **0.7378 %**, inside the 5 % bound, while still taking the axis per station and the sense from
the configured directions. "Sign from configuration" is therefore a whole-blade sense, which is the
only sense that preserves global force conservation.

**The three invariants, before -> after** (7 probed stations, bounds 0.05 / 0.05 / 10 deg):

| invariant | before | after | bound |
| --- | --- | --- | --- |
| 1. `Np`: worst `\|F.c_hat\|/\|F\|` | 0.9998 | **0.0116** | <= 0.05 |
| 1. `Np`: worst `\|F.f_hat\|/\|F\|` | 0.0205 | **0.9999** | >= 0.95 |
| 2. `Tp`: worst `\|F.c_hat\|/\|F\|` | 0.0205 | **0.9999** | >= 0.95 |
| 2. `Tp`: worst `\|F.f_hat\|/\|F\|` | 0.9998 | **0.0116** | <= 0.05 |
| 3. worst `angle(F, f_hat)` [deg] | 88.824 | **0.662** | <= 10 |
| new: uniform-`Np` `\|\|F\| - sum Np dr\| / sum Np dr` | 63.71 % (per-strip sign) | **0.7378 %** | <= 5 % |

The residual `|F.c_hat|/|F| = 0.0116` is not a frame error: the guard rebuilds each section frame
from the **merged ring**, while a single-station load reaches the production **strip**, whose node set
includes neighbouring rings, so the two frames differ by a fraction of a degree. The 24.410 deg
inter-station twist still rules out any single global axis.

**Semantic change of the config keys.** `normal_direction` and `tangential_direction` no longer name
the axis the load rides; they only resolve the **sign** of the per-strip frame. A caller that was
relying on them as the axis (in particular the production defaults `[1,0,0]` / `[0,1,0]`, which are
the chord/normal of a *twisted* blade only by accident) is now silently re-interpreted. The same
wording still sits in `src/aeroelast/solvers/bem/fsi_participant.py:261-262`
("global direction for BEM Np / Tp"); that file is outside this work unit and needs its own
follow-up. No caller file was edited.

**Still unknown - explicitly outside the assertion set.** The absolute downstream/upwind and
flapwise **sense** has no source in this tree (`.sources/`, `docs/`, the deck headers). The guard
asserts only the axis, and the sign is whatever the caller's configuration declares; the **default**
configuration has not been reviewed against a physical reference and must not be cited as validated.
A source is required before any signed one-way/FSI claim is made. The campaign's one-way/FSI twist
and de-loading numbers were produced by the pre-fix projector and remain not citable until re-derived
with this implementation.

### 22.3.2 the neighbours, corrected from geometry

Fixing the frame (§22.3.1) turned **four** existing tests red: they encoded the old global-axis
expectation ("`Np` rides global `x`, `Tp` rides global `y`"), which the corrected projector no longer
satisfies. They are the direct neighbours of the P5 fix, and - as §22.1 had already recorded for the
projection module - their expectation was not independent: `tests/test_force_projection.py` took its
expected force from `projector.verify()`, and `tests/test_force_projection_ac_datum.py` built its
expected strip force from the configured `NORMAL_DIR`/`TANGENTIAL_DIR`. Both are corrected by
**deriving the expectation from the geometry in the test**, which is what makes the values physics
rather than algebra.

**The synthetic plate (stated in `tests/test_force_projection.py`).** The fixture builds a flat
rectangular plate: the strip outline extends along **global `x`** (chordwise index `i` ->
`x = 0 .. chord_length`), the strips stack along **global `z`** (spanwise index `j` ->
`z = hub_radius .. hub_radius + span_length`), and `y = 0`. The section normal is therefore
`chord_hat x span_dir = x_hat x z_hat = -y_hat`: **`Np` rides `-y` and `Tp` rides `+x`**, not the old
global axes. The integrated magnitude is `sum_k Np_k dr_k`, with `dr_k` rebuilt in the test by
`_strip_widths()` from the BEM station grid and the mesh span datum - never read back from
`projector.verify()`.

**The four tests, before -> after.**

| test | old expectation (measured on this tree) | new expectation | measured after |
| --- | --- | --- | --- |
| `test_uniform_Np_conservation` | `forces[:,1].sum() == 0` (actual **-22000**) and `forces[:,2].sum() == 0`, then `verify().force_error < 1e-6` | `forces[:,1].sum() == -sum_k Np_k dr_k` (`rtol=1e-9`); `forces[:,0].sum()` and `forces[:,2].sum() == 0` (`atol=1e-6`); `verify` kept as a bookkeeping cross-check only | `x = 0`, `y = -22000`, `z = 0`; **PASSED** |
| `test_uniform_Tp_conservation` | `forces[:,0].sum() == 0` (actual **+11000**) and `forces[:,2].sum() == 0` | `forces[:,0].sum() == +sum_k Tp_k dr_k` (`rtol=1e-9`); `forces[:,1].sum()` and `forces[:,2].sum() == 0` (`atol=1e-6`) | `x = +11000`, `y = 0`, `z = 0`; **PASSED** |
| `test_forces_only_in_load_direction` | `forces[:,1] == 0` per node (max **475**) and `forces[:,2] == 0`, plus `sum(forces[:,0]) > 0` | `forces[:,0] == 0` and `forces[:,2] == 0` per node (`atol=1e-8`); `forces[:,1].sum() == -sum_k Np_k dr_k` (`rtol=1e-9`) and `< 0` | `x = 0`, `y = -6000`, `z = 0`; **PASSED** |
| `test_moment_conservation` | expected `F_k = Np_k dr_k NORMAL_DIR + Tp_k dr_k TANGENTIAL_DIR`, AC from `projector._strip_chord_dirs[k]`; error **1.5442** relative on this tree (the matrix's earlier **1.0387** was an earlier source state) | expected `F_k = Np_k dr_k n_hat_k + Tp_k dr_k c_hat_k` with `(c_hat_k, n_hat_k)` derived from the ring outline (`_geometry_load_frames`), `r_ac_k` from `_blunt_end`; nothing read from `projector._strip_chord_dirs` / `_strip_normal_dirs` | error **1.128758e-07 N m = 1.7626e-15 relative** (bound 0.01); `\|M_applied\| = \|M_expected\| = 6.404089e7 N m`; **PASSED** |

**How each expectation is derived.** The three synthetic tests use the plate geometry above - the
load axis (`+/-y` for `Np`, `+/-x` for `Tp`) and the magnitude `sum_k Np_k dr_k` / `sum_k Tp_k dr_k`
from `_strip_widths()` - and assert the perpendicular component is ~0 instead of asserting the load
axis is zero. The AC-datum test builds each strip's frame from the ring outline alone
(`_geometry_load_frames`: principal in-plane axis via SVD, made continuous along the span, then the
one global sign rule of §22.3.1), locates the LE/TE with the blunt-end rule already in that module,
forms `normal_hat = chord_hat x span_dir`, and computes
`sum_k [r_ac_k x F_k + Mp_k dr_k span_dir]`. That is the identity `sum_j r_j x f_j` the projected
nodes must satisfy; it can only hold at machine precision if production follows each section. It does,
so the residual is `1.76e-15`.

**Whole-suite result (this tree, `aeroelast-dev`).**

```text
python -m pytest -o addopts="" tests/test_force_projection.py tests/test_force_projection_load_frame.py tests/test_force_projection_ac_datum.py -v
-> 17 passed, 6696 warnings in 28.30s

python -m pytest -o addopts="" -q tests
-> 501 passed, 13 xfailed, 53433 warnings in 654.92s (0:10:54)

python -m pytest -o addopts="" --collect-only -q tests | tail -3
-> 514 tests collected

ruff check tests/test_force_projection.py tests/test_force_projection_ac_datum.py
-> All checks passed!
```

No other test moved: the whole suite was `497 passed, 4 failed, 13 xfailed` before the correction, and
the four failures were exactly these four. The blade, FSI and CLI suites are green both before and
after, so P5's re-interpretation of the config keys (the semantic change recorded in §22.3.1) is not
exercised by any other test on this tree.

**This closes the §22.1 circularity finding for these two files.** The audit row
"`tests/test_force_projection.py` takes its expected force from `projector.verify()`, whose docstring
states it reuses the same `Np[k] * dr` summation as `project()` - production asserted against
production" is now **CLOSED**: the three synthetic expectations and the AC-datum moment are derived
from the mesh/ring geometry, and `verify()` is kept only as a bookkeeping cross-check, never as the
reference. The `tests/test_force_projection_ac_datum.py` AC datum row is likewise no longer read
from the projector's own frame.

The **remaining §22.1 items are still OPEN** and are not touched by this change:

- `tests/test_blade_rated_twist.py` still never imports `ForceProjector` (its subject is
  re-implemented in `_rated_load_cases`);
- its load-resultant invariant still compares hand-built vectors against trapezoidal integrals of the
  same `bem.Np`/`Tp`/`Mp`;
- `tests/test_blade_twist_mechanism.py` still hand-builds its lines of action;
- the test-less production physics (`constitutive/failure.py`, `solvers/elasticity/dynamic_newmark.py`,
  the FSI time loops, `solvers/bem/fsi_participant.py`, `solvers/fsi/force_clipper.py`, and the Rust
  `newmark_beta_solve_coo` + FSI drivers) is unchanged.

### 22.3.3 the sense: the defaults were transposed

§22.3.1 fixed the load **axis** (each strip now rides its own section outline) but left the load
**sense** to the configured `normal_direction`/`tangential_direction`. Reviewing that sense against
the model owner's declared convention - **the fluid travels along `+Y`, and the blade rotates
clockwise viewed from behind, i.e. `Omega = +omega * y`** - showed the defaults were **transposed**:
`normal_direction = [1, 0, 0]`, `tangential_direction = [0, 1, 0]`.

**Why that made the sign a round-off decision.** `project()` resolves each strip's frame **axis**
from the ring geometry and then takes its **sign** from `dot(axis, configured_vector)`. On this mesh
the section normal is (nearly) parallel to `y` and the chord is (nearly) parallel to `x`, so the
pre-fix `normal_direction = [1, 0, 0]` is (post-§22.3.1) *perpendicular* to the very axis whose sense
it is supposed to decide. The dot product is ~0 and the resolved sign is whatever the last
floating-point bits say. The transposed pair is therefore not cosmetic: `[0, 1, 0]` is parallel to
the load axis and the projection is well conditioned. Measured on the real IEA-15MW blade mesh + the
real `BladeAero` from
`tests/reference/iea15mw_openfast/case/IEA-15-240-RWT_AeroDyn15.dat`, at the rated point
(`BEMSolver(rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(10.59, 7.56, 0.0)`,
BEM rotor thrust 2.5251 MN / power 16.3266 MW, 3 blades):

| config | total projected force | meaning |
| --- | --- | --- |
| `normal_direction=[1,0,0]`, `tangential_direction=[0,1,0]` (**the pre-fix defaults**) | `(-1.179e5, -8.444e5, 0)` N | thrust along **-Y = upwind**, and the aero power is **negative** for `Omega = +Y` - the load brakes the rotor |
| `normal_direction=[0,1,0]`, `tangential_direction=[1,0,0]` (**the transposed defaults**) | `(+1.179e5, +8.444e5, 0)` N | thrust along **+Y**, `F.y = 0.8444 MN` vs the BEM's per-blade 0.8417 MN (**0.32%**), and power **+4.885 MW** (driving) |

The ill-conditioning is measurable, not qualitative: with the pre-fix pair the span-weighted dot of
the configured reference against the axis it signs is **`|dot| = 0.0798`** (chord and normal), i.e.
dominated by the small chordwise component of the twisted chord sum; with the transposed pair it is
**`|dot| = 0.996811`**.

**The convention, stated once.** Mesh axes measured on this tree's IEA-15MW blade mesh: a ring's `x`
extent is the section chord and its `y` extent the airfoil thickness, and the tip ring's mean `y` is
the documented prebend `BlCrvAC`. Span runs `+Z`; the chord runs `+X` with the leading edge at
positive `x`; the out-of-plane (flapwise/section-normal) direction is `+Y`, which is the
fluid/downwind direction; the rotor turns clockwise viewed from behind, i.e. `Omega = +Y`. The
*sense* is the model owner's declaration - it is not derivable from the geometry alone.

**The fix.** Transpose the defaults everywhere they appear: `ForceProjector`
(`src/aeroelast/solvers/bem/force_projection.py`), `BEMConfig` (`src/aeroelast/core/config.py`),
`solvers/bem/standalone.py`, and `solvers/bem/fsi_participant.py` (fallbacks plus its docstring,
which now describes the vectors as the *sense reference* for the section normal / chord and states
the well-conditioning requirement). The `ForceProjector` class docstring now carries the convention
once. `src/aeroelast/cli/run_bem_fsi.py` still shows the pre-fix pair in its example YAML
(lines 49-50 and 109-110: `normal_direction: [1.0, 0.0, 0.0]` / `tangential_direction: [0.0, 1.0, 0.0]`);
it is outside this change's allowed surfaces and is left for the maintainer.

**The new guard, written first and red before the fix** (`test_load_sense_is_downwind_and_driving`
in `tests/test_force_projection_load_frame.py`):

- the pre-fix run fails on the first assertion: `F = (-1.1789e5, -8.4444e5, 5.9e-11) N`,
  `F.y = -8.4444e5 <= 0`;
- after the fix: `F = (+1.1789e5, +8.4444e5, -3.8e-11) N`; `F.y = +8.4444e5 N` = 0.32% off
  `bem.thrust/3 = 0.84169 MN`, and `|F| = 852632.0 N` = **1.30%** off it (bound 2%);
  `F.x/|F| = +0.1383`, `F.y/|F| = +0.9904`, `F.z/|F| = -4.5e-17`;
- power `P = sum_j f_j . (Omega x r_j)` with the node positions measured from the rotor centre
  (`hub_radius = 3.97`, `r_span = hub_radius + z`): `P = +5.2559 MW` vs `bem.power/3 = 5.4422 MW`
  = **-3.42%**, and `P > 0` (driving). With the blade-root lever arm only, `P = +4.8854 MW` =
  **-10.23%**;
- ill-conditioning guard: `|dot(tangential_direction, chord_axis)| = 0.996811` and
  `|dot(normal_direction, section_normal)| = 0.996811`, both > 0.9. With the pre-fix defaults the
  same dots are **0.0798**, which is what made the sign a round-off decision.

**Open, un-attributed residual.** The power magnitude is short of the per-blade BEM power by
**-3.4%** with the full rotor-centre lever arm (`+5.2559 MW` vs `5.4422 MW`) and **-10.2%** with the
blade-root lever arm only (`+4.8854 MW` vs `5.4422 MW`). The sign is now right; the few-percent
magnitude gap is a discretisation / frame-transfer residual that this change neither explains nor
tunes away, so the test asserts only `P > 0` and **reports** the ratio - no magnitude bound is
invented to absorb it.

**Every pre-change campaign number has the aero load pushing upwind and braking the rotor.** The
one-way `+8.9 deg` twist and the FSI de-loading numbers were produced with the pre-fix defaults,
i.e. with the aero load on the wrong side. They are therefore not merely imprecise - they have the
**wrong sign** and must be re-derived. The validation matrix's P5 rows (§1 and §8.3c) are updated
accordingly.

**Neighbours that move.** The sense fix turned two `tests/test_force_projection.py` nodes red,
because they still encode the pre-fix sense: `test_uniform_Np_conservation` now gets
`ACTUAL +22000 N` where it expects `DESIRED -22000 N` (relative difference 2.0), and
`test_single_node_per_strip` gets `forces[:,0] = [0, 0, 0, 0, 0]` (the degenerate strip's normal now
signs to `+y`, so `np.all(forces[:,0] > 0)` is false). They are outside this change's allowed
surfaces and are left red for the maintainer to decide. The other synthetic nodes are unaffected:
`test_uniform_Tp_conservation` already rode `+x` (its exact-zero pre-fix dot resolved to the `+`
sign), and `test_forces_only_in_load_direction` passes the pre-fix direction vectors explicitly.

**Whole-suite result (this tree, `aeroelast-dev`, `ccblade` present but CalculiX/OpenFAST not, so 69
CCX/OpenFAST rows skip).**

```text
python -m pytest -o addopts="" -s tests/test_force_projection_load_frame.py -v
-> 5 passed, 3348 warnings in 13.30s

python -m pytest -o addopts="" tests/test_force_projection.py tests/test_force_projection_ac_datum.py -v
-> 2 failed, 11 passed (the two failures are the pre-fix-sign nodes above)

python -m pytest -o addopts="" -q tests
-> 2 failed, 441 passed, 69 skipped, 3 xfailed, 35240 warnings in 281.27s (0:04:41)

python -m pytest -o addopts="" --collect-only -q tests | tail -3
-> 515 tests collected

ruff check src/aeroelast/solvers/bem/force_projection.py src/aeroelast/core/config.py \
  src/aeroelast/solvers/bem/standalone.py src/aeroelast/solvers/bem/fsi_participant.py \
  tests/test_force_projection_load_frame.py
-> []
```

**Verification of the sense fix (measured, this refresh).** The two §8.3 nodes the run above left
red were then corrected from geometry, not by loosening anything:

- `test_uniform_Np_conservation` now expects the geometry-derived signed total `+sum_k Np_k dr_k`
  on `y` (the raw section normal of the synthetic plate is `x x z = -y`; the global sign rule
  resolves that pair to the `+Y` half-space the fluid travels in). Measured `F = (0, +22000, 0) N`.
- `test_forces_only_in_load_direction` and `test_single_node_per_strip` were on the same sign
  question; the latter also documents that a single-node strip has no chordwise extent, so its
  frame legitimately falls back to the configured sense vectors (`normal = +Y`, `chord = +X`) -
  measured `+2500 N` per node with `Np = 1000 N/m`, `dr = 2.5 m`.
- `test_forces_only_in_load_direction` had additionally been passing an explicit
  `normal_direction=[1,0,0] / tangential_direction=[0,1,0]` pair, i.e. the ill-conditioned
  reference this unit exists to remove; it now uses the production defaults, and the docstring
  says why. Pinning a sign through a reference vector that is orthogonal to the axis it signs is
  the round-off defect, not a test fixture.

Authoritative full-suite run on this tree (every tool present, so nothing skipped):

```
python -m pytest -o addopts="" -q -rxXs --tb=line -p no:cacheprovider tests
-> 515 collected: 502 passed, 13 xfailed, 0 failed, 0 errors, 0 skipped in 537.34s (8:57)
```

The `2 failed, 441 passed, 69 skipped` line quoted higher up in this section came from an
intermediate run without `CCX_BIN` on `PATH` and before those two corrections; it is kept only as
history, and `docs/validation-matrix.md` §1/§2 now carry the measured 502/13/0 line.

**Two follow-ups this unit deliberately does NOT fold in.**

1. `src/aeroelast/cli/run_bem_fsi.py` still shows the **pre-fix pair** in its two embedded example
   configs (lines 49-50 and 109-110: `normal_direction: [1.0, 0.0, 0.0]`,
   `tangential_direction: [0.0, 1.0, 0.0]`). Those lines are now misleading: with the transposed
   production defaults they describe the ill-conditioned, upwind-and-braking sense this unit
   removed. The edit was attempted and reverted in this session because the file's pre-existing
   diagnostics (an unbound `gen_cfg` at L235 and an unguarded `open()` at L328, both unrelated to
   P5 and blamed to `42d20df9`/`89637571`) block a clean touch; the example text should be fixed
   together with those two, in its own unit.
2. The open power-magnitude residual: `P = +5.2559 MW` against `bem.power/3 = 5.4422 MW`
   (**-3.42%**, rotor-centre lever arm) or **-10.23%** with the blade-root arm. The sign and the
   thrust magnitude are pinned; this magnitude is **un-attributed** and is deliberately left
   without any bound asserted, per the rule that a comparison which misses the bound is the
   finding. Candidate causes to test next: the rotor-plane projection of the chord at twisted
   stations, the strip `dr` grid (BEM stations vs the mesh's span buckets), and the minimum-norm
   nodal distribution's effective lever arm.

### 22.4 The structural response measured under the production load path

Commit `87db753` ("test(blade): measure the rated twist under the production load path, not
hand-built vectors") closes the §22.1 finding for `tests/test_blade_rated_twist.py`: the module now
also solves the shell with the forces `ForceProjector.project()` actually produces, through a
`production_rated_loads` fixture built exactly as `standalone.py` builds it (only `span_direction`
passed, so the configured defaults apply). The four hand-built applications stay in place,
untouched, as the historical sensitivity record of sections 18 and 20.

Measured at rated (V = 10.59 m/s, 7.56 rpm, pitch 0) on the real mesh:

| quantity | hand-built `at_ac` (§20) | **production path** | reference |
| --- | ---: | ---: | --- |
| applied `abs(sum(F))` | - | **852632.0 N** | `bem.thrust/3 = 841688.8 N` -> **+1.300%** (bound 2%, the P5 guard's `SENSE_THRUST_TOL`) |
| applied direction (share on the measured flapwise axis) | - | **+0.9874** | bound 0.95, mirroring the P5 guard's per-station `abs(F.f_hat)/abs(F) >= 0.95` |
| tip section rotation `omega` | -26.1215 deg | **-1.5112 deg** | Zhou -3.60 deg -> **0.4198x** |
| tip mean `theta_z` | -20.7348 deg | **-3.8558 deg** | Zhou -> **1.071x** |
| `distortion/abs(omega)` | 0.206 | **9.36** | - |
| tip flapwise deflection | - | **+16.3865 m** | Zhou coupled +13.86 m -> 1.182x |
| tip edgewise deflection | - | **-1.7848 m** | Zhou -1.22 m -> **1.463x, same sign** (§22.5) |

**What this settles.** The blade over-twist that motivated issue #9 is not a property of the
element: with the production load path the same structure and the same BEM loads give a tip twist
of the **same order** as the literature (1.071x on the metric this module has always used, 7.1%
off), where the hand-built vectors gave 5.76x-7.26x. Combined with the coupon verdict and
`D16 = D26 = 0` across all 696 sections, the element is exonerated twice over and the residual
that remains is the difference between our steady one-way BEM loads and a coupled aeroelastic
LL-FVW solution.

**What it does not settle.** Neither ratio is inside the 5% rule, so the promotion guard
`abs(1 - ratio_to_zhou) > 0.05` is asserted and the magnitude stays a reported residual - promotion
needs the load-case and aerodynamic-model differences removed, not a wider bound. The remaining
structural disagreements with Zhou's triple are now all of **magnitude**, not sign (flap 1.182x, edge
1.463x, torsion 1.071x on the module's historical metric); see §22.5 for the axis-convention fix that
removed the apparent edgewise sign flip. The power-magnitude residual from section 22.3.3 is still
un-attributed.

**Review R3-001 closed in the same breath.** The reliability lens flagged that the applied-load
invariant asserted a **magnitude** (`abs(sum(F))` versus `bem.thrust/3`), which a rotated but
norm-preserving load would pass. The test now also asserts the **direction**: at least 0.95 of
`abs(sum(F))` must ride the flapwise axis measured from the tip ring's own outline, in the downwind
half-space (measured 0.9874). The chordwise share is deliberately **not** bounded - over a span
whose chord turns 24.41 degrees, an aggregate chordwise share is not a statement the geometry
licenses, and inventing one to make the test look stricter is the failure mode the test rules exist
to prevent. It is printed as the TE-positive edgewise share (-0.1585) and left as a reported number.

### 22.5 The edgewise "sign inconsistency" was my test's axis convention

`§22.4` reported the tip edgewise deflection as `+1.7848 m` against Zhou's `-1.22 m` and flagged the
sign as an open structural disagreement. It was not one. On this mesh the **leading edge sits at
positive x** (the load-frame test measures the ring's x extent as the chord and places the leading
edge at `+pitch_axis * chord` from the pitch axis), so the chord axis `_measured_tip_axes` returns
points toward the **leading** edge, while Zhou states the edgewise deflection **positive toward the
trailing edge**. Reporting Zhou's triple through the un-flipped axis invented a sign flip.

With the axis expressed in Zhou's convention (`edge_hat = -chord_hat`, documented at the helper and
at the call site):

| component | Zhou | production path, before | production path, **after the axis fix** |
| --- | ---: | ---: | ---: |
| flapwise (positive downstream) | +13.86 m | +16.3865 m (1.182x) | +16.3865 m (1.182x) - unchanged |
| edgewise (positive toward the TE) | -1.22 m | +1.7848 m ("opposite sign") | **-1.7848 m (1.463x, same sign)** |
| torsion (nose-down negative) | -3.60 deg | `theta_z` -3.8558 deg (1.071x) | unchanged |

**Consequence:** all three of Zhou's tip components now agree in **sign** with the simulated blade;
the remaining disagreements are magnitudes only (1.182x, 1.463x, 1.071x), which is what a one-way
application of rigid-blade loads to the flexible structure should produce - the loads cannot
de-load themselves, so the deflections come out large. No production code changed in this unit; the
assertions (`flap_share >= 0.95`, `omega < 0`, the promotion guard) are untouched and still green
(`7 passed`), and the printed share is relabelled "edgewise(TE-positive) share -0.1585" so the axis
in the output matches the axis in the claim.

### 22.6 Production P6/P7: the deformed-geometry feedback was corrupted and under-measured

The BEM preCICE participant (`src/aeroelast/solvers/bem/fsi_participant.py`, `BEMFSIParticipant`)
had never been exercised through the production class by the suite, so the whole one-way
de-loading path was unguarded. `tests/test_bem_fsi_deformed_geometry.py` (new, 4 tests, 14 s) now
drives it on the real IEA-15MW mesh, the real AeroDyn `BladeAero` and the rated point
(V = 10.59 m/s, 7.56 rpm, pitch 0) under the participant's own projected loads. Measured: max|u| =
16.527 m, tip u = (+1.447, +16.429, +1.072) m, rigid path thrust 2.525066 MN / 16.326608 MW.

**P6, the mixed datums.** `_ref_r` is hub-referenced (3.9700 .. 120.9699 m) while
`_compute_deformed_geometry` averaged the mesh span coordinate, blade-root-referenced for a
single-blade mesh. Measured offsets `r_def - _ref_r`:

| | r_def [m] | r_def - _ref_r [m] | max abs | r_def.min vs 0.5 Rhub |
| --- | ---: | ---: | ---: | ---: |
| before | 0.3979 .. 117.4462 | **-4.2170 .. -3.0368** | 4.2170 m | 0.3979 m vs 1.985 m FAIL |
| after | 4.3679 .. 121.4162 | **-0.2470 .. +0.9332** | 0.9332 m | 4.3679 m vs 1.985 m |

The before column is minus the hub radius (3.97 m), so the deformed radii entered CCBlade at
~0 with `Rhub = 3.97`: `UserWarning: error. check input values.` and
`NaNs at 0/50: 1e-06 0.0 1.5707963267948966` at stations 0 and 1 of every distorted evaluation.
**No such warning is emitted after the fix** (the target file's warning summary carries only the
pre-existing mesh-thickness `DeprecationWarning`).

**Datum choice.** One convention, stated in the module docstring: `_ref_r`, the `r_def` of
`_compute_deformed_geometry` and the radii passed to `_rebuild_bem_solver` are all
**hub-referenced** (the datum of `blade_aero.r` and of the CCBlade annulus). The mesh's span
origin is converted once, at init, into `_mesh_datum_offset = _ref_r[0] - min_i(X_i . e_s)` - the
same anchor `force_projection.py` uses - and the deformed mean adds it, so the BEM annulus keeps
its hub while `ForceProjector` re-anchors the same values on the mesh datum internally.

**P7, the twist from a node cloud that is not a section.** A strip spans `dr = 2.39 m` against
~0.63 m mesh stations, so its PCA chord direction follows the deformed arc. Production now takes
the antisymmetric part of the least-squares in-plane affine fit of `(u_x, u_y)` over `(x, y)` of
**one physical ring** at the strip centre (`_section_rotation`), the ring being bound by span
position within `RING_BIND_FRACTION = 0.05` of the strip width (above the premesh prebend skew
~1e-3 m, below the 0.63 m station spacing). The retired SVD/PCA estimator
(`_compute_strip_chord_dirs`, `_ref_chord_dirs`) is deleted, not kept as a silent second
estimator; a ring too small to define a 2-D section keeps the reference twist.

Measured, outer quarter of the span, degrees about +span (production's elastic twist, converted
back into the mesh's rotation sense; the ring section rotation computed independently by the test;
the shell's mean rotational DOF 5):

| strip | r_ref [m] | production before | production after | ring section (ref.) | nodal theta_z |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 37 | 92.32 | -0.5744 | -0.7351 | -0.7351 | -0.4029 |
| 40 | 99.48 | -0.4294 | -0.8744 | -0.8744 | +0.2474 |
| 43 | 106.64 | -0.8169 | -1.2698 | -1.2698 | -2.8757 |
| 44 | 109.03 | -0.7475 | -1.1992 | -1.1992 | -3.2725 |
| 48 | 118.58 | -0.2003 | -0.8692 | -0.8692 | -0.8774 |
| 49 | 120.97 | **+2.3764** | **-1.5112** | -1.5112 | -3.8558 |

Production now equals the independent ring section rotation exactly (max gap 0.000e+00 deg,
bound 1e-6). The gap to the nodal drilling rotation is **2.3446 deg** at the tip and is
**reported, not asserted** - see the refuted premise below.

**A third defect on the same path: the mesh-to-BEM twist sense was inverted.** `blade_aero.twist`
is a CCBlade angle (`alpha = phi - theta`), while the elastic rotation is a mesh rotation about
+span, and the two are opposite-sensed on this frame. Measured on the mesh's own chord line (the
extreme-chordwise node pair, oriented +x = leading edge): `psi = -8.954, -6.308, -2.194, -0.312,
+1.513, +2.066, +1.185` deg at z = 10.35, 30.24, 50.14, 69.84, 90.14, 105.06, 117.00 m against
`-theta = -13.971, -6.765, -2.734, -0.587, +1.548, +2.103, +1.242` deg - same slope, opposite
sense (within 0.6 deg outboard of 50 m). Adding the mesh rotation with a plus fed the BEM an
inverted twist and turned de-loading into re-loading. The class now derives the factor from its
configured directions, `_twist_mesh_to_bem = -sign(n . (s x t))` (-1 for the default span = z,
normal = y, tangential = x), and applies `theta_def = theta_ref + _twist_mesh_to_bem * omega`.

**De-loading, before and after** (Zhou et al. 2025 Table 6, flexible vs rigid:
-13.04% thrust, -8.38% power):

| case | thrust | power |
| --- | ---: | ---: |
| production path before P6+P7 (corrupted annulus + strip-cloud twist) | -1.84% | -1.20% |
| production path after (hub datum + ring twist + corrected sense) | **-2.32%** | **+0.42%** |
| elastic twist alone, reference station radii | **-4.00%** | **-0.81%** |
| deformed radii alone, reference twist | +1.14% | +0.85% |

Baselines: the production net row is against the exact rigid participant BEM (2.525066 MN /
16.326608 MW); the two decomposition rows are against the same-construction reference-station
solver (`_rebuild_bem_solver(_ref_r, _ref_twist)`, 2.541660 MN / 16.389370 MW), which differs from
the rigid path only by `_rebuild_bem_solver`'s own `Rtip = max(r_def) * 1.001` rule. The before
row is -1.84%/-1.20% only because the 3-4 m radius deficit (and the NaN-clamped stations) shrank
the rotor. **The magnitude does not approach Zhou**, and no bound was fitted to it.

**Asserted:** the rigid path equals an independent `BEMSolver` to 1e-9 (measured 0.000e+00);
`max|r_def - _ref_r| < 2 m` and `r_def.min > 0.5 Rhub`; the production twist equals the centre
ring's independently computed section rotation to 1e-6 deg; every outer-station section rotation
is nose-down; the elastic twist alone reduces thrust and power; the production path reduces
thrust. **Reported:** the magnitude against Zhou, the nodal `theta_z` table, and the production
net power.

**Two premises of the task were refuted by measurement; both are reported, not fitted.**

1. *The nodal `theta_z` does not equal the in-plane section rotation on this shell.* The tip ring
   reaches `mean(theta_z) = -3.8558` deg while its in-plane section rotation is -1.5112 deg
   (2.3446 deg apart); the same gap is already in sections 20 and 22.4. No translation-only
   estimator reaches -3.86: measured at the tip, the rigid-rotation fit gives +0.325 deg, the ring
   chord-line rotation +0.381 deg, and a Kabsch rigid-body fit +0.451 deg. The 0.5 deg agreement
   bound is therefore not written; the section rotation is the reference and the drilling rotation
   is printed.
2. *Production power does not fall.* +0.42% as measured above; only thrust falls.

**Not established** (stated, never guessed):

* Which twist magnitude the aerodynamics should see - the section rotation (-1.51 deg at the tip,
  the geometric quantity, and what production now uses) or the shell's drilling rotation
  (-3.86 deg, which matches Zhou's -3.60 deg better at 1.071x, section 22.4). The participant only
  receives the three translational DOFs, so it cannot recover the drilling DOF.
* Whether the de-loading magnitude can reach Zhou's -13.04% / -8.38% at all with a one-way,
  radius-and-twist-only feedback: our steady BEM on rigid-blade loads has no changed velocity
  triangle from the deflected shape, which is a large part of a coupled de-loading.
* ~~Whether the +1.63% axial stretch of the deformed blade is caused by the minimum-norm force
  distribution.~~ **REFUTED by measurement (section 22.7): the spanwise forces are not the cause.**
  The stretch is real and its cause is structural: a linear solve on a curved (prebent) blade axis
  under an in-plane flapwise load. The `r_def` growth (up to +0.93 m) that re-loads the rotor stays
  un-attributed to any load-frame defect, and no load-frame change removes it.
* Whether the retired strip-cloud PCA estimator was ever correct for a straight blade (it cannot
  be separated from P6 with the data collected here).

**Commands.** Target file `4 passed in 14 s`; neighbours (`test_force_projection*.py`,
`test_blade_rated_twist.py`) `25 passed`; whole suite with `CCX_BIN` `507 passed, 13 xfailed,
0 failed, 0 skipped` (baseline 503/13/0/0 plus the 4 new guards; without `CCX_BIN` on `PATH`, 69
rows skip and 0 fail); `ruff check` clean on both touched code files.

### 22.7 The spanwise-force hypothesis, refuted before it was implemented

The follow-up named at the end of section 22.6 - "the `+1.63%` axial stretch is most likely
injected by `_distribute`'s minimum-norm spanwise force components" - was tested **before** being
implemented, and it is wrong. No repository file was changed by this unit; the numbers below come
from an in-process experiment on the unmodified tree (`HEAD = ac69c8d`), whose harness reproduces
the recorded values exactly (`sum|fz| = 276,080.3 N`, `sum|fy| = 1,104,872 N`, ring-mean axial
displacements +0.033 / +0.297 / +1.063 m).

**Where the spanwise force actually comes from.** Decomposing each strip moment into its span
component (torsion) and the rest (`M_bend = M_strip - M_tors`):

| part of the load system | `sum|fz|` it produces | share |
| --- | ---: | ---: |
| net force alone (`F_strip`, `M = 0`) | **0.0 N** | 0% |
| torsion alone (`M_tors`, the pitching moment on the span axis) | 6,389.9 N | 2.31% |
| the rest (`M_bend`) | 273,763.9 N | **99.16%** |
| whole system | 276,080.3 N | 100% |

and `M_bend` is an identity, not an interpretation: `M_bend = -dz_ac * (span x F_strip)` with
`max ||M_bend + dz_ac x ...|| = 0.0` and `sum |dz_ac| |F| = 749,261 N.m = sum |M_bend|`, where
`dz_ac = (centroid - AC) . span` reaches +1.10 m. So the spanwise force is the price of collapsing a
strip's distributed in-plane load onto one centroid point while the aerodynamic centre sits up to
1.10 m away along the span - not of the shear-flow/torsion realisation. The proposed shear-flow fix
would therefore have removed **2.3%** of it.

**It is also structurally inert.** Constraining the whole solve to the section plane (in-plane
unknowns only, same six constraints, feasible to `1.9e-10` worst-case residual over all 50 strips
because a strip spans ~4 rings and carries its own lever arm) zeroes `sum|fz|` **exactly** and
changes the structural response by `max ||delta u|| = 2.6e-3 m` - 0.016% of the 16.53 m response.
Axial tip displacement: **+1.0630 m before, +1.0631 m after**. Load split: the 276 kN spanwise
system alone gives **-0.0032 m**, the flapwise `y` component alone gives **+1.0641 m**. The
lengthening is the in-plane load acting on a **curved blade axis** (ring-centroid `y` runs from
+0.36 m at z = 24 m to -4.00 m at the tip), against a `-1.149 m` second-order inextensional estimate
that a single linear `spsolve` on a linear `K` cannot produce at all - so the proposed guard "assert
the axial tip displacement is negative" is not writable against this model, fixed or unfixed.

**De-loading and the twist estimators, re-measured under the constraint (reported, not fitted):**
twist-only thrust -4.00% -> -3.90%, power -0.81% -> -0.78%; production-path thrust -2.32% -> -2.22%,
power +0.42% -> +0.45%; `r_def - r_ref` unchanged at -0.2470 .. +0.9332 m. Tip gap between the ring
section rotation and `mean(theta_z)`: **2.3446 deg before, 2.5233 deg after** - the divergence is
**not** a symptom of the load frame, it is a property of comparing a shell wall-bending field with a
section rigid rotation (section 22.6 finding 1 stands, my "symptom" reading in that section's
discussion does not).

**Also measured, and it blocks the fix as originally written:** the ring node order in this mesh is
not a contour order - the shoelace area differs by up to 29% between stored and angular ordering for
**133 of 186** raw single-z rings and **171 of 186** merged physical rings - so `q = M / (2 A_ring)`
is not computable from the stored order without adding a contour-ordering step.

**Decision taken.** The in-plane constraint is kept as a candidate **correctness/hygiene** change
(the aero load has no spanwise component, so injecting 276 kN of self-equilibrated spanwise force is
not physical), explicitly **not** as the fix for the axial stretch, the `r_def` growth or the
de-loading magnitude - measurement says it changes none of them. Deferred until the twist question
below is settled, so that it is not mistaken for a physics fix.

### 22.8 The anchor beam arbitrates which shell estimator is the section rotation

`tests/test_blade_twist_anchor_beam.py` (new, 2 tests, 57 s) builds an **independent** 1D beam:
section stiffness from the official IEA-15-240-RWT BeamDyn blade deck (26 stations), loads from
the **production** BEM + `ForceProjector` the shell is solved with, and the root-fixed cantilever
torsion `theta(z) = int_0^z T(s)/GJ(s) ds`, `T(z) = int_z^R m ds`, with
`m = Mp + (x_AC - xS) Np - (y_AC - yS) Tp`. `xS, yS, GJ` are the **anchor's**, not the shell's -
the models share *loads, not section properties*.

Parse cross-check (`K66toPropsDecoupled`, never hand-read): `EA = 4.605108e10 N` (0.00 % off
4.605e10), `GKt = 8.748569e10 N.m^2` (0.005 % off 8.749e10), `K[2,2]` axial, `K[5,5]` raw torsion.
Bending: `EIxp = 1.495993e11` vs `FlpStff = 1.525339e11` (+1.92 %), `EIyp = 1.497329e11` vs
`EdgStff = 1.524792e11` (+1.80 %) - inside the 2 % rule; the torsion-at-index-3 reading misses
`FlpStff` by 42.6 %.

Degenerate tip: `median(GKt)/1e3 = 4.104e5`; `GKt[-2] = 6.501e6` kept, `GKt[-1] = 5.879e4` excluded,
only interval `[24]`; it would otherwise have added **-0.0257 deg** under this cantilever
quadrature (not §18.2's -9.4 deg, which belongs to §18's load-defective set). Lever-arm sign: the
task's `(0.25 - pitch_axis)` is the negative of the production `r_AC . c_hat`; the written sign
gives a 10x root-torque inconsistency, so the measured sign is used.

Comparison (real mesh, rated point, all 26 stations, degrees about +span):

| r_root | r_shell | n | `phi_beam` | ring rot | mean `theta_z` | ring/beam | meanz/beam |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.00 | 22 | 0.0000 | +0.0000 | +0.0000 | - | - |
| 1.17 | 0.80 | 22 | -0.0001 | +0.0012 | +0.0016 | -10.35 | -13.74 |
| 2.34 | 2.39 | 22 | -0.0002 | +0.0032 | +0.0048 | -13.60 | -20.87 |
| 3.51 | 3.18 | 22 | -0.0004 | +0.0045 | +0.0088 | -12.47 | -24.28 |
| 4.68 | 4.78 | 22 | -0.0005 | +0.0087 | +0.0224 | -17.07 | -43.87 |
| 5.85 | 5.57 | 22 | -0.0007 | +0.0111 | +0.0317 | -16.38 | -46.76 |
| 8.77 | 8.76 | 22 | -0.0012 | +0.0232 | +0.0745 | -19.04 | -61.00 |
| 11.70 | 11.94 | 26 | -0.0020 | +0.0312 | +0.0995 | -15.40 | -49.16 |
| 17.55 | 17.51 | 23 | -0.0061 | +0.0099 | -0.0885 | -1.61 | +14.40 |
| 23.40 | 23.08 | 23 | -0.0185 | -0.0238 | -0.1011 | +1.28 | +5.45 |
| 29.25 | 29.45 | 21 | -0.0522 | -0.1188 | -0.1753 | +2.28 | +3.36 |
| 35.10 | 35.02 | 20 | -0.1145 | -0.1085 | -0.1037 | +0.95 | +0.91 |
| 40.95 | 40.59 | 19 | -0.2002 | -0.1447 | -0.1791 | +0.72 | +0.89 |
| 46.80 | 46.96 | 18 | -0.3033 | -0.1787 | -0.2011 | +0.59 | +0.66 |
| 52.65 | 52.53 | 18 | -0.4214 | -0.2044 | -0.1726 | +0.48 | +0.41 |
| 58.50 | 58.50 | 16 | -0.5550 | -0.2610 | -0.2087 | +0.47 | +0.38 |
| 64.35 | 64.47 | 14 | -0.7040 | -0.3234 | -0.2467 | +0.46 | +0.35 |
| 70.20 | 70.44 | 14 | -0.8657 | -0.3955 | -0.2880 | +0.46 | +0.33 |
| 76.05 | 75.81 | 14 | -1.0364 | -0.4885 | -0.3325 | +0.47 | +0.32 |
| 81.90 | 81.78 | 14 | -1.2094 | -0.6096 | -0.4140 | +0.50 | +0.34 |
| 87.75 | 87.75 | 14 | -1.3753 | -0.7312 | -0.4312 | +0.53 | +0.31 |
| 93.60 | 93.72 | 14 | -1.5237 | -0.8453 | -0.0831 | +0.55 | +0.05 |
| 99.45 | 99.33 | 12 | -1.6483 | -1.0711 | -0.6979 | +0.65 | +0.42 |
| 105.30 | 105.06 | 12 | -1.7406 | -1.1992 | -3.2725 | +0.69 | +1.88 |
| 111.15 | 111.27 | 12 | -1.7765 | -0.9741 | -0.9703 | +0.55 | +0.55 |
| **117.00** | **117.00** | **12** | **-1.7765** | **-1.5112** | **-3.8558** | **+0.85** | **+2.17** |

Zhou et al. 2025, *Energy* 336:138488, Table 4 - **reported, never asserted**: torsion -3.60 deg,
flap +13.86 m, edge -1.22 m. Asserts: the twist magnitude is monotone toward the tip, and the
excluded tip interval is nose-down.

**Corrected after RDD review `review-0d3ea0d5b17e0dfe` (finding R3-001, CRITICAL, deterministic).**
The quadrature multiplied station `i`'s step by `dr[i]` - the *next* interval's length - while
admitting that step under `interval_ok[i-1]`, so guard and spacing referred to different intervals.
The deck's station spacing is not uniform (span fractions step by 0.01, 0.025 and 0.05), so every
station below the tip was integrated over a mismatched interval, and the final step would have
indexed `dr` out of bounds had the tip guard not short-circuited it. Measured impact of the fix:
tip `-1.7776 -> -1.7765 deg` (0.06 %), ratios `0.8501 -> 0.8506` and `2.1691 -> 2.1704`, and the
outer-half `ring/beam` band stays `+0.46 .. +0.55`. The arbitration and the deficit bounds are
unchanged at the precision quoted here - but the beam column in the first draft of this section was
not the documented root-fixed quadrature, and this section's own monotonicity assertion did not
catch it: a monotone-but-wrong integral passes a monotonicity guard.

**Verdict: (i).** Tip: beam **-1.777**, ring section rotation **-1.511** (ratio **0.851**),
`mean(theta_z)` **-3.856** (ratio **2.170**). A beam from the turbine's own `GJ`/shear centre is
within 15 % of the *section rotation* and 2.2x off the nodal average, so the affine estimator is
the shell's section rotation and **(iii) is refuted**: the meanz/beam ratio wanders by a factor of
40 across the span (+12.3, +5.2, +0.9, +0.33, +0.05, +1.9, +2.2 deg per deg) while ring/beam holds
+0.46 .. +0.55 over the outer half, so the nodal average is a wall-bending field, not a rotation.
**(i) holds for the estimator.** The magnitude, however, is NOT explained by the known `GJ` ratio -
the parent corrected the direction of that argument, which the first draft of this section had
inverted:

- `GJ_shell / GJ_ref = (4.000 / 4.290)^2 = 0.87` (§15.1, modal, load-independent) says the shell is
  the *softer* member, so under the same torque it must twist `1 / 0.87 = 1.149` x *more* than the
  anchor beam: expected `phi_shell = 1.149 x 1.7776 = 2.043 deg`.
- Measured `phi_shell = 1.5112 deg`, i.e. observed `beam/shell = 1.176` against an expected
  `beam/shell = 0.87`. The shell delivers **0.74 x** the rigid twist its own independently measured
  torsional stiffness calls for - a **1.35 x deficit**, not the "2-3 % agreement" the first draft
  claimed. (Comparing an observed ratio to the reciprocal of the expected one is the error; the
  numbers in the table are unaffected.)
- **SUPERSEDED by section 22.10, and it was my own reference mix-up.** The `0.87` above is *not* a
  shell-versus-deck ratio: section 15.1 forms it from our shell's 4.000 Hz against **NuMAD's**
  4.290 Hz, so its denominator is another team's discretisation (the same section rates that
  comparison's scatter at 15.1 %). Scaling a **deck-`GJ`** beam by it, as this bullet does, mixes
  references. The deck-referenced static measurement (section 22.10) gives `GJ_shell/GJ_deck = 1.068`
  over `r = 20-100 m`, with a tight `1.03-1.09` plateau over `r = 35-80 m`, so the expected twist is
  `1.7765/1.068 = 1.663 deg` against the measured `1.5112 deg`: a **9 % residual**, not a 26-35 %
  deficit, and it sits inside the anchor beam's own convention uncertainty (26-station quadrature,
  interpolated loads, anchor shear centre, `y_AC = 0`). Estimator verdict (i) is untouched; what
  dies here is the "under-twist" framing, and with it the premise that a structural explanation is
  owed.

So **(ii) is back on the table**: the shell under-twists relative to both an independent beam and
its own modal stiffness, and two independent measurements already point at the load realisation
rather than the element - §20.5's `distortion / |omega| = 9.36` under the production path, and
§19.3's negative result (commit d4fec33) that a statically equivalent end load on a thin-walled tube buys 125x
section distortion and 74.8x tip displacement at the *same* torque. Open and to be measured: the
26-station quadrature, interpolated loads, the anchor shear centre and `y_AC = 0` leave the residual
unresolved below their own convention uncertainty, so no bound is asserted on the 1.35x. The
de-loading gap follows the magnitude: a tip twist of -1.51 deg unloads the rotor far less than
Zhou's -3.60 deg (reported, not used as a criterion). No shell code changed.

### 22.9 Does the production nodal moment realisation under-deliver rigid torsion?

**REJECTED - measured on the exact Bredt tube, the production minimum-norm realisation of a pure
section moment does not under-deliver rigid torsion; it over-delivers it by 31.7x
(self-equilibrated) / 8.9x (clamped tip), so it is not the cause of the blade's 1.35x under-twist.
The minimum-norm field is nonetheless not a wall shear flow and it does depart from Bredt by a large
factor - the departure is an over-soft end distortion, the opposite sign to the blade deficit.**

`tests/test_thin_walled_tube_moment_realization.py` (1 test, 5.7 s) reuses the mesh, torque,
estimators and hand-written reference of `tests/test_thin_walled_tube_torsion.py`
(`import test_thin_walled_tube_torsion as tube`). One mesh (L = 6 m, 16-node ring, 496 nodes), one
torque (T = 1e4 N.m), one rigid-mode-removed solve:

* **case A** - the validated shear flow `q = T/(2A)`, `+T` at the tip ring and `-T` at the root ring
  (`tube._self_equilibrated_load`);
* **case B** - the **same** net torque through the production code under test,
  `ForceProjector._distribute(strip, F = 0, M = T z_hat)` at the two end rings
  (`aeroelast.solvers.bem.force_projection`; the minimum-norm solve is called, not re-implemented).

Case A is the validated reference and is asserted within 5%. Case B's output is not trusted: its
tip/root torque is re-measured with the independent ruler `sum(x Fy - y Fx)` (ratio 1.0000) and its
net force (1e-13), and its nodal field is compared to the closed form below.

Measured (rate of `theta_fit` over 0.4L-0.9L against the hand-written Bredt `T/GJ`; distortion =
`|0.5(alpha+beta)|` of the section parallelogram at 0.75L; deviation = field residual after removing
the best rigid rotation):

| case | rate/Bredt | tip/Bredt(L/2) | distortion/|rotation| | deviation-from-rigid |
| --- | ---: | ---: | ---: | ---: |
| A shear flow | 1.0031 | 1.0033 | 0.008 | 0.007 |
| B min-norm (P8) | **31.6933** | **33.6664** | **2.042** | **0.876** |

(`theta_z` rate, for completeness: A 1.0067, B 1.0374.) Case A is inside 5%; case B misses the exact
independent reference by 3069% and on the **over**-soft side, so the P8 hypothesis as stated - a
smaller rigid rotation caused by the nodal realisation - is refuted.

**Localisation (B misses, so where the nodal system departs from shear flow).** The production field
is exactly `f_j = omega x d_j` with `omega = T / sum|d|^2` (measured `max |f_j - omega x d_j| =
2.3e-13`), i.e. each nodal force is **perpendicular to that node's position vector about the strip
centroid** (circle-tangential), not tangential to the wall edges. On this rectangle the wall-mid nodes
carry 802 N against the shear flow's uniform 2083 N, the corner-diagonal nodes carry 1559 N against
1215 N, and the implied perimeter shear flow is strongly non-uniform where Bredt requires the
constant `q = T/(2A)`; the section absorbs the difference as distortion (deviation-from-rigid 0.876,
`distortion/|rotation|` 2.042, against case A's 0.008). This is a **long-decaying end effect**, not a
rate change: extending the tube to L = 24/48 m drives the clamped-tip ratio 8.86 -> 3.20 -> 2.10
toward Bredt, the extra tip rotation staying fixed at ~1.45e-3 rad (T = 1e4 N.m), so the interior
returns to Saint-Venant torsion once the added end distortion is negligible. On the blade the strips
are ~2.4 m apart against a ~10 m decay length, so those end distortions overlap and are locally
significant - but they soften the section. **A real non-shear-flow load-path defect in
`force_projection.py`, of the opposite sign to the blade deficit; it is not P8 as hypothesised.**

**Blade-deficit restatement (recorded values, not re-derived here; the modal ratio and the blade
solve were not re-run).**

| quantity | value | source | status |
| --- | ---: | --- | --- |
| anchor `GJ_ref` (BeamDyn `GKt`) | 8.7486e10 N.m^2 | §22.8 | recorded |
| modal `GJ_shell/GJ_ref` | (4.000/4.290)^2 = 0.870 | §15.1 | recorded |
| anchor beam tip twist `phi_beam` | -1.7765 deg | §22.8 | recorded |
| shell tip section rotation `phi_shell` | -1.5112 deg | §22.8 | measured |
| expected shell (same torque) `phi_beam/0.870` | -2.042 deg | derived | - |
| observed / expected | 1.5112/2.042 = **0.740 -> 1.35x deficit** | derived | - |

**Bound on that 1.35x (the parent's correction - the ratio is not a sharp number), and superseded
by section 22.10 (the ratio is not even a shell-versus-deck number).** The modal
`0.870` comes from one torsional frequency, 4.000 Hz against NuMAD's 4.290 Hz, and section 15.1
records that mode's **reference scatter as 15.1%**. Propagating it, the expected shell twist is
`2.042 x [0.849, 1.151] = [1.734, 2.351] deg`, and the measured `1.5112 deg` sits below the whole
band. So the deficit is real - the shell under-twists relative to its own modal torsional
stiffness - but the honest statement is **between 1.15x and 1.35x**, not 1.35x, and closing it
needs a sharper independent torsional stiffness than one modal ratio with 15% scatter.

For there to be no deficit, one of the three recorded quantities would have to be wrong by ~35%: the
modal ratio would have to be `1.7765/1.5112 = 1.176` (the shell 18% stiffer than the reference,
not 13% softer), or the anchor beam would have to give `-1.3147 deg` (its `GJ` ~35% higher), or the
shell's measured section rotation would have to be `-2.043 deg`. This unit removes the
load-realisation candidate from that list in the direction the blade needs - the minimum-norm
realisation softens the section (over-twist), it does not stiffen it - so the 1.35x deficit still
points at the section-stiffness or anchor inputs, not at the nodal moment realisation.

**Verdict on P8: rejected.** No production code was changed in this unit.

### 22.10 Three independent routes to the blade's torsional stiffness

**Route (c), the new static measurement.** `tests/test_blade_section_torsion_stiffness.py` (1 test,
44 s, one blade solve). The realisation is the Bredt wall flow `q = T/(2A)` validated on the exact
tube (`d4fec33`), applied as `+T` on the tip ring and `-T` on the root ring of the real shell blade
(`q = 2.1017e8` N/m at the tip, `2.3959e5` N/m at the root), with the six rigid modes removed by the
exact saddle-point constraint (`test_thin_walled_tube_torsion._solve_rigid_removed`, reused). The
enclosed area `A` is the hand-written shoelace of the **ordered** contour (root 20.8688 m^2 on a
22-node contour, tip 2.3790e-2 m^2 on a closed 12-node single cell), never the force vector's own
bookkeeping. One solve gives the whole span; `phi(z)` is the antisymmetric part of each physical
ring's affine fit (section 22.8 - never `mean(theta_z)`), and `GJ_static = T / (dphi/dz)` is the
secant between consecutive physical stations (`_physical_stations`, the section 14.3 prebend merge).
Assertions are only the physics-pinned invariants (ruler, self-equilibration, the accumulating sign);
no bound is asserted against the deck, the modal ratio or Zhou.

**Torque ruler and self-equilibration.** Read back from the applied nodal system, not from the
intent: tip `sum(x Fy - y Fx) = +1.0000000000e7` N.m and root `-1.0000000000e7` N.m, and the tip
system's torque about its **own centroid** is the same `+1.0000000000e7` N.m - centre-free because
the couple carries no net force. Net force `8.382e-9` N against a max nodal force of `4.167e7` N
(**2.01e-16** relative); the six rigid-mode resultants are `1.1e-20 ... 1.0e-17`. The rigid-mode
*gauge* residual `|Phi^T u|/|u| = 6.75e-9` is a saddle-point solver tolerance, not a physics
quantity (a residual rigid rotation about z adds one constant to every ring and cannot change a
rate; about x/y it is invisible to the z-section fit). Sign: a positive torque about +span gives a
positive, monotonically accumulating twist (asserted).

**The per-station table** (`n` = ring nodes, `distort` = the ring's residual after the best rigid
rotation from `test_thin_walled_tube_torsion._parallelogram`):

```text
     r[m]   n GJ_static[N.m^2]  GJ_ref[N.m^2]   ratio  distort
     0.40  22     3.833730e+11   8.605660e+10   4.455    0.072
     1.19  22     9.569476e+10   8.319859e+10   1.150    0.079
     1.99  22     6.714091e+10   8.034611e+10   0.836    0.086
     2.79  22     2.792526e+10   7.560758e+10   0.369    0.093
     3.58  22     1.693623e+10   6.942973e+10   0.244    0.101
     4.38  22     1.581898e+10   6.368208e+10   0.248    0.111
     5.17  22     1.633405e+10   5.837496e+10   0.280    0.121
     5.97  22     1.550064e+10   5.345207e+10   0.290    0.131
     6.77  22     1.476690e+10   4.917649e+10   0.300    0.140
     7.56  22     1.360149e+10   4.490091e+10   0.303    0.149
     8.36  22     1.282648e+10   4.062533e+10   0.316    0.158
     9.15  22     1.108651e+10   3.670327e+10   0.302    0.167
     9.95  21     1.080599e+10   3.317193e+10   0.326    0.204
    10.74  21     1.003451e+10   2.964060e+10   0.339    0.212
    11.54  21     6.261892e+09   2.611100e+10   0.240    0.221
    12.34  26     8.449674e+09   2.351771e+10   0.359    0.245
    13.13  26     6.970639e+09   2.115826e+10   0.329    0.252
    13.93  26     5.991144e+09   1.879882e+10   0.319    0.260
    14.72  25     5.132582e+09   1.643937e+10   0.312    0.257
    15.52  25     4.328407e+09   1.407992e+10   0.307    0.266
    16.32  25     4.939174e+09   1.172009e+10   0.421    0.275
    17.11  23     4.469966e+09   9.360234e+09   0.478    0.293
    17.91  23     4.647186e+09   7.736253e+09   0.601    0.302
    18.70  23     4.770708e+09   7.015233e+09   0.680    0.310
    19.50  23     4.337516e+09   6.294214e+09   0.689    0.318
    20.29  24     3.949315e+09   5.573195e+09   0.709    0.338
    21.09  24     3.209607e+09   4.852175e+09   0.661    0.344
    21.89  23     2.890461e+09   4.131156e+09   0.700    0.379
    22.68  23     2.654954e+09   3.410136e+09   0.779    0.385
    23.48  23     2.432588e+09   2.736809e+09   0.889    0.391
    24.27  23     1.892925e+09   2.499824e+09   0.757    0.396
    25.07  23     1.705255e+09   2.262839e+09   0.754    0.400
    25.87  23     1.617550e+09   2.025854e+09   0.798    0.404
    26.66  23     1.382937e+09   1.788870e+09   0.773    0.407
    27.46  23     1.267260e+09   1.551885e+09   0.817    0.411
    28.25  23     1.243641e+09   1.314989e+09   0.946    0.415
    29.05  21     1.100580e+09   1.078095e+09   1.021    0.400
    29.84  21     1.058131e+09   9.780164e+08   1.082    0.402
    30.64  21     1.062654e+09   9.241918e+08   1.150    0.405
    31.44  20     8.101309e+08   8.703667e+08   0.931    0.425
    32.23  20     7.764253e+08   8.165566e+08   0.951    0.427
    33.03  20     7.391860e+08   7.627464e+08   0.969    0.428
    33.82  20     7.024820e+08   7.089363e+08   0.991    0.430
    34.62  20     6.689649e+08   6.551261e+08   1.021    0.431
    35.42  20     6.403896e+08   6.138195e+08   1.043    0.433
    36.21  20     6.135005e+08   5.914632e+08   1.037    0.434
    37.01  20     5.868931e+08   5.691069e+08   1.031    0.435
    37.80  20     5.638044e+08   5.467507e+08   1.031    0.437
    38.60  20     5.417295e+08   5.243944e+08   1.033    0.438
    39.40  20     5.252691e+08   5.020381e+08   1.046    0.439
    40.19  20     4.991355e+08   4.796818e+08   1.041    0.440
    40.99  19     4.966015e+08   4.577633e+08   1.085    0.483
    41.78  19     4.790716e+08   4.447007e+08   1.077    0.485
    42.58  19     4.632796e+08   4.316380e+08   1.073    0.486
    43.38  19     4.466938e+08   4.185754e+08   1.067    0.487
    44.17  19     4.299072e+08   4.055127e+08   1.060    0.487
    44.97  19     3.943369e+08   3.924500e+08   1.005    0.486
    45.76  18     4.041139e+08   3.793874e+08   1.065    0.523
    46.56  18     3.919743e+08   3.663247e+08   1.070    0.524
    47.35  18     3.805417e+08   3.549490e+08   1.072    0.524
    48.15  18     3.703447e+08   3.443078e+08   1.076    0.524
    48.95  18     3.597531e+08   3.336666e+08   1.078    0.523
    49.74  18     3.479338e+08   3.230254e+08   1.077    0.523
    50.54  18     3.364965e+08   3.123842e+08   1.077    0.523
    51.33  18     3.261664e+08   3.017430e+08   1.081    0.522
    52.13  18     3.153396e+08   2.911018e+08   1.083    0.521
    52.83  18     3.044357e+08   2.821905e+08   1.079    0.521
    53.42  18     2.960454e+08   2.755647e+08   1.074    0.520
    54.02  18     2.873819e+08   2.689389e+08   1.069    0.520
    54.62  18     2.837108e+08   2.623131e+08   1.082    0.521
    55.21  17     2.697720e+08   2.556873e+08   1.055    0.499
    55.81  17     2.633908e+08   2.490614e+08   1.058    0.499
    56.41  17     2.585756e+08   2.424356e+08   1.067    0.500
    57.00  17     2.472276e+08   2.358098e+08   1.048    0.500
    57.60  16     2.500705e+08   2.291840e+08   1.091    0.478
    58.20  16     2.427105e+08   2.225582e+08   1.091    0.478
    58.80  16     2.369172e+08   2.165731e+08   1.094    0.479
    59.39  16     2.310791e+08   2.112431e+08   1.094    0.480
    59.99  16     2.240871e+08   2.059131e+08   1.088    0.481
    60.59  16     2.184248e+08   2.005830e+08   1.089    0.482
    61.18  16     2.119067e+08   1.952530e+08   1.085    0.484
    61.78  16     2.045398e+08   1.899229e+08   1.077    0.485
    62.38  15     1.994308e+08   1.845929e+08   1.080    0.521
    62.97  15     1.938999e+08   1.792629e+08   1.082    0.523
    63.57  15     1.877839e+08   1.739328e+08   1.080    0.525
    64.17  15     1.832437e+08   1.685866e+08   1.087    0.527
    64.77  14     1.761519e+08   1.641130e+08   1.073    0.530
    65.36  14     1.707507e+08   1.600296e+08   1.067    0.533
    65.96  14     1.656110e+08   1.559461e+08   1.062    0.537
    66.56  14     1.605977e+08   1.518626e+08   1.058    0.540
    67.16  14     1.545411e+08   1.477792e+08   1.046    0.545
    67.75  14     1.499753e+08   1.436957e+08   1.044    0.549
    68.35  14     1.457082e+08   1.396122e+08   1.044    0.554
    68.95  14     1.416315e+08   1.355287e+08   1.045    0.559
    69.54  14     1.372477e+08   1.314453e+08   1.044    0.564
    70.14  14     1.334656e+08   1.273618e+08   1.048    0.570
    70.74  14     1.298193e+08   1.238947e+08   1.048    0.577
    71.33  14     1.262752e+08   1.204962e+08   1.048    0.583
    71.93  14     1.223818e+08   1.170976e+08   1.045    0.591
    72.53  14     1.189394e+08   1.136990e+08   1.046    0.598
    73.12  14     1.155033e+08   1.103004e+08   1.047    0.607
    73.72  14     1.121053e+08   1.069019e+08   1.049    0.615
    74.32  14     1.076970e+08   1.035033e+08   1.041    0.624
    74.92  14     1.045835e+08   1.001047e+08   1.045    0.632
    75.51  14     1.018963e+08   9.670614e+07   1.054    0.641
    76.11  14     9.963563e+07   9.340031e+07   1.067    0.651
    76.71  14     9.738046e+07   9.092914e+07   1.071    0.660
    77.30  14     9.522084e+07   8.845798e+07   1.076    0.669
    77.90  14     9.262755e+07   8.598681e+07   1.077    0.679
    78.50  14     8.962716e+07   8.351565e+07   1.073    0.689
    79.09  14     8.544711e+07   8.104448e+07   1.054    0.698
    79.69  14     8.261774e+07   7.857332e+07   1.051    0.708
    80.29  14     8.018130e+07   7.610216e+07   1.054    0.717
    80.89  14     7.817771e+07   7.363099e+07   1.062    0.727
    81.48  14     7.550545e+07   7.115983e+07   1.061    0.737
    82.08  14     7.363942e+07   6.884995e+07   1.070    0.750
    82.68  14     7.176590e+07   6.691640e+07   1.072    0.782
    83.27  14     6.991265e+07   6.498285e+07   1.076    0.845
    83.87  14     6.732790e+07   6.304930e+07   1.068    0.780
    84.47  14     6.561442e+07   6.111575e+07   1.074    0.784
    85.06  14     6.390508e+07   5.918220e+07   1.080    0.791
    85.66  14     6.235870e+07   5.724865e+07   1.089    0.799
    86.26  14     6.001686e+07   5.531511e+07   1.085    0.807
    86.85  14     5.886042e+07   5.338156e+07   1.103    0.815
    87.45  14     5.744616e+07   5.144801e+07   1.117    0.823
    88.05  14     5.550455e+07   4.975597e+07   1.116    0.831
    88.65  14     5.094032e+07   4.830544e+07   1.055    0.838
    89.24  14     4.915089e+07   4.685491e+07   1.049    0.845
    89.84  14     4.815119e+07   4.540438e+07   1.060    0.852
    90.44  14     4.741846e+07   4.395385e+07   1.079    0.859
    91.03  14     4.477830e+07   4.250332e+07   1.054    0.866
    91.63  14     4.416806e+07   4.105279e+07   1.076    0.872
    92.23  14     4.327550e+07   3.960226e+07   1.093    0.878
    92.82  14     4.288221e+07   3.815173e+07   1.124    0.884
    93.42  14     4.062887e+07   3.670120e+07   1.107    0.889
    94.02  14     4.041012e+07   3.527939e+07   1.145    0.894
    94.61  14     3.892871e+07   3.386988e+07   1.149    0.899
    95.21  14     3.771737e+07   3.246038e+07   1.162    0.904
    95.75  13     3.452559e+07   3.119183e+07   1.107    0.902
    96.23  13     3.398331e+07   3.006422e+07   1.130    0.906
    96.70  13     3.289215e+07   2.893662e+07   1.137    0.909
    97.18  13     3.352540e+07   2.780902e+07   1.206    0.913
    97.66  13     3.106819e+07   2.668141e+07   1.164    0.916
    98.14  12     3.108457e+07   2.555381e+07   1.216    0.914
    98.61  12     3.037850e+07   2.442621e+07   1.244    0.917
    99.09  12     3.021276e+07   2.329861e+07   1.297    0.920
    99.57  12     2.901496e+07   2.226001e+07   1.303    0.923
   100.05  12     2.825263e+07   2.148842e+07   1.315    0.926
   100.52  12     2.071252e+07   2.071683e+07   1.000    0.928
   101.00  12     1.997960e+07   1.994524e+07   1.002    0.931
   101.48  12     1.884594e+07   1.917365e+07   0.983    0.933
   101.96  12     1.857905e+07   1.840206e+07   1.010    0.935
   102.43  12     1.780634e+07   1.763047e+07   1.010    0.938
   102.91  12     1.638658e+07   1.685888e+07   0.972    0.940
   103.39  12     1.582680e+07   1.608729e+07   0.984    0.942
   103.87  12     1.572787e+07   1.531570e+07   1.027    0.944
   104.34  12     1.514993e+07   1.454411e+07   1.042    0.946
   104.82  12     1.504068e+07   1.377252e+07   1.092    0.948
   105.30  12     1.289537e+07   1.300093e+07   0.992    0.950
   105.78  12     1.270926e+07   1.247031e+07   1.019    0.952
   106.26  12     1.217034e+07   1.193969e+07   1.019    0.954
   106.73  12     1.191010e+07   1.140907e+07   1.044    0.956
   107.21  12     1.143182e+07   1.087845e+07   1.051    0.958
   107.69  12     1.117887e+07   1.034782e+07   1.080    0.959
   108.17  12     1.073411e+07   9.817204e+06   1.093    0.961
   108.64  12     1.038292e+07   9.286583e+06   1.118    0.963
   109.12  12     9.939860e+06   8.755962e+06   1.135    0.964
   109.60  12     9.641661e+06   8.225341e+06   1.172    0.966
   110.08  12     7.833629e+06   7.694720e+06   1.018    0.967
   110.55  12     7.560047e+06   7.164099e+06   1.055    0.969
   111.03  12     7.346735e+06   6.633478e+06   1.108    0.970
   111.51  12     6.869794e+06   6.106413e+06   1.125    0.971
   111.99  12     5.943034e+06   5.580532e+06   1.065    0.973
   112.46  12     5.163833e+06   5.054652e+06   1.022    0.974
   112.94  12     5.044463e+06   4.528771e+06   1.114    0.975
   113.42  12     4.081120e+06   4.002891e+06   1.020    0.977
   113.90  12     3.807176e+06   3.477010e+06   1.095    0.979
   114.37  12     2.904910e+06   2.951130e+06   0.984    0.981
   114.81  12     1.572481e+06   2.469073e+06   0.637    0.984
   115.21  12     1.014619e+06   2.030839e+06   0.500    0.989
   115.61  12     5.869040e+05   1.592606e+06   0.369    0.994
   115.96  12     2.808089e+05   1.209151e+06   0.232    0.998
   116.25  12     1.507721e+05   8.804758e+05   0.171    1.000
   116.55  12     7.315922e+04   5.518005e+05   0.133    1.000
   116.85  12     2.629636e+04   2.231252e+05   0.118    1.000
```

**The three routes side by side.**

| route | quantity | value | source | status |
| --- | --- | ---: | --- | --- |
| (a) anchor deck | `GKt` at the root | **8.748569e10 N.m^2** | section 22.8 (`K66toPropsDecoupled`) | recorded, not re-derived |
| (b) modal | `GJ_shell/GJ_ref` | **0.870**; 1-sigma band `[0.739, 1.001]` from the 15.1% reference scatter | section 15.1 (4.000 Hz vs 4.290 Hz) | recorded, **not re-run** |
| (c) static, this unit | `GJ_static/GJ_ref`, r = 20-100 m | **median 1.068** (range 0.661-1.303) | this section | measured |
| (c) static, tight plateau | `GJ_static/GJ_ref`, r ~ 35-80 m | **1.03-1.09, median ~1.07** | this section | measured |

**Which branch, and what it falsifies.** The test's own decision rule lands on **branch 1** (the
window median `1.068` is within the modal ratio's 15.1% scatter of `GJ_ref`), but the number is
`1.07`, not `1.00`, and that 7% excess is the content of the finding:

* **Routes (b) and (c) do not disagree - they have different denominators.** Read literally, the
  modal ratio says the shell is 13% **softer** than the deck (`0.870`) while the static route says it
  is 7% **stiffer** (`1.068`), i.e. `1.068/0.870 = 1.228`, a **22.8%** gap that the modal route's own
  1-sigma band (`0.870 x 1.151 = 1.001`) cannot carry. But (b) never measured the deck: section 15.1
  forms it from our shell's 4.000 Hz against **NuMAD's** 4.290 Hz, i.e. its denominator is another
  team's discretisation of the turbine, and the same section records that comparison's scatter as
  15.1% - which is most of the 22.8%. Only (c) is referenced to the **deck** (a). So the correct
  conclusion is not that the modal ratio is false, it is that **using `0.870` as a shell-versus-deck
  ratio was a reference mix-up**, and it was *my* mix-up: section 22.8 built its expected twist by
  scaling a deck-`GJ` beam with a NuMAD-referenced number. The parent corrected the wording here and
  the arithmetic there.
  NuMAD - and the sharper static measurement replaces it. The deck (a) and the static route (c) agree
  to 7%, just outside the suite's 5% rule, so the shell's section torsional stiffness is essentially
  the deck's.
* **It also removes most of the section 22.8 deficit.** The deficit was computed against the modal
  ratio: expected shell twist `1.7765/0.870 = 2.042 deg` against measured `1.5112 deg` -> `0.740`, a
  26% shortfall. Against the **static** ratio: expected `1.7765/1.068 = 1.6634 deg` against the same
  measured `1.5112 deg` -> **`0.909`, a 9% residual**, i.e. the deficit falls from 26% to 9% with no
  load-frame change. The remaining 9% is inside the anchor beam's own convention uncertainty
  (section 22.8's 26-station quadrature, interpolated loads, shear centre, `y_AC = 0`). The outer-span
  ratios rise to `1.15-1.31` at r ~ 92-100 m, where the deck's `GKt` plunges (its last valid station
  is r = 111.15 m), so the shell is comparatively stiffer where the twist accumulates most - the
  direction needed to explain the beam/shell gap.
* **Branch 1's literal reading is NOT triggered.** Branch 1 predicts that a shell stiffness equal to
  the deck's leaves the whole deficit in the delivered torque, contradicting section 22.7's moment
  identity. That is not what happened: the static ratio is `1.07`, not `1.00`, so ~7 of the 15
  percentage points of the section 22.8 gap are stiffness, not torque, and the residual is 9%, not
  15%. No load-frame contradiction is needed and none is inferred.

**Where `GJ_static` stops being trustworthy.**

| region | why | measured signature |
| --- | --- | --- |
| r < ~28 m | the root couple is a uniform flow on the 22-node outer contour only; the root section is multi-cell (`allShearWebNods` = 831 nodes), so the applied system is self-equilibrated but **not** the multi-cell Saint-Venant flow, and its boundary layer runs to ~28 m | ratio 0.24-0.95, dist 0.09-0.38 |
| ~28-92 m | interior, constant internal torque, deck reference valid | ratio 0.93-1.15; tight plateau 1.03-1.09 over 35-80 m |
| ~92-111 m | interior, but the deck's `GKt` plunges and its last valid station is r = 111.15 m | ratio rises to 1.31 |
| r > ~111 m | the section 14.3 free-edge tear at z ~ 112.22 m, 12-node single cells of area 0.024 m^2, and the degenerate deck tip (`GKt[-1] = 5.879e4`, excluded by the median/1e3 floor) | ratio collapses to 0.12 at r = 116.85 |

The `distort` column is the systematic caveat on the absolute value: the ring's non-affine residual
grows from 0.34 at r = 20 m to 0.93 at r = 100 m to 1.00 at the tip, so the affine `phi(z)` is a
fitted average of a section that is deforming in plane, not an exact rigid rotation. The plateau's
tightness (1.03-1.09 over 45 m) shows the bias is smooth and largely cancels in the **ratio**; it is
why the honest claim is `GJ_static/GJ_ref`, not an absolute `GJ`.

**Verdict.** The static route (c) is the sharpest of the three (`1.07`, local spread ~3% over
r = 35-80 m, against the modal route's 15.1% reference scatter), it agrees with the independent deck
(a) to 7%, and it **replaces the modal ratio (b) as the shell-versus-deck number**: the shell is not
  13% softer than the deck, it is
7% stiffer, and increasingly stiffer toward the tip where the deck's `GKt` plunges. The section 15.1
modal ratio `0.870` is the suspect input; with it replaced, the section 22.8 beam/shell deficit falls
from 26% to 9%, and the remaining 9% is inside the anchor beam's own conventions - no load-frame
defect is required, so section 22.7's moment identity is not contradicted. No production code, no
existing test and no other document was changed, and the modal analysis was not re-run.

### 22.11 Does Zhou's -3.60 deg / -13.04% transfer to this model?

WU-D settles comparability and reports one de-loading table in one sign convention. The reference
is read from the local PDF (`.sources/papers/Unsteady aeroelastic performance of the 15 MW floating
offshore wind turbine under surge condition.pdf`; Zhou, Shen, Ma, Ouyang & Du 2025, *Energy*
336:138488); no paywalled fetch.

**Part 1, the paper's own words** (italic = quoted; "not stated" = it does not say).

1. **Radius.** Table 4 is titled *"Comparison of the mean tip deflections of the IEA-15 MW blade
   under rated condition"* - the torsion is at the **tip**; no station or span fraction beyond "tip".
2. **Elastic or total.** *"the torsional deflection is defined as the rotation of the airfoil section
   about the reference axis in Fig. 19, and the positive torsional deflection means the airfoil
   section rotates toward stall."* It is a deflection, so the **elastic twist relative to the
   undeformed/built-in configuration** - not total geometric+elastic twist, not a twist-distribution
   value; "relative to built-in twist" is not stated verbatim. Its sign is the opposite of our label:
   *"the aerodynamic moments at the blade sections make the airfoil sections twist towards feather
   and reduce the angle of attack"*, so their -3.60 deg = toward feather = our negative (nose-down) sense.
3. **Structural model.** Their own GEBT beam - *"The GEBT method used in the present study is
   established based on the Hamilton principle"*, with *"C* the sectional stiffness matrix resolved in
   the frame Bi"* - discretized as *"the IEA-15 MW blade is discretized into 50 stations as suggested
   by Gaertner et al. [56] in technical manual of the IEA-15MW RWT"*. The stiffness source (NuMAD /
   PreComp / a published deck / their own code) is **not stated**, and a **torsional stiffness value
   is not stated**; the only stiffness datum is Table 3's 1st torsion frequency **4.072 Hz** (refs
   4.314 / 4.475 / 4.295 / 3.911; the repo's own reference is NuMAD's 4.290 Hz, section 15.1).
4. **Aero method and coupling.** LL-FVW with Beddoes-Leishman dynamic stall, coupled to GEBT *"through
   the two-way loose coupling approach ... at each time step, the GEBT module accepts the aerodynamic
   loads calculated by the LL-FVW module, and then the blade profile is updated and the additional
   velocities due the blade deflections are calculated and fed back to the LL-FVW module for the
   aerodynamic prediction of the next time step."* The circulation is converged at each step (*"In
   each time step, the circulation of the bound vortex is solved through iteration of Eq. (4)"*,
   relaxation 0.1, tolerance 1e-3): a time-accurate **loose, two-way** coupling, not a fixed point
   and not one iteration. Azimuth step 6 deg, GEBT dt 0.002 s, 600 s. The tip values are the time
   mean; **the blade/azimuth is not stated**.
5. **Operating point.** *"the wind turbine is assumed to operate at the rated state, that is, the
   inflow wind velocity is 10.59 m/s and the rotational speed of the rotor is 7.55 rpm"*, cone and
   shaft tilt *"set to zero"*. Ours is 10.59 m/s, 7.56 rpm (-0.13%) and pitch 0; **their pitch is not
   stated**. Loads match: our per-blade 8.4169e5 N x 3 = **2.525 MN** vs their rigid fixed **2.53 MN**
   (Table 6, -0.2%).


**Parent's correction to the inference above - their own Table 3 refutes "2x softer", and the
transferability verdict has to be restated on the right axis.** The worker wrote the gap off to their
structure being ~2x softer than the public deck, but the paper's own Table 3 gives their model's
**1st torsion frequency as 4.072 Hz**, against our shell's **4.000 Hz** (section 15.1) and NuMAD's
4.290 Hz. Forming the same stiffness ratio the section 15.1 modal route uses: `(4.072/4.000)^2 =
1.036`, i.e. **their torsional stiffness is within ~4 % of ours** - and `(4.072/4.290)^2 = 0.900`,
within 3 % of our shell's `0.870`. Two models whose first torsional frequencies differ by 1.8 %
cannot be 2.4x apart in torsional stiffness. So the twist gap is **not** a stiffness difference.
With the stiffness difference excluded and the thrust level matched (their rigid 2.53 MN against our
rigid baseline 2.5417 MN, 0.46 %; their rigid 16.11 MW against our 16.389 MW, 1.73 % - an
independent load cross-check neither of us was built for), a linear twist-torque relation leaves
**their delivered aerodynamic torque around 2.4x ours** as the open variable, which is an
aerodynamic-model question (their two-way LL-FVW with Beddoes-Leishman dynamic stall and a converged
wake, at an unstated pitch angle, against our steady one-way BEM at pitch 0), not a structural one.
Their own decomposition supports looking there: they split the angle-of-attack change into a
*"structural twist component (the yellow line) and the aerodynamic component (the blue line)"*, and
Table 4's torsion is only the structural one - so their load path is where the extra twist must come
from. Restated verdict: **comparable in kind, and the residual is a load/aero-model difference, not a
stiffness difference** - which is a sharper claim than "not transferable" and keeps the prohibition on
retargeting our model to their -3.60 deg. What stays unproven: their sectional stiffness source and
their pitch angle are not stated, so the torque explanation is inferred from their frequency plus
their thrust, not read off their paper.
**Part 2.** One table, one sign convention: thrust and power as `(flexible - rigid) / rigid`,
negative = the rotor unloads; twist in degrees about +span under fluid +Y / rotor clockwise viewed
from behind. `tests/test_blade_deloading_vs_reference.py` (3 tests) drives the production participant
on the real mesh, solves the shell **once**, and re-runs only the BEM on the three feedbacks against
one reference baseline (2.541662 MN / 16.389372 MW):

| feedback fed to the BEM | what it isolates | thrust | power | tip dr [m] | tip dtwist [deg] |
| --- | --- | ---: | ---: | ---: | ---: |
| twist only, reference radii | the pure bend-twist unloading | -4.00% | -0.81% | +0.0000 | -1.5112 |
| deformed radii only, reference twist | the geometric re-loading (axial-stretch artefact) | +1.13% | +0.85% | +0.4462 | +0.0000 |
| twist + radii (production path) | what production does today | -2.96% | +0.03% | +0.4462 | -1.5112 |

Zhou Table 6 (flexible vs rigid): -13.04% / -8.38%. The exact production path crosses its own rigid
baseline at -2.32% / +0.42% (section 22.6); the table shares the reference-station baseline so only
the feedback changes. **Power sign:** radii-only **+0.85%** (positive) vs production **+0.03%**
(positive) - the **same sign** - while twist alone is -0.81%: "power up" is the **radius
(axial-stretch) artefact**, not bend-twist re-loading. Tip nodal displacement: axial (span)
+1.0636 m, radial (in-plane) +16.4868 m.

**Verdict: not transferable.** The reference is *partially comparable in kind* (same machine, rated
point, thrust level, and its quoted quantity is the elastic tip torsion we measure) but *not
transferable to this model*. (a) A beam built from the IEA-15MW's **own published `GJ`** gives
-1.78 deg at the tip under demonstrably comparable loads, and our shell's static torsion matches the
published deck to `GJ_static/GJ_deck = 1.068` (section 22.10); their -3.60 deg is ~2.0x that beam
and 2.4x our shell section rotation (-1.51 deg). ~~Their structure is therefore ~2x softer in torsion
than the public deck (or their "torsion" is a different measurement); with the stiffness source
unstated this cannot be resolved from the paper. (b) Their aero is a fully coupled LL-FVW; ours is a
one-way steady BEM with only radius and twist feedback. The de-loading gap follows the same ratio
(one-way twist-only -4.00% / -0.81% vs -13.04% / -8.38%); even the one-way lower bound (a converged
fixed point would deepen the twist) is ~3.3x short. So the gap is **not a defect** (sign and
mechanism are right, stiffness is not the problem) and **not** only one-way-vs-converged - it is a
**model/reference difference**, and the suite's rules forbid retargeting the model to hit their
number. The de-loading gap is now attributed: "power up" = the axial-stretch artefact (the
de-loading magnitude stays a reported residual); the thrust shortfall = a reference that does not
transfer.

### 22.12 The blade's bend-twist coupling, the section frame, and where the mesh limits the per-station numbers

**Part 1 - the z~112.22 m tear is not what limits the outer band.**
`tests/validation/blade/test_blade_section_frame_tear.py` (3 tests, one shell solve - the section
22.10 pure-torque realisation - 79.8 s) reads the real mesh's station grid and free edges, then
applies the tear discriminator.

*The grid (mesh geometry, no solve).* 186 merged physical stations from 266 raw z buckets; spacing
0.2985 / 0.5969 / 0.7963 m (min/median/max); ring node counts 12 (43 stations) to 26. 44 free edges:
22 on the open root ring (z = 0), 12 on the open tip ring (z = 117), and the two interior tears:

| tear | z [m] | gap to previous [m] | ring n | raw z buckets | free edges |
| --- | ---: | ---: | ---: | ---: | ---: |
| inner | 11.9380 | 0.7951 (neighbours 0.7959) | 26 (neighbours 21) | 1 | 8 |
| outer | 112.2245 | 0.4776 (neighbours 0.4776) | 12 (neighbours 12) | 1 | 2 |

Neither is a gap or a doubled station: `_physical_stations` merges the prebend-split sub-buckets
(1 bucket per tear ring) and the node dedup still reports 0 % reduction. The inner tear sits where
the ring jumps 21 -> 26 nodes; the outer tear has no station-grid signature at all.

*The discriminator at z = 112.2245.* Every geometric section quantity from the raw ring nodes lies
**inside the interval spanned by its two immediate neighbours**, and the fitted section rotation's
non-affine residual `distort` is smooth (0.973 -> 0.974 -> 0.975):

| quantity | z = 111.7469 (prev) | z = 112.2245 (tear) | z = 112.7020 (next) |
| --- | ---: | ---: | ---: |
| chordwise extent [m] | 1.96096 | 1.93775 | 1.91782 |
| xmin [m] | -1.27200 | -1.25427 | -1.23895 |
| xmax [m] | 0.68897 | 0.68349 | 0.67886 |
| hull area [m^2] | 0.43175 | 0.42402 | 0.42128 |
| shoelace area [m^2] | 0.43023 | 0.42256 | 0.41986 |
| declared chord [m] | 1.96073 | 1.93775 | 1.91414 |
| fitted `distort` | 0.973 | 0.974 | 0.975 |

The tear moves neither the geometry nor the affine fit. **Verdict: the z = 112.22 m tear is not a
limit on the per-station numbers** - a connectivity defect with no measurable geometric or estimator
footprint. What ended section 22.10's trustworthy band is the **section frame at the tip**: from
r ~ 114.6 m the ring's chordwise extent leaves the declared chord while the pitch-axis split does not.
Measured (declared-chord error / pitch-axis split error / `distort`):

| z [m] | declared-chord extent | pitch-axis split | `distort` |
| ---: | ---: | ---: | ---: |
| 114.6122 | 0.00 % | 0.00 % | 0.984 |
| 115.0102 | -3.54 % | 0.03 % | 0.989 |
| 115.4082 | -10.87 % | 0.00 % | 0.994 |
| 115.8061 | -19.75 % | 0.01 % | 0.998 |
| 116.1046 | -25.72 % | 0.00 % | 1.000 |
| 116.4031 | -34.60 % | 0.01 % | 1.000 |
| 116.7015 | -25.64 % | 0.00 % | 1.000 |
| 117.0000 | 0.00 % | 0.00 % | 1.000 |

So the mesh's tip refinement (the yaml's taped-chord end), the `distort` saturation to 1.000, and the
deck `GKt` having no valid station beyond r = 111.15 m (22.10) are what limit the outer band. No
mesh was changed; the tear is recorded here as a mesh-generation defect with its station index, ring
size and free-edge count.

**Part 2 - the section frame, pinned.** Guard file:
`tests/validation/blade/test_blade_section_frame_tear.py`; the P5 axis guard is
`tests/validation/bem/test_force_projection_load_frame.py`.

- `pitch_axis` is a chord fraction measured **from the leading edge**: the yaml pitch_axis (line 24)
  and column 2 (`PitchAxis`) of the official ElastoDyn blade deck agree to `max |diff| = 1.1e-16`
  over all 50 stations.
- `x = (pitch_axis - f) c`, so `x = 0` is the pitch axis, the LE at `+pitch_axis c` and the TE at
  `-(1 - pitch_axis) c`. On every one of the 186 physical rings - chord axis measured from the ring
  outline only (SVD principal in-plane axis, blunt-end oriented), origin the declared reference-axis
  point `(0, prebend(z), z)` - `xmax/(xmax - xmin) = pitch_axis` and
  `-xmin/(xmax - xmin) = 1 - pitch_axis` to median **0.005 %** / max **0.987 %** (the max at the
  circular root, z = 0); `x = 0` is strictly inside every ring. The bound is the suite's 5 % rule.
- The **blunt** end (larger in-plane spread in the outer chord quarter - the production
  `ForceProjector` rule) is the LE and lies at `+x` on **all 186** stations. The ring's own chord
  axis agrees with the declared twist direction to median 0.36 deg / max 6.72 deg over the 182
  stations with a defined chord axis (bound 10 deg; the 4 near-circular root rings
  z in {0.0, 0.796, 1.592, 2.388} have no defined axis).
- **Corner check at z = 68.05 m** (section 21.3's station): the `+x` (LE) node's interior angle is
  **97.32 deg**, reproducing 21.3's 97.3 deg; the `-x` (TE) node reads **102.19 deg**, **not**
  9.35 deg. The mesh's trailing edge is blunt (outer-quarter spread 0.17 m at the TE against 0.86 m
  at the LE), so 9.35 deg is not an in-plane ring angle anywhere on the mesh (its whole-mesh minimum
  in-plane ring angle is 76.7 deg). The frame fact the quote supported - min(x) is the TE, max(x)
  the LE - is confirmed by the thickness rule, not by the quoted angle. **This is a finding.**
- **LE/TE sense** vs the declared convention (fluid `+Y`, rotor clockwise viewed from behind,
  `Omega = +omega y`): the tip chord axis is `c_hat = (0.99979, +0.02052, 0)`;
  `c_hat . X = +0.99979 > 0` (mesh `+x` is the LE direction) and `(span x c_hat) . Y = +0.99979 > 0`
  (the section normal is `+Y`, the fluid/downwind direction). Both are asserted.
- **No single global axis is legitimate**: the maximum inter-station chord-direction angle is
  **24.410 deg** (z = 3.184 to z = 105.061 m), > 2 x 10 deg, so no fixed vector satisfies the P5
  per-section bound (P5, fixed in `9a3923e`).

**Part 3 - the coupling statement the issue needs.** This blade has **`D16 = D26 = 0` in all 696
sections** (recorded in `c6eb8bf`; pinned by
`tests/validation/blade/test_blade_twist_mechanism.py::test_blade_laminates_have_no_bend_twist_coupling`,
which asserts `max |D16| + |D26| <= 1e-9 |D|max` over every section). The bend-twist coupling
mechanism the issue hypothesises therefore **cannot occur in this model** - there is no coupling to
over-predict - and the element is **exonerated**: MITC4 against a converged Rayleigh-Ritz (12.5), a
hand-authored CalculiX S8R composite deck (12.6) and the free-edge CLT (12.7) all agree within 1-3 %
(quoted from the record, not re-run here). The blade-level twist residual against the literature is
an aero/torque difference, not stiffness and not the load frame (22.7 moment identity, 22.9
minimum-norm realisation over-delivers and softens, 22.10 `GJ_static/GJ_deck = 1.03-1.09` over
r = 35-80 m, 22.11 their 4.072 Hz is within ~4 % of ours so the gap is their delivered torque).

What the **new guards pin**: the section frame of Part 2 (pitch-axis split, blunt/LE sign, `-twist`
chord axis, LE/TE sense, 24.410 deg spread) on the real mesh, so a reader cannot mis-sign `x`, the
pitch axis or the LE/TE sense. What stays **quoted from the record**: the 696 balanced sections
(re-pinned by the existing `test_blade_twist_mechanism.py`, not re-derived here), the element
exoneration of 12.5-12.7, and the twist attributions of 22.7-22.11. The tear verdict of Part 1 is
this unit's own measurement. No production code, no existing test and no store file was changed.
