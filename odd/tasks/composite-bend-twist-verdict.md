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
