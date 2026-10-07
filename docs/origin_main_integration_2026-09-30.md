# Integrating `origin/main` into the IEA 15 MW validation line — measured result

**Date**: 2026-09-30
**Branch**: `integrate/origin-main-2026-09-30` (from `integrate/main-into-rust` @ `d571917`)
**Merge commit**: `f5f92f7` (merges `origin/main` @ `dffacab`, 169 commits after the previous
merge `e2bc081` @ `7a2196f`)
**Feature log**: `odd/tasks/integrate-origin-main.md`

> **Correction pointer (2026-10-05) — do not carry forward the coupled-FSI columns.**
> Commit `1146265` later fixed three linked sign errors in the blade's torsional load
> chain, arbitrated against the deck's own aerofoil geometry by
> `tools/diagnose_sign_chain.py`. The convention it settles: the deck puts the leading
> edge at **+x** and the load-frame downwind (thrust) direction at **+y**, so a rigid
> **+z** section rotation is **nose-down** and nose-down is **`omega > 0`**. Every number
> in this document was measured before that fix. In particular, a reader must NOT carry
> forward the **`flap mean`**, **`thrust`**, **`power`** and **`CP`** columns of the
> "Yaw sweep, after vs `_mitc3fix`" table, the `yaw loads (flap/thrust/power)` row of
> the "Verdict, item by item" table, or the absolute `flap`/`thrust`/`power` levels of
> the `h`/`dt` convergence tables: those are coupled-FSI campaigns whose elastic twist
> and load feedback ride the pre-fix sign chain. The corrected production-path
> de-loading (`-25.31% / -14.77%` thrust/power, against Zhou Table 6's
> `-13.04% / -8.38%`) supersedes the yaw-sweep de-loading reading. The static blade
> matrix, the CCX convergence tables, the element deltas and the rigid-BEM results do
> not depend on that sign chain and stand as recorded. The already-recorded config
> confound (mesh + hub radius) is a second, separate reason the coupled columns are not
> an element A/B.

## Why

`origin/main` carries the reviewed MITC4+/D element, the span-relative ply-angle fix
(`9de3731`), the "uncorrected section shear" material fix (`8cbfc0b`), the CCX writer
id-scheme fix (`930d055`), the modal-solver fixes (`b40561d`, `7e395ef`) and the OpenFAST
AeroDyn comparison (`d10d26d`). We wanted those on the branch where the IEA 15 MW
validation suite runs, and a measured before/after.

## Method

Same command, same environment, same node, before and after:

```bash
module load gcc/14.2.0_sequana
export LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64:/scratch/app_sequana/gcc/14.2.0/lib:$LD_LIBRARY_PATH
python -m pytest tests/ -q --tb=short \
    --ignore=tests/test_blade_mesh.py --ignore=tests/test_rotor_inertial.py \
    --ignore=tests/test_vol_mesh.py
```

Rust rebuilt with `maturin develop --release` at each measurement point.
Logs: `$SCRATCH/tmp/odd-integrate-origin-main/{before-full,after-full3}.log`.

| | collected | passed | failed | errors | skipped | xfailed |
|---|---|---|---|---|---|---|
| **before** (`cce8165`) | 974 | 933 | 6 | 0 | 26 | 9 |
| after, first merge pass (`28339a4`) | 987 | 905 | 15 | 22 | 38 | 6 |
| after, all of `main`, CCX unblocked, ignores dropped (`0cfa523`) | 974 | 901 | 18 | 0 | 38 | 17 |
| after the K_G unit fix and the `airfoil_spacing` pins (`4ff9a6c`) | 974 | 906 | 13 | 0 | 38 | 17 |
| after the delta investigation (box, D-Tube, elastica) | 974 | 911 | 7 | 0 | 38 | 17 |
| **after S-7, V-02, the composite FRD fix and the blade span direction** (`69ec5d2`) | 974 | **917** | **4** | **0** | 38 | 14 |

The four that remain are named, not hidden:

| test | why it is red |
|---|---|
| `test_corotational_is_frame_objective_tl_is_not` | deliberate: the MITC4 path of `assemble_kt_corotational` uses upstream's total-Lagrangian tangent, so the property the test asserts no longer differs |
| `test_static_nonlinear_ul_elastica` | converged ~4% stiffer than the exact elastica (the load-stepping table is in the docstring); an accuracy question about the UL path |
| `test_dtube_camarena_bending` | converged 1.07% outside the 1% bound, ~10x the transverse-shear estimate (table in the module docstring) |
| `test_iea15mw_s2[LC1_gravity-utd]` | the UTD blade input, whose layup divergence is already documented; the `official` variant passes |

The last row is the honest one to compare against `before` only with a caveat: the two
later merges (`8f3deca`, `0cfa523`) brought upstream's rewritten test files, which
replaced several weak or tautological tests with fewer, stronger ones -- upstream's
`test_shell_validation_fixed.py` and `test_shell_comprehensive.py` load the strip to
~0.17 L now and compare against the closed-form Bisshopp-Drucker elastica instead of
asserting a loose O(L) window.  Seven fewer passes with the same eighteen failures is
that change of instrument, not a regression.  The two `--ignore`s for
`test_blade_mesh.py` and `test_rotor_inertial.py` are gone: upstream fixed both.

The suite grew by 13 tests; `origin/main` deleted three solid/volumetric test files and
added its own CalculiX and AeroDyn parity suites (~150 tests, of which some skip).

## What the merge fixed

1. **The CalculiX writer's id scheme.** Before, `write_ccx_mesh` labelled the `.msh` nodes
   by 1-based index but wrote connectivity and `.nam` sets from raw entity ids, so with
   drifted ids the deck was corrupt (CalculiX reported 696 743 nodes for a 3 333-element
   mesh) and `_build_angle_bucket_sets` raised `IndexError: index 7500 is out of bounds for
   axis 0 with size 125`. After `930d055` (plus the id→index mapping applied to
   `_build_angle_bucket_sets`, which upstream had missed) the deck is consistent and
   `test_composite_ccx_parity` passes its linear cases.

2. **The composite gravity load (found while measuring).** `assemble_f_body` handed
   `material_rho` to the shell body-load kernel, which multiplies by the element thickness,
   but for a composite that function returns mass per *area*. The load came out a thickness
   too small: on the IEA 15 MW blade the root reaction was 51 440 N against an assembled
   weight of 693 771 N (92.6% off) and the gravity deflection collapsed to 5% of the beam
   value. Fixed with `body_load_rho`; `test_iea15mw_v03_static_gravity` and the S-2
   `LC1_gravity` cases pass again (12 passed).

3. **Three of our interfaces were dropped by the module split and are now restored**:
   `PyMeshAssembler::assemble_kt_corotational`, plus the core `centrifugal_load`,
   `assemble_geometric_k_from_disp` and `update_node_coordinates` that
   `src/aeroelast/core/assembler.py` and the S-0/S-4 tests call.

4. **Four line-level union artifacts repaired** (merge damage, not physics):
   `force_projection.py` lost `hub_r = hub_radius if hub_radius is not None else
   blade_aero.hub_radius` (10 failures), `test_composite_beam_parity::_run_ccx` lost its
   `return` (3 failures), `test_mitc3_benchmarks.py` mixed ours' load helper with theirs'
   sign assertion (5 failures), and `generators.py` lost the `ElementSet` import
   (1 failure). All four are now green.

## What the merge did not fix, and what it moved

### A. One blade-mesh defect blocked every CalculiX parity test (28 tests) — **RESOLVED 2026-09-30**

**Root cause: our own doing, not the element's.**  Commit `e705820` (2026-05-12,
"default airfoil_spacing to constant to reduce tip element AR") left two workarounds
in `BladeMesh`, both aimed at the aspect ratio of the *previous* element:

- `_refine_high_gradient_sections` skipped a refinement when
  `gap / 2 < element_size * min_ar_ratio` (4.0), leaving the tip coarser;
- `airfoil_spacing` defaulted to `"constant"` instead of the library's `"cosine"`.

At `element_size = 0.5 m` that produced **9251 nodes where upstream's own test pins
9277**, and the missing tip refinement left a sliver element CalculiX cannot
integrate.  Removing both workarounds restores the reference mesh exactly (9277
nodes) and CalculiX finishes (`Job finished`, no `e_c3d`).  Fixed in `2561aad`.

Recovered: `test_blade_mesh.py` (2 passed), `test_ccx_shell_element_types_parity.py`
(4 passed, after its element-block assertion was relaxed from a pinned
`ELSET=Eall` to the emitted TYPE set — this writer names blocks per region so the
per-set `SHELL SECTION` cards can reference them), upstream's
`test_blade_iea15mw_validation.py` (10 passed / 8 xfailed, as upstream reports),
and our five `test_blade_ccx_parity.py` cases now *execute* instead of dying at the
deck.

Those five are now a measurement, not a blocker: they report 78-89% displacement
gaps and a +70 deg CCX twist against -4.9 deg for the shell.  That is far outside
an element tolerance, so the deck the test builds is the suspect, not the element.
Unresolved.

### A-bis. The structural anchors did not move with the mesh fix

Re-measured after `2561aad`: the eight structural failures are identical, so they
are element behaviour and not the tip mesh.

### B. Our validation anchors moved (element-behaviour deltas, need re-baselining)

`origin/main` is fully merged as of `8f3deca` (`5bfa2b2`).  No band was touched:
per the user's decision each delta gets investigated against an independent
reference first.

| Anchor | Before | After | Reading |
|---|---|---|---|
| V-02 2nd flap | 1.6946 Hz (+2.15%) | 1.5545 Hz (−6.3%, tol 5%) | new element ~8% softer in 2F |
| S-4 rotating 1F vs OpenFAST MBC3 | 0.5698 Hz (+0.6%) | 2.1375 Hz (+277%) | mode identification or K_G scale — see below |
| S-7 shell/beam torsion ratio | twist 0.926 (`GJ` 1.080) | twist 1.406 (band ≤1.3) | **superseded 2026-10-07**: the two columns compared a `GJ` ratio against a twist ratio, which are reciprocals, and 1.406 does not reproduce at HEAD at any mesh size (0.9396 measured, 0.9336 converged, i.e. *stiffer*). See `odd/tasks/s7-torsion-ratio.md`, issue #20, and the § S-7 closure |
| S-6 one-way dynamic OoP mean/std | in band | `nan` | a NaN enters the S-6 path |
| Box EI vs analytic | <2% | 3.27% | |
| D-Tube tip vs beam | +0.2% | 4.75% | |
| UL elastica α=1 | −0.5% | 3.63% | |
| Composite bend-twist vs CCX | 0.4% (2026-09-25) | AE −0.4499 vs CCX −0.1584 (2.84x) | the `_build_angle_bucket_sets` crash had been hiding this; it is the 3.23x question again |

Two candidate causes for the S-4 jump, both unresolved:

- ours' `assemble_geometric_k` filtered the membrane stress to its **tensile part**
  (`tensile_part_membrane`) before assembling K_σ; upstream's does not. The merge dropped
  the filter, so K_G is now assembled from the full (partly compressive) membrane state.
- upstream's modal fixes (`b40561d`, `7e395ef`) changed which eigenpairs the solver returns
  and in what order. If the S-4 test classifies "1st flap" by index, a different mode is
  being labelled: 2.1375 Hz is plausible as a higher flap or edge mode.

### C-bis. The K_G unit convention: MITC3 and MITC4 disagreed (found 2026-09-30)

`assemble_geometric_k(sigma)` is one API fed by one `sigma` array, so both shell
families must read the same quantity.  MITC3's `compute_k_sigma_local` multiplies the
membrane stress by the element thickness to obtain the force resultant N [N/m].
The reviewed MITC4 did not: `geometric_stiffness_contribution` used `sigma` directly as
if it were already N/m.  On a mixed mesh the quad part of K_G therefore came out
**1/h too stiff** while the triangle part was right, and the IEA 15 MW blade is mostly
quads.

Measured per term at 7.56 rpm (`$SCRATCH/tmp/odd-integrate-origin-main/diag_s4.py`):

```
K               1F = 0.5336 Hz   (parked; the pre-stress chain is correct:
K + K_G         1F = 2.1298 Hz   |u|inf = 0.3759 m, exactly the documented
K + K_SP        1F = 0.5328 Hz   pre-bend straightening; radial load ~1.2 MN)
K + K_G + K_SP  1F = 2.1293 Hz   ref 0.5666 Hz
```

With the factor of `h` restored (`4ff9a6c`):

```
K + K_G         1F = 0.5565 Hz
K + K_G + K_SP  1F = 0.5557 Hz   (-1.9% vs 0.5666)
```

and the centrifugal stiffening contributes **+4.3%** (0.5336 -> 0.5565) against the
**+4.2%** OpenFAST measures (0.544 -> 0.5666) -- the physical cross-check that the fix
is the right one, not a tolerance fit.

This single fix closed **S-4** (3/3) and **S-6** (3/3, the NaN was the same inflated
K_G in the one-way dynamic path).  `tensile_part_membrane`, the other candidate, is a
no-op for uniaxial tension (a rotating blade) and was not the cause.

**Why upstream never caught it:** their K_G unit tests assert zero-stress-is-zero,
finite entries, symmetry, and "local equals the local transform" -- none pins the
magnitude or the units, so the mixed convention passes all four.  This is worth
reporting upstream.

### D-bis. Our own blade CCX test: a metric bug, and one real finding (2026-09-30)

`tests/test_blade_ccx_parity.py` is ours, and it was failing on all five cases with
78-89% gaps.  Two defects of its own, plus one real question:

1. **The metric compared two different node sets.**  `_tip_metrics` derived its
   index list as `sorted(disp_nodes)`, and while the CalculiX side passes only the
   tip nodes, the AeroElast side passed a dict covering every node.  AeroElast's
   "tip" displacement was therefore the mean over the whole blade (1.796 m) against
   CalculiX's real tip (8.469 m).  With the tip set named on both sides the flap
   case is **6.3% apart** and the two axial cases 17.0% -- inside the 20% band.
2. **The twist was noise.**  A least-squares rotation over the thin tip ring has a
   tiny denominator; the CalculiX side of the flap case returns +70 deg, which no
   1e5 N flap load produces.  The twist is now reported, not asserted.
3. **`tip_torsion` is not comparable.**  The CalculiX writer has no validated
   couple path: sending the moment as a `*CLOAD` on the sixth DOF makes CalculiX
   return 418 m and a -42826 deg "twist".  Removed from the case list with the
   reason recorded; it belongs with S-7/S-8, which use a force couple.

**Still red, on purpose: `tip_edge`.**  AeroElast gives 1.229 m where CalculiX
gives 3.450 m (64%), while flap matches to 6%.  The flap/edge stiffness ratio is
6.5 for AeroElast and 2.5 for CalculiX, so the disagreement is in the section's
**edgewise** stiffness -- a question about the laminate mapping, and the first
item on the delta list below.

The stronger evidence sits in upstream's own suite: on this tree
`test_blade_iea15mw_validation.py` passes its flapwise static tip against CCX and
all five modal frequencies against CCX S8R (11 passed, 7 xfailed by declared
limits).  The blade pipeline and the element are verified there.

### D-ter. Blade vs CalculiX: the two formulations do converge to the same value

MITC4 (linear) and S8R (quadratic) are different formulations, so a fixed-mesh gap
says nothing on its own.  `tools/blade_ccx_convergence.py` refines the mesh and
feeds both sides the same one:

| element_size | nodes | flap | edge | axial |
|---|---|---|---|---|
| 2.0 | 1 460 | 10.56% | 11.60% | 12.43% |
| 1.0 | 3 043 | 6.22% | 5.15% | 7.43% |
| 0.5 | 9 277 | **3.68%** | **1.56%** | **4.90%** |

The gap shrinks monotonically for every load case, so the residual is
discretisation and not a model difference.  Recorded caveat: the metric is the mean
over the tip nodes and moves with the tip node distribution (flap reads 7.55 / 7.94
/ 7.12 m across the three meshes), so the trend of the gap is the signal; a
mesh-independent version compares the work done by the load.

### D-quater. `tip_edge`: our test omitted the span direction (resolved)

`tip_edge` was the last red case at 64%, while flap matched to 6%.  The cause was
this test's AeroElast side: it built the assembler with the direct
`PyMeshAssembler(...)` constructor, which has **no `span_direction` argument**, so
the element ply angles were read as element-local instead of span-relative.  Flap
is insensitive to that; edge is not:

```
flap  mine (no span dir) 7.934 m | from_model(span dir) 7.942 m | CCX 8.469 m
edge  mine (no span dir) 1.229 m | from_model(span dir) 3.272 m | CCX 3.450 m
```

The edgewise stiffness came out **2.7x too high** without it, which confirms that
`9de3731` ("supply the span direction so ply angles are span-relative") is
load-bearing and names the failure mode when it is omitted.  The test now runs both
sides through upstream's pipeline (`Blade` + `get_element_properties()` +
`from_model(..., SPAN_DIRECTION, ...)`); all four cases pass: flap 6.22%, edge 5.15%,
axial tension and compression 7.43%.

## The blade structural matrix, before and after (0.25 m mesh)

`tools/run_blade_structural_matrix.py --element-size 0.25`, run on `cce8165` (before)
and on `063af4d`-equivalent HEAD (after).  The before run reproduces the recorded
baseline exactly -- 32 321 nodes, 33 454 elements, 684 element sets, B1 1.42%,
B2 9.12%, B3 3.85%, B4/B5 7.49%, B7 1.70% -- which validates the tool and the
comparison.

| case | before shell [m] | after shell [m] | change | before vs beam | after vs beam |
|---|---|---|---|---|---|
| B1 tip flap | 67.441 | 67.395 | −0.07% | 1.42% | 1.61% |
| B2 tip edge | 35.369 | 35.161 | −0.59% | 9.12% | 8.37% |
| B3 distributed | 7.271 | 7.256 | −0.21% | 3.85% | 4.43% |
| B4 traction axial | 0.455169 | 0.455945 | +0.17% | 7.49% | 7.66% |
| B5 compression axial | −0.455169 | −0.455945 | +0.17% | 7.49% | 7.66% |
| B7 gravity | −1.935772 | −1.929050 | +0.35% | 1.70% | 3.11% |

The mesh moved with the `e705820` fix (32 321 → 32 336 nodes, 684 → 696 element
sets), and the beam column is sampled at the shell's span stations, so it moves
too -- B7's beam value goes from −1.9693 to −1.9910 m (+1.1%), which is most of
that row's error change.

**Verdict: on the blade's static matrix the reviewed element changes the shell's
own response by at most 0.6%.**  That is a very different picture from the unit
cases (box, D-Tube, UL elastica, 3-5% each) and it is the useful headline: the
blade-level static behaviour is essentially unchanged, while the anchors that
moved are the small-geometry unit cases.

## G4 FSI campaigns: the finished before/after (yaw sweep and h/dt)

The first two G4 campaigns finished, so the recorded snapshot can be diffed:

```bash
python tools/campaign_metrics.py collect \
    $SCRATCH/frontiersin_results_corotational_100s \
    $SCRATCH/convergence_b1_b2_rerun --csv /tmp/after_yaw_hdt.csv
python tools/campaign_metrics.py compare /tmp/after_yaw_hdt.csv \
    --alias frontiersin_results_corotational_100s=frontiersin_results_corotational_100s_mitc3fix \
    --alias convergence_b1_b2_rerun=convergence_b1_b2_results_official
```

A re-run directory does not carry the recorded campaign's name, so `compare` now
takes an explicit `--alias NEW=OLD` (and refuses to guess when two recorded states
exist).  The snapshot holds **two** recorded yaw states: `_mitc3fix` (the immediate
predecessor, flap 8.02 m) and `_twistfix` (flap 8.23 m).  The primary column below
is `_mitc3fix`; `_twistfix` is quoted once as the secondary.

### The confound: before and after are not the same configuration

**The G4 delta cannot be attributed to the element.**  The after campaigns started
2026-09-30 22:55, after `8e6488e` (20:30) fixed two pre-flight defects on 30 solid
case YAMLs:

- `airfoil_spacing: constant -> cosine` (the old `e705820` workaround), so the yaw
  case now builds the reference 32 336-node / 33 473-element mesh;
- `force_projection.py` no longer honours an explicit `hub_radius: 0.0` when the
  aero has a hub, removing a 3.97 m shift of every strip.

The before campaigns (`_mitc3fix`, 2026-09-23 and earlier) ran the old config.

`Thrust/CT` is the effective dynamic pressure `q·A`; it dropped a **uniform 4.76%
across all five yaw cases** (effective inflow -2.4%).  A yaw-independent,
global shift is what a hub-radius or mesh change produces and what an element
change cannot.  The element's own static effect is already measured at ≤0.6%, and
the campaign flap moved ~8%, thirteen times that.

### Yaw sweep, after vs `_mitc3fix`

| yaw | flap mean [m] | flap Δ | thrust [MN] | thrust Δ | power [MW] | power Δ | CP |
|---|---|---|---|---|---|---|---|
| 0  | 8.022 → 7.341 | -8.48% | 1.726 → 1.639 | -5.02% | 12.964 → 12.342 | -4.80% | 0.3681 → 0.3680 |
| 10 | 7.957 → 7.285 | -8.44% | 1.699 → 1.614 | -5.02% | 12.485 → 11.889 | -4.77% | 0.3545 → 0.3545 |
| 20 | 7.757 → 7.112 | -8.32% | 1.620 → 1.539 | -5.01% | 11.111 → 10.589 | -4.69% | 0.3155 → 0.3158 |
| 30 | 7.425 → 6.825 | -8.09% | 1.490 → 1.416 | -4.97% | 9.008 → 8.604 | -4.48% | 0.2559 → 0.2566 |
| 40 | 6.957 → 6.422 | -7.70% | 1.312 → 1.249 | -4.86% | 6.407 → 6.156 | -3.92% | 0.1820 → 0.1836 |

Edge mean moves -7.83% to -6.26%.  Against `_twistfix` the same re-run reads flap
-10.8% to -9.6% and power -6.2% to -4.7% (the old `_twistfix → _mitc3fix` step plus
this one).  The coefficients are flat because they are normalized by the same
effective inflow that shifted; the dimensional loads are what moved.  Versus the
recorded literature bands (power -12% to -34%, flap -42%) the ~5% drop deepens the
deficit -- an aero-loading consequence of the corrected aero geometry, not a
measured element effect.

### Convergence, after (read within the re-run)

The `convergence_b1_b2_results_official` column is flap **17.84 m** -- a different
geometry, so it does not read as a convergence comparison.  The re-run does:

| pair | flap mean [m] | gap | order | thrust [MN] | power [MW] |
|---|---|---|---|---|---|
| h 0.5 → 0.25 m | 7.7344 → 7.3805 | 4.58% | -- | 1.659 → 1.641 | 12.500 → 12.355 |
| h 0.25 → 0.125 m | 7.3805 → 7.3187 | **0.84%** | ~2.4 | 1.641 → 1.639 | 12.355 → 12.332 |
| dt 0.02 → 0.005 s (0.25 m) | 7.3806 → 7.3804 | **0.003%** | -- | 1.6414 vs 1.6414 | 12.3546 vs 12.3545 |

**Temporal convergence is closed** at the output resolution.  **Spatial converges
second-order out to 0.25 m**: the observed order is ~2.4, and a Richardson fit puts
the 0.125 m flap within ~0.2% of the limit (~7.307 m), so 0.125 m is the usable
reference mesh and 0.5 m is the outlier.  (The earlier "h not converged" reading,
7.734 vs 7.053, was a truncated `h_fine`; the complete run settles at 7.319.)

### Verdict, item by item

| item | before -> after | read |
|---|---|---|
| element, static (blade matrix) | ≤0.6% | **unchanged**; the FSI deltas cannot come from here |
| yaw loads (flap/thrust/power) | -8.5/-5.0/-4.8% at yaw 0 | **not an element verdict**: the campaigns moved onto the reference mesh and correct hub radius; the uniform `q·A` -4.76% is that config fix |
| yaw coefficients (CP/CT) | ~flat (<=+0.9%) | consistent with an inflow-normalized shift, not a new aero efficiency |
| yaw vs literature | deficit ~5 pp deeper | worse vs the band, but a consequence of the corrected geometry; needs a physical reassessment, not a re-baseline |
| temporal convergence (dt) | 0.003% | **closed** |
| spatial convergence (h) | 0.84% at 0.25→0.125, order ~2.4 | **healthy**; 0.125 m is the reference |

A clean element A/B still needs one campaign re-run with the *other* config held
fixed; the recorded before/after is a **config** before/after.

## What is left, and its cost

| cost | item |
|---|---|
| done | blade matrix before/after (above), the counter reset, the merge of `034ba81` |
| done | G4 yaw sweep and B1/B2 h-dt (verdict above; config-confounded, temporal closed) |
| ~30 min | the work-based (global scalar) version of the blade convergence check |
| cancelled | every remaining G4 job: the V-06 (`11604992/93`, queued, normal-operation and biased by the twist over-deloop, issue #9), ch6 (`11605077`) and parked (`11605081`) -- analysed, then stopped; the parked run diverges (issue #10) |
| report upstream | the MITC4 `K_G` unit mismatch; the two over-specified assertions in `test_blade_iea15mw_mesh_convergence.py`; the composite `0.819` known-behaviour back-out |

### D-quinquies. The work-based convergence scalar, and the issue filed

**The MITC4 `K_G` unit mismatch is reported upstream**: `efirvida/AeroElast#7`, with
a reproducer that uses only this repository's public API -- one element, one stress,
two thicknesses -- so it does not depend on the FSI path this branch adds.  The
reproducer against `2fd847d`: MITC4's `|K_sigma|` is the same for h = 0.10 and
h = 0.01 (ratio 1.000 where the theory wants 10), while MITC3 scales correctly
(10.000).  The FSI numbers (S-4's +276%, the S-6 NaN, the +4.3%-against-+4.2%
cross-check) are carried in the issue as context, labelled as coming from a
downstream branch.

The node-mean tip displacement only had a readable *trend* (flap reads 7.55 / 7.94
/ 7.12 m across the meshes), so `tools/blade_ccx_convergence.py` now also reports
the work done by the load, `W = 1/2 sum f_i . u_i`, a global scalar whose two load
representations converge to the same uniform traction:

| element_size | nodes | flap | edge | axial |
|---|---|---|---|---|
| 2.0 | 1 460 | 10.70% | 11.65% | 14.55% |
| 1.0 | 3 043 | 6.38% | 5.17% | 9.52% |
| 0.5 | 9 277 | **3.87%** | **1.55%** | **7.73%** |

Every column falls with the mesh, so the convergence claim now rests on the global
scalar rather than on a node-count-dependent average.  Axial settles slowest, as
expected: the resultant is spread over the tip section, so its convergence carries
the section's local distortion too.

### C. A deliberate semantic change

`assemble_kt_corotational`'s MITC4 path now uses upstream's total-Lagrangian tangent
(`mitc4::compute_kt_global`) because the reviewed element does not expose a corotational
variant. That is the honest adaptation, and it invalidates
`test_corotational_large_rotation_validation::test_corotational_is_frame_objective_tl_is_not`,
which asserts the very property that no longer differs. The test needs a decision: port our
`Mitc4Precomputed::compute_kt_corotational` into the new element, or restate the test
around the TL tangent.

## Inherited limitations (not regressions)

- 3D solid elements are gone, in Rust (`7b295f3`) and Python (`9230ea2`). The upstream
  deletions took `tests/test_solid_elements.py`, `tests/test_vol_mesh.py` and
  `tests/test_beam_4cases_parity.py` with them; the surviving CCX tests read
  `tests/_ccx_io.py` instead.
- The volume-mesh STL boundary-face path in `write_meshio` is gone with them; the STL
  writer now caps every open boundary loop (upstream capped only the tip loop, which the
  rotor export test rejects).
- The `--ignore=tests/test_vol_mesh.py` in the documented command is now a no-op.

## Next steps (decisions needed)

1. ~~**Blade mesh element 2790**~~ — resolved in `2561aad`; see section A above.
   The remaining CCX item is our own `test_blade_ccx_parity.py` deck (78-89% gaps,
   +70 deg twist), which has to be checked against a deck built by the S-4/S-5
   tooling before its numbers mean anything.
2. **S-4**: restore or justify the tensile-part filter, and check the mode classification
   against the new modal filtering. This is the one anchor whose failure could indicate a
   real regression rather than a re-baselining.
3. **S-6 NaN**: trace where the one-way dynamic path produces NaN.
4. **Investigate each delta against an independent reference** (user decision, no band
   edits): V-02 2F 1.5539 vs 1.6590 Hz, S-7 1.406 (superseded 2026-10-07:
   HEAD measures 0.9396 and converges to 0.9336), box EI 4.8194e8 vs 4.6667e8,
   D-Tube 23.359 vs 24.525 m, UL elastica 0.29075 vs 0.30172.
5. **K_T corotational**: port `Mitc4Precomputed::compute_kt_corotational` or restate the
   frame-objectivity test.
6. **G4 FSI campaigns**: the yaw sweep and the B1/B2 h-dt re-runs are finished and
   verdicted (section above).  The rest were cancelled after analysis: the
   normal-operation cases carry the twist over-deloop, and the parked run diverges
   (issues #9 and #10).  The one open methodological item is the config confound: a
   clean element A/B needs a re-run with the case config held fixed.
7. **`origin/main` moved again** (checked 2026-10-01): `5f189fa -> 31565b3`, 13
   commits, and a new `validation-2026-09` tag.  Contents: nine `numad` refactor lint/
   typing commits (the xlsx/legacy blade-import path), two **unwired** MITC3 fixes
   (`089cd54` `union_rotation` rotational block; `9d394de` the strain-smoothed MITC3+
   relative frame, closing upstream `#2`, both in the smoothed/union path the MITC4
   blade does not use), one CI workflow plus an env pin (`bdbc6c0`), and one docs
   reconciliation (`2446335`).  `git merge-tree --write-tree HEAD origin/main` exits 0
   with no conflicts.  It is low-risk for the recorded campaigns, but it needs a Rust
   rebuild and a suite re-run -- so **do not merge/rebuild while the G4 jobs are
   queued or in flight**, or the G4 set ends up on two different binaries.
