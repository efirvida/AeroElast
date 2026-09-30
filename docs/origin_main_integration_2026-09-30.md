# Integrating `origin/main` into the IEA 15 MW validation line — measured result

**Date**: 2026-09-30
**Branch**: `integrate/origin-main-2026-09-30` (from `integrate/main-into-rust` @ `d571917`)
**Merge commit**: `f5f92f7` (merges `origin/main` @ `dffacab`, 169 commits after the previous
merge `e2bc081` @ `7a2196f`)
**Feature log**: `odd/tasks/integrate-origin-main.md`

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
| **after** (`28339a4`) | 987 | 905 | 15 | 22 | 38 | 6 |

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
| S-7 shell/beam torsion ratio | 1.080 | 1.406 (band ≤1.3) | softer in torsion, the expected direction of `8cbfc0b` |
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
   edits): V-02 2F 1.5539 vs 1.6590 Hz, S-7 1.406, box EI 4.8194e8 vs 4.6667e8,
   D-Tube 23.359 vs 24.525 m, UL elastica 0.29075 vs 0.30172.
5. **K_T corotational**: port `Mitc4Precomputed::compute_kt_corotational` or restate the
   frame-objectivity test.
6. **G4 FSI campaigns**: still to be relaunched on the cluster, cheapest first.
