# ODD — force projection: realise a section moment as a wall shear flow (#11)

Base: `integrate/origin-main-2026-09-30` after merging `origin/main` (`26ffe6e`).
Issue: https://github.com/efirvida/AeroElast/issues/11 (bug).

## State update (2026-10-05) — the sign chain is fixed

Commit `1146265` fixed the three linked sign errors this plan tracked and arbitrated the
convention against the deck's own aerofoil geometry (`tools/diagnose_sign_chain.py`: every
real station ring matched to its WindIO aerofoil, residual 0.5-3.6% of chord). The
convention, stated once: the deck puts the leading edge at **+x** and the load-frame
downwind (thrust) direction at **+y**, so a rigid **+z** rotation moves the leading edge
downwind and **nose-down is `omega > 0`**. Consequences: `ForceProjector` records
`_strip_moment_axis_sign` and applies `Mp` on it; `fsi_participant._twist_mesh_to_bem` is
**+1**; the multi-cell realisation (skin `q_i`, shared webs `q_i - q_j`, total moment
exact, net force zero) is validated against an independent hand-assembled Bredt-Batho
system, and the two realisations of the same rated moment agree in sign.

Corrected rated numbers (same point throughout):

| quantity | pre-fix | corrected |
| --- | ---: | ---: |
| tip section rotation, minimum-norm | -0.8696 deg | **+8.1048 deg** |
| tip section rotation, with-properties multi-cell | +0.0773 deg | **+9.6669 deg** |
| de-loading, twist only (thrust/power) | -4.00% / -0.81% | **-26.01% / -15.61%** |
| de-loading, radii only | +1.13% / +0.85% | +1.24% / +0.95% |
| de-loading, production path | -2.96% / +0.03% | **-25.31% / -14.77%** |
| Zhou Table 6 reference | -13.04% / -8.38% | (unchanged) |

The pre-fix numbers that follow in the merge-reconciliation and T2c sections are that
record, not the current state. T2c-2b and T3 are closed by `1146265`; the tube's Bredt
validation is untouched at `1.0031x` (0.3091%), and the sign is now pinned by
`tests/test_multicell_shear_flow.py::test_min_norm_and_with_properties_realisations_of_the_rated_moment_agree_in_sign`.

## Problem

`ForceProjector._distribute` (`src/aeroelast/solvers/bem/force_projection.py`) realises a
strip's section moment as a constrained **minimum-norm** nodal force system. On the one
geometry where the exact answer is independent of the code under test (the closed thin-wall
tube), that field is **circle-tangential about the section centroid**, `f_j = omega x d_j`
(residual 2.3e-13), not the constant wall shear flow `q = T/(2A)` a closed thin-walled section
carries. Measured on the validated tube, same torque and solver: shear-flow realisation =
1.0031x Bredt (passes the 5% rule), production `_distribute` = 31.6933x Bredt at the ends.

Two consequences the issue pins:
- it is an **end / Saint-Venant** effect (clamped-tip ratio decays 8.86 -> 3.20 -> 2.10 with
  L = 6 -> 24 -> 48 at fixed ~1.45e-3 rad extra tip rotation); never quote it at blade scale;
- the same `_distribute` injects structurally inert spanwise nodal force
  (`sum|fz|` = 24.99% of `sum|fy|`, net zero) that the aero load does not have. Fixing the
  frame is hygiene, not a physics fix, and must not be sold as one.

Ordering is a prerequisite, not a refinement: `q = M/(2A)` is not computable on the stored
node order (the stored ring's shoelace area differs from the angle-sorted one by up to 29% for
133 of 186 raw single-z rings and 171 of 186 merged physical rings).

## Already settled — do not re-open

- MITC4 is not implicated: the coupon discrepancy was a reference bug and the element matches
  converged Rayleigh-Ritz, CCX S8R and free-edge CLT within 1-3% (#9).
- `GJ_static/GJ_deck = 1.068`, plateau 1.03-1.09 over r = 35-80 m. Not a stiffness defect.
- The aero-visible twist is the section's rigid in-plane rotation; `mean(theta_z)` is a
  wall-bending field (verdict 22.8).
- The earlier 0.74x/1.35x torsional deficit was a reference mix-up and is retracted (22.10).

## Merge reconciliation findings (the base is not neutral)

Merging `origin/main` into this branch produced a hybrid where main's guards and this line's
numerical-core variants disagree. Two are measured, and fixing one exposes the other:

- **Mesh frame (fixed, `680cf81`).** This line forced NuMAD `rotorspin = -1` in
  `models/blade/numad/io/yaml_to_blade.py`. Measured on the IEA-15MW blade it mirrors the mesh in
  x (`sum(x)` +1986.78 m against -1986.78 m) and flips the section normal's sense, so
  `ForceProjector` reads it the wrong way: main's guard
  `test_force_projection_load_frame.py::test_load_sense_is_downwind_and_driving` fails with
  `F.y = -8.48e5` N (upwind) where main gives `+8.44e5` N. Removing the override restores main's
  frame; the force-projection / tube / ordering set goes 37 passed.
- **Mesh winding (open).** This line canonicalises element windings
  (`core/mesh/winding.py`, called from `BladeMesh.generate`); main has neither the module nor the
  call. `test_blade_rated_twist.py::test_rated_aero_loads_reproduce_the_bem_resultants` fails with
  the in-plane edge construction (`16.8070%` against a `0.5%` bound). Disabling the
  canonicalisation in-process makes that test pass but breaks
  `test_rated_twist_under_production_loads`, and on pure `origin/main` both pass (7/7). So the two
  lines disagree about the mesh's own connectivity, not about one call site: this is a core
  divergence, not a one-line defect.
- **BEM polars (open).** This line resamples every polar onto a shared `[-180, 180]` degree grid
  (`solvers/bem/engine.py::_build_ccblade_polar_columns`, imported by
  `tests/validation/bem/test_bem_polars.py`); main passes the polar columns through. `bem.Np`
  differs by ~0.4% at rated. Adopting main's engine alone does not fix the winding failure.

Consequence: continuing on this merge means, for every core file both lines rewrote, choosing one
line. main's guards only pass with main's core; this line's core variants (winding canonicalisation,
rotorspin, polar resampling) were built against this line's mesher and BEM. That decision is the
user's, and it is the true prerequisite for #11 T2b/T3 and for #12.

### The structural numbers, measured (user decision: main wins in the core)

`test_blade_rated_twist.py` prints the two quantities #11 and #12 quote. Measured with `-s`:

| tree | tip section `omega` | `distortion/|omega|` | `tip mean theta_z` |
| --- | ---: | ---: | ---: |
| pure `origin/main` | **-1.5112 deg** | **9.3633** | -3.8558 deg |
| this branch (canonicalisation on) | -0.3895 deg | 4.8062 | +1.4813 deg |
| this branch (canonicalisation off) | **+0.1194 deg** | 15.4966 | +2.0090 deg |

These are **pre-fix** (before `1146265`) and measured under the old convention. After the
sign fix the properties-less production path measures `omega = +8.1048 deg` (nose-down) and
the with-properties multi-cell path `+9.6669 deg`, both reported in the state update above.

main's numbers are the ones the issue text quotes, so main's structural core is the reference and
this branch's is not.

The two main guards that matter here trade off against each other on this tree:

- canonicalisation **on**: `test_rated_twist_under_production_loads` passes (`omega < 0`) and
  `test_rated_aero_loads_reproduce_the_bem_resultants` fails (16.8070% against 0.5%);
- canonicalisation **off**: the loads test passes and the twist test fails on
  `omega < 0` (`+0.1194 deg`), the same assertion main satisfies with `-1.5112 deg`.

The `omega < 0` assertions above were the inverted convention and are now `omega > 0`
(issue #11, `1146265`).

So this line's mesh needs the winding canonicalisation for the composite ply-angle mapping to give
the nose-down twist, while main's mesh needs none: the two meshers do not produce the same
connectivity. Swapping main's `models/blade/numad/*` files and main's `core/assembler.py` does not
change it (`omega` stays `+0.1194 deg`), and the blade fixture assembles through the Rust
`PyMeshAssembler` directly, so the divergence is in the mesh-building path this line rewrote
(`core/mesh/generators.py::BladeMesh`, `_deduplicate_and_create_mesh`) and in
`core/mesh/model.py`, not in the element kernels (`mitc4.rs`'s body is identical to main's).

### Mesh/BEM alignment (done)

It did not need `generators.py`'s additive surface re-added: `BladeMesh.generate` was already
identical to main's except for the canonicalisation call, so the alignment was three targeted
changes - drop the winding canonicalisation, take main's BEM polar path, keep the `rotorspin`
removal - and the additive CLI surface (`--export-mesh`, `RotorHubMesh`, `airfoil_spacing`, the CCX
deck) was never at risk. Acceptance met exactly: `omega = -1.5112 deg`, `distortion/|omega| =
9.3633`, 7 passed. (Pre-fix values under the pre-fix convention; the corrected numbers are in the
state update above.)

Still divergent from main, measured and recorded, not needed by #11/#12 and not covered by its
guards: `crates/aeroelast-core/src/elements/mitc3.rs`, `solvers/elasticity/*`, `solvers/fsi/*`,
`core/config.py`, `core/assembler.py`, `models/blade/aerodynamics.py` (this line's prebend/sweep)
and the numad extras.

## Tasks

- [x] T0 — Environment: rebuilt the Rust extension after the merge (`maturin develop
      --release`, GCC 14 on `LD_LIBRARY_PATH`), `_aeroelast` imports and both tube modules
      collect. Evidence: `assemble_m_lumped` present, tube baseline 5 passed, case B = 31.6933x
      Bredt (matches the issue). Commit `26ffe6e` carries the merge.
- [x] T1 — Contour ordering: `src/aeroelast/solvers/bem/section_contour.py`, pure and
      numpy-only — `section_plane_axes` (`u x v = +span`), `order_ring` (angular order about the
      in-plane centroid, deterministic tie-break), `signed_area` (shoelace in that frame),
      `is_simple` (exact orientation predicates, no epsilon), `contour_report` (frozen dataclass:
      closed / simple / area -> `usable` + `reason` + `witness`). Tests:
      `tests/test_section_contour.py`, 9 exact-property tests (no store rows: these are
      properties of our own data, not references). Outcome: 9 passed; tube regression 5 passed.
- [x] T2prev — Mesh/BEM core aligned to `origin/main`, the blocking prerequisite found by
      measurement: the winding canonicalisation (`core/mesh/winding.py`, absent from main)
      removed, main's BEM polar path taken, the `rotorspin` override dropped earlier. Acceptance
      met exactly: `tests/validation/blade/test_blade_rated_twist.py` prints tip section rotation
      -1.5112 deg and distortion/|omega| 9.3633 (main's own numbers), 7 passed. Commits `680cf81`,
      `22d3ccc`, `5520d77`.
- [x] T2a — Torsional realisation as a wall shear flow, pinned on the tube: rate/Bredt 1.0031
      (0.3091%) for the production entry point `realise_section_load`, on top of the case-A guard,
      with the tolerance-free invariances (net force, independent torque ruler, sign). Corrected
      after measurement: the ring-gap tolerance must ride the mesh's span extent (a prebent ring
      was being sheared into partial arcs), and the realisation is gated to an
      exactly-one-usable-ring strip because the multi-ring equal-share flow inverts the blade's tip
      rotation (+0.1273 deg against the minimum-norm -1.4843 deg; pre-fix, before `1146265`).
      Commits `f9d4144`, `7c82e9c`.
- [x] T2b — Store refreshed for group 31: both comparisons are real Bredt rows at 0.3091%, the
      adjudication points at the new row id, the two diagnostic prints are declared unasserted, and
      `gaps.yaml`'s open defect now describes the remaining multi-ring case. `check` 0 errors,
      `coherence` 29/29, `regression --group 31 --write` recorded 2 baselines. Commit `bb19fc4`.
- [x] Verification — main's whole validation suite is green on this tree: `tests/validation`
      434 passed, 18 skipped, 12 xfailed, 0 failed (parity+element 199 passed / 1 xfailed;
      bem+benchmarks+blade+rotor 235 passed / 18 skipped / 11 xfailed).
- [x] T3 — Blade re-measurement: closed by `1146265`. The blade no longer sits at main's
      pre-fix numbers: the properties-less production projector measures tip section rotation
      **`+8.1048 deg`** (nose-down, `omega > 0`) and the with-properties multi-cell path
      **`+9.6669 deg`**; the de-loading table is twist only **`-26.01% / -15.61%`**, radii only
      **`+1.24% / +0.95%`**, production path **`-25.31% / -14.77%`** (thrust/power) against Zhou
      Table 6's `-13.04% / -8.38%`. The multi-ring realisation is its own work unit now, with its
      own independent reference (the hand-assembled Bredt-Batho system in
      `tests/test_multicell_shear_flow.py`). The remaining gap is magnitude (~1.9x over Zhou),
      reported not tuned.

## Structural closure — the convergence analysis (2026-10-05)

The two FE codes must not be compared on the same mesh. MITC4 is a linear 4-node element and
CalculiX's S8R is quadratic, so a same-mesh comparison conflates **element order** with **mesh
size** — the trap `tests/validation/blade/test_blade_iea15mw_mesh_convergence.py` already
documents ("the 4-node-vs-8-node element-order difference ... not discretisation"). Each code is
converged on its own sequence instead, and compared where each has settled.

`tools/ccx_blade_twist_arbitration.py --element-size <h>` runs both codes on the **identical
nodal force vector**, root clamped, with the same section-rotation estimator (`_ring_kinematics`):

| h [m] | AeroElast (MITC4) | CalculiX (S8R) | difference |
| --- | ---: | ---: | ---: |
| 1.00 | +9.6669 | +11.9649 | 23.77% |
| 0.50 | +8.3937 | +10.7368 | 27.91% |
| 0.25 | +8.0980 | +8.4396 | **4.22%** |

At 0.25 m, where the shell's own sequence has settled — increments of 13.2% then 3.5% for a
halving of h, about second order, Richardson limit near **+8.0 deg** — the two independent codes
agree to **4.22%** under the same load. The 23.77-27.91% at the coarser meshes was
discretisation: neither the element formulation nor our load path, and not a number that can be
quoted from a same-mesh run.

**Structural side closed.** Both independent meshes amplify the tip rotation far above the deck's
decoupled sectional beam (1.7765 deg), so the extra twist over a sectional model is real 3D
behaviour; the shell matches an independent FE code to 4% where it is converged; and the tube
(1.0031x Bredt), the two-realisation invariant and the mesh convergence of the gap between the
realisations (1.56 -> 0.78 deg over 1.0 -> 0.25 m) all hold. The residual magnitude against Zhou
is therefore the **aero** side, which is the next phase.

## Constraints

- `src/aeroelast/solvers/bem/force_projection.py` is production; this issue authorises changing
  it, with a test that exercises the production entry point (not a test-local copy of the rule).
- Tests as cheap as the physics allows: one tube solve is seconds; one blade solve is ~40-90 s
  and ~450 MB. Run the full suite **once, serially, alone** (peak 4.26 GB).
- Store rules (`docs/validation-workflow.md`): a bound must be a module-level literal (an
  imported alias reads `UNREADABLE`); one source file per group; one print shape per group;
  adjudications go in `docs/validation/adjudications/<group>-<slug>.yaml`; `extract --write`
  regenerates the row's `justification` and `notes` from the code.
- Test-first by default: observe RED before GREEN for each behaviour, and never invent evidence.

## Evidence pointers

- `tests/validation/parity/test_thin_walled_tube_moment_realization.py` (two-realisation
  comparison, group 31, 0.3091% vs Bredt for the valid realisation).
- `tests/validation/parity/test_thin_walled_tube_torsion.py` (`_ring_shear_flow`,
  `_self_equilibrated_load`, `_realised_torque`, `_bredt_isotropic`, `_theta_fit`, `_theta_z`).
- `odd/tasks/composite-bend-twist-verdict.md` sections 22.7, 22.9, 19.3; commits `40bf36a`,
  `d4fec33`, `e1b747a`; store: `gaps.yaml` id `moment_realization_over_delivers`.

## T2c — the section is multi-cell, and the AC was on the wrong datum

- **T2c-0 (measured).** The blade's section is not a tube: **159 of 186 physical rings form 3
  cells** (outer skin plus two webs); only the 27 circular root rings are single-cell. `q = T/(2A)`
  is the single-cell flow, so on most of the span it is the wrong field.
- **T2c-1 (`0082b75`).** `section_cells` + `cell_adjacency` in `section_contour.py`: the bounded
  faces of the wall graph and which cells each edge bounds. 8 property tests (ring, two-cell box,
  three-cell box = the blade's own topology, dangling stub, non-z span, input order, degenerate
  inputs, immutability).
- **T2c-2a (`d9e8f70`).** The section frame and the acoustic centre now come from each strip's
  **station ring**, not from the whole BEM band. Measured: **12 of 50 strips chose their LE and TE
  on different physical rings** (24%), so the band's "chord" was not any section's; against the
  AeroDyn chord the band datum errs 1.517%/5.538% (median/p90) where the ring datum errs
  0.648%/1.736%. The two guards move to the same datum **without mirroring the implementation**
  (the ring is selected from the deck's own radius), and their bounds are unchanged. Measured
  consequence: the rated tip rotation moved **-1.5112 -> -0.8696 deg** (pre-fix convention) and
  the twist-only de-loading deepened **-4.00% -> -5.13%** of thrust toward Zhou's -13.04%, while
  the magnitude against Zhou's -3.60 deg moved away (the one-way-versus-coupled comparison the
  module already documents). After `1146265` the corrected picture is a nose-down `+8.1048 deg`
  and a twist-only de-loading of **`-26.01% / -15.61%`**, i.e. now over Zhou. 42 + 10 passed;
  no bound widened.
- **T2c-2b (done, `1146265`).** The multi-cell Bredt-Batho realisation: the skin carries
  `q_i`, a shared web carries `q_i - q_j`, the total moment is exact and the net force is zero.
  Validated on an **asymmetric** two-cell box against an independently hand-assembled `(n+1)`
  Bredt-Batho system solved with `numpy.linalg.solve`, with the three-cell blade topology as a
  second case and the production `project()` path reading the shared web back as `q_i - q_j`
  (`tests/test_multicell_shear_flow.py`). `test_box_multicell_torsion.py` already documents that a
  symmetric box cannot detect the error because its web carries zero net flow, which is why the
  pin is asymmetric. The old "single-cell flow inverts the twist" observation was the sign chain,
  not the flow: with the moment applied on `_strip_moment_axis_sign` both realisations are
  nose-down (`+8.1048` vs `+9.6669` deg).
- **T2c-2c (closed, reframed).** CCX arbitration of the blade's load path is replaced by the
  arbitrated sign chain: `tools/diagnose_sign_chain.py` matches every real station ring to its
  WindIO aerofoil (residual 0.5-3.6% of chord) and settles the leading edge, the `_strip_chord_dirs`
  sense and the rotation convention. The two CCX limits measured earlier still hold:
  `write_ccx_mesh` spreads a load **uniformly** over the nodeset (`_write_ccx_cload` divides by the
  node count), so it cannot represent a per-node wall flow; and sending the torque as a moment on
  the 6th DOF blows up through the drilling DOF (418 m), documented in
  `tests/test_blade_ccx_parity.py` - which is this line's own file, not main's. CCX remains the
  **structural** arbiter; the load-path sign is arbitrated by the deck's geometry.
