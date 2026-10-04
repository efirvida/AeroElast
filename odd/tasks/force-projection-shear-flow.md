# ODD — force projection: realise a section moment as a wall shear flow (#11)

Base: `integrate/origin-main-2026-09-30` after merging `origin/main` (`26ffe6e`).
Issue: https://github.com/efirvida/AeroElast/issues/11 (bug).

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

main's numbers are the ones the issue text quotes, so main's structural core is the reference and
this branch's is not.

The two main guards that matter here trade off against each other on this tree:

- canonicalisation **on**: `test_rated_twist_under_production_loads` passes (`omega < 0`) and
  `test_rated_aero_loads_reproduce_the_bem_resultants` fails (16.8070% against 0.5%);
- canonicalisation **off**: the loads test passes and the twist test fails on
  `omega < 0` (`+0.1194 deg`), the same assertion main satisfies with `-1.5112 deg`.

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
9.3633`, 7 passed.

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
      rotation (+0.1273 deg against the minimum-norm -1.4843 deg). Commits `f9d4144`, `7c82e9c`.
- [x] T2b — Store refreshed for group 31: both comparisons are real Bredt rows at 0.3091%, the
      adjudication points at the new row id, the two diagnostic prints are declared unasserted, and
      `gaps.yaml`'s open defect now describes the remaining multi-ring case. `check` 0 errors,
      `coherence` 29/29, `regression --group 31 --write` recorded 2 baselines. Commit `bb19fc4`.
- [x] Verification — main's whole validation suite is green on this tree: `tests/validation`
      434 passed, 18 skipped, 12 xfailed, 0 failed (parity+element 199 passed / 1 xfailed;
      bem+benchmarks+blade+rotor 235 passed / 18 skipped / 11 xfailed).
- [ ] T3 — Blade re-measurement: the change does **not** reach the blade (gated), so the blade sits
      at main's numbers. T3 owes the recorded finding plus what the issue asked for: tip section
      rotation, spanwise `sum|fz|`, `distortion/|omega|` (9.3633 under production), the de-loading
      table, and the multi-ring realisation as its own work unit with its own reference. A result
      outside the bound is the finding.

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
