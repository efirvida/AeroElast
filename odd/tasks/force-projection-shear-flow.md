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
- [ ] T2 — Torsional realisation: realise the section moment's span component as a wall shear
      flow in the production path (`_distribute` or a new function the projector calls), on top
      of the existing case-A guard. Pin it on the tube where Bredt is exact: rate/Bredt within
      5% for the production path. Assert the invariances without tolerance (ruler read back from
      the applied forces, self-equilibration, sign) and never widen a bound.
- [ ] T3 — Blade re-measurement: re-measure with the change and report without fitting — tip
      section rotation, spanwise `sum|fz|` (in-plane part expected ~0), `distortion/|omega|`
      (was 9.36 under production, 0.206 in the older hand-built case), and the de-loading table.
      A result outside the bound is the finding.

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
