# Radii feedback: linear kinematics with the wrong-sign second-order term (#33 L5)

Store gap: `radii_linear_kinematics_axial_term` (`docs/validation/gaps.yaml`, level L5).
Predecessor, closed 2026-10-11: `radii_datum_definition_bias` (`d61993a`, `b155c95`; closure
section in `docs/validation_closures.md`).

## State (2026-10-11)

`_compute_deformed_geometry` feeds `r_def[k] = _ref_r[k] + mean(u[strip k] . s)`. On the
de-loading fixture (`tests/validation/blade/test_blade_deloading_vs_reference.py`: official
IEA-15MW mesh, `assemble_k()` only, aero load only, no `K_G`, no centrifugal), tip strip:

| term | value [m] |
| --- | ---: |
| measured mean `u_z` (what is fed) | +1.1056 |
| linear prebend/sweep kinematics `-int y0' u_y' dz` | +1.0973 (99.2 %) |
| axial strain residual | +0.0083 |
| omitted second order `-1/2 int u_y'^2 dz` | -2.1534 |
| flap `u_y` | +17.007 |

Ratio second-order / linear = 1.96, so the fed term has the wrong physical sign **on this
state**. The verdict holds only while the production deflection is above about half of this one.

## Next session: T1 (decides whether the defect is real)

Repeat the decomposition on the **production structural state** at rated: stress-stiffened `K`
(`K + K_G` from the static centrifugal solve, `assemble_geometric_stiffness` with `free_dofs`)
and the centrifugal body load plus the participant's aero load, rotor at 7.56 rpm.

- Reuse the fixture's integration (trapezoid along strip centroids: `lin = -cum(y0' uy' + x0' ux')`,
  `nl = -1/2 cum(uy'^2 + ux'^2)`), measured against `uz`.
- Report tip `u_y`, `uz`, `lin`, `nl`, and the ratio `|nl| / lin`.
- Outcome A (`|nl| > lin`): the defect is confirmed, go to T2.
- Outcome B (`|nl| < lin`): record the number, downgrade the gap to `bounded` with the ratio as
  its bound, and stop.

## T2 (only if T1 confirms)

Pick a fix and its reference:
- second-order kinematic correction of `r_def` (add `-1/2 int |u_perp'|^2`), plus the
  `Rtip = max(r_def) * 1.001` rule in `_rebuild_bem_solver`, which inherits the same term; or
- a geometrically nonlinear structural state for the one-way path.
Reference: BeamDyn/OpenFAST projected tip position at the same operating point.

## Also open, not in this task

- #19's 500/500 coupled convergence was measured with the old `r_def`; it needs a SLURM re-run.
- The L5 arbiter that matters: a coupled rated run against OpenFAST coupled (roadmap
  `odd/tasks/validation-review-roadmap.md`, L5 rows). The twist, not the radii, drives the
  -26 % vs Zhou's -13 %, and Zhou is non-transferable (#14).
- Store drift, not bisected: tip `dr_proj` 1.106156 -> 1.105628 m and `dr_path` 2.112061 ->
  2.116631 m between `7bd497e` and `db181e2`.
