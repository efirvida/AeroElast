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

## T1 result (2026-10-11): Outcome A, the defect is real on the production state

`tools/diagnose_radii_kinematics.py` (revision `e6b971c`; `_aeroelast` built 2026-10-04 13:24,
after the last `crates/` commit at 13:12, so not stale):
`scripts/aeroenv.sh python -W ignore tools/diagnose_radii_kinematics.py` (~24 s).
Production state: `(K + K_G + K_SP) u = F_aero + F_cf(X0)` at 7.56 rpm, axis `+y`, centre at the
origin as in every case YAML, `K_G` from the static centrifugal solve with `free_dofs`, `K_SP`
lumped as `tools/run_s4_rotating_modal.py`. Rigid-blade rated aero load. Tip strip
(`z0 = 116.403 m`), values in m:

| state | `u_y` | `uz` fed | `lin` | `nl` | `\|nl\|/lin` |
| --- | ---: | ---: | ---: | ---: | ---: |
| control: `K`, aero only (the fixture) | 17.007 | +1.1056 | +1.0973 | -2.1534 | 1.962 |
| **prod: `K+K_G+K_SP`, aero + `F_cf`** | 15.931 | +1.0389 | +1.0278 | -1.8931 | **1.842** |
| prod operator, aero only | 15.568 | +1.0114 | +1.0034 | -1.8063 | 1.800 |
| prod operator, `F_cf` only | 0.364 | +0.0275 | +0.0243 | -0.0010 | 0.042 |

- The control reproduces the fixture to four decimals, so the integration is the same one.
- The ratio is above 1 along the whole span on the production state (z0 57 / 88 / 105 m:
  2.64 / 1.72 / 1.76).
- Stress stiffening lowers the tip flap only 6 %. The centrifugal stretch adds +0.0028 m of
  axial residual, which is negligible.
- The ratio scales with the flap deflection (1.962 x 15.93/17.01 = 1.838 against the measured
  1.842). The verdict therefore reduces to one number: **the fed term has the wrong sign whenever
  tip `u_y` > ~8.7 m** (15.93 / 1.84).
- Independent anchors sit far above that crossover. S-5's ElastoDyn rated OoP is 16.10 m, and the
  shell under the AeroDyn loads gives 15.38 m (`tests/test_iea15mw_s5_oneway.py` docstring).
  Gravity and the twist feedback's de-loading are omitted. Gravity is mostly edgewise, so it
  only raises `u_x'^2`. To bring `u_y` under 8.7 m the de-loading would have to cut the flap
  load by ~45 %. Neither omission can flip the outcome.
- The gap stays `open_defect`; its evidence now carries the production number.

## Next session: T2 (T1 confirmed)

Pick a fix and its reference:
- second-order kinematic correction of `r_def` (add `-1/2 int |u_perp'|^2`), plus the
  `Rtip = max(r_def) * 1.001` rule in `_rebuild_bem_solver`, which inherits the same term; or
- a geometrically nonlinear structural state for the one-way path.
Reference: BeamDyn/OpenFAST projected tip position at the same operating point.

## Also open, not in this task

- Hub radius: blade-mesh `z = 0` is the blade root and every case YAML puts `rotation_center`
  at the origin, so `F_cf` and `K_G` omit the hub-radius arm (production does the same, so T1
  is faithful to it). That under-stiffens slightly; S-4 still closes at +0.6 % on 1st flap.
  This is a separate item, not this gap.
- #19's 500/500 coupled convergence was measured with the old `r_def`; it needs a SLURM re-run.
- The L5 arbiter that matters: a coupled rated run against OpenFAST coupled (roadmap
  `odd/tasks/validation-review-roadmap.md`, L5 rows). The twist, not the radii, drives the
  -26 % vs Zhou's -13 %, and Zhou is non-transferable (#14).
- Store drift, not bisected: tip `dr_proj` 1.106156 -> 1.105628 m and `dr_path` 2.112061 ->
  2.116631 m between `7bd497e` and `db181e2`.
