# Feature: the coupled rotor path diverges at t ≈ 1.9 s

Status: open
Owner: this session (2026-10-07)
Blocks: every roadmap item of issue #18 that needs a production FSI run — the
campaign relaunch, the #13 production measurement, and the smoke gate recorded in
`odd/tasks/force-projection-shear-flow.md` ("Campaign re-run — the gates", gate 1).

## Why this exists

The roadmap says item 1 (#13, the axial artefact) first, and #13's first step is a
measurement on production outputs. There is no healthy production run to measure:
the gate itself is failing.

## Evidence as of 2026-10-07

### E1 — the smoke diverges at t ≈ 1.9 s

Job `11609965` (HEAD with Fix A/C, `tests/run_step1b_smoke.srm`, started 12:42,
cancelled 13:47 after this reading). Outputs in
`$SCRATCH/smoke_fix_results/base_fix/`, `fluid/bem_report.csv`:

| t [s] | Max Disp [m] | Thrust | Power |
| ---: | ---: | ---: | ---: |
| 0.01 | 2.05 | 2.66 MN | 16.5 MW |
| 0.41 | 14.6 | 2.50 MN | 16.5 MW |
| 0.81 | 21.9 | 0.27 MN | −30.5 MW |
| 1.61 | 13.7 | 3.23 MN | 13.6 MW |
| 2.01 | 7.5e14 | 0 | 0 |
| 5.52 (last row) | 5.4e14 | 0 | 0 |

`logs/fluid.log` keeps "completing" windows in 1 iteration because `||F|| = 0`; the
run would have burned its 6 h wall producing nothing. Cancelled.

This **refutes** the recorded conclusion `fsi-coupling-ab-verdict-startup-transient`
("the non-convergence is the initial transient, not a broken coupling"): every A/B
of that session ran 10 windows = 0.1 s, and the process that breaks needs ~1.9 s.
The comparison that matters is "the 2026-09-30 campaign ran 100 s clean" against
"HEAD blows up at 1.9 s".

### E2 — the case is not the difference

`diff` against the campaign's own case file
(`tests/IEA15MW/frontiersin_2025_yaw/solid_corotational_yaw_0.yaml`) shows ONE line:
`airfoil_spacing: constant` (smoke) vs `cosine` (campaign). The fluid yamls differ in
one line: `velocity_data: null` vs `Velocity`, and that is inert — nothing in the
Python rotor path or in `crates/aeroelast-solvers/src/petsc/fsi/` writes a `Velocity`
field to preCICE (grep: only `AngularVelocity` and `Displacement`).

Same `element_size: 0.25`, same `zeta: 0.03`, same `omega_ramp_time: 0` and
`force_ramp_time: 0`, same precice XML except `max-time`. The campaign's own output
mesh (`frontiersin_results_corotational_100s/yaw_0/fluid_mesh.vtu`, 1 065 820 B)
matches the smoke's (1 057 988 B).

So the divergence is either a **code regression** between 2026-09-30 and HEAD, or the
`airfoil_spacing` case difference.

### E3 — one free clue: the extracted outboard twist collapsed

`bem_sectional.csv` per sub-iteration carries the production `r[m]` and `twist[deg]`
that the radii/twist feedback hands to CCBlade. Reference (AeroDyn deck, root→tip):
twist `15.59° → −1.24°`. Comparing the settled window of the campaign against the
first windows of the smoke at the same strip index (53 strips):

| source | t [s] | r(48) [m] | twist(48) [deg] | elastic Δθ (48) |
| --- | ---: | ---: | ---: | ---: |
| campaign (pre-fix, ran 100 s) | 0.5 | 118.65 | +6.76 | +8.0 |
| campaign (settled) | 100 | 117.46 | +6.75 | +8.0 |
| smoke (HEAD) | 0.01 | 118.24 | −3.47 | −2.2 |
| smoke (HEAD) | 0.2 | 118.80 | −4.67 | −3.4 |
| smoke (HEAD) | 0.5 | 118.27 | −0.32 | +0.9 |
| smoke (HEAD) | 1.4 | 109.22 | −16.40 | −15.2 |
| smoke (HEAD) | 1.7 | 116.74 | −11.45 | −10.2 |

`Δθ` here is `twist − ref_twist`, so it carries the flip of `_twist_mesh_to_bem`
(−1 pre-fix, +1 after `1146265`); the low-magnitude rows are comparable, the violent
ones are not a convention artefact.

Reading: the corrected production loop runs with an outboard elastic twist that starts
at zero/negative and swings to −15°, where the corrected **static fixture**
(`tests/validation/blade/test_blade_deloading_vs_reference.py`, group 27) measures
`+8.1048°` nose-down and de-loading. The loop keeps the rotor loaded (2.5–2.7 MN
against the 2.541 MN rigid anchor) instead of de-loading, drives the blade past 21 m
of tip deflection at t = 0.8 s, and then diverges.

## The defect (found 2026-10-07, same session)

`BEMFSIParticipant` derived its mesh-rotation -> BEM-twist factor from
`normal . (span x tangential)`.  `tangential_direction` configures the sense of the
tangential **force** `Tp` (`F = Np*n + Tp*t`), and it was doing double duty as a proxy
for the chord axis.  The static fixtures never pass the key (default `+x`), every
production case YAML (62 files) sets `[-1, 0, 0]`, so the same expression evaluated to
`+1` in the fixtures and `-1` in production.

`1146265` ("three linked sign errors") inverted the polarity of that condition while
keeping the condition, so it took production from `+1` (the value its own new comment
calls correct) to `-1` (the value its own new comment calls the `+21.88%` re-loading).
The fix was validated only in the fixtures' frame, which is why its assertions passed.

Measured, on the real blade mesh and the real AeroDyn deck, with the same rigid nose-down
section rotation applied to both frames (`tools/diagnose_twist_sense_frame.py`):

| frame | `_twist_mesh_to_bem` | rigid +2 deg nose-down | thrust |
| --- | ---: | ---: | ---: |
| fixtures (no `tangential_direction`) | +1 | `twist_def - ref = +2.0 deg` | **-13.05%** |
| production (`tangential_direction = -x`) | -1 | `twist_def - ref = -2.0 deg` | **+10.63%** |

The physical indicator is identical in both frames: the leading edge sits at `+x` and the
load frame puts downwind at `+y`, so a positive rotation about `+span` is nose-down and
CCBlade's `alpha = phi - theta` requires `theta` to RISE.  The factor is `+1` in both.

The coupled consequence, read off `bem_sectional.csv` (the production `twist[deg]` handed
to CCBlade, less the deck's reference twist): the smoke's outboard elastic twist ran
`-2.2 -> -15.2 deg` (nose-up = loading) where the 2026-09-30 campaign - which de-loaded and
ran 100 s - held `+8.0 deg`.  `+10.63%` of thrust per 2 deg of twist is the runaway, and it
is why the path broke at `t ~ 1.9 s`, one flap period in.

Why the pre-fix campaign was healthy: it carried the inverted projecting moment axis
*and* this frame, so the two cancelled.  `1146265` corrected the moment axis and left the
frame inverting, so only the error remains.

## The fix

`_resolve_twist_sense()` reads the leading edge from what `ForceProjector` already
measures: `_section_ends` identifies it on each strip's station ring and
`_strip_moment_axis_sign` records whether the strip's chord axis runs leading-to-trailing,
arbitrated against the deck's own WindIO aerofoils.  The leading-edge direction is
`-_strip_moment_axis_sign * _strip_chord_dirs`, width-weighted, and the factor is `+1`
when it has a positive component along `normal x span`.  That is the same class of fix as
the projecting-moment one: measure the section frame, do not infer it from a load sense.
The degenerate case (no measurable chord axis) keeps the configured-frame expression,
documented as a fallback rather than an arbitrated answer.

One stale assertion fell out of the same tangle:
`test_bem_fsi_deformed_geometry.py::test_elastic_twist_matches_the_nodal_section_rotation`
asserted `nose-down = negative about +span` and was **failing on HEAD** next to the
`+8.1048 deg` it measures, contradicting both the deck's geometry and
`test_blade_deloading_vs_reference.py`'s declared convention.  Corrected with the
reasoning recorded in its docstring.

Also corrected, and **unrelated to this defect**: the `# pyright: ignore` on the
`meshio.Mesh(...)` call in `_write_sections_vtu` sat on the closing paren while pyright
reports the two argument lines, so it suppressed nothing (verified identical at HEAD,
lines 1055/1056).  Moving it onto the reported lines clears the file.

## T3 result — the 5 s run (job `11610007`, 2026-10-07)

The sign fix is confirmed live; **a second, independent defect is now the blocker.**

What the fix changed (same case, `max-time 5.0`, past the 1.9 s breakdown):

| | factor `-1` (`base_fix`) | factor `+1` (`twistfix_5s`) |
| --- | ---: | ---: |
| windows completed | 500 | 500 |
| `max disp` at `t ~ 2 s` | 7.5e14 m (permanent runaway) | 21 m, returns to 11.5 m |
| elastic tip twist, last second | never settles | **+4.5 deg, nose-down** |

What did **not** change: the implicit coupling does not contract.

| run | windows | converged | iters mean | max |
| --- | ---: | ---: | ---: | ---: |
| campaign 2026-09-30 (`frontiersin.../yaw_0`) | 500 | **500 (100%)** | **4.1** | 18 |
| `twistfix_5s` (today, factor `+1`) | 500 | 51 (10%) | 28.5 | 30 |
| `base_fix` (today, factor `-1`) | 500 | 328 (66%) | 12.1 | 30 |

The `base_fix` 66% is an artefact: once the structure has run away `||F|| = 0` and every window "converges" in 2 iterations. Its own fault, already fixed.

The `twistfix_5s` transient is still violent: `max disp` peaks at 2783 m at `t = 3.0`,
thrust swings `-0.12 -> +3.47 MN`, and the first 20 windows all sit at the 30-iteration
ceiling.  The outboard twist reads `-4.8, -9.7, +0.1, -82.4, +139.9, +4.5 deg` along the
run: the wild values are mid-iteration states of windows that never converged, and the
settled ones carry the correct nose-down sign.

So: one defect down, one to go, and the second one is a **regression** by construction.
The campaign that ran 100 s at 100% convergence predates the 2026-10-04 merge and the
#11 batch, and its own case file differs from the smoke by exactly one line
(`airfoil_spacing` on the SOLID mesh: `cosine` vs `constant`).

## Tasks

- [x] T1 Cancel the diverged smoke (11609965) and record why.
- [x] T2 Offline: name the defect and prove it, `tools/diagnose_twist_sense_frame.py`
      (`-1` / `+10.63%` production against `+1` / `-13.05%` fixtures).
- [x] T2b Fix `_resolve_twist_sense`, RED first: the new guard test failed with
      `factor = -1, +10.63%` before the change and passes with `+1, -13.05%` after.
- [x] T3 Live 5 s (job `11610007`): runaway gone, twist sense correct, coupling still
      does not contract.
- [x] T4 Separate case from code for the non-contraction: job `11610151` runs the
      campaign's own case (`airfoil_spacing: cosine` on the solid mesh, `max-time 5.0`).
      **Done 2026-10-07: it does not converge.** 34 of its 493 windows closed under the
      30-iteration ceiling (6.9%) against the gate smoke's 51 of 500 (10%), the rest sat on the
      ceiling, `thrust ~= -0.7 MN` and not settling at `t = 4.93 s`. The case is ruled out, so
      the non-contraction is a code regression, and the next rung is a worktree at `552566d`
      (last pre-merge state) on the same case.
- [ ] T5 Re-run the 30 s gate (`tests/run_step1b_smoke.srm`) once T4 is clean, and
      re-anchor the campaign baseline.

## Follow-ups, out of scope here

- The 62 production case YAMLs still carry `tangential_direction = [-1, 0, 0]` for the
  force sense.  That is correct for `Tp` and must stay; the twist sense no longer reads
  it.
- The `_twist_mesh_to_bem` name is now a misnomer for a measurement, but renaming it
  touches the deloading fixtures' assertions; own change.
- `tests/validation/blade/test_blade_deloading_vs_reference.py` builds the participant
  without `normal_direction`/`tangential_direction`, so it validates the default frame.
  It should run the production frame too, which is what the new guard does for the sign.
- The `# pyright: ignore` on `meshio.Mesh(...)` in `_write_sections_vtu` sat on the
  closing paren while pyright reports the two argument lines (identical at HEAD, 1055/1056);
  moved onto the reported lines, which clears the file.

