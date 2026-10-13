# Feature: the coupled rotor path diverges at t ≈ 1.9 s

> **Handoff 2026-10-08 — #19 and #26 are CLOSED and pushed.** Read this block first;
> the rest of the file is the ladder and stays as the record.
>
> **Closed.** #19 (the coupling did not contract) and #26 (the BEM pair rode the
> section, not the rotor plane). Fixes: `6765633` (the projection keeps the
> reference geometry; the BEM loads still follow the deformation) and `7a84da1`
> (the load frame); guards `e3278ba` and the rewritten
> `tests/validation/bem/test_force_projection_load_frame.py`. Verified on the
> campaign's own case at `max-time 5.0`: **500/500 windows, mean 2.60 sub-iterations**
> against the 2026-09-30 campaign's 4.08. Ladders, the consequence table and the
> known trade-off are in `docs/validation_closures.md`, section "#19 y #26"; the
> store and the roadmap (#18) are reconciled.
>
> **In flight when this handoff was written.** The 30 s gate,
> `tests/run_step1b_smoke.srm`, job `11611263`, writing to
> `$SCRATCH/smoke_fix_results/base_fix/` (`fluid/bem_report.csv` carries thrust and
> tip disp; `precice-Solid-iterations.log` carries the windows). At 277 windows it
> read 277 converged (100%), mean 3.23. When it lands: record the steady thrust
> against the rigid anchor 2.541 MN and the window statistics, then close T7.
>
> **Left open on purpose.** (a) The campaign re-anchor — unblocked, but it is a long
> cluster run (5 cases x 100 s) and the user's decision; the new baseline carries
> two declared moves (load frame, projection geometry) versus the stored snapshots.
> (b) `reference_projection_geometry_approximation`, the residual gap the fix
> introduces: the within-strip distribution and the AC arm are read from the
> reference geometry. The fidelity follow-up is to make the `#11` derivation
> deformation-robust (hysteresis-free ring selection, LE/TE tie-break, chord-axis
> continuity) instead of frozen — the campaign's own projector followed the
> deformation and contracted, so a smooth derivation exists.
>
> **Instrument notes that cost time, so they are not rediscovered.** The rung case
> lives in `tests/smoke_fix/frame_ab{,_abs}/` (the second one spells `config_file`
> absolutely because the `origin/main` side resolves it against the workdir while
> HEAD resolves it against the yaml's directory). Shadowing an old Python against
> the venv `_aeroelast` (2026-10-04) needs two shims: the five K_G-deformed /
> omega-rebuild kwargs (`rotor.py`) and a `len(rotor_perf_tuple) >= 6` fallback. The
> worktrees' shims are labelled `DIAGNOSTIC` and are not for commit. Startup varies
> 5-45 min with node contention: budget 3.5 h and cut a run yourself as soon as
> about 20 windows show the pattern (broken = 30 from window 2; healthy = 3-15).

Status: #19 and #26 closed 2026-10-08; the 30 s gate was in flight
Owner: this session (2026-10-07, resumed 2026-10-08)
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
`#11` batch, and its own case file differs from the smoke by exactly one line
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

---

## The regression is the projector's load frame (2026-10-08)

**The issue's ladder could never have found it.** All three named rungs
(`f9d4144` 2026-10-04, `5520d77` 2026-10-04, `638965f` 2026-10-05) *postdate* the
suspect; only rung 1 (`552566d`) brackets it, and by three days.

**The campaign revision is `b5d369e`, not `552566d`.** Job `11604984` (the
`yaw_batch_corotational` batch that produced
`$SCRATCH/frontiersin_results_corotational_100s`) started 2026-09-30 23:37.
The reflog puts HEAD at `b5d369e` (23:16) then: `552566d` is its ancestor two
commits back, and both intervening commits are campaigns tooling. The 100 s
campaign therefore ran on `b5d369e`'s Python, and it converged 500/500.

**The two earliest post-campaign commits on the BEM/FSI path are `6064aa8`
(10-02 09:40, AC datum and moment arm) and `9a3923e` (10-02 15:31).**
`git log --reverse b5d369e..HEAD -- src/aeroelast/solvers/bem/` lists them
before every rung the issue names.

### `9a3923e` is the defect, and it is #26

> **Status 2026-10-08, afternoon: partially refuted.** The fix below landed and is
> unit-verified, and it is the correct physics for #26 - but it does **not**
> restore the contraction. Job `11610819` (`headfix`, the repo's fixed
> `ForceProjector` on the campaign's own case at `max-time 5.0`) reads window 1
> = 15 iterations converged, window 2 onward pinned at the 30-iteration ceiling
> (43 windows in, all saturated), against HEAD's `34/493` (job `11610151`) and
> the campaign's `15 11 11 13 13 ...`. So `9a3923e`'s load frame was **not** the
> coupling regression; the load frame is still wrong (see the measurements
> below) and is fixed, but #19 has another cause. The rungs that localise it are
> in flight (see T6).
>
> The topology also corrects the framing: `b5d369e` (the campaign) is **not** an
> ancestor of `ac2e9e8`/`9a3923e`. They are the two sides of the 2026-10-04 merge
> `26ffe6e` - the campaign is the `integrate/...` side (with
> `552566d`), the BEM batch is the `origin/main` side (`5f189fa`...). So the
> regression is *what the merge imported*, exactly as the issue said, and both
> sides carried their own implementation of the same BEM deformation feedback.

Its own message asserts the premise that #26 later refuted:

> `Np`/`Tp` come from ccblade's `distributedAeroLoads` and are normal/tangential
> to the section chord at that azimuth (`engine.py:210`)

and on that premise it moved the applied pair from the configured fixed frame

    F_strip = F_n * _normal_dir + F_t * _tangential_dir      # b5d369e

to each strip's own section frame

    F_strip = F_n * _strip_normal_dirs[k] + F_t * _strip_chord_dirs[k]   # 9a3923e

Those two frames differ by the local section twist. The section frame is the
rotor frame rotated by `theta` about the span, so from `9a3923e` on, every
applied blade load was **the physical force vector rotated by the local twist**
- a load direction that tracks the deformation the coupled loop iterates on.
That is positive feedback in the partitioned fixed point: it changes the
interface Jacobian the IQN-ILS scheme factorises, which is exactly what a
contraction regression looks like from the outside, while every static fixture
(the untwisted tube of group 35, the band-wide frames of group 27) is blind to
it.

Issue #26 had already settled the reading, from ccblade's source and from our own
output: `cn = cl*cos(phi) + cd*sin(phi)` with `phi` measured **from the rotor
plane**, and the identity
`atan2(Tp, Np) + atan2(cd, cl) - twist = alpha` holds to `max |residual| =
1.4e-14 deg` over the 50 stations (`tools/diagnose_zhou_tp_frame.py`). The
section reading is off by up to `15.6 deg` - the root twist.

### The fix

`ForceProjector` builds the rotor-plane frame once from the configured pair and
applies the BEM pair on it:

    self._rotor_tangential_dir, self._rotor_normal_dir = self._load_frame(self._tangential_dir)
    F_strip = F_n * self._rotor_normal_dir + F_t * self._rotor_tangential_dir

and `verify()` recomputes with the same frame. The per-strip section frame stays
what it always was physically: the datum for the aerodynamic centre and the
pitching-moment axis. `_load_frame` already orthonormalised the configured pair
against the span and kept the blade-wide sense, so no new geometry rule was
invented.

### RED, measured on the real mesh (element_size 1.0, rated point)

`tests/validation/bem/test_force_projection_load_frame.py` rewritten to the
rotor-plane reading (it used to assert the section one). With the fix reverted
the new guards fail; the printed witness is the section normal against the rotor
normal, per station:

| frac | angle(f_hat, z_rotor) | pre-fix `|F.t_rotor|/|F|`, `Np` only | post-fix |
| ---: | ---: | ---: | ---: |
| 0.03 | 22.245 deg | 0.3918 | 0.0000 |
| 0.15 | 11.647 deg | 0.2096 | 0.0000 |
| 0.25 | 7.954 deg | 0.1425 | 0.0000 |
| 0.50 | 2.257 deg | 0.0367 | 0.0000 |
| 0.75 | 1.256 deg | 0.0234 | 0.0000 |
| 0.90 | 2.166 deg | 0.0378 | 0.0000 |
| 1.00 | 1.176 deg | 0.0205 | 0.0000 |

Inter-station section-normal spread `24.410 deg` (frac 0.03 vs 0.90): the mesh
discriminates the two readings. Uniform `Np` conservation: `0.7547%` pre-fix
(the section normals fan out), `0.0000%` post-fix. `Tp` mirrors it exactly.

### The consequence at rated (same mesh and deck, issue #26 step 3)

| application | rotor-axis force | in-plane force | edgewise root moment |
| --- | ---: | ---: | ---: |
| rotor plane (HEAD + fix) | 851 437.8 N | 103 983.5 N | 6.9503 MN.m |
| section frame (pre-fix) | 844 453.9 N | 117 306.3 N | 6.6847 MN.m |
| delta | **+0.83%** | **-11.36%** | **+3.97%** |

Measured with both applications decomposed on the *same* rotor frame and the
`blade_aero.r` lever arm, so the comparison is physical. The signs of the
deltas disagree with #26's estimate (which set the sign convention aside); the
magnitudes agree, and this is the measured pair.

### Tasks

- [ ] T5 Establish the instrument: worktrees at `ac2e9e8` (`9a3923e^`) and
      `9a3923e`, plus `shadow-headfix`; all three on the venv `_aeroelast` build
      so the Rust side is constant across rungs (the window's only Rust changes
      are MITC3-only - `crates/aeroelast-core/src/elements/mitc3.rs` and
      `smoothing.rs`, and the blade is MITC4).
- [x] T6 The anchor holds. `campaign` (`b5d369e`, job `11610862`) on the
      consensus case `tests/smoke_fix/frame_ab/` at `max-time 5.0`: **423 of 423
      windows converged**, first 12 = `15 11 11 13 13 12 11 10 11 10 11 10` -
      the campaign's own sequence. So the 5 s instrument reproduces the healthy
      baseline and the regression is in the code between the campaign and HEAD.
      Case construction: the campaign's own yaml with the 4 rotor keys whose
      values equal HEAD's defaults removed, so every revision in the window
      parses it; `config_file` and the mesh paths are the case's own.
- [ ] T6a **The `origin/main`-side rungs cannot run against the venv Rust.**
      `ac2e9e8`, `9a3923e`, `pre4fee` and `fee690` all die at start-up with
      `run_rotor_fsi_solver() missing 5 required positional arguments:
      omega_rebuild_rel_high, omega_rebuild_rel_low, kg_use_deformed_coords,
      kg_deflection_rebuild_rel_high, kg_deflection_rebuild_rel_low`. Those five
      are the K_G-deformed / omega-rebuild plumbing: the **local** line has it
      (`b5d369e` passes 9 kwargs, HEAD passes 9) and `origin/main` at 10-02 does
      not (0 kwargs matched). The venv `_aeroelast` is the 2026-10-04 build, so
      it wants the newer call. Bisecting that line needs one Rust rebuild per
      rung; not attempted. The four jobs were cancelled (`11610844`, `11610845`,
      `11610863`, `11610864`) because the Python shadow cannot drive the newer
      Rust.
- [x] T6b **The deformation feedback is not the cause either.** `nofeed` (job
      `11610901`), HEAD with both feedbacks frozen, reads **144 of 500 windows
      (28.8%)**, first 12 = `12 30 30 30 30 30 30 30 30 30 30 30`: window 1 is
      close to the campaign's 15, then the ceiling. Freezing the feedback buys
      `28.8%` against HEAD's `7.8%` and `headfix`'s `8.8%`, so it *contributes*
      but three quarters of the windows still saturate. The defect is upstream
      of the geometry feedback. (The `noradii`/`notwist` attribution rungs are
      therefore not worth a slot yet.)
- [x] T6f **Both sides of the merge are healthy on their own.** The same
      instrument, the same (absolute-path) case, 1 h budget:

      | tree | windows | converged | first 8 |
      | --- | ---: | ---: | --- |
      | `campaign` = `b5d369e` (local side) | 500 | **500 (100%)** | (job 11610862), run to 5 s |
      | `origmain` = `5f22f51` (merge's 2nd parent) | 206 | 185 (89.8%) | 30 8 4 4 5 4 5 5 |
      | `ac2e9e8` (origin/main @ 10-02) | 243 | 232 (95.5%) | 30 9 10 8 6 4 3 3 |
      | `fee690` = `4fee690` | 217 | 192 (88.5%) | 30 9 7 7 7 5 5 8 |
      | HEAD (`7a84da1`) | 500 | 39 (7.8%) | 23 30 30 30 30 30 30 30 |
      | `headfix` (was `7a84da1`) | 159 | 14 (8.8%) | 15 30 30 30 30 30 30 30 |
      | `nofeed` (feedback frozen) | 500 | 144 (28.8%) | 12 30 30 30 30 30 30 30 |

      So **the regression is the merge's own combination**, not a single commit
      on either side: each parent contracts, HEAD does not. The origin/main side
      does saturate its *first* window (30) and then settles to 3-12, while the
      campaign opens at 15 and never saturates.
- [x] T6g Where the search is now, from file-level ancestry at the merge. The
      merge's **own** resolutions in the coupled path are only
      `core/config.py`, `core/mesh/__init__.py`, `core/mesh/generators.py` and
      two MITC3 Rust files. `config.py`'s only substantive resolution is the
      `normal_direction`/`tangential_direction` **default** pair (origin/main's
      transposed one, `[0,1,0]`/`[1,0,0]`) - and the case sets both explicitly,
      so it is inert here. The mesh is **identical** across every rung by count
      (32336 nodes / 33473 elements solid, 32325 / 33462 fluid), so the
      merge-own `generators.py`/`__init__.py` differences are formatting and
      dead code, not geometry. And `solvers/fsi/rotor.py`, `corotational.py`,
      `bem/engine.py`, `mesh/model.py`, `mesh/winding.py` are **byte-identical
      to the campaign's** in the merge, so they are out. What is left is exactly
      two files, both taken from origin/main and both rewritten relative to the
      campaign: `solvers/bem/fsi_participant.py` (+389/-179) and
      `solvers/bem/force_projection.py` (+151/-285), plus `standalone.py`
      (2 lines).
- [x] T6h **Found it: the lever is the projector rebuilt on the deformed mesh.**
      Job `11611087` (`projfrozen`, HEAD but keeping the **reference** projector
      instead of rebuilding it per sub-iteration) reads **94 of 94 windows
      (100%)**, first 12 = `11 7 8 11 10 8 7 9 9 8 7 8` - as healthy as the
      campaign. Job `11611088` (`bemhead`, the healthy `origmain` tree carrying
      HEAD's `fsi_participant.py` + `force_projection.py` + `standalone.py`)
      reads **2 of 24 (8.3%)**, i.e. it fails exactly like HEAD. So the defect
      lives in HEAD's fluid files, and the mechanism is the applied load's
      geometry following the deformation.
      But the rebuild itself is not the bug: the **campaign rebuilds the
      projector on the deformed mesh too**
      (`b5d369e:fsi_participant.py:1045`, same call, same arguments) and it
      converges 500/500. What breaks is *this* projector's constructor reacting
      to the deformation - the `#11` rewrite of `force_projection.py`
      (+151/-285 against the campaign's).
- [x] T6i Attribution done - **no single geometric output is the lever; the
      sum is.** All rungs are HEAD at `7a84da1` with the rebuild kept and one
      attribute pinned to the reference projector afterwards:

      | pin | windows | converged | first 12 |
      | --- | ---: | ---: | --- |
      | none (HEAD) | 500 | 39 (7.8%) | 23 30 30 30 30 30 30 30 |
      | `_strip_moment_axis_sign` (`11611101`) | 6 | 1 (17%) | 15 30 30 30 30 30 |
      | `_strip_ac_offsets` (`11611102`) | 18 | 8 (44%) | 20 20 30 30 30 30 24 17 18 12 30 30 |
      | `_strips` (`11611103`) | 18 | 2 (11%) | 18 12 30 30 30 30 30 30 |
      | `_strip_ring_groups` + `_strip_ring_sections` (`11611116`) | 24 | 1 (4%) | 15 30 30 30 30 30 |
      | **all of them** (`11611119`) | 39 | **39 (100%)** | 11 7 8 11 10 8 7 9 9 8 7 8 |

      So the coupling has no margin to spare: each deformation-dependent piece
      of the projected geometry contributes a small gain, none of them alone
      crosses the threshold, and their sum does. `frzall` (100%) is exactly
      `projfrozen` (100%) in effect, which also says the rebuild's *other* side
      effects are irrelevant - only the geometry values matter.
      Note the two main pieces are the *distribution* geometry: the strip grid
      and the ring groups/sections (`realise_section_load`'s inputs), which is
      why pinning either one alone leaves the other free. Job `11611231`
      (`frzdist`) tests exactly that pair, leaving the AC arm and the moment
      sign free, to fix the minimal shape of the fix.
- [ ] T6e Those three (`11610978`-`11610980`) returned **zero windows in 3 h**,
      and it was not the coupling: the **fluid** participant died with
      `RuntimeError: XML parser was unable to open configuration file
      "<RUN_DIR>/precice-config.xml"` and the solid sat in preCICE's
      "Setting up primary communication" handshake until the timeout. The
      cause is a path-resolution difference on the origin/main side:
      `run_bem_fsi.py:331-333` resolves `config_file` against the **workdir**
      (`cfg["config_file"] = str(workdir / cfg["config_file"])`), while HEAD
      resolves it against the **yaml's directory** (`_resolve`).
      `tests/smoke_fix/frame_ab_abs/` is the same case with `config_file`
      spelled absolutely - verified that only that line differs from
      `frame_ab/` - and the three rungs were resubmitted with it
      (`11611050` `origmain`, `11611051` `ac2e9e8`, `11611052` `fee690`). The
      healthy rungs already recorded (`campaign`, `nofeed`, `headfix`) ran with
      the relative spelling; the XML file they loaded is the same one, so
      nothing but the path spelling changes.
- [ ] T6c **The `exp_feedback_off` yamls do not work at HEAD** (checked
      2026-10-08): `deformed_twist`/`deformed_radius` were declared in
      `BEMConfig` by `d571917` on the *local* line, and the merge kept the
      `origin/main` `_compute_deformed_geometry`, which applies both feedbacks
      **unconditionally** - nothing under `src/` reads either key, so the yaml
      values are inert. The working substitute is the diagnostic worktree
      `~/fem-shell-nofeed` (at `7a84da1`) plus `~/shadow-nofeed`: it returns
      `(_ref_r, _ref_twist)` from `_compute_deformed_geometry`. Validated - under
      that shadow `tests/validation/rotor/test_bem_fsi_deformed_geometry.py` goes
      RED (the loads stop moving with the deformation), which is exactly the
      mutation. Run it with
      `sbatch -p sequana_cpu --time=03:30:00 -J frame_ab --export=ALL,LABEL=nofeed,SHADOW=$SCRATCH/shadow-nofeed,RUN_DIR=$SCRATCH/frame_ab/nofeed,FLUID_YAML=<repo>/tests/smoke_fix/frame_ab/fluid_yaw_0.yaml,SOLID_YAML=<repo>/tests/smoke_fix/frame_ab/solid_corotational_yaw_0.yaml,INNER_TIMEOUT=10800 tests/run_coupling_probe.srm`.
      If `nofeed` contracts, the non-contraction is the *strength* of the
      deformation feedback interacting with the coupling, not a geometry bug;
      if it still saturates, the defect is upstream of the feedback.
- [x] T6j **Fixed and verified.** Commit `6765633` keeps the reference projector in
      `_compute_forces` (the BEM geometry still follows the deformation) and removes
      the now-unused `_rebuild_projector`. Job `11611246` (`fixfull`, the repo tree
      at that commit, campaign's own case, `max-time 5.0`) reads **90 of 90 windows
      (100%)**, first 12 = `11 7 8 11 10 8 7 9 9 8 7 8` - the campaign's own healthy
      sequence. Guard: `tests/validation/rotor/test_bem_fsi_deformed_geometry.py::
      test_force_projection_geometry_stays_on_the_reference` (commit `e3278ba`) -
      for the same deformed BEM result the nodal forces must be the reference
      projector's. RED against `7a84da1`: all 9129 elements differ, max absolute
      difference 1372.75 N; green after. Suite: the participant + multicell guards
      23 passed, the de-loading table and the rated-twist fixture 15 passed, ruff
      clean.
- [ ] T6k **Known trade-off, for the reviewer.** The within-strip load distribution
      no longer sees the deformed offsets. It is small (the offsets are chord-scale
      while the deformation is blade-scale) but it is a fidelity loss, and the
      campaign's *own* projector followed the deformation and contracted - so the
      smoother alternative stays open: make the `#11` geometry derivation
      deformation-robust (hysteresis-free ring selection, LE/TE tie-break and
      chord-axis continuity) instead of frozen. That is the fidelity follow-up, not
      part of this closure.
- [ ] T7 Re-run the 30 s gate (`tests/run_step1b_smoke.srm`) and re-anchor the
      campaign baseline now that the coupled path contracts.
      Resume (2026-10-11): #33 L5 (`fee2cc8`) changed `r_def`/`Rtip`; the 5 s case still
      reads 500/500 (job `11613014`) but is not stationary, so the steady thrust needs
      this gate. The script always writes `$SCRATCH/smoke_fix_results/base_fix/`: copy it
      to `base_fix_pre_radii` first (job `11611263`, the before side), then
      `sbatch --time=06:00:00 tests/run_step1b_smoke.srm` (~2.5 h; do not edit `src/`
      while it runs). Check stationarity before averaging, then thrust vs the rigid
      2.541 MN and vs `base_fix_pre_radii`.
      **Result of the before side (2026-10-11): T7 is RED, it does not close.** Job
      `11611263` (tree of 2026-10-08, after `6765633`/`7a84da1`, before `fee2cc8`), kept
      in `$SCRATCH/smoke_fix_results/base_fix_pre_radii/` with its `.out`/`.err` in
      `logs/`. The coupling contracts (3000/3000 windows, mean 2.14 sub-iterations, none
      at 30), but the response is a full-amplitude **1P limit cycle**: peak spacing
      7.95 s = 2π/Ω, peak at a fixed azimuth (Ωt mod 2π = 306°), and every whole period
      is identical, so it is periodic-stationary, not drifting. Whole-period statistics
      (T = 7.947 s, last 3 periods): thrust **2.726 MN (1.92-3.52)**, tip **17.99 m
      (11.4-24.3)**, power 13.89 MW, dropping to 8.7 MW at the thrust peak. The script's
      "steady 2.680 MN" averages [15, 30] = 1.89 periods and is biased; the script's
      own criterion (flexible thrust below the rigid 2.541 MN) fails either way.
      The same operating point on the 2026-10-01 tree
      (`$SCRATCH/frontiersin_results_corotational_100s/yaw_0`, case
      `tests/IEA15MW/frontiersin_2025_yaw/`, which differs from `smoke_fix/base_fix` only in
      `velocity_data` and `airfoil_spacing`) reads 1.640 MN (1.60-1.66), tip 7.44 m
      (7.0-7.9), 12.35 MW: a small 1P. Neither matches G4's recorded 2.20 MN /
      12.7-12.8 m / 14.71 MW, so G4's source run is a third, unidentified state.
      **The blind spot:** every #19/#26 verification (500/500, `11613014`) ran at
      `max-time 5.0`, shorter than one rotor period (7.95 s), so it could not see a 1P
      mechanism. Contraction says the coupling converges, not that the load is right.
      Mechanism candidates, none tested yet: the solid is **one blade** and the BEM
      scales by `n_blades: 3`, so a per-blade 1P shows up in "rotor" thrust; gravity is
      the only 1P forcing in the case (uniform inflow, no shear/tilt/yaw) but is mostly
      edgewise and implausible for +/-0.8 MN and +/-6 m of flap; prime suspect is
      `6765633` (projector frozen on the reference geometry). Check its frame against
      the azimuth: a projector frozen at one azimuth and applied to a rotating geometry
      modulates the load at exactly 1P.
      After side: job `11613031` at `cecf020` (includes `fee2cc8`), writes `base_fix/`; it is a second sample
      of the cycle amplitude, **not** the T7 closure. T8/T9 and the campaign re-anchor
      are blocked until the 1P is explained.
- [ ] T7a The 1P is a frame mismatch at the BEM participant's preCICE boundary
      (found 2026-10-11, branch `fix/bem-participant-rotor-frame`). `rotor_fsi.rs`
      writes `disps_to_inertial(u, theta_target)` = R(θ)·u and reads the forces as
      global (`forces_to_rotating`, Rᵀ). `BEMFSIParticipant` used the read
      displacements as rotating-frame and wrote rotating-frame forces, so the BEM saw
      the in-plane and span components (`r_def`, ring twist) mixed by sin θ/cos θ.
      Evidence on `base_fix_pre_radii`, last 3 periods: raw tip z swings −2.2..1.5 m,
      while z_rot = sin θ·x + cos θ·z has corr(z_rot, flap y) = 1.0000 (raw: 0.717) and
      z_rot/y = 0.0650 ± 0.0007, the ratio at t = 0.01 s. The Rust convention dates from
      2026-04-28 (`54b25ee`), so this is not a post-10-01 regression: the amplitude grew
      once #11/#33 began reading exactly the components the rotation scrambles.
      Route: delegated writer (2+ non-trivial files: participant + test + case key).
      Job `11613031` (after side) could not be cancelled (no permission), so it runs to
      completion on the corrupted state; the work happens in a worktree.
  - [x] T7a.1 RED/GREEN unit test: round trip against Rust's Rodrigues sense.
        Evidence: `tests/validation/bem/test_bem_participant_rotor_frame.py` RED 53 failed
        (AttributeError, helpers absent) before the change, 53 passed after.
  - [x] T7a.2 Transform at the IO boundary (end-of-window θ, explicit axis and
        partner-frame keys); the smoke case declares them.
        Evidence: `bem.partner_frame` (default `rotating`) + `bem.rotation_axis` in
        `fsi_participant.py`; frame angle mirrors Rust θ from 0 (not `bem.azimuth`);
        `base_fix/fluid_yaw_0.yaml` declares `inertial`/[0,1,0]. Existing participant
        tests (deformed geometry, wall flow, rated twist): 23 passed; ruff clean.
  - [ ] T7a.3 30 s gate on the fix: the 1P is gone (≥ 2 periods), then the mean
        against 2.541 MN, reported separately.
        Ran: job `11613040` (submitted 2026-10-11 from this worktree,
        `RUN_NAME=frame_fix_30s`), output in
        `runs/smoke_fix/frame_fix_30s/fluid/bem_report.csv` (worktree-relative); the
        pre-radii before side `11613031` was cancelled on purpose. The transform is live:
        over the first 3 s the tip reads z/y = 0.065 constant while θ sweeps 23°→136°
        (raw columns would swing). Early warning, not a closure: thrust 2.67-2.69 MN and
        flap tip ≈ 18 m at t ≤ 3 s, still above the rigid 2.541 MN and far from the
        2026-10-01 tree's 7.4 m - if the 1P is gone but the mean stays there, that is a
        separate finding.
        **Result (2026-10-12, job `11613040` COMPLETED, 30 s):** the limit cycle is
        gone, a clean 1P is not. Least-squares fit (mean + trend + 1P, T = 7.947 s) over
        [10, 30] s, after vs `base_fix_pre_radii`:
        | quantity | frame_fix_30s | base_fix_pre_radii |
        |---|---|---|
        | thrust 1P amplitude | 42.4 kN (1.57 %), resid std 303 N | 853 kN (31.4 %), resid std 64 kN |
        | flap tip Y 1P amplitude | 0.228 m (1.25 %) | 7.04 m (39.4 %) |
        | edge tip X 1P amplitude | 1.24 m | 1.54 m |
        | trend | −1.9 N/s (steady) | −319 N/s |
        The residual is deterministic (resid std ≪ amplitude), so the criterion "1P gone"
        is **not met**; no tolerance is set after the fact. Phases: tip X 0°, thrust
        and Mb −178°, torque +4°, i.e. locked to the gravity edge cycle (case: uniform
        inflow, shear/tilt/yaw/precone 0, no tower key; `gravity:` on in the solid).
        Thrust is reported as B × one replicated blade (`BladeAero(n_blades=…)` on a
        single deformed blade, `fsi_participant.py:999`), so a per-blade gravity 1P
        reaches the rotor total instead of cancelling to 3P; read the residual per
        blade (1.57 %, Mb ±1.29 MN·m), not as a rotor load.
        Confound: the after side carries `fee2cc8` (radii) and the frame fix; the
        one-variable side `11613031` (radii, no frame fix) stopped at 8.33 s, under 2
        periods, so the 20× is attributed through the corr = 1.0000 diagnosis above,
        not through a single-variable pair.
        **Mean, a separate finding (not chased here):** thrust 2.709 MN, tip Y 18.15 m,
        essentially unchanged from `base_fix_pre_radii` (2.726 MN / 17.99 m); +6.6 %
        above rigid 2.541 MN despite an 18 m flap (sign-check first, discipline #4),
        and far from 10-01 (1.640 MN / 7.44 m). The frame fix did not move the mean.
  - [x] T7a.3b Discriminate the residual 1P: gravity-off run of the same case (zero
        `gravity:` in `solid_corotational_yaw_0.yaml`, ~20 s ≈ 1.7 h SLURM). Any 1P
        that survives is an artifact; prime suspect the frozen projector `6765633`
        (feeds T7a.4). Then T7a.4, then add `partner_frame: inertial` to the other ~60
        fluid YAMLs paired with `rotor_fsi` (the `solid_inertial_*` pairs need their
        own frame analysis).
        In flight: job `11613217` (submitted 2026-10-12 at `ae5cc26`, 30 s, same command
        as `11613040` with `SOLID_YAML=solid_corotational_yaw_0_nograv.yaml`,
        `RUN_NAME=nograv_30s`). **Resume recipe:** `sacct -j 11613217`; run the 1P fit
        (mean + trend + 1P, T = 7.947 s, [10, 30] s) on
        `runs/smoke_fix/nograv_30s/fluid/bem_report.csv`; 1P amplitude ≈ 0 (≪ 42 kN,
        resid-std level) ⇒ the residual is gravity; a surviving 1P ⇒ artifact, go to
        T7a.4. The mean without gravity is also a new datum for the +6.6 % finding.
        **Result (2026-10-12, job `11613217` COMPLETED, 30 s): the residual 1P is
        gravity.** Same fit (mean + trend + 1P, T = 7.947 s, [10, 30] s) on
        `runs/smoke_fix/nograv_30s/fluid/bem_report.csv` vs `frame_fix_30s`:
        | quantity | nograv_30s | frame_fix_30s |
        |---|---|---|
        | thrust 1P amplitude | 3.9e-10 N (resid std 4.7e-10 N) | 42.4 kN (1.57 %), resid std 303 N |
        | Mb 1P amplitude | 6.8e-9 N·m | 1.29 MN·m (1.76 %) |
        | edge tip X 1P amplitude | 5.7e-10 m | 1.24 m |
        | flap tip Y 1P amplitude | 7.4e-15 m | 0.228 m |
        | mean thrust / tip Y / power | 2.7086 MN / 18.142 m / 15.948 MW | 2.7082 MN / 18.141 m / 15.942 MW |
        Without gravity every channel is constant to round-off over [10, 30] s: the
        coupled solution is a true fixed point, so no frame, projector or replication
        artifact survives at 1P. The residual in `frame_fix_30s` is the per-blade
        gravity response, which reaches the "rotor" thrust only because the BEM
        replicates one deformed blade (`n_blades: 3`); three phased blades would
        cancel it to 3P. T7a.4 is no longer a 1P suspect; it stays as the fidelity
        follow-up of T6k. Gravity does not move the mean (Δ 0.4 kN), so the +6.6 %
        above rigid 2.541 MN is not gravity.
  - [x] T7a.4 Retest the deformation-following projector: the rationale for `6765633`
        was measured on rotated displacements.
        **Decision (2026-10-13): keep the frozen projector.** Grounds: cost, not
        contraction. The following projector rebuilds at ~63 s per sub-iteration on the
        27.6 k-node smoke mesh (~40× the whole frozen sub-iteration), so a 30 s gate would
        take ~114 h; it cannot be the production path whatever Gate A says. The job
        `11614191` was cancelled by the user before window 0 converged (sub-iteration 8 at
        05:37:56; the baseline's window 1 took 11), so **Gate A is untested: whether the
        old 39/500 non-contraction was the frame bug stays open**, and Gate B was never
        reachable. Not a null result, an unanswered one.
        What survives as a finding is the static conservation bias of the frozen
        projector (+1.23 % of ∫Np dr at max |u| 17.6 m, following +0.15 %; numbers
        below). It is a candidate for a cheap, targeted change (strip widths from the
        deformed radii on the frozen projector, no rebuild), to be tested on its own
        (follow-up T7a.5), not folded into this decision. The comparison worktree
        `../fem-shell-worktrees/t7a4-follow` (branch `exp/t7a4-follow-projector`,
        `b0f01e4`) is throwaway; keep it until T7a.5 is scoped, then remove it.
  - [x] T7a.5 Does the frozen projector conserve the BEM load
        on the coupled run? Metric: ΣF_normal vs ∫Np dr per window (and `bem.thrust/B`);
        reference: the BEM integral itself; static bias to explain: +1.23 % at 17.6 m.
        First prove the mechanism (strip widths vs shortened radii, discipline #4), RED
        test on the `test_bem_fsi_deformed_geometry.py` fixture, then the fix and the 30 s
        gate against `t7c_30s` with `tools/coupled_run_stats.py`. Run pytest in a worktree
        with `PYTHONPATH=$PWD/src` (incident below).
        **Setup (2026-10-13).** Comparison tree: worktree `../fem-shell-worktrees/t7a4-follow`,
        branch `exp/t7a4-follow-projector` = HEAD `f418fb4` + `git revert 6765633` + one
        patch so the rebuilt projector receives exactly the reference's constructor
        arguments (`tangential_direction=self._configured_tangential`,
        `element_properties`). The plain revert was **not** one variable: it passed a
        hard-coded `[1,0,0]` and dropped `element_properties` (`81a1e15`, `d42ddda` came
        later). Only the coordinates and the strip radii follow the deformation. Frame:
        `_compute_forces` receives rotating-frame `u` (`_displacements_from_partner`,
        `fsi_participant.py:667`); the case reads `partner_frame: inertial`.
        Wall-flow realisation is off in the smoke case (`fluid.log:10`), so both sides use
        the minimum-norm strip distribution.
        Liveness (fixture of `test_bem_fsi_deformed_geometry.py`, max |u| 17.6 m): the guard
        `test_force_projection_geometry_stays_on_the_reference` goes RED in the comparison
        tree (max nodal diff 1383 N; `e3278ba` read 1372.75 N), the other 7 pass.
        **Incident:** pytest in a worktree without `PYTHONPATH=$PWD/src` imports the
        *main checkout's* `aeroelast` (editable install; at `chore/cite-by-reference`), so
        the guard passed and the `-x` refusal failed on the wrong tree. With the variable
        set, HEAD reads 7/7 and `test_bem_participant_rotor_frame.py` 53/53. The SLURM
        srm sets it (`a48ca20`), so the coupled runs were not affected.
        **Static finding (same fixture, same BEM result), the frozen projector does not
        conserve the load:** ΣF_normal reference 624.31 kN (**+1.23 %** of ∫Np dr =
        616.75 kN), following 617.66 kN (+0.15 %), `bem.thrust/B` 614.84 kN. Hypothesis,
        not proven: the reference strip widths carry Np evaluated on the shortened
        deformed radii (`fee2cc8`). Cost: projector rebuild 2.68 s on the 9 k-node fixture
        (BEM rebuild 0.13 s); the smoke fluid mesh has 27 609 nodes.
        **Pre-registered (before launch).** Same srm, same case, only the tree differs;
        before side `t7c_30s` (job `11613577`, `src/` identical to HEAD).
        Metrics: `tools/coupled_run_stats.py` on both run dirs, with `--t-end` = the
        common end time (truncating the deterministic frozen run is exact).
        - Gate A, the rationale of `6765633`: every window converges, 0 windows at the
          30-iteration limit. Reference `t7c_30s`: 3000/3000, 0 at limit, mean 2.11,
          max 11. Fail ⇒ the rationale still holds ⇒ **keep frozen**.
        - Gate B (only if A passes): means over the last 2 whole periods
          (T = 7.947 s) of thrust, flap tip Y, power, Mb. Reference `t7c_30s`
          (14.106, 30]: 1.9746 MN / 10.3366 m / 14.3384 MW / 49.1957 MN·m (identical over
          (7.946, 23.840]: stationary). Decision threshold 1 % on thrust and flap,
          justified as: well under the 9-point gap to Zhou (-22.3 % vs -13 %) that
          motivates the fidelity work, and above the 0.4 % flap spread between two
          physically identical runs (`tplus_30s`/`t7c_30s`). |Δ| < 1 % on both ⇒ null:
          **keep frozen** (no material change at 20× the projector cost), record the
          static conservation bias as bounded by the measured Δ. |Δ| ≥ 1 % ⇒ the choice
          is material: **change**, toward whichever variant conserves ∫Np dr (statically:
          the following one), RED conservation test first, then the 30 s gate.
        Price: ~6300 sub-iterations × ~8 s rebuild (27 k nodes, linear scaling
        assumed) + 2.75 h base ≈ 17 h on 8 cores; submitted with `--time=24:00:00`.
        Gate A is readable from the iterations log within the first hours.
        **Launched:** job `11614191` (2026-10-13 05:28, from the comparison worktree,
        `RUN_NAME=t7a4_follow_30s`, output `../t7a4-follow/runs/smoke_fix/t7a4_follow_30s/`).
        **Price miss, recorded:** the measured cost is ~65 s per sub-iteration (fluid log,
        window 0, iterations 0-7), against ~1.6 s for the frozen baseline (2 h 45 min /
        6330). The rebuild scales far worse than linearly in nodes (2.7 s at 9 k, ~63 s at
        27.6 k). Within the 24 h cap the job reaches ~1300 sub-iterations, about the
        rationale's own 5 s / 500 windows, so **Gate A is testable and Gate B is not**
        (30 s ≈ 114 h). Cost finding regardless of A: a full rebuild per sub-iteration is
        ~40× the whole frozen sub-iteration, so "change to the following projector" is not
        viable as-is for production runs even if it contracts.
        **Mechanism proven (2026-10-13, step 1-2 of T7a.5).** Code path: `project()`
        computes `F_n = Np[k]*strip.dr` (`force_projection.py:711-787`), pairing reference
        strip k with deformed BEM station k by index, no interpolation; the widths are the
        reference `s.dr` (`fsi_participant.py:497`). The rotor-normal axis is the same in
        both projectors and the AC arm / `Mp` enter only moments, so neither can move
        ΣF_normal. `ForceProjector.verify` checks against the projector's own `dr`, so it
        cannot see this. Metric correction: `bem.thrust/B` is not the conservation target
        (CCBlade integrates it azimuth-averaged with its own end rule, 0.31 % below
        trapz); the target is ∫Np dr on the deformed span.
        Diagnostic (`/tmp/t7a5_diag.py`, throwaway; same fixture, one `bem_result`,
        `aeroelast` from this worktree, max |u| 17.63 m, 50 strips), kN:
        | quantity | deformed | rigid |
        |---|---|---|
        | (a) ΣF_normal, frozen nodes | 624.3 | 851.4 |
        | (b) Σ Np·dr_ref | 624.3 ((a)−(b) = −1e-10 N) | 851.4 |
        | (c) Σ Np·dr_def (same strip rule on r_def) = rebuilt projector | 617.7 | — |
        | (d) trapz(Np, r_def) | 616.8 | 847.5 (on r_ref) |
        | `bem.thrust/B` | 614.8 | 841.7 |
        Split of (a)−(d) as % of (d): distribution loss 0.0000 + width term D **+1.079**
        + end half-cells E **+0.147** = +1.225. D is the reference widths on radii pulled
        inward by the second-order shortening (tip r_ref−r_def = +1.093 m; largest width
        loss 0.102 m at k=47, 4.3 % of the strip). E is the mirrored end half-cells,
        Np₀Δr₀/2 + Np_NΔr_N/2 (matched to 1e-10 N); it is not a constant offset: rigid
        +0.467 % (tip Np 3203 N/m) vs deformed +0.147 % (787 N/m). So deformed widths
        remove D only; E is a load-dependent quadrature bias already present in the rigid
        path. Side finding, unverified: the deformed root station sits 1.6e-5 m inboard of
        r_ref[0], CCBlade reports NaN at station 0 and zeroes it (deformed Np₀ = 0 vs
        115 N/m rigid); suspected slip inside Rhub, effect on E ≈ 0.016 %.
        **Scope decision (user, 2026-10-13): fix D only.** `project()` takes the strip
        widths from `bem_result.r` with the constructor's strip rule; node assignment, AC
        arm, axes and distribution stay frozen (no rebuild). Rigid path is a no-op. RED
        target: ΣF_normal = Σ Np·dr_rule(r_def) at rel 1e-9 (what the rebuilt projector
        delivers). E (end half-cells, also in the rigid path, moves every BEM-FSI anchor)
        is deferred to its own follow-up, T7a.6, together with the root-station NaN.
        **Fix: `6e5d6a2`.** RED `test_frozen_projection_conserves_the_deformed_strip_load`
        624.3097 vs 617.6559 kN (+1.077 %) → GREEN; deformed-geometry 8/8, rotor-frame
        53/53, force projection 15/15, tube parity 1/1, `check.sh quick` OK (only the
        known-red). Rigid path: widths vs `strip.dr` 7e-15 m, ΣF 2.2e-16 rel. Left on
        reference widths, not load integration: VTU `Force_section`/`dr_m`
        (`_write_sections_vtu`, diagnostic, ~1 % off when deformed), `_strip_widths` and
        the LE weighting in `_resolve_twist_sense` (geometry weights).
        **Pre-registered 30 s gate (before launch).** One variable: `git diff f418fb4
        6e5d6a2~1 -- src tests/smoke_fix tests/run_step1b_smoke.srm` is empty, so the
        before side is `t7c_30s` (job `11613577`). Same srm, `RUN_NAME=t7a5_widths_30s`.
        Metrics: `tools/coupled_run_stats.py` on both run dirs, common `--t-end`.
        - Gate A (the 6765633 risk: a projector-geometry change can break contraction):
          every window converges, 0 at the 30-iteration limit; reference 3000/3000,
          mean 2.11, max 11. Fail ⇒ finding, the fix does not ship as-is.
        - Gate B (record, not decide): Δ of the last-2-period means (T = 7.947 s) of
          thrust, flap tip Y, power, Mb vs 1.9746 MN / 10.3366 m / 14.3384 MW /
          49.1957 MN·m. Prediction: the width term scales with the shortening ∝ u², so at
          the coupled flap (10.3 m vs the fixture's 17.6 m) it is ≈ 1.08 % × (10.3/17.6)²
          ≈ 0.37 % of the applied normal load; expected flap Δ ≈ −0.4 %, Mb ≈ −0.4 %,
          BEM thrust slightly up (less de-loading). This is at the 0.4 % flap spread
          between `tplus_30s`/`t7c_30s`, so a null is the expected outcome. The fix is kept
          on conservation grounds at zero cost whatever B reads; B only bounds how much
          the bias mattered to the coupled solution. A |Δ flap| > 1 % would contradict the
          u² scaling and is itself a finding.
        **Launched:** job `11614215` (2026-10-13, from this worktree at `f939653`,
        `RUN_NAME=t7a5_widths_30s`, ~2 h 45 min). **Resume recipe:** `sacct -j 11614215`;
        then `PYTHONPATH=$PWD/src scripts/aeroenv.sh python tools/coupled_run_stats.py
        runs/smoke_fix/t7c_30s runs/smoke_fix/t7a5_widths_30s --t-end 30` and read Gates
        A/B above.
        **Result (2026-10-13, job `11614215` COMPLETED in 2 h 44 min, read with
        `tools/coupled_run_stats.py`, common `--t-end 30`, window (14.106, 30] s).**
        Gate A **passes**: 3000/3000 windows converge, 0 at the 30-iteration limit,
        mean 2.17 iterations (reference 2.11), max 11 (reference 11).
        Gate B (record): thrust 1.9752 MN (+0.03 %), power 14.342 MW (+0.03 %), Mb
        49.219 MN m (+0.05 %), flap tip Y 10.290 m (-0.45 %), edge tip X 1.457 m (+0.3 %).
        Against the prediction: the flap Δ (-0.45 % vs expected -0.4 %) matches; Mb has
        the opposite sign to the expected -0.4 % but moves 0.05 %, inside the 0.4 % spread
        between `tplus_30s` and `t7c_30s`, so it neither confirms nor contradicts the u²
        scaling. |Δ flap| is far below the 1 % that would have been a finding. Verdict:
        the fix stays on conservation grounds; the width bias barely moved the coupled
        solution. Whole-period means were not separated from the last-window means here.
  - Known gap: Rust `with_initial_state` restart θ is not mirrored on the fluid side.
- [x] T7b The +6.6 % mean thrust is a reversed tangential force: every production
      fluid YAML (65 of 65 carrying the key) sets `tangential_direction: -x`, which
      pushes the edgewise load toward the trailing edge.
      Evidence (2026-10-12, all local, no HPC):
      - The coupled structure twists **nose-up**: BEM tip elastic twist −3.18°, AoA
        +2.3°, tip Np +28 % (`nograv_30s`, t = 30 s vs 0.01 s); the solid's own
        translation fit at z = 80 m reads −1.22° against the participant's −1.25°, so
        the participant transfers it faithfully (factor +1). Zhou is +3.6° nose-down,
        BeamDyn +2.08°, #14's static shell +5.78°.
      - Linear split on the fixture mesh (K + K_G + K_SP, rigid rated aero): aero alone
        +8.25° nose-down at the tip, centrifugal alone −0.13°. The centrifugal / K_G /
        K_SP path is not the flip. The VTU's `F_INERTIAL` is written in the inertial
        frame (R_y(Ωt)·f); rotated back it is the constant (−18 kN, 0, +1.209 MN).
        VTU `points` are deformed coordinates.
      - Same mesh, polars and K, only `tangential_direction` changed: +7.78° → −1.87° at
        r = 118.6 m. Both fields apply nose-down torque about the section centroid
        (outboard of 80 m: +524 vs +256 kN·m); the sign is set by the edgewise force
        (outboard Fx +34 vs −34 kN) through the shell's bend-twist coupling.
      - Which sense is physical: the projector measures the leading edge at +x
        (`-moment_axis_sign · chord_dir` = (0.999, 0.032, 0)) and the solid's ω × r at
        the tip is (+92.6, 0, 0.25) m/s. A blade moves toward its leading edge, and
        positive Tp (the BEM's positive torque and power) acts along that motion, so
        the force must point +x. afac31c's docstring ("sets it to -x so that Tp drives
        the rotor at positive torque") has it backwards: z × (−x) = −y opposes Ω = +y.
        The code's own original design agrees: `_resolve_twist_sense`'s fallback
        `normal · (span × tangential)` is +1 only for tangential +x (z × x = +y). The
        defaults are +x too (`core/config.py:467`, `fsi_participant.py:429`), so a YAML
        without the key is already physical and the 65 are the complete reversed set.
      - Static FSI fixed point, `tools/diagnose_tangential_sense.py` (14 iterations,
        relaxation 0.5, rigid 2.554 MN): −x → **+5.55 %**, flap 16.6 m, tip twist
        −0.86°; +x → **−21.6 %**, flap 10.0 m, +6.28°. The −x row reproduces the
        coupled anomaly (+6.6 %, 18.1 m). The +x overshoot against Zhou (−13 %,
        13.9 m) is a magnitude gap (here +6.28° vs +3.6° ≈ 1.7×). **Correction
        2026-10-13:** the "#14 shell over-twist, 1.61×" label is not a validated cause. The 1.61× is
        the figure of a diagnostic with foreign loads, not the production rated case
        (`docs/validation_closures.md`, S-7 block), and the shell is stiffer than the BeamDyn deck
        (#20). See T7d/T7e and `shell-twist-mixed-mesh.md`. Separate finding, not tuned.
        Both rows converged to 4 digits over the last 3 iterations.
      **Pre-registered expectation for `11613528`** (written before the result, from the
      static +x row shifted by the −x pair's static→coupled offset, +5.55 → +6.6 %):
      mean thrust ≈ −20 % of rigid (≈ 2.03 MN), flap tip ≈ 10-11 m, BEM tip elastic
      twist ≈ +6° (nose-down), coupling still contracting (≈ 2 sub-iterations per
      window), gravity 1P ≈ 1.5 % per blade. A large miss on any of these is a finding,
      not a reason to move the expectation.
      Blast radius if confirmed, to list in T8, not edited here: the 65 fluid YAMLs; the
      tests and tools that encode −x as the production frame
      (`tests/validation/rotor/test_bem_fsi_deformed_geometry.py`,
      `tests/test_multicell_shear_flow.py`, `tests/test_rotor_mesh_export.py`,
      `tests/validation/rotor/test_rotor_inertial.py`, `tools/diagnose_fsi_deformed.py`,
      `tools/diagnose_twist_sense_frame.py`, afac31c's docstring); every coupled or
      BEM-loaded anchor produced with −x (S-5, S-6, the V-0x and campaign rows) is
      expected to move - report the move, never re-band.
      **Orientation confirmed with the maintainer (2026-10-12).** Maintainer: "la
      orientacion correcta es la de Aeroelas estructural de winio" - wind +Y, rotation
      axis +Y (right hand), blade along +Z with prebend toward −Y. Measured on the
      aeroelast mesh (`solid_mesh.vtu`): the leading edge is +x at every station
      (blunt end, nearer the pitch axis), so ω × r = +x moves the blade toward its
      leading edge. The `turbinesFoam/geometry/stl/iea15mw_blade.stl` is the x-mirror
      of this mesh (leading edge −x); the maintainer handles that side.
      **Decision: the BEM follows the structure.** Fix plan (T7c, after `11613528`):
      1. RED tests first: a configured `tangential_direction` that contradicts the
         measured leading edge is refused; F_t points toward the measured leading
         edge; the aerodynamic torque on the structure about +y is positive.
      2. `ForceProjector` derives the Tp sense from the measured leading edge (the
         afac31c pattern) instead of the config; a contradicting config raises with a
         message naming the fix.
      3. Drop the key from the 65 fluid YAMLs; fix afac31c's docstring and the tests
         and tools listed above that encode −x.
      4. Independent verifier (high-risk item 3), then the 30 s gate on the fixed tree.
      In flight: job `11613528` at `d233e38`, `RUN_NAME=tplus_30s`,
      `FLUID_YAML=fluid_yaw_0_tplus.yaml` (only `tangential_direction` differs from the
      base case; gravity on), the one-variable pair of `frame_fix_30s`.
      **Result of `11613528` (2026-10-12, COMPLETED 02:18:45, `tplus_30s`, window
      [10, 30] s, 2001 samples):** mean thrust **1.976 MN = -22.3 %** of rigid 2.541
      (expected ~ -20 %, hit); flap tip mean **10.31 m** (expected 10-11, hit; the
      -x pair is 18.1 m); trend +0.36 kN/s (flat); 1P amplitude **2.7 %** of mean
      (expected ~1.5 % per blade; the single-blade figure is not separated here, so
      this row neither confirms nor refutes it); 2 sub-iterations per window
      (iteration 0 and 1, dF_iter ~ 22 N vs |F| ~ 7 kN: contracting, hit). BEM tip
      elastic twist was **not measured** (the log does not carry it); the sign is
      inferred from the thrust and flap moving together with the static +x row.
      The reversed tangential force is the cause; T7c proceeds as decided.
      **Resume recipe:** `sacct -j 11613528`; run the mean + trend + 1P fit (T = 7.947 s,
      [10, 30] s) on `runs/smoke_fix/tplus_30s/fluid/bem_report.csv`; expect thrust
      below the rigid 2.541 MN, flap tip well below 18 m and a positive BEM tip
      elastic twist. If it holds, the decision left to the maintainer is the fix's
      form: flip the 65 YAMLs, or have the projector derive the Tp sense from the
      measured leading edge (as afac31c did for the twist) and refuse a configured
      sense that contradicts it. Either way it is a contract change (high-risk item
      3) and invalidates every coupled result produced with `-x`.
- [ ] T7c The projector derives the Tp sense from the measured leading edge
      (`d42ddda`, `0c96b41`, `685b3c7`, `1877b76`; 65 YAMLs lose `-x`, a
      contradicting `tangential_direction` raises `ValueError`).
  - [x] T7c.1-3 RED/GREEN, tests and tools, YAMLs. RED: only the refusal test failed
        (`DID NOT RAISE`); the default was already +x.
  - [x] T7c.4a Independent verifier (2026-10-12): PASS. 161 passed, 12 skipped
        (OpenFAST absent), 1 xfailed (NeuralFoil 5.9 % vs 5 %, not new); no bound
        widened; 685b3c7 is deletions only. Open observations: (A) the in-plane
        leading-edge guard is `norm <= 1e-12`, so a near-cancelling edge would set the
        sense from noise (unproven on real blades); (B) `test_multicell_shear_flow.py:1150`
        still passes `[-1,0,0]` and relies on the guard not firing on its ring mesh.
  - [x] T7c.4b 30 s gate on the fixed tree: job `11613577`, `RUN_NAME=t7c_30s`,
        COMPLETED 02:45:01 (2026-10-13). Window [10, 30] s, fit T = 7.947 s, same
        script on both sides:
        | quantity | tplus_30s (YAML +x) | t7c_30s (derived sense) |
        |---|---|---|
        | mean thrust | 1.976 MN (-22.3 % of rigid 2.541) | 1.976 MN |
        | flap tip Y mean | 10.31 m | 10.35 m |
        | thrust 1P / mean | 2.69 % | 2.69 % |
        | power mean / Mb mean | 14.34 MW / 49.23 MN m | 14.34 MW / 49.23 MN m |
        The projector-derived sense reproduces the one-variable +x run to 4 digits
        (flap +0.4 %), so the fix is the same physics without the YAML key. Not
        measured: BEM tip elastic twist (the log does not carry it).
  - [x] T7d BEM tip elastic twist on the fixed tree (2026-10-13, one measurement session,
        no tuning; precondition T7a.4 decided, Engram obs 357). Metric: last-strip
        `twist_deg` in `fluid/<t>/bem_sectional.csv` minus the production reference
        `load_blade_aero(tests/IEA-15-240-RWT.yaml).twist`, whole periods
        (T = 7.947 s, [14.11, 30] s); positive = nose-down, Zhou's sense. Instrument
        checks: root strip +0.000 deg, tip reference twist -1.242 deg (official sheet).
        Tolerance, written after the first read and derived from the references:
        consistent with Zhou iff |M - 3.6| / 3.6 <= 0.42 (Zhou to BeamDyn +2.08 deg).
        **Result:** `t7c_30s` and `tplus_30s` identical to 4 digits, tip **+5.589 deg**
        (min +5.243, max +5.934); profile peak +6.09 deg at r = 113.7 m. Zhou 1.55x
        (+55 %, outside the band); static +x fixture +6.28 deg 0.89x; BeamDyn 2.7x.
        Thrust per degree of twist: ours 4.0 %/deg vs Zhou 3.6 %/deg (weak check).
        **Outcome (b): no model defect with a named cause.** The "stiffness over-twist"
        reading is not supported by earlier measurements (shell stiffer than BeamDyn by
        6.4 % global / 23.9 % sectional, #20; independent-code response 4.22 %; Zhou
        comparison non-transferable because of section distortion). **Superseded by T7e:**
        the load side is ruled out as the main source (same load, shell 2.82x the beam); the
        open candidates are in the shell's own response, see `shell-twist-mixed-mesh.md`. Recorded in `docs/validation_closures.md` as an open magnitude gap
        that does not block re-anchoring. Proposed, not run: Zhou's Fig. 11 load on the
        coupled shell to separate load from structure. Script: not committed (scratch).
  - [x] T7e Zhou's published loads on our shell and beam (2026-10-13, existing
        `tools/diagnose_zhou_loads_reverse.py` run unmodified; needed a git-ignored
        `.sources` symlink into the worktree). Same load, tip twist nose-down positive:
        shell +5.855 deg (1.63x Zhou; 0.25 m, the 0.5 m figure +5.955 is not citable), standard beam +2.079 (0.58x), Zhou +3.60. The coupled
        shell under BEM loads gives +5.59, so the load source is not what drives the 1.55x;
        T7d's "load-side candidates" is withdrawn. Open: shell vs beam differ 2.82x under
        identical load (section distortion the beam cannot represent, vs an application or
        mesh effect in the shell). Not a named defect; still an open magnitude gap that does
        not block re-anchoring. Next, not run: a closed section with an analytic torsion
        solution, shell vs beam.
  - Leftovers for T8: `docs/model_parity_audit.md:51` (stale MATCH claim),
    `tools/diagnose_tangential_sense.py:84`, `tools/diagnose_twist_sense_frame.py:11`.
- [x] T8 Reconcile the store and the closure log. Finding: both ids named here were
      already gone (`1b5256c`, 2026-10-09: `force_projection_sense_p5` removed,
      `rotor_coupling_noncontraction` replaced by
      `reference_projection_geometry_approximation`), and groups 27/35 carry no
      `tangential_direction` claim, so the store needed no edit
      (`validation_matrix.py check`: 211 rows, 286 comparisons, 0 errors). Done: new
      closure section in `docs/validation_closures.md`, `docs/model_parity_audit.md`
      A23 corrected, `tools/diagnose_tangential_sense.py` no longer builds a `-x`
      participant (it would raise) and reports the refusal. `diagnose_twist_sense_frame.py:11`
      already says "used to". Not run: the diagnose script itself (needs the full env).
- [x] T9 Comment and close #19 (done 2026-10-13; #26 was already closed), and fold #26's remaining closure criteria into
      it.

### Corrections to the issue's own text

- The window is not "the 2026-10-04 merge plus the #11 batch": `9a3923e` is
  2026-10-02, before the merge. The merge merely carries it forward.
- "Rung 1 is `552566d` (the last pre-merge state)" understates it: the campaign
  ran at `b5d369e`, two tooling commits later.
- "The campaign, pre-merge" is true of the 2026-10-04 merge, not of `552566d`.
