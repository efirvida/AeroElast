# Feature: the coupled rotor path diverges at t ≈ 1.9 s

Status: projector load-frame fix **committed as `7a84da1`** and unit-verified
2026-10-08; it is NOT the #19 cause; the anchor holds (`campaign` 500/500) and
the origin/main-side ladder runs (jobs 11611050-11611052)
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
- [ ] T6i Attribution among the constructor's deformation-dependent outputs
      (jobs `11611101` `frzsign`, `11611102` `frzac`, `11611103` `frzgrid`), each
      one keeping HEAD's rebuild and pinning a single attribute to the reference
      projector afterwards: `_strip_moment_axis_sign` (a discrete +/-1 from the
      LE/TE blunt rule - the one candidate that can *jump*), `_strip_ac_offsets`
      (the AC->centroid moment arm) and `_strips` (node-to-strip assignment plus
      the strip widths and offsets). A pin that restores contraction names the
      physical quantity the fix has to keep deformation-independent; note the
      code already states the intent for its own neighbour - `_strip_node_indices
      is intentionally kept as the *reference* assignment and is NOT updated from
      the deformed projector here` - and the projector reassigns internally
      anyway.
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
- [ ] T7 Re-run the 30 s gate (`tests/run_step1b_smoke.srm`) and re-anchor the
      campaign baseline once T6 is clean.
- [ ] T8 Reconcile the store (`rotor_coupling_noncontraction`,
      `force_projection_sense_p5`), the group 27/35 rows and
      `docs/validation_closures.md`; the consequence numbers above are the
      record.
- [ ] T9 Comment and close #19, and fold #26's remaining closure criteria into
      it.

### Corrections to the issue's own text

- The window is not "the 2026-10-04 merge plus the #11 batch": `9a3923e` is
  2026-10-02, before the merge. The merge merely carries it forward.
- "Rung 1 is `552566d` (the last pre-merge state)" understates it: the campaign
  ran at `b5d369e`, two tooling commits later.
- "The campaign, pre-merge" is true of the 2026-10-04 merge, not of `552566d`.
