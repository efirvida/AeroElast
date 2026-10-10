# Feature: the coupled arbiter for the wall-flow activation decision (P2b, feeds #16)

**Trigger.** Issue **#16** stays open on one question: whether
`bem.wall_flow_realisation` should default to `true`. The realisation is implemented, guarded
(group 37) and parked off; the mesh defect it denounced is fixed (0/671 -> 671/671) and the
coupled movement is attributed (#30: structural response to the pattern, not a load or a
transfer effect). On 2026-10-10 the user chose the third branch of the decision:
**do not activate and do not close #16 until the coupled response has an independent arbiter.**

**That arbiter is already named in the record**: `odd/tasks/production-path-and-independent-arbiters.md`
§P2c assigns *"is our coupled response right?"* to **P2b — OpenFAST's flexible rated response**,
still marked `[ ]` / *"not yet run"* (§P2b). This document is the P2b work unit.

**Scope of this unit.** Produce the arbiter and define what it decides. It does **not** flip
`bem.wall_flow_realisation` and does **not** close #16; those are consequences of the number.

**Boundary.** No `src/` or `crates/` change. No validation-store row until the comparison is
actually measured against our side. No change to `docs/validation_closures.md` until then.

---

## T1 — inventory and feasibility (DONE 2026-10-10, all evidence below is measured)

### What exists

| piece | path | state |
|---|---|---|
| OpenFAST binary | `/scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast` | present, 43.4 MB; the run that produced the reference `.out` reports `v5.0.0`, built Apr 11 2026, single precision; `aerodyn_driver`, `beamdyn_driver` sit beside it |
| converted v5.0 glue deck | `/scratch/leahk/eduardo.donestevez/ofruns/OpenFAST/IEA-15-240-RWT-Monopile/` | complete case: `.fst`, `_ElastoDyn.dat`, `_ElastoDyn_tower.dat`, `_AeroDyn15.dat`, `_ServoDyn.dat` + `_DISCON.IN` + `_ROSCO.yaml`, `_SeaState.dat`, `_HydroDyn.dat`, `_SubDyn.dat`, and a finished `IEA-15-240-RWT-Monopile.out` |
| shared model dir | `/scratch/leahk/eduardo.donestevez/IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT/` | airfoils, `_AeroDyn15_blade.dat`, `_ElastoDyn_blade.dat`, **`_BeamDyn.dat`**, `_BeamDyn_blade.dat`, `_InflowFile.dat`, `ServoData/`, `Wind/`; the case dirs reach it as `../IEA-15-240-RWT/` |
| BeamDyn run (diagnostic) | `/scratch/leahk/eduardo.donestevez/tmp/opencode/ofrun-s6-bd/IEA-15-240-RWT-Monopile/` | ran 100 s on 2026-09-15: 33.8 MB `.out`, 6.2 MB `.outb`, three `*.BD.R1.B{1,2,3}.ech`, and (corrected in T3) it **does** carry the 51 `B1N###_RDxr` torsion channels; **`tmp/`**, i.e. the location the standing preference forbids for anything needed |
| in-repo deck | `tests/reference/iea15mw_openfast/` | **AeroDyn-only**: `case/IEA-15-240-RWT_AeroDyn15.dat`, `IEA-15-240-RWT/IEA-15-240-RWT_AeroDyn15_blade.dat`, 50 polars, `NOTICE`. No `.fst`, no ElastoDyn/ServoDyn/InflowWind. `tests/support/openfast_bem.py:47` resolves it from `tests/support/paths.py:24` (`DATA_DIR = tests/`) |

### The finding: the deck in use has no blade torsion

`ofruns/OpenFAST/IEA-15-240-RWT-Monopile/IEA-15-240-RWT-Monopile.fst:18` is
`CompElast = 1` — **ElastoDyn**, `CompAero = 2`, `CompServo = 0`, `TMax = 100.0`,
`DT = 0.005`. Its `_ElastoDyn.dat` sets `FlapDOF1/2 = True`, `EdgeDOF = True`
(`:8-10`), `DrTrDOF`, `GenDOF`, `TwFADOF1/2` false (`:13,14,16,17`), `RotSpeed = 7.56`,
`TipRad = 120.97`, `HubRad = 3.97`, `BldFile{1,2,3} = ../IEA-15-240-RWT/IEA-15-240-RWT_ElastoDyn_blade.dat`.

ElastoDyn's blade is a **linear modal beam with two flapwise and one edgewise mode and no
torsion degree of freedom**; blade torsion in OpenFAST needs **BeamDyn** (`CompElast = 2`).
So the S-5 reference run (`TipDxc1 = 16.03 m`, recorded in Engram
`openfast-5-0-0-corriendo-iea-15-mw-deck-v4-2-v5-0-convertido`) is a **flap/edge** number:
**no OpenFAST number produced by this deck can arbitrate a twist**, in either direction.
That is why P2b has never run: it is not a scheduling gap, the deck lacks the DOF.

### The BeamDyn run that does exist was not usable as it stood

`tmp/opencode/ofrun-s6-bd/.../IEA-15-240-RWT-Monopile.fst:18` is `CompElast = 2` with
`BDBldFile(1..3) = ../IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat` (`:42-44`) — the flexible
torsion model — and it ran (`bd_smoke_11594900`). Its three defects are *not* the ones this
section first claimed; T3 measured them:

1. **It lives in `tmp/`.** The only flexible-torsion reference on the machine sat in the
   diagnostic location the standing preference forbids for anything needed.
2. **Its model tree is not the named one.** The case reaches `../IEA-15-240-RWT/`; that
   sibling carries `HWindSpeed = 10.59` and matches `$SCRATCH/bfs16/case5s`, while the
   upstream clone `/scratch/leahk/.../IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT/` carries
   `10.0`. Taking the named clone would have moved the operating point off rated — silently.
3. **Its AeroDyn deck is stale at source.** The case's
   `IEA-15-240-RWT-Monopile_AeroDyn15.dat` (mtime 2026-09-15 23:04, *after* the 16:28 run)
   lists `"B1Mp"`, which this build rejects: `SetOutParam:B1Mp is not an available output
   channel`. A re-run from the source as found fails in setup.

**Correction to this document's first reading.** The `NNodeOuts = 0` claim was measured on the
*upstream clone's* `IEA-15-240-RWT_BeamDyn.dat:87`, not on the case's own sibling, which
already had `NNodeOuts = 1` and the nodal `OutList`. The torsion channels were therefore
present in the 2026-09-15 `.outb` all along, and the defect was **durability and provenance**,
not missing output. Recorded because the wrong diagnosis would have sent the next reader to
edit a deck that did not need editing.

### What T1 left for T2/T3

A durable, channel-carrying flexible-torsion reference:
1. rebuild the BeamDyn case in `$SCRATCH` (not `tmp/`, not the repo) with the shared model dir
   reachable, and `NNodeOuts >= 1` (`1..9` allowed) plus the torsion `OutList` channels on the
   BeamDyn blade outputs;
2. re-run rated to `TMax = 100 s` and keep the `.outb` and the exact command;
3. record the operating point and confirm it matches the case #16's A/B used
   (`$SCRATCH/bfs16/case5s`), because a comparison across different operating points is not a
   comparison.

## T2 — the arbitrable quantity, fixed BEFORE either profile is read

Rules, from `odd/tasks/production-path-and-independent-arbiters.md` §P2c: a number without its
arbiter is not a number, and the definition is fixed while blind to the result.

* **Our side.** `tools/run_s7_torsion.py::section_twists_deg(mesh_model, n_slices=40)` — the
  section twist of the shell, band-limited to `0.3*117 .. 0.9*117` m, as the S-7 test and the
  tool use it. The *coupled* variant, when the campaigns side is measured, is the same
  extraction applied to the coupled run's fields.
* **OpenFAST side.** The BeamDyn blade-1 nodal rotation about the local blade axis at the same
  span stations — `B1N###RDxr` (BeamDyn "rotation" channel) with `NNodeOuts >= 1`, time-averaged
  over a converged window at the end of the run, in the same blade-root frame.
* **Sign and reference.** Frozen by `odd/tasks/rated-twist-sign-convention.md` (E1: `omega` and
  `theta_z` are the same signed quantity; E2: the Rodrigues frame derivation fixes the sense at
  `omega > 0`). The conversion from the BeamDyn channel to the same signed section twist, and
  the span stations, are written into the probe **before** the run's numbers are read.
* **Explicitly NOT decided by this arbiter.** (a) the sectional stiffness — that is S-7 /
  group 34; (b) the load application pattern — group 35; (c) comparability with Zhou 2025 —
  nothing, per #14; (d) the shell's section distortion: OpenFAST's blade is a beam, so the
  arbiter bounds the **beam-comparable part** of the coupled twist, not the distortion that a
  beam cannot represent (`distortion/|omega| = 1.6728` on the rated path).

## T3 — build and run the durable reference (DONE 2026-10-10)

Regenerator: `tools/openfast_flexible_rated_reference.py` (new, 557 lines, `ruff check` and
`ruff format --check` clean). It materialises the durable case, patches the case-local decks
only, runs the binary, gates on the channels being present in the written file before reading a
number, and extracts the profile.

| item | value |
|---|---|
| durable root | `$SCRATCH/bfs16/openfast-flexible/` (`IEA-15-240-RWT-Monopile/`, the `IEA-15-240-RWT/` sibling it needs, `runs/{short,full}/`, logs) |
| command | `cd <repo> && scripts/aeroenv.sh python tools/openfast_flexible_rated_reference.py both` |
| the run itself | `cd $SCRATCH/bfs16/openfast-flexible/IEA-15-240-RWT-Monopile && /scratch/leahk/eduardo.donestevez/conda-envs/openfast/bin/openfast IEA-15-240-RWT-Monopile.fst` |
| profile | `$SCRATCH/bfs16/openfast-flexible/blade1_torsion_profile.csv` (51 rows) |
| channel | `B1N###_RDxr`, nodes `B1N001..B1N051`, BeamDyn "rotational displacement in X, rad"; the `.out` unit row prints `-` and the `rad` comes from the BeamDyn registry |
| stations | the `BD_Blade_R1B1_Reference.vtp` output-node polyline, 51 points, `span_m` = cumulative arc length from `B1N001`, total **117.1487 m** against the summary's `Length: 117.149 m` |
| window | `90 .. 100 s`, `DT_Out = 0.05 s`, 201 samples, `TMax = 100 s` |
| rotor speed (mean) | **7.5599997 rpm** (`case5s` omega 0.7906341464750989 rad/s = 7.5486 rpm) |
| pitch | **0.0 deg** (`BldPitch1`) |
| power | **15.065 MW** rotor aero power; there is **no `GenPwr`** in this deck (`CompServo = 0`, `ServoFile = "unused"`) |
| wall time | 18.1 s at `TMax = 3 s`; **307.2 s = 5.12 min** at `TMax = 100 s`, both on the login node, no queue |
| reproducibility | `--skip-run` re-extraction rewrote the CSV to an identical md5; against the archived 2026-09-15 run (same window) `max |dRDxr| = 6.2e-8 rad`, max relative `4.1e-6` (single precision) |

The profile is monotone to ~0.01578 rad at `span = 102.46 m` and flat to the tip
(`117.15 m`: 0.0157281 rad = **0.901 deg**), raw — no sign, frame or unit conversion applied.

**Sanity checks that the deck can and cannot give.** `TipDxc1` is `INVALID` under
`CompElast = 2`, so the S-5 cross-check came from the ElastoDyn deck instead: full-run mean
15.948 m, 90-100 s mean 15.995 m, max 16.80 m — consistent with the recorded 16.03 m as a
late-run mean, not a peak. Flagged, not forced.

**No ROSCO was needed**: `CompServo = 0`; no `libdiscon.*` exists under
`/scratch/leahk/eduardo.donestevez` within depth 6, so turning ServoDyn on later is a missing
input to hunt, not a step to assume.

## T4 — register the arbiter (DONE 2026-10-10)

Registered in `odd/tasks/production-path-and-independent-arbiters.md`: §P2b is now `[x]` with the
reference's location, its regenerator and the ElastoDyn-has-no-torsion finding (which sharpens the S-5
claim: that deck's flexibility is flap/edge only); §P2c's coupled row and rule 4 now carry the T5
outcome — the reference exists, the *comparand* is what is missing. No store row and no
`validation_closures.md` section: there is no comparison to register until a comparand exists and our
side is settled.

## T5 — the cheap rung, measured: it cannot decide (DONE 2026-10-10)

T5 asked whether already-saved coupled fields (or the frozen-field static rung) could substitute for a
campaign-scale run. Every number below comes from the runs' own saved `fields.vtu` /
`rotor_performance.csv` (production path, no test-local construction) and from T3's reference. Probes
live outside the repo, as usual: `$SCRATCH/bfs16/wf_p2b_diag/probe_*.py`
(`twist_vs_openfast`, `run_timeseries`, `cycle_mean_vs_openfast`, `two_estimators_vs_openfast`).

### What exists to compare

| piece | where | state |
|---|---|---|
| our side, default (min-norm, pre-activation) | `$SCRATCH/smoke_fix_results/base_fix` | 60 steps at 0.5 s = **30 s**; identical to the #16 `A` run over its first 5 s (same values to 6 digits) |
| our side, activated (wall flow) | `$SCRATCH/bfs16/B` | **5 s only** |
| reference | `$SCRATCH/bfs16/openfast-flexible/blade1_torsion_profile.csv` | 90-100 s window mean of `B1N###_RDxr` |

### 1. Our coupled answer is a sustained cycle, and the A/B runs end on its maximum

`base_fix`, 3000 windows, every window "converged" in 2 solid iterations: `sum F_y` **0.653 .. 1.160 MN**,
`Q_y` **-6.65 .. -3.35 MN.m**, `max|U|` **11.4 .. 24.3 m**, period **~5 s**, with no decay between t = 15 s
and t = 30 s. Both `A` and `B` stop at **t = 5 s, the cycle's load maximum**: a snapshot is a phase of the
cycle, not the response. (`frontiersin_results_corotational_100s/yaw_0`, 100 s, older line and lower load
0.546 MN/blade, settles by ~5 s with a 0.4 m / 1 % residual ripple, so a settled coupled answer *is*
reachable in this code family; the case5s family's cycle is a property of that case or code state and
wants its own diagnosis.)

### 2. The reference oscillates too

`B1N046_RDxr` window means [rad]: `[0,1] 0.01094`, `[2,3] 0.02995`, `[4,5] 0.00411`, `[5,10] 0.01292`,
`[20,30] 0.01227`, `[50,60] 0.01646`, `[90,100] 0.01578` - a **+-0.6 deg torsional limit cycle** about the
+0.75 deg the profile carries. So the reference is usable only as a window mean, and any
snapshot-to-mean comparison carries an error bar the size of the signal.

### 3. The comparand: the S-7 construction survives, and it took a control to see it (corrected 2026-10-10)

**Correction, and it is mine.** The first pass of this section reported that "the frozen estimator and the
chord-line rotation of the same ring disagree in sign", which would have made the comparand unquotable.
That was a **bug in the chord-line probe**, not a property of the field: the implementation computed
`atan2(c' x c, c' . c)`, the rotation *from the deformed chord to the reference one*, which returns
`-theta`. The arbiter's control exposes it -- see §T6: on the S-7 case, where the answer is known, the
four estimators give `frozen 0.9396`, `affine 0.9932`, `rigid 0.9395` and (with `c x c'` the right way
round) `chord 0.9274` of the beam. With the sign fixed the coupled field gives `frozen -0.676` and
`chord -0.611` at z = 70.9 m: **same sign, 10 % apart.** No residual was ever quoted from the inverted
numbers, and the two figures that did escape into §4 below are corrected here.

What survives as a real measurement is smaller and different:

* the estimators agree in the mid band (within ~4 % at 70.9 m) and **diverge in the outer span**, where the
  section distorts most: at t = 15 s, z = 98.0 m, `frozen -1.92` against `chord -0.14`, and at 103.8 m
  `frozen -1.50` against `chord -1.90`;
* the ring's motion is strongly non-affine (affine-fit residual RMS 5.7x the rotation term at 70.9 m), so
  the single scalar is a model, not a direct reading -- but on the control, where that same machinery is
  exercised against a known beam answer, all four constructions land within 7 % of it. The non-affinity
  therefore does not by itself invalidate the S-7 construction.

### 4. Window means, for the record (not verdicts)

`base_fix` 10-30 s window mean band value: frozen estimator **-0.747 deg**, chord-line **-0.413 deg**
(the +0.413 of the first pass was the sign bug), reference **+0.750 deg**. Per-station chord-line means
(sign corrected): -0.10 deg at 45 m, -0.37 at 60, -0.60 at 70.9, -0.87 at 80, -1.16 at 90, then
oscillating (-0.45 at 100, -0.02 at 103.8 with a +-1.9 spread, -2.51 at 110) against the reference's
monotone 0.59 -> 0.90. The activated side has no window at all: its chord-line twist runs -0.81 -> -3.69
deg at 70.9 m and -2.75 -> -12.82 at 100 m over 0.5 -> 5 s (sign corrected), i.e. pre-steady.

### 5. The two sides are not at the same loading

| | reference (OpenFAST) | our coupled case (case5s) |
|---|---|---|
| rotor speed / pitch | 7.560 rpm / 0.0 deg | 7.55 rpm / 0.0 deg |
| thrust per blade | 2.436 MN rotor / 3 = **0.812 MN** | `sum F_y` = **0.900 MN** mean (0.653-1.160) |
| torque per blade | 19.03 MN.m rotor / 3 = **6.34 MN.m** | `Q_y` = **5.55 MN.m** mean (3.35-6.65) |
| rotor aero power | **15.065 MW** | **13.17 MW** (single blade x3) |

The reference is at the official rated point; our side sits at ~1.11x thrust and ~0.88x torque per blade
with a +-28 % swing, so a torsion residual measured there would confound the realization question with a
different aeroelastic state.

### Verdict: the hand-back to #16

The cheap rung still cannot decide the activation question, but the list is shorter and better aimed
after the arbiter's first run:

1. **Our coupled answer is not a state.** `A` and `B` stop on the load maximum of a sustained ~5 s cycle,
   and the activated side has no window at all (5 s only).
2. **The reference is itself a window mean** with a +-0.6 deg torsional limit cycle about its band value.
3. **The two sides are not at the same loading** (~1.11x thrust and ~0.88x torque per blade, 13.17 against
   15.065 MW, +-28 % swing).
4. **The one beam reference available for the coupled field disagrees with our shell in sign and by ~8x**:
   with the field's own internal torque `M_z(z)` (nose-down positive, +0.56 MN.m at the root at t = 15 s)
   and the deck's `GJ`, a beam twists `+5.09 deg` (band mean, `+12.7` at 103.8 m) while our shell reads
   `-0.65` (frozen) / `-0.28` (chord). This is **not** an estimator problem (the control validates the
   estimators against the beam at 0.94) and it is **not** a stiffness problem (the same control says
   0.94): the aero internal torque is nose-down and our shell twists nose-up, so the missing term is the
   one the beam reference does not carry either -- the **centrifugal/inertial torsion**, which by this
   reading dominates the coupled blade's torsion path. That is the next thing to measure, not the
   comparand.

So the arbiter's order changed: **(a)** quantify the centrifugal torsional moment on this mesh and check
its sign and size against the aero torque above -- if it dominates, the comparison is a different
question than the one T2 posed; **(b)** only then re-open the comparand, whose only remaining live issue
is the outer-span estimator divergence (98-104 m); **(c)** the settled activated run at campaign scale is
still needed for any of it to become a coupled number; **(d)** the reference's +-0.6 deg cycle stays the
tolerance floor. The activation decision still sits where the 2026-10-10 branch put it: activation off,
#16 open.

## T6 — the arbitrer's first run: a control, and a beam reference that refuses (2026-10-10)

Probe: `$SCRATCH/bfs16/wf_p2b_diag/probe_estimator_arbiter.py`; log `logs/tmp_p2b_arbiter_fixed.txt`.

**The control, where the answer is known.** `tools/run_s7_torsion.py::run_torsion_case` on the campaign
mesh (0.25 m) applies a prescribed tip couple of `1000 N.m`; the beam answer is `M * integral(1/GJ)` from
the BeamDyn blade deck. Four estimators of the section rotation, ratio to that beam over the 0.3-0.9 band:

| estimator | mean estimator/beam | min .. max |
|---|---|---|
| `section_twists_deg` (the S-7 / T2 construction) | **+0.9396** | +0.79 .. +1.24 |
| affine-fit omega (`_ring_kinematics`) | +0.9932 | +0.80 .. +1.40 |
| rigid-only fit (`rigid_rotation`) | +0.9395 | +0.79 .. +1.24 |
| chord-line rotation (unit-tested against the same control) | +0.9274 | +0.73 .. +1.22 |

That reproduces the recorded S-7 anchor (0.939603) exactly, and it is what caught the sign bug of the
chord implementation described in §T5-3. The instrument is now controlled.

**The coupled field, with the same beam machinery.** The internal torque profile comes from the run's own
aero field (`M_z(z)`, cumulative outboard moment about the span axis, nose-down positive per the sign
convention); the beam twist is `integral_0^z M(z')/GJ(z') dz'`, the S-7 construction with a varying
moment. At t = 15 s: `M_z` +0.5565 MN.m at the root, +0.0076 near the tip; beam band value +5.0946 deg
(+12.69 at 103.8 m). Our shell: frozen -0.6545 deg, chord -0.2761 deg band mean. At t = 25 s: beam
+2.0025 deg, frozen -0.8369, chord -0.6804.

Read with care, because three things are wrong with this comparison and only one of them is a defect:
the beam reference carries **no centrifugal or inertial torsion** (so it is only the aero path), the
`GJ` of the deck's outermost stations is degenerate (the S-7 helper already excludes the last station,
and the remaining taper still weights the outer span heavily), and the shell's own model is a shell.
What it establishes is narrower than it looks and still useful: **the aero internal torque is nose-down
and our shell answers nose-up**, so on this load case the shell's section rotation is not produced by the
aero torque at all -- something bigger and opposite acts on the torsion path. The candidate is the
centrifugal term (`include_centrifugal: true` in the case), which a beam comparison of this form cannot
arbitrate because it is missing on the beam side too.


## Traps

- `scripts/aeroenv.sh` for anything importing `aeroelast`; OpenFAST itself needs its own env
  (`conda-envs/openfast`) and ROSCO's `libdiscon.so` for the rated controller.
- The case dirs reach the model tree through `../IEA-15-240-RWT/` — a copy that forgets that
  relative path fails in setup, not in the solver.
- `NNodeOuts = 0` is the default in the shared `_BeamDyn.dat`: a re-run without editing it
  produces a big file with no torsion in it, which reads exactly like a successful run.
- `openfast_toolbox.io.FASTOutputFile` reads `.outb`; use its `channels`/`data` attributes, not
  `channel_names` (not present in the installed version).
- `rtk` hides and repeats lines: verify listings with counts and read from a file, as used here.
- The tree is shared with another session: `git add -- <paths>` **and** `git commit -- <paths>`.
- `smoke_fix_results/base_fix` and `bfs16/A` are **the same run**: the first 5 s agree to 6 digits. Do not
  count them as two measurements of the min-norm realization.
- On a distorting shell, "the section twist" is not unique: `section_twists_deg`, the affine fit's
  rotation and the chord-line rotation are three different numbers on the same ring, and here they
  disagree in sign. Name the estimator, never quote one as *the* twist (the same discipline as "never
  widen a pass band").
- The A and B runs stop at t = 5 s, which is the load maximum of the ~5 s cycle: a phase, not a mean.
- `RotThrust` is kN and `RotTorq` kN-m in this deck, and `TipDxc1`/`OoPDefl1` are `INVALID` under
  `CompElast = 2`; the blade-1 nodal **force** channels are not written, so the reference's per-blade
  load split is the rotor channel divided by 3, and it is stated as such.

## Commit plan (work units)

`docs(odd):` T1+T2 (this file), then one commit per task as it lands. Branch
`integrate/origin-main-2026-09-30`.
