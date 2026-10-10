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

### 3. The comparand itself is not quotable on our fields

Same ring, same field (`base_fix`, z ~ 70.9 m); all three are positive for nose-down per E2:

| estimator | t = 15 s | t = 30 s |
|---|---|---|
| `section_twists_deg` (the T2-frozen S-7 construction) | -0.676 deg | -0.801 deg |
| rotation of the least-squares affine fit (`_ring_kinematics` omega) | -0.685 deg | -0.833 deg |
| chord-line rotation LE->TE (what a beam reference reports) | **+0.611 deg** | **+0.734 deg** |
| affine-fit residual RMS / (\|theta\| * r_rms) | 5.7 | 4.7 |

At z = 90 m the same three read -1.194 / -1.092 / **+0.841**. Restricting the frozen estimator to the
outer shell (`allOuterShellNods`) moves it by 0.02 deg, and the ring's node asymmetry
(`sum rx ry / sum r^2` = 0.02 .. 0.03) cannot explain the flip: **the non-affine in-plane deformation of
the section dominates the fitted rotation by ~5x**, so "the section twist" is a small difference of large
cancelling components and two standard constructions of it disagree in sign, by more than the whole
signal (~1.4 deg at 70.9 m). This is `distortion/|omega| = 1.6728` (P2c) met at the level of the
comparand. No residual against the reference may be quoted until the comparand is arbitrated.

What *is* estimator-robust: the A/B amplification measured at t = 5 s is ~6x in the band under **both**
estimators, so #30's "the amplification is structural" conclusion stands; only the absolute comparand is
blocked.

### 4. Window means, for the record (not verdicts)

`base_fix` 10-30 s window mean band value: frozen estimator **-0.747 deg**, chord-line **+0.413 deg**,
reference **+0.750 deg**. Per-station chord-line means: 0.10 deg at 45 m, 0.37 at 60, 0.60 at 70.9, 0.87
at 80, 1.16 at 90, then oscillating (0.45 at 100, 0.02 at 103.8 with a +-1.9 spread, 2.51 at 110)
against the reference's monotone 0.59 -> 0.90. The activated side has no window at all: its chord-line
twist runs +0.81 -> +3.69 deg at 70.9 m and +2.75 -> +12.82 at 100 m over 0.5 -> 5 s, i.e. pre-steady.

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

The cheap rung cannot decide the activation question, for four measured reasons: our side is not steady
at 5 s and the activated side has no window at all; the reference is itself a window mean with a +-0.6 deg
cycle; the comparand is estimator-dependent at the size of the signal, sign included; and the two sides
run at different loading. What the decision now needs, in order:

1. **an arbiter for the comparand** - a case with a known section twist under a comparable distorting
   load (a prescribed moment on this mesh, or the `beamdyn_driver` that ships in the same OpenFAST env);
   the S-7 static case is the closest existing one and does not distort enough to test the estimator;
2. **a settled activated coupled run at campaign scale** - >= 30 s with 0.5 s output, reporting the
   window mean *and* the envelope, at the reference's per-blade loading (0.812 MN / 6.34 MN.m);
3. **a diagnosis of the case5s ~5 s cycle** before any coupled number of ours is re-baselined onto it
   (#19-adjacent; `frontiersin/yaw_0` settles, so it is a case or code-state property);
4. the reference's +-0.6 deg cycle stated as the tolerance floor of any comparison.

Until (1) is closed, "activate or park" cannot be decided by a twist comparison, so the decision stays
where the 2026-10-10 branch put it: activation off, #16 open.

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
