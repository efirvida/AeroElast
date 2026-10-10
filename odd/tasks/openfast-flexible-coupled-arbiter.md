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

## T4 — register the arbiter

A note in `odd/tasks/production-path-and-independent-arbiters.md` §P2b/§P2c: the arbiter now
exists, where it lives, how it is regenerated, and the ElastoDyn-has-no-torsion finding (which
sharpens the S-5 claim: that deck's flexibility is flap/edge only). No store row, no
`validation_closures.md` section yet: there is no comparison to register until our side is
measured.

## T5 — hand back to #16

State, in one place, what #16's activation decision now needs: either our coupled side (a
campaign-scale run with the activation on) or the frozen-field static rung already instrumented
in `tools/diagnose_wall_flow_static_ab.py`, compared profile by profile against T3's reference.

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

## Commit plan (work units)

`docs(odd):` T1+T2 (this file), then one commit per task as it lands. Branch
`integrate/origin-main-2026-09-30`.
