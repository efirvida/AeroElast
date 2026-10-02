# Feature: refresh the validation matrix against the new `main` line

**Status**: in progress
**Opened**: 2026-09-30
**Branch**: `integrate/origin-main-2026-09-30`
**Engram mirror**: `odd/validation-matrix-update/tasks`
**Predecessor**: `odd/tasks/integrate-origin-main.md` (the first merge + the measured before/after)

## Goal

Have **every** change of `origin/main` in this tree, then refresh the recorded
validation matrix (G1 elements, G2 structural solvers, G3 BEM, G4 FSI campaigns)
against the new element line, and report per item whether the behaviour improved or
worsened.  Launch order is by computational cost, cheapest first.

## User decisions (2026-09-30, second round)

1. **First**: make sure we have all of `main`.  The ten commits between `dffacab`
   (already merged) and `5bfa2b2` are still missing.
2. **Then**: start launching the cheapest items first.
3. **On the G1/G2 deltas** (box 3.27%, D-Tube 4.75%, UL elastica 3.63%, V-02 2F
   −6.3%, S-7 1.406): **investigate each delta against an independent reference
   before touching any band.**  No re-baselining by tolerance widening.

## Known cost ladder (to be refined with real numbers)

| Cost | Item | Instrument | Needs |
|---|---|---|---|
| minutes | G1 unit anchors (cantilever, box, D-Tube, elastica, MITC4+) | `pytest` | nothing |
| minutes | G2 anchors (V-02, V-03, V-05, S-1, S-4, S-6, S-7, S-8c, Campbell) | `pytest` | nothing |
| ~40-50 min | Blade structural matrix B1..B7 @ 0.25 m, **both sides** | `tools/run_blade_structural_matrix.py` | a pre-merge build |
| ~1 h | S-1/S-7 sectional + GJ matrix on both sides | `tools/run_s1_sectional.py`, `tools/run_s7_torsion.py` | nothing |
| hours | G3 BEM V-01 vs AeroDyn (6 points) | `tools/run_s5_oneway.py`, BEM engine | CCBlade |
| hours-days | G4 FSI: yaw sweep 0-40°, V-06 15 cases, parked DLC, B1/B2 h/dt | Slurm + preCICE + OpenFOAM | cluster |

## The ten commits being merged (`dffacab..5bfa2b2`)

All test hygiene, 2026-09-30, none touches elements, mesh generators or the CCX writer:

- `5bfa2b2` lands ten measured mismatches as documented `pytest.xfail`
  (reports "Blade: 10 passed, 8 xfailed" — i.e. CalculiX integrates the blade
  **on their side**; see the open lead in the predecessor document)
- `53925ca` deletes rotor tests that cannot fail; `38c88f3` replaces two
  assertions that could not fail; `7aa6948` asserts M-orthogonality instead of a
  non-zero norm; `ff78b9a` gives `test_polar_cl_not_constant` a real lift check
- `dabbe36` tightens previously widened tolerances to the real 5% bound
- `c904dc20`→`904dc20` adds the test→reference map; `e3c9724` reverts `c9f8fe8`
- `20229fa` documentation

## Risk of this merge (accepted, recorded)

`5bfa2b2` turns measured differences into `xfail`: if it lands before S-4/S-6/2790
are fixed, part of the signal becomes an expected failure instead of a red test.
The user asked for all of `main` first, so it lands now and the matrix must carry
the distinction explicitly (`xfail` = "we know and have named the gap").

## Tasks

- [x] T1 — Fetch and merge `origin/main` @ `5bfa2b2`; resolve conflicts; rebuild the
      Rust extension; verify the import surface.
- [x] T2 — Refresh the AFTER column of the G1/G2 anchors (pytest) and record which
      items became `xfail`.
- [x] T3 — Launch the cheapest matrix items in cost order, recording real durations.
- [x] T4 — Investigate the G1/G2 deltas against independent references (no band edits).
- [ ] T5 — Refresh `docs/origin_main_integration_2026-09-30.md` (or a successor) with
      the per-group verdict and the campaign launch plan.

## Progress

### 2026-10-01 — G4 yaw + h/dt verdict, and its confound

The first two G4 campaigns finished (yaw 0-40, 100 s each; h-coarse/fine and
dt-coarse/fine, 30 s each) and the snapshot was diffed with
`tools/campaign_metrics.py`.  To make that possible the tool gained `--alias NEW=OLD`
because the re-run directories do not carry the recorded campaign names, and its
delta was switched to be relative to the *before* value; a unit test pins both
(`tests/test_campaign_metrics.py`, 8 passed).

**Finding (the verdict's key caveat): the before/after is config-confounded.**
The after campaigns started 2026-09-30 22:55, after `8e6488e` (20:30) changed 30
solid case YAMLs (`airfoil_spacing: constant -> cosine`) and stopped honouring an
explicit `hub_radius: 0.0` when the aero has a hub (a 3.97 m strip shift).  The
recorded `_mitc3fix` campaigns ran the old config.  `Thrust/CT` (= `q·A`) dropped a
uniform 4.76% across all five yaw cases -- a global, yaw-independent shift, not an
element effect; the element's static effect is already measured at <=0.6% on the
blade matrix.  So the yaw delta reads flap -8.5% to -7.7%, thrust ~-5%, power -4.8%
to -3.9%, CP/CT flat: a *config* before/after, not an element before/after.

Convergence, read within the re-run: dt 0.02 vs 0.005 closes at 0.003%; h converges
second-order (0.5 -> 0.25 = 4.6%, 0.25 -> 0.125 = 0.84%, observed order ~2.4, 0.125 m
within ~0.2% of the Richardson limit).  The `convergence_b1_b2_results_official`
column (flap 17.84 m) is a different geometry and does not read as convergence.

Full write-up: `docs/origin_main_integration_2026-09-30.md`, section "G4 FSI
campaigns: the finished before/after".  Open follow-up: a clean element A/B needs one
campaign re-run with the case config held fixed.

Upstream moved again while this ran: `origin/main` `5f189fa -> 31565b3` (13 commits,
new `validation-2026-09` tag); merge-tree is clean, but do not rebuild until the
in-flight G4 jobs finish.

### 2026-10-01 (later) — the FSI bias traced to the shell twist; issue #9 and the
### campaign decision

Chasing why the FSI over-de-loads, the mechanism is the shell's elastic twist feeding
the BEM's deformed-geometry feedback: rigid BEM 16.51 MW / 2.54 MN (matches Zhou rigid)
vs FSI 12.96 MW / 1.73 MN, i.e. **-21%/-32% against Zhou's -8%/-13%**. The extraction
is clean (a controlled field test: pure flap -> 0.3 deg, real torsion -> exact), so the
twist is a real FEM response: ~8.9 deg vs BeamDyn's 0.98 deg (`docs/twist_distortion_diagnosis.md`,
open since 2026-09-16).

A controlled flat coupon (`tests/test_laminate_bend_twist.py`, `xfail`) gives shell/CLT
**3.4-4.8x** (0 when `D16 = 0`), but the same shell matches **CalculiX S8R to 0.4%** on
the same strip -- so the reference is in dispute, not confirmed as an element bug.
Filed upstream: **issue #9** (`efirvida/AeroElast`) plus a comment with the exact layup
sweep reproducer.

**Campaign decision** (the defect scales with aerodynamic load): normal-operation cases
(yaw, ch6 1.1, V-06) have `CT ~ 0.5` and are biased; the feathered/extreme ones (ch6
6.1/6.3, parked) have `CT ~ 0` and are usable. Cancelled **every** remaining G4 job:
the queued V-06 (`11604992/93`, no compute spent) and, after analysing them, `ch6`
(`11605077`) and `parked` (`11605081`). The parked run diverges (`Tip Disp Y` > 2000 m)
where the recorded response settles at 4.705 m -- filed as **issue #10**. Convergence
h/dt is unaffected (discretisation).
