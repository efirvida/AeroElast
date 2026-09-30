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

- [ ] T1 — Fetch and merge `origin/main` @ `5bfa2b2`; resolve conflicts; rebuild the
      Rust extension; verify the import surface.
- [ ] T2 — Refresh the AFTER column of the G1/G2 anchors (pytest) and record which
      items became `xfail`.
- [ ] T3 — Launch the cheapest matrix items in cost order, recording real durations.
- [ ] T4 — Investigate the G1/G2 deltas against independent references (no band edits).
- [ ] T5 — Refresh `docs/origin_main_integration_2026-09-30.md` (or a successor) with
      the per-group verdict and the campaign launch plan.
