# Feature: close the tip-twist discrepancies with a finite, strict plan

Written 2026-10-13 after the beam-vs-ring result (`shell-twist-mixed-mesh.md`, baf38a3). Tracked
upstream as an L3 item under roadmap #33 (#42). Goal: stop
the probe-by-probe loop. The plan is finite: six phases at most, two of them conditional, one
phase per session, each with an entry, a pre-registration, a control, a decision table and an exit
state. Nothing is tuned toward Zhou. Blade mesh 0.15-0.25 m.

Issue: efirvida/AeroElast#42 (linked from #33, section C, L3 row).

## Why we kept iterating

Three discrepancies were chased as one, and no step had an independent judge:

| id | discrepancy | numbers (tip, nose-down +) | question |
|---|---|---|---|
| D1 | shell vs beam, same load | shell +5.855 vs decoupled beam +2.079 (2.82x; ~2.4x per segment from the root, ring rigid at 75 %) | structural model or load path |
| D2 | our beam vs Zhou's GEBT, Zhou's load | +2.079 vs +3.60 (Energy 336:138488, Table 4) | beam class, geometric nonlinearity, or the reconstructed Mp |
| D3 | our coupled FSI vs Zhou | +5.589 vs +3.60 | follows from D1 and D2: the shell under Zhou's load alone gives +5.855 |

D3 is arithmetic once D1 and D2 are attributed. Issue #10 (parked divergence) is a different
question (dynamic stability) and is out of this plan.

## Facts the plan rests on (verified 2026-10-13)

- Deck provenance: the IEA-15 report v1.0 (Gaertner et al. 2020, `.sources/papers/75698.pdf`,
  Fig. 2-11) computed the blade properties with PreComp and plots a shear centre; the repository
  release notes say the elastic properties later "consistently originate from SONATA/ANBA4"
  (`IEA-15-240-RWT/ReleaseNotes.md:25`). Which one produced our deck file is not yet pinned (P0).
- K46 = K56 = 0 exactly at all 26 deck stations (parsed); K16/K26 (shear-torsion) are non-zero,
  max |K16|/K66 = 1.49. The exact zero suggests the bend-twist terms were dropped on export
  (inference). The deck can judge GJ and the shear centre, not bend-twist coupling.
- WU-C (`odd/tasks/shear-centre-arbitration.md`, issue #25): the deck's shear centre is within
  0.018c of its reference axis; the mesh's 0.477c claimed in a docstring is unexplained.
- S-7 is disputed (#20, open): kappa-corrected ratio 1.047 at 0.25 m (softer) vs the earlier
  0.9396 (stiffer). It is therefore NOT used as a control below.
- Zhou's GEBT uses a generic sectional stiffness matrix `C*` (their Eq. 23) and cites the IEA-15
  repository; their stiffness source, pitch and torsion definition are not stated (#14, closed as
  non-transferable). Their first torsional frequency 4.072 Hz matches ours 4.000 Hz.
- Shell implementation = CalculiX S8R within 1.46 % on the same nodal force (Line A); #34, #38
  ruled out (E3, Line B); distortion ruled out (beam-vs-ring).
- The shell load (`force_projection.py:784-805`: rotor-plane Np/Tp directions, full 3-D arm about
  `LE + 0.25(TE-LE)`, Tp included) and the beam torque (`m = Mp + ((pitch-0.25)c - x_shear) Np`,
  Tp arm zeroed, no twist projection) are not shown to be the same torque.
- `beamdyn_driver` v5 works (`tools/diagnose_radii_beamdyn.py`: static, one point load per strip
  binned from the shell's own nodal forces about the reference axis; reads displacements only).
- Papers in `.sources/papers/`: `Unsteady aeroelastic performance of the 15 MW ... surge
  condition.pdf` (Zhou, Energy 336:138488), `75698.pdf` (IEA-15 report), `jmse-13-00116.pdf`
  (ALM-GEBT, structural nonlinearity on IEA-15), `fenrg-13-1571567.pdf` (LL-FVW+GEBT, yaw).
  `AIAA2023Paper_EscaleraMendoza.pdf` models the UTD blade (disowned input, AGENTS rule 7): not a
  judge.
- Conda: `module load anaconda3/2020.02_sequana`; environments live under
  `$SCRATCH/conda-envs/` (the `openfast` env is there).

## Acceptance (from #33, binding)

- A physics comparison passes at <= 5 % against an independent reference; a wider band only with
  a number from an independent source fixed before measuring.
- Numerical references settle only in pairs. One code agreeing is not validation.
- Bands in P1/P2 below are TRIAGE (which phase runs next), not acceptance. Only P3/P4 produce
  verdicts.

## Phases

Order note: #33 attacks C items once their level's audit (B3, L3) has reached them. P0 is
documentation only and can run before; P1 onward run when the maintainer schedules them.

### [ ] P0. Provenance and pre-registration (docs only; no solver run; one session)

Entry: this file; #33 acceptance rule; #14, #20, #25.
1. Pin the deck's provenance: IEA-15 repository version in use, `ReleaseNotes.md` entry, and
   whether the ANBA4 run produced non-zero K46/K56 that the export dropped (report, release notes,
   the WindIO `elastic_properties_mb.six_x_six`). Record the answer with page/line cites.
2. Zhou (Energy 336:138488): their `C*` source, linear or nonlinear GEBT, reference axis, torsion
   definition, pitch. Record what is stated and what is not, with page cites.
3. ALM-GEBT (`jmse-13-00116.pdf`): the linear-vs-nonlinear tip-twist factor for this blade, if
   reported (the external number for P1's geometric factor).
4. Write and commit the P1 pre-registration (the table below, filled with exact commands, outputs
   and file names).
Closes when: items 1-3 recorded with cites (an "not stated" is a valid answer) and the P1
pre-registration committed; a comment on the issue lists them.

### [ ] P1. Attribution run with BeamDyn (one job, ~10 min, `sequana_cpu_dev`; one session)

Entry: P0 record. Tool: extend `tools/diagnose_radii_beamdyn.py` (add `RDzr` to the OutList,
Wiener-Milenkovic to angle; strip loads from Np/Tp/Mp) and reuse the Zhou-load shell solve.
Runs (all on the official deck, K16/K26 kept, no gravity, no rotation):
- `BD(shell_load)`: the shell's own nodal forces binned to strips about the reference axis.
- `BD(zhou_formula)`: Np, Tp, Mp per strip with the beam's arm.
- `BD(zhou_formula + Tp arm)`, `BD(zhou_formula + twist projection)`,
  `BD(zhou_formula, arm about the AC point used by the shell)`: one term toggled each.
- Each at full load (`nl`) and at 1 % load rescaled x100 (`lin`).
Control C1 (must pass before any factor is read): BeamDyn under a pure tip couple equals our own
`GJ(z)` integration of the same deck within 1 % at 50/75/90 % and tip. It tests frames, units and
the Wiener-Milenkovic conversion, without the shell. Fail -> report, stop; fixing the harness is
the next session.
Factors at the rigid 75 % ring and the tip (same estimator as `_ring_kinematics` on the shell):
| factor | definition | triage |
|---|---|---|
| structural | `omega_shell / BD_lin(shell_load)` | <= 1.05: D1 is load path, P2 skipped. > 1.05: P2 |
| load path | `BD_lin(shell_load) / BD_lin(zhou_formula)` | the toggled runs name the arm term(s) whose sum reproduces it within 5 %; otherwise "unattributed", recorded, not chased |
| geometric | `BD_nl / BD_lin` | reported beside the ALM-GEBT number from P0 |
| beam class | `BD_lin(zhou_formula) / decoupled beam` | reported; > 5 % means K16/K26 matter and the decoupled beam is retired as a reference |
D2: `BD_nl(zhou_formula)` within 5 % of 3.60 -> Zhou reproduced on its own model class, D2 closed.
Otherwise D2 is closed as non-transferable, consistent with #14 (no further work on D2).
Also reported, not gated: the tip-couple shell/BeamDyn ratio, as evidence for #20 (not closing it).
Closes when: C1 passed and the four factors plus the D2 line are committed and commented.

### [ ] P2. Name the section property (only if structural > 1.05; one session)

Entry: P1 record. Shell and BeamDyn, tip load cases, 0.25 m:
(a) pure tip torque; (b) pure tip flap couple (twist from bend-twist coupling only; the deck gives
0); (c) tip flap force at two chordwise points (the shell's shear centre near the tip vs the deck
x_shear, evidence for #25).
Prediction check: the named property, put into a COPY of the deck (labelled attribution-only, never
a reference), must reproduce the P1 structural factor within 5 %. Otherwise it is not the cause:
record "not attributed by (a)-(c)" and go to P3 with the unattributed residual.
Closes when: one property named with its prediction check, or the null recorded.

### [ ] P3a. Install the section judges (one session)

SONATA + ANBA4 (github NLRWindSystems/SONATA) and OpenSG (github wenbinyugroup/OpenSG), each in
its own conda env under `$SCRATCH/conda-envs/`.
Install gate per tool: on the official WindIO yaml, EA, EI_flap, EI_edge and GJ at three stations
(r/R 0.3, 0.6, 0.9) within 5 % of the deck. Proves both tools read the same blade.
Closes when: both gates pass, or the failing one is recorded with its error (no workaround that
changes the blade input).

### [ ] P3b. Judge the named property with the pair (one session)

Compute the property P2 named (or the full 6x6 if P2 was null) with SONATA/ANBA4 and OpenSG.
Verdict per #33: shell within 5 % of both -> the shell is right and the deck/decoupled beam is
incomplete for this blade. The two judges disagreeing by > 5 % -> not validated, recorded. Shell
off both by > 5 % -> a shell finding: propose the fix, wait for the maintainer's OK.

### [ ] P4. Close (one session)

One attribution table for D1/D2/D3. Ships, per #33: the verdict in `contract.yaml`, the closure
section in `docs/validation_closures.md`, the comment on this issue; the #33 row updated.

## Strict rules (apply to every phase)

1. One phase per session, from a clean context. Resume recipe: this file's first unchecked box, the
   issue's last comment, one `mem_search "twist-discrepancy-closure"`.
2. The phase's pre-registration (numbers, bands, commands) is committed before any tool is edited.
3. A failed control stops the phase. No reinterpretation, no second hypothesis in the same session.
4. A surprise becomes a new issue under #33 (section C), never absorbed into the phase.
5. Nothing is tuned toward Zhou. Deck edits exist only in P2's prediction check, labelled.
6. Every phase ends with: commit(s), this file updated, Engram mirror updated, a comment on the
   issue with the phase verdict and the next phase's entry.
7. Budget: at most two job submissions per phase (one run, one harness-fix rerun); context
   ~200 k -> hand off.
8. Out of scope: parked FSI aeroelast vs ccx_preCICE (#10); `k_drill` units (#38); wall-flow (#16);
   0.5 m re-calibration; any Rust change.
9. Environment: `scripts/aeroenv.sh` with `PYTHONPATH=$PWD/src`, never `python -I`; `.sources`
   linked; mesh 0.15-0.25 m.

## Why the parked FSI comparison is not in this plan

At pitch 90 the load is near zero (CT ~ 0.008), so bend-twist barely activates; CCX S8R shares the
CLT and the geometry writer with our shell and already agrees to 1.46 %; D3 is explained by
structure, not by the coupling layer. It is the right tool for #10 (settle vs diverge), on its own
track (`ccx-fsi-parity-parked-v50.md` T5-T6), after its blocker: the `*SHELL COMPOSITE` deck
diverges under load; move to `*SHELL GENERAL SECTION` (A/B/D) and re-gate. Cost today ~9 h per 40
windows at 0.5 m, 0 windows at 0.25 m.
