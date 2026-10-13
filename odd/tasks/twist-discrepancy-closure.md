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

- Deck provenance and version remain bounded by the P0 evidence below: the copied repository
  describes itself as v1.1 and contains release notes through v1.1.6, but the test deck has no
  version or raw section-solver identifier. Do not identify the report's PreComp/VABS results with
  this particular exported matrix without a trace.
- K46 = K56 = 0 exactly at all 26 deck stations (parsed); K16/K26 (shear-torsion) are non-zero,
  max |K16|/K66 = 1.49. Export loss of non-zero ANBA4 terms is **unproven**, not an established
  cause. The exported deck specifies no direct bend-twist stiffness through K46/K56.
- WU-C (`odd/tasks/shear-centre-arbitration.md`, issue #25): the deck's shear centre is within
  0.018c of its reference axis; the mesh's 0.477c claimed in a docstring is unexplained.
- S-7 is disputed (#20, open): kappa-corrected ratio 1.047 at 0.25 m (softer) vs the earlier
  0.9396 (stiffer). It is therefore NOT used as a control below.
- Zhou's nonlinear GEBT uses a sectional `C*` (Eq. 23); its *source*, axis location and rated
  pitch are not specified. Its torsion definition **is** stated: section rotation about the
  reference axis, positive toward stall (Fig. 19, PDF p. 12). The first torsional frequency is
  4.072 Hz (Table 3, PDF p. 8), near our 4.000 Hz; see P0 for limits of comparison.
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

### [x] P0. Provenance and pre-registration (docs only; no solver run; one session)

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
Closes when: items 1-3 recorded with cites (a "not stated" is a valid answer) and the P1
pre-registration committed; a comment on the issue lists them.

#### P0 record: provenance and external model (PDF pages are PDF page indices)

| Question | Evidence and bounded answer |
|---|---|
| Which deck/version? | The vendored `tests/IEA15MW/validation papers/github-IEA-15-240-RWT/README.md:3` calls the repository **v1.1**; its `ReleaseNotes.md:3-25` extends through **v1.1.6**, with the SONATA/ANBA4 statement specifically under **v1.1.4** (line 25). `tests/IEA-15-240-RWT.yaml:1,333-339` has no revision or solver label; git first added it at `6e5722b` (2026-04-26, “Tests update”). Its root K11/K16 match the vendored WindIO `WT_Ontology/IEA-15-240-RWT.yaml:599-606,611` and the BeamDyn blade deck at `OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn_blade.dat:12-17`. **Exact upstream release/tag and generating solver for this particular copy: not established.** Do not equate a repository README version with a pinned deck revision. |
| Does the report identify its solver? | Gaertner et al., `.sources/papers/75698.pdf`, **PDF p. 23 (printed p. 13), §2.2**: section stiffness was computed with PreComp **and VABS** in NREL's tool chain, or BECAS in DTU's; Figs. 2-10 and 2-11 (both PDF p. 26) show resulting blade properties including shear centre. That is report-era provenance, **not** proof of this test YAML's export path. |
| Was non-zero K46/K56 lost? | The vendored WindIO `elastic_properties.stiffness_matrix.K46/K56` contains 26 exact zeros (`WT_Ontology/IEA-15-240-RWT.yaml:618,620`). The test WindIO-style `elastic_properties_mb.six_x_six.stiff_matrix` has 26 upper-triangular 6×6 rows (`tests/IEA-15-240-RWT.yaml:333-365`); K46/K56 were parsed as zero at all stations, while WindIO K16/K26 are non-zero (`WT_Ontology/.../IEA-15-240-RWT.yaml:605,612`). The BeamDyn file stores a full 6×6 matrix, including zero [4,6]/[5,6] at its root (`OpenFAST/.../IEA-15-240-RWT_BeamDyn_blade.dat:12-17`). `ReleaseNotes.md:25` says properties consistently originate from SONATA/ANBA4 as of v1.1.4; it neither includes raw ANBA4 matrices nor specifies an export operation. **Whether the raw run had non-zero K46/K56, and whether an export discarded them: not stated/unverified.** Zeros in a delivered matrix cannot establish loss. |
| Zhou `C*` and nonlinearity | Zhou et al., `.sources/papers/Unsteady aeroelastic performance of the 15 MW floating offshore wind turbine under surge condition.pdf`, **PDF p. 5, §3.2, Eq. (23)**: `(F*, T*) = C* (γ*, κ*)`, with forces, strains and `C*` resolved in the **deformed-section frame Bi**. The model is geometrically exact and handles nonlinear large deflection (§3.2, PDF p. 5); §2, PDF p. 3 points to the IEA report/repository for turbine data, but **does not name a generator/version of `C*` or say whether the particular stiffness includes K46/K56**. Its reference axis is described kinematically (§3.2, Fig. 5, PDF p. 5); its **chordwise position is not given**. |
| Zhou twist and operating point | Zhou **PDF p. 12, §4.2/Fig. 19** explicitly defines torsional deflection as rotation of the airfoil section about the reference axis, **positive toward stall** (negative is feather). **PDF p. 10, Table 4**: rated, fixed-platform mean tip torsion **−3.60° (Present LL-FVW+GEBT)** vs **−3.77° (ALM-GEBT)**; do not silently replace their signed convention with ours. Rated wind **10.59 m/s**, speed **7.55 rpm** (§3.3, PDF p. 6); Table 1, PDF p. 3. The paper gives no explicit **rated blade pitch for Table 4** nor a numerical reference-axis position. **PDF p. 8, Table 3** lists their first torsional frequency **4.072 Hz**. These are comparisons, not proof of transferable section stiffness (#14). |
| ALM-GEBT linear-vs-nonlinear tip-twist factor | Chuang et al., `.sources/papers/jmse-13-00116.pdf`, **PDF p. 15, §4.2**, compare nonlinear GALM/GEBT with *linear Euler–Bernoulli* EALM under otherwise identical computational settings; **PDF p. 16, §4.2.1** explicitly says EALM does **not** account for blade torsional deformation. GALM tip torsion is **−4.1° to −3.3°** in Fig. 18 (PDF p. 18), mean **−3.768°** in Table 5 (PDF p. 14), rated **10.59 m/s, 7.55 rpm** (PDF p. 9). **No comparable linear tip twist or tip-twist ratio is reported; the external P1 geometric-twist factor is unavailable.** Their flapwise tips **15.920 m linear vs 14.10 m nonlinear** (PDF p. 17, §4.2.2; ratio 0.886) reflect changed aerodynamic loads *and* a different beam class, **not** a fixed-load structural torsion factor. |

#### P1 pre-registration — frozen before modifying the diagnostic tool

**Scope and inputs.** One `sequana_cpu_dev` job, maximum one harness-fix rerun; official
`tests/IEA-15-240-RWT.yaml` at **0.25 m** mesh (never the existing radii script's 1.0 m
fixture), and the associated vendored `OpenFAST/IEA-15-240-RWT/IEA-15-240-RWT_BeamDyn.dat`
and `..._BeamDyn_blade.dat`. Use `docs/validation_data/zhou_2025_fig11_loads.csv` and
`zhou_2025_fig10_aoa.csv` without changing their values; obtain the shell nodal force
from the unmodified `tools/diagnose_zhou_loads_reverse.py` formulation at 0.25 m
(`:66-96,99,192-270` for deck, mesh, reconstructed moment and projection). No gravity, rotation or deck-property edits. Do not use S-7 (#20) as a
control or infer a pitch that Zhou did not state. All results refer to the same blade and load.

**Future CLI contract, not a command that works today.** P1 must extend
`tools/diagnose_radii_beamdyn.py` with `--twist-p1` and `--output-dir` (retain its existing
`--deck`, `--beamdyn`, `--work`); add `RDzr` to the BeamDyn OutList and convert its
Wiener–Milenkovic rotations to section angles. Its current `read_beamdyn` reads translations
only (`:179-194`) and its current `shell_state` uses 1.0 m (`:71-108`): running it as-is is
**not P1**. Implement and inspect the interface before the one job; no numerical result is
registered until C1 passes. P1's job command, from the `rotor-frame` worktree root, is:

```bash
# BEAMDYN_DRIVER must point to the existing v5 beamdyn_driver executable.
# In the sequana_cpu_dev job (one submission; log from the scheduler):
export PYTHONPATH="$PWD/src"
: "${BEAMDYN_DRIVER:?set the v5 beamdyn_driver path before submission}"
scripts/aeroenv.sh python -W ignore tools/diagnose_radii_beamdyn.py \
  --twist-p1 --deck 'tests/IEA15MW/validation papers/github-IEA-15-240-RWT/OpenFAST/IEA-15-240-RWT' \
  --beamdyn "$BEAMDYN_DRIVER" --work logs/check/twist_p1_beamdyn \
  --output-dir logs/check/twist_p1
```

Submit that command using a single `sbatch` script on `sequana_cpu_dev` with a 20-minute
limit, and scheduler stdout/stderr `logs/check/twist_p1_%j.out` and
`logs/check/twist_p1_%j.err` (`%j` is the job ID); the job script is
`logs/check/twist_p1.srm` (to be written in P1, not executed in P0). Save machine-readable
`logs/check/twist_p1/results.csv` (one row per case and station with signed torsion, force,
moment and units), `logs/check/twist_p1/summary.json` (C1 residuals, four signed factors,
D2 and branch decision, source revision/input hashes), and each generated driver/deck in
`logs/check/twist_p1_beamdyn/<case>/` (including `bd_driver.inp` and `bd_driver.out`). The
script must fail non-zero before factor calculation if C1 fails; preserve logs and commit
**a small summary table** into this feature document in P1 (ignored logs are not the record).

**Order and decision, fixed now.** Run C1 first: apply one pure tip couple, independently
integrate the section twist using the **same deck's full 6×6 compliance**, not `T/K66`:
transform the constant global applied moment into each deck section's local frame, compute
`(K(z)^{-1} [0, 0, 0, Mx, My, Mz])_6` and integrate over the section arclength in radians;
apply the documented BeamDyn frame/orientation conversion before comparing angles at
50/75/90 % and tip. The nonzero K16/K26 mean `T/K66` is *not* an independent 1 %
reference. Each station must be within **1 %** with consistent signed frames, units and
station mappings. If geometry/section-frame transport makes that 1-D reference
incommensurate with the static BeamDyn solution, report a **failed control**; do not
reinterpret a factor. A C1 failure halts the job and P1 (one harness-fix rerun allowed
only in a *new* session).
Then, for each of `shell_load`, `zhou_formula`, `zhou_formula+Tp_arm`,
`zhou_formula+twist_projection`, `zhou_formula+shell_AC_arm`, evaluate **both** the full
load (`nl`) and 1 % load divided by 0.01 (`lin`) using the same undeformed strips. Fix
conversion to nose-down-positive for the comparison, but keep signed raw BeamDyn angles
in the CSV. Record four factors at the rigid 75 % ring and tip using the definitions below;
negative or near-zero denominator is an explicit **unattributed/undefined** result, not a
reason to switch to absolute values or widen a band. Attribute a load-path arm only if
its toggled sum reproduces the shell-load factor within **5 %**; structural ratio **<= 1.05**
skips P2, **> 1.05** enters P2; beam-class difference **> 5 %** retires the decoupled beam
as a reference. For D2 report signed `BD_nl(zhou_formula)` at the tip beside Zhou's
**−3.60° toward feather**, and the numerical difference after a documented frame
conversion; a relative difference <= **5 %** is only a diagnostic coincidence, **not**
reproduction or closure (#14 applies in either case). Report `BD_nl/BD_lin`
as an internal fixed-load geometric factor; **no comparable ALM-GEBT tip-twist factor exists**
(P0), so no external pass claim. Report tip-couple shell/BeamDyn for #20, not as a gate.
No other hypotheses, jobs or thresholds are admitted inside P1.

**P0 closure (2026-10-13).** Evidence and P1 pre-registration: local work-unit commit
`21cb8d9a1597` on `fix/bem-participant-rotor-frame` (not pushed). Issue comment:
https://github.com/efirvida/AeroElast/issues/42#issuecomment-6148712619.
Docs-only, no meaningful RED or solver run; `scripts/check.sh quick` exited 0
(66 passed, 1 xfailed, 1 documented known-red; 697 warnings), and `git diff --check`
exited 0. This smoke gate says nothing about twist physics. Next entry: P1 only when the
maintainer schedules it under #33 after its B3/L3 audit; use this frozen pre-registration
before modifying the BeamDyn diagnostic.

### [ ] P1. Attribution run with BeamDyn (one job, ~10 min, `sequana_cpu_dev`; one session)

Entry: P0 record. Tool: extend `tools/diagnose_radii_beamdyn.py` (add `RDzr` to the OutList,
Wiener-Milenkovic to angle; strip loads from Np/Tp/Mp) and reuse the Zhou-load shell solve.
Runs (all on the official deck, K16/K26 kept, no gravity, no rotation):
- `BD(shell_load)`: the shell's own nodal forces binned to strips about the reference axis.
- `BD(zhou_formula)`: Np, Tp, Mp per strip with the beam's arm.
- `BD(zhou_formula + Tp arm)`, `BD(zhou_formula + twist projection)`,
  `BD(zhou_formula, arm about the AC point used by the shell)`: one term toggled each.
- Each at full load (`nl`) and at 1 % load rescaled x100 (`lin`).
Control C1 (must pass before any factor is read): BeamDyn under a pure tip couple equals an
independent full-6×6-compliance integration of the same deck (section-frame moment transport,
not `T/K66`) within 1 % at 50/75/90 % and tip. It tests frames, units and the
Wiener–Milenkovic conversion, without the shell. Fail -> report, stop; fixing the harness is
the next session.
Factors at the rigid 75 % ring and the tip (same estimator as `_ring_kinematics` on the shell):
| factor | definition | triage |
|---|---|---|
| structural | `omega_shell / BD_lin(shell_load)` | <= 1.05: D1 is load path, P2 skipped. > 1.05: P2 |
| load path | `BD_lin(shell_load) / BD_lin(zhou_formula)` | the toggled runs name the arm term(s) whose sum reproduces it within 5 %; otherwise "unattributed", recorded, not chased |
| geometric | `BD_nl / BD_lin` | internal fixed-load factor only; no comparable ALM-GEBT tip-twist number (P0) |
| beam class | `BD_lin(zhou_formula) / decoupled beam` | reported; > 5 % means K16/K26 matter and the decoupled beam is retired as a reference |
D2 diagnostic: compare `BD_nl(zhou_formula)` with Zhou's signed −3.60° (5 % is a
**triage band only**), but label both agreement and disagreement **non-transferable**:
our static dead-load, non-rotating, no-gravity beam with a reconstructed Mp and unstated
rated pitch cannot reproduce Zhou's mean rotating, gravity-affected, two-way aeroelastic
solution. Keep #14's non-transferability verdict; **do not close D2 on this numerical
comparison** or treat a coincident value as an independent reference.
Also reported, not gated: the tip-couple shell/BeamDyn ratio, as evidence for #20 (not closing it).
Closes when: C1 passed and the four factors plus the D2 diagnostic (explicitly
non-transferable per #14) are committed and commented.

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
