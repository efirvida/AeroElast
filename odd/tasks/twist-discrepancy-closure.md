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

#### P1 re-pre-registration (2026-10-23), frozen before any tool edit or run

User (verbatim): "Re-preregistrar P1 en odd/tasks/twist-discrepancy-closure.md: convención
nose-up, factor geométrico con rotación, el control C1 primero, y T5.d absorbido. Hacerlo en un
commit propio, antes de cualquier corrida."

This block **supersedes the P0 pre-registration above and the P1 box below where they differ**.
Everything it does not mention stays as frozen in P0: the official deck, the 0.25 m mesh, the
unchanged Zhou CSVs, the K16/K26 kept, no gravity, no deck edits, S-7 not a control, and the 5 %
triage meaning of the bands.

**Why the re-registration.** `4ee98ab` restored the aerodynamic torque sense against OpenFAST v5,
and `7189314` carried the convention into the anchor beam, the tools and the docs
(`odd/tasks/l3-physical-validation.md` L34-L36). Every number in "Why we kept iterating" was taken
with the inverted torque sense and the old nose-down convention, so D1 and D3 are restated, not
compared:
- rated production tip twist at `4ee98ab`, es 0.25: **-3.13°** (OpenFAST v5 AeroDyn+BeamDyn -3.40°,
  Zhou -3.60°; all nose-up positive);
- anchor beam under the production load (SLURM 11622758): beam -1.7765°, shell ring -4.0462°
  (ratio 2.28), mean θz -2.9720°;
- shell vs BeamDyn full 6×6 under the same rated nodal load, linear: 0.1125 vs 0.1228 rad (-8.4 %;
  L33, load taken before the sign fix, sign-independent ratio);
- T5.d first reading (L36): a static, non-rotating BeamDyn run flaps 19.0 m against OpenFAST's
  14.08 m with rotation, so its nonlinear/linear twist factor is not the rotor's.

**R1. Convention: nose-up positive.** Every angle is a right-hand rotation about the mesh +z (span)
axis, with the LE at +x and the suction side downwind at +y. +z is nose-up, toward stall (L34 mesh
probe), which is Zhou's "positive toward stall" (Fig. 19, PDF p. 12). Zhou's -3.60° and
OpenFAST's -3.40° therefore enter **unchanged**. P0's "conversion to nose-down-positive" is
withdrawn. BeamDyn's blade frame (x_BD = mesh +y, y_BD = -mesh x, z_BD = mesh z) is a proper
rotation about z, so `RDzr` keeps the sense of the mesh θz. C1 checks this numerically and does
not assume it. The CSV keeps the raw Wiener-Milenkovic components `RDxr/RDyr/RDzr`.

**R2. Estimators.** Primary BeamDyn section angle: φz = 4·atan(c_z/4) of the `RDzr` component,
the definition behind the OpenFAST -3.40° (`logs/check/l3-bd-openfast/analyze2.py:16-17`).
Secondary, reported and not gated: the twist about the deformed tangent, from the swing-twist
split of the full rotation (R = R_twist(t)·R_swing). The two differ at second order in the flap
slope, so their gap is printed beside every nonlinear angle. Shell: `_ring_kinematics` ring
rotation (primary) and mean nodal θz (secondary, L37 N2). Output format `OutFmt "ES16.8E3"`
in every generated BeamDyn file: the default lost digits at small scale (L36 item 2).

**R3. Controls first, in their own job.** Job J1 runs only the controls. Job J2 (the case matrix)
is submitted only after J1 has passed and its numbers are committed in this file. A failed control
stops P1, as P0 fixed (rule 3). The harness fix is the next session.
- **C1 (frames, units, WM conversion).** As P0: a pure tip couple T_z = 1.0e4 N·m (predicted tip
  twist 1e-3 to 1e-2 rad, the linear range). Reference: the same deck's full 6×6 compliance,
  `(K(z)^-1 [0,0,0,M_local])_6`, with the constant global couple transported into each section's
  local frame and integrated over the arc in radians. Pass: every one of 50 %, 75 %, 90 % and tip
  within **1 %** and with the same sign. The BeamDyn run is the static solve
  (`DynamicSolve False`).
- **C2 (rotation; new, because the geometric factor now carries rotation).** BeamDyn driver with
  `DynamicSolve True`, `QuasiStaticInit True`, root angular velocity Ω = 7.56 rpm
  (0.791681 rad/s) about the rotor axis x_BD (`RootVel(4)`), blade root at the hub radius
  `GlbPos = (0, 0, 3.97 m)` (ElastoDyn `HubRad`, `bd_mu30/ed.dat:47`), no precone and no tilt.
  The read state is the first output row, the quasi-static initialisation; its time is printed.
  - C2a, no applied load: root axial force against ∫ m(s) Ω² ρ(s) ds, with m(s) the deck's mass
    per length and ρ(s) the distance from the rotor axis of the deck key-point curve. Pass: within
    **1 %**.
  - C2b, Ω = 0 with the P-shell load at scale 1: the quasi-static path against the static solve.
    Pass: tip φz and tip flap within **0.1 %**.
  If C2 fails and C1 passes, J2 runs without the rotating column, and the rotating geometric
  factor is reported as "not measured, C2 failed". No other substitute is admitted.

**R4. T5.d absorbed.** The decomposition planned in `l3-physical-validation.md` L36 step 3 moves
into P1 as the P load family. It is no longer a separate L3 task. Its scratch directory
`logs/check/l3-twist-residual/` is not created. P1's own directories are used.

**R5. Cases.** All cases use the same 0.25 m strips and the same binning about the reference axis
(the existing `shell_state` construction at 0.25 m, never its 1.0 m default).

| family | case | load | origin |
|---|---|---|---|
| Z (D1, D2) | `Z_shell` | shell nodal force from the Zhou loads, binned | P0 |
| Z | `Z_formula`, `Z_formula+Tp_arm`, `Z_formula+twist_projection`, `Z_formula+shell_AC_arm` | Np, Tp, Mp per strip, one term toggled each | P0 |
| P (T5.d) | `P_shell` | production rated rigid-blade nodal load at HEAD (`ForceProjector` after `4ee98ab`, the de-loading fixture), binned | T5.d (ii) |
| P | `P_BlSwpAC` | the same Np, Tp, Mp with the AC lever on the AeroDyn `BlSwpAC` line instead of the mesh ring | T5.d (iii) |

Each case runs in three columns: `lin` (scale 1e-2, divided by 1e-2; P0's 1 %), `nl` (scale 1,
static, no rotation) and `nlrot` (scale 1, C2 configuration, T5.d (iv)). Shell side: one linear
one-way solve per family, `omega_shell(Z)` (P0 route) and `omega_shell(P)` (T5.d (i)). The
production structural solve is linear, so the shell has no nl column.

**R6. Factors** at the rigid 75 % ring and the tip, with the primary estimators (R2):

| factor | definition | triage |
|---|---|---|
| structural | `omega_shell(F) / BD_lin(F_shell)`, F = Z and P | <= 1.05 on Z: D1 is load path, P2 skipped. > 1.05: P2 |
| load path | `BD_lin(Z_shell) / BD_lin(Z_formula)` | as P0: toggles name the arm term(s) within 5 %, else "unattributed" |
| AC route | `BD_lin(P_BlSwpAC) / BD_lin(P_shell)` | reported; L36 predicts about 1.16 from the torque alone |
| geometric | **`BD_nlrot / BD_lin`** per case | replaces P0's static `BD_nl / BD_lin` as the factor. The static ratio is still printed, as the T5.d over-flap control, beside the `nlrot` and `nl` tip flaps |
| beam class | `BD_lin(Z_formula) / decoupled beam` | as P0 |

Negative or near-zero denominators are reported as undefined (P0 rule).

**R7. Readings that are not verdicts.**
- D2: `BD_nlrot(Z_formula)` and `BD_nl(Z_formula)` at the tip, beside Zhou's -3.60°. These are
  non-transferable per #14, and P0's limits apply unchanged.
- Rated residual (T5.d): shell -3.13° vs OpenFAST -3.40°, R = 1.086. The P cases use a rigid-blade
  one-way load. OpenFAST's load is two-way (thrust -13.2 %), so the comparison is triage only.
  Classification, fixed now:
  - (a) structural, if `1/structural(P)` alone carries at least half of ln R with the same sign;
  - (b) AC route, if the AC factor alone does;
  - (c) geometric with rotation, if `geometric(P_BlSwpAC)` alone does;
  - (d) unattributed, if none does or two disagree in sign.
  `BD_nlrot(P_BlSwpAC)` is also printed beside OpenFAST's -3.40° (tip twist) and 14.08 m (tip
  flap).
- The tip-couple shell/BeamDyn ratio, for #20 (not a gate).

**R8. Commands, outputs, budget.** The CLI contract of P0 stays (`--twist-p1`, `--output-dir`,
`--deck`, `--beamdyn`, `--work`), with a required stage, `--twist-p1 controls|matrix`. Logic
that makes `diagnose_radii_beamdyn.py` unreadable may live in a sibling module the flag
dispatches to. The job script `logs/check/twist_p1.srm <stage>` runs on `sequana_cpu_dev`
(20 min, `--ntasks=4`), with stdout/stderr in `logs/check/twist_p1_%j.{out,err}` and
`BEAMDYN_DRIVER=$SCRATCH/conda-envs/openfast/bin/beamdyn_driver`. Outputs:
- `logs/check/twist_p1/controls.json` (C1 and C2 residuals, revision, input hashes);
- `results.csv` (one row per case, column and station: raw WM components, φz, the tangent twist
  and flap, in SI units);
- `summary.json` (factors, classification, D2);
- each BeamDyn deck and output under `logs/check/twist_p1_beamdyn/<case>/`.

The controls stage exits non-zero on a failed control. The matrix stage refuses to run unless
`controls.json` records C1 passed. The pure functions (WM to angle, swing-twist split, frame map,
6×6 compliance integration) get unit tests before the harness is used. Budget, amending rule 7
for P1 only because C1 now has its own job: J1 + J2 + at most **one** harness-fix rerun, in a new
session, in total.

### [ ] P1. Attribution run with BeamDyn (one job, ~10 min, `sequana_cpu_dev`; one session)

> Superseded where it differs by the 2026-10-23 re-pre-registration above: nose-up, `nlrot` geometric
> factor, the controls in their own job, and the P (T5.d) family.

#### P1 J1 result (2026-10-23): C1 failed, P1 stopped

Harness `5978b73` (`tools/twist_p1.py`, 24 unit tests; `--twist-p1 controls`). J1 = SLURM
**11623050** on `sequana_cpu_dev`, revision 5978b73, 3 min 43 s, exit 3 (C1 failure), log
`logs/check/twist_p1_11623050.out`, `logs/check/twist_p1/controls.json`. A development run of the
same stage on the login node gave identical numbers (`logs/check/twist_p1-loginnode-dev/`). It is
not the record.

| control | station | BeamDyn φz (rad) | reference (rad) | rel. | band | result |
|---|---|---|---|---|---|---|
| C1, T_z = 1e4 N·m | 50 % | 9.37e-4 | 7.93e-4 | +18.1 % | 1 % | FAIL |
| | 75 % | 3.864e-3 | 3.739e-3 | +3.4 % | 1 % | FAIL |
| | 90 % | 9.339e-3 | 1.050e-2 | -11.0 % | 1 % | FAIL |
| | tip | 4.491e-2 | 5.904e-2 | -23.9 % | 1 % | FAIL |
| C2a, Ω = 7.56 rpm, no load | root | `RootFzr` 1.31224e6 N | 1.31404e6 N | -0.14 % | 1 % | PASS |
| C2b, Ω = 0, P_shell scale 1 | tip | QS = static to 1e-7 (φz -0.077095 rad, flap 19.228 m) | | | 0.1 % | PASS |

The signs agree at every station. The tip tangent twist (secondary) is 0.04508 rad. The
reference is converged: with a quarter of the subdivisions it reads 0.05917 at the tip (0.2 %).
The non-gating prismatic diagnostic, with the root section at every station, matches BeamDyn
within 0.15 %, and the root Mz reaction is 9999.99 N·m.

**Verdict: C1 FAILED.** Per R3 and rule 3, J2 was not submitted and no factor was computed. No
cause is adopted in this session. Observations for the next session, all unverified:
- K66 falls about 190× over the last 5 % of span (1.34e7 to 7e4). In the reference, 82 % of
  the tip twist accrues outboard of 90 % span.
- The reference interpolates K linearly between deck stations. BeamDyn's own interpolation and
  quadrature (one member, the deck's element order, trapezoidal with refine 2) were not checked.
  Linear-in-K and linear-in-compliance differ strongly over such a drop.
- The OpenFAST -3.40° arbiter ran this same deck discretisation, so whatever resolves C1 also
  bears on that number.

**Next entry (new session, the one harness-fix rerun of R8).** Decide what C1 compares before
changing anything: the deck's continuum (refine BeamDyn until it converges) or BeamDyn's discrete
interpolation (match the reference to it). Pre-register that choice, then rerun J1. J2 follows
only after C1 passes.

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
