# cases/l5_twoway_rated: L5 two-way, steady rated (thesis plan step 5, T5)

PREPARED, NOT RUN. Step 4 (one-way, SLURM 11741547) is not closed, so this case claims no numbers
(strict order, odd/tasks/thesis-validation-plan.md). Spec: plan S2, S3 (5 % ceiling, every number an independent reference).

Operating point: 10.59 m/s uniform (no shear), 7.56 rpm fixed, pitch 0, tilt -6, precone -4 (deck values), tower off, wake off,
gravity off, air 1.225 kg/m3. Mesh: gmsh chord-scaled N20, AR 1, mid skin (reference mesh; same as cases/l5_oneway_rated).

## Files

| file | what |
|---|---|
| `fluid/fluid.yaml`, `solid/solid.yaml`, `precice-config.xml` | full 100 s run (dt 0.01 s). Copy of the one-way case; the ONLY functional change is `bem.deformed_twist: true`, `bem.deformed_radius: true` (the BEM twist and radii follow the shell). tilt/precone explicit and labelled |
| `smoke/` | 0.5 s timing smoke, same files with `max-time 0.5` and `*_smoke` output folders |
| `openfast_deck_manifest.json` | copy of the generated OpenFAST deck manifest (sha256 of every source and generated file, options) |
| `../../tools/l5_twoway_openfast_deck.py` | reproducible generator of the OpenFAST decks (`--smoke S` runs each deck S seconds) |
| `../../logs/check/l5_twoway/run-twoway.sh` | sbatch launcher (sequana_cpu_long, `--ntasks=8`, TMPDIR under `$SCRATCH/tmp`); `full` or `smoke`. Output `$SCRATCH/l5_twoway/<mode>-<jobid>/` |
| `../../logs/check/l5_twoway/run-openfast.sh` | runs both OpenFAST decks (serial) |
| `$SCRATCH/l5_twoway/openfast/{flex,rigid,common}` | generated decks: ElastoDyn (rotor fixed 7.56 rpm) + BeamDyn + AeroDyn + InflowWind. `rigid` = same deck with BeamDyn 6x6 stiffness x1e3 (same geometry and mass), the flexible-minus-rigid twin |

## Our rigid twin (to be run, not yet configured)

Our delta needs a rigid-blade run of the same code: `deformed_*` false does not make the shell rigid (the structure still
deflects; that is the one-way case). Use the L4 rigid BEM (`BEMSolver`, tilt -6 / precone -4, deck-polars or YAML polars as in
l4-bem-rigid-closure.md) as our rigid reference, or a second case with a stiffened shell. Decision pending (open question).

## Quantities to extract (settled window = whole revolutions, T_rev = 7.937 s, tail of the run)

Per code (ours / OpenFAST BeamDyn+AeroDyn flex / OpenFAST rigid twin):
- mean rotor torque, power (torque x 0.79168 rad/s), thrust. OpenFAST: `RtAeroFxh`, `RtAeroMxh`, `RtAeroPwr` (AeroDyn, hub frame) in `flex/run.out`, `rigid/run.out`; ours: BEM log / `bem_fsi_results`.
- tip deflection and tip twist: OpenFAST `B1TipTDxr/yr/zr`, `B1TipRDzr` (BeamDyn, root frame) and `B1N0nn_RDzr` for spanwise twist; ours: the tip watch-point and the ring section rotation (`section_contour.ring_section_rotation`; the production BEM estimator has equal nodal weights, S8 pending).
- deltas: flexible-minus-rigid in torque, power and thrust, for each code; compare the two deltas (the thesis quantity, 5 %).
- gate: 5 % relative on mean torque/power/thrust and tip twist per code pair. Sign conventions of twist per odd/tasks/rated-twist-sign-convention.md and plan l5-one-way L10 (twist sense arbitration): report as measured.

## Secondary reference: Zhou et al. 2025 (Energy 336:138488)

Rigid 16.11 MW / 2.53 MN; flexible 14.76 MW / 2.20 MN (reductions 8.38 % / 13.04 %); tip (LL-FVW & GEBT) flap 13.86 m, edge -1.22 m, torsion -3.60 deg
(ALM-GEBT 14.10 / -1.27 / -3.77). Zhou ran cone and tilt zero, 7.55 rpm, pitch not stated: a Zhou comparison is a separate 0/0 case,
labelled, never mixed with the -6/-4 OpenFAST case (L4: -6/-4 vs 0/0 is +2.44 % torque). Source: odd/tasks/composite-bend-twist-verdict.md:688.
