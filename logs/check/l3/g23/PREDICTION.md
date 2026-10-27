# Group 23 prediction (written before the run, 2026-10-27, branch l3/g23-bernardi-p20)

Run: whole `tests/validation/blade/test_blade_iea15mw_validation.py`, es 0.25 mid-surface, ccx 2.20 binary.
Signed error = shell / reference - 1. Basis: S-4 parked vs p20 (1F +3.42, 1E +0.20, 2F +1.80 %), S-1 GJ -10.7 %
(shell torsion soft, but warping restraint stiffens it), Bernardi comparison of 2026-10-22 (worst 11.4 % at mode 7),
CCX parity of 2026-10-24 (0.01-0.12 % on the *best* five of ten).

## BeamDyn p20 rows (rows `beamd.3`-`beamd.7`, bound 5 %)
| row | mode | p20 Hz | predicted signed | range I would accept as "predicted" | red if |
|---|---|---|---|---|---|
| .3 | 2E | 2.136560 | +1.5 % | -1 .. +4 % | > 5 % |
| .4 | 3F | 2.935456 | +1.0 % | -2 .. +4 % | > 5 % |
| .5 | 1T | 4.093666 | -3 % | -7 .. +3 % | > 5 % (GJ -10.7 % says low; warping says up) |
| .6 | 3E | 4.379791 | +2 % | -2 .. +8 % | > 5 % (most likely red: Bernardi gap 11.4 % was at mode 7) |
| .7 | 4F | 4.844770 | +2 % | -3 .. +6 % | > 5 % |
Expectation: 0-2 of five red. A red row is reported as a finding (no xfail, no band change).

## CCX rows (`ccx.0`-`ccx.4`, CCX modes 1-5 now fixed; bound 5 %)
Both are shell FE on the same mesh: predicted gap 0.01-0.2 % per mode, all four-fold inside the bound; the old 0.01-0.12 %
were the five best of ten, the true modes 1-5 are expected only slightly larger (worst <= 0.3 %). Refuted if any mode 1-5 > 1 %.
Static row `static_tip`: 0.3-0.5 % (0.37 % on 2026-10-24).
Mode identity: shell modes 1-5 = 1F, 1E, 2F, 2E, 3F (no local or torsion mode below 3F); the Hungarian pairing against the
eight p20 modes then gives the same partners. Refuted if the p20 pairing returns a non-monotone partner list.

## Other tests in the file (unchanged, no row)
article / NREL / NuMAD bounds: same outcomes as 2026-10-22 (2F vs NuMAD xfail at ~12 %). Mass and DLC 1.4 tests: unchanged.
The log prints `ccx binary: ... This is Version 2.20`.

## Job
SLURM 11625842 (sequana_cpu_bigmem: 0 of the user's jobs running there vs 3 on sequana_cpu), `--ntasks=8 --time=01:30:00`,
log `logs/check/l3/g23/job-11625842.out`, script `logs/check/l3/g23/run.sbatch`. Submitted 2026-10-27, not awaited.

## After the job (parent)
```
export PYTHONPATH=$PWD/src; E=/scratch/leahk/eduardo.donestevez/fem-shell/scripts/aeroenv.sh
$E python tools/validation_matrix.py regression --group 23 --write    # run on a compute node (it re-runs the file, ~25 min)
$E python tools/validation_matrix.py extract --group 23 --write       # only if the job log changed a row id (none expected)
# then replace the 11 `undecided` entries of group 23 in docs/validation/contract.yaml by verdicts (K3: > 5 % by >= 1 pp rejected,
# within 1 pp of 5 % without es 0.125 undecided), set audited_rev, and gate:
$E python tools/validation_matrix.py check && $E python tools/validation_matrix.py contract --level L3 --rows | grep blade_iea15mw
```
`regression --write` re-executes the test file, so it may be submitted instead of reading the job log (the job 11625842
only gives the numbers to compare with this prediction and to write the verdicts).
