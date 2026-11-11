# 63 v2_blade_sectional_distribution (T4, v2 counterpart of group 38)

Spec coverage (odd/tasks/v2-validation-store.md): S1 (same tests as v1 on v2), S2 (own files only), S4 (metric,
independent reference, tolerance on every row; es 0.25; negative results reported), S5 (nothing outside the allowed
surfaces; the store was edited in a copy, `logs/check/v2_s1/store`). S3: prepared and measured once (SLURM 11641763,
190 s); no campaign.

## Same as v1 (`tests/test_iea15mw_s1_sectional.py`)
Reference: BeamDyn 6x6 deck (`ei_*_ec` elastic-centre EI, mass, GJ, EA). Metric: mean |rel err| over 0.2-0.85 of the tip z,
five tests. Tolerance: 5 % in every test (`_TOL`, read from the v1 TEST CODE). Extractions are the v1 code, unchanged:
tangent-frame curvature EI (`ei_from_static_response(frame="tangent")`), lumped-mass partition, static torsion GJ with
the end-ring wall shear flow, direct CLT EA (`tools/run_s1_sectional.ea_gj_from_static_response`). Skin on the
mid-surface. 49 evaluation stations (bay centres of the 50-point reference-axis z grid, identical in the v1 and v2 yaml).

## What differs, and why
* Mesh: `build_v2_blade(0.25, skin_reference="mid")` (27 960 nodes, 29 143 elements, 497 eta rings, mass 67.106 t).
* Station membership: v1 reads element sets `{cc}_{ss}_...`, absent in v2. `tests/support/v2_sectional.py`
  (`V2SectionalExtractor`, subclass of the v1 `SectionalExtractor`) overrides ONLY `station_membership`: station = the
  nodes of the elements between ring r and r+1 (one bay, shared rings between neighbours, root ring = `RootNodes`).
  The ring of each node is recovered by exact coordinate match against the builder's raw arrays (the mesh is rcm
  renumbered). EA/GJ reuse the v1 function through a module-attribute swap (`v2_extractor_substituted`), not a copy.
* `extract()` computes only the mass path by default (the v1 test reads `sec.z` and `sec.mass_dens` only; the v1
  docstring flags the slice stiffness experimental); the other `SectionProperties` fields are NaN.
* v2 bays are 0.087-0.25 m (non-uniform: 0.087 m near z = 95.3), v1 bays are those of its own mesher. Stations are the
  bays nearest to the 49 grid centres (`select_stations`, as v1), so the sampled z differs by up to half a bay (~0.12 m).
* Assertions use `assert_residual_below` (canonical form the extractor reads) instead of a bare `assert err < 0.05`;
  same bound, and the printed residual feeds `regression`. Group 38 has no residual pattern (its rows were
  transcribed); this group has one (`residual.json`, same regex as group 57/58).
* Known inherited weakness (not changed): the v1 EA-CLT width uses `dz = z_st - z_(st-1)` of neighbouring bay mean
  z's, exact only for uniform bays. With v2's non-uniform bays this is an approximation next to the 0.087 m bays
  (stations at z = 94.4 and 96.6 sit there); effect on the mid-span mean not quantified.

## Measured (SLURM 11641763, 190.9 s, rev 0da03be; reproduced bit-for-bit by `regression --write`, SLURM 11641839)
mean |rel err| (signed) over 0.2-0.85, 31 stations, vs the v1 store rows (group 38, es 0.25 mid, f1c9e91):

| channel | v2 | v1 store | bound | v2 verdict |
|---|---|---|---|---|
| EI_edge | 8.65 % (-6.72) | 8.01 % (-1.89) | 5 % | FAIL |
| EI_flap | 4.53 % (+0.64) | 4.53 % (+3.87) | 5 % | pass |
| mass | 3.38 % (+2.08) | 3.54 % (+2.26) | 5 % | pass |
| GJ | 9.38 % (-9.38) | 10.74 % (-10.74) | 5 % | FAIL |
| EA (CLT) | 7.11 % (+6.96) | 6.68 % (+6.56) | 5 % | FAIL |

pytest: 3 failed, 2 passed (exit 1). Same pass/fail pattern as v1 and same order of magnitude; the three reds are the
v1 reds, left red (not xfailed, like the v1 file; no band widened). EI_edge signed error moved (-1.9 -> -6.7 %): the
v2 mesh/extraction differs, the abs error is similar. The v1 module docstring numbers (EI_edge 8.82 %, GJ 5.05 %,
EA 7.59 %, 2026-10-22) differ from the v1 registry rows above (8.01 / 10.74 / 6.68); the registry was not touched.
No registry/test-code tolerance disagreement for group 38 (all 5 %).

## Store
`check` on the copy (unscoped): 312 rows, 0 errors, 3 warnings (group 63 on no level: the parent places it).
`check --group 63` alone prints 300 cross-reference errors that also appear for any scoped group on the shared store
(e.g. `--group 34`); use the unscoped check. Rows carry `run: 0da03be`.

## Unit tests (`tests/support/test_v2_sectional.py`, 8 passed, 27 s)
Ring order and bay membership on an eta patch (0.30-0.32), mass-partition identity, total lumped mass vs an
independent layup integral (rel 2e-3), station selection, z-gap ring split on a synthetic stack, 49 bay centres,
bad-input rejection, substitution restore.

## Rerun
`sbatch -p sequana_cpu_dev --time=00:20:00 -J v2s1 --ntasks=1 logs/check/v2_s1/run_s1.sh` (~3.5 min).
