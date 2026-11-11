# T5 — group 64 `v2_blade_static_prescribed` (v2 counterpart of group 39)

Spec coverage (odd/tasks/v2-validation-store.md): S1, S2, S4, S5 (S3: prepared and measured once, cheap enough: 5 min 50 s on SLURM).
Test: `tests/validation/blade/v2/test_v2_iea15mw_s2_static_prescribed.py` (8 tests). Run: SLURM 11641783, `sequana_cpu_dev`, exit 0, 8 passed,
log `logs/check/v2_s2/run-11641783.out`. Tree 0da03be + uncommitted worktree (v2 mesher edits); es 0.25, `skin_reference="mid"`, no strip collapse.

## Same as v1
Official blade only; LC1 gravity, LC2/LC3 1 MN tip, LC4 5 kN/m (60 slices); loads through `run_s2_static_cases.run_shell_case`
(tip `equal` for the tip tests, `arc` for the station tests; LC4 `equal` for the tip test, `area` for the station test); BeamDyn p60 references
(tip uy -2.231483 / 82.273278 / 9.062631, LC3 ux 34.469441; station nodes `tests/IEA15MW/bd_static`, scale gate 1e-3); metrics (tip node = argmax z;
S12 arc ring mean interpolated at the predeclared BeamDyn node, inner of the two nodes >= 2 tip chords inboard); tolerance K3 = 5 % everywhere (`_TOL`);
asserted: tip LC1, LC4 and station LC1, LC4; printed: tip LC2, LC3 ("K11.3"), station LC2, LC3 (S13). Test ids mirror v1 (`test_dominant_component_matches_beam[...]`).

## What differs, and why
* Mesh: `build_v2_blade(0.25, skin_reference="mid")`.
* No `utd` parametrisation (v1 UTD is a disowned input-sensitivity record) and no in-house EB beam (printed only in v1; built from z-sliced v1 geometry).
* v2 rings are constant-eta rows, not constant-z planes (tip ring: 30 nodes spanning ~0.1 m in z; z.max = 117.005 m). Adapted in the test file:
  - tip load nodes = `TipNodes` (v1 `run_shell_case` takes `z >= zmax - 1e-6`, which is one node on v2) with v1 `ring_arc_weights`;
  - ring profile = eta rows (497 rows), weights = v1 rule (half of each element side with both corners on the row, skin and web), ring z = arc-weighted mean z.
    The raw eta of each node is recovered by matching raw to mesh elements by centroid (rcm renumbering reorders elements); bijection asserted;
  - tip chord = largest in-plane extent of `TipNodes` (v1: z-plane tip ring); z_tip = z.max() as in v1.
* Distributed loads (LC1, LC4) use the unmodified v1 z-slice code (`z.min..z.max`, 60 slices, `area` weights).

## Measured (v2, SLURM 11641783) next to v1
Tip node vs BeamDyn p60 (v1 = registry rows of group 39, run d5d0fd2, 2026-10-24, es 0.25 mid-surface):

| case | v2 shell | BeamDyn | v2 | v1 (registry) |
|---|---|---|---|---|
| LC1 uy | -2.190924 | -2.231483 | -1.82 % PASS | -5.26 % FAIL |
| LC2 uy (not asserted) | 90.829308 | 82.273278 | +10.40 % | -2.54 % |
| LC3 ux (not asserted) | 38.637278 (ring mean 38.794224, node vs ring -0.40 %) | 34.469441 | +12.09 % | +10.26 % |
| LC4 uy | 8.829434 | 9.062631 | -2.57 % PASS | -5.48 % FAIL |

Interior station (node 56, z 114.845 m, 4.35 tip chords; v1 has no station rows in the registry):

| case | v2 shell | BeamDyn | v2 signed |
|---|---|---|---|
| LC1 uy (asserted) | -2.124655 | -2.158325 | -1.56 % PASS |
| LC4 uy (asserted) | 8.446527 | 8.632019 | -2.15 % PASS |
| LC2 uy (printed) | 74.085445 | 73.611629 | +0.64 % |
| LC3 ux (printed) | 35.188603 | 31.903817 | +10.30 % |

Prototype (`logs/check/t16/fem_anchor2.py`, equal tip / slice loads) expected LC1 -1.6 %, LC4 -2.4 %, LC3 +8.8 %: LC1 matches, LC4 -2.15 % (area vs slice load), LC3 +10.3 % (arc tip load vs equal split).
The v1 prototype json (`logs/check/t16/fem-anchor-v1.json`, a different run, not the registry) shows station -5.1 / -5.2 % for LC1 / LC4 on the v1 mesh.

## Findings
* Nothing xfailed; both asserted tip cases and both asserted station cases pass at 5 % on v2 (v1 tip rows are red at -5.3 / -5.5 %).
* Not asserted, large: tip LC2 +10.4 % and LC3 +12.1 % (tip point loads without MPC, unlike v1 LC2 -2.5 %), station LC3 +10.3 %. These are printed as in v1; they are diagnostics of the load introduction, not a pass.
* Single mesh point (es 0.25); no convergence trend measured here (0.15 would be a second point).

## Store
`group.yaml` (uses `source_files`, as group 62), `rows.yaml` (extract --write generated the 4 station rows only: the tip tests compare against a dict lookup, which the extractor reads as `no_comparison`;
the 4 tip rows are hand-written in the v1 row shape, with `asserted`/`info_only` as v1, source line 209; measured blocks filled from the run above, run id `0da03be+worktree`).
No `residual.json`: v1 group 39 has no entry in `residual-patterns.json`. `check --group 64` on the store copy (`logs/check/v2_s2/store`) shows no error naming group 64 rows; the only differences from `check --group 39` are
4 contract.yaml entries of the v1 rows (artifact of loading group 64 only) and "group 64 is on no level" (parent to place it in levels.yaml next to group 39). The 297 baseline store errors are pre-existing.
v1 registry vs test code: station tests exist in the v1 file but have no v1 rows; otherwise no disagreement found (tolerance 0.05 everywhere).
