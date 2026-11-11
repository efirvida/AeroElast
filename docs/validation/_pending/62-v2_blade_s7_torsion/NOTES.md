# T3 — group 62 `v2_blade_s7_torsion` (v1 counterpart: group 34)

Spec coverage: S1, S2, S4, S5 (S3: nothing beyond the one measurement below was run).

## Result (RED — a finding, not a pass)
`tests/validation/blade/v2/test_v2_iea15mw_s7_torsion.py::test_v2_s7_global_torsion_ratio_within_band`
- v2 (es 0.25, mid-surface, 27 960 nodes, 40 slices, 24 in the 0.3–0.9 band): **ratio 1.07028**, error 7.0281 % vs bound 5 % → fails.
  Per-station ratio min 1.023, max 1.364; mid-span stations 1.02–1.04, rising toward the tip (1.057 at z 83, 1.093 at z 95).
- Shell twists MORE than the BeamDyn integral (shell GJ ~7 % below the deck), the opposite side to the v1 OML reading.
- Measured: SLURM 11641755 (pytest, exit 1, 44.6 s), 11641759 (per-station, exit 0), 11641829 (`regression --write`, exit 0). No Traceback in any log. Revision 0da03be (working tree dirty with other workers' edits).
- No xfail: v1 file is a plain red-allowed test too; left red.

## Same as v1
Load (tip ring Saint-Venant wall shear flow, 1000 N.m, `section_torque_wall_forces`), section rotation (`section_twists_deg`, 40 slices, all nodes), reference (`beam_twist_profile_deg`, BeamDyn K[5,5] integrated to z), metric (mean ratio over 0.3–0.9 of 117 m), tolerance `_RATIO_TOL = 0.05`, `assert_relative_error` text (v1 residual pattern for group 34 reused as `residual.json` key "62"), `ELEMENT_SIZE = 0.25`, `skin_reference="mid"`. Helpers imported from `tools/run_s7_torsion.py` unchanged.

## Differences
1. Mesh: `build_v2_blade(0.25, skin_reference="mid")`.
2. `run_torsion_case_v2` (copy of `run_torsion_case` in the test file): the tip ring comes from the `TipNodes` set (30 nodes). The v1 mask `z >= z.max()-1e-6` selects ONE node on the v2 mesh (tip z spans 116.995–117.005, root z min −0.017), so `ring_section` raised "needs at least 3 ring nodes" (first SLURM run). Nothing else changed.
3. Not repeated: v1's two machinery tests (`..._beam_reference_integrates_to_the_query_point`, `..._applied_moment_matches_request`) — mesher-independent, `non_validation_tests` in v1; group 62 therefore has no `non_validation_tests` block.

## v1 numbers vs registry disagreement (recorded, v1 not touched)
- v1 test docstring: ratio 1.0009 on mid at 7dfe78d (0.8980 OML). v1 registry row 34: 1.07171 (7.1712 %) at f1c9e91, 2026-10-27, and its justification text says "0.909 on the OML skin" while the test code uses `mid` and tol 0.05 (registry tol 0.05 agrees). I did NOT re-run v1; the v2 1.0703 is within 0.1 % of the registry's 1.0717, so the v2 mesh behaves like the current v1 HEAD, and the docstring's 1.0009 looks stale. Worth a v1 re-measure (parent decision).
- Tolerance: registry 0.05 == code `_RATIO_TOL = 0.05`; no disagreement. (The registry provenance_note still says ±30 % — stale prose.)

## Store work
Copy at `logs/check/v2_s7/store`; group entry appended, `extract --group 62 --write`, `set` for validates/justification/notes, `regression --group 62 --write` for the measurement. `check --group 62` reports 301 errors / 3 warnings, but the same store gives 300 errors for group 34: they are pre-existing contract/adjudication/gaps dangling refs, none mention group 62. Only new warning: "group 62 is on no level" (parent adds it to levels.yaml). The extractor's own `validates` was a mangled test name; replaced via `set`. Note `residual-patterns.json` in the copy was re-indented by json.dump — merge only `residual.json`'s "62" key.

## Files
`group.yaml`, `rows.yaml`, `residual.json`, this file; test `tests/validation/blade/v2/test_v2_iea15mw_s7_torsion.py`; scratch scripts in `logs/check/v2_s7/` (probe.py, measure.py, *.sh, logs, store/).
