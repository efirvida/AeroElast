# T8 group 67 — v2 counterpart of group 42 (V-05)

Covers S1, S2, S4, S5 (S3: prepared, nothing heavy run).

Test: `tests/validation/blade/v2/test_v2_iea15mw_v05_structural_properties.py` (3 tests).

## Same as v1
- Reference: integral of BeamDyn v1.1 deck M11 over the 117.149 m reference axis (`beamdyn_deck_mass`) = 66,996.9 kg.
- Metric: row-sum of assembled M over x-DOF rows; relative error. Tolerance 0.05 (same as the v1 test code; the v1 registry row also says 0.05, no disagreement).
- Guards: K/M assemble; K size multiple of 6; n_nodes >= 50 (declared `non_validation_tests`, as in group 42).
- es 0.25, mid-surface skin.

## Differs
- Mesh: `build_v2_blade(0.25, skin_reference="mid")` instead of `BladeMesh(...).generate("rcm")` + `build_rust_properties`.
- The mass assertion uses `assert_relative_error(..., tol=0.05, reference_name=..., kind="code")` (same `rel < tol` metric) because the extractor only reads that canonical call; a bare `assert rel_err < 0.05` yields no row.
- NOT repeated: v1 `TestBeamDynK55K66Consistency` (deck-vs-deck BeamDyn K55 vs ElastoDyn FlpStff, xfail record; no aeroelast code, mesher independent; already a non_validation_test in group 42).
- No `residual.json`: the v1 group 42 has no pattern in residual-patterns.json.

## Measured (login node, cached mesh, 12 s, exit 0, 3 passed)
| | mass kg | vs deck 66,996.9 | bound |
|---|---|---|---|
| v1 (registry, d5d0fd2) | 67,150.0 | +0.2285 % | 5 % |
| v1 (task note, 2026-10 mid mesh) | 67,105 | +0.25 % vs 66,934.6 | |
| **v2** | **67,105.6** | **+0.1623 %** | 5 % |

No failures, no xfail. The v2 mass is nearly identical to the v1 mid-surface mass (67,105 kg).

## Store
- Entry added to a copy (`logs/check/v2_v05/store`); `extract --group 67 --write` → 1 row; `check --group 67` reports 301 errors, all pre-existing store-wide ones (the baseline `check --group 42` on docs/validation gives 300; none mention group 67) plus the warning "group 67 is on no level" (the parent's levels.yaml merge).
- rows.yaml: the `measured` block and a reworded justification were set by hand on the generated row (`set` cannot index into the comparisons list). `run` is `<HEAD>+worktree` (uncommitted tree); restamp via `regression` after merge.
