# T2 — group 61 `v2_blade_v02_natural_frequencies` (counterpart of group 32)

Covers S1, S2, S4, S5 (S3: prepared, one cheap SLURM run only).

## Same as v1
Test body, references (NREL/TP-5000-75698 Table 5-2: 0.5585 / 0.6406 / 1.659 Hz), metric (closest of 10 modes), tolerances from the v1 TEST CODE
(1F 5 %, 1E 5 %, 2F xfail when err >= 5 %), modal solver, clamp of `RootNodes`, es 0.25, `skin_reference="mid"`. Residual pattern copied from group 32.

## Differs
Mesh/properties from `build_v2_blade(0.25, skin_reference="mid")` instead of `BladeMesh(...).generate("rcm")` + `build_rust_properties`; the
`iea_blade_yaml` fixture is gone (builder reads the windIO v2 YAML). Class renamed `TestV2BladeNaturalFrequencies`. Nothing else needed from the builder.

## Measured (SLURM 11641748, sequana_cpu_dev, exit 0, 104 s test time; no Traceback; log logs/check/v2_v02/pytest.log)
| mode | ref Hz | v2 Hz | v2 err | v1 (registry row, 19ff1a3) Hz | v1 err | bound |
|---|---|---|---|---|---|---|
| 1st flap | 0.5585 | 0.514244 | -7.92 % | 0.518922 | -7.09 % | 5 % |
| 1st edge | 0.6406 | 0.684892 | +6.91 % | 0.697231 | +8.84 % | 5 % |
| 2nd flap (closest mode to 1.659) | 1.659 | 1.4691 | -11.4 % | ~1.554 (es 0.25, docstring) | -6.3 % | xfail, 5 % |

Result: pytest `2 failed, 1 xfailed` — **v2 misses both 5 % bands** (as v1 does; v1 is also red by design on both). Not widened. v2 is
slightly further on 1F and 2F (softer), closer on 1E. The 2F gap (-11.4 %) is larger than v1's -6.3 %: unexplained, worth a follow-up
(possible mode-identification effect: closest-mode pick, or v2 mesh / OML->mid inset differences). Not investigated here.

## Disagreements v1 registry vs v1 code
None on tolerances (registry rows already 5 %/5 %). The v1 group-32 `provenance_note` still says the 1E bound is 10 % (stale, ADJ-1034 withdrew it);
the group-61 note states 5 %. v1 registry not touched.

## Store work
Copy at `logs/check/v2_v02/store`. `extract --group 61 --write`: 2 rows, 1 non-comparison (2F xfail, as v1). `check --group 61`: no error mentions group 61;
301 errors are store-wide (contract.yaml/adjudications/gaps references, present on the base copy too) and a levels.yaml warning "group 61 is on no level".
Rows have `measured: not_measured` (no `regression --write`; the numbers above are from the plain pytest run). For the parent: add group 61 to levels.yaml
next to 32; `regression --group 61 --write` should be submitted via SLURM (~2 min).
Residual pattern: `residual.json` has key `"61"`, to merge under `patterns`.
Note: the SLURM account hit AssocMaxSubmitJobLimit while other workers' jobs were queued; resubmitted in a loop until accepted.
