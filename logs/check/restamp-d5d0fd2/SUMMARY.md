# Restamp at d5d0fd2 (SLURM 11623360, worktree wt/restamp-a35bba3)

Rows carried from the worktree with `git diff -- docs/validation | git apply` (clean); G40/G42 by hand. contract.yaml untouched. Not committed.
Old run = the run stamped in the row before (mostly 9f72c6b / 73743de / ee50e3a). `in/OUT` = |margin| vs the row bound.

## Flagged
- **G34 S-7 torsion ratio: CROSSED its bound in -> OUT**: 1.00095 (0.0947 %) -> 1.0725 (7.2496 %), bound 5 %.
- G32 V-02 1st flap: stays OUT, 6.07 % -> 7.0865 % (bound 5 %); 1st edge stays OUT, 8.5046 % -> 8.8403 %. Both moved further out.
- G23 berna.3: had no prior measurement; now 4.4589 % (in, bound 5 %). Nothing crossed OUT -> in.
- No anchor crossed OUT -> in.
- Not measured (stay stale, old run kept): G5 (1 comparison pair), G20 (4 comparisons), G23 static_deflection (1) and berna.4-7 (4), G11 (11 unmapped, 2 mapped).

## Per group (json verdicts)
| G | mapped/written | verdicts | non-same/changed/new_baseline |
|---|---|---|---|
| 5 | 0 | unmapped 1 | unmapped (see cause) |
| 11 | 2 | same 2, unmapped 11, declared_non_validation 2 | 11 unmapped (aero-only; BEM tests print no residual the mapper accepts) |
| 20 | 0 | unmapped 1, dnv 1 | unmapped (see cause) |
| 23 | 11 | changed 10, new_baseline 1, unmapped 5, unclaimed 10 | 5 unmapped, 10 unclaimed (test_blade_iea15mw_validation first_modes/parked_modes) |
| 27 | 7 | same 2, changed 5, unclaimed 7 | 7 unclaimed (rated-twist tests, no row) |
| 32 | 2 | changed 2, unclaimed 1 | 1 unclaimed (2nd flap) |
| 34 | 1 | changed 1, dnv 2 | none |
| 43 | 2 | changed 2, dnv 2 | none |
| 58 | 13 | same 13, dnv 3 | none |
No `regression`/`not_measured`/error verdicts anywhere.

## Values (row comparison: old -> new, bound)
- G11 | test_viterna_extension_matches_official.single | rtol at line 371 | 3c02781 None (0.2333%) | None (0.2333%) | rtol 0.01 | in->in
- G11 | test_viterna_extension_matches_official.single | rtol at line 378 | 3c02781 None (0.1275%) | None (0.1275%) | rtol 0.01 | in->in
- G23 | test_blade_mass_matches_published_models.singl | rtol at line 328 | 73743de None (3.7406%) | None (1.3599%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_ccx.0 | rtol at line 355 | 73743de None (0.4423%) | None (0.0096%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_ccx.1 | rtol at line 355 | 73743de None (0.7076%) | None (0.0465%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_ccx.2 | rtol at line 355 | 73743de None (0.7789%) | None (0.0473%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_ccx.3 | rtol at line 355 | 73743de None (0.8358%) | None (0.0475%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_ccx.4 | rtol at line 355 | 73743de None (1.6478%) | None (0.1202%) | rtol 0.05 | in->in
- G23 | test_blade_static_tip_deflection_matches.singl | rtol at line 464 | 73743de None (1.8603%) | None (0.3674%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_berna.0 | rtol at line 516 | 73743de None (1.9404%) | None (1.3791%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_berna.1 | rtol at line 516 | 73743de None (3.3516%) | None (3.3485%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_berna.2 | rtol at line 516 | 73743de None (4.3692%) | None (4.0552%) | rtol 0.05 | in->in
- G23 | test_blade_modal_frequencies_match_berna.3 | rtol at line 516 | None None (None%) | None (4.4589%) | rtol 0.05 | ?->in  <<CROSS
- G27 | test_rated_bem_matches_the_literature_ri.singl | rtol at line 106 | ee50e3a None (0.195%) | None (0.195%) | rtol 0.03 | in->in
- G27 | test_rated_bem_matches_the_literature_ri.singl | rtol at line 113 | ee50e3a None (1.3446%) | None (1.3446%) | rtol 0.03 | in->in
- G27 | test_rated_aero_loads_reproduce_the_bem.single | rtol at line 1020 | ee50e3a None (0.1676%) | None (0.1179%) | rtol 0.005 | in->in
- G27 | test_rated_aero_loads_reproduce_the_bem.single | rtol at line 1027 | ee50e3a None (0.1933%) | None (0.1234%) | rtol 0.005 | in->in
- G27 | test_rated_aero_loads_reproduce_the_bem.single | rtol at line 1034 | ee50e3a None (0.0564%) | None (0.011%) | rtol 0.005 | in->in
- G27 | test_rated_aero_loads_reproduce_the_bem.single | rtol at line 1041 | ee50e3a None (0.0076%) | None (0.0015%) | rtol 0.005 | in->in
- G27 | test_rated_twist_under_production_loads.single | rtol at line 1287 | 4ee98ab None (1.9099%) | None (1.9099%) | rtol 0.02 | in->in
- G32 | testbladenaturalfrequencies_test_1st_fla.singl | rtol at line 133 | 9f72c6b 0.524599 (6.07%) | 0.518922 (7.0865%) | rtol 0.05 | OUT->OUT
- G32 | testbladenaturalfrequencies_test_1st_edg.singl | rtol at line 158 | 9f72c6b 0.695081 (8.5046%) | 0.697231 (8.8403%) | rtol 0.05 | OUT->OUT
- G34 | test_s7_global_torsion_ratio_within_band.singl | rtol at line 107 | 9f72c6b 1.00095 (0.0947%) | 1.0725 (7.2496%) | rtol 0.05 | in->OUT  <<CROSS
- G43 | s3_undamped_free_vibration_f1 | First free-vibration frequency relative  | a0e9279 6.64195 (2.9385%) | 6.65021 (2.8178%) | rel_err 0.05 | in->in
- G43 | s3_newmark_dispersion | Newmark frequency vs trapezoidal dispers | a0e9279 6.64195 (0.0105%) | 6.65021 (0.0062%) | rel_err 0.0005 | in->in
- G58 | test_exact_closed_section_mass.single | rtol at line 129 | 9f72c6b None (0.0%) | None (0.0%) | rtol 1e-10 | in->in
- G58 | test_exact_closed_section_mass.single | rtol at line 133 | 9f72c6b None (0.0%) | None (0.0%) | rtol 1e-10 | in->in
- G58 | test_parked_closed_section_modes.flap_0 | rtol at line 155 | 9f72c6b None (0.4053%) | None (0.4053%) | rtol 0.02 | in->in
- G58 | test_parked_closed_section_modes.flap_1 | rtol at line 155 | 9f72c6b None (0.4368%) | None (0.4368%) | rtol 0.02 | in->in
- G58 | test_parked_closed_section_modes.edge_0 | rtol at line 155 | 9f72c6b None (0.4043%) | None (0.4043%) | rtol 0.02 | in->in
- G58 | test_parked_closed_section_modes.edge_1 | rtol at line 155 | 9f72c6b None (0.4785%) | None (0.4785%) | rtol 0.02 | in->in
- G58 | test_parked_closed_section_modes.torsion_0 | rtol at line 155 | 9f72c6b None (0.0173%) | None (0.0173%) | rtol 0.02 | in->in
- G58 | test_rotating_southwell_same_s4_path.flap_0 | rtol at line 181 | 9f72c6b None (0.1883%) | None (0.1883%) | rtol 0.03 | in->in
- G58 | test_rotating_southwell_same_s4_path.edge_1 | rtol at line 181 | 9f72c6b None (0.6039%) | None (0.6039%) | rtol 0.03 | in->in
- G58 | test_explicit_reference_surface_inertia.mid | rtol at line 202 | 9f72c6b None (0.3413%) | None (0.3413%) | rtol 0.005 | in->in
- G58 | test_explicit_reference_surface_inertia.outer | rtol at line 202 | 9f72c6b None (0.3235%) | None (0.3235%) | rtol 0.005 | in->in
- G58 | test_composite_closed_box_effective_clt.ea_0_1 | rtol at line 217 | 9f72c6b None (0.0%) | None (0.0%) | rtol 1e-05 | in->in
- G58 | test_composite_closed_box_effective_clt.ei_1_0 | rtol at line 217 | 9f72c6b None (0.3323%) | None (0.3323%) | rtol 0.005 | in->in

(`raw` is None for rtol rows; the margin is the number shown. G27/G58/G11 rows mostly re-stamped with identical margins: they do not depend on the blade geometry changed in the audit, or the margin is unchanged to 4 digits.)

## G5 / G20 cause
Not an environment or collection error. Stderr empty; the batch log shows only an OpenFabrics warning. Both runs exit 1 because the mapper returned `unmapped` ("0 printed residual(s) for N asserted comparison(s)"), which is in REGRESSION_FAILING_VERDICTS. G20's 33 s is just the test's runtime: a login-node run of `tests/validation/blade/test_blade_twist_anchor_beam.py` passes both tests and prints a table, but no `signed error ... bound` residual lines the mapper pairs with the 4 asserted comparisons. G5 (blade_mesh_convergence, 2 comparisons) is the same mapper pattern; I did not rerun it (long). Rows 5 and 20 keep their old measurements (G20 at e1b7da7): still stale; the fix is a test or mapper change, outside this task.

## Hand groups
S-4 (G40) from pytest-s4.log (the same printed lines as capture-s4.pkl; the pickle holds the shell frequencies, Hz, in test order 1F, 1E, 2F per rpm and was used as cross-check, values agree). Mapped by label and test node; rated 1E uses the driver reference as the row does.
| comparison | old (4d64ea8 / 2eab2c6) | new d5d0fd2 | ref Hz | bound |
|---|---|---|---|---|
| parked 1F | +3.4237 % | 0.518922 Hz, +2.3348 % | 0.507083 | 5 % in |
| parked 1E | +0.2043 % | 0.697231, +0.5559 % | 0.693376 | in |
| parked 2F | +1.7977 % | 1.483824, +0.1974 % | 1.4809 | in |
| rated 1F | +2.8319 % | 0.541819, +1.8451 % | 0.532003 | in |
| rated 1E | +0.1210 % | 0.702508, +0.4424 % | 0.699413 | in |
| rated 2F | +1.6196 % | 1.508960, +0.0598 % | 1.508057 | in |
The old texts' es 0.125 numbers and the ElastoDyn comparison were dropped from the text (not re-measured). Note the ElastoDyn parked comparison printed in the log reads -4.61 % / -5.65 % (comparison only). The test centrifugal-stiffening 1F rise: 0.5189 -> 0.5418 Hz (passes).
V-05 (G42): the capture-v05.pkl is empty (the plugin only wraps the S-4 helper and picks fixtures not used by this test); the value comes from the printed line in pytest-v05.log: V-05 assembled mass 67,150.0 kg vs 66,996.9, +0.2285 % (old -0.4695 %), bound 5 %, in. 6 passed, 1 xfailed, 3 xpassed (S-6 included, no S-6 row touched).
Not mapped: the S-6 tests (no V-05/S-6 row edited), the K55 XPASS cases, and G42 beyond the single mass comparison.

## Checks
- `tools/validation_matrix.py check`: `check: 313 row(s), 432 comparison(s), 0 error(s), 1 warning(s)` (warning: levels.yaml group 27 is on no level).
- `scripts/check.sh quick`: 75 passed, 1 xfailed; OK (log logs/check/quick-20261024T103523.log).
