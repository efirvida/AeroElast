# Blade mesh size: the rule and the evidence behind it

## The rule

**The blade is meshed at `element_size` 0.15-0.25 m in every test, diagnostic, campaign case and
documentation example. 0.5 m is not a converged mesh for this blade and its results are not
citable as evidence of anything.** It exists in the repository only inside the mesh-convergence
studies already done (below), where it is the coarse point of the sequence.

- A new test or tool that meshes the blade at `element_size >= 0.3` is wrong unless it IS a
  convergence study, says so in its docstring, and reports at least two finer points.
- 0.125 m is the reference mesh when a convergence claim is made (Richardson fit within ~0.2 %
  of the limit, see 1). 0.25 m is the working mesh.
- A coarser mesh needs a written convergence statement for the quantity in question, in the same
  commit. "It ran faster" is not one.
- Never re-band an anchor to absorb a mesh change: re-measure it on the finer mesh and report
  the move (AGENTS.md, "Physics discipline").

## Evidence (all already in the repository; nothing here is new measurement except 5)

1. **Static flap tip, coupled convergence study** (`docs/origin_main_integration_2026-09-30.md`,
   section "Convergence", `tests/IEA15MW/convergence_b1_b2`): h 0.5 -> 0.25 m moves the flap tip
   7.7344 -> 7.3805 m (**4.58 %**); h 0.25 -> 0.125 m moves it 7.3805 -> 7.3187 m (**0.84 %**),
   observed order ~2.4, Richardson limit ~7.307 m. The document's own verdict: "0.125 m is the
   usable reference mesh and 0.5 m is the outlier".
2. **Blade torsion** (`docs/validation_closures.md`, 2026-09-19 MITC3 batch): `GJ_beam/GJ_shell`
   0.932 at 0.5 m vs 0.954 at 0.25 m (CalculiX 0.930 / 0.953).
3. **Aeroelastic twist** (same document, item 4): about 9.9 deg at 0.5 m vs 5-6 deg at 0.25 m on
   the pre-fix element (so the size of the gap is not a current number, its direction and order
   are); `docs/twist_distortion_diagnosis.md`: S-8c at 0.125 m gives mid +1.57 deg, tip
   +6.43 deg, i.e. the excess twist persists under refinement.
4. **Modal convergence** (`docs/validation_data/generated/s2e_mesh_convergence.csv`, 1.0 / 0.5 /
   0.25 m): the first edgewise frequency is already within 0.03 % at 0.5 m (0.697503 vs
   0.697671 Hz) but `tx_2e` still moves 11.22 -> 11.12 %. **Frequencies converge earlier than
   displacement and torsion**, which is why a coarse mesh can look fine in a modal test and be
   wrong in a static or coupled one. The rule does not rely on modal agreement.
5. **Zhou-load shell twist** (2026-10-13, `tools/diagnose_zhou_loads_reverse.py`): +5.955 deg at
   0.5 m, **+5.855 deg at 0.25 m** (-1.7 %). A small move for this metric; the 0.5 m figure was
   discarded anyway because the metric was not known to be insensitive beforehand.

## Known debt: places that still mesh the blade at 0.5 m (inventory 2026-10-13)

Not edited, because each one carries a measured anchor or an exact-count pin that has to be
re-measured on the finer mesh and reported, not silently moved. Re-baseline them one at a time,
each with its own before/after.

- Tests with anchors: `tests/test_iea15mw_v03_static_gravity.py` (3 meshes; its own comment
  admits the body-load integral error at 0.5 m), `tests/test_iea15mw_s2_static_prescribed.py`,
  `tests/test_iea15mw_v05_structural_properties.py`, `tests/test_iea15mw_s1_sectional.py`,
  `tests/test_iea15mw_s4_rotating_modal.py` (`_build_assembly`),
  `tests/software/mesh/test_blade_mesh.py` (exact node/element counts).
- Anchor-producing tools: `tools/run_s5_oneway.py`, `tools/run_s6_oneway_dynamic.py`,
  `tools/run_sx_campbell.py`, `tools/run_s8a_gj_circumferential.py`,
  `tools/diagnose_multicell_convergence.py` (docstring about refining to 0.5 m).
- Intentional convergence points (keep, they are the coarse end of the study):
  `tests/IEA15MW/convergence_b1_b2/*`, `tests/IEA15MW/_short_hcoarse*/solid_h_coarse.yaml`.
- Factual remark kept: `src/aeroelast/core/mesh/generators.py:1177` (upstream's own test pins
  9277 nodes at 0.5 m).
- Historical measurements quoted at 0.5 m (`docs/twist_distortion_diagnosis.md`,
  `docs/origin_main_integration_2026-09-30.md`, `docs/validation_closures.md` item 4): records,
  not instructions; left as written.

Already moved to 0.25 m on 2026-10-13: the CLI/doc examples (`cli/run_fsi.py`, `run_bem_fsi.py`,
`run_aero_fsi.py`, `docs/cli-reference.md`, `docs/mesh-cli.md`) and the diagnostics without an
anchor (`diagnose_newmark_only`, `diagnose_load_point_sign`, `diagnose_rotating_assembly`,
`diagnose_damping_matrix`, `debug_kg_magnitudes`, `diagnose_zhou_loads_reverse`,
`diagnose_mixed_mesh_consistency`).

## Re-calibration of the anchors: deferred until it is mandatory (decision 2026-10-13)

Moving the debt list above to 0.15-0.25 m changes recorded anchors, so each one needs its own
re-measurement and a reported before/after. That work is **deliberately postponed** and becomes
mandatory when any of these happens:

- an anchor from the debt list is cited as evidence in a report, paper or review answer;
- a change touches the code path one of those tests/tools exercises (element, K_G, projector,
  mass, loads) and the result has to be believed rather than just kept green;
- a new test is added next to them and would inherit their mesh;
- the full suite is re-recorded (`scripts/check.sh full --record`).

Until then: those numbers are labelled "0.5 m, not converged" wherever they are quoted, and no new
result is built on top of them. Do not "fix" the debt by editing a band or a mesh size without
re-measuring.

## Correction to the inventory above, and the guard (2026-10-13)

The inventory above searched 0.5 m only. A scan of every `element_size >= 0.3` found **45 files**
(`scripts/mesh-size-debt.txt` is the authoritative list): 7 convergence-study points (keep), 7
geometry-only (turbine, tower, aero surface, fluid meshes; no structural blade result), 10 at
0.5 m, and **21 at 1.0 m** - about fourteen validation tests under `tests/validation/` and
`tests/test_multicell_shear_flow.py`, and ten diagnostics under `tools/`, all coarser than the
0.5 m this document is about. Their numbers are not converged either: in `docs/validation_data/generated/s2e_mesh_convergence.csv`
the 1.0 m row is further from the 0.25 m one than the 0.5 m row (`f1e` 0.696398 vs 0.697503 vs
0.697671 Hz; `tx_2e` 11.32 vs 11.22 vs 11.12 %). They are
debt of the same kind and follow the same deferral and the same triggers.

`scripts/check_mesh_size.py`, called by `scripts/check.sh` (quick and full), fails on any file or
value not in the registry: the debt can only shrink. Removing a file from the registry requires
having re-measured it on 0.15-0.25 m. Adding one needs a reason written in this document.
