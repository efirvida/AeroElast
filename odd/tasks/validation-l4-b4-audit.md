# Feature: B4 audit of L4 (BEM aerodynamics against OpenFAST AeroDyn)

Issue: efirvida/AeroElast#33, item B4. Branch `audit/b4-l4` from `origin/main` 2b95e5f, worktree `logs/worktrees/b4-l4`
(separate from `audit/b3-1-l3-inventory`, where B3 K7 jobs were running at entry). L4 `requires: []`, so it does not wait on B3.
Audit work: produces measured verdict proposals; the maintainer closes or reopens. No band is widened, no verdict precedes its criterion.

## Specs

- S1 (user, verbatim): "vamos a arrancar con L4"
- S2 (roadmap #33, B4, verbatim): "| B4 | Audit L4 | as B1 |", where B1 reads: "every L1 row has a verdict at one revision, with
  its arbiters named; the documented reds (...) are each an exception with a source, or rejected; every `atol` comparison (23 in
  the store) is classed or given a relative equivalent, because the ceiling cannot judge an absolute bound"
- S3 (roadmap acceptance rule, verbatim): "A **physics comparison** passes at **≤ 5 % relative error** against an independent
  reference (...)" and "**Numerical references settle it only in pairs:** agreement with a single model is not validation, and at
  least two independent ones are needed. With no arbiter available, the quantity is **not validated**, and the record says so."
- S5 (user, verbatim, amends S3 for IEA 15 MW): "para la turbina IEA 15mw aunque no tenemos arbitro analitico o pulicaciones
  tenemos el repositorio oficial de OpenFAST por lo que ese arbitro es valido"; scope (L4): "Todo OpenFAST oficial (L3 y L4)".
- S4 (roadmap, verbatim): "A closure ships three things in one push: 1. the verdict in `contract.yaml`; 2. the closure section;
  3. the comment on the item." Push and the issue comment are the maintainer's.

## Tasks

- [x] T1 (S2) Entry measurement and environment at 2b95e5f - inline - this commit
- [x] T2 (S2, S3) Per-row inventory and source findings (below), no verdict - inline
- [ ] T3 (S2, D2) Import V-01 and `test_viterna_matches_aerodyn_theory_manual` as manual rows; `bem_load_frame` to out_of_scope - inline
- [x] T4 (S3, S5) Amend the arbiter rule: `reference_model` (official OpenFAST model on the official deck) settles alone, L3 and L4;
  tool + test-first in `tools/validation_matrix.py`, roadmap mirror amended, issue text proposed - inline
- [-] T4-old (D1) second-arbiter search: dropped by L3/L4 (the official model is the arbiter)
- [x] T5 (S2, S3) Predeclared verdict criterion, written before any re-measurement - inline
- [ ] T6 (S2, D3) Adapt the AeroDyn driver writer to OpenFAST v5.0.0 (test-first: v5 setup fails today) - inline
- [ ] T7 (S2) Re-measure group 11 + imported rows at the audited revision with v5 - SLURM dev
- [ ] T8 (S2, S4) Verdicts in `contract.yaml` + closure section; issue comment text proposed, not posted

## Entry measurement (T1, 2026-10-20, 2b95e5f, read-only)

- `contract --level L4 --rows` (`logs/check/b4/contract-l4-rows-entry.log`): L4 open; unaudited 12, absent 3, over ceiling 0,
  single or no arbiter 15, gaps 0. `check`: 289 rows, 399 comparisons, 0 errors, 0 warnings.
- Store: all 12 rows of group 11 measured at `1d9a665` (2026-10-04), 655 commits before 2b95e5f; since then `solvers/bem/engine.py`
  changed (4 lines) along with force projection / participant files. Stored margins are not a measurement at the audited revision.
- Environment: the store and the test labels name **OpenFAST 4.2.1**. On this cluster the only binaries are
  `conda-envs/openfast/bin/openfast` and `tools/openfast-env/bin/openfast`, both **v5.0.0** (GCC 14.3). `find_openfast_bin`
  finds neither (it looks at `<venv>/../openfast` and `~/miniconda3/envs`), so without `OPENFAST_BIN` all 12 rows skip.
  `docs/validation-environment.md:30-36` records that v5 refuses the 4.x driver input this test writes.

## T2 inventory (source `tests/validation/bem/test_bem_openfast_parity.py` at 2b95e5f)

| rows | quantity | bound | stored | reference as labelled | finding |
|---|---|---|---|---|---|
| `test_rotor_performance_matches_aerodyn.{0,1,2}` (A1, V 8/9/12, 7.56 rpm, pitch 0) | thrust, torque | rtol 0.015 | 0.05-0.79 % | AeroDyn 4.2.1 | physics row; single numerical arbiter |
| `test_yaw_and_shear_match_aerodyn.{0,1}` (A2, yaw 10 deg / shear 0.2) | azimuth-avg thrust, torque | rtol 0.015 | 0.33 / 0.27 % (row .0) | AeroDyn 4.2.1 | physics row; single numerical arbiter |
| `test_spanwise_loads_match_aerodyn.{0,1,2}` | interior max abs d-alpha; mean abs d-Cn | atol 2.0 deg; rtol 0.025 | - | "the official AeroDyn polar tables, interior stations" | F1: the label is wrong, the reference is the AeroDyn run (`ref.alpha_deg`, `ref.cn`); F2: the 3 `atol` of L4 are these d-alpha bounds |
| `test_bem_with_repo_default_polars_matche.{0,1,2}` (A3) | thrust, torque with the WindIO polars | rtol 0.03 | up to 2.12 % | AeroDyn 4.2.1 | F3: its own docstring says "a **sensitivity**, not a parity claim": different inputs on the two sides, so not a physics comparison of the BEM |
| `test_viterna_extension_matches_official.single` | worst Cl, Cd over 48-88 deg | rtol 0.01 | - | official polar tables | F4: AR is fitted from the reference's own Cd(90), so Cd near 90 deg agrees by construction; F5: module docstring says 0.23 % / 0.13 %, test docstring and comment say 0.48 % / 0.18 % |

Non-validation by `groups.yaml`: `test_neuralfoil_generated_polar_differs_from_official`, `test_reference_polars_are_the_official_aerodyn_tables`.

## T5 predeclared verdict criterion (written BEFORE any v5 number, 2026-10-20; T6 RED seen: v5 fails parsing `a1.dvr`)

The 4.2.1 numbers in the store (1d9a665) are known; they are history, not the baseline, and no v5 number exists yet.
- K1 measurement: each row's own printed metric, at the branch head after T6, `aerodyn_driver` of OpenFAST **v5.0.0** on the
  vendored official deck with the `write_bem_primary` overrides (quasi-steady, tower-free, BEM_Mod 1); BEM side unchanged.
- K2 arbiter: `reference_model` "OpenFAST v5.0.0 AeroDyn, official IEA 15 MW deck" (S5) for A1, A2, spanwise and the Viterna
  like-for-like rows. The verdict note states it is verification against the turbine's definition, not physics.
- K3 ceiling: |e| <= 5 % -> `trusted`; |e| > 5 % -> `rejected`. A row that is <= 5 % but fails its own tighter source bound
  (1.5 %, 2.5 %, 1 %) is a moved anchor: reported as a finding and left `undecided`; no bound is touched.
- K4 the `atol` (spanwise max |d alpha| <= 2 deg, the 3 L4 atols): relative equivalent fixed now as
  e_alpha = max_i |d alpha_i| / max_i |alpha_AeroDyn,i| over the same interior stations, judged by K3. Station-wise
  |d alpha|/|alpha| is not used (ill-conditioned where alpha -> 0).
- K5 AeroDyn time convergence: one control run with dt halved (A1: 0.5 s, A2: 0.125 s), same tmax. |delta e| <= 0.1 pp per
  comparison -> converged; otherwise the row is `undecided` whatever its e.
- K6 class: A3 `test_bem_with_repo_default_polars_matche.{0,1,2}` compare different inputs (its docstring: "a **sensitivity**,
  not a parity claim") -> class `regression_pin`, outside the ceiling. Imported Viterna closed form (T3) -> `analytical`.
  Imported V-01 vs the WISDEM tabular (CCBlade, the same BEM code as ours, not OpenFAST) -> `published_numerical`, single and
  not independent -> `unarbitrated` whatever e. F1 label is corrected to the AeroDyn run before the verdict.

## Decisions (asked and answered 2026-10-20, all as recommended; see L2)

- D1 arbiter: every row rests on one numerical code (AeroDyn). By S3 that leaves them unarbitrated unless a second independent code
  or a closed form exists.
- D2 absent: import the three absent L4 tests as rows now (A6c slice) or keep them absent with B4 open.
- D3 environment: re-measure with v5.0.0 (adapt the driver writer; the arbiter version changes) or install 4.2.1 (dependency
  install, item 5) to reproduce the stored numbers first.

## Log

- L1 (user, 2026-10-20): "vamos a arrancar con L4"
- L2 (user, 2026-10-20, answers): D1 "Unarbitrated ahora + buscar 2º código (Recomendado)"; D2 "Importar V-01 y Viterna;
  load_frame fuera (Recom.)"; D3 "Adaptar el driver a v5 y medir con v5 (Recomendado)".
- L3 (user, 2026-10-20, mid-turn, revises D1): "para la turbina IEA 15mw aunque no tenemos arbitro analitico o pulicaciones
  tenemos el repositorio oficial de OpenFAST por lo que ese arbitro es valido"
- L4 (user, 2026-10-20, scope of L3): "Todo OpenFAST oficial (L3 y L4) (Recomendado)" - ElastoDyn/BeamDyn/AeroDyn on the official
  deck count as a sole arbiter; L3 verdicts are re-adjudicated on the B3 branch, not here. Agent caveat recorded: with identical
  inputs on both sides this is verification against the reference implementation, not validation against physics.
- L5 (agent, T4 evidence): RED `test_contract_the_official_reference_model_arbitrates_alone` (contract refused kind
  `reference_model`), GREEN after adding it to `SOLE_ARBITERS`. Non-slow `tools/tests`: 185 passed, 5 failed
  (`logs/check/b4/tools-tests-T4.log`); the same 5 fail at clean 2b95e5f (`logs/check/b4-base-5fail.log`: coherence x2,
  real store coherent, no pending migration, depends_on covers L1/L2 sources) -> pre-existing store line drift, not this change.

