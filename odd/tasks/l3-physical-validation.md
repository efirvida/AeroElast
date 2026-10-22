# L3 physical solver validation from the matrix store
Reproduce physical residuals, diagnose their cause, correct demonstrated implementation errors, and record measured evidence in the validation store.
Branch: `audit/b3-s2-beamdyn-static`; baseline `3e381af`; delivery ask-on-risk, local commits/no push; forecast depends on demonstrated defect (initial execution/evidence slice <400 lines); runner `scripts/aeroenv.sh python -u` on SLURM for blade solves.

## Specs
- S1: "yo estoy validando un solver fisico, por favor haz las validaciones, utiliza el store de la matriz de validacion, identifica errores en los modelos fisicos implementados en el proceso, y corrigelos". Use the store as authority; reproduction, independently referenced physical quantities and accepted tolerances decide outcomes, not software smoke.
- S2: "K14 VOID; K15 condicional; LC3 undecided; D1-D3 pendientes." "No tocar contract.yaml ni ejecutar E/F. No repetir K13 completo." Retain these boundaries; a new physical measurement does not silently promote a verdict.
- S3: "Yaml oficial; malla 0.15–0.25 m, 0.125 solo convergencia." "Python con scripts/aeroenv.sh y -u." "SLURM solo --ntasks=N; M1 0.125 en sequana_cpu, 16 tasks, nunca dev." "Esperas background sin polling. Sin push." Baseline reproduction uses existing official es 0.25 fixture; no M1/.125 launch is needed initially.
- S4: "Preservar todos los directorios untracked de referencias BeamDyn." Never overwrite or rerun in those directories. No reinstall is needed.
- S5: "no tengo idea de que me estas hablando". Do not ask the human to decide technical prerequisites; communicate physical results, demonstrated errors, corrections and remaining validity limits in plain language.
- Assumption: continue the current L3/B3-3 scope first, using store group 39 (official S-2 cases) rather than open unrelated solver campaigns. Never infer that a shell/beam discrepancy by itself proves a coding defect or change a pass band to hide it.

## Tasks
- T1 | S1–S5 | parent execution + delegated diagnosis (context-backstop) | done | Baseline3e381af source pinned, old extension hash frozen unchanged; physical run launched as background job4; triangle defect isolated, compliant RED required in T2.
- T2 | S1–S5 | delegated bounded writer (context-backstop) | done | RED 4 failed/8 passed (91.6667%); fix in `assembler.rs`; release rebuild; compiled GREEN 26 passed; 12 rows registered in group 28 (not_measured until recorded at a revision); commit: this work unit.
- T3 | S1–S5 | parent execution + independent high-risk verification | pending | Collect immutable baseline, run corrected official es0.25 four-case physics and matching element anchors; full-suite measurement on SLURM after Rust rebuild; quantify blade impact without claiming entire gap fixed.
- T4 | S1–S5 | parent store reconciliation | pending | Reconcile physical measurements and new analytical regression with store, check it, record remaining validity limits/local commits; no contract.yaml verdict edits.

## Log
### L1 — User correction (verbatim)
```text
no tengo idea de que me estas hablando, yo estoy validando un solver fisico, por favor haz las validaciones, utiliza el store de la matriz de validacion, identifica errores en los modelos fisicos implementados en el proceso, y corrigelos
```
### L2 — Reconciled starting evidence
- Baseline HEAD `3e381af`, branch as above, tracked clean with the same six untracked BeamDyn reference directories. K16 input-only tooling is not physical evidence and does not close any solver validation.
- `validation_matrix check`: 301 rows, 420 comparisons, 0 errors, 1 warning (group27 has no level), exit0. Store group39 contains official LC1/LC4 asserted failures at 5%, LC2/LC3 info_only due non-clean node/ring metrics. Existing measurements are from earlier revisions, not this run.
- Read-only mapper identified existing physical test `tests/test_iea15mw_s2_static_prescribed.py`, fixture es0.25 official mid-surface, all four static load cases; known LC3 discrepancy ~17%. Tip loads and measurement selection need to be distinguished from actual solver constitutive errors. No concrete solver defect is established yet.
- One batched global status/help/scheduler query timed out at 30 s; stdout/stderr were not returned by the tool exception. This is a failed query, not physical evidence; do not retry that command. Continue from the scoped store row and launch scheduler interactions as background events.
- Execution cost: prior S-2 accounting now observed (exit0): job11621347 used8tasks,3m41s; job11621034 used12tasks,5m21s. Both regular CPU partitions currently have0 user running jobs. Use sequana_cpu with --ntasks=12 to reproduce the demonstrated allocation rather than underallocate memory; no other resource-shape flags. Wait in background without polling.

### L3 — Physical defect isolated; baseline protection before rebuild
- Read-only diagnosis observed triangle uniform-strain energy0.0503778337531J versus independent orthotropic plane-stress/CLT0.604534005038J (-91.6667%, exactly E2/E1=1/12) for requested fibres+Z. Quad1.20906801007J versus1.20906801008J agrees. It implicates production material-frame mapping, not gravity: gravity resultants matched rho*h*A*g.
- Source: `crates/aeroelast-py/src/assembler.rs:224–244` corrects composite quads only; `crates/aeroelast-core/src/elements/mitc3.rs:162–190,316` uses a projected global-X reference and subsequently rotates to element-local. Do not reuse the quad-local correction blindly. Blade impact remains unmeasured.
- Existing tiny physical anchors: cantilever3.5115%/5%, axial extension2.0103%/5%, direct MITC3 composite bending0.6667%/5%, direct quad axial3.4578%/5%, exact mass and laminate shear controls pass. These bypass the defective production span-direction path.
- Diagnostic agent used noncompliant /tmp scratch; do not access/remove it. Reproduce regression under required managed SCRATCH/tmp before correction; source remains unchanged. Freeze baseline extension before rebuilding; observed old extension SHA256f1d3ef01382fb25ba2d823dab988c029c81c70233efe4a633dd9cda19d0c3d2e, path `src/_aeroelast/_aeroelast.cpython-312-x86_64-linux-gnu.so`.
- Risk: item3 (production assembler material preparation is shared by existing element paths; independent verification and preservation anchors required). Rust change additionally requires release rebuild and full-suite measurement before any merge. No merge/push is authorized.

### L4 — Correction compiled; baseline job invalidated
- Release rebuild (`maturin develop --release`) exit 0 in 1m30s. Imported extension is `venv/lib/python3.12/site-packages/_aeroelast/_aeroelast.cpython-312-x86_64-linux-gnu.so`, SHA256 `9e183ce64a3f3de4482d2882e35881ab66b460a54bea3bf1be3ca1199ea565b8`.
- Compiled GREEN (parent, `logs/check/l3-triangle-frame-green.log`): 26 passed in 4.68 s, exit 0. Production triangle +Z fibre energy 0.604534 J vs CLT 0.604534 J (0.0000 %, bound 5 %); four triangle frame/projection cases 0.0000 %; default/normal-span/isotropic tri+quad unchanged 0.0000 %; cantilever 3.6958 %, axial 1.8580 %, direct MITC3 ply0 bending 0.6667 %, direct quad ply90 axial 3.4453 %, exact mass 0 %, laminate shear coupling 8 passed.
- **Baseline job 11621827 is INVALID as a before-measurement.** It loaded the ignored workspace `src/_aeroelast/*.so` (mtime 2026-09-25, SHA256 f1d3ef…), which predates the Rust fixes `41e9a16` and `c0fdeea` (2026-10-19). Its numbers (LC1 −7.3691 %, LC2 −3.3180 %, LC3 +16.0884 %, LC4 −6.1805 %) are a stale-binary run and are not cited. `logs/check/l3-physical-baseline-extension/` holds that stale binary, not a valid baseline.
- Valid before-measurement: the store's group-39 values (es 0.25, runs 44ebffa/6318955, 2026-10-21). `git diff 44ebffa HEAD -- src/ crates/` is empty; only test printing and the in-house beam comparison changed. The after-run therefore differs from that baseline only by `assembler.rs`.
- Independent verifier (Risk item 3) did not complete: its 17 G workspace copy hit the 30 s timeout (not retried), then the subagent errored. The store-registration worker also errored. Both were runtime failures, not findings; registration was done inline instead. A 99 M partial copy remains at `$SCRATCH/tmp/l3-t2-verify.VkWjf0`.
- Store: extractor `--write` on a copy would also relabel nine existing group-28 comparisons and add an unrelated row (pre-existing extractor drift, out of scope). Only the 12 new production-frame rows were appended; existing 24 rows byte-identical. `validation_matrix check`: 313 rows, 432 comparisons, 0 errors, 1 warning (group 27 on no level).
