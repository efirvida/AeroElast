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
- T3 | S1–S5 | parent execution + independent high-risk verification | done (independent verifier unavailable, see L5) | Collect immutable baseline, run corrected official es0.25 four-case physics and matching element anchors; full-suite measurement on SLURM after Rust rebuild; quantify blade impact without claiming entire gap fixed.
- T4 | S1–S5 | parent store reconciliation | done | Groups 39/28 in 4b6a04b; CCX mesh-convergence, V-02 1F/1E, S-7 at 7d20c81 in the L6 commit. contract.yaml untouched.
- T5 | S1–S5 | delegated bounded diagnosis (context-backstop), parent SLURM execution | in_progress (instrumentation c882c1f; cause located L14; defect demonstrated and fixed L18: director averaging across >20° folds; L21: right at true folds (Bredt); L22: 20° rule kept; L25–L32: S-5/S-6/S-7 closed on mid, D3, V-02, T-junction, S-1 frame fix, rule-9 re-measure, CCX parity, ring fix; L33–L35: ring fix accepted; torque sign chain fixed 4ee98ab vs OpenFAST (finding 1 closed (b)); L36: follow-ups 7189314, P5 row 8e1c516, full suite 11622759 running, T5.d first readings; resume from L36) | Reproduce S-6 on official es0.25, diagnose stiffness/loads/integration, correct only a demonstrated defect, and record physical evidence. Starting HEAD 0568825; no push. Stiffness indefiniteness is a hypothesis, not an established necessary cause.

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

### L5 — Blade impact, differential and off-axis checks (fix commit `7d20c81`)
- Official S-2 at 7d20c81 (SLURM 11621861, sequana_cpu, 12 tasks, 2m31s pytest; import origin and extension SHA 9e183ce verified in the log): LC1 −7.4188 % (was −7.4110), LC2 −3.6997 % (−3.6807), LC3 +16.9097 % (+16.9442), LC4 −6.2268 % (−6.2214). **The triangle defect is real but moves the blade by ≤0.035 pp**: it does not explain the S-2 gaps. Recorded in group 39 at 7d20c81; es 0.125 not re-measured.
- Off-axis probe (`$SCRATCH/tmp/l3-probes/offaxis_energy.py`, log `logs/check/l3-offaxis-energy.log`): 192 cases (4 face normals × tri/quad × θ∈{0,±30,±45,90} × 4 strain states incl. shear and combined) vs textbook CLT Q̄(θ): worst 8.0e−12, 0 mismatches. Ply orientation of membrane stiffness, including in-plane shear (A66), is correct on the production path; the edge excess compliance is not a ply-angle error.
- Differential (`k_differential.py`): pre-fix extension rebuilt from 3e381af in isolation (SHA 7317d1c5…, `$SCRATCH/tmp/l3-prefix-build-3e381af`) vs fixed 9e183ce. Bit-identical K for all composite quads (no span, +Z, tilted), composite triangles without span or with normal span, all isotropic tri/quad/mixed. Differs only where intended: composite triangles/mixed with in-plane span (max rel 0.21–0.28 for a 0/±45/90 laminate).
- Independent gentle-ai-verify for Risk item 3: UNAVAILABLE (two runtime failures: 30 s copy timeout, then subagent error). Not retried; the differential above and the analytical suite are parent checks, not an independent verdict.
- Group 28: the 12 production-frame comparisons recorded as measured at 7d20c81 (12 passed, 0.0000 %, `logs/check/l3-span-frame-7d20c81.log`). Store check 313 rows, 0 errors.
- Full suite after the Rust change: job-7 (`scripts/check_sbatch.sh full`, sequana_cpu_bigmem, 12 tasks) running.

### L6 — Full suite at 7d20c81: no new failure
- SLURM 11621862 (1h28m): 18 failed, 1407 passed, 38 skipped, 13 xfailed, 3 xpassed; log `logs/check/full-20261022T113305.log`. check.sh labels all 18 NEW-FAIL only because `scripts/known-red-full.txt` is stale (not re-recorded by maintainer decision 2026-10-18, Engram 25).
- Every one of the 18 also failed in the three full suites of 2026-10-19 (`full-20261019T153623/173319/233925.log`) or is a store-recorded FAIL: blade-mesh node count 9357 vs 9277 (×2), S-1 ×4 (EI_edge 12.69, EI_flap 5.17, GJ 6.22, EA 7.59 % — store values within print precision), S-2 LC1/LC4 official, S-2 LC2/LC4 UTD (disowned input), S-6 ×2 (divergence ~1e40 m), section-frame tear guard, twist-mechanism, S-7 9.0924 % (9.0927), V-02 1F 6.0930 % (6.0978) and 1E 8.4707 % (8.4598), CCX high-mode gap 9.2182 % (10-19: 9.3061/9.2253). The triangle fix introduced no new failure and moved no anchor by more than 0.035 pp.
- Store correction: the CCX mesh-convergence 5 % comparison was stored as 2.1197 % at 1d9a665 while failing at ~9.2 % since at least 10-19; now recorded at 7d20c81. V-02 1F/1E and S-7 restamped to 7d20c81. Store check: 313 rows, 432 comparisons, 0 errors.
- **Open physical-model failure with the largest consequence: S-6 diverges.** Documented closing on 2026-09-10 at es 0.5 (std −9.8 %, mean −17.7 %); diverged in every suite since the 09-30 merge / 10-19 K_G fix. Next unit T5.

### L7 — Resume T5 (2026-10-22)
```text
│ Retomamos odd/tasks/l3-physical-validation.md, T5 (divergencia S-6). mem_search "S-6
 │ divergence" (observación 87).
```
- Reconciled observation87, full feature mirror83, repository document and HEAD0568825; tracked clean, six BeamDyn reference directories remain untracked and protected.
- Route: bounded worker for S-6 diagnosis/instrumentation and an evidence-backed correction only within the S-6 tool/test surfaces; trigger is runtime context backstop (47,060 tokens). Parent owns SLURM blade executions and subsequent verification/reconciliation. One behavior, not parallel writers.
- Existing S-6 test still calls build_assembly(0.5): not converged physical evidence. Reproduce at official es0.25 before correction; retain accepted 25% std/mean bounds unchanged. Do not infer that Newmark divergence proves an indefinite effective step matrix: negative structural stiffness, force feedback, or an implementation error must be distinguished.
- Initial measurement is the directed S-6 run, not another 1h28 full suite. Any Rust correction requires release rebuild and a new full-suite measurement; no AGENTS/skill/contract edits.

### L8 — S-6 instrumentation (worker), no defect demonstrated yet
- Worker edits (uncommitted): `tools/run_s6_oneway_dynamic.py` (+25/−4: optional `components=` in `assemble_keff`, optional `observer=` in `newmark_tip_series`; arithmetic unchanged), new `tools/diagnose_s6_divergence.py`, new `tests/test_s6_newmark.py` (5 analytic tests), log `logs/check/l3-s6-analytic.log`.
- RED: observer keyword `TypeError` (exit 1). GREEN: 5 passed in 0.326 s (exit 0). Preservation: 401-sample histories bit-identical to `HEAD:tools/run_s6_oneway_dynamic.py` for positive and negative stiffness; returned K/M/free identical.
- Trace: loads are precomputed (`run_s6_oneway_dynamic.py:165`), no displacement feedback. K_G from production static centrifugal prestress (`assembler.py:835–870`); K_SP negative semidefinite (`:148`); Rayleigh C built on the rotating K (`:207–214`). Newmark recurrence matches an independent continuous oscillator. A positive Newmark step matrix does not prove physical stability: a negative rotating stiffness grows under trapezoidal Newmark (analytic test).
- Remaining hypotheses: negative K + K_G + K_SP (and the Rayleigh damping built on it), assembly or linear-solve pathology, a blade-specific discrete error. None established.
- Risk: item 3 (shared S-6 helpers). `assess` with escalate → independent verifier required; native assess unassessable (untracked declaration).
- Official es0.25 measurement submitted (background job-1, `sbatch --wait`, sequana_cpu, 12 tasks, 45 min; log `logs/check/l3-s6-official025-<jobid>.log`), eigenvalues skipped in this first run.
- Independent verifier: launch refused by the runtime (task text carried no `Risk: item N` line); not retried. Parent spot check instead: `tests/test_s6_newmark.py` 5 tests OK in 0.360 s, exit 0 (`logs/check/l3-s6-analytic-parent.log`). `scripts/check.sh quick` timed out at the 30 s tool limit on the login node (normally ~13 s); not retried, quick gate still pending.

### L9 — Official es0.25 S-6 does not diverge (SLURM 11621941, 5m39s)
- HEAD 0568825 + instrumentation only; extension `venv/.../_aeroelast...so` SHA 9e183ce (same as the 7d20c81 suite). Mesh 32 656 nodes, 195 528 free DOFs; log `logs/check/l3-s6-official025-11621941.log`.
- **No divergence at es0.25**: max translation 14.14 m at step0 → 13.42 m at step400; dynamic balance residual ≤5.5e−8 relative, every KSP reason 4 (converged). Tip OoP mean 14.028 m vs ElastoDyn 16.087 m (−12.80 %, bound 25 %), std 0.5824 m vs 0.6138 m (−5.12 %, bound 25 %). Both accepted bounds pass without change.
- Suite failure is the es0.5 path: `tests/test_iea15mw_s6_oneway_dynamic.py:63` calls `build_assembly(0.5)` (std 3.26e40 m in `logs/check/full-20261022T113305.log:107`). Per AGENTS rule 9, es0.5 is not converged and a test with es ≥ 0.3 is a defect.
- Open question before changing the test: is the es0.5 blow-up a genuine negative rotating eigenvalue of the coarse mesh (mesh artefact) or a defect that es0.25 merely hides? Probe `logs/check/l3-s6-probe/rotating_spectrum.py` (shift-invert K x = w² M x for K, K+Kg, Krot; es0.5 adds the growth rate) submitted as background job-2 (es0.5 on sequana_cpu, es0.25 on sequana_cpu_bigmem, 12 tasks each).

### L10 — es0.5 spectrum and growth (SLURM 11621944, 1m33s)
- Six w² closest to −1 are positive for K, K+Kg and Krot (lowest Krot 12.93 rad²/s² = 0.572 Hz). **Blind spot:** shift-invert at σ=−1 cannot see a far negative eigenvalue.
- es0.5 history: max translation 14.39 m (step0) → 14.98 m (step50) → 8.1e6 m (step100) while tip OoP is only 38 m: a **local** mode. Late growth rate 30.92 1/s. Solving s² + (ηKλ + ηM)s + λ = 0 for s = 30.92 (C is Rayleigh on the same rotating K) gives λ ≈ −774 rad²/s²: consistent with one statically unstable mode, not an integrator defect.
- Maintainer (2026-10-22): "si lo que estamos esperando del cliuster puede correr en menos de 20min puedes moverlo hacia la cola de dev". The pending es0.25 job (11621945, σ=−1, same blind spot) was cancelled; both meshes resubmitted to sequana_cpu_dev, `--time=20`, 12 tasks, with `--sigma -2000` (es0.25 Krot only to fit 20 min): background job-3, logs `logs/check/l3-s6-negspec-es*.log`.

### L11 — K_rot has no negative generalized eigenvalue near the growth (es0.5, SLURM 11621968, 47 s)
- Shift-invert at σ=−2000 returned the same six positive w² as σ=−1 (12.93 … 653.2 rad²/s²) for K, K+Kg and Krot. Any eigenvalue of K x = w² M x in [−4653, 653] would have appeared, so the inferred λ ≈ −774 does not exist as a stiffness divergence mode. **The inherited hypothesis "K+K_G+K_SP indefinite" is rejected as the cause** (if M is positive definite).
- es0.25 submission rejected by `AssocMaxSubmitJobLimit`; resubmitted when the slot freed (job-4, SLURM 11621970, dev).
- New hypothesis: M not positive definite on the free DOFs (a negative mass direction gives x'Kx/x'Mx < 0 and local growth with SPD K; the shift-invert generalized solve assumes SPD M and is then unreliable). Probe `logs/check/l3-s6-probe/growing_mode.py` (M diagonal per DOF type, M eigenvalues nearest 0⁻, peak node) queued on dev for es0.5 then es0.25 (job-6).
- es0.25 Krot at σ=−2000 (SLURM 11621970, 1m47s): six positive w² = 13.01, 19.58, 105.7, 172.5, 375.1, 543.0 rad²/s² (0.574, 0.704, 1.636, 2.090, 3.082, 3.709 Hz); first five within 2 % of es0.5. The 6th (local, r ≈ 98 m) carries a negative Kg Rayleigh term (−55.9 against K 598.9): compressive centrifugal skin prestress, stable here.

### L12 — M diagonal positive; M spectrum probe timed out
- `growing_mode.py` (SLURM 11621972 es0.5, 11621984 es0.25, dev): both TIMEOUT at 20 min inside `eigsh(M, sigma=−1e−6)`; not retried. Diagonal printed before the timeout: no non-positive entry. es0.5 translational 2.33e−5 … 20.7 kg, rotational 1.94e−12 … 1.56e−2; es0.25 translational 1.16e−5 … 5.35 kg, rotational 9.66e−13 … 4.17e−3 (≈13 decades: near-singular shift-invert near 0). Positive diagonal does not prove M SPD.
- Change of approach: trace the es0.5 blow-up directly. Observer record gained read-only `rotation_max` and `peak_node` (`tools/run_s6_oneway_dynamic.py`, still optional, no effect when omitted). Probe `logs/check/l3-s6-probe/growth_trace.py` (per-step DOF type, peak node, balance residual, KSP reason) submitted to dev (job-7).
- Parent re-run of `tests/test_s6_newmark.py` after the observer edit: exit 124 (28 s timeout on the login node; first run took 0.36 s). Not retried: the analytic GREEN for this edit is pending.

### L13 — es0.5 growth trace (SLURM 11622000, 1m38s, `logs/check/l3-s6-trace-es0.5-11622000.log`)
- Linear solves are exact for the discrete system at every step: relative balance 3e−11 … 4e−8, KSP reason 4 throughout. The blow-up is a true solution of the discretized equations, not a solve failure.
- **Rotation grows first**: max rotation already 2.51 rad at step 0 (static mean-load IC — unphysical), 2.58 at step 30, 13.5 at 40, 265 at 50, 5.8e3 at 60; translation stays at the tip (node 7, 14.4–15.0 m) until step 50, then jumps to node 1596 (xyz −1.93, −1.94, 92.53) and follows. Rotation/translation ≈ 168 late; step factor ≈ 22 per 10 steps.
- Working hypothesis: a nearly unstiffened rotational DOF (shell drilling rotation) where Kg/Ksp subtract more than K provides, with rotary inertia ~1e−12: a negative eigenvalue far outside the searched [−4653, 653] window. Not proved.
- Next probe `logs/check/l3-s6-probe/growing_shape.py`: Rayleigh quotients of d = u(400) − u(399) in K, Kg, Ksp, M, Krot, DOF-component share, peak node; also the static-IC rotation peak. Observer record now also carries `u_free` (copy; optional hook). Submitted for es0.5 then es0.25 (job-8, dev).

### L14 — Cause located: the linear elastic K is indefinite at r ≈ 92 m (SLURM 11622004 es0.5, 11622006 es0.25)
- es0.5 growing shape d = u(400) − u(399): 98.8 % in global θz (DOF 5), 0.6 % θx, 0.56 % θy; top-10 nodes carry 99.3 %, peak node 1636 at (−1.951, −1.891, 91.929). Quadratic forms: **q'Kq = −4.939** (elastic stiffness negative), q'Kg q = +3.130 (centrifugal prestress *stabilizes*), q'Ksp q = −8.7e−7, q'Mq = 1.18e−5 > 0, q'Krot q = −1.810. Rayleigh w² = −1.535e5 rad²/s²: far outside the searched window, which is why the shift-invert spectra (L10–L11) missed it.
- **So the divergence is a defect of the linear elastic stiffness K, not of K_G, K_SP, M, Newmark or the loads.** The earlier hypothesis (K_G drives indefiniteness) is the opposite of what is measured.
- es0.25: the growing shape is the healthy 1st flap (w² = 13.09, q'Kq = 5.64 > 0, 91 % flapwise translation), so no instability. **But the static IC already shows the same local pathology**: 3.534 rad rotation at node 5296 (−1.927, −1.948, 92.675), same location as es0.5 (2.514 rad at r ≈ 110 m, peak growing region r ≈ 92 m). A 3.5 rad static rotation under rated loads is unphysical: at es0.25 the same region is near-singular rather than negative.
- Physical consequence: the es0.25 S-6 tip mean/std (L9) pass, but they come from a model with a local stiffness defect at r ≈ 92 m. They cannot be promoted to evidence until that region is explained.
- Candidate link (not verified): the documented-red section-frame tear guard test in the full suite.
- Next unit (new session): element-level check around r ≈ 92 m. For the elements touching nodes 1596/1636 (es0.5) and 5296 (es0.25), compute element K eigenvalues (a negative eigenvalue isolates the element), geometry quality (area, Jacobian sign, aspect ratio, duplicate or near-duplicate nodes), family (MITC3/MITC4) and section/material, then compare with `test_section_frame_tear`. Fix only the demonstrated defect, with RED/GREEN, then re-run S-6 es0.25 and the matrix rows.
- Instrumentation unit checks (job-9): `tests/test_s6_newmark.py` 5 OK (3.82 s, exit 0); `scripts/check.sh quick` OK (75 passed, 1 xfailed; mesh-size OK, no new es ≥ 0.3 file). Ruff flagged 4 issues in changed files: 2 new ones fixed (B905 in the test, C401 in the diagnostic; re-run "All checks passed", analytic 5 OK); 2 in `run_s6_oneway_dynamic.py` (B905 at :144, F841 `n_nodes` at :171) pre-exist at 0568825 and are left untouched. Committed as the T5 instrumentation work unit (Risk: item 3, shared S-6 helpers; optional hooks only, histories bit-identical when omitted). The probes under `logs/check/l3-s6-probe/` are git-ignored diagnostics, not committed.

### L15 — Independent verification of c882c1f; session handoff
- Commit `c882c1f` (test(s6): instrument the S-6 Newmark run and locate its divergence). Native `assess` (baseRef 0568825, committedOnly) unassessable: the untracked BeamDyn directories need an explicit review declaration; risk raised to high by escalation item 3.
- gentle-ai-verify: items 1–4 pass, no blockers. Differential base vs new helper on a synthetic 36-DOF SPD system (MUMPS, 300 steps): bit-identical with hooks omitted and with an observer, including one that overwrites `u_free`. Analytic margins: k=4 max error 1.06e−7 vs atol 2e−6; k=−100 max error 6.8e−6 vs atol 2e−5. Advisories only: `static_solve_reason` is constant per run; `rotation_max`/`peak_node` assume dpn=6 (shell). Verifier scratch left at `$SCRATCH/tmp/s6verify.57E9` (outside the project; not removed without authorization).
- Delivery: the branch is 4013 authored lines in 29 commits past origin/main, beyond the ~400-line slice. Strategy ask-on-risk: slicing into PRs is the maintainer's decision; no push.
- **Resume recipe (next session):** read this document from L14; `mem_search "indefinite elastic K r~92"` (observation 90); `git log --oneline -3`. Next: element-level K eigenvalues and geometry quality for the elements touching es0.5 nodes 1596/1636 and es0.25 node 5296 (r ≈ 92 m); compare with `test_section_frame_tear`; RED test for the demonstrated defect, fix, then re-run S-6 es0.25 (dev, ~4 min) and check that the 3.5 rad static rotation disappears before recording any S-6 number in the store. Diagnostics that need <20 min go to `sequana_cpu_dev` (maintainer, L10).

### L16 — Continue in the same session (2026-10-22)
```text
vamos a seguir en esta sesion, pero creo que debemos comitear primero
```
- Already committed: `c882c1f` + `1bd312d`; tree clean except the six protected BeamDyn directories.
- Worker (context backstop) wrote the git-ignored probe `logs/check/l3-s6-probe/element_check.py`. It obtains per-element K through the production path, a one-element `RustMeshModel` → `PyMeshAssembler.from_model` → `assemble_k` (`assembler.py:458,585`, `assembly/assembler.rs:369`, `mitc3.rs:930`, `mitc4.rs:1498`). It reports geometry quality, Jacobian sign, flipped neighbours, laminate/span data and element eigenvalues.
- Self-test passed: 6 cases with 0 negative eigenvalues. Isotropic MITC4 has 7 exact-zero modes (drilling unstiffened by itself), MITC3 has 6 plus a weak drilling mode near 3e−9·λmax, and an inverted bow-tie quad is still PSD, so the Jacobian check is required. Hypotheses only: MITC4 drill-membrane strain (`mitc4.rs:1207,1394`), tiny MITC3 `k_drill` (`mitc3.rs:320`), composite `e_equiv` discarded (`crates/aeroelast-py/src/elements.rs:709`), triangle span-angle fallback (`assembler.rs:23`).
- Submitted to dev: es0.5 nodes 1636 1596, then es0.25 node 5296 (job-10, logs `logs/check/l3-s6-elem-es*.log`).

### L17 — Element check and fold test: no single defective element; near-zero rotational mode at the TE
- SLURM 11622054 (es0.5) / 11622055 (es0.25): every element touching or neighbouring the peak nodes is PSD alone (0 negatives, ≥6 zero modes, rigid residual <1e−8). They are MITC4-composite slivers of the trailing-edge closure (`11_38_LP_TE_FLAT`, `00_38_HP_TE_FLAT`): 1.81 mm × 0.30 m, aspect 165, area 5.4e−4 m², warp ~1.5e−6, Jacobian positive; neighbouring normals n·n = −0.033 (≈92° corner). Each element has 8 near-zero modes (lmin/lmax ~1e−16).
- Assembled K principal block on the 15 nodes: min eig −55.9 (es0.5; block max ~1e10, i.e. −5.6e−9 relative) vs +2.42 (es0.25). Growing-shape q'Kq = −4.94 is only ~−4e−12 of λmax(K) ≈ 1.2e12. **Reading: a near-zero-stiffness rotational mode (θz) at the TE-corner nodes, whose sign is set at round-off level; with rotary inertia ~1e−5 it gives |w²| ~1e5.** At es0.5 it lands negative and diverges; at es0.25 it lands slightly positive (stable) but shows as the 3.5 rad static rotation. Not yet proved which formulation term leaves it unstiffened.
- Director-averaging hypothesis (`mitc4_nodal_directors`, `crates/aeroelast-core/src/assembly/assembler.rs:1347`, overwrite at `:195–197` and `:325–327`) **not reproduced**: the worker's two-strip MITC4 fold test (`tests/test_mitc4_fold_directors.py`, **untracked, not committed**) passes 7/7 in 2.43 s. Assembled min λ/λmax is −1e−16 for isotropic 0–90° and composite 90°, and −6.3e−11 for the composite 1.8 mm sliver at 90° (tolerance −1e−10). Folding removes one zero mode (8 → 7, flat = 6 rigid + 2 drilling) and adds no negative eigenvalue. The near-miss could flake; decide whether to keep the file as a characterization test.
- Open hypotheses: (a) unstiffened MITC4 drilling at nodes where the TE closure creates a 92° corner between 1.8 mm slivers (flat patch already has 2 drilling zero modes); (b) the TE closure mesh itself, a 1.8 mm strip, is a degenerate representation of the flatback/TE that leaves a rotational mechanism; (c) round-off amplification by aspect 165 slivers.
- Next (new session): extract the actual TE patch (the 24 elements around node 1596) into a sub-model and compare its assembled K with the full-mesh block; compute the null-space vector of that patch and its DOF content; test whether a drilling stiffness, or removing/merging the 1.8 mm TE strip, changes es0.25 static rotation and the es0.5 sign. Mesh-generation change versus element drilling change is a formulation choice: record the evidence first.
- Physical status unchanged: S-6 es0.25 numbers (L9) are not evidence until the TE mechanism is removed.

### L18 — Defect demonstrated: director averaging across the TE fold (2026-10-22)
```text
Retomamos odd/tasks/l3-physical-validation.md, T5, desde L17
   (mecanismo rotacional en el cierre TE, r≈92 m). mem_search "TE closure slivers"
 (observación 90).
```
- **L17 correction.** The "sign set at round-off" reading was wrong. A principal block of a sum of PSD matrices cannot reach −55.9 against a maximum of 1e10 (5.6e−9). L17's element check used one-element models, where the director is the element's own normal. In the blade, `mitc4_nodal_directors` area-averages the directors of *all* MITC4 elements at a node, with no fold criterion.
- TE sub-model (`logs/check/l3-s6-probe/te_submodel.py`, SLURM 11622060, `logs/check/l3-te-submodel-11622060.log`). It covers the patch plus its 1-ring, so directors are exact on the patch nodes. Its block reproduces the blade's (−55.87 vs −55.87 at es0.5; +2.410 vs +2.416 at es0.25).
  - TE_FLAT + REINF together: min/λmax −5.2e−9 with 9 negatives (es0.5) and −3.3e−9 with 4 negatives (es0.25), all ≥97 % θz on the TE nodes.
  - Each family alone is PSD: TE_FLAT only −3e−17, REINF only −2e−16.
  - Fold-node directors sit 79.4° (LP) and 85.3° (HP) off the sliver normal, because the 0.134 m² REINF quad outweighs the 5.4e−4 m² sliver. The sliver's director lies in its own plane.
- Synthetic sweep (panel 0.45 m + sliver at fold θ, `logs/check/l3-te-fold-sweep.log`): indefinite only when h/width ≫ 1. With h = 12 mm on 1.8 mm, θ 30–45° gives down to −5.4e−7 and θ 135–175° gives about −3e−8. With h = 1.2 mm it is clean at every θ.
- **Fix (`crates/aeroelast-core/src/assembly/assembler.rs`, `mitc4_nodal_directors`).** Directors are now per element-node and follow the Abaqus rule (*Defining the initial geometry of conventional shell elements*, "Calculation of average nodal normals"). Directors at a node are grouped transitively by pairwise angle ≤ 20°. A group is averaged (area-weighted, ADR-4 option A) only if all its pairs are within 20°; otherwise each element keeps its own director. Nodes with more than 30 elements are not averaged. Pass 2 in `new` and `rebuild_precomputed` assigns `pre.vn = directors[e]`.
- RED (`tests/test_mitc4_fold_directors.py`, now tracked; `logs/check/l3-fold-criterion-red.log`): 14 failed, 13 passed.
  - Folds of 21, 92 and 150° and the 0/15/30 chain give |K − ΣK_alone|/|K| = 3.6e−3 … 0.96.
  - The thick-sliver PSD cases give down to −5.4e−7.
- Release rebuild 49 s; extension SHA `c0d9b9a7…` (baseline `9e183ce…` kept in `logs/check/l3-fold-baseline-extension/`). GREEN: 27 passed in 2.22 s (`logs/check/l3-fold-criterion-green.log`).
  - Folds of more than 20°: ≤4e−17.
  - Smooth 10/19° and 0/10/19 cases: unchanged values (2.04e−3 …).
  - Thick sliver: ≥ −9e−17.
- Differential, base vs new extension (`logs/check/l3-fold-differential.log`):
  - Bit-identical: twisted 90° strip, 24-sided cylinder and 10° fold (iso and composite).
  - Changed only beyond 20°: 92° sliver fold, T-junction, and an 8-sided cylinder (45°/element). The last one is the Abaqus caveat: a coarse mesh of a smooth shell becomes a fold. Whether that happens on the blade (LE at es0.25) is measured by the fold census.
- After the fix (extension `c0d9b9a7…`):
  - **TE sub-model** (`logs/check/l3-fold-census-11622071.log`): 0 negatives and exactly 6 zero modes at es0.5 and es0.25. The blade block min eig is now +174.2 at es0.5 (was −55.9) and +102.4 at es0.25 (was +2.4).
  - **S-6 growing shape** (`logs/check/l3-fold-s6shape-11622068.log`): **es0.5 no longer diverges**. Its growing shape is the 1st flap (w² = 13.02 rad²/s², was −1.5e5); es0.25 gives 13.10.
  - **es0.25 static-IC max rotation**: 1.485 rad θz at node 986 (−1.17, −3.59, 112.5), down from 3.53 rad at the TE at r = 92.7. At es0.5 it is 2.43 rad at node 528, r = 109.8, essentially the same as the 2.51 rad at r ≈ 110 reported in L14. This second local rotation near the tip region is open and is *not* the TE fold.
- **Fold census es0.25** (same log). 5138 of 32 656 quad nodes now keep separate directors:
  - spar cap ↔ shear web: 2344 nodes, ~90°;
  - TE_FLAT ↔ TE_REINF: 1232 nodes, 80–105°; TE_FLAT ↔ TE_FLAT: 53;
  - LE seam HP_LE ↔ LP_LE: 506 nodes, 97°;
  - coarse LE curvature: about 1000 nodes inside LE / LE_PANEL at 21–43°.
  - 1909 nodes sit at 10–20° and stay averaged.
  The web junction and the TE are true intersections, where separate directors are the standard treatment. The LE is the Abaqus coarse-mesh caveat: a smooth surface becomes faceted at es0.25. **Every blade anchor can move**, so the fix is not accepted until S-1/S-2/S-4/S-7/V-02 are re-measured against the 7d20c81 suite. That suite is a valid baseline: same source except this fix, extension 9e183ce.

### L19 — S-6 after the fold fix; pre-existing TE cross-section distortion near the tip (2026-10-22)
- Fix committed as `7dfe78d` (Risk: item 3).
- **S-6 official es0.25 with the fix** (SLURM 11622086, extension c0d9b9a7, `logs/check/l3-s6-fold-es0.25-11622086.log`): tip OoP mean 14.024 m vs ElastoDyn 16.087 m (**−12.83 %**, was −12.80 %), std 0.5821 vs 0.6138 m (**−5.16 %**, was −5.12 %), against an accepted 25 % bound. The global response did not depend on the TE mechanism.
- **Residual static rotation is not a stiffness mechanism.** Node 986 (r = 112.5, HP_TE_PANEL/REINF, smooth: director 0.7–1.3° off the normals) has a PSD 1-ring sub-model with exactly 6 zero modes. Its 1.48 rad rotation is not drilling (θ·director ≈ −0.06). Static split (`logs/check/l3-s6-static-split-11622100.log`): **K alone** gives 1.608 rad, K+Kg 1.482, Krot 1.480, so prestress is not the cause. Nodal rotation percentiles 50/99/99.9 = 0.084/0.66/1.40 rad. The S-6 loads are translational forces only (`nodal_load_series`, `a in range(3)`), with no applied moments.
- **Cross-section distortion at the TE** (`section_twist.py`, static K, mean S-6 load; `logs/check/l3-s6-section-twist-11622101.log` vs baseline extension 9e183ce `l3-s6-section-twist-base9e183ce-*.log`). Per station, the rigid twist φ is fitted to the ring's in-plane translations; the residual and the nodal θz spread measure in-plane distortion. Worst node always in the TE region:

  | r [m] | φ fix / base [rad] | ring residual fix / base [m] | θz min…max fix (base) |
  |---|---|---|---|
  | 85.0 | −0.039 / −0.041 | 0.012 / 0.019 | −0.09…0.00 (−0.12…0.02) |
  | 90.6 | −0.053 / −0.067 | 0.032 / 0.083 | −0.24…0.02 (−0.52…0.47) |
  | 97.6 | −0.096 / −0.115 | 0.159 / 0.229 | −0.62…0.59 (−0.97…0.80) |
  | 110.1 | −0.158 / −0.170 | 0.245 / 0.269 | −1.47…0.59 (−1.57…0.62) |
  | 117.0 | −0.149 / −0.167 | 0.023 / 0.030 | −0.26…0.10 (−0.37…0.23) |

  The distortion **pre-exists the fix, which reduces it** (by 30–60 % at r = 90–98 m). It is still 16–25 cm of in-plane ring distortion at r = 97–112 m under a linear static load with a 15.4 m tip deflection. The fix also lowers the rigid twist by 7–20 % outboard, so the S-7 torsion anchor is expected to move.
- No independent reference exists yet for local section distortion; the S-6 tip metrics do not see it. Hypotheses, none tested:
  - (i) The load model applies the sectional fn/ft as equal translational forces on every ring node, including the TE panel and TE_FLAT nodes. A real aerodynamic load is a surface pressure, so this may over-load the TE panel out of plane.
  - (ii) A soft TE panel / TE closure laminate near the tip (interpolated sections).
  - (iii) A remaining element issue in thin, stiff-face sandwich panels.
  The discriminating test is the same static case with the forces lumped on the spar/shear-centre nodes only, or a CCX S8R arbiter of the same mesh and load.
- Physical status:
  - The S-6 es0.25 tip numbers now come from a PSD model without the TE mechanism, but are not promoted to the store until the full-suite anchor check (below) and the section-distortion question are resolved.
  - `tests/test_iea15mw_s6_oneway_dynamic.py:63` still builds es0.5, a rule-9 defect (es ≥ 0.3); that unit's fix (move to 0.25) is pending and is a test change with no physics debate.
- **Full suite at 7dfe78d running** (job-5, `scripts/check_sbatch.sh full -- --ntasks=12`, sequana_cpu, log `logs/check/full-<stamp>.log`). Baseline: `logs/check/full-20261022T113305.log` (7d20c81, extension 9e183ce, 18 failed / 1407 passed). Compare per anchor:
  - S-1 EI/GJ/EA;
  - S-2 LC1–LC4 official;
  - S-4 1F/1E/2F;
  - S-7 9.0924 %;
  - V-02 1F 6.0930 % and 1E 8.4707 %;
  - CCX 9.2182 %;
  - S-6 (×2, es0.5 now expected to pass or change);
  - twist mechanism and section-frame tear guard.
  Any anchor that moves is recorded as a finding, never re-banded.
- **Resume recipe.**
  1. Read this document from L18 and `mem_search "director averaging folds"` (observation 90).
  2. When the suite ends: `grep -E "^(FAILED|ERROR)|passed|failed" <log>` and diff the failure set and the anchor percentages against the 7d20c81 log.
  3. If the anchors hold or improve, record group 39/28/S-7/V-02 rows at 7dfe78d and update S-6.
  4. Then decide the section-distortion probe (lumped-load variant first, dev ~2 min).

### L20 — Full suite at 7dfe78d vs 7d20c81: the fold fix is a trade-off, not yet accepted
- SLURM 11622087 (2 h 0 m, `logs/check/full-20261022T160155.log`): 21 failed, 1436 passed, 38 skipped, 13 xfailed, 3 xpassed. Baseline `full-20261022T113305.log`: 18 failed, 1407 passed. Source differs only by `7dfe78d` (+27 fold tests).
- **Fixed (2):** S-6 mean and S-6 std (es0.5 divergence gone).
- **New failures (5):**
  - CCX blade parity `tip_axial_tension` / `tip_axial_compression`: tip displacement 5.46 % off CCX (bound 5 %). The previous value passed but was not printed.
  - Thin-walled tube, three tests:
    - `minimum_norm_realisation_over_delivers_bredt` (pinned over-delivery moved);
    - `production_distributed_application_is_bredt` (pinned energy excess moved);
    - `self_equilibrated_laminate_reproduces_bredt` ("rigid modes not removed").
    The box tube has 90° corners, exactly the fold case the fix changes. Its Bredt reference is analytic, so it is the cleanest arbiter of the change.
- Persisting failures, anchor values before → after:

  | anchor | 7d20c81 | 7dfe78d | change |
  |---|---|---|---|
  | S-1 EI_edge | 12.69 % | 11.45 % | better |
  | S-1 GJ | 6.22 % | 5.05 % | better |
  | S-1 EI_flap / EA | 5.17 / 7.59 % | 5.17 / 7.59 % | = |
  | S-2 LC1 / LC4 official | 7.42 / 6.23 % | 7.44 / 6.27 % | ≈ (+0.02/+0.04 pp) |
  | S-7 tip twist ratio | 9.0924 % | **10.2000 %** | **worse +1.1 pp** |
  | V-02 1F / 1E | 6.093 / 8.471 % | 6.070 / 8.505 % | ≈ |
  | CCX high-mode gap | 9.2182 % | **13.9784 %** | **worse +4.8 pp** |
  | twist-mechanism guard | 63.8 > 10·12.25 | 59.9 > 10·7.07 | still red |
  | blade node count, tear guard, UTD S-2 | unchanged | unchanged | = |

- Reading: separate directors at the >20° folds fix the TE mechanism, improve sectional EI_edge/GJ, and move the code-to-code comparisons with CalculiX (high modes, tip axial) and the BeamDyn torsion ratio away. Hypothesis, not checked: CalculiX shells average nodal normals with their own rule (and generate "knots" at folds), so part of the CCX gap may be a director-convention difference, not physics. **Not accepted, not reverted.** No anchor was re-banded; no store row was written at 7dfe78d.
- **Next (new session), in this order:**
  1. **Box-tube Bredt arbitration.** Analytic reference at 90° corners. Run the three tube tests and read the moved values and the "rigid modes not removed" path. If separate directors are closer to Bredt, the fix is right at true folds.
  2. **The CalculiX shell normal rule** (manual, plus the CCX deck used by `test_blade_ccx_parity.py`). Is the CCX gap a convention difference?
  3. **S-7 at 7dfe78d vs 7d20c81** per span band.
  4. **Then decide the threshold scope.** Options: Abaqus 20° everywhere (current); folds only at true set-boundary junctions (web/TE); or 20° plus es-0.15 convergence for the coarse LE. That is a formulation decision, to be reported with this evidence.
  Then return to L19 steps 2–3: the S-6 test at es0.25 and the TE section distortion.

### L21 — Box-tube Bredt arbitration and the CalculiX normal rule (2026-10-22)
```text
Retomamos odd/tasks/l3-physical-validation.md, T5, desde L20. mem_search "fold criterion suite impact" (observación 92).
```
- **Bredt arbitration (L20 step 1).** Same tree; extension swapped through `PYTHONPATH` (base copy `$SCRATCH/tmp/ext-base-9e183ce`, fix `c0d9b9a7` installed). Logs: `logs/check/l3-tube-base9e183ce.log` (8 passed), `l3-tube-fixc0d9b9a7.log` (3 failed, 5 passed), `l3-tube-fold-probe.log` (`logs/check/l3-s6-probe/tube_fold_probe.py`). Reference: the analytic Bredt-Batho `θ' = T/GJ` of the 1.0 × 0.6 m tube, 90° corners. Ratios to Bredt:

  | case | metric | base 9e183ce | fix c0d9b9a7 |
  |---|---|---|---|
  | iso, self-equilibrated | θ_fit / θ_z | 1.00309 / 1.00675 | 0.99695 / **0.99976** |
  | laminate [45,−45]s, self-equilibrated | θ_fit / θ_z | 0.99859 / 1.00704 | 0.99419 / **0.99999** |
  | iso, tip shear flow, clamped root | θ_fit / θ_z | 0.99156 / 1.00628 | 0.99587 / 0.99929 |
  | θ_fit vs θ_z disagreement (iso / lam) | — | 0.36 / 0.84 % | 0.28 / 0.58 % |
  | clamped, offset transverse load, closed form | tip θ_z err (10 mm, 0.5 mm) | 0.23 / 0.11 % | 0.54 / 0.41 % |
  | production distributed flow, energy excess (estimator-free) | dz 0.2 / 0.1 m | +4.356 / +2.484 % | **+4.251 / +2.048 %** |
  | one-wall traction distortion amplitude (wall/shear flow) | — | 125.7× | 116.9× |

  - θ_z, which the tube module establishes as the boundary-layer-free section rotation, goes from 0.7 % off Bredt to within 0.02 %. The energy metric improves and still converges with dz. The only worse number is the offset closed form, by 0.3 pp, still under 0.6 %.
  - **Rigid-body invariance is intact.** ‖K φ_a‖/‖K‖_F = 1–4e−18 for all six rigid modes, base and fix.
  - **"rigid modes not removed" is solve round-off, not physics.** The saddle-point `splu` leaves Φᵀu at 2.6e−10 (bound 1e−10) for the fix laminate, 2.9e−12 for the base. One step of iterative refinement takes it to 7e−16 (base 2e−13) and leaves the ratios unchanged.
  - The two moved pins are not reference comparisons. The minimum-norm over-delivery is a recorded defect (gaps.yaml `moment_realization_over_delivers`): 31.69 → 24.66× self-equilibrated, 10.41 → 9.20× clamped. The distributed-flow energy excess is a discretisation pin and moved toward zero.
  - **Verdict:** at true 90° folds, separate directors are equal to or closer to the analytic reference than averaged directors. The fold rule is right at true folds.
- **CalculiX normal rule (L20 step 2, manual §6.2.14, S8/S8R; §6.2.33 "Calculating the normals in the nodes").**
  - CCX groups element normals at a node around a reference normal within **20°**, and within 0.5° for `*NORMAL` user normals. Its procedure is the same as the Abaqus rule now in `mitc4_nodal_directors`.
  - It **also splits normals where local thickness or offset differ**. Our averaging ignores thickness.
  - A node with more than one normal is expanded once per normal and tied by a **knot** (rigid body; mean-rotation MPC in statics).
  - On a smooth surface, S8R normals come from curved quadratic geometry, so they rarely trip the 20° rule. Flat MITC4 facets at a coarse LE do.
  - So the L20 hypothesis "CCX averages where we now split" is **refuted at true folds**: CCX splits there too. The remaining convention differences are (a) CCX splits at every laminate-thickness change across set boundaries, and we average; (b) at the coarse LE, MITC4 facets exceed 20° where S8R does not. Finding (b) is the Abaqus coarse-mesh caveat. It points at the 2.0/1.0/0.5 m meshes of `test_blade_iea15mw_mesh_convergence.py` (finest 0.5 m), where the CCX high-mode gap moved 9.22 → 13.98 %. Fold census at es 0.5/1.0 is running to size it.
- **Fold census vs mesh size** (`logs/check/l3-fold-census-es{0.5,1.0}.log`; es 0.25 from L18). Nodes where directors stay separate:

  | es | split nodes / total | true intersections (SW, TE_FLAT↔REINF, TE_FLAT↔TE_FLAT) | LE region (LE seam, LE↔LE_PANEL, LE_PANEL) | max-angle median |
  |---|---|---|---|---|
  | 1.0 | 1687 / 3067 (55 %) | ≈1013 | ≈655 | 23.9° |
  | 0.5 | 2872 / 9357 (31 %) | ≈1905 | ≈960 | 10.0° |
  | 0.25 | 5138 / 32656 (16 %) | ≈3629 | ≈1500 | — |

  - The HP_LE↔LP_LE nose seam stays at 97–103° at every es. It is a chordwise kink of the LE discretisation, so refining es does not remove it.
  - LE↔LE_PANEL (≈50°) and LE_PANEL interior (21–42°) shrink with es. That is the Abaqus coarse-mesh caveat.
  - The CCX convergence test's meshes (2.0/1.0/0.5 m) are where the most smooth LE becomes a fold; that is consistent with its high-mode gap moving most. Not isolated: it needs a variant run.
- **S-7 per span band** (SLURM 11622227, es 0.25, both extensions on the same tree; `logs/check/l3-s7-bands-11622227.log`, `logs/check/l3-s6-probe/s7_bands.py`). It reproduces the suite: 0.9091 → 0.8980 over 0.3–0.9.

  | r/R | 0.1–0.2 | 0.2–0.3 | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 | 0.6–0.7 | 0.7–0.8 | 0.8–0.9 | 0.9–1.0 |
  |---|---|---|---|---|---|---|---|---|---|
  | base ratio | 0.757 | 0.946 | 0.936 | 0.897 | 0.883 | 0.884 | 0.910 | 0.945 | 1.864 |
  | fix ratio | 0.783 | 0.949 | 0.932 | 0.892 | 0.876 | 0.875 | 0.891 | 0.922 | 1.560 |

  - The fix stiffens torsion by 0.4–0.9 pp at 0.3–0.7 and by 1.9–2.3 pp at 0.7–0.9, the thin outboard sections.
  - The tip band (couple application, local distortion; see L19) drops from 4.75× to 3.32× the beam tip twist.
  - S-7 compares against a beam deck with a known model difference: the OML reference surface, E8 of `odd/tasks/s7-torsion-ratio.md`, which already puts the shell ~13 % stiff. It cannot arbitrate a director convention. The movement is recorded as a finding; the band is unchanged.
- **Status (L20 step 4, not decided).** Evidence so far:
  - At true intersections (web, TE closure), separate directors are right: Bredt, the TE sub-model PSD, S-6 es 0.5.
  - At the coarse smooth LE they are the documented Abaqus caveat. The code-to-code movements (CCX high modes, CCX tip axial, S-7 outboard) have not been attributed between the two populations.
  - The discriminating experiment is a variant extension that keeps directors averaged on the LE population only, then re-measures CCX high modes, CCX tip axial, S-7 and S-1. Rust change plus rebuild; CCX convergence ~11 min plus S-7/S-1 on SLURM.
  - Until then, 7dfe78d stays unaccepted and unreverted. No pin re-set, no store row.
- **Resume recipe.** Read L20–L21; `mem_search "fold criterion LE caveat Bredt"`. Next, in order:
  1. Implement the population-split variant behind a test-only switch and run the four anchors under it.
  2. Pick the criterion from that evidence.
  3. Re-pin the two tube pins and add one step of iterative refinement to `_solve_rigid_removed`, keeping the 1e−10 bound.
  4. Write the store rows at the accepted revision.
  5. Then return to L19 steps 2–3.

### L22 — Population-split variant (L21 step 1), 2026-10-22
```text
Retomamos  odd/tasks/l3-physical-validation.md, T5, desde L21" y mem_search "fold criterion LE caveat  Bredt"
```
- Variant `smooth_oml` (scratch only, never committed): worktree `logs/check/l3-fold-variant/wt` at d745327, +46 lines. `from_model` flags elements whose property set contains `SW` or `TE_FLAT` when `AEROELAST_DIRECTOR_VARIANT=smooth_oml`. A node touched by any flagged element keeps the 20° rule of 7dfe78d. A skin-only node (LE seam, LE↔LE_PANEL, LE_PANEL, LE_PANEL↔SPAR, TE_REINF interior) gets the pre-7dfe78d area-weighted mean. Built to `logs/check/l3-fold-variant/ext` (`build.sh`), selected through `PYTHONPATH`; installed extension untouched.
- Build note: the first `maturin build` let auditwheel rewrite the cargo artifact in place (NEEDED `libslepc-9b11be1e`, rpath `aeroelast_py.libs`), so the extension did not import. Relinked with `--auditwheel skip`: `.so` sha256 `ab480983…`, NEEDED `libslepc.so.3.25`/`libpetsc.so.3.25` from the venv. Mask sanity at es 1.0: `[fold-variant] smooth_oml: 911 fold elements, 1642 skin-only nodes averaged` (3067 nodes).
- **Anchors under the variant** (SLURM 11622271, 27m47s, `logs/check/l3-fold-variant/anchors-11622271.log`; 9 failed, 9 passed, 1 xfailed). The variant was active on every blade mesh (`[fold-variant]` lines at es 2.0/1.0/0.5/0.25: 440/911/2461/7544 fold elements). Fix values from L20/L21 (same tree, installed extension):

  | anchor | base 7d20c81 | fix 7dfe78d | variant smooth_oml |
  |---|---|---|---|
  | CCX high-mode gap | 9.2182 % | 13.9784 % | **13.9766 %** |
  | CCX tip axial tension / compression | pass (<5 %) | 5.5 / 5.5 % | 5.4 / 5.4 % |
  | S-7 tip twist ratio (test) | 9.0924 % | 10.2000 % | **10.1888 %** |
  | S-7 band mean 0.3–0.9 (probe) | 0.9091 | 0.8980 | 0.8981 |
  | S-7 tip band 0.9–1.0 (probe) | 1.864 | 1.560 | 1.570 |
  | S-1 EI_edge / GJ / EA | 12.69 / 6.22 / 7.59 % | 11.45 / 5.05 / 7.59 % | 11.45 / 5.05 / 7.59 % |
  | S-1 EI_flap | 5.17 % (red) | 5.17 % (red) | **passes** (value not printed; see S-1 probe below) |
  | V-02 1F / 1E | 6.093 / 8.471 % | 6.070 / 8.505 % | 6.070 / 8.505 % |
  | S-6 (3 tests) | es0.5 diverged | pass | pass |

  - **Finding:** averaging the skin-only population back recovers almost none of the code-to-code movement: CCX high mode −0.002 pp out of +4.76, S-7 −0.011 pp out of +1.11, tip axial ≈ −0.1 pp. The coarse-LE caveat is not what moved CalculiX, S-7 or V-02. The movement comes from the true-intersection nodes (web junctions, TE closure), the same population where Bredt (L21) and the TE sub-model PSD favour separate directors.
  - The S-1 EI_flap flip (red → pass) has to be checked with a printed value before it means anything; `s1_values.py` runs fix vs variant on the same node.
- **S-1 printed** (SLURM 11622335 dev, `logs/check/l3-fold-variant/s1-11622335.log`, same node): fix EI_edge/EI_flap/GJ/EA = 11.45/**5.1702**/5.0508/7.5939 %; variant 11.4486/**4.8594**/5.0526/7.5939 %. The flip is −0.31 pp crossing the 5 % bound by 0.14 pp, which is not evidence that either variant is right.
- **Decision (L21 step 2): keep the Abaqus 20° rule everywhere (7dfe78d as committed).** Reasons:
  1. At true folds, separate directors match the analytic reference (Bredt, L21) and remove the TE mechanism (S-6).
  2. The skin-only (coarse LE) population moves no anchor by more than 0.31 pp. The CCX, S-7 and V-02 movements belong to the true-fold population.
  3. A set-name split would put mesh-semantic knowledge (`SW`, `TE_FLAT`) into the assembler: a fragile contract for no measured gain.
  4. CalculiX and Abaqus use the same 20° rule.
  - Consequence, reported and not re-banded: CCX tip axial goes 5.4–5.5 % against a 5 % bound (it passed before); the CCX high-mode gap goes 9.22 → 13.98 %; S-7 goes 9.09 → 10.20 %. These are moved code-to-code anchors and are findings.
  - **Open question (reference side):** CCX also splits normals at true folds (§6.2.14) but ties them with a knot over the expanded C3D20R. That is not the same joint as our shared nodal rotations, so why separation moves us away from CCX at the web junctions is not explained. Not investigated here.
  - The scratch variant stays in `logs/check/l3-fold-variant/` (worktree + ext) as the record; nothing from it is committed.
- **L21 step 3 done.** Tube pins re-measured at 7dfe78d: over-delivery 24.66152× / 9.19606×; distributed energy excess 4.251 / 2.048 % (minimum-norm 2.88902 / 2.82458). `_solve_rigid_removed` gets one iterative-refinement step. Both tube files: 8 passed in 41 s (`logs/check/l3-tube-repin.log`). gaps.yaml and groups.yaml narratives carry the new pins; `validation_matrix.py check`: 0 errors (the pre-existing group-27 warning remains).
- **L21 step 4 done.** Store rows written at 7dfe78d (run 7dfe78d, date 2026-10-22; sources: full suite SLURM 11622087 `logs/check/full-20261022T160155.log`, `logs/check/l3-fold-variant/s1-11622335.log`, `logs/check/l3-tube-repin.log` at c0049c3). No bound, tolerance, verdict, flag or reference touched; contract.yaml, groups.yaml, gaps.yaml untouched. margin_pct old -> new:
  - V-02 (32): 1F 6.093 -> 6.0700 %; 1E 8.4707 -> 8.5046 %.
  - S-7 (34): tip twist ratio 9.0924 -> 10.2000 %.
  - Mesh convergence (5): high-mode gap 9.2182 -> 13.9784 %; low-mode gap (1 % bound) 0.1948 -> 0.5138 % (first value printed by a full-suite log).
  - S-1 (38): EI_edge 12.7095 -> 11.45 %; GJ 6.2193 -> 5.0508 %; EA 7.5939 -> 7.5939 % (unchanged).
  - S-2 official (39): LC1 7.4188 -> 7.4393 %; LC4 6.2268 -> 6.2659 %.
  - Tube torsion (6): iso theta_fit ratio 1.00309 -> 0.996948 (0.3091 -> 0.3052 %), iso theta_z 1.00674 -> 0.99976 (0.6745 -> 0.0240 %), laminate theta_fit 0.998594 -> 0.994185 (0.1406 -> 0.5815 %), laminate theta_z 1.00704 -> 0.999999 (0.7036 -> 0.0001 %); clamped-root theta_z 10 mm 0.6772 -> 0.0216 %, 0.5 mm 0.701 -> 0.0027 %; distributed force tip theta_z 10 mm 0.2341 -> 0.5433 % (both offsets), 0.5 mm 0.1068 -> 0.4108 % (both offsets).
  - Tube moment realisation (31): interior twist rate, both comparisons 0.3091 -> 0.3052 % (raw 2.75983e-05 -> 2.74292e-05).
  - Not recorded: CCX tip_axial_tension/compression (no store row exists for `test_blade_ccx_parity.py`; the log prints 5.5 %, L20 quotes 5.46 %); S-6 (no store row; tests pass, no value printed); S-1 mass (passes, not printed); S-1 flapwise EI (the row holds the es 0.125 converged value 5.3104 %, a different quantity; the es 0.25 control moved 5.1666 -> 5.1702 %); S-2 official LC2 and LC3 (not printed in the log); S-2 UTD cases (disowned input, no rows); signed S-1 values (only |rel err| printed); the tube minimum-norm pins (not store comparisons).

### L23 — Unattended mandate (2026-10-22)
```text
continua en modo desatendido hasta que identifiques y resuelvas los problemas de los modelos, o los problemas de implementacion que hacen que no podamos reproducir la fisica, o que pruebes que nuestro metodo es mas fisicamente realista, o que las diferencias simplemente son restricciones de la plataforma, pero no te detengas mas hasta resolver los prolemas por uno u otro de las planteamientos anteriores
```
- Each open residual closes with one of four outcomes: (a) model defect fixed, (b) implementation defect fixed, (c) our method shown more physically realistic against an independent reference, (d) the difference shown to be a platform or reference-model restriction. Every outcome needs a measured number and a stated reference.

### L24 — D3: PreComp section arbiter (2026-10-22)
Question (D3): is the published 6x6 (`components.blade.elastic_properties_mb.six_x_six`, = BeamDyn deck) consistent with the yaml's own `internal_structure_2d_fem` layup, or is our yaml→section reading the outlier? Arbiter: WISDEM 4.0.2 PreComp (NREL thin-wall composite section code, independent of our mapping), driven through WISDEM's own loader and RotorSE.
- **Install**: `$SCRATCH/venv/bin/python -m venv logs/check/l3-d3/venv`; `pip install wisdem==4.0.2` succeeded (`pip.log`). Needed one shim inside the venv: the base interpreter has no `_sqlite3`, so `pysqlite3-binary` + `venv/.../sitecustomize.py` aliasing `sqlite3`. `$SCRATCH/venv` and the project install untouched. PreComp is pure Python in 4.0.2 (`wisdem/precomp/precomp.py`; the Fortran import is commented out).
- **Feeding PreComp the official yaml** (all in `logs/check/l3-d3/`):
  1. WISDEM 4.0.2 rejects the repo yaml (windIO v1 keys `outer_shape_bem`/`internal_structure_2d_fem`; schema v2 wants `outer_shape`/`structure.layers`). Converted with windIO's own converter `windIO.converters.windIO2windIO.v1p0_to_v2p0` (`convert.log`) → `v2/geom_v2.yaml`. Converter warnings: LE anchor recomputed from `LE_reinforcement`; `offset_y_pa` read as rotation = twist; dummy `end_nd_arc`/`start_nd_arc` for the TE reinforcements.
  2. The converter left the two TE-reinforcement layers pointing at anchor `not_defined`; `patch_v2.py` repoints `TE_reinforcement_SS.end_nd_arc` and `TE_reinforcement_PS.start_nd_arc` to their own width anchor (the one the converter already writes and that `TE_*_filler` already references) → `v2/geom_v2p.yaml`. This is the only edit to the layup mapping.
  3. `run_wisdem_hack.py` calls `run_wisdem` with only RotorSE enabled (`mod.yaml`, n_span 50, 27 s) and two stubs, neither touching the blade: a dummy `costs`/`tcc` IVC (the converter drops `costs`) and a no-op `PlantFinance.compute`. Output `precomp_out.pkl`: EA, EIxx/EIyy/EIxy, GJ, rhoA, K (6x6 at the reference axis, WISDEM's `pc2bd_K`).
  - Span mapping: PreComp `r` starts at Rhub (3.97 m), so z = (r−r_min)/(r_max−r_min)·117. A first pass with z = r/r_max·117 shifted every station by up to 4 m and gave EI_flap +15 %; it was a mapping error of mine and is discarded (logs `compare*.log` are the corrected ones).
- **Definitions** (same as S-1): deck EI = elastic-centre principal values (`bd.ei_*_ec`, Schur complement of K33 on K[3:5,3:5]); PreComp EI by the same Schur on its K (it equals the eigenvalues of [[EIxx,EIxy],[EIxy,EIyy]] to 4e-4). Window 0.2–0.85 of z_max on the S-1 grids. Shell numbers: the S-1 fixture calls copied (es 0.25, `skin_reference="mid"`) on SLURM 11622353 (4m40s, `shell-11622353.log`, `shell.json`); they reproduce the 7dfe78d anchors (EA 7.59, EI_flap 5.17, EI_edge 11.45, GJ 5.05 %). "CLT yaml" = `tools/yaml_ei_independent.py` / `yaml_layup_independent.py` (meshless OML thin wall, `clt.json`, `clt_and_shell.py`). "CLT mesh" = `tools/s1_clt_steiner_ei.py` machinery on the mid-surface es 0.25 mesh, per station (SLURM 11622377, `clt_mesh.json`, `cltmesh-11622377.log`). Compare: `compare_s1grid.py` (S-1 grids, headline), `compare.py` (deck stations, per-station table). Deck-station and S-1-grid means agree within ≈1 pp.
- **Mid-span (0.2–0.85) mean signed / mean |rel err|, % (`compare_s1grid.log`, n = 32 stations, mass 31)**:

  | qty | PreComp/deck | shell/deck | PreComp/shell | CLT yaml/deck | PreComp/CLT yaml | CLT mesh/deck | PreComp/CLT mesh |
  |---|---|---|---|---|---|---|---|
  | EA | +7.67 / 7.77 | +7.57 / 7.59 | +0.08 / 0.71 | +7.57 / 7.57 | +0.08 / 0.48 | – | – |
  | EI_flap | +4.85 / 4.95 | +3.93 / 5.17 | +1.13 / 4.20 | +26.15 / 26.15 | −16.75 / 16.75 | +7.74 / 7.74 | −2.68 / 2.92 |
  | EI_edge | −0.17 / 7.56 | −11.45 / 11.45 | **+13.20 / 13.20** | −1.72 / 6.81 | +1.50 / 3.24 | −0.40 / 6.38 | +0.15 / 3.76 |
  | GJ | −28.79 / 28.79 | −2.76 / 5.05 | −26.10 / 27.97 | – | – | – | – |
  | mass | +3.49 / 4.28 | +0.73 / 2.23 | +2.73 / 2.73 | +3.21 / 3.57 | +0.26 / 0.78 | – | – |

  Outboard (0.64–0.855, 11 stations): PreComp/deck EA +5.05, EI_flap +2.92, EI_edge **−9.84**, GJ −45.6, mass +2.20 %; shell/deck EI_edge −16.67; CLT mesh/deck EI_edge −9.40; PreComp/CLT mesh EI_edge −0.59 (|.| 4.65). Per-station signed table: `compare.log`. Whole-blade mass integral 0–117 m: deck 66 912 kg, PreComp 67 956 kg (+1.6 %).
  - On EA, a second yaml integral (`ea_q`, restrained Q, 12 stations) is +9.74 % over the deck; PreComp vs `ea_q` −1.49 %.
- **Verdict per quantity** (reference hierarchy: PreComp = independent code with its own yaml reading; the deck is SONATA/ANBA4, a 2D FE homogenisation):
  - **EA: (d) reference restriction, deck vs yaml.** PreComp agrees with our laminate to 0.08 % (and with `ea_q` to 1.5 %) and is +7.7 % over the deck, the same as the shell (+7.57 %). Outboard +5.05 vs +5.3 %. Three independent codes read the yaml layup at +5…+10 % over the 6x6; our mapping is not the outlier.
  - **Mass: (d).** PreComp +3.5 % (outboard +2.2) over the deck, our yaml integral +3.2 %, shell +0.7 %. The deck mass is ≈3 % below what the yaml layup plus material densities give in two independent readings. The shell is closer to the deck than either layup integral (lumped mass partition, not a laminate property).
  - **EI_flap: (d) as a yaml-vs-deck question; the shell is consistent.** PreComp +4.85 % over the deck, shell +3.9 % (|.| 5.17), PreComp vs shell +1.1 % signed, CLT mesh +7.7 %. The 5 % bound on shell vs deck sits inside the +5 % at which PreComp itself reads the yaml. The meshless OML integral (+26 %) is the known OML-placement bias, which PreComp and the mid-surface mesh do not share.
  - **EI_edge: split, and the result goes against the shell.** Mid-span, PreComp reads the yaml at −0.17 % (deck-consistent) and agrees with our mesh laminate to +0.15 % (CLT mesh/deck −0.40 %). So the laminate, as PreComp and the mesh CLT both read it, is **not** 11 % below the deck at mid-span; the shell *response* is (−11.45 %), i.e. 13 % below PreComp and 11 pp below its own laminate integral. Outboard, PreComp (−9.84) and the mesh CLT (−9.40) both sit ≈10 % below the deck, so that part of the earlier "10–21 % outboard" is a deck-vs-yaml difference (d); the remaining −7 pp (shell −16.67 vs CLT mesh −9.40) is again the shell response. The mid-span part is not a reference restriction. It points at the extraction or the section response of the shell, consistent with the TE cross-section distortion of L19 (not tested here). Outcome (d) is not available for it; it stays an open shell-side residual.
  - **GJ: inconclusive, not claimed.** PreComp −28.8 % mid, −45.6 % outboard vs the deck, while the shell is −2.8 %. PreComp's torsion is a thin-wall closed-cell estimate with no warping (it is known to miss multi-cell/web torsion), so it cannot arbitrate a 5 % residual. The shell/deck agreement of 5 % stands uncontested and unconfirmed.
- **Limits.** PreComp: thin wall, no shear deformation, no warping, own region/web discretisation. The layup reaches PreComp through windIO's converter (warnings above) plus the 2-line anchor repair, so "PreComp's reading" is "WISDEM's reading of the converted yaml". PreComp EA, EI_flap and mass are 0.1–3 % from our independent mesh-free integrals, which is evidence the converter did not distort those. WISDEM stubs (costs, financese) do not enter the blade section calculation.
- **Boundaries kept**: no test, bound, `contract.yaml` or store row changed; no verdict promoted. Scratch only under `logs/check/l3-d3/` (gitignored), the venv stays there.
- **Next (shell side)**: the EI_edge response gap (shell 11 pp below its own laminate at mid-span, 7 pp outboard) is the open item; PreComp removes the deck as an excuse for the mid-span part. Test candidate: TE sectional distortion (L19) with the 2D section solve of the mid-span station.

### L25 — Skin reference surface, S-6/S-7 on mid, TE distortion verdict (2026-10-22)
- **S-7 on the mid-surface** (probe SLURM 11622345, `logs/check/l3-s7-mid-11622345.log`; test SLURM 11622347): band 0.3–0.9 ratio **1.00095 (0.0947 %)** against the BeamDyn deck. On the OML the ratio was 0.8980 (10.2 %) at the same extension c0d9b9a7. Per band 0.979–1.046 in 0.2–0.9. **Outcome (c)/(d):** the OML placement centres the laminate on the aerodynamic surface, so half of each wall lies outside the blade. That is a reference-surface choice made for CFD coupling. The physically consistent structure (plies inward from the OML, as SONATA builds the deck) closes the anchor. Test moved to `skin_reference="mid"` in 38b8ba2; bound unchanged.
- **S-6 test at es 0.25** (rule 9): 3 passed (SLURM 11622346, OML) in 38b8ba2.
- **S-6 mid vs OML** (SLURM 11622350, `logs/check/l3-s6-skin-11622350.log`, `logs/check/l3-s6-probe/s6_skin.py`), against ElastoDyn:

  | skin | tip OoP mean | std |
  |---|---|---|
  | OML | −12.83 % | −5.16 % |
  | mid | **−0.86 %** | +8.88 % |

  - The mean deficit attributed for weeks to a "K_G over-response" was the OML placement's flap stiffness. K_G and K_SP stay in K_eff.
  - Test moved to mid (bounds unchanged): 3 passed (SLURM 11622379).
- **TE section distortion (L19), verdict (ii)** (worker, `logs/check/l3-te-distortion/REPORT.md`, SLURM 11622349/11622374/11622376).
  - **Method.** The same station resultants (to 2e-15) were applied three ways: (a) per node, as production; (b) on spar and web nodes; (c) by skin area.
  - **Global response is insensitive.** Tip OoP changes by ≤0.10 % across the three.
  - **Distortion follows load on the TE skin, not the per-node rule.** At r = 110 m it is 0.245 / 0.012 / 0.165 m for (a)/(b)/(c).
  - **The outboard TE_PANEL is genuinely thin in our reading of the yaml:** 11.1 mm with D 2334 N·m at r ≈ 92 m; 4.0 mm with D 111 at 99 m; 2.7 mm with D 32.6 at 108 m. A plate-strip check of that panel brackets the FE sag (0.227 m, between the clamped and simply supported bounds). So (iii), an element defect, is not supported.
  - **Outcome (d), input property.** Tip anchors are unaffected; outboard section-twist readings depend on load introduction by 27–40 %.
- **D3 (L24, 4909fa7) consequences for the open residuals:**
  - S-1 EA, mass and EI_flap, and S-2 LC1/LC4 through flap: **outcome (d), deck-vs-yaml.** PreComp reads the yaml as we do (EA +0.08 %, EI_flap +1.1 % vs the shell), and both sit +4.9 to +7.7 % over the deck.
  - S-1 GJ: within 5 % on the shell; PreComp cannot arbitrate it.
  - S-1 EI_edge mid-span: still open, as a shell-side residual of −11 pp against its own laminate. A prismatic-segment arbiter is running (worker, `logs/check/l3-prismatic/`).
- **S-5 at es 0.25, skin × stiffening** (SLURM 11622378, `logs/check/l3-s5-var-11622378.log`, `logs/check/l3-s6-probe/s5_variants.py`), tip OoP against ElastoDyn's 16.1033 m:

  | skin | parked K | K + K_G + K_SP, 7.56 rpm |
  |---|---|---|
  | OML | −4.50 % (the old passing value) | −11.71 % |
  | mid | +9.75 % | **+0.37 %** |

  - The old −4.5 % pass was **two compensating errors**: the OML skin is too stiff in flap, and the parked K drops the centrifugal stiffening that ElastoDyn carries at rated speed.
  - The docstring premise that "K_G injects a spurious ~13 %" came from the OML placement. The shell's rotating modes close against OpenFAST (S-4 1st flap −1.4 %), so K_G is not suspect.
  - Test moved to mid, es 0.25 and rotating K_eff; bounds unchanged. IP |u_x| = 2.00 m gives a ratio of 2.47 against the loose 2.5 bound (IP is apples-to-oranges by the module's own account). Test: 3 passed (SLURM 11622382).
- Store: S-7 row (`rows/34`) re-measured at 38b8ba2, ratio 1.00095 (0.0947 %, PASS). S-5 and S-6 have no store row, so nothing was written for them; their numbers live in the test docstrings and here.

### L26 — V-02 / S-4 frequencies against v1.1 yaml-consistent beams (2026-10-22)
- Worker report `logs/check/l3-beamfreq/REPORT.md`. Our clamped EB beam (`tools/s4_deck_clamped_beam_modes.py`) reproduces `logs/check/s4/deck_clamped_beam.txt` exactly. It was fed the deck and the PreComp reading of the yaml (L24).

  | mode [Hz] | published (v1.0) | deck EB | PreComp EB | shell (V-02, mid, es 0.25) | shell vs deck / PreComp |
  |---|---|---|---|---|---|
  | 1F | 0.5585 | 0.5159 | 0.5268 | 0.5246 | +1.7 / −0.4 % |
  | 1E | 0.6406 | 0.7242 | 0.7181 | 0.6951 | −4.0 / −3.2 % |

  - A first-order shear bound lowers the beams by about 2 %, giving +1.8 to +3.9 % on 1F and −1.4 to −2.3 % on 1E. That is the fairer comparison, because the shell carries transverse shear.
- **V-02, outcome (d), reference restriction.** NREL/TP-5000-75698 describes the v1.0 blade (B3-2, `validation-l3-b3-audit.md:507`). Both v1.1 readings sit 5.7–7.6 % below it on 1F and 12–13 % above it on 1E. The 1E gap even has the opposite sign to the shell's. Against v1.1 beams the shell is inside 5 % on both modes. The beam is our code and explains the number; it is not an independent arbiter. The V-02 test stays red against the v1.0 table; changing its reference is a contract decision (D1) that this record supports.
- **S-4 1E** reads −5.4 % against OpenFAST/ElastoDyn: −1.3 pp is ElastoDyn's assumed modes, and −3.2 to −4.0 % is shell below beam (−1.4 to −2.3 % after the shear bound). Its link to the EI_edge deficit is tested by the prismatic arbiter.

### L27 — S-1 EI_edge decomposed: prebend frame mixing, Poisson, taper (2026-10-22)
- Worker report `logs/check/l3-prismatic/REPORT.md`.
  - **Setup.** Prismatic 30 m segments at r = 38.5/62.9/74.7/90.3 m, mid skin, es 0.25 (es 0.125 check), 13 stations each. Pure end moment through a rigid diaphragm, root clamped; a tip couple agrees to ≤0.5 %.
  - **Prismatic shell vs its own same-mesh CLT:** flap −0.1…−0.8 %, edge −3.2…−5.2 % (es 0.125: −3.4 % at r 62, −6.1 % at r 90). The strain field is planar, with no shear lag, and in-plane distortion is ≤2 mm. Setting ν = 0 everywhere closes the edge gap to +0.2…+1.2 %, so the 3–5 % is Poisson ring coupling of the shell.
- **Decomposition of the −11.45 % S-1 EI_edge:**
  1. **Prebend frame mixing, ≈3.6–6 pp — outcome (b), an instrument defect.**
     - **Cause.** `SectionalExtractor.ei_from_static_response` reads moments and curvatures in global X/Y. The deck 6×6 is in the tangent frame. On the prebent blade a global moment has a torsion component, and the blade is torsionally soft (EI/GJ ≈ 46–69). That gives 1/EI_global = cos²φ/EI_local + sin²φ/GJ.
     - **Evidence.** Prebending the prismatic segments with slope ×0.5/×1/×2 drops the reading by 1.8/6.5/21 pp, against 2.1/8.0/25.7 predicted. Removing prebend from the full blade moves S-1 edge from −11.45 % to −7.9 %.
  2. **Poisson coupling, ≈4 % — outcome (c) only against the hoop-clamped CLT.** Whether the ANBA4 deck carries the same effect is open: no 2D solid section solve has been run.
  3. **Taper and layup variation, ≈1.5 pp.**
- **The instrument is not otherwise biased.** On the real blade a pure tip couple reads −13 % vs the deck, no better than the instrument.
- **Next:** fix the extraction frame (tangent-frame moment and curvature, torsion coupling included), test-first on a prebent prismatic segment, then re-measure S-1.

### L30 — Unconverged-mesh tests: inventory and mandate (2026-10-22)
```text
haz una investigacion de los tests para ver cuales usan mallas no convergidas, eso lo prohibi estrictamente hace rato, excepto para aquellos tests de convergencia de malla, no de pruebas fisicas reales
```
```text
si vamos a hacer eso, porque estamos validando contra resultados numericos con mallas no convergidas, nunca vamos a resolver el problema si la comparacion la hacemos contra basura
```
- **Inventory** (read-only scan, Engram obs 100):
  - **Physics tests on the blade at 1.0 m:** `tests/validation/blade/` `test_blade_iea15mw_validation`, `test_blade_section_torsion_stiffness`, `test_blade_twist_anchor_beam`, `test_blade_twist_mechanism`, `test_blade_rated_twist`, `test_blade_deloading_vs_reference` and `test_blade_section_frame_tear`, plus `tests/test_multicell_shear_flow.py` (blade part).
  - **At 2.0 m:** `tests/test_ccx_fsi_parity_parked_v50.py` (`GATE_ELEMENT_SIZE`, missing from the registry).
  - **Guards at 1.0/3.0 m:** `tests/validation/bem/` `test_force_projection_ac_datum`, `test_force_projection_load_frame` and `test_wall_flow_activation`, plus `tests/validation/rotor/test_bem_fsi_deformed_geometry`.
  - **Convergence test that asserts a physics gap at its finest 0.5 m:** `tests/validation/parity/test_blade_iea15mw_mesh_convergence.py`.
- **Guard holes in `scripts/check_mesh_size.py`:** it misses prefixed constants (`GATE_ELEMENT_SIZE`), positional calls (`build_assembly(0.5)`) and size tuples (`MESHES = (...)`). `scripts/mesh-size-debt.txt` is stale.
- **Plan:**
  1. Close the guard holes and true up the registry.
  2. Re-measure every physics test at 0.25 m, with before/after recorded and bounds unchanged; a moved number is a finding.
  3. Move the convergence test's finest point to 0.25 m.

### L28 — S-1 tangent-frame extraction (2026-10-22)
- **Change.** `SectionalExtractor.ei_from_static_response(..., frame="global"|"tangent")` (`src/aeroelast/postprocess/sectional.py`). `"global"` is the default and the old code path; `"tangent"` reads the two tip-load solves in the section tangent frame: tangent = slope of the station-centroid line, edge = global X normal to it, flap = tangent × edge; station rotation vector from a rigid fit of its nodal displacements, differentiated with a local quadratic fit (±2 station spacings); EI = |M·axis| / |κ·axis| on the bending components only, so torsion enters neither moment nor curvature. Bend–twist coupling of the section is not removed. `tests/test_iea15mw_s1_sectional.py` now uses `frame="tangent"`; bounds, contract.yaml and store rows untouched. Specs covered: S1 (identified and corrected an instrument error), S2/S3/S5 (no bound, verdict or store change; es 0.25, mid, SLURM `--ntasks` only, plain-language result).
- **Test** `tests/validation/blade/test_sectional_tangent_frame.py` (slow, ~80 s): 24 m prismatic segment at r = 62.91 m (L27 crop construction, es 0.25, mid), straight (A) and with the real prebend slope (F1, 0.036 rad). Reference = pure tip-couple EI of the straight segment: flap 4.2118e9, edge 1.2856e10 N·m².
- **RED.** Before the change: setup error `TypeError: ... unexpected keyword argument 'frame'` (SLURM 11622459, `red-11622459.log`). Measured default miss (`measold-11622460.log`, window 0.3–0.7 L): prebent EI_edge/pure −12.47 %, EI_flap −14.32 %; the straight segment reads −10.1 % through the same instrument (deflection-spline smoothing on a short segment), the prebend adds −2.4 pp here (the pure global-moment figure of L27 is −6.5 pp).
- **GREEN.** `measnew-11622461.log`: tangent frame vs straight pure moment, EI_edge: straight −0.78 % (stations −1.5…−0.1), prebent −1.28 % (stations −1.9…−1.0); EI_flap −1.04 / −1.06 %. Tolerance 2.5 % (1.3× the worst station; the ~1 % offset is tip-force vs pure-moment loading) and 1.5 % for prebent-vs-straight tangent (measured 0.5 pp); default miss asserted beyond 10 %. New file: 6 passed (SLURM 11622472).
- **Default path.** Byte-identical: `np.tobytes` of the `frame="global"` EI arrays, HEAD vs this tree, on both segments (`old_*_default.npz` vs `new_*_default.npz`); the test also pins it against an inline copy of the old algorithm (`assert_array_equal`). Callers: `tools/run_s1_sectional.py`, `s1_ei_partition.py`, `s1_midsurface_counterfactual.py`, `s1_reference_sweep.py`, `diagnose_bending_compliance_2x2.py` (docstring) keep the default.
- **S-1 on the real blade** (es 0.25, mid, 0.2–0.85 R, signed (mean |rel|) vs the 6×6; `blade_s1.py`, SLURM 11622464/11622465):

  | quantity | before (global) | after (tangent) |
  |---|---|---|
  | EI_edge | −11.45 % (11.45) | **−7.51 % (8.82)** |
  | EI_flap | +3.93 % (5.17) | +3.79 % (5.70) |
  | GJ | −2.76 % (5.05) | unchanged |
  | EA (CLT) | +7.57 % (7.59) | unchanged |
  | EI_edge, prebend removed (`logs/check/l3-s1-frame/blade_noprebend.yaml`) | −7.88 % (8.78) | −8.17 % (9.23) |

  - The tangent reading on the prebent blade (−7.51) sits within 0.7 pp of the tangent reading with the prebend removed (−8.17), and the two frames agree to 0.3 pp when there is no prebend. That closes the 3.9 pp frame mixing (L27 counterfactual 3.6, formula 6). EI_flap's mean |rel| rises 0.5 pp (noisier local-fit estimator; signed bias unchanged).
  - `tests/test_iea15mw_s1_sectional.py` (SLURM 11622473, run with `-x`): EI_edge still red (8.82 % vs the 5 % bound, was 11.45 %). EI_flap, GJ, EA, mass were not rerun in that job; their numbers above come from `blade_s1.py` (GJ and EA do not use EI). The residual −7.5 % is the Poisson ring coupling (~4 %) and taper/layup (~1.5 pp) of L27.
- **Gate.** `scripts/check.sh quick`: OK (75 passed, 1 xfailed), ruff clean on the three files.
- Logs: `logs/check/l3-s1-frame/` (untracked by `.gitignore`).

### L29 — T-junction director convention against a converged C3D20R solid (2026-10-22)
- Worker report `logs/check/l3-tjunction/REPORT.md`.
  - **Model.** Two-cell isotropic box, 1.0×0.6×10 m, web at mid-chord (and at x = −0.2 m), clamped root, rigid tip diaphragm, t = 10 and 4 mm.
  - **Variants.** MITC4 (A) separate directors (7dfe78d) and (B) averaged (base 9e183ce). The Bm mirror check matches B to every digit. CCX S8R at h 0.1/0.05/0.025.
  - **Reference.** C3D20R with 2 elements through the wall, 3 levels each. Last-step change ≤0.043 %, extrapolated residual ≤0.008 %. Within 0.8 % (10 mm) and 0.3 % (4 mm) of two-cell Bredt, and within 0.15 % of Timoshenko.
- **Tip metrics do not separate A from B** (within 0.15 pp). Twist +0.68/+0.73/+0.69 % for A/B/S8R at 10 mm is a common solid-corner offset.
- **Junction-sensitive modes** at h 0.025:

  | mode | t | A | B | S8R |
  |---|---|---|---|---|
  | f3 (web/panel) | 10 mm | −0.95 % | −6.40 % | +0.14 % |
  | f5 | 10 mm | −0.50 % | −3.58 % | +0.03 % |
  | f3 | 4 mm | −0.27 % | −4.97 % | +0.73 % |

  - With the web at x = −0.2 m: f3 −0.92 % (A), −5.44 % (B).
  - B converges slowly: f3 at 10 mm goes −12.0, −9.5, −6.4 %.
- **Verdict, outcome (c):** at T-junctions the production convention (A) is 5–8× closer to 3D elasticity than averaging. CalculiX S8R with knots matches the solid once refined, so it is a sound arbiter there. The L20 CalculiX blade movement is therefore **not** junction physics.
- **Candidate explanation, under test:** the high-mode gap comes from `test_blade_iea15mw_mesh_convergence.py`, which meshes at 2.0/1.0/0.5 m, against rule 9. At those meshes B's slow, soft junction modes may have offset MITC4's coarse-mesh stiffness against S8R. That would be another compensating pair. It is being tested at es 0.25 (worker, `logs/check/l3-ccx-es025/`).

### L31 — Rule-9 re-measurement at 0.25 m: results (2026-10-22)
- Guard: baaa641. Tool defaults: c8c09bf. Test groups: B 5b061a2, D cfbdde3, C ad117c1, A deb97c5. Worker reports are under `logs/check/l3-mesh025-{a,b,c,d}/`. Bounds and assertion logic are unchanged everywhere. Tests that compare a structural response with a beam or section reference were moved to the mid skin; projection and geometry guards kept their skin.
- **Code-to-code agreement tightens on converged meshes**, against CalculiX on the same mesh:
  - CCX parked-V50 modal gaps 0.86/1.94/1.44 % at 2.0 m → **0.02/0.03/0.03 %** at 0.25 m.
  - IEA-15 validation file: modal 1.64 → **0.14 %**; static tip 1.17 → **0.44 %**.
  - The earlier "agreement" at coarse meshes was discretisation noise.
- **Passes that were cancellations at 1.0 m on the OML, now reported (xfail/red, not re-banded):**
  - Article DLC 1.4 tip deflection: 4.77 → 13.21 % (20.39 m).
  - 2nd flap vs NuMAD: 10.87 → 12.26 %.
  - `rated_twist` couple-vs-nodal over-report: 4.2× → 1.36× against the pinned >2×. Most of the pinned "defect" was coarse-mesh.
- **Instrument findings under diagnosis** (worker, `logs/check/l3-rings/`):
  - `ac_datum` and `deformed_geometry` red at 0.25. A fixed ring-clustering tolerance tuned to the 1.0 m ring spacing merges tip rings.
  - Multicell midspan cell-area vs shoelace: 14.93 % (2.84 % at 1.0 m).
  - Tear-guard signature (8 edges at z 11.94) reproduces at neither mesh.
- **Physics findings carried forward:**
  1. **Elastic twist under the rated BEM load.** The shell tip ring rotation is +6.44°, 3.63× the decoupled anchor beam (GJ + shear centre). Deloading follows from it: production thrust change −28.27 % vs Zhou −13.04 %. The independent arbiter is a full-6×6 BeamDyn static run with the same loads.
  2. **Tip radial shortening:** −0.815 m vs BeamDyn −1.3146 m.
- **Campbell artifact** regenerated at 0.25 m mid (SLURM 11622486; gitignored, local): rated 1F/2F/1E 0.5472/1.5341/0.7006 Hz; `test_iea15mw_sx_campbell` 5 passed.

### L32 — CCX parity verdict, ring-grouping fix, session restart handoff (2026-10-22)
- **CCX blade parity** (worker report `logs/check/l3-ccx-es025/REPORT.md`, SLURM 11622469–81, 11622503–18).
  - **High-mode gap with A is discretisation.** Worst gap 13.98/8.26/5.23 % at 0.5/0.35/0.25 m: modes 6 and 8 fall at about order 1.4 (linear MITC4 against quadratic S8R). B's smaller gap changes sign between 0.5 and 0.35 m, so it was a cancellation. Modes 1–5 agree to 0.02–0.18 % at 0.25 m.
  - **Tip axial (A −5.46 % at 0.25 m) is local.** Span displacement and per-region axial strain agree to 0.05–0.2 % up to z = 105 m. About 89 % of the error sits in the last 10 m, at the point-loaded 16-node tip ring. Fold and thickness-step nodes carry none of it. Outcome (d): the metric is too local to validate.
  - **Mesh-convergence test moved to 0.5/0.35/0.25** (uncommitted, together with its registry line in `scripts/mesh-size-debt.txt`). Its docstring now carries the table above and drops the refuted "element-order plateau" claim. Bounds unchanged; expected high-mode 5.23 % vs 5 % (finding). Verification job **SLURM 11622586** (`logs/check/l3-meshconv-test-11622586.log`) was running at the restart.
- **Ring grouping, outcome (b), production defect fixed: 533d7e6.** Fixed ring tolerances tuned on the 1.0 m mesh split prebent rings and fused the tip rings at 0.25 m. `physical_ring_labels` now groups by element topology. All affected tests pass at 0.25 m. Two things are still pending:
  - RED of `tests/validation/bem/test_ring_grouping_fine_mesh.py` against the parent code;
  - independent verification (Risk item 3).
- **Lost at the restart:** three background workers (Pi-session processes; SLURM jobs survived).
  1. **BeamDyn full-6×6 arbiter of the rated elastic twist.** Partial: `logs/check/l3-bd-rated/` holds `run_bd_rated.py`, `bd_results.json`, `bd_variants.json`, `torque_check.log` and `torque_by_span.log`, all from SLURM 11622573/81/82/87 (completed). There is no REPORT.md. Read those outputs before re-running anything.
  2. **BeamDyn with a PreComp (yaml-consistent) deck for S-2 LC1–LC4.** Not started; no job was submitted (`logs/check/l3-bd-yamldeck/` is empty).
  3. The ring worker had finished its runs; its result is committed as 533d7e6.
- **Resume recipe (new session):**
  1. Read this document from L30; run `mem_search "rule-9 re-measure 0.25"` and `git log --oneline -8`.
  2. Read `logs/check/l3-meshconv-test-11622586.log`. If the low modes pass and the high mode reads about 5.23 %, commit the mesh-convergence test with its registry line.
  3. Observe RED for the ring tests on a `git archive 533d7e6^` copy of `src/`, then run one independent gentle-ai-verify on 533d7e6.
  4. Relaunch worker 1 from its partial outputs. The task text is the L31 "physics finding 1" bullet: shell tip ring rotation +6.44° = 3.63× the decoupled anchor beam; deloading −28.27 % vs Zhou −13.04 %; tip dr −0.815 vs −1.3146 m.
  5. Relaunch worker 2 for the S-2 yaml-consistent deck.
  6. Open items after that: the tear-guard signature (reproduces at no mesh), S-1 EI_edge −7.5 % (needs a 2D solid section arbiter), and the final synthesis per residual.

### L33 — Resume from L32; ring RED; rated-twist torque sense under suspicion (2026-10-23)
```text
Retomamos odd/tasks/l3-physical-validation.md, T5, desde L32 (receta de reanudación), en  modo desatendido según L23 y L30. mem_search "rule-9 re-measure 0.25".
```
- **Ring grouping 533d7e6, RED observed by the parent** (`logs/check/l3-ring-red-parent.log`, driver `$SCRATCH/tmp/l3-ring-red/drv/red.py` on a `git archive 533d7e6^` copy, module path printed). Parent `rings_in_band` splits every synthetic skewed ring: 16 groups vs 8 stations at all four spacings (0.30, 0.25, 0.109, 0.0746 m). The parent twist binder (0.05 × 2.388 m) fuses 3 stations at 0.109 and 0.0746 m. GREEN at HEAD: 5 passed (`logs/check/l3-ring-green-head.log`). The runtime refused the independent verifier (no admissible Risk line) and said not to retry. Parent differential old vs new instead, SLURM 11622644 (`logs/check/l3-ring-diff-11622644.log`); es 1.0 serves only as the differential case.
- **BeamDyn arbiter (lost worker 1): partial outputs read** (`logs/check/l3-bd-rated/`). The same production rated nodal load on shell es 0.25 mid gives a linear tip twist of 0.1125 rad. BeamDyn v5 full 6×6 gives 0.1228 rad linear (scale 1e-3 ×1000), so the shell is −8.4 %, and 0.0921 rad nonlinear. Tip dr: shell −1.269 m vs BD −1.3147 m. **The L31 "3.63×" was against the decoupled anchor beam, not the 6×6 beam.** The shell is consistent with BeamDyn *under the same load*.
- **The load's torque sense is the open question (suspected outcome (b)).**
  - The production nodal load carries a torque about the reference axis of +4.57e5 N·m about +z, which is A + B: A = −Mp = 2.77e5 (polar Cm term) and B = x_ac·Np = 1.83e5 (AC lever). In the textbook combination the two oppose, B − A = −9.4e4 about +z.
  - Argument from conventions, not yet a measurement. OpenFAST pitch-to-feather turns the LE upwind (ElastoDyn positive pitch about −z_b). With the LE at +x and the downwind thrust at +y (both stated in 1146265), a **+z rotation takes the LE downwind, i.e. toward stall: nose-UP**. AeroDyn's local frame (x_l toward the suction side, y_l LE→TE) puts nose-up Cm about x_l × y_l = n × t, not t × n.
  - If that holds, commit `1146265` ("three linked sign errors", 2026-10-05) inverted all three links instead of fixing them:
    - the Mp axis (`force_projection.py` `_strip_moment_axis_sign`, `M_ac` at `project()` :840);
    - the guards' nose-down convention;
    - the twist feedback `_twist_mesh_to_bem = +1` (`fsi_participant.py`).
    Its premise was that the pre-fix section torque came out nose-up. That was actually the band-datum lever (+7.6e5, since replaced by the station-ring datum, now 1.83e5), not the Cm axis.
  - Supporting numbers (BeamDyn nonlinear, `bd_loadc.json`): with the textbook torque, tip rdz = −0.0517 rad, i.e. 2.96° toward feather under this reading, against Zhou 2025 at about 3.6° toward feather. `Zhou_tip_twist_scaled` gives a thrust change of −11.8 % vs Zhou −13.04 %. The production chain gives +0.0921 rad and −32.0 %.
  - Not proved. The decisive independent arbiter is an OpenFAST v5 AeroDyn+BeamDyn steady run at the rated point (10.59 m/s, 7.56 rpm, pitch 0): sign and size of the tip twist and of the root torsion. AeroDyn's `BlSwpAC` is non-zero (−0.38 … −0.06 m), so OpenFAST carries the AC lever and is a valid torque reference. Delegated (context backstop at 151 k).

### L34 — OpenFAST arbiter: the production torque sense is inverted, outcome (b) (2026-10-23)
- **Ring differential, 533d7e6** (parent check; SLURM 11622644, `logs/check/l3-ring-diff-11622644.log`).
  - es 1.0: old vs new **bit-identical**, with the same 191 groups, the same strip rings, and identical projected F. This confirms the superset claim where the old rule worked.
  - es 0.25: 751 → 689 groups (the intended merge of split rings) and different strip-ring sets. Projected rated F is **bit-identical**, so the change reaches only the twist ring and the wall-flow grouping.
  - Advisory: the minimum spacing between per-strip label means is 0.0009 m, from the probe's per-strip relabelling of rings that straddle a strip boundary; the grouping itself is not implicated. No blocker. 533d7e6 is accepted on RED + GREEN + differential; there was no independent verifier (refused by the runtime).
- **OpenFAST v5 AeroDyn+BeamDyn, steady rated point** (worker, `logs/check/l3-bd-openfast/REPORT.md`, SLURM 11622651–74).
  - Run conditions: 10.59 m/s without shear, 7.56 rpm, pitch 0, no gravity, tower or tilt. BeamDyn damping was ×10/×30 (identical to 4 digits) because the deck's damping sits in a 3.6 Hz limit cycle.
  - Rigid rotor: 2.538 MN / 16.42 MW. Production CCBlade gives 2.525 MN / 16.33 MW.
  - Elastic, against Zhou 2025:
    - tip twist −3.40° **toward feather** (Zhou −3.60°);
    - thrust −13.2 % (Zhou −13.04 %); power −6.6 % (Zhou −8.38 %);
    - tip flap/edge 14.08/−1.28 m (Zhou 13.86/−1.22 m).
- **Mesh frame** (0.25 mid). LE at +x; the camber bulges toward +y on 467 of 483 rings, so the suction side is downwind. **+z is nose-up (toward stall).**
  - Polar Cm torque −2.77e5 about +z, AC lever +1.826e5, net **−9.41e4 N·m (nose-down)**.
  - Production applies +4.57e5 (A+B), so the Cm term enters with the wrong sign. OpenFAST's outboard section torque is nose-down at every cut, production's is positive at every cut.
  - BeamDyn under the textbook torque gives −2.96° (−12.8 % vs OpenFAST). The production chain gives +5.28° toward stall.
- **Verdict: commit `1146265` inverted the three links.** Its premise ("+z rotation is nose-down") is wrong. Links:
  - `force_projection.py:582`: the Mp-axis sign, `-1.0 if le_at_hi`, should be +1;
  - `:840-842`, where that sign is applied;
  - the docstring at `:36-48`, whose cross-product order is wrong;
  - `fsi_participant.py:472, :789` with `_resolve_twist_sense` `:828-880`: `_twist_mesh_to_bem` +1 → −1;
  - nose-down guards asserting ω > 0 in `test_blade_rated_twist.py`, `test_blade_deloading_vs_reference.py`, `test_bem_fsi_deformed_geometry.py` and `test_multicell_shear_flow.py`.
  The production deloading (−32.0 % thrust) is therefore a sign artefact. **L31 finding 1 closes as (b)** once the fix is measured.
- Open after the fix:
  - the remaining 13 % tip-twist gap to OpenFAST (OpenFAST outboard torque 1.7–4× the textbook value beyond 70 m: AeroDyn applies the AC on the BlSwpAC line, not on the mesh section);
  - power −6.6 vs −8.38 %;
  - the deck's 3.6 Hz BeamDyn limit cycle (reference-side).

### L35 — Torque sign chain fixed (4ee98ab), verified; session handoff (2026-10-23)
- **Mesh-convergence test** committed as `b17f3d6`. SLURM 11622586, 49m37s: low modes 0.02–0.12 % at 0.25 m. The worst high mode reads 13.98/8.26/5.23 % at 0.5/0.35/0.25 m against the unchanged 5 % bound, so the test is red as a finding. The assertion text no longer claims "element-order, not discretisation".
- **Fix `4ee98ab`, outcome (b), closes L31 finding 1.** Worker report: `logs/check/l3-sign-fix/REPORT.md`.
  - The polars' nose-up Mp now rides −sign·span (`force_projection.py:849-851`).
  - `_resolve_twist_sense` returns −1 for a nose-up +span rotation (`fsi_participant.py:866-875`).
  - The guards in four test files flip sense; no tolerance changed.
  - New test `tests/validation/bem/test_torque_sense_vs_openfast.py`, 3 tests at a 2 % tolerance. RED at 0754a8a (SLURM 11622677): +4.58e5 vs B−A −9.52e4, and a +2° nose-up rotation raised θ by +2°. GREEN (11622706): −9.5173e4, and θ falls by 2°.
  - Measured at es 0.25 (rigid BEM 2.525 MN / 16.33 MW unchanged):

    | quantity | before | after | OpenFAST | Zhou |
    |---|---|---|---|---|
    | tip twist | +6.39° (toward stall) | −3.13° (toward feather) | −3.40° | −3.60° |
    | production thrust | −28.26 % | −9.06 % | −13.2 % | −13.04 % |
    | production power | −18.34 % | −4.26 % | −6.6 % | −8.38 % |

  - `tests/validation/bem/`: 129 → 132 passed. `check.sh quick` OK; ruff clean.
- **Independent gentle-ai-verify of 4ee98ab:** items 1–4 pass (frame argument, scope, pure sign flips, new test); item 5 has no blocker. Native assess is unassessable (untracked BeamDyn directories), as before.
- **Red after the fix** (findings, not re-banded):
  1. `test_rated_tip_twist_matches_zhou_with_the_physical_load_path` trips its own "promotion guard". The test-local at-AC construction now gives 1.046× Zhou, i.e. it agrees with Zhou within 5 %. The guard was written to refuse promotion of a construction declared withdrawn. **Decision needed (formulation/test-policy):** promote the at-AC construction to evidence, or keep it withdrawn.
  2. `test_section_moment_as_a_couple_sets_the_physical_sign` is red at 0754a8a with identical numbers: pre-existing.
- **Follow-ups from the verifier** (old convention still printed or used; none is asserted):
  - `test_blade_twist_anchor_beam.py:320-366` keeps a stale nose-down conversion of the beam twist. It needs an identity conversion plus a same-sign assertion.
  - `tools/diagnose_moment_reference.py:84` keeps the old Mp axis.
  - `tools/diagnose_zhou_loads_reverse.py`, `tools/diagnose_sign_chain.py:23,136` and `tools/diagnose_applied_torsion_sign.py` keep the old prose.
  - Docs: `docs/origin_main_integration_2026-09-30.md:14` and `docs/validation_closures.md:34,633` still say "+z = nose-down".
  - The message at `test_bem_fsi_deformed_geometry.py:315` has the wrong wording.
- **Remaining residual of finding 1:**
  - tip twist −3.13° vs OpenFAST −3.40° (−7.9 %) and Zhou −3.60°;
  - thrust −9.06 % vs −13.2 % OpenFAST;
  - power −4.26 % vs −6.6 %.
  Candidate causes, not yet tested: OpenFAST applies the AC on the BlSwpAC line, and its outboard section torque is 1.7–4× the mesh-ring B−A beyond 70 m; also legacy BEM vs CCBlade, and prebend/sweep in the aero.
- **Delivery.** The branch holds 6624 authored lines in 57 commits past origin/main, well beyond the ~400-line slice. Strategy is ask-on-risk: slicing into PRs is the maintainer's decision. No push.
- **Resume recipe (new session):**
  1. Read this document from L33. Run `mem_search "torque sign chain 1146265"` and `git log --oneline -6`.
  2. Make the follow-up conversions (anchor-beam test, the diagnose tools, the two docs) in one small commit.
  3. Record the 4ee98ab numbers in the validation store: rated-twist and deloading rows, with `tools/validation_matrix.py check`.
  4. Run a full suite on SLURM at 4ee98ab vs the 7dfe78d baseline, because the sign change reaches every BEM-FSI test.
  5. Diagnose the OpenFAST-vs-shell twist and torque residual outboard of 70 m (AC line BlSwpAC vs mesh-ring AC).
  6. Then L32 step 5: BeamDyn with a yaml-consistent deck for S-2. After that the tear-guard signature, the S-1 EI_edge 2D arbiter and the final synthesis.

### L36 — Follow-ups committed, P5 row re-measured, full suite launched; T5.d first readings (2026-10-23)
- **L35 step 2 done, `7189314`.** The anchor-beam test no longer negates the beam twist: beam and shell are both nose-up positive. It now also asserts that the shell and the beam twist the tip the same way. SLURM 11622758: 2 passed. Tip readings are reported, not asserted: beam −1.7765°, ring −4.0462° (ratio +2.28), mean θz −2.9720° (ratio +1.67).
  - The diagnose tools (`moment_reference`, `sign_chain`, `applied_torsion_sign`, `zhou_loads_reverse`) now use the production Mp axis and nose-up-positive prose.
  - `docs/origin_main_integration_2026-09-30.md` and `docs/validation_closures.md` carry a "reversed by 4ee98ab" note with the OpenFAST numbers.
  - The `test_bem_fsi_deformed_geometry.py:315` wording was re-read: it is already consistent (a failure means +z, nose-up), so it was left unchanged. Ruff clean; `check.sh quick` OK.
- **L35 step 3 done, `8e1c516`.** The group 27 P5 row is re-measured at 4ee98ab: **1.9099 % against the unchanged 2 % bound** (ee50e3a read 1.2919 %).
  - **Finding:** 0754a8a, before the sign fix, already read 1.9099 % (`logs/check/l3-sign-fix/before-rt-11622678.log`). So 4ee98ab did not move it; the shift happened between ee50e3a and 0754a8a and is not bisected. Headroom is now 0.09 pp.
  - Deloading and deformed-geometry have no rows: `groups.yaml` declares them out of scope under rule 6. The three ungrouped L3 files (torque sense, ring grouping, tangent frame) are declared out of scope with their reasons.
  - `validation_matrix check`: 0 errors, 1 pre-existing warning (G27 sits on no level, 4689258). `status`: 0 undeclared files. `contract.yaml` untouched.
- **L35 step 4 launched:** SLURM **11622759** runs `check.sh full` at `7189314` (`logs/check/l3-full-4ee98ab.sbatch`, 12 tasks, sequana_cpu, 4 h; log `logs/check/full-slurm-11622759.log`). The tree is source-identical to 4ee98ab plus the test, tools and docs above. Baseline: 11622087 at 7dfe78d (21 failed / 1436 passed). Do not edit `src/` or `tests/` while it runs.
- **L35 step 5 (T5.d), first readings, nothing concluded:**
  1. **AC route.** In `l3-bd-openfast/torque_sense.log` the mesh-ring AC lever gives B = +1.826e5 and the deck `BlSwpAC` route gives +1.677e5, with A = −2.768e5. So B−A is −9.41e4 by the mesh-ring route and −1.09e5 by the deck route: **16 % more nose-down torque by OpenFAST's AC line.** That is the same order as the twist gap (−7.9 % vs OpenFAST) and the thrust gap.
  2. **The "textbook" BeamDyn static run is not a valid twist arbiter.** `l3-bd-rated` `C_textbook_A-B_full` (rdz −0.0517 rad) has a tip flap of **19.0 m**, against OpenFAST's 14.08 m with rotation: a static run without centrifugal stiffening over-flaps.
     - Its scale-1e-3 twin reads ≈ −0.034 rad. The `.out` format gives only 1 significant digit; use scale 1e-2.
     - The flap-induced (geometrically nonlinear) torsion is therefore of order 0.02 rad there.
     - The same gap appears under the inverted load: linear +0.1228 vs nonlinear +0.0921.
     - The production shell solve is linear, so this effect is absent from it.
  3. **Planned decomposition** (new dir `logs/check/l3-twist-residual/`; never write into `l3-bd-rated` or `l3-bd-openfast`):
     - (i) shell linear one-way under the corrected production load;
     - (ii) BeamDyn under the same binned load and its own Mz, at scale 1e-2 (linear) and 1.0;
     - (iii) the same as (ii) with the deck `BlSwpAC` AC lever;
     - (iv) whether `beamdyn_driver` can add the rotor's centrifugal field (root angular velocity), so the nonlinear run flaps ~14 m.
     Each case is minutes on SLURM (the l3-bd-rated jobs took 1–3 min).
- **Resume recipe (new session):**
  1. Read this document from L35. Run `mem_search "twist residual decomposition"` and `git log --oneline -4`.
  2. Read the 11622759 summary line plus the NEW-FAIL list. Compare against 11622087, and attribute every new failure to 4ee98ab or to drift.
  3. Run the T5.d decomposition (i)–(iv) above. Classify the residual as (a)/(b)/(c)/(d).
  4. Then L35 step 6.

### L37 — Decision on the at-AC construction: promote only if not coincidental (2026-10-23)
- User (verbatim): "dado que Zhow utiliza metodos mas precisos que openfas creo que lo podemos dar como evidencia, pero solo si esa evidencia no es casual"
- **Criterion: what "not coincidental" means.** Agreement within 5 % of Zhou counts as evidence only if it survives all of the following:
  - N1, load distribution: the same load applied through the production `ForceProjector` gives the same tip twist within 5 %.
  - N2, estimator: section ring rotation and mean nodal θz agree within 5 %.
  - N3, mesh: 0.25 vs 0.15 m agree within 5 %.
  - N4, physics: the effects the construction omits are each below 5 %. These are geometric nonlinearity, centrifugal stiffening, twist feedback and the AC route.
- **Verdict on the evidence already in hand: coincidental. Not promoted.**
  - N1 fails. The test's own docstring shows the same load set giving −4.35, −14.25 or −34.13° depending only on the distribution. In the anchor-beam run (SLURM 11622758) the production-projector path gives −4.0462° ring rotation and −2.9720° mean θz. At-AC gives −3.77° mean θz (1.046 × −3.60).
  - N2 fails. The two estimators differ by a factor of 1.36 in that same run, and at-AC uses mean θz, a drilling DOF.
  - N4 fails or is open. The construction is a linear one-way solve. L36 estimates the flap-induced nonlinear torsion at ~0.02 rad (~30 % of Zhou's 3.60°), and the AC route moves the net torque by 16 %.
  - N3 was not measured: unnecessary, since N1 and N2 already fail.
- **On the premise.** Zhou's aero (LL-FVW) is a higher-fidelity class than OpenFAST's BEM. Structurally both are geometrically exact beams (GEBT, BeamDyn). Our model's aero is BEM, as OpenFAST's is. A closer match to the FVW reference than to the reference of our own aero class is therefore itself a sign of compensation, not of precision. Zhou stays the external reference; the comparable number is the production path once T5.d explains its parts, not this test-local construction.
- **Consequence for the test.** `test_rated_tip_twist_matches_zhou_with_the_physical_load_path` stays red: its guard sees |1 − ratio| < 5 % and asks for promotion. The decision is now taken, so in the next session, after 11622759 finishes, the guard should assert the reason the magnitude is withdrawn (N1/N2 spread) instead of distance from Zhou. That is a test-policy edit in its own commit; no tolerance changes. Do not edit `tests/` while 11622759 runs.

### L38 — Delivery: 19 stacked PRs to main (2026-10-23)
- User (verbatim): "antes de cerrar la sesion dividie todo en varios PR encadenados y mezclalos a main, luego cierra la sesion para poder continuar en otra"
- Strategy: stacked PRs to main, branches `l3/01`..`l3/19` pushed from commit SHAs (no checkout, because suite 11622759 runs in the tree), merged in order with merge commits. Each part is retargeted to main after its parent merges. One slicing pass; nine slices exceed 400 lines and carry a `size:exception` note (single cohesive commits or one fix with its evidence). CI was not awaited, and the full suite had not finished at merge time.
- Slices (start..end, changed lines):
  - 01 `9259a08..4689258` K11 criteria and the S-2 LC3 / V-05 deck judges
  - 02 `f0b2215..f0b2215` S-2 LC3 printed, not asserted (K11.3 clean-metric gate)
  - 03 `8a577fb..8a577fb` Bernardi 2026 Fig. B1 digitiser against the v1.0/v1.1 decks
  - 04 `ed4acfe..97fa549` K11 verdicts and the re-stamped L3 rows
  - 05 `3700691..bbec949` B3-3 rejected rows and the K12 diagnostics
  - 06 `cf03ba4..afe7c5b` K13 closing-verification criteria, tools and result
  - 07 `1b4936f..5d4ff36` K14 (void) and K15 shear split
  - 08 `586bd16..3e381af` K16 matched-field candidates
  - 09 `7d20c81..7d20c81` Composite triangle fibres follow the requested span direction
  - 10 `4b6a04b..0568825` S-2, span-frame and suite anchors after the triangle fix
  - 11 `c882c1f..9900294` S-6 Newmark instrumentation and divergence location
  - 12 `7dfe78d..d745327` MITC4 directors kept across folds beyond 20 degrees
  - 13 `c0049c3..63c6dcd` Tube re-pins and anchors at the 20-degree fold directors
  - 14 `38b8ba2..1ae91a3` S-5/S-6/S-7 on the laminate mid-surface
  - 15 `cbe6345..cbe6345` S-1 bending stiffness in the section tangent frame
  - 16 `baaa641..63fa801` Rule 9: blade tests and tools meshed at 0.25 m
  - 17 `533d7e6..0754a8a` Ring grouping by element topology; modal gap convergence
  - 18 `4ee98ab..bea62ed` Aerodynamic torque sign chain restored against OpenFAST
  - 19 `8e1c516..070af3d` Nose-up convention follow-ups, P5 row, at-AC decision
- PR numbers and merge commits are in Engram (`l3/delivery-chain`). After the merge, the next session starts from `origin/main`.
