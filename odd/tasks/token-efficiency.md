# Token efficiency for the fem-shell workflow

Date: 2026-10-09 · Branch: `integrate/origin-main-2026-09-30`

Trigger: the project burns ~1.1 B tokens across 13 recorded sessions while the actual
artefacts (tool results, assistant output) sum to well under 10 M tokens. The gap is
context re-read per turn, plus a handful of recurring rediscoveries that a local
artefact can remove.

Scope: this repository only. Everything here is implemented inside `fem-shell/`
(plus one Engram mirror). No global Pi setting, no `~/.pi`, no `rtk`, no shared
`gentle-*` component is touched. Global opportunities are listed as proposals at the
end and are **not** implemented.

## Baseline (measured, not estimated)

Source: Pi session logs under `~/.pi/agent/sessions/--scratch-leahk-eduardo.donestevez-fem-shell--/`,
parsed by `scripts/token_audit.py` (13 sessions, 2026-09-30 → 2026-10-09). Snapshot:
`logs/check/token-baseline-2026-10-09.json`. The numbers grow while the reading session
is live, so they are a snapshot, not a target.

| metric | value |
|---|---|
| assistant turns | 3 027 |
| input (uncached) | 80 066 855 |
| cache read | 1 011 400 832 |
| cache write | 0 |
| output | 4 132 544 |
| reasoning | 2 492 978 |
| **total tokens** | **1 095 600 231** |
| mean context per turn | 334 032 |
| peak context | 982 936 |

Tool-result payload actually produced: 6 611 191 chars ≈ 1.65 M tokens, i.e. **0.15 % of
the billed volume**. The other 99.85 % is repeated context, not new information.

Where the repeated context comes from:

1. **Session length.** Cost ≈ Σ context per turn, which grows roughly linearly inside a
   session, so a session is ≈ quadratic in its turn count. The two marathons
   (944 and 570 turns) alone account for 649 M of the 1 011 M cache reads. A 150-turn
   session costs ~30 M; six of them cost less than one 944-turn session doing the same work.
2. **Fixed per-turn floor.** The system prompt is 18 246 tokens/turn:
   addendum 6 392 · skills 3 721 · `AGENTS.md` 3 757 · rules 2 962 · tools 837 · docs 365.
3. **Memory recall is expensive per call.** `mem_search` 20 calls → 176 701 tokens
   (8 835/call); `mem_context` 10 calls → 20 546. Broad queries return whole observations.
4. **Repeated rediscovery.** `bash` 2 307 calls / 918 049 tokens and `read` 294 calls /
   352 402 tokens are not individually large; they are numerous, and a large fraction of
   them re-derive facts that are stable (environment bootstrap, test command, landmark
   locations, known-red tests).

How the 1 016 M cache reads divide up:

| component | tokens | share |
|---|---|---|
| system-prompt floor (≈19 370 tok × 3 081 turns) | ≈60 M | ≈5.9 % |
| accumulated conversation history | ≈956 M | ≈94.1 % |
| tool payload produced (once) | 1.7 M | 0.15 % |

So trimming the always-on prompt is a ~6 % game at best, while the length of the context
that history builds is the whole game. That is why the plan spends its per-turn budget
carefully and puts the detail in a skill that is loaded on demand.

## Findings that shape the plan (verified 2026-10-09)

- **The documented test command does not run in a fresh shell.** `AGENTS.md` gives
  `python -m pytest tests/ ...`; without the GCC 14 runtime the collection fails with
  `ImportError: /lib64/libstdc++.so.6: version 'CXXABI_1.3.15' not found`. Reproduced:
  3 collection errors before the bootstrap, 66 passed after. Shell state does not survive
  between `bash` calls, so **every** aeroelast command needs the bootstrap, not just the
  first one.
- **All three `--ignore=` flags in `AGENTS.md` are dead.** `tests/test_blade_mesh.py`,
  `tests/test_rotor_inertial.py` and `tests/test_vol_mesh.py` do not exist. The stale
  flags are harmless but they hide that the command has not been re-derived since the
  2026-09-30 merge, which is exactly what `docs/origin_main_integration_2026-09-30.md`
  already says about the `vol_mesh` one.
- **The full suite is not a smoke test.** 20–26 min wall time on 2026-09-30/10-01, and
  192 000–267 000 warnings. A re-run on the same tree on 2026-10-09 reached only 47 % of
  the collected tests after 33 min (≈70 min end to end): the integration and validation
  directories dominate the first half, so the figure depends on the node and on which
  heavy cases run. The failure count also moves between runs (4 → 13 failures between
  `after-full9` and `after-full10`). A "run everything and read the log" gate is both slow
  and unstable.
- **`ruff check` on the repo emits 161 801 bytes (313 pre-existing errors).** With
  `--statistics` the same information about *which rules* fail is 1 134 bytes (143×).
  pi-lens already runs Ruff per edited file, so the global run is almost pure duplication.
- **The four known-red tests in the integration doc are still red (three re-verified
  here):** `test_corotational_is_frame_objective_tl_is_not`,
  `test_dtube_tip_deflection_matches_beam_theory`,
  `test_dominant_component_matches_beam[LC1_gravity-utd]`, and the UL-elastica accuracy
  case. They are documented decisions, not regressions, but nothing in the run output says
  so — an agent re-diagnoses them every time.
- **`pi-rtk-optimizer` already filters tool output** (project scope: 3 450 commands,
  7.4 M tokens saved, 72.8 %). The remaining wins are not about making individual outputs
  shorter; they are about not producing or re-reading them at all.

## Decisions

- **D1 — one bootstrap script, used by everything.** `scripts/aeroenv.sh` resolves the
  GCC 14 runtime, the venv and `LD_LIBRARY_PATH`, and either exports into the current
  shell or execs a command. Hard-coded paths are only a fallback after `module load`.
- **D2 — a two-mode gate.** `scripts/check.sh quick` (≈13 s: Ruff on changed files plus a
  fast physical subset) is the pre-report check; `scripts/check.sh full` runs the
  documented suite once, writes the log to `logs/check/`, and prints only the summary plus
  non-known-red failures. Quick is a smoke gate and **never** a substitute for the
  physical validation anchors (S-2 / S-4 / S-5 / S-7).
- **D3 — know your red.** `scripts/known-red-<mode>.txt` records the current known-red set with
  the revision and date it was verified at; `quick` and `full` keep separate registries so a
  measurement run cannot clobber the gate. `check.sh` exits non-zero only for failures outside
  the mode's file.
- **D4 — machine-readable pytest output.** The gate runs pytest with
  `-o addopts="" -q --tb=line --disable-warnings`, because the project `addopts` (`-v`)
  triples the output of a passing run (measured: 7 558 → 1 684 bytes on the same seven
  files) and the warning summary adds nothing a gate acts on. `pyproject.toml` is left
  alone so human runs stay verbose.
- **D5 — progressive disclosure instead of a fatter `AGENTS.md`.** The detail lives in the
  Pi project skill `.pi/skills/fem-shell/SKILL.md`, which is advertised by name and loaded
  only when the task matches. `AGENTS.md` gets ~15 lines: the bootstrap rule, the two
  canonical commands, the handoff rule, and pointers.
- **D6 — measurement before and after.** `scripts/token_audit.py` reads the Pi session
  logs of the current working directory and prints the same table as the baseline, so the
  effect of these changes is checkable with the same instrument.

## Tasks

1. [x] **WU-1** `scripts/aeroenv.sh` + verified import in a clean shell.
2. [x] **WU-2** `scripts/check.sh` quick/full + the mode-scoped `scripts/known-red-*.txt`; quick
   baseline recorded (66 passed, 1 known-red, 1 xfailed in 12.45 s).
3. [x] **WU-3** `scripts/token_audit.py` reproducing the baseline table.
4. [x] **WU-4** `.pi/skills/fem-shell/SKILL.md`.
5. [x] **WU-5** `AGENTS.md` section: bootstrap rule, canonical commands, session hygiene.
6. [x] **WU-6** Evidence recorded here and in Engram; global proposals listed, not applied.
7. [x] **WU-7 (second pass)** Renamed the context file to `AGENTS.md` (the cross-tool
   convention, so Copilot/Cursor/Codex read the same file) and rewrote it: 264 → 173 lines,
   15 030 → 12 550 chars of `project_context`, while *adding* the physics-discipline section,
   the validation-anchor table and the traps table. Verified factual corrections in the same
   pass: `rotor.py:1827` → `:2021`, the `dynamic_newmark.rs:492-496` KSP claim →
   `solve_with_cached_ksp`, the solver class names (`StaticLinearSolver` → `StaticLinearSolver`,
   `DynamicNewmarkSolver`/`LinearDynamicSolver`, `LinearDynamicFSISolver`), the removed `SOLID`
   family and the wrong "max stride 6" rule (real rule: `MeshAssembler._FAMILY_PROPERTIES`,
   mixed meshes infer from `dofs_per_node == 6`). References to the old filename updated in
   `tests/test_iea15mw_v05_structural_properties.py` and
   `docs/blade_input_divergence_utd_vs_official.md`.
8. [x] **WU-8 (make the rules passive)** Turned the usage tips into mechanisms instead of
   advice, at three layers: `.pi/extensions/session-guard.ts` (project extension, zero prompt
   cost — status-bar budget meter at 150 k / 250 k / 500 k plus one-off notifications, and a
   `/handoff` command), `## Session economy` in `AGENTS.md` (six agent-side rules), and
   "The session contract" in the skill (open / during / close / resume templates).

## Validity and invalidation

Reused artefacts are only useful while they are still true. Each one carries its own
invalidation trigger:

| artefact | valid while | invalidated by |
|---|---|---|
| `scripts/known-red-quick.txt` | the recorded revision and `.so` are unchanged | any `crates/` change, a `maturin develop` rebuild, or a suite run whose new failures appear outside the file — re-derive with `scripts/check.sh quick --record` |
| `.pi/extensions/session-guard.ts` | the thresholds are a session-length heuristic | `ctx.getContextUsage()` stops returning `tokens`, or the Pi extension API is renamed |
| `logs/check/*.log` | same revision | same as above |
| `.pi/skills/fem-shell/SKILL.md` landmark map | paths exist | a rename/move on the referenced paths (checked by the quick gate) |
| `scripts/token_audit.py` output | Pi keeps the session-log format (v3) | a Pi major version change |
| the baseline table above | sessions recorded before 2026-10-09 | by design: it is a baseline, not a target |

Rules that keep the local artefacts from rotting:

- Never edit `known-red-quick.txt` by hand to make a run green. Re-derive it with `--record` and
  say which run produced it.
- A red test that is not in the registry is treated as a regression, always.
- The quick gate is expected to run in one `bash` call. If it needs a second call to
  interpret, the gate is wrong, not the run.
- Evidence that only exists in a session transcript is not evidence. Verified numbers go
  into this file, with the command that produced them.

## Evidence

```bash
# 66 passed, 1 failed (known-red), 1 xfailed — 12.45 s, 635 bytes of output
scripts/check.sh quick

# the documented full suite, quiet, log kept under logs/check/
scripts/check.sh full

# the baseline table at the top of this file
python3 scripts/token_audit.py

# the environment failure the bootstrap prevents
bash -c 'python -c "import aeroelast"'   # ImportError: CXXABI_1.3.15 not found
scripts/aeroenv.sh python -c "import aeroelast; print('ok')"
```

### Verified, with the command that verified it

| claim | check | result |
|---|---|---|
| bootstrap fixes the import | `scripts/aeroenv.sh --check` | `ok`; `import aeroelast` succeeds in a clean shell |
| bootstrap is idempotent and pipefail-safe | `bash -c 'set -uo pipefail; source scripts/aeroenv.sh'` | rc=0, no duplicated `LD_LIBRARY_PATH` |
| the gate passes with only documented-red failures | `scripts/check.sh quick` | `1 failed, 66 passed, 1 xfailed in 12.7 s` → `OK`, exit 0 |
| the gate fails on an undocumented failure | same run with `known-red-quick.txt` moved away | `NEW-FAIL` + exit 1 |
| the registry is re-derivable | `scripts/check.sh quick --record` | 1 entry, stamped with revision `bd4cc05` |
| quick and full registries stay separate | `scripts/check.sh --list` | each mode reads only its own `known-red-<mode>.txt` |
| Ruff-on-diff branch works | `scripts/check.sh quick` with an untracked `scripts/*.py` | `ruff: clean` |
| the parser sees the structure I assumed | `scripts/token_audit.py` | reproduces the baseline table above |
| snapshot/compare round trip | `--save` then `--compare` | 0 delta on unchanged data |
| the skill is actually advertised | headless `pi -p` run, then the session log's `skills` section | `fem-shell` present with name, description and location |
| `AGENTS.md` is loaded once, not twice | same headless run, `project_context` section | 12 550 chars, `# AGENTS.md` occurs exactly once |
| the 2nd-pass corrections hold | symbol/line greps against the tree | `_step_cb` at `rotor.py:2021`, `solve_with_cached_ksp`, `_FAMILY_PROPERTIES`, no `SOLID` family |
| the AGENTS.md pointers resolve | existence scan of every backticked repo path in `AGENTS.md` | 0 dead pointers outside the explicitly historical ones |

### Measured effect on the always-on prompt

The system prompt is re-sent every turn, so additions to it are a per-turn tax. The
skill body deliberately lives outside it:

| stage | `project_context` | `skills` (advertisement) | total prompt |
|---|---|---|---|
| original `CLAUDE.md` | 15 030 | 14 887 | 72 985 |
| 1st pass | 16 640 | 15 542 | 74 121 |
| 2nd pass (`AGENTS.md`) | 12 550 | 15 542 | 70 031 |
| + `## Session economy` + `session-guard` extension | 14 065 | 15 542 | **71 546** |

The `session-guard` extension adds **zero** characters to the prompt: the 1 515-char delta is
exactly the `project_context` growth, and the extension is code plus UI, not context. Net
against the untouched original the file is ~6 % smaller while carrying the physics-discipline
section, the traps table, the validation-anchor table and the session rules. The skill body
(9.1 kB) would cost ~2 300 tokens/turn if it were pasted into the context file instead of being
loaded on demand.

### Cost of the change itself

One headless `pi -p` run (~20 k input tokens) was spent to prove the skill is discovered.
That is the only model call this work made; everything else is local execution.

## Not yet measured

The honest gap: there is **no before/after on real work yet**. `scripts/token_audit.py`
now exists to produce it. To close it, record a snapshot before a normal piece of work and
compare afterwards:

```bash
python3 scripts/token_audit.py --save logs/check/snap-before.json
# ... do the work ...
python3 scripts/token_audit.py --compare logs/check/snap-before.json
```

Watch `cacheRead/turn`, not the totals: the totals grow simply because the current session
is included.

The `check.sh full` mode was validated as far as launch mechanics (Ruff-on-diff, the quiet
pytest invocation, the log under `logs/check/`), but no run has completed yet: on 2026-10-09 a
login-node run was killed at ~17 %, and a second attempt started by the same invocation died the
same way. Anything involving the full suite on this cluster belongs in a SLURM job. Re-derive
`known-red-full.txt` from that job's log.

## Observations with limits — what the logs do and do not show

- **Prompt caching already works.** cacheRead is 92 % of billed volume, cacheWrite is 0.
  The lever left is volume, not hit rate; do not go hunting for cache-miss waste.
- **Retries are a wall-time cost, not a demonstrated token cost.** 41 of 3 116 assistant
  messages (1.3 %) ended with `stopReason: "error"`, and they carry no `usage` block at
  all — 0 cacheRead attributed to them. Median assistant call is 11.9 s, p90 25.8 s, max
  58 s, so the 41 dead turns cost up to ~40 minutes of wall time, but the log cannot show
  a token charge for them. Treat "retries waste tokens" here as unproven.
- **Tool errors are rare and concentrated.** 109 of 3 518 tool results are `isError`
  (3.1 %), and the largest class is the environment failure the bootstrap removes.
- **Observability limits.** Pi logs `input`, `output`, `cacheRead`, `cacheWrite`,
  `reasoning` and `totalTokens` per assistant message, plus per-tool payload sizes. It does
  not log token prices, provider-side cache write charges, or the cost of aborted requests,
  so nothing here converts to money, and the char/4 token estimates for tool payloads are
  estimates, not measurements.

## Out of scope — global proposals (needs explicit authorization)

These would help this project but live outside it, so they are **not** implemented:

1. `rtk init -g` — install the rtk hook for automatic interception (currently every
   `bash` call prints a 20-token "No hook installed" banner; 2 307 calls in the baseline).
   Global tool configuration; also affects every other project.
2. Pi compaction tuning (`compaction.keepRecentTokens`, `reserveTokens`) — the marathons
   reached 983 k tokens before compaction. A project `.pi/settings.json` override is
   technically local, but it changes the reliability of every other session in this repo,
   so it needs the maintainer's decision first.
3. A cheaper model profile for mechanical steps (finding files, running the gate) — a
   provider/settings decision, not a repository one.
