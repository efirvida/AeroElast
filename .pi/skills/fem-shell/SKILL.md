---
name: fem-shell
description: "Work inside the fem-shell / aeroelast repository: run the aeroelast test suite, a single test, or the Ruff gate; import or rebuild the Rust extension with maturin; diagnose an import failure such as `CXXABI_1.3.15 not found` or a stale `_aeroelast.so`; interpret the known-red tests; find a solver, config, validation anchor, or theory document; resume interrupted work from an `odd/tasks/` feature document. Use whenever a task runs code, tests, or builds in this repository."
---

# fem-shell workflow

This repository is expensive to rediscover. These are the rules that stop the
rediscovery; everything here was measured, not assumed (`odd/tasks/token-efficiency.md`).

## 0. Bootstrap before any command that imports the package

The Rust extension is built against GCC 14. The system `libstdc++` has no `CXXABI_1.3.15`,
so a bare import fails with:

```
ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
```

Shell state does **not** survive between tool calls, so this is needed before *every*
aeroelast command — not just the first one in a session.

```bash
scripts/aeroenv.sh python -c "import aeroelast"     # one-shot: run a command
scripts/aeroenv.sh --check                          # one line: is the env sane
source scripts/aeroenv.sh                           # export into the current shell
```

Treat a bare `python -m pytest` or a bare `import aeroelast` failure as *your* mistake,
not a broken tree. Do not start a diagnosis loop over it.

## 1. Canonical commands

```bash
scripts/check.sh quick                     # ~13 s gate: Ruff on changed files + fast subset
scripts/check.sh full                      # the documented suite (~25 min), log kept
scripts/check.sh --list                    # the selections and the known-red registry
scripts/check.sh quick --record            # re-derive known-red-quick.txt from a real run

ruff check --statistics                    # 1 kB of rule counts (bare `ruff check` is 161 kB)
ruff check path/to/changed.py              # what pi-lens already does on every edit

cd crates/aeroelast-py && maturin develop --release   # after any crates/ change
```

Rebuild rules: `pip install -e .` does **not** recompile Rust; use `maturin develop
--release` from `crates/aeroelast-py/` (`scripts/aeroenv.sh` exports the `VIRTUAL_ENV`
maturin requires). On the SDumont cluster the module/libstdc++ handling is in
`.github/skills/build-aeroelast/SKILL.md`.

## 2. What "verified" means here

`scripts/check.sh quick` is a smoke gate, not validation, and it says nothing about
physics. Escalate by what changed:

| change | minimum evidence |
|---|---|
| docs, comments, pure refactor | `check.sh quick` |
| element / DOF / convention / BC | quick + the matching unit test + one CalculiX or beam reference |
| solver formulation (rotor, corotational, inertial) | the relevant anchor below, with numbers, not "it runs" |
| anything touching `crates/` | `maturin develop --release` first, then a full-suite measurement before a merge |

Physical anchors: **S-2** static prescribed, **S-4** rotating modal vs OpenFAST v5 MBC3,
**S-5** one-way 100 s rated run, **S-7** torsion ratio, plus V-02 natural frequencies.
An anchor that moves is a finding to report, never a tolerance to widen.

## 3. Known-red tests

`scripts/known-red-quick.txt` lists tests that are red by documented decision. `check.sh` prints
them as `KNOWN-RED` and exits non-zero only for anything else. The current entry is
`test_corotational_is_frame_objective_tl_is_not` — deliberate: the MITC4 path of
`assemble_kt_corotational` uses the upstream total-Lagrangian tangent, so the property the
test asserts no longer differs (`docs/origin_main_integration_2026-09-30.md` §C).

Rules: never hand-edit the registry to make a run green; re-derive with `--record` and say
which run produced it. A red test outside the registry is a regression. A registry entry
that now passes is stale — re-record.

## 4. Landmark map

Read the file before the directory. Whole-doc reads of `odd/tasks/*.md` (up to 156 kB) are
the single worst read in this repo; use the section you need.

| need | path |
|---|---|
| solver dispatch, YAML → solver | `src/aeroelast/solvers/fsi/runner.py` |
| production rotor solver | `src/aeroelast/solvers/fsi/rotor.py` |
| inertial variant (WIP) | `src/aeroelast/solvers/fsi/rotor_inertial.py` |
| config dataclasses, `SolverType` | `src/aeroelast/core/config.py` |
| Rust per-iteration FSI loop | `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` |
| PyO3 crate (maturin root) | `crates/aeroelast-py/Cargo.toml` |
| stress / ply recovery | `src/aeroelast/postprocess/stress_recovery.py` |
| theory: FSI rotor formulation | `docs/formulations/teoria_formulacion_fsi_rotor.md` |
| theory: frames, inertial forces | `docs/formulations/solvers.md` |
| shell/element formulation | `docs/formulations/shell-elements.md`, `mitc4plusd-2025-extract.md` |
| before/after of the 2026-09-30 merge, anchors, open decisions | `docs/origin_main_integration_2026-09-30.md` |
| validation closures and re-run triggers | `docs/validation_closures.md` |
| consolidated rotor validation report | `docs/validation_corotational_rotor_solver.md` |
| validation matrix (must stay at 0 errors) | `tools/validation_matrix.py check` |
| official IEA 15 MW blade input | `tests/IEA-15-240-RWT.yaml` (never the UTD xlsx in new work) |

## 5. Session hygiene — the largest lever

Cost is roughly the sum of the per-turn context, so it grows with the **square** of the
session length. Measured here: 13 sessions, 3 000 turns, 1.0 B cache-read tokens, of which
tool payload was 0.15 %. The marathons (944 and 570 turns) were 64 % of it.

- Long work goes in bounded sessions. Before the context passes ~200 k, write the state
  into the feature document and start a fresh session.
- The handoff is `odd/tasks/<feature>.md` (trigger, scope, findings with `path:line`,
  tasks as checkboxes, evidence). Resume from it plus one narrow `mem_search` — not from a
  transcript.
- Do not re-read what a subagent already summarized; spot-check at most one detail.
- Bounded tool output: `--stat`, `tail`, `wc -c`, `grep -c`, `-q`, and redirect long logs
  to a file, then read the summary line. `grep -q` over `grep` when you only need a verdict.
- Full-suite runs belong in `logs/check/`, once per measurement point, not once per edit.

### The session contract

`.pi/extensions/session-guard.ts` puts the budget in the status bar every session
(`ctx 132k · 41t · ok` → `cerrá pronto` → `CERRÁ + /handoff` → `CARO: CERRÁ YA`, at 150 k /
250 k / 500 k). It costs no tokens: it is UI, not context. Use it as the cue, and shape the
session like this:

**Open** with one objective, the evidence level, and the boundaries:

> Una sola cosa: `<objetivo>`. Alcanza con `scripts/check.sh quick` + `<test o ancla>`. No toques `<X>`.

**During**: bounded tool output, no re-reading, and stop to propose a split if a second
workstream appears.

**Close** at a work-unit boundary, before the context passes ~200 k:

```text
/handoff <feature>
```

That writes the handoff and returns the resume recipe. If you close by hand, the next session
needs exactly four things:

1. objective and current state;
2. decisions taken and why (including the rejected options);
3. verified evidence — command plus result, never narrative;
4. the next concrete step.

**Resume** in a fresh session from `odd/tasks/<feature>.md` + one narrow `mem_search` +
`git log --oneline -5`. Never resume from a transcript.

```bash
python3 scripts/token_audit.py                      # this project, from the Pi session logs
python3 scripts/token_audit.py --save logs/check/snap.json
python3 scripts/token_audit.py --compare logs/check/snap.json
```

## 6. Memory protocol

- One narrow query beats one broad query: `mem_search` has cost up to 8 800 tokens per call
  here. Prefer `topic_key`, a distinctive identifier (`S-4`, `K_G`, a file name), and a
  `limit`.
- Stable, evolving knowledge gets a `topic_key` so recall updates instead of duplicating:
  `fem-shell/<subsystem>/<decision>` (e.g. `fem-shell/rotor/k-g-assembly`).
- Verified numbers and closed decisions go into `odd/tasks/<feature>.md`; memory mirrors it.
  Evidence that lives only in a transcript is not evidence.

## 7. Validity of reused results

| artefact | invalidated by |
|---|---|
| `scripts/known-red-quick.txt`, `logs/check/*.log` | any `crates/` change, a `maturin` rebuild, or a new failure outside the list |
| `.pi/skills/fem-shell/SKILL.md` landmark map | a rename or move of the referenced paths |
| `scripts/token_audit.py` output | a Pi session-format change |
| cached validation numbers in `docs/` | a change to the mesh, the blade input, the solver formulation, or `K_G` assembly |

When one of those changes, re-derive the artefact instead of trusting it.

## 8. Out of scope here

Global or cross-project changes (`rtk init -g`, Pi compaction settings, provider/model
choices) are proposals in `odd/tasks/token-efficiency.md`, not this repository's business.
Never apply them from a task in this repository.
