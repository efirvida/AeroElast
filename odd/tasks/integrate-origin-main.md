# Feature: integrate origin/main into the validation line

**Status**: in progress
**Opened**: 2026-09-30
**Branch**: `integrate/origin-main-2026-09-30` (cut from `integrate/main-into-rust` @ d571917)
**Engram mirror**: `odd/integrate-origin-main/tasks`

## Goal

Bring the 169 commits that `origin/main` (dffacab, Eduardo Firvida) added after our
last merge (`e2bc081`, which took `origin/main` @ `7a2196f`) into the IEA 15 MW
validation line, and measure the effect on the validation suite with a documented
before/after.

## Why

`origin/main` carries the reviewed MITC4+/D element, the span-relative ply-angle
fix, the "uncorrected section shear" material fix, the CCX writer id-scheme fix,
and its own CalculiX parity suite. Our line carries the IEA 15 MW S/V-series
validation suite, the campaign tooling, the static K_G path and the CCX blade
column. The user wants the structural improvements on the branch where the
validation runs.

## Topology (verified 2026-09-30)

```
main 2cf3df5 (2026-02-19)  ── common ancestor of both lines
  ├── rust_implementation 7df69dc   +114 commits  → our validation line
  └── origin/main        dffacab    +435 commits  → reviewed element line

HEAD integrate/main-into-rust d571917 = rust_implementation line
  + merge e2bc081 (origin/main @ 7a2196f)
  + 7 commits (CCX column, blade structural matrix)
merge-base(HEAD, origin/main) = 7a2196f ; origin/main advanced 169 commits since.
```

## Decisions taken by the user (2026-09-30)

1. **Strategy**: new branch from HEAD, then merge `origin/main` into it. Keep
   `integrate/main-into-rust` untouched.
2. **3D solids**: accept `origin/main`'s removal (Rust `7b295f3`, Python
   `9230ea2`). Solid-dependent tests become documented skips/xfails.
3. **Execution**: login node (sdumont17). Full suite first, then the IEA 15 MW
   subset, for both the before and the after measurement.

## Non-goals

- No push, no PR, no merge to `main`. Delivery stays a user decision.
- No re-derivation of campaign results; this feature measures the suite, not the
  article numbers.
- No attempt to keep 3D solid support alive on our side.

## Environment recipe (verified)

```bash
module load gcc/14.2.0_sequana          # do NOT pipe this call: the pipe subshell loses it
export LD_LIBRARY_PATH=/scratch/app_sequana/gcc/14.2.0/lib64:/scratch/app_sequana/gcc/14.2.0/lib:$LD_LIBRARY_PATH
/scratch/leahk/eduardo.donestevez/venv/bin/python -c "import aeroelast"
# Rust rebuild: maturin develop --release in crates/aeroelast-py (VIRTUAL_ENV set); pip install -e does NOT recompile
```

Suite command (`docs/validation_closures.md`):

```bash
python -m pytest tests/ -q --tb=short \
    --ignore=tests/test_blade_mesh.py \
    --ignore=tests/test_rotor_inertial.py \
    --ignore=tests/test_vol_mesh.py
```

Historical baseline: 923 passed / 26 skipped / 4 xfailed / 0 failed.

## Tasks

- [x] T1 — Cut the integration branch and commit the pending validation work.
      Commits: `df75313` (untrack the .atl tool state), `940ab31` (frame-fix
      erratum), `2a348f8` (blade CCX parity test), `5d179a5` (this plan).
- [x] T2 — Rebuild the Rust extension at pre-merge HEAD and capture the BEFORE
      measurement (full suite + IEA 15 MW subset).
      Rebuild: `maturin develop --release` OK in 1m59s (11 warnings, 0 errors).
      Full suite: **933 passed, 26 skipped, 9 xfailed, 6 failed** in 17m45s.
      Subset (S0-S7, sx, V-01..V-05, composite CCX, blade CCX):
      **68 passed, 6 failed, 5 xfailed** in 10m06s.
      Logs: `$SCRATCH/tmp/odd-integrate-origin-main/{before-full,before-subset}.log`.
      The 6 failures are all on the CalculiX parity path and share one root
      cause in `write_ccx_mesh`:
        - `tests/test_composite_ccx_parity.py[unbalanced_45_0s_bend-...]` ->
          `writers.py:122 _build_angle_bucket_sets` ->
          `IndexError: index 7500 is out of bounds for axis 0 with size 125`
          (the writer mixes node-id spaces).
        - the 5 `tests/test_blade_ccx_parity.py` cases -> CalculiX exits 201 with
          `*INFO in gen3dnor: in some nodes opposite normals are defined`, on a
          deck CCX reads as 696743 nodes for a 3333-element mesh.
      This is exactly the area `origin/main` touched (`930d055 fix(ccx): label
      nodes, elements and sets consistently for any id scheme`), so it is the
      headline candidate for the AFTER comparison.
- [x] T3 — `git merge --no-commit origin/main` and resolve the 42 conflict hunks
      in 19 files per the resolution table below.
      Commit: `f5f92f7`.  42/42 hunks resolved with a per-hunk, auditable choice
      (`$SCRATCH/tmp/odd-integrate-origin-main/resolve.py`); the four union
      artifacts it left were repaired in `f93fc5a` and `28339a4`.
- [x] T4 — Apply the solid-removal consequence: skip/xfail the solid-dependent
      tests with a documented reason.
      No skip was needed: upstream's `7e13d1e` deleted the solid test paths
      (`test_solid_elements.py`, `test_vol_mesh.py`, `test_beam_4cases_parity.py`)
      along with the Rust and Python support.  The two `ccx_bin_or_skip` helpers
      from upstream's `conftest.py` were kept (a union) because four test files
      import `ccx_bin_or_skip` from it.  The surviving CCX tests read
      `tests/_ccx_io.py`.
- [x] T5 — Rebuild the Rust extension post-merge and fix any build breakage.
      Commits: `f5f92f7` (core `centrifugal_load`, `assemble_geometric_k_from_disp`,
      `update_node_coordinates`, `assemble_kt_corotational`; py `fsi.rs` rebuilt from
      our five runner signatures; `alpha_prev`/`velocity_write_data` restored) and
      `28339a4` (py binding for `assemble_kt_corotational`).  Build green.
- [x] T6 — Capture the AFTER measurement (full suite + IEA 15 MW subset).
      Full suite: **905 passed, 15 failed, 22 errors, 38 skipped, 6 xfailed**
      (`after-full3.log`); the earlier passes are in `after-full.log` (41 failed,
      22 errors) and `after-full2.log` (22 failed, 22 errors) before the union
      artifacts and the gravity-load unit bug were fixed.
- [x] T7 — Write the before/after comparison on the validation subset, including
      which previously-failing or xfail items closed.
      Commit: `docs/origin_main_integration_2026-09-30.md`.
- [x] T8 — Update `docs/validation_closures.md` and `CLAUDE.md` with the new
      state and the inherited limitations.
      Both updated in the same commit as T7.

## Outcome

See `docs/origin_main_integration_2026-09-30.md`.  In one line: the merge lands the
reviewed element and the CCX writer fix, and while measuring it we also found and fixed a
real composite gravity-load unit bug and four merge artifacts — but 28 CalculiX tests are
blocked by one degenerate blade-mesh element, and six of our validation anchors moved
(V-02 2F, S-4 rotating 1F, S-7 torsion, S-6 NaN, box EI, D-Tube, elastica).

## Conflict resolution policy (dry-run inventory)

`git merge-tree --write-tree HEAD origin/main` → 42 hunks, 19 files + 2 modify/delete.

| Area | Files | Policy |
|---|---|---|
| Element kernels | `mitc3.rs` (3), `mitc4.rs` (6) | **Theirs**. Upstream is the reviewed MITC4+/D line; our mitc3 hunks are comments/whitespace only. |
| Rust bindings/FFI | `aeroelast-py/src/lib.rs` (1) | **Theirs**, then re-add only what our Python layer still calls after the solid removal. |
| FSI Rust loop | `rotor_fsi.rs` (1) | **Theirs** unless it drops our K_G/K_SP call sites — verify. |
| Python assembler | `assembler.py` (2) | **Per hunk**: theirs for structure, but our orthotropic→`RustLaminate` path and the static K_G call sites must be checked before discarding. |
| Mesh generators/writers | `generators.py` (2), `writers.py` (5) | **Per hunk**: upstream for the CCX id-scheme and region ELSETs; ours for anything only our validation suite consumes. |
| CLI / package surface | `aeroelast.py` (1), `__init__.py` (1) | **Theirs** for style; keep our GCC-runtime guard if it disappears. |
| Monitor / rotor | `power.py` (1), `rotor.py` (1) | **Theirs** unless it drops a K_G/`free_dofs` call site. |
| Tests | `conftest.py`, `test_composite_beam_parity`, `test_mass_matrix_validation` (3), `test_mitc3_benchmarks`, `test_rotor_physical_consistency` (10), `test_shell_comprehensive` | **Per hunk**: prefer upstream expectations where the element behavior changed; keep our reference artifacts where upstream has none. |
| Docs | `.github/copilot-instructions.md` (2) | **Theirs**. |
| Tool state | `.atl/skill-registry*.md` | **Deletion** (upstream untracked them and `.gitignore` ignores `.atl/`). |

## Inherited limitations to record (not regressions)

- 3D solid elements gone (Rust + Python): `tests/test_solid_elements.py`,
  `tests/test_vol_mesh.py` lose their backend.
- Upstream has no `assemble_geometric_k_from_disp` / `centrifugal_load` Rust
  bindings; our static K_G path must survive the merge or S-4 loses its meaning.
- The bend-twist 3.23x and the "3-12x over-twist" hypotheses get a second chance
  via `9de3731` (span-relative ply angles) and `8cbfc0b` (uncorrected section
  shear) — verify, do not assume.
