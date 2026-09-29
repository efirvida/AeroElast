# Handoff: the MITC4+/D element (SDD change `mitc4plusd-faithful`)

Written at the end of a very long session, before it degrades. **Read this file first, then
`openspec/changes/mitc4plusd-faithful/` and the two extracts.** Everything below is measured;
nothing is remembered.

## The objective (the user's words)

> *"cuando tengamos el elemento bien implementado y estable segun una formulacion concreta,
> trabajar sobre él para implementar la versión de 2025. La idea es hacer FEM en palas de
> aerogeneradores, las cuales son flexibles y geometrías curvas. Pero primero necesitamos que
> el elemento pase los test teóricos de su formulación estándar."*

So: a clean, paper-faithful **MITC4+/D** — the 2017 MITC4+ of Ko, Lee & Bathe (C&S
182:404-418) plus the 2025 penalty-free drilling DOF of Ko, Bathe & Zhang (C&S 308:107622) —
passing its own papers' basic tests, as the base for the 2025/continuum work. The final
target is wind-turbine blades: flexible, curved, warped.

## Where things stand

| | |
| --- | --- |
| **PRODUCTION** — the hybrid, `crates/aeroelast-core/src/elements/mitc4.rs` | **Works.** Twisted beam: 0.9976 / 0.9982 / 0.9986 / 0.9984 / 0.9990 against the published 0.9959 / 0.9975 / 0.9980 / 0.9972 / 0.9972. It is a **hybrid of four papers**: the Winkler & Plakomytis ERC, a `beta_w` warping penalty, selective reduced integration, a 2-DOF rotation bubble, and a shear correction factor — none of which any single paper asks for. Its plane normal was fixed this session to the paper's Eq. (10). |
| **THE NEW ELEMENT** — `crates/aeroelast-core/src/elements/mitc4_plusd.rs` | **Tier 1 GREEN** (both papers' own theoretical tests), **3 of 5** twisted-beam cells, **not production**. The flip was measured and **reverted**. |
| Rust | `cd crates && cargo test -p aeroelast-core` = **168 passed / 0 failed** |
| Python | 344 passed / 3 failed / 2 skipped. Two pre-existing; the third is a **real finding** exposed by correcting mis-sourced cells (see below). |
| The SDD change | `openspec/changes/mitc4plusd-faithful/` — explore, proposal (rev 2), spec (**rev 8**, 16 requirements / 39 scenarios), design (rev 2), tasks (**47/59**), apply-progress. Store `openspec` + Engram mirror. |

## The single open technical question (why the new element is not production)

**The three thin twisted-beam cells match; the two thick ones are 2-3x too soft.**

| case | published | hybrid | new element |
| --- | --- | --- | --- |
| thin N=8 in | 0.9959 | 0.9976 | **0.9958** |
| thin N=16 in | 0.9975 | 0.9982 | **0.9978** |
| thin N=16 out | 0.9980 | 0.9986 | **0.9987** |
| thick N=16 in (`t/L = 0.02667`) | 0.9972 | 0.9984 | **1.9623** |
| thick N=16 out | 0.9972 | 0.9990 | **2.9933** |

What has been **measured and excluded**, so do not redo it:

- **Not an `h`-power error.** The drill block scales exactly `h¹`, the same as the hybrid's,
  with an `h`-independent ratio (2.984e-4); the benchmark error goes as `h²` and the absolute
  deficit doubles as `h` halves. A constant factor, not a power.
- **Not the drill magnitude.** No drill magnitude reaches the published cells: at infinite
  drill stiffness the cells go to thick-in `1.00983`, thick-out `0.90127`, thin-in `0.96841` —
  wrong in both directions. The per-cell multipliers they would need are `≈0.5` (thin in),
  `≈30` (thick out) and `>1e4` (thick in): **mutually inconsistent**.
- **Not** the transverse-shear construction, the ADR-1 uncorrected-shear factor, the
  `ẽ^m + e^md` cross term, Eq. (21)'s centre metric, the `j0/j` ratio, the Eq. (19c) edge
  coefficients, the `t¹`/`t²` placements, `s2`, or `W_22` — each was measured and is inert or
  already correct.
- **Not a paper factor error.** An earlier note in `mitc4plusd-2025-extract.md` claimed the
  paper's Eq. (16a) dropped a `1/‖x_m^I‖` that Eq. (15a) carries. **That was my error and it is
  corrected in the extract**: the factor is cancelled exactly by the `L_l/8` of Eq. (12d)/(14c),
  so `(15a) = (16a)` term for term and the code follows the paper correctly.

The two structural facts that remain, and they point at the **thick regime itself**:

1. the paper's drill stiffness scales as `h¹` while the rotational stiffness scales as `h³`;
2. **the paper's own convergence tests are all thin (`t/L ≤ 1/100`)** while this repository's
   thick case is `t/L = 0.02667 ≈ 1/37.5`, i.e. ~2.7x thicker than anything the paper validated.

**The decision this needs is the maintainer's**, not the agent's: is the MITC4+/D's thick-cell
behaviour acceptable here, must the benchmark expectation be restated, or does the thick-regime
discrepancy point at something in the element still to be found? **Do not add a drill scale or
a penalty to close it** — neither paper has one.

## What this session established, and it is the real yield

1. **Neither the benchmark nor the paper's theory tests is sufficient alone.**
   The twisted beam **cannot see** a 2.4-degree frame error (all five numbers identical before
   and after fixing it); the paper's own basic tests **cannot see** a 13x locking. You need
   both.
2. **Three times an internal cross-check missed the bug**, because the "independent" reference
   shared the error: the frame convention (WU4's reference shared the frame), the identity lock
   in WU9b (`ke_ref` had the same double scaling), and the drill oracle (whose variant had
   become vacuous). **Every bug was found by something external**: a benchmark, the rigid-body
   oracle, or a vision read of the paper.
3. **The papers have errors of their own, and they are found by measuring, not by trusting the
   equations.** One is confirmed and relevant: the 2017 paper's Eq. (21) as printed omits the
   leading `e_rs^m|bil` term (proved by Eq. (27c)'s `(1 + a_E·rs)` coefficient) — the
   implementation must use the corrected form.
4. **A mis-sourced benchmark cell can hide a formulation error.** Correcting four of them
   (task 10.6) exposed that the distorted pinched cylinder is **6.676%** off the published
   MITC4+ cell (0.9321 vs a measured 0.9943) for the **hybrid**; the expectation had been
   0.9892, which hid it.

## The method rules that must survive

- **Vision reads only** for equations: `pdftoppm -png -r 300` (or 400) + the read tool, then
  crop with PIL. **Never** `pdftotext` for equations. The extracts record what was read.
- **Never run the twisted-beam benchmarks at N ≥ 32** (SuperLU fill-in needs tens of GB; N ≤ 16
  is ~47 s per case).
- **Never weaken a test**, never widen a tolerance, and **never add an ingredient the papers do
  not have** (no penalty, no numerical factor, no drill scale). If the only way to pass is an
  ingredient the paper lacks, **stop and report**.
- **Every requirement corrected from measurement**, not from reasoning: the zero-energy
  requirement was reasoned wrong twice and fixed from a measured table the third time.
- The Cargo workspace root is **`crates/`**, so `cd crates && cargo test -p aeroelast-core`.
  `crates/aeroelast-solvers` cannot build here (preCICE).
- Project config and per-phase rules: `openspec/config.yaml`.

## The decisions that are open

1. **The thick-cell disposition** (above) — the blocker for the new element.
2. **Delivery**: `single-pr` with an accepted `size:exception` (currently ~5000 code lines), or
   `feature-branch-chain`. The natural cut, now that Tier 1 is green for both papers:
   PR1 = WU0-WU4 (the element + its identity tests), PR2 = WU5-WU7 (Tier 1), PR3 = WU8-WU11.
3. **The rename**: the user wants the implementation named **`mitc4`** with the "plusd" only in
   the documentation — and chose to **rename only after the element passes the benchmarks**. That
   requires swapping `mitc4.rs` (the hybrid, currently production) and `mitc4_plusd.rs`, and the
   types `Mitc4Precomputed` / `Mitc4PlusDPrecomputed`.
4. **The hybrid's retirement** (WU10/S4) is gated on the new element passing **Tier 1 and Tier
   2**, and must not change the PyO3 surface.

## Where to look

| what | where |
| --- | --- |
| the new element + its 30-odd tests | `crates/aeroelast-core/src/elements/mitc4_plusd.rs` |
| the production hybrid (read-only reference) | `crates/aeroelast-core/src/elements/mitc4.rs` |
| the 2017 formulation, transcribed with page citations | `docs/formulations/mitc4plus-2017-extract.md` |
| the 2025 formulation, Eqs. (1)-(25), the Fig. 7 patch data, and the corrections | `docs/formulations/mitc4plusd-2025-extract.md` |
| the SDD artifacts | `openspec/changes/mitc4plusd-faithful/{explore,proposal,design,tasks,apply-progress}.md` and `…/specs/mitc4plusd-element/spec.md` |
| the papers | `.sources/papers/A_new_MITC4+_shell_element.pdf` (2017 formulation); `1-s2.0-S0045794917309550-main.pdf` (2017 benchmarks); `1-s2.0-S0045794924003511-main.pdf` (2025) |
| the canonical bibliography | `docs/references.md` |
| the benchmark probes | `/tmp/probe_drill.py <N> quad <in\|out> <scales…>`, `/tmp/probe_thick.py <N> <t_over_L> <in\|out> <scales…>` (rebuild them if `/tmp` was cleared) |
| the extension build | `source ~/miniconda3/etc/profile.d/conda.sh && conda activate aeroelast-dev && export PKG_CONFIG_PATH="$CONDA_PKG…"` — the full command is in `openspec/config.yaml` under `testing.build_command` |

## Session commits (newest first)

`142811d` the (15)→(16a) note (since corrected) · `5e90b4b` the t-rule double-scaling fix ·
`426249d` the production Eq. (10) normal · `e39f7c2` four benchmark cells corrected ·
`3dfa7f6` WU8/S2 the Tier-2 tests moved · `af8cf78` Tier 1b green · `f47b7f7` Tier 1a green ·
`acfca1d`, `de78c15`, `2beaef7`, `5ab74ac`, `bb0d671`, `82c7dc4`, `f7b5eb3`, `8f52502`,
`8eec86c`, `c5e800d`, `d518cb0`, `ee13216`, `b50d1f1`, `2f6cdf0`, `7f5370c`, `e0107f2`.

## The first three things to do in the new session

1. Read this file, then `apply-progress.md`'s `## WU9c` and `## WU9d` sections — they carry the
   thick-cell evidence and the excluded candidates in full.
2. Decide the thick-cell disposition (maintainer decision), because everything downstream —
   the re-flip, the rename, the retirement — is gated on it.
3. If the decision is "find it", start from the two structural facts above (the `h¹`-vs-`h³`
   stiffness mismatch and the paper's thin-only validation), **not** from a drill scale.
