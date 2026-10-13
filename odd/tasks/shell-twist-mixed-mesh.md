# Feature: is the shell's 2.86x twist (vs the beam, same load) a mixed MITC3/MITC4 effect?

Spun off from `coupled-divergence.md` T7d/T7e (2026-10-13). Worktree `rotor-frame`, branch
`fix/bem-participant-rotor-frame`. Read-only investigation: **no element code changes without
the maintainer's OK.**

## Question

Under Zhou's published loads the production shell turns +5.96 deg, the standard beam +2.08 deg
(2.86x). T7e left two readings open: (i) the shell is right and the beam cannot represent section
distortion; (ii) the shell over-twists through an effect the beam lacks. The blade mesh mixes
MITC4 quads with MITC3 triangles (root, trailing edge, tip), and two open defects live exactly
there: #34 (quads and triangles use different fibre-axis conventions) and #38 (MITC3+ drilling
penalty dimensionally inconsistent). They are candidates for (ii).

## What the code does (read 2026-10-13, `mitc3.rs:166`, `assembler.rs:234`)

- Quads: `is_composite && node_ids.len() == 4` gets the span-relative angle offset `phi`; MITC4
  consumes the laminate as given (axis 1 = span).
- Triangles: no `phi`; `rotate_constitutive_to_local` sets axis 1 = global X projected on the
  element, falling back to global Z when X is normal to it.
- Consequence on a blade (span = Z): a triangle on a cap/skin panel (normal ~ Y) has axis 1
  across the chord, 90 deg from the quads around it.

## Design correction (before any run)

The planned thin-walled **tube in torsion** cannot see #34: `Qbar66(theta)` depends on
`sin^2(2 theta)`, so a 90 deg axis swap leaves the Bredt rate unchanged (checked numerically:
ratio 1.0000 for [+-45]s and [+-30]s). The discriminants are axial, bending and bend-twist
coupling (E1 vs E2; `D16` changes sign).

## Pre-registered measurements

- **E1 element-level** (`tools/diagnose_mixed_mesh_consistency.py --element`): cantilever
  composite strip, span along Y, quad axis 1 = span. Meshes: all-quad (44), all-tri (33), rows
  alternating. Laminates: UD [0]4 and unbalanced [45,0]s; isotropic control (4 vs 3).
  Metrics: tip axial uy, tip bending uz, tip twist about y. Predictions from CLT only:
  all-quad matches the closed form (axial) within 5 %; if #34 is active, tri/quad equals the
  CLT ratio of the laminate rotated 90 deg (printed by the tool before solving); isotropic
  tri/quad within 1 %. tri/quad ~ 1 on the composite means #34 is not active in this path.
- **E2 blade geometry** (`--blade`): triangles by r/R, area share, and the angle between the
  triangle's effective axis 1 and the span-projected axis (the quad contract). Reading rule:
  area-weighted error and share of tri area inside 0.3-0.9 R decide whether E3 is worth asking.
- **E3 (not run, needs OK)**: Zhou-load shell run with triangles forced to the span contract.

## Tasks

- [x] T1 Tool written, predictions printed before solves (`tools/diagnose_mixed_mesh_consistency.py`)
- [x] T2 E1 run, log in `logs/check/mixed_mesh_e1.log` (git-ignored; numbers below)
- [x] T3 E2 run (blade mesh), log in `logs/check/mixed_mesh_e2.log`
- [x] T4 Handoff for the analysis session (below)

## Results (2026-10-13, one run each, nothing tuned)

**E1: #34 is real at element level and the CLT prediction, written before solving, holds.**
Triangle/quad tip-response ratio, laminate axis 1 = span for quads:

| laminate | metric | predicted | measured | delta |
|---|---|---|---|---|
| UD [0]4 | axial | 12.000 | 11.842 | -1.3 % |
| UD [0]4 | bending | 12.000 | 12.002 | 0.0 % |
| [45,0]s | axial | 4.938 | 5.082 | +2.9 % |
| [45,0]s | bending | 1.950 | 1.866 | -4.3 % |
| [45,0]s | bend-twist | -2.438 | -2.295 | -5.9 % (sign flips, as predicted) |
| isotropic control (3 vs 4) | axial / bending | 1.000 | 0.996 / 0.999 | within 1 % |

All-quad matches the closed form (UD axial 1.015, bend 0.999; [45,0]s 0.972, 0.957). Alternating
rows land between the two (axial 6.42 vs 6.50 predicted for UD; bending is not a plain average
once coupling is present, [45,0]s measured 1.065 vs 1.475 naive). The twist column of the tool
is meaningless where the quad twist is ~0 (UD, isotropic): a one-diagonal triangulation has its
own chirality (isotropic tri strip twists 2.1 deg at 45 m deflection), not a #34 effect.

**E2: on the production blade the defect is real but touches very little.** 9277 nodes: 9824
quads, **43 triangles = 0.44 % of the elements, 0.22 % of the area** (0.17 % inside 0.3-0.9 R),
spread evenly along the span, not only root/tip/TE. Their effective axis 1 is **~88 deg off the
span in every bin, 100 % of them above 45 deg** (the Z fallback of `rotate_constitutive_to_local`
fires only when X is exactly normal, which never holds once the blade is twisted). Where they sit
(sets with a triangle): LE/TE panels (~25), shear webs 14, TE reinforcement 2, spar cap 2
(one HP, one LP).

## Reading, and what is NOT shown

- Shown: the quad/triangle fibre-axis contract differs by ~90 deg (E1), and every blade triangle
  is on the wrong side of it (E2).
- Not shown: that this explains the 2.86x. Area share is 0.22 %; a single transverse-axis
  element in a 12x-stiffer UD cap is a local soft spot whose effect on torsion is not
  estimated by area. Plausible mechanism for an effect that is not small: the two spar-cap
  triangles and the 14 web triangles. Plausible reason it is small: 0.22 %.
- Not tested: #38 (drilling penalty on composite triangles). The isotropic control agrees to
  1 %, which says nothing about the composite penalty (`e_equiv = tr(A)/3h`).

## Handoff for the analysis session

Resume: `mem_search "shell-twist-mixed-mesh"`, then this file. Env: worktree
`/scratch/leahk/eduardo.donestevez/fem-shell-worktrees/rotor-frame`, `source scripts/aeroenv.sh`,
`PYTHONPATH=$PWD/src`, never `python -I`; the Zhou-load tool needs `.sources` (symlink to the main
checkout's, git-ignored). Logs: `logs/check/mixed_mesh_e1.log`, `mixed_mesh_e2.log`,
`zhou_loads_reverse_t7e.log` (ignored by git; rerun the tools, they take ~1-3 min each).

Next measurement, decided but NOT run (needs the maintainer's OK because it is a model change in
a sensitivity copy, not a fix): **E3, the decisive one.** Run `tools/diagnose_zhou_loads_reverse.py`
twice on the same mesh and load: (a) as is (shell +5.955 deg), (b) with the 43 triangles' composite
properties pre-rotated in Python by minus the angle `rotate_constitutive_to_local` applies, so
their effective axis 1 is the span (the #34 contract). No Rust change, no parameter moved.
Pre-register before running: the tolerance for "explains it" is that (b) moves the tip twist by
more than the 0.42 band width used in T7d (|delta| > 0.42 x 3.6 deg = 1.5 deg); anything smaller
closes mixed-mesh #34 as a non-cause of the 2.86x. If it moves less, the next candidates are #38
and section distortion (reading (i)). Do NOT propose the #34 fix as a remedy for the twist gap
unless E3 shows an effect; #34 is a correctness defect on its own and can be reported on its own.
Possible outward-facing step (maintainer's decision): comment the E1/E2 numbers on issue #34.
