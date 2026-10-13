# Feature: is the shell's 2.82x twist (vs the beam, same load) a mixed MITC3/MITC4 effect?

Spun off from `coupled-divergence.md` T7d/T7e (2026-10-13). Worktree `rotor-frame`, branch
`fix/bem-participant-rotor-frame`. Read-only investigation: **no element code changes without
the maintainer's OK.**

## Question

Under Zhou's published loads the production shell turns +5.86 deg, the standard beam +2.08 deg
(2.82x). T7e left two readings open: (i) the shell is right and the beam cannot represent section
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

**E2: on the production blade (0.25 m) the defect is real but touches very little.** 32336 nodes:
33385 quads, **88 triangles = 0.26 % of the elements, 0.12 % of the area** (0.11 % inside 0.3-0.9
R), spread along the whole span (log `logs/check/mixed_mesh_e2_025.log`). Their effective axis 1
is **~82-89 deg off the span in every bin** (the Z fallback of `rotate_constitutive_to_local` fires
only when X is exactly normal, which never holds once the blade is twisted). Sets with a triangle:
LE/TE panels 52, shear webs 26, TE reinforcement 4, spar cap 4, leading edge 2. (An earlier pass
at 0.5 m, 43 triangles, was discarded: 0.5 m is not a converged blade mesh,
`docs/blade_mesh_size_rule.md`. The strip of E1 does not depend on the blade mesh.)

## Reading, and what is NOT shown

- Shown: the quad/triangle fibre-axis contract differs by ~90 deg (E1), and every blade triangle
  is on the wrong side of it (E2).
- Not shown: that this explains the 2.82x. Area share is 0.12 %; a single transverse-axis
  element in a 12x-stiffer UD cap is a local soft spot whose effect on torsion is not
  estimated by area. Plausible mechanism for an effect that is not small: the four spar-cap
  triangles and the 26 web triangles. Plausible reason it is small: 0.12 %.
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
twice on the same mesh and load: (a) as is (shell +5.855 deg (0.25 m; 0.5 m gave +5.955, not citable)), (b) with the 88 triangles' composite
properties pre-rotated in Python by minus the angle `rotate_constitutive_to_local` applies, so
their effective axis 1 is the span (the #34 contract). No Rust change, no parameter moved.
Pre-register before running: the tolerance for "explains it" is that (b) moves the tip twist by
more than the 0.42 band width used in T7d (|delta| > 0.42 x 3.6 deg = 1.5 deg); anything smaller
closes mixed-mesh #34 as a non-cause of the 2.82x. If it moves less, the next candidates are #38
and section distortion (reading (i)). Do NOT propose the #34 fix as a remedy for the twist gap
unless E3 shows an effect; #34 is a correctness defect on its own and can be reported on its own.
Possible outward-facing step (maintainer's decision): comment the E1/E2 numbers on issue #34.

**Mesh rule (2026-10-13):** every blade run in this feature uses 0.25 m or finer; see `docs/blade_mesh_size_rule.md`. Re-run E3 on the same mesh as the baseline (a) above.

## E3 result (2026-10-13, 0.25 m, one run per mode, nothing tuned)

Pre-registered rule (written and committed in `e498163` before running): #34 is a relevant cause
of the shell/beam gap only if moving the triangles onto the quad contract changes the Zhou-load
tip twist by more than 1.5 deg (0.42 x 3.6 deg).

| run (`ZHOU_TRIANGLE_AXIS`) | tip twist | vs production |
|---|---|---|
| `asis` (production) | +5.8550 deg | - |
| `split` (88 triangles in their own sets, laminates unchanged; control) | +5.8550 deg | 0.0000 |
| `span` (same, ply angles shifted by psi so axis 1 = projected span) | **+5.8511 deg** | **-0.0039 deg (-0.07 %)** |

- Sign and method of the shift were calibrated first on a strip turned 30 deg
  (`--calibrate`): +psi gives tri/quad 1.033 axial, 0.996 bending; -psi gives 1.647 / 0.482;
  unshifted 6.46 / 2.63. psi over the 88 blade triangles: min/median/max |psi| 20.7 / 90.1 /
  121.4 deg, so the shift is large and the null is not a no-op.
- **Verdict: #34 is NOT a cause of the 2.82x.** The effect is 0.07 % against a 1.5 deg bar
  (385x below it). Logs: `logs/check/zhou_e3_split.log`, `zhou_e3_span.log`.
- #34 stays a correctness defect on its own (E1: tri/quad up to 12x on UD), with negligible
  weight on this blade's torsion because 88 triangles are 0.12 % of the area. Reporting it does
  not depend on the twist gap.
- Limits: one metric (tip ring rotation under Zhou's loads, linear static); the in-plane
  triangle fibre axis only. #38 (drilling penalty of composite triangles) was not touched by E3.
- Still open for the 2.82x: (i) the shell is right and the beam cannot represent section
  distortion (`distortion/|omega|` = 1.67 in the rated case); (ii-bis) #38. A discriminating
  next step for (i), not run: shell vs beam on a section with an exact distortion-free answer
  (closed tube, torsion: shell 0.998 of Bredt in `test_thin_walled_tube_torsion.py`, so the
  shell torsion itself is validated where distortion is absent).

## Next session (decided 2026-10-13): two open lines, in this order

Both are about the same open question: why the production shell turns 2.82x the standard beam
under the same load (Zhou's, +5.855 vs +2.079 deg at 0.25 m). #34 is ruled out by E3.

**Line A - section distortion (reading (i): shell right, beam cannot represent it).**
- Idea: find a case where shell and beam have an independent reference AND the section distorts.
  The closed rectangular tube is distortion-free (shell 0.998 of Bredt), so it validates torsion
  but not this. Candidates to check first, in the repo: `tests/test_box_multicell_torsion.py`,
  `docs/twist_distortion_diagnosis.md` (S-8c, `distortion/|omega|`), `tools/diagnose_zhou_loads_reverse.py`
  with the section-distortion split already printed by `test_blade_rated_twist._ring_kinematics`.
- Pre-register before running: the metric (`omega` vs `theta_fit` of the tip ring, i.e. the
  distortion share), the reference (an independent 3D/shell code, CalculiX S8R, or an analytic
  open-section/warping result) and the tolerance (5 %, the suite's rule).
- Done when: either the shell matches an independent reference where the beam cannot (reading (i)
  gets evidence) or it does not (then the shell over-twists and the cause is next to find).

**Line B - #38, drilling penalty of composite triangles.**
- `mitc3.rs:320`: `k_drill = e_modulus * t^2 * 0.15 * drilling_scale`, dimensionally inconsistent;
  for composites `e_modulus = tr(A)/(3h)`.
- Cheap probe, same style as E3 and no Rust change: the Zhou-load shell with `drilling_scale`
  varied is NOT possible for composites (hard-wired 1.0), so instead use the strip of
  `tools/diagnose_mixed_mesh_consistency.py` with a flat patch under a rigid rotation about its
  normal: spurious energy should be zero and must not scale with h (the issue's own "closes when").
- Pre-register: spurious energy / membrane energy at two thicknesses (h and 2h); a ratio
  different from 1 within 5 % confirms the h-scaling. Expected weight on the blade is small for the
  same reason as #34 (triangles are 0.12 % of the area); say so if it is.

Conventions that must carry over: blade mesh 0.15-0.25 m only (`docs/blade_mesh_size_rule.md`);
bootstrap `scripts/aeroenv.sh` with `PYTHONPATH=$PWD/src`, never `python -I`; the Zhou-load tool
needs the git-ignored `.sources` symlink to the main checkout; thresholds are committed before the
run; nothing is tuned toward Zhou. Work in this worktree/branch; commits are local, no push.
