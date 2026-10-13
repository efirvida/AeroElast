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

## Closed before handing off (2026-10-13)

- Issue #34 commented with the E1/E2/E3 numbers (https://github.com/efirvida/AeroElast/issues/34#issuecomment-6145882312).
- Mesh-size guard: `scripts/check_mesh_size.py`, called by `scripts/check.sh`; shrink-only registry
  `scripts/mesh-size-debt.txt`. Found 45 files with `element_size >= 0.3`, **21 of them at 1.0 m**
  (validation tests and diagnostics), which the first 0.5 m inventory had missed. Debt, same
  deferral as the 0.5 m anchors (`docs/blade_mesh_size_rule.md`).

## Line A pre-registration (2026-10-13, committed BEFORE the tool is written or run)

**Question.** Does our MITC4/MITC3 shell reproduce an independent shell code on the exact
Zhou-load case where it turns 2.82x the beam? This separates "the shell implementation or its
discretisation is off" from "the shell model gives this and the beam cannot".

**Case.** Production blade, `element_size = 0.25` m (`docs/blade_mesh_size_rule.md`), root
clamped, the SAME nodal force vector `f` that `tools/diagnose_zhou_loads_reverse.py` builds with
`ForceProjector` from Zhou's loads (translations only, no nodal moments), applied to:
- ours: `PyMeshAssembler` MITC4/MITC3, `ZHOU_TRIANGLE_AXIS=asis` (production);
- reference: CalculiX `ccx`, quadratic mesh, `S8R` / `S6`, deck from `write_ccx_mesh(...,
  shell_element_type="S8R", load_field=...)`, the same composite layup export the blade parity
  test uses (`tests/test_blade_ccx_parity.py`).

**Primary metric.** Tip-ring section rotation `omega` (deg), computed by the SAME function
(`test_blade_rated_twist._ring_kinematics`, as the Zhou tool does) on the corner-node displacement
field of each code; CalculiX corner nodes matched to ours by coordinates (tolerance 1e-4 m).

**Reference and tolerance.** Reference = CalculiX S8R. Gate: `|omega_ours - omega_ccx| /
|omega_ccx| <= 5 %` (the suite's rule). Secondary, same 5 % gate: total work `W = 1/2 f.u` over
the free DOFs (global, node-count independent; AGENTS.md rule 3). Reported, not gated: tip
distortion share (`distortion/|omega|`) and `omega` at three further rings (50, 75, 90 % span).

**Precondition (checked before reading any result, not a tuning knob).** The resultant of the
`*CLOAD` block equals the resultant of `f` to 1e-6 relative; the deck's clamped node set equals
ours; one run per code, nothing re-run with a different setting.

**Reading, fixed in advance.**
- Both gates inside 5 %: the shell implementation and its discretisation agree with an independent
  shell code on this load. The 2.82x is then a property of the plate/shell MODEL (shared by both
  codes) against the beam's section properties. Reading (i) gains support but is NOT proven:
  S8R shares the CLT/thin-wall assumptions, so a 3D solid or VABS/BECAS judge is still the only
  closure of "the shell is physically right". Say that in the close.
- Either gate outside 5 %: the shell differs from the independent shell code on the case that
  matters. Report the number, name the cause if one is found, propose the fix and wait for OK.
- Not shown by this experiment in any outcome: that Zhou's +3.6 deg is matched, or that the beam
  is wrong.

## Line A result (2026-10-13, 0.25 m, one measurement, nothing tuned)

Job 11614310 (`sequana_cpu_dev`, 5 min 44 s), log `logs/check/zhou_ccx_s8r_11614310.out`; tool
`tools/zhou_ccx_compare.py` via `ZHOU_CCX_DIR` in `diagnose_zhou_loads_reverse.py`. Earlier
attempts (11614305, 11614307) died before any result was read, on two harness bugs of mine (a list
minus a float; an FRD coordinate match tolerance of 1e-4 m against a 5-digit FRD, ~5e-4 m
resolution at z ~ 100 m). Precondition held: `*CLOAD` resultant equals `f` (rel diff 0.0).

| metric | ours (MITC4/3) | CalculiX S8R/S6 | diff | gate |
|---|---|---|---|---|
| tip ring `omega` | +5.8550 deg | +5.9421 deg | -1.46 % | 5 % PASS |
| work `1/2 f.u` | 1.510489e6 | 1.509266e6 | +0.08 % | 5 % PASS |
| `omega` 50 % span | 1.4389 | 1.4446 | -0.40 % | reported |
| `omega` 75 % span | 3.6048 | 3.6140 | -0.25 % | reported |
| `omega` 90 % span | 5.7896 | 5.8567 | -1.15 % | reported |
| tip `distortion/|omega|` | 0.636 | 0.562 | - | reported |

**Verdict (by the reading fixed in advance).** The shell implementation and its discretisation
agree with an independent shell code on the very load where it turns 2.82x the beam. The 2.82x is
therefore a property of the shell MODEL (shared by both codes) against the beam properties, not a
MITC4/MITC3 defect. Together with E3 (#34 negligible) the mixed-mesh hypothesis has no support.

**Not shown.** (1) That the shell is physically right: S8R shares the CLT/thin-wall assumptions, so
a 3D solid or VABS/BECAS judge is still the only closure of reading (i). (2) That Zhou's +3.6 deg
is matched (shell 1.63x, beam 0.58x of it). (3) That the beam is wrong.

**Observation that weakens reading (i) and is NOT yet a finding.** The tip ring is mostly
distortion (`distortion/|omega|` 0.64), but at 75 % span the ring is nearly rigid (0.019) and the
shell already rotates 3.60 deg there. If the beam is well below 3.6 deg at 75 % span, the gap
exists where distortion is ~0 and distortion cannot be its explanation. The beam's `theta(r)` is
already computed in the Zhou tool but only the tip value is printed. Candidate next measurement
(needs pre-registration and an OK): print the beam twist at the same 50/75/90 % stations and
compare it with the shell's `omega` there.

## Line B pre-registration (2026-10-13, committed BEFORE the tool is written or run)

**What #38 is, from the code (no run yet).** Quads (MITC4+/D) use the variational drill-membrane
strain of Ko-Bathe-Zhang 2025; the composite triangle (`mitc3.rs:320,560`) uses a plain penalty
`k = 0.15 * E_eq * t^2` on `gamma = (dv/dx - du/dy)/2 - theta_z`, with `E_eq = tr(A)/(3h)` and
`drilling_scale` hard-wired to 1.0 for composites (`assembler.rs:1270`). Units: `E t^2` is Pa.m^2,
so the energy is N.m^2, not N.m; the penalty is one length unit off a stiffness per area.

**Change to the plan in "Next session".** The planned probe (flat patch, rigid rotation about the
normal, spurious energy vs h) is zero BY ALGEBRA: `u=-phi y, v=phi x, theta_z=phi` gives
`gamma=0`. It cannot discriminate, so it is kept only as a numerical check of that algebra (B0),
and the discriminating measurement is B1 below, which asks whether the penalty changes a
PHYSICAL response where drilling couples to bending (the folds of a closed section).

**B0 (element, numeric check of the algebra).** Single flat composite triangle (code 33), rigid
rotation about its normal. Metric: `u^T K u / (u_s^T K u_s)` with `u_s` a uniform stretch of
strain 1e-3, at thickness h and 2h. Gate: ratio <= 1e-10 at both thicknesses.

**B1 (the discriminating measurement).** Self-equilibrated torsion of the closed rectangular
`[45,-45]s` tube of `tests/validation/parity/test_thin_walled_tube_torsion.py` (same mesh, load,
rigid-mode removal, `theta_fit` rate), in three meshes: all quads (44), all triangles (each quad
split in two, 33), and mixed (alternate z-rows triangle/quad). The triangle drilling stiffness is
scaled through the property's `e_equiv` by s = 0.1, 1, 10, 100 (the composite triangle reads it only
for `k_drill`; a Rust change is NOT needed). Reference: the hand-written Bredt rate
`T/GJ` of that file (`_bredt_laminate`, independent CLT).
- Gate 1 (correctness at production setting, s=1): `theta_fit` rate / Bredt within 5 % for the
  triangle and the mixed mesh (the suite's rule). The quad run is reported beside it (known 0.998).
- Gate 2 (is the penalty invisible?): relative change of the triangle ratio between s=0.1 and
  s=10 (two decades). <= 1 % means the dimensional inconsistency has no consequence on this
  response; > 5 % means #38 changes a physical result and has a named mechanism.
- Reported, not gated: the rate at each s for each mesh, and the s=100 point.

**Reading, fixed in advance.**
- Gate 1 and Gate 2 inside: #38 is a dimensional defect without measurable consequence on torsion
  of a closed section. Say so, and that the blade weight is bounded by the triangle area share
  (0.12 %), which is NOT a blade measurement.
- Gate 1 fails, or Gate 2 > 5 %: #38 changes a physical result; name the mechanism, propose the
  fix (e.g. a penalty scaled by `G t`, not `E t^2`) and STOP for OK; no Rust edit without it.
- Not shown in any outcome: that #38 explains the blade's 2.82x (the blade has 88 triangles and
  its drilling scale cannot be varied for composites without a Rust change); torsion is one
  response only.

## Line B result (2026-10-13, one run, nothing tuned)

Tool `tools/diagnose_triangle_drilling.py`, log `logs/check/triangle_drilling.log` (seconds, login
node). Control before B1: the `e_equiv` knob scales the triangle's `K[theta_z,theta_z]` by exactly
10 at s=10, so the sensitivity run is not void.

| check | result | gate |
|---|---|---|
| B0 rigid rotation about the normal, energy / stretch energy | 4.0e-17 (h), 8.1e-17 (2h) | <= 1e-10 PASS |
| B1 tube torsion rate / Bredt, quad, s=0.1 / 1 / 10 / 100 | 0.99417 / 0.99417 / 0.99409 / 0.99370 | reported |
| B1 same, all triangles | 0.99375 / 0.99300 / 0.99267 / 0.99262 | s=1 within 5 %: PASS |
| B1 same, mixed (alternate z-rows) | 0.99326 / 0.99335 / 0.99316 / 0.99248 | s=1 within 5 %: PASS |
| Gate 2: triangle ratio change s=0.1 -> 10 | **0.109 %** | <= 1 % = invisible |

**Verdict (by the reading fixed in advance).** #38 is a dimensional defect (`E t^2` penalty, energy
in N.m^2) WITHOUT a measurable consequence on closed-section torsion: two decades of drilling
stiffness move the triangle result by 0.11 %, and triangle, mixed and quad meshes all sit at
0.993-0.994 of Bredt. The plan's rigid-rotation probe was zero by algebra (B0 confirms it to
1e-17), so it never could have shown the defect.

**Not shown.** (1) Not a blade measurement: the blade's composite `drilling_scale` is hard-wired
to 1.0, so varying it needs a Rust change (not made). The blade has 88 triangles = 0.12 % of the
area, so even a visible sensitivity would be weighted by that. (2) One response only (torsion of a
convex rectangular tube, flat walls, 90 deg corners); a loaded curved/tapered section or a
web-flange junction under bending is not covered. (3) That the penalty is *right*: only that it is
*insensitive* where tested. Correcting the units stays a cleanup, not a fix for the 2.82x.

## Where the 2.82x stands after Lines A and B

Ruled out with measurements: #34 triangle fibre axis (E3, -0.07 %), #38 drilling penalty (B,
0.11 % on torsion, 0.12 % of blade area), and the MITC4/3 implementation and discretisation (A,
-1.46 % vs CalculiX S8R on the very load). What remains: the shell MODEL (shared by both codes)
versus the beam's section properties, i.e. reading (i) not proven, and the load path (applying
the force at the nodes vs the beam's shear-centre arm).
Next candidates, each needs its own pre-registration and an OK:
1. Print the beam twist at 50/75/90 % span beside the shell's `omega` (Line A observation: at
   75 % the ring is rigid, distortion share 0.019, and the shell already rotates 3.60 deg; if the
   beam is far below there, distortion cannot be the whole explanation). ~3 min, no CCX.
2. External judge for "the shell is physically right": 3D solid segment or VABS/BECAS on an
   outboard section (z ~ 95-117), compared on the chord rotation.

## Next session (decided 2026-10-13, after Lines A and B): ONE task

**Task: where along the span does the shell/beam gap open, and does distortion explain it?**
Print the beam twist next to the shell ring rotation at the same stations, in the Zhou tool,
no CCX, no new physics. ~3 min on `sequana_cpu_dev` (or the login node if it fits).

Why this and not the 3D judge first: it is cheap and it can already kill reading (i). Line A
showed that at 75 % span the ring is rigid (`distortion/|omega|` 0.019) and the shell already
rotates 3.60 deg. If the beam is far below 3.6 deg at 75 %, the gap exists where distortion is ~0
and distortion cannot be its explanation. If the beam matches there and the gap only opens near the
tip, reading (i) gains real support.

**Where to hook.** `tools/diagnose_zhou_loads_reverse.py`: the beam `theta` array and `r_b`
(`hub + frac * (R - hub)`) are computed after `# the beam, on the same completed load`; only
`theta[-1]` is printed. The shell side is `t._ring_kinematics(coords, u, ring)["omega"]`, rings at
50/75/90 % already built in the `ZHOU_CCX_DIR` block (`_ring_at`). Reuse, do not rewrite.

**Things to settle BEFORE running (these are the traps):**
1. Station mapping: the shell `z` is measured from the blade root, the beam `r_b` from the rotor
   axis (includes `hub`). Interpolate `theta` at `r = hub + z`, not at `z`. A wrong offset moves the
   comparison by a hub radius and would look like a physical result.
2. Sign: the tool prints the beam as `-phi` (its `m` is nose-up positive); compare in the
   same convention as the shell (nose-down positive). Print both raw and converted once.
3. The beam only twists where `ok` is true (stations with `gkt` above the floor); report which
   stations were skipped, because a skipped station is a hole, not a zero.

**Pre-register (draft, confirm and commit before the run):**
- Metric: `theta_beam(z) / omega_shell(z)` at 50, 75, 90 % and the tip.
- Reference: the shell ring rotation on rings with distortion share < 0.05 (rigid), because there
  the two models are measuring the same kinematic quantity.
- Tolerance: "agree" = within 5 % (the suite's rule). "The gap exists without distortion" =
  ratio < 0.75 at the rigid ring (the literature shell-beam spread is 14-24 %, so < 0.75 is outside
  it). In between: report, do not call it.
- Reading: ratio within 5 % at 50/75 %, gap only at 90 %/tip -> reading (i) supported; ratio
  < 0.75 at 75 % -> distortion is not the explanation, the cause is elsewhere (load path: nodal
  force vs the beam's shear-centre arm; or the beam properties themselves); do not tune anything.

**Not in this session:** 3D solid or VABS/BECAS judge (separate task, bigger, needs its own
pre-registration); the unit cleanup of `k_drill` (#38); any Rust change.

**Conventions that carry over:** blade mesh 0.15-0.25 m only; `scripts/aeroenv.sh` with
`PYTHONPATH=$PWD/src`, never `python -I`; the Zhou tool needs the git-ignored `.sources` symlink;
job script `logs/check/zhou_ccx_s8r.srm` (ignored by git, edit the env line to drop
`ZHOU_CCX_DIR` if CCX is not needed); thresholds committed before the run; nothing tuned toward
Zhou; if a named cause appears, propose the fix and wait for the OK.
