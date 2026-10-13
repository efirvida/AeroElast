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

- [ ] T1 Tool written, predictions printed before solves, committed
- [ ] T2 E1 run, log in `logs/check/mixed_mesh_e1.log`
- [ ] T3 E2 run (blade mesh), log in `logs/check/mixed_mesh_e2.log`
- [ ] T4 Handoff for the analysis session (below)
