# Core torsion reliability (shell FE torsional stiffness)

Gate before continuing L3 (T6/T7/T8) and L4/L5: decide whether the shell FE (MITC4) or the twist measure is soft
in torsion on thin, twisted or tapered sections, and fix it if it is the core. Origin: `odd/tasks/l3-closure.md` T9,
Handoff 4. Branch `l3/core-torsion` from main 3a2d222.

## Specs

- S1 (user, verbatim): "recuerda que no es que yo quiera arreglar el modelo de una u otra pala, yo quero que aeroelast sea cientificamente confiable para el proposito que te di"
- S2 (user, verbatim, after the reliability summary): "entonces sabes ien lo que tienes que corregir antes de poder seguir"
- S3 (repository physics discipline): every number with a metric, an independent reference and a written tolerance;
  prediction written before each run; no pass band widened; a moved anchor is a finding.
- S4 Done means: the IEA-3.4 outboard softness (FE twist 1.2-1.5x the integral of T/GJ of its own mesh section at
  eta 0.65-0.85) is attributed to element or twist measure with closed-form evidence, and, if element, fixed with
  the element benchmarks (D-Tube, KO2017 3.5 twisted beam, thin tubes) and S-7 / S-4 re-measured.

- S5 (user, verbatim): "es que no entiendo, yo siempre lo unico que he intercambiado desde el BEM hacia la estructura son fuerzas nodales, no se como es que ahora se intercambiamomento, la idea es recibir fuerzas nodales, porq es lo mismo que va a recibir cuando cambie el solver BEM por openfoam q entrega los vectores de fuerza en cada nodo"
- S6 (user, verbatim): "si, pero eso que me cuentas es trabajo que tiene que hacer el solver bem a partir de los desplazamientos que recibe del solver estructural, la responsabilidad de adaptarlos a su problema es del solver aerodinamico, no del solver estrucutral, ademas ahora me viene a la zabeza que tambien hay que revisar q los marcos de referencia esten bien aplicados al hacer esas transformaciones"

- S7 (user, verbatim, on how to fix D1): "el solver inercial va a ser descontinuado, el solver estructural siempre va a ser corrotacional"

## Tasks

- C1 | S3 | parent | resubmitted | Baseline: full suite at d20911c (main 3a2d222 + tools/docs only) via `scripts/check_sbatch.sh full` (first try 11625395 failed: no git on the node), compare with known-red-full at a9ca498.
- C2 | S3 S4 | parent | done (D-Tube, KO2017 green) | Collect the three symptoms with their current numbers: D-Tube torsion vs Bredt, KO2017 3.5 twisted beam (MITC4 thin, 8.5 % vs 1 %), 3.4 outboard 1.3x (`tools/l3_twist_integral.py`); IEA-15 S-7 1.072. Prediction of what an element defect vs a distortion measure would show in each.
- C3 | S3 S4 | parent/worker | done e60b7ca: measure | Split element vs measure on the 3.4 at 0.65-0.85: rigid-fit residual per station; twist fitted on the box nodes (caps + webs) vs all nodes; element quality of the outboard rows (aspect, warp) and es 0.15 vs 0.10.
- C4 | S3 S4 | parent/worker | done 4809427 as prismatic extrusion: core right | Tapered closed-form tube through MITC4 at es 0.15 with the 3.4 outboard perimeter, wall G t and taper (extend `tools/diagnose_s1_gj_tube.py` / `tests/test_s1_gj_tube.py`), RED first if it fails.
- C7 | S6 | parent/explore | pending | Frame audit of the coupling transformations: structure-published displacements and received forces (frame, units), `BEMFSIParticipant._displacements_from_partner` / `_forces_to_partner` / `_frame_theta_end`, radius/twist extraction, `ForceProjector` rotor-plane frame; each step mapped to a test with an independent reference, gaps listed.
- C8 | S6 S7 | worker (after C6) | pending | D1 fix: the BEM participant always exchanges with a corotational structure that publishes R u and reads inertial forces: remove the `rotating` partner frame (a YAML declaring it fails at load with the reason), the inertial transform is always on; RED test first (participant round trip against the Rust convention at theta = 90 / 180 deg); migrate in-tree BEM YAMLs; then D3, D2, D4 with RED tests; D5 recorded for the OpenFOAM interface.
- C5 | S4 | parent | revised: twist-definition sensitivity done d20911c; user decision pending | If element: fix (test first), then element benchmarks, S-7, S-4 on both blades; if measure: fix the twist extraction and re-measure S-1 GJ / S-7.

## C2 result (main 3a2d222, login, 54 s; `logs/check/core-torsion/c2-benchmarks.log`, 14 passed)

| benchmark | measured | bound |
|---|---|---|
| D-Tube GJ vs Bredt-Batho closed cell (`test_dtube_torsion_matches_bredt`) | 1.33336e11 vs 1.34554e11, -0.91 % | 2 % |
| KO2017 3.5 twisted beam MITC4, t 0.02667 / 0.0002667, in / out of plane | 0.10 / 0.25 / 0.06 / 0.01 % | 1 % |
| Thin twisted beam vs Simo, Fox & Rifai 1989 | 0.24 % | 5 % |
| S-1 tube truth tests (rectangle, circle, oblique, S-7 slices, torque helper) | 7 passed | per test |

**Correction to L2:** the "element reds" (D-Tube, KO2017 twisted beam 8.5 %) are not red on main; the AGENTS.md
note dates from 890d4e6. The element is right in torsion on prismatic thin closed sections and on a twisted
thin strip. What none of them has: **taper** (perimeter and wall changing along the span), a multi-cell section
with webs, and the distorted outboard quads of a real blade. The 3.4 softness has to come from one of those or
from the twist measure (C3, C4).

## C3 prediction (written before the run)

On the IEA-3.4 at eta 0.65-0.85, same solve as `tools/l3_twist_integral.py` (T40 shear-flow tip torque):
- **Measure (section distortion)**: the twist fitted on the box nodes only (elements of the spar caps `*SPAR*` and
  webs `*_SW`) gives an increment within 0.95-1.05 of the mesh-Bredt integral while the all-node fit stays at
  ~1.3; the rigid-fit residual (RMS of u - (u0 + theta x d) over |theta| RMS |d|) grows outboard above 5 %.
- **Element**: the box-only ratio is also >= 1.2 and the residual stays below 2 %.
- Mixed (box 1.05-1.2): both; then C4 decides the element share.
Control IEA-15 at 0.3-0.85: box and all-node ratios within 0.05 of each other, residual < 2 %.

## C3 result (worker; jobs 11625399 3.4, 11625400 15; `logs/check/core-torsion/c3-*.log|json`; `tools/l3_twist_integral.py --c3`)

Box = element sets `*SPAR_A/B/C` (caps) and `*_SW` (webs). FE delta-theta / mesh-Bredt integral, all nodes -> box:

| blade | 0.3-0.6 | 0.6-0.85 | 0.7-0.85 |
|---|---|---|---|
| IEA-3.4 | 1.014 -> 1.083 | 1.314 -> **0.926** | 1.380 -> **0.911** |
| IEA-15 | 1.002 -> 0.971 | 0.952 -> 1.018 | 0.932 -> 1.032 |

Rigid-fit residual (all / box), 3.4: 0.7 0.13 / 0.045, 0.75 0.14 / 0.034, 0.8 0.24 / 0.020, 0.85 0.36 / 0.051
(inboard 0.46-0.69 / 0.04-0.08, small |theta| in the denominator); 15: 0.04-0.10 / 0.009-0.04 at 0.7-0.85.
Elements benign (median aspect 1.2-2.1, warp < 1 deg) except the thin TE reinf/flat strips (max aspect 31-84).
Verdict: **measure** (prediction held in kind; box 0.91-0.93 is below the 0.95-1.05 band, and the 15 control missed
its 0.05 / 2 % limits: distortion is present on both blades, larger on the 3.4). The FE box is not soft; the
section distorts in its plane, and an all-node rigid fit reads that as twist. Open: is the distortion physical
(non-prismatic: taper, pretwist, prebend; no diaphragms) or an artefact; and which twist the BEM feedback uses
(`fsi_participant.py _compute_deformed_geometry`).

## C4 (revised) prediction (written before the run)

Prismatic extrusion of the real IEA-3.4 section at eta 0.8 (a copy of the yaml, under `logs/check/core-torsion/`,
with chord, airfoil/relative thickness, pitch axis and every layer and web value held at their eta-0.8 value over
the whole span; twist, prebend and sweep zero), same mesher, es 0.15, same solve. Saint-Venant/Bredt is exact for
a prismatic tube away from the ends. Predicted, at eta 0.3-0.7 (away from root and tip):
- all-node and box ratios both within 0.97-1.03 of the mesh Bredt, residual all-node < 3 %: then the 3.4
  outboard distortion is a **non-prismatic physical effect**, element and measure are right, and the open item is
  which twist definition the aerodynamics and the beam comparison must use;
- all-node ratio >= 1.1 or residual >= 5 % on the prismatic: the distortion is produced by the section itself
  (shell model of the TE strips, load path), an artefact to fix in the model or element.

## C4 result (worker; job 11625404; `logs/check/core-torsion/prismatic/`; `tools/core_torsion_prismatic_yaml.py`)

Prismatic copy of the 3.4 at eta 0.8 (chord 1.9, pitch axis 0.308, twist 0, no prebend/sweep, DU91/DU08 blend
w 0.5576 -> t/c 0.2277, all 17 layers and both webs held at their eta-0.8 values): mesh Bredt GJ 9.561e6 at all 100
stations; FE / mesh Bredt **0.982 (all nodes) = 0.982 (box)** at every slice 0.25-0.85 (0.989 at 0.9, tip);
residual 0.15-0.45 % (all and box). **First reading holds:** element and twist measure are right on the real
multi-cell section; the 3.4 outboard distortion (C3) is a non-prismatic physical effect (taper, pretwist, prebend)
that a rigid-section beam cannot represent. The uniform 1.8 % stiffer-than-Bredt is the thick-wall / open-wall term
(T27 measured +2.0 pt on the IEA-15), not distortion. TE strips aspect 49-51 do no harm.
**Consequences.** (1) The core is reliable in torsion: no element fix (C5 as written does not apply). (2) Against
the beam references (S-1 GJ, S-7) the beam-equivalent twist is the box twist; the all-node twist carries physical
section distortion. (3) Production BEM feedback (`fsi_participant._section_rotation`: antisymmetric part of an
affine fit over the whole ring) includes the non-affine distortion; the aerodynamic twist is the chord-line
rotation. Which one the BEM must take is a modelling decision, to be made on its measured torque effect.

## C5 (revised): twist definition sensitivity; prediction (written before the run)

IEA-15, rated, one-way, same solve as `tools/l3_torque_sensitivity.py` (rigid production loads, es 0.25): per BEM
strip, elastic twist by (a) the production ring-affine rotation, (b) the chord line (LE and TE nodes of the ring),
(c) the box rigid fit; rigid-geometry BEM torque/thrust with each. Prediction: (b) and (c) within 5 % of (a) in
tip elastic twist (the 15's all-node residual is 4-10 % at 0.7-0.85), so torque within 0.15 % of (a) (twist alone
moves torque -3.2 %). Refuted if any definition moves torque > 0.3 % or tip twist > 10 %: then the definition is a
thesis-level choice and goes to the user with the numbers.

## C5 result (worker; job 11625406; `logs/check/core-torsion/c5/`; `tools/l3_torque_sensitivity.py --twist-definitions`)

IEA-15 rated one-way, rigid-geometry BEM with the deformed radius; default lines reproduce T4 (tip +3.758 deg, torque
-4.404 %, twist alone -3.175 %). (a) production `_section_rotation` (called, matches to 0.0 rad); (b) chord line,
LE / TE from `ForceProjector._section_ends` (flat-back: midpoint of the corner nodes within 1 % chord); (c) box
rigid fit (`*SPAR*`, `*_SW`, >= 6 nodes per ring).

| def | thrust vs (a) | torque = power vs (a) | twist 0.7 / 0.85 / 0.95 R / tip [deg] |
|---|---|---|---|
| a production | 0 | 0 | +0.547 / +1.850 / +4.163 / +3.758 |
| b chord line | +1.711 % | **+0.457 %** | +0.478 / +0.625 / +4.849 / +6.780 |
| c box | +0.768 % | **+0.249 %** | +0.635 / +1.684 / +1.684 / +10.058 |

Profile: at 0.25-0.8 R the definitions agree within 0.23 deg; at r 99-109 m (b) is up to 1.4 deg below (a); the
last three strips (r 113-121 m) are erratic in all three (a 4.27 / 2.82 / 3.76). **Prediction refuted** for (b)
(torque 0.46 % > 0.3 %), (c) between bands; tip twist differences 80-170 % (tip strips not credible). Scale: the
definition is worth 0.25-0.46 % of rated torque (5 % acceptance; twist alone 3.2 %), thrust 0.8-1.7 %.
Side findings, not chased: the BEM prints `NaNs at 0/50` lines in this run; the tip-strip rings give erratic
rotations in every definition. **Goes to the user** (thesis-level choice) with these numbers.

## C7 result (explore, read-only, code + YAML; consequences inferred, not measured)

Path: Rust rotor (rotating frame, X0 = inertial at theta 0) reads Force as **inertial** and applies R^T(theta_target)
(`rotor_fsi.rs:799-815`), writes Displacement as R(theta_target) u, elastic only (`:955-958`); BEM participant reads u,
un-rotates only if `partner_frame == "inertial"`, **default "rotating" = no transform** (`fsi_participant.py:275-280`,
`:557-567`, theta_end same time level as Rust); ForceProjector Np +y, Tp +x toward the LE, torque about +y positive
(`test_force_projection_tangential_sense.py:75-82`); pitch/cone/tilt only on the BEM side, all 0 in-tree.
Suspected defects, by consequence for torque (parent spot-checked D1 in the code):
- **D1** only 6 YAMLs declare `partner_frame`; the ~38 BEM cases paired with the Rust rotor solver (e.g.
  `IEA15MW/frontiersin_2025_yaw/fluid_yaw_*.yaml`) run with the silent default: the structure receives
  R^T(theta) F_rot (torque ~cos theta, reversed at 180 deg) and the BEM sees u_x/u_z mixed (radius, twist). No guard.
- **D2** power/CP/CQ sign stale since T7c (Tp +x, tau_aero > 0): `rotor_physics.rs:836-845`, `rotor.py:1846-1853`,
  `:2246-2250` negate it.
- **D3** Python step callback rotates the already-rotating X0 again (`rotor.py:2164-2171`, `:2185`): reported torques
  modulated ~cos theta, disagree with the Rust tau that drives omega.
- **D4** theta integrated separately on both sides (preCICE exchanges omega, not theta): drift on restart, without an
  omega mesh, and in computed-omega mode at iteration 0.
- **D5** written displacement is R u, not the total R(X+u) - X: right for the BEM, wrong for a partner (OpenFOAM)
  expecting total displacement; the inertial solver's total mode pairs wrongly with the BEM `inertial` mode.
Gaps (no independent check): theta agreement across participants; partner_frame declared for Rust pairs; Python vs
Rust torque and the CP/power sign; BEM integrated torque vs torque of the nodal forces; CCBlade azimuth sense
(yaw cases); no precurve to CCBlade (prebend ignored); rotation-axis defaults disagree (`config.py:344` [1,0,0],
`rotor.py:223` (0,0,1), `fsi_participant.py:285` [0,1,0]).

## C6 result (worker; jobs 11625411 base, 11625421 final; `logs/check/core-torsion/tip/`)

**NaN: defect, fixed.** Deformed root strip (r 3.97, chord 5.2) came back at Rhub - 8e-6 m (second-order span
shortening > first-order stretch); CCBlade hub loss flips sign, residual NaN, strip 0 load silently zeroed (Np 0 vs
115 N/m rigid). Fix: `_rebuild_bem_solver` holds r at Rhub. RED -> GREEN `tests/validation/rotor/test_bem_nan_inputs.py`
(2 tests); `test_bem_fsi_deformed_geometry.py` + `test_bem_participant_rotor_frame.py` 61 passed; quick OK. Deformed
torque -0.0008 %, T4 numbers unchanged, 0 NaN lines. **Tip strips: benign.** Right rings (tip row z 117.0, 16 nodes;
114.61, 26; 112.22, 30); the production rotation steps 0.1-0.34 deg where the ring node count changes (least-squares
fit on a non-rigid section); smoothing the last 3 strips moves torque -0.05..-0.13 %. Open: a node-count-independent
ring fit (arc-length weighting) would shift T4; not started. C5 table unchanged by the fix.
Note: the src fix landed while the full suite 11625408 was running (after collection); C1 measures 4809427+docs.

## C8 D1 result (worker; not a measurement of past runs)

`BEMFSIParticipant` always un-rotates read displacements and rotates written forces with R(theta_end) (the Rust
corotational convention); `partner_frame: inertial` accepted and default; `rotating` rejected at load with the reason.
RED 6 failed / 4 passed -> GREEN (`tests/validation/rotor/test_bem_partner_frame.py`: Rust R^T recovers F_rot and u at
90 / 180 deg); participant + NaN tests 65 passed; deformed-geometry, torque-sense-vs-OpenFAST, wall-flow 18 passed;
quick OK. Updated `tests/validation/bem/test_bem_participant_rotor_frame.py` (3 places relied on the rotating default).
**Inventory: 67 BEM fluid YAMLs; 65 omitted `partner_frame` and so ran untransformed against the corotational
solver** (29 dirs Rotor-only, e.g. `IEA15MW/frontiersin_2025_yaw`, `convergence_b1_b2`, `parked_v50`, `ab_test`,
`smoke_fix/*`, `exp_*`, `zeta_sweep`; 37 dirs mixed with the inertial solver, e.g. `IEA15MW/ch6_prepared_cases`).
**Every coupled BEM run made from those cases before this commit is suspect** (torque ~cos theta, BEM geometry
from mixed u_x/u_z) and needs re-running before it is evidence.

## C8 D3 / D2 result (worker; commits ab46009 D3, next D2)

**D3 confirmed, fixed:** the step callback rotated the rotating-frame coordinates again for the aero, inertial, gravity
and total torques; now `_window_torques` uses them as Rust does. RED 3 failed / 1 passed (aero 198 % off at 180 deg,
inertial 1063 % at 90 deg) -> GREEN vs the closed form sum (X+u) x F . axis (no Rust binding exposes compute_torque).
**D2 confirmed, fixed:** CP, CQ, power_aero, power_total now +tau_aero omega (Tp drives the rotor about +axis since T7c);
RED -> GREEN Python `test_rotor_power_sign.py` 2/2 and Rust `driving_torque_reports_positive_cp_and_cq` (was cp -0.289);
`cargo test -p aeroelast-solvers --lib rotor_physics` 48, `-p aeroelast-core` 170; maturin release rebuilt; rotor tests 94
passed, inertial/wiring 78 passed 7 skipped; quick OK. Computed-omega mode unaffected (omega used tau, never power).
Open: inertial solver power sign left inconsistent (to be discontinued, S7); `shaft_torque` sign not reviewed. Past
rotor CSVs: CP, CQ, power_aero and callback torques at theta != 0 are suspect; the Rust-driven omega was not affected.
Rust test needs RUSTFLAGS="-L /scratch/app/openmpi/4.1.4_gnu/lib -C link-arg=-lmpi" to link on this cluster.

## Log

- L1 (2026-10-26) User: S1 verbatim (during the PR of l3/twist-p1, merged as #68 / 3a2d222).
- L2 (2026-10-26) Parent reliability summary: bending and frequencies within 2.4 % (IEA-15) and 4 % (IEA-3.4) of the
  BeamDyn 6x6; torsion not reliable (IEA-15 S-7 7 % soft; 3.4 outboard FE 1.3x its own section; KO2017 twisted beam
  and D-Tube element reds); mode shapes unmeasured (T8). Recommendation: core torsion first, then T8.
- L3 (2026-10-26) User: S2 verbatim.
- L4 (2026-10-26) User: "continua en esta sesion". Parent: C2 (benchmarks green), C3 (measure), C4 (prismatic: core right).
- L5 (2026-10-26) User chose, on the twist definition for the BEM: "Primero arreglar punta y NaN, después decidir" (over keep (a) and declare, chord line, box).
- L6 (2026-10-26) User: S5. Parent: BEM -> structure is nodal forces only (`ForceProjector`, the pitching moment as a force couple); the twist definition is the structure -> BEM geometric reduction, needed only by a strip model; OpenFOAM takes nodal displacements.
- L7 (2026-10-26) User: S6. Parent: the reduction already lives in the BEM participant (`fsi_participant.py:297`, reads `Displacement` :628, writes `Force` :688); the structure only publishes displacements; C7 added (frame audit).
- L8 (2026-10-26) User: S7 (custom answer to the D1 fix options). Parent: single structural partner -> BEM always inertial, `rotating` mode removed; C8.
