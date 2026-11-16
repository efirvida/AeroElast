# L5 one-way: BEM loads -> shell structure, no feedback (thesis plan step 4)

Step 4 of odd/tasks/thesis-validation-plan.md. Opened 2026-11-15 after L4 closed (plan L5). Closes when every row
below is trusted against its independent reference or declared with cause and torque/power sensitivity.

## Specs

- S1 (plan step 4, verbatim): "BEM loads -> structure, no feedback | OpenFAST BeamDyn under the same loads (S-2 type), S-5 rated run | 5 % on deflection and twist; S-6 NaN closed"
- S2 (plan S3, verbatim): "every compared number has an independent reference and a written tolerance (5 % relative ceiling); a red that stays red is reported with its cause and, where it can move the thesis quantity, with its measured effect on rotor torque and power (sensitivity), never hidden by widening."
- S3 (handed over by L4, odd/tasks/l4-bem-rigid-closure.md:60,144, verbatim): "production tilt/precone must match the deck (-6 / -4) or the difference be declared (+2.44 % torque)".
- S4 (handed over by L3, plan L4 / AGENTS.md): the production-twist estimator (S8, `BEMFSIParticipant._section_rotation`, equal nodal weights) "stays deferred to L5".

## Tasks

- T1 | S1 S2 | parent | pending | BeamDyn under the S-5 AeroDyn loads: deflection OoP/IP and tip twist vs our shell (replaces the ElastoDyn IP apples-to-oranges bound)
- T2 | S1 | parent | done (L12: S-6 3 passed on HEAD, no NaN) | S-6 NaN: reproduce on HEAD (mid-surface vs default OML skin), root cause, close or declare
- T3 | S1 S3 | parent | pending | our rigid BEM loads (deck tilt -6 / precone -4, T7 model) through the production load mapping onto the shell, the same loads onto BeamDyn; load conservation check; 5 % on deflection and twist
- T4 | S4 | parent | pending | twist metric for the one-way deliverable: equal-weight vs skin-weighted estimator on T1/T3, decide or carry to two-way
- T5 | S1 S2 | parent | pending | store rows + verdicts, thesis plan T4 closed or declared
- T6 | S1 S3 | worker | in progress (smoke 11741409 on dev, L15) | S-5 rated run with the production code: `aeroelast-bem-fsi` one-way (`deformed_twist: false`, `deformed_radius: false`), gmsh N20 AR 1, deck tilt -6 / precone -4, rated 10.59 m/s 7.56 rpm, ~100 s; its settled mean tip deflection and twist must match the T3 static solve (same loads, same K_eff) within 1 %
- T7 | S3 | worker (with T6) | done (L15: L5 cases explicit; legacy cases untouched by user decision) | tilt/precone explicit in every L5/L5-two-way case file: deck -6/-4 for OpenFAST comparisons, 0/0 only for Zhou cases, labelled; no case relies on the 0/0 default
- T9 | S1 S2 | worker | pending | flap-edge cross-coupling arbiter (L9): per-station bending compliance of shell and BeamDyn under the same loads, deck 6x6 cross terms and principal angle, and a third reference independent of both
- T10 | S1 S2 | parent | pending | twist-sense arbitration (L10): is the shell pretwist mirrored or a consistent mirror world; then the correct mesh->BeamDyn map and a rerun of the T1/T3 twist and edge rows
- T11 | S1 S2 | worker + parent | in progress (fix committed, L14; cache rebuild + L3 re-measure pending) | fix the windIO-v2 pretwist sign (geometry.py:538), RED-first; mesh cache invalidation; then re-measure the orientation-dependent L3 rows (before/after, same reference) and rerun T1/T3
- T8 | S2 | worker | done (ce50b2e, L6) | torque sensitivity harness: a measured spanwise twist difference (shell - BeamDyn) -> delta torque / power / thrust at rated with the production BEM (reuse tools/l3_torque_sensitivity.py), ready for any red of T1/T3

## Log

- L1 (2026-11-15) User: "si" (to: "¿Abro L5 one-way ya, en paralelo con T7, o esperamos a que T7 confirme? Recomiendo abrirlo:
  el primer trabajo es armar el documento de L5 y resolver tilt/precone, que no dependen de T7."). Parent: step 4 opened.
  L4 T7 blade-test check (SLURM 11739091) still running; moved from sequana_cpu_long to sequana_cpu.
- L2 (2026-11-15) Parent, inventory (explorer): S-5 = tests/test_iea15mw_s5_oneway.py + tools/run_s5_oneway.py, AeroDyn
  azimuth-mean blade-1 loads (tests/IEA15MW/reference/s5_ad_blade1_loads.csv) on the shell at 7.56 rpm (K+K_G+K_SP, gravity);
  ElastoDyn ref OoP 16.1033 m, IP -0.8103 m; es 0.25 mid-surface OoP 16.1630 (+0.37 %), IP 2.0000 m (bound 2.5x, not a
  comparison). S-6 = tests/test_iea15mw_s6_oneway_dynamic.py + tools/run_s6_oneway_dynamic.py: mid-surface mean -0.86 %,
  std +8.88 % (bound 25 %); divergence notes (~1e40-1e74 m) are on the default OML skin, no recorded root cause
  (tools/diagnose_s6_divergence.py exists). No "our BEM -> shell, no feedback" runner and no BeamDyn-same-loads row exist;
  no store group for L5. Production BEM tilt/precone default 0/0 (fsi_participant.py:468-469; gen_v5_s5_fluid.yaml:42-43).
  S8 is consumed only in the feedback loop (fsi_participant.py:823,878): a two-way concern unless one-way twist is measured with it.
- L3 (2026-11-15) Parent, T1/T3 design (one tool, `tools/l5_oneway_bd_parity.py`, `--loads ad|bem`):
  * mesh: reference gmsh chord-scaled N20 AR 1, mid-surface (`build_v2_blade(0.25, mid, gmsh, (0.05,0.25), 20)`), not the S-5
    numad mesh (AGENTS rule 9);
  * "same loads" = the shell's own nodal force vector f, binned in 0.5 m span strips into resultant force and moment about the
    BeamDyn reference axis (key-point curve mapped to the mesh frame) and fed to `beamdyn_driver` as point loads: both models see
    identical station resultants whatever the shell distribution (AD case: equal-per-node split, no Mm; BEM case: production
    ForceProjector with arm and Cm);
  * rotation 7.56 rpm on both, hub radius 3.97 m on both (shell K + K_G + K_SP with rotation centre (0,0,-3.97); BD GlbPos 3.97),
    no precone on the structure, gravity off on both (gravity is S-2's anchor, not the coupling's);
  * aero-induced response only: shell u = K_eff^-1 f (linear about the rotating prestress); BD delta = run(load) - run(no load),
    at load scale 1 (`nl`, geometrically exact) and 0.01 x100 (`lin`): lin is the like-for-like reference, nl/lin is the
    geometric-nonlinearity factor, reported;
  * metrics: deflection at the reference axis (ring affine fit evaluated at the axis point) flap/edge along span and at tip; elastic
    twist along span (BD phi_z = 4 atan(c_z/4); shell production equal-weight estimator and skin-weighted); W = 1/2 f.u per model;
    load conservation (total F and root M of the rows vs f).
  Prior art: twist_p1 C2b (P_shell -> BD at Omega 0) and the L3 unit-force twist scatter (shell/beam 0.28-1.5 by load case,
  Engram discovery/unit-force-twist-shell-vs-beamdyn-iea15): twist may not close at 5 %; a red goes through S2 (cause + torque
  sensitivity), not a wider band.
- L4 (2026-11-15) User: "eso va a demorar un tiempo podemos ir planificando otros analisis?"; chose all four proposals: A (pre-register
  T1/T3 criteria), B (rated one-way with the production code -> T6), C (tilt/precone explicit in the cases -> T7), D (torque
  sensitivity harness -> T8). Jobs in flight: 11739134 (T1 ad), 11739135 (T3 bem), 11739136 (S-5+S-6 tests on HEAD 2c7a47f, T2),
  11739091 (L4 T7).
- L5 (2026-11-15) Parent, PRE-REGISTRATION of T1/T3 acceptance, written before any number of 11739134/11739135 was seen:
  * reference = BeamDyn `lin` (load x0.01 rescaled), the like-for-like counterpart of the linear shell; `nl` is reported, not gated;
  * deflection at the reference axis, flap and edge: |shell - BD_lin| / |BD_lin| <= 5 % at 0.7R, 0.9R and tip (R = 117 m arc z);
    0.3R-0.6R reported only;
  * elastic twist, gated estimator = skin-weighted ring fit (the blade guards' estimator); the production equal-weight value is
    reported beside it (T4 input): <= 5 % relative at 0.9R and tip; at 0.5R-0.8R the band is 5 % of the BD_lin tip twist (absolute),
    because the inboard twist tends to zero and a relative band would be meaningless there;
  * twist sense: R1 of tools/twist_p1.py (both nose-up positive about +z); an opposite sign is a finding reported as measured,
    never flipped to agree;
  * geometric nonlinearity: nl/lin tip flap reported; if it exceeds 5 % it is declared as a linear-shell one-way limitation with its
    torque sensitivity (T8), not used to move the gate;
  * load identity: the binned rows reproduce the nodal total force and root moment (relative difference < 1e-9, exact by
    construction); a failure invalidates the run, not the gate;
  * any red: cause + measured torque/power/thrust sensitivity (T8), per S2.
- L6 (2026-11-15) Parent:
  * T8 done (ce50b2e): tools/l5_twist_torque_sensitivity.py, 7 unit tests; sense d theta_BEM = -omega_mesh
    (fsi_participant.py:824,863); demo 0/-0.5/-1.0 deg at z 60/90/117 -> torque/power -0.22 %, thrust -2.48 %.
  * T6 code (a8130d0): `deformed_twist`/`deformed_radius` were INERT in production (declared in BEMConfig, read nowhere;
    coupled-divergence.md T6c): now read in `_compute_deformed_geometry`, false = reference twist/radius, defaults bit-identical;
    the exp_feedback_off cases ran two-way until today. gmsh sizing keys (chord_divisions, span_aspect, h_range, min_band_*) are
    now forwarded by the CLI and the runner (silently dropped before). check.sh quick OK. Case cases/l5_oneway_rated (a8130d0),
    tilt -6 / precone -4 explicit; full run NOT launched (smoke first to size it: mesh + modal analysis > 3 min on login node).
  * T7 (read-only): every existing fluid YAML is 0/0 by default (32 ch6_prepared_cases incl. S-5/S-6, convergence_b1_b2,
    parked_v50 -> proposed -6/-4; Frontiers/Zhou yaw cases and their diagnostics -> 0/0 written explicitly). None edited.
  * Load identity of the parity tool measured on the login node: ad 236 rows, total F (76 735, 827 942, 0) N exact, root M rel
    diff 1e-15; bem (deck tilt/precone) total F (104 656, 849 884, 0) N, 3e-16.
  * INCIDENT: every SLURM job of the user NODE_FAILs at start and is requeued (11739091 ~340 restarts since 11-14, on many nodes;
    no job completed since 2026-11-14). 11739134/11739135 cancelled; T1/T3 run on the login node instead (MUMPS solve).
- L7 (2026-11-15) Parent, T1/T3 RESULT (HEAD a8cb874, login node, 5 min both; logs/check/l5_oneway/{ad,bem}.json, gmsh N20 AR 1
  38 500 nodes, 236 point loads), shell vs BeamDyn lin, judged against the L5 pre-registration:
  | quantity | ad (S-5 AeroDyn loads) | bem (production BEM, deck tilt/precone) | gate | verdict |
  |---|---|---|---|---|
  | flap 0.7R / 0.9R / tip | -2.9 / -3.3 / -2.2 % | -2.7 / -3.2 / -2.1 % | 5 % | PASS |
  | edge 0.7R / 0.9R / tip | -41 / -37 / -27 % (tip 1.170 vs 1.611 m) | -45 / -43 / -34 % (0.921 vs 1.394 m) | 5 % | RED |
  | twist tip (skin) | -4.51 vs -4.64 deg, -2.8 % | -2.24 vs -2.55 deg, -12.0 % | 5 % | ad PASS, bem RED |
  | twist 0.9R (z 105.3, skin) | -3.93 vs -3.50, +12 % | -3.54 vs -1.80, +97 % | 5 % | RED |
  | twist 0.5-0.8R (band 5 % of tip) | up to 0.81 deg off (band 0.23) | up to 0.70 deg off, shell sign + inboard (band 0.13) | abs | RED |
  BD geometric nonlinearity (nl/lin, reported): tip flap -1.3 % (ad) / -1.6 % (bem); tip twist +30 % (ad -4.64 -> -6.02) and
  +25 % (bem -2.55 -> -3.20): > 5 %, so declared as a linear-shell one-way limitation.
  Torque sensitivity (T8, production rated BEM, twist only, % of rigid; shell - BD lin): ad skin torque +0.07 % thrust +1.64 %;
  ad eq +0.66 % / +2.47 %; bem skin -1.60 % / -1.14 %; bem eq -0.89 % / -0.02 %. BD nl - lin twist (bem): torque -0.51 %,
  thrust -4.23 %. Edge red: cause not yet known (S-4 1st edge reads the shell SOFTER, -5.6 % frequency, here it deflects LESS);
  next diagnostic: the same loads parked (Omega 0, shell K only vs BD static) to separate the rotating terms from the section.
- L8 (2026-11-15) Parent, edge-red attribution (ad loads, parked, lin; `--rpm 0`, `--keep x|y`; logs/check/l5_oneway/ad-parked*.json):
  * parked reproduces the rotating edge red (-29 % tip, -34..-43 % inboard): NOT the rotating terms (K_G, K_SP, hub radius);
  * edgewise (mesh x) load only: edge 0.308/0.309 m at 0.7R, 0.634/0.625 at 0.8R, 1.144/1.091 at tip (shell/BD, <= +4.9 %):
    the direct edgewise stiffness AGREES (consistent with L3: EI_edge CLT mean +0.4 %);
  * flapwise (mesh y) load only: edge (cross) 0.075/0.357, 0.064/0.525, 0.048/0.591 m; and reciprocally under the edge load the
    flap cross term is 0.026 vs 0.081 m at tip: the red is the FLAP-EDGE BENDING CROSS-COUPLING, 5-12x larger in the BeamDyn deck
    6x6 than in the shell. Flap direct agrees (-2..-3 %);
  * twist split: edge load +1.34 vs +0.95 deg tip (+40 %), flap load -6.63 vs -6.32 (+5 %); they partly cancel in the full load.
  Open: which side is right about the cross-coupling needs an arbiter independent of both (the deck 6x6 cross terms in the
  principal/section frame vs the shell's recovered section compliance at the same stations; candidates: the S-1 sectional tools,
  or the deck's PreComp/SONATA provenance). Not chased in this session.
- L9 (2026-11-15) User: "claro que 1, vamos a diagnositicar primero para entender lo que pasa" (option 1: section arbiter for the
  flap-edge cross-coupling before any declaration). Parent: T9 opened. Plan, three independent views at z = 0.5-0.95R:
  (a) model-level: from the existing --keep x / --keep y parked runs, the bending curvatures kappa_flap, kappa_edge per station of
      both models (shell rings: rotations about x and y from the affine fit of u_z; BD: RDx/RDy) against the statically exact
      internal moments M_flap, M_edge of the applied rows -> the 2x2 bending compliance per station, cross term included;
  (b) deck: the BD 6x6 bending block per station, its frame (BeamDyn convention: which frame the 6x6 is expressed in, how the
      key-point twist enters), principal-axis angle vs the chord/twist angle;
  (c) third reference: an in-repo section analysis of the windIO layup (CLT/PreComp-type EI_flap, EI_edge, EI_cross, principal
      angle) independent of both the shell FE and the deck.
  The model whose 2x2 compliance matches (c) within 5 % is right; the other is the finding.
- L10 (2026-11-15) Parent, T9 RESULT (worker; tools/l5_section_coupling.py, 13 unit tests; logs/check/l5_oneway/t9/summary.md):
  * the two models ORIENT THE SECTION OPPOSITELY about +z: BeamDyn negates the file's kp twist on read (BeamDyn_IO.f90:898) and
    rotates the section-frame 6x6 at run time (BeamDyn.f90:2692-2710, Stif = RR0 Stif0 RR0^T); the shell mesh chord sits at
    +twist (+15.2 deg at z 6, +1.70 at 58.5, -2.1 at 105.3). Harness validated: a uniform BD blade matches the deck compliance at
    Rz(-tw) to 5e-4 and is off 0.8-8x at Rz(+tw);
  * principal-axis angle: the deck in the MESH twist sense matches the independent yaml thin-wall reference (c) within 0.82 deg at
    all six stations; the deck in BeamDyn's sense does not (up to 4.6 deg). A 1D bending prediction with the deck in BeamDyn's sense
    reproduces the L8 BeamDyn cross terms (0.70 / 0.082 vs measured 0.591 / 0.081 m); the reference (c) gives 0.042 vs the shell's
    0.048 m. So the L8 edge red is a mirrored section orientation between the two models, not a stiffness error;
  * which sense is physical is NOT closed: AeroDyn BlTwist is "positive to feather, leading edge upwind"; the mesh has LE at +x and
    +y downwind, which the worker reads as physical chord angle -tw, i.e. the shell pretwist mirrored. The repo itself disagrees:
    tests/test_iea15mw_s5_oneway.py:38 "BlTwist, negated - the shell is built mirrored in the twist" vs tools/run_s5_oneway.py:240
    "follow the WindIO twist sign (rotorspin fix)". A mirrored blade spinning the mirrored way is a consistent mirror world for the
    aero loads, but then MESH_TO_BD (a proper rotation) is the wrong map to BeamDyn (needs det -1) and bend-twist coupling signs
    of the two codes are not comparable as computed; the twist rows of L7 are suspect until this closes;
  * secondary: (c) itself is not a 5 % reference (mid-surface cut vs OML: C_ff +1..+29 %); BD's own solve with one order-10
    element does not reproduce its deck compliance outboard (C_ff up to 73 %, as twist_p1 C1); shell curvature fit unresolved at
    0.9R/0.95R (window spread 30-160 %).
  Next (T10, new session): arbitrate the twist SENSE with one sign chain from physical definitions (rotor axis and spin sense,
  upwind, LE, AeroDyn BlTwist, windIO twist, mesh chord) -> decide mirror world vs mirrored-pretwist defect; see
  odd/tasks/rated-twist-sign-convention.md and tools/diagnose_sign_chain.py.
- L11 (2026-11-15) User: "vamos a ahcer todo lo que se pueda hacer sin el hpc". Parent: continue on the login node: T10 now,
  S-6 test (T2) on the login node (background, 50 min cap); T6 needs hours of coupled run, stays for SLURM. Existing case YAMLs
  are not edited (tilt/precone, T7) until the user decides: they are the provenance of past campaign results.
- L12 (2026-11-15) Parent:
  * T2: `tests/test_iea15mw_s6_oneway_dynamic.py` on HEAD f7f3c35+ (login node, numad mid-surface 0.25): 3 passed in 300 s, no
    NaN. The recorded divergence was on the default OML skin (L2); the mid-surface S-6 anchor runs clean.
  * T10 RESULT (worker; tools/l5_twist_sense_chain.py, 13 unit tests; logs/check/l5_oneway/t10/summary.md): verdict (a), the
    windIO-v2 shell blade carries a MIRRORED PRETWIST; not a mirror world (mesh frame proper: LE +x, suction +y downwind, span +z,
    det +1, prebend upwind, spin +y drives LE-first, production BEM torque drives the rotor). Evidence: the mesh chord sits at
    psi = +BlTwist (0.37 deg); with the BEM inflow the +twist chord sees an AoA off alpha_BEM by up to 39.6 deg, the -twist chord
    within 0.046 deg. Defect line: src/aeroelast/models/blade/numad/objects/geometry.py:538 `angle = deg2rad(twist) + correction`
    then Rz(+angle) (:539-541); windIO twist is positive to feather, i.e. Rz(-twist) in the windIO frame (parent re-derived: LE at
    -y must move to -x upwind -> negative rotation about +z). The v1 NuMAD path places it at -twist (geometry.py:577-579).
    Falsifier: a negated-twist scratch deck (no code change) gives tip edge-under-flap 0.831 m (stock 0.048, BD 0.591) and
    flap-under-edge 0.108 (stock 0.026, BD 0.081): sign and magnitude of the cross-coupling now follow BeamDyn; direct flap/edge
    barely move (-4.1 % / +4.7 %); cross terms still 33-41 % above BD (not closed).
  * Blast radius: every windIO-v2 mesh (gmsh_skin / windio2_mesh / build_v2_blade): the reference gmsh N20 mesh, the v2 S-1, S-2,
    S-4, S-7, V-02, V-03, V-05 tests, the L3 closure rows measured on it, L7 T1/T3, cases/l5_oneway_rated;
    tests/software/mesh/test_windio2_frames.py encodes the + sign. Stale comments: tools/run_s5_oneway.py:239-244 (+theta frame),
    tests/test_iea15mw_s5_oneway.py:38. The correction term's sign under the flipped twist is unverified (0 for a straight axis).
  * Consequence for the strict order: L3 was closed on this mesh; a mesher fix re-opens the L3 rows that depend on the section
    orientation (bend-twist, flap-edge coupling, modal 1E/1F) -> user decision.
- L13 (2026-11-15) User: "q hacemos?" -> parent recommended: fix the mesher sign, re-measure the orientation-dependent L3 rows
  on the login node, rerun T1/T3, SLURM-only items later. User: "arreglalo aqui". Parent: T11 opened.
- L14 (2026-11-16) Parent: T11 fix verified. geometry.py `angle = correction - twist` (chord at psi = -twist), v2_blade
  CACHE_VERSION 4 with geometry.py in the digest, apply_ad_loads frame at psi = -BlTwist, docs/windio2-input.md updated.
  tests/software/mesh/test_windio2_frames.py: GREEN 17 passed; against HEAD geometry.py RED 8 failed / 9 passed
  (logs/check/l5_oneway/t11/frames_{green,red}.log). check.sh quick OK (75 passed, 1 xfailed); the one ruff F541 in the
  touched files is pre-existing. Every pre-fix v2/gmsh cache entry is now unreachable; L3 orientation rows re-open until
  re-measured. Status (user question): T2 passed (L12), T8 harness committed (ce50b2e), both await the T5 store rows.
- L15 (2026-11-16) User: "que podemos ir haciendo mientras esperamos?" -> parent proposed D1-D4 triage, doc status, T6 launch,
  T7 with new case files instead of editing legacy ones. User: "si". Parent:
  * T7 closed: the only L5 case (cases/l5_oneway_rated) already writes tilt -6 / precone -4; legacy cases stay as written
    (provenance); step-5 two-way cases are created with explicit, labelled values.
  * T6: smoke (0.5 s) submitted, SLURM 11741409 on sequana_cpu_dev, on the fixed mesh (GmshBladeMesh -> geometry.py); full run
    after the smoke sizes it. run-t6.sh: stale "uncommitted edits" banner replaced by HEAD + dirty count.
  * D1-D4 orientation triage (l3-outboard-section-stiffness.md:9-12): D1 Bredt GJ and D4 (S-7 metric) are section-intrinsic,
    invariant under an in-plane rotation of the section -> unaffected; D2 (LC3 edge station ux under a global-x load) depends on
    the section angle -> re-measured by the v2 S-2 rerun (SLURM 11741381); D3 (tip FE local GJ) can move through the twist-rate
    sign, weakly -> re-check with T5 of l3-outboard if S-7 moves.
  * T11 T1/T3 RERUN (HEAD 0c04c36, login node, gmsh N20 AR 1, same tool and BeamDyn lin x0.01 reference as L7;
    logs/check/l5_oneway/t11/{ad,bem,ad-parked-x,ad-parked-y}.json; z = 35.1 58.5 70.2 81.9 93.6 105.3 111.2 117.0):
    | row | before (L7/L8) | after | gate | verdict |
    |---|---|---|---|---|
    | ad flap | -1.5..-3.4 % | -6.4 -5.4 -5.6 -5.6 -5.1 -4.4 -4.1 -3.6 % | 5 % | RED inboard z<=93.6 (was PASS) |
    | ad edge | -27..-41 % | -1.6 +8.6 +11.3 +13.3 +15.2 +16.9 +17.3 +18.8 % | 5 % | RED, sign flipped, halved |
    | ad twist tip (skin) | -2.8 % | +10.5 % (-5.63 vs -5.09 deg) | 5 % | RED (was PASS) |
    | bem flap | -1.3..-3.3 % | -6.5..-3.6 % (same shape as ad) | 5 % | RED inboard |
    | bem edge | -34..-45 % | -1.8 +8.0 .. +17.0 % | 5 % | RED |
    | bem twist tip / 0.9R | -12 % / +97 % | -11.2 % / +140 % | 5 % | RED |
    | parked flap load: edge-under-flap (cross) | tip 0.048 vs 0.591 m (-92 %) | 0.834 vs 0.591 m (+41 %) | - | sign now BeamDyn's, magnitude +41 % |
    | parked flap load: flap (direct) | -0.3..-3.0 % | -4.1..-6.8 % | 5 % | RED inboard |
    | parked edge load: flap-under-edge (cross) | -68..-90 % | +6..+26 % | - | sign right, +25 % tip |
    | parked edge load: edge (direct) | -4.1..+4.9 % | +2.7..+4.5 % | 5 % | PASS |
    Load totals: ad F_x 76 735 -> 97 486 N (the AD Fn/Ft frame now at -BlTwist, so the AD loads themselves were mis-rotated
    before); bem identical (104 656, 849 884 N). Reading: the mirror was the dominant edge error (edge -35 % -> +15 %), and
    the remaining edge excess tracks the cross-term excess (+25..+41 %), i.e. the shell's flap-edge coupling is stronger
    than the deck's; T9 (cross-coupling arbiter) is the open question now. New: the direct flap compliance moved from
    -2..-3 % to -4..-7 % (shell stiffer), red inboard; the L3 S-2/S-4 reruns (11741380/1) will say whether this is
    blade-intrinsic or load-path.
