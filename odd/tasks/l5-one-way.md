# L5 one-way: BEM loads -> shell structure, no feedback (thesis plan step 4)

Step 4 of odd/tasks/thesis-validation-plan.md. Opened 2026-11-15 after L4 closed (plan L5). Closes when every row
below is trusted against its independent reference or declared with cause and torque/power sensitivity.

## Specs

- S1 (plan step 4, verbatim): "BEM loads -> structure, no feedback | OpenFAST BeamDyn under the same loads (S-2 type), S-5 rated run | 5 % on deflection and twist; S-6 NaN closed"
- S2 (plan S3, verbatim): "every compared number has an independent reference and a written tolerance (5 % relative ceiling); a red that stays red is reported with its cause and, where it can move the thesis quantity, with its measured effect on rotor torque and power (sensitivity), never hidden by widening."
- S3 (handed over by L4, odd/tasks/l4-bem-rigid-closure.md:60,144, verbatim): "production tilt/precone must match the deck (-6 / -4) or the difference be declared (+2.44 % torque)".
- S4 (handed over by L3, plan L4 / AGENTS.md): the production-twist estimator (S8, `BEMFSIParticipant._section_rotation`, equal nodal weights) "stays deferred to L5".

- S5 (user L32, verbatim): "Claro que B, si vamos por A el codigo sera cada vez mas inmanejable" -- B = the Rust RotorFsi K_G rebuild uses static-solve centrifugal membrane stresses (K u = f_cf on the free DOFs, stress recovery, K_G(u) on deformed coords), the same physics as the Python `assemble_geometric_stiffness`; the local sigma = rho w^2 r l_char estimate is removed, not kept as a fallback.

## Tasks

- T1 | S1 S2 | parent | pending | BeamDyn under the S-5 AeroDyn loads: deflection OoP/IP and tip twist vs our shell (replaces the ElastoDyn IP apples-to-oranges bound)
- T2 | S1 | parent | done (L12: S-6 3 passed on HEAD, no NaN) | S-6 NaN: reproduce on HEAD (mid-surface vs default OML skin), root cause, close or declare
- T3 | S1 S3 | parent | rerun on the production load (L42): flap pass, edge and twist red | our rigid BEM loads (deck tilt -6 / precone -4, T7 model) through the production load mapping onto the shell, the same loads onto BeamDyn; load conservation check; 5 % on deflection and twist
- T4 | S4 | parent | pending | twist metric for the one-way deliverable: equal-weight vs skin-weighted estimator on T1/T3, decide or carry to two-way
- T5 | S1 S2 | parent | pending | store rows + verdicts, thesis plan T4 closed or declared
- T6 | S1 S3 | worker | done (L41: dyn mean vs static_t6 +0.30 % flap, -0.34 % edge, +0.61 % twist on 048bf2f) | S-5 rated run with the production code: `aeroelast-bem-fsi` one-way (`deformed_twist: false`, `deformed_radius: false`), gmsh N20 AR 1, deck tilt -6 / precone -4, rated 10.59 m/s 7.56 rpm, ~100 s; its settled mean tip deflection and twist must match the T3 static solve (same loads, same K_eff) within 1 %
- T7 | S3 | worker (with T6) | done (L15: L5 cases explicit; legacy cases untouched by user decision) | tilt/precone explicit in every L5/L5-two-way case file: deck -6/-4 for OpenFAST comparisons, 0/0 only for Zhou cases, labelled; no case relies on the 0/0 default
- T9 | S1 S2 | worker | partial (L38, 3c293a5) | flap-edge cross-coupling arbiter (L9): per-station bending compliance of shell and BeamDyn under the same loads, deck 6x6 cross terms and principal angle, and a third reference independent of both
- T10 | S1 S2 | parent | pending | twist-sense arbitration (L10): is the shell pretwist mirrored or a consistent mirror world; then the correct mesh->BeamDyn map and a rerun of the T1/T3 twist and edge rows
- T11 | S1 S2 | worker + parent | in progress (fix committed, L14; cache rebuild + L3 re-measure pending) | fix the windIO-v2 pretwist sign (geometry.py:538), RED-first; mesh cache invalidation; then re-measure the orientation-dependent L3 rows (before/after, same reference) and rerun T1/T3
- T12 | S1 S2 | parent | done (L17) | measured orientation audit (tools/l5_orientation_audit.py): physical frame, mesh rings vs yaml vs BeamDyn deck (LE, suction side, chord angle, prebend, sweep, chord), physics-only AoA, response-angle reading of the T1/T3 differences
- T13 | S1 S2 | worker | in progress (L18-L19) | pre-fix baseline of the v2 L3 rows (worktree ../fem-shell-prefix at ad1aecf, same commands, gmsh N20), for the before/after table
- T14 | S1 | worker | done for mass/EI (L23); EA/GJ extraction on gmsh open | S-1 v2 sectional extractor on the gmsh mesh: m(r) 72 %, GJ 20 %, EA 1.99e10 % are orientation-free -> diagnose the extraction (login node), real EI_edge/GJ/mass values or a declared defect
- T15 | S1 | worker | done (L24, 93beb2d) | T6 cost: 11.7 s per preCICE window with a constant K_eff; profile the smoke on the login node, find where the time goes (refactorization, Python callback, BEM, I/O), report before any code change
- T16 | S2 | worker | done (L22) | torque/power/thrust sensitivity of the current T1/T3 reds (t11 jsons: twist shell - BeamDyn spanwise) with tools/l5_twist_torque_sensitivity.py, ad and bem
- T17 | S5 S1 | parent | implemented (L34), T6 passes (L41); full suite pending (L41) | Rust K_G rebuild with static-solve prestress (option B): RED test (Rust K_G == Python K_G at X0), port, maturin, quick + rotor tests, full suite on SLURM, rerun T6
- T18 | S1 S2 | worker | done (L37, 77d7257); prebend-drop fixed (L39) | aero torsional load: T6 time-mean F_AERO vs T3 `_compute_forces(zeros)` give the same flap but twist -3.20 vs -2.25 deg on the same K_eff (L31); per-strip force and torsional moment about the reference axis, cause (Cm/arm projection, azimuth mean with tilt, BEM config fluid.yaml vs bem_config), fix or declare
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
- L16 (2026-11-16) User: "lo de la pala espejada lo deduciste a partir de un comentario viejo escrito, que tal si en vez de
  leer mides todo en la malla y en la geometria para asegurarse que este todas las orientaciones y los prebend bien
  orientados. Ademas de esto no entiendo porque las diferencias tan grandes". Parent: T12 opened, every orientation measured
  from coordinates (mesh rings, yaml, BeamDyn key points), AoA from kinematics only (no BEM), and the T1/T3 edge differences
  read as a rotation of the deflection vector (edge is ~1/12 of flap, so 1 deg of response direction = ~20 % of edge).
- L17 (2026-11-16) Parent, T12 RESULT (HEAD 0c04c36 + tool, login node 100 s; logs/check/l5_oneway/t12/audit.json; 32 airfoil
  rings z 10-115 of 50, gmsh N20 AR 1 mid):
  * frame measured: rotor axis +y, Omega > 0 -> blade velocity +x at +z; production BEM thrust +849 884 N along +y -> wind +y,
    upwind -y.
  * mesh vs physics: LE advances with the blade velocity on every ring; suction (cambered) side faces downwind on every ring.
  * pretwist: |psi_mesh + twist_yaml| <= 0.35 deg (max at z 26, transition; <= 0.12 deg from z 45 out), |psi - twist| up to
    16 deg -> chord at -twist; BeamDyn deck initial_twist == yaml twist (3e-5 deg) and the deck section maps to the same psi.
  * prebend: tip ring reference point mesh y -3.988 m, yaml x -3.995 m, deck kp_xr -> mesh y -3.995 m (all upwind);
    max |mesh - yaml| 0.012 m, |mesh - deck| 0.012 m along span; sweep |x| <= 0.004 m (yaml 0); chord within 0.47 %.
  * physics-only AoA (V(1-a), a = 1/3, Omega r, no BEM): measured chord 5.8-8.2 deg along the span (near-uniform, as a designed
    rotor); mirrored chord 2.6-24.6 deg. Independent of any twist-sign convention: the fixed mesh is the physical blade.
  * response direction (T11 jsons): shell deflection vector rotated toward +edge vs BeamDyn by 0.5-1.3 deg (ad, bem), 0.1-1.1
    deg parked pure flap; the pretwist placement (<= 0.12 deg outboard) cannot produce it -> the remaining edge difference is
    stiffness coupling (principal axes / 6x6 cross terms, T9), not orientation.
- L18 (2026-11-16) Parent, T11 v2 L3 rerun on the fixed mesh (SLURM 11741381, gmsh N20 AR 1 via V2_BLADE_MESHER=gmsh
  V2_BLADE_CHORD_DIVISIONS=20, 17 min; logs/check/l5_oneway/t11/static-11741381.out): 16 passed, 7 failed, 1 xfailed.
  * PASS: S-2 all (LC1-LC4 dominant component, LC3 edge tip ux, interior station LC1-LC4 -> D2 station now inside 5 %),
    S-1 EI_flap, V-03 (mass, gravity), V-05.
  * FAIL: S-1 EI_edge 9.0 %, S-7 torsion ratio 8.30 % (pre-fix gmsh N20 +5.90 %, blade_mesh_size_rule.md:31-33),
    V-02 1F +7.2 %, 1E +6.8 % (5 % bound, vs NREL/TP-5000-75698 Table 5-2).
  * S-1 m(r) 72 %, GJ 20 %, EA 1.99e10 %: orientation-independent quantities -> extraction artifact of the S-1 sectional
    slicer on the gmsh mesh, not physics; no gmsh pre-fix value recorded.
  * None of these v2 ids are in scripts/known-red-full.txt; no pre-fix run of this file set on gmsh N20 exists except S-7,
    so the before/after needs a same-command pre-fix run.
- L19 (2026-11-16) Parent:
  * S-4 v2 on the fixed mesh (SLURM 11741380, gmsh N20, 36 min): 4 passed. Parked 1F +2.19 %, 1E -1.29 %, 2F -0.02 %; rated
    7.56 rpm 1F +1.72 %, 1E -1.38 %, 2F -0.13 % vs OpenFAST v5 BeamDyn p20 (5 % bound). Pre-fix gmsh N20 1E -1.55 %
    (blade_mesh_size_rule.md:30). V-02's 1F +7.2 % / 1E +6.8 % are against the NREL report table (ElastoDyn-type values; the
    same S-4 output puts ElastoDyn 4.7-7.4 % above BeamDyn), not against BeamDyn.
  * T6 smoke (SLURM 11741409, 14 min, both participants exit 0): mesh 163 s, 50 windows (0.5 s) + solver setup 583 s ->
    <= 11.7 s per window; the full 100 s (10 000 windows at dt 0.01) is ~20-32 h, beyond run-t6.sh's 4 h -> user decision on
    duration / restart chaining.
  * Pre-fix baseline worker: first cache job 11741507 failed (worktree src/_aeroelast stub shadows the built extension);
    worker steered to a shim PYTHONPATH.
- L20 (2026-11-16) User: "creo que podemos usar las colas que permiten jobs de hasta 96h". Parent: T6 full 100 s submitted
  as-is (dt 0.01, 10 000 windows, est. 20-32 h): `sbatch -p sequana_cpu_long --time=96:00:00 logs/check/l5_oneway/run-t6.sh
  full` -> SLURM 11741547 (PENDING Priority). Output $SCRATCH/l5_oneway/full-11741547/.
- L21 (2026-11-16) User: "y no hay nada que podamos hacer mientras tanto 30h+tiempo de espera en cola" -> parent listed T9,
  S-1 extractor, T6 cost, T8 on the reds, step-5 prep. User: "arranca todo lo que puedas arrancar en el nodo login que nuestra
  cola en el HPC ya tiene 45 trabajos y no estoy seguro siquiera que pueda manejar tantos jobs". Parent: login-only workers
  for T14, T15, T16 and the step-5 two-way preparation (thesis-validation-plan.md T5, no numbers claimed); T9 (CalculiX,
  SLURM) waits for the queue; no new SLURM submissions.
- L22 (2026-11-16) Worker T16 (login; T8 harness on the t11 fixed-mesh tables, rated, tilt -6/precone -4; logs/check/l5_oneway/t16/
  summary.{json,md}, run_t16.py). Delta in % of rigid, power == torque (fixed rpm):
  | case | ad thrust | ad torque | bem thrust | bem torque |
  |---|---|---|---|---|
  | shell skin - BD lin | -2.34 | -2.57 | -3.85 | -2.71 |
  | shell eq - BD lin | -1.67 | -1.96 | -2.80 | -1.94 |
  | shell skin - BD nl (x1) | +0.77 | -0.41 | +0.71 | -0.80 |
  | shell skin vs rigid | -16.95 | -8.42 | -9.62 | -3.85 |
  | BD lin vs rigid | -14.60 | -5.84 | -5.78 | -1.14 |
  | BD nl vs rigid | -17.71 | -8.01 | -10.33 | -3.05 |
  The linear twist difference moves torque 2-2.7 % (> 1 %: cannot be declared negligible); against nonlinear BeamDyn the shell
  is within 0.4-0.8 %. Flap/edge deflection: no harness path; first-order radial foreshortening estimate dR ~ -0.11 m ->
  dP <~ -0.2 % (indicative only). Open question for T4/T9: why a LINEAR shell tracks NONLINEAR BeamDyn twist (coincidence
  of compensating errors vs the shell capturing a coupling BD-linear misses) -> needs the T9 arbiter before any claim.
- L23 (2026-11-16) Worker T14 (login): S-1 v2 extractor defect on gmsh. tests/support/v2_sectional.py ring_ids_for_v2_mesh
  (np.unique(eta), assumes the structured windio2 shared eta row) produced 27 621 one-to-64-node "rings" on gmsh -> 815
  arbitrary stations, 99.6 % centroid fallback, mass / tiny dz. Fix: unstructured meshes (element without exactly 2 eta
  values) use uniform 0.25 m z strips (GMSH_STRIP_M); windio2 unchanged; no tolerance or src touched. Independent strip-mass
  cross-check (rho t A per element, no extractor): total 68 191 kg gmsh / 68 193 kg windio2 vs BeamDyn 66 912 kg (+1.9 %);
  strip m / BD 0.94-1.09 z 1-102, tip 1.63 on both meshers. S-1 gmsh N20 after the fix (mid-span mean |rel err|, 5 % bound;
  logs/check/l5_oneway/t14/s1_gmsh_fixed.log): m 71.97 -> 4.35 % PASS, EI_flap 4.14 % PASS, EI_edge 9.00 -> 8.83 % FAIL
  (genuine, T9 family), GJ 20.16 -> 12.10 % FAIL, EA(CLT) 1.99e10 -> 14.48 % FAIL (not verified as physics: strip-width /
  edge effect suspected, windio2 EA/GJ pipeline not compared). The T13 pre-fix baseline runs the OLD extractor, so its S-1
  gmsh values are artifacts; only S-2/S-4/S-7/V-02 are comparable before/after.
  Parent check: tests/support/test_v2_sectional.py 11 passed (windio2 default); under the gmsh env 5 passed / 6 errors both
  with and without the fix (fixtures pass a collapse ratio the gmsh mesher rejects: pre-existing, not a regression). RED/GREEN
  evidence for the fix is the S-1 gmsh m(r) row (71.97 % -> 4.35 %); no synthetic unit test for the unstructured branch yet.
- L24 (2026-11-16) Worker T15 (login profile; logs/check/l5_oneway/t15/summary.md): ~8 s per window = ~6-7 s Python stress
  recovery in the step callback every window (_compute_stress_fields 6.9 s/call; stress_output_interval default 1,
  rotor.py:577, :2132-2141, against its own docstring :2503) + ~2 s BEM wait (4 IQN-ILS sub-iterations/window) + newmark 50 ms;
  cached KSP works; refactorization only after a K_G rebuild (~+15 s per event, 2 in 50 windows). The YAML key was silently
  dropped (SolverConfig/_build_model_config whitelist): worker 2 plumbed it, RED 4 failed -> GREEN 4 passed
  (tests/software/test_stress_output_interval_option.py), quick OK; L5 cases set 50 (93beb2d). Estimated full T6 ~5.6 h + setup.
  T6 11741547 (15 min in, every-window stress) cancelled, resubmitted as 11741636. Fix-2 (fewer sub-iterations) not taken:
  changes coupling semantics.
  T13: worker stalled (never applied the shim steer); 11741508/9 cancelled; parent re-ran provenance (pyshim -> prefix src,
  `angle = deg2rad(twist) + correction`, CACHE_VERSION 3, _aeroelast from the venv) and runs the static set + S-4 on the
  LOGIN node (user L21), separate cache dir; logs/check/l5_oneway/t13/{static,s4}-login.out.
- L25 (2026-11-16) User: "si, todo lo que puedas correr en el login tienes permiso para hacerlo". Parent: login runs started:
  step-5 OpenFAST flex/rigid 60 s (logs/check/l5_twoway/openfast-login.out; reference data only, step 4 still open), T14b
  worker (S-1 GJ/EA on gmsh: windio2 comparison, strip-width scan, independent EA/GJ), T13 baseline continues on login.
- L26 (2026-11-16) Parent:
  * Login node limit: `ulimit -t` = 1800 CPU-s per process. OpenFAST flex 60 s was killed at t = 22.75 s (rc 137, threaded
    BLAS); rerun single-threaded (OMP/OPENBLAS = 1). Rigid twin 60 s finished normally. Long pytest runs on login must be split
    per test node and single-threaded (S-4 rerun per node).
  * T13 BEFORE/AFTER, static set, same command, gmsh N20 AR 1, BeamDyn p60 / NREL report references (before = ad1aecf via the
    pyshim, logs/check/l5_oneway/t13/static-login.out; after = 0c04c36, t11/static-11741381.out). Both 16 passed / 7 failed.
    | row | before (mirrored) | after (fixed) | bound |
    |---|---|---|---|
    | S-2 LC2 tip flap uy | +6.32 % FAIL | +4.75 % | 5 % |
    | S-2 LC4 distributed uy | -2.93 % | -4.97 % | 5 % |
    | S-2 LC1 gravity station | -1.48 % | -4.57 % | 5 % |
    | S-2 LC2 / LC4 station | -0.50 / -2.52 % | -1.81 / -4.59 % | 5 % |
    | S-2 LC3 tip edge ux | +12.14 % | +11.78 % | (diagnostic tip; station test passes) |
    | S-7 torsion ratio | 5.9 % FAIL | 8.30 % FAIL | 5 % |
    | V-02 1F / 1E (NREL report) | 8.02 / 6.76 % | 7.22 / 6.84 % | 5 % |
    | S-1 EI_edge (old extractor both) | 8.72 % | 9.00 % | 5 % |
    Reading: with the physical pretwist the shell is 2-3 % stiffer in global flap (every flap row moves toward -5 %),
    consistent with T1/T3 (flap -2..-3 % -> -4..-6.5 %): bending about the global axis now engages the stiff edge axis at
    the correct angles. The global flap rows now sit at the 5 % edge; S-7 moved away (5.9 -> 8.3 %).
- L27 (2026-11-16) Worker T14b (login): S-1 EA/GJ on gmsh.
  * EA(CLT): the v1 CLT loop (tools/run_s1_sectional.py:287-300, all-nodes-in-strip) aliases on gmsh (0.67-1.33x at z 72-92;
    moves with strip width; 0.125 m leaves empty stations). Fix: plane-cut EA (`ea_clt_plane_cut`, tests/support/v2_sectional.py,
    unstructured meshes only). Independent cut (EA = sum A11_eff L over the z = const cut, own CLT): +6.9 % on BOTH meshers vs
    deck K33 -> S-1 EA(CLT) gmsh 14.48 -> 7.00 % (windio2 7.12 %): a genuine +7 % red.
  * GJ: windio2 12.13 % = gmsh 12.10 %; closed-form multi-cell Bredt on the same cut (tools/l3_bredt_mesh.py graph_bredt,
    G t = A66) already -9.3..-9.9 % vs K66 (D1 family), the FE shell adds ~-2.5 %: physics, not extraction.
  * EI_edge 8.8-9.0 % on both meshers (genuine). tests/support/test_v2_sectional.py 11 passed after the change (parent).
  * S-1 file on gmsh: 3 failed (EI_edge, GJ, EA), 2 passed (m, EI_flap).
  T13 S-4 BEFORE (login, per node, single thread; t13/s4-login.out) vs AFTER (11741380), BeamDyn p20, 5 % bound:
    parked 1F +1.31 -> +2.19 %, 1E -1.37 -> -1.29 %, 2F -0.79 -> -0.02 %; rated 1F +1.02 -> +1.72 %, 1E -1.55 -> -1.38 %,
    2F -0.87 -> -0.13 %. All pass on both; the fix stiffens 1F by ~0.7-0.9 % (same direction as the S-2 flap rows).
- L28 (2026-11-16) User: "metele" (prepare the T9 arbiter on the login node). Parent: worker prepares the CalculiX solid
  segment arbiter (z ~60/80/95: flap moment, edge moment, torque -> section compliance incl. flap-edge cross term, principal
  angle, torsion) for shell, BeamDyn deck 6x6 and solid on the same cut; login smoke only if < 1800 CPU-s per process.
- L29 (2026-11-16) Worker T9 prep (tools/l5_t9_arbiter.py, 13 unit tests RED->GREEN; reuses tools/l3_solid_section_ccx.py and
  l5_section_coupling.deck_at). One fit (fit_section_kinematics) -> 3x3 compliance [k_flap, k_edge, twist rate] per [M_flap,
  M_edge, T] in the mesh frame for shell (6-chord segment of the N20 mesh), deck 6x6 and layered solid (ccx). SMOKE (coarse,
  not evidence; logs/check/l5_oneway/t9b/smoke/table.md):
  | z | model | psi deg | C_fe/C_ff | GJ/deck | EI_flap/deck | EI_edge/deck |
  |---|---|---|---|---|---|---|
  | 60 | shell | -2.56 (deck -1.49) | -0.029 (deck -0.017) | 0.90 | 1.04 | 1.05 |
  | 80 | shell | 0.52 (deck 1.23) | -0.003 (deck +0.010) | 0.82 | 1.00 | 0.94 |
  | 95 | shell | -0.68 (deck 2.91) | -0.032 (deck +0.026) | 0.84 | 0.99 | 0.86 |
  | 95 | solid A | 2.52 | +0.022 | 0.80 | 0.99 | 0.94 |
  First sign (z 95 only, coarse): the solid sides with the deck on the principal angle / cross term and with the shell on GJ
  (deck TE filler, D1). Risks: shell psi moves +-1 deg with segment length and +-2 deg along span (window decision needed);
  shell GJ depends on fit weights; the deck carries a TE filler the shell/solid lack (fill job); level B settings unvalidated.
  Login runtimes: solid A 4.5 min wall (3 ccx x ~230 CPU-s). sbatch prepared, not submitted (logs/check/l5_oneway/t9b/).
- L30 (2026-11-16) Parent, T9 solid runs on the LOGIN node (logs/check/l5_oneway/t9b/{z80_A_n4,z80_B_n4,z95_B_n4,
  z95_A_n4_fill0.10}, smoke/solid_z95_A; z95_A_n1 is the worker's 1-chord script test, meaningless, ignore its table row).
  Solid converged A -> B: z 80 psi 0.391 -> 0.390 deg, GJ 6.790 -> 6.794e7, EI_edge 3.043 -> 3.043e9; z 95 psi 2.52 -> 2.516,
  GJ 2.64 -> 2.6465e7, EI_edge 7.87 -> 7.872e8 (< 0.3 %). z 60 fails at 4, 3 and 2 chords: "tip spans two stacks" (a layup
  stack boundary near z 60; the solid builder needs one stack) -> pick another inboard z next session.
  | z | model | psi deg | GJ | EI_flap | EI_edge |
  |---|---|---|---|---|---|
  | 80 | deck | 1.23 | 7.77e7 | 1.53e9 | 2.93e9 |
  | 80 | solid B | 0.39 | 6.79e7 | 1.55e9 | 3.04e9 |
  | 80 | shell (smoke, 6 chords) | 0.52 | 6.38e7 | 1.52e9 | 2.74e9 |
  | 95 | deck | 2.91 | 3.32e7 | 4.06e8 | 8.42e8 |
  | 95 | solid B | 2.52 | 2.65e7 | 4.01e8 | 7.87e8 |
  | 95 | shell (smoke, 6 chords) | -0.68 | 2.77e7 | 4.02e8 | 7.21e8 |
  Reading: EI_flap agrees in all three (<= 2 %). EI_edge: shell/solid 0.90 (z 80), 0.92 (z 95), solid/deck 1.04 / 0.93 -> the
  S-1 EI_edge -9 % is a SHELL deficit. GJ: solid/deck 0.87 / 0.80, shell/solid 0.94 / 1.05 -> the DECK is stiff in torsion
  (D1); the TE filler 0.10 m adds +12.6 % solid GJ (2.98e7) without moving psi, still -10 % vs deck. Principal angle: z 80
  shell ~ solid (deck +0.8 deg off), z 95 deck ~ solid (shell -3.2 deg off); the shell psi moves +-1-2 deg with segment
  length/window -> NOT decided until the shell measurement window is fixed.

- L31 (2026-11-16) Parent, T6 verdict (tools/l5_t6_settled.py, logs/check/l5_oneway/t6/settled.json; SLURM 11741636 done,
  exit 0, 4 h 24 min). Mean over 8 whole revolutions (t 36.5-100 s, 128 dumps, forces rotated to the rotating frame with each
  dump's theta; VTU points are X0+U, same node order). Tip z 116.95, ring-axis metric of the parity tool:
  | solve | flap m | edge m | twist eq deg |
  |---|---|---|---|
  | T3 aero only (t11/bem.json) | 16.911 | 1.630 | -2.25 |
  | T6 mean F_AERO, static, Python K_G | 16.919 | 1.610 | -3.20 |
  | centrifugal only (T6 F_INERTIAL = calc f_cf) | 0.388 | 0.019 | -0.22 |
  | static like-for-like (Python K_G) | 17.306 | 1.628 | -3.42 |
  | static with the Rust K_G | 19.041 | 1.808 | -4.29 |
  | static with no K_G | 19.042 | 1.809 | -4.29 |
  | T6 dynamic mean | 19.041 | 1.808 | -4.29 |
  Criterion (1 %) FAILS: +10.0 % flap, +11 % edge, -0.87 deg twist. Root cause measured: after step 1 the Rust loop rebuilds K_G
  (rotor_fsi.rs ~:590-665, kg_update_interval 0 -> every step, rebuilds logged at steps 1 and 4) with
  assemble_centrifugal_k (aeroelast-core assembler.rs:1170, local sigma = rho w^2 r l_char, the estimate AGENTS.md records as
  2-3 orders too small) and replaces the Python static-solve K_G: ||K_G||_F free DOFs Python 8.53e7 vs Rust 1.21e4 (x7000).
  Production rotor FSI therefore runs WITHOUT centrifugal stiffening (K_SP and f_cf still on). The Newmark integration itself is
  consistent (dyn mean = its own static to 1e-5). Exposed: every RotorFsi run with omega != 0 regardless of
  kg_use_deformed_coords (both branches call the same kernel), incl. L5 two-way. Not exposed: Python-side statics/modal (S-4,
  T1/T3 parity tool). Second finding, open: same K_eff, T6 mean aero load vs T3 _compute_forces(zeros) agree in flap (0.05 %)
  but not in twist (-3.20 vs -2.25 deg; F_x 100.3 vs 104.7 kN) -> production torsional load path / azimuth mean differs.
- L32 (2026-11-16) User, on the T6 fix (A: keep the Python K_G scaled by w^2; B: port the static-solve stress recovery to Rust): "Claro que B, si vamos por A el codigo sera cada vez mas inmanejable". -> S5, T17.
- L33 (2026-11-16) Parent, T17 design (option B, one implementation). Map: the element kernels already exist in aeroelast-core
  (`centrifugal_load` assembler.rs:645, lumped, kg/m2 per element; `compute_prestress_field` :746, centroid (A eps+B kappa)/h in
  Pa; `assemble_geometric_k` :762, Pa, x h inside). The Python `assemble_geometric_stiffness` (core/assembler.py:811-863) does
  load -> PETSc MUMPS solve of K_red on free_dofs -> `assemble_geometric_k_from_disp`; the Rust RotorFsi rebuild
  (rotor_fsi.rs ~:590-665) and rotor_inertial.rs (:765-790) instead call `assemble_centrifugal_k` (local sigma, also fed
  volumetric rho for isotropic materials). Design:
  1. core: `MeshAssembler::mass_per_area_per_elem()` (isotropic/plane-stress rho*thickness, composite mass_per_area); DELETE
     `assemble_centrifugal_k` (core and its PyO3 binding); no fallback.
  2. solvers: new `petsc/elasticity/prestress.rs`, `CentrifugalPrestress`: reduced linear K from `assemble_k` on the free DOFs,
     factorized once (preonly LU, cached like `setup_ksp`); `geometric_k(asm, omega, axis, center, deformed_coords: Option)`:
     f_cf at X0 -> solve -> u_full -> sigma = compute_prestress_field(u) at X0 -> assemble_geometric_k(sigma), on X0+u only when
     deformed coords are requested (K stays linear, "Option A discipline"; prestress always at X0, = the Python reference).
  3. RotorFsi and InertialRotorFsi rebuilds use it (build it once in `new`, keep it in the struct).
  4. PyO3: `PyMeshAssembler.assemble_centrifugal_geometric_k(omega, axis, center, free_dofs)` and `mass_per_area()` call the same
     code; Python `assemble_geometric_stiffness` centrifugal branch calls the binding (its PETSc orchestration deleted);
     `_rho_per_elem` from `mass_per_area()` (the isotropic rescale loop deleted).
  RED first: K_G of the new binding == the current Python orchestration (load -> MUMPS -> from_disp, written inside the test as
  the reference) entrywise on the composite plate of tests/test_kg_prestress_coupling.py and on an isotropic clamped plate with
  a hub offset; a Rust test that the rotor rebuild K_G equals the helper's at X0. Acceptance: T6 tool `static_t6_kgrust` row
  (now through the production path) == `static_t6` within 1e-6 relative.
- L34 (2026-11-16) Worker, T17 implemented per L33 (S5). Commits: 2da4bc1 (Rust: core `mass_per_area_per_elem`, `assemble_centrifugal_k`
  deleted, solvers `petsc/elasticity/prestress.rs` `CentrifugalPrestress`, RotorFsi + InertialRotorFsi rebuilds, PyO3
  `assemble_centrifugal_geometric_k` + `mass_per_area`, Rust tests), 8fe1260 (Python `assemble_geometric_stiffness` -> binding, `_rho_per_elem`
  from `mass_per_area()`, Python tests, tools/l5_t6_settled.py `Kg_rs` from the new binding), d2d70ac (docs/formulations/solvers.md).
  Deviations/notes: (1) RotorFsi builds the helper lazily at the first K_G rebuild (`new` returns Self, PETSc errors become FsiError);
  (2) InertialRotorFsi refreshes the reduced K at X_rotated(theta) on every K_G rebuild (its geometry rotates; a K factorized once at X0 would be
  inconsistent with f_cf(theta)); frame-invariance K(theta)=R K R^T, f_cf(theta)=R f_cf; (3) the old Rust feed used volumetric rho for isotropic
  elements (rho, not rho*t) - gone with the kernel; (4) `rho_per_elem()` (volumetric) binding kept, docs point to `mass_per_area()`.
  RED (observed before the port): Python `pytest tests/test_kg_prestress_coupling.py` -> 3 failed / 7 passed (AttributeError
  `assemble_centrifugal_geometric_k` / `mass_per_area`; the reference orchestration itself ran, ||K_G||>0). Rust: tests added to a scratch copy of
  HEAD (old rotor_fsi.rs + helper) `cargo test -p aeroelast-solvers --lib rotor_kg_rebuild -- --test-threads=1` -> 2 failed, e.g.
  "31109.43 vs 134756.99" (old local-sigma K_G entry vs helper) and "31144.10 vs 134717.93" (deformed variant).
  GREEN: same two tests pass on the new tree; Python 10 passed incl. binding == Python orchestration rtol 1e-10 on the composite plate and an isotropic clamped hub-offset plate.
  Checks (OMP_NUM_THREADS=1, TMPDIR $SCRATCH/tmp/t17): `cargo test -p aeroelast-core -q` 170 passed; `cargo test -p aeroelast-solvers -- --test-threads=1`
  lib 119 passed, fsi_mock_loop 10 passed 2 ignored (the 2 doctests of dynamic_newmark.rs fail identically on HEAD, untouched); PETSc tests need
  `--test-threads=1` and `RUSTFLAGS="-C link-arg=-L<openmpi lib> -C link-arg=-lmpi"` (parallel threads race MPI_Init; no libmpi on the link line);
  `maturin develop --release` Installed aeroelast-py-0.1.0 (56 s); pytest test_kg_prestress_coupling + test_centrifugal_load_density + test_s6_newmark +
  test_inertial_rotor_rust_parity + element/test_prestress_centrifugal_analytic + test_plate_buckling_analytic + test_stress_stiffened_solver
  87 passed 6 skipped; test_iea15mw_s0_omega_zero_consistency 5 passed; test_rotor_centrifugal_branch/inertial_factorization/inertial_fixes 110 passed 1 skipped;
  `scripts/check.sh quick` OK (75 passed, 1 xfailed). Not run: full suite (SLURM, parent), tools/l5_t6_settled.py (parent). Pending: T6 rerun,
  acceptance `static_t6_kgrust` == `static_t6` within 1e-6 relative, and a dynamic mean within 1 % of the T3 static solve.
- L35 (2026-11-16) Parent, T17 acceptance (tools/l5_t6_settled.py on the 11741636 dumps, logs/check/l5_oneway/t6/settled_t17.json):
  ||K_G||_F free DOFs production Rust path 8.53293996901178e7 = Python 8.53293996901182e7; static_t6_kgrust tip flap 17.3063 /
  edge 1.6281 / twist -3.424 deg = static_t6 (identical to 4 decimals). Submitted: full suite SLURM 11743282 (sequana_cpu_bigmem,
  snapshot f8a1639) and the T6 rerun SLURM 11743283 (sequana_cpu_long, 8 h; criterion: dyn mean within 1 % of static_t6
  17.306 m flap, 1.628 m edge, twist -3.42 deg).
- L36 (2026-11-16) User: "arranca con los dos" (T18 aero twist load difference and T9) and "no preguntes mas sobre el estado del git mezclalo todo a main". Parent: main fast-forwarded to acefb5f locally; T18 and T9 to two workers while SLURM 11743282/11743283 run.
- L37 (2026-11-16) Worker T18 (77d7257, tools/l5_t18_*.py, logs/check/l5_oneway/t18/report.md, twist.json). The T6 twist -3.204 deg
  is reproduced from the production config alone (participant built as run_bem_fsi.main: allOuterShellNods coupling set,
  36-azimuth mean, u != 0 branch): -3.205 deg, F 100.31/842.29 kN. T3 -> T6 chain on the same K_eff: wall-flow realisation off in
  production (projector gets element_properties=None) -0.311; AeroDyn -> windIO polars + station grid (53 vs 50) -0.690; azimuth
  mean (tilt -6) -0.042; u != 0 branch +0.057; coupling node set (no web nodes) +0.035 deg. Ruled out: Cm on/off (both apply it,
  totals equal 0.1 %), AC fraction swap (< 0.001 deg), BEMSolver kwargs, deformed feedback. FINDING (production bug):
  `_rebuild_bem_solver` (fsi_participant.py ~:1033) rebuilds AeroStations WITHOUT prebend; production always takes that branch
  (max|u| > 1e-12): thrust -1.1 %, torque -1.9 %, power -1.9 % vs the prebend-carrying solver. Also: `include_pitching_moment`
  is never read (Cm always on). Intended load for T6: the production load; the parity tool's `shell_loads("bem")` must be
  re-pointed at the production participant (T3/T1 BeamDyn rows use binned resultants; T6 uses the subset).
- L38 (2026-11-16) Worker T9 (3c293a5, logs/check/l5_oneway/t9c/summary.md), PARTIAL. Shell measurement: bias-corrected
  (deck read through the same fit), segment-length criterion psi <= 0.5 deg, stiffness <= 2 % -> 8 chords, window L/4 (z80
  still +0.45 deg 8->10 chords); shell psi oscillates +-1.5 deg along z (wavelength ~6 m, cause open, absent in the deck read).
  Inboard station z45 (one TE stack z 35-62). Model - solid (level B, 4 chords):
  | z | shell dpsi | deck dpsi | shell C_fe/C_ff (solid) | GJ shell/deck vs solid | EI_edge shell/deck vs solid |
  |---|---|---|---|---|---|
  | 45 | +1.02 | +4.00 | -0.059 (-0.062) | solid GJ not length-converged (-11 % 4->6 chords) | 1.06 / 1.12 |
  | 80 | +1.45 | +0.84 | +0.016 (+0.003) | 0.95 / 1.14 | 0.89 / 0.96 |
  | 95 | -2.60 | +0.40 | -0.024 (+0.022, sign wrong) | 1.08 / 1.25 (1.11 TE filler) | 0.93 / 1.07 |
  Reading: no single winner on flap-edge coupling (inboard shell ~ solid, deck 4 deg off; outboard deck ~ solid, shell 2.6 deg off
  with the wrong cross-term sign at z95). EI_edge shell deficit -11/-7.5 % outboard, +6 % at z45 (not uniform). Deck GJ excess D1
  confirmed outboard (+14/+25 %). Open: z45 solid torsion convergence; cause of the shell psi oscillation; band 1 deg set after data.
- L39 (2026-11-16) Parent, T18 fix: `_rebuild_bem_solver` now carries `prebend` and `sweep` of every station
  (fsi_participant.py). RED first: tests/validation/rotor/test_bem_rebuild_keeps_prebend.py (rebuild at the reference geometry
  must keep station prebend/sweep and reproduce the construction-time Np/Tp at rtol 1e-12) failed 2/2, then GREEN. Callers of
  the rebuild (test_bem_fsi_deformed_geometry, test_blade_deloading_vs_reference, test_torque_sense_vs_openfast, nan inputs)
  24 passed; check.sh quick OK. RDD push left for later (user: "si no puedes dejamos el rdd para despues y cerramos la fisica
  primero"); the earlier SLURM 11743282/11743283 were cancelled while pending and are resubmitted on this tree (L40).
- L40 (2026-11-16) Parent. Resubmitted on 75d1913 (K_G T17 + prebend fix): full suite SLURM 11743498 (sequana_cpu_bigmem,
  snapshot logs/check/git-snapshot-20261116T222722, pending AssocMaxJobsLimit) and T6 rerun SLURM 11743499 (moved to
  sequana_cpu, RUNNING since ~22:30, 8 h limit, output $SCRATCH/l5_oneway/full-11743499/). Independent verify of T17: pass with
  notes, no defect (18 Python tests pass; dense-reference probe 6.5e-15; X0 restore on all paths; one implementation). Notes for
  after the suite: delete the dead l_char kernels `compute_centrifugal_prestress` (mitc3.rs:1596, mitc4.rs:3012, own tests
  only); add regression tests (dense reference + omega^2 scaling with a non-unit axis; out-of-range free_dofs -> ValueError);
  Python K_G now uses PETSc default LU instead of MUMPS (production runs keep MUMPS via PETSC_OPTIONS). RDD: push of main to
  origin blocked by the RDD gate (review due, 18.9k lines); `gentle-ai` is not on PATH in either shell; user: leave RDD for later.
- L41 (2026-11-17) Parent. T6 rerun SLURM 11743499 (75d1913: K_G T17 + prebend fix) COMPLETED in 4h06;
  `tools/l5_t6_settled.py` -> logs/check/l5_oneway/t6/settled_rerun.json. Tip, run's own mean loads, same K_eff:
  | row | flap m | edge m | twist eq deg | twist skin deg |
  |---|---|---|---|---|
  | static_t6 | 17.1470 | 1.6039 | -3.433 | -3.463 |
  | dyn_mean | 17.1990 | 1.5984 | -3.454 | -3.482 |
  | delta | +0.30 % | -0.34 % | +0.61 % | +0.55 % |
  | static_t6_nokg | 18.8669 | 1.7822 | -4.295 | -4.329 |
  Criterion (1 %, T6) met on all four: the production run now carries the centrifugal stiffening (L31's +10 % flap
  is gone; static_t6_kgrust == static_t6). T6 done. Full suite SLURM 11743498 TIMEOUT at 4 h (68 %, 2 F, ids not in
  the -q log); resubmitted 11744840 with --time=600 on sequana_cpu, snapshot logs/check/git-snapshot-20261117T085709
  (048bf2f). Handoff step 3 started: tools/l5_oneway_bd_parity.py gains `--loads prod` (`production_loads`: the
  participant built as run_bem_fsi.main, rebuilt branch, 36-azimuth mean, coupling-mesh forces written onto the full
  mesh by vertex match, tilt/precone checked against the config); `bem` keeps the pre-T18 T3 load for the earlier
  records and its other callers (t6_settled, orientation_audit, twist_sense_chain); production_cfg/FLUID_YAML/N_AZ
  moved there from l5_t18_loads (re-exported). T3 rerun: logs/check/l5_oneway/t3prod/.
- L42 (2026-11-17) Parent. T3 rerun with the production load: `tools/l5_oneway_bd_parity.py --loads prod`
  (048bf2f + uncommitted tool change) -> logs/check/l5_oneway/t3prod/prod.json (5.7 min, login). Conservation nodal = rows
  (F 99.16/836.90 kN, M_root exact). Independent verify of `production_loads`: pass (recipe parity 6e-13; `bem` callers
  unchanged); t18/subset.npz is stale (pre-L39, ~1 % Fx), re-record before using it as a reference.
  Shell vs BeamDyn (x1, nonlinear), same station resultants:
  | z m | 35.1 | 58.5 | 70.2 | 81.9 | 93.6 | 105.3 | 111.2 | 117.0 |
  |---|---|---|---|---|---|---|---|---|
  | flap d% | -6.0 | -4.8 | -4.9 | -4.7 | -4.1 | -3.0 | -2.6 | -2.0 |
  | edge d% | -1.3 | +8.6 | +10.9 | +12.1 | +12.6 | +12.1 | +11.1 | +11.0 |
  | twist eq shell-BD deg | 0.00 | +0.53 | +0.89 | +1.13 | +0.80 | +0.23 | +0.05 | +0.21 |
  Tip: flap 16.760 vs 17.102 m, edge 1.585 vs 1.428 m, twist -3.214 vs -3.428 deg (skin -3.244). The shell tip twist
  equals T6 aero_t6 (-3.213): T3 and T6 now carry the same load. Against the pre-T18 load (t11/bem.json): flap and edge
  ratios unchanged to 0.5 %; the outboard twist difference loses its sign flip (z105-117: -0.97/-0.53/+0.99 -> +0.23/+0.05/
  +0.21 deg); mid-span z58-94 stays +0.5 to +1.1 deg. Verdict against 5 % (S3): flap pass at z >= 58.5 (-6.0 % at z35, where
  the deflection is 0.6 m); edge red (+11-13 %, the flap-edge cross-coupling of T9, not load-caused); twist red mid-span,
  tip within 0.21 deg (6 %). Remaining T3 reds belong to T9 (coupling) and T4/T10 (twist metric and sense).
- L43 (2026-11-19) Parent. [L5 paused after this entry: user chose to close the gmsh-official-mesher gates and retire NuMAD first, gmsh-official-mesher L-S37.] Full suite 11744840 (048bf2f, 9h24): 45 failed, 2097 passed, 8 errors; against
  known-red-full.txt (recorded 2026-10-28 at bd3e73f, 236 commits back) 27 NEW-FAIL. Triage (explorer + parent):
  (i) 11 x test_ring_kinematics_estimator: TypeError, `_ring_section_rotation` gained `edges` in T19 (6a6ad35)
  and the catalog adapters were not updated -> fixed (closed-ring edges for the synthetic rings), 14 passed 2 xfailed.
  (ii) v2 copies of the known-red v1 accuracy tests (S-1 EI_edge/GJ/EA, S-7, V-02 1st flap/edge): new ids, same
  documented reds (v2 GJ 12.1 % vs v1 9.7 %) -> re-record, not a regression. (iii) mesher group (test_blade_mesh
  counts, red since 2026-10-17; gmsh_chord_scaled 26015 vs 26018 elements; windio2 oblique y flip; web nd_arc offset
  1.4-1.8 cm; quality gate h0.25 web0/web1 coincident at z 9.53-9.58): mesher workstream, not L5. (iv) ERROR x4
  frame_tear/torsion_stiffness fixture `float(None)` (logs/check/triage-frame-tear.log): a 7-node station ring
  (nodes 2809-2813, 2844-2845) whose skin edges leave < 3 positively weighted nodes, so `ring_section_rotation`
  returns None (section_contour.py:600). The numad mesh now yields a partial station ring; likely the same mesher
  change as the web0/web1 coincidence at z 9.55 (z of the ring not printed) -> mesher group, not L5. (v) ERROR x4 S-2 utd
  `Laminate must have at least one ply`: disowned input. (vi) S-2 official LC1 interior 5.01 % vs 5 %.
  None of these touches the T17 K_G change (rotating path only).
  T9d(b) psi oscillation, logs/check/l5_oneway/t9d/osc/analyze.md (login, 2026-11-17): the bias-corrected shell psi
  oscillation is mesh-converged (z80 cd 15/20/30: 1.52/1.44/1.48 deg at half 1, wavelength 8-9 m), independent of the
  segment length (nchord 6/8/10: 1.56/1.44/1.45 deg), absent from the deck same-fit (<= 0.05 deg), uncorrelated with
  the mean layup thickness (|corr| <= 0.16), shrinks with the fit window (half 1 -> 4: 1.44 -> 0.48 deg) and grows
  outboard (z45 0.33, z80 1.44, z95 4.49 deg). Reading: real shell content at the 4-10 m scale that the deck's
  station-wise 6x6 does not carry; not a discretization or fit artifact. Cause not yet named.
  T9d(a) solid z45: 11744952-55 OOM (ntasks 3 = 24 GB, three ccx cases at once); resubmitted ntasks 12 as
  11750509-12 (nchord 8/10/12/16), old outputs kept as z45_A_n*.oom.

- L44 (2026-11-19) T9d(a) solid z45 READ (SLURM 11750509-12, ntasks 12, 19-42 min, all ccx cases ran;
  logs/check/l5_oneway/t9d/solid/slurm-1175050{9..12}.out). Level A, NQ 813, chord 4.786 m, segment length L:
  | nchord | L m | psi deg | fe/ff | EIf N m2 | EIe N m2 | GJ N m2 |
  |---|---|---|---|---|---|---|
  | 8 | 38.29 | -7.537 | -0.0616 | 9.9755e9 | 1.8364e10 | 2.3908e8 |
  | 10 | 47.86 | -7.532 | -0.0615 | 9.9757e9 | 1.8368e10 | 2.2448e8 |
  | 12 | 57.43 | -7.529 | -0.0615 | 9.9757e9 | 1.8369e10 | 2.1335e8 |
  | 16 | 76.57 | -7.529 | -0.0615 | 9.9757e9 | 1.8369e10 | 1.9960e8 |
  Bending, principal angle and flap-edge coupling are length-converged (<= 0.03 %, 0.008 deg). GJ is NOT: -16.5 %
  from L 38 to 77 m, still falling; GJ = a + b/L fits all four within 0.3 % and extrapolates to 1.60e8,
  1/GJ = a + b/L gives 1.71e8 (fit within 1.3 %); the 2D St-Venant value of the same section is GJ_st2d 1.891e8.
  Reading: the solid's torsion carries an end (warping-restraint) stiffening that decays slowly with L; the solid
  GJ at z45 cannot arbitrate torsion until its limit is pinned (longer segment or free-warping ends). The T9
  coupling question (edge +11-13 %) is answered on the converged quantities: fe/ff -0.0615 at z45.
## Handoff (2026-11-16 night, end of session)
Resume: this file to `## Log`, then L31-L40; one narrow mem_search "T17 K_G rebuild prestress".
State: branch feat/windio2-support = local main = HEAD (origin/main still at 5acbbd8; push blocked by RDD, see L40).
Login node: ulimit -t 1800 CPU-s per process; OMP/OPENBLAS=1. Do not edit src/ or crates/ while 11743498 is pending/running.
Next, in order:
 1. T6 verdict on SLURM 11743499: `scripts/aeroenv.sh python tools/l5_t6_settled.py --results $SCRATCH/l5_oneway/full-11743499/results
    --out logs/check/l5_oneway/t6/settled_rerun.json`; criterion: `dyn_mean` within 1 % of `static_t6` (flap, edge, twist) with
    the run's own mean loads (prebend fix changes the loads, so compare to the new static_t6, not 17.306).
 2. Full suite 11743498: compare with scripts/known-red-full.txt; if green, T17 cleanup (L40 notes), maturin, quick.
 3. Re-point tools/l5_oneway_bd_parity.py `shell_loads("bem")` at the production participant (T18, L37), rerun T1/T3 rows.
 4. T9 open items (L38): z45 solid torsion length convergence, cause of the shell psi oscillation (~6 m wavelength).
 5. L3 re-closure decision (user), then T4 (twist estimator), T5 (store rows), close step 4.
 6. RDD: user locates `gentle-ai` and runs `gentle-ai review mode disable --scope clone`, or reviews T17 range 59bc420..f8a1639; then push main.
