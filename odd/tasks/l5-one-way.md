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
- T2 | S1 | parent | pending | S-6 NaN: reproduce on HEAD (mid-surface vs default OML skin), root cause, close or declare
- T3 | S1 S3 | parent | pending | our rigid BEM loads (deck tilt -6 / precone -4, T7 model) through the production load mapping onto the shell, the same loads onto BeamDyn; load conservation check; 5 % on deflection and twist
- T4 | S4 | parent | pending | twist metric for the one-way deliverable: equal-weight vs skin-weighted estimator on T1/T3, decide or carry to two-way
- T5 | S1 S2 | parent | pending | store rows + verdicts, thesis plan T4 closed or declared
- T6 | S1 S3 | worker | pending | S-5 rated run with the production code: `aeroelast-bem-fsi` one-way (`deformed_twist: false`, `deformed_radius: false`), gmsh N20 AR 1, deck tilt -6 / precone -4, rated 10.59 m/s 7.56 rpm, ~100 s; its settled mean tip deflection and twist must match the T3 static solve (same loads, same K_eff) within 1 %
- T7 | S3 | worker (with T6) | pending | tilt/precone explicit in every L5/L5-two-way case file: deck -6/-4 for OpenFAST comparisons, 0/0 only for Zhou cases, labelled; no case relies on the 0/0 default
- T8 | S2 | worker | pending | torque sensitivity harness: a measured spanwise twist difference (shell - BeamDyn) -> delta torque / power / thrust at rated with the production BEM (reuse tools/l3_torque_sensitivity.py), ready for any red of T1/T3

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
