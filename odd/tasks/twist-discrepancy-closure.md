# Feature: close the tip-twist discrepancies with a finite plan

Written 2026-10-13 after the beam-vs-ring result (`shell-twist-mixed-mesh.md`, baf38a3). Goal: stop
the probe-by-probe loop. Every step below ends in a decision, has a control that must pass before
its number is read, and has a stop rule. Nothing is tuned toward Zhou. Blade mesh 0.15-0.25 m.

## Why we kept iterating

Three different discrepancies were being chased as one, and no step had an independent judge:

| id | discrepancy | numbers (tip, nose-down +) | class of question |
|---|---|---|---|
| D1 | shell vs beam, same load | shell +5.855 vs decoupled beam +2.079 (2.82x; ~2.4x per segment, uniform from the root) | structural model or load path |
| D2 | our beam vs Zhou's GEBT, Zhou's load | +2.079 vs +3.60 (1.73x) | beam class, geometric nonlinearity, or the reconstructed Mp |
| D3 | our coupled FSI vs Zhou | +5.589 vs +3.60 (1.55x) | follows from D1 and D2: the shell under Zhou's load gives +5.855, so the load is not the source |

D3 is not a separate question: once D1 and D2 are attributed, D3 is arithmetic. Issue #10 (parked
divergence) is a fourth, unrelated question (dynamic stability), see "Not in this plan".

## Facts this plan rests on (verified 2026-10-13)

- The BeamDyn deck comes from SONATA/ANBA4 (`IEA-15-240-RWT/ReleaseNotes.md:25`). ANBA4 is a 2D
  FE cross-section solver (VABS class), so the deck's GJ and its shear-centre offsets (K16/K26,
  max |K16|/K66 = 1.49) are an independent section judge we already hold.
- K46 = K56 = 0 exactly at all 26 deck stations (parsed). An exact zero in a composite blade
  suggests the bend-twist terms were dropped in export, not that they are physically zero
  (inference, unverified). The deck therefore cannot judge bend-twist coupling.
- Global GJ: S-7 after the kappa fix (72f64de) gives shell/beam twist ratio 1.047 at 0.25 m, i.e.
  shell GJ ~4 % below the deck under a pure tip couple (`docs/validation_closures.md:982-988`).
  S-1's +23.9 % is at 0.5 m (not converged, rule 9). Torsional stiffness cannot carry a 2.4x.
- Shell implementation = CalculiX S8R within 1.46 % on the same nodal force (Line A); #34 and #38
  ruled out (E3, Line B); the gap is not section distortion (beam-vs-ring, rigid 75 % ring rho 0.406).
- The shell load (`force_projection.py:784-805`) and the beam's torque
  (`m = Mp + ((pitch-0.25)c - x_shear) Np`, Tp arm zeroed, no twist projection) are not shown to be
  the same torque. D1 therefore still mixes "structure" with "load path".
- `beamdyn_driver` v5 works (`tools/diagnose_radii_beamdyn.py`): static, point loads per strip
  binned from the shell's own nodal forces about the reference axis; it reads displacements only.

## Tasks

- [ ] **P0. Pre-register P1** (doc commit before any tool edit): the factors, the control, the bands.
- [ ] **P1. Split D1 into a load-path factor and a structural factor, and get D2's beam-class
      answer, in one run (~10 min, `sequana_cpu_dev`).**
      On the existing Zhou-load shell solve:
      1. Bin the shell's nodal force vector into strip loads about the deck reference axis
         (reuse `diagnose_radii_beamdyn.shell_state`), and run BeamDyn with the full deck (K16/K26
         kept): `BD_shell_load`. Add `RDzr` to the OutList, convert Wiener-Milenkovic to an angle.
      2. Run BeamDyn on the Zhou-formula load (Np, Tp, Mp per strip, the beam's own arm):
         `BD_zhou_formula`.
      3. Both at full load and at 1 % load rescaled x100: `BD_lin` vs `BD_nl` gives the geometric
         nonlinearity factor (the shell solve is linear, Zhou's GEBT is not).
      Factors at the rigid 75 % ring and at the tip:
      - structural = `omega_shell / BD_lin(shell_load)` (same force, two structural models);
      - load path = `BD_lin(shell_load) / BD_lin(zhou_formula)` (same model, two loads);
      - geometric = `BD_nl / BD_lin`;
      - beam class = `BD_lin(zhou_formula) / our decoupled beam` (must be ~1 if our beam is right).
      **Control, must pass before any factor is read:** the S-7 pure tip couple through the same
      BeamDyn pipeline reproduces the shell/beam ratio 1.047 within 3 %. It tests frames, units and
      the WM conversion. If it fails, the run is a harness bug and nothing is reported.
      **Decision:** structural factor <= 1.25 (inside the 14-24 % literature spread) -> D1 is load
      path; name which arm term (Tp, twist projection, AC point vs pitch axis) and stop D1.
      Structural > 1.5 -> P2. In between -> report, P2 only with an OK.
      D2: if `BD_nl(zhou_formula)` is within 10 % of 3.60 -> Zhou reproduced on its own model class,
      D2 closed. If not -> the residual is Zhou's load reconstruction or properties we do not hold;
      record it as not resolvable here (Mp error is unquantified) and close D2 as such.
- [ ] **P2 (only if P1 says structural). Name the section property (~15 min).** Tip load cases on
      the shell, reusing `run_s7_torsion`/`run_s1_sectional` (kappa-corrected couple), each also
      through BeamDyn: (a) pure tip torque (control, known 1.047); (b) pure tip flap couple -> twist
      from bend-twist coupling only (deck gives 0); (c) tip flap force at two chordwise points ->
      the shell's shear centre near the tip vs the deck (ANBA4) x_shear. Prediction check: the
      named property must reproduce the P1 structural factor within 15 % when put into the beam,
      otherwise it is not the cause.
- [ ] **P3 (only if P2 names it). Judge that one property.**
      Shear centre or GJ -> the deck (ANBA4) is the judge; a shell mismatch is a shell finding,
      propose the fix and wait for an OK. Bend-twist coupling -> the deck cannot judge (K46/K56
      zeroed). Judges in order of cost: the full ANBA4 6x6 if published (SONATA outputs, IEA-15
      history); a layered 3D solid segment in CalculiX at z ~ 88 m; installing a cross-section
      solver (SONATA/ANBA4 or OpenSG, hours). Pick only the one needed.
- [ ] **P4. Close.** One attribution table for D1/D2/D3 into `docs/validation_closures.md`, each
      line either matched by a judge within its band or a named model-class difference with its
      bound. After P4 the question is closed: reopen only if a re-run trigger in that doc fires.

## Stop rules

- A step whose control fails reports the failure and stops; no new hypothesis is opened from it.
- No step adds a hypothesis that is not already one of P1-P3. A surprise is written down as a
  follow-up, not chased in the same session.
- One step per session; the next session starts from this file's first unchecked box.

## Not in this plan, and why

**Parked FSI aeroelast vs ccx_preCICE.** It cannot decide D1-D3: at pitch 90 the load is near
zero (CT ~ 0.008) so bend-twist barely activates; CCX S8R shares the CLT and the geometry writer
with our shell and already agrees to 1.46 %; D3 is explained by structure (the shell under Zhou's
load alone gives +5.855), so the coupling layer is not where the twist gap lives. It is the right
tool for issue #10 (settle vs diverge), as its own track (`ccx-fsi-parity-parked-v50.md` T5-T6),
with its blocker first: the `*SHELL COMPOSITE` deck diverges under load, move to `*SHELL GENERAL
SECTION` (A/B/D) and re-gate; cost today ~9 h per 40 windows at 0.5 m, 0 windows at 0.25 m.

## External inputs (articles and software), ranked by what they decide

Located 2026-10-13 (web search); none is installed or read yet unless stated.

| # | input | decides | needed at | cost / access |
|---|---|---|---|---|
| 1 | **Zhou et al. 2025 (Energy), structural-model section** | which sectional properties their GEBT used (this deck? a coupled 6x6? their own VABS?), linear or nonlinear, and their twist along the span (not only the tip) | P1 decision on D2 | paywalled; ask the user for the PDF. Without it D2 can only be closed as "not resolvable here" |
| 2 | **Gaertner et al. 2020, NREL/TP-5000-75698** (IEA-15 definition), blade structural section | whether K46/K56 were computed by ANBA4 and dropped on export, and the shear-centre/elastic-axis data | P3 if coupling | open: docs.nlr.gov/docs/fy20osti/75698.pdf |
| 3 | **Feil et al., SONATA paper** (osti.gov 1660247) | ANBA4 validated against VABS: tells how far the deck can be trusted as a section judge | P3 | open |
| 4 | **SONATA + ANBA4** (github NLRWindSystems/SONATA; gitlab.lrz.de HTMWTUM/SONATA) | regenerates the full 6x6 (K46/K56 included) from the same WindIO yaml: the judge for bend-twist coupling and shear centre with the same input as the shell | P3 if coupling or shear centre | open source; ANBA4 needs legacy FEniCS/dolfin (conda), likely hours to install on the cluster |
| 5 | **OpenSG** (github wenbinyugroup/OpenSG; Sandia/Purdue, MSG with FEniCS, verified vs VABS/SwiftComp on a wind blade) | same as 4, independent implementation (Yu's group, the same theory family as VABS and GEBT) | P3, alternative or cross-check to 4 | open source; FEniCSx |
| 6 | **ALM-GEBT paper, J. Mar. Sci. Eng. 13(1):116 (2025)**, structural nonlinearity effects on IEA-15 | the literature's own linear-vs-nonlinear twist factor for this blade (cross-check for the P1 geometric factor) | P1 | open access (MDPI) |
| 7 | Frontiers Energy Res. 2025 (fenrg.2025.1571567 yaw; 1591837 swept), LL-FVW + GEBT on IEA-15 | more twist values and, possibly, the structural input used by the same model family | P1/P4 context | open access |

Not proposed: VABS, BECAS (licence, and no MATLAB on the cluster: `module avail matlab` empty);
GEBT itself (BeamDyn is the same geometrically exact beam class and is already installed).
