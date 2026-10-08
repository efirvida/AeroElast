# Feature: the +31% tangential load (roadmap item 3, issue #15)

Status: T1-T4 closed 2026-10-08; T5 in progress (verification done, the close is pending)
Owner: this session
Related: issue #15 (roadmap item 3 of #18); issue #14 (item 2, closed as documented
non-transferability and building on the same Zhou publication); `docs/validation_results_v01_v02_v04.md`
§5.9 and §5.11; `docs/validation_plots/plot_spanwise_loads.py`; `tools/diagnose_zhou_loads_reverse.py`.

## Why this exists

Issue #15 states that our tangential load integral runs `+31%` above Zhou et al. 2025's while the
normal force and the pitching moment match to 1%, and that the difference "has not been attributed:
neither to the polar used, nor to the tangential induction, nor to the definition of the tangential
direction in their frame". Its own closure path is explicit:

> Reproduce the tangential integral from their published angle of attack and the official polars
> exactly as the pitching moment was, then compare term by term: polar `Cd`/`Cl` projection, the
> tangential induction factor, and the sign/direction convention of the tangential axis.

This feature runs exactly that path and lands the attribution. It is offline and deterministic:
`BEMSolver.compute` is production (`src/aeroelast/solvers/bem/engine.py`), the polars are the official
AeroDyn deck, and the reference is the digitized Fig. 11 already in the tree.

## What is already established (read-only, 2026-10-08)

Measured with the production BEM at Zhou's own rated point (`V = 10.59 m/s`, `Omega = 7.55 rpm`,
`pitch = 0`, official AeroDyn deck, 50 stations):

```
thrust 2.481 MN   power 15.723 MW   torque 19.861 MNm
rotor-plane frame: intNp 922.2 kN   intTp 127.8 kN   int(Tp r) 8.786 MNm
section frame    : intNp 925.4 kN   intTp 117.3 kN   int(Tp r) 9.251 MNm
Zhou digitised   : intNp 732.2 kN   intTp  70.0 kN   int(Tp r) 4.923 MNm
```

1. **The `+31%` is a pointwise statistic, not an integral.** In the §5.9 write-up the number is the
   spanwise `Tp` comparison (peak `1.05` against `0.80 kN/m` = `1.31`, and the mean pointwise ratio over
   the digitized stations is `1.300`, measured on the only surviving local campaign, which is
   low-load). The *integral* ratios differ by how you pick them: `1.24`
   (that campaign, `40-100 s`), `1.49` (§5.9's stored, non-reproducible numbers), `1.83` (the rigid
   production BEM). The issue's title overstates a pointwise deviation as an integral.
2. **Zhou's Fig. 11 is the section (chord) frame, not the rotor plane.** Taking their published pair
   (`Fig. 11` `Np`, `Tp`) and the official polar at their published `alpha` (`Fig. 10`), the identity
   `psi := atan2(Tp, Np) = alpha - gamma` (with `gamma := atan2(Cd, Cl)`) reproduces `Fig. 10`'s
   `alpha` to **`+/-0.4 deg` over `r/R in [0.26, 0.80]`**. The rotor-plane reading
   (`alpha = psi + gamma - theta`, with `theta` the blade twist) is wrong by **up to 9 deg**. Our
   `BEMSolver` output is rotor-plane **by construction** (`theta_implied` reproduces the deck twist
   exactly), so §5.9 compares a rotor-plane curve against a section-frame curve.
3. **The digitised curve does not carry the paper's own integrals.** `intNp dr x 3 = 2.20 MN`, which
   equals Zhou's reported **flexible** thrust `2.20 MN` (`Table 6`) and not the rigid `2.53 MN`, so
   `Fig. 11` is the flexible case. Its `int(Tp r) x 3 = 14.8 MNm` against a torque implied by their
   `14.76 MW` of `20.4 MNm`, i.e. the digitised tangential torque is **~27% low**. Integral
   comparisons through the digitization therefore carry a `~10-30%` quadrature band.
4. **§5.9's own numbers are not reproducible from the tree.** The campaign it cites
   (`frontiersin_results_corotational`) is gone; the surviving `frontiersin_results_corotational_100s`
   is a different, low-load run (`thrust 1.64 MN`, `power 12.34 MW`, tip `7.4 m` against the article's
   `14.71 MW` / `12.73-12.79 m`), so the committed figure and its parentheses are stale.

## Tasks

- [x] **T1 - The frame of Zhou's Fig. 11, established by measurement, and the polar term.**
  **Done 2026-10-08** (`51f9718`): `tools/diagnose_zhou_tp_frame.py`; the production BEM's
  rotor-plane identity holds to `1.4e-14 deg`, the section reading of their published pair
  reproduces `Fig. 10`'s `alpha` to `0.37 deg` over the blade core (`r/R 0.26-0.80`, 10 stations)
  where the rotor-plane reading is off by `6.83 deg`; guard
  `tests/validation/bem/test_bem_load_frame.py`, declared `out_of_scope` in `groups.yaml`. The
  polar enters both readings identically, so it is not the cause.
- [x] **T2 - The digitization quadrature band, and which case Fig. 11 is.** **Done 2026-10-08**:
  the digitized `intNp dr x 3 = 2.20 MN` equals `Table 6`'s **flexible** thrust (not the rigid
  `2.53 MN`), and its tangential trapezoid carries `11.68 MW` of the paper's own `14.76 MW`
  (`-21%`). Pointwise digitization band `+/-0.05 kN/m` (`+/-3-5%`), `Fig. 10` band `+/-0.5 deg`.
- [x] **T3 - Re-do the §5.9 comparison in one declared frame, on the production path.** **Done
  2026-10-08**: the frame is declared in `docs/validation_plots/plot_spanwise_loads.py` and the
  rotor-plane pair is rotated into the section frame before the overlay; §5.9 and the `V-09` entry
  carry the attribution and keep the old numbers as non-reproducible record. The residual is our
  inflow-angle bias (the issue's second term), measured in one frame and one `qc` with the polar
  re-evaluated at each angle: feeding our `alpha` into their load formula multiplies `Tp` by
  `1.19-2.24` and `Np` by `1.06-1.34` over `r/R 0.26-0.80` (direction ratio `Tp/Np` `1.12-1.67`).
  `Tp` is the small difference of two large terms, so a 1-5 deg inflow shift moves it 3-7x more in
  relative terms than `Np`; the absolute bias is the 0.9-to-4.8 deg that §5.11 already documents.
- [x] **T4 - Land the attribution: write-up, store, closures.** **Done 2026-10-08** (`5bfad2a`):
  §5.9 and the `V-09` entry rewritten with the frame declared, `docs/validation_plots/plot_spanwise_loads.py`
  rotates before overlaying, `gaps.yaml` gains `zhou_spanwise_load_frame` (`not_validated`,
  `citations_forbidden: true`), and `docs/validation_closures.md` carries the #15 section and the
  re-run trigger. The report's old numbers stay as non-reproducible record. `check` 198 / 259 / 0.
- [ ] **T5 - Independent read-only verification, then close.**
  A separate verifier re-derives the frame identity and the quadrature band from the raw CSVs and the
  deck.constants, checks every printed number against the live script, and checks the store. Then
  comment on #15 with the attribution and close it, and update the #18 map.
  **Verification done 2026-10-08** (independent read-only, `gentle-ai-verify`): frames PASS
  (`0.366 deg` core section reading vs `6.826 deg` rotor, re-derived from the CSVs and the deck; the
  `1.421e-14 deg` identity re-checked against ccblade's `:734-739`), quadrature PASS (`2.197 MN`,
  `11.676 MW`, and the conclusion is not `Rtip`-mapping-sensitive over 117-122.5 m), gates PASS
  (`check` 0 errors, `status` 0 undeclared, 2 tests pass). Two corrections it forced and this unit
  applied: the `Tp`/`Np` pair mixed two polar variants (now one variant: `Tp` `1.19-2.24`, `Np`
  `1.06-1.34`, direction `1.12-1.67`), and the frame effect's `-7%` is range-dependent (full span
  `0.20-0.98`; `-15%` over the core `0.26-0.80`, so the reported `-7%` is tip-dominated) - both
  ranges are now stated in the closure. The verifier also noted that the issue's own "Np matches to
  1%" comes from the reconstructed-`Mp` comparison, not from the inflow experiment, which the prose
  now says explicitly.

## Follow-ups this unit opened

- **The `ForceProjector` load frame, on our side.** The plan documents `ccblade`'s `Np`/`Tp` as
  "normal and tangential to the section chord"; the engine emits the *rotor-plane* pair
  (`cn = cl cos(phi) + cd sin(phi)`), so `ForceProjector.project` puts rotor-plane components on the
  section's own chord-normal / chord axes, a rotation by the local twist (up to 15.6 deg inboard).
  The net effect on the integrated force is small and of both signs along the span, but it is the
  same class of implicit-convention defect as the P5 sense chain, and it should be arbitrated with a
  case whose answer is known before attributing anything to it. Not this unit's scope: the reported
  numbers compared `bem_sectional.csv`, which is upstream of the projector.
- **Their `Tf`/`Np` digitization is the flexible case.** Any future integral comparison must match
  case and quadrature, or be stated as a bound.

## Constraints

- Measure on the **production path** (`BEMSolver`, the official AeroDyn deck), never a test-local
  reconstruction of the rule.
- One work-unit commit per task, Conventional Commits, on `integrate/origin-main-2026-09-30`, no
  amend (shared branch). No push unless asked.
- Diagnostics outside the repo (`$SCRATCH/tp_diag/`); cited evidence in the repo.
- Test-first by default where a runnable deterministic check exists (the frame identity and the
  quadrature identity are both assertable); otherwise state why and run structural verification.
- Never invent metadata (their pitch, stiffness source, or torsion definition are unpublished): declare
  the gap.
- `tools/tests` must be run with `-m "not slow"` (issue #24: `coherence` is destructive).

## Evidence pointers

- `docs/validation_data/zhou_2025_fig11_loads.csv`, `docs/validation_data/zhou_2025_fig10_aoa.csv`
- `docs/validation_results_v01_v02_v04.md` §5.9 (`:769`), §5.11 (`:800`), §5.11.1 (`:811`)
- `docs/validation_plots/plot_spanwise_loads.py`, `tools/diagnose_zhou_loads_reverse.py`
- `src/aeroelast/solvers/bem/engine.py` (`Np`, `Tp`, `qc`, `Mp`), `tests/support/openfast_bem.py`
- Probes (outside the repo): `$SCRATCH/tp_diag/{probe.py,mapping.py,frame.py,frame2.py}` and their logs
