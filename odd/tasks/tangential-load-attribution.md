# Feature: the +31% tangential load (roadmap item 3, issue #15)

Status: started 2026-10-08 (T1-T5)
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
   the digitized stations is `1.300`). The *integral* ratios differ by how you pick them: `1.24`
   (the only surviving local campaign, `40-100 s`), `1.49` (§5.9's own numbers), `1.83` (the rigid
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

- [ ] **T1 - The frame of Zhou's Fig. 11, established by measurement, and the polar term.**
  Turn the read-only probe into a repo artifact (`tools/diagnose_zhou_tp_frame.py` or an extension of
  `tools/diagnose_zhou_loads_reverse.py`) that prints the term-by-term table for both frame
  hypotheses; pin the section-frame verdict and the polar `gamma`; state what remains unmeasured (their
  `a'`, their twist/pitch) rather than inferring it. RED first: a check that the probe's own identity
  holds on our BEM output before it is trusted on theirs.
- [ ] **T2 - The digitization quadrature band, and which case Fig. 11 is.**
  Compute each digitized curve's implied rotor thrust and torque against the paper's own `Table 6`
  (`rigid 16.11 MW / 2.53 MN`, `flexible 14.76 MW / 2.20 MN`); pin that Fig. 11 is the flexible case
  and carry the resulting `~10-30%` band explicitly into every integral claim.
- [ ] **T3 - Re-do the §5.9 comparison in one declared frame, on the production path.**
  `plot_spanwise_loads.py` (and the prose it feeds) declares the frame of each curve and converts to a
  common one (or reports both). Decide the honest, reproducible statistic: the rigid BEM at the rated
  point is not the comparable case for a flexible `Fig. 11`, so either use a one-way deformed geometry
  or state the comparison as rigid-vs-flexible and bound it. The stale §5.9 numbers are replaced or
  withdrawn, not widened.
- [ ] **T4 - Land the attribution: write-up, store, closures.**
  `docs/validation_results_v01_v02_v04.md` §5.9/§5.11, `docs/iea15mw_validation_reference.md` `V-09`,
  the article draft if it carries the number; and the store (`docs/validation/gaps.yaml` if a residual
  survives, else a recorded note) plus the `docs/validation_closures.md` re-run trigger. Never widen a
  bound or delete a comparison: a number that moves IS the finding.
- [ ] **T5 - Independent read-only verification, then close.**
  A separate verifier re-derives the frame identity and the quadrature band from the raw CSVs and the
  deck.constants, checks every printed number against the live script, and checks the store. Then
  comment on #15 with the attribution and close it, and update the #18 map.

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
