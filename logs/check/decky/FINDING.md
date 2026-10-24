# DECKY: does the deck-era root TE reinforcement move S-2 LC3? (T26)

- Metric: S-2 shell tip-node dominant component (LC1/LC2/LC4 uy, LC3 ux), `s2_results` fixture of
  `tests/test_iea15mw_s2_static_prescribed.py`. Reference: OpenFAST v5 BeamDyn p60 (LC1 -2.231483, LC2 82.273278,
  LC4 9.062631, LC3 ux 34.469441 m). Tolerance: 5 % (`rel_tol_dominant`); LC3 is printed, not asserted.
- Revision: both captures are ad5f3ac (T14m = clean tree; DECKY = worktree `wt/input-fix-DECKY`, HEAD ad5f3ac, only
  `tests/IEA-15-240-RWT.yaml` modified, 262+/467-). Captures read, nothing rerun. Scripts: `table.py`, `te.py`, `other.py`.

## Numbers (`table.txt`; err = shell/BD - 1; delta = DECKY - T14m)

| case | comp | DECKY [m] | T14m [m] | BeamDyn | err DECKY | err T14m | delta [m] | delta [pp] |
|---|---|---|---|---|---|---|---|---|
| LC1 | uy | -2.101038 | -2.107695 | -2.231483 | -5.846 % | -5.547 % | +0.006658 | -0.298 |
| LC2 | uy | 80.970808 | 80.194963 | 82.273278 | -1.583 % | -2.526 % | +0.775845 | +0.943 |
| LC3 | ux | 38.514987 | 38.352309 | 34.469441 | +11.737 % | +11.265 % | +0.162678 | +0.472 |
| LC4 | uy | 8.565501 | 8.566503 | 9.062631 | -5.485 % | -5.474 % | -0.001002 | -0.011 |

Tip-ring mean (LC3 ux): DECKY 38.339, T14m 37.924 (+1.1 %). No other input variant is in these pickles.

## YAML difference (bounded; `te.txt`, `yaml.diff`)

TE_reinforcement_SS and _PS (glass_uni, width 0.75, same material and width as the official yaml): the deck-era
thickness grid adds 2.98 mm at r/R 0.061 and 8.39 mm at r/R 0.082 (official: 0 there); from r/R 0.102 on the
thickness is identical (13.97 mm, 30 mm plateau). The width is 0 up to r/R 0.082 and ramps to 0.75 at 0.102 in
both files, so the added material acts only in the ramp 0.082-0.102 (about 10-12 m of the span), at partial width.
The root TE change is therefore tiny.

## Confounders (DECKY is not a pure root-TE probe)

Git: DECKY HEAD = ad5f3ac, same as T14m, so no code difference. The whole yaml is swapped, though, and it also
differs in: (1) the tip chord grid: official has points at r/R 0.985/0.99/0.995 (1.78/1.708/1.472 m), deck-era
goes linearly 2.053 (0.939) -> 1.82 (0.98) -> 0.5 (1.0) (chord at 0.99: 1.708 vs 1.16 m, at 0.995: 1.472 vs 0.83);
(2) web `offset_y_pa` grid resampled (2 -> many points); (3) the width/arc grids of the other layers resampled
(same values); (4) `elastic_properties_mb` 6x6 tables (not read by the shell mesh path; S-2 beam reference comes
from the ElastoDyn/BeamDyn files); (5) Gelcoat Xc/Xt units (strength only). Twist, reference axis, airfoil
positions and all other thickness grids are identical.

## Verdict

- The root TE reinforcement does **not** explain LC3. LC3 moves +0.47 pp (+0.42 % in ux), i.e. the deck-era shell is
  more flexible, +11.74 % against +11.27 %. Added root TE material would stiffen (smaller ux), so the sign is
  opposite to what the root TE could cause, and the size is 0.47 pp against a ~10-11 pp gap.
- LC2 moves +0.94 pp (flapwise, insensitive to TE) and LC1/LC4 almost nothing: the shift is not specific to
  edgewise stiffness. It is most plausibly the tip chord grid (confounder 1), which changes the tip geometry and
  the tip-node response; this is an inference, not measured (it would need a yaml with only one change).
- Baseline note: T14m LC3 is +11.265 %, not +10.26 %; +10.26 % is the T15m/a35bba3 row (T15 web corners).
  The gap that remains after T15 (+10.26 %) is not addressed by the deck-era root TE either.
- Not tested: root TE alone, tip chord alone. If wanted, two single-change yamls on ad5f3ac (SLURM, ~S-2 case time).
