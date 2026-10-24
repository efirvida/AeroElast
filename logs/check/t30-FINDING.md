# T30 — TE band extent ["d","e"] (T24a) measured vs d5d0fd2

Commit 5edd82c (fix/input-T30, off ae483fa; blade code d5d0fd2 + one-line extent change). es 0.25 mid-surface, official yaml.
SLURM 11623402 (S-1, S-2, V-03; `capture-A.pkl`, `pytest-A.log`) and 11623403 (S-4, V-02; `pytest-B.log`, `capture-B.pkl`).
Baseline = T15m capture (11623286) / S-4, V-02 values of the task. Prediction is in `measure.sbatch`.

| metric | d5d0fd2 | T30 | delta pp | bound | verdict |
|---|---|---|---|---|---|
| V-03 mass | 67,150 kg (-1.14 %) | 66,783 kg | -0.55 % mass | 5 % | pass; as predicted |
| S-1 EA mean abs | 6.68 | 6.32 | -0.36 | 5 | red before and after |
| S-1 EI_edge mean abs (signed) | 8.01 (-1.89) | 8.14 (-5.08) | +0.13 (-3.2 signed) | 5 | red before and after; signed -5 % as predicted |
| S-1 EI_flap mean abs | 4.53 | 4.53 | 0.00 | 5 | pass |
| S-1 GJ mean abs (signed) | 10.18 (-9.78) | **100.36 (+69.7)** | **+90.2** | 5 | red before; now absurd |
| S-1 mass mean abs | 3.54 | 3.28 | -0.26 | 5 | pass |
| S-2 LC1 | -5.26 | -5.56 | -0.30 | 5 | red before and after |
| S-2 LC2 | -2.54 | -2.53 | +0.01 | 5 | pass |
| S-2 LC3 (edge) | +10.26 | +15.67 | +5.4 | none asserted (ratio test passes) | worse |
| S-2 LC4 | -5.48 | -5.48 | 0.00 | 5 | red before and after |
| S-4 parked 1F / 1E / 2F | +2.33 / +0.56 / +0.20 | +2.26 / -3.09 / **-7.53** | -0.07 / -3.65 / **-7.73** | 5 | 2F newly red |
| S-4 rated 1F / 1E / 2F | +1.85 / +0.44 / +0.06 | +1.77 / -3.00 / **-7.39** | -0.08 / -3.44 / **-7.45** | 5 | 2F newly red |
| V-02 1F / 1E | 7.09 / 8.84 | 7.15 / 4.90 | +0.06 / -3.94 | 5 | 1F red before and after; 1E now passes |

## Against the prediction
Matched: mass -0.55 %, EI_edge signed -1.9 -> -5.1 %, EA down, EI_flap and 1F unchanged, LC3 deflection up, 1E down.
**Not predicted, unexplained:** (1) S-1 GJ mean abs 10 -> 100 % (signed -9.8 -> +69.7); (2) S-4 2nd flap -7.5 % (moved 7.7 pp
while 1F did not move). Removing glass from a thin flat face should not do either. Possible causes: GJ post-processing breaking when
TE_FLAT loses its glass (thin 22-ply flat face with a different laminate), or a mesh/ply-grouping artifact. Not investigated here; no
diagnostic was run (per-station GJ, TE_FLAT section) — that is the next step before any judgement of the physics.

## Recommendation: RECORD, do not merge
The change is a fidelity gain against the yaml (366.7 kg of glass_uni the yaml does not place; mass and EI_edge moved as predicted), so the
defect stays documented as T24a and must not be kept silently. But it worsens anchors: S-4 2F newly outside 5 % (parked and rated),
S-1 GJ collapses to 100 % error, S-2 LC3 +5.4 pp; the offsetting gain is S-4 1E (+0.5 -> -3.1 %, in band) and V-02 1E. Anchors are not moved
to fit. Before re-proposing: diagnose the GJ and 2F jumps (they look like an artifact, not physics); if they are an artifact the
candidate may still be merged as the fidelity fix, and then S-1 EI_edge/EA re-judged.

## Tests (observed)
- RED: `tests/test_blade_te_flat_glass.py` on the unfixed extent: 1 failed (`_band_stacks(stacks, "_TE_FLAT") == []` failed, TE_reinforcement_* plies on `*_TE_FLAT`).
- GREEN: after the extent change `["d","te"]` -> `["d","e"]` (4 lines, yaml_to_blade.py:449-460): 1 passed. `ruff check` on both files: all checks passed.
- Neighbouring `tests/test_blade_*.py`: **not completed** (login-node runs hit the tool timeout; the first run showed `test_blade_bay_length` 2 passed, then
  `test_blade_ccx_parity` with an F, unconfirmed whether pre-existing). Needs a SLURM run.
