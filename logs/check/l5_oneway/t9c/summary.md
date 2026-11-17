## Convergence of the shell measurement (bias-corrected, window half = L/4) vs segment length

| z | chords | L m | psi deg | C_fe/C_ff | GJ/deck | EI_flap/deck | EI_edge/deck |
|---|---|---|---|---|---|---|---|
| 45 | 4 | 19.0 | -6.73 | -0.0637 | 0.951 | 1.065 | 0.958 |
| 45 | 6 | 28.6 | -6.57 | -0.0611 | 0.935 | 1.065 | 0.953 |
| 45 | 8 | 38.2 | -6.37 | -0.0586 | 0.937 | 1.064 | 0.950 |
| 45 | 10 | 47.7 | -6.18 | -0.0565 | 0.940 | 1.062 | 0.949 |
| 45 | 12 | 57.1 | -6.01 | -0.0548 | 0.942 | 1.060 | 0.949 |
| 80 | 4 | 13.2 | -0.02 | -0.0106 | 0.866 | 1.005 | 0.956 |
| 80 | 6 | 19.9 | 1.01 | +0.0042 | 0.831 | 1.003 | 0.939 |
| 80 | 8 | 26.3 | 1.84 | +0.0161 | 0.833 | 1.001 | 0.927 |
| 80 | 10 | 33.0 | 2.30 | +0.0228 | 0.844 | 1.000 | 0.922 |
| 80 | 12 | 39.6 | 2.46 | +0.0254 | 0.853 | 1.001 | 0.921 |
| 95 | 4 | 10.8 | -0.56 | -0.0302 | 0.865 | 1.000 | 0.868 |
| 95 | 6 | 16.3 | 0.02 | -0.0220 | 0.852 | 0.997 | 0.860 |
| 95 | 8 | 21.8 | -0.09 | -0.0239 | 0.860 | 0.998 | 0.864 |
| 95 | 10 | 27.2 | 0.10 | -0.0213 | 0.870 | 0.999 | 0.861 |
| 95 | 12 | 32.6 | 0.40 | -0.0172 | 0.875 | 0.997 | 0.854 |

## Fit window at 8 chords (half width / L)

| z | half/L | psi deg | C_fe/C_ff | GJ/deck | EI_flap/deck | EI_edge/deck |
|---|---|---|---|---|---|---|
| 45 | 0.167 | -6.64 | -0.0622 | 0.940 | 1.066 | 0.955 |
| 45 | 0.250 | -6.37 | -0.0586 | 0.937 | 1.064 | 0.950 |
| 45 | 0.333 | -6.13 | -0.0559 | 0.938 | 1.061 | 0.949 |
| 45 | 0.400 | -5.94 | -0.0541 | 0.942 | 1.059 | 0.950 |
| 80 | 0.167 | 1.30 | +0.0083 | 0.831 | 1.003 | 0.934 |
| 80 | 0.250 | 1.84 | +0.0161 | 0.833 | 1.001 | 0.927 |
| 80 | 0.333 | 1.99 | +0.0186 | 0.841 | 1.000 | 0.927 |
| 80 | 0.400 | 1.98 | +0.0187 | 0.850 | 1.001 | 0.928 |
| 95 | 0.167 | 0.42 | -0.0167 | 0.846 | 0.997 | 0.851 |
| 95 | 0.250 | -0.09 | -0.0239 | 0.860 | 0.998 | 0.864 |
| 95 | 0.333 | -0.03 | -0.0233 | 0.867 | 0.999 | 0.865 |
| 95 | 0.400 | 0.44 | -0.0164 | 0.870 | 0.997 | 0.855 |

(z 80 and z 95 also in study_table.md; raw and bias-corrected columns in study2/study_table.md.)

## Chosen shell setting and criterion
Criterion (tools/l5_t9_arbiter.py `converged_setting`, stated before choosing): per station, the smallest segment length (chords) whose every later
step to the next studied length (4, 6, 8, 10, 12) at window half = L/4 moves the bias-corrected psi by <= 0.5 deg and GJ, EI_flap, EI_edge by <= 2 %.
Result (study2/converged.json): z45 -> 4, z80 -> 8, z95 -> 6 chords; common setting = the maximum, **8 chords, window half width L/4 (2 chords)**.
Sensitivity: with 0.3 deg / 1 % only z45 converges (6 chords), z80 and z95 do not; with 1 deg / 3 %: 4 / 6 / 4.
"Bias corrected" = the same fit applied to a beam with the deck compliance on the SAME shell nodes, shell - (deck_samefit - deck_exact) for psi and ratios,
shell * deck_exact / deck_samefit for stiffnesses: raw shell GJ/EI drift by 1-3 % (taper) with the window, the corrected ones by < 1 %.
Reading band of the shell psi at this setting: about +-0.5 deg (n = 8 -> 10 step, spread over windows L/6 ... L/3 at n = 10:
z45 0.5, z80 0.5, z95 0.3 deg).  Not fully converged: z80 still moves +0.45 (8->10) and +0.16 (10->12) deg, z95 +0.2 deg per step.
Cause of the +-1.5 deg local oscillation of the shell psi along z (wavelength ~6 m, profile at 2.5 m half window; $SCRATCH exploration, not kept):
NOT established (the deck read through the same fit shows none, so it is a property of the shell model, not of the fit).

## Inboard station
Layup stack scan (stack_scan.txt, cells and TE stack thickness every 2.5 m): the TE stack is constant (0.033 m) for z 35-62, thinning from z 65;
z 60 fails the solid builder ("tip spans two stacks"); z 45 builds at 4 chords (z 35.4-54.6, L 19.1 m, inside the constant stack).
Shell cut at z 45 (8 chords, 38 m) is longer than the solid cut and crosses stack steps; that is a property of the production mesh, not of the criterion.

## Solid z45 checks
A (h 0.036 chord) -> B (h 0.018): psi -7.525 -> -7.389 deg (0.14), EI_flap 9.9695 -> 9.9693e9, EI_edge 1.8345 -> 1.8658e10 (1.7 %): bending converged.
Length 4 -> 6 chords (A): psi -7.525 -> -7.545, EI_flap/EI_edge unchanged (< 0.1 %), **GJ 2.904 -> 2.591e8 (-11 %)** and the weight-free 2D St Venant GJ of the
same section is 1.891e8: the solid GJ at z45 is NOT length-converged (warping restraint) and disagrees with its own 2D limit by 37 %
(z80, z95: within 2 %).  The z45 GJ comparison against the solid is therefore not evidence; psi, C_fe/C_ff, EI_flap, EI_edge at z45 are.

## Reproduce
  scripts/aeroenv.sh env PYTHONPATH=$PWD OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python tools/l5_t9_arbiter.py study --z 80 --nchords 4 6 8 10 12 --fracs 0.1667 0.25 0.3333 0.4 --out DIR   (per z; ~5.5 min each)
  ... study-table --out DIR
  ... deck --z 45 80 95 --out stations ; ... shell --z 45 80 95 --out stations (8 chords, 200 s)
  ... solid --z 45 --level B --out SOLIDDIR   (~14 min wall, 3 ccx x ~830 s)  [--cmd build|run|analyze]
  ... table --out stations --solid-dirs ...   ;   python tools/l5_t9_verdict.py logs/check/l5_oneway/t9c
