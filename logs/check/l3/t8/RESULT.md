# T8 RESULT (CORRECTED, current Rust build): shell vs BeamDyn p20 mode shapes

**The first run of this file used a stale `_aeroelast.so` (Sep 25) that `PYTHONPATH=<repo>/src` selected in the main checkout; its rated numbers are
superseded** (kept in `stale-binary/`). Rerun with the Oct 26 build (`rerun-*.log`): rated shell 1F/1E/2F = 0.5418/0.7025/1.5090 Hz, which reproduces the stored S-4
(+1.85/+0.44/+0.06 % vs the driver). Parked is unchanged within 0.2 %. Tables below (RESULT-tables.md) are from the current build; the narrative section
after them is the stale-build text and is retained only for the method.

Verdicts on the current build (combined-vector MAC, fit A = all ring nodes):
- P1 HELD: 1F 1E 2F 2E 3F parked bending MAC 0.995 0.992 0.999 0.993 0.995 (rated >= 0.989). P2 HELD for these (change <= 0.004).
- P3 HELD: 1T twist-shape MAC 0.990 (box 0.999) but combined 0.848 parked / 0.810 rated, bending MAC 0.02: the small translation part differs.
  **The 1T frequency by shape is 3.8044 Hz vs 4.0937 Hz, -7.07 % parked, -6.54 % rated** (the group 23 +1.59 % paired BeamDyn 1T with the shell's 3E).
  f ~ sqrt(GJ/Ip): -7 % is about -14 % effective GJ with Ip within a few %, in line with S-1 GJ -10.7 % and S-7 +7.2 %: the shell IS torsionally soft;
  the retraction of the earlier "torsion is the measurement" stands, and T7 (TE closure) is justified by an extractor-free dynamic measure.
- P4: 3E parked 0.963 combined (bending 0.973) at -5.05 %; 4F REFUTED: BeamDyn's single 4F (4.8448 Hz) is split over shell modes at 4.32/4.52/4.63 Hz
  (best: shell 5F 4.6340 Hz, combined 0.870, -4.35 %); rated pairing is not stable (3E/4E/4F swap). Higher modes are not identifiable by shape at 0.95.
- Rated BeamDyn reference here is the MBC3 collective (hub 3.97 m, precone), shell hub 0: rated gaps are not the S-4 driver gaps (1E -1.31 % vs the driver's +0.44 %).

# T8 RESULT: mode shapes, shell vs BeamDyn p20 (l3-closure T8; S6, S10, S4)

Tool `tools/l3_mode_shapes.py` (steps `bd`, `shell-parked`, `shell-rated`, `compare`), tests `tests/software/test_mode_shapes.py`.
BeamDyn: existing `.lin` (no OpenFAST run), frequency gate vs `modes.json`: max rel diff 6e-12 (parked), 5e-12 (rated), 12 modes each.
Shell: es 0.25, skin mid, `_build_assembly` + `assemble_rotating` (n_modes=12), eigenvectors in `shell-{parked,rated}.npz`.
Wall: bd 18 s; shell parked 151 s; shell rated 721 s; compare 13 s. Only the IEA-15 blade (the 3.4 has no p20 .lin; S10 not covered here).

Check of the shell side: the real S-4 rated test run at this HEAD (`s4-rated-crosscheck.log`, 12 min) gives the same frequencies as the
mode solve here (1F 0.555924, 1E 0.711703, 2F 1.527232 Hz), so the saved eigenvectors are the S-4 model's. **They do not match the stored
S-4 rated numbers** (`logs/check/l3/store-pass1/pyt-s4.txt`, 2026-10-27 07:14: 1F 0.541819, 1E 0.702508, 2F 1.508960 Hz; +1.85/+0.44/+0.06 % vs BeamDyn
there, +4.50/+1.76/+1.27 % now). Parked agrees (1F 0.5187 vs 0.5189). Cause not investigated (no src change to K_G since, `git log`).

## Verdicts (combined-vector MAC, twist from all ring nodes = fit A unless stated; tables below)

- P1 (parked 1F 1E 2F 2E 3F MAC >= 0.98): HELD on bending and combined: bending 0.995 0.992 0.999 0.993 0.995, combined 0.995 0.992 0.999 0.993 0.994.
  Caveat: the twist-only MAC of 2F and 3F is 0.917 / 0.898 (A), 0.704 / 0.804 (box), below 0.98 (twist amplitude of flap modes is tiny). The
  parenthetical "frequencies within 2.4 %" is slightly off: 3F -2.64 %.
- P2 (rated within 0.02 of parked): HELD for 1F 1E 2F 2E 3F (largest change 0.004, 1F combined 0.995 -> 0.991). REFUTED for the higher modes
  if they are meant: 1T 0.829 -> 0.893 (+0.064), 3E 0.938 -> 0.986 (+0.048), 4F 0.449 -> 0.974 (+0.525; the parked 4F row pairs a different
  shell mode, see P4).
- P3 (1T combined MAC weakest of the set, 0.80-0.95): HELD. 1T parked 0.829 (A) / 0.822 (box), rated 0.893 / 0.891; the weakest among
  1F 1E 2F 2E 3F 1T (others >= 0.992). 1T is not below 0.80, so the "not an arbiter of GJ" declaration is not triggered. Mechanism claim
  (shell 1T carries more flap/edge translation): direction held, size small: translation energy share shell 11.4 % (A) / 13.6 % (box) vs
  BeamDyn 7.2 % parked; shell 8.2 % / 8.5 % vs 6.5 % rated. The twist shape itself matches (twist MAC 0.983 A, 0.999 box parked; 0.997 / 0.998
  rated); the bending MAC of the 1T pair is 0.014 / 0.076 (the small translation parts are unrelated shapes).
- P4 (parked 3E and 4F MAC 0.90-0.98): 3E HELD (0.938 A, 0.949 box; bending 0.954). 4F REFUTED: 8x8 pairing 0.449 (A) / 0.501 (box), shell 4F
  at 4.2399 Hz (-12.5 %); with the 12x8 pairing BeamDyn 4F goes to shell 6F (4.6524 Hz, -3.97 %): combined 0.818 (A) / 0.881 (box), bending
  0.918. Rated 4F: 0.974 (A) / 0.959 (box), shell 4F 4.6422 Hz (-4.61 %), so rated would pass.

Contradicting or not predicted:
- Parked, the shell has three flap-like modes between 4.24 and 4.65 Hz (shell 4F, 5F, 6F; BeamDyn has one, 4F at 4.84 Hz, and its next mode is at
  6.79 Hz). The 4F shape is split across them: bending MAC of BeamDyn 4F with shell 8, 9, 10 = 0.573, 0.589, 0.918.
- Shell 3E (7) in the parked row is at -5.30 % from BeamDyn 3E though the shape is matched (0.938); rated -4.81 %.
- The 8x8 Hungarian pairing is the identity (frequency order) in all four tables (parked/rated x A/box).

## Method (details in the tool docstring)

- BeamDyn states: `BD_1 finite element node N ... translational / rotational displacement in X|Y|Z` (derivative order 2, nodes 2..21), i.e. the
  displacement part of the state vector `[q; qdot]` in the blade-root axes; MBC3 collective triplet column 0. Rotations are Wiener-Milenkovic
  parameters, rotation vector to first order. Node z = 117 (1 + xi_GLL)/2 (reproduces the 9 tip nodes of the p60 BeamDyn VTK to 1e-4 m, unit test).
- Axes BeamDyn (x, y, z) -> shell (-y, x, z); phase of a BeamDyn mode: rotation to the dominant real phase of the combined vector, real part.
- Shell: arc-weighted ring mean of (u_x, u_y) (`ring_profile` weights: half of each in-plane element side, webs included, as in S12), twist
  from `rigid_rotation_fit` on all ring nodes (A) or `box_node_mask` nodes (B); linear in z between the bracketing rings; NaN stations 0 (none occurred).
- Combined vector: rotation scaled by chord/2, chord from `tests/IEA-15-240-RWT.yaml` at z/117.
- Pairing on the combined MAC (`pair_modes_by_mac`); 8x8 (spec) and 12x8 (all saved shell modes) both given.
- The shell ring mean is the arc-weighted centroid, BeamDyn's node is its reference axis: a section rotation moves the two by theta x offset, unremoved.

# T8 RESULT: shell vs BeamDyn p20 mode shapes (es 0.25 mid skin; BeamDyn p20 from the existing .lin)

Generated by `tools/l3_mode_shapes.py compare`. Vectors: bending [u_x(z_i); u_y(z_i)], twist [theta_z(z_i)], combined [u_x; u_y; theta c/2] (c from the official yaml chord at z/117). Shell ring means arc-weighted (`ring_profile`), twist from `rigid_rotation_fit` on all ring nodes (A) or box nodes (B), linear in z between rings. Pairing: Hungarian on the combined MAC of shell modes 1..8 x BeamDyn 1..8 (square).

## parked (0 rpm)

BeamDyn p20 modes (family by station energy, E = u_x, F = u_y, T = theta c/2): 1:1F 0.5071 Hz, 2:1E 0.6934 Hz, 3:2F 1.4809 Hz, 4:2E 2.1366 Hz, 5:3F 2.9355 Hz, 6:1T 4.0937 Hz, 7:3E 4.3798 Hz, 8:4F 4.8448 Hz

Shell modes: 1:1F 0.5187 Hz, 2:1E 0.6967 Hz, 3:2F 1.4830 Hz, 4:2E 2.0888 Hz, 5:3F 2.8581 Hz, 6:1T 3.7945 Hz, 7:3E 4.1475 Hz, 8:4F 4.2399 Hz, 9:5F 4.5188 Hz, 10:6F 4.6524 Hz, 11:2T 4.7975 Hz, 12:3T 5.0899 Hz

### parked, A: twist from all ring nodes, pairing 8x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5071 | 1 1F | 0.5187 | +2.29 | 0.995 | 0.989 | 0.995 |
| 2 1E | 0.6934 | 2 1E | 0.6967 | +0.48 | 0.992 | 0.992 | 0.992 |
| 3 2F | 1.4809 | 3 2F | 1.4830 | +0.14 | 0.999 | 0.917 | 0.999 |
| 4 2E | 2.1366 | 4 2E | 2.0888 | -2.24 | 0.993 | 0.989 | 0.993 |
| 5 3F | 2.9355 | 5 3F | 2.8581 | -2.64 | 0.995 | 0.898 | 0.994 |
| 6 1T | 4.0937 | 6 1T | 3.7945 | -7.31 | 0.014 | 0.983 | 0.829 |
| 7 3E | 4.3798 | 7 3E | 4.1475 | -5.30 | 0.954 | 0.556 | 0.938 |
| 8 4F | 4.8448 | 8 4F | 4.2399 | -12.49 | 0.573 | 0.189 | 0.449 |

### parked, A: twist from all ring nodes, pairing 12x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5071 | 1 1F | 0.5187 | +2.29 | 0.995 | 0.989 | 0.995 |
| 2 1E | 0.6934 | 2 1E | 0.6967 | +0.48 | 0.992 | 0.992 | 0.992 |
| 3 2F | 1.4809 | 3 2F | 1.4830 | +0.14 | 0.999 | 0.917 | 0.999 |
| 4 2E | 2.1366 | 4 2E | 2.0888 | -2.24 | 0.993 | 0.989 | 0.993 |
| 5 3F | 2.9355 | 5 3F | 2.8581 | -2.64 | 0.995 | 0.898 | 0.994 |
| 6 1T | 4.0937 | 6 1T | 3.7945 | -7.31 | 0.014 | 0.983 | 0.829 |
| 7 3E | 4.3798 | 7 3E | 4.1475 | -5.30 | 0.954 | 0.556 | 0.938 |
| 8 4F | 4.8448 | 10 6F | 4.6524 | -3.97 | 0.918 | 0.114 | 0.818 |

Full MAC rows (shell modes 1..12 x BeamDyn 1..8), A: twist from all ring nodes; shell modes within 3.7-4.9 Hz and all BeamDyn columns (bend | twist | comb):

- shell 6 1T 3.7945 Hz  bend [0.000 0.272 0.061 0.433 0.199 0.014 0.460 0.176]  twist [0.963 0.977 0.226 0.982 0.673 0.983 0.292 0.965]  comb [0.000 0.024 0.010 0.017 0.014 0.829 0.052 0.074]
- shell 7 3E 4.1475 Hz  bend [0.047 0.504 0.098 0.746 0.102 0.138 0.954 0.069]  twist [0.435 0.325 0.025 0.380 0.598 0.361 0.556 0.216]  comb [0.046 0.494 0.096 0.732 0.102 0.003 0.938 0.064]
- shell 8 4F 4.2399 Hz  bend [0.137 0.191 0.235 0.252 0.378 0.281 0.190 0.573]  twist [0.094 0.068 0.000 0.082 0.007 0.078 0.070 0.189]  comb [0.101 0.142 0.172 0.194 0.280 0.001 0.149 0.449]
- shell 9 5F 4.5188 Hz  bend [0.245 0.133 0.435 0.205 0.558 0.000 0.343 0.589]  twist [0.115 0.104 0.010 0.111 0.040 0.108 0.045 0.186]  comb [0.181 0.097 0.324 0.137 0.407 0.032 0.249 0.472]
- shell 10 6F 4.6524 Hz  bend [0.428 0.009 0.623 0.019 0.775 0.094 0.004 0.918]  twist [0.091 0.067 0.001 0.080 0.027 0.077 0.086 0.114]  comb [0.397 0.008 0.577 0.017 0.714 0.012 0.003 0.818]
- shell 11 2T 4.7975 Hz  bend [0.071 0.022 0.119 0.033 0.162 0.049 0.028 0.180]  twist [0.192 0.151 0.000 0.173 0.123 0.165 0.164 0.213]  comb [0.034 0.009 0.056 0.010 0.080 0.060 0.010 0.063]

BeamDyn columns 1..8 are 1:1F, 2:1E, 3:2F, 4:2E, 5:3F, 6:1T, 7:3E, 8:4F

BeamDyn modes within 3.7-4.9 Hz, column MAC over all shell modes (comb):

- BD 6 1T 4.0937 Hz  comb over shell 1..12 [0.003 0.017 0.000 0.003 0.009 0.829 0.003 0.001 0.032 0.012 0.060 0.075]
- BD 7 3E 4.3798 Hz  comb over shell 1..12 [0.002 0.577 0.019 0.848 0.008 0.052 0.938 0.149 0.249 0.003 0.010 0.009]
- BD 8 4F 4.8448 Hz  comb over shell 1..12 [0.481 0.002 0.710 0.000 0.864 0.074 0.064 0.449 0.472 0.818 0.063 0.120]

### parked, B: twist from box nodes, pairing 8x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5071 | 1 1F | 0.5187 | +2.29 | 0.995 | 0.965 | 0.995 |
| 2 1E | 0.6934 | 2 1E | 0.6967 | +0.48 | 0.992 | 0.991 | 0.992 |
| 3 2F | 1.4809 | 3 2F | 1.4830 | +0.14 | 0.999 | 0.704 | 0.999 |
| 4 2E | 2.1366 | 4 2E | 2.0888 | -2.24 | 0.993 | 0.983 | 0.993 |
| 5 3F | 2.9355 | 5 3F | 2.8581 | -2.64 | 0.995 | 0.804 | 0.991 |
| 6 1T | 4.0937 | 6 1T | 3.7945 | -7.31 | 0.014 | 0.999 | 0.822 |
| 7 3E | 4.3798 | 7 3E | 4.1475 | -5.30 | 0.954 | 0.643 | 0.949 |
| 8 4F | 4.8448 | 8 4F | 4.2399 | -12.49 | 0.573 | 0.468 | 0.501 |

### parked, B: twist from box nodes, pairing 12x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5071 | 1 1F | 0.5187 | +2.29 | 0.995 | 0.965 | 0.995 |
| 2 1E | 0.6934 | 2 1E | 0.6967 | +0.48 | 0.992 | 0.991 | 0.992 |
| 3 2F | 1.4809 | 3 2F | 1.4830 | +0.14 | 0.999 | 0.704 | 0.999 |
| 4 2E | 2.1366 | 4 2E | 2.0888 | -2.24 | 0.993 | 0.983 | 0.993 |
| 5 3F | 2.9355 | 5 3F | 2.8581 | -2.64 | 0.995 | 0.804 | 0.991 |
| 6 1T | 4.0937 | 6 1T | 3.7945 | -7.31 | 0.014 | 0.999 | 0.822 |
| 7 3E | 4.3798 | 7 3E | 4.1475 | -5.30 | 0.954 | 0.643 | 0.949 |
| 8 4F | 4.8448 | 10 6F | 4.6524 | -3.97 | 0.918 | 0.002 | 0.881 |

Full MAC rows (shell modes 1..12 x BeamDyn 1..8), B: twist from box nodes; shell modes within 3.7-4.9 Hz and all BeamDyn columns (bend | twist | comb):

- shell 6 1T 3.7945 Hz  bend [0.000 0.272 0.061 0.433 0.199 0.014 0.460 0.176]  twist [0.976 0.995 0.231 0.997 0.731 0.999 0.297 0.940]  comb [0.000 0.029 0.012 0.023 0.018 0.822 0.064 0.080]
- shell 7 3E 4.1475 Hz  bend [0.047 0.504 0.098 0.746 0.102 0.138 0.954 0.069]  twist [0.624 0.493 0.006 0.561 0.692 0.538 0.643 0.385]  comb [0.047 0.501 0.097 0.743 0.103 0.002 0.949 0.063]
- shell 8 4F 4.2399 Hz  bend [0.137 0.191 0.235 0.252 0.378 0.281 0.190 0.573]  twist [0.583 0.652 0.249 0.623 0.593 0.632 0.108 0.468]  comb [0.131 0.177 0.219 0.220 0.369 0.101 0.183 0.501]
- shell 9 5F 4.5188 Hz  bend [0.245 0.133 0.435 0.205 0.558 0.000 0.343 0.589]  twist [0.310 0.224 0.031 0.266 0.462 0.251 0.434 0.138]  comb [0.239 0.131 0.425 0.202 0.548 0.001 0.338 0.559]
- shell 10 6F 4.6524 Hz  bend [0.428 0.009 0.623 0.019 0.775 0.094 0.004 0.918]  twist [0.070 0.028 0.138 0.047 0.213 0.040 0.303 0.002]  comb [0.411 0.009 0.601 0.019 0.742 0.002 0.003 0.881]
- shell 11 2T 4.7975 Hz  bend [0.071 0.022 0.119 0.033 0.162 0.049 0.028 0.180]  twist [0.110 0.030 0.351 0.064 0.237 0.051 0.618 0.010]  comb [0.070 0.022 0.118 0.034 0.156 0.000 0.025 0.180]

BeamDyn columns 1..8 are 1:1F, 2:1E, 3:2F, 4:2E, 5:3F, 6:1T, 7:3E, 8:4F

BeamDyn modes within 3.7-4.9 Hz, column MAC over all shell modes (comb):

- BD 6 1T 4.0937 Hz  comb over shell 1..12 [0.003 0.017 0.000 0.003 0.011 0.822 0.002 0.101 0.001 0.002 0.000 0.009]
- BD 7 3E 4.3798 Hz  comb over shell 1..12 [0.002 0.577 0.019 0.849 0.008 0.064 0.949 0.183 0.338 0.003 0.025 0.024]
- BD 8 4F 4.8448 Hz  comb over shell 1..12 [0.481 0.002 0.709 0.000 0.858 0.080 0.063 0.501 0.559 0.881 0.180 0.174]

## rated (7.56 rpm)

BeamDyn p20 modes (family by station energy, E = u_x, F = u_y, T = theta c/2): 1:1F 0.5320 Hz, 2:1E 0.7118 Hz, 3:2F 1.5081 Hz, 4:2E 2.1668 Hz, 5:3F 2.9604 Hz, 6:1T 4.0676 Hz, 7:3E 4.4059 Hz, 8:4F 4.8667 Hz

Shell modes: 1:1F 0.5559 Hz, 2:1E 0.7117 Hz, 3:2F 1.5272 Hz, 4:2E 2.1183 Hz, 5:3F 2.9063 Hz, 6:1T 3.9084 Hz, 7:3E 4.1939 Hz, 8:4F 4.6422 Hz, 9:2T 5.5116 Hz, 10:3T 6.1560 Hz, 11:5F 6.2937 Hz, 12:4T 6.6705 Hz

### rated, A: twist from all ring nodes, pairing 8x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5320 | 1 1F | 0.5559 | +4.50 | 0.991 | 0.993 | 0.991 |
| 2 1E | 0.7118 | 2 1E | 0.7117 | -0.02 | 0.988 | 0.996 | 0.988 |
| 3 2F | 1.5081 | 3 2F | 1.5272 | +1.27 | 0.999 | 0.960 | 0.999 |
| 4 2E | 2.1668 | 4 2E | 2.1183 | -2.24 | 0.993 | 0.994 | 0.993 |
| 5 3F | 2.9604 | 5 3F | 2.9063 | -1.83 | 0.995 | 0.940 | 0.995 |
| 6 1T | 4.0676 | 6 1T | 3.9084 | -3.91 | 0.076 | 0.997 | 0.893 |
| 7 3E | 4.4059 | 7 3E | 4.1939 | -4.81 | 0.989 | 0.835 | 0.986 |
| 8 4F | 4.8667 | 8 4F | 4.6422 | -4.61 | 0.979 | 0.401 | 0.974 |

### rated, A: twist from all ring nodes, pairing 12x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5320 | 1 1F | 0.5559 | +4.50 | 0.991 | 0.993 | 0.991 |
| 2 1E | 0.7118 | 2 1E | 0.7117 | -0.02 | 0.988 | 0.996 | 0.988 |
| 3 2F | 1.5081 | 3 2F | 1.5272 | +1.27 | 0.999 | 0.960 | 0.999 |
| 4 2E | 2.1668 | 4 2E | 2.1183 | -2.24 | 0.993 | 0.994 | 0.993 |
| 5 3F | 2.9604 | 5 3F | 2.9063 | -1.83 | 0.995 | 0.940 | 0.995 |
| 6 1T | 4.0676 | 6 1T | 3.9084 | -3.91 | 0.076 | 0.997 | 0.893 |
| 7 3E | 4.4059 | 7 3E | 4.1939 | -4.81 | 0.989 | 0.835 | 0.986 |
| 8 4F | 4.8667 | 8 4F | 4.6422 | -4.61 | 0.979 | 0.401 | 0.974 |

Full MAC rows (shell modes 1..12 x BeamDyn 1..8), A: twist from all ring nodes; shell modes within 3.7-4.9 Hz and all BeamDyn columns (bend | twist | comb):

- shell 6 1T 3.9084 Hz  bend [0.046 0.385 0.000 0.586 0.032 0.076 0.517 0.008]  twist [0.966 0.997 0.276 0.992 0.759 0.997 0.235 0.931]  comb [0.005 0.025 0.000 0.018 0.000 0.893 0.045 0.022]
- shell 7 3E 4.1939 Hz  bend [0.024 0.542 0.056 0.798 0.054 0.178 0.989 0.026]  twist [0.587 0.417 0.056 0.505 0.671 0.474 0.835 0.325]  comb [0.024 0.538 0.056 0.791 0.054 0.007 0.986 0.024]
- shell 8 4F 4.6422 Hz  bend [0.391 0.009 0.613 0.016 0.807 0.042 0.001 0.979]  twist [0.166 0.198 0.178 0.182 0.025 0.190 0.000 0.401]  comb [0.384 0.010 0.605 0.018 0.787 0.001 0.001 0.974]

BeamDyn columns 1..8 are 1:1F, 2:1E, 3:2F, 4:2E, 5:3F, 6:1T, 7:3E, 8:4F

BeamDyn modes within 3.7-4.9 Hz, column MAC over all shell modes (comb):

- BD 6 1T 4.0676 Hz  comb over shell 1..12 [0.006 0.014 0.001 0.002 0.003 0.893 0.007 0.001 0.204 0.047 0.170 0.061]
- BD 7 3E 4.4059 Hz  comb over shell 1..12 [0.001 0.573 0.017 0.849 0.008 0.045 0.986 0.001 0.013 0.047 0.161 0.001]
- BD 8 4F 4.8667 Hz  comb over shell 1..12 [0.476 0.004 0.712 0.001 0.868 0.022 0.024 0.974 0.000 0.326 0.336 0.199]

### rated, B: twist from box nodes, pairing 8x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5320 | 1 1F | 0.5559 | +4.50 | 0.991 | 0.964 | 0.991 |
| 2 1E | 0.7118 | 2 1E | 0.7117 | -0.02 | 0.988 | 0.992 | 0.988 |
| 3 2F | 1.5081 | 3 2F | 1.5272 | +1.27 | 0.999 | 0.755 | 0.999 |
| 4 2E | 2.1668 | 4 2E | 2.1183 | -2.24 | 0.993 | 0.984 | 0.993 |
| 5 3F | 2.9604 | 5 3F | 2.9063 | -1.83 | 0.995 | 0.835 | 0.993 |
| 6 1T | 4.0676 | 6 1T | 3.9084 | -3.91 | 0.076 | 0.998 | 0.891 |
| 7 3E | 4.4059 | 7 3E | 4.1939 | -4.81 | 0.989 | 0.706 | 0.987 |
| 8 4F | 4.8667 | 8 4F | 4.6422 | -4.61 | 0.979 | 0.002 | 0.959 |

### rated, B: twist from box nodes, pairing 12x8 on the combined MAC

| BD mode | BD Hz | shell mode | shell Hz | gap % | MAC bend | MAC twist | MAC comb |
|---|---|---|---|---|---|---|---|
| 1 1F | 0.5320 | 1 1F | 0.5559 | +4.50 | 0.991 | 0.964 | 0.991 |
| 2 1E | 0.7118 | 2 1E | 0.7117 | -0.02 | 0.988 | 0.992 | 0.988 |
| 3 2F | 1.5081 | 3 2F | 1.5272 | +1.27 | 0.999 | 0.755 | 0.999 |
| 4 2E | 2.1668 | 4 2E | 2.1183 | -2.24 | 0.993 | 0.984 | 0.993 |
| 5 3F | 2.9604 | 5 3F | 2.9063 | -1.83 | 0.995 | 0.835 | 0.993 |
| 6 1T | 4.0676 | 6 1T | 3.9084 | -3.91 | 0.076 | 0.998 | 0.891 |
| 7 3E | 4.4059 | 7 3E | 4.1939 | -4.81 | 0.989 | 0.706 | 0.987 |
| 8 4F | 4.8667 | 8 4F | 4.6422 | -4.61 | 0.979 | 0.002 | 0.959 |

Full MAC rows (shell modes 1..12 x BeamDyn 1..8), B: twist from box nodes; shell modes within 3.7-4.9 Hz and all BeamDyn columns (bend | twist | comb):

- shell 6 1T 3.9084 Hz  bend [0.046 0.385 0.000 0.586 0.032 0.076 0.517 0.008]  twist [0.973 0.997 0.257 0.995 0.785 0.998 0.252 0.921]  comb [0.005 0.026 0.000 0.019 0.000 0.891 0.046 0.022]
- shell 7 3E 4.1939 Hz  bend [0.024 0.542 0.056 0.798 0.054 0.178 0.989 0.026]  twist [0.509 0.366 0.040 0.440 0.580 0.415 0.706 0.272]  comb [0.024 0.539 0.056 0.792 0.053 0.007 0.987 0.025]
- shell 8 4F 4.6422 Hz  bend [0.391 0.009 0.613 0.016 0.807 0.042 0.001 0.979]  twist [0.044 0.017 0.119 0.029 0.172 0.024 0.259 0.002]  comb [0.383 0.009 0.603 0.017 0.789 0.000 0.001 0.959]

BeamDyn columns 1..8 are 1:1F, 2:1E, 3:2F, 4:2E, 5:3F, 6:1T, 7:3E, 8:4F

BeamDyn modes within 3.7-4.9 Hz, column MAC over all shell modes (comb):

- BD 6 1T 4.0676 Hz  comb over shell 1..12 [0.006 0.015 0.001 0.002 0.004 0.891 0.007 0.000 0.014 0.022 0.122 0.044]
- BD 7 3E 4.4059 Hz  comb over shell 1..12 [0.001 0.573 0.017 0.850 0.008 0.046 0.987 0.001 0.085 0.100 0.175 0.003]
- BD 8 4F 4.8667 Hz  comb over shell 1..12 [0.476 0.004 0.712 0.001 0.865 0.022 0.025 0.959 0.050 0.518 0.386 0.402]

