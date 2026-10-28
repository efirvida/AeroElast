# T8 RESULT (CORRECTED, current Rust build): shell vs BeamDyn p20 mode shapes

**The first run of this file used a stale `_aeroelast.so` (Sep 25) that `PYTHONPATH=<repo>/src` selected in the main checkout; its rated numbers are
superseded** (kept in `stale-binary/`). Rerun with the Oct 26 build (`rerun-*.log`): rated shell 1F/1E/2F = 0.5418/0.7025/1.5090 Hz, which reproduces the stored S-4
(+1.85/+0.44/+0.06 % vs the driver). Parked modes 1-7 moved by <= 0.27 %, the higher ones did not (shell 4F 4.2399 -> 4.3221, 5F 4.5188 -> 4.6340 Hz). The tables are in RESULT-tables.md
(current build); the stale-build narrative and tables are only in `stale-binary/RESULT.md`.

Verdicts on the current build (combined-vector MAC, fit A = all ring nodes):
- P1 HELD: 1F 1E 2F 2E 3F parked bending MAC 0.995 0.992 0.999 0.993 0.995 (rated >= 0.989). P2 HELD for these (change <= 0.004).
- P3 HELD: 1T twist-shape MAC 0.990 (box 0.999) but combined 0.848 parked / 0.810 rated, bending MAC 0.02: the small translation part differs.
  **The 1T frequency by shape is 3.8044 Hz vs 4.0937 Hz, -7.07 % parked, -6.54 % rated** (the group 23 +1.59 % paired BeamDyn 1T with the shell's 3E).
  f ~ sqrt(GJ/Ip): -7 % is about -14 % effective GJ with Ip within a few %, in line with S-1 GJ -10.7 % and S-7 +7.2 %: the shell IS torsionally soft;
  the retraction of the earlier "torsion is the measurement" stands, and T7 (TE closure) is justified by an extractor-free dynamic measure.
- P4: 3E parked 0.963 combined (bending 0.973) at -5.05 %; 4F REFUTED: BeamDyn's single 4F (4.8448 Hz) is split over shell modes at 4.32/4.63/4.70 Hz (bending MAC 0.766/0.924/0.901)
  (best: shell 5F 4.6340 Hz, combined 0.870, -4.35 %); rated pairing is not stable (3E/4E/4F swap). Higher modes are not identifiable by shape at 0.95.
- Rated BeamDyn reference here is the MBC3 collective (hub 3.97 m, precone), shell hub 0: rated gaps are not the S-4 driver gaps (1E -1.31 % vs the driver's +0.44 %).

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
