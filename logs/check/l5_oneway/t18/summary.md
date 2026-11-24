# T18 summary (tip twist equal-weight, flap, same K_eff, aero only; HEAD 84d7019 tree)

| load | F_x kN | F_y kN | tip flap m | edge m | twist eq deg |
|---|---|---|---|---|---|
| T6 production dump mean (127 dumps, 8 rev) | 100.34 | 842.25 | 16.919 | 1.610 | -3.204 |
| production config only: coupling subset, windIO, wall-flow off, u!=0 branch, 36-azimuth mean (l5_t18_subset.py) | 100.31 | 842.29 | 16.920 | 1.609 | -3.205 |
| T3 (AeroDyn polars, wall-flow on, u=0 branch, az 0) | 104.66 | 849.88 | 16.911 | 1.630 | -2.253 |

Chain T3 -> T6 (twist, deg): wall-flow realisation off -0.311 (-2.564); AeroDyn -> windIO polars/station grid -0.690 (-3.254, az 0);
azimuth mean -0.042 (-3.296); u!=0 branch (prebend dropped by `_rebuild_bem_solver`) +0.057 (-3.239);
coupling subset (no web nodes) +0.035 (-3.204). Sum -0.951. Flap: polars -0.10 m, azimuth mean -0.23 m, prebend drop +0.34 m, rest ~0.
Ruled out (numbers in twist.json/acswap.json/loads.json): Cm on/off (both always on; removing Cm flips twist to +3.47 / +2.43 deg, the Cm
total is equal AeroDyn vs windIO to 0.1 %), AC fraction (swap both ways: twist -3.254 vs -3.254, -2.564 vs -2.564), BEMSolver kwargs
(identical dicts), omega/wind/shear/yaw/tilt/precone/hub/tip radius (identical), deformed_twist/radius (false in production).
