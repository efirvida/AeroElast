# L5 T10 sign chain: measured links

- mesh rings (good) 30/50: LE on +x True, suction on +y True, handedness(LE, suction, span) +1, tip ring centre y -3.98 m (prebend upwind = -y)
- chord psi vs BlTwist: max |psi - BlTwist| 0.37 deg, max |psi + BlTwist| 15.31 deg
- alpha(mesh chord, BEM inflow) - alpha_BEM: max 39.60 deg; with the chord at -BlTwist: 0.046 deg; production _twist_mesh_to_bem -1

| z [m] | psi mesh [deg] | BlTwist [deg] | alpha BEM | alpha mesh chord | alpha chord at -tw | d alpha / +1 deg about +z |
|---:|---:|---:|---:|---:|---:|---:|
| -0.0 | +21.22 | +15.59 | +53.38 | -94.69 | +57.87 | -1.000 |
| 9.5 | +14.26 | +14.26 | +26.71 | +55.34 | +26.83 | +1.000 |
| 19.2 | +10.28 | +10.37 | +11.91 | +32.51 | +11.86 | +1.000 |
| 28.6 | +6.83 | +7.19 | +8.09 | +22.12 | +8.11 | +1.000 |
| 38.2 | +4.79 | +4.80 | +7.20 | +16.79 | +7.20 | +1.000 |
| 47.7 | +3.08 | +3.08 | +6.82 | +12.99 | +6.83 | +1.000 |
| 57.3 | +1.85 | +1.83 | +6.58 | +10.26 | +6.58 | +1.000 |
| 66.9 | +0.70 | +0.84 | +6.44 | +7.99 | +6.45 | +1.000 |
| 76.4 | +0.03 | +0.04 | +6.38 | +6.46 | +6.39 | +1.000 |
| 86.0 | -1.05 | -1.00 | +6.67 | +4.64 | +6.68 | +1.000 |
| 95.5 | -2.11 | -2.09 | +7.11 | +2.93 | +7.13 | +1.000 |
| 105.0 | -2.10 | -2.10 | +6.82 | +2.64 | +6.84 | +1.000 |
| 114.6 | -1.49 | -1.51 | +5.08 | +2.10 | +5.10 | +1.000 |

- loads bem_deck_tilt-6_precone-4: F = (+104656, +849884, +0) N, torque about +y +7.0077e+06 N m, spin +y -> power/blade +5.548e+06 W
- loads bem_tilt0_precone0: F = (+104160, +851603, +0) N, torque about +y +6.9618e+06 N m, spin +y -> power/blade +5.512e+06 W
- AD Fn/Ft physical_minus_theta: Fx +97428 N, Fy +825041 N, torque 6.3837e+06 N m, 3-blade power 15.16 MW
- AD Fn/Ft apply_ad_loads_plus_theta: Fx +76795 N, Fy +827929 N, torque 7.1447e+06 N m, 3-blade power 16.97 MW
- v1 NuMAD (rotorspin 1): psi = -twist True; v2 gmsh/windIO frames: psi = +twist True
