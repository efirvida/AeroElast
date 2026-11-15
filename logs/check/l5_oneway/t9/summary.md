# L5 T9 flap-edge cross-coupling arbiter

Frame: mesh x = edge, mesh y = flap, not rotated by twist; C = d(kappa_x, kappa_y)/d(Mx, My) [1/(N m^2)];
psi = angle of the largest-compliance (flap-like) principal axis from mesh +x toward +y.
Window (a): +-3.5 m (shell slab fit, order 2; BeamDyn linear slope). Tip couple 10000 N m.

| z [m] | view | C_ff | C_ee | C_fe | C_ef | C_fe/C_ff | psi [deg] |
|---|---|---|---|---|---|---|---|
| 58.50 | shell (a) | 1.9377e-10 | 6.9960e-11 | 1.075e-12 | 8.462e-13 | +0.0056 | +0.44 |
| 58.50 | BeamDyn (a) | 1.8880e-10 | 7.8779e-11 | -8.012e-12 | 3.110e-13 | -0.0424 | -2.00 |
| 58.50 | deck (b) R=Rz(-tw) [BeamDyn sense] | 2.0355e-10 | 7.4562e-11 | -3.750e-12 | -3.750e-12 | -0.0184 | -1.66 |
| 58.50 | deck (b) R=Rz(+tw) [mesh sense] | 2.0353e-10 | 7.4582e-11 | 4.087e-12 | 4.087e-12 | +0.0201 | +1.81 |
| 58.50 | reference (c) yaml thin-wall | 1.5759e-10 | 6.6120e-11 | 1.588e-12 | 1.588e-12 | +0.0101 | +0.99 |
| 58.50 | reference (c2) mesh-cut mid | 1.8727e-10 | 6.4607e-11 | 7.516e-13 | 7.516e-13 | +0.0040 | +0.35 |
| 58.50 | reference (c2) mesh-cut OML | 1.5669e-10 | 6.2736e-11 | 5.010e-13 | 5.010e-13 | +0.0032 | +0.31 |
| 70.20 | shell (a) | 3.6109e-10 | 1.4023e-10 | 1.299e-12 | 7.508e-13 | +0.0036 | +0.27 |
| 70.20 | BeamDyn (a) | 1.4095e-10 | 6.7266e-11 | -5.574e-12 | -4.743e-12 | -0.0395 | -3.98 |
| 70.20 | deck (b) R=Rz(-tw) [BeamDyn sense] | 3.7607e-10 | 1.4777e-10 | -1.088e-12 | -1.088e-12 | -0.0029 | -0.27 |
| 70.20 | deck (b) R=Rz(+tw) [mesh sense] | 3.7602e-10 | 1.4782e-10 | 3.717e-12 | 3.717e-12 | +0.0099 | +0.93 |
| 70.20 | reference (c) yaml thin-wall | 2.8519e-10 | 1.1813e-10 | 1.511e-12 | 1.511e-12 | +0.0053 | +0.52 |
| 70.20 | reference (c2) mesh-cut mid | 3.4907e-10 | 1.1527e-10 | -3.875e-13 | -3.875e-13 | -0.0011 | -0.09 |
| 70.20 | reference (c2) mesh-cut OML | 2.8454e-10 | 1.1279e-10 | -8.726e-14 | -8.726e-14 | -0.0003 | -0.03 |
| 81.90 | shell (a) | 7.4929e-10 | 4.2141e-10 | -1.860e-11 | -7.150e-13 | -0.0248 | -1.69 |
| 81.90 | BeamDyn (a) | 1.2854e-09 | 5.7521e-10 | 2.545e-11 | 5.348e-12 | +0.0198 | +1.24 |
| 81.90 | deck (b) R=Rz(-tw) [BeamDyn sense] | 7.4372e-10 | 4.0488e-10 | 7.971e-12 | 7.971e-12 | +0.0107 | +1.35 |
| 81.90 | deck (b) R=Rz(+tw) [mesh sense] | 7.4388e-10 | 4.0473e-10 | 1.300e-12 | 1.300e-12 | +0.0017 | +0.22 |
| 81.90 | reference (c) yaml thin-wall | 5.6317e-10 | 3.7644e-10 | 3.310e-12 | 3.310e-12 | +0.0059 | +1.02 |
| 81.90 | reference (c2) mesh-cut mid | 7.2661e-10 | 3.6345e-10 | -2.711e-12 | -2.711e-12 | -0.0037 | -0.43 |
| 81.90 | reference (c2) mesh-cut OML | 5.8279e-10 | 3.5759e-10 | -1.402e-13 | -1.402e-13 | -0.0002 | -0.04 |
| 93.60 | shell (a) | 2.1430e-09 | 1.2727e-09 | -1.235e-10 | -9.591e-12 | -0.0576 | -4.35 |
| 93.60 | BeamDyn (a) | 1.5988e-09 | 1.0560e-09 | 2.971e-11 | 3.195e-11 | +0.0186 | +3.24 |
| 93.60 | deck (b) R=Rz(-tw) [BeamDyn sense] | 2.1748e-09 | 1.0712e-09 | 5.436e-11 | 5.436e-11 | +0.0250 | +2.81 |
| 93.60 | deck (b) R=Rz(+tw) [mesh sense] | 2.1765e-09 | 1.0695e-09 | -2.976e-11 | -2.976e-11 | -0.0137 | -1.54 |
| 93.60 | reference (c) yaml thin-wall | 1.8304e-09 | 1.0622e-09 | -1.452e-11 | -1.452e-11 | -0.0079 | -1.08 |
| 93.60 | reference (c2) mesh-cut mid | 2.0831e-09 | 1.0096e-09 | -3.140e-11 | -3.140e-11 | -0.0151 | -1.67 |
| 93.60 | reference (c2) mesh-cut OML | 1.8090e-09 | 9.9600e-10 | -1.969e-11 | -1.969e-11 | -0.0109 | -1.39 |
| 105.30 | shell (a) (UNRESOLVED: window spread / plane-fit distortion) | 8.9779e-09 | 4.9230e-09 | 1.040e-10 | 2.404e-09 | +0.0116 | +15.87 |
| 105.30 | BeamDyn (a) | 9.2876e-09 | 2.4018e-09 | 2.595e-10 | 3.105e-10 | +0.0279 | +2.37 |
| 105.30 | deck (b) R=Rz(-tw) [BeamDyn sense] | 8.9784e-09 | 3.0301e-09 | 2.684e-10 | 2.684e-10 | +0.0299 | +2.58 |
| 105.30 | deck (b) R=Rz(+tw) [mesh sense] | 8.9829e-09 | 3.0256e-09 | -2.071e-10 | -2.071e-10 | -0.0231 | -1.99 |
| 105.30 | reference (c) yaml thin-wall | 8.2470e-09 | 3.6033e-09 | -1.616e-10 | -1.616e-10 | -0.0196 | -1.99 |
| 105.30 | reference (c2) mesh-cut mid | 8.7068e-09 | 3.3204e-09 | -2.043e-10 | -2.043e-10 | -0.0235 | -2.17 |
| 105.30 | reference (c2) mesh-cut OML | 8.1486e-09 | 3.2954e-09 | -1.785e-10 | -1.785e-10 | -0.0219 | -2.10 |
| 111.15 | shell (a) (UNRESOLVED: window spread / plane-fit distortion) | 1.0196e-06 | 4.0990e-09 | -1.015e-08 | -4.685e-10 | -0.0100 | -0.30 |
| 111.15 | BeamDyn (a) | 4.2983e-08 | 9.6996e-09 | 1.005e-09 | 1.091e-09 | +0.0234 | +1.80 |
| 111.15 | deck (b) R=Rz(-tw) [BeamDyn sense] | 2.8027e-08 | 5.9769e-09 | 7.435e-10 | 7.435e-10 | +0.0265 | +1.93 |
| 111.15 | deck (b) R=Rz(+tw) [mesh sense] | 2.8028e-08 | 5.9755e-09 | -7.213e-10 | -7.213e-10 | -0.0257 | -1.87 |
| 111.15 | reference (c) yaml thin-wall | 2.6472e-08 | 7.9762e-09 | -7.036e-10 | -7.036e-10 | -0.0266 | -2.18 |
| 111.15 | reference (c2) mesh-cut mid | 2.6782e-08 | 7.1637e-09 | -7.641e-10 | -7.641e-10 | -0.0285 | -2.23 |
| 111.15 | reference (c2) mesh-cut OML | 2.5697e-08 | 7.1268e-09 | -7.180e-10 | -7.180e-10 | -0.0279 | -2.21 |

## Shell (a) fit quality

| z [m] | spread C_ff over windows | spread C_ee | plane-fit distortion (rms/rotation signal) | resolved |
|---|---|---|---|---|
| 58.50 | 0.005 | 0.005 | 0.008 | True |
| 70.20 | 0.006 | 0.024 | 0.011 | True |
| 81.90 | 0.015 | 0.010 | 0.019 | True |
| 93.60 | 0.033 | 0.016 | 0.057 | True |
| 105.30 | 0.298 | 0.007 | 0.860 | False |
| 111.15 | 0.468 | 1.600 | 1.173 | False |

## Chord/twist angle per station

| z [m] | mesh chord angle [deg] | windIO twist [deg] | deck kp twist [deg] | yaml-frame principal angle [deg] | alignment: reflect / 2nd-best residual ratio |
|---|---|---|---|---|---|
| 58.50 | +1.70 | +1.69 | +1.69 | +1.11 | True / 1.3x |
| 70.20 | +0.46 | +0.56 | +0.56 | +0.30 | True / 1.5x |
| 81.90 | -0.36 | -0.48 | -0.49 | -1.07 | True / 1.6x |
| 93.60 | -1.92 | -1.93 | -1.92 | -0.78 | True / 2.0x |
| 105.30 | -2.10 | -2.10 | -2.09 | -0.05 | True / 1.5x |
| 111.15 | -1.82 | -1.81 | -1.80 | +0.41 | True / 1.3x |

## Verdicts (5 % on C_ff and C_ee, 1 deg on psi; per-station detail in summary.json)

| comparison | stations | all ok | worst dC_ff | worst dC_ee | worst dpsi [deg] |
|---|---|---|---|---|---|
| BD(a) vs deck(b) bd_sense_minus | 6 | False | +0.728 | +0.623 | 3.71 |
| BD(a) vs deck(b) harness_sense_plus | 6 | False | +0.728 | +0.623 | 4.92 |
| shell (a) vs (c) | 4 | False | +0.330 | +0.198 | 3.26 |
| BeamDyn (a) vs (c) | 6 | False | +1.282 | +0.528 | 4.50 |
| deck (b) BeamDyn sense vs (c) | 6 | False | +0.321 | +0.251 | 4.57 |
| deck (b) mesh sense vs (c) | 6 | False | +0.321 | +0.251 | 0.82 |
| shell (a) vs deck (b) mesh sense | 4 | False | +0.048 | +0.190 | 2.81 |
| shell (a) vs deck (b) BeamDyn sense | 4 | False | +0.048 | +0.188 | 7.16 |

(c) is a thin-wall estimate with its own spread: mesh-cut/yaml - 1 on C_ff (mid, OML) = +0.19/-0.01, +0.22/-0.00, +0.29/+0.03, +0.14/-0.01, +0.06/-0.01, +0.01/-0.03; on C_ee = -0.02/-0.05, -0.02/-0.05, -0.03/-0.05, -0.05/-0.06, -0.08/-0.09, -0.10/-0.11

## Harness check: BD on a uniform blade (K of the station, constant twist) vs view (b)

| z [m] | kp twist [deg] | BD C_fe/C_ff | psi_BD | max rel diff (b) Rz(-tw) | max rel diff (b) Rz(+tw) |
|---|---|---|---|---|---|
| 58.50 | +1.69 | -0.0184 | -1.62 | 4.9e-04 | 2.1e+00 |
| 70.20 | +0.56 | -0.0029 | -0.25 | 4.9e-04 | 4.4e+00 |
| 81.90 | -0.49 | +0.0107 | +1.17 | 5.0e-04 | 8.4e-01 |
| 93.60 | -1.92 | +0.0251 | +2.49 | 5.0e-04 | 1.5e+00 |
| 105.30 | -2.09 | +0.0300 | +2.36 | 5.0e-04 | 1.8e+00 |
| 111.15 | -1.80 | +0.0266 | +1.83 | 5.0e-04 | 2.0e+00 |

## 1D prediction of the L8 tip deflections [m] (parked, AeroDyn rows, deck/yaml compliance)

| compliance | flap-only load: u_edge / u_flap | edge-only load: u_edge / u_flap |
|---|---|---|
| bd_sense_minus | +0.702 / +17.881 | +1.101 / +0.082 |
| harness_sense_plus | -0.155 / +18.103 | +1.075 / -0.020 |
| reference_c_yaml_0.25-0.97R_deck_plus_outside | +0.042 / +15.341 | +1.048 / +0.004 |
| L8 shell (3-D, magnitudes) | 0.048 / - | - / 0.026 |
| L8 BeamDyn (magnitudes) | 0.591 / - | - / 0.081 |
