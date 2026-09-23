# Torque-deformation frequency overlay

Source: `frontiersin_results_corotational_100s/yaw_*/fluid/bem_report.csv`, window `40-100 s`.
Rotation: `1P = 0.12530 Hz`; for a three-blade ideal rotor, `3P`, `6P`, `9P`, ... survive exact phase summation, while non-multiples of 3 are local/single-blade harmonics unless symmetry is broken.
Frequency-bin resolution is approximately `0.01667 Hz`; zero-padding is used only to make peak locations easier to read in the figure.

## Dominant detected peaks at yaw=0

| Signal | f [Hz] | f/1P | nearest P | modal match | rel. amp. | 3-blade reading |
|---|---:|---:|---:|---|---:|---|
| Tip edgewise X | 0.1251 | 1.00 | 1P | - | 1.000 | single-blade/local; cancels in ideal 3-blade sum |
| Tip edgewise X | 0.2518 | 2.01 | 2P | - | 0.848 | single-blade/local; cancels in ideal 3-blade sum |
| Tip edgewise X | 0.3784 | 3.02 | 3P | - | 0.048 | survives ideal 3-blade sum |
| Tip flapwise Y | 0.5035 | 4.02 | 4P | - | 1.000 | single-blade/local; cancels in ideal 3-blade sum |
| Tip flapwise Y | 0.1251 | 1.00 | 1P | - | 0.986 | single-blade/local; cancels in ideal 3-blade sum |
| Tip flapwise Y | 0.2518 | 2.01 | 2P | - | 0.568 | single-blade/local; cancels in ideal 3-blade sum |
| Torque Q | 0.1251 | 1.00 | 1P | - | 1.000 | single-blade/local; cancels in ideal 3-blade sum |
| Torque Q | 0.2518 | 2.01 | 2P | - | 0.580 | single-blade/local; cancels in ideal 3-blade sum |
| Torque Q | 0.3769 | 3.01 | 3P | - | 0.177 | survives ideal 3-blade sum |

## Marker amplitudes

| Marker | Torque rel. amp. | Flapwise rel. amp. | Edgewise rel. amp. | Reading |
|---|---:|---:|---:|---|
| 1P | 1.000 | 0.986 | 1.000 | single-blade/local; cancels in ideal 3-blade sum |
| 2P | 0.577 | 0.565 | 0.844 | single-blade/local; cancels in ideal 3-blade sum |
| 3P | 0.175 | 0.410 | 0.047 | survives ideal 3-blade sum |
| 4P | 0.160 | 0.981 | 0.025 | single-blade/local; cancels in ideal 3-blade sum |
| 1F | 0.007 | 0.178 | 0.001 | structural modal marker |
| 1E | 0.033 | 0.085 | 0.047 | structural modal marker |

## Yaw-sweep marker hierarchy

| Yaw [deg] | Torque dominant marker | Torque 1P/2P/3P | Flapwise dominant marker | Flapwise 1P/3P/4P | Edgewise dominant marker | Edgewise 1P/2P |
|---:|---|---:|---|---:|---|---:|
| 0 | 1P (1.00) | 1.00/0.58/0.18 | 1P (0.99) | 0.99/0.41/0.98 | 1P (1.00) | 1.00/0.84 |
| 10 | 1P (1.00) | 1.00/0.65/0.19 | 4P (0.98) | 0.24/0.41/0.98 | 1P (1.00) | 1.00/0.84 |
| 20 | 1P (1.00) | 1.00/0.80/0.22 | 1P (1.00) | 1.00/0.24/0.55 | 1P (1.00) | 1.00/0.85 |
| 30 | 2P (1.00) | 0.98/1.00/0.20 | 1P (1.00) | 1.00/0.11/0.25 | 1P (1.00) | 1.00/0.86 |
| 40 | 1P (1.00) | 1.00/0.39/0.04 | 1P (1.00) | 1.00/0.05/0.13 | 1P (1.00) | 1.00/0.91 |

## Reading

- At yaw=0, torque and edgewise deformation share the strongest low-order rotational content near 1P, with visible 2P and weaker 3P content.
- Across yaw, the behaviour is not identical: torque remains dominated by low-order rotational content, but flapwise dominance moves between 1P and 4P, with significant 3P content at low yaw.
- The yaw dependence is physical information, not a failure of the method: yaw reorganizes inflow azimuthally, changes phase, and changes how local harmonics survive in the rotor-equivalent signal.
- The strongest yaw=0 flapwise peak falls near 4P, not on the static 1F marker; it should be read as rotational-order content in the modal neighbourhood, not as a closed modal-identification claim.
- Peaks near integer `mP` but not at multiples of 3 should not be interpreted as true rotor-global harmonics without a phase-summed three-blade structural response.
- The clearly modal-matched peaks near 1E are weak in all three signals; they are much smaller than the low-order rotational content.
