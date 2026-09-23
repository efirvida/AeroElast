# Sub-band spectral variance of rotor torque

Source: `bem_report.csv` column `Torque [N.m]` (units: N·m, converted to kN·m).
Window: t ∈ [40.0, 100.0] s. Hann window, Parseval-coherent PSD.
The standalone corotational report uses the corotational rows; optional inertial
rows are retained only for cross-solver diagnostics when available.

## Bands

| Label | Range [Hz] | Physical content |
|---|---|---|
| `B_low`  | 0.05–1.0 | 1P (0.125), 3P (0.376), f1_flap (0.554), f1_edge (0.629) |
| `B_mid`  | 1.0–5.0  | 2nd flap (1.69), 2nd edge (1.98), 1st torsional (4.54) |
| `B_high` | 5.0–25.0 | sin contraparte modal estructural conocida |

## σ_Q por banda (corotacional y diagnostico opcional)

| Yaw [°] | Banda | σ corot [kN·m] | σ inert [kN·m] | ratio i/c |
|---:|---|---:|---:|---:|
| 0 | B_low | 238.1551 | 259.0073 | 1.09 |
| 0 | B_mid | 5.1677 | 13.1203 | 2.54 |
| 0 | B_high | 0.7410 | 16.6582 | 22.48 |
| 10 | B_low | 202.6007 | 221.9926 | 1.10 |
| 10 | B_mid | 4.9914 | 12.0539 | 2.41 |
| 10 | B_high | 0.5843 | 15.1601 | 25.95 |
| 20 | B_low | 143.1476 | 158.2142 | 1.11 |
| 20 | B_mid | 4.3601 | 9.0776 | 2.08 |
| 20 | B_high | 0.4567 | 10.6691 | 23.36 |
| 30 | B_low | 83.6104 | 91.4791 | 1.09 |
| 30 | B_mid | 3.2209 | 4.1253 | 1.28 |
| 30 | B_high | 0.2604 | 2.2127 | 8.50 |
| 40 | B_low | 65.1988 | 65.4623 | 1.00 |
| 40 | B_mid | 2.1064 | 6.3335 | 3.01 |
| 40 | B_high | 0.2010 | 7.9733 | 39.67 |

## Fracción de varianza total por banda (corotacional)

| Yaw [°] | B_low | B_mid | B_high |
|---:|---:|---:|---:|
| 0 | 0.9995 | 0.0005 | 0.0000 |
| 10 | 0.9994 | 0.0006 | 0.0000 |
| 20 | 0.9991 | 0.0009 | 0.0000 |
| 30 | 0.9985 | 0.0015 | 0.0000 |
| 40 | 0.9989 | 0.0010 | 0.0000 |

## Fracción de varianza total por banda (inercial)

| Yaw [°] | B_low | B_mid | B_high |
|---:|---:|---:|---:|
| 0 | 0.9933 | 0.0025 | 0.0041 |
| 10 | 0.9924 | 0.0029 | 0.0046 |
| 20 | 0.9922 | 0.0033 | 0.0045 |
| 30 | 0.9974 | 0.0020 | 0.0006 |
| 40 | 0.9764 | 0.0091 | 0.0145 |

## Parseval sanity check

Ratio var_freq / var_time ∈ [0.988, 1.057] — expected ~1.0.
Pequeñas desviaciones (<10 %) son aceptables y provienen de la ventana Hann + Δf finito.
