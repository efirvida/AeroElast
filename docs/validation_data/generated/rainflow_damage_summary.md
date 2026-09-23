# A4 — Rainflow / DEL analysis on $M_b$ root (flapwise root bending moment)

Source: `bem_report.csv` column `Mb [N.m]` (units: N·m, converted to MN·m here).
Window: t ∈ [40.0, 100.0] s (target span ≈ 7.5 revolutions at Ω = 7.518 RPM).
Wöhler slope m = 10 (typical glass-fibre composite, IEC 61400-1 / DNV-GL).
DEL normalised by N_eq = T_sim, so units are MN·m at 1 Hz equivalent cycle rate.

## Comparison table

| Solver | Yaw [°] | n samples | span [s] | Mb mean [MN·m] | Mb std [MN·m] | Mb p-p [MN·m] | n cycles | DEL_m10 [MN·m] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| corotational | 0 | 6001 | 60.00 | 50.8581 | 1.5144 | 5.5876 | 214 | 4.3274 |
| corotational | 10 | 6001 | 60.00 | 50.2701 | 1.3335 | 4.9734 | 191 | 3.8070 |
| corotational | 20 | 6001 | 60.00 | 48.4679 | 1.0775 | 4.0860 | 168 | 3.0751 |
| corotational | 30 | 6001 | 60.00 | 45.2603 | 0.7868 | 3.0416 | 135 | 2.3376 |
| corotational | 40 | 6001 | 60.00 | 40.3234 | 0.5844 | 2.0867 | 86 | 1.6410 |
| inertial | 0 | 3001 | 30.00 | 50.6089 | 1.5842 | 6.4813 | 459 | 5.0358 |
| inertial | 10 | 3001 | 30.00 | 50.0315 | 1.4112 | 5.8494 | 456 | 4.5260 |
| inertial | 20 | 3001 | 30.00 | 48.2449 | 1.1603 | 4.8485 | 460 | 3.7600 |
| inertial | 30 | 3001 | 30.00 | 45.0633 | 0.8730 | 3.6835 | 461 | 2.8327 |
| inertial | 40 | 3185 | 31.84 | 40.1815 | 0.6601 | 2.6211 | 483 | 1.9969 |

## Lectura

- **Tendencia con yaw**: en la rama corrotacional, el DEL del Mb raíz decrece
  con yaw junto con la carga media y la amplitud cíclica absoluta. La métrica
  captura una tendencia comparativa regular, sin picos aislados que dominen el daño.
- **Alcance de las filas adicionales**: si existen filas de otro solver, deben leerse
  solo como diagnóstico auxiliar y con atención a `span [s]`; el reporte corrotacional
  usa las filas `corotational` de la campaña de 100 s.
- **Limitación de tiempo simulado**: 60 s ≈ 7.5 revoluciones es corto para una
  reclamación de fatiga absoluta. La DEL aquí es comparativa, no absoluta. Para
  certificación IEC 61400-1 se requiere DLC 1.1 (NTM, 600 s × 6 semillas) que
  excede el alcance de este dossier.
- **Publicabilidad**: la combinación Rainflow + DEL sobre el barrido yaw de la
  IEA 15 MW con BEM+shell es una caracterización dinámica original, pero todavía
  no reemplaza una campaña DLC turbulenta.
