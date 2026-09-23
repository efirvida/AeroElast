# A7 — Cross-baseline BEM-only at V-01 and V-04 operating points

Same `BEMSolver` + same blade YAML (`IEA-15-240-RWT.yaml`) evaluated at both
the V-01 (NREL rated) and V-04 (Frontiers) operating points. This isolates
the operating-point contribution from the flexible-vs-rigid deficit.

| Operating point | V_inf [m/s] | Ω [RPM] | Torque [MN·m] | Power [MW] | Thrust [MN] |
|---|---:|---:|---:|---:|---:|
| V-01 (NREL rated) | 10.659 | 7.518 | 21.371 | 16.825 | 2.551 |
| V-04 (Frontiers) | 10.590 | 7.550 | 20.884 | 16.511 | 2.541 |

**Delta torque V-04 vs V-01 (rigid BEM only)**: -2.28%
(-0.487 MN·m).

## Decomposition of the §3.7.1 "flexible vs rigid" claim

- Dossier rigid baseline (used in §3.7.1, V-04 point): 20.081 MN·m
- This run rigid at V-04 point: 20.884 MN·m
- This run rigid at V-01 point: 21.371 MN·m
- Dossier FSI flexible (corotational at V-04 point): 18.599 MN·m

Deficit flexible vs rigid (V-04 point, §3.7.1 number): -7.38%
Deficit flexible vs rigid (V-01 point, hypothetical): -12.97%

Operating-point contribution relative to the V-01 rigid baseline: -2.28%.
This is why the aeroelastic deficit must be computed against the V-04 rigid
baseline, not against the NREL-rated V-01 rigid point.
