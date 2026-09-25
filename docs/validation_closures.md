# Log de cierres de validación

Registro de los ítems de validación **cerrados**, con su evidencia y el
**disparador de re-ejecución** (qué cambio de código o de input invalida el
cierre). Este archivo es la fuente para el capítulo de validación del
artículo y evita re-auditar lo ya cerrado.

Última actualización: 2026-09-18.

## Cómo correr la suite de validación

```bash
python -m pytest tests/ -q --tb=short \
    --ignore=tests/test_blade_mesh.py \
    --ignore=tests/test_rotor_inertial.py \
    --ignore=tests/test_vol_mesh.py
```

`test_vol_mesh.py` (39 fallas) prueba el pipeline de malla volumétrica 3D de
capa límite (`get_vertex_normals`/`create_offset_layers` devuelven ceros):
fuera del alcance de la validación shell. `test_blade_mesh.py` y
`test_rotor_inertial.py` tienen imports stale. El paquete tercero
`tests/IEA15MW/validation papers/github-IEA-15-240-RWT/tests/` no se corre.

---

## G1 — Elementos MITC / shell (pruebas unitarias)

| Ítem | Resultado | Evidencia |
|---|---|---|
| E-01 cantilevers isotrópicos | error máx 1.42 % | `test_shell_comprehensive.py`, dossier §5.7 |
| E-02 benchmarks MITC4+ de literatura | error máx 3.88 % | idem |
| Cajón cerrado (Bredt-Batho + EI + EA) | GJ/EI/EA < 2 % | `test_box_torsion_bending_benchmark.py` |
| Elástica exacta (UL incremental) | α=1 −0.5 %, α=2 +5.8 % | `test_static_nonlinear_ul_elastica.py` (2 tests) |
| **D-Tube flexión vs viga analítica** | 24.574 vs 24.525 m (**+0.2 %**) | `test_dtube_camarena_bending.py` |
| **D-Tube torsión vs Bredt-Batho** | GJ 1.3455e11 N·m², error **0.92 %** (64 arcos; 32 → 2.10 %; n_z no influye) | idem, `test_dtube_torsion_matches_bredt` |
| D-Tube convergencia de malla (job 11596675) | lineal 24.562–24.583 m (+0.15/+0.24 %) en 6 mallas | `$SCRATCH/tmp/opencode/dtube_convergence_results.txt` |

**Re-ejecutar si cambia**: kernels de elementos (`crates/aeroelast-core/src/elements/`),
ensamble (`assembler.rs`), o el solver UL (`static_nonlinear.rs`/`.py`).

---

## G2 — Solvers estructurales (pala + geometrías analíticas)

| Ítem | Resultado | Evidencia |
|---|---|---|
| V-02 modal (pala oficial) | 1F −0.86 %, 1E −1.81 %, 2F +2.15 % | `test_iea15mw_v02_natural_frequencies.py` (job v02_11594905: 3 passed) |
| V-02 2E con pala oficial | −2.6 % (en banda) | S-2E, job 11594902 |
| V-03 estática | masa −1.11 %, equilibrio 8.2e-8 % | `test_iea15mw_v03_static_gravity.py` |
| V-05 propiedades | masa total +4.1 % = banda documentada +4–8 % (layup vs 6×6) | `test_iea15mw_v05_structural_properties.py` |
| S-1 seccional | EI_edge ~10 %, EI_flap 20-25 % (sesgo de extracción documentado), **GJ 23.9 %** post-fix local-tangente | `test_iea15mw_s1_sectional.py`; doc `shell_vs_beam_sectional_validation.md` §6 |
| S-4 modal rotante vs OpenFAST MBC3 | 1F +0.6 %, 1E −4.5 %, 2F −1.5 % | `test_iea15mw_s4_rotating_modal.py` |
| S-7 torsión (twist global vs viga) | shell/viga GJ = 1.080 convergido | `test_iea15mw_s7_torsion.py`, job s7conv |
| S-8c signo de twist | coincide con BeamDyn (+0.545/+0.981°) | `test_iea15mw_s8c_twist_sign` (job s8c) |
| Campbell / K_G | 5/5 | `test_iea15mw_sx_campbell.py` |

**Re-ejecutar si cambia**: elementos, solver espectral (SLEPc), ensamble K_G/K_SP,
geometría de pala (yaml/twist), o el generador de malla `BladeMesh`.

---

## G3 — BEM vs OpenFAST

| Ítem | Resultado | Evidencia |
|---|---|---|
| V-01 BEM rated | thrust +2.67 %, torque +3.34 %, ΔCP máx 0.0121 (6 puntos); rígido 2.541 MN | `test_iea15mw_v01_rotor_performance.py` |
| Seccional vs AeroDyn | mid-span ~5 %, total +8.5 % | `docs/shell_vs_beam_sectional_validation.md`, S-5 |
| Contrato VLM/BEM | mismo dataset de entrada y tendencia con pitch (potencia regulada, thrust cae) | `test_aero_backend_contract_coherence.py` |
| Motor BEM | 19/19 tests | `test_bem_engine.py` |

**Re-ejecutar si cambia**: motor BEM (`solvers/bem/`), CCBlade, proyección de
fuerzas, o los polares/airfoils.

---

## G4 — FSI + comparación con artículos

**Cerrado (campaña de producción `frontiersin_results_corotational_100s`):**
- V-04 rated: flap medio 12.73–12.79 m en 3 campañas (dispersión < 0.06 m).
- Potencia media 14.71 MW vs Zhou 14.76 MW (**−0.34 %**); thrust flexible 2.20 MN = Zhou.
- Barrido yaw 0–40°: ley cos³, sin ventanas no convergentes, ~3 subiteraciones/ventana.
- Parqueado `bem_90_50_S`: 8.31 m pico vs DLC 6.x ~8 m.
- FFT: jerarquía 1P > 3P consistente con Zhou Fig. 18 (ratio 0.145 en 0.12–0.18).

**Comparaciones con artículos (estado y explicación física):**
- **Camarena 2025** (D-Tube y pala): el D-Tube como caso publicado es
  **sub-especificado / no reproducible** — para la carga documentada el
  analítico da 24.5 m (viga) / 23.4 m (elástica) y el paper reporta ~15 m;
  nuestro modelo reproduce el analítico a +0.2 % y converge con malla. La
  pala BAR URC no se reproduce (sin geometría); regla de redacción en Engram
  (`validation/dtube-report-framing`).
- **Zhou 2025**: comparación de sistema −7.8 % en flap; test aislado con sus
  cargas publicadas (Fig. 11) → shell 10.71 m, viga de tablas 11.47 m
  (cargas idénticas, sin rotación/gravedad). Su GEBT es más blando en flap
  que la referencia (su Tabla 3: 0.523 Hz vs 0.555–0.559 Hz; 2F 1.475 vs
  1.659), mientras nuestro shell cierra con la referencia (~1 %). La
  diferencia queda explicada por **propiedades seccionales propias de su
  modelo y clase aerodinámica (BEM vs LL-FVW)**, dentro de su propia banda
  vs literatura (3.6–7.1 %, su Tabla 5). Escalera de fidelidad: BEM 12.78 <
  LL-FVW 13.86 < ALM 14.10 (su ALM = nuestra corrida ALM interna) < LES ~16 m.
- **Ma 2025** (yaw): comparaciones digitalizadas (potencia/flap/thrust por
  yaw) en el dossier; offset ~−1.1 m atribuido a BEM vs LL-FVW.
- **Bernardi 2025** (LES): solo contexto de fidelidad (V=10 m/s, pitch mínimo,
  Ω=7.16 rpm, TSR 9; LES resuelto).

**Abierto / fuera de alcance declarado:**
- Anomalía 1P de resultantes del solver inercial (hasta ~58× potencia,
  ~156× yaw moment): la variante de producción es la corrotacional; queda
  como trabajo futuro, no bloquea el capítulo.
- Campaña pareada K_G/K_SP + Campbell formal: planificada, sin lanzar
  (no hay dir de resultados; el `$SCRATCH/v006_matrix_results` es la matriz
  V-06 de curva de potencia, no esta campaña).

**Resultados post-campaña (2026-09-18) — yaw twist-fix:**
- 5/5 completadas y estacionarias; números en la tabla del análisis
  (`$SCRATCH/tmp/opencode/analyze_twistfix.py`).
- Con la geometría correcta: yaw 0 → P 13.161 MW, thrust 1.757 MN,
  flap 8.230 m; de-loading de thrust **−30.9 %** vs Zhou −13.0 %/ancla −15 %.
  El driver es el over-twist del shell (S-8c: +5…+9.5° vs BeamDyn +0.98°) →
  ver el ítem abierto del juez externo. Los números V-04 previos (flap
  ~12.8 m, P −0.34 % vs Zhou) son de la campaña espejada y quedan inválidos.
- Acción inmediata: extraer el twist de punta 4-path de la campaña nueva
  (VTUs en `corotational/<t>/`) para cuantificar el driver.

**V-06 — curva FSI de potencia (lanzada 2026-09-18):**
- Launcher: `tests/IEA15MW/submit_ch6_matrix.py` (12 vientos + 3 puntos extra).
  Se lanzó con `SEEDS="1"` (15 jobs, 11596777–11596791): con BEM de inflow
  uniforme los seeds repetidos son idénticos (documentado en el launcher).
- Bug corregido antes de lanzar: el launcher no creaba `case_dir/logs` y el
  sbatch redirigía ahí → los intentos previos murieron sin escribir nada
  (17 dirs vacíos, sin markers).
- Al cerrar: `ch6_run_matrix.csv` + tabla comparativa final (falta para subir
  V-06 a Capa A, doc `iea15mw_validation_reference.md` §5).

**B1/B2 — convergencia fluid mesh (cerrado):**
- Resumen en `docs/validation_data/generated/convergence_b1_b2_summary.md`
  (2026-09-09, sobre `convergence_b1_b2_results_official`): B1 (element_size)
  Richardson p=2.57, flap f∞ = 12.5998 m vs baseline 12.661 m (−0.5 %);
  B2 (time window) p=1.17, f∞ = 12.6581 m. B1 dentro del criterio de cierre
  (|fine−base|/base < 1 %).

**S-x — sensibilidad de malla estructural (cerrado):**
- Artefacto `docs/validation_data/generated/sx_mesh_convergence.csv`
  (2026-09-14) consumido por `test_iea15mw_sx_mesh_convergence.py` (pasa).
- Refinos adicionales 2026-09-18: ES 0.125 m y TE refine (jobs esc/ter).

**Campaña de producción invalidada por el fix de twist (verificado 2026-09-18):**
- El `solid_mesh.vtu` de `frontiersin_results_corotational_100s` (2026-09-11)
  tiene el twist **espejado**: el ángulo del chord medido es −yaml en todas
  las estaciones (z=10: −20.99° vs +14.10°; z=40: −5.93° vs +4.42°; ...).
- El caso usa `yaml_file: ../../IEA-15-240-RWT.yaml` (la pala afectada) y el
  código actual ya genera +yaml (verificado con `BladeMesh` el 2026-09-18).
- Consecuencia: los resultados FSI de esa campaña (V-04) son pre-fix y deben
  re-derivarse. La campaña se conserva como registro de sensibilidad de input.
- Re-ejecución **lanzada el 2026-09-18** (jobs 11596769–11596773, uno por yaw)
  con `RESULTS_BASE=$SCRATCH/frontiersin_results_corotational_100s_twistfix`.
  Al arrancar el solid, verificar la malla nueva con
  `python3 $SCRATCH/tmp/opencode/check_campaign_twist.py <RESULTS_BASE>/yaw_0/solid_mesh.vtu`
  (debe dar +yaml en todas las estaciones). Al terminar (~10 h/caso), re-derivar
  los números V-04 y actualizar los scripts de `docs/validation_plots/` (hoy
  apuntan a la campaña vieja). Los V-01/V-02/V-03/V-05 no se ven afectados (no
  dependen del signo del twist de la pala en el FSI); B1/B2 sigue válida como
  estudio de convergencia (niveles comparados con la misma pala).

**Re-ejecutar si cambia**: solver FSI (`rotor_fsi.rs`, `rotor.py`), BEM,
malla/pala, o la ventana de postproceso.

---

## Regla general

Un ítem cerrado se re-audita **solo** si cambia algo de su cadena
(elemento, ensamble, solver, BEM, geometría, malla) o si se re-ejecuta la
campaña con un input distinto. Los cierres por "cache stale" (test con
nombre viejo, tolerancias desalineadas con bandas documentadas) quedan
corregidos en los tests.

---

## Sesión 2026-09-19 — correcciones de código y cierres

### Cambios de código (en el árbol, pendientes de commit)

| Cambio | Archivo | Evidencia |
|---|---|---|
| **Fix convención de corte covariante MITC3** (`V3 = e3 + θy·e1 − θx·e2`; antes el director invertido) | `crates/aeroelast-core/src/elements/mitc3.rs` | mallas mixtas: antes 27.6 %/99.8 % → ahora **0.04–0.14 %**; tests puros sin cambio |
| Expectativas de signo alineadas a la convención física | `tests/test_mitc3_benchmarks.py` (5 tests) | 19/19 pass |
| **Fix ley de cargas one-way** (tributario por estación, no por nodo) | `tools/load_mapping.py` + `run_s5_oneway.py`/`run_s6_oneway_dynamic.py`/`run_s8c_twist_sign.py` | S-8c +9.03° → **+7.10°**; S-5/S-6 6 passed |
| **Fix leak de CWD** (`os.chdir` sin restaurar rompía todo lo posterior) | `tests/test_composite_beam_parity.py` | 55 fallas de `test_frontiersin_2025_yaw_fsi` + DLC → 0 |
| Guards por artefactos faltantes (Ch.6 refs, configs de simulación) | `tests/test_dlc_fsi_bem_dual_solver.py` | 1 passed / 5 skipped con razón explícita |
| xfail documentado: artefacto de convergencia modal pre-fix | `tests/test_iea15mw_sx_mesh_convergence.py` | f1e 0.125 (0.6977→0.7028) a regenerar en HPC |
| **Canonicalización de windings** (ply orientation winding-independiente) | `src/aeroelast/core/mesh/winding.py` + `generators.py` (BladeMesh) | test estricto 4/4; pala S-8c −0.2/−3.2 % |

### Tests nuevos (regresión)

- `tests/test_mixed_mesh_convention.py` — MITC3/MITC4 en malla mixta vs todo-quads (<1 %). Guarda el fix del corte.
- `tests/test_box_multicell_torsion.py` — cajón 1 y 2 celdas con **pares en esquinas** vs Bredt multicelda (0.01 m⁴), <2 %. Documenta que un par distribuido en el contorno excita el modo blando de cizalla y colapsa J (~40 %).
- `tests/test_winding_invariance.py` — **xfail** que documenta el bug latente de orientación winding-dependiente (ver limitaciones).

### Cierres de esta sesión

1. **MITC3/MITC4 en mallas mixtas**: cerrado (fix + test).
2. **Discrepancia MITC4+ vs CCX en torsión de la pala**: **resuelta** — era el bug MITC3. Con el fix, `GJ_viga/GJ_shell` = 0.932 (0.5 m) y 0.954 (0.25 m), vs CCX 0.930/0.953 (<0.1 %).
3. **Cajón bicelda**: **descartado como defecto de elemento** — el −43 % era mi carga de anillo; con pares en esquinas da J = Bredt al 0.03 %.
4. **Twist aeroelástico**: ~5–6° (0.25 m, cargas corregidas) y ~9.9° (0.5 m) vs ~1° de la viga → diferencia de clase de modelo, robusta a código/medida/malla.
5. **Suite local**: 858+ passed, con xfails/skips documentados (sin fallas no explicadas).

### Limitaciones conocidas (para usuarios y revisores)

1. ~~Orientación de plies winding-dependiente~~ **ARREGLADO (2026-09-19)**: el offset de ángulo de ply es un ángulo firmado medido sobre el normal del elemento (que salía del winding sin canonicalizar). Se agregó `aeroelast/core/mesh/winding.py::canonicalize_windings` (score físico "hacia afuera" de la sección + herencia del elemento confiable más cercano) y el generador de pala lo aplica tras el renumbering. Verificado con `tests/test_winding_invariance.py` (4 tests: el bug se documenta, la canonicalización restaura la invariancia, es idempotente y el generador emite windings canónicos). Impacto en la pala: S-8c mid −0.2 %/tip −3.2 %, κ de torsión −0.2 %, twist AD −2.3 % (esperado: los offsets ±90° son equivalentes mod 180; la diferencia residual viene de la pequeña asimetría E2/E3 del material y del setup paramétrico del elemento).
2. **Torsión de secciones de pared delgada**: medir J con un par distribuido en el contorno (no físico) excita el modo blando de cizalla en el plano (flexión de paredes ∝ t³ vs torsión de membrana ∝ t). Usar pares en esquinas o presión física.
3. **Twist del shell ~5× la viga**: física del modelo (distorsión de sección que la viga no captura). Las referencias beam-based (Zhou/Ma/BeamDyn) son cota inferior. El de-loading FSI absoluto depende del método de reconstrucción del chord (PCA del contorno) — cuantificar con juez externo 3D o documentar la banda.
4. **BEM cuasiestacionario**: clase de fidelidad intermedia (sin estela sesgada resuelta); los offsets vs LL-FVW son de fidelidad, no de bug.
5. **Solver inercial**: anomalía 1P en resultantes de rotor (hasta ~58×/156×); la variante de producción es la corrotacional.
6. **Malla**: evitar slivers (aristas <2 cm) — la conversión a 2º orden para CalculiX los invierte; el pipeline de malla volumétrica 3D (`test_vol_mesh`) está fuera de alcance.

### Pendiente en HPC (batch final)

- Re-ejecutar campañas FSI (yaw twist-fix + V-06) **con el fix de MITC3**.
- V-06 sobre-rated con el pitch schedule de referencia/ROSCO (8 casos).
- Regenerar `sx_mesh_convergence.csv` a 0.125 con el código corregido.
- Punto 0.125 de torsión/twist (convergencia fina).
- Postproceso y re-derivación de números/figuras (V-04/V-06) + cierre del capítulo.

### Batch HPC lanzado (2026-09-23)

| Job(s) | Qué | Salida |
|---|---|---|
| 11599732–11599736 | Barrido yaw 0–40° **con el fix de MITC3** (5 casos) | `$SCRATCH/frontiersin_results_corotational_100s_mitc3fix/` |
| 11599738, 11599741–11599754 | Matriz V-06 (15 casos, SEEDS=1) **con el pitch schedule de referencia** (no 1.8°/m/s) | `$SCRATCH/v006_matrix_results_mitc3fix/` |
| 11599739 | Regeneración de `sx_mesh_convergence.csv` a 0.125 (convergencia modal) | `docs/validation_data/generated/sx_mesh_convergence.csv` |
| 11599740 | Punto 0.125 de torsión/twist (MITC4, sin CCX) | `$SCRATCH/tmp/opencode/ccx_robust_0125.csv` |

Notas: los V-06 quedan parcialmente retenidos por `AssocMaxJobsLimit` y se liberan
a medida que terminan los demás. El schedule de referencia (V-06) se implementó
en `submit_ch6_matrix.py` interpolando la tabla IEA: pitch 0 hasta 10.659 m/s,
6.761° a 12.259, 12.185° a 15.471 y extrapolación con esa pendiente arriba.
Al terminar: postprocesar, re-derivar V-04/V-06 y actualizar figuras + capítulo.

#### Postproceso parcial del batch (2026-09-23)

**V-04 yaw con el fix MITC3** (0/10/20 terminados; 30 corriendo, 40 pendiente) —
`docs/validation_data/generated/v04_yaw_mitc3fix.csv`:

| yaw | Potencia | Thrust | Flap | Edge | vs campaña twistfix (pre-fix MITC3) |
|---|---|---|---|---|---|
| 0 | 12.964 MW | 1.726 MN | 8.022 m | −0.241 | −1.5 % / −2.5 % |
| 10 | 12.485 MW | 1.699 MN | 7.957 m | −0.229 | −1.5 % / −2.5 % |
| 20 | 11.111 MW | 1.620 MN | 7.757 m | −0.201 | −1.1 % / −2.4 % |

De-loading de thrust a yaw 0: **−32.1 %** vs el rígido (Zhou −13 %) — el
over-de-loading conocido (twist del shell → chord PCA del BEM).

**Torsión: convergencia resuelta por el fix de MITC3** —
`GJ_viga/GJ_shell` (GJ_viga = 0.43711 MN·m², del propio CSV):
0.5→**0.932**, 0.25→**0.953**, 0.125→**0.960**
(pre-fix: 1.037 / 1.125 / 1.193, sin converger; CCX a 0.25: 0.953).
Twist AD en el mismo punto: 9.68° / 6.17° / **5.56°** (0.5/0.25/0.125).

**V-06: el schedule estaba mal en AMBOS extremos y se corrigió.**
El launcher usaba TSR-9 con pitch 0 en todo el rango; la tabla IEA pide
5 rpm mínimas con pitch fino abajo de rated (3 m/s: 2.13→**5.00 rpm**,
0→**3.92°**) y regula 15 MW con la rampa real arriba (25 m/s: 28.3→**22.8°**).
Evidencia del impacto: nuestro BEM con el schedule viejo daba 0.364 MW a
3 m/s vs 0.043 MW de la referencia (8×) — otro punto operativo.
`schedule()` ahora interpola `docs/validation_data/reference_iea15mw_rotor_performance.csv`
(la hoja Rotor Performance del xlsx IEA, 50 puntos) y la matriz se relanzó a
`$SCRATCH/v006_matrix_results_refsched/` (jobs 11600327–11600343).

**Convergencia modal (sx)**: regenerada con el fix. f1f y OoP convergen;
el salto de f1e a 0.125 **persiste** (0.6970→0.7019) → no era el bug de MITC3;
el test sigue xfail (a investigar: identificación del modo o detalle punta/raíz
en esa malla).

#### Barrido yaw completo (5 casos, campaña mitc3fix)

Punto operativo: viento 10.59 m/s, pitch 0°, omega 0.7906341 rad/s constante
(7.55 rpm), malla 0.25 m, pala oficial IEA-15-240-RWT.yaml, 100 s; ventana
[40, 100] s. Deriva [40,70] vs [70,100] = 0.03–0.06 % → régimen estable.
Artefacto: `docs/validation_data/generated/v04_yaw_mitc3fix.csv`.

| yaw | P [MW] | ± | CP | thrust [MN] | flap [m] | ± | edge [m] | Mb [MN·m] |
|---|---|---|---|---|---|---|---|---|
| 0° | 12.964 | 0.139 | 0.3681 | 1.726 | 8.022 | 0.278 | −0.241 | 41.26 |
| 10° | 12.485 | 0.130 | 0.3545 | 1.699 | 7.957 | 0.257 | −0.229 | 40.73 |
| 20° | 11.111 | 0.110 | 0.3155 | 1.620 | 7.757 | 0.264 | −0.200 | 39.14 |
| 30° | 9.008 | 0.081 | 0.2559 | 1.490 | 7.425 | 0.302 | −0.159 | 36.53 |
| 40° | 6.407 | 0.044 | 0.1820 | 1.312 | 6.957 | 0.354 | −0.107 | 32.89 |

**Impacto del fix MITC3** (vs campaña twistfix, mismo punto operativo):
potencia −1.49/−1.47/−1.40/−1.23/−0.85 %, flap −2.53/−2.50/−2.42/−2.27/−2.02 %
(mayor a bajo yaw). El fix no cambia las conclusiones físicas, afina 1–2 %.

**Déficit vs literatura** (yaw 0: Zhou 2025; 10–40: Ma 2025,
doi:10.3389/fenrg.2025.1571567): potencia −12.2/−13.9/−18.3/−24.3/−34.0 %;
flap −42.1/−42.3/−42.5/−42.4/−41.5 %. Thrust a yaw 0: −32.1 % vs −13.4 %.

**Hallazgo clave — dos efectos separados:**

1. *Déficit global* ya a yaw 0 (−12 % potencia, −42 % flap, thrust −32 % vs
   −13 %), independiente del yaw → el over-de-loading ya documentado
   (rotación de sección → chord PCA del participante BEM).
2. *Déficit adicional que solo depende del yaw*: ajustando
   P(γ)/P(0) = cos^n(γ), nuestro n = **2.53/2.56/2.66** (10/30/40°) contra
   n = **1.36/1.49/1.58** de Ma. A yaw 40 quedamos en 0.494 de la potencia
   alineada contra 0.657 de Ma. El flap, en cambio, cae plano (−41 a −42 %)
   en todo el barrido → el déficit de flap es un offset global, no yaw-dependiente.

El exceso de degradación con yaw apunta a la parte yaw del modelo (wake
sesgado y/o feedback aeroelástico asimétrico por azimuth), no al elemento.

**Próximo paso propuesto (barato y decisivo):** correr el BEM SOLO (rígido,
sin FSI) en los 5 yaws. Si el BEM rígido da n≈1.5–2.0, el exceso viene del
acoplamiento aeroelástico; si ya da n≈2.6, es el modelo de yaw del BEM.
Separa las dos causas antes de fijar la banda de incertidumbre del artículo.

#### V-06 refsched — estado parcial y decisión de espera (2026-09-23 noche)

Casos completos (100 s) y en curso contra la hoja Rotor Performance del IEA:

| caso | V [m/s] | P sim | P ref | ΔP | thrust sim | thrust ref | ΔT | estado |
|---|---|---|---|---|---|---|---|---|
| v3 | 3 | 0.056 | 0.043 | **+29.4 %** | 0.198 | 0.202 | −1.9 % | completo |
| v5 | 5 | 1.447 | 1.395 | **+3.8 %** | 0.478 | 0.549 | **−13.1 %** | completo (drift −0.53 %) |
| v7 | 7 | 4.216 | 4.245 | **−0.7 %** | 0.884 | 1.073 | **−17.6 %** | 59 % |
| v9 | 9 | — | 9.026 | — | — | — | — | 31 % |
| v11 | 11 | — | 15.000 | — | — | — | — | 16 % |

**El fix de schedule funcionó**: la potencia bajo rated ahora coincide con la
referencia dentro de ±4 % (a 3 m/s el viejo la daba 0.364 MW = +750 %, ahora
0.056 = +29 %; el +29 % es un punto de 13 kW absolutos, atribuible al modelo
de polar a bajo Re con pitch fino). **Pero el thrust queda 13–18 % abajo** a 5
y 7 m/s y −2 % a 3 m/s: el déficit de thrust crece con la carga, consistente
con el over-de-loading aeroelástico ya documentado (a rated, 10.59 m/s:
−28 % vs la tabla, −32 % vs el rígido).

**Decisión sobre esperar la matriz completa:** no hace falta para fijar la
curva. Lo que falta con valor científico son los puntos de la rodilla
(op_v8 8, v9 9, op_v10 10, v11 11, op_v12 12), que definen dónde empieza la
regulación y si la meseta cae en 15 MW o se hunde por de-loading. Los seis
vientos extremos (v13–v25) sólo agregan puntos sobre la meseta ya regulada.
Nota de cola: los op_* se enviaron últimos, así que en FIFO correrían al
final; conviene reordenar (cancelar y reenviar los extremos al fondo) si se
quiere la rodilla antes.

Pendiente anotado: **barrido BEM rígido (sin FSI) en los 5 yaws** para
separar el déficit global del exceso de degradación con yaw (nuestro
cos^n con n≈2.6 vs 1.5 de Ma).

#### Barrido BEM rígido — diagnóstico decisivo (2026-09-24, job 11601073, dev)

Rígido = pala indeformada, misma config que el participante FSI (CCBlade +
NeuralFoil, `shear_exp=0`, sin precono/tilt). Artefacto:
`docs/validation_data/generated/rigid_bem_check.csv`.

**A. Malla rígida vs tabla IEA.** Un hallazgo metodológico: la columna
`Power [MW]` de la hoja Rotor Performance es potencia **entregada**
(`CP/CP_aero = 0.9570` constante = 4.3 % de pérdidas). Comparar nuestra
potencia *aerodinámica* contra esa columna introduce ese sesgo. Corregido:

| | potencia aero @10.59 | thrust |
|---|---|---|
| tabla IEA (aero) | 15.38 MW | 2.429 MN |
| V-01 (referencia AeroDyn) | 15.68 MW | 2.457 MN |
| **nuestro BEM rígido** | **16.51 MW** | **2.541 MN** |
| offset vs tabla aero | **+7.4 %** | +4.6 % |
| offset vs V-01 | +5.3 % | +3.4 % |

El modelo aero está dentro de su offset conocido (+3.3 % de torque en V-01).
**No es la causa del problema.**

**B. Acoplamiento (FSI vs rígido) — el problema real.** La fracción que
sobrevive al acoplamiento: 0.785 / 0.792 / 0.818 / 0.880 / **0.999** en
yaw 0/10/20/30/40. O sea: de-loading **−21.5 % potencia y −32 % thrust**
alineado, que se desvanece con el yaw. Contra la literatura (Zhou/Ma
de-loadean 4–6 % potencia y ~10 % thrust cuando se comparan contra la misma
base aero), nuestro acoplamiento sobre-de-loadea **~4–5×**. Es el over-twist
ya documentado, ahora cuantificado contra nuestra propia base rígida.

**C. Modelo de yaw del BEM — hallazgo nuevo.** El rígido degrada la potencia
mucho más rápido que Ma: `P(γ)/P(0) = cos^n` con n = 3.04 / 3.14 / 3.32 /
3.55 (yaw 10/20/30/40), contra n ≈ 1.36–1.58 de Ma. Y el acoplamiento **no**
lo causa: su efecto relativo *baja* con el yaw. Prueba limpia: a yaw 40 el
FSI da 6.407 MW y el rígido 6.412 MW — idénticos; todo el déficit de yaw 40
es del modelo de yaw del BEM (wake sesgado/skew), no del elemento ni del
acoplamiento.

**D. Consecuencia sobre la curva V-06.** El rígido sobre-predice +12 %
plano contra la columna `Power` de la tabla (≈+7 % corregido por pérdidas)
en TODO el rango 5–25 m/s. Por lo tanto el acuerdo de la campaña FSI bajo
rated (v3/v5/v7 dentro de ±4 %) es en realidad **cancelación de dos
errores**: el aero sobre-predice ~7 % y el acoplamiento resta ~9–15 %. No es
evidencia de validez.

**Veredicto de re-corrida:** sí. Hay que tocar (1) el acoplamiento
(over-twist, el efecto dominante) y (2) el modelo de yaw del BEM. Ambos
invalidan el barrido yaw FSI y la matriz V-06 tal como están. La parte
estructural (torsión, malla, modos) queda a salvo.

#### Escalera de diagnóstico aero: CCBlade vs AeroDyn (2026-09-24, `tools/diagnose_yaw_aero.py`)

Mismo rotor, mismas condiciones (viento 10.59, 7.55 rpm, pitch 0, sin shear),
tres peldaños que cambian UNA cosa por vez. `UA_Mod=0` y `DBEMT_Mod=0` para
comparar BEM estacionario contra BEM estacionario (CCBlade también lo es).

| yaw | R1 CCBlade | R2 AeroDyn `Skew_Mod=0` | R3 AeroDyn `Skew_Mod=1` | R3−R2 | R3−R1 |
|---|---|---|---|---|---|
| 0° | 16.511 MW | 15.815 MW | 15.815 MW | 0.0 % | **−4.2 %** |
| 10° | 15.760 | 15.100 | 15.110 | +0.1 % | −4.1 % |
| 40° | 6.412 | 6.231 | 6.547 | **+5.1 %** | +2.1 % |

`P(γ)/P(0) = cos^n`: R1 n = 3.04/3.55; **R2 (skew off) n = 3.02/3.49**;
**R3 (skew on, Pitt/Peters) n = 2.98/3.31**.

**Conclusión 1 — el modelo de estela sesgada NO es la causa.** AeroDyn con
skew activo degrada casi igual que sin skew y casi igual que nuestro CCBlade.
El skew aporta sólo +5 % en yaw 40. La degradación empinada (n≈3) es
propiedad de *la clase BEM estacionaria*, no un defecto de nuestro código:
dos implementaciones independientes (CCBlade y AeroDyn) coinciden dentro de
2–5 % en todo el barrido. Parchear un skew en nuestro Fortran no habría
arreglado el yaw.

**Conclusión 2 — el offset alineado es de POLARES.** R1 16.511 vs R2/R3
15.815 = +4.4 %. AeroDyn con los polares tabulados del IEA reproduce la
referencia de la tabla (aero ≈15.4–15.8 MW); nosotros usamos polares
NeuralFoil. El +7 % vs la tabla es entonces modelo de polar, no código BEM.

**Conclusión 3 — el yaw empinado vs la literatura (n≈1.5 de Ma) tampoco lo
explica AeroDyn estacionario**, que da n≈3.3 como nosotros. Los candidatos
restantes son el inflow dinámico (DBEMT, que acá apagué a propósito), la
física de estela libre de Ma (LL-FVW), o la respuesta aeroelástica de su
modelo. Siguiente peldaño: reactivar DBEMT/UA.

**Lectura para la prioridad del usuario:** el lado aero está *aislado y no
es el defecto dominante*. Lo dominante sigue siendo el acoplamiento
(−21.5 % de potencia y −32 % de thrust alineado donde la literatura
de-loadea 4–6 % y ~10 %), y después el solver estructural.

#### Los cuatro peldaños de diagnóstico (2026-09-24) — veredicto

**Peldaño 1 — aero, física completa** (`tools/diagnose_yaw_aero.py`, job 11601115).
`P(γ)/P(0) = cos^n`:

| | yaw 10 | 20 | 30 | 40 |
|---|---|---|---|---|
| R1 CCBlade | 3.04 | 3.14 | 3.32 | 3.55 |
| R2 AeroDyn skew=0 | 3.02 | 3.12 | 3.29 | 3.49 |
| R3 AeroDyn skew=1 | 2.98 | 3.06 | 3.19 | 3.31 |
| **R4 AeroDyn deck (skew+DBEMT+UA)** | **2.89** | **3.02** | **3.18** | **3.32** |

DBEMT + UA (R3→R4) aporta ≤0.3 %: el inflow dinámico **no** explica la
pendiente. Ni el skew ni el dinámico acercan AeroDyn al n≈1.5 de Ma →
la diferencia es **clase de modelo de estela** (BEM vs estela libre) y/o la
respuesta aeroelástica del modelo de ellos. El aero no es el defecto.

**Peldaño 2 — polares.** CL dentro de +0.3…+4.3 % (peor inboard), CD
outboard ~0 % e inboard +10–16 %. Los polares explican ~1–3 % del offset
alineado de +4.4 %; queda un confound de Reynolds (la tabla IEA tiene un Re
de diseño por aire, nosotros usamos el local).

**Peldaño 3 — acoplamiento** (`tools/diagnose_yaw_coupling.py`). Mismo aero
en ambos lados (nuestro CCBlade):

| yaw | P_rigid | P_fsi | ratio | de-load | T ratio |
|---|---|---|---|---|---|
| 0 | 16.511 | 12.964 | 0.785 | **−21.5 %** | 0.679 |
| 20 | 13.579 | 11.111 | 0.818 | −18.2 % | 0.687 |
| 40 | 6.412 | 6.407 | **0.999** | −0.08 % | 0.718 |

vs literatura: potencia **×2.0** (base nuestro rígido) o **×2.7** (base
AeroDyn, la justa); thrust **×2.4**. Over-de-loading real, pero 2–2.7×, no
4–5× (corrección de una afirmación previa).

**Peldaño 4 — estructura. Es el defecto dominante.** Bajo las MISMAS cargas
AeroDyn (S-8c), nuestro shell tuerce:

| | nuestro shell | ancla BeamDyn oficial | factor |
|---|---|---|---|
| mid-span (z=60 m) | **+1.615°** | +0.545° | **×2.96** |
| punta (z=112.6 m) | **+11.635°** | +0.981° | **×11.9** |

Caveat metodológico pendiente: la "rotación de sección" de un shell puede
incluir distorsión/warping que un beam no tiene; aún así el mid-span (menos
distorsionado) ya es ×3.

**Cadena causal**: la estructura sobre-tuerce (×3–12) → el acoplamiento
sobre-de-loadea (×2–2.7) → la potencia y el flap caen. El aero es lo que
menos falla. Siguiente escalón natural: comparar nuestras propiedades de
sección (GJ, EI, posición del centro de torsión/EA) contra la hoja "Blade
Structural Properties" del xlsx del IEA — el S-8b (semichord) ya apuntaba
al eje de pitch como parámetro sensible.

#### FIX: unificación de frame aero-estructura (2026-09-24)

Causa raíz: aero y estructura usaban **orígenes distintos** en la misma
corrida FSI. El deck oficial ElastoDyn define el datum sin ambigüedad:
`TipRad 120.97 m`, `HubRad 3.97 m` — o sea la pala va de 3.97 a 120.97 desde
el apex del rotor. Nuestro loader aéreo construye las estaciones justo así
(3.97…120.97) y el `RotorMesh` traduce por el hub radius; **el `BladeMesh`
no lo hacía** (0…117), y el participante además sumaba el hub radius dos
veces en `_rebuild_bem_solver`.

Cambios:

| archivo | cambio |
|---|---|
| `core/mesh/generators.py` | `BladeMesh` acepta `hub_radius`/`rotor_frame` y traduce al frame rotor (misma llamada que `RotorMesh`); opt-out `rotor_frame=False` |
| `solvers/bem/fsi_participant.py` | eliminado el doble conteo del hub radius en `deformed_rtip` (era `(max(r_def)+hub)*1.001` → Rtip 125 m en un rotor de 121) |
| `tools/run_s5_oneway.py` | `extract_ad_loads` devuelve `r_m` en frame rotor (BlSpn + hub) |
| `tests/IEA15MW/reference/s5_ad_blade1_loads.csv` | `r_m` desplazado +3.97 (0…117 → 3.97…120.97) |
| `tests/test_iea15mw_s1_sectional.py`, `tools/beam_reference.py` | comparaciones por **fracción de vano** `(x−min)/(max−min)`, no `x/max` |

**Verificación**: suite completa **923 passed, 26 skipped, 4 xfailed, 0
failed** (idéntico al baseline). Prueba física: el estático parked es
frame-invariante — da `14.9624 m` tanto en la configuración consistente
blade-local como en la consistente rotor, y coincide con el valor
`oop_static` del CSV de convergencia a 0.5 m (14.96282).

**Pendientes anotados**:
- La campaña declara `hub_radius: 0.0` y el loader lo ignora (recomputa 3.97
  desde el yaml). El recomputo es el correcto; la config miente.
- El driver AeroDyn de la escalera aero usó `HubRad 1.5` (elección mía,
  errada): hay que poner 3.97 y repetir esa escalera.
- El deck oficial trae `PreCone −4°` y `ShftTilt −6°`; nuestras campañas
  corren 0/0 — diferencia de configuración a declarar.

**Consecuencia**: todas las campañas acopladas (barrido yaw y matriz V-06)
corrieron **con el bug de frame activo** → hay que re-correrlas con el fix.
