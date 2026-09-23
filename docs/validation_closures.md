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
