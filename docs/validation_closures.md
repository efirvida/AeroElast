# Log de cierres de validación

Registro de los ítems de validación **cerrados**, con su evidencia y el
**disparador de re-ejecución** (qué cambio de código o de input invalida el
cierre). Este archivo es la fuente para el capítulo de validación del
artículo y evita re-auditar lo ya cerrado.

Última actualización: 2026-10-08.

> **El orden de trabajo lo rige el roadmap: issue #18** ("roadmap: what is still unvalidated or
> unexplained (order of record)"). Es la única fuente del orden. Cada ítem abierto de #18 trae su
> *Entry* (qué leer primero) y su *closes when*; cada ítem cerrado apunta a su sección de este
> archivo. Un ítem se ataca **por separado, en una sesión limpia sin contexto**, y la sesión que lo
> cierra deja en el mismo push: su sección acá, el cambio en el store y el comentario en el issue.
> Este log es el espejo local de los **cierres**; `docs/validation/gaps.yaml` es el espejo de lo que
> queda **sin validar** (con `citations_forbidden: true` donde no se puede citar).

> **2026-09-30 — la línea de elementos cambió.** Se integró `origin/main` (el elemento
> MITC4+/D revisado, el fix de ángulos de ply span-relative, el fix de corte no corregido,
> el fix del writer CCX y los fixes del solver modal). **Todos los cierres de G1 y G2 de
> abajo quedaron invalidados**: son el disparador de re-ejecución que este documento define.
> Estado medido antes/después, regresiones y próximos pasos:
> `docs/origin_main_integration_2026-09-30.md`.  Resumen: 933→905 passed, 6→15 failed,
> 0→22 errors; 28 tests de CalculiX bloqueados por un elemento degenerado de la malla de
> pala; V-02 2F, S-4 rotante, S-7, S-6 (NaN), box EI, D-Tube y elástica se movieron.

> **2026-10-05 — la cadena de signo del momento torsor de la pala cambió.** El commit
> `1146265` corrige tres errores de signo enlazados en la carga torsional de la pala,
> arbitrados contra la geometría de aerofoil del propio deck por
> `tools/diagnose_sign_chain.py` (cada anillo de estación real emparejado con su aerofoil
> WindIO, residuo 0.5–3.6 % de cuerda). La convención queda enunciada una sola vez aquí: el
> deck pone el borde de ataque en **+x** y el empuje del marco de carga (downwind) en
> **+y**, así que una rotación rígida **+z** mueve el borde de ataque aguas abajo y
> **nose-down es `omega > 0`**. Consecuencia medida en el punto rated: giro de sección en la
> punta **`+8.1048°`** (mínima norma) y **`+9.6669°`** (multi-celda con propiedades);
> de-loading (flexible−rígido)/rígido: **twist solo `−26.01 % / −15.61 %`**, **radios solo
> `+1.24 % / +0.95 %`**, **producción `−25.31 % / −14.77 %`** (thrust/potencia), contra Zhou
> Table 6 **`−13.04 % / −8.38 %`**. Las cifras de twist y de de-loading **pre-fix** de las
> secciones de abajo (barrido yaw `−30.9 %`/`−32.1 %`, over-twist ×3–12, S-8c `+5…+9.5°`)
> son registro histórico: no se arrastran como estado actual. Detalle en "Fix de la cadena
> de signo (2026-10-05)" al final.

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

Desde el merge de `origin/main` (2026-09-30) el soporte de sólidos 3D ya no
existe en Rust ni en Python, y upstream borró `test_solid_elements.py`,
`test_vol_mesh.py` y `test_beam_4cases_parity.py`; el `--ignore` de
`test_vol_mesh.py` quedó como no-op por compatibilidad.

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
| S-7 torsión (twist global vs viga) | twist **0.9336** convergido ⇒ shell/viga `GJ` **1.071** (re-medido 2026-10-07; el `1.080` del job `s7conv` es pre-merge y quedó invalidado por el merge del 30-09) | `tests/validation/blade/test_iea15mw_s7_torsion.py`; § S-7 (2026-10-07) |
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
flap 8.230 m; de-loading de thrust **−30.9 %** vs Zhou −13.0 %/ancla −15 %. *(pre-fix del
signo, ver fix de la cadena de signo 2026-10-05.)*
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
- ~~Regenerar `sx_mesh_convergence.csv` a 0.125 con el código corregido.~~
  **ENTREGADO (2026-09-23)**: job **11599739** COMPLETED en 2h27m; el CSV
  regenerado está en `docs/validation_data/generated/sx_mesh_convergence.csv`
  (4 tamaños: 1.0 / 0.5 / 0.25 / 0.125 m; log
  `$SCRATCH/tmp/opencode/sxc_11599739.out`, que termina con
  `wrote docs/validation_data/generated/sx_mesh_convergence.csv`). El artefacto
  queda, sin embargo, **obsoleto por otra causa** (ver la corrección del batch,
  abajo): **no** es un artefacto pre-fix de MITC3.
- ~~Punto 0.125 de torsión/twist (convergencia fina).~~ **ENTREGADO (2026-09-23)**:
  job **11599740** COMPLETED en 9m34s; salida
  `$SCRATCH/tmp/opencode/ccx_robust_0125.csv`. **Parcial**: solo MITC4-vs-viga
  (`kappa_ccx = nan`, `ad_ccx_089 = nan`), sin valor de CCX en ese refinamiento.
- Postproceso y re-derivación de números/figuras (V-04/V-06) + cierre del capítulo.

### Batch HPC lanzado (2026-09-23)

| Job(s) | Qué | Salida |
|---|---|---|
| 11599732–11599736 | Barrido yaw 0–40° **con el fix de MITC3** (5 casos) | `$SCRATCH/frontiersin_results_corotational_100s_mitc3fix/` |
| 11599738, 11599741–11599754 | Matriz V-06 (15 casos, SEEDS=1) **con el pitch schedule de referencia** (no 1.8°/m/s) | `$SCRATCH/v006_matrix_results_mitc3fix/` |
| 11599739 | Regeneración de `sx_mesh_convergence.csv` a 0.125 (convergencia modal) — **COMPLETED 2026-09-23, 2h27m** | `docs/validation_data/generated/sx_mesh_convergence.csv` (log `$SCRATCH/tmp/opencode/sxc_11599739.out`) |
| 11599740 | Punto 0.125 de torsión/twist (MITC4, sin CCX) — **COMPLETED 2026-09-23, 9m34s** | `$SCRATCH/tmp/opencode/ccx_robust_0125.csv` (parcial: `kappa_ccx = nan`) |

Notas: los V-06 quedan parcialmente retenidos por `AssocMaxJobsLimit` y se liberan
a medida que terminan los demás. El schedule de referencia (V-06) se implementó
en `submit_ch6_matrix.py` interpolando la tabla IEA: pitch 0 hasta 10.659 m/s,
6.761° a 12.259, 12.185° a 15.471 y extrapolación con esa pendiente arriba.
Al terminar: postprocesar, re-derivar V-04/V-06 y actualizar figuras + capítulo.

#### Entrega y corrección de los dos ítems 0.125 (actualizado 2026-10-05)

Los dos ítems de la tabla ligados a la malla fina **ya se entregaron**: jobs
**11599739** (convergencia modal) y **11599740** (torsión/twist), ambos
**COMPLETED** el 2026-09-23. Lo que sigue **corrige el diagnóstico** con el que
quedaron registrados; no borra la entrada previa.

- **`sx_mesh_convergence.csv` — el diagnóstico "predata el fix de MITC3" queda
  falsado.** El artefacto **sí** es el regenerado: el job 11599739 corrió la
  malla de 0.125 m (`120352 nodes, 122587 elements`, `1E parked 0.7019`) y
  escribió el CSV al final. Sus filas de f1e park son
  **0.69569824** (1.0 m) / **0.69661868** (0.5) / **0.69698433** (0.25) /
  **0.70194314** (0.125), con incrementos consecutivos 0.13 % / 0.05 % /
  **0.71 %**: el salto **persiste después del fix de MITC3** (~14× el incremento
  anterior), así que no era ese artefacto. El resumen del propio job **no**
  imprimía f1e (solo `1F parked`, `1F rotating`, `OoP static`:
  `1.86 % | 0.40 % | 0.11 %`, `1.28 % | 0.28 % | 0.08 %`,
  `8.59 % | 1.82 % | 0.46 %`), y por eso el salto pasó desapercibido.
- **Causa correcta de obsolescencia del artefacto**: se generó con la
  canonicalización de winding todavía en su lugar (el log dice
  `Canonicalised 121305 element winding(s)`, y `core/mesh/winding.py` ya se
  eliminó para alinear con `origin/main`) y **antes** de la integración del
  elemento de `origin/main`. La instrucción de regenerar sigue en pie, con esa
  causa.
- **Punto 0.125 de torsión/twist — parcial**: `ccx_robust_0125.csv` trae
  `ES 0.125`, `kappa_mitc4 0.4551740622`, `kappa_beam 0.4371100257`,
  `gj_ratio_mitc4 0.960314`, `ad_mitc4_089 5.5587`; el lado CCX **no** produjo
  valor en ese refinamiento (`kappa_ccx = nan`, `ad_ccx_089 = nan`) → el punto es
  MITC4-vs-viga únicamente, sin referencia CCX a 0.125.
- **Hallazgo abierto**: el incremento de f1e a 0.125 m es **no monótono**
  (0.71 % contra 0.05 % previo) y **sobrevive al fix de MITC3**; el artefacto
  regenerado es ahora obsoleto por la remoción de la canonicalización de winding
  y la integración del elemento de `origin/main`. El test
  `TestSxMeshConvergence::test_f1e_frequency_converges` sigue en
  `xfail(strict=True)` hasta regenerar contra la línea de elemento actual.

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

> **Actualizado (2026-10-05)**: ver "Entrega y corrección de los dos ítems 0.125"
> más arriba. El "a investigar" queda acotado: el artefacto regenerado por el job
> 11599739 es obsoleto por la remoción de la canonicalización de winding
> (`core/mesh/winding.py`) y por la integración del elemento de `origin/main` —
> no por el fix de MITC3, que el salto de f1e sobrevive.

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
menos falla. *(Los números de de-loading de esta cadena son pre-fix del signo; ver el
fix de la cadena de signo, 2026-10-05. La estructura sobre-tuerce sigue siendo el hallazgo,
pero el de-loading corregido es `−25.31 % / −14.77 %`, sobre Zhou en vez de corto.)*
Siguiente escalón natural: comparar nuestras propiedades de
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

### Fix de la cadena de signo del momento torsor (2026-10-05)

Tres errores de signo enlazados escondían el signo del twist y del de-loading FSI, y el
tercero compensaba al primero. Convención arbitrada una sola vez: borde de ataque en **+x**,
downwind del marco de carga en **+y**, rotación rígida **+z** = **nose-down**, o sea
**`omega > 0`**. La fija `tools/diagnose_sign_chain.py` contra la geometría del deck
(residuo 0.5–3.6 % de cuerda en cada estación real).

1. El momento de cabeceo aerodinámico se aplicaba sobre el eje equivocado: la dirección
   leading-to-trailing del deck corre **contra** `_strip_chord_dirs` en toda estación real
   (`LE . c_hat = -1.0`). `ForceProjector` guarda `_strip_moment_axis_sign` por strip y
   aplica `Mp` sobre él.
2. Los dos guards afirmaban `omega < 0` (nose-up); ahora afirman `omega > 0`.
3. `fsi_participant._twist_mesh_to_bem` era `-1` y convertía la rotación nose-up pre-fix en
   un descenso aparente de `-4.00 %`; ahora es **+1**.

Evidencia, mismo punto rated: `-0.8696 -> +8.1048 deg` (mínima norma),
`+0.0773 -> +9.6669 deg` (multi-celda con propiedades); de-loading twist solo
`-4.00 %/-0.81 % -> -26.01 %/-15.61 %`; radios solo `+1.13 %/+0.85 % -> +1.24 %/+0.95 %`;
producción `-2.96 %/+0.03 % -> -25.31 %/-14.77 %`; Zhou Table 6 `-13.04 %/-8.38 %`
(sin cambio). La validación Bredt del tubo queda intacta en `1.0031x` (0.3091 %). El guard que
habría atrapado todo esto es la invariancia sin convención de
`tests/test_multicell_shear_flow.py`: las dos realizaciones fieles del mismo momento rated
deben rotar la sección en el mismo sentido. Lo que queda es una cuestión de **magnitud**
(~1.9× por encima de Zhou), no de signo (acoplamiento one-way vs convergido, o nivel de
carga), más la prohibición de citar números pre-fix. Ver `docs/validation/gaps.yaml`
(`moment_realization_over_delivers`, `force_projection_sense_p5`).

### El artefacto axial de la realimentación de radios: REFUTADO (2026-10-07)

Cierra la premisa de issue #13 registrada en el gap `force_projection_axial_extension`: que
parte del `+1.24 %` de thrust de la fila de radios fuera un artefacto de medir un *path length*
sobre el eje curvado de la pala, y no una extensión axial física. La premisa se mide y se
refuta. Es la segunda y última causa registrada para ese gap: la primera, la fuerza spanwise
que inyecta la distribución nodal, ya se había medido y refutado (sección 22.7 del veredicto,
`odd/tasks/composite-bend-twist-verdict.md`), así que el gap se cierra como **refutado**, no
como pendiente.

Tres hechos de código, en `src/aeroelast/solvers/bem/fsi_participant.py`:

1. La definición deformada es una **proyección de span**
   (`r_def[k] = mean((X + u)[nodos del strip] . s) + offset`, `:632`), con `s` un vector
   unitario fijo. No existe ninguna cantidad de longitud de camino en producción: un `grep` de
   `arclen|path_len|path length|curve length|cumsum` sobre `src/aeroelast/solvers/` y
   `crates/aeroelast-solvers/src/` no devuelve nada, y el bucle rotor en Rust no recibe ningún
   campo de radios.
2. El camino rígido hace **cortocircuito** antes de medir geometría deformada alguna
   (`:910-912`, `disp_max < 1e-12` → solver y proyector de referencia preconstruidos), así que
   el guard del camino rígido compara el cortocircuito consigo mismo.
3. El datum de la malla se ancla en un solo punto (`:402`,
   `_mesh_datum_offset = _ref_r[0] - _mesh_span.min()`).

Medición (fixture `tests/validation/blade/test_blade_deloading_vs_reference.py`, HEAD
2026-10-07: `BEMFSIParticipant` de producción, malla IEA-15MW real y
`BladeAero` real de AeroDyn en el punto rated, un solo solve de cáscara bajo las cargas
proyectadas del participante: **5 passed**). Todas las filas contra la misma baseline de construcción
`_rebuild_bem_solver(_ref_r, _ref_twist)` = 2.541662 MN / 16.389372 MW:

| fila alimentada al BEM | d thrust | d power | tip dr [m] |
| --- | ---: | ---: | ---: |
| radios, solo sesgo de definición (desplazamiento cero) | **-0.38 %** | **-0.17 %** | **-0.5969** |
| radios, solo incremento por deformación (métrica de proyección) | **+1.89 %** | **+1.86 %** | **+1.1062** |
| radios, solo incremento por deformación (métrica de camino) | **+3.63 %** | **+3.56 %** | **+2.1121** |
| radios como producción los alimenta (sesgo + incremento) | **+1.24 %** | **+0.95 %** | **+0.5093** |

Los deltas de las tres filas de radios son tres evaluaciones BEM separadas y **no son
aditivos** (`-0.38 % + 1.89 % != +1.24 %`).

**Conclusión**: el artefacto aporta **cero** al `+1.24 %` real. Si producción midiera un path,
la fila daría `+3.63 % / +3.56 %`; el incremento por deformación medido con la métrica de
producción es `+1.89 % / +1.86 %`, con `dr_proj` = `+1.106156 m` contra `dr_path` =
`+2.112061 m` (`1.909x`, término transversal de segundo orden `47.6 %` del incremento de
camino; longitudes `L_ref` = 116.181290 m, `L_def` = 118.293351 m, `+1.8179 %`).

**Identidad pineada** (assert a `1e-9`): `dr_proj[k] == mean(displacements[strip k] .
span_dir)`. La respuesta de `r_def` a la deformación es el desplazamiento axial del propio
strip, a primer orden: una respuesta axial de un solve lineal, no un artefacto de camino.

**Defecto nuevo que la reemplaza**: `radii_datum_definition_bias` (entra en
`docs/validation/gaps.yaml`, y `force_projection_axial_extension` sale de ese archivo). El
camino rígido usa los radios de estación del deck; el deformado reconstruye las estaciones con
la proyección de span de la malla (`:632`) y su datum anclado en un solo punto (`:402`). A
desplazamiento cero las dos definiciones difieren `+0.397959 m` en el strip 0 y `-0.596870 m`
en el strip de punta, de modo que la primera iteración deformada mueve cada radio de estación
hasta 0.6 m por sí sola y la fila de radios es la **suma neta** de un sesgo de definición
`-0.38 % / -0.17 %` y un incremento por deformación `+1.89 % / +1.86 %`. Ningún guard lo ve:
el guard del camino rígido compara el cortocircuito consigo mismo. Citar una cifra de la
familia de radios como "la recarga geométrica" sobreestima la parte de la deformación; el
resto de los números de la tabla conservan su significado, y los números acoplados pre-fix
siguen sin ser citables.

**Disparador de re-ejecución:** invalida este cierre cualquier cambio en
`_compute_deformed_geometry`, en `_mesh_datum_offset`, en el cortocircuito rígido
(`fsi_participant.py:910-912`), en la asignación de nodos por strip, o en el fixture
`tests/validation/blade/test_blade_deloading_vs_reference.py`.

---

### S-7: la tabla de convergencia contradecía a la medición, y el lado "más blanda" es el artefacto (2026-10-07)

Cierra el issue #20 (ítem 8 del roadmap de #18). S-7
(`tests/validation/blade/test_iea15mw_s7_torsion.py`) sostenía las dos direcciones a la vez: la
aserción medía la cáscara **más rígida** que el deck BeamDyn y la tabla de convergencia
commiteada en el mismo archivo decía **más blanda**. El roadmap lo puso antes de los ítems 2
(#14) y 3 (#15) porque la dirección de la rigidez torsional decide la cota estructural del
`1.61x` contra Zhou.

**La definición no era el problema.** El test y `tools/run_s7_torsion.py` calculan lo mismo:
banda `0.3·117 … 0.9·117`, `theta_s[mid]/theta_beam[mid]`,
`section_twists_deg(n_slices=40)` y `beam_twist_profile_deg`. La expresión del ratio no cambió
desde que el módulo se creó (`5234b48`, 2026-09-23).

**La tabla no se reproduce en HEAD, con ninguna convención de winding.** Registro medido el
2026-10-07 (misma malla, mismo `MOMENT_NM = 1000 N·m`):

| es | nodos | HEAD | HEAD + winding canonicalizado | tabla commiteada |
| ---: | ---: | ---: | ---: | ---: |
| 2.000 | 1460 | 1.0237 | 1.1530 | 1.386 |
| 1.000 | 3040 | 1.0219 | 1.1640 | 1.406 |
| 0.500 | 9271 | 0.9790 | 1.0763 | 1.303 |
| 0.250 | 32325 | **0.939603** | 1.014223 | 1.286 |
| 0.125 | 120352 | **0.933641** | — | 1.273 |

La columna "canonicalizado" restaura el `canonicalize_windings(mesh_model, span_axis=2)` que
`22d3ccc` (2026-10-04) quitó de `BladeMesh.generate`: voltea **1732 de 33462** elementos y mueve
el ratio `0.939603 → 1.014223`, factor **1.079** — la misma dirección que la tabla, pero solo un
cuarto del camino a `1.286`. ⇒ La tabla es un registro irrecuperable de la ventana
**2026-09-30 → 2026-10-04** (merge de `origin/main`, `22d3ccc`, y los cambios de geometría nuMAD
`680cf81`/`6b2cbd6`). El residuo no está bisecado.

**El fixture sí converge, y la dirección es una sola.** `0.250` queda **0.64 %** por encima de
`0.125`, así que `ELEMENT_SIZE = 0.25` se sostiene (el costo de `0.125` es 42 s de malla +
230 s de solve y **12 GB de RSS pico**). La dirección "más rígida" la sostienen cuatro fuentes
independientes (este test, `1/0.939603 = 1.064`; el `GJ` seccional de S-1, `+23.9 %`;
`docs/model_parity_audit.md`, `0.84`; y el `GJ = 1.080` pre-merge, `1/1.080 = 0.926`) y solo la
tabla disiente.

**Qué se quemó:** la tabla de convergencia del módulo (reemplazada por las filas medidas
arriba), el "expected 0.84" del docstring, y la fila de
`docs/origin_main_integration_2026-09-30.md:156`, que además comparaba un ratio de `GJ` contra un
ratio de twist (recíprocos). La cota `_RATIO_TOL = 0.30` **no** se tocó: no es el defecto y el
`6.04 %` medido la pasa por 24 puntos. La aserción del par aplicado sigue sin referente
independiente y va declarada como tal en el store.

**Disparador de re-ejecución:** invalida este cierre cualquier cambio en
`crates/aeroelast-core/src/elements/`, en el generador de malla de la pala
(`BladeMesh.generate`, sus opciones de spacing y su convención de winding), en
`tools/beam_reference.py::load_beamdyn_blade`, en el deck
`tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat`, o en
`tools/run_s7_torsion.py::section_twists_deg` (la construcción del twist de sección).

---

### Los dos defectos de signo de la carga torsional rated (2026-10-07)

Arbitrados al abrir #14 (roadmap item 2): el tool que compara contra Zhou imprimía nuestra
cáscara y nuestra viga con signos opuestos bajo la misma carga. La convención queda fijada
**midiendo** dónde está el borde de ataque, no eligiendo: `tools/diagnose_leading_edge.py` lo da
en el extremo de **x alto** en 8 de las 9 estaciones muestreadas (dos métodos independientes
coinciden; la que difiere es el anillo degenerado de raíz), así que
con el LE en `+x` y el empuje del marco de carga (downwind) en `+y`, **nariz-abajo es `omega > 0`**.
Lo corroboran el signo del de-loading (`-25.31 % / -14.77 %`) y las dos realizaciones rated
(`+8.1048` mínima norma, `+9.6669` multi-celda con propiedades).

| defecto | dónde | medición | arreglo |
| --- | --- | --- | --- |
| `Mp` se aplicaba con el signo del BEM, no del frame | `_section_couples` y el término `Mp` de `_rated_load_cases` | la misma `Mp` nariz-abajo (integral `-1.0029e6 N.m`): `+29.0724` por el proyector de producción contra `-38.6420` por `_section_couples` | `M_z = -Mp`, con un check permanente entre caminos |
| el centro aerodinámico estaba a 0.25 c del **borde de fuga** | `_rated_load_cases`, `x_ac = xs.min() + 0.25 c` | el par de transferencia salía `-1.302 … -0.130` donde debe ser `+1.288 … +0.120` | `x_ac = xs.max() - 0.25 c` |

**El árbitro del signo fue el proyector de producción, y el lado defectuoso era el test-local.**
La cadena del proyector la fija geometría medida fuera del módulo (LE en `+x`, empuje aguas abajo
en `+y`, `+z` mueve el borde de ataque aguas abajo: `tools/diagnose_leading_edge.py`), así que el
signo que carga `_strip_moment_axis_sign` es el arbitrado; `_section_couples` y el término `Mp` de
`_rated_load_cases` —las construcciones **test-local**— se corrigieron para coincidir con él. El
test cruzado `test_the_two_moment_applications_agree_in_sign` lo dice ahora explícito en su
docstring.

Consecuencia medida tras el arreglo: `mp_only` `+4.8337` (era `-4.8337`), `at_ac` `+15.4921`
(era `-20.7348`), `ratio_omega` `2.7327` y `ratio_theta_z` `4.3034` contra Zhou (positivos y en
la convención del frame), `spread_omega` `6.170` (el spread de las cuatro aplicaciones
**test-local** de `_rated_load_cases`, no una propiedad del camino de producción: a este lo acota
el grupo 35), `distortion[at_ac]/distortion[mp_only]`
`18.367`. La tabla del anchor beam pasa de ratios negativos a `ring/beam 4.5622`. Los dos
invariantes de resultantes de `at_ac` se re-midieron (0.0767 %→0.1676 % y 0.3623 %→0.1933 %,
cota 0.5 %). Dos asertos que afirmaban "toda aplicación es nariz-abajo" se acotaron: las
aplicaciones que ponen la fuerza en el centroide del perímetro dan el sentido opuesto, que es la
diferencia de línea de acción que el caso `at_ac` existe para medir.

**Disparador de re-ejecución:** invalida este cierre cualquier cambio en
`_section_couples`, en `_rated_load_cases` (su `Mp`, su `x_ac`, sus pesos tributarios), en
`_ring_kinematics`, en `ForceProjector.project`, en el generador de malla de la pala (la
posición del borde de ataque) o en `tools/diagnose_leading_edge.py`.

---

### #14: el twist contra Zhou 2025 no es transferible, y el árbitro que lo sustituye (2026-10-08)

Cierra el issue #14 (ítem 2 del roadmap de #18). El issue pedía una de tres cosas: (a) obtener los
tres inputs que Zhou no publica (pitch, fuente de rigidez, definición exacta de su torsión) y
rehacer la comparación; (b) cerrarlo como **no-transferibilidad documentada** y sacar el claim de
magnitud del paper; (c) comparar contra una referencia **cuyo modelo sí esté publicado**, p. ej. un
shell CalculiX S8R de la misma pala con las mismas cargas. **Se cierra por (b), con (c) ya
satisfecho.**

**(a) no es posible.** El texto de Zhou no declara pitch en ninguna parte ni la fuente de su rigidez
seccional; su cantidad torsional es "rotation of the airfoil section about the reference axis ...
positive toward stall", que plausiblemente no es nuestra rotación de sección.

**(c) está medido y registrado.** `tools/ccx_blade_twist_arbitration.py` alimenta CalculiX S8R con
**el mismo vector nodal** que produce el `ForceProjector` de producción, con el mismo empotramiento,
y converge cada código en su propia secuencia de malla (comparar MITC4 contra S8R en la *misma*
malla es inválido: conflaciona orden de elemento con tamaño):

| h [m] | AeroElast (MITC4) | CalculiX (S8R) | diferencia |
| ---: | ---: | ---: | ---: |
| 1.00 | +9.6669 | +11.9649 | 23.77 % |
| 0.50 | +8.3937 | +10.7368 | 27.91 % |
| **0.25** | **+8.0980** | **+8.4396** | **4.22 %** |

⇒ La respuesta estructural de nuestro shell, bajo una carga dada, la confirma un código
independiente al **4.22 %**. Ese es el árbitro que sustituye a la comparación con el paper.

**Por qué la comparación con Zhou no es transferible.** Zhou es una **viga** (LL-FVW + GEBT); lo
nuestro es un **shell**, y el rated de producción lleva una distorsión de sección que la viga no
puede representar (`distortion / |omega| = 1.6728`). El par arbitrable es **viga-viga**: bajo las
cargas de su Fig. 11 nuestra viga da `-2.0790°` contra su `-3.6000°`, un **`0.58x`** (déficit de
torsión, no exceso), medido en `tools/diagnose_zhou_loads_reverse.py`. El `1.61x` que el issue cita
es la figura del *shell* de ese mismo tool con cargas ajenas, **no** el rated de producción, y no
debe citarse como tal.

**Números vivos del camino de producción** (2026-10-08, los tres tests de producción del módulo
rated): el camino que usan las campañas acopladas (`BEMFSIParticipant`, sin `element_properties` ⇒
realización de mínima norma) da `omega = +8.1048°`; `standalone.py` con `elements.properties` (flujo
multi-celda) da `+9.6669°`. Contra Zhou eso es `2.2513x` / `2.6853x` **de una magnitud no
transferible**. El patrón de la aplicación está acotado aparte: el tubo cerrado por el
`ForceProjector` de producción reproduce la respuesta discreta exacta al `0.22 % / 0.09 %` (store
grupo **35**, `tests/validation/parity/test_thin_walled_tube_projection.py`).

**El paper no lleva el claim.** `docs/article_draft_wind_energy.md` compara contra Zhou flapwise
(12.78 vs 13.86 m), potencia (14.71 vs 14.76 MW, `-0.34 %`), thrust y frecuencias (Tabla 3); **no**
hay ningún claim de magnitud torsional contra él.

**Asignación de árbitros** (lo que hace citable cada número): `odd/tasks/production-path-and-independent-arbiters.md`,
sección "P2c" — patrón de aplicación → Bredt en el tubo (grupo 35); rigidez seccional → OpenFAST MBC3
/ BeamDyn (S-7, grupo 34); respuesta estructural → CalculiX S8R (4.22 %); respuesta acoplada →
OpenFAST flexible (pendiente, P2b); comparabilidad con el paper → Zhou no decide nada de nuestro
modelo.

**Disparador de re-ejecución:** invalida este cierre cualquier cambio en
`src/aeroelast/solvers/bem/force_projection.py` (el frame de carga, `_strip_moment_axis_sign`, la
realización multi-celda), en `tools/ccx_blade_twist_arbitration.py` o en el escritor de decks CCX, en
el generador de malla de la pala, o la publicación de los tres inputs por parte de Zhou.

### #15: el `+31%` de carga tangencial es frame y sesgo de ángulo, no un defecto de polar (2026-10-08)

Cierra el issue #15 (ítem 3 del roadmap de #18). El issue pedía reproducir el integral tangencial
desde el ángulo de ataque publicado y los polares oficiales —igual que se había hecho con el momento—
y comparar término por término: proyección polar `Cd`/`Cl`, factor de inducción tangencial, y la
convención de signo/dirección del eje tangencial. **Se cierra con los tres términos medidos**, y el
resultado es que no es un defecto de polar: es el frame de la referencia más el sesgo de ángulo de
ataque que el propio reporte ya lleva en §5.11.

**Qué es el `+31%`, y qué no es.** No es un integral: es un estadístico pointwise de `Tp` (el pico
`1.05` contra `0.80 kN/m`, y la media pointwise `1.300` en la única campaña que sobrevive, que es de
baja carga). Los cocientes *integrales* dependen de qué se elija: `1.24` (la misma campaña, `t ∈
[40, 100] s`), `1.49` (derivado de los números de §5.9, no reproducibles), `1.83` (el BEM rígido de
producción). El título del issue presenta una desviación pointwise como un integral.

**Término 1 — la proyección polar: descartada.** El polar es el mismo (el deck oficial) en las dos
lecturas; entra en ambas idénticamente y no puede por sí solo producir una diferencia con forma de
frame.

**Término 3 — la convención del eje tangencial: real y medida.** Las dos curvas no están en el mismo
frame:

- Zhou Fig. 11 es un par del **frame de sección (cuerda)**: leído con
  `alpha = atan2(Tp, Np) + atan2(Cd, Cl)` reproduce su propia Fig. 10 a **0.37°** en `r/R` 0.26–0.80
  (rms `1.31°` en todo el span digitalizado, peor en la punta digitalizada); en el plano del rotor
  falla por **6.83°** en el núcleo.
- `BEMSolver` emite el par del **plano del rotor** por construcción:
  `atan2(Tp, Np) + atan2(Cd, Cl) − twist == alpha` a **1.4e-14°**.

La rotación entre los dos es el twist local. Su efecto neto sobre el cociente medio pointwise es
chico (`1.82 → 1.69`, −7%) sobre todo el span digitalizado (`r/R` 0.20–0.98) y `1.79 → 1.51` (−15%)
sobre el núcleo `r/R` 0.26–0.80, o sea que el −7% está dominado por la punta digitalizada; el twist
cambia de signo a lo largo del span, así que el frame es **fundamental para la definición** pero no
domina la magnitud. Medido por `tools/diagnose_zhou_tp_frame.py`
y guardado por `tests/validation/bem/test_bem_load_frame.py`.

**Término 2 — la inducción tangencial: es el término que queda, y ya estaba reportado.** Medido en el
**mismo frame y con el mismo `qc`**, alimentando el ángulo de ataque de nuestro BEM en la fórmula de
carga de ellos (polar oficial **re-evaluado en cada `alpha`**, `qc` de su propio `|F|`): los cocientes
por componente van de **1.19 a 2.24** en `Tp` y de **1.06 a 1.34** en `Np` sobre `r/R` 0.26–0.80, o
sea que el cociente de dirección `Tp/Np` va de **1.12 a 1.67**. La misma diferencia de ángulo mueve
`Tp` entre 3 y 7 veces más en términos relativos que `Np` (×1.19 vs ×1.06 en 0.26; ×2.24 vs ×1.34 en
0.80). Nuestro ángulo de ataque es mayor que el suyo en
`+0.9°` en `r/R = 0.26`, `+1.2°` a 0.5 y `+4.8°` a la punta. `Tp = qc(CL sin(phi) − CD cos(phi))` es
la **diferencia de dos términos grandes**, así que 1–5° de ángulo de flujo aparecen como decenas de
por ciento en `Tp` y ~1% en `Np`. Ésa es la firma reportada, y es el sesgo BEM-vs-LL-FVW que §5.11 ya
documentaba (`+3.09°` medio).

**Banda del comparador, dicha antes de usarla.** La Fig. 11 digitalizada es el caso **flexible**:
`∫Np·3 = 2.20 MN` = el thrust flexible de la Tabla 6 (`14.76 MW / 2.20 MN`), no el rígido
(`16.11 MW / 2.53 MN`). Y su trapecio **no lleva el torque del propio paper**: `∫(Tp·r)·3·Omega` da
`11.68 MW` contra los `14.76 MW` reportados (−21%). Con eso, ningún claim *integral* se toma a través
de la digitalización; la comparación utilizable es la de dirección (`Tp/Np`) punto a punto. La banda
pointwise de la digitalización es `±0.05 kN/m` (±3–5%) y la de la Fig. 10 `±0.5°`.

**Números del camino de producción** (BEM en el punto rated de Zhou: `V = 10.59 m/s`, `Omega = 7.55 rpm`,
pitch 0, deck AeroDyn oficial): `thrust 2.478 MN`, `power 15.724 MW`, `torque 19.888 MN·m`;
`∫Np dr = 922.2 kN`, `∫Tp dr = 127.8 kN` en el plano del rotor. Contra la curva flexible de ellos la
diferencia de nivel es de caso (rígido-vs-flexible) más el sesgo de ángulo, y por eso no se cita como
magnitud.

**Números anteriores, conservados y no reproducibles.** §5.9 reportaba `∫Np dr = 669.1 kN`,
`∫Tp dr = 101.7 kN`, `∫(Tp·r) = 6.21 MN·m`, pico `Np` 10.65 vs 9.95 y pico `Tp` 1.05 vs 0.80 kN/m.
La campaña que los produjo (`frontiersin_results_corotational`) **ya no existe en disco**; la única
superviviente (`..._100s`) es otra corrida de baja carga (`thrust` 1.64 MN, `power` 12.34 MW, tip
7.4 m contra los 14.71 MW / 12.73–12.79 m del artículo). Quedan como registro, no como evidencia.

**Store.** `docs/validation/gaps.yaml` id `zhou_spanwise_load_frame` (`not_validated`,
`citations_forbidden: true`): la referencia no puede arbitrar una magnitud de `Tp`. V-09 sigue en
Capa B, como comparación contextual. El guard del frame es
`tests/validation/bem/test_bem_load_frame.py`, declarado `out_of_scope` en `groups.yaml` (su sujeto
es en qué frame están dibujadas las dos curvas, no una magnitud física contra una referencia).

**Disparador de re-ejecución:** invalida este cierre cualquier cambio en
`src/aeroelast/solvers/bem/engine.py` (el frame de `Np`/`Tp` que emite `BEMSolver`), en
`tools/diagnose_zhou_tp_frame.py`, en `docs/validation_data/zhou_2025_fig11_loads.csv` o
`..._fig10_aoa.csv` (una re-digitalización), en `tests/support/openfast_bem.py`, o que Zhou publique
la curva tabular, el pitch, la fuente de rigidez o la definición torsional. Reabre además la pregunta
de si `ForceProjector` aplica `Np`/`Tp` en el frame correcto (ver el follow-up anotado en
`odd/tasks/tangential-load-attribution.md`).

### Cierres de tooling del store (2026-10-08)

No son ítems de validación, pero el roadmap los sigue porque un gate roto o destructivo esconde
defectos de datos. Se registran acá para que el contrato de cierre valga también para ellos.

**#21 — `regression` no honraba `non_validation_tests` (cerrado 2026-10-08).** El comando armaba su
lista de sobrantes comparando los nodos colectados contra `store.rows` y **nunca** consultaba las
declaraciones del grupo (`tools/validation_matrix.py`), así que un grupo con las declaraciones
correctas se reportaba `unclaimed` — y `unclaimed` está en el set que falla, o sea salía non-zero
justamente en los grupos que habían declarado bien. Arreglado en `2c309ed`: el split pasa por el
mismo `classify_unclaimed` que usa `extract` (`validation_matrix.py:2902`), así que los dos comandos
no pueden discrepar. Verdicts: `declared_non_validation` (se imprime, no falla), `unclaimed` (sigue
fallando), `stale_declaration` (se muestra); `REGRESSION_FAILING_VERDICTS` hace el contrato
testeable, y el test unitario se escribió test-first (falló con `AttributeError` antes de existir el
helper). Aceptación medida sin `--write`: `--group 34` pasa de `1 same, 1 unclaimed` (exit != 0) a
**exit 0** con `1 declared_non_validation, 1 same`; `--group 30` pasa de `5 same, 29 unclaimed` a
**exit 0** con `29 declared_non_validation, 5 same`; la suite de tools da `121 passed` con
`-m "not slow"` y `ruff` limpio.

**El defecto que el arreglo reveló**: `test_compute_performance_coefficients_uses_given_radius`
(grupo 30) estaba contado entre esos 29 falsos `unclaimed`, así que su falta de declaración era
invisible. Es una propiedad (desigualdades, sin referencia independiente) ⇒ se declara, no se le da
fila; ahora está en `docs/validation/groups.yaml`. Un defecto de herramienta puede **enmascarar** un
defecto de datos, y arreglar la herramienta es lo que lo expone.

**Abierto con registro: #24 — `coherence` es destructivo.** `test_every_group_re_derives_to_the_rows_on_disk`
corre `coherence` **con write** por diseño y esa corrida reescribe los archivos de fila enteros,
borrando la prosa escrita a mano y los bloques `measured` (grupos 6, 20, 34 y 35 el 2026-10-08).
Mientras siga abierto: **correr `tools/tests` con `-m "not slow"` y verificar `git status --short`
después de cualquier corrida de tools.** El detalle está en
`odd/tasks/production-path-and-independent-arbiters.md`.


### #19 y #26: la proyección conserva la geometría de referencia, y el par del BEM rueda en el plano del rotor (2026-10-08)

Los dos ítems P0 del roadmap se cierran juntos porque comparten archivo, aunque son independientes:
el #26 es un error de *frame* con la lectura arbitrada desde el fuente de `ccblade`, y el #19 era la
no-contracción del acople implícito.

**#26 — `Np`/`Tp` ruedan en el plano del rotor, no sobre los ejes de la sección.** `ccblade`
descompone la fuerza de sección **en el plano del rotor** (`cn = cl·cos(phi) + cd·sin(phi)` con
`phi` medido desde ese plano), así que `Np` va por el normal del eje del rotor y `Tp` por la
tangencial en el plano — un solo frame para la pala entera. `9a3923e` montaba el par sobre los ejes
de cada sección, que son ese frame **rotado por el twist local**, así que cada carga aplicada era el
vector físico rotado por `theta`. Medido en la malla real y el deck real a rated, con las dos
aplicaciones descompuestas en el mismo frame y el brazo `blade_aero.r`: in-plane `+11.36%`, eje del
rotor `-0.83%`, momento edgewise de raíz `-3.97%` contra la aplicación correcta. Por estación el
normal de sección está a `22.245°` del normal del rotor en `frac 0.03` (39.18% de `|F|` sobre la
tangencial) y baja a `1.176°` en la punta; la dispersión entre estaciones es `24.410°`. La identidad
sobre nuestra propia salida, `atan2(Tp,Np) + atan2(cd,cl) - twist = alpha`, cierra a
`max |residuo| = 1.4e-14°` sobre las 50 estaciones (`tools/diagnose_zhou_tp_frame.py`).
**Arreglado en `7a84da1`**: la carga se aplica sobre `_rotor_normal_dir`/`_rotor_tangential_dir`,
construidos una vez desde el par configurado con `_load_frame`; el frame de sección queda como datum
del centro aerodinámico y del eje del pitching. El guard es
`tests/validation/bem/test_force_projection_load_frame.py`, reescrito al reading del rotor (**4 de
sus 5 guards fallan antes del cambio**) y con la malla twisted como caso de respuesta conocida: el
tubo sin twist del grupo 35 no puede discriminar los dos frames. 122 tests verdes en el conjunto
afectado (bem 58; projector + ac datum + load frame 22; blade + parity + rotor + multicell 42).
**Re-ejecutar si cambia**: `solvers/bem/force_projection.py`, los pares de dirección configurados, o
el deck.

**#19 — la no-contración era la geometría de la proyección siguiendo la deformación.** El
participante BEM reconstruía el `ForceProjector` sobre la malla deformada en **cada sub-iteración**,
así que la grilla nodo→franja, el brazo al AC y el eje del momento de cabeceo se movían con la
deformación mientras las cargas también. Esa ganancia extra lleva al acople particionado (IQN-ILS)
por encima de su umbral de contracción. La escalera, sobre el caso propio de la campaña a
`max-time 5.0`, con el Rust constante y el caso consenso `tests/smoke_fix/frame_ab_abs/`:

| árbol / experimento | ventanas | convergidas |
| --- | ---: | ---: |
| `campaign` = `b5d369e` (lado local del merge) | 500 | 500 (100%) |
| `origmain` = `5f22f51` (2º padre del merge) | 206 | 185 (89.8%) |
| `ac2e9e8` / `fee690` = `4fee690` (lado origin/main) | 243 / 217 | 232 (95.5%) / 192 (88.5%) |
| **HEAD** (`7a84da1`, con el fix del #26) | 500 | 39 (7.8%) |
| `headfix` / `nofeed` (radios+twist congelados) | 159 / 500 | 14 (8.8%) / 144 (28.8%) |
| `projfrozen` (projector de referencia) | 94 | **94 (100%)** |
| `bemhead` (árbol `origmain` + fluido de HEAD) | 24 | 2 (8.3%) |
| `frzsign` / `frzac` / `frzgrid` / `frzrings` / `frzdist` | — | 17% / 44% / 11% / 4% / 9.5% |
| `frzall` (signo + brazo + grilla + anillos) | 39 | **39 (100%)** |
| **`fixfull` (el fix en el árbol)** | 500 | **500 (100%)**, media 2.60, máx 11 |

Dos lecturas valen más que el número final: **los dos lados del merge son sanos por separado**
(`b5d369e` 100% y `5f22f51` 89.8%), así que la regresión es la **combinación** que armó `26ffe6e`,
no un commit de ninguno de los dos; y **ninguna salida geométrica sola es la palanca** — fijar el
signo del momento, el brazo al AC, la grilla de franjas o los grupos/anillos deja la contracción
rota, mientras fijarlos todos la restaura. El camino de anillos ni siquiera está activo en el fluido
(no hay `element_properties`, los anillos quedan no-físicos y se cae a `_distribute`), lo que
confirma que el efecto es la suma y no una pieza discreta. De paso, el fix deja el acople mejor que
el baseline sano: **media 2.60 sub-iteraciones por ventana contra 4.08 de la campaña**.
**Arreglado en `6765633`**: `_compute_forces` usa el projector de referencia; la geometría del BEM
sigue la deformación intacta (`_compute_deformed_geometry` y `_rebuild_bem_solver` sin tocar) y
`_rebuild_projector` queda sin llamadores, así que se elimina. Guard `e3278ba`
(`test_force_projection_geometry_stays_on_the_reference`): para el mismo resultado BEM **deformado** —
que se assertea distinto del rígido, para que se vea que el feedback corre — las fuerzas nodales
tienen que ser las del projector de referencia. **RED contra `7a84da1`: 9129/9129 elementos
difieren, diferencia absoluta máxima 1372.75 N**; verde con el fix. 38 tests verdes (participante +
multicell 23; de-loading + rated-twist 15).
**Contra conocida, para el revisor**: la distribución *dentro* de la franja ya no ve los offsets
deformados — escala de cuerda contra una deformación de escala de pala. El projector de la propia
campaña *sí* seguía la deformación y contraía, así que el camino fino queda abierto: hacer la
derivación del lote #11 robusta a la deformación (selección de anillo sin switch, tie-break LE/TE y
continuidad del eje de cuerda sin histéresis) en vez de congelarla. Es fidelidad, no este cierre.
**Re-ejecutar si cambia**: `solvers/bem/fsi_participant.py` (`_compute_forces`), el projector, o el
esquema de acople en `precice-config.xml`. La verificación larga es el gate de 30 s
(`tests/run_step1b_smoke.srm`, job `11611263`) y, después, el re-anclaje de la campaña.

### #27: la mitad que faltaba éramos nosotros — la recuperación de tensiones de compuestos no daba tensión de ply (2026-10-09)

Cierra el issue #27 (ítem 4 del roadmap de #18, P1, sucesor de #3). El issue asumía que la referencia
no existía porque `*SHELL SECTION, COMPOSITE` ignora `OUTPUT=3D`. Medido, las dos mitades de esa
premisa son falsas y en direcciones opuestas: CalculiX **sí** juzga la fibra externa (su valor en el
FRD sigue la ply de cierre —3.6820/3.6853 MPa en `[0/90]s`, 0.3150/0.3148 MPa en `[90/0]s` a 8x2/16x4—
mientras el nodo de midsurface da exactamente 0.00000 MPa: la tarjeta expande la sección), y
`OUTPUT=3D` es inerte para compuestos (idéntico valor a valor en el deck de la fila; solo desaparecen
los 233 nodos de referencia sin tensión). Lo que faltaba era nuestro: la recuperación devolvía la
media homogeneizada de espesor en TOP/MIDDLE/BOTTOM (2.0000 MPa en los dos stacks, contra
3.696 / 0.304 MPa de la fibra real: −45.9 % / +558 %), con el `Cm()/h` de
`src/aeroelast/postprocess/stress_recovery.py` como mecanismo, no un error de escala.

El fix (`208f220`) resuelve la ply que contiene `z` y evalúa `σ_ply = Qbar(θ_ply)·ε(z)`; el camino
isótropo queda bit-idéntico (0 de 108 arrays distintos) y un dict ABD crudo conserva la media
avisando una vez por elemento set. El defecto *smeared* que el ítem expuso se publicó como issue
propio (**#28**, cerrado el mismo día) en vez de absorberse en #27, según la regla del roadmap.

**Árbitro doble, en una fila.** CalculiX 2.20 S8R con `*SHELL SECTION, COMPOSITE` (`kind="code"`) y
una CLT reimplementada dentro del test (`kind="analytical"`), convergido-vs-convergido con malla
propia por código, leído en **campo libre** porque el máximo global es una singularidad de
introducción de carga en el **borde libre cargado** (5.8838 → 6.5803 MPa al refinar; el
empotramiento no es singular, 3.69 → 3.71 MPa). Medido a 16x4: nuestro 3.6959 / 0.3041 MPa, la forma
cerrada igual a ≤0.0006 %, CalculiX 3.6853 / 0.3148 MPa → 0.29 % / 3.41 %.

**Dos cosas para llevarse, no para festejar.** `TOL_CCX = 0.05` es **cota, no convergencia**: la
brecha de `[90/0]s` casi no se mueve al refinar (3.48 % → 3.41 %) y todo el gap está del lado de
CalculiX. Y este deck reproduce el **mecanismo** de la sonda E2 pero no sus valores (0.07 % / 1.6 %):
ninguna elección de nodo sobre él da 3.6984/0.2991 MPa, o sea que la sonda usó otro deck, ya
irrecuperable (`/tmp/ccx_composite_probe/` no existe). Queda registrado como **ancla movida**.

**Store.** Grupo 36 `composite_ply_stress` (13 filas / 18 comparaciones, 12 medidas; 3 filas sin medir
con motivo escrito: el extractor deduplica aserciones dentro de loops). `gaps.yaml`
`composite_stress_recovery` **afinado, no borrado**: solo quedan los dos residuales sin evidencia
citable — un dict ABD crudo conserva la media de espesor, y el caso de flexión `[0/90/90/0]`
(616.10 vs 595.59 MPa) no tiene fila. Los dos defectos de herramienta que aparecieron corriendo la
gate de este ítem se publicaron como **#29**.

**Re-ejecutar si cambia**: `src/aeroelast/postprocess/stress_recovery.py` (stack normalizado o
`Qbar`), el escritor de decks compuestos (`core/mesh/io/writers.py`), o la versión de CalculiX
(aquí 2.20). Commits: `95bfbab` (RED), `208f220` (fix), `bd4cc05` (cobertura), `62d806a` (fila),
`5260233` (store).

### #16: la realización wall-flow llega a producción y queda **apagada por defecto** — su efecto acoplado no se pudo atribuir (2026-10-09)

Ítem P2 del roadmap de #18, rama A aprobada por el usuario ("activar skin wall-flow + medir el
movimiento"). El ítem se **cierra a medias a propósito**: el defecto que denunciaba está arreglado y
medido, pero la activación queda opt-in porque el movimiento que produce no tiene árbitro. #16 queda
**abierto y bloqueado** por el issue de seguimiento, según la regla del roadmap (un ítem bloqueado lo
dice en su propio comentario y no se mueve).

**El defecto era real y era de la malla, no de la realización.** `_build_mesh` filtraba la malla a
`allOuterShellNods` con `MeshModel(nodes=...)` y tiraba **todos** los elementos, así que `ring_section`
no encontraba ninguna arista de pared: **0 de 671** anillos usables a `element_size: 0.25`, y cada
strip caía al mínimo-norma. Se arregla en dos piezas: `MeshModel.subset_to_nodes` (nodos exactos del
caller, ids preservados, elementos totalmente contenidos, sets restringidos) y el mapa de propiedades
del deck restringido a los sets que la malla conserva (`_element_properties_for`, que lo **retiene y
avisa** si algún elemento de acople queda sin cobertura). Medido: **671 de 671** anillos realizables,
todos de una celda en la malla de campaña, `S = G·t` resuelto desde los laminados, cobertura 100 %.

**Un segundo hop, y lo cazé yo mismo.** Mi primer A/B salió **idéntico a 10 dígitos** (rel. L2
2.9e-15 en el campo de fuerzas): las campañas no usan `aeroelast-bem-fsi` sino `aeroelast`, que
despacha los configs con `bem:` a `run_aero_fsi` → `build_aero_participant_from_config` (participante
legacy). La malla llevaba los elementos (está en el log), pero el mapa nunca llegaba al projector.
Arreglado en el mismo ítem, con test que falla si se cae el hop.

**El movimiento acoplado (caso gate a `max-time 5.0`, A = fluido `bd8237d`, B = activado, C = malla
con elementos sin mapa, sólido idéntico en las tres):**

| cantidad (media t∈[2.5,5] s) | A | C | B |
| --- | --- | --- | --- |
| empuje | 2.492465 MN | 2.492465 MN | 2.486703 MN |
| torque / potencia | 19.2082 / 15.1867 | 19.2082 / 15.1867 | 19.1036 / 15.1040 |
| tip disp X (edgewise) | -1.445670 m | -1.445670 m | -1.327938 m |
| tip disp Y (flapwise) | 16.324849 m | 16.324849 m | 16.354168 m |
| fuerza nodal máxima | 168.460 N | 168.460 N | 194.999 N |
| von Mises TOP/MID/BOT | 2683 / 954 / 155 MPa | 2683 / 954 / 155 MPa | 2165 / 746 / 187 MPa |
| L2 de la fuerza aero | 9.936e3 N | 9.936e3 N | 1.653e4 N |
| rotación de sección (ROTZ) | +2.583° | +2.583° | **-27.404°** |
| rotación de sección (best-fit) | +1.208° | +1.208° | **-9.270°** |

**C ≡ A a 8-9 dígitos y B ≠ ambas**: el mecanismo es la realización y nada más (queda refutada la
hipótesis del feedback de twist por la malla con elementos). La diferencia de twist es *distribuida*
(~0 en la raíz → ~-10° en la punta), con los resultantes por strip idénticos a 1e-15 y el flapwise
casi igual: **no es un refinamiento, es otro estado físico**.

**El árbitro que se pidió antes de decidir, y que empató.** Momento torsor puro de 1e6 N·m sobre un
segmento skin-only de la malla de producción (empotrado hasta z = 40.95 m, cargado en z = 63.28 m,
L = 22.33 m; `$SCRATCH/bfs16/probe_section_gj.py`): ambas vías aplican el mismo resultante exacto
(1.000000e6 N·m) y el giro de sección queda **+5.002° (wall-flow) vs +5.454° (mínimo-norma) contra
+8.582° de Bredt** del propio anillo (`multi_cell_shear_flow` con sus `S = G·t`). 9 % entre ellas,
~40 % las dos por debajo de la referencia: **sin veredicto**. Dato extra que corrige al store: los
anillos no son todos de una celda — a `0.25` los 671 son de una celda, a `1.0` hay 94 de una, 4 de
dos y 88 de tres, así que los tests en malla gruesa **sí** ejercitan la rama multi-celda.

**Decisión del usuario (2026-10-09): parquear.** `bem.wall_flow_realisation` (nuevo, default `false`)
gatea la entrega del mapa en los **dos** puntos de entrada de producción; el camino de producción y
toda cifra de campaña siguen sobre el mínimo-norma y el baseline viejo sigue en pie. El movimiento no
atribuido se publica como issue propio (con las tres mediciones) y #16 queda bloqueado por él.

**Store.** Grupo nuevo **37** `wall_flow_activation` (7 tests: los 6 guards de activación +
el default opt-in), 0 filas y `non_validation_tests` declarados con motivo. `gaps.yaml`
`moment_realization_over_delivers` sigue **acotado** y ahora dice explícitamente que la realización
existe pero está apagada, que encenderla **no** es un refinamiento, y que una corrida activada debe
declararlo y no puede llamar validado su twist. `check`: 211 filas / 277 comparaciones / 0 errores / 0
warnings; `status`: 0 archivos sin agrupar. Resultado nulo registrado: la tabla one-way contra
Zhou et al. 2025 **no** se re-corrió (la producción acopla sólo el skin; esa tabla usa la malla
completa), así que la predicción del doc de feature sigue sin confirmar.

**Re-ejecutar si cambia**: `_effective_element_properties` / `bem.wall_flow_realisation`,
`MeshModel.subset_to_nodes`, `_element_properties_for`, o el camino
`run_aero_fsi` → `build_aero_participant_from_config`. Commits: `a0b9a37` (RED), `fe6ca30`
(`subset_to_nodes`), `756a3b2` (elementos en la malla), `81a1e15` (mapa al projector), `790b80e`
(el hop del dispatcher), `2028dc8` (store), y el commit de parqueo (config key + guard).

### #30: el movimiento acoplado del wall-flow está atribuido — es respuesta estructural al patrón, no la carga (2026-10-10)

Ítem P2 del roadmap de #18. #30 nacía de la parqueada de #16: con la realización wall-flow activada,
la rotación de sección de punta del caso gate se movía `+2.583 → -27.404°` (`ROTZ`) / `+1.208 →
-9.270°` (best-fit) **con el resultante por strip idéntico**, y nada lo atribuía. Se cierra con
cuatro mediciones, bajo la misma regla con la que se escribió el ítem: cada hipótesis primero recibe
una sonda estática o de un solo participante, y **el A/B acoplado de 5 s (5 h de cola) no se volvió a
correr**.

**T1 — el Transfer no es el mecanismo (resultado nulo).** Las tres corridas imprimen el mismo bloque
de preCICE: `nearest-neighbor`, `Mapping distance min:0 max:0 avg: 0 var: 0 cnt: 27609` en **ambas**
direcciones. Verificado por fuera con `meshio`: `fluid_mesh.vtu` y `solid_mesh.vtu` son el mismo
conjunto de 32325 puntos (`max = 0.0`), y no hay conectividad declarada en ningún lado
(`set_mesh_vertices` solamente). Un mapeo de orden superior devolvería el mismo vector salvo
redondeo. El prior del usuario (que `nearest-neighbor` pierde el patrón) sigue siendo válido **donde
los meshes no coinciden** — sus casos OpenFOAM — y eso es otra observación, no #30.

**T2 — el wall flow ES Bredt sobre la sección del propio deck; el que sobre-entrega es el
mínimo-norma.** Sonda `tools/diagnose_section_moment_realization.py`, control = el rectángulo del
grupo 31, tratamiento = el anillo de piel del deck en `z ≈ 40 m` (51 nodos), extruido a `L = 30.268 m`,
ventana interior `(0.4L, 0.9L)`, `T = 1e4 N·m`, referencia Bredt con el `A66` propio de cada pared:

| forma | realización | BC | rate/Bredt | distorsión rms/dim |
| --- | --- | --- | --- | --- |
| rectángulo (control) | wall flow | self-eq | **1.00309** | 1.270e-7 |
| rectángulo | mínimo-norma | self-eq | **31.69331** | 1.0649e-3 |
| rectángulo | wall flow / mínimo-norma | clamped | 0.99156 / **10.40654** | 5.204e-6 / 6.0596e-4 |
| **aire del deck (51 nodos)** | wall flow | self-eq | **0.95943** | 4.1267e-6 |
| aire del deck | mínimo-norma | self-eq | **1.65658** (sube con el conteo de nodos) | 7.4327e-5 |
| aire del deck | wall flow / mínimo-norma | clamped | 0.96089 / 1.56841 | 4.0818e-6 / 7.7050e-5 |

Referencia del aire: `A = 5.0655 m²`, `GJ = 4.623668e8 N·m²`, `θ' = 2.162785e-5 rad/m` (suma armónica
del `A66` por pared, que varía 16.8x en el perímetro). El wall flow queda en Bredt en **las dos**
formas (0.96–1.00, distorsión 4e-6): la lectura "es un artefacto de los DOF de drilling" queda
**refutada en el camino de torsión**. El mínimo-norma sobre-entrega en toda sección no circular. La
hipótesis de elongación de la nota de feature queda refutada: el aire es más elongado que el
rectángulo y diverge mucho menos; el driver es el desajuste radio-vs-tangente en los nodos de
discretización.

**T3b — el par aplicado no invierte el signo; la amplificación es estructural.** Dos sondas, ninguna
de ellas una campaña.

*Parte 1* (`tools/diagnose_wall_flow_pattern_ab.py`, lee los `fields.vtu` de A y B, no resuelve nada).
Primero, dos hechos que hay que verificar antes de tocar estos archivos: `points == solid_mesh.vtu +
U` **index a index** (7e-15) en los dos runs, y `A/solid_mesh.vtu == B/solid_mesh.vtu` bit a bit
(32325 nodos) — o sea que la geometría de referencia es compartida y comparar `points` entre runs es
comparar deformadas (0.73 m de diferencia), no mallas. Y el `F_AERO` guardado está en el frame de
referencia/rotante, no inercial (axial share vs el eje de pala `5.3e-3`; rotarlo por `R(-θ)` sobre
`Y` lo sube a `1.3e-1`).

| cantidad | A (mínimo-norma) | B (wall flow) |
| --- | --- | --- |
| `\|F_AERO\|` L2 | 9935.9 N | 16528.3 N (**1.664x**) |
| campo diferencia `F_B - F_A` | L2 `1.326e4` N = **133.4 % de A**, mueve `27609/27609` nodos | resultante global difiere `0.45 %` |
| momento torsor acumulado (lo que tuerce una sección empotrada) | `5.941e5` N·m | `6.497e5` N·m (**+9.4 %**, mismo signo) |

El exceso se entrega fuera de `z ≈ 78 m`. La respuesta no acompaña: con `1.09x` el par, B gira
`2.6x` más que A (`-15.76` vs `-6.05°`, brazos de referencia, en el anillo de punta), con distorsión
normalizada `7.94e-2` vs `5.35e-2` y RMS de rotación nodal sobre el eje de sección `1.04` vs `0.245`
rad en `z = 109 m`.

*Parte 2* (`tools/diagnose_wall_flow_static_ab.py`): el campo `F_AERO` **congelado** de cada corrida,
aplicado a la malla de producción en un solve **estático lineal** — sin dinámica y sin feedback aero.
La malla se regenera con los parámetros del propio caso y sale **index a index idéntica** al
`solid_mesh.vtu` de los runs (drift `0.0`), así que el campo se mapea nodo a nodo; `RootNodes`
empotrado como en el `dirichlet` del caso. Control: la deflexión estática reproduce la acoplada
(`max|u|` 23.63 / 23.35 m contra 24.10 / 23.63 m). Con el par aplicado en el anillo de punta
diferiendo `+4.5 %`, en la ventana de pala externa `z ∈ [93.6, 112.3] m` la rotación afín de sección
de B es **mediana 1.98x** la de A (rango 1.02–3.20), la rotación rígida **3.31x** (1.41–4.43) y la
distorsión `1.77x`. El `2.6x` acoplado cae adentro de ese rango: **el mecanismo es la respuesta
estructural al patrón y no necesita el lazo acoplado**. Estaciones testigo: `z = 104.52 m` → `2.52x`
(afín) y `3.57x` (rígida).

**Caveat que toca al titular del issue.** El anillo de punta tiene 16 nodos y radio medio `0.176 m`, y
el centroide deformado queda a ~24 m del de referencia, así que cualquier best-fit que mezcle brazos
de referencia con centroide deformado degenera (`-0.0004°`). El mismo punto da `-6.05 / -15.76`
(brazos de referencia), `-5.98 / -15.01` (brazos deformados), `+5.71 / -9.58` y `+2.33 / -3.88`
(outer 5 % del span), y el `+1.208 / -9.270` del issue sale de otra construcción (media del perfil
interpolado de 40 slices). **El valor absoluto de la punta y hasta su signo dependen del estimador**;
el signo del cambio y el DOF crudo (`ROTZ` medio) no. Toda cifra de rotación de sección de este issue
debe nombrar su estimador, y el anillo final no decide ningún ratio (distorsión 0.65–1.0, totalmente
distorsionado).

**T3 — el over-delivery del mínimo-norma queda pineado en el repo.**
`tests/validation/parity/test_thin_walled_tube_moment_realization.py::test_minimum_norm_realisation_over_delivers_bredt_by_the_recorded_factor`
ejecuta el `31.69` de la nota T6 sobre el tubo cerrado validado: `31.69331x` Bredt auto-equilibrado
(`8.719852e-4` rad/m) y `10.40654x` empotrado (`2.863174e-4`), con `distortion/|rotation| = 2.042`
contra `0.008` del wall flow. El `1.65658x` del anillo del deck queda como medición de sonda:
promoverlo pediría construir el tubo extrusionado dentro de un test, y T3 fijaba el tubo cerrado
validado como árbitro, no una construcción nueva.

**Veredicto.** El movimiento acoplado **no** es un cambio de signo del estado de carga ni un artefacto
del Transfer: es la **respuesta estructural al patrón** de la realización. El resultante es el mismo
y el par torsor aplicado es `+9 %` mayor y del mismo signo; lo que cambia es que el patrón wall-flow
— tangencial y concentrado en las paredes — es respondido con `2–3.3x` la rotación de sección
(estático, campo congelado, sin dinámica ni feedback aero) y con más distorsión de sección. El `ROTZ`
crudo que titulaba el issue (`+2.583 → -27.404°`) sí invierte el signo, pero es un DOF **nodal**
(`src/aeroelast/solvers/checkpoint.py:34-35`), no una rotación de sección: es el contenido de
drilling que el patrón tangencial excita `4–10x` más, y es la lectura 2 del issue ganando evidencia
en el camino acoplado — donde T2 la había refutado en el de torsión pura. La lectura 1 se sostiene en
el camino de torsión: el wall flow es Bredt en las dos secciones con referencia y el mínimo-norma
es el que sobre-entrega.

**Lo que no cambia.** Las reglas del gap siguen: una corrida con la realización activada debe
declararlo en sus propios resultados, **no** se puede comparar con una de mínimo-norma como si fuera
el mismo estado, y no puede llamar validado su twist. Ahora tienen mecanismo medido, no solo
magnitud. El A/B de 5 s no se re-corrió; el costo de una campaña activada sigue siendo una decisión
de #16.

**Sin fila nueva en el store, a propósito.** El pin del mínimo-norma es una aserción **desnuda** y no
un `assert_relative_error`: las comparaciones del store necesitan una referencia independiente que la
medición deba *cumplir*, y Bredt es la que este campo **falla**, así que registrarlo vestiría un
defecto de validación. `check`: 211 filas / 277 comparaciones / 0 errores / 0 warnings, sin cambios.
El `provenance_note` del grupo 31 se corrigió: describía una aserción mínimo-norma que no existía en
el archivo.

**Re-ejecutar si cambia**: `ForceProjector._distribute`, `realise_section_load`,
`_realise_multi_cell_section_load`, la malla de acople (`element_size`), o el camino de proyección.
Sondas: `tools/diagnose_section_moment_realization.py` (T2), `tools/diagnose_wall_flow_pattern_ab.py`
(T3b parte 1), `tools/diagnose_wall_flow_static_ab.py` (T3b parte 2). Commits: `dccf9e2`, `a890b30`
(T1/T2), `5b80b7c`, `d376d4c` (T3b-1), `113b135`, `4b5ea82` (T3b-2), `9152e2b`, `40f8978` (T3).
