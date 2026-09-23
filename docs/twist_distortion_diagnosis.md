# Diagnóstico: twist elástico y distorsión de sección del shell

**Fecha**: 2026-09-16
**Estado**: diagnóstico abierto — documentación de los hallazgos y evidencia hasta la fecha.

## 1. Contexto

Tras el fix del signo del twist (`rotorspin = -1` en el path YAML, ver
`docs/validation_theory/` y Engram "Twist sign fix applied and verified"),
el smoke FSI mostró de-loading excesivo: **thrust 1.736 MN = −31.7%** vs el
ancla rígida (2.541 MN, shear 0). El ancla oficial flexible (OpenFAST
BeamDyn, deck CompElast=2 del paquete IEA-15-240-RWT) da **2.159 MN =
−15.0%** de-loading y torque 19.02 MN·m (−8.9%) — consistente con Zhou 2025
(−13.0% thrust, −8.4% power).

El exceso de de-loading del shell (~2×) se origina en el **twist elástico**:
el shell tuerce +8.9° en la punta bajo cargas AD rated, vs +0.98° del ancla
BeamDyn. Este twist alimenta el feedback de torsión del acoplamiento BEM y
produce el de-loading excesivo.

## 2. Evidencia del twist excesivo (mediciones limpias)

Medición limpia: anillo fino (±0.3 m) + rotación del vector LE→TE, cargas AD
distribuidas (S-5 rated), malla 0.5 m:

| z [m] | uy deflexión [m] | pendiente dy/dz [°] | twist shell [°] |
|---|---|---|---|
| 60 | 2.5 | 6.5 | +1.5 |
| 80 | 5.5 | 10.4 | +3.8 |
| 100 | 10.1 | 15.5 | +6.4 |
| 110 | 12.7 | 17.6 | +8.9 |

Ancla BeamDyn: +0.545° mid-span, +0.981° tip (rotación Wiener-Milenkovic
RDxr, estado estacionario t>80 s del run rated).

Convergencia de malla (S-8c a 0.125 m, 120k nodos): mid +1.57°, tip +6.43°
— **persiste con refinamiento**, no es discretización spanwise.

## 3. Causas eliminadas

| Hipótesis | Test | Resultado |
|---|---|---|
| Rigidez torsional (GJ) | S-7 (par de fuerzas en punta) | Shell ≈ beam (±14%) — OK |
| Acoplamiento material del layup | Inspección del YAML | fiber_orientation = 0° en todas las capas — sin acoplamiento material |
| Artefacto de medición (best-fit) | Anillo fino + vector LE→TE | Twist real, no artefacto |
| Refinamiento spanwise | 0.5 → 0.25 → 0.125 m | Persiste (~-13% en mid) |
| Conectividad de webs | Duplicados de nodos coincidentes | 0 pares — malla estanca |
| Brazo fuerza–centro de corte | Fuerza en centroid vs mid-chord | Mismo twist (−5.2 vs −5.5°) — no es el brazo |
| Pendiente de la línea elástica | dy/dz vs twist | Correlacionada pero no 1:1 (ratio 0.4–0.5) |

## 4. Hipótesis activa: distorsión de sección

El desglose de la deformación en el plano de la sección (tras quitar la
traslación rígida y la mejor rotación) muestra que el residuo (distorsión)
es ~98% de la señal: **la sección no rota como cuerpo rígido, se distorsiona
en su plano** (diferencial LE–TE ~0.2 m sobre una cuerda de 1.5 m en la
punta). El vector LE→TE (y el SVD de la cuerda que usa el feedback BEM)
lee esa distorsión como twist.

Candidatos del mecanismo:
1. **Drilling stiffness del MITC3/MITC4**: `k_drill = 0.15·E·t²` — para los
   paneles delgados (t ~ 2–5 mm, E ~ 40 GPa) el penalty es
   ~54–340 N·m, diminuto. La rotación en el plano de los paneles queda casi
   libre → la sección puede distorsionar con energía casi nula.
2. **Paneles del TE/LE muy delgados** con pocos elementos en el contorno.
3. **Respuesta real de sección** (ovalización Brazier + shear lag) que el
   shell captura y el beam no — pero la magnitud (+8.9° ≫ 1–2° físicos)
   sugiere que está sobre-predicha.

## 5. Test de referencia: el cajón de pared delgada

El cajón rectangular (w=h=1 m, t variable, L=10 m, E=70 GPa, nu=0.3) del
benchmark `tests/test_box_torsion_bending_benchmark.py` bajo carga
distribuida (flap):

| t [m] | q [N/m] | distorsión mid-span [mm] |
|---|---|---|
| 0.010 | 200000 | 15.9 |
| 0.005 | 100000 | 62.9 |
| 0.002 | 40000 | 387 |
| 0.001 | 20000 | 1530 |

**La distorsión escala ~1/t²** — consistente con la ovalización de Brazier de
secciones cerradas de pared delgada (la rigidez anti-ovalización escala con
la flexión local de los paneles, ~t³, bajo un momento flector que escala
~t). **No es un artefacto numérico del MITC4.**

- `drilling_scale` 1→1000×: la distorsión no cambia (1529.6 → 1529.4 mm en
  t=0.001) — el drilling DOF NO es el mecanismo.
- En el cajón simétrico la ovalización no produce rotación LE→TE aparente
  (0.0000°): la distorsión es simétrica. En la pala (perfil asimétrico) la
  ovalización asimétrica se proyecta como rotación aparente del vector de
  cuerda.

**Conclusión**: el twist del shell de la pala (+8.9° tip) refleja la
ovalización REAL de las secciones de pared delgada (t/c ~ 0.001-0.003
outboard — exactamente el régimen donde el cajón muestra distorsión fuerte).
La viga de sección rígida (ancla BeamDyn/ElastoDyn) NO puede representar
esa deformación — **subestima el twist aerodinámico efectivo**. El shell
puede estar más cerca de la física real; sin un juez externo (solid
resuelto o BECAS/VABS sobre el layup) la banda honesta del de-loading es
[shell −32%, viga −15%].

## 6. Consecuencia para el FSI y el paper

- La señal de torque del rotor (objetivo del usuario) depende del de-loading
  aeroelástico: con −31.7% en vez de −15%, la deformación modula el torque
  ~2× más que el ancla — el análisis de la señal (armónicos, rainflow,
  sub-bandas) quedaría contaminado. Hay que cerrar la distorsión antes de
  fijar la campaña FSI.
- El ancla BeamDyn (−15.0%/−8.9%, consistente con Zhou 2025) es la
  referencia de cierre para el paper.

## 7. Artefactos de este diagnóstico

- `tools/run_s8c_twist_sign.py` — twist elástico bajo cargas AD (malla 0.5).
- `docs/validation_data/generated/s8c_twist_sign.csv` — perfil spanwise.
- Test de regresión: `tests/test_iea15mw_s8c_twist_sign.py`.

## 8. Cierre de la vía "artefacto de extracción / malla" (2026-09-18)

Se investigaron sistemáticamente las dos hipótesis restantes para el twist
excesivo (extracción y malla). **Ambas quedaron descartadas**:

**Extracción — el twist es robusto a la medida.** Todas las construcciones
dan el mismo resultado outboard (~8-10°):
- SVD del contorno global-z (feedback BEM actual): +8.8° tip.
- SVD proyectado sobre la tangente local deformada: +9.7°.
- Rotación del vector LE→TE (frame Camarena/Almeida) sobre la tangente
  local: +8.8° tip (con el fallback corregido, ver abajo).
- `section_twists_deg` (best-fit de rotación de sección del S-8c oficial):
  plateau 7.7-9.0° en z=104-113, 5.0° en la última estación (borde libre).

**Bug corregido en el participante**: el fallback de la tangente para los
strips del tip usaba una banda *por debajo* de la sección invertida, lo que
producía tangente invertida y twist de signo opuesto en los últimos dos
strips (−8.7°, −10.7°). Con el fallback corregido (la ventana ±2 m se
desliza hacia abajo manteniendo hi sobre lo) todos los strips del tip dan
+8..+10.7°. El promedio tip(−6) de +2.37° que aparecía en los reportes era
una cancelación de signos de ese bug, no una medida física.

**Malla — el twist no depende de la malla.** `n_samples` 120/240/300/480
(9251-9262 nodos, misma topología efectiva) dan el mismo perfil; el
espaciado del contorno `cosine` (más puntos en LE/TE) *aumenta* el twist
(+10.5° tip). El refinamiento del contorno no lo reduce.

**Conclusión actualizada**: el twist de ~8-10° del shell bajo la carga
one-way S-8c es una propiedad del modelo shell (sección de pared delgada
del layup oficial) bajo flexión flap — la distorsión Brazier de los paneles
del TE. La pregunta abierta es por qué nuestro shell distorsiona ~8× más
que la viga cuando los papers (Almeida/Camarena) reportan diferencias
shell-viga de solo 14-24% en twist. Próximos candidatos: mapeo del layup
del TE al mesh (refuerzo), o reproducir el benchmark del paper con su
geometría/layup para aislar el pipeline.

**Números de referencia consolidados** (cargas S-8c AD rated, one-way):
ancla BeamDyn +0.55 mid / +0.98 tip; S-8c oficial +1.81 mid / +5.01 tip
(plateau outboard 7.7-9.0); participante BEM +2.2 mid / +8.8 tip; de-loading
FSI smoke −31.7% (feedback leyendo el twist del SVD) vs ancla −15.0%.

## 9. Vías del layup, de las cargas y del frame (2026-09-18, continuación)

**El layup del TE está bien mapeado.** Inspección directa del modelo NuMAD
(672 element sets): la familia `HP/LP_TE_REINF` reproduce el yaml oficial —
`glass_uni` 30 mm en z=13-85 m, decreciendo 25→5 mm hasta z≈87, y ausente
outboard (z>89), donde queda solo la piel `glass_triax` 1+1 mm + gelcoat,
exactamente como el diseño oficial. No hay un refuerzo perdido.

**La rotación de cargas no es la causa.** El S-8c transforma fn/ft al
sistema global con el twist del archivo AD (`fx=c·ft−s·fn`). Repetido el
caso con 3 transformaciones (actual, sin rotación, signo opuesto): el twist
del tip cambia solo ±0.4° (8.83 / 8.59 / 8.35) y la deflexión del tip queda
igual (14.7-14.8 m). La transformación no amplifica el twist.

**La rotación del frame es real (4-path, Camarena).** Medida la rotación
de los dos vectores de la sección — cuerda LE→TE y espesor LP→HP —
alrededor de la tangente local deformada, con posiciones deformadas reales:
ambos paths rotan 6-13° en el tip (promedio 7-11°), consistente con el
LE-TE del participante (6-10.7°). No es distorsión antisimétrica (que
cancelaría en el promedio): la sección rota como frame.

**Estado**: descartadas extracción, malla (n_samples y spacing), layup y
rotación de cargas. La única vía abierta es la física del modelo shell —
la distorsión Brazier de la sección de pared delgada oficial bajo flexión.
En cola: convergencia con `element_size` 0.5/0.25/0.125 (job 11596582) —
el test preliminar con 0.125 m daba tip ≈6.4° vs 8.8° con 0.5 m, señal de
que refinar reduce el twist ~20%. Cierre pendiente: convergencia + un juez
externo (solid resuelto / BECAS) sobre el layup oficial.

## 10. Convergencia final y veredicto (2026-09-18)

**La malla está completamente descartada.** Dos barridos:

1. `element_size` 0.5 → 0.25 → 0.125 m (job 11596582): el twist **no
   converge hacia abajo**: mid estable ~2.2°, tip +8.8/+11.3/+9.7°,
   **máximo outboard +10.7/+16.1/+16.7°**. El refinamiento lo aumenta.
2. Contorno fino combinado (job 11596588, `n_samples` 600 + es 0.25/0.125):
   **idéntico a n_samples=300** (diferencias < 0.05°): tip +11.346 vs
   +11.344, max +16.113 vs +16.081. La discretización del TE no influye.

**El layup completo está verificado**: `HP_SPAR` (CarbonUD) reproduce el
yaml — 98 mm en z=40-50, decreciendo a 5.2 mm en z=111, 0 en z=115.8.
`HP_TE_REINF` (glass_uni) 30 mm en z=13-85, muere en z≈89. Ambos mapeados
correctamente desde el yaml oficial.

**Veredicto**: bajo las cargas AD rated one-way (deflexión del tip
14.3-14.7 m), el shell MITC4 + CLT con la sección oficial de la IEA-15
produce un twist de frame outboard de **~9-17°** (creciente con el
refinamiento), robusto a extracción, malla, contorno, layup y rotación de
cargas. La viga BeamDyn con el mismo layup da **+0.98°**. La diferencia es
la distorsión de sección que el modelo de placa resuelve y la viga no.
La literatura (Almeida/Camarena, DTU 10MW) reporta diferencias shell-viga
de solo 14-24% — nuestro modelo excede eso en ~un orden de magnitud, lo que
deja dos explicaciones: (a) el CLT/MITC4 exagera la ovalización de esta
sección, o (b) la pala IEA-15 outboard (t/c≈0.001-0.002 efectivo) es un
caso más extremo que el DTU 10MW. **Cierre pendiente: juez externo** —
solid 3D layer-wise de un segmento outboard (z≈95-117) con el layup
oficial, o BECAS/VABS — comparando la rotación de la cuerda con el shell.

**Implicancia FSI**: el de-loading −31.7% del smoke es la consecuencia
directa de este twist leído por el feedback BEM. Sin el juez externo no se
puede afirmar si el de-loading correcto es ≈−15% (viga/ancla) o mayor.

## 11. La lectura de Camarena y el candidato final: el solver lineal

Leído Camarena & Anderson 2025 (Sandia, *Composite Structures* 360, "A
critical verification of beam and shell models of wind turbine blades",
pala de 100 m, `tests/IEA15MW/validation papers/shell_vs_beam/
1-s2.0-S0263822325001643-main.pdf`). Puntos clave:

1. **Su setup es "static, large-deflection, flap-wise load"** — análisis
   geométricamente NO LINEAL. Nuestro S-8c y todos los tests del twist
   corrieron con `StaticLinearSolver` (teoría lineal) con 14.7 m de
   deflexión del tip y rotaciones de sección de 8-17° — **la teoría lineal
   no es válida en ese régimen**.
2. Su método de twist es el frame de 4 paths (LE/TE/HP/LP) con la base
   proyectada sobre la tangente — idéntico a lo implementado aquí.
3. Sus resultados: beam ≈ solid (twist +5.6%), shell −24% (1.5° de
   déficit) vs solid. **Su shell, viga y solid concuerdan dentro de ±25%**
   — un orden de magnitud menos que nuestra discrepancia shell-viga.
4. Su artefacto de punta del shell ("twist slightly reduces towards the
   tip") se corrigió fijando el número de elementos en la dirección
   circumferencial (hoop) — en nuestro caso n_samples 300→600 no cambió
   nada, así que nuestro exceso no es ese artefacto.
5. Advertencia literal: "Shell models should be used with caution in
   aeroelastic studies, where the torsional response could have a
   significant impact."

**Hipótesis final**: el twist de 8-17° es un artefacto de la cinemática
lineal (rotaciones grandes sin stress-stiffening de membrana, que en la
realidad limita la ovalización). Test en cola: S-8c con
`StaticNonlinearSolver` (SNES Newton-Raphson con continuación de carga,
job 11596594). Si el twist colapsa a ~1-3°, el cierre es: rehacer el S-8c
no lineal y reevaluar el feedback del FSI.

## 12. Cierre: solver no lineal terminado + resultado (2026-09-18)

**El solver no lineal estaba incompleto y se terminó.** Verificación del
código: el camino SNES del Rust usa `F_int` Green-Lagrange completo
(`assemble_fint(u, true)`) con la tangente `assemble_kt(u)`, pero la
**tangente es inconsistente con el F_int a grandes desplazamientos**
(test de diferencias finitas direccionales: ratio 6.6e-4 con deformación
chica → 8.0e-2 con la deformación real del S-8c). El job 11596594
divergió por line search exactamente por eso. El diseño original
(documentado en el docstring del Rust) era el **Updated-Lagrangian
incremental** con F_int linealizado y `update_reference` tras cada paso —
las piezas ya existían en el Rust. Se completó:

- `StaticNonlinearSolver.solve()` ahora ejecuta el loop UL: en cada paso
  resuelve `K_e(x_ref)·du = dλ·F_ext` con el ensamble lineal y avanza la
  referencia con `update_reference(du)` (32 pasos recomendado).
- Fallback GMRES+LU agregado en `linear_static_solve` (Rust, recompilado):
  el CG falla cuando K_e se vuelve indefinida por la compresión de
  membrana en estados muy curvados.
- Validación contra la elástica exacta (Bisshopp & Drucker): α=1
  (−0.5%), α=2 (+5.8%) — tests nuevos en
  `tests/test_static_nonlinear_ul_elastica.py`; los 2 tests de
  "divergencia" antiguos convertidos en tests de convergencia física.
  4/4 pasan; 23/23 con los vecinos.

**Resultado del S-8c no lineal**: twist mid +1.94, tip +9.69, max +12.63°
(deflexión del tip 14.7 m, igual que el lineal). **La linealidad NO era la
causa del exceso** — el shell no lineal tuerce igual o más que el lineal.
También descartado el momento torsor de brazo de las cargas nodales del
contorno (2.6 kN·m integrados, despreciable vs los ~200 kN·m del Mp).

**Cierre del diagnóstico**: el twist de ~9-13° del shell outboard bajo
las cargas AD one-way es la física del modelo (distorsión de sección que
la viga no representa) — robusto a extracción, malla, contorno, layup,
cargas, linealidad y momento de brazo. La comparación cuantitativa con
Camarena no procede (él usa una tracción uniforme arbitraria, pala de
100 m del IEA Task 37). La discrepancia del paper — de-loading FSI −31.7%
vs ancla −15.0% — sigue abierta y solo la resuelve un juez externo:
solid 3D layer-wise del segmento outboard, o datos experimentales.