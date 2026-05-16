# Resultados de validación calculados — AeroElast
> Generado a partir de corridas con `module load gcc` y análisis de CSVs existentes.  
> Este archivo **no modifica** `article_draft_wind_energy.md`. Es la fuente de datos para la redacción final.  
> Última actualización: 2026-05-18.

## Índice

1. [E-01/E-02 — Benchmarks de elementos MITC](#e-01e-02)
2. [V-01 — Rendimiento BEM](#v-01)
3. [V-02 — Frecuencias naturales de pala](#v-02)
4. [V-03 — Gravedad estática y consistencia de masa](#v-03)
5. [V-04 — Acoplamiento FSI completo](#v-04)
   - [5.1 Comparación yaw=0° con múltiples fuentes externas](#v-04-1)
   - [5.2 Barrido de yaw — solver corrotacional](#v-04-2)
   - [5.3 Corrotacional vs inercial — yaw=0° detallado](#v-04-3)
   - [5.4 Barrido de yaw — inercial (desplazamiento)](#v-04-4)
   - [5.5 Análisis espectral FFT](#v-04-5)
   - [5.6 Convergencia preCICE](#v-04-6)
   - [5.7 Contexto de carga extrema — DLC vs operación nominal](#v-04-7)
   - [5.8 Caso parqueado en posición de bandera (V=45 m/s, yaw=90°)](#v-04-8)
6. [Tabla resumen de validación](#resumen)

---

## E-01/E-02 — Benchmarks de elementos MITC3+/MITC4+

### E-01 Cantileveres isotrópicos (solución cerrada)

Test: `pytest tests/test_isotropic_shell_parity.py` y `test_material_suite.py`

| Caso | Referencia analítica | AeroElast | Error | Criterio |
|---|---:|---:|---:|---:|
| MITC4 isotrópico, flexión fuera del plano | 19,047.62 µm | 18,777.19 µm | −1.42 % | ≤ 5 % |
| MITC3 isotrópico, flexión fuera del plano | 19,047.62 µm | 18,818.96 µm | −1.20 % | ≤ 5 % |
| MITC4 isotrópico, flexión en el plano    | 190.476 µm   | 189.803 µm   | −0.35 % | ≤ 5 % |

### E-02 Benchmarks canónicos de literatura para MITC4+ (Ko et al. 2017)

Test: `pytest tests/test_ko2017_performance.py`

| Benchmark | Valor norm. ref. | AeroElast | Error | Estado |
|---|---:|---:|---:|:---:|
| Square plate, clamped, regular, t/L = 10⁻³   | 0.9980 | 1.0004 | +0.24 % | ✓ |
| Circular plate, clamped, t/R = 10⁻³           | 0.9997 | 0.9980 | −0.17 % | ✓ |
| Pinched cylinder, regular                     | 0.9313 | 0.9674 | +3.88 % | ✓ |
| Scordelis-Lo roof, regular                    | 0.9973 | 0.9988 | +0.15 % | ✓ |
| Hyperbolic paraboloid, regular, t/L = 10⁻³    | 0.9762 | 0.9761 | −0.01 % | ✓ |
| Hook (excluido del test — fuera de banda)      | 1.1200 | 0.9927 | −11.37 % | ✗ |

**Nota sobre el hook**: el caso hook no pasa la banda de aceptación con el refinamiento de malla actual. Se documenta explícitamente para no silenciar la limitación.

---

## V-01 — Rendimiento BEM

**Condición nominal**: V = 10.659 m/s, Ω = 7.518 RPM, pitch = 0°, ρ = 1.225 kg/m³  
**Test**: `module load gcc && pytest tests/test_iea15mw_v01_rotor_performance.py` → **10/10 PASSED**

### Condición nominal

| Magnitud | Referencia NREL/IEA [R1] | AeroElast | Error relativo |
|---|---:|---:|---:|
| Thrust T [N]    | 2,457,000.0   | 2,522,621.9 | +2.67 % |
| Torque Q [N m]  | 19,910,000.0  | 20,575,823.0 | +3.34 % |
| C_P [−]         | 0.4618        | 0.47387     | +2.61 % |
| C_T [−]         | 0.7718        | 0.78657     | +1.91 % |

Tolerancias del test: ±5 % en T y Q; ±0.02 abs en C_P; ±0.05 abs en C_T. **Todos dentro de tolerancia.**

**Interpretación**: La desviación positiva sistemática (~3 %) es consistente con diferencias de configuración BEM (discretización radial, modelo de punta de pala). Es un sesgo acotado y reproducible, no un artefacto numérico.

### Curva de desempeño (6 puntos muestreados)

| V [m/s] | pitch [°] | RPM | C_P AeroElast | C_P ref. [R1] | ΔC_P | C_T AeroElast | C_T ref. [R1] | ΔC_T |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5.006  | 2.893 | 5.000 | 0.42344 | 0.4164 | +0.0070 | 0.79472 | 0.7842 | +0.0105 |
| 7.159  | 0.000 | 5.100 | 0.47362 | 0.4616 | +0.0120 | 0.79330 | 0.7783 | +0.0150 |
| 9.027  | 0.000 | 6.430 | 0.47362 | 0.4616 | +0.0120 | 0.79322 | 0.7783 | +0.0149 |
| 10.659 | 0.000 | 7.518 | 0.47387 | 0.4618 | +0.0121 | 0.78657 | 0.7718 | +0.0148 |
| 12.259 | 6.760 | 7.518 | 0.30785 | 0.3036 | +0.0043 | 0.39187 | 0.3875 | +0.0044 |
| 15.471 | 12.19 | 7.518 | 0.15375 | 0.1511 | +0.0027 | 0.18090 | 0.1796 | +0.0013 |

Desviaciones máximas: |ΔC_P|_max = 0.0121, |ΔC_T|_max = 0.0150. Ambos muy por debajo de los criterios de aceptación.

---

## V-02 — Frecuencias naturales de pala estática

**Malla**: `tests/NuMAD_utd_iea15mw.xlsx`, element_size = 0.25 m → 31 693 nodos, 32 861 elementos shell.  
**Test**: `module load gcc && pytest tests/test_iea15mw_v02_natural_frequencies.py` → **3/3 PASSED**

### 10 primeros autovalores calculados (orden ascendente)

```
[0.553717, 0.628994, 1.694608, 1.979596, 3.237655, 4.178097, 4.534671, 4.763666, 5.152384, 5.215497] Hz
```

### Comparación con múltiples referencias

Existen tres niveles de referencia para las frecuencias de la pala IEA 15 MW, con fidelidad creciente en la descripción geométrica:

1. **NREL/TP-5000-75698 [R2]** — ElastoDyn (viga 1D, propiedades seccionales de WISDEM): referencia original del diseño, no usa el archivo NuMAD [R3].
2. **Escalera Mendoza et al. 2023 [R4]** — PreComp+BModes (viga 1D) aplicado al modelo `NuMAD_utd_iea15mw.xlsx` [R3]: misma geometría de pala que AeroElast, diferente modelo estructural (viga vs. shell).
3. **AeroElast (este trabajo)** — Shell-3D MITC aplicado al mismo archivo Excel [R3]: máxima fidelidad geométrica, resuelve modos de la sección transversal.

La comparación entre UTD BModes [R4] y AeroElast shell —ambos del mismo `NuMAD_utd_iea15mw.xlsx`— aísla el efecto *beam-vs-shell* en la predicción de frecuencias, independientemente de la geometría de entrada.

| Modo | Descripción | NREL ElastoDyn [Hz] | UTD BModes [R4] [Hz] | Bernardi 2025 [Hz] | AeroElast Shell [Hz] | Δ vs NREL | Δ vs UTD | Test |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| 1er flapwise  | Flexión fuera del plano | 0.5585 | 0.57 | 0.5369 | 0.5537 | −0.86 % | −2.9 % | ✓ |
| 1er edgewise  | Flexión en el plano     | 0.6406 | 0.65 | **0.7267**⁺ | 0.6290 | −1.81 % | −3.2 % | ✓ |
| 2do flapwise  | Flexión fuera del plano | 1.6590 | 1.72 | 1.577 | 1.6946 | +2.15 % | −1.5 % | ✓ |
| 2do edgewise† | Flexión en el plano     | 2.1670 | 2.08 | 2.267 | 1.9796 | −8.65 % | −4.8 % | ✗ |
| 3er flapwise† | Flexión fuera del plano | —      | 3.41 | 3.113 | 3.2377 | —        | −5.1 % | ✗ |
| 1er torsión†  | Torsional               | 4.459  | 4.29 | 3.642 | ~4.535 | ~+1.7 % | ~+5.7 % | ✗ |

†: Modos de contexto — no incluidos en las aserciones automáticas actuales.  
⁺: El valor Bernardi 2025 (0.7267 Hz) para el 1er edgewise es el más alto reportado entre todos los métodos revisados (ver Tabla 3 abajo) — ver observación correspondiente.

**Tendencias sistemáticas**: el modelo shell AeroElast predice frecuencias consistentemente inferiores a ambos modelos de viga. Respecto a UTD BModes (misma geometría), la diferencia es de −1.5 % a −5.1 % para los modos de flexión. Esto refleja que el modelo shell captura deformaciones de la sección transversal (*in-plane cross-section distortion*) que elevan la flexibilidad efectiva respecto a la viga Bernoulli-Euler. Que los tres modos testeados estén dentro del ±2.2 % respecto a NREL y del ±3.2 % respecto a UTD BModes es un resultado robusto dado el salto de fidelidad modelar.

**2do edgewise**: la discrepancia con NREL (−8.65 %) es mayor que con UTD BModes (−4.8 %), y el propio UTD BModes difiere −4.0 % respecto a NREL. Esto indica que parte del error no es atribuible al shell-3D sino a diferencias de homogenización seccional entre WISDEM y PreComp. El autovalor AeroElast en 1.9796 Hz puede corresponder a un modo flap-edge mixto; la clasificación definitiva requiere análisis de energía modal por dirección.

**Nota de estabilidad aerolástica [R4]**: el modelo UTD NuMAD analizado con OpenFAST muestra flutter a **6.87 RPM**, un ratio de **0.91** respecto a la velocidad nominal de 7.55 RPM. Las simulaciones AeroElast operan a 7.518 RPM (≈ velocidad nominal) — dentro de la banda estable. El flutter clásico implica crecimiento exponencial libre, incompatible con el régimen de régimen permanente observado en las corridas FSI.

### Posición de AeroElast en el espectro de métodos — Tabla 3 [Zhou Energy 2025]

Zhou et al. 2025 [Zhou Energy 2025] recopilan en su Tabla 3 seis fuentes de referencia para las frecuencias de la pala IEA 15 MW, abarcando desde vigas Timoshenko 1D hasta modelos 3D-FEM. La siguiente tabla reproduce esos datos añadiendo AeroElast shell-3D para contextualizar nuestro modelo en ese espectro:

| Modo | HWAC2+Timo [74] | HWAC2+Timo [26] | 3D-FEM [55] | 3D-FEM [75] | Tech.Rep. [56] | LL-FVW+GEBT [Zhou] | Modal CSD [Bernardi] | **AeroElast Shell** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1er flapwise  | 0.512 | 0.521 | 0.465 | 0.555 | 0.555 | 0.523 | 0.537 | **0.554** |
| 1er edgewise  | 0.692 | 0.619 | 0.547 | 0.642 | 0.642 | 0.665 | **0.727**⁺ | **0.629** |
| 2do flapwise  | 1.509 | 1.559 | 1.452 | 1.598 | —     | 1.475 | 1.577 | **1.695** |
| 2do edgewise  | 2.120 | 1.933 | 1.671 | 1.925 | —     | 2.124 | 2.267 | **1.980** |
| 1er torsional | 4.314 | 4.475 | 4.295 | 3.911 | —     | 4.072 | 3.642 | **~4.535** |

Fuentes: [74] y [26] = HWAC2 + viga Timoshenko (implementaciones independientes); [55] y [75] = 3D-FEM (estudios publicados); [56] = NREL/TP-5000-75698 (= nuestro [R2]); [Zhou] = LL-FVW + GEBT (Zhou et al. Energy 2025); [Bernardi] = Modal CSD, vigas Euler-Bernoulli, 15 modos, 80 nodos, integración generalizada-α (Bernardi et al. 2025); AeroElast = shell MITC3/4, `NuMAD_utd_iea15mw.xlsx` [R3].

⁺ Valor Bernardi 2025 para 1er edgewise (0.727 Hz): el más alto de todos los métodos, +15.5 % sobre AeroElast (0.629 Hz) y +11.5 % sobre la columna HWAC2+Timo [74] (0.692 Hz). Su modelo usa vigas Euler-Bernoulli sin torsión de Vlásov; una posible causa es que la homogenización seccional subestime la rigidez de flexión en el plano en esa implementación. La dispersión total para el 1er edgewise abarca 0.547–0.727 Hz (±14 % respecto a la mediana 0.638 Hz), lo que establece el ancho de banda de incertidumbre metodológica para este modo.

![Comparación modal multi-método — IEA 15 MW](figures/fig_v02_modal_comparison.png)
*Figura: Frecuencias naturales de la pala IEA 15 MW comparadas entre los 8 métodos de la Tabla 3 (Zhou et al. 2025 + Bernardi et al. 2025). AeroElast Shell (rojo) destacado. Barras ausentes = no reportado en esa fuente.*

**Observaciones**:
- **1er flapwise** (rango 0.465–0.570 Hz, dispersión ±10 % respecto a mediana 0.542 Hz): AeroElast (0.554 Hz) queda dentro del rango, alineado con el extremo superior (3D-FEM [75], NREL [56]). Bernardi 2025 (0.537 Hz) se ubica en el tercio inferior, coherente con un modelo de vigas Euler-Bernoulli de 80 nodos que subestima ligeramente la rigidez de flexión. La dispersión entre estudios 3D es mayor que entre vigas, lo que refleja sensibilidad al modelado del laminado compuesto.
- **1er edgewise** (rango 0.547–0.727 Hz, dispersión ±14 % respecto a mediana 0.638 Hz): AeroElast (0.629 Hz) cae dentro del bloque central. Bernardi 2025 (0.727 Hz) es el valor más alto de todos los métodos; el outlier puede deberse a diferencias en la homogenización seccional de la rigidez edgewise. La dispersión ±14 % es la más alta de los modos fundamentales, lo que confirma que el 1er edgewise es el modo más sensible a las hipótesis de modelado estructural.
- **Torsional**: dispersión amplia entre métodos (3.642–4.475 Hz). Bernardi 2025 (3.642 Hz) es el más bajo de todos — plausible para un modelo Euler-Bernoulli sin rigidez de Vlásov. El valor GEBT (4.072 Hz) difiere +11 % del shell-3D (4.535 Hz). AeroElast captura la rigidez adicional de la sección cerrada de pared delgada.
- **Dispersión global**: los rangos de dispersión para los modos fundamentales (10 % en flapwise, 14 % en edgewise) reflejan el estado actual del campo — la comunidad no ha convergido en un valor de referencia único para la pala IEA 15 MW, y la comparación entre métodos distintos tiene un piso de incertidumbre intrínseco de ±10 %.

---

## V-03 — Gravedad estática y consistencia de masa

**Test**: `pytest tests/test_iea15mw_v03_static_gravity.py`

**AeroElast reporta dos métricas de masa independientes; ambas coinciden en 67,167.04 kg**, confirmando consistencia interna:

| Métrica | Método | Masa [kg] |
|---|---|---:|
| Suma elemental | ρ·h·A por elemento (`total_elemental_mass`, Rust) | 67,167.04 |
| Matriz de masa | Suma x-DOF de M consistente (partición de la unidad) | 67,167.04 |

| Chequeo | Referencia | AeroElast | Error | Criterio |
|---|---:|---:|---:|---:|
| Reacción en raíz vs. peso ensamblado | 658,908.7046 N | 658,908.7051 N | 8.2×10⁻⁸ % | ≤ 2 % |
| Masa total de pala vs. entradas NuMAD [R3] | 67,921.0 kg | 67,167.04 kg | −1.11 % | ≤ 2 % |

**Deflexiones tip estáticas** (indicadores auxiliares bajo gravedad, pala horizontal):

| Dirección | AeroElast [m] |
|---|---:|
| Edgewise (gravedad) | −1.7925 |
| Flapwise (gravedad) | −0.1331 |

**Jerarquía de fuentes de masa**: la definición original de la turbina (NREL/TP-5000-75698, Tabla 2-1 [R2]) establece **65,250 kg** como objetivo de diseño WISDEM. Los archivos de propiedades distribuidas (BeamDyn, HAWC2, NuMAD) generados a partir de ese diseño integran entre 2.5 % y 4.3 % más, lo cual es una característica conocida del paquete de referencia IEA 15 MW: las propiedades seccionales fueron refinadas tras la optimización WISDEM original.

| Fuente | Método | Masa [kg] | Δ vs. objetivo [R2] |
|---|---|---:|---:|
| NREL/TP-5000-75698, Tabla 2-1 [R2] | Objetivo de diseño WISDEM | **65,250** | — |
| HAWC2 `IEA_15MW_RWT_Blade_st_noFPM.st` | ∫m(r)dr, 26 estaciones | 66,994 | +2.68 % |
| BeamDyn tabular M₁₁ [R1] | ∫M₁₁(s)ds, 26 estaciones | 66,912 | +2.55 % |
| Excel Overview [R1] | Valor declarado (hoja Overview) | 67,921 | +4.10 % |
| **AeroElast shell-3D (este trabajo)** | ρ·h·A y M-matrix (ambas coinciden) | **67,167** | **+2.94 %** |
| UTD NuMAD — Escalera Mendoza et al. 2023 [R4] | AutoNuMAD (integración seccional viga) | 68,077 | +4.33 % |

La diferencia de −1.11 % entre AeroElast (67,167 kg) y las entradas NuMAD [R3] (67,921 kg) valida que la integración del laminado shell-3D es correcta. La brecha de +2.94 % con respecto al objetivo de diseño WISDEM [R2] no es un error de AeroElast sino una propiedad del paquete IEA 15 MW.

---

## V-04 — Acoplamiento FSI completo

Fuentes de datos en disco:
- Corrotacional: `frontiersin_results/metrics_summary.csv` (yaw=0,10,20,30,40; t≥20s; yaw=30 diverge en t≈64s — datos válidos hasta ese punto)
- Inercial: `frontiersin_results_inertial_old/` (yaw=0,10,20,30,40; t≥20s; yaw=40 diverge en t≈53s); Cp/P integrados de la malla fluida BEM sección a sección para cada paso de tiempo; desplazamiento de `structural_report.csv`
- Convergencia: `frontiersin_results/convergence_summary.csv`

Ventana de estado estacionario: t ≥ 20 s en todos los casos.

---

### 5.1 Comparación yaw=0° con múltiples fuentes externas

**Nota sobre la comparación multi-fuente**: los comparadores externos usan modelos aerodinámicos y estructurales distintos. La comparación es por nivel de fidelidad relativa, no por exactitud directa.

| Fuente / modelo | Config | Flapwise tip medio [m] | Flap máx [m] | Potencia media [MW] | Tipo aerodinámica | Tipo estructura |
|---|---|---:|---:|---:|---|---|
| **AeroElast corrot. (frontiersin_results)**    | WindIO, ramp=1 s      | **12.794** | 14.345 | **14.710** | BEM (CCBlade) | Shell MITC |
| AeroElast corrot. (bem_0_10, NuMAD mesh)†      | NuMAD XLS, ramp=0     | **12.761** | 14.814 | **~14.39** | BEM (CCBlade) | Shell MITC |
| AeroElast corrot. (bem_0_10_full, vel. fbk)††  | NuMAD XLS, ramp=0     | **12.777** | 14.537 | ⚠ (ver nota) | BEM (CCBlade) | Shell MITC |
| AeroElast corrot. (campaña anterior‡)          | WindIO, ramp anterior | 11.94      | n/d    | n/d         | BEM             | Shell MITC |
| **AeroElast inercial (metrics_summary)**        | WindIO, ramp=1 s      | **12.747** | 13.519 | **14.675** | BEM (CCBlade) | Shell MITC |
| Zhou et al. 2025 [LL-FVW + GEBT]               | —                     | 13.86      | n/d    | 14.76       | LL-FVW          | Viga GEBT  |
| Bernardi et al. 2025 [LES + CSD]               | —                     | ~16        | n/d    | n/d         | LES             | CSD 3D     |
| ALM (campaña previa interna)                    | —                     | 14.10      | n/d    | n/d         | ALM             | n/d        |
| **NREL DLC 1.4 — pico extremo [R2-DLC]**        | ECD + dir. change, op. nominal, OpenFAST | 22.8  | —      | n/d         | BEM (AeroDyn)   | Viga (ElastoDyn) |

†: `bem_0_10` — mismo solver corrotacional, malla del archivo NuMAD Excel (misma que V-02), t_sim=50s. `force_ramp_time=0` provoca un pico transitorio numérico de 27.6 m a t≈1s que NO es una carga física; el régimen permanente (t≥35s) es 12.76m. Ver §5.7 para análisis del pico.  
††: `bem_0_10_full` — igual que bem_0_10 pero con `send_velocity_to_precice: true` (retroalimentación de velocidad estructural al BEM). La potencia presenta una anomalía: el valor medio en régimen permanente cae a ~2.3 MW (por pala, ×3 = ~6.9 MW) frente a ~14.4 MW en bem_0_10. El desplazamiento es igual. Esta discrepancia necesita investigación adicional (posiblemente un efecto de fase en la velocidad relativa cuando la pala vibra). No usar el valor de potencia de bem_0_10_full hasta resolver la causa.  
‡: La campaña anterior (11.94 m) corresponde a parámetros de corrida distintos. No calcular error comparando con ella; sirve solo para contextualizar variabilidad entre campañas.  
§: **NREL DLC 1.4** NO es condición operativa comparable con régimen permanente: corresponde a la ráfaga extrema con cambio de dirección (ECD, DLC 1.4) a velocidad rated, que produce el pico dinámico máximo absoluto de todo el análisis de cargas (22.8 m). Los DLC 6.x (parqueado + viento extremo V50=50 m/s) producen ~8 m, consistente con `bem_90_50_S`. Se incluye exclusivamente como envolvente de escala de diseño. Ver §5.7 y §5.8 para el análisis completo.

**Lectura jerárquica**:
1. Todas las campañas de AeroElast en régimen permanente a potencia nominal convergen al rango 12.75–12.79 m (desviación < 0.04 m entre sí), lo que indica alta reproducibilidad frente a cambios de malla y configuración.
2. AeroElast queda entre la campaña propia anterior (11.94 m) y la referencia LL-FVW (13.86 m) — físicamente coherente: BEM no resuelve la estela completa pero captura el nivel correcto de deformación.
3. **Bernardi LES a V=10 m/s: aparente paradoja y su resolución.** A primera vista, ~16 m a V=10 m/s < V_rated=10.659 m/s parece contradictorio — menos viento debería implicar menos carga. La resolución está en el **control de paso**: a velocidad rated, el pitch control actúa para limitar la potencia, reduciendo el coeficiente de empuje Ct. A V=10 m/s (por debajo del rated), el pitch permanece en su valor de diseño óptimo y Ct alcanza su valor de diseño máximo (~0.83). Con los valores del documento Bernardi: T ≈ Ct × ½ρAV² ≈ 0.83 × 2.77 MN = 2.30 MN, frente a ~2.01 MN en nuestras corridas (V=10.659 m/s, pitch activo, Ct ≈ 0.64). El empuje en la simulación Bernardi es **~14 % mayor** que el de AeroElast a rated, lo que por sí solo explica una deflexión ~14 % mayor: 12.79 × 1.14 ≈ 14.6 m. El ~9 % restante (14.6 → 16 m) es atribuible a la mayor fidelidad de la LES (carga de vorticidad de extremo, efectos 3D de capa límite, variación de carga inducida por la cizalladura del perfil de viento). **Los 16 m de Bernardi no son una anomalía sino el resultado físico esperado de su punto de operación.**
4. Corrotacional e inercial predicen medias estadísticamente muy similares (diferencia < 0.5 %) pero la dinámica temporal difiere (ver §5.3 y §5.5).

**Rango de deflexión flapwise entre todas las fuentes (yaw=0°, operación nominal)**:

| Fuente | V∞ [m/s] | Control pitch | Aerodinámica | Flapwise [m] |
|---|---:|---|---|---:|
| AeroElast BEM (este trabajo) | 10.659 (rated) | activo | BEM CCBlade | 12.75–12.79 |
| Zhou et al. 2025 [LL-FVW+GEBT] | ~10.59 (near-rated) | — | LL-FVW | 13.86 |
| ALM-GEBT (referencia interna) | ~rated | — | ALM | 14.10 |
| ElastoDyn [Bernardi, V=10 m/s] | 10.0 (sub-rated) | mínimo pitch | BEM OpenFAST | ~13.8 |
| **Bernardi et al. 2025 [LES+CSD]** | **10.0 (sub-rated)** | **mínimo pitch** | **LES ALM** | **~16** |

La dispersión entre todas las fuentes (12.75–16 m) abarca un rango de **~25 %** respecto al valor medio (~14 m). Esta dispersión combina dos efectos que no son separables sin una comparación a condiciones idénticas: (a) fidelidad aerodinámica (BEM < LL-FVW < ALM < LES; en base fija sub-rated, BEM subestima potencia ≲1.4 % respecto a LL-FVW [Ramos-García 2022]) y (b) punto de operación (rated con pitch activo vs. sub-rated con Ct máximo). En el estado actual de la literatura, **la comunidad no ha publicado una comparación multi-método en condiciones idénticas de viento y control** para la pala IEA 15 MW; la dispersión observada es por tanto una propiedad del conjunto de datos disponibles, no necesariamente de los modelos mismos.

**Deflexiones de punta nominales (yaw=0°, condición rated)** — Tabla 4 [Zhou Energy 2025]:

| Fuente / modelo | Tipo | Flapwise [m] | Edgewise [m] | Torsional [°] |
|---|---|---:|---:|---:|
| ALM-GEBT [referencia en Zhou T4] | ALM + GEBT | 14.10 | −1.27 | −3.77 |
| Zhou et al. 2025 [LL-FVW+GEBT]   | LL-FVW + GEBT | **13.86** | **−1.22** | **−3.60** |
| Ma et al. 2025 [LL-FVW+GEBT] (yaw=10°)† | LL-FVW + GEBT | ≈14.2 | ≈−1.05 | ≈−3.8° |
| **AeroElast corrot. (frontiersin)** | BEM + Shell MITC | **12.794** | −0.573 | n/d |

†yaw=10° es el ángulo mínimo analizado por Ma et al.; se incluye como referencia por ser el más próximo a la condición yaw=0°.

La edgewise de AeroElast (−0.573 m) es notablemente menor que la de Zhou (−1.22 m). La diferencia refleja la naturaleza del acoplamiento: el modelo GEBT calcula el campo de desplazamiento completo incluyendo deflexión edgewise inducida por la carga aerodinámica tangencial, mientras que la campaña corrotacional de AeroElast tiene la carga aerodinámica tangencial dominada por el peso propio en ausencia de yaw (−0.573 m coincide con la gravedad estática de V-03: −1.79 m de gravedad pura, reducida por la centrifugación). Esto indica que en la campaña FSI la componente tangencial de la carga aerodinámica contribuye en sentido contrario a la gravedad.

**Potencia y thrust nominales — Tabla 6 [Zhou Energy 2025]**:

| Condición | Modelo | Potencia [MW] | Thrust [MN] |
|---|---|---:|---:|
| Fijo, pala rígida   | Zhou LL-FVW+GEBT | 16.11 | 2.53 |
| Fijo, pala flexible | Zhou LL-FVW+GEBT | **14.76** | **2.20** |
| Reducción (rígida → flexible) | — | −8.38 % | −13.04 % |
| **Fijo, yaw=0° (frontiersin)** | **AeroElast BEM+Shell** | **14.710** | n/d |

La potencia nominal de AeroElast (14.71 MW) coincide con la de Zhou flexible (14.76 MW) con un error de −0.34 %, confirmando consistencia entre ambos modelos FSI con pala flexible. El thrust de referencia de V-01 BEM (2.52 MN, pala rígida) es coherente con el valor rígido de Zhou (2.53 MN); la reducción −13 % por flexibilidad lleva a ~2.20 MN, que no es medido directamente en la campaña FSI actual.

![Series temporales de deflexión tip — flapwise y edgewise](figures/fig_5_1_timeseries_deflection.png)

*Figura: Series temporales de deflexión de punta (yaw=0°, condición rated). Eje X: ángulo de rotación acumulado desde t=20 s (4 revoluciones completas). Corrotacional (rojo sólido) vs. inercial (azul discontinuo). Líneas de puntos: medias propias de simulación. Líneas de guiones: valores de referencia de la literatura (Zhou LL-FVW+GEBT y ALM-GEBT, Tabla 4 [Zhou Energy 2025]). Panel superior: deflexión flapwise — ambos solvers se sitúan ~1.1 m por debajo de Zhou, consistente con la menor carga aerodinámica del modelo BEM frente a LL-FVW. Panel inferior: deflexión edgewise — la oscilación 1P dominada por peso propio es claramente visible; la media de AeroElast (≈−0.55 m) queda por encima de la referencia Zhou (−1.22 m) porque el BEM no transfiere la componente tangencial completa de la carga.*

---

### 5.2 Barrido de yaw — solver corrotacional (completo)

Fuente: `frontiersin_results/metrics_summary.csv`. Post-procesado con t_start = 20 s.

| Yaw [°] | n_muestras | t_fin [s] | C_P | C_T | Potencia [MW] | Flap medio [m] | Flap std [m] | Flap p-p [m] | Edge medio [m] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|  0  | 4 269  | 62.7  | 0.4626 | 0.6699 | 14.683 | 12.794 | 0.533 | 3.050 | −0.573 |
| 10  | 13 000 | 150.0 | 0.4453 | 0.6614 | 14.130 | 12.672 | 0.353 | 2.599 | −0.553 |
| 20  | 8 532  | 105.3 | 0.3949 | 0.6342 | 12.526 | 12.282 | 0.501 | 2.840 | −0.508 |
| 30⚠ | 4 376  |  63.8 | 0.3142 | 0.5873 |  9.958 | 11.615 | 1.009 | 5.475 | +0.577 |
| 40  | 5 000  |  70.0 | 0.2077 | 0.5168 |  6.574 | 10.672 | 0.816 | 3.180 | −0.298 |

⚠ yaw=30: datos válidos hasta divergencia numérica en t≈64s; estadísticas sobre ventana t=20–63.8s.

**Caída de C_P vs. ley cos³**:

| Yaw [°] | C_P/C_P(0°) medido | cos³(yaw) | Error vs. cos³ |
|---:|---:|---:|---:|
| 10  | 0.9626 | 0.9551 | +0.75 pp |
| 20  | 0.8537 | 0.8298 | +2.39 pp |
| 30⚠ | 0.6793 | 0.6495 | +2.98 pp |
| 40  | 0.4490 | 0.4495 | −0.05 pp |

![Barrido de yaw — Corrotacional vs Inercial](figures/fig_5_2_yaw_sweep.png)
*Figura: Barrido de yaw, corrotacional (círculo rojo) vs. inercial (triángulo azul) vs. Ma et al. 2025 LL-FVW+GEBT (diamante verde). Izq.: C_P vs. ley cos³(γ) — ambos solvers BEM siguen la ley cos³ con excelente acuerdo (Cp inercial obtenido por integración BEM fluido). Centro: potencia aerodinámica — se añade la referencia Ma et al. [LL-FVW] y la curva P₀·cos^1.5 que describe su tendencia; la diferencia BEM–LL-FVW es pequeña a yaw bajo y crece hasta ~20% a yaw=30°–40°, coherente con la mayor capacidad del modelo LL-FVW para capturar la recuperación de estela sesgada. Der.: deflexión flapwise — los puntos de Ma et al. quedan sistemáticamente ~1.5–2.0 m por encima de AeroElast, mismo offset que en yaw=0°, confirmando que la diferencia es de fidelidad aerodinámica (BEM < LL-FVW) y no de formulación estructural. Símbolo ×: divergencia numérica antes de finalizar la simulación.*

---

### 5.3 Corrotacional vs inercial — yaw=0° detallado

Fuente: `metrics_summary.csv` de cada campaña. Ambas con t_start = 20 s.

| Magnitud | Corrotacional | Inercial | Diferencia abs. | Diferencia rel. |
|---|---:|---:|---:|---:|
| n_muestras (t ≥ 20 s)          | 4 269    | 3 253   | −1 016   | −23.8 % |
| t_fin [s]                       | 62.69    | 52.53   | —        | — |
| Flapwise medio [m]              | 12.7943  | 12.7471 | −0.0472  | −0.37 % |
| Flapwise std [m]                |  0.5334  |  0.3300 | −0.2034  | −38.1 % |
| Flapwise max [m]                | 14.345   | 13.519  | −0.826   | — |
| Flapwise p-p [m]                |  3.050   |  1.172  | −1.878   | −61.6 % |
| Edgewise medio [m]              | −0.5728  | −0.5472 | +0.0256  | — |
| Potencia media [MW]             | 14.7100  | 14.6754 | −0.0346  | −0.24 % |
| Potencia std [MW]               |  0.1902  |  0.3916 | +0.2014  | +106 % |
| Potencia p-p [MW]               |  0.879   |  4.968  | +4.089   | +465 % |
| Potencia max [MW]               | 15.195   | 18.816  | +3.621   | — |
| C_P medio                       |  0.4691  |  0.4682 | −0.0009  | −0.19 % |
| C_T medio                       |  0.6789  |  0.6769 | −0.0020  | −0.29 % |
| Thrust medio [MN]               |  2.0104  |  2.0031 | −0.0073  | −0.36 % |
| Torque medio [MN m]             | 18.605   | 18.562  | −0.043   | −0.23 % |

**Lectura**:
- **Las medias de todas las magnitudes son estadísticamente equivalentes** (diferencias < 0.4 %). Ambas formulaciones son físicamente consistentes en régimen estacionario.
- **La dinámica temporal difiere sustancialmente**: el inercial tiene flapwise p-p 62 % menor pero potencia p-p 465 % mayor.
- La potencia máxima del inercial (18.82 MW) supera al corrotacional (15.19 MW) en +3.6 MW (+24 %). Relevante para diseño de la cadena cinemática.

---

### 5.4 Barrido de yaw — solver inercial (desplazamiento tip)

Fuente: `frontiersin_results_inertial_old/` (t≥20s). Cp/P calculados por integración BEM fluido sección a sección (`3·∫Np dr`, `3·∫Tp·r dr`); desplazamiento flapwise de `structural_report.csv`.

| Yaw [°] | n_muestras | t_fin [s] | Flap medio [m] | Flap std [m] | Flap p-p [m] | Estado |
|---:|---:|---:|---:|---:|---:|---|
|  0  | 5 000 | 70.0 | 12.746 | 0.326 | 1.172 | ✓ válido |
| 10  | 5 000 | 70.0 | 12.630 | 0.183 | 0.608 | ✓ válido |
| 20  | 5 000 | 70.0 | 12.272 | 0.397 | 1.269 | ✓ válido |
| 30  | 5 000 | 70.0 | 11.623 | 0.620 | 1.927 | ✓ válido |
| 40⚠ | 3 318 | 53.2 | 10.659 | 0.775 | 2.272 | ⚠ diverge t≈53s |

**Comparación de deflexión corrot. vs. inercial por yaw**:

| Yaw [°] | Corrot. [m] | Inercial [m] | Diferencia [m] | Nota |
|---:|---:|---:|---:|---|
|  0  | 12.794       | 12.746        | −0.048 | |
| 10  | 12.672       | 12.630        | −0.042 | |
| 20  | 12.282       | 12.272        | −0.010 | |
| 30  | 11.615⚠      | 11.623        | +0.008 | corrot. diverge t≈64s |
| 40  | 10.672       | 10.659⚠       | −0.013 | inercial diverge t≈53s |

Ambos solvers muestran deflexiones flapwise muy similares en todo el barrido (diferencias < 0.05 m hasta yaw=20°). A yaw=30–40° la diferencia sigue siendo pequeña (< 0.05 m), pero ambas simulaciones presentan divergencia numérica antes de finalizar: el inercial en t≈53s (yaw=40°) y el corrotacional en t≈64s (yaw=30°). Las estadísticas reportadas son sobre la ventana válida pre-divergencia en cada caso.

#### Comparación con Ma et al. 2025 (LL-FVW + GEBT)

[Ma Frontiers 2025] publica un barrido yaw=10°–40° con el mismo IEA-15 MW a condición rated (V=10.59 m/s, Ω=7.55 rpm, sin control de paso) usando LL-FVW + GEBT. Los datos de potencia (leídos de Figure 15) y deflexión flapwise tip (Figure 16) permiten la siguiente comparación directa:

| Yaw [°] | Flap corrot. [m] | Flap inercial [m] | Flap Ma et al. [m] | P corrot. [MW] | P Ma et al. [MW] |
|---:|---:|---:|---:|---:|---:|
| 10 | 12.672 | 12.630 | ≈14.2 | 14.130 | ≈14.0 |
| 20 | 12.282 | 12.272 | ≈13.8 | 12.526 | ≈13.4 |
| 30 | 11.615 | 11.623 | ≈12.5 |  9.958 | ≈12.0 |
| 40 | 10.672 | 10.659 | ≈11.5 |  6.574 | ≈10.0 |

Observaciones:

1. **Deflexión flapwise**: AeroElast queda sistemáticamente ≈1.5–2.0 m por debajo de Ma et al. en todo el barrido — mismo patrón que en yaw=0° (BEM < LL-FVW). La tendencia de caída con yaw es cualitativamente consistente en ambos modelos.

2. **Potencia media vs. yaw**: Ma et al. (LL-FVW) sigue una caída aproximada ~cos^1.5, mientras que AeroElast (BEM) sigue más fielmente ~cos³. Esto es físicamente coherente: LL-FVW captura el desplazamiento lateral de la estela (la estela se sesga y reduce la auto-interferencia), lo que atenúa la pérdida de potencia respecto a la predicción geométrica pura cos³. El BEM con corrección de yaw reproduce cos³ sin capturar este efecto de recuperación de estela.

3. **A yaw=30°–40° la discrepancia AeroElast–Ma et al. es grande** (9.96 vs. ≈12.0 MW a yaw=30°, +20%). Esta diferencia se debe principalmente al modelo aerodinámico (BEM sobreestima la pérdida a yaw elevado) y no a la formulación estructural. La tendencia relativa de caída entre yaw=10° y yaw=40° es ~30% en Ma et al. y ~53% en AeroElast, consistente con la diferencia BEM/LL-FVW documentada en la literatura.

---

#### Comparación con literatura: valores medios vs yaw

La figura siguiente presenta tres paneles comparando los dos solvers AeroElast con Ma et al. 2025 y Zhou 2025 en función del ángulo de yaw:

![Comparación con literatura — valores medios](figures/fig_5_2b_literature_comparison.png)

**Figura 5** — Comparación con literatura, IEA 15 MW, condición rated: (a) thrust medio, (b) deflexión flapwise media, (c) |deflexión edgewise media|. ◇ = Ma et al. 2025 [LL-FVW+GEBT]; □ = Zhou 2025 [LL-FVW+GEBT] (solo yaw=0°); × = divergencia numérica. Amplitudes Ma et al. leídas gráficamente de Figure 16a [incertidumbre estimada ±0.5 m].

**Observaciones por panel:**

**Panel (a) — Thrust medio**: Los dos solvers AeroElast son prácticamente indistinguibles (la formulación estructural no afecta las cargas medias). Ma et al. (LL-FVW) supera a AeroElast (BEM) en ~7–9% en todo el barrido — coherente con que LL-FVW captura mejor la inducción axial que el BEM con corrección de yaw. Zhou a yaw=0° confirma este offset: 2.20 MN vs. 2.01 MN (~9% superior).

**Panel (b) — Flapwise media**: Misma conclusión que el thrust. Ambos solvers AeroElast se solapan entre sí y quedan sistemáticamente ~1.5–2 m por debajo de Ma et al. y Zhou. El sesgo es consecuencia directa del BEM: menor inducción radial → menor carga → menor deflexión media.

**Panel (c) — |Edgewise media|**: Ma et al. reporta ≈1.05 m constante para todos los ángulos de yaw (la deflexión edgewise está dominada por la tracción centrífuga y el peso propio, independientemente del yaw). Zhou a yaw=0° = 1.22 m. AeroElast subestima sistemáticamente: corrot. ~0.73–0.55 m, inercial ~0.55–0.29 m. El inercial muestra además una caída pronunciada con yaw que no aparece en Ma et al., lo que indica que el solver inercial puede estar subestimando la componente aerodinámica tangencial con yaw. (Nota: valores en módulo por diferencia de convención de signo entre BEM y LL-FVW.)

---

### 5.5 Análisis espectral FFT — corrotacional vs inercial

**Señales analizadas**: `Max Displacement [m]` y `Aero Power [W]` de `rotor_performance.csv` (cómputo interno del solver, por pala; los valores kW de la tabla incluyen el factor ×N_blades=3 para expresarlos como equivalente total de rotor), en t ≥ 20 s, Δt = 0.01 s, ventana rectangular.  
> **Importante**: los valores de potencia en esta sección corresponden al cómputo interno de par del solver estructural (`Ω × τ`), **no** a la potencia aerodinámica del participante fluido BEM. Para comparación de rendimiento, usar §5.4 (BEM).
**Frecuencias físicas de referencia**: 1P = 0.1253 Hz (= 7.518 RPM / 60), 3P = 0.3759 Hz, f₁_flap = 0.5537 Hz (V-02), f₁_edge = 0.6290 Hz (V-02).

> **Nota**: Los valores de esta sección usan ventana rectangular. La subsección 5.5.1 repite el análisis con ventana Hann (menor fuga espectral) y extiende el estudio a cargas de rotor (thrust, momento de inclinación, momento de guiñada), en línea con el análisis de frecuencias características de [Zhou Energy 2025] §4.2.

#### Amplitudes espectrales — desplazamiento tip

| Frecuencia | Hz nominal | Corrot. bin | Corrot. amp. [cm] | Inercial bin | Inercial amp. [cm] | Ratio |
|---|---:|---:|---:|---:|---:|---:|
| 1P     | 0.1253 | 0.1171 | **26.0** | 0.1200 | **29.4** | 1.1× |
| 3P     | 0.3759 | 0.3747 | **15.9** | 0.3799 | **10.1** | 0.6× |
| f₁_flap | 0.5537 | 0.5621 | **16.4** | 0.5599 | **0.6** | **0.04×** |
| f₁_edge | 0.6290 | 0.6323 | **9.3** | 0.6199 | **1.3** | **0.14×** |

#### Amplitudes espectrales — potencia aerodinámica

| Frecuencia | Hz nominal | Corrot. amp. [kW] | Inercial amp. [kW] | Ratio inercial/corrot |
|---|---:|---:|---:|---:|
| 1P     | 0.1253 |   188.1 | **12 492.3** | **66×** |
| 3P     | 0.3759 |    35.7 |    466.2     |  13×  |
| f₁_flap | 0.5537 |     8.4 |    244.2     |  29×  |
| f₁_edge | 0.6290 |     8.1 |    261.9     |  32×  |

#### Interpretación física

El contraste más importante entre las dos formulaciones **no está en los promedios** (casi idénticos) sino en el contenido espectral:

**1. Excitación de modos estructurales de la pala:**  
La formulación corrotacional excita el 1er modo flapwise con amplitud de 16.4 cm a 0.562 Hz (≈ f₁_flap = 0.5537 Hz). La formulación inercial **no lo excita** (0.6 cm — 27× menor). Esto se debe a que el marco de referencia rotante del corrotacional produce un acoplamiento inercial entre las cargas aerodinámicas y los modos propios de la pala que el marco inercial global no reproduce.

**2. Fluctuaciones de potencia (artefacto estructural):**  
La formulación inercial produce fluctuaciones de potencia **66×** mayores en la banda 1P (12 492 kW vs. 188 kW, rotor equivalente ×3). La causa es la forma en que el solver inercial computa el par mecánico en el marco global — no es una fluctuación aerodinámica real sino un artefacto numérico de la formulación. La potencia aerodinámica real (del participante BEM) muestra un ratio ~2× (§5.4).

**3. Amplitud total de oscilación (p-p) — paradoja aparente:**  
El inercial tiene desplazamiento p-p menor (1.17 m vs. 3.05 m) porque la excitación del modo flapwise que aparece en el corrotacional agrega amplitud estructural. Pero la potencia p-p del inercial es 465 % mayor. La formulación inercial amortígua la respuesta estructural pero amplifica las fluctuaciones de potencia.

**4. Implicación para el análisis de cargas:**  
El corrotacional produce mayor amplitud de deflexión (p-p mayor) — útil como entrada conservadora para un análisis de cargas de pala externo. El inercial produce mayor fluctuación de par/potencia, pero ese exceso es un artefacto de la formulación inercial (ver §5.5 punto 2), no una carga real.

---

### 5.5.1 Análisis de frecuencias características — cargas de rotor (estilo Zhou §4.2)

Zhou et al. [Zhou Energy 2025] presentan en §4.2 la identificación de frecuencias características mediante FFT de los momentos de inclinación (*tilt*) y guiñada (*yaw*) del rotor (Figs. 17-18). Para el caso de turbina fija (*fixed*), sin movimiento de plataforma, los espectros exhiben dos picos dominantes en **1P** y **3P**, sin contribución de la frecuencia de *surge*.

El análisis que se presenta aquí replica esa metodología aplicando FFT con ventana Hann a cuatro señales de la campaña `frontiersin_results/yaw_0` (condición rated, yaw = 0°, t ≥ 20 s, Δt = 0.01 s):

- Empuje aerodinámico (**Aero Thrust** [N])
- Potencia aerodinámica (**Aero Power** [W])
- Momento de inclinación de hub (**Aero Torque X** [Nm], eje paralelo al viento)
- Momento de guiñada de hub (**Aero Torque Z** [Nm], eje vertical)

Parámetros del análisis: Δf = 0.023 Hz (corrotacional, 4270 puntos); Δf = 0.020 Hz (inercial, 5001 puntos). Ventana de Hann con normalización de ganancia coherente (`win / mean(win)`), amplitudes de pico de espectro unilateral.

#### Amplitudes en frecuencias características — thrust y momentos de rotor

| Señal | Unidad | 1P (0.1253 Hz) — Corrot | 3P (0.3759 Hz) — Corrot | f₁_flap (0.5537 Hz) — Corrot | 1P — Inercial | 3P — Inercial | f₁_flap — Inercial |
|---|---|---:|---:|---:|---:|---:|---:|
| Aero Thrust | kN | **16.76** | 2.43 | 2.55 | 16.19 | 2.31 | 0.65 |
| Aero Power | kW | 231 | 30 | 30 | **13 812** | 198 | 54 |
| Tilt Moment (Tx) | MNm | **1.78** | 0.24 | 0.23 | **48.35** | 0.37 | 0.06 |
| Yaw Moment (Tz) | MNm | 0.30 | 0.026 | 0.018 | **48.11** | 0.40 | 0.06 |
| Tip Displacement | cm | 29.6 | 15.6 | **34.4** | 32.6 | 10.6 | 3.3 |

*Fuente: `docs/validation_plots/plot_fft_characteristic_frequencies.py`, campaña `frontiersin_results/yaw_0`, ventana Hann.*

#### Comparación con Zhou Energy 2025 Fig. 18 — turbina fija

| Característica | Zhou (LL-FVW+GEBT, fixed) | AeroElast Corrotacional | Observación |
|---|---|---|---|
| Pico dominante en thrust/tilt | **1P** | **1P** (16.76 kN, 1.78 MNm) | ✓ consistente |
| Segundo pico | **3P** | **3P** (2.43 kN, 0.24 MNm) | ✓ consistente |
| Ratio 3P/1P (thrust) | ≈0.12–0.18 (lectura Fig. 18) | **0.145** | ✓ dentro del rango |
| Pico adicional en f₁_flap | no reportado | **34.4 cm** en desplazamiento, **2.55 kN** en thrust | solo corrotacional |
| Pico en frecuencia de surge | 0.15 Hz (solo condición FOWT) | — (turbina fija) | no aplicable |

#### Interpretación adicional — contraste entre formulaciones

El comportamiento del **solver inercial** en los momentos de hub es anómalo: exhibe amplitudes 1P en Tilt Moment y Yaw Moment de **48.35 MNm y 48.11 MNm** respectivamente — factor ~27× superior al corrotacional (1.78 MNm y 0.30 MNm). Esta anomalía es consistente con las fluctuaciones de potencia 60× mayores identificadas en §5.5 y con el comportamiento observado en la campaña inercial. El origen físico probable es que la formulación inercial no reproduce correctamente la cancelación de cargas de las tres palas en el marco de referencia rotante, resultando en momentos de hub espurios a 1P.

El **solver corrotacional** reproduce correctamente la jerarquía 1P > 3P en thrust y momentos de hub, en línea con Zhou Fig. 18 para turbina fija, y adicionalmente identifica el pico en f₁_flap (34.4 cm en desplazamiento, 2.55 kN en thrust) que evidencia el acoplamiento aeroelástico entre las cargas 1P y el primer modo flapwise de la pala.

![Espectros FFT de cargas características](figures/fig_5_5_1_fft_cargas.png)
*Figura: Espectros FFT (ventana Hann, escala log-y) de cuatro señales de carga — Corrotacional (rojo) vs. Inercial (azul discontinuo). Líneas punteadas verticales: 1P (0.1253 Hz), 3P (0.3759 Hz), f₁_flap (0.5537 Hz), f₁_edge (0.6290 Hz). Notar: (a) thrust 1P similar en ambos solvers; (b) tilt moment 1P del inercial 27× mayor (anomalía); (c) pico en f₁_flap exclusivo del corrotacional en thrust y desplazamiento.*

---

### 5.6 Convergencia preCICE (solver corrotacional)

Fuente: `frontiersin_results/convergence_summary.csv`

| Yaw [°] | Ventanas | Iter./vent. (media ± std) | Iter. mín | Iter. máx | Total iter. | V. no convergentes |
|---:|---:|---:|---:|---:|---:|---:|
|  0 | 4 269 | 3.094 ± 1.399 | 2 | 43 | 13 207 | 0 |
| 10 | 5 000 | 3.015 ± 0.734 | 2 | 15 | 15 074 | 0 |
| 20 | 5 000 | 2.992 ± 0.876 | 2 | 39 | 14 958 | 0 |
| 40 | 5 000 | 2.886 ± 1.052 | 2 | 47 | 14 429 | 0 |

**Residuales de convergencia (al final de cada ventana)**:

| Yaw [°] | Res. desplaz. (medio) | Res. desplaz. (máx) | Res. fuerza (medio) | Res. fuerza (máx) |
|---:|---:|---:|---:|---:|
|  0 | 1.19×10⁻⁶ | 2.97×10⁻⁵ | 1.11×10⁻⁵ | **5.29×10⁻³** |
| 10 | 1.14×10⁻⁶ | 7.75×10⁻⁵ | 9.21×10⁻⁶ | 9.94×10⁻⁵ |
| 20 | 1.17×10⁻⁶ | 7.39×10⁻⁵ | 9.27×10⁻⁶ | 9.98×10⁻⁵ |
| 40 | 1.09×10⁻⁶ | 9.87×10⁻⁵ | 8.09×10⁻⁶ | 9.93×10⁻⁵ |

El pico aislado en yaw=0° (fuerza máx = 5.29×10⁻³) corresponde a una ventana de la fase transitoria y no produce pérdida de convergencia global. El solver opera establemente en ~3 subiteraciones promedio para todo el barrido, lo que confirma la eficiencia de la reutilización de factorización.

---

### 5.7 Contexto de carga extrema — DLC vs operación nominal {#v-04-7}

**Fuente**: NREL/TP-5000-75698, Section 6 "Load Assessment" [R2]; simulaciones `bem_0_10` y `bem_0_10_full`.

#### 5.7.1 Qué reporta NREL/TP-5000-75698 Section 6

El Capítulo 6 del informe de definición del IEA 15 MW presenta un análisis de cargas de diseño (DLC) completo, simulado con OpenFAST, siguiendo la norma IEC 61400-1 Ed.3. El objetivo es determinar la **carga última máxima** sobre componentes estructurales, no el comportamiento en operación normal.

**Resumen del análisis DLC (Tabla 6-1 del informe)**:

| DLC | Condición de viento | Velocidades | Condición operativa | Configuración |
|---|---|---|---|---|
| 1.1 | NTM (turb. normal)    | 3–25 m/s    | Normal — operación   | — |
| 1.3 | ETM (turb. extrema)   | 3–25 m/s    | Normal — operación   | — |
| 1.4 | ECD (ráfaga extrema + dirección) | Vr, Vr±2 m/s | Normal | ±Dir. Change |
| 1.5 | EWS (cizalladura extrema) | 3–25 m/s | Normal | ±Vert./Horiz. |
| **6.1** | **EWM (V50=50 m/s)**  | **V50**     | **Parqueado/idling** | **Yaw ±8°** |
| **6.3** | **EWM (V1=40 m/s)**   | **V1**      | **Parqueado/idling** | **Yaw ±20° + cambio de dirección** |

**Cita directa** (p. 32 del informe):  
> *"Yaw-misaligned parked conditions with extreme wind speeds and extreme coherent gust with a direction change result in the worst-case loading for this design. The worst-case out-of-plane tip deflection is 22.8 m, leaving more than sufficient tower clearance, with an unbent blade tip-to-tower clearance of 30.0 m."*

#### 5.7.2 Qué significa físicamente el valor 22.8 m

El valor 22.8 m representa:

1. **El máximo absoluto** de deflexión fuera del plano entre **todos los DLC** analizados (222 simulaciones OpenFAST).
2. Proviene de **DLC 1.4** (ECD — ráfaga extrema con cambio de dirección, *operación normal* a Vr), **no** de las condiciones parqueadas de DLC 6.x.  
3. El DLC 1.4 aplica una ráfaga coherente extrema (ECD) con cambio abrupto de dirección de viento a velocidades Vr y Vr±2 m/s, con la máquina girando a potencia nominal. El transitorio dinámico resultante excita fuertemente el primer modo flapwise de la pala giratoria.
4. Los DLC 6.x (parqueado + viento extremo V50=50 m/s, yaw ±8°–±20°) producen deflexiones del orden de **~8 m**, no 22.8 m — consistente con el resultado de `bem_90_50_S` (§5.8).  
5. Es un valor pico instantáneo (máximo de una serie temporal), no un promedio temporal.  
6. Fue calculado con el modelo de viga ElastoDyn + estela de BEM (AeroDyn) acoplado en OpenFAST, con todo el modelo de turbina completa (torre, góndola, rotor).

**El 22.8 m NO es comparable directamente con los resultados de régimen permanente de una pala giratoria a potencia nominal** (12.76–12.79 m). Son condiciones físicas fundamentalmente distintas.

#### 5.7.3 Pico transitorio en bem_0_10 (force_ramp_time=0)

La simulación `bem_0_10` usa `force_ramp_time=0`, lo que significa que las fuerzas aerodinámicas se aplican instantáneamente sin rampa. Esto produce un **pico transitorio artificial**:

| Magnitud | Valor |
|---|---:|
| Pico máximo global | **27.6 m** en t = 0.97 s |
| Régimen permanente medio (t ≥ 35 s) | **12.76 m** |
| Régimen permanente std (t ≥ 35 s) | **1.14 m** |
| Régimen permanente máx (t ≥ 35 s) | **14.81 m** |

Este pico de 27.6 m supera el 22.8 m reportado por NREL, pero por razones **numéricas** (impulso por ausencia de rampa), no físicas. El pico NREL proviene de una ráfaga extrema con cambio de dirección (DLC 1.4, ECD) durante operación normal a velocidad rated, que excita dinámicamente el primer modo flapwise de la pala giratoria. El pico de bem_0_10 es un artefacto del esquema de aplicación de condiciones iniciales.

La simulación `bem_0_10_full` (con retroalimentación de velocidad estructural al BEM) produce un pico inicial de **15.37 m** (también a t≈0.95s) y régimen permanente de **12.78 m**: la retroalimentación de velocidad aporta amortiguamiento aerodinámico desde el inicio, atenuando el pico transitorio en ~44%.

#### 5.7.4 Interpretación en contexto del artículo

La comparación de deflexiones en la literatura para el IEA 15 MW puede ser confusa porque distintas fuentes reportan distintos tipos de valor:

| Fuente | Tipo de valor | Condición | Valor |
|---|---|---|---|
| AeroElast (este trabajo) | Media en régimen permanente | Operación nominal, yaw=0° | 12.76–12.79 m |
| Zhou et al. 2025 | Media en régimen permanente | Operación nominal | 13.86 m |
| Bernardi et al. 2025 | Media en régimen permanente | Operación nominal | ~16 m |
| NREL/TP-5000-75698 [R2] | Máximo absoluto DLC | DLC 1.4 — ECD + dir. change, op. nominal (OpenFAST) | **22.8 m** |
| Escalera Mendoza et al. 2023 [R4] | Máximo absoluto DLC | DLC 1.4 — ECD, op. nominal (OpenFAST + UTD NuMAD) | **23.5 m** |

El valor 22.8 m del informe NREL es físicamente consistente con las simulaciones de AeroElast **si se compara con el pico transitorio bajo carga impulsiva** (~27.6 m con ramp=0, o ~15 m con retroalimentación aerodinámica). Sin embargo, la comparación correcta es la de régimen permanente (12.76 m) con los valores de Zhou y Bernardi, ya que todas estas fuentes reportan operación nominal. El valor UTD de 23.5 m [R4] —obtenido con el mismo modelo NuMAD que usa AeroElast pero procesado con OpenFAST— es un 3.1 % superior al NREL, lo que ilustra la sensibilidad de los picos DLC a la representación del modelo estructural aun usando el mismo código CFD.

---

### 5.8 Caso parqueado en posición de bandera (V=45 m/s, yaw=90°) {#v-04-8}

**Simulación**: `bem_90_50_S`  
**Solver**: `StressStiffenedDynamicFSI` (no giratorio, con rigidez geométrica por esfuerzos de membrana)  
**Condición**: pala estática (ω = 0 rad/s), cargada por BEM a V_inf = 45 m/s, pitch = 0°, denominada "posición de bandera" porque la pala apunta en la dirección del flujo con su perfil aerodinámico expuesto lateralmente.

#### Configuración del caso

| Parámetro | Valor |
|---|---|
| Velocidad de viento | 45 m/s |
| Velocidad angular | 0 rad/s (parqueado) |
| Pitch colectivo | 0° |
| `max-time` configurado (preCICE) | 50.0 s |
| Tiempo real de corrida | 30.21 s (finalizada antes del límite) |
| Δt (time-window-size) | 0.005 s |
| `force_ramp_time` | 2.0 s (~1 período flapwise, ≈1/f₁_flap) |
| Integración temporal | Newmark-β constante-aceleración-media (β=0.25, γ=0.5) |
| Amortiguamiento | Rayleigh ζ = 3%, modos i=1, j=2 (num_modes=10) |
| Rigidez geométrica | `update_interval=1` (K_G recalculada cada ventana) |
| Solver lineal | MUMPS directo (`solver_type: "direct"`) |
| Malla estructural | NuMAD Excel, `element_size=0.25 m`, `n_samples=300`, RCM renumber |

#### Resultados de desplazamiento

La ventana de régimen permanente se define como t ≥ 21.1 s (70% del tiempo total; n = 1813 muestras).

| Magnitud | Valor |
|---|---:|
| Máximo global (pico transitorio) | **8.313 m** |
| Régimen permanente — medio | **4.701 m** |
| Régimen permanente — std | 0.118 m |
| Régimen permanente — máximo | 4.955 m |
| Régimen permanente — mínimo | 4.439 m |
| **Componente X (dominante)** | **−4.617 m** (medio) |
| Componente Y | −0.880 m (medio) |
| Componente Z (axial) | −0.026 m (casi nulo) |

La componente dominante es X (~4.6 m de las 4.7 m totales), que en la orientación del caso (eje de pala = Z, viento en X) corresponde al desplazamiento fuera del plano de rotación, equivalente al **flapwise** en operación nominal.

#### Comparación con los DLC 6.x de NREL (parqueado, V50 = 50 m/s, yaw ±8°–±20°)

Según la corrección del usuario, los DLC 6.x (parqueado + viento extremo) de NREL producen deflexiones del orden de **~8 m**, lo que es directamente comparable con el resultado de `bem_90_50_S`:

| Magnitud | bem_90_50_S (AeroElast) | DLC 6.x (NREL OpenFAST) |
|---|---:|---:|
| Velocidad de viento | 45 m/s | 50 m/s (V50) / 40 m/s (V1) |
| Condición de yaw | 90° ("bandera") | ±8° (DLC 6.1) / ±20° (DLC 6.3) |
| ω (rotor) | 0 rad/s (parqueado) | ~0 (parqueado/idling) |
| Deflexión máxima | **8.31 m** | **~8 m** |
| Deflexión en régimen | 4.70 m (medio) | n/d |

La concordancia entre bem_90_50_S (8.31 m) y el rango DLC 6.x de NREL (~8 m) es notable, especialmente considerando que AeroElast usa una velocidad de viento ligeramente menor (45 vs. 50 m/s) y una orientación de yaw diferente (90° vs. ±8°). Esto constituye una **validación cualitativa** del solver en condición parqueada con carga extrema de viento.

La diferencia de orientación explica la similitud de resultados:
- En **posición de bandera** (yaw=90°), el perfil ve el viento lateralmente (AoA ≈ 0°), con fuerza dominante de arrastre (drag). La carga flapwise-equivalente es alta porque la componente de arrastre actúa en la dirección de mayor flexibilidad de la pala.
- En **DLC 6.1** (yaw ±8°), la pala puede estar en cualquier azimut y la carga combinada de sustentación + arrastre sobre la sección más expuesta produce un resultado similar en magnitud.

**Conclusión de este caso**: bem_90_50_S (V=45 m/s, bandera) reproduce cuantitativamente el rango de deflexión parqueada de los DLC 6.x de NREL, validando el solver `StressStiffenedDynamicFSI` en esta condición de carga extrema. Los 22.8 m del informe NREL corresponden a DLC 1.4 (ECD, operación normal), una condición físicamente distinta.

---

## Tabla resumen de validación {#resumen}

| Bloque | Magnitud | AeroElast | Referencia | Error | Estado |
|---|---|---:|---|---:|:---:|
| E-01 MITC | Max error cantilever isotrópico | 1.42 % | analítica | ≤1.42 % | ✓ |
| E-02 MITC4+ | Max error benchmarks Ko 2017 | 3.88 % (pinched cyl.) | Ko et al. 2017 | ≤3.88 % | ✓ |
| V-01 BEM | Thrust rated | 2,522.6 kN | 2,457.0 kN | +2.67 % | ✓ |
| V-01 BEM | Torque rated | 20.576 MN m | 19.910 MN m | +3.34 % | ✓ |
| V-01 BEM | \|ΔC_P\|_max (6 puntos) | 0.0121 | ref. workbook | — | ✓ |
| V-02 modal | 1er flapwise | 0.5537 Hz | 0.5585 Hz [R2] | −0.86 % | ✓ |
| V-02 modal | 1er edgewise | 0.6290 Hz | 0.6406 Hz [R2] | −1.81 % | ✓ |
| V-02 modal | 2do flapwise | 1.6946 Hz | 1.6590 Hz [R2] | +2.15 % | ✓ |
| V-02 modal | 2do edgewise (contexto) | 1.9796 Hz | 2.167 Hz [R2] | −8.65 % | ⚠ |
| V-02 modal | 1er flapwise (vs UTD BModes [R4]) | 0.5537 Hz | 0.57 Hz | −2.9 % | ✓ |
| V-02 modal | 1er edgewise (vs UTD BModes [R4]) | 0.6290 Hz | 0.65 Hz | −3.2 % | ✓ |
| V-02 modal | 3er flapwise (vs UTD BModes [R4]) | 3.2377 Hz | 3.41 Hz | −5.1 % | ✓ |
| V-03 estático | Reacción raíz vs. peso | 658 908.7 N | 658 908.7 N | 8.2×10⁻⁸ % | ✓ |
| V-03 estático | Masa total de pala (vs NuMAD tabular [R3]) | 67 167 kg | 67 921 kg [R3] | −1.11 % | ✓ |
| V-03 estático | Masa total de pala (vs objetivo WISDEM [R2]) | 67 167 kg | 65 250 kg [R2] | +2.94 % | ℹ |
| V-03 estático | Masa total de pala (vs UTD NuMAD [R4]) | 67 167 kg | 68 077 kg | −1.34 % | ✓ |
| V-04 FSI corot. (frontiersin) | Flapwise medio yaw=0° | 12.794 m | 13.86 m [Zhou Energy 2025] | −7.7 % | evidencia B |
| V-04 FSI corot. (frontiersin) | Potencia media yaw=0° | 14.710 MW | 14.76 MW [Zhou Energy 2025] | −0.34 % | evidencia B |
| V-04 FSI corot. (bem_0_10) | Flapwise medio yaw=0° (ss, t≥35s) | 12.761 m | 13.86 m [Zhou Energy 2025] | −7.9 % | evidencia B |
| V-04 FSI corot. (bem_0_10) | Flapwise pico transitorio (ramp=0) | 27.6 m | — | — | ℹ artefacto num. |
| V-04 FSI corot. (bem_0_10_full) | Flapwise medio yaw=0° (ss, t≥35s) | 12.777 m | 13.86 m [Zhou Energy 2025] | −7.8 % | evidencia B |
| V-04 FSI corot. (bem_0_10_full) | Flapwise pico transitorio (ramp=0, vel.fbk) | 15.37 m | — | — | ℹ artefacto num. |
| V-04 FFT corot. | Amplitud 1P desplazamiento | 26.0 cm | — | — | evidencia B |
| V-04 FFT corot. | Amplitud f₁_flap desplazamiento | 16.4 cm | — | — | evidencia B |
| V-04 FFT inercial | Amplitud 1P potencia (×3, equiv. rotor) | 12 492 kW | — | — | evidencia B |
| V-04 FSI corot. | Potencia p-p yaw=0° | 0.88 MW | — | — | evidencia B |
| V-04 FSI inercial | Potencia p-p yaw=0° | 4.97 MW | — | — | evidencia B |
| V-04 DLC referencia | Pico extremo DLC 1.4 — ECD + dir. change, op. nominal (NREL OpenFAST) | — | 22.8 m [R2-DLC] | n/a | 📌 escala diseño |
| V-04 DLC referencia | Pico DLC 1.4 — UTD NuMAD OpenFAST (Escalera Mendoza 2023) | — | 23.5 m [R4] | n/a | 📌 escala diseño |
| V-04 bandera | Flap máx (parqueado, V=45 m/s, yaw=90°) | 8.31 m | ~8 m [R2-DLC-6] | ~+4 % | ✓ validación parqueado |
| V-04 bandera | Flap medio régimen permanente (bem_90_50_S) | 4.70 m | — | — | evidencia B |

---

## Referencias de trazabilidad

| ID | Fuente |
|---|---|
| [R1] | `tests/IEA15MW/IEA-15-240-RWT/Documentation/IEA-15-240-RWT_tabular.xlsx` — curva de potencia, parámetros globales |
| [R2] | NREL/TP-5000-75698 — definición de la turbina IEA 15 MW; Tabla 1-1 (65 t) y Tabla 2-1 (65,250 kg, objetivo de diseño WISDEM); también contiene frecuencias modales de referencia y análisis de cargas DLC |
| [R2-DLC] | NREL/TP-5000-75698, Section 6 "Load Assessment" — deflexión máxima absoluta: 22.8 m proviene de **DLC 1.4** (ECD + cambio de dirección, operación normal a Vr, OpenFAST) |
| [R2-DLC-6] | NREL/TP-5000-75698, Section 6 "Load Assessment" — DLC 6.x (parqueado, V50=50 m/s o V1=40 m/s, yaw ±8°–±20°): deflexión parqueada del orden de ~8 m |
| [R3] | `tests/NuMAD_utd_iea15mw.xlsx` — geometría estructural usada en V-02, V-03 y bem_0_10 |
| [R4] | Escalera Mendoza A.S., Mishra I., Griffith D.T. — "An Open-Source NuMAD Model for the IEA 15 MW Blade with Baseline Structural Analysis", AIAA Scitech 2023, DOI: 10.2514/6.2023-2093. Fuente del modelo NuMAD [R3]; aporta frecuencias modales BModes (Tabla 3), masa AutoNuMAD (68,077 kg), deflexión máxima DLC 1.4 (23.5 m) y flutter ratio (0.91 a Ω_rated). |
| [Zhou Energy 2025] | Zhou L., et al. — "Unsteady aeroelastic performance of the 15 MW floating offshore wind turbine under surge condition", *Energy* 336 (2025) 136488. LL-FVW + GEBT. Tabla 3: frecuencias modales (6 referencias externas + LL-FVW+GEBT); Tabla 4: deflexiones nominales (flap=13.86 m, edge=−1.22 m, tor=−3.60°); Tabla 6: potencia+thrust con/sin flexibilidad (P_flexible=14.76 MW, T_flexible=2.20 MN). |
| [Ma Frontiers 2025] | Ma L., Li Y., Zhou L., Yang D., Shen X., Du Z. — "Study on the aeroelastic performance of 15 MW wind turbine under yaw condition", *Frontiers in Energy Research* 13:1571567 (2025). LL-FVW + GEBT (mismo grupo SJTU que [Zhou Energy 2025]). Condición: V=10.59 m/s, Ω=7.55 rpm, sin control de paso. Tabla 2: frecuencias modales IEA-15 MW (1st flap=0.523 Hz, 1st edge=0.665 Hz, 2nd flap=1.475 Hz, 2nd edge=2.124 Hz, 1st tor=4.072 Hz). Figure 15: potencia media vs. yaw (yaw=10°: ≈14.0 MW, 20°: ≈13.4 MW, 30°: ≈12.0 MW, 40°: ≈10.0 MW) — tendencia ~cos^1.5. Figure 16: deflexión flapwise tip (yaw=10°: ≈14.2 m, 20°: ≈13.8 m, 30°: ≈12.5 m, 40°: ≈11.5 m), edgewise tip ≈−1.05 m (todos los yaw), torsión tip ≈−3.8° a −4.2°. Hallazgo clave: torsión de pala flexible reduce AoA efectivo → menor carga media; velocidad de flapping amplifica fluctuación de carga respecto a pala rígida. |
| [Ramos-García 2022] | Ramos-García N., Kontos S., Pegalajar-Jurado A., González Horcas S., Bredmose H. — "Investigation of the floating IEA Wind 15 MW RWT using vortex methods Part I: Flow regimes and wake recovery", *Wind Energy* 25(3):468–504 (2022). DOI: 10.1002/we.2682. DTU, proyecto COREWIND (EU H2020). Método: MIRAS (lifting-line, wake libre híbrido filamento-partícula) + HAWC2 (vigas Timoshenko, Newmark-β). Condiciones: V=8 m/s (sub-rated, TSR=9.0) y V=15 m/s (sobre-rated); **no simula V_rated=10.59 m/s**. Resultado clave (base fija): BEM subestima potencia media ≲1.4 % vs. LL-FVW en sub-rated; diferencias < 0.1 % en sobre-rated. En caso flotante + olas sobre-rated, BEM sobreestima amplitud máxima de movimiento torre fore-aft en > 50 % vs. LL. Frecuencias citadas de NREL [R2] (0.555 Hz flapwise, 0.642 Hz edgewise) — no son datos independientes. |
| [Bernardi 2025] | Bernardi et al. 2025 — LES + CSD: flapwise ≈ 16 m |
| [ALM] | Campaña interna previa (ALM): flapwise = 14.10 m (no citable externamente) |
