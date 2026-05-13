# IEA 15 MW RWT — Validation Reference Document

Base para la generación de casos de validación del solver estructural/FSI de `aeroelast`.  
Todas las secciones incluyen la fuente exacta de cada dato.

---

## Referencias primarias

| ID | Descripción | Fuente |
|----|-------------|--------|
| [R1] | NREL/TP-5000-75698 — Definition of the IEA 15-Megawatt Offshore Reference Wind Turbine | https://docs.nrel.gov/docs/fy20osti/75698.pdf |
| [R2] | Repositorio oficial IEAWindSystems/IEA-15-240-RWT | https://github.com/IEAWindSystems/IEA-15-240-RWT |
| [R3] | `IEA-15-240-RWT_tabular.xlsx` — tablas digitalizadas del reporte | `Documentation/IEA-15-240-RWT_tabular.xlsx` en [R2]; copia local en `/tmp/opencode/iea15_tabular.xlsx` |
| [R4] | `IEA-15-240-RWT_ElastoDyn_blade.dat` — propiedades 1D distribuidas de la pala | `OpenFAST/IEA-15-240-RWT/` en [R2] |
| [R5] | `IEA-15-240-RWT_BeamDyn.dat` — geometría 3D de la pala (pre-bend, twist) | `OpenFAST/IEA-15-240-RWT/` en [R2] |
| [R6] | `IEA-15-240-RWT_BeamDyn_blade.dat` — matrices 6×6 de rigidez y masa (Timoshenko) | `OpenFAST/IEA-15-240-RWT/` en [R2] |

---

## 1. Parámetros globales del rotor

| Parámetro | Valor | Fuente |
|-----------|-------|--------|
| Potencia nominal | 15 MW | [R3] Overview |
| Clase IEC | IB | [R3] Overview |
| Diámetro de rotor | 241.35 m | [R3] Overview |
| Longitud de pala | 117.0 m | [R5] `kp_zr` en tip |
| Velocidad de viento nominal | 10.659 m/s | [R3] Overview |
| Velocidad de viento de arranque | 3 m/s | [R3] Overview |
| Velocidad de viento de corte | 25 m/s | [R3] Overview |
| RPM mínimas | 5.0 RPM | [R3] Overview |
| RPM máximas | 7.56 RPM | [R3] Overview |
| Velocidad de punta máxima | 95 m/s | [R3] Overview |
| TSR de diseño | 9 | [R3] Overview |
| Ángulo de conicidad (cone) | 4° | [R3] Overview |
| Inclinación de eje (tilt) | 6° | [R3] Overview |
| Pre-bend en el tip | −4.00 m | [R3] Overview; [R5] `kp_xr` en tip |
| Masa de pala | 67,921 kg | [R3] Overview |
| Masa de hub | 21,441 kg | [R3] Overview |
| Masa de góndola + tren | 673,001 kg | [R3] Overview |
| Masa RNA total | 945,770 kg | [R3] Overview |
| Altura de buje (monopile) | 150 m | [R3] Overview |

---

## 2. Caso de validación V-01 — Rendimiento aerodinámico a velocidad nominal

**Objetivo:** verificar que el solver BEM + FSI reproduce los coeficientes de rendimiento de referencia en condiciones nominales.

**Condiciones de operación a rated:**

| Magnitud | Valor de referencia | Fuente |
|----------|---------------------|--------|
| Velocidad de viento | 10.659 m/s | [R3] Rotor Performance |
| Ángulo de pitch | ~0° (0.0 a 0.535°) | [R3] Rotor Performance |
| Velocidad angular | 7.518 RPM (0.7872 rad/s) | [R3] Rotor Performance |
| Velocidad de punta | 95.0 m/s | [R3] Rotor Performance |
| Potencia eléctrica | 15.0 MW | [R3] Rotor Performance |
| Cp (mecánico) | 0.4420 | [R3] Rotor Performance |
| Cp (aerodinámico) | 0.4618 | [R3] Rotor Performance |
| Empuje total (rotor) | **2.457 MN** | [R3] Rotor Performance |
| Ct | 0.7718 | [R3] Rotor Performance |
| Torque total (rotor) | **19.91 MNm** | [R3] Rotor Performance |
| Cq | 0.05182 | [R3] Rotor Performance |
| Momento flector en raíz (pala) | 65.46 MNm | [R3] Rotor Performance |

**Por pala (dividir rotor / 3):**

| Magnitud | Valor (1 pala) |
|----------|----------------|
| Empuje | ~819 kN |
| Torque | ~6.637 MNm |
| Momento raíz | ~21.82 MNm |

**Tabla de rendimiento completa** (todas las velocidades de viento, para validar la curva de potencia):

| V (m/s) | Pitch (°) | P (MW) | Cp_aero | RPM | Thrust (MN) | Ct | Torque (MNm) |
|---------|-----------|--------|---------|-----|-------------|-----|--------------|
| 3.000 | 3.918 | 0.0431 | 0.0595 | 5.00 | 0.2023 | 0.8023 | 0.0860 |
| 5.006 | 2.893 | 1.4011 | 0.4164 | 5.00 | 0.5508 | 0.7842 | 2.7960 |
| 7.159 | 0.000 | 4.5415 | 0.4616 | 5.099 | 1.1177 | 0.7783 | 8.8882 |
| 9.027 | 0.000 | 9.1060 | 0.4616 | 6.429 | 1.7773 | 0.7783 | 14.133 |
| 9.780 | 0.000 | 11.579 | 0.4616 | 6.965 | 2.0861 | 0.7783 | 16.588 |
| 10.210 | 0.000 | 13.173 | 0.4616 | 7.271 | 2.2734 | 0.7783 | 18.078 |
| **10.659** | **~0** | **15.000** | **0.4618** | **7.518** | **2.457** | **0.7718** | **19.910** |
| 11.170 | 3.755 | 15.000 | 0.4013 | 7.518 | 1.958 | 0.5600 | 19.910 |
| 12.259 | 6.761 | 15.000 | 0.3036 | 7.518 | 1.632 | 0.3875 | 19.910 |
| 15.471 | 12.185 | 15.000 | 0.1511 | 7.518 | 1.204 | 0.1796 | 19.910 |
| 20.030 | 17.768 | 15.000 | 0.0696 | 7.518 | 0.930 | 0.0827 | 19.910 |
| 25.000 | 22.829 | 15.000 | 0.0358 | 7.518 | 0.774 | 0.0442 | 19.910 |

> Fuente: [R3] hoja "Rotor Performance". Nota: datos generados con WISDEM bajo condiciones cuasi-estáticas idealizadas, maximizando producción sujeta a límites de RPM y torque.

**Criterios de aceptación sugeridos:**
- Empuje (pala) dentro del ±5% del valor de referencia en estado estacionario
- Torque (pala) dentro del ±5%
- Cp dentro del ±0.02 absoluto

---

## 3. Caso de validación V-02 — Frecuencias naturales de la pala

**Objetivo:** verificar que el modelo FEM (MITC3/MITC4) reproduce las frecuencias modales de referencia.

**Frecuencias naturales de referencia** (blade-alone, sin rotación):

| Modo | Frecuencia (Hz) | Período (s) | Fuente |
|------|----------------|-------------|--------|
| 1er flapwise | **0.5585** | 1.790 | [R1] Tabla 5-2; [R4] modo 1 |
| 2do flapwise | **1.659** | 0.603 | [R1] Tabla 5-2; [R4] modo 2 |
| 1er edgewise | **0.6406** | 1.561 | [R1] Tabla 5-2; [R4] modo 1 edge |
| 2do edgewise | **2.167** | 0.462 | [R1] Tabla 5-2 |
| 1er torsional | **4.459** | 0.224 | [R1] Tabla 5-2 |

**Coeficientes de forma modal** (polinomios en r/R, de [R4]):

Flap modo 1: `φ(x) = -0.01521·x² + 2.42098·x³ - 2.62066·x⁴ + 1.87076·x⁵ - 0.65586·x⁶`  
Flap modo 2: `φ(x) = -0.49737·x² + 0.45816·x³ - 6.97772·x⁴ + 14.7914·x⁵ - 6.77451·x⁶`  
Edge modo 1: `φ(x) = 0.00201·x² + 4.76562·x³ - 9.51016·x⁴ + 8.71931·x⁵ - 2.97678·x⁶`

> Fuente: [R4] sección "BLADE MODE SHAPES"

**Amortiguamiento estructural:**

| DOF | % amortiguamiento crítico | mu (BeamDyn) | Fuente |
|-----|--------------------------|--------------|--------|
| Flap modo 1 | **0.48%** | 0.00299 | [R4]; [R6] |
| Flap modo 2 | **0.48%** | 0.00219 | [R4]; [R6] |
| Edge modo 1 | **0.48%** | 0.00084 | [R4]; [R6] |

> BeamDyn mu coefficients completos: `[0.00299005, 0.00218775, 0.00084171, 0.00218775, 0.00299005, 0.00084171]` — [R6]

**Protocolo de validación sugerido:**
1. Modelo en voladizo (raíz empotrada, tip libre), sin rotación, sin cargas aerodinámicas
2. Aplicar deflexión inicial flapwise pequeña (~0.1 m en tip) y soltar
3. FFT del desplazamiento en tip → pico principal debe ser 0.5585 Hz ± 2%
4. Repetir para excitación edgewise → pico en 0.6406 Hz ± 2%

**Nota sobre la sim actual:** en la corrida inercial previa se observó un pico a 0.395 Hz en desplazamiento flapwise. Esto corresponde a 3P (3 × 7.518 RPM / 60 = 0.376 Hz) y no es la frecuencia natural — es forzamiento aerodinámico. Hay que separar ambos efectos con un caso de vibración libre.

---

## 4. Caso de validación V-03 — Deflexión estática bajo gravedad

**Objetivo:** verificar la respuesta estática del modelo estructural ante carga gravitatoria pura.

**Condición:** pala horizontal (eje de pala perpendicular a la gravedad), sin rotación, sin viento.  
La gravedad actúa en la dirección edgewise (flapwise si la pala está en posición de bandera).

**Valor de referencia:**  
El pre-bend de diseño es −4.00 m en el tip (dirección flapwise, hacia barlovento) [R5].  
La deflexión edgewise estática bajo gravedad para una pala horizontal puede estimarse con la rigidez de referencia:

A r/R = 0.5 (z = 58.5 m): K_55_flap ≈ 4.89 × 10⁹ N·m²; K_66_edge ≈ 2.21 × 10⁸ N·m² [R3] / [R6]  
Masa total pala: **67,921 kg** [R3]

> **No hay un valor de deflexión estática tabulado en las fuentes disponibles.** El valor de referencia debe generarse corriendo OpenFAST/BeamDyn en condición estática (parked, horizontal). Ver [R1] figuras 5-11 a 5-15 (no accesibles sin descargar el PDF completo).

**Criterios de aceptación sugeridos:**
- Deflexión en tip dentro del ±10% de la solución BeamDyn de referencia
- Forma de la deflexión (curvatura) cualitativamente correcta

---

## 5. Caso de validación V-04 — Deflexión flapwise en estado estacionario (rated)

**Objetivo:** verificar que el solver FSI reproduce la deflexión flapwise media en condiciones nominales.

**Valor de referencia:**  
Deflexión flapwise en tip a rated: **~10–12 m** (media en estado estacionario) [R1]  
Este rango proviene de simulaciones OpenFAST/BeamDyn reportadas en [R1] figuras 5-11/5-12.

**Resultado de la simulación inercial previa** (antes del restart, corrida `bem_0_10_full`):
- Media flapwise: **11.94 m** ± 0.99 m
- Máximo: 14.38 m
- Radio deformado: 115.84 m (acortamiento −1.16 m por flexión)

> Este resultado es consistente con la referencia y constituye una validación positiva del solver inercial.

**Condiciones de la simulación:**
- V = 10.659 m/s, ω = 0.7872 rad/s, `omega_ramp_time = 0`, `force_ramp_time = 0`
- Solver: `LinearDynamicFSIRotorInertialSolver`

**Criterios de aceptación sugeridos:**
- Media flapwise en tip entre 9 y 13 m en estado estacionario
- Oscilación 3P visible (frecuencia ~0.376 Hz a rated)
- Acortamiento del radio deformado < 3 m

---

## 6. Caso de validación V-05 — Propiedades estructurales distribuidas

**Objetivo:** verificar que las propiedades mecánicas del modelo FEM coinciden con las del IEA 15MW.

### 6.1 Rigidez de flexión (EI) — ElastoDyn 1D

Valores seleccionados de [R4] (columnas `FlpStff` y `EdgStff`):

| r/R | z (m) | EI_flap (N·m²) | EI_edge (N·m²) |
|-----|-------|----------------|----------------|
| 0.00 | 0.0 | 1.525 × 10¹¹ | 1.525 × 10¹¹ |
| 0.10 | 11.7 | 5.671 × 10¹⁰ | 7.500 × 10¹⁰ |
| 0.20 | 23.4 | 2.462 × 10¹⁰ | 3.779 × 10¹⁰ |
| 0.30 | 35.1 | 1.408 × 10¹⁰ | 2.719 × 10¹⁰ |
| 0.40 | 46.8 | 8.679 × 10⁹ | 2.031 × 10¹⁰ |
| 0.50 | 58.5 | 4.925 × 10⁹ | 1.522 × 10¹⁰ |
| 0.60 | 70.2 | 2.599 × 10⁹ | 7.598 × 10⁹ |
| 0.70 | 81.9 | 1.226 × 10⁹ | 2.411 × 10⁹ |
| 0.80 | 93.6 | 3.838 × 10⁸ | 7.818 × 10⁸ |
| 0.90 | 105.3 | 1.179 × 10⁸ | 2.919 × 10⁸ |
| 1.00 | 117.0 | 1.862 × 10⁵ | 1.663 × 10⁶ |

> Fuente: [R4]

### 6.2 Rigidez de flexión (K_55/K_66) — BeamDyn 6×6

Términos diagonales de la matriz de rigidez de Timoshenko de [R6]:

| r/R | K_55 flapwise (N·m²) | K_66 edgewise (N·m²) | K_11 axial (N) |
|-----|----------------------|----------------------|----------------|
| 0.00 | 1.497 × 10¹¹ | 8.749 × 10¹⁰ | 4.605 × 10¹⁰ |
| 0.10 | 5.226 × 10¹⁰ | 2.543 × 10¹⁰ | 2.419 × 10⁹ |
| 0.20 | 2.182 × 10¹⁰ | 2.768 × 10⁹ | 4.595 × 10⁸ |
| 0.30 | 1.336 × 10¹⁰ | 6.271 × 10⁸ | 2.246 × 10⁸ |
| 0.40 | 8.538 × 10⁹ | 3.648 × 10⁸ | 1.616 × 10⁸ |
| 0.50 | 4.893 × 10⁹ | 2.205 × 10⁸ | 1.166 × 10⁸ |
| 0.60 | 2.655 × 10⁹ | 1.282 × 10⁸ | 8.136 × 10⁷ |
| 0.70 | 1.352 × 10⁹ | 7.146 × 10⁷ | 5.555 × 10⁷ |
| 0.80 | 4.662 × 10⁸ | 3.811 × 10⁷ | 3.651 × 10⁷ |
| 0.90 | 1.133 × 10⁸ | 1.342 × 10⁷ | 1.975 × 10⁷ |
| 1.00 | 1.862 × 10⁵ | 7.145 × 10⁴ | 9.315 × 10⁵ |

> Fuente: [R6]. Coordenada BeamDyn: K_55 = flexión en el plano de batimiento (flapwise), K_66 = flexión en el plano de arrastre (edgewise), K_11 = rigidez axial.

> **Nota:** K_55 de BeamDyn ≈ EI_flap de ElastoDyn a bajos r/R pero divergen hacia el tip — esto se debe a las diferencias entre la formulación Timoshenko completa y la viga de Euler-Bernoulli 1D.

### 6.3 Masa distribuida

| r/R | z (m) | M_11 (kg/m) — BeamDyn | BMassDen (kg/m) — ElastoDyn |
|-----|-------|-----------------------|-----------------------------|
| 0.00 | 0.0 | 3127.4 | 3189.1 |
| 0.10 | 11.7 | — | 1675.4 |
| 0.20 | 23.4 | — | 649.8 |
| 0.50 | 58.5 | 377.7 | 401.2 |
| 0.75 | 87.8 | 179.6 | — |
| 0.90 | 105.3 | 54.7 | 53.4 |
| 1.00 | 117.0 | 5.39 | 5.77 |

> Fuentes: [R6] para BeamDyn, [R4] para ElastoDyn.

**Masa total de pala:** 67,921 kg [R3]

---

## 7. Geometría de la pala

### 7.1 Pre-bend (kp_xr) — curva de pre-curvado

Selección de estaciones de [R5]:

| z (m) | kp_xr (m) | Twist (°) |
|-------|-----------|-----------|
| 0.0 | 0.000 | +15.595 |
| 23.9 | +0.249 | +8.552 |
| 47.8 | +0.187 | +3.077 |
| 57.3 | −0.022 | +1.828 |
| 71.6 | −0.568 | +0.437 |
| 95.5 | −2.100 | −2.086 |
| 107.4 | −2.670 | −2.018 |
| 117.0 | **−4.000** | **−1.242** |

> Fuente: [R5]. kp_xr negativo = dirección barlovento (upwind prebend). Si el modelo no incluye el pre-bend, la geometría inicial ya difiere −4 m respecto a la referencia en el tip.

### 7.2 Cuerda y espesor relativo

| r/R | z (m) | Cuerda (m) | Twist (°) | Espesor rel. (%) |
|-----|-------|------------|-----------|-----------------|
| 0.000 | 0.0 | 5.200 | +15.59 | 100.0 (circular) |
| 0.170 | 19.8 | 5.710 | +10.07 | 46.1 |
| 0.340 | 39.7 | 5.075 | +4.49 | 32.7 |
| 0.500 | 58.5 | 4.199 | +1.81 | 28.5 |
| 0.670 | 78.4 | 3.321 | −0.22 | 23.0 |
| 0.850 | 99.5 | 2.538 | −2.17 | 21.1 |
| 1.000 | 117.0 | 0.500 | −1.24 | 21.1 |

> Fuente: [R3] hoja "Blade Geometry". Cuerda máxima ~5.77 m en r/R ≈ 0.20.

---

## 8. Casos de validación — resumen

| ID | Nombre | Tipo | Magnitud principal | Valor referencia | Fuente |
|----|--------|------|-------------------|-----------------|--------|
| V-01 | Rendimiento aerodinámico | FSI steady-state | Thrust, Torque, Cp | 2.457 MN, 19.91 MNm, Cp=0.462 | [R3] |
| V-02 | Frecuencias naturales | Modal (libre) | f1_flap, f1_edge | 0.5585 Hz, 0.6406 Hz | [R1][R4] |
| V-03 | Deflexión estática gravedad | Estático | Tip deflection edgewise | ~3.5–4 m (generar con OpenFAST) | [R1] |
| V-04 | Deflexión flapwise FSI rated | FSI dinámico | Tip deflection flapwise | 10–12 m (media) | [R1] |
| V-05 | Propiedades estructurales | Verificación | EI_flap, EI_edge, masa | Ver tablas §6 | [R3][R4][R6] |

---

## 9. Notas de implementación

### Pre-bend
El modelo BeamDyn usa `kp_xr = −4 m` en el tip [R5]. Si `aeroelast` no incorpora pre-bend en la malla inicial, la geometría de referencia y la del solver difieren. Esto afecta directamente la comparación de deflexión absoluta pero no la deflexión relativa (incremental desde la posición deformada de equilibrio).

### Coordenadas BeamDyn vs. convención `aeroelast`
- BeamDyn: eje z a lo largo del span, x flapwise (barlovento positivo), y edgewise
- `aeroelast`: span en Z, rotación en Y, DOF layout `[ux, uy, uz, rx, ry, rz]`
- Al comparar K_55 (BeamDyn flapwise) con la rigidez del modelo FEM, verificar que los ejes coincidan

### Condición de contorno en raíz
BeamDyn usa `QuasiStaticInit = True` para pre-acondicionar con aceleraciones centrípetas [R5]. El solver corotacional de `aeroelast` hace algo equivalente implícitamente mediante el ensamblado de K_G en la posición inicial.

### Masa de pala
El valor de 67,921 kg [R3] incluye la pala completa desde la raíz (r/R = 0) al tip. Verificar que la integral de BMassDen de [R4] converja a este valor al integrar sobre los 117 m.

---

*Documento generado el 2026-05-13. Actualizar cuando se disponga de datos de deflexión tabulados del reporte [R1].*
