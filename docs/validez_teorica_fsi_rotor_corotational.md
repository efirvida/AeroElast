# Informe de Validez Teórica
## `LinearDynamicFSIRotorCorotationalSolver`

**Fecha:** 2026-05-06
**Alcance:** Análisis crítico de la formulación teórica del solver corotacional de FSI para rotores de turbinas eólicas. Basado en la revisión del código fuente en `src/aeroelast/solvers/fsi/rotor.py`, `src/aeroelast/solvers/fsi/corotational.py` y la documentación de diseño existente.
**Perspectiva:** Análisis aeroelástico de turbinas eólicas, mecánica de marcos no inerciales, FEM de cáscaras rotantes.

---

## 1. Ecuación de Gobierno

El solver resuelve, en el marco rotante:

```
[M]{ü} + [C]{u̇} + ([K] + [K_G] + [K_SP]){u} = {F_aero} + {F_cf} + {F_cor} + {F_euler} + {F_g}
```

Esta ecuación es **correcta** y corresponde a la formulación estándar de mecánica analítica para marcos no inerciales (ANSYS Theory Reference §14.4.1, Eq. 14-57). La descomposición LHS/RHS de los efectos físicos es uno de los puntos más delicados de la formulación. A continuación se evalúa cada término.

---

## 2. Marco Corotacional — Transformaciones de Coordenadas

**Validez: ✅ CORRECTA**

La transformación usada es la fórmula de Rodrigues (`corotational.py:123`):

```
R = I + sin(θ)·K + (1 − cos(θ))·K²
```

Esto es exacto para cualquier ángulo (no es aproximación de pequeño ángulo). Las convenciones de transformación son correctas:

- Fuerzas CFD → marco rotante: `F_local = R^T · F_global` ✅
- Desplazamientos → marco inercial: `u_global = R · u_local` ✅

**Observación menor:** El método `to_rotating` para arrays usa `vec_global @ R` que es equivalente a `(R^T · vec_global.T).T`. La implementación es correcta pero la notación puede confundir a quien compare contra derivaciones row-vector vs column-vector.

---

## 3. Rigidez Geométrica K_G (Stress Stiffening)

**Validez: ⚠️ CORRECTA EN FORMULACIÓN, CON RIESGO EN LA PRESTRESS**

### 3.1 Formulación

K_G se ensambla elemento a elemento:

```
K_G = ∫ B_G^T · S̃ · B_G dA
```

donde S̃ es el tensor de tensión en el plano generado por la carga centrífuga. Esta es la formulación estándar de rigidez geométrica para elementos de cáscara (`rotor.py:20-23`). **Correcto.**

K_G aparece en el LHS → **aumenta** las frecuencias naturales (efecto de tensado de cuerda, stiffening flap). **Correcto.**

### 3.2 El problema: la prestress se estima, no se resuelve

El código utiliza `σ_cf ≈ ρ·ω²·r·L_char`, donde `L_char` es una longitud característica del elemento (`rotor.py:22-23`). Esto es una **aproximación geométrica** que reemplaza el problema estático correcto:

```
Proceso riguroso:
  Paso 1: Resolver [K]{u_cf} = {F_cf}  → obtener σ₀ exacta
  Paso 2: Ensamblar K_G(σ₀)
```

Para una pala con sección transversal variable y torsión geométrica, la distribución real de tensiones puede diferir del estimado `ρ·ω²·r·L_char` por un factor que depende de la geometría local. **El error en K_G se propaga directamente a los cruces del diagrama de Campbell**, donde más importa la precisión de las frecuencias naturales.

**Recomendación:** Para validación modal, implementar el paso estático previo. Para análisis de cargas aeroelásticas en operación nominal, la aproximación es prácticamente aceptable.

### 3.3 K_G con ω variable (ComputedOmega)

El parámetro `kg_update_interval` tiene default `0` (nunca se reconstruye después de la primera vez). Para ω dinámico (arranque, ráfagas), K_G calculado con `ω_inicial` se vuelve obsoleto. El error en K_G es O((Δω/ω)²·K_G). Durante una rampa de arranque donde ω varía 0→12 rpm, K_G está **sistemáticamente subestimado** durante toda la rampa.

**Recomendación:** Habilitar `kg_update_interval > 0` para simulaciones con ω variable.

---

## 4. Spin Softening K_SP

**Validez: ✅ CORRECTA EN FORMULACIÓN, ⚠️ INCOMPLETA EN DOFs ROTACIONALES**

### 4.1 Formulación

```
K_SP = -ω² · M · (I - n̂⊗n̂)
```

Captura la variación de la fuerza centrífuga con el desplazamiento elástico (ANSYS Eq. 3-74). Aparece en el LHS y **reduce** la rigidez en el plano perpendicular al eje de rotación (softening lead-lag). **Formulación correcta.**

### 4.2 La relación F_cf @ X₀ + K_SP·u = F_cf @ X₀+u

El código evalúa F_cf en coordenadas no deformadas X₀ y K_SP está en el LHS. El comentario en `rotor.py:34-37` explica correctamente que:

```
F_cf(X₀ + u) ≈ F_cf(X₀) + ω²·M·(I - n̂⊗n̂)·u = F_cf(X₀) - K_SP·u
```

Por lo tanto, llevar K_SP al LHS con F_cf en X₀ no duplica la fuerza centrífuga. **Matemáticamente impecable** — es la linearización correcta de la fuerza centrífuga con el desplazamiento.

**Advertencia:** Esto es una linearización. Para deflexiones de punta > 5% de la longitud de pala (e.g., una pala de 60 m con 3+ m de deflexión), el residuo de orden superior `O(u²)` no está capturado.

### 4.3 DOFs rotacionales — problema con lumped mass

La mass matrix lumped por row-sum produce **ceros en la diagonal para los DOFs rotacionales** de los elementos de cáscara MITC3/MITC4. Por lo tanto:

```
K_SP = -ω² · M_lumped · (I - n̂⊗n̂)
```

produce **K_SP = 0 para los DOFs θx, θy, θz** de cada nodo. El spin softening actúa únicamente sobre los DOFs traslacionales. Para una cáscara con 5-6 DOFs/nodo, esto significa que aproximadamente la mitad de los términos de K_SP son nulos. Este es un problema estructural de la formulación que no se documenta explícitamente como limitación.

**Impacto práctico:** Subestimación del softening para modos que involucran rotaciones nodales (torsión, modos de cáscara de orden superior).

---

## 5. Fuerzas Inerciales

### 5.1 Fuerza Centrífuga F_cf

**Validez: ✅ CORRECTA**

```
F_cf = m · ω² · r_perp   (corotational.py:452)
```

Evaluada en X₀. Correcta. Ya cubierta en §4.2.

### 5.2 Fuerza de Coriolis F_cor — LA LIMITACIÓN MÁS IMPORTANTE

**Validez: ⚠️ TEÓRICAMENTE INCOMPLETA — RIESGO PARA DIAGRAMAS DE CAMPBELL**

El código trata la fuerza de Coriolis como una **fuerza explícita en el RHS** con velocidad retrasada:

```python
F_cor = -2 · m · (ω × v_{n})   # rotor.py:1914
```

La formulación **exacta** de la ecuación de movimiento en marco rotante coloca el término de Coriolis como una **matriz giroscópica antisimétrica [G] en el LHS**:

```
[M]{ü} + ([C] + [G]){u̇} + ([K]+[K_G]+[K_SP]){u} = {F_ext}

donde [G] = 2·[M]·[Ω_skew]  (Ω_skew = skew-symmetric de ω·n̂)
```

La diferencia teórica tiene dos consecuencias:

#### a) Splitting giroscópico — NO capturado

Con [G] en el LHS, los pares de modos degenerados en ω=0 se **separan linealmente** con la velocidad angular:

```
ω_forward  = ω_natural + ω_rotor · Γ_i
ω_backward = ω_natural - ω_rotor · Γ_i
```

donde Γ_i es el factor giroscópico del modo i. Este es el efecto de **forward/backward whirl**. Con Coriolis explícito en RHS, el solver produce diagramas de Campbell **simétricos** — no predice la separación de modos. Para análisis de resonancia (cruces de Campbell a n·rpm), esto puede llevar a errores en la predicción de velocidades críticas.

**El solver subestimará el margen de estabilidad de modos lead-lag**, que son los más afectados por el efecto giroscópico.

#### b) Estabilidad numérica

[G] antisimétrico en LHS no contribuye al trabajo por ciclo (cero work rate). Con Coriolis explícito en RHS, la fuerza sí hace trabajo:

```
δW_Cor ≈ O((ω·Δt)² · KE) por paso de tiempo
```

Para una turbina eólica típica (ω ≈ 1.5 rad/s, dt ≈ 0.01 s):

```
ω·Δt ≈ 0.015  →  error ≈ 0.02% de KE por paso
```

**Esto es estable y aceptable** para turbinas de eje horizontal (HAWT) a velocidades nominales. La condición de estabilidad práctica es `ω·Δt ≪ 1`, que se cumple para HAWT con dt ≤ 0.05 s. **NO es seguro para rotores de alta RPM** (helicópteros, turbinas pequeñas, etc.).

### 5.3 Fuerza de Euler F_euler

**Validez: ✅ CORRECTA**

```python
F_euler = -m · (α × (X₀ + u))   # rotor.py:1920-1923
```

Evaluada en coordenadas deformadas. Es correcto — a diferencia de F_cf, no existe una corrección LHS para F_euler, por lo que evaluarla en deformadas es la aproximación más consistente disponible. **Bien documentado en `rotor.py:39-41`.**

### 5.4 Gravedad en Marco Rotante

**Validez: ✅ CORRECTA**

```python
gravity_local = R^T(θ) · g_global   # rotor.py:1928
```

La gravedad es un vector constante en el marco inercial. Al rotarla al marco local produce la **carga 1P** (once-per-revolution) que es fundamental en análisis de fatiga de palas. La implementación es correcta.

---

## 6. Integración Temporal — Newmark-β

**Validez: ✅ CORRECTA**

```
β = 0.25, γ = 0.5  (regla trapezoidal / aceleración promedio)
```

Esta es la variante incondicionalmente estable de Newmark para sistemas lineales. La formulación:

```
K_eff = [K] + [K_G] + [K_SP] + a₀·[M] + a₁·[C]
```

es la estándar. Para sistemas con [G] explícito en RHS, el Newmark trapezoidal añade una disipación numérica artificial proporcional a `(ω·Δt)²` en los modos giroscópicos — esto es bajo para turbinas HAWT como se argumenta en §5.2.

**Observación:** La integración de ω (ComputedOmega) usa **Euler explícito** (`rotor.py:107-108`), que es de primer orden. Para rampa lineal de ω durante arranque, Euler reproduce exactamente la rampa, sin error. Para ω variable no lineal (respuesta a ráfaga), el error acumula como O(Δt). Con Δt = 0.01s y variaciones de torque típicas, esto es aceptable.

---

## 7. Acoplamiento FSI — preCICE IQN-ILS

**Validez: ✅ ESTADO DEL ARTE**

El acoplamiento implícito con aceleración IQN-ILS (quasi-Newton) es el método de referencia actual para FSI particionado. La arquitectura es correcta:

- Cada ventana temporal mantiene `ω̄` constante durante las sub-iteraciones → K_SP y K_G no cambian dentro de la ventana. Correcto y eficiente.
- La convergencia del punto fijo es manejada por preCICE, no por iteraciones internas. Correcto.

**Riesgo menor:** Para ω fuertemente variable entre ventanas, el salto de K_SP entre ventanas puede afectar la convergencia del IQN-ILS (cambia la función de mapeo fluido-estructura). El rebuild threshold `_OMEGA_CHANGE_THRESHOLD = 1e-4` parece conservador — podría ajustarse para mejorar performance.

---

## 8. Coeficientes de Rendimiento

**Validez: ✅ FÍSICAMENTE CORRECTOS, ⚠️ COMPARACIÓN CON ESTÁNDAR IEC**

Los coeficientes se calculan usando solo fuerzas aerodinámicas (`rotor.py:150-154`), consistente con IEC 61400-1 y BEM estándar. **Correcto.**

**Problema de trazabilidad:** El radio R usado es el deformado (max distancia perpendicular al eje en configuración deformada). Esto introduce variaciones de Cp/Ct debidas a deflexión estructural, no a aerodinámica. Los mapas de curva de potencia industriales usan R nominal. Para validación contra ensayos o contra diseños de referencia, esta diferencia puede ser confusa. **No es incorrecto físicamente, pero rompe comparabilidad directa.**

---

## 9. Dinámica de ω — OmegaProvider

**Validez: ✅ CORRECTO EN PRINCIPIO**

El balance de momento (torque solo de fuerzas externas, no inerciales — `rotor.py:121-123`) es **termodinámicamente correcto**. Las fuerzas ficticias del marco rotante no producen aceleración angular neta del rotor.

El cálculo de inercia `I = Σᵢ mᵢ · r_⊥,ᵢ²` usando masa lumped es aceptable para una estimación de inercia de la pala. Para el rotor completo, habría que incluir inercias de buje y tren de potencia — si no están incluidas, I subestimado → sobreestimación de α.

---

## 10. Tabla Resumen de Validez

| Componente | Estado | Impacto |
|---|---|---|
| Marco corotacional (Rodrigues) | ✅ Correcto | — |
| Newmark-β (β=0.25, γ=0.5) | ✅ Correcto | — |
| K_G formulación | ✅ Correcto | — |
| K_SP formulación y partición LHS/RHS | ✅ Correcto | — |
| F_cf @ X₀ + K_SP·u (no doble conteo) | ✅ Correcto | — |
| F_euler @ X₀+u | ✅ Correcto | — |
| Gravedad → marco rotante | ✅ Correcto | — |
| Torque solo fuerzas externas | ✅ Correcto | — |
| preCICE IQN-ILS | ✅ Estado del arte | — |
| **K_G prestress ≈ ρω²rL_char** | ⚠️ Aproximación | Frecuencias naturales, Campbell |
| **K_SP nulo en DOFs rotacionales** | ⚠️ Incompleto | Modos torsionales, cáscara |
| **Coriolis explícito (RHS)** | ⚠️ Limitación significativa | **Splitting giroscópico no capturado** |
| **ω por Forward Euler** | ⚠️ 1er orden | ω variable rápida |
| **K_G estático cuando ω varía** | ⚠️ Default peligroso | Simulaciones transitorias |
| **R deformado en Cp/Ct/Cq** | ⚠️ Trazabilidad | Comparación vs IEC |

---

## 11. Veredicto Final

**Para turbinas eólicas de eje horizontal a condiciones nominales (7–15 RPM, palas de rigidez moderada): el solver es teóricamente sólido para análisis de cargas aeroelásticas y respuesta dinámica.**

La implementación corotacional es limpia. El splitting LHS/RHS de efectos físicos está cuidadosamente documentado y es correcto. Newmark trapezoidal es el integrador correcto.

**La limitación dominante** es la Coriolis explícita: los **diagramas de Campbell generados por este solver no predicen el splitting giroscópico** (forward/backward whirl). Para análisis modal aeroelástico o identificación de velocidades críticas, los resultados son conservadoramente incorrectos — las frecuencias aparecerán como degeneradas cuando físicamente se habrían separado.

**Para validación modal:** implementar [G] en el LHS es el paso más importante que está faltando. El propio código lo reconoce como "reservado para el futuro" (`rotor.py:289`).

**Para análisis de cargas en operación estacionaria:** el solver produce resultados correctos dentro de las hipótesis de pequeña deformación en el marco rotante.

---

## 12. Recomendaciones Priorizadas

### Alta prioridad (validez teórica)
1. **Implementar [G] giroscópico en el LHS** — habilita análisis modal correcto y predicción de whirl modes. Requiere solver no simétrico.
2. **K_G desde prestress estática real** — paso previo de equilibrio centrífugo antes del ensamble de K_G. Mejora frecuencias naturales para palas tapered/twisted.

### Media prioridad (precisión numérica)
3. **K_SP con masa consistente o lumping HRZ** — recuperar términos rotacionales de K_SP que actualmente son nulos.
4. **Default `kg_update_interval > 0`** para simulaciones con ω variable.
5. **Integrador RK4 o trapezoidal implícito para ω** — reduce error acumulado en transitorios.

### Baja prioridad (trazabilidad)
6. **Opción de usar R nominal en Cp/Ct/Cq** — facilita comparación con curvas industriales y normativa IEC.

---

## 13. Referencias

- ANSYS Mechanical APDL Theory Reference §3.4 (Stress Stiffening), §3.5 (Spin Softening), §14.4.1 (Rotating Frame Dynamics)
- Shabana, A. A. — *Computational Dynamics*, 3rd ed., §6 (Reference Frames)
- Hansen, M. H. — *Aeroelastic stability analysis of wind turbines using an eigenvalue approach*, Wind Energy, 2004
- Bathe, K. J. — *Finite Element Procedures*, §10.2 (Newmark Method), §6.4 (Geometric Stiffness)
- Bauchau, O. A. — *Flexible Multibody Dynamics*, §17 (Rotating Beams)
- IEC 61400-1 ed. 4 — *Wind turbines – Part 1: Design requirements*

---

*Documento generado a partir de la revisión de:*
- `src/aeroelast/solvers/fsi/rotor.py`
- `src/aeroelast/solvers/fsi/corotational.py`
- `src/aeroelast/solvers/fsi/linear_dynamic.py`
- `docs/teoria_formulacion_fsi_rotor.md`
- `docs/rotor_inertial_solver_design.md`
