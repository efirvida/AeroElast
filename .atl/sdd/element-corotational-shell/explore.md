# Exploration: element-corotational-shell

**Change ID**: element-corotational-shell  
**Date**: 2026-05-12  
**Scope**: Formulación corotacional de elemento para MITC3/MITC4 con soporte en el solver dinámico

---

## Current State

### Lo que hay hoy (evidencia de código)

**MITC4** (`crates/aeroelast-core/src/elements/mitc4.rs`):
- `Mitc4Precomputed` tiene `gp_initial_frames: [GpLocalFrame; N_GAUSS]` (línea 119)
- `update_corotational_frame()` implementado (líneas 2708–2750): calcula el frame deformado desde las coordenadas actuales via tangentes covariant `g_r`, `g_s`, produce `{e1, e2, e3}` actualizados
- `frame_incremental_rotation()` implementado (líneas 2756–2766): `R_inc = R_new · R_old^T`
- `polar_decomposition()` implementado (líneas 2646–2671): **con bug conocido** (ver análisis B)
- **NO existe** `compute_kt_corotational()`
- `compute_kt_global()` usa Total Lagrangian: `K_T = T^T·(K0 + K_L + K_σ)·T` con `T` fijo al frame de referencia (líneas 1330–1374)
- `compute_fint_global()` implementado (línea 1459), nonlinear flag usa B_L + B_NL GL

**MITC3** (`crates/aeroelast-core/src/elements/mitc3.rs`):
- `Mitc3Precomputed` tiene `t3: Matrix3<f64>` (línea 75) — frame local fijo al estado de referencia
- `compute_kt_global()` (línea 878): `K_T = T^T·(K0 + K_L + K_σ)·T` con `T` = `t3` (fijo)
- `compute_fint_global()` implementado (línea 920)
- **NO HAY** ninguna estructura corotacional (ni `gp_initial_frames`, ni `update_corotational_frame`, ni frame deformado)

**Assembler** (`crates/aeroelast-core/src/assembly/assembler.rs`):
- `assemble_k()` (línea 413): K lineal, llama `compute_ke_global()` — ignora desplazamientos
- `assemble_kt()` (línea 848): K_T Total Lagrangian, llama `compute_kt_global(pre, &ue)` — desplazamientos en frame fijo
- `assemble_fint()` (línea 945): F_int, llama `compute_fint_global(pre, &ue, nonlinear)`
- `update_reference()` (línea 276): actualiza las coordenadas de referencia (UL incremental)
- **NO HAY** `assemble_kt_corotational()` ni similares

**NewmarkStepper** (`dynamic_newmark.rs`):
- Almacena `k_vals`, `kg_vals`, `ksp_diag` como COO (líneas 542–562)
- Factoriza `K_eff = K + a0·M + a1·C` UNA VEZ, reutiliza (líneas 570–578)
- `refactorize_count` para diagnóstico (línea 644)
- El stepper NO conoce el assembler: recibe COO valores, no el assembler

**rotor_fsi.rs**:
- `update_kg_if_needed()` (línea 490): actualiza K_G con frecuencia configurable (`should_update_kg_on_step`)
- Llama `refactorize()` del stepper cuando K_G cambia
- Patrón: triggear rebuild de K_G por cambio de ω o cambio de deflexión (`kg_deflection_rebuild_rel_high`)

---

## Análisis A: Formulación corotacional — elección correcta

### Iteración 1

Para palas con >20° de rotación flap, hay tres candidatos. Pensemos desde primeros principios:

**TL + K_G** (estado actual): el frame de referencia NO se mueve. El operador lineal `B_L` se calcula siempre respecto al estado sin deformar. El término K_G = B_NL^T·S·B_NL captura el efecto de pretensión pero **no corrige el frame de medición de strains**. Cuando la rotación real supera ~10–15°, el tensor de deformación Green-Lagrange tiene componentes de rotación rígida que generan "strain falso". El error es O(θ²) donde θ es el ángulo de rotación. Para θ=20°, el error de strain es ~12%.

**Updated Lagrangian (UL) incremental**: actualiza la referencia en cada paso, así los strains de cada incremento son pequeños aún cuando el total sea grande. Es correcto para grandes rotaciones + grandes strains. Sin embargo requiere actualizar la geometría de referencia (llamar `update_reference()`), y dentro de cada paso sigue siendo lineal — la consistencia tangente es solo primera orden en el paso.

**Corotacional puro (Crisfield/Felippa)**: el frame LOCAL sigue la rotación rígida del elemento. K_T_global = T^T(θ)·K_L·T(θ) + K_σ donde T(θ) se actualiza en cada paso a partir de la posición deformada. Los strains se miden siempre en el frame rotado → siempre son pequeños. Correcto para grandes rotaciones + pequeñas strains. Es exactamente el régimen de palas (flap/edge < 3% strain, pero rotaciones > 30°).

### Iteración 2

¿UL no es suficiente? Para UL incremental hay que elegir el tamaño de paso cuidadosamente — si Δθ_paso > ~5° el error lineal dentro del paso acumula. En FSI dinámico los pasos son `dt = precice_time_window`, tipicamente 0.001–0.01s con ω = 1–3 rev/s → Δθ_paso ≈ 0.004–0.1 rad (0.2° – 6°). UL puede funcionar PERO requiere llamar `update_reference()` en cada paso convergido, lo que recalcula `Mitc3Precomputed::new()` y `Mitc4Precomputed::new()` para todos los elementos — esto es caro y destruye el caché del stepper.

El corotacional puro NO requiere actualizar la geometría de referencia. Sólo actualiza T(θ) (una matriz 3×3 por Gauss point, a partir del frame deformado actual) y recomputa `K_T = T^T·K_L·T + K_σ`. K_L es la rigidez local **que ya existe** y es constante (la geometría de referencia no cambia).

### Iteración 3

¿Qué dice el código ya implementado? `update_corotational_frame()` calcula el frame deformado a partir de `current_coords` (las coordenadas actuales de los 4 nodos). `frame_incremental_rotation()` da R_inc. El pattern es correcto para el corotacional puro: las herramientas están listas para MITC4, solo falta conectarlas en `compute_kt_corotational()`.

El UL sería más trabajo y costoso en performance (recalcular todo el Precomputed). El TL + K_G es insuficiente para > 15°.

**Conclusión A**: La opción correcta es **corotacional puro de elemento** (Crisfield/Felippa). Es el match exacto para el régimen de palas:
- Grandes rotaciones (30°+ flap) + pequeñas strains (<3%)
- El código de MITC4 ya tiene las primitivas (`update_corotational_frame`, `frame_incremental_rotation`, `gp_initial_frames`)
- No requiere actualizar la geometría de referencia → compatible con el stepper actual
- La formulación matricial es `K_T_global = T^T·K_L·T + K_σ(σ_local)`

**Riesgo**: La consistencia tangente de K_σ (el término de fuerzas internas al girar el frame) es no trivial para shells MITC — si se omite, la convergencia Newton-Raphson falla, pero como este change usa Newmark lineal (no NR), basta con K_σ consistente para la dinámica.

---

## Análisis B: Polar decomposition — corrección del bug

### Iteración 1

El algoritmo de Denman-Beavers (DB) para `sqrt(C)` es:
```
Y_{k+1} = (Y_k + Z_k^{-1}) / 2
Z_{k+1} = (Z_k + Y_k^{-1}) / 2
```
con `Y_0 = C`, `Z_0 = I`. La secuencia `Y_k → sqrt(C)`, `Z_k → (sqrt(C))^{-1}`.

El código actual (línea 2654–2661):
```rust
let mut u = 0.5 * (&ct + Matrix3::identity());  // Y_0 = (C + I)/2 — INCORRECTO: debería ser C
for _ in 0..3 {
    u = 0.5 * (&u + &ct * &u_inv);  // u = (Y + C·Y^{-1})/2
```

Hay DOS problemas:
1. El estado inicial `u = (C+I)/2` no es el correcto para DB (debería inicializarse en `C`)
2. La iteración `u = (Y + C·Y^{-1})/2` no es DB estándar — en DB el segundo término es `Z_k^{-1}`, no `C·Y^{-1}`. Son equivalentes matemáticamente sólo en la primera iteración cuando `Z_0 = I`, pero luego divergen.

### Iteración 2

¿Hay una alternativa directa? Sí: **SVD**. Para `F = U·Σ·V^T` (SVD), la descomposición polar es `F = R·U_stretch` donde `R = U·V^T` y `U_stretch = V·Σ·V^T`. Esto da R exacto sin iteración y funciona para cualquier magnitud de deformación.

nalgebra tiene `nalgebra::SVD::new(f, true, true)` que devuelve `(u, singular_values, v_t)`.

Caso edge: si `det(F) < 0` (elemento invertido), `det(U·V^T)` puede ser -1. En ese caso R no es una rotación sino una reflexión impropia. La convención es verificar `det(R)` y, si es -1, reemplazar la columna de V correspondiente al singular value mínimo por su negativo (Umeyama correction). Esto sucede cuando el elemento se degrada (colapso nodal o carga extrema que invierte el elemento).

### Iteración 3

¿Para el régimen de palas (pequeñas strains), cuántas iteraciones DB serían suficientes si se corrigiera? Con el DB correcto (Y_0 = C, Z_0 = I), la convergencia es cuadrática. Para strains < 5%: 2 iteraciones dan error de orden (ε/2)^4 ≈ 1.5e-6 — suficiente. Para el DB malo actual: la convergencia es lenta y puede divergir para C con valores propios muy diferentes.

**Conclusión B**: El bug tiene dos componentes:
1. Inicialización incorrecta: `(C+I)/2` en lugar de `C`
2. Iteración incorrecta: no sigue DB estándar correctamente para 3D

**La corrección recomendada es usar SVD** (`nalgebra::SVD`) porque:
- Exacto, sin iteración
- Maneja grandes rotaciones directamente
- Maneja el caso edge (det=-1) con la corrección de columna
- Performance: nalgebra SVD 3×3 es ~50ns — negligible comparado con la integración Gauss

El test `test_polar_decomposition` usa tolerancia `< 0.05` (5%) — tapa el bug. Con SVD, la tolerancia debería ser `< 1e-12`.

**Riesgo**: El caso `det(F) < 0` (elemento invertido) debe detectarse explícitamente para retornar un error en lugar de una rotación silenciosa incorrecta.

---

## Análisis C: Frame corotacional para shells MITC

### El problema del director

El shell MITC tiene 6 DOFs/nodo: `[u, v, w, θx, θy, θz]`. Las rotaciones `(θx, θy)` son rotaciones del **director nodal** (vector normal) alrededor de los ejes locales. `θz` es la rotación drill (rotación del plano de la sección).

El frame corotacional de elemento típicamente sigue la rotación del **plano medio**. Pero los directores pueden rotar independientemente del plano medio (flapping/twisting de la sección transversal). Este es el punto crítico.

### Lo que hace el código actual

`update_corotational_frame()` (línea 2708–2750):
- Computa `g_r = ∂x/∂ξ`, `g_s = ∂x/∂η` desde las posiciones actuales de los 4 nodos (solo traslaciones `u,v,w`)
- `e3 = normalize(g_r × g_s)` — normal al plano medio deformado
- `e1 = normalize(g_r - (g_r·e3)·e3)` — proyección de g_r al plano
- `e2 = e3 × e1`

**No usa los directores (θx, θy) nodales**. Sigue el plano medio.

### ¿Es correcto?

Para el régimen de palas de rotor:
- La deformación dominante es flap (w grande, θx, θy moderados) y edge (u/v) 
- Los directores se rotan por coherencia con el plano medio (no hay deformación transversal de sección)
- La hipótesis de Kirchhoff-Love/Mindlin es que el director sigue el plano medio

El frame que sigue el plano medio deformado (lo que hace el código) es la elección correcta según Crisfield (1990) y Felippa & Haugen (2005) para shells delgados. Para shells gruesos o con gran shear transversal sería necesario usar el director promedio, pero palas de rotor son shells delgados.

Sin embargo, hay una sutileza: al transformar las rotaciones `(θx, θy, θz)` al nuevo frame corotacional, se debe transformar **también las rotaciones**, no solo las traslaciones. El `build_t24()` actual usa `t3` fijo (no el frame deformado).

**Conclusión C**: El frame corotacional siguiendo el plano medio deformado es correcto para el régimen de palas. Pero la transformación `T24` debe actualizarse para usar el frame deformado `{e1_def, e2_def, e3_def}` en lugar del fijo `t3`. Esto afecta tanto las traslaciones como las rotaciones en el vector de DOF.

**Riesgo**: Si se actualiza `T24` para rotaciones también, hay que asegurarse que la transformación de `(θx, θy)` al nuevo frame sea consistente. El error tipico es transformar sólo las traslaciones y dejar las rotaciones en el frame viejo.

---

## Análisis D: K_T corotacional — implicaciones para el assembler

### Lo que necesita el assembler

El `K_T_corotacional_global = T(θ)^T · K_L_local · T(θ) + K_σ`

donde:
- `T(θ)`: transformación 24×24 (o 18×18) construida desde el frame deformado actual `{e1_def, e2_def, e3_def}` — NOT el `t3` de referencia
- `K_L_local`: rigidez lineal en el frame deformado — es `compute_ke_local(pre)` evaluada con la geometría de referencia pero expresada en el frame actual
- `K_σ`: `compute_geometric_stiffness_local()` usando las tensiones actuales

Para llamar `K_T_corotacional` el assembler necesita:
1. Las coordenadas actuales de los nodos (para actualizar el frame)
2. El vector de desplazamiento actual `u` (para calcular las tensiones actuales y K_σ)

Hoy `assemble_kt(&self, u: &[f64])` ya recibe `u`. Lo que falta es:
- Un método en `Mitc4Precomputed`/`Mitc3Precomputed` que dado `u_global` y las coordenadas actuales (= `initial_coords_3d + u_traslaciones`) calcule el frame deformado y produzca `K_T_coro`
- Un nuevo método en el assembler: `assemble_kt_corotational(&self, u: &[f64]) -> (Vec<i64>, Vec<i64>, Vec<f64>)`

**Conclusión D**: El assembler necesita un nuevo método `assemble_kt_corotational()`. La interfaz puede ser idéntica a `assemble_kt()` (solo recibe `u`) porque las coordenadas actuales se derivan de `initial_coords_3d + u_traslaciones`. No hay que cambiar la API hacia afuera — el stepper y el solver pasan `u` igual que antes.

**Riesgo**: El `assemble_kt_corotational()` es inherentemente más caro que `assemble_kt()` porque recalcula el frame para todos los elementos en cada llamada. Sin embargo, este costo ya existe para la opción K_G con deformed coords.

---

## Análisis E: Consistencia con los solvers FSI

### Arquitectura actual del stepper

`NewmarkStepper` almacena `k_vals` como COO fijo (línea 542). La factorización `K_eff = K + a0·M + a1·C` se hace en `refactorize()`, que es caro (factorización LU).

El patrón actual en `rotor_fsi.rs` para K_G: actualiza `kg_vals` y llama `refactorize()` solo cuando se cumplen los predicados (ω-change o deflexión-change). Esto se hace cada N pasos configurable.

### ¿Hay que refactorizar en cada paso con K_T corotacional?

Primer pensamiento: K_T cambia en cada paso porque el frame cambia → hay que refactorizar en cada paso. Eso es enormemente costoso.

Segundo pensamiento: En el régimen de palas de rotor con Newmark implícito, la rigidez efectiva `K_eff = K_T + a0·M + a1·C`. Si K_T cambia suavemente (rotaciones lentas), el cambio relativo entre pasos es O(ω·dt). Para ω=2π rad/s y dt=0.005s → Δθ=0.03 rad → ΔK/K ≈ θ·ΔK_rel ≈ 2% por paso. Se puede usar la misma estrategia que K_G: actualizar K_T corotacional cada N pasos (configurable). 

Tercer pensamiento: Pero hay una diferencia con K_G: K_G es el efecto de prestress (que varía lento con ω), mientras K_T corotacional cambia la **orientación** del frame — si se usa K_T obsoleto en el Newmark, los residuos de fuerza pueden ser incorrectos. La clave es que en el FSI de palas el problema es cuasiestático en el sentido de que la rigidez varía lentamente. Para el límite conservador: refactorizar cada paso (máxima precisión, máximo costo). Para el límite eficiente: refactorizar cada K_update_freq pasos (igual que K_G hoy).

**Conclusión E**: Se puede reutilizar el mecanismo existente de `should_update_kg_on_step()` para controlar también la frecuencia de actualización de K_T corotacional. El impacto en performance depende de la frecuencia elegida:
- Cada paso: ~10x más caro para la factorización (igual que NR completo)
- Cada 10 pasos: ~1.1x más caro — aceptable

La estrategia recomendada: una flag `use_corotational_kt` en `RotorFsiConfig` con un `kt_update_freq` separado (default 1, actualizable). Esto permite al usuario elegir entre precisión y performance.

**Riesgo principal**: Si K_T se actualiza con lag (cada N pasos), la integración puede ser ligeramente inconsistente. Para análisis de frecuencias (flutter) esto podría afectar los resultados. Debe documentarse.

---

## Análisis F: MITC3 — qué necesita específicamente

### Diferencias respecto a MITC4

**Geometría**: El triángulo tiene Jacobiano constante (no varía con ξ, η). El frame local `t3` se calcula en `new()` una sola vez. Los tangentes `g_r`, `g_s` son los vectores lado `p2-p1`, `p3-p1`.

**Corotacional en triángulo**: Para actualizar el frame deformado de un MITC3, el cálculo es más simple que en MITC4 porque el Jacobiano del triángulo es constante — hay un solo frame por elemento (no uno por Gauss point). Esto simplifica la implementación.

La polar decomposition del triángulo usa la misma F = I + H, H = dU/dx. El gradiente de desplazamiento H para el triángulo es:
```
H = [[dUx/dx, dUx/dy, 0], ...]
```
donde dUx/dx se calcula con las derivadas de forma `dh`. En un triángulo plano (shell), H es 2D (en el plano del elemento), no 3D completo. La polar decomposition 3D aplica igual si se completa con la normal.

**El problema**: `Mitc3Precomputed` no tiene `gp_initial_frames` ni ninguna estructura de frame corotacional. Para MITC3 habría que:
1. Agregar campos análogos (pero un solo frame, no N_GAUSS frames porque Jacobiano constante)
2. Implementar `update_corotational_frame_tri()` 
3. Implementar `compute_kt_corotational()` para MITC3

El enfoque es análogo al MITC4 pero simplificado por el Jacobiano constante.

**Conclusión F**: El MITC3 puede seguir el mismo approach que MITC4. La complicación triangular es menor que la cuadrilateral porque el frame es único por elemento. El mayor costo es la implementación desde cero (mientras MITC4 ya tiene las primitivas). Se puede secuenciar: primero MITC4, luego MITC3 en el mismo change o en uno subsiguiente.

---

## Affected Areas

- `crates/aeroelast-core/src/elements/mitc4.rs` — bug polar_decomp (SVD fix) + nuevo `compute_kt_corotational()`
- `crates/aeroelast-core/src/elements/mitc3.rs` — agregar structs corotacionales + `update_corotational_frame_tri()` + `compute_kt_corotational()`
- `crates/aeroelast-core/src/assembly/assembler.rs` — nuevo `assemble_kt_corotational()`
- `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs` — flag `use_corotational_kt`, llamar `assemble_kt_corotational()` en lugar de `assemble_kt()`
- `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs` — potencialmente `update_tangent_stiffness()` para reemplazar k_vals en caliente

---

## Scope definitivo del change

### ENTRA:

1. **Fix bug polar_decomp MITC4** — reemplazar iteración Denman-Beavers por SVD. Fix test con tolerancia realista.
2. **`compute_kt_corotational()` en MITC4** — usa `update_corotational_frame()` + nuevo `build_t24_deformed()` + K_L + K_σ
3. **Primitivas corotacionales para MITC3** — struct de frame, `update_corotational_frame_tri()`, `compute_kt_corotational()`
4. **`assemble_kt_corotational()` en MeshAssembler** — dispatch a los nuevos métodos de elemento
5. **Integración en rotor_fsi** — flag `use_corotational_kt` + `kt_update_freq` + llamada al nuevo assembler
6. **Tests**: test polar_decomp grande rotación (30°, 45°), test K_T corotacional vs TL (debe coincidir a 0° y diferir a 30°)

### NO ENTRA (fuera de scope):

- Newton-Raphson completo (no es el objetivo)
- Updated Lagrangian para grandes strains
- Corotacional para elementos sólidos (HEXA, TETRA)
- Actualización del director independiente del plano medio (shells gruesos)
- Integración en el rotor_inertial solver (puede venir después)

---

## Decisiones de diseño clave para el spec

1. **SVD para polar decomp**: la corrección de bug debe ser SVD directo, no Denman-Beavers corregido. Razón: exactitud y manejo del caso det=-1.

2. **Frame = plano medio deformado**: no el director promedio. Correcto para shells delgados. Documentar en código.

3. **T24_deformed**: el `build_t24()` actual usa `pre.t3` fijo. Hay que crear `build_t24_deformed(frame_def: &GpLocalFrame)` que use el frame actualizado. Las rotaciones (θx, θy, θz) deben transformarse con el mismo frame que las traslaciones.

4. **API assembler**: `assemble_kt_corotational(&self, u: &[f64])` — misma firma que `assemble_kt()`, derivar coords actuales internamente.

5. **Frecuencia de update**: parámetro `kt_update_freq` en `RotorFsiConfig`, default = 1 (actualizar cada paso), con comentario de tradeoff.

6. **MITC3 en el mismo change**: un solo frame por elemento (no N_GAUSS). Simplifica la implementación.

---

## Risks

- **Bug silencioso en la transformación de rotaciones**: Si `build_t24_deformed` transforma solo traslaciones y no rotaciones, los resultados serán silenciosamente incorrectos. Necesita test específico con rotación pura.
- **det(F) < 0**: elementos muy deformados o degenerados. La corrección SVD debe manejar esto explícitamente y loguear advertencia.
- **Inconsistencia K_T / f_int**: Para que Newton-Raphson (si algún día se agrega) sea consistente, `compute_fint_global` debe usar el mismo frame que `compute_kt_corotational`. Si este change agrega solo K_T corotacional sin actualizar f_int, hay inconsistencia. Para Newmark lineal (sin NR) esto es aceptable en este change pero debe documentarse.
- **MITC4 tiene N_GAUSS=4 frames, MITC3 tiene 1**: el diseño del assembler debe ser uniforme. Recomendado: usar el frame del centroide para K_T (un solo frame por elemento es suficiente para la transformación de rigidez).
- **Performance**: `assemble_kt_corotational()` recalcula frames en cada llamada. Con rayon paralelo (ya implementado en `assemble_kt`) esto es aceptable, pero hay que benchmarkear.

---

## Estimación de complejidad

| Tarea | Complejidad | Riesgo |
|-------|-------------|--------|
| Fix polar_decomp → SVD + test | Bajo | Bajo |
| `compute_kt_corotational` MITC4 | Medio | Medio (transformación rotaciones) |
| Primitivas corotacionales MITC3 | Medio | Bajo (Jacobiano constante) |
| `assemble_kt_corotational` | Bajo | Bajo |
| Integración rotor_fsi | Medio | Medio (refactorization trigger) |
| Tests de integración | Medio | Bajo |
| **Total** | **Medium-High** | |

Estimación de implementación: 3–4 días de trabajo enfocado. El riesgo más alto es la corrección en la transformación de rotaciones (punto 3 del scope).

---

## Ready for Proposal

**Sí.** Las decisiones técnicas principales están claras:
- Formulación: corotacional puro de elemento (Crisfield)
- Corrección bug: SVD en polar decomp
- Frame: plano medio deformado
- Scope: MITC3+MITC4 + assembler + rotor_fsi integration
- Performance: update_freq configurable

El próximo paso es `sdd-propose` para formalizar el intent, scope y approach.
