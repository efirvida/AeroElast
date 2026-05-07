# Informe de Mejoras de Rendimiento
## `LinearDynamicFSIRotorCorotationalSolver`

**Fecha:** 2026-05-06
**Alcance:** Identificación de oportunidades de mejora de rendimiento en el solver corotacional de FSI para rotores.
**Disclaimer epistémico importante:** Este informe distingue entre tres categorías:
- ✅ **Verificado en código** — el comportamiento se confirmó leyendo la fuente.
- 📊 **Probable cuello de botella** — análisis estructural sugiere que es costoso, pero no se midió.
- 💡 **Code smell / oportunidad teórica** — vale la pena considerar pero el impacto real depende del perfil.

> **NO HAY DATOS DE PROFILING EN ESTE INFORME.** Las recomendaciones están priorizadas por análisis estático del código, no por mediciones reales. La recomendación #0 es ejecutar un profiler antes de invertir esfuerzo significativo en cualquier otra optimización.

---

## 0. Recomendación #0 — Medir antes de optimizar

**ANTES de cualquier cambio, ejecutar profiling sobre una simulación representativa:**

```bash
# Profiling Python (callback _step_cb)
py-spy record --output fsi_profile.svg --duration 120 --rate 250 -- \
    python -m aeroelast.cli run_fsi --config caso_test.yaml

# O cProfile para detalles más finos
python -m cProfile -o fsi.prof -m aeroelast.cli run_fsi --config caso_test.yaml
snakeviz fsi.prof
```

**Métricas a verificar:**
1. ¿Cuánto tiempo (%) se gasta en `_step_cb` vs solver Rust vs preCICE I/O?
2. ¿Cuántas veces se invoca `rotation_matrix(theta)` por callback?
3. ¿Cuántas refactorizaciones de K_eff ocurren por simulación (cuenta los `refactorize` en logs)?
4. ¿Cuánto tarda el `_compute_stress_fields` cuando se llama?
5. Tamaño de la malla (n_dofs) y mesh de interfaz (n_iface_nodes).

**Sin estos datos, los porcentajes de impacto en este informe son estimaciones de orden de magnitud.**

---

## 1. Arquitectura — Lo que ya está bien

Antes de listar problemas, conviene reconocer las decisiones de diseño que ya capturan los wins más importantes:

### ✅ Factorización de K_eff cacheada en PETSc

Verificado en `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs:492-496`:

> *"K_eff is constant for fixed Δt, so this pays the O(n^α) factorization cost ONCE and reuses only O(n) back-substitution on every step() call."*

El `NewmarkDynamicStepper` mantiene una `KSP` PETSc cacheada con `KSPSetOperators` y solo refactoriza cuando K_G, K_SP o dt cambian. **Este es el win individual más grande posible para sistemas implícitos** y ya está implementado.

### ✅ Histéresis adaptativa para K_G

Verificado en `rotor_fsi.rs:425-499`. El rebuild de K_G usa doble threshold:
- `THRESHOLD_REBUILD = 0.005` (alto, dentro de los primeros 10 pasos)
- `THRESHOLD_SKIP = 0.003` (bajo, después)

Esto evita "chattering" de refactorizaciones cuando ω está cerca de un umbral. Diseño cuidadoso.

### ✅ Rust hace el trabajo pesado, no el callback Python

El loop interno (sub-iteraciones preCICE, ensamblaje de fuerzas, paso Newmark) corre en Rust. El callback Python `_step_cb` solo se invoca **una vez por ventana convergida**, no por sub-iteración. Esto evita el GIL en el hot loop.

### ✅ Stress recovery diferido

`_compute_stress_fields` se llama solo cuando hay checkpoint o probes activos (`rotor.py:1876-1885`). Si se mantiene `stress_output_interval > 1`, no afecta el time-stepping.

---

## 2. Optimizaciones de Bajo Riesgo — Quick Wins (1-3 días)

Estas son cambios localizados con bajo riesgo de regresión. ROI individual modesto pero acumulan.

### 2.1 Cachear R(θ) dentro de `_step_cb` 💡

**Verificado:** En `_step_cb` se invocan **8-9 transformaciones de coordenadas** que recomputan internamente la matriz `R(θ)` con la fórmula de Rodrigues:

```python
# rotor.py — todas estas llamadas son con el MISMO theta dentro de un callback:
self._coord_transforms.to_rotating(f_aero_nodes_global, theta)        # línea 1896
self._coord_transforms.to_rotating(iface_force_global, theta)         # línea 1902
self._coord_transforms.to_rotating(self._gravity, theta)              # línea 1928
self._coord_transforms.to_inertial(f_inertial_nodes_local, theta)     # línea 1987
self._coord_transforms.to_inertial(f_gravity_nodes_local, theta)      # línea 1989
self._coord_transforms.to_inertial(f_total_nodes_local, theta)        # línea 1992
# + 2-3 más dentro de _compute_axis_torque
```

Cada llamada ejecuta `rotation_matrix()` (`corotational.py:123`) que hace `np.cos`, `np.sin`, y dos sumas de matrices 3×3.

**Cambio propuesto:** Computar R(θ) y R^T(θ) una vez al inicio del callback, y exponer un método `to_rotating_with_R(vec, R)` que reciba la matriz precomputada.

**Impacto estimado:** 📊 ~10-30% del tiempo del callback (depende de proporción transforms vs otros costos). El impacto absoluto depende de cuánto pesa el callback en el total — por eso #0 es profiling.

**Riesgo:** Bajo. Los cambios son localizados.

### 2.2 Vectorizar extracción de masa nodal ✅

**Verificado** en `rotor.py:1737-1743`:

```python
all_node_masses = np.array([
    float(diag_arr[i * dofs]) if i * dofs < len(diag_arr) else 0.0
    for i in range(n_nodes)
], dtype=np.float64)
```

Loop Python sobre todos los nodos. Se ejecuta una vez al inicio, no en hot loop, pero es trivial de mejorar:

```python
expected = n_nodes * dofs
if len(diag_arr) >= expected:
    all_node_masses = diag_arr[::dofs][:n_nodes].copy().astype(np.float64)
else:
    all_node_masses = np.zeros(n_nodes, dtype=np.float64)
    valid_n = len(diag_arr) // dofs
    all_node_masses[:valid_n] = diag_arr[:valid_n * dofs:dofs]
```

**Impacto:** Negligible para tiempo total (es startup). **Por qué hacerlo igual:** Es código limpio y el patrón actual sugiere que el autor desconoce el slicing strided de numpy — vale la pena alinear.

### 2.3 Vectorizar scatter de DOFs en interfaz 💡

**Verificado** en `rotor.py:1871-1874`:

```python
iface_u_local = np.zeros_like(interface_coords_nodes)
for local_idx in range(min(len(iface_dofs_flat), iface_u_local.size)):
    gdof = int(iface_dofs_flat[local_idx])
    if gdof < n_total:
        iface_u_local.flat[local_idx] = u_full[gdof]
```

Loop Python sobre `n_iface_nodes * 3`. Reemplazable con:

```python
mask = iface_dofs_flat < n_total
valid_dofs = iface_dofs_flat[mask].astype(np.int64)
iface_u_local.flat[:len(valid_dofs)] = u_full[valid_dofs]
```

**Impacto:** Pequeño en mallas con interfaz pequeña (<1000 nodos). Crece linealmente con tamaño de interfaz.

### 2.4 Pre-allocar buffers del callback 💡

**Verificado:** `_step_cb` allocata varios arrays en cada llamada (líneas 1855, 1860, 1863, 1870, 1888, 1906, 1926). Con el callback ejecutándose una vez por ventana convergida, esto NO causa "TB de allocations" como sugería el análisis preliminar — el GC de Python libera los buffers entre callbacks.

**Pero:** las allocations repetidas presionan al allocator, fragmentan memoria y disparan GC. Para una ventana convergida pesada (n_total ~ 1M DOFs), `np.zeros(n_total)` es ~8 MB.

**Cambio propuesto:** Allocar buffers persistentes en `__init__` y reusarlos:

```python
def __init__(self, ...):
    ...
    self._cb_buf_u_full = None  # Lazy alloc on first callback when n_total known
    self._cb_buf_v_full = None
    self._cb_buf_a_full = None
    # ...

def _step_cb(self, ...):
    if self._cb_buf_u_full is None:
        self._cb_buf_u_full = np.zeros(n_total)
        # ...
    self._cb_buf_u_full.fill(0)
    self._cb_buf_u_full[_free_dofs_arr] = u_red
    # ...
```

**Impacto:** Marginal por callback, acumula en simulaciones largas (>10k pasos). **Riesgo:** Hay que cuidar thread-safety si se paraleliza el callback en el futuro.

---

## 3. Optimizaciones de Alto Impacto — Requieren verificación

Estas tienen potencial de impacto significativo pero **requieren medición o cambio arquitectural**.

### 3.1 ⚠️ DEFINITIVAMENTE VALE LA PENA: Eliminar recomputación de fuerzas en Python

**Verificado:** El bloque de líneas 1906-1933 en `_step_cb` recomputa F_centrifuga, F_coriolis, F_euler y F_gravedad **en Python**, con el único propósito de logging y reporte de torques:

```python
f_inertial_nodes_local = np.zeros_like(all_node_coords_nodes)
if self._include_centrifugal:
    f_inertial_nodes_local += self._inertial_calculator.compute_centrifugal_force(...)
if self._include_coriolis:
    f_inertial_nodes_local += self._inertial_calculator.compute_coriolis_force(...)
if self._include_euler and abs(alpha_window) > 1e-14:
    f_inertial_nodes_local += self._inertial_calculator.compute_euler_force(...)

f_gravity_nodes_local = np.zeros_like(all_node_coords_nodes)
if self._include_gravity:
    gravity_local = self._coord_transforms.to_rotating(self._gravity, theta)
    f_gravity_nodes_local = all_node_masses[:, np.newaxis] * gravity_local
```

**El problema:** Rust ya computó estas mismas fuerzas en `rotor_fsi.rs:625-680` para resolver el paso estructural. Recomputarlas en Python es duplicación de trabajo.

**Cambio propuesto:** Modificar la signatura del callback de Rust para devolver los arrays de fuerzas (`f_centrifugal`, `f_coriolis`, `f_euler`, `f_gravity`) ya en el frame rotante. El callback Python simplemente las consume.

**Esto requiere:**
1. Modificar `crates/aeroelast-py/src/lib.rs:3534+` (la función `run_rotor_fsi_solver`).
2. Que Rust mantenga los buffers de fuerzas accesibles (no destructivos).
3. Modificar la signatura `_step_cb` para aceptar las nuevas fuerzas.

**Por qué es la recomendación más sólida del informe:** Es un win arquitectural claro independiente de profiling. Estás ejecutando el mismo cálculo dos veces. **Eso siempre es desperdicio.**

**Impacto estimado:** 📊 ~15-25% del tiempo del callback. Mayor en mallas grandes (cálculos vectoriales O(n_nodes)).

**Riesgo:** Medio. Toca interfaz Python ↔ Rust. Necesita tests de regresión cuidadosos para verificar que las fuerzas reportadas coinciden con la versión actual.

### 3.2 ⚠️ Verificar frecuencia real de refactorización 📊

**Lo que sabemos verificado:**
- Cada update de K_G dispara `refactorize()` → O(n^α)
- Cada update de K_SP dispara `refactorize()` → O(n^α)
- `kg_update_interval` default = 0 → normalizado a 1 → **rebuild evaluado cada paso**
- La histéresis omega_sq (`THRESHOLD_REBUILD = 0.005`, 0.5%) bloquea el rebuild si ω no cambia >0.5%
- `ksp_omega_threshold` default = `1e-4 rad/s`

**Lo que NO sabemos sin profiling:** Cuántas refactorizaciones realmente ocurren en una corrida típica.

**Casos a evaluar empíricamente:**

| Escenario | ω comportamiento | Predicción de refactorizaciones |
|---|---|---|
| Operación nominal (ConstantOmega) | ω = const | Solo refactoriza al inicio. **Óptimo.** |
| Arranque (RampedOmega) | ω rampa de 0 → ω_target | K_G: probablemente cada paso (Δω/ω grande al principio). K_SP: cada paso (threshold 1e-4 muy bajo). **Costoso.** |
| Operación dinámica (ComputedOmega) | ω varía con torque | Variable. Depende de magnitud de oscilación. |

**Acción:** Añadir un contador de refactorizaciones en `dynamic_newmark.rs:refactorize()` y loggear el total al final de la simulación. Si el conteo es alto en operación nominal, hay un bug; si es alto durante rampa, considerar:

1. Subir `ksp_omega_threshold` (e.g. a `1e-3` o `1e-2 rad/s`) durante rampa.
2. Usar la misma estrategia de histéresis adaptativa que ya tiene K_G, también para K_SP.

**Impacto estimado:** Si las refactorizaciones son frecuentes, el ahorro puede ser **dramático** (la factorización suele ser >50% del costo del paso para sistemas medianos a grandes). Si no son frecuentes, este punto es moot. **Por eso requiere medir primero.**

### 3.3 Reducir conversiones PETSc → COO al inicio 💡

**Verificado** en `rotor.py:1709-1721`:

```python
k_rows, k_cols, k_vals = self._petsc_to_coo(self.K)
m_rows, m_cols, m_vals = self._petsc_to_coo(self.M)
if K_G is not None:
    kg0_rows, kg0_cols, kg0_vals = self._petsc_to_coo(K_G)
```

Cada `_petsc_to_coo` hace `getValuesCSR` + expansión `np.repeat` (CSR → COO) y luego dtype casting. Es O(nnz). Se ejecuta una sola vez al inicio.

**Impacto:** Solo afecta tiempo de startup. Para sistemas grandes (~10⁷ nnz), puede ser segundos. **Para una simulación de 1000+ pasos, irrelevante.** No optimizar a menos que se necesite arranque rápido (e.g., baterías de tests).

---

## 4. Mejoras Arquitecturales — Largo plazo

### 4.1 Mover parte del callback a Rust

**Hipótesis:** Si el profiling muestra que `_step_cb` consume una fracción significativa del tiempo (>20%), considerar mover toda la lógica de logging y métricas a Rust:

- Cálculo de torques (`_compute_axis_torque`)
- Cálculo de coeficientes (`_compute_performance_coefficients`)
- Escritura de CSV de performance

**Beneficios:**
- Elimina el roundtrip GIL Python ↔ Rust
- Permite paralelizar el callback con el siguiente paso de tiempo (asincronía)

**Costos:**
- Duplicación de lógica si el código Python sigue existiendo
- Reducción de flexibilidad para usuarios que quieran customizar logging

**Decisión:** Solo justificable si el profiling lo demanda.

### 4.2 Eliminar `_extract_nodal_translation_field` redundante

`u_full → u_nodes_local` y `v_full → v_nodes_local` se hacen en cada callback. Para mallas con DOFs rotacionales (cáscaras MITC4 con 5-6 DOF/nodo), esto descarta información. Si el callback solo necesita traducciones, podría exponerse Rust pasando ya el campo nodal sin la necesidad de "expandir luego contraer".

### 4.3 GPU acceleration

**No recomendado** para este tamaño de problema típico (n_dofs ~ 10⁴-10⁵). El overhead de transferencia GPU no se amortiza. Solo evaluar si se mueve a problemas de >10⁶ DOFs y >10⁴ nodos en interfaz.

---

## 5. Configuración del Usuario — Impacto Inmediato

Estas son palancas que el usuario puede ajustar **sin tocar código** y deberían documentarse en `cli-reference.md`:

| Parámetro | Default | Recomendación performance | Cuándo aplica |
|---|---|---|---|
| `solver.rotor.kg_update_interval` | `0` (cada paso) | `5-10` para ω casi constante | Operación nominal |
| `solver.rotor.ksp_omega_threshold` | `1e-4 rad/s` | `1e-3` durante rampa | Arranque/transitorios |
| `solver.stress_output_interval` | varía | `≥ checkpoint_interval` | Siempre |
| `coupling.preCICE.acceleration` | IQN-ILS | Verificar `max-used-iterations` | FSI en general |

**Especial atención a:** En operación nominal con ω casi constante, **subir `kg_update_interval` a 5-10 puede ahorrar refactorizaciones significativas** (ver §3.2). El default `0` es seguro pero conservador.

---

## 6. Tabla Resumen de Mejoras

Ordenadas por **ratio impacto estimado / riesgo**:

| # | Mejora | Categoría | Riesgo | Impacto estimado | Requiere medición previa |
|---|---|---|---|---|---|
| 0 | **Ejecutar profiler primero** | Diagnóstico | Nulo | Habilita el resto | — |
| 1 | Eliminar recomputación de fuerzas en `_step_cb` | Arquitectural | Medio | 📊 Alto | No (es waste claro) |
| 2 | Cachear R(θ) dentro de `_step_cb` | Quick win | Bajo | 📊 Medio | Sí (cuánto pesa el callback) |
| 3 | Subir `kg_update_interval` para ω constante | Config | Nulo | 📊 Alto en op. nominal | Recomendado |
| 4 | Histéresis adaptativa para K_SP | Rust internal | Medio | 📊 Alto en transitorios | Sí (contar refactorizaciones) |
| 5 | Pre-allocar buffers en callback | Quick win | Bajo | 📊 Bajo | No |
| 6 | Vectorizar scatter de DOFs en interfaz | Quick win | Bajo | 📊 Bajo | No |
| 7 | Vectorizar extracción de masa nodal | Quick win | Bajo | 📊 Negligible (startup) | No |
| 8 | Mover callback completo a Rust | Arquitectural | Alto | 📊 Variable | **Sí, esencial** |

---

## 7. Plan de Acción Sugerido

### Semana 1: Diagnóstico
- Ejecutar `py-spy record` y `cProfile` sobre 3 casos representativos (ConstantOmega, RampedOmega, ComputedOmega).
- Añadir contadores de refactorización en `dynamic_newmark.rs` y loggear totales.
- Documentar el perfil real.

### Semana 2: Quick wins
- Implementar #2 (cache R(θ)) con tests de regresión que verifiquen invariancia numérica.
- Implementar #5 (pre-alloc) y #6 (vectorize scatter).
- Documentar el parámetro `kg_update_interval` con guía de uso (#3).

### Semana 3-4: Win arquitectural
- Implementar #1 (eliminar recomputación de fuerzas en Python).
- Re-perfilar y verificar mejora.

### Más adelante (condicional al profiling)
- #4 si las refactorizaciones de K_SP resultan frecuentes.
- #8 solo si el callback Python sigue siendo >20% del tiempo después de #1-#7.

---

## 8. Conclusión Honesta

El solver tiene una **arquitectura razonable** con las decisiones más impactantes ya correctas (factorización K_eff cacheada, histéresis adaptativa para K_G, callback Python solo en convergencia). Los quick wins disponibles son reales pero modestos individualmente.

**La única optimización con impacto claro y certero independiente de profiling** es la #1 (eliminar la recomputación de fuerzas en Python para logging) — porque es trabajo duplicado verificado, no especulación sobre cuellos de botella.

**Todo lo demás depende de medir primero.** Sin un profiler, optimizar es apostar.

---

## 9. Referencias al código

- `src/aeroelast/solvers/fsi/rotor.py:1827-2090` — `_step_cb` callback
- `src/aeroelast/solvers/fsi/corotational.py:123-141` — `rotation_matrix`
- `crates/aeroelast-solvers/src/petsc/fsi/rotor_fsi.rs:425-499` — K_G rebuild logic con histéresis
- `crates/aeroelast-solvers/src/petsc/elasticity/dynamic_newmark.rs:441-820` — `NewmarkDynamicStepper` con KSP cacheada
- `crates/aeroelast-py/src/lib.rs:3534-3680` — interfaz Python ↔ Rust del solver
- `src/aeroelast/core/config.py:359` — defaults de parámetros del rotor
