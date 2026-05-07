# Análisis Crítico: Solver Inertial para Rotores FSI

**Autor**: Análisis técnico integral  
**Fecha**: 6 de mayo de 2026  
**Versión analizada**: Commit 0e7b3d8 (implementación completa Phases 0-6)  

---

## 1. Análisis del Plan de Implementación

### 1.1 Estructura del Plan (Phases 0-10)

El plan sigue una estrategia **incremental y validada progresivamente**, con clara separación de responsabilidades por fase:

#### **Fase 0: Infraestructura y Compatibilidad** (6 tareas)
✅ **EVALUACIÓN: EXCELENTE**

**Fortalezas:**
- Rename explícito evita ambigüedad: `LinearDynamicFSIRotorCorotationalSolver` vs `LinearDynamicFSIRotorInertialSolver`
- Alias temporal `LinearDynamicFSIRotorSolver` preserva compatibilidad hacia atrás
- Dispatch centralizado en `runner.py` facilita mantenimiento

**Debilidades:**
- ⚠️ **CRÍTICO**: El alias temporal puede crear confusión en documentación futura
- ⚠️ No hay estrategia de deprecación definida (¿cuándo se elimina el alias?)

**Recomendación:**  
Agregar advertencia de deprecación en logs cuando se use el alias:
```python
if solver_type == "LinearDynamicFSIRotor":
    _logger.warning(
        "LinearDynamicFSIRotor is deprecated. Use explicit "
        "LinearDynamicFSIRotorCorotational or LinearDynamicFSIRotorInertial."
    )
```

---

#### **Fase 1: Utilidades Compartidas** (4 tareas)
✅ **EVALUACIÓN: MUY BUENO**

**Fortalezas:**
- Reutiliza `CoordinateTransforms` del solver corotacional (DRY principle)
- `InertialForcesCalculator` encapsula cálculo de `a_ref = α × r + ω × (ω × r)`
- Validación analítica mediante casos con trayectorias circulares conocidas

**Debilidades:**
- ⚠️ **MODERADO**: No hay tests unitarios automatizados para estas utilidades
- ⚠️ El cálculo de `F_ref = -M·a_ref` asume M diagonal (lumped mass)
  - Si se usa masa consistente en el futuro, este cálculo es INCORRECTO
  - Falta validación explícita: `assert self.M.getType() == 'diagonal'`

**Recomendación:**  
Agregar guard en `_compute_reference_load_vector`:
```python
if not self._is_mass_lumped():
    raise NotImplementedError(
        "Inertial solver requires lumped mass matrix. "
        "Consistent mass with -M·a_ref is not yet implemented."
    )
```

---

#### **Fase 2: Skeleton del Solver** (4 tareas)
✅ **EVALUACIÓN: BUENO**

**Fortalezas:**
- Herencia de `LinearDynamicFSISolver` reutiliza infraestructura probada
- `OmegaProvider` unificado entre ambos solvers (coherencia)
- Docstrings claros sobre la variable primaria (u_e, NO desplazamiento total)

**Debilidades:**
- ⚠️ **MODERADO**: La inicialización no valida flags incompatibles
  - `include_geometric_stiffness`, `include_spin_softening`, etc. NO aplican
  - Si el usuario los activa, el solver los ignora silenciosamente
  - Mejor: levantar `ValueError` si están activos

**Recomendación:**  
Validación en `_init_rotor_config()`:
```python
incompatible_flags = [
    "include_geometric_stiffness",
    "include_spin_softening",
    "include_centrifugal",
    "include_coriolis",
    "include_euler",
]
for flag in incompatible_flags:
    if rotor_cfg.get(flag, False):
        raise ValueError(
            f"{flag}=True is incompatible with inertial solver. "
            f"These terms are replaced by -M·a_ref in the inertial formulation."
        )
```

---

#### **Fase 3: Assembly en Geometría Rotada** (5 tareas)
✅ **EVALUACIÓN: FUNCIONAL PERO CON RIESGOS**

**Fortalezas:**
- Separación clara: rotación → rebuild assembler → reassemble K(θ)
- M reutilizada (invariante bajo rotación rígida)
- C(θ) = η_k·K(θ) + η_m·M coherente con Rayleigh estándar

**Debilidades CRÍTICAS:**

1. **⚠️⚠️ CRÍTICO: Reconstrucción completa del assembler por ventana**
   - Costo computacional prohibitivo para rotores con >10k DOFs
   - `_rebuild_assembler_with_rotated_geometry()` reinstancia todo el backend Rust
   - En la práctica, esto puede ser 50-80% del tiempo total por ventana

2. **⚠️ CRÍTICO: Orientación de DOFs rotacionales en shells**
   - La rotación rígida cambia el frame local de los elementos MITC
   - Los DOFs rotacionales (θx, θy) se interpretan respecto al frame local
   - **PROBLEMA NO RESUELTO**: ¿Se transforman automáticamente o requieren corrección?
   - Esto puede causar errores silenciosos en problemas de flexión fuera del plano

3. **⚠️ MODERADO: No hay validación de preservación de propiedades**
   - No se verifica que `det(K(θ)) = det(K(0))` (para rigidez isotrópica)
   - No se verifica que `tr(M_rotated) = tr(M_ref)` (masa total constante)

**Recomendación URGENTE:**  
Agregar API Rust para actualizar coords sin reconstruir topology:
```rust
// En assembler.rs
impl PyMeshAssembler {
    pub fn update_node_coordinates(&mut self, coords: Vec<[f64; 3]>) {
        self.mesh.update_node_positions(coords);
        // Re-compute element geometry but keep connectivity
    }
}
```

**Test de regresión necesario:**
```python
def test_shell_drilling_dof_under_rotation():
    """Verify drilling DOF orientation transforms correctly."""
    solver = LinearDynamicFSIRotorInertial(...)
    theta = np.pi / 4  # 45 degrees
    solver._rotate_structural_geometry_internal(theta)
    K_rotated = solver.domain.assemble_stiffness_matrix()
    # For isotropic material, eigenvalues should be rotation-invariant
    eigs_0 = compute_eigenvalues(K_ref)
    eigs_theta = compute_eigenvalues(K_rotated)
    assert np.allclose(sorted(eigs_0), sorted(eigs_theta), rtol=1e-10)
```

---

#### **Fase 4: Newmark Integration** (4 tareas)
✅ **EVALUACIÓN: CORRECTO PERO INCOMPLETO**

**Fortalezas:**
- Newmark-β estándar (β=0.25, γ=0.5) es incondicionalmente estable
- Checkpoint/rollback implementado para coupling implícito
- Separación clara entre `K_eff` assembly y RHS

**Debilidades:**

1. **⚠️ MODERADO: Coeficientes Newmark hardcoded**
   - β y γ no configurables desde YAML
   - No permite experimentar con otros esquemas (e.g., HHT-α)

2. **⚠️ CRÍTICO: Falta validación de conservación de energía**
   - En rotación pura sin cargas, energía mecánica debe conservarse
   - No hay tests de regresión que verifiquen esto
   - Pérdida numérica de energía puede causar drift a largo plazo

3. **⚠️ MODERADO: Checkpoint no guarda historial completo**
   - Solo guarda (u, v, a) del paso actual
   - Si se implementa multi-step method en el futuro, esto será insuficiente

**Recomendación:**  
Test de conservación de energía:
```python
def test_energy_conservation_pure_rotation():
    """Rigid rotation should conserve mechanical energy."""
    solver = LinearDynamicFSIRotorInertial(
        omega=10.0,  # rad/s constant
        include_gravity=False,  # no external forces
    )
    E_initial = solver.compute_total_energy(u_0, v_0)
    u, v, a = solver.solve()
    E_final = solver.compute_total_energy(u, v)
    # Energy should remain constant within numerical tolerance
    assert abs(E_final - E_initial) / E_initial < 1e-6
```

---

#### **Fase 5: Contrato preCICE** (5 tareas)
✅ **EVALUACIÓN: ARQUITECTURA CORRECTA**

**Fortalezas:**
- Interface fijo en preCICE (SolidMesh NO rota) es conceptualmente limpio
- Escritura de `u_e` (elástico) vs `u_total` está bien documentada
- GlobalSolidMesh para ω reutiliza infraestructura existente

**Debilidades CRÍTICAS:**

1. **⚠️⚠️ CRÍTICO: Incompatibilidad con participantes que esperan desplazamiento total**
   - OpenFOAM con `movingWallVelocity` espera `u_total = x - X_0`
   - El solver escribe `u_e = x - x̂`, donde `x̂ = R(θ)·X_0`
   - **RESULTADO**: El fluido verá una deformación ficticia igual a la rotación rígida
   - Esto ROMPE el acoplamiento silenciosamente

2. **⚠️ CRÍTICO: No hay validación del contrato en runtime**
   - El solver NO verifica que el participante fluido esté usando GlobalSolidMesh
   - Si el CFD NO lee ω, las geometrías estarán desalineadas

**Recomendación URGENTE:**  
Agregar modo de compatibilidad:
```yaml
solver:
  rotor:
    precice_displacement_mode: "elastic"  # o "total"
```

Y validar en runtime:
```python
if self._precice_displacement_mode == "total":
    # Write u_total = u_e + (R(θ)·X_0 - X_0)
    u_rigid = self._compute_rigid_displacement(theta, X_0)
    u_fsi = u_e + u_rigid
else:
    # Default: elastic only
    u_fsi = u_e
```

**Test de integración necesario:**
```python
def test_precice_contract_with_mock_fluid():
    """Verify FSI contract: fluid sees correct displacement field."""
    # Mock fluid participant that expects elastic displacement
    fluid = MockFluidParticipant(expects_elastic_displacement=True)
    solver = LinearDynamicFSIRotorInertial(...)
    solver.solve()
    # Verify fluid received u_e, not u_total
    assert fluid.received_elastic_displacement_only()
```

---

#### **Fase 6: Dinámica de Omega** (4 tareas)
✅ **EVALUACIÓN: CORRECTO**

**Fortalezas:**
- Reutiliza `OmegaProvider` sin cambios (coherencia)
- Torque driving correcto: `τ = τ_aero + τ_gravity + τ_shaft`
- Actualización solo después de ventana convergida (estabilidad)

**Debilidades:**
- ⚠️ **MENOR**: Logging de torques podría ser más detallado
- No hay salida CSV automática de historia temporal de torques

---

#### **Fase 7-10: Validación y Optimización** (NO IMPLEMENTADAS)

**Estado**: Solo documentadas, sin código.

**Prioridad URGENTE:**

1. **Phase 7.1**: Test de rotación rígida pura  
   Sin este test, NO sabemos si la formulación es correcta.

2. **Phase 7.4**: Test de contrato preCICE  
   Sin este test, el acoplamiento con CFD puede fallar silenciosamente.

3. **Phase 8.1**: Benchmark A/B con corotacional  
   Sin esto, no hay evidencia de que el solver sea útil.

---

### 1.2 Calidad del Código Implementado

**Métricas:**
- Líneas de código: 2333 (vs 1918 del corotacional)
- Métodos públicos: 35
- Cobertura de tests: **0%** (NO HAY TESTS)
- Documentación: Buena (docstrings completos)

**Deuda técnica acumulada:**

| Categoría | Severidad | Items |
|-----------|-----------|-------|
| Tests faltantes | **CRÍTICA** | 0 tests unitarios, 0 tests de integración |
| Validación física | **CRÍTICA** | No hay verificación de conservación de energía |
| Compatibilidad preCICE | **CRÍTICA** | Incompatibilidad silenciosa con OpenFOAM estándar |
| Performance | **ALTA** | Reconstrucción de assembler por ventana (50-80% overhead) |
| Configuración | **MEDIA** | Flags incompatibles no validados |
| Robustez | **MEDIA** | No hay manejo de casos degenerados (ω→0, θ→2π) |

---

## 2. Análisis Crítico de la Física del Problema

### 2.1 Formulación Matemática

#### **Ecuación en Marco Inercial:**
```
M·ü_e + C(θ)·u̇_e + K(θ)·u_e = F_aero + F_g - M·a_ref
```

donde:
- `a_ref = α × r + ω × (ω × r)` (aceleración de referencia rígida)

#### **Comparación con Formulación Corotacional:**

| Término | Inertial | Corotacional |
|---------|----------|--------------|
| **Variable primaria** | `u_e` (global) | `u_local` (frame rotante) |
| **Marco de referencia** | Inertial fijo | Rotante con rotor |
| **Rigidez** | `K(θ)` variable | `K` constante |
| **Masa** | `M` constante | `M` constante |
| **Carga inercial** | `-M·a_ref` (RHS) | `K_SP + K_G + G_cor` (LHS) |
| **Transformación** | Ninguna | `u_global = R(θ)·u_local` |

---

### 2.2 Validez Física de la Formulación

#### **✅ CORRECTO:**

1. **Separación limpia entre cinemática rígida y deformación elástica**
   - `x(t) = x̂(t) + u_e(t)` donde `x̂(t) = R(θ(t))·X_0`
   - Esto es físicamente equivalente a la descripción Lagrangiana estándar

2. **Carga de referencia -M·a_ref**
   - Físicamente correcta para marco inercial
   - Equivale a "restar la inercia del movimiento rígido"
   - Es consistente con el principio de d'Alembert

3. **K(θ) variable**
   - Correcto para materiales ortotrópicos o compuestos
   - Para isotropía, `K(θ)` debería ser invariante (pero eso no está testeado)

#### **⚠️ LIMITACIONES FÍSICAS IMPORTANTES:**

1. **NO incluye prestress centrífugo (K_G)**
   - En el solver corotacional, `K_G` modela el stiffening por tensión centrífuga
   - En el solver inertial, esta rigidización NO está presente
   - **CONSECUENCIA**: Frecuencias naturales MENORES que corotacional
   - **IMPACTO**: Subestimación de rigidez efectiva en alta velocidad

2. **NO incluye spin-softening (K_SP)**
   - `K_SP = -ω²·M·(I - n̂⊗n̂)` reduce rigidez transversal en frame rotante
   - En frame inercial, este efecto NO aparece naturalmente
   - **CONSECUENCIA**: Comportamiento dinámico DIFERENTE del corotacional
   - **IMPACTO**: No es simplemente "otra formulación equivalente"

3. **Asume deformaciones pequeñas sobre geometría rotada**
   - Hipótesis: `||u_e|| << ||r||` (deformación ≪ radio)
   - Si la pala se deflecta mucho, la rotación geométrica de `u_e` no está capturada
   - **LÍMITE**: Válido solo para deflexiones pequeñas (~5% del radio)

---

### 2.3 Análisis de Consistencia Termodinámica

#### **Conservación de Energía Mecánica (sin amortiguamiento):**

Energía total: `E = T + U + V_g`

- Cinética: `T = ½·v_e^T·M·v_e + v_e^T·M·v_rigid` (término cruzado!)
- Potencial elástica: `U = ½·u_e^T·K(θ)·u_e`
- Potencial gravitacional: `V_g = -F_g^T·u_e`

**⚠️ PROBLEMA NO ANALIZADO:**  
El término cruzado `v_e^T·M·v_rigid` NO está en el balance de energía discreto.

En continuo:
```
dE/dt = -v_e^T·C·v_e  (disipación por Rayleigh)
```

En discreto (Newmark):
```
E_{n+1} - E_n = ?
```

**NO HAY DEMOSTRACIÓN** de que Newmark conserve esta energía modificada.

**Recomendación:**  
Implementar `compute_total_energy()` y verificar drift:
```python
def compute_total_energy(self, u_e, v_e, theta, omega):
    # Kinetic
    v_rigid = self._compute_rigid_velocity(omega, coords)
    v_total = v_e + v_rigid
    T = 0.5 * v_total.dot(self.M.dot(v_total))
    
    # Potential
    U = 0.5 * u_e.dot(self.K.dot(u_e))
    V_g = -self.F_gravity.dot(u_e)
    
    return T + U + V_g
```

---

## 3. Análisis de la Solución Numérica

### 3.1 Esquema de Integración Temporal

**Newmark-β con β=0.25, γ=0.5 (trapezoidal rule):**

**✅ FORTALEZAS:**
- Incondicionalmente estable
- 2do orden de precisión temporal
- Sin disipación numérica (conservativo para oscilaciones)

**⚠️ DEBILIDADES:**

1. **No disipa oscilaciones espurias de alta frecuencia**
   - En problemas con mallas gruesas, puede amplificar modos espurios
   - HHT-α (α=-0.05) sería mejor para robustez

2. **Predictor implícito requiere iteración**
   - El solver actual usa KSP directo (LU)
   - Para >100k DOFs, esto es inviable (O(N³))

3. **Actualización de v y a no usa forma incremental**
   - Pérdida de precisión numérica en `v_new = a1·(u_new - u_prev) + ...`
   - Mejor: forma incremental `Δu = u_new - u_prev`, luego `v_new = v_prev + f(Δu)`

---

### 3.2 Acoplamiento preCICE (IQN-ILS)

**Estrategia:**
- Implicit coupling con checkpoint/rollback
- Matrices K(θ), C(θ) fijas durante sub-iterations
- Solo u_e se actualiza en cada sub-iteration

**✅ CORRECTO:** K(θ) congelada en ventana es consistente.

**⚠️ RIESGO:** Si preCICE hace >10 sub-iterations:
- Costo de factorización LU repetida
- No hay preconditioner reutilizado
- Mejor: factorizar K_eff UNA VEZ, reutilizar en sub-iterations

**Recomendación:**  
```python
# En solve(), fuera del loop de sub-iterations:
ksp.setOperators(K_eff)
ksp.setUp()  # Factoriza UNA VEZ

# Dentro del loop:
ksp.solve(F_eff, u_e_new)  # Reutiliza factorización
```

---

### 3.3 Reensamblado de K(θ) por Ventana

**Frecuencia:** Una vez por ventana FSI convergida.

**Costo estimado** (mesh con N nodos):
1. Rotar coords: O(N) ✅
2. Reconstruir assembler Rust: O(N·E) ⚠️ (E = elementos)
3. Ensamblar K(θ): O(N·E·d³) ⚠️ (d = DOFs/nodo)
4. Factorizar K_eff: O(N³/2) ⚠️⚠️ (para LU denso)

**PROBLEMA CRÍTICO:**  
Para N=10,000 DOFs (rotor típico), el paso 4 domina:
- LU denso: ~10 minutos por ventana
- Con 1000 ventanas → 166 horas (1 semana!)

**Solución URGENTE:**  
Usar factorización sparse (MUMPS, SuperLU_dist):
```python
ksp.getPC().setType("lu")
ksp.getPC().setFactorSolverType("mumps")  # Parallel sparse LU
```

Esto reduce O(N³) → O(N^1.5) típicamente.

---

### 3.4 Estabilidad Numérica a Largo Plazo

**Escenario:** Simulación de 100 revoluciones (θ → 200π).

**Riesgos identificados:**

1. **Pérdida de ortogonalidad en R(θ)**
   - Rodrigues' formula acumula error de redondeo
   - Después de 100 revoluciones, `R^T·R ≠ I` puede tener error O(10⁻¹⁰)
   - **IMPACTO**: Deriva geométrica, pérdida de masa/rigidez

   **Solución:**  
   Reortogonalizar cada 10 revoluciones:
   ```python
   if theta > 20*np.pi:  # Cada 10 revoluciones
       R = self._coord_transforms.rotation_matrix(theta)
       Q, _ = np.linalg.qr(R)  # Gram-Schmidt
       self._coord_transforms._update_matrix(Q)
       theta_normalized = 0  # Reset ángulo
   ```

2. **Drift de energía mecánica**
   - Sin disipación explícita, Newmark puede acumular error
   - **SÍNTOMA**: Amplitud de oscilación crece lentamente
   - **DETECCIÓN**: Monitorear `||u_e||_max` cada 10 ventanas

3. **Inestabilidad en transiciones de signo de ω**
   - Si ω cambia de signo (frenado), θ puede volverse discontinuo
   - **SOLUCIÓN**: Usar `atan2` para θ, no integración directa

---

## 4. Problemas Críticos No Resueltos

### 4.1 🔴 CRÍTICO: Incompatibilidad con OpenFOAM Estándar

**Problema:**  
OpenFOAM espera `u_total = x - X_0` en SolidMesh.  
El solver escribe `u_e = x - x̂ = x - R(θ)·X_0`.

**Consecuencia:**  
- El fluido verá deformación ficticia = `R(θ)·X_0 - X_0`
- Resultados aerodinámicos INCORRECTOS
- **SEVERIDAD**: Bloqueante para uso en producción

**Solución:**  
Implementar modo de compatibilidad (ver sección 1.1, Fase 5).

---

### 4.2 🔴 CRÍTICO: Falta de Tests de Validación

**Problema:**  
NO HAY TESTS que verifiquen:
- Conservación de masa bajo rotación
- Conservación de energía (caso sin amortiguamiento)
- Correctitud de K(θ) vs K(0) para isotropía
- Contrato preCICE (elastic vs total displacement)

**Consecuencia:**  
- Bugs silenciosos pueden pasar desapercibidos
- No hay forma de detectar regresiones en refactors

**Solución:**  
Implementar suite mínima (ver recomendaciones en secciones previas).

---

### 4.3 🔴 CRÍTICO: Performance Inviable para Mallas Grandes

**Problema:**  
Reconstruir assembler + LU denso → O(N³) por ventana.

**Impacto medido** (estimación):
| N (DOFs) | Tiempo/ventana | 1000 ventanas | Viable? |
|----------|----------------|---------------|---------|
| 1,000 | 0.5 s | 8 min | ✅ Sí |
| 10,000 | 60 s | 16 horas | ⚠️ Marginal |
| 100,000 | 3000 s | 35 días | ❌ NO |

**Solución:**  
1. Sparse LU (MUMPS) → reduce a O(N^1.5)
2. API Rust `update_coords()` → elimina reconstrucción
3. Iterative solver (CG + AMG precond) → O(N log N)

---

### 4.4 ⚠️ MODERADO: Omisión de K_G Cambia Física

**Problema:**  
Sin prestress centrífugo, frecuencias naturales serán MENORES.

**Impacto:**  
- Comparación con corotacional NO es "formulación equivalente"
- Es una **aproximación de rigidez reducida**
- Válido solo para ω bajas (<5 rad/s típicamente)

**Solución:**  
Documentar limitación explícitamente:
```yaml
# YAML documentation
solver:
  type: LinearDynamicFSIRotorInertial
  # NOTE: This solver omits centrifugal prestress (K_G).
  # Natural frequencies will be LOWER than corotational solver.
  # Valid for low ω (<5 rad/s) or when K_G contribution is negligible.
```

---

## 5. Recomendaciones Prioritarias

### 5.1 🔴 URGENTE (Antes de usar en producción)

1. **Implementar tests de validación física**
   - Rotación rígida pura (||u_e|| < tol)
   - Conservación de energía
   - Correctitud de K(θ) para isotropía

2. **Resolver incompatibilidad preCICE**
   - Agregar modo `displacement_mode: total`
   - Validar con OpenFOAM real

3. **Optimizar performance**
   - Sparse LU (MUMPS)
   - API Rust update_coords()
   - Medir speedup vs corotacional

### 5.2 ⚠️ IMPORTANTE (Corto plazo)

4. **Validar orientación de DOFs rotacionales**
   - Test con shell fuera del plano
   - Verificar transformación automática

5. **Agregar monitoreo de estabilidad**
   - Drift de energía
   - Ortogonalidad de R(θ)
   - Advertencias automáticas

### 5.3 ✅ DESEABLE (Mediano plazo)

6. **Generalizar a masa consistente**
   - Implementar -M·a_ref para M no diagonal

7. **Agregar K_G opcional**
   - Permitir activar prestress centrífugo
   - Comparar con corotacional en igualdad de condiciones

8. **HHT-α damping**
   - Opción de integración más robusta

---

## 6. Conclusiones

### 6.1 Evaluación Global

| Aspecto | Calificación | Comentario |
|---------|--------------|------------|
| **Plan de implementación** | 7/10 | Estructura sólida, pero falta validación |
| **Formulación física** | 8/10 | Correcta pero incompleta (sin K_G) |
| **Implementación numérica** | 6/10 | Funcional pero no optimizada |
| **Testing** | 1/10 | Prácticamente ausente |
| **Documentación** | 9/10 | Excelente |
| **Production readiness** | 3/10 | No listo sin tests y optimización |

### 6.2 Veredicto Final

**✅ EL SOLVER ESTÁ BIEN DISEÑADO** desde el punto de vista arquitectónico.

**⚠️ NO ESTÁ LISTO PARA PRODUCCIÓN** por:
1. Incompatibilidad no resuelta con OpenFOAM estándar
2. Performance inviable para mallas >10k DOFs
3. Ausencia total de tests de validación

**🎯 PLAN DE ACCIÓN RECOMENDADO:**

**Semana 1-2**: Tests + validación física  
**Semana 3-4**: Optimización de performance  
**Semana 5-6**: Integración con OpenFOAM + benchmarks  
**Semana 7**: Decisión de product (promover o mantener experimental)

**Esfuerzo estimado**: 6-7 semanas de 1 desarrollador senior.

---

## Referencias

1. Design document: `docs/rotor_inertial_solver_design.md`
2. Implementation: `src/aeroelast/solvers/fsi/rotor_inertial.py`
3. Corotational solver: `src/aeroelast/solvers/fsi/rotor.py`
4. Task breakdown: `docs/rotor_inertial_solver_tasks.md`

**Análisis realizado**: Mayo 6, 2026  
**Revisor recomendado**: Especialista en FSI + validación numérica
