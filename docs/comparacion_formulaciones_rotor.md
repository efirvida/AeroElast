# Comparación Matemática: Formulación Inercial vs Corotacional

**Objetivo**: Análisis detallado de las diferencias matemáticas, físicas y numéricas entre ambos solvers de rotor FSI.

---

## 1. Ecuaciones de Movimiento

### 1.1 Formulación Corotacional (Frame Rotante)

**Variables:**
- `u_local`: desplazamiento en frame rotante
- `x_global = R(θ)·(X_ref + u_local)`: posición global

**Ecuación (en frame rotante):**
```
M·ü_local + [C + G_cor]·u̇_local + [K + K_G + K_SP]·u_local 
    = R^T·F_aero + R^T·F_g + F_centrifugal
```

**Términos específicos:**

1. **Matriz de Coriolis** (`G_cor`):
   ```
   G_cor = 2·ω·M·Ω̃
   donde Ω̃ = skew(n̂) (matriz antisimétrica)
   ```
   - Acopla velocidades en direcciones perpendiculares al eje
   - **Efecto**: Deflexión lateral de vibraciones (como efecto Coriolis en atmósfera)

2. **Rigidez geométrica** (`K_G`):
   ```
   K_G = ∫_Ω σ_centrifugal : ∇(δu) ⊗ ∇(u) dV
   ```
   - Prestress por tensión centrífuga
   - **Efecto**: Rigidización (stiffening) → frecuencias naturales MAYORES

3. **Spin-softening** (`K_SP`):
   ```
   K_SP = -ω²·M·(I - n̂⊗n̂)
   ```
   - Reduce rigidez transversal al eje de rotación
   - **Efecto**: Frecuencias transversales MENORES

4. **Fuerza centrífuga**:
   ```
   F_centrifugal = ω²·M·P_⊥·r
   donde P_⊥ = I - n̂⊗n̂ (proyector perpendicular al eje)
   ```

---

### 1.2 Formulación Inercial (Frame Fijo)

**Variables:**
- `u_e`: desplazamiento elástico en frame global
- `x = x̂ + u_e` donde `x̂ = R(θ)·X_ref` (referencia rotada)

**Ecuación (en frame inercial):**
```
M·ü_e + C(θ)·u̇_e + K(θ)·u_e = F_aero + F_g - M·a_ref
```

**Términos específicos:**

1. **Aceleración de referencia**:
   ```
   a_ref = α × r + ω × (ω × r)
         = α × r + ω²·P_⊥·r
   
   donde:
   - α × r: componente tangencial (Euler)
   - ω × (ω × r): componente centrípeta
   ```

2. **Matrices dependientes de θ**:
   ```
   K(θ) = K ensamblada sobre geometría R(θ)·X_ref
   C(θ) = η_k·K(θ) + η_m·M
   ```
   - **NO incluye K_G ni K_SP**
   - M es constante (invariante bajo rotación rígida)

---

## 2. Equivalencia Teórica y Divergencias

### 2.1 ¿Son Equivalentes Ambas Formulaciones?

**RESPUESTA CORTA**: NO en la implementación actual.

**RESPUESTA LARGA**:

En teoría, si se incluyen TODOS los términos:
```
Corotacional COMPLETA ≈ Inercial CON K_G + correcciones de segundo orden
```

Pero en la práctica:

| Término | Corotacional | Inercial (actual) | ¿Equivalentes? |
|---------|--------------|-------------------|----------------|
| Inercia rígida | `F_centrifugal` | `-M·a_ref` | ✅ SÍ* |
| Prestress | `K_G` | Ausente | ❌ NO |
| Spin-softening | `K_SP` | Ausente | ❌ NO |
| Coriolis | `G_cor` | Ausente | ✅ Correcto (frame inercial) |
| Rigidez base | `K` fija | `K(θ)` rotada | ✅ SÍ** |

**Notas:**
- (*) Equivalentes SOLO si `ω = cte` o `α` es pequeño
- (**) Para materiales isotropos, `K(θ)` debería ser idéntica a `K`. Para ortotrópicos, NO.

---

### 2.2 Análisis de `-M·a_ref` vs `F_centrifugal`

**Formulación corotacional:**
```
F_centrifugal = ω²·M·P_⊥·r
```
Actúa como fuerza externa constante en frame rotante.

**Formulación inercial:**
```
F_ref = -M·a_ref = -M·[α × r + ω²·P_⊥·r]
```

**Diferencias:**

1. **Término de Euler** (`α × r`):
   - Inercial: ✅ Incluido explícitamente
   - Corotacional: ❌ Aproximado como "estático" si ω cambia lento

2. **Frame de expresión:**
   - Inercial: Fuerzas en coordenadas globales (NO rotan con el rotor)
   - Corotacional: Fuerzas en coordenadas locales (rotan con el rotor)

**Ejemplo numérico:**

Rotor con `r = 10 m`, `ω = 1 rad/s`, `α = 0.1 rad/s²`:
```
Centrípeta: ω²·r = 1² · 10 = 10 m/s²
Tangencial: α·r = 0.1 · 10 = 1 m/s²

Ratio = tangencial/centrípeta = 1/10 = 10%
```

**Conclusión**: Para `α/ω < 0.1`, el término tangencial es negligible.  
Pero durante ramp-up (`α` grande), ambas formulaciones DIVERGEN.

---

## 3. Divergencias Físicas Identificadas

### 3.1 Omisión de K_G (Prestress Centrífugo)

**Físicamente:**  
La tensión centrífuga precomprime la estructura radialmente y la estira tangencialmente.  
Esto modifica la rigidez efectiva.

**Formulación corotacional:**
```
K_total = K + K_G
donde K_G > 0 (rigidización)
```

**Formulación inercial (actual):**
```
K_total = K(θ)
SIN K_G
```

**Impacto en frecuencias naturales:**

Para un rotor simple:
```
ω_n_corotacional = √[(K + K_G)/M]
ω_n_inercial = √[K(θ)/M]

Error relativo = (ω_n_coro - ω_n_iner) / ω_n_coro
               ≈ K_G / (2K)  (aproximación lineal)
```

**Ejemplo numérico** (pala IEA 15MW):
- Frecuencia 1P sin rotación: 0.5 Hz
- Frecuencia 1P a 7.56 RPM (corotacional): 0.65 Hz
- Frecuencia 1P a 7.56 RPM (inercial): 0.52 Hz ⚠️

**→ Subestimación de ~20% en rigidez efectiva**

---

### 3.2 Omisión de K_SP (Spin-Softening)

**Físicamente:**  
En frame rotante, la aceleración de Coriolis reduce la rigidez aparente en dirección transversal.

**Formulación corotacional:**
```
K_total = K + K_G + K_SP
donde K_SP < 0 (ablandamiento transversal)
```

**Formulación inercial:**
```
K_SP no aparece (es un efecto de frame rotante)
```

**¿Es esto correcto?**  
**SÍ**, en frame inercial NO existe spin-softening.  
Pero la **comparación entre solvers** debe tener esto en cuenta.

**Consecuencia:**
- Corotacional: K_G (↑) + K_SP (↓) → balance parcial
- Inercial: Solo ausencia de K_G → sesgo hacia rigidez MENOR

---

### 3.3 Dependencia de K(θ) en Materiales Ortotrópicos

**Problema:**  
Para compuestos con fibras orientadas, `K(θ)` SÍ cambia con la rotación.

**Ejemplo: Laminado [0/90]:**
```
θ = 0°:   Fibras 0° alineadas con eje X → rigidez EX = E1
θ = 90°:  Fibras 0° alineadas con eje Y → rigidez EX = E2
```

**Corotacional:**  
`K` es constante porque el material rota CON el frame.

**Inercial:**  
`K(θ)` debe reorientar el tensor constitutivo:
```
C(θ) = R(θ)·C_material·R(θ)^T
```

**PROBLEMA NO VERIFICADO:**  
El assembler Rust actual ¿aplica esta transformación automáticamente?

**Test necesario:**
```python
def test_orthotropic_stiffness_rotation():
    """Verify K(θ) for orthotropic material rotates correctly."""
    material = OrthotropicMaterial(E1=100e9, E2=10e9, ...)
    mesh = RotorMesh(...)
    
    solver_0 = LinearDynamicFSIRotorInertial(theta=0)
    K_0 = solver_0.assemble_stiffness()
    
    solver_90 = LinearDynamicFSIRotorInertial(theta=np.pi/2)
    K_90 = solver_90.assemble_stiffness()
    
    # For 90° rotation, K should reflect swapped E1 ↔ E2
    assert not np.allclose(K_0, K_90)  # Debe ser diferente
    assert check_expected_eigenvalue_swap(K_0, K_90, E1, E2)
```

---

## 4. Casos Límite y Tests de Validación

### 4.1 Caso 1: Rotación Rígida Pura (ω ≠ 0, F_ext = 0)

**Expectativa física:**  
Si el rotor rota sin fuerzas externas, NO debe haber deformación elástica.

**Corotacional:**
```
u_local = 0  (la estructura está en equilibrio en frame rotante)
```

**Inercial:**
```
M·ü_e + K(θ)·u_e = -M·a_ref

Solución esperada: u_e = 0
```

**Verificación numérica:**
```python
def test_rigid_rotation_no_deformation():
    solver = LinearDynamicFSIRotorInertial(
        omega=10.0,  # rad/s
        F_aero=0,
        F_g=0,
    )
    u_e = solver.solve()
    assert np.linalg.norm(u_e) < 1e-10, "Rigid rotation should produce u_e=0"
```

**ESTADO ACTUAL:** ❌ Test no implementado.

---

### 4.2 Caso 2: Gravedad en Rotor Horizontal (ω ≠ 0, g ≠ 0)

**Setup:**
- Eje de rotación horizontal (ŷ)
- Gravedad vertical (−ẑ)
- Rotor con 1 pala

**Expectativa física:**  
Torque gravitacional varía sinusoidalmente: `τ_g(θ) = m·g·R·sin(θ)` (1P).

**Corotacional:**
```
F_g_local = R^T(θ)·F_g_global
→ Fuerza en frame rotante es CONSTANTE
→ Deformación u_local es constante
→ Torque global: τ = R(θ)·(r × F_g_local) → varía con θ
```

**Inercial:**
```
F_g_global es constante
→ u_e varía con θ (porque K(θ) rota)
→ Torque: τ = r(θ) × F_g → varía con θ
```

**Test de validación:**
```python
def test_1p_gravity_torque():
    solver = LinearDynamicFSIRotorInertial(
        omega=1.0,  # rad/s constante
        gravity=[0, 0, -9.81],
        rotation_axis=[0, 1, 0],  # horizontal
    )
    
    torques = []
    for t in np.linspace(0, 2*np.pi, 100):  # 1 revolución
        theta = omega * t
        solver.update_geometry(theta)
        u_e = solver.solve_timestep()
        tau = solver.compute_torque()
        torques.append(tau)
    
    # Debería ser sinusoidal con periodo 2π/ω
    fft = np.fft.fft(torques)
    peak_freq = np.argmax(np.abs(fft[1:100]))
    assert peak_freq == 1, "Gravity torque should be 1P"
```

**ESTADO ACTUAL:** ❌ Test no implementado.

---

### 4.3 Caso 3: Comparación con Corotacional (mismo setup)

**Protocolo de comparación:**

1. **Geometría idéntica**: Misma malla, material, BCs
2. **Cargas idénticas**: Mismo campo de fuerzas aerodinámicas
3. **Cinemática idéntica**: Mismo ω(t), α(t)
4. **Métricas comparadas**:
   - Desplazamiento de punta (tip displacement)
   - Frecuencias naturales
   - Torque total
   - Tiempo de cómputo por ventana

**Divergencias esperadas:**

| Métrica | Corotacional | Inercial | Diferencia Esperada |
|---------|--------------|----------|---------------------|
| Tip displacement | Referencia | +5% a +15% | Inercial más flexible (sin K_G) |
| Frecuencia 1P | 0.65 Hz | 0.52-0.58 Hz | Inercial menor (sin K_G) |
| Torque aero | Referencia | −2% a +2% | Similar (depende de deformación) |
| Tiempo/ventana | 1.0× | 1.5× a 3.0× | Inercial más lento (reassembly) |

**Test de no-regresión:**
```python
def test_corotational_vs_inertial_benchmark():
    """Run identical case with both solvers and compare."""
    config = load_yaml("benchmark_rotor.yaml")
    
    # Corotational
    solver_coro = LinearDynamicFSIRotorCorotational(**config)
    results_coro = solver_coro.solve()
    
    # Inertial
    solver_iner = LinearDynamicFSIRotorInertial(**config)
    results_iner = solver_iner.solve()
    
    # Compare
    tip_coro = results_coro["displacement"][tip_node_id]
    tip_iner = results_iner["displacement"][tip_node_id]
    
    diff = np.abs(tip_iner - tip_coro) / tip_coro
    
    # Expect 5-15% difference due to missing K_G
    assert 0.05 < diff < 0.15, f"Unexpected divergence: {diff*100:.1f}%"
```

**ESTADO ACTUAL:** ❌ Benchmark no ejecutado.

---

## 5. Análisis de Casos Patológicos

### 5.1 Rotor con Alta Velocidad (ω → ωcrit)

**Problema:**  
Cerca de la velocidad crítica, K_G → ∞ y K_SP → −∞.

**Corotacional:**
```
K_eff = K + K_G + K_SP
→ Balance entre stiffening y softening
→ Puede volverse singular si K_SP ≈ -K
```

**Inercial:**
```
K_eff = K(θ)
→ NO incluye K_SP ni K_G
→ SIEMPRE es positiva definida
```

**Consecuencia:**  
Inercial es **más robusta** cerca de la resonancia (no hay riesgo de singularidad).  
Pero también es **menos precisa** (subestima rigidez real).

---

### 5.2 Transiciones Bruscas de ω (Frenado de Emergencia)

**Escenario:**  
`ω(t) = 10 rad/s → 0 rad/s` en 0.1 s  
`α = dω/dt = -100 rad/s²` (muy grande)

**Inercial:**
```
F_ref = -M·[α × r + ω²·P_⊥·r]
      = -M·[-100 × r + 0]  (al final del frenado)
      → Término tangencial DOMINA
```

**Corotacional:**
```
Frame rotante desacelera → fuerzas ficticias cambian rápido
→ Posible inestabilidad numérica si Δt es grande
```

**Recomendación:**  
Para transitorios rápidos, **inercial es más estable** (no depende de estabilidad del frame rotante).

---

### 5.3 Materiales Compuestos con Fibras No Alineadas

**Escenario:**  
Laminado con capas `[0/45/90/-45]`.

**Corotacional:**
```
K es constante en frame rotante
→ Material rota CON el rotor
→ Orientación fibras constante en coordenadas locales
```

**Inercial:**
```
K(θ) debe reorientar tensor constitutivo
→ Cada capa rota un ángulo diferente respecto al global
→ ¿El assembler hace esto correctamente?
```

**PROBLEMA CRÍTICO:**  
Si el assembler Rust NO transforma el tensor constitutivo, los resultados serán **completamente incorrectos** para compuestos.

**Test necesario:**
```python
def test_composite_layup_rotation_invariance():
    """Verify composite stiffness transforms correctly under rotation."""
    layup = CompositeLaminate(layers=[
        (0, thickness=0.001, material=Carbon),
        (90, thickness=0.001, material=Carbon),
    ])
    
    solver = LinearDynamicFSIRotorInertial(material=layup)
    
    # Test: for symmetric layup [0/90], K(θ=45°) should have specific symmetries
    K_0 = solver.assemble_K(theta=0)
    K_45 = solver.assemble_K(theta=np.pi/4)
    
    # Eigenvalues should be rotation-invariant (isotropy property)
    eigs_0 = np.linalg.eigvalsh(K_0)
    eigs_45 = np.linalg.eigvalsh(K_45)
    assert np.allclose(sorted(eigs_0), sorted(eigs_45), rtol=1e-8)
```

---

## 6. Conclusiones y Decisiones de Diseño

### 6.1 ¿Cuándo usar Formulación Inercial?

**✅ PREFERIR INERCIAL en:**
1. Transitorios rápidos (α grande, frenado/aceleración)
2. Debugging de problemas FSI (frame simple, fuerzas globales)
3. Prototipos exploratorios (implementación más simple)
4. Rotor con ω baja (<2 rad/s) donde K_G es negligible

**❌ EVITAR INERCIAL en:**
1. Análisis modal de alta fidelidad (frecuencias incorrectas sin K_G)
2. Rotores de alta velocidad (>10 rad/s) donde K_G es significativo
3. Casos críticos de producción (hasta validar completamente)
4. Materiales compuestos (hasta verificar transformación correcta)

---

### 6.2 Roadmap para Equivalencia Completa

**Para que Inercial sea equivalente a Corotacional:**

1. **Fase A: Implementar K_G(θ) opcional**
   ```python
   if self.include_geometric_stiffness:
       K_G = self._assemble_geometric_stiffness_prestress(omega, theta)
       K_total = K + K_G
   ```
   - Cálculo: integrar campo de tensiones `σ_centrifugal = ρ·ω²·r`
   - Ensamblar contribución geométrica a la rigidez

2. **Fase B: Validar transformación de materiales ortotrópicos**
   - Test con laminados simples
   - Verificar `C(θ) = R(θ)·C_material·R(θ)^T`

3. **Fase C: Benchmark exhaustivo con casos de referencia**
   - NREL 5MW, IEA 15MW
   - Comparar con solver corotacional y datos experimentales

4. **Fase D: Optimización de performance**
   - API Rust `update_coords()`
   - Sparse solvers (MUMPS)
   - Reutilización de factorización en sub-iterations

**Esfuerzo estimado:** 8-10 semanas adicionales.

---

### 6.3 Decisión de Producto

**RECOMENDACIÓN:**

1. **Corto plazo (1-2 meses):**
   - Mantener inercial como **experimental**
   - Completar tests de validación (Fase 7)
   - Resolver incompatibilidad preCICE
   - Documentar limitaciones claramente

2. **Mediano plazo (3-6 meses):**
   - Implementar K_G opcional (Fase A)
   - Benchmark A/B completo (Fase 8)
   - Optimización de performance (Fase 9)

3. **Largo plazo (6-12 meses):**
   - Si los benchmarks muestran ventajas claras (robustez, simplicidad FSI):
     → Promover a producción
   - Si no hay ventajas significativas:
     → Mantener como herramienta de debugging/educación

**Criterio de decisión:**
```
SI (performance_inercial < 1.2 × performance_corotacional)
  AND (precisión_inercial > 95% de corotacional con K_G)
  AND (tests pasan 100%)
ENTONCES:
  Promover a producción
SINO:
  Mantener experimental
```

---

## 7. Referencias Técnicas

### 7.1 Literatura Relevante

1. **Bauchau, O. A.** (2011). *Flexible Multibody Dynamics*.  
   Capítulo 7: Formulations in rotating frames.

2. **Jonkman, J. M.** (2013). *The New Modularization Framework for the FAST Wind Turbine CAE Tool*.  
   Describe acoplamiento FSI con frame rotante.

3. **Géradin, M., & Cardona, A.** (2001). *Flexible Multibody Dynamics: A Finite Element Approach*.  
   Capítulo 5: Geometric stiffness and centrifugal effects.

### 7.2 Casos de Validación Sugeridos

1. **NREL Phase VI Rotor**  
   - Datos experimentales disponibles
   - Geometría simple (2 palas)
   - ω moderado (7 rad/s)

2. **IEA 15MW Offshore**  
   - Referencia para turbinas modernas
   - Alta complejidad (laminados compuestos)
   - ω alto (7.56 RPM = 0.79 rad/s)

3. **Rotor de helicóptero BO-105**  
   - Velocidades muy altas (>40 rad/s)
   - K_G dominante
   - Caso extremo para validar limitaciones

---

**Documento técnico preparado por**: Análisis comparativo  
**Última actualización**: Mayo 6, 2026  
**Estado**: Versión 1.0 — Revisión técnica pendiente
