# Proposal: element-corotational-shell

## Intent

Los elementos MITC3/MITC4 operan en Total Lagrangian con frame fijo. Para grandes rotaciones
flap/edgewise (θ > 15°) en palas de rotor, el error en la deformación es O(θ²):
a 20° el error alcanza ~12%, lo que degrada la fidelidad del modelo FSI sin aumentar
la complejidad del problema no-lineal.

La formulación **corotacional de elemento** (Crisfield/Felippa) captura grandes rotaciones
con pequeñas strains locales, sin requerir Newton-Raphson completo. Esto la hace
compatible con el solver Newmark lineal-por-paso ya existente.

**Por qué ahora**: el solver `LinearDynamicFSIRotorCorotationalSolver` ya aplica
un frame corotacional a nivel de estructura (rotación del rotor), pero el kernel
de elemento sigue siendo TL. El error O(θ²) acumulado en cada paso de tiempo
contamina los resultados flap a velocidades nominales (ω ≥ 5 rpm en palas > 40 m).

## Approach: Formulación Corotacional de Elemento

### Marco teórico

Para cada elemento, la rigidez tangente global corotacional es:

```
K_T^coro = Tᵀ · K_L · T + K_σ^coro
```

donde:
- **T (24×24 o 18×18)**: matriz de transformación de coordenadas globales al frame
  corotacional deformado. Depende de u_n (configuración actual).
- **K_L**: rigidez local evaluada en el frame corotacional ≡ compute_ke_local() existente
  (los strains locales son pequeños → K_L es la rigidez lineal del elemento referenciada
  al frame deformado).
- **K_σ^coro**: contribución geométrica del spin del frame T(u) al diferencial de
  fuerzas internas.

### Construcción de T — `build_t24_deformed(frame_def)`

El frame corotacional deformado {ê₁_def, ê₂_def, ê₃_def} se extrae del plano medio
deformado:
1. Coordenadas deformadas: x_def = x₀ + u_{trans} (solo traslaciones nodales).
2. Normal deformada: ê₃_def = normalize(g_r_def × g_s_def) (en el centroide del elemento).
3. ê₁_def = normalize(g_r_def − (g_r_def·ê₃_def)·ê₃_def).
4. ê₂_def = ê₃_def × ê₁_def.

La matriz T24_deformed transforma **tanto traslaciones como rotaciones**:
```
T24_def[3i:3i+3, 3i:3i+3] = R_frame  ∀ i ∈ {0..7} (6 tripletes)
```
donde R_frame = [ê₁_def | ê₂_def | ê₃_def]ᵀ es la misma rotación para traslaciones
y rotaciones (Mindlin shell, rotaciones no-independientes en el frame corotacional).

### Polar decomposition via SVD — `polar_decomposition_svd(F)`

F = R · S  donde R ortogonal, S simétrico positivo definida.
Implementación: F = U Σ Vᵀ (SVD), luego R = U Vᵀ, S = V Σ Vᵀ.
Ventajas sobre Denman-Beavers: exacto, O(1) iteraciones, maneja det(F) ≈ −1 sin
divergencia, tolerancia numérica 1e-14.

**Bug actual**: el código usa iteración DB con fórmula incorrecta
`U_new = 0.5*(U + C·U⁻¹)` que NO converge a sqrt(C) — converge solo si se usa
`0.5*(U + U^{-T})` sobre F directamente. El test con tol=0.05 lo tapaba.

### K_σ^coro (rigidez geométrica corotacional)

En la formulación Crisfield (Vol 2, Cap 17), K_σ en el marco corotacional incluye:
1. K_σ_material = Bᵀ σ B (igual al caso TL, calculado en frame deformado).
2. K_spin = ∂(Tᵀ f_int)/∂u|_{spin} = contribución del cambio de T con u.

Para el solver Newmark lineal (sin NR), **K_spin ≈ 0** es justificado:
- La linearización en un paso es consistente a primer orden en Δu.
- K_spin es O(||f_int||·||ΔT/Δu||) ≈ O(σ·θ̇·Δt) — pequeño para Δt típicos FSI.
- Incluir K_spin sin NR puede desestabilizar: agrega términos asimétricos sin
  iteración de equilibrio que los corrija.
- **Decisión**: K_σ^coro ≡ K_σ_TL(frame_def) = Bᵀ σ B evaluado en coords deformadas.
  Error: O(θ²·σ·Δt) — aceptable para el régimen objetivo (strains < 2%).

## Scope

### In
1. Fix `polar_decomposition` → SVD (nalgebra `svd()`)
2. `build_t24_deformed(frame_def)` en mitc4.rs — nuevo método
3. `compute_kt_corotational(u_elem, coords_elem)` en MITC4
4. Primitivas corotacionales en MITC3 + `compute_kt_corotational()`
5. `assemble_kt_corotational(u: &[f64])` en MeshAssembler
6. Flag `use_corotational_kt: bool` + `kt_update_freq: usize` en rotor_fsi.rs y rotor_inertial.rs
7. Tests: rotación 30°/45°/90°, patch test, regresión TL, comparación lineal

### Out
- Newton-Raphson completo
- Updated Lagrangian
- Elementos sólidos
- Director independiente por nodo
- Elementos QUAD planos (plane stress)

## Risks y mitigaciones

| Riesgo | Probabilidad | Mitigación |
|--------|-------------|------------|
| Bug en SVD sign convention (R no SO(3)) | Media | Test det(R)=+1 para todos los ángulos |
| build_t24_deformed inconsistente con traslaciones vs rotaciones | Alta | Spec SC-012: invariancia bajo rotación rígida |
| kt_update_freq > 1 introduce deriva temporal | Media | Default=1; advertencia en log si freq>5 |
| Regresión en casos existentes con use_corotational_kt=false | Baja | SC-031: test de regresión explícito |
| Inestabilidad Newmark con K_σ geométrico en el frame deformado | Baja | Filtro tensile_part existente aplica |

## Success Criteria

1. `polar_decomposition_svd`: error en R < 1e-12 para θ = 0°, 30°, 45°, 90°, 180°
2. `compute_kt_corotational`: rigidez invariante bajo rotación rígida pura (K_T_coro ≡ K_T_TL)
3. Patch test corotacional: solución exacta para deformación uniforme
4. Regresión: con `use_corotational_kt=false` todos los tests existentes pasan sin cambio
5. Comparación TL: para θ < 5° la diferencia entre K_T_coro y K_T_TL < 0.1%
6. Error corotacional: para θ = 20° la diferencia con TL > 5% (confirma el beneficio)
