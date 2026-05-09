# Design: LinearDynamicFSIRotor Inertial Solver

## Technical Approach

El objetivo es introducir un segundo solver de rotor que resuelva la respuesta estructural en el marco inercial, manteniendo la rotacion rigid-body solo como cinematica interna del solver estructural y preservando la interfaz FSI con malla estatica en preCICE.

La estrategia propuesta es:

1. Renombrar el solver actual a `LinearDynamicFSIRotorCorotationalSolver` y mantener compatibilidad temporal con el nombre publico actual.
2. Crear un nuevo `LinearDynamicFSIRotorInertialSolver` que use una formulacion con referencia rigidamente rotada y desplazamiento elastico expresado en coordenadas globales inerciales.
3. Mantener `SolidMesh` fija en preCICE y enviar unicamente desplazamiento elastico global, mientras la rotacion rigid-body del rotor sigue viviendo del lado del solver estructural y del canal `GlobalSolidMesh` para la velocidad angular.
4. Reensamblar la rigidez estructural solo cuando cambia la configuracion rigidamente rotada de la ventana temporal. Dentro de una misma ventana FSI, la cinematica angular se congela y la rigidez se reutiliza en todas las subiteraciones de preCICE.

La variable primaria del nuevo solver no es el desplazamiento total absoluto del nodo, sino el desplazamiento elastico respecto de una configuracion de referencia que rota rigidamente:

$$
\mathbf{x}(t) = \hat{\mathbf{x}}(t) + \mathbf{u}_e(t),
\qquad
\hat{\mathbf{x}}(t) = \mathbf{R}(\theta(t))\,\mathbf{X}_0
$$

donde $\hat{\mathbf{x}}(t)$ representa la geometria rigida del rotor y $\mathbf{u}_e$ la deformacion elastica pequena superpuesta a esa configuracion.

En esta formulacion, si $\mathbf{u}_e$ se expresa en coordenadas globales inerciales, no aparecen fuerzas ficticias de marco rotante. La ecuacion estructural propuesta para la version inicial es:

$$
[\mathbf{M}]\{\ddot{\mathbf{u}}_e\}
+ [\mathbf{C}(\theta)]\{\dot{\mathbf{u}}_e\}
+ [\mathbf{K}(\theta)]\{\mathbf{u}_e\}
= \{\mathbf{F}_{aero}^{glob}\} + \{\mathbf{F}_g^{glob}\} - [\mathbf{M}]\{\ddot{\hat{\mathbf{x}}}\}
$$

con la aceleracion de referencia rigid-body

$$
\ddot{\hat{\mathbf{x}}}_i =
\boldsymbol{\alpha} \times \mathbf{r}_i
+ \boldsymbol{\omega} \times (\boldsymbol{\omega} \times \mathbf{r}_i)
$$

donde $\mathbf{r}_i = \hat{\mathbf{x}}_i - \mathbf{c}$ es la posicion del nodo respecto al centro de rotacion.

Esta carga inercial equivale, en el marco inercial, al efecto de imponer que la referencia del rotor siga una rotacion rigid-body prescrita o calculada. En la version inicial no se incluyen $[\mathbf{K}_G]$, $[\mathbf{K}_{SP}]$ ni $[\mathbf{G}_{cor}]$; esos terminos pertenecen a la formulacion corotacional actual y no deben mezclarse en una primera implementacion del solver inercial.

## Architecture Decisions

| Decision | Choice | Alternatives | Rationale |
|----------|--------|--------------|-----------|
| Nomenclatura publica | Introducir `LinearDynamicFSIRotorCorotational` y `LinearDynamicFSIRotorInertial`; mantener `LinearDynamicFSIRotor` apuntando temporalmente al corotacional | Renombrado inmediato y swap del alias por defecto | Reduce riesgo de romper casos existentes y permite comparacion A/B controlada |
| Variable primaria | Resolver desplazamiento elastico $\mathbf{u}_e$ sobre referencia rigidamente rotada | Resolver desplazamiento total absoluto con BC temporales en la raiz | La referencia movil evita condiciones de borde temporales y separa claramente rotacion rigida de deformacion elastica |
| Marco de expresion | Expresar $\mathbf{u}_e$ en coordenadas globales inerciales | Expresar $\mathbf{u}_e$ en coordenadas locales rotantes | En coordenadas globales no aparecen terminos ficticios de Coriolis/Euler en la ecuacion elastica |
| Contrato FSI | Mantener `SolidMesh` estatica y escribir solo desplazamiento elastico global | Rotar la malla de preCICE o escribir desplazamiento total | Evita doble rotacion del lado fluido y sigue la recomendacion de preCICE para mallas geometricamente fijas |
| Reensamblado estructural | Reensamblar una vez por ventana temporal convergida, no por subiteracion FSI | Reensamblar en cada subiteracion de preCICE | Dentro de la ventana se congela $\bar{\omega}$ y $\theta_{target}$, por lo que $[\mathbf{K}(\theta)]$ no cambia |
| Matriz de masa | Reutilizar $[\mathbf{M}]$ fija | Reensamblar $[\mathbf{M}]$ junto con $[\mathbf{K}]$ | La masa consistente o lumped es invariante bajo rotacion rigid-body |
| Amortiguamiento | Reconstruir $[\mathbf{C}(\theta)] = \eta_m[\mathbf{M}] + \eta_k[\mathbf{K}(\theta)]$ cuando cambie $[\mathbf{K}(\theta)]$ | Congelar $[\mathbf{C}]$ en $t=0$ | Mantiene consistencia con el modelo de Rayleigh ya usado en el repo |
| Implementacion inicial | Prototipo en Python/PETSc reutilizando el ensamblador Rust ya expuesto | Nuevo solver Rust completo desde el dia uno | Permite validar fisica y contrato FSI antes de optimizar la capa de orquestacion |
| Version 2 de performance | Agregar API de actualizacion geometrica en el ensamblador Rust o usar transformacion global $\mathbf{K}(\theta)=\mathbf{T}^T\mathbf{K}_0\mathbf{T}$ | Mantener reconstruccion total del ensamblador | La version inicial prioriza consistencia fisica; la optimizacion viene despues de validar resultados |

## Data Flow

```text
OmegaProvider -> (theta_target, omega_window, alpha_window)
              -> rotacion rigida interna de la geometria estructural
              -> ensamblado de K(theta) y C(theta)

preCICE SolidMesh (fija)
    -> lee F_aero global
    -> NO rota la geometria de interfaz

Solver estructural inercial
    -> calcula a_ref = alpha x r + omega x (omega x r)
    -> arma F_eff = F_aero + F_g - M a_ref + terminos Newmark
    -> resuelve u_e global
    -> escribe a preCICE solo u_e global
    -> escribe omega en GlobalSolidMesh

CFD / BEM
    -> usa omega para mover su propia malla o su cinematica del rotor
    -> aplica u_e sobre una geometria de interfaz que permanece fija en preCICE
```

## Solver Loop

Para cada ventana temporal $[t^n, t^{n+1}]$:

1. Guardar checkpoint estructural y de cinemativa rotacional.
2. Obtener de `OmegaProvider` la cinemativa representativa de ventana: $\theta_{target}$, $\bar{\omega}$, $\bar{\alpha}$.
3. Construir internamente la geometria rigidamente rotada $\hat{\mathbf{x}}(\theta_{target})$.
4. Reinstanciar o actualizar el ensamblador estructural con esa geometria y ensamblar $[\mathbf{K}(\theta)]$.
5. Construir $[\mathbf{C}(\theta)]$ a partir de Rayleigh y reutilizar $[\mathbf{M}]$ fija.
6. Leer de preCICE las fuerzas aerodinamicas en marco global y aplicar rampa/cap si corresponde.
7. Calcular la carga de referencia $\mathbf{F}_{ref} = -[\mathbf{M}]\ddot{\hat{\mathbf{x}}}$ y sumar gravedad global.
8. Resolver el paso Newmark para $\mathbf{u}_e^{n+1}$ en coordenadas globales.
9. Escribir en `SolidMesh` unicamente $\mathbf{u}_e^{glob}$.
10. Escribir $\bar{\omega}$ en `GlobalSolidMesh` si la configuracion lo requiere.
11. Si preCICE pide rollback, restaurar estado y reutilizar la misma geometria de ventana.
12. Si la ventana converge, actualizar $\omega$ desde el torque externo convergido y avanzar al siguiente paso.

## Interface Contract

### SolidMesh

- La geometria registrada en preCICE permanece igual a la de referencia inicial.
- El solver inercial escribe solo desplazamiento elastico global:

$$
\mathbf{u}_{fsi} = \mathbf{x} - \hat{\mathbf{x}} = \mathbf{u}_e^{glob}
$$

- El solver NO debe escribir desplazamiento total absoluto $\mathbf{x} - \mathbf{X}_0$, porque eso duplicaria la rotacion si el participante fluido ya usa `GlobalSolidMesh` para mover el rotor.

### GlobalSolidMesh

- Se mantiene sin cambios respecto del solver actual.
- Sigue transmitiendo la velocidad angular representativa de la ventana FSI convergida.

## File Changes

| File | Action | Description |
|------|--------|-------------|
| `src/aeroelast/core/config.py` | Modify | Agregar nuevos `SolverType` y documentar la separacion entre solver corotacional e inercial |
| `src/aeroelast/solvers/fsi/rotor.py` | Modify | Renombrar la implementacion actual a `LinearDynamicFSIRotorCorotationalSolver` y dejar alias de compatibilidad |
| `src/aeroelast/solvers/fsi/rotor_inertial.py` | Create | Nuevo solver inercial con referencia rigidamente rotada y desplazamiento elastico global |
| `src/aeroelast/solvers/fsi/__init__.py` | Modify | Exportar ambos solvers y ajustar imports perezosos |
| `src/aeroelast/solvers/fsi/runner.py` | Modify | Despachar `LinearDynamicFSIRotorCorotational` y `LinearDynamicFSIRotorInertial` |
| `docs/cli-reference.md` | Modify | Documentar los dos tipos de solver y el contrato de desplazamiento del solver inercial |
| `docs/teoria_formulacion_fsi_rotor.md` | Modify | Dejar explicito que el documento actual describe el solver corotacional y referenciar el nuevo diseño inercial |
| `crates/aeroelast-py/src/lib.rs` | Modify later | Fase 2: exponer helper o binding especifico para el solver inercial si se decide mover el loop a Rust |
| `crates/aeroelast-core/src/assembly/assembler.rs` | Modify later | Fase 2: API para refrescar coordenadas nodales o construir ensamblador rotado sin reprocesar topologia |

## Phase Plan

### Phase 0: Rename and Compatibility Layer

- Renombrar la clase Python actual a `LinearDynamicFSIRotorCorotationalSolver`.
- Mantener `LinearDynamicFSIRotorSolver` como alias temporal al corotacional.
- Agregar nuevos valores de `SolverType` sin cambiar el comportamiento por defecto.

Go/no-go: todos los casos existentes siguen ejecutando sin cambios de YAML.

### Phase 1: Inertial Prototype in Python

- Implementar `rotor_inertial.py` reutilizando PETSc y el flujo de `LinearDynamicFSISolver`.
- Reconstruir el ensamblador interno con geometria rigidamente rotada una vez por ventana temporal.
- Mantener el contrato FSI de malla fija y escritura de desplazamiento elastico global.

Go/no-go: el nuevo solver reproduce casos de control con fisica consistente y sin contaminar la interfaz preCICE.

### Phase 2: Performance Path

- Medir costo total por ventana: reconstruccion del ensamblador, ensamblado de $K$, factorizacion de $K_{eff}$ y solve.
- Si el costo de reconstruccion domina, agregar una API Rust de actualizacion de coordenadas o una transformacion global de matrices.

Go/no-go: la optimizacion conserva resultados dentro de tolerancia respecto al prototipo de Fase 1.

### Phase 3: Comparative Validation

- Ejecutar benchmark A/B entre corotacional e inercial.
- Comparar torque total, torque aero, 1P gravitacional, desplazamiento de punta, velocidad angular y costo computacional.
- Decidir si `LinearDynamicFSIRotor` pasa a apuntar al nuevo solver o si ambos conviven de forma permanente.

Go/no-go: la decision de producto se toma con datos de precision y costo, no solo con intuicion arquitectonica.

## Testing Strategy

| Layer | What to Test | Approach |
|-------|--------------|----------|
| Unit | Rotacion rigida de coordenadas internas | Verificar que `rotate_mesh` o helper equivalente preserve distancias, masa y topologia |
| Unit | Aceleracion de referencia $\ddot{\hat{x}}$ | Casos analiticos simples con eje fijo, $\omega$ constante y $\alpha \neq 0$ |
| Unit | Contrato FSI | Verificar que el solver inercial escribe `u_e` y no desplazamiento total |
| Integration | Rotor sin cargas externas con rotacion prescrita | La respuesta elastica debe permanecer aproximadamente nula |
| Integration | Rotor con gravedad y $\omega$ constante | Debe aparecer la modulacion 1P esperada en torque/desplazamiento |
| Integration | Caso aeroelastico comparable | Comparar desplazamiento de punta, torque aero y torque total contra el corotacional |
| Integration | `ComputedOmega` / `RampedComputedOmega` | Verificar que la actualizacion de $\omega$ usa solo torque aero + gravedad + shaft |
| Performance | Benchmark por ventana | Reportar tiempo de ensamblado, refactorizacion y solve en ambos enfoques |

## Rejected Alternatives

### Resolver desplazamiento total absoluto

Se descarta como version inicial porque obliga a imponer condiciones de borde temporales dependientes de la rotacion rigid-body en la raiz y mezcla cinematica rigida con respuesta elastica en la misma variable primaria.

### Rotar `SolidMesh` en preCICE

Se descarta porque rompe la hipotesis de malla de interfaz estatica usada hoy en el acoplamiento y complica innecesariamente el contrato con OpenFOAM/BEM.

### Portar directamente a Rust antes de validar la fisica

Se descarta para la primera iteracion porque combina dos riesgos a la vez: cambio de formulacion fisica y cambio de infraestructura de ejecucion.

## Risks and Open Questions

- **Prestress centrifugo en la version inicial:** la primera version inercial no incluira un equivalente directo de $[\mathbf{K}_G]$. Eso puede cambiar frecuencias efectivas y deflexiones respecto al corotacional; la comparacion debe explicitar este alcance.
- **Orientacion de DOFs rotacionales de shell:** hay que confirmar que al reensamblar sobre geometria rigidamente rotada, los DOFs rotacionales de los elementos MITC siguen interpretandose como incrementales respecto del frame local del elemento y no requieren una transformacion adicional.
- **Costo de reconstruccion del ensamblador:** la API actual del ensamblador Rust no expone `update_node_coords`, por lo que la Fase 1 probablemente reconstruya el ensamblador completo por ventana.
- **Contrato con el participante fluido:** el solver inercial asume que el participante fluido usa `GlobalSolidMesh` para la rotacion rigid-body y `SolidMesh` solo para deformacion elastica. Si existe un participante que espere desplazamiento total, habra que definir un modo de compatibilidad explicito.
- **Transformacion global de matrices:** la optimizacion $\mathbf{K}(\theta)=\mathbf{T}^T\mathbf{K}_0\mathbf{T}$ puede ser valida para rotacion rigid-body uniforme, pero debe verificarse cuidadosamente con shells compuestos y orientaciones materiales locales.

## Acceptance Criteria

1. El solver actual sigue disponible como corotacional sin romper YAMLs existentes.
2. El nuevo solver inercial ejecuta el loop FSI con `SolidMesh` fija y `GlobalSolidMesh` activa.
3. El nuevo solver no usa $[\mathbf{K}_G]$, $[\mathbf{K}_{SP}]$ ni $[\mathbf{G}_{cor}]$ en la fase inicial.
4. La carga inercial rigid-body $-[\mathbf{M}]\ddot{\hat{\mathbf{x}}}$ queda implementada y validada con casos simples.
5. Existe un benchmark reproducible de comparacion entre `Corotational` e `Inertial` en terminos de torque, desplazamiento y costo.

---

## Addendum: K_SP and K_G Hysteresis Triggers (inertial-spin-softening-kg)

This section documents the changes introduced by the `inertial-spin-softening-kg` SDD change, which adds per-matrix Δω-driven rebuild triggers for $[\mathbf{K}_{SP}]$ and $[\mathbf{K}_G]$ in the inertial rotor FSI solver.

### Context

K_G and K_SP were already implemented in the Rust inertial solver (`rotor_inertial.rs`) before this change. The original design doc above stated they would not be included in the initial version; they were added in a later iteration. This addendum documents the configuration surface and rebuild semantics introduced to make them controllable.

### New YAML Configuration Keys

The following keys are accepted under the `rotor:` block in the simulation YAML:

```yaml
rotor:
  # ... existing keys ...

  # Enable spin-softening stiffness K_SP (default: true).
  # K_SP = -ω² · M_lump · (I - n̂⊗n̂), added to effective stiffness K_eff.
  include_spin_softening: true

  # Enable geometric stiffness K_G from centrifugal prestress (default: false).
  # Replaces the deprecated kg0_rows/cols/vals sentinel pattern.
  include_geometric_stiffness: false

  # Hysteresis thresholds for K_SP rebuild (relative Δω² / ω²):
  #   rebuild when |Δ(ω²)| / ω² > ksp_omega_rebuild_high  (when not currently rebuilt)
  #   suppress rebuild when |Δ(ω²)| / ω² < ksp_omega_rebuild_low  (when currently rebuilt)
  ksp_omega_rebuild_high: 0.005
  ksp_omega_rebuild_low: 0.003

  # Hysteresis thresholds for K_G rebuild (same semantics as K_SP thresholds):
  kg_omega_rebuild_high: 0.005
  kg_omega_rebuild_low: 0.003
```

**Defaults summary:**

| Key | Default | Notes |
|-----|---------|-------|
| `include_spin_softening` | `true` | Active by default; set `false` to disable entirely |
| `include_geometric_stiffness` | `false` | Opt-in; replaces `kg0_rows/cols/vals` sentinel |
| `ksp_omega_rebuild_high` | `0.005` | Relative Δω² threshold to trigger K_SP rebuild |
| `ksp_omega_rebuild_low` | `0.003` | Relative Δω² threshold to suppress K_SP rebuild (hysteresis) |
| `kg_omega_rebuild_high` | `0.005` | Relative Δω² threshold to trigger K_G rebuild |
| `kg_omega_rebuild_low` | `0.003` | Relative Δω² threshold to suppress K_G rebuild (hysteresis) |

### Rebuild Trigger Semantics

Each matrix (K_SP, K_G) has an independent two-threshold hysteresis predicate:

```
omega_changed_significantly(omega_new, omega_sq_at_last_rebuild,
                            threshold_high, threshold_low, currently_rebuilt) → bool
```

The predicate returns `true` (trigger rebuild) when:
- ω = 0 and the matrix was previously built (needs to go to zero), OR
- `currently_rebuilt = false` and relative change > `threshold_high`, OR
- `currently_rebuilt = true` and relative change > `threshold_low`.

This prevents oscillation near a threshold: once a rebuild is triggered, the bar to suppress it is lower than the bar to trigger the next one.

State variables `omega_sq_at_last_ksp_rebuild` and `omega_sq_at_last_kg_rebuild` track the ω² value at the last rebuild for each matrix independently. Both are initialized to `f64::NEG_INFINITY` so the first window always triggers a rebuild.

### Single Refactorization Invariant (Spec Scenario R5)

When K(θ), K_G, and K_SP all trigger in the same coupling window, the solver calls `update_tangent_terms(k_vals, kg_vals, ksp_vals)` exactly once. This results in a single `refactorize()` call regardless of how many matrices changed. The new `NewmarkStepper::refactorize_count()` accessor exposes this counter for testing.

If only K_G changes (K_SP and K(θ) unchanged), the solver reuses the stored `stepper.k_vals()` and `stepper.ksp_vals()` slices and still calls `update_tangent_terms` once, passing the unchanged slices alongside the new K_G values.

### Checkpoint / Restore for Tracker State

Two new fields are captured in `InertialKgKspCheckpoint` and saved alongside the structural checkpoint:

- `omega_sq_at_last_ksp_rebuild`
- `omega_sq_at_last_kg_rebuild`
- `last_ksp_rebuild_step`
- `last_kg_rebuild_step`

Without this, a preCICE rollback would restore ω to its pre-window value but leave the rebuild-step trackers pointing to a future step, causing the first post-rollback rebuild to be suppressed incorrectly.

### Deprecated: kg0_rows/cols/vals Sentinel

The previous pattern of passing non-empty `kg0_rows`, `kg0_cols`, `kg0_vals` to `run_inertial_rotor_fsi_solver` as a sentinel to enable K_G is deprecated. Use `include_geometric_stiffness: true` instead. A deprecation warning is emitted at runtime if the sentinel is detected. The sentinel will be removed in a future version.