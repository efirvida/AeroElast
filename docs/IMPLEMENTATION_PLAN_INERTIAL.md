# Plan de Implementación: Solver Inertial — Corrección y Optimización

**Fecha**: 6 de mayo de 2026  
**Rama**: `rust_implementation`  
**Estado actual**: 14 commits, prototipo funcional, sin tests, 5 bugs críticos  

**Criterio de priorización**: Correctitud física > Robustez numérica > Performance

---

## Resumen de Issues

### Grupo A — Inconsistencia física (bloquean uso correcto)

| ID | Issue | Severidad | Esfuerzo |
|----|-------|-----------|----------|
| A1 | Contrato preCICE: escribe `u_e` en lugar de `u_total` cuando el CFD lo espera | 🔴 Bloqueante | ~50 LOC |
| A2 | K(θ) para materiales ortotrópicos: tensor constitutivo no se transforma al rotar | 🔴 Bloqueante | ~30 LOC test |
| A3 | Masa consistente no detectada: si M no es diagonal, `F_ref = -M·a_ref` es silenciosamente incorrecto | 🔴 Bloqueante | ~10 LOC |
| A4 | Flags incompatibles (`include_geometric_stiffness`, `include_spin_softening`) ignorados sin aviso | ⚠️ Importante | ~20 LOC |
| A5 | K_G ausente (prestress centrífugo): subestima rigidez efectiva ~20% para ω > 2 rad/s | ⚠️ Limitación de diseño | ~300 LOC |

### Grupo B — Robustez numérica (afectan correctitud a largo plazo)

| ID | Issue | Severidad | Esfuerzo |
|----|-------|-----------|----------|
| B1 | Checkpoint no guarda K(θ)/C(θ): después de rollback FSI las matrices son inconsistentes | 🔴 Bloqueante | ~15 LOC |
| B2 | Rampa de fuerzas no aplicada: divergencia numérica probable en t=0 | ⚠️ Importante | ~20 LOC |
| B3 | Overflow de θ para simulaciones largas (>100 revoluciones) | ⚠️ Importante | ~5 LOC |
| B4 | Pérdida de ortogonalidad en R(θ) acumulada en tiempo largo | ℹ️ Menor | ~15 LOC |

### Grupo C — Performance (no afectan correctitud)

| ID | Issue | Severidad | Esfuerzo |
|----|-------|-----------|----------|
| C1 | KSP factoriza K_eff en cada sub-iteration (debería ser 1×/ventana) | 🔴 5-10× lento | ~30 LOC |
| C2 | Assembler Rust reconstruye topología completa cada ventana | ⚠️ 2-4× lento | ~150 LOC Rust + 20 Python |
| C3 | Normas MPI locales en logging (incorrecto en paralelo) | ℹ️ Cosmético | ~10 LOC |

### Grupo D — Tests (bloquean validación)

| ID | Test | Prioridad |
|----|------|-----------|
| D1 | Rotación rígida pura: debe dar `‖u_e‖ < ε` | 🔴 Primero |
| D2 | Conservación de energía mecánica (sin amortiguamiento, sin cargas) | 🔴 |
| D3 | Contrato preCICE: escritura de `u_e` vs `u_total` | 🔴 |
| D4 | K(θ) para material isótropo: eigenvalores invariantes bajo rotación | ⚠️ |
| D5 | Torque gravitacional 1P (rotor horizontal, ω constante) | ⚠️ |
| D6 | Checkpoint/rollback: estado consistente después de rollback | ⚠️ |
| D7 | MPI: normas globales correctas en multirank | ℹ️ |

---

## Sprint 1 — Correctitud física esencial

**Objetivo**: El solver nunca produce resultados silenciosamente incorrectos.  
**Duración estimada**: 1 semana  
**Criterio de salida**: D1 + D2 + D3 pasan.

### Tarea A3: Guard para masa no-lumped

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `_assemble_system_matrices()` o `_init_rotor_config()`

```python
# Después de ensamblar M, antes de continuar:
m_type = self.M.getType()
if m_type not in ("diagonal", "seqdense"):  # PETSc types for diagonal/lumped
    raise NotImplementedError(
        "LinearDynamicFSIRotorInertialSolver requires a lumped (diagonal) mass "
        "matrix. The current assembler uses consistent mass. Either switch to "
        "lumped mass or implement full M·a_ref matrix-vector product."
    )
```

---

### Tarea A4: Validar flags incompatibles

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `_init_rotor_config()`

```python
_INERTIAL_INCOMPATIBLE_FLAGS = [
    "include_geometric_stiffness",
    "include_spin_softening",
    "include_centrifugal",
    "include_coriolis",
    "include_euler",
]

for flag in _INERTIAL_INCOMPATIBLE_FLAGS:
    if rotor_cfg.get(flag, False):
        _logger.warning(
            "%s=True has no effect in the inertial solver. "
            "These terms are replaced by -M·a_ref. "
            "To suppress this warning, set %s=False explicitly.",
            flag, flag,
        )
```

---

### Tarea A1: Modo de desplazamiento preCICE configurable

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `_write_elastic_displacement_to_precice()`

Agregar soporte para `precice_displacement_mode: "elastic" | "total"` en el YAML del rotor.

```python
# En _init_rotor_config():
self._precice_displacement_mode = rotor_cfg.get("precice_displacement_mode", "elastic")
if self._precice_displacement_mode not in ("elastic", "total"):
    raise ValueError(
        f"precice_displacement_mode must be 'elastic' or 'total', "
        f"got {self._precice_displacement_mode!r}"
    )

# En _write_elastic_displacement_to_precice():
if self._precice_displacement_mode == "total":
    # u_total = u_e + (R(θ)·X₀ - X₀)  [rigid displacement over reference]
    X_0 = self._interface_coords_reference  # shape (N, 3)
    R = self._coord_transforms.rotation_matrix(theta)
    u_rigid = ((R @ X_0.T).T - X_0).ravel()  # global DOF order
    u_fsi = u_e_flat + u_rigid
else:
    u_fsi = u_e_flat  # default: elastic only
adapter.write_data(mesh_name, data_name, u_fsi)
```

**YAML**:
```yaml
solver:
  rotor:
    precice_displacement_mode: "total"   # para OpenFOAM estándar
    # precice_displacement_mode: "elastic"  # default (inertial-aware CFD)
```

También agregar `precice_displacement_mode` a `RotorConfig` en `config.py`.

---

### Tarea B1: Checkpoint completo (incluye K y C)

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `_checkpoint_elastic_state()` / `solve()`

```python
# En solve(), al hacer checkpoint:
if adapter.requires_writing_checkpoint:
    checkpoint = self._checkpoint_elastic_state(u_e, v_e, a_e, theta, omega, alpha)
    checkpoint["K_current"] = K_current  # guardar referencia (no copy — no cambia en sub-iters)
    checkpoint["C_current"] = C_current

# En solve(), al hacer rollback:
if adapter.requires_reading_checkpoint:
    self._rollback_elastic_state(checkpoint, u_e, v_e, a_e)
    theta   = float(checkpoint["theta"])
    omega   = float(checkpoint["omega"])
    alpha   = float(checkpoint["alpha"])
    K_current = checkpoint["K_current"]   # restaurar
    C_current = checkpoint["C_current"]   # restaurar
```

---

### Tarea D1: Test de rotación rígida pura

**Archivo**: `tests/solvers/fsi/test_rotor_inertial_physics.py` (nuevo)

```python
def test_rigid_rotation_no_elastic_deformation():
    """
    A rotor spinning at constant ω with no external forces should produce
    zero elastic displacement. This validates -M·a_ref cancels the inertial
    load exactly.
    """
    solver = build_minimal_inertial_solver(omega=10.0, F_aero=None, gravity=None)
    u_e, v_e, a_e = solver.solve_one_window(dt=0.01, theta=np.pi / 6)
    assert np.linalg.norm(u_e) < 1e-10, (
        f"Rigid rotation should give ‖u_e‖ ≈ 0, got {np.linalg.norm(u_e):.3e}"
    )
```

---

### Tarea D3: Test de contrato preCICE

```python
def test_precice_writes_elastic_displacement_by_default():
    """Verify solver writes u_e (elastic) not u_total when mode='elastic'."""
    written = []
    mock_adapter = MockAdapter(on_write=written.append)

    solver = build_minimal_inertial_solver(precice_displacement_mode="elastic")
    solver.solve_with_mock_adapter(mock_adapter, theta=np.pi / 4)

    u_e_expected = solver.last_u_e_interface
    assert np.allclose(written[-1], u_e_expected, atol=1e-12), (
        "Solver must write u_e (elastic displacement) not u_total"
    )


def test_precice_writes_total_displacement_when_configured():
    """Verify solver writes u_total when mode='total'."""
    written = []
    mock_adapter = MockAdapter(on_write=written.append)
    theta = np.pi / 4

    solver = build_minimal_inertial_solver(precice_displacement_mode="total")
    solver.solve_with_mock_adapter(mock_adapter, theta=theta)

    R = solver._coord_transforms.rotation_matrix(theta)
    X_0 = solver._interface_coords_reference
    u_rigid = ((R @ X_0.T).T - X_0).ravel()
    u_total_expected = solver.last_u_e_interface + u_rigid
    assert np.allclose(written[-1], u_total_expected, atol=1e-12)
```

---

## Sprint 2 — Robustez numérica + test de energía

**Objetivo**: El solver es estable a largo plazo y tiene tests de regresión física.  
**Duración estimada**: 1 semana  
**Criterio de salida**: D2 + D5 + D6 pasan; B2 y B3 corregidos.

### Tarea B2: Rampa de fuerzas

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `solve()` — en el bloque de lectura de fuerzas

```python
# Después de leer F_aero, antes de _solve_fsi_step():
if self._force_ramp_time > 0.0 and t < self._force_ramp_time:
    ramp = t / self._force_ramp_time
    F_aero = F_aero * ramp

if self._force_max_magnitude is not None:
    f_norm = np.linalg.norm(F_aero)
    if f_norm > self._force_max_magnitude:
        F_aero = F_aero * (self._force_max_magnitude / f_norm)
```

Estos parámetros ya existen en `RotorConfig` (`force_ramp_time`, `force_max_magnitude`) y en `_init_rotor_config()` del solver corotacional — solo falta propagarlos al solver inertial.

---

### Tarea B3: Normalización de θ

**Archivo**: `src/aeroelast/solvers/fsi/rotor_inertial.py`  
**Método**: `solve()` — en el bloque post-convergencia de ventana

```python
theta += omega * dt

# Prevent unbounded accumulation: normalize every 10 full rotations.
# sin/cos precision degrades for large arguments; beyond 20π the roundoff
# in Rodrigues' formula starts to matter.
_TWO_PI = 2.0 * np.pi
if theta > 20.0 * _TWO_PI:
    theta = theta % _TWO_PI
```

---

### Tarea D2: Test de conservación de energía

```python
def test_energy_conservation_undamped_free_vibration():
    """
    Undamped free vibration (F_aero=0, gravity=0, no Rayleigh damping)
    should conserve total mechanical energy within Newmark truncation error.
    For β=0.25, γ=0.5 (trapezoidal rule) energy is exactly conserved in
    the linear case.
    """
    solver = build_minimal_inertial_solver(omega=0.0, damping=None)
    energies = []
    for _ in range(100):
        u, v, a = solver.solve_one_step(dt=0.001)
        T = 0.5 * v @ solver.M_diag * v   # kinetic
        U = 0.5 * u @ (solver.K @ u)      # strain
        energies.append(T + U)

    E0 = energies[0]
    drift = max(abs(E - E0) / (E0 + 1e-30) for E in energies[1:])
    assert drift < 1e-6, f"Energy drift {drift:.2e} exceeds tolerance for undamped Newmark"
```

---

### Tarea D5: Test de torque gravitacional 1P

```python
def test_gravity_torque_is_1p():
    """
    For a horizontal-axis rotor at constant ω, gravity produces a torque
    τ(θ) = M_total · g · R_cm · sin(θ) — a pure 1P signal.
    """
    solver = build_minimal_inertial_solver(
        omega=1.0,
        gravity=[0.0, 0.0, -9.81],
        rotation_axis=[0.0, 1.0, 0.0],
        F_aero=None,
    )
    torques = []
    for k in range(64):   # one full revolution at Δθ = 2π/64
        theta_k = k * 2 * np.pi / 64
        tau = solver.compute_gravity_torque_at(theta_k)
        torques.append(tau)

    spectrum = np.abs(np.fft.rfft(torques))
    dominant = int(np.argmax(spectrum[1:])) + 1   # skip DC
    assert dominant == 1, f"Expected 1P gravity torque, dominant frequency is {dominant}P"
```

---

### Tarea D6: Test de checkpoint/rollback

```python
def test_rollback_restores_consistent_state():
    """After rollback, K and C must match the geometry at checkpoint theta."""
    solver = build_minimal_inertial_solver(omega=5.0)

    # Advance one step to get a non-trivial state
    solver.step(dt=0.01)
    theta_before = solver._theta

    # Simulate checkpoint + two steps forward
    ckpt = solver.checkpoint()
    solver.step(dt=0.01)
    solver.step(dt=0.01)
    solver.rollback(ckpt)

    # K must correspond to theta_before, not theta_before + 2·Δθ
    K_restored = solver.K_current
    solver._rotate_structural_geometry_internal(theta_before)
    solver._rebuild_assembler_with_rotated_geometry()
    K_expected = solver.domain.assemble_stiffness_matrix()

    diff = (K_restored - K_expected).norm(PETSc.NormType.FROBENIUS)
    assert diff < 1e-10, f"K after rollback differs from expected: ‖ΔK‖ = {diff:.3e}"
```

---

## Sprint 3 — Performance crítica (factorización K_eff)

**Objetivo**: Eliminar el 80% de tiempo desperdiciado en factorizaciones repetidas.  
**Duración estimada**: 1 semana  
**Criterio de salida**: Benchmark 10k DOFs muestra speedup ≥ 4× vs baseline.

### Tarea C1: Cache de factorización K_eff por ventana

**Diseño**:
- `K_eff` y su factorización solo cambian al inicio de cada ventana (cuando cambia θ).
- Dentro de la ventana (sub-iterations de preCICE), solo cambia el RHS (`F_eff`).
- La factorización LU se realiza **una sola vez** al inicio de la ventana.

**Cambios en `solve()`**:

```python
# Fuera del while — cache vacío al inicio
_ksp_cache: Optional["PETSc.KSP"] = None

while adapter.is_coupling_ongoing:
    if adapter.requires_writing_checkpoint:
        checkpoint = ...
        _ksp_cache = None   # Nueva ventana → invalidar cache

    u_e_new, v_e_new, a_e_new = self._solve_fsi_step(
        ...,
        K_current, C_current,
        ksp_cache=_ksp_cache,
    )
    if _ksp_cache is None:
        _ksp_cache = self._last_ksp   # captura la primera factorización
```

**Cambios en `_solve_fsi_step()`**:

```python
def _solve_fsi_step(self, ..., ksp_cache=None):
    # ...
    F_eff = self._assemble_inertial_rhs(...)

    if ksp_cache is None:
        # Primera sub-iteration: ensamblar y factorizar
        K_eff = self._assemble_inertial_effective_system(K, C, M, a0, a1)
        bc_manager.apply_dirichlet_to_system(K_eff, F_eff)
        ksp = PETSc.KSP().create(comm=self.comm)
        ksp.setOperators(K_eff)
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        ksp.getPC().setFactorSolverType("mumps")   # sparse LU
        ksp.setFromOptions()
        ksp.setUp()   # factoriza aquí
        self._last_ksp = ksp
    else:
        # Sub-iterations siguientes: reutilizar factorización
        ksp = ksp_cache
        bc_manager.apply_dirichlet_to_rhs(F_eff)   # solo RHS

    u_e_new = PETSc.Vec().createMPI(self.domain.dofs_count, comm=self.comm)
    ksp.solve(F_eff, u_e_new)
    # ...
```

**Ganancia esperada**: 5-10× por ventana cuando preCICE usa coupling implícito.

---

## Sprint 4 — Performance mayor (Rust API)

**Objetivo**: Eliminar reconstrucción de topología en el assembler.  
**Duración estimada**: 2 semanas  
**Criterio de salida**: Benchmark 10k DOFs muestra speedup adicional ≥ 2× sobre Sprint 3.

### Tarea C2: API Rust `update_node_coordinates()`

**Archivo Rust**: `crates/aeroelast-core/src/assembly/assembler.rs`  
**Binding Python**: `crates/aeroelast-py/src/lib.rs`

```rust
/// Update node positions without rebuilding topology or material state.
/// Only element Jacobians, area vectors, and derived geometric quantities
/// are recomputed. Connectivity, DOF mapping, and material tensors are preserved.
pub fn update_node_coordinates(&mut self, coords: Vec<[f64; 3]>) -> PyResult<()> {
    if coords.len() != self.mesh.n_nodes() {
        return Err(PyValueError::new_err(format!(
            "Expected {} node coordinates, got {}",
            self.mesh.n_nodes(),
            coords.len()
        )));
    }
    self.mesh.update_positions(&coords);
    self.refresh_element_geometry();   // Jacobians, normals, areas
    Ok(())
}
```

**En Python** (`rotor_inertial.py`):

```python
def _rotate_structural_geometry_internal(self, theta: float) -> None:
    R = self._coord_transforms.rotation_matrix(theta)
    X_rotated = (R @ self._reference_coords.T).T   # (N, 3)
    
    if hasattr(self.domain.assembler, "update_node_coordinates"):
        # Fast path: update coords only, preserve topology
        self.domain.assembler.update_node_coordinates(X_rotated.tolist())
    else:
        # Fallback: full rebuild (Phase 1 behavior, for backward compatibility)
        self.domain.mesh.set_node_coordinates(X_rotated)
        self._rebuild_assembler_with_rotated_geometry()
```

**Ganancia esperada**: 2-4× adicional sobre Sprint 3.

---

## Sprint 5 — Validación completa y decisión de producto

**Objetivo**: Datos suficientes para decidir si promover a producción.  
**Duración estimada**: 1-2 semanas

### Tarea A2: Test de K(θ) para ortotrópicos

```python
def test_orthotropic_stiffness_transforms_under_rotation():
    """
    For an orthotropic material with E1 ≠ E2, K(θ=90°) must reflect
    the swapped principal directions. Eigenvalues should be preserved
    (rotation-invariant), but eigenvectors must rotate.
    """
    # Material: E1=100 GPa, E2=10 GPa
    solver_0   = build_orthotropic_inertial_solver(theta=0.0)
    solver_90  = build_orthotropic_inertial_solver(theta=np.pi / 2)

    K_0  = solver_0.domain.assemble_stiffness_matrix().getDenseArray()
    K_90 = solver_90.domain.assemble_stiffness_matrix().getDenseArray()

    eigs_0  = np.sort(np.linalg.eigvalsh(K_0))
    eigs_90 = np.sort(np.linalg.eigvalsh(K_90))

    # Eigenvalues invariant under rotation (isospectral property for SO(3))
    assert np.allclose(eigs_0, eigs_90, rtol=1e-8), (
        "K eigenvalues must be rotation-invariant for orthotropic material"
    )
    # Stiffness matrices must differ (otherwise rotation has no effect)
    assert not np.allclose(K_0, K_90, atol=1e-6), (
        "K(0°) and K(90°) must differ for anisotropic material"
    )
```

### Benchmark A/B: inertial vs corotacional

Script `benchmarks/compare_solvers.py`:
- Misma malla, material, BCs, cargas
- Comparar: tip displacement, torque, frecuencias naturales, tiempo/ventana
- Output: tabla Markdown + plots de convergencia

**Criterio de decisión**:
```
SI  speedup_inertial < 1.5× corotacional (con Sprints 3-4)
AND error_tip_displacement < 5% vs corotacional (con K_G opcional)
AND todos los tests de D1-D7 pasan
ENTONCES promover a producción
SINO mantener como experimental
```

---

## Cronograma

| Sprint | Foco | Duración | Issues resueltos |
|--------|------|----------|-----------------|
| **1** | Correctitud física + tests críticos | 1 semana | A1, A3, A4, B1, D1, D3 |
| **2** | Robustez numérica + tests de física | 1 semana | B2, B3, B4, D2, D5, D6 |
| **3** | Performance: cache factorización | 1 semana | C1 |
| **4** | Performance: API Rust | 2 semanas | C2 |
| **5** | Validación completa + decisión | 1-2 semanas | A2, D4, D7, benchmark |

**Total**: 6-7 semanas, 1 desarrollador senior.

---

## Estado para Iniciar

**Rama base**: `rust_implementation` (commit `9b30537`)  
**Rama de trabajo sugerida**: `fix/inertial-physics-correctness`

```bash
git checkout -b fix/inertial-physics-correctness
```

**Orden de commits sugerido por sprint**:

```
Sprint 1:
  fix(fsi): validate mass matrix is lumped in inertial solver            [A3]
  fix(fsi): warn on incompatible rotor flags in inertial solver          [A4]
  feat(fsi): add precice_displacement_mode config (elastic|total)        [A1]
  fix(fsi): include K and C matrices in FSI checkpoint state             [B1]
  test(fsi): rigid rotation produces zero elastic displacement           [D1]
  test(fsi): precice displacement mode contract elastic and total        [D3]

Sprint 2:
  fix(fsi): apply force ramp and cap in inertial solver                  [B2]
  fix(fsi): normalize theta accumulation to prevent float overflow       [B3]
  fix(fsi): reorthogonalize rotation matrix every 10 revolutions         [B4]
  test(fsi): energy conservation for undamped free vibration             [D2]
  test(fsi): gravity torque is 1P for horizontal-axis rotor              [D5]
  test(fsi): rollback restores consistent K and C matrices               [D6]

Sprint 3:
  perf(fsi): cache K_eff factorization across FSI sub-iterations        [C1]

Sprint 4:
  feat(core/rust): add update_node_coordinates to PyMeshAssembler        [C2]
  perf(fsi): use update_node_coordinates fast path for geometry rotation [C2]

Sprint 5:
  test(fsi): K(theta) eigenvalues invariant for orthotropic material     [A2]
  bench(fsi): add corotational vs inertial benchmark script
```
