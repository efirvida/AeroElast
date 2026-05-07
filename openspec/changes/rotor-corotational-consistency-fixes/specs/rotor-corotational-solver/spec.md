# Delta for Rotor Corotational Solver

## Purpose

This specification defines what MUST be true after applying the `rotor-corotational-consistency-fixes` change to `LinearDynamicFSIRotorCorotationalSolver` and its Rust hot-loop counterpart. The change corrects five inconsistencies between the documented theoretical contract and the Rust implementation: centrifugal evaluation coordinates, Newmark RHS history term for the gyroscopic matrix, ω-rebuild policy coherence, Coriolis docstring, and per-window inertial-force recomputation.

---

## MODIFIED Requirements

### Requirement: Centrifugal Force Evaluation Coordinates

The centrifugal force computation in the FSI rotor hot-loop MUST select the nodal coordinate set based on the `include_ksp` flag:

- When `include_ksp = true`: centrifugal force MUST be evaluated at undeformed nodal coordinates $X_0$. The spin-softening stiffness $K_{SP}$ already accounts for the displacement-dependent linearisation on the LHS; evaluating the force at $X_0 + u$ simultaneously would double-count the $O(\omega^2 |u|)$ term.
- When `include_ksp = false`: centrifugal force MUST be evaluated at deformed nodal coordinates $X_0 + u$. No LHS counterpart exists, so the full geometrically-nonlinear correction must appear in the force.

The Euler force block has no LHS counterpart regardless of `include_ksp`; it MUST continue to be evaluated at deformed coordinates in both cases.

#### Scenario: K_SP on — centrifugal evaluated at undeformed coordinates

- GIVEN a rotating beam at constant ω with non-zero steady-state deflection u, assembled with `include_ksp = true`
- WHEN the centrifugal force vector is computed by the Rust hot-loop for any sub-iteration
- THEN each nodal centrifugal force component MUST equal $m_i \cdot \omega^2 \cdot r_{0,i,\perp}$, where $r_{0,i,\perp}$ is the perpendicular distance from the rotation axis at the UNDEFORMED position of node $i$
- AND the force MUST NOT include any contribution from the current displacement $u_i$

#### Scenario: K_SP off — centrifugal evaluated at deformed coordinates

- GIVEN the same rotating beam with the same ω and deflection u, but assembled with `include_ksp = false`
- WHEN the centrifugal force vector is computed
- THEN each nodal centrifugal force component MUST equal $m_i \cdot \omega^2 \cdot r_{i,\perp}(X_0 + u)$, including the deformation correction
- AND the result MUST differ from the `include_ksp = true` result by an amount proportional to $\omega^2 |u|$

#### Scenario: Regression — spin-up tip displacement matches reference

- GIVEN a validated reference case (blade-only static spin-up with `include_ksp = true`) whose expected tip displacement was computed with the correct formulation
- WHEN the corrected solver runs the same case
- THEN the steady-state tip displacement MUST match the reference within the tolerance agreed in the design phase
- AND the result MUST NOT shift by an amount attributable to the prior double-counting of $K_{SP}$ and deformed-coordinate centrifugal force

---

### Requirement: Newmark RHS History Term for Gyroscopic Matrix

(Note: design phase MUST line-verify `step()` end-to-end before this requirement is locked. If design confirms the term is already present under a different form, this requirement is replaced by a documentation-only requirement.)

When the gyroscopic (Coriolis) matrix $G_{cor}$ is included on the LHS of the effective stiffness $K_{eff}$ via the $a_1 G_{cor}$ contribution, the RHS history vector in `step()` MUST include the corresponding Newmark history term $G_{cor} \cdot (a_1 u_n + a_4 v_n + a_5 a_n)$, mirroring the treatment of the structural damping matrix $C$.

The Newmark coefficients $a_1$, $a_4$, $a_5$ are defined by the trapezoidal (average acceleration) rule and MUST be consistent with the values used for the $C$ history term.

#### Scenario: Pure-Coriolis system converges at second order

- GIVEN a synthetic rotating system where only the gyroscopic term is non-zero (no stiffness, no damping, no centrifugal), with a known analytical solution
- WHEN the Newmark stepper integrates the system at successively halved time steps
- THEN the global displacement error MUST decrease at $O(\Delta t^2)$ (second-order trapezoidal convergence)
- AND a run without the $G_{cor}$ history term on the RHS MUST show $O(\Delta t)$ error, demonstrating the correction is load-bearing

---

### Requirement: Omega-Rebuild Policy for K_G and K_SP

The rebuild triggers for $K_G$ (geometric stiffness) and $K_{SP}$ (spin-softening stiffness) MUST use the same form: a relative change predicate on $\omega^2$, defined as $|\Delta(\omega^2)| / \omega^2 > \text{threshold}$, with hysteresis. It is NOT acceptable for one matrix to use an absolute $\Delta\omega$ test while the other uses a relative $\Delta(\omega^2)/\omega^2$ test.

A single internal helper (`omega_changed_significantly` or equivalent) MUST gate both rebuild paths, with the same threshold pair.

#### Scenario: Low-ω small Δω — coherent non-rebuild

- GIVEN a solver running at $\omega_0 = 0.5$ rad/s with a per-step increment $\Delta\omega$ smaller than the rebuild threshold in relative terms for both matrices
- WHEN the step advances
- THEN neither $K_G$ nor $K_{SP}$ MUST be rebuilt
- AND the pre-fix behaviour (where the absolute-$\Delta\omega$ gate for $K_{SP}$ would trigger a rebuild at this small $\omega$) MUST no longer occur

#### Scenario: High-ω steady state — coherent stale

- GIVEN a solver running at steady $\omega = 50$ rad/s where $\Delta\omega = 0$ between steps
- WHEN the step advances
- THEN both $K_G$ and $K_{SP}$ MUST remain stale (no rebuild)
- AND a step with $\Delta(\omega^2)/\omega^2 > \text{omega\_rebuild\_rel\_high}$ MUST trigger a rebuild of both matrices on the same step

---

### Requirement: Coriolis Treatment Docstring in rotor.py

The docstring of `LinearDynamicFSIRotorCorotationalSolver` in `rotor.py` — specifically the "Theoretical Limitations & Risks" section — MUST describe the actual implementation: $G_{cor}$ is placed implicitly on the LHS of $K_{eff}$ (via the $a_1 G_{cor}$ contribution). The prior description of Coriolis as an "explicit lagged force on RHS" MUST be removed or corrected.

The updated docstring MUST cross-reference the Newmark RHS history term for $G_{cor}$ (Finding #2 above) once that term is confirmed or added.

---

## ADDED Requirements

### Requirement: Configurable Omega-Rebuild Thresholds on RotorConfig

`RotorConfig` MUST expose two new fields for the relative $\omega^2$ rebuild thresholds:

- `omega_rebuild_rel_high: float` — rebuild trigger threshold (default `0.005`)
- `omega_rebuild_rel_low: float` — hysteresis release threshold (default `0.003`)

These fields MUST be forwarded through `RotorConfig.to_dict()` to the Rust `RotorFSIConfig` struct. The defaults MUST match the current $K_G$ behaviour so that existing cases without explicit configuration observe no change in rebuild frequency.

### Requirement: Per-Window Inertial Force Field in Rust Result Struct

The per-converged-window result struct returned from the Rust FSI hot-loop MUST include the last-applied combined inertial force vector (centrifugal + Coriolis + Euler, or equivalent decomposition agreed in design). This vector is already computed each sub-iteration; it MUST be retained and returned rather than discarded.

The new field MUST be additive in the PyO3 binding: existing Python callers that do not reference the new field MUST continue to work without modification. No existing field may be removed or renamed.

The Python `_step_cb` in `rotor.py` MUST consume the inertial force from this result field instead of recomputing it independently. The recomputation block currently at `_step_cb` MUST be removed once the passthrough field is available.

#### Scenario: Diagnostic CSV columns are equivalent before and after

- GIVEN an FSI run producing diagnostic output (inertial-force CSV columns or equivalent) before the change
- WHEN the same run is executed after the change with the passthrough field in use
- THEN the inertial-force values in the diagnostic output MUST be numerically equivalent within $10^{-12}$ (absolute), accounting only for floating-point order-of-operations differences between the Rust and Python evaluation paths

### Requirement: Validated Rotor Benchmarks Must Not Regress

All validated reference cases currently passing in the test suite (blade-only static spin-up, dynamic FSI rotor benchmarks) MUST continue to pass after this change. Numerical drift attributable to the correction of the centrifugal evaluation coordinates (Finding #1) is expected and accepted; the tolerance for that drift MUST be agreed in the design phase and documented as the new reference baseline.

No previously passing test MUST transition to a failing state for any reason other than the intentional correction of Finding #1.

### Requirement: Wall-Clock Performance Must Not Regress

The benchmark suite execution wall-clock time MUST NOT increase by more than 5% relative to the pre-change baseline, measured on the same hardware. The unified ω-rebuild predicate MUST NOT introduce additional refactorizations on existing cases whose `RotorConfig` does not override the new threshold fields.

---

## Requirements Summary Table

| Requirement | Strength | Scenarios |
|-------------|----------|-----------|
| Centrifugal force evaluation coordinates | MUST | 3 |
| Newmark RHS history term for G_cor (conditional) | MUST (if confirmed by design) | 1 |
| Omega-rebuild policy coherence (K_G and K_SP) | MUST | 2 |
| Coriolis docstring in rotor.py | MUST | 0 (manual review) |
| Configurable omega-rebuild thresholds on RotorConfig | MUST | 0 |
| Per-window inertial force passthrough (additive PyO3 field) | MUST | 1 |
| Validated rotor benchmarks must not regress | MUST | 0 (existing test suite) |
| Wall-clock performance must not regress | MUST NOT exceed +5% | 0 (benchmark suite) |
