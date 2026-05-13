# Tasks: element-corotational-shell

## Review Workload Forecast

| Field | Value |
|-------|-------|
| Estimated changed lines | 600–900 |
| 400-line budget risk | High |
| Chained PRs recommended | Yes |
| Suggested split | PR 1 (Rust core) → PR 2 (assembler + solver wiring + Python) |
| Delivery strategy | ask-on-risk |
| Chain strategy | pending |

Decision needed before apply: Yes
Chained PRs recommended: Yes
Chain strategy: pending
400-line budget risk: High

### Suggested Work Units

| Unit | Goal | Likely PR | Notes |
|------|------|-----------|-------|
| 1 | SVD + corotational primitives in MITC4 and MITC3 | PR 1 | Base: main; self-contained, all element-level tests included |
| 2 | Assembler + solver wiring + Python config flags | PR 2 | Base: PR 1 branch; depends on Unit 1 exports |

---

## Phase 1: Foundation — SVD polar decomposition (MITC4)

- [ ] 1.1 **[RED]** `test_mitc4_corotational.rs`: add tests SC-001..SC-006 (SVD cases + old DB removal check)
- [ ] 1.2 **[GREEN]** `mitc4.rs`: replace `polar_decomposition()` with `polar_decomposition_svd()` using `nalgebra::Matrix3::svd`; fix det(R) sign convention
- [ ] 1.3 **[GREEN]** `mitc4.rs`: delete or gate-out the old Denman-Beavers code path
- [ ] 1.4 **[VERIFY]** `cargo test -p aeroelast-core` — SC-001..SC-006 green

## Phase 2: Core — Corotational frame and T matrix (MITC4)

- [ ] 2.1 **[RED]** `test_mitc4_corotational.rs`: add SC-010, SC-011 (T=I for zero/translation u)
- [ ] 2.2 **[GREEN]** `mitc4.rs`: implement `build_t24_deformed(coords: &[f64; 12], u_trans: &[f64; 12]) → Matrix24`; degenerate fallback to identity
- [ ] 2.3 **[RED]** add SC-012 (rigid rotation invariance of K_T)
- [ ] 2.4 **[GREEN]** `mitc4.rs`: implement `compute_kt_corotational(u_elem: &[f64; 24], coords: &[f64; 12]) → Mat24`
- [ ] 2.5 **[GREEN]** `mitc4.rs`: implement K_σ in deformed frame (reuse `compute_kg_local()` at deformed coords)
- [ ] 2.6 **[VERIFY]** SC-010..SC-012 green

## Phase 3: Core — Corotational primitives for MITC3

- [ ] 3.1 **[RED]** `test_mitc3_corotational.rs`: add SC-013 (rigid rotation invariance)
- [ ] 3.2 **[GREEN]** `mitc3.rs`: implement `build_t18_deformed(coords: &[f64; 9], u_trans: &[f64; 9]) → Matrix18`
- [ ] 3.3 **[GREEN]** `mitc3.rs`: implement `compute_kt_corotational(u_elem: &[f64; 18], coords: &[f64; 9]) → Mat18`
- [ ] 3.4 **[VERIFY]** SC-013 green

## Phase 4: Integration — Assembler

- [ ] 4.1 **[RED]** add assembler-level test: SC-020 (patch test MITC4), SC-021 (patch test MITC3)
- [ ] 4.2 **[GREEN]** `assembler.rs`: add `assemble_kt_corotational(u: &[f64]) → CooMatrix`; dispatch to element `compute_kt_corotational` in rayon loop
- [ ] 4.3 **[GREEN]** `assembler.rs`: ensure zero-displacement path returns same result as `assemble_kt()`
- [ ] 4.4 **[VERIFY]** SC-020, SC-021 green

## Phase 5: Integration — Solver and config wiring

- [ ] 5.1 **[RED]** `tests/test_config_corotational.py`: add SC-040, SC-041, SC-042
- [ ] 5.2 **[GREEN]** `src/aeroelast/core/config.py`: add `use_corotational_kt: bool = False` and `kt_update_freq: int = 1` to `RotorConfig`; update `to_dict()`
- [ ] 5.3 **[GREEN]** `rotor_fsi.rs`: add `use_corotational_kt: bool`, `kt_update_freq: usize` to `RotorFsiConfig`; conditional call in assembly loop (every `kt_update_freq` windows)
- [ ] 5.4 **[GREEN]** `rotor_inertial.rs`: same wiring as 5.3
- [ ] 5.5 **[VERIFY]** SC-040..SC-042 green; Python tests pass

## Phase 6: Verification and regression

- [ ] 6.1 **[RED]** add SC-030 (small-angle agreement < 0.1%), SC-031 (regression TL), SC-032 (large-angle benefit > 5%)
- [ ] 6.2 **[GREEN]** run full Rust test suite: `cargo test --manifest-path crates/Cargo.toml`
- [ ] 6.3 **[GREEN]** run full Python suite: `python -m pytest tests/ -q --tb=short --ignore=tests/test_blade_mesh.py --ignore=tests/test_rotor_inertial.py`
- [ ] 6.4 **[VERIFY]** all SC-001..SC-042 pass; no pre-existing test regressions
