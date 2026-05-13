# Testing Capabilities — fem-shell

**Strict TDD Mode**: enabled
**Detected**: 2026-05-12

## Test Runner
- Command: `python -m pytest tests/ -q --tb=short --ignore=tests/test_blade_mesh.py --ignore=tests/test_rotor_inertial.py`
- Framework: pytest >= 9.0.3

## Test Layers
| Layer | Available | Tool |
|-------|-----------|------|
| Unit | ✅ | pytest |
| Integration | ✅ | pytest (parity/validation tests cover Python↔Rust boundary) |
| E2E | ❌ | — |

## Coverage
- Available: ✅
- Command: `python -m pytest tests/ --cov=src/aeroelast --cov-report=term-missing`
- Config: `pyproject.toml [tool.coverage.*]`, branch=true, source=src/aeroelast

## Quality Tools
| Tool | Available | Command |
|------|-----------|---------|
| Linter | ✅ | `ruff check` |
| Type checker | ❌ | — (no mypy/pyright config found) |
| Formatter | ✅ | `ruff format` (via ruff) |
| Rust lint | ✅ | `cargo clippy` (crates/) |
| Rust test | ✅ | `cargo test` (crates/) |

## Known Exclusions
- `tests/test_blade_mesh.py` — stale imports, excluded by default
- `tests/test_rotor_inertial.py` — stale imports, excluded by default
- `tests/test_ko2017_performance.py` — 8 pre-existing failures on coarse meshes, not regressions
- Slow/benchmark: deselect with `-m "not slow and not benchmark"`

## Strict TDD Rationale
Test runner (pytest) is present and configured. No explicit `strict_tdd: false` marker found. Default: `strict_tdd: true`.
