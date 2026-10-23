"""The rule-9 mesh-size guard (scripts/check_mesh_size.py) sees every way a blade size is written.

2026-10-22 inventory (odd/tasks/l3-physical-validation.md L30): the guard only matched a bare
``element_size = <v>`` and let three spellings through - a prefixed constant
(``GATE_ELEMENT_SIZE = 2.0``), a positional size passed to a mesh builder
(``build_assembly(0.5)``, how S-5 and S-6 stayed at 0.5 m) and a tuple of sizes
(``MESHES = (2.0, 1.0, 0.5)``).
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "check_mesh_size.py"
_spec = importlib.util.spec_from_file_location("check_mesh_size", SCRIPT)
guard = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(guard)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("blade = Blade(yaml, element_size=1.0)", {1.0}),
        ("element_size: 0.5\n", {0.5}),
        ("GATE_ELEMENT_SIZE = 2.0", {2.0}),
        ("_ELEMENT_SIZE = 1.0", {1.0}),
        ("asm, mesh = build_assembly(0.5)", {0.5}),
        ("mesh, props = build_mesh_and_properties(1.0, 'mid')", {1.0}),
        ("asm, mesh = _build_assembly(0.5)", {0.5}),
        ("MESHES = (2.0, 1.0, 0.5)", {2.0, 1.0, 0.5}),
        ("ELEMENT_SIZES = [0.5, 0.25]", {0.5}),
    ],
)
def test_guard_sees_every_coarse_spelling(text, expected):
    assert guard.coarse_values(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "element_size=0.25",
        "ELEMENT_SIZE = 0.15",
        "build_assembly(0.25, skin_reference='mid')",
        "MESHES = (0.25, 0.125)",
        "rayleigh_damping=(0.5, 1e-3)",
        "MESHES = [(40, 2), (80, 4), (160, 8)]",
        "CONVERGENCE_MESHES = [(4, 2), (8, 4)]",
        "n_blades = 3",
    ],
)
def test_guard_ignores_converged_sizes_and_unrelated_numbers(text):
    assert guard.coarse_values(text) == set()
