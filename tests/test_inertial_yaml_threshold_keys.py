"""T4.1 — YAML key parsing for new spin-softening / K_G threshold fields.

Covers spec scenario R6:
  - New keys are parsed and defaults are applied when keys are absent.
  - Type errors (non-numeric threshold values) raise ``ValueError``.
  - ``include_spin_softening`` defaults to ``True``.
  - ``include_geometric_stiffness`` defaults to ``False``.
"""

from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# Optional imports — collection must not crash without PETSc / Rust
# ---------------------------------------------------------------------------

try:
    from aeroelast.solvers.fsi.rotor_inertial import LinearDynamicFSIRotorInertialSolver

    _HAS_WRAPPER = True
except (ImportError, OSError):
    _HAS_WRAPPER = False
    LinearDynamicFSIRotorInertialSolver = None  # type: ignore[assignment,misc]

_skip_wrapper = pytest.mark.skipif(
    not _HAS_WRAPPER,
    reason="LinearDynamicFSIRotorInertialSolver not importable (missing PETSc/deps)",
)


# ---------------------------------------------------------------------------
# Helper — build a minimal model_properties dict understood by _init_rotor_config
# ---------------------------------------------------------------------------

def _make_model_props(rotor_overrides: dict | None = None) -> dict:
    """Minimal model_properties that satisfies _init_rotor_config without a real mesh."""
    rotor = {
        "omega": 10.0,
        "rotation_axis": [0.0, 0.0, 1.0],
        "rotation_center": [0.0, 0.0, 0.0],
    }
    if rotor_overrides:
        rotor.update(rotor_overrides)
    return {
        "solver": {
            "rotor": rotor,
            "coupling": {
                "participant": "Solid",
                "config_file": "precice-config.xml",
                "coupling_mesh": "Solid-Mesh",
                "write_data": "Displacement",
                "read_data": "Force",
            },
        },
    }


def _make_solver_instance(rotor_overrides: dict | None = None):
    """Instantiate the solver wrapper bypassing any mesh / PETSc state."""
    props = _make_model_props(rotor_overrides)
    # Call _init_rotor_config directly, bypassing __init__ which requires a mesh.
    # We construct a bare object and inject the minimum attributes it needs.
    import types

    obj = object.__new__(LinearDynamicFSIRotorInertialSolver)
    obj.model_properties = props
    obj.solver_params = props["solver"]
    # _init_rotor_config also accesses self.domain and self.comm if it reaches
    # far enough; the tests below only exercise the first part that reads YAML
    # keys, so we inject stubs.
    obj.domain = types.SimpleNamespace(nodes=[], mesh=types.SimpleNamespace(node_id_to_index={}))
    obj._init_rotor_config()
    return obj


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@_skip_wrapper
class TestSpinSofteningDefault:
    """include_spin_softening defaults to True when absent from YAML."""

    def test_default_is_true(self):
        s = _make_solver_instance()
        assert s._include_spin_softening is True, (
            "include_spin_softening must default to True to preserve HEAD behavior"
        )

    def test_explicit_true_is_stored(self):
        s = _make_solver_instance({"include_spin_softening": True})
        assert s._include_spin_softening is True

    def test_explicit_false_is_stored(self):
        s = _make_solver_instance({"include_spin_softening": False})
        assert s._include_spin_softening is False


@_skip_wrapper
class TestGeometricStiffnessDefault:
    """include_geometric_stiffness defaults to False when absent from YAML."""

    def test_default_is_false(self):
        s = _make_solver_instance()
        assert s._include_geometric_stiffness is False

    def test_explicit_true_is_stored(self):
        s = _make_solver_instance({"include_geometric_stiffness": True})
        assert s._include_geometric_stiffness is True


@_skip_wrapper
class TestKspThresholdDefaults:
    """ksp_omega_rebuild_high/low default to 0.005 / 0.003."""

    def test_defaults(self):
        s = _make_solver_instance()
        assert s._ksp_omega_rebuild_high == pytest.approx(0.005)
        assert s._ksp_omega_rebuild_low == pytest.approx(0.003)

    def test_explicit_values_stored(self):
        s = _make_solver_instance({
            "ksp_omega_rebuild_high": 0.01,
            "ksp_omega_rebuild_low": 0.007,
        })
        assert s._ksp_omega_rebuild_high == pytest.approx(0.01)
        assert s._ksp_omega_rebuild_low == pytest.approx(0.007)

    def test_invalid_type_raises_value_error(self):
        with pytest.raises(ValueError, match="ksp_omega_rebuild_high"):
            _make_solver_instance({"ksp_omega_rebuild_high": "bad"})

    def test_invalid_low_type_raises_value_error(self):
        with pytest.raises(ValueError, match="ksp_omega_rebuild_low"):
            _make_solver_instance({"ksp_omega_rebuild_low": [0.003]})


@_skip_wrapper
class TestKgThresholdDefaults:
    """kg_omega_rebuild_high/low default to 0.005 / 0.003."""

    def test_defaults(self):
        s = _make_solver_instance()
        assert s._kg_omega_rebuild_high == pytest.approx(0.005)
        assert s._kg_omega_rebuild_low == pytest.approx(0.003)

    def test_explicit_values_stored(self):
        s = _make_solver_instance({
            "kg_omega_rebuild_high": 0.02,
            "kg_omega_rebuild_low": 0.01,
        })
        assert s._kg_omega_rebuild_high == pytest.approx(0.02)
        assert s._kg_omega_rebuild_low == pytest.approx(0.01)

    def test_invalid_type_raises_value_error(self):
        with pytest.raises(ValueError, match="kg_omega_rebuild_high"):
            _make_solver_instance({"kg_omega_rebuild_high": "wrong"})

    def test_invalid_low_type_raises_value_error(self):
        with pytest.raises(ValueError, match="kg_omega_rebuild_low"):
            _make_solver_instance({"kg_omega_rebuild_low": None})
