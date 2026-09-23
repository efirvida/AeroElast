from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aeroelast.solvers.aero.sharpy_backend import SharpyAeroBackend
from aeroelast.solvers.aero.types import AeroFSIRuntimeContext, RotorAeroLoads, RotorAeroState


class _FakeSharpySolver:
    instances: list["_FakeSharpySolver"] = []

    def __init__(self):
        self.data = None
        self.initialise_calls = []
        self.run_calls = []
        type(self).instances.append(self)

    def initialise(self, data, custom_settings=None, restart=False):
        self.data = data
        self.initialise_calls.append((data, dict(custom_settings or {}), bool(restart)))

    def run(self, **kwargs):
        self.run_calls.append(dict(kwargs))
        self.data.ran = True
        return self.data


class _FakeCaseAdapter:
    def __init__(self):
        self.build_calls = []
        self.update_calls = []
        self.extract_calls = []
        self.export_calls = []
        self.restored_states = []
        self.finalized = 0

    def build_data(self, context, state, settings):
        self.build_calls.append((context, state, dict(settings)))
        return SimpleNamespace(ran=False)

    def update_state(self, data, context, state):
        self.update_calls.append((data, context, state))

    def build_run_kwargs(self, data, context, state):
        return {"custom": state.extra.get("custom", "value")}

    def extract_loads(self, data, context, state):
        self.extract_calls.append((data, context, state))
        return RotorAeroLoads(
            nodal_forces=np.zeros((4, 3), dtype=float),
            metadata={"ran": data.ran},
        )

    def capture_state(self, data):
        return {"ran": data.ran}

    def restore_state(self, data, checkpoint_state):
        self.restored_states.append(checkpoint_state)
        if data is not None and isinstance(checkpoint_state, dict):
            data.ran = bool(checkpoint_state.get("ran", False))

    def finalize_time_window(self, data):
        self.finalized += 1

    def export_aero_mesh(self, data, output_file):
        self.export_calls.append((data, Path(output_file)))


class _FakeInitSettingsAdapter(_FakeCaseAdapter):
    def build_initialise_settings(self, data, context, state):
        return {"horseshoe": False, "num_cores": 3}


def _runtime_context() -> AeroFSIRuntimeContext:
    return AeroFSIRuntimeContext(
        backend="sharpy",
        participant_name="Fluid",
        coupling_mesh_name="Fluid-Mesh",
        blade_file="blade.yaml",
        rotor={"n_blades": 2},
        normalized_config={"aero": {"backend": "sharpy"}},
    )


def _state() -> RotorAeroState:
    return RotorAeroState(
        time=1.25,
        dt=0.1,
        omega_rad_s=2.0,
        wind_velocity=np.array([10.0, 0.0, 0.0]),
        extra={"custom": "kept"},
    )


def test_sharpy_backend_rejects_unknown_model():
    with pytest.raises(ValueError, match="Unsupported SHARPy aerodynamic model"):
        SharpyAeroBackend(model="panel")


def test_sharpy_backend_requires_case_adapter_before_importing_sharpy():
    backend = SharpyAeroBackend(model="uvlm")

    with pytest.raises(RuntimeError, match="case_adapter"):
        backend.compute_loads(_runtime_context(), _state())


def test_sharpy_uvlm_backend_delegates_to_step_solver(monkeypatch):
    _FakeSharpySolver.instances.clear()

    def fake_import_module(name: str):
        assert name == "sharpy.solvers.stepuvlm"
        return SimpleNamespace(StepUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(
        model="uvlm",
        solver_settings={"convection_scheme": 2},
        case_adapter=adapter,
    )

    loads = backend.compute_loads(_runtime_context(), _state())

    solver = _FakeSharpySolver.instances[-1]
    assert loads.metadata["ran"] is True
    assert adapter.build_calls[0][2]["convection_scheme"] == 2
    assert adapter.build_calls[0][2]["dt"] == 0.1
    assert solver.initialise_calls[0][1]["n_time_steps"] == 1
    assert solver.run_calls[0]["custom"] == "kept"
    assert solver.run_calls[0]["dt"] == 0.1
    assert solver.run_calls[0]["t"] == 1.25
    assert solver.run_calls[0]["convect_wake"] is True


def test_sharpy_vlm_backend_delegates_to_static_solver(monkeypatch):
    _FakeSharpySolver.instances.clear()

    def fake_import_module(name: str):
        assert name == "sharpy.solvers.staticuvlm"
        return SimpleNamespace(StaticUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    backend = SharpyAeroBackend(model="vlm", case_adapter=_FakeCaseAdapter())

    backend.compute_loads(_runtime_context(), _state())

    solver = _FakeSharpySolver.instances[-1]
    assert "dt" not in solver.initialise_calls[0][1]
    assert "convect_wake" not in solver.run_calls[0]


def test_sharpy_backend_forwards_checkpoint_hooks(monkeypatch):
    def fake_import_module(name: str):
        return SimpleNamespace(StepUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(model="uvlm", case_adapter=adapter)

    backend.compute_loads(_runtime_context(), _state())
    checkpoint = backend.capture_state()
    backend.restore_state(checkpoint)
    backend.finalize_time_window()

    assert checkpoint == {"data": {"ran": True}}
    assert adapter.restored_states == [{"ran": True}]
    assert adapter.finalized == 1


def test_sharpy_backend_checkpoint_before_first_compute_is_safe(monkeypatch):
    def fake_import_module(name: str):
        return SimpleNamespace(StepUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(model="uvlm", case_adapter=adapter)

    checkpoint = backend.capture_state()
    backend.restore_state({"data": {"ran": True}})
    backend.finalize_time_window()

    assert checkpoint == {"data": None}
    assert adapter.restored_states == []
    assert adapter.finalized == 0


def test_sharpy_backend_rebuild_per_step_reinitialises_solver(monkeypatch):
    _FakeSharpySolver.instances.clear()

    def fake_import_module(name: str):
        assert name == "sharpy.solvers.staticuvlm"
        return SimpleNamespace(StaticUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(
        model="vlm",
        case_adapter=adapter,
        rebuild_per_step=True,
    )

    backend.compute_loads(_runtime_context(), _state())
    backend.compute_loads(_runtime_context(), _state())

    assert len(_FakeSharpySolver.instances) == 2
    assert len(adapter.build_calls) == 2


def test_sharpy_backend_uses_adapter_initialise_settings(monkeypatch):
    _FakeSharpySolver.instances.clear()

    def fake_import_module(name: str):
        assert name == "sharpy.solvers.staticuvlm"
        return SimpleNamespace(StaticUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeInitSettingsAdapter()
    backend = SharpyAeroBackend(
        model="vlm",
        solver_settings={"horseshoe": True},
        case_adapter=adapter,
    )

    backend.compute_loads(_runtime_context(), _state())

    solver = _FakeSharpySolver.instances[-1]
    assert solver.initialise_calls[0][1]["horseshoe"] is False
    assert solver.initialise_calls[0][1]["num_cores"] == 3


def test_sharpy_backend_exports_aerodynamic_mesh_via_adapter(monkeypatch, tmp_path):
    def fake_import_module(name: str):
        assert name == "sharpy.solvers.staticuvlm"
        return SimpleNamespace(StaticUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(model="vlm", case_adapter=adapter)

    backend.compute_loads(_runtime_context(), _state())
    output_file = tmp_path / "fluid_mesh.vtu"
    backend.export_aero_mesh(output_file)

    assert adapter.export_calls == [(backend._data, output_file)]


def test_sharpy_backend_primes_aerodynamic_mesh_before_first_compute(monkeypatch, tmp_path):
    _FakeSharpySolver.instances.clear()

    def fake_import_module(name: str):
        assert name == "sharpy.solvers.staticuvlm"
        return SimpleNamespace(StaticUvlm=_FakeSharpySolver)

    from aeroelast.solvers.aero import sharpy_backend

    monkeypatch.setattr(sharpy_backend.importlib, "import_module", fake_import_module)
    adapter = _FakeCaseAdapter()
    backend = SharpyAeroBackend(model="vlm", case_adapter=adapter)
    context = _runtime_context()
    state = _state()
    output_file = tmp_path / "fluid_mesh.vtu"

    backend.prime_aero_mesh_export(context, state, output_file)
    backend.compute_loads(context, state)

    assert len(_FakeSharpySolver.instances) == 1
    assert len(adapter.build_calls) == 1
    assert adapter.export_calls == [(backend._data, output_file)]
