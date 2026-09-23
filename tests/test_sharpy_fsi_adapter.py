from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from aeroelast.solvers.aero.sharpy_fsi_adapter import SharpyRotorFSICaseAdapter
from aeroelast.solvers.aero.types import AeroFSIRuntimeContext, RotorAeroState


def _build_surface_zeta(z_level: float) -> np.ndarray:
    zeta = np.zeros((3, 2, 3), dtype=float)
    for chord_index in range(2):
        for span_index in range(3):
            zeta[:, chord_index, span_index] = [float(chord_index), float(span_index), z_level]
    return zeta


def test_export_aero_mesh_writes_sharpy_aerogrid_quads(tmp_path):
    meshio = pytest.importorskip("meshio")
    adapter = SharpyRotorFSICaseAdapter(blade_file="blade.yaml", route=tmp_path)
    data = SimpleNamespace(
        aero=SimpleNamespace(
            timestep_info=[
                SimpleNamespace(
                    zeta=[
                        _build_surface_zeta(0.0),
                        _build_surface_zeta(1.0),
                    ]
                )
            ]
        )
    )
    output_file = tmp_path / "fluid_mesh.vtu"

    adapter.export_aero_mesh(data, output_file)
    preview_file = output_file.with_suffix(".vtk")

    mesh = meshio.read(output_file)
    preview_mesh = meshio.read(preview_file)
    triangle_blocks = [block.data for block in mesh.cells if block.type == "triangle"]
    preview_triangle_blocks = [
        block.data for block in preview_mesh.cells if block.type == "triangle"
    ]

    assert output_file.exists()
    assert preview_file.exists()
    assert mesh.points.shape == (12, 3)
    assert len(triangle_blocks) == 1
    assert triangle_blocks[0].shape == (8, 3)
    assert len(preview_triangle_blocks) == 1
    assert preview_triangle_blocks[0].shape == (8, 3)
    assert set(np.unique(mesh.points[:, 2]).tolist()) == {0.0, 1.0}
    assert np.array_equal(
        np.unique(mesh.cell_data_dict["surface_id"]["triangle"]),
        np.array([1, 2]),
    )
    assert np.array_equal(
        np.unique(preview_mesh.cell_data_dict["surface_id"]["triangle"]),
        np.array([1, 2]),
    )


def test_build_data_forwards_configured_collective_pitch(monkeypatch, tmp_path):
    captured: dict[str, float] = {}

    class _FakeRotor:
        def __init__(self):
            self.StructuralInformation = SimpleNamespace(num_node=3)

        def generate_h5_files(self, route, case):
            captured["generated_case"] = 1.0

    class _FakeSimulation:
        solvers = {"SHARPy": {"route": str(tmp_path), "case": "fsi_sharpy_vlm_fsi"}}

        def generate_solver_file(self):
            captured["solver_file"] = 1.0

        def generate_dyn_file(self, steps):
            captured["dyn_steps"] = float(steps)

    class _FakeNuMAD:
        def __init__(self, **kwargs):
            captured["numad_init"] = 1.0

        def _build_rotor(self, pitch_deg, **kwargs):
            captured["pitch_deg"] = float(pitch_deg)
            return _FakeRotor()

        def _build_simulation(self, **kwargs):
            return _FakeSimulation()

    fake_gc = SimpleNamespace(clean_test_files=lambda route, case: None)

    monkeypatch.setattr(
        "aeroelast.solvers.aero.sharpy_fsi_adapter.load_blade_aero",
        lambda *args, **kwargs: SimpleNamespace(),
    )
    monkeypatch.setattr(
        "aeroelast.solvers.aero.sharpy_fsi_adapter.NuMADSharpyRotorCaseAdapter",
        _FakeNuMAD,
    )
    monkeypatch.setattr(
        "aeroelast.solvers.aero.sharpy_fsi_adapter._load_sharpy_modules",
        lambda: (fake_gc, None, None),
    )

    runtime_context = AeroFSIRuntimeContext(
        backend="sharpy",
        participant_name="Fluid",
        coupling_mesh_name="Fluid-Mesh",
        blade_file="blade.yaml",
        rotor={"hub_radius": 3.0, "n_blades": 3, "rotation_axis": [0.0, 1.0, 0.0]},
        normalized_config={"aero": {"backend": "sharpy"}},
    )

    adapter = SharpyRotorFSICaseAdapter.from_runtime_context(
        runtime_context,
        {"route": str(tmp_path), "pitch_deg": 7.25},
    )
    monkeypatch.setattr(adapter, "_build_presharpy_data", lambda case_file: SimpleNamespace())
    monkeypatch.setattr(adapter, "_prime_sharpy_data", lambda data: None)
    monkeypatch.setattr(
        adapter,
        "_build_cache",
        lambda data, context: SimpleNamespace(blade_maps=()),
    )

    state = RotorAeroState(
        time=0.0,
        dt=0.1,
        omega_rad_s=1.0,
        wind_velocity=np.array([0.0, 10.0, 0.0], dtype=float),
    )
    adapter.build_data(runtime_context, state, {})

    assert captured["pitch_deg"] == pytest.approx(7.25)
