from __future__ import annotations

import csv
import json
import logging
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

from aeroelast.cli.aero_report_info import main as aero_report_info_main
from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.solvers.aero import (
    AeroFSIRuntimeContext,
    AeroReportSchemaManifest,
    QuasiSteadyVLMBackend,
    RotorAeroFSIParticipant,
    RotorAeroLoads,
    RotorAeroState,
    build_aero_runtime_context,
    build_rotor_aero_participant,
    load_aero_report_schema,
    load_aero_report_table,
)


def _build_toy_rotor_mesh() -> MeshModel:
    mesh = MeshModel()

    blade_1_nodes = [
        Node([0.0, 0.0, 1.0]),
        Node([1.0, 0.0, 1.0]),
        Node([1.0, 1.0, 1.0]),
        Node([0.0, 1.0, 1.0]),
    ]
    blade_2_nodes = [
        Node([0.0, 0.0, -1.0]),
        Node([-1.0, 0.0, -1.0]),
        Node([-1.0, 1.0, -1.0]),
        Node([0.0, 1.0, -1.0]),
    ]

    for node in blade_1_nodes + blade_2_nodes:
        mesh.add_node(node)

    blade_1_elem = MeshElement(blade_1_nodes, ElementType.quad)
    blade_2_elem = MeshElement(blade_2_nodes, ElementType.quad)
    mesh.add_element(blade_1_elem)
    mesh.add_element(blade_2_elem)

    mesh.add_node_set(NodeSet("rotor_blade_1", set(blade_1_nodes)))
    mesh.add_node_set(NodeSet("rotor_blade_2", set(blade_2_nodes)))
    mesh.add_element_set(ElementSet("rotor_blade_1", {blade_1_elem}))
    mesh.add_element_set(ElementSet("rotor_blade_2", {blade_2_elem}))

    return mesh


def _build_runtime_context(
    *, backend: str = "bem", backend_config_overrides: dict | None = None
) -> AeroFSIRuntimeContext:
    mesh = _build_toy_rotor_mesh()
    backend_config = {
        "wind_speed": 12.0,
        "wind_direction": [0.0, 0.0, -1.0],
        "omega": {
            "initial": 1.5,
            "mesh": "OmegaMesh",
            "data": "AngularVelocity",
            "vertex": [1.0, 2.0, 3.0],
        },
        "azimuth": 15.0,
    }
    if backend_config_overrides:
        backend_config.update(backend_config_overrides)

    cfg = {
        "participant": "Fluid",
        "config_file": "precice-config.xml",
        "coupling_mesh": "Fluid-Mesh",
        "displacement_data": "Displacement",
        "force_data": "Force",
        "velocity_data": "Velocity",
        "aero": {
            "backend": backend,
            "blade_file": "blade.yaml",
            "rotor": {
                "n_blades": 2,
                "hub_radius": 3.5,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": backend_config,
        },
    }
    return build_aero_runtime_context(mesh, cfg)


def _build_panel_metadata(
    runtime_context: AeroFSIRuntimeContext,
    *,
    panel_forces: np.ndarray | None = None,
    panel_wind_velocities: np.ndarray | None = None,
    wake_row_count: int = 0,
    wake_corrector_iterations_used: int = 0,
) -> dict[str, np.ndarray | int]:
    surface = runtime_context.surface
    assert surface is not None

    collocation_points = np.vstack(
        [blade.panel_collocation_points for blade in surface.blades],
    )
    panel_count = collocation_points.shape[0]
    if panel_forces is None:
        panel_forces = np.zeros((panel_count, 3), dtype=float)
    if panel_wind_velocities is None:
        panel_wind_velocities = np.zeros((panel_count, 3), dtype=float)

    return {
        "panel_forces": np.asarray(panel_forces, dtype=float).reshape(panel_count, 3),
        "panel_collocation_points": collocation_points,
        "panel_wind_velocities": np.asarray(panel_wind_velocities, dtype=float).reshape(
            panel_count, 3
        ),
        "wake_row_count": int(wake_row_count),
        "wake_corrector_iterations_used": int(wake_corrector_iterations_used),
    }


class _RecordingBackend:
    name = "recording"

    def __init__(self, returned_forces: np.ndarray | None = None):
        self.calls: list[tuple[AeroFSIRuntimeContext, RotorAeroState]] = []
        self._returned_forces = returned_forces

    def compute_loads(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        self.calls.append((context, state))
        if self._returned_forces is None:
            nodal_forces = np.zeros((context.geometry.rotor_mesh.node_count, 3), dtype=float)
        else:
            nodal_forces = self._returned_forces
        return RotorAeroLoads(nodal_forces=nodal_forces)


class _PrimingExportBackend(_RecordingBackend):
    name = "priming-export"

    def __init__(self):
        super().__init__()
        self.prime_calls: list[tuple[AeroFSIRuntimeContext, RotorAeroState, Path]] = []
        self.export_calls: list[Path] = []

    def prime_aero_mesh_export(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
        output_file: str | Path,
    ) -> None:
        output_path = Path(output_file)
        self.prime_calls.append((context, state, output_path))
        output_path.write_text("primed mesh")

    def export_aero_mesh(self, output_file: str | Path) -> None:
        self.export_calls.append(Path(output_file))


class _StateTrackingBackend(_RecordingBackend):
    name = "state-tracking"

    def __init__(self):
        super().__init__()
        self.checkpoint_token = 0
        self.finalized_windows = 0

    def capture_state(self) -> dict[str, int]:
        return {"token": self.checkpoint_token}

    def restore_state(self, checkpoint_state: object | None) -> None:
        if isinstance(checkpoint_state, dict):
            self.checkpoint_token = int(checkpoint_state["token"])

    def finalize_time_window(self) -> None:
        self.finalized_windows += 1
        self.checkpoint_token += 1


class _RunLoopBackend(_RecordingBackend):
    name = "run-loop"

    def __init__(self):
        super().__init__()
        self.finalized_windows = 0

    def compute_loads(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        self.calls.append((context, state))
        nodal_forces = np.zeros((context.geometry.rotor_mesh.node_count, 3), dtype=float)
        nodal_forces[:, 2] = -1.5
        metadata = _build_panel_metadata(
            context,
            panel_forces=np.array(
                [
                    [0.0, 0.0, -3.0],
                    [0.0, 0.0, -3.0],
                ],
                dtype=float,
            ),
            panel_wind_velocities=np.array(
                [
                    [10.0, -2.0, -1.5],
                    [4.0, -1.0, -0.5],
                ],
                dtype=float,
            ),
            wake_row_count=1,
            wake_corrector_iterations_used=0,
        )
        return RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.array([1.0, 2.0, 3.0], dtype=float),
            integrated_moment=np.array([4.0, 5.0, 6.0], dtype=float),
            metadata=metadata,
        )

    def finalize_time_window(self) -> None:
        self.finalized_windows += 1


class _RollbackRunLoopBackend(_RecordingBackend):
    name = "run-loop-rollback"

    def __init__(self):
        super().__init__()
        self.finalized_windows = 0
        self.restore_tokens: list[int] = []
        self.checkpoint_token = 0

    def capture_state(self) -> dict[str, int]:
        return {"token": self.checkpoint_token}

    def restore_state(self, checkpoint_state: object | None) -> None:
        token = 0
        if isinstance(checkpoint_state, dict):
            token = int(checkpoint_state.get("token", 0))
        self.checkpoint_token = token
        self.restore_tokens.append(token)

    def compute_loads(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        self.calls.append((context, state))
        magnitude = float(self.checkpoint_token + 1)
        nodal_forces = np.zeros((context.geometry.rotor_mesh.node_count, 3), dtype=float)
        nodal_forces[:, 2] = -magnitude
        self.checkpoint_token = 9

        metadata = _build_panel_metadata(
            context,
            panel_forces=np.array(
                [
                    [0.0, 0.0, -2.0 * magnitude],
                    [0.0, 0.0, -2.0 * magnitude],
                ],
                dtype=float,
            ),
            panel_wind_velocities=np.array(
                [
                    [10.0, -2.0, -1.5],
                    [4.0, -1.0, -0.5],
                ],
                dtype=float,
            ),
        )
        return RotorAeroLoads(nodal_forces=nodal_forces, metadata=metadata)

    def finalize_time_window(self) -> None:
        self.finalized_windows += 1
        self.checkpoint_token += 1


class _StubAdapter:
    last_instance = None

    def __init__(self, participant: str, config_file: str, coupling_meshes: dict[str, np.ndarray]):
        self.participant = participant
        self.config_file = config_file
        self.coupling_meshes = coupling_meshes
        self._dt = 0.25
        self._is_coupling_ongoing = True
        self._is_time_window_complete = False
        self._requires_reading_checkpoint = False
        self._requires_writing_checkpoint = False
        self.initialized = False
        self.finalized = False
        self.read_requests: list[tuple[str, str]] = []
        self.write_requests: list[tuple[str, str, np.ndarray]] = []
        self.advance_calls: list[float] = []
        self.stored_checkpoints: list[tuple] = []

        n_nodes = coupling_meshes["Fluid-Mesh"].shape[0]
        self._displacements = np.zeros((n_nodes, 3), dtype=float)
        self._displacements[1] = [0.1, -0.2, 0.3]
        self._velocities = np.zeros((n_nodes, 3), dtype=float)
        self._velocities[5] = [0.0, 0.4, 0.3]
        self._omega = np.array([2.25], dtype=float)

        type(self).last_instance = self

    def initialize(self) -> float:
        self.initialized = True
        return self._dt

    def read_data(self, mesh_name: str, data_name: str) -> np.ndarray:
        self.read_requests.append((mesh_name, data_name))
        if mesh_name == "Fluid-Mesh" and data_name == "Displacement":
            return self._displacements.copy()
        if mesh_name == "Fluid-Mesh" and data_name == "Velocity":
            return self._velocities.copy()
        if mesh_name == "OmegaMesh" and data_name == "AngularVelocity":
            return self._omega.copy()
        raise AssertionError(f"Unexpected read_data request: {(mesh_name, data_name)!r}")

    def write_data(self, mesh_name: str, data_name: str, data: np.ndarray) -> None:
        self.write_requests.append((mesh_name, data_name, np.asarray(data, dtype=float).copy()))

    def advance(self, dt: float) -> float:
        self.advance_calls.append(float(dt))
        self._is_time_window_complete = True
        self._is_coupling_ongoing = False
        return self._dt

    def store_checkpoint(self, states: tuple) -> None:
        self.stored_checkpoints.append(states)

    def retrieve_checkpoint(self):
        raise AssertionError("retrieve_checkpoint should not be called in this test")

    def finalize(self) -> None:
        self.finalized = True

    @property
    def is_coupling_ongoing(self) -> bool:
        return self._is_coupling_ongoing

    @property
    def is_time_window_complete(self) -> bool:
        return self._is_time_window_complete

    @property
    def requires_reading_checkpoint(self) -> bool:
        return self._requires_reading_checkpoint

    @property
    def requires_writing_checkpoint(self) -> bool:
        return self._requires_writing_checkpoint

    @property
    def dt(self) -> float:
        return self._dt


class _RollbackStubAdapter(_StubAdapter):
    def __init__(self, participant: str, config_file: str, coupling_meshes: dict[str, np.ndarray]):
        super().__init__(participant, config_file, coupling_meshes)
        self._requires_writing_checkpoint = True
        self._phase = 0
        self._checkpoint_state: tuple | None = None

    def advance(self, dt: float) -> float:
        self.advance_calls.append(float(dt))
        if self._phase == 0:
            self._phase = 1
            self._requires_writing_checkpoint = False
            self._requires_reading_checkpoint = True
            self._is_time_window_complete = False
            self._is_coupling_ongoing = True
        else:
            self._phase = 2
            self._requires_reading_checkpoint = False
            self._is_time_window_complete = True
            self._is_coupling_ongoing = False
        return self._dt

    def store_checkpoint(self, states: tuple) -> None:
        super().store_checkpoint(states)
        self._checkpoint_state = states

    def retrieve_checkpoint(self):
        assert self._checkpoint_state is not None
        return self._checkpoint_state


def test_build_rotor_aero_participant_uses_runtime_context_defaults():
    runtime_context = _build_runtime_context()
    backend = _RecordingBackend()

    participant = build_rotor_aero_participant(runtime_context, backend)

    assert participant.runtime_context is runtime_context
    assert participant.backend is backend
    assert participant._participant_name == "Fluid"
    assert participant._coupling_mesh == "Fluid-Mesh"
    assert participant._omega_mesh == "OmegaMesh"
    assert participant._omega_data == "AngularVelocity"
    np.testing.assert_allclose(participant._omega_vertex, [1.0, 2.0, 3.0])
    np.testing.assert_allclose(participant._wind_velocity, [0.0, 0.0, -12.0])
    assert participant._current_omega == pytest.approx(1.5)
    assert participant._azimuth == pytest.approx(15.0)


def test_rotor_aero_participant_builds_backend_state_from_full_rotor_fields():
    runtime_context = _build_runtime_context()
    backend = _RecordingBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)

    displacements = np.full((runtime_context.geometry.rotor_mesh.node_count, 3), 0.25)
    velocities = np.full((runtime_context.geometry.rotor_mesh.node_count, 3), -0.5)

    loads = participant.compute_iteration_loads(
        displacements,
        velocities=velocities,
        current_time=2.0,
        dt=0.1,
    )

    assert loads.nodal_forces.shape == (runtime_context.geometry.rotor_mesh.node_count, 3)
    assert len(backend.calls) == 1
    recorded_context, recorded_state = backend.calls[0]
    assert recorded_context is runtime_context
    assert recorded_state.time == pytest.approx(2.0)
    assert recorded_state.dt == pytest.approx(0.1)
    assert recorded_state.omega_rad_s == pytest.approx(1.5)
    assert recorded_state.azimuth_deg == pytest.approx(15.0)
    np.testing.assert_allclose(recorded_state.wind_velocity, [0.0, 0.0, -12.0])
    np.testing.assert_allclose(recorded_state.nodal_displacements, displacements)
    np.testing.assert_allclose(recorded_state.nodal_velocities, velocities)
    deformed_surface = recorded_state.extra["deformed_surface"]
    np.testing.assert_allclose(
        deformed_surface.get_blade(0).panel_centroids[0],
        runtime_context.surface.get_blade(0).panel_centroids[0] + [0.25, 0.25, 0.25],
    )
    assert recorded_state.extra["backend"] == "bem"
    assert recorded_state.extra["window_index"] == 0
    assert recorded_state.extra["iteration_index"] == 0
    assert recorded_state.extra["wind_velocity_provider"] is participant._wind_velocity_provider


def test_rotor_aero_participant_rejects_invalid_backend_force_shape():
    runtime_context = _build_runtime_context()
    backend = _RecordingBackend(returned_forces=np.zeros((3, 3), dtype=float))
    participant = build_rotor_aero_participant(runtime_context, backend)

    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    with pytest.raises(ValueError, match="nodal_forces must have shape"):
        participant.compute_iteration_loads(displacements)


def test_rotor_aero_participant_requires_rotor_geometry():
    runtime_context = AeroFSIRuntimeContext(
        backend="bem",
        participant_name="Fluid",
        coupling_mesh_name="Fluid-Mesh",
        blade_file="blade.yaml",
        rotor={"n_blades": 3},
        normalized_config={"bem": {"wind_speed": 10.0}},
        geometry=None,
        surface=None,
    )

    with pytest.raises(ValueError, match="full-rotor geometry"):
        RotorAeroFSIParticipant(runtime_context, _RecordingBackend())


def test_rotor_aero_participant_advances_azimuth_from_omega():
    runtime_context = _build_runtime_context()
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    participant._advance_azimuth(0.5)

    assert participant._azimuth == pytest.approx(57.9718346348)


def test_rotor_aero_participant_delegates_backend_checkpoint_hooks():
    runtime_context = _build_runtime_context()
    backend = _StateTrackingBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)

    assert participant._capture_backend_checkpoint_state() == {"token": 0}

    backend.checkpoint_token = 7
    participant._restore_backend_checkpoint_state({"token": 3})
    participant._finalize_backend_time_window()

    assert backend.checkpoint_token == 4
    assert backend.finalized_windows == 1


def test_rotor_aero_participant_interpolates_tabulated_wind_velocity():
    runtime_context = _build_runtime_context(
        backend_config_overrides={
            "wind": {
                "type": "table",
                "times": [0.0, 1.0, 2.0],
                "velocities": [
                    [0.0, 0.0, -12.0],
                    [0.0, 0.0, -18.0],
                    [0.0, 0.0, -6.0],
                ],
            }
        }
    )
    backend = _RecordingBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)
    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    participant.compute_iteration_loads(displacements, current_time=0.5, dt=0.1)
    participant.compute_iteration_loads(displacements, current_time=3.0, dt=0.1)

    first_state = backend.calls[0][1]
    second_state = backend.calls[1][1]
    np.testing.assert_allclose(first_state.wind_velocity, [0.0, 0.0, -15.0])
    np.testing.assert_allclose(second_state.wind_velocity, [0.0, 0.0, -6.0])


def test_rotor_aero_participant_loads_tabulated_wind_velocity_from_csv(tmp_path):
    wind_csv = tmp_path / "wind.csv"
    wind_csv.write_text("time,u,v,w\n0.0,0.0,0.0,-12.0\n1.0,2.0,0.0,-18.0\n2.0,-1.0,0.5,-6.0\n")

    runtime_context = _build_runtime_context(
        backend_config_overrides={
            "wind": {
                "type": "table",
                "file": str(wind_csv),
            }
        }
    )
    backend = _RecordingBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)
    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    participant.compute_iteration_loads(displacements, current_time=0.5, dt=0.1)
    participant.compute_iteration_loads(displacements, current_time=1.5, dt=0.1)

    first_state = backend.calls[0][1]
    second_state = backend.calls[1][1]
    np.testing.assert_allclose(first_state.wind_velocity, [1.0, 0.0, -15.0])
    np.testing.assert_allclose(second_state.wind_velocity, [0.5, 0.25, -12.0])


def test_rotor_aero_participant_accumulates_physical_time_independently_of_dt():
    runtime_context = _build_runtime_context()
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    participant._current_time = 1.25
    participant._window_count = 4
    participant._iteration_count = 3
    participant._advance_time_window(0.2)

    assert participant._current_time == pytest.approx(1.45)
    assert participant._window_count == 5
    assert participant._iteration_count == 0


def test_rotor_aero_participant_builds_power_law_shear_wind_field():
    runtime_context = _build_runtime_context(
        backend_config_overrides={
            "wind_velocity": [0.0, 0.0, -12.0],
            "hub_height": 10.0,
            "shear_exp": 1.0,
            "vertical_axis": [0.0, 0.0, 1.0],
        }
    )
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    velocities = participant._wind_velocity_provider.evaluate_at_points(
        np.array(
            [
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, -2.0],
            ],
            dtype=float,
        ),
        time=0.0,
    )

    np.testing.assert_allclose(velocities[1], [0.0, 0.0, -12.0])
    assert np.linalg.norm(velocities[0]) > np.linalg.norm(velocities[1])
    assert np.linalg.norm(velocities[2]) < np.linalg.norm(velocities[1])


def test_rotor_aero_participant_writes_aero_report_csv(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    displacements[2] = [0.1, -0.2, 0.3]
    velocities = np.zeros_like(displacements)
    velocities[5] = [0.0, 0.4, 0.3]
    nodal_forces = np.zeros_like(displacements)
    nodal_forces[:, 2] = -1.0
    metadata = _build_panel_metadata(
        runtime_context,
        panel_wind_velocities=np.array(
            [
                [10.0, -2.0, -1.5],
                [4.0, -1.0, -0.5],
            ],
            dtype=float,
        ),
        wake_row_count=2,
        wake_corrector_iterations_used=1,
    )

    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.array([1.0, 2.0, 3.0]),
            integrated_moment=np.array([4.0, 5.0, 6.0]),
            metadata=metadata,
        ),
        displacements=displacements,
        velocities=velocities,
        time=0.25,
        window_index=2,
    )

    assert participant._output_folder == tmp_path
    assert participant._log_interval == 3

    csv_path = tmp_path / "aero_report.csv"
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader)

    assert reader.fieldnames is not None
    assert "Blade 1 Fx [N]" in reader.fieldnames
    assert "Blade 2 Moment Mag [N.m]" in reader.fieldnames
    assert "Blade 1 Max Disp [m]" in reader.fieldnames
    assert "Blade 2 Mean Panel Wind Speed [m/s]" in reader.fieldnames
    assert "Mean Panel Wind Speed [m/s]" in reader.fieldnames
    assert row["Schema Version"] == "1.0"
    assert row["Report Type"] == "aero_report"
    assert row["Backend"] == "vlm"
    assert float(row["Time [s]"]) == pytest.approx(0.25)
    assert int(row["Window"]) == 2
    assert float(row["Fx [N]"]) == pytest.approx(1.0)
    assert float(row["Fy [N]"]) == pytest.approx(2.0)
    assert float(row["Fz [N]"]) == pytest.approx(3.0)
    assert float(row["Mx [N.m]"]) == pytest.approx(4.0)
    assert float(row["My [N.m]"]) == pytest.approx(5.0)
    assert float(row["Mz [N.m]"]) == pytest.approx(6.0)
    assert float(row["Force Mag [N]"]) == pytest.approx(np.linalg.norm([1.0, 2.0, 3.0]))
    assert float(row["Max Disp [m]"]) == pytest.approx(np.linalg.norm([0.1, -0.2, 0.3]))
    assert float(row["Max Vel [m/s]"]) == pytest.approx(0.5)
    assert float(row["Wind Z [m/s]"]) == pytest.approx(-12.0)
    assert float(row["Mean Panel Wind Speed [m/s]"]) == pytest.approx(
        np.mean([
            np.linalg.norm([10.0, -2.0, -1.5]),
            np.linalg.norm([4.0, -1.0, -0.5]),
        ])
    )
    assert float(row["Blade 1 Fz [N]"]) == pytest.approx(-4.0)
    assert float(row["Blade 2 Fz [N]"]) == pytest.approx(-4.0)
    assert float(row["Blade 1 Force Mag [N]"]) == pytest.approx(4.0)
    assert float(row["Blade 2 Force Mag [N]"]) == pytest.approx(4.0)
    assert float(row["Blade 1 Mx [N.m]"]) == pytest.approx(-2.0)
    assert float(row["Blade 2 Mx [N.m]"]) == pytest.approx(-2.0)
    assert float(row["Blade 1 My [N.m]"]) == pytest.approx(2.0)
    assert float(row["Blade 2 My [N.m]"]) == pytest.approx(-2.0)
    assert float(row["Blade 1 Moment Mag [N.m]"]) == pytest.approx(2.0 * np.sqrt(2.0))
    assert float(row["Blade 2 Moment Mag [N.m]"]) == pytest.approx(2.0 * np.sqrt(2.0))
    assert float(row["Blade 1 Max Disp [m]"]) == pytest.approx(np.linalg.norm([0.1, -0.2, 0.3]))
    assert float(row["Blade 2 Max Disp [m]"]) == pytest.approx(0.0)
    assert float(row["Blade 1 Max Vel [m/s]"]) == pytest.approx(0.0)
    assert float(row["Blade 2 Max Vel [m/s]"]) == pytest.approx(0.5)
    assert float(row["Blade 1 Mean Panel Wind Speed [m/s]"]) == pytest.approx(
        np.linalg.norm([10.0, -2.0, -1.5])
    )
    assert float(row["Blade 2 Mean Panel Wind Speed [m/s]"]) == pytest.approx(
        np.linalg.norm([4.0, -1.0, -0.5])
    )
    assert int(row["Wake Rows"]) == 2
    assert int(row["Wake Corrector Iterations"]) == 1

    manifest_path = tmp_path / "aero_report_schema.json"
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

    assert manifest["manifest_version"] == "1.0"
    assert manifest["backend"] == "vlm"
    assert manifest["n_blades"] == 2
    assert manifest["sectional_bins"] == 8
    assert manifest["reports"]["aero_report"]["file"] == "aero_report.csv"
    assert manifest["reports"]["aero_spanwise_report"]["fieldnames"][11] == "Radius Fraction [-]"
    assert (
        "Radius Center Fraction [-]" in manifest["reports"]["aero_sectional_report"]["fieldnames"]
    )


def test_rotor_aero_participant_writes_aero_spanwise_report_csv(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    metadata = _build_panel_metadata(
        runtime_context,
        panel_forces=np.array(
            [
                [1.0, 2.0, -3.0],
                [-1.0, 0.5, -4.0],
            ],
            dtype=float,
        ),
        panel_wind_velocities=np.array(
            [
                [10.0, -2.0, -1.5],
                [4.0, -1.0, -0.5],
            ],
            dtype=float,
        ),
        wake_row_count=2,
        wake_corrector_iterations_used=1,
    )

    participant._append_aero_spanwise_report(
        RotorAeroLoads(
            nodal_forces=np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float),
            metadata=metadata,
        ),
        time=0.25,
        window_index=2,
    )

    csv_path = tmp_path / "aero_spanwise_report.csv"
    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    assert rows[0]["Schema Version"] == "1.0"
    assert rows[0]["Report Type"] == "aero_spanwise_report"
    assert rows[0]["Backend"] == "vlm"
    assert int(rows[0]["Blade"]) == 1
    assert int(rows[1]["Blade"]) == 2
    assert int(rows[0]["Panel In Blade"]) == 1
    assert int(rows[1]["Panel In Blade"]) == 1
    assert float(rows[0]["Fx [N]"]) == pytest.approx(1.0)
    assert float(rows[1]["Fz [N]"]) == pytest.approx(-4.0)
    assert float(rows[0]["Wind Speed [m/s]"]) == pytest.approx(np.linalg.norm([10.0, -2.0, -1.5]))
    assert float(rows[0]["Radius [m]"]) == pytest.approx(
        runtime_context.surface.get_blade(0).panel_reference_radii[0]
    )
    assert float(rows[0]["Radius Fraction [-]"]) == pytest.approx(1.0)
    assert float(rows[1]["Spanwise [m]"]) == pytest.approx(
        runtime_context.surface.get_blade(1).panel_spanwise_coordinates[0]
    )
    assert float(rows[1]["Blade Span Fraction [-]"]) == pytest.approx(0.0)
    assert float(rows[0]["Area [m2]"]) == pytest.approx(
        runtime_context.surface.get_blade(0).panel_areas[0]
    )
    assert int(rows[0]["Wake Rows"]) == 2
    assert int(rows[0]["Wake Corrector Iterations"]) == 1


def test_rotor_aero_participant_writes_aero_sectional_report_csv(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
        "sectional_bins": 4,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    metadata = _build_panel_metadata(
        runtime_context,
        panel_forces=np.array(
            [
                [1.0, 2.0, -3.0],
                [-1.0, 0.5, -4.0],
            ],
            dtype=float,
        ),
        panel_wind_velocities=np.array(
            [
                [10.0, -2.0, -1.5],
                [4.0, -1.0, -0.5],
            ],
            dtype=float,
        ),
        wake_row_count=2,
        wake_corrector_iterations_used=1,
    )

    participant._append_aero_sectional_report(
        RotorAeroLoads(
            nodal_forces=np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float),
            metadata=metadata,
        ),
        time=0.25,
        window_index=2,
    )

    csv_path = tmp_path / "aero_sectional_report.csv"
    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    assert rows[0]["Schema Version"] == "1.0"
    assert rows[0]["Report Type"] == "aero_sectional_report"
    assert int(rows[0]["Blade"]) == 1
    assert int(rows[1]["Blade"]) == 2
    assert int(rows[0]["Section Bin"]) == 1
    assert int(rows[0]["Panel Count"]) == 1
    assert float(rows[0]["Fx [N]"]) == pytest.approx(1.0)
    assert float(rows[1]["Fz [N]"]) == pytest.approx(-4.0)
    assert float(rows[0]["Area [m2]"]) == pytest.approx(
        runtime_context.surface.get_blade(0).panel_areas[0]
    )
    assert float(rows[0]["Radius Center Fraction [-]"]) == pytest.approx(1.0)
    assert float(rows[0]["Blade Span Center Fraction [-]"]) == pytest.approx(0.0)
    assert float(rows[0]["Wind Speed [m/s]"]) == pytest.approx(np.linalg.norm([10.0, -2.0, -1.5]))
    assert int(rows[0]["Wake Rows"]) == 2


def test_rotor_aero_report_csv_header_is_stable_when_panel_wind_metadata_appears_later(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    nodal_forces = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata={},
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.0,
        window_index=1,
    )

    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata=_build_panel_metadata(
                runtime_context,
                panel_wind_velocities=np.array(
                    [
                        [10.0, -2.0, -1.5],
                        [4.0, -1.0, -0.5],
                    ],
                    dtype=float,
                ),
            ),
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.25,
        window_index=2,
    )

    csv_path = tmp_path / "aero_report.csv"
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    assert reader.fieldnames is not None
    assert len(rows) == 2
    assert "Blade 1 Mean Panel Wind Speed [m/s]" in reader.fieldnames
    assert rows[0]["Blade 1 Mean Panel Wind Speed [m/s]"] == ""
    assert float(rows[1]["Blade 1 Mean Panel Wind Speed [m/s]"]) == pytest.approx(
        np.linalg.norm([10.0, -2.0, -1.5])
    )


def test_load_aero_report_schema_returns_typed_report_paths(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
        "sectional_bins": 6,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    nodal_forces = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata={},
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.0,
        window_index=1,
    )

    schema = load_aero_report_schema(tmp_path)

    assert isinstance(schema, AeroReportSchemaManifest)
    assert schema.backend == "vlm"
    assert schema.n_blades == 2
    assert schema.sectional_bins == 6
    assert schema.manifest_path == tmp_path / "aero_report_schema.json"
    assert schema.aero_report.path == tmp_path / "aero_report.csv"
    assert schema.aero_spanwise_report.path == tmp_path / "aero_spanwise_report.csv"
    assert schema.aero_sectional_report.path == tmp_path / "aero_sectional_report.csv"
    assert schema.aero_report.fieldnames[0] == "Schema Version"
    assert schema.aero_spanwise_report.fieldnames[11] == "Radius Fraction [-]"
    assert "Radius Center Fraction [-]" in schema.aero_sectional_report.fieldnames


def test_load_aero_report_table_validates_manifest_field_order(tmp_path):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    nodal_forces = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata={},
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.0,
        window_index=1,
    )

    table = load_aero_report_table(tmp_path, "aero_report")
    schema = load_aero_report_schema(tmp_path)

    assert table.height == 1
    assert tuple(table.columns) == schema.aero_report.fieldnames


def test_aero_report_info_cli_prints_manifest_summary(tmp_path, capsys):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    nodal_forces = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata={},
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.0,
        window_index=1,
    )

    exit_code = aero_report_info_main([str(tmp_path)])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Backend: vlm" in captured.out
    assert "aero_report.csv" in captured.out
    assert "aero_sectional_report.csv" in captured.out


def test_aero_report_info_cli_prints_report_fields(tmp_path, capsys):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 3,
    }
    participant = build_rotor_aero_participant(runtime_context, _RecordingBackend())

    nodal_forces = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    participant._append_aero_report(
        RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=np.zeros(3, dtype=float),
            integrated_moment=np.zeros(3, dtype=float),
            metadata={},
        ),
        displacements=np.zeros_like(nodal_forces),
        velocities=np.zeros_like(nodal_forces),
        time=0.0,
        window_index=1,
    )

    exit_code = aero_report_info_main([str(tmp_path), "--report", "aero_report", "--fields"])
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "Schema Version" in captured.out
    assert "Wake Corrector Iterations" in captured.out


def test_rotor_aero_participant_run_emits_report_for_converged_window(
    tmp_path, monkeypatch, caplog
):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 1,
    }
    backend = _RunLoopBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)

    fake_base_module = ModuleType("aeroelast.solvers.fsi.base")
    fake_base_module.Adapter = _StubAdapter
    monkeypatch.setitem(sys.modules, "aeroelast.solvers.fsi.base", fake_base_module)

    with caplog.at_level(logging.INFO):
        participant.run()

    adapter = _StubAdapter.last_instance
    assert adapter is not None
    assert adapter.initialized is True
    assert adapter.finalized is True
    assert adapter.advance_calls == [0.25]
    assert len(adapter.write_requests) == 1

    mesh_name, data_name, written_forces = adapter.write_requests[0]
    assert mesh_name == "Fluid-Mesh"
    assert data_name == "Force"
    np.testing.assert_allclose(written_forces[:, 2], -1.5)

    assert backend.finalized_windows == 1
    assert len(backend.calls) == 1
    _, recorded_state = backend.calls[0]
    assert recorded_state.time == pytest.approx(0.0)
    assert recorded_state.dt == pytest.approx(0.25)
    assert recorded_state.omega_rad_s == pytest.approx(2.25)

    assert participant._current_time == pytest.approx(0.25)
    assert participant._window_count == 1
    assert participant._current_omega == pytest.approx(2.25)

    csv_path = tmp_path / "aero_report.csv"
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        row = next(reader)

    spanwise_csv_path = tmp_path / "aero_spanwise_report.csv"
    with spanwise_csv_path.open(newline="") as handle:
        spanwise_rows = list(csv.DictReader(handle))

    sectional_csv_path = tmp_path / "aero_sectional_report.csv"
    with sectional_csv_path.open(newline="") as handle:
        sectional_rows = list(csv.DictReader(handle))

    manifest_path = tmp_path / "aero_report_schema.json"
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

    assert row["Backend"] == "vlm"
    assert float(row["Time [s]"]) == pytest.approx(0.25)
    assert int(row["Window"]) == 1
    assert float(row["Omega [rad/s]"]) == pytest.approx(2.25)
    assert float(row["Max Disp [m]"]) == pytest.approx(np.linalg.norm([0.1, -0.2, 0.3]))
    assert float(row["Max Vel [m/s]"]) == pytest.approx(0.5)
    assert int(row["Wake Rows"]) == 1
    assert len(spanwise_rows) == 2
    assert int(spanwise_rows[0]["Blade"]) == 1
    assert int(spanwise_rows[1]["Blade"]) == 2
    assert len(sectional_rows) == 2
    assert sectional_rows[0]["Report Type"] == "aero_sectional_report"
    assert manifest["reports"]["aero_sectional_report"]["file"] == "aero_sectional_report.csv"
    assert "Window=1" in caplog.text
    assert "wake_rows=1" in caplog.text


def test_rotor_aero_participant_run_restores_checkpoint_before_converged_window(
    tmp_path, monkeypatch
):
    runtime_context = _build_runtime_context(backend="vlm")
    runtime_context.normalized_config["output"] = {
        "folder": str(tmp_path),
        "log_interval": 1,
    }
    backend = _RollbackRunLoopBackend()
    participant = build_rotor_aero_participant(runtime_context, backend)

    fake_base_module = ModuleType("aeroelast.solvers.fsi.base")
    fake_base_module.Adapter = _RollbackStubAdapter
    monkeypatch.setitem(sys.modules, "aeroelast.solvers.fsi.base", fake_base_module)

    participant.run()

    adapter = _RollbackStubAdapter.last_instance
    assert adapter is not None
    assert adapter.initialized is True
    assert adapter.finalized is True
    assert len(adapter.stored_checkpoints) == 1
    assert len(adapter.write_requests) == 2

    assert backend.restore_tokens == [0]
    assert backend.finalized_windows == 1
    assert len(backend.calls) == 2
    assert backend.calls[0][1].extra["iteration_index"] == 0
    assert backend.calls[1][1].extra["iteration_index"] == 1

    assert participant._current_time == pytest.approx(0.25)
    assert participant._window_count == 1
    assert participant._last_forces.shape == (
        runtime_context.geometry.rotor_mesh.node_count,
        3,
    )
    np.testing.assert_allclose(participant._last_forces[:, 2], -1.0)

    csv_path = tmp_path / "aero_report.csv"
    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    spanwise_csv_path = tmp_path / "aero_spanwise_report.csv"
    with spanwise_csv_path.open(newline="") as handle:
        spanwise_rows = list(csv.DictReader(handle))

    sectional_csv_path = tmp_path / "aero_sectional_report.csv"
    with sectional_csv_path.open(newline="") as handle:
        sectional_rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    assert float(rows[0]["Fz [N]"]) == pytest.approx(-8.0)
    assert int(rows[0]["Window"]) == 1
    assert len(spanwise_rows) == 2
    assert float(spanwise_rows[0]["Fz [N]"]) == pytest.approx(-2.0)
    assert len(sectional_rows) == 2
    assert float(sectional_rows[0]["Fz [N]"]) == pytest.approx(-2.0)


def test_vlm_participant_uses_shear_to_break_blade_symmetry():
    runtime_context = _build_runtime_context(
        backend="vlm",
        backend_config_overrides={
            "wind_velocity": [10.0, -2.0, -1.5],
            "hub_height": 10.0,
            "shear_exp": 1.0,
            "vertical_axis": [0.0, 0.0, 1.0],
        },
    )
    participant = build_rotor_aero_participant(
        runtime_context,
        QuasiSteadyVLMBackend(diagonal_regularization=1.0e-6),
    )
    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    loads = participant.compute_iteration_loads(displacements, current_time=0.0, dt=0.1)

    panel_wind_speeds = np.linalg.norm(loads.metadata["panel_wind_velocities"], axis=1)
    upper_blade_force = np.linalg.norm(loads.nodal_forces[:4], axis=1).sum()
    lower_blade_force = np.linalg.norm(loads.nodal_forces[4:], axis=1).sum()

    assert panel_wind_speeds[0] > panel_wind_speeds[1]
    assert upper_blade_force > lower_blade_force


def test_participant_primes_aerodynamic_mesh_before_first_compute(tmp_path):
    runtime_context = _build_runtime_context(backend="sharpy")
    backend = _PrimingExportBackend()
    output_file = tmp_path / "fluid_mesh.vtu"

    participant = RotorAeroFSIParticipant(
        runtime_context,
        backend,
        aero_mesh_file=output_file,
    )

    assert participant._aero_mesh_written is True
    assert output_file.read_text() == "primed mesh"
    assert backend.export_calls == []
    assert len(backend.prime_calls) == 1

    _, primed_state, primed_output = backend.prime_calls[0]
    assert primed_output == output_file
    assert primed_state.time == pytest.approx(0.0)
    assert primed_state.dt == pytest.approx(0.0)
    np.testing.assert_allclose(primed_state.nodal_displacements, 0.0)

    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)
    participant.compute_iteration_loads(displacements, current_time=0.0, dt=0.1)

    assert backend.export_calls == []


def test_vlm_participant_combines_csv_shear_and_persistent_wake(tmp_path):
    wind_csv = tmp_path / "transient_wind.csv"
    wind_csv.write_text("time,u,v,w\n0.0,9.0,-2.0,-1.0\n0.5,14.0,-3.0,-1.8\n1.0,8.0,-1.5,-0.7\n")

    runtime_context = _build_runtime_context(
        backend="vlm",
        backend_config_overrides={
            "wind": {
                "type": "table",
                "file": str(wind_csv),
                "hub_height": 10.0,
                "shear_exp": 1.0,
                "vertical_axis": [0.0, 0.0, 1.0],
            }
        },
    )
    displacements = np.zeros((runtime_context.geometry.rotor_mesh.node_count, 3), dtype=float)

    def _run_participant(wake_history_steps: int, wake_corrector_iterations: int):
        participant = build_rotor_aero_participant(
            runtime_context,
            QuasiSteadyVLMBackend(
                diagonal_regularization=1.0e-6,
                wake_history_steps=wake_history_steps,
                wake_corrector_iterations=wake_corrector_iterations,
            ),
        )
        loads_history = []
        for _ in range(4):
            loads = participant.compute_iteration_loads(
                displacements,
                current_time=participant._current_time,
                dt=0.25,
            )
            loads_history.append(loads)
            participant._finalize_backend_time_window()
            participant._advance_time_window(0.25)
        return loads_history

    memoryless_history = _run_participant(0, 0)
    persistent_history = _run_participant(2, 1)

    wake_row_counts = [loads.metadata["wake_row_count"] for loads in persistent_history]
    assert wake_row_counts == [1, 2, 2, 2]

    second_step_diff = np.linalg.norm(
        persistent_history[1].integrated_force - memoryless_history[1].integrated_force
    )
    assert second_step_diff > 1.0

    panel_wind_speeds = np.linalg.norm(
        persistent_history[1].metadata["panel_wind_velocities"], axis=1
    )
    assert panel_wind_speeds[0] > panel_wind_speeds[1]
    assert persistent_history[1].metadata["wake_corrector_iterations_used"] == 1
