"""Backend-agnostic aerodynamic preCICE participant for full rotors.

This participant owns the generic coupling loop for full-rotor aerodynamic
backends. It reads structural kinematics from preCICE, builds a
``RotorAeroState`` with the full-rotor ordering preserved, delegates the load
computation to a ``RotorAerodynamicBackend``, and writes the assembled nodal
forces back to preCICE.
"""

from __future__ import annotations

import csv
import json
import logging
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

from .types import AeroFSIRuntimeContext, RotorAerodynamicBackend, RotorAeroLoads, RotorAeroState

logger = logging.getLogger(__name__)

AERO_REPORT_SCHEMA_VERSION = "1.0"
AERO_SPANWISE_REPORT_SCHEMA_VERSION = "1.0"
AERO_SECTIONAL_REPORT_SCHEMA_VERSION = "1.0"
AERO_REPORT_SCHEMA_MANIFEST_VERSION = "1.0"
AERO_REPORT_SCHEMA_MANIFEST_FILENAME = "aero_report_schema.json"

AERO_REPORT_FIELDNAMES_BASE = (
    "Schema Version",
    "Report Type",
    "Time [s]",
    "Window",
    "Backend",
    "Omega [rad/s]",
    "Azimuth [deg]",
    "Fx [N]",
    "Fy [N]",
    "Fz [N]",
    "Mx [N.m]",
    "My [N.m]",
    "Mz [N.m]",
    "Force Mag [N]",
    "Max Disp [m]",
    "Max Vel [m/s]",
    "Wind X [m/s]",
    "Wind Y [m/s]",
    "Wind Z [m/s]",
    "Mean Panel Wind X [m/s]",
    "Mean Panel Wind Y [m/s]",
    "Mean Panel Wind Z [m/s]",
    "Mean Panel Wind Speed [m/s]",
    "Wake Rows",
    "Wake Corrector Iterations",
)

AERO_REPORT_BLADE_FIELD_TEMPLATES = (
    "Blade {blade} Fx [N]",
    "Blade {blade} Fy [N]",
    "Blade {blade} Fz [N]",
    "Blade {blade} Force Mag [N]",
    "Blade {blade} Mx [N.m]",
    "Blade {blade} My [N.m]",
    "Blade {blade} Mz [N.m]",
    "Blade {blade} Moment Mag [N.m]",
    "Blade {blade} Max Disp [m]",
    "Blade {blade} Max Vel [m/s]",
    "Blade {blade} Mean Panel Wind X [m/s]",
    "Blade {blade} Mean Panel Wind Y [m/s]",
    "Blade {blade} Mean Panel Wind Z [m/s]",
    "Blade {blade} Mean Panel Wind Speed [m/s]",
)

AERO_SPANWISE_REPORT_FIELDNAMES = (
    "Schema Version",
    "Report Type",
    "Time [s]",
    "Window",
    "Backend",
    "Omega [rad/s]",
    "Azimuth [deg]",
    "Blade",
    "Blade Azimuth [deg]",
    "Panel In Blade",
    "Radius [m]",
    "Radius Fraction [-]",
    "Spanwise [m]",
    "Blade Span Fraction [-]",
    "Chordwise [m]",
    "Area [m2]",
    "Chord Length [m]",
    "Span Length [m]",
    "Collocation X [m]",
    "Collocation Y [m]",
    "Collocation Z [m]",
    "Fx [N]",
    "Fy [N]",
    "Fz [N]",
    "Force Mag [N]",
    "Wind X [m/s]",
    "Wind Y [m/s]",
    "Wind Z [m/s]",
    "Wind Speed [m/s]",
    "Wake Rows",
    "Wake Corrector Iterations",
)

AERO_SECTIONAL_REPORT_FIELDNAMES = (
    "Schema Version",
    "Report Type",
    "Time [s]",
    "Window",
    "Backend",
    "Omega [rad/s]",
    "Azimuth [deg]",
    "Blade",
    "Blade Azimuth [deg]",
    "Section Bin",
    "Panel Count",
    "Spanwise Start [m]",
    "Spanwise End [m]",
    "Spanwise Center [m]",
    "Blade Span Start Fraction [-]",
    "Blade Span End Fraction [-]",
    "Blade Span Center Fraction [-]",
    "Radius Start [m]",
    "Radius End [m]",
    "Radius Center [m]",
    "Radius Start Fraction [-]",
    "Radius End Fraction [-]",
    "Radius Center Fraction [-]",
    "Area [m2]",
    "Mean Chord Length [m]",
    "Mean Span Length [m]",
    "Mean Chordwise [m]",
    "Fx [N]",
    "Fy [N]",
    "Fz [N]",
    "Force Mag [N]",
    "Wind X [m/s]",
    "Wind Y [m/s]",
    "Wind Z [m/s]",
    "Wind Speed [m/s]",
    "Wake Rows",
    "Wake Corrector Iterations",
)


def _normalize_vector(vector: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    if arr.shape != (3,):
        raise ValueError(f"{name} must be a 3-component vector")
    norm = np.linalg.norm(arr)
    if np.isclose(norm, 0.0):
        raise ValueError(f"{name} cannot be the zero vector")
    return arr / norm


def _coerce_velocity(vector: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    if arr.shape != (3,):
        raise ValueError(f"{name} must be a 3-component vector")
    return arr


def _normalize_fraction(value: float, *, lower: float, upper: float) -> float:
    if np.isclose(upper, lower):
        return 0.0
    return float((float(value) - float(lower)) / (float(upper) - float(lower)))


def _normalize_radius_fraction(radius: float, *, tip_radius: float) -> float:
    if np.isclose(tip_radius, 0.0):
        return 0.0
    return float(float(radius) / float(tip_radius))


def _resolve_speed_direction_velocity(
    *,
    speed: float,
    direction: Sequence[float],
    name: str,
) -> np.ndarray:
    return float(speed) * _normalize_vector(direction, name=name)


class _ConstantWindVelocity:
    def __init__(self, velocity: Sequence[float]) -> None:
        self._velocity = _coerce_velocity(velocity, name="wind_velocity")

    def evaluate(self, time: float) -> np.ndarray:
        return self._velocity.copy()

    def evaluate_at_points(self, points: np.ndarray, time: float) -> np.ndarray:
        points_arr = np.asarray(points, dtype=float).reshape(-1, 3)
        return np.broadcast_to(self._velocity, (points_arr.shape[0], 3)).copy()


class _TabulatedWindVelocity:
    def __init__(self, times: Sequence[float], velocities: Sequence[Sequence[float]]) -> None:
        self._times = np.asarray(times, dtype=float).reshape(-1)
        self._velocities = np.asarray(velocities, dtype=float)

        if self._times.size == 0:
            raise ValueError("wind.times cannot be empty")
        if self._velocities.shape != (self._times.size, 3):
            raise ValueError("wind.velocities must have shape (n_samples, 3) matching wind.times")
        if self._times.size > 1 and np.any(np.diff(self._times) < 0.0):
            raise ValueError("wind.times must be monotonically increasing")

    def evaluate(self, time: float) -> np.ndarray:
        if self._times.size == 1:
            return self._velocities[0].copy()

        return np.array(
            [
                np.interp(time, self._times, self._velocities[:, component])
                for component in range(3)
            ],
            dtype=float,
        )

    def evaluate_at_points(self, points: np.ndarray, time: float) -> np.ndarray:
        points_arr = np.asarray(points, dtype=float).reshape(-1, 3)
        velocity = self.evaluate(time)
        return np.broadcast_to(velocity, (points_arr.shape[0], 3)).copy()


class _PowerLawShearWindVelocity:
    def __init__(
        self,
        base_provider,
        *,
        reference_height: float,
        shear_exp: float,
        vertical_axis: Sequence[float],
        reference_point: Sequence[float],
    ) -> None:
        if reference_height <= 0.0:
            raise ValueError("hub_height must be positive when using shear_exp")

        self._base_provider = base_provider
        self._reference_height = float(reference_height)
        self._shear_exp = float(shear_exp)
        self._vertical_axis = _normalize_vector(vertical_axis, name="vertical_axis")
        self._reference_point = np.asarray(reference_point, dtype=float).reshape(3)
        self._min_height = max(1.0e-6, 1.0e-6 * self._reference_height)

    def evaluate(self, time: float) -> np.ndarray:
        return self._base_provider.evaluate(time)

    def evaluate_at_points(self, points: np.ndarray, time: float) -> np.ndarray:
        points_arr = np.asarray(points, dtype=float).reshape(-1, 3)
        base_velocity = np.asarray(self._base_provider.evaluate(time), dtype=float).reshape(3)
        if np.allclose(base_velocity, 0.0):
            return np.zeros((points_arr.shape[0], 3), dtype=float)

        heights = (
            self._reference_height + (points_arr - self._reference_point) @ self._vertical_axis
        )
        effective_heights = np.maximum(heights, self._min_height)
        scale = np.power(effective_heights / self._reference_height, self._shear_exp)
        return scale[:, None] * base_velocity.reshape(1, 3)


def _extract_backend_cfg(runtime_context: AeroFSIRuntimeContext) -> dict:
    backend_cfg = runtime_context.normalized_config.get(runtime_context.backend, {})
    return dict(backend_cfg) if isinstance(backend_cfg, Mapping) else {}


def _extract_initial_omega_cfg(
    runtime_context: AeroFSIRuntimeContext,
) -> tuple[float, str | None, str, list[float]]:
    backend_cfg = _extract_backend_cfg(runtime_context)
    omega_cfg = backend_cfg.get("omega", 0.0)
    if isinstance(omega_cfg, Mapping):
        omega_initial = float(omega_cfg.get("initial", 0.0))
        omega_mesh = omega_cfg.get("mesh")
        omega_data = str(omega_cfg.get("data", "AngularVelocity"))
        omega_vertex = list(omega_cfg.get("vertex", [0.0, 0.0, 0.0]))
        return omega_initial, omega_mesh, omega_data, omega_vertex
    return float(omega_cfg), None, "AngularVelocity", [0.0, 0.0, 0.0]


def _extract_wind_velocity(runtime_context: AeroFSIRuntimeContext) -> np.ndarray:
    backend_cfg = _extract_backend_cfg(runtime_context)

    wind_velocity = backend_cfg.get("wind_velocity")
    if wind_velocity is not None:
        return np.asarray(wind_velocity, dtype=float).reshape(3)

    wind_speed = float(backend_cfg.get("wind_speed", 0.0))
    wind_direction = backend_cfg.get("wind_direction", [1.0, 0.0, 0.0])
    return _resolve_speed_direction_velocity(
        speed=wind_speed,
        direction=wind_direction,
        name="wind_direction",
    )


def _build_tabulated_wind_velocity(wind_cfg: Mapping[str, object]) -> _TabulatedWindVelocity:
    file_path = wind_cfg.get("file")
    if file_path is not None:
        path = Path(str(file_path))
        if not path.exists():
            raise FileNotFoundError(f"wind.file not found: {path}")

        delimiter = str(wind_cfg.get("delimiter", ","))
        time_column = str(wind_cfg.get("time_column", "time"))
        velocity_columns = wind_cfg.get("velocity_columns")
        if velocity_columns is None:
            velocity_column_names = ["u", "v", "w"]
        else:
            velocity_column_names = [str(name) for name in velocity_columns]
            if len(velocity_column_names) != 3:
                raise ValueError("wind.velocity_columns must contain exactly three column names")

        speed_column = wind_cfg.get("speed_column")
        direction = wind_cfg.get("direction", [1.0, 0.0, 0.0])

        times: list[float] = []
        velocities: list[np.ndarray] = []
        with path.open(newline="") as csv_file:
            reader = csv.DictReader(csv_file, delimiter=delimiter)
            if reader.fieldnames is None:
                raise ValueError(f"wind.file has no header row: {path}")

            for row_index, row in enumerate(reader, start=2):
                if time_column not in row:
                    raise ValueError(f"wind.file missing time column '{time_column}': {path}")
                times.append(float(row[time_column]))

                if speed_column is not None:
                    speed_column_name = str(speed_column)
                    if speed_column_name not in row:
                        raise ValueError(
                            f"wind.file missing speed column '{speed_column_name}': {path}"
                        )
                    velocities.append(
                        _resolve_speed_direction_velocity(
                            speed=float(row[speed_column_name]),
                            direction=direction,
                            name=f"wind.direction[row {row_index}]",
                        )
                    )
                    continue

                missing_columns = [name for name in velocity_column_names if name not in row]
                if missing_columns:
                    missing = ", ".join(missing_columns)
                    raise ValueError(f"wind.file missing velocity columns [{missing}]: {path}")

                velocities.append(
                    _coerce_velocity(
                        [row[name] for name in velocity_column_names],
                        name=f"wind.file velocity row {row_index}",
                    )
                )

        return _TabulatedWindVelocity(times, velocities)

    samples = wind_cfg.get("samples")
    if isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)):
        times: list[float] = []
        velocities: list[np.ndarray] = []
        for index, sample in enumerate(samples):
            if not isinstance(sample, Mapping):
                raise ValueError("wind.samples entries must be mappings")
            if "time" not in sample:
                raise ValueError("wind.samples entries must include time")
            times.append(float(sample["time"]))
            velocity = sample.get("velocity")
            if velocity is not None:
                velocities.append(
                    _coerce_velocity(velocity, name=f"wind.samples[{index}].velocity")
                )
                continue

            if "speed" not in sample:
                raise ValueError("wind.samples entries must provide either velocity or speed")
            direction = sample.get("direction", wind_cfg.get("direction", [1.0, 0.0, 0.0]))
            velocities.append(
                _resolve_speed_direction_velocity(
                    speed=float(sample["speed"]),
                    direction=direction,
                    name=f"wind.samples[{index}].direction",
                )
            )

        return _TabulatedWindVelocity(times, velocities)

    times = wind_cfg.get("times")
    if times is None:
        raise ValueError("wind table requires either samples or times")

    velocities = wind_cfg.get("velocities")
    if velocities is not None:
        return _TabulatedWindVelocity(times, velocities)

    speeds = wind_cfg.get("speeds")
    if speeds is None:
        raise ValueError("wind table requires velocities or speeds")

    direction = wind_cfg.get("direction", [1.0, 0.0, 0.0])
    return _TabulatedWindVelocity(
        times,
        [
            _resolve_speed_direction_velocity(
                speed=float(speed),
                direction=direction,
                name="wind.direction",
            )
            for speed in speeds
        ],
    )


def _build_wind_velocity_provider(runtime_context: AeroFSIRuntimeContext):
    backend_cfg = _extract_backend_cfg(runtime_context)
    wind_cfg = backend_cfg.get("wind")
    if runtime_context.geometry is not None:
        reference_point = runtime_context.geometry.rotation_center
    else:
        reference_point = runtime_context.rotor.get("rotation_center", [0.0, 0.0, 0.0])

    if isinstance(wind_cfg, Mapping):
        wind_type = str(wind_cfg.get("type", "constant")).strip().lower()
        if wind_type == "constant":
            velocity = wind_cfg.get("velocity")
            if velocity is not None:
                base_provider = _ConstantWindVelocity(velocity)
            else:
                base_provider = _ConstantWindVelocity(
                    _resolve_speed_direction_velocity(
                        speed=float(wind_cfg.get("speed", 0.0)),
                        direction=wind_cfg.get("direction", [1.0, 0.0, 0.0]),
                        name="wind.direction",
                    )
                )

        elif wind_type in {"table", "tabulated", "time_series"}:
            base_provider = _build_tabulated_wind_velocity(wind_cfg)
        else:
            raise ValueError(
                f"Unsupported wind.type '{wind_type}'. Supported values: constant, table"
            )

        shear_exp = float(wind_cfg.get("shear_exp", backend_cfg.get("shear_exp", 0.0)))
        if not np.isclose(shear_exp, 0.0):
            return _PowerLawShearWindVelocity(
                base_provider,
                reference_height=float(
                    wind_cfg.get("hub_height", backend_cfg.get("hub_height", 0.0))
                ),
                shear_exp=shear_exp,
                vertical_axis=wind_cfg.get(
                    "vertical_axis", backend_cfg.get("vertical_axis", [0.0, 0.0, 1.0])
                ),
                reference_point=reference_point,
            )

        return base_provider

    base_provider = _ConstantWindVelocity(_extract_wind_velocity(runtime_context))
    shear_exp = float(backend_cfg.get("shear_exp", 0.0))
    if np.isclose(shear_exp, 0.0):
        return base_provider

    return _PowerLawShearWindVelocity(
        base_provider,
        reference_height=float(backend_cfg.get("hub_height", 0.0)),
        shear_exp=shear_exp,
        vertical_axis=backend_cfg.get("vertical_axis", [0.0, 0.0, 1.0]),
        reference_point=reference_point,
    )


class RotorAeroFSIParticipant:
    """Generic preCICE participant for full-rotor aerodynamic backends."""

    def __init__(
        self,
        runtime_context: AeroFSIRuntimeContext,
        backend: RotorAerodynamicBackend,
        *,
        participant: str | None = None,
        config_file: str = "precice-config.xml",
        coupling_mesh: str | None = None,
        displacement_data: str = "Displacement",
        force_data: str = "Force",
        velocity_data: str | None = "Velocity",
        omega_mesh: str | None = None,
        omega_data: str = "AngularVelocity",
        omega_vertex: Sequence[float] = (0.0, 0.0, 0.0),
        omega_rad_s: float = 0.0,
        azimuth_deg: float = 0.0,
        wind_velocity: Sequence[float] = (0.0, 0.0, 0.0),
        wind_velocity_provider=None,
        state_extra: Mapping[str, object] | None = None,
        output_folder: str | Path | None = "aero_fsi_results",
        aero_mesh_file: str | Path | None = None,
        log_interval: int = 10,
        sectional_bins: int = 8,
    ) -> None:
        if runtime_context.geometry is None:
            raise ValueError(
                "RotorAeroFSIParticipant requires a runtime context with full-rotor geometry"
            )

        self._runtime_context = runtime_context
        self._backend = backend
        self._participant_name = participant or runtime_context.participant_name
        self._config_file = str(config_file)
        self._coupling_mesh = coupling_mesh or runtime_context.coupling_mesh_name
        self._displacement_data = displacement_data
        self._force_data = force_data
        self._velocity_data = velocity_data
        self._omega_mesh = omega_mesh
        self._omega_data = omega_data
        self._omega_vertex = np.asarray(omega_vertex, dtype=float).reshape(3)
        self._current_omega = float(omega_rad_s)
        self._azimuth = float(azimuth_deg)
        self._wind_velocity = np.asarray(wind_velocity, dtype=float).reshape(3)
        self._wind_velocity_provider = wind_velocity_provider
        self._state_extra = dict(state_extra or {})
        self._output_folder = None if output_folder is None else Path(output_folder)
        self._aero_mesh_file = None if aero_mesh_file is None else Path(aero_mesh_file)
        self._log_interval = max(1, int(log_interval))
        self._sectional_bins = max(1, int(sectional_bins))
        self._schema_manifest_written = False
        self._current_time = 0.0
        self._aero_mesh_written = False

        self._ref_coords = runtime_context.geometry.rotor_mesh.coords_array.copy()
        self._n_nodes = self._ref_coords.shape[0]

        self._last_forces = np.zeros((self._n_nodes, 3), dtype=float)
        self._window_count = 0
        self._iteration_count = 0

        self._maybe_prime_aero_mesh()

    @property
    def runtime_context(self) -> AeroFSIRuntimeContext:
        return self._runtime_context

    @property
    def backend(self) -> RotorAerodynamicBackend:
        return self._backend

    def _coerce_nodal_field(self, field: np.ndarray, *, name: str) -> np.ndarray:
        arr = np.asarray(field, dtype=float)
        if arr.shape != (self._n_nodes, 3):
            raise ValueError(f"{name} must have shape ({self._n_nodes}, 3), got {arr.shape}")
        return arr

    def _read_optional_nodal_data(self, adapter, data_name: str | None) -> np.ndarray | None:
        if data_name is None:
            return None

        raw = adapter.read_data(self._coupling_mesh, data_name)
        if raw is None or np.asarray(raw).size == 0:
            return None
        return self._coerce_nodal_field(raw, name=data_name)

    def _update_live_omega(self, adapter) -> None:
        if self._omega_mesh is None:
            return

        omega_arr = adapter.read_data(self._omega_mesh, self._omega_data)
        if omega_arr is None or np.asarray(omega_arr).size == 0:
            return

        self._current_omega = float(np.asarray(omega_arr, dtype=float).ravel()[0])

    def _build_backend_state(
        self,
        *,
        current_time: float,
        dt: float,
        displacements: np.ndarray,
        velocities: np.ndarray | None,
    ) -> RotorAeroState:
        wind_velocity = self._resolve_wind_velocity(current_time)
        extra = dict(self._state_extra)
        extra.setdefault("window_index", self._window_count)
        extra.setdefault("iteration_index", self._iteration_count)
        extra.setdefault("backend", self._runtime_context.backend)
        extra.setdefault("wind_velocity_provider", self._wind_velocity_provider)
        if self._runtime_context.surface is not None:
            extra.setdefault(
                "deformed_surface",
                self._runtime_context.surface.build_deformed_state(displacements),
            )

        return RotorAeroState(
            time=current_time,
            dt=dt,
            omega_rad_s=self._current_omega,
            azimuth_deg=self._azimuth,
            wind_velocity=wind_velocity,
            nodal_displacements=displacements.copy(),
            nodal_velocities=None if velocities is None else velocities.copy(),
            extra=extra,
        )

    def _resolve_wind_velocity(self, current_time: float) -> np.ndarray:
        if self._wind_velocity_provider is None:
            return self._wind_velocity.copy()

        self._wind_velocity = np.asarray(
            self._wind_velocity_provider.evaluate(current_time),
            dtype=float,
        ).reshape(3)
        return self._wind_velocity.copy()

    def compute_iteration_loads(
        self,
        displacements: np.ndarray,
        *,
        velocities: np.ndarray | None = None,
        current_time: float = 0.0,
        dt: float = 0.0,
    ) -> RotorAeroLoads:
        displacements = self._coerce_nodal_field(displacements, name="displacements")
        if velocities is not None:
            velocities = self._coerce_nodal_field(velocities, name="velocities")

        state = self._build_backend_state(
            current_time=current_time,
            dt=dt,
            displacements=displacements,
            velocities=velocities,
        )
        loads = self._backend.compute_loads(self._runtime_context, state)
        if not isinstance(loads, RotorAeroLoads):
            raise TypeError(
                f"Aerodynamic backend '{self._backend.name}' must return RotorAeroLoads"
            )

        self._coerce_nodal_field(loads.nodal_forces, name="nodal_forces")
        self._maybe_export_aero_mesh()
        return loads

    def _maybe_prime_aero_mesh(self) -> None:
        if self._aero_mesh_written or self._aero_mesh_file is None:
            return

        prime_aero_mesh_export = getattr(self._backend, "prime_aero_mesh_export", None)
        if not callable(prime_aero_mesh_export):
            return

        try:
            zero_displacements = np.zeros((self._n_nodes, 3), dtype=float)
            initial_state = self._build_backend_state(
                current_time=self._current_time,
                dt=0.0,
                displacements=zero_displacements,
                velocities=None,
            )
            prime_aero_mesh_export(
                self._runtime_context,
                initial_state,
                self._aero_mesh_file,
            )
        except Exception as exc:
            logger.warning(
                "[AERO-FSI] Failed to prime aerodynamic mesh export to %s: %s",
                self._aero_mesh_file,
                exc,
            )
            return

        self._aero_mesh_written = True

    def _maybe_export_aero_mesh(self) -> None:
        if self._aero_mesh_written or self._aero_mesh_file is None:
            return

        export_aero_mesh = getattr(self._backend, "export_aero_mesh", None)
        if not callable(export_aero_mesh):
            return

        try:
            export_aero_mesh(self._aero_mesh_file)
        except Exception as exc:
            logger.warning(
                "[AERO-FSI] Failed to export aerodynamic mesh to %s: %s",
                self._aero_mesh_file,
                exc,
            )
            return

        self._aero_mesh_written = True

    def _capture_backend_checkpoint_state(self) -> object | None:
        capture_state = getattr(self._backend, "capture_state", None)
        if callable(capture_state):
            return capture_state()
        return None

    def _restore_backend_checkpoint_state(self, checkpoint_state: object | None) -> None:
        restore_state = getattr(self._backend, "restore_state", None)
        if callable(restore_state):
            restore_state(checkpoint_state)

    def _finalize_backend_time_window(self) -> None:
        finalize_time_window = getattr(self._backend, "finalize_time_window", None)
        if callable(finalize_time_window):
            finalize_time_window()

    def _advance_azimuth(self, dt: float) -> None:
        self._azimuth += np.degrees(self._current_omega * dt)
        self._azimuth %= 360.0

    def _advance_time_window(self, dt: float) -> None:
        self._current_time += float(dt)
        self._window_count += 1
        self._iteration_count = 0
        self._advance_azimuth(dt)

    def _iter_blade_panel_slices(self) -> list[tuple[int, object, slice]]:
        surface = self._runtime_context.surface
        if surface is None:
            return []

        panel_slices: list[tuple[int, object, slice]] = []
        cursor = 0
        for blade_number, blade_surface in enumerate(surface.blades, start=1):
            count = blade_surface.n_panels
            panel_slices.append((blade_number, blade_surface, slice(cursor, cursor + count)))
            cursor += count
        return panel_slices

    def _aero_report_fieldnames(self) -> list[str]:
        fieldnames = list(AERO_REPORT_FIELDNAMES_BASE)
        for blade_number in range(1, self._runtime_context.geometry.n_blades + 1):
            fieldnames.extend(
                field.format(blade=blade_number) for field in AERO_REPORT_BLADE_FIELD_TEMPLATES
            )
        return fieldnames

    def _build_aero_report_schema_manifest(self) -> dict[str, object]:
        return {
            "manifest_version": AERO_REPORT_SCHEMA_MANIFEST_VERSION,
            "backend": str(self._runtime_context.backend),
            "n_blades": int(self._runtime_context.geometry.n_blades),
            "sectional_bins": int(self._sectional_bins),
            "reports": {
                "aero_report": {
                    "file": "aero_report.csv",
                    "schema_version": AERO_REPORT_SCHEMA_VERSION,
                    "granularity": "one row per converged time window",
                    "summary": "Rotor-global loads plus blade-integrated diagnostics.",
                    "fieldnames": self._aero_report_fieldnames(),
                },
                "aero_spanwise_report": {
                    "file": "aero_spanwise_report.csv",
                    "schema_version": AERO_SPANWISE_REPORT_SCHEMA_VERSION,
                    "granularity": "one row per aerodynamic panel",
                    "summary": "Per-panel aerodynamic loads, collocation coordinates, and local wind.",
                    "fieldnames": list(AERO_SPANWISE_REPORT_FIELDNAMES),
                },
                "aero_sectional_report": {
                    "file": "aero_sectional_report.csv",
                    "schema_version": AERO_SECTIONAL_REPORT_SCHEMA_VERSION,
                    "granularity": "one row per blade radial bin",
                    "summary": "Radial-bin aggregates derived from panel metadata.",
                    "fieldnames": list(AERO_SECTIONAL_REPORT_FIELDNAMES),
                },
            },
        }

    def _ensure_output_schema_manifest(self) -> None:
        if self._output_folder is None:
            return

        self._output_folder.mkdir(parents=True, exist_ok=True)
        if self._schema_manifest_written:
            return

        manifest_path = self._output_folder / AERO_REPORT_SCHEMA_MANIFEST_FILENAME
        with manifest_path.open("w", encoding="utf-8") as handle:
            json.dump(self._build_aero_report_schema_manifest(), handle, indent=2)
            handle.write("\n")
        self._schema_manifest_written = True

    def _build_aero_report_row(
        self,
        loads: RotorAeroLoads,
        *,
        displacements: np.ndarray,
        velocities: np.ndarray | None,
        time: float,
        window_index: int,
    ) -> dict[str, float | int | str]:
        integrated_force = (
            np.asarray(loads.integrated_force, dtype=float).reshape(3)
            if loads.integrated_force is not None
            else np.sum(loads.nodal_forces, axis=0)
        )
        integrated_moment = loads.integrated_moment
        if integrated_moment is None:
            lever_arms = self._ref_coords - self._runtime_context.geometry.rotation_center.reshape(
                1, 3
            )
            integrated_moment = np.cross(lever_arms, loads.nodal_forces).sum(axis=0)
        integrated_moment = np.asarray(integrated_moment, dtype=float).reshape(3)

        disp_mag = np.linalg.norm(displacements, axis=1)
        max_disp = float(np.max(disp_mag)) if disp_mag.size else 0.0

        if velocities is None:
            max_vel = 0.0
        else:
            vel_mag = np.linalg.norm(velocities, axis=1)
            max_vel = float(np.max(vel_mag)) if vel_mag.size else 0.0

        mean_panel_wind = self._wind_velocity.copy()
        mean_panel_wind_speed = float(np.linalg.norm(mean_panel_wind))
        panel_wind_velocities = loads.metadata.get("panel_wind_velocities")
        if panel_wind_velocities is not None and np.asarray(panel_wind_velocities).size:
            panel_wind_arr = np.asarray(panel_wind_velocities, dtype=float).reshape(-1, 3)
            mean_panel_wind = np.mean(panel_wind_arr, axis=0)
            mean_panel_wind_speed = float(np.mean(np.linalg.norm(panel_wind_arr, axis=1)))

        blade_panel_wind_stats: dict[int, tuple[np.ndarray, float]] = {}
        blade_panel_slices = self._iter_blade_panel_slices()
        if (
            panel_wind_velocities is not None
            and np.asarray(panel_wind_velocities).size
            and blade_panel_slices
            and np.asarray(panel_wind_velocities).reshape(-1, 3).shape[0]
            == blade_panel_slices[-1][2].stop
        ):
            panel_wind_arr = np.asarray(panel_wind_velocities, dtype=float).reshape(-1, 3)
            for blade_number, _blade_surface, blade_slice in blade_panel_slices:
                blade_panel_wind = panel_wind_arr[blade_slice]
                if blade_panel_wind.size == 0:
                    continue
                blade_panel_wind_stats[blade_number] = (
                    np.mean(blade_panel_wind, axis=0),
                    float(np.mean(np.linalg.norm(blade_panel_wind, axis=1))),
                )

        row = {
            "Schema Version": AERO_REPORT_SCHEMA_VERSION,
            "Report Type": "aero_report",
            "Time [s]": float(time),
            "Window": int(window_index),
            "Backend": str(self._runtime_context.backend),
            "Omega [rad/s]": float(self._current_omega),
            "Azimuth [deg]": float(self._azimuth),
            "Fx [N]": float(integrated_force[0]),
            "Fy [N]": float(integrated_force[1]),
            "Fz [N]": float(integrated_force[2]),
            "Mx [N.m]": float(integrated_moment[0]),
            "My [N.m]": float(integrated_moment[1]),
            "Mz [N.m]": float(integrated_moment[2]),
            "Force Mag [N]": float(np.linalg.norm(integrated_force)),
            "Max Disp [m]": max_disp,
            "Max Vel [m/s]": max_vel,
            "Wind X [m/s]": float(self._wind_velocity[0]),
            "Wind Y [m/s]": float(self._wind_velocity[1]),
            "Wind Z [m/s]": float(self._wind_velocity[2]),
            "Mean Panel Wind X [m/s]": float(mean_panel_wind[0]),
            "Mean Panel Wind Y [m/s]": float(mean_panel_wind[1]),
            "Mean Panel Wind Z [m/s]": float(mean_panel_wind[2]),
            "Mean Panel Wind Speed [m/s]": mean_panel_wind_speed,
            "Wake Rows": int(loads.metadata.get("wake_row_count", 0)),
            "Wake Corrector Iterations": int(
                loads.metadata.get("wake_corrector_iterations_used", 0)
            ),
        }

        rotation_center = self._runtime_context.geometry.rotation_center.reshape(1, 3)
        for blade_number, blade in enumerate(self._runtime_context.geometry.blades, start=1):
            blade_forces = loads.nodal_forces[blade.parent_node_indices]
            blade_integrated_force = np.sum(blade_forces, axis=0)
            blade_lever_arms = self._ref_coords[blade.parent_node_indices] - rotation_center
            blade_integrated_moment = np.cross(blade_lever_arms, blade_forces).sum(axis=0)
            blade_displacements = displacements[blade.parent_node_indices]
            blade_disp_mag = np.linalg.norm(blade_displacements, axis=1)
            blade_max_disp = float(np.max(blade_disp_mag)) if blade_disp_mag.size else 0.0
            if velocities is None:
                blade_max_vel = 0.0
            else:
                blade_velocities = velocities[blade.parent_node_indices]
                blade_vel_mag = np.linalg.norm(blade_velocities, axis=1)
                blade_max_vel = float(np.max(blade_vel_mag)) if blade_vel_mag.size else 0.0
            row.update({
                f"Blade {blade_number} Fx [N]": float(blade_integrated_force[0]),
                f"Blade {blade_number} Fy [N]": float(blade_integrated_force[1]),
                f"Blade {blade_number} Fz [N]": float(blade_integrated_force[2]),
                f"Blade {blade_number} Force Mag [N]": float(
                    np.linalg.norm(blade_integrated_force)
                ),
                f"Blade {blade_number} Mx [N.m]": float(blade_integrated_moment[0]),
                f"Blade {blade_number} My [N.m]": float(blade_integrated_moment[1]),
                f"Blade {blade_number} Mz [N.m]": float(blade_integrated_moment[2]),
                f"Blade {blade_number} Moment Mag [N.m]": float(
                    np.linalg.norm(blade_integrated_moment)
                ),
                f"Blade {blade_number} Max Disp [m]": blade_max_disp,
                f"Blade {blade_number} Max Vel [m/s]": blade_max_vel,
            })

            if blade_number in blade_panel_wind_stats:
                blade_mean_panel_wind, blade_mean_panel_wind_speed = blade_panel_wind_stats[
                    blade_number
                ]
                row.update({
                    f"Blade {blade_number} Mean Panel Wind X [m/s]": float(
                        blade_mean_panel_wind[0]
                    ),
                    f"Blade {blade_number} Mean Panel Wind Y [m/s]": float(
                        blade_mean_panel_wind[1]
                    ),
                    f"Blade {blade_number} Mean Panel Wind Z [m/s]": float(
                        blade_mean_panel_wind[2]
                    ),
                    f"Blade {blade_number} Mean Panel Wind Speed [m/s]": float(
                        blade_mean_panel_wind_speed
                    ),
                })

        return row

    def _build_aero_spanwise_report_rows(
        self,
        loads: RotorAeroLoads,
        *,
        time: float,
        window_index: int,
    ) -> list[dict[str, float | int | str]]:
        surface = self._runtime_context.surface
        if surface is None:
            return []

        panel_forces = loads.metadata.get("panel_forces")
        panel_collocation_points = loads.metadata.get("panel_collocation_points")
        panel_wind_velocities = loads.metadata.get("panel_wind_velocities")
        if (
            panel_forces is None
            or panel_collocation_points is None
            or panel_wind_velocities is None
        ):
            return []

        panel_forces_arr = np.asarray(panel_forces, dtype=float).reshape(-1, 3)
        panel_collocation_arr = np.asarray(panel_collocation_points, dtype=float).reshape(-1, 3)
        panel_wind_arr = np.asarray(panel_wind_velocities, dtype=float).reshape(-1, 3)

        blade_panel_slices = self._iter_blade_panel_slices()
        if not blade_panel_slices:
            return []

        total_panel_count = blade_panel_slices[-1][2].stop
        if (
            panel_forces_arr.shape[0] != total_panel_count
            or panel_collocation_arr.shape[0] != total_panel_count
            or panel_wind_arr.shape[0] != total_panel_count
        ):
            return []

        rows: list[dict[str, float | int | str]] = []
        for blade_number, blade_surface, blade_slice in blade_panel_slices:
            blade_azimuth = float(
                (
                    self._azimuth
                    + self._runtime_context.geometry.get_blade(blade_number - 1).azimuth_deg
                )
                % 360.0
            )
            blade_tip_radius = float(np.max(blade_surface.panel_reference_radii))
            blade_span_min = float(np.min(blade_surface.panel_spanwise_coordinates))
            blade_span_max = float(np.max(blade_surface.panel_spanwise_coordinates))
            blade_forces = panel_forces_arr[blade_slice]
            blade_collocation = panel_collocation_arr[blade_slice]
            blade_wind = panel_wind_arr[blade_slice]

            for panel_number, (
                force,
                collocation,
                wind_velocity,
                radius,
                spanwise,
                chordwise,
                area,
                chord_length,
                span_length,
            ) in enumerate(
                zip(
                    blade_forces,
                    blade_collocation,
                    blade_wind,
                    blade_surface.panel_reference_radii,
                    blade_surface.panel_spanwise_coordinates,
                    blade_surface.panel_chordwise_coordinates,
                    blade_surface.panel_areas,
                    blade_surface.panel_chord_lengths,
                    blade_surface.panel_span_lengths,
                ),
                start=1,
            ):
                rows.append({
                    "Schema Version": AERO_SPANWISE_REPORT_SCHEMA_VERSION,
                    "Report Type": "aero_spanwise_report",
                    "Time [s]": float(time),
                    "Window": int(window_index),
                    "Backend": str(self._runtime_context.backend),
                    "Omega [rad/s]": float(self._current_omega),
                    "Azimuth [deg]": float(self._azimuth),
                    "Blade": int(blade_number),
                    "Blade Azimuth [deg]": blade_azimuth,
                    "Panel In Blade": int(panel_number),
                    "Radius [m]": float(radius),
                    "Radius Fraction [-]": _normalize_radius_fraction(
                        radius,
                        tip_radius=blade_tip_radius,
                    ),
                    "Spanwise [m]": float(spanwise),
                    "Blade Span Fraction [-]": _normalize_fraction(
                        spanwise,
                        lower=blade_span_min,
                        upper=blade_span_max,
                    ),
                    "Chordwise [m]": float(chordwise),
                    "Area [m2]": float(area),
                    "Chord Length [m]": float(chord_length),
                    "Span Length [m]": float(span_length),
                    "Collocation X [m]": float(collocation[0]),
                    "Collocation Y [m]": float(collocation[1]),
                    "Collocation Z [m]": float(collocation[2]),
                    "Fx [N]": float(force[0]),
                    "Fy [N]": float(force[1]),
                    "Fz [N]": float(force[2]),
                    "Force Mag [N]": float(np.linalg.norm(force)),
                    "Wind X [m/s]": float(wind_velocity[0]),
                    "Wind Y [m/s]": float(wind_velocity[1]),
                    "Wind Z [m/s]": float(wind_velocity[2]),
                    "Wind Speed [m/s]": float(np.linalg.norm(wind_velocity)),
                    "Wake Rows": int(loads.metadata.get("wake_row_count", 0)),
                    "Wake Corrector Iterations": int(
                        loads.metadata.get("wake_corrector_iterations_used", 0)
                    ),
                })

        return rows

    def _build_aero_sectional_report_rows(
        self,
        loads: RotorAeroLoads,
        *,
        time: float,
        window_index: int,
    ) -> list[dict[str, float | int | str]]:
        surface = self._runtime_context.surface
        if surface is None:
            return []

        panel_forces = loads.metadata.get("panel_forces")
        panel_wind_velocities = loads.metadata.get("panel_wind_velocities")
        if panel_forces is None or panel_wind_velocities is None:
            return []

        panel_forces_arr = np.asarray(panel_forces, dtype=float).reshape(-1, 3)
        panel_wind_arr = np.asarray(panel_wind_velocities, dtype=float).reshape(-1, 3)
        blade_panel_slices = self._iter_blade_panel_slices()
        if not blade_panel_slices:
            return []

        total_panel_count = blade_panel_slices[-1][2].stop
        if (
            panel_forces_arr.shape[0] != total_panel_count
            or panel_wind_arr.shape[0] != total_panel_count
        ):
            return []

        rows: list[dict[str, float | int | str]] = []
        for blade_number, blade_surface, blade_slice in blade_panel_slices:
            blade_forces = panel_forces_arr[blade_slice]
            blade_wind = panel_wind_arr[blade_slice]
            spanwise = np.asarray(blade_surface.panel_spanwise_coordinates, dtype=float)
            radii = np.asarray(blade_surface.panel_reference_radii, dtype=float)
            chordwise = np.asarray(blade_surface.panel_chordwise_coordinates, dtype=float)
            areas = np.asarray(blade_surface.panel_areas, dtype=float)
            chord_lengths = np.asarray(blade_surface.panel_chord_lengths, dtype=float)
            span_lengths = np.asarray(blade_surface.panel_span_lengths, dtype=float)

            if spanwise.size == 0:
                continue

            span_min = float(np.min(spanwise))
            span_max = float(np.max(spanwise))
            radius_min = float(np.min(radii))
            radius_max = float(np.max(radii))
            blade_tip_radius = radius_max

            if np.isclose(span_min, span_max):
                span_edges = np.array([span_min, span_max], dtype=float)
                radius_edges = np.array([radius_min, radius_max], dtype=float)
                bin_indices = np.zeros(spanwise.shape[0], dtype=int)
                active_bin_count = 1
            else:
                span_edges = np.linspace(span_min, span_max, self._sectional_bins + 1)
                radius_edges = np.linspace(radius_min, radius_max, self._sectional_bins + 1)
                bin_indices = np.clip(
                    np.searchsorted(span_edges[1:], spanwise, side="right"),
                    0,
                    self._sectional_bins - 1,
                )
                active_bin_count = self._sectional_bins

            blade_azimuth = float(
                (
                    self._azimuth
                    + self._runtime_context.geometry.get_blade(blade_number - 1).azimuth_deg
                )
                % 360.0
            )

            for bin_number in range(active_bin_count):
                mask = bin_indices == bin_number
                if not np.any(mask):
                    continue

                bin_force = blade_forces[mask].sum(axis=0)
                bin_wind = blade_wind[mask]
                bin_area = float(np.sum(areas[mask]))
                rows.append({
                    "Schema Version": AERO_SECTIONAL_REPORT_SCHEMA_VERSION,
                    "Report Type": "aero_sectional_report",
                    "Time [s]": float(time),
                    "Window": int(window_index),
                    "Backend": str(self._runtime_context.backend),
                    "Omega [rad/s]": float(self._current_omega),
                    "Azimuth [deg]": float(self._azimuth),
                    "Blade": int(blade_number),
                    "Blade Azimuth [deg]": blade_azimuth,
                    "Section Bin": int(bin_number + 1),
                    "Panel Count": int(np.count_nonzero(mask)),
                    "Spanwise Start [m]": float(span_edges[bin_number]),
                    "Spanwise End [m]": float(span_edges[bin_number + 1]),
                    "Spanwise Center [m]": float(np.mean(spanwise[mask])),
                    "Blade Span Start Fraction [-]": _normalize_fraction(
                        span_edges[bin_number],
                        lower=span_min,
                        upper=span_max,
                    ),
                    "Blade Span End Fraction [-]": _normalize_fraction(
                        span_edges[bin_number + 1],
                        lower=span_min,
                        upper=span_max,
                    ),
                    "Blade Span Center Fraction [-]": _normalize_fraction(
                        float(np.mean(spanwise[mask])),
                        lower=span_min,
                        upper=span_max,
                    ),
                    "Radius Start [m]": float(radius_edges[bin_number]),
                    "Radius End [m]": float(radius_edges[bin_number + 1]),
                    "Radius Center [m]": float(np.mean(radii[mask])),
                    "Radius Start Fraction [-]": _normalize_radius_fraction(
                        radius_edges[bin_number],
                        tip_radius=blade_tip_radius,
                    ),
                    "Radius End Fraction [-]": _normalize_radius_fraction(
                        radius_edges[bin_number + 1],
                        tip_radius=blade_tip_radius,
                    ),
                    "Radius Center Fraction [-]": _normalize_radius_fraction(
                        float(np.mean(radii[mask])),
                        tip_radius=blade_tip_radius,
                    ),
                    "Area [m2]": bin_area,
                    "Mean Chord Length [m]": float(np.mean(chord_lengths[mask])),
                    "Mean Span Length [m]": float(np.mean(span_lengths[mask])),
                    "Mean Chordwise [m]": float(np.mean(chordwise[mask])),
                    "Fx [N]": float(bin_force[0]),
                    "Fy [N]": float(bin_force[1]),
                    "Fz [N]": float(bin_force[2]),
                    "Force Mag [N]": float(np.linalg.norm(bin_force)),
                    "Wind X [m/s]": float(np.mean(bin_wind[:, 0])),
                    "Wind Y [m/s]": float(np.mean(bin_wind[:, 1])),
                    "Wind Z [m/s]": float(np.mean(bin_wind[:, 2])),
                    "Wind Speed [m/s]": float(np.mean(np.linalg.norm(bin_wind, axis=1))),
                    "Wake Rows": int(loads.metadata.get("wake_row_count", 0)),
                    "Wake Corrector Iterations": int(
                        loads.metadata.get("wake_corrector_iterations_used", 0)
                    ),
                })

        return rows

    def _append_aero_report(
        self,
        loads: RotorAeroLoads,
        *,
        displacements: np.ndarray,
        velocities: np.ndarray | None,
        time: float,
        window_index: int,
    ) -> None:
        if self._output_folder is None:
            return

        self._ensure_output_schema_manifest()
        csv_path = self._output_folder / "aero_report.csv"
        write_header = not csv_path.exists()
        row = self._build_aero_report_row(
            loads,
            displacements=displacements,
            velocities=velocities,
            time=time,
            window_index=window_index,
        )

        with csv_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self._aero_report_fieldnames())
            if write_header:
                writer.writeheader()
            writer.writerow(row)

    def _append_aero_spanwise_report(
        self,
        loads: RotorAeroLoads,
        *,
        time: float,
        window_index: int,
    ) -> None:
        if self._output_folder is None:
            return

        self._ensure_output_schema_manifest()

        rows = self._build_aero_spanwise_report_rows(
            loads,
            time=time,
            window_index=window_index,
        )
        if not rows:
            return

        csv_path = self._output_folder / "aero_spanwise_report.csv"
        write_header = not csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(AERO_SPANWISE_REPORT_FIELDNAMES))
            if write_header:
                writer.writeheader()
            writer.writerows(rows)

    def _append_aero_sectional_report(
        self,
        loads: RotorAeroLoads,
        *,
        time: float,
        window_index: int,
    ) -> None:
        if self._output_folder is None:
            return

        self._ensure_output_schema_manifest()

        rows = self._build_aero_sectional_report_rows(
            loads,
            time=time,
            window_index=window_index,
        )
        if not rows:
            return

        csv_path = self._output_folder / "aero_sectional_report.csv"
        write_header = not csv_path.exists()
        with csv_path.open("a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(AERO_SECTIONAL_REPORT_FIELDNAMES))
            if write_header:
                writer.writeheader()
            writer.writerows(rows)

    def _log_window_summary(self, loads: RotorAeroLoads, *, time: float, window_index: int) -> None:
        if window_index % self._log_interval != 0:
            return

        integrated_force = (
            np.asarray(loads.integrated_force, dtype=float).reshape(3)
            if loads.integrated_force is not None
            else np.sum(loads.nodal_forces, axis=0)
        )
        logger.info(
            "[AERO-FSI] Window=%d time=%.6f force=[%.3e, %.3e, %.3e] wake_rows=%d",
            window_index,
            time,
            integrated_force[0],
            integrated_force[1],
            integrated_force[2],
            int(loads.metadata.get("wake_row_count", 0)),
        )

    def run(self) -> None:
        """Execute the generic preCICE coupling loop."""
        from aeroelast.solvers.fsi.base import Adapter

        coupling_meshes = {self._coupling_mesh: self._ref_coords}
        if self._omega_mesh is not None:
            coupling_meshes[self._omega_mesh] = self._omega_vertex.reshape(1, 3)

        logger.info(
            "[AERO-FSI] Preparing preCICE adapter: participant=%s config=%s",
            self._participant_name,
            self._config_file,
        )
        adapter = Adapter(
            participant=self._participant_name,
            config_file=self._config_file,
            coupling_meshes=coupling_meshes,
        )
        logger.info("[AERO-FSI] Initializing preCICE handshake...")
        adapter.initialize()

        logger.info(
            "[AERO-FSI] Initialized. Participant=%s mesh=%s nodes=%d backend=%s",
            self._participant_name,
            self._coupling_mesh,
            self._n_nodes,
            self._backend.name,
        )
        logger.info("[AERO-FSI] Entered coupling loop.")

        while adapter.is_coupling_ongoing:
            if adapter.requires_writing_checkpoint:
                adapter.store_checkpoint((
                    self._last_forces.copy(),
                    self._current_omega,
                    self._azimuth,
                    self._current_time,
                    self._capture_backend_checkpoint_state(),
                ))

            if adapter.requires_reading_checkpoint:
                states = adapter.retrieve_checkpoint()
                self._last_forces = states[0]
                self._current_omega = float(states[1])
                self._azimuth = float(states[2])
                if len(states) > 4:
                    self._current_time = float(states[3])
                    backend_state = states[4]
                else:
                    backend_state = states[3] if len(states) > 3 else None
                self._restore_backend_checkpoint_state(backend_state)
                self._iteration_count += 1

            dt = adapter.dt
            current_time = self._current_time

            displacements = self._coerce_nodal_field(
                adapter.read_data(self._coupling_mesh, self._displacement_data),
                name=self._displacement_data,
            )
            velocities = self._read_optional_nodal_data(adapter, self._velocity_data)
            self._update_live_omega(adapter)

            loads = self.compute_iteration_loads(
                displacements,
                velocities=velocities,
                current_time=current_time,
                dt=dt,
            )
            adapter.write_data(self._coupling_mesh, self._force_data, loads.nodal_forces)
            adapter.advance(dt)

            if adapter.is_time_window_complete:
                self._finalize_backend_time_window()
                self._last_forces = loads.nodal_forces.copy()
                completed_window = self._window_count + 1
                completed_time = current_time + float(dt)
                self._append_aero_report(
                    loads,
                    displacements=displacements,
                    velocities=velocities,
                    time=completed_time,
                    window_index=completed_window,
                )
                self._append_aero_spanwise_report(
                    loads,
                    time=completed_time,
                    window_index=completed_window,
                )
                self._append_aero_sectional_report(
                    loads,
                    time=completed_time,
                    window_index=completed_window,
                )
                self._log_window_summary(
                    loads,
                    time=completed_time,
                    window_index=completed_window,
                )
                self._advance_time_window(dt)

        logger.info("[AERO-FSI] Finalizing preCICE adapter...")
        adapter.finalize()


def build_rotor_aero_participant(
    runtime_context: AeroFSIRuntimeContext,
    backend: RotorAerodynamicBackend,
    *,
    config_file: str = "precice-config.xml",
) -> RotorAeroFSIParticipant:
    """Build the generic full-rotor participant from normalized config defaults."""
    normalized_cfg = runtime_context.normalized_config
    output_cfg_raw = normalized_cfg.get("output")
    output_cfg = dict(output_cfg_raw) if isinstance(output_cfg_raw, Mapping) else {}
    omega_initial, omega_mesh, omega_data, omega_vertex = _extract_initial_omega_cfg(
        runtime_context
    )
    wind_velocity_provider = _build_wind_velocity_provider(runtime_context)
    initial_wind_velocity = wind_velocity_provider.evaluate(0.0)

    return RotorAeroFSIParticipant(
        runtime_context=runtime_context,
        backend=backend,
        participant=str(normalized_cfg.get("participant", runtime_context.participant_name)),
        config_file=str(normalized_cfg.get("config_file", config_file)),
        coupling_mesh=str(normalized_cfg.get("coupling_mesh", runtime_context.coupling_mesh_name)),
        displacement_data=str(normalized_cfg.get("displacement_data", "Displacement")),
        force_data=str(normalized_cfg.get("force_data", "Force")),
        velocity_data=normalized_cfg.get("velocity_data", "Velocity"),
        omega_mesh=omega_mesh,
        omega_data=omega_data,
        omega_vertex=omega_vertex,
        omega_rad_s=omega_initial,
        azimuth_deg=float(_extract_backend_cfg(runtime_context).get("azimuth", 0.0)),
        wind_velocity=initial_wind_velocity,
        wind_velocity_provider=wind_velocity_provider,
        output_folder=output_cfg.get("folder", "aero_fsi_results"),
        aero_mesh_file=output_cfg.get("aero_mesh_file"),
        log_interval=int(output_cfg.get("log_interval", 10)),
        sectional_bins=int(output_cfg.get("sectional_bins", 8)),
    )
