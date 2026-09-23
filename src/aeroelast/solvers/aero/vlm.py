"""Quasi-steady vortex-lattice backend for full-rotor aerodynamic participants."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .surface import FullRotorAeroSurface, FullRotorAeroSurfaceState
from .types import AeroFSIRuntimeContext, RotorAeroLoads, RotorAeroState


def _normalize_vector(vector: np.ndarray, *, name: str) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(3)
    norm = np.linalg.norm(arr)
    if np.isclose(norm, 0.0):
        raise ValueError(f"{name} cannot be the zero vector")
    return arr / norm


@dataclass
class _WakeRow:
    points: np.ndarray
    convection_velocities: np.ndarray


def _copy_wake_rows(rows: list[_WakeRow]) -> list[_WakeRow]:
    return [
        _WakeRow(
            points=row.points.copy(),
            convection_velocities=row.convection_velocities.copy(),
        )
        for row in rows
    ]


def _segment_induced_velocity(
    point: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    core_radius: float,
) -> np.ndarray:
    r1 = point - start
    r2 = point - end
    r0 = end - start

    cross = np.cross(r1, r2)
    cross_sq = float(np.dot(cross, cross) + core_radius**2)
    r1_norm = float(np.linalg.norm(r1))
    r2_norm = float(np.linalg.norm(r2))
    if r1_norm < 1e-12 or r2_norm < 1e-12:
        return np.zeros(3, dtype=float)

    coeff = np.dot(r0, r1 / r1_norm - r2 / r2_norm) / (4.0 * np.pi * cross_sq)
    return coeff * cross


def _horseshoe_induced_velocity(
    point: np.ndarray,
    bound_root: np.ndarray,
    bound_tip: np.ndarray,
    wake_direction: np.ndarray,
    wake_length: float,
    core_radius: float,
) -> np.ndarray:
    root_far = bound_root + wake_length * wake_direction
    tip_far = bound_tip + wake_length * wake_direction
    return (
        _segment_induced_velocity(point, bound_root, bound_tip, core_radius)
        + _segment_induced_velocity(point, bound_tip, tip_far, core_radius)
        + _segment_induced_velocity(point, tip_far, root_far, core_radius)
        + _segment_induced_velocity(point, root_far, bound_root, core_radius)
    )


def _panel_mean_velocity(
    nodal_velocities: np.ndarray | None,
    parent_node_indices: np.ndarray,
) -> np.ndarray:
    if nodal_velocities is None:
        return np.zeros(3, dtype=float)
    return nodal_velocities[parent_node_indices].mean(axis=0)


def _evaluate_wind_velocity_points(
    state: RotorAeroState,
    points: np.ndarray,
) -> np.ndarray:
    points_arr = np.asarray(points, dtype=float).reshape(-1, 3)
    provider = state.extra.get("wind_velocity_provider")
    if provider is None or not hasattr(provider, "evaluate_at_points"):
        return np.broadcast_to(state.wind_velocity, (points_arr.shape[0], 3)).copy()

    velocities = np.asarray(provider.evaluate_at_points(points_arr, state.time), dtype=float)
    if velocities.shape != (points_arr.shape[0], 3):
        raise ValueError(
            "wind_velocity_provider.evaluate_at_points must return shape (n_points, 3)"
        )
    return velocities


def _select_surface_state(
    context: AeroFSIRuntimeContext,
    state: RotorAeroState,
) -> tuple[FullRotorAeroSurface, FullRotorAeroSurfaceState]:
    if context.surface is None:
        raise ValueError("VLM backend requires a runtime context with aerodynamic surface data")

    deformed_surface = state.extra.get("deformed_surface")
    if isinstance(deformed_surface, FullRotorAeroSurfaceState):
        return context.surface, deformed_surface

    if state.nodal_displacements is not None:
        return context.surface, context.surface.build_deformed_state(state.nodal_displacements)

    zero_displacements = np.zeros((context.surface.total_parent_nodes, 3), dtype=float)
    return context.surface, context.surface.build_deformed_state(zero_displacements)


@dataclass
class _PanelData:
    blade_index: int
    panel_index: int
    collocation_point: np.ndarray
    normal: np.ndarray
    centroid: np.ndarray
    bound_root: np.ndarray
    bound_tip: np.ndarray
    bound_vector: np.ndarray
    area: float
    parent_node_indices: np.ndarray
    chord_length: float
    wake_direction: np.ndarray
    collocation_wind_velocity: np.ndarray
    bound_wind_velocities: np.ndarray
    material_velocity: np.ndarray
    relative_flow: np.ndarray


def _build_shed_wake_row(
    panels: list[_PanelData],
    dt: float,
) -> _WakeRow | None:
    if dt <= 0.0 or not panels:
        return None

    convection_velocities = np.asarray(
        [panel.bound_wind_velocities for panel in panels],
        dtype=float,
    )
    if np.allclose(convection_velocities, 0.0):
        return None

    step = dt * convection_velocities
    if np.allclose(step, 0.0):
        return None

    points = (
        np.asarray(
            [[panel.bound_root, panel.bound_tip] for panel in panels],
            dtype=float,
        )
        + step
    )
    return _WakeRow(points=points, convection_velocities=convection_velocities)


def _advect_wake_rows(rows: list[_WakeRow], dt: float) -> list[_WakeRow]:
    advected_rows: list[_WakeRow] = []
    for row in rows:
        points = row.points.copy()
        if dt > 0.0:
            points += dt * row.convection_velocities
        advected_rows.append(
            _WakeRow(
                points=points,
                convection_velocities=row.convection_velocities.copy(),
            )
        )
    return advected_rows


def _total_induced_velocity(
    point: np.ndarray,
    panels: list[_PanelData],
    circulation: np.ndarray,
    wake_rows: list[_WakeRow],
    wake_length: float,
    core_radius: float,
) -> np.ndarray:
    velocity = np.zeros(3, dtype=float)
    for panel_index, (gamma, panel) in enumerate(zip(circulation, panels)):
        velocity += gamma * _wake_induced_velocity(
            point,
            panel,
            panel_index,
            wake_rows,
            wake_length,
            core_radius,
        )
    return velocity


def _refresh_wake_rows_convection_velocities(
    rows: list[_WakeRow],
    panels: list[_PanelData],
    circulation: np.ndarray,
    state: RotorAeroState,
    wake_length: float,
    core_radius: float,
) -> list[_WakeRow]:
    refreshed_rows: list[_WakeRow] = []

    for row in rows:
        free_stream = _evaluate_wind_velocity_points(state, row.points.reshape(-1, 3)).reshape(
            row.points.shape
        )
        convection_velocities = np.zeros_like(row.points)
        for panel_index in range(row.points.shape[0]):
            for edge_index in range(row.points.shape[1]):
                point = row.points[panel_index, edge_index]
                convection_velocities[panel_index, edge_index] = free_stream[
                    panel_index, edge_index
                ] + _total_induced_velocity(
                    point,
                    panels,
                    circulation,
                    rows,
                    wake_length,
                    core_radius,
                )

        refreshed_rows.append(
            _WakeRow(
                points=row.points.copy(),
                convection_velocities=convection_velocities,
            )
        )

    return refreshed_rows


def _correct_wake_rows_geometry(
    predictor_rows: list[_WakeRow],
    refreshed_rows: list[_WakeRow],
    dt: float,
) -> list[_WakeRow]:
    if dt <= 0.0 or not predictor_rows:
        return _copy_wake_rows(refreshed_rows)

    corrected_rows: list[_WakeRow] = []
    for predictor_row, refreshed_row in zip(predictor_rows, refreshed_rows):
        corrected_points = predictor_row.points + 0.5 * dt * (
            refreshed_row.convection_velocities - predictor_row.convection_velocities
        )
        corrected_rows.append(
            _WakeRow(
                points=corrected_points,
                convection_velocities=refreshed_row.convection_velocities.copy(),
            )
        )

    return corrected_rows


def _wake_induced_velocity(
    point: np.ndarray,
    panel: _PanelData,
    panel_index: int,
    wake_rows: list[_WakeRow],
    wake_length: float,
    core_radius: float,
) -> np.ndarray:
    if not wake_rows:
        return _horseshoe_induced_velocity(
            point,
            panel.bound_root,
            panel.bound_tip,
            panel.wake_direction,
            wake_length,
            core_radius,
        )

    velocity = _segment_induced_velocity(point, panel.bound_root, panel.bound_tip, core_radius)
    previous_tip = panel.bound_tip
    previous_root = panel.bound_root

    for row in wake_rows:
        tip_point = row.points[panel_index, 1]
        root_point = row.points[panel_index, 0]
        velocity += _segment_induced_velocity(point, previous_tip, tip_point, core_radius)
        velocity += _segment_induced_velocity(point, root_point, previous_root, core_radius)
        previous_tip = tip_point
        previous_root = root_point

    velocity += _segment_induced_velocity(point, previous_tip, previous_root, core_radius)
    return velocity


def _flatten_surface_panels(
    surface_state: FullRotorAeroSurfaceState,
    state: RotorAeroState,
) -> list[_PanelData]:
    omega_vector = state.omega_rad_s * surface_state.rotation_axis
    panels: list[_PanelData] = []

    for blade_state in surface_state.blades:
        for panel_index in range(blade_state.n_panels):
            collocation = blade_state.panel_collocation_points[panel_index]
            bound_points = blade_state.panel_bound_vortex_points[panel_index]
            collocation_wind_velocity = _evaluate_wind_velocity_points(
                state,
                collocation.reshape(1, 3),
            )[0]
            bound_wind_velocities = _evaluate_wind_velocity_points(state, bound_points)
            bound_root = bound_points[0]
            bound_tip = bound_points[1]
            bound_vector = bound_tip - bound_root
            rotation_velocity = np.cross(omega_vector, collocation - surface_state.rotation_center)
            elastic_velocity = _panel_mean_velocity(
                state.nodal_velocities,
                blade_state.panel_parent_node_indices[panel_index],
            )
            material_velocity = rotation_velocity + elastic_velocity
            relative_flow = collocation_wind_velocity - material_velocity
            flow_norm = np.linalg.norm(relative_flow)
            if np.isclose(flow_norm, 0.0):
                wake_direction = blade_state.chord_direction.copy()
            else:
                wake_direction = relative_flow / flow_norm

            panels.append(
                _PanelData(
                    blade_index=blade_state.blade_index,
                    panel_index=panel_index,
                    collocation_point=collocation,
                    normal=blade_state.panel_normals[panel_index],
                    centroid=blade_state.panel_centroids[panel_index],
                    bound_root=bound_root,
                    bound_tip=bound_tip,
                    bound_vector=bound_vector,
                    area=float(blade_state.panel_areas[panel_index]),
                    parent_node_indices=blade_state.panel_parent_node_indices[panel_index],
                    chord_length=float(blade_state.panel_chord_lengths[panel_index]),
                    wake_direction=wake_direction,
                    collocation_wind_velocity=collocation_wind_velocity,
                    bound_wind_velocities=bound_wind_velocities,
                    material_velocity=material_velocity,
                    relative_flow=relative_flow,
                )
            )

    return panels


class QuasiSteadyVLMBackend:
    """Minimal horseshoe-vortex backend built on the reduced aerodynamic surface."""

    name = "vlm"

    def __init__(
        self,
        *,
        air_density: float = 1.225,
        wake_length: float | None = None,
        wake_length_scale: float = 25.0,
        wake_history_steps: int = 0,
        wake_corrector_iterations: int = 0,
        diagonal_regularization: float = 1e-8,
        core_radius: float = 1e-8,
    ) -> None:
        self._air_density = float(air_density)
        self._wake_length = None if wake_length is None else float(wake_length)
        self._wake_length_scale = float(wake_length_scale)
        self._wake_history_steps = max(0, int(wake_history_steps))
        self._wake_corrector_iterations = max(0, int(wake_corrector_iterations))
        self._diagonal_regularization = float(diagonal_regularization)
        self._core_radius = float(core_radius)
        self._committed_wake_rows: list[_WakeRow] = []
        self._last_working_wake_rows: list[_WakeRow] = []

    def capture_state(self) -> dict[str, Any]:
        return {
            "committed_wake_rows": _copy_wake_rows(self._committed_wake_rows),
            "last_working_wake_rows": _copy_wake_rows(self._last_working_wake_rows),
        }

    def restore_state(self, checkpoint_state: Any) -> None:
        if not isinstance(checkpoint_state, dict):
            self._committed_wake_rows = []
            self._last_working_wake_rows = []
            return

        committed_rows = checkpoint_state.get("committed_wake_rows", [])
        self._committed_wake_rows = _copy_wake_rows(list(committed_rows))

        last_working_rows = checkpoint_state.get("last_working_wake_rows")
        if isinstance(last_working_rows, list):
            self._last_working_wake_rows = _copy_wake_rows(last_working_rows)
        else:
            self._last_working_wake_rows = _copy_wake_rows(self._committed_wake_rows)

    def finalize_time_window(self) -> None:
        if self._wake_history_steps <= 0:
            self._committed_wake_rows = []
            self._last_working_wake_rows = []
            return

        self._committed_wake_rows = _copy_wake_rows(
            self._last_working_wake_rows[: self._wake_history_steps]
        )

    def _build_working_wake_rows(
        self,
        panels: list[_PanelData],
        dt: float,
    ) -> list[_WakeRow]:
        if self._wake_history_steps <= 0:
            return []

        working_rows = _advect_wake_rows(self._committed_wake_rows, dt)
        shed_row = _build_shed_wake_row(panels, dt)
        if shed_row is not None:
            working_rows.insert(0, shed_row)

        return working_rows[: self._wake_history_steps]

    def _resolve_wake_length(self, panels: list[_PanelData]) -> float:
        if self._wake_length is not None:
            return self._wake_length
        mean_chord = max(1.0e-6, float(np.mean([panel.chord_length for panel in panels])))
        return max(mean_chord * self._wake_length_scale, mean_chord)

    def _assemble_influence_matrix(
        self,
        panels: list[_PanelData],
        wake_rows: list[_WakeRow],
        wake_length: float,
    ) -> np.ndarray:
        panel_count = len(panels)
        influence = np.zeros((panel_count, panel_count), dtype=float)
        for i, target in enumerate(panels):
            for j, source in enumerate(panels):
                induced = _wake_induced_velocity(
                    target.collocation_point,
                    source,
                    j,
                    wake_rows,
                    wake_length,
                    self._core_radius,
                )
                influence[i, j] = np.dot(induced, target.normal)

        if self._diagonal_regularization > 0.0:
            influence += self._diagonal_regularization * np.eye(panel_count)
        return influence

    def _solve_circulation(self, influence: np.ndarray, panels: list[_PanelData]) -> np.ndarray:
        rhs = -np.array(
            [np.dot(panel.relative_flow, panel.normal) for panel in panels], dtype=float
        )
        try:
            return np.linalg.solve(influence, rhs)
        except np.linalg.LinAlgError:
            return np.linalg.lstsq(influence, rhs, rcond=None)[0]

    def _induced_velocity_on_panels(
        self,
        panels: list[_PanelData],
        circulation: np.ndarray,
        wake_rows: list[_WakeRow],
        wake_length: float,
    ) -> np.ndarray:
        induced = np.zeros((len(panels), 3), dtype=float)
        for i, target in enumerate(panels):
            velocity = np.zeros(3, dtype=float)
            for j, (gamma, source) in enumerate(zip(circulation, panels)):
                velocity += gamma * _wake_induced_velocity(
                    target.collocation_point,
                    source,
                    j,
                    wake_rows,
                    wake_length,
                    self._core_radius,
                )
            induced[i] = velocity
        return induced

    def _panel_forces(
        self,
        panels: list[_PanelData],
        circulation: np.ndarray,
        induced_velocity: np.ndarray,
    ) -> np.ndarray:
        panel_forces = np.zeros((len(panels), 3), dtype=float)
        for index, panel in enumerate(panels):
            total_velocity = panel.relative_flow + induced_velocity[index]
            raw_force = (
                self._air_density
                * circulation[index]
                * np.cross(total_velocity, panel.bound_vector)
            )
            normal_component = np.dot(raw_force, panel.normal)
            panel_forces[index] = normal_component * panel.normal
        return panel_forces

    def compute_loads(
        self,
        context: AeroFSIRuntimeContext,
        state: RotorAeroState,
    ) -> RotorAeroLoads:
        reference_surface, surface_state = _select_surface_state(context, state)
        panels = _flatten_surface_panels(surface_state, state)
        if not panels:
            nodal_forces = np.zeros((reference_surface.total_parent_nodes, 3), dtype=float)
            return RotorAeroLoads(nodal_forces=nodal_forces)

        wake_length = self._resolve_wake_length(panels)
        wake_rows = self._build_working_wake_rows(
            panels,
            state.dt,
        )
        predictor_wake_rows = _copy_wake_rows(wake_rows)
        corrector_iterations_used = 0

        for iteration in range(self._wake_corrector_iterations + 1):
            influence = self._assemble_influence_matrix(panels, wake_rows, wake_length)
            circulation = self._solve_circulation(influence, panels)
            induced_velocity = self._induced_velocity_on_panels(
                panels,
                circulation,
                wake_rows,
                wake_length,
            )
            panel_forces = self._panel_forces(panels, circulation, induced_velocity)
            refreshed_wake_rows = _refresh_wake_rows_convection_velocities(
                wake_rows,
                panels,
                circulation,
                state,
                wake_length,
                self._core_radius,
            )
            if iteration >= self._wake_corrector_iterations:
                wake_rows = refreshed_wake_rows
                break

            wake_rows = _correct_wake_rows_geometry(
                wake_rows,
                refreshed_wake_rows,
                state.dt,
            )
            corrector_iterations_used = iteration + 1

        self._last_working_wake_rows = _copy_wake_rows(wake_rows)

        blade_panel_forces: list[np.ndarray] = []
        cursor = 0
        for blade in reference_surface.blades:
            count = blade.n_panels
            blade_panel_forces.append(panel_forces[cursor : cursor + count])
            cursor += count

        nodal_forces = reference_surface.assemble_global_nodal_forces(blade_panel_forces)
        integrated_force = panel_forces.sum(axis=0)
        integrated_moment = np.zeros(3, dtype=float)
        for panel, force in zip(panels, panel_forces):
            integrated_moment += np.cross(panel.centroid - surface_state.rotation_center, force)

        metadata: dict[str, Any] = {
            "circulation": circulation,
            "panel_forces": panel_forces,
            "panel_relative_flows": np.asarray(
                [panel.relative_flow for panel in panels],
                dtype=float,
            ),
            "panel_induced_velocities": induced_velocity,
            "panel_incident_velocities": np.asarray(
                [panel.relative_flow for panel in panels],
                dtype=float,
            )
            + induced_velocity
            - np.asarray(
                [
                    circulation[index]
                    * _wake_induced_velocity(
                        panel.collocation_point,
                        panel,
                        index,
                        wake_rows,
                        wake_length,
                        self._core_radius,
                    )
                    for index, panel in enumerate(panels)
                ],
                dtype=float,
            ),
            "panel_total_velocities": np.asarray(
                [panel.relative_flow for panel in panels],
                dtype=float,
            )
            + induced_velocity,
            "panel_collocation_points": np.asarray([panel.collocation_point for panel in panels]),
            "panel_wind_velocities": np.asarray(
                [panel.collocation_wind_velocity for panel in panels],
                dtype=float,
            ),
            "wake_length": wake_length,
            "wake_row_count": len(wake_rows),
            "wake_predictor_rows": np.asarray(
                [row.points for row in predictor_wake_rows], dtype=float
            ),
            "wake_rows": np.asarray([row.points for row in wake_rows], dtype=float),
            "wake_convection_velocities": np.asarray(
                [row.convection_velocities for row in wake_rows],
                dtype=float,
            ),
            "wake_corrector_iterations_used": corrector_iterations_used,
        }
        return RotorAeroLoads(
            nodal_forces=nodal_forces,
            integrated_force=integrated_force,
            integrated_moment=integrated_moment,
            metadata=metadata,
        )


def build_vlm_backend(runtime_context: AeroFSIRuntimeContext) -> QuasiSteadyVLMBackend:
    """Build the quasi-steady VLM backend from the normalized aerodynamic config."""
    aero_cfg = runtime_context.normalized_config.get("aero", {})
    backend_cfg_raw = aero_cfg.get("backend_config", {}) if isinstance(aero_cfg, dict) else {}
    legacy_cfg_raw = runtime_context.normalized_config.get("vlm", {})
    backend_cfg = dict(legacy_cfg_raw) if isinstance(legacy_cfg_raw, dict) else {}
    if isinstance(backend_cfg_raw, dict):
        backend_cfg.update(backend_cfg_raw)

    return QuasiSteadyVLMBackend(
        air_density=float(backend_cfg.get("air_density", 1.225)),
        wake_length=backend_cfg.get("wake_length"),
        wake_length_scale=float(backend_cfg.get("wake_length_scale", 25.0)),
        wake_history_steps=int(backend_cfg.get("wake_history_steps", 0)),
        wake_corrector_iterations=int(backend_cfg.get("wake_corrector_iterations", 0)),
        diagonal_regularization=float(backend_cfg.get("diagonal_regularization", 1.0e-8)),
        core_radius=float(backend_cfg.get("core_radius", 1.0e-8)),
    )
