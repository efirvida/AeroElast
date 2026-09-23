"""Validation helpers for reduced-order VLM versus BEM rotor comparisons.

These utilities build a lightweight full-rotor lifting-surface mesh directly
from the blade aerodynamic definition so the current Python VLM backend can be
compared against the BEM solver across operating-point sweeps without relying
on the much heavier structural shell mesh.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Sequence

import numpy as np

from aeroelast.core.mesh import ElementSet, ElementType, MeshElement, MeshModel, Node, NodeSet
from aeroelast.models.blade.aerodynamics import BladeAero, load_blade_aero

from .participant import build_aero_runtime_context
from .types import RotorAeroState
from .vlm import _flatten_surface_panels, _wake_induced_velocity, build_vlm_backend


@dataclass(frozen=True)
class RotorOperatingPoint:
    """Single rigid-rotor operating point used for aerodynamic comparisons."""

    name: str
    wind_speed: float
    omega_rpm: float
    pitch_deg: float = 0.0


@dataclass(frozen=True)
class RotorPerformanceMetrics:
    """Integrated rotor performance metrics for one aerodynamic model."""

    thrust: float
    torque: float
    power: float
    cp: float
    ct: float
    cq: float


@dataclass(frozen=True)
class RotorModelComparison:
    """Flattened VLM versus BEM comparison at one operating point."""

    operating_point: RotorOperatingPoint
    vlm: RotorPerformanceMetrics
    bem: RotorPerformanceMetrics
    panel_count: int
    thrust_ratio: float
    torque_ratio: float
    power_ratio: float
    cp_ratio: float
    ct_ratio: float
    cq_ratio: float
    cp_delta: float
    ct_delta: float


@dataclass(frozen=True)
class RotorSectionalDiagnostics:
    """Per-strip reduced-VLM diagnostics for one rotor operating point."""

    operating_point: RotorOperatingPoint
    radii: np.ndarray
    alpha_incident_deg: np.ndarray
    alpha_total_deg: np.ndarray
    cl_polar: np.ndarray
    cd_polar: np.ndarray
    cl_vlm_trailing: np.ndarray
    gamma_leading: np.ndarray
    gamma_trailing: np.ndarray
    strip_base_thrust: np.ndarray
    strip_base_torque: np.ndarray


@dataclass(frozen=True)
class RotorValidationObjectiveMetrics:
    """Objective error metrics for reduced-order model-vs-BEM sweeps."""

    mean_abs_thrust_ratio_error: float
    mean_abs_torque_ratio_error: float
    mean_abs_power_ratio_error: float
    mean_abs_cp_error: float
    zero_pitch_mean_abs_torque_ratio_error: float
    cp_curve_rmse: float
    cp_peak_rpm_error: float
    cp_peak_value_error: float
    bem_cp_peak_rpm: float
    model_cp_peak_rpm: float
    bem_cp_peak: float
    model_cp_peak: float
    reference_sample_count: int
    tsr_sample_count: int


IEA15MW_REFERENCE_OPERATING_POINTS: tuple[RotorOperatingPoint, ...] = (
    RotorOperatingPoint("v5.006", wind_speed=5.006, omega_rpm=5.0, pitch_deg=2.893),
    RotorOperatingPoint("v7.159", wind_speed=7.159, omega_rpm=5.10, pitch_deg=0.0),
    RotorOperatingPoint("v9.027", wind_speed=9.027, omega_rpm=6.43, pitch_deg=0.0),
    RotorOperatingPoint("v10.659", wind_speed=10.659, omega_rpm=7.518, pitch_deg=0.0),
    RotorOperatingPoint("v12.259", wind_speed=12.259, omega_rpm=7.518, pitch_deg=6.76),
    RotorOperatingPoint("v15.471", wind_speed=15.471, omega_rpm=7.518, pitch_deg=12.19),
)

IEA15MW_ZERO_PITCH_TSR_SWEEP: tuple[RotorOperatingPoint, ...] = tuple(
    RotorOperatingPoint(
        name=f"rpm{rpm:04.1f}",
        wind_speed=10.659,
        omega_rpm=rpm,
        pitch_deg=0.0,
    )
    for rpm in (2.5, 4.0, 5.5, 7.0, 8.5, 10.0, 11.5)
)


PolarCorrectionLiftMode = Literal[
    "drag_only",
    "delta_cl_trailing",
    "delta_cl_trailing_positive_gated",
    "delta_cl_trailing_positive_gated_tsr_weighted",
]


def _coerce_blade_aero(blade: BladeAero | str | Path) -> tuple[BladeAero, str]:
    if isinstance(blade, BladeAero):
        return blade, "<in-memory-blade-aero>"
    blade_path = Path(blade).expanduser().resolve()
    return load_blade_aero(str(blade_path)), str(blade_path)


def _interpolate_station_field(
    blade_aero: BladeAero,
    span_grid: np.ndarray,
    values: Sequence[float],
) -> np.ndarray:
    station_span = np.array([station.span_fraction for station in blade_aero.stations], dtype=float)
    return np.interp(span_grid, station_span, np.asarray(values, dtype=float))


def build_reduced_vlm_rotor_mesh(
    blade_aero: BladeAero,
    *,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    collective_pitch_deg: float = 0.0,
) -> MeshModel:
    """Build a lightweight full-rotor mesh from aerodynamic blade stations.

    The reduced mesh is a twisted lifting surface with real spanwise chord,
    twist, pitch-axis location, prebend, sweep, hub offset, and full-rotor
    azimuth placement.
    It is intended for validation sweeps where the full shell mesh would make
    the current Python VLM backend prohibitively expensive.
    """

    if n_spanwise_panels < 1:
        raise ValueError("n_spanwise_panels must be at least 1")
    if n_chordwise_panels < 1:
        raise ValueError("n_chordwise_panels must be at least 1")

    mesh = MeshModel()
    rotation_axis = np.array([0.0, 1.0, 0.0], dtype=float)
    span_grid = np.linspace(
        blade_aero.stations[0].span_fraction,
        blade_aero.stations[-1].span_fraction,
        n_spanwise_panels + 1,
        dtype=float,
    )
    chord_grid = np.linspace(0.0, 1.0, n_chordwise_panels + 1, dtype=float)

    radii = _interpolate_station_field(blade_aero, span_grid, blade_aero.r)
    chords = _interpolate_station_field(blade_aero, span_grid, blade_aero.chord)
    twists = _interpolate_station_field(blade_aero, span_grid, blade_aero.twist)
    prebend = _interpolate_station_field(blade_aero, span_grid, blade_aero.prebend)
    sweep = _interpolate_station_field(blade_aero, span_grid, blade_aero.sweep)
    pitch_axes = _interpolate_station_field(
        blade_aero,
        span_grid,
        [station.pitch_axis for station in blade_aero.stations],
    )
    # Keep spanwise shape while anchoring root at the hub frame.
    prebend = prebend - prebend[0]
    sweep = sweep - sweep[0]
    collective_pitch_rad = math.radians(collective_pitch_deg)
    twists = twists + collective_pitch_rad

    for blade_index in range(blade_aero.n_blades):
        azimuth = blade_index * 2.0 * math.pi / blade_aero.n_blades
        radial_direction = np.array([math.sin(azimuth), 0.0, math.cos(azimuth)], dtype=float)
        tangential_direction = np.cross(rotation_axis, radial_direction)

        blade_station_nodes: list[list[Node]] = []
        for radius, chord, twist, pitch_axis, prebend_i, sweep_i in zip(
            radii,
            chords,
            twists,
            pitch_axes,
            prebend,
            sweep,
        ):
            chord_direction = (
                -math.cos(twist) * tangential_direction - math.sin(twist) * rotation_axis
            )
            chord_direction = chord_direction / np.linalg.norm(chord_direction)

            station_nodes: list[Node] = []
            for chord_fraction in chord_grid:
                point = (
                    radial_direction * radius
                    + rotation_axis * prebend_i
                    + tangential_direction * sweep_i
                    + (chord_fraction - pitch_axis) * chord * chord_direction
                )
                node = Node(point)
                mesh.add_node(node)
                station_nodes.append(node)
            blade_station_nodes.append(station_nodes)

        blade_elements: set[MeshElement] = set()
        blade_nodes = {node for station_nodes in blade_station_nodes for node in station_nodes}
        for span_index in range(n_spanwise_panels):
            for chord_index in range(n_chordwise_panels):
                element = MeshElement(
                    [
                        blade_station_nodes[span_index][chord_index],
                        blade_station_nodes[span_index][chord_index + 1],
                        blade_station_nodes[span_index + 1][chord_index + 1],
                        blade_station_nodes[span_index + 1][chord_index],
                    ],
                    ElementType.quad,
                )
                mesh.add_element(element)
                blade_elements.add(element)

        set_name = f"rotor_blade_{blade_index + 1}"
        mesh.add_node_set(NodeSet(set_name, blade_nodes))
        mesh.add_element_set(ElementSet(set_name, blade_elements))

    mesh.renumber_mesh("simple")
    return mesh


def _build_vlm_runtime_context(
    blade_aero: BladeAero,
    blade_file: str,
    mesh: MeshModel,
    wind_speed: float,
    *,
    rho: float,
    wake_length_scale: float,
    wake_history_steps: int,
    wake_corrector_iterations: int,
    diagonal_regularization: float,
    core_radius: float,
):
    cfg = {
        "participant": "Fluid",
        "coupling_mesh": "Fluid-Mesh",
        "mesh": {
            "source": "generator",
            "generator": {
                "type": "RotorMesh",
                "params": {"yaml_file": blade_file},
            },
        },
        "aero": {
            "backend": "vlm",
            "blade_file": blade_file,
            "rotor": {
                "n_blades": blade_aero.n_blades,
                "hub_radius": blade_aero.hub_radius,
                "rotation_axis": [0.0, 1.0, 0.0],
                "rotation_center": [0.0, 0.0, 0.0],
            },
            "backend_config": {
                "air_density": rho,
                "wind_velocity": [0.0, -wind_speed, 0.0],
                "wake_length_scale": wake_length_scale,
                "wake_history_steps": wake_history_steps,
                "wake_corrector_iterations": wake_corrector_iterations,
                "diagonal_regularization": diagonal_regularization,
                "core_radius": core_radius,
            },
        },
    }
    return build_aero_runtime_context(mesh, cfg)


def _compute_rotor_performance_metrics(
    blade_aero: BladeAero,
    operating_point: RotorOperatingPoint,
    integrated_force: np.ndarray,
    integrated_moment: np.ndarray,
    *,
    rho: float,
) -> RotorPerformanceMetrics:
    rotation_axis = np.array([0.0, 1.0, 0.0], dtype=float)
    thrust = abs(float(np.dot(integrated_force, rotation_axis)))
    torque = abs(float(np.dot(integrated_moment, rotation_axis)))
    power = torque * operating_point.omega_rpm * 2.0 * math.pi / 60.0

    rotor_area = math.pi * blade_aero.rotor_radius**2
    dynamic_pressure_area = 0.5 * rho * rotor_area * operating_point.wind_speed**2
    cp = power / (dynamic_pressure_area * operating_point.wind_speed)
    ct = thrust / dynamic_pressure_area
    cq = torque / (dynamic_pressure_area * blade_aero.rotor_radius)

    return RotorPerformanceMetrics(
        thrust=thrust,
        torque=torque,
        power=power,
        cp=cp,
        ct=ct,
        cq=cq,
    )


def _build_reduced_vlm_case(
    blade: BladeAero | str | Path,
    operating_point: RotorOperatingPoint,
    *,
    rho: float,
    n_spanwise_panels: int,
    n_chordwise_panels: int,
    wake_length_scale: float,
    wake_history_steps: int,
    wake_corrector_iterations: int,
    diagonal_regularization: float,
    core_radius: float,
) -> tuple[BladeAero, MeshModel, object, object]:
    blade_aero, blade_file = _coerce_blade_aero(blade)
    mesh = build_reduced_vlm_rotor_mesh(
        blade_aero,
        n_spanwise_panels=n_spanwise_panels,
        n_chordwise_panels=n_chordwise_panels,
        collective_pitch_deg=operating_point.pitch_deg,
    )
    runtime_context = _build_vlm_runtime_context(
        blade_aero,
        blade_file,
        mesh,
        operating_point.wind_speed,
        rho=rho,
        wake_length_scale=wake_length_scale,
        wake_history_steps=wake_history_steps,
        wake_corrector_iterations=wake_corrector_iterations,
        diagonal_regularization=diagonal_regularization,
        core_radius=core_radius,
    )
    backend = build_vlm_backend(runtime_context)

    nodal_displacements = np.zeros((mesh.node_count, 3), dtype=float)
    state = RotorAeroState(
        dt=0.0,
        omega_rad_s=operating_point.omega_rpm * 2.0 * math.pi / 60.0,
        wind_velocity=np.array([0.0, -operating_point.wind_speed, 0.0], dtype=float),
        nodal_displacements=nodal_displacements,
    )
    loads = backend.compute_loads(runtime_context, state)

    if runtime_context.surface is None:
        raise RuntimeError("Reduced VLM validation requires an aerodynamic surface")
    surface_state = runtime_context.surface.build_deformed_state(nodal_displacements)

    return blade_aero, mesh, surface_state, loads


def _local_panel_frame(
    blade_state: object, panel_index: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    node_indices = blade_state.panel_node_indices[panel_index]
    if len(node_indices) != 4:
        raise ValueError("Reduced VLM+polar correction expects quadrilateral panels")

    vertices = blade_state.node_coordinates[node_indices]
    chord_direction = 0.5 * ((vertices[1] - vertices[0]) + (vertices[2] - vertices[3]))
    chord_norm = np.linalg.norm(chord_direction)
    if np.isclose(chord_norm, 0.0):
        raise ValueError("Panel-local chord direction cannot be zero")
    chord_direction = chord_direction / chord_norm

    span_direction = 0.5 * ((vertices[3] - vertices[0]) + (vertices[2] - vertices[1]))
    span_direction = span_direction - np.dot(span_direction, chord_direction) * chord_direction
    span_norm = np.linalg.norm(span_direction)
    if np.isclose(span_norm, 0.0):
        raise ValueError("Panel-local span direction cannot be zero")
    span_direction = span_direction / span_norm

    thickness_direction = np.cross(span_direction, chord_direction)
    thickness_norm = np.linalg.norm(thickness_direction)
    if np.isclose(thickness_norm, 0.0):
        raise ValueError("Panel-local thickness direction cannot be zero")
    thickness_direction = thickness_direction / thickness_norm

    if np.dot(thickness_direction, blade_state.panel_normals[panel_index]) < 0.0:
        thickness_direction = -thickness_direction

    return chord_direction, span_direction, thickness_direction


def compute_vlm_sectional_diagnostics(
    blade: BladeAero | str | Path,
    operating_point: RotorOperatingPoint,
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> RotorSectionalDiagnostics:
    """Compute stripwise diagnostics for the reduced VLM validation mesh.

    The diagnostics are reported on blade 1 only. For the reduced rigid mesh,
    all blades are azimuth copies, so this captures the relevant strip behavior.

    Section-level force/moment diagnostics are reconstructed from local
    incident-flow lift/drag directions with polar Cl/Cd, which gives a more
    physically interpretable tangential loading trace than raw panel-normal
    pressure forces.
    """

    blade_aero, blade_file = _coerce_blade_aero(blade)
    mesh = build_reduced_vlm_rotor_mesh(
        blade_aero,
        n_spanwise_panels=n_spanwise_panels,
        n_chordwise_panels=n_chordwise_panels,
        collective_pitch_deg=operating_point.pitch_deg,
    )
    runtime_context = _build_vlm_runtime_context(
        blade_aero,
        blade_file,
        mesh,
        operating_point.wind_speed,
        rho=rho,
        wake_length_scale=wake_length_scale,
        wake_history_steps=wake_history_steps,
        wake_corrector_iterations=wake_corrector_iterations,
        diagonal_regularization=diagonal_regularization,
        core_radius=core_radius,
    )
    backend = build_vlm_backend(runtime_context)

    nodal_displacements = np.zeros((mesh.node_count, 3), dtype=float)
    state = RotorAeroState(
        dt=0.0,
        omega_rad_s=operating_point.omega_rpm * 2.0 * math.pi / 60.0,
        wind_velocity=np.array([0.0, -operating_point.wind_speed, 0.0], dtype=float),
        nodal_displacements=nodal_displacements,
    )

    if runtime_context.surface is None:
        raise RuntimeError("Reduced VLM diagnostics require an aerodynamic surface")
    surface_state = runtime_context.surface.build_deformed_state(nodal_displacements)
    panels = _flatten_surface_panels(surface_state, state)
    wake_length = backend._resolve_wake_length(panels)
    wake_rows = backend._build_working_wake_rows(panels, state.dt)
    influence = backend._assemble_influence_matrix(panels, wake_rows, wake_length)
    circulation = backend._solve_circulation(influence, panels)
    induced_velocity = backend._induced_velocity_on_panels(
        panels,
        circulation,
        wake_rows,
        wake_length,
    )
    blade_state = surface_state.blades[0]
    strip_count = blade_state.n_panels // n_chordwise_panels

    radii = np.zeros(strip_count, dtype=float)
    alpha_incident_deg = np.zeros(strip_count, dtype=float)
    alpha_total_deg = np.zeros(strip_count, dtype=float)
    cl_polar = np.zeros(strip_count, dtype=float)
    cd_polar = np.zeros(strip_count, dtype=float)
    cl_vlm_trailing = np.zeros(strip_count, dtype=float)
    gamma_leading = np.zeros(strip_count, dtype=float)
    gamma_trailing = np.zeros(strip_count, dtype=float)
    strip_base_thrust = np.zeros(strip_count, dtype=float)
    strip_base_torque = np.zeros(strip_count, dtype=float)

    for strip_index in range(strip_count):
        local_leading = strip_index * n_chordwise_panels
        local_trailing = local_leading + (n_chordwise_panels - 1)
        global_leading = local_leading
        global_trailing = local_trailing

        chord_direction, span_direction, thickness_direction = _local_panel_frame(
            blade_state,
            local_trailing,
        )

        total_velocity = panels[global_trailing].relative_flow + induced_velocity[global_trailing]
        self_induced = circulation[global_trailing] * _wake_induced_velocity(
            panels[global_trailing].collocation_point,
            panels[global_trailing],
            global_trailing,
            wake_rows,
            wake_length,
            backend._core_radius,
        )
        incident_velocity = total_velocity - self_induced

        section_velocity_total = (
            total_velocity - np.dot(total_velocity, span_direction) * span_direction
        )
        section_velocity_incident = (
            incident_velocity - np.dot(incident_velocity, span_direction) * span_direction
        )
        speed_incident = np.linalg.norm(section_velocity_incident)

        alpha_total = -math.atan2(
            np.dot(section_velocity_total, thickness_direction),
            np.dot(section_velocity_total, chord_direction),
        )
        alpha_incident = -math.atan2(
            np.dot(section_velocity_incident, thickness_direction),
            np.dot(section_velocity_incident, chord_direction),
        )

        strip_radii = [
            float(blade_state.panel_radii[local_leading + offset])
            for offset in range(n_chordwise_panels)
        ]
        strip_chord = sum(
            float(blade_state.panel_chord_lengths[local_leading + offset])
            for offset in range(n_chordwise_panels)
        )
        station_index = int(np.argmin(np.abs(blade_aero.r - np.mean(strip_radii))))
        reynolds = max(rho * speed_incident * strip_chord / mu, 1.0)
        polar = blade_aero.stations[station_index].airfoil.get_polar(reynolds)
        cl_eval, cd_eval, _cm_eval = polar.evaluate(np.array([alpha_incident], dtype=float))

        strip_force = np.zeros(3, dtype=float)
        strip_moment = np.zeros(3, dtype=float)
        for offset in range(n_chordwise_panels):
            panel_local_index = local_leading + offset
            panel_global_index = panel_local_index
            panel_total_velocity = (
                panels[panel_global_index].relative_flow + induced_velocity[panel_global_index]
            )
            panel_self_induced = circulation[panel_global_index] * _wake_induced_velocity(
                panels[panel_global_index].collocation_point,
                panels[panel_global_index],
                panel_global_index,
                wake_rows,
                wake_length,
                backend._core_radius,
            )
            panel_incident_velocity = panel_total_velocity - panel_self_induced

            panel_chord_direction, panel_span_direction, panel_thickness_direction = (
                _local_panel_frame(blade_state, panel_local_index)
            )
            panel_section_velocity = panel_incident_velocity - (
                np.dot(panel_incident_velocity, panel_span_direction) * panel_span_direction
            )
            panel_section_speed = np.linalg.norm(panel_section_velocity)
            if np.isclose(panel_section_speed, 0.0):
                continue

            panel_drag_direction = panel_section_velocity / panel_section_speed
            panel_alpha = -math.atan2(
                np.dot(panel_section_velocity, panel_thickness_direction),
                np.dot(panel_section_velocity, panel_chord_direction),
            )
            panel_radius = float(blade_state.panel_radii[panel_local_index])
            panel_station_index = int(np.argmin(np.abs(blade_aero.r - panel_radius)))
            panel_chord = float(blade_state.panel_chord_lengths[panel_local_index])
            panel_area = float(blade_state.panel_areas[panel_local_index])
            panel_reynolds = max(rho * panel_section_speed * panel_chord / mu, 1.0)
            panel_polar = blade_aero.stations[panel_station_index].airfoil.get_polar(panel_reynolds)
            panel_cl, panel_cd, _panel_cm = panel_polar.evaluate(
                np.array([panel_alpha], dtype=float)
            )

            panel_lift_direction = np.cross(panel_span_direction, panel_drag_direction)
            panel_lift_norm = np.linalg.norm(panel_lift_direction)
            if np.isclose(panel_lift_norm, 0.0):
                continue
            panel_lift_direction = panel_lift_direction / panel_lift_norm

            panel_dynamic = 0.5 * rho * panel_section_speed**2 * panel_area
            panel_force = (
                panel_dynamic * float(panel_cl[0]) * panel_lift_direction
                + panel_dynamic * float(panel_cd[0]) * panel_drag_direction
            )
            panel_centroid = blade_state.panel_centroids[panel_local_index]

            strip_force += panel_force
            strip_moment += np.cross(panel_centroid - surface_state.rotation_center, panel_force)

        radii[strip_index] = float(np.mean(strip_radii))
        alpha_incident_deg[strip_index] = math.degrees(alpha_incident)
        alpha_total_deg[strip_index] = math.degrees(alpha_total)
        cl_polar[strip_index] = float(cl_eval[0])
        cd_polar[strip_index] = float(cd_eval[0])
        cl_vlm_trailing[strip_index] = (
            2.0
            * float(circulation[global_trailing])
            / max(speed_incident * float(blade_state.panel_chord_lengths[local_trailing]), 1.0e-9)
        )
        gamma_leading[strip_index] = float(circulation[global_leading])
        gamma_trailing[strip_index] = float(circulation[global_trailing])
        strip_base_thrust[strip_index] = float(
            np.dot(strip_force, np.array([0.0, 1.0, 0.0], dtype=float))
        )
        strip_base_torque[strip_index] = float(
            np.dot(strip_moment, np.array([0.0, 1.0, 0.0], dtype=float))
        )

    return RotorSectionalDiagnostics(
        operating_point=operating_point,
        radii=radii,
        alpha_incident_deg=alpha_incident_deg,
        alpha_total_deg=alpha_total_deg,
        cl_polar=cl_polar,
        cd_polar=cd_polar,
        cl_vlm_trailing=cl_vlm_trailing,
        gamma_leading=gamma_leading,
        gamma_trailing=gamma_trailing,
        strip_base_thrust=strip_base_thrust,
        strip_base_torque=strip_base_torque,
    )


def compute_vlm_reduced_rotor_metrics(
    blade: BladeAero | str | Path,
    operating_point: RotorOperatingPoint,
    *,
    rho: float = 1.225,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> tuple[RotorPerformanceMetrics, int]:
    """Evaluate the current VLM backend on a reduced full-rotor lattice."""

    blade_aero, mesh, _surface_state, loads = _build_reduced_vlm_case(
        blade,
        operating_point,
        rho=rho,
        n_spanwise_panels=n_spanwise_panels,
        n_chordwise_panels=n_chordwise_panels,
        wake_length_scale=wake_length_scale,
        wake_history_steps=wake_history_steps,
        wake_corrector_iterations=wake_corrector_iterations,
        diagonal_regularization=diagonal_regularization,
        core_radius=core_radius,
    )

    return (
        _compute_rotor_performance_metrics(
            blade_aero,
            operating_point,
            loads.integrated_force,
            loads.integrated_moment,
            rho=rho,
        ),
        mesh.elements_count,
    )


def compute_vlm_polar_corrected_reduced_rotor_metrics(
    blade: BladeAero | str | Path,
    operating_point: RotorOperatingPoint,
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    lift_mode: PolarCorrectionLiftMode = "drag_only",
    experimental_lift_gain: float = 0.15,
    experimental_alpha_gate_deg: float = 6.0,
    experimental_tsr_cutoff: float = 9.5,
    experimental_tsr_falloff: float = 3.0,
    experimental_radial_tip_scale: float = 1.0,
    experimental_radial_power: float = 1.0,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> tuple[RotorPerformanceMetrics, int]:
    """Apply a local VLM-plus-polar force correction on the reduced rotor.

    This validation-only prototype preserves the baseline VLM panel force and
    adds a sectional profile-drag correction from BladeAero polars using the
    panel-local incident flow with self-induction removed.

    `lift_mode` controls optional lift closure experiments:
    - `drag_only`: current stable baseline (default)
    - `delta_cl_trailing`: non-default experimental branch that adds a
      trailing-panel lift term proportional to (Cl_polar - Cl_vlm_trailing)
    - `delta_cl_trailing_positive_gated`: non-default experimental branch
      that adds lift only when Cl_polar > Cl_vlm_trailing and |alpha| exceeds
      a gate threshold, to reduce high-TSR over-correction.
        - `delta_cl_trailing_positive_gated_tsr_weighted`: gated positive-delta-Cl
            branch with an additional global TSR damping term above
            `experimental_tsr_cutoff`; optional radial weighting can be applied using
            `experimental_radial_tip_scale` and `experimental_radial_power`.
    """

    if lift_mode not in (
        "drag_only",
        "delta_cl_trailing",
        "delta_cl_trailing_positive_gated",
        "delta_cl_trailing_positive_gated_tsr_weighted",
    ):
        raise ValueError(f"Unsupported lift_mode '{lift_mode}'")
    if experimental_lift_gain < 0.0:
        raise ValueError("experimental_lift_gain must be >= 0")
    if experimental_alpha_gate_deg <= 0.0:
        raise ValueError("experimental_alpha_gate_deg must be > 0")
    if experimental_tsr_cutoff <= 0.0:
        raise ValueError("experimental_tsr_cutoff must be > 0")
    if experimental_tsr_falloff <= 0.0:
        raise ValueError("experimental_tsr_falloff must be > 0")
    if not (0.0 < experimental_radial_tip_scale <= 1.0):
        raise ValueError("experimental_radial_tip_scale must be in (0, 1]")
    if experimental_radial_power <= 0.0:
        raise ValueError("experimental_radial_power must be > 0")

    blade_aero, mesh, surface_state, loads = _build_reduced_vlm_case(
        blade,
        operating_point,
        rho=rho,
        n_spanwise_panels=n_spanwise_panels,
        n_chordwise_panels=n_chordwise_panels,
        wake_length_scale=wake_length_scale,
        wake_history_steps=wake_history_steps,
        wake_corrector_iterations=wake_corrector_iterations,
        diagonal_regularization=diagonal_regularization,
        core_radius=core_radius,
    )

    panel_incident_velocities_raw = loads.metadata.get("panel_incident_velocities")
    panel_forces_raw = loads.metadata.get("panel_forces")
    circulation_raw = loads.metadata.get("circulation")
    if panel_incident_velocities_raw is None or panel_forces_raw is None or circulation_raw is None:
        raise RuntimeError(
            "VLM validation metadata is missing panel_incident_velocities, panel_forces, or circulation"
        )

    panel_incident_velocities = np.asarray(panel_incident_velocities_raw, dtype=float).reshape(
        -1, 3
    )
    panel_forces = np.asarray(panel_forces_raw, dtype=float).reshape(-1, 3)
    circulation = np.asarray(circulation_raw, dtype=float).reshape(-1)
    station_radii = blade_aero.r

    integrated_force = np.zeros(3, dtype=float)
    integrated_moment = np.zeros(3, dtype=float)
    panel_cursor = 0
    alpha_gate_rad = math.radians(experimental_alpha_gate_deg)
    omega = operating_point.omega_rpm * 2.0 * math.pi / 60.0
    rotor_tsr = omega * blade_aero.rotor_radius / max(operating_point.wind_speed, 1.0e-9)
    tsr_excess = max(rotor_tsr - experimental_tsr_cutoff, 0.0)
    tsr_weight = 1.0 / (1.0 + (tsr_excess / experimental_tsr_falloff) ** 2)
    radius_min = float(np.min(blade_aero.r))
    radius_span = float(max(np.max(blade_aero.r) - radius_min, 1.0e-9))

    for blade_state in surface_state.blades:
        for panel_index in range(blade_state.n_panels):
            incident_velocity = panel_incident_velocities[panel_cursor]
            chord_direction, span_direction, thickness_direction = _local_panel_frame(
                blade_state,
                panel_index,
            )
            section_velocity = (
                incident_velocity - np.dot(incident_velocity, span_direction) * span_direction
            )
            section_speed = np.linalg.norm(section_velocity)
            if np.isclose(section_speed, 0.0):
                integrated_force += panel_forces[panel_cursor]
                integrated_moment += np.cross(
                    blade_state.panel_centroids[panel_index] - surface_state.rotation_center,
                    panel_forces[panel_cursor],
                )
                panel_cursor += 1
                continue

            drag_direction = section_velocity / section_speed
            alpha = -math.atan2(
                np.dot(section_velocity, thickness_direction),
                np.dot(section_velocity, chord_direction),
            )
            panel_radius = float(blade_state.panel_radii[panel_index])
            station_index = int(np.argmin(np.abs(station_radii - panel_radius)))
            panel_chord = float(blade_state.panel_chord_lengths[panel_index])
            panel_area = float(blade_state.panel_areas[panel_index])
            reynolds = max(rho * section_speed * panel_chord / mu, 1.0)
            polar = blade_aero.stations[station_index].airfoil.get_polar(reynolds)
            cl, cd, _cm = polar.evaluate(np.array([alpha], dtype=float))
            drag_magnitude = 0.5 * rho * section_speed**2 * panel_area * float(cd[0])
            corrected_panel_force = panel_forces[panel_cursor] + drag_magnitude * drag_direction

            if (
                lift_mode
                in (
                    "delta_cl_trailing",
                    "delta_cl_trailing_positive_gated",
                    "delta_cl_trailing_positive_gated_tsr_weighted",
                )
                and n_chordwise_panels > 1
                and (panel_index + 1) % n_chordwise_panels == 0
            ):
                cl_vlm_trailing = (
                    2.0
                    * float(circulation[panel_cursor])
                    / max(
                        section_speed * panel_chord,
                        1.0e-9,
                    )
                )
                delta_cl = float(cl[0]) - cl_vlm_trailing
                gate = 1.0
                if lift_mode == "delta_cl_trailing_positive_gated":
                    delta_cl = max(delta_cl, 0.0)
                    if abs(alpha) > alpha_gate_rad:
                        gate = min((abs(alpha) - alpha_gate_rad) / alpha_gate_rad, 1.0)
                    else:
                        gate = 0.0
                elif lift_mode == "delta_cl_trailing_positive_gated_tsr_weighted":
                    delta_cl = max(delta_cl, 0.0)
                    if abs(alpha) > alpha_gate_rad:
                        gate = min((abs(alpha) - alpha_gate_rad) / alpha_gate_rad, 1.0)
                    else:
                        gate = 0.0
                    radius_norm = min(max((panel_radius - radius_min) / radius_span, 0.0), 1.0)
                    radial_weight = (
                        experimental_radial_tip_scale
                        + (1.0 - experimental_radial_tip_scale)
                        * (1.0 - radius_norm) ** experimental_radial_power
                    )
                    gate *= tsr_weight
                    gate *= radial_weight
                lift_direction = np.cross(span_direction, drag_direction)
                lift_norm = np.linalg.norm(lift_direction)
                if not np.isclose(lift_norm, 0.0) and not np.isclose(gate, 0.0):
                    lift_direction = lift_direction / lift_norm
                    if np.dot(lift_direction, thickness_direction) < 0.0:
                        lift_direction = -lift_direction
                    lift_magnitude = (
                        gate
                        * experimental_lift_gain
                        * 0.5
                        * rho
                        * section_speed**2
                        * panel_area
                        * delta_cl
                    )
                    corrected_panel_force += lift_magnitude * lift_direction

            integrated_force += corrected_panel_force
            integrated_moment += np.cross(
                blade_state.panel_centroids[panel_index] - surface_state.rotation_center,
                corrected_panel_force,
            )
            panel_cursor += 1

    return (
        _compute_rotor_performance_metrics(
            blade_aero,
            operating_point,
            integrated_force,
            integrated_moment,
            rho=rho,
        ),
        mesh.elements_count,
    )


def compute_bem_metrics(
    blade: BladeAero | str | Path,
    operating_point: RotorOperatingPoint,
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
) -> RotorPerformanceMetrics:
    """Evaluate the BEM solver at one operating point."""

    blade_aero, _ = _coerce_blade_aero(blade)
    from aeroelast.solvers.bem.engine import BEMSolver

    bem = BEMSolver(blade_aero, rho=rho, mu=mu)
    result = bem.compute(
        v_inf=operating_point.wind_speed,
        omega=operating_point.omega_rpm,
        pitch=operating_point.pitch_deg,
    )

    rotor_area = math.pi * blade_aero.rotor_radius**2
    dynamic_pressure_area = 0.5 * rho * rotor_area * operating_point.wind_speed**2
    omega_rad_s = operating_point.omega_rpm * 2.0 * math.pi / 60.0
    cp = (
        result.CP
        if result.CP is not None
        else result.power / (dynamic_pressure_area * operating_point.wind_speed)
    )
    ct = result.CT if result.CT is not None else result.thrust / dynamic_pressure_area
    cq = (
        result.CQ
        if result.CQ is not None
        else result.torque / (dynamic_pressure_area * blade_aero.rotor_radius)
    )

    return RotorPerformanceMetrics(
        thrust=float(result.thrust),
        torque=float(result.torque),
        power=float(result.power if result.power is not None else result.torque * omega_rad_s),
        cp=float(cp),
        ct=float(ct),
        cq=float(cq),
    )


def _safe_ratio(numerator: float, denominator: float) -> float:
    if math.isclose(denominator, 0.0, abs_tol=1.0e-12):
        return math.nan
    return numerator / denominator


def _build_model_comparison(
    operating_point: RotorOperatingPoint,
    model_metrics: RotorPerformanceMetrics,
    bem_metrics: RotorPerformanceMetrics,
    panel_count: int,
) -> RotorModelComparison:
    return RotorModelComparison(
        operating_point=operating_point,
        vlm=model_metrics,
        bem=bem_metrics,
        panel_count=panel_count,
        thrust_ratio=_safe_ratio(model_metrics.thrust, bem_metrics.thrust),
        torque_ratio=_safe_ratio(model_metrics.torque, bem_metrics.torque),
        power_ratio=_safe_ratio(model_metrics.power, bem_metrics.power),
        cp_ratio=_safe_ratio(model_metrics.cp, bem_metrics.cp),
        ct_ratio=_safe_ratio(model_metrics.ct, bem_metrics.ct),
        cq_ratio=_safe_ratio(model_metrics.cq, bem_metrics.cq),
        cp_delta=model_metrics.cp - bem_metrics.cp,
        ct_delta=model_metrics.ct - bem_metrics.ct,
    )


def compare_vlm_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> list[RotorModelComparison]:
    """Compare reduced-order VLM and BEM over a set of rigid operating points."""

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        vlm_metrics, panel_count = compute_vlm_reduced_rotor_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            n_spanwise_panels=n_spanwise_panels,
            n_chordwise_panels=n_chordwise_panels,
            wake_length_scale=wake_length_scale,
            wake_history_steps=wake_history_steps,
            wake_corrector_iterations=wake_corrector_iterations,
            diagonal_regularization=diagonal_regularization,
            core_radius=core_radius,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(operating_point, vlm_metrics, bem_metrics, panel_count)
        )
    return comparisons


def compare_sharpy_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    model: Literal["vlm", "uvlm"],
    route: str | Path,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    solver_settings: dict[str, object] | None = None,
) -> list[RotorModelComparison]:
    """Compare a SHARPy rotor case against BEM over rigid operating points.

    The same ``BladeAero`` definition feeds both branches. SHARPy may use
    solver-specific settings, but geometry, operating point and coefficient
    normalisation remain shared through this contract.
    """

    from .sharpy_numad import NuMADSharpyRotorCaseAdapter

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    adapter = NuMADSharpyRotorCaseAdapter(blade_aero=blade_aero, route=Path(route), rho=rho)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        sharpy_result = adapter.run(
            operating_point,
            model=model,
            solver_settings=solver_settings,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(
                operating_point,
                sharpy_result.metrics,
                bem_metrics,
                sharpy_result.panel_count,
            )
        )
    return comparisons


def compare_vlm_polar_corrected_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> list[RotorModelComparison]:
    """Compare the reduced VLM+polar correction against BEM over a sweep."""

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        corrected_metrics, panel_count = compute_vlm_polar_corrected_reduced_rotor_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
            n_spanwise_panels=n_spanwise_panels,
            n_chordwise_panels=n_chordwise_panels,
            wake_length_scale=wake_length_scale,
            wake_history_steps=wake_history_steps,
            wake_corrector_iterations=wake_corrector_iterations,
            diagonal_regularization=diagonal_regularization,
            core_radius=core_radius,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(
                operating_point,
                corrected_metrics,
                bem_metrics,
                panel_count,
            )
        )
    return comparisons


def compare_vlm_experimental_lift_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    experimental_lift_gain: float = 0.15,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> list[RotorModelComparison]:
    """Compare experimental delta-Cl trailing lift closure against BEM."""

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        corrected_metrics, panel_count = compute_vlm_polar_corrected_reduced_rotor_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
            lift_mode="delta_cl_trailing",
            experimental_lift_gain=experimental_lift_gain,
            n_spanwise_panels=n_spanwise_panels,
            n_chordwise_panels=n_chordwise_panels,
            wake_length_scale=wake_length_scale,
            wake_history_steps=wake_history_steps,
            wake_corrector_iterations=wake_corrector_iterations,
            diagonal_regularization=diagonal_regularization,
            core_radius=core_radius,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(
                operating_point,
                corrected_metrics,
                bem_metrics,
                panel_count,
            )
        )
    return comparisons


def compare_vlm_experimental_gated_lift_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    experimental_lift_gain: float = 1.0,
    experimental_alpha_gate_deg: float = 6.0,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> list[RotorModelComparison]:
    """Compare gated positive-delta-Cl trailing lift closure against BEM."""

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        corrected_metrics, panel_count = compute_vlm_polar_corrected_reduced_rotor_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
            lift_mode="delta_cl_trailing_positive_gated",
            experimental_lift_gain=experimental_lift_gain,
            experimental_alpha_gate_deg=experimental_alpha_gate_deg,
            n_spanwise_panels=n_spanwise_panels,
            n_chordwise_panels=n_chordwise_panels,
            wake_length_scale=wake_length_scale,
            wake_history_steps=wake_history_steps,
            wake_corrector_iterations=wake_corrector_iterations,
            diagonal_regularization=diagonal_regularization,
            core_radius=core_radius,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(
                operating_point,
                corrected_metrics,
                bem_metrics,
                panel_count,
            )
        )
    return comparisons


def compare_vlm_experimental_gated_tsr_lift_to_bem_sweep(
    blade: BladeAero | str | Path,
    operating_points: Iterable[RotorOperatingPoint],
    *,
    rho: float = 1.225,
    mu: float = 1.81206e-5,
    experimental_lift_gain: float = 1.0,
    experimental_alpha_gate_deg: float = 6.0,
    experimental_tsr_cutoff: float = 9.5,
    experimental_tsr_falloff: float = 3.0,
    experimental_radial_tip_scale: float = 1.0,
    experimental_radial_power: float = 1.0,
    n_spanwise_panels: int = 12,
    n_chordwise_panels: int = 2,
    wake_length_scale: float = 20.0,
    wake_history_steps: int = 0,
    wake_corrector_iterations: int = 0,
    diagonal_regularization: float = 1.0e-6,
    core_radius: float = 1.0e-4,
) -> list[RotorModelComparison]:
    """Compare gated positive-delta-Cl branch with TSR damping against BEM."""

    blade_aero, _blade_file = _coerce_blade_aero(blade)
    comparisons: list[RotorModelComparison] = []
    for operating_point in operating_points:
        corrected_metrics, panel_count = compute_vlm_polar_corrected_reduced_rotor_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
            lift_mode="delta_cl_trailing_positive_gated_tsr_weighted",
            experimental_lift_gain=experimental_lift_gain,
            experimental_alpha_gate_deg=experimental_alpha_gate_deg,
            experimental_tsr_cutoff=experimental_tsr_cutoff,
            experimental_tsr_falloff=experimental_tsr_falloff,
            experimental_radial_tip_scale=experimental_radial_tip_scale,
            experimental_radial_power=experimental_radial_power,
            n_spanwise_panels=n_spanwise_panels,
            n_chordwise_panels=n_chordwise_panels,
            wake_length_scale=wake_length_scale,
            wake_history_steps=wake_history_steps,
            wake_corrector_iterations=wake_corrector_iterations,
            diagonal_regularization=diagonal_regularization,
            core_radius=core_radius,
        )
        bem_metrics = compute_bem_metrics(
            blade_aero,
            operating_point,
            rho=rho,
            mu=mu,
        )
        comparisons.append(
            _build_model_comparison(
                operating_point,
                corrected_metrics,
                bem_metrics,
                panel_count,
            )
        )
    return comparisons


def comparison_as_dict(comparison: RotorModelComparison) -> dict[str, float | str]:
    """Convert a comparison row to a flat mapping for CSV or tabular reporting."""

    point = comparison.operating_point
    return {
        "Case": point.name,
        "Wind Speed [m/s]": point.wind_speed,
        "Omega [rpm]": point.omega_rpm,
        "Pitch [deg]": point.pitch_deg,
        "VLM Panels [-]": comparison.panel_count,
        "BEM Thrust [N]": comparison.bem.thrust,
        "VLM Thrust [N]": comparison.vlm.thrust,
        "Thrust Ratio [-]": comparison.thrust_ratio,
        "BEM Torque [N m]": comparison.bem.torque,
        "VLM Torque [N m]": comparison.vlm.torque,
        "Torque Ratio [-]": comparison.torque_ratio,
        "BEM Power [W]": comparison.bem.power,
        "VLM Power [W]": comparison.vlm.power,
        "Power Ratio [-]": comparison.power_ratio,
        "BEM CP [-]": comparison.bem.cp,
        "VLM CP [-]": comparison.vlm.cp,
        "CP Ratio [-]": comparison.cp_ratio,
        "CP Delta [-]": comparison.cp_delta,
        "BEM CT [-]": comparison.bem.ct,
        "VLM CT [-]": comparison.vlm.ct,
        "CT Ratio [-]": comparison.ct_ratio,
        "CT Delta [-]": comparison.ct_delta,
        "BEM CQ [-]": comparison.bem.cq,
        "VLM CQ [-]": comparison.vlm.cq,
        "CQ Ratio [-]": comparison.cq_ratio,
    }


def _mean_abs_unit_ratio_error(values: Iterable[float]) -> float:
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return math.nan
    return float(np.mean([abs(1.0 - value) for value in finite_values]))


def compute_validation_objective_metrics(
    reference_rows: Sequence[RotorModelComparison],
    tsr_rows: Sequence[RotorModelComparison],
) -> RotorValidationObjectiveMetrics:
    """Compute compact objective metrics for model-vs-BEM validation sweeps."""

    thrust_ratio_error = _mean_abs_unit_ratio_error(row.thrust_ratio for row in reference_rows)
    torque_ratio_error = _mean_abs_unit_ratio_error(row.torque_ratio for row in reference_rows)
    power_ratio_error = _mean_abs_unit_ratio_error(row.power_ratio for row in reference_rows)

    cp_errors = [abs(row.cp_delta) for row in reference_rows]
    mean_abs_cp_error = float(np.mean(cp_errors)) if cp_errors else math.nan

    zero_pitch_rows = [
        row for row in reference_rows if abs(row.operating_point.pitch_deg) < 1.0e-12
    ]
    zero_pitch_torque_error = _mean_abs_unit_ratio_error(
        row.torque_ratio for row in zero_pitch_rows
    )

    cp_deltas = [row.vlm.cp - row.bem.cp for row in tsr_rows]
    cp_curve_rmse = float(np.sqrt(np.mean(np.square(cp_deltas)))) if cp_deltas else math.nan

    bem_cp_peak = max(tsr_rows, key=lambda row: row.bem.cp)
    model_cp_peak = max(tsr_rows, key=lambda row: row.vlm.cp)

    return RotorValidationObjectiveMetrics(
        mean_abs_thrust_ratio_error=thrust_ratio_error,
        mean_abs_torque_ratio_error=torque_ratio_error,
        mean_abs_power_ratio_error=power_ratio_error,
        mean_abs_cp_error=mean_abs_cp_error,
        zero_pitch_mean_abs_torque_ratio_error=zero_pitch_torque_error,
        cp_curve_rmse=cp_curve_rmse,
        cp_peak_rpm_error=abs(
            model_cp_peak.operating_point.omega_rpm - bem_cp_peak.operating_point.omega_rpm
        ),
        cp_peak_value_error=abs(model_cp_peak.vlm.cp - bem_cp_peak.bem.cp),
        bem_cp_peak_rpm=bem_cp_peak.operating_point.omega_rpm,
        model_cp_peak_rpm=model_cp_peak.operating_point.omega_rpm,
        bem_cp_peak=bem_cp_peak.bem.cp,
        model_cp_peak=model_cp_peak.vlm.cp,
        reference_sample_count=len(reference_rows),
        tsr_sample_count=len(tsr_rows),
    )


__all__ = [
    "IEA15MW_REFERENCE_OPERATING_POINTS",
    "IEA15MW_ZERO_PITCH_TSR_SWEEP",
    "RotorModelComparison",
    "RotorOperatingPoint",
    "RotorPerformanceMetrics",
    "RotorSectionalDiagnostics",
    "RotorValidationObjectiveMetrics",
    "build_reduced_vlm_rotor_mesh",
    "compare_vlm_experimental_lift_to_bem_sweep",
    "compare_vlm_experimental_gated_lift_to_bem_sweep",
    "compare_vlm_experimental_gated_tsr_lift_to_bem_sweep",
    "compare_vlm_polar_corrected_to_bem_sweep",
    "compare_vlm_to_bem_sweep",
    "compare_sharpy_to_bem_sweep",
    "compute_vlm_sectional_diagnostics",
    "comparison_as_dict",
    "compute_validation_objective_metrics",
    "compute_bem_metrics",
    "compute_vlm_polar_corrected_reduced_rotor_metrics",
    "compute_vlm_reduced_rotor_metrics",
]
