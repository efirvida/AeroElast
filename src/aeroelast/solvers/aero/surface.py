"""Mesh-derived aerodynamic surface helpers for full-rotor participants."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aeroelast.core.mesh.entities import ElementType, MeshElement

from .rotor_geometry import FullRotorAeroGeometry, RotorBladeAeroGeometry

_SURFACE_ELEMENT_TYPES = {
    ElementType.triangle,
    ElementType.triangle6,
    ElementType.quad,
    ElementType.quad8,
    ElementType.quad9,
}


def _normalize_vector(vector: np.ndarray, label: str) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if np.isclose(norm, 0.0):
        raise ValueError(f"{label} cannot be the zero vector")
    return vector / norm


def _polygon_area_normal(vertices: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    centroid = vertices.mean(axis=0)
    area_vector = np.zeros(3, dtype=float)
    for idx in range(len(vertices)):
        p0 = vertices[idx] - centroid
        p1 = vertices[(idx + 1) % len(vertices)] - centroid
        area_vector += np.cross(p0, p1)
    area_vector *= 0.5

    area = float(np.linalg.norm(area_vector))
    if np.isclose(area, 0.0):
        raise ValueError("Degenerate surface element with zero area")

    return area, area_vector / area, centroid


def _surface_node_indices(element: MeshElement, node_index_map: dict[int, int]) -> np.ndarray:
    geometric_nodes = [node for node in element.nodes if node.geometric_node]
    if len(geometric_nodes) < 3:
        geometric_nodes = list(element.nodes)
    return np.array([node_index_map[node.id] for node in geometric_nodes], dtype=int)


def _estimate_surface_frame(
    coords: np.ndarray,
    rotation_axis: np.ndarray,
    rotation_center: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    relative_coords = coords - rotation_center
    radial_components = relative_coords - np.outer(relative_coords @ rotation_axis, rotation_axis)
    radial_distance = np.linalg.norm(radial_components, axis=1)
    order = np.argsort(radial_distance)
    sample_size = min(len(order), max(2, len(order) // 10 or 1))

    root_point = coords[order[:sample_size]].mean(axis=0)
    tip_point = coords[order[-sample_size:]].mean(axis=0)

    span_direction = tip_point - root_point
    span_direction -= np.dot(span_direction, rotation_axis) * rotation_axis
    if np.isclose(np.linalg.norm(span_direction), 0.0):
        span_direction = radial_components.mean(axis=0)
    span_direction = _normalize_vector(span_direction, "span_direction")

    centered = coords - coords.mean(axis=0)
    off_span = centered - np.outer(centered @ span_direction, span_direction)
    _, _, vh = np.linalg.svd(off_span, full_matrices=False)
    chord_direction = vh[0]
    chord_direction -= np.dot(chord_direction, span_direction) * span_direction
    chord_direction = _normalize_vector(chord_direction, "chord_direction")

    reference_chord = np.cross(rotation_axis, span_direction)
    if not np.isclose(np.linalg.norm(reference_chord), 0.0):
        reference_chord = reference_chord / np.linalg.norm(reference_chord)
        if np.dot(chord_direction, reference_chord) < 0.0:
            chord_direction = -chord_direction

    thickness_direction = _normalize_vector(
        np.cross(span_direction, chord_direction),
        "thickness_direction",
    )

    return root_point, span_direction, chord_direction, thickness_direction


def _build_panel_lattice(
    node_coordinates: np.ndarray,
    panel_node_indices: tuple[np.ndarray, ...],
    rotation_center: np.ndarray,
    root_point: np.ndarray,
    root_reference_radius: float,
    span_direction: np.ndarray,
    chord_direction: np.ndarray,
    thickness_direction: np.ndarray,
) -> dict[str, np.ndarray]:
    panel_centroids: list[np.ndarray] = []
    panel_normals: list[np.ndarray] = []
    panel_areas: list[float] = []
    panel_span_lengths: list[float] = []
    panel_chord_lengths: list[float] = []
    panel_collocation_points: list[np.ndarray] = []
    panel_bound_vortex_points: list[np.ndarray] = []

    for node_indices in panel_node_indices:
        panel_vertices = node_coordinates[node_indices]
        area, normal, centroid = _polygon_area_normal(panel_vertices)
        if np.dot(normal, thickness_direction) < 0.0:
            normal = -normal

        panel_span_coords = (panel_vertices - rotation_center) @ span_direction
        panel_chord_coords = (panel_vertices - root_point) @ chord_direction

        span_min = float(np.min(panel_span_coords))
        span_max = float(np.max(panel_span_coords))
        chord_min = float(np.min(panel_chord_coords))
        chord_max = float(np.max(panel_chord_coords))

        span_length = span_max - span_min
        chord_length = chord_max - chord_min
        if np.isclose(span_length, 0.0):
            raise ValueError("Surface panel has zero spanwise extent")
        if np.isclose(chord_length, 0.0):
            raise ValueError("Surface panel has zero chordwise extent")

        quarter_chord = chord_min + 0.25 * chord_length
        three_quarter_chord = chord_min + 0.75 * chord_length
        span_min_local = span_min - root_reference_radius
        span_max_local = span_max - root_reference_radius
        span_mid_local = 0.5 * (span_min_local + span_max_local)

        collocation_point = (
            root_point + span_mid_local * span_direction + three_quarter_chord * chord_direction
        )
        bound_vortex_root = (
            root_point + span_min_local * span_direction + quarter_chord * chord_direction
        )
        bound_vortex_tip = (
            root_point + span_max_local * span_direction + quarter_chord * chord_direction
        )

        panel_centroids.append(centroid)
        panel_normals.append(normal)
        panel_areas.append(area)
        panel_span_lengths.append(span_length)
        panel_chord_lengths.append(chord_length)
        panel_collocation_points.append(collocation_point)
        panel_bound_vortex_points.append(np.vstack((bound_vortex_root, bound_vortex_tip)))

    centroid_array = np.asarray(panel_centroids, dtype=float)
    node_radii = (node_coordinates - rotation_center) @ span_direction
    panel_radii = (centroid_array - rotation_center) @ span_direction

    return {
        "node_radii": node_radii,
        "panel_centroids": centroid_array,
        "panel_normals": np.asarray(panel_normals, dtype=float),
        "panel_areas": np.asarray(panel_areas, dtype=float),
        "panel_radii": panel_radii,
        "panel_spanwise_coordinates": panel_radii - root_reference_radius,
        "panel_chordwise_coordinates": (centroid_array - root_point) @ chord_direction,
        "panel_span_lengths": np.asarray(panel_span_lengths, dtype=float),
        "panel_chord_lengths": np.asarray(panel_chord_lengths, dtype=float),
        "panel_collocation_points": np.asarray(panel_collocation_points, dtype=float),
        "panel_bound_vortex_points": np.asarray(panel_bound_vortex_points, dtype=float),
    }


def _coerce_blade_displacements(displacements: np.ndarray, node_count: int) -> np.ndarray:
    arr = np.asarray(displacements, dtype=float)
    if arr.shape != (node_count, 3):
        raise ValueError(f"displacements must have shape ({node_count}, 3), got {arr.shape}")
    return arr


def _coerce_global_displacements(displacements: np.ndarray, total_parent_nodes: int) -> np.ndarray:
    arr = np.asarray(displacements, dtype=float)
    if arr.shape != (total_parent_nodes, 3):
        raise ValueError(
            f"displacements must have shape ({total_parent_nodes}, 3), got {arr.shape}"
        )
    return arr


@dataclass
class BladeAeroSurface:
    """Reduced aerodynamic surface built from one structural blade mesh."""

    blade_index: int
    azimuth_deg: float
    element_set_name: str
    reference_node_coordinates: np.ndarray
    parent_node_indices: np.ndarray
    parent_element_indices: np.ndarray
    total_parent_nodes: int
    rotation_axis: np.ndarray
    rotation_center: np.ndarray
    span_direction: np.ndarray
    chord_direction: np.ndarray
    thickness_direction: np.ndarray
    root_point: np.ndarray
    root_reference_radius: float
    node_reference_radii: np.ndarray
    panel_node_indices: tuple[np.ndarray, ...]
    panel_parent_node_indices: tuple[np.ndarray, ...]
    panel_centroids: np.ndarray
    panel_normals: np.ndarray
    panel_areas: np.ndarray
    panel_reference_radii: np.ndarray
    panel_spanwise_coordinates: np.ndarray
    panel_chordwise_coordinates: np.ndarray
    panel_span_lengths: np.ndarray
    panel_chord_lengths: np.ndarray
    panel_collocation_points: np.ndarray
    panel_bound_vortex_points: np.ndarray

    @property
    def n_panels(self) -> int:
        return len(self.panel_areas)

    @property
    def node_count(self) -> int:
        return len(self.parent_node_indices)

    def build_deformed_state(
        self,
        displacements: np.ndarray,
        *,
        scale: float = 1.0,
    ) -> "BladeAeroSurfaceState":
        """Build the blade aerodynamic surface state on the current deformed geometry."""
        blade_displacements = _coerce_blade_displacements(displacements, self.node_count)
        node_coordinates = self.reference_node_coordinates + scale * blade_displacements

        root_point, span_direction, chord_direction, thickness_direction = _estimate_surface_frame(
            node_coordinates,
            rotation_axis=self.rotation_axis,
            rotation_center=self.rotation_center,
        )
        if np.dot(span_direction, self.span_direction) < 0.0:
            span_direction = -span_direction
        if np.dot(chord_direction, self.chord_direction) < 0.0:
            chord_direction = -chord_direction

        thickness_direction = _normalize_vector(
            np.cross(span_direction, chord_direction),
            "thickness_direction",
        )
        if np.dot(thickness_direction, self.thickness_direction) < 0.0:
            thickness_direction = -thickness_direction
            chord_direction = -chord_direction

        lattice = _build_panel_lattice(
            node_coordinates,
            self.panel_node_indices,
            self.rotation_center,
            root_point,
            self.root_reference_radius,
            span_direction,
            chord_direction,
            thickness_direction,
        )

        return BladeAeroSurfaceState(
            blade_index=self.blade_index,
            azimuth_deg=self.azimuth_deg,
            element_set_name=self.element_set_name,
            parent_node_indices=self.parent_node_indices,
            parent_element_indices=self.parent_element_indices,
            total_parent_nodes=self.total_parent_nodes,
            node_coordinates=node_coordinates,
            rotation_axis=self.rotation_axis.copy(),
            rotation_center=self.rotation_center.copy(),
            span_direction=span_direction,
            chord_direction=chord_direction,
            thickness_direction=thickness_direction,
            root_point=root_point,
            root_reference_radius=self.root_reference_radius,
            node_radii=lattice["node_radii"],
            panel_node_indices=self.panel_node_indices,
            panel_parent_node_indices=self.panel_parent_node_indices,
            panel_centroids=lattice["panel_centroids"],
            panel_normals=lattice["panel_normals"],
            panel_areas=lattice["panel_areas"],
            panel_radii=lattice["panel_radii"],
            panel_spanwise_coordinates=lattice["panel_spanwise_coordinates"],
            panel_chordwise_coordinates=lattice["panel_chordwise_coordinates"],
            panel_span_lengths=lattice["panel_span_lengths"],
            panel_chord_lengths=lattice["panel_chord_lengths"],
            panel_collocation_points=lattice["panel_collocation_points"],
            panel_bound_vortex_points=lattice["panel_bound_vortex_points"],
        )

    def assemble_nodal_forces(self, panel_forces: np.ndarray) -> np.ndarray:
        """Assemble panel-integrated forces into blade-local nodal forces."""
        panel_forces = np.asarray(panel_forces, dtype=float)
        if panel_forces.shape != (self.n_panels, 3):
            raise ValueError(
                f"panel_forces must have shape ({self.n_panels}, 3), got {panel_forces.shape}"
            )

        nodal_forces = np.zeros((self.node_count, 3), dtype=float)
        for panel_index, node_indices in enumerate(self.panel_node_indices):
            share = panel_forces[panel_index] / len(node_indices)
            nodal_forces[node_indices] += share
        return nodal_forces

    def assemble_nodal_forces_from_pressure(self, panel_pressures: np.ndarray) -> np.ndarray:
        """Convert scalar pressure values into blade-local nodal forces."""
        panel_pressures = np.asarray(panel_pressures, dtype=float)
        if panel_pressures.shape != (self.n_panels,):
            raise ValueError(
                f"panel_pressures must have shape ({self.n_panels},), got {panel_pressures.shape}"
            )

        panel_forces = panel_pressures[:, None] * self.panel_areas[:, None] * self.panel_normals
        return self.assemble_nodal_forces(panel_forces)


@dataclass
class BladeAeroSurfaceState:
    """Aerodynamic blade surface state on the current deformed geometry."""

    blade_index: int
    azimuth_deg: float
    element_set_name: str
    parent_node_indices: np.ndarray
    parent_element_indices: np.ndarray
    total_parent_nodes: int
    node_coordinates: np.ndarray
    rotation_axis: np.ndarray
    rotation_center: np.ndarray
    span_direction: np.ndarray
    chord_direction: np.ndarray
    thickness_direction: np.ndarray
    root_point: np.ndarray
    root_reference_radius: float
    node_radii: np.ndarray
    panel_node_indices: tuple[np.ndarray, ...]
    panel_parent_node_indices: tuple[np.ndarray, ...]
    panel_centroids: np.ndarray
    panel_normals: np.ndarray
    panel_areas: np.ndarray
    panel_radii: np.ndarray
    panel_spanwise_coordinates: np.ndarray
    panel_chordwise_coordinates: np.ndarray
    panel_span_lengths: np.ndarray
    panel_chord_lengths: np.ndarray
    panel_collocation_points: np.ndarray
    panel_bound_vortex_points: np.ndarray

    @property
    def n_panels(self) -> int:
        return len(self.panel_areas)


@dataclass
class FullRotorAeroSurface:
    """Reduced aerodynamic surface for the full rotor."""

    rotation_axis: np.ndarray
    rotation_center: np.ndarray
    blades: tuple[BladeAeroSurface, ...]

    @property
    def n_blades(self) -> int:
        return len(self.blades)

    @property
    def n_panels(self) -> int:
        return sum(blade.n_panels for blade in self.blades)

    def get_blade(self, blade_index: int) -> BladeAeroSurface:
        return self.blades[blade_index]

    @property
    def total_parent_nodes(self) -> int:
        return max(blade.total_parent_nodes for blade in self.blades)

    def build_deformed_state(
        self,
        displacements: np.ndarray,
        *,
        scale: float = 1.0,
    ) -> "FullRotorAeroSurfaceState":
        """Build the full-rotor aerodynamic surface state on the deformed geometry."""
        global_displacements = _coerce_global_displacements(displacements, self.total_parent_nodes)
        blades = tuple(
            blade.build_deformed_state(global_displacements[blade.parent_node_indices], scale=scale)
            for blade in self.blades
        )
        return FullRotorAeroSurfaceState(
            rotation_axis=self.rotation_axis.copy(),
            rotation_center=self.rotation_center.copy(),
            blades=blades,
        )

    def assemble_global_nodal_forces(
        self,
        blade_panel_forces: tuple[np.ndarray, ...] | list[np.ndarray],
    ) -> np.ndarray:
        """Assemble blade-panel force arrays into the full-rotor node ordering."""
        if len(blade_panel_forces) != self.n_blades:
            raise ValueError(
                f"Expected {self.n_blades} blade force blocks, got {len(blade_panel_forces)}"
            )

        global_forces = np.zeros((self.total_parent_nodes, 3), dtype=float)
        for blade, panel_forces in zip(self.blades, blade_panel_forces):
            blade_nodal_forces = blade.assemble_nodal_forces(panel_forces)
            global_forces[blade.parent_node_indices] += blade_nodal_forces
        return global_forces

    def assemble_global_nodal_forces_from_pressure(
        self,
        blade_panel_pressures: tuple[np.ndarray, ...] | list[np.ndarray],
    ) -> np.ndarray:
        """Assemble blade-panel pressure arrays into the full-rotor node ordering."""
        if len(blade_panel_pressures) != self.n_blades:
            raise ValueError(
                f"Expected {self.n_blades} blade pressure blocks, got {len(blade_panel_pressures)}"
            )

        global_forces = np.zeros((self.total_parent_nodes, 3), dtype=float)
        for blade, panel_pressures in zip(self.blades, blade_panel_pressures):
            blade_nodal_forces = blade.assemble_nodal_forces_from_pressure(panel_pressures)
            global_forces[blade.parent_node_indices] += blade_nodal_forces
        return global_forces


@dataclass
class FullRotorAeroSurfaceState:
    """Aerodynamic surface state for the full rotor on the deformed geometry."""

    rotation_axis: np.ndarray
    rotation_center: np.ndarray
    blades: tuple[BladeAeroSurfaceState, ...]

    @property
    def n_blades(self) -> int:
        return len(self.blades)

    @property
    def n_panels(self) -> int:
        return sum(blade.n_panels for blade in self.blades)

    def get_blade(self, blade_index: int) -> BladeAeroSurfaceState:
        return self.blades[blade_index]


def build_blade_aero_surface(
    blade: RotorBladeAeroGeometry,
    rotation_axis: np.ndarray,
    rotation_center: np.ndarray,
    hub_radius: float | None = None,
) -> BladeAeroSurface:
    """Build a reduced aerodynamic surface from one blade mesh."""
    node_index_map = blade.mesh.node_id_to_index
    panel_node_indices: list[np.ndarray] = []
    panel_parent_node_indices: list[np.ndarray] = []

    if blade.element_count != blade.mesh.elements_count:
        raise RuntimeError(
            f"Blade '{blade.element_set_name}' lost element ordering between parent mesh and submesh"
        )

    for element in blade.mesh.elements:
        if element.element_type not in _SURFACE_ELEMENT_TYPES:
            continue

        local_node_indices = _surface_node_indices(element, node_index_map)
        panel_vertices = blade.reference_coords[local_node_indices]
        _polygon_area_normal(panel_vertices)
        panel_node_indices.append(local_node_indices)
        panel_parent_node_indices.append(blade.parent_node_indices[local_node_indices])

    if not panel_node_indices:
        raise ValueError(
            f"Blade '{blade.element_set_name}' does not contain supported surface elements"
        )

    root_point, span_direction, chord_direction, thickness_direction = _estimate_surface_frame(
        blade.reference_coords,
        rotation_axis,
        rotation_center,
    )
    root_reference_radius = (
        float(hub_radius)
        if hub_radius is not None
        else float(np.min((blade.reference_coords - rotation_center) @ span_direction))
    )
    root_reference_radius = max(0.0, root_reference_radius)

    lattice = _build_panel_lattice(
        blade.reference_coords,
        tuple(panel_node_indices),
        rotation_center,
        root_point,
        root_reference_radius,
        span_direction,
        chord_direction,
        thickness_direction,
    )

    return BladeAeroSurface(
        blade_index=blade.blade_index,
        azimuth_deg=blade.azimuth_deg,
        element_set_name=blade.element_set_name,
        reference_node_coordinates=blade.reference_coords.copy(),
        parent_node_indices=blade.parent_node_indices,
        parent_element_indices=blade.parent_element_indices,
        total_parent_nodes=blade.total_parent_nodes,
        rotation_axis=rotation_axis.copy(),
        rotation_center=rotation_center.copy(),
        span_direction=span_direction,
        chord_direction=chord_direction,
        thickness_direction=thickness_direction,
        root_point=root_point,
        root_reference_radius=root_reference_radius,
        node_reference_radii=lattice["node_radii"],
        panel_node_indices=tuple(panel_node_indices),
        panel_parent_node_indices=tuple(panel_parent_node_indices),
        panel_centroids=lattice["panel_centroids"],
        panel_normals=lattice["panel_normals"],
        panel_areas=lattice["panel_areas"],
        panel_reference_radii=lattice["panel_radii"],
        panel_spanwise_coordinates=lattice["panel_spanwise_coordinates"],
        panel_chordwise_coordinates=lattice["panel_chordwise_coordinates"],
        panel_span_lengths=lattice["panel_span_lengths"],
        panel_chord_lengths=lattice["panel_chord_lengths"],
        panel_collocation_points=lattice["panel_collocation_points"],
        panel_bound_vortex_points=lattice["panel_bound_vortex_points"],
    )


def build_full_rotor_aero_surface(geometry: FullRotorAeroGeometry) -> FullRotorAeroSurface:
    """Build reduced aerodynamic surfaces for all blades in a rotor."""
    blades = tuple(
        build_blade_aero_surface(
            blade,
            rotation_axis=geometry.rotation_axis,
            rotation_center=geometry.rotation_center,
            hub_radius=geometry.hub_radius,
        )
        for blade in geometry.blades
    )
    return FullRotorAeroSurface(
        rotation_axis=geometry.rotation_axis.copy(),
        rotation_center=geometry.rotation_center.copy(),
        blades=blades,
    )
