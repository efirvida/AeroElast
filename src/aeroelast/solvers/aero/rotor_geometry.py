"""Full-rotor aerodynamic geometry helpers.

This module extracts a reproducible blade-by-blade view of a rotor mesh so a
future aerodynamic backend can work on each blade explicitly while still
mapping fields back to the parent full-rotor ordering.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np

from aeroelast.core.mesh.model import MeshModel

ROTOR_BLADE_SET_PREFIX = "rotor_blade_"


def _normalize_vector(vector: Sequence[float]) -> np.ndarray:
    arr = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(arr)
    if np.isclose(norm, 0.0):
        raise ValueError("rotation_axis cannot be the zero vector")
    return arr / norm


def _coerce_vector(vector: Sequence[float], *, name: str) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    if arr.shape != (3,):
        raise ValueError(f"{name} must be a 3-component vector")
    return arr


def _coerce_nodal_field(
    field: np.ndarray,
    total_nodes: int,
    dofs_per_node: int | None = None,
) -> tuple[np.ndarray, int]:
    field = np.asarray(getattr(field, "array", field), dtype=float)

    if field.ndim == 1:
        if dofs_per_node is None:
            if field.size % total_nodes != 0:
                raise ValueError(
                    "Field size is not divisible by the number of parent nodes; "
                    "specify dofs_per_node explicitly."
                )
            dofs_per_node = field.size // total_nodes
        field = field.reshape(total_nodes, dofs_per_node)
    elif field.ndim == 2:
        if field.shape[0] != total_nodes:
            raise ValueError(f"Field has {field.shape[0]} rows, but mesh has {total_nodes} nodes.")
        if dofs_per_node is None:
            dofs_per_node = field.shape[1]
        elif field.shape[1] != dofs_per_node:
            raise ValueError(f"Field has {field.shape[1]} columns, expected {dofs_per_node}.")
    else:
        raise ValueError(f"Field must be 1D or 2D, got {field.ndim}D")

    assert dofs_per_node is not None
    return field, dofs_per_node


def discover_rotor_blade_set_names(
    mesh: MeshModel,
    prefix: str = ROTOR_BLADE_SET_PREFIX,
) -> list[str]:
    """Discover and sort aggregate rotor blade element-set names."""
    names = [name for name in mesh.element_sets_names if name.startswith(prefix)]
    if not names:
        raise ValueError(
            "No aggregate rotor blade element sets were found. "
            f"Expected names starting with '{prefix}'."
        )

    def _sort_key(name: str) -> tuple[int, int | str]:
        suffix = name[len(prefix) :]
        if suffix.isdigit():
            return (0, int(suffix))
        return (1, suffix)

    return sorted(names, key=_sort_key)


@dataclass
class RotorBladeAeroGeometry:
    """Aerodynamic geometry slice for one blade in a full rotor."""

    blade_index: int
    azimuth_deg: float
    element_set_name: str
    node_set_name: str
    mesh: MeshModel
    parent_node_indices: np.ndarray
    parent_element_indices: np.ndarray
    total_parent_nodes: int

    @property
    def node_count(self) -> int:
        return self.mesh.node_count

    @property
    def element_count(self) -> int:
        return self.mesh.elements_count

    @property
    def reference_coords(self) -> np.ndarray:
        return self.mesh.coords_array

    def extract_parent_field(
        self,
        field: np.ndarray,
        dofs_per_node: int | None = None,
    ) -> np.ndarray:
        """Extract a parent-ordered nodal field for this blade only."""
        parent_field, _ = _coerce_nodal_field(
            field,
            total_nodes=self.total_parent_nodes,
            dofs_per_node=dofs_per_node,
        )
        return parent_field[self.parent_node_indices]

    def scatter_blade_field(
        self,
        field: np.ndarray,
        dofs_per_node: int | None = None,
    ) -> np.ndarray:
        """Scatter a blade-local nodal field into the full-rotor ordering."""
        blade_field, blade_dofs = _coerce_nodal_field(
            field,
            total_nodes=self.node_count,
            dofs_per_node=dofs_per_node,
        )
        parent_field = np.zeros((self.total_parent_nodes, blade_dofs), dtype=float)
        parent_field[self.parent_node_indices] = blade_field
        return parent_field

    def get_deformed_mesh(
        self,
        displacements: np.ndarray,
        scale: float = 1.0,
        dofs_per_node: int | None = None,
        displacement_indices: tuple[int, ...] | None = None,
    ) -> MeshModel:
        """Return a deformed blade mesh extracted from a full-rotor field."""
        blade_displacements, blade_dofs = _coerce_nodal_field(
            displacements,
            total_nodes=self.total_parent_nodes,
            dofs_per_node=dofs_per_node,
        )
        blade_displacements = blade_displacements[self.parent_node_indices]
        if displacement_indices is None:
            return self.mesh.get_deformed_copy(
                blade_displacements,
                scale=scale,
                dofs_per_node=blade_dofs,
            )
        return self.mesh.get_deformed_copy(
            blade_displacements,
            scale=scale,
            dofs_per_node=blade_dofs,
            displacement_indices=displacement_indices,
        )


@dataclass
class FullRotorAeroGeometry:
    """Full-rotor aerodynamic geometry resolved blade by blade."""

    rotor_mesh: MeshModel
    blades: tuple[RotorBladeAeroGeometry, ...]
    rotation_axis: np.ndarray
    rotation_center: np.ndarray
    hub_radius: float | None = None
    azimuth_offset_deg: float = 0.0

    @property
    def n_blades(self) -> int:
        return len(self.blades)

    def get_blade(self, blade_index: int) -> RotorBladeAeroGeometry:
        return self.blades[blade_index]

    def get_deformed_blade_meshes(
        self,
        displacements: np.ndarray,
        scale: float = 1.0,
        dofs_per_node: int | None = None,
        displacement_indices: tuple[int, ...] | None = None,
    ) -> tuple[MeshModel, ...]:
        """Return deformed meshes for all blades from a full-rotor field."""
        total_parent_nodes = self.rotor_mesh.node_count
        _, blade_dofs = _coerce_nodal_field(
            displacements,
            total_nodes=total_parent_nodes,
            dofs_per_node=dofs_per_node,
        )
        return tuple(
            blade.get_deformed_mesh(
                displacements,
                scale=scale,
                dofs_per_node=blade_dofs,
                displacement_indices=displacement_indices,
            )
            for blade in self.blades
        )

    def extract_blade_fields(
        self,
        field: np.ndarray,
        dofs_per_node: int | None = None,
    ) -> tuple[np.ndarray, ...]:
        """Split a full-rotor nodal field into bladewise blocks."""
        total_parent_nodes = self.rotor_mesh.node_count
        parent_field, blade_dofs = _coerce_nodal_field(
            field,
            total_nodes=total_parent_nodes,
            dofs_per_node=dofs_per_node,
        )
        return tuple(
            blade.extract_parent_field(parent_field, dofs_per_node=blade_dofs)
            for blade in self.blades
        )

    def assemble_global_nodal_field(
        self,
        blade_fields: Sequence[np.ndarray],
        dofs_per_node: int | None = None,
    ) -> np.ndarray:
        """Assemble blade-local nodal fields into the full-rotor ordering."""
        if len(blade_fields) != self.n_blades:
            raise ValueError(f"Expected {self.n_blades} blade fields, got {len(blade_fields)}.")

        assembled: np.ndarray | None = None
        field_width = dofs_per_node

        for blade, blade_field in zip(self.blades, blade_fields):
            scattered = blade.scatter_blade_field(blade_field, dofs_per_node=field_width)
            field_width = scattered.shape[1]
            if assembled is None:
                assembled = scattered
            else:
                assembled += scattered

        if assembled is None:
            width = 1 if dofs_per_node is None else dofs_per_node
            return np.zeros((self.rotor_mesh.node_count, width), dtype=float)

        return assembled


def build_full_rotor_aero_geometry(
    mesh: MeshModel,
    blade_set_names: Iterable[str] | None = None,
    rotation_axis: Sequence[float] = (0.0, 1.0, 0.0),
    rotation_center: Sequence[float] = (0.0, 0.0, 0.0),
    hub_radius: float | None = None,
    azimuth_offset_deg: float = 0.0,
) -> FullRotorAeroGeometry:
    """Build a blade-resolved aerodynamic view of a full rotor mesh."""
    if hub_radius is not None and hub_radius < 0.0:
        raise ValueError(f"hub_radius must be non-negative, got {hub_radius}")

    if blade_set_names is None:
        ordered_blade_sets = discover_rotor_blade_set_names(mesh)
    else:
        ordered_blade_sets = list(blade_set_names)

    if not ordered_blade_sets:
        raise ValueError("At least one blade set is required to build full-rotor geometry")

    n_blades = len(ordered_blade_sets)
    rotor_coords = mesh.coords_array
    blades: list[RotorBladeAeroGeometry] = []

    for blade_index, set_name in enumerate(ordered_blade_sets):
        if set_name not in mesh.element_sets:
            raise ValueError(f"Element set '{set_name}' not found in the rotor mesh")

        node_set_name = set_name
        selected_element_ids = mesh.element_sets[set_name].element_ids
        if node_set_name in mesh.node_sets:
            parent_node_ids = mesh.node_sets[node_set_name].node_ids
        else:
            parent_node_ids = set()
            for element in mesh.element_sets[set_name].elements:
                parent_node_ids.update(element.node_ids)

        parent_element_indices = np.array(
            [
                mesh.get_element_index(element.id)
                for element in mesh.elements
                if element.id in selected_element_ids
            ],
            dtype=int,
        )

        parent_node_indices = np.array(
            sorted(mesh.get_node_index(node_id) for node_id in parent_node_ids),
            dtype=int,
        )

        blade_mesh = mesh.extract_submesh(set_name)
        expected_coords = rotor_coords[parent_node_indices]
        if blade_mesh.node_count != len(parent_node_indices):
            raise RuntimeError(
                f"Blade submesh '{set_name}' has inconsistent node count: "
                f"{blade_mesh.node_count} vs {len(parent_node_indices)}"
            )
        if not np.allclose(blade_mesh.coords_array, expected_coords):
            raise RuntimeError(f"Blade submesh '{set_name}' does not preserve parent node ordering")

        azimuth_deg = (azimuth_offset_deg + blade_index * 360.0 / n_blades) % 360.0
        blades.append(
            RotorBladeAeroGeometry(
                blade_index=blade_index,
                azimuth_deg=azimuth_deg,
                element_set_name=set_name,
                node_set_name=node_set_name,
                mesh=blade_mesh,
                parent_node_indices=parent_node_indices,
                parent_element_indices=parent_element_indices,
                total_parent_nodes=mesh.node_count,
            )
        )

    return FullRotorAeroGeometry(
        rotor_mesh=mesh,
        blades=tuple(blades),
        rotation_axis=_normalize_vector(rotation_axis),
        rotation_center=_coerce_vector(rotation_center, name="rotation_center"),
        hub_radius=None if hub_radius is None else float(hub_radius),
        azimuth_offset_deg=azimuth_offset_deg,
    )
