"""
Turbine mesh assembly: blades + hub/nacelle + tower as separate meshes.

:class:`TurbineMesh` reads one WindIO turbine YAML, generates every component
(the blade mesh through the numad shell mesher, the axisymmetric components
through :mod:`aeroelast.core.mesh.components`) and places them in a single
turbine coordinate frame.  :class:`TurbineMeshes` keeps the components as
separate :class:`~aeroelast.core.mesh.model.MeshModel` objects and can write
them all to a directory.

Coordinate convention
---------------------
* Tower base at ``z = tower_base_z`` (default ``0``), tower axis **+Z**.
* Hub centre at ``(overhang, 0, hub_height + base_offset)``, where
  ``base_offset`` is the rigid shift that moves the imported tower base to
  ``tower_base_z``.
* Rotor axis tilted out of **+X** by the nacelle uptilt.
* Blade span is the base-blade **+Z** direction; each blade is coned by the hub
  cone angle, azimuthally distributed about the rotor axis and rooted at the
  hub surface (``hub_diameter / 2``).

The base blade mesh has span **+Z**, chord **+X** and thickness **+Y** (the
convention :class:`~aeroelast.core.mesh.generators.RotorMesh` rotates about Y).
The remap **Rz(-90 deg)** turns that into span **+Z**, rotor axis **+X**, chord
**-Y**, which is the frame used here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from aeroelast.core.mesh.components import (
    HubMesh,
    NacelleMesh,
    TowerMesh,
    read_windio_components,
)
from aeroelast.core.mesh.entities import ElementSet, MeshElement, Node, NodeSet
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.core.mesh.io.writers import write_mesh
from aeroelast.core.mesh.model import MeshModel


def _unit(vector: np.ndarray | Sequence[float]) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(3)
    norm = float(np.linalg.norm(arr))
    if norm <= 0.0:
        raise ValueError("direction must be a non-zero vector")
    return arr / norm


def _rotation(axis: np.ndarray | Sequence[float], angle: float) -> np.ndarray:
    """Right-handed rotation matrix of ``angle`` radians about ``axis``."""
    k = _unit(axis)
    cross = np.array(
        [
            [0.0, -k[2], k[1]],
            [k[2], 0.0, -k[0]],
            [-k[1], k[0], 0.0],
        ]
    )
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * (cross @ cross)


def _copy_mesh(mesh: MeshModel) -> tuple[MeshModel, dict[int, Node], dict[int, MeshElement]]:
    """Deep-copy a mesh, returning the copy and its node/element id maps."""
    clone = MeshModel()
    node_map: dict[int, Node] = {}
    for node in mesh.nodes:
        new_node = Node(node.coords.copy())
        clone.add_node(new_node)
        node_map[node.id] = new_node

    element_map: dict[int, MeshElement] = {}
    for element in mesh.elements:
        new_element = MeshElement(
            [node_map[node.id] for node in element.nodes], element.element_type
        )
        clone.add_element(new_element)
        element_map[element.id] = new_element

    for name, node_set in mesh.node_sets.items():
        clone.add_node_set(NodeSet(name, {node_map[nid] for nid in node_set.node_ids}))
    for name, element_set in mesh.element_sets.items():
        clone.add_element_set(
            ElementSet(name, {element_map[eid] for eid in element_set.element_ids})
        )
    return clone, node_map, element_map


def _transform_mesh(
    mesh: MeshModel,
    matrix: np.ndarray,
    offset: np.ndarray,
    set_suffix: Optional[str] = None,
) -> MeshModel:
    """Return a transformed copy: ``x -> matrix @ x + offset``.

    With ``set_suffix`` every node/element set is renamed ``<name>_<suffix>``.
    """
    clone, _node_map, _element_map = _copy_mesh(mesh)

    for node in clone.nodes:
        coords = matrix @ node.coords + offset
        node.coords = coords
        node.x, node.y, node.z = (float(v) for v in coords)

    if set_suffix:
        clone.node_sets = {
            f"{name}_{set_suffix}": node_set for name, node_set in clone.node_sets.items()
        }
        clone.element_sets = {
            f"{name}_{set_suffix}": element_set for name, element_set in clone.element_sets.items()
        }

    return clone


@dataclass
class TurbineMeshes:
    """Separate component meshes in one turbine frame."""

    blade_meshes: dict[str, MeshModel] = field(default_factory=dict)
    hub: Optional[MeshModel] = None
    nacelle: Optional[MeshModel] = None
    tower: Optional[MeshModel] = None
    hub_center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    rotor_axis: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.0, 0.0]))
    hub_radius: float = 0.0
    tower_top: np.ndarray = field(default_factory=lambda: np.zeros(3))

    @property
    def meshes(self) -> dict[str, MeshModel]:
        """All component meshes: blades first, then hub, nacelle and tower."""
        ordered: dict[str, MeshModel] = dict(self.blade_meshes)
        if self.hub is not None:
            ordered["hub"] = self.hub
        if self.nacelle is not None:
            ordered["nacelle"] = self.nacelle
        if self.tower is not None:
            ordered["tower"] = self.tower
        return ordered

    def write(self, directory: str | Path, format: str = "stl", prefix: str = "") -> dict[str, str]:
        """Write every component to ``directory``, returning ``{name: path}``."""
        out_dir = Path(directory)
        out_dir.mkdir(parents=True, exist_ok=True)
        suffix = format if format.startswith(".") else f".{format}"
        written: dict[str, str] = {}
        for name, mesh in self.meshes.items():
            path = out_dir / f"{prefix}{name}{suffix}"
            write_mesh(mesh, str(path))
            written[name] = str(path)
        return written


class TurbineMesh:
    """Generate a full turbine as separate blade/hub/nacelle/tower meshes."""

    def __init__(
        self,
        yaml_file: str,
        n_blades: Optional[int] = None,
        element_size: float = 0.5,
        n_samples: int = 300,
        include_webs: bool = False,
        span_grading: str = "chord",
        tower_element_size: Optional[float] = None,
        tower_n_axial: int = 20,
        tower_n_circ: int = 32,
        hub_n_circ: int = 32,
        hub_n_merid: int = 16,
        nacelle_n_circ: int = 32,
        nacelle_n_axial: int = 12,
        nacelle_n_nose: int = 6,
        tower_base_z: float = 0.0,
        hub_height: Optional[float] = None,
        overhang: Optional[float] = None,
        hub_diameter: Optional[float] = None,
        tower_height: Optional[float] = None,
        tower_base_diameter: Optional[float] = None,
        tower_top_diameter: Optional[float] = None,
        nacelle_length: Optional[float] = None,
        nacelle_body_diameter: Optional[float] = None,
        nacelle_nose_diameter: Optional[float] = None,
    ) -> None:
        if tower_height is not None and tower_base_diameter is None:
            raise ValueError("tower_height requires tower_base_diameter")
        if tower_base_diameter is not None and tower_height is None:
            raise ValueError("tower_base_diameter requires tower_height")
        if nacelle_body_diameter is not None and nacelle_length is None:
            raise ValueError("nacelle_body_diameter requires nacelle_length")

        self.yaml_file = yaml_file
        self.n_blades = n_blades
        self.element_size = element_size
        self.n_samples = n_samples
        self.include_webs = include_webs
        self.span_grading = span_grading

        self.tower_element_size = tower_element_size
        self.tower_n_axial = tower_n_axial
        self.tower_n_circ = tower_n_circ
        self.hub_n_circ = hub_n_circ
        self.hub_n_merid = hub_n_merid
        self.nacelle_n_circ = nacelle_n_circ
        self.nacelle_n_axial = nacelle_n_axial
        self.nacelle_n_nose = nacelle_n_nose

        self.tower_base_z = tower_base_z
        self.hub_height = hub_height
        self.overhang = overhang
        self.hub_diameter = hub_diameter
        self.tower_height = tower_height
        self.tower_base_diameter = tower_base_diameter
        self.tower_top_diameter = tower_top_diameter
        self.nacelle_length = nacelle_length
        self.nacelle_body_diameter = nacelle_body_diameter
        self.nacelle_nose_diameter = nacelle_nose_diameter

        self._blade_generator: Optional[BladeMesh] = None

    # -- components --------------------------------------------------------

    def _tower_builder(self) -> TowerMesh:
        base_diameter = self.tower_base_diameter
        if base_diameter is not None:
            height = self.tower_height
            if height is None:
                raise ValueError("tower_height is required with tower_base_diameter")
            return TowerMesh.from_params(
                height=height,
                base_diameter=base_diameter,
                top_diameter=self.tower_top_diameter,
                n_axial=self.tower_n_axial,
                n_circ=self.tower_n_circ,
            )
        return TowerMesh.from_windio(
            self.yaml_file,
            n_circ=self.tower_n_circ,
            n_axial=self.tower_n_axial,
            element_size=self.tower_element_size,
        )

    def _build_blade(self, verbose: bool) -> MeshModel:
        self._blade_generator = BladeMesh(
            yaml_file=self.yaml_file,
            element_size=self.element_size,
            n_samples=self.n_samples,
            include_webs=self.include_webs,
            span_grading=self.span_grading,
        )
        return self._blade_generator.generate(verbose=verbose)

    def _build_hub(self) -> MeshModel:
        if self.hub_diameter is not None:
            hub = HubMesh(self.hub_diameter, n_circ=self.hub_n_circ, n_merid=self.hub_n_merid)
        else:
            hub = HubMesh.from_windio(
                self.yaml_file, n_circ=self.hub_n_circ, n_merid=self.hub_n_merid
            )
        return hub.generate()

    # -- assembly ----------------------------------------------------------

    def generate(self, renumber: str | None = None, verbose: bool = False) -> TurbineMeshes:
        """Generate and assemble every component."""
        definition = read_windio_components(self.yaml_file)

        n_blades = self.n_blades if self.n_blades is not None else definition.number_of_blades
        if n_blades < 1:
            raise ValueError(f"n_blades must be >= 1, got {n_blades}")

        base_blade = self._build_blade(verbose=verbose)
        hub_radius = self._hub_radius(definition)

        tower_builder = self._tower_builder()
        tower_stations = tower_builder.stations
        base_offset = float(self.tower_base_z) - float(tower_stations[0, 2])
        shift = np.array([0.0, 0.0, base_offset])

        tower_mesh = _transform_mesh(tower_builder.generate(), np.eye(3), shift)
        tower_top = tower_stations[-1, :3] + shift

        hub_height = self._resolve_hub_height(definition, tower_top, base_offset)
        overhang = self._resolve_overhang(definition)
        hub_center = np.array([overhang, 0.0, hub_height])

        uptilt = float(definition.nacelle.uptilt) if definition.nacelle is not None else 0.0
        rotor_axis = _rotation((0.0, 1.0, 0.0), -uptilt) @ np.array([1.0, 0.0, 0.0])

        hub_mesh = _transform_mesh(self._build_hub(), np.eye(3), hub_center)

        nacelle_mesh = self._assemble_nacelle(definition, tower_top, hub_center, hub_radius)

        cone_angle = float(definition.hub.cone_angle) if definition.hub is not None else 0.0
        remap = _rotation((0.0, 0.0, 1.0), -np.pi / 2.0)
        cone = _rotation((0.0, 1.0, 0.0), cone_angle)
        tilt = _rotation((0.0, 1.0, 0.0), -uptilt)

        blade_meshes: dict[str, MeshModel] = {}
        for index in range(n_blades):
            azimuth = _rotation((1.0, 0.0, 0.0), 2.0 * np.pi * index / n_blades)
            matrix = tilt @ azimuth @ cone @ remap
            offset = hub_center + matrix @ np.array([0.0, 0.0, hub_radius])
            name = f"blade_{index + 1}"
            blade_meshes[name] = _transform_mesh(base_blade, matrix, offset, set_suffix=name)

        result = TurbineMeshes(
            blade_meshes=blade_meshes,
            hub=hub_mesh,
            nacelle=nacelle_mesh,
            tower=tower_mesh,
            hub_center=hub_center,
            rotor_axis=rotor_axis,
            hub_radius=hub_radius,
            tower_top=tower_top,
        )
        if renumber is not None:
            for mesh in result.meshes.values():
                mesh.renumber_mesh(algorithm=renumber)
        return result

    def _assemble_nacelle(
        self,
        definition,
        tower_top: np.ndarray,
        hub_center: np.ndarray,
        hub_radius: float,
    ) -> Optional[MeshModel]:
        if self.nacelle_body_diameter is not None:
            body_diameter = self.nacelle_body_diameter
            nose_diameter = self.nacelle_nose_diameter
        elif definition.nacelle is not None:
            body_diameter = float(definition.nacelle.body_diameter)
            nose_diameter = float(definition.nacelle.nose_diameter)
        else:
            return None

        # Span the nacelle from the yaw axis (tower top) to the hub.
        span = hub_center - tower_top
        length = float(np.linalg.norm(span))
        if length <= 0.0:
            return None
        return NacelleMesh(
            length=length,
            body_diameter=body_diameter,
            nose_diameter=nose_diameter,
            n_circ=self.nacelle_n_circ,
            n_axial=self.nacelle_n_axial,
            n_nose=self.nacelle_n_nose,
            center=tower_top,
            axis=span / length,
        ).generate()

    # -- helpers -----------------------------------------------------------

    def _hub_radius(self, definition) -> float:
        if self.hub_diameter is not None:
            return float(self.hub_diameter) / 2.0
        if self._blade_generator is not None and self._blade_generator.numad_blade is not None:
            hub_diameter = self._blade_generator.numad_blade.definition.hub_diameter
            if hub_diameter:
                return float(hub_diameter) / 2.0
        if definition.hub is not None:
            return float(definition.hub.diameter) / 2.0
        raise ValueError("cannot resolve hub diameter; pass hub_diameter explicitly")

    def _resolve_hub_height(self, definition, tower_top: np.ndarray, base_offset: float) -> float:
        if self.hub_height is not None:
            return float(self.hub_height) + base_offset
        if definition.hub_height is not None:
            return float(definition.hub_height) + base_offset
        if definition.nacelle is not None and definition.nacelle.distance_tt_hub:
            return float(tower_top[2]) + float(definition.nacelle.distance_tt_hub)
        raise ValueError("cannot resolve hub height; pass hub_height explicitly")

    def _resolve_overhang(self, definition) -> float:
        if self.overhang is not None:
            return float(self.overhang)
        if definition.nacelle is not None and definition.nacelle.overhang:
            return float(definition.nacelle.overhang)
        return 0.0
