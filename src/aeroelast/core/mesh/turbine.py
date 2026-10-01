"""
Turbine mesh assembly: blade(s) + hub/nacelle body + tower as separate meshes.

:class:`TurbineMesh` reads one WindIO turbine YAML, generates every component
(the blade mesh through the numad shell mesher, the axisymmetric components
through :mod:`aeroelast.core.mesh.components`) and places them in a single
turbine coordinate frame.

Coordinate convention
---------------------
* The **rotor is centred at the origin** ``(0, 0, 0)``.
* The **rotor axis is configurable** (``rotor_axis``, default **+Y**, the same
  convention :class:`~aeroelast.core.mesh.generators.RotorMesh` uses).  The
  rotor plane is therefore the ``X-Z`` plane by default.
* The blade base mesh already has span **+Z**, chord **+X** and thickness
  **+Y**, so with ``rotor_axis=+Y`` no remap is needed.  For any other axis the
  whole blade is rotated by the minimal rotation taking **+Y** to ``rotor_axis``.
* The **tower is the displaced component**: it is vertical along ``+Z`` with its
  top placed at ``tower_offset`` relative to the rotor centre.  When the YAML
  carries the data, the default offset is ``-overhang * rotor_axis -
  distance_tt_hub * Z`` (yaw axis behind the rotor plane, tower top below the
  hub); otherwise ``tower_offset`` and/or ``tower_base_z`` must be given.
* The **hub and the nacelle are a single body** (:class:`HubNacelleMesh`): a
  constant-radius cylinder with a hemispherical tip, running from the tower top
  to the rotor centre.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from aeroelast.core.mesh.components import (
    HubNacelleMesh,
    TowerMesh,
    read_windio_components,
)
from aeroelast.core.mesh.entities import ElementSet, MeshElement, Node, NodeSet
from aeroelast.core.mesh.generators import BladeMesh
from aeroelast.core.mesh.io.writers import write_mesh
from aeroelast.core.mesh.model import MeshModel

_DEFAULT_ROTOR_AXIS = (0.0, 1.0, 0.0)


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


def _rotation_between(
    source: np.ndarray | Sequence[float], target: np.ndarray | Sequence[float]
) -> np.ndarray:
    """Minimal rotation taking ``source`` onto ``target``."""
    a = _unit(source)
    b = _unit(target)
    cross = np.cross(a, b)
    sin_angle = float(np.linalg.norm(cross))
    cos_angle = float(np.dot(a, b))
    if sin_angle <= 1e-12:
        if cos_angle > 0.0:
            return np.eye(3)
        helper = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        return _rotation(_unit(np.cross(a, helper)), np.pi)
    return _rotation(cross / sin_angle, float(np.arctan2(sin_angle, cos_angle)))


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
    """Separate component meshes in one turbine frame (rotor centred at origin)."""

    blade_meshes: dict[str, MeshModel] = field(default_factory=dict)
    hub: Optional[MeshModel] = None
    tower: Optional[MeshModel] = None
    rotor_axis: np.ndarray = field(default_factory=lambda: np.array(_DEFAULT_ROTOR_AXIS))
    hub_radius: float = 0.0
    tower_top: np.ndarray = field(default_factory=lambda: np.zeros(3))

    @property
    def meshes(self) -> dict[str, MeshModel]:
        """All component meshes: blades first, then the hub/nacelle body, tower."""
        ordered: dict[str, MeshModel] = dict(self.blade_meshes)
        if self.hub is not None:
            ordered["hub"] = self.hub
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
    """Generate a full turbine as separate blade/hub/tower meshes."""

    def __init__(
        self,
        yaml_file: str,
        n_blades: Optional[int] = None,
        element_size: float = 0.5,
        n_samples: int = 300,
        include_webs: bool = False,
        span_grading: str = "chord",
        rotor_axis: Sequence[float] = _DEFAULT_ROTOR_AXIS,
        tower_element_size: Optional[float] = None,
        tower_n_axial: int = 20,
        tower_n_circ: int = 32,
        tower_height: Optional[float] = None,
        tower_base_diameter: Optional[float] = None,
        tower_top_diameter: Optional[float] = None,
        hub_n_circ: int = 32,
        hub_n_axial: int = 12,
        hub_n_tip: int = 6,
        hub_diameter: Optional[float] = None,
        hub_length: Optional[float] = None,
        tower_offset: Optional[Sequence[float]] = None,
        tower_base_z: Optional[float] = None,
        overhang: Optional[float] = None,
        distance_tt_hub: Optional[float] = None,
    ) -> None:
        if (tower_height is None) != (tower_base_diameter is None):
            raise ValueError("tower_height and tower_base_diameter must be given together")

        self.yaml_file = yaml_file
        self.n_blades = n_blades
        self.element_size = element_size
        self.n_samples = n_samples
        self.include_webs = include_webs
        self.span_grading = span_grading
        self.rotor_axis = _unit(rotor_axis)

        self.tower_element_size = tower_element_size
        self.tower_n_axial = tower_n_axial
        self.tower_n_circ = tower_n_circ
        self.tower_height = tower_height
        self.tower_base_diameter = tower_base_diameter
        self.tower_top_diameter = tower_top_diameter

        self.hub_n_circ = hub_n_circ
        self.hub_n_axial = hub_n_axial
        self.hub_n_tip = hub_n_tip
        self.hub_diameter = hub_diameter
        self.hub_length = hub_length

        self.tower_offset = (
            np.asarray(tower_offset, dtype=float).reshape(3) if tower_offset is not None else None
        )
        self.tower_base_z = tower_base_z
        self.overhang = overhang
        self.distance_tt_hub = distance_tt_hub

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

    def _resolve_tower_top(self, definition) -> np.ndarray:
        """Tower-top position relative to the rotor centre (origin)."""
        if self.tower_offset is not None:
            return self.tower_offset.copy()

        overhang = self.overhang
        if overhang is None and definition.nacelle is not None:
            overhang = definition.nacelle.overhang
        distance_tt_hub = self.distance_tt_hub
        if distance_tt_hub is None and definition.nacelle is not None:
            distance_tt_hub = definition.nacelle.distance_tt_hub

        if overhang is None and distance_tt_hub is None:
            raise ValueError(
                "cannot place the tower: the input file has no nacelle overhang/"
                "distance_tt_hub; pass tower_offset (and/or overhang, distance_tt_hub)"
            )
        return np.array(
            [
                -float(overhang or 0.0) * self.rotor_axis[0],
                -float(overhang or 0.0) * self.rotor_axis[1],
                -float(distance_tt_hub or 0.0),
            ]
        )

    def _body_radius(self, definition) -> float:
        if self.hub_diameter is not None:
            return float(self.hub_diameter) / 2.0
        if self._blade_generator is not None and self._blade_generator.numad_blade is not None:
            hub_diameter = self._blade_generator.numad_blade.definition.hub_diameter
            if hub_diameter:
                return float(hub_diameter) / 2.0
        if definition.hub is not None:
            return float(definition.hub.diameter) / 2.0
        raise ValueError("cannot resolve the hub diameter; pass hub_diameter explicitly")

    def _body_length(self, definition, tower_top: np.ndarray) -> float:
        if self.hub_length is not None:
            return float(self.hub_length)
        return float(np.linalg.norm(tower_top))

    # -- assembly ----------------------------------------------------------

    def generate(self, renumber: str | None = None, verbose: bool = False) -> TurbineMeshes:
        """Generate and assemble every component (rotor centred at the origin)."""
        definition = read_windio_components(self.yaml_file)

        n_blades = self.n_blades if self.n_blades is not None else definition.number_of_blades
        if n_blades < 1:
            raise ValueError(f"n_blades must be >= 1, got {n_blades}")

        base_blade = self._build_blade(verbose=verbose)
        hub_radius = self._body_radius(definition)

        # --- tower: the displaced component -------------------------------
        tower_top = self._resolve_tower_top(definition)
        tower_builder = self._tower_builder()
        source_top = tower_builder.stations[-1, :3]
        shift = tower_top - source_top
        if self.tower_base_z is not None:
            shift[2] += float(self.tower_base_z) - (tower_builder.stations[0, 2] + shift[2])
        tower_mesh = _transform_mesh(tower_builder.generate(), np.eye(3), shift)

        # --- hub + nacelle: one constant-radius body with a round tip -----
        body_span = -tower_top  # from the tower top to the rotor centre
        body_length = self._body_length(definition, body_span)
        if body_length <= 0.0:
            raise ValueError("hub/nacelle body length must be positive")
        body_axis = _unit(body_span) if np.linalg.norm(body_span) > 0.0 else self.rotor_axis
        hub_mesh = HubNacelleMesh(
            length=body_length,
            radius=hub_radius,
            n_circ=self.hub_n_circ,
            n_axial=self.hub_n_axial,
            n_tip=self.hub_n_tip,
            center=tower_top,
            axis=body_axis,
        ).generate()

        # --- blades -------------------------------------------------------
        cone_angle = float(definition.hub.cone_angle) if definition.hub is not None else 0.0

        # The base blade already spins about +Y; R0 maps that onto rotor_axis.
        # ``rotor_axis`` is used exactly as given: no hidden tilt is applied.
        # Coning tilts the blade out of the rotor plane away from the tower.
        r0 = _rotation_between((0.0, 1.0, 0.0), self.rotor_axis)
        cone = _rotation((1.0, 0.0, 0.0), -cone_angle)

        blade_meshes: dict[str, MeshModel] = {}
        for index in range(n_blades):
            azimuth = _rotation((0.0, 1.0, 0.0), 2.0 * np.pi * index / n_blades)
            matrix = r0 @ azimuth @ cone
            offset = matrix @ np.array([0.0, 0.0, hub_radius])
            name = f"blade_{index + 1}"
            blade_meshes[name] = _transform_mesh(base_blade, matrix, offset, set_suffix=name)

        result = TurbineMeshes(
            blade_meshes=blade_meshes,
            hub=hub_mesh,
            tower=tower_mesh,
            rotor_axis=self.rotor_axis,
            hub_radius=hub_radius,
            tower_top=tower_top,
        )
        if renumber is not None:
            for mesh in result.meshes.values():
                mesh.renumber_mesh(algorithm=renumber)
        return result
