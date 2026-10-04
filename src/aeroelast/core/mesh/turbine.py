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
* The **tower is the displaced component**: vertical along ``+Z``, with its
  top at ``tower_offset`` relative to the rotor centre.  When the YAML carries
  the data the default offset is ``+overhang * rotor_axis`` (the tower sits
  **behind** the rotor plane and rises to meet the horizontal nacelle);
  ``distance_tt_hub`` adds an optional vertical drop.  ``tower_offset`` /
  ``overhang`` / ``distance_tt_hub`` / ``tower_base_z`` are all parameters, so a
  file with no nacelle data can still be placed.
* The **hub and the nacelle are a single body** (:class:`NacelleMesh`): a
  **hemispherical hub centred on the rotor origin**, a constant-radius cylinder
  extending **backwards along the rotor axis** (always horizontal) past the
  tower, and a **rear hemisphere** at the tail.  Its length defaults to
  ``2 * overhang + radius``, so the body is symmetric about the tower axis and
  the tail clears the tower by more than the tower diameter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from aeroelast.core.mesh.components import (
    NacelleMesh,
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
    nacelle: Optional[MeshModel] = None
    tower: Optional[MeshModel] = None
    rotor_axis: np.ndarray = field(default_factory=lambda: np.array(_DEFAULT_ROTOR_AXIS))
    nacelle_radius: float = 0.0
    tower_top: np.ndarray = field(default_factory=lambda: np.zeros(3))

    @property
    def meshes(self) -> dict[str, MeshModel]:
        """All component meshes: blades first, then the hub/nacelle body, tower."""
        ordered: dict[str, MeshModel] = dict(self.blade_meshes)
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
        nacelle_n_circ: int = 32,
        nacelle_n_axial: int = 12,
        nacelle_n_tip: int = 6,
        nacelle_n_tail: int = 6,
        nacelle_rear_tip: bool = True,
        nacelle_diameter: Optional[float] = None,
        nacelle_length: Optional[float] = None,
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

        self.nacelle_n_circ = nacelle_n_circ
        self.nacelle_n_axial = nacelle_n_axial
        self.nacelle_n_tip = nacelle_n_tip
        self.nacelle_n_tail = nacelle_n_tail
        self.nacelle_rear_tip = nacelle_rear_tip
        self.nacelle_diameter = nacelle_diameter
        self.nacelle_length = nacelle_length

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
        """Tower-top position relative to the rotor centre (origin).

        The nacelle is always horizontal along the rotor axis, so the tower top
        sits **behind** the rotor at ``+overhang * rotor_axis`` and rises to meet
        the nacelle.  ``distance_tt_hub``, when given, adds a vertical drop.
        """
        if self.tower_offset is not None:
            return self.tower_offset.copy()

        overhang = self.overhang
        if overhang is None and definition.nacelle is not None:
            overhang = definition.nacelle.overhang
        if overhang is None:
            raise ValueError(
                "cannot place the tower: the input file has no nacelle overhang; "
                "pass tower_offset (and/or overhang)"
            )

        drop = float(self.distance_tt_hub) if self.distance_tt_hub is not None else 0.0
        return float(overhang) * self.rotor_axis + np.array([0.0, 0.0, -drop])

    def _nacelle_radius(self, definition) -> float:
        if self.nacelle_diameter is not None:
            return float(self.nacelle_diameter) / 2.0
        if self._blade_generator is not None and self._blade_generator.numad_blade is not None:
            hub_diameter = self._blade_generator.numad_blade.definition.hub_diameter
            if hub_diameter:
                return float(hub_diameter) / 2.0
        if definition.hub is not None:
            return float(definition.hub.diameter) / 2.0
        raise ValueError("cannot resolve the hub diameter; pass hub_diameter explicitly")

    def _nacelle_length(self, definition, tower_top: np.ndarray, radius: float) -> float:
        """Axial extent from the hub centre to the rear tip.

        Defaults to ``2 * |tower_top| + radius`` so the body (front pole at
        ``-radius``) is **symmetric about the tower axis**: the tail sticks out
        behind the tower by the same distance the hub sticks out in front.
        ``hub_length`` overrides it.
        """
        if self.nacelle_length is not None:
            return float(self.nacelle_length)
        return 2.0 * float(np.linalg.norm(tower_top)) + float(radius)

    # -- assembly ----------------------------------------------------------

    def generate(self, renumber: str | None = None, verbose: bool = False) -> TurbineMeshes:
        """Generate and assemble every component (rotor centred at the origin)."""
        definition = read_windio_components(self.yaml_file)

        n_blades = self.n_blades if self.n_blades is not None else definition.number_of_blades
        if n_blades < 1:
            raise ValueError(f"n_blades must be >= 1, got {n_blades}")

        base_blade = self._build_blade(verbose=verbose)
        nacelle_radius = self._nacelle_radius(definition)

        # --- tower: the displaced component -------------------------------
        tower_top = self._resolve_tower_top(definition)
        tower_builder = self._tower_builder()
        source_top = tower_builder.stations[-1, :3]
        shift = tower_top - source_top
        if self.tower_base_z is not None:
            shift[2] += float(self.tower_base_z) - (tower_builder.stations[0, 2] + shift[2])
        tower_mesh = _transform_mesh(tower_builder.generate(), np.eye(3), shift)

        # --- hub + nacelle: one horizontal body, hub sphere on the rotor ----
        body_length = self._nacelle_length(definition, tower_top, nacelle_radius)
        if body_length <= 0.0:
            raise ValueError("hub/nacelle body length must be positive")
        nacelle_mesh = NacelleMesh(
            length=body_length,
            radius=nacelle_radius,
            n_circ=self.nacelle_n_circ,
            n_axial=self.nacelle_n_axial,
            n_tip=self.nacelle_n_tip,
            n_tail=self.nacelle_n_tail,
            center=(0.0, 0.0, 0.0),
            axis=self.rotor_axis,
            rear_tip=self.nacelle_rear_tip,
        ).generate()

        # --- blades -------------------------------------------------------
        cone_angle = float(definition.hub.cone_angle) if definition.hub is not None else 0.0

        # The base blade already spins about +Y; R0 maps that onto rotor_axis.
        # ``rotor_axis`` is used exactly as given: no hidden tilt is applied.
        # Coning tilts the blade out of the rotor plane, away from the tower
        # (the tower sits behind the rotor at +rotor_axis).
        r0 = _rotation_between((0.0, 1.0, 0.0), self.rotor_axis)
        cone = _rotation((1.0, 0.0, 0.0), cone_angle)

        blade_meshes: dict[str, MeshModel] = {}
        for index in range(n_blades):
            azimuth = _rotation((0.0, 1.0, 0.0), 2.0 * np.pi * index / n_blades)
            matrix = r0 @ azimuth @ cone
            offset = matrix @ np.array([0.0, 0.0, nacelle_radius])
            name = f"blade_{index + 1}"
            blade_meshes[name] = _transform_mesh(base_blade, matrix, offset, set_suffix=name)

        result = TurbineMeshes(
            blade_meshes=blade_meshes,
            nacelle=nacelle_mesh,
            tower=tower_mesh,
            rotor_axis=self.rotor_axis,
            nacelle_radius=nacelle_radius,
            tower_top=tower_top,
        )
        if renumber is not None:
            for mesh in result.meshes.values():
                mesh.renumber_mesh(algorithm=renumber)
        return result
