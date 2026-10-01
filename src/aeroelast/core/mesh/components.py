"""
Axisymmetric turbine component meshes: tower and hub/nacelle body.

These components are *surface* meshes (2D elements) suitable for CFD walls and
for STL / VTK / OBJ export.  Every mesh is generated analytically from a
meridian profile revolved around an axis, so the element counts are exact and
deterministic — no mesher call is involved.

WindIO input
------------
:func:`read_windio_components` reads ``components.tower``, ``components.hub``
and ``components.nacelle`` (plus ``assembly``) from a WindIO turbine YAML.
``TowerMesh.from_windio`` / ``HubNacelleMesh.from_windio`` use that definition;
``from_params`` is the fallback when the input file does not carry the
component.

Coordinate convention for a standalone component
------------------------------------------------
Each generator is built along ``axis`` through ``center`` in its own local
frame.  :class:`~aeroelast.core.mesh.turbine.TurbineMesh` is responsible for
placing them in the shared turbine frame (tower ``+Z``, configurable rotor
axis, rotor centred at the origin).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import yaml

from aeroelast.core.mesh.entities import (
    ElementSet,
    ElementType,
    MeshElement,
    Node,
    NodeSet,
)
from aeroelast.core.mesh.model import MeshModel

DEFAULT_N_CIRC = 32

# Radius below which a station is treated as a degenerate point (pole).
_POLE_TOL = 1e-9


# ============================================================================
# Generic revolved-shell builder
# ============================================================================


def _unit(vector: np.ndarray | Sequence[float]) -> np.ndarray:
    arr = np.asarray(vector, dtype=float).reshape(3)
    norm = float(np.linalg.norm(arr))
    if norm <= 0.0:
        raise ValueError("axis must be a non-zero vector")
    return arr / norm


def _frame(axis: np.ndarray | Sequence[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return an orthonormal frame ``(u, v, w)`` with ``w`` along ``axis``.

    ``u`` and ``v`` span the plane normal to ``w`` and are chosen so that a
    circle traversed as ``cos(theta) * u + sin(theta) * v`` winds
    counter-clockwise when looking down ``-w``.
    """
    w = _unit(axis)
    ref = np.array([0.0, 0.0, 1.0])
    if abs(float(w @ ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    u = np.cross(ref, w)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return u, v, w


def build_revolved_shell(
    centers: np.ndarray | Sequence[Sequence[float]],
    radii: np.ndarray | Sequence[float],
    n_circ: int = DEFAULT_N_CIRC,
    axis: np.ndarray | Sequence[float] = (0.0, 0.0, 1.0),
    cap_start: bool = False,
    cap_end: bool = False,
    start_name: str = "start",
    end_name: str = "end",
) -> MeshModel:
    """Revolve a meridian profile into a closed-oriented surface mesh.

    Parameters
    ----------
    centers:
        ``(n_stations, 3)`` centreline points, ordered along ``axis``.
    radii:
        ``(n_stations,)`` radius at each station.  A radius of zero is a
        degenerate pole: it becomes a single node and the adjacent faces
        become triangles.
    n_circ:
        Number of circumferential elements.  Must be at least 3.
    axis:
        Revolution axis direction in the frame of ``centers``.
    cap_start / cap_end:
        Close the open ends with a flat fan of triangles.
    start_name / end_name:
        Names of the node sets holding the first and last station.

    Returns
    -------
    MeshModel
        Element sets ``surface`` (and ``cap_start`` / ``cap_end`` when capped)
        and node sets ``start_name``, ``end_name`` and ``all``.
    """
    centers_arr = np.asarray(centers, dtype=float).reshape(-1, 3)
    radii_arr = np.asarray(radii, dtype=float).reshape(-1)
    if centers_arr.shape[0] != radii_arr.shape[0]:
        raise ValueError(
            f"centers and radii must have the same length, got "
            f"{centers_arr.shape[0]} and {radii_arr.shape[0]}"
        )
    if centers_arr.shape[0] < 2:
        raise ValueError("a revolved shell needs at least two stations")
    if n_circ < 3:
        raise ValueError(f"n_circ must be >= 3, got {n_circ}")
    if np.any(radii_arr < 0.0):
        raise ValueError("radii must be non-negative")

    u, v, w = _frame(axis)
    theta = 2.0 * np.pi * np.arange(n_circ) / n_circ
    circle = np.outer(np.cos(theta), u) + np.outer(np.sin(theta), v)

    mesh = MeshModel()

    # Stations: a ring of n_circ nodes, or a single pole node.
    rings: list[list[Node] | Node] = []
    for center, radius in zip(centers_arr, radii_arr, strict=True):
        if radius <= _POLE_TOL:
            pole = Node(center)
            mesh.add_node(pole)
            rings.append(pole)
        else:
            ring = [Node(center + radius * offset) for offset in circle]
            for node in ring:
                mesh.add_node(node)
            rings.append(ring)

    surface: list[MeshElement] = []
    for lower, upper in zip(rings[:-1], rings[1:], strict=True):
        if isinstance(lower, list) and isinstance(upper, list):
            for k in range(n_circ):
                k2 = (k + 1) % n_circ
                elem = MeshElement([lower[k], lower[k2], upper[k2], upper[k]], ElementType.quad)
                surface.append(elem)
        elif isinstance(lower, Node) and isinstance(upper, list):
            for k in range(n_circ):
                k2 = (k + 1) % n_circ
                surface.append(MeshElement([lower, upper[k2], upper[k]], ElementType.triangle))
        elif isinstance(lower, list) and isinstance(upper, Node):
            for k in range(n_circ):
                k2 = (k + 1) % n_circ
                surface.append(MeshElement([lower[k], lower[k2], upper], ElementType.triangle))
        else:
            raise ValueError("two consecutive degenerate stations do not form a surface")
    for elem in surface:
        mesh.add_element(elem)
    mesh.add_element_set(ElementSet("surface", set(surface)))

    caps: dict[str, list[MeshElement]] = {}
    if cap_start:
        ring = rings[0]
        if not isinstance(ring, list):
            raise ValueError("cannot cap a degenerate start station")
        cap_center = Node(centers_arr[0])
        mesh.add_node(cap_center)
        caps["cap_start"] = [
            MeshElement([cap_center, ring[(k + 1) % n_circ], ring[k]], ElementType.triangle)
            for k in range(n_circ)
        ]
    if cap_end:
        ring = rings[-1]
        if not isinstance(ring, list):
            raise ValueError("cannot cap a degenerate end station")
        cap_center = Node(centers_arr[-1])
        mesh.add_node(cap_center)
        caps["cap_end"] = [
            MeshElement([cap_center, ring[k], ring[(k + 1) % n_circ]], ElementType.triangle)
            for k in range(n_circ)
        ]
    for name, elems in caps.items():
        for elem in elems:
            mesh.add_element(elem)
        mesh.add_element_set(ElementSet(name, set(elems)))

    def _nodes_of(station: list[Node] | Node) -> set[Node]:
        return set(station) if isinstance(station, list) else {station}

    mesh.add_node_set(NodeSet(start_name, _nodes_of(rings[0])))
    mesh.add_node_set(NodeSet(end_name, _nodes_of(rings[-1])))
    mesh.add_node_set(NodeSet("all", set(mesh.nodes)))

    return mesh


# ============================================================================
# WindIO component definitions
# ============================================================================


@dataclass
class TowerDefinition:
    """Tower outer shape: centreline offset and outer diameter per station."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    diameter: np.ndarray

    @property
    def height(self) -> float:
        return float(self.z[-1] - self.z[0])


@dataclass
class HubDefinition:
    """Hub outer shape: sphere diameter and rotor cone angle (radians)."""

    diameter: float
    cone_angle: float = 0.0


@dataclass
class NacelleDefinition:
    """Nacelle envelope derived from the drivetrain block."""

    body_diameter: float
    nose_diameter: float
    length: float
    uptilt: float = 0.0
    overhang: float = 0.0
    distance_tt_hub: float = 0.0


@dataclass
class TurbineDefinition:
    """Everything the turbine mesher needs from the imported model."""

    tower: Optional[TowerDefinition] = None
    hub: Optional[HubDefinition] = None
    nacelle: Optional[NacelleDefinition] = None
    number_of_blades: int = 3
    hub_height: Optional[float] = None
    rotor_diameter: Optional[float] = None


def load_windio(path: str | Path) -> dict:
    """Load a WindIO turbine YAML as a plain nested dictionary."""
    try:
        with open(path) as handle:
            data = yaml.safe_load(handle)
    except FileNotFoundError as exc:
        raise ValueError(f"{path}: WindIO input file not found") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a WindIO mapping at the top level")
    return data


def _read_tower(block: dict | None) -> Optional[TowerDefinition]:
    if not block:
        return None
    outer = block.get("outer_shape_bem", block)
    reference = outer.get("reference_axis")
    diameter = outer.get("outer_diameter")
    if not reference or not diameter:
        return None

    grid = np.asarray(reference["z"]["grid"], dtype=float)
    z = np.asarray(reference["z"]["values"], dtype=float)
    x = np.interp(grid, reference["x"]["grid"], reference["x"]["values"])
    y = np.interp(grid, reference["y"]["grid"], reference["y"]["values"])
    d = np.interp(grid, diameter["grid"], diameter["values"])

    order = np.argsort(z)
    return TowerDefinition(x=x[order], y=y[order], z=z[order], diameter=d[order])


def _read_hub(block: dict | None) -> Optional[HubDefinition]:
    if not block:
        return None
    outer = block.get("outer_shape_bem", block)
    diameter = outer.get("diameter", block.get("diameter"))
    if diameter is None:
        return None
    return HubDefinition(
        diameter=float(diameter),
        cone_angle=float(block.get("cone_angle", 0.0)),
    )


def _read_nacelle(block: dict | None) -> Optional[NacelleDefinition]:
    if not block:
        return None
    drivetrain = block.get("drivetrain", block)
    lss = drivetrain.get("lss_diameter", [3.0, 3.0])
    nose = drivetrain.get("nose_diameter", lss)
    # The nacelle spans from the yaw axis (tower top) forward to the hub, which
    # is what ``overhang`` measures.  ``distance_tt_hub`` is the vertical
    # tower-top-to-hub distance, kept separately for the assembly.
    overhang = float(drivetrain.get("overhang", drivetrain.get("distance_tt_hub", 5.0)))
    length = float(drivetrain.get("overhang", drivetrain.get("distance_tt_hub", 5.0)))
    return NacelleDefinition(
        body_diameter=float(np.mean(lss)),
        nose_diameter=float(np.mean(nose)),
        length=length,
        uptilt=float(drivetrain.get("uptilt", 0.0)),
        overhang=overhang,
        distance_tt_hub=float(drivetrain.get("distance_tt_hub", 0.0)),
    )


def read_windio_components(path: str | Path) -> TurbineDefinition:
    """Extract tower, hub, nacelle and assembly data from a WindIO YAML."""
    data = load_windio(path)
    assembly = data.get("assembly", {}) or {}
    components = data.get("components", {}) or {}

    hub_height = assembly.get("hub_height")
    rotor_diameter = assembly.get("rotor_diameter")

    return TurbineDefinition(
        tower=_read_tower(components.get("tower")),
        hub=_read_hub(components.get("hub")),
        nacelle=_read_nacelle(components.get("nacelle")),
        number_of_blades=int(assembly.get("number_of_blades", 3)),
        hub_height=float(hub_height) if hub_height is not None else None,
        rotor_diameter=float(rotor_diameter) if rotor_diameter is not None else None,
    )


# ============================================================================
# Tower
# ============================================================================


class TowerMesh:
    """Tapered tower surface mesh built from a centreline and a diameter profile.

    ``stations`` is an ``(n, 4)`` sequence of ``(x, y, z, diameter)`` items.
    Use :meth:`from_windio` / :meth:`from_params` to build it.
    """

    def __init__(
        self,
        stations: np.ndarray | Sequence[Sequence[float]],
        n_circ: int = DEFAULT_N_CIRC,
        cap_start: bool = True,
        cap_end: bool = True,
    ) -> None:
        arr = np.asarray(stations, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 4:
            raise ValueError("stations must be an (n, 4) sequence of (x, y, z, diameter)")
        if np.any(arr[:, 3] <= 0.0):
            raise ValueError("tower diameters must be positive")
        self.stations = arr
        self.n_circ = int(n_circ)
        self.cap_start = cap_start
        self.cap_end = cap_end

    @classmethod
    def from_windio(
        cls,
        path: str | Path,
        n_circ: int = DEFAULT_N_CIRC,
        n_axial: int = 20,
        element_size: Optional[float] = None,
        cap_start: bool = True,
        cap_end: bool = True,
    ) -> "TowerMesh":
        """Build the tower from ``components.tower`` in a WindIO YAML."""
        definition = read_windio_components(path).tower
        if definition is None:
            raise ValueError(
                f"{path}: no tower definition found "
                "(components.tower.outer_shape_bem); use TowerMesh.from_params instead"
            )
        return cls(
            _resample_tower(definition, n_axial, element_size),
            n_circ=n_circ,
            cap_start=cap_start,
            cap_end=cap_end,
        )

    @classmethod
    def from_params(
        cls,
        height: float,
        base_diameter: float,
        top_diameter: Optional[float] = None,
        n_axial: int = 20,
        n_circ: int = DEFAULT_N_CIRC,
        base_center: np.ndarray | Sequence[float] = (0.0, 0.0, 0.0),
        cap_start: bool = True,
        cap_end: bool = True,
    ) -> "TowerMesh":
        """Build a linear-taper tower when the input file has no tower block."""
        if height <= 0.0:
            raise ValueError("tower height must be positive")
        if base_diameter <= 0.0:
            raise ValueError("tower base_diameter must be positive")
        if top_diameter is None:
            top_diameter = 0.5 * base_diameter
        if top_diameter <= 0.0:
            raise ValueError("tower top_diameter must be positive")
        if n_axial < 1:
            raise ValueError("n_axial must be >= 1")

        z = np.linspace(0.0, float(height), int(n_axial) + 1)
        diameter = np.linspace(float(base_diameter), float(top_diameter), int(n_axial) + 1)
        center = np.asarray(base_center, dtype=float).reshape(3)
        stations = np.column_stack([center + np.outer(z, [0.0, 0.0, 1.0]), diameter])
        return cls(stations, n_circ=n_circ, cap_start=cap_start, cap_end=cap_end)

    def generate(self) -> MeshModel:
        return build_revolved_shell(
            self.stations[:, :3],
            self.stations[:, 3] / 2.0,
            n_circ=self.n_circ,
            axis=(0.0, 0.0, 1.0),
            cap_start=self.cap_start,
            cap_end=self.cap_end,
            start_name="base",
            end_name="top",
        )


def _resample_tower(
    definition: TowerDefinition, n_axial: int, element_size: Optional[float]
) -> np.ndarray:
    z = definition.z
    if element_size is not None:
        if element_size <= 0.0:
            raise ValueError("element_size must be positive")
        n_axial = max(2, int(np.ceil(definition.height / element_size)))
    if n_axial < 1:
        raise ValueError("n_axial must be >= 1")

    z_new = np.linspace(z[0], z[-1], int(n_axial) + 1)
    x_new = np.interp(z_new, z, definition.x)
    y_new = np.interp(z_new, z, definition.y)
    d_new = np.interp(z_new, z, definition.diameter)
    return np.column_stack([x_new, y_new, z_new, d_new])


# ============================================================================
# Hub + nacelle (single rotor body)
# ============================================================================


class HubNacelleMesh:
    """Single rotor hub/nacelle body.

    A constant-radius cylinder along ``axis``, closed by a flat cap at the tail
    and a hemispherical tip at the nose.  The hub and the nacelle are the same
    body, so a turbine exports exactly one ``hub`` mesh (one STL).
    """

    def __init__(
        self,
        length: float,
        radius: Optional[float] = None,
        diameter: Optional[float] = None,
        n_circ: int = DEFAULT_N_CIRC,
        n_axial: int = 12,
        n_tip: int = 6,
        center: np.ndarray | Sequence[float] = (0.0, 0.0, 0.0),
        axis: np.ndarray | Sequence[float] = (0.0, 1.0, 0.0),
        cap_tail: bool = True,
    ) -> None:
        if length <= 0.0:
            raise ValueError("hub/nacelle length must be positive")
        if radius is not None and diameter is not None:
            raise ValueError("give radius or diameter, not both")
        if radius is None:
            if diameter is None:
                raise ValueError("hub/nacelle needs either radius or diameter")
            resolved = float(diameter) / 2.0
        else:
            resolved = float(radius)
        if resolved <= 0.0:
            raise ValueError("hub/nacelle radius must be positive")
        if resolved >= length:
            raise ValueError(
                f"radius ({resolved}) must be smaller than length ({length}) "
                "so the hemispherical tip fits"
            )
        if n_tip < 1:
            raise ValueError("n_tip must be >= 1")

        self.length = float(length)
        self.radius = resolved
        self.diameter = 2.0 * resolved
        self.n_circ = int(n_circ)
        self.n_axial = int(n_axial)
        self.n_tip = int(n_tip)
        self.center = np.asarray(center, dtype=float).reshape(3)
        self.axis = _unit(axis)
        self.cap_tail = cap_tail

    @classmethod
    def from_windio(
        cls,
        path: str | Path,
        n_circ: int = DEFAULT_N_CIRC,
        n_axial: int = 12,
        n_tip: int = 6,
        length: Optional[float] = None,
        radius: Optional[float] = None,
        diameter: Optional[float] = None,
        center: np.ndarray | Sequence[float] = (0.0, 0.0, 0.0),
        axis: np.ndarray | Sequence[float] = (0.0, 1.0, 0.0),
        cap_tail: bool = True,
    ) -> "HubNacelleMesh":
        """Build the body from ``components.hub`` / ``components.nacelle``.

        The radius defaults to ``components.hub.diameter / 2`` and the length to
        the nacelle overhang (yaw axis -> hub).
        """
        definition = read_windio_components(path)
        if radius is None and diameter is None:
            if definition.hub is None:
                raise ValueError(
                    f"{path}: no hub definition found (components.hub.diameter); "
                    "use HubNacelleMesh.from_params instead"
                )
            radius = definition.hub.diameter / 2.0
        if length is None:
            if definition.nacelle is None or not definition.nacelle.length:
                raise ValueError(
                    f"{path}: no nacelle length found (components.nacelle.drivetrain); "
                    "use HubNacelleMesh.from_params instead"
                )
            length = definition.nacelle.length
        return cls(
            length=length,
            radius=radius,
            diameter=diameter,
            n_circ=n_circ,
            n_axial=n_axial,
            n_tip=n_tip,
            center=center,
            axis=axis,
            cap_tail=cap_tail,
        )

    @classmethod
    def from_params(
        cls,
        length: float,
        radius: Optional[float] = None,
        diameter: Optional[float] = None,
        n_circ: int = DEFAULT_N_CIRC,
        n_axial: int = 12,
        n_tip: int = 6,
        center: np.ndarray | Sequence[float] = (0.0, 0.0, 0.0),
        axis: np.ndarray | Sequence[float] = (0.0, 1.0, 0.0),
        cap_tail: bool = True,
    ) -> "HubNacelleMesh":
        """Build the body from explicit parameters when the YAML has none."""
        return cls(
            length=length,
            radius=radius,
            diameter=diameter,
            n_circ=n_circ,
            n_axial=n_axial,
            n_tip=n_tip,
            center=center,
            axis=axis,
            cap_tail=cap_tail,
        )

    def generate(self) -> MeshModel:
        radius = self.radius
        body_end = self.length - radius

        # Constant-radius cylinder, then a hemispherical tip.
        s_body = np.linspace(0.0, body_end, self.n_axial + 1)
        r_body = np.full_like(s_body, radius)
        phi = np.linspace(0.0, np.pi / 2.0, self.n_tip + 1)[1:]
        s_all = np.concatenate([s_body, body_end + radius * np.sin(phi)])
        r_all = np.concatenate([r_body, radius * np.cos(phi)])

        centers = self.center + np.outer(s_all, self.axis)
        return build_revolved_shell(
            centers,
            r_all,
            n_circ=self.n_circ,
            axis=self.axis,
            cap_start=self.cap_tail,
            cap_end=False,
            start_name="tail",
            end_name="nose",
        )
