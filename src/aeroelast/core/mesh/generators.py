"""
Mesh generators module.

This module contains classes for generating various types of structured meshes:
- SquareShapeMesh: 2D rectangular meshes
- BoxSurfaceMesh: 3D box surface meshes
- MultiFlapMesh: Multi-flap structures for FSI simulations
"""

from __future__ import annotations

import importlib
import math
from typing import TYPE_CHECKING, Any, List, Tuple

import numpy as np

from aeroelast.core.mesh.entities import (
    ELEMENT_NODES_MAP,
    ElementSet,
    ElementType,
    MeshElement,
    Node,
    NodeSet,
)


# Populated on first use by ``_import_gmsh``; kept as an explicit module-level
# name so that the helper methods which reference the bare ``gmsh`` are
# statically resolvable instead of relying on a dynamic injection into this
# module's dict.
gmsh = None


def _import_gmsh():
    """Lazy import of gmsh to avoid loading libGLU on headless systems."""
    global gmsh
    if gmsh is None:
        gmsh = importlib.import_module("gmsh")
    return gmsh


if TYPE_CHECKING:
    from aeroelast.core.mesh.model import MeshModel


class SquareShapeMesh:
    """
    Generates structured 2D meshes using Gmsh and returns a MeshModel instance.

    Attributes
    ----------
    width : float
        Domain width (x-direction)
    height : float
        Domain height (y-direction)
    nx : int
        Number of divisions in x-direction
    ny : int
        Number of divisions in y-direction
    quadratic : bool
        Use quadratic elements
    triangular : bool
        Use triangular elements
    distorted : bool
        Apply Ko2017 ratio-based mesh distortion (default: False)
    """

    def __init__(
        self,
        width: float,
        height: float,
        nx: int,
        ny: int,
        quadratic: bool = False,
        triangular: bool = False,
        distorted: bool = False,
    ):
        self.width = width
        self.height = height
        self.nx = nx
        self.ny = ny
        self.quadratic = quadratic
        self.triangular = triangular
        self.distorted = distorted

    def _apply_ko2017_distortion(self, mesh_model: "MeshModel") -> None:
        """Apply Ko2017 ratio-based mesh distortion in-place."""
        if not self.distorted:
            return

        def _ratio_positions(n: int) -> np.ndarray:
            seg = np.arange(1, n + 1, dtype=float)
            cum = np.concatenate([[0.0], np.cumsum(seg)])
            return cum / cum[-1]

        # Domain is [-width/2, width/2] x [0, height]
        x0, x1 = -self.width / 2, self.width / 2
        y0, y1 = 0.0, self.height

        x_uniform = np.linspace(x0, x1, self.nx + 1)
        y_uniform = np.linspace(y0, y1, self.ny + 1)
        x_dist = x0 + _ratio_positions(self.nx) * (x1 - x0)
        y_dist = y0 + _ratio_positions(self.ny) * (y1 - y0)

        for node in mesh_model.nodes:
            x, y, z = node.coords
            i_x = int(np.argmin(np.abs(x_uniform - x)))
            i_y = int(np.argmin(np.abs(y_uniform - y)))
            node.coords[0] = float(x_dist[i_x])
            node.coords[1] = float(y_dist[i_y])

    def generate(self) -> "MeshModel":
        """Generates and returns a MeshModel with the structured mesh"""
        # Import here to avoid circular imports
        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("rectangle")

            self._create_geometry()
            self._configure_mesh()

            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.optimize()

            mesh_model = self._create_mesh_model(MeshModel)
            self._apply_ko2017_distortion(mesh_model)
            return mesh_model
        finally:
            gmsh.finalize()

    def _create_geometry(self):
        """Create geometric entities with proper boundary handling"""
        x0 = -self.width / 2
        x1 = self.width / 2
        y0 = 0.0
        y1 = self.height

        self.p1 = gmsh.model.geo.addPoint(x0, y0, 0)
        self.p2 = gmsh.model.geo.addPoint(x1, y0, 0)
        self.p3 = gmsh.model.geo.addPoint(x1, y1, 0)
        self.p4 = gmsh.model.geo.addPoint(x0, y1, 0)

        self.bottom = gmsh.model.geo.addLine(self.p1, self.p2)
        self.right = gmsh.model.geo.addLine(self.p2, self.p3)
        self.top = gmsh.model.geo.addLine(self.p3, self.p4)
        self.left = gmsh.model.geo.addLine(self.p4, self.p1)

        loop = gmsh.model.geo.addCurveLoop([self.bottom, self.right, self.top, self.left])
        self.surface = gmsh.model.geo.addPlaneSurface([loop])

        self._add_physical_groups()

    def _add_physical_groups(self):
        """Add physical groups with proper corner handling"""
        gmsh.model.addPhysicalGroup(1, [self.top], name="top")
        gmsh.model.addPhysicalGroup(1, [self.bottom], name="bottom")
        gmsh.model.addPhysicalGroup(1, [self.left], name="left")
        gmsh.model.addPhysicalGroup(1, [self.right], name="right")

        gmsh.model.addPhysicalGroup(0, [self.p1], name="corner_p1")
        gmsh.model.addPhysicalGroup(0, [self.p2], name="corner_p2")
        gmsh.model.addPhysicalGroup(0, [self.p3], name="corner_p3")
        gmsh.model.addPhysicalGroup(0, [self.p4], name="corner_p4")

    def _configure_mesh(self):
        """Configure meshing parameters"""
        for curve in [self.bottom, self.top]:
            gmsh.model.geo.mesh.setTransfiniteCurve(curve, self.nx + 1)
        for curve in [self.left, self.right]:
            gmsh.model.geo.mesh.setTransfiniteCurve(curve, self.ny + 1)

        gmsh.model.geo.mesh.setTransfiniteSurface(
            self.surface, "Right", [self.p1, self.p2, self.p3, self.p4]
        )

        if self.triangular:
            gmsh.option.setNumber("Mesh.Algorithm", 6)
        else:
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
            gmsh.option.setNumber("Mesh.Algorithm", 8)

        if self.quadratic:
            gmsh.option.setNumber("Mesh.ElementOrder", 2)

    def _create_mesh_model(self, MeshModelClass) -> "MeshModel":
        """Converts Gmsh mesh to MeshModel instance"""
        mesh_model = MeshModelClass()
        self.tag_map = {}
        # Reset Node ID counter
        from aeroelast.core.mesh.entities import Node

        Node._id_counter = 0

        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(2)

        geometric_node_tags = set()
        for et, _conn in zip(elementTypes, nodeTags, strict=False):
            props = gmsh.model.mesh.getElementProperties(et)
            total_nodes = props[3]
            e_type = ELEMENT_NODES_MAP[total_nodes]
            if e_type in (ElementType.quad, ElementType.quad8, ElementType.quad9):
                num_corners = 4
            elif e_type in (ElementType.triangle, ElementType.triangle6):
                num_corners = 3
            else:
                raise ValueError(f"Unknown element type with {total_nodes} nodes")
            geometric_node_tags.update(
                nodeTags[0].reshape(-1, total_nodes)[:, :num_corners].flatten()
            )

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)

        # Sort by tag to fix ID mapping assumption
        p = np.argsort(node_tags)
        node_tags = node_tags[p]
        coords = coords[p]

        for tag, coord in zip(node_tags, coords, strict=False):
            if tag in geometric_node_tags:
                n = Node(coord, geometric_node=True)
            else:
                n = Node(coord, geometric_node=False)
            mesh_model.add_node(n)
            self.tag_map[tag] = n

        self._add_elements(mesh_model)
        self._create_node_sets(mesh_model)

        return mesh_model

    def _add_elements(self, mesh_model: "MeshModel"):
        """Add elements with corrected node indices"""
        elem_types = gmsh.model.mesh.getElementTypes()
        for elem_type in elem_types:
            elem_props = gmsh.model.mesh.getElementProperties(elem_type)
            if elem_props[1] == 2:
                _, elem_node_tags = gmsh.model.mesh.getElementsByType(elem_type)
                num_nodes = elem_props[3]
                connectivity = elem_node_tags.reshape(-1, num_nodes)

                for nodes in connectivity:
                    node_objs = [self.tag_map[nt] for nt in nodes]
                    e_type = ELEMENT_NODES_MAP.get(len(nodes), ElementType.quad)
                    mesh_model.add_element(MeshElement(nodes=node_objs, element_type=e_type))

    def _create_node_sets(self, mesh_model: "MeshModel"):
        """Creates node sets including boundary corners"""
        physical_groups = gmsh.model.getPhysicalGroups()

        # Store GMsh Tags in sets
        boundary_sets = {"top": set(), "bottom": set(), "left": set(), "right": set()}

        for dim, tag in physical_groups:
            name = gmsh.model.getPhysicalName(dim, tag)
            node_tags_list = []

            if dim == 0:
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
                for e in entities:
                    nt, _, _ = gmsh.model.mesh.getNodes(dim=0, tag=e)
                    node_tags_list.extend([int(n) for n in nt])

                if "corner" in name:
                    if "p1" in name:
                        boundary_sets["bottom"].update(node_tags_list)
                        boundary_sets["left"].update(node_tags_list)
                    elif "p2" in name:
                        boundary_sets["bottom"].update(node_tags_list)
                        boundary_sets["right"].update(node_tags_list)
                    elif "p3" in name:
                        boundary_sets["top"].update(node_tags_list)
                        boundary_sets["right"].update(node_tags_list)
                    elif "p4" in name:
                        boundary_sets["top"].update(node_tags_list)
                        boundary_sets["left"].update(node_tags_list)

            elif dim == 1:
                entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
                for e in entities:
                    nt, _, _ = gmsh.model.mesh.getNodes(dim=1, tag=e)
                    node_tags_list.extend([int(n) for n in nt])

                if name in boundary_sets:
                    boundary_sets[name].update(node_tags_list)

        boundary_nodes_objs = set()

        for name, tag_list in boundary_sets.items():
            if tag_list:
                node_objs = {self.tag_map[t] for t in tag_list if t in self.tag_map}
                mesh_model.add_node_set(NodeSet(name=name, nodes=node_objs))
                boundary_nodes_objs.update(node_objs)

        all_nodes = set(mesh_model.nodes)
        mesh_model.add_node_set(NodeSet(name="all", nodes=all_nodes))

        surface_nodes = all_nodes - boundary_nodes_objs
        mesh_model.add_node_set(NodeSet(name="surface", nodes=surface_nodes))

    @classmethod
    def create_rectangle(
        cls,
        width: float,
        height: float,
        nx: int,
        ny: int,
        quadratic: bool = False,
        triangular: bool = False,
    ) -> "MeshModel":
        """Helper method to create rectangular mesh"""
        return cls(width, height, int(nx), int(ny), quadratic, triangular).generate()

    @classmethod
    def create_unit_square(
        cls,
        nx: int,
        ny: int,
        quadratic: bool = False,
        triangular: bool = False,
    ) -> "MeshModel":
        """Helper method to create unit square mesh"""
        return cls(1.0, 1.0, nx, ny, quadratic, triangular).generate()


class BoxSurfaceMesh:
    """
    Generates structured 3D box surface meshes using Gmsh's classic geo API.

    Parameters
    ----------
    center : Tuple[float, float, float]
        Center coordinates of the box (x, y, z)
    dims : Tuple[float, float, float]
        Total box dimensions (dx, dy, dz)
    nx : int
        Number of divisions in x-direction
    ny : int
        Number of divisions in y-direction
    nz : int
        Number of divisions in z-direction
    quadratic : bool, optional
        Use quadratic elements (default: False)
    triangular : bool, optional
        Generate triangular mesh (default: False)
    """

    def __init__(
        self,
        center: Tuple[float, float, float],
        dims: Tuple[float, float, float],
        nx: int,
        ny: int,
        nz: int,
        quadratic: bool = False,
        triangular: bool = False,
    ):
        self.center = center
        self.dims = dims
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.quadratic = quadratic
        self.triangular = triangular

    def generate(self) -> "MeshModel":
        """Generates and returns the MeshModel instance"""
        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("box_surface")

            self._create_geometry()
            gmsh.model.geo.synchronize()

            self._configure_mesh()
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.optimize()

            return self._create_mesh_model(MeshModel)
        finally:
            gmsh.finalize()

    def _create_geometry(self):
        cx, cy, cz = self.center
        dx, dy, dz = self.dims

        x = (cx - dx / 2, cx + dx / 2)
        y = (cy - dy / 2, cy + dy / 2)
        z = (cz - dz / 2, cz + dz / 2)

        self.points = {
            "p1": gmsh.model.geo.addPoint(x[0], y[0], z[0]),
            "p2": gmsh.model.geo.addPoint(x[1], y[0], z[0]),
            "p3": gmsh.model.geo.addPoint(x[1], y[0], z[1]),
            "p4": gmsh.model.geo.addPoint(x[0], y[0], z[1]),
            "p5": gmsh.model.geo.addPoint(x[0], y[1], z[0]),
            "p6": gmsh.model.geo.addPoint(x[1], y[1], z[0]),
            "p7": gmsh.model.geo.addPoint(x[1], y[1], z[1]),
            "p8": gmsh.model.geo.addPoint(x[0], y[1], z[1]),
        }

        self._create_edges()
        self._create_faces()

        gmsh.model.geo.synchronize()

    def _create_edges(self):
        p = self.points
        self.edges = {
            "l1": gmsh.model.geo.addLine(p["p1"], p["p2"]),
            "l2": gmsh.model.geo.addLine(p["p2"], p["p3"]),
            "l3": gmsh.model.geo.addLine(p["p3"], p["p4"]),
            "l4": gmsh.model.geo.addLine(p["p4"], p["p1"]),
            "l5": gmsh.model.geo.addLine(p["p5"], p["p6"]),
            "l6": gmsh.model.geo.addLine(p["p6"], p["p7"]),
            "l7": gmsh.model.geo.addLine(p["p7"], p["p8"]),
            "l8": gmsh.model.geo.addLine(p["p8"], p["p5"]),
            "l9": gmsh.model.geo.addLine(p["p1"], p["p5"]),
            "l10": gmsh.model.geo.addLine(p["p2"], p["p6"]),
            "l11": gmsh.model.geo.addLine(p["p3"], p["p7"]),
            "l12": gmsh.model.geo.addLine(p["p4"], p["p8"]),
        }

    def _create_faces(self):
        """Creates box faces with proper orientation and physical groups"""
        edges = self.edges
        self.faces = {
            "bottom": self._create_face_loop([edges["l1"], edges["l2"], edges["l3"], edges["l4"]]),
            "top": self._create_face_loop([edges["l5"], edges["l6"], edges["l7"], edges["l8"]]),
            "front": self._create_face_loop([
                edges["l3"],
                edges["l12"],
                -edges["l7"],
                -edges["l11"],
            ]),
            "back": self._create_face_loop([edges["l1"], edges["l10"], -edges["l5"], -edges["l9"]]),
            "left": self._create_face_loop([-edges["l4"], edges["l12"], edges["l8"], -edges["l9"]]),
            "right": self._create_face_loop([
                edges["l2"],
                edges["l11"],
                -edges["l6"],
                -edges["l10"],
            ]),
        }

        for name, face_tag in self.faces.items():
            gmsh.model.addPhysicalGroup(2, [face_tag], name=name)

    def _create_face_loop(self, curves):
        """Helper to create face from oriented curves"""
        loop = gmsh.model.geo.addCurveLoop(curves)
        return gmsh.model.geo.addPlaneSurface([loop])

    def _configure_mesh(self):
        """Configures mesh parameters and transfinite settings"""
        self._set_transfinite_curves()
        self._set_transfinite_surfaces()

        if self.triangular:
            gmsh.option.setNumber("Mesh.Algorithm", 6)
        else:
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.option.setNumber("Mesh.RecombineAll", 1)

        if self.quadratic:
            gmsh.option.setNumber("Mesh.ElementOrder", 2)

    def _set_transfinite_curves(self):
        """Configure transfinite curves"""
        for e in ["l1", "l3", "l5", "l7"]:
            gmsh.model.mesh.setTransfiniteCurve(self.edges[e], self.nx + 1)

        for e in ["l2", "l4", "l6", "l8"]:
            gmsh.model.mesh.setTransfiniteCurve(self.edges[e], self.nz + 1)

        for e in ["l9", "l10", "l11", "l12"]:
            gmsh.model.mesh.setTransfiniteCurve(self.edges[e], self.ny + 1)

    def _set_transfinite_surfaces(self):
        """Configures transfinite surfaces"""
        for face in self.faces.values():
            gmsh.model.mesh.setTransfiniteSurface(face)
            gmsh.model.mesh.setRecombine(2, face)

    def _create_mesh_model(self, MeshModelClass) -> "MeshModel":
        """Converts Gmsh mesh to MeshModel instance"""
        mesh_model = MeshModelClass()

        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(2)

        geometric_node_tags = set()
        for et, _conn in zip(elementTypes, nodeTags, strict=False):
            props = gmsh.model.mesh.getElementProperties(et)
            total_nodes = props[3]
            e_type = ELEMENT_NODES_MAP[total_nodes]
            if e_type in (ElementType.quad, ElementType.quad8, ElementType.quad9):
                num_corners = 4
            elif e_type in (ElementType.triangle, ElementType.triangle6):
                num_corners = 3
            else:
                raise ValueError(f"Unknown element type with {total_nodes} nodes")
            geometric_node_tags.update(
                nodeTags[0].reshape(-1, total_nodes)[:, :num_corners].flatten()
            )

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)

        for tag, coord in zip(node_tags, coords, strict=False):
            if tag in geometric_node_tags:
                mesh_model.add_node(Node(coord, geometric_node=True))
            else:
                mesh_model.add_node(Node(coord, geometric_node=False))

        self._add_elements(mesh_model)
        self._create_node_sets(mesh_model)

        return mesh_model

    def _add_elements(self, mesh_model: "MeshModel"):
        """Add elements with corrected node indices"""
        elem_types = gmsh.model.mesh.getElementTypes()
        for elem_type in elem_types:
            elem_props = gmsh.model.mesh.getElementProperties(elem_type)
            if elem_props[1] == 2:
                _, elem_node_tags = gmsh.model.mesh.getElementsByType(elem_type)
                num_nodes = elem_props[3]
                connectivity = elem_node_tags.reshape(-1, num_nodes)

                for nodes in connectivity:
                    node_objs = [mesh_model.get_node_by_id(int(nt - 1)) for nt in nodes]
                    e_type = ELEMENT_NODES_MAP.get(len(nodes), ElementType.quad)
                    mesh_model.add_element(MeshElement(nodes=node_objs, element_type=e_type))

    def _create_node_sets(self, mesh_model: "MeshModel"):
        """Creates node sets for all boundary faces"""
        physical_groups = gmsh.model.getPhysicalGroups()

        face_sets = {
            "top": set(),
            "bottom": set(),
            "front": set(),
            "back": set(),
            "left": set(),
            "right": set(),
        }

        for dim, tag in physical_groups:
            if dim != 2:
                continue

            name = gmsh.model.getPhysicalName(dim, tag)
            if name not in face_sets:
                continue

            node_tags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
            face_sets[name].update((tag - 1 for tag in node_tags))

        for name, node_ids in face_sets.items():
            if node_ids:
                node_objs = {mesh_model.get_node_by_id(nid) for nid in node_ids}
                mesh_model.add_node_set(NodeSet(name=name, nodes=node_objs))

        all_nodes = set(mesh_model.nodes)
        mesh_model.add_node_set(NodeSet(name="all", nodes=all_nodes))

    @classmethod
    def create_box(
        cls,
        center: Tuple[float, float, float],
        dims: Tuple[float, float, float],
        nx: int,
        ny: int,
        nz: int,
        quadratic: bool = False,
        triangular: bool = False,
    ) -> "MeshModel":
        """Creates a box with specified dimensions"""
        return cls(center, dims, nx, ny, nz, quadratic, triangular).generate()

    @classmethod
    def create_unit_box(
        cls,
        divisions: int = 10,
        quadratic: bool = False,
        triangular: bool = False,
    ) -> "MeshModel":
        """Creates a unit cube (1x1x1) centered at origin"""
        return cls(
            (0, 0, 0), (1, 1, 1), divisions, divisions, divisions, quadratic, triangular
        ).generate()


class MultiFlapMesh:
    """
    Generates a structured 2D mesh for multiple flaps connected to a common base.

    The geometry consists of:
    - A horizontal base (rectangle) at y=0 with fixed boundary condition
    - Multiple vertical flaps extending upward from the base

    Parameters
    ----------
    n_flaps : int
        Number of flaps
    flap_width : float
        Width of each flap (X direction)
    flap_height : float
        Height of each flap (Y direction)
    x_spacing : float
        Spacing between consecutive flaps (edge to edge)
    base_height : float
        Height of the base strip (Y direction)
    nx_flap : int
        Number of divisions in X direction for each flap
    ny_flap : int
        Number of divisions in Y direction for each flap
    nx_base_segment : int
        Number of divisions in X for base segments between flaps
    ny_base : int
        Number of divisions in Y for the base
    quadratic : bool
        Use quadratic elements (default: False)
    """

    def __init__(
        self,
        n_flaps: int,
        flap_width: float,
        flap_height: float,
        x_spacing: float,
        base_height: float = 0.05,
        nx_flap: int = 4,
        ny_flap: int = 20,
        nx_base_segment: int = 10,
        ny_base: int = 2,
        quadratic: bool = False,
    ):
        self.n_flaps = n_flaps
        self.flap_width = flap_width
        self.flap_height = flap_height
        self.x_spacing = x_spacing
        self.base_height = base_height
        self.nx_flap = nx_flap
        self.ny_flap = ny_flap
        self.nx_base_segment = nx_base_segment
        self.ny_base = ny_base
        self.quadratic = quadratic

        self.points = {}
        self.lines = {}
        self.surfaces = {}
        self.flap_surfaces = []
        self.base_surfaces = []

    def _calculate_positions(self) -> List[dict]:
        """Calculate X positions of all flaps centered around x=0"""
        total_width = self.n_flaps * self.flap_width + (self.n_flaps - 1) * self.x_spacing
        x_start = -total_width / 2

        positions = []
        for i in range(self.n_flaps):
            x_left = x_start + i * (self.flap_width + self.x_spacing)
            positions.append({
                "index": i + 1,
                "x_left": x_left,
                "x_right": x_left + self.flap_width,
            })
        return positions

    def generate(self) -> "MeshModel":
        """Generates and returns the MeshModel instance"""
        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("multiflap")

            self._create_geometry()
            gmsh.model.geo.synchronize()

            self._configure_mesh()
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.optimize()

            return self._create_mesh_model(MeshModel)
        finally:
            gmsh.finalize()

    def _create_geometry(self):
        """Create the multi-flap geometry with base"""
        flap_positions = self._calculate_positions()

        y_base_bottom = -self.base_height
        y_base_top = 0.0
        y_flap_top = self.flap_height

        x_left_total = flap_positions[0]["x_left"]
        x_right_total = flap_positions[-1]["x_right"]

        point_id = 1

        base_bottom_points = []
        x_coords = [x_left_total]
        for fp in flap_positions:
            x_coords.extend([fp["x_left"], fp["x_right"]])
        x_coords.append(x_right_total)
        x_coords = sorted(set(x_coords))

        for x in x_coords:
            self.points[f"bb_{point_id}"] = gmsh.model.geo.addPoint(x, y_base_bottom, 0)
            base_bottom_points.append((x, self.points[f"bb_{point_id}"]))
            point_id += 1

        base_top_points = []
        for x in x_coords:
            self.points[f"bt_{point_id}"] = gmsh.model.geo.addPoint(x, y_base_top, 0)
            base_top_points.append((x, self.points[f"bt_{point_id}"]))
            point_id += 1

        flap_top_points = []
        for fp in flap_positions:
            p_left = gmsh.model.geo.addPoint(fp["x_left"], y_flap_top, 0)
            p_right = gmsh.model.geo.addPoint(fp["x_right"], y_flap_top, 0)
            self.points[f"ft_{fp['index']}_left"] = p_left
            self.points[f"ft_{fp['index']}_right"] = p_right
            flap_top_points.append({
                "index": fp["index"],
                "x_left": fp["x_left"],
                "x_right": fp["x_right"],
                "p_left": p_left,
                "p_right": p_right,
            })

        self._create_base_structure(base_bottom_points, base_top_points, flap_positions)
        self._create_flaps(base_top_points, flap_top_points, flap_positions)
        self._add_physical_groups(flap_positions)

    def _create_base_structure(self, base_bottom_points, base_top_points, flap_positions):
        """Create the base horizontal strip"""
        bottom_lines = []
        top_lines = []

        for i in range(len(base_bottom_points) - 1):
            l_bottom = gmsh.model.geo.addLine(
                base_bottom_points[i][1], base_bottom_points[i + 1][1]
            )
            l_top = gmsh.model.geo.addLine(base_top_points[i][1], base_top_points[i + 1][1])
            bottom_lines.append(l_bottom)
            top_lines.append(l_top)

        vertical_lines = []
        for i in range(len(base_bottom_points)):
            l_vert = gmsh.model.geo.addLine(base_bottom_points[i][1], base_top_points[i][1])
            vertical_lines.append(l_vert)

        for i in range(len(bottom_lines)):
            loop = gmsh.model.geo.addCurveLoop([
                bottom_lines[i],
                vertical_lines[i + 1],
                -top_lines[i],
                -vertical_lines[i],
            ])
            surf = gmsh.model.geo.addPlaneSurface([loop])
            self.base_surfaces.append(surf)

        self.lines["base_bottom"] = bottom_lines
        self.lines["base_top"] = top_lines
        self.lines["base_vertical"] = vertical_lines

    def _create_flaps(self, base_top_points, flap_top_points, flap_positions):
        """Create the vertical flaps extending from the base"""
        base_top_lookup = {round(x, 10): p for x, p in base_top_points}

        for ftp in flap_top_points:
            idx = ftp["index"]
            x_left = round(ftp["x_left"], 10)
            x_right = round(ftp["x_right"], 10)

            p_base_left = base_top_lookup[x_left]
            p_base_right = base_top_lookup[x_right]

            p_top_left = ftp["p_left"]
            p_top_right = ftp["p_right"]

            l_left = gmsh.model.geo.addLine(p_base_left, p_top_left)
            l_top = gmsh.model.geo.addLine(p_top_left, p_top_right)
            l_right = gmsh.model.geo.addLine(p_top_right, p_base_right)
            l_bottom = gmsh.model.geo.addLine(p_base_right, p_base_left)

            loop = gmsh.model.geo.addCurveLoop([l_left, l_top, l_right, l_bottom])
            surf = gmsh.model.geo.addPlaneSurface([loop])

            self.flap_surfaces.append({
                "index": idx,
                "surface": surf,
                "l_left": l_left,
                "l_top": l_top,
                "l_right": l_right,
                "l_bottom": l_bottom,
            })

    def _add_physical_groups(self, flap_positions):
        """Add physical groups for boundary conditions"""
        all_bottom_lines = self.lines["base_bottom"]
        gmsh.model.addPhysicalGroup(1, all_bottom_lines, name="base_bottom")

        gmsh.model.addPhysicalGroup(1, [self.lines["base_vertical"][0]], name="base_left")
        gmsh.model.addPhysicalGroup(1, [self.lines["base_vertical"][-1]], name="base_right")

        for flap_data in self.flap_surfaces:
            idx = flap_data["index"]
            gmsh.model.addPhysicalGroup(1, [flap_data["l_left"]], name=f"flap{idx}_left")
            gmsh.model.addPhysicalGroup(1, [flap_data["l_top"]], name=f"flap{idx}_top")
            gmsh.model.addPhysicalGroup(1, [flap_data["l_right"]], name=f"flap{idx}_right")
            gmsh.model.addPhysicalGroup(2, [flap_data["surface"]], name=f"flap{idx}")

        all_flap_left = [f["l_left"] for f in self.flap_surfaces]
        all_flap_top = [f["l_top"] for f in self.flap_surfaces]
        all_flap_right = [f["l_right"] for f in self.flap_surfaces]
        gmsh.model.addPhysicalGroup(1, all_flap_left, name="flaps_left")
        gmsh.model.addPhysicalGroup(1, all_flap_top, name="flaps_top")
        gmsh.model.addPhysicalGroup(1, all_flap_right, name="flaps_right")

        all_surfaces = self.base_surfaces + [f["surface"] for f in self.flap_surfaces]
        gmsh.model.addPhysicalGroup(2, all_surfaces, name="all_surfaces")

    def _configure_mesh(self):
        """Configure mesh parameters"""
        for line in self.lines["base_bottom"]:
            gmsh.model.mesh.setTransfiniteCurve(line, self.nx_base_segment + 1)

        for line in self.lines["base_top"]:
            gmsh.model.mesh.setTransfiniteCurve(line, self.nx_base_segment + 1)

        for line in self.lines["base_vertical"]:
            gmsh.model.mesh.setTransfiniteCurve(line, self.ny_base + 1)

        for surf in self.base_surfaces:
            gmsh.model.mesh.setTransfiniteSurface(surf)
            gmsh.model.mesh.setRecombine(2, surf)

        for flap_data in self.flap_surfaces:
            gmsh.model.mesh.setTransfiniteCurve(flap_data["l_left"], self.ny_flap + 1)
            gmsh.model.mesh.setTransfiniteCurve(flap_data["l_right"], self.ny_flap + 1)
            gmsh.model.mesh.setTransfiniteCurve(flap_data["l_top"], self.nx_flap + 1)
            gmsh.model.mesh.setTransfiniteCurve(flap_data["l_bottom"], self.nx_flap + 1)

            gmsh.model.mesh.setTransfiniteSurface(flap_data["surface"])
            gmsh.model.mesh.setRecombine(2, flap_data["surface"])

        gmsh.option.setNumber("Mesh.Algorithm", 8)
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

        if self.quadratic:
            gmsh.option.setNumber("Mesh.ElementOrder", 2)

    def _create_mesh_model(self, MeshModelClass) -> "MeshModel":
        """Converts Gmsh mesh to MeshModel instance"""
        mesh_model = MeshModelClass()

        elementTypes, elementTags, nodeTags = gmsh.model.mesh.getElements(2)

        geometric_node_tags = set()
        for et, _conn in zip(elementTypes, nodeTags, strict=False):
            props = gmsh.model.mesh.getElementProperties(et)
            total_nodes = props[3]
            e_type = ELEMENT_NODES_MAP[total_nodes]
            if e_type in (ElementType.quad, ElementType.quad8, ElementType.quad9):
                num_corners = 4
            elif e_type in (ElementType.triangle, ElementType.triangle6):
                num_corners = 3
            else:
                raise ValueError(f"Unsupported element type: {e_type}")
            geometric_node_tags.update(
                nodeTags[0].reshape(-1, total_nodes)[:, :num_corners].flatten()
            )

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)

        self._gmsh_tag_to_id = {}
        for i, (tag, coord) in enumerate(zip(node_tags, coords, strict=False)):
            is_geometric = tag in geometric_node_tags
            mesh_model.add_node(Node(coord, geometric_node=is_geometric))
            self._gmsh_tag_to_id[int(tag)] = i

        self._add_elements(mesh_model)
        self._create_node_sets(mesh_model)

        return mesh_model

    def _add_elements(self, mesh_model: "MeshModel"):
        """Add elements to the mesh model"""
        elem_types = gmsh.model.mesh.getElementTypes()
        for elem_type in elem_types:
            elem_props = gmsh.model.mesh.getElementProperties(elem_type)
            if elem_props[1] == 2:
                _, elem_node_tags = gmsh.model.mesh.getElementsByType(elem_type)
                num_nodes = elem_props[3]
                connectivity = elem_node_tags.reshape(-1, num_nodes)

                for nodes in connectivity:
                    node_objs = [
                        mesh_model.get_node_by_id(self._gmsh_tag_to_id[int(nt)]) for nt in nodes
                    ]

                    if len(node_objs) >= 4:
                        coords = [n.coords[:2] for n in node_objs[:4]]
                        area = 0.0
                        for i in range(4):
                            j = (i + 1) % 4
                            area += coords[i][0] * coords[j][1]
                            area -= coords[j][0] * coords[i][1]
                        if area < 0:
                            if len(node_objs) == 4:
                                node_objs = [node_objs[0], node_objs[3], node_objs[2], node_objs[1]]
                            elif len(node_objs) == 8:
                                node_objs = [
                                    node_objs[0],
                                    node_objs[3],
                                    node_objs[2],
                                    node_objs[1],
                                    node_objs[7],
                                    node_objs[6],
                                    node_objs[5],
                                    node_objs[4],
                                ]
                            elif len(node_objs) == 9:
                                node_objs = [
                                    node_objs[0],
                                    node_objs[3],
                                    node_objs[2],
                                    node_objs[1],
                                    node_objs[7],
                                    node_objs[6],
                                    node_objs[5],
                                    node_objs[4],
                                    node_objs[8],
                                ]

                    e_type = ELEMENT_NODES_MAP.get(len(nodes), ElementType.quad)
                    mesh_model.add_element(MeshElement(nodes=node_objs, element_type=e_type))

    def _create_node_sets(self, mesh_model: "MeshModel"):
        """Create node sets for boundaries"""
        physical_groups = gmsh.model.getPhysicalGroups()

        for dim, tag in physical_groups:
            name = gmsh.model.getPhysicalName(dim, tag)
            if dim == 1 and name:
                node_tags = gmsh.model.mesh.getNodesForPhysicalGroup(dim, tag)[0]
                node_ids = [self._gmsh_tag_to_id[int(nt)] for nt in node_tags]
                node_objs = {mesh_model.get_node_by_id(nid) for nid in node_ids}
                mesh_model.add_node_set(NodeSet(name=name, nodes=node_objs))

        all_nodes = set(mesh_model.nodes)
        mesh_model.add_node_set(NodeSet(name="all", nodes=all_nodes))

        if "base_bottom" in mesh_model.node_sets:
            bottom_nodes = mesh_model.node_sets["base_bottom"].nodes
            mesh_model.add_node_set(NodeSet(name="bottom", nodes=set(bottom_nodes.values())))

    @classmethod
    def create_from_config(
        cls,
        n_flaps: int,
        flap_width: float,
        flap_height: float,
        x_spacing: float,
        base_height: float = 0.05,
        nx_flap: int = 4,
        ny_flap: int = 20,
        nx_base_segment: int = 10,
        ny_base: int = 2,
        quadratic: bool = False,
    ) -> "MeshModel":
        """Create multi-flap mesh from parameters"""
        return cls(
            n_flaps=n_flaps,
            flap_width=flap_width,
            flap_height=flap_height,
            x_spacing=x_spacing,
            base_height=base_height,
            nx_flap=nx_flap,
            ny_flap=ny_flap,
            nx_base_segment=nx_base_segment,
            ny_base=ny_base,
            quadratic=quadratic,
        ).generate()


class BladeMesh:
    """
    Generates shell meshes for wind turbine blades from YAML definitions.

    This generator uses the numad library to create parametric blade geometries
    and mesh them with shell elements. The resulting mesh includes element sets
    for different blade regions and sections.

    Attributes
    ----------
    yaml_file : str
        Path to the blade YAML definition file
    element_size : float
        Target element size for meshing
    n_samples : int
        Number of samples for airfoil discretization
    include_webs : bool
        When *False* the shear webs are omitted and only the blade outer shell
        (outer-mold-line surface) is meshed.  Use this for CFD wall surfaces.

    Examples
    --------
    >>> blade_mesh = BladeMesh("blade_definition.yaml", element_size=0.05)
    >>> mesh = blade_mesh.generate(renumber="rcm")
    >>> mesh.write_mesh("blade.vtk")

    >>> cfd = BladeMesh("blade_definition.yaml", include_webs=False)
    >>> cfd.generate().write_mesh("blade_surface.stl")
    """

    def __init__(
        self,
        yaml_file: str = None,
        element_size: float = 0.1,
        n_samples: int = 300,
        excel_file: str = None,
        airfoil_dir: str = None,
        refine_tip: bool = True,
        span_grading: str = "chord",
        airfoil_spacing: str = "cosine",
        include_webs: bool = True,
    ):
        self.yaml_file = yaml_file
        self.excel_file = excel_file
        self.airfoil_dir = airfoil_dir
        self.element_size = element_size
        self.n_samples = n_samples
        self.refine_tip = refine_tip
        self.span_grading = span_grading
        self.airfoil_spacing = airfoil_spacing
        self.include_webs = include_webs
        # The numad objects come from the untyped ported mesher (``nuMAD``): their
        # attributes are optional in the stubs but populated here, so ``Any`` is the honest
        # annotation and keeps the Optional out of every downstream subscript.
        self._numad_blade: Any = None
        self._numad_mesh: Any = None

    def generate(self, renumber: str | None = None, verbose: bool = True) -> "MeshModel":
        """
        Generate the blade mesh from the YAML definition.

        Parameters
        ----------
        renumber : str, optional
            Renumbering algorithm to apply after mesh generation.
            - None (default): No renumbering
            - "simple": Direct index assignment
            - "rcm": Reverse Cuthill-McKee for bandwidth reduction
        verbose : bool, optional
            Whether to print progress messages (default: True)

        Returns
        -------
        MeshModel
            The generated mesh model with nodes, elements, and sets
        """
        # Import here to avoid circular imports and optional dependency
        from aeroelast.core.mesh.model import MeshModel

        try:
            from aeroelast.models.blade.numad import Blade as numadBlade
            from aeroelast.models.blade.numad.mesh_gen import get_shell_mesh
        except ImportError as e:
            raise ImportError(
                "BladeMesh requires the numad module. "
                "Make sure aeroelast.models.blade.numad is available."
            ) from e

        mesh_model = MeshModel()

        # Initialize numad blade
        if verbose:
            print("      Loading blade definition...")
        self._numad_blade = numadBlade()
        if self.excel_file:
            self._numad_blade.read_excel(self.excel_file, airfoil_dir=self.airfoil_dir)
        elif self.yaml_file:
            self._numad_blade.read_yaml(self.yaml_file)
        else:
            raise ValueError("BladeMesh requires either yaml_file or excel_file")

        # Resample airfoils
        for stat in self._numad_blade.definition.stations:
            stat.airfoil.resample(n_samples=self.n_samples, spacing=self.airfoil_spacing)

        # Update geometry
        if verbose:
            print("      Updating blade geometry...")
        self._numad_blade.update_blade()

        # Refine cross-sections where chord changes rapidly.  The default
        # YAML grid is uniform in span, so the last section can jump from
        # chord ≈ 1.8 m to 0.5 m in a single 2.4 m span gap.  We subdivide
        # sections whose chord ratio exceeds a threshold so that element
        # aspect ratios stay reasonable.
        if self.refine_tip:
            self._refine_high_gradient_sections(verbose=verbose)

        # Expand trailing edge
        n_stations = self._numad_blade.geometry.coordinates.shape[2]
        min_TE_lengths = 0.001 * np.ones(n_stations)
        self._numad_blade.expand_blade_geometry_te(min_TE_lengths)

        # Generate shell mesh using numad
        if verbose:
            print("      Generating shell mesh...")
        self._numad_mesh = get_shell_mesh(
            self._numad_blade,
            self.element_size,
            spanGrading=self.span_grading,
            include_webs=self.include_webs,
        )

        num_raw_nodes = len(self._numad_mesh["nodes"])
        num_elements = len(self._numad_mesh["elements"])
        if verbose:
            print(f"      Raw mesh: {num_raw_nodes} nodes, {num_elements} elements")

        # Deduplicate nodes
        self._deduplicate_and_create_mesh(mesh_model, verbose)

        # Renumber mesh if requested
        if renumber is not None:
            if verbose:
                print(f"      Renumbering mesh using {renumber} algorithm...")
            mesh_model.renumber_mesh(algorithm=renumber)

        if verbose:
            print(
                f"      Blade mesh generated: {mesh_model.node_count} nodes, {mesh_model.elements_count} elements"
            )

        return mesh_model

    def _refine_high_gradient_sections(
        self,
        chord_ratio_threshold=0.5,
        verbose=True,
    ):
        """Add interpolated stations where the chord ratio between adjacent
        sections is below *chord_ratio_threshold* (i.e. chord shrinks by
        more than 50 %).  Stations are subdivided until the ratio of the
        smaller chord to the larger chord in every pair exceeds the
        threshold, up to a maximum of 3 subdivisions per original gap.

        The refinement is unconditional, matching the reference mesh.  An
        earlier local workaround skipped a refinement when
        ``gap / 2 < element_size * 4``: that dropped 26 nodes at the tip and
        left a sliver element CalculiX cannot integrate
        (``*ERROR in e_c3d: nonpositive jacobian``, element 2790).  The
        workaround targeted the aspect ratio of the previous element, while
        the reviewed element and upstream's own mesh test both assume the
        refined mesh (9277 nodes at element_size = 0.5 m).
        """
        blade = self._numad_blade
        ichord = blade.geometry.ichord
        ispan = blade.ispan.copy()
        added = 0

        for _ in range(3):  # iterate because adding a station shifts indices
            ichord = blade.geometry.ichord
            ispan = blade.ispan.copy()
            n = len(ichord)
            inserted_any = False
            for i in range(n - 1):
                c_lo = min(ichord[i], ichord[i + 1])
                c_hi = max(ichord[i], ichord[i + 1])
                if c_hi > 0 and c_lo / c_hi < chord_ratio_threshold:
                    mid = 0.5 * (ispan[i] + ispan[i + 1])
                    blade.add_interpolated_station(mid)
                    added += 1
                    inserted_any = True
                    break  # restart scan — indices shifted
            if not inserted_any:
                break

        if verbose and added > 0:
            print(f"      Added {added} interpolated section(s) for tip refinement")

    def _deduplicate_and_create_mesh(self, mesh_model: "MeshModel", verbose: bool = True):
        """
        Deduplicate nodes and create mesh elements.

        Numad may create duplicate nodes at region boundaries, so we need to
        deduplicate them using spatial tolerance matching.
        """
        from scipy.spatial import cKDTree

        from aeroelast.core.mesh.entities import ElementSet

        node_coords = self._numad_mesh["nodes"]
        num_raw_nodes = len(node_coords)
        num_elements = len(self._numad_mesh["elements"])

        if verbose:
            print("      Deduplicating nodes...")

        tolerance = 1e-6
        numad_to_unique = {}  # Maps numad node index to unique node index
        unique_nodes = []

        # Build KDTree for efficient nearest neighbor search
        all_coords = np.array(node_coords)
        tree = cKDTree(all_coords)

        # Find duplicate groups using the tree
        processed = set()
        for idx in range(num_raw_nodes):
            if idx in processed:
                continue

            neighbors = tree.query_ball_point(all_coords[idx], r=tolerance)
            unique_idx = len(unique_nodes)
            unique_nodes.append(all_coords[idx])

            for neighbor_idx in neighbors:
                numad_to_unique[neighbor_idx] = unique_idx
                processed.add(neighbor_idx)

        num_unique_nodes = len(unique_nodes)
        if verbose:
            reduction_pct = 100 * (num_raw_nodes - num_unique_nodes) // max(1, num_raw_nodes)
            print(
                f"      Deduplicated: {num_raw_nodes} -> {num_unique_nodes} nodes ({reduction_pct}% reduction)"
            )

        # Create nodes
        node_index_to_node = {}
        for unique_idx, unique_coord in enumerate(unique_nodes):
            node = Node(unique_coord)
            mesh_model.add_node(node)
            node_index_to_node[unique_idx] = node

        # Create elements (skip degenerate elements with repeated nodes)
        if verbose:
            print(f"      Creating {num_elements} elements...")
        skipped = 0
        raw_to_mesh_id = {}  # maps raw element index → mesh element id
        for elem_idx, node_ids in enumerate(self._numad_mesh["elements"]):
            if node_ids[3] == -1:
                element_type = ElementType.triangle
                node_ids = node_ids[:3]
            else:
                element_type = ElementType.quad
            node_objs = [node_index_to_node[numad_to_unique[n_id]] for n_id in node_ids]
            # Skip elements where node deduplication collapsed nodes together
            unique_node_ids = {id(n) for n in node_objs}
            if len(unique_node_ids) < 3:
                skipped += 1
                continue
            elem = MeshElement(nodes=node_objs, element_type=element_type)
            mesh_model.add_element(elem)
            raw_to_mesh_id[elem_idx] = elem.id
        if skipped and verbose:
            print(f"      Skipped {skipped} degenerate element(s)")

        # Create element sets
        if verbose:
            print("      Creating element and node sets...")
        for element_set in self._numad_mesh["sets"]["element"]:
            name = element_set["name"]
            elements = set()
            for raw_idx in element_set["labels"]:
                if raw_idx in raw_to_mesh_id:
                    elements.add(mesh_model.get_element_by_id(raw_to_mesh_id[raw_idx]))
            if elements:
                mesh_model.add_element_set(ElementSet(name=name, elements=elements))

        # Create node sets
        for node_set in self._numad_mesh["sets"]["node"]:
            name = node_set["name"]
            nodes = {node_index_to_node[numad_to_unique[i]] for i in node_set["labels"]}
            mesh_model.add_node_set(NodeSet(name=name, nodes=nodes))

        # Create node sets from element sets with "all" in name
        for element_set in [
            eset for eset in mesh_model.element_sets.values() if "all" in eset.name.lower()
        ]:
            name = element_set.name.replace("Els", "Nods")
            nodes = {node for element in element_set.elements for node in element.nodes}
            mesh_model.add_node_set(NodeSet(name=name, nodes=nodes))

        # Assign per-element thickness from section layup data
        if self._numad_mesh.get("sections"):
            section_thickness = {}
            for section in self._numad_mesh["sections"]:
                set_name = section["elementSet"]
                total = sum(layer[1] for layer in section["layup"])
                section_thickness[set_name] = total
            for set_name, thickness in section_thickness.items():
                if set_name in mesh_model.element_sets:
                    for elem in mesh_model.element_sets[set_name].elements:
                        elem.thickness = thickness

    @property
    def numad_mesh_data(self):
        """
        Access the raw numad mesh data dictionary.

        This includes materials, sections, and other metadata from the blade definition.
        Only available after calling generate().

        Returns
        -------
        dict
            Raw mesh data from numad including nodes, elements, sets, materials, sections
        """
        if self._numad_mesh is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate() first.")
        return self._numad_mesh

    @property
    def numad_blade(self):
        """
        Access the numad blade object.

        Only available after calling generate().

        Returns
        -------
        numadBlade
            The numad Blade object with geometry and definition data
        """
        if self._numad_blade is None:
            raise RuntimeError("Blade has not been loaded yet. Call generate() first.")
        return self._numad_blade

    @classmethod
    def create_from_yaml(
        cls,
        yaml_file: str,
        element_size: float = 0.1,
        n_samples: int = 300,
        renumber: str = None,
        include_webs: bool = True,
    ) -> "MeshModel":
        """
        Convenience method to create blade mesh directly from YAML file.

        Parameters
        ----------
        yaml_file : str
            Path to the blade YAML definition file
        element_size : float, optional
            Target element size (default: 0.1)
        n_samples : int, optional
            Number of samples for airfoil discretization (default: 300)
        renumber : str, optional
            Renumbering algorithm ("simple", "rcm", or None)
        include_webs : bool, optional
            Set *False* to mesh only the outer shell (CFD surface)

        Returns
        -------
        MeshModel
            The generated mesh model
        """
        return cls(
            yaml_file=yaml_file,
            element_size=element_size,
            n_samples=n_samples,
            include_webs=include_webs,
        ).generate(renumber=renumber)


class RotorMesh:
    """
    Generates shell meshes for wind turbine rotors from blade YAML definitions.

    This generator creates a multi-blade rotor by generating a base blade mesh
    and then copying, translating, and rotating it for each blade position.

    The coordinate system follows wind turbine conventions:
    - Y-axis: Rotor rotation axis (blade span direction)
    - Z-axis: Radial direction (blades are offset in this direction)
    - X-axis: Tangential direction

    Attributes
    ----------
    yaml_file : str
        Path to the blade YAML definition file
    n_blades : int
        Number of blades in the rotor
    hub_radius : float or None
        Radial distance from rotation axis to blade root
    element_size : float
        Target element size for meshing
    n_samples : int
        Number of samples for airfoil discretization

    Examples
    --------
    >>> rotor = RotorMesh(yaml_file="blade.yaml", n_blades=3, hub_radius=1.5)
    >>> mesh = rotor.generate(renumber="rcm")
    >>> mesh.write_mesh("rotor.vtk")
    """

    BLADE_PART_SET_PREFIX = "rotor_blade_"

    def __init__(
        self,
        yaml_file: str | None = None,
        n_blades: int = 3,
        hub_radius: float | None = None,
        hub_diameter: float | None = None,
        rotor_diameter: float | None = None,
        element_size: float = 0.1,
        n_samples: int = 300,
        excel_file: str | None = None,
        airfoil_dir: str | None = None,
        airfoil_spacing: str = "cosine",
        include_webs: bool = True,
    ):
        self.yaml_file = yaml_file
        self.excel_file = excel_file
        self.airfoil_dir = airfoil_dir
        self.n_blades = n_blades
        self.hub_radius = hub_radius
        self.hub_diameter = hub_diameter
        self.rotor_diameter = rotor_diameter
        self.element_size = element_size
        self.n_samples = n_samples
        self.airfoil_spacing = airfoil_spacing
        self.include_webs = include_webs
        self._blade_generator: BladeMesh = None

    def generate(self, renumber: str | None = None, verbose: bool = True) -> "MeshModel":
        """
        Generate the rotor mesh by creating and transforming blade meshes.

        Parameters
        ----------
        renumber : str, optional
            Renumbering algorithm to apply after mesh generation.
            - None (default): No renumbering
            - "simple": Direct index assignment
            - "rcm": Reverse Cuthill-McKee for bandwidth reduction
        verbose : bool, optional
            Whether to print progress messages (default: True)

        Returns
        -------
        MeshModel
            The generated rotor mesh model
        """
        from aeroelast.core.mesh.model import MeshModel

        rotor_mesh = MeshModel()

        # Generate base blade mesh
        if verbose:
            print("\nGenerating base blade mesh...")
        self._blade_generator = BladeMesh(
            yaml_file=self.yaml_file,
            excel_file=self.excel_file,
            airfoil_dir=self.airfoil_dir,
            element_size=self.element_size,
            n_samples=self.n_samples,
            airfoil_spacing=self.airfoil_spacing,
            include_webs=self.include_webs,
        )
        base_mesh = self._blade_generator.generate(renumber=None, verbose=verbose)

        if verbose:
            print(
                f"  Base blade: {base_mesh.node_count} nodes, {base_mesh.elements_count} elements"
            )

        actual_hub_radius, hub_radius_source = self._resolve_hub_radius()
        if verbose:
            if hub_radius_source == "explicit":
                print(f"  Using explicit hub radius: {actual_hub_radius}")
            else:
                print(
                    "  Using hub radius from blade definition "
                    f"({hub_radius_source}): {actual_hub_radius}"
                )

        if verbose:
            print(f"\nGenerating rotor mesh with {self.n_blades} blades...")

        for i in range(self.n_blades):
            if verbose:
                print(f"  Processing blade {i + 1}/{self.n_blades}...")

            # Create a copy of the base blade mesh with renamed sets
            blade_mesh = self._create_blade_copy(base_mesh, i)

            # Apply radial translation (along Z axis) to offset blade from center
            if actual_hub_radius > 0:
                blade_mesh.translate_mesh(vector=(0, 0, 1), distance=actual_hub_radius)
                if verbose:
                    print(f"    Translated by Z={actual_hub_radius}")

            # Apply rotation around Y axis (rotor rotation axis) to distribute blades
            angle_rad = i * 2 * np.pi / self.n_blades
            if angle_rad > 0:
                blade_mesh.rotate_mesh(axis=(0, 1, 0), angle=angle_rad)
                if verbose:
                    print(f"    Rotated by {np.degrees(angle_rad):.1f}°")

            # Merge transformed blade mesh into rotor mesh
            self._merge_mesh(rotor_mesh, blade_mesh)

            if verbose:
                print(
                    f"    Rotor now has: {rotor_mesh.node_count} nodes, {rotor_mesh.elements_count} elements"
                )

        # Renumber mesh if requested
        if renumber is not None:
            if verbose:
                print(f"\nRenumbering mesh using {renumber} algorithm...")
            rotor_mesh.renumber_mesh(algorithm=renumber)

        if verbose:
            print(
                f"\nRotor mesh generated: {rotor_mesh.node_count} nodes, {rotor_mesh.elements_count} elements"
            )

        return rotor_mesh

    def _resolve_hub_radius(self) -> tuple[float, str]:
        """Resolve the radial blade offset using explicit or NuMAD geometry."""
        if self._blade_generator is None:
            raise RuntimeError(
                "Rotor mesh has not been initialized. Generate the blade mesh first."
            )

        definition = self._blade_generator.numad_blade.definition
        if self.hub_diameter is not None:
            definition.hub_diameter = float(self.hub_diameter)
        if self.rotor_diameter is not None:
            definition.rotor_diameter = float(self.rotor_diameter)

        try:
            return definition.resolve_hub_radius(override=self.hub_radius)
        except ValueError as exc:
            if (
                self.excel_file is not None
                and self.hub_radius is None
                and self.hub_diameter is None
                and self.rotor_diameter is None
            ):
                raise ValueError(
                    "Unable to determine hub radius from Excel blade definition. "
                    "Pass hub_radius explicitly or provide hub_diameter/rotor_diameter "
                    "when exporting a RotorMesh from Excel."
                ) from exc
            raise

    @classmethod
    def blade_part_set_name(cls, blade_index: int) -> str:
        """Return the aggregate set name used for a full rotor blade."""
        return f"{cls.BLADE_PART_SET_PREFIX}{blade_index + 1}"

    def _create_blade_copy(self, base_mesh: "MeshModel", blade_index: int) -> "MeshModel":
        """Create a copy of the base blade mesh with renamed sets."""
        from aeroelast.core.mesh.entities import ElementSet
        from aeroelast.core.mesh.model import MeshModel

        blade_mesh = MeshModel()
        node_map = {}

        # Copy nodes
        for node in base_mesh.nodes:
            new_node = Node([node.x, node.y, node.z])
            blade_mesh.add_node(new_node)
            node_map[node.id] = new_node

        # Copy elements
        element_map = {}
        for element in base_mesh.elements:
            new_nodes = [node_map[n.id] for n in element.nodes]
            new_element = MeshElement(new_nodes, element.element_type)
            blade_mesh.add_element(new_element)
            element_map[element.id] = new_element

        # Copy node sets with blade identifier
        for name, node_set in base_mesh.node_sets.items():
            new_set_name = f"{name}_blade_{blade_index + 1}"
            new_set_nodes = {node_map[n_id] for n_id in node_set.node_ids}
            blade_mesh.add_node_set(NodeSet(new_set_name, new_set_nodes))

        # Copy element sets with blade identifier
        for name, element_set in base_mesh.element_sets.items():
            new_set_name = f"{name}_blade_{blade_index + 1}"
            new_set_elements = {element_map[e.id] for e in element_set.elements}
            blade_mesh.add_element_set(ElementSet(new_set_name, new_set_elements))

        # Aggregate sets allow exporting an entire blade without reverse-engineering
        # the fine-grained material/region set naming.
        blade_part_name = self.blade_part_set_name(blade_index)
        blade_mesh.add_node_set(NodeSet(blade_part_name, set(node_map.values())))
        blade_mesh.add_element_set(ElementSet(blade_part_name, set(element_map.values())))

        return blade_mesh

    def _merge_mesh(self, target_mesh: "MeshModel", source_mesh: "MeshModel") -> None:
        """Merge a source mesh into the target mesh."""
        from aeroelast.core.mesh.entities import ElementSet

        node_map = {}

        for node in source_mesh.nodes:
            new_node = Node([node.x, node.y, node.z])
            target_mesh.add_node(new_node)
            node_map[node.id] = new_node

        element_map = {}
        for element in source_mesh.elements:
            new_nodes = [node_map[n.id] for n in element.nodes]
            new_element = MeshElement(new_nodes, element.element_type)
            target_mesh.add_element(new_element)
            element_map[element.id] = new_element

        for name, node_set in source_mesh.node_sets.items():
            new_set_nodes = {node_map[n_id] for n_id in node_set.node_ids}
            target_mesh.add_node_set(NodeSet(name, new_set_nodes))

        for name, element_set in source_mesh.element_sets.items():
            new_set_elements = {element_map[e.id] for e in element_set.elements}
            target_mesh.add_element_set(ElementSet(name, new_set_elements))

    @property
    def numad_mesh_data(self):
        """Access raw numad mesh data (only available after generate())."""
        if self._blade_generator is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate() first.")
        return self._blade_generator.numad_mesh_data

    @property
    def numad_blade(self):
        """Access the numad blade object (only available after generate())."""
        if self._blade_generator is None:
            raise RuntimeError("Mesh has not been generated yet. Call generate() first.")
        return self._blade_generator.numad_blade

    @classmethod
    def create_from_yaml(
        cls,
        yaml_file: str,
        n_blades: int,
        hub_radius: float = None,
        element_size: float = 0.1,
        n_samples: int = 300,
        renumber: str = None,
        include_webs: bool = True,
    ) -> "MeshModel":
        """
        Convenience method to create rotor mesh directly from YAML file.

        Parameters
        ----------
        yaml_file : str
            Path to the blade YAML definition file
        n_blades : int
            Number of blades in the rotor
        hub_radius : float, optional
            Radial distance from rotation axis to blade root
        element_size : float, optional
            Target element size (default: 0.1)
        n_samples : int, optional
            Number of samples for airfoil discretization (default: 300)
        renumber : str, optional
            Renumbering algorithm ("simple", "rcm", or None)

        Returns
        -------
        MeshModel
            The generated rotor mesh model
        """
        return cls(
            yaml_file=yaml_file,
            n_blades=n_blades,
            hub_radius=hub_radius,
            element_size=element_size,
            n_samples=n_samples,
            include_webs=include_webs,
        ).generate(renumber=renumber)


class RotorHubMesh:
    """
    Generates a CFD-oriented hub surface mesh for a wind turbine rotor.

    The generated geometry uses Gmsh OCC booleans to fuse:
    - a central hub cylinder aligned with the rotor axis (Y),
    - one radial connector cylinder per blade, and
    - a spherical nose on the front face so the hub is not flat.

    The hub dimensions are derived from the blade definition when not given
    explicitly. In particular, the default hub length is slightly larger than
    the diameter of the first blade profile.
    """

    SURFACE_SET_NAME = "hub_surface"

    def __init__(
        self,
        yaml_file: str | None = None,
        n_blades: int = 3,
        hub_radius: float | None = None,
        hub_diameter: float | None = None,
        rotor_diameter: float | None = None,
        element_size: float = 0.1,
        n_samples: int = 300,
        excel_file: str | None = None,
        airfoil_dir: str | None = None,
        airfoil_spacing: str = "cosine",
        hub_length: float | None = None,
        connector_radius: float | None = None,
        nose_radius: float | None = None,
        hub_length_factor: float = 1.1,
        connector_radius_factor: float = 0.5,
        connector_length_factor: float = 1.01,
    ):
        self.yaml_file = yaml_file
        self.excel_file = excel_file
        self.airfoil_dir = airfoil_dir
        self.n_blades = n_blades
        self.hub_radius = hub_radius
        self.hub_diameter = hub_diameter
        self.rotor_diameter = rotor_diameter
        self.element_size = element_size
        self.n_samples = n_samples
        self.airfoil_spacing = airfoil_spacing

        self.hub_length = hub_length
        self.connector_radius = connector_radius
        self.nose_radius = nose_radius
        self.hub_length_factor = hub_length_factor
        self.connector_radius_factor = connector_radius_factor
        self.connector_length_factor = connector_length_factor

        self._numad_blade = None

    def generate(self, renumber: str | None = None, verbose: bool = True) -> "MeshModel":
        """Generate the fused hub surface mesh."""
        from aeroelast.core.mesh.model import MeshModel

        self._load_blade_definition(verbose=verbose)

        actual_hub_radius, hub_radius_source = self._resolve_hub_radius()
        root_profile = self._first_profile_coords()
        root_diameter = self._estimate_first_profile_diameter(root_profile)
        hub_length, connector_length, connector_radius, nose_radius = (
            self._resolve_hub_geometry_parameters(actual_hub_radius, root_diameter)
        )

        if verbose:
            print("\nGenerating rotor hub mesh...")
            if hub_radius_source == "explicit":
                print(f"  Using explicit hub radius: {actual_hub_radius}")
            else:
                print(
                    "  Using hub radius from blade definition "
                    f"({hub_radius_source}): {actual_hub_radius}"
                )
            print(f"  First profile diameter: {root_diameter:.3f} m")
            print(f"  Hub length: {hub_length:.3f} m")
            print(f"  Connector length: {connector_length:.3f} m")
            print(f"  Connector radius: {connector_radius:.3f} m")
            print(f"  Nose radius: {nose_radius:.3f} m")

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("rotor_hub")

            self._create_geometry(
                hub_radius=actual_hub_radius,
                hub_length=hub_length,
                connector_length=connector_length,
                connector_radius=connector_radius,
                nose_radius=nose_radius,
            )
            self._configure_mesh(root_diameter=root_diameter)

            gmsh.model.occ.synchronize()

            surface_tags = [tag for dim, tag in gmsh.model.getEntities(2)]
            if surface_tags:
                gmsh.model.addPhysicalGroup(2, surface_tags, name=self.SURFACE_SET_NAME)

            gmsh.model.mesh.generate(2)

            mesh_model = self._create_mesh_model(MeshModel)
            if renumber is not None:
                mesh_model.renumber_mesh(algorithm=renumber)

            if verbose:
                print(
                    f"  Rotor hub mesh generated: {mesh_model.node_count} nodes, "
                    f"{mesh_model.elements_count} elements"
                )

            return mesh_model
        finally:
            gmsh.finalize()

    def _load_blade_definition(self, verbose: bool = True) -> None:
        """Load the blade definition and update its geometry without meshing."""
        try:
            from aeroelast.models.blade.numad import Blade as numadBlade
        except ImportError as exc:
            raise ImportError(
                "RotorHubMesh requires the numad module. "
                "Make sure aeroelast.models.blade.numad is available."
            ) from exc

        if verbose:
            print("  Loading blade definition for hub geometry...")

        self._numad_blade = numadBlade()
        if self.excel_file:
            self._numad_blade.read_excel(self.excel_file, airfoil_dir=self.airfoil_dir)
        elif self.yaml_file:
            self._numad_blade.read_yaml(self.yaml_file)
        else:
            raise ValueError("RotorHubMesh requires either yaml_file or excel_file")

        for station in self._numad_blade.definition.stations:
            station.airfoil.resample(n_samples=self.n_samples, spacing=self.airfoil_spacing)

        self._numad_blade.update_blade()

    def _resolve_hub_radius(self) -> tuple[float, str]:
        """Resolve the hub radius using the same priority as RotorMesh."""
        if self._numad_blade is None:
            raise RuntimeError("Blade definition has not been initialized yet.")

        definition = self._numad_blade.definition
        if self.hub_diameter is not None:
            definition.hub_diameter = float(self.hub_diameter)
        if self.rotor_diameter is not None:
            definition.rotor_diameter = float(self.rotor_diameter)

        try:
            return definition.resolve_hub_radius(override=self.hub_radius)
        except ValueError as exc:
            if (
                self.excel_file is not None
                and self.hub_radius is None
                and self.hub_diameter is None
                and self.rotor_diameter is None
            ):
                raise ValueError(
                    "Unable to determine hub radius from Excel blade definition. "
                    "Pass hub_radius explicitly or provide hub_diameter/rotor_diameter "
                    "when exporting a RotorHubMesh from Excel."
                ) from exc
            raise

    def _first_profile_coords(self) -> np.ndarray:
        """Return XYZ coordinates of the first blade profile from the updated geometry."""
        if self._numad_blade is None:
            raise RuntimeError("Blade definition has not been initialized yet.")

        coords = np.asarray(self._numad_blade.geometry.coordinates, dtype=float)
        if coords.ndim != 3 or coords.shape[2] == 0:
            raise ValueError("Blade geometry does not contain any spanwise stations.")
        return coords[:, :, 0]

    @staticmethod
    def _estimate_first_profile_diameter(profile_coords: np.ndarray) -> float:
        """Estimate the first-profile diameter from its XY footprint."""
        profile_coords = np.asarray(profile_coords, dtype=float)
        if profile_coords.ndim != 2 or profile_coords.shape[0] < 3 or profile_coords.shape[1] < 2:
            raise ValueError("At least 3 profile points with XY coordinates are required.")

        profile_xy = profile_coords[:, :2]
        centroid = profile_xy.mean(axis=0)
        radii = np.linalg.norm(profile_xy - centroid, axis=1)
        diameter = 2.0 * float(radii.max())
        if diameter <= 0.0:
            raise ValueError("First profile diameter must be positive.")
        return diameter

    def _resolve_hub_geometry_parameters(
        self,
        hub_radius: float,
        root_diameter: float,
    ) -> tuple[float, float, float, float]:
        """Resolve the fused hub geometry dimensions from defaults and overrides."""
        hub_length = (
            float(self.hub_length)
            if self.hub_length is not None
            else float(self.hub_length_factor * root_diameter)
        )
        # The connector cylinder starts at the rotor center, so its total length
        # should only slightly exceed the hub radius to guarantee overlap.
        connector_length = float(self.connector_length_factor * hub_radius)
        connector_radius = (
            float(self.connector_radius)
            if self.connector_radius is not None
            else float(min(self.connector_radius_factor * root_diameter, 0.95 * hub_radius))
        )
        nose_radius = float(self.nose_radius) if self.nose_radius is not None else float(hub_radius)

        if hub_length <= 0.0:
            raise ValueError("Hub length must be positive.")
        if connector_length <= 0.0:
            raise ValueError("Connector length must be positive.")
        if connector_radius <= 0.0:
            raise ValueError("Connector radius must be positive.")
        if nose_radius <= 0.0:
            raise ValueError("Nose radius must be positive.")

        return hub_length, connector_length, connector_radius, nose_radius

    def _create_geometry(
        self,
        hub_radius: float,
        hub_length: float,
        connector_length: float,
        connector_radius: float,
        nose_radius: float,
    ) -> None:
        """Create and fuse the central hub, blade connectors, and front nose."""
        hub_cylinder = gmsh.model.occ.addCylinder(
            0.0,
            -0.5 * hub_length,
            0.0,
            0.0,
            hub_length,
            0.0,
            hub_radius,
        )
        front_nose = gmsh.model.occ.addSphere(0.0, -0.5 * hub_length, 0.0, nose_radius)

        tools = [(3, front_nose)]
        for blade_idx in range(self.n_blades):
            angle_rad = blade_idx * 2.0 * np.pi / self.n_blades
            direction = np.array([np.sin(angle_rad), 0.0, np.cos(angle_rad)], dtype=float)
            connector = gmsh.model.occ.addCylinder(
                0.0,
                0.0,
                0.0,
                *(direction * connector_length),
                connector_radius,
            )
            tools.append((3, connector))

        gmsh.model.occ.fuse([(3, hub_cylinder)], tools, removeObject=True, removeTool=True)

    def _configure_mesh(self, root_diameter: float) -> None:
        """Configure a reasonable surface mesh size for the hub geometry."""
        surface_size = min(float(self.element_size), 0.25 * root_diameter)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", surface_size)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", surface_size)
        gmsh.option.setNumber("Mesh.Algorithm", 6)

    def _create_mesh_model(self, MeshModelClass) -> "MeshModel":
        """Convert the meshed Gmsh hub surfaces into a MeshModel."""
        mesh_model = MeshModelClass()

        elem_types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(2, -1)
        used_node_tags = set()
        parsed_elements = []

        for i, elem_type in enumerate(elem_types):
            props = gmsh.model.mesh.getElementProperties(elem_type)
            num_nodes = props[3]
            current_elem_tags = np.asarray(elem_tags_list[i], dtype=np.int64)
            current_elem_node_tags = np.asarray(elem_node_tags_list[i], dtype=np.int64).reshape(
                -1, num_nodes
            )
            parsed_elements.append((current_elem_tags, current_elem_node_tags))
            for tags in current_elem_node_tags:
                used_node_tags.update(tags)

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.asarray(coords, dtype=float).reshape(-1, 3)

        node_map = {}
        for tag, coord in zip(node_tags, coords):
            if tag not in used_node_tags:
                continue
            node = Node(coord, geometric_node=True)
            mesh_model.add_node(node)
            node_map[int(tag)] = node

        for elem_tags, current_elem_node_tags in parsed_elements:
            del elem_tags
            for nodes_gmsh in current_elem_node_tags:
                nodes = [node_map[int(tag)] for tag in nodes_gmsh]
                element_type = ELEMENT_NODES_MAP.get(len(nodes))
                if element_type is None:
                    continue
                mesh_model.add_element(MeshElement(nodes=nodes, element_type=element_type))

        all_nodes = set(mesh_model.nodes)
        all_elements = set(mesh_model.elements)
        mesh_model.add_node_set(NodeSet(name="all", nodes=all_nodes))
        mesh_model.add_node_set(NodeSet(name=self.SURFACE_SET_NAME, nodes=all_nodes))
        mesh_model.add_element_set(ElementSet(name=self.SURFACE_SET_NAME, elements=all_elements))

        return mesh_model


class CylindricalSurfaceMesh:
    """
    Generates a structured cylindrical shell mesh using Gmsh.

    Coordinate system:
    - Cylinder axis aligns with Z-axis.
    - Base center at (0, 0, 0).
    - Surface generated at radius R.
    - Angle measures from X-axis towards Y-axis.

    Attributes
    ----------
    radius : float
        Cylinder radius
    length : float
        Cylinder axial length (Z direction)
    angle_deg : float
        Circumferential span in degrees (default: 90.0)
    nx : int
        Number of circumferential elements
    ny : int
        Number of axial elements
    quadratic : bool
        Use quadratic elements (default: False)
    triangular : bool
        Use triangular elements instead of quads (default: False)
    distorted : bool
        Apply Ko2017 ratio-based mesh distortion (default: False)
    """

    def __init__(
        self,
        radius: float,
        length: float,
        nx: int,
        ny: int,
        angle_deg: float = 90.0,
        quadratic: bool = False,
        triangular: bool = False,
        distorted: bool = False,
    ):
        self.radius = radius
        self.length = length
        self.nx = nx
        self.ny = ny
        self.angle_deg = angle_deg
        self.quadratic = quadratic
        self.triangular = triangular
        self.distorted = distorted

    def generate(self) -> "MeshModel":
        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("cylinder")

            self._create_geometry()
            self._configure_mesh()

            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)

            return self._create_mesh_model(MeshModel)
        finally:
            gmsh.finalize()

    def _create_geometry(self):
        """Create cylindrical patch geometry"""
        R = self.radius
        L = self.length
        theta = np.radians(self.angle_deg)

        # Determine number of segments needed (arcs < 180 deg, prefer <= 90)
        n_segments = int(np.ceil(self.angle_deg / 90.0))
        d_theta = theta / n_segments

        # Center points at different Z levels (for axis)
        c_bottom = gmsh.model.geo.addPoint(0, 0, 0)
        c_top = gmsh.model.geo.addPoint(0, 0, L)

        points_bottom = []
        points_top = []

        # Create points along the arc
        for i in range(n_segments + 1):
            ang = i * d_theta
            x = R * np.cos(ang)
            y = R * np.sin(ang)

            points_bottom.append(gmsh.model.geo.addPoint(x, y, 0))
            points_top.append(gmsh.model.geo.addPoint(x, y, L))

        self.surfaces = []
        self.curves_bottom = []
        self.curves_top = []
        self.lines_vertical = []

        # First vertical line
        last_vertical = gmsh.model.geo.addLine(points_bottom[0], points_top[0])
        self.lines_vertical.append(last_vertical)

        # Create surfaces segment by segment
        for i in range(n_segments):
            # Arcs
            arc_bot = gmsh.model.geo.addCircleArc(points_bottom[i], c_bottom, points_bottom[i + 1])
            arc_top = gmsh.model.geo.addCircleArc(points_top[i], c_top, points_top[i + 1])

            # Next vertical line
            next_vertical = gmsh.model.geo.addLine(points_bottom[i + 1], points_top[i + 1])

            # Loop
            # Orientation: bot -> right -> -top -> -left
            loop = gmsh.model.geo.addCurveLoop([arc_bot, next_vertical, -arc_top, -last_vertical])
            surf = gmsh.model.geo.addSurfaceFilling([loop])

            self.surfaces.append(surf)
            self.curves_bottom.append(arc_bot)
            self.curves_top.append(arc_top)
            self.lines_vertical.append(next_vertical)

            last_vertical = next_vertical

        # Store boundaries for physical groups
        self.edge_bottom = self.curves_bottom
        self.edge_top = self.curves_top
        self.edge_left = [self.lines_vertical[0]]
        self.edge_right = [self.lines_vertical[-1]]

        self._add_physical_groups()

    def _add_physical_groups(self):
        gmsh.model.addPhysicalGroup(1, self.edge_bottom, name="bottom")
        gmsh.model.addPhysicalGroup(1, self.edge_top, name="top")
        gmsh.model.addPhysicalGroup(1, self.edge_left, name="left")
        gmsh.model.addPhysicalGroup(1, self.edge_right, name="right")
        gmsh.model.addPhysicalGroup(2, self.surfaces, name="surface")

    def _configure_mesh(self):
        """Set transfinite meshing"""
        # Distribute elements per segment
        n_segments = len(self.surfaces)
        nx_seg = max(1, self.nx // n_segments)
        # Remainder distribution handled simply here (might lose 1-2 elements if not divisible)
        # For benchmarks usually power of 2, so it's fine.

        for i, surf in enumerate(self.surfaces):
            # Curves order in loop: bot, right, top, left
            # Transfinite curve needs specific points if lines are used
            # Here we just set counts on curves

            gmsh.model.geo.mesh.setTransfiniteCurve(self.curves_bottom[i], nx_seg + 1)
            gmsh.model.geo.mesh.setTransfiniteCurve(self.curves_top[i], nx_seg + 1)

            # Vertical lines are shared, set only once effectively
            # Left of current is self.lines_vertical[i]
            # Right of current is self.lines_vertical[i+1]
            gmsh.model.geo.mesh.setTransfiniteCurve(self.lines_vertical[i], self.ny + 1)
            gmsh.model.geo.mesh.setTransfiniteCurve(self.lines_vertical[i + 1], self.ny + 1)

            gmsh.model.geo.mesh.setTransfiniteSurface(surf)

        # Configure element type (triangular or quad)
        if self.triangular:
            gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for triangles
        else:
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
            gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for Quads

        if self.quadratic:
            gmsh.option.setNumber("Mesh.ElementOrder", 2)

    def _apply_ko2017_distortion(self, mesh_model: "MeshModel") -> None:
        """Apply Ko2017 ratio-based mesh distortion in-place.

        The distortion pattern uses monotonically increasing segment lengths
        with ratios 1:2:3:...:n for both parametric directions.
        """
        if not self.distorted:
            return

        def _ratio_positions(n: int) -> np.ndarray:
            """Monotone node positions in [0,1] with segment lengths 1:2:...:n."""
            seg = np.arange(1, n + 1, dtype=float)
            cum = np.concatenate([[0.0], np.cumsum(seg)])
            return cum / cum[-1]

        # Original uniform parametric positions
        thetas_uniform = np.linspace(0.0, np.radians(self.angle_deg), self.nx + 1)
        zs_uniform = np.linspace(0.0, self.length, self.ny + 1)

        # Distorted positions using Ko2017 ratio pattern
        theta_dist = _ratio_positions(self.nx) * np.radians(self.angle_deg)
        z_dist = _ratio_positions(self.ny) * self.length

        R = self.radius
        for node in mesh_model.nodes:
            x, y, z = node.coords
            # Convert to parametric (theta, z)
            th = float(np.arctan2(y, x))
            if th < 0:
                th += 2 * np.pi
            # Snap to nearest uniform level
            i_th = int(np.argmin(np.abs(thetas_uniform - th)))
            i_z = int(np.argmin(np.abs(zs_uniform - z)))
            # Apply distorted positions
            new_th = theta_dist[i_th]
            new_z = z_dist[i_z]
            node.coords[0] = float(R * np.cos(new_th))
            node.coords[1] = float(R * np.sin(new_th))
            node.coords[2] = float(new_z)

    def _create_mesh_model(self, MeshModelClass) -> "MeshModel":
        # Reuse logic from SquareShapeMesh or similar base
        # But since I cannot inherit easily from here due to structure, I'll copy-paste adaptation
        # Ideally this should be a mixin or base class method.
        # For now, minimal implementation to get nodes/elements.

        mesh_model = MeshModelClass()

        # 1. Get all elements first to find used nodes
        # Retrieve all 2D elements (shells) irrespective of the entity tag (-1)
        elem_types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(2, -1)
        used_node_tags = set()

        parsed_elements = []  # List of (et, current_elem_node_tags)

        for i, et in enumerate(elem_types):
            props = gmsh.model.mesh.getElementProperties(et)
            num_nodes = props[3]

            # elem_node_tags_list[i] is a flat array of node tags
            current_elem_node_tags = elem_node_tags_list[i].reshape(-1, num_nodes)
            parsed_elements.append((et, current_elem_node_tags))

            for tags in current_elem_node_tags:
                used_node_tags.update(tags)

        # 2. Get all nodes and filter
        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)
        node_map_gmsh = {}  # gmsh_tag -> Node object

        for tag, coord in zip(node_tags, coords, strict=False):
            if tag in used_node_tags:
                # Basic node creation
                node = Node(coord, geometric_node=False)
                mesh_model.add_node(node)
                node_map_gmsh[tag] = node

        # 3. Create elements
        for _et, current_elem_node_tags in parsed_elements:
            for nodes_gmsh in current_elem_node_tags:
                nodes = [node_map_gmsh[tag] for tag in nodes_gmsh]
                e_type = ELEMENT_NODES_MAP.get(len(nodes), ElementType.quad)
                mesh_model.add_element(MeshElement(nodes=nodes, element_type=e_type))

        # Create node sets from physical groups
        physical_groups = gmsh.model.getPhysicalGroups(1)
        for dim, tag in physical_groups:
            name = gmsh.model.getPhysicalName(dim, tag)
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            nodes_set = set()
            for e_tag in entities:
                nt, _, _ = gmsh.model.mesh.getNodes(dim, e_tag, includeBoundary=True)
                for n_tag in nt:
                    nodes_set.add(node_map_gmsh[n_tag])

            if nodes_set:
                mesh_model.add_node_set(NodeSet(name=name, nodes=nodes_set))

        # Apply distortion if requested
        self._apply_ko2017_distortion(mesh_model)

        return mesh_model


class HyperbolicParaboloidMesh(SquareShapeMesh):
    """
    Hyperbolic Paraboloid Shell Mesh Generator.
    Surface defined by z = c * (x^2 - y^2) over a square domain.

    Parameters
    ----------
    length : float
        Side length of the square domain
    c : float
        Curvature parameter for z = c * (x^2 - y^2)
    nx : int
        Number of elements in x-direction
    ny : int
        Number of elements in y-direction
    thickness : float, optional
        Shell thickness (stored for reference)
    triangular : bool
        Use triangular elements instead of quads (default: False)
    distorted : bool
        Apply Ko2017 ratio-based mesh distortion (default: False)
    """

    def __init__(self, length, c, nx, ny, thickness=None, triangular=False, distorted=False):
        super().__init__(width=length, height=length, nx=nx, ny=ny, triangular=triangular)
        self.c = c
        self.thickness = thickness
        self.distorted = distorted

    def _create_geometry(self):
        # Create base flat square first
        super()._create_geometry()
        pass

    def _apply_ko2017_distortion(self, mesh_model: "MeshModel") -> None:
        """Apply Ko2017 ratio-based mesh distortion in-place."""
        if not self.distorted:
            return

        def _ratio_positions(n: int) -> np.ndarray:
            seg = np.arange(1, n + 1, dtype=float)
            cum = np.concatenate([[0.0], np.cumsum(seg)])
            return cum / cum[-1]

        # Domain is [-width/2, width/2] x [0, height]
        x0, x1 = -self.width / 2, self.width / 2
        y0, y1 = 0.0, self.height

        x_uniform = np.linspace(x0, x1, self.nx + 1)
        y_uniform = np.linspace(y0, y1, self.ny + 1)
        x_dist = x0 + _ratio_positions(self.nx) * (x1 - x0)
        y_dist = y0 + _ratio_positions(self.ny) * (y1 - y0)

        for node in mesh_model.nodes:
            x, y, z = node.coords
            i_x = int(np.argmin(np.abs(x_uniform - x)))
            i_y = int(np.argmin(np.abs(y_uniform - y)))
            new_x = x_dist[i_x]
            new_y = y_dist[i_y]
            # Recalculate z with distorted x, y
            y_shifted = new_y - self.height / 2.0
            new_z = self.c * (new_x**2 - y_shifted**2)
            node.coords[0] = float(new_x)
            node.coords[1] = float(new_y)
            node.coords[2] = float(new_z)

    def generate(self) -> "MeshModel":
        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("hyperbolic_paraboloid")

            # 1. Create flat mesh
            super()._create_geometry()
            self._configure_mesh()

            gmsh.model.geo.synchronize()
            gmsh.model.mesh.generate(2)

            # 2. Warp nodes to z = c * (x^2 - y^2)
            node_tags, coords_flat, _ = gmsh.model.mesh.getNodes()

            for i in range(0, len(coords_flat), 3):
                x = coords_flat[i]
                y = coords_flat[i + 1]

                y_shifted = y - self.height / 2.0
                z_new = self.c * (x**2 - y_shifted**2)

                gmsh.model.mesh.setNode(node_tags[i // 3], [x, y, z_new], [])

            mesh_model = self._create_mesh_model(MeshModel)

            # Apply distortion if requested
            self._apply_ko2017_distortion(mesh_model)

            return mesh_model

        finally:
            gmsh.finalize()


class RaaschHookMesh:
    """
    Raasch Challenge Hook Mesh.

    Based on Knight (1997) and Ko et al. geometry:
    - Two circular arc segments with G1 continuity
    - R2 is the MAIN arc (large radius) starting from the clamped end
    - R1 is the TIP arc (small radius) at the free end
    - theta2 is the angle of the main arc (typically 150°)
    - theta1 is the angle of the tip arc (typically 60°)

    Geometry (in X-Z plane, extruded in Y):
    - Clamped end at origin A=(0,0,0)
    - Hook curves toward negative X (inward)
    - Main arc center C2 at (-R2, 0, 0)
    - Tip arc center C1 positioned for G1 continuity

    Parameters from Knight (1997):
    - R1 = 14 (tip arc radius)
    - R2 = 46 (main arc radius)
    - theta1 = 60° (tip arc angle)
    - theta2 = 150° (main arc angle)
    - width = 20 (extrusion in Y)
    """

    def __init__(self, width, R1, R2, nx, ny, angle1_deg=60, angle2_deg=150, triangular=False):
        self.width = width
        self.R1 = R1  # Tip arc radius (small)
        self.R2 = R2  # Main arc radius (large)
        self.nx = nx  # Elements along width (Y direction)
        self.ny = ny  # Elements along curve (total for both arcs)
        self.angle1 = math.radians(angle1_deg)  # Tip arc angle
        self.angle2 = math.radians(angle2_deg)  # Main arc angle
        self.triangular = triangular  # Use triangular elements

    def generate(self) -> "MeshModel":
        import math

        from aeroelast.core.mesh.model import MeshModel

        gmsh = _import_gmsh()
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("raasch_hook")

            # ============================================================
            # Geometry: Hook in X-Z plane, extruded in Y
            # Based on the figure from Ko et al. / Knight (1997)
            # ============================================================

            # Point A: Origin (clamped end)
            # Tangent at A is vertical (along +Z axis)
            p0 = gmsh.model.geo.addPoint(0, 0, 0)

            # ---------------------------------------------------------
            # Main arc (Arc 2): R2, angle theta2
            # Center C2 at (-R2, 0, 0) so the arc starts at origin
            # and curves toward negative X
            # ---------------------------------------------------------
            c2 = gmsh.model.geo.addPoint(-self.R2, 0, 0)

            # End point of main arc P1
            # Arc sweeps counterclockwise from angle 0 to angle theta2
            # Point at angle 0 relative to C2: (0, 0, 0) ✓
            # Point at angle theta2: (-R2 + R2*cos(theta2), 0, R2*sin(theta2))
            x1 = -self.R2 + self.R2 * math.cos(self.angle2)
            z1 = self.R2 * math.sin(self.angle2)
            p1 = gmsh.model.geo.addPoint(x1, 0, z1)

            # Create main arc
            arc2 = gmsh.model.geo.addCircleArc(p0, c2, p1)

            # ---------------------------------------------------------
            # Tip arc (Arc 1): R1, angle theta1
            # For G1 continuity, C1 lies on the line C2-P1 extended
            # Direction from C2 to P1 (outward radial)
            # ---------------------------------------------------------
            radial_x = x1 - (-self.R2)  # = x1 + R2 = R2*cos(theta2)
            radial_z = z1 - 0  # = R2*sin(theta2)
            norm = math.sqrt(radial_x**2 + radial_z**2)
            dir_x = radial_x / norm
            dir_z = radial_z / norm

            # Center C1: For OPPOSITE curvature direction (true hook - curves back),
            # C1 = P1 + dir * R1 (center is "outside" the curve, opposite to main arc)
            cx1 = x1 + dir_x * self.R1
            cz1 = z1 + dir_z * self.R1
            c1 = gmsh.model.geo.addPoint(cx1, 0, cz1)

            # End point P2 (tip of hook)
            # Since we reversed curvature direction, we need to sweep in CLOCKWISE direction
            # Starting angle at P1 relative to C1 is theta2 + π (opposite direction)
            # End angle = (theta2 + π) - theta1 = theta2 + π - theta1
            theta_start_at_p1 = self.angle2 + math.pi  # P1 is now opposite from main arc
            theta_end = theta_start_at_p1 - self.angle1  # Sweep clockwise (negative)
            x2 = cx1 + self.R1 * math.cos(theta_end)
            z2 = cz1 + self.R1 * math.sin(theta_end)
            p2 = gmsh.model.geo.addPoint(x2, 0, z2)

            # Create tip arc
            arc1 = gmsh.model.geo.addCircleArc(p1, c1, p2)

            # ---------------------------------------------------------
            # Transfinite mesh: split ny proportionally by arc length
            # ---------------------------------------------------------
            L2 = self.R2 * self.angle2  # Main arc length
            L1 = self.R1 * self.angle1  # Tip arc length
            L_tot = L1 + L2

            # Distribute elements proportionally
            ny2 = max(2, int(self.ny * (L2 / L_tot)) + 1)
            ny1 = max(2, self.ny - ny2 + 2)

            gmsh.model.geo.mesh.setTransfiniteCurve(arc2, ny2)  # Main arc (from origin)
            gmsh.model.geo.mesh.setTransfiniteCurve(arc1, ny1)  # Tip arc

            # Extrude in Y with transfinite layers
            # Order: main arc first (from clamped end), then tip arc
            # recombine=True generates quads, recombine=False generates triangles
            recombine = not self.triangular
            extrusion = gmsh.model.geo.extrude(
                [(1, arc2), (1, arc1)], 0, self.width, 0, numElements=[self.nx], recombine=recombine
            )

            # Collect surfaces
            surfaces = []
            for e in extrusion:
                if e[0] == 2:
                    surfaces.append(e[1])

            gmsh.model.geo.synchronize()

            # Physical groups
            gmsh.model.addPhysicalGroup(2, surfaces, 1, "surface")

            gmsh.model.mesh.generate(2)

            from aeroelast.core.mesh.model import MeshModel

            return self._extract_mesh(MeshModel)
        finally:
            gmsh.finalize()

    def _extract_mesh(self, MeshModelClass):
        """Duplicated extraction logic with Raasch-specific sets"""

        # Reset Node ID counter
        from aeroelast.core.mesh.entities import Node

        Node._id_counter = 0
        mesh_model = MeshModelClass()

        elem_types, elem_tags_list, elem_node_tags_list = gmsh.model.mesh.getElements(2, -1)
        used_node_tags = set()
        parsed_elements = []
        for i, et in enumerate(elem_types):
            props = gmsh.model.mesh.getElementProperties(et)
            num_nodes = props[3]
            # Filter Quad (4 nodes) or Triangle (3 nodes)
            if num_nodes not in (3, 4):
                continue

            current_elem_node_tags = elem_node_tags_list[i].reshape(-1, num_nodes)
            parsed_elements.append((et, current_elem_node_tags))
            for tags in current_elem_node_tags:
                used_node_tags.update(tags)

        node_tags, coords, _ = gmsh.model.mesh.getNodes()
        coords = np.array(coords).reshape(-1, 3)
        node_map = {}
        for tag, coord in zip(node_tags, coords, strict=False):
            if tag in used_node_tags:
                node = Node(coord, geometric_node=False)
                mesh_model.add_node(node)
                node_map[tag] = node

        for _et, current_elem_node_tags in parsed_elements:
            for nodes_gmsh in current_elem_node_tags:
                nodes = [node_map[tag] for tag in nodes_gmsh]
                e_type = ELEMENT_NODES_MAP.get(len(nodes), ElementType.quad)
                mesh_model.add_element(MeshElement(nodes=nodes, element_type=e_type))

        # Create node sets (boundaries)
        root_nodes = []
        tip_nodes = []

        # Find tip and root by geometry:
        # - Tip (load applied): at the origin (x=0, z=0), varies in Y
        # - Root (clamped): at the end of the hook (minimum X), varies in Y
        tol = 1e-3

        # Find min X to identify root (clamped end)
        all_nodes = list(mesh_model.nodes)
        if all_nodes:
            min_x = min(n.x for n in all_nodes)

            for node in all_nodes:
                # Tip: near (0, *, 0) - origin, extrusion in Y - where load is applied
                if abs(node.x) < tol and abs(node.z) < tol:
                    tip_nodes.append(node)
                # Root: at minimum X (clamped end of hook)
                if abs(node.x - min_x) < tol:
                    root_nodes.append(node)

        mesh_model.add_node_set(NodeSet("root", set(root_nodes)))
        mesh_model.add_node_set(NodeSet("tip", set(tip_nodes)))

        return mesh_model


class SphericalSurfaceMesh:
    """
    Generates a structured spherical shell mesh (hemisphere or partial sphere).

    Coordinate system:
    - Sphere center at origin (0, 0, 0).
    - Z-axis is the polar axis.
    - θ (theta): polar angle from Z-axis (0° = north pole, 90° = equator)
    - φ (phi): azimuthal angle in X-Y plane from X-axis

    Parameters
    ----------
    radius : float
        Sphere radius
    theta_min_deg : float
        Minimum polar angle (hole at apex), default 0.0
    theta_max_deg : float
        Maximum polar angle, default 90.0 (equator)
    phi_max_deg : float
        Azimuthal span in degrees (default: 90.0 for 1/4 sphere)
    n_theta : int
        Number of elements in polar direction
    n_phi : int
        Number of elements in azimuthal direction
    triangular : bool
        Use triangular elements instead of quads (default: False)
    distorted : bool
        Apply Ko2017 ratio-based mesh distortion (default: False)

    Example
    -------
    # 1/4 hemisphere with 18° hole (MacNeal-Harder pinched hemisphere)
    gen = SphericalSurfaceMesh(
        radius=10.0,
        theta_min_deg=18.0,
        theta_max_deg=90.0,
        phi_max_deg=90.0,
        n_theta=16,
        n_phi=16
    )
    mesh = gen.generate()
    """

    def __init__(
        self,
        radius: float,
        n_theta: int,
        n_phi: int,
        theta_min_deg: float = 0.0,
        theta_max_deg: float = 90.0,
        phi_max_deg: float = 90.0,
        triangular: bool = False,
        distorted: bool = False,
    ):
        self.radius = radius
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.theta_min = np.radians(theta_min_deg)
        self.theta_max = np.radians(theta_max_deg)
        self.phi_max = np.radians(phi_max_deg)
        self.theta_min_deg = theta_min_deg
        self.theta_max_deg = theta_max_deg
        self.phi_max_deg = phi_max_deg
        self.triangular = triangular
        self.distorted = distorted

    def _apply_ko2017_distortion(self, mesh_model: "MeshModel") -> None:
        """Apply Ko2017 ratio-based mesh distortion in-place."""
        if not self.distorted:
            return

        def _ratio_positions(n: int) -> np.ndarray:
            seg = np.arange(1, n + 1, dtype=float)
            cum = np.concatenate([[0.0], np.cumsum(seg)])
            return cum / cum[-1]

        R = self.radius
        thetas_uniform = np.linspace(self.theta_min, self.theta_max, self.n_theta + 1)
        phis_uniform = np.linspace(0.0, self.phi_max, self.n_phi + 1)

        theta_dist = self.theta_min + _ratio_positions(self.n_theta) * (
            self.theta_max - self.theta_min
        )
        phi_dist = _ratio_positions(self.n_phi) * self.phi_max

        for node in mesh_model.nodes:
            x, y, z = node.coords
            # Convert to spherical
            th = float(np.arccos(np.clip(z / R, -1.0, 1.0)))
            ph = float(np.arctan2(y, x))
            if ph < 0:
                ph += 2 * np.pi

            # Snap to nearest uniform level
            i_th = int(np.argmin(np.abs(thetas_uniform - th)))
            i_ph = int(np.argmin(np.abs(phis_uniform - ph)))

            # Apply distorted positions
            new_th = theta_dist[i_th]
            new_ph = phi_dist[i_ph]

            node.coords[0] = float(R * np.sin(new_th) * np.cos(new_ph))
            node.coords[1] = float(R * np.sin(new_th) * np.sin(new_ph))
            node.coords[2] = float(R * np.cos(new_th))

    def generate(self) -> "MeshModel":
        from aeroelast.core.mesh.entities import ElementType, MeshElement, Node, NodeSet
        from aeroelast.core.mesh.model import MeshModel

        Node._id_counter = 0
        mesh_model = MeshModel()

        R = self.radius
        n_theta = self.n_theta
        n_phi = self.n_phi

        # Create nodes
        node_grid = {}  # (j, i) -> Node

        for j in range(n_theta + 1):
            theta = self.theta_min + j * (self.theta_max - self.theta_min) / n_theta
            for i in range(n_phi + 1):
                phi = i * self.phi_max / n_phi

                x = R * np.sin(theta) * np.cos(phi)
                y = R * np.sin(theta) * np.sin(phi)
                z = R * np.cos(theta)

                node = Node([x, y, z], geometric_node=False)
                mesh_model.add_node(node)
                node_grid[(j, i)] = node

        # Create elements
        for j in range(n_theta):
            for i in range(n_phi):
                n0 = node_grid[(j, i)]
                n1 = node_grid[(j, i + 1)]
                n2 = node_grid[(j + 1, i + 1)]
                n3 = node_grid[(j + 1, i)]

                if self.triangular:
                    # Split quad into two triangles
                    elem1 = MeshElement(nodes=[n0, n1, n2], element_type=ElementType.triangle)
                    elem2 = MeshElement(nodes=[n0, n2, n3], element_type=ElementType.triangle)
                    mesh_model.add_element(elem1)
                    mesh_model.add_element(elem2)
                else:
                    elem = MeshElement(nodes=[n0, n1, n2, n3], element_type=ElementType.quad)
                    mesh_model.add_element(elem)

        # Create boundary node sets
        # theta_min edge (hole/apex)
        theta_min_nodes = [node_grid[(0, i)] for i in range(n_phi + 1)]
        mesh_model.add_node_set(NodeSet("theta_min", theta_min_nodes))

        # theta_max edge (equator)
        theta_max_nodes = [node_grid[(n_theta, i)] for i in range(n_phi + 1)]
        mesh_model.add_node_set(NodeSet("equator", theta_max_nodes))

        # phi=0 edge
        phi_0_nodes = [node_grid[(j, 0)] for j in range(n_theta + 1)]
        mesh_model.add_node_set(NodeSet("phi_0", phi_0_nodes))

        # phi=phi_max edge
        phi_max_nodes = [node_grid[(j, n_phi)] for j in range(n_theta + 1)]
        mesh_model.add_node_set(NodeSet("phi_max", phi_max_nodes))

        # Apply distortion if requested
        self._apply_ko2017_distortion(mesh_model)

        return mesh_model
