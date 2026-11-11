"""The gmsh skin mesher as a production blade mesh generator.

``GmshBladeMesh`` has the contract of ``BladeMesh`` that the framework relies on (``generate()`` returns a ``MeshModel``;
``numad_mesh_data`` carries the materials and sections for ``build_rust_properties``) and builds the blade with
``numad.mesh_gen.gmsh_skin.generate_mesh_gmsh`` from a windIO 2.x YAML.  The mesh arrays -> ``MeshModel`` conversion is the one the
validation tests used (``tests/support/v2_blade.py``), moved here so production and tests share it.

Mesh sets: nodes ``RootNodes`` / ``TipNodes`` (minimum / maximum chart ``eta``), ``allOuterShellNods``, ``allWebNods``; elements one
``SEC<i>`` per laminate, ``allOuterShellEls``, ``allWebEls`` and one set per web.  The shell is built on the mid-surface
(``skin_reference="mid"``, the default and the configuration every structural validation anchor is measured on) or left on the outer
mould line (``"oml"``).
"""

from __future__ import annotations

import json
from typing import Any, Dict, Optional, Tuple

import numpy as np

SKIN_REFERENCES = ("mid", "oml")


def element_areas(x: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """3-D area of every element (quads padded with -1 for triangles)."""
    areas = np.empty(len(elements))
    for k, el in enumerate(elements):
        p = x[el[el >= 0]]
        cross = np.cross(p[2] - p[0], p[-1] - p[1]) / 2 if len(p) == 4 else np.cross(p[1] - p[0], p[2] - p[0]) / 2
        areas[k] = np.linalg.norm(cross)
    return areas


def raw_to_arrays(raw) -> Dict[str, Any]:
    """Plain arrays and JSON-able tables from a ``RawMesh`` with material stacks (what ``generate_mesh_gmsh`` returns)."""
    ne = len(raw.elements)
    params = np.full((ne, 4, 2), np.nan)
    for k in range(ne):
        p = np.asarray(raw.element_parameters[k], float)
        params[k, : len(p)] = p
    surfaces = [str(s) for s in raw.surfaces]
    layups, layup_index = [], {}
    layup_id = np.empty(ne, dtype=np.int64)
    normals = np.empty((ne, 3))
    offsets = np.empty(ne)
    materials: Dict[str, Any] = {}
    for k, stack in enumerate(raw.material_stacks):
        layers = []
        for layer in stack.layers:
            s = layer.source
            layers.append([s.material_name, round(float(s.thickness_m), 14), float(s.fiber_angle_deg)])
            if s.material_name not in materials:
                m = s.material
                tup = lambda v: [float(t) for t in np.atleast_1d(np.asarray(v, float))]  # noqa: E731
                materials[s.material_name] = {
                    "E": tup(m["E"]) if np.ndim(m["E"]) else float(m["E"]),
                    "nu": tup(m["nu"]) if np.ndim(m["nu"]) else float(m["nu"]),
                    "G": (tup(m["G"]) if np.ndim(m["G"]) else float(m["G"])) if "G" in m else None,
                    "rho": float(np.mean(np.asarray(m["rho"], float))),
                }
        key = json.dumps(layers)
        if key not in layup_index:
            layup_index[key] = len(layups)
            layups.append(layers)
        layup_id[k] = layup_index[key]
        normals[k] = np.asarray(stack.physical_normal, float)
        offsets[k] = float(stack.midplane_offset_m)
    nodes = np.asarray(raw.nodes, float)
    elements = np.asarray(raw.elements, dtype=np.int64)
    tables = {"layups": layups, "materials": materials, "surfaces": sorted(set(surfaces), key=surfaces.index)}
    return {
        "nodes": nodes,
        "elements": elements,
        "surfaces": np.array(surfaces),
        "params": params,
        "layup_id": layup_id,
        "normals": normals,
        "offsets": offsets,
        "areas": element_areas(nodes, elements),
        "tables_json": np.frombuffer(json.dumps(tables).encode(), dtype=np.uint8),
    }


def inset(arrays: Dict[str, Any]) -> np.ndarray:
    """Outer mould line -> mid-surface node coordinates (area-weighted mean mid-plane offset along the area-weighted normal)."""
    x = np.array(arrays["nodes"], float)
    nsum, osum, wsum = np.zeros_like(x), np.zeros(len(x)), np.zeros(len(x))
    for k, el in enumerate(arrays["elements"]):
        if arrays["surfaces"][k] != "skin":
            continue
        ids = el[el >= 0]
        a = arrays["areas"][k]
        nsum[ids] += a * arrays["normals"][k]
        osum[ids] += a * arrays["offsets"][k]
        wsum[ids] += a
    moved = wsum > 0
    norm = np.maximum(np.linalg.norm(nsum[moved], axis=1), 1e-300)[:, None]
    x[moved] += (osum[moved] / wsum[moved])[:, None] * nsum[moved] / norm
    return x


def numad_data(tables: Dict[str, Any], layup_names: list) -> Dict[str, Any]:
    """Materials and sections in the dictionary layout ``build_rust_properties`` reads."""
    materials = []
    for name, m in tables["materials"].items():
        elastic = {"E": m["E"], "nu": m["nu"]}
        if m["G"] is not None:
            elastic["G"] = m["G"]
        materials.append({"name": name, "density": m["rho"], "elastic": elastic})
    sections = [{"elementSet": name, "layup": layup} for name, layup in zip(layup_names, tables["layups"], strict=True)]
    return {"materials": materials, "sections": sections}


def mesh_from_arrays(arrays: Dict[str, Any], skin_reference: str = "mid", renumber: Optional[str] = "rcm"):
    """``(MeshModel, rust_properties)`` from the arrays of :func:`raw_to_arrays`."""
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node, NodeSet
    from aeroelast.core.mesh.model import MeshModel
    from aeroelast.models.blade.model import build_rust_properties

    if skin_reference not in SKIN_REFERENCES:
        raise ValueError(f"skin_reference must be one of {SKIN_REFERENCES}, got {skin_reference!r}")
    tables = json.loads(bytes(arrays["tables_json"]).decode())
    x = inset(arrays) if skin_reference == "mid" else np.array(arrays["nodes"], float)
    mesh = MeshModel()
    for p in x:
        mesh.add_node(Node(p))
    nodes = list(mesh.nodes)
    surfaces = [str(s) for s in arrays["surfaces"]]
    elements = []
    for el in arrays["elements"]:
        ids = [int(i) for i in el if i >= 0]
        e = MeshElement([nodes[i] for i in ids], ElementType.quad if len(ids) == 4 else ElementType.triangle)
        mesh.add_element(e)
        elements.append(e)

    layup_names = [f"SEC{i}" for i in range(len(tables["layups"]))]
    by_layup: Dict[int, set] = {}
    by_surface: Dict[str, set] = {}
    for k, e in enumerate(elements):
        by_layup.setdefault(int(arrays["layup_id"][k]), set()).add(e)
        by_surface.setdefault(surfaces[k], set()).add(e)
    for i, es in sorted(by_layup.items()):
        mesh.add_element_set(ElementSet(name=layup_names[i], elements=es))
    skin_els = by_surface.get("skin", set())
    web_els = set().union(*[v for s, v in by_surface.items() if s != "skin"]) if len(by_surface) > 1 else set()
    mesh.add_element_set(ElementSet(name="allOuterShellEls", elements=skin_els))
    if web_els:
        mesh.add_element_set(ElementSet(name="allWebEls", elements=web_els))
        for s, v in by_surface.items():
            if s != "skin":
                mesh.add_element_set(ElementSet(name=s, elements=v))

    eta = np.full(len(x), np.nan)
    for k, el in enumerate(arrays["elements"]):
        for j, i in enumerate(el[el >= 0]):
            eta[i] = arrays["params"][k, j, 0]
    lo, hi = np.nanmin(eta), np.nanmax(eta)
    mesh.add_node_set(NodeSet(name="RootNodes", nodes={nodes[i] for i in np.nonzero(eta == lo)[0]}))
    mesh.add_node_set(NodeSet(name="TipNodes", nodes={nodes[i] for i in np.nonzero(eta == hi)[0]}))
    mesh.add_node_set(NodeSet(name="allOuterShellNods", nodes={n for e in skin_els for n in e.nodes}))
    if web_els:
        mesh.add_node_set(NodeSet(name="allWebNods", nodes={n for e in web_els for n in e.nodes}))

    numad = numad_data(tables, layup_names)
    properties = build_rust_properties(numad)
    if renumber:
        mesh.renumber_mesh(algorithm=renumber)
    mesh.numad_mesh_data = numad
    return mesh, properties


class GmshBladeMesh:
    """Blade shell mesh from a windIO 2.x YAML with the structured gmsh skin mesher.

    Parameters
    ----------
    yaml_file : str
        windIO 2.x turbine/blade YAML.
    element_size : float
        Target element size [m]; the mesher accepts ``h_range`` (default 0.05-0.25 m; the converged size is per blade:
        IEA-15 0.25, IEA-22 0.177, see odd/tasks/gmsh-validation-15-22.md).
    skin_reference : {"mid", "oml"}
        Surface the shell is built on (default ``"mid"``).
    chord_divisions : int, optional
        N elements per chord (>= 4).  When given, ``element_size`` and ``h_range`` are not used: the chordwise size is c(eta)/N,
        the spanwise size ``span_aspect`` * c(eta)/N and every material band gets at least ``min_band_elements`` elements across.
        None (default) keeps the uniform ``element_size`` path.
    min_band_fraction : float, optional
        With ``chord_divisions``: the ``min_band_elements`` floor applies only to a band at least this fraction of c/N wide
        (default 1.0, the S3 rule; narrower bands get their rounded count, >= 1).  None = the floor on every band.
    """

    def __init__(self, yaml_file: str, element_size: float = 0.25, skin_reference: str = "mid",
                 h_range: Tuple[float, float] = (0.05, 0.25), chord_divisions: Optional[float] = None,
                 span_aspect: float = 1.0, min_band_elements: int = 2,
                 min_band_fraction: Optional[float] = 1.0) -> None:
        if skin_reference not in SKIN_REFERENCES:
            raise ValueError(f"skin_reference must be one of {SKIN_REFERENCES}, got {skin_reference!r}")
        self.yaml_file = yaml_file
        self.element_size = float(element_size)
        self.skin_reference = skin_reference
        self.h_range = tuple(h_range)
        self.chord_divisions = None if chord_divisions is None else float(chord_divisions)
        self.span_aspect = float(span_aspect)
        self.min_band_elements = int(min_band_elements)
        self.min_band_fraction = None if min_band_fraction is None else float(min_band_fraction)
        self.mesh = None
        self.numad_mesh_data: Dict[str, Any] = {}
        self.numad_blade = None

    def _arrays(self) -> Dict[str, Any]:
        """Run the gmsh mesher (minutes) and return the arrays of :func:`raw_to_arrays`."""
        from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import generate_mesh_gmsh
        from aeroelast.models.blade.numad.objects.blade import Blade

        blade = Blade()
        blade.read_yaml(str(self.yaml_file))
        self.numad_blade = blade
        chord = {} if self.chord_divisions is None else {
            "chord_divisions": self.chord_divisions, "span_aspect": self.span_aspect, "min_band_elements": self.min_band_elements}
        if chord:
            chord["min_band_fraction"] = self.min_band_fraction
        raw = generate_mesh_gmsh(blade.windio2_topology, self.element_size, with_materials=True, h_range=self.h_range, **chord)
        return raw_to_arrays(raw)

    def generate(self, renumber: Optional[str] = None, verbose: bool = False):
        """Build the ``MeshModel``; ``renumber`` (None or ``"rcm"``) as in ``BladeMesh.generate``."""
        arrays = self._arrays()
        mesh, _ = mesh_from_arrays(arrays, self.skin_reference, renumber=renumber)
        self.mesh = mesh
        self.numad_mesh_data = mesh.numad_mesh_data
        return mesh
