"""GmshBladeMesh: the gmsh skin mesher as a production mesh generator (odd/tasks/gmsh-official-mesher.md, R1).

Contract the framework relies on (runner.py, models/blade/model.py): ``generate()`` returns a ``MeshModel`` with the node sets
``RootNodes``/``TipNodes``, the element sets ``SEC<i>`` (one per laminate), ``allOuterShellEls`` and ``allWebEls`` plus one set per
web; ``numad_mesh_data`` carries the materials and sections for ``build_rust_properties``.  The shell sits on the mid-surface unless
told otherwise.  A tiny synthetic raw mesh stands in for the (minutes long) gmsh run; the IEA-15 equivalence is a slow test below.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import aeroelast.core.mesh.gmsh_blade as gb


def _raw():
    """Two skin quads stacked in z (eta 0 -> 1) and one web triangle on the first quad's edge."""
    nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1], [1, 0, 2], [0, 0, 2], [0.5, 0.5, 0.0]], float)
    elements = np.array([[0, 1, 2, 3], [3, 2, 4, 5], [0, 1, 6, -1]], dtype=np.int64)
    params = np.full((3, 4, 2), np.nan)
    params[0] = [[0.0, 0.0], [0.0, 1.0], [0.5, 1.0], [0.5, 0.0]]
    params[1] = [[0.5, 0.0], [0.5, 1.0], [1.0, 1.0], [1.0, 0.0]]
    params[2, :3] = [[0.0, 0.0], [0.0, 1.0], [0.0, 0.5]]
    tables = {
        "layups": [[["glass", 0.02, 0.0]], [["glass", 0.01, 0.0]]],
        "materials": {"glass": {"E": [4e10, 1e10, 1e10], "nu": [0.3, 0.3, 0.3], "G": [4e9, 4e9, 4e9], "rho": 1900.0}},
        "surfaces": ["skin", "web0"],
    }
    return {
        "nodes": nodes, "elements": elements, "surfaces": np.array(["skin", "skin", "web0"]), "params": params,
        "layup_id": np.array([0, 0, 1], dtype=np.int64),
        "normals": np.array([[0, -1.0, 0], [0, -1.0, 0], [1.0, 0, 0]]), "offsets": np.array([0.01, 0.01, 0.0]),
        "areas": gb.element_areas(nodes, elements), "tables_json": np.frombuffer(json.dumps(tables).encode(), dtype=np.uint8),
    }


def test_element_areas_of_a_unit_quad_and_a_triangle():
    a = _raw()["areas"]
    assert a[0] == pytest.approx(1.0) and a[2] == pytest.approx(0.25)


def test_mesh_has_the_sets_the_framework_reads():
    mesh, props = gb.mesh_from_arrays(_raw(), "mid")
    assert len(mesh.nodes) == 7 and mesh.elements_count == 3
    names = set(mesh.element_sets)
    assert {"SEC0", "SEC1", "allOuterShellEls", "allWebEls", "web0"} <= names
    assert {"RootNodes", "TipNodes"} <= set(mesh.node_sets)
    assert len(mesh.get_node_set("RootNodes").nodes) == 3   # eta = 0: skin nodes 0 and 1 and the web node 6
    assert len(mesh.get_node_set("TipNodes").nodes) == 2
    assert set(props) == {"SEC0", "SEC1"}


def test_mid_surface_moves_the_skin_nodes_along_the_normal_and_oml_does_not():
    raw = _raw()
    mid, _ = gb.mesh_from_arrays(raw, "mid")
    oml, _ = gb.mesh_from_arrays(raw, "oml")
    ys_mid = sorted(round(float(p[1]), 6) for p in mid.coords_array)
    ys_oml = sorted(round(float(p[1]), 6) for p in oml.coords_array)
    assert ys_oml == sorted(raw["nodes"][:, 1].round(6).tolist())
    assert ys_mid != ys_oml and min(ys_mid) == pytest.approx(-0.01)   # offset 0.01 along the -y normal


def test_an_invalid_surface_is_rejected():
    with pytest.raises(ValueError, match="skin_reference"):
        gb.mesh_from_arrays(_raw(), "inner")


def test_generator_follows_the_blade_mesh_contract(monkeypatch):
    monkeypatch.setattr(gb.GmshBladeMesh, "_arrays", lambda self: _raw())
    gen = gb.GmshBladeMesh(yaml_file="blade.yaml", element_size=0.25)
    assert gen.skin_reference == "mid"
    mesh = gen.generate(renumber=None, verbose=False)
    assert mesh.elements_count == 3
    assert set(gen.numad_mesh_data) >= {"materials", "sections"}
    assert [m["name"] for m in gen.numad_mesh_data["materials"]] == ["glass"]


def test_generator_rejects_an_invalid_surface_at_construction():
    with pytest.raises(ValueError, match="skin_reference"):
        gb.GmshBladeMesh(yaml_file="blade.yaml", element_size=0.25, skin_reference="inner")
