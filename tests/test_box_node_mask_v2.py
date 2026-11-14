"""``sectional.box_node_mask`` on windIO v2 meshes (SEC<n> sets) and the unchanged v1 reading (l3-outboard T8, D4)."""

import numpy as np

from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
from aeroelast.core.mesh.model import MeshModel
from aeroelast.postprocess.sectional import box_node_mask


def _strip(names):
    """Four quads in a row (x = 0..4, z = 0..1) of one named element set each; returns (mesh, elements)."""
    mesh = MeshModel()
    for k in (0, 1):
        for i in range(5):
            mesh.add_node(Node([float(i), 0.0, float(k)]))
    nodes = list(mesh.nodes)
    els = []
    for i in range(4):
        e = MeshElement([nodes[i], nodes[i + 1], nodes[6 + i], nodes[5 + i]], ElementType.quad)
        mesh.add_element(e)
        els.append(e)
    for name, e in zip(names, els, strict=True):
        if name:
            mesh.add_element_set(ElementSet(name=name, elements={e}))
    return mesh, els


def _nodes_of(mesh, mask):
    return sorted(int(i) for i in np.nonzero(mask)[0])


def test_v2_box_is_caps_and_webs_only():
    # e0 skin glass (SEC0), e1 skin carbon cap (SEC1), e2 skin glass_uni TE (SEC2), e3 web.
    mesh, els = _strip(["SEC0", "SEC1", "SEC2", None])
    mesh.add_element_set(ElementSet(name="allWebEls", elements={els[3]}))
    mesh.numad_mesh_data = {
        "materials": [
            {"name": "glass_triax", "elastic": {"E": [2.8e10, 1.3e10, 1.3e10]}},
            {"name": "CarbonUD", "elastic": {"E": [1.14e11, 8e9, 8e9]}},
            {"name": "glass_uni", "elastic": {"E": [4.2e10, 1.4e10, 1.4e10]}},
        ],
        "sections": [
            {"elementSet": "SEC0", "layup": [["glass_triax", 0.01, 0.0]]},
            {"elementSet": "SEC1", "layup": [["glass_triax", 0.01, 0.0], ["CarbonUD", 0.03, 0.0]]},
            {"elementSet": "SEC2", "layup": [["glass_uni", 0.01, 0.0]]},
        ],
    }
    mask = box_node_mask(mesh)
    expected = {mesh.node_id_to_index[n] for e in (els[1], els[3]) for n in e.node_ids}
    assert mask.any()
    assert set(_nodes_of(mesh, mask)) == expected  # cap element + web element nodes, nothing else


def test_v2_mesh_without_webs_or_data_gives_empty_mask():
    mesh, _ = _strip(["SEC0", "SEC1", "SEC2", "SEC3"])
    assert not box_node_mask(mesh).any()


def test_v1_names_keep_selecting_spar_and_sw_sets():
    mesh, els = _strip(["00_01_TE_SKIN", "00_01_SPAR_CAP", "00_01_TE_FLAT", "00_01_SW"])
    mask = box_node_mask(mesh)
    expected = {mesh.node_id_to_index[n] for e in (els[1], els[3]) for n in e.node_ids}
    assert set(_nodes_of(mesh, mask)) == expected
