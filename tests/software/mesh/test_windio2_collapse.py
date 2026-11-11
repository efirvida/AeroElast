"""Collapse of thin strips in windIO v2 raw meshes (synthetic flat strips, no loft, no solver)."""

import numpy as np
import pytest

from aeroelast.models.blade.numad.mesh_gen.windio2_collapse import collapse_thin_strips

ES = 0.25


def strip(widths, rows_above=1, dx=0.24):
    """Lines y=0 and y=w_i (thin band), then ``rows_above`` regular rows of height ES above the band line."""
    n = len(widths)
    nodes, ids = [], {}
    for i in range(n):
        for j, y in enumerate([0.0, widths[i]] + [widths[i] + ES * (r + 1) for r in range(rows_above)]):
            ids[i, j] = len(nodes)
            nodes.append((i * dx, y, 0.0))
    elems = []
    for i in range(n - 1):
        for j in range(1 + rows_above):
            elems.append([ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]])
    return np.array(nodes), np.array(elems), ids


def areas(nodes, elems):
    out = []
    for el in elems:
        p = nodes[el[el >= 0]]
        out.append(np.linalg.norm(np.cross(p[2] - p[0], p[-1] - p[1])) / 2 if len(p) == 4
                   else np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) / 2)
    return np.array(out)


def aspect(nodes, elems):
    out = []
    for el in elems:
        p = nodes[el[el >= 0]]
        ln = np.linalg.norm(np.roll(p, -1, 0) - p, axis=1)
        out.append(ln.max() / ln.min())
    return np.array(out)


def test_constant_width_band_is_removed_and_the_rest_stays_regular():
    nodes, elems, _ = strip([0.01] * 11)
    new_nodes, new_elems, src, rep = collapse_thin_strips(nodes, elems, max_edge_len=0.15 * ES)
    assert rep["dropped"] == 10 and rep["to_triangle"] == 0
    assert len(new_elems) == 10 and (new_elems >= 0).all()
    assert aspect(new_nodes, new_elems).max() < 1.1
    assert rep["max_node_shift_m"] <= 0.005 + 1e-12
    # the 10 surviving quads are the regular ones; the band (w=0.01) is removed
    assert np.allclose(areas(new_nodes, new_elems).sum(), 10 * 0.24 * (ES + 0.005), rtol=1e-9)  # merged line sits mid-band
    assert rep["removed_area_m2"] == pytest.approx(10 * 0.24 * 0.01, rel=1e-9)
    assert list(src) == [i for i in range(len(elems)) if i % 2 == 1]


def test_tapering_band_ends_in_triangles_and_nothing_inverts():
    w = [0.20, 0.16, 0.12, 0.08, 0.05, 0.03, 0.02, 0.01, 0.004, 0.001, 0.0005]
    nodes, elems, _ = strip(w)
    new_nodes, new_elems, src, rep = collapse_thin_strips(nodes, elems, max_edge_len=0.15 * ES)
    assert rep["to_triangle"] >= 1 and rep["dropped"] >= 1
    assert (areas(new_nodes, new_elems) > 0).all()
    normals = [np.cross(new_nodes[e[e >= 0]][1] - new_nodes[e[e >= 0]][0], new_nodes[e[e >= 0]][2] - new_nodes[e[e >= 0]][1]) for e in new_elems]
    assert all(n[2] > 0 for n in normals)
    assert aspect(new_nodes, new_elems).max() < aspect(nodes, elems).max()
    # elements that were not touched are bit-identical in shape
    touched = rep["to_triangle"] + rep["dropped"]
    assert len(new_elems) == len(elems) - rep["dropped"] and touched >= 2


def test_fixed_nodes_never_move():
    nodes, elems, ids = strip([0.01] * 11)
    fixed = [ids[i, 1] for i in range(11)]
    new_nodes, new_elems, _, rep = collapse_thin_strips(nodes, elems, max_edge_len=0.15 * ES, fixed=fixed)
    for i in range(11):
        assert any(np.allclose(p, nodes[ids[i, 1]]) for p in new_nodes)
    assert rep["max_node_shift_m"] == 0.0 or rep["max_node_shift_m"] <= 0.01 + 1e-12


def test_mesh_without_thin_edges_is_returned_unchanged():
    nodes, elems, _ = strip([ES] * 6)
    new_nodes, new_elems, src, rep = collapse_thin_strips(nodes, elems, max_edge_len=0.15 * ES)
    assert np.array_equal(new_nodes, nodes) and np.array_equal(new_elems, elems)
    assert list(src) == list(range(len(elems))) and rep["merged"] == 0


def test_padded_triangles_are_accepted():
    nodes, elems, _ = strip([ES] * 4)
    tri = np.array([[elems[0][0], elems[0][1], elems[0][2], -1]])
    new_nodes, new_elems, _, rep = collapse_thin_strips(nodes, np.vstack([elems, tri]), max_edge_len=0.15 * ES)
    assert new_elems.shape[1] == 4 and rep["merged"] == 0


def test_generate_mesh_collapse_is_opt_in_and_cuts_slivers_on_the_official_patch():
    from aeroelast.models.blade.numad.mesh_gen.windio2_mesh import generate_mesh
    from aeroelast.models.blade.numad.objects.blade import Blade
    from pathlib import Path

    yaml = Path(__file__).resolve().parents[2] / "IEA15MW" / "reference" / "upstream" / "IEA-15-240-RWT.yaml"
    blade = Blade()
    blade.read_yaml(str(yaml))
    topo = blade.windio2_topology
    base = generate_mesh(topo, ES, (0.85, 0.87), min_strip_m=1e-6)
    cut = generate_mesh(topo, ES, (0.85, 0.87), min_strip_m=1e-6, collapse_ratio=0.15)
    again = generate_mesh(topo, ES, (0.85, 0.87), min_strip_m=1e-6)
    assert np.array_equal(base.elements, again.elements)  # default path is untouched
    assert aspect(np.asarray(cut.nodes), np.asarray(cut.elements)).max() < aspect(np.asarray(base.nodes), np.asarray(base.elements)).max()
    assert aspect(np.asarray(base.nodes, float), np.asarray(base.elements)).max() > 10  # the TE band is there
    assert aspect(np.asarray(cut.nodes, float), np.asarray(cut.elements)).max() < 10
    assert len(cut.elements) < len(base.elements)
    # a patch without slivers comes back identical
    clean = generate_mesh(topo, ES, (0.30, 0.32), min_strip_m=1e-6)
    same = generate_mesh(topo, ES, (0.30, 0.32), min_strip_m=1e-6, collapse_ratio=0.15)
    assert np.array_equal(clean.elements, same.elements) and np.array_equal(clean.nodes, same.nodes)
    n = len(cut.elements)
    assert len(cut.cells) == len(cut.surfaces) == len(cut.element_parameters) == len(cut.bases) == n
    assert all(0 <= i < len(cut.nodes) for chain in cut.attachments.values() for i in chain)
