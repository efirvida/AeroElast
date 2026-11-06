"""Shear-web meshing of the gmsh skin mesher on a synthetic cylinder with one diametral web (no blade, no PETSc)."""
import numpy as np
import pytest

pytest.importorskip("gmsh")
gs = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")

H, L, PER = 0.25, 3.0, 2.0
OPTS = {"patch_rows": 1, "extra_options": {"Mesh.Smoothing": 10.0, "Mesh.RecombineMinimumQuality": 0.7}}


def cylinder(eta, nd):
    r = PER / (2 * np.pi)
    phi = -2 * np.pi * np.asarray(nd, float)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), np.asarray(eta, float) * L])


def skin_with_web():
    per = lambda e: np.full_like(np.asarray(e, float), PER)  # noqa: E731
    lines = [gs.Line("web/w/start_nd_arc", np.array([0.0, 1.0]), np.array([0.25, 0.25]), 0.0, 1.0, True),
             gs.Line("web/w/end_nd_arc", np.array([0.0, 1.0]), np.array([0.75, 0.75]), 0.0, 1.0, True)]
    src = gs.SkinSource(lines, L, per, [], cylinder)
    res = gs.generate_skin_mesh(src, H, enforce_gate=False, structured=True, h_range=(0.1, 0.25), **OPTS)
    el = res.elements
    used = np.unique(res.master[el[el >= 0]])
    remap = np.full(len(res.master), -1)
    remap[used] = np.arange(len(used))
    skin_el = np.where(el >= 0, remap[res.master[np.where(el >= 0, el, 0)]], -1)
    rn = res.info["line_row_nodes"]
    f = lambda a: np.where(a >= 0, remap[res.master[np.where(a >= 0, a, 0)]], -1)  # noqa: E731
    return res, res.xyz[used], skin_el, f(rn["web/w/start_nd_arc"]), f(rn["web/w/end_nd_arc"]), res.info["rows"]


def web_of(active_rows=None):
    res, xyz, skin_el, ns, ne, rows = skin_with_web()
    active = np.ones(len(rows) - 1, bool) if active_rows is None else active_rows(len(rows) - 1)
    wm = gs.mesh_web_chart(xyz, ns, ne, rows, active, H, lambda e: np.full_like(np.asarray(e, float), H), L, **OPTS)
    return xyz, skin_el, ns, ne, rows, active, wm


def test_web_edges_are_skin_nodes_and_junction_is_conforming():
    xyz, skin_el, ns, ne, rows, active, wm = web_of()
    allxyz = np.vstack([xyz, wm["xyz_new"]])
    wel = wm["elements"]
    edge_ids = set(ns[ns >= 0].tolist()) | set(ne[ne >= 0].tolist())
    used = set(wel[wel >= 0].tolist())
    assert len(ns) == len(rows) and (ns >= 0).all() and (ne >= 0).all()
    assert edge_ids <= used  # every attachment node is a web corner (one division per row segment)
    assert min(used & set(range(len(xyz)))) >= 0 and (used & set(range(len(xyz)))) == edge_ids  # only skin ids on the two edges
    assert not gs.conformity_report(allxyz, np.vstack([skin_el, wel]))["coincident_distinct_id_pairs"]
    # the web edge nodes lie on the skin attachment lines
    assert np.allclose(allxyz[sorted(edge_ids)][:, 0] ** 2 + allxyz[sorted(edge_ids)][:, 1] ** 2, (PER / (2 * np.pi)) ** 2)
    # every skin attachment edge is the edge of exactly one web element (junction closed)
    edges = {}
    for e in wel:
        ids = e[e >= 0]
        for a, b in zip(ids, np.roll(ids, -1), strict=True):
            edges[(min(a, b), max(a, b))] = edges.get((min(a, b), max(a, b)), 0) + 1
    for chain in (ns, ne):
        for a, b in zip(chain[:-1], chain[1:], strict=True):
            assert edges[(min(a, b), max(a, b))] == 1


def test_web_is_in_the_plane_has_positive_area_and_good_quality():
    xyz, _, _, _, _, _, wm = web_of()
    allxyz = np.vstack([xyz, wm["xyz_new"]])
    wel = wm["elements"]
    m = gs.element_metrics(allxyz, wel, np.tile([1.0, 0.0, 0.0], (len(wel), 1)))
    assert np.allclose(allxyz[np.unique(wel[wel >= 0])][:, 0], 0.0, atol=1e-9)  # the diametral web is the plane x = 0
    assert abs(m["area"].sum() - L * 2 * PER / (2 * np.pi)) < 1e-9
    assert m["aspect"].max() <= 3.0 and m["min_angle"].min() >= 30.0
    assert gs.quality_gate({**m, "jacobian": np.abs(m["jacobian"])}, raise_on_fail=False)["passed"]


def test_web_only_on_active_intervals():
    xyz, _, ns, ne, rows, active, wm = web_of(lambda n: np.arange(n) % 12 >= 4)
    wel = wm["elements"]
    allxyz = np.vstack([xyz, wm["xyz_new"]])
    z = allxyz[np.unique(wel[wel >= 0])][:, 2]
    first = np.nonzero(active)[0].min()
    assert z.min() == pytest.approx(rows[first]) and z.max() == pytest.approx(rows[np.nonzero(active)[0].max() + 1])


def test_web_rejects_an_active_interval_without_skin_node():
    res, xyz, _, ns, ne, rows = skin_with_web()
    ns = ns.copy()
    ns[len(rows) // 2] = -1
    ns[len(rows) // 2 + 1] = -1
    active = np.ones(len(rows) - 1, bool)
    with pytest.raises(gs.MeshQualityError):
        gs.mesh_web_chart(xyz, ns, ne, rows, active, H, lambda e: np.full_like(np.asarray(e, float), H), L, **OPTS)


def test_orient_positive_reverses_only_negative_charts_and_keeps_triangles_padded():
    import numpy as np

    from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import MeshQualityError, _orient_positive

    nan = 0.0
    quad_ok = np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]])
    tri_neg = np.array([[0, 0], [0, 1], [1, 0], [nan, nan]])
    params = np.stack([quad_ok, tri_neg])
    elements = np.array([[0, 1, 2, 3], [4, 5, 6, -1]])
    el, par = _orient_positive(elements, params)
    assert el[0].tolist() == [0, 1, 2, 3]
    assert el[1].tolist() == [4, 6, 5, -1]
    assert np.isfinite(par).all() and par[1, 1].tolist() == [1, 0]  # padding repeats a corner, never NaN
    degenerate = np.stack([quad_ok, np.zeros((4, 2))])
    with pytest.raises(MeshQualityError):
        _orient_positive(elements, degenerate)
