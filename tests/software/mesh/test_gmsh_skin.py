"""Cheap tests of the gmsh skin prototype on a synthetic ring (cylinder): no blade, no PETSc."""
import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
gs = pytest.importorskip("aeroelast.models.blade.numad.mesh_gen.gmsh_skin")

H, L, PER = 0.25, 3.0, 2.0  # element size, span length, constant perimeter


def cylinder(eta, nd):
    r = PER / (2 * np.pi)
    phi = -2 * np.pi * np.asarray(nd, float)  # clockwise nd -> outward normal for positive (eta, nd)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), np.asarray(eta, float) * L])


def source(lines, drops=()):
    return gs.SkinSource(lines, L, lambda e: np.full_like(np.asarray(e, float), PER), list(drops), cylinder)


def straight(name, nd, fixed=False, lo=0.0, hi=1.0):
    return gs.Line(name, np.array([0.0, 1.0]), np.array([nd, nd]), lo, hi, fixed)


def mesh_of(lines, **kw):
    return gs.generate_skin_mesh(source(lines), H, enforce_gate=False, **kw)


def line_edges_cover(res, x_target, tol=1e-6):
    """True when mesh edges lying on the chart line x = x_target chain from y = 0 to y = L."""
    chart = res.chart
    on = np.abs(chart[:, 0] - x_target) < tol
    total = 0.0
    for a, b in gs.unique_edges(res.elements):
        if on[a] and on[b]:
            total += abs(chart[a, 1] - chart[b, 1])
    return abs(total - L) < 1e-6


def test_event_lines_are_conforming_mesh_edges():
    res = mesh_of([straight("a", 0.25), straight("b", 0.5), straight("c", 0.75)])
    for nd in (0.25, 0.5, 0.75):
        assert line_edges_cover(res, nd * PER), nd
    assert res.info["clean"]["n_bands_merged"] == 0
    assert res.metrics["jacobian"].min() > 0  # every element has the outward orientation


def test_near_coincident_lines_are_merged_and_reported():
    gap = 0.01  # chart metres, far below tau * h = 0.0625
    res = mesh_of([straight("a", 0.5), straight("b", 0.5 + gap / PER)])
    rep = res.info["clean"]
    assert rep["n_bands_merged"] == 1
    assert rep["total_band_chart_area_m2"] == pytest.approx(gap * L, rel=1e-6)
    assert rep["max_material_boundary_displacement_m"] == pytest.approx(gap / 2, rel=1e-6)
    # exactly one line survives: nothing in the open band between the two original positions
    x = res.chart[:, 0]
    mid = 0.5 * PER + gap / 2
    assert not np.any((x > 0.5 * PER + 1e-6) & (x < 0.5 * PER + gap - 1e-6) & (np.abs(x - mid) > 1e-6))
    assert line_edges_cover(res, 0.5 * PER + gap / 2)
    snap = mesh_of([straight("a", 0.5), straight("b", 0.5 + gap / PER)], mode="snap").info["clean"]
    assert snap["max_material_boundary_displacement_m"] == pytest.approx(gap, rel=1e-6)


def test_wedge_is_closed_where_the_gap_reaches_tau_h():
    # gap grows from 0 to 0.2 m over the span: closes at s = tau*h / 0.2 * L
    b = gs.Line("b", np.array([0.0, 1.0]), np.array([0.5, 0.5 + 0.2 / PER]))
    src = source([straight("a", 0.5), b])
    cleaned = gs.clean_source(src, H, tau=0.25)
    rep = cleaned.report
    assert rep["n_bands_merged"] == 1
    band = rep["bands"][0]
    assert band["s_from_m"] == pytest.approx(0.0, abs=1e-9)
    assert band["s_to_m"] == pytest.approx(0.25 * H / 0.2 * L, rel=1e-6)
    assert band["max_gap_m"] <= 0.25 * H + 1e-9
    assert len(cleaned.polylines["b"]) > 2  # the closing step is part of the line


def test_periodic_seam_has_identical_nodes():
    res = mesh_of([straight("a", 0.5)])
    slave, _master = res.master != np.arange(len(res.master)), None
    pairs = np.column_stack([np.nonzero(slave)[0], res.master[slave]])
    assert len(pairs) > 0
    assert np.allclose(res.chart[pairs[:, 0], 1], res.chart[pairs[:, 1], 1], atol=1e-9)
    assert np.allclose(res.chart[pairs[:, 0], 0], PER, atol=1e-6) and np.allclose(res.chart[pairs[:, 1], 0], 0.0)
    assert np.allclose(res.xyz[pairs[:, 0]], res.xyz[pairs[:, 1]])
    # identified nodes close the ring: unique positions == merged node count
    merged = np.unique(res.master)
    assert len(np.unique(np.round(res.xyz[merged], 9), axis=0)) == len(merged)


def test_deterministic():
    a = mesh_of([straight("a", 0.3), straight("b", 0.7)])
    b = mesh_of([straight("a", 0.3), straight("b", 0.7)])
    assert np.array_equal(a.elements, b.elements) and np.array_equal(a.chart, b.chart)


def test_gate_rejects_and_never_repairs():
    res = mesh_of([straight("a", 0.5)])
    assert res.gate["passed"], res.gate
    xyz = res.xyz.copy()
    el = res.elements[~res.metrics["triangle"]][:1]
    xyz[el[0, 1]] = xyz[el[0, 3]]  # collapse a corner onto the opposite one: degenerate quad
    bad = gs.element_metrics(xyz, el)
    with pytest.raises(gs.MeshQualityError):
        gs.quality_gate(bad, {"max_ar": 3.0})
    # raising leaves the metrics untouched (no repair)
    assert gs.quality_gate(bad, {"max_ar": 3.0}, raise_on_fail=False)["passed"] is False
    with pytest.raises(ValueError):
        gs.generate_skin_mesh(source([straight("a", 0.5)]), 0.5)


def kinked(name, nd, eta_k, d=0.0):
    """A line with a 0.5-nd kink vertex at eta_k (an event station)."""
    return gs.Line(name, np.array([0.0, eta_k, 1.0]), np.array([nd, nd + d, nd]))


def test_snap_coalesces_close_stations_and_default_is_unchanged():
    # two kinks 4 mm apart in span: without snap they are two stations (sliver band), with snap one
    lines = [kinked("a", 0.3, 0.5), kinked("b", 0.7, 0.5 + 0.004 / L)]
    base = gs.clean_source(source(lines), H)
    snapped = gs.clean_source(source(lines), H, snap=0.4)
    assert base.report["n_stations_snapped"] == 0
    assert snapped.report["n_stations_snapped"] == 1
    assert snapped.report["n_stations"] == base.report["n_stations"] - 1
    assert snapped.stations[0] == 0.0 and snapped.stations[-1] == 1.0  # the end stations are kept
    # snap = 0 reproduces the previous stations exactly
    assert np.array_equal(gs.clean_source(source(lines), H, snap=0.0).stations, base.stations)


def test_parity_full_quad_runs_and_has_no_triangles():
    with pytest.raises(Exception, match="divided by 2"):
        mesh_of([straight("a", 0.5)], recombination=3)  # odd divisions: the documented failure
    res = mesh_of([straight("a", 0.5)], recombination=3, parity=1, mesh_scale=1.4)
    assert int(res.metrics["triangle"].sum()) == 0
    assert res.metrics["jacobian"].min() > 0
    assert line_edges_cover(res, 0.5 * PER)


def test_post_optimisation_is_deterministic_and_keeps_conformity():
    a = mesh_of([straight("a", 0.3)], optimize=("Laplace2D", "Relocate2D"), optimize_iter=2)
    b = mesh_of([straight("a", 0.3)], optimize=("Laplace2D", "Relocate2D"), optimize_iter=2)
    assert np.array_equal(a.chart, b.chart) and np.array_equal(a.elements, b.elements)
    assert line_edges_cover(a, 0.3 * PER) and a.metrics["jacobian"].min() > 0


class _Layer:
    def __init__(self, name, t):
        self.name, self.thickness_m, self.material = name, t, {"rho": [1000.0]}


class _FakeLaminate:
    def cells(self, eta, surface, pos, bounds, xyz):
        # stack name by the cell's lower bound, thickness linear in eta
        return [type("C", (), {"layers": [_Layer(f"s{lo:.2f}", 0.01 * (1 + eta))]})() for lo, hi in bounds]


class _FakeTopology:
    _laminate = _FakeLaminate()

    def _positions(self, eta, surface):
        return {"x": 0.4, "y": 0.8}


def test_areal_mass_buckets_match_exact_within_slope_bound():
    rng = np.random.default_rng(0)
    eta, nd = np.sort(rng.uniform(0, 1, 300)), rng.uniform(0, 1, 300)
    top = _FakeTopology()
    exact, k_exact = gs.areal_mass(top, eta, nd)
    fast, k_fast = gs.areal_mass(top, eta, nd, L, 0.05)
    assert k_exact == k_fast  # the cell structure does not depend on eta here
    # areal mass = 10 (1 + eta) kg/m2 is linear in eta: error <= slope * dy (a point is within dy of its bucket mean)
    assert np.max(np.abs(fast - exact)) <= 10.0 * 0.05 / L + 1e-12
    with pytest.raises(ValueError):
        gs.areal_mass(top, eta, nd, None, 0.05)


def structured_of(lines, **kw):
    return gs.generate_skin_mesh(source(lines), H, enforce_gate=False, structured=True, patch_rows=1, **kw)


def test_structured_parallel_strips_are_pure_regular_quads_on_the_rows():
    res = structured_of([straight("a", 0.25), straight("b", 0.5), straight("c", 0.75)])
    assert int(res.metrics["triangle"].sum()) == 0
    assert res.info["n_patches"] == 0 and res.info["n_patch_cells"] == 0
    assert res.metrics["jacobian"].min() > 0 and res.gate["passed"], res.gate
    for nd in (0.25, 0.5, 0.75):
        assert line_edges_cover(res, nd * PER), nd
    # strips are PER/4 = 0.5 m wide: 2 elements of 0.25 m; rows every 0.25 m (L = 3 -> 12 intervals): 4 strips x 2 x 12
    assert len(res.elements) == 4 * 2 * 12
    assert np.allclose(np.unique(np.round(res.chart[:, 1], 9)), np.linspace(0, L, 13))  # every node is on a row
    assert np.max(res.metrics["aspect"]) < 1.05 and np.allclose(res.metrics["max_angle"], 90.0)  # chords of the cylinder


def test_structured_seam_has_identical_nodes_and_closes_the_ring():
    res = structured_of([straight("a", 0.5)])
    slave = res.master != np.arange(len(res.master))
    pairs = np.column_stack([np.nonzero(slave)[0], res.master[slave]])
    assert len(pairs) == 12 + 1  # one node per row
    assert np.allclose(res.chart[pairs[:, 0], 1], res.chart[pairs[:, 1], 1], atol=1e-12)
    assert np.allclose(res.xyz[pairs[:, 0]], res.xyz[pairs[:, 1]])


def test_structured_count_change_makes_a_small_patch_and_keeps_conformity():
    # strip between a (fixed at 0.25) and b: 0.5 m wide at the root, 1.5 m at the tip -> the count goes 2 -> 6
    b = gs.Line("b", np.array([0.0, 1.0]), np.array([0.5, 0.5 + 1.0 / PER]))
    res = structured_of([straight("a", 0.25), b])
    info = res.info
    assert info["n_count_changes"] >= 3 and info["n_patches"] >= 3
    assert info["patch_cell_fraction"] < 0.5
    assert bool(np.all(info["in_patch"][res.metrics["triangle"]]))  # triangles only inside transition patches
    assert res.metrics["jacobian"].min() > 0
    assert line_edges_cover(res, 0.25 * PER)
    # b is slanted: every mesh node on it must lie on the straight line x = 0.5 PER + y / L
    x, y = res.chart[:, 0], res.chart[:, 1]
    on_b = np.abs(x - (0.5 * PER + y / L)) < 1e-9
    assert on_b.sum() >= 13  # one node per row at least
    assert res.metrics["max_angle"][~res.metrics["triangle"]].max() < 150.0  # harsh taper (2 -> 6 elements in 3 m)


def test_structured_thin_band_is_closed_by_the_cleaning_and_reported():
    gap = 0.01
    res = structured_of([straight("a", 0.5), straight("b", 0.5 + gap / PER)], snap=0.4)
    assert res.info["clean"]["n_bands_merged"] == 1
    assert res.info["clean"]["max_material_boundary_displacement_m"] == pytest.approx(gap / 2, rel=1e-6)
    assert res.metrics["jacobian"].min() > 0 and int(res.metrics["triangle"].sum()) == 0
    assert line_edges_cover(res, 0.5 * PER + gap / 2)


def test_web_chart_shares_rows_has_distinct_edges_and_only_active_intervals():
    # a web chart: q in [0, 1] of a width that shrinks 1.2 -> 0.6 m; rows shared with another chart; only rows 2..9 active
    rows = np.linspace(0.0, 3.0, 13)
    width = np.linspace(1.2, 0.6, 13)
    eta = rows / L
    src = gs.SkinSource([], L, lambda e: np.interp(np.asarray(e, float), eta, width), [], None, 0.0, 1.0)
    active = np.zeros(12, bool)
    active[2:10] = True
    cm = gs.mesh_chart_structured(src, gs.Cleaned(eta, {}, [], {}), H, rows=rows, active=active, periodic=False, patch_rows=1)
    y = cm.chart[:, 1]
    assert np.array_equal(cm.master, np.arange(len(cm.master))) and len(cm.seam_pairs) == 0  # q = 0 and q = 1 stay distinct
    assert y.min() == pytest.approx(rows[2]) and y.max() == pytest.approx(rows[10])  # nothing outside the active rows
    left, right = np.isclose(cm.chart[:, 0], 0.0), np.isclose(cm.chart[:, 0], np.interp(y / L, eta, width))
    # edge nodes only on rows, exactly one per active row boundary (1 division per row segment)
    assert np.allclose(np.sort(y[left]), rows[2:11]) and np.allclose(np.sort(y[right]), rows[2:11])
    rn = cm.info["line_row_nodes"]
    assert set(rn) == {"seam0", "seam1"} and (rn["seam0"][2:11] >= 0).all() and (rn["seam0"][:2] < 0).all()
    assert np.allclose(cm.chart[rn["seam1"][2:11], 0], width[2:11])
    # the count follows the width: 4 -> 2 or 3 elements, so transition patches exist
    assert cm.info["n_count_changes"] >= 1 and cm.info["n_patches"] >= 1


def test_lookup_cells_matches_one_call_and_web_surface_uses_its_own_positions():
    top = _FakeTopology()
    eta = np.array([0.1, 0.1, 0.8])
    q = np.array([0.1, 0.6, 0.9])
    cells = gs.lookup_cells(top, "skin", eta, q)
    assert [c.layers[0].name for c in cells] == ["s0.00", "s0.40", "s0.80"]  # cells bounded by the positions 0.4 / 0.8
    web = gs.lookup_cells(top, "webX", eta, q, L, 0.5)  # bucketed: same stack keys, one call per bucket
    assert [c.layers[0].name for c in web] == [c.layers[0].name for c in cells]
    with pytest.raises(ValueError):
        gs.lookup_cells(top, "skin", eta, q, None, 0.5)


def test_conformity_report_counts_coincident_distinct_ids():
    nodes = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 0]], float)  # id 4 duplicates id 1
    el = np.array([[0, 1, 2, 3], [1, 4, 2, -1]])
    rep = gs.conformity_report(nodes, el)
    assert rep["coincident_distinct_id_pairs"] == 1 and rep["unused_nodes"] == 0 and rep["elements_with_repeated_node"] == 0
    assert gs.conformity_report(nodes[:4], el[:1])["coincident_distinct_id_pairs"] == 0


def test_le_landmark_line_from_loft_stations():
    class Loft:
        rows = np.array([0.0, 0.5, 1.0])
        stations = [(None, -0.5 * 4.0, 4.0), (None, -0.49 * 4.0, 4.0), (None, -0.48 * 2.0, 2.0)]  # a0 = -nd_arc_fraction * perimeter
    line = gs.le_landmark_line(Loft)
    assert line.fixed and line.name == "le/landmark"
    assert np.allclose(line.values, [0.5, 0.51, 0.52]) and line.lo == 0.0 and line.hi == 1.0


def test_count_rule_proportional_divides_a_strip_in_equal_parts_and_never_adds_changes():
    b = gs.Line("b", np.array([0.0, 1.0]), np.array([0.5, 0.5 + 1.0 / PER]))  # strip 0.5 -> 1.5 m wide
    both = {}
    for rule in ("round", "proportional"):
        res = structured_of([straight("a", 0.25), b], count_rule=rule)
        both[rule] = res
        assert res.metrics["jacobian"].min() > 0
        # the strip's cells in the block next to the line a have EQUAL chordwise widths on every row they span
        x, y = res.chart[:, 0], res.chart[:, 1]
        for yr in np.unique(np.round(y, 9))[:3]:
            row = np.sort(x[(np.abs(y - yr) < 1e-9) & (x > 0.25 * PER + 1e-9) & (x < 0.5 * PER + yr / L - 1e-9)])
            d = np.diff(np.concatenate([[0.25 * PER], row, [0.5 * PER + yr / L]]))
            assert d.max() / d.min() < 1.0 + 1e-6, (rule, yr, d)
    assert both["proportional"].info["n_count_changes"] <= both["round"].info["n_count_changes"]
    assert 0.8 <= np.mean(both["proportional"].metrics["aspect"] < 3) <= 1.0


def test_row_blocks_keep_every_regular_node_exactly_on_a_row_between_slanted_lines():
    # two slanted lines with different slopes (a multi-row Coons block would tilt its interior grid lines)
    a = gs.Line("a", np.array([0.0, 1.0]), np.array([0.2, 0.3]))
    b = gs.Line("b", np.array([0.0, 1.0]), np.array([0.5, 0.9]))
    on_rows = {}
    for rb in (True,):
        res = structured_of([a, b], row_blocks=rb)
        rows = res.info["rows"]
        patch = res.info["in_patch"]
        reg = np.unique(res.elements[~patch][res.elements[~patch] >= 0])
        pat = set(np.unique(res.elements[patch][res.elements[patch] >= 0]).tolist())
        reg = np.array([i for i in reg if i not in pat])
        off = np.min(np.abs(res.chart[reg, 1][:, None] - rows[None, :]), axis=1)
        on_rows[rb] = float(off.max())
        assert res.metrics["jacobian"].min() > 0
    assert on_rows[True] < 1e-9  # (on the blade a multi-row block tilted its grid lines by up to 0.17 m: row_blocks=False)


def _cell(nd, *layers):
    from aeroelast.models.blade.numad.objects.windio2_laminate import LaminateCell, LayerSample

    zero = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    return LaminateCell("skin", nd, zero, tuple(LayerSample(n, n, m, t, 0.0, None) for n, m, t in layers))


GELCOAT = ("UV_protection", {"orth": 0, "E": 3.44e9, "rho": 1235.0}, 0.001)


def _triax(name):
    return (name, {"orth": 1, "E": [28.2e9, 16.2e9, 15.2e9], "rho": 1940.0}, 0.003)


def _strip(nd_edges):
    """One chordwise strip of quads between vertical edges at chart nd ``nd_edges`` (bottom node i, top node i + m)."""
    m = len(nd_edges)
    elements = np.array([[i, i + 1, i + 1 + m, i + m] for i in range(m - 1)])
    nd = np.asarray(nd_edges, float)
    return elements, np.column_stack([nd[:-1], nd[1:], nd[1:], nd[:-1]])


def test_te_coating_fill_gives_a_te_coating_band_the_laminate_of_its_own_side():
    """T5b L34-L36 / T5c L40 option B: coating-only elements edge-connected to the TE seam (a corner at nd 0 or 1) take the
    cell of an edge neighbour with a laminate, each from its own side; the band is found geometrically, not from the cell
    bounds, which other layer breakpoints split (IEA-15: (0.99, 0.9987) / (0.9987, 1)).  A coating-only element away from the
    TE and the laminate cells keep theirs."""
    # e0 interior coating | e1 laminate A | e2, e3 TE coating up to the seam nd 1 | e4 TE coating from nd 0 | e5 laminate B
    elements, corner_nd = _strip([0.40, 0.50, 0.99, 0.9987, 1.0, 0.01, 0.50])
    corner_nd[4] = [0.0, 0.01, 0.01, 0.0]  # the element after the seam starts again at nd 0 (its left edge is the seam)
    a, b = _cell((0.5, 0.99), GELCOAT, _triax("A")), _cell((0.01, 0.5), GELCOAT, _triax("B"))
    interior = _cell((0.4, 0.5), GELCOAT)
    te = [_cell((0.99, 0.9987), GELCOAT), _cell((0.9987, 1.0), GELCOAT), _cell((0.0, 0.01), GELCOAT)]
    cells, n = gs.fill_te_coating_cells(elements, [interior, a, *te, b], corner_nd)
    assert n == 3
    assert cells[2] is a and cells[3] is a and cells[4] is b
    assert cells[0] is interior and cells[1] is a and cells[5] is b


def test_te_coating_fill_leaves_an_all_coating_ring_and_a_laminate_only_mesh_alone():
    elements, corner_nd = _strip([0.0, 0.5, 1.0])
    only = [_cell((0.0, 1.0), GELCOAT)] * 2
    cells, n = gs.fill_te_coating_cells(elements, only, corner_nd)
    assert n == 0 and all(c is only[0] for c in cells)
    lam = [_cell((0.0, 0.5), GELCOAT, _triax("A")), _cell((0.5, 1.0), GELCOAT, _triax("B"))]
    cells, n = gs.fill_te_coating_cells(elements, lam, corner_nd)
    assert n == 0 and cells == lam


def test_te_coating_fill_finds_a_band_touching_only_one_seam_side():
    elements, corner_nd = _strip([0.50, 0.99, 0.9987, 1.0])  # band reaches nd 1 only
    a = _cell((0.5, 0.99), GELCOAT, _triax("A"))
    cells, n = gs.fill_te_coating_cells(elements, [a, _cell((0.99, 0.9987), GELCOAT), _cell((0.9987, 1.0), GELCOAT)], corner_nd)
    assert n == 2 and cells[1] is a and cells[2] is a
    elements, corner_nd = _strip([0.0, 0.0013, 0.01, 0.5])  # band starts at nd 0 only
    b = _cell((0.01, 0.5), GELCOAT, _triax("B"))
    cells, n = gs.fill_te_coating_cells(elements, [_cell((0.0, 0.0013), GELCOAT), _cell((0.0013, 0.01), GELCOAT), b], corner_nd)
    assert n == 2 and cells[0] is b and cells[1] is b
