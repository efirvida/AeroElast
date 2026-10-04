"""Exact-property tests for :mod:`aeroelast.solvers.bem.section_contour`.

Pure geometry: no external reference, no solver, no ``_aeroelast`` import.
Every assertion is an exact property of the inputs, except the translation
test, which uses ``rel=1e-12`` to absorb float addition only.
"""

import numpy as np
import pytest

from aeroelast.solvers.bem import section_contour as sc

SPAN_Z = np.array([0.0, 0.0, 1.0])
SPAN_Y = np.array([0.0, 1.0, 0.0])

# A convex pentagon in the z = 0 plane, used for the invariance tests.
PENTA = np.array(
    [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 1.0, 0.0], [1.5, 3.0, 0.0], [0.0, 1.0, 0.0]]
)


def _orientation(area: float) -> int:
    return 1 if area > 0.0 else (-1 if area < 0.0 else 0)


def test_shuffled_square_orders_ccw_and_is_usable():
    # Raw corners are CCW; deliver them in a shuffled order.
    square = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    shuffled = square[[2, 0, 3, 1]]

    order = sc.order_ring(shuffled, SPAN_Z)
    ordered = shuffled[order]
    assert sc.is_simple(ordered, SPAN_Z)[0] is True
    assert sc.signed_area(ordered, SPAN_Z) == 1.0

    report = sc.contour_report(shuffled, SPAN_Z)
    assert report.orientation == 1
    assert report.simple is True
    assert report.closed is True
    assert report.usable is True
    assert report.area == 1.0
    assert sorted(report.order.tolist()) == [0, 1, 2, 3]


def test_clockwise_input_is_normalised_but_raw_sequence_keeps_sign():
    # Clockwise square: the ordered report flips it to the +span-CCW convention.
    clockwise = np.array([[0.0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]])

    report = sc.contour_report(clockwise, SPAN_Z)
    assert report.orientation == 1
    assert report.usable is True

    # The raw sequence is still evaluated with its own (negative) sign.
    assert _orientation(sc.signed_area(clockwise, SPAN_Z)) == -1
    assert sc.is_simple(clockwise, SPAN_Z)[0] is True


def test_self_intersecting_and_self_overlapping_rings_are_not_usable():
    # A raw bow-tie order: edges (0,1) and (2,3) cross.
    bowtie = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])
    simple, witness = sc.is_simple(bowtie, SPAN_Z)
    assert simple is False
    assert witness is not None
    assert len(witness) == 2

    # A ring that reordering cannot repair: two vertices share one in-plane
    # projection (they differ only along span), so the ordered ring keeps a
    # zero-length wall edge and must be reported as not simple / not usable.
    degenerate = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    report = sc.contour_report(degenerate, SPAN_Z)
    assert report.usable is False
    assert report.simple is False
    assert "simple" in report.reason


def test_repeated_consecutive_point_is_not_closed():
    points = np.array([[0.0, 0, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    report = sc.contour_report(points, SPAN_Z)
    assert report.closed is False
    assert report.usable is False
    assert "closed" in report.reason


def test_degenerate_inputs_are_not_usable():
    fewer_than_three = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert sc.contour_report(fewer_than_three, SPAN_Z).usable is False

    collinear = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    report = sc.contour_report(collinear, SPAN_Z)
    assert report.usable is False
    assert report.area == 0.0


def test_non_z_span_frame_chirality_is_pinned():
    # Square in the xz plane; span is +y.
    square = np.array([[0.0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]])

    report = sc.contour_report(square, SPAN_Y)
    # Exact, signed value (not abs): pins u x v == +span.
    assert report.area == 1.0
    assert report.orientation == 1
    assert report.usable is True

    u, v = sc.section_plane_axes(SPAN_Y)
    assert u @ u == 1.0
    assert v @ v == 1.0
    assert u @ v == 0.0
    assert np.array_equal(np.cross(u, v), SPAN_Y)


def test_ordering_is_deterministic_under_shuffles():
    first = PENTA[[0, 1, 2, 3, 4]]
    second = PENTA[[3, 4, 0, 2, 1]]

    order_first = sc.order_ring(first, SPAN_Z)
    order_second = sc.order_ring(second, SPAN_Z)
    # Same ordered point sequence, not necessarily the same index values.
    assert np.array_equal(first[order_first], second[order_second])


def test_in_plane_translation_preserves_area_and_perimeter():
    shift = np.array([3.0, 4.0, 0.0])

    base = sc.contour_report(PENTA, SPAN_Z)
    moved = sc.contour_report(PENTA + shift, SPAN_Z)
    assert moved.area == pytest.approx(base.area, rel=1e-12)
    assert moved.perimeter == pytest.approx(base.perimeter, rel=1e-12)


def test_zero_span_is_rejected():
    with pytest.raises(ValueError):
        sc.section_plane_axes(np.zeros(3))


# ---------------------------------------------------------------------------
# Span-station grouping (merge_span_stations / rings_in_band)
#
# The projector's strips contain several *physical* rings: prebend splits one
# ring's nodes over ~1e-3 m of span while neighbouring BEM stations are ~0.44 m
# apart. Plain properties, not store comparisons.
# ---------------------------------------------------------------------------


def test_merge_span_stations_merges_close_values_as_means():
    # Two physical rings, each split by a prebend of ~1e-3 m, 0.44 m apart.
    values = np.array([0.0, 0.0004, 0.0009, 0.44, 0.4403, 0.4408])
    merged = sc.merge_span_stations(values, 0.01)
    assert merged == pytest.approx([0.0013 / 3.0, (0.44 + 0.4403 + 0.4408) / 3.0])


def test_merge_span_stations_is_increasing_and_keeps_exact_duplicates():
    # Zero tolerance is the planar section: all nodes share one span, they must
    # still come back as one station (exact duplicates always merge).
    assert sc.merge_span_stations(np.array([2.0, 2.0, 5.0]), 0.0) == [2.0, 5.0]
    merged = sc.merge_span_stations(np.array([3.0, 1.0, 2.0, 1.0]), 0.01)
    assert merged == [1.0, 2.0, 3.0]


def test_rings_in_band_covers_every_index_once_in_span_order():
    spans = np.array([10.0, 0.0, 0.0005, 20.0, 0.0, 20.0004])
    groups = sc.rings_in_band(spans, 0.01)
    assert len(groups) == 3
    flat = np.concatenate(groups)
    assert sorted(flat.tolist()) == [0, 1, 2, 3, 4, 5]
    means = [float(spans[g].mean()) for g in groups]
    assert means == sorted(means)


def test_rings_in_band_merges_a_split_prebent_ring():
    # One physical ring: same chordwise station, prebend spreads it by ~1e-3 m.
    spans = np.array([5.0, 5.0004, 4.9996, 5.0002])
    groups = sc.rings_in_band(spans, 0.01)
    assert len(groups) == 1
    assert sorted(groups[0].tolist()) == [0, 1, 2, 3]


def test_rings_in_band_keeps_well_separated_rings_separate():
    spans = np.array([0.0, 0.44, 0.88, 0.4401, 0.0002])
    groups = sc.rings_in_band(spans, 0.01)
    assert [sorted(g.tolist()) for g in groups] == [[0, 4], [1, 3], [2]]


# ---------------------------------------------------------------------------
# Cells of a planar wall graph (SectionCell / section_cells / cell_adjacency)
#
# A thin-walled blade section is measurably multi-cell: the outer skin plus
# internal webs bound several closed cells, so the single-cell q = T/(2A) is the
# wrong shear-flow field there.  These tests pin the bounded-face extraction on
# synthetic planar wall graphs whose cell areas are known exactly.  Plain exact
# properties of our own data, not store comparisons.
# ---------------------------------------------------------------------------


def _cell_signature(cells):
    """Ordered ``(boundary, area)`` view with the boundary compared as a set."""
    return [(sorted(c.boundary.tolist()), c.area) for c in cells]


def test_square_ring_is_one_cell():
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

    cells = sc.section_cells(points, edges, SPAN_Z)
    assert len(cells) == 1
    cell = cells[0]
    assert cell.area == 1.0
    # The bounded face is normalised CCW starting at its lowest node index.
    assert cell.boundary.tolist() == [0, 1, 2, 3]
    assert set(cell.edges) == {(0, 1), (1, 2), (2, 3), (3, 0)}
    assert len(cell.edges) == 4


def test_two_cell_rectangle_reports_the_shared_web():
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0 bottom-left
            [2.0, 0.0, 0.0],  # 1 bottom-right
            [2.0, 1.0, 0.0],  # 2 top-right
            [0.0, 1.0, 0.0],  # 3 top-left
            [1.0, 0.0, 0.0],  # 4 bottom-middle
            [1.0, 1.0, 0.0],  # 5 top-middle
        ]
    )
    edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]

    cells = sc.section_cells(points, edges, SPAN_Z)
    assert len(cells) == 2
    assert sorted(c.area for c in cells) == [1.0, 1.0]

    adjacency = sc.cell_adjacency(cells)
    assert sorted(adjacency[(4, 5)]) == [0, 1]  # the web bounds both cells
    skin = [key for key in adjacency if key != (4, 5)]
    assert len(skin) == 6
    assert all(len(adjacency[key]) == 1 for key in skin)


def test_three_cell_rectangle_pins_blade_topology():
    # This is the blade's own topology: outer skin plus two parallel webs.
    points = np.array(
        [
            [0.0, 0.0, 0.0],  # 0 bottom-left
            [3.0, 0.0, 0.0],  # 1 bottom-right
            [3.0, 1.0, 0.0],  # 2 top-right
            [0.0, 1.0, 0.0],  # 3 top-left
            [1.0, 0.0, 0.0],  # 4 bottom of web 1
            [2.0, 0.0, 0.0],  # 5 bottom of web 2
            [1.0, 1.0, 0.0],  # 6 top of web 1
            [2.0, 1.0, 0.0],  # 7 top of web 2
        ]
    )
    edges = [
        (0, 4),
        (4, 5),
        (5, 1),
        (1, 2),
        (2, 7),
        (7, 6),
        (6, 3),
        (3, 0),
        (4, 6),  # web 1
        (5, 7),  # web 2
    ]

    cells = sc.section_cells(points, edges, SPAN_Z)
    assert len(cells) == 3
    assert sorted(c.area for c in cells) == [1.0, 1.0, 1.0]

    adjacency = sc.cell_adjacency(cells)
    assert set(adjacency) == {(min(a, b), max(a, b)) for a, b in edges}
    assert len(adjacency[(4, 6)]) == 2
    assert len(adjacency[(5, 7)]) == 2
    histogram: dict[int, int] = {}
    for owners in adjacency.values():
        histogram[len(owners)] = histogram.get(len(owners), 0) + 1
    assert histogram == {1: 8, 2: 2}


def test_dangling_stub_does_not_create_or_corrupt_a_cell():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, -0.5, 0.0],  # 4: free end of the stub, off any wall ray
        ]
    )
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)]

    cells = sc.section_cells(points, edges, SPAN_Z)
    assert len(cells) == 1
    assert cells[0].area == 1.0
    assert 4 not in cells[0].boundary.tolist()
    assert set(cells[0].boundary.tolist()) == {0, 1, 2, 3}

    adjacency = sc.cell_adjacency(cells)
    assert (0, 4) not in adjacency
    assert (4, 0) not in adjacency


def test_non_z_span_axis_yields_the_same_cells():
    # Same physical square expressed in two planes; the cells and their areas
    # must not depend on which span axis the section frame is built from.
    square_z = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    square_y = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.0, 0.0, 1.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

    cells_z = sc.section_cells(square_z, edges, SPAN_Z)
    cells_y = sc.section_cells(square_y, edges, SPAN_Y)
    assert _cell_signature(cells_y) == _cell_signature(cells_z)
    assert len(cells_y) == 1
    assert cells_y[0].area == 1.0


def test_cells_are_independent_of_edge_order_and_endpoint_order():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ]
    )
    edges = np.array(
        [
            [0, 4],
            [4, 5],
            [5, 1],
            [1, 2],
            [2, 7],
            [7, 6],
            [6, 3],
            [3, 0],
            [4, 6],
            [5, 7],
        ],
        dtype=np.intp,
    )
    base = sc.section_cells(points, edges, SPAN_Z)
    base_signature = [(c.boundary.tolist(), c.area) for c in base]
    assert len(base) == 3

    rng = np.random.default_rng(20240517)
    for _ in range(5):
        shuffled = edges[rng.permutation(edges.shape[0])]
        flip = rng.random(shuffled.shape[0]) < 0.5
        shuffled[flip] = shuffled[flip][:, ::-1]
        cells = sc.section_cells(points, shuffled, SPAN_Z)
        assert [(c.boundary.tolist(), c.area) for c in cells] == base_signature


def test_graphs_without_a_closed_loop_have_no_cells():
    chain = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]])

    assert sc.section_cells(np.zeros((0, 3)), [], SPAN_Z) == []
    assert sc.section_cells(chain[:3], [], SPAN_Z) == []
    assert sc.section_cells(chain[:2], [(0, 1)], SPAN_Z) == []
    assert sc.section_cells(chain, [(0, 1), (1, 2), (2, 3)], SPAN_Z) == []


def test_section_cells_does_not_mutate_its_inputs():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ]
    )
    edges = np.array([[0, 4], [4, 1], [1, 2], [2, 5], [5, 3], [3, 0], [4, 5]], dtype=np.intp)
    points_before = points.copy()
    edges_before = edges.copy()

    cells = sc.section_cells(points, edges, SPAN_Z)
    assert len(cells) == 2
    assert np.array_equal(points, points_before)
    assert np.array_equal(edges, edges_before)
