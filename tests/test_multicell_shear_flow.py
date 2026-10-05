"""Multi-cell Bredt-Batho shear-flow solver: closed-form and independent tests.

:func:`aeroelast.solvers.bem.force_projection.multi_cell_shear_flow` solves the
thin-walled multi-cell torsion system

    cell i:  sum over the walls of cell i of  q_wall * ell_wall / t_wall = 2*A_i*phi
    torque:  sum_i 2*A_i*q_i = T

where a wall bounding only cell i carries ``q_i`` and a wall shared by cells i
and j carries ``q_i - q_j``.  ``phi = G*theta'`` is the common twist stiffness
returned by the solver (units N/m^3 = stress/length; the caller keeps ``G`` and
recovers ``theta' = phi / G``).  ``G`` cancels from the linear system because
every compatibility equation is homogeneous in it.

References used below (all stated in closed form in the test docstrings):

* single cell: ``q = T / (2 A)`` and ``phi = T P / (4 A^2 t)`` (Bredt-Batho);
* symmetric two-cell box: the central web is inert, so the result equals the
  single-cell value for the same outer dimensions and thickness;
* asymmetric two-cell box: an independent 3x3 system assembled here by hand
  from the wall lengths and solved with ``numpy.linalg.solve``.

Cells always come from the production extractor
(:func:`section_contour.section_cells` / :func:`section_contour.cell_adjacency`)
so the solver is fed real ``SectionCell`` objects, never hand-rolled ones
(except the deliberate degenerate cases in the last tests).
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.solvers.bem.force_projection import multi_cell_shear_flow
from aeroelast.solvers.bem.section_contour import SectionCell, cell_adjacency, section_cells

SPAN_Z = np.array([0.0, 0.0, 1.0])


def _wall_graph(points, edges, thickness):
    """Build ``(cells, adjacency, edge_length, edge_thickness)`` from a wall graph.

    ``thickness`` is either a scalar (applied to every edge) or a mapping from the
    undirected edge key ``(min(a, b), max(a, b))`` to a thickness.
    """
    pts = np.asarray(points, dtype=float)
    cells = section_cells(pts, edges, SPAN_Z)
    adjacency = {key: tuple(value) for key, value in cell_adjacency(cells).items()}
    edge_length: dict[tuple[int, int], float] = {}
    edge_thickness: dict[tuple[int, int], float] = {}
    for a, b in edges:
        key = (int(a), int(b)) if a < b else (int(b), int(a))
        edge_length[key] = float(np.linalg.norm(pts[a] - pts[b]))
        if hasattr(thickness, "keys"):
            edge_thickness[key] = float(thickness[key])
        else:
            edge_thickness[key] = float(thickness)
    return cells, adjacency, edge_length, edge_thickness


def _cell_with_node(cells, node):
    return next(i for i, c in enumerate(cells) if node in c.boundary.tolist())


def _compatibility_residual(cells, adjacency, edge_length, edge_thickness, index, flows, phi):
    """``sum_wall q_wall * ell / t - 2 * A_i * phi`` for cell ``index``.

    The wall flow is ``q_i`` on a sole-owner wall and ``q_i - q_j`` on a wall
    shared by cells i and j: the sign convention the multi-cell system uses.
    """
    cell = cells[index]
    total = 0.0
    for a, b in cell.edges:
        key = (a, b) if a < b else (b, a)
        wall_flow = flows[index]
        for other in adjacency[key]:
            if other != index:
                wall_flow -= flows[other]
        total += wall_flow * edge_length[key] / edge_thickness[key]
    return total - 2.0 * cell.area * phi


def test_single_cell_ring_matches_bredt():
    """One closed cell must reduce exactly to Bredt-Batho.

    Closed form: ``q = T / (2 A)`` and ``phi = G theta' = T P / (4 A^2 t)``,
    with ``A`` the enclosed area and ``P`` the perimeter of the same ring.
    """
    points = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    thickness, torsion = 0.01, 5.0

    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, thickness)
    assert len(cells) == 1
    area = cells[0].area
    perimeter = sum(edge_length.values())

    flows, phi = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    assert flows[0] == pytest.approx(torsion / (2.0 * area), rel=1e-12)
    assert phi == pytest.approx(torsion * perimeter / (4.0 * area**2 * thickness), rel=1e-12)


def test_symmetric_two_cell_box_web_carries_no_net_flow():
    """Equal cells sharing one web: the web is inert.

    Bredt-Batho on the symmetric two-cell box gives equal cell flows, so the
    shared web carries ``q_1 - q_2 = 0`` and ``phi`` collapses to the
    single-cell value ``T P_outer / (4 A_outer^2 t)`` with ``A_outer`` the whole
    enclosed rectangle and ``P_outer`` its outer perimeter (this is the case
    ``tests/test_box_multicell_torsion.py`` documents).
    """
    width, height, thickness, torsion = 2.0, 1.0, 0.01, 5.0
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [width, 0.0, 0.0],
            [width, height, 0.0],
            [0.0, height, 0.0],
            [width / 2.0, 0.0, 0.0],
            [width / 2.0, height, 0.0],
        ]
    )
    edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, thickness)
    assert len(cells) == 2

    flows, phi = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    assert flows[0] == pytest.approx(flows[1], rel=1e-12)
    left, right = adjacency[(4, 5)]
    assert flows[left] - flows[right] == pytest.approx(0.0, abs=1e-12)

    outer_area = width * height
    outer_perimeter = 2.0 * (width + height)
    assert phi == pytest.approx(
        torsion * outer_perimeter / (4.0 * outer_area**2 * thickness), rel=1e-12
    )


def test_asymmetric_two_cell_box_matches_independent_hand_assembly():
    """Independent 3x3 reference for the multi-cell system.

    Geometry: a 3 x 2 box with a web at x = 1, uniform thickness t; cell areas
    A0 = 2 (x in [0, 1]) and A1 = 4 (x in [1, 3]).  The reference assembles

        [ D0     -d_web   -2 A0 ] [ q0 ]   [ 0 ]
        [ -d_web  D1      -2 A1 ] [ q1 ] = [ 0 ]
        [ 2 A0   2 A1      0    ] [ phi]   [ T ]

    from the wall lengths alone, with D0 = 6/t, D1 = 8/t and d_web = 2/t, and
    solves it with ``numpy.linalg.solve``.
    """
    thickness, torsion = 0.01, 7.0
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
        ]
    )
    edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, thickness)
    assert len(cells) == 2
    area0, area1 = cells[0].area, cells[1].area
    assert area0 == 2.0 and area1 == 4.0

    flows, phi = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    # Independent reference: D0 = 6/t, D1 = 8/t, web = 2/t (from the lengths).
    reference = np.array(
        [
            [6.0 / thickness, -2.0 / thickness, -2.0 * area0],
            [-2.0 / thickness, 8.0 / thickness, -2.0 * area1],
            [2.0 * area0, 2.0 * area1, 0.0],
        ]
    )
    reference_solution = np.linalg.solve(reference, np.array([0.0, 0.0, torsion]))
    reference_q = reference_solution[:2]
    reference_phi = reference_solution[2]

    assert flows[0] == pytest.approx(reference_q[0], rel=1e-10)
    assert flows[1] == pytest.approx(reference_q[1], rel=1e-10)
    assert phi == pytest.approx(reference_phi, rel=1e-10)

    # Exact invariants, independent of the reference matrix.
    assert 2.0 * (area0 * flows[0] + area1 * flows[1]) == pytest.approx(torsion, rel=1e-10)
    for index in (0, 1):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_thickness, index, flows, phi
        ) == pytest.approx(0.0, abs=1e-8)
    left, right = adjacency[(4, 5)]  # the shared web, owners in cell order
    assert left != right


def test_three_cell_blade_topology_shares_one_phi():
    """Blade topology: outer skin plus two parallel webs, three equal cells.

    For the parallel-web rectangle the reflection symmetry only forces the two
    outer cells to be equal (``q_left = q_right``); the middle cell is coupled to
    a neighbour by *two* webs while each outer cell has one web plus a free side
    wall, so ``q_middle`` differs.  The physical invariants that must hold
    exactly are the torque equilibrium and one common ``phi`` per cell.
    """
    thickness, torsion = 0.01, 9.0
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
    edges = [
        (0, 4),
        (4, 5),
        (5, 1),
        (1, 2),
        (2, 7),
        (7, 6),
        (6, 3),
        (3, 0),
        (4, 6),
        (5, 7),
    ]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, thickness)
    assert len(cells) == 3
    assert [c.area for c in cells] == [1.0, 1.0, 1.0]

    flows, phi = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    torque = 2.0 * sum(cell.area * flow for cell, flow in zip(cells, flows, strict=True))
    assert torque == pytest.approx(torsion, rel=1e-12)
    for index in range(3):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_thickness, index, flows, phi
        ) == pytest.approx(0.0, abs=1e-9)

    left = _cell_with_node(cells, 0)
    right = _cell_with_node(cells, 1)
    assert flows[left] == pytest.approx(flows[right], rel=1e-12)


def test_three_cell_fan_is_symmetric_and_webs_carry_no_net_flow():
    """A 3-fold symmetric fan has equal cell flows, hence inert internal webs.

    Three congruent cells around a central node (angles 90/210/330 deg, equipped
    with a 3-fold rotation symmetry) must all carry the same flow, so every
    internal wall's net flow ``q_i - q_j`` is zero.  This is the configuration
    in which the "middle cell equals its neighbours" property actually holds.
    """
    radius, thickness, torsion = 1.0, 0.01, 9.0
    angles = np.deg2rad([90.0, 210.0, 330.0])
    points = np.zeros((4, 3))
    points[1:, 0] = radius * np.cos(angles)
    points[1:, 1] = radius * np.sin(angles)
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, thickness)
    assert len(cells) == 3
    assert cells[0].area == pytest.approx(cells[1].area, rel=1e-12)
    assert cells[1].area == pytest.approx(cells[2].area, rel=1e-12)

    flows, phi = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    assert flows[0] == pytest.approx(flows[1], rel=1e-12)
    assert flows[1] == pytest.approx(flows[2], rel=1e-12)
    for key in [(0, 1), (0, 2), (0, 3)]:
        i, j = adjacency[key]
        assert flows[i] - flows[j] == pytest.approx(0.0, abs=1e-12)
    for index in range(3):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_thickness, index, flows, phi
        ) == pytest.approx(0.0, abs=1e-9)


def test_thicker_walled_cell_carries_the_larger_share():
    """Two equal-area cells; the thicker-walled one must carry more flow.

    For equal areas and a common ``phi``, a cell with thicker walls has the
    smaller compliance ``sum ell / t``, so it needs the larger flow to reach the
    same ``2 A phi``.  Only the direction is pinned, no magnitude.
    """
    width, height, thickness, torsion = 2.0, 1.0, 0.01, 5.0
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [width, 0.0, 0.0],
            [width, height, 0.0],
            [0.0, height, 0.0],
            [width / 2.0, 0.0, 0.0],
            [width / 2.0, height, 0.0],
        ]
    )
    edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]

    thin = thickness
    thick = 2.0 * thickness

    def key(a, b):
        return (a, b) if a < b else (b, a)

    by_edge = {key(a, b): thin for a, b in edges}
    # Right cell (bounded by nodes 1, 2, 5, 4) gets the thicker walls; the web
    # stays thin so the two cells differ only in their outer skin thickness.
    by_edge[key(4, 1)] = thick
    by_edge[key(1, 2)] = thick
    by_edge[key(2, 5)] = thick

    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, by_edge)
    assert len(cells) == 2

    flows, _ = multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, torsion)

    thin_cell = _cell_with_node(cells, 0)  # left cell: thin outer skin
    thick_cell = _cell_with_node(cells, 1)  # right cell: thick outer skin
    assert flows[thin_cell] != pytest.approx(flows[thick_cell])
    assert flows[thick_cell] > flows[thin_cell]


def test_degenerate_inputs_raise_value_error():
    """The solver rejects empty, non-positive, missing and singular inputs."""
    rectangle = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(rectangle, edges, 0.01)

    with pytest.raises(ValueError):
        multi_cell_shear_flow([], {}, {}, {}, 1.0)

    zero_area = SectionCell(
        boundary=np.array([0, 1, 2]),
        area=0.0,
        edges=((0, 1), (1, 2), (2, 0)),
    )
    zero_area_edges = {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0}
    with pytest.raises(ValueError):
        multi_cell_shear_flow([zero_area], {}, zero_area_edges, zero_area_edges, 1.0)

    zero_thickness = dict(edge_thickness)
    zero_thickness[(0, 1)] = 0.0
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, zero_thickness, 1.0)

    missing_length = {key: value for key, value in edge_length.items() if key != (0, 1)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, missing_length, edge_thickness, 1.0)

    missing_thickness = {key: value for key, value in edge_thickness.items() if key != (0, 1)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, missing_thickness, 1.0)

    # Singular system: two positive-area cells with no walls.
    wall_less = [
        SectionCell(boundary=np.array([0]), area=1.0, edges=()),
        SectionCell(boundary=np.array([1]), area=1.0, edges=()),
    ]
    with pytest.raises(ValueError):
        multi_cell_shear_flow(wall_less, {}, {}, {}, 1.0)


def test_non_manifold_shared_edge_is_rejected():
    """A wall owned by three cells is not a thin-walled section wall."""
    cell = SectionCell(boundary=np.array([0, 1, 2]), area=1.0, edges=((0, 1), (1, 2), (2, 0)))
    cells = [cell, cell, cell]
    edge_length = {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0}
    edge_thickness = {(0, 1): 0.01, (1, 2): 0.01, (0, 2): 0.01}
    adjacency = {(0, 1): (0, 1, 2), (1, 2): (0, 1, 2), (0, 2): (0, 1, 2)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, 1.0)


def test_inputs_are_not_mutated():
    """The solver is pure: cells, adjacency, lengths and thicknesses survive."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 2.0, 0.0],
        ]
    )
    edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]
    cells, adjacency, edge_length, edge_thickness = _wall_graph(points, edges, 0.01)

    adjacency_before = {key: list(value) for key, value in adjacency.items()}
    length_before = dict(edge_length)
    thickness_before = dict(edge_thickness)
    boundaries_before = [cell.boundary.copy() for cell in cells]

    multi_cell_shear_flow(cells, adjacency, edge_length, edge_thickness, 7.0)

    assert {key: list(value) for key, value in adjacency.items()} == adjacency_before
    assert edge_length == length_before
    assert edge_thickness == thickness_before
    for cell, boundary in zip(cells, boundaries_before, strict=True):
        assert np.array_equal(cell.boundary, boundary)
