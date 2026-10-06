"""Multi-cell Bredt-Batho shear-flow solver: closed-form and independent tests.

:func:`aeroelast.solvers.bem.force_projection.multi_cell_shear_flow` solves the
thin-walled multi-cell torsion system

    cell i:  sum over the walls of cell i of  q_wall * ell_wall / S_wall = 2*A_i*theta'
    torque:  sum_i 2*A_i*q_i = T

where a wall bounding only cell i carries ``q_i`` and a wall shared by cells i
and j carries ``q_i - q_j``.  ``S = G*t`` is the wall's membrane shear stiffness
in N/m and ``theta'`` is the twist rate in rad/m returned by the solver.  Using
``S`` per wall (rather than one ``t`` for the whole section) is what lets
spar-cap, skin and web laminates split the flow; for a uniform ``S`` the system
reduces exactly to the old thickness-only one, with ``theta' = phi / G``.

References used below (all stated in closed form in the test docstrings):

* single cell: ``q = T / (2 A)`` and ``theta' = T P / (4 A^2 S)`` with ``S = G*t``
  (Bredt-Batho);
* symmetric two-cell box: the central web is inert, so the result equals the
  single-cell value for the same outer dimensions and thickness;
* asymmetric two-cell box: an independent 3x3 system assembled here by hand
  from the wall lengths and solved with ``numpy.linalg.solve``.

Cells always come from the production extractor
(:func:`section_contour.section_cells` / :func:`section_contour.cell_adjacency`)
so the solver is fed real ``SectionCell`` objects, never hand-rolled ones
(except the deliberate degenerate cases in the last tests).

The file also carries the sign-chain guard for the applied pitching moment (issue #11): a
cheap synthetic-section test that ``_strip_moment_axis_sign`` follows ``_section_ends`` and
flips with the chord axis, and one slow test that the minimum-norm projector and the
with-properties multi-cell projector - two faithful realisations of the SAME rated BEM
moment - rotate the real IEA-15MW blade tip the same way.  The slow test needs the real
blade mesh, the AeroDyn polars and **two** shell solves, so it is marked ``slow``.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.solvers.bem.engine import BEMResult
from aeroelast.solvers.bem.force_projection import (
    ForceProjector,
    multi_cell_shear_flow,
    ring_section,
)
from aeroelast.solvers.bem.section_contour import (
    SectionCell,
    cell_adjacency,
    order_ring,
    section_cells,
    signed_area,
)

SPAN_Z = np.array([0.0, 0.0, 1.0])


def _wall_graph(points, edges, shear_stiffness):
    """Build ``(cells, adjacency, edge_length, edge_shear_stiffness)`` from a wall graph.

    ``shear_stiffness`` is either a scalar (applied to every edge) or a mapping from
    the undirected edge key ``(min(a, b), max(a, b))`` to a membrane shear stiffness
    ``S = G*t`` in N/m.
    """
    pts = np.asarray(points, dtype=float)
    cells = section_cells(pts, edges, SPAN_Z)
    adjacency = {key: tuple(value) for key, value in cell_adjacency(cells).items()}
    edge_length: dict[tuple[int, int], float] = {}
    edge_shear_stiffness: dict[tuple[int, int], float] = {}
    for a, b in edges:
        key = (int(a), int(b)) if a < b else (int(b), int(a))
        edge_length[key] = float(np.linalg.norm(pts[a] - pts[b]))
        if hasattr(shear_stiffness, "keys"):
            edge_shear_stiffness[key] = float(shear_stiffness[key])
        else:
            edge_shear_stiffness[key] = float(shear_stiffness)
    return cells, adjacency, edge_length, edge_shear_stiffness


def _cell_with_node(cells, node):
    return next(i for i, c in enumerate(cells) if node in c.boundary.tolist())


def _compatibility_residual(
    cells, adjacency, edge_length, edge_shear_stiffness, index, flows, theta_rate
):
    """``sum_wall q_wall * ell / S - 2 * A_i * theta'`` for cell ``index``.

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
        total += wall_flow * edge_length[key] / edge_shear_stiffness[key]
    return total - 2.0 * cell.area * theta_rate


def test_single_cell_ring_matches_bredt():
    """One closed cell must reduce exactly to the uniform-S Bredt formula.

    Closed form: ``q = T / (2 A)`` and ``theta' = T P / (4 A^2 S)``, with ``A``
    the enclosed area, ``P`` the perimeter of the same ring and ``S = G*t`` the
    uniform membrane shear stiffness.  With one ``S`` per wall this is the old
    thickness-only ``phi = T P / (4 A^2 t)`` after dividing by ``G``.
    """
    points = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    shear_stiffness, torsion = 0.01, 5.0  # S = G*t; uniform, so flows match the old ones

    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(
        points, edges, shear_stiffness
    )
    assert len(cells) == 1
    area = cells[0].area
    perimeter = sum(edge_length.values())

    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )

    assert flows[0] == pytest.approx(torsion / (2.0 * area), rel=1e-12)
    assert theta_rate == pytest.approx(
        torsion * perimeter / (4.0 * area**2 * shear_stiffness), rel=1e-12
    )


def test_symmetric_two_cell_box_web_carries_no_net_flow():
    """Equal cells sharing one web: the web is inert.

    Bredt-Batho on the symmetric two-cell box gives equal cell flows, so the
    shared web carries ``q_1 - q_2 = 0`` and ``theta'`` collapses to the
    single-cell value ``T P_outer / (4 A_outer^2 S)`` with ``A_outer`` the whole
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
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(points, edges, thickness)
    assert len(cells) == 2

    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )

    assert flows[0] == pytest.approx(flows[1], rel=1e-12)
    left, right = adjacency[(4, 5)]
    assert flows[left] - flows[right] == pytest.approx(0.0, abs=1e-12)

    outer_area = width * height
    outer_perimeter = 2.0 * (width + height)
    assert theta_rate == pytest.approx(
        torsion * outer_perimeter / (4.0 * outer_area**2 * thickness), rel=1e-12
    )


def test_asymmetric_two_cell_box_matches_independent_hand_assembly():
    """Independent 3x3 reference for the multi-cell system.

    Geometry: a 3 x 2 box with a web at x = 1, uniform shear stiffness S; cell
    areas A0 = 2 (x in [0, 1]) and A1 = 4 (x in [1, 3]).  The reference assembles

        [ D0     -d_web   -2 A0 ] [ q0    ]   [ 0 ]
        [ -d_web  D1      -2 A1 ] [ q1    ] = [ 0 ]
        [ 2 A0   2 A1      0    ] [ theta']   [ T ]

    from the wall lengths alone, with D0 = 6/S, D1 = 8/S and d_web = 2/S, and
    solves it with ``numpy.linalg.solve``.
    """
    shear_stiffness, torsion = 0.01, 7.0
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
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(
        points, edges, shear_stiffness
    )
    assert len(cells) == 2
    area0, area1 = cells[0].area, cells[1].area
    assert area0 == 2.0 and area1 == 4.0

    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )

    # Independent reference: D0 = 6/S, D1 = 8/S, web = 2/S (from the lengths).
    reference = np.array(
        [
            [6.0 / shear_stiffness, -2.0 / shear_stiffness, -2.0 * area0],
            [-2.0 / shear_stiffness, 8.0 / shear_stiffness, -2.0 * area1],
            [2.0 * area0, 2.0 * area1, 0.0],
        ]
    )
    reference_solution = np.linalg.solve(reference, np.array([0.0, 0.0, torsion]))
    reference_q = reference_solution[:2]
    reference_theta_rate = reference_solution[2]

    assert flows[0] == pytest.approx(reference_q[0], rel=1e-10)
    assert flows[1] == pytest.approx(reference_q[1], rel=1e-10)
    assert theta_rate == pytest.approx(reference_theta_rate, rel=1e-10)

    # Exact invariants, independent of the reference matrix.
    assert 2.0 * (area0 * flows[0] + area1 * flows[1]) == pytest.approx(torsion, rel=1e-10)
    for index in (0, 1):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_shear_stiffness, index, flows, theta_rate
        ) == pytest.approx(0.0, abs=1e-8)
    left, right = adjacency[(4, 5)]  # the shared web, owners in cell order
    assert left != right


def test_three_cell_blade_topology_shares_one_theta_rate():
    """Blade topology: outer skin plus two parallel webs, three equal cells.

    For the parallel-web rectangle the reflection symmetry only forces the two
    outer cells to be equal (``q_left = q_right``); the middle cell is coupled to
    a neighbour by *two* webs while each outer cell has one web plus a free side
    wall, so ``q_middle`` differs.  The physical invariants that must hold
    exactly are the torque equilibrium and one common ``theta'`` per cell.
    """
    shear_stiffness, torsion = 0.01, 9.0
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
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(
        points, edges, shear_stiffness
    )
    assert len(cells) == 3
    assert [c.area for c in cells] == [1.0, 1.0, 1.0]

    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )

    torque = 2.0 * sum(cell.area * flow for cell, flow in zip(cells, flows, strict=True))
    assert torque == pytest.approx(torsion, rel=1e-12)
    for index in range(3):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_shear_stiffness, index, flows, theta_rate
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
    radius, shear_stiffness, torsion = 1.0, 0.01, 9.0
    angles = np.deg2rad([90.0, 210.0, 330.0])
    points = np.zeros((4, 3))
    points[1:, 0] = radius * np.cos(angles)
    points[1:, 1] = radius * np.sin(angles)
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1)]
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(
        points, edges, shear_stiffness
    )
    assert len(cells) == 3
    assert cells[0].area == pytest.approx(cells[1].area, rel=1e-12)
    assert cells[1].area == pytest.approx(cells[2].area, rel=1e-12)

    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )

    assert flows[0] == pytest.approx(flows[1], rel=1e-12)
    assert flows[1] == pytest.approx(flows[2], rel=1e-12)
    for key in [(0, 1), (0, 2), (0, 3)]:
        i, j = adjacency[key]
        assert flows[i] - flows[j] == pytest.approx(0.0, abs=1e-12)
    for index in range(3):
        assert _compatibility_residual(
            cells, adjacency, edge_length, edge_shear_stiffness, index, flows, theta_rate
        ) == pytest.approx(0.0, abs=1e-9)


def test_thicker_walled_cell_carries_the_larger_share():
    """Two equal-area cells; the thicker-walled one must carry more flow.

    For equal areas and a common ``theta'``, a cell with thicker walls has the
    larger ``S = G*t`` and so the smaller compliance ``sum ell / S``, so it needs
    the larger flow to reach the same ``2 A theta'``.  Only the direction is
    pinned, no magnitude.
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

    thin = thickness  # S = G*t with G = 1, so the numbers are the old thicknesses
    thick = 2.0 * thickness

    def key(a, b):
        return (a, b) if a < b else (b, a)

    by_edge = {key(a, b): thin for a, b in edges}
    # Right cell (bounded by nodes 1, 2, 5, 4) gets the thicker walls; the web
    # stays thin so the two cells differ only in their outer skin stiffness.
    by_edge[key(4, 1)] = thick
    by_edge[key(1, 2)] = thick
    by_edge[key(2, 5)] = thick

    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(points, edges, by_edge)
    assert len(cells) == 2

    flows, _ = multi_cell_shear_flow(cells, adjacency, edge_length, edge_shear_stiffness, torsion)

    thin_cell = _cell_with_node(cells, 0)  # left cell: thin outer skin
    thick_cell = _cell_with_node(cells, 1)  # right cell: thick outer skin
    assert flows[thin_cell] != pytest.approx(flows[thick_cell])
    assert flows[thick_cell] > flows[thin_cell]


def test_degenerate_inputs_raise_value_error():
    """The solver rejects empty, non-positive, missing and singular inputs."""
    rectangle = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(rectangle, edges, 0.01)

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

    zero_stiffness = dict(edge_shear_stiffness)
    zero_stiffness[(0, 1)] = 0.0
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, zero_stiffness, 1.0)

    missing_length = {key: value for key, value in edge_length.items() if key != (0, 1)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, missing_length, edge_shear_stiffness, 1.0)

    missing_stiffness = {key: value for key, value in edge_shear_stiffness.items() if key != (0, 1)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, missing_stiffness, 1.0)

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
    edge_shear_stiffness = {(0, 1): 0.01, (1, 2): 0.01, (0, 2): 0.01}
    adjacency = {(0, 1): (0, 1, 2), (1, 2): (0, 1, 2), (0, 2): (0, 1, 2)}
    with pytest.raises(ValueError):
        multi_cell_shear_flow(cells, adjacency, edge_length, edge_shear_stiffness, 1.0)


def test_inputs_are_not_mutated():
    """The solver is pure: cells, adjacency, lengths and stiffnesses survive."""
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
    cells, adjacency, edge_length, edge_shear_stiffness = _wall_graph(points, edges, 0.01)

    adjacency_before = {key: list(value) for key, value in adjacency.items()}
    length_before = dict(edge_length)
    stiffness_before = dict(edge_shear_stiffness)
    boundaries_before = [cell.boundary.copy() for cell in cells]

    multi_cell_shear_flow(cells, adjacency, edge_length, edge_shear_stiffness, 7.0)

    assert {key: list(value) for key, value in adjacency.items()} == adjacency_before
    assert edge_length == length_before
    assert edge_shear_stiffness == stiffness_before
    for cell, boundary in zip(cells, boundaries_before, strict=True):
        assert np.array_equal(cell.boundary, boundary)


# -----------------------------------------------------------------------------
# Part B/C: the per-ring wall graph (``ring_section``) and its stiffness
# -----------------------------------------------------------------------------


def _box_ring_mesh():
    """A one-ring, two-cell box section built as quad shell elements.

    Ring A (z = 0) and ring B (z = 1) share the same six in-plane nodes.  Each of
    the box graph's seven wall edges becomes one quad ``[A_a, A_b, B_b, B_a]``,
    so the ring's chordwise edges are exactly the box's wall graph (six outer
    skin segments plus the central web) and ``ring_section`` can recover it.

    Returns ``(mesh, ring_nodes, wall_edges)`` with ``ring_nodes`` the positional
    indices of ring A.
    """
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
    from aeroelast.core.mesh.model import MeshModel

    Node._id_counter = 0
    MeshElement._id_counter = 0
    outline = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 1.0),
        (0.0, 1.0),
        (1.0, 0.0),
        (1.0, 1.0),
    ]
    mesh = MeshModel()
    ring_a = [Node([x, y, 0.0]) for x, y in outline]
    ring_b = [Node([x, y, 1.0]) for x, y in outline]
    for node in ring_a + ring_b:
        mesh.add_node(node)

    wall_edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]
    element_set = ElementSet("box")
    for a, b in wall_edges:
        element = MeshElement([ring_a[a], ring_a[b], ring_b[b], ring_b[a]], ElementType.quad)
        mesh.add_element(element)
        element_set.add_element(element)
    mesh.add_element_set(element_set)
    return mesh, np.arange(len(ring_a), dtype=np.intp), wall_edges


def test_ring_section_two_cell_box_recovers_the_wall_graph_and_closed_form():
    """The ring's own wall graph and the symmetric two-cell closed form.

    With a uniform ``S`` the central web is inert, so both cells carry
    ``q = T / (2 * A_outer)`` and ``theta' = T * P_outer / (4 * A_outer**2 * S)``
    (the single-cell Bredt formula on the outer boundary).  ``ring_section`` must
    recover the two cells and their edge data from the mesh alone.
    """
    mesh, ring_nodes, _ = _box_ring_mesh()
    cells, adjacency, edge_length, edge_shear_stiffness = ring_section(mesh, ring_nodes, SPAN_Z)

    assert len(cells) == 2
    outer_area = 2.0 * 1.0
    outer_perimeter = 2.0 * (2.0 + 1.0)
    assert sum(cell.area for cell in cells) == pytest.approx(outer_area, rel=1e-12)
    assert set(edge_shear_stiffness.values()) == {1.0}

    torsion = 5.0
    flows, theta_rate = multi_cell_shear_flow(
        cells, adjacency, edge_length, edge_shear_stiffness, torsion
    )
    assert flows[0] == pytest.approx(torsion / (2.0 * outer_area), rel=1e-12)
    assert flows[1] == pytest.approx(torsion / (2.0 * outer_area), rel=1e-12)
    assert theta_rate == pytest.approx(
        torsion * outer_perimeter / (4.0 * outer_area**2 * 1.0), rel=1e-12
    )
    shared = next(key for key, owners in adjacency.items() if len(owners) == 2)
    left, right = adjacency[shared]
    assert flows[left] - flows[right] == pytest.approx(0.0, abs=1e-12)


def test_ring_section_without_properties_uses_uniform_unit_stiffness():
    """``element_properties=None`` (and missing entries) fall back to ``S = 1.0``.

    The cells and the edge lengths are unchanged by the property source; only
    the per-edge stiffness differs.  A supplied isotropic property gives the
    expected ``A66 = G * t`` through the ``(2, 2)`` entry of the merged membrane
    stiffness.
    """
    mesh, ring_nodes, _ = _box_ring_mesh()
    cells, adjacency, edge_length, edge_shear = ring_section(mesh, ring_nodes, SPAN_Z)
    assert all(value == 1.0 for value in edge_shear.values())

    e_modulus, nu, thickness = 70.0e9, 0.33, 0.02
    expected = e_modulus / (2.0 * (1.0 + nu)) * thickness
    props = {
        "box": {
            "type": "isotropic",
            "e": e_modulus,
            "nu": nu,
            "rho": 2700.0,
            "thickness": thickness,
            "shear_correction": 5.0 / 6.0,
        }
    }
    cells_prop, adjacency_prop, length_prop, shear_prop = ring_section(
        mesh, ring_nodes, SPAN_Z, props
    )

    assert [cell.area for cell in cells_prop] == [cell.area for cell in cells]
    assert length_prop == edge_length
    assert adjacency_prop == adjacency
    assert all(value == pytest.approx(expected, rel=1e-12) for value in shear_prop.values())

    # A property map that does not name the ring's element set is the documented
    # "no entry" fallback, not an error.
    _, _, _, shear_missing = ring_section(mesh, ring_nodes, SPAN_Z, {"other": props["box"]})
    assert all(value == 1.0 for value in shear_missing.values())


def test_ring_section_rejects_degenerate_rings():
    """Fewer than three ring nodes, or a wall graph with no closed cell, raise."""
    from aeroelast.core.mesh.entities import ElementType, MeshElement, Node
    from aeroelast.core.mesh.model import MeshModel

    Node._id_counter = 0
    MeshElement._id_counter = 0
    mesh = MeshModel()
    nodes = [
        Node([0.0, 0.0, 0.0]),
        Node([1.0, 0.0, 0.0]),
        Node([1.0, 1.0, 0.0]),
        Node([0.0, 0.0, 1.0]),
    ]
    for node in nodes:
        mesh.add_node(node)

    with pytest.raises(ValueError):
        ring_section(mesh, np.array([0, 1], dtype=np.intp), SPAN_Z)

    # One triangle spans the ring: only its chordwise edge (0, 1) lies in the ring,
    # so the wall graph is an open chain and bounds no cell.
    mesh.add_element(MeshElement([nodes[0], nodes[1], nodes[3]], ElementType.triangle))
    with pytest.raises(ValueError):
        ring_section(mesh, np.arange(3, dtype=np.intp), SPAN_Z)


def test_ring_section_on_iea15mw_midspan_ring():
    """The repo's own IEA-15MW blade: a mid-span ring has three cells.

    Building the blade fixture costs ~20 s.  The count and the positivity of each
    cell area are asserted; the area sum against the ring's own shoelace area and
    the per-edge stiffness spread are printed only (the sum is a few percent,
    not exact, because the web chords cut the outer shoelace polygon).
    """
    pytest.importorskip("_aeroelast", reason="Rust backend not available")
    from aeroelast.core.mesh.entities import MeshElement, Node
    from aeroelast.models.blade.model import Blade
    from aeroelast.solvers.bem.force_projection import RING_GAP_FRACTION
    from aeroelast.solvers.bem.section_contour import rings_in_band
    from tests.support.paths import DATA_DIR

    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(DATA_DIR / "IEA-15-240-RWT.yaml"), element_size=1.0)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None
    properties = model.get_element_properties()

    coords = mesh.coords_array
    span = coords @ SPAN_Z
    tolerance = RING_GAP_FRACTION * float(span.max() - span.min())
    groups = rings_in_band(span, tolerance)
    ring_nodes = groups[len(groups) // 2]

    cells, _, edge_length, edge_shear = ring_section(mesh, ring_nodes, SPAN_Z, properties)
    assert len(cells) == 3
    assert all(cell.area > 0.0 for cell in cells)

    ring_points = coords[ring_nodes]
    shoelace = signed_area(ring_points[order_ring(ring_points, SPAN_Z)], SPAN_Z)
    cell_sum = sum(cell.area for cell in cells)
    print(
        f"\n[midspan ring] {len(ring_nodes)} nodes, {len(cells)} cells; "
        f"sum(cell area) = {cell_sum:.6e} m^2 vs ring shoelace = {shoelace:.6e} m^2 "
        f"(relative difference {abs(cell_sum - shoelace) / abs(shoelace):.4%})"
    )
    shear_values = np.asarray([edge_shear[key] for key in sorted(edge_shear)])
    print(
        f"[midspan ring] per-edge shear stiffness S = G*t: "
        f"min {shear_values.min():.6e} N/m, max {shear_values.max():.6e} N/m "
        f"(spread {shear_values.max() / shear_values.min():.3f}x)"
    )
    assert len(edge_length) == len(edge_shear) >= 3


# -----------------------------------------------------------------------------
# Part D: the production project() path realises the multi-cell flow
# (issue #11, ODD T2c-2b-iii)
# -----------------------------------------------------------------------------

#: Round-off guard for the tolerance-free invariances below.  The net force and the
#: realised span moment are exact in real arithmetic; the per-cell walk, the solver's
#: equilibrium and the independent reference agree to floating-point round-off.
ROUNDOFF = 1e-9


def _asymmetric_two_cell_ring_mesh():
    """One ring of an **asymmetric** two-cell box, built as quad shell elements.

    Same construction as :func:`_box_ring_mesh` (a ring A at z = 0 and a ring B at
    z = 1 sharing six in-plane nodes, one quad per wall edge), but the web sits at
    x = 1 in a width-3 rectangle instead of at the centre of a width-2 one.  The
    symmetry matters: a symmetric box's shared web carries zero net flow, so the
    single-cell outer-boundary flow reproduces the multi-cell one exactly and the
    test would be blind (the ODD T2c-2b note says so).  The two cells here have
    areas 2 and 4 m^2 and different flows.

    Returns ``(mesh, wall_edges)`` with ``wall_edges`` the ring's local wall graph.
    """
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
    from aeroelast.core.mesh.model import MeshModel

    Node._id_counter = 0
    MeshElement._id_counter = 0
    outline = [
        (0.0, 0.0),
        (3.0, 0.0),
        (3.0, 2.0),
        (0.0, 2.0),
        (1.0, 0.0),
        (1.0, 2.0),
    ]
    mesh = MeshModel()
    ring_a = [Node([x, y, 0.0]) for x, y in outline]
    ring_b = [Node([x, y, 1.0]) for x, y in outline]
    for node in ring_a + ring_b:
        mesh.add_node(node)

    wall_edges = [(0, 4), (4, 1), (1, 2), (2, 5), (5, 3), (3, 0), (4, 5)]
    element_set = ElementSet("box")
    for a, b in wall_edges:
        element = MeshElement([ring_a[a], ring_a[b], ring_b[b], ring_b[a]], ElementType.quad)
        mesh.add_element(element)
        element_set.add_element(element)
    mesh.add_element_set(element_set)
    return mesh, wall_edges


def _independent_multicell_flows(cells, adjacency, edge_length, stiffness, torsion):
    """The ``(n + 1)`` Bredt-Batho system assembled by hand and solved with numpy.

    This is the independent reference: the ``numpy.linalg.solve`` of the matrix the
    solver's own docstring writes down, built here from the wall geometry alone.
    It never calls :func:`multi_cell_shear_flow`.
    """
    n = len(cells)
    matrix = np.zeros((n + 1, n + 1))
    rhs = np.zeros(n + 1)
    for i, cell in enumerate(cells):
        for a, b in cell.edges:
            key = (a, b) if a < b else (b, a)
            ratio = edge_length[key] / stiffness[key]
            matrix[i, i] += ratio
            for owner in adjacency[key]:
                if owner != i:
                    matrix[i, owner] -= ratio
        matrix[i, n] = -2.0 * cell.area
        matrix[n, i] = 2.0 * cell.area
    rhs[n] = torsion
    solution = np.linalg.solve(matrix, rhs)
    return list(solution[:n])


def _per_cell_wall_field(cells, flows, offsets):
    """The nodal field of a set of cell flows, by walking each cell's own boundary.

    A shared wall is traversed by the two cells in opposite directions, so it
    collects ``q_i - q_j`` with no special case.  This is the expectation the
    production realisation must reproduce.
    """
    field = np.zeros_like(offsets)
    for cell, flow in zip(cells, flows, strict=True):
        for a, b in cell.edges:
            edge = offsets[b] - offsets[a]
            ell = float(np.linalg.norm(edge))
            if ell == 0.0:
                continue
            half = 0.5 * (flow * ell) * (edge / ell)
            field[a] = field[a] + half
            field[b] = field[b] + half
    return field


def _minimum_norm_span_field(offsets, span_hat, moment_span):
    """The constrained minimum-norm nodal field for a pure span moment.

    ``_distribute`` spreads a pure span moment as the circle-tangential field
    ``f_j = omega x d_j`` (the tube test documents the residual is 2.3e-13), with
    ``omega`` set by ``sum_j d_j x f_j = M``; for a planar ring
    ``omega = M / sum_j |d_perp,j|^2``.  Derived here from the geometry, never from
    the projector.
    """
    offsets = np.asarray(offsets, dtype=float)
    in_plane = offsets - np.outer(offsets @ span_hat, span_hat)
    inertia = float(np.sum(in_plane**2))
    omega = moment_span / inertia
    return omega * np.cross(span_hat, offsets)


def test_production_projector_realises_multicell_flow_on_a_two_cell_ring():
    """``project()`` carries a pure span moment as the multi-cell wall flow.

    The mesh is the asymmetric two-cell box ring above; a two-station ``BladeAero``
    makes each span ring its own BEM strip, so the strip's one physical ring is the
    two-cell section and its requested span moment is ``Mp * dr``.  The assertions
    are the tolerance-free invariances of the realisation (net force, realised
    moment) plus the decisive one: the nodal field equals the multi-cell wall flow
    built independently from :func:`section_cells` and a hand-assembled Bredt-Batho
    system, and is **not** the minimum-norm circle-tangential field.  That last pair
    is what fails if the projection falls back to the old single-cell/min-norm path.
    """
    from tests.validation.bem.test_force_projection import (
        _make_simple_blade_aero,
        _strip_widths,
    )

    mesh, wall_edges = _asymmetric_two_cell_ring_mesh()
    blade_aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=1.0)
    span_hat = SPAN_Z
    moment_per_length = 7.0
    bem = BEMResult(
        r=blade_aero.r,
        Np=np.zeros(2),
        Tp=np.zeros(2),
        alpha=np.zeros(2),
        cl=np.zeros(2),
        cd=np.zeros(2),
        a=np.zeros(2),
        ap=np.zeros(2),
        thrust=0.0,
        torque=0.0,
        power=0.0,
        Mp=np.full(2, moment_per_length),
    )
    props = {
        "box": {
            "type": "isotropic",
            "e": 70.0e9,
            "nu": 0.33,
            "rho": 2700.0,
            "thickness": 0.02,
            "shear_correction": 5.0 / 6.0,
        }
    }

    projector = ForceProjector(mesh, blade_aero, element_properties=props)
    forces = projector.project(bem)
    coords = mesh.coords_array

    widths = _strip_widths(blade_aero, mesh, span_hat)
    assert widths.shape == (2,)

    # Each ring is its own strip: ring A is nodes 0..5 (z = 0), ring B nodes 6..11
    # (z = 1).  The requested span moment of a strip is Mp * dr, derived from the
    # station grid, never from the projector.
    for k, strip_nodes in enumerate((np.arange(6, dtype=np.intp), np.arange(6, 12, dtype=np.intp))):
        # The box's section is symmetric, so ``_section_ends`` keeps the leading edge at
        # the HIGH end and ``_strip_moment_axis_sign`` is -1: the polars' nose-up Cm rides
        # ``-span_hat`` here (issue #11).  The requested span moment carries that sign; the
        # realisation mechanics below are unchanged.
        requested_span_moment = float(
            moment_per_length * widths[k] * projector._strip_moment_axis_sign[k]
        )
        ring_points = coords[strip_nodes]
        centroid = ring_points.mean(axis=0)
        offsets = ring_points - centroid
        strip_forces = forces[strip_nodes]

        # 1. Self-equilibrated: the wall flow has no net force.
        assert np.abs(strip_forces.sum(axis=0)).max() <= ROUNDOFF, (
            f"net force of the realised multi-cell field is not zero: {strip_forces.sum(axis=0)}"
        )

        # 2. The realised moment about the strip centroid is the requested span
        #    moment, with no transverse component.
        realised_moment = np.cross(offsets, strip_forces).sum(axis=0)
        assert abs(float(realised_moment @ span_hat) - requested_span_moment) <= (
            ROUNDOFF * abs(requested_span_moment)
        ), (
            f"realised span moment {float(realised_moment @ span_hat):.6e} != "
            f"requested {requested_span_moment:.6e}"
        )
        assert np.abs(realised_moment - (realised_moment @ span_hat) * span_hat).max() <= (
            ROUNDOFF * abs(requested_span_moment)
        )

        # 3. The decisive expectation, built without the projector: the multi-cell
        #    wall flow from the hand-assembled system and the per-cell walk.
        cells, adjacency, edge_length, stiffness = _wall_graph(ring_points, wall_edges, 1.0)
        assert len(cells) == 2
        flows = _independent_multicell_flows(
            cells, adjacency, edge_length, stiffness, requested_span_moment
        )
        expected_field = _per_cell_wall_field(cells, flows, offsets)
        scale = float(np.abs(expected_field).max())
        np.testing.assert_allclose(strip_forces, expected_field, rtol=0.0, atol=ROUNDOFF * scale)

        # 3b. The multi-cell split is readable on the shared web: its net flow is
        #     ``q_i - q_j``.  The realisation pushes each wall's ``0.5*q*ell*t`` onto
        #     BOTH endpoints, so summing the two web nodes' forces and projecting on
        #     the web tangent recovers the net flow.  The two skin walls meeting the
        #     web at its nodes run perpendicular to it on this rectangular box, so
        #     they drop out of that projection.
        shared_key = next(key for key, owners in adjacency.items() if len(owners) == 2)
        wa, wb = shared_key
        owner_i, owner_j = adjacency[shared_key]
        mc_flows, _ = multi_cell_shear_flow(
            cells, adjacency, edge_length, stiffness, requested_span_moment
        )
        traverses_ab = any((int(u), int(v)) == (wa, wb) for u, v in cells[owner_i].edges)
        orientation = 1.0 if traverses_ab else -1.0
        expected_net_flow = orientation * (mc_flows[owner_i] - mc_flows[owner_j])
        web_edge = offsets[wb] - offsets[wa]
        web_length = float(np.linalg.norm(web_edge))
        web_tangent = web_edge / web_length
        read_back_flow = float((strip_forces[wa] + strip_forces[wb]) @ web_tangent) / web_length
        assert read_back_flow == pytest.approx(expected_net_flow, rel=ROUNDOFF), (
            f"shared-web net flow read back as {read_back_flow:.6e} N/m != q_i - q_j = "
            f"{expected_net_flow:.6e} N/m from multi_cell_shear_flow"
        )

        # 4. ... and it is not the minimum-norm field, i.e. the multi-cell path ran.
        minimum_norm = _minimum_norm_span_field(offsets, span_hat, requested_span_moment)
        assert np.abs(strip_forces - minimum_norm).max() > 1e-3 * scale, (
            "the projected field matches the minimum-norm circle-tangential field; "
            "the multi-cell realisation did not run"
        )


def _three_ring_single_cell_mesh():
    """Three single-cell square rings at z = 0, 0.4, 1.0, joined by quads.

    The non-uniform span positions are the point: ``0.5*(z1-z0)``,
    ``0.5*(z2-z0)`` and ``0.5*(z2-z1)`` are the trapezoidal tributary lengths and
    they are all different, so an equal share cannot be mistaken for the correct
    split.  Every ring bounds exactly one cell.
    """
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
    from aeroelast.core.mesh.model import MeshModel

    Node._id_counter = 0
    MeshElement._id_counter = 0
    square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    spans = [0.0, 0.4, 1.0]
    mesh = MeshModel()
    rings = []
    for z in spans:
        ring = [Node([x, y, z]) for x, y in square]
        for node in ring:
            mesh.add_node(node)
        rings.append(ring)

    element_set = ElementSet("tube")
    walls = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for lower, upper in zip(rings[:-1], rings[1:], strict=True):
        for a, b in walls:
            element = MeshElement([lower[a], lower[b], upper[b], upper[a]], ElementType.quad)
            mesh.add_element(element)
            element_set.add_element(element)
    mesh.add_element_set(element_set)
    return mesh, rings


def test_projector_splits_one_strip_over_three_single_cell_rings_by_tributary_weight():
    """A BEM strip wider than one physical ring shares its torsion by tributary span.

    The three single-cell square rings sit at z = 0, 0.4, 1.0 inside ONE strip, so
    the old single-ring gate would fall back and an equal share would hand each ring
    ``torsion / 3``.  The trapezoidal tributary rule gives the shares
    ``(0.2, 0.5, 0.3)`` from the rings' own span positions; each ring's realised span
    moment read back from the production nodal field must equal its share.
    """
    from tests.validation.bem.test_force_projection import (
        _make_simple_blade_aero,
        _strip_widths,
    )

    mesh, _ = _three_ring_single_cell_mesh()
    blade_aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=4.0)
    moment_per_length = 7.0
    bem = BEMResult(
        r=blade_aero.r,
        Np=np.zeros(2),
        Tp=np.zeros(2),
        alpha=np.zeros(2),
        cl=np.zeros(2),
        cd=np.zeros(2),
        a=np.zeros(2),
        ap=np.zeros(2),
        thrust=0.0,
        torque=0.0,
        power=0.0,
        Mp=np.full(2, moment_per_length),
    )
    props = {
        "tube": {
            "type": "isotropic",
            "e": 70.0e9,
            "nu": 0.33,
            "rho": 2700.0,
            "thickness": 0.02,
            "shear_correction": 5.0 / 6.0,
        }
    }
    projector = ForceProjector(mesh, blade_aero, element_properties=props)

    # One physical strip holds all three rings; each ring is a usable single cell.
    assert len(projector._strips[0].node_indices) == 12
    assert len(projector._strip_ring_groups[0]) == 3
    raw_sections = projector._strip_ring_sections[0]
    assert len(raw_sections) == 3
    assert all(section is not None and len(section.cells) == 1 for section in raw_sections)
    sections = [section for section in raw_sections if section is not None]
    assert len(sections) == 3

    forces = projector.project(bem)
    widths = _strip_widths(blade_aero, mesh, SPAN_Z)
    # Symmetric square section -> ``_section_ends`` keeps the leading edge at the HIGH
    # end, so ``_strip_moment_axis_sign`` is -1 and the requested torsion carries that
    # sign (issue #11); the tributary split below is unchanged.
    torsion = float(moment_per_length * widths[0] * projector._strip_moment_axis_sign[0])

    strip = projector._strips[0]
    positions = np.array(
        [
            float((strip.centroid + strip.offsets[section.local_nodes]).mean(axis=0) @ SPAN_Z)
            for section in sections
        ]
    )
    np.testing.assert_allclose(positions, [0.0, 0.4, 1.0])
    # Trapezoidal tributary lengths, half a cell at each end, normalised to one.
    shares = np.array(
        [
            0.5 * (positions[1] - positions[0]),
            0.5 * (positions[2] - positions[0]),
            0.5 * (positions[2] - positions[1]),
        ]
    )
    shares = shares / shares.sum()
    np.testing.assert_allclose(shares, [0.2, 0.5, 0.3])

    for section, share in zip(sections, shares, strict=True):
        local = section.local_nodes
        ring_forces = forces[strip.node_indices[local]]
        ring_offsets = strip.offsets[local]
        realised = float(np.cross(ring_offsets, ring_forces).sum(axis=0) @ SPAN_Z)
        assert realised == pytest.approx(torsion * share, rel=ROUNDOFF), (
            f"ring near z={float((strip.centroid + ring_offsets.mean(axis=0)) @ SPAN_Z):.3f} "
            f"realised {realised:.6e} N.m, not its tributary share {torsion * share:.6e} N.m"
        )


def test_projector_without_element_properties_falls_back_to_minimum_norm():
    """The multi-cell path is gated on real element properties, explicitly.

    ``ring_section`` returns a uniform ``S = 1.0`` when no properties are supplied,
    so a geometric-only split would look plausible but is not a physical split (the
    blade's per-wall ``G*t`` spans 18.8x).  The gate is therefore the *presence* of
    the property map, not an inference from the uniform stiffness: a projector built
    without ``element_properties`` routes the whole strip through ``_distribute``,
    and the realised field is the minimum-norm circle-tangential one, not the
    geometric multi-cell one.
    """
    from tests.validation.bem.test_force_projection import (
        _make_simple_blade_aero,
        _strip_widths,
    )

    mesh, wall_edges = _asymmetric_two_cell_ring_mesh()
    blade_aero = _make_simple_blade_aero(n_stations=2, hub_radius=3.0, span_length=1.0)
    moment_per_length = 7.0
    bem = BEMResult(
        r=blade_aero.r,
        Np=np.zeros(2),
        Tp=np.zeros(2),
        alpha=np.zeros(2),
        cl=np.zeros(2),
        cd=np.zeros(2),
        a=np.zeros(2),
        ap=np.zeros(2),
        thrust=0.0,
        torque=0.0,
        power=0.0,
        Mp=np.full(2, moment_per_length),
    )
    projector = ForceProjector(mesh, blade_aero)  # deliberately NO element_properties
    assert all(
        section is None or not getattr(section, "from_element_properties", False)
        for sections in projector._strip_ring_sections
        for section in sections
    ), "no-properties rings must be marked non-physical so the gate reroutes them"

    forces = projector.project(bem)
    coords = mesh.coords_array
    widths = _strip_widths(blade_aero, mesh, SPAN_Z)
    for k, strip_nodes in enumerate((np.arange(6, dtype=np.intp), np.arange(6, 12, dtype=np.intp))):
        # The same signed requested moment as the properties-bearing test above (issue
        # #11): the fallback realisation must reproduce whatever sign the projector applies.
        requested = float(moment_per_length * widths[k] * projector._strip_moment_axis_sign[k])
        ring_points = coords[strip_nodes]
        offsets = ring_points - ring_points.mean(axis=0)
        minimum_norm = _minimum_norm_span_field(offsets, SPAN_Z, requested)
        scale = float(np.abs(minimum_norm).max())
        np.testing.assert_allclose(
            forces[strip_nodes], minimum_norm, rtol=0.0, atol=ROUNDOFF * scale
        )

        cells, adjacency, edge_length, stiffness = _wall_graph(ring_points, wall_edges, 1.0)
        flows, _ = multi_cell_shear_flow(cells, adjacency, edge_length, stiffness, requested)
        multicell = _per_cell_wall_field(cells, flows, offsets)
        assert np.abs(forces[strip_nodes] - multicell).max() > 1e-3 * np.abs(multicell).max(), (
            "the no-properties field still matches the geometric multi-cell field; "
            "the missing-properties gate did not reroute it to _distribute"
        )


# -----------------------------------------------------------------------------
# Part E: the applied pitching-moment axis (issue #11, the sign chain)
#
# The BEM polars define Cm (and hence Mp) about the aerofoil chord **leading to
# trailing edge**, positive nose-up, so the aerodynamic moment axis is
# ``chord_le_to_te x n_hat``.  ``ForceProjector`` resolves ``_strip_chord_dirs`` from the
# tangential-force sense, which need not run leading-to-trailing, so it records the
# per-strip conversion ``_strip_moment_axis_sign``.  ``tools/diagnose_sign_chain.py``
# measures the IEA-15MW deck's real stations (every one has LE . c_hat = -1); the cheap
# test below pins the wiring on a synthetic asymmetric section, and the slow one pins
# the convention-free invariant on the real blade.
# -----------------------------------------------------------------------------


def _asymmetric_ring_mesh():
    """A closed section blunt at x=0 and sharp at x=1, repeated at z=0 and z=1.

    ``_section_ends`` reads the leading edge from the blunt end - the larger in-plane
    spread inside the outer quarter of the chord.  This polygon spreads +/-0.20 at x=0
    and only +/-0.02 at x=1, so with the chord resolved along +x the leading edge is the
    LOW end and the aerodynamic moment axis is ``+span``; reversing the chord axis through
    ``tangential_direction`` moves the blunt end to the HIGH end and flips the sign.
    """
    from aeroelast.core.mesh.entities import ElementSet, ElementType, MeshElement, Node
    from aeroelast.core.mesh.model import MeshModel

    Node._id_counter = 0
    MeshElement._id_counter = 0
    outline = [
        (0.0, 0.20),
        (0.5, 0.05),
        (1.0, 0.02),
        (1.0, -0.02),
        (0.5, -0.05),
        (0.0, -0.20),
    ]
    mesh = MeshModel()
    rings = []
    for z in (0.0, 1.0):
        ring = [Node([x, y, z]) for x, y in outline]
        for node in ring:
            mesh.add_node(node)
        rings.append(ring)
    element_set = ElementSet("skin")
    for a in range(len(outline)):
        b = (a + 1) % len(outline)
        element = MeshElement(
            [rings[0][a], rings[0][b], rings[1][b], rings[1][a]], ElementType.quad
        )
        mesh.add_element(element)
        element_set.add_element(element)
    mesh.add_element_set(element_set)
    return mesh


def test_strip_moment_axis_sign_follows_section_ends_and_reverses_with_the_chord():
    """The stored moment-axis sign is ``_section_ends``, and the applied ``Mp`` term flips.

    No finite-element solve: the sign lives entirely in ``ForceProjector.__init__``.  Two
    checks, both convention-free:

    1. for every strip, ``_strip_moment_axis_sign[k]`` is exactly the sign implied by
       ``_section_ends`` (``-1`` when the leading edge is at the HIGH end of the chord
       projection, ``+1`` when it is at the low end);
    2. the applied pitching-moment term flips with the chord axis: with ``Np = Tp = 0``
       the only applied moment is ``M_ac = Mp * dr * sign * span_hat`` and the realisation
       preserves it, so the moment read back about each strip's own centroid must reverse.
    """
    from tests.validation.bem.test_force_projection import _make_simple_blade_aero

    mesh = _asymmetric_ring_mesh()
    blade_aero = _make_simple_blade_aero(n_stations=2, hub_radius=0.0, span_length=1.0)
    forward = ForceProjector(mesh, blade_aero, span_direction=SPAN_Z)
    reversed_chord = ForceProjector(
        mesh, blade_aero, span_direction=SPAN_Z, tangential_direction=[-1.0, 0.0, 0.0]
    )

    for projector in (forward, reversed_chord):
        for k, strip in enumerate(projector._strips):
            ring_pts = projector._station_ring_points(k, strip, SPAN_Z)
            chord_hat = projector._strip_chord_dirs[k]
            le_i, _te_i = projector._section_ends(ring_pts, chord_hat, SPAN_Z)
            le_at_hi = le_i == int(np.argmax(ring_pts @ chord_hat))
            expected = -1.0 if le_at_hi else 1.0
            assert projector._strip_moment_axis_sign[k] == expected, (
                f"strip {k}: stored sign {projector._strip_moment_axis_sign[k]} but "
                f"_section_ends puts the leading edge at the "
                f"{'high' if le_at_hi else 'low'} end (expected {expected})"
            )

    # The blunt (leading) end is the low end of the +x chord, so the aerodynamic axis is
    # +span; reversing the chord moves it to the high end and flips the sign.
    assert forward._strip_moment_axis_sign == [1.0, 1.0]
    assert reversed_chord._strip_moment_axis_sign == [-1.0, -1.0]

    moment_per_length = 7.0
    bem = BEMResult(
        r=blade_aero.r,
        Np=np.zeros(2),
        Tp=np.zeros(2),
        alpha=np.zeros(2),
        cl=np.zeros(2),
        cd=np.zeros(2),
        a=np.zeros(2),
        ap=np.zeros(2),
        thrust=0.0,
        torque=0.0,
        power=0.0,
        Mp=np.full(2, moment_per_length),
    )

    def realised_span_moments(projector):
        forces = projector.project(bem)
        return [
            float(np.cross(strip.offsets, forces[strip.node_indices]).sum(axis=0) @ SPAN_Z)
            for strip in projector._strips
        ]

    forward_moments = realised_span_moments(forward)
    reversed_moments = realised_span_moments(reversed_chord)
    for forward_moment, reversed_moment in zip(forward_moments, reversed_moments, strict=True):
        assert forward_moment > 0.0 and reversed_moment < 0.0, (
            f"the Mp term did not flip with the chord axis: forward {forward_moment:.6e} "
            f"N.m, reversed {reversed_moment:.6e} N.m"
        )
        assert forward_moment == pytest.approx(-reversed_moment, rel=ROUNDOFF)


@pytest.mark.slow
def test_min_norm_and_with_properties_realisations_of_the_rated_moment_agree_in_sign():
    """Two faithful realisations of the SAME rated moment must rotate the tip the same way.

    The minimum-norm projector (no ``element_properties``) and the with-properties
    multi-cell projector are two different realisations of the same requested section
    moment.  A rotation is linear in the applied load, so two faithful realisations of one
    moment cannot give tip rotations of opposite sign; before the pitching-moment axis fix
    they did (``-0.8696`` vs ``+0.0773`` deg, ``tools/diagnose_pitching_moment_sense.py``),
    which is exactly the whole defect.  With ``_strip_moment_axis_sign`` both read
    ``+8.1048`` vs ``+9.6669`` deg.

    Cost: the real IEA-15MW blade mesh, the AeroDyn polars and **two** shell solves, so it
    is marked ``slow`` and keeps its own blade/mesh construction (no heavy import at
    module scope, so the cheap tests above run without ccblade).
    """
    pytest.importorskip("ccblade", reason="ccblade not installed (pip install -e '.[bem]')")
    pytest.importorskip("_aeroelast", reason="Rust backend not available")

    from _aeroelast import PyMeshAssembler
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve

    from aeroelast.models.blade.model import Blade
    from aeroelast.solvers.bem.engine import BEMSolver
    import tests.validation.blade.test_blade_iea15mw_validation as blade_validation
    import tests.validation.blade.test_blade_rated_twist as rated
    from tests.support.openfast_bem import build_blade_aero_from_aerodyn

    from aeroelast.core.mesh.entities import MeshElement, Node

    Node._id_counter = 0
    MeshElement._id_counter = 0
    model = Blade(str(rated.YAML), element_size=1.0)
    model.generate_mesh()
    mesh = model.mesh
    assert mesh is not None
    props = model.get_element_properties()
    assembler = PyMeshAssembler.from_model(
        blade_validation._to_rust_mesh(mesh, props),
        props,
        list(blade_validation.SPAN_DIRECTION),
        None,
    )
    n = assembler.dofs_count
    rows, cols, vals = assembler.assemble_k()
    K = coo_matrix((np.asarray(vals), (np.asarray(rows), np.asarray(cols))), shape=(n, n)).tocsr()
    root = {mesh.node_id_to_index[nid] for nid in mesh.get_node_set("RootNodes").node_ids}
    free = np.array(
        [i for i in range(n) if i not in {6 * r + d for r in root for d in range(6)}],
        dtype=np.int64,
    )
    coords = mesh.coords_array
    phys = rated._physical_stations(coords)
    tip = np.where(np.abs(coords[:, 2] - phys[-1]) < rated.STATION_GAP_TOLERANCE)[0]

    aero = build_blade_aero_from_aerodyn(rated.AD_PRIMARY)
    bem = BEMSolver(aero, rho=1.225, mu=1.81206e-5, hub_height=150.0, shear_exp=0.0).compute(
        rated.V_RATED, rated.RPM_RATED, rated.PITCH_RATED
    )
    assert bem.Mp is not None

    plain = ForceProjector(mesh, aero, span_direction=SPAN_Z)
    with_props = ForceProjector(mesh, aero, span_direction=SPAN_Z, element_properties=props)

    def tip_rotation(projector):
        forces = projector.project(bem)
        force = np.zeros(n)
        force[0::6] = forces[:, 0]
        force[1::6] = forces[:, 1]
        force[2::6] = forces[:, 2]
        u = np.zeros(n)
        u[free] = np.asarray(spsolve(K[np.ix_(free, free)], force[free])).ravel()
        return float(np.rad2deg(rated._ring_kinematics(coords, u, tip)["omega"]))

    minimum_norm = tip_rotation(plain)
    multi_cell = tip_rotation(with_props)
    print(
        f"\nrated pitching-moment realisations, same BEM result: minimum-norm "
        f"{minimum_norm:+.4f} deg, with-properties multi-cell {multi_cell:+.4f} deg"
    )
    assert minimum_norm * multi_cell > 0.0, (
        f"two faithful realisations of the same requested section moment rotate the tip "
        f"the opposite way ({minimum_norm:+.4f} deg vs {multi_cell:+.4f} deg): the applied "
        f"moment's sign chain is broken"
    )
