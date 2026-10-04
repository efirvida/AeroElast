"""Pure ordering of a blade/tube section ring for thin-walled shear flow.

A closed thin-walled section under a torque ``M`` about the span axis carries the
Bredt wall shear flow ``q = M / (2A)`` where ``A`` is the enclosed area.  That
formula is only meaningful on a ring that is *closed*, *simple* (no
self-intersections) and *consistently oriented*.  Stored FE node rings are none of
these in general, so this module:

* builds the section frame ``(u, v)`` perpendicular to the span axis,
* orders the ring by increasing polar angle about its in-plane centroid,
* checks closure / simplicity / area with exact predicates (no tolerances),
* reports whether the ordering is usable for ``q = M / (2A)``.

Sign convention
---------------
``section_plane_axes`` returns ``(u, v)`` with ``u x v = +span``.  A positive
shoelace area therefore means the ring runs counter-clockwise **as seen from
``+span``**, and the reported ``orientation`` is ``+1`` for that case.

Limitations
-----------
The angle sort is exact for a star-shaped/convex section and is only a
*candidate* for a concave one: a naive angular order can connect across a notch
that a hand-built wall contour would not.  Utilisation is therefore decided by
the checks below, never assumed.  Note also that for a set of distinct points the
angular order about the centroid is always a simple polygon, so a convex bow-tie
*input* is repaired by ordering; the checks flag rings that ordering cannot
repair (self-overlap / zero-length wall edges) and degenerate inputs.

Every function is pure and never mutates its inputs.  numpy only; no I/O, no
class state.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "ContourReport",
    "SectionCell",
    "section_plane_axes",
    "order_ring",
    "signed_area",
    "is_simple",
    "contour_report",
    "section_cells",
    "cell_adjacency",
    "merge_span_stations",
    "rings_in_band",
]


def _as_points(points: np.ndarray) -> np.ndarray:
    """Validate and copy an ``(n, 3)`` point array as float."""
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError(f"points must have shape (n, 3), got {arr.shape}")
    return arr


def _reference_axis(span: np.ndarray) -> np.ndarray:
    """Deterministic coordinate axis least aligned with ``span`` (lowest index wins)."""
    axis = np.zeros(3)
    axis[int(np.argmin(np.abs(span)))] = 1.0
    return axis


def section_plane_axes(span_dir: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return orthonormal ``(u, v)`` spanning the plane perpendicular to ``span_dir``.

    The frame is built so that ``u x v = +span``; consequently a positive
    shoelace area in ``(u, v)`` is counter-clockwise as seen from ``+span``.
    The choice is deterministic: the in-plane ``u`` is the coordinate axis least
    aligned with the span (ties broken by lowest index), orthogonalised against
    the span, and ``v = span x u``.

    Raises ``ValueError`` if ``span_dir`` is zero-length.
    """
    span = np.asarray(span_dir, dtype=float).reshape(3)
    norm = float(np.linalg.norm(span))
    if norm == 0.0:
        raise ValueError("span_dir must be non-zero")
    s = span / norm
    ref = _reference_axis(s)
    u = ref - s * float(ref @ s)
    u = u / float(np.linalg.norm(u))
    v = np.cross(s, u)
    return u, v


def _in_plane(points: np.ndarray, span_dir: np.ndarray) -> np.ndarray:
    """Project ``(n, 3)`` points into the section frame, returning ``(n, 2)``."""
    u, v = section_plane_axes(span_dir)
    return np.column_stack((points @ u, points @ v))


def order_ring(points: np.ndarray, span_dir: np.ndarray) -> np.ndarray:
    """Return an int permutation visiting ``points`` by increasing polar angle.

    The angle is measured about the points' in-plane centroid using
    :func:`section_plane_axes`, so the resulting ring is always the ``+span``-CCW
    one.  The tie-break is deterministic: equal angle -> smaller in-plane radius
    -> smaller original index.  The result is stable under a shuffle of the
    input (the ordered *point sequence* is unchanged).  A permutation is also
    returned for fewer than three points; the caller decides usability.

    As documented in the module, the angle sort is exact for a star-shaped or
    convex section and only a *candidate* for a concave one.
    """
    pts = _as_points(points)
    n = pts.shape[0]
    if n == 0:
        return np.empty(0, dtype=np.intp)
    u, v = section_plane_axes(span_dir)
    rel = pts - pts.mean(axis=0)
    pu = rel @ u
    pv = rel @ v
    angle = np.arctan2(pv, pu)
    radius = np.hypot(pu, pv)
    index = np.arange(n)
    return np.lexsort((index, radius, angle)).astype(np.intp)


def signed_area(points: np.ndarray, span_dir: np.ndarray) -> float:
    """Shoelace area of the closed polygon through ``points`` in their given order.

    Positive means counter-clockwise as seen from ``+span`` (see
    :func:`section_plane_axes`).  Fewer than three points yields ``0.0``.
    """
    pts = _as_points(points)
    if pts.shape[0] < 3:
        return 0.0
    plane = _in_plane(pts, span_dir)
    x = plane[:, 0]
    y = plane[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def _orient(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> float:
    """Twice the signed area of triangle ``(p, q, r)`` (exact zero test)."""
    return float((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))


def _on_segment(p: np.ndarray, q: np.ndarray, r: np.ndarray) -> bool:
    """Is ``r``, known collinear with ``pq``, inside segment ``pq``?"""
    return min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1])


def _segments_intersect(p1, p2, p3, p4) -> bool:
    """Exact inclusive segment intersection (proper crossings, touches, overlaps)."""
    d1 = _orient(p3, p4, p1)
    d2 = _orient(p3, p4, p2)
    d3 = _orient(p1, p2, p3)
    d4 = _orient(p1, p2, p4)
    if ((d1 > 0.0 and d2 < 0.0) or (d1 < 0.0 and d2 > 0.0)) and (
        (d3 > 0.0 and d4 < 0.0) or (d3 < 0.0 and d4 > 0.0)
    ):
        return True
    if d1 == 0.0 and _on_segment(p3, p4, p1):
        return True
    if d2 == 0.0 and _on_segment(p3, p4, p2):
        return True
    if d3 == 0.0 and _on_segment(p1, p2, p3):
        return True
    if d4 == 0.0 and _on_segment(p1, p2, p4):
        return True
    return False


def is_simple(points: np.ndarray, span_dir: np.ndarray) -> tuple[bool, tuple[int, int] | None]:
    """Is the closed polygon through ``points`` a simple ring?

    Returns ``(False, (i, j))`` for the first offending edge pair when two
    non-adjacent edges intersect (proper crossings plus collinear overlap or
    touching) or when an edge has zero in-plane length.  Adjacent edges sharing
    one endpoint are allowed.  Predicates are exact: comparisons are against
    ``0.0``, never an epsilon.  Fewer than three points is not a ring and returns
    ``(False, None)``.
    """
    pts = _as_points(points)
    n = pts.shape[0]
    if n < 3:
        return False, None
    plane = _in_plane(pts, span_dir)

    for i in range(n):
        j = (i + 1) % n
        if plane[i, 0] == plane[j, 0] and plane[i, 1] == plane[j, 1]:
            return False, (i, j)

    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                continue
            if _segments_intersect(plane[i], plane[(i + 1) % n], plane[j], plane[(j + 1) % n]):
                return False, (i, j)
    return True, None


def _clustered_spans(sorted_values: np.ndarray, gap_tolerance: float) -> list[list[float]]:
    """Single-linkage clusters of sorted, adjacent span values.

    A new cluster starts when the gap to the previous value is at least
    ``gap_tolerance``, so a chain of values each within ``gap_tolerance`` of its
    neighbour merges into one physical station.  Exact duplicates always merge,
    even at ``gap_tolerance == 0`` (a planar section, where every ring node has
    the same span coordinate).
    """
    if sorted_values.size == 0:
        return []
    clusters: list[list[float]] = [[float(sorted_values[0])]]
    for value in sorted_values[1:]:
        gap = float(value) - clusters[-1][-1]
        if gap < gap_tolerance or gap == 0.0:
            clusters[-1].append(float(value))
        else:
            clusters.append([float(value)])
    return clusters


def merge_span_stations(span_coords: np.ndarray, gap_tolerance: float) -> list[float]:
    """Merge adjacent span coordinates into physical section stations, as means.

    ``span_coords`` is a flat set of scalar span positions (one per section
    node).  They are sorted and clustered by single linkage: consecutive values
    closer than ``gap_tolerance`` belong to the same physical station, and each
    station is reported as the mean of its members.  The result is in increasing
    span order.  On this repo's prebent IEA-15MW blade one physical ring spreads
    over ~1e-3 m of span while neighbouring BEM stations sit ~0.44 m apart, so a
    scale-relative ``gap_tolerance`` separates rings without splitting one.

    Deterministic and pure: the input is never mutated.
    """
    values = np.asarray(span_coords, dtype=float).ravel()
    if values.size == 0:
        return []
    clusters = _clustered_spans(np.sort(values), gap_tolerance)
    return [float(np.mean(cluster)) for cluster in clusters]


def rings_in_band(span_coords: np.ndarray, gap_tolerance: float) -> list[np.ndarray]:
    """Group section-node indices into physical rings by their span coordinate.

    Uses the same single-linkage rule as :func:`merge_span_stations`.  Returns
    one int index array per merged station, in increasing span order, and every
    input index appears in exactly one group.  The indices index into
    ``span_coords`` itself (for a projector strip, into ``strip.offsets``); no
    mesh connectivity is used or needed.

    Deterministic and pure: the input is never mutated.
    """
    values = np.asarray(span_coords, dtype=float).ravel()
    if values.size == 0:
        return []
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    clusters = _clustered_spans(sorted_values, gap_tolerance)

    groups: list[np.ndarray] = []
    start = 0
    for cluster in clusters:
        groups.append(order[start : start + len(cluster)].astype(np.intp))
        start += len(cluster)
    return groups


@dataclass(frozen=True)
class ContourReport:
    """Outcome of ordering and validating one section ring."""

    order: np.ndarray
    area: float
    orientation: int
    perimeter: float
    closed: bool
    simple: bool
    usable: bool
    reason: str
    witness: tuple[int, int] | None


def _perimeter(ordered: np.ndarray) -> float:
    """Sum of Euclidean edge lengths around the closed ordered ring."""
    if ordered.shape[0] < 2:
        return 0.0
    edges = ordered - np.roll(ordered, -1, axis=0)
    return float(np.sum(np.linalg.norm(edges, axis=1)))


def _is_closed(ordered: np.ndarray) -> bool:
    """At least three points and no two consecutive points coincident (3-D exact)."""
    n = ordered.shape[0]
    if n < 3:
        return False
    for i in range(n):
        if np.array_equal(ordered[i], ordered[(i + 1) % n]):
            return False
    return True


def contour_report(points: np.ndarray, span_dir: np.ndarray) -> ContourReport:
    """Order ``points`` and report whether the ordered ring can carry ``q = M / (2A)``.

    ``order`` is :func:`order_ring`; ``area``/``orientation``/``perimeter`` come
    from the ordered ring; ``closed`` and ``simple`` are the exact predicates
    above applied to that ring.  ``usable`` is ``closed and simple and
    abs(area) > 0``.  ``reason`` names the first failing condition in plain words
    (``closed``, ``simple`` or ``area``), or ``"ok"`` when usable; ``witness``
    carries the offending edge pair when the ring is not simple.
    """
    pts = _as_points(points)
    order = order_ring(pts, span_dir)
    ordered = pts[order]

    area = signed_area(ordered, span_dir)
    orientation = 1 if area > 0.0 else (-1 if area < 0.0 else 0)
    perimeter = _perimeter(ordered)
    closed = _is_closed(ordered)
    simple, witness = is_simple(ordered, span_dir)

    if not closed:
        return ContourReport(
            order,
            area,
            orientation,
            perimeter,
            closed,
            simple,
            False,
            "ring is not closed: fewer than three points or coincident consecutive points",
            witness,
        )
    if not simple:
        return ContourReport(
            order,
            area,
            orientation,
            perimeter,
            closed,
            simple,
            False,
            "ring is not simple: intersecting edges or a zero-length wall edge",
            witness,
        )
    if not abs(area) > 0.0:
        return ContourReport(
            order,
            area,
            orientation,
            perimeter,
            closed,
            simple,
            False,
            "ring encloses zero area",
            witness,
        )
    return ContourReport(order, area, orientation, perimeter, closed, simple, True, "ok", witness)


@dataclass(frozen=True)
class SectionCell:
    """One bounded face (cell) of a thin-walled section's planar wall graph.

    A thin-walled section bounded by an outer skin plus internal webs splits
    into several closed cells; each one carries its own Bredt shear flow in the
    multi-cell torsional problem.  ``boundary`` is the cell's node indices in
    traversal order, a closed loop whose first index is **not** repeated at the
    end and is the smallest index in the loop, so the dataclass is stable and
    independent of the input edge order.  ``area`` is the positive shoelace area
    of that loop in the section frame (see :func:`section_plane_axes`): every
    bounded face is stored counter-clockwise as seen from ``+span``.
    ``edges`` lists the loop's boundary edges as ordered ``(a, b)`` pairs in the
    same traversal direction (``b = (a + 1)`` around the cell).
    """

    boundary: np.ndarray
    area: float
    edges: tuple[tuple[int, int], ...]


def _wall_edges(edges, n_points: int) -> set[tuple[int, int]]:
    """Normalise the wall graph to a set of undirected ``(min, max)`` pairs.

    Self-loops are dropped (a wall always joins two distinct nodes) and repeated
    edges are collapsed.  Out-of-range indices raise ``ValueError``.
    """
    normalised: set[tuple[int, int]] = set()
    for pair in edges:
        a, b = pair
        a = int(a)
        b = int(b)
        if a == b:
            continue
        if not (0 <= a < n_points) or not (0 <= b < n_points):
            raise ValueError(f"edge ({a}, {b}) is out of range for {n_points} points")
        normalised.add((a, b) if a < b else (b, a))
    return normalised


def _clockwise_successor(
    plane: np.ndarray, wall_edges: set[tuple[int, int]], n_points: int
) -> dict[tuple[int, int], int]:
    """Map each directed wall edge ``u -> v`` to the next edge ``v -> w``.

    Around every node the incident neighbours are sorted by increasing polar
    angle in the section plane; neighbours on the same ray (not angularly
    distinguishable) are ordered farthest-first, then by node index.  That
    deterministic tie-break treats a collinear wall as perturbed slightly
    clockwise, which keeps a dangling stub on the same ray as an existing wall
    in the outer face instead of cutting into a cell.  The successor of
    ``u -> v`` is the neighbour immediately *clockwise* from ``u`` about ``v`` --
    the previous element of the counter-clockwise angular order.  This keeps the
    cell interior on the left of each directed edge, so every bounded face is
    traversed counter-clockwise and has positive shoelace area.  A node whose
    incident edges are collinear is not an error: the tie-break above just picks
    one of the possible planar embeddings of the degenerate drawing.
    """
    neighbours: list[list[int]] = [[] for _ in range(n_points)]
    for a, b in wall_edges:
        neighbours[a].append(b)
        neighbours[b].append(a)

    successor: dict[tuple[int, int], int] = {}
    for node, incident in enumerate(neighbours):
        if not incident:
            continue
        rel = plane[incident] - plane[node]
        angle = np.arctan2(rel[:, 1], rel[:, 0])
        radius = np.hypot(rel[:, 0], rel[:, 1])
        # Primary key: angle.  Ties (same ray): farthest first (``-radius``),
        # then the node index, so the embedding is deterministic and collinear
        # stubs stay in the outer face.
        order = np.lexsort((np.asarray(incident, dtype=np.intp), -radius, angle))
        ring = [incident[i] for i in order]
        for i, incoming in enumerate(ring):
            successor[(incoming, node)] = ring[i - 1]
    return successor


def _face_loops(
    wall_edges: set[tuple[int, int]], successor: dict[tuple[int, int], int]
) -> list[list[int]]:
    """Traverse the planar embedding, returning one node loop per face.

    Every directed edge belongs to exactly one face walk.  The order in which
    walks are started is made deterministic (sorted undirected edges, each with
    both directions); the final cell order does not depend on it.
    """
    directed: list[tuple[int, int]] = []
    for a, b in sorted(wall_edges):
        directed.append((a, b))
        directed.append((b, a))

    visited: set[tuple[int, int]] = set()
    loops: list[list[int]] = []
    for start in directed:
        if start in visited:
            continue
        loop: list[int] = []
        current = start
        while current not in visited:
            visited.add(current)
            loop.append(current[0])
            current = (current[1], successor[current])
        loops.append(loop)
    return loops


def section_cells(points: np.ndarray, edges, span_dir: np.ndarray) -> list[SectionCell]:
    """Extract the bounded faces (cells) of a section's planar wall graph.

    ``points`` is an ``(n, 3)`` node array and ``edges`` a sequence of undirected
    index pairs: the in-plane wall graph of the section, made of skin segments
    and web segments.  ``span_dir`` is the span axis; the section frame comes
    from :func:`section_plane_axes`, so ``+span`` counter-clockwise is positive.

    The graph is embedded in that frame and its faces are traversed with the
    successor rule of :func:`_clockwise_successor`: from ``u -> v`` the next edge
    is ``v -> w`` with ``w`` the neighbour immediately clockwise from ``u`` about
    ``v``.  That keeps the interior on the left, so **bounded faces are
    counter-clockwise (positive shoelace) and only they are returned**; the
    single unbounded outer face of a connected graph is clockwise and is dropped
    by the ``area > 0`` filter (a graph with no closed loop yields only
    non-positive faces and returns ``[]``).

    Every returned cell's ``boundary`` starts at its lowest node index and runs
    counter-clockwise, and cells are ordered by that index, so the result is
    deterministic and independent of the input edge order.  Dangling walls (a
    stub that leads nowhere) are traversed into and back out of as part of the
    dropped outer face and never appear in a returned cell.

    The graph is assumed planar and drawn without crossing walls; a crossing
    graph is **not** detected here (no full geometric validation is performed).
    Degenerate input (no edges, fewer than three nodes, no closed loop) returns
    ``[]``.  Pure: ``points`` and ``edges`` are never mutated.
    """
    pts = _as_points(points)
    n_points = pts.shape[0]
    wall_edges = _wall_edges(edges, n_points)
    if n_points < 3 or not wall_edges:
        return []

    plane = _in_plane(pts, span_dir)
    successor = _clockwise_successor(plane, wall_edges, n_points)

    cells: list[SectionCell] = []
    for loop in _face_loops(wall_edges, successor):
        boundary = np.asarray(loop, dtype=np.intp)
        cell_plane = plane[boundary]
        x = cell_plane[:, 0]
        y = cell_plane[:, 1]
        area = 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
        if not area > 0.0:
            continue
        # Normalise the loop to start at its lowest node index.  All kept cells
        # are already counter-clockwise, so the direction needs no reversal.
        boundary = np.roll(boundary, -int(np.argmin(boundary)))
        boundary_edges = tuple(
            (int(boundary[i]), int(boundary[(i + 1) % boundary.shape[0]]))
            for i in range(boundary.shape[0])
        )
        cells.append(SectionCell(boundary, area, boundary_edges))

    cells.sort(key=lambda cell: int(cell.boundary[0]))
    return cells


def cell_adjacency(cells: list[SectionCell]) -> dict[tuple[int, int], list[int]]:
    """Map each bounded wall edge to the cell indices that use it.

    The key is the undirected edge ``(min, max)`` and the value lists the indices
    of the cells whose boundary contains that edge, in :func:`section_cells`
    output order.  A skin edge bounds one cell, a web shared by two cells appears
    with two owners; this is exactly the incidence the multi-cell Bredt shear
    flow system needs.  Edges that bound no cell (dangling stubs) are absent.
    Deterministic and pure.
    """
    adjacency: dict[tuple[int, int], list[int]] = {}
    for cell_index, cell in enumerate(cells):
        for a, b in cell.edges:
            key = (a, b) if a < b else (b, a)
            adjacency.setdefault(key, []).append(cell_index)
    return adjacency
