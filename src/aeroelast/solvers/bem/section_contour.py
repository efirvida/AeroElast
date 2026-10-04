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
    "section_plane_axes",
    "order_ring",
    "signed_area",
    "is_simple",
    "contour_report",
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
