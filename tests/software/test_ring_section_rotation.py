"""Arc-length-weighted section rotation of one ring (l3-closure T19, closed forms).

The section rotation is the antisymmetric part of the ring's in-plane affine fit.  Weighted by
each node's tributary wall length it converges to a wall-length integral, so it does not depend on
where the mesher puts the nodes; with equal nodal weights it follows the node density.  The
closed form is a circle of radius ``R`` carrying the tangential field
``u_t = R theta0 (1 + a cos s)``: its arc-weighted affine fit is exactly ``theta0`` (the
``cos s`` term projects onto ``sin 2s``, which the affine fit cannot see), while nodes clustered
at ``s = 0`` pull an equal-weight fit toward ``theta0 (1 + a)``.
"""

from __future__ import annotations

import numpy as np
import pytest

from aeroelast.solvers.bem.section_contour import ring_section_rotation, ring_tributary_lengths

SPAN = np.array([0.0, 0.0, 1.0])
R = 2.0
THETA0 = 1.0e-3
A = 0.5


def _graded_circle(n: int, grading: float):
    """``n`` nodes on a circle, spacing ratio ~(1 + g)/(1 - g) with the dense end at ``s = 0``."""
    t = np.arange(n) / n
    s = 2.0 * np.pi * t - grading * np.sin(2.0 * np.pi * t)
    points = np.column_stack([R * np.cos(s), R * np.sin(s), np.full(n, 30.0)])
    edges = np.column_stack([np.arange(n), (np.arange(n) + 1) % n])
    return s, points, edges


def _tangential_field(s: np.ndarray, a: float) -> np.ndarray:
    u_t = R * THETA0 * (1.0 + a * np.cos(s))
    return np.column_stack([-u_t * np.sin(s), u_t * np.cos(s), np.zeros_like(s)])


def test_tributary_lengths_of_a_closed_ring_with_a_web():
    """Unit square (side 1) with a web from node 0 to node 2: half of every incident wall."""
    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0], [0, 2]])
    w = ring_tributary_lengths(np.arange(4), edges, points)
    web = np.sqrt(2.0)
    np.testing.assert_allclose(w, [1.0 + web / 2, 1.0, 1.0 + web / 2, 1.0], rtol=0, atol=1e-15)
    assert w.sum() == pytest.approx(4.0 + web, abs=1e-14)


def test_tributary_lengths_ignore_edges_leaving_the_ring():
    """Global indices; an edge with one end outside the ring (a spanwise edge) carries no wall."""
    points = np.array([[0.0, 0, 0], [9.0, 9, 9], [1.0, 0, 0], [1.0, 1, 0], [0.0, 1, 0]])
    ring = np.array([0, 2, 3, 4])
    edges = np.array([[0, 2], [2, 3], [3, 4], [4, 0], [0, 1], [1, 3]])
    np.testing.assert_allclose(ring_tributary_lengths(ring, edges, points), [1.0, 1.0, 1.0, 1.0])


@pytest.mark.parametrize("weighted", [True, False])
def test_rigid_motion_is_exact_with_or_without_weights(weighted):
    s, points, edges = _graded_circle(40, 0.6)
    theta = 2.5e-3
    disp = theta * np.cross(np.broadcast_to(SPAN, points.shape), points) + [0.3, -0.2, 0.1]
    weights = ring_tributary_lengths(np.arange(len(s)), edges, points) if weighted else None
    assert ring_section_rotation(points, disp, SPAN, weights) == pytest.approx(theta, rel=1e-12)


def test_arc_weighted_rotation_does_not_depend_on_the_node_layout():
    """P1: theta0 within 1e-3 at n >= 48 on two layouts; equal weights miss it by > 5 %."""
    results = {}
    for n, grading in ((48, 0.6), (96, 0.3)):
        s, points, edges = _graded_circle(n, grading)
        weights = ring_tributary_lengths(np.arange(n), edges, points)
        disp = _tangential_field(s, A)
        results[n] = (ring_section_rotation(points, disp, SPAN, weights),
                      ring_section_rotation(points, disp, SPAN))
    for n, (arc, equal) in results.items():
        print(f"n={n}: arc {arc / THETA0:.6f} theta0, equal weights {equal / THETA0:.6f} theta0")
        assert arc == pytest.approx(THETA0, rel=1e-3)
    assert results[48][0] == pytest.approx(results[96][0], rel=1e-3)
    assert abs(results[48][1] / THETA0 - 1.0) > 0.05


def test_degenerate_rings_have_no_rotation():
    """Fewer than three weighted nodes, or collinear ones, define no 2-D section."""
    points = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.0, 0, 0], [0.0, 1, 0]])
    disp = np.zeros_like(points)
    assert ring_section_rotation(points[:2], disp[:2], SPAN) is None
    assert ring_section_rotation(points[:3], disp[:3], SPAN) is None
    assert ring_section_rotation(points, disp, SPAN, np.array([1.0, 1.0, 1.0, 0.0])) is None
