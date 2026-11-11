"""Second-phase chart smoother for windIO v2 raw meshes (synthetic evaluators)."""

import numpy as np
import pytest

from aeroelast.models.blade.numad.mesh_gen.windio2_smooth import (
    element_quality,
    node_jacobians,
    smooth_chart,
)

N = 11


def flat(p):
    p = np.asarray(p, float)
    return np.column_stack([2.0 * p[:, 0], p[:, 1], 0.0 * p[:, 0]])


def cylinder(p):
    p = np.asarray(p, float)
    r, ln = 0.5, 3.0
    th = 2.0 * np.pi * 0.25 * p[:, 1] / (N - 1)
    return np.column_stack([r * np.cos(th), r * np.sin(th), ln * p[:, 0] / (N - 1)])


def grid(n=N):
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    par = np.column_stack([ii.ravel(), jj.ravel()]).astype(float)
    idx = lambda i, j: i * n + j
    el = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)]
          for i in range(n - 1) for j in range(n - 1)]
    return par, np.array(el), np.isin(par, [0, n - 1]).any(axis=1)


def to_tris(el):
    out = []
    for q in el:
        out.append([q[0], q[1], q[2], -1])
        out.append([q[0], q[2], q[3], -1])
    return np.array(out)


def skew(par, node=5 * N + 5, d=(0.9, 0.0)):
    p = par.copy()
    p[node] += d
    return p


def bounds_for(par):
    b = np.empty((len(par), 4))
    b[:, 0], b[:, 1], b[:, 2], b[:, 3] = 0, N - 1, 0, N - 1
    return b


def signed_areas(par, el):
    q = el[:, :4]
    x, y = par[:, 0], par[:, 1]
    a = np.zeros(len(el))
    n = (q >= 0).sum(1)
    for k in range(4):
        i = q[:, k]
        j = np.where(k + 1 >= n, q[:, 0], q[:, (k + 1) % 4])
        ok = k < n
        a += np.where(ok, x[i] * y[j] - x[j] * y[i], 0.0)
    return a


@pytest.mark.parametrize("ev", [flat, cylinder])
@pytest.mark.parametrize("tri", [False, True])
def test_smooth_improves_skewed_patch(ev, tri):
    par, el, fixed = grid()
    if tri:
        el = to_tris(el)
    bad = skew(par)
    J = node_jacobians(ev, bad)
    b = bounds_for(bad)
    new, rep = smooth_chart(bad, el, J, fixed, bounds=b)
    assert rep["min_angle_after"] > rep["min_angle_before"] + 1.0
    assert rep["max_aspect_after"] <= rep["max_aspect_before"]
    assert rep["moved"] > 0 and rep["iterations"] >= 1
    # fixed bit-identical
    assert np.array_equal(new[fixed], bad[fixed])
    # bounds
    assert (new[:, 0] >= b[:, 0]).all() and (new[:, 0] <= b[:, 1]).all()
    assert (new[:, 1] >= b[:, 2]).all() and (new[:, 1] <= b[:, 3]).all()
    # no inversion
    assert (signed_areas(new, el) * signed_areas(bad, el) > 0).all()
    # independent recomputation of the reported quality
    q = element_quality(new, el, J)
    assert np.isclose(np.degrees(q["min_angle"]).min(), rep["min_angle_after"])


def test_only_bad_neighbourhood_moves():
    par, el, fixed = grid()
    bad = skew(par)
    new, _ = smooth_chart(bad, el, node_jacobians(flat, bad), fixed)
    moved = np.flatnonzero((new != bad).any(1))
    assert 0 < len(moved) < 60
    far = np.flatnonzero(np.abs(par - 5).max(1) > 3)
    assert not np.isin(far, moved).any()


def test_good_mesh_unchanged():
    par, el, fixed = grid()
    new, rep = smooth_chart(par, el, node_jacobians(flat, par), fixed)
    assert np.array_equal(new, par)
    assert rep["moved"] == 0 and rep["iterations"] == 0


def test_sliding_node_stays_on_line():
    par, el, fixed = grid()
    node = 5 * N + 5
    bad = skew(par, node, (0.9, 0.0))
    d = np.array([1.0, 0.0])
    J = node_jacobians(flat, bad)
    new, rep = smooth_chart(bad, el, J, fixed, slide_lines={node: d})
    off = new[node] - bad[node]
    assert off[1] == 0.0 or abs(off[1]) < 1e-12
    assert abs(off[0]) > 0.1
    assert rep["max_aspect_after"] < rep["max_aspect_before"]


def test_sliding_oblique_line():
    par, el, fixed = grid()
    node = 5 * N + 5
    bad = skew(par, node, (0.6, 0.6))
    d = np.array([1.0, 1.0]) / np.sqrt(2)
    new, _ = smooth_chart(bad, el, node_jacobians(flat, bad), fixed, slide_lines={node: d})
    off = new[node] - bad[node]
    assert abs(off[0] * d[1] - off[1] * d[0]) < 1e-12


def test_sliding_respects_bounds():
    par, el, fixed = grid()
    node = 5 * N + 5
    bad = skew(par, node, (0.9, 0.0))
    b = bounds_for(bad)
    b[node] = (5.5, 6.0, 0, N - 1)
    new, _ = smooth_chart(bad, el, node_jacobians(flat, bad), fixed,
                          slide_lines={node: np.array([1.0, 0.0])}, bounds=b)
    assert new[node, 0] >= 5.5 - 1e-15


def test_node_jacobians_cylinder():
    p = np.array([[3.0, 4.0], [7.0, 1.0]])
    J = node_jacobians(cylinder, p)
    assert J.shape == (2, 3, 2)
    assert np.allclose(np.linalg.norm(J[:, :, 0], axis=1), 3.0 / (N - 1), rtol=1e-6)
    assert np.allclose(np.linalg.norm(J[:, :, 1], axis=1), 0.5 * 2 * np.pi * 0.25 / (N - 1), rtol=1e-6)


def test_callable_evaluator_accepted():
    par, el, fixed = grid()
    bad = skew(par)
    a, _ = smooth_chart(bad, el, flat, fixed)
    b, _ = smooth_chart(bad, el, node_jacobians(flat, bad), fixed)
    assert np.allclose(a, b, atol=1e-9)
