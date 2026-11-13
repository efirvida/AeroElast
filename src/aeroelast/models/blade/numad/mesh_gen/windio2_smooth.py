"""Second-phase smoother for windIO v2 raw skin/web meshes.

Works in the native chart ``(eta, nd_arc)`` with a per-node physical metric
``G = J^T J`` (``J = dX/d(eta, nd)``), so element quality is measured in metres
without re-evaluating the loft inside the iteration loop.  Only nodes of bad
elements and their one-ring neighbours move; fixed nodes never move, sliding
nodes move along one chart line, and every node stays inside its ``bounds`` box
so material-cell ownership cannot change.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

__all__ = ["node_jacobians", "element_quality", "smooth_chart"]


def node_jacobians(evaluate: Callable[[np.ndarray], np.ndarray],
                   parameters: np.ndarray, step: float = 1e-6) -> np.ndarray:
    """Central finite-difference Jacobians ``dX/d(eta, nd)``, shape ``(N, 3, 2)``.

    ``evaluate(parameters (N,2)) -> xyz (N,3)``.
    """
    p = np.asarray(parameters, dtype=float)
    J = np.empty((len(p), 3, 2))
    for k in range(2):
        dp = np.zeros(2)
        dp[k] = step
        J[:, :, k] = (np.asarray(evaluate(p + dp), float)
                      - np.asarray(evaluate(p - dp), float)) / (2.0 * step)
    return J


def _node_metric(J: np.ndarray) -> np.ndarray:
    return np.einsum("nji,njk->nik", J, J)


def _corners(elements):
    el = np.asarray(elements)
    valid = el >= 0
    n = valid.sum(axis=1)
    k = np.arange(4)[None, :]
    nxt = np.where(k + 1 >= n[:, None], 0, k + 1)
    prv = np.where(k == 0, n[:, None] - 1, k - 1)
    rows = np.arange(len(el))[:, None]
    cur = np.where(valid, el, 0)
    return valid, cur, el[rows, nxt], el[rows, prv]


def _element_metric(G, el, valid):
    idx = np.where(valid, el, 0)
    g = (G[idx] * valid[:, :, None, None]).sum(axis=1)
    return g / valid.sum(axis=1)[:, None, None]


def _quality(par, el, G, corners=None):
    valid, cur, nxt, prv = corners if corners is not None else _corners(el)
    Ge = _element_metric(G, el, valid)
    a = par[nxt] - par[cur]
    b = par[prv] - par[cur]
    gaa = np.einsum("mki,mij,mkj->mk", a, Ge, a)
    gbb = np.einsum("mki,mij,mkj->mk", b, Ge, b)
    gab = np.einsum("mki,mij,mkj->mk", a, Ge, b)
    cross = np.sqrt(np.maximum(gaa * gbb - gab ** 2, 0.0))
    ang = np.where(valid, np.arctan2(cross, gab), np.inf)
    ln = np.where(valid, np.sqrt(np.maximum(gaa, 0.0)), np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        ar = np.nanmax(ln, axis=1) / np.nanmin(ln, axis=1)
    ar = np.where(np.isfinite(ar), ar, np.inf)
    area = a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]  # chart signed corner area
    return {"min_angle": ang.min(axis=1), "aspect": ar,
            "corner_area": np.where(valid, area, 0.0)}


def element_quality(parameters, elements, jacobians) -> dict:
    """Per-element ``min_angle`` (rad) and ``aspect`` ratio measured in the node metric."""
    G = _node_metric(np.asarray(jacobians, float))
    q = _quality(np.asarray(parameters, float), np.asarray(elements), G)
    return {"min_angle": q["min_angle"], "aspect": q["aspect"]}


def _edges(el):
    valid, cur, nxt, _ = _corners(el)
    i = cur[valid]
    j = nxt[valid]
    return np.concatenate([i, j]), np.concatenate([j, i])


def smooth_chart(parameters, elements, jacobians_or_evaluator, fixed, slide_lines=None,
                 iterations: int = 30, bounds=None, min_angle_deg: float = 20.0,
                 max_aspect: float = 5.0, relax: float = 1.0, max_backtrack: int = 8):
    """Smooth bad elements in chart space; returns ``(new_parameters, report)``.

    Parameters
    ----------
    parameters : (N, 2) chart coordinates.
    elements : (M, 4) int, triangles padded with -1.
    jacobians_or_evaluator : (N, 3, 2) array, or ``evaluate(parameters) -> xyz``
        (Jacobians then come from :func:`node_jacobians`).
    fixed : (N,) bool mask, or index array; these rows are returned bit-identical.
    slide_lines : ``{node: direction (2,)}``; the node moves only along that chart
        direction through its initial position.
    bounds : (N, 4) ``(eta_lo, eta_hi, nd_lo, nd_hi)`` per node (optional).

    A step is accepted only if no element's chart orientation flips and neither
    the worst min-angle gets smaller nor the worst aspect ratio larger
    (step halved up to ``max_backtrack`` times).  ``report`` keys: ``moved``,
    ``min_angle_before/after`` (deg), ``max_aspect_before/after``,
    ``bad_before/after``, ``iterations``.
    """
    par0 = np.array(parameters, dtype=float)
    el = np.asarray(elements)
    n = len(par0)
    J = (np.asarray(jacobians_or_evaluator, float)
         if not callable(jacobians_or_evaluator)
         else node_jacobians(jacobians_or_evaluator, par0))
    G = _node_metric(J)
    fx = np.zeros(n, bool)
    fixed = np.asarray(fixed)
    if fixed.dtype == bool:
        fx[:] = fixed
    else:
        fx[fixed] = True

    corners = _corners(el)
    lim = np.radians(min_angle_deg)

    def evaluate(p):
        q = _quality(p, el, G, corners)
        return q["min_angle"], q["aspect"], q["corner_area"]

    ang, ar, area = evaluate(par0)
    sgn = np.sign(area.sum(axis=1))
    sgn[sgn == 0] = 1.0
    ok0 = area * sgn[:, None] > 0
    bad_el = (ang < lim) | (ar > max_aspect)

    def report(par, iters):
        a, r, _ = evaluate(par)
        return {"moved": int((par != par0).any(axis=1).sum()),
                "min_angle_before": float(np.degrees(ang.min())) if len(ang) else np.nan,
                "min_angle_after": float(np.degrees(a.min())) if len(a) else np.nan,
                "max_aspect_before": float(ar.max()) if len(ar) else np.nan,
                "max_aspect_after": float(r.max()) if len(r) else np.nan,
                "bad_before": int(bad_el.sum()),
                "bad_after": int(((a < lim) | (r > max_aspect)).sum()),
                "iterations": iters}

    if not bad_el.any():
        return par0.copy(), report(par0, 0)

    ei, ej = _edges(el)
    mark = np.zeros(n, bool)
    mark[el[bad_el][el[bad_el] >= 0]] = True
    ring = mark.copy()
    ring[ej[mark[ei]]] = True
    movable = ring & ~fx
    mv = np.flatnonzero(movable)
    if len(mv) == 0:
        return par0.copy(), report(par0, 0)

    # edge list restricted to edges ending at movable nodes (i = movable, j = neighbour)
    sel = movable[ei]
    ei, ej = ei[sel], ej[sel]
    Gij = 0.5 * (G[ei] + G[ej])

    slide = np.zeros(n, bool)
    sdir = np.zeros((n, 2))
    for node, d in (slide_lines or {}).items():
        d = np.asarray(d, float)
        sdir[node] = d / np.linalg.norm(d)
        slide[node] = True
    slide &= movable

    if bounds is None:
        lo = np.full((n, 2), -np.inf)
        hi = np.full((n, 2), np.inf)
    else:
        b = np.asarray(bounds, float)
        lo = np.minimum(b[:, [0, 2]], par0)
        hi = np.maximum(b[:, [1, 3]], par0)

    # parameter range of t along each sliding line that stays in the bounds box
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = (lo - par0) / sdir
        t2 = (hi - par0) / sdir
    tlo = np.where(sdir != 0, np.minimum(t1, t2), -np.inf).max(axis=1)
    thi = np.where(sdir != 0, np.maximum(t1, t2), np.inf).min(axis=1)
    tlo, thi = np.minimum(tlo, 0.0), np.maximum(thi, 0.0)

    par = par0.copy()
    done = 0
    for _it in range(iterations):
        cur_ang, cur_ar, _ = evaluate(par)
        if not ((cur_ang < lim) | (cur_ar > max_aspect)).any():
            break
        A = np.zeros((n, 2, 2))
        for r in range(2):
            for c in range(2):
                A[:, r, c] = np.bincount(ei, Gij[:, r, c], minlength=n)
        rhs = np.stack([np.bincount(ei, np.einsum("ek,ek->e", Gij[:, r, :], par[ej]),
                                    minlength=n) for r in range(2)], axis=1)
        Am, rm = A[mv], rhs[mv]
        det = Am[:, 0, 0] * Am[:, 1, 1] - Am[:, 0, 1] * Am[:, 1, 0]
        good = np.abs(det) > 1e-300
        tgt = par[mv].copy()
        inv = np.stack([np.stack([Am[:, 1, 1], -Am[:, 0, 1]], -1),
                        np.stack([-Am[:, 1, 0], Am[:, 0, 0]], -1)], axis=1)
        sol = np.einsum("nij,nj->ni", inv, rm) / np.where(good, det, 1.0)[:, None]
        tgt[good] = sol[good]
        delta = np.zeros((n, 2))
        delta[mv] = tgt - par[mv]
        # sliding: project the (metric-weighted) step onto the line
        sm = slide[mv]
        t = np.zeros(n)
        if sm.any():
            s = mv[sm]
            Ad = np.einsum("nij,nj->ni", A[s], sdir[s])
            t[s] = np.einsum("ni,ni->n", Ad, delta[s]) / np.einsum("ni,ni->n", Ad, sdir[s])
        accepted = False
        alpha = relax
        for _ in range(max_backtrack):
            trial = par.copy()
            free = movable & ~slide
            trial[free] = np.clip(par[free] + alpha * delta[free], lo[free], hi[free])
            sl = np.flatnonzero(slide)
            if len(sl):
                tt = np.clip(alpha * t[sl], tlo[sl], thi[sl])
                trial[sl] = par[sl] + tt[:, None] * sdir[sl]
            a, r, area_t = evaluate(trial)
            flip = ((area_t * sgn[:, None] <= 0) & ok0).any()
            if (not flip and a.min() >= cur_ang.min() - 1e-14
                    and r.max() <= cur_ar.max() + 1e-12
                    and (a.min() > cur_ang.min() or r.max() < cur_ar.max())):
                par, accepted = trial, True
                break
            alpha *= 0.5
        if not accepted:
            break
        done += 1
    return par, report(par, done)
