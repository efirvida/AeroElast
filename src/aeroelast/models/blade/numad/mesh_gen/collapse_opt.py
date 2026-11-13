"""Shape optimizer for the thin-strip collapse of windIO v2 raw shell meshes.

``windio2_collapse.collapse_thin_strips`` merges the ends of short cross edges.  Where a thin band
opens up again the last quad keeps one merged edge and becomes a triangle with sides ``(w, l, l)``
(``w`` = band width, ``l`` = element length): an isosceles sliver of min angle ``asin(w / 2l)``.
It is equilateral only if the merged apex sits ``sqrt(3)/2 * w`` from the base, i.e. the apex has to
*slide along its merged line towards the base*; the next interval of that line then grows and the
quads behind it get an aspect ratio up to about 2.

This module is deliberately outside the ``windio2_*`` glob (that glob feeds the mesh cache key of
``tests/support/v2_blade.py``) and only *reads* ``windio2_collapse``.  It provides

* ``collapse_capped``     the same merge as ``collapse_thin_strips`` plus a bound ``max_shift`` on
                          the cumulative displacement of every merged node (identical result when
                          ``max_shift`` is infinite);
* ``optimize_collapsed``  local 1D search per node, a few sweeps, maximising the worst element quality
                          around the node; refuses any move that inverts an element, moves a fixed
                          node, breaks the slide / shift budgets, or leaves an incident quad below a
                          quality floor (or worse than it already was);
* ``collapse_and_optimize_raw_mesh``  the ``RawMesh`` wrapper (same remapping as ``collapse_raw_mesh``).

Quality (``element_quality``) is 1 for the ideal element and falls to 0:
  triangle  min_angle / 60 deg          (equilateral = 1)
  quad      min(min_angle / 60, 2 / AR, (180 - max_angle) / 60), clipped to 1
            (a rectangle up to AR 2 is "good": it is what the neighbours of a slid apex become).
Nothing here knows about a particular blade.
"""

from __future__ import annotations

import numpy as np

from .windio2_collapse import _normal, _pad

__all__ = ["collapse_capped", "cumulative_shift", "element_quality", "mesh_stats",
           "optimize_collapsed", "collapse_and_optimize_raw_mesh"]


# --------------------------------------------------------------------------------------------------
# quality
# --------------------------------------------------------------------------------------------------
def _angles_lengths(p):
    """``p`` (..., k, 3) -> (angles in deg (..., k) at each vertex, edge lengths (..., k))."""
    e = np.roll(p, -1, axis=-2) - p                       # edge j: v_j -> v_{j+1}
    ln = np.linalg.norm(e, axis=-1)
    prev = -np.roll(e, 1, axis=-2)                         # v_j -> v_{j-1}
    cosv = np.einsum("...i,...i->...", prev, e) / np.maximum(np.roll(ln, 1, axis=-1) * ln, 1e-300)
    return np.degrees(np.arccos(np.clip(cosv, -1.0, 1.0))), ln


def element_quality(p):
    """Quality in [0, 1] of polygons ``p`` (..., k, 3), k = 3 or 4 (see module docstring)."""
    ang, ln = _angles_lengths(p)
    if p.shape[-2] == 3:
        return np.clip(ang.min(axis=-1) / 60.0, 0.0, 1.0)
    ar = ln.max(axis=-1) / np.maximum(ln.min(axis=-1), 1e-300)
    q = np.minimum.reduce([ang.min(axis=-1) / 60.0, 2.0 / ar, (180.0 - ang.max(axis=-1)) / 60.0])
    return np.clip(q, 0.0, 1.0)


def _corner_dots(p, nref):
    """Per-corner cross(e_{j-1}, e_j) . nref / |nref| for polygons (..., k, 3); positive = not inverted."""
    e = np.roll(p, -1, axis=-2) - p
    cr = np.cross(np.roll(e, 1, axis=-2), e)
    return np.einsum("...i,i->...", cr, nref / np.linalg.norm(nref))


def mesh_stats(nodes, elements):
    """Shape statistics of a (padded) mesh: AR counts, triangle min-angle distribution, quad max angle."""
    x = np.asarray(nodes, float)
    el = _pad(elements)
    tri = (el < 0).any(axis=1)
    ar = np.empty(len(el))
    amin = np.empty(len(el))
    amax = np.empty(len(el))
    for sel, k in ((tri, 3), (~tri, 4)):
        if sel.any():
            ang, ln = _angles_lengths(x[el[sel][:, :k]])
            ar[sel] = ln.max(-1) / np.maximum(ln.min(-1), 1e-300)
            amin[sel] = ang.min(-1)
            amax[sel] = ang.max(-1)
    out = {"elements": int(len(el)), "triangles": int(tri.sum()), "ar_gt5": int((ar > 5).sum()),
           "ar_gt10": int((ar > 10).sum()), "max_ar": float(ar.max()) if len(el) else 0.0,
           "quad_ar_gt5": int(((ar > 5) & ~tri).sum()), "quad_max_angle_gt120": int(((amax > 120) & ~tri).sum()),
           "quad_min_angle_p5": float(np.percentile(amin[~tri], 5)) if (~tri).any() else float("nan")}
    if tri.any():
        a = amin[tri]
        p5, p25, p50 = np.percentile(a, [5, 25, 50])
        out.update(tri_min_angle_p5=float(p5), tri_min_angle_p25=float(p25), tri_min_angle_p50=float(p50),
                   tri_min_angle_lt20=int((a < 20).sum()), tri_min_angle_lt30=int((a < 30).sum()),
                   tri_min_angle_ge50=int((a >= 50).sum()))
    return out


# --------------------------------------------------------------------------------------------------
# capped collapse
# --------------------------------------------------------------------------------------------------
def cumulative_shift(nodes0, nodes_new, node_map):
    """Max over original nodes of |final position - original position| (the true cumulative shift)."""
    d = np.linalg.norm(np.asarray(nodes_new, float)[np.asarray(node_map)] - np.asarray(nodes0, float), axis=1)
    return float(d.max()) if len(d) else 0.0


def collapse_capped(nodes, elements, *, max_edge_len, fixed=(), min_aspect=3.0, max_shift=np.inf):
    """``collapse_thin_strips`` with a bound on the cumulative shift of every merged node.

    A merge is refused when its midpoint would lie farther than ``max_shift`` from any original node of
    either group (merges cascade, so the per-merge bound ``max_edge_len`` alone lets lines drift tens of
    mm).  With ``max_shift = inf`` the result equals ``collapse_thin_strips``.  Same return signature;
    the report adds ``max_cumulative_shift_m``, ``refused_by_cap``.
    """
    x0 = np.asarray(nodes, dtype=float)
    el0 = _pad(elements)
    pos = x0.copy()
    rep = np.arange(len(x0))
    group = {}
    is_fixed = np.zeros(len(x0), bool)
    is_fixed[list(fixed)] = True
    members = [el0[k][el0[k] >= 0].tolist() for k in range(len(el0))]

    def find(i):
        while rep[i] != i:
            rep[i] = rep[rep[i]]
            i = rep[i]
        return i

    def current(k):
        ids = [find(i) for i in members[k]]
        out = [i for j, i in enumerate(ids) if i != ids[j - 1]]
        if len(out) > 1 and out[0] == out[-1]:
            out.pop()
        return out

    n0 = [_normal(x0[m]) for m in members]
    area0 = np.array([np.linalg.norm(v) for v in n0])
    incident = [set() for _ in range(len(x0))]
    for k, m in enumerate(members):
        for i in m:
            incident[i].add(k)

    def aspect(k):
        p = x0[members[k]]
        ln = np.linalg.norm(np.roll(p, -1, 0) - p, axis=1)
        return ln.max() / max(ln.min(), 1e-300)

    cand = {}
    for k, m in enumerate(members):
        if aspect(k) < min_aspect:
            continue
        for j in range(len(m)):
            a, b = m[j], m[(j + 1) % len(m)]
            length = float(np.linalg.norm(x0[a] - x0[b]))
            if length < max_edge_len:
                cand[(min(a, b), max(a, b))] = length
    merged = refused = 0
    for (a, b), _ in sorted(cand.items(), key=lambda kv: kv[1]):
        ra, rb = find(a), find(b)
        if ra == rb or (is_fixed[ra] and is_fixed[rb]):
            continue
        keep, gone = (ra, rb) if (is_fixed[ra] or (not is_fixed[rb] and ra < rb)) else (rb, ra)
        new = pos[keep] if is_fixed[keep] else 0.5 * (pos[ra] + pos[rb])
        if np.linalg.norm(pos[gone] - new) > max_edge_len or np.linalg.norm(pos[keep] - new) > max_edge_len:
            continue
        gm = group.get(keep, [keep]) + group.get(gone, [gone])
        if np.isfinite(max_shift) and np.linalg.norm(x0[gm] - new, axis=1).max() > max_shift:
            refused += 1
            continue
        old_keep = pos[keep].copy()
        pos[keep] = new
        rep[gone] = keep
        ok = True
        for k in incident[keep] | incident[gone]:
            ids = current(k)
            if len(ids) >= 3 and float(np.dot(_normal(pos[ids]), n0[k])) <= 0.0:
                ok = False
                break
        if not ok:
            rep[gone] = gone
            pos[keep] = old_keep
            continue
        group[keep] = gm
        group.pop(gone, None)
        incident[keep] |= incident[gone]
        merged += 1

    out_el, source = [], []
    dropped = to_tri = 0
    removed = 0.0
    for k in range(len(el0)):
        ids = current(k)
        if len(ids) < 3:
            dropped += 1
            removed += area0[k]
            continue
        if len(ids) < len(members[k]):
            to_tri += 1
        out_el.append(ids + [-1] * (4 - len(ids)))
        source.append(k)
    if merged == 0:
        return x0, el0, np.arange(len(el0)), {"merged": 0, "dropped": 0, "to_triangle": 0, "removed_area_m2": 0.0,
                                              "max_node_shift_m": 0.0, "max_cumulative_shift_m": 0.0,
                                              "refused_by_cap": refused, "node_map": np.arange(len(x0))}
    out_el = np.array(out_el, dtype=int).reshape(-1, 4)
    used = np.unique(out_el[out_el >= 0])
    remap = -np.ones(len(x0), dtype=int)
    remap[used] = np.arange(len(used))
    out_el = np.where(out_el >= 0, remap[np.where(out_el >= 0, out_el, 0)], -1)
    nmap = np.array([remap[find(i)] for i in range(len(x0))], dtype=int)
    new_nodes = pos[used]
    report = {"merged": merged, "dropped": dropped, "to_triangle": to_tri, "removed_area_m2": float(removed),
              "nodes_before": len(x0), "nodes_after": len(used), "refused_by_cap": refused,
              "max_cumulative_shift_m": cumulative_shift(x0, new_nodes, nmap), "node_map": nmap}
    return new_nodes, out_el, np.array(source, dtype=int), report


# --------------------------------------------------------------------------------------------------
# optimizer
# --------------------------------------------------------------------------------------------------
def _span_field(apex, nb, hops):
    """Reference line direction for every node within ``hops`` edges of an apex (nearest apex wins)."""
    ref = dict(apex)
    frontier = list(apex)
    for _ in range(hops):
        nxt = []
        for i in frontier:
            for j in nb[i]:
                if j not in ref:
                    ref[j] = ref[i]
                    nxt.append(j)
        frontier = nxt
    return ref


def _directions(x, p, nbrs, apex_dir=None, ref=None, tol_deg=25.0):
    """Slide directions at node ``p`` that keep it on its line.

    A line through ``p`` is a pair of neighbours ``(a, b)`` seen at more than ~135 deg; the direction is
    ``b - a``.  Edge directions are NOT used: a cross edge would drag the node across the band and move
    a material boundary.  An apex (the vertex of a triangle opposite its shortest edge) is the end of a
    merged line: its direction is base midpoint -> apex, given as ``apex_dir``.  Neighbours of an apex
    (``ref`` = the apex direction, ``hops`` edges around it) may slide only along the line family that is
    within ``tol_deg`` of ``ref``: in a structured patch both families are straight chains and only the
    apex tells which one is the material line.
    """
    if apex_dir is not None:
        return [apex_dir]
    if ref is None:                                          # not near any merged line end: nothing to slide
        return []
    vec = [x[q] - x[p] for q in nbrs]
    vec = [v / np.linalg.norm(v) for v in vec if np.linalg.norm(v) > 0]
    out = []
    for i in range(len(vec)):
        for j in range(i + 1, len(vec)):
            if vec[i] @ vec[j] < -0.7:                      # angle at p > ~135 deg: a line through p
                d = vec[j] - vec[i]
                d = d / np.linalg.norm(d)
                if abs(d @ ref) >= np.cos(np.radians(tol_deg)) and all(abs(d @ o) < 0.9995 for o in out):
                    out.append(d)
    return out


def _apex_directions(x, el, cnt, sharp=0.9):
    """``{apex node: unit direction base-midpoint -> apex}`` for triangles whose shortest edge is < sharp * the others."""
    out = {}
    for k in np.nonzero(cnt == 3)[0]:
        ids = el[k, :3]
        p = x[ids]
        ln = np.linalg.norm(np.roll(p, -1, 0) - p, axis=1)      # edge j: v_j -> v_{j+1}; opposite vertex of edge j is j+2
        j = int(np.argmin(ln))
        others = np.delete(ln, j)
        if ln[j] >= sharp * others.min():
            continue
        apex = int(ids[(j + 2) % 3])
        mid = 0.5 * (p[j] + p[(j + 1) % 3])
        d = x[apex] - mid
        if np.linalg.norm(d) > 0:
            out.setdefault(apex, d / np.linalg.norm(d))
    return out


def optimize_collapsed(nodes, elements, *, fixed=(), nodes0=None, node_map=None, max_slide=0.15,
                       max_shift=np.inf, quad_floor=0.5, max_quad_angle=120.0, max_normal=0.005, target=0.95, sweeps=6,
                       n_grid=49, hops=2):
    """Slide nodes along their lines to make the worst incident element as good as possible.

    ``nodes`` (n, 3), ``elements`` (m, 4) padded.  ``fixed`` ids never move.  ``nodes0`` / ``node_map``
    (original positions and original->current node map, as returned by the collapse) define the
    *cumulative* shift of a node: the largest distance from its position to any original node merged
    into it; a move is refused when that would exceed ``max(max_shift, current shift)``.  ``max_slide``:
    largest distance a node may end from the position it had on entry.  ``quad_floor``: an incident quad
    may not end below ``min(quad_floor, its quality before the move)`` nor with a larger max angle than
    ``max(max_quad_angle, its max angle before)`` degrees.  ``max_normal``: a straight slide leaves the
    curved skin, so no node may end farther than this (m) from where it entered, measured along the mean
    normal of its incident elements at entry.  Only nodes of elements with
    quality below ``target`` and within ``hops`` edges of a triangle apex are visited; ``sweeps`` bounds the passes.

    Returns ``(new_nodes, report)``; the topology is unchanged.
    """
    x = np.array(nodes, dtype=float)
    xin = x.copy()
    el = _pad(elements)
    nn = len(x)
    is_fixed = np.zeros(nn, bool)
    is_fixed[list(fixed)] = True
    cnt = (el >= 0).sum(1)
    inc = [[] for _ in range(nn)]
    nb = [set() for _ in range(nn)]
    for k in range(len(el)):
        ids = el[k, :cnt[k]]
        for j, i in enumerate(ids):
            inc[i].append(k)
            nb[i].update((int(ids[j - 1]), int(ids[(j + 1) % cnt[k]])))
    if nodes0 is not None and node_map is not None:
        x0 = np.asarray(nodes0, float)
        grp = [[] for _ in range(nn)]
        for o, c in enumerate(np.asarray(node_map)):
            grp[c].append(o)
        grp = [x0[g] if g else x[i:i + 1] for i, g in enumerate(grp)]
    else:
        grp = [x[i:i + 1] for i in range(nn)]

    def shift_of(i, pts):                                   # pts (T, 3) -> (T,)
        return np.linalg.norm(grp[i][None, :, :] - pts[:, None, :], axis=2).max(axis=1)

    def qual(k, xs):
        return element_quality(xs[el[k, :cnt[k]]])

    q = np.array([qual(k, x) for k in range(len(el))])
    q_before = q.copy()
    moved_total = 0
    sweeps_done = 0
    for _ in range(sweeps):
        active = np.unique(el[(q < target)[:, None] & (el >= 0)])
        active = [int(i) for i in active if not is_fixed[i]]
        apex = _apex_directions(x, el, cnt)
        ref = _span_field(apex, nb, hops)
        active.sort(key=lambda i: (i not in apex, min(q[k] for k in inc[i])))      # apexes first
        changed = 0
        for p in active:
            ks = inc[p]
            nrm = [_normal(x[el[k, :cnt[k]]]) for k in ks]
            nrm = [v if np.linalg.norm(v) > 0 else np.array([0.0, 0.0, 1.0]) for v in nrm]
            cmin0 = [_corner_dots(x[el[k, :cnt[k]]], nrm[a]).min() for a, k in enumerate(ks)]
            q0 = np.array([q[k] for k in ks])
            floors = np.array([min(quad_floor, q[k]) if cnt[k] == 4 else -1.0 for k in ks])
            n_in = sum(_normal(xin[el[k, :cnt[k]]]) for k in ks)
            n_in = n_in / np.linalg.norm(n_in) if np.linalg.norm(n_in) > 0 else np.array([0.0, 0.0, 1.0])
            amax0 = [_angles_lengths(x[el[k, :4]])[0].max() if cnt[k] == 4 else 0.0 for k in ks]
            s0 = float(shift_of(p, x[p:p + 1])[0])
            lim = max(max_shift, s0)

            def evaluate(cands):
                """Score (T,) for candidate positions (T, 3); -inf where a constraint fails."""
                T = len(cands)
                ok = np.linalg.norm(cands - xin[p], axis=1) <= max_slide + 1e-12
                ok &= np.abs((cands - xin[p]) @ n_in) <= max_normal + 1e-12
                if np.isfinite(lim):
                    ok &= shift_of(p, cands) <= lim + 1e-12
                qs = np.empty((T, len(ks)))
                for a, k in enumerate(ks):
                    ids = el[k, :cnt[k]]
                    pts = np.broadcast_to(x[ids], (T, len(ids), 3)).copy()
                    pts[:, list(ids).index(p), :] = cands
                    qs[:, a] = element_quality(pts)
                    cm = _corner_dots(pts, nrm[a]).min(axis=1)
                    ok &= (cm > 1e-9) | (cm >= cmin0[a] - 1e-12)
                    ok &= qs[:, a] >= floors[a] - 1e-12
                    if cnt[k] == 4:
                        ok &= _angles_lengths(pts)[0].max(axis=1) <= max(max_quad_angle, amax0[a]) + 1e-9
                score = qs.min(axis=1) + 0.01 * qs.mean(axis=1)
                return np.where(ok, score, -np.inf), qs

            base_score = float(q0.min() + 0.01 * q0.mean())
            best = (base_score + 1e-4, None, None)
            for d in _directions(x, p, nb[p], apex.get(p), ref.get(p)):
                center, half = 0.0, max_slide
                for rounds in range(3):
                    ts = center + np.linspace(-half, half, n_grid if rounds == 0 else 21)
                    cands = x[p][None, :] + ts[:, None] * d[None, :]
                    sc, qs = evaluate(cands)
                    pen = sc - 1e-3 * np.abs(ts) / max(max_slide, 1e-12)     # prefer the shortest slide among ties
                    j = int(np.argmax(pen))
                    if not np.isfinite(pen[j]):
                        break
                    if pen[j] > best[0]:
                        best = (pen[j], cands[j].copy(), qs[j].copy())
                    center, half = float(ts[j]), 2.0 * (ts[1] - ts[0])
            if best[1] is not None:
                x[p] = best[1]
                for a, k in enumerate(ks):
                    q[k] = best[2][a]
                changed += 1
        moved_total += changed
        sweeps_done += 1
        if changed == 0:
            break
    moved = np.linalg.norm(x - xin, axis=1)
    # displacement out of the local surface: component of each move along the mean normal at entry
    out_of_surface = 0.0
    for i in np.nonzero(moved > 0)[0]:
        n = sum(_normal(xin[el[k, :cnt[k]]]) for k in inc[i])
        if np.linalg.norm(n) > 0:
            out_of_surface = max(out_of_surface, abs(float((x[i] - xin[i]) @ n) / np.linalg.norm(n)))
    cum = np.array([shift_of(i, x[i:i + 1])[0] for i in range(nn)])
    rep = {"nodes_moved": int((moved > 0).sum()), "moves": moved_total, "sweeps": sweeps_done,
           "max_slide_m": float(moved.max()) if nn else 0.0, "max_out_of_surface_m": out_of_surface,
           "max_cumulative_shift_m": float(cum.max()) if nn else 0.0,
           "quality_min_before": float(q_before.min()), "quality_min_after": float(q.min()),
           "inverted": int(sum(1 for k in range(len(el)) if np.dot(_normal(x[el[k, :cnt[k]]]), _normal(xin[el[k, :cnt[k]]])) <= 0))}
    return x, rep


# --------------------------------------------------------------------------------------------------
# RawMesh wrapper
# --------------------------------------------------------------------------------------------------
def collapse_and_optimize_raw_mesh(raw, element_size, ratio=0.15, min_aspect=3.0, *, max_shift=0.05,
                                   max_slide=None, quad_floor=0.5, max_quad_angle=120.0, max_normal=0.005, optimize=True):
    """Drop-in for ``collapse_raw_mesh`` using ``collapse_capped`` and, if ``optimize``, ``optimize_collapsed``.

    ``max_shift`` bounds the cumulative shift of a merged node during the collapse; the optimizer may slide
    a node up to ``max_slide`` (default ``0.6 * element_size``) along its line, bounded in total by
    ``max_shift + max_slide`` from the original positions.  Returns ``(new_raw, report)``.
    """
    from dataclasses import replace

    from .windio2_basis import shell_bases

    if max_slide is None:
        max_slide = 0.6 * element_size
    fixed = {int(i) for chain in raw.attachments.values() for i in chain}
    nodes, elements, src, rep = collapse_capped(
        raw.nodes, raw.elements, max_edge_len=ratio * element_size, fixed=fixed, min_aspect=min_aspect,
        max_shift=max_shift)
    if rep["merged"] == 0:
        return raw, rep
    nmap = rep["node_map"]
    if optimize:
        fixed_new = {int(nmap[i]) for i in fixed}
        nodes, orep = optimize_collapsed(nodes, elements, fixed=fixed_new, nodes0=raw.nodes, node_map=nmap,
                                         max_slide=max_slide, max_shift=max_shift + max_slide, quad_floor=quad_floor,
                                         max_quad_angle=max_quad_angle, max_normal=max_normal)
        rep = {**rep, "optimizer": orep}
    stacks = None if raw.material_stacks is None else tuple(raw.material_stacks[k] for k in src)
    new = replace(
        raw, nodes=nodes, elements=elements,
        cells=tuple(raw.cells[k] for k in src), surfaces=tuple(raw.surfaces[k] for k in src),
        element_parameters=np.asarray(raw.element_parameters)[src],
        bases=shell_bases(nodes, elements),
        attachments={key: tuple(int(nmap[i]) for i in chain) for key, chain in raw.attachments.items()},
        material_stacks=stacks)
    return new, rep
