"""Cut collapsing strips out of a windIO v2 raw shell mesh.

Where two source lines run closer than a fraction of the element size (the TE seam and its first
event line, a layer edge that creeps towards another), the strip between them is meshed as a row of
quads whose short side is millimetres and whose long side is an element: aspect ratio 15-40.  No
smoothing can fix that, because both lines are source events and cannot move apart.

``collapse_thin_strips`` merges the two ends of every short *cross* edge.  Inside a thin band both
long sides end up on one line and the band's quads degenerate to a segment (dropped).  Where the band
opens up again, the last quad keeps one short edge merged and one long edge, so it becomes a
triangle: the collapsing line is cut with triangles.  Everything is local and general; nothing knows
about a particular blade.

What it changes, stated: the material boundary of a removed band moves by at most half the merge
length (``max_edge_len / 2``) and the removed band's area is reported.  A merge that would invert an
incident element, or move a ``fixed`` node, is refused.
"""

from __future__ import annotations

import numpy as np

__all__ = ["collapse_thin_strips", "collapse_raw_mesh"]


def _pad(elements):
    el = np.asarray(elements, dtype=int)
    if el.ndim != 2 or el.shape[1] not in (3, 4):
        raise ValueError("elements must be (n, 3) or (n, 4); triangles in a quad array are padded with -1")
    if el.shape[1] == 3:
        el = np.column_stack([el, -np.ones(len(el), dtype=int)])
    return el


def _normal(p):
    """Newell normal (area-weighted) of a polygon given as (k, 3)."""
    q = np.roll(p, -1, axis=0)
    return 0.5 * np.array([
        np.sum((p[:, 1] - q[:, 1]) * (p[:, 2] + q[:, 2])),
        np.sum((p[:, 2] - q[:, 2]) * (p[:, 0] + q[:, 0])),
        np.sum((p[:, 0] - q[:, 0]) * (p[:, 1] + q[:, 1])),
    ])


def collapse_thin_strips(nodes, elements, *, max_edge_len, fixed=(), min_aspect=3.0):
    """Merge the ends of short cross edges of sliver elements.

    ``nodes`` (n, 3); ``elements`` (m, 4) with -1 padding for triangles (or (m, 3)); ``max_edge_len``
    metres: an edge is a collapse candidate when it is shorter than this and belongs to an element of
    aspect ratio at least ``min_aspect``; ``fixed``: node ids that never move (they absorb their
    partner).

    Returns ``(nodes, elements, source, report)``: compacted node array, compacted padded elements,
    ``source[k]`` = index of the input element that became output element ``k``, and a report dict
    (merged, dropped, to_triangle, removed_area_m2, max_node_shift_m, ...).
    """
    x0 = np.asarray(nodes, dtype=float)
    el0 = _pad(elements)
    pos = x0.copy()
    rep = np.arange(len(x0))
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
    merged = 0
    shift = 0.0
    for (a, b), _ in sorted(cand.items(), key=lambda kv: kv[1]):
        ra, rb = find(a), find(b)
        if ra == rb or (is_fixed[ra] and is_fixed[rb]):
            continue
        keep, gone = (ra, rb) if (is_fixed[ra] or (not is_fixed[rb] and ra < rb)) else (rb, ra)
        new = pos[keep] if is_fixed[keep] else 0.5 * (pos[ra] + pos[rb])
        if np.linalg.norm(pos[gone] - new) > max_edge_len or np.linalg.norm(pos[keep] - new) > max_edge_len:
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
        incident[keep] |= incident[gone]
        merged += 1
        shift = max(shift, float(np.linalg.norm(new - x0[keep])), float(np.linalg.norm(new - x0[gone])))

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
        return x0, el0, np.arange(len(el0)), {"merged": 0, "dropped": 0, "to_triangle": 0,
                                              "removed_area_m2": 0.0, "max_node_shift_m": 0.0,
                                              "node_map": np.arange(len(x0))}
    out_el = np.array(out_el, dtype=int).reshape(-1, 4)
    used = np.unique(out_el[out_el >= 0])
    remap = -np.ones(len(x0), dtype=int)
    remap[used] = np.arange(len(used))
    out_el = np.where(out_el >= 0, remap[np.where(out_el >= 0, out_el, 0)], -1)
    report = {"merged": merged, "dropped": dropped, "to_triangle": to_tri,
              "removed_area_m2": float(removed), "max_node_shift_m": shift,
              "nodes_before": len(x0), "nodes_after": len(used),
              "node_map": np.array([remap[find(i)] for i in range(len(x0))], dtype=int)}
    return pos[used], out_el, np.array(source, dtype=int), report


def collapse_raw_mesh(raw, element_size, ratio=0.15, min_aspect=3.0):
    """``collapse_thin_strips`` on a ``RawMesh``; returns ``(new_raw, report)``.

    Skin nodes that carry a web attachment chain are fixed.  Cells, surfaces, chart parameters and
    material stacks follow their elements; bases are recomputed on the new connectivity and the
    attachment chains are remapped onto the merged node ids.
    """
    from dataclasses import replace

    from .windio2_basis import shell_bases

    fixed = {int(i) for chain in raw.attachments.values() for i in chain}
    nodes, elements, src, rep = collapse_thin_strips(
        raw.nodes, raw.elements, max_edge_len=ratio * element_size, fixed=fixed, min_aspect=min_aspect)
    if rep["merged"] == 0:
        return raw, rep
    nmap = rep["node_map"]
    stacks = None if raw.material_stacks is None else tuple(raw.material_stacks[k] for k in src)
    new = replace(
        raw, nodes=nodes, elements=elements,
        cells=tuple(raw.cells[k] for k in src), surfaces=tuple(raw.surfaces[k] for k in src),
        element_parameters=np.asarray(raw.element_parameters)[src],
        bases=shell_bases(nodes, elements),
        attachments={key: tuple(int(nmap[i]) for i in chain) for key, chain in raw.attachments.items()},
        material_stacks=stacks)
    return new, rep
