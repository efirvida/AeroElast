"""Mesh-free Bredt and BeamDyn-deck helpers shared by the L3 torsion tools.

Kept from ``tools/l3_bredt_yaml.py`` when that windIO 1.x section reader was deleted with
``tools/audit_blade_input.py`` (odd/tasks/gmsh-official-mesher.md, S9, L-S66). Geometry arrays in,
numbers out; no yaml schema and no mesh. BeamDyn 6x6 order is [Fx Fy Fz Mx My Mz]: 0,1 shear,
2 axial, 3,4 bending, 5 torsion.
"""

from __future__ import annotations

import numpy as np


def _cross_sum(p: np.ndarray, q: np.ndarray) -> float:
    return float(np.sum(p[:, 0] * q[:, 1] - p[:, 1] * q[:, 0]))


def multicell_bredt(pts, seg_gt, seg_cell, webs, seg_t=None) -> dict:
    """Multi-cell Bredt GJ for unit twist rate.

    ``pts`` closed loop vertices (M, 2); segment i joins vertex i to i+1 (mod M), shear stiffness
    ``seg_gt[i]`` and cell ``seg_cell[i]``. ``webs`` sorted TE to LE, each ``{i_s, i_p, gt}`` with the
    suction and pressure vertex indices; web j separates cell j (TE side) and j+1, and is traversed
    suction -> pressure for cell j.
    """
    pts = np.asarray(pts, float)
    pts = pts - pts.mean(axis=0)
    seg_gt, seg_cell = np.asarray(seg_gt, float), np.asarray(seg_cell, int)
    nc = len(webs) + 1
    nxt = np.roll(pts, -1, axis=0)
    length = np.linalg.norm(nxt - pts, axis=1)
    sigma = np.sign(_cross_sum(pts, nxt))
    area = np.zeros(nc)
    flex_ext = np.zeros(nc)
    len_ext = np.zeros(nc)
    for c in range(nc):
        sel = seg_cell == c
        area[c] += 0.5 * _cross_sum(pts[sel], nxt[sel])
        flex_ext[c] = float(np.sum(length[sel] / seg_gt[sel]))
        len_ext[c] = float(length[sel].sum())
    f = np.zeros((nc, nc))
    web_out = []
    for j, w in enumerate(webs):
        a, b = pts[w["i_s"]], pts[w["i_p"]]
        cr = a[0] * b[1] - a[1] * b[0]
        area[j] += 0.5 * cr
        area[j + 1] -= 0.5 * cr
        lw = float(np.linalg.norm(b - a))
        fw = lw / w["gt"]
        f[j, j + 1] -= fw
        f[j + 1, j] -= fw
        f[j, j] += fw
        f[j + 1, j + 1] += fw
        web_out.append({"length": lw, "flex": fw})
    area = sigma * area
    f[np.arange(nc), np.arange(nc)] += flex_ext
    q = np.linalg.solve(f, 2.0 * area)
    open_term = float(np.sum(seg_gt * (seg_t if seg_t is not None else 0.0) ** 2 * length) / 3.0)
    for w, wo in zip(webs, web_out, strict=True):
        open_term += w.get("t", 0.0) ** 2 * w["gt"] * wo["length"] / 3.0
    return {"gj": float(2.0 * area @ q), "q": q, "areas": area, "n_cells": nc, "F": f,
            "flex_ext": flex_ext, "len_ext": len_ext, "webs": web_out, "open": open_term}


def condense_shear(k) -> float:
    """K66 condensed over the two shear strains (indices 0, 1) of a 6x6 stiffness matrix."""
    k = np.asarray(k, float)
    s = k[np.ix_([0, 1], [0, 1])]
    c = k[[0, 1], 5]
    return float(k[5, 5] - c @ np.linalg.solve(s, c))


# =============================================================================================
# yaml section
# =============================================================================================


# =============================================================================================
# BeamDyn deck
# =============================================================================================
def load_beamdyn_6x6(path) -> tuple[np.ndarray, np.ndarray]:
    """(etas, K) with K (n, 6, 6): each station is an eta line, 6 lines of K, 6 lines of M."""
    lines = [ln.strip() for ln in open(path)]
    i = next(k for k, ln in enumerate(lines) if "DISTRIBUTED PROPERTIES" in ln) + 1
    etas, ks = [], []
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) == 1:
            try:
                eta = float(parts[0])
            except ValueError:
                i += 1
                continue
            rows = []
            j = i + 1
            while len(rows) < 6:
                if lines[j]:
                    rows.append([float(x) for x in lines[j].split()])
                j += 1
            etas.append(eta)
            ks.append(rows)
            i = j
            seen = 0
            while seen < 6:  # skip the mass block
                if lines[i]:
                    seen += 1
                i += 1
        else:
            i += 1
    k = np.asarray(ks, float)
    if not (k[:, 2, 2] > k[:, 0, 0]).all():
        raise ValueError("BeamDyn K index convention check failed: expected K[2,2] (EA) > K[0,0] (shear)")
    return np.asarray(etas), k


def interp_deck(etas, ks, eta) -> np.ndarray:
    etas, ks = np.asarray(etas, float), np.asarray(ks, float)
    return np.array([[np.interp(eta, etas, ks[:, r, c]) for c in range(6)] for r in range(6)])
