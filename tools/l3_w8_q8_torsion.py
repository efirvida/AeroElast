#!/usr/bin/env python
"""W8 / L3 T11: 2D St-Venant torsion on the SAME 2D mesh and mid-side nodes as the 3D C3D20R model, with 8-node serendipity quads (the in-plane
element of the solid) instead of the bilinear Q4 of ``torsion_2d``. Separates 2D-Q4 discretisation from a 3D/2D model difference.
Run: python tools/l3_w8_q8_torsion.py ETA ARC H TMAX [GAUSS=3]"""
import math
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

sys.path.insert(0, str(Path(__file__).parent))
import l3_solid_section_ccx as T  # noqa: E402

# serendipity Q8, nodes: corners (-1,-1) (1,-1) (1,1) (-1,1), mids (0,-1) (1,0) (0,1) (-1,0)
XI = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1], [0, -1], [1, 0], [0, 1], [-1, 0]], float)


def shape(a, b):
    N = np.zeros(8)
    dN = np.zeros((2, 8))
    for i, (xi, eta) in enumerate(XI):
        if i < 4:
            N[i] = 0.25 * (1 + xi * a) * (1 + eta * b) * (xi * a + eta * b - 1)
            dN[0, i] = 0.25 * xi * (1 + eta * b) * (2 * xi * a + eta * b)
            dN[1, i] = 0.25 * eta * (1 + xi * a) * (xi * a + 2 * eta * b)
        elif xi == 0:
            N[i] = 0.5 * (1 - a * a) * (1 + eta * b)
            dN[0, i] = -a * (1 + eta * b)
            dN[1, i] = 0.5 * (1 - a * a) * eta
        else:
            N[i] = 0.5 * (1 + xi * a) * (1 - b * b)
            dN[0, i] = 0.5 * xi * (1 - b * b)
            dN[1, i] = -b * (1 + xi * a)
    return N, dN


def torsion_q8(sec, m2, ng=3):
    xy, quads = m2["xy"], m2["quads"]
    mid, nid = {}, len(xy)
    coords = list(xy)
    for ids, *_ in quads:
        for a, b in zip(ids, ids[1:] + ids[:1], strict=False):
            key = (min(a, b), max(a, b))
            if key not in mid:
                mid[key] = nid
                nid += 1
                coords.append(T.edge_mid(sec, m2, *key))
    coords = np.array(coords)
    gp = {2: [-1 / math.sqrt(3), 1 / math.sqrt(3)], 3: [-math.sqrt(0.6), 0.0, math.sqrt(0.6)]}[ng]
    gw = {2: [1.0, 1.0], 3: [5 / 9, 8 / 9, 5 / 9]}[ng]
    I, J, V = [], [], []
    f = np.zeros(nid)
    const = 0.0
    for ids, mat, tv, _a, _tag in quads:
        m = sec.materials[mat]
        G12, G13 = (m["G"][0], m["G"][1]) if np.ndim(m["E"]) else (m["G"], m["G"])
        tv = np.asarray(tv) / np.linalg.norm(tv)
        nv = np.array([-tv[1], tv[0]])
        Gt = G12 * np.outer(tv, tv) + G13 * np.outer(nv, nv)
        # edge order of ids: (0-1),(1-2),(2-3),(3-0) = serendipity mids 4..7
        nodes = list(ids) + [mid[(min(a, b), max(a, b))] for a, b in zip(ids, ids[1:] + ids[:1], strict=False)]
        X = coords[nodes]
        Ke, fe = np.zeros((8, 8)), np.zeros(8)
        for ia, a in enumerate(gp):
            for ib, b in enumerate(gp):
                N, dN = shape(a, b)
                Jm = dN @ X
                dets = np.linalg.det(Jm) * gw[ia] * gw[ib]
                B = np.linalg.solve(Jm, dN)
                g0 = np.array([-(N @ X[:, 1]), N @ X[:, 0]])
                Ke += B.T @ Gt @ B * dets
                fe += -(B.T @ Gt @ g0) * dets
                const += g0 @ Gt @ g0 * dets
        for i in range(8):
            f[nodes[i]] += fe[i]
            for j in range(8):
                I.append(nodes[i])
                J.append(nodes[j])
                V.append(Ke[i, j])
    K = sp.csr_matrix((V, (I, J)), shape=(nid, nid))
    free = np.arange(1, nid)
    phi = np.zeros(nid)
    phi[free] = spl.spsolve(K[free][:, free].tocsc(), f[free])
    return float(const - f @ phi)


def main():
    eta, arc, h, tmax = (float(x) for x in sys.argv[1:5])
    ng = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    sec = T.blade_section(eta)
    m2 = T.build_2d(sec, h, tmax, "bonded", {"mat": "adhesive", "arc": arc} if arc > 0 else None)
    sf = T.deck_metrics(T.deck_station(eta))[1]["GJ"]
    g4, g8 = T.torsion_2d(sec, m2), torsion_q8(sec, m2, ng)
    print(f"eta {eta} arc {arc} h {h} tmax {tmax}: GJ 2D-Q4 {g4 / sf:.4f}  2D-Q8({ng}x{ng}) {g8 / sf:.4f} of deck GJ_sc")


if __name__ == "__main__":
    main()
