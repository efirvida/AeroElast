"""W4/T5: mesh-cut multi-cell Bredt GJ (copy of the T7b recipe logs/chord_scaled/t7b/bredt_stations.py, gj_at unchanged) as a function.
make(raw, deck_path) -> (gj_at, eta_deck, K66_deck).  Mid-surface cut (``_inset``)."""

import json

import numpy as np

from tests.support.v2_blade import _inset


def make(raw, deck_path):
    X = _inset(raw)
    E = np.asarray(raw["elements"])
    tables = json.loads(bytes(raw["tables_json"]).decode())

    def qbar66(m, deg):
        g = lambda v, i: float(np.atleast_1d(v)[min(i, np.size(v) - 1)])  # noqa: E731
        E1, E2, G12, nu12 = g(m["E"], 0), g(m["E"], 1), g(m["G"], 0), g(m["nu"], 0)
        nu21 = nu12 * E2 / E1
        d = 1 - nu12 * nu21
        q11, q22, q12, q66 = E1 / d, E2 / d, nu12 * E2 / d, G12
        c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
        return (q11 + q22 - 2 * q12 - 2 * q66) * s * s * c * c + q66 * (s**4 + c**4)

    a66_layup = np.array(
        [
            sum(qbar66(tables["materials"][n], a) * th for n, th, a in lay)
            for lay in tables["layups"]
        ]
    )
    a66 = a66_layup[raw["layup_id"]]

    deck = open(deck_path).read().split("\n")
    i0 = [i for i, ln in enumerate(deck) if "DISTRIBUTED PROPERTIES" in ln][0] + 1
    n = int([ln for ln in deck if "station_total" in ln][0].split()[0])
    eta, K66, i = [], [], i0
    while i < len(deck) and len(eta) < n:
        t = deck[i].split()
        if len(t) == 1:
            eta.append(float(t[0]))
            K66.append(
                np.array([[float(x) for x in deck[i + 1 + r].split()] for r in range(6)])[5, 5]
            )
            i += 14
        else:
            i += 1
    eta, K66 = np.array(eta), np.array(K66)

    zn = X[:, 2]
    zmin = np.array([zn[e[e >= 0]].min() for e in E])
    zmax = np.array([zn[e[e >= 0]].max() for e in E])

    def gj_at(zc):
        cand = np.nonzero((zmin < zc) & (zmax > zc))[0]
        pts, segs = {}, []
        for k in cand:
            ids = E[k][E[k] >= 0]
            cut = []
            for a, b in zip(ids, np.roll(ids, -1), strict=True):
                za, zb = zn[a] - zc, zn[b] - zc
                if za * zb < 0:
                    key = (min(a, b), max(a, b))
                    w = za / (za - zb)
                    pts.setdefault(key, X[a, :2] + w * (X[b, :2] - X[a, :2]))
                    cut.append(key)
            if len(cut) == 2 and cut[0] != cut[1]:
                segs.append((cut[0], cut[1], a66[k]))
        keys = list(pts)
        ix = {k: j for j, k in enumerate(keys)}
        P = np.array([pts[k] for k in keys])
        adj = {j: [] for j in range(len(keys))}
        for s, (u, v, _g) in enumerate(segs):
            adj[ix[u]].append((ix[v], s))
            adj[ix[v]].append((ix[u], s))
        # half-edge faces
        ang = {}
        for u, lst in adj.items():
            ang[u] = sorted(lst, key=lambda t: np.arctan2(*(P[t[0]] - P[u])[::-1]))
        visited, faces = set(), []
        for u in adj:
            for v, s in adj[u]:
                if (u, v, s) in visited:
                    continue
                face, cu, cv, cs = [], u, v, s
                while (cu, cv, cs) not in visited:
                    visited.add((cu, cv, cs))
                    face.append((cu, cv, cs))
                    lst = ang[cv]
                    j = [t for t, (w, ss) in enumerate(lst) if w == cu and ss == cs][0]
                    nw, ns = lst[(j - 1) % len(lst)]
                    cu, cv, cs = cv, nw, ns
                faces.append(face)
        area = [0.5 * sum(P[a, 0] * P[b, 1] - P[b, 0] * P[a, 1] for a, b, _ in f) for f in faces]
        cells = [f for f, ar in zip(faces, area, strict=True) if ar > 1e-6]
        ca = [ar for ar in area if ar > 1e-6]
        owner = {}
        for c, f in enumerate(cells):
            for _a, _b, s in f:
                owner.setdefault(s, []).append(c)
        nc = len(cells)
        M, rhs = np.zeros((nc, nc)), 2 * np.array(ca)
        for c, f in enumerate(cells):
            for a, b, s in f:
                ds = np.linalg.norm(P[a] - P[b]) / segs[s][2]
                M[c, c] += ds
                others = [o for o in owner[s] if o != c]
                if others:
                    M[c, others[0]] -= ds
        q = np.linalg.solve(M, rhs)
        return float(np.sum(2 * np.array(ca) * q)), nc, float(sum(ca))

    return gj_at, eta, K66
