"""FEM-free multi-cell St-Venant torsion constant GJ of a cached v2 mesh, at planes z = const, against the BeamDyn deck K[5,5] and (optionally)
against the shell/CalculiX twist profile of tools/ccx_s7_torsion_gmsh.py.

Model: every shell element cut by the plane gives one wall segment of length ds and shear stiffness A66 = sum_plies Qbar66(theta) t (ply G12 and
fibre angle from the span); the segments form a planar graph whose bounded faces are the cells (skin ring split by the webs). Compatibility
in each cell, sum_over_boundary (q_i - q_j) ds / A66 = 2 A_i * (twist rate = 1), is solved for the shear flows; GJ = sum 2 A_i q_i.
usage: python tools/bredt_batho_bands.py MESHER ELEMENT_SIZE BEAMDYN_BLADE_DAT L [twist.json]
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
mesher, es, deck_path, L = sys.argv[1], float(sys.argv[2]), sys.argv[3], float(sys.argv[4])
os.environ["V2_BLADE_MESHER"] = mesher
from tests.support.v2_blade import _inset, load_raw  # noqa: E402

raw = load_raw(es, mesher=mesher, generate=False, es_range=(0.05, 0.25))
X = np.asarray(raw["nodes"], float) if os.environ.get("BREDT_REF") == "oml" else _inset(raw)   # BREDT_REF=oml: outer mould line, no inset to the mid-surface
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


a66_layup = np.array([sum(qbar66(tables["materials"][n], a) * th for n, th, a in lay) for lay in tables["layups"]])
a66 = a66_layup[raw["layup_id"]]

deck = open(deck_path).read().split("\n")
i0 = [i for i, ln in enumerate(deck) if "DISTRIBUTED PROPERTIES" in ln][0] + 1
n = int([ln for ln in deck if "station_total" in ln][0].split()[0])
eta, K66, i = [], [], i0
while i < len(deck) and len(eta) < n:
    t = deck[i].split()
    if len(t) == 1:
        eta.append(float(t[0]))
        K66.append(np.array([[float(x) for x in deck[i + 1 + r].split()] for r in range(6)])[5, 5])
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
    for s, (u, v, g) in enumerate(segs):
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
        for a, b, s in f:
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


zs = np.linspace(0.05 * L, 0.95 * L, 19)
res = []
print(" z[m]  cells  area[m2]   GJ_bredt[N m2]   K66_deck      bredt/K66")
for zc in zs:
    g, nc, ar = gj_at(zc + 0.0137)
    kd = np.interp(zc / L, eta, K66)
    res.append((zc, nc, g, kd))
    print("%6.1f  %3d   %7.2f   %.4e   %.4e   %7.3f" % (zc, nc, ar, g, kd, g / kd))
if len(sys.argv) > 5:
    d = json.load(open(sys.argv[5]))
    z, th = np.array(d["z"]), np.radians(np.array(d["shell_deg"]))
    thc = np.radians(np.array(d["ccx_deg"]))
    print("\nFE torsion stiffness over z intervals (T=1000 N m; GJ = T dz / dtheta) vs Bredt at the midpoint and the deck:")
    for a, b in ((20, 40), (40, 60), (60, 80), (80, 100)):
        sa, sb = np.argmin(abs(z - a)), np.argmin(abs(z - b))
        g_sh = 1000.0 * (z[sb] - z[sa]) / (th[sb] - th[sa])
        g_cx = 1000.0 * (z[sb] - z[sa]) / (thc[sb] - thc[sa])
        zm = 0.5 * (z[sa] + z[sb])
        gb = np.mean([1.0 / gj_at(zz + 0.0137)[0] for zz in np.linspace(z[sa], z[sb], 9)]) ** -1   # harmonic: the intervals are springs in series
        kd = np.mean(np.interp(np.linspace(z[sa], z[sb], 50) / L, eta, K66) ** -1) ** -1
        print("z %3d-%3d  FEshell %.3e  FEccx %.3e  Bredt %.3e  deck(harm.) %.3e | shell/deck %.3f ccx/deck %.3f bredt/deck %.3f bredt/shell %.3f"
              % (a, b, g_sh, g_cx, gb, kd, g_sh / kd, g_cx / kd, gb / kd, gb / g_sh))
