#!/usr/bin/env python
"""L3 / T4: CalculiX solid arbiter for ONE prismatic wind-turbine-blade section (odd/tasks/l3-outboard-section-stiffness.md, L15).

A 2D section (OML x chord, layers grown INWARD from the OML by their thickness, webs as straight solid walls with
their own layers, centred on the web chart line) is meshed with a structured strip mesh (normal columns, one element
row per layer sub-interval) and extruded into a C3D20R solid of length L = n_chord * chord. Orthotropic plies carry the
local frame 1 = span (z), 2 = in-plane wall tangent, 3 = wall normal (windIO E1/E2/E3, G12/G13/G23, nu12/nu13/nu23).
Not SONATA: the layer geometry comes from the windIO v2 topology reader (``WindIO2Topology.section``) only.

Loads: end A clamped; end B rigid body (reference node at the section centroid); four step load cases (axial force,
torque, flap-axis moment Mx, My). Generalised strains (axial strain, twist rate, two curvatures) are least-squares plane
fits of the displacement of the section nodes at two planes of the middle third; K4 = inv(compliance 4x4) is the
shear-free section stiffness for (F, T, Mx, My) about the reference point.

Cases: ``tube_iso`` / ``tube_ortho`` / ``tube_web`` (pipeline proof with closed forms) and ``blade`` (--eta).

  scripts/aeroenv.sh env PYTHONPATH=$PWD python tools/l3_solid_section_ccx.py run --case blade --eta 0.5 --out DIR
  ... build | analyze  (build writes DIR/mesh.inp, case_{F,T,MX,MY}.inp, meta.json; analyze reads case_*.dat)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
YAML = REPO / "tests/IEA15MW/reference/upstream/IEA-15-240-RWT.yaml"
DECK = REPO / "tests/IEA15MW/reference/IEA-15-240-RWT_BeamDyn_blade.dat"
CCX = os.environ.get("CCX_BIN", "/scratch/leahk/eduardo.donestevez/venv/bin/ccx")
TOL_LEVEL = 1e-6  # m, depth levels closer than this are one level
_trapz = getattr(np, "trapezoid", None) or np.trapz


# ------------------------------------------------------------------------------------------------ section model
class Section:
    """P(s) -> OML xy for s in [0,1]; cells (s0, s1, [(mat, t)] OML->inward); webs; materials."""

    def __init__(self, name, P, cells, webs, materials, closed=False, le_point=None, te_point=None):
        self.name, self.P_, self.cells, self.webs, self.materials = name, P, cells, webs, materials
        self.closed, self.le_point, self.te_point = closed, le_point, te_point
        s = np.linspace(0.0, 1.0, 20001)
        xy = P(s)
        self.s_dense, self.xy_dense = s, xy
        area2 = np.sum(xy[:-1, 0] * xy[1:, 1] - xy[1:, 0] * xy[:-1, 1])
        self.ccw = area2 > 0
        from scipy.interpolate import CubicSpline

        self._cs = CubicSpline(s, xy, bc_type="periodic" if closed else "not-a-knot")
        seg = np.linalg.norm(np.diff(xy, axis=0), axis=1)
        self.arc = np.concatenate([[0.0], np.cumsum(seg)])

    def P(self, s):
        return self._cs(s)

    def normal_in(self, s):
        t = self._cs(s, 1)
        t = t / np.linalg.norm(t, axis=-1, keepdims=True)
        n = np.stack([-t[..., 1], t[..., 0]], axis=-1)  # left normal
        return n if self.ccw else -n

    def tangent(self, s):
        t = self._cs(s, 1)
        return t / np.linalg.norm(t, axis=-1, keepdims=True)

    def cell_at(self, s):
        for c in self.cells:
            if c[0] - 1e-12 <= s <= c[1] + 1e-12:
                return c
        raise ValueError(f"no cell at s={s}")

    def depth_total(self, s):
        return sum(t for _, t in self.cell_at(s)[2])

    def arc_at(self, s):
        return float(np.interp(s, self.s_dense, self.arc))


def _circle(R, theta0=0.5 * math.pi):
    return lambda s: np.stack([R * np.cos(theta0 + 2 * math.pi * np.asarray(s)),
                               R * np.sin(theta0 + 2 * math.pi * np.asarray(s))], axis=-1)


def tube_section(kind, Ro=0.5, t=0.02):
    """Circular tube (CCW). iso: E=70 GPa; ortho: transversely isotropic E1=100 GPa along the span; web: iso tube + diametral web."""
    if kind in ("tube_ortho", "tube_ortho_g13"):
        g13 = 2e9 if kind == "tube_ortho_g13" else 5e9  # G13 != G12 discriminates the local-axis mapping (GJ must stay G12*J)
        mats = {"UD": {"E": (100e9, 10e9, 10e9), "nu": (0.3, 0.3, 0.4), "G": (5e9, g13, 3.5714e9)}}
        layers = [("UD", t / 2), ("UD", t / 2)]
    else:
        mats = {"ISO": {"E": 70e9, "nu": 0.3, "G": 70e9 / 2.6}}
        layers = [("ISO", t / 2), ("ISO", t / 2)]
    webs = []
    if kind == "tube_web":
        tw = 0.03
        # horizontal web along x; circle starts at 90 deg, CCW: s=0.25 -> 180 deg (-Ro,0), s=0.75 -> 360 deg (Ro,0)
        webs = [{"name": "web", "A": np.array([-Ro, 0.0]), "B": np.array([Ro, 0.0]), "sA": 0.25, "sB": 0.75,
                     "layers": [("ISO", tw / 3), ("ISO", tw / 3), ("ISO", tw / 3)], "le_dir": np.array([0.0, 1.0])}]
    return Section(kind, _circle(Ro), [(0.0, 1.0, layers)], webs, mats, closed=True)


def tube_exact(kind, Ro=0.5, t=0.02):
    Ri = Ro - t
    A = math.pi * (Ro**2 - Ri**2)
    I = math.pi / 4 * (Ro**4 - Ri**4)
    J = 2 * I
    if kind in ("tube_ortho", "tube_ortho_g13"):
        E1, G12 = 100e9, 5e9
        return {"EA": E1 * A, "GJ": G12 * J, "EI_min": E1 * I, "EI_max": E1 * I}
    E, G = 70e9, 70e9 / 2.6
    ex = {"EA": E * A, "GJ": G * J, "EI_min": E * I, "EI_max": E * I}
    if kind == "tube_web":
        tw = 0.03
        y = np.linspace(-tw / 2, tw / 2, 200001)
        L = 2 * np.sqrt(np.maximum(Ri**2 - y**2, 0))
        Aw = _trapz(L, y)
        Ix = _trapz(L * y**2, y)  # about x axis (the web lies along x)
        Iy = _trapz(2 * (np.sqrt(np.maximum(Ri**2 - y**2, 0)) ** 3) / 3, y)  # integral of x^2 over the strip
        ex = {"EA": E * (A + Aw), "GJ": None, "EI_min": E * min(I + Ix, I + Iy), "EI_max": E * max(I + Ix, I + Iy)}
    return ex


def blade_section(eta, yaml_path=YAML):
    from aeroelast.models.blade.numad.objects.blade import Blade

    b = Blade()
    b.read_yaml(str(yaml_path))
    topo = b.windio2_topology
    sec = topo.section(eta)
    sd = np.linspace(0.0, 1.0, 20001)
    xyd = np.asarray(topo.skin_points(eta, sd), float)[:, :2]
    from scipy.interpolate import CubicSpline

    cs = CubicSpline(sd, xyd)
    P = lambda s: cs(np.asarray(s))  # noqa: E731
    mats = {}

    def mprop(m):
        E = m["E"]
        if np.ndim(E) == 0:
            return {"E": float(E), "nu": float(m["nu"]), "G": float(m["G"]) if "G" in m else float(E) / (2 * (1 + float(m["nu"])))}
        return {"E": tuple(float(v) for v in E), "nu": tuple(float(v) for v in m["nu"]), "G": tuple(float(v) for v in m["G"])}

    cells = []
    for c in sec.skin_cells:
        lay = []
        for l in c.layers:
            mats.setdefault(l.material_name, mprop(l.material))
            lay.append((l.material_name, float(l.thickness_m)))
        if cells and cells[-1][2] == lay and abs(cells[-1][1] - c.nd_bounds[0]) < 1e-12:
            cells[-1] = (cells[-1][0], float(c.nd_bounds[1]), lay)
        else:
            cells.append((float(c.nd_bounds[0]), float(c.nd_bounds[1]), lay))
    xy0 = P(0.0)
    d = np.linalg.norm(xyd - xy0, axis=1)
    le = xyd[np.argmax(d)]
    webs = []
    wc = {c.surface: c for c in sec.web_cells}
    for name, (pa, pb) in sec.webs:
        lay = []
        for l in wc[name].layers:
            mats.setdefault(l.material_name, mprop(l.material))
            lay.append((l.material_name, float(l.thickness_m)))
        webs.append({"name": name, "A": np.array(pa.xyz[:2]), "B": np.array(pb.xyz[:2]), "sA": float(pa.nd_arc), "sB": float(pb.nd_arc),
                         "layers": lay, "le_dir": le - xy0})
    return Section(f"blade_eta{eta}", P, cells, webs, mats, closed=False, le_point=le, te_point=xy0)


# ------------------------------------------------------------------------------------------------ 2D strip mesh
def refine_levels(levels, tmax):
    out = [levels[0]]
    for a, b in zip(levels[:-1], levels[1:], strict=False):
        n = max(1, math.ceil((b - a) / tmax - 1e-9))
        out += [a + (b - a) * i / n for i in range(1, n + 1)]
    return out


def unify(vals, tol=TOL_LEVEL):
    vals = sorted(vals)
    out = [vals[0]]
    for v in vals[1:]:
        if v - out[-1] > tol:
            out.append(v)
    return out


def valid_range(sec, ds_check=None):
    """TE-tip cut: first/last fine sample where the inward layer stack folds or overlaps (distance of the inner point to the OML
    boundary < 0.97 * depth). Returns (s_lo, s_hi) and the removed region polygon area."""
    from shapely.geometry import LinearRing, Point, Polygon

    s = np.linspace(0, 1, 4001)
    xy = sec.P(s)
    poly = Polygon(xy)
    ring = LinearRing(xy[:-1] if not sec.closed else xy[:-1])
    n = sec.normal_in(s)
    bad = np.zeros(len(s), bool)
    for i, si in enumerate(s):
        D = sec.depth_total(min(max(si, 0), 1))
        q = xy[i] + D * n[i]
        pt = Point(q)
        if not poly.contains(pt) or ring.distance(pt) < 0.97 * D:
            bad[i] = True
    if sec.closed:
        assert not bad.any(), "tube wall folds"
        return 0.0, 1.0, 0.0
    bi = np.where(bad)[0]
    lo_bad = bi[s[bi] < 0.5]
    hi_bad = bi[s[bi] > 0.5]
    assert (len(lo_bad) == 0 or s[lo_bad.max()] < 0.2) and (len(hi_bad) == 0 or s[hi_bad.min()] > 0.8), "fold away from the TE"
    s_lo = float(s[lo_bad.max() + 1]) if len(lo_bad) else 0.0
    s_hi = float(s[hi_bad.min() - 1]) if len(hi_bad) else 1.0
    # removed tip polygon: P(s<s_lo) ... P(s>s_hi) closed by the chord
    pts = np.vstack([sec.P(np.linspace(s_hi, 1.0, 400)), sec.P(np.linspace(0.0, s_lo, 400))])
    area = Polygon(pts).buffer(0).area if (s_lo > 0 or s_hi < 1) else 0.0
    return s_lo, s_hi, float(area)


def tip_setup(sec, h_s, f_min=0.25):
    """Trailing-edge tip closure. Inward layer stacks of the two sides overlap (or fold) near the sharp TE. Regular columns stop at s_c, the
    last SS station whose normal ray still finds room (ray length to the opposite surface >= D_ss + D_ps, no fold). Beyond it, tip stations
    s_j < s_c are rays along the SS normal: the SS layers are laid from the SS point, the PS layers from the ray's hit point Q_j along -u,
    all depth levels scaled by f_j = t_j / (D_ss + D_ps) when f_j < 1, so the innermost SS and PS nodes coincide (union of the wedge, no
    double counting). Stations stop where f_j would fall below f_min (wedge thinner than f_min of the two stacks); the face there is free (area beyond it is lost)."""
    if sec.closed:
        return None
    from shapely.geometry import LinearRing, LineString, Point, Polygon

    ring = LinearRing(sec.xy_dense)

    def partner(s):
        P, u = sec.P(s), sec.normal_in(s)
        o = P + 1e-6 * u
        inter = ring.intersection(LineString([o, P + 6.0 * u]))
        pts = [inter] if inter.geom_type == "Point" else [g for g in getattr(inter, "geoms", []) if g.geom_type == "Point"]
        if not pts:
            return None
        q = min(pts, key=lambda g: g.distance(Point(*o)))
        Q = np.array([q.x, q.y])
        return P, u, Q, float(np.linalg.norm(Q - P)), float(np.interp(ring.project(q), sec.arc, sec.s_dense))

    def room(s):
        r = partner(s)
        if r is None:
            return False, None
        P, u, Q, tt, sp = r
        Dss, Dps = sec.depth_total(s), sec.depth_total(sp)
        fold = ring.distance(Point(*(P + Dss * u))) < 0.97 * Dss
        return (tt >= (Dss + Dps) * (1 + 1e-3)) and not fold, r

    grid = np.linspace(0.0, 0.2, 2001)
    ok = [room(s)[0] for s in grid]
    assert ok[-1], "no room for the layer stacks at s=0.2 (TE closure not supported here)"
    ic = len(grid) - 1
    while ic > 0 and ok[ic - 1]:
        ic -= 1
    s_c = float(grid[ic])
    _, (P0, u0, Q0, t0, sp0) = room(s_c)
    t_min = f_min * (sec.depth_total(s_c) + sec.depth_total(sp0))
    ie = ic
    while ie > 0:
        r = partner(grid[ie - 1])
        if r is None or r[3] < t_min:
            break
        ie -= 1
    s_end = float(grid[ie])
    n_t = max(1, math.ceil((sec.arc_at(s_c) - sec.arc_at(s_end)) / min(h_s, 0.01) - 1e-9)) if ie < ic else 0
    st = []
    for j in range(n_t + 1):
        s = s_c if j == 0 else float(np.interp(sec.arc_at(s_c) + (sec.arc_at(s_end) - sec.arc_at(s_c)) * j / n_t, sec.arc, sec.s_dense))
        P, u, Q, tt, sp = partner(s)
        Dss, Dps = sec.depth_total(s), sec.depth_total(sp)
        f = 1.0 if j == 0 else min(1.0, tt / (Dss + Dps))
        st.append({"s": s, "P": P, "u": u, "Q": Q, "t": tt, "sp": sp, "Dss": Dss, "Dps": Dps, "f": f})
    sp_end = st[-1]["sp"]
    poly = Polygon(np.vstack([sec.P(np.linspace(st[-1]["s"], 0.0, 200)), sec.P(np.linspace(1.0, sp_end, 200))]))
    return {"s_c": s_c, "sp_c": sp0, "stations": st, "area": float(poly.buffer(0).area)}


def build_2d(sec, h_s, t_max, web_mode="bonded"):
    from scipy.optimize import brentq

    tip = tip_setup(sec, h_s)
    if tip:
        s_lo, s_hi, tip_area = tip["s_c"], tip["sp_c"], tip["area"]
    else:
        s_lo, s_hi, tip_area = 0.0, 1.0, 0.0
    cells = [(max(a, s_lo), min(b, s_hi), lay) for a, b, lay in sec.cells if b > s_lo and a < s_hi]
    # global depth levels
    lv = {0.0}
    for _, _, lay in cells:
        c = 0.0
        for _, t in lay:
            c += t
            lv.add(c)
    for _w in sec.webs:
        pass
    G = unify(list(lv))
    Dn = refine_levels(G, t_max)
    Dn = np.array(Dn)

    def kmax_of(D):
        return int(np.argmin(np.abs(Dn - D)))

    def inner(s, Dt):
        return sec.P(s) + Dt * sec.normal_in(s)

    # mandatory columns
    must = {s_lo, s_hi}
    for a, b, _ in cells:
        must.update([a, b])
    for w in sec.webs:
        tl = [t for _, t in w["layers"]]  # (LE skin, core, TE skin)
        u = (w["B"] - w["A"])
        w["len"] = float(np.linalg.norm(u))
        u = u / w["len"]
        wd = np.array([-u[1], u[0]])
        if wd @ w["le_dir"] < 0:
            wd = -wd
        w["u"], w["w"] = u, wd
        # across lines from the TE side to the LE side: TE skin, core, LE skin
        order = [2, 1, 0]
        bounds = [-(tl[1] / 2 + tl[2]), -tl[1] / 2, tl[1] / 2, tl[1] / 2 + tl[0]]
        lam, lay_of = [bounds[0]], []
        for qi, li in enumerate(order):
            n = max(1, math.ceil((bounds[qi + 1] - bounds[qi]) / t_max - 1e-9))
            for i in range(1, n + 1):
                lam.append(bounds[qi] + (bounds[qi + 1] - bounds[qi]) * i / n)
                lay_of.append(li)
        w["lam"], w["lay_of"] = lam, lay_of
        w["foot"] = {"A": [], "B": []}
        for end, s0 in (("A", w["sA"]), ("B", w["sB"])):
            for l in lam:
                def g(s):
                    D = 0.0 if web_mode == "through" else sec.depth_total(s)  # footprint on the OML (through) or the inner skin surface (bonded)
                    return float((inner(s, D) - w["A"]) @ wd) - l

                win = np.linspace(max(s_lo, s0 - 0.04), min(s_hi, s0 + 0.04), 801)
                gv = np.array([g(x) for x in win])
                sc = np.where(np.sign(gv[:-1]) * np.sign(gv[1:]) <= 0)[0]
                assert len(sc), f"web {w['name']} footprint not found"
                j = sc[np.argmin(np.abs(win[sc] - s0))]
                sf = brentq(g, win[j], win[j + 1], xtol=1e-14)
                w["foot"][end].append(sf)
                must.add(sf)
    must = sorted(must)
    # column list: subdivide by arc length
    cols = [must[0]]
    for a, b in zip(must[:-1], must[1:], strict=False):
        if b - a < 1e-10:
            continue
        L = sec.arc_at(b) - sec.arc_at(a)
        n = max(1, math.ceil(L / h_s - 1e-9))
        aa, bb = sec.arc_at(a), sec.arc_at(b)
        for i in range(1, n + 1):
            target = aa + (bb - aa) * i / n
            cols.append(float(np.interp(target, sec.arc, sec.s_dense)) if i < n else b)
    cols = np.array(cols)
    # intervals -> cell
    nc = len(cols)
    web_ivl = {}  # web_mode == 'through': skin intervals under the web footprint take the web layer (web has priority over the skin)
    if web_mode == "through":
        for w in sec.webs:
            for end in ("A", "B"):
                fc = [int(np.argmin(np.abs(cols - sf))) for sf in w["foot"][end]]
                for q in range(len(fc) - 1):
                    for i in range(fc[q], fc[q + 1]):
                        web_ivl[i] = (w["layers"][w["lay_of"][q]][0], w["u"], w["name"])
    ivl_cell = []
    for i in range(nc - 1):
        ivl_cell.append(sec.cell_at(0.5 * (cols[i] + cols[i + 1])) if True else None)
    # restrict to the trimmed range
    # node ids
    xy2, node_id, node_sk = [], {}, {}

    def node(key, xy, sk=None):
        if key not in node_id:
            node_id[key] = len(xy2)
            xy2.append(np.asarray(xy, float))
            if sk is not None:
                node_sk[node_id[key]] = sk
        return node_id[key]

    col_key = lambda c: 0 if (sec.closed and c == nc - 1) else c  # noqa: E731
    quads = []  # (n0..n3, mat, tangent)
    ktot_col = {}
    for i in range(nc - 1):
        cell = sec.cell_at(0.5 * (cols[i] + cols[i + 1]))
        tot = sum(t for _, t in cell[2])
        kt = kmax_of(tot)
        cum = np.cumsum([t for _, t in cell[2]])
        for c in (i, i + 1):
            ktot_col[c] = max(ktot_col.get(c, 0), kt)
        tvec = sec.P(cols[i + 1]) - sec.P(cols[i])
        tvec = tvec / np.linalg.norm(tvec)
        for k in range(kt):
            dmid = 0.5 * (Dn[k] + Dn[k + 1])
            li = int(np.searchsorted(cum, dmid))
            mat = cell[2][li][0]
            tv_k, tag_k = tvec, "skin"
            if i in web_ivl:
                mat, tv_k, tag_k = web_ivl[i]
            ids = []
            for c, kk in ((i, k), (i + 1, k), (i + 1, k + 1), (i, k + 1)):
                s = cols[c]
                if tip and c == nc - 1:  # PS end column: laid along -u from the ray's hit point (shared with the tip stations)
                    ids.append(node((col_key(c), kk), tip["stations"][0]["Q"] - Dn[kk] * tip["stations"][0]["u"]))
                else:
                    ids.append(node((col_key(c), kk), sec.P(s) + Dn[kk] * sec.normal_in(s), (float(s), float(Dn[kk]))))
            quads.append((ids, mat, tv_k, tag_k if tag_k == "skin" else "skin_" + tag_k))
    # trailing-edge tip stations (closure of the aft cell)
    if tip and len(tip["stations"]) > 1:
        st = tip["stations"]
        cellS = sec.cell_at(0.5 * (st[0]["s"] + st[-1]["s"]))
        cellP = sec.cell_at(0.5 * (st[0]["sp"] + st[-1]["sp"]))
        assert abs(sum(x for _, x in cellS[2]) - st[0]["Dss"]) < 1e-9 and abs(sum(x for _, x in cellP[2]) - st[0]["Dps"]) < 1e-9, "tip spans two stacks"
        ktS, ktP = kmax_of(st[0]["Dss"]), kmax_of(st[0]["Dps"])

        def tnode(side, j, k):
            sj = st[j]
            if side == "S":
                if j == 0:
                    return node_id[(col_key(0), k)]
                return node(("tipS", j, k), sj["P"] + sj["f"] * Dn[k] * sj["u"])
            if j == 0:
                return node_id[(col_key(nc - 1), k)]
            if sj["f"] < 1.0 and k == ktP:
                return tnode("S", j, ktS) if False else node(("tipS", j, ktS), sj["P"] + sj["f"] * Dn[ktS] * sj["u"])
            return node(("tipP", j, k), sj["Q"] - sj["f"] * Dn[k] * sj["u"])

        for j in range(len(st) - 1):
            for side, kt, cell, key in (("S", ktS, cellS, "P"), ("P", ktP, cellP, "Q")):
                a, b = st[j][key], st[j + 1][key]
                tv = b - a
                tv = tv / np.linalg.norm(tv)
                cum = np.cumsum([x for _, x in cell[2]])
                for k in range(kt):
                    mat = cell[2][int(np.searchsorted(cum, 0.5 * (Dn[k] + Dn[k + 1])))][0]
                    ids = [tnode(side, j, k), tnode(side, j + 1, k), tnode(side, j + 1, k + 1), tnode(side, j, k + 1)]
                    quads.append((ids, mat, tv, "tip"))
    # webs
    for w in sec.webs:
        cA = [int(np.argmin(np.abs(cols - sf))) for sf in w["foot"]["A"]]
        cB = [int(np.argmin(np.abs(cols - sf))) for sf in w["foot"]["B"]]
        for c, sf in zip(cA + cB, w["foot"]["A"] + w["foot"]["B"], strict=False):
            assert abs(cols[c] - sf) < 1e-9
        Q = len(w["lam"])
        endA = [node_id[(col_key(c), ktot_col[c])] for c in cA]
        endB = [node_id[(col_key(c), ktot_col[c])] for c in cB]
        FA = np.array([xy2[i] for i in endA])
        FB = np.array([xy2[i] for i in endB])
        J = max(1, math.ceil(float(np.linalg.norm(FB[0] - FA[0])) / h_s))
        grid = {}
        for q in range(Q):
            for j in range(J + 1):
                if j == 0:
                    grid[(j, q)] = endA[q]
                elif j == J:
                    grid[(j, q)] = endB[q]
                else:
                    grid[(j, q)] = node(("web", w["name"], j, q), FA[q] + (FB[q] - FA[q]) * j / J)
        for q in range(Q - 1):
            mat = w["layers"][w["lay_of"][q]][0]
            for j in range(J):
                quads.append(([grid[(j, q)], grid[(j + 1, q)], grid[(j + 1, q + 1)], grid[(j, q + 1)]], mat, w["u"], w["name"]))
    xy2 = np.array(xy2)
    # CCW enforcement
    out = []
    for ids, mat, tv, tag in quads:
        p = xy2[ids]
        a = 0.5 * np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
        if a < 0:
            ids = ids[::-1]
            a = -a
        assert a > 1e-14, f"degenerate quad area {a} tag {tag}"
        out.append((ids, mat, tv, a, tag))
    return {"xy": xy2, "quads": out, "node_sk": node_sk, "sec": sec, "s_lo": s_lo, "s_hi": s_hi, "tip_area": tip_area, "ncol": nc, "Dn": Dn}


# ------------------------------------------------------------------------------------------------ solid + input deck
def eng_const_line(m):
    E, nu, G = m["E"], m["nu"], m["G"]
    E1, E2, E3 = E
    S = np.zeros((6, 6))
    S[0, 0], S[1, 1], S[2, 2] = 1 / E1, 1 / E2, 1 / E3
    S[0, 1] = S[1, 0] = -nu[0] / E1
    S[0, 2] = S[2, 0] = -nu[1] / E1
    S[1, 2] = S[2, 1] = -nu[2] / E2
    S[3, 3], S[4, 4], S[5, 5] = 1 / G[0], 1 / G[1], 1 / G[2]
    ev = np.linalg.eigvalsh(S)
    assert ev.min() > 0, f"material compliance not positive definite: {m}"
    return (f"{E1:.8g}, {E2:.8g}, {E3:.8g}, {nu[0]:.8g}, {nu[1]:.8g}, {nu[2]:.8g}, {G[0]:.8g}, {G[1]:.8g}\n{G[2]:.8g}\n")


def write_model(sec, m2, out_dir, n_chord, hz, nmid, planes_frac):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    xy = m2["xy"]
    quads = m2["quads"]
    allxy = np.vstack([sec.P(np.linspace(0, 1, 2001))])
    chord = float(np.ptp(allxy[:, 0]) if sec.name.startswith("blade") else 2 * 0.5)  # tube: diameter
    if sec.name.startswith("blade"):
        chord = float(np.linalg.norm(sec.le_point - sec.te_point))
    L = n_chord * chord
    n1 = max(1, round(L / 3 / hz))
    n2 = nmid  # elements in the middle third (quadratic in z: the St-Venant fields there are polynomial of degree <= 2)
    zs = np.concatenate([np.linspace(0, L / 3, n1 + 1), np.linspace(L / 3, 2 * L / 3, n2 + 1)[1:], np.linspace(2 * L / 3, L, n1 + 1)[1:]])
    nz = len(zs) - 1
    N2, NQ = len(xy), len(quads)
    # in-plane edges
    edge_id = {}
    qedges = []
    for ids, *_ in quads:
        es = []
        for a, b in zip(ids, ids[1:] + ids[:1], strict=False):
            key = (min(a, b), max(a, b))
            if key not in edge_id:
                edge_id[key] = len(edge_id)
            es.append(edge_id[key])
        qedges.append(es)
    NE = len(edge_id)
    ekeys = sorted(edge_id, key=edge_id.get)
    base_e = N2 * (nz + 1)  # in-plane edge mid: 1 + base_e + k*NE + e
    base_v = base_e + NE * (nz + 1)  # vertical mid: 1 + base_v + k*N2 + i  (k = 0..nz-1)
    nnode = base_v + N2 * nz
    ref_id, rot_id = nnode + 1, nnode + 2
    corner = lambda i, k: 1 + k * N2 + i  # noqa: E731
    emid = lambda e, k: 1 + base_e + k * NE + e  # noqa: E731
    vmid = lambda i, k: 1 + base_v + k * N2 + i  # noqa: E731

    # lumped corner-node areas (for the plane fits)
    wnode = np.zeros(N2)
    for ids, _, _, a, _ in quads:
        for i in ids:
            wnode[i] += a / 4
    cx = float(np.sum(wnode * xy[:, 0]) / wnode.sum())
    cy = float(np.sum(wnode * xy[:, 1]) / wnode.sum())

    # orientations
    ori_key, ori_list = {}, []
    sets = {}
    for q, (_ids, mat, tv, _a, _tag) in enumerate(quads):
        key = (round(float(tv[0]), 7), round(float(tv[1]), 7))
        if key not in ori_key:
            ori_key[key] = len(ori_list)
            ori_list.append(key)
        sets.setdefault((mat, ori_key[key]), []).append(q)

    f = open(out / "mesh.inp", "w")
    f.write(f"** l3_solid_section_ccx {sec.name} NQ={NQ} N2={N2} nz={nz} L={L:.4f}\n*NODE\n")
    for k in range(nz + 1):
        for i in range(N2):
            f.write(f"{corner(i, k)},{xy[i,0]:.10g},{xy[i,1]:.10g},{zs[k]:.10g}\n")
    for k in range(nz + 1):
        for e, (a, b) in enumerate(ekeys):
            m = 0.5 * (xy[a] + xy[b])
            sa, sb = m2["node_sk"].get(a), m2["node_sk"].get(b)
            if sa and sb and abs(sa[1] - sb[1]) < 1e-12:  # along-wall edge: midside node on the true offset curve
                s_a, s_b = sa[0], sb[0]
                if abs(s_a - s_b) > 0.5:  # closed seam: the node is shared by s=0 and s=1
                    s_a, s_b = (s_a + 1, s_b) if s_a < s_b else (s_a, s_b + 1)
                sm = (0.5 * (s_a + s_b)) % 1.0
                m = sec.P(sm) + sa[1] * sec.normal_in(sm)
            f.write(f"{emid(e, k)},{m[0]:.10g},{m[1]:.10g},{zs[k]:.10g}\n")
    for k in range(nz):
        zm = 0.5 * (zs[k] + zs[k + 1])
        for i in range(N2):
            f.write(f"{vmid(i, k)},{xy[i,0]:.10g},{xy[i,1]:.10g},{zm:.10g}\n")
    f.write(f"{ref_id},{cx:.10g},{cy:.10g},{L:.10g}\n{rot_id},{cx:.10g},{cy:.10g},{L:.10g}\n")
    f.write("*ELEMENT, TYPE=C3D20R\n")
    for k in range(nz):
        for q, (ids, *_r) in enumerate(quads):
            eid = 1 + k * NQ + q
            e = qedges[q]
            n = ([corner(i, k) for i in ids] + [corner(i, k + 1) for i in ids]
                 + [emid(e[j], k) for j in range(4)] + [emid(e[j], k + 1) for j in range(4)]
                 + [vmid(i, k) for i in ids])
            # CCX order: 1-4 bottom, 5-8 top, 9-12 bottom mids, 13-16 top mids, 17-20 vertical mids
            f.write(f"{eid}, " + ", ".join(map(str, n[:15])) + ",\n" + ", ".join(map(str, n[15:])) + "\n")
    # sets and sections
    for (mat, oi), ql in sorted(sets.items()):
        name = f"S_{mat}_{oi}"
        ids = [1 + k * NQ + q for k in range(nz) for q in ql]
        f.write(f"*ELSET, ELSET={name}\n")
        for j in range(0, len(ids), 16):
            f.write(", ".join(map(str, ids[j:j + 16])) + ",\n")
    for mat, p in sec.materials.items():
        f.write(f"*MATERIAL, NAME={mat}\n")
        if np.ndim(p["E"]) == 0:
            f.write(f"*ELASTIC\n{p['E']:.8g}, {p['nu']:.8g}\n")
        else:
            f.write("*ELASTIC, TYPE=ENGINEERING CONSTANTS\n" + eng_const_line(p))
    for oi, (tx, ty) in enumerate(ori_list):
        f.write(f"*ORIENTATION, NAME=O{oi}, SYSTEM=RECTANGULAR\n0., 0., 1., {tx:.8f}, {ty:.8f}, 0.\n")
    for (mat, oi) in sorted(sets):
        f.write(f"*SOLID SECTION, ELSET=S_{mat}_{oi}, MATERIAL={mat}, ORIENTATION=O{oi}\n")
    # node sets
    def nset(name, ids):
        f.write(f"*NSET, NSET={name}\n")
        ids = list(ids)
        for j in range(0, len(ids), 16):
            f.write(", ".join(map(str, ids[j:j + 16])) + ",\n")

    def plane_nodes(k):
        return ([corner(i, k) for i in range(N2)] + [emid(e, k) for e in range(NE)])

    nset("ENDA", plane_nodes(0))
    nset("ENDB", plane_nodes(nz))
    levels = []
    for fr in planes_frac:
        k = int(np.argmin(np.abs(zs - fr * L)))
        levels.append(k)
    levels = sorted(set(levels))
    nset("PLANES", [corner(i, k) for k in levels for i in range(N2)])
    f.write(f"*RIGID BODY, NSET=ENDB, REF NODE={ref_id}, ROT NODE={rot_id}\n")
    f.close()
    cases = [("F", ref_id, 3, 1.0e6), ("T", rot_id, 3, 1.0e6), ("MX", rot_id, 1, 1.0e6), ("MY", rot_id, 2, 1.0e6)]
    for name, nd, dof, val in cases:  # one input per load case (the four run as parallel single-thread ccx processes)
        with open(out / f"case_{name}.inp", "w") as g:
            g.write(f"*INCLUDE, INPUT=mesh.inp\n*STEP\n*STATIC\n*BOUNDARY\nENDA, 1, 3, 0.\n*CLOAD\n{nd}, {dof}, {val:.8g}\n")
            g.write("*NODE PRINT, NSET=PLANES\nU\n*END STEP\n")
    meta = {"case": sec.name, "N2": N2, "NQ": NQ, "nz": nz, "nelem": NQ * nz, "nnode": nnode + 2, "L": L, "chord": chord, "zs_levels": [float(zs[k]) for k in levels],
                "levels": levels, "cx": cx, "cy": cy, "cases": [(c[0], c[2], c[3]) for c in cases], "s_lo": m2["s_lo"], "s_hi": m2["s_hi"],
                "tip_area": m2["tip_area"], "ncol": m2["ncol"], "depth_levels": len(m2["Dn"]), "area": float(wnode.sum()),
                "corner_xy": xy.tolist(), "wnode": wnode.tolist(), "n_orient": len(ori_list), "n_sets": len(sets), "N2_corner_id_base": 1}
    json.dump(meta, open(out / "meta.json", "w"))
    return meta


# ------------------------------------------------------------------------------------------------ analysis
def parse_dat(dat, nsteps=4):
    blocks, cur = [], None
    for line in open(dat):
        if "displacements" in line and "set" in line:
            cur = {}
            blocks.append(cur)
            continue
        if cur is not None:
            p = line.split()
            if len(p) == 4:
                try:
                    cur[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]))
                except ValueError:
                    pass
    assert len(blocks) == nsteps, f"expected {nsteps} displacement blocks, got {len(blocks)}"
    return blocks


def fit_plane(xy, w, ux, uy, uz, cx, cy):
    X = np.column_stack([np.ones(len(xy)), xy[:, 0] - cx, xy[:, 1] - cy])
    W = w[:, None]
    A = X.T @ (W * X)
    cw = np.linalg.solve(A, X.T @ (w * uz))  # a, b, c
    cu = np.linalg.solve(A, X.T @ (w * ux))
    cv = np.linalg.solve(A, X.T @ (w * uy))
    theta = 0.5 * (cv[1] - cu[2])
    return cw, theta


def strain_pair(meta, blk, k1, k2):
    xy = np.array(meta["corner_xy"])
    w = np.array(meta["wnode"])
    N2 = meta["N2"]
    res = []
    for k in (k1, k2):
        ids = [1 + k * N2 + i for i in range(N2)]
        U = np.array([blk[i] for i in ids])
        cw, th = fit_plane(xy, w, U[:, 0], U[:, 1], U[:, 2], meta["cx"], meta["cy"])
        res.append((cw, th))
    z1, z2 = meta["zs_k"][k1], meta["zs_k"][k2]
    dz = z2 - z1
    (c1, t1), (c2, t2) = res
    a_p = (c2[0] - c1[0]) / dz
    th_p = (t2 - t1) / dz
    kx = (c2[2] - c1[2]) / dz  # phi_x = c
    ky = -(c2[1] - c1[1]) / dz  # phi_y = -b
    return np.array([a_p, th_p, kx, ky])


def stiffness(meta, blocks, win):
    k1, k2 = win
    loads = np.array([c[2] for c in meta["cases"]])
    C = np.column_stack([strain_pair(meta, b, k1, k2) / l for b, l in zip(blocks, loads, strict=False)])  # strains per unit load
    sym = float(np.max(np.abs(C - C.T) / np.sqrt(np.outer(np.abs(np.diag(C)), np.abs(np.diag(C))))))
    Cs = 0.5 * (C + C.T)
    K = np.linalg.inv(Cs)
    return K, sym


def invariants(K4):
    """K4 ordered (F, T, Mx, My). Returns EA, GJ, bending eigenvalues after axial condensation."""
    EA, GJ = K4[0, 0], K4[1, 1]
    S = K4[np.ix_([0, 2, 3], [0, 2, 3])]
    B = S[1:, 1:] - np.outer(S[1:, 0], S[0, 1:]) / S[0, 0]
    ev = np.sort(np.linalg.eigvalsh(B))
    return {"EA": float(EA), "GJ": float(GJ), "EI_min": float(ev[0]), "EI_max": float(ev[1])}


def analyze(out_dir):
    out = Path(out_dir)
    meta = json.load(open(out / "meta.json"))
    zs_lv = meta["zs_levels"]
    levels = meta["levels"]
    meta["zs_k"] = dict(zip(levels, zs_lv, strict=False))
    blocks = []
    for c in meta["cases"]:
        b = parse_dat(out / f"case_{c[0]}.dat", 1)
        blocks.append(b[0])
    L = meta["L"]
    near = lambda fr: levels[int(np.argmin([abs(z - fr * L) for z in zs_lv]))]  # noqa: E731
    wins = {"mid_third": (near(1 / 3), near(2 / 3)), "mid_fifth": (near(0.4), near(0.6)), "mid_half": (near(0.25), near(0.75))}
    res = {}
    for nm, win in wins.items():
        K, sym = stiffness(meta, blocks, win)
        inv = invariants(K)
        inv.update(sym_err=sym, win_z=[meta["zs_k"][win[0]], meta["zs_k"][win[1]]], K4=K.tolist())
        res[nm] = inv
    res["meta"] = {k: v for k, v in meta.items() if k not in ("corner_xy", "wnode")}
    json.dump(res, open(out / "result.json", "w"), indent=1)
    return res


# ------------------------------------------------------------------------------------------------ 2D St-Venant cross-check
def torsion_2d(sec, m2):
    """Independent 2D St-Venant torsion (Prandtl warping function, bilinear quads, tangent/normal shear moduli G12/G13) on the SAME
    2D mesh and materials as the solid: checks the 3D path (rigid caps, extraction) against its own continuum limit. Returns GJ."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spl

    xy, quads = m2["xy"], m2["quads"]
    n = len(xy)
    gp = [(-1 / math.sqrt(3), 1 / math.sqrt(3))]
    pts = [(a, b) for a in gp[0] for b in gp[0]]
    I, J, V = [], [], []
    f = np.zeros(n)
    const = 0.0
    for ids, mat, tv, _area, _tag in quads:
        m = sec.materials[mat]
        G12, G13 = (m["G"][0], m["G"][1]) if np.ndim(m["E"]) else (m["G"], m["G"])
        tv = np.asarray(tv) / np.linalg.norm(tv)
        nv = np.array([-tv[1], tv[0]])
        Gt = G12 * np.outer(tv, tv) + G13 * np.outer(nv, nv)
        X = xy[ids]
        Ke = np.zeros((4, 4))
        fe = np.zeros(4)
        for a, b in pts:
            dN = 0.25 * np.array([[-(1 - b), (1 - b), (1 + b), -(1 + b)], [-(1 - a), -(1 + a), (1 + a), (1 - a)]])
            Jm = dN @ X
            dets = np.linalg.det(Jm)
            B = np.linalg.solve(Jm, dN)  # rows: d/dx, d/dy
            N = 0.25 * np.array([(1 - a) * (1 - b), (1 + a) * (1 - b), (1 + a) * (1 + b), (1 - a) * (1 + b)])
            xp, yp = N @ X[:, 0], N @ X[:, 1]
            g0 = np.array([-yp, xp])
            Ke += B.T @ Gt @ B * dets
            fe += -(B.T @ Gt @ g0) * dets
            const += g0 @ Gt @ g0 * dets
        for i in range(4):
            f[ids[i]] += fe[i]
            for j in range(4):
                I.append(ids[i])
                J.append(ids[j])
                V.append(Ke[i, j])
    K = sp.csr_matrix((V, (I, J)), shape=(n, n))
    # fix the rigid constant mode: pin node 0
    free = np.arange(1, n)
    phi = np.zeros(n)
    phi[free] = spl.spsolve(K[free][:, free].tocsc(), f[free])
    return float(const - f @ phi)  # GJ = int g0.G.g0 - f.phi, minimum of the warping functional


# ------------------------------------------------------------------------------------------------ deck comparison
def deck_station(eta, deck=DECK):
    lines = [l.split() for l in open(deck).read().split("DISTRIBUTED PROPERTIES")[1].split("\n")[1:]]
    nums = [l for l in lines if l]
    i = 0
    while i < len(nums):
        if len(nums[i]) == 1:
            e = float(nums[i][0])
            rows = nums[i + 1:i + 7]
            if abs(e - eta) < 1e-4 and all(len(r) == 6 for r in rows):
                return np.array(rows, float)
            i += 13
        else:
            i += 1
    raise KeyError(eta)


def deck_metrics(D):
    """Like cmp.py: GJ condensed over the two shear strains, EA raw, bending eigenvalues after condensing axial (raw 6x6).
    Plus the shear-free (zero shear force) variant K4_sf, the like-for-like of the solid."""
    cond = lambda K, t, s: K[t, t] - K[t, s] @ np.linalg.solve(K[np.ix_(s, s)], K[s, t])  # noqa: E731
    raw = {"GJ_sc": float(cond(D, 5, [0, 1])), "EA": float(D[2, 2])}
    S = D[np.ix_([2, 3, 4], [2, 3, 4])]
    B = S[1:, 1:] - np.outer(S[1:, 0], S[0, 1:]) / S[0, 0]
    ev = np.sort(np.linalg.eigvalsh(B))
    raw.update(EI_min=float(ev[0]), EI_max=float(ev[1]))
    idx = [2, 5, 3, 4]
    C6 = np.linalg.inv(D)
    K4 = np.linalg.inv(C6[np.ix_(idx, idx)])
    sf = invariants(K4)
    return raw, sf


# ------------------------------------------------------------------------------------------------ driver
def build(args):
    if args.case == "blade":
        sec = blade_section(args.eta)
    else:
        sec = tube_section(args.case)
    m2 = build_2d(sec, args.h, args.tmax, args.web_mode)
    meta = write_model(sec, m2, args.out, args.nchord, args.hz, args.nmid, [0.25, 1 / 3, 0.4, 0.5, 0.6, 2 / 3, 0.75])
    if args.st2d:
        meta["GJ_st2d"] = torsion_2d(sec, m2)
        json.dump(meta, open(Path(args.out) / "meta.json", "w"))
        print(f"2D St-Venant GJ on the same mesh: {meta['GJ_st2d']:.6e}")
    print(f"built {sec.name}: NQ={meta['NQ']} nz={meta['nz']} elements={meta['nelem']} nodes={meta['nnode']} L={meta['L']:.3f} "
          f"chord={meta['chord']:.3f} s_range=[{meta['s_lo']:.5f},{meta['s_hi']:.5f}] tip_area={meta['tip_area']:.3e}")
    return meta


def run_ccx(out_dir, threads=None):
    """The four load cases run as four parallel single-thread ccx processes (SPOOLES here is serial)."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    procs = []
    for c in ("F", "T", "MX", "MY"):
        lg = open(Path(out_dir) / f"ccx_{c}.log", "w")
        procs.append((c, subprocess.Popen([CCX, f"case_{c}"], cwd=out_dir, stdout=lg, stderr=subprocess.STDOUT, env=env)))
    bad = [c for c, p in procs if p.wait() != 0]
    if bad:
        raise SystemExit(f"ccx failed for load cases {bad}; see {out_dir}/ccx_*.log")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["build", "run", "analyze"])
    ap.add_argument("--case", default="blade", choices=["tube_iso", "tube_ortho", "tube_ortho_g13", "tube_web", "blade"])
    ap.add_argument("--eta", type=float, default=0.5)
    ap.add_argument("--h", type=float, default=0.1, help="in-plane (arc) element size, m")
    ap.add_argument("--tmax", type=float, default=0.03, help="max element height through a depth interval / web layer, m")
    ap.add_argument("--hz", type=float, default=0.3, help="span element size in the end thirds, m")
    ap.add_argument("--nmid", type=int, default=5, help="span elements in the middle third (5 puts levels at 0.4 L and 0.6 L)")
    ap.add_argument("--st2d", action="store_true", help="also run the 2D St-Venant torsion cross-check on the same mesh")
    ap.add_argument("--nchord", type=float, default=4.0)
    ap.add_argument("--web-mode", dest="web_mode", default="bonded", choices=["bonded", "through"],
                    help="bonded: web ends on the inner skin surface (cap stays whole); through: web (OML chart line) replaces the skin under its footprint")
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_NTASKS", "4")))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.cmd in ("build", "run"):
        build(args)
    if args.cmd == "run":
        run_ccx(args.out, args.threads)
    if args.cmd in ("run", "analyze"):
        res = analyze(args.out)
        r = res["mid_third"]
        print({k: r[k] for k in ("EA", "GJ", "EI_min", "EI_max", "sym_err")})
        if args.case == "blade":
            raw, sf = deck_metrics(deck_station(args.eta))
            print("deck raw", raw, "\ndeck shear-free", {k: sf[k] for k in ('EA', 'GJ', 'EI_min', 'EI_max')})
            print("solid/deck(raw): GJ %.3f EA %.3f flap %.3f edge %.3f" % (r["GJ"] / raw["GJ_sc"], r["EA"] / raw["EA"],
                  r["EI_min"] / raw["EI_min"], r["EI_max"] / raw["EI_max"]))
        else:
            ex = tube_exact(args.case)
            print("exact", ex)
            print({k: (r[k] / v - 1 if v else None) for k, v in ex.items()})


if __name__ == "__main__":
    main()
