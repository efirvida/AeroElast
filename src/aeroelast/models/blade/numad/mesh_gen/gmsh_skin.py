"""Gmsh skin mesher PROTOTYPE for the windIO v2 blade (skin only, no webs, not wired into generate_mesh).

Pipeline (see logs/check/gmsh_proto/NOTE.md for the measured outcome)

1. chart      x = u = nd_arc * P(eta)   (P = loft perimeter), y = s = eta * L   (L = reference span length).
              The skin is a ring, so the chart is a strip whose two vertical sides (nd = 0 and nd = 1) are the TE seam.
2. cleaning   ``clean_source``: the source event lines (laminate / web-attachment boundaries nd = f(eta)) are examined in
              pairs of neighbours.  Where two neighbours are closer than ``tau * h`` they are replaced by ONE line
              (``mode='mid'``: their mean; ``'snap'``: the lower one; a fixed line - seam, web attachment - always
              wins).  The thin band is closed by a short chart-transverse step at the station where the separation
              reaches ``tau * h``, so no band thinner than ``tau * h`` reaches the mesher.  An open stretch shorter than
              ``tau * h`` between two merged stretches of the same pair is merged too.  The report counts the bands,
              their chart area and the largest displacement of a material boundary.  Nothing is repaired silently.
3. gmsh       the cleaned linework is noded with shapely, simplified (Douglas-Peucker) and meshed with the built-in
              ``geo`` kernel: boundary curves + embedded event curves (conforming mesh edges), a background size field
              and Frontal-Delaunay for quads with recombination; the seam is made periodic.
4. 3D         nodes are evaluated on the loft (``RowLoft``: ``Geometry.generate`` is ~0.08 s per station whatever the
              batch size, so a few hundred stations are cubic-interpolated across eta at fixed nd).
5. gate       ``element_metrics`` / ``quality_gate`` REJECT (``MeshQualityError``) - they never repair.

The module name does not start with ``windio2_`` on purpose: the on-disk mesh cache key of tests/support/v2_blade.py
hashes ``windio2_*.py`` only.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.spatial import cKDTree

__all__ = ["generate_mesh_gmsh", "map_mesh_materials_gmsh", "le_landmark_line", "lookup_cells", "areal_mass", "conformity_report", "Line", "SkinSource", "Cleaned", "ChartMesh", "MeshQualityError", "RowLoft", "clean_source",
           "build_linework", "mesh_chart", "mesh_chart_structured", "add_le_lines", "map_to_3d", "element_metrics", "quality_gate", "generate_skin_mesh",
           "source_from_topology", "ring_size", "unique_edges", "mesh_web_chart"]

EPS_GAP = 1e-12  # nd tolerance for "equal"
EPS_G = 1e-9  # chart-metre tolerance of the gap-vs-threshold sign (stations are rounded to 1e-13 in eta ~ 1e-11 m)


class MeshQualityError(RuntimeError):
    """The gate found elements outside the accepted limits; ``.report`` carries the numbers."""

    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report


# --------------------------------------------------------------------------------------------------- source model
@dataclass
class Line:
    """Event line nd(eta), piecewise linear on ``grid``; drawn on ``[lo, hi]``. ``fixed`` lines never move."""

    name: str
    grid: np.ndarray
    values: np.ndarray
    lo: float = 0.0
    hi: float = 1.0
    fixed: bool = False

    def at(self, eta):
        return np.interp(eta, self.grid, self.values)


@dataclass
class SkinSource:
    lines: list
    length: float  # reference span length L [m]
    perimeter: object  # eta -> P(eta) [m], vectorised
    drops: list = field(default_factory=list)  # (eta, line_a, line_b, 'left'|'right'): layer drop-off steps
    loft: object = None  # (eta array, nd array) -> (n, 3)
    eta0: float = 0.0  # chart window [eta0, eta1] (a patch of the blade); y = (eta - eta0) * L
    eta1: float = 1.0
    extra_stations: list = field(default_factory=list)  # eta stations that must be rows (topology events: web/layer spans)


@dataclass
class Cleaned:
    stations: np.ndarray
    polylines: dict  # name -> (m, 2) [eta, nd]
    closures: list  # [(m, 2) [eta, nd]] drop-off steps
    report: dict
    row_limits: object = None  # None, or per station interval the longest row step [m] the kept strips allow (inf = free)


def ring_size(h, perimeter, min_ring_divisions=0):
    """Local element size h(eta) = min(h, P(eta)/min_ring_divisions) (0 disables), the guard the v2 mesher has too."""
    if not min_ring_divisions:
        return lambda eta: np.full_like(np.asarray(eta, float), h)
    return lambda eta: np.minimum(h, np.asarray(perimeter(eta), float) / min_ring_divisions)


# --------------------------------------------------------------------------------------------------- 2D cleaning
def _interval_rows(a, b, size, eta0, length, max_step=np.inf):
    """Rows (chart y, ``a`` excluded, ``b`` included) of the station interval [a, b]: equal steps of the integral of
    1/size, one count per interval (raised until no step exceeds ``max_step``, the bound a kept narrow strip imposes)."""
    yy = np.linspace(a, b, max(2, int((b - a) / 0.05) + 2))
    inv = 1.0 / np.asarray(size(eta0 + yy / length), float)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(yy))])
    m = max(1, int(round(cum[-1])))
    while True:
        rows = np.interp(np.linspace(0.0, cum[-1], m + 1), cum, yy)
        if not np.isfinite(max_step) or np.max(np.diff(rows)) <= max_step * (1.0 + 1e-9):
            return rows[1:]
        m += 1


def _uniq(values):
    return np.unique(np.round(np.asarray(sorted(values), float), 13))


def clean_source(src, h, tau=0.25, mode="mid", size=None, max_station_gap=1.0, snap=0.0, max_strip_ar=5.0, min_row_fraction=0.5):
    """Merge event lines closer than ``tau * size(eta)`` (chart metres). Returns :class:`Cleaned`.

    ``max_strip_ar`` (default 5 = half the gate limit ``DEFAULT_LIMITS['max_ar']``; None/0 = merge every thin band, the
    previous behaviour) keeps a thin stretch of a pair of lines as a one-element-wide STRIP and bounds its aspect ratio by
    shortening the rows over it: every row of the stretch is at most ``max_strip_ar * gap`` long (``Cleaned.row_limits``,
    honoured by ``mesh_chart_structured``).  The stretch is kept only while that costs at most a factor 1/``min_row_fraction``
    (default 0.5) in row density, i.e. ``max_strip_ar * gap >= min_row_fraction * size`` on every interval; a stretch that
    gets narrower anywhere (a wedge closing to zero, a pinch) is merged as before, over its whole thin extent.

    ``snap`` > 0 coalesces stations closer than ``snap * size(eta)`` (chart metres) into the lowest one of the chain
    (the first/last station and the layer drop-off stations are kept).  Event stations that sit a few millimetres apart
    (closing steps of several bands at the ring-size cap, source vertices) otherwise force 1e-2 m elements, AR 10-18."""
    if mode not in ("mid", "snap"):
        raise ValueError("mode must be 'mid' or 'snap'")
    L, P = src.length, src.perimeter
    e0, e1 = src.eta0, src.eta1

    def uniq(values):
        a = _uniq(values)
        return a[(a >= e0 - 1e-13) & (a <= e1 + 1e-13)]

    size = size or (lambda e: np.full_like(np.asarray(e, float), h))

    def thr(e):
        return tau * np.asarray(size(e), float)

    lines = [Line("seam0", np.array([0., 1.]), np.array([0., 0.]), 0., 1., True),
             Line("seam1", np.array([0., 1.]), np.array([1., 1.]), 0., 1., True), *src.lines]
    index = {ln.name: i for i, ln in enumerate(lines)}
    n = len(lines)

    def matrix(stations):
        out = np.full((n, len(stations)), np.nan)
        for k, ln in enumerate(lines):
            alive = (stations >= ln.lo - 1e-12) & (stations <= ln.hi + 1e-12)
            out[k, alive] = np.interp(stations[alive], ln.grid, ln.values)
        return out

    s = {e0, e1}
    for ln in lines:
        s.update(float(g) for g in ln.grid if ln.lo < g < ln.hi)
        s.update((ln.lo, ln.hi))
    s.update(float(d[0]) for d in src.drops)
    s.update(float(e) for e in src.extra_stations)
    S = uniq(s)
    extra = []
    for a, b in zip(S[:-1], S[1:], strict=False):
        m = int(math.ceil((b - a) * L / max_station_gap))
        extra.extend(np.linspace(a, b, m + 1)[1:-1] if m > 1 else [])
    S = uniq([*S, *extra])

    # crossings of any two lines are stations of their own
    M = matrix(S)
    crossings = []
    for k in range(n - 1):
        D = M[k][None, :] - M[k + 1:]
        sg = np.where(np.abs(D) < EPS_GAP, 0., np.sign(D))
        with np.errstate(invalid="ignore"):
            prod = sg[:, :-1] * sg[:, 1:]
        for a, j in zip(*np.nonzero(prod < 0), strict=False):
            d0, d1 = D[a, j], D[a, j + 1]
            crossings.append(S[j] + (S[j + 1] - S[j]) * d0 / (d0 - d1))
    S = uniq([*S, *crossings])

    # the station where a neighbour gap reaches tau*h (and sign dips inside an interval)
    for _ in range(8):
        M = matrix(S)
        added = []
        mid = 0.5 * (S[:-1] + S[1:])
        Pl, Pm, Pr = P(S[:-1]), P(mid), P(S[1:])
        Tl, Tm, Tr = thr(S[:-1]), thr(mid), thr(S[1:])
        for j in range(len(S) - 1):
            alive = np.nonzero(np.isfinite(M[:, j]) & np.isfinite(M[:, j + 1]))[0]
            if len(alive) < 2:
                continue
            order = alive[np.argsort(M[alive, j] + M[alive, j + 1], kind="stable")]
            nl, nr = M[order, j], M[order, j + 1]
            nm = 0.5 * (nl + nr)
            gl, gm, gr = np.diff(nl) * Pl[j] - Tl[j], np.diff(nm) * Pm[j] - Tm[j], np.diff(nr) * Pr[j] - Tr[j]
            for i in range(len(gl)):
                a, b = order[i], order[i + 1]
                x, y, z = gl[i], gm[i], gr[i]
                x, y, z = (0. if abs(v) < EPS_G else v for v in (x, y, z))
                if x * z < 0:
                    root = brentq(lambda e: (lines[b].at(e) - lines[a].at(e)) * float(P(e)) - float(thr(e)),
                                  S[j], S[j + 1], xtol=1e-15)
                    added.append(root)
                elif x * y < 0 or z * y < 0:
                    added.append(mid[j])
        if not added:
            break
        S = uniq([*S, *added])
    else:
        raise RuntimeError(f"clean_source: station refinement did not converge (last additions {sorted(added)[:6]}, "
                           f"{len(added)} of them)")

    snapped = 0
    if snap > 0 and len(S) > 2:
        hard = {round(float(e0), 13), round(float(e1), 13), *(round(float(d[0]), 13) for d in src.drops),
                *(round(float(e), 13) for e in src.extra_stations)}
        keep = [S[0]]
        for v in S[1:]:
            close = (v - keep[-1]) * L < snap * float(size(keep[-1]))
            if close and round(float(v), 13) not in hard:
                snapped += 1
                continue
            if close and round(float(keep[-1]), 13) not in hard and len(keep) > 1:  # a hard station replaces the last one
                keep.pop()
                snapped += 1
            keep.append(v)
        S = np.asarray(keep)

    M = matrix(S)
    m_int = len(S) - 1
    mid = 0.5 * (S[:-1] + S[1:])
    Pe, Pm = P(S), P(mid)
    Tm = thr(mid)
    orders, thin = [], []
    for j in range(m_int):
        alive = np.nonzero(np.isfinite(M[:, j]) & np.isfinite(M[:, j + 1]))[0]
        order = alive[np.argsort(M[alive, j] + M[alive, j + 1], kind="stable")]
        nm = 0.5 * (M[order, j] + M[order, j + 1])
        orders.append(order)
        thin.append(list((np.diff(nm) * Pm[j] - Tm[j]) < -EPS_G) if len(order) > 1 else [])

    # per-pair runs of consecutive intervals; a short open lens between two thin stretches is thin as well
    def pair_runs(flag_of):
        runs = {}
        for j in range(m_int):
            for i, flag in enumerate(flag_of(j)):
                key = (orders[j][i], orders[j][i + 1])
                lst = runs.setdefault(key, [])
                if lst and lst[-1][-1][0] == j - 1:
                    lst[-1].append((j, flag))
                else:
                    lst.append([(j, flag)])
        return runs

    flips = 0
    for key, chains in pair_runs(lambda j: thin[j]).items():
        for chain in chains:
            j = 0
            while j < len(chain):
                if chain[j][1]:
                    j += 1
                    continue
                e = j
                while e < len(chain) and not chain[e][1]:
                    e += 1
                if 0 < j and e < len(chain):
                    a, b = chain[j][0], chain[e - 1][0]
                    if (S[b + 1] - S[a]) * L < tau * float(np.min(size(mid[a:b + 1]))):
                        for jj, _ in chain[j:e]:
                            thin[jj][list(zip(orders[jj][:-1], orders[jj][1:], strict=False)).index(key)] = True
                            flips += 1
                j = e

    # a thin stretch whose gap stays above size / max_strip_ar is a persistent band: keep it as a strip
    kept, row_limits = [], np.full(m_int, np.inf)
    if max_strip_ar:
        for key, chains in pair_runs(lambda j: thin[j]).items():
            a, b = key
            for chain in chains:
                groups, cur = [], []
                for jj, flag in chain:
                    if flag:
                        cur.append(jj)
                    elif cur:
                        groups.append(cur)
                        cur = []
                if cur:
                    groups.append(cur)
                for g in groups:
                    gl_ = (M[b, g] - M[a, g]) * Pe[g]
                    gr_ = (M[b, [x + 1 for x in g]] - M[a, [x + 1 for x in g]]) * Pe[[x + 1 for x in g]]
                    allowed = max_strip_ar * np.minimum(gl_, gr_)  # longest row step that keeps the strip's AR <= max_strip_ar
                    if np.all(allowed >= min_row_fraction * np.asarray(size(mid[g]), float) - EPS_G):
                        for jj in g:
                            thin[jj][list(zip(orders[jj][:-1], orders[jj][1:], strict=False)).index(key)] = False
                        row_limits[g] = np.minimum(row_limits[g], allowed)
                        kept.append({"lines": [lines[a].name, lines[b].name], "s_from_m": float(S[g[0]] * L),
                                     "s_to_m": float(S[g[-1] + 1] * L), "min_gap_m": float(min(gl_.min(), gr_.min())),
                                     "max_gap_m": float(max(gl_.max(), gr_.max())), "n_intervals": len(g)})

    # clusters and replacement values at both ends of every interval
    W, conflicts = {}, []
    for j in range(m_int):
        order = orders[j]
        i = 0
        while i < len(order):
            e = i
            while e < len(thin[j]) and thin[j][e]:
                e += 1
            comp = list(order[i:e + 1])
            if len(comp) > 1:
                fixed = [k for k in comp if lines[k].fixed]
                if len({round(M[k, j], 12) for k in fixed} | {round(M[k, j + 1], 12) for k in fixed}) > 2:
                    conflicts.append({"interval": [float(S[j]), float(S[j + 1])],
                                      "lines": [lines[k].name for k in comp]})
                else:
                    ends = []
                    for col in (j, j + 1):
                        if fixed:
                            ends.append(M[fixed[0], col])
                        elif mode == "mid":
                            ends.append(float(np.mean(M[comp, col])))
                        else:
                            ends.append(M[comp[0], col])
                    for k in comp:
                        W[(k, j)] = (float(ends[0]), float(ends[1]))
            i = e + 1

    # report: bands, area, displacement
    bands, coincident = [], 0
    for key, chains in pair_runs(lambda j: thin[j]).items():
        a, b = key
        for chain in chains:
            js = [jj for jj, flag in chain if flag]
            groups, cur = [], []
            for jj in js:
                if cur and jj != cur[-1] + 1:
                    groups.append(cur)
                    cur = []
                cur.append(jj)
            if cur:
                groups.append(cur)
            for g in groups:
                gl = (M[b, g] - M[a, g]) * Pe[g]
                gr = (M[b, [x + 1 for x in g]] - M[a, [x + 1 for x in g]]) * Pe[[x + 1 for x in g]]
                ds = (S[[x + 1 for x in g]] - S[g]) * L
                area = float(np.sum(0.5 * (gl + gr) * ds))
                disp = 0.0
                for k in (a, b):
                    if lines[k].fixed:
                        continue
                    for jj in g:
                        wl, wr = W.get((k, jj), (M[k, jj], M[k, jj + 1]))
                        disp = max(disp, abs(wl - M[k, jj]) * Pe[jj], abs(wr - M[k, jj + 1]) * Pe[jj + 1])
                rec = {"lines": [lines[a].name, lines[b].name], "s_from_m": float(S[g[0]] * L),
                       "s_to_m": float(S[g[-1] + 1] * L), "max_gap_m": float(max(gl.max(), gr.max())),
                       "area_m2": area, "max_displacement_m": float(disp)}
                if rec["max_gap_m"] < 1e-9:
                    coincident += 1
                else:
                    bands.append(rec)
    max_disp = 0.0
    for (k, j), (wl, wr) in W.items():
        if not lines[k].fixed:
            max_disp = max(max_disp, abs(wl - M[k, j]) * Pe[j], abs(wr - M[k, j + 1]) * Pe[j + 1])

    # polylines per drawn line
    polylines = {}
    for k, ln in enumerate(lines):
        if k < 2:
            continue
        pts = []
        for j in range(m_int):
            if not (np.isfinite(M[k, j]) and np.isfinite(M[k, j + 1])):
                continue
            wl, wr = W.get((k, j), (M[k, j], M[k, j + 1]))
            for p in ((S[j], wl), (S[j + 1], wr)):
                if not pts or abs(pts[-1][0] - p[0]) > 1e-14 or abs(pts[-1][1] - p[1]) > 1e-14:
                    pts.append(p)
        if len(pts) > 1:
            polylines[ln.name] = np.asarray(pts)

    closures = []
    for eta_d, na, nb, side in src.drops:
        j = int(np.searchsorted(S, round(eta_d, 13)))
        jj, col = (j - 1, 1) if side == "left" else (j, 0)
        vals = []
        for name in (na, nb):
            k = index[name]
            vals.append(W.get((k, jj), (M[k, jj], M[k, jj + 1]))[col])
        if abs(vals[0] - vals[1]) * float(P(eta_d)) > 1e-9:
            closures.append(np.array([[eta_d, vals[0]], [eta_d, vals[1]]]))

    report = {
        "tau": tau, "h": h, "mode": mode, "snap": snap, "n_stations_snapped": snapped, "n_stations": int(len(S)), "n_lines": n - 2,
        "max_strip_ar": max_strip_ar, "min_row_fraction": min_row_fraction, "n_bands_kept": len(kept), "kept": kept,
        "n_bands_merged": len(bands), "n_coincident_runs": coincident, "n_short_lenses_merged": flips,
        "total_band_chart_area_m2": float(sum(b["area_m2"] for b in bands)),
        "max_band_gap_m": float(max((b["max_gap_m"] for b in bands), default=0.0)),
        "max_material_boundary_displacement_m": float(max_disp),
        "unresolved_fixed_conflicts": conflicts, "bands": bands,
    }
    return Cleaned(S, polylines, closures, report, row_limits if np.isfinite(row_limits).any() else None)


# --------------------------------------------------------------------------------------------------- linework
def _chart_xy(src, eta, nd):
    eta = np.asarray(eta, float)
    return np.column_stack([np.asarray(nd, float) * src.perimeter(eta), (eta - src.eta0) * src.length])


def build_linework(src, cleaned, snap_tol=1e-7, simplify_tol=0.005):
    """Cleaned lines -> noded, simplified chart edges (list of (m, 2) arrays). Seam/root/tip pieces are dropped
    here (the boundary is rebuilt in ``mesh_chart``), interior edges that END on the boundary are kept."""
    import shapely
    from shapely.geometry import LineString

    L, P = src.length, src.perimeter
    Ly = (src.eta1 - src.eta0) * L
    pls = [_chart_xy(src, pl[:, 0], pl[:, 1]) for pl in cleaned.polylines.values()]
    pls += [_chart_xy(src, c[:, 0], c[:, 1]) for c in cleaned.closures]

    def on(p, kind):
        x, y = p
        return {"l": abs(x) < snap_tol, "r": abs(x - float(P(src.eta0 + y / L))) < snap_tol,
                "b": abs(y) < snap_tol, "t": abs(y - Ly) < snap_tol}[kind]

    runs = []
    for pl in pls:
        cur = [pl[0]]
        for a, b in zip(pl[:-1], pl[1:], strict=False):
            if any(on(a, k) and on(b, k) for k in "lrbt") or np.hypot(*(b - a)) < snap_tol:
                if len(cur) > 1:
                    runs.append(np.asarray(cur))
                cur = [b]
            else:
                cur.append(b)
        if len(cur) > 1:
            runs.append(np.asarray(cur))
    if not runs:
        return []
    allp = np.vstack(runs)
    tree = cKDTree(allp)
    parent = np.arange(len(allp))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in tree.query_pairs(snap_tol):
        parent[find(max(i, j))] = find(min(i, j))
    rep = np.array([find(i) for i in range(len(allp))])
    snapped = allp[rep]
    out, pos = [], 0
    for r in runs:
        q = snapped[pos:pos + len(r)]
        pos += len(r)
        keep = [0, *[i for i in range(1, len(q)) if np.hypot(*(q[i] - q[i - 1])) > 0]]
        if len(keep) > 1:
            out.append(LineString(q[keep]))
    merged = shapely.unary_union(out)
    geoms = list(merged.geoms) if hasattr(merged, "geoms") else [merged]
    edges = []
    for g in geoms:
        g = g.simplify(simplify_tol, preserve_topology=False) if simplify_tol else g
        c = np.asarray(g.coords)
        if len(c) >= 2 and np.hypot(*(c[-1] - c[0])) > 0:
            edges.append(c)
    return edges


# --------------------------------------------------------------------------------------------------- gmsh
@dataclass
class ChartMesh:
    chart: np.ndarray  # (n, 2) x = u, y = s   (unmerged: the right seam keeps its own nodes)
    elements: np.ndarray  # (ne, 4), -1 pads triangles; positive (eta, nd) orientation
    master: np.ndarray  # (n,) index of the node each node is identified with (seam), else itself
    seam_pairs: np.ndarray  # (k, 2) [slave, master]
    info: dict


def mesh_chart(src, edges, h, size=None, algorithm=8, recombination=1, optimize_topology=5, tol=1e-7,
               recombine=True, extra_options=None, parity=0, optimize=(), optimize_iter=1, field_data=None):
    """Mesh the chart with gmsh (geo kernel). Deterministic: one thread, fixed seed."""
    import gmsh

    Lspan, P = src.length, src.perimeter
    L = (src.eta1 - src.eta0) * Lspan  # chart height
    size = size or (lambda e: np.full_like(np.asarray(e, float), h))
    Pmax = float(np.max(P(np.linspace(0, 1, 400))))

    def key(x, y):
        return (int(round(x / tol)), int(round(y / tol)))

    ys_seam, xs_root, xs_tip = {0.0, L}, {0.0, float(P(src.eta0))}, {0.0, float(P(src.eta1))}
    for e in edges:
        for p in (e[0], e[-1]):
            x, y = p
            if abs(x) < tol or abs(x - float(P(src.eta0 + y / Lspan))) < tol:
                ys_seam.add(float(y))
            if abs(y) < tol:
                xs_root.add(float(x))
            if abs(y - L) < tol:
                xs_tip.add(float(x))
    junction_y = [0.0]
    for y in sorted(ys_seam):
        if y - junction_y[-1] > tol:
            junction_y.append(y)

    def seam_ys(a, b):
        """Node heights of the seam between two junctions: equal steps of the integral of 1/size (one element each)."""
        yy = np.linspace(a, b, max(2, int((b - a) / 0.05) + 2))
        inv = 1.0 / np.asarray(size(src.eta0 + yy / Lspan), float)
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(yy))])
        n = max(1, int(round(cum[-1])))
        return np.interp(np.linspace(0.0, cum[-1], n + 1), cum, yy)

    gmsh.initialize()
    t0 = time.perf_counter()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.logger.start()
        gmsh.option.setNumber("General.NumThreads", 1)
        for opt, val in (("Mesh.MaxNumThreads1D", 1), ("Mesh.MaxNumThreads2D", 1), ("Mesh.MaxNumThreads3D", 1),
                         ("Mesh.RandomSeed", 1), ("Mesh.RandomFactor", 1e-9)):
            try:
                gmsh.option.setNumber(opt, val)
            except Exception:  # option name differs between gmsh versions
                pass
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(np.max(size(src.eta0 + np.linspace(0, L, 801) / Lspan))))
        gmsh.option.setNumber("Mesh.Algorithm", algorithm)
        gmsh.option.setNumber("Mesh.RecombineAll", 1 if recombine else 0)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", recombination)
        gmsh.option.setNumber("Mesh.RecombineOptimizeTopology", optimize_topology)
        for opt, val in (extra_options or {}).items():  # last word: any gmsh option, e.g. Mesh.Smoothing
            gmsh.option.setNumber(opt, val)
        gmsh.model.add("skin")
        geo = gmsh.model.geo
        points = {}

        def pt(x, y):
            k = key(x, y)
            if k not in points:
                points[k] = geo.addPoint(x, y, 0.0, h)
            return points[k]

        def curve(c):
            ids = [pt(x, y) for x, y in c]
            return geo.addLine(ids[0], ids[1])

        def chain(c):
            # straight segments, never a spline: an event polyline has corners (closure steps, kinks) and a spline
            # through a corner overshoots across the boundary
            return [curve(c[i:i + 2]) for i in range(len(c) - 1)]

        def pieces(vals, make):
            return [make(a, b) for a, b in zip(vals[:-1], vals[1:], strict=False)]

        root = pieces(sorted(xs_root), lambda a, b: curve(np.array([[a, 0.0], [b, 0.0]])))
        tip = pieces(sorted(xs_tip), lambda a, b: curve(np.array([[a, L], [b, L]])))
        # The period P(eta) varies along the span, so a gmsh periodic constraint (one affine map) cannot express the
        # seam.  Both seams get the SAME node heights explicitly: one element per piece, identical y on both sides.
        left, right = [], []
        for a, b in zip(junction_y[:-1], junction_y[1:], strict=False):
            yk = seam_ys(a, b)
            for y0, y1 in zip(yk[:-1], yk[1:], strict=False):
                left.append(curve(np.array([[0.0, y0], [0.0, y1]])))
                right.append(curve(np.array([[float(P(src.eta0 + y0 / Lspan)), y0], [float(P(src.eta0 + y1 / Lspan)), y1]])))
        loop = [*root, *right, *[-t for t in reversed(tip)], *[-t for t in reversed(left)]]
        cl = geo.addCurveLoop(loop)
        surf = geo.addPlaneSurface([cl])
        interior = [t for e in edges for t in chain(e)]
        geo.synchronize()
        if interior:
            gmsh.model.mesh.embed(1, interior, 2, surf)
        if parity:
            # full-quad recombination (RecombinationAlgorithm 2/3) meshes coarse and subdivides: gmsh requires an EVEN
            # number of 1D elements on EVERY curve ("1D mesh cannot be divided by 2" otherwise).  The final element
            # size is half the size field, so pass size = 2h.  Seam pieces get exactly 2 divisions (mirrored on both
            # sides); root/tip pieces and embedded event edges get an even count from their length.
            def ncount(c):
                (x0, y0), (x1, y1) = (np.asarray(gmsh.model.getValue(0, p_, [])[:2]) for p_ in
                                       [v for _, v in gmsh.model.getBoundary([(1, c)], oriented=False)])
                n = math.hypot(x1 - x0, y1 - y0) / float(size(src.eta0 + 0.5 * (y0 + y1) / Lspan))
                return 2 * max(1, int(round(n / 2)))
            for c in (*left, *right):
                gmsh.model.mesh.setTransfiniteCurve(c, 3)
            for c in (*root, *tip, *interior):
                gmsh.model.mesh.setTransfiniteCurve(c, ncount(c) + 1)
        else:
            for c in (*left, *right):
                gmsh.model.mesh.setTransfiniteCurve(c, 2)
        ys_grid = np.linspace(0, L, 801)
        sizes = np.asarray(size(src.eta0 + ys_grid / Lspan), float)
        if field_data is not None:  # per-region size field (see ``region_size_data``)
            view = gmsh.view.add("size")
            gmsh.view.addListData(view, "ST", len(field_data) // 12, field_data)
            f = gmsh.model.mesh.field.add("PostView")
            gmsh.model.mesh.field.setNumber(f, "ViewTag", view)
        elif np.ptp(sizes) < 1e-12:
            f = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(f, "F", repr(float(sizes[0])))
        else:
            view = gmsh.view.add("size")
            data = []
            x0, x1 = -1.0, Pmax + 1.0
            for i in range(len(ys_grid) - 1):
                ya, yb, va, vb = ys_grid[i], ys_grid[i + 1], sizes[i], sizes[i + 1]
                data += [x0, x1, x1, ya, ya, yb, 0, 0, 0, va, va, vb]
                data += [x0, x1, x0, ya, yb, yb, 0, 0, 0, va, vb, vb]
            gmsh.view.addListData(view, "ST", 2 * (len(ys_grid) - 1), data)
            f = gmsh.model.mesh.field.add("PostView")
            gmsh.model.mesh.field.setNumber(f, "ViewTag", view)
        gmsh.model.mesh.field.setAsBackgroundMesh(f)
        gmsh.model.mesh.generate(1)
        gmsh.model.mesh.generate(2)
        for method in optimize:  # post-optimisation of the 2D mesh: 'Laplace2D', 'Relocate2D', 'Netgen', ...
            gmsh.model.mesh.optimize(method, niter=int(optimize_iter))

        tags, coords, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        tags, xyz = tags[order], coords.reshape(-1, 3)[order]
        pos = {int(t): i for i, t in enumerate(tags)}
        elements = []
        for etype, k in ((2, 3), (3, 4)):
            _, nodes = gmsh.model.mesh.getElementsByType(etype)
            for row in np.asarray(nodes, dtype=np.int64).reshape(-1, k):
                elements.append([pos[int(t)] for t in row] + [-1] * (4 - k))
        elements = np.asarray(elements, dtype=int).reshape(-1, 4)
        if len(elements) == 0:
            tail = [m for m in gmsh.logger.get() if "Error" in m or "Warning" in m][:12]
            raise MeshQualityError("gmsh produced no 2D elements: " + " | ".join(tail))

        def seam_nodes(curves):
            ids = set()
            for c in curves:
                t, _, _ = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)
                ids.update(pos[int(v)] for v in t)
            return sorted(ids)

        ln, rn = seam_nodes(left), seam_nodes(right)
        chart = xyz[:, :2].copy()
    finally:
        gmsh.finalize()
    elapsed = time.perf_counter() - t0

    # gmsh also meshes the control points of splines as isolated nodes: keep only nodes an element uses
    used = np.unique(elements[elements >= 0])
    compact = np.full(len(chart), -1, dtype=int)
    compact[used] = np.arange(len(used))
    elements = np.where(elements >= 0, compact[np.where(elements >= 0, elements, 0)], -1)
    chart = chart[used]
    ln = [int(compact[i]) for i in ln if compact[i] >= 0]
    rn = [int(compact[i]) for i in rn if compact[i] >= 0]
    lo = np.array(ln)[np.argsort(chart[ln, 1], kind="stable")]
    ro = np.array(rn)[np.argsort(chart[rn, 1], kind="stable")]
    if len(lo) != len(ro):
        yl = chart[lo, 1]
        miss_r = [(round(float(chart[i, 1]), 5), round(float(chart[i, 0] - float(P(src.eta0 + chart[i, 1] / Lspan))), 6))
                  for i in ro if np.min(np.abs(yl - chart[i, 1])) > 1e-9][:8]
        yr = chart[ro, 1]
        miss_l = [round(float(chart[i, 1]), 5) for i in lo if np.min(np.abs(yr - chart[i, 1])) > 1e-9][:8]
        raise MeshQualityError(f"periodic seam not mirrored: {len(lo)} vs {len(ro)} nodes; right-only (y, x-P) {miss_r}; "
                               f"left-only y {miss_l}")
    # gmsh mirrors the 1D mesh by curve parameter; the right seam is a curve x = P(y), so its nodes sit a hair off the
    # master's y.  Slide each slave along the seam to the master's y (reported; refused when it is not a hair).
    dy = chart[ro, 1] - chart[lo, 1]
    seam_dy = float(np.max(np.abs(dy))) if len(dy) else 0.0
    if seam_dy > 0.1 * h:
        raise MeshQualityError(f"periodic seam nodes differ in y by {seam_dy:.4g} m (> 0.1 h)")
    chart[ro, 1] = chart[lo, 1]
    chart[ro, 0] = np.asarray(P(src.eta0 + chart[lo, 1] / Lspan), float)
    master = np.arange(len(chart))
    master[ro] = lo
    # orientation: positive in (eta, nd) == clockwise in (x = u, y = s)
    for i, e in enumerate(elements):
        ids = e[e >= 0]
        a = 0.5 * np.sum(chart[ids, 0] * np.roll(chart[ids, 1], -1) - np.roll(chart[ids, 0], -1) * chart[ids, 1])
        if a > 0:
            elements[i, :len(ids)] = ids[::-1]
    return ChartMesh(chart, elements, master, np.column_stack([ro, lo]),
                     {"gmsh_seconds": elapsed, "n_embedded_edges": len(edges), "algorithm": algorithm,
                      "recombination": recombination, "seam_nodes": int(len(lo)), "seam_slide_max_m": seam_dy})


def region_size_data(src, cleaned, h, size=None, dy=0.25, lo=0.8, hi=1.25):
    """Background size field for gmsh built from the STRIP WIDTHS: in each region between two adjacent cleaned lines
    the target is ``w / n`` with ``n = max(1, round(w / size))`` (clipped to [lo, hi] * size when n >= 2), so that an
    integer number of elements fits the region and the chordwise width tracks the spanwise length.  Returns the flat
    triangle list for a gmsh 'ST' post-processing view (x, y in chart metres)."""
    L, P = src.length, src.perimeter
    size = size or (lambda e: np.full_like(np.asarray(e, float), h))
    ys = [(cleaned.stations - src.eta0) * L]
    Y = np.unique(np.round(np.concatenate([ys[0], np.arange(0.0, (src.eta1 - src.eta0) * L, dy), [(src.eta1 - src.eta0) * L]]), 9))
    pls = list(cleaned.polylines.values())

    def xs_at(eta, tiny):
        """Sorted distinct line x at ``eta`` (nd * P), including the seam lines, using the polylines alive across it."""
        out = [0.0, float(P(eta))]
        for pl in pls:
            if pl[0, 0] - 1e-12 <= eta + tiny and pl[-1, 0] + 1e-12 >= eta + tiny:
                out.append(float(np.interp(eta + tiny, pl[:, 0], pl[:, 1])) * float(P(eta)))
        return np.sort(np.asarray(out))

    data = []
    for ya, yb in zip(Y[:-1], Y[1:], strict=False):
        ea, eb = src.eta0 + ya / L, src.eta0 + yb / L
        em = 0.5 * (ea + eb)
        tiny = (eb - ea) * 1e-6
        xa, xb = xs_at(ea, tiny), xs_at(eb, -tiny)
        xm = xs_at(em, 0.0)
        if not (len(xa) == len(xb) == len(xm)):  # topology changes inside the interval: fall back to the plain size
            xa = xb = xm
        else:
            pass
        for i in range(len(xm) - 1):
            wa, wb = xa[i + 1] - xa[i], xb[i + 1] - xb[i]
            if wa < 1e-9 and wb < 1e-9:
                continue

            def tgt(w, e):
                sz = float(size(e))
                n = max(1, int(round(w / sz)))
                return float(np.clip(w / n, lo * sz, hi * sz)) if n >= 2 else float(np.clip(w, lo * sz, hi * sz))
            ha, hb = tgt(max(wa, 1e-9), ea), tgt(max(wb, 1e-9), eb)
            c = [(xa[i], ya), (xa[i + 1], ya), (xb[i + 1], yb), (xb[i], yb)]
            v = [ha, ha, hb, hb]
            for tri in ((0, 1, 2), (0, 2, 3)):
                data += [c[t][0] for t in tri] + [c[t][1] for t in tri] + [0.0, 0.0, 0.0] + [v[t] for t in tri]
    return data


# --------------------------------------------------------------------------------------------------- structured chart
@dataclass
class _Cell:
    k: int  # row interval [rows[k], rows[k+1]]
    a: int  # left / right line index
    b: int
    xb0: float
    xb1: float
    xt0: float
    xt1: float
    h: float
    n: int = 1
    patch: bool = False


def _line_row_values(pl, eta_rows):
    """nd of a cleaned polyline at the global rows.  A closing/opening STEP (two points with the same eta) is replaced
    by a slant over the next row interval: at the station row the value BELOW the step wins, the next row takes the
    value of the piece above, so no transverse edge shorter than a row exists (the material boundary moves by at most
    the step, i.e. tau*h/2, and only over that one row)."""
    pieces = [pc for pc in np.split(pl, np.nonzero(np.abs(np.diff(pl[:, 0])) < 1e-13)[0] + 1) if len(pc) >= 2]
    out = np.full(len(eta_rows), np.nan)
    for i, pc in enumerate(pieces):
        lo, hi = pc[0, 0], pc[-1, 0]
        sel = (eta_rows >= lo - 1e-12) & (eta_rows <= hi + 1e-12)
        v = np.interp(eta_rows[sel], pc[:, 0], pc[:, 1])
        if i > 0:
            v[np.abs(eta_rows[sel] - lo) < 1e-12] = pieces[i - 1][-1, 1]
        out[sel] = v
    return out


def mesh_chart_structured(src, cleaned, h, size=None, patch_rows=2, band=(0.8, 1.25), algorithm=8, recombination=1,
                          optimize_topology=5, extra_options=None, tol=1e-7, optimize=(), optimize_iter=3, strip_scale=None, optimize_scope="patches", rows=None, active=None, periodic=True, count_rule="round", row_blocks=False):
    """Structured chart mesh with gmsh as the engine.

    1. rows: the event stations plus uniform subdivisions with one count per interval (equal steps of the integral of
       1/size, as the v2 mesher does); every nd line is straight between two rows.
    2. cells: (strip between two consecutive distinct lines) x (row interval).  Chordwise count
       n = max(1, round(width/size)) per strip, held (hysteresis ``band`` on width/(n size)) while the strip is
       regular, so n changes only when the element width leaves the band.
    3. REGULAR blocks: maximal runs of cells with the same strip and count -> one transfinite recombined surface whose
       nodes lie on the rows.  Shared edges carry one count (a conflict raises ``MeshQualityError``).
    4. TRANSITION patches: where the count changes, a strip appears/disappears or two lines converge, the cells within
       ``patch_rows`` rows of the mismatch form a small unstructured Frontal-Delaunay + recombination patch whose
       boundary divisions are imposed by the neighbouring blocks.
    5. the two TE seams use the same node heights (the rows).

    ``optimize`` runs gmsh's own 2D optimisers (``gmsh.model.mesh.optimize``) after meshing on ``optimize_scope``:
    'patches' (default), 'neighbors' (patches + the blocks sharing a curve with one) or 'all'.  Nodes classified on
    curves (event lines, seams, root/tip rows, embedded lines) are not moved by a surface optimiser; the displacement of
    every node and of the curve nodes is measured and reported in ``info['optimize_report']``.

    ``count_rule`` 'proportional': a strip is divided into n EQUAL parts (a transfinite division of the strip's
    own width, so the grid lines follow the strip when its width varies along the span); n is the integer whose cell
    width stays in ``band`` * h and is closest to h (floor when the cell would drop below the band, ceil when it would
    exceed it), judged also by the cell aspect ratio against the row length, keeping the previous count when it is still
    in the band or tied (fewer count changes).  'round' (default): n = round(width / h) outside the band.  Measured on
    the blade (h 0.25 / 0.15) both give the same uniformity (83.3 vs 83.1 % / 89.5 vs 89.4 % of cells within 20 % of h)
    and 'proportional' has slightly more AR > 5 (7 vs 5 / 13 vs 10), so 'round' stays the default.

    ``row_blocks`` (default False): one transfinite surface per cell (one row high) instead of per run of equal cells, so
    every chordwise grid line is EXACTLY a row (gmsh's Coons map of a multi-row block between slanted sides tilts the
    interior grid lines, up to 0.17 m off the row on the blade).  Measured on the blade it costs uniformity (cells within
    20 % of h: 78.6 -> 73.2 % at h 0.25, 88.5 -> 86.3 % at 0.15; same q120, AR, triangles), so it is opt-in.

    ``rows`` (chart y of every row, default: built from the stations) lets a second chart (a web) share the rows of
    the skin; ``active`` (bool per row interval) restricts the cells to part of the span; ``periodic=False`` skips the
    TE seam identification (a web has two distinct edges).

    ``strip_scale`` {(line_a, line_b): factor} scales the chordwise target of one strip (e.g. a leading-edge strip
    bounded by two mesh-control lines, see ``add_le_lines``); the spanwise row spacing is global."""
    import gmsh
    from shapely.geometry import Point, Polygon
    from shapely.ops import unary_union

    Lspan, P = src.length, src.perimeter
    Ly = (src.eta1 - src.eta0) * Lspan
    size = size or (lambda e: np.full_like(np.asarray(e, float), h))
    t0 = time.perf_counter()

    # ---- rows
    if rows is None:
        st = np.clip((cleaned.stations - src.eta0) * Lspan, 0.0, Ly)
        rows = [st[0]]
        lim = getattr(cleaned, "row_limits", None)
        for i_, (a, b) in enumerate(zip(st[:-1], st[1:], strict=False)):
            rows.extend(_interval_rows(a, b, size, src.eta0, Lspan, np.inf if lim is None else lim[i_]))
    rows = np.asarray(rows)
    nr = len(rows)
    eta_rows = src.eta0 + rows / Lspan
    Pr = np.asarray(P(eta_rows), float)

    # ---- line table  X[l, k]  (chart x), seams first
    names = ["seam0", "seam1", *cleaned.polylines.keys()]
    X = np.full((len(names), nr), np.nan)
    X[0], X[1] = 0.0, Pr
    for i, nm in enumerate(names[2:], 2):
        X[i] = _line_row_values(cleaned.polylines[nm], eta_rows) * Pr
    srcline = {ln.name: ln for ln in src.lines}

    # ---- cells
    cells = {}
    by_k = [[] for _ in range(nr - 1)]
    for k in range(nr - 1):
        if active is not None and not active[k]:
            continue
        alive = [i for i in range(len(names)) if np.isfinite(X[i, k]) and np.isfinite(X[i, k + 1])]
        emid = src.eta0 + 0.5 * (rows[k] + rows[k + 1]) / Lspan
        # order by the SOURCE nd at the interval middle (ties between merged lines keep the source order)
        def okey(i, emid=emid):
            if i == 0:
                return (-1e9, 0)
            if i == 1:
                return (1e9, 0)
            return (float(srcline[names[i]].at(emid)), i)
        alive.sort(key=okey)
        xm = 0.5 * (X[:, k] + X[:, k + 1])
        if any(xm[b] - xm[a] < -1e-9 for a, b in zip(alive[:-1], alive[1:], strict=False)):
            # the cleaned lines swapped order inside the interval (a closing step displaced one of two near-coincident
            # source lines past the other): the CLEANED geometry decides, the source order only breaks ties
            alive.sort(key=lambda i: (float(xm[i]), okey(i)))
        for a, b in zip(alive[:-1], alive[1:], strict=False):
            wb, wt = X[b, k] - X[a, k], X[b, k + 1] - X[a, k + 1]
            if max(abs(wb), abs(wt)) < 1e-9:
                continue
            if min(wb, wt) < -1e-9:
                raise MeshQualityError(f"structured mesher: lines {names[a]} / {names[b]} cross inside row interval {k} "
                                       f"(widths {wb:.3g}, {wt:.3g} m)")
            c = _Cell(k, a, b, X[a, k], X[b, k], X[a, k + 1], X[b, k + 1],
                      float(size(emid)) * float((strip_scale or {}).get((names[a], names[b]), 1.0)))
            c.patch = min(wb, wt) < 1e-9  # a pinched cell (a triangle) is never a regular block
            cells[(k, a, b)] = c
            by_k[k].append(c)

    # ---- chordwise count with hysteresis, per strip run
    runs = {}
    for (_k, a, b), c in sorted(cells.items()):
        runs.setdefault((a, b), []).append(c)
    for lst in runs.values():
        prev = None
        for c in lst:
            w = 0.5 * ((c.xb1 - c.xb0) + (c.xt1 - c.xt0))
            before = prev.n if (prev is not None and prev.k == c.k - 1) else None
            if before is not None and band[0] * c.h <= w / before <= band[1] * c.h:
                c.n = before
            elif count_rule == "round":
                c.n = max(1, int(round(w / c.h)))
            else:
                rl = float(rows[c.k + 1] - rows[c.k])
                lo_n = max(1, int(math.floor(w / (band[1] * c.h))))
                best = None
                for n_ in range(lo_n, int(math.ceil(w / (band[0] * c.h))) + 2):
                    ew = w / n_
                    key = (not band[0] * c.h <= ew <= band[1] * c.h,  # inside the band first
                           round(abs(math.log(ew / c.h)), 3),  # closest to h
                           round(max(ew, rl) / min(ew, rl), 3),  # cell aspect against the row length
                           n_ != before)  # tie: keep the previous count
                    if best is None or key < best[0]:
                        best = (key, n_)
                c.n = best[1]
            prev = c

    # ---- layer drop-off closures: transverse segments on a row (a material step).  The strip set changes there, so the
    #      cells on both sides become one transition patch and the segment is embedded in it as a chain of mesh edges
    #      (split at every line vertex that lies on it), see ``closure_chains`` / ``closure_report``.
    closures = []
    for cl in cleaned.closures:
        cl = np.asarray(cl, float)
        if len(cl) != 2 or abs(cl[0, 0] - cl[1, 0]) > 1e-12:
            raise NotImplementedError("a layer drop-off closure must be one transverse (constant-eta) segment")
        y = (cl[0, 0] - src.eta0) * Lspan
        kc = int(np.argmin(np.abs(rows - y)))
        if abs(rows[kc] - y) > 1e-6:
            raise MeshQualityError(f"structured mesher: drop-off closure at y={y:.9g} m is not on a row (nearest {rows[kc]:.9g})")
        xa, xb = sorted(float(v) * Pr[kc] for v in cl[:, 1])
        closures.append({"k": kc, "eta": float(cl[0, 0]), "y": float(rows[kc]), "x0": xa, "x1": xb})

    # ---- mismatches -> patch cells
    def overlap(c1, c2):
        return min(c1.xt1, c2.xb1) - max(c1.xt0, c2.xb0)

    seeds = []
    for cl in closures:  # cells on both sides of the segment (a boundary row needs no patch: its edges are boundary edges)
        c = cl["k"]
        if 1 <= c <= nr - 2:
            seeds.extend((c, u, "up") for u in by_k[c] if min(u.xb1, cl["x1"]) - max(u.xb0, cl["x0"]) > 1e-9)
            seeds.extend((c, l, "low") for l in by_k[c - 1] if min(l.xt1, cl["x1"]) - max(l.xt0, cl["x0"]) > 1e-9)
    for c in range(1, nr - 1):  # row boundary c between interval c-1 (lower) and c (upper)
        low, up = by_k[c - 1], by_k[c]
        upmap = {(u.a, u.b): u for u in up}
        for l in low:
            u = upmap.get((l.a, l.b))
            same = (u is not None and u.n == l.n and abs(l.xt0 - u.xb0) < 1e-9 and abs(l.xt1 - u.xb1) < 1e-9)
            if not same:
                seeds.append((c, l, "low"))
                seeds.extend((c, v, "up") for v in up if overlap(l, v) > 1e-9)
        lowmap = {(l.a, l.b) for l in low}
        for u in up:
            if (u.a, u.b) not in lowmap:
                seeds.append((c, u, "up"))
                seeds.extend((c, l, "low") for l in low if overlap(l, u) > 1e-9)
    for c, cell, side in seeds:
        rng = range(c - patch_rows, c) if side == "low" else range(c, c + patch_rows)
        for k in rng:
            hit = cells.get((k, cell.a, cell.b))
            if hit is not None:
                hit.patch = True
        cell.patch = True

    # ---- components of patch cells (vertical overlap)
    pcells = [c for c in cells.values() if c.patch]
    parent = {id(c): id(c) for c in pcells}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for c in range(1, nr - 1):
        low = [x for x in by_k[c - 1] if x.patch]
        up = [x for x in by_k[c] if x.patch]
        for l in low:
            for u in up:
                if overlap(l, u) > 1e-9 or (l.a, l.b) == (u.a, u.b):
                    parent[find(id(l))] = find(id(u))
    comps = {}
    for c in pcells:
        comps.setdefault(find(id(c)), []).append(c)

    # ---- gmsh geometry
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.logger.start()
        gmsh.option.setNumber("General.NumThreads", 1)
        for opt, val in (("Mesh.MaxNumThreads1D", 1), ("Mesh.MaxNumThreads2D", 1), ("Mesh.RandomSeed", 1), ("Mesh.RandomFactor", 1e-9)):
            try:
                gmsh.option.setNumber(opt, val)
            except Exception:
                pass
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 1)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.Algorithm", algorithm)
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", recombination)
        gmsh.option.setNumber("Mesh.RecombineOptimizeTopology", optimize_topology)
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(h))
        for opt, val in (extra_options or {}).items():
            gmsh.option.setNumber(opt, val)
        gmsh.model.add("skin_structured")
        geo = gmsh.model.geo
        pts, curves, _counts = {}, {}, {}
        conflicts = []

        def pt(x, y, lc):
            key = (int(round(x / tol)), int(round(y / tol)))
            if key not in pts:
                pts[key] = geo.addPoint(x, y, 0.0, lc)
            return pts[key]

        def seg(p0, p1, n):
            """signed curve tag p0 -> p1 with ``n`` divisions; a second request must agree (shared-edge assertion)."""
            key = (min(p0, p1), max(p0, p1))
            if key not in curves:
                curves[key] = (geo.addLine(*key), n)
                geo.mesh.setTransfiniteCurve(curves[key][0], n + 1)
            elif curves[key][1] != n:
                conflicts.append({"points": key, "counts": [curves[key][1], n]})
            tag = curves[key][0]
            return tag if p0 < p1 else -tag

        def hrow(x0, x1, k, n):  # horizontal edge on row k (0..nr-1)
            lc = float(size(eta_rows[k]))
            return seg(pt(x0, rows[k], lc), pt(x1, rows[k], lc), n)

        line_pts = {}

        def vseg(i, k):  # line i, row k -> k+1, upward
            p0, p1 = pt(X[i, k], rows[k], float(size(eta_rows[k]))), pt(X[i, k + 1], rows[k + 1], float(size(eta_rows[k + 1])))
            line_pts[(i, k)], line_pts[(i, k + 1)] = p0, p1
            return seg(p0, p1, 1)

        seam_curves = {0: [], 1: []}
        if periodic:
            for i in (0, 1):
                for k in range(nr - 1):
                    seam_curves[i].append(abs(vseg(i, k)))

        # blocks: maximal runs of non-patch cells of one strip with the same count
        blocks, _nsurf_block = [], 0
        for lst in runs.values():
            cur = []
            for c in lst + [None]:
                if c is not None and not c.patch and (not cur or (c.k == cur[-1].k + 1 and c.n == cur[-1].n)):
                    cur.append(c)
                    continue
                if cur:
                    blocks.extend([[b_] for b_ in cur] if row_blocks else [cur])
                cur = [c] if (c is not None and not c.patch) else []
        surfaces = []  # (tag, kind)
        surf_curves = {}
        for blk in blocks:
            c0, c1 = blk[0], blk[-1]
            n = c0.n
            bottom = hrow(c0.xb0, c0.xb1, c0.k, n)
            top = hrow(c1.xt0, c1.xt1, c1.k + 1, n)
            right = [vseg(c0.b, c.k) for c in blk]
            left = [vseg(c0.a, c.k) for c in blk]
            loop = geo.addCurveLoop([bottom, *right, -top, *[-t for t in reversed(left)]])
            s = geo.addPlaneSurface([loop])
            lc0 = float(size(eta_rows[c0.k]))
            corners = [pt(c0.xb0, rows[c0.k], lc0), pt(c0.xb1, rows[c0.k], lc0),
                       pt(c1.xt1, rows[c1.k + 1], lc0), pt(c1.xt0, rows[c1.k + 1], lc0)]
            geo.mesh.setTransfiniteSurface(s, "Left", corners)
            geo.mesh.setRecombine(2, s)
            surfaces.append((s, "block", []))
            surf_curves[s] = {abs(t) for t in (bottom, top, *right, *left)}

        # patches
        n_patch_bad = 0
        for members in comps.values():
            polys = [Polygon(_dedupe([(c.xb0, rows[c.k]), (c.xb1, rows[c.k]), (c.xt1, rows[c.k + 1]), (c.xt0, rows[c.k + 1])]))
                     for c in members]
            U = unary_union([p.buffer(0) for p in polys])
            geoms = list(U.geoms) if hasattr(U, "geoms") else [U]
            for poly in geoms:
                if poly.area <= 0:
                    continue
                if len(poly.interiors):
                    n_patch_bad += 1
                    raise MeshQualityError("structured mesher: a transition patch has a hole")
                bnd = poly.boundary
                mine = [c for c, pg in zip(members, polys, strict=False) if poly.buffer(1e-8).contains(pg.representative_point())]
                edges_b, edges_i = [], []
                for c in mine:
                    yb, yt = rows[c.k], rows[c.k + 1]
                    cand = [((c.xb0, yb), (c.xb1, yb), "h", c.n), ((c.xt0, yt), (c.xt1, yt), "h", c.n),
                            ((c.xb0, yb), (c.xt0, yt), "v", c.a), ((c.xb1, yb), (c.xt1, yt), "v", c.b)]
                    for p, q, kind, tag in cand:
                        if math.hypot(q[0] - p[0], q[1] - p[1]) < 1e-9:
                            continue
                        mid = Point(0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]))
                        (edges_b if bnd.distance(mid) < 1e-8 else edges_i).append((p, q, kind, tag, c))
                # boundary loop
                lcp = float(size(eta_rows[mine[0].k]))
                # chain the undirected boundary edges into one loop
                und = {}
                for p, q, kind, tag, c in edges_b:
                    ip, iq = pt(p[0], p[1], lcp), pt(q[0], q[1], lcp)
                    nn = tag if kind == "h" else 1
                    if kind == "v":  # bottom / top vertex of line ``tag`` at rows k / k + 1
                        line_pts[(tag, c.k)], line_pts[(tag, c.k + 1)] = ip, iq
                    und.setdefault(ip, []).append((iq, ip, iq, nn))
                    und.setdefault(iq, []).append((ip, ip, iq, nn))
                used, loop_tags = set(), []
                start = next(iter(und))
                cur = start
                while True:
                    nxt = next(((o, a_, b_, nn) for (o, a_, b_, nn) in und[cur] if (min(a_, b_), max(a_, b_)) not in used), None)
                    if nxt is None:
                        break
                    o, a_, b_, nn = nxt
                    used.add((min(a_, b_), max(a_, b_)))
                    loop_tags.append(seg(cur, o, nn))
                    cur = o
                    if cur == start:
                        break
                if len(used) != len({(min(a_, b_), max(a_, b_)) for lst in und.values() for (_, a_, b_, _) in lst}):
                    raise MeshQualityError("structured mesher: a transition patch boundary is not a single loop")
                loop = geo.addCurveLoop(loop_tags)
                s = geo.addPlaneSurface([loop])
                geo.mesh.setRecombine(2, s)
                emb = []
                seen = set()
                for p, q, kind, tag, c in edges_i:
                    if kind != "v":
                        continue
                    ip, iq = pt(p[0], p[1], lcp), pt(q[0], q[1], lcp)
                    line_pts[(tag, c.k)], line_pts[(tag, c.k + 1)] = ip, iq
                    key = (min(ip, iq), max(ip, iq))
                    if key in seen:
                        continue
                    seen.add(key)
                    emb.append(abs(seg(ip, iq, 1)))
                for cl in closures:
                    xs = {cl["x0"], cl["x1"]}
                    for i in range(len(names)):  # every line vertex on the row inside the segment splits it
                        if np.isfinite(X[i, cl["k"]]) and cl["x0"] - 1e-7 <= X[i, cl["k"]] <= cl["x1"] + 1e-7:
                            xs.add(float(X[i, cl["k"]]))
                    xs = sorted(xs)
                    chain = [xs[0]]
                    for v in xs[1:]:
                        if v - chain[-1] > 1e-7:
                            chain.append(v)
                    for xp, xq in zip(chain[:-1], chain[1:], strict=False):
                        mid = Point(0.5 * (xp + xq), cl["y"])
                        if not poly.buffer(1e-8).contains(mid) or bnd.distance(mid) < 1e-8:
                            continue  # outside this patch, or already a boundary edge of it
                        ip, iq = pt(xp, cl["y"], lcp), pt(xq, cl["y"], lcp)
                        key = (min(ip, iq), max(ip, iq))
                        if key in seen or key in curves:
                            continue
                        seen.add(key)
                        emb.append(abs(seg(ip, iq, max(1, int(round((xq - xp) / float(size(eta_rows[cl["k"]]))))))))
                surfaces.append((s, "patch", emb))
                surf_curves[s] = {abs(t) for t in loop_tags}
        if conflicts:
            raise MeshQualityError(f"structured mesher: {len(conflicts)} shared edges with conflicting counts, e.g. {conflicts[:3]}")
        geo.synchronize()
        for s, _kind, emb in surfaces:
            if emb:
                gmsh.model.mesh.embed(1, emb, 2, s)
        gmsh.model.mesh.generate(2)
        patch_tags = [(2, s) for s, kind, _ in surfaces if kind == "patch"]
        opt_report = {}
        if optimize and patch_tags:
            t_opt = time.perf_counter()
            if optimize_scope == "patches":
                scope = list(patch_tags)
            elif optimize_scope == "all":
                scope = [(2, s) for s, _, _ in surfaces]
            elif optimize_scope == "neighbors":
                pc = set().union(*(surf_curves[s] for _, s in patch_tags))
                scope = [(2, s) for s, kind, _ in surfaces if kind == "patch" or (surf_curves[s] & pc)]
            else:
                raise ValueError("optimize_scope must be 'patches', 'neighbors' or 'all'")
            t0n, c0n, _ = gmsh.model.mesh.getNodes()
            before = dict(zip(t0n.tolist(), c0n.reshape(-1, 3), strict=False))
            curve_nodes = set()
            for _, c in gmsh.model.getEntities(1):
                tt, _, _ = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)
                curve_nodes.update(int(v) for v in tt)
            for method in optimize:  # 'Relocate2D', 'Laplace2D', 'Netgen', 'QuadQuasiStructured', ...
                gmsh.model.mesh.optimize(method, niter=int(optimize_iter), dimTags=scope)
            t1n, c1n, _ = gmsh.model.mesh.getNodes()
            disp = {int(t): float(np.linalg.norm(c - before[int(t)])) for t, c in zip(t1n, c1n.reshape(-1, 3), strict=False) if int(t) in before}
            opt_report = {"methods": list(optimize), "scope": optimize_scope, "iter": int(optimize_iter), "surfaces": len(scope),
                          "seconds": time.perf_counter() - t_opt,
                          "max_node_displacement_m": max(disp.values(), default=0.0),
                          "max_curve_node_displacement_m": max((v for k, v in disp.items() if k in curve_nodes), default=0.0),
                          "n_nodes_moved": int(sum(v > 1e-12 for v in disp.values())), "nodes_added": int(len(t1n) - len(t0n))}

        tags, coords, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        tags, xyz = tags[order], coords.reshape(-1, 3)[order]
        pos = {int(t): i for i, t in enumerate(tags)}
        elements, patch_tri, patch_el, in_patch = [], 0, 0, []
        for s, kind, _ in surfaces:
            for etype, kk in ((2, 3), (3, 4)):
                _, nodes = gmsh.model.mesh.getElementsByType(etype, tag=s)
                arr = np.asarray(nodes, dtype=np.int64).reshape(-1, kk)
                for row in arr:
                    elements.append([pos[int(t)] for t in row] + [-1] * (4 - kk))
                    in_patch.append(kind == 'patch')
                if kind == "patch":
                    patch_el += len(arr)
                    patch_tri += len(arr) if kk == 3 else 0
        elements = np.asarray(elements, dtype=int).reshape(-1, 4)
        if len(elements) == 0:
            tail = [m for m in gmsh.logger.get() if "Error" in m or "Warning" in m][:12]
            raise MeshQualityError("gmsh produced no 2D elements: " + " | ".join(tail))

        def seam_nodes(ts):
            ids = set()
            for c in ts:
                t, _, _ = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)
                ids.update(pos[int(v)] for v in t)
            return sorted(ids)

        ln, rn = (seam_nodes(seam_curves[0]), seam_nodes(seam_curves[1])) if periodic else ([], [])
        chart = xyz[:, :2].copy()
        vertex_node = {}
        for key, ptag in line_pts.items():
            vt, _, _ = gmsh.model.mesh.getNodes(0, ptag)
            if len(vt):
                vertex_node[key] = pos[int(vt[0])]
    finally:
        gmsh.finalize()
    elapsed = time.perf_counter() - t0

    used = np.unique(elements[elements >= 0])
    compact = np.full(len(chart), -1, dtype=int)
    compact[used] = np.arange(len(used))
    elements = np.where(elements >= 0, compact[np.where(elements >= 0, elements, 0)], -1)
    chart = chart[used]
    ln = [int(compact[i]) for i in ln if compact[i] >= 0]
    rn = [int(compact[i]) for i in rn if compact[i] >= 0]
    lo = np.array(ln, dtype=int)[np.argsort(chart[ln, 1], kind="stable")] if ln else np.zeros(0, dtype=int)
    ro = np.array(rn, dtype=int)[np.argsort(chart[rn, 1], kind="stable")] if rn else np.zeros(0, dtype=int)
    if len(lo) != len(ro) or not np.allclose(chart[lo, 1], chart[ro, 1], atol=1e-9):
        raise MeshQualityError(f"structured seam: {len(lo)} vs {len(ro)} nodes or different heights")
    master = np.arange(len(chart))
    master[ro] = lo
    line_row_nodes = {}
    for (i, k), nd_ in vertex_node.items():
        if compact[nd_] >= 0:
            line_row_nodes.setdefault(names[i], np.full(nr, -1, dtype=int))[k] = int(compact[nd_])
    for i, e in enumerate(elements):
        ids = e[e >= 0]
        a = 0.5 * np.sum(chart[ids, 0] * np.roll(chart[ids, 1], -1) - np.roll(chart[ids, 0], -1) * chart[ids, 1])
        if a > 0:
            elements[i, :len(ids)] = ids[::-1]
    closure_info = closure_report(chart, elements, closures)
    bad = [c for c in closure_info if not c["covered"]]
    if bad:
        raise MeshQualityError(f"structured mesher: {len(bad)} drop-off closure(s) are not a chain of mesh edges: {bad[:3]}")
    ncells = len(cells)
    info = {"closures": closure_info, "gmsh_seconds": elapsed, "structured": True, "n_rows": int(nr), "n_cells": ncells,
            "n_patch_cells": len(pcells), "patch_cell_fraction": len(pcells) / max(ncells, 1), "n_blocks": len(blocks),
            "in_patch": np.asarray(in_patch), "optimize_report": opt_report, "rows": rows, "line_row_nodes": line_row_nodes,
            "line_names": names, "n_patches": int(sum(1 for s in surfaces if s[1] == "patch")), "patch_elements": int(patch_el),
            "patch_triangles": int(patch_tri), "count_conflicts": 0, "seam_nodes": int(len(lo)),
            "n_count_changes": int(sum(1 for lst in runs.values() for c0, c1 in zip(lst[:-1], lst[1:], strict=False) if c1.k == c0.k + 1 and c1.n != c0.n))}
    return ChartMesh(chart, elements, master, np.column_stack([ro, lo]), info)


def closure_report(chart, elements, closures, tol=1e-6):
    """For every closure ``{y, x0, x1}`` (chart metres): the chart nodes on the segment, sorted by x, and whether
    consecutive ones are joined by a mesh edge from ``x0`` to ``x1`` (then no element can straddle the step)."""
    if not closures:
        return []
    edges = {(int(a), int(b)) for a, b in unique_edges(elements)}
    out = []
    for cl in closures:
        on = np.nonzero((np.abs(chart[:, 1] - cl["y"]) < tol) & (chart[:, 0] > cl["x0"] - tol) & (chart[:, 0] < cl["x1"] + tol))[0]
        on = on[np.argsort(chart[on, 0], kind="stable")]
        joined = all((min(int(a), int(b)), max(int(a), int(b))) in edges for a, b in zip(on[:-1], on[1:], strict=False))
        ends = len(on) >= 2 and abs(chart[on[0], 0] - cl["x0"]) < tol and abs(chart[on[-1], 0] - cl["x1"]) < tol
        out.append({"eta": cl["eta"], "y": cl["y"], "x0": cl["x0"], "x1": cl["x1"], "n_nodes": int(len(on)),
                    "covered": bool(joined and ends)})
    return out


def add_le_lines(src, nd_le=0.4968, half_width=0.15, grid_step=1.0):
    """Two fixed mesh-control lines ``LE-`` / ``LE+`` at +- ``half_width`` metres (chart) from the leading edge, so the
    structured mesher can give the LE strip its own, finer chordwise count.  Returns the new source (shallow copy)."""
    from copy import copy

    out = copy(src)
    g = np.unique(np.concatenate([np.arange(src.eta0, src.eta1, grid_step / src.length), [src.eta1]]))
    P = np.asarray(src.perimeter(g), float)
    out.lines = [*src.lines, Line("LE-", g, nd_le - half_width / P, src.eta0, src.eta1, True),
                 Line("LE+", g, nd_le + half_width / P, src.eta0, src.eta1, True)]
    return out


def _dedupe(points, eps=1e-9):
    out = []
    for p in points:
        if not out or math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > eps:
            out.append(p)
    if len(out) > 1 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) <= eps:
        out.pop()
    return out


# --------------------------------------------------------------------------------------------------- 3D + metrics
def map_to_3d(src, mesh):
    """Evaluate every chart node on the loft; seam duplicates share the master's position (nd = 0)."""
    eta = src.eta0 + mesh.chart[:, 1] / src.length
    nd = np.clip(mesh.chart[:, 0] / np.asarray(src.perimeter(eta), float), 0.0, 1.0)
    nd[mesh.seam_pairs[:, 0]] = 0.0
    xyz = np.asarray(src.loft(eta, nd), float)
    xyz[mesh.seam_pairs[:, 0]] = xyz[mesh.seam_pairs[:, 1]]
    return xyz, eta, nd


def unique_edges(elements):
    out = set()
    for e in np.asarray(elements):
        ids = e[e >= 0]
        for a, b in zip(ids, np.roll(ids, -1), strict=False):
            out.add((int(min(a, b)), int(max(a, b))))
    return np.array(sorted(out), dtype=int).reshape(-1, 2)


def element_metrics(xyz, elements, nref=None):
    """Per-element AR (longest/shortest edge), min/max corner angle [deg], min signed corner sine (Jacobian, > 0 is
    good, measured against ``nref`` (ne, 3) when given), 3D area. ``elements`` padded with -1 for triangles."""
    el = np.asarray(elements)
    ne = len(el)
    ar, amin, amax, jac, area = (np.empty(ne) for _ in range(5))
    tri = el[:, 3] < 0
    for sel, k in ((tri, 3), (~tri, 4)):
        if not sel.any():
            continue
        p = xyz[el[sel][:, :k]]
        e = np.roll(p, -1, axis=1) - p
        ln = np.linalg.norm(e, axis=-1)
        prev = -np.roll(e, 1, axis=1)
        den = np.maximum(np.roll(ln, 1, axis=1) * ln, 1e-300)
        cosv = np.einsum("nki,nki->nk", prev, e) / den
        ang = np.degrees(np.arccos(np.clip(cosv, -1, 1)))
        cr = np.cross(np.roll(e, 1, axis=1), e)
        if nref is None:
            nn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]) if k == 3 else np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1])
        else:
            nn = nref[sel]
        nn = nn / np.maximum(np.linalg.norm(nn, axis=1, keepdims=True), 1e-300)
        sine = np.einsum("nki,ni->nk", cr, nn) / den
        ar[sel], amin[sel], amax[sel], jac[sel] = ln.max(1) / np.maximum(ln.min(1), 1e-300), ang.min(1), ang.max(1), sine.min(1)
        cross = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]) if k == 3 else np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 1]) / 2
        area[sel] = np.linalg.norm(cross, axis=1) / (2 if k == 3 else 1)
    return {"aspect": ar, "min_angle": amin, "max_angle": amax, "jacobian": jac, "area": area, "triangle": tri}


DEFAULT_LIMITS = {"min_jacobian": 0.0, "max_ar": 10.0, "max_quad_angle": 165.0, "min_angle": 5.0}


def quality_gate(metrics, limits=None, raise_on_fail=True):
    """Count violations; REJECT (raise) when any exist. Nothing is repaired."""
    lim = {**DEFAULT_LIMITS, **(limits or {})}
    tri = metrics["triangle"]
    v = {"non_positive_jacobian": int((metrics["jacobian"] <= lim["min_jacobian"]).sum()),
         "ar_above_limit": int((metrics["aspect"] > lim["max_ar"]).sum()),
         "quad_max_angle_above_limit": int(((metrics["max_angle"] > lim["max_quad_angle"]) & ~tri).sum()),
         "min_angle_below_limit": int((metrics["min_angle"] < lim["min_angle"]).sum())}
    report = {"limits": lim, "violations": v, "passed": not any(v.values())}
    if raise_on_fail and not report["passed"]:
        raise MeshQualityError(f"quality gate rejected the mesh: {v}", report)
    return report


# --------------------------------------------------------------------------------------------------- driver
@dataclass
class SkinResult:
    xyz: np.ndarray
    elements: np.ndarray
    master: np.ndarray
    chart: np.ndarray
    eta: np.ndarray
    nd: np.ndarray
    metrics: dict
    gate: dict
    cleaned: Cleaned
    info: dict


def generate_skin_mesh(src, h, tau=0.25, mode="mid", size=None, limits=None, enforce_gate=True, simplify_tol=0.005,
                       snap=0.0, mesh_scale=1.0, h_range=(0.15, 0.25), structured=False, max_strip_ar=5.0, **gmsh_options):
    """clean -> linework -> gmsh -> 3D -> gate. ``h`` must lie in ``h_range`` ([0.15, 0.25] m, the project rule for the
    blade; a finer-mesh study of the prototype widens it explicitly)."""
    if not (np.isfinite(h) and h_range[0] <= h <= h_range[1]):
        raise ValueError(f"element size must be in [{h_range[0]}, {h_range[1]}] m")
    t = {}
    t0 = time.perf_counter()
    cleaned = clean_source(src, h, tau, mode, size, snap=snap, max_strip_ar=max_strip_ar)
    t["clean_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    size_m = size and (lambda e: mesh_scale * np.asarray(size(e), float)) if mesh_scale != 1.0 else size
    if structured:
        cm = mesh_chart_structured(src, cleaned, h, size_m, **gmsh_options)
        t["linework_s"] = 0.0
    else:
        edges = build_linework(src, cleaned, simplify_tol=simplify_tol)
        t["linework_s"] = time.perf_counter() - t0
        cm = mesh_chart(src, edges, h, size_m, **gmsh_options)
    t["gmsh_s"] = cm.info["gmsh_seconds"]
    t0 = time.perf_counter()
    xyz, eta, nd = map_to_3d(src, cm)
    t["map3d_s"] = time.perf_counter() - t0
    # reference normal: d(loft)/d(eta) x d(loft)/d(nd) at the element centroid (positive (eta, nd) orientation)
    ids = np.where(cm.elements >= 0, cm.elements, cm.elements[:, :1])
    ce = src.eta0 + cm.chart[ids, 1].mean(1) / src.length  # chart centroid (the seam duplicates keep nd = 1 here)
    cn = cm.chart[ids, 0].mean(1) / np.asarray(src.perimeter(ce), float)
    d1, d2 = 1e-4, 1e-4
    base = np.asarray(src.loft(ce, cn), float)
    nref = np.cross(np.asarray(src.loft(np.clip(ce + d1, 0, 1), cn), float) - base,
                    np.asarray(src.loft(ce, np.clip(cn + d2, 0, 1)), float) - base)
    metrics = element_metrics(xyz, cm.elements, nref)
    gate = quality_gate(metrics, limits, raise_on_fail=False)
    info = {**cm.info, **t, "n_elements": int(len(cm.elements)), "n_triangles": int(metrics["triangle"].sum()),
            "n_nodes": int(len(np.unique(cm.master))), "clean": cleaned.report}
    if enforce_gate and not gate["passed"]:
        raise MeshQualityError(f"quality gate rejected the mesh: {gate['violations']}", gate)
    return SkinResult(xyz, cm.elements, cm.master, cm.chart, eta, nd, metrics, gate, cleaned, info)


# --------------------------------------------------------------------------------------------------- real blade
class RowLoft:
    """Loft evaluation from a few hundred ``Geometry.generate`` stations (cost ~0.08 s each, independent of batch
    size), cubic-interpolated across eta at fixed nd; exact Pchip in nd on each station. Float64 path."""

    def __init__(self, topology, spacing=1.0, tip_spacing=0.5, tip_length=12.0, eta_bounds=(0.0, 1.0), end_spacing=0.0,
                 end_length=1.0):
        from copy import copy

        from aeroelast.models.blade.numad.objects.geometry import Geometry

        d = copy(topology._definition)
        self.length = L = float(d.span[-1])
        if end_spacing:  # the loft closes in the last metre: a finer station set avoids interpolation overshoot there
            s = np.concatenate([np.arange(0, L - tip_length, spacing), np.arange(L - tip_length, L - end_length, tip_spacing),
                                np.arange(L - end_length, L, end_spacing), [L]])
        else:
            s = np.concatenate([np.arange(0, L - tip_length, spacing), np.arange(L - tip_length, L, tip_spacing), [L]])
        e = s / L
        first = max(int(np.searchsorted(e, eta_bounds[0], side="right")) - 2, 0)
        last = min(int(np.searchsorted(e, eta_bounds[1], side="left")) + 2, len(e))
        s = s[first:last]
        self.rows = s / L
        d.ispan = s
        g = Geometry().generate(d, windio2_eta=self.rows)
        self.stations = []
        per = []
        for i in range(len(self.rows)):
            arc = np.asarray(g.arclength[:, i], float)
            dist, idx = np.unique(arc, return_index=True)
            f = PchipInterpolator(dist, np.asarray(g.coordinates[idx, :, i], float))
            self.stations.append((f, arc[0], arc[-1] - arc[0]))
            per.append(arc[-1] - arc[0])
        self.P = PchipInterpolator(self.rows, per)

    def perimeter(self, eta):
        return self.P(np.asarray(eta, float))

    def station(self, i, nd):
        f, a0, per = self.stations[i]
        return f(a0 + (1.0 - nd) * per)

    def __call__(self, eta, nd):
        eta, nd = np.broadcast_arrays(np.asarray(eta, float), np.asarray(nd, float))
        j = np.clip(np.searchsorted(self.rows, eta, side="right") - 1, 1, len(self.rows) - 3)
        X = np.empty((len(eta), 4, 3))
        for r in range(4):
            idx = j - 1 + r
            for i in np.unique(idx):
                m = idx == i
                X[m, r] = self.station(int(i), nd[m])
        t = self.rows[j[:, None] - 1 + np.arange(4)[None, :]]
        out = np.zeros((len(eta), 3))
        for r in range(4):
            w = np.ones(len(eta))
            for q in range(4):
                if q != r:
                    w *= (eta - t[:, q]) / (t[:, r] - t[:, q])
            out += w[:, None] * X[:, r]
        return out


def le_landmark_line(loft):
    """Fixed line nd_LE(eta) of the loft's leading-edge landmark (arc 0: the highest-curvature point of every section).
    The v2 mesher makes it a mandatory ring cut; without it one element wraps the apex and the edge chord misses it."""
    nd = np.array([1.0 + a0 / per for _, a0, per in loft.stations])
    return Line("le/landmark", np.asarray(loft.rows, float), nd, float(loft.rows[0]), float(loft.rows[-1]), True)


def source_from_topology(topology, loft, include_attachments=True, eta_bounds=(0.0, 1.0), le_landmark=False):
    """:class:`SkinSource` of the windIO v2 skin: layer boundary curves (+ web attachment lines, fixed), de-duplicated
    when identical. Lines that sit on the seam (nd == 0 or 1 everywhere) are dropped (they are the boundary)."""
    from aeroelast.models.blade.numad.objects.windio2_laminate import span_domain

    curves = {c.identity: c for c in topology.curves if c.surface == "skin"}
    lam = topology._laminate
    dom, drops = {}, []
    for layer, _ in lam.layers:
        if layer.get("web", "skin") != "skin":
            continue
        ids = [("layer", layer["name"], h) for h in ("start_nd_arc", "end_nd_arc")]
        if not all(i in curves for i in ids):
            continue
        lo, hi = span_domain(layer)
        lo = max(lo, *(curves[i].grid[0] for i in ids))
        hi = min(hi, *(curves[i].grid[-1] for i in ids))
        for i in ids:
            a, b = dom.get(i, (lo, hi))
            dom[i] = (min(a, lo), max(b, hi))
        if lo > 0:
            drops.append((lo, ids, "right"))
        if hi < 1:
            drops.append((hi, ids, "left"))
    fixed_ids = set()
    if include_attachments:
        for web in topology.source["components"]["blade"]["structure"].get("webs", ()):
            lo, hi = span_domain(web)
            for h in ("start_nd_arc", "end_nd_arc"):
                i = ("web", web["name"], h)
                if i in curves:
                    dom[i] = (max(lo, curves[i].grid[0]), min(hi, curves[i].grid[-1]))
                    fixed_ids.add(i)
    lines, alias = [], {}
    for ident, (lo, hi) in dom.items():
        c = curves[ident]
        v = np.asarray(c.values, float)
        if np.all(v < 1e-12) or np.all(v > 1 - 1e-12):
            alias[ident] = None
            continue
        dup = next((ln for ln in lines if ln.fixed == (ident in fixed_ids) and (lo, hi) == (ln.lo, ln.hi)
                    and np.array_equal(ln.grid, np.asarray(c.grid)) and np.allclose(ln.values, v, atol=1e-12, rtol=0)), None)
        if dup:
            alias[ident] = dup.name
            continue
        name = "/".join(ident)
        lines.append(Line(name, np.asarray(c.grid, float), v, float(lo), float(hi), ident in fixed_ids))
        alias[ident] = name
    if le_landmark:
        lines.append(le_landmark_line(loft))
    out_drops = []
    for eta, ids, side in drops:
        a, b = (alias.get(i) for i in ids)
        if a and b and a != b and eta_bounds[0] <= eta <= eta_bounds[1]:
            out_drops.append((eta, a, b, side))
    extra = [float(e) for e in topology.event_stations if eta_bounds[0] <= e <= eta_bounds[1]]
    return SkinSource(lines, loft.length, loft.perimeter, out_drops, loft, float(eta_bounds[0]), float(eta_bounds[1]), extra)


def lookup_cells(topology, surface, eta, q, length=None, dy=0.0):
    """Source :class:`LaminateCell` of ``surface`` ('skin' or a web name) at every chart point ``(eta, q)``.

    The cells are the ones ``topology.section(eta)`` would return (same ``_positions`` / ``_laminate.cells`` calls, without
    its per-eta loft rebuild).  ``dy`` > 0 (chart metres, needs ``length``) evaluates the laminate once per span bucket
    of that height at the bucket's mean eta, for all the cells of the bucket in ONE call; ``dy = 0`` is exact.
    ``LaminateCell.xyz`` is a placeholder (zeros): no consumer of the raw mesh reads it."""
    eta, q = np.asarray(eta, float), np.asarray(q, float)
    out = [None] * len(eta)

    def positions(e):
        try:
            return topology._positions(float(e), surface)
        except ValueError:
            return topology._positions(float(e) + 1e-12, surface)

    if dy > 0:
        if not length:
            raise ValueError("dy > 0 needs the span length")
        bucket = np.floor(eta * length / dy).astype(int)
        groups = [np.nonzero(bucket == b)[0] for b in np.unique(bucket)]
    else:
        groups = [np.array([i]) for i in range(len(eta))]
    zero = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    for idx in groups:
        e = float(np.mean(eta[idx])) if dy > 0 else float(eta[idx[0]])
        pos = positions(e)
        vals = np.array(sorted({0.0, 1.0, *pos.values()}))
        lo_i = np.searchsorted(vals, q[idx], side="right") - 1
        lo = vals[lo_i]
        nxt = np.minimum(lo_i + 1, len(vals) - 1)
        hi = np.where(vals[nxt] > lo, vals[nxt], 1.0)
        pairs = sorted(set(zip(lo.tolist(), hi.tolist(), strict=True)))
        cells = topology._laminate.cells(e, surface, pos, tuple(pairs), tuple(zero for _ in pairs))
        look = dict(zip(pairs, cells, strict=True))
        for i, l_, h_ in zip(idx, lo.tolist(), hi.tolist(), strict=True):
            out[i] = look[(l_, h_)]
    return out


def areal_mass(topology, eta, nd, length=None, dy=0.0):
    """Per-point (kg/m2, stack key) from the source laminate cell containing (eta, nd); same lookup as the v2 mesher
    (``lookup_cells``).  ``dy`` > 0 buckets the span (see there): the thickness fields are piecewise linear in eta, so the
    difference to the exact per-point lookup is below ``slope * dy``."""
    cells = lookup_cells(topology, "skin", eta, nd, length, dy)
    return (np.array([sum(float(np.mean(np.asarray(s.material["rho"], float))) * s.thickness_m for s in c.layers) for c in cells]),
            [tuple(s.name for s in c.layers) for c in cells])


# --------------------------------------------------------------------------------------------------- skin + webs
def _pad_corners(values, elements):
    """(ne, 4, k) corner values; a triangle repeats its first corner (the convention of the v2 mesher)."""
    idx = np.where(elements >= 0, elements, elements[:, :1])
    return values[idx]


def conformity_report(nodes, elements, tol=1e-9, surfaces=None):
    """Distinct node ids at coincident positions (must be 0 between skin and webs), unused nodes, duplicate-corner
    elements.  With ``surfaces`` (one name per element) the coincident pairs are split by the surfaces that own them."""
    nodes = np.asarray(nodes, float)
    pairs = cKDTree(nodes).query_pairs(tol, output_type="ndarray")
    used = np.zeros(len(nodes), bool)
    used[elements[elements >= 0]] = True
    dup = sum(len(set(r[r >= 0].tolist())) != int((r >= 0).sum()) for r in elements)
    rep = {"coincident_distinct_id_pairs": int(len(pairs)), "unused_nodes": int((~used).sum()), "elements_with_repeated_node": int(dup)}
    if surfaces is not None:
        own = {}
        for i, r in enumerate(elements):
            for v in r[r >= 0]:
                own.setdefault(int(v), set()).add(surfaces[i])
        kinds = {}
        for a, b in pairs:
            key = "|".join(sorted(own.get(int(a), set()))) + " ~ " + "|".join(sorted(own.get(int(b), set())))
            kinds[key] = kinds.get(key, 0) + 1
        rep["coincident_pairs_by_surfaces"] = kinds
        rep["skin_web_coincident_pairs"] = int(sum(v for k, v in kinds.items() if ("skin" in k.split(" ~ ")[0]) != ("skin" in k.split(" ~ ")[1])
                                                    or ("skin" in k.split(" ~ ")[0] and "skin" in k.split(" ~ ")[1] and False)))
    return rep


def mesh_web_chart(xyz, start_nodes, end_nodes, rows, active, h, size, length, eta0=0.0, eta1=1.0, *, n_new_first=None, **mesh_opts):
    """Mesh ONE shear web as a chart (span y, height x = q * width) between two skin attachment node chains.

    ``start_nodes`` / ``end_nodes``: node id of the skin on the web's two attachment lines at every row (-1 where the
    skin has none); ``xyz`` the coordinates of those ids.  ``rows`` (span positions [m]) are the rows of the skin chart, so
    every web edge node IS a skin node (shared id, no hanging node at the junction); ``active`` (bool per row interval)
    keeps only the intervals where the web has material.  The interior is meshed by :func:`mesh_chart_structured`
    (quads, transition patches of triangles/quads where the element count across the height changes), element size ``size(eta)``.
    Returns ``dict(xyz_new, elements, params, info)``; new node ids start at ``n_new_first`` (default ``len(xyz)``), the elements
    reference skin ids on the two edges.  Raises :class:`MeshQualityError` if an active interval lacks a skin node."""
    rows = np.asarray(rows, float)
    nr = len(rows)
    n0 = len(xyz) if n_new_first is None else n_new_first
    ns, ne_ = np.asarray(start_nodes), np.asarray(end_nodes)
    sx, ex = xyz[np.where(ns >= 0, ns, 0)], xyz[np.where(ne_ >= 0, ne_, 0)]
    width = np.where((ns >= 0) & (ne_ >= 0), np.linalg.norm(ex - sx, axis=1), np.nan)
    ok = np.isfinite(width)
    wfill = np.interp(rows, rows[ok], width[ok])
    eta_rows = eta0 + rows / length
    web_src = SkinSource([], length, lambda e, w=wfill: np.interp(np.asarray(e, float), eta_rows, w), [], None, eta0, eta1)
    cm = mesh_chart_structured(web_src, Cleaned(eta_rows, {}, [], {}), h, size, rows=rows, active=active, periodic=False,
                               **{k: v for k, v in mesh_opts.items() if k != "strip_scale"})
    ch = cm.chart
    k = np.clip(np.searchsorted(rows, ch[:, 1] + 1e-9) - 1, 0, nr - 2)
    t = np.clip((ch[:, 1] - rows[k]) / (rows[k + 1] - rows[k]), 0.0, 1.0)
    S = sx[k] + (sx[k + 1] - sx[k]) * t[:, None]
    E = ex[k] + (ex[k + 1] - ex[k]) * t[:, None]
    w = wfill[k] + (wfill[k + 1] - wfill[k]) * t
    q = np.clip(ch[:, 0] / w, 0.0, 1.0)
    wxyz = S + q[:, None] * (E - S)
    r = np.clip(np.searchsorted(rows, ch[:, 1] - 1e-9), 0, nr - 1)  # nearest row at or above y (an edge node sits ON a row)
    on_row = np.abs(ch[:, 1] - rows[r]) < 1e-9
    start_n = np.abs(ch[:, 0]) < 1e-9
    end_n = np.abs(ch[:, 0] - w) < 1e-9 * max(1.0, w.max())
    if np.any((start_n | end_n) & ~on_row):
        raise MeshQualityError("web: an edge node is not on a row")
    new_id = np.full(len(ch), -1, dtype=int)
    side_id = np.where(start_n, ns[r], ne_[r])
    if np.any((start_n | end_n) & (side_id < 0)):
        raise MeshQualityError("web: an active row interval has no skin node on its attachment line")
    new_id[start_n | end_n] = side_id[start_n | end_n]
    interior = new_id < 0
    new_id[interior] = n0 + np.arange(int(interior.sum()))
    wel = np.where(cm.elements >= 0, new_id[np.where(cm.elements >= 0, cm.elements, 0)], -1)
    par = _pad_corners(np.column_stack([eta0 + ch[:, 1] / length, q]), cm.elements)
    return {"xyz_new": wxyz[interior], "elements": wel, "params": par, "info": cm.info}


def _chart_jacobian(params, elements):
    """Signed chart Jacobian per element (centred differences, triangles and quads); ``params`` is (ne, 4, 2) with NaN padding."""
    uv = np.nan_to_num(params)
    tri = elements[:, 3] < 0
    du = np.where(tri[:, None], uv[:, 1] - uv[:, 0], uv[:, 1] - uv[:, 0] + uv[:, 2] - uv[:, 3])
    dv = np.where(tri[:, None], uv[:, 2] - uv[:, 0], uv[:, 3] - uv[:, 0] + uv[:, 2] - uv[:, 1])
    return du[:, 0] * dv[:, 1] - du[:, 1] * dv[:, 0]


def _orient_positive(elements, params):
    """Reverse every element whose chart orientation is negative, so the connectivity normal follows the chart (what
    ``windio2_mesh.map_mesh_materials`` assumes).  A zero Jacobian is not repaired: it raises."""
    elements, params = elements.copy(), params.copy()
    det = _chart_jacobian(params, elements)
    if np.any(det == 0):
        raise MeshQualityError(f"{int((det == 0).sum())} elements have a zero chart Jacobian")
    tri = elements[:, 3] < 0
    params[tri, 3] = params[tri, 0]  # triangles repeat one corner (windio2_mesh._centroid sums params * live-mask: NaN would poison it)
    for i in np.nonzero(det < 0)[0]:
        order = [0, 2, 1] if tri[i] else [0, 3, 2, 1]
        elements[i, : len(order)] = elements[i, order]
        params[i, : len(order)] = params[i, order]
    return elements, params


def map_mesh_materials_gmsh(topology, mesh):
    """``windio2_mesh.map_mesh_materials`` for the gmsh mesh: same frames, normals and web-to-LE orientation, with ONE change.

    An element whose frame is degenerate with the beam axis (``map_stack`` raises "Degenerate beam-axis") takes the element's first
    shell-basis axis as the reference direction instead.  Only the MAPPED fibre angle depends on it (the SOURCE angle, which the FEM
    reads, does not); the number of such elements is recorded in ``mesh``-independent form by the caller through ``degenerate_axis``.
    """
    from dataclasses import replace

    from ..objects.geometry import Geometry
    from .windio2_basis import map_stack, shell_bases
    from .windio2_mesh import Pchip, _centroid

    eta = _centroid(mesh.element_parameters, mesh.elements)[:, 0]
    definition = topology._definition
    frames = Geometry()
    frames.ispan = eta * definition.span[-1]
    twist = definition.windio2_outer_shape["twist"]
    frames._windio2_eta = eta
    frames._v2_pchip = Pchip
    frames.idegreestwist = Pchip(twist["grid"], twist["values"])(eta)
    frames._build_windio2_frames(definition)
    canonical = shell_bases(mesh.nodes, mesh.elements)[:, 2]
    det = _chart_jacobian(np.asarray(mesh.element_parameters, float), mesh.elements)
    if np.any(det == 0):
        raise ValueError("Degenerate material mapping chart")
    normals = canonical * np.sign(det)[:, None]
    is_web = np.array([c.surface != "skin" for c in mesh.cells])
    # A web's physical normal points to the LE side.  ``map_stack`` needs the connectivity normal to BE that normal, so a web whose
    # connectivity normal points the other way is reversed (as the windio2 mesher builds its webs); the chart normal is unchanged.
    to_le = -np.array([f[:, 1] for f in frames._windio2_frames])
    alignment = np.einsum("ij,ij->i", normals, to_le)
    if np.any(is_web & (np.abs(alignment) <= 1e-12)):
        raise ValueError("Web normal has no unambiguous LE side")
    normals = np.where(is_web[:, None], normals * np.sign(alignment)[:, None], normals)
    reverse = is_web & (np.einsum("ij,ij->i", canonical, normals) < 0)
    if reverse.any():
        elements, params = mesh.elements.copy(), np.array(mesh.element_parameters, float)
        for i in np.nonzero(reverse)[0]:
            order = [0, 2, 1] if elements[i, 3] < 0 else [0, 3, 2, 1]
            elements[i, : len(order)] = elements[i, order]
            params[i, : len(order)] = params[i, order]
        mesh = replace(mesh, elements=elements, element_parameters=params, bases=shell_bases(mesh.nodes, elements))
    stacks = []
    for cell, frame, normal, basis in zip(mesh.cells, frames._windio2_frames, normals, mesh.bases, strict=True):
        try:
            stacks.append(map_stack(cell, frame[:, 2], normal, basis))
        except ValueError as exc:
            if "Degenerate beam-axis" not in str(exc):
                raise
            stacks.append(map_stack(cell, basis[0], normal, basis))
    return replace(mesh, material_stacks=tuple(stacks))


def generate_mesh_gmsh(topology, element_size=0.25, *, tau=0.25, snap=0.4, ring=16, structured_options=None, le_landmark=True,
                       with_materials=False, material_dy=0.0, enforce_gate=True, h_range=(0.15, 0.25), loft=None,
                       return_info=False, max_strip_ar=5.0):
    """Skin + webs raw mesh with the STRUCTURED gmsh mesher; a drop-in for ``windio2_mesh.generate_mesh`` outputs.

    skin   ``generate_skin_mesh(structured=True)`` on the official topology (web attachment lines are fixed lines of the
           skin, so the skin carries one node per row on each of them);
    webs   one structured chart (eta, q) per web with the SAME rows as the skin, only on the row intervals where the web
           layup has material (``thickness > 0``), q = distance between its two attachment nodes at the row;
           the q = 0 / 1 edges are the skin nodes (1 division per row segment, shared ids);
    out    :class:`windio2_mesh.RawMesh` (nodes, elements with -1 padded triangles, cells, surfaces ['skin', web...],
           element_parameters (eta, nd) / (eta, q), shell bases, attachments chains) and, with ``with_materials``,
           ``map_mesh_materials``.  Nodes are float64 (the v2 mesher uses extended precision)."""
    from .windio2_basis import shell_bases
    from .windio2_mesh import RawMesh

    t0 = time.perf_counter()
    h = float(element_size)
    loft = loft or RowLoft(topology, end_spacing=0.1)
    src = source_from_topology(topology, loft, le_landmark=le_landmark)
    L = loft.length
    size = ring_size(h, loft.perimeter, ring)
    opts = {"patch_rows": 1, "extra_options": {"Mesh.Smoothing": 10.0, "Mesh.RecombineMinimumQuality": 0.7},
                "optimize": ("Relocate2D",), "optimize_iter": 5, "optimize_scope": "patches"}
    opts.update(structured_options or {})
    res = generate_skin_mesh(src, h, tau, "mid", size, enforce_gate=enforce_gate, snap=snap, structured=True, h_range=h_range,
                             max_strip_ar=max_strip_ar, **opts)
    inf = res.info
    rows, _names, line_nodes = inf["rows"], inf["line_names"], inf["line_row_nodes"]
    nr = len(rows)
    eta_rows = src.eta0 + rows / L
    t_skin = time.perf_counter() - t0

    # ---- skin
    el = res.elements
    used = np.unique(res.master[el[el >= 0]])
    remap = np.full(len(res.master), -1, dtype=int)
    remap[used] = np.arange(len(used))

    def merged(i):
        return remap[res.master[i]]

    skin_el = np.where(el >= 0, merged(np.where(el >= 0, el, 0)), -1)
    nodes = [res.xyz[used]]
    nd1 = np.clip(res.chart[:, 0] / np.asarray(src.perimeter(res.eta), float), 0.0, 1.0)  # the right seam keeps nd = 1
    skin_par = _pad_corners(np.column_stack([res.eta, nd1]), el)
    live = (el >= 0)[:, :, None]
    c_eta, c_nd = ((skin_par * live).sum(1) / live.sum(1)).T
    t0 = time.perf_counter()
    skin_cells = lookup_cells(topology, "skin", c_eta, c_nd, L, material_dy)
    t_cells = time.perf_counter() - t0
    elements, params, cells, surfaces = [skin_el], [skin_par], list(skin_cells), ["skin"] * len(skin_el)

    # ---- webs
    attachments, web_info = {}, {}
    n_nodes = len(used)
    web_list = list(topology.source["components"]["blade"]["structure"].get("webs", ()))
    from aeroelast.models.blade.numad.objects.windio2_laminate import span_domain

    for web in web_list:
        nm = web["name"]
        keys = [f"web/{nm}/{hd}" for hd in ("start_nd_arc", "end_nd_arc")]
        if not all(k in line_nodes for k in keys):
            continue
        ns, ne_ = line_nodes[keys[0]], line_nodes[keys[1]]
        lo_, hi_ = span_domain(web)

        def has_material(eta, nm=nm):
            try:
                pos = topology._positions(float(eta), nm)
            except ValueError:
                pos = topology._positions(float(eta) + 1e-12, nm)
            vals = sorted({0.0, 1.0, *pos.values()})
            bounds = tuple(zip(vals[:-1], vals[1:], strict=True))
            zero = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
            return any(c.thickness_m > 0 for c in topology._laminate.cells(float(eta), nm, pos, bounds, tuple(zero for _ in bounds)))

        mid = 0.5 * (eta_rows[:-1] + eta_rows[1:])
        material = np.array([lo_ <= mid[k] <= hi_ and has_material(mid[k]) for k in range(nr - 1)])
        nodes_ok = np.array([ns[k] >= 0 and ns[k + 1] >= 0 and ne_[k] >= 0 and ne_[k + 1] >= 0 for k in range(nr - 1)])
        if np.any(material & ~nodes_ok):
            raise MeshQualityError(f"web {nm}: {int((material & ~nodes_ok).sum())} row intervals carry web material but the skin "
                                   f"has no node on its attachment line at one of the rows")
        active = material
        web_info[nm] = {"active_row_intervals": int(active.sum()), "of": int(nr - 1)}
        if not active.any():
            continue
        ns_f = np.where(ns >= 0, merged(np.where(ns >= 0, ns, 0)), -1)
        ne_f = np.where(ne_ >= 0, merged(np.where(ne_ >= 0, ne_, 0)), -1)
        wm = mesh_web_chart(nodes[0], ns_f, ne_f, rows, active, h, size, L, src.eta0, src.eta1, n_new_first=n_nodes, **opts)
        nodes.append(wm["xyz_new"])
        n_nodes += len(wm["xyz_new"])
        wel, wpar = wm["elements"], wm["params"]
        wlive = (wel >= 0)[:, :, None]
        we, wq = ((wpar * wlive).sum(1) / wlive.sum(1)).T
        t0 = time.perf_counter()
        wcells = lookup_cells(topology, nm, we, wq, L, material_dy)
        t_cells += time.perf_counter() - t0
        elements.append(wel)
        params.append(wpar)
        cells.extend(wcells)
        surfaces.extend([nm] * len(wel))
        act_rows = sorted(set(np.nonzero(active)[0]) | {kk + 1 for kk in np.nonzero(active)[0]})
        for hd, arr in zip(("start_nd_arc", "end_nd_arc"), (ns_f, ne_f), strict=True):
            attachments[(nm, hd)] = tuple(int(arr[r]) for r in act_rows)
        web_info[nm].update(elements=int(len(wel)), triangles=int((wel[:, 3] < 0).sum()), patches=wm["info"]["n_patches"],
                            count_changes=wm["info"]["n_count_changes"], nodes_new=int(len(wm["xyz_new"])), gmsh_seconds=wm["info"]["gmsh_seconds"])

    nodes = np.vstack(nodes)
    elements, params_all = _orient_positive(np.vstack(elements), np.vstack(params))
    mesh = RawMesh(nodes, elements, tuple(cells), tuple(surfaces), params_all, shell_bases(nodes, elements), attachments)
    t0 = time.perf_counter()
    if with_materials:
        mesh = map_mesh_materials_gmsh(topology, mesh)
    info = {"skin_result": res, "webs": web_info, "n_skin_elements": len(skin_el), "gate": res.gate, "seconds_skin_mesh": t_skin,
            "seconds_cells": t_cells, "seconds_materials_map": time.perf_counter() - t0, "rows": rows,
            "conformity": conformity_report(nodes, elements, surfaces=surfaces)}
    return (mesh, info) if return_info else mesh
