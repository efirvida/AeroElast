"""Thin-wall multi-cell Bredt GJ of a windIO v1 yaml section, no mesh (odd/tasks/l3-closure.md, T9).

Independent of the shell FE: the section comes only from the yaml-side reader of
``tools/audit_blade_input.py`` (``YamlBlade``), and is compared with the BeamDyn deck K66 and its
shear-condensed value. Conventions:

* wall line: the open airfoil polyline x chord, nd_arc 0 = suction TE, 1 = pressure TE, resampled by
  arc length with the web attachment arcs as exact samples; the TE is closed by a straight segment
  whose laminate is the one active at nd_arc 0 (TE reinforcement included);
* Gt of a segment = sum over active layers of G12 * t at the segment's arc midpoint; the mid-line is
  the polyline offset inward by half the local laminate thickness (mitre, clamped), or the outer
  mould line with ``--no-offset``;
* webs: straight chords between the mid-line points at the web arcs, Gt of ``web_layers``; a web with
  no layer at the station does not exist;
* cells are split by the webs ordered by their suction-side arc (cell 0 = TE cell, last = LE cell);
* BeamDyn 6x6 order is [Fx Fy Fz Mx My Mz]: 0,1 shear, 2 axial, 3,4 bending, 5 torsion.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_blade_input import YamlBlade  # noqa: E402

MITRE_MAX = 3.0


# =============================================================================================
# Bredt core (geometry arrays in, numbers out)
# =============================================================================================
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
def _wall(yb: YamlBlade, frac: float, a: float, perimeter: float, cache: dict):
    layers = yb.active_layers(frac, a, perimeter)
    gt = t = 0.0
    for lay in layers:
        key = lay["material"]
        if key not in cache:
            cache[key] = yb.material(key)["g12"]
        gt += cache[key] * lay["thickness"]
        t += lay["thickness"]
    return gt, t


def _web(yb: YamlBlade, k: int, frac: float, cache: dict):
    gt = t = 0.0
    for lay in yb.web_layers(k, frac):
        if lay["material"] not in cache:
            cache[lay["material"]] = yb.material(lay["material"])["g12"]
        gt += cache[lay["material"]] * lay["thickness"]
        t += lay["thickness"]
    return gt, t


def _offset_loop(pts: np.ndarray, h: np.ndarray) -> np.ndarray:
    """Move each vertex inward by h (per vertex) along the mitre direction of the closed loop."""
    nxt = np.roll(pts, -1, axis=0)
    seg = nxt - pts
    ln = np.linalg.norm(seg, axis=1)
    tan = seg / ln[:, None]
    sigma = np.sign(_cross_sum(pts, nxt))
    nrm = sigma * np.c_[-tan[:, 1], tan[:, 0]]  # inward for either orientation
    n_in, n_out = np.roll(nrm, 1, axis=0), nrm
    bis = n_in + n_out
    nb = np.linalg.norm(bis, axis=1)
    t_in, t_out = np.roll(tan, 1, axis=0), tan
    alt = sigma * np.c_[-(t_out - t_in)[:, 1], (t_out - t_in)[:, 0]]
    alt = np.where(np.linalg.norm(alt, axis=1)[:, None] > 0, alt / np.maximum(np.linalg.norm(alt, axis=1), 1e-300)[:, None], nrm)
    direction = np.where((nb > 1e-6)[:, None], bis / np.maximum(nb, 1e-300)[:, None], alt)
    cos_half = np.clip(np.sum(direction * n_out, axis=1), 1.0 / MITRE_MAX, 1.0)
    return pts + (h / cos_half)[:, None] * direction


def yaml_section(yb: YamlBlade, frac: float, n: int = 800, offset: bool = True) -> dict:
    chord = yb.outer_shape(frac)["chord"]
    poly, _, scheme = yb.airfoil_at(frac)
    poly = chord * np.asarray(poly, float)
    perimeter = yb.perimeter(frac)
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(poly, axis=0), axis=1))]
    s /= s[-1]

    webs_raw = []
    for k in range(len(yb.struct["webs"])):
        gt, t = _web(yb, k, frac, {})
        if gt > 0.0:
            ws = yb.web_state(k, frac)
            a, b = sorted((ws["start"], ws["end"]))
            webs_raw.append({"k": k, "name": ws["name"], "a_s": a, "a_p": b, "gt": gt, "t": t})
    webs_raw.sort(key=lambda w: w["a_s"])

    u = np.linspace(0.0, 1.0, n)
    extra = np.array([w[key] for w in webs_raw for key in ("a_s", "a_p")], float)
    u = np.unique(np.r_[u, extra])
    pts = np.c_[np.interp(u, s, poly[:, 0]), np.interp(u, s, poly[:, 1])]
    if np.linalg.norm(pts[-1] - pts[0]) < 1e-9 * chord:  # sharp TE: the loop closes on its own vertex
        pts, u = pts[:-1], u[:-1]
    m = len(pts)
    cache: dict = {}
    mid = np.r_[0.5 * (u[:-1] + u[1:]), 0.0 if m == len(u) else 0.5 * (u[-1] + 1.0)]
    gt_t = [_wall(yb, frac, a, perimeter, cache) for a in mid]
    seg_gt = np.array([g for g, _ in gt_t])
    seg_t = np.array([t for _, t in gt_t])
    if np.any(seg_gt <= 0.0):
        raise ValueError(f"no active laminate at nd_arc {mid[np.argmin(seg_gt)]:.4f}, frac {frac}")
    closing = m == len(u)  # the closing segment exists only when the TE is open
    idx = {key: [int(np.argmin(abs(u - w[key]))) for w in webs_raw] for key in ("a_s", "a_p")}
    webs = [{"i_s": i, "i_p": j, "gt": w["gt"], "t": w["t"]}
            for w, i, j in zip(webs_raw, idx["a_s"], idx["a_p"], strict=True)]
    a_s = np.array([w["a_s"] for w in webs_raw])
    a_p = np.array([w["a_p"] for w in webs_raw])
    seg_cell = np.array([int(np.sum(a_s < x) - np.sum(a_p < x)) for x in mid])
    mid_pts = _offset_loop(pts, 0.5 * 0.5 * (seg_t + np.roll(seg_t, 1))) if offset else pts
    out = multicell_bredt(mid_pts, seg_gt, seg_cell, webs, seg_t=seg_t)
    out.update(_regions(out, seg_cell, webs_raw))
    out.update({"frac": frac, "chord": chord, "scheme": scheme, "offset": offset, "n_pts": m,
                "te_closure_gt_at_nd_arc0": float(seg_gt[-1]) if closing else None,
                "web_names": [w["name"] for w in webs_raw]})
    for key in ("F", "q"):
        out[key] = np.asarray(out[key]).tolist()
    out["areas"] = np.asarray(out["areas"]).tolist()
    out.pop("flex_ext"), out.pop("len_ext")
    return out


def _regions(out: dict, seg_cell, webs_raw) -> dict:
    fe, le = out["flex_ext"], out["len_ext"]
    nc = out["n_cells"]
    reg = {"te_cell_exterior": (fe[0], le[0])}
    if nc > 2:
        reg["box_exterior"] = (fe[1:-1].sum(), le[1:-1].sum())
    if nc > 1:
        reg["le_cell_exterior"] = (fe[-1], le[-1])
    for w, wo in zip(webs_raw, out["webs"], strict=True):
        reg[f"web:{w['name']}"] = (wo["flex"], wo["length"])
    total = sum(v[0] for v in reg.values())
    regions = {k: {"flex": v[0], "length": v[1], "share": v[0] / total} for k, v in reg.items()}
    return {"regions": regions, "flex_total": total, "te_share": regions["te_cell_exterior"]["share"]}


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


# =============================================================================================
# CLI
# =============================================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--yaml", required=True)
    ap.add_argument("--beamdyn", required=True)
    ap.add_argument("--stations", type=float, nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--no-offset", action="store_true", help="report only the outer mould line")
    args = ap.parse_args(argv)

    yb = YamlBlade.from_file(args.yaml)
    etas, ks = load_beamdyn_6x6(args.beamdyn)
    print(f"blade length {yb.length:.3f} m; deck stations {len(etas)}; K order [Fx Fy Fz Mx My Mz], "
          f"torsion = K[5,5] condensed over shear 0,1; TE closure laminate = active at nd_arc 0")
    print(f"{'eta':>5} {'Bredt_off':>11} {'Bredt_omt':>11} {'open':>10} {'K66':>11} {'GJ_sc':>11} "
          f"{'%off':>7} {'%omt':>7} {'TE share':>8}")
    rows = []
    for eta in args.stations:
        k = interp_deck(etas, ks, eta)
        gjsc = condense_shear(k)
        a = yaml_section(yb, eta, args.n, True)
        b = yaml_section(yb, eta, args.n, False)
        pa, pb = 100 * (a["gj"] / gjsc - 1), 100 * (b["gj"] / gjsc - 1)
        print(f"{eta:5.2f} {a['gj']:11.4e} {b['gj']:11.4e} {a['open']:10.3e} {k[5, 5]:11.4e} {gjsc:11.4e} "
              f"{pa:+7.1f} {pb:+7.1f} {a['te_share']:8.3f}")
        rows.append({"eta": eta, "k66": k[5, 5], "gj_sc": gjsc, "pct_offset": pa, "pct_no_offset": pb,
                     "offset": a, "no_offset": b})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps({"yaml": args.yaml, "beamdyn": args.beamdyn, "n": args.n,
                                          "stations": rows}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
