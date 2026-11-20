"""Mesh-quality measurement and VTU export for the IEA-15 blade shell meshes.

Measurement only: nothing here changes a mesh, a tolerance or the physics.

Usage (always through ``scripts/aeroenv.sh python``)::

    tools/mesh_quality.py v2 <yaml> <out_dir> [--element-size 0.25] [--tag NAME]   # gmsh on a windIO 2.x yaml; default sizing N20 AR 1
    tools/mesh_quality.py v2 --from-npz <raw.npz> <out_dir> [--tag NAME]
    tools/mesh_quality.py summary <out_dir> <tag> [<tag> ...]

Writes ``<out_dir>/<tag>.vtu`` (ParaView) and ``quality_<tag>.json``.

Metrics, per element in float64 (quad corners p0..p3 counter-clockwise, edge
``ei = p[i+1]-p[i]``; triangles use the three edges and corners):

* edge lengths ``|ei|``
* aspect ratio        = longest edge / shortest edge
* stretch (quad)      = max(r, 1/r), r = mean(|e0|,|e2|) / mean(|e1|,|e3|)
* min / max angle     = extreme interior corner angles (3D vectors), degrees
* warping (quad)      = for each diagonal split, angle between the unit normals of the
  two triangles; the larger of the two diagonals, degrees (triangles: 0)
* scaled Jacobian     = min over corners of ((a x b).n)/(|a||b|), a/b the next/previous
  edge vectors from the corner, n the unit element normal (quad: cross of diagonals;
  triangle: cross of two edges)
* taper (quad, Verdict) = max(|X12|/min(|X1|,|X2|)...) with X1 = p1+p2-p0-p3,
  X2 = p2+p3-p0-p1, X12 = p0-p1+p2-p3, taper = |X12| / min(|X1|,|X2|)

Span position uses the axis with the largest node extent (the blade axis), 10 equal bins.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

METRICS = ("aspect", "min_angle", "max_angle", "warping", "scaled_jacobian")
EXTRA_QUAD = ("stretch", "taper")


def _unit(v):
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.where(n == 0, 1.0, n)


def _corner_metrics(P):
    """P: (m,k,3) polygon corners. Returns edge lengths, angles (deg), scaled jac per corner, normal."""
    k = P.shape[1]
    nxt = np.roll(P, -1, axis=1) - P  # edge i: p[i+1]-p[i]
    prv = np.roll(P, 1, axis=1) - P  # to previous corner
    if k == 4:
        nrm = _unit(np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 1]))
    else:
        nrm = _unit(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]))
    ln, lp = np.linalg.norm(nxt, axis=-1), np.linalg.norm(prv, axis=-1)
    denom = np.where(ln * lp == 0, 1.0, ln * lp)
    cr = np.cross(nxt, prv)
    sj = np.einsum("mkc,mc->mk", cr, nrm) / denom
    cosang = np.clip(np.einsum("mkc,mkc->mk", nxt, prv) / denom, -1, 1)
    return ln, np.degrees(np.arccos(cosang)), sj, nrm


def _tri_normal(a, b, c):
    return _unit(np.cross(b - a, c - a))


def quad_quality(P):
    """Quality arrays for quads, P (m,4,3)."""
    P = np.asarray(P, float)
    ln, ang, sj, _ = _corner_metrics(P)
    r = 0.5 * (ln[:, 0] + ln[:, 2]) / (0.5 * (ln[:, 1] + ln[:, 3]))
    p0, p1, p2, p3 = (P[:, i] for i in range(4))

    def dang(n1, n2):
        return np.degrees(np.arccos(np.clip(np.einsum("mc,mc->m", n1, n2), -1, 1)))

    w1 = dang(_tri_normal(p0, p1, p2), _tri_normal(p0, p2, p3))
    w2 = dang(_tri_normal(p0, p1, p3), _tri_normal(p1, p2, p3))
    x1 = np.linalg.norm(p1 + p2 - p0 - p3, axis=1)
    x2 = np.linalg.norm(p2 + p3 - p0 - p1, axis=1)
    x12 = np.linalg.norm(p0 - p1 + p2 - p3, axis=1)
    return {
        "aspect": ln.max(1) / ln.min(1),
        "stretch": np.maximum(r, 1 / r),
        "min_angle": ang.min(1),
        "max_angle": ang.max(1),
        "warping": np.maximum(w1, w2),
        "scaled_jacobian": sj.min(1),
        "taper": x12 / np.minimum(x1, x2),
        "edges": ln,
    }


def tri_quality(P):
    """Quality arrays for triangles, P (m,3,3)."""
    P = np.asarray(P, float)
    ln, ang, sj, _ = _corner_metrics(P)
    return {
        "aspect": ln.max(1) / ln.min(1),
        "min_angle": ang.min(1),
        "max_angle": ang.max(1),
        "warping": np.zeros(len(P)),
        "scaled_jacobian": sj.min(1),
        "edges": ln,
    }


THRESHOLDS = (
    ("aspect", ">", (5, 10, 20)),
    ("min_angle", "<", (30, 15)),
    ("max_angle", ">", (150, 170)),
    ("warping", ">", (10,)),
    ("scaled_jacobian", "<", (0.3, 0.0)),
)


def _stats(x):
    if len(x) == 0:
        return None
    q = np.percentile(x, [50, 95, 99])
    return {"min": float(x.min()), "median": float(q[0]), "p95": float(q[1]),
            "p99": float(q[2]), "max": float(x.max())}


def _summarise(q):
    out = {"count": int(len(q["aspect"]))}
    for m in METRICS + tuple(k for k in EXTRA_QUAD if k in q):
        out[m] = _stats(q[m])
    out["beyond"] = {f"{m}{op}{t}": int((q[m] > t if op == ">" else q[m] < t).sum())
                     for m, op, ts in THRESHOLDS for t in ts}
    e = q["edges"].ravel()
    out["edge_length"] = {"median": float(np.median(e)), "p5": float(np.percentile(e, 5)),
                          "p95": float(np.percentile(e, 95))}
    return out


def analyse(nodes, quads, tris, qlabels, tlabels, legend):
    """Return (report dict, per-block quality dicts). labels are int arrays aligned with quads/tris."""
    nodes = np.asarray(nodes, float)
    quads = np.asarray(quads, int).reshape(-1, 4)
    tris = np.asarray(tris, int).reshape(-1, 3)
    qq = quad_quality(nodes[quads]) if len(quads) else None
    tq = tri_quality(nodes[tris]) if len(tris) else None
    qlabels, tlabels = np.asarray(qlabels, int), np.asarray(tlabels, int)

    def sel(q, mask):
        return {k: v[mask] for k, v in q.items()}

    rep = {"n_nodes": int(len(nodes)), "n_quads": int(len(quads)), "n_tris": int(len(tris)),
           "tri_to_quad_ratio": (len(tris) / len(quads)) if len(quads) else None,
           "surface_legend": legend, "overall": {}, "by_surface": {}}
    both = [q for q in (qq, tq) if q is not None]
    allq = {k: np.concatenate([b[k].ravel() if k == "edges" else b[k] for b in both])
            for k in METRICS + ("edges",)}
    rep["overall"]["all"] = _summarise(allq)
    if qq is not None:
        rep["overall"]["quad"] = _summarise(qq)
    if tq is not None:
        rep["overall"]["tri"] = _summarise(tq)
    for sid, name in legend.items():
        sid = int(sid)
        d = {}
        parts = []
        if qq is not None and (qlabels == sid).any():
            d["quad"] = _summarise(sel(qq, qlabels == sid))
            parts.append(sel(qq, qlabels == sid))
        if tq is not None and (tlabels == sid).any():
            d["tri"] = _summarise(sel(tq, tlabels == sid))
            parts.append(sel(tq, tlabels == sid))
        if parts:
            d["n_quads"] = int((qlabels == sid).sum()) if qq is not None else 0
            d["n_tris"] = int((tlabels == sid).sum()) if tq is not None else 0
            d["tri_to_quad_ratio"] = d["n_tris"] / d["n_quads"] if d["n_quads"] else None
            d["all"] = _summarise({k: np.concatenate([p[k].ravel() if k == "edges" else p[k] for p in parts])
                                   for k in METRICS + ("edges",)})
            rep["by_surface"][name] = d
    # span distribution
    ext = nodes.max(0) - nodes.min(0)
    ax = int(np.argmax(ext))
    cent = np.concatenate([nodes[quads].mean(1)[:, ax] if len(quads) else np.empty(0),
                           nodes[tris].mean(1)[:, ax] if len(tris) else np.empty(0)])
    edges = np.linspace(nodes[:, ax].min(), nodes[:, ax].max(), 11)
    cnt, _ = np.histogram(cent, edges)
    rep["span"] = {"axis": "xyz"[ax], "bin_edges": edges.tolist(), "element_count": cnt.tolist()}
    return rep, (qq, tq)


def write_vtu(path, nodes, quads, tris, qlabels, tlabels, qq, tq):
    import meshio

    cells, data = [], {k: [] for k in ("surface",) + METRICS + EXTRA_QUAD}
    for conn, lab, q, kind in ((quads, qlabels, qq, "quad"), (tris, tlabels, tq, "triangle")):
        if q is None or len(conn) == 0:
            continue
        cells.append((kind, np.asarray(conn, np.int64)))
        data["surface"].append(np.asarray(lab, np.int32))
        for k in METRICS + EXTRA_QUAD:
            data[k].append((q[k] if k in q else np.full(len(conn), np.nan)).astype(np.float32))
    data = {k: v for k, v in data.items() if v}
    meshio.Mesh(np.asarray(nodes, float), cells, cell_data=data).write(str(path))


def _finish(tag, out, nodes, quads, tris, ql, tl, legend, extra=None):
    rep, (qq, tq) = analyse(nodes, quads, tris, ql, tl, legend)
    if extra:
        rep.update(extra)
    out.mkdir(parents=True, exist_ok=True)
    write_vtu(out / f"{tag}.vtu", nodes, quads, tris, ql, tl, qq, tq)
    (out / f"quality_{tag}.json").write_text(json.dumps(rep, indent=1))
    print(f"{tag}: nodes={rep['n_nodes']} quads={rep['n_quads']} tris={rep['n_tris']} "
          f"ratio={rep['tri_to_quad_ratio']} -> {out / (tag + '.vtu')}")
    return rep


def run_v2(args):
    out = Path(args.out_dir)
    if args.from_npz:
        z = np.load(args.from_npz, allow_pickle=False)
        nodes, elems, surf = z["nodes_f64"], z["elements"], z["surfaces"].astype(str)
    else:
        from aeroelast.models.blade.numad.mesh_gen.gmsh_skin import generate_mesh_gmsh
        from aeroelast.models.blade.numad.objects.blade import Blade
        from v2_blade_mesh import REFERENCE_CHORD_DIVISIONS, REFERENCE_SPAN_ASPECT

        if args.element_size is not None and args.element_size >= 0.3:
            raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")
        blade = Blade()
        blade.read_yaml(args.yaml)
        sizing = ({"element_size": args.element_size, "h_range": (0.05, 0.25)} if args.element_size is not None
                  else {"chord_divisions": REFERENCE_CHORD_DIVISIONS, "span_aspect": REFERENCE_SPAN_ASPECT})
        t0 = time.time()
        raw = generate_mesh_gmsh(blade.windio2_topology, **sizing)
        print(f"generate_mesh_gmsh: {time.time() - t0:.0f} s")
        nodes = np.asarray(raw.nodes, np.float64)
        elems = np.asarray(raw.elements)
        surf = np.asarray(raw.surfaces).astype(str)
        npz = out / (args.tag + "_raw.npz")
        out.mkdir(parents=True, exist_ok=True)
        np.savez(npz, nodes_f64=nodes, nodes_longdouble=np.asarray(raw.nodes),
                 elements=elems, surfaces=surf)
        print(f"raw saved: {npz}")
    names = sorted(set(surf.tolist()))
    legend = dict(enumerate(names))
    lab = np.array([names.index(s) for s in surf], np.int32)
    return _finish(args.tag, out, nodes, elems, np.empty((0, 3), int), lab, np.empty(0, int), legend)


def write_summary(out_dir, tags):
    """Markdown comparison table of quality_<tag>.json files (overall + per surface)."""
    rows = ["| mesh | group | quads | tris | AR med/p99/max | minang min/med | maxang med/max | warp p99/max | sjac min | AR>5 | AR>10 | minang<30 | sjac<0.3 |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag in tags:
        r = json.loads((Path(out_dir) / f"quality_{tag}.json").read_text())
        groups = [("ALL", r["overall"]["all"], r["n_quads"], r["n_tris"])]
        groups += [(k, v["all"], v["n_quads"], v["n_tris"]) for k, v in r["by_surface"].items()]
        for name, g, nq, nt in groups:
            b = g["beyond"]
            rows.append(
                f"| {tag} | {name} | {nq} | {nt} | {g['aspect']['median']:.2f}/{g['aspect']['p99']:.2f}/{g['aspect']['max']:.2f}"
                f" | {g['min_angle']['min']:.1f}/{g['min_angle']['median']:.1f} | {g['max_angle']['median']:.1f}/{g['max_angle']['max']:.1f}"
                f" | {g['warping']['p99']:.2f}/{g['warping']['max']:.2f} | {g['scaled_jacobian']['min']:.3f}"
                f" | {b['aspect>5']} | {b['aspect>10']} | {b['min_angle<30']} | {b['scaled_jacobian<0.3']} |")
    rows.append("")
    for tag in tags:
        r = json.loads((Path(out_dir) / f"quality_{tag}.json").read_text())
        e = r["overall"]["all"]["edge_length"]
        rows.append(f"{tag}: edge length median/p5/p95 = {e['median']:.3f}/{e['p5']:.3f}/{e['p95']:.3f} m; "
                    f"tri/quad = {r['tri_to_quad_ratio']:.4f}; span axis {r['span']['axis']} bin counts {r['span']['element_count']}")
    (Path(out_dir) / "quality_summary.md").write_text("\n".join(rows) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="ver", required=True)
    p = sub.add_parser("v2")
    p.add_argument("yaml", nargs="?")
    p.add_argument("out_dir")
    p.add_argument("--element-size", type=float, default=None,
                   help="uniform element size [m, 0.05-0.25]; default: reference mesh chord_divisions 20, span_aspect 1")
    p.add_argument("--from-npz")
    p.add_argument("--tag", default="iea15_v2_full")
    ps = sub.add_parser("summary")
    ps.add_argument("out_dir")
    ps.add_argument("tags", nargs="+")
    args = ap.parse_args(argv)
    if args.ver == "summary":
        write_summary(args.out_dir, args.tags)
    else:
        run_v2(args)


if __name__ == "__main__":
    sys.exit(main())
