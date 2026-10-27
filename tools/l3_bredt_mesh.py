#!/usr/bin/env python
"""Thin-wall multi-cell Bredt GJ of our shell mesh section, vs the yaml Bredt, the S-1 FE GJ and the deck (odd/tasks/l3-closure.md, T9).

Discriminates where the S-1 GJ deficit of a blade enters: the mesh geometry/laminate (mesh Bredt below
the yaml Bredt of ``tools/l3_bredt_yaml.py``) or the element / extraction (mesh Bredt ~ yaml Bredt, FE
below it). Generalises ``model_bredt`` of ``tools/bredt_gj.py`` (branch ``fix/cap-arcs``, IEA-15 only)
to any windIO v1 yaml. Conventions:

* mesh: ``BladeMesh(yaml, element_size, skin_reference='mid')`` as ``tools/l3_s1_blade.py``; the section
  is the element row whose centroid z is nearest ``eta * blade_length``, cut at that z (each element
  polygon gives one wall segment on its mid-surface);
* wall G t = laminate A66 of the element's section layup (Qbar66 per ply; G12 t for 0 deg plies);
* cells = bounded faces of the planar wall graph, ranked by chordwise position from the TE; regions as
  in ``l3_bredt_yaml._regions``: exterior wall per cell (TE cell, box, LE cell) and each shared wall
  (web), with their share of the summed integral ds/(G t);
* deck GJ_sc = K66 condensed over the two shear strains (``l3_bredt_yaml.condense_shear``).

  python tools/l3_bredt_mesh.py --yaml tests/IEA-3.4-130-RWT.yaml \\
      --beamdyn tests/IEA34/reference/IEA-3.4-130-RWT_BeamDyn_blade.dat --blade-length 63.0 \\
      --element-size 0.15 --stations 0.3 0.5 0.7 --s1-json logs/check/l3/t9/s1-iea34-es0.15.json \\
      --yaml-bredt-json logs/check/l3/t9/bredt/iea34.json --out logs/check/l3/t9/bredt/mesh-iea34.json
  python tools/l3_bredt_mesh.py --selftest
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from l3_bredt_yaml import condense_shear, interp_deck, load_beamdyn_6x6, multicell_bredt  # noqa: E402

ROW_TOL = 1e-5  # [m] elements of one span row share their centroid z to this tolerance


# =============================================================================================
# planar wall graph -> multi-cell Bredt (ported from tools/bredt_gj.py, fix/cap-arcs)
# =============================================================================================
def polygon_area(p: np.ndarray) -> float:
    return 0.5 * float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))


def graph_faces(xy, edges):
    """Faces of a planar straight-line graph: list of (half-edge list, signed area); bounded faces have area > 0."""
    nbr = defaultdict(list)
    for u, v in edges:
        nbr[u].append(v)
        nbr[v].append(u)
    order = {u: sorted(vs, key=lambda v: np.arctan2(*(xy[v] - xy[u])[::-1])) for u, vs in nbr.items()}
    seen, faces = set(), []
    for u, v in list(edges) + [(v, u) for u, v in edges]:
        if (u, v) in seen:
            continue
        face, (a, b) = [], (u, v)
        while (a, b) not in seen:
            seen.add((a, b))
            face.append((a, b))
            ring = order[b]
            a, b = b, ring[(ring.index(a) - 1) % len(ring)]  # next edge clockwise: face on the left
        faces.append((face, polygon_area(xy[[h[0] for h in face]])))
    return faces


def graph_bredt(xy, edges, gt) -> dict:
    """Multi-cell Bredt GJ (unit twist rate) of a wall graph; cells ranked by mean x (cell 0 = lowest x)."""
    xy, gt = np.asarray(xy, float), np.asarray(gt, float)
    bounded = [(f, a) for f, a in graph_faces(xy, edges) if a > 1e-9]
    if not bounded:
        raise ValueError("the wall graph closes no cell")
    cx = [float(np.mean(xy[[h[0] for h in f]][:, 0])) for f, _ in bounded]
    rank = {c: r for r, c in enumerate(np.argsort(cx))}
    n = len(bounded)
    index = {frozenset(e): k for k, e in enumerate(edges)}
    owner = defaultdict(list)
    for c, (face, _) in enumerate(bounded):
        for a, b in face:
            owner[index[frozenset((a, b))]].append(rank[c])
    length = np.array([np.linalg.norm(xy[u] - xy[v]) for u, v in edges])
    flex = length / gt
    f = np.zeros((n, n))
    ext = np.zeros(n)
    shared = defaultdict(lambda: [0.0, 0.0])
    open_edges = []
    for e in range(len(edges)):
        cells = owner.get(e, [])
        if len(cells) == 1:
            f[cells[0], cells[0]] += flex[e]
            ext[cells[0]] += flex[e]
        elif len(cells) == 2 and cells[0] != cells[1]:
            i, j = sorted(cells)
            f[i, i] += flex[e]
            f[j, j] += flex[e]
            f[i, j] -= flex[e]
            f[j, i] -= flex[e]
            shared[(i, j)][0] += flex[e]
            shared[(i, j)][1] += length[e]
        else:
            open_edges.append(e)
    area = np.zeros(n)
    for c, (_, a) in enumerate(bounded):
        area[rank[c]] = a
    q = np.linalg.solve(f, 2 * area)
    reg = {"te_cell_exterior": ext[0]}
    if n > 2:
        reg["box_exterior"] = float(ext[1:-1].sum())
    if n > 1:
        reg["le_cell_exterior"] = float(ext[-1])
    for (i, j), (fl, _) in sorted(shared.items()):
        reg[f"shared:{i}-{j}"] = fl
    total = sum(reg.values())
    return {"gj": float(2 * area @ q), "areas": area.tolist(), "n_cells": n, "n_edges": len(edges),
            "n_open": len(open_edges), "open_length": float(length[open_edges].sum()) if open_edges else 0.0,
            "regions": {k: {"flex": float(v), "share": float(v / total)} for k, v in reg.items()},
            "te_share": float(reg["te_cell_exterior"] / total),
            "shared_length": {f"{i}-{j}": v[1] for (i, j), v in sorted(shared.items())}}


def laminate_a66(layup, mats) -> float:
    """A66 [N/m] of a layup [(material, t, angle_deg)]."""
    a66 = 0.0
    for material, t, angle in layup:
        el = mats[material]["elastic"]
        e = el["E"] if isinstance(el["E"], (list, tuple)) else [el["E"]] * 3
        nu = el["nu"] if isinstance(el["nu"], (list, tuple)) else [el["nu"]] * 3
        g = el["G"] if isinstance(el["G"], (list, tuple)) else [el["G"]] * 3
        e1, e2, nu12 = e[0], e[1], nu[0]
        g12 = g[0] if g[0] else e1 / (2 * (1 + nu12))
        den = 1 - nu12 * nu12 * e2 / e1
        q11, q22, q12 = e1 / den, e2 / den, nu12 * e2 / den
        c, s = np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))
        a66 += ((q11 + q22 - 2 * q12 - 2 * g12) * s * s * c * c + g12 * (s ** 4 + c ** 4)) * t
    return a66


# =============================================================================================
# mesh section
# =============================================================================================
def _section_segment(points: np.ndarray, z: float) -> np.ndarray:
    hits = []
    for a, b in zip(points, np.roll(points, -1, axis=0), strict=True):
        if min(a[2], b[2]) <= z <= max(a[2], b[2]) and abs(a[2] - b[2]) > 1e-10:
            hits.append(a + (b - a) * (z - a[2]) / (b[2] - a[2]))
    if len(hits) != 2:
        raise ValueError(f"expected two row intersections, got {len(hits)}")
    return np.asarray(hits)


def _chord_frame(blade, z: float):
    """TE point and (chord, normal) basis of the NuMAD section at span z (same as t37_section_compare.frame_at)."""
    geo = blade.geometry
    span = np.asarray(geo.ispan)
    k = int(np.clip(np.searchsorted(span, z) - 1, 0, len(span) - 2))
    f = (z - span[k]) / (span[k + 1] - span[k])
    oml = (1 - f) * geo.coordinates[:, :, k] + f * geo.coordinates[:, :, k + 1]
    key = (1 - f) * blade.keypoints.key_points[:, :, k] + f * blade.keypoints.key_points[:, :, k + 1]
    te = (key[0, :2] + key[9, :2]) / 2
    chord = oml[geo.LEindex, :2] - te
    chord /= np.linalg.norm(chord)
    return te, np.array([chord, [-chord[1], chord[0]]]).T


def build_mesh(yaml: str, element_size: float):
    from aeroelast.core.mesh.generators import BladeMesh

    gen = BladeMesh(yaml_file=yaml, element_size=element_size, skin_reference="mid")
    mesh = gen.generate(renumber="rcm")
    data = gen.numad_mesh_data
    sections = {s["elementSet"]: s for s in data["sections"]}
    mats = {m["name"]: m for m in data["materials"]}
    elem_section = {el.id: name for name, eset in mesh.element_sets.items() if name in sections
                    for el in eset.elements}
    coords, idx = mesh.coords_array, mesh.node_id_to_index
    elements = []
    for el in mesh.elements:
        p = coords[[idx[n] for n in el.node_ids]]
        elements.append((elem_section[el.id], p, float(p[:, 2].mean())))
    return gen, sections, mats, elements


def mesh_section(gen, sections, mats, elements, rows, z_target: float) -> dict:
    z = float(rows[np.argmin(abs(rows - z_target))])
    te, basis = _chord_frame(gen.numad_blade, z)
    # every element the plane z crosses: node rings of a pre-bent blade are not at constant z, so the
    # plane through a row's centroids also cuts elements of the neighbouring rows (counted as off-row)
    cut = [(name, p, zm) for name, p, zm in elements if p[:, 2].min() < z < p[:, 2].max()]
    off_row = Counter(name for name, _, zm in cut if abs(zm - z) >= ROW_TOL)
    segs = [((_section_segment(p, z)[:, :2] - te) @ basis, laminate_a66(sections[name]["layup"], mats))
            for name, p, _ in cut]
    pts = np.concatenate([s[0] for s in segs])
    keys, inv = np.unique(np.round(pts / 1e-6).astype(np.int64), axis=0, return_inverse=True)
    inv = inv.ravel()
    edges, gt = [], []
    for k, (_, a66) in enumerate(segs):
        u, v = int(inv[2 * k]), int(inv[2 * k + 1])
        if u != v:
            edges.append((u, v))
            gt.append(a66)
    out = graph_bredt(keys * 1e-6, edges, gt)
    out.update({"z_row": z, "n_row_elements": len(segs), "n_off_row": sum(off_row.values()),
                "off_row_sections": dict(off_row)})
    return out


# =============================================================================================
# selftest: closed forms and the yaml-side multi-cell solver
# =============================================================================================
def selftest() -> int:
    ok = True
    r, gt, n = 1.3, 2.0e7, 2000
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xy = np.c_[r * np.cos(th), r * np.sin(th)]
    res = graph_bredt(xy, [(i, (i + 1) % n) for i in range(n)], [gt] * n)
    ref = 2 * np.pi * r ** 3 * gt  # thin circular tube: 4 A^2 G t / P
    err = res["gj"] / ref - 1
    ok &= abs(err) < 1e-5 and res["n_cells"] == 1 and res["n_open"] == 0
    print(f"selftest circle: GJ {res['gj']:.6e} vs 2 pi r^3 G t {ref:.6e} ({100 * err:+.5f} %)")
    # two cells (4 x 2 rectangle, web at x = 1.5), plus a dangling open strip, vs multicell_bredt
    # (its loop starts on the first side walked from the TE, so the web runs vertex 1 -> 4)
    xy = np.array([[0, -1], [1.5, -1], [4, -1], [4, 1], [1.5, 1], [0, 1], [4, 1.5]], float)
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (1, 4), (3, 6)]
    gts = [1e7, 2e7, 3e7, 2e7, 1e7, 4e7, 5e7, 9e9]
    res = graph_bredt(xy, edges, gts)
    loop = xy[[0, 1, 2, 3, 4, 5]]
    ref = multicell_bredt(loop, [1e7, 2e7, 3e7, 2e7, 1e7, 4e7], [0, 1, 1, 1, 0, 0],
                          [{"i_s": 1, "i_p": 4, "gt": 5e7}])
    err = res["gj"] / ref["gj"] - 1
    ok &= abs(err) < 1e-12 and res["n_cells"] == 2 and res["n_open"] == 1
    print(f"selftest two cells: graph {res['gj']:.6e} vs multicell_bredt {ref['gj']:.6e} ({100 * err:+.2e} %), "
          f"open walls {res['n_open']}, TE share {res['te_share']:.3f} vs {ref['flex_ext'][0] / (ref['flex_ext'].sum() + ref['webs'][0]['flex']):.3f}")
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# =============================================================================================
# CLI
# =============================================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--yaml")
    ap.add_argument("--beamdyn")
    ap.add_argument("--blade-length", type=float, help="length that scales the deck eta to the mesh z")
    ap.add_argument("--element-size", type=float, default=0.15)
    ap.add_argument("--stations", type=float, nargs="+")
    ap.add_argument("--s1-json", help="tools/l3_s1_blade.py output: FE GJ interpolated at the row z")
    ap.add_argument("--yaml-bredt-json", help="tools/l3_bredt_yaml.py output at the same stations")
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    if args.selftest:
        return selftest()
    if not all((args.yaml, args.beamdyn, args.blade_length, args.stations, args.out)):
        ap.error("--yaml, --beamdyn, --blade-length, --stations and --out are required")
    if args.element_size >= 0.3:
        raise SystemExit("element_size >= 0.3 is not converged (docs/blade_mesh_size_rule.md)")

    etas, ks = load_beamdyn_6x6(args.beamdyn)
    fe = None
    if args.s1_json:
        g = json.loads(Path(args.s1_json).read_text())["metrics"]["GJ"]
        fe = (np.asarray(g["z"], float), np.asarray(g["ours"], float))
    yb = {}
    if args.yaml_bredt_json:
        yb = {round(s["eta"], 4): s["offset"] for s in json.loads(Path(args.yaml_bredt_json).read_text())["stations"]}
    gen, sections, mats, elements = build_mesh(args.yaml, args.element_size)
    rows = np.unique(np.round([zm for _, _, zm in elements], 6))
    print(f"mesh {Path(args.yaml).name} es {args.element_size}: {len(elements)} elements, {len(rows)} rows; "
          "% vs deck GJ_sc; FE = S-1 GJ linearly interpolated at the row z")
    print(f"{'eta':>5} {'z_row':>7} {'cells':>5} {'open':>4} {'mesh':>11} {'%mesh':>7} {'%yaml':>7} {'%FE':>7} "
          f"{'mesh/yaml':>9} {'FE/mesh':>8} {'TE m/y':>11}")
    out = []
    for eta in args.stations:
        sec = mesh_section(gen, sections, mats, elements, rows, eta * args.blade_length)
        gjsc = condense_shear(interp_deck(etas, ks, sec["z_row"] / args.blade_length))
        y = yb.get(round(eta, 4))
        fe_gj = float(np.interp(sec["z_row"], *fe)) if fe is not None else float("nan")
        y_gj = y["gj"] if y else float("nan")
        rec = {"eta": eta, "gj_sc": gjsc, "mesh": sec, "yaml_bredt": y_gj, "fe_gj": fe_gj,
               "pct_mesh": 100 * (sec["gj"] / gjsc - 1), "pct_yaml": 100 * (y_gj / gjsc - 1),
               "pct_fe": 100 * (fe_gj / gjsc - 1), "mesh_over_yaml": 100 * (sec["gj"] / y_gj - 1),
               "fe_over_mesh": 100 * (fe_gj / sec["gj"] - 1)}
        te_y = y["te_share"] if y else float("nan")
        print(f"{eta:5.2f} {sec['z_row']:7.3f} {sec['n_cells']:5d} {sec['n_open']:4d} {sec['gj']:11.4e} "
              f"{rec['pct_mesh']:+7.1f} {rec['pct_yaml']:+7.1f} {rec['pct_fe']:+7.1f} {rec['mesh_over_yaml']:+9.1f} "
              f"{rec['fe_over_mesh']:+8.1f} {sec['te_share']:5.3f}/{te_y:5.3f}"
              + (f"  off-row cut {sec['n_off_row']}: {sec['off_row_sections']}" if sec["n_off_row"] else ""))
        out.append(rec)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps({"yaml": args.yaml, "beamdyn": args.beamdyn, "element_size": args.element_size,
                                          "blade_length": args.blade_length, "stations": out}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
